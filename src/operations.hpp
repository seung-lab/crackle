#ifndef __CRACKLE_OPERATIONS_HPP__
#define __CRACKLE_OPERATIONS_HPP__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <span>
#include <type_traits>
#include <unordered_map>
#include <mutex>

#include "robin_hood.hpp"

#include "cc3d.hpp"
#include "crc.hpp"
#include "header.hpp"
#include "crackcodes.hpp"
#include "lib.hpp"
#include "labels.hpp"
#include "pins.hpp"
#include "markov.hpp"
#include "dual_graph.hpp"
#include "threadpool.hpp"

namespace crackle {
namespace operations {

template <typename LABEL>
std::unordered_map<uint64_t, std::vector<uint16_t>>
point_cloud(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const std::optional<uint64_t> label = std::nullopt,
	size_t parallel = 1
) {

	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	z_start = std::max(std::min(z_start, static_cast<int64_t>(header.sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(header.sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(header.sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const int64_t szr = z_end - z_start;

	const uint64_t voxels = (
		static_cast<uint64_t>(header.sx) 
		* static_cast<uint64_t>(header.sy) 
		* static_cast<uint64_t>(szr)
	);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, std::vector<uint16_t>>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, std::vector<uint16_t>> ptc;
	ptc.reserve(num_labels / (header.sz / szr));

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(szr));

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel);
	std::vector<std::vector<uint32_t>> ccls(parallel);

	for (size_t t = 0; t < parallel; t++) {
		vcgs[t].resize(header.sx * header.sy);
		ccls[t].resize(header.sx * header.sy);
	}

	std::mutex mtx;

	for (size_t z = z_start, i = 0; i < crack_codes.size(); z++, i++) {
			pool.enqueue([&,z,i](size_t t){
				std::vector<uint8_t>& vcg = vcgs[t];
				std::vector<uint32_t>& ccl = ccls[t];
				auto& crack_code = crack_codes[i];

				crack_code_to_vcg(
					/*code=*/crack_code,
					/*sx=*/header.sx, /*sy=*/header.sy,
					/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
					/*markov_model=*/markov_model,
					/*slice_edges=*/vcg.data()
				);

				uint64_t N = 0;
				crackle::cc3d::color_connectivity_graph<uint32_t>(
					vcg, header.sx, header.sy, 1, ccl.data(), N
				);

				std::vector<LABEL> label_map = decode_label_map<LABEL>(
					header, binary, ccl.data(), N, z, z+1
				);

				uint64_t label_i = 0;

				auto ptc_ccls = crackle::dual_graph::extract_contours(vcg, ccl, N, header.sx, header.sy);

				std::unique_lock<std::mutex> lock(mtx);

				for (auto ptc_ccl : ptc_ccls) {
					uint64_t current_label = label_map[label_i];
					
					if (label.has_value() && current_label != *label) {
						label_i++;
						continue;
					}

					std::vector<uint16_t>& label_points = ptc[current_label];

					for (uint32_t loc : ptc_ccl) {
						uint16_t y = loc / header.sx;
						uint16_t x = loc - (header.sx * y);

						label_points.push_back(x);
						label_points.push_back(y);
						label_points.push_back(static_cast<uint16_t>(z));
					}

					label_i++;
				}
		});
	}

	pool.join();

	return ptc;
}

auto point_cloud(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const std::optional<uint64_t> label = std::nullopt,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return point_cloud<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, label, parallel
		);
	}
	else if (header.data_width == 2) {
		return point_cloud<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, label, parallel
		);
	}
	else if (header.data_width == 4) {
		return point_cloud<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, label, parallel
		);
	}
	else {
		return point_cloud<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, label, parallel
		);
	}
}

auto point_cloud(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	const std::optional<uint64_t> label = std::nullopt,
	size_t parallel = 1
) {
	return point_cloud(
		buffer.data(),
		buffer.size(),
		z_start, z_end, label, parallel
	);
}

template <typename LABEL>
std::unordered_map<uint64_t, uint64_t>
voxel_counts(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {

	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	z_start = std::max(std::min(z_start, static_cast<int64_t>(header.sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(header.sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(header.sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const int64_t szr = z_end - z_start;

	const uint64_t voxels = (
		static_cast<uint64_t>(header.sx) 
		* static_cast<uint64_t>(header.sy) 
		* static_cast<uint64_t>(szr)
	);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, uint64_t>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, uint64_t> cts;
	cts.reserve(num_labels / (header.sz / szr));

	const uint64_t sxy = header.sx * header.sy;

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(szr));

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel);
	std::vector<std::vector<uint32_t>> ccls(parallel);

	for (size_t t = 0; t < parallel; t++) {
		vcgs[t].resize(sxy);
		ccls[t].resize(sxy);
	}

	std::mutex mtx;

	for (size_t z = z_start, i = 0; i < crack_codes.size(); z++, i++) {
			pool.enqueue([&,z,i](size_t t){
				std::vector<uint8_t>& vcg = vcgs[t];
				std::vector<uint32_t>& ccl = ccls[t];
				auto& crack_code = crack_codes[i];

				crack_code_to_vcg(
					/*code=*/crack_code,
					/*sx=*/header.sx, /*sy=*/header.sy,
					/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
					/*markov_model=*/markov_model,
					/*slice_edges=*/vcg.data()
				);

				uint64_t N = 0;
				crackle::cc3d::color_connectivity_graph<uint32_t>(
					vcg, header.sx, header.sy, 1, ccl.data(), N
				);

				std::vector<LABEL> label_map = decode_label_map<LABEL>(
					header, binary, ccl.data(), N, z, z+1
				);
				std::vector<uint64_t> subcounts(N);

				for (uint64_t i = 0; i < sxy; i++) {
					subcounts[ccl[i]]++;
				}

				std::unique_lock<std::mutex> lock(mtx);
				for (uint64_t i = 0; i < N; i++) {
						cts[label_map[i]] += subcounts[i];
				}
		});
	}

	pool.join();

	return cts;
}

auto voxel_counts(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return voxel_counts<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 2) {
		return voxel_counts<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 4) {
		return voxel_counts<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else {
		return voxel_counts<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
}

auto voxel_counts(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1
) {
	return voxel_counts(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}

};
};

#endif