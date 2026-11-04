#ifndef __CRACKLE_OPERATIONS_HPP__
#define __CRACKLE_OPERATIONS_HPP__

#include <concepts>
#include <algorithm>
#include <array>
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

crackle::CrackleHeader get_header(
	const unsigned char* buffer, 
	const size_t num_bytes
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

	return header;
}

int64_t get_szr(
	const crackle::CrackleHeader& header,
	int64_t z_start, 
	int64_t z_end
) {
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

	return z_end - z_start;
}

uint64_t get_voxels(
	const crackle::CrackleHeader& header,
	int64_t z_start, 
	int64_t z_end
) {
	const int64_t szr = get_szr(header, z_start, z_end);

	return (
		static_cast<uint64_t>(header.sx) 
		* static_cast<uint64_t>(header.sy) 
		* static_cast<uint64_t>(szr)
	);
}

template <typename LABEL>
void for_each_z_parallel(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start,
	int64_t z_end,
	size_t parallel,
	const std::function<void(
		std::vector<uint8_t>&, 
		std::vector<uint32_t>&,
		uint64_t,
		std::vector<LABEL>&,
		int64_t
	)>& fn
) {
	const CrackleHeader header = get_header(buffer, num_bytes);

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
		return;
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

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

			fn(vcg, ccl, N, label_map, z);
		});
	}

	pool.join();

	return;
}

template <typename LABEL>
void for_each_z_vcg_parallel(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start,
	int64_t z_end,
	size_t parallel,
	const std::function<void(
		std::vector<uint8_t>&, 
		std::vector<LABEL>&,
		int64_t
	)>& fn
) {
	const CrackleHeader header = get_header(buffer, num_bytes);

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
		return;
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(szr));

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel);

	for (size_t t = 0; t < parallel; t++) {
		vcgs[t].resize(header.sx * header.sy);
	}

	for (size_t z = z_start, i = 0; i < crack_codes.size(); z++, i++) {
		pool.enqueue([&,z,i](size_t t){
			std::vector<uint8_t>& vcg = vcgs[t];
			auto& crack_code = crack_codes[i];

			crack_code_to_vcg(
				/*code=*/crack_code,
				/*sx=*/header.sx, /*sy=*/header.sy,
				/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model,
				/*slice_edges=*/vcg.data()
			);

			std::vector<LABEL> label_map = decode_label_map<LABEL>(
				header, binary, nullptr, 0, z, z+1
			);

			fn(vcg, label_map, z);
		});
	}

	pool.join();

	return;
}

typedef std::unordered_map<uint64_t, std::vector<uint16_t>> pointcloud_t;

template <typename LABEL>
pointcloud_t _point_cloud_binary_image(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start,
	int64_t z_end,
	const std::optional<uint64_t> label = std::nullopt,
	const bool skip_background = false,
	const bool skip_last_2d = false
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);

	if (voxels == 0) {
		return pointcloud_t();
	}

	const int64_t sx = header.sx;
	const int64_t sy = header.sy;
	const int64_t sxy = header.sx * header.sy;

	const bool selective = label.has_value();
	std::span<const unsigned char> binary(buffer, num_bytes);	
	std::vector<std::vector<uint16_t>> ptc(2);

	std::vector<std::vector<uint8_t>> vcgs(2);
	vcgs[0].resize(sxy);
	vcgs[1].resize(sxy);
	std::vector<std::vector<uint64_t>> label_maps(2);

	std::vector<uint64_t> uniq = crackle::labels::unique(binary);

	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	crack_code_to_vcg(
		/*code=*/crack_codes[0],
		/*sx=*/sx, /*sy=*/sy,
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
		/*markov_model=*/markov_model,
		/*slice_edges=*/vcgs[0].data()
	);

	auto ptc_ccls = crackle::dual_graph::extract_contours_binary_image(vcgs[0], header.sx, header.sy);
	label_maps[0] = decode_label_map<uint64_t>(
		header, binary, nullptr, 0, z_start, z_start+1
	);

	auto copy_contour_fn = [&](
		const std::vector<std::vector<uint32_t>>& ptc_ccls,
		const std::vector<uint64_t>& label_map,
		const uint16_t z
	) {
		uint64_t label_i = (label_map[0] == uniq[1])
			? 1
			: 0;

		for (auto ptc_ccl : ptc_ccls) {
			uint64_t current_label = uniq[label_i];

			if ((skip_background && current_label == 0) 
				|| (selective && current_label != *label)
			) {
				label_i = !label_i;
				continue;
			}

			std::vector<uint16_t>& label_points = ptc[label_i];

			for (uint32_t loc : ptc_ccl) {
				uint16_t y = loc / header.sx;
				uint16_t x = loc - (header.sx * y);

				label_points.push_back(x);
				label_points.push_back(y);
				label_points.push_back(z);
			}

			label_i = !label_i;
		}
	};

	copy_contour_fn(ptc_ccls, label_maps[0], z_start);

	for (int64_t z = z_start + 1, i = 1; i < (int64_t)crack_codes.size(); z++, i++) {
		int64_t bottom_j = i & 0b1;
		int64_t top_j = (i-1) & 0b1;

		crack_code_to_vcg(
			/*code=*/crack_codes[i],
			/*sx=*/sx, /*sy=*/sy,
			/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
			/*markov_model=*/markov_model,
			/*slice_edges=*/vcgs[bottom_j].data()
		);

		const auto& top_vcg = vcgs[top_j];
		auto& bottom_vcg = vcgs[bottom_j];

		label_maps[bottom_j] = decode_label_map<uint64_t>(
			header, binary, nullptr, 0, z, z+1
		);

		const auto& top_label_map = label_maps[top_j];
		const auto& bottom_label_map = label_maps[bottom_j];

		if (!skip_last_2d || (z < z_end - 1)) {
			copy_contour_fn(
				crackle::dual_graph::extract_contours_binary_image(bottom_vcg, header.sx, header.sy),
				bottom_label_map,
				z
			);
		}

		const uint8_t top_offset = (top_label_map[0] == uniq[1]);
		const uint8_t bottom_offset = (bottom_label_map[0] == uniq[1]);

		for (int64_t y = 0; y < sy; y++) {
			for (int64_t x = 0; x < sx; x++) {
				int64_t loc = x + sx * y;

				const uint8_t top_val = (top_vcg[loc] - top_offset);
				const uint8_t bottom_val = (bottom_vcg[loc] - bottom_offset);

				uint64_t current_label = uniq[bottom_val];
				if (
					(skip_background && current_label == 0) 
					|| (selective && *label != current_label)
				) {
					continue;
				}

				if (top_val != bottom_val) {
					std::vector<uint16_t>& label_points = ptc[current_label];
					label_points.push_back(x);
					label_points.push_back(y);
					label_points.push_back(static_cast<uint16_t>(z));
				}
			}
		}
	}

	pointcloud_t return_ptc;
	return_ptc[uniq[0]] = std::move(ptc[0]);
	return_ptc[uniq[1]] = std::move(ptc[1]);
	return return_ptc;
}

template <typename LABEL>
pointcloud_t _point_cloud_multivalued_image(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start,
	int64_t z_end,
	const std::optional<std::vector<uint64_t>> labels = std::nullopt,
	const bool skip_background = false,
	const bool skip_last_2d = false
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);

	if (voxels == 0) {
		return pointcloud_t();
	}

	const int64_t sx = header.sx;
	const int64_t sy = header.sy;
	const int64_t sxy = header.sx * header.sy;

	const bool selective = labels.has_value();
	robin_hood::unordered_flat_set<uint64_t> labels_set;
	if (selective) {
		labels_set.reserve(labels->size());
		labels_set.insert(labels->begin(), labels->end());
	}

	std::vector<std::vector<uint8_t>> vcgs(2);
	std::vector<std::vector<uint32_t>> ccls(2);
	std::vector<std::vector<uint64_t>> label_maps(2);
	for (size_t t = 0; t < ccls.size(); t++) {
		vcgs[t].resize(sxy);
		ccls[t].resize(sxy);
	}

	std::span<const unsigned char> binary(buffer, num_bytes);	
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	crack_code_to_vcg(
		/*code=*/crack_codes[0],
		/*sx=*/sx, /*sy=*/sy,
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
		/*markov_model=*/markov_model,
		/*slice_edges=*/vcgs[0].data()
	);

	uint64_t N = 0;
	crackle::cc3d::color_connectivity_graph<uint32_t>(
		vcgs[0], sx, sy, 1, ccls[0].data(), N
	);

	label_maps[0] = decode_label_map<uint64_t>(
		header, binary, ccls[0].data(), N, z_start, z_start+1
	);

	pointcloud_t ptcs;

	auto copy_contour_fn = [&](
		const std::vector<std::vector<uint32_t>>& ptc_ccls,
		const std::vector<uint64_t>& label_map,
		const uint16_t z
	) {
		uint64_t label_i = 0;
		for (auto ptc_ccl : ptc_ccls) {
			const uint64_t current_label = label_map[label_i];

			if ((skip_background && current_label == 0) 
				|| (selective && !labels_set.contains(current_label))
			) {
				label_i++;
				continue;
			}

			std::vector<uint16_t>& label_points = ptcs[current_label];

			for (uint32_t loc : ptc_ccl) {
				uint16_t y = loc / header.sx;
				uint16_t x = loc - (header.sx * y);

				label_points.push_back(x);
				label_points.push_back(y);
				label_points.push_back(z);
			}

			label_i++;
		}
	};

	auto ptc_ccls = crackle::dual_graph::extract_contours(vcgs[0], ccls[0], N, header.sx, header.sy);
	copy_contour_fn(ptc_ccls, label_maps[0], z_start);

	for (int64_t z = z_start + 1, i = 1; i < (int64_t)crack_codes.size(); z++, i++) {
		int64_t bottom_j = i & 0b1;
		int64_t top_j = (i-1) & 0b1;

		const auto& top_vcg = vcgs[top_j];
		auto& bottom_vcg = vcgs[bottom_j];
		const auto& top_ccl = ccls[top_j];
		auto& bottom_ccl = ccls[bottom_j];

		crack_code_to_vcg(
			/*code=*/crack_codes[i],
			/*sx=*/sx, /*sy=*/sy,
			/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
			/*markov_model=*/markov_model,
			/*slice_edges=*/bottom_vcg.data()
		);

		uint64_t N = 0;
		crackle::cc3d::color_connectivity_graph<uint32_t>(
			bottom_vcg, sx, sy, 1, bottom_ccl.data(), N
		);

		label_maps[bottom_j] = decode_label_map<uint64_t>(
			header, binary, bottom_ccl.data(), N, z, z+1
		);

		const auto& top_label_map = label_maps[top_j];
		const auto& bottom_label_map = label_maps[bottom_j];

		if (!skip_last_2d || ((int)z < ((int)z_end - 1))) {
			copy_contour_fn(
				crackle::dual_graph::extract_contours(bottom_vcg, bottom_ccl, N, header.sx, header.sy),
				bottom_label_map,
				z
			);
		}

		for (int64_t y = 0; y < sy; y++) {
			for (int64_t x = 0; x < sx; x++) {
				const int64_t loc = x + sx * y;

				const uint64_t current_label = bottom_label_map[bottom_ccl[loc]];

				if (
					(skip_background && current_label == 0) 
					|| (selective && !labels_set.contains(current_label))
				) {
					continue;
				}

				const uint64_t top_label = top_label_map[top_ccl[loc]];

				if (top_label != current_label) {
					std::vector<uint16_t>& label_points = ptcs[current_label];
					label_points.push_back(x);
					label_points.push_back(y);
					label_points.push_back(static_cast<uint16_t>(z));
				}
			}
		}
	}

	return ptcs;
}

template <typename LABEL>
std::unordered_map<uint64_t, std::vector<uint16_t>>
point_cloud(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const std::optional<std::vector<uint64_t>> labels = std::nullopt,
	const bool skip_background = false,
	int parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, std::vector<uint16_t>>();
	}

	const int64_t sz = header.sz;
	const int64_t szr = get_szr(header, z_start, z_end);
	z_start = std::max(std::min(z_start, static_cast<int64_t>(sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const bool selective = labels.has_value();
	robin_hood::unordered_flat_set<uint64_t> labels_set;
	if (selective) {
		labels_set.insert(labels->begin(), labels->end());
	}

	std::span<const unsigned char> binary(buffer, num_bytes);
	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, std::vector<uint16_t>> ptc;
	ptc.reserve(num_labels / (header.sz / szr));

	const bool binary_image_optimization = (
		num_labels == 2
		&& header.is_remapped == false
	);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, (int)(szr / 2));
	parallel = std::max(parallel, (int)1);

	int64_t cz = szr / parallel;
	cz = std::max(cz, (int64_t)2);
	cz = std::min(cz, (int64_t)szr);

	const int64_t num_chunks = std::max((int64_t)((szr + cz - 1) / cz), (int64_t)1);

	ThreadPool pool(parallel);
	std::vector<pointcloud_t> ptcs(num_chunks);

	for (int64_t z = z_start, i = 0; z < (int64_t)z_end; z += cz, i++) {
		pool.enqueue([&,i,z](size_t t) {
			int64_t z_start_i = (z == z_start) ? z : z - 1;
			int64_t z_end_i = std::min(z + cz, z_end);

			if (binary_image_optimization) {
				ptcs[i] = _point_cloud_binary_image<LABEL>(
					buffer, num_bytes,
					z_start_i, z_end_i,
					/*label=*/( 
						(labels.has_value() && labels->size()) > 0
							? std::optional<uint64_t>{(*labels)[0]}
							: std::nullopt
					),
					skip_background,
					/*skip_last_2d=*/(z_end_i < z_end)
				);
			}
			else {
				ptcs[i] = _point_cloud_multivalued_image<LABEL>(
					buffer, num_bytes,
					z_start_i, z_end_i,
					labels,
					skip_background,
					/*skip_last_2d=*/(z_end_i < z_end)
				);
			}
		});
	}

	pool.join();

	size_t size = 0;
	for (int64_t i = 0; i < num_chunks; i++) {
		size += ptcs[i].size();
	}
	
	pointcloud_t accumulated_ptc = std::move(ptcs[0]);
	accumulated_ptc.reserve(size);

	for (int64_t i = 1; i < num_chunks; i++) {
		for (auto& [label, pts] : ptcs[i]) {
			auto& vec = accumulated_ptc[label];
			vec.insert(vec.end(), pts.begin(), pts.end());
		}
	}

	return accumulated_ptc;
}

auto point_cloud(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const std::optional<std::vector<uint64_t>> labels = std::nullopt,
	const bool skip_background = false,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return point_cloud<uint8_t>(
			buffer, num_bytes,
			z_start, z_end,
			labels, skip_background, parallel
		);
	}
	else if (header.data_width == 2) {
		return point_cloud<uint16_t>(
			buffer, num_bytes,
			z_start, z_end,
			labels, skip_background, parallel
		);
	}
	else if (header.data_width == 4) {
		return point_cloud<uint32_t>(
			buffer, num_bytes,
			z_start, z_end,
			labels, skip_background, parallel
		);
	}
	else {
		return point_cloud<uint64_t>(
			buffer, num_bytes,
			z_start, z_end,
			labels, skip_background, parallel
		);
	}
}

auto point_cloud(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	const std::optional<std::vector<uint64_t>> labels = std::nullopt,
	const bool skip_background = false,
	size_t parallel = 1
) {
	return point_cloud(
		buffer.data(),
		buffer.size(),
		z_start, z_end, 
		labels, skip_background, parallel
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
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);
	const int64_t szr = get_szr(header, z_start, z_end);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, uint64_t>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);
	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, uint64_t> cts;
	cts.reserve(num_labels / (header.sz / szr));

	const uint64_t sxy = header.sx * header.sy;

	std::mutex mtx;

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			std::vector<uint64_t> subcounts(N);

			for (uint64_t i = 0; i < sxy; i++) {
				subcounts[ccl[i]]++;
			}

			std::unique_lock<std::mutex> lock(mtx);
			for (uint64_t i = 0; i < N; i++) {
				cts[label_map[i]] += subcounts[i];
			}
		}
	);

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

template <typename LABEL>
std::unordered_map<uint64_t, std::array<double, 3>>
centroids(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);
	const int64_t szr = get_szr(header, z_start, z_end);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, std::array<double,3>>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);
	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, std::array<uint64_t, 4>> cts; // label: x,y,z,N
	cts.reserve(num_labels / (header.sz / szr));

	std::mutex mtx;

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			std::vector<std::array<uint64_t, 4>> subcounts(N);

			for (uint64_t y = 0; y < header.sy; y++) {
				for (uint64_t x = 0; x < header.sx; x++) {
					uint32_t ccl_label = ccl[x + header.sx * y];

					subcounts[ccl_label][0] += x;
					subcounts[ccl_label][1] += y;
					subcounts[ccl_label][2] += z;
					subcounts[ccl_label][3]++;
				}
			}

			std::unique_lock<std::mutex> lock(mtx);
			for (uint64_t i = 0; i < N; i++) {
				auto& arr = cts[label_map[i]];
				arr[0] += subcounts[i][0];
				arr[1] += subcounts[i][1];
				arr[2] += subcounts[i][2];
				arr[3] += subcounts[i][3];
			}
		}
	);

	std::unordered_map<uint64_t, std::array<double, 3>> ctsf;
	ctsf.reserve(cts.size());

	for (auto& [label, ctarr] : cts) {
		ctsf[label] = {
			static_cast<double>(ctarr[0]) / static_cast<double>(ctarr[3]),
			static_cast<double>(ctarr[1]) / static_cast<double>(ctarr[3]),
			static_cast<double>(ctarr[2]) / static_cast<double>(ctarr[3]),
		};
	}

	return ctsf;
}

auto centroids(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return centroids<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 2) {
		return centroids<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 4) {
		return centroids<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else {
		return centroids<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
}

auto centroids(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1
) {
	return centroids(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}

template <typename LABEL>
std::unordered_map<uint64_t, std::array<uint32_t, 6>>
bounding_boxes(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {

	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);

	if (voxels == 0) {
		return std::unordered_map<uint64_t, std::array<uint32_t,6>>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);
	uint64_t num_labels = crackle::labels::num_labels(binary);
	std::unordered_map<uint64_t, std::array<uint32_t, 6>> bbxes; // label: xmin,ymin,zmin,xmax,ymax,zmax
	bbxes.reserve(num_labels);

	std::vector<uint64_t> uniq = crackle::labels::unique(binary);
	for (auto lbl : uniq) {
		auto& bbx = bbxes[lbl];
		bbx[0] = std::numeric_limits<uint32_t>::max();
		bbx[1] = std::numeric_limits<uint32_t>::max();
		bbx[2] = std::numeric_limits<uint32_t>::max();
	}

	std::mutex mtx;

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			std::vector<std::array<uint32_t, 4>> slice_bbxs(N);

			for (uint64_t i = 0; i < N; i++) {
				slice_bbxs[i][0] = std::numeric_limits<uint32_t>::max(); // xmin
				slice_bbxs[i][1] = std::numeric_limits<uint32_t>::max(); // ymin
			}

			for (uint64_t y = 0; y < header.sy; y++) {
				for (uint64_t x = 0; x < header.sx; x++) {
					uint32_t ccl_label = ccl[x + header.sx * y];

					auto& bbx = slice_bbxs[ccl_label];

					bbx[0] = std::min(bbx[0], static_cast<uint32_t>(x));
					bbx[1] = std::min(bbx[1], static_cast<uint32_t>(y));
					bbx[2] = std::max(bbx[2], static_cast<uint32_t>(x));
					bbx[3] = std::max(bbx[3], static_cast<uint32_t>(y));
				}
			}

			std::unique_lock<std::mutex> lock(mtx);
			for (uint64_t i = 0; i < N; i++) {
				auto& bbx = bbxes[label_map[i]];
				bbx[0] = std::min(bbx[0], slice_bbxs[i][0]);
				bbx[1] = std::min(bbx[1], slice_bbxs[i][1]);
				bbx[2] = std::min(bbx[2], static_cast<uint32_t>(z));
				bbx[3] = std::max(bbx[3], slice_bbxs[i][2]);
				bbx[4] = std::max(bbx[4], slice_bbxs[i][3]);
				bbx[5] = std::max(bbx[5], static_cast<uint32_t>(z));
			}
		}
	);

	return bbxes;
}

auto bounding_boxes(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return bounding_boxes<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 2) {
		return bounding_boxes<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 4) {
		return bounding_boxes<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else {
		return bounding_boxes<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
}

auto bounding_boxes(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1
) {
	return bounding_boxes(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}

uint8_t* voxel_connectivity_graph(
	const unsigned char* buffer,
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1,
	const int connectivity = 4
) {
	if (connectivity != 4 && connectivity != 6) {
		std::string err = "crackle: voxel_connectivity_graph: only connectivity 4 and 6 are currently supported.";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	const uint64_t sz = header.sz;
	const uint64_t sxy = header.sx * header.sy;

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	z_start = std::max(std::min(z_start, static_cast<int64_t>(sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const int64_t szr = z_end - z_start;

	const uint64_t voxels = (
		static_cast<uint64_t>(sx) 
		* static_cast<uint64_t>(sy) 
		* static_cast<uint64_t>(szr)
	);

	uint8_t* vcg = new uint8_t[voxels];

	if (voxels == 0) {
		return vcg;
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(szr));

	ThreadPool pool(parallel);

	for (uint64_t i = 0; i < crack_codes.size(); i++) {
		pool.enqueue([&,i](size_t t){
			auto& crack_code = crack_codes[i];

			crack_code_to_vcg(
				/*code=*/crack_code,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model,
				/*slice_edges=*/(vcg + i * sxy)
			);
		});
	}

	pool.join();

	if (sz == 1 || connectivity == 4) {
		return vcg;
	}

	std::vector<std::vector<uint32_t>> ccls(2);
	std::vector<std::vector<uint64_t>> label_maps(2);
	for (size_t t = 0; t < ccls.size(); t++) {
		ccls[t].resize(sxy);
	}

	std::span<const uint8_t> init_sub_vcg(vcg, vcg + sxy);

	uint64_t N = 0;
	crackle::cc3d::color_connectivity_graph<uint32_t>(
		init_sub_vcg, sx, sy, 1, ccls[0].data(), N
	);

	label_maps[0] = decode_label_map<uint64_t>(
		header, binary, ccls[0].data(), N, z_start, z_start+1
	);

	for (size_t z = z_start + 1, i = 1; i < crack_codes.size(); z++, i++) {
		size_t bottom_j = i & 0b1;
		size_t top_j = (i-1) & 0b1;

		std::vector<uint32_t>& top_ccl = ccls[top_j];
		std::vector<uint32_t>& bottom_ccl = ccls[bottom_j];

		const std::span<uint8_t> sub_vcg(
			vcg + i * sxy, 
			vcg + (i+1) * sxy
		);

		uint64_t N = 0;
		crackle::cc3d::color_connectivity_graph<uint32_t>(
			sub_vcg, sx, sy, 1, bottom_ccl.data(), N
		);

		label_maps[bottom_j] = decode_label_map<uint64_t>(
			header, binary, bottom_ccl.data(), N, z, z+1
		);

		const std::vector<uint64_t>& top_labels = label_maps[top_j];
		const std::vector<uint64_t>& bottom_labels = label_maps[bottom_j];

		for (size_t y = 0; y < sy; y++) {
			for (size_t x = 0; x < sx; x++) {
				const int64_t loc = x + sx * y;
				const uint64_t bottom_label = bottom_labels[bottom_ccl[loc]];
				const uint64_t top_label = top_labels[top_ccl[loc]];

				if (bottom_label == top_label) {
					vcg[loc + (z-1) * sxy] |= 0b010000;
					vcg[loc + z * sxy] |= 0b100000;
				}
			}
		}
	}

	// set boundary to passable for compatibility with
	// certain operations like autapse elimination
	for (size_t y = 0; y < sy; y++) {
		for (size_t x = 0; x < sx; x++) {
			const int64_t loc = x + sx * y;
			vcg[loc] |= 0b100000;
			vcg[loc + (sz-1) * sxy] |= 0b010000;
		}
	}

	return vcg;
}

auto voxel_connectivity_graph(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1,
	const int connectivity = 4
) {
	return voxel_connectivity_graph(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel,
		connectivity
	);
}

struct pair_hash {
	inline std::size_t operator()(const std::pair<uint64_t,uint64_t> & v) const {
		return v.first * 31 + v.second; // arbitrary hash fn
	}
};

typedef std::unordered_map<std::pair<uint64_t,uint64_t>, float, pair_hash> contacts_t;

// single threaded version
contacts_t _contacts(
	const unsigned char* buffer,
	const size_t num_bytes,
	int64_t z_start,
	int64_t z_end,
	const float wx = 1.0,
	const float wy = 1.0,
	const float wz = 1.0,
	const bool skip_last_2d = false
) {
	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	const uint64_t sz = header.sz;
	const uint64_t sxy = header.sx * header.sy;

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	z_start = std::max(std::min(z_start, static_cast<int64_t>(sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const int64_t szr = z_end - z_start;

	const uint64_t voxels = (
		static_cast<uint64_t>(sx) 
		* static_cast<uint64_t>(sy) 
		* static_cast<uint64_t>(szr)
	);

	contacts_t edges;

	if (voxels == 0) {
		return edges;
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	std::vector<uint8_t> vcg(sxy);
	std::vector<std::vector<uint32_t>> ccls(2);
	std::vector<std::vector<uint64_t>> label_maps(2);
	for (size_t t = 0; t < ccls.size(); t++) {
		ccls[t].resize(sxy);
	}

	crack_code_to_vcg(
		/*code=*/crack_codes[0],
		/*sx=*/sx, /*sy=*/sy,
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
		/*markov_model=*/markov_model,
		/*slice_edges=*/vcg.data()
	);

	uint64_t N = 0;
	crackle::cc3d::color_connectivity_graph<uint32_t>(
		vcg, sx, sy, 1, ccls[0].data(), N
	);

	label_maps[0] = decode_label_map<uint64_t>(
		header, binary, ccls[0].data(), N, z_start, z_start+1
	);

	const float area_x = wy * wz;
	const float area_y = wx * wz;
	const float area_z = wx * wy;

	auto add_edge = [&edges](uint64_t a, uint64_t b, float area) {
		if (a != 0 && b != 0) {
			const std::pair<uint64_t,uint64_t> p = (a <= b) 
				? std::pair<uint64_t,uint64_t>(a,b)
				: std::pair<uint64_t,uint64_t>(b,a);
				
			edges[p] += area;
		}
	};

	for (size_t y = 0; y < sy; y++) {
		for (size_t x = 0; x < sx; x++) {
			const int64_t loc = x + sx * y;

			if (x < sx - 1 && ccls[0][loc] != ccls[0][loc+1]) {
				add_edge(label_maps[0][ccls[0][loc]], label_maps[0][ccls[0][loc+1]], area_x);
			}
			if (y < sy - 1 && ccls[0][loc] != ccls[0][loc+sx]) {
				add_edge(label_maps[0][ccls[0][loc]], label_maps[0][ccls[0][loc+sx]], area_y);
			}
		}
	}

	for (size_t z = z_start + 1, i = 1; i < crack_codes.size(); z++, i++) {
		size_t bottom_j = i & 0b1;
		size_t top_j = (i-1) & 0b1;

		crack_code_to_vcg(
			/*code=*/crack_codes[i],
			/*sx=*/sx, /*sy=*/sy,
			/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
			/*markov_model=*/markov_model,
			/*slice_edges=*/vcg.data()
		);

		std::vector<uint32_t>& top_ccl = ccls[top_j];
		std::vector<uint32_t>& bottom_ccl = ccls[bottom_j];

		uint64_t N = 0;
		crackle::cc3d::color_connectivity_graph<uint32_t>(
			vcg, sx, sy, 1, bottom_ccl.data(), N
		);

		label_maps[bottom_j] = decode_label_map<uint64_t>(
			header, binary, bottom_ccl.data(), N, z, z+1
		);

		const std::vector<uint64_t>& top_labels = label_maps[top_j];
		const std::vector<uint64_t>& bottom_labels = label_maps[bottom_j];

		for (size_t y = 0; y < sy; y++) {
			for (size_t x = 0; x < sx; x++) {
				const int64_t loc = x + sx * y;

				if (!skip_last_2d || i < crack_codes.size() - 1) {
					if (x < sx - 1 && bottom_ccl[loc] != bottom_ccl[loc+1]) {
						add_edge(
							bottom_labels[bottom_ccl[loc]], 
							bottom_labels[bottom_ccl[loc+1]],
							area_x
						);
					}
					if (y < sy - 1 && bottom_ccl[loc] != bottom_ccl[loc+sx]) {
						add_edge(
							bottom_labels[bottom_ccl[loc]],
							bottom_labels[bottom_ccl[loc+sx]],
							area_y
						);
					}
				}
				if (bottom_labels[bottom_ccl[loc]] != top_labels[top_ccl[loc]]) {
					add_edge(
						bottom_labels[bottom_ccl[loc]],
						top_labels[top_ccl[loc]],
						area_z
					);
				}
			}
		}
	}

	return edges;
}

auto contacts(
	const unsigned char* buffer,
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const float wx = 1,
	const float wy = 1,
	const float wz = 1,
	int parallel = 1
) {
	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	const uint64_t sz = header.sz;
	const uint64_t sxy = header.sx * header.sy;

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	z_start = std::max(std::min(z_start, static_cast<int64_t>(sz - 1)), static_cast<int64_t>(0));
	z_end = z_end < 0 ? static_cast<int64_t>(sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(sz)), static_cast<int64_t>(0));

	if (z_start >= z_end) {
		std::string err = "crackle: Invalid range: ";
		err += std::to_string(z_start);
		err += std::string(" - ");
		err += std::to_string(z_end);
		throw std::runtime_error(err);
	}

	const uint64_t szr = z_end - z_start;

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, (int)szr / 2);
	parallel = std::max(parallel, (int)1);

	int64_t cz = szr / parallel;
	cz = std::max(cz, (int64_t)2);
	cz = std::min(cz, (int64_t)szr);

	const uint64_t num_chunks = std::max((szr + cz - 1) / cz, (uint64_t)1);

	ThreadPool pool(parallel);
	std::vector<contacts_t> edges(num_chunks);

	for (int64_t z = z_start, i = 0; z < (int64_t)z_end; z += cz, i++) {
		pool.enqueue([&,i,z](size_t t) {
			int64_t z_start_i = (z == z_start) ? z : z - 1;
			int64_t z_end_i = std::min(z + cz, z_end);
			edges[i] = _contacts(
				buffer, num_bytes,
				z_start_i, z_end_i,
				wx, wy, wz,
				/*skip_last_2d=*/(z_end_i < z_end)
			);
		});
	}

	pool.join();

	contacts_t accumulated_edges = std::move(edges[0]);
	for (int64_t i = 1; i < num_chunks; i++) {
		for (auto& [edge, area] : edges[i]) {
			accumulated_edges[edge] += area;
		}
	}

	return accumulated_edges;
}

auto contacts(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	const float wx = 1,
	const float wy = 1,
	const float wz = 1,
	int parallel = 1
) {
	return contacts(
		buffer.data(),
		buffer.size(),
		z_start, z_end,
		wx, wy, wz,
		parallel
	);
}

bool array_equal(
	const unsigned char* buffer1,
	const size_t num_bytes1,
	const unsigned char* buffer2,
	const size_t num_bytes2,
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);
	const CrackleHeader header2 = get_header(buffer2, num_bytes2);

	const uint64_t voxels1 = get_voxels(header1, 0, -1);
	const uint64_t voxels2 = get_voxels(header2, 0, -1);

	if (voxels1 == 0 || voxels2 == 0) {
		return voxels1 == voxels2;
	}

	if (
		header1.sx != header2.sx 
		|| header1.sy != header2.sy 
		|| header1.sz != header2.sz
	) {
		return false;
	}

	const size_t sx = header1.sx;
	const size_t sy = header1.sy;
	const size_t sz = header1.sz;

	std::span<const unsigned char> binary1(buffer1, num_bytes1);
	std::span<const unsigned char> binary2(buffer2, num_bytes2);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model1 = decode_markov_model(header1, binary1);
	std::vector<std::vector<uint8_t>> markov_model2 = decode_markov_model(header2, binary2);
	
	auto crack_codes1 = get_crack_codes(header1, binary1, 0, sz);
	auto crack_codes2 = get_crack_codes(header2, binary2, 0, sz);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, sz);

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel * 2);
	std::vector<std::vector<uint32_t>> ccls(parallel * 2);

	for (size_t t = 0; t < parallel * 2; t++) {
		vcgs[t].resize(header1.sx * header1.sy);
		ccls[t].resize(header1.sx * header1.sy);
	}
	
	std::atomic<bool> is_equal{true};

	for (size_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			if (!is_equal) {
				return;
			}

			t *= 2;

			std::vector<uint8_t>& vcg1 = vcgs[t];
			std::vector<uint8_t>& vcg2 = vcgs[t+1];
			
			std::vector<uint32_t>& ccl1 = ccls[t];
			std::vector<uint32_t>& ccl2 = ccls[t+1];

			auto& crack_code1 = crack_codes1[z];
			auto& crack_code2 = crack_codes2[z];

			crack_code_to_vcg(
				/*code=*/crack_code1,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header1.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model1,
				/*slice_edges=*/vcg1.data()
			);

			if (!is_equal) {
				return;
			}

			crack_code_to_vcg(
				/*code=*/crack_code2,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header2.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model2,
				/*slice_edges=*/vcg2.data()
			);

			if (!is_equal) {
				return;
			}

			uint64_t N1 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg1, sx, sy, 1, ccl1.data(), N1
			);

			if (!is_equal) {
				return;
			}

			uint64_t N2 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg2, sx, sy, 1, ccl2.data(), N2
			);

			if (!is_equal) {
				return;
			}

			if (N1 != N2) {
    			is_equal.store(false, std::memory_order_relaxed);
				return;
			}

			std::vector<uint64_t> label_map1 = decode_label_map<uint64_t>(
				header1, binary1, ccl1.data(), N1, z, z+1
			);
			std::vector<uint64_t> label_map2 = decode_label_map<uint64_t>(
				header2, binary2, ccl2.data(), N2, z, z+1
			);

			for (size_t y = 0; y < sy; y++) {
				for (size_t x = 0; x < sx; x++) {
					const int64_t loc = x + sx * y;
					const uint64_t a = label_map1[ccl1[loc]];
					const uint64_t b = label_map2[ccl2[loc]];

					if (a != b) {
						is_equal.store(false, std::memory_order_relaxed);
						return;
					}
				}
			}
		});
	}

	pool.join();

	return is_equal;
}

bool array_equal(
	const std::span<const unsigned char>& buffer1,
	const std::span<const unsigned char>& buffer2,
	const int parallel = 1
) {
	return array_equal(
		buffer1.data(),
		buffer1.size(),
		buffer2.data(),
		buffer2.size(),
		parallel
	);
}

std::unordered_map<uint64_t,uint64_t>
component_map(
	const unsigned char* buffer1,
	const size_t num_bytes1,
	const unsigned char* buffer2,
	const size_t num_bytes2,
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);
	const CrackleHeader header2 = get_header(buffer2, num_bytes2);

	const uint64_t voxels1 = get_voxels(header1, 0, -1);
	const uint64_t voxels2 = get_voxels(header2, 0, -1);

	if (
		header1.sx != header2.sx 
		|| header1.sy != header2.sy 
		|| header1.sz != header2.sz
	) {
		throw std::invalid_argument("The shapes of the input arrays must match.");
	}

	if (voxels1 == 0 || voxels2 == 0) {
		return {};
	}

	const size_t sx = header1.sx;
	const size_t sy = header1.sy;
	const size_t sz = header1.sz;

	std::span<const unsigned char> binary1(buffer1, num_bytes1);
	std::span<const unsigned char> binary2(buffer2, num_bytes2);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model1 = decode_markov_model(header1, binary1);
	std::vector<std::vector<uint8_t>> markov_model2 = decode_markov_model(header2, binary2);
	
	auto crack_codes1 = get_crack_codes(header1, binary1, 0, sz);
	auto crack_codes2 = get_crack_codes(header2, binary2, 0, sz);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, sz);

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel * 2);
	std::vector<std::vector<uint32_t>> ccls(parallel * 2);

	for (size_t t = 0; t < parallel * 2; t++) {
		vcgs[t].resize(header1.sx * header1.sy);
		ccls[t].resize(header1.sx * header1.sy);
	}

	std::mutex mtx;
	
	const uint64_t num_component_labels = crackle::labels::num_labels(binary1);

	std::unordered_map<uint64_t, uint64_t> mapping;
	mapping.reserve(num_component_labels);

	for (size_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			t *= 2;

			std::vector<uint8_t>& vcg1 = vcgs[t];
			std::vector<uint8_t>& vcg2 = vcgs[t+1];
			
			std::vector<uint32_t>& ccl1 = ccls[t];
			std::vector<uint32_t>& ccl2 = ccls[t+1];

			auto& crack_code1 = crack_codes1[z];
			auto& crack_code2 = crack_codes2[z];

			crack_code_to_vcg(
				/*code=*/crack_code1,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header1.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model1,
				/*slice_edges=*/vcg1.data()
			);

			crack_code_to_vcg(
				/*code=*/crack_code2,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header2.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model2,
				/*slice_edges=*/vcg2.data()
			);

			uint64_t N1 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg1, sx, sy, 1, ccl1.data(), N1
			);

			uint64_t N2 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg2, sx, sy, 1, ccl2.data(), N2
			);

			std::vector<uint64_t> label_map1 = decode_label_map<uint64_t>(
				header1, binary1, ccl1.data(), N1, z, z+1
			);
			std::vector<uint64_t> label_map2 = decode_label_map<uint64_t>(
				header2, binary2, ccl2.data(), N2, z, z+1
			);

			std::unordered_map<uint64_t, uint64_t> submap;

			uint64_t last_label = label_map1[ccl1[0]];
			submap[last_label] = label_map2[ccl2[0]];
			for (size_t y = 0; y < sy; y++) {
				for (size_t x = 0; x < sx; x++) {
					const int64_t loc = x + sx * y;
					const uint64_t child_label = label_map1[ccl1[loc]];
					const uint64_t parent_label = label_map2[ccl2[loc]];

					if (last_label == child_label) {
						continue;
					}
					else {
						submap[child_label] = parent_label;
						last_label = child_label;
					}
				}
			}

			std::unique_lock<std::mutex> lock(mtx);
			mapping.merge(submap);
		});
	}

	pool.join();

	return mapping;
}

auto component_map(
	const std::span<const unsigned char>& buffer1,
	const std::span<const unsigned char>& buffer2,
	const int parallel = 1
) {
	return component_map(
		buffer1.data(),
		buffer1.size(),
		buffer2.data(),
		buffer2.size(),
		parallel
	);
}

template <typename LABEL>
std::vector<std::vector<unsigned char>>
mode_pooling_2x2x1(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const bool sparse = false,
	size_t parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);
	const int64_t szr = get_szr(header, z_start, z_end);

	if (voxels == 0) {
		return std::vector<std::vector<unsigned char>>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	std::vector<std::vector<unsigned char>> ds_binaries(szr);

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	
	const uint64_t osx = ((sx + 1) >> 1);
	const uint64_t osy = ((sy + 1) >> 1);
	const uint64_t osxy = osx * osy;

	const size_t xodd = (sx & 0x01);
	const size_t yodd = (sy & 0x01);

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			std::vector<LABEL> oimg(osxy);

			auto idx = [&](size_t x, size_t y) {
				return x + sx * y;
			};

			auto oidx = [&](size_t ox, size_t oy) {
				return ox + osx * oy;
			};

			uint64_t oy = 0;

			for (uint64_t y = 0; y < sy - yodd; y += 2) {
			  uint64_t ox = 0;

			  for (uint64_t x = 0; x < sx - xodd; x += 2) {
				const LABEL a = label_map[ccl[idx(x,y)]];
				const LABEL b = label_map[ccl[idx(x+1,y)]];
				const LABEL c = label_map[ccl[idx(x,y+1)]];
				const LABEL d = label_map[ccl[idx(x+1,y+1)]];

				uint64_t out = oidx(ox, oy);

				if (sparse) {
					if (a > 0 && a == b) {
						oimg[out] = a;
					}
					else if (a > 0 && a == c) {
						oimg[out] = a;
					}
					else if (b > 0 && b == c) {
						oimg[out] = b;
					}
					else {
						oimg[out] = d ? d : (a ? a : (b ? b : c));
					}					
				}
				else {
					if (a == b) {
						oimg[out] = a;
					}
					else if (a == c) {
						oimg[out] = a;
					}
					else if (b == c) {
						oimg[out] = b;
					}
					else {
						oimg[out] = d;
					}
				}

				ox++;
			  }
			  if (xodd) {
			    uint64_t out = oidx(osx - 1, oy);
			    oimg[out] = label_map[ccl[idx(sx - 1, y)]];
			  }
			  oy++;
			}

			if (yodd) {
			  for (uint64_t x = 0; x < osx - xodd; x++) {
			    oimg[oidx(x, osy - 1)] = label_map[ccl[idx(2*x, sy - 1)]];
			  }
			  if (xodd) {
			    oimg[oidx(osx - 1, osy - 1)] = label_map[ccl[idx(sx - 1, sy - 1)]];
			  }
			}

			ds_binaries[z-z_start] = crackle::compress<LABEL>(
				oimg.data(),
				osx, osy, 1,
				/*allow_pins=*/false,
				/*fortran_order=*/true,
				/*markov_model_order=*/0,
				/*optimize_pins=*/false,
				/*auto_bgcolor=*/false,
				/*manual_bgcolor=*/0,
				/*parallel=*/1,
				/*crack_format=*/crackle::CrackFormat::IMPERMISSIBLE
			);
		}
	);

	return ds_binaries;
}

auto mode_pooling_2x2x1(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	const bool sparse = false,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return mode_pooling_2x2x1<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, sparse, parallel
		);
	}
	else if (header.data_width == 2) {
		return mode_pooling_2x2x1<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, sparse, parallel
		);
	}
	else if (header.data_width == 4) {
		return mode_pooling_2x2x1<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, sparse, parallel
		);
	}
	else {
		return mode_pooling_2x2x1<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, sparse, parallel
		);
	}
}

auto mode_pooling_2x2x1(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	const bool sparse = false,
	size_t parallel = 1
) {
	return mode_pooling_2x2x1(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}

template <typename LABEL>
robin_hood::unordered_flat_set<uint64_t>
edge_labels(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);
	const int64_t szr = get_szr(header, z_start, z_end);

	if (voxels == 0) {
		return {};
	}

	std::span<const unsigned char> binary(buffer, num_bytes);
	uint64_t num_labels = crackle::labels::num_labels(binary);
	robin_hood::unordered_flat_set<uint64_t> boundary_labels;
	boundary_labels.reserve(num_labels / (header.sz / szr));

	const int64_t sx = (int64_t)header.sx;
	const int64_t sy = (int64_t)header.sy;
	const int64_t sz = (int64_t)header.sz;

	std::mutex mtx;

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			robin_hood::unordered_flat_set<uint64_t> boundary_labels_subset;

			if (z == 0 || z == sz - 1) {
				uint64_t last = ccl[0];
				boundary_labels_subset.emplace(label_map[last]);
				for (int64_t i = 1; i < sx * sy; i++) {
					if (ccl[i] != last) {
						boundary_labels_subset.emplace(label_map[ccl[i]]);
						last = ccl[i];
					}
				}
			}
			else {
				uint64_t last_a = ccl[0];
				uint64_t last_b = ccl[sx * (sy-1)];
				boundary_labels_subset.emplace(label_map[last_a]);
				boundary_labels_subset.emplace(label_map[last_b]);
				for (int64_t x = 0; x < sx; x++) {
					if (ccl[x] != last_a) {
						boundary_labels_subset.emplace(label_map[ccl[x]]);
						last_a = ccl[x];
					}
					if (ccl[x + sx * (sy-1)] != last_b) {
						boundary_labels_subset.emplace(label_map[ccl[x + sx * (sy-1)]]);
						last_b = ccl[x + sx * (sy-1)];
					}
				}
				last_a = ccl[0];
				last_b = ccl[(sx-1)];
				for (int64_t y = 1; y < sy - 1; y++) {
					if (ccl[sx * y] != last_a) {
						boundary_labels_subset.emplace(label_map[ccl[sx * y]]);
						last_a = ccl[sx * y];
					}
					if (ccl[(sx - 1) + sx * y] != last_b) {
						boundary_labels_subset.emplace(label_map[ccl[(sx-1) + sx * y]]);
						last_b = ccl[(sx - 1) + sx * y];
					}
				}
			}

			std::unique_lock<std::mutex> lock(mtx);
			boundary_labels.insert(
				boundary_labels_subset.begin(), boundary_labels_subset.end()
			);
		}
	);

	return boundary_labels;
}

auto edge_labels(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return edge_labels<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 2) {
		return edge_labels<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 4) {
		return edge_labels<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else {
		return edge_labels<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
}

auto edge_labels(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1
) {
	return edge_labels(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}


template <typename LABEL>
std::vector<std::vector<unsigned char>>
apply_binary_mask(
	const unsigned char* buffer1, // input image
	const size_t num_bytes1,
	const unsigned char* buffer2,
	const size_t num_bytes2, // mask buffer
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);
	const CrackleHeader header2 = get_header(buffer2, num_bytes2);

	const uint64_t voxels1 = get_voxels(header1, 0, -1);
	const uint64_t voxels2 = get_voxels(header2, 0, -1);

	if (voxels1 == 0 || voxels2 == 0) {
		return std::vector<std::vector<unsigned char>>();
	}

	if (
		header1.sx != header2.sx 
		|| header1.sy != header2.sy 
		|| header1.sz != header2.sz
	) {
		return std::vector<std::vector<unsigned char>>();
	}

	const size_t sx = header1.sx;
	const size_t sy = header1.sy;
	const size_t sz = header1.sz;

	std::span<const unsigned char> binary1(buffer1, num_bytes1);
	std::span<const unsigned char> binary2(buffer2, num_bytes2);

	std::vector<std::vector<unsigned char>> masked_binaries(sz);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model1 = decode_markov_model(header1, binary1);
	std::vector<std::vector<uint8_t>> markov_model2 = decode_markov_model(header2, binary2);
	
	auto crack_codes1 = get_crack_codes(header1, binary1, 0, sz);
	auto crack_codes2 = get_crack_codes(header2, binary2, 0, sz);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, sz);

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel * 2);
	std::vector<std::vector<uint32_t>> ccls(parallel * 2);
	std::vector<std::vector<LABEL>> new_slices(parallel);

	for (size_t t = 0; t < parallel * 2; t++) {
		vcgs[t].resize(header1.sx * header1.sy);
		ccls[t].resize(header1.sx * header1.sy);
	}

	for (size_t t = 0; t < parallel; t++) {
		new_slices[t].resize(header1.sx * header1.sy);
	}
	
	for (size_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			const size_t tt = t * 2;

			std::vector<uint8_t>& vcg1 = vcgs[tt];
			std::vector<uint8_t>& vcg2 = vcgs[tt+1];
			
			std::vector<uint32_t>& ccl1 = ccls[tt];
			std::vector<uint32_t>& ccl2 = ccls[tt+1];

			std::vector<LABEL>& ablated = new_slices[t];
			std::fill(ablated.begin(), ablated.end(), static_cast<LABEL>(0));

			auto& crack_code1 = crack_codes1[z];
			auto& crack_code2 = crack_codes2[z];

			crack_code_to_vcg(
				/*code=*/crack_code1,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header1.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model1,
				/*slice_edges=*/vcg1.data()
			);

			crack_code_to_vcg(
				/*code=*/crack_code2,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header2.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model2,
				/*slice_edges=*/vcg2.data()
			);

			uint64_t N1 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg1, sx, sy, 1, ccl1.data(), N1
			);

			uint64_t N2 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg2, sx, sy, 1, ccl2.data(), N2
			);

			std::vector<uint64_t> label_map1 = decode_label_map<uint64_t>(
				header1, binary1, ccl1.data(), N1, z, z+1
			);
			std::vector<uint64_t> label_map2 = decode_label_map<uint64_t>(
				header2, binary2, ccl2.data(), N2, z, z+1
			);

			for (size_t y = 0; y < sy; y++) {
				for (size_t x = 0; x < sx; x++) {
					const int64_t loc = x + sx * y;
					const uint64_t mask_value = label_map2[ccl2[loc]];

					if (mask_value) {
						ablated[loc] = label_map1[ccl1[loc]];
					}
				}
			}

			masked_binaries[z] = crackle::compress<LABEL>(
				ablated.data(),
				sx, sy, 1,
				/*allow_pins=*/false,
				/*fortran_order=*/true,
				/*markov_model_order=*/0,
				/*optimize_pins=*/false,
				/*auto_bgcolor=*/false,
				/*manual_bgcolor=*/0,
				/*parallel=*/1,
				/*crack_format=*/crackle::CrackFormat::IMPERMISSIBLE
			);
		});
	}

	pool.join();

	return masked_binaries;
}

auto apply_binary_mask(
	const unsigned char* buffer1, // input image
	const size_t num_bytes1,
	const unsigned char* buffer2,
	const size_t num_bytes2, // mask buffer
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);

	if (header1.data_width == 1) {
		return apply_binary_mask<uint8_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2,
			parallel
		);
	}
	else if (header1.data_width == 2) {
		return apply_binary_mask<uint16_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2,
			parallel
		);
	}
	else if (header1.data_width == 4) {
		return apply_binary_mask<uint32_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2,
			parallel
		);
	}
	else {
		return apply_binary_mask<uint64_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2,
			parallel
		);
	}
}

auto apply_binary_mask(
	const std::span<const unsigned char>& buffer1,
	const std::span<const unsigned char>& buffer2,
	size_t parallel = 1
) {
	return apply_binary_mask(
		buffer1.data(), buffer1.size(),
		buffer2.data(), buffer2.size(),
		parallel
	);
}

template <typename LABEL>
std::vector<std::vector<unsigned char>>
recompress(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);
	const uint64_t voxels = get_voxels(header, z_start, z_end);
	const int64_t szr = get_szr(header, z_start, z_end);

	if (voxels == 0) {
		return std::vector<std::vector<unsigned char>>();
	}

	std::vector<std::vector<unsigned char>> recompressed_binaries(szr);

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	const uint64_t sxy = sx * sy;

	std::span<const unsigned char> binary(buffer, num_bytes);

	for_each_z_parallel<LABEL>(
		buffer, num_bytes,
		z_start, z_end, parallel,
		[&](
			std::vector<uint8_t>& vcg, 
			std::vector<uint32_t>& ccl, 
			uint64_t N,
			std::vector<LABEL>& label_map,
			int64_t z
		){
			std::unique_ptr<LABEL[]> rendered(new LABEL[sx*sy]);
			for (uint64_t i = 0; i < sxy; i++) {
				rendered[i] = label_map[ccl[i]];
			}

			recompressed_binaries[z-z_start] = crackle::compress<LABEL>(
				rendered.get(),
				sx, sy, 1,
				/*allow_pins=*/false,
				/*fortran_order=*/true,
				/*markov_model_order=*/0,
				/*optimize_pins=*/false,
				/*auto_bgcolor=*/false,
				/*manual_bgcolor=*/0,
				/*parallel=*/1,
				/*crack_format=*/crackle::CrackFormat::IMPERMISSIBLE
			);
		}
	);

	return recompressed_binaries;
}

auto recompress(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return recompress<uint8_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 2) {
		return recompress<uint16_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else if (header.data_width == 4) {
		return recompress<uint32_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
	else {
		return recompress<uint64_t>(
			buffer, num_bytes,
			z_start, z_end, parallel
		);
	}
}

auto recompress(
	const std::span<const unsigned char>& buffer,
	const int64_t z_start = -1, 
	const int64_t z_end = -1,
	size_t parallel = 1
) {
	return recompress(
		buffer.data(),
		buffer.size(),
		z_start, z_end, parallel
	);
}

enum BinaryOperator : uint8_t {
	EQ = 0,
	NEQ = 1,
	GT = 2,
	GTE = 3,
	LT = 4,
	LTE = 5,
	ADD = 6,
	SUB = 7,
	MUL = 8,
	DIV = 9,
	MOD = 10,
	XOR = 11,
	LSHIFT = 12,
	RSHIFT = 13,
	BITWISE_OR = 14,
	BITWISE_AND = 15
};

template <typename LABEL>
std::vector<std::vector<unsigned char>>
binary_operator(
	const unsigned char* buffer1, // input image
	const size_t num_bytes1,
	const unsigned char* buffer2,
	const size_t num_bytes2, // mask buffer,
	const BinaryOperator oper,
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);
	const CrackleHeader header2 = get_header(buffer2, num_bytes2);

	const uint64_t voxels1 = get_voxels(header1, 0, -1);
	const uint64_t voxels2 = get_voxels(header2, 0, -1);

	if (voxels1 == 0 || voxels2 == 0) {
		return std::vector<std::vector<unsigned char>>();
	}

	if (
		header1.sx != header2.sx 
		|| header1.sy != header2.sy 
		|| header1.sz != header2.sz
	) {
		return std::vector<std::vector<unsigned char>>();
	}

	const size_t sx = header1.sx;
	const size_t sy = header1.sy;
	const size_t sz = header1.sz;

	std::span<const unsigned char> binary1(buffer1, num_bytes1);
	std::span<const unsigned char> binary2(buffer2, num_bytes2);

	std::vector<std::vector<unsigned char>> result_binaries(sz);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model1 = decode_markov_model(header1, binary1);
	std::vector<std::vector<uint8_t>> markov_model2 = decode_markov_model(header2, binary2);
	
	auto crack_codes1 = get_crack_codes(header1, binary1, 0, sz);
	auto crack_codes2 = get_crack_codes(header2, binary2, 0, sz);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, sz);

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel * 2);
	std::vector<std::vector<uint32_t>> ccls(parallel * 2);
	std::vector<std::vector<LABEL>> new_slices(parallel);

	for (size_t t = 0; t < parallel * 2; t++) {
		vcgs[t].resize(header1.sx * header1.sy);
		ccls[t].resize(header1.sx * header1.sy);
	}

	for (size_t t = 0; t < parallel; t++) {
		new_slices[t].resize(header1.sx * header1.sy);
	}
	
	for (size_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			const size_t tt = t * 2;

			std::vector<uint8_t>& vcg1 = vcgs[tt];
			std::vector<uint8_t>& vcg2 = vcgs[tt+1];
			
			std::vector<uint32_t>& ccl1 = ccls[tt];
			std::vector<uint32_t>& ccl2 = ccls[tt+1];

			std::vector<LABEL>& result_map = new_slices[t];

			auto& crack_code1 = crack_codes1[z];
			auto& crack_code2 = crack_codes2[z];

			crack_code_to_vcg(
				/*code=*/crack_code1,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header1.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model1,
				/*slice_edges=*/vcg1.data()
			);

			crack_code_to_vcg(
				/*code=*/crack_code2,
				/*sx=*/sx, /*sy=*/sy,
				/*permissible=*/(header2.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model2,
				/*slice_edges=*/vcg2.data()
			);

			uint64_t N1 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg1, sx, sy, 1, ccl1.data(), N1
			);

			uint64_t N2 = 0;
			crackle::cc3d::color_connectivity_graph<uint32_t>(
				vcg2, sx, sy, 1, ccl2.data(), N2
			);

			std::vector<uint64_t> label_map1 = decode_label_map<uint64_t>(
				header1, binary1, ccl1.data(), N1, z, z+1
			);
			std::vector<uint64_t> label_map2 = decode_label_map<uint64_t>(
				header2, binary2, ccl2.data(), N2, z, z+1
			);

			for (size_t y = 0; y < sy; y++) {
				for (size_t x = 0; x < sx; x++) {
					const int64_t loc = x + sx * y;
					const uint64_t lvalue = label_map1[ccl1[loc]];
					const uint64_t rvalue = label_map2[ccl2[loc]];

					if (oper == BinaryOperator::EQ) {
						result_map[loc] = lvalue == rvalue;
					}
					else if (oper == BinaryOperator::NEQ) {
						result_map[loc] = lvalue != rvalue;
					}
					else if (oper == BinaryOperator::GT) {
						result_map[loc] = lvalue > rvalue;
					}
					else if (oper == BinaryOperator::GTE) {
						result_map[loc] = lvalue >= rvalue;
					}
					else if (oper == BinaryOperator::LT) {
						result_map[loc] = lvalue < rvalue;
					}
					else if (oper == BinaryOperator::LTE) {
						result_map[loc] = lvalue <= rvalue;
					}
					else if (oper == BinaryOperator::ADD) {
						result_map[loc] = lvalue + rvalue;
					}
					else if (oper == BinaryOperator::SUB) {
						result_map[loc] = lvalue - rvalue;
					}
					else if (oper == BinaryOperator::MUL) {
						result_map[loc] = lvalue * rvalue;
					}
					else if (oper == BinaryOperator::DIV) {
						result_map[loc] = lvalue / rvalue;
					}
					else if (oper == BinaryOperator::MOD) {
						result_map[loc] = lvalue % rvalue;
					}
					else if (oper == BinaryOperator::LSHIFT) {
						result_map[loc] = lvalue << rvalue;
					}
					else if (oper == BinaryOperator::RSHIFT) {
						result_map[loc] = lvalue >> rvalue;
					}
					else if (oper == BinaryOperator::BITWISE_OR) {
						result_map[loc] = lvalue | rvalue;
					}
					else if (oper == BinaryOperator::BITWISE_AND) {
						result_map[loc] = lvalue & rvalue;
					}
					else if (oper == BinaryOperator::XOR) {
						result_map[loc] = lvalue ^ rvalue;
					}
					else {
						throw std::invalid_argument("Invalid operator.");
					}
				}
			}

			result_binaries[z] = crackle::compress<LABEL>(
				result_map.data(),
				sx, sy, 1,
				/*allow_pins=*/false,
				/*fortran_order=*/true,
				/*markov_model_order=*/0,
				/*optimize_pins=*/false,
				/*auto_bgcolor=*/false,
				/*manual_bgcolor=*/0,
				/*parallel=*/1,
				/*crack_format=*/crackle::CrackFormat::IMPERMISSIBLE
			);
		});
	}

	pool.join();

	return result_binaries;
}

auto binary_operator(
	const unsigned char* buffer1, // lvalue buf
	const size_t num_bytes1,
	const unsigned char* buffer2, // rvalue buf
	const size_t num_bytes2, 
	const BinaryOperator oper,
	size_t parallel = 1
) {
	const CrackleHeader header1 = get_header(buffer1, num_bytes1);
	const CrackleHeader header2 = get_header(buffer2, num_bytes2);

	bool is_boolean_operator = (
		oper == BinaryOperator::EQ
		|| oper == BinaryOperator::NEQ
		|| oper == BinaryOperator::GT
		|| oper == BinaryOperator::GTE
		|| oper == BinaryOperator::LT
		|| oper == BinaryOperator::LTE
	);

	const int data_width = std::max(
		header1.data_width,
		header2.data_width
	);

	if (is_boolean_operator || data_width == 1) {
		return binary_operator<uint8_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2, oper,
			parallel
		);
	}
	else if (data_width == 2) {
		return binary_operator<uint16_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2, oper,
			parallel
		);
	}
	else if (data_width == 4) {
		return binary_operator<uint32_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2, oper,
			parallel
		);
	}
	else {
		return binary_operator<uint64_t>(
			buffer1, num_bytes1,
			buffer2, num_bytes2, oper,
			parallel
		);
	}
}

auto binary_operator(
	const std::span<const unsigned char>& buffer1,
	const std::span<const unsigned char>& buffer2,
	const BinaryOperator oper,
	size_t parallel = 1
) {
	return binary_operator(
		buffer1.data(), buffer1.size(),
		buffer2.data(), buffer2.size(),
		oper,
		parallel
	);
}

template <std::integral LABEL, std::integral POINT_T>
std::vector<LABEL> scattered_points(
	const unsigned char* buffer,
	const size_t num_bytes,
	const std::span<const POINT_T> points,
	size_t parallel = 1
) {
	const CrackleHeader header = get_header(buffer, num_bytes);

	std::vector<LABEL> values(points.size() / 3);

	if (header.voxels() == 0) {
		return values;
	}

	if (points.size() == 0) {
		return values;
	}

	if (points.size() % 3 != 0) {
		throw std::invalid_argument("points must contain xyz triples.");
	}

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;
	const uint64_t sz = header.sz;

	std::vector<std::vector<uint64_t>> clustered_by_z(header.sz);

	for (uint64_t i = 0, j = 0; i < points.size(); i += 3, j++) {
		if (
			points[i+0] < 0 || points[i+0] >= sx ||
			points[i+1] < 0 || points[i+1] >= sy ||
			points[i+2] < 0 || points[i+2] >= sz

		) {
			throw std::invalid_argument("points must lie inside the volume");
		}

		clustered_by_z[points[i+2]].push_back(j);
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, 0, sz);

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(sz));

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcgs(parallel);
	std::vector<std::vector<uint32_t>> ccls(parallel);

	const uint64_t num_labels = crackle::labels::num_labels(binary);

	const bool binary_image_optimization = (
		num_labels == 2 
		&& header.label_format == LabelFormat::FLAT
		&& header.is_remapped == false
	);

	for (size_t t = 0; t < parallel; t++) {
		vcgs[t].resize(header.sx * header.sy);
		if (!binary_image_optimization) {
			ccls[t].resize(header.sx * header.sy);
		}
	}

	for (size_t z = 0; z < sz; z++) {
		if (clustered_by_z[z].size() == 0) {
			continue;
		}

		pool.enqueue([&,z](size_t t) {
			std::vector<uint8_t>& vcg = vcgs[t];
			std::vector<uint32_t>& ccl = ccls[t];
			auto& crack_code = crack_codes[z];

			crack_code_to_vcg(
				/*code=*/crack_code,
				/*sx=*/header.sx, /*sy=*/header.sy,
				/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
				/*markov_model=*/markov_model,
				/*slice_edges=*/vcg.data()
			);

			std::vector<LABEL> label_map;

			if (binary_image_optimization) {
				label_map = decode_label_map<LABEL>(
					header, binary, nullptr, 0, z, z+1
				);
				cc3d::color_binary_connectivity_graph_2d(vcg, sx, sy);
			}
			else {
				uint64_t N = 0;
				crackle::cc3d::color_connectivity_graph<uint32_t>(
					vcg, header.sx, header.sy, 1, ccl.data(), N
				);

				label_map = decode_label_map<LABEL>(
					header, binary, ccl.data(), N, z, z+1
				);
			}

			for (uint64_t idx : clustered_by_z[z]) {
				const POINT_T x = points[idx*3];
				const POINT_T y = points[idx*3 + 1];
				const uint64_t loc = (uint64_t)x + sx * (uint64_t)y;

				values[idx] = binary_image_optimization
					? label_map[vcg[loc]]
					: label_map[ccl[loc]];
			}
		});
	}

	pool.join();

	return values;
}
 
template <std::integral POINT_T>
auto scattered_points(
	const unsigned char* buffer, 
	const size_t num_bytes,
	const std::span<const POINT_T> points,
	size_t parallel = 1
) {
	CrackleHeader header(buffer);

	if (header.data_width == 1) {
		return scattered_points<uint8_t, POINT_T>(
			buffer, num_bytes,
			points, parallel
		);
	}
	else if (header.data_width == 2) {
		return scattered_points<uint16_t, POINT_T>(
			buffer, num_bytes,
			points, parallel
		);
	}
	else if (header.data_width == 4) {
		return scattered_points<uint32_t, POINT_T>(
			buffer, num_bytes,
			points, parallel
		);
	}
	else {
		return scattered_points<uint64_t, POINT_T>(
			buffer, num_bytes,
			points, parallel
		);
	}
}

template <std::integral POINT_T>
auto scattered_points(
	const std::span<const unsigned char>& buffer,
	const std::span<const POINT_T> points,
	size_t parallel = 1
) {
	return scattered_points<POINT_T>(
		buffer.data(),
		buffer.size(),
		points, parallel
	);
}

};
};

#endif