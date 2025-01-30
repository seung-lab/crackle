#ifndef __CRACKLE_HXX__
#define __CRACKLE_HXX__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
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

// COMPRESSION CODE STARTS HERE

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> compress_helper(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0,
	const bool optimize_pins = false,
	const bool auto_bgcolor = true,
	const bool manual_bgcolor = 0,
	size_t parallel = 1
) {
	const int64_t voxels = sx * sy * sz;

	int64_t num_pairs = crackle::lib::pixel_pairs(labels, voxels);

	CrackFormat crack_format = CrackFormat::IMPERMISSIBLE;
	LabelFormat label_format = LabelFormat::PINS_VARIABLE_WIDTH;
	if (num_pairs < voxels / 2) {
		crack_format = CrackFormat::PERMISSIBLE;
		label_format = LabelFormat::FLAT;
	}

	constexpr bool is_signed = std::numeric_limits<LABEL>::is_signed;

	// Full support for signed integers in pins is not yet implemented
	// but is not impossible. The is_signed condition can be removed
	// at that point.
	if (sz == 1 || !allow_pins || is_signed) {
		label_format = LabelFormat::FLAT;
	}

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(sz));

	CrackleHeader header(
		/*format_version=*/CrackleHeader::current_version,
		/*label_format=*/label_format,
		/*crack_format=*/crack_format,
		/*is_signed=*/is_signed,
		/*data_width=*/sizeof(LABEL),
		/*stored_data_width=*/sizeof(STORED_LABEL),
		/*sx=*/sx,
		/*sy=*/sy,
		/*sz=*/sz,
		
		// grid size is not yet supported, but will be
		// used for within-slice random access.
		// very large value so that way when we start 
		// using random access decoders, old formats
		// will be backwards compatible (each slice is 1 grid).
		/*grid_size*/2147483648,
		
		/*num_label_bytes=*/0,
		/*fortran_order*/fortran_order,
		/*markov_model_order=*/markov_model_order, // 0 is disabled
		/*is_sorted=*/true,
		/*crc=*/0xFF // will be recalculated on saving
	);

	if (voxels == 0) {
		return header.tobytes();
	}

	std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> 
		crack_codepoints = crackle::crackcodes::encode_boundaries(
			labels, sx, sy, sz, 
			/*permissible=*/(crack_format == CrackFormat::PERMISSIBLE),
			/*parallel=*/parallel
		);

	if (header.markov_model_order > 0) {
		bool empty_cracks = true;
		for (auto& code : crack_codepoints) {
			if (code.size()) {
				empty_cracks = false;
				break;
			}
		}
		if (empty_cracks) {
			header.markov_model_order = 0;
		}
	}

	std::vector<unsigned char> stored_model; // only needed for markov
	std::vector<std::vector<unsigned char>> crack_codes(crack_codepoints.size());

	ThreadPool pool(parallel);

	if (header.markov_model_order > 0) {
		auto stats = crackle::markov::gather_statistics(
			crack_codepoints, header.markov_model_order, parallel
		);
		auto model = crackle::markov::stats_to_model(stats);
		stored_model = crackle::markov::to_stored_model(model);

		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			pool.enqueue([&,z](size_t t){
				crack_codes[z] = crackle::markov::compress(
					crack_codepoints[z], model, header.markov_model_order,
					sx, sy
				);
			});
		}
		// pool.join must be inside to ensure the lifetime of model
		pool.join(); 
	}
	else {
		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			pool.enqueue([&,z](size_t t){
				crack_codes[z] = crackle::crackcodes::pack_codepoints(crack_codepoints[z], sx, sy);
			});
		}
		pool.join();
	}
	
	std::vector<unsigned char> labels_binary;
	std::vector<uint32_t> crack_crcs;
	if (label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		auto [all_pins, num_components_per_slice, num_components, crack_crcs_tmp] = crackle::pins::compute(labels, sx, sy, sz, optimize_pins);
		labels_binary = crackle::labels::encode_condensed_pins<LABEL, STORED_LABEL>(
			all_pins,
			sx, sy, sz,
			header.pin_index_width(),
			num_components_per_slice, num_components,
			auto_bgcolor, manual_bgcolor
		);
		crack_crcs = std::move(crack_crcs_tmp);
	}
	else {
		auto [labels_binary_tmp, crack_crcs_tmp] = crackle::labels::encode_flat<LABEL, STORED_LABEL>(labels, sx, sy, sz, parallel);
		labels_binary = std::move(labels_binary_tmp);
		crack_crcs = std::move(crack_crcs_tmp);
	}

	header.num_label_bytes = labels_binary.size();

	std::vector<unsigned char> z_index_binary;
	z_index_binary.reserve(sizeof(uint32_t) * (sz + 1));
	z_index_binary.resize(sizeof(uint32_t) * sz);

	int64_t i = 0, z = 0;
	uint64_t code_size = 0;
	for (; z < sz; z++) {
		i += crackle::lib::itoc(static_cast<uint32_t>(crack_codes[z].size()), z_index_binary, i);
		code_size += crack_codes[z].size();
	}
	const uint32_t z_index_crc = crackle::crc::crc32c(z_index_binary);
	z_index_binary.resize(sizeof(uint32_t) * (sz + 1));
	i += crackle::lib::itoc(z_index_crc, z_index_binary, i);

	const uint32_t labels_binary_crc = crackle::crc::crc32c(labels_binary);

	std::vector<unsigned char> final_binary;
	final_binary.reserve(
		header.header_bytes()
		+ z_index_binary.size()
		+ labels_binary.size()
		+ stored_model.size()
		+ code_size
		+ sizeof(uint32_t) // labels crc32c
		+ header.sz * sizeof(uint32_t) // crack codes crc32cs
	);

	std::vector<unsigned char> header_binary = header.tobytes();
	final_binary.insert(final_binary.end(), header_binary.begin(), header_binary.end());
	final_binary.insert(final_binary.end(), z_index_binary.begin(), z_index_binary.end());
	final_binary.insert(final_binary.end(), labels_binary.begin(), labels_binary.end());
	if (header.markov_model_order > 0) {
		final_binary.insert(final_binary.end(), stored_model.begin(), stored_model.end());
	}
	for (auto& code : crack_codes) {
		final_binary.insert(final_binary.end(), code.begin(), code.end());
	}

	crackle::lib::itoc_push_back(labels_binary_crc, final_binary);
	for (uint32_t crack_crc : crack_crcs) {
		crackle::lib::itoc_push_back(crack_crc, final_binary);
	}

	return final_binary;
}

// note: labels expected to be in fortran order
template <typename LABEL>
std::vector<unsigned char> compress(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0,
	const bool optimize_pins = false,
	const bool auto_bgcolor = true,
	const int64_t manual_bgcolor = 0,
	size_t parallel = 1
) {
	const int64_t voxels = sx * sy * sz;
	uint8_t stored_data_width = crackle::lib::compute_byte_width(
		crackle::lib::max_label(labels, voxels)
	);

#define CALL_COMPRESS_HELPER(STORED_T) compress_helper<LABEL, STORED_T>(\
		labels, sx, sy, sz,\
		allow_pins, fortran_order, markov_model_order,\
		optimize_pins, auto_bgcolor, manual_bgcolor, parallel\
	)

	if (stored_data_width == 1) {
		return CALL_COMPRESS_HELPER(uint8_t);
	}
	else if (stored_data_width == 2) {
		return CALL_COMPRESS_HELPER(uint16_t);
	}
	else if (stored_data_width == 4) {
		return CALL_COMPRESS_HELPER(uint32_t);
	}
	else {
		return CALL_COMPRESS_HELPER(uint64_t);
	}

#undef CALL_COMPRESS_HELPER
}


// DECOMPRESSION CODE STARTS HERE

std::vector<uint64_t> get_crack_code_offsets(
	const CrackleHeader &header,
	const std::span<const unsigned char> &binary
) {
	uint64_t offset = header.header_bytes();

	const uint64_t z_width = sizeof(uint32_t);

	if (offset + header.grid_index_bytes() > binary.size()) {
		throw std::runtime_error("crackle: get_crack_code_offsets: Unable to read past end of buffer.");
	}

	const unsigned char* buf = binary.data();

	if (header.format_version > 0) {
		const uint32_t stored_crc32c = crackle::lib::ctoi<uint32_t>(
			buf, offset + z_width * header.sz
		);
		const uint32_t computed_crc32c = crackle::crc::crc32c(
			buf + offset, header.grid_index_bytes() - sizeof(uint32_t)
		);

		if (stored_crc32c != computed_crc32c) {
			std::string err = "crackle: grid index crc32c did not match. stored: ";
			err += std::to_string(stored_crc32c);
			err += " computed: ";
			err += std::to_string(computed_crc32c);
			throw std::runtime_error(err);
		}
	}

	std::vector<uint64_t> z_index(header.sz + 1);
	for (uint64_t z = 0; z < header.sz; z++) {
		z_index[z+1] = lib::ctoi<uint32_t>(buf, offset + z_width * z);
	}
	for (uint64_t z = 0; z < header.sz; z++) {
		z_index[z+1] += z_index[z];
	}

	uint64_t markov_model_offset = 0;
	if (header.markov_model_order > 0) {
		markov_model_offset = header.markov_model_bytes();
	}

	for (uint64_t i = 0; i < header.sz + 1; i++) {
		z_index[i] += (
			offset + header.grid_index_bytes() + 
			header.num_label_bytes + markov_model_offset
		);
	}
	return z_index;
}

std::vector<std::span<const unsigned char>> get_crack_codes(
	const CrackleHeader &header,
	const std::span<const unsigned char> &binary,
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<uint64_t> z_index = get_crack_code_offsets(header, binary);

	if (z_index.back() > binary.size()) {
		throw std::runtime_error("crackle: get_crack_codes: Unable to read past end of buffer.");
	}

	std::vector<std::span<const unsigned char>> crack_codes(z_end - z_start);

	for (uint64_t z = z_start; z < z_end; z++) {
		uint64_t code_size = z_index[z+1] - z_index[z];
		crack_codes[z - z_start] = std::span<const unsigned char>(
			&binary[z_index[z]], code_size
		);
	}

	return crack_codes;
}

uint32_t get_labels_crc(
	const CrackleHeader &header,
	const std::span<const unsigned char> &binary
) {
	uint64_t offset = (header.sz + 1) * sizeof(uint32_t);
	return crackle::lib::ctoi<uint32_t>(binary.data(), binary.size() - offset);
}

uint32_t reset_labels_crc(std::span<unsigned char> &binary) {
	crackle::CrackleHeader header(binary);
	std::span<const unsigned char> labels_binary = crackle::labels::raw_labels(binary);
	const uint32_t labels_binary_crc = crackle::crc::crc32c(labels_binary);

	uint64_t offset = (header.sz + 1) * sizeof(uint32_t);
	crackle::lib::itoc(labels_binary_crc, binary, binary.size() - offset);
	return labels_binary_crc;
}

std::span<const uint32_t> get_crack_code_crcs(
	const CrackleHeader &header,
	const std::span<const unsigned char> &binary
) {
  if (binary.size() < header.sz * sizeof(uint32_t)) {
  	throw std::out_of_range("Insufficient binary data for crack code CRCs.");
  }

	// Compute the start of the uint32_t array
	const uint32_t* start = reinterpret_cast<const uint32_t*>(
		binary.data() + (binary.size() - header.sz * sizeof(uint32_t))
	);

	return std::span<const uint32_t>(start, header.sz);
}

std::vector<std::vector<uint8_t>> decode_markov_model(
	const CrackleHeader &header,
	const std::span<const unsigned char> &binary
) {
	if (header.markov_model_order == 0) {
		return std::vector<std::vector<uint8_t>>();
	}

	uint64_t model_offset = (
		header.header_bytes() 
		+ header.grid_index_bytes()
		+ header.num_label_bytes
	);

	std::vector<unsigned char> stored_model(
		binary.begin() + model_offset,
		binary.begin() + model_offset + header.markov_model_bytes()
	);
	return crackle::markov::from_stored_model(stored_model, header.markov_model_order);
}

std::vector<std::pair<uint64_t, std::vector<unsigned char>> >
crack_code_to_symbols(
  const std::span<const unsigned char>& code,
  const uint64_t sx, const uint64_t sy,
  const std::vector<std::vector<uint8_t>>& markov_model
) {
	std::vector<uint64_t> nodes = crackle::crackcodes::read_boc_index(code, sx, sy);

	std::vector<uint8_t> codepoints;
	if (markov_model.size() == 0) {
		codepoints = crackle::crackcodes::unpack_codepoints(code, sx, sy);
	}
	else {
		uint32_t index_size = 4 + crackle::lib::ctoid(code, 0, 4);
		std::span<const uint8_t> markov_stream(code.begin() + index_size, code.size() - index_size);
		codepoints = crackle::markov::decode_codepoints(markov_stream, markov_model);
	}

	return crackle::crackcodes::codepoints_to_symbols(nodes, codepoints);
}

// vcg: voxel connectivity graph
void crack_code_to_vcg(
  const std::span<const unsigned char>& code,
  const uint64_t sx, const uint64_t sy,
  const bool permissible, 
  const std::vector<std::vector<uint8_t>>& markov_model,
  uint8_t* slice_edges
) {
	auto symbol_stream = crack_code_to_symbols(code, sx, sy, markov_model);
	crackle::crackcodes::decode_crack_code(
		symbol_stream, sx, sy, permissible, slice_edges
	);
}

template <typename CCL>
CCL* crack_codes_to_cc_labels(
  std::span<const unsigned char>& crack_codes,
  const uint64_t sx, const uint64_t sy,
  const bool permissible, uint64_t &N,
  const std::vector<std::vector<uint8_t>>& markov_model,
  std::vector<uint8_t>& vcg,
  CCL* out = NULL
) {
	crack_code_to_vcg(
		crack_codes, sx, sy,
		permissible, markov_model,
		vcg.data()
	);

	return crackle::cc3d::color_connectivity_graph<CCL>(
		vcg, sx, sy, 1, out, N
	);
}

template <typename LABEL>
std::vector<LABEL> decode_label_map(
	const CrackleHeader &header,
	const std::span<const unsigned char>& binary,
	const uint32_t* cc_labels,
	uint64_t N,
	int64_t z_start,
	int64_t z_end
) {
	if (header.is_signed) {
		if (header.stored_data_width == 1) {
			return crackle::labels::decode_label_map<LABEL, int8_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 2) {
			return crackle::labels::decode_label_map<LABEL, int16_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 4) {
			return crackle::labels::decode_label_map<LABEL, int32_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else {
			return crackle::labels::decode_label_map<LABEL, int64_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
	}
	else {
		if (header.stored_data_width == 1) {
			return crackle::labels::decode_label_map<LABEL, uint8_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 2) {
			return crackle::labels::decode_label_map<LABEL, uint16_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 4) {
			return crackle::labels::decode_label_map<LABEL, uint32_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else {
			return crackle::labels::decode_label_map<LABEL, uint64_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
	}
}


template <typename LABEL, typename OUT>
OUT* decompress(
	const unsigned char* buffer, 
	const size_t num_bytes,
	OUT* output = NULL,
	int64_t z_start = -1,
	int64_t z_end = -1,
	size_t parallel = 1,
	const std::optional<uint64_t> label = std::nullopt
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
		return output;
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);

	if (output == NULL) {
		output = new OUT[voxels]();
	}

	const uint64_t sxy = header.sx * header.sy;

	std::span<const uint32_t> crack_code_crcs;

	if (header.format_version > 0) {
		crack_code_crcs = get_crack_code_crcs(header, binary);
	}

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(szr));

	ThreadPool pool(parallel);

	std::vector<std::vector<uint8_t>> vcg_scratch(parallel);
	std::vector<std::vector<uint32_t>> cc_labels_scratch(parallel);
	for (size_t t = 0; t < parallel; t++) {
		vcg_scratch[t].resize(sxy);
		cc_labels_scratch[t].resize(sxy);
	}

	for (uint64_t z = 0; z < static_cast<uint64_t>(szr); z++) {
		pool.enqueue([&,z](size_t t) {
			std::vector<uint8_t>& vcg = vcg_scratch[t];
			std::vector<uint32_t>& cc_labels = cc_labels_scratch[t];

			uint64_t N = 0;
			crack_codes_to_cc_labels<uint32_t>(
				crack_codes[z], header.sx, header.sy,
				/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE), 
				/*N=*/N,
				/*markov_model=*/markov_model,
				/*vcg=*/vcg,
				/*output=*/cc_labels.data()
			);

			if (header.format_version > 0) {
				const uint32_t computed_crc = crackle::crc::crc32c(cc_labels.data(), sxy);

				if (crack_code_crcs[z_start + z] != computed_crc) {
					std::string err = "crackle: crack code crc mismatch on z=";
					err += std::to_string(z_start + z);
					err += " computed: ";
					err += std::to_string(computed_crc);
					err += " stored: ";
					err += std::to_string(crack_code_crcs[z_start + z]);
					throw std::runtime_error(err);
				}
			}

			const std::vector<LABEL> label_map = decode_label_map<LABEL>(
				header, binary, cc_labels.data(), N, z_start+z, z_start+z+1
			);

			if (label.has_value()) {
				const uint64_t label_v = *label;

				if (header.fortran_order) {
					for (uint64_t i = 0; i < sxy; i++) {
						output[i + z * sxy] = (label_map[cc_labels[i]] == label_v);
					}
				}
				else {
					uint64_t i = 0;
					for (uint64_t y = 0; y < header.sy; y++) {
						for (uint64_t x = 0; x < header.sx; x++, i++) {
							output[z + szr * (y + header.sy * x)] = (label_map[cc_labels[i]] == label_v);
						}
					}
				}
			}
			else {
				if (header.fortran_order) {
					for (uint64_t i = 0; i < sxy; i++) {
						output[i + z * sxy] = label_map[cc_labels[i]];
					}
				}
				else {
					uint64_t i = 0;
					for (uint64_t y = 0; y < header.sy; y++) {
						for (uint64_t x = 0; x < header.sx; x++, i++) {
							output[z + szr * (y + header.sy * x)] = label_map[cc_labels[i]];
						}
					}
				}
			}
		});
	}

	pool.join();

	return output;
}

template <typename LABEL, typename OUT>
LABEL* decompress(
	const std::vector<unsigned char>& buffer,
	LABEL* output = NULL,
	const int64_t z_start = -1, const int64_t z_end = -1,
	size_t parallel = 1,
	const std::optional<uint64_t> label = std::nullopt
) {
	return decompress<LABEL, OUT>(
		buffer.data(),
		buffer.size(),
		output,
		z_start, z_end,
		parallel,
		label
	);
}

template <typename LABEL, typename OUT>
LABEL* decompress(
	const std::span<const unsigned char>& buffer,
	LABEL* output = NULL,
	const int64_t z_start = -1, const int64_t z_end = -1,
	size_t parallel = 1, 	
	const std::optional<uint64_t> label = std::nullopt
) {
	return decompress<LABEL, OUT>(
		buffer.data(),
		buffer.size(),
		output,
		z_start, z_end,
		parallel,
		label
	);
}


template <typename LABEL, typename OUT>
LABEL* decompress(
	const std::string &buffer,
	LABEL* output = NULL,
	const int64_t z_start = -1, const int64_t z_end = -1,
	size_t parallel = 1,
	const std::optional<uint64_t> label = std::nullopt
) {
	return decompress<LABEL, OUT>(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size(),
		output,
		z_start, z_end,
		parallel,
		label
	);
}

void decompress(
	const unsigned char* buffer, 
	const size_t num_bytes,
	unsigned char* output, const int64_t out_num_bytes,
	const int64_t z_start = -1, const int64_t z_end = -1,
	size_t parallel = 1,
	const std::optional<uint64_t> label = std::nullopt
) {
	CrackleHeader header(buffer);

	if (header.nbytes() > static_cast<uint64_t>(out_num_bytes)) {
		throw new std::runtime_error("Output buffer too small.");
	}

	if (header.data_width == 1) {
		decompress<uint8_t, uint8_t>(
			buffer, num_bytes, reinterpret_cast<uint8_t*>(output), 
			z_start, z_end,
			parallel, label
		);
	}
	else if (header.data_width == 2) {
		if (label.has_value()) {
			decompress<uint16_t, uint8_t>(
				buffer, num_bytes, reinterpret_cast<uint8_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
		else {
			decompress<uint16_t, uint16_t>(
				buffer, num_bytes, reinterpret_cast<uint16_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
	}
	else if (header.data_width == 4) {
		if (label.has_value()) {
			decompress<uint32_t, uint8_t>(
				buffer, num_bytes, reinterpret_cast<uint8_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
		else {
			decompress<uint32_t, uint32_t>(
				buffer, num_bytes, reinterpret_cast<uint32_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
	}
	else {
		if (label.has_value()) {
			decompress<uint32_t, uint8_t>(
				buffer, num_bytes, reinterpret_cast<uint8_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
		else {
			decompress<uint32_t, uint32_t>(
				buffer, num_bytes, reinterpret_cast<uint32_t*>(output), 
				z_start, z_end,
				parallel, label
			);
		}
	}
}

std::vector<uint8_t>
decode_slice_vcg(
	const unsigned char* buffer, 
	const size_t num_bytes,
	int64_t z
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
	if (z >= header.sz || z < 0) {
		std::string err = "crackle: Invalid z: ";
		err += std::to_string(z);
		throw std::runtime_error(err);
	}

	const uint64_t voxels = (
		static_cast<uint64_t>(header.sx) 
		* static_cast<uint64_t>(header.sy) 
	);

	if (voxels == 0) {
		return std::vector<uint8_t>();
	}

	std::span<const unsigned char> binary(buffer, num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z, z+1);
	std::vector<uint8_t> vcg(header.sx * header.sy);

	for (auto crack_code : crack_codes) {
		crack_code_to_vcg(
			/*code=*/crack_code,
			/*sx=*/header.sx, /*sy=*/header.sy,
			/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE),
			/*markov_model=*/markov_model,
			/*slice_edges=*/vcg.data()
		);
	}

	return vcg;
}

auto decode_slice_vcg(
	const std::string &buffer,
	const int64_t z
) {
	return decode_slice_vcg(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size(),
		z
	);
}

// take an existing crackle stream and change the markov order
// returning a reencoded copy
std::vector<unsigned char> reencode_with_markov_order(
	const unsigned char* buffer, 
	const size_t num_bytes,
	const int markov_model_order,
	size_t parallel = 1
) { 
	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	CrackleHeader header(buffer);

	if (header.format_version > CrackleHeader::current_version) {
		std::string err = "crackle: Invalid format version.";
		err += std::to_string(header.format_version);
		throw std::runtime_error(err);
	}

	if (parallel == 0) {
		parallel = std::thread::hardware_concurrency();
	}
	parallel = std::min(parallel, static_cast<size_t>(header.sz));

	std::span<const unsigned char> binary(buffer, num_bytes);

	if (header.markov_model_order == markov_model_order) {
		return std::vector<unsigned char>(binary.begin(), binary.end());
	}

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto existing_crack_codes = get_crack_codes(header, binary, 0, header.sz);

	std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> crack_codepoints;
	for (auto& crack_code : existing_crack_codes) {
		auto symbol_stream = crack_code_to_symbols(crack_code, header.sx, header.sy, markov_model);
		auto unpacked_codepoints = crackle::crackcodes::symbols_to_codepoints(symbol_stream);
		crack_codepoints.push_back(unpacked_codepoints);
	}

	existing_crack_codes = std::vector<std::span<const unsigned char>>();

	header.markov_model_order = markov_model_order;

	std::vector<unsigned char> stored_model; // only needed for markov
	std::vector<std::vector<unsigned char>> crack_codes(crack_codepoints.size());

	ThreadPool pool(parallel);

	if (header.markov_model_order > 0) {
		auto stats = crackle::markov::gather_statistics(
			crack_codepoints, header.markov_model_order, parallel
		);
		auto model = crackle::markov::stats_to_model(stats);
		stored_model = crackle::markov::to_stored_model(model);

		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			pool.enqueue([&,z](size_t t){
				crack_codes[z] = crackle::markov::compress(
					crack_codepoints[z], model, header.markov_model_order,
					header.sx, header.sy
				);
			});
		}
		// pool.join must be inside to ensure the lifetime of model
		pool.join(); 
	}
	else {
		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			pool.enqueue([&,z](size_t t){
				crack_codes[z] = crackle::crackcodes::pack_codepoints(crack_codepoints[z], header.sx, header.sy);
			});
		}
		pool.join();
	}

	uint64_t code_size = 0;
	std::vector<unsigned char> z_index_binary(sizeof(uint32_t) * header.sz);
	int64_t i = 0, z = 0;
	for (; z < header.sz; z++) {
		code_size += crack_codes[z].size();
		i += crackle::lib::itoc(static_cast<uint32_t>(crack_codes[z].size()), z_index_binary, i);
	}

	if (header.format_version > 0) {
		const uint32_t z_index_crc = crackle::crc::crc32c(z_index_binary.data(), header.sz * sizeof(uint32_t));
		z_index_binary.resize(z_index_binary.size() + sizeof(z_index_crc));
		i += crackle::lib::itoc(z_index_crc, z_index_binary, i);
	}

	std::span<const unsigned char> labels_binary = crackle::labels::raw_labels(binary);

	std::vector<unsigned char> final_binary = header.tobytes();
	final_binary.reserve(
		final_binary.size()
		+ z_index_binary.size() + labels_binary.size() 
		+ stored_model.size() + code_size
	);

	final_binary.insert(final_binary.end(), z_index_binary.begin(), z_index_binary.end());
	final_binary.insert(final_binary.end(), labels_binary.begin(), labels_binary.end());
	if (header.markov_model_order > 0) {
		final_binary.insert(final_binary.end(), stored_model.begin(), stored_model.end());
	}
	for (auto& code : crack_codes) {
		final_binary.insert(final_binary.end(), code.begin(), code.end());
	}

	if (header.format_version == 0) {
		return final_binary;
	}

	uint32_t labels_binary_crc = get_labels_crc(header, binary);
	std::span<const uint32_t> crack_crcs = get_crack_code_crcs(header, binary);

	crackle::lib::itoc_push_back(labels_binary_crc, final_binary);
	for (uint32_t crack_crc : crack_crcs) {
		crackle::lib::itoc_push_back(crack_crc, final_binary);
	}

	return final_binary;
}

auto reencode_with_markov_order(
	const std::string &buffer,
	const int markov_model_order,
	size_t parallel = 1
) {
	return reencode_with_markov_order(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size(),
		markov_model_order, parallel
	);
}

// MAPTYPE b/c maybe someone will want to use the C++ interface and use 
// std::unordered_map or a wrapper around py::dict directly to avoid copies.
template <typename STORED_LABEL, typename MAPTYPE = robin_hood::unordered_flat_map<uint64_t,uint64_t>>
void remap(
	unsigned char* buffer,
	const uint64_t num_bytes,
	const MAPTYPE& mapping,
	const bool preserve_missing_labels = false,
	size_t parallel = 0
) {

	std::span<unsigned char> binary(buffer, num_bytes);

	crackle::CrackleHeader header(binary);

	std::span<const unsigned char> labels_binary = crackle::labels::raw_labels(binary);

	if (header.format_version > 0) {
		uint32_t crc_computed = crackle::crc::crc32c(labels_binary);
		uint32_t crc_stored = get_labels_crc(header, binary);
		
		if (crc_computed != crc_stored) {
			throw std::runtime_error("crackle::remap: crc mismatch on labels binary.");
		}
	}

	std::span<const STORED_LABEL> unique_span = crackle::labels::decode_uniq<STORED_LABEL>(
		header, labels_binary
	);

	const uint64_t target_size = unique_span.size() + (header.label_format == LabelFormat::PINS_VARIABLE_WIDTH);

	if (!preserve_missing_labels && mapping.size() < target_size) {
		throw std::runtime_error("mapping must be at least as large as unique labels when preserve_missing_labels is false.");
	}

	STORED_LABEL* unique_ptr = const_cast<STORED_LABEL*>(unique_span.data());
	std::span<STORED_LABEL> unique(unique_ptr, unique_span.size());

	const uint64_t target_factor = static_cast<uint64_t>(ceil(log2(target_size)));

	if (header.is_sorted && preserve_missing_labels && mapping.size() < target_size / target_factor) {
		std::vector<typename std::span<STORED_LABEL>::iterator> iters;
		iters.reserve(mapping.size());

		for (auto& [k,v] : mapping) {
    	auto it = std::lower_bound(unique.begin(), unique.end(), k);
	    if (it != unique.end() && *it == k) {
	    	iters.push_back(it);
			}
		}
		for (auto& it : iters) {
			auto itm = mapping.find(*it);
			*it = itm->second;
		}
	}
	else {
		if (parallel == 0) {
			parallel = std::thread::hardware_concurrency();
		}

		const uint64_t step = 1e5;

		parallel = std::min(parallel, static_cast<size_t>(unique.size() / step));
		parallel = std::max(parallel, static_cast<size_t>(1));

		ThreadPool pool(parallel);

		for (uint64_t i = 0; i < unique.size(); i += step) {
			pool.enqueue([&,i](size_t t) {
				const uint64_t limit = std::min(i + step, static_cast<uint64_t>(unique.size()));
				for (uint64_t j = i; j < limit; j++) {
					const uint64_t uniq = static_cast<uint64_t>(unique[j]);
					auto it = mapping.find(uniq);
					if (it != mapping.end()) {
						unique[j] = it->second;
					}
					else if (preserve_missing_labels == false) {
						throw std::runtime_error("Label was missing.");
					}
				}
			});
		}

		pool.join();
	}

	if (header.label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		STORED_LABEL bgcolor = crackle::labels::background_color(binary);
		auto it = mapping.find(bgcolor);
		if (it != mapping.end()) {
			bgcolor = it->second;
			unsigned char* lbl_binary = const_cast<unsigned char*>(labels_binary.data());
			std::span<unsigned char> lbl_binary_span(lbl_binary, sizeof(STORED_LABEL));
			crackle::lib::itocd(bgcolor, lbl_binary_span, 0, sizeof(STORED_LABEL));
		}
	}

	header.is_sorted = true;
	for (uint64_t i = 1; i < unique.size(); i++) {
		if (unique[i] < unique[i-1]) {
			header.is_sorted = false;
			break;
		}
	}

	header.tochars(binary);

	if (header.format_version > 0) {
		reset_labels_crc(binary);
	}
}

};

#endif
