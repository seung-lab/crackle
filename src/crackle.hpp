#ifndef __CRACKLE_HXX__
#define __CRACKLE_HXX__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <type_traits>

#include "robin_hood.hpp"

#include "cc3d.hpp"
#include "header.hpp"
#include "crackcodes.hpp"
#include "lib.hpp"
#include "labels.hpp"
#include "pins.hpp"
#include "markov.hpp"

namespace crackle {

// COMPRESSION CODE STARTS HERE

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> compress_helper(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0
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

	CrackleHeader header(
		/*format_version=*/0,
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
		/*markov_model_order=*/markov_model_order // 0 is disabled
	);

	if (voxels == 0) {
		return header.tobytes();
	}

	std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> 
		crack_codepoints = crackle::crackcodes::encode_boundaries(
			labels, sx, sy, sz, 
			/*permissible=*/(crack_format == CrackFormat::PERMISSIBLE)
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
	if (header.markov_model_order > 0) {
		auto stats = crackle::markov::gather_statistics(crack_codepoints, header.markov_model_order);
		auto model = crackle::markov::stats_to_model(stats);
		stored_model = crackle::markov::to_stored_model(model);

		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			crack_codes[z] = crackle::markov::compress(
				crack_codepoints[z], model, header.markov_model_order,
				sx, sy
			);
		}
	}
	else {
		for (uint64_t z = 0; z < crack_codepoints.size(); z++) {
			crack_codes[z] = crackle::crackcodes::pack_codepoints(crack_codepoints[z], sx, sy);
		}
	}
	
	std::vector<unsigned char> labels_binary;
	if (label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		auto [all_pins, num_components_per_slice, num_components] = crackle::pins::compute(labels, sx, sy, sz);
		labels_binary = crackle::labels::encode_condensed_pins<LABEL, STORED_LABEL>(
			all_pins,
			sx, sy, sz,
			header.pin_index_width(),
			num_components_per_slice, num_components
		);
	}
	// else if (label_format == LabelFormat::PINS_FIXED_WIDTH) {
	// 	auto [all_pins, num_components_per_slice, num_components] = crackle::pins::compute(labels, sx, sy, sz);
	// 	labels_binary = crackle::labels::encode_fixed_width_pins<LABEL, STORED_LABEL>(
	// 		all_pins,
	// 		sx, sy, sz,
	// 		header.pin_index_width(), header.depth_width()
	// 	);		
	// }
	else {
		labels_binary = crackle::labels::encode_flat<LABEL, STORED_LABEL>(labels, sx, sy, sz);
	}

	header.num_label_bytes = labels_binary.size();

	std::vector<unsigned char> z_index_binary(sizeof(uint32_t) * sz);
	for (int64_t i = 0, z = 0; z < sz; z++) {
		i += crackle::lib::itoc(static_cast<uint32_t>(crack_codes[z].size()), z_index_binary, i);
	}

	std::vector<unsigned char> final_binary;
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
	return final_binary;
}

// note: labels expected to be in fortran order
template <typename LABEL>
std::vector<unsigned char> compress(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0
) {
	const int64_t voxels = sx * sy * sz;
	uint8_t stored_data_width = crackle::lib::compute_byte_width(
		crackle::lib::max_label(labels, voxels)
	);

	if (stored_data_width == 1) {
		return compress_helper<LABEL, uint8_t>(labels, sx, sy, sz, allow_pins, fortran_order, markov_model_order);
	}
	else if (stored_data_width == 2) {
		return compress_helper<LABEL, uint16_t>(labels, sx, sy, sz, allow_pins, fortran_order, markov_model_order);
	}
	else if (stored_data_width == 4) {
		return compress_helper<LABEL, uint32_t>(labels, sx, sy, sz, allow_pins, fortran_order, markov_model_order);
	}
	else {
		return compress_helper<LABEL, uint64_t>(labels, sx, sy, sz, allow_pins, fortran_order, markov_model_order);
	}
}


// DECOMPRESSION CODE STARTS HERE

std::vector<uint64_t> get_crack_code_offsets(
	const CrackleHeader &header,
	const std::vector<unsigned char> &binary
) {
	uint64_t offset = CrackleHeader::header_size;

	const uint64_t z_width = sizeof(uint32_t);
	const uint64_t zindex_bytes = z_width * header.sz;

	if (offset + zindex_bytes > binary.size()) {
		throw std::runtime_error("crackle: get_crack_code_offsets: Unable to read past end of buffer.");
	}

	const unsigned char* buf = binary.data();

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
			offset + zindex_bytes + 
			header.num_label_bytes + markov_model_offset
		);
	}
	return z_index;
}

std::vector<std::vector<unsigned char>> get_crack_codes(
	const CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<uint64_t> z_index = get_crack_code_offsets(header, binary);

	if (z_index.back() > binary.size()) {
		throw std::runtime_error("crackle: get_crack_codes: Unable to read past end of buffer.");
	}

	std::vector<std::vector<unsigned char>> crack_codes(z_end - z_start);

	for (uint64_t z = z_start; z < z_end; z++) {
		uint64_t code_size = z_index[z+1] - z_index[z];
		std::vector<unsigned char> code;
		code.reserve(code_size);
		for (uint64_t i = z_index[z]; i < z_index[z+1]; i++) {
			code.push_back(binary[i]);
		}
		crack_codes[z - z_start] = std::move(code);
	}

	return crack_codes;
}

std::vector<std::vector<uint8_t>> decode_markov_model(
	const CrackleHeader &header,
	const std::vector<unsigned char> &binary
) {
	if (header.markov_model_order == 0) {
		return std::vector<std::vector<uint8_t>>();
	}

	uint64_t model_offset = header.header_size + header.num_label_bytes;
	const uint64_t z_width = sizeof(uint32_t);
	const uint64_t zindex_bytes = z_width * header.sz;
	model_offset += zindex_bytes;

	std::vector<unsigned char> stored_model(
		binary.begin() + model_offset,
		binary.begin() + model_offset + header.markov_model_bytes()
	);
	return crackle::markov::from_stored_model(stored_model, header.markov_model_order);
}

// vcg: voxel connectivity graph
std::vector<uint8_t> crack_code_to_vcg(
  const std::vector<unsigned char>& code,
  const uint64_t sx, const uint64_t sy,
  const bool permissible, 
  const std::vector<std::vector<uint8_t>>& markov_model
) {
	std::vector<uint64_t> nodes = crackle::crackcodes::read_boc_index(code, sx, sy);

	std::vector<uint8_t> codepoints;
	if (markov_model.size() == 0) {
		codepoints = crackle::crackcodes::unpack_codepoints(code, sx, sy);
	}
	else {
		uint32_t index_size = 4 + crackle::lib::ctoid(code, 0, 4);
		std::vector<uint8_t> markov_stream(code.begin() + index_size, code.end());
		codepoints = crackle::markov::decode_codepoints(markov_stream, markov_model);
	}

	auto symbol_stream = crackle::crackcodes::codepoints_to_symbols(nodes, codepoints);
	return crackle::crackcodes::decode_crack_code(
		symbol_stream, sx, sy, permissible
	);
}

template <typename CCL>
std::vector<CCL> crack_codes_to_cc_labels(
  const std::vector<std::vector<unsigned char>>& crack_codes,
  const uint64_t sx, const uint64_t sy, const uint64_t sz,
  const bool permissible, uint64_t &N,
  const std::vector<std::vector<uint8_t>>& markov_model
) {
	const uint64_t sxy = sx * sy;

	std::vector<uint8_t> edges(sx*sy*sz);

	for (uint64_t z = 0; z < crack_codes.size(); z++) {
		if (crack_codes[z].size() == 0) {
			continue;
		}

		auto code = crack_codes[z];
		std::vector<uint8_t> slice_edges = crack_code_to_vcg(
			code, sx, sy,
			permissible, markov_model
		);
		for (uint64_t i = 0; i < sxy; i++) {
			edges[i + sxy*z] = slice_edges[i];
		}
	}

	return crackle::cc3d::color_connectivity_graph<CCL>(
		edges, sx, sy, sz, N
	);
}

template <typename LABEL>
LABEL* decompress(
	const unsigned char* buffer, 
	const size_t num_bytes,
	LABEL* output = NULL,
	int64_t z_start = -1,
	int64_t z_end = -1
) {
	if (num_bytes < CrackleHeader::header_size) {
		std::string err = "crackle: Input too small to be a valid stream. Bytes: ";
		err += std::to_string(num_bytes);
		throw std::runtime_error(err);
	}

	const CrackleHeader header(buffer);

	if (header.format_version > 0) {
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

	std::vector<unsigned char> binary(buffer, buffer + num_bytes);

	// only used for markov compressed streams
	std::vector<std::vector<uint8_t>> markov_model = decode_markov_model(header, binary);
	
	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);
	uint64_t N = 0;
	std::vector<uint32_t> cc_labels = crack_codes_to_cc_labels<uint32_t>(
		crack_codes, header.sx, header.sy, szr, 
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE), 
		/*N=*/N,
		/*markov_model*/markov_model
	);

	std::vector<LABEL> label_map;
	if (header.is_signed) {
		if (header.stored_data_width == 1) {
			label_map = crackle::labels::decode_label_map<LABEL, int8_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 2) {
			label_map = crackle::labels::decode_label_map<LABEL, int16_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 4) {
			label_map = crackle::labels::decode_label_map<LABEL, int32_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else {
			label_map = crackle::labels::decode_label_map<LABEL, int64_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
	}
	else {
		if (header.stored_data_width == 1) {
			label_map = crackle::labels::decode_label_map<LABEL, uint8_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 2) {
			label_map = crackle::labels::decode_label_map<LABEL, uint16_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else if (header.stored_data_width == 4) {
			label_map = crackle::labels::decode_label_map<LABEL, uint32_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
		else {
			label_map = crackle::labels::decode_label_map<LABEL, uint64_t>(
				header, binary, cc_labels, N, z_start, z_end
			);
		}
	}

	if (output == NULL) {
		output = new LABEL[voxels]();
	}

	if (header.fortran_order) {
		for (uint64_t i = 0; i < voxels; i++) {
			output[i] = label_map[cc_labels[i]];
		}
	}
	else { // cc_labels is in fortran order so transpose it
		uint64_t i = 0;
		for (uint64_t z = 0; z < static_cast<uint64_t>(szr); z++) {
			for (uint64_t y = 0; y < header.sy; y++) {
				for (uint64_t x = 0; x < header.sx; x++, i++) {
					output[z + szr * (y + header.sy * x)] = label_map[cc_labels[i]];
				}
			}
		}
	}

	return output;
}

template <typename LABEL>
LABEL* decompress(
	const std::vector<unsigned char>& buffer,
	LABEL* output = NULL,
	const int64_t z_start = -1, const int64_t z_end = -1
) {
	return decompress<LABEL>(
		buffer.data(),
		buffer.size(),
		output,
		z_start, z_end
	);
}

template <typename LABEL>
LABEL* decompress(
	const std::string &buffer,
	LABEL* output = NULL,
	const int64_t z_start = -1, const int64_t z_end = -1
) {
	return decompress<LABEL>(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size(),
		output,
		z_start, z_end
	);
}

};

#endif