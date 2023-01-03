#ifndef __CRACKLE_HXX__
#define __CRACKLE_HXX__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "robin_hood.hpp"

#include "cc3d.hpp"
#include "header.hpp"
#include "crackcodes.hpp"
#include "lib.hpp"
#include "labels.hpp"
#include "pins.hpp"

namespace crackle {

// COMPRESSION CODE STARTS HERE

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> compress_helper(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool allow_pins = false,
	const bool fortran_order = true
) {
	const int64_t voxels = sx * sy * sz;

	int64_t num_pairs = crackle::lib::pixel_pairs(labels, voxels);

	CrackFormat crack_format = CrackFormat::IMPERMISSIBLE;
	LabelFormat label_format = LabelFormat::PINS_FIXED_WIDTH;
	if (num_pairs < voxels / 2) {
		crack_format = CrackFormat::PERMISSIBLE;
		label_format = LabelFormat::FLAT;
	}

	if (sz == 1 || !allow_pins) {
		label_format = LabelFormat::FLAT;
	}

	CrackleHeader header(
		/*format_version=*/0,
		/*label_format=*/label_format,
		/*crack_format=*/crack_format,
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
		/*fortran_order*/fortran_order
	);
	std::vector<std::vector<unsigned char>> 
		crack_codes = crackle::crackcodes::encode_boundaries(
			labels, sx, sy, sz, 
			/*permissible=*/(crack_format == CrackFormat::PERMISSIBLE)
		);
	
	std::vector<unsigned char> labels_binary;
	if (label_format == LabelFormat::PINS_FIXED_WIDTH) {
		std::unordered_map<uint64_t, std::vector<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>>
			all_pins = crackle::pins::compute(labels, sx, sy, sz);
		labels_binary = crackle::labels::encode_fixed_width_pins<LABEL, STORED_LABEL>(
			all_pins,
			sx, sy, sz,
			header.pin_index_width(),
			header.depth_width()
		);
	}
	else if (label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		std::unordered_map<uint64_t, std::vector<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>>
			all_pins = crackle::pins::compute(labels, sx, sy, sz);
		labels_binary = crackle::labels::encode_condensed_pins<LABEL, STORED_LABEL>(
			all_pins,
			sx, sy, sz,
			header.pin_index_width()
		);
	}
	else {
		labels_binary = crackle::labels::encode_flat<LABEL, STORED_LABEL>(labels, sx, sy, sz);
	}

	header.num_label_bytes = labels_binary.size();

	std::vector<unsigned char> z_index_binary(header.z_index_width() * sz);
	for (int64_t i = 0, z = 0; z < sz; z++) {
		i += crackle::lib::itocd(crack_codes[z].size(), z_index_binary, i, header.z_index_width());
	}

	std::vector<unsigned char> final_binary;
	std::vector<unsigned char> header_binary = header.tobytes();
	final_binary.insert(final_binary.end(), header_binary.begin(), header_binary.end());
	final_binary.insert(final_binary.end(), labels_binary.begin(), labels_binary.end());
	final_binary.insert(final_binary.end(), z_index_binary.begin(), z_index_binary.end());
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
	const bool fortran_order = true
) {
	const int64_t voxels = sx * sy * sz;
	uint8_t stored_data_width = crackle::lib::compute_byte_width(
		crackle::lib::max_label(labels, voxels)
	);

	if (stored_data_width == 1) {
		return compress_helper<LABEL, uint8_t>(labels, sx, sy, sz, allow_pins, fortran_order);
	}
	else if (stored_data_width == 2) {
		return compress_helper<LABEL, uint16_t>(labels, sx, sy, sz, allow_pins, fortran_order);
	}
	else if (stored_data_width == 4) {
		return compress_helper<LABEL, uint32_t>(labels, sx, sy, sz, allow_pins, fortran_order);
	}
	else {
		return compress_helper<LABEL, uint64_t>(labels, sx, sy, sz, allow_pins, fortran_order);
	}
}


// DECOMPRESSION CODE STARTS HERE

std::vector<uint64_t> get_crack_code_offsets(
	const CrackleHeader &header,
	const std::vector<unsigned char> &binary
) {
	uint64_t offset = CrackleHeader::header_size + header.num_label_bytes;

	const uint64_t z_width = header.z_index_width();
	const uint64_t zindex_bytes = z_width * header.sz;

	if (offset + zindex_bytes > binary.size()) {
		throw std::runtime_error("crackle: get_crack_code_offsets: Unable to read past end of buffer.");
	}

	const unsigned char* buf = binary.data();

	std::vector<uint64_t> z_index(header.sz + 1);
	for (uint64_t z = 0; z < header.sz; z++) {
		if (z_width == 1) {
			z_index[z+1] = lib::ctoi<uint8_t>(buf, offset + z_width * z);
		}
		else if (z_width == 2) {
			z_index[z+1] = lib::ctoi<uint16_t>(buf, offset + z_width * z);
		}
		else if (z_width == 4) {
			z_index[z+1] = lib::ctoi<uint32_t>(buf, offset + z_width * z);
		}
		else if (z_width == 8) {
			z_index[z+1] = lib::ctoi<uint64_t>(buf, offset + z_width * z);
		}
	}
	for (uint64_t z = 0; z < header.sz; z++) {
		z_index[z+1] += z_index[z];
	}
	for (uint64_t i = 0; i < header.sz + 1; i++) {
		z_index[i] += offset + zindex_bytes;
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

template <typename CCL>
std::vector<CCL> crack_codes_to_cc_labels(
  const std::vector<std::vector<unsigned char>>& crack_codes,
  const uint64_t sx, const uint64_t sy, const uint64_t sz,
  const bool permissible, uint64_t &N
) {
	const uint64_t sxy = sx * sy;

	std::vector<uint8_t> edges(sx*sy*sz);

	for (uint64_t z = 0; z < crack_codes.size(); z++) {
		auto code = crackle::crackcodes::unpack_binary(crack_codes[z], sx, sy);
		std::vector<uint8_t> slice_edges = crackle::crackcodes::decode_crack_code(
			code, sx, sy, permissible
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

	z_start = std::max(std::min(z_start, static_cast<int64_t>(header.sz - 1)), 0LL);
	z_end = z_end < 0 ? static_cast<int64_t>(header.sz) : z_end;
	z_end = std::max(std::min(z_end, static_cast<int64_t>(header.sz)), 0LL);

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

	auto crack_codes = get_crack_codes(header, binary, z_start, z_end);
	uint64_t N = 0;
	std::vector<uint32_t> cc_labels = crack_codes_to_cc_labels<uint32_t>(
		crack_codes, header.sx, header.sy, szr, 
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE), 
		/*N=*/N
	);

	std::vector<LABEL> label_map;
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