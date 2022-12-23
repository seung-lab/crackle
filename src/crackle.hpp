#ifndef __CRACKLE_HXX__
#define __CRACKLE_HXX__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
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

namespace crackle {

// COMPRESSION CODE STARTS HERE

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> encode_flat_labels(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz
) {

	const int64_t voxels = sx * sy * sz;

	std::vector<uint64_t> num_components_per_slice(sz);
	uint64_t N = 0;
	std::unique_ptr<uint32_t[]> cc_labels(crackle::cc3d::connected_components<LABEL, uint32_t>(
		labels, sx, sy, sz,
		num_components_per_slice,
		NULL, N
	));

	robin_hood::unordered_flat_map<uint32_t, LABEL> mapping;
	LABEL last = cc_labels[0];
	mapping[cc_labels[0]] = labels[0];
	for (int64_t i = 1; i < voxels; i++) {
		if (cc_labels[i] != last) {
			mapping[cc_labels[i]] = labels[i];
			last = cc_labels[i];
		}
	}

	cc_labels.release();

	robin_hood::unordered_flat_set<LABEL> uniq;
	for (auto& pair : mapping) {
		uniq.emplace(pair.second);
	}

	std::vector<LABEL> vecuniq(uniq.begin(), uniq.end());
	std::sort(vecuniq.begin(), vecuniq.end());

	robin_hood::unordered_flat_map<LABEL, STORED_LABEL> remapping;
	for (STORED_LABEL i = 0; i < vecuniq.size(); i++) {
		remapping[vecuniq[i]] = i;
	}

	std::vector<STORED_LABEL> stored_labels(N);

	for (auto [ccid, label] : mapping) {
		stored_labels[ccid] = remapping[label];
	}

	std::vector<unsigned char> binary(
		8 + sizeof(LABEL) * vecuniq.size() 
		  + sizeof(uint32_t) * num_components_per_slice.size()
		  + sizeof(STORED_LABEL) * stored_labels.size()
	);

	int64_t i = 0;
	i = crackle::lib::itoc(
		static_cast<uint64_t>(stored_labels.size()), binary, i
	);
	for (auto val : vecuniq) {
		i = crackle::lib::itoc(
			static_cast<STORED_LABEL>(val), binary, i
		);		
	}
	for (auto val : num_components_per_slice) {
		i = crackle::lib::itoc(
			static_cast<uint32_t>(val), binary, i
		);		
	}

	int key_width = crackle::lib::compute_byte_width(uniq.size());

	for (auto val : stored_labels) {
		i = crackle::lib::itocd(
			val, binary, i, sizeof(STORED_LABEL), key_width
		);		
	}

	return binary;
}

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> compress_helper(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz
) {
	const int64_t voxels = sx * sy * sz;
	const int64_t sxy = sx * sy;

	int64_t num_pairs = crackle::lib::pixel_pairs(labels, sx, sy, sz);

	CrackFormat crack_format = CrackFormat::IMPERMISSIBLE;
	LabelFormat label_format = LabelFormat::PINS_FIXED_WIDTH;
	if (num_pairs < voxels / 2) {
		crack_format = CrackFormat::PERMISSIBLE;
		label_format = LabelFormat::FLAT;
	}

	if (sz == 1) {
		label_format = LabelFormat::FLAT;
	}

	CrackleHeader header(
		/*format_version=*/0,
		/*label_format=*/label_format,
		/*crack_format=*/crack_format,
		/*data_width=*/sizeof(LABEL),
		/*stored_data_width=*/stored_data_width,
		/*sx=*/sx,
		/*sy=*/sy,
		/*sz=*/sz,
		/*num_label_bytes=*/0,
	);
	std::vector<std::vector<unsigned char>> 
		crack_codes = crackle::crackcodes::encode_boundaries(
			labels, sx, sy, sz, 
			/*permissible=*/(crack_format == CrackFormat::PERMISSIBLE)
		);
	std::vector<int64_t> z_index;
	z_index.reserve(sz);
	for (auto& code : crack_codes) {
		z_index.push_back(code.size());
	}

	if (label_format == LabelFormat::PINS_FIXED_WIDTH) {
		auto all_pins = crackle::pins::compute(labels, sx, sy, sz);

	}
	else {
		labels_binary = encode_flat_labels<LABEL, STORED_LABEL>(labels, sx, sy, sz);
	}

}

// note: labels expected to be in fortran order
template <typename LABEL>
std::vector<unsigned char> compress(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz
) {
	const int64_t voxels = sx * sy * sz;
	uint8_t stored_data_width = crackle::lib::compute_byte_width(
		crackle::lib::max_label(labels, voxels)
	);

	if (stored_data_width == 1) {
		return compress_helper<LABEL, uint8_t>(labels, sx, sy ,sz);
	}
	else if (stored_data_width == 2) {
		return compress_helper<LABEL, uint16_t>(labels, sx, sy ,sz);
	}
	else if (stored_data_width == 4) {
		return compress_helper<LABEL, uint32_t>(labels, sx, sy ,sz);
	}
	else {
		return compress_helper<LABEL, uint64_t>(labels, sx, sy ,sz);
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
	const std::vector<unsigned char> &binary
) {
	std::vector<uint64_t> z_index = get_crack_code_offsets(header, binary);

	if (z_index.back() > binary.size()) {
		throw std::runtime_error("crackle: get_crack_codes: Unable to read past end of buffer.");
	}

	std::vector<std::vector<unsigned char>> crack_codes(header.sz);

	for (uint64_t z = 0; z < header.sz; z++) {
		uint64_t code_size = z_index[z+1] - z_index[z];
		std::vector<unsigned char> code;
		code.reserve(code_size);
		for (uint64_t i = z_index[z]; i < z_index[z+1]; i++) {
			code.push_back(binary[i]);
		}
		crack_codes[z] = std::move(code);
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
	const bool fortran_order = true,
	LABEL* output = NULL
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

	const uint64_t voxels = (
		static_cast<uint64_t>(header.sx) 
		* static_cast<uint64_t>(header.sy) 
		* static_cast<uint64_t>(header.sz)
	);

	std::vector<unsigned char> binary(buffer, buffer + num_bytes);

	auto crack_codes = get_crack_codes(header, binary);
	uint64_t N = 0;
	std::vector<uint32_t> cc_labels = crack_codes_to_cc_labels<uint32_t>(
		crack_codes, header.sx, header.sy, header.sz, 
		/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE), 
		/*N=*/N
	);

	std::vector<LABEL> label_map;
	if (header.stored_data_width == 1) {
		label_map = crackle::labels::decode_label_map<LABEL, uint8_t>(
			header, binary, cc_labels, N
		);
	}
	else if (header.stored_data_width == 2) {
		label_map = crackle::labels::decode_label_map<LABEL, uint16_t>(
			header, binary, cc_labels, N
		);
	}
	else if (header.stored_data_width == 4) {
		label_map = crackle::labels::decode_label_map<LABEL, uint32_t>(
			header, binary, cc_labels, N
		);
	}
	else {
		label_map = crackle::labels::decode_label_map<LABEL, uint64_t>(
			header, binary, cc_labels, N
		);
	}

	if (output == NULL) {
		output = new LABEL[voxels]();
	}

	if (fortran_order) {
		for (uint64_t i = 0; i < voxels; i++) {
			output[i] = label_map[cc_labels[i]];
		}
	}
	else { // cc_labels is in fortran order so transpose it
		uint64_t i = 0;
		for (uint64_t z = 0; z < header.sz; z++) {
			for (uint64_t y = 0; y < header.sy; y++) {
				for (uint64_t x = 0; x < header.sx; x++, i++) {
					output[z + header.sz * (y + header.sy * x)] = label_map[cc_labels[i]];
				}
			}
		}
	}

	return output;
}

template <typename LABEL>
LABEL* decompress(
	const std::string &buffer, 
	const bool fortran_order = true,
	LABEL* output = NULL
) {
	return decompress<LABEL>(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size(),
		fortran_order,
		output
	);
}

};

#endif