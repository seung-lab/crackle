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

#include "cc3d.hpp"
#include "header.hpp"
#include "crackcodes.hpp"
#include "lib.hpp"
#include "labels.hpp"

namespace crackle {

std::vector<uint64_t> get_crack_code_offsets(
	const CrackleHeader &header,
	const std::vector<unsigned char> &binary
) {
	uint64_t offset = CrackleHeader::header_size + header.num_label_bytes;

	const uint64_t z_width = header.z_index_width();
	const uint64_t zindex_bytes = z_width * header.sz;

	if (offset + zindex_bytes >= binary.size()) {
		throw std::runtime_error("crackle: get_crack_code_offsets: Unable to read past end of buffer.");
	}

	const unsigned char* buf = binary.data();

	std::vector<uint64_t> z_index(header.sz + 1);
	for (uint64_t i = 0; i < zindex_bytes; i++) {
		if (z_width == 1) {
			z_index[i+1] = lib::ctoi<uint8_t>(buf, offset + z_width * i);
		}
		else if (z_width == 2) {
			z_index[i+1] = lib::ctoi<uint16_t>(buf, offset + z_width * i);
		}
		else if (z_width == 4) {
			z_index[i+1] = lib::ctoi<uint32_t>(buf, offset + z_width * i);
		}
		else if (z_width == 8) {
			z_index[i+1] = lib::ctoi<uint64_t>(buf, offset + z_width * i);
		}
	}
	for (uint64_t z = 1; z < header.sz + 1; z++) {
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
		std::vector<unsigned char> code(code_size);
		for (uint64_t i = z_index[z]; i < z_index[z+1]; i++) {
			code.push_back(binary[i]);
		}
	}

	return crack_codes;
}

template <typename CCL>
CCL* crack_codes_to_cc_labels(
  const std::vector<std::vector<unsigned char>>& crack_codes,
  const uint64_t sx, const uint64_t sy, const uint64_t sz,
  const bool permissible, uint64_t &N
) {
	const uint64_t sxy = sx * sy;

	std::unique_ptr<uint8_t[]> edges(new uint8_t[sx*sy*sz]());

	for (uint64_t z = 0; z < crack_codes.size(); z++) {
		auto code = crackle::crackcodes::unpack_binary(crack_codes[z], sx, sy);
		uint8_t* slice_edges = crackle::crackcodes::decode_crack_code(
			code, sx, sy, permissible
		);
		for (uint64_t i = 0; i < sxy; i++) {
			edges[i + sxy*z] = slice_edges[i];
		}
	}

	return crackle::cc3d::color_connectivity_graph<CCL>(
		edges.get(), sx, sy, sz, N
	);
}

template <typename LABEL>
LABEL* decompress(
	const unsigned char* buffer, 
	size_t num_bytes,
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
	if (output == NULL) {
		output = new LABEL[voxels]();
	}

	std::vector<unsigned char> binary(buffer, buffer + num_bytes);

	auto crack_codes = get_crack_codes(header, binary);
	uint64_t N = 0;
	std::unique_ptr<uint32_t[]> cc_labels(
		crack_codes_to_cc_labels<uint32_t>(
			crack_codes, header.sx, header.sy, header.sz, 
			/*permissible=*/(header.crack_format == CrackFormat::PERMISSIBLE), 
			/*N=*/N
		)
	);

	std::vector<LABEL> label_map;
	if (header.stored_data_width == 1) {
		label_map = crackle::labels::decode_label_map<LABEL, uint8_t>(
			header, binary, cc_labels.get(), N
		);
	}
	else if (header.stored_data_width == 2) {
		label_map = crackle::labels::decode_label_map<LABEL, uint16_t>(
			header, binary, cc_labels.get(), N
		);
	}
	else if (header.stored_data_width == 4) {
		label_map = crackle::labels::decode_label_map<LABEL, uint32_t>(
			header, binary, cc_labels.get(), N
		);
	}
	else {
		label_map = crackle::labels::decode_label_map<LABEL, uint64_t>(
			header, binary, cc_labels.get(), N
		);
	}

	for (uint64_t i = 0; i < voxels; i++) {
		output[i] = label_map[cc_labels[i]];
	}

	return output;
}

template <typename LABEL>
LABEL* decompress(const std::string &buffer) {
	return decompress<LABEL>(
		reinterpret_cast<const unsigned char*>(buffer.c_str()),
		buffer.size()
	);
}

};

#endif