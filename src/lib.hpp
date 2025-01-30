#ifndef __CRACKLE_LIB_HXX__
#define __CRACKLE_LIB_HXX__

#include <vector>
#include <span>

namespace crackle {
namespace lib {

// d is for dynamic
inline uint64_t itocd(uint64_t x, std::span<unsigned char> buf, uint64_t idx, int byte_width) { 
	for (int i = 0; i < byte_width; i++) {
		buf[idx + i] = static_cast<unsigned char>(
			(x >> (8*i)) & 0xFF
		);
	}
	return byte_width;
}

// little endian serialization of integers to chars
// returns bytes written
inline uint64_t itoc(uint8_t x, std::span<unsigned char> buf, uint64_t idx) {
	buf[idx] = x;
	return 1;
}

inline uint64_t itoc(uint16_t x, std::span<unsigned char> buf, uint64_t idx) {
	buf[idx + 0] = x & 0xFF;
	buf[idx + 1] = (x >> 8) & 0xFF;
	return 2;
}

inline uint64_t itoc(uint32_t x, std::span<unsigned char> buf, uint64_t idx) {
	buf[idx + 0] = x & 0xFF;
	buf[idx + 1] = (x >> 8) & 0xFF;
	buf[idx + 2] = (x >> 16) & 0xFF;
	buf[idx + 3] = (x >> 24) & 0xFF;
	return 4;
}

inline uint64_t itoc(uint64_t x, std::span<unsigned char> buf, uint64_t idx) {
	buf[idx + 0] = x & 0xFF;
	buf[idx + 1] = (x >> 8) & 0xFF;
	buf[idx + 2] = (x >> 16) & 0xFF;
	buf[idx + 3] = (x >> 24) & 0xFF;
	buf[idx + 4] = (x >> 32) & 0xFF;
	buf[idx + 5] = (x >> 40) & 0xFF;
	buf[idx + 6] = (x >> 48) & 0xFF;
	buf[idx + 7] = (x >> 56) & 0xFF;
	return 8;
}

void itoc_push_back(uint32_t x, std::vector<unsigned char> &buf) {
	buf.push_back(x & 0xFF);
	buf.push_back((x >> 8) & 0xFF);
	buf.push_back((x >> 16) & 0xFF);
	buf.push_back((x >> 24) & 0xFF);
}

template <typename T>
T ctoi(const unsigned char* buf, const uint64_t idx = 0);

template <>
int64_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int64_t x = 0;
	x |= static_cast<uint64_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint64_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint64_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint64_t>(buf[idx + 3]) << 24;
	x |= static_cast<uint64_t>(buf[idx + 4]) << 32;
	x |= static_cast<uint64_t>(buf[idx + 5]) << 40;
	x |= static_cast<uint64_t>(buf[idx + 6]) << 48;
	x |= static_cast<uint64_t>(buf[idx + 7]) << 56;
	return x;
}

template <>
uint64_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint64_t x = 0;
	x |= static_cast<uint64_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint64_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint64_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint64_t>(buf[idx + 3]) << 24;
	x |= static_cast<uint64_t>(buf[idx + 4]) << 32;
	x |= static_cast<uint64_t>(buf[idx + 5]) << 40;
	x |= static_cast<uint64_t>(buf[idx + 6]) << 48;
	x |= static_cast<uint64_t>(buf[idx + 7]) << 56;
	return x;
}

template <>
int32_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int32_t x = 0;
	x |= static_cast<uint32_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint32_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint32_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint32_t>(buf[idx + 3]) << 24;
	return x;
}

template <>
uint32_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint32_t x = 0;
	x |= static_cast<uint32_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint32_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint32_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint32_t>(buf[idx + 3]) << 24;
	return x;
}

template <>
int16_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int16_t x = 0;
	x |= static_cast<uint16_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint16_t>(buf[idx + 1]) << 8;
	return x;
}

template <>
uint16_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint16_t x = 0;
	x |= static_cast<uint16_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint16_t>(buf[idx + 1]) << 8;
	return x;
}

template <>
uint8_t ctoi(const unsigned char* buf, const uint64_t idx) {
	return static_cast<uint8_t>(buf[idx]);
}

template <>
int8_t ctoi(const unsigned char* buf, const uint64_t idx) {
	return static_cast<int8_t>(buf[idx]);
}

uint64_t ctoid(
	const unsigned char* buf, const uint64_t idx, const int byte_width
) {
	uint64_t val = 0;
	for (int i = 0; i < byte_width; i++) {
		val |= (buf[idx + i] << (i*8));
	}
	return val;
}

uint64_t ctoid(
	const std::vector<unsigned char>& buf,
	const uint64_t idx, const int byte_width
) {
	return ctoid(buf.data(), idx, byte_width);
}

uint64_t ctoid(
	const std::span<const unsigned char>& buf,
	const uint64_t idx, const int byte_width
) {
	return ctoid(buf.data(), idx, byte_width);
}

// bits must be < 8
// data must be <= 1 byte
std::vector<unsigned char> write_packed_bitstream(
	const std::vector<uint8_t>& data, const uint64_t bits
) {
	std::vector<unsigned char> output((data.size() * 8 + bits) / bits);

	uint64_t pos = 0;
	uint64_t j = 0;
	for (uint64_t i = 0; i < data.size(); i++) {
		if (pos+bits < 8) {
			output[j] |= data[i] << (8-pos-bits);
			pos += bits;
		}
		else if (pos+bits == 8) {
			output[j] |= data[i] << (8-pos-bits);
			pos = 0;
			j++;
		}
		else {
			int wrote = 8 - pos;
			output[j] |= data[i] << wrote;
			j++;
			output[j] |= data[i] >> wrote;
			pos = bits - wrote;
		}
	}
	return output;
}

std::vector<uint8_t> read_packed_bitstream(
	const std::vector<unsigned char>& bitstream, 
	const uint64_t bits, const uint64_t n_fields
) {

	std::vector<uint8_t> output(n_fields);

	uint64_t pos = 0;
	uint64_t j = 0;
	for (uint64_t i = 0; i < n_fields && i < bitstream.size(); i++) {
		if (pos+bits < 8) {
			output[i] |= (bitstream[j] << (8-pos-bits) >> (8-bits));
			pos += bits;
		}
		else if (pos+bits == 8) {
			output[i] |= (bitstream[j] >> pos);
			pos = 0;
			j++;
		}
		else {
			uint64_t read = 8 - pos;
			output[i] |= (bitstream[j] >> pos);
			j++;
			output[i] |= (bitstream[j] << (8-(bits-read)) >> (8-(bits-read)));
			pos = bits - read;
		}
	}
	
	return output;
}



template <typename LABEL>
LABEL max_label(const LABEL* labels, const uint64_t voxels) {
	LABEL mx = 0;
	if (voxels > 0) {
		mx = labels[0];
	}
	for (uint64_t i = 1; i < voxels; i++) {
		mx = std::max(mx, labels[i]);
	}
	return mx;
}

int compute_byte_width(const uint64_t x) {
	if (x <= std::numeric_limits<uint8_t>::max()) {
		return sizeof(uint8_t);
	}
	else if (x <= std::numeric_limits<uint16_t>::max()) {
		return sizeof(uint16_t);
	}
	else if (x <= std::numeric_limits<uint32_t>::max()) {
		return sizeof(uint32_t);
	}
	return sizeof(uint64_t);
}

template <typename LABEL>
uint64_t pixel_pairs (const LABEL* labels, const uint64_t voxels) {
	uint64_t pairs = 0;
	for (uint64_t i = 1; i < voxels; i++) {
		pairs += (labels[i] == labels[i-1]);
	}
	return pairs;
}

};
};

#endif
