#ifndef __CRACKLE_LIB_HXX__
#define __CRACKLE_LIB_HXX__

namespace crackle {
namespace lib {

// d is for dynamic
inline uint64_t itocd(uint64_t x, std::vector<unsigned char> &buf, uint64_t idx, int byte_width) { 
	for (int i = 0; i < byte_width; i++) {
		buf[idx + i] = static_cast<unsigned char>(
			(x >> (8*i)) & 0xFF
		);
	}
	return byte_width;
}

// little endian serialization of integers to chars
// returns bytes written
inline uint64_t itoc(uint8_t x, std::vector<unsigned char> &buf, uint64_t idx) {
	buf[idx] = x;
	return 1;
}

inline uint64_t itoc(uint16_t x, std::vector<unsigned char> &buf, uint64_t idx) {
	buf[idx + 0] = x & 0xFF;
	buf[idx + 1] = (x >> 8) & 0xFF;
	return 2;
}

inline uint64_t itoc(uint32_t x, std::vector<unsigned char> &buf, uint64_t idx) {
	buf[idx + 0] = x & 0xFF;
	buf[idx + 1] = (x >> 8) & 0xFF;
	buf[idx + 2] = (x >> 16) & 0xFF;
	buf[idx + 3] = (x >> 24) & 0xFF;
	return 4;
}

inline uint64_t itoc(uint64_t x, std::vector<unsigned char> &buf, uint64_t idx) {
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

template <typename T>
T ctoi(const unsigned char* buf, const uint64_t idx = 0);

template <>
uint64_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint64_t x = 0;
	x += static_cast<uint64_t>(buf[idx + 0]) << 0;
	x += static_cast<uint64_t>(buf[idx + 1]) << 8;
	x += static_cast<uint64_t>(buf[idx + 2]) << 16;
	x += static_cast<uint64_t>(buf[idx + 3]) << 24;
	x += static_cast<uint64_t>(buf[idx + 4]) << 32;
	x += static_cast<uint64_t>(buf[idx + 5]) << 40;
	x += static_cast<uint64_t>(buf[idx + 6]) << 48;
	x += static_cast<uint64_t>(buf[idx + 7]) << 56;
	return x;
}

template <>
uint32_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint32_t x = 0;
	x += static_cast<uint32_t>(buf[idx + 0]) << 0;
	x += static_cast<uint32_t>(buf[idx + 1]) << 8;
	x += static_cast<uint32_t>(buf[idx + 2]) << 16;
	x += static_cast<uint32_t>(buf[idx + 3]) << 24;
	return x;
}

template <>
uint16_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint16_t x = 0;
	x += static_cast<uint16_t>(buf[idx + 0]) << 0;
	x += static_cast<uint16_t>(buf[idx + 1]) << 8;
	return x;
}

template <>
uint8_t ctoi(const unsigned char* buf, const uint64_t idx) {
	return static_cast<uint8_t>(buf[idx]);
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
	if (x < std::numeric_limits<uint8_t>::max()) {
		return sizeof(uint8_t);
	}
	else if (x < std::numeric_limits<uint16_t>::max()) {
		return sizeof(uint16_t);
	}
	else if (x < std::numeric_limits<uint32_t>::max()) {
		return sizeof(uint32_t);
	}
	return sizeof(uint64_t);
}

template <typename LABEL>
uint64_t pixel_pairs (const LABEL* labels, const uint64_t voxels) {
	uint64_t pairs = 0;
	LABEL label = labels[0];

	for (uint64_t i = 1; i < voxels; i++) {
		if (label == labels[i]) {
			pairs++;
		}
		else {
			label = labels[i];
		}
	}
	return pairs;
}

};
};

#endif
