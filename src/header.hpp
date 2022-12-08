#ifndef __CRACKLE_LABELS_HXX__
#define __CRACKLE_LABELS_HXX__

#include "lib.hpp"

namespace crackle {

enum LabelFormat {
	LABEL_FMT_FLAT = 0,
	LABEL_FMT_PINS_FIXED_WIDTH = 1,
	LABEL_FMT_PINS_VARIABLE_WIDTH = 2
}

enum CrackFormat {
	CRACK_FMT_IMPERMISSIBLE = 0,
	CRACK_FMT_PERMISSIBLE = 1
}

/* Header: 
 *   'crkl'            : magic number (4 bytes)
 *   format version    : unsigned integer (1 byte) 
 *   format byte       : unsigned integer (1 byte) 
 *     bits 0-1: data width (2 ** dw == byte width)
 *     bits 2-3: storate data width (2 ** sdw == byte width)
 *     bit  4:   crack format (impermissible or permissible)
 *     bits 5-6: label format
 *   sx, sy, sz        : size of each dimension (4 bytes x3)
 *   num_label_bytes   : number of label format bytes
 */
struct CrackleHeader {
public:
	static constexpr size_t header_size{22};

	static constexpr char magic[4]{ 'c', 'r', 'k', 'l' }; 
	uint8_t format_version; 
	LabelFormat label_format;
	CrackFormat crack_format;
	uint8_t data_width;
	uint8_t stored_data_width;
	uint32_t sx;
	uint32_t sy;
	uint32_t sz;
	uint32_t num_label_bytes;

	CrackleHeader() :
		format_version(0),
		label_format(LABEL_FMT_FLAT),
		crack_format(CRACK_FMT_IMPERMISSIBLE),
		data_width(1), stored_data_width(1),
		sx(1), sy(1), sz(1), 
		num_label_bytes(0)
	{}

	CrackleHeader(
		const uint8_t _format_version, 
		const LabelFormat _label_fmt,
		const CrackFormat _crack_fmt,
		const uint8_t _data_width,
		const uint8_t _stored_data_width,
		const uint32_t _sx, const uint32_t _sy, const uint32_t _sz,
		const uint32_t _num_label_bytes
	) : 
		format_version(_format_version),
		label_format(_label_fmt),
		crack_format(_crack_fmt),
		data_width(_data_width), stored_data_width(_stored_data_width),
		sx(_sx), sy(_sy), sz(_sz),
		num_label_bytes(_num_label_bytes)
	{}

	CrackleHeader(unsigned char* buf) {
		bool valid_magic = (buf[0] == 'c' && buf[1] == 'r' && buf[2] == 'k' && buf[3] == 'l');
		format_version = buf[4];

		if (!valid_magic || format_version > 0) {
			throw std::runtime_error("crackle: Data stream is not valid. Unable to decompress.");
		}

		format_byte = ctoi<uint8_t>(buf, 5);
		sx = ctoi<uint32_t>(buf, 6); 
		sy = ctoi<uint32_t>(buf, 10); 
		sz = ctoi<uint32_t>(buf, 14);
		num_label_bytes = ctoi<uint32_t>(buf, 18);

		data_width = pow(2, (format_byte & 0b00000011));
		stored_data_width = pow(2, (format_byte & 0b00001100) >> 2);
		crack_format = (format_byte & 0b00010000) >> 4;
		label_format = (label_format & 0b01100000) >> 5;
	}

	CrackleHeader(const std::vector<unsigned char> &buf) {
		return CrackleHeader(buf.data());
	}

	int z_index_width() const {
		return crackle::lib::compute_byte_width(2 * sx * sy);
	}

	int pin_index_width() const {
		return crackle::lib::compute_byte_width(sx * sy * sz);
	}

	int depth_width() const {
		return crackle::lib::compute_byte_width(sx * sy);	
	}

	size_t tochars(std::vector<unsigned char> &buf, size_t idx = 0) const {
		if ((idx + CrackleHeader::header_size) > buf.size()) {
			throw std::runtime_error("crackle: Unable to write past end of buffer.");
		}

		size_t i = idx;
		for (int j = 0; j < 4; j++, i++) {
			buf[i] = magic[j];
		}

		uint8_t format_byte = 0;
		format_byte |= static_cast<uint8_t>(log2(data_width));
		format_byte |= static_cast<uint8_t>(log2(stored_data_width)) << 2;
		format_byte |= static_cast<uint8_t>(crack_format) << 4;
		format_byte |= static_cast<uint8_t>(label_format) << 5;

		i += itoc(format_version, buf, i);
		i += itoc(format_byte, buf, i);
		i += itoc(sx, buf, i);
		i += itoc(sy, buf, i);
		i += itoc(sz, buf, i);
		i += itoc(num_label_bytes, buf, i);

		return i - idx;
	}

	static bool valid_header(unsigned char* buf) {
		bool valid_magic = (buf[0] == 'c' && buf[1] == 'r' && buf[2] == 'k' && buf[3] == 'l');
		uint8_t format_version = buf[4];
		return valid_magic && (format_version == 0);
	}

	static CrackleHeader fromchars(unsigned char* buf) {
		return CrackleHeader(buf);
	}
};

	
};

#endif

