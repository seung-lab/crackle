#ifndef __CRACKLE_HEADER_HXX__
#define __CRACKLE_HEADER_HXX__

#include "lib.hpp"
#include <cstdint>

namespace crackle {

enum LabelFormat {
	FLAT = 0,
	PINS_FIXED_WIDTH = 1,
	PINS_VARIABLE_WIDTH = 2
};

enum CrackFormat {
	IMPERMISSIBLE = 0,
	PERMISSIBLE = 1
};

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
	static constexpr size_t header_size{24};

	static constexpr char magic[4]{ 'c', 'r', 'k', 'l' }; 
	uint8_t format_version; 
	LabelFormat label_format;
	CrackFormat crack_format;
	bool is_signed;
	uint8_t data_width;
	uint8_t stored_data_width;
	uint32_t sx;
	uint32_t sy;
	uint32_t sz;
	uint32_t grid_size;
	uint32_t num_label_bytes;
	bool fortran_order;

	CrackleHeader() :
		format_version(0),
		label_format(LabelFormat::FLAT),
		crack_format(CrackFormat::IMPERMISSIBLE),
		is_signed(false),
		data_width(1), stored_data_width(1),
		sx(1), sy(1), sz(1), grid_size(2147483648),
		num_label_bytes(0), fortran_order(true)
	{}

	CrackleHeader(
		const uint8_t _format_version, 
		const LabelFormat _label_fmt,
		const CrackFormat _crack_fmt,
		const bool _is_signed,
		const uint8_t _data_width,
		const uint8_t _stored_data_width,
		const uint32_t _sx, const uint32_t _sy, const uint32_t _sz,
		const uint32_t _grid_size,
		const uint32_t _num_label_bytes,
		const bool _fortran_order
	) : 
		format_version(_format_version),
		label_format(_label_fmt),
		crack_format(_crack_fmt),
		is_signed(_is_signed),
		data_width(_data_width), stored_data_width(_stored_data_width),
		sx(_sx), sy(_sy), sz(_sz),
		grid_size(_grid_size),
		num_label_bytes(_num_label_bytes), 
		fortran_order(_fortran_order)
	{}

	void assign_from_buffer(const unsigned char* buf) {
		bool valid_magic = (buf[0] == 'c' && buf[1] == 'r' && buf[2] == 'k' && buf[3] == 'l');
		format_version = buf[4];

		if (!valid_magic || format_version > 0) {
			throw std::runtime_error("crackle: Data stream is not valid. Unable to decompress.");
		}

		uint16_t format_bytes = lib::ctoi<uint16_t>(buf, 5);
		sx = lib::ctoi<uint32_t>(buf, 7); 
		sy = lib::ctoi<uint32_t>(buf, 11); 
		sz = lib::ctoi<uint32_t>(buf, 15);
		grid_size = static_cast<uint32_t>(
			pow(2, lib::ctoi<uint8_t>(buf, 19))
		);
		num_label_bytes = lib::ctoi<uint32_t>(buf, 20);

		data_width = pow(2, (format_bytes & 0b00000011));
		stored_data_width = pow(2, (format_bytes & 0b00001100) >> 2);
		crack_format = static_cast<CrackFormat>((format_bytes & 0b00010000) >> 4);
		label_format = static_cast<LabelFormat>((format_bytes & 0b01100000) >> 5);
		fortran_order = static_cast<bool>((format_bytes & 0b10000000) >> 7);
		is_signed = static_cast<bool>((format_bytes >> 8) & 0b1);
	}

	CrackleHeader(const unsigned char* buf) {
		assign_from_buffer(buf);
	}

	CrackleHeader(const std::string &buf) {
		assign_from_buffer(reinterpret_cast<const unsigned char*>(buf.c_str()));
	}

	CrackleHeader(const std::vector<unsigned char> &buf) {
		assign_from_buffer(buf.data());
	}

	uint64_t voxels() const {
		return static_cast<uint64_t>(sx) * static_cast<uint64_t>(sy) * static_cast<uint64_t>(sz);
	}

	int pin_index_width() const {
		return crackle::lib::compute_byte_width(sx * sy * sz);
	}

	int depth_width() const {
		return crackle::lib::compute_byte_width(sz == 0 ? 0 : sz - 1);	
	}

	uint64_t num_grids() const {
		uint64_t gsize = std::min(grid_size, std::max(sx, sy));
		uint64_t ngrids = ((sx + gsize - 1) / gsize) * ((sy + gsize - 1) / gsize);
		ngrids = std::max(ngrids, static_cast<uint64_t>(1));
		ngrids *= sz;
		return ngrids;
	}

	size_t tochars(std::vector<unsigned char> &buf, size_t idx = 0) const {
		if ((idx + CrackleHeader::header_size) > buf.size()) {
			throw std::runtime_error("crackle: Unable to write past end of buffer.");
		}

		size_t i = idx;
		for (int j = 0; j < 4; j++, i++) {
			buf[i] = magic[j];
		}

		uint16_t format_bytes = 0;
		format_bytes |= static_cast<uint16_t>(log2(data_width));
		format_bytes |= static_cast<uint16_t>(log2(stored_data_width)) << 2;
		format_bytes |= static_cast<uint16_t>(crack_format) << 4;
		format_bytes |= static_cast<uint16_t>(label_format) << 5;
		format_bytes |= static_cast<uint16_t>(fortran_order) << 7;
		format_bytes |= static_cast<uint16_t>(is_signed) << 8;

		i += lib::itoc(format_version, buf, i);
		i += lib::itoc(format_bytes, buf, i);
		i += lib::itoc(sx, buf, i);
		i += lib::itoc(sy, buf, i);
		i += lib::itoc(sz, buf, i);
		i += lib::itoc(static_cast<uint8_t>(log2(grid_size)), buf, i);
		i += lib::itoc(num_label_bytes, buf, i);

		return i - idx;
	}

	std::vector<unsigned char> tobytes() const {
		std::vector<unsigned char> buf(header_size);
		tochars(buf);
		return buf;
	}

	uint64_t nbytes() const {
		return (
			  static_cast<uint64_t>(sx) 
			* static_cast<uint64_t>(sy) 
			* static_cast<uint64_t>(sz) 
			* static_cast<uint64_t>(data_width)
		);
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

