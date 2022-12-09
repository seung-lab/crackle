#ifndef __CRACKLE_LABELS_HXX__
#define __CRACKLE_LABELS_HXX__

#include <vector>

#include "header.hpp"
#include "lib.hpp"

namespace crackle {
namespace labels {

std::vector<unsigned char> raw_labels(
	const std::vector<unsigned char> &binary
) {
	crackle::CrackleHeader header(binary);
	uint64_t hb = crackle::CrackleHeader::header_size;
	std::vector<unsigned char> labels_binary(
		binary.begin() + hb,
		binary.begin() + hb + header.num_label_bytes
	);
	return labels_binary;
}

uint64_t decode_num_labels(
	const CrackleHeader &header,
	const std::vector<unsigned char> &labels_binary
) {
	if (header.label_format == LabelFormat::FLAT) {
		return crackle::lib::ctoi<uint64_t>(labels_binary.data(), 0);
	}
	else {
		return crackle::lib::ctoi<uint64_t>(labels_binary.data(), header.stored_data_width);
	}
}

template <typename STORED_LABEL>
std::vector<STORED_LABEL> decode_uniq(
	const CrackleHeader &header,
	const std::vector<unsigned char> &labels_binary
) {
	const uint64_t num_labels = decode_num_labels(header, labels_binary);
	std::vector<STORED_LABEL> uniq(num_labels);

	uint64_t idx = header.label_format == LabelFormat::FLAT
		? 8 // num labels
		: header.stored_data_width + 8; // bgcolor + numlabels for pins

	const unsigned char* buf = labels_binary.data();
	for (uint64_t i = 0; i < num_labels; i++, idx += sizeof(STORED_LABEL)) {
		uniq[i] = crackle::lib::ctoi<STORED_LABEL>(buf, idx);
	}

	return uniq;
}

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_flat(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary
) {
  std::vector<unsigned char> labels_binary = raw_labels(binary);
  const unsigned char* buf = labels_binary.data();

  const uint64_t num_labels = decode_num_labels(header, labels_binary);
  std::vector<STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

  const int cc_label_width = crackle::lib::compute_byte_width(num_labels);
  uint64_t offset = 8 + sizeof(STORED_LABEL) * num_labels + 4 * header.sz;

  uint64_t num_fields = (labels_binary.size() - offset) / cc_label_width;
  std::vector<LABEL> label_map(num_fields);

  for (uint64_t i = 0, j = offset; i < num_fields; i++, j += cc_label_width) {
  	if (cc_label_width == 1) {
		label_map[i] = static_cast<LABEL>(
			uniq[crackle::lib::ctoi<uint8_t>(buf, j)]
		);
	}
	else if (cc_label_width == 2) {
		label_map[i] = static_cast<LABEL>(
			uniq[crackle::lib::ctoi<uint16_t>(buf, j)]
		);
	}
	else if (cc_label_width == 4) {
		label_map[i] = static_cast<LABEL>(
			uniq[crackle::lib::ctoi<uint32_t>(buf, j)]
		);
	}
	else {
		label_map[i] = static_cast<LABEL>(
			uniq[crackle::lib::ctoi<uint64_t>(buf, j)]
		);
	}
  }
  return label_map;
}

template <typename LABEL, typename INDEX, typename DEPTH>
struct Pin {
	LABEL label;
	INDEX index;
	DEPTH depth;

	Pin() {
		label = 0;
		index = 0;
		depth = 0;
	}

	Pin(LABEL _lbl, INDEX _idx, DEPTH _depth) 
		: label(_lbl), index(_idx), depth(_depth)
	{} 
	
	uint64_t decode_buffer(const unsigned char* buf, const uint64_t idx) {
		label = crackle::lib::ctoi<LABEL>(buf, idx);
		index = crackle::lib::ctoi<INDEX>(buf, idx + sizeof(LABEL));
		depth = crackle::lib::ctoi<DEPTH>(buf, idx + sizeof(LABEL) + sizeof(INDEX));
		return sizeof(LABEL) + sizeof(INDEX) + sizeof(DEPTH);
	}
};

template <
	typename LABEL, typename STORED_LABEL, 
	typename RENUM_LABEL, typename INDEX, 
	typename DEPTH
>
std::vector<LABEL> decode_pins_helper3(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &labels_binary,
	const std::vector<STORED_LABEL> &uniq,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N,
	const LABEL bgcolor
) {
	typedef Pin<RENUM_LABEL, INDEX, DEPTH> PinType;
	const unsigned char* buf = labels_binary.data();

	const uint64_t pin_size = sizeof(RENUM_LABEL) + sizeof(INDEX) + sizeof(DEPTH);

	uint64_t offset = 8 + sizeof(STORED_LABEL) * (uniq.size() + 1);
	const uint64_t num_pins = (labels_binary.size() - offset) / pin_size;

	std::vector<PinType> pins(num_pins);
	for (uint64_t i = 0, j = offset; i < num_pins; i++) {
		j += pins[i].decode_buffer(buf, j);
	}

	const uint64_t sx = header.sx;
	const uint64_t sy = header.sy;

	const uint64_t sxy = sx * sy;

	std::vector<LABEL> label_map(N, bgcolor);
	for (uint64_t i = 0; i < num_pins; i++) {
		PinType pin = pins[i];
		for (uint64_t z = 0; z <= pin.depth; z++) {
			auto cc_id = cc_labels[pin.index + sxy * z];
			label_map[cc_id] = uniq[pin.label];
		}
	}

	return label_map;
}

template <typename LABEL, typename STORED_LABEL, typename RENUM_LABEL, typename INDEX>
std::vector<LABEL> decode_pins_helper2(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &labels_binary,
	const std::vector<STORED_LABEL> &uniq,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N,
	const LABEL bgcolor
) {
	int depth = header.depth_width();
	if (depth == 1) {
		return decode_pins_helper3<LABEL, STORED_LABEL, RENUM_LABEL, INDEX, uint8_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else if (depth == 2) {
		return decode_pins_helper3<LABEL, STORED_LABEL, RENUM_LABEL, INDEX, uint16_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else if (depth == 4) {
		return decode_pins_helper3<LABEL, STORED_LABEL, RENUM_LABEL, INDEX, uint32_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else {
		return decode_pins_helper3<LABEL, STORED_LABEL, RENUM_LABEL, INDEX, uint64_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
}

template <typename LABEL, typename STORED_LABEL, typename RENUM_LABEL>
std::vector<LABEL> decode_pins_helper(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &labels_binary,
	const std::vector<STORED_LABEL> &uniq,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N,
	const LABEL bgcolor
) {
	int width = header.pin_index_width();
	if (width == 1) {
		return decode_pins_helper2<LABEL, STORED_LABEL, RENUM_LABEL, uint8_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else if (width == 2) {
		return decode_pins_helper2<LABEL, STORED_LABEL, RENUM_LABEL, uint16_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else if (width == 4) {
		return decode_pins_helper2<LABEL, STORED_LABEL, RENUM_LABEL, uint32_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
	else {
		return decode_pins_helper2<LABEL, STORED_LABEL, RENUM_LABEL, uint64_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
	}
}

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_fixed_width_pins(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N
) {
  std::vector<unsigned char> labels_binary = raw_labels(binary);
  const LABEL bgcolor = static_cast<LABEL>(
  	crackle::lib::ctoi<STORED_LABEL>(
  		labels_binary.data(), 0
  	)
  );
  const uint64_t num_labels = decode_num_labels(header, labels_binary);
  std::vector<STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

  // bgcolor, num labels (u64), N labels, pins
  const int renum_width = crackle::lib::compute_byte_width(num_labels);

  if (renum_width == 1) {
		return decode_pins_helper<LABEL, STORED_LABEL, uint8_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
  }
  else if (renum_width == 2) {
		return decode_pins_helper<LABEL, STORED_LABEL, uint16_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
  }
  else if (renum_width == 4) {
		return decode_pins_helper<LABEL, STORED_LABEL, uint32_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
  }
  else {
		return decode_pins_helper<LABEL, STORED_LABEL, uint64_t>(
			header, labels_binary, uniq, cc_labels, N, bgcolor
		);
  }
}

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_label_map(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N
) {
	if (header.label_format == LabelFormat::FLAT) {
		return decode_flat<LABEL, STORED_LABEL>(header, binary);
	}
	else if (header.label_format == LabelFormat::PINS_FIXED_WIDTH) {
		return decode_fixed_width_pins<LABEL, STORED_LABEL>(
			header, binary, cc_labels, N
		);
	}
	else {
		std::string err = "crackle: Unsupported label format. Got: ";
		err += std::to_string(header.label_format);
		throw std::runtime_error(err);
	}
}

};
};

#endif