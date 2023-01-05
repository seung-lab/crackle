#ifndef __CRACKLE_LABELS_HXX__
#define __CRACKLE_LABELS_HXX__

#include <vector>

#include "robin_hood.hpp"

#include "header.hpp"
#include "lib.hpp"
#include "pins.hpp"

namespace crackle {
namespace labels {

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> encode_flat(
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

	std::vector<STORED_LABEL> mapping(N);

	uint32_t last = cc_labels[0];
	mapping[cc_labels[0]] = labels[0];
	for (int64_t i = 1; i < voxels; i++) {
		if (cc_labels[i] != last) {
			mapping[cc_labels[i]] = labels[i];
			last = cc_labels[i];
		}
	}

	cc_labels.release();

	std::vector<STORED_LABEL> uniq(mapping.begin(), mapping.end());
	std::sort(uniq.begin(), uniq.end());

	uint64_t j = 1;
	for (uint64_t i = 1; i < N; i++) {
		if (uniq[i] != uniq[i-1]) {
			uniq[j] = uniq[i];
			j++;
		}
	}
	uniq.resize(j);

	robin_hood::unordered_flat_map<STORED_LABEL, STORED_LABEL> remapping;
	for (STORED_LABEL i = 0; i < uniq.size(); i++) {
		remapping[uniq[i]] = i;
	}

	std::vector<STORED_LABEL> stored_labels(N);

	for (uint64_t i = 0; i < N; i++) {
		stored_labels[i] = remapping[mapping[i]];
	}

	int key_width = crackle::lib::compute_byte_width(uniq.size());
	int component_width = crackle::lib::compute_byte_width(sx * sy);

	std::vector<unsigned char> binary(
		8 
		+ sizeof(STORED_LABEL) * uniq.size() 
		+ component_width * num_components_per_slice.size()
		+ key_width * stored_labels.size()
	);

	int64_t i = 0;
	i += crackle::lib::itoc(
		static_cast<uint64_t>(uniq.size()), binary, i
	);
	for (auto val : uniq) {
		i += crackle::lib::itoc(
			static_cast<STORED_LABEL>(val), binary, i
		);		
	}
	for (auto val : num_components_per_slice) {
		i += crackle::lib::itocd(
			val, binary, i, component_width
		);		
	}

	for (auto val : stored_labels) {
		i += crackle::lib::itocd(
			val, binary, i, key_width
		);		
	}

	return binary;
}

template <typename STORED_LABEL>
STORED_LABEL find_bgcolor(
	std::unordered_map<uint64_t, std::vector<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>>& all_pins,
	const int64_t sz
) {
	// find bg color, pick the most pins
	// first, and then the one with the most
	// pin depth (so less decoding work later)
	STORED_LABEL bgcolor = 0;
	uint64_t max_pins = 0;
	uint64_t max_pins_depth = sz;
	for (auto& [label, pins] : all_pins) {
		if (pins.size() > max_pins) {
			bgcolor = static_cast<STORED_LABEL>(label);
			max_pins = pins.size();
			max_pins_depth = 0;
			for (auto& pin : pins) {
				max_pins_depth += pin.depth;
			}
		} 
		else if (pins.size() == max_pins) {
			uint64_t candidate_max_depth = 0;
			for (auto& pin : pins) {
				candidate_max_depth += pin.depth;
			}
			if (candidate_max_depth > max_pins_depth) {
				bgcolor = static_cast<STORED_LABEL>(label);
				max_pins_depth = candidate_max_depth;
			}
		}
	}

	return bgcolor;
}

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> encode_fixed_width_pins(
	std::unordered_map<uint64_t, std::vector<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>>& all_pins,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const int64_t index_width, const int64_t z_width
) {

	STORED_LABEL bgcolor = find_bgcolor<STORED_LABEL>(all_pins, sz);
	all_pins.erase(bgcolor);

	std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> linear;
	for (auto& [label, pins] : all_pins) {
		for (auto& pin : pins) {
			linear.emplace_back(
				pin.label, pin.index, pin.depth
			);
		}
	}

	struct {
		bool operator()(
			std::tuple<uint64_t,uint64_t,uint64_t>& a, std::tuple<uint64_t,uint64_t,uint64_t>& b
		) const { 
			return std::get<1>(a) < std::get<1>(b); 
		}
	} CmpIndex;

	std::sort(linear.begin(), linear.end(), CmpIndex);

	std::vector<STORED_LABEL> all_labels;
	all_labels.reserve(all_pins.size());
	for (auto& [label, pins] : all_pins) {
		all_labels.push_back(label);
	}
	std::sort(all_labels.begin(), all_labels.end());

	robin_hood::unordered_flat_map<STORED_LABEL, STORED_LABEL> renumbering;
	renumbering.reserve(all_labels.size());
	for (uint64_t i = 0; i < all_labels.size(); i++) {
		renumbering[all_labels[i]] = static_cast<STORED_LABEL>(i);
	}

	int renum_data_width = crackle::lib::compute_byte_width(all_labels.size());

	std::vector<unsigned char> binary(
		sizeof(STORED_LABEL) // bgcolor
		+ 8 // num labels
		+ sizeof(STORED_LABEL) * all_labels.size()
		+ (renum_data_width + index_width + z_width) * linear.size()
	);

	int64_t i = 0;
	i += crackle::lib::itoc(bgcolor, binary, i);
	i += crackle::lib::itoc(static_cast<uint64_t>(all_labels.size()), binary, i);
	for (auto label : all_labels) {
		i += crackle::lib::itoc(label, binary, i);
	}

	for (auto& pin : linear) {
		i += crackle::lib::itocd(renumbering[std::get<0>(pin)], binary, i, renum_data_width);
		i += crackle::lib::itocd(std::get<1>(pin), binary, i, index_width);
		i += crackle::lib::itocd(std::get<2>(pin), binary, i, z_width);
	}

	return binary;
}

template <typename LABEL, typename STORED_LABEL>
std::vector<unsigned char> encode_condensed_pins(
	std::unordered_map<uint64_t, std::vector<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>>& all_pins,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const int64_t index_width
) {
	STORED_LABEL bgcolor = find_bgcolor<STORED_LABEL>(all_pins, sz);
	all_pins.erase(bgcolor);

	uint64_t max_pins = 0;
	uint64_t max_depth = 0;
	uint64_t total_pins = 0;
	for (auto& [label, pins] : all_pins) {
		max_pins = std::max(static_cast<uint64_t>(pins.size()), max_pins);
		total_pins += pins.size();
		for (auto& pin : pins) {
			max_depth = std::max(max_depth, pin.depth);
		}
	}

	std::vector<STORED_LABEL> all_labels;
	all_labels.reserve(all_pins.size());
	for (auto& [label, pins] : all_pins) {
		all_labels.push_back(label);
	}
	std::sort(all_labels.begin(), all_labels.end());

	uint8_t num_pins_width = crackle::lib::compute_byte_width(max_pins);
	uint8_t depth_width = crackle::lib::compute_byte_width(max_depth);

	uint8_t combined_width = static_cast<uint8_t>(log2(num_pins_width)) | (static_cast<uint8_t>(log2(depth_width)) << 2);

	struct {
		bool operator()(
			crackle::pins::Pin<uint64_t, uint64_t, uint64_t>& a, 
			crackle::pins::Pin<uint64_t, uint64_t, uint64_t>& b
		) const { 
			return a.index < b.index; 
		}
	} CmpIndex;

	std::vector<unsigned char> binary(
		sizeof(STORED_LABEL) // bgcolor
		+ 8 // num labels
		+ sizeof(STORED_LABEL) * all_labels.size() // unique
		+ 1 // depth size, num_pins_size
		+ (num_pins_width * all_labels.size())
		+ (index_width + depth_width) * total_pins
	);

	int64_t i = 0;
	i += crackle::lib::itoc(bgcolor, binary, i);
	i += crackle::lib::itoc(static_cast<uint64_t>(all_labels.size()), binary, i);
	for (auto label : all_labels) {
		i += crackle::lib::itoc(label, binary, i); // STORED_LABEL size
	}
	i += crackle::lib::itoc(combined_width, binary, i);

	for (uint64_t label = 0; label < all_labels.size(); label++) {
		auto& pins = all_pins[all_labels[label]];
		std::sort(pins.begin(), pins.end(), CmpIndex);
		if (pins.size() > 1) {
			for (uint64_t j = pins.size() - 1; j >= 1; j--) {
				pins[j].index -= pins[j-1].index;
			}
		}

		i += crackle::lib::itocd(pins.size(), binary, i, num_pins_width);
		for (auto& pin : pins) {
			i += crackle::lib::itocd(pin.index, binary, i, index_width);
		}
		for (auto& pin : pins) {
			i += crackle::lib::itocd(pin.depth, binary, i, depth_width);
		}
	}

	return binary;
}


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

std::tuple<
	std::vector<uint64_t>,
	uint64_t,
	uint64_t
>
decode_components(
	const crackle::CrackleHeader &header,
	const unsigned char *buf,
	const uint64_t offset,
	const uint64_t num_grids, 
	const uint64_t component_width,
	const uint64_t z_start,
	const uint64_t z_end
) {
	std::vector<uint64_t> components(num_grids);
	for (uint64_t i = 0, j = offset; i < num_grids; i++, j += component_width) {
		components[i] = crackle::lib::ctoid(buf, j, component_width);
	}
	uint64_t component_left_offset = 0;
	uint64_t component_right_offset = 0;
	for (uint64_t z = 0; z < z_start; z++) {
		component_left_offset += components[z];
	}
	for (uint64_t z = header.sz - 1; z >= z_end; z--) {
		component_right_offset += components[z];
	}
	return std::make_tuple(components, component_left_offset, component_right_offset);
}

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_flat(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<unsigned char> labels_binary = raw_labels(binary);
	const unsigned char* buf = labels_binary.data();

	const uint64_t num_labels = decode_num_labels(header, labels_binary);
	std::vector<STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

	const int cc_label_width = crackle::lib::compute_byte_width(num_labels);

	const uint64_t num_grids = header.num_grids();
	uint64_t component_width = crackle::lib::compute_byte_width(header.sx * header.sy);

	uint64_t offset = 8 + sizeof(STORED_LABEL) * num_labels;
	auto [components, component_left_offset, component_right_offset] = decode_components(
		header, labels_binary.data(), offset, num_grids, component_width,
		z_start, z_end
	);
	offset += component_width * num_grids + component_left_offset * cc_label_width;
	uint64_t num_fields = (
		labels_binary.size() 
		- offset 
		- (component_right_offset * cc_label_width)
	) / cc_label_width;
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
				uniq[crackle::lib::ctoid(buf, j, component_width)]
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
	const uint64_t renum_width = crackle::lib::compute_byte_width(num_labels);
	const uint64_t index_width = header.pin_index_width();
	const uint64_t depth_width = header.depth_width();

	typedef crackle::pins::Pin<uint64_t, uint64_t, uint64_t> PinType;
	const unsigned char* buf = labels_binary.data();

	const uint64_t pin_size = renum_width + index_width + depth_width;

	uint64_t offset = 8 + sizeof(STORED_LABEL) * (uniq.size() + 1);
	const uint64_t num_pins = (labels_binary.size() - offset) / pin_size;

	std::vector<PinType> pins(num_pins);
	for (uint64_t i = 0, j = offset; i < num_pins; i++) {
		j += pins[i].dynamic_decode_buffer(
			buf, j,
			renum_width, index_width, depth_width
		);
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

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_condensed_pins(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N, 
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<unsigned char> labels_binary = raw_labels(binary);
	const LABEL bgcolor = static_cast<LABEL>(
		crackle::lib::ctoi<STORED_LABEL>(
			labels_binary.data(), 0
		)
	);
	std::vector<STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

	// bgcolor, num labels (u64), N labels, fmt depth num_pins, 
	// [num_pins][idx_1][depth_1]...[idx_n][depth_n]
	const uint64_t index_width = header.pin_index_width();

	typedef crackle::pins::Pin<uint64_t, int64_t, int64_t> PinType;
	const unsigned char* buf = labels_binary.data();

	uint64_t offset = 8 + sizeof(STORED_LABEL) * (uniq.size() + 1);
	uint8_t combined_width = crackle::lib::ctoi<uint8_t>(buf, offset);
	offset += 1;

	const uint8_t num_pins_width = pow(2, (combined_width & 0b11));
	const uint8_t depth_width = pow(2, (combined_width >> 2) & 0b11);

	std::vector<PinType> pins;
	for (uint64_t i = offset, label = 0; label < uniq.size(); label++) {
		uint64_t num_pins = crackle::lib::ctoid(buf, i, num_pins_width);
		i += num_pins_width;
		for (uint64_t j = 0; j < num_pins; j++) {
			uint64_t index = crackle::lib::ctoid(buf, i + (j * index_width), index_width);
			uint64_t depth = crackle::lib::ctoid(buf, i + (num_pins * index_width) + (j * depth_width), depth_width);
			pins.emplace_back(label, index, depth);
		}
		if (num_pins > 1) {
			for (uint64_t j = pins.size() - (num_pins-1); j < pins.size(); j++) {
				pins[j].index += pins[j-1].index;
			}
		}
		i += num_pins * (index_width + depth_width);
	}

	const int64_t sx = header.sx;
	const int64_t sy = header.sy;

	const int64_t sxy = sx * sy;

	std::vector<LABEL> label_map(N, bgcolor);
	for (auto& pin : pins) {
		int64_t pin_z = pin.index / sxy;
		int64_t loc = pin.index - (pin_z * sxy);
		int64_t pin_z_start = std::max(pin_z, static_cast<int64_t>(z_start));
		int64_t pin_z_end = pin_z + pin.depth + 1;
		pin_z_end = std::min(pin_z_end, static_cast<int64_t>(z_end));

		pin_z_start -= static_cast<int64_t>(z_start);
		pin_z_end -= static_cast<int64_t>(z_start);

		for (int64_t z = pin_z_start; z < pin_z_end; z++) {
			auto cc_id = cc_labels[loc + sxy * z];
			label_map[cc_id] = uniq[pin.label];
		}
	}

	return label_map;
}

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_label_map(
	const crackle::CrackleHeader &header,
	const std::vector<unsigned char> &binary,
	const std::vector<uint32_t> &cc_labels,
	const uint64_t N,
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<LABEL> label_map;
	if (header.label_format == LabelFormat::FLAT) {
		return decode_flat<LABEL, STORED_LABEL>(header, binary, z_start, z_end);
	}
	else if (header.label_format == LabelFormat::PINS_FIXED_WIDTH) {
		label_map = decode_fixed_width_pins<LABEL, STORED_LABEL>(
			header, binary, cc_labels, N
		);
	}
	else if (header.label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		label_map = decode_condensed_pins<LABEL, STORED_LABEL>(
			header, binary, cc_labels, N, z_start, z_end
		);
	}
	else {
		std::string err = "crackle: Unsupported label format. Got: ";
		err += std::to_string(header.label_format);
		throw std::runtime_error(err);
	}

	return label_map;
}

};
};

#endif