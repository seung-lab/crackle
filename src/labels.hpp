#ifndef __CRACKLE_LABELS_HXX__
#define __CRACKLE_LABELS_HXX__

#include <span>
#include <vector>
#include <mutex>

#include "robin_hood.hpp"

#include "crc.hpp"
#include "header.hpp"
#include "lib.hpp"
#include "pins.hpp"
#include "threadpool.hpp"

namespace crackle {
namespace labels {

// For pin encodings only, extract the background color.
uint64_t background_color(std::span<unsigned char> binary) {
	crackle::CrackleHeader header(binary);

	if (header.label_format == LabelFormat::FLAT) {
		throw std::runtime_error("Background color can only be extracted from pin encoded streams.");
	}
	uint64_t offset = header.header_bytes() + header.grid_index_bytes();
	return crackle::lib::ctoid(binary, offset, header.stored_data_width);
}

template <typename LABEL, typename STORED_LABEL>
std::tuple<
	std::vector<unsigned char>,
	std::vector<uint32_t>
> encode_flat(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const size_t parallel
) {
	const int64_t sxy = sx * sy;

	std::vector<uint64_t> num_components_per_slice(sz);
	std::vector<uint32_t> crcs(sz);

	std::vector<std::vector<uint32_t>> cc_labels_scratch(parallel);
	std::vector<std::vector<STORED_LABEL>> mapping_scratch(sz);

	for (size_t t = 0; t < parallel; t++) {
		cc_labels_scratch[t].resize(sxy);
	}

	ThreadPool pool(parallel);

	std::mutex mtx;
	uint64_t N = 0;

	for (int64_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			std::vector<uint32_t>& cc_labels = cc_labels_scratch[t];
			std::vector<STORED_LABEL>& mapping = mapping_scratch[z];

			uint64_t tmp_N = 0;
			crackle::cc3d::connected_components2d_4<LABEL, uint32_t>(
				(labels + sxy * z), 
				sx, sy, 1, 
				cc_labels.data(),
				0, tmp_N
			);

			mapping.resize(tmp_N);

			uint64_t last = cc_labels[0];
			mapping[cc_labels[0]] = labels[sxy * z];
			for (int64_t i = 1; i < sxy; i++) {
				if (cc_labels[i] != last) {
					mapping[cc_labels[i]] = labels[sxy * z + i];
					last = cc_labels[i];
				}
			}

			num_components_per_slice[z] = tmp_N;
			crcs[z] = crackle::crc::crc32c(cc_labels.data(), sxy);
			
			std::unique_lock<std::mutex> lock(mtx);
			N += tmp_N;
		});
	}

	pool.join();

	cc_labels_scratch = std::vector<std::vector<uint32_t>>(); // deallocate memory

	std::vector<STORED_LABEL> mapping;
	mapping.reserve(N);
	for (int64_t z = 0; z < sz; z++) {
		mapping.insert(mapping.end(), mapping_scratch[z].begin(), mapping_scratch[z].end());
	}
	mapping_scratch = std::vector<std::vector<STORED_LABEL>>(); // deallocate memory

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
	remapping.reserve(uniq.size());
	for (uint64_t i = 0; i < uniq.size(); i++) {
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

	return std::make_tuple(binary, crcs);
}

template <typename STORED_LABEL>
STORED_LABEL find_bgcolor(
	std::unordered_map<uint64_t, std::vector<crackle::pins::CandidatePin>>& all_pins,
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
				max_pins_depth += pin.depth();
			}
		} 
		else if (pins.size() == max_pins) {
			uint64_t candidate_max_depth = 0;
			for (auto& pin : pins) {
				candidate_max_depth += pin.depth();
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
std::vector<unsigned char> encode_condensed_pins(
	std::unordered_map<uint64_t, std::vector<crackle::pins::CandidatePin>>& all_pins,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const int64_t index_width, 
	const std::vector<uint64_t>& num_components_per_slice,
	const uint64_t num_components,
	const bool auto_bgcolor = true,
	const STORED_LABEL manual_bgcolor = 0
) {
	STORED_LABEL bgcolor = manual_bgcolor;
	if (auto_bgcolor) {
		bgcolor = find_bgcolor<STORED_LABEL>(all_pins, sz);
	}

	all_pins.erase(bgcolor);

	uint64_t max_pins = 0;
	uint64_t max_depth = 0;
	uint64_t total_pins = 0;
	for (auto& [label, pins] : all_pins) {
		max_pins = std::max(static_cast<uint64_t>(pins.size()), max_pins);
		total_pins += pins.size();
		for (auto& pin : pins) {
			max_depth = std::max(max_depth, pin.depth());
		}
	}

	std::vector<STORED_LABEL> all_labels;
	all_labels.reserve(all_pins.size());
	for (auto& [label, pins] : all_pins) {
		all_labels.push_back(label);
	}
	std::sort(all_labels.begin(), all_labels.end());

	const uint8_t num_pins_width = crackle::lib::compute_byte_width(max_pins);
	const uint8_t depth_width = crackle::lib::compute_byte_width(max_depth);
	const uint8_t cc_label_width = crackle::lib::compute_byte_width(num_components);
	const uint8_t component_width = crackle::lib::compute_byte_width(sx * sy);

	const uint8_t pin_bytes = index_width + depth_width;
	const uint8_t cc_efficient_threshold = pin_bytes / cc_label_width;

	uint8_t combined_width = (
		static_cast<uint8_t>(log2(num_pins_width))
		| (static_cast<uint8_t>(log2(depth_width)) << 2)
		| (static_cast<uint8_t>(log2(cc_label_width)) << 4)
	);

	struct CmpIndex {
		uint64_t sx;
		uint64_t sy;

		CmpIndex(uint64_t _sx, uint64_t _sy) 
			: sx(_sx), sy(_sy)
		{}

		bool operator()(
			crackle::pins::CandidatePin& a, 
			crackle::pins::CandidatePin& b
		) const { 
			return a.start_idx(sx, sy) < b.start_idx(sx, sy); 
		}
	};

	const CmpIndex cmp(sx,sy);

	// overestimate size using more expensive pins
	// and then resize at end
	std::vector<unsigned char> binary(
		sizeof(STORED_LABEL) // bgcolor
		+ 8 // num labels
		+ sizeof(STORED_LABEL) * all_labels.size() // unique
		+ (component_width * num_components_per_slice.size())
		+ 1 // depth size, num_pins_size
		+ (2 * num_pins_width * all_labels.size())
		+ (index_width + depth_width) * total_pins
	);

	int64_t i = 0;
	i += crackle::lib::itoc(bgcolor, binary, i);
	i += crackle::lib::itoc(static_cast<uint64_t>(all_labels.size()), binary, i);
	for (auto label : all_labels) {
		i += crackle::lib::itoc(label, binary, i); // STORED_LABEL size
	}
	for (auto val : num_components_per_slice) {
		i += crackle::lib::itocd(
			val, binary, i, component_width
		);		
	}
	i += crackle::lib::itoc(combined_width, binary, i);

	for (uint64_t label = 0; label < all_labels.size(); label++) {
		auto& pins = all_pins[all_labels[label]];
		std::sort(pins.begin(), pins.end(), cmp);

		std::vector<uint64_t> pin_repr;
		std::vector<uint64_t> cc_repr;
		
		for (uint64_t j = 0; j < pins.size(); j++) {
			auto& pin = pins[j];
			if (pin.depth() < cc_efficient_threshold) {
				cc_repr.push_back(j);
			}
			else {
				pin_repr.push_back(j);
			}
		}

		std::vector<uint64_t> pin_index;
		pin_index.reserve(pins.size());
		for (uint64_t j : pin_repr) {
			pin_index.push_back(pins[j].start_idx(sx, sy));
		}

		if (pin_index.size() > 1) {
			for (uint64_t j = pin_index.size() - 1; j >= 1; j--) {
				pin_index[j] -= pin_index[j-1];
			}
		}

		i += crackle::lib::itocd(pin_repr.size(), binary, i, num_pins_width);
		for (uint64_t index : pin_index) {
			i += crackle::lib::itocd(index, binary, i, index_width);
		}
		for (uint64_t j : pin_repr) {
			i += crackle::lib::itocd(pins[j].depth(), binary, i, depth_width);
		}

		std::vector<uint32_t> cc_ids;
		cc_ids.reserve(cc_repr.size() * cc_efficient_threshold);
		for (uint64_t j : cc_repr) {
			auto& pin = pins[j];
			for (uint32_t ccid : pin.ccids) {
				cc_ids.push_back(ccid);
			}
		}
		std::sort(cc_ids.begin(), cc_ids.end());
		if (cc_ids.size() > 1) {
			for (uint64_t j = cc_ids.size() - 1; j >= 1; j--) {
				cc_ids[j] -= cc_ids[j-1];
			}
		}

		i += crackle::lib::itocd(cc_ids.size(), binary, i, num_pins_width);
		for (uint32_t ccid : cc_ids) {
			i += crackle::lib::itocd(ccid, binary, i, cc_label_width);
		}
	}

	binary.resize(i);
	return binary;
}


std::span<const unsigned char> raw_labels(
	const std::span<const unsigned char> &binary
) {
	crackle::CrackleHeader header(binary);
	return std::span<const unsigned char>(
		(binary.begin() + header.header_bytes() + header.grid_index_bytes()),
		header.num_label_bytes
	);
}

uint64_t decode_num_labels(
	const CrackleHeader &header,
	const std::span<const unsigned char> &labels_binary
) {
	if (header.label_format == LabelFormat::FLAT) {
		return crackle::lib::ctoi<uint64_t>(labels_binary.data(), 0);
	}
	else {
		return crackle::lib::ctoi<uint64_t>(labels_binary.data(), header.stored_data_width);
	}
}

uint64_t num_labels(const std::span<const unsigned char> &binary) {
	const CrackleHeader header(binary);
	std::span<const unsigned char> labels_binary = raw_labels(binary);
	return decode_num_labels(header, labels_binary);
}

template <typename STORED_LABEL>
std::span<const STORED_LABEL> decode_uniq(
	const CrackleHeader &header,
	const std::span<const unsigned char> &labels_binary
) {
	const uint64_t num_labels = decode_num_labels(header, labels_binary);

	uint64_t idx = header.label_format == LabelFormat::FLAT
		? 8 // num labels
		: header.stored_data_width + 8; // bgcolor + numlabels for pins

    const unsigned char* buf = labels_binary.data();
	return std::span<const STORED_LABEL>(
		 reinterpret_cast<const STORED_LABEL*>(buf + idx),
		num_labels
	);
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
	const std::span<const unsigned char> &binary,
	const uint64_t z_start, const uint64_t z_end
) {
	std::span<const unsigned char> labels_binary = raw_labels(binary);
	const unsigned char* buf = labels_binary.data();

	const uint64_t num_labels = decode_num_labels(header, labels_binary);
	std::span<const STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

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

template <typename LABEL, typename STORED_LABEL>
std::vector<LABEL> decode_condensed_pins(
	const crackle::CrackleHeader &header,
	const std::span<const unsigned char> &binary,
	const uint32_t* cc_labels,
	const uint64_t N, 
	const uint64_t z_start, const uint64_t z_end
) {
	std::span<const unsigned char> labels_binary = raw_labels(binary);
	const LABEL bgcolor = static_cast<LABEL>(
		crackle::lib::ctoi<STORED_LABEL>(
			labels_binary.data(), 0
		)
	);
	std::span<const STORED_LABEL> uniq = decode_uniq<STORED_LABEL>(header, labels_binary);

	// bgcolor, num labels (u64), N labels, fmt depth num_pins, 
	// [num_pins][idx_1][depth_1]...[idx_n][depth_n][num_cc][cc_1][cc_2]...[cc_n]
	const uint64_t index_width = header.pin_index_width();
	const uint64_t component_width = crackle::lib::compute_byte_width(header.sx * header.sy);

	typedef crackle::pins::Pin<uint64_t, int64_t, int64_t> PinType;
	const unsigned char* buf = labels_binary.data();

	uint64_t offset = 8 + sizeof(STORED_LABEL) * (uniq.size() + 1);

	auto [components, component_left_offset, component_right_offset] = decode_components(
		header, labels_binary.data(), offset, header.num_grids(), component_width,
		z_start, z_end
	);

	uint64_t N_all = 0;
	for (uint64_t j = 0; j < components.size(); j++) {
		N_all += components[j];
	}

	component_right_offset = N_all - component_right_offset;
	offset += component_width * header.num_grids();

	uint8_t combined_width = crackle::lib::ctoi<uint8_t>(buf, offset);
	offset += 1;

	const uint8_t num_pins_width = pow(2, (combined_width & 0b11));
	const uint8_t depth_width = pow(2, (combined_width >> 2) & 0b11);
	const uint8_t cc_label_width = pow(2, (combined_width >> 4) & 0b11);

	std::vector<LABEL> label_map(N, bgcolor);

	std::vector<PinType> pins;
	
	for (uint64_t i = offset, label = 0; label < uniq.size(); label++) {
		if (i >= labels_binary.size()) {
			throw std::runtime_error("crackle: pin section is malformed or corrupted.");
		}

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

		uint64_t num_cc_labels = crackle::lib::ctoid(buf, i, num_pins_width);
		i += num_pins_width;
		std::vector<uint32_t> cc_labels(num_cc_labels);
		for (uint64_t j = 0; j < num_cc_labels; j++) {
			cc_labels[j] = crackle::lib::ctoid(buf, i, cc_label_width);
			i += cc_label_width;
		}
		for (uint64_t j = 1; j < num_cc_labels; j++) {
			cc_labels[j] += cc_labels[j-1];
		}
		for (uint64_t j = 0; j < num_cc_labels; j++) {
			if (cc_labels[j] < component_left_offset || cc_labels[j] >= component_right_offset) {
				continue;
			}
			label_map[cc_labels[j] - component_left_offset] = uniq[label];
		}
	}

	const int64_t sx = header.sx;
	const int64_t sy = header.sy;
	const int64_t sxy = sx * sy;

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
	const std::span<const unsigned char> &binary,
	const uint32_t* cc_labels,
	const uint64_t N,
	const uint64_t z_start, const uint64_t z_end
) {
	std::vector<LABEL> label_map;
	if (header.label_format == LabelFormat::FLAT) {
		return decode_flat<LABEL, STORED_LABEL>(header, binary, z_start, z_end);
	}
	else if (header.label_format == LabelFormat::PINS_VARIABLE_WIDTH) {
		if (cc_labels == NULL) {
			std::string err = "crackle: cc_labels must not be null.";
			throw std::runtime_error(err);
		}
		
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