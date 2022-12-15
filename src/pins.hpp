#ifndef __CRACKLE_PINS_HXX__
#define __CRACKLE_PINS_HXX__

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "cc3d.hpp"
#include "lib.hpp"

namespace crackle {
namespace pins {

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

struct CandidatePin {
	uint32_t x;
	uint32_t y;
	uint32_t z_s;
	uint32_t z_e;
	std::unordered_set<uint32_t> ccids;

	CandidatePin() {
		x = 0;
		y = 0;
		z_s = 0;
		z_e = 0;
	}

	CandidatePin(
		const uint32_t _x,
		const uint32_t _y,
		const uint32_t _z_s,
		const uint32_t _z_e,
		std::unordered_set<uint32_t> _ccids) 
		: x(_x), y(_y), z_s(_z_s), 
		  z_e(_z_e), ccids(_ccids)
	{}

	uint64_t start_idx(uint64_t sx, uint64_t sy) {
		return static_cast<uint64_t>(x) + sx * (static_cast<uint64_t>(y) + sy * static_cast<uint64_t>(z_s));
	}
	uint64_t depth() {
		return static_cast<uint64_t>(z_e - z_s);
	}
};


template <typename LABEL>
void add_pin(
	std::unordered_map<LABEL, std::vector<CandidatePin>> &pinsets,
	LABEL label,
	const uint64_t z_start,
	const uint64_t x, const uint64_t y, const uint64_t z,
	const uint64_t sx, const uint64_t sy, const uint64_t sz,
	const std::unordered_set<uint32_t> &cc_set
) {
	if (pinsets[label].size() == 0) {
		pinsets[label].emplace_back(x, y, z_start, z, cc_set);
		return;
	}

	CandidatePin last_pin = pinsets[label].back();

	if (last_pin.x == x - 1 && last_pin.y == y) {
		if (last_pin.z_s <= z_start && last_pin.z_e >= z) {
			return;
		}
		else if (last_pin.z_s >= z_start && last_pin.z_e <= z) {
			pinsets[label].back() = CandidatePin(x, y, z_start, cc_set);
		}
		else {
			pinsets[label].emplace_back(x, y, z_start, z, cc_set);
		}
	}
	else {
		pinsets[label].emplace_back(x, y, z_start, z, cc_set);
	}
}

template <typename LABEL>
std::unordered_map<LABEL, std::vector<CandidatePin>> extract_columns(
	LABEL* labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {
	std::vector<uint64_t> num_components_per_slice(sz);
	uint64_t N_total = 0;

	uint32_t* cc_labels = crackle::cc3d::connected_components<uint32_t>(
		labels, sx, sy, sz, 
		num_components_per_slice,
		NULL, N_total
	);

	std::unordered_map<LABEL, std::vector<CandidatePin>> pinsets;

	for (uint64_t y = 0; y < sy; y++) {
		for (uint64_t x = 0; x < sx; x++) {
			uint64_t loc = x + sx * y;
			LABEL label = labels[loc];

			std::unordered_set<uint32_t> label_set;
			label_set.insert(cc_labels[loc]);
			uint64_t z_start = 0;
			for (uint64_t z = 0; z < sz; z++) {
				uint64_t zoff = sx * sy * z;
				LABEL cur = labels[loc + zoff];
				if (label != cur) {
					add_pin(pinsets, label, z_start, x, y, z-1, label_set);
					label = cur;
					z_start = z;
					label_set.clear();
				}
				label_set.insert(cc_labels[loc + zoff]);
			}
			add_pin(pinsets, label, z_start, x, y, z, label_set);
		}
	}

	return pinsets;
}

std::vector<CandidatePin> find_optimal_pins(
	std::vector<CandidatePin> &pinsets,
	const uint64_t sx, const uint64_t sy
) {	
	std::unordered_set<uint32_t> universe;
	for (auto candidate : pinsets) {
		universe.merge(candidate.ccids);
	}

	std::vector<CandidatePin> final_pins;
	final_pins.reserve(final_pins.size() / 10);

	std::unordered_set<uint64_t> isets(pinsets.size());
	for (uint64_t i = 0; i < isets.size(); i++) {
		isets.emplace(i);
	}

	std::vector<uint64_t> sizes;
	sizes.reserve(pinsets.size());
	while (universe.size()) {
		uint64_t idx = 0;
		uint64_t ct = 0;
		for (auto i : isets) {
			if (pinsets[i].ccids.size() > ct) {
				ct = pinsets[i].size();
				idx = i;
			}
		}
		uint64_t j = *isets.begin();
		isets.erase(j);

		CandidatePin cur = pinsets[j];
		for (uint32_t ccid : cur.ccids) {
			universe.erase(ccid);
		}

		for (auto i : isets) {
			auto& tmp = pinsets[i].ccids;
			for (uint32_t ccid : cur.ccids) {
				tmp.erase(ccid);
			}

			if (tmp.size() == 0) {
				isets.erase(i);
			}
		}

		final_pins.emplace_back(cur);
	}

	return final_pins;
}

template <typename LABEL, typename INDEX, typename DEPTH>
std::vector<Pin<LABEL, INDEX, DEPTH>> compute(
	LABEL* labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {
	typedef Pin<LABEL, INDEX, DEPTH> PinType;

	auto pinsets = extract_columns(labels, sx, sy, sz);
	std::unordered_map<LABEL, std::vector<PinType>> all_pins;
	for (auto [label, pins] : pinsets) {
		std::vector<CandidatePin> solution = find_optimal_pins(pins, sx, sy);
		std::vector<PinType> encoded_pins(solution.size());
		for (auto pin : solution) {
			encoded_pins.emplace_back(
				label, 
				pin.start_index(sx, sy),
				pin.depth()
			);
		}
		all_pins[label] = encoded_pins;
	}

	return all_pins;
}

};
};

#endif