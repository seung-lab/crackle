#ifndef __CRACKLE_PINS_HXX__
#define __CRACKLE_PINS_HXX__

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "robin_hood.hpp"

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
	robin_hood::unordered_flat_set<uint32_t> ccids;

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
		const std::vector<uint32_t> &_ccids) 
		: x(_x), y(_y), z_s(_z_s), 
		  z_e(_z_e)
	{
		ccids.insert(_ccids.begin(), _ccids.end());

	}

	uint64_t start_idx(uint64_t sx, uint64_t sy) {
		return static_cast<uint64_t>(x) + sx * (static_cast<uint64_t>(y) + sy * static_cast<uint64_t>(z_s));
	}
	uint64_t depth() {
		return static_cast<uint64_t>(z_e - z_s);
	}
};

template <typename LABEL>
void add_pin(
	robin_hood::unordered_node_map<LABEL, std::vector<CandidatePin>> &pinsets,
	LABEL label,
	const uint64_t z_start,
	const uint64_t x, const uint64_t y, const uint64_t z,
	const std::vector<uint32_t> &cc_set
) {
	if (pinsets[label].size() == 0) {
		pinsets[label].emplace_back(x, y, z_start, z, cc_set );
		return;
	}

	CandidatePin last_pin = pinsets[label].back();

	if (last_pin.x == x - 1 && last_pin.y == y) {
		if (last_pin.z_s <= z_start && last_pin.z_e >= z) {
			return;
		}
		else if (last_pin.z_s >= z_start && last_pin.z_e <= z) {
			pinsets[label].back() = CandidatePin(x, y, z_start, z, cc_set);
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
robin_hood::unordered_node_map<LABEL, std::vector<CandidatePin>> extract_columns(
	const LABEL* labels,
	const uint32_t* cc_labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {
	robin_hood::unordered_node_map<LABEL, std::vector<CandidatePin>> pinsets;

	for (uint64_t y = 0; y < sy; y++) {
		for (uint64_t x = 0; x < sx; x++) {
			uint64_t loc = x + sx * y;
			LABEL label = labels[loc];

			std::vector<uint32_t> label_set;
			label_set.reserve(sz);
			label_set.push_back(cc_labels[loc]);
			uint64_t z_start = 0;
			uint64_t z = 1;
			for (; z < sz; z++) {
				uint64_t zoff = sx * sy * z;
				LABEL cur = labels[loc + zoff];
				if (label != cur) {
					add_pin(pinsets, label, z_start, x, y, z-1, label_set);
					label = cur;
					z_start = z;
					label_set.clear();
				}
				label_set.push_back(cc_labels[loc + zoff]);
			}
			if (sz == 1) {
				z = 0;
			}
			add_pin(pinsets, label, z_start, x, y, z-1, label_set);
		}
	}

	return pinsets;
}

template <typename LABEL>
std::unordered_map<LABEL, robin_hood::unordered_flat_set<uint32_t>> 
compute_multiverse(
	const LABEL* labels,
	const uint32_t* cc_labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz,
	const uint64_t N
) {

	std::unordered_map<
		LABEL, robin_hood::unordered_flat_set<uint32_t>
	> multiverse;
	multiverse.reserve(N);

	const uint64_t voxels = sx * sy * sz;

	if (voxels == 0) {
		return multiverse;
	}

	uint32_t last = cc_labels[0];
	multiverse[labels[0]].emplace(cc_labels[0]);

	for (uint64_t i = 1; i < voxels; i++) {
		if (labels[i-1] != last) {
			multiverse[labels[i]].emplace(cc_labels[i]);
			last = cc_labels[i];
		}
	}

	multiverse[last].emplace(cc_labels[voxels - 1]);

	return multiverse;
}


std::vector<CandidatePin> find_optimal_pins(
	std::vector<CandidatePin> &pinsets,
	robin_hood::unordered_flat_set<uint32_t> &universe,
	const uint64_t sx, const uint64_t sy
) {	
	std::vector<CandidatePin> final_pins;
	final_pins.reserve(final_pins.size() / 10);

	std::unordered_set<uint64_t> isets;
	isets.reserve(pinsets.size());
	for (uint64_t i = 0; i < pinsets.size(); i++) {
		isets.emplace(i);
	}

	std::vector<uint64_t> sizes;
	sizes.reserve(pinsets.size());
	while (universe.size()) {
		uint64_t idx = 0;
		uint64_t ct = 0;
		for (auto i : isets) {
			if (pinsets[i].ccids.size() > ct) {
				ct = pinsets[i].ccids.size();
				idx = i;
			}
		}

		isets.erase(idx);

		CandidatePin cur = pinsets[idx];
		for (uint32_t ccid : cur.ccids) {
			universe.erase(ccid);
		}

		if (universe.size() == 0) {
			final_pins.emplace_back(cur);
			break;
		}

		std::vector<uint64_t> to_erase;
		for (auto i : isets) {
			auto& tmp = pinsets[i].ccids;
			for (uint32_t ccid : cur.ccids) {
				tmp.erase(ccid);
			}

			if (tmp.size() == 0) {
				to_erase.push_back(i);
			}
		}
		for (uint64_t i : to_erase) {
			isets.erase(i);
		}

		final_pins.emplace_back(cur);
	}

	return final_pins;
}

template <typename LABEL>
std::unordered_map<uint64_t, std::vector<Pin<uint64_t, uint64_t, uint64_t>>> 
compute(
	const LABEL* labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {
	typedef Pin<uint64_t, uint64_t, uint64_t> PinType;

	std::vector<uint64_t> num_components_per_slice(sz);
	uint64_t N_total = 0;

	std::unique_ptr<uint32_t[]> cc_labels(
		crackle::cc3d::connected_components<LABEL, uint32_t>(
			labels, sx, sy, sz,
			num_components_per_slice,
			/*out=*/NULL, N_total
		)
	);

	auto pinsets = extract_columns(labels, cc_labels.get(), sx, sy, sz);
	std::unordered_map<uint64_t, std::vector<PinType>> all_pins;
	all_pins.reserve(128);

	auto multiverse = compute_multiverse<LABEL>(
		labels, cc_labels.get(), sx, sy, sz, N_total
	);

	for (auto [label, pins] : pinsets) {
		std::vector<CandidatePin> solution = find_optimal_pins(
			pins, multiverse[label], sx, sy
		);
		std::vector<PinType> encoded_pins;
		encoded_pins.reserve(solution.size());
		for (auto pin : solution) {
			encoded_pins.emplace_back(
				label, 
				pin.start_idx(sx, sy),
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