#ifndef __CRACKLE_PINS_HXX__
#define __CRACKLE_PINS_HXX__

#include <cstdint>
#include <unordered_map>

#include "robin_hood.hpp"

#include "cc3d.hpp"
#include "crc.hpp"
#include "lib.hpp"
#include "pairing_heap.hpp"

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

	uint64_t dynamic_decode_buffer(
		const unsigned char* buf, const uint64_t idx,
		const uint64_t label_width, const uint64_t index_width, const uint64_t depth_width  
	) {
		label = crackle::lib::ctoid(buf, idx, label_width);
		index = crackle::lib::ctoid(buf, idx + label_width, index_width);
		depth = crackle::lib::ctoid(buf, idx + label_width + index_width, depth_width);
		return label_width + index_width + depth_width;
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
		ccids.reserve(_ccids.size());
		ccids.insert(_ccids.begin(), _ccids.end());

	}

	uint64_t start_idx(uint64_t sx, uint64_t sy) const {
		return static_cast<uint64_t>(x) + sx * (static_cast<uint64_t>(y) + sy * static_cast<uint64_t>(z_s));
	}
	uint64_t depth() const {
		return static_cast<uint64_t>(z_e - z_s);
	}

	Pin<uint64_t,uint64_t,uint64_t> to_pin(
		const uint64_t label, const uint64_t sx, const uint64_t sy
	) const {
		return Pin<uint64_t,uint64_t,uint64_t>(
			label, start_idx(sx, sy), depth()
		);
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
		if (cc_labels[i] != last) {
			multiverse[labels[i]].emplace(cc_labels[i]);
			last = cc_labels[i];
		}
	}

	multiverse[labels[voxels - 1]].emplace(cc_labels[voxels - 1]);

	return multiverse;
}

void shrink_pin_to_fit(
	CandidatePin& pin,
	const std::unique_ptr<uint32_t[]> &cc_labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {
	const uint64_t sxy = sx * sy;
	const uint64_t voxels = sx * sy * sz;

	uint32_t min_id = cc_labels[voxels - 1];
	uint32_t max_id = 0;
	for (uint32_t ccid : pin.ccids) {
		min_id = std::min(min_id, ccid);
		max_id = std::max(max_id, ccid);
	}

	for (uint64_t z = pin.z_s; z <= pin.z_e; z++) {
		uint64_t idx = pin.x + sx * pin.y + sxy * z;
		if (cc_labels[idx] == min_id) {
			pin.z_s = z;
		}
		if (cc_labels[idx] == max_id) {
			pin.z_e = z;
			break;
		}
	}
}

std::vector<CandidatePin> find_optimal_pins(
	std::vector<CandidatePin> &pinsets,
	robin_hood::unordered_flat_set<uint32_t> &universe,
	const std::unique_ptr<uint32_t[]> &cc_labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {	
	std::vector<CandidatePin> final_pins;

	if (pinsets.size() == 0) {
		return final_pins;
	}

	final_pins.reserve(pinsets.size() / 10);

	crackle::pairing_heap::MinHeap<int64_t, int64_t> heap;

	std::vector<crackle::pairing_heap::PHNode<int64_t, int64_t>*> 
		heap_ptrs(pinsets.size());

	robin_hood::unordered_flat_set<int64_t> isets;
	isets.reserve(pinsets.size());

	for (int64_t i = 0; i < static_cast<int64_t>(pinsets.size()); i++) {
		heap_ptrs[i] = heap.emplace(-1 * pinsets[i].ccids.size(), i);
		isets.emplace(i);
	}

	while (universe.size()) {
		int64_t idx = heap.min_value();
		heap.pop();
		isets.erase(idx);

		CandidatePin& cur = pinsets[idx];
		for (uint32_t ccid : cur.ccids) {
			universe.erase(ccid);
		}

		shrink_pin_to_fit(cur, cc_labels, sx, sy, sz);

		if (universe.size() == 0) {
			final_pins.emplace_back(cur);
			break;
		}

		std::vector<uint64_t> to_erase;
		for (auto i : isets) {
			if (pinsets[i].z_s > cur.z_e || pinsets[i].z_e < cur.z_s) {
				continue;
			}

			auto& tmp = pinsets[i].ccids;
			for (uint32_t ccid : cur.ccids) {
				tmp.erase(ccid);
			}

			if (tmp.size() == 0) {
				to_erase.push_back(i);
				heap.erase(heap_ptrs[i]);
			}
			else {
				heap.update_key(heap_ptrs[i], -1 * tmp.size());
			}
		}
		for (uint64_t i : to_erase) {
			isets.erase(i);
		}

		final_pins.emplace_back(cur);
	}

	return final_pins;
}

std::vector<CandidatePin> find_suboptimal_pins(
	std::vector<CandidatePin> &pinsets,
	robin_hood::unordered_flat_set<uint32_t> &universe,
	const std::unique_ptr<uint32_t[]> &cc_labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz
) {	
	std::vector<CandidatePin> final_pins;

	if (pinsets.size() == 0) {
		return final_pins;
	}

	final_pins.reserve(pinsets.size() / 10);

	robin_hood::unordered_node_map<int64_t, std::vector<int64_t>> component_to_pins;
	component_to_pins.reserve(universe.size());

	for (uint64_t i = 0; i < pinsets.size(); i++) {
		CandidatePin& pin = pinsets[i];
		for (uint64_t label : pin.ccids) {
			component_to_pins[label].push_back(i);
		}
	}

	while (universe.size()) {
		uint32_t picked_ccid = *universe.begin();
		auto& pins = component_to_pins[picked_ccid];

		CandidatePin max_pin = pinsets[pins[0]];
		int max_depth = max_pin.z_e - max_pin.z_s;
		for (uint64_t i = 1; i < pins.size(); i++) {
			CandidatePin cur = pinsets[pins[i]];
			int depth = cur.z_e - cur.z_s;
			if (depth > max_depth) {
				max_pin = cur;
			}
		}

		for (uint32_t ccid : max_pin.ccids) {
			universe.erase(ccid);
		}

		final_pins.push_back(max_pin);
	}

	return final_pins;
}

template <typename LABEL>
std::tuple<
	std::unordered_map<uint64_t, std::vector<CandidatePin>>,
	std::vector<uint64_t>,
	uint64_t,
	std::vector<uint32_t>
>
compute(
	const LABEL* labels,
	const uint64_t sx, const uint64_t sy, const uint64_t sz,
	const bool optimize
) {
	const uint64_t sxy = sx * sy;

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
	std::unordered_map<uint64_t, std::vector<CandidatePin>> all_pins;
	all_pins.reserve(128);

	auto multiverse = compute_multiverse<LABEL>(
		labels, cc_labels.get(), sx, sy, sz, N_total
	);

	auto find_pins_fn = optimize ? find_optimal_pins : find_suboptimal_pins;

	for (auto [label, pins] : pinsets) {
		all_pins[label] = find_pins_fn(
			pins, multiverse[label], cc_labels, 
			sx, sy, sz
		);
	}

	std::vector<uint32_t> crcs(sz);
	uint64_t tmp_N = 0;
	for (uint64_t z = 0; z < sz; z++) {
		if (tmp_N > 0) {
			for (uint64_t i = 0; i < sxy; i++) {
				cc_labels[sxy * z + i] -= tmp_N;
			}
		}
		crcs[z] = crackle::crc::crc32c(cc_labels.get() + sxy * z, sxy);
		tmp_N += num_components_per_slice[z];
	}

	return std::make_tuple(all_pins, num_components_per_slice, N_total, crcs);
}

};
};

#endif