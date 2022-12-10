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
	uint64_t startpt;
	uint64_t endpt;
	std::unordered_set<uint32_t> ccids;

	CandidatePin() {
		startpt = 0;
		endpt = 0;
	}

	CandidatePin(uint64_t _startpt, uint64_t _endpt, std::unordered_set<uint32_t> _ccids) 
		: startpt(_startpt), endpt(_endpt), ccids(_ccids)
	{}
};

template <typename LABEL>
std::vector<CandidatePin> extract_columns(
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
			for (uint64_t z = 0; z < sz; z++) {
				uint64_t zoff = sx * sy * z;
				LABEL cur = labels[loc + zoff];
				if (label != cur) {
					pinsets[label].emplace_back(loc, loc+zoff, label_set);
					label = cur;
					label_set.clear();
				}
				label_set.insert(cc_labels[loc + zoff]);
			}
		}
	}

	return pinsets;
}


};
};

#endif