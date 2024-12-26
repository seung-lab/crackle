#ifndef __CRACKLE_DUALGRAPH_HXX__
#define __CRACKLE_DUALGRAPH_HXX__

#include <vector>
#include <cstdint>
#include <algorithm>

#include "crackcodes.hpp"
#include "builtins.hpp"
#include "cc3d.hpp"

namespace crackle {
namespace dual_graph {

const uint8_t VISITED_BIT = 0b10000;

// voxel connectivity matches cc3d_graphs.hpp 4 connected
// four bits: -y+y-x+x true is passable
enum VCGDirectionCode {
	NONE = 0b0000,
	LEFT = 0b0010,
	RIGHT = 0b0001,
	UP = 0b1000,
	DOWN = 0b0100,
	ANY = 0b1111
};

void print_bits(uint8_t val) {
	printf("vcg: 0b");
	for (int i = 4; i >= 0; i--) {
		printf("%d", (val >> i) & 0b1);
	}
}

struct VCGGraph {
	std::vector<uint8_t>& vcg;
	
	int64_t sx;
	int64_t sy;

	VCGGraph(std::vector<uint8_t>& _vcg, int64_t _sx, int64_t _sy) 
		: vcg(_vcg), sx(_sx), sy(_sy) {}

	bool next_contour(uint32_t& barriers, int64_t& idx, int64_t& y) {
		int64_t x = idx - sx * y;

		// important for barrier to be after
		// if even numbers of barriers are contours and odd numbers
		// are holes, assume the first contour starting at 0,0 is not
		// a hole. If two contours are separated by more than one space,
		// there will be a contour and a hole and then a contour. If
		// they are separated by only one space, there will be no hole
		// contour (because the contour and hole are the same).
		// So add the sum of the left and right barriers, but after
		// b/c otherwise the first contour will frequently be considered
		// a hole and that throws everything off.

		for (; y < sy; y++) {
			for (; x < sx; x++, idx++) {
				// condensing this conditional seems to save 5% in one speed test
				// if (((vcg[idx] & 0b11) < 0b11) && (vcg[idx] & VISITED_BIT) == 0) {
				if ((vcg[idx] & 0b110011) < 0b11 || vcg[idx] == 0b11100) {
					return true;
				}
				barriers += static_cast<uint32_t>(popcount((~vcg[idx]) & 0b11));
			}
			barriers = 0;
			x = 0;
		}

		return false;
	}
};

#define TRY_LEFT if (allowed_dirs & VCGDirectionCode::LEFT) {return VCGDirectionCode::LEFT;}
#define TRY_RIGHT if (allowed_dirs & VCGDirectionCode::RIGHT) {return VCGDirectionCode::RIGHT;}
#define TRY_UP if (allowed_dirs & VCGDirectionCode::UP) {return VCGDirectionCode::UP;}
#define TRY_DOWN if (allowed_dirs & VCGDirectionCode::DOWN) {return VCGDirectionCode::DOWN;}


uint8_t compute_next_move(
	const bool clockwise,
	const uint8_t last_move,
	const uint8_t allowed_dirs
) {
	if (clockwise) {
		if (last_move == VCGDirectionCode::RIGHT) {
			TRY_DOWN
			else TRY_RIGHT
			else TRY_UP
			else TRY_LEFT
		}
		else if (last_move == VCGDirectionCode::LEFT) {
			TRY_UP
			else TRY_LEFT
			else TRY_DOWN
			else TRY_RIGHT
		}
		else if (last_move == VCGDirectionCode::UP) {
			TRY_RIGHT
			else TRY_UP
			else TRY_LEFT
			else TRY_DOWN
		}
		else { // last_move == 'd'
			TRY_LEFT
			else TRY_DOWN
			else TRY_RIGHT
			else TRY_UP
		}
	}
	else {
		if (last_move == VCGDirectionCode::RIGHT) {
			TRY_UP
			else TRY_RIGHT
			else TRY_DOWN
			else TRY_LEFT
		}
		if (last_move == VCGDirectionCode::LEFT) {
			TRY_DOWN
			else TRY_LEFT
			else TRY_UP
			else TRY_RIGHT
		}
		if (last_move == VCGDirectionCode::UP) {
			TRY_LEFT
			else TRY_UP
			else TRY_RIGHT
			else TRY_DOWN
		}
		else { // last_move == 'd'
			TRY_RIGHT
			else TRY_DOWN
			else TRY_LEFT
			else TRY_UP
		}
	}

	return VCGDirectionCode::NONE;
}


#undef TRY_UP
#undef TRY_DOWN
#undef TRY_LEFT
#undef TRY_RIGHT

void extract_contours_helper(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy,
	std::vector<std::vector<uint32_t>>& contours,
	std::vector<std::vector<uint32_t>>& hole_contours
) {

	VCGGraph G(vcg, sx, sy);
	for (uint64_t i = 0; i < sx; i++) {
		vcg[i] = vcg[i] & ~VCGDirectionCode::UP;
		int idx = i + sx * (sy-1);
		vcg[idx] = vcg[idx] & ~VCGDirectionCode::DOWN;
	}
	for (uint64_t i = 0; i < sy; i++) {
		int idx = sx * i;
		vcg[idx] = vcg[idx] & ~VCGDirectionCode::LEFT;
		idx = (sx-1) + sx * i;
		vcg[idx] = vcg[idx] & ~VCGDirectionCode::RIGHT;
	}

	// clockwise for outer boundaries
	// counterclockwise for inner boundaries
	// bool is_hole = false; 
	bool clockwise = true;
	int64_t start_node = 0;
	uint32_t barriers = 0;

	// corresponds to VCGDirectionCodes
	int64_t move_amt[9];
	move_amt[VCGDirectionCode::NONE] = 0;
	move_amt[VCGDirectionCode::RIGHT] = 1;
	move_amt[VCGDirectionCode::LEFT] = -1;
	move_amt[VCGDirectionCode::DOWN] = static_cast<int64_t>(sx);
	move_amt[VCGDirectionCode::UP] = -static_cast<int64_t>(sx);

	// Moore Neighbor Tracing variation
	int64_t y = 0; // breaking abstraction to save a frequent division
	while (G.next_contour(barriers, start_node, y)) {

		// is_hole = (barriers & 0b1) == 1;
		std::vector<uint32_t> connected_component;

		int64_t node = start_node;
		uint8_t allowed_dirs = vcg[node] & 0b1111;
		uint8_t next_move, ending_orientation;

		if (allowed_dirs == VCGDirectionCode::NONE) {
			vcg[node] |= VISITED_BIT;
			connected_component.push_back(node);
		}
		else {
			connected_component.reserve(100);
			connected_component.push_back(start_node);

			// go counterclockwise for |x  vs clockwise for x|
			next_move = VCGDirectionCode::UP;
			clockwise = (vcg[start_node] & 0b1) == 0;

			ending_orientation = compute_next_move(
				clockwise, next_move, allowed_dirs
			);
			next_move = ending_orientation;

			do {
				node += move_amt[next_move];
				connected_component.push_back(node);
				uint8_t visit_count = vcg[node] >> 4;
				vcg[node] = ((visit_count + 1) << 4) | (vcg[node] & 0b1111);
				allowed_dirs = vcg[node] & 0b1111;
				next_move = compute_next_move(
					clockwise, next_move, allowed_dirs
				);
			} while (
				!(node == start_node && next_move == ending_orientation)
			);
		}

		if (connected_component.size() == 0) {
			continue;
		}

		// Rotate the contour so that the characteristic min element
		// is at the front.
		std::vector<uint32_t>::iterator it = std::min_element(connected_component.begin(), connected_component.end());
		std::vector<uint32_t> rotated;
		rotated.reserve(connected_component.size());

		rotated.insert(rotated.end(), it, connected_component.end());
		rotated.insert(rotated.end(), connected_component.begin(), it);

		// if (is_hole) {
		// 	hole_contours.push_back(std::move(rotated));
		// }
		// else {
		contours.push_back(std::move(rotated));
		// }
	}
}

bool polygonContainsPoint(
	const std::vector<uint32_t>& poly, 
	const uint32_t idx, 
	const uint64_t sx
) {

	uint32_t pt_y = idx / sx;
	uint32_t contacts = 0;

	for (uint64_t i = 0; i < poly.size(); i++) {
		uint32_t elem_y = poly[i] / sx;
		contacts += (elem_y > pt_y);	
	}

	return contacts & 0b1;
}


std::vector<std::vector<uint32_t>>
merge_holes(
	std::vector<std::vector<uint32_t>>& candidate_contours,
	const uint64_t sx
) {

	std::vector<std::pair<uint32_t, uint32_t>> bboxes(candidate_contours.size());
	std::vector<uint16_t> merge_ct(candidate_contours.size());
	crackle::cc3d::DisjointSet<uint32_t> equivalences(candidate_contours.size() + 1);

	for (uint64_t i = 0; i < candidate_contours.size(); i++) {
		auto& vec = candidate_contours[i];
		std::vector<uint32_t>::iterator it = std::max_element(vec.begin(), vec.end());
		bboxes[i].first = vec[0];
		bboxes[i].second = *it;
		equivalences.add(i+1); // +1 bc 0 is "null" in DisjointSet
	}

	// printf("3\n");

	for (uint64_t i = 0; i < candidate_contours.size(); i++) {
		for (uint64_t j = i + 1; j < candidate_contours.size(); j++) {
			auto& bbx1 = bboxes[i];
			auto& bbx2 = bboxes[j];

			// non-intersecting bounding boxes
			// Without parsing these into x,y coords first
			// this test mainly resticts the y axis, so performance can be
			// improved a lot.
			if (
				!(bbx2.first >= bbx1.first && bbx2.second <= bbx1.second)
			) {
				continue;
			}

			if (polygonContainsPoint(candidate_contours[i], candidate_contours[j][0], sx)) {
				merge_ct[j]++;
				// if an odd number of merges, then it's a "hole"
				// otherwise it's its own connected component.
				if ((merge_ct[j] & 0b1) == 0) {
					equivalences.ids[j] = j+1;
				}
				else {
					equivalences.unify(i+1, j+1); // +1 bc 0 is null in DisjointSet
				}
			}
		}
	}

	std::vector<std::vector<uint32_t>> merged_contours(candidate_contours.size());

	for (uint64_t i = 0; i < merged_contours.size(); i++) {
		uint32_t idx = equivalences.root(i + 1) - 1;

		merged_contours[idx].insert(
			merged_contours[idx].end(), 
			candidate_contours[i].begin(), 
			candidate_contours[i].end()
		);
	}

  	merged_contours.erase(
        std::remove_if(merged_contours.begin(), merged_contours.end(), [](const std::vector<uint32_t>& item) { return item.empty(); }),
        merged_contours.end()
    );


	std::sort(merged_contours.begin(), merged_contours.end(),
		[](const auto& a, const auto& b) {
			return a[0] < b[0];
		});

	return merged_contours;
}


std::vector<std::vector<uint32_t>> 
extract_contours(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy 
) {

	std::vector<std::vector<uint32_t>> contours;
	std::vector<std::vector<uint32_t>> hole_contours;

	extract_contours_helper(vcg, sx, sy, contours, hole_contours);

	// return merge_holes(contours, sx);



	// std::sort(contours.begin(), contours.end(),
	// 	[](const auto& a, const auto& b) {
	// 		return a[0] < b[0];
	// 	});

	// std::sort(hole_contours.begin(), hole_contours.end(),
	// 	[](const auto& a, const auto& b) {
	// 		return a[0] < b[0];
	// 	});

	// // merge holes with parent contours
	// for (uint64_t j = 0; j < hole_contours.size(); j++) {
	// 	auto& hc = hole_contours[j];
	// 	auto min_element_hole = hc[0];

	// 	for (uint64_t i = 0; i < contours.size(); i++) {
	// 		auto& contour = contours[i];
	// 		auto min_element_ct = contour[0];
			
	// 		if (min_element_ct > min_element_hole) {
	// 			break;
	// 		}
			
	// 		contour.insert(contour.end(), hc.begin(), hc.end());
	// 		break;
	// 	}
	// }

	return contours;
}

};
};

#endif