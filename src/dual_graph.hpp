#ifndef __CRACKLE_DUALGRAPH_HXX__
#define __CRACKLE_DUALGRAPH_HXX__

#include <vector>
#include <cstdint>
#include <algorithm>

#include "builtins.hpp"
#include "cc3d.hpp"
#include "libdivide.hpp"

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

	bool next_contour(int64_t& idx, int64_t& y) {
		int64_t x = idx - sx * y;

		for (; y < sy; y++) {
			for (; x < sx; x++, idx++) {
				// condensing this conditional seems to save 5% in one speed test
				// if (((vcg[idx] & 0b11) < 0b11) && (vcg[idx] & VISIT_COUNT) == 0) {
				// -----
				// check that the next voxel isn't visited and is a barrier

				if ((vcg[idx] & 0b110011) < 0b11 
					|| (x < sx - 1 && (vcg[idx+1] & 0b11110010) == 0b0)) {
					return true;
				}
			}
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

std::vector<std::vector<uint32_t>> 
extract_contours_helper(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy
) {
	std::vector<std::vector<uint32_t>> contours;

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
	bool clockwise = true;
	int64_t start_node = 0;

	// corresponds to VCGDirectionCodes
	int64_t move_amt[9];
	move_amt[VCGDirectionCode::NONE] = 0;
	move_amt[VCGDirectionCode::RIGHT] = 1;
	move_amt[VCGDirectionCode::LEFT] = -1;
	move_amt[VCGDirectionCode::DOWN] = static_cast<int64_t>(sx);
	move_amt[VCGDirectionCode::UP] = -static_cast<int64_t>(sx);

	// Moore Neighbor Tracing variation
	int64_t y = 0; // breaking abstraction to save a frequent division
	while (G.next_contour(start_node, y)) {

		std::vector<uint32_t> connected_component;

		int64_t node = start_node;
		uint8_t allowed_dirs = vcg[node] & 0b1111;
		uint8_t next_move, ending_orientation;

		uint64_t nodes_already_visited = (vcg[node] >> 4) > 0;

		if (allowed_dirs == VCGDirectionCode::NONE) {
			vcg[node] |= VISITED_BIT;
			connected_component.push_back(node);
		}
		else {
			connected_component.reserve(100);
			connected_component.push_back(start_node);

			// go counterclockwise for |x  vs clockwise for x|
			next_move = VCGDirectionCode::UP;
			clockwise = ((vcg[start_node] & 0b1) == 0) || (vcg[start_node] == 0b11100);

			ending_orientation = compute_next_move(
				clockwise, next_move, allowed_dirs
			);
			next_move = ending_orientation;

			do {
				node += move_amt[next_move];
				connected_component.push_back(node);
				uint8_t is_visited = (vcg[node] >> 4) > 0;
				nodes_already_visited += is_visited;
				vcg[node] |= VISITED_BIT;
				allowed_dirs = vcg[node] & 0b1111;

				next_move = compute_next_move(
					clockwise, next_move, allowed_dirs
				);
			} while (
				!(node == start_node && next_move == ending_orientation)
			);
		}

		start_node++;

		if (connected_component.size() == 0 || connected_component.size() == nodes_already_visited) {
			continue;
		}

		// Rotate the contour so that the characteristic min element
		// is at the front.
		std::vector<uint32_t>::iterator it = std::min_element(connected_component.begin(), connected_component.end());
		std::vector<uint32_t> rotated;
		rotated.reserve(connected_component.size());

		rotated.insert(rotated.end(), it, connected_component.end());
		rotated.insert(rotated.end(), connected_component.begin(), it);

		contours.push_back(std::move(rotated));
	}

	return contours;
}

// an alternative is the complex computational geometry strategy above...
// which is much lower memory but hairy to implement and much slower so
// long as we have that dang quadratic loop.
std::vector<std::vector<uint32_t>> merge_contours_via_vcg_coloring(
	const std::vector<std::vector<uint32_t>>& contours,
	const std::vector<uint8_t>& vcg,
	std::unique_ptr<uint32_t[]>& cc_labels,
	const uint64_t N,
	const uint64_t sx, const uint64_t sy
) {
	std::vector<std::vector<uint32_t>> merged_contours(N);
	for (uint64_t i = 0; i < contours.size(); i++) {
		auto& contour = contours[i];
		uint32_t cc_label = cc_labels[contour[0]];

		auto insertion_point_it = merged_contours[cc_label].end();
		if (
			merged_contours[cc_label].size() > 0 
			&& merged_contours[cc_label][0] > contour[0]
		) {
			insertion_point_it = merged_contours[cc_label].begin();
		}

		merged_contours[cc_label].insert(
			insertion_point_it, 
			contour.begin(), 
			contour.end()
		);
	}

	return merged_contours;
}


std::vector<std::vector<uint32_t>> 
extract_contours(
	std::vector<uint8_t>& vcg,
	std::unique_ptr<uint32_t[]>& cc_labels,
	const uint64_t N,
	const uint64_t sx, const uint64_t sy 
) {
	std::vector<std::vector<uint32_t>> contours = extract_contours_helper(vcg, sx, sy);
	return merge_contours_via_vcg_coloring(contours, vcg, cc_labels, N, sx, sy);
}

};
};

#endif