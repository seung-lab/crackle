#ifndef __CRACKLE_DUALGRAPH_HXX__
#define __CRACKLE_DUALGRAPH_HXX__

#include <vector>
#include <cstdint>
#include <algorithm>

#include "robin_hood.hpp"
#include "crackcodes.hpp"
#include "builtins.hpp"

namespace crackle {
namespace dual_graph {

const uint8_t VISITED_BIT = 0b10000;

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

	// returns clockwise and next id
	bool next_contour(uint32_t& barriers, int64_t& idx) {
		int64_t y = idx / sx;
	
		idx = sx * y;

		for (; y < sy; y++) {
			barriers = 0;
			for (int64_t x = 0; x < sx; x++, idx++) {
				barriers += static_cast<uint32_t>((vcg[idx] & 0b11) < 0b11);
				if (((vcg[idx] & 0b11) < 0b11) && (vcg[idx] & VISITED_BIT) == 0) {
					return true;
				}
			}
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
extract_contours(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy 
) {

	// for (int y = 0; y < sy; y++) {
	// 	for (int x = 0; x < sx; x++) {
	// 		print_bits(vcg[x + sx * y]);
	// 		printf(" ");
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

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

	// for (int y = 0; y < sy; y++) {
	// 	for (int x = 0; x < sx; x++) {
	// 		print_bits(vcg[x + sx * y]);
	// 		printf(" ");
	// 	}
	// 	printf("\n");
	// }

	std::vector<std::vector<uint32_t>> contours;
	std::vector<std::vector<uint32_t>> hole_contours;

	// clockwise for outer boundaries
	// counterclockwise for inner boundaries
	bool is_hole = false; 
	bool clockwise = true;
	int64_t start_node = 0;
	uint32_t barriers = 0;

	// corresponds to VCGDirectionCodes
	int64_t move_amt[9] = { 
		0, // none
		1, // right
		-1, // left
		0, // n/a
		static_cast<int64_t>(sx), // down
		0, // n/a 
		0, // n/a
		0, // n/a 
		-static_cast<int64_t>(sx) // up
	};

	char translate[9];
	translate[VCGDirectionCode::NONE] = 'x';
	translate[VCGDirectionCode::LEFT] = 'l';
	translate[VCGDirectionCode::RIGHT] = 'r';
	translate[VCGDirectionCode::UP] = 'u';
	translate[VCGDirectionCode::DOWN] = 'd';

	// Moore Neighbor Tracing variation
	while (G.next_contour(barriers, start_node)) {

		is_hole = (barriers & 0b1) == 0;
		clockwise = is_hole;

		// printf("clockwise %d start_node %d\n", clockwise, start_node);
		std::vector<uint32_t> connected_component;
		int64_t node = start_node;
		uint8_t allowed_dirs = vcg[node] & 0b1111;
		uint8_t next_move, ending_orientation;

		// int y = node / sx;
		// int x = node - sx * y;

		if (allowed_dirs == VCGDirectionCode::NONE) {
			vcg[node] |= VISITED_BIT;
			connected_component.push_back(node);
		}
		else {
			next_move = (clockwise ? VCGDirectionCode::UP : VCGDirectionCode::DOWN); 
			ending_orientation = compute_next_move(clockwise, next_move, allowed_dirs);
			next_move = ending_orientation;

			// printf("x %d y %d next: %c ", x, y, translate[next_move]);
			// print_bits(vcg[node]);
			// printf(" node %d hole %d clockwise %d start %d\n", node, is_hole, clockwise, start_node);


			// printf("barriers: %d\n", barriers);

			do {
				node += move_amt[next_move];
				// y = node / sx;
				// x = node - sx * y;

				connected_component.push_back(node);
				vcg[node] |= VISITED_BIT;
				next_move = compute_next_move(clockwise, next_move, (vcg[node] & 0b1111));

				// printf("x %d y %d next: %c ", x, y, translate[next_move]);
				// print_bits(vcg[node]);
				// printf(" node %d hole %d clockwise %d start %d\n", node, is_hole, clockwise, start_node);
			} while (!(node == start_node && (vcg[node] & VISITED_BIT) && next_move == ending_orientation));
		}
		// printf("NEW COMPONENT\n");

		if (connected_component.size() == 0) {
			continue;
		}

		// make it slightly easier to parse contours after merging
		// as each contour starts and ends with the same node
		connected_component.push_back(start_node); 

		// Rotate the contour so that the characteristic min element
		// is at the front.
		std::vector<uint32_t>::iterator it = std::min_element(connected_component.begin(), connected_component.end());
		std::vector<uint32_t> rotated;
		rotated.reserve(connected_component.size());

		rotated.insert(rotated.end(), it, connected_component.end());
		rotated.insert(rotated.end(), connected_component.begin(), it);

		if (is_hole) {
			hole_contours.push_back(std::move(rotated));
		}
		else {
			contours.push_back(std::move(rotated));
		}
	}

	std::sort(contours.begin(), contours.end(),
		[](const auto& a, const auto& b) {
			return a[0] < b[0];
		});

	std::sort(hole_contours.begin(), hole_contours.end(),
		[](const auto& a, const auto& b) {
			return a[0] < b[0];
		});

	// printf("MERGING\n");

	// merge holes with parent contours
	for (uint64_t j = 0; j < hole_contours.size(); j++) {
		auto& hc = hole_contours[j];
		auto min_element_hole = hc[0];

		for (uint64_t i = 0; i < contours.size(); i++) {
			auto& contour = contours[i];
			auto min_element_ct = contour[0];
			
			if (min_element_ct > min_element_hole) {
				break;
			}
			// printf("Merged c%d hc%d %d %d\n", i, j, min_element_ct, min_element_hole);
			contour.insert(contour.end(), hc.begin(), hc.end());
			break;
		}
	}

	return contours;
}

};
};

#endif