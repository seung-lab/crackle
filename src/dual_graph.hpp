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

// this array specifies the allowed
// directions when you hit a wall.
// ->| means you can go up or down
// to follow the wall contour

// -y +y -x +x
const uint8_t contour_lookup[16] = {
	0b0000, // 0b0000
	0b0001, // 0b0001
	0b0010, // 0b0010
	0b0011, // 0b0011
	0b0100, // 0b0100
	0b0101, // 0b0101
	0b0110, // 0b0110
	0b0011, // 0b0111
	0b1000, // 0b1000
	0b1001, // 0b1001
	0b1010, // 0b1010
	0b0011, // 0b1011
	0b1100, // 0b1100
	0b1100, // 0b1101
	0b1100, // 0b1110
	0b1111  // 0b1111
};

const uint8_t VISITED_BIT = 0b10000;

enum VCGDirectionCode {
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

	for (int y = 0; y < sy; y++) {
		for (int x = 0; x < sx; x++) {
			print_bits(vcg[x + sx * y]);
			printf(" ");
		}
		printf("\n");
	}

	std::vector<std::vector<uint32_t>> contours;
	std::vector<std::vector<uint32_t>> hole_contours;

	// clockwise for outer boundaries
	// counterclockwise for inner boundaries
	bool is_hole = false; 
	bool clockwise = true;
	int64_t start_node = 0;
	uint32_t barriers = 0;


	// Moore Neighbor Tracing variation
	while (G.next_contour(barriers, start_node)) {

		is_hole = (barriers & 0b1) == 0;
		clockwise = is_hole;

		// printf("clockwise %d start_node %d\n", clockwise, start_node);
		std::vector<uint32_t> connected_component;
		uint32_t node = start_node;

		char last_move = (clockwise ? 'u' : 'd'); 

		// printf("barriers: %d\n", barriers);


#define TRY_LEFT if (allowed_dirs & VCGDirectionCode::LEFT) {\
						node -= 1;\
						last_move = 'l';\
					}

#define TRY_RIGHT if (allowed_dirs & VCGDirectionCode::RIGHT) {\
						node += 1;\
						last_move = 'r';\
					}

#define TRY_UP if (allowed_dirs & VCGDirectionCode::UP) {\
						node -= sx;\
						last_move = 'u';\
					}

# define TRY_DOWN if (allowed_dirs & VCGDirectionCode::DOWN) {\
						node += sx;\
						last_move = 'd';\
					}

		do {
			// int y = node / sx;
			// int x = node - sx * y;
			// printf("x %d y %d last: %c ", x, y, last_move);
			// print_bits(vcg[node]);
			// printf(" node %d hole %d clockwise %d \n", node, is_hole, clockwise);

			connected_component.push_back(node);
			uint8_t allowed_dirs = contour_lookup[vcg[node]];
			vcg[node] = allowed_dirs | VISITED_BIT;

			if (clockwise) {
				if (last_move == 'r') {
					TRY_DOWN
					else TRY_RIGHT
					else TRY_UP
					else TRY_LEFT
				}
				else if (last_move == 'l') {
					TRY_UP
					else TRY_LEFT
					else TRY_DOWN
					else TRY_RIGHT
				}
				else if (last_move == 'u') {
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
				if (last_move == 'r') {
					TRY_UP
					else TRY_RIGHT
					else TRY_DOWN
					else TRY_LEFT
				}
				else if (last_move == 'l') {
					TRY_DOWN
					else TRY_LEFT
					else TRY_UP
					else TRY_RIGHT
				}
				else if (last_move == 'u') {
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
		} while (node != start_node); // defective version of moore algorithm

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

#undef TRY_UP
#undef TRY_DOWN
#undef TRY_LEFT
#undef TRY_RIGHT

#endif