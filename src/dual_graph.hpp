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
	bool next_contour(bool& clockwise, int64_t& idx) {
		int64_t y = idx / sx;
		int64_t x = idx - sx * y;
		// printf("<%d %d %d %d>\n", x, y, sx, sy);
		for (; y < sy; y++) {
			int barriers = 0;
			for (; x < sx; x++, idx++) {
				barriers += static_cast<int>((vcg[idx] & 0b11) < 0b11);
				if (((vcg[idx] & 0b11) < 0b11) && (vcg[idx] & VISITED_BIT) == 0) {
					// odd clockwise, even counterclockwise
					// this is bc the outer border has a barrier already at 0
					clockwise = true;//(barriers & 0b1) == 1; 
					// printf("barriers: %d\n", barriers);
					return true;
				}
			}
			x = 0;
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

	// for (int y = 0; y < sy; y++) {
	// 	for (int x = 0; x < sx; x++) {
	// 		print_bits(vcg[x + sx * y]);
	// 		printf(" ");
	// 	}
	// 	printf("\n");
	// }

	std::vector<std::vector<uint32_t>> connected_components;

	// clockwise for outer boundaries
	// counterclockwise for inner boundaries
	bool clockwise = true; 
	int64_t start_node = 0;

	// Moore Neighbor Tracing variation

	while (G.next_contour(clockwise, start_node)) {

		// printf("clockwise %d start_node %d\n", clockwise, start_node);
		std::vector<uint32_t> connected_component;
		uint32_t node = start_node;
	
		char last_move = 'r'; 

		while ((vcg[node] & VISITED_BIT) == 0) {
			// int y = node / sx;
			// int x = node - sx * y;
			// printf("x %d y %d last: %c ", x, y, last_move);
			// print_bits(vcg[node]);
			// printf("\n");

			connected_component.push_back(node);
			uint8_t allowed_dirs = contour_lookup[vcg[node]];
			vcg[node] = allowed_dirs | VISITED_BIT;

			if (clockwise) {
				if (last_move == 'r') {
					if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
					else if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
					else if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
				}
				else if (last_move == 'l') {
					if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
					else if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
					else if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
				}
				else if (last_move == 'u') {
					if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
					else if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
					else if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
				}
				else { // last_move == 'd'
					if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
					else if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
					else if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
				}
			}
			else {
				if (last_move == 'r') {
					if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
					else if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
					else if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
				}
				else if (last_move == 'l') {
					if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
					else if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
					else if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
				}
				else if (last_move == 'u') {
					if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
					
					else if (allowed_dirs & VCGDirectionCode::UP) {
						node -= sx;
						last_move = 'u';
					}
					else if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
				}
				else { // last_move == 'd'
					if (allowed_dirs & VCGDirectionCode::LEFT) {
						node -= 1;
						last_move = 'l';
					}
					else if (allowed_dirs & VCGDirectionCode::DOWN) {
						node += sx;
						last_move = 'd';
					}
					else if (allowed_dirs & VCGDirectionCode::RIGHT) {
						node += 1;
						last_move = 'r';
					}
				}
			}
		}

		// printf("NEW COMPONENT\n");
		connected_components.push_back(std::move(connected_component));
	}

	// printf("final clockwise %d start_node %d\n", clockwise, start_node);

	std::vector<std::pair<uint32_t, const std::vector<uint32_t>*>> min_with_vectors;
	for (const auto& vec : connected_components) {
	    if (!vec.empty()) {
	        min_with_vectors.emplace_back(*std::min_element(vec.begin(), vec.end()), &vec);
	    }
	}

	std::sort(min_with_vectors.begin(), min_with_vectors.end(),
	          [](const auto& a, const auto& b) {
	              return a.first < b.first;
	          });

	std::vector<std::vector<uint32_t>> sorted_components;
	for (const auto& pair : min_with_vectors) {
	    sorted_components.push_back(*pair.second);
	}
	connected_components = std::move(sorted_components);

	return connected_components;
}

};
};

#endif