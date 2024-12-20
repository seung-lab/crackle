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

struct VCGGraph {
	std::vector<uint8_t>& vcg;
	
	int64_t sx;
	int64_t sy;

	VCGGraph(std::vector<uint8_t>& _vcg, int64_t _sx, int64_t _sy) 
		: vcg(_vcg), sx(_sx), sy(_sy) {}

	int64_t next_contour(int64_t idx) {
		for (int64_t i = idx; i < sx * sy; i++) {
			if (vcg[i] != 0b1111 && (vcg[i] & 0b10000) == 0) {
				return i;
			}
		}

		return -1;
	}

};

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
	0b0001  // 0b1111
};

enum VCGDirectionCode {
	LEFT = 0b0010,
	RIGHT = 0b0001,
	UP = 0b1000,
	DOWN = 0b0100,
	ANY = 0b1111,
	VISITED = 0b10000
};

void print_bits(uint8_t val) {
	printf("vcg: 0b");
	for (int i = 4; i >= 0; i--) {
		printf("%d", (val >> i) & 0b1);
	}
}

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

	int64_t start_node = 0;

	std::vector<std::vector<uint32_t>> connected_components;

	const uint8_t visited_bit = 0b10000;

	while ((start_node = G.next_contour(start_node)) != -1) {

		std::vector<uint32_t> connected_component;
		uint32_t node = start_node;
		
		while ((vcg[node] & visited_bit) == 0) {
			// int y = node / sx;
			// int x = node - sx * y;
			// printf("x %d y %d ", x, y);
			// print_bits(vcg[node]);
			// printf("\n");

			connected_component.push_back(node);
			uint8_t allowed_dirs = contour_lookup[vcg[node]];
			vcg[node] = allowed_dirs | visited_bit;
			// print_bits(vcg[node]);
			// printf("\n");

			// circulate clockwise
			if (allowed_dirs & VCGDirectionCode::RIGHT) {
				node += 1;
				// printf("r\n");
			}
			else if (allowed_dirs & VCGDirectionCode::DOWN) {
				node += sx;
				// printf("d\n");
			}
			else if (allowed_dirs & VCGDirectionCode::LEFT) {
				node -= 1;
				// printf("l\n");
			}
			else if (allowed_dirs & VCGDirectionCode::UP) {
				node -= sx;
				// printf("u\n");
			}
		}

		connected_components.push_back(connected_component);
	}

    std::sort(connected_components.begin(), connected_components.end(), 
    	[](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
	        uint32_t minA = *std::min_element(a.begin(), a.end());
	        uint32_t minB = *std::min_element(b.begin(), b.end());
	        return minA < minB;
	    });

	return connected_components;
}

};
};

#endif