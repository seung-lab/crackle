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
			int y = node / sx;
			int x = node - sx * y;
			printf("x %d y %d ", x, y);
			print_bits(vcg[node]);
			printf("\n");

			connected_component.push_back(node);
			uint8_t allowed_dirs = contour_lookup[vcg[node]];
			vcg[node] = allowed_dirs | visited_bit;
			print_bits(vcg[node]);
			printf("\n");

			// circulate clockwise
			if (allowed_dirs & VCGDirectionCode::RIGHT) {
				node += 1;
				// vcg[node] = vcg[node] & ~VCGDirectionCode::RIGHT;
				printf("r\n");
			}
			else if (allowed_dirs & VCGDirectionCode::DOWN) {
				node += sx;
				// vcg[node] = vcg[node] & ~VCGDirectionCode::DOWN;
				printf("d\n");
			}
			else if (allowed_dirs & VCGDirectionCode::LEFT) {
				node -= 1;
				// vcg[node] = vcg[node] & ~VCGDirectionCode::LEFT;
				printf("l\n");
			}
			else if (allowed_dirs & VCGDirectionCode::UP) {
				node -= sx;
				// vcg[node] = vcg[node] & ~VCGDirectionCode::UP;
				printf("u\n");
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



// template <typename T>
// class MapDisjointSet {
// public:
//   robin_hood::unordered_flat_map<T,T> ids;

//   MapDisjointSet () {}

//   MapDisjointSet (const MapDisjointSet &cpy) {
//   	ids.insert(ids.end(), cpy.ids.begin(), cpy.ids.end());
//   }

//   T root (T n) {
//     T i = ids[n];
//     while (i != ids[i]) {
//       ids[i] = ids[ids[i]]; // path compression
//       i = ids[i];
//     }

//     return i;
//   }

//   bool find (T p, T q) {
//     return root(p) == root(q);
//   }

//   void add(T p) {
//     if (ids[p] == 0) {
//       ids[p] = p;
//     }
//   }

//   void unify (T p, T q) {
//     if (p == q || q < 0) {
//       return;
//     }

//     T i = root(p);
//     T j = root(q);

//     if (i == 0) {
//       add(p);
//       i = p;
//     }

//     if (j == 0) {
//       add(q);
//       j = q;
//     }

//     ids[i] = j;
//   }

//   T renumber() {
//   	std::vector<T> uniq;
//   	for (auto pair : ids) {
//   		uniq.push_back(root(pair.second));
//   	}

//   	std::sort(uniq.begin(), uniq.end());

//   	robin_hood::unordered_flat_map<T, T> relabel;

//   	T N = 1;
//   	relabel[uniq[0]] = N;

//   	for (uint64_t i = 1; i < uniq.size(); i++) {
//   		if (uniq[i] != uniq[i-1]) {
//   			N++;
//   			relabel[uniq[i]] = N;
//   		}
//   	}

//   	for (auto pair : ids) {
//   		pair.second = relabel[pair.second];
//   	}

//   	return N;
//   }
// };

// // returns lists of contours sorted by connected component
// std::vector<std::vector<uint32_t>> 
// crack_codes_to_dual_graph(
// 	const std::vector<std::pair<uint64_t, std::vector<unsigned char>> >& chains,
// 	const uint64_t sx, const uint64_t sy
// ) {

// 	int64_t sxe = sx + 1;

// 	const int64_t pixels = (sx+1) * (sy+1);

// 	MapDisjointSet<int32_t> equivalences;

// 	// graph is of corners and edges
// 	// origin is located at top left
// 	// corner of the image
// 	for (auto& [node, symbols]: chains) {
// 		int64_t y = node / sxe;
// 		int64_t x = node - (sxe * y);
// 		int64_t loc = x + sx * y;

// 		unsigned char last_symbol = symbols[0];

// 		std::stack<int64_t> revisit;
// 		for (unsigned char symbol : symbols) {
// 			if (loc >= pixels) {
// 				throw std::runtime_error("crackle: crack_codes_to_dual_graph: index out of range.");
// 			}

// 			if (symbol == 'u') {
// 				if (last_symbol == 'u') {
// 					equivalences.unify(loc-sx, loc);
// 					equivalences.unify(loc-1-sx, loc-1);
// 				}
// 				else if (last_symbol == 'l') {
// 					equivalences.unify(loc-1, loc);
// 					equivalences.unify(loc-1-sx, loc - 1);
// 				}
// 				else {
// 					equivalences.unify(loc, loc - 1);
// 					equivalences.unify(loc-sx, loc);
// 				}

// 				y--;
// 				loc -= sx;
// 			}
// 			else if (symbol == 'd') {
// 				if (last_symbol == 'd') {
// 					equivalences.unify(loc, loc - sx);
// 					equivalences.unify(loc-1, loc - 1 - sx);
// 				}
// 				else if (last_symbol == 'l') {
// 					equivalences.unify(loc-1, loc - 1 - sx);
// 					equivalences.unify(loc-1-sx, loc - sx);
// 				}
// 				else {
// 					equivalences.unify(loc-sx, loc - sx - 1);
// 					equivalences.unify(loc, loc - sx);
// 				}
// 				y++;
// 				loc += sx;
// 			}
// 			else if (symbol == 'l') {
// 				if (last_symbol == 'l') {
// 					equivalences.unify(loc-1, loc);
// 					equivalences.unify(loc-1-sx, loc - sx);
// 				}
// 				else if (last_symbol == 'd') {
// 					equivalences.unify(loc-1, loc);
// 					equivalences.unify(loc, loc-sx);
// 				}
// 				else {
// 					equivalences.unify(loc-sx, loc);
// 					equivalences.unify(loc-1-sx, loc - sx);
// 				}
// 				x--;
// 				loc--;
// 			}
// 			else if (symbol == 'r') {
// 				if (last_symbol == 'r') {
// 					equivalences.unify(loc, loc - 1);
// 					equivalences.unify(loc-sx, loc - sx - 1);
// 				}
// 				else if (last_symbol == 'd') {
// 					equivalences.unify(loc, loc - 1);
// 					equivalences.unify(loc-1, loc - sx - 1);
// 				}
// 				else {
// 					equivalences.unify(loc-sx, loc - sx - 1);
// 					equivalences.unify(loc-sx-1, loc - 1);
// 				}
// 				x++;
// 				loc++;
// 			}
// 			else if (symbol == 'b') {
// 				revisit.push(loc);
// 			}
// 			else if (symbol =='t') {
// 				if (!revisit.empty()) {
// 					loc = revisit.top();
// 					revisit.pop();
// 					y = loc / sx;
// 					x = loc - (sx * y);
// 				}
// 			}

// 			last_symbol = symbol;
// 		}
// 	}

// 	uint32_t N = equivalences.renumber();
// 	printf("%d\n", N);

// 	std::vector<std::vector<uint32_t>> connected_components(N-1);

// 	for (const auto& pair : equivalences.ids) {
// 		if (pair.first == 0) {
// 			continue;
// 		}
// 		connected_components[equivalences.root(pair.second) - 1].push_back(pair.first);
// 	}

//     std::sort(connected_components.begin(), connected_components.end(), 
//     	[](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
// 	        uint32_t minA = *std::min_element(a.begin(), a.end());
// 	        uint32_t minB = *std::min_element(b.begin(), b.end());
// 	        return minA < minB;
// 	    });

// 	return connected_components;
// }

};
};

#endif