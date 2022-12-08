#ifndef __CRACKLE_CRACKCODE_HXX__
#define __CRACKLE_CRACKCODE_HXX__

#include <vector>
#include <unordered_map>
#include <stack>
#include <cstdint>

#include "lib.hpp"

namespace crackle {
namespace crackcodes {

template <typename T>
std::unordered_map<uint64_t, std::vector<unsigned char>> 
unpack_binary_helper(
	const std::vector<unsigned char> &code, 
	const uint64_t sx, const uint64_t sy
) {
	std::unordered_map<uint64_t, std::vector<unsigned char>> chains;
	const std::vector<T> int_code(code.begin(), code.end());

	std::vector<unsigned char> symbols;
	symbols.reserve(code.size() * 4);

	uint64_t branches_taken = 0;
	T node = 0;

	char remap[4] = { 'u', 'r', 'l', 'd' };

	uint64_t num_moves = sizeof(T) * 8 / 2;

	for (T moveset : int_code) {
		if (branches_taken == 0) {
			node = moveset;
			branches_taken = 1;
			continue;
		}

		for (uint64_t i = 0; i < num_moves; i++) {
			uint8_t move = static_cast<uint8_t>((moveset >> (2*i)) & 0b11);

			if (symbols.size()) {
				if (
					(move == 0 && symbols.back() == 'd')
					|| (move == 2 && symbols.back() == 'r')
				) {
					symbols.back() = 't';
				}
				else if (
					(move == 3 && symbols.back() == 'u')
					|| (move == 1 && symbols.back() == 'l')
				) {
					symbols.back() = 'b';
				}
				else {
					symbols.push_back(remap[move]);
				}
			}
			else {
				symbols.push_back(remap[move]);
			}
		}

		if (branches_taken == 0) {
			auto vec = chains[node];
			vec.insert(vec.end(), symbols.begin(), symbols.end());
		}
	}

	return chains;
}

// decodes to symbols
std::unordered_map<uint64_t, std::vector<unsigned char>> unpack_binary(
	const std::vector<unsigned char> &code, 
	const uint64_t sx, const uint64_t sy
) {
	std::unordered_map<uint64_t, std::vector<unsigned char>> chains;

	if (code.size() == 0) {
		return chains;
	}

	uint64_t byte_width = crackle::lib::compute_byte_width((sx+1) * (sy+1));
	
	if (byte_width == 1) {
		return unpack_binary_helper<uint8_t>(code, sx, sy);
	}
	else if (byte_width == 2) {
		return unpack_binary_helper<uint16_t>(code, sx, sy);
	}
	else if (byte_width == 4) {
		return unpack_binary_helper<uint32_t>(code, sx, sy);
	}
	else {
		return unpack_binary_helper<uint64_t>(code, sx, sy);
	}
}

std::vector<uint8_t> decode_permissible_crack_code(
	const std::unordered_map<uint64_t, std::vector<unsigned char>> &chains,
	const uint64_t sx, const uint64_t sy
	) {
	// voxel connectivity
	// four bits: -y-x+y+x true is passable
	std::vector<uint8_t> edges(sx * sy);

	uint64_t sxe = sx + 1;

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		uint64_t y = node / sxe;
		uint64_t x = node - (sxe * y);

		std::stack<uint64_t> revisit;
		for (unsigned char symbol : symbols) {
			uint64_t loc = x + sx * y;
			if (symbol == 'u') {
				if ((x-1) < sx && (y-1) < sy) {
					edges[loc - 1 - sx] |= 0b0001;
				}
				if (y-1 < sy) {
					edges[loc - sx] |= 0b0010;
				}
				y--;
			}
			else if (symbol == 'd') {
				if ((x-1) < sx) {
					edges[loc-1] |= 0b0001;
				}
				edges[loc] |= 0b0010;
				y++;
			}
			else if (symbol == 'l') {
				if ((x-1) < sx && (y-1) < sy) {
					edges[loc-1-sx] |= 0b0100;
				}
				if (x-1 < sx) {
					edges[loc-1] |= 0b1000;
				}
				x--;
			}
			else if (symbol == 'r') {
				if (y-1 < sy) {
					edges[loc-sx] |= 0b0100;
				}
				edges[loc] |= 0b1000;
				x++;
			}
			else if (symbol == 'b') {
				revisit.push(loc);
			}
			else if (symbol =='t') {
				if (!revisit.empty()) {
					loc = revisit.top();
					revisit.pop();
					y = loc / sxe;
					x = loc - (sxe * y);
				}
			}
		}
	}

	return edges;
}

std::vector<uint8_t> decode_impermissible_crack_code(
	const std::unordered_map<uint64_t, std::vector<unsigned char>> &chains,
	const uint64_t sx, const uint64_t sy
	) {
	// voxel connectivity
	// four bits: -y-x+y+x true is passable
	std::vector<uint8_t> edges(sx * sy);
	for (uint64_t i = 0; i < sx * sy; i++) {
		edges[i] = 0b1111;
	}

	uint64_t sxe = sx + 1;

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		uint64_t y = node / sxe;
		uint64_t x = node - (sxe * y);

		std::stack<uint64_t> revisit;
		for (unsigned char symbol : symbols) {
			uint64_t loc = x + sx * y;
			if (symbol == 'u') {
				if ((x-1) < sx && (y-1) < sy) {
					edges[loc - 1 - sx] &= 0b1110;
				}
				if (y-1 < sy) {
					edges[loc - sx] &= 0b1101;
				}
				y--;
			}
			else if (symbol == 'd') {
				if ((x-1) < sx) {
					edges[loc-1] &= 0b1110;
				}
				edges[loc] &= 0b1101;
				y++;
			}
			else if (symbol == 'l') {
				if ((x-1) < sx && (y-1) < sy) {
					edges[loc-1-sx] &= 0b1011;
				}
				if (x-1 < sx) {
					edges[loc-1] &= 0b0111;
				}
				x--;
			}
			else if (symbol == 'r') {
				if (y-1 < sy) {
					edges[loc-sx] &= 0b1011;
				}
				edges[loc] &= 0b0111;
				x++;
			}
			else if (symbol == 'b') {
				revisit.push(loc);
			}
			else if (symbol =='t') {
				if (!revisit.empty()) {
					loc = revisit.top();
					revisit.pop();
					y = loc / sxe;
					x = loc - (sxe * y);
				}
			}
		}
	}

	return edges;
}

std::vector<uint8_t> decode_crack_code(
	const std::unordered_map<uint64_t, std::vector<unsigned char>> &chains,
	const uint64_t sx, const uint64_t sy,
	const bool permissible
) {
	if (permissible) {
		return decode_permissible_crack_code(chains, sx, sy);
	}
	else {
		return decode_impermissible_crack_code(chains, sx, sy);
	}
}

};
};

#endif