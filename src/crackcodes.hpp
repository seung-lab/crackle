#ifndef __CRACKLE_CRACKCODE_HXX__
#define __CRACKLE_CRACKCODE_HXX__

#include <vector>
#include <unordered_map>
#include <stack>
#include <span>
#include <cstdint>
#include <mutex>

#include "robin_hood.hpp"
#include "lib.hpp"
#include "cc3d.hpp"
#include "builtins.hpp"
#include "threadpool.hpp"

namespace crackle {
namespace crackcodes {

enum DirectionCode {
	LEFT = 0b11,
	RIGHT = 0b01,
	UP = 0b00,
	DOWN = 0b10,
	NONE = 255
};

inline std::pair<int64_t, int64_t> mkedge(int64_t a, int64_t b){
	if (a < b) {
		return std::make_pair(a,b);
	}
	return std::make_pair(b,a);
};

struct Graph {
	std::vector<uint8_t> adjacency;
	
	int64_t sxe;
	int64_t sye;

	int64_t next_cluster(int64_t idx) {
		for (int64_t i = idx; i < sxe * sye; i++) {
			if (adjacency[i]) {
				return i;
			}
		}

		return -1;
	}

	void erase_edge(const std::pair<int64_t, int64_t> &edge) {
		int dir = edge.second - edge.first;
		// b/c second is always > first,
		// the direction will either be pointing 
		// right or down.
		if (dir == 1) {
			adjacency[edge.first] &= 0b1110;
			adjacency[edge.second] &= 0b1101;
		}
		else {
			adjacency[edge.first] &= 0b1011;
			adjacency[edge.second] &= 0b0111;	
		}
	}

	template <typename LABEL>
	bool init(
		const LABEL* labels,
		const int64_t sx, const int64_t sy,
		const bool permissible
	) {
		sxe = sx + 1; // sx edges
		sye = sy + 1; // sy edges

		adjacency.resize(sxe * sye);

		bool any_edges = false;

		if (permissible) {
			for (int64_t y = 0; y < sy; y++) {
				for (int64_t x = 0; x < sx; x++) {
					// assign vertical edges
					if (x > 0 && labels[x + sx * y] == labels[(x-1) + sx * y]) {
						int64_t node_up = x + sxe * y;
						int64_t node_down = x + sxe * (y + 1);
						adjacency[node_up] |= 0b0100;
						adjacency[node_down] |= 0b1000;
						any_edges = true;
					}
					// assign horizontal edges
					if (y > 0 && labels[x + sx * y] == labels[x + sx * (y-1)]) {
						int64_t node_left = x + sxe * y;
						int64_t node_right = (x+1) + sxe * y;
						adjacency[node_left] |= 0b0001;
						adjacency[node_right] |= 0b0010;
						any_edges = true;
					}
				}
			}
		}
		else {
			for (int64_t y = 0; y < sy; y++) {
				for (int64_t x = 0; x < sx; x++) {
					// assign vertical edges
					if (x > 0 && labels[x + sx * y] != labels[(x-1) + sx * y]) {
						int64_t node_up = x + sxe * y;
						int64_t node_down = x + sxe * (y + 1);
						adjacency[node_up] |= 0b0100;
						adjacency[node_down] |= 0b1000;
						any_edges = true;
					}
					// assign horizontal edges
					if (y > 0 && labels[x + sx * y] != labels[x + sx * (y-1)]) {
						int64_t node_left = x + sxe * y;
						int64_t node_right = (x+1) + sxe * y;
						adjacency[node_left] |= 0b0001;
						adjacency[node_right] |= 0b0010;
						any_edges = true;
					}
				}
			}
		}

		return any_edges;
	}
};

robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>
symbols_to_codepoints(
	std::vector<std::pair<uint64_t, std::vector<unsigned char>>> &chains
) {
	robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>> encoded_chains;

	const uint64_t BRANCH[2] = { DirectionCode::UP, DirectionCode::DOWN };
	const uint64_t BRANCH2[2] = { DirectionCode::LEFT, DirectionCode::RIGHT };
	const uint64_t TERM[2] = { DirectionCode::DOWN, DirectionCode::UP };
	const uint64_t TERM2[2] = { DirectionCode::RIGHT, DirectionCode::LEFT };

	constexpr uint8_t lookup_size = static_cast<uint8_t>('u') + 1;
	uint8_t lookup[lookup_size];
	lookup[static_cast<uint8_t>('u')] = DirectionCode::UP;
	lookup[static_cast<uint8_t>('d')] = DirectionCode::DOWN;
	lookup[static_cast<uint8_t>('l')] = DirectionCode::LEFT;
	lookup[static_cast<uint8_t>('r')] = DirectionCode::RIGHT;

	for (auto [node, chain] : chains) {
		std::vector<uint8_t> code;
		code.reserve(chain.size() * 3 / 2); // account for b and t

		for (uint64_t i = 0; i < chain.size(); i++) {
			char symbol = chain[i];
			if (symbol == 's') {
				continue;
			}
			else if (symbol == 'b') {
				if (i > 0 && code.back() != BRANCH[1]) {
					code.push_back(BRANCH[0]);
					code.push_back(BRANCH[1]);
				}
				else {
					code.push_back(BRANCH2[0]);
					code.push_back(BRANCH2[1]);	
				}
			}
			else if (symbol == 't') {
				if (i > 0 && code.back() != TERM[1]) {
					code.push_back(TERM[0]);
					code.push_back(TERM[1]);
				}
				else {
					code.push_back(TERM2[0]);
					code.push_back(TERM2[1]);	
				}
			}
			else {
				code.push_back(lookup[static_cast<uint8_t>(symbol)]);
			}
		}
		encoded_chains[node] = std::move(code);
	}

	return encoded_chains;
}

int64_t remove_initial_branch(
	int64_t node,
	std::vector<unsigned char>& code,
	const int64_t sx, const int64_t sy
) {
	if (code.empty()) {
		return node;
	}
	else if (code[0] != 'b') {
		return node;
	}

	int64_t i = 1;
	while (code[i] != 't') {
		if (code[i] == 'b') {
			return node;
		}
		i++;
	}
	
	int64_t sxe = sx + 1;
	int64_t y = node / sxe;
	int64_t x = node - (sxe * y);

	robin_hood::unordered_flat_map<char,char> flip({
		{'u', 'd'},
		{'d', 'u'},
		{'l', 'r'},
		{'r', 'l'},
		{'s', 's'},
	});

	robin_hood::unordered_flat_map<char,std::pair<int,int>> mvmt({
		{'u', std::make_pair(0,-1)},
		{'d', std::make_pair(0,+1)},
		{'l', std::make_pair(-1,0)},
		{'r', std::make_pair(+1,0)},
		{'s', std::make_pair(0,0)},
	});

	code[0] = 's';
	int64_t pos_x = x;
	int64_t pos_y = y;

	for (i = 1; code[i] != 't'; i++) {
		pos_x += mvmt[code[i]].first;
		pos_y += mvmt[code[i]].second;
		code[i] = flip[code[i]];
	}
	code[i] = 's';
	const int64_t last = i - 1;
	i--;
	for (; i > last/2; i--) {
		std::swap(code[i], code[last - i + 1]);
	}

	return pos_x + sxe * pos_y;
}

// during code creation, if the path loops and passed
// through a previously marked branch point, a spurious
// set of b/t is created. Previously, we tried removing
// these during creation, but this required adding hash
// access that was slower and higher memory. So let's just
// strip them afterwards.
void remove_spurious_branches(
	std::vector<unsigned char>& code
) {
	std::vector<int64_t> branch_stack;
	branch_stack.push_back(-1);

	std::vector<uint32_t> branch_lens(code.size() + 1);

	std::vector<std::pair<int64_t,int64_t>> to_erase;

	int64_t current_branch = -1;
	for (int64_t i = 0; i < static_cast<int64_t>(code.size()); i++) {
		if (code[i] == 'b') {
			branch_stack.push_back(i);
		}
		else if (code[i] == 't') {
			if (current_branch >= 0 && branch_lens[current_branch+1] == 0) {
				to_erase.emplace_back(current_branch, i);
			}
			current_branch = branch_stack.back();
			branch_stack.pop_back();
		}
		else {
			branch_lens[current_branch+1]++;
		}
	}

	for (auto pair : to_erase) {
		code[pair.first] = 's';
		code[pair.second] = 's';
	}
}

std::vector<uint64_t> read_boc_index(
	const std::span<const unsigned char>& binary,
	const uint64_t sx, const uint64_t sy
) {
	std::vector<uint64_t> nodes;

	const uint64_t sxe = sx + 1;

	const uint64_t x_width = crackle::lib::compute_byte_width(sx+1);
	const uint64_t y_width = crackle::lib::compute_byte_width(sy+1);

	uint64_t idx = 4; // skip over index size
	uint64_t num_y = crackle::lib::ctoid(binary, idx, y_width);
	idx += y_width;

	uint64_t y = 0; 

	for (uint64_t yi = 0; yi < num_y; yi++) {
		y += crackle::lib::ctoid(binary, idx, y_width);
		idx += y_width;

		uint64_t num_x = crackle::lib::ctoid(binary, idx, x_width);
		idx += x_width;

		uint64_t x = 0;
		for (uint64_t xi = 0; xi < num_x; xi++) {
			x += crackle::lib::ctoid(binary, idx, x_width);
			idx += x_width;
			nodes.push_back(x + sxe * y);
		}
	}

	return nodes;
}

std::vector<unsigned char> write_boc_index(
	const std::vector<uint64_t>& sorted_nodes,
	const uint64_t sx, const uint64_t sy
) {
	const uint64_t sxe = sx + 1;
	const uint64_t sye = sy + 1;

	const uint64_t x_width = crackle::lib::compute_byte_width(sx+1);
	const uint64_t y_width = crackle::lib::compute_byte_width(sy+1);

	// beginning of chain index
	robin_hood::unordered_node_map<uint64_t, std::vector<uint64_t>> boc(sye);

	for (uint64_t node : sorted_nodes) {
		uint64_t y = node / sxe;
		uint64_t x = node - sxe * y;
		boc[y].push_back(x);
	}

	std::vector<uint64_t> all_y;
	for (auto pair : boc) {
		all_y.push_back(pair.first);
	}
	std::sort(all_y.begin(), all_y.end());

	uint64_t index_size = y_width;
	for (uint64_t y : all_y) {
		index_size += y_width;
		index_size += (boc[y].size() + 1) * x_width;
	}

	std::vector<unsigned char> binary(4 + index_size);

	uint64_t idx = 0;
	idx += crackle::lib::itocd(index_size, binary, idx, 4);
	idx += crackle::lib::itocd(all_y.size(), binary, idx, y_width);
	
	for (uint64_t i = 0; i < all_y.size(); i++) {
		if (i == 0) {
			idx += crackle::lib::itocd(all_y[0], binary, idx, y_width);
		}
		else {
			idx += crackle::lib::itocd(all_y[i] - all_y[i-1], binary, idx, y_width);
		}
		uint64_t y = all_y[i];
		idx += crackle::lib::itocd(boc[y].size(), binary, idx, x_width);
		uint64_t last_x = 0;
		for (uint64_t x : boc[y]) {
			idx += crackle::lib::itocd(x - last_x, binary, idx, x_width);
			last_x = x;
		}
	}

	return binary;
}

template <typename LABEL>
robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>
create_crack_codes(
	const LABEL* labels,
	const int64_t sx, const int64_t sy,
	bool permissible
) {
	Graph G;
	const bool any_edges = G.init(labels, sx, sy, permissible);

	std::vector<std::pair<uint64_t, std::vector<unsigned char>>> chains;

	if (!any_edges) {
		return symbols_to_codepoints(chains);
	}

	std::vector<int64_t> revisit;
	revisit.reserve(sx);
	std::vector<uint8_t> revisit_ct((sx+1)*(sy+1));

	int64_t start_node = 0;

	const int64_t lookup_dir[4] = { 1, -1, sx+1, -(sx+1) };
	constexpr char lookup_symbol[4] = { 'r', 'l', 'd', 'u' };

	while ((start_node = G.next_cluster(start_node)) != -1) {
		int64_t node = start_node;

		std::vector<unsigned char> code;
		code.reserve(sx >> 2);

		int64_t branches_taken = 1;

		while (G.adjacency[node] || !revisit.empty()) {
			if (!G.adjacency[node]) {
				code.push_back('t');
				branches_taken--;
				while (!revisit.empty()) {
					node = revisit.back();
					revisit.pop_back();
					revisit_ct[node]--;
					break;
				}
				continue;
			}
			else if (popcount(G.adjacency[node]) > 1) {
				code.push_back('b');
				revisit.push_back(node);
				revisit_ct[node]++;
				branches_taken++;
			}

			const int next_neighbor_idx = ctz(G.adjacency[node]);
			const int64_t dir_taken = lookup_dir[next_neighbor_idx];
			const int64_t next_node = node + dir_taken;
			code.push_back(lookup_symbol[next_neighbor_idx]);

			auto edge = mkedge(node, next_node);
			G.erase_edge(edge);
			node = next_node;
		}

		while (branches_taken > 0) {
			code.push_back('t');
			branches_taken--;
		}

		const int64_t adjusted_start_node = remove_initial_branch(
			start_node, code, sx, sy
		);

		remove_spurious_branches(code);

		chains.push_back(
			std::make_pair(static_cast<uint64_t>(adjusted_start_node), code)
		);
	}

	return symbols_to_codepoints(chains);
}

std::vector<unsigned char> pack_codepoints(
	robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>& chains,
	const uint64_t sx, const uint64_t sy
) {

	std::vector<uint64_t> nodes;
	for (auto& [node, code] : chains) {
		nodes.push_back(node);
	}
	std::sort(nodes.begin(), nodes.end());

	std::vector<unsigned char> binary = write_boc_index(nodes, sx, sy);

	std::vector<uint8_t> codepoints;
	uint8_t last = 0;
	for (uint64_t node : nodes) {
		auto chain = chains[node];
		for (uint8_t codepoint : chain) {
			uint8_t diffcoded = (codepoint - last) & 0b11;
			last = codepoint;
			codepoints.push_back(diffcoded);
		}
	}

	uint8_t encoded = 0;
	int pos = 0;
	for (uint64_t i = 0; i < codepoints.size(); i++) {
		encoded |= (codepoints[i] << pos);
		pos += 2;
		if (pos == 8) {
			binary.push_back(encoded);
			encoded = 0;
			pos = 0;
		}
	}

	if (pos > 0) {
		binary.push_back(encoded);
	}

	return binary;
}

template <typename LABEL>
std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> 
encode_boundaries(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool permissible,
	const size_t parallel
) {
	std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> binary_components(sz);

	const int64_t sxy = sx * sy;

	ThreadPool pool(parallel);

	for (int64_t z = 0; z < sz; z++) {
		pool.enqueue([&,z](size_t t){
			binary_components[z] = create_crack_codes(labels + z * sxy, sx, sy, permissible);
		});
	}

	pool.join();

	return binary_components;
}

std::vector<std::pair<uint64_t, std::vector<unsigned char>> >
codepoints_to_symbols(
	const std::vector<uint64_t>& sorted_nodes,
	const std::vector<uint8_t>& codepoints
) {

	std::vector<std::pair<uint64_t, std::vector<unsigned char>> > chains;

	std::vector<unsigned char> symbols;
	symbols.reserve(codepoints.size() * 4 * 2);

	uint64_t branches_taken = 0;
	uint64_t node = 0;

	constexpr char remap[4] = { 'u', 'r', 'd', 'l' };

	uint64_t node_i = 0;

	uint8_t last_move = DirectionCode::NONE;

	for (uint64_t i = 0; i < codepoints.size(); i++) {
		if (branches_taken == 0) {
			if (node_i >= sorted_nodes.size()) {
				break;
			}
			node = sorted_nodes[node_i];
			node_i++;
			i--; // b/c i will be incremented
			branches_taken = 1;
			continue;
		}

		uint8_t move = codepoints[i];
	
		// by chance, up^down and left^right 
		// both evaluate to 0b10
		if ((move ^ last_move) != 0b10) {
			symbols.push_back(remap[move]);
			last_move = move;
			continue;
		}
		else if (
			// equivalent to:
			// move == DirectionCode::UP || move == DirectionCode::LEFT
			// 
			// which is equivalent to (because we already check 
			// against last_move in move ^ last_move = 0b10) which
			// means last move is guaranteed to be its opposite.
			//
			// (move == DirectionCode::UP && last_move == DirectionCode::DOWN)
			// || (move == DirectionCode::LEFT && last_move == DirectionCode::RIGHT)
			popcount(move) != 1 // 00 (LEFT) or 11 (UP), 7 operations -> 2
		) {
			symbols.back() = 't';
			branches_taken--;
			last_move = DirectionCode::NONE;
		}
		else { // the code here is DOWN+UP or RIGHT+LEFT
			symbols.back() = 'b';
			branches_taken++;
			last_move = DirectionCode::NONE;
		}

		if (branches_taken == 0) {
			chains.push_back(std::make_pair(node, symbols));
			symbols.clear();
		}
	}

	return chains;
}

std::vector<uint8_t> unpack_codepoints(
	const std::span<const unsigned char> &code, 
	const uint64_t sx, const uint64_t sy
) {
	if (code.size() == 0) {
		return std::vector<uint8_t>();
	}

	uint32_t index_size = 4 + crackle::lib::ctoid(code, 0, 4);

	std::vector<uint8_t> codepoints;
	codepoints.reserve(4 * (code.size() - index_size));

	uint8_t last = 0;

	for (uint64_t i = index_size; i < code.size(); i++) {
		for (uint64_t j = 0; j < 4; j++) {
			uint8_t codepoint = static_cast<uint8_t>((code[i] >> (2*j)) & 0b11);
			codepoint += last;
			codepoint &= 0b11;
			last = codepoint;
			codepoints.push_back(codepoint);
		}
	}

	return codepoints;
}

void decode_permissible_crack_code(
	const std::vector<std::pair<uint64_t, std::vector<unsigned char>> > &chains,
	const int64_t sx, const int64_t sy,
	uint8_t* edges
) {
	// voxel connectivity matches cc3d_graphs.hpp 4 connected
	// four bits: -y+y-x+x true is passable
	std::fill(edges, edges + sx * sy, 0);

	const int64_t sxe = sx + 1;

	const uint64_t pixels = (sx+1) * (sy+1);

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		int64_t y = node / sxe;
		int64_t x = node - (sxe * y);
		uint64_t loc = x + sx * y;

		std::stack<int64_t> revisit;
		for (unsigned char symbol : symbols) {
			if (loc >= pixels) {
				std::string err = std::string("crackle: decode_permissible_crack_code: index out of range. loc: ");
				err += std::to_string(loc);
				throw std::runtime_error(err);
			}

			if (symbol == 'u') {
				if (x > 0 && y > 0) {
					edges[loc - 1 - sx] |= 0b0001;
				}
				if (y > 0) {
					edges[loc - sx] |= 0b0010;
				}
				y--;
				loc -= sx;
			}
			else if (symbol == 'd') {
				if (x > 0) {
					edges[loc - 1] |= 0b0001;
				}
				edges[loc] |= 0b0010;
				y++;
				loc += sx;
			}
			else if (symbol == 'l') {
				if (x > 0 && y > 0) {
					edges[loc-1-sx] |= 0b0100;
				}
				if (x > 0) {
					edges[loc-1] |= 0b1000;
				}
				x--;
				loc--;
			}
			else if (symbol == 'r') {
				if (y > 0) {
					edges[loc-sx] |= 0b0100;
				}
				edges[loc] |= 0b1000;
				x++;
				loc++;
			}
			else if (symbol == 'b') {
				revisit.push(loc);
			}
			else if (symbol =='t') {
				if (!revisit.empty()) {
					loc = revisit.top();
					revisit.pop();
					y = loc / sx;
					x = loc - (sx * y);
				}
			}
		}
	}
}

void decode_impermissible_crack_code(
	const std::vector<std::pair<uint64_t, std::vector<unsigned char>> > &chains,
	const int64_t sx, const int64_t sy,
	uint8_t* edges
) {
	// voxel connectivity matches cc3d_graphs.hpp 4 connected
	// four bits: -y+y-x+x true is passable
	std::fill(edges, edges + sx * sy, 0b1111);

	const int64_t sxe = sx + 1;

	const uint64_t pixels = (sx+1) * (sy+1);

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		int64_t y = node / sxe;
		int64_t x = node - (sxe * y);
		uint64_t loc = x + sx * y;

		std::stack<int64_t> revisit;
		for (unsigned char symbol : symbols) {
			if (loc >= pixels) {
				throw std::runtime_error("crackle: decode_impermissible_crack_code: index out of range.");
			}

			if (symbol == 'u') {
				if (x > 0 && y > 0) {
					edges[loc - 1 - sx] &= 0b1110;
				}
				if (y > 0) {
					edges[loc - sx] &= 0b1101;
				}
				y--;
				loc -= sx;
			}
			else if (symbol == 'd') {
				if (x > 0) {
					edges[loc - 1] &= 0b1110;
				}
				edges[loc] &= 0b1101;
				y++;
				loc += sx;
			}
			else if (symbol == 'l') {
				if (x > 0 && y > 0) {
					edges[loc - 1 - sx] &= 0b1011;
				}
				if (x > 0) {
					edges[loc-1] &= 0b0111;
				}
				x--;
				loc--;
			}
			else if (symbol == 'r') {
				if (y > 0) {
					edges[loc - sx] &= 0b1011;
				}
				edges[loc] &= 0b0111;
				x++;
				loc++;
			}
			else if (symbol == 'b') {
				revisit.push(loc);
			}
			else if (symbol =='t') {
				if (!revisit.empty()) {
					loc = revisit.top();
					revisit.pop();
					y = loc / sx;
					x = loc - (sx * y);
				}
			}
		}
	}
}

void decode_crack_code(
	const std::vector<std::pair<uint64_t, std::vector<unsigned char>> > &chains,
	const uint64_t sx, const uint64_t sy,
	const bool permissible, 
	uint8_t* slice_edges
) {
	if (permissible) {
		decode_permissible_crack_code(chains, sx, sy, slice_edges);
	}
	else {
		decode_impermissible_crack_code(chains, sx, sy, slice_edges);
	}
}

};
};

#endif