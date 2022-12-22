#ifndef __CRACKLE_CRACKCODE_HXX__
#define __CRACKLE_CRACKCODE_HXX__

#include <vector>
#include <unordered_map>
#include <stack>
#include <cstdint>

#include "robin_hood.hpp"
#include "lib.hpp"
#include "cc3d.hpp"

namespace crackle {
namespace crackcodes {

enum DirectionCode {
	LEFT = 0b10,
	RIGHT = 0b01,
	UP = 0b00,
	DOWN = 0b11
};

std::pair<int64_t, int64_t> mkedge(int64_t a, int64_t b){
	if (a < b) {
		return std::make_pair(a,b);
	}
	return std::make_pair(b,a);
};

struct Graph {
	std::vector<uint8_t> adjacency;
	std::vector<std::vector<std::pair<int64_t, int64_t>>> component_edge_list;
	int64_t sxe;
	int64_t sye;

	std::vector<int64_t> neighbors(int64_t node) {
		std::vector<int64_t> nbrs;
		if (adjacency[node] & 0b0001) {
			nbrs.push_back(node + 1);
		}
		if (adjacency[node] & 0b0010) {
			nbrs.push_back(node - 1);
		}
		if (adjacency[node] & 0b0100) {
			nbrs.push_back(node + sxe);
		}
		if (adjacency[node] & 0b1000) {
			nbrs.push_back(node - sxe);
		}
		return nbrs;
	}

	int64_t num_components() {
		return component_edge_list.size();
	}

	template <typename LABEL>
	void init(
		const LABEL* labels,
		const int64_t sx, const int64_t sy,
		const bool permissible
	) {
		sxe = sx + 1; // sx edges
		sye = sy + 1; // sy edges

		adjacency.resize(sxe * sye);

		crackle::cc3d::DisjointSet<uint32_t> equivalences(sxe * sye);
		for (uint64_t i = 0; i < sxe * sye; i++) {
			equivalences.ids[i] = i;
		}

		std::vector<std::pair<int64_t, int64_t>> all_edges;

		auto check = [](int a, int b) { return a != b; };
		if (permissible) {
			check = [](int a, int b) { return a == b; };
		}

		// assign vertical edges
		for (uint64_t y = 0; y < sy; y++) {
			for (uint64_t x = 1; x < sx; x++) {
				if (check(labels[x + sx * y], labels[(x-1) + sx * y])) {
					int64_t node_up = x + sxe * y;
					int64_t node_down = x + sxe * (y + 1);
					adjacency[node_up] |= 0b0100;
					adjacency[node_down] |= 0b1000;
					equivalences.unify(node_up, node_down);
					all_edges.push_back(mkedge(node_up, node_down));
				}
			}
		}

		// assign horizontal edges
		for (uint64_t y = 1; y < sy; y++) {
			for (uint64_t x = 0; x < sx; x++) {
				if (check(labels[x + sx * y], labels[x + sx * (y-1)])) {
					int64_t node_left = x + sxe * y;
					int64_t node_right = (x+1) + sxe * y;
					adjacency[node_left] |= 0b0001;
					adjacency[node_right] |= 0b0010;
					equivalences.unify(node_left, node_right);
					all_edges.push_back(mkedge(node_left, node_right));
				}
			}
		}

		component_edge_list.resize(all_edges.size() * 2);

		robin_hood::unordered_flat_map<int64_t,int64_t> renumber;
		int64_t next_label = 1;
		int64_t label = 0;
		int64_t membership = 0;
		for (auto pair : all_edges) {
			label = equivalences.root(pair.first);
			if (renumber[label] == 0) {
				renumber[label] = next_label;
				membership = next_label;
				next_label++;
			}
			else {
				membership = renumber[label];
			}

			component_edge_list[membership].push_back(pair);
		}

		component_edge_list.resize(next_label);
	}
};

std::vector<std::vector<uint64_t>>
symbols_to_integers(
	std::vector<std::pair<int64_t, std::vector<char>>> &chains
) {
	std::vector<std::vector<uint64_t>> encoded_chains;

	const uint64_t BRANCH[2] = { DirectionCode::UP, DirectionCode::DOWN };
	const uint64_t BRANCH2[2] = { DirectionCode::LEFT, DirectionCode::RIGHT };
	const uint64_t TERM[2] = { DirectionCode::DOWN, DirectionCode::UP };
	const uint64_t TERM2[2] = { DirectionCode::RIGHT, DirectionCode::LEFT };

	for (auto [node, chain] : chains) {
		std::vector<uint64_t> code;
		code.reserve(chain.size());
		code.push_back(node);
		for (int64_t i = 0; i < chain.size(); i++) {
			char symbol = chain[i];
			if (symbol == 's') {
				continue;
			}
			else if (symbol == 'u') {
				code.push_back(DirectionCode::UP);
			}
			else if (symbol == 'd') {
				code.push_back(DirectionCode::DOWN);
			}
			else if (symbol == 'l') {
				code.push_back(DirectionCode::LEFT);
			}
			else if (symbol == 'r') {
				code.push_back(DirectionCode::RIGHT);
			}
			else if (symbol == 'b') {
				if (i > 0 && chain[i-1] != BRANCH[1]) {
					code.push_back(BRANCH[0]);
					code.push_back(BRANCH[1]);
				}
				else {
					code.push_back(BRANCH2[0]);
					code.push_back(BRANCH2[1]);	
				}
			}
			else if (symbol == 't') {
				if (i > 0 && chain[i-1] != TERM[1]) {
					code.push_back(TERM[0]);
					code.push_back(TERM[1]);
				}
				else {
					code.push_back(TERM2[0]);
					code.push_back(TERM2[1]);	
				}
			}
			else {
				throw std::runtime_error("Invalid symbol.");
			}
		}
		encoded_chains.push_back(code);
	}

	return encoded_chains;
}

void remove_initial_branch(
	int64_t& node,
	std::vector<char>& code,
	const int64_t sx, const int64_t sy
) {
	if (code.empty()) {
		return;
	}
	else if (code[0] != 'b') {
		return;
	}

	int64_t i = 1;
	while (code[i] != 't') {
		if (code[i] == 'b') {
			return;
		}
		i++;
	}
	
	int64_t sxe = sx + 1;
	int64_t y = node / sxe;
	int64_t x = node - (sxe * y);
	// int64_t pos 

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
	const int64_t len = i;
	for (; i >= i/2; i--) {
		std::swap(code[i], code[len - i]);
	}

	node = pos_x + sxe * pos_y;
}

template <typename LABEL>
std::vector<std::vector<uint64_t>>
create_crack_codes(
	const LABEL* labels,
	const int64_t sx, const int64_t sy,
	bool permissible
) {
  Graph G;
  G.init(labels, sx, sy, permissible);

  const int64_t n_clusters = G.num_components();

  std::unordered_map<int64_t, char> dirmap = {
    {1, 'r'},
    {-1, 'l'},
    {sx+1, 'd'},
    {-(sx+1), 'u'}
  };

  std::vector<std::pair<int64_t, std::vector<char>>> chains;
  std::vector<int64_t> revisit;

  if (n_clusters == 0) {
    return symbols_to_integers(chains);
  }

  for (int64_t cluster = 0; cluster < n_clusters; cluster++) {
  	robin_hood::unordered_flat_set<std::pair<int64_t, int64_t>>
  		remaining;

  	auto& cluster_edges = G.component_edge_list[cluster];
  	remaining.insert(cluster_edges.begin(), cluster_edges.end());
  	
  	std::pair<int64_t, int64_t> start_edge = *remaining.begin();

  	int64_t node = start_edge.first;
  	int64_t start_node = start_edge.first;
  	remaining.erase(start_edge);

    std::vector<char> code;
    std::vector<int64_t> revisit;
    std::unordered_map<int64_t, std::vector<int64_t>> branch_nodes;
    int64_t branches_taken = 1;

    while (!remaining.empty() || !revisit.empty()) {
    	auto neighbors = G.neighbors(node);

    	if (neighbors.empty()) {
    		code.push_back('t');
    		branches_taken--;
    		if (!revisit.empty()) {
    			node = revisit.back();
    			revisit.pop_back();
    		}
    		else if (!remaining.empty()) {
    			node = (*remaining.begin()).first;
    		}
    		continue;
    	}
    	else if (neighbors.size() > 1) {
    		code.push_back('b');
    		revisit.push_back(node);
    		branch_nodes[node].push_back(code.size() - 1);
    		branches_taken++;
    	}

    	int64_t next_node = neighbors[0];
    	int64_t dir_taken = dirmap[next_node - node];
    	code.push_back(dir_taken);
    	remaining.erase(mkedge(node, next_node));
    	node = next_node;

    	if (revisit.count(node)) {
    		int64_t pos = 0;
				for (pos = revisit.size() - 1; pos >= 0; pos--) {
					if (revisit[pos] == node) {
						break;
					}
				}
				revisit.erase(pos);
				branches_taken--;
				code[branch_nodes[node].back()] = 's';
				branch_nodes[node].pop_back();
    	}
    }

    while (branches_taken > 0) {
    	code.push_back('t');
    	branches_taken--;
    }

    remove_initial_branch(start_node, code, sx, sy);
    chains.push_back(
    	std::make_pair(start_node, code)
    );
  }

  return symbols_to_integers(chains);
}

template <typename T>
std::unordered_map<uint64_t, std::vector<unsigned char>> 
unpack_binary_helper(
	const std::vector<unsigned char> &code, 
	const uint64_t sx, const uint64_t sy
) {
	std::unordered_map<uint64_t, std::vector<unsigned char>> chains;
	uint64_t code_size = code.size() / sizeof(T);
	std::vector<T> int_code(code_size);
	for (uint64_t i = 0; i < code_size; i++) {
		int_code[i] = crackle::lib::ctoi<T>(code.data(), i * sizeof(T));
	}

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
					branches_taken--;
					if (branches_taken == 0) {
						break;
					}
				}
				else if (
					(move == 3 && symbols.back() == 'u')
					|| (move == 1 && symbols.back() == 'l')
				) {
					symbols.back() = 'b';
					branches_taken++;
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
			chains[node] = vec;
			symbols.clear();
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
	const int64_t sx, const int64_t sy
) {
	// voxel connectivity
	// four bits: -y-x+y+x true is passable
	std::vector<uint8_t> edges(sx * sy);

	int64_t sxe = sx + 1;

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		int64_t y = node / sxe;
		int64_t x = node - (sxe * y);

		std::stack<int64_t> revisit;
		for (unsigned char symbol : symbols) {
			int64_t loc = x + sx * y;
			if (loc < 0 || loc >= (sx+1) * (sy+1)) {
				throw std::runtime_error("crackle: decode_permissible_crack_code: index out of range.");
			}

			if (symbol == 'u') {
				if (x > 0 && y > 0) {
					edges[loc - 1 - sx] |= 0b0001;
				}
				if (y > 0) {
					edges[loc - sx] |= 0b0010;
				}
				y--;
			}
			else if (symbol == 'd') {
				if (x > 0) {
					edges[loc - 1] |= 0b0001;
				}
				edges[loc] |= 0b0010;
				y++;
			}
			else if (symbol == 'l') {
				if (x > 0 && y > 0) {
					edges[loc-1-sx] |= 0b0100;
				}
				if (x > 0) {
					edges[loc-1] |= 0b1000;
				}
				x--;
			}
			else if (symbol == 'r') {
				if (y > 0) {
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
					y = loc / sx;
					x = loc - (sx * y);
				}
			}
		}
	}

	return edges;
}

std::vector<uint8_t> decode_impermissible_crack_code(
	const std::unordered_map<uint64_t, std::vector<unsigned char>> &chains,
	const int64_t sx, const int64_t sy
) {
	// voxel connectivity
	// four bits: -y-x+y+x true is passable
	std::vector<uint8_t> edges(sx * sy);
	for (int64_t i = 0; i < sx * sy; i++) {
		edges[i] = 0b1111;
	}

	int64_t sxe = sx + 1;

	// graph is of corners and edges
	// origin is located at top left
	// corner of the image
	for (auto& [node, symbols]: chains) {
		int64_t y = node / sxe;
		int64_t x = node - (sxe * y);

		std::stack<int64_t> revisit;
		for (unsigned char symbol : symbols) {
			int64_t loc = x + sx * y;

			if (loc < 0 || loc >= (sx+1) * (sy+1)) {
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
			}
			else if (symbol == 'd') {
				if (x > 0) {
					edges[loc - 1] &= 0b1110;
				}
				edges[loc] &= 0b1101;
				y++;
			}
			else if (symbol == 'l') {
				if (x > 0 && y > 0) {
					edges[loc - 1 - sx] &= 0b1011;
				}
				if (x > 0) {
					edges[loc-1] &= 0b0111;
				}
				x--;
			}
			else if (symbol == 'r') {
				if (y > 0) {
					edges[loc - sx] &= 0b1011;
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
					y = loc / sx;
					x = loc - (sx * y);
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