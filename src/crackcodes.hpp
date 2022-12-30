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

inline std::pair<int64_t, int64_t> mkedge(int64_t a, int64_t b){
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

	void neighbors(int64_t node, std::vector<int64_t> &nbrs) {
		nbrs.clear();
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
	}

	std::vector<int64_t> neighbors(int64_t node) {
		std::vector<int64_t> nbrs;
		nbrs.reserve(4);
		neighbors(node, nbrs);
		return nbrs;
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
		for (int64_t i = 0; i < sxe * sye; i++) {
			equivalences.ids[i] = i;
		}

		std::vector<std::pair<int64_t, int64_t>> all_edges;
		all_edges.reserve(sxe * sye / 10);

		if (permissible) {
			// assign vertical edges
			for (int64_t y = 0; y < sy; y++) {
				for (int64_t x = 1; x < sx; x++) {
					if (labels[x + sx * y] == labels[(x-1) + sx * y]) {
						int64_t node_up = x + sxe * y;
						int64_t node_down = x + sxe * (y + 1);
						adjacency[node_up] |= 0b0100;
						adjacency[node_down] |= 0b1000;
						equivalences.ids[node_down] = node_up;
						all_edges.emplace_back(node_up,node_down);
					}
				}
			}

			// assign horizontal edges
			for (int64_t y = 1; y < sy; y++) {
				for (int64_t x = 0; x < sx; x++) {
					if (labels[x + sx * y] == labels[x + sx * (y-1)]) {
						int64_t node_left = x + sxe * y;
						int64_t node_right = (x+1) + sxe * y;
						adjacency[node_left] |= 0b0001;
						adjacency[node_right] |= 0b0010;
						equivalences.unify(node_left, node_right);
						all_edges.emplace_back(node_left,node_right);
					}
				}
			}
		}
		else {
			// assign vertical edges
			for (int64_t y = 0; y < sy; y++) {
				for (int64_t x = 1; x < sx; x++) {
					if (labels[x + sx * y] != labels[(x-1) + sx * y]) {
						int64_t node_up = x + sxe * y;
						int64_t node_down = x + sxe * (y + 1);
						adjacency[node_up] |= 0b0100;
						adjacency[node_down] |= 0b1000;
						equivalences.ids[node_down] = node_up;
						all_edges.emplace_back(node_up,node_down);
					}
				}
			}

			// assign horizontal edges
			for (int64_t y = 1; y < sy; y++) {
				for (int64_t x = 0; x < sx; x++) {
					if (labels[x + sx * y] != labels[x + sx * (y-1)]) {
						int64_t node_left = x + sxe * y;
						int64_t node_right = (x+1) + sxe * y;
						adjacency[node_left] |= 0b0001;
						adjacency[node_right] |= 0b0010;
						equivalences.unify(node_left, node_right);
						all_edges.emplace_back(node_left, node_right);
					}
				}
			}			
		}

		component_edge_list.resize(all_edges.size() * 2);

		std::vector<int64_t> renumber(sxe*sye);
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

			component_edge_list[membership - 1].push_back(mkedge(pair.first, pair.second));
		}

		component_edge_list.resize(next_label - 1);
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
		for (uint64_t i = 0; i < chain.size(); i++) {
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

	node = pos_x + sxe * pos_y;
}

struct pair_hash {
	size_t operator()(const std::pair<int64_t, int64_t>& p) const {
		return p.first + 31 * p.second;
	}
};


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

  std::vector<std::pair<int64_t, std::vector<char>>> chains;
  std::vector<int64_t> revisit;
  revisit.reserve(sx);
  std::vector<uint8_t> revisit_ct((sx+1)*(sy+1));

  if (n_clusters == 0) {
    return symbols_to_integers(chains);
  }

  std::vector<int64_t> neighbors;
  for (int64_t cluster = 0; cluster < n_clusters; cluster++) {
  	robin_hood::unordered_flat_set<std::pair<int64_t, int64_t>, pair_hash>
  		remaining;

  	auto& cluster_edges = G.component_edge_list[cluster];
  	remaining.insert(cluster_edges.begin(), cluster_edges.end());
  	
  	std::pair<int64_t, int64_t> start_edge = *remaining.begin();

  	int64_t node = start_edge.first;
  	int64_t start_node = start_edge.first;

    std::vector<char> code;
    std::unordered_map<int64_t, std::vector<int64_t>> branch_nodes;
    int64_t branches_taken = 1;

    while (!remaining.empty() || !revisit.empty()) {
    	G.neighbors(node, neighbors);

    	if (!G.adjacency[node]) {
    		code.push_back('t');
    		branches_taken--;
    		if (!revisit.empty()) {
	    		while (!revisit.empty()) {
	    			node = revisit.back();
	    			revisit.pop_back();
	    			if (node > -1) {
	    				revisit_ct[node]--;
	    				break;
	    			}
	    		}
	    	}
    		else if (!remaining.empty()) {
    			node = (*remaining.begin()).first;
    		}
    		continue;
    	}
    	else if (neighbors.size() > 1) {
    		code.push_back('b');
    		revisit.push_back(node);
    		revisit_ct[node]++;
    		branch_nodes[node].push_back(code.size() - 1);
    		branches_taken++;
    	}

    	int64_t next_node = neighbors[0];
    	int64_t dir_taken = next_node - node;

    	if (dir_taken == 1) {
    		code.push_back('r');
    	}
    	else if (dir_taken == -1) {
    		code.push_back('l');
    	}
    	else if (dir_taken > 1) {
    		code.push_back('d');
    	}
    	else {
    		code.push_back('u');
    	}

    	auto edge = mkedge(node, next_node);
    	remaining.erase(edge);
    	G.erase_edge(edge);
    	node = next_node;

    	// if we reencounter a node we've already visited,
    	// remove it from revisit and replace the branch. 
    	// with a skip.
    	if (revisit_ct[node]) {
				int64_t pos = revisit.size() - 1;
				for (; pos >= 0; pos--) {
					if (revisit[pos] == node) {
						break;
					}
				}
				revisit[pos] = -1;
				branches_taken--;
				code[branch_nodes[node].back()] = 's';
				branch_nodes[node].pop_back();
				revisit_ct[node]--;
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


std::vector<uint64_t> read_boc_index(
	const std::vector<unsigned char>& binary,
	const uint64_t sx, const uint64_t sy
) {
	std::vector<uint64_t> nodes;

	const uint64_t sxe = sx + 1;
	const uint64_t sye = sy + 1;

	const uint64_t x_width = crackle::lib::compute_byte_width(sx+1);
	const uint64_t y_width = crackle::lib::compute_byte_width(sy+1);

	uint64_t idx = 0;
	uint32_t index_size = crackle::lib::ctoid(binary, idx, 4);

	uint64_t num_y = crackle::lib::ctoid(binary, idx, y_width);
	idx += y_width;

	uint64_t y = 0; 

	for (uint64_t yi = 0; yi < num_y; yi++) {
		y += crackle::lib::ctoid(binary, idx, y_width);
		idx += y_width;

		uint64_t num_x = crackle::lib::ctoid(binary, idx, x_width);
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
	const std::vector<std::vector<uint64_t>>& chains,
	const uint64_t sx, const uint64_t sy
) {
	const uint64_t sxe = sx + 1;
	const uint64_t sye = sy + 1;

	const uint64_t x_width = crackle::lib::compute_byte_width(sx+1);
	const uint64_t y_width = crackle::lib::compute_byte_width(sy+1);

	// beginning of chain index
	robin_hood::unordered_node_map<uint64_t, std::vector<uint64_t>> boc(sye);

	uint64_t start_y = sye;
	for (auto& chain : chains) {
		uint64_t node = chain[0];
		uint64_t y = node / sxe;
		uint64_t x = node - sxe * y;
		start_y = std::min(start_y, y);
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

std::vector<unsigned char> pack_codes(
	std::vector<std::vector<uint64_t>>& chains,
	const uint64_t sx, const uint64_t sy
) {
	struct {
		bool operator()(
			const std::vector<uint64_t> &a, const std::vector<uint64_t> &b 
		) const { 
			return a[0] < b[0];
		}
	} CmpNodeIndex;

	std::sort(chains.begin(), chains.end(), CmpNodeIndex);

	std::vector<unsigned char> binary = write_boc_index(chains, sx, sy);

	for (auto& chain : chains) {
		uint64_t all_moves = chain.size() - 1;
		all_moves -= (all_moves % 4);

		uint64_t i = 1;
		uint8_t encoded = 0;
		while (i < all_moves) {
			encoded = 0;
			for (uint64_t j = 0; j < 4; j++, i++) {
				encoded |= (chain[i] << (2*j));
			}
			binary.push_back(encoded);
		}
		if (i < chain.size()) {
			encoded = 0;
			for (uint64_t j = 0; i < chain.size(); j++, i++) {
				encoded |= (chain[i] << (2*j));
			}
			binary.push_back(encoded);
		}
	}

	return binary;
}

template <typename LABEL>
std::vector<std::vector<unsigned char>> 
encode_boundaries(
	const LABEL* labels,
	const int64_t sx, const int64_t sy, const int64_t sz,
	const bool permissible
) {
	std::vector<std::vector<unsigned char>> binary_components;

	const int64_t sxy = sx * sy;

	for (int64_t z = 0; z < sz; z++) {

		auto codes = create_crack_codes(labels + z * sxy, sx, sy, permissible);
		binary_components.push_back(
			pack_codes(codes, sx, sy)
		);
	}

	return binary_components;
}

std::unordered_map<uint64_t, std::vector<unsigned char>> 
unpack_binary(
	const std::vector<unsigned char> &code, 
	const uint64_t sx, const uint64_t sy
) {
	std::unordered_map<uint64_t, std::vector<unsigned char>> chains;

	if (code.size() == 0) {
		return chains;
	}

	std::vector<unsigned char> symbols;
	symbols.reserve(code.size() * 4 * 2);

	uint64_t branches_taken = 0;
	uint64_t node = 0;

	char remap[4] = { 'u', 'r', 'l', 'd' };

	std::vector<uint64_t> nodes = read_boc_index(code, sx, sy);
	int node_i = 0;

	for (uint64_t i = 0; i < code.size(); i++) {
		if (branches_taken == 0) {
			node = nodes[node_i];
			node_i++;
			if (node_i >= nodes.size()) {
				throw std::runtime_error("corrupted crack code.");
			}
			branches_taken = 1;
			continue;
		}

		for (uint64_t j = 0; j < 4; j++) {
			uint8_t move = static_cast<uint8_t>((code[i] >> (2*j)) & 0b11);

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
	std::vector<uint8_t> edges(sx * sy, 0b1111);

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