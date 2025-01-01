#ifndef __CRACKLE_DUALGRAPH_HXX__
#define __CRACKLE_DUALGRAPH_HXX__

#include <vector>
#include <cstdint>
#include <algorithm>

#include "crackcodes.hpp"
#include "builtins.hpp"
#include "cc3d.hpp"
#include "libdivide.hpp"
#include "robin_hood.hpp"

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
				// vcg[idx] == 0b11100 means we are in a pinch that has been visited
				// exactly once (it will need to be visited twice). 

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

struct TreeNode {
public:
	TreeNode* parent;
	uint32_t value;
	std::vector<TreeNode*> children;

	TreeNode() : parent(NULL), value(0) {}

	void setParent(TreeNode* parent_) {
		if (parent_ != NULL) {
			parent_->children.push_back(this);
		}
		if (parent != NULL && parent != parent_) {
			parent->children.erase(
				std::remove(parent->children.begin(), parent->children.end(), this), 
				parent->children.end()
			);
		}
		this->parent = parent_;
	}

	bool isRoot() const {
		return parent == NULL;
	}

	const TreeNode* root() const {
		const TreeNode* cur = this;
		while (cur->parent != NULL) {
			cur = cur->parent;
		}
		return cur;
	}

	int depth() const {
		int depth = 0;
		const TreeNode* cur = this;
		while (cur->parent != NULL) {
			cur = cur->parent;
			depth++;
		}
		return depth;
	}

private:
	void _allValuesHelper (std::vector<uint32_t>& all_vals) const {
		all_vals.push_back(value);

		for (TreeNode* node : children) {
			node->_allValuesHelper(all_vals);
		}
	}

public:
	std::vector<uint32_t> allValues () const {
		std::vector<uint32_t> all_vals;
		_allValuesHelper(all_vals);
		return all_vals;
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

void extract_contours_helper(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy,
	std::vector<std::vector<uint32_t>>& contours
) {

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
	// bool is_hole = false; 
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
				uint8_t visit_count = vcg[node] >> 4;
				nodes_already_visited += (visit_count > 0);
				vcg[node] = ((visit_count + 1) << 4) | (vcg[node] & 0b1111);
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
}

bool polygonContainsPoint(
	const std::vector<uint32_t>& poly,
	const std::vector<uint8_t>& vcg,
	const uint32_t pt,
	const uint64_t sx
) {

	const libdivide::divider<uint32_t> fast_sx(sx); 

	uint32_t pt_y = pt / fast_sx;
	uint32_t pt_x = pt - pt_y * sx;
	uint32_t contacts_y = 0;
	uint32_t contacts_x = 0;

	robin_hood::unordered_set<uint32_t> seen;
	seen.reserve(poly.size());

	for (uint64_t i = 0; i < poly.size(); i++) {
		if (seen.count(poly[i])) {
			continue;
		}

		uint32_t elem_y = poly[i] / fast_sx;
		uint32_t elem_x = poly[i] - elem_y * sx;

		// need to check that a vertical chain actually touches
		// a vertical boundary in the vcg
		bool barrier = (vcg[poly[i]] & 0b1100) != 0b1100; // at least one y barrier
		uint32_t incr = (pt_x == elem_x) && (elem_y <= pt_y) && barrier;
		contacts_y += incr;

		barrier = (vcg[poly[i]] & 0b11) != 0b11; // at least one x barrier
		incr = (pt_y == elem_y) && (elem_x <= pt_x) && barrier;
		contacts_x += incr;

		seen.emplace(poly[i]);
	}

	// printf("poly %d pt %d contacts x%d y%d\n", poly[0], pt, contacts_x, contacts_y);

	return (contacts_y & 0b1) && (contacts_x & 0b1);
}


std::vector<std::vector<uint32_t>>
merge_holes(
	std::vector<std::vector<uint32_t>>& candidate_contours,
	const std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy
) {

	const libdivide::divider<uint32_t> fast_sx(sx);
	
	std::vector<
		std::pair<std::pair<uint32_t, uint32_t>, std::pair<uint32_t, uint32_t>>
	> bboxes(candidate_contours.size());
	std::vector<TreeNode> links(candidate_contours.size());

	for (uint64_t i = 0; i < candidate_contours.size(); i++) {
		auto& vec = candidate_contours[i];
		uint32_t minx = sx - 1;
		uint32_t miny = sy - 1;
		uint32_t maxx = 0;
		uint32_t maxy = 0;

		for (auto pt : vec) {
			uint32_t y = pt / fast_sx;
			uint32_t x = pt - y * sx;
			minx = std::min(minx, x);
			maxx = std::max(maxx, x);
			miny = std::min(miny, y);
			maxy = std::max(maxy, y);
		}

		bboxes[i].first = std::make_pair(minx, miny);
		bboxes[i].second = std::make_pair(maxx, maxy);
		links[i].value = i;
	}

	// printf("sz: %d\n", candidate_contours.size());

	for (uint64_t i = 0; i < candidate_contours.size(); i++) {
		if (candidate_contours[i].size() == 1) {
			continue;
		}

		for (uint64_t j = i + 1; j < candidate_contours.size(); j++) {
			auto& bbx1 = bboxes[i];
			auto& bbx2 = bboxes[j];

			// non-intersecting bounding boxes
			if (bbx2.first.second > bbx1.second.second) {
				break;
			}
			else if (!(
				(bbx2.first.first >= bbx1.first.first && bbx2.first.second >= bbx1.first.second)
			 && (bbx2.second.first <= bbx1.second.first && bbx2.second.second <= bbx1.second.second)
			)) {
				continue;
			}

			if (polygonContainsPoint(candidate_contours[i], vcg, candidate_contours[j][0], sx)) {
				// printf("link i %d j %d (idx %d %d)\n", i, j, candidate_contours[i][0], candidate_contours[j][0]);
				links[j].setParent(&links[i]);
			}
		}
	}

	robin_hood::unordered_set<uint32_t> roots;
	roots.reserve(candidate_contours.size() / 10);

	for (uint64_t i = 0; i < candidate_contours.size(); i++) {
		uint32_t depth = links[i].depth();
		uint32_t root = links[i].root()->value;

		// printf("i=%d depth=%d, root=%d, nchild=%d\n", 
		// 	i, depth, links[i].root()->value, links[i].children.size());
		// uint64_t picked = 666;
		
		if (links[i].children.size() == 0) {
			links[i].setParent(NULL);
			roots.emplace(i);
			// picked = i;
		}
		else if (depth == 0) {
			roots.emplace(root);
			// picked = root;
		}
		else if (depth == 1) {
			uint32_t minpt = candidate_contours[i][0];

			// if not blocked on up and left
			// it'll be the same connected component
			// (a hole)
			if (vcg[minpt] & 0b1010) {
				roots.emplace(root);
				// picked = root;
			}
			else {
				links[i].setParent(NULL);
				roots.emplace(i);
				// picked = i;				
			}
		}
		else {
			auto& vec = candidate_contours[root];
			auto it = std::find(vec.begin(), vec.end(), candidate_contours[i][0]);
			if (it != vec.end()) {
				links[i].setParent(&links[root]);
				roots.emplace(root);
				// picked = root;
			}
			else if ((depth & 0b0) == 0) {
				links[i].setParent(NULL);
				roots.emplace(i);
				// picked = i;
			}
			else {
				links[i].parent->setParent(NULL);
				uint32_t root = links[i].root()->value;
				roots.emplace(root);
				// picked = root;
			}
		}

		// printf("i=%d depth=%d, root=%d, picked=%d\n\n", i, depth, links[i].root()->value, picked);
	}

	std::vector<std::vector<uint32_t>> merged_contours(roots.size());

	uint64_t i = 0;
	for (auto root : roots) {
		auto vals = links[root].allValues();
		for (uint32_t val : vals) {

			auto insertion_point_it = merged_contours[i].end();
			if (
				merged_contours[i].size() > 0 
				&& merged_contours[i][0] > candidate_contours[val][0]
			) {
				insertion_point_it = merged_contours[i].begin();
			}

			merged_contours[i].insert(
				insertion_point_it, 
				candidate_contours[val].begin(), 
				candidate_contours[val].end()
			);
		}
		i++;
	}

	return merged_contours;
}


std::vector<std::vector<uint32_t>> 
extract_contours(
	std::vector<uint8_t>& vcg,
	const uint64_t sx, const uint64_t sy 
) {

	std::vector<std::vector<uint32_t>> contours;

	extract_contours_helper(vcg, sx, sy, contours);

	std::sort(contours.begin(), contours.end(),
		[](const auto& a, const auto& b) {
			return a[0] < b[0];
		});

	// int i = 0;
	// for (auto ct : contours) {
	// 	printf("new contour %d\n", i++);
	// 	for (auto val : ct) {
	// 		printf("%d, ", val);
	// 	}		
	// 	printf("\n");
	// }


	std::vector<std::vector<uint32_t>> merged = merge_holes(contours, vcg, sx);

	std::sort(merged.begin(), merged.end(),
		[](const auto& a, const auto& b) {
			return a[0] < b[0];
		});

	// i = 0;
	// for (auto ct : merged) {
	// 	printf("contour %d %d\n", i++, ct[0]);
	// 	printf("\n");
	// }

	return merged;
}

};
};

#endif