#ifndef __CRACKLE_PARTITION_HPP__
#define __CRACKLE_PARTITION_HPP__

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <span>
#include <type_traits>
#include <unordered_map>
#include <mutex>

#include "robin_hood.hpp"

#include "cc3d.hpp"
#include "crc.hpp"
#include "header.hpp"
#include "crackcodes.hpp"
#include "lib.hpp"
#include "labels.hpp"
#include "pins.hpp"
#include "markov.hpp"
#include "dual_graph.hpp"
#include "threadpool.hpp"

namespace crackle {
namespace partition {

struct DecodedFlatLabels {
	uint64_t m_num_labels;
	std::vector<uint64_t> m_unique;
	std::vector<uint64_t> m_components_per_grid;
	std::vector<uint64_t> m_cc_map;

	DecodedFlatLabels(const std::span<unsigned char>& binary) {
		CrackleHeader head(binary);
		if (head.label_format != LabelFormat::FLAT) {
			throw std::runtime_error("Can only decode flat labels.");
		}

		std::span<const unsigned char> labels_binary = crackle::labels::raw_labels(binary);
		m_num_labels = crackle::labels::num_labels(binary);

		uint64_t offset = 8;
		uint64_t uniq_bytes = m_num_labels * head.stored_data_width;
		m_unique = crackle::labels::unique(binary);

		offset += uniq_bytes;

		const uint64_t num_grids = head.num_grids();
		
		m_components_per_grid.resize(num_grids);

		const uint64_t component_width = head.component_width();
		std::span<const unsigned char> components_per_grid_binary(
			labels_binary.begin() + offset, 
			labels_binary.begin() + offset + component_width * num_grids
		);

		for (uint64_t i = 0; i < num_grids; i++) {
			m_components_per_grid[i] = crackle::lib::ctoid(
				components_per_grid_binary, 
				i * component_width, 
				component_width
			);
		}

		offset += components_per_grid_binary.size();

		uint64_t cc_label_width = crackle::lib::compute_byte_width(num_labels);

		uint64_t num_ccs = (labels_binary.size() - offset) / cc_label_width;

		m_cc_map.resize(num_ccs);

		std::span<const unsigned char> cc_map_binary(
			labels_binary.begin() + offset, 
			labels_binary.end()
		);

		for (uint64_t i = 0; i < num_ccs; i++) {
			m_cc_map[i] = crackle::lib::ctoid(cc_map_binary, i * cc_label_width, cc_label_width);
		}
	}
};

std::span<unsigned char> _zstack_flat_labels(
  std::span<uint64_t> uniq,
  std::vector<std::span<unsigned char>> binaries
) {
	robin_hood::unordered_flat_map<uint64_t, uint64_t> uniq_map;
	uniq_map.reserve(uniq.size());

	for (uint64_t i = 0; i < uniq.size(); i++) {
		uniq_map[uniq[i]] = i;
	}

	CrackleHeader first_head(binaries[0]);

	first_head.stored_data_width = crackle::lib::compute_byte_width(uniq.back());
	const uint64_t key_width = crackle::lib::compute_byte_width(uniq.size());

	std::vector<uint64_t> component_index;
	std::vector<uint64_t> all_keys;

	for (auto& binary : binaries) {
		if (binary.size() == 0) {
			continue;
		}

		CrackleHeader head(binary);
		DecodedFlatLabels decoded(binary);

		component_index.insert(
			component_index.begin(), 
			decoded.m_components_per_grid.begin(),
			decoded.m_components_per_grid.end()
		);

		auto& cc_map = decoded.m_cc_map;

		for (uint64_t i = 0; i < cc_map.size(); i++) {
			cc_map[i] = uniq_map[cc_map[i]];
		}

		all_keys.insert(all_keys.end(), cc_map.begin(), cc_map.end());
	}

	std::vector<unsigned char> stacked_label_binary;
	crackle::lib::itocd_insert(stacked_label_binary, uniq, first_head.stored_data_width);
	crackle::lib::itocd_insert(stacked_label_binary, component_index, first_head.component_width());
	crackle::lib::itocd_insert(stacked_label_binary, all_keys, key_width);
	return stacked_label_binary;
}

template <typename Iterator>
std::span<const unsigned char>
zstack(Iterator& begin, const Iterator& end) {
	Iterator it = begin;

	if (it == end) {
		return {};
	}

	std::vector<std::span<unsigned char>> binaries;

	int64_t sz = 0;
	const CrackleHeader first_head(*it);
	int data_width = 1;

	if (first_head.label_format != LabelFormat::FLAT) {
		throw std::runtime_error("Pins are not supported in this implementation.");
	}

	for (Iterator it = begin; it != end; it++) {
		std::span<unsigned char> binary = crackle::reencode_with_markov_order(
			*it, /*markov_model_order=*/0
		);
		CrackleHeader head(binary);
		data_width = std::max(data_width, head.data_width);

		head.fortran_order = first_head.fortran_order;
		std::vector<unsigned char> head_bytes = head.tobytes();
		std::copy(head_bytes.begin(), head_bytes.end(), data.begin());

		if (first_head.sx != head.sx || first_head.sy != head.sy) {
			throw std::runtime_error("All images must have the same width and height.");
		}
		if (head.grid_size != first_head.grid_size) {
			throw std::runtime_error("Grid sizes must match.");
		}
		if (head.crack_format != first_head.crack_format) {
			throw std::runtime_error("All crack formats must match.");
		}
		if (head.signed != first_head.signed) {
			throw std::runtime_error("All binaries must have the same sign.");
		}

		sz += head.sz;
		binaries.push_back(binary);
	}

	if (binaries.size() == 1) {
		return binaries[0];
	}

	first_head.sz = sz;
	first_head.data_width = data_width;

	std::vector<uint64_t> unique;
	for (Iterator it = begin; it != end; it++) {
		auto uniq = crackle::labels::unique(*it);
		unique.insert(unique.end(), uniq.begin(), uniq.end());
	}

	std::sort(unique.begin(), unique.end());
	auto uiter = std::unique(unique.begin(), unique.end());
	unique.erase(uiter, unique.end());

	first_head.stored_data_width = crackle::lib::compute_byte_width(unique.back());

	auto& labels_binary = _zstack_flat_labels(unique, binaries);

	uint64_t z = 0;
	std::vector<uint32_t> zindex(sz);
	std::vector<std::span<unsigned char>> all_crack_codes_binary;

	for (auto& binary : binaries) {
		CrackleHeader head(binary);
		auto& crack_codes = crackle::get_crack_codes(
			head, binary, 0, head.sz
		);
		for (auto& cc : crack_codes) {
			zindex[z] = cc.size();
			all_crack_codes_binary.insert(
				all_crack_codes_binary.end(),
				cc.begin(),
				cc.end()
			);
		}
		z++;
	}

	std::vector<unsigned char> zindex_binary;
	crackle::lib::itocd_insert(zindex_binary, zindex, 4);
	uint32_t grid_crc = crackle::crc::crc32c(zindex_binary);
	crackle::lib::itoc_push_back(grid_crc, zindex_binary);

	std::vector<unsigned char> crack_crcs_binary;

  	if (first_head.format_version > 0) {
		for (auto& binary : binaries) {
			CrackleHeader head(binary);
			auto& crack_code_crcs = crackle::get_crack_code_crcs(
				head, binary
			);
			for (uint32_t crc : crack_code_crcs) {
				crackle::lib::itoc_push_back(crc, crack_crcs_binary);
			}
		}
	}
    
    first_head.is_sorted = true;
    first_head.num_label_bytes = labels_binary.size();

	std::vector<unsigned char> stacked_binary;
	auto& header_bytes = first_head.tobytes();
	stacked_binary.insert(stacked_binary.end(), head_bytes.begin(), head_bytes.end());
	stacked_binary.insert(stacked_binary.end(), zindex_binary.begin(), zindex_binary.end());
	stacked_binary.insert(stacked_binary.end(), labels_binary.begin(), labels_binary.end());
	stacked_binary.insert(stacked_binary.end(), all_crack_codes_binary.begin(), all_crack_codes_binary.end());

	if (first_head.format_version > 0) {
		uint32_t labels_crc = crackle::crc::crc32c(labels_binary);
		crackle::lib::itoc_push_back(labels_crc, stacked_binary);
		stacked_binary.insert(stacked_binary.end(), crack_crcs_binary.begin(), crack_crcs_binary.end());
	}

	return stacked_binary;
}


};
};

#endif