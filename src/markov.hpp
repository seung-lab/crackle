#ifndef __MARKOV_HXX__
#define __MARKOV_HXX__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <type_traits>

#include "robin_hood.hpp"

namespace crackle {
namespace markov {

	struct CircularBuf {
		uint8_t* data;
		int length;
		int idx;

		CircularBuf(int model_order) {
			data = new uint8_t[model_order]();
			length = model_order;
			idx = 0;
		}

		~CircularBuf() {
			delete [] data;
		}

		void push_back(uint8_t elem) {
			data[idx] = elem;
			idx++;
			if (idx >= length) {
				idx = 0;
			}
		}

		std::vector<uint8_t> read() const {
			std::vector<uint8_t> out;
			out.reserve(length);
			for (int i = 0, j = idx; i < length; i++, j++) {
				if (j >= length) {
					j = 0;
				}
				out.push_back(data[j]);
			}
			return out;
		}

		int change_to_base_10() const {
			int base_10 = 0;
			for (int i = 0, j = idx; i < length; i++, j++) {
				if (j >= length) {
					j = 0;
				}
				base_10 += pow(4, i) * static_cast<int>(data[j]);
			}
			return base_10;
		}
	};

	std::tuple<std::vector<uint64_t>, std::vector<uint8_t>> difference_codepoints(
		robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>& chains
	) {
		std::vector<uint64_t> nodes;
		for (auto& [node, code] : chains) {
			nodes.push_back(node);
		}
		std::sort(nodes.begin(), nodes.end());

		std::vector<uint8_t> codepoints;
		for (uint64_t node : nodes) {
			auto chain = chains[node];
			for (uint8_t codepoint : chain) {
				codepoints.push_back(codepoint);
			}
		}
		if (codepoints.size() > 0) {
			for (uint64_t i = codepoints.size() - 1; i >= 1; i--) {
				codepoints[i] -= codepoints[i-1];
				if (codepoints[i] > 3) {
					codepoints[i] += 4;
				}
			}
		}
		return std::make_tuple(nodes, codepoints);
	}

	std::vector<robin_hood::unordered_flat_map<uint8_t,int>> 
	gather_statistics(
		const std::vector<robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>> &crack_codes,
		const uint64_t model_order
	) {
		std::vector<robin_hood::unordered_flat_map<uint8_t,int>> stats(
			pow(4, model_order)
		);

		for (auto slice : crack_codes) {
			for (auto& [node, code] : slice) {
				CircularBuf buf(model_order);
				for (uint64_t i = 0; i < code.size() - model_order; i++) {
					buf.push_back(static_cast<uint8_t>(code[i]));
					int idx = buf.change_to_base_10();
					stats[idx][code[i]]++;
				}
			}
		}

		return stats;
	}

	std::vector<std::vector<uint8_t>> stats_to_model(
		std::vector<robin_hood::unordered_flat_map<uint8_t,int>>& stats
	) {
		struct {
			bool operator()(
				robin_hood::pair<uint8_t,int>& a, robin_hood::pair<uint8_t,int>& b
			) const { 
				return a.second >= b.second;
			}
		} CmpIndex;

		// model is: index is direction, value is which 
		// codepoint to use

		std::vector<std::vector<uint8_t>> model(stats.size());
		for (uint64_t i = 0; i < model.size(); i++) {
			std::vector<robin_hood::pair<uint8_t,int>> pair_row;
			pair_row.reserve(4);
			for (auto& pair : stats[i]) {
				pair_row.push_back(pair);
			}
			// most frequent in lowest index
			std::sort(pair_row.begin(), pair_row.end(), CmpIndex);
			std::vector<uint8_t> row(4);
			std::vector<bool> marked(4);
			uint64_t j = 0;
			for (j = 0; j < pair_row.size(); j++) {
				row[pair_row[j].first] = j;
				marked[pair_row[j].first] = true;
			}
			// handle sparse statistics
			if (j < 4) {
				for (uint64_t k = 0; k < 4; k++) {
					if (marked[k]) {
						continue;
					}
					row[k] = j;
					j++;
				}
			}
			model[i] = std::move(row);
		}

		return model;
	}

	std::vector<uint8_t> decode_codepoints(
		std::vector<unsigned char>& crack_code,
		const std::vector<std::vector<uint8_t>>& model
	) {
		std::vector<uint8_t> data_stream;

		int model_order = static_cast<int>(log2(model.size()) / log2(4));
		CircularBuf buf(model_order);

		int pos = 2;
		uint8_t start_dir = crack_code[0] & 0b11;
		data_stream.push_back(start_dir);
		buf.push_back(start_dir);

		for (uint64_t i = 0; i < crack_code.size(); i++) {
			uint16_t byte = crack_code[i]; 
			if (i < crack_code.size() - 1) {
				byte |= (crack_code[i+1] << 8);
			}

			while (pos < 8) {
				uint8_t codepoint = (byte >> pos) & 0b111;
				uint64_t model_row = buf.change_to_base_10();

				if ((codepoint & 0b1) == 0) {
					data_stream.push_back(model[model_row][0]);
					pos++;
				}
				else if ((codepoint & 0b10) == 0) {
					data_stream.push_back(model[model_row][1]);
					pos += 2;
				}
				else if ((codepoint & 0b100) == 0) {
					data_stream.push_back(model[model_row][2]);	
					pos += 3;
				}
				else {
					data_stream.push_back(model[model_row][3]);
					pos += 3;
				}

				buf.push_back(data_stream.back());
			}

			pos -= 8;
		}

		for (uint64_t i = 1; i < data_stream.size(); i++) {
			data_stream[i] += data_stream[i-1];
			if (data_stream[i] > 3) {
				data_stream[i] -= 4;
			}
		}

		return data_stream;
	}

	std::vector<unsigned char> to_stored_model(
		const std::vector<std::vector<uint8_t>>& model
	) {
		std::vector<unsigned char> stored_model(model.size());

		// invert keys and values for model to make decoding faster.
		// assumption: reading occurs more often than writing

		struct {
			bool operator()(
				robin_hood::pair<uint8_t,uint8_t>& a, robin_hood::pair<uint8_t,uint8_t>& b
			) const { 
				return a.second < b.second;
			}
		} CmpValue;

		for (uint64_t i = 0; i < model.size(); i++) {
			std::vector<robin_hood::pair<uint8_t, uint8_t>> decode_row;
			decode_row.reserve(4);
			for (int j = 0; j < 4; j++) {
				decode_row.emplace_back(j, model[i][j]);
			}
			std::sort(decode_row.begin(), decode_row.end(), CmpValue);

			stored_model[i] = (
				  (decode_row[0].first & 0b11)
				| ((decode_row[1].first & 0b11) << 2)
				| ((decode_row[2].first & 0b11) << 4)
				| ((decode_row[3].first & 0b11) << 6)
			);
		}

		return stored_model;
	}

	std::vector<std::vector<uint8_t>> from_stored_model(
		const std::vector<unsigned char>& model_stream
	) {
		std::vector<std::vector<uint8_t>> model(model_stream.size());

		for (uint64_t i = 0; i < model_stream.size(); i++) {
			model[i].resize(4);
			model[i][0] = model_stream[i] & 0b11;
			model[i][1] = (model_stream[i] >> 2) & 0b11;
			model[i][2] = (model_stream[i] >> 4) & 0b11;
			model[i][3] = (model_stream[i] >> 6) & 0b11;
		}

		return model;
	}

	std::vector<unsigned char> encode_markov(
		const std::vector<uint8_t> &codepoints,
		const std::vector<std::vector<uint8_t>>& model,
		int64_t model_order
	) {
		std::vector<unsigned char> bitstream;

		if (codepoints.size() == 0) {
			return bitstream;
		}

		CircularBuf buf(model_order);

		int pos = 2;
		uint16_t byte = codepoints[0];

		buf.push_back(codepoints[0]);
		for (uint64_t i = 1; i < codepoints.size(); i++) {
			uint8_t idx = model[buf.change_to_base_10()][codepoints[i]];

			if (idx == 0) {
				pos++;
			}
			else if (idx == 1) {
				byte |= (0b01 << pos);
				pos += 2;
			}
			else if (idx == 2) {
				byte |= (0b011 << pos);
				pos += 3;
			}
			else {
				byte |= (0b111 << pos);
				pos += 3;
			}

			if (pos >= 8) {
				bitstream.push_back(static_cast<uint8_t>(byte));
				pos -= 8;
				byte >>= 8;
			}

			buf.push_back(codepoints[i]);
		}
		if (pos > 0) {
			bitstream.push_back(static_cast<uint8_t>(byte));
		}

		return bitstream;
	}

	std::vector<unsigned char>
	compress(
		robin_hood::unordered_node_map<uint64_t, std::vector<uint8_t>>& chains,
		const std::vector<std::vector<uint8_t>>& model,
		const int64_t model_order,
		const uint64_t sx, const uint64_t sy
	) {
		std::vector<uint64_t> nodes;
		for (auto& [node, code] : chains) {
			nodes.push_back(node);
		}
		std::sort(nodes.begin(), nodes.end());

		std::vector<uint8_t> codepoints;
		for (uint64_t node : nodes) {
			auto chain = chains[node];
			for (uint8_t codepoint : chain) {
				codepoints.push_back(codepoint);
			}
		}
		if (codepoints.size() > 0) {
			for (uint64_t i = codepoints.size() - 1; i >= 1; i--) {
				codepoints[i] -= codepoints[i-1];
				if (codepoints[i] > 3) {
					codepoints[i] += 4;
				}
			}
		}

		std::vector<unsigned char> binary = crackle::crackcodes::write_boc_index(nodes, sx, sy);
		std::vector<unsigned char> bitstream = encode_markov(
			codepoints, model, model_order
		);
		binary.insert(binary.end(), bitstream.begin(), bitstream.end());
		return binary;
	}
};
};

#endif
