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

		std::vector<uint8_t> read() {
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

		int change_to_base_10() {
			int base_10 = 0;
			for (int i = 0, j = idx; i < length; i++, j++) {
				if (j >= length) {
					j = 0;
				}
				base_10 += pow(4, i) * static_cast<int>(input[j]);
			}
			return base_10;
		}
	}

	void apply_difference_code(
		std::vector<std::vector<uint64_t>>& crack_codes,
	) {
		for (auto code : crack_codes) {
			for (uint64_t i = code.size() - 1; i >= 1; i--) {
				code[i] -= code[i-1];
				if (code[i] > 3) {
					code += 4;
				}
			}
		}
	}

	std::vector<robin_hood::unordered_flat_map<uint8_t,int>> 
	gather_statistics(
		std::vector<std::vector<unsigned char>> &crack_codes,
		int64_t model_order
	) {
		std::vector<robin_hood::unordered_flat_map<uint8_t,int>> stats(
			pow(4, model_order)
		);
		CircularBuf buf(model_order);

		// might be an off-by-one here
		for (auto codes : crack_codes) {
			for (int64_t i = 0; i < codes.size() - model_order; i++) {
				buf.push_back(static_cast<uint8_t>(code[i]));
				int idx = buf.change_to_base_10();
				stats[idx][code[i]]++;
			}
		}

		return stats;
	}

	std::vector<std::vector<uint8_t>> stats_to_model(
		std::vector<robin_hood::unordered_flat_map<int,int>>& stats,
	) {
		struct {
			bool operator()(
				std::pair<int,int>& a, std::pair<int,int>& b
			) const { 
				return a.second >= b.second;
			}
		} CmpIndex;

		std::vector<std::vector<uint8_t>> model(stats.size());
		for (uint64_t i = 0; i < model.size(); i++) {
			std::vector<std::pair<int,int>> pair_row(stats[i].begin(), stats[i].end());
			// most frequent in lowest index
			std::sort(pair_row.begin(), pair_row.end(), CmpIndex);
			std::vector<uint8_t> row(4);
			for (int j = 0; j < 4; j++) {
				row[pair_row[j].first] = j;
			}
			model[i] = std::move(row);
		}

		return model;
	}

	std::vector<uint8_t> decode_markov(
		std::vector<unsigned char>& crack_code,
		std::vector<std::vector<uint8_t>>& model
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

		for (uint64_t i = 0; i < model.size(); i++) {
			stored_model[i] = (
				  (model[i] & 0b11)
				| ((model[i] >> 2) & 0b11)
				| ((model[i] >> 4) & 0b11)
				| ((model[i] >> 6) & 0b11)
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
		std::vector<std::vector<unsigned char>> &crack_codes,
		std::vector<robin_hood::unordered_flat_map<int,int>>& stats,
		int64_t model_order
	) {
		std::vector<unsigned char> bitstream;
		auto model = stats_to_model(stats);

		CircularBuf buf(model_order);

		for (auto code : crack_codes) {
			int pos = 2;
			uint16_t byte = code[0];

			buf.push_back(code[0]);
			for (uint64_t i = 1; i < code.size(); i++) {
				uint8_t idx = model[buf.change_to_base_10()][code[i]];

				if (idx == 0) {
					pos++;
				}
				else if (idx == 1) {
					byte |= (0b10 << pos);
					pos += 2;
				}
				else if (idx == 2) {
					byte |= (0b110 << pos);
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

				buf.push_back(code[i]);
			}
			if (pos > 0) {
				bitstream.push_back(static_cast<uint8_t>(byte));
			}
		}

		return bitstream;
	}
};
};

#endif