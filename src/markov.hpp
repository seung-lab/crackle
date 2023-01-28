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

	std::vector<robin_hood::unordered_flat_map<int,int>> 
	gather_statistics(
		std::vector<std::vector<unsigned char>> &crack_codes,
		int64_t model_order
	) {
		std::vector<robin_hood::unordered_flat_map<int,int>> stats(
			pow(4, model_order)
		);
		CircularBuf buf(model_order);

		// might be an off-by-one here
		for (auto codes : crack_codes) {
			for (int64_t i = 0; i < codes.size() - model_order; i++) {
				buf.push_back(code[i]);
				int idx = buf.change_to_base_10();
				stats[idx][code[i]]++;
			}
		}

		return stats;
	}

	

};
};

#endif