#ifndef CRACKLE_PACKED_ARRAY_HPP
#define CRACKLE_PACKED_ARRAY_HPP

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace crackle::packed_array {

struct OneBitArray {
	std::vector<uint8_t> m_data;
	uint64_t m_num_entries;

	void resize(uint64_t size) {
		m_num_entries = size;
		m_data.resize(
			(size + 7) >> 3
		);	
	}

	OneBitArray() { 
		m_num_entries = 0;
	}

	OneBitArray(const uint64_t size) {
		resize(size);
	}

	const uint8_t* data() const {
		return m_data.data();
	}

	uint64_t size() const {
		return m_num_entries;
	}

	uint64_t nbytes() const { 
		return m_data.size();
	}

	void set(const uint64_t index, const bool value) {
		if (index >= m_num_entries) {
			throw std::out_of_range("index >= size");
		}

		const uint64_t offset = index >> 3;
		const uint64_t remainder = index & 0b111;
		const uint8_t mask = (uint8_t(1) << remainder);
		std::atomic_ref<uint8_t> ref(m_data[offset]);

		if (value) {
			ref.fetch_or(mask, std::memory_order_relaxed);
		}
		else {
			ref.fetch_and(uint8_t(~mask), std::memory_order_relaxed);
		}
	}

	bool operator[](const uint64_t index) const {
		if (index >= m_num_entries) { 
			throw std::runtime_error("index greater than table size.");
		}
		const uint64_t offset = index >> 3;
		const uint64_t remainder = index & 0b111;
		return (m_data[offset] >> remainder) & 1;
	}

	void fill(const bool value) {
		const uint8_t packed_value = value
			? 0xff
			: 0x00;

		std::fill(m_data.begin(), m_data.end(), packed_value);
	}
};

};
#endif