#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>

#include "crackle.hpp"

namespace py = pybind11;

template <typename LABEL>
py::array decompress_helper(
	const crackle::CrackleHeader& head, 
	const py::bytes &buffer
) {
	py::array arr = py::array_t<LABEL>(head.voxels());
	crackle::decompress<LABEL>(
		buffer, 
		/*fortran_order=*/true,
		reinterpret_cast<LABEL*>(const_cast<void*>(arr.data()))
	);
	return arr;
}

py::array decompress(const py::bytes &buffer) {
	crackle::CrackleHeader head(buffer);

	if (head.data_width == 1) {
		return decompress_helper<uint8_t>(head, buffer);
	}
	else if (head.data_width == 2) {
		return decompress_helper<uint16_t>(head, buffer);
	}
	else if (head.data_width == 4) {
		return decompress_helper<uint32_t>(head, buffer);	
	}
	else {
		return decompress_helper<uint64_t>(head, buffer);
	}
}

PYBIND11_MODULE(fastcrackle, m) {
	m.doc() = "Accelerated crackle functions."; 
	m.def("decompress", &decompress, "Decompress a crackle file into a numpy array.");
}
