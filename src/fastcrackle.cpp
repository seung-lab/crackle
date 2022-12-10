#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>

#include "crackle.hpp"
#include "cc3d.hpp"

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

py::tuple connected_components(const py::array &labels) {
	int width = labels.dtype().itemsize();

	const uint64_t sx = labels.shape()[0];
	const uint64_t sy = labels.shape()[1];
	const uint64_t sz = labels.shape()[2];

	std::vector<uint64_t> num_components_per_slice(sz);
	py::array cc_labels = py::array_t<uint32_t>(sx * sy * sz);
	uint64_t N = 0;

	if (width == 1) {
		crackle::cc3d::connected_components<uint8_t, uint32_t>(
			reinterpret_cast<uint8_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, 
			num_components_per_slice,
			reinterpret_cast<uint32_t*>(const_cast<void*>(cc_labels.data())), 
			N 
		);
	}
	else if (width == 2) {
		crackle::cc3d::connected_components<uint16_t, uint32_t>(
			reinterpret_cast<uint16_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, 
			num_components_per_slice,
			reinterpret_cast<uint32_t*>(const_cast<void*>(cc_labels.data())), 
			N 
		);
	}
	else if (width == 4) {
		crackle::cc3d::connected_components<uint32_t, uint32_t>(
			reinterpret_cast<uint32_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, 
			num_components_per_slice,
			reinterpret_cast<uint32_t*>(const_cast<void*>(cc_labels.data())), 
			N 
		);
	}
	else {
		crackle::cc3d::connected_components<uint64_t, uint32_t>(
			reinterpret_cast<uint64_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, 
			num_components_per_slice,
			reinterpret_cast<uint32_t*>(const_cast<void*>(cc_labels.data())), 
			N 
		);
	}

	return py::make_tuple(cc_labels, num_components_per_slice, N);
}


PYBIND11_MODULE(fastcrackle, m) {
	m.doc() = "Accelerated crackle functions."; 
	m.def("decompress", &decompress, "Decompress a crackle file into a numpy array.");
	m.def(
		"connected_components", 
		&connected_components,
		"Perform 4-connected components in layers on a 3D array."
	);
}
