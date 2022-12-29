#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <vector>

#include "crackle.hpp"
#include "cc3d.hpp"
#include "pins.hpp"

namespace py = pybind11;

template <typename LABEL>
py::array decompress_helper(
	const crackle::CrackleHeader& head, 
	const py::bytes &buffer
) {
	py::array arr = py::array_t<LABEL>(head.voxels());
	crackle::decompress<LABEL>(
		buffer,
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

template <typename LABEL>
py::bytes compress_helper(
	const py::array &labels, 
	const bool force_flat = false,
	const bool fortran_order = true
) {
	const uint64_t sx = labels.shape()[0];
	const uint64_t sy = labels.shape()[1];
	const uint64_t sz = labels.shape()[2];

	std::vector<unsigned char> buf = crackle::compress<LABEL>(
		reinterpret_cast<LABEL*>(const_cast<void*>(labels.data())),
		sx, sy, sz,
		force_flat,
		fortran_order
	);
	return py::bytes(reinterpret_cast<char*>(buf.data()), buf.size());
}

py::bytes compress(
	const py::array &labels, 
	const bool force_flat = false,
	const bool fortran_order = true
) {
	int width = labels.dtype().itemsize();

	if (width == 1) {
		return compress_helper<uint8_t>(labels, force_flat, fortran_order);
	}
	else if (width == 2) {
		return compress_helper<uint16_t>(labels, force_flat, fortran_order);
	}
	else if (width == 4) {
		return compress_helper<uint32_t>(labels, force_flat, fortran_order);
	}
	else {
		return compress_helper<uint64_t>(labels, force_flat, fortran_order);
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

auto compute_pins(const py::array &labels) {
	int width = labels.dtype().itemsize();

	const uint64_t sx = labels.shape()[0];
	const uint64_t sy = labels.shape()[1];
	const uint64_t sz = labels.shape()[2];

	py::dict results;

	if (width == 1) {
		return crackle::pins::compute<uint8_t>(
			reinterpret_cast<uint8_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz
		);
	}
	else if (width == 2) {
		return crackle::pins::compute<uint16_t>(
			reinterpret_cast<uint16_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz
		);
	}
	else if (width == 4) {
		return crackle::pins::compute<uint32_t>(
			reinterpret_cast<uint32_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz
		);
	}
	else {
		return crackle::pins::compute<uint64_t>(
			reinterpret_cast<uint64_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz
		);		
	}
}

PYBIND11_MODULE(fastcrackle, m) {
	m.doc() = "Accelerated crackle functions."; 
	m.def("decompress", &decompress, "Decompress a crackle file into a numpy array.");
	m.def("compress", &compress, "Compress a numpy array into a binary crackle file returned as bytes.");
	m.def(
		"connected_components", 
		&connected_components,
		"Perform 4-connected components in layers on a 3D array."
	);
	m.def("compute_pins", &compute_pins, "Compute a pinset.");

	py::class_<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>(m, "CppPin")
	    .def(py::init<>())
	    .def_readwrite("index", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::index)
	    .def_readwrite("depth", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::depth)
	    .def_readwrite("label", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::label);
}
