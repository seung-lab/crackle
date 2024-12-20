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
	const py::bytes &buffer,
	int64_t z_start, int64_t z_end
) {
	int64_t voxels = head.sx * head.sy;
	z_start = std::max(z_start, static_cast<int64_t>(0));
	if (z_end == -1) {
		z_end = head.sz;
	}
	z_end = std::min(
		std::max(z_end, static_cast<int64_t>(0)), 
		static_cast<int64_t>(head.sz)
	);

	voxels *= z_end - z_start;

	py::array arr = py::array_t<LABEL>(voxels);
	crackle::decompress<LABEL>(
		buffer,
		reinterpret_cast<LABEL*>(const_cast<void*>(arr.data())),
		z_start, z_end
	);
	return arr;
}

py::array decompress(
	const py::bytes &buffer, 
	const int64_t z_start = 0, const int64_t z_end = -1
) {
	crackle::CrackleHeader head(buffer);

	py::array labels;

	if (head.data_width == 1) {
		labels = decompress_helper<uint8_t>(head, buffer, z_start, z_end);
	}
	else if (head.data_width == 2) {
		labels = decompress_helper<uint16_t>(head, buffer, z_start, z_end);
	}
	else if (head.data_width == 4) {
		labels = decompress_helper<uint32_t>(head, buffer, z_start, z_end);	
	}
	else {
		labels = decompress_helper<uint64_t>(head, buffer, z_start, z_end);
	}
	
	return labels;
}

template <typename LABEL>
py::bytes compress_helper(
	const py::array &labels, 
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0
) {
	const uint64_t sx = labels.shape()[0];
	const uint64_t sy = labels.ndim() < 2
		? 1 
		: labels.shape()[1];
	const uint64_t sz = labels.ndim() < 3 
		? 1 
		: labels.shape()[2];

	std::vector<unsigned char> buf = crackle::compress<LABEL>(
		reinterpret_cast<LABEL*>(const_cast<void*>(labels.data())),
		sx, sy, sz,
		allow_pins,
		fortran_order,
		markov_model_order
	);
	return py::bytes(reinterpret_cast<char*>(buf.data()), buf.size());
}

py::bytes compress(
	const py::array &labels, 
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0
) {
	int width = labels.dtype().itemsize();
	bool is_signed = labels.dtype().kind() == 'i';

	if (is_signed) {
		if (width == 1) {
			return compress_helper<int8_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else if (width == 2) {
			return compress_helper<int16_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else if (width == 4) {
			return compress_helper<int32_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else {
			return compress_helper<int64_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
	}
	else {
		if (width == 1) {
			return compress_helper<uint8_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else if (width == 2) {
			return compress_helper<uint16_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else if (width == 4) {
			return compress_helper<uint32_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
		else {
			return compress_helper<uint64_t>(
				labels, allow_pins, fortran_order, markov_model_order
			);
		}
	}
}

py::bytes reencode_markov(
	const py::bytes &buffer, const int markov_model_order
) {
	auto buf = crackle::reencode_with_markov_order(buffer, markov_model_order);
	return py::bytes(reinterpret_cast<char*>(buf.data()), buf.size());
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

auto point_cloud(	
	const py::bytes &buffer, 
	const int64_t z_start = 0, 
	const int64_t z_end = -1,
	const int64_t label = -1
) {
	return crackle::point_cloud(buffer, z_start, z_end, label);
}


PYBIND11_MODULE(fastcrackle, m) {
	m.doc() = "Accelerated crackle functions."; 
	m.def("decompress", &decompress, "Decompress a crackle file into a numpy array.");
	m.def("compress", &compress, "Compress a numpy array into a binary crackle file returned as bytes.");
	m.def("reencode_markov", &reencode_markov, "Change the markov order of an existing crackle binary.");
	m.def(
		"connected_components", 
		&connected_components,
		"Perform 4-connected components in layers on a 3D array."
	);
	m.def("compute_pins", &compute_pins, "Compute a pinset.");
	m.def("point_cloud", &point_cloud, "Extract one or more point clouds without decompressing.");

	py::class_<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>(m, "CppPin")
	    .def(py::init<>())
	    .def_readwrite("index", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::index)
	    .def_readwrite("depth", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::depth)
	    .def_readwrite("label", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::label);
}
