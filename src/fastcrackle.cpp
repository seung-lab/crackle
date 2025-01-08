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
	const uint8_t* buffer,
	const uint64_t num_bytes,
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
		buffer, num_bytes,
		reinterpret_cast<LABEL*>(const_cast<void*>(arr.data())),
		z_start, z_end
	);
	return arr;
}

py::array decompress(
	const py::buffer buffer, 
	const int64_t z_start = 0, const int64_t z_end = -1
) {
	py::buffer_info info = buffer.request();

	if (info.ndim != 1) {
		throw std::runtime_error("Expected a 1D buffer");
	}

	uint8_t* data = static_cast<uint8_t*>(info.ptr);

	crackle::CrackleHeader head(data);

	py::array labels;

	if (head.data_width == 1) {
		labels = decompress_helper<uint8_t>(
			head, data, info.size, z_start, z_end
		);
	}
	else if (head.data_width == 2) {
		labels = decompress_helper<uint16_t>(
			head, data, info.size, z_start, z_end
		);
	}
	else if (head.data_width == 4) {
		labels = decompress_helper<uint32_t>(
			head, data, info.size, z_start, z_end
		);	
	}
	else {
		labels = decompress_helper<uint64_t>(
			head, data, info.size, z_start, z_end
		);
	}
	
	return labels;
}

template <typename LABEL>
py::bytes compress_helper(
	const py::array &labels, 
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0,
	const bool optimize_pins = false,
	const bool auto_bgcolor = true,
	const int64_t manual_bgcolor = 0
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
		markov_model_order,
		optimize_pins,
		auto_bgcolor,
		manual_bgcolor
	);
	return py::bytes(reinterpret_cast<char*>(buf.data()), buf.size());
}

py::bytes compress(
	const py::array &labels, 
	const bool allow_pins = false,
	const bool fortran_order = true,
	const uint64_t markov_model_order = 0,
	const bool optimize_pins = false,
	const bool auto_bgcolor = true,
	const int64_t manual_bgcolor = 0
) {
	int width = labels.dtype().itemsize();
	bool is_signed = labels.dtype().kind() == 'i';

#define CALL_COMPRESS_HELPER(LABELS_T) compress_helper<LABELS_T>(\
		labels, allow_pins, fortran_order, markov_model_order,\
		optimize_pins, auto_bgcolor, manual_bgcolor\
	)

	if (is_signed) {
		if (width == 1) {
			return CALL_COMPRESS_HELPER(int8_t);
		}
		else if (width == 2) {
			return CALL_COMPRESS_HELPER(int16_t);
		}
		else if (width == 4) {
			return CALL_COMPRESS_HELPER(int32_t);
		}
		else {
			return CALL_COMPRESS_HELPER(int64_t);
		}
	}
	else {
		if (width == 1) {
			return CALL_COMPRESS_HELPER(uint8_t);
		}
		else if (width == 2) {
			return CALL_COMPRESS_HELPER(uint16_t);
		}
		else if (width == 4) {
			return CALL_COMPRESS_HELPER(uint32_t);
		}
		else {
			return CALL_COMPRESS_HELPER(uint64_t);
		}
	}
#undef CALL_COMPRESS_HELPER
}

py::bytes reencode_markov(
	const py::buffer buffer, const int markov_model_order
) {
	py::buffer_info info = buffer.request();

	if (info.ndim != 1) {
		throw std::runtime_error("Expected a 1D buffer");
	}

	uint8_t* data = static_cast<uint8_t*>(info.ptr);
	auto buf = crackle::reencode_with_markov_order(data, info.size, markov_model_order);
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
			sx, sy, sz, false
		);
	}
	else if (width == 2) {
		return crackle::pins::compute<uint16_t>(
			reinterpret_cast<uint16_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, false
		);
	}
	else if (width == 4) {
		return crackle::pins::compute<uint32_t>(
			reinterpret_cast<uint32_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, false
		);
	}
	else {
		return crackle::pins::compute<uint64_t>(
			reinterpret_cast<uint64_t*>(const_cast<void*>(labels.data())),
			sx, sy, sz, false
		);
	}
}

py::dict point_cloud(	
	const py::buffer buffer, 
	const int64_t z_start = 0, 
	const int64_t z_end = -1,
	const int64_t label = -1
) {
	py::buffer_info info = buffer.request();

	if (info.ndim != 1) {
		throw std::runtime_error("Expected a 1D buffer");
	}

	uint8_t* data = static_cast<uint8_t*>(info.ptr);
	auto ptc = crackle::point_cloud(data, info.size, z_start, z_end, label);

	py::dict py_ptc;
	for (const auto& [key, vec] : ptc) {
		py::array array = py::array_t<uint16_t>(vec.size());
		std::memcpy(array.mutable_data(), vec.data(), vec.size() * sizeof(uint16_t));
		py_ptc[py::int_(key)] = array;
	}

	return py_ptc;
}

py::array get_slice_vcg(
	const py::bytes &buffer, 
	const int64_t z
) {
	std::vector<uint8_t> vcg = crackle::decode_slice_vcg(buffer, z);
	py::array array = py::array_t<uint8_t>(vcg.size());
	std::memcpy(array.mutable_data(), vcg.data(), vcg.size() * sizeof(uint8_t));
	return array;
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
	m.def("get_slice_vcg", &get_slice_vcg, "Debugging tool for examining the voxel connectivity graph of a slice.");

	py::class_<crackle::pins::Pin<uint64_t, uint64_t, uint64_t>>(m, "CppPin")
		.def(py::init<>())
		.def_readwrite("index", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::index)
		.def_readwrite("depth", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::depth)
		.def_readwrite("label", &crackle::pins::Pin<uint64_t, uint64_t, uint64_t>::label);
}
