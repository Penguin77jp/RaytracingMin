#include <pybind11/pybind11.h>

int hoge() {
	return 1;
}

PYBIND11_MODULE(cppmodule, m) {
	m.doc() = "pybind11";
	m.def("hoge", &hoge, "hoge function");
}