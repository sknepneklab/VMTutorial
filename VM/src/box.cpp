/*!
 * \file box.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Simulation Box
 */

#include "box.hpp"

namespace VMTutorial
{
	void export_Box(py::module &m)
	{
		py::class_<Box, shared_ptr<Box>>(m, "Box")
			.def(py::init<double, double>())
			.def(py::init<double, double, double, double>())
			.def("a", [](const Box &b)
				 { return vector<double>({b.h._mxx, b.h._myx}); })
			.def("b", [](const Box &b)
				 { return vector<double>({b.h._mxy, b.h._myy}); });
	}
}