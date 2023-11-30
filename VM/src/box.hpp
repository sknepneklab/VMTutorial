/*!
 * \file box.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Simulation Box
 */

#ifndef __BOX_HPP__
#define __BOX_HPP__

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "matrix.hpp"

using std::shared_ptr;

namespace py = pybind11;

namespace VMTutorial
{

	struct Box
	{
		Box(double lx, double ly) : h{lx, 0.0, 0.0, ly},
									inv_h{h.inv()},
									s_mat{1.0, 0.0, 0.0, 1.0}
		{
		}
		Box(double ax, double ay, double bx, double by) : h{ax, bx, ay, by},
														  inv_h{h.inv()},
														  s_mat{1.0, 0.0, 0.0, 1.0}
		{
		}
		double area() { return h.det(); }
		vector<double> get() { return vector<double>({h._mxx, h._mxy, h._myx, h._myy}); }

		Matrix h;
		Matrix inv_h;
		Matrix s_mat;
	};

	void export_Box(py::module &m);

}

#endif