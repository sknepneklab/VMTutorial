/*!
 * \file force_compute.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceCompute class
 */

#include "force_compute.hpp"

namespace VMTutorial
{
	void export_ForceCompute(py::module &m)
	{
		py::class_<ForceCompute>(m, "Force")
			.def(py::init<System &>())
			.def("set_params", &ForceCompute::set_params)
			.def("set_vec_params", &ForceCompute::set_vec_params)
			.def("set_flag", &ForceCompute::set_flag)
			.def("add", &ForceCompute::add_force)
			.def("compute", &ForceCompute::compute_forces)
			.def("energy", &ForceCompute::total_energy);
	}
}