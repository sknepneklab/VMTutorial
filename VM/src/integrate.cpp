/*!
 * \file integrate.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-May-2019
 * \brief Integrate class 
 */ 

#include "integrate.hpp"

namespace VMTutorial
{
  void export_Integrate(py::module& m)
  {
    py::class_<Integrate>(m, "Integrate")
      .def(py::init<System&, ForceCompute&, int>())
      .def("set_params", &Integrate::set_params)
      .def("set_type_params", &Integrate::set_type_params)
      .def("set_string_params", &Integrate::set_string_params)
      .def("set_external_force", &Integrate::set_external_force)
      .def("set_radial_force", &Integrate::set_radial_force)
      .def("set_flag", &Integrate::set_flag)
      .def("set_affine_vel", &Integrate::set_affine_vel)
      .def("enable", &Integrate::enable)
      .def("disable", &Integrate::disable)
      .def("set_dt", &Integrate::set_dt)
      .def("enable_constraint", &Integrate::enable_constraint)
      .def("converged", &Integrate::converged)
      .def("add", &Integrate::add_integrator);
  }
}