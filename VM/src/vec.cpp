/*!
 * \file vec.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief 2D vectors 
 */ 


#include "vec.hpp"

namespace VMTutorial
{
  void export_Vec(py::module& m)
  {
    py::class_<Vec>(m, "Vec")
          .def(py::init<double,double>())
          .def(py::init<double,double,const shared_ptr<Box>&>())
          .def_readwrite("x", &Vec::x)
          .def_readwrite("y", &Vec::y)
          .def("box", [](const Vec& v) { if (v.box) return v.box->get(); else return vector<double>({0,0,0,0}); })
          .def("__add__", [](const Vec& v1, const Vec& v2) { return Vec(v1 + v2); })
          .def("__sub__", [](const Vec& v1, const Vec& v2) { return Vec(v1 - v2); })
          .def("__mul__", [](const Vec& v, double c) { return c*v; })
          .def("__rmul__", [](const Vec& v, double c) { return c*v; })
          .def("__repr__", [](const Vec& v) { return "("+std::to_string(v.x)+","+std::to_string(v.y)+")";  })
          .def("unit", [](const Vec& v) { return Vec(v.x/v.len(), v.y/v.len(), v.box); });

  }
}