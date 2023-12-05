/*!
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 30-Nov-2023
 * \brief Topology class 
*/

#include "topology.hpp"

namespace VMTutorial
{
  
  void Topology::T1()
  {
    for (auto e : _sys.mesh().edges())
    {
      if (!(e.he()->from()->data().constraint == "fixed" || e.he()->to()->data().constraint == "fixed"))
        if (!e.erased && _sys.mesh().len(e) < _min_edge_len)
          if (_sys.mesh().T1(e, _new_edge_len))
            _sys.set_topology_change(true);
    }
  }

  
  void export_Topology(py::module& m)
  {
    py::class_<Topology>(m, "Topology")
      .def(py::init<System&, ForceCompute&, int>(), py::arg("sys"), py::arg("force_compute"), py::arg("seed") = 0)
      .def("set_params", &Topology::set_params)
      .def("set_type_params", &Topology::set_type_params)
      .def("set_flag", &Topology::set_flag);
  }
  
}
