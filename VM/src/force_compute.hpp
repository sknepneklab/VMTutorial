/*!
 * \file force_compute.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceCompute class 
 */ 

#ifndef __FORCE_COMPUTE_HPP__
#define __FORCE_COMPUTE_HPP__

#include <exception>
#include <algorithm>
#include <map>
#include <string>

#include "system.hpp"
#include "class_factory.hpp"
#include "force.hpp"
#include "force_area.hpp"
#include "force_perimeter.hpp"
#include "force_self_propulsion.hpp"

 
using std::runtime_error;
using std::transform;
using std::map;
using std::string;


namespace VMTutorial
{

  class ForceCompute : public ClassFactory<Force>
  {
    public:

      ForceCompute(System& sys) : _sys{sys}
      { 

      }

      ~ForceCompute() = default; 

      ForceCompute(const ForceCompute&) = delete;

      void compute_forces()
      {
        for (auto v : _sys.mesh().vertices())
          if (!v.erased)
            this->compute(v);
      }
      
      void compute(Vertex<Property> &v)
      {
        v.data().force = Vec(0.0,0.0);
        for (auto he : v.circulator())
          v.data().force += this->compute(v, he);
      }

      Vec compute(Vertex<Property> &v, const HalfEdge<Property> &he)
      {
        Vec ftot(0,0);
        for (auto& f : this->factory_map)
          ftot += f.second->compute(v, he);
        return ftot;
      }

      double tension(HalfEdge<Property>& he)
      {
        double T = 0.0;
        for (auto& f : this->factory_map)
          T += f.second->tension(he);
        return T;
      }

      double energy(const Face<Property>& face)
      {
        double E = 0.0;
        for (auto& f : this->factory_map)
          E += f.second->energy(face);
        return E;
      }

      double total_energy()
      {
        double E = 0.0;
        for (auto f : _sys.mesh().faces())  
          E += this->energy(f);
        return E;
      }

      void set_params(const string& fname, const string& type, const params_type& params)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_params(type, params);
        else
          throw runtime_error("set_params: Force type " + fname + " is not used in this simulation.");
      }

      void set_vec_params(const string& fname, const string& type, const vec_params_type& params)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_vec_params(type, params);
        else
          throw runtime_error("set_vec_params: Force type " + fname + " is not used in this simulation.");
      }

      void set_flag(const string& fname, const string& flag)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_flag(flag);
        else
          throw runtime_error("set_flag: Force type " + fname + " is not used in this simulation.");
      }

      void add_force(const string& fname)
      {
        string name = fname; 
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (name == "area")
          this->add<ForceArea,System&>(name, _sys);
        else if (name == "perimeter")
          this->add<ForcePerimeter,System&>(name, _sys);
        else if (name == "self-propulsion")
          this->add<ForceSelfPropulsion,System&>(name, _sys);
        else 
          throw runtime_error("Unknown force type : " + name + ".");
      }

    private: 

      System& _sys;
 
  };

  void export_ForceCompute(py::module&);

}
#endif
