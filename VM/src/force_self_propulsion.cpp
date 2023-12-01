/*!
 * \file force_self_propulsion.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceSelfPropulsion class 
*/ 

#include "force_self_propulsion.hpp"

namespace VMTutorial
{
  Vec ForceSelfPropulsion::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  {
    
    const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    
    double v0 = (f.outer)   ? 0.0 : _v0[f.data().face_type];
    Vec  n  = (f.outer)   ? Vec(0,0) : _n[f.data().face_type];
    
    return (v0/v.coordination)*n;
  }

}
