/*!
 * \file force_perimeter.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class 
*/ 

#include "force_perimeter.hpp"

namespace VMTutorial
{
  Vec ForcePerimeter::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  {
    Vec l = he.to()->r - v.r;                    // vector along the junction pointing away from the vertex
    Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    Face<Property>& fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    double P1 = _sys.mesh().perim(f);
    double P2 = _sys.mesh().perim(fp);
    double lambda_1, lambda_2;
    
    double  gamma_1 = (f.outer)    ? 0.0 : _gamma[f.data().face_type];
    double  gamma_2 = (fp.outer)   ? 0.0 : _gamma[fp.data().face_type];
    if (!_lambda_P0)
    {
      lambda_1 = (f.outer)   ? 0.0 : _lambda[f.data().face_type];
      lambda_2 = (fp.outer)  ? 0.0 : _lambda[fp.data().face_type];
    }
    else 
    {
      lambda_1 = gamma_1*f.data().P0;
      lambda_2 = gamma_2*fp.data().P0;
    }

    double lambda = lambda_1 + lambda_2;
    double fedges = gamma_1*P1 + gamma_2*P2 - lambda;

    return fedges*l.unit();
  }

  double ForcePerimeter::tension(const HalfEdge<Property>& he)
  {
    Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    Face<Property>& fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    
    double lambda_1, lambda_2;
    
    double gamma_1 = (f.outer)    ? 0.0 : _gamma[f.data().face_type];
    double gamma_2 = (fp.outer)   ? 0.0 : _gamma[fp.data().face_type];
    if (!_lambda_P0)
    {
      lambda_1 = (fh->outer)   ? 0.0 : _lambda[fh->data().face_type];
      lambda_2 = (fh_p->outer) ? 0.0 : _lambda[fh_p->data().face_type];
    }
    else
    {
      lambda_1 = gamma_1*f.data().P0;
      lambda_2 = gamma_2*fp.data().P0;
    }

    double lambda = lambda_1 + lambda_2;
  
    return gamma_1*_sys.mesh().perim(f) + gamma_2*_sys.mesh().perim(fp) - lambda;
  }

  double ForcePerimeter::energy(const Face<Property>& f)
  {
    if (f.outer || f.erased)
      return 0.0;

    double P = _sys.mesh().perim(f);
    
    double lambda = 0.0;
    double gamma = _gamma[f.data().face_type];
    if (!_lambda_P0)
        lambda = _lambda[f.data().face_type];
    }
    
    double P0;
    if (_lambda_P0)
      P0 = f.data().P0;
    else
      P0 = lambda/gamma;
    
    double dP = P - P0;

    return 0.5*gamma*dP*dP;
    
  }

}
