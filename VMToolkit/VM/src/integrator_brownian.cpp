/*!
 * \file integrator.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 30-Nov-2023
 * \brief IntegratorBrownian class 
*/

#include "integrator_brownian.hpp"

namespace VMTutorial
{
  void IntegratorBrownian::step()
  {
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    // Compute force on each vertex
    for (auto& v : _sys.mesh().vertices())
      if (!v.erased)
        _force_compute.compute(v);
    
    // This is actual integrator 
    for (auto& v : _sys.mesh().vertices())
    {
      if (!v.erased)
      {
        // add external force 
        Vec f = v.data().force + _constant_force[v.data().vert_type];
        
        // apply constraint
        if (_constraint_enabled)
          f = _constrainer->apply_vertex(v, f);

        Vec rold = v.r;
        v.r += _dt*mu*f;  // deterministic part of the integrator step
        if (_T > 0.0)
        {
          Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
          Vec fr = ffr;
          if (_constraint_enabled) 
            fr = _constrainer->apply_vector(v, ffr);
          v.r += sqrt_dt*fr;  // update vertex position due to noise
        }
        v.data().vel = (1.0 / _dt) * (v.r - rold);  
      }
    }

    // Update direction of cell director using simple Brownian dynamics
    if (_update_n)
    {
      double stoch_coeff = sqrt(_Dr*_dt);
      for (auto& f : _sys.mesh().faces())
      {
        if (!(f.erased || f.outer))
        {
          double dtheta = stoch_coeff*_rng.gauss_rng();
          Vec n = f.data().n;
          // Sine and cosine of the rotation angle
          double c = cos(dtheta), s = sin(dtheta);
          // Rotation matrix around z axis
          double Rxx = c, Rxy = -s;
          double Ryx = s, Ryy = c;
          // Apply rotation matrix
          double nx = Rxx*n.x + Rxy*n.y;
          double ny = Ryx*n.x + Ryy*n.y;
          double len = sqrt(nx*nx + ny*ny);
          f.data().n = Vec(nx/len, ny/len);  // Rotation does not change length of a vector, but numerical errors can accumulate, so we normalize it
        }
      }
    }

  }

}