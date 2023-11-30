/* ***************************************************************************
 *
 *  Copyright (C) 2017 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of AJM (Active Junction Model) program.
 *
 *  AJM is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  AJM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file integrator.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 20-May-2019
 * \brief IntegratorBrownian class 
*/

#include "integrator_brownian.hpp"

namespace AJM
{
  void IntegratorBrownian::step()
  {
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    // Compute force on each vertex
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);


    // If we are computing stress, add the term due to friction
    double SXX = 0.0, SYY = 0.0, SXY = 0.0, SYX = 0.0;
    if (_force_compute.stress_compute())
    {
      int Nc = 0;
      for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
      {
        if (!(fh->erased || fh->outer)) 
        {
          _force_compute.compute_stress(fh);
          if (_add_friction_stress)
          {
            HEHandle<Property> he = fh->he();
            HEHandle<Property> first = fh->he();
            Vec rc = _sys.mesh().get_face_centre(fh);
            vector<double> stress(4,0.0);
            do
            {
              double z = he->from()->coordination;
              Vec ri = he->from()->r - rc;
              Vec vi = mu*he->from()->data().force;
              stress[0] -= ri.x*vi.x/z;
              stress[1] -= 0.5*(ri.x*vi.y + ri.y*vi.x)/z;
              stress[2] -= 0.5*(ri.x*vi.y + ri.y*vi.x)/z;
              stress[3] -= ri.y*vi.y/z;
              he = he->next();
            } while (he != first);
            double fact = _gamma/_sys.mesh().area(fh);
            transform(stress.begin(),stress.end(),stress.begin(),[fact](double v) -> double {return fact*v; });
            transform(fh->data().stress.begin(),fh->data().stress.end(), stress.begin(), fh->data().stress.begin(), std::plus<double>());
            transform(fh->data().stress_v.begin(),fh->data().stress_v.end(), stress.begin(), fh->data().stress_v.begin(), std::plus<double>());
          }
          Nc++;
          SXX += fh->data().stress[0];
          SXY += fh->data().stress[1];
          SYX += fh->data().stress[2];
          SYY += fh->data().stress[3];
        }
      }
      SXX /= Nc;
      SXY /= Nc;
      SYX /= Nc;
      SYY /= Nc;
    }
    
    // This is actual integrator 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        // add external force 
        Vec f = vh->data().force + _constant_force[vh->data().vert_type];
        vh->data().f_type["fpull"] =  _constant_force[vh->data().vert_type];

        if (_has_radial_force)
          f = f + _radial_force_magnitude[vh->data().vert_type] * vh->r.unit();

        // apply constraint
        if (_constraint_enabled)
          f = _constrainer->apply_vertex(vh, f);

        Vec rold = vh->r;
        vh->r += _dt*mu*f;  // deterministic part of the integrator step
        if (_T > 0.0)
        {
          Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
          Vec fr = ffr;
          if (_constraint_enabled) 
            fr = _constrainer->apply_vector(vh, ffr);
          vh->r += sqrt_dt*fr;  // update vertex position due to noise
        }
        vh->data().vel = (1.0 / _dt) * (vh->r - rold);  
      } 
    }
    
    if (_apply_pressure)
    {
      double eta_xx = 1.0 - _beta*_dt*(_Pxx - SXX);
      double eta_xy = _beta*_dt*SXY;
      double eta_yx = _beta*_dt*SYX;
      double eta_yy = 1.0 - _beta*_dt*(_Pyy - SYY);
      Matrix inv_h0 = _sys.mesh().box()->inv_h;
      _sys.mesh().box()->transform(eta_xx, eta_xy, eta_yx, eta_yy, false, true);
      Matrix h = _sys.mesh().box()->h;
      Matrix h_inv_h0 = h * inv_h0;
      for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
        vh->r = h_inv_h0*vh->r;
    }

  }

}