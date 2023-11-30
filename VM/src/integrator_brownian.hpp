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
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief IntegratorBrownian class 
*/

#ifndef __INTEGRATOR_BROWNIAN_HPP__
#define __INTEGRATOR_BROWNIAN_HPP__

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"
#include "constraint_x.hpp"
#include "constraint_y.hpp"
#include "constraint_ux.hpp"
#include "constraint_uy.hpp"
#include "constraint_radial.hpp"

#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>

using namespace std::chrono;
using std::map;
using std::make_unique;

namespace AJM
{

  using ConstrainerType = unique_ptr<Constrainer>;

  class IntegratorBrownian : public Integrator 
  {

    public:

      IntegratorBrownian(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc, seed},
                                                                    _T{0.0},
                                                                    _gamma{1.0}
                                                                    
      { 
        map<string,int>& vert_types = _sys.vert_types();
        map<string,int>& cell_types = _sys.cell_types();
        for (int i = 0; i < vert_types.size(); i++)  _constant_force.push_back(Vec(0.0,0.0));

        _constrainer = make_unique<Constrainer>();
        _constrainer->add<ConstraintNone>("none");
        _constrainer->add<ConstraintFixed>("fixed");
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "T")
            _T = p.second;
          if (p.first == "gamma")
            _gamma = p.second;
        }
      };
      void set_type_params(const string& type, const params_type& params) override { }
      void set_string_params(const string_params_type& params) override { }
      void set_external_force(const string& vtype, const Vec& f) override 
      { 
        _constant_force[_sys.vert_types()[vtype]] = f;
      }
      
      void set_flag(const string& flag) override 
      {  
        
      }

      
    private:

      ConstrainerType _constrainer; // Apply various constraints
      double _T;                 // temperature
      double _gamma;             // friction 
      vector<Vec> _constant_force;
      
    };

}

#endif

