/*!
 * \file force_area.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceArea class 
 */ 

#ifndef __FORCE_AREA_HPP__
#define __FORCE_AREA_HPP__

#include "force.hpp"


namespace VMTutorial
{

  // Force on a vertex
  class ForceArea : public Force
  {
    public:
                                                                                                                                
      ForceArea(System& sys) : Force(sys) 
      { 
        _kappa.resize(_sys.cell_types().size(), 0.0);
      }
      virtual ~ForceArea() { }
        
      // computes force on vertex by a given edge
      Vec compute(const Vertex<Property>&, const HalfEdge<Property>&) override;  
      
      // Are term does not generate tension
      double tension(const HalfEdge<Property>& he) override
      {
        return 0.0;
      }


      // Energy calculation 
      double energy(const Face<Property>&) override;
      
      // set all parameters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "kappa")
            throw runtime_error("Unknown parameter "+p.first+".");
            
        if (params.find("kappa") == params.end())
          throw runtime_error("Area force requires parameter kappa.");

        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore area: Cell type " + cell_type + " is not defined.");
          if (_kappa.size() < _sys.get_num_cell_types())
            _kappa.resize(_sys.get_num_cell_types(), 0.0);
          int ct = _sys.cell_types()[cell_type];
          _kappa[ct] = params.at("kappa");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting area force parameters. Exception: " << e.what() << '\n';
          throw;
        }
      }

      
      // set all vector-valued parameters for a given type
      void set_vec_params(const string& cell_type, const vec_params_type& params) override { };


      void set_flag(const string& flag) override   {   }

    
    private:

      vector<double> _kappa; 
      
      
  };

  
}

#endif
