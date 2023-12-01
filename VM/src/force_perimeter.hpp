/*!
 * \file force_perimeter.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class
 */

#ifndef __FORCE_PERIMETER_HPP__
#define __FORCE_PERIMETER_HPP__

#include "force.hpp"

namespace VMTutorial
{

	// Force on a vertex
	class ForcePerimeter : public Force
	{
	public:
		ForcePerimeter(System &sys) : Force{sys},
									  _lambda_P0{false}
		{
			_gamma.resize(_sys.cell_types().size(), 0.0);
			_lambda.resize(_sys.cell_types().size(), 0.0);
		}
		virtual ~ForcePerimeter() {}

		// computes force on vertex by a given edge
		Vec compute(const Vertex<Property> &, const HalfEdge<Property> &) override;

		double tension(const HalfEdge<Property> &) override;

		// Energy calculation
		double energy(const Face<Property> &) override;

		// set all parameters for a given type
		void set_params(const string &cell_type, const params_type &params) override
		{
			for (auto p : params)
				if (p.first != "gamma" && p.first != "lambda")
					throw runtime_error("Unknown parameter " + p.first + ".");

			if (params.find("gamma") == params.end())
				throw runtime_error("Perimeter force requires parameter gamma.");
			if (params.find("lambda") == params.end())
				throw runtime_error("Perimeter force requires parameter lambda.");

			try
			{
				if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
					throw runtime_error("Fore perimeter: Cell type " + cell_type + " is not defined.");
				if (_gamma.size() < _sys.get_num_cell_types())
					_gamma.resize(_sys.get_num_cell_types(), 0.0);
				if (_lambda.size() < _sys.get_num_cell_types())
					_lambda.resize(_sys.get_num_cell_types(), 0.0);
				int ct = _sys.cell_types()[cell_type];
				_gamma[ct] = params.at("gamma");
				_lambda[ct] = params.at("lambda");
			}
			catch (const exception &e)
			{
				cerr << "Problem with setting perimeter force parameters. Exception: " << e.what() << '\n';
				throw;
			}
		}

		// set all vector-valued parameters for a given type
		void set_vec_params(const string &cell_type, const vec_params_type &params) override{};

		void set_flag(const string &flag) override
		{
			if (flag == "use_P0")
				_lambda_P0 = true;
			else
				throw runtime_error("Unknown flag : " + flag + ".");
		}

		
	private:
		vector<double> _gamma;
		vector<double> _lambda;
		bool _lambda_P0; // If true, lambda will be computed as gamma*P0 where P0 is read from the input configuration
	};

}

#endif
