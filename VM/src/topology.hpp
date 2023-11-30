/*!
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Topology class
 */

#include "system.hpp"
#include "force_compute.hpp"


#include <exception>


using std::exception;

namespace VMTutorial
{
	class Topology
	{
	public:
		Topology(System &sys, ForceCompute &fc, int seed) : _sys{sys},
															_force_compute{fc},
															_min_edge_len{0.02},
															_new_edge_len{0.022},
		{
		}
		~Topology() = default;

		void set_params(const params_type &params)
		{
			for (auto &p : params)
			{
				if (p.first == "min_edge_len")
					_min_edge_len = p.second;
				else if (p.first == "new_edge_len")
					_new_edge_len = p.second;
				else
					throw runtime_error("Unknown topoloy flag.");
			}
		};

		void set_type_params(const string &type, const params_type &params)
		{
		}

		void set_flag(const string &flag)
		{
		}

		void T1();

	private:
		System &_sys;
		ForceCompute &_force_compute;
		double _min_edge_len;
		double _new_edge_len;
	};

	void export_Topology(py::module &);

}
