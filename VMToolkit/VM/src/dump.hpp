/*!
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#ifndef __DUMP_HPP__
#define __DUMP_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <regex>

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

#include "json.hpp"

#include "system.hpp"
#include "force_compute.hpp"

using std::back_inserter;
using std::cerr;
using std::copy;
using std::ifstream;
using std::istream_iterator;
using std::istringstream;
using std::map;
using std::move;
using std::ofstream;
using std::runtime_error;
using std::setprecision;
using std::setw;
using std::stod;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::tolower;
using std::transform;
using std::unique_ptr;
using std::vector;

namespace VMTutorial
{

	class Dump
	{
	public:
		Dump(System &sys, ForceCompute &fc) : _sys{sys}, _force_compute{fc}, _sfc(0.95) {}

		void dump_cells(const string &, bool = false, bool = false);
		void dump_junctions(const string &, bool = false);
		void dump_mesh(const string &, bool = false);
		void dump_json(const string &);
		void dump_cell_directors(const string&, bool = false, bool = false);

		void set_sfc(double sfc) { _sfc = sfc; }

	private:
		System &_sys;
		ForceCompute &_force_compute;
		double _sfc; // Scaling factor for junction output
	};

	void to_json(json &, const HalfEdge<Property> &);
	void to_json(json &, const Edge<Property> &);
	void to_json(json &, const Vertex<Property> &);
	void to_json(json &, const Face<Property> &);
	void to_json(json &, const Box &);

	void from_json(const json &, HalfEdge<Property> &);
	void from_json(const json &, Edge<Property> &);
	void from_json(const json &, Vertex<Property> &);
	void from_json(const json &, Face<Property> &);

	vector<string> split(const std::string &, char);

	void export_Dump(py::module &);

}

#endif