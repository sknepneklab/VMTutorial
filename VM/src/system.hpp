/*!
 * \file system.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief System class 
 */ 

#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/detail/xml_parser_writer_settings.hpp>

#include "json.hpp"

#include "rng.hpp"
#include "mesh.hpp"
#include "myostore.hpp"
#include "property.hpp"


using std::string;
using std::istringstream;
using std::istream_iterator;
using std::copy;
using std::back_inserter;
using std::stoi;
using std::stod;
using std::ifstream;
using std::map;
using std::cerr;
using std::runtime_error;
using std::setw;
using std::to_string;
using std::shared_ptr;
using std::make_shared;
using std::ofstream;
using std::sqrt;

namespace pt = boost::property_tree;
using json = nlohmann::json;

namespace AJM
{

  typedef Mesh<Property> MyMesh;
  typedef map<string,int> type_data;
  typedef map<int, string> type_name;
  
  enum SPLIT_TYPE { FORWARD, BACKWARD };

  typedef map<string, double> params_type;         // Used when we set a numerical value to a parameter, e.g., kappa = 1.0
  typedef map<string, Vec> vec_params_type;        // Used when we set a vectorial value to a parameter, e.g., n = Vec(1,0)
  typedef map<string, vector<double>> multi_params_type;   // Used when we need to at least two values to set a parameter, e.g., a parameters is drawn from a random distribution 
  typedef map<string, string> string_params_type;         // Used when we set a string value to a parameter, e.g., force_type = "nematic_vm"
  typedef map<string, RNGState> rng_state;    // Used to store state of various random number generators

  bool operator<(const VertexHandle<Property>&, const VertexHandle<Property>&);

  class System
  {
    public:

      System(MyMesh& mesh) : _mesh{mesh}, 
                             _time_step{0},
                             _simulation_time{0.0},
                             _symmetric_myosin{false},
                             _has_spokes{false},
                             _num_cell_types{0},
                             _num_vert_types{0},
                             _num_junction_types{0},
                             _mesh_set{false},
                             _topology_changed{true},
                             _max_unique_id{0},
                             _has_rng_state{false}
                             { 
                               
                             }

      MyMesh& mesh()  { return _mesh; }

      void set_box(const shared_ptr<Box>& box) { _mesh.set_box(box); }
      const shared_ptr<Box> &box() const { return _mesh.box(); }
      bool periodic() { return (_mesh.box() != nullptr); }

      void read_input(const string&, bool = false);

      void set_simulation_time_step(int time_step) { _time_step = time_step; }
      void increment_max_unique_id() { _max_unique_id++; }
      int& time_step() { return _time_step; }
      double& simulation_time() { return _simulation_time; }
      
      void add_myostore(VertexHandle<Property>& vh, const MyoStore<Property>& ms) { _myostore[vh] = ms; }
      MyoStore<Property>& get_myostore(VertexHandle<Property>& vh) { return _myostore[vh]; }

      type_data& cell_types() { return _cell_types; }
      type_data& vert_types() { return _vert_types; }
      type_data& junction_types() { return _junction_types; }
      const string get_cell_type_name(const int type_id) const { return _cell_types_map.at(type_id); }
      const string get_vert_type_name(const int type_id) const { return _vert_types_map.at(type_id); }
      const string get_junction_type_name(const int type_id) const { return _junction_types_map.at(type_id); }
      int get_num_cell_types() const { return _num_cell_types; }
      int get_num_vert_types() const { return _num_vert_types; }
      int get_max_unique_id() const { return _max_unique_id; }
      void add_cell_type(const string& cell_type)
      {
        if (_cell_types.find(cell_type) == _cell_types.end())
        {
          _cell_types_map[_num_cell_types] = cell_type;
          _cell_types[cell_type] = _num_cell_types++;
        }
      }
      void set_cell_type(int id, const string& type)
      {
        this->add_cell_type(type);
        _mesh.get_face(id).data().face_type = this->_cell_types[type];
        _mesh.get_face(id).data().type_name = type;
      }

      void add_vert_type(const string& vert_type)
      {
        if (_vert_types.find(vert_type) == _vert_types.end())
        {
          _vert_types_map[_num_vert_types] = vert_type;
          _vert_types[vert_type] = _num_vert_types++;
        }
      }
      void add_junction_type(const string& junction_type)
      {
        if (_junction_types.find(junction_type) == _junction_types.end())
        {
          _junction_types_map[_num_junction_types] = junction_type;
          _junction_types[junction_type] = _num_junction_types++;
        }
      }
      void set_junction_type(int id, const string& type)
      {
        this->add_junction_type(type);
        _mesh.get_edge(id).data().edge_type = this->_junction_types[type];
      }

      void displace_vertices(const string &, const Vec&);

      void set_topology_change(bool flag) { _topology_changed = flag; }
      bool topology_change() { return _topology_changed;  }

      void save_rng_state(const string& name, RNGState state)
      {
        _rng_state[name] = state;
        _has_rng_state = true;
      }

      RNGState get_rng_state(const string& name)
      {
        if (_rng_state.find(name) != _rng_state.end())
          return _rng_state.at(name);
        return RNGState{"none", "none", "none"};
      }

      rng_state& get_all_rng_states() { return _rng_state; }

      bool has_rng_state() { return _has_rng_state; }

    private:
      MyMesh &_mesh;
      type_data _cell_types;
      type_name _cell_types_map;
      type_data _vert_types;
      type_name _vert_types_map;
      type_data _junction_types;
      type_name _junction_types_map;
      myostore_map _myostore;
      int _time_step;
      double _simulation_time;
      int _num_cell_types;
      int _num_vert_types; 
      int _num_junction_types; 
      bool _mesh_set;
      bool _topology_changed;  // If true, mesh topology has changed
      int _max_unique_id;
      rng_state _rng_state;  // State of various random number generators
      bool _has_rng_state;   // If true, RNG state is set
  };

  void export_T1_stats(py::module&);
  void export_VertexProperty(py::module&);
  void export_EdgeProperty(py::module&);
  void export_HEProperty(py::module&);
  void export_FaceProperty(py::module&);
  void export_Spoke(py::module &);
  void export_Vertex(py::module &);
  void export_Edge(py::module&);
  void export_HalfEdge(py::module&);
  void export_Face(py::module&);
  void export_Mesh(py::module&);
  void export_System(py::module&);

}

#endif
