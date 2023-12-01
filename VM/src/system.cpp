/*!
 * \file system.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief System class 
 */ 

#include "system.hpp"

namespace VMTutorial
{
  
  // Read input from a JSON file
  void System::read_input(const string& json_file, bool read_params)
  {
    if (_mesh_set)
    {
      cout << "Warning! Mesh has already been set. Overwriting it." << endl;
      _mesh.wipe();
    }
    ifstream inp(json_file.c_str());
    json j;
    inp >> j;
    inp.close();
    // Check if simulation box exists
    if (j["mesh"].find("box") != j["mesh"].end())
      if (j["mesh"]["box"]["periodic"])
      {
        if (j["mesh"]["box"].find("a") != j["mesh"]["box"].end() && j["mesh"]["box"].find("b") != j["mesh"]["box"].end())
        {
          vector<double> a = j["mesh"]["box"]["a"];
          vector<double> b = j["mesh"]["box"]["b"];
          this->set_box(make_shared<Box>(a[0], a[1], b[0], b[1]));
        }
        else if (j["mesh"]["box"].find("lx") != j["mesh"]["box"].end() && j["mesh"]["box"].find("ly") != j["mesh"]["box"].end())
        {
          double lx = j["mesh"]["box"]["lx"];
          double ly = j["mesh"]["box"]["ly"];
          this->set_box(make_shared<Box>(lx, ly));
        }
        else 
          throw runtime_error("There is a problem with the box field in the input JSON file.");
        cout << "Setting periodic simulation box." << endl;
      }

    // Check for time step
    if (j["mesh"].find("time_step") != j["mesh"].end())
      _time_step = j["mesh"]["time_step"];

    // Populate vertices
    for (int i = 0; i < j["mesh"]["vertices"].size(); i++)
    {
      int id = j["mesh"]["vertices"][i]["id"];
      double x = j["mesh"]["vertices"][i]["r"][0];
      double y = j["mesh"]["vertices"][i]["r"][1];
      bool boundary = j["mesh"]["vertices"][i]["boundary"];
      _mesh.add_vertex(Vertex<Property>(id, Vec(x, y, _mesh.box()), boundary, _mesh));
      Vertex<Property> &v = _mesh.vertices().back();
      this->add_vert_type(j["mesh"]["vertices"][i]["type"]);
      v.erased = j["mesh"]["vertices"][i]["erased"];
      v.data().vert_type = _vert_types[j["mesh"]["vertices"][i]["type"]];
      v.data().type_name = get_vert_type_name(v.data().vert_type);
      v.data().constraint = j["mesh"]["vertices"][i]["constraint"];
      if (j["mesh"]["vertices"][i].find("velocity") != j["mesh"]["vertices"][i].end())
      {
        v.data().vel.x = j["mesh"]["vertices"][i]["velocity"][0];
        v.data().vel.y = j["mesh"]["vertices"][i]["velocity"][1];
      }
    }
    cout << "Finished reading vertices." << endl;
    // Populate faces
    bool erased_face = false;
    for (int i = 0; i < j["mesh"]["faces"].size(); i++)
    {
      if (j["mesh"]["faces"][i].find("erased") != j["mesh"]["faces"][i].end())
      {
        erased_face = j["mesh"]["faces"][i]["erased"];
        if (erased_face)
          _mesh.add_face(vector<int>(),  true);
        else
          _mesh.add_face(j["mesh"]["faces"][i]["vertices"], false);
      }
      else 
        _mesh.add_face(j["mesh"]["faces"][i]["vertices"]);
      Face<Property>& f = _mesh.faces().back();
      f.data().unique_id = f.id;
      this->increment_max_unique_id();
      if (!erased_face)
      {
        this->add_cell_type(j["mesh"]["faces"][i]["type"]);
        f.data().face_type = _cell_types[j["mesh"]["faces"][i]["type"]];
        f.data().type_name = get_cell_type_name(f.data().face_type);
        f.outer = j["mesh"]["faces"][i]["outer"];
        if (j["mesh"]["faces"][i].find("A0") != j["mesh"]["faces"][i].end())
          f.data().A0 = j["mesh"]["faces"][i]["A0"];
        if (j["mesh"]["faces"][i].find("P0") != j["mesh"]["faces"][i].end())
          f.data().P0 = j["mesh"]["faces"][i]["P0"];
      }
    }
    cout << "Finished reading faces." << endl;
    _mesh.tidyup();
    cout << "Finished mesh setup." << endl;
    cout << "Mesh has " << _mesh.vertices().size() << " vertices " << _mesh.edges().size() << " edges and " << _mesh.faces().size() << " faces." << endl;
    
    _mesh_set = true;
    
    cout << "Finished reading input configuration." << endl;
  }

  

  // Displace vertices of a given type by a given amount
  void System::displace_vertices(const string & vtype, const Vec& dr)
  {
    if (_vert_types.find(vtype) == _vert_types.end())
      throw runtime_error("Unknown vertex type " + vtype + " in displace vertices.");
    int vert_type = _vert_types[vtype];
    for (VertexHandle<Property> vh = _mesh.vertices().begin(); vh != _mesh.vertices().end(); vh++)
      if (vh->data().vert_type == vert_type)
        _mesh.move_vertex(vh->id, dr);
  }

  // Used to be able to make a map of MyoStore
  bool operator<(const VertexHandle<Property>& lv, const VertexHandle<Property>& rv) 
  {
    return (lv->id < rv->id);
  }

  // Python exports

  
  void export_VertexProperty(py::module& m)
  {
    py::class_<Property::VertexProperty>(m, "VertexProperty")
      .def_readonly("type", &Property::VertexProperty::vert_type)
      .def_readonly("type_name", &Property::VertexProperty::type_name)
      .def_readonly("vel", &Property::VertexProperty::vel)
      .def_readonly("force", &Property::VertexProperty::force);
  }

  void export_EdgeProperty(py::module& m)
  {
    py::class_<Property::EdgeProperty>(m, "EdgeProperty")
      .def_readwrite("tension", &Property::EdgeProperty::tension)
      .def_readwrite("l0", &Property::EdgeProperty::l0)
      .def_readwrite("type", &Property::EdgeProperty::edge_type);
  }

  void export_HEProperty(py::module& m)
  {
    py::class_<Property::HEProperty>(m, "HEProperty")
      .def_readonly("tension", &Property::HEProperty::tension)
      .def_readonly("l0", &Property::HEProperty::l0)
      .def_readonly("force_type", &Property::HEProperty::force_type);
  }

  
  void export_FaceProperty(py::module& m)
  {
    py::class_<Property::FaceProperty>(m, "CellProperty")
      .def_readonly("type", &Property::FaceProperty::face_type)
      .def_readonly("type_name", &Property::FaceProperty::type_name)
      .def_readonly("unique_id", &Property::FaceProperty::unique_id)
      .def_readwrite("A0", &Property::FaceProperty::A0)
      .def_readwrite("P0", &Property::FaceProperty::P0)
      .def_readwrite("cell_type", &Property::FaceProperty::face_type)
      .def_readwrite("n", &Property::FaceProperty::n);
  }


  void export_Vertex(py::module& m)
  {
    py::class_<Vertex<Property>>(m, "Vertex")
      .def(py::init<Mesh<Property>&>())
      .def_readwrite("r", &Vertex<Property>::r)
      .def_readonly("id", &Vertex<Property>::id)
      .def_readonly("erased", &Vertex<Property>::erased)
      .def_readwrite("boundary", &Vertex<Property>::boundary)
      .def_readonly("coordination", &Vertex<Property>::coordination)
      .def("he", [](Vertex<Property>& v) { return *(v.he()); })
      .def("force", [](Vertex<Property>& v) { return v.data().force; })
      .def("vel", [](Vertex<Property>& v) { return v.data().vel; })
      .def("property", (Property::VertexProperty& (Vertex<Property>::*)()) &Vertex<Property>::data, py::return_value_policy::reference);
  }

  void export_Edge(py::module& m)
  {
    py::class_<Edge<Property>>(m, "Edge")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("i", &Edge<Property>::i)
      .def_readonly("j", &Edge<Property>::j)
      .def_readonly("erased", &Edge<Property>::erased)
      .def_readonly("boundary", &Edge<Property>::boundary)
      .def_property_readonly("id", &Edge<Property>::idx)
      .def("he", [](Edge<Property>& e) { return *(e.he()); })
      .def("type", [](Edge<Property>& e) { return e.data().edge_type; })
      .def("property", (Property::EdgeProperty& (Edge<Property>::*)()) &Edge<Property>::data, py::return_value_policy::reference);
  }

  void export_HalfEdge(py::module& m)
  {
    py::class_<HalfEdge<Property>>(m, "Junction")
      .def(py::init<Mesh<Property>&>())
      .def("vfrom", [](HalfEdge<Property>& he) { return *(he.from()); })
      .def("vto", [](HalfEdge<Property>& he) { return *(he.to()); })
      .def("edge", [](HalfEdge<Property>& he) { return *(he.edge()); })
      .def("id", &HalfEdge<Property>::idx)
      .def("next", [](HalfEdge<Property>& he) { return *(he.next()); })
      .def("prev", [](HalfEdge<Property>& he) { return *(he.prev()); })
      .def("pair", [](HalfEdge<Property>& he) { return *(he.pair()); })
      .def("face", [](HalfEdge<Property>& he) { return *(he.face()); })
      .def("property", (Property::HEProperty& (HalfEdge<Property>::*)()) &HalfEdge<Property>::data, py::return_value_policy::reference);
  }

  void export_Face(py::module& m)
  {
    py::class_<Face<Property>>(m, "Cell")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("id", &Face<Property>::id)
      .def_readonly("neighbours", &Face<Property>::nsides)
      .def_readonly("outer", &Face<Property>::outer)
      .def_readonly("erased", &Face<Property>::erased)
      .def("type", [](Face<Property>& f) { return f.data().face_type; })
      .def("A0", [](Face<Property>& f) { return f.data().A0; })
      .def("P0", [](Face<Property>& f) { return f.data().P0; })
      .def("he", [](Face<Property>& f) { return *(f.he()); })
      .def("property", (Property::FaceProperty& (Face<Property>::*)()) &Face<Property>::data, py::return_value_policy::reference);
  }

  void export_Mesh(py::module& m)
  {
    py::class_<Mesh<Property>>(m, "Tissue")
      .def(py::init<>())
      .def("num_vert", &Mesh<Property>::num_vert)
      .def("num_cells", &Mesh<Property>::num_faces)
      .def("tidyup", &Mesh<Property>::tidyup)
      .def("set_cell_type", [](Mesh<Property>& m, int i, int type) { m.get_face(i).data().face_type = type; })
      .def("set_junction_type", [](Mesh<Property>& m, int i, int type) { m.get_edge(i).data().edge_type = type; })
      .def("set_cell_A0", [](Mesh<Property>& m, int i, double A0) { m.get_face(i).data().A0 = A0;  })
      .def("set_cell_P0", [](Mesh<Property>& m, int i, double P0) { m.get_face(i).data().P0 = P0;  })
      .def("get_vertex", &Mesh<Property>::get_vertex , py::return_value_policy::reference)
      .def("get_junction", &Mesh<Property>::get_halfedge, py::return_value_policy::reference)
      .def("get_cell", &Mesh<Property>::get_face, py::return_value_policy::reference)
      .def("vertices", &Mesh<Property>::vertices, py::return_value_policy::reference)
      .def("junctions", &Mesh<Property>::edges, py::return_value_policy::reference)
      .def("halfedges", &Mesh<Property>::halfedges, py::return_value_policy::reference)
      .def("cells", &Mesh<Property>::faces, py::return_value_policy::reference)
      .def("get_cell_centre", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centre(*(m.get_mesh_face(i))); })
      .def("get_cell_centroid", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centroid(*(m.get_mesh_face(i))); })
      .def("get_centre", &Mesh<Property>::get_centre)
      .def("cell_area", [](Mesh<Property> &m, int i) -> double { return m.area(*(m.get_mesh_face(i))); })
      .def("cell_perim", [](Mesh<Property> &m, int i) -> double { return m.perim(*(m.get_mesh_face(i))); });
  }

  void export_System(py::module& m)
  {
    py::class_<System>(m, "System")
      .def(py::init<MyMesh&>())
      .def("read_input", &System::read_input, py::arg("input_file"), py::arg("read_params") = false)
      .def("add_cell_type", &System::add_cell_type)
      .def("set_cell_type", &System::set_cell_type)
      .def("mesh", &System::mesh)
      .def("cell_types", &System::cell_types)
      .def("time_step", &System::time_step)
      .def("simulation_time", &System::simulation_time)
      .def("set_simulation_time_step", &System::set_simulation_time_step)
      .def("get_cell_type_name", &System::get_cell_type_name)
      .def("get_vert_type_name", &System::get_vert_type_name)
      .def("displace_vertices", &System::displace_vertices);
  }

  
}

