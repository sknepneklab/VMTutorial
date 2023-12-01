/*!
 * \file dump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#include "dump.hpp"

namespace VMTutorial
{
	void Dump::dump_cells(const string &vtk_file, bool binary_output, bool draw_periodic)
	{
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkDoubleArray> forces = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> face_ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> vert_type = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkDoubleArray> areas = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> perims = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkIntArray> cell_types = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkDoubleArray> num_neigh = vtkSmartPointer<vtkDoubleArray>::New();
		
		ids->SetName("Id");
		ids->SetNumberOfComponents(1);
		areas->SetName("Area");
		areas->SetNumberOfComponents(1);
		perims->SetName("Perimeter");
		perims->SetNumberOfComponents(1);
		face_ids->SetName("FaceId");
		face_ids->SetNumberOfComponents(1);
		vert_type->SetName("VertType");
		vert_type->SetNumberOfComponents(1);
		forces->SetName("Force");
		forces->SetNumberOfComponents(3);
		num_neigh->SetName("NumNeigh");
		num_neigh->SetNumberOfComponents(1);
		cell_types->SetName("CellTypes");
		cell_types->SetNumberOfComponents(1);

		int id = 0;
		map<int, int> id_map;
		for (auto v : _sys.mesh().vertices())
		{
			if (!v.erased)
			{
				points->InsertNextPoint(v.r.x, v.r.y, 0.0);
				ids->InsertNextValue(id);
				vert_type->InsertNextValue(v.data().vert_type);
				double f[3] = {v.data().force.x, v.data().force.y, 0.0};
				forces->InsertNextTuple(f);
				id_map[v.id] = id;
				id++;
			}
		}

		for (auto f : _sys.mesh().faces())
		{
			if (!f.erased)
			{
				double alpha = 0;
				bool omit = false;
				bool crosses_pbc = false;
				vector<int> face_verts;
				if (_sys.periodic())
				{
					Vec r0 = _sys.mesh().get_face_centroid(f);
					Vec s0 = r0.box->inv_h * r0;
					for (auto he : f.circulator())
					{
						Vec s = r0.box->inv_h * he.from()->r;
						double sx = s.x - s0.x, sy = s.y - s0.y;
						if (fabs(rint(sx)) > 1e-6 || fabs(rint(sy)) > 1e-6)
						{
							if (!draw_periodic)
							{
								omit = true;
								break;
							}
							else
							{
								crosses_pbc = true;
								Vec r = r0.box->h * Vec(s.x - rint(sx), s.y - rint(sy));
								points->InsertNextPoint(r.x, r.y, 0.0);
								ids->InsertNextValue(id);
								vert_type->InsertNextValue(he.from()->data().vert_type);
								double f[3] = {0.0, 0.0, 0.0};
								forces->InsertNextTuple(f);
								face_verts.push_back(id);
								id++;
							}
						}
						else
							face_verts.push_back(id_map[he.from()->id]);
					} 
				}
				if (!(f.outer || omit))
				{
					face->GetPointIds()->SetNumberOfIds(_sys.mesh().face_sides(f));
					int i = 0;
					for (auto he : f.circulator())
					{
						double le = (he.from()->r - he.to()->r).len();
						if (!crosses_pbc)
							face->GetPointIds()->SetId(i++, id_map[he.from()->id]);
						else
						{
							face->GetPointIds()->SetId(i, face_verts[i]);
							i++;
						}
					} 
					faces->InsertNextCell(face);
					areas->InsertNextValue(_sys.mesh().area(f));
					perims->InsertNextValue(_sys.mesh().perim(f));
					face_ids->InsertNextValue(f.id);
					cell_types->InsertNextValue(f.data().face_type);
					num_neigh->InsertNextValue(f.nsides);
				}
			}
		}

		polydata->SetPoints(points);
		polydata->GetPointData()->AddArray(ids);
		polydata->GetPointData()->AddArray(vert_type);
		polydata->GetPointData()->AddArray(forces);


		polydata->SetPolys(faces);
		polydata->GetCellData()->AddArray(areas);
		polydata->GetCellData()->AddArray(perims);
		polydata->GetCellData()->AddArray(face_ids);
		polydata->GetCellData()->AddArray(cell_types);
		polydata->GetCellData()->AddArray(num_neigh);

		// Write the file
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetFileName(vtk_file.c_str());
#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(polydata);
#else
		writer->SetInputData(polydata);
#endif
		if (binary_output)
		{
			writer->SetDataModeToBinary();
			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
			compressor->SetCompressionLevel(9); // max compression level
			writer->SetCompressor(compressor);
		}
		else
			writer->SetDataModeToAscii();
		writer->Write();
	}

	void Dump::dump_junctions(const string &vtk_file, bool binary_output)
	{
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

		vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();
		vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> he_id = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> edge_id = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkDoubleArray> lens = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> tension = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> l0 = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkIntArray> true_ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> edge_type = vtkSmartPointer<vtkIntArray>::New();
				

		ids->SetName("Id");
		ids->SetNumberOfComponents(1);
		he_id->SetName("HE_Id");
		he_id->SetNumberOfComponents(1);
		edge_id->SetName("Edge_Id");
		edge_id->SetNumberOfComponents(1);
		lens->SetName("Length");
		lens->SetNumberOfComponents(1);
		tension->SetName("Tension");
		tension->SetNumberOfComponents(1);
		l0->SetName("l0");
		l0->SetNumberOfComponents(1);
		true_ids->SetName("TrueId");
		true_ids->SetNumberOfComponents(1);
		edge_type->SetName("Type");
		edge_type->SetNumberOfComponents(1);
		

		int id = 0;
		map<int, int> id_map;
		map<int, bool> omit;

		// build list of points belonging to faces, with shift (for display purpose)
		for (auto f : _sys.mesh().faces())
		{
			if (!f.erased)
			{
				Vec fc = _sys.mesh().get_face_centre(f); // face centre coordinates

				for (auto he : f.circulator())
				{
					Edge<Property>& e = *(he.edge());
					if (!e.erased)
					{
						Vertex<Property>& v1 = *(he.from());
						Vertex<Property>& v2 = *(he.to());
						Vec r1_shift = fc + _sfc * (v1.r - fc);
						Vec r2_shift = fc + _sfc * (v2.r - fc);
						if (f.outer)
						{
							r1_shift = v1.r;
							r2_shift = v2.r;
						}
						omit[he.idx()] = false;
						if (_sys.periodic())
						{
							Vec s1 = r1_shift.box->inv_h * r1_shift;
							Vec s2 = r2_shift.box->inv_h * r2_shift;
							double sx = s1.x - s2.x, sy = s1.y - s2.y;
							if (fabs(rint(sx)) != 0.0 || fabs(rint(sy)) != 0.0)
								omit[he.idx()] = true;
						}
						if (!(omit[he.idx()]))
						{
							points->InsertNextPoint(r1_shift.x, r1_shift.y, 0.0);
							ids->InsertNextValue(id);
							true_ids->InsertNextValue(v1.id);
							id_map[v1.id] = id;
							id++;

							points->InsertNextPoint(r2_shift.x, r2_shift.y, 0.0);
							ids->InsertNextValue(id);
							true_ids->InsertNextValue(v2.id);
							id_map[v2.id] = id;
							id++;
						}
					}
				} 
			}
		}

		polydata->SetPoints(points);
		polydata->GetPointData()->AddArray(ids);
		polydata->GetPointData()->AddArray(true_ids);

		// create connectivity
		id = 0;
		for (auto f :  _sys.mesh().faces())
		{
			if (!f.erased)
			{
				for (auto he : f.circulator())
				{
					Edge<Property>& e = *(he.edge());
					if (!(e.erased || omit[he.idx()]))
					{
						edge->GetPointIds()->SetId(0, id);
						id++;
						edge->GetPointIds()->SetId(1, id);
						id++;
						he_id->InsertNextValue(he.idx());
						edge_id->InsertNextValue(e.idx());
						lines->InsertNextCell(edge);
						lens->InsertNextValue(_sys.mesh().len(e));
						double he_tension = he.data().tension;
						if (fabs(he_tension) < 1e-10)
							he_tension = _force_compute.tension(he);
						tension->InsertNextValue(he_tension);
						l0->InsertNextValue(e.data().l0);
						edge_type->InsertNextValue(e.data().edge_type);
					}
				} 
			}
		}

		polydata->SetLines(lines);
		polydata->GetCellData()->AddArray(lens);
		polydata->GetCellData()->AddArray(he_id);
		polydata->GetCellData()->AddArray(edge_id);
		polydata->GetCellData()->AddArray(tension);
		polydata->GetCellData()->AddArray(l0);
		polydata->GetCellData()->AddArray(edge_type);

		// Write the file
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetFileName(vtk_file.c_str());
#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(polydata);
#else
		writer->SetInputData(polydata);
#endif
		if (binary_output)
		{
			writer->SetDataModeToBinary();
			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
			compressor->SetCompressionLevel(9); // max compression level
			writer->SetCompressor(compressor);
		}
		else
			writer->SetDataModeToAscii();
		writer->Write();
	}

	
	void Dump::dump_mesh(const string &mesh_file, bool copy_params)
	{
		vector<string> strs = split(mesh_file, '.');
		string ext = strs[strs.size() - 1];
		transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c)
				  { return tolower(c); });
		pt::ptree out;
		pt::ptree rng;
		pt::ptree mesh;
		pt::ptree vertices;
		pt::ptree faces;
		if (_sys.periodic())
		{
			pt::ptree box;
			pt::ptree lx, ly;
			pt::ptree A, B;
			pt::ptree periodic;
			pt::ptree mxx, mxy, myx, myy;
			mxx.put("", _sys.mesh().box()->h._mxx);
			myx.put("", _sys.mesh().box()->h._myx);
			A.push_back(std::make_pair("", mxx));
			A.push_back(std::make_pair("", myx));
			mxy.put("", _sys.mesh().box()->h._mxy);
			myy.put("", _sys.mesh().box()->h._myy);
			B.push_back(std::make_pair("", mxy));
			B.push_back(std::make_pair("", myy));
			lx.put("", _sys.mesh().box()->h._mxx);
			ly.put("", _sys.mesh().box()->h._myy);
			if (fabs(_sys.mesh().box()->h._mxy) >= 1e-6 || fabs(_sys.mesh().box()->h._myx) >= 1e-6)
			{
				box.add_child("a", A);
				box.add_child("b", B);
			}
			else
			{
				box.add_child("lx", lx);
				box.add_child("ly", ly);
			}
			periodic.put("", "true");
			box.add_child("periodic", periodic);
			mesh.add_child("box", box);
		}
		for (auto v : _sys.mesh().vertices())
		{
			pt::ptree vertex;
			pt::ptree id;
			pt::ptree boundary;
			pt::ptree x;
			pt::ptree y;
			pt::ptree type;
			pt::ptree r;
			pt::ptree fx;
			pt::ptree fy;
			pt::ptree force;
			pt::ptree vx;
			pt::ptree vy;
			pt::ptree velocity;
			pt::ptree constraint;
			pt::ptree erased;
			id.put("", v.id);
			boundary.put("", v.boundary);
			type.put("", v.data().type_name);
			x.put("", v.r.x);
			y.put("", v.r.y);
			r.push_back(std::make_pair("", x));
			r.push_back(std::make_pair("", y));
			force.push_back(std::make_pair("", fx));
			force.push_back(std::make_pair("", fy));
			vx.put("", v.data().vel.x);
			vy.put("", v.data().vel.y);
			velocity.push_back(std::make_pair("", vx));
			velocity.push_back(std::make_pair("", vy));
			constraint.put("", v.data().constraint);
			erased.put("", v.erased);
			vertex.add_child("id", id);
			vertex.add_child("boundary", boundary);
			vertex.add_child("type", type);
			vertex.add_child("r", r);
			vertex.add_child("force", force);
			vertex.add_child("velocity", velocity);
			vertex.add_child("constraint", constraint);
			vertex.add_child("erased", erased);
			vertices.push_back(std::make_pair("", vertex));
		}
		mesh.add_child("vertices", vertices);
		int face_id = 0;
		for (auto f : _sys.mesh().faces())
		{
			vector<int> verts;
			vector<double> tension;
			if (!f.erased)
			{
				for (auto he : f.circulator())
				{
					verts.push_back(he.from()->id);
					tension.push_back(he.data().tension);
				} 
			}
			
			pt::ptree face;
			pt::ptree id;
			pt::ptree original_id;
			pt::ptree outer;
			pt::ptree nsides;
			pt::ptree type;
			pt::ptree erased;
			pt::ptree A0;
			pt::ptree P0;
			pt::ptree press;
			pt::ptree fverts;
			pt::ptree fneighs;
			pt::ptree rc;
			pt::ptree rc_x, rc_y;
			pt::ptree kappa, gamma, lambda;
			pt::ptree beta, beta_a, alpha;
			pt::ptree cell_myo;
			pt::ptree active_myo;
			pt::ptree k;
			pt::ptree fl0;
			pt::ptree fmyo;
			pt::ptree ftension;
			pt::ptree unique_id;
			pt::ptree division_time;
			pt::ptree mother_unique_id;

			id.put("", face_id++);
			original_id.put("", f.id);
			outer.put("", f.outer);
			nsides.put("", f.nsides);
			type.put("", f.data().type_name);
			erased.put("", f.erased);
			A0.put("", f.data().A0);
			P0.put("", f.data().P0);
			for (int vid : verts)
			{
				pt::ptree pvid;
				pvid.put("", vid);
				fverts.push_back(std::make_pair("", pvid));
			}
			
			face.add_child("id", id);
			face.add_child("outer", outer);
			face.add_child("type", type);
			face.add_child("erased", erased);
			face.add_child("A0", A0);
			face.add_child("P0", P0);
			face.add_child("vertices", fverts);
			for (double t : tension)
			{
				pt::ptree pt;
				pt.put("", t);
				ftension.push_back(std::make_pair("", pt));
			}
			face.add_child("tension", ftension);
			faces.push_back(std::make_pair("", face));
			//}
		}
		mesh.add_child("faces", faces);

		// save some system properties in the json file
		out.add_child("mesh", mesh);

		if (ext == "json")
		{
			std::ostringstream oss;
			pt::write_json(oss, out);
			std::regex reg("\\\"([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))\\\"");
			// std::regex reg("\\\"(\\-{0,1}[0-9]+(\\.[0-9]+){0,1})\\\"");
			std::regex reg_bool("\\\"(true|false)\\\"");
			std::string result = std::regex_replace(oss.str(), reg, "$1");
			result = std::regex_replace(result, reg_bool, "$1");

			std::ofstream file;
			file.open(mesh_file);
			file << result;
			file.close();
		}
		else if (ext == "xml")
			pt::write_xml(mesh_file, out, std::locale(), pt::xml_writer_make_settings<string>(' ', 4));
		else
			pt::write_json(mesh_file, out);
	}

	// JSON output
	void Dump::dump_json(const string &json_file)
	{
		ofstream jout(json_file.c_str());
		json j;
		if (_sys.periodic())
			j["mesh"]["box"] = *(_sys.mesh().box());
		j["mesh"]["time_step"] = _sys.time_step();
		j["mesh"]["vertices"] = _sys.mesh().vertices();
		j["mesh"]["faces"] = _sys.mesh().faces();
		jout << setw(2) << j << endl;
		jout.close();
	}

	// JSON conversions
	// Conversation to and from JSON

	// Half-edge
	void to_json(json &j, const HalfEdge<Property> &he)
	{
		j = json{{"from", he.from()->id}, {"to", he.to()->id}};
	}

	void from_json(const json &j, HalfEdge<Property> &he)
	{
		
	}

	// Edge
	void to_json(json &j, const Edge<Property> &e)
	{
		j = json{{"i", e.i}, {"j", e.j}, {"boundary", e.boundary}, {"l0", e.data().l0}};
	}

	void from_json(const json &j, Edge<Property> &e)
	{
		e.i = j.at("i").get<int>();
		e.j = j.at("j").get<int>();
		e.boundary = j.at("boundary").get<bool>();
	}

	// Vertex
	void to_json(json &j, const Vertex<Property> &v)
	{
		vector<int> neigh;
		for (auto he : v.circulator())
			neigh.push_back(he.to()->id);
		 
		string vert_type;
		j = json{
			{"id", v.id},
			{"r", {v.r.x, v.r.y}},
			{"type", v.data().type_name},
			{"erased", v.erased},
			{"boundary", v.boundary},
			{"constraint", v.data().constraint},
			{"force", {v.data().force.x, v.data().force.y}}
			};
	}

	void from_json(const json &j, Vertex<Property> &v)
	{
		v.id = j.at("id").get<int>();
		v.r.x = j.at("r").get<vector<double>>()[0];
		v.r.y = j.at("r").get<vector<double>>()[1];
		v.data().vert_type = j.at("type").get<int>();
		v.erased = j.at("erased").get<bool>();
		v.boundary = j.at("boundary").get<bool>();
		v.data().constraint = j.at("constraint").get<string>();
		v.data().force.x = j.at("force").get<vector<double>>()[0];
		v.data().force.y = j.at("force").get<vector<double>>()[1];
	}

	// Face
	void to_json(json &j, const Face<Property> &f)
	{
		vector<int> verts;
		vector<double> lx, ly, ndx, ndy, l0, myo, tension;
		for (const auto he : f.circulator())
		{
			verts.push_back(he.from()->id);
			tension.push_back(he.data().tension);
		}
		j = json{
			{"id", f.id},
			{"outer", f.outer},
			{"nsides", f.nsides},
			{"type", f.data().type_name},
			{"A0", f.data().A0},
			{"P0", f.data().P0},
			{"vertices", verts},
			{"neighbours", f.data().neighs},
			{"tension", tension}};
	}

	void from_json(const json &j, Face<Property> &f)
	{
		f.id = j.at("id").get<int>();
		f.outer = j.at("outer").get<bool>();
		f.nsides = j.at("nsides").get<int>();
		f.data().face_type = j.at("type").get<int>();
	}

	// Box
	void to_json(json &j, const Box &b)
	{
		vector<double> A = {b.h._mxx, b.h._myx};
		vector<double> B = {b.h._mxy, b.h._myy};
		j = json{
			{"lx", b.h._mxx},
			{"ly", b.h._myy},
			{"a", A},
			{"b", B},
			{"periodic", true}};
	}

	vector<string> split(const string &s, char delim)
	{
		stringstream ss{s};
		string item;
		vector<string> elems;
		while (std::getline(ss, item, delim))
			elems.push_back(move(item));
		return elems;
	}

	void export_Dump(py::module &m)
	{
		py::class_<Dump>(m, "Dump")
			.def(py::init<System &, ForceCompute &>())
			.def("dump_cells", &Dump::dump_cells, py::arg("vtk_file"), py::arg("binary_output") = false, py::arg("draw_periodic") = false)
			.def("dump_junctions", &Dump::dump_junctions, py::arg("vtk_file"), py::arg("binary_output") = false)
			.def("dump_mesh", &Dump::dump_mesh, py::arg("mesh_file"), py::arg("copy_params") = false)
			.def("dump_json", &Dump::dump_json)
			.def("set_sfc", &Dump::set_sfc);
	}
}
