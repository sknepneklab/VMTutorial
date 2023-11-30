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
		std::map<std::string, vtkSmartPointer<vtkDoubleArray>> force_comp;
		for (auto &fc : _force_compute.factory())
			force_comp[fc.first] = vtkSmartPointer<vtkDoubleArray>::New();
		force_comp["fpull"] = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> face_ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> true_ids = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkIntArray> vert_type = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkDoubleArray> areas = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> perims = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkIntArray> cell_types = vtkSmartPointer<vtkIntArray>::New();
		vtkSmartPointer<vtkDoubleArray> stress = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> num_neigh = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> press = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> sigma_a = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> sigma_p = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> p0 = vtkSmartPointer<vtkDoubleArray>::New();

		ids->SetName("Id");
		ids->SetNumberOfComponents(1);
		true_ids->SetName("TrueId");
		true_ids->SetNumberOfComponents(1);
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
		stress->SetName("Stress");
		stress->SetNumberOfComponents(9);
		time->SetName("time");
		time->SetNumberOfComponents(1);
		time->InsertNextValue(_sys.simulation_time());
		num_neigh->SetName("NumNeigh");
		num_neigh->SetNumberOfComponents(1);
		press->SetName("TotalHydroStress");
		press->SetNumberOfComponents(1);
		sigma_a->SetName("HydroStressArea");
		sigma_a->SetNumberOfComponents(1);
		sigma_p->SetName("HydroStressPerim");
		sigma_p->SetNumberOfComponents(1);
		p0->SetName("realised_p0");
		p0->SetNumberOfComponents(1);

		for (auto &fc : _force_compute.factory())
		{
			force_comp[fc.first]->SetName(("force_" + fc.first).c_str());
			force_comp[fc.first]->SetNumberOfComponents(3);
		}
		force_comp["fpull"]->SetName("force_fpull");
		force_comp["fpull"]->SetNumberOfComponents(3);

		cell_types->SetName("CellTypes");
		cell_types->SetNumberOfComponents(1);

		int id = 0;
		map<int, int> id_map;
		for (auto vh : _sys.mesh().vertices())
		{
			if (!vh.erased)
			{
				points->InsertNextPoint(vh.r.x, vh.r.y, 0.0);
				ids->InsertNextValue(id);
				true_ids->InsertNextValue(vh.id);
				vert_type->InsertNextValue(vh.data().vert_type);
				double f[3] = {vh.data().force.x + vh.data().f_type["fpull"].x, vh.data().force.y + vh.data().f_type["fpull"].y, 0.0};
				forces->InsertNextTuple(f);
				for (auto &fc : _force_compute.factory())
				{
					Vec fcomp = vh.data().f_type[fc.first];
					double fcomponent[3] = {fcomp.x, fcomp.y, 0.0};
					force_comp[fc.first]->InsertNextTuple(fcomponent);
				}
				Vec fcomp = vh.data().f_type["fpull"];
				double fcomponent[3] = {fcomp.x, fcomp.y, 0.0};
				force_comp["fpull"]->InsertNextTuple(fcomponent);
				id_map[vh.id] = id;
				id++;
			}
		}

		for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
		{
			if (!fh->erased)
			{
				double alpha = 0;
				bool omit = false;
				bool crosses_pbc = false;
				vector<int> face_verts;
				if (_sys.periodic())
				{
					HEHandle<Property> he = fh->he();
					HEHandle<Property> first = he;
					Vec r0 = _sys.mesh().get_face_centroid(fh);
					Vec s0 = r0.box->inv_h * r0;
					do
					{
						Vec s = r0.box->inv_h * he->from()->r;
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
								true_ids->InsertNextValue(he->from()->id);
								vert_type->InsertNextValue(he->from()->data().vert_type);
								double f[3] = {0.0, 0.0, 0.0};
								forces->InsertNextTuple(f);
								for (auto &fc : _force_compute.factory())
								{
									double fcomponent[3] = {0.0, 0.0, 0.0};
									force_comp[fc.first]->InsertNextTuple(fcomponent);
								}
								double fcomponent[3] = {0.0, 0.0, 0.0};
								force_comp["fpull"]->InsertNextTuple(fcomponent);
								face_verts.push_back(id);
								id++;
							}
						}
						else
							face_verts.push_back(id_map[he->from()->id]);
						he = he->next();
					} while (he != first);
				}
				if (!(fh->outer || omit))
				{
					face->GetPointIds()->SetNumberOfIds(_sys.mesh().face_sides(fh));
					HEHandle<Property> he = fh->he();
					HEHandle<Property> first = he;

					int i = 0;
					do
					{
						double le = (he->from()->r - he->to()->r).len();
						if (!crosses_pbc)
							face->GetPointIds()->SetId(i++, id_map[he->from()->id]);
						else
						{
							face->GetPointIds()->SetId(i, face_verts[i]);
							i++;
						}
						he = he->next();
					} while (he != first);
					faces->InsertNextCell(face);
					areas->InsertNextValue(_sys.mesh().area(fh));
					perims->InsertNextValue(_sys.mesh().perim(fh));
					face_ids->InsertNextValue(fh->id);
					cell_types->InsertNextValue(fh->data().face_type);
					if (_force_compute.stress_compute())
					{
						stress->InsertNextTuple9(fh->data().stress[0], fh->data().stress[1], 0.0, fh->data().stress[2], fh->data().stress[3], 0.0, 0.0, 0.0, 0.0);
						press->InsertNextValue(fh->data().stress_a[0] + 0.5 * (fh->data().stress_p[0] + fh->data().stress_p[3]));
						sigma_a->InsertNextValue(fh->data().stress_a[0]);
						sigma_p->InsertNextValue(0.5 * (fh->data().stress_p[0] + fh->data().stress_p[3]));
					}
					num_neigh->InsertNextValue(fh->nsides);
					p0->InsertNextValue(_sys.mesh().perim(fh) / sqrt(_sys.mesh().area(fh)));
				}
			}
		}

		polydata->SetPoints(points);
		polydata->GetPointData()->AddArray(ids);
		polydata->GetPointData()->AddArray(true_ids);
		polydata->GetPointData()->AddArray(vert_type);
		polydata->GetPointData()->AddArray(forces);

		for (auto &fc : _force_compute.factory())
			polydata->GetPointData()->AddArray(force_comp[fc.first]);
		polydata->GetPointData()->AddArray(force_comp["fpull"]);

		polydata->SetPolys(faces);
		polydata->GetCellData()->AddArray(areas);
		polydata->GetCellData()->AddArray(perims);
		polydata->GetCellData()->AddArray(face_ids);
		polydata->GetCellData()->AddArray(cell_types);
		polydata->GetCellData()->AddArray(num_neigh);
		polydata->GetCellData()->AddArray(p0);
		if (_force_compute.stress_compute())
		{
			polydata->GetCellData()->AddArray(stress);
			polydata->GetCellData()->AddArray(press);
			polydata->GetCellData()->AddArray(sigma_a);
			polydata->GetCellData()->AddArray(sigma_p);
		}
		polydata->GetFieldData()->AddArray(time);

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
		vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
		std::map<std::string, vtkSmartPointer<vtkDoubleArray>> tension_type;
		for (auto &fc : _force_compute.factory())
			tension_type[fc.first] = vtkSmartPointer<vtkDoubleArray>::New();

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
		time->SetName("time");
		time->SetNumberOfComponents(1);
		time->InsertNextValue(_sys.simulation_time());
		for (auto &fc : _force_compute.factory())
		{
			tension_type[fc.first]->SetName(("tension_" + fc.first).c_str());
			tension_type[fc.first]->SetNumberOfComponents(1);
		}

		int id = 0;
		map<int, int> id_map;
		map<int, bool> omit;

		// build list of points belonging to faces, with shift (for display purpose)
		for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
		{
			if (!fh->erased)
			{
				HEHandle<Property> he = fh->he();
				HEHandle<Property> first = he;
				Vec fc = _sys.mesh().get_face_centre(fh); // face centre coordinates

				do
				{
					EdgeHandle<Property> eh = he->edge();
					if (!eh->erased)
					{
						VertexHandle<Property> vh1 = he->from();
						VertexHandle<Property> vh2 = he->to();
						Vec r1_shift = fc + _sfc * (vh1->r - fc);
						Vec r2_shift = fc + _sfc * (vh2->r - fc);
						if (fh->outer)
						{
							r1_shift = vh1->r;
							r2_shift = vh2->r;
						}
						omit[he->idx()] = false;
						if (_sys.periodic())
						{
							Vec s1 = r1_shift.box->inv_h * r1_shift;
							Vec s2 = r2_shift.box->inv_h * r2_shift;
							double sx = s1.x - s2.x, sy = s1.y - s2.y;
							if (fabs(rint(sx)) != 0.0 || fabs(rint(sy)) != 0.0)
								omit[he->idx()] = true;
						}
						if (!(omit[he->idx()]))
						{
							points->InsertNextPoint(r1_shift.x, r1_shift.y, 0.0);
							ids->InsertNextValue(id);
							true_ids->InsertNextValue(vh1->id);
							id_map[vh1->id] = id;
							id++;

							points->InsertNextPoint(r2_shift.x, r2_shift.y, 0.0);
							ids->InsertNextValue(id);
							true_ids->InsertNextValue(vh2->id);
							id_map[vh2->id] = id;
							id++;
						}
					}
					he = he->next();
				} while (he != first);
			}
		}

		polydata->SetPoints(points);
		polydata->GetPointData()->AddArray(ids);
		polydata->GetPointData()->AddArray(true_ids);

		// create connectivity
		id = 0;
		for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
		{
			if (!fh->erased)
			{
				HEHandle<Property> he = fh->he();
				HEHandle<Property> first = he;
				do
				{
					EdgeHandle<Property> eh = he->edge();
					if (!(eh->erased || omit[he->idx()]))
					{
						edge->GetPointIds()->SetId(0, id);
						id++;
						edge->GetPointIds()->SetId(1, id);
						id++;
						he_id->InsertNextValue(he->idx());
						edge_id->InsertNextValue(eh->idx());
						lines->InsertNextCell(edge);
						lens->InsertNextValue(_sys.mesh().len(eh));
						double he_tension = he->data().tension;
						if (fabs(he_tension) < 1e-10)
							he_tension = _force_compute.tension(he, he->data().myo, eh->data().l0);
						tension->InsertNextValue(he_tension);
						l0->InsertNextValue(eh->data().l0);
						edge_type->InsertNextValue(eh->data().edge_type);
						for (auto &fc : _force_compute.factory())
							tension_type[fc.first]->InsertNextValue(_force_compute.tension_of_type(fc.first, he, he->data().myo, eh->data().l0));
					}
					he = he->next();
				} while (he != first);
			}
		}

		polydata->SetLines(lines);
		polydata->GetCellData()->AddArray(lens);
		polydata->GetCellData()->AddArray(he_id);
		polydata->GetCellData()->AddArray(edge_id);
		polydata->GetCellData()->AddArray(tension);
		// polydata->GetCellData()->AddArray(myo);
		polydata->GetCellData()->AddArray(l0);
		polydata->GetCellData()->AddArray(edge_type);
		polydata->GetFieldData()->AddArray(time);
		for (auto &fc : _force_compute.factory())
			polydata->GetCellData()->AddArray(tension_type[fc.first]);

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

	void Dump::dump_box(const string &vtk_file, bool binary_output)
	{
		if (!_sys.box())
			throw runtime_error("Simulation box not defined.");

		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

		vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();

		Vec r1 = _sys.box()->h * Vec(1, 0, nullptr);
		Vec r2 = _sys.box()->h * Vec(0, 1, nullptr);
		double xc = 0.5 * (r1.x + r2.x), yc = 0.5 * (r1.y + r2.y);

		points->InsertNextPoint(-xc, -yc, 0.0);
		points->InsertNextPoint(-xc + r1.x, -yc + r1.y, 0.0);
		points->InsertNextPoint(-xc + r1.x + r2.x, -yc + r1.y + r2.y, 0.0);
		points->InsertNextPoint(-xc + r2.x, -yc + r2.y, 0.0);

		polydata->SetPoints(points);

		edge->GetPointIds()->SetId(0, 0);
		edge->GetPointIds()->SetId(1, 1);
		lines->InsertNextCell(edge);
		edge->GetPointIds()->SetId(0, 1);
		edge->GetPointIds()->SetId(1, 2);
		lines->InsertNextCell(edge);
		edge->GetPointIds()->SetId(0, 2);
		edge->GetPointIds()->SetId(1, 3);
		lines->InsertNextCell(edge);
		edge->GetPointIds()->SetId(0, 3);
		edge->GetPointIds()->SetId(1, 0);
		lines->InsertNextCell(edge);

		polydata->SetLines(lines);
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
		if (copy_params)
			_force_compute.copy_type_param_to_cell();
		vector<string> strs = split(mesh_file, '.');
		string ext = strs[strs.size() - 1];
		transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c)
				  { return tolower(c); });
		pt::ptree out;
		pt::ptree rng;
		if (_sys.has_rng_state())
			for (auto &rng_state : _sys.get_all_rng_states())
			{
				pt::ptree rng_param;
				pt::ptree name;
				pt::ptree seed;
				pt::ptree u_dist;
				pt::ptree n_dist;
				name.put("", rng_state.first);
				seed.put("", rng_state.second.seed);
				u_dist.put("", rng_state.second.u_dist);
				n_dist.put("", rng_state.second.n_dist);
				rng_param.add_child("name", name);
				rng_param.add_child("seed", seed);
				rng_param.add_child("u_dist", u_dist);
				rng_param.add_child("n_dist", n_dist);
				rng.push_back(std::make_pair("", rng_param));
			}
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
		for (auto vh : _sys.mesh().vertices())
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
			id.put("", vh.id);
			boundary.put("", vh.boundary);
			type.put("", vh.data().type_name);
			x.put("", vh.r.x);
			y.put("", vh.r.y);
			r.push_back(std::make_pair("", x));
			r.push_back(std::make_pair("", y));
			fx.put("", vh.data().force.x + vh.data().f_type["fpull"].x);
			fy.put("", vh.data().force.y + vh.data().f_type["fpull"].y);
			force.push_back(std::make_pair("", fx));
			force.push_back(std::make_pair("", fy));
			vx.put("", vh.data().vel.x);
			vy.put("", vh.data().vel.y);
			velocity.push_back(std::make_pair("", vx));
			velocity.push_back(std::make_pair("", vy));
			constraint.put("", vh.data().constraint);
			erased.put("", vh.erased);
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
		_sys.mesh().populate_face_neighbours();
		int face_id = 0;
		for (auto fh : _sys.mesh().faces())
		{
			// if (!fh.erased)
			//{
			vector<int> verts;
			vector<double> l0, myo, tension;
			if (!fh.erased)
			{
				HECHandle<Property> he = fh.he();
				HECHandle<Property> first = he;
				do
				{
					verts.push_back(he->from()->id);
					l0.push_back(he->edge()->data().l0);
					myo.push_back(he->data().myo);
					Vec er = he->to()->r - he->from()->r;
					Vec fr = he->pair()->face()->data().rc - fh.data().rc;
					tension.push_back(he->data().tension);
				} while ((he = he->next()) != first);
			}
			/*
			else
			{
			  verts.push_back(-1);
			  l0.push_back(-1);
			  myo.push_back(-1);
			  tension.push_back(-1);
			}
			*/
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
			original_id.put("", fh.id);
			outer.put("", fh.outer);
			nsides.put("", fh.nsides);
			type.put("", fh.data().type_name);
			erased.put("", fh.erased);
			A0.put("", fh.data().A0);
			P0.put("", fh.data().P0);
			press.put("", fh.data().stress_a[0]);
			for (int vid : verts)
			{
				pt::ptree pvid;
				pvid.put("", vid);
				fverts.push_back(std::make_pair("", pvid));
			}
			for (int nid : fh.data().neighs)
			{
				pt::ptree pnid;
				pnid.put("", nid);
				fneighs.push_back(std::make_pair("", pnid));
			}
			rc_x.put("", fh.data().rc.x);
			rc_y.put("", fh.data().rc.y);
			rc.push_back(std::make_pair("", rc_x));
			rc.push_back(std::make_pair("", rc_y));
			kappa.put("", fh.data().kappa);
			gamma.put("", fh.data().gamma);
			lambda.put("", fh.data().lambda);
			beta.put("", fh.data().beta);
			beta_a.put("", fh.data().beta_a);
			alpha.put("", fh.data().alpha);
			cell_myo.put("", fh.data().cell_myo);
			active_myo.put("", fh.data().active_myo);
			k.put("", fh.data().k);
			unique_id.put("", fh.data().unique_id);
			division_time.put("", fh.data().division_time);
			mother_unique_id.put("", fh.data().mother_unique_id);

			face.add_child("id", id);
			face.add_child("original_id", original_id);
			face.add_child("outer", outer);
			face.add_child("nsides", nsides);
			face.add_child("type", type);
			face.add_child("erased", erased);
			face.add_child("A0", A0);
			face.add_child("P0", P0);
			face.add_child("pressure", press);
			face.add_child("vertices", fverts);
			face.add_child("neighbours", fneighs);
			face.add_child("rc", rc);
			face.add_child("kappa", kappa);
			face.add_child("gamma", gamma);
			face.add_child("lambda", lambda);
			face.add_child("beta", beta);
			face.add_child("beta_a", beta_a);
			face.add_child("alpha", alpha);
			face.add_child("cell_myo", cell_myo);
			face.add_child("active_myo", active_myo);
			face.add_child("k", k);
			face.add_child("unique_id", unique_id);
			face.add_child("division_time", division_time);
			face.add_child("mother_unique_id", mother_unique_id);
			for (double ll0 : l0)
			{
				pt::ptree pl0;
				pl0.put("", ll0);
				fl0.push_back(std::make_pair("", pl0));
			}
			face.add_child("l0", fl0);
			for (double m : myo)
			{
				pt::ptree pm;
				pm.put("", m);
				fmyo.push_back(std::make_pair("", pm));
			}
			face.add_child("myo", fmyo);
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
		pt::ptree system;
		pt::ptree time;
		pt::ptree time_step;
		pt::ptree symmetric_myosin;
		pt::ptree has_spokes;
		time.put("", _sys.simulation_time());
		system.add_child("time", time);
		time_step.put("", _sys.time_step());
		system.add_child("time_step", time_step);
		symmetric_myosin.put("", _sys.symmetric_myosin());
		system.add_child("symmetric_myosin", symmetric_myosin);
		has_spokes.put("", _sys.has_spokes());
		system.add_child("has_spokes", has_spokes);
		out.add_child("system", system);
		if (_sys.has_rng_state())
			out.add_child("RNG", rng);
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
		_sys.mesh().populate_face_neighbours();
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
		j = json{{"from", he.from()->id}, {"to", he.to()->id}, {"myo", he.data().myo}};
	}

	void from_json(const json &j, HalfEdge<Property> &he)
	{
		he.data().myo = j.at("myo").get<double>();
	}

	// Edge
	void to_json(json &j, const Edge<Property> &e)
	{
		j = json{{"i", e.i}, {"j", e.j}, {"boundary", e.boundary}, {"myo", 0.5 * (e.he()->data().myo + e.he()->pair()->data().myo)}, {"l0", e.data().l0}};
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
		HECHandle<Property> he = v.he();
		HECHandle<Property> first = v.he();
		do
		{
			neigh.push_back(he->to()->id);
			he = he->pair()->next();
		} while (he != first);
		string vert_type;
		j = json{
			{"id", v.id},
			{"r", {v.r.x, v.r.y}},
			{"type", v.data().type_name},
			{"erased", v.erased},
			{"boundary", v.boundary},
			{"constraint", v.data().constraint},
			{"coordination", v.coordination},
			{"force", {v.data().force.x, v.data().force.y}},
			{"neighbours", neigh}};
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
		v.coordination = j.at("coordination").get<int>();
		v.data().force.x = j.at("force").get<vector<double>>()[0];
		v.data().force.y = j.at("force").get<vector<double>>()[1];
	}

	// Face
	void to_json(json &j, const Face<Property> &f)
	{
		vector<int> verts;
		HECHandle<Property> he = f.he();
		HECHandle<Property> first = he;
		vector<double> lx, ly, ndx, ndy, l0, myo, tension;
		do
		{
			verts.push_back(he->from()->id);
			l0.push_back(he->edge()->data().l0);
			myo.push_back(he->data().myo);
			Vec er = he->to()->r - he->from()->r;
			Vec fr = he->pair()->face()->data().rc - f.data().rc;
			lx.push_back(er.x);
			ly.push_back(er.y);
			ndx.push_back(fr.x);
			ndy.push_back(fr.y);
			tension.push_back(he->data().tension);
		} while ((he = he->next()) != first);
		j = json{
			{"id", f.id},
			{"outer", f.outer},
			{"nsides", f.nsides},
			{"type", f.data().type_name},
			{"A0", f.data().A0},
			{"P0", f.data().P0},
			{"vertices", verts},
			{"neighbours", f.data().neighs},
			{"rc", {f.data().rc.x, f.data().rc.y}},
			/*
			{"lx", lx},
			{"ly", ly},
			{"ndx", ndx},
			{"ndy", ndy},
			*/
			{"kappa", f.data().kappa},
			{"gamma", f.data().gamma},
			{"lambda", f.data().lambda},
			{"beta", f.data().beta},
			{"beta_a", f.data().beta_a},
			{"alpha", f.data().alpha},
			{"cell_myo", f.data().cell_myo},
			{"active_myo", f.data().active_myo},
			{"k", f.data().k},
			{"l0", l0},
			{"myo", myo},
			{"tension", tension}};
	}

	void from_json(const json &j, Face<Property> &f)
	{
		f.id = j.at("id").get<int>();
		f.outer = j.at("outer").get<bool>();
		f.nsides = j.at("nsides").get<int>();
		f.data().face_type = j.at("type").get<int>();
		f.data().kappa = j.at("kappa").get<double>();
		f.data().gamma = j.at("gamma").get<double>();
		f.data().lambda = j.at("lambda").get<double>();
		f.data().beta = j.at("beta").get<double>();
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
			.def("dump_box", &Dump::dump_box, py::arg("vtk_file"), py::arg("binary_output") = false)
			.def("dump_mesh", &Dump::dump_mesh, py::arg("mesh_file"), py::arg("copy_params") = false)
			.def("dump_json", &Dump::dump_json)
			.def("set_sfc", &Dump::set_sfc);
	}
}
