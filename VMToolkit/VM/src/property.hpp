/*!
 * \file property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Nov-2023
 * \brief Property class
 */

#ifndef __PROPERTY_HPP__
#define __PROPERTY_HPP__

#include <map>
#include <vector>
#include <string>

#include "vec.hpp"

#include "base_property.hpp"

using std::map;
using std::string;
using std::vector;

namespace VMTutorial
{

	struct Property : public BaseProperty
	{
		struct HEProperty : public BaseProperty::HEProperty
		{
			double tension = 0.0;
			double l0;					 // Native length of an edge
			int old_face_id;			 // the way to distinguish if a vertex split was actual T1 or a bounce back
			map<string, Vec> force_type; // Force of a given type on the he->from() vertex due to along this half-edge
		};
		struct VertexProperty : public BaseProperty::VertexProperty
		{
			Vec vel;
			Vec force;
			string constraint;		 // if "x" move only along x-axis, if "y" move only along y axis, if "radial", move along radius; otherwise ignore
			map<string, Vec> f_type; // Force from a given interaction type (e.g., area term, perimeter term, etc)
			string type_name;		 // String with the actual name of the the vertex type
			VertexProperty &operator=(const VertexProperty &p)
			{
				if (this == &p)
					return *this;
				this->vert_type = p.vert_type;
				this->type_name = p.type_name;
				return *this;
			}
		};
		struct EdgeProperty : public BaseProperty::EdgeProperty
		{
			double tension = 0.0;
			double l0; // Native length of an edge
		};
		struct FaceProperty : public BaseProperty::FaceProperty
		{
			double A0;
			double P0;
			int unique_id;
			int original_face;	// one of the faces collapsed edge belonged to
			double kappa;		// area modulus
			double gamma;		// perimeter modulus
			double lambda;		// line tension
			Vec n;				// Self-propulsion direction
			Vec rc;				// centre of the face
			vector<int> neighs; // indices of neighbouring faces
			string type_name;	// String with the actual name of the the cell type
			Vec v;				// Cell velocity (average velocity of its vertices)
			FaceProperty &operator=(const FaceProperty &p)
			{
				if (this == &p)
					return *this;
				this->face_type = p.face_type;
				this->A0 = p.A0;
				this->P0 = p.P0;
				this->type_name = p.type_name;
				return *this;
			}
		};
	};

} // namespace AJM

#endif
