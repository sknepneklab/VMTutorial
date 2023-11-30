
/*!
 * \file base_property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief BaseProperty class
 */

#ifndef __BASE_PROPERTY_HPP__
#define __BASE_PROPERTY_HPP__

namespace VMTutorial
{
	struct BaseProperty
	{
		struct HEProperty
		{
			int he_type = -1;
		};
		struct VertexProperty
		{
			int vert_type = -1;
			VertexProperty &operator=(const VertexProperty &p)
			{
				if (this == &p)
					return *this;
				this->vert_type = p.vert_type;
				return *this;
			}
		};
		struct EdgeProperty
		{
			int edge_type = -1;
		};
		struct FaceProperty
		{
			int face_type = -1;
			FaceProperty &operator=(const FaceProperty &p)
			{
				if (this == &p)
					return *this;
				this->face_type = p.face_type;
				return *this;
			}
		};
	};
}

#endif
