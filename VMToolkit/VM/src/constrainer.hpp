/*!
 * \file constrainer.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Applies constraint to a vertex
 */


#ifndef __CONSTRINER_HPP__
#define __CONSTRINER_HPP__

#include "class_factory.hpp"

#include "class_factory.hpp"
#include "constraint.hpp"

namespace VMTutorial
{

	class Constrainer : public ClassFactory<Constraint>
	{
	public:
		Constrainer() = default;
		~Constrainer() = default;

		Constrainer(const Constrainer &) = delete;

		Vec apply_vertex(Vertex<Property> &v, const Vec &fc)
		{
			if (this->factory_map.find(v.data().constraint) != this->factory_map.end())
				return this->factory_map[v.data().constraint]->apply(v, fc);
			else
				return Vec(fc);
		}

		Vec apply_vector(Vertex<Property> &v, const Vec &fc)
		{
			if (this->factory_map.find(v.data().constraint) != this->factory_map.end())
				return this->factory_map[v.data().constraint]->apply(fc);
			else
				return Vec(fc);
		}
	};

}
#endif
