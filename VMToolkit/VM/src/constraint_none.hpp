/*!
 * \file constraint_none.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief No constraint
 */

#ifndef __CONSTRAINT_NONE_HPP__
#define __CONSTRAINT_NONE_HPP__

#include "constraint.hpp"

namespace VMTutorial
{

	class ConstraintNone : public Constraint
	{
	public:
		ConstraintNone() {}
		Vec apply(const Vertex<Property> &v, const Vec &f) override
		{
			return Vec(f);
		}
		Vec apply(const Vec &f) override
		{
			return Vec(f);
		}
	};

}
#endif