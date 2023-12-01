/*!
 * \file constraint_fixed.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Fixed vertex constraint
 */

#ifndef __CONSTRAINT_FIXED_HPP__
#define __CONSTRAINT_FIXED_HPP__

#include "constraint.hpp"

namespace VMTutorial
{
	class ConstraintFixed : public Constraint
	{
	public:
		ConstraintFixed() {}
		Vec apply(const Vertex<Property> &v, const Vec &f) override
		{
			return Vec(0.0, 0.0);
		}
		Vec apply(const Vec &f) override
		{
			return Vec(0.0, 0.0);
		}
	};

}
#endif