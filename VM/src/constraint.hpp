/*!
 * \file constraint.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Abstract constraints class
 */

#ifndef __CONSTRAINT_HPP__
#define __CONSTRAINT_HPP__

#include "system.hpp"

namespace VMTutorial
{

	class Constraint
	{
	public:
		Constraint() {}
		virtual ~Constraint() {}
		virtual Vec apply(const Vertex<Property> &, const Vec &) = 0;
		virtual Vec apply(const Vec &) = 0;
	};

}
#endif