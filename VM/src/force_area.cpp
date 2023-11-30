/* ***************************************************************************
 *
 *  Copyright (C) 2017 University of Dundee
 *  All rights reserved.
 *
 *  This file is part of AJM (Active Junction Model) program.
 *
 *  AJM is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  AJM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class
 */

#include "force_area.hpp"

namespace AJM
{
	Vec ForceArea::compute(const Vertex<Property> &v, const HalfEdge<Property> &he)
	{
		Vec l = he.to()->r - v.r;					// vector along the junction pointing away from the vertex
		Face<Property> &f = *(he.face());			// cell to the right of the half edge
		Face<Property> &fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
		double A1 = _sys.mesh().area(f);
		double A2 = _sys.mesh().area(fp);
		double A0_1 = f.data().A0;
		double A0_2 = fp.data().A0;

		double kappa_1 = (f.outer) ? 0.0 : _kappa[f.data().face_type];
		double kappa_2 = (fp.->outer) ? 0.0 : _kappa[f.data().face_type];

		Vec farea_vec = 0.5 * (kappa_1 * (A1 - A0_1) - kappa_2 * (A2 - A0_2)) * l.ez_cross_v();

		return farea_vec;
	}

	double ForceArea::energy(const Face<Property> &f)
	{
		double A = _sys.mesh().area(f.id);
		double A0 = f.data().A0;
		double kappa;

		if (f.outer || f.erased)
			return 0.0;

		double kappa = _kappa[f.data().face_type];

		double dA = A - A0;
		return 0.5 * kappa * dA * dA;
	}
}
