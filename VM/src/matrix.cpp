/*!
 * \file matrix.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Nov-2023
 * \brief 2x2 matrix manipulations
 */

#include "matrix.hpp"

namespace VMTutorial
{
	Matrix Matrix::inv()
	{
		double D = this->det();
		if (fabs(D) <= 1e-7)
			throw runtime_error("Encountered singular matrix.");
		return Matrix(_myy / D, -_mxy / D, -_myx / D, _mxx / D);
	}

	Matrix Matrix::operator*(const Matrix &m)
	{
		return Matrix(_mxx * m._mxx + _mxy * m._myx,
					  _mxx * m._mxy + _mxy * m._myy,
					  _myx * m._mxx + _myy * m._myx,
					  _myx * m._mxy + _myy * m._myy);
	}

}
