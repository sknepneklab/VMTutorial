/*!
 * \file matrix.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Nov-2023
 * \brief 2x2 matrix manipulations
 */

#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <vector>
#include <stdexcept>
#include <cmath>

using std::fabs;
using std::runtime_error;
using std::vector;

namespace VMTutorial
{

	class Matrix
	{
	public:
		Matrix() : _mxx{0}, 
				   _mxy{0}, 
				   _myx{0}, 
				   _myy{0} 
		{

		}
		Matrix(double mxx, double mxy, double myx, double myy) : _mxx{mxx}, 
																 _mxy{mxy}, 
																 _myx{myx}, 
																 _myy{myy} 
		{
			
		}
		Matrix(const vector<vector<double>> &m)
		{
			if (m.size() != 2)
				throw runtime_error("Matrix m has have 2 rows.");
			if (m[0].size() != 2 || m[1].size() != 2)
				throw runtime_error("Each row of matrix m has to have 2 columns.");
			_mxx = m[0][0];
			_mxy = m[0][1];
			_myx = m[1][0];
			_myy = m[1][1];
		}

		Matrix(const vector<double> &r1, const vector<double> &r2)
		{
			if (r1.size() != 2 || r2.size() != 2)
				throw runtime_error("Each row of matrix m has to have 2 columns.");
			_mxx = r1[0];
			_mxy = r1[1];
			_myx = r2[0];
			_myy = r2[1];
		}
		Matrix(const Matrix &m)
		{
			_mxx = m._mxx;
			_mxy = m._mxy;
			_myx = m._myx;
			_myy = m._myy;
		}

		double trace() { return (_mxx + _myy); }
		double trace() const { return (_mxx + _myy); }
		double det() { return (_mxx * _myy - _mxy * _myx); }
		double det() const { return (_mxx * _myy - _mxy * _myx); }
		Matrix inv();
		Matrix T()
		{
			return Matrix(_mxx, _myx, _mxy, _myy);
		}

		Matrix operator*(const Matrix &);

		double _mxx, _mxy, _myx, _myy;
	};

}
#endif