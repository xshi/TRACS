// GNU General Public License Agreement
// Copyright (C) 2004-2010 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by CodeCogs. 
// You must retain a copy of this licence in all copies. 
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
// ---------------------------------------------------------------------------------
//! Linearly interpolates a given set of points.

#ifndef MATHS_INTERPOLATION_LINEAR_H
#define MATHS_INTERPOLATION_LINEAR_H

namespace Maths
{

namespace Interpolation
{

//! Linearly interpolates a given set of points.

class Linear
{
public:

//! Class constructor

Linear(int n, double *x, double *y)
{

            m_x = new double[n];
            m_y = new double[n];

            for (int i = 0; i < n; ++i) {
                m_x[i] = x[i];
                m_y[i] = y[i];
            }

        }

//! Class destructor

~Linear()
{

            delete [] m_x;
            delete [] m_y;

        }

//! Returns an interpolated value.

double getValue(double x)
{

            int i = 0;
            while (x > m_x[++i]);

            double a = (x - m_x[i - 1]) / (m_x[i] - m_x[i - 1]);
            return m_y[i - 1] + a * (m_y[i] - m_y[i - 1]);

        }

private:

double *m_x, *m_y;
};


//! A static function implementing the Linear Class for one off calculations

double Linear_once(int N, double *x, double *y, double a )
{
  // This function is created to enable an Instant Calculator on CodeCogs. 
  // You probably shouldn't be using this function otherwise. 

   Maths::Interpolation:: Linear A(N, x, y);
   return A.getValue(a);
}

}

}

#endif

