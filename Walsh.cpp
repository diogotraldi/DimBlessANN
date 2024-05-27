#include "stdafx.h"
# include <cmath>
#include "Walsh.h"

namespace DimBless
{
	void Walsh::Init(int N)
	{
		n = N;
		y = new float[n];

		n2 = n / 2;
		m = i4_log_2(n);
	}

	void Walsh::Close()
	{
		if (y != nullptr)
		{
			delete[] y;
			y = nullptr;
		}
	}

	void Walsh::fwt(float x[])

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    FWT performs a fast Walsh transform.
		//
		//  Discussion:
		//
		//    This routine performs a fast Walsh transform on an input series X
		//    leaving the transformed results in X. 
		//    X is dimensioned N, which must be a power of 2.
		//    The results of this Walsh transform are in sequency order.
		//
		//    The output sequence could be normalized by dividing by N.
		//
		//    Note that the program text in the reference included the line
		//      y(jd) = abs ( x(j) - x(j2) )
		//    which has been corrected to:
		//      y(jd) = x(j) - x(j2)
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    16 March 2011
		//
		//  Author:
		//
		//    Ken Beauchamp
		//
		//  Reference:
		//
		//    Ken Beauchamp,
		//    Walsh functions and their applications,
		//    Academic Press, 1975,
		//    ISBN: 0-12-084050-2,
		//    LC: QA404.5.B33.
		//
		//  Parameters:
		//
		//    Input, int N, the number of items in X.
		//    N must be a power of 2.
		//
		//    Input/output, double X[N], the data to be transformed.
		//
	{
		int i;
		int j;
		int j2;
		int jd;
		int js;
		int l;
		int nx;
		int ny;
		int nz;
		int nzi;
		int nzn;

		for (l = 1; l <= m; l++)
		{
			ny = 0;
			nz = i4_power(2, l - 1);
			nzi = 2 * nz;
			nzn = n / nzi;
			for (i = 1; i <= nzn; i++)
			{
				nx = ny + 1;
				ny = ny + nz;
				js = (i - 1) * nzi;
				jd = js + nzi + 1;
				for (j = nx; j <= ny; j++)
				{
					js = js + 1;
					j2 = j + n2;
					y[js - 1] = x[j - 1] + x[j2 - 1];
					jd = jd - 1;
					y[jd - 1] = x[j - 1] - x[j2 - 1];
				}
			}
			r8vec_copy(n, y, x);
		}

		return;
	}
	//****************************************************************************80

	int Walsh::i4_log_2(int i)

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
		//
		//  Example:
		//
		//        I  I4_LOG_10
		//    -----  --------
		//        0    0
		//        1    0
		//        2    1
		//        3    1
		//        4    2
		//        5    2
		//        7    2
		//        8    3
		//        9    3
		//     1000    9
		//     1024   10
		//
		//  Discussion:
		//
		//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    04 January 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int I, the number whose logarithm base 2 is desired.
		//
		//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
		//    the absolute value of X.
		//
	{
		int i_abs;
		int two_pow;
		int value;

		if (i == 0)
		{
			value = 0;
		}
		else
		{
			value = 0;
			two_pow = 2;

			i_abs = abs(i);

			while (two_pow <= i_abs)
			{
				value = value + 1;
				two_pow = two_pow * 2;
			}
		}

		return value;
	}
	//****************************************************************************80

	int Walsh::i4_power(int i, int j)

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    I4_POWER returns the value of I^J.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    01 April 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int I, J, the base and the power.  J should be nonnegative.
		//
		//    Output, int I4_POWER, the value of I^J.
		//
	{
		int k;
		int value;

		if (j < 0)
		{
			if (i == 1)
			{
				value = 1;
			}
			else if (i == 0)
			{
				return 1;
			}
			else
			{
				value = 0;
			}
		}
		else if (j == 0)
		{
			if (i == 0)
			{
				return 1;
			}
			else
			{
				value = 1;
			}
		}
		else if (j == 1)
		{
			value = i;
		}
		else
		{
			value = 1;
			for (k = 1; k <= j; k++)
			{
				value = value * i;
			}
		}
		return value;
	}
	//****************************************************************************80

	void Walsh::r8vec_copy(int n, float a1[], float a2[])

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    R8VEC_COPY copies an R8VEC.
		//
		//  Discussion:
		//
		//    An R8VEC is a vector of R8's.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    03 July 2005
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int N, the number of entries in the vectors.
		//
		//    Input, double A1[N], the vector to be copied.
		//
		//    Output, double A2[N], the copy of A1.
		//
	{
		int i;

		for (i = 0; i < n; i++)
		{
			a2[i] = a1[i];
		}
		return;
	}
	//**
}