// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1985, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include <cmath>
#include <limits>

double MACHEP = 1.11022302462515654042E-16;    /* 2**-53 */

static double A[8] = {
    4.65128586073990045278E-5,
    7.31589045238094711071E-3,
    1.33847639578309018650E-1,
    8.79691311754530315341E-1,
    2.71149851196553469920E0,
    4.25697156008121755724E0,
    3.29771340985225106936E0,
    1.00000000000000000126E0,
};

static double B[8] = {
    6.90990488912553276999E-4,
    2.54043763932544379113E-2,
    2.82974860602568089943E-1,
    1.41172597751831069617E0,
    3.63800533345137075418E0,
    5.03278880143316990390E0,
    3.54771340985225096217E0,
    9.99999999999999998740E-1,
};

double polevl(double x, double coef[], int N)
{
    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;

    do
        ans = ans * x + *p++;
    while (--i);

    return (ans);
}

double spence(double x)
{
    double w, y, z;
    int flag;

    if (x < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    if (x == 1.0)
        return (0.0);

    if (x == 0.0)
        return (M_PI * M_PI / 6.0);

    flag = 0;

    if (x > 2.0) {
        x = 1.0 / x;
        flag |= 2;
    }

    if (x > 1.5) {
        w = (1.0 / x) - 1.0;
        flag |= 2;
    }

    else if (x < 0.5) {
        w = -x;
        flag |= 1;
    }

    else
        w = x - 1.0;


    y = -w * polevl(w, A, 7) / polevl(w, B, 7);

    if (flag & 1)
        y = (M_PI * M_PI) / 6.0 - std::log(x) * std::log(1.0 - x) - y;

    if (flag & 2) {
        z = std::log(x);
        y = -0.5 * z * z - y;
    }

    return (y);
}
