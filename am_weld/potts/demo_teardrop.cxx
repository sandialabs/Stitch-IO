/* Copyright 2019 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525
 * with NTESS, the U.S. Government retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * For more information, contact Jay Lofstead (gflofst@sandeia.gov) or
 * John Mitchell (jamitch@sandia.gov) for more information.
 */ 
#include <vector>
#include <tuple>
#include "am_teardrop.h"

using std::vector;
using RASTER::pool_shape::AM_TEARDROP::Teardrop2D;
using RASTER::pool_shape::AM_TEARDROP::BezierCurve;


vector< vector<double> > get_demo_control_points() {
    int n=9;
    vector< vector<double> > cp(n);
    double a[]={-3./2.,0.0};
    double b[]={a[0],0.5};
    double c[]={0,3*b[1]};
    double d[]={1,2*b[1]};
    double e[]={d[0],0.0};
    double dd[]={d[0],-d[1]};
    double cc[]={c[0],-c[1]};
    double bb[]={b[0],-b[1]};
    double aa[]={a[0], a[1]};
    for(int i=0;i<n;i++)
       cp[i]=vector<double>(2);

    cp[0][0]= a[0];  cp[0][1]= a[1];
    cp[1][0]= b[0];  cp[1][1]= b[1];
    cp[2][0]= c[0];  cp[2][1]= c[1];
    cp[3][0]= d[0];  cp[3][1]= d[1];
    cp[4][0]= e[0];  cp[4][1]= e[1];
    cp[5][0]=dd[0];  cp[5][1]=dd[1];
    cp[6][0]=cc[0];  cp[6][1]=cc[1];
    cp[7][0]=bb[0];  cp[7][1]=bb[1];
    cp[8][0]=aa[0];  cp[8][1]=aa[1];
    return cp;
}

int main(int nargs, char* argv[]){
   vector< vector<double> > cp=get_demo_control_points();
   BezierCurve curve(cp);
   double input_pool_width=100.0;
   Teardrop2D t2d(cp,input_pool_width,true);
   double L0=t2d.get_pool_length();
   double W0=t2d.get_pool_width();
   double x0=t2d.get_pool_position();
   std::cout << "Pool length = " << L0  << std::endl;
   std::cout << "Pool width = " << W0 << std::endl;
   std::cout << "Pool L/W ratio = " << L0/W0 << std::endl;
   std::cout << "Pool position along x-axis = " << x0 << std::endl;
   std::cout << "Parametric coordinate at maximum pool width, u = " 
	   << t2d.get_parametric_coordinate_at_maximum_pool_width() << std::endl;
   return 0;
}
