#include <vector>
#include <tuple>
#include "teardrop.h"

using std::vector;
using weld::pool_shape::teardrop::Teardrop2D;
using weld::pool_shape::BezierCurve;
using weld::pool_shape::teardrop::get_pool_position;
using weld::pool_shape::teardrop::get_parametric_coordinate_at_maximum_pool_width;
using weld::pool_shape::teardrop::predict_parametric_coordinate;


vector< vector<double> > get_demo_control_points() {
    int n=9;
    vector< vector<double> > cp(n);
    double a[]={0.0,-3./2.};
    double b[]={0.5,a[1]};
    double c[]={3*b[0],0};
    double d[]={2*b[0],1};
    double e[]={0.0,d[1]};
    double dd[]={-d[0],d[1]};
    double cc[]={-c[0],c[1]};
    double bb[]={-b[0],b[1]};
    double aa[]={a[0],a[1]};
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
   Teardrop2D t2d(cp);
   double L0=t2d.get_pool_length();
   double W0=t2d.get_pool_width();
   double y0=t2d.get_pool_position();
   std::cout << "Pool length = " << L0  << std::endl;
   std::cout << "Pool width = " << W0 << std::endl;
   std::cout << "Pool L/W ratio = " << L0/W0 << std::endl;
   std::cout << "Pool position along y-axis = " << y0 << std::endl;
   std::cout << "Parametric coordinate at maximum pool width, u = " 
	   << get_parametric_coordinate_at_maximum_pool_width(curve) << std::endl;
   double L=1.0;
   double alpha=L/L0;
   double dy=0.0;
   dy-=y0;
   // Move pool returns new 'position'
   double y=t2d.move_pool(dy);
   t2d.scale_pool_size(alpha);
   L=t2d.get_pool_length();
   double W=t2d.get_pool_width();
   std::cout << "UPDATED Pool length = " << L  << std::endl;
   std::cout << "UPDATED Pool width = " << W << std::endl;
   std::cout << "UPDATED Pool L/W ratio = " << L/W << std::endl;
   std::cout << "UPDATE Pool position along y-axis = " << y << std::endl;
   std::cout << "Parametric coordinate at maximum pool width, u = " 
	   << get_parametric_coordinate_at_maximum_pool_width(curve) << std::endl;
   //double x0=0.5, y0=-0.0;
   //double u, umin, umax;
   //std::tie(u,umin,umax)=predict_parametric_coordinate(curve,x0,y0);
   //std::cout << "Predict parametric coordinate: u=" << u << "; umin=" << umin << "; umax=" << umax << std::endl;
   return 0;
}
