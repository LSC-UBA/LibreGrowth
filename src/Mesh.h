#ifndef MESH_H_
#define MESH_H_
#include <math.h> 
#include <iostream>
using namespace std;

class Mesh{
	private:

		long double x_min;
		long double x_max;
		long double y_min;
		long double y_max;
		long double z_min;
		long double z_max;
		int ii; // Number of divisions in x axis
		int jj; // Number of divisions in y axis
		int kk; // Number of divisions in z axis
		long double dx;  // Uniform spatial step in x axis
		long double dy;  // Uniform spatial step in y axis
		long double dz;  // Uniform spatial step in z axis

	public:

		Mesh(const long double x_min0, const long double x_max0,
		 const long double y_min0, const long double y_max0,
		 const long double z_min0, const long double z_max0,
		 const int ii0, const int jj0, const int kk0)
		{
			x_min=x_min0;
			x_max=x_max0;
			y_min=y_min0;
			y_max=y_max0;
			z_min=z_min0;
			z_max=z_max0;
			ii=ii0;
			jj=jj0;
			kk=kk0;
			dx=(x_max-x_min)/(ii-1);
			dy=(y_max-y_min)/(jj-1);
			dz=(z_max-z_min)/(kk-1);
		}

		long double getXMax() const
		{
			return x_max;
		}
	    
		long double getYMax() const
		{
			return y_max;
		}
	    
		long double getZMax() const
		{
			return z_max;
		}
	    
		long double getXMin() const
		{
			return x_min;
		}
	    
		long double getYMin() const
		{
			return y_min;
		}
	    
		long double getZMin() const
		{
			return z_min;
		}
	    
		long double getDx() const
		{
			return dx;
		}
	    
		long double getDy() const
		{
			return dy;
		}

		long double getDz() const
		{
			return dz;
		}

		int getII() const
		{
			return ii;
		}

		int getJJ() const
		{
			return jj;
		}

		int getKK() const
		{
			return kk;
		}

		long double getX(const int i) const
		{
			return x_min+i*dx;
		}

		long double getY(const int j) const
		{
			return y_min+j*dy;
		}

		long double getZ(const int k) const
		{
			return z_min+k*dz;
		}

		~Mesh()
		{
		}

};

#endif
