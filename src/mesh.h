#ifndef MESH_H_
#define MESH_H_
#include <math.h> 
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

class Mesh
{

    private:

        /* Domain limits */
        long double x_min;
        long double x_max;
        long double y_min;
        long double y_max;
        long double z_min;
        long double z_max;
        /* Number of divisions in x, y, and z axis */
        int ii;
        int jj;
        int kk;
        /* Uniform spatial step in x, y, and z axis */
        long double dx;
        long double dy;
        long double dz;

    public:

        Mesh(const long double x_min0, const long double x_max0,
             const long double y_min0, const long double y_max0,
             const long double z_min0, const long double z_max0,
             const int ii0, const int jj0, const int kk0);

        long double get_xmax() const;

        long double get_ymax() const;

        long double get_zmax() const;

        long double get_xmin() const;

        long double get_ymin() const;

        long double get_zmin() const;

        long double get_dx() const;

        long double get_dy() const;

        long double get_dz() const;

        int get_ii() const;

        int get_jj() const;

        int get_kk() const;

        long double get_x(const int i) const;

        long double get_y(const int j) const;

        long double get_z(const int k) const;

        ~Mesh();

};

#endif
