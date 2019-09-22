#include "mesh.h"

/* Mesh class keeps information about geometry and domain discretization */

Mesh::Mesh(const long double x_min0, const long double x_max0,
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

long double Mesh::get_xmax() const
{
    return x_max;
}

long double Mesh::get_ymax() const
{
    return y_max;
}

long double Mesh::get_zmax() const
{
    return z_max;
}

long double Mesh::get_xmin() const
{
    return x_min;
}

long double Mesh::get_ymin() const
{
    return y_min;
}

long double Mesh::get_zmin() const
{
    return z_min;
}

long double Mesh::get_dx() const
{
    return dx;
}

long double Mesh::get_dy() const
{
    return dy;
}

long double Mesh::get_dz() const
{
    return dz;
}

int Mesh::get_ii() const
{
    return ii;
}

int Mesh::get_jj() const
{
    return jj;
}

int Mesh::get_kk() const
{
    return kk;
}

long double Mesh::get_x(const int i) const
{
    return x_min+i*dx;
}

long double Mesh::get_y(const int j) const
{
    return y_min+j*dy;
}

long double Mesh::get_z(const int k) const
{
    return z_min+k*dz;
}

Mesh::~Mesh()
{
}


