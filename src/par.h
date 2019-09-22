/* Params.h */

#ifndef PAR_H_
#define PAR_H_

#include "scalar_field.h"
#include "mesh.h"

/* --- Here you will find the variables and parameters of the simulation ----*/

namespace par
{

/* --------------- Biological and mathematical parameters ------------------ */

/* Note: The following parameters are fitted with respect to to the in-vitro
   experiment: D_inv_dim, vc_dim, p_coef_dim, s_dim, filter_dim */

/* PI */
const long double pi = 3.14159265359;
/* Domain limits. Units: um */
const long double x_min_dim = -450.;
const long double x_max_dim = 450.;
const long double y_min_dim = -450;
const long double y_max_dim = 450;
const long double z_min_dim = -92;
const long double z_max_dim = 92;
/* Max. concentration of tumoral cells. Units: cells/um^3. */
const long double u_max_dim = 0.00073;
/* Max. concentration of core cells. Units: cells/um^3. */
const long double u_core_max_dim = 0.0007; 
/* Invasion diffusion coefficient. Units: um^2/h */
const long double D_inv_dim = 70.; 
/* Proliferation coefficient. Units: 1/h */
const long double p_coef_dim = 0.03;
/* Source coefficient. Units: cell/(um^3 h) */
const long double s_dim = 0.00009; 
/* Total simulation time (bening + invasion): 15 days. Units: h */
const long double t_total_dim = 24. * 15.; 
/* Max time of benign stage. 10 days. Units: h */
const long double t_max_benign_dim = 24. * 10.; 
/* Max time after invasion begins: 5 days. Units: h */
const long double t_max_inv_dim = t_total_dim - t_max_benign_dim + 0.01; 
/* Initial core radius. Units: um */
const long double r_min_core_dim = 5.; 
/* Critical invasion radio. Units: um */
const long double r_inv_dim = 45.08; 
/* Core velocity. Units: um/h */
const long double vc_dim = ( r_inv_dim - r_min_core_dim ) / t_max_benign_dim; 
/* Core radius, changes over time. */
const long double r_core_dim = r_min_core_dim; 


/* ----------------------- Nondimensionalization --------------------------- */

/* Reference parameters (for non dimensional parameters) */
const long double u0 = 1.;
const long double r0 = 1.;
const long double t0 = 1.;
const long double D0 = ( r0 * r0 ) / t0;

/* Non dimensional parameters */
const long double r_min_core = r_min_core_dim / r0;
const long double x_min = x_min_dim / r0;
const long double x_max = x_max_dim / r0;
const long double y_min = y_min_dim / r0;
const long double y_max = y_max_dim / r0;
const long double z_min = z_min_dim / r0;
const long double z_max = z_max_dim / r0;
const long double u_max = u_max_dim / u0;
const long double u_core_max = u_core_max_dim / u0;
const long double D_inv = D_inv_dim / D0;
const long double vc = vc_dim * t0 / r0;
const long double p_coef = p_coef_dim * t0;
const long double s = s_dim * t0 / u0;
const long double t_max_benign = t_max_benign_dim / t0;
const long double t_max_inv = t_max_inv_dim / t0;

/* ------------------------ Numerical parameters ------------------------- */
 
/* Invasion mesh */
/* Number of divisions in x-axis. */
const int ii = 100;
/* Number of divisions in y-axis. */
const int jj = 100;
/* Number of divisions in z-axis. */
const int kk = 31; 
/* Mesh */ 
const Mesh inv_mesh(x_min, x_max, y_min, y_max, z_min, z_max, ii, jj, kk);
/* Distance increment. Units: um */
const long double dx = inv_mesh.get_dx();
const long double dy = inv_mesh.get_dy();
const long double dz = inv_mesh.get_dz();
/* Time increment. Units: h */
const double dt_dim = 0.005 * 0.25; 
const double dt = dt_dim / t0; 
/* Save every 6 hours */
const int save_step = 6. / dt; 
/* Only save invasive concentrations greater that filter_dim. Units: cells/um^3. */
const long double filter_dim = 0;
/* Non dimesional filter_dim */
const long double filter = filter_dim / u0;
/* File names */
const std::string filename_D = "D";
/* Save format:  "vtk" or "csv" */
const std::string save_format = "vtk";

}

#endif
