// Params.h
 
#ifndef SCALAR_FIELD
#define SCALAR_FIELD
#include "ScalarField.h"
#endif
 
#ifndef MESH
#define MESH
#include "Mesh.h"
#endif
 
/////////////////////////////////////////////
// Here you will find the variables  
// and parameters of the simulation
////////////////////////////////////////////
namespace Params
{
   
	/////////////////////////////////////////////
	// Biological (and mathematical) parameters
	////////////////////////////////////////////
   
		// Note: The following parameters are fitted with respect
		// to to the in-vitro experiment: 
		// D_inv_dim, vc_dim, p_coef_dim, s_dim, filter_dim

		const long double pi = 3.14159265359;

		// Dimensional parameter values
		const long double x_min_dim = -450.;
		const long double x_max_dim = 450.;
		const long double y_min_dim = -450;
		const long double y_max_dim = 450;
		const long double z_min_dim = -92;
		const long double z_max_dim = 92;
		const long double u_max_dim = 0.00073; // Max. concentration of tumoral cells. Units: cells/um^3.
		const long double u_core_max_dim = 0.0007; //  Max. concentration of core cells. Units: cells/um^3.
		const long double D_inv_dim = 70.;//1.7; // Invasion diffusion coefficient. Units: um^2/h
		const long double p_coef_dim = 0.03; // Proliferation coefficient. Units: 1/h
		const long double s_dim = 0.00009; // Source coefficient. Units: cell/(um^3 h)
		const long double t_total_dim = 24. * 15.; // Total simulation time (bening + invasion): 15 days. Units: h
		const long double t_max_benign_dim = 24. * 10.; // Max time of benign stage. 10 days. Units: h
		const long double t_max_inv_dim = t_total_dim - t_max_benign_dim + 0.01; // Max time after invasion begins: 5 days. Units: h
		const long double r_min_core_dim = 5.; // Initial core radius. Units: um
		const long double r_inv_dim = 45.08; // Critical invasion radio. Units: um
		const long double vc_dim = ( r_inv_dim - r_min_core_dim ) / t_max_benign_dim; // Core velocity. Units: um/h
		const long double r_core_dim = r_min_core_dim; // Core radius, changes over time.
		const long double curr_time_dim = 0.; // Units: h

		// Reference parameters (for non dimensional parameters)
		const long double u0 = 1.;
		const long double r0 = 1.;
		const long double t0 = 1.;
		const long double D0 = ( r0 * r0 ) / t0;

		// Non dimensional parameters
		const long double r_min_core = r_min_core_dim / r0;
		long double r_core = r_core_dim / r0;
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
		long double curr_time = curr_time_dim / t0;
		
		
	///////////////////////////////
	// Numerical parameters
	///////////////////////////////
 
		// Invasion mesh
		const int ii = 100; // Number of divisions in x-axis.
		const int jj = 100; // Number of divisions in y-axis.
		const int kk = 31; // Number of divisions in z-axis.
		
		int iterationNumber = 0;
		 
		const Mesh inv_mesh(x_min, x_max, y_min, y_max, z_min, z_max, ii, jj, kk);
		ScalarField u_new(inv_mesh); // Tumor Concentration in time N+1
		ScalarField u_old(inv_mesh); // Tumor Concentration in time N
		 
		ScalarField D(inv_mesh); // Diffusion coeficient
		  
		const long double dx = inv_mesh.getDx();
		const long double dy = inv_mesh.getDy();
		const long double dz = inv_mesh.getDz();
 
		const double dt_dim = 0.005 * 0.25; // Units: h
		const double dt = dt_dim / t0;  
		const int save_step = 6. / dt; // Save each 6 hours
	
		const long double filter_dim = 0; // Only save invasive concentrations greater that filter_dim. Units: cells/um^3.
		const long double filter = filter_dim / u0; // Non dimesional filter_dim
		 
		 
	///////////////////////////////
	// Auxiliar variables
	///////////////////////////////
	 
		long double t_max_core_growth = 0.; // Time needed for the core to get radii r_min_dim (it will be updated in the simulation, at the end of stage 1)
		long double curr_u,r,r2,vr,Dif,u,du_dr,d2u_dr2,dD_dr;
		long double curr_max = 0.;
		long double curr_min = 0.;	
		long double v_phi_aux[kk];
		long double aux;
		int i,j,k;
		int core_index;

		ostringstream filename("");
		ostringstream diffusionFile("../D.vtk");
	
	/////////////////////////////////////
	// Print parameters
	////////////////////////////////////
	
	void print_parameters() {

		cout << "DIMENSIONAL PARAMETERS:"<< endl;

		cout << "  r_min_core_dim:"<< r_min_core_dim << endl;
		cout << "  x_min_dim" << x_min_dim << endl;
		cout << "  x_max_dim" << x_min_dim << endl;
		cout << "  y_min_dim" << x_min_dim << endl;
		cout << "  y_max_dim" << x_min_dim << endl;
		cout << "  z_min_dim" << z_min_dim << endl;
		cout << "  z_max_dim" << z_max_dim << endl;
		cout << "  u_max_dim:"<< u_max_dim  << endl;
		cout << "  u_core_max_dim:"<< u_core_max_dim  << endl;
		cout << "  D_inv_dim:"<< D_inv_dim  << endl;
		cout << "  vc_dim:"<< vc_dim << endl;
		cout << "  p_coef_dim:"<<p_coef_dim << endl;
		cout << "  s_dim:"<< s_dim << endl;
		cout << "  t_max_inv_dim:"<< t_max_inv_dim  << endl;
		cout <<  endl;

		cout << "REFERENCE PARAMETERS FOR NON-DIMENSIONAL ANALYSIS:"<< endl;
		cout << "  u0:"<< u0 << endl;
		cout << "  r0:"<< r0 << endl;
		cout << "  t0:"<< t0  << endl;
		cout << "  D0:"<< D0  << endl;
		cout <<  endl;

		cout << "NON-DIMENSIONAL PARAMETERS:"<< endl;
		cout << "  r_min_core:"<< r_min_core << endl;
		cout << "  x_min" << x_min << endl;
		cout << "  x_max" << x_min << endl;
		cout << "  y_min" << x_min << endl;
		cout << "  y_max" << x_min << endl;
		cout << "  z_min" << z_min << endl;
		cout << "  z_max" << z_max << endl;
		cout << "  u_max:"<< u_max  << endl;
		cout << "  u_core_max:"<< u_core_max  << endl;
		cout << "  D_inv:"<< D_inv  << endl;
		cout << "  vc:"<< vc << endl;
		cout << "  p_coef:"<<p_coef << endl;
		cout << "  s:"<< s << endl;
		cout << "  t_max_inv:"<< t_max_inv  << endl;
		cout <<  endl;

		cout << "NUMERICAL PARAMETERS:"<< endl;
		cout << "  dt_dim:" << dt_dim << endl;
		cout << "  dx: "<< dx << ", dy: "<< dy << ", dz: "<< dz << endl;  
		cout << "  ii: " << ii << ", jj: "<< jj << ", kk: "<< kk << endl;
		cout << "  filter_dim:" << filter_dim << endl;
		cout << "  filter:" << filter << endl;
		cout << "  save_step:" << save_step << endl;
		
		cout <<  endl;
		cout <<  endl;
 
	}
	
 
}