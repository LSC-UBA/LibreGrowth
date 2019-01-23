#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>	
#include <unistd.h>

#include <assert.h>

#ifndef SCALAR_FIELD
#define SCALAR_FIELD
#include "ScalarField.h"
#endif

#ifndef MESH
#define MESH
#include "Mesh.h"
#endif

#include "Params.h"

using namespace std;
using namespace Params;

void setDiffusionCoefficient( ScalarField & D )
{
	
	double D0, D1;
	D1 = D_inv;
	D0 = 0.;
		
	int particlesNbr = 2000;
	double r;
	int i, j ,k ,in ,jn ,kn;
	srand (time(NULL));
		
	//DLA

	D.init(D0);
		
	int bottom = kk / 2 - 3;
	int top = kk / 2 + 3;
	
	for (int k = bottom; k <= top; k++) 
	{
		D(ii/2,jj/2,k) = D1;
	}

	int half;
	half = int ( kk / 2 );
	for(int p = 0; p < particlesNbr; p++)
	{
		i = rand() % ii;
		j = rand() % jj;
		k = half;
		
		bool attached = false;
							  
		while( !attached ) // While the particle is not attached
		{
			r = ( rand() % 1000 ) / 1000.;
			in = i;
			jn = j;
			kn = k;

			// Move particle
			
			double prob = 1. / 6.;
			if ( r < prob )
			{
				in = i - 1;
			}
			else if ( r < 2 * prob )
			{
				kn = k - 1;
			}
			else if ( r < 3 * prob )
			{
				kn = k + 1;
			}
			else if ( r < 4 * prob )
			{
				jn = j + 1;
			}
			else if ( r < 5 * prob )
			{
				jn = j - 1;
			}
			else
			{
				in = i + 1;
			}
			  
			// Check boundaries
			if ( in != -1 && in != ii  && jn != -1 && jn != jj && kn != bottom -1  && kn != top + 1)
			{
				i = in;
				j = jn;
				k = kn;
			}

			attached = ( ( i > 0 ) && D(i-1,j,k) != D0 ) ||
			( ( i < ii - 1 ) && D(i+1,j,k) != D0 ) ||
			( ( j > 0 ) && D(i,j-1,k) != D0 ) ||
			( ( j < jj - 1 ) && D(i,j+1,k) != D0 ) ||
			( ( k > 0 ) && D(i,j,k-1) != D0 ) ||
			( ( k < kk - 1 ) && D(i,j,k+1) != D0 );
						
		}
		 
		r = ( rand() % 1000  ) / 1000.;
		if( r < 1./6. )
			 D(i+1,j,k) = D1;
		else if( r < 2./6. )
			 D(i-1,j,k) = D1;
		else if( r < 3./6. )
			 D(i,j+1,k) = D1;
		else if( r < 4./6. )
			 D(i,j-1,k) = D1;
		else if( r < 5./6. )
			 D(i,j,k+1) = D1;
		else
			 D(i,j,k-1) = D1;

	}

	// Save
	std::ostringstream os;
	os << "D.vtk";
	D.saveVtk(os);

}

void save()
{
	// Set tumor concetration file name and save it
	filename.clear();
	filename.str("");
	filename << "tumor_concentration_" << curr_time << "_hs.vtk";
	u_new.saveVtk(filename); 

	// filename.clear();
	// filename.str("");
	// filename << "tumor_concentration_" << curr_time << "_hs.csv";
	// u_new.saveCsv(filename);
}

// Print iteration data
void print()
{
	cout << "  Time: " << curr_time * t0 << ", core radio: "<< r_core * r0 << ", iteration: " << iterationNumber << endl;
}

// Simulation
int main()
{
	/////////////////////////////////////
	// Parameters: See Params.h
	////////////////////////////////////
	
	print_parameters();

	/////////////////////////////////////
	// Initial and boundary conditions
	////////////////////////////////////
	 
	// Set tumor concentration
	u_old.init(0.);
	u_old(ii/2.,jj/2.,kk/2.) = u_max;
	u_new = u_old;
		
	double x, y, z;

	///////////////////////////////
	// Diffusion coefficient
	///////////////////////////////

	// Generate diffusion coefficient
	setDiffusionCoefficient(D);

	// Load diffusion coefficient
	// D.readVtk(diffusionFile);

	///////////////////////////////
	// Simulation
	///////////////////////////////
		
	cout << "SIMULATION LOG:"<< endl;
	cout << "  Time (hs), core radio (um), time iteration" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// First stage or benign stage: only the core grows, from radii r_min_core to radii r_min
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Temporal iteration
	while( curr_time <= t_max_benign )
	{

		// Save & print
		if( iterationNumber % save_step == 0 ) 
		{
			
				// Update core
				r_core = r_min_core + vc * curr_time;
				for (int k = 0; k < kk; k++) 
				{
					for (int j = 0; j < jj; j++)
					{
						for (int i = 0; i < ii; i++)
						{
							x = x_min + dx * i;
							y = y_min + dy * j;
							z = z_min + dz * k;
							long double r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
							if ( r <= r_core )
							{
								u_new(i ,j , k) = u_max;
							}
						}
					}
				}
				
				// Save
				save();
				
				// Print iteration data
				print();
				
		}
			
		// Increment time
		curr_time += dt;
		iterationNumber++;
	}
	
	u_old = u_new;

	// Update core
	r_core = r_min_core + vc * curr_time;
	  
	cout << "Invasion begins"<< endl;
	
	/////////////////////////////////////////////////////////////////////////////////////
	// Second stage or malignant stage: core and invasion concentration grows
	/////////////////////////////////////////////////////////////////////////////////////
	
	// Temporal iteration
	while( curr_time <= ( t_max_benign + t_max_inv ) )
	{
		// Solving linear equation system
		
				#pragma omp parallel for collapse(3) private(aux,x,y,z) 
				for(int k = 1; k < kk - 1; k++)
				{
						for(int j = 1; j < jj-1; j++)
						{
								for(int i = 1; i < ii-1; i++)
								{
										x = x_min + dx * i;
										y = y_min + dy * j;
										z = z_min + dz * k;
										if( pow(x,2.) + pow(y,2.) + pow(z,2.) > pow(r_core,2.) )
										{

												long  double D_i_12 = ( D(i+1,j,k) + D(i,j,k) ) / 2.;
												long  double D_i_m12 = ( D(i,j,k) + D(i-1,j,k) ) / 2.;
												long  double du_i_12 = D_i_12  * (u_old(i+1,j,k) - u_old(i,j,k)) / dx;
												long  double du_i_m12 = D_i_m12  * (u_old(i,j,k) - u_old(i-1,j,k)) / dx;
												long  double d2u_i  = ( D_i_12 * du_i_12 - D_i_m12 * du_i_m12  ) / dx;

												long  double D_j_12 = ( D(i,j+1,k) + D(i,j,k) ) / 2.;
												long  double D_j_m12 = ( D(i,j,k) + D(i,j-1,k) ) / 2.;
												long  double du_j_12 = D_j_12  * (u_old(i,j+1,k) - u_old(i,j,k)) / dy;
												long  double du_j_m12 = D_j_m12  * (u_old(i,j,k) - u_old(i,j-1,k)) / dy;
												long  double d2u_j  = ( D_j_12 * du_j_12 - D_j_m12 * du_j_m12  ) / dy;

												long  double D_k_12 = ( D(i,j,k+1) + D(i,j,k) ) / 2.;
												long  double D_k_m12 = ( D(i,j,k) + D(i,j,k-1) ) / 2.;
												long  double du_k_12 = D_k_12  * (u_old(i,j,k+1) - u_old(i,j,k)) / dz;
												long  double du_k_m12 = D_k_m12  * (u_old(i,j,k) - u_old(i,j,k-1)) / dz;
												long  double d2u_k  = ( D_k_12 * du_k_12 - D_k_m12 * du_k_m12  ) / dz;

												long double prolif_ijk = p_coef * u_old(i,j,k) * (1. - u_old(i,j,k) / u_max );

												aux = ( ( d2u_i + d2u_j + d2u_k ) + prolif_ijk ) * dt + u_old(i,j,k);

												u_new(i,j,k) = aux;

												if(isnan(aux))
												{
														cout << "  Time: " << curr_time << ", core radio: "<< r_core*r0 << ", iteration: " << iterationNumber << ", i:" << i << ", j:" << j << ", k:" << k <<  ", Curr.conc.:" << aux << endl;
														exit(EXIT_FAILURE);
												}	

										}

								}
						}
				}
				
					
				// Main computation: neumann boundary conditions
				i = ii - 1;
				#pragma omp parallel for collapse(2) schedule(static)
				for(int j = 0; j < jj; j++)
				{
						for(int k = 0; k < kk; k++)
						{
								u_new(i,j,k) = u_new(i-1,j,k);
						}
				}
						
				i = 0;
				#pragma omp parallel for collapse(2) schedule(static)
				for(int j = 0; j < jj; j++)
				{
						for(int k = 0; k < kk; k++)
						{
								u_new(i,j,k) = u_new(i+1,j,k);
						}
				}
						
				j= jj - 1;
				#pragma omp parallel for collapse(2) schedule(static)
				for(int i = 1; i < ii; i++)
				{
						for(int k = 0; k < kk; k++)
						{
								u_new(i,j,k) = u_new(i,j-1,k);
						}
				}
						
				j = 0;
				#pragma omp parallel for collapse(2) schedule(static)
				for(int i = 1; i < ii; i++)
				{
						for(int k = 0; k < kk; k++)
						{
								u_new(i,j,k) = u_new(i,j+1,k);
						}
				}

				// Source term. Using Lie-Trotter splitting method for Dirac Delta term
				#pragma omp parallel for collapse(3) private(x,y,z) schedule(static)
				for(int k = 0; k < kk; k++)
				{
						for(int j = 0; j < jj; j++)
						{
								for(int i = 0; i < ii; i++)
								{
										x = x_min + dx * i;
										y = y_min + dy * j;
										z = z_min + dz * k;
										long double dr = sqrt( pow(dx,2.) + pow(dy,2.)+pow(dz,2.) );
										long double r = sqrt( pow(x,2.) + pow(y,2.)+pow(z,2.) );
										if( fabs( r - r_core ) < dr / 2. ) // r == r_core
										{
											u_new(i,j,k) = s / vc;
										}
										else if ( r < r_core )
										{
											u_new(i,j,k) = u_max;
										}
								}	
						}
				}
	
		// Save & print
		if(iterationNumber % save_step == 0)
		{
			// Save
			save();

			// Print iteration data
			print();
		}
			
		// Update invasion concentration
		u_old = u_new;
			
		// Update core
		r_core = r_min_core + vc * curr_time;
			
		// Increment time
		curr_time += dt;
		iterationNumber++;
			
	}
	   
	return 0;
}
