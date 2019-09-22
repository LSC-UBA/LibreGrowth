#include "par.h"
#include "mesh.h"
#include "scalar_field.h"
#include <iomanip>

using namespace std;
using namespace par;

void set_diffusion_coefficient( ScalarField & D )
{
    
    double D0, D1;
    D1 = par::D_inv;
    D0 = 0.;
        
    int particlesNbr = 2000;
    double r;
    int i, j ,k ,in ,jn ,kn;
    srand (time(NULL));
        
    /* DLA */

    D.init(D0);
        
    int bottom = par::kk / 2 - 3;
    int top = par::kk / 2 + 3;
    
    for (int k = bottom; k <= top; k++) 
    {
        D(par::ii/2,par::jj/2,k) = D1;
    }

    int half;
    half = int ( par::kk / 2 );
    for(int p = 0; p < particlesNbr; p++)
    {
        i = rand() % par::ii;
        j = rand() % par::jj;
        k = half;
        
        bool attached = false;
                           
        /* While the particle is not attached */   
        while( !attached ) 
        {
            r = ( rand() % 1000 ) / 1000.;
            in = i;
            jn = j;
            kn = k;

            /* Move particle */
            
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
              
            /* Check boundaries */
            if ( in != -1 && in != par::ii  && jn != -1 && jn != par::jj && kn != bottom -1  && kn != top + 1)
            {
                i = in;
                j = jn;
                k = kn;
            }

            attached = ( ( i > 0 ) && D(i-1,j,k) != D0 ) ||
            ( ( i < par::ii - 1 ) && D(i+1,j,k) != D0 ) ||
            ( ( j > 0 ) && D(i,j-1,k) != D0 ) ||
            ( ( j < par::jj - 1 ) && D(i,j+1,k) != D0 ) ||
            ( ( k > 0 ) && D(i,j,k-1) != D0 ) ||
            ( ( k < par::kk - 1 ) && D(i,j,k+1) != D0 );
                        
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

    /* Save */
    D.save(par::filename_D,par::save_format);

}

void save(long double curr_time, ScalarField & u_new)
{
    /*Set tumor concetration file name and save it */

    std::ostringstream stream;
    stream  << "tumor_concentration_hs_"
            << std::setfill('0')
            << std::setw(10) << std::fixed 
            << std::setprecision(0) << curr_time;
    u_new.save(stream.str(),par::save_format);
}

/* Print iteration data */
void print(long double curr_time, long double iterationNumber, long double r_core)
{
    cout << "  Time: " << curr_time * par::t0 << ", core radio: "<< r_core * par::r0 <<
            ", iteration: " << iterationNumber << endl;
}

/* Print parameters */

void print_parameters() {

    cout << "DIMENSIONAL PARAMETERS:"<< endl;

    cout << "  r_min_core_dim:"<< par::r_min_core_dim << endl;
    cout << "  x_min_dim" << par::x_min_dim << endl;
    cout << "  x_max_dim" << par::x_min_dim << endl;
    cout << "  y_min_dim" << par::x_min_dim << endl;
    cout << "  y_max_dim" << par::x_min_dim << endl;
    cout << "  z_min_dim" << par::z_min_dim << endl;
    cout << "  z_max_dim" << par::z_max_dim << endl;
    cout << "  u_max_dim:"<< par::u_max_dim  << endl;
    cout << "  u_core_max_dim:"<< par::u_core_max_dim  << endl;
    cout << "  D_inv_dim:"<< par::D_inv_dim  << endl;
    cout << "  vc_dim:"<< par::vc_dim << endl;
    cout << "  p_coef_dim:"<< par::p_coef_dim << endl;
    cout << "  s_dim:"<< par::s_dim << endl;
    cout << "  t_max_inv_dim:"<< par::t_max_inv_dim  << endl;
    cout <<  endl;

    cout << "REFERENCE PARAMETERS FOR NON-DIMENSIONAL ANALYSIS:"<< endl;
    cout << "  u0:"<< par::u0 << endl;
    cout << "  r0:"<< par::r0 << endl;
    cout << "  t0:"<< par::t0  << endl;
    cout << "  D0:"<< par::D0  << endl;
    cout <<  endl;

    cout << "NON-DIMENSIONAL PARAMETERS:"<< endl;
    cout << "  r_min_core:"<< par::r_min_core << endl;
    cout << "  x_min" << par::x_min << endl;
    cout << "  x_max" << par::x_min << endl;
    cout << "  y_min" << par::x_min << endl;
    cout << "  y_max" << par::x_min << endl;
    cout << "  z_min" << par::z_min << endl;
    cout << "  z_max" << par::z_max << endl;
    cout << "  u_max:"<< par::u_max  << endl;
    cout << "  u_core_max:"<< par::u_core_max  << endl;
    cout << "  D_inv:"<< par::D_inv  << endl;
    cout << "  vc:"<< par::vc << endl;
    cout << "  p_coef:"<< par::p_coef << endl;
    cout << "  s:"<< par::s << endl;
    cout << "  t_max_inv:"<< par::t_max_inv  << endl;
    cout <<  endl;

    cout << "NUMERICAL PARAMETERS:"<< endl;
    cout << "  dt_dim:" << par::dt_dim << endl;
    cout << "  par::dx: "<< par::dx << ", par::dy: "<< par::dy << ", par::dz: "<< par::dz << endl;  
    cout << "  par::ii: " << par::ii << ", par::jj: "<< par::jj << ", par::kk: "<< par::kk << endl;
    cout << "  filter_dim:" << par::filter_dim << endl;
    cout << "  filter:" << par::filter << endl;
    cout << "  save_step:" << par::save_step << endl;

    cout <<  endl;
    cout <<  endl;

}


/* Simulation */
int main()
{

    /* ------------------------- Print parameters -------------------------- */
    print_parameters();

    /* -------------------- Declare dependent variables -------------------- */
    /* Tumor Concentration in time N+1 */
    ScalarField u_new(par::inv_mesh);
    /* Tumor Concentration in time N */
    ScalarField u_old(par::inv_mesh);
    /* Diffusion coeficient */
    ScalarField D(par::inv_mesh); 

    /* ----------------- Declare state and aux. variables --------------------*/
    /* Time needed for the core to get radpar::ii r_min_dim */
    /* (it will be updated in the simulation, at the end of stage 1) */
    long double aux;
    long double curr_time = 0.;
    long double r_core = r_core_dim / r0;
    long double x, y, z;
    int i,j;
    int iterationNumber = 0;


    /* ------------------ Initial and boundary conditions ------------------ */
    /* Set tumor concentration */
    u_old.init(0.);
    u_old(par::ii/2.,par::jj/2.,par::kk/2.) = u_max;
    u_new = u_old;
    /* Generate diffusion coefficient */
    set_diffusion_coefficient(D);
    /* Load diffusion coefficient */
    /*D.readVtk(diffusionFile); */


    /* ----------------------------- Simulation ---------------------------- */

    cout << "SIMULATION LOG:"<< endl;
    cout << "  Time (hs), core radio (um), time iteration" << endl;

    /* First stage or benign stage: only the core grows, 
    from radpar::ii r_min_core to radpar::ii r_min */

    /* Temporal iteration */
    while( curr_time <= par::t_max_benign )
    {

        /* Save & print */
        if( iterationNumber % par::save_step == 0 ) 
        {
            
                /* Update core */
                r_core = par::r_min_core + par::vc * curr_time;
                for (int k = 0; k < par::kk; k++) 
                {
                    for (int j = 0; j < par::jj; j++)
                    {
                        for (int i = 0; i < par::ii; i++)
                        {
                            x = par::x_min + par::dx * i;
                            y = par::y_min + par::dy * j;
                            z = par::z_min + par::dz * k;
                            long double r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
                            if ( r <= r_core )
                            {
                                u_new(i ,j , k) = par::u_max;
                            }
                        }
                    }
                }
                
                /* Save */
                save(curr_time,u_new);
                
                /* Print iteration data */
                print(curr_time, iterationNumber, r_core);

        }
            
        /* Increment time  */
        curr_time += dt;
        iterationNumber++;
    }
    
    u_old = u_new;

    /* Update core */
    r_core = par::r_min_core + par::vc * curr_time;
      
    cout << "Invasion begins"<< endl;
    
    /* Second stage or malignant stage: core and invasion concentration grows */


    /* Temporal iteration */
    while( curr_time <= ( par::t_max_benign + par::t_max_inv ) )
    {
        /* Solving linear equation system */
        
        #pragma omp parallel for collapse(3) private(aux,x,y,z) 
        for(int k = 1; k < par::kk-1; k++)
        {
            for(int j = 1; j < par::jj-1; j++)
            {
                for(int i = 1; i < par::ii-1; i++)
                {
                    x = x_min + par::dx * i;
                    y = y_min + par::dy * j;
                    z = z_min + par::dz * k;
                    if( pow(x,2.) + pow(y,2.) + pow(z,2.) > pow(r_core,2.) )
                    {

                        long double D_i_12 = ( D(i+1,j,k) + D(i,j,k) ) / 2.;
                        long double D_i_m12 = ( D(i,j,k) + D(i-1,j,k) ) / 2.;
                        long double du_i_12 = D_i_12  * (u_old(i+1,j,k) - u_old(i,j,k)) / par::dx;
                        long double du_i_m12 = D_i_m12  * (u_old(i,j,k) - u_old(i-1,j,k)) / par::dx;
                        long double d2u_i  = ( D_i_12 * du_i_12 - D_i_m12 * du_i_m12  ) / par::dx;

                        long double D_j_12 = ( D(i,j+1,k) + D(i,j,k) ) / 2.;
                        long double D_j_m12 = ( D(i,j,k) + D(i,j-1,k) ) / 2.;
                        long double du_j_12 = D_j_12  * (u_old(i,j+1,k) - u_old(i,j,k)) / par::dy;
                        long double du_j_m12 = D_j_m12  * (u_old(i,j,k) - u_old(i,j-1,k)) / par::dy;
                        long double d2u_j  = ( D_j_12 * du_j_12 - D_j_m12 * du_j_m12  ) / par::dy;

                        long double D_k_12 = ( D(i,j,k+1) + D(i,j,k) ) / 2.;
                        long double D_k_m12 = ( D(i,j,k) + D(i,j,k-1) ) / 2.;
                        long double du_k_12 = D_k_12  * (u_old(i,j,k+1) - u_old(i,j,k)) / par::dz;
                        long double du_k_m12 = D_k_m12  * (u_old(i,j,k) - u_old(i,j,k-1)) / par::dz;
                        long double d2u_k  = ( D_k_12 * du_k_12 - D_k_m12 * du_k_m12  ) / par::dz;

                        long double prolif_ijk = par::p_coef * u_old(i,j,k) * (1. - u_old(i,j,k) / par::u_max );

                        aux = ( ( d2u_i + d2u_j + d2u_k ) + prolif_ijk ) * par::dt + u_old(i,j,k);

                        u_new(i,j,k) = aux;

                        if(isnan(aux))
                        {
                                cout << "  Time: " << curr_time << ", core radio: "<< r_core*r0 <<
                                        ", iteration: " << iterationNumber << ", i:" << i << 
                                        ", j:" << j << ", k:" << k <<  ", Curr.conc.:" << aux << endl;
                                exit(EXIT_FAILURE);
                        }

                    }

                }
            }
        }

        /* Main computation: neumann boundary conditions */
        i = par::ii - 1;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < par::jj; j++)
        {
            for(int k = 0; k < par::kk; k++)
            {
                u_new(i,j,k) = u_new(i-1,j,k);
            }
        }

        i = 0;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < par::jj; j++)
        {
            for(int k = 0; k < par::kk; k++)
            {
                u_new(i,j,k) = u_new(i+1,j,k);
            }
        }

        j= par::jj - 1;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < par::ii; i++)
        {
            for(int k = 0; k < par::kk; k++)
            {
                u_new(i,j,k) = u_new(i,j-1,k);
            }
        }

        j = 0;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < par::ii; i++)
        {
            for(int k = 0; k < par::kk; k++)
            {
                u_new(i,j,k) = u_new(i,j+1,k);
            }
        }

        /* Source term. Using Lie-Trotter splitting method for Dirac Delta term */
        #pragma omp parallel for collapse(3) private(x,y,z) schedule(static)
        for(int k = 0; k < par::kk; k++)
        {
            for(int j = 0; j < par::jj; j++)
            {
                for(int i = 0; i < par::ii; i++)
                {
                    x = par::x_min + par::dx * i;
                    y = par::y_min + par::dy * j;
                    z = par::z_min + par::dz * k;
                    long double dr = sqrt( pow(par::dx,2.) + pow(par::dy,2.)+pow(par::dz,2.) );
                    long double r = sqrt( pow(x,2.) + pow(y,2.)+pow(z,2.) );
                    /* r == r_core */
                    if( fabs( r - r_core ) < dr / 2. ) 
                    {
                        u_new(i,j,k) = s / par::vc;
                    }
                    else if ( r < r_core )
                    {
                        u_new(i,j,k) = par::u_max;
                    }
                }
            }
        }

        /* Save & print */
        if(iterationNumber % save_step == 0)
        {
            /* Save */
            save(curr_time,u_new);

            /* Print iteration data */
            print(curr_time, iterationNumber, r_core);
        }
            
        /* Update invasion concentration */
        u_old = u_new;
            
        /* Update core */
        r_core = par::r_min_core + par::vc * curr_time;
            
        /* Increment time */
        curr_time += dt;
        iterationNumber++;

    }

    return 0;
}
