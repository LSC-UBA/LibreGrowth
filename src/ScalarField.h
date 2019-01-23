#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


#include "Mesh.h"

using namespace std;

class ScalarField
{

	private:

		int ii;
		int jj;
		int kk;
		int n;
		long double * p_mat;
		const Mesh * mesh; // The mesh can be shared with other scalar field objects

	public:

		ScalarField(const Mesh &mesh0)
		{
			mesh = &mesh0;
			ii = mesh->getII();
			jj = mesh->getJJ();
			kk = mesh->getKK();
			n = ii * jj * kk;
			p_mat = new long double[n];
			init(0);
		}

		void init(long double val)
		{
			for (int l = 0; l < n; l++)
				p_mat[l] = val;
		}

		void set(const int i, const int j, const int k, long double val)
		{
			p_mat[i + j * ii + ii * jj * k] = val;
		}

		long double get(const int i, const int j, const int k) const
		{
			return p_mat[i + j * ii + ii * jj * k];

		}

		const Mesh& getMesh() const
		{
			return *mesh;
		}

		void print() const
		{
			for (int i = 0; i < ii; i++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int k = 0; k < kk; k++)
					{
						cout << mesh->getX(i)	<< "," << mesh->getY(j) << ","
									<< mesh->getZ(k) << ","
									<< p_mat[i + j * ii + ii * jj * k] << endl;
					}
				}
			}
		}

		void saveCsv(const ostringstream &filename) const
		{
			ofstream myfile1;
			long double aux;
			// Open file
			myfile1.open(filename.str().c_str());

			// Write content
			myfile1 << "X,Y,Z,S" << endl;
			for (int k = 0; k < kk; k++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int i = 0; i < ii; i++)
					{
						aux = p_mat[i + j * ii + ii * jj * k];
						myfile1 << mesh->getX(i) << "," << mesh->getY(j) << ","
							<< mesh->getZ(k) << "," << aux << endl;
					}
				}
			}
			// Close file name
			myfile1.close();

		}
		
		void saveVtk(const ostringstream &filename) const
		{

			ofstream myfile1;

			// Open file
			myfile1.open(filename.str().c_str());

			// Write content
			myfile1 << "# vtk DataFile Version 2.0\n";
			myfile1 << "Comment goes here\n";
			myfile1 << "ASCII\n";
			myfile1 << "\n";
			myfile1 << "DATASET STRUCTURED_POINTS\n";
			myfile1 << "DIMENSIONS  " << ii << " " << jj << " " << kk << endl;
			myfile1 << "\n";
			myfile1 << "ORIGIN	0.000   0.000   0.000\n";
			myfile1 << "SPACING	"<< mesh->getDx() <<" "<< mesh->getDy() << " " << mesh->getDz() << "\n";
			myfile1 << "\n";
			myfile1 << "POINT_DATA " << ii * jj * kk << endl;
			myfile1 << "SCALARS scalars float\n";
			myfile1 << "LOOKUP_TABLE default\n";
			myfile1 << "\n";
		
			for (int k = 0; k < kk; k++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int i = 0; i < ii; i++)
					{
						myfile1 << p_mat[i + j * ii + ii * jj * k] << " ";
					}
					myfile1 << endl;
				}
			}

			// Close file name
			myfile1.close();
		}

		// It is assumed that the mesh is correct
		void readVtk(const ostringstream & filename)
		{

			ifstream myfile1;
			string line;
			vector<std::string> result;
			istringstream * iss;

			// Open file
			myfile1.open(filename.str().c_str());

			// Read file and print
			if (myfile1.is_open())
			{

				// Drop 5 lines
				for (int i = 0; i < 5; i++)
					getline(myfile1,line);

				// Extract ii, jj, kk
				getline(myfile1,line);

				// Drop 2 lines
				for (int i = 0; i < 2; i++)
					getline(myfile1,line);

				// Extract origin: x_min, y_min, z_min
				getline(myfile1,line);

				// Extract dx, dy, dz
				getline(myfile1,line);

				// Drop 5 lines
				for (int i = 0; i < 4; i++)
					getline(myfile1,line);
				
				// Extract data
				for (int k = 0; k < kk; k++)
				{
					for (int j = 0; j < jj; j++)
					{

						getline(myfile1,line);

						iss = new istringstream(line);
										
						for (int i = 0; i < ii; i++)
						{
							string elem; *iss >> elem;
							p_mat[i + j * ii + ii * jj * k] = stod( elem );
						}

						delete iss;
					}
				}

			}
			else cout << "Unable to open file"; 

			// Close file name
			myfile1.close();

		}
		
		long double getNorm() const
		{
			long double max = p_mat[0];
			for (int l = 1; l < n; l++)
			{
				if (max < p_mat[l])
					max = p_mat[l];
			}
			return max;
		}

		// It assums that ii, jj and kk doesn't change
		ScalarField& operator =(const ScalarField & other)
		{
			mesh = &other.getMesh();
			ii = mesh->getII();
			jj = mesh->getJJ();
			kk = mesh->getKK();
			if (p_mat == NULL)
			{
				int n = ii * jj * kk;
				p_mat = new long double[n];
			}
			for (int i = 0; i < ii; i++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int k = 0; k < kk; k++)
					{
						p_mat[i + j * ii + ii * jj * k] = other.get(i, j, k);
					}
				}
			}
			return *this;
		}

		long double & operator()(const int i, const int j, const int k) const
		{
			return p_mat[i+j*ii+ii*jj*k];
		}

		ScalarField operator +(const ScalarField & other) const
		{
			ScalarField * add = new ScalarField(*mesh);
			for (int i = 0; i < ii; i++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int k = 0; k < kk; k++)
					{
						add->set(i, j, k,
							p_mat[i + j * ii + ii * jj * k]
							- other.get(i, j, k));
					}
				}
			}
			return *add;
		}

		ScalarField operator -(const ScalarField & other) const
		{
			ScalarField * diff = new ScalarField(*mesh);
			for (int i = 0; i < ii; i++)
			{
				for (int j = 0; j < jj; j++)
				{
					for (int k = 0; k < kk; k++)
					{
						diff->set(i, j, k,
							p_mat[i + j * ii + ii * jj * k]
							- other.get(i, j, k));
					}
				}
			}
			return *diff;
		}

		~ScalarField()
		{
			if (p_mat != NULL)
				delete[] p_mat;
			// The mesh can be shared with other scalar field objects
			// ergo I will not delete the Mesh
		}

};

