// Import standard libraries
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>
#include <dirent.h>
#include <string> // For string operations
#include <sstream> // For converting int to string
#include <sys/types.h> // For mkdir
#include <sys/stat.h>   // For mkdir
#include <iomanip> // For setprecision

#include "param.c" // Include parameters from an external file

using namespace std;

// Constants and initial values
clock_t start = clock(); // Start the clock to measure execution time
constexpr double radius = 1; // Particle radius

// Parameters from the included param.c file
int N1 = par_N;
double phi = par_phi;
double Dr = par_Dr;
double v0 = par_Pe * Dr * radius;
double k1 = par_K1;
double k2 = par_K2;
int NN1 = par_N1;
int NN2 = par_N2;
double k;

// Scaled radius
constexpr double pi = 22.000 / 7;
constexpr double TAU = 1;
constexpr double dt = 0.001 * TAU; // Time step
double MAXIT = par_maxIT * TAU / dt; // Maximum iterations

// Calculate parameters for cell list
double rcut = 2 * radius; // Cutoff radius
int n = sqrt(N1); // Number of particles in one dimension
int N = n * n; // Total number of particles
double L = sqrt(N * pi * radius * radius / phi); // System size
int M = int(L / rcut); // Number of cells in one dimension
double alat = 2 * L / (2 * n + 3); // Lattice spacing

double BOXL = M * rcut; // Box length
int ncell = M * M; // Total number of cells
int mapsiz = 4 * ncell; // Size of the map
double A = BOXL * BOXL; // Area of the box

// Constants for diffusion
double DrDtConstat = sqrt(2 * Dr * dt); 

// Vectors for storing particle properties
vector<int> list_(N);
vector<double> Fx(N), Fy(N);
vector<int> head(ncell + 1);
vector<int> map_(mapsiz + 2);
int ix, iy, iz, imap;

vector<double> theta(N);
vector<double> x(N), y(N);

vector<double> vx(N);
vector<double> vy(N);
vector<double> spring_const(N);
vector<int> mi(N);

// Random number generators
mt19937 generator(time(NULL));
normal_distribution<double> gau_dist(0.0, 1.0);
uniform_real_distribution<double> uni_dist(-pi / 4, pi / 4);

// Function prototypes
void initialize_square_lattice();
void identify_particles();
int icell(int ix, int iy);
void maps();
void links(int it);
void force();
void update_position_angle();
void write_VMD_data(ostream &os);

// Main function
int main(int argc, char *argv[])
{
    // Initialization and output file streams
    ifstream list;
    ofstream o("vmd_data.xyz");
    ofstream out("data.xyz");

    // Write initial parameters to output file
    out << "phi: " << phi << endl;
    out << "N: " << N << endl;
    out << "System Size: " << BOXL << endl;
    out << "v0: " << v0 << endl;
    out << "Pe: " << par_Pe << endl;
    out << "k: " << k1 << " k: " << k2 << endl;

    int frame = par_maxFrame;
    int tn = MAXIT / frame;

    // Initialize lattice and particles
    initialize_square_lattice();
    identify_particles();
    maps();

    // Time-stepping loop
    for (int it = 0; it < MAXIT; it++)
    {
        links(it);
        force();
        if (it % tn == 0)
            write_VMD_data(o);
        update_position_angle();
    }

    // Write final parameters to output file
    out << "k: " << k << endl;
    out << "Time elapsed: " << ((double(clock()) - start) / CLOCKS_PER_SEC) << endl;
}

/******************* Sub functions *******************/

// Function to initialize square lattice
void initialize_square_lattice()
{
    int nptl = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j % 2 == 0)
                x[nptl] = -L / 2.0 + i * alat;
            else
                x[nptl] = -L / 2.0 + (2.0 * i + 1.0) * alat / 2.0;
            if (i == 0)
                y[nptl] = -L / 2.0 + (2.0 * j + 1.0) * alat / 2.0;
            else
                y[nptl] = -L / 2.0 + j * alat;

            theta[nptl] = 0;
            nptl++;
        }
    }
}

// Function to identify particles based on some criteria
void identify_particles()
{
    for (int i = 0; i < N; i++)
    {
        if (i >= int(NN2 * N / 100))
        {
            mi[i] = 1;
            spring_const[i] = k1;
        }
        else
        {
            mi[i] = 2;
            spring_const[i] = k2;
        }
    }
}

// Function to calculate cell index based on x, y positions
int icell(int ix, int iy)
{
    int a = 1 + (ix - 1 + M) % M + ((iy - 1 + M) % M) * M;
    return a;
}

// Function to create maps for cells
void maps()
{
    for (int iy = 1; iy <= M; iy++)
    {
        for (int ix = 1; ix <= M; ix++)
        {
            imap = (icell(ix, iy) - 1) * 4;
            map_[imap + 1] = icell(ix + 1, iy);
            map_[imap + 2] = icell(ix + 1, iy + 1);
            map_[imap + 3] = icell(ix, iy + 1);
            map_[imap + 4] = icell(ix - 1, iy + 1);
        }
    }
}

// Function to create links between cells
void links(int it)
{
    // Zero head of chain array
    for (int icell = 1; icell <= ncell; icell++)
    {
        head[icell] = -1;
    }

    double celli = M;
    double cell = BOXL / celli;
    if (cell < rcut)
        cout << "Cell size too small for cutoff" << endl;

    // Sort all atoms
    for (int i = 0; i < N; i++)
    {
        int icel = 1 + int((x[i] / BOXL + 0.5) * celli) + int((y[i] / BOXL + 0.5) * celli) * M;
        list_[i] = head[icel];
        head[icel] = i;
        if (icel > M * M)
            cout << "Error in icel at iteration " << it << endl;
    }
}

// Function to calculate forces
void force()
{
    double rcutsq = rcut * rcut;
    double xi, yi, Fxi, Fyi;
    double xij, yij, rijsq, Fxij, Fyij;
    double rij;

    // Zero forces
    for (int i = 0; i < N; i++)
    {
        Fx[i] = 0.0;
        Fy[i] = 0.0;
    }

    // Loop over all cells
    for (int icell = 1; icell <= ncell; icell++)
    {
        int i = head[icell];
        // Loop over all molecules in the cell
        while (i > -1)
        {
            xi = x[i];
            yi = y[i];
            Fxi = Fx[i];
            Fyi = Fy[i];

            // Loop over all molecules below i in the current cell
            int j = list_[i];
            while (j > -1)
            {
                xij = x[i] - x[j];
                yij = y[i] - y[j];

                xij = xij - BOXL * rint(xij / BOXL);
                yij = yij - BOXL * rint(yij / BOXL);
                rijsq = xij * xij + yij * yij;
                rij = sqrt(rijsq);

                if (rijsq < rcutsq)
                {
                    k = spring_const[i] * spring_const[j] / (spring_const[i] + spring_const[j]);
                    Fxij = k * (2 * radius - rij) * xij / rij;
                    Fyij = k * (2 * radius - rij) * yij / rij;
                    Fxi += Fxij;
                    Fyi += Fyij;
					Fx[j] -= Fxij;
					Fy[j] -= Fyij;
				}
				j = list_[j];
			}
			int jcell0, jcell;
			jcell0 = 4 * (icell - 1);
			for (int nabor = 1; nabor <= 4; nabor++)
			{
				jcell = map_[jcell0 + nabor];
				// Loop over all molecules in neighbouring cells
				j = head[jcell];
				while (j > -1)
				{
					xij = x[i] - x[j];
					yij = y[i] - y[j];

					xij = xij - BOXL * rint(xij / BOXL);
					yij = yij - BOXL * rint(yij / BOXL);
					rijsq = xij * xij + yij * yij;
					rij = sqrt(rijsq);

					if (rijsq < rcutsq)
					{
						k = spring_const[i] * spring_const[j] / (spring_const[i] + spring_const[j]);
						Fxij = k * (2 * radius - rij) * xij / rij;
						Fyij = k * (2 * radius - rij) * yij / rij;
						Fxi += Fxij;
						Fyi += Fyij;
						Fx[j] -= Fxij;
						Fy[j] -= Fyij;
					}
					j = list_[j];
				}
			}
			Fx[i] = Fxi;
			Fy[i] = Fyi;
			i = list_[i];
		}
	}
}

// Function to update particle positions and angles
void update_position_angle()
{
	for (int i = 0; i < N; i++)
	{
		vx[i] = (Fx[i] + v0 * cos(theta[i]));
		vy[i] = (Fy[i] + v0 * sin(theta[i]));
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
		x[i] = x[i] - BOXL * rint(x[i] / BOXL);
		y[i] = y[i] - BOXL * rint(y[i] / BOXL);
		theta[i] += DrDtConstat * gau_dist(generator);
	}
}

// Function to write data in VMD compatible format
void write_VMD_data(ostream &os)
{
	os << N << endl;
	os << "LM"<< endl;
	os.setf(ios::fixed, ios::floatfield);
	os.precision(5);
	for (int i = 0; i < N; i++)
	{
		os << int(mi[i]) << " " << x[i] << " " << y[i] << " " << vx[i] << " " << vy[i] << " " << theta[i] << endl;
	}
	os.flush();
}