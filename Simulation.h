#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

//Pressure function can't be solved independently - how solve?

typedef struct
{
	double u;
	double v;
	double p;
} datumPoint;

class Simulation
{
	public:
		Simulation();
		double buildUpB(int xLocation, int yLocation);
		double pressurePreSolve(int xLocation, int yLocation);
		double pressureSolve(int xLocation, int yLocation);
		double xMomentumSolve(int xLocation, int yLocation);
		double yMomentumSolve(int xLocation, int yLocation);
		void iterate(void);

	private:
		//Size values
		const static int nx = 101;
		const static int ny = 101;
		const static int nt = 100;
		const static int nit = nt;
		double dx = 2/(nx-1);
		double dy = 2/(nx-1);
		double dt = 0.0001;

		//Physical Values
		int rho;
		int nu; 
		int F; //needed?

		int numCells;
		const static int problemSize = nx*ny*nt;
		int leftOvers;
		int startingLocation;
		int xLocation, yLocation;

		//Solved data storage
		datumPoint localVelData[];
		datumPoint solvedVelData[nt][problemSize];
		double localPrePresData[];
		double solvedPrePresData[nt][problemSize];
//		double bPressureData[numCells];

		//MPI data
		int myRank, size;
		int counter = 1;
		int subCounter = 1;

};
