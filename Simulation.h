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
		void buildUpB;
		float pressurePreSolve(int xLocation, int yLocation);
		float xMomentumSolve(int xLocation, int yLocation);
		float yMomentumSolve(int xLocation, int yLocation);
		void iterate(void);

	private:
		//Size values
		int nx = 101;
		int ny = 101;
		int nt = 100;
		int nit = nt;
		int dx = 2/(nx-1);
		int dy = 2/(nx-1);
		int dt = 0.0001;

		//Physical Values
		int rho;
		int nu; 
		int F; //needed?

		int numCells;
		int leftOvers;
		int startingLocation;

		//Solved data storage
		datumPoint localVelData[numCells];
		datumPoint solvedVelData[nt][problemSize];
		double localPrePresData[numCells];
		double solvedPrePresData[nt][problemSize];
//		double bPressureData[numCells];

		//MPI data
		int myRank, size;
		int counter = 1;
		int subCounter = 1;
		int recCounts[size];
		int displs[size];

};
