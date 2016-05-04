#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

//Pressure function can't be solved independently - how solve?

typedef struct
{
	double p;
	double u;
	double v;
} datumPoint;

class Simulation
{
	public:
		Simulation();
		void buildUpB;
		float pressureSolve(int xLocation, int yLocation);
		float xMomentumSolve(int xLocation, int yLocation);
		float yMomentumSolve(int xLocation, int yLocation);
		void iterate(void);

	private:
		//Size values
		int nx;
		int ny;
		int nt;
		int dx;
		int dy;
		int dt;

		//Physical Values
		int rho;
		int mu; 
		int F; //needed?

		int numCells;
		int leftOvers;
		int startingLocation;

		//Solved data storage
/*		datumPoint nGlobalSolvedData[problemSize];
		datumPoint localData[numCells];
		datumPoint allSolvedData[nt][problemSize];
*/
		datumPoint* nGlobalSolvedData;
		datumPoint* localData;
		datumPoint* allSolvedData;
		//MPI data
		int myRank, size;
		int counter;

};
