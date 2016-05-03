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

void buildUpB;

void pressureSolve;

void xMomentumSolve;

void yMomentumSolve;

int main(int argc, char* argv[])
{
	int myRank, size, i, j;

	int rho = 
	int mu = 
	int F = //needed? 

	int nx = 101;
	int ny = 101;
	int nt = 100;
	int dx = 2/(nx-1);
	int dy = 2/(ny-1);
	int dt = 0.0001;
	int counter = 0;

	int* pressureGradient = (int*) malloc( sizeof(double) * nx * ny );
	int* xVelocity =  (int*) malloc( sizeof(double) * nx * ny );
	int* yVelocity =  (int*) malloc( sizeof(double) * nx * ny );
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists
//	int* pointMatrix =  (int*) malloc( sizeof(int) * nx * ny );


	//Various Vars needed -- init in MPI_init?
	int xLocation, yLocation, startingLocation;

	int problemSize = nx*ny;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//Workload distribution
	int numCells = problemSize/size;
	int leftOvers = problemSize % size;
	for (i=1; i<leftOvers+1; i++)
		if (i == myRank)
			numCells += 1;

	startingLocation = numCells * rank;

	datumPoint globalSolvedData[problemSize];
	datumPoint localData[numCells];

	//U and V populating
	while (counter < nt)
	{
		//Runs over local workload
		for (i=0; i<problemSize; i++)
		{
			xLocation = startingLocation % nx;
			yLocation = startingLocation / ny;

			localData[i].p = pressureFind(xLocation, yLocation);
			localData[i].u = xMomentumFind(xLocation, yLocation);
			localData[i].v = yMomentumFind(xLocation, yLocation);
			
		}

//		MPI_allgather
	}



	MPI_Finalize();

	free(pressureGradient);
	free(xVelocity);
	free(yVelocity);

	return 0;
}
