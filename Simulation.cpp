#include "Simulation.h"

Simulation::Simulation(void)
{
	nx = 101;
	ny = 101;
	nt = 100;
	dx = 2/(nx-1);
	dy = 2/(ny-1);
	dt = 0.0001;

	//Physical Values
	rho = 
	mu = 
	F = //needed? 

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	counter = 0;

	//Workload distribution
	numCells = problemSize/size;
	leftOvers = problemSize % size;
	for (i=1; i<leftOvers+1; i++)
		if (i == myRank)
			numCells += 1;

	startingLocation = numCells * rank;

	nGlobalSolvedData = (datumPoint*) malloc(sizeof(datumPoint)*nx*ny);
	localData = (datumPoint*) malloc(sizeof(datumPoint)*numCells);
	allSolvedData = (datumPoint*) malloc(sizeof(datumPoint)*nx*ny*nt);
}

void Simulation::buildUpB(int xLocation, int yLocation)
{
	//Aliases
	datumPoints lastgbl[problemSize] = allSolvedData[counter-1];
	datumPoint gbl[problemSize] = globalSolvedData;
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double ui-1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double ui+1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uij-1n = gbl[xLocation-1+yLocation*ny].u;
	double uij+1n = glb[xLocation+1+yLocation*ny].u;
	double vi-1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vi+1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vij-1n = gbl[xLocation-1+yLocation*ny].v;
	double vij+1n = glb[xLocation+1+yLocation*ny].v;

	b=rho*(1/dt*((ui+1jn - ui-1jn)/(2*dx) + (vij+1n - vij-1n)/(2*dy)) -\
	pow((ui+1jn - ui-1j)/(2*dx),2) - 2*(uij+1n - uij-1n)/(2*dy) * \
	(vi+1jn - vi-1jn)/(2*dx) - pow((vij+1n - vij-1n)/(2*dy),2));


}

float Simulation::pressureSolve(int xLocation, int yLocation)
{
	int pseudoTime;
	//For distributed workload collective needs to be done at end of each ptimestep
	//Will square the communication necessary
	//Conversely, could have one dedicated proc for solving pressure, while all others
	//do velocity calculations
		//Isn't this still an issue? The dedicated proc will need to do nx*ny*nt
		//calculations
		//All others will have to do ny*nx/P. Smaller than atleast a factor of 100
		//Using OMP could only reduce this time by potentially 16.
	for(pseudoTime=0; pseudoTime < nt; pseudoTime++)
	{
		for(i=0; i<ny; i++)
		{
			for(j=0; j<nx; j++)
			{
				p=((pi+1jn+pi-1jn)*dy*dy+(pij+1n-pij-1n)*dx*dx)/(2*(dx*dx+dy*dy)) - dx*dx*dy*dy/(dx*dx+dy*dy)*b
			}
		}
	}
};

float Simulation::xMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	datumPoint gbl[problemSize] = globalSolvedData;
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double ui-1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double ui+1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uij-1n = gbl[xLocation-1+yLocation*ny].u;
	double uij+1n = glb[xLocation+1+yLocation*ny].u;
	double vi-1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vi+1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vij-1n = gbl[xLocation-1+yLocation*ny].v;
	double vij+1n = glb[xLocation+1+yLocation*ny].v;
	double pij+1n = glb[xLocation+1+yLocation*ny].p;
	double pij-1n = glb[xLocation-1+yLocation*ny].p;
	double pi+1jn = glb[xLocation+(yLocation+1)*ny].p;
	double pi-1jn = glb[xLocation+(yLocation-1)*ny].p;

	return  vijn - uijn*dt/dx*(vijn - vi-1jn) - vijn*dt/dy*(vijn - vij-1n) -\
	dt/(rho*2*dy)*(pij+1n - pij+1n) + nu*(dt/(dx*dx)*(vi+1jn - 2*vijn + vi-1jn)\
	dt/(rho*dy*dy)*(vij+1n - 2*vijn + vij-1n));
}

float Simulation::yMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	datumPoint gbl[problemSize] = globalSolvedData;
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double ui-1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double ui+1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uij-1n = gbl[xLocation-1+yLocation*ny].u;
	double uij+1n = glb[xLocation+1+yLocation*ny].u;
	double vi-1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vi+1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vij-1n = gbl[xLocation-1+yLocation*ny].v;
	double vij+1n = glb[xLocation+1+yLocation*ny].v;
	double pij+1n = glb[xLocation+1+yLocation*ny].p;
	double pij-1n = glb[xLocation-1+yLocation*ny].p;
	double pi+1jn = glb[xLocation+(yLocation+1)*ny].p;
	double pi-1jn = glb[xLocation+(yLocation-1)*ny].p;

	return  uijn - uijn*dt/dx*(uijn - ui-1jn) - vijn*dt/dy*(uijn - uij-1n) -\
	dt/(rho*2*dx)*(pi+1jn - pi-1jn) + nu*(dt/(dx*dx)*(ui+1jn - 2*uijn + ui-1jn)\
	dt/(rho*dy*dy)*(uij+1n - 2*uijn + uij-1n)) + F*dt;
}

void Simulation::iterate(void)
{
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
		counter += 1;
	}
}
