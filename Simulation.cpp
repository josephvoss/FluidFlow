#include "Simulation.h"

Simulation::Simulation(void)
{
	nx = 10;
	ny = 10;
	nt = 10;
	dx = 2/(nx-1);
	dy = 2/(ny-1);
	dt = 0.0001;

	//Physical Values
	rho = 1; 
	nu = .1;
	F = //needed? 

	myRank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();

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

	//Populate recCount and displs arrays
	MPI::COMM_WORLD.Allgather(numCells, 1, MPI::INT, recCounts, size, MPI::INT);
	int sum = 0;
	for (i=0; i<size; i++)
	{
		sum += recCount[i]*sizeof(int);
		displs[i] = sum;
	}		
}

void Simulation::buildUpB(int xLocation, int yLocation)
{
	//Aliases
	datumPoints lastgbl[problemSize] = allSolvedData[counter-1];
	datumPoint gbl[problemSize] = globalSolvedData;
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double uim1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double uip1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uijm1n = gbl[xLocation-1+yLocation*ny].u;
	double uijp1n = glb[xLocation+1+yLocation*ny].u;
	double vim1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vip1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vijm1n = gbl[xLocation-1+yLocation*ny].v;
	double vijp1n = glb[xLocation+1+yLocation*ny].v;

	return rho*(1/dt*((uip1jn - uim1jn)/(2*dx) + (vijp1n - vijm1n)/(2*dy)) -\
	pow((uip1jn - uim1j)/(2*dx),2) - 2*(uijp1n - uijm1n)/(2*dy) * \
	(vip1jn - vim1jn)/(2*dx) - pow((vijp1n - vijm1n)/(2*dy),2));

}

float Simulation::pressurePreSolve(int xLocation, int yLocation)
{
	int pseudoTime;
	//Aliases
	datumPoint gbl[problemSize] = solvedVelData[counter];
	double pre[problemSize] = solvedPrePresData[subCounter-1]
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double uim1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double uip1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uijm1n = gbl[xLocation-1+yLocation*ny].u;
	double uijp1n = glb[xLocation+1+yLocation*ny].u;
	double vim1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vip1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vijm1n = gbl[xLocation-1+yLocation*ny].v;
	double vijp1n = glb[xLocation+1+yLocation*ny].v;
	double pijp1n = glb[xLocation+1+yLocation*ny].p;
	double pijm1n = glb[xLocation-1+yLocation*ny].p;
	double pip1jn = pre[xLocation+(yLocation+1)*ny];
	double pim1jn = pre[xLocation+(yLocation-1)*ny];

	//For distributed workload collective needs to be done at end of each ptimestep
	//Will square the communication necessary
	//Conversely, could have one dedicated proc for solving pressure, while all others
	//do velocity calculations
		//Isn't this still an issue? The dedicated proc will need to do nx*ny*nt
		//calculations
		//All others will have to do ny*nx/P. Smaller than atleast a factor of 100
		//Using OMP could only reduce this time by potentially 16.
	//Treat pressure solving as velocity solving
	return ((pip1jn+pim1jn)*dy*dy+(pijp1n-pijm1n)*dx*dx)/(2*(dx*dx+dy*dy));// - dx*dx*dy*dy/(dx*dx+dy*dy)*b
};

float Simulation::xMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	datumPoint gbl[problemSize] = solvedVelData[counter-1];
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double uim1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double uip1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uijm1n = gbl[xLocation-1+yLocation*ny].u;
	double uijp1n = glb[xLocation+1+yLocation*ny].u;
	double vim1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vip1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vijm1n = gbl[xLocation-1+yLocation*ny].v;
	double vijp1n = glb[xLocation+1+yLocation*ny].v;
	double pijp1n = glb[xLocation+1+yLocation*ny].p;
	double pijm1n = glb[xLocation-1+yLocation*ny].p;
	double pip1jn = glb[xLocation+(yLocation+1)*ny].p;
	double pim1jn = glb[xLocation+(yLocation-1)*ny].p;

	return  vijn - uijn*dt/dx*(vijn - vim1jn) - vijn*dt/dy*(vijn - vijm1n) -\
	dt/(rho*2*dy)*(pijp1n - pijp1n) + nu*(dt/(dx*dx)*(vip1jn - 2*vijn + vim1jn)\
	dt/(rho*dy*dy)*(vijp1n - 2*vijn + vijm1n));
}

float Simulation::yMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	datumPoint gbl[problemSize] = globalSolvedData;
	double uijn = gbl[xLocation+yLocation*ny].u;
	double vijn = gbl[xLocation+yLocation*ny].v;
	double uim1jn = gbl[xLocation+(yLocation-1)*ny].u;
	double uip1jn = gbl[xLocation+(yLocation+1)*ny].u;
	double uijm1n = gbl[xLocation-1+yLocation*ny].u;
	double uijp1n = glb[xLocation+1+yLocation*ny].u;
	double vim1jn = gbl[xLocation+(yLocation-1)*ny].v;
	double vip1jn = gbl[xLocation+(yLocation+1)*ny].v;
	double vijm1n = gbl[xLocation-1+yLocation*ny].v;
	double vijp1n = glb[xLocation+1+yLocation*ny].v;
	double pijp1n = glb[xLocation+1+yLocation*ny].p;
	double pijm1n = glb[xLocation-1+yLocation*ny].p;
	double pip1jn = glb[xLocation+(yLocation+1)*ny].p;
	double pim1jn = glb[xLocation+(yLocation-1)*ny].p;

	return  uijn - uijn*dt/dx*(uijn - uim1jn) - vijn*dt/dy*(uijn - uijm1n) -\
	dt/(rho*2*dx)*(pip1jn - pim1jn) + nu*(dt/(dx*dx)*(uip1jn - 2*uijn + uim1jn)\
	dt/(rho*dy*dy)*(uijp1n - 2*uijn + uijm1n)) + F*dt;
}

double pressureSolve(int xLocation, int yLocation)
{
	return localPresData[xLocation+yLocation*ny] - dx*dx*dy*dy/(dx*dx+dy*dy)*buildUpB(xLocation, yLocation); 
}

void Simulation::iterate(void)
{
	//U and V populating
	while (counter < nt)
	{
		//Pressure Populating
		//Solve pressure for n+1
		while (subCounter < nit)
		{
			for (i=0; i<numCells; i++)
			{
				xLocation = startingLocation % nx;
				yLocation = startingLocation / ny;

				localPresData[i] = pressurePreSolve(xLocation, yLocation); //needs to be for n-1
			}
			MPI::COMM_WORLD.Allgatherv(localPresData, numCells, MPI::DOUBLE, solvedPrePresData[subCounter], recCounts, disps, MPI::DOUBLE);
			subCounter += 1;
		}
		
		//Runs over local workload
		for (i=0; i<numCells; i++)
		{
			xLocation = startingLocation % nx;
			yLocation = startingLocation / ny;

			localVelData[i].u = xMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].v = yMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].p = pressureSolve(xLocation, yLocation); //needs to be for n
	
		}

//		MPI_allgather
		counter += 1;
	}
}
