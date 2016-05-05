#include "Simulation.h"

//Simulation::Simulation(void)

void Simulation::init(void)
{
	myRank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();

	//Workload distribution
	numCells = problemSize/size;
	leftOvers = problemSize % size;
	int i;
	for (i=1; i<leftOvers+1; i++)
		if (i == myRank)
			numCells += 1;

	startingLocation = numCells * myRank;

	for (i=0; i<problemSize; i++)
	{
		solvedVelData[0][i].p = 0;
		solvedVelData[0][i].u = 0;
		solvedVelData[0][i].v = 0;
		solvedPrePresData[0][i] = 0;
	}

	for (i=0; i<numCells; i++)
	{
		localVelData[i].p = 0;
		localVelData[i].u = 0;
		localVelData[i].v = 0;
		localPrePresData[i] = 0;
	}
}

double Simulation::buildUpB(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].v;

	return rho*(1/dt*((uip1jn - uim1jn)/(2*dx) + (vijp1n - vijm1n)/(2*dy)) -\
	pow((uip1jn - uim1jn)/(2*dx),2) - 2*(uijp1n - uijm1n)/(2*dy) * \
	(vip1jn - vim1jn)/(2*dx) - pow((vijp1n - vijm1n)/(2*dy),2));

}

double Simulation::pressurePreSolve(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].v;
	double pijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].p;
	double pijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].p;
	double pip1jn = solvedPrePresData[subCounter][xLocation+(yLocation+1)*ny];
	double pim1jn = solvedPrePresData[subCounter][xLocation+(yLocation-1)*ny];

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

double Simulation::xMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].v;
	double pijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].p;
	double pijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].p;
	double pip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].p;
	double pim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].p;

	return vijn - uijn*dt/dx*(vijn - vim1jn) - vijn*dt/dy*(vijn - vijm1n) -
	dt/(rho*2*dy)*(pijp1n - pijp1n) + nu*(dt/(dx*dx)*(vip1jn - 2*vijn + vim1jn) + 
	dt/(rho*dy*dy)*(vijp1n - 2*vijn + vijm1n));
}

double Simulation::yMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].v;
	double pijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].p;
	double pijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].p;
	double pip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].p;
	double pim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].p;

	return  uijn - uijn*dt/dx*(uijn - uim1jn) - vijn*dt/dy*(uijn - uijm1n) -
	dt/(rho*2*dx)*(pip1jn - pim1jn) + nu*(dt/(dx*dx)*(uip1jn - 2*uijn + uim1jn) +
	dt/(rho*dy*dy)*(uijp1n - 2*uijn + uijm1n)) + F*dt;
}

double Simulation::pressureSolve(int xLocation, int yLocation)
{
	return localPrePresData[xLocation+yLocation*ny] - dx*dx*dy*dy/(dx*dx+dy*dy)*buildUpB(xLocation, yLocation); 
}

void Simulation::iterate(void)
{
	int recCounts[size];
	int displs[size];

	//Populate recCount and displs arrays
//	MPI::COMM_WORLD.Allgather(&numCells[0], 1, MPI::INT, (void*) &recCounts[0], size, MPI::INT);
	int sum = 0;
	int i;
	for (i=0; i<size; i++)
	{
		sum += recCounts[i]*sizeof(int);
		displs[i] = sum;
	}		

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

				//localPrePresData[i] = pressurePreSolve(xLocation, yLocation); //needs to be for n-1
			}
//			MPI::COMM_WORLD.Allgatherv(&localPrePresData[0], numCells, MPI::DOUBLE, solvedPrePresData[subCounter], (const int*) &recCounts[0], displs, MPI::DOUBLE);
			subCounter += 1;
		}
		
		//Runs over local workload
		for (i=0; i<numCells; i++)
		{
			xLocation = startingLocation % nx;
			yLocation = startingLocation / ny;

		//	localVelData[i].u = xMomentumSolve(xLocation, yLocation); //needs to be for n
		//	localVelData[i].v = yMomentumSolve(xLocation, yLocation); //needs to be for n
		//	localVelData[i].p = pressureSolve(xLocation, yLocation); //needs to be for n
	
		}

//		MPI_allgather
		counter += 1;
	}
}
