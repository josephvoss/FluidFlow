#include "Simulation.h"

//Simulation::Simulation(void)

Simulation::Simulation()
{
	dx = (float) 2/(nx-1);
	dy = (float) 2/(nx-1);
	dt = 0.0001;

	//Workload distribution
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	numCells = problemSize/size;
	leftOvers = problemSize % size;
	for (i=1; i<leftOvers+1; i++)
		if (i == myRank)
			numCells += 1;

	//init local data Arrays
	recCounts = (int*) malloc(sizeof(int)*size);
	for (i=0; i<size; i++)
		recCounts[i] = 0;

	localPrePresData = (double*) malloc(sizeof(double)*numCells);
	localB = (double*) malloc(sizeof(double)*numCells);
	localVelData = (datumPoint*) malloc(sizeof(datumPoint)*numCells);
	for (i=0; i<nit; i++)
		solvedPrePresData[i] = (double*) malloc(sizeof(double)*problemSize);

	startingLocation = numCells * myRank;

	for (i=0; i<problemSize; i++)
	{
		solvedVelData[0][i].p = 1;
		solvedVelData[0][i].u = 0;
		solvedVelData[0][i].v = 0;
		solvedPrePresData[0][i] = 1;
	}

	for (i=0; i<numCells; i++)
	{
		localVelData[i].p = 1;
		localVelData[i].u = 0;
		localVelData[i].v = 0;
		localPrePresData[i] = 1;
		localB[i] = 0;
	}
}

double Simulation::buildUpB(int xLocation, int yLocation)
{
	//Questionable BC
	if (yLocation == ny-1)
		return 0;
	if (yLocation == 0)
		return 0;
	//Periodic
	if (xLocation == nx-1)
	{
		uijp1n = solvedVelData[counter][0+yLocation*ny].u;
		vijp1n = solvedVelData[counter][0+yLocation*ny].v;
	}
	else
	{
		uijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].u;
		vijp1n = solvedVelData[counter][xLocation+1+yLocation*ny].v;
	}

	if (xLocation == 0)
	{
		uijm1n = solvedVelData[counter][nx-1+yLocation*ny].u;
		vijm1n = solvedVelData[counter][nx-1+yLocation*ny].v;
	}
	else
	{
		uijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].u;
		vijm1n = solvedVelData[counter][xLocation-1+yLocation*ny].v;
	}

	//Aliases
	uijn = solvedVelData[counter][xLocation+yLocation*ny].u;
	vijn = solvedVelData[counter][xLocation+yLocation*ny].v;

	uip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].u;
	vip1jn = solvedVelData[counter][xLocation+(yLocation+1)*ny].v;
	uim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].u;
	vim1jn = solvedVelData[counter][xLocation+(yLocation-1)*ny].v;

	double b = rho*(1/dt*((uijp1n - uijm1n)/(2*dx) + (vip1jn - vim1jn)/(2*dy)) - pow((uijp1n - uijm1n)/(2*dx),2) - 2*(uip1jn - uim1jn)/(2*dy) * (vijp1n - vijm1n)/(2*dx) - pow((vip1jn - vim1jn)/(2*dy),2));
	return b;

}

double Simulation::pressureSolve(int xLocation, int yLocation, int i)
{
	//Questionable BC
	if (yLocation == ny-1)
		return 1;
	if (yLocation == 0)
		return 1;
	//Periodic
	if (xLocation == nx-1)
		pijp1n = solvedPrePresData[subCounter-1][1+yLocation*ny];
	else
		pijp1n = solvedPrePresData[subCounter-1][xLocation+1+yLocation*ny];

	if (xLocation == 0)
		pijm1n = solvedPrePresData[subCounter-1][nx-2+yLocation*ny];
	else
		pijm1n = solvedPrePresData[subCounter-1][(xLocation-1)+yLocation*ny];

	//Aliases
	pim1jn = solvedPrePresData[subCounter-1][xLocation+(yLocation-1)*ny];
	pip1jn = solvedPrePresData[subCounter-1][xLocation+(yLocation+1)*ny];

	//For distributed workload collective needs to be done at end of each ptimestep
	//Will square the communication necessary
	//Conversely, could have one dedicated proc for solving pressure, while all others
	//do velocity calculations
		//Isn't this still an issue? The dedicated proc will need to do nx*ny*nt
		//calculations
		//All others will have to do ny*nx/P. Smaller than atleast a factor of 100
		//Using OMP could only reduce this time by potentially 16.
	//Treat pressure solving as velocity solving

	double x = ((pip1jn+pim1jn)*dy*dy+(pijp1n+pijm1n)*dx*dx)/(2*(dx*dx+dy*dy)) - dx*dx*dy*dy/(dx*dx+dy*dy)*localB[i];
	return x;
};

double Simulation::yMomentumSolve(int xLocation, int yLocation)
{
	if (xLocation == nx-1)
		return 0;
	if (xLocation == 0)
		return 0;
	//Periodic
	if (xLocation == nx-1)
	{
		pijp1n = solvedPrePresData[nit-1][1+yLocation*ny];
		uijp1n = solvedVelData[counter-1][1+yLocation*ny].u;
		vijp1n = solvedVelData[counter-1][1+yLocation*ny].v;
	}
	else
	{
		uijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
		vijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;
		pijp1n = solvedPrePresData[nit-1][xLocation+1+yLocation*ny];
	}
	if (xLocation == 0)
	{
		pijm1n = solvedPrePresData[nit-1][nx-2+yLocation*ny];
		uijm1n = solvedVelData[counter-1][nx-2+yLocation*ny].u;
		vijm1n = solvedVelData[counter-1][nx-2+yLocation*ny].v;
	}
	else
	{
		uijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
		vijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
		pijm1n = solvedPrePresData[nit-1][xLocation-1+yLocation*ny];
	}

	//Aliases
	uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	pip1jn = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	pim1jn = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];
	uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;

	return vijn - uijn*dt/dx*(vijn - vim1jn) - vijn*dt/dy*(vijn - vijm1n) -
	dt/(rho*2*dy)*(pijp1n - pijp1n) + nu*(dt/(dx*dx)*(vip1jn - 2*vijn + vim1jn) + 
	dt/(dy*dy)*(vijp1n - 2*vijn + vijm1n));
}

double Simulation::xMomentumSolve(int xLocation, int yLocation)
{
	if (yLocation >= ny-1)
		return 0;
	if (yLocation == 0)
		return 0;
	if (xLocation == nx-1)
	{
		pijm1n = solvedPrePresData[nit-1][0+yLocation*ny];
		uijm1n = solvedVelData[counter-1][0+yLocation*ny].u;
		vijm1n = solvedVelData[counter-1][0+yLocation*ny].v;
	}
	else
	{
		uijm1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
		vijm1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;
		pijm1n = solvedPrePresData[nit-1][xLocation+1+yLocation*ny];
	}
	if (xLocation == 0)
	{
		pijp1n = solvedPrePresData[nit-1][nx-1+yLocation*ny];
		uijp1n = solvedVelData[counter-1][nx-1+yLocation*ny].u;
		vijp1n = solvedVelData[counter-1][nx-1+yLocation*ny].v;
	}
	else
	{
		uijp1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
		vijp1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
		pijp1n = solvedPrePresData[nit-1][xLocation-1+yLocation*ny];
	}
	//Aliases
	uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	pip1jn = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	pim1jn = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];
	uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;

	double x = uijn - uijn*dt/dx*(uijn - uim1jn) - vijn*dt/dy*(uijn - uijm1n) -
	dt/(rho*2*dx)*(pip1jn - pim1jn) + nu*(dt/(dx*dx)*(uip1jn - 2*uijn + uim1jn) +
	dt/(dy*dy)*(uijp1n - 2*uijn + uijm1n)) + F*dt;
	return x;
}

void Simulation::iterate(void)
{
	//Populate recCount and displs arrays
	int displs[size];
	int* pNumCells = &numCells;
	MPI_Allgather(pNumCells, 1, MPI_INT, recCounts, 1, MPI_INT, MPI_COMM_WORLD);
	int sum = 0;
	int i;

	//Defining MPI data type
	MPI_Datatype newType;
	MPI_Type_contiguous(3,MPI_DOUBLE, &newType);
	MPI_Type_commit(&newType);

	//Length from receive buffer to put data
	for (i=0; i<size; i++)
	{
		displs[i] = sum;
		sum += recCounts[i];
	}		

	//U and V populating
	while (counter < nt)
	{
		//Pressure Populating
		for (i=0; i<numCells; i++)
		{
			xLocation = (startingLocation+i)% nx;
			yLocation = (startingLocation+i)/ nx;
			solvedPrePresData[0][i] = solvedVelData[counter-1][i].p; //pressurePreSolve uses this array, iterates in ptime
			localB[i] = buildUpB(xLocation, yLocation);
		}
		//Solve pressure for n
		while (subCounter < nit)
		{
			for (i=0; i<numCells; i++)
			{
				xLocation = (startingLocation+i)% nx;
				yLocation = (startingLocation+i)/ nx;

				localPrePresData[i] = pressureSolve(xLocation, yLocation,i);
			}
			MPI_Allgatherv(localPrePresData, numCells, MPI_DOUBLE, solvedPrePresData[subCounter], recCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
			subCounter += 1;

			// dP/dy  bc needs to be satified after all values collected
			for (i=0; i<problemSize; i++)
			{
				xLocation = i% nx;
				yLocation = i/ nx;

				if (yLocation == 0)
					solvedPrePresData[subCounter-1][xLocation] = solvedPrePresData[subCounter-1][xLocation+(1)*ny];
				if (yLocation == ny-1)
					solvedPrePresData[subCounter-1][xLocation+(ny-1)*ny] = solvedPrePresData[subCounter-1][xLocation+(ny-2)*ny];
			}
		}
		
		//Runs over local workload
		for (i=0; i<numCells; i++)
		{
			xLocation = (startingLocation+i) % nx;
			yLocation = (startingLocation+i) / ny;

			localVelData[i].u = xMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].v = yMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].p = solvedPrePresData[nit-1][startingLocation+i]; //needs to be for n
	
		}

		MPI_Allgatherv(localVelData, numCells, newType, &(solvedVelData[counter][0]), recCounts, displs, newType, MPI_COMM_WORLD);

		subCounter = 1;
		counter += 1;
	}

	MPI_Type_free(&newType);
}
