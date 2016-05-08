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
	//Aliases
	double uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;

	if (yLocation == ny-1)
		return 0;
	if (yLocation == 0)
		return 0;
	if (xLocation == ny-1)
		return 0;
	if (xLocation == 0)
		return 0;

	double b = rho*(1/dt*((uip1jn - uim1jn)/(2*dx) + (vijp1n - vijm1n)/(2*dy)) - pow((uip1jn - uim1jn)/(2*dx),2) - 2*(uijp1n - uijm1n)/(2*dy) * (vip1jn - vim1jn)/(2*dx) - pow((vijp1n - vijm1n)/(2*dy),2));
	printf("%d: (%d, %d) %f\n", counter, xLocation, yLocation, b);
	return b;

}

double Simulation::pressureSolve(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;
	double pip1jn = solvedPrePresData[subCounter-1][xLocation+(yLocation+1)*ny];
	double pim1jn = solvedPrePresData[subCounter-1][xLocation+(yLocation-1)*ny];
	double pijp1n = solvedPrePresData[subCounter-1][xLocation+(yLocation+1)*ny];
	double pijm1n = solvedPrePresData[subCounter-1][xLocation+(yLocation-1)*ny];

	//For distributed workload collective needs to be done at end of each ptimestep
	//Will square the communication necessary
	//Conversely, could have one dedicated proc for solving pressure, while all others
	//do velocity calculations
		//Isn't this still an issue? The dedicated proc will need to do nx*ny*nt
		//calculations
		//All others will have to do ny*nx/P. Smaller than atleast a factor of 100
		//Using OMP could only reduce this time by potentially 16.
	//Treat pressure solving as velocity solving

	if (yLocation == ny-1)
		return 1;
	if (yLocation == 0)
		return 1;
	if (xLocation == ny-1)
		return 1;
	if (xLocation == 0)
		return 1;

	double x = ((pip1jn+pim1jn)*dy*dy+(pijp1n+pijm1n)*dx*dx)/(2*(dx*dx+dy*dy)) - dx*dx*dy*dy/(dx*dx+dy*dy)*localB[xLocation+yLocation*ny];
//	printf("%d.%d: (%d, %d) %f\n", counter, subCounter, xLocation, yLocation, x);
	return x;
};

double Simulation::yMomentumSolve(int xLocation, int yLocation)
{
	if (yLocation == ny-1)
		return 0;
	if (yLocation == 0)
		return 0;
	if (xLocation == nx-1)
		return 0;
	if (xLocation == 0)
		return 0;

	//Aliases
	double uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;
	double pip1jn = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	double pim1jn = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];
	double pijp1n = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	double pijm1n = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];

	return vijn - uijn*dt/dx*(vijn - vim1jn) - vijn*dt/dy*(vijn - vijm1n) -
	dt/(rho*2*dy)*(pijp1n - pijp1n) + nu*(dt/(dx*dx)*(vip1jn - 2*vijn + vim1jn) + 
	dt/(dy*dy)*(vijp1n - 2*vijn + vijm1n));
}

double Simulation::xMomentumSolve(int xLocation, int yLocation)
{
	//Aliases
	double uijn = solvedVelData[counter-1][xLocation+yLocation*ny].u;
	double vijn = solvedVelData[counter-1][xLocation+yLocation*ny].v;
	double uim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].u;
	double uip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].u;
	double uijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].u;
	double uijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].u;
	double vim1jn = solvedVelData[counter-1][xLocation+(yLocation-1)*ny].v;
	double vip1jn = solvedVelData[counter-1][xLocation+(yLocation+1)*ny].v;
	double vijm1n = solvedVelData[counter-1][xLocation-1+yLocation*ny].v;
	double vijp1n = solvedVelData[counter-1][xLocation+1+yLocation*ny].v;
	double pip1jn = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	double pim1jn = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];
	double pijp1n = solvedPrePresData[nit-1][xLocation+(yLocation+1)*ny];
	double pijm1n = solvedPrePresData[nit-1][xLocation+(yLocation-1)*ny];

	if (yLocation == ny-1)
		return 0;
	if (yLocation == 0)
		return 0;
	if (xLocation == ny-1)
		return 0;
	if (xLocation == 0)
		return 0;

	return  uijn - uijn*dt/dx*(uijn - uim1jn) - vijn*dt/dy*(uijn - uijm1n) -
	dt/(rho*2*dx)*(pip1jn - pim1jn) + nu*(dt/(dx*dx)*(uip1jn - 2*uijn + uim1jn) +
	dt/(dy*dy)*(uijp1n - 2*uijn + uijm1n)) + F*dt;
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
			solvedPrePresData[0][i] = solvedVelData[counter-1][i].p; //pressurePreSolve uses this array
			localB[i] = buildUpB(xLocation, yLocation);
		}
		i=numCells;
		//Solve pressure for n+1
		while (subCounter < nit)
		{
			for (i=0; i<numCells; i++)
			{
				xLocation = (startingLocation+i)% nx;
				yLocation = (startingLocation+i)/ nx;

				localPrePresData[i] = pressureSolve(xLocation, yLocation); //needs to be for n-1
				//Uses B, from n-1
			}
			MPI_Allgatherv(localPrePresData, numCells, MPI_DOUBLE, solvedPrePresData[subCounter], recCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
			subCounter += 1;
		}
		
		//Runs over local workload
		for (i=0; i<numCells; i++)
		{
			xLocation = (startingLocation+i) % nx;
			yLocation = (startingLocation+i) / ny;

			localVelData[i].u = xMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].v = yMomentumSolve(xLocation, yLocation); //needs to be for n
			localVelData[i].p = localPrePresData[i]; //needs to be for n
	
		}

//		MPI_allgather
		MPI_Allgatherv(localVelData, numCells, newType, solvedVelData[counter], recCounts, displs, newType, MPI_COMM_WORLD);
/*		int x,y;
		for (i=0; i<problemSize; i++)
		{
			x = i/nx; y = i%nx;
			//swap x and y for better hdf5 formatting, unnecessary
			solvedPMat[counter][x][y] = i;//solvedVelData[counter][startingLocation+i].p;
			solvedUMat[counter][x][y] = i;// solvedVelData[counter][i].u;
			solvedVMat[counter][x][y] = i;//solvedVelData[counter][i].v;
		}
*/
		subCounter = 1;
		counter += 1;
	}

	MPI_Type_free(&newType);
}
