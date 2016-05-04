#include <"Simulation.h">

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
	Simulation workBench;
	workBench.iterate();

	MPI_Finalize();
	return 0;
}
