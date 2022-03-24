#include <mpi.h>
#include <vector>
#include <iostream>

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
MPI_Comm comm = MPI_COMM_WORLD;

int rank;
MPI_Comm_rank(comm, &rank);

// scatter
int seed;
std::vector<int> randseeds(3);
if (rank==0){
	randseeds = {10,20,30};
}
MPI_Scatter(&randseeds[0], 1, MPI_INTEGER, &seed, 1, MPI_INTEGER, 0, comm);

std::cout<<"this is rank "<<rank<<" and I received seed "<<seed<<std::endl;



// reduce
/*std::vector <double> Tconf(3);
std::vector <double> Tconf_avg(3);

for(int i=0; i<Tconf.size(); ++i){
	Tconf[i]=i*rank;
}

for(int i=0; i<Tconf_avg.size(); ++i){
	MPI_Reduce(&Tconf[i], &Tconf_avg[i], 1, MPI_DOUBLE, MPI_SUM, 0, comm);
}
if(rank==0){
	std::cout<<"Results:"<<"\n";
	for(int i=0; i<Tconf_avg.size(); ++i){
			std::cout<< "i="<<i<<" :  "<<Tconf_avg[i]<<std::endl;
	}
}*/


// send a vector
/*MPI_Status status;
std::vector<double> xvec(10);
std::vector<double> resultvec(10);
if(rank == 1){
	for(int i=0; i<xvec.size(); ++i){
		xvec[i] = i;	
	}
	MPI_Ssend(&xvec[0], xvec.size(), MPI_DOUBLE, 0, 0, comm);
}
else if(rank==0){
	MPI_Recv(&resultvec[0], resultvec.size(), MPI_DOUBLE, 1, 0, comm, &status);
	std::cout<<"received!"<<std::endl;
	for(int i=0; i<resultvec.size(); ++i){std::cout<<resultvec[i]<<"\n"<<std::endl;}
}
*/


MPI_Finalize();

return 0;

}