/* 
This script samples the posterior density of the means of two independent Gaussians used to build a 1D Gaussian mixture model.
It uses Metropolized Langevin integrators (given by the "OBABO" scheme) and uses full gradients or stochastic gradients obtained by data subsampling.
Using K processors, it draws K trajectories simultaneously, and averages over them (also in time) before printing out the results.

This version tracks the two means and the Metropolis acceptance probability in time.
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <chrono>
#include <math.h> 
#include <iomanip>
#include <algorithm>
#include<limits>
#include <numeric>
#include <mpi.h>


using namespace std;

const double PI = 3.141592653589793;


mt19937 twister;

struct params {
double mu1,mu2,p1,p2;
};
struct forces{
double fmu1, fmu2;
};

struct measurement{
vector <double> mu1, mu2, acceptP;
};

/*struct measurement{
vector <double> mu1, mu2;
};*/


double likelihood_GM(const params& theta, const double sig1, const double sig2, const double a1, const double a2, const double x);
double U_pot_GM( const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  const vector <double>& Xdata, const double sig0);
forces get_noisy_force_GM(const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  			  vector <double>& Xdata, const size_t B, const double sig0);
measurement OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	
									const size_t B, const double sig1, const double sig2, const double a1, const double a2, const size_t n_meas, const double sig0);
measurement MOBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	const size_t B, 
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0);
measurement OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0);
vector <double> read_dataset(string datafile);


int main(int argc, char *argv[]){

double sig1 = 3;			// GM params
double sig2 = 0.5;
double a1 = 0.8;
double a2 = 0.2;

double sig0 = 5;			// Gauss. prior std.dev.

params param0{-4,3,0,0};	// initial conditions


double T = 1;
double gamma = 1;
string datafile = "GM_data_5000.csv";
bool tavg = true;								// time average after ensemble average?
int n = 5e5; 									// time average over last n values
int ndist = 10000;							// write out any ndist-th result entry to file 
//int ndist = 1;

int method = atoi(argv[1]);   // must be 1 (OBABO), 2 (MOBABO), or 3 (OMBABO)
size_t L = atoi(argv[2]);
string SF = argv[3];				// must be "A", "R", or "0"
double h = atof(argv[4]);
size_t N_0 = atol(argv[5]);
size_t B = atof(argv[6]);
int n_meas = atoi(argv[7]);	// store sample any n_meas steps (governs RAM used during execution, as opposed to ndist which governs size of output file)

//size_t N = N_0/(L+2);
size_t N = N_0;

vector <double> Xdata = read_dataset(datafile);

MPI_Init(&argc, &argv);				// initialize MPI, use rank as random seed
MPI_Comm comm = MPI_COMM_WORLD;
int rank, nr_proc;
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &nr_proc);
int randomseed = rank;

seed_seq seq{1,20,3200,403,5*randomseed+1,12000,73667,9474+randomseed,19151-randomseed};
vector<std::uint32_t> seeds(1);
seq.generate(seeds.begin(), seeds.end());
twister.seed(seeds.at(0)); 

string label;
measurement results; // stores mus and acceptance probs.
if ( method == 1 ){
	results = OBABO_simu(param0, N_0, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, n_meas, sig0);
	label = "OBABO";
}
else if ( method == 2 ){
	results = MOBABO_simu(param0, N, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas, sig0);
	label = "MOBABO_SF" + SF +"_L" + to_string(L);
}
else if ( method == 3 ){
	results = OMBABO_simu(param0, N, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas, sig0);
	label = "OMBABO_SF" + SF +"_L" + to_string(L);
}
else {
	cout << "No valid method number" << endl;
	return 0;
}

cout<<"rank "<<rank<<" reached barrier"<<endl;
MPI_Barrier(comm);

// average over results of different processors and print out file
measurement results_avg;
if( rank==0 ){
	 results_avg.mu1.resize(results.mu1.size());  // do only on rank 0 to save RAM 
	 results_avg.mu2.resize(results.mu2.size());
	 results_avg.acceptP.resize(results.acceptP.size());
	 
}	 

MPI_Reduce(&results.mu1[0], &results_avg.mu1[0], results.mu1.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
MPI_Reduce(&results.mu2[0], &results_avg.mu2[0], results.mu2.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
MPI_Reduce(&results.acceptP[0], &results_avg.acceptP[0], results.acceptP.size(), MPI_DOUBLE, MPI_SUM, 0, comm);

for (int i=0; i<results_avg.mu1.size(); ++i){
	results_avg.mu1[i] /= nr_proc;
	results_avg.mu2[i] /= nr_proc;
	results_avg.acceptP[i] /= nr_proc;
}

if( rank==0 ){	

	// time-average results
	if( tavg == true ){	
		cout << "Time averaging...\n"; 
		
		for(int i=results_avg.mu1.size()-1; i>=0; i-=ndist){
			if(i<=n-1) n=i;
			for(int j=i-n; j<i; ++j){
				results_avg.mu1[i] += results_avg.mu1[j];
				results_avg.mu2[i] += results_avg.mu2[j];
				results_avg.acceptP[i] += results_avg.acceptP[j];					
			}
				results_avg.mu1[i] /= n+1 ;
				results_avg.mu2[i] /= n+1 ;  		
				results_avg.acceptP[i] /= n+1 ; 				
		}	
		
	}			

	
	// Print results
	stringstream stream, stream2;
	string final_label;
	stream << std::fixed << std::setprecision(3) << h;
	if(B==Xdata.size()) final_label = "mus_" + label + "_h"+stream.str()+"_avg"+to_string(nr_proc);
	else 					  final_label = "mus_" + label + "_h"+stream.str()+"_gradnoiseB"+to_string(B)+"_avg"+to_string(nr_proc);			 

	ofstream file {final_label};
	cout<<"Writing to file...\n";
	for(int i=0; i<results_avg.mu1.size(); i+=ndist){
		file << i*n_meas << " " << results_avg.mu1.at(i) << " " << results_avg.mu2.at(i) << " " << results_avg.acceptP.at(i) << "\n";
//		file << i*n_meas << " " << results_avg.mu1.at(i) << " " << results_avg.mu2.at(i) << "\n";
	}
	file.close();
}

MPI_Finalize();

return 0;

}


measurement OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	
									const size_t B, const double sig1, const double sig2, const double a1, const double a2, const size_t n_meas, const double sig0){
    
   cout<<"Starting OBABO simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta = param0; 
	forces force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B, sig0);    // HERE FORCES!!!!

	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1, 1)};		
//	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};
	meas.mu1[0] = theta.mu1;
	meas.mu2[0] = theta.mu2;
	int k=1;

   double a = exp(-1*gamma*h);      

	double Rn1, Rn2;
	normal_distribution<> normal{0,1}; 	
	for(size_t i=1; i<=N; ++i){
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1 + 0.5*h*force.fmu1;  // O+B step		
		theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2 + 0.5*h*force.fmu2;		
														 
		theta.mu1 += h*theta.p1;															// A step							
		theta.mu2 += h*theta.p2;		

		force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B, sig0);   // HERE FORCES!!!!

		theta.p1 += 0.5*h*force.fmu1;														 // B step
		theta.p2 += 0.5*h*force.fmu2;		

		Rn1 = normal(twister);
		Rn2 = normal(twister);								 
		theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;							 // O step
		theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;

		if(i%n_meas == 0 ) {
			meas.mu1[k] = theta.mu1;			
			meas.mu2[k] = theta.mu2;	
			++k;		
		}
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;	
	}

	
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;
	    
   return meas;

}



/*vector <double> MOBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
									 const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0){
   
   cout<<"Starting MOBABO + SF" + SF + " simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 
	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B, sig0);	    // HERE FORCES!!!!

	double Tconfig = -1*(force_curr.fmu1 * theta_curr.mu1  +  force_curr.fmu2 * theta_curr.mu2);
	vector <double> Tconfigs(0);
	Tconfigs.push_back(Tconfig);
	double Tconfig_sum = Tconfig;

   double a = exp(-1*gamma*h);      
	
	double Rn1, Rn2;
	normal_distribution<> normal{0,1};
	uniform_real_distribution<> uniform(0, 1); 
	size_t ctr = 0; 	
	double kin_energy, MH, U1, U0;   // kin. energies are necessary for MH criterion
	params theta;
	forces force; 

	for(size_t i=1; i<N; ++i){
		theta = theta_curr;
		force = force_curr;		
		
		kin_energy = 0;
		
		cout<<"\nIter: "<<i<<endl;
		cout<<" F0: "<<force.fmu1<<"  "<<force.fmu2<<endl;		
		cout<<"Params0: (mu1,p1)=("<<theta.mu1<<", "<<theta.p1<<"),  (mu2,p2)=("<<theta.mu2<<", "<<theta.p2<<")"<<endl;

		for(size_t j=0; j<L; ++j){
			// L OBABO steps			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			
//			cout<<"Rns: "<<Rn1<<"  "<<Rn2<<endl;

			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;

			cout<<"after O-step:\n";
			cout<<"p1, p2 "<<theta.p1<<"  "<<theta.p2<<endl;		
		
			kin_energy -= 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);
		
			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;
			
			cout<<"after B-step:\n";
			cout<<"p1, p2:  "<<theta.p1<<"  "<<theta.p2<<endl;		

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			cout<<"after A-step:\n";
			cout<<"mu1, mu2:  "<<theta.mu1<<"  "<<theta.mu2<<endl;	
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B, sig0);   // HERE FORCES!!!!
//			cout<<"new F: "<<force.fmu1<< "  "<<force.fmu2<<endl;

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			
			
			cout<<"after B-step:\n";
			cout<<"p1, p2:  "<<theta.p1<<"  "<<theta.p2<<endl;	

			kin_energy += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);			
			
			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
//			cout<<"Rns: "<<Rn1<<"  "<<Rn2<<endl;
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  		// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;
			
			cout<<"after O-step:\n";
			cout<<"p1, p2 "<<theta.p1<<"  "<<theta.p2<<endl;				       
		}
		
		// MH criterion
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata, sig0);					//HERE UPOTS!!!
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata, sig0);
		MH = exp( (-1/T) * (U1 - U0 + kin_energy) );					
		cout<<"MH CRITERION:\n";
		cout<<"U0, U1: "<<U0<<"  "<<U1<<endl;
		cout<<"kin energy  "<<kin_energy<<endl;
		cout<<"MH:  "<<MH<<endl;		
		
		if( uniform(twister) < min(1., MH) ){ 					// ACCEPT SAMPLE
				
			cout<<"accepted"<<endl;		
			Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2);
			Tconfig_sum += Tconfig;
			cout<<i<<" accepted"<<endl;
			cout<<"T "<<Tconfig<<endl;	
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	

			theta_curr = theta;
			force_curr = force;
			theta_curr.p1 = SF=="A" ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="A" ? -1*theta_curr.p2 : theta_curr.p2;			

         ctr += 1;
        
		}
		else{  // REJECT SAMPLE
			cout<<"rejected"<<endl;
			Tconfig_sum += Tconfig;
			cout<<i<<" rejected"<<endl;
			cout<<"T "<<Tconfig<<endl;	
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));
			
			theta_curr.p1 = SF=="R" ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="R" ? -1*theta_curr.p2 : theta_curr.p2;					
		}
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return Tconfigs;

}*/


measurement MOBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0){
   
   cout<<"Starting MOBABO + SF" + SF + " simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 

	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B, sig0);	    // HERE FORCES!!!!

	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};		
//	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};
	meas.mu1[0] = theta_curr.mu1;
	meas.mu2[0] = theta_curr.mu2;
	int k=1;

   double a = exp(-1*gamma*h);      
	
	double Rn1, Rn2;
	normal_distribution<> normal{0,1};
	uniform_real_distribution<> uniform(0, 1); 
	size_t ctr = 0; 	
	double kin_energy, MH, U1, U0;   // kin. energies are necessary for MH criterion
	params theta;
	forces force; 

	for(size_t i=1; i<=N; ++i){
		theta = theta_curr;
		force = force_curr;		
		
		kin_energy = 0;

		for(size_t j=0; j<L; ++j){
			// L OBABO steps			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;
		
			kin_energy -= 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);
		
			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B, sig0);   // HERE FORCES!!!!

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			kin_energy += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);			
			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  		// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;				       
		}
		
		// MH criterion
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata, sig0);					//HERE UPOTS!!!
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata, sig0);
		MH = exp( (-1/T) * (U1 - U0 + kin_energy) );					
		
		if( uniform(twister) < min(1., MH) ){ 												// ACCEPT SAMPLE
			if(i%n_meas == 0 ) {
				meas.mu1[k] = theta.mu1;			
				meas.mu2[k] = theta.mu2;		
				meas.acceptP[k] = min(1., MH);	
				++k;		
			}			

			theta_curr = theta;
			force_curr = force;
			theta_curr.p1 = SF=="A" ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="A" ? -1*theta_curr.p2 : theta_curr.p2;			

         ctr += 1;
        
		}
		else{ 																				 // REJECT SAMPLE		

			if(i%n_meas == 0 ) {
				meas.mu1[k] = theta_curr.mu1;			
				meas.mu2[k] = theta_curr.mu2;	
				meas.acceptP[k] = min(1., MH);	
				++k;		
			}	
			
			theta_curr.p1 = SF=="R" ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="R" ? -1*theta_curr.p2 : theta_curr.p2;					
		}
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return meas;

}



measurement OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0){
   
   cout<<"Starting OMBABO + SF" + SF + " simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 
	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B, sig0);	    // HERE FORCES!!!!

	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};		
//	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};
	meas.mu1[0] = theta_curr.mu1;
	meas.mu2[0] = theta_curr.mu2;
	int k=1;

   double a = exp(-1*gamma*h);      
	
	double Rn1, Rn2;
	normal_distribution<> normal{0,1};
	uniform_real_distribution<> uniform(0, 1); 
	size_t ctr = 0; 	
	double MH, U1, U0, K1, K0;   // pot. and kin. energies for MH criterion
	params theta;
	forces force; 

	for(size_t i=1; i<=N; ++i){
		
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta_curr.p1 = sqrt(a)*theta_curr.p1 + sqrt((1-a)*T)*Rn1;  				// O step
		theta_curr.p2 = sqrt(a)*theta_curr.p2 + sqrt((1-a)*T)*Rn2;			
		
		theta = theta_curr;
		force = force_curr;

		for(size_t j=0; j<L; ++j){
			// L BAB steps					
			
			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B, sig0);   // HERE FORCES!!!!

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			
			    
		
		}
		
		// MH criterion
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata, sig0);	//HERE UPOTS!!!
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata, sig0);
		K0 = 0.5*(theta_curr.p1*theta_curr.p1 + theta_curr.p2*theta_curr.p2);	
		K1 = 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2); 	
		
		MH = exp( (-1/T) * (U1+K1 - (U0+K0)) );					
		
		if( uniform(twister) < min(1., MH) ){ 						// ACCEPT SAMPLE

			if(i%n_meas == 0 ) {
				meas.mu1[k] = theta.mu1;			
				meas.mu2[k] = theta.mu2;		
				meas.acceptP[k] = min(1., MH);	
				++k;			
			}			
	
			theta_curr = theta;
			force_curr = force;			
		
			theta_curr.p1 = SF=="A"? -1*theta_curr.p1 : theta_curr.p1;  // sign flip (SF) 
			theta_curr.p2 = SF=="A"? -1*theta_curr.p2 : theta_curr.p2;

         ctr += 1;
		
		}
		else{  // REJECT SAMPLE
			
			if(i%n_meas == 0 ) {
				meas.mu1[k] = theta_curr.mu1;			
				meas.mu2[k] = theta_curr.mu2;	
				meas.acceptP[k] = min(1., MH);	
				++k;		
			}	
						
			theta_curr.p1 = SF=="R" ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF)  
			theta_curr.p2 = SF=="R" ? -1*theta_curr.p2 : theta_curr.p2;
			 
		
		}
		
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta_curr.p1 = sqrt(a)*theta_curr.p1 + sqrt((1-a)*T)*Rn1;  				// O step
		theta_curr.p2 = sqrt(a)*theta_curr.p2 + sqrt((1-a)*T)*Rn2;		
		
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return meas;

}


double likelihood_GM(const params& theta, const double sig1, const double sig2, const double a1, const double a2, const double x){
	double e1 = exp( -1*(x-theta.mu1)*(x-theta.mu1)/(2*sig1*sig1) );
	double e2 = exp( -1*(x-theta.mu2)*(x-theta.mu2)/(2*sig2*sig2) );
	double p = a1/(sqrt(2*PI)*sig1) * e1 + a2/(sqrt(2*PI)*sig2) * e2;	
	return p;
}

double U_pot_GM( const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  const vector <double>& Xdata, const double sig0){
	double U = 0;
	for(int i=0; i<Xdata.size(); ++i){
		U += log( likelihood_GM(theta, sig1, sig2, a1, a2, Xdata[i]) );
	}
	U -= (theta.mu1*theta.mu1 + theta.mu2*theta.mu2)/(2*sig0*sig0) + log(2*PI*sig0*sig0);
	return -1*U;					  
}


forces get_noisy_force_GM(const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  			  vector <double>& Xdata, const size_t B, const double sig0){
	forces F{0,0};
	double P, e1, e2, x, scale;
	scale = Xdata.size()/double(B);

	if(B<Xdata.size()) shuffle( Xdata.begin(), Xdata.end(), twister );
	
	for(int i=0; i<B; ++i){
			x = Xdata[i];
			P = likelihood_GM(theta, sig1, sig2, a1, a2, x);		
			e1 = exp( -1*(x-theta.mu1)*(x-theta.mu1)/(2*sig1*sig1) );
			e2 = exp( -1*(x-theta.mu2)*(x-theta.mu2)/(2*sig2*sig2) );
			F.fmu1 += 1/P * e1 * (x-theta.mu1);
			F.fmu2 += 1/P * e2 * (x-theta.mu2);
	}


	F.fmu1 *= a1/(sqrt(2*PI)*sig1*sig1*sig1) * scale;
	F.fmu2 *= a2/(sqrt(2*PI)*sig2*sig2*sig2) * scale;
	
	F.fmu1 -= theta.mu1/(sig0*sig0);   // prior part
	F.fmu2 -= theta.mu2/(sig0*sig0);

	return F;
}

vector <double> read_dataset(string datafile){
	ifstream datasource(datafile);
	vector <double> Xdata;
	string row;
	while (getline(datasource, row)){
		Xdata.push_back(stod(row));
	}
return Xdata;
}

