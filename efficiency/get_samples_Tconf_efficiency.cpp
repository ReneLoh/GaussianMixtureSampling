/* 
This script samples the posterior density of the means of two independent Gaussians used to build a 1D Gaussian mixture model.
It uses Metropolized Langevin integrators (given by the "OBABO" scheme) and uses full gradients or stochastic gradients obtained by data subsampling.
Using K processors, it draws K trajectories simultaneously, and averages over them (also in time) before printing out the results.

This version tracks the configurational temperature as an observable, and the Metropolis acceptance probability in time.
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
vector <double> Tconf, acceptP;
};

double likelihood_GM(const params& theta, const double x, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2);
double U_pot_GM(const params& theta, const vector <double>& Xdata, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2, 
					 const double two_sigsig0, const double log_PI_two_sigsig0);
forces get_noisy_force_GM(const params& theta, vector <double>& Xdata, const size_t B, vector <int>& idx_arr, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2, 
					  			  const double F_scale_1, const double F_scale_2, const double sigsig0); 
measurement OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	
								const size_t B, const double sig1, const double sig2, const double a1, const double a2, const size_t n_meas, const double sig0);
measurement MOBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	const size_t B, 
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0);
measurement OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	const size_t B, 
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0);

//double costfactor_pot_vs_grad(const size_t Ntrials, const double sig1, const double sig2, const double a1, const double a2, const double sig0, vector <double>& Xdata);

vector <double> read_dataset(string datafile);


int main(int argc, char *argv[]){

double sig1 = 3;			// GM params (maybe make them global as 
double sig2 = 0.5;		// method signatures are a mess)
double a1 = 0.8;
double a2 = 0.2;

double sig0 = 5;			// Gauss. prior std.dev.

params param0{-4,3,0,0};	// initial conditions


double T = 5;
double gamma = 1;
string datafile = "GM_data_500.csv";
bool tavg = true;								// time average after ensemble average?
int n = 5e5; 									// time average over last n values
int ndist = 10000;							// write out any ndist-th result entry to file 


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

/*// measure the factor between computing time for potential energy and full gradient evaluation
// only needed for Metropolized schemes
double q_avg, time_scale;
int M;
if ( method != 1 ){
	if ( rank == 0 ) cout << "Measuring q..." << endl;
	size_t N_trials = 10000000;	
	double q = costfactor_pot_vs_grad(N_trials, sig1, sig2, a1, a2, sig0, Xdata);
	cout<<"THIS IS RANK "<<rank<<" AND MY Q IS: "<<q<<endl;
	// average q and compute no. of samples to obtain for a run time equivalent to N_0 (input argument) full gradient evals.
	MPI_Reduce(&q, &q_avg, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	if ( rank == 0){
		q_avg /= nr_proc;
		time_scale = q_avg + float(B)/Xdata.size();	// execution time estimate per sample 
																	// in units of a single full gradient evaluation
		M = min( (int)(N_0 / time_scale), (int)2e7 ); // minimum ensures that sample size M does not lead to RAM issues
	}	
	MPI_Bcast(&M, 1, MPI_INT, 0, comm);

	if ( rank == 0) cout << "q = " << q_avg << ", time_scale = " << time_scale << ", M = " << M << endl; 
}*/
int M = N;	// delete if costfactor active or if it was completely removed

string label;
measurement results; // stores Tconf and acceptance probs.

if ( method == 1 ){
	results = OBABO_simu(param0, N_0, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, n_meas, sig0);
	label = "OBABO";
//	time_scale = float(B)/Xdata.size(); 			
}
else if ( method == 2 ){
	results = MOBABO_simu(param0, M, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas, sig0);
	label = "MOBABO_SF" + SF +"_L" + to_string(L);
}
else if ( method == 3 ){
	results = OMBABO_simu(param0, M, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas, sig0);
	label = "OMBABO_SF" + SF +"_L" + to_string(L);
}
else {
	cout << "No valid method number" << endl;
	return 0;
}

cout << "Rank " << rank << " reached barrier" << endl;
MPI_Barrier(comm);


// average over results of different processors and print out file
measurement results_avg;
if( rank==0 ){
	 results_avg.Tconf.resize(results.Tconf.size());  // do only on rank 0 to save RAM 
	 results_avg.acceptP.resize(results.acceptP.size());	 
}	 

MPI_Reduce(&results.Tconf[0], &results_avg.Tconf[0], results.Tconf.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
MPI_Reduce(&results.acceptP[0], &results_avg.acceptP[0], results.acceptP.size(), MPI_DOUBLE, MPI_SUM, 0, comm);

for (int i=0; i<results_avg.Tconf.size(); ++i){
	results_avg.Tconf[i] /= nr_proc;
	results_avg.acceptP[i] /= nr_proc;
}

if( rank==0 ){	

	// time-average results
	if( tavg == true ){	
		cout << "Time averaging...\n"; 
		
		for(int i=results_avg.Tconf.size()-1; i>=0; i-=ndist){
			if(i<=n-1) n=i;
			for(int j=i-n; j<i; ++j){
				results_avg.Tconf[i] += results_avg.Tconf[j];
				results_avg.acceptP[i] += results_avg.acceptP[j];					
			}
				results_avg.Tconf[i] /= n+1 ;
				results_avg.acceptP[i] /= n+1 ; 				
		}	
		
	}			

	
	// Print results
	stringstream stream, stream2;
	string final_label;
	stream << std::fixed << std::setprecision(3) << h;
	if(B==Xdata.size()) final_label = "Tconf_" + label + "_h"+stream.str()+"_avg"+to_string(nr_proc);
	else 					  final_label = "Tconf_" + label + "_h"+stream.str()+"_gradnoiseB"+to_string(B)+"_avg"+to_string(nr_proc);			 

	ofstream file {final_label};
	cout<<"Writing to file...\n";
	for(int i=0; i<results_avg.Tconf.size(); i+=ndist){
		file << i*n_meas << " " << results_avg.Tconf.at(i) << " " << results_avg.acceptP.at(i) << "\n";	
//		file << i*n_meas*time_scale << " " << results_avg.Tconf.at(i) << " " << results_avg.acceptP.at(i) << "\n";
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

	vector <int> idx_arr(Xdata.size());
	for (int i=0; i<Xdata.size(); ++i){		// list of indices, used for subsampling in stochastic gradient comp.
		idx_arr[i] = i;
	}


	// some constants for the integrator, force, potential...
   const double a = exp(-1*gamma*h);    // used for integrator
   const double sqrt_a = sqrt(a);
   const double sqrt_aT = sqrt((1-a)*T);
   const double h_half = 0.5*h;      
   
   const double two_sigsig1 = 2*sig1*sig1;    // used in likelihood
   const double two_sigsig2 = 2*sig2*sig2;
   const double pref_exp1 = a1/(sqrt(2*PI)*sig1);
   const double pref_exp2 = a2/(sqrt(2*PI)*sig2);
				
	const double F_scale_1 = a1/(sqrt(2*PI)*sig1*sig1*sig1);  // used in force
	const double F_scale_2 = a2/(sqrt(2*PI)*sig2*sig2*sig2);
	const double sigsig0 = sig0*sig0;	
	const size_t Xdata_size = Xdata.size();

   params theta = param0; 
	forces force = get_noisy_force_GM(theta, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 					// HERE FORCES!!!!
					  			  				F_scale_1, F_scale_2, sigsig0);    					
	forces full_force = get_noisy_force_GM(theta, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	
					  			  						F_scale_1, F_scale_2, sigsig0);																		// to compute Tconf				
	
	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1, 1)};		

	meas.Tconf[0] = -0.5 * (theta.mu1*full_force.fmu1 + theta.mu2*full_force.fmu2);

	int k=1;
    
	double Rn1, Rn2;
	normal_distribution<> normal{0,1}; 	
	for(size_t i=1; i<=N; ++i){
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta.p1 = sqrt_a*theta.p1 + sqrt_aT*Rn1 + h_half*force.fmu1;  // O+B step		
		theta.p2 = sqrt_a*theta.p2 + sqrt_aT*Rn2 + h_half*force.fmu2;		
														 
		theta.mu1 += h*theta.p1;															// A step							
		theta.mu2 += h*theta.p2;		
  
		force = get_noisy_force_GM(theta, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 					// HERE FORCES!!!!
					  			  				F_scale_1, F_scale_2, sigsig0); 

		theta.p1 += h_half*force.fmu1;													// B step
		theta.p2 += h_half*force.fmu2;		

		Rn1 = normal(twister);
		Rn2 = normal(twister);								 
		theta.p1 = sqrt_a*theta.p1 + sqrt_aT*Rn1;							 // O step
		theta.p2 = sqrt_a*theta.p2 + sqrt_aT*Rn2;

		if(i%n_meas == 0 ) {
			full_force = get_noisy_force_GM(theta, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	
					  			  						F_scale_1, F_scale_2, sigsig0); 	// HERE FORCES!!!!
			meas.Tconf[k] = -0.5 * (theta.mu1*full_force.fmu1 + theta.mu2*full_force.fmu2);				
			++k;		
		}
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;	
	}

	
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;
	    
   return meas;

}



measurement MOBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0){
   
   cout<<"Starting MOBABO + SF" + SF + " simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

	vector <int> idx_arr(Xdata.size());
	for (int i=0; i<Xdata.size(); ++i){		// list of indices, used for subsampling in stochastic gradient comp.
		idx_arr[i] = i;
	}


	// some constants for the integrator, force, potential...
   const double a = exp(-1*gamma*h);    // used for integrator
   const double sqrt_a = sqrt(a);
   const double sqrt_aT = sqrt((1-a)*T);
   const double h_half = 0.5*h;      
   
   const double two_sigsig1 = 2*sig1*sig1;    // used in likelihood
   const double two_sigsig2 = 2*sig2*sig2;
   const double pref_exp1 = a1/(sqrt(2*PI)*sig1);
   const double pref_exp2 = a2/(sqrt(2*PI)*sig2);
	
	const double two_sigsig0 = 2*sig0*sig0;    					// used in U
	const double log_PI_two_sigsig0 = log(PI*two_sigsig0);
				
	const double F_scale_1 = a1/(sqrt(2*PI)*sig1*sig1*sig1);  // used in force
	const double F_scale_2 = a2/(sqrt(2*PI)*sig2*sig2*sig2);
	const double sigsig0 = sig0*sig0;	
	const size_t Xdata_size = Xdata.size();


   params theta_curr = param0; 

	forces force_curr = get_noisy_force_GM(theta_curr, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	// HERE FORCES!!!!
					  			  						F_scale_1, F_scale_2, sigsig0);	    													
	forces full_force = get_noisy_force_GM(theta_curr, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	
					  			  						F_scale_1, F_scale_2, sigsig0);																				// to compute Tconf		
	
	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1, 1)};		

	meas.Tconf[0] = -0.5 * (theta_curr.mu1*full_force.fmu1 + theta_curr.mu2*full_force.fmu2);
	int k=1;

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
			theta.p1 = sqrt_a*theta.p1 + sqrt_aT*Rn1;  				// O step
			theta.p2 = sqrt_a*theta.p2 + sqrt_aT*Rn2;
		
			kin_energy -= 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);
		
			theta.p1 += h_half*force.fmu1;						// B step
			theta.p2 += h_half*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(theta, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	// HERE FORCES!!!!
					  			  				F_scale_1, F_scale_2, sigsig0);   

			theta.p1 += h_half*force.fmu1;						// B step
			theta.p2 += h_half*force.fmu2;			

			kin_energy += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);			
			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt_a*theta.p1 + sqrt_aT*Rn1;  		// O step
			theta.p2 = sqrt_a*theta.p2 + sqrt_aT*Rn2;				       
		}
		
		// MH criterion
		U1 = U_pot_GM(theta, Xdata, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, two_sigsig0, log_PI_two_sigsig0);					//HERE UPOTS!!!  		
		U0 = U_pot_GM(theta_curr, Xdata, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, two_sigsig0, log_PI_two_sigsig0);
		MH = exp( (-1/T) * (U1 - U0 + kin_energy) );					
		
		if( uniform(twister) < min(1., MH) ){ 												// ACCEPT SAMPLE
			if(i%n_meas == 0 ) {
				full_force = get_noisy_force_GM(theta, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	// HERE FORCES!!!!
					  			  						   F_scale_1, F_scale_2, sigsig0);					
				meas.Tconf[k] = -0.5 * (theta.mu1*full_force.fmu1 + theta.mu2*full_force.fmu2);				
				meas.acceptP[k] = min(1., MH);	
				++k;		
			}			

			theta_curr = theta;
			force_curr = force;
			theta_curr.p1 = SF=="A" ? -theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="A" ? -theta_curr.p2 : theta_curr.p2;			

         ctr += 1;
        
		}
		else{ 																				 // REJECT SAMPLE		

			if(i%n_meas == 0 ) {
				meas.Tconf[k] = meas.Tconf[k-1];			
				meas.acceptP[k] = min(1., MH);	
				++k;		
			}	
			
			theta_curr.p1 = SF=="R" ? -theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF=="R" ? -theta_curr.p2 : theta_curr.p2;					
		}
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return meas;

}



measurement OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	const size_t B, 
								const double sig1, const double sig2, const double a1, const double a2, const size_t L, const string SF, const size_t n_meas, const double sig0){
   
   cout<<"Starting OMBABO + SF" + SF + " simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

	vector <int> idx_arr(Xdata.size());
	for (int i=0; i<Xdata.size(); ++i){		// list of indices, used for subsampling in stochastic gradient comp.
		idx_arr[i] = i;
	}


	// some constants for the integrator, force, potential...
   const double a = exp(-1*gamma*h);    // used for integrator
   const double sqrt_a = sqrt(a);
   const double sqrt_aT = sqrt((1-a)*T);
   const double h_half = 0.5*h;      
   
   const double two_sigsig1 = 2*sig1*sig1;    // used in likelihood
   const double two_sigsig2 = 2*sig2*sig2;
   const double pref_exp1 = a1/(sqrt(2*PI)*sig1);
   const double pref_exp2 = a2/(sqrt(2*PI)*sig2);
	
	const double two_sigsig0 = 2*sig0*sig0;    					// used in U
	const double log_PI_two_sigsig0 = log(PI*two_sigsig0);
				
	const double F_scale_1 = a1/(sqrt(2*PI)*sig1*sig1*sig1);  // used in force
	const double F_scale_2 = a2/(sqrt(2*PI)*sig2*sig2*sig2);
	const double sigsig0 = sig0*sig0;	
	const size_t Xdata_size = Xdata.size();



   params theta_curr = param0; 
	
	forces force_curr = get_noisy_force_GM(theta_curr, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 				// HERE FORCES!!!!
					  			  						F_scale_1, F_scale_2, sigsig0);	    													
	forces full_force = get_noisy_force_GM(theta_curr, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	
					  			  						F_scale_1, F_scale_2, sigsig0);																				// to compute Tconf		
	
	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1, 1)};		

	meas.Tconf[0] = -0.5 * (theta_curr.mu1*full_force.fmu1 + theta_curr.mu2*full_force.fmu2);
	int k=1;

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
		theta_curr.p1 = sqrt_a*theta_curr.p1 + sqrt_aT*Rn1;  				// O step
		theta_curr.p2 = sqrt_a*theta_curr.p2 + sqrt_aT*Rn2;
		
		theta = theta_curr;
		force = force_curr;

		for(size_t j=0; j<L; ++j){
			// L BAB steps					
			
			theta.p1 += h_half*force.fmu1;						// B step
			theta.p2 += h_half*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(theta, Xdata, B, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	// HERE FORCES!!!!
					  			  				F_scale_1, F_scale_2, sigsig0);   


			theta.p1 += h_half*force.fmu1;						// B step
			theta.p2 += h_half*force.fmu2;			
			    
		
		}
		
		// MH criterion
		U1 = U_pot_GM(theta, Xdata, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, two_sigsig0, log_PI_two_sigsig0);					//HERE UPOTS!!!  		
		U0 = U_pot_GM(theta_curr, Xdata, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, two_sigsig0, log_PI_two_sigsig0);
		K0 = 0.5*(theta_curr.p1*theta_curr.p1 + theta_curr.p2*theta_curr.p2);	
		K1 = 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2); 	
		
		MH = exp( (-1/T) * (U1+K1 - (U0+K0)) );					
		
		if( uniform(twister) < min(1., MH) ){ 						// ACCEPT SAMPLE

			if(i%n_meas == 0 ) {
				full_force = get_noisy_force_GM(theta, Xdata, Xdata_size, idx_arr, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2, 	// HERE FORCES!!!!
					  			  						   F_scale_1, F_scale_2, sigsig0);						
				meas.Tconf[k] = -0.5 * (theta.mu1*full_force.fmu1 + theta.mu2*full_force.fmu2);				
				meas.acceptP[k] = min(1., MH);	
				++k;			
			}			
	
			theta_curr = theta;
			force_curr = force;			
		
			theta_curr.p1 = SF=="A"? -theta_curr.p1 : theta_curr.p1;  // sign flip (SF) 
			theta_curr.p2 = SF=="A"? -theta_curr.p2 : theta_curr.p2;

         ctr += 1;
		
		}
		else{  // REJECT SAMPLE
			
			if(i%n_meas == 0 ) {
				meas.Tconf[k] = meas.Tconf[k-1];			
				meas.acceptP[k] = min(1., MH);	
				++k;			
			}	
						
			theta_curr.p1 = SF=="R" ? -theta_curr.p1 : theta_curr.p1;		// sign flip (SF)  
			theta_curr.p2 = SF=="R" ? -theta_curr.p2 : theta_curr.p2;
			 
		
		}
		
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta_curr.p1 = sqrt_a*theta_curr.p1 + sqrt_aT*Rn1;  				// O step
		theta_curr.p2 = sqrt_a*theta_curr.p2 + sqrt_aT*Rn2;
		
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return meas;

}


// This is the likelihood evaluation of a given data point x. It is used by both force and energy.
double likelihood_GM(const params& theta, const double x, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2){
	double e1 = exp( -1*(x-theta.mu1)*(x-theta.mu1)/(two_sigsig1) );
	double e2 = exp( -1*(x-theta.mu2)*(x-theta.mu2)/(two_sigsig2) );
	double p = pref_exp1 * e1  +  pref_exp2 * e2;	
	return p;
}



// This is the energy routine, it needs to process the whole vector Xdata
double U_pot_GM(const params& theta, const vector <double>& Xdata, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2, 
					 const double two_sigsig0, const double log_PI_two_sigsig0){
	
	double U = 0;
	
	for(int i=0; i<Xdata.size(); ++i){

		U += log( likelihood_GM(theta, Xdata[i], two_sigsig1, two_sigsig2, pref_exp1, pref_exp2) );	// sum over log-likelihoods

	}

	U -= (theta.mu1*theta.mu1 + theta.mu2*theta.mu2)/two_sigsig0 + log_PI_two_sigsig0;  // log-prior part.

	return -U;					  
}



// This is the force routine, it needs to process B random elements of the vector Xdata
forces get_noisy_force_GM(const params& theta, vector <double>& Xdata, const size_t B, vector <int>& idx_arr, const double two_sigsig1, const double two_sigsig2, const double pref_exp1, const double pref_exp2, 
					  			  const double F_scale_1, const double F_scale_2, const double sigsig0){
	forces F{0,0};
	double P, e1, e2, x, scale;
	int help_int, idx;
	int size_minus_B = Xdata.size()-B;
	scale = Xdata.size() / double(B);
	
	if(Xdata.size() != B){	// in case of subsampling...
		
		// idx_arr stores the possible indices of the vector Xdata
		// this loop randomly chooses B of them and stores them
		// at the end of idx_arr.
		for(int i=Xdata.size()-1; i>=size_minus_B; --i){ 

			uniform_int_distribution<> distrib(0, i);		// recreates this in every iter... is there a better way?
		
			idx = distrib(twister);
			help_int = idx_arr[i];
			idx_arr[i] = idx_arr[idx];
			idx_arr[idx] = help_int; 

		}

	}

	for(int i=idx_arr.size()-1; i>=size_minus_B; --i){		// actual force evaluation. 
																			// the B data points to be considered are given 
																			// by the B last indices stored in idx_arr.
	
			x = Xdata[ idx_arr[i] ];
			P = likelihood_GM(theta, x, two_sigsig1, two_sigsig2, pref_exp1, pref_exp2);		// likelihood of a single data point
			e1 = exp( -(x-theta.mu1)*(x-theta.mu1)/(two_sigsig1) );
			e2 = exp( -(x-theta.mu2)*(x-theta.mu2)/(two_sigsig2) );
			F.fmu1 += 1/P * e1 * (x-theta.mu1);
			F.fmu2 += 1/P * e2 * (x-theta.mu2);
				
	}


	F.fmu1 *= F_scale_1 * scale;
	F.fmu2 *= F_scale_2 * scale;
	
	F.fmu1 -= theta.mu1/(sigsig0);   // prior part
	F.fmu2 -= theta.mu2/(sigsig0);

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



/*double costfactor_pot_vs_grad(const size_t N_trials, const double sig1, const double sig2, const double a1, const double a2, const double sig0, vector <double>& Xdata){

	auto T_exe1 = chrono::high_resolution_clock::now();

	vector <int> idx_arr(Xdata.size());		// needed to pass to force functions, but bears no meaning here
	double q;
	params theta{1,1,1,1};  // dummy param
	
	
	// measure time for full gradient computations
	auto TG1 = chrono::high_resolution_clock::now();
	forces F1, F2{0,0};
	for (size_t k=1; k<N_trials; ++k){
		F1 = get_noisy_force_GM(theta, sig1, sig2, a1, a2, Xdata, Xdata.size(), idx_arr, sig0);
//		F2.fmu1 += 1./k * F1.fmu1;
//		theta.mu1 *= 0.9999;		
	} 

	auto TG2 = chrono::high_resolution_clock::now();
	auto ms1 = chrono::duration_cast<chrono::milliseconds>(TG2 - TG1);
	double TG = ms1.count();
	cout<<"TGRADIENT: "<<TG<<endl;
//	cout<<"Necessary printout of F2.fmu1 in q measurement: "<<F2.fmu1<<endl;	
	
	
	// measure time for potential energy computation
	auto TU1 = chrono::high_resolution_clock::now();
	double U = 0;
	for (size_t k=1; k<N_trials; ++k){
//		U += 1./k * U_pot_GM(theta, sig1, sig2, a1, a2, Xdata, sig0);
		U = U_pot_GM(theta, sig1, sig2, a1, a2, Xdata, sig0);	
//		theta.mu1 *= 0.9999;	
	} 
	auto TU2 = chrono::high_resolution_clock::now();
	auto ms2 = chrono::duration_cast<chrono::milliseconds>(TU2 - TU1);
	double TU = ms2.count(); 
	cout<< "TPOTENTIAL: "<<TU<<endl;
//	cout<<"Necessary printout of U in q measurement: "<<U<<endl;	
	
	auto T_exe2 = chrono::high_resolution_clock::now();
	auto ms3 = chrono::duration_cast<chrono::milliseconds>(T_exe2 - T_exe1);
	cout << "q measurement took " << ms3.count() << " ms"<<endl; 

return TU / TG; 
}
*/