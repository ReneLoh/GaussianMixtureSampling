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


using namespace std;

const double PI = 3.141592653589793;


mt19937 twister;

struct params {
double mu1,mu2,p1,p2;
};
struct forces{
double fmu1, fmu2;
};

double likelihood_GM(const params& theta, const double sig1, const double sig2, const double a1, const double a2, const double x);
double U_pot_GM( const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  const vector <double>& Xdata);
forces get_noisy_force_GM(const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  			  vector <double>& Xdata, const size_t B);
/*vector <vector <size_t>> OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, 
									const double xrange, const double yrange, size_t nrbins, const vector <double>& Xdata,	const double sig1, 
									const double sig2, const double a1, const double a2, const double biasx=0, const double biasy=0);*/
vector <double> OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	
									const size_t B, const double sig1, const double sig2, const double a1, const double a2, const size_t n_meas);
vector <double> OBABO_simu_MH(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	const size_t B, 
										const double sig1, const double sig2, const double a1, const double a2, const size_t L, const bool SF, const size_t n_meas);
vector <double> OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
									 const double sig1, const double sig2, const double a1, const double a2, const size_t L, const bool SF, const size_t n_meas);
void calculate_bins(vector <vector <size_t>>& bin_ctr, const params& theta, const double xrange, const double deltax, 
						  const double yrange, const double deltay, const size_t nrbins);						
vector <double> read_dataset(string datafile);


int main(int argc, char *argv[]){

double sig1 = 3;
double sig2 = 0.5;
double a1 = 0.8;
double a2 = 0.2;

params param0{2,3,0,0};

double xrange = 1;
double biasx = 5;
double yrange = 1;
double biasy = 0;

double T=0.1;
double gamma = 1;
size_t nrbins = 500;
string datafile = "GM_data.csv";

int method = atoi(argv[1]);   // must be in {1,2,3}
double h = atof(argv[2]);
size_t N_0 = atol(argv[3]);
size_t B = atof(argv[4]);
size_t L = atoi(argv[5]);
bool SF = atoi(argv[6]);
int n_meas = atoi(argv[7]);	// any n_meas steps, store avg. Tconfig until then
int randomseed = atoi(argv[8]);

size_t N = N_0/(L+2);

vector <double> Xdata = read_dataset(datafile);

seed_seq seq{1,20,3200,403,5*randomseed+1,12000,73667,9474+randomseed,19151-randomseed};
vector<std::uint32_t> seeds(1);
seq.generate(seeds.begin(), seeds.end());
twister.seed(seeds.at(0)); 

string label;
//vector <vector <size_t>> bin_ctr;
vector <double> Tconfigs;
if(method==1){
//	bin_ctr = OBABO_simu(param0, N_0, h, T, gamma, xrange, yrange, nrbins, Xdata, sig1, sig2, a1, a2, biasx);
	Tconfigs = OBABO_simu(param0, N_0, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, n_meas);
	label = "OBABO";
}
else if(method==2){
	Tconfigs = OBABO_simu_MH(param0, N, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas);
	label = "OBABO_MH_SF" + to_string(SF) + "_L" + to_string(L);
}
else if (method==3){
	Tconfigs = OMBABO_simu(param0, N, h, T, gamma, Xdata, B, sig1, sig2, a1, a2, L, SF, n_meas);
	label = "OMBABO" + string("_L") + to_string(L);
}
else {
	cout<<"No valid method number"<<endl;
	return 0;
}



/*stringstream stream, stream2;

stream << std::fixed << std::setprecision(3) << h;
ofstream file {"histo_" + label + "_h"+stream.str()};

cout<<"Writing to file...\n";
for(int i=0; i<bin_ctr.size(); ++i){
	for(int j=0; j<bin_ctr.at(0).size(); ++j){	
		file << bin_ctr.at(i).at(j) << " ";
	}
	file << "\n";
}
file.close();*/


stringstream stream, stream2;
string final_label;

stream << std::fixed << std::setprecision(3) << h;
if(B==Xdata.size()) final_label = "TEST_Tconf_" + label + "_h"+stream.str();
else 					  final_label = "Tconf_" + label + "_h"+stream.str()+"_gradnoiseB"+to_string(B);			 

ofstream file {final_label};

cout<<"Writing to file...\n";
for(int i=0; i<Tconfigs.size(); ++i){
	file << i*n_meas << " " << Tconfigs.at(i) << "\n";
}
file.close();

return 0;

}




vector <double> OBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata,	
									const size_t B, const double sig1, const double sig2, const double a1, const double a2, const size_t n_meas){
    
   cout<<"Starting OBABO simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta = param0; 
	forces force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B);    // HERE FORCES!!!!

	double Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2) ;
	vector <double> Tconfigs(0);
	Tconfigs.push_back(Tconfig);
	double Tconfig_sum = Tconfig;

   double a = exp(-1*gamma*h);      

	double Rn1, Rn2;
	normal_distribution<> normal{0,1}; 	
	for(size_t i=1; i<N; ++i){
		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1 + 0.5*h*force.fmu1;  // O+B step		
		theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2 + 0.5*h*force.fmu2;		
														 
		theta.mu1 += h*theta.p1;															// A step							
		theta.mu2 += h*theta.p2;		

		force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B);   // HERE FORCES!!!!

		theta.p1 += 0.5*h*force.fmu1;														 // B step
		theta.p2 += 0.5*h*force.fmu2;		

		Rn1 = normal(twister);
		Rn2 = normal(twister);								 
		theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;							 // O step
		theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;
		
		
		Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2);
		Tconfig_sum += Tconfig;		

		if(i%n_meas ==0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;	
	}
	 
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;
	    
   return Tconfigs;

}


vector <double> OBABO_simu_MH(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
										const double sig1, const double sig2, const double a1, const double a2, const size_t L, const bool SF, const size_t n_meas){
   
   cout<<"Starting OBABO w. MH + SF(" + to_string(SF) + ") simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 
	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B);	    // HERE FORCES!!!!

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
		kin_energy = 0;
		
		force = force_curr;

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
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B);   // HERE FORCES!!!!

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			kin_energy += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);			
			
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;				       
		}
		
		// MH criterion
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata);
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata);
		MH = exp( (-1/T) * (U1 - U0 + kin_energy) );					//HERE UPOTS!!!
		
		if( uniform(twister) < min(1., MH) ){ 								// ACCEPT SAMPLE
			
			theta_curr = theta;
			force_curr = force;
			theta_curr.p1 = SF==true ? -1*theta_curr.p1 : theta_curr.p1;		// sign flip (SF) 		
			theta_curr.p2 = SF==true ? -1*theta_curr.p2 : theta_curr.p2;			

         ctr += 1;
         
			Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2);
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
		}
		else{  // REJECT SAMPLE
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		}
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return Tconfigs;

}


/*vector <double> OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
									 const double sig1, const double sig2, const double a1, const double a2, const size_t L, const bool SF, const size_t n_meas){
   
   cout<<"Starting OMBABO simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 
	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B);	    // HERE FORCES!!!!

	double Tconfig = -1*(force_curr.fmu1 * theta_curr.mu1  +  force_curr.fmu2 * theta_curr.mu2);
	vector <double> Tconfigs(0);
	Tconfigs.push_back(Tconfig);
	double Tconfig_sum = Tconfig;

   double a = exp(-1*gamma*h);      
	
	double Rn1, Rn2;
	normal_distribution<> normal{0,1};
	uniform_real_distribution<> uniform(0, 1); 
	size_t ctr = 0; 	
	double MH, U1, U0, K1, K0;   // pot. and kin. energies for MH criterion
	params theta;
	forces force; 

	for(size_t i=1; i<N; ++i){
		
		theta = theta_curr;
		force = force_curr;

		Rn1 = normal(twister);
		Rn2 = normal(twister);
		theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
		theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;

		for(size_t j=0; j<L; ++j){
			// L BAB steps			
		
			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B);   // HERE FORCES!!!!

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			
			       
		}
		
		// MH criterion
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata);	//HERE UPOTS!!!
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata);
		K0 = 0.5*(theta_curr.p1*theta_curr.p1 + theta_curr.p2*theta_curr.p2);	
		K1 = 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2); 	
		MH = exp( (-1/T) * (U1+K1 - (U0+K0)) );					
		
		if( uniform(twister) < min(1., MH) ){ 						// ACCEPT SAMPLE
			
			theta_curr = theta;
			force_curr = force;			

         ctr += 1;
         
			Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2);
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
		}
		else{  // REJECT SAMPLE
			
			theta_curr.p1 *= -1;		// reverse velocity  
			theta_curr.p2 *= -1;
			 
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
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
	    
   return Tconfigs;

}*/

vector <double> OMBABO_simu(const params param0, const size_t N, const double h, const double T, const double gamma, vector <double>& Xdata, const size_t B,
									 const double sig1, const double sig2, const double a1, const double a2, const size_t L, const bool SF, const size_t n_meas){
   
   cout<<"Starting OMBABO simulation..."<<endl;
	auto t1 = chrono::high_resolution_clock::now();

   params theta_curr = param0; 
	forces force_curr = get_noisy_force_GM(T, theta_curr, sig1, sig2, a1, a2, Xdata, B);	    // HERE FORCES!!!!

	double Tconfig = -1*(force_curr.fmu1 * theta_curr.mu1  +  force_curr.fmu2 * theta_curr.mu2);
	vector <double> Tconfigs(0);
	Tconfigs.push_back(Tconfig);
	double Tconfig_sum = Tconfig;

   double a = exp(-1*gamma*h);      
	
	double Rn1, Rn2;
	normal_distribution<> normal{0,1};
	uniform_real_distribution<> uniform(0, 1); 
	size_t ctr = 0; 	
	double MH, U1, U0, K1, K0;   // pot. and kin. energies for MH criterion
	params theta;
	forces force; 

	for(size_t i=1; i<N; ++i){
		
		theta = theta_curr;
		force = force_curr;
		K0 = 0;
		K1 = 0;

		for(size_t j=0; j<L; ++j){
			// L BAB steps			
		
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;				

//			K0 += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);			
			
			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			

			theta.mu1 += h*theta.p1;							// A step
			theta.mu2 += h*theta.p2;
			
			force = get_noisy_force_GM(T, theta, sig1, sig2, a1, a2, Xdata, B);   // HERE FORCES!!!!

			theta.p1 += 0.5*h*force.fmu1;						// B step
			theta.p2 += 0.5*h*force.fmu2;			
			       
//			K1 += 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2);    
			    
			Rn1 = normal(twister);
			Rn2 = normal(twister);
			theta.p1 = sqrt(a)*theta.p1 + sqrt((1-a)*T)*Rn1;  				// O step
			theta.p2 = sqrt(a)*theta.p2 + sqrt((1-a)*T)*Rn2;		
		
		}
		
		// MH criterion
		U0 = U_pot_GM( T, theta_curr, sig1, sig2, a1, a2, Xdata);	//HERE UPOTS!!!
		U1 = U_pot_GM( T, theta, sig1, sig2, a1, a2, Xdata);
		K0 = 0.5*(theta_curr.p1*theta_curr.p1 + theta_curr.p2*theta_curr.p2);	
		K1 = 0.5*(theta.p1*theta.p1 + theta.p2*theta.p2); 	
		MH = exp( (-1/T) * (U1+K1 - (U0+K0)) );					
		
		if( uniform(twister) < min(1., MH) ){ 						// ACCEPT SAMPLE
			
			theta_curr = theta;
			force_curr = force;			

         ctr += 1;
         
			Tconfig = -1*(force.fmu1 * theta.mu1  +  force.fmu2 * theta.mu2);
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
		}
		else{  // REJECT SAMPLE
			
			theta_curr.p1 *= -1;		// reverse velocity  
			theta_curr.p2 *= -1;
			 
			Tconfig_sum += Tconfig;
			if(i%n_meas == 0 ) Tconfigs.push_back(Tconfig_sum / (i+1));	
		
		}
		
		
		
	
		if(i%int(1e6)==0) cout<<"Iteration "<<i<<" done!"<<endl;
		
	}

	cout <<"Acceptance probability was "<<float(ctr)/N<<endl;
	auto t2 = chrono::high_resolution_clock::now();
	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;	
	    
   return Tconfigs;

}

double likelihood_GM(const params& theta, const double sig1, const double sig2, const double a1, const double a2, const double x){
	double e1 = exp( -1*(x-theta.mu1)*(x-theta.mu1)/(2*sig1*sig1) );
	double e2 = exp( -1*(x-theta.mu2)*(x-theta.mu2)/(2*sig2*sig2) );
	double p = a1/(sqrt(2*PI)*sig1) * e1 + a2/(sqrt(2*PI)*sig2) * e2;	
	return p;
}

double U_pot_GM( const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  const vector <double>& Xdata){
	double U = 0;
	for(int i=0; i<Xdata.size(); ++i){
		U += log( likelihood_GM(theta, sig1, sig2, a1, a2, Xdata[i]) );
	}
	return -1*T*U;					  
}


forces get_noisy_force_GM(const double T, const params& theta, const double sig1, const double sig2, const double a1, const double a2, 
					  			  vector <double>& Xdata, const size_t B){
	forces F{0,0};
	double P, e1, e2, x, scale;
	scale = Xdata.size()/double(B);

	if(B<Xdata.size()) shuffle( Xdata.begin(), Xdata.end(), twister );

	for(int i=0; i<B; ++i){
			x = Xdata[i];			
			P = likelihood_GM(theta, sig1, sig2, a1, a2, x);		
			e1 = exp( -1*(x-theta.mu1)*(x-theta.mu1)/(2*sig1*sig1) );
			e2 = exp( -1*(x-theta.mu2)*(x-theta.mu2)/(2*sig2*sig2) );
			F.fmu1 += 1/P * e1 * (x-theta.mu1)/(sig1*sig1);
			F.fmu2 += 1/P * e2 * (x-theta.mu2)/(sig2*sig2);
	}
	F.fmu1 *= T*a1/(sqrt(2*PI)*sig1) * scale;
	F.fmu2 *= T*a2/(sqrt(2*PI)*sig2) * scale;

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


void calculate_bins(vector <vector <size_t>>& bin_ctr, const params& theta, const double xrange, const double deltax, 
						  const double yrange, const double deltay, const size_t nrbins, const double biasx, const double biasy){
	size_t k1,k2;		
	k1 = 1 + floor( (theta.mu1-biasy + yrange)/deltay );
	k1 = theta.mu1-biasy > yrange ? nrbins+1 : k1;
	k1 = theta.mu1-biasy < -1*yrange ? 0 : k1;		
	
	k2 = 1 + floor( (theta.mu2-biasx + xrange)/deltax );
	k2 = theta.mu2-biasx > xrange ? nrbins+1 : k2;
	k2 = theta.mu2-biasx < -1*xrange ? 0 : k2;		
	bin_ctr[k1][k2]++;
}
