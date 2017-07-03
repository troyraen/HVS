//
//	findBeta.cc
//
//	written by Troy Joseph Raen on 8/5/2015
//
//	returns array of betas representing shells of roughly constant mass
//	for BT potential dM
//  .
//  z as short axis
//  std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<

#include <vector>
#include "/Users/troyraen/dropbox/research/headers/constants.h"


using namespace std;

int N=0;
vector<double> InC, sigC, quuC, LC;


// FUNCTIONS
void setMGEf(void);
double deltaMBT(const vector<double> &qObl, const double &rho, const double &dBeta);
double calcDensity(const vector<double> &qObl);


// MAIN
int main(){
	clock_t t1,t2;
	t1=clock();
	const double halfMNSC= 1.4e7*Msun, a= 4.2*pc, (1-q*q), dMconst= halfMNSC/1.e4, deltaB= 1.e-6;
	double mass= 0., beta= deltaB, shellM= 0., ai= beta*a, rhoai= 0.;
	vector<double> qObl(6,0.); qObl[4]= ecc; qObl[5]= a*ecc; // u, v, phi, ai, e, ae (BT)
	setMGEf();

	while (ai/pc<10.) {
		ai= beta*a;
		qObl[0]= beta; qObl[1]= pi/4.; qObl[2]= pi/4.; qObl[3]= ai;
		
		rhoai= calcDensity(qObl);
		double deltaM= deltaMBT(qObl,rhoai,deltaB);
		mass+= deltaM; shellM+= deltaM;
		
		if (shellM > dMconst/beta) {
			cout << beta << "\t" << shellM/Msun << "\t" << mass/Msun << "\t" << ai/pc << "\t" << rhoai/Msun*pc*pc*pc << endl;
			shellM= 0.;
		}
		
		beta+= deltaB;
	}
	
	cout << beta << "\t" << shellM/Msun << "\t" << mass/Msun << "\t" << ai/pc << "\t" << rhoai/Msun*pc*pc*pc << endl;
	// Runtime
	t2=clock();
	cerr << "\truntime: " << (t2-t1)/CLOCKS_PER_SEC/60. << " minutes (" << (t2-t1)/CLOCKS_PER_SEC << " seconds)" << endl;
	return 0;
}



// ------ deltaM from Binney & Tremaine
double deltaMBT(const vector<double> &qObl, const double &rho, const double &dBeta){
	const double beta= qObl[0], ai= qObl[3], e= qObl[4];
	return 4.*pi*rho*ai*ai*ai*sqrt(1.-e*e)*beta*beta*dBeta;
}




// ------ density from Emsellem 1994, luminosity from Fieldmeier 14
double calcDensity(const vector<double> &qObl){
	double density= 0.;
	const double ae= qObl[3]*qObl[4], R= ae*cosh(qObl[0])*sin(qObl[1]), z= ae*sinh(qObl[0])*cos(qObl[1]), R2= R*R, z2= z*z;
	
	for (int ii=0; ii<N; ii++) { // NSC MGE density
		const double sig2 = sigC[ii]*sigC[ii];
		density += LC[ii]/pow(sig2*2.*pi,3./2.)/quuC[ii] * exp(-.5/sig2*(R2+z2/quuC[ii]/quuC[ii]));
	}
	return density;
}





// ------- set MGE parameters, from Feldmeier 14 -------
void setMGEf(void){
	InC.clear(); sigC.clear(); quuC.clear(); LC.clear();
	const double max=10000.*pc, eC, MoL = 0.56*Msun/Lsun;
	
	double cgs = 1.e5*Lsun/pc/pc;
	InC.push_back(112.*cgs);
	InC.push_back(46.2*cgs);
	InC.push_back(16.0*cgs);
	InC.push_back(20.4*cgs);
	InC.push_back(7.48*cgs);
	InC.push_back(4.53*cgs);
	InC.push_back(0.77*cgs);
	InC.push_back(0.73*cgs);
	InC.push_back(0.47*cgs);
	InC.push_back(0.17*cgs);
	InC.push_back(0.31*cgs);
	N=InC.size();
	
	sigC.push_back(0.1*asec);
	sigC.push_back(2.1*asec);
	sigC.push_back(8.4*asec);
	sigC.push_back(11.5*asec);
	sigC.push_back(22.8*asec);
	sigC.push_back(66.4*asec);
	sigC.push_back(143.*asec);
	sigC.push_back(184.*asec);
	sigC.push_back(581.*asec);
	sigC.push_back(2705.*asec);
	sigC.push_back(2705.*asec);
	
	quuC.push_back(0.9);
	quuC.push_back(1.0);
	quuC.push_back(0.6);
	quuC.push_back(1.0);
	quuC.push_back(0.7);
	quuC.push_back(0.7);
	quuC.push_back(1.0);
	quuC.push_back(0.2);
	quuC.push_back(0.4);
	quuC.push_back(1.0);
	quuC.push_back(0.2);
	
	for (int ii=0; ii<N; ii++) {
		LC.push_back(MoL*2.*pi*InC[ii]*sigC[ii]*sigC[ii]*Q*(1-exp(-0.5/sigC[ii]/sigC[ii]*max*max)));
		
	}
}

