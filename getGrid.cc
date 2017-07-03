//
//	getGrid.cc
//
//	written by Troy Joseph Raen on 8/12/2015
//
//	returns array of [R, z, rho, phi (potential), acc (acceleration)]
//	for MGE model (including NSC and bulge, see plotDensity folder
//	for more info on the combination).
//
//  z as short axis
//  std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<

#include <vector>
#include <complex>
#include "/Users/troyraen/dropbox/research/headers/constants.h"


using namespace std;


vector< vector<int> > jGroups;
vector<double> In, sig, quu, gam;


// FUNCTIONS
void setMGE(void);
void transCtoObl(const vector<double> &qCart, vector<double> &qObl, const double &ae);
double calcPhiAcc(const vector<double> &cart, const double &dM, vector<double> &acc, const double &beta);
double calcDensity(const vector<double> &coords, const int &j);


// MAIN
int main(){
	clock_t t1,t2;
	t1=clock();
	setMGE();
	
	double max= 5.*pc, min= 1.e-2*pc, deltaCoords= 0.1*pc; // cartesian coordinates
	vector<double> coords(7,0.); // [0:R, 1:z, 2:R^2, 3:z^2, 4:qj, 5:a, 6:ecc]
	double a= 4.2*pc, dBeta= 1.e-1; coords[5]= a;
	
	for (coords[0]= min; coords[0]<max; coords[0]+= deltaCoords) {
		coords[2]= coords[0]*coords[0];
		for (coords[1]= min; coords[1]<max; coords[1]+= deltaCoords) {
			coords[3]= coords[1]*coords[1];
			double rhoTot= 0., phiTot= 0.; vector<double> accTot(2,0.);
			
			for (int j=0; j<jGroups.size(); j++) {
				int qInd= jGroups[j][0]; double qj= quu[qInd]; coords[4]= qj;
				cout << qj << endl;
				double ecc= sqrt(1-qj*qj); coords[6]= ecc;
				
				double beta= dBeta;
				while (beta<2.) {
					
					double rhok= calcDensity(coords,j);
					double dMass= 4.*pi*rhok* a*a*a*sqrt(1-ecc*ecc)* beta*beta*dBeta; // deltaM from Binney & Tremaine
					
					vector<double> acck(2,0.);
					double phik= calcPhiAcc(coords,dMass,acck,beta);
					
					phiTot+= phik; accTot[0]+= acck[0]; accTot[1]+= acck[1];
					
					beta+= dBeta;
				}
			}
			
			// check density
			for (int j=0; j<jGroups.size(); j++) {
				int qInd= jGroups[j][0]; double qj= quu[qInd]; coords[4]= qj;
				rhoTot+= calcDensity(coords,j);
			}
			
			cout << coords[0]/pc << "\t" << coords[1]/pc << "\t" << rhoTot/Msun*pc*pc*pc << "\t" << phiTot << "\t" << accTot[0]/1.e5 << "\t" << accTot[1]/1.e5 << endl;
		}
	}
	

	
	
	// Runtime
	t2=clock();
	cerr << "\truntime: " << (t2-t1)/CLOCKS_PER_SEC/60. << " minutes (" << (t2-t1)/CLOCKS_PER_SEC << " seconds)\n" << endl;
	return 0;
}




//-------- potential and acceleration from BT
double calcPhiAcc(const vector<double> &cart, const double &dM, vector<double> &acc, const double &beta){
	double phi=0.;
	const double u= acosh(beta/cart[6]), ae= cart[5]*cart[6], pre= -G*dM/ae;
	vector<double> obl(2,0.); transCtoObl(cart,obl,ae);
	
	if (u< obl[0]){
		double sechu0= 1./cosh(obl[0]), coshu02= cosh(obl[0])*cosh(obl[0]),
		cosv02= cos(obl[1])*cos(obl[1]);
		
		// calculate potential
		phi= pre*asin(sechu0);
		cout << "inside: " << phi << endl;
		
		// calculate acceleration
		double accObl= pre*sechu0*tanh(obl[0])/ sqrt((1-sechu0*sechu0)* (coshu02+cosv02-1)); // in u direction, then translate to R,z:
		double C= 1./sqrt(sinh(obl[0])*sinh(obl[0])+ sin(obl[1])*sin(obl[1]));
		acc[0]= accObl*C* sinh(obl[0])*cos(obl[1]);
		acc[1]= accObl*C* cosh(obl[0])*sin(obl[1]);
	}
	
	else { phi= pre*asin(1./cosh(u));
		cout << "outside: " << phi << "\t" << cart[6] << endl;
	}

	return phi;
}



// ------ cartesian to oblate spheroidal coord transform (short axis is z)
void transCtoObl(const vector<double> &qCart, vector<double> &qObl, const double &ae){
	complex<double> arg(qCart[0]/ae, qCart[1]/ae);
	qObl[0] = real(acosh(arg)); // u
	qObl[1] = imag(acosh(arg)); // v
}


// ------ density from Emsellem 1994
double calcDensity(const vector<double> &coords, const int &j){
	int size= jGroups[j].size();
	double density= 0.;
	
	for (int ii= 0; ii<size; ii++) {
		int ind= jGroups[j][ii];
		const double sig2 = sig[ind]*sig[ind];
		density += gam[ind]/pow(sig2*2.*pi,3./2.)/coords[4] * exp(-.5/sig2*(coords[2]+coords[3]/coords[4]/coords[4]));
	}
	
	return density;
}



// ------- set MGE parameters, from Feldmeier 14 + Hernquist->MGE for bulge -------
void setMGE(void){
	In.clear(); sig.clear(); quu.clear(); gam.clear();
	const double max=10000.*pc, eps=0.29, Q=sqrt(1-eps*eps), MoL = 0.56*Msun/Lsun;
	
	double cgs = 1.e5*Lsun/pc/pc;
	In.push_back(112.*cgs);
	In.push_back(46.2*cgs);
	In.push_back(16.0*cgs);
	In.push_back(20.4*cgs);
	In.push_back(7.48*cgs);
	In.push_back(4.53*cgs);
	In.push_back(0.77*cgs);
	In.push_back(0.73*cgs);
	In.push_back(0.47*cgs);
	In.push_back(0.17*cgs);
	In.push_back(0.31*cgs);
	
	sig.push_back(0.1*asec);
	sig.push_back(2.1*asec);
	sig.push_back(8.4*asec);
	sig.push_back(11.5*asec);
	sig.push_back(22.8*asec);
	sig.push_back(66.4*asec);
	sig.push_back(143.*asec);
	sig.push_back(184.*asec);
	sig.push_back(581.*asec);
	sig.push_back(2705.*asec);
	sig.push_back(2705.*asec);
	
	quu.push_back(0.9);
	quu.push_back(1.0);
	quu.push_back(0.6);
	quu.push_back(1.0);
	quu.push_back(0.7);
	quu.push_back(0.7);
	quu.push_back(1.0);
	quu.push_back(0.2);
	quu.push_back(0.4);
	quu.push_back(1.0);
	quu.push_back(0.2);
	
	for (int ii=0; ii<In.size(); ii++) {
		gam.push_back(MoL*2.*pi*In[ii]*sig[ii]*sig[ii]*Q*(1-exp(-0.5/sig[ii]/sig[ii]*max*max)));
	}
	
	vector<int> j02, j04, j06, j07, j09, j10;
	j02.push_back(7); j02.push_back(10);
	j04.push_back(8);
	j06.push_back(2);
	j07.push_back(4); j07.push_back(5);
	j09.push_back(0);
	j10.push_back(1); j10.push_back(3); j10.push_back(6); j10.push_back(9);
	
	jGroups.push_back(j02); jGroups.push_back(j04); jGroups.push_back(j06); jGroups.push_back(j07); jGroups.push_back(j09); jGroups.push_back(j10);
	
}
