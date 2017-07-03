
//
//  main.cpp
//  HighVelocityStars
//
//  Created by Troy Joseph Raen on 1/27/15.
//
//  short axis is z

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>	// srand 48, rand
#include "constants.h"
#include "ctemplates.h"
#include "randx.h"
#include "star.h"
#include "nsc_potential.h"


using namespace std;


// load beta (const mass) vector for NSC potential
vector<double> betaShells; double mShell; int bsize; void loadBeta();

// other functions
void vec3Reset(vector<double> &reset, const vector<double> &val);
void leap(vector<double> &q, vector<double> &vel, vector<double> &acc, const double dt);
inline double doubleDt(double dt);
bool errCheck(vector<double> &q, vector<double> &qhi, vector<double> &vel, const double startE);
void calcAcc(const vector<double> &q, vector<double> &acc);
double calcE(const vector<double> &q, const vector<double> &vel);
double calcR(const vector<double> &q, const bool ind);
	
const int potential = 3;

//.................................................................................||
int main(){
	clock_t t1,t2;
	t1=clock();
//	srand48(1095343394); //potential 3 works but orbits in plane
//	srand48(109534339); //potential 3 gets stuck with energy error
	srand48(109534331345);
	if(potential==4) { loadBeta(); }
	
	//				......... CONSTRUCTORS ..........
	Star star; // default contructor (IRS9)
	star.randomize();
	vector<double> q(3), vel(3), acc(3);
	q[0] = star.IC[0]; q[1] = star.IC[1]; q[2] = star.IC[2];
	vel[0] = star.IC[3]; vel[1] = star.IC[4]; vel[2] = star.IC[5];
	// USER INPUT (include which potential to tests)
	double time = 0, tmax=1e7*secyr, dt=10*secyr, n=1;
	double startE = calcE(q, vel); if(startE==1) { exit (EXIT_FAILURE); }
	star.storePosition(q);
	
	
	//----------------------------- While ---------------------------------|
	while (time < tmax) {
		vector<double> qtmp(q), vtmp(vel), atmp(acc),
		qlow(q), vlow(vel), qhi(q), vhi(vel);
		
		leap(qlow, vlow, acc, dt);
		leap(qhi, vhi, acc, 0.5*dt); leap(qhi, vhi, acc, 0.5*dt);	//	leap twice with dt/2

		for (int ii=0; ii<3; ii++) {	//	set q, v best
			q[ii] = 4./3.*qhi[ii] - 1./3.*qlow[ii];
			vel[ii] = 4./3.*vhi[ii] - 1./3*vlow[ii];
		}
		
		
		//		...... check error ......
		bool dtInd=0;
		if (dt>1*secyr) { dtInd = errCheck(q, qhi, vel, startE); }
		if (dtInd==1) { // reset variables, 1/2dt, repeat the step
			vec3Reset(q,qtmp); vec3Reset(vel,vtmp); vec3Reset(acc,atmp);
			dt = 0.5*dt;
		}
		else {// if (dtInd == 0) { //	dtInd == 0, error is fine, move to next step
			if (time>n*500.*secyr){ star.storePosition(q); n++;}
			time += dt;
			dt = doubleDt(dt);
		}
	}
	//---------------------------------------------------------------------|
	
	star.printOrbit(0);
	
	// Runtime
	t2=clock();
	cerr << "runtime: " << (t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}//................................................................................||



void leap(vector<double> &q, vector<double> &vel, vector<double> &acc, const double dt){
	for (int ii=0; ii<3; ii++){ q[ii] += 0.5*dt*vel[ii]; }
	calcAcc(q, acc);
	for (int ii=0; ii<3; ii++) { vel[ii] += dt*acc[ii]; q[ii] += 0.5*dt*vel[ii]; }
}


// check error in q, step back and 0.5dt if needed
bool errCheck(vector<double> &q, vector<double> &qhi, vector<double> &vel, const double startE){
	bool tmpReturn;
	double rhi = calcR(qhi,1);// sqrt(qhi[0]*qhi[0] + qhi[1]*qhi[1] + qhi[2]*qhi[2]);
	double rbest = calcR(q,1); //sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
	double deltaR = abs((rbest - rhi)/rbest);
	double energy = calcE(q, vel);
	double efrac = abs((startE-energy)/startE);
//	if (efrac > eFracMax) { eFracMax = efrac; }
	if ( (deltaR < 1e-5) && (efrac < 1e-5) ) { tmpReturn = 0; }
	else { tmpReturn = 1; }
	return tmpReturn;
}



//-----------------------------------ENERGY-----------------------------------//
double calcE(const vector<double> &q, const vector<double> &vel){
	double r = calcR(q, 1), R = calcR(q, 0);
	double vtot2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
	double Mstar = Msun, Mbh = 3.5e6*Msun;
	double potBH = -G*Mbh/r;
	
	//----------Spherically Symmetric
	if (potential==0) {
		double C = 1.396e4*Msun/pow(pc,3.), rc = 8*pc, rc2 = rc*rc, // Constants from Kenyon 2008
		potSS = -2*pi*G*C*rc2*(2*rc*atan(r/rc)/r+log(1+r*r/rc2)) + potBH;
		return (0.5*vtot2 + potSS)*Mstar;
	}
	
	//----------3 Component from Kenyon 2008
	else if (potential==3) {
		double Mb = 3.76e9*Msun, Md = 4.0e10*Msun, Mh = 1.0e12*Msun, ab = 100*pc, rh = 2.0e4*pc, ad = 2.0e3*pc, bd = 300*pc,	// Constants from Kenyon 2008
		bulge = -G*Mb/(r+ab),
		disk = -G*Md/sqrt(R*R+pow(ad+sqrt(q[2]*q[2]+bd*bd),2.)),
		halo = -G*Mh*log(1+r/rh)/r,
		pot3comp = bulge + disk + halo + potBH;
		return (0.5*vtot2 + pot3comp)*Mstar;
	}
	
	//----------4 Component from Kenyon 2008 plus NSC
	else if (potential==4) {
		double Mb = 3.76e9*Msun, Md = 4.0e10*Msun, Mh = 1.0e12*Msun, ab = 100*pc, rh = 2.0e4*pc, ad = 2.0e3*pc, bd = 300*pc,	// Constants from Kenyon 2008
		bulge = -G*Mb/(r+ab),
		disk = -G*Md/sqrt(R*R+pow(ad+sqrt(q[2]*q[2]+bd*bd),2.)),
		halo = -G*Mh*log(1+r/rh)/r,
		NSC = calcNSCpot(q),
		pot4comp = NSC + bulge + disk + halo + potBH;
		return (0.5*vtot2 + pot4comp)*Mstar;
	}
	
	else {	cerr << "Which potential do you want to test?\n0 = Spherically Symmetric\n3 = 3 Component\n4 = 4 Component\n";		return 1.;	}
	
}



//---------------------------------ACCELERATION---------------------------------//
void calcAcc(const vector<double> &q, vector<double> &acc){
	double r = calcR(q, 1);
	const double Mbh = 3.5e6*Msun,	accBH = -G*Mbh/r/r;
	
	//----------Spherically Symmetric
	if (potential==0) {
		double C = 1.396e4*Msun/pow(pc, 3.), rc = 8*pc, rc2 = rc*rc, // Constants from Kenyon 2008
		accSph = -4*pi*G*C*rc2*(r-rc*atan(r/rc))/r/r,
		a = accSph + accBH; //in r direction, spherical coords
		acc[0] = a*q[0]/r;
		acc[1] = a*q[1]/r;
		acc[2] = a*q[2]/r;
	}
	
	//----------3 Component
	else if (potential==3) {
		double Mb = 3.76e9*Msun, Md = 4.0e10*Msun, Mh = 1.0e12*Msun, ab = 100*pc, rh = 2.0e4*pc, ad = 2.0e3*pc, bd = 300*pc,	// Constants from Kenyon 2008
		
		bulge = -G*Mb/pow(r+ab,2.), // r direction
		halo = -G*Mh/r*(log(1+r/rh)/r - 1/(rh+r)), // r direction
		sqrtzb = sqrt(q[2]*q[2]+bd*bd), denom = pow(r*r+bd*bd+ad*ad+2*ad*sqrtzb,1.5),
		disk_xy = -G*Md/denom, // x and y directions
		disk_z = -G*Md*(ad+sqrtzb)/denom/sqrtzb, // z direction
		
		a_r = bulge+halo+accBH;
		
		acc[0] = q[0]*(a_r/r + disk_xy);
		acc[1] = q[1]*(a_r/r + disk_xy);
		acc[2] = q[2]*(a_r/r + disk_z);
	}
	
	//----------4 Component
	else if (potential==4){
		double Mb = 3.76e9*Msun, Md = 4.0e10*Msun, Mh = 1.0e12*Msun, ab = 100*pc, rh = 2.0e4*pc, ad = 2.0e3*pc, bd = 300*pc,	// Constants from Kenyon 2008
		bulge = -G*Mb/pow(r+ab,2.), // r direction
		halo = -G*Mh/r*(log(1+r/rh)/r - 1/(rh+r)), // r direction
		sqrtzb = sqrt(q[2]*q[2]+bd*bd), denom = pow(r*r+bd*bd+ad*ad+2*ad*sqrtzb,1.5),
		disk_xy = -G*Md/denom, // x and y directions
		disk_z = -G*Md*(ad+sqrtzb)/denom/sqrtzb, // z direction
		a_r = bulge+halo+accBH;
		vector<double> nscAcc(3);	calcNSCacc(q,nscAcc);
		
		acc[0] = q[0]*(a_r/r + disk_xy)+nscAcc[0];
		acc[1] = q[1]*(a_r/r + disk_xy)+nscAcc[1];
		acc[2] = q[2]*(a_r/r + disk_z)+nscAcc[2];
	}
}


// randomly change dt to 2*dt
inline double doubleDt(double dt){
	double ind = drand48();
	if(ind > 0.95) { dt = 2.*dt; }
	return dt;
}


void vec3Reset(vector<double> &reset, const vector<double> &value){
	reset.clear();		reset.resize(3);
	for (int ii=0; ii<3; ii++){ reset[ii] = value[ii]; }
}


// Calculate R - ind ==1 --> r (spherical coords), ind ==0 --> R (cylindrical coords)
double calcR(const vector<double> &q, const bool ind){
	if (ind == 1) {	return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);	}
	else { return sqrt(q[0]*q[0] + q[1]*q[1]); }
}


// load beta (const mass) vector for NSC potential
void loadBeta(){
	bool dmInd = 0;
	ifstream bfile("/Users/troyraen/dropbox/research/hvs/code/nscPotential/betaShells.dat");
	if (!bfile.is_open()) {	cerr << "Beta vector file could not be opened.\n";	 }
	double phold;
		if (dmInd==0){mShell = phold;	dmInd=1;}
		else {betaShells.push_back(phold);}
	}
	bfile.close();
	bsize = betaShells.size();
	if (bsize==0){cerr << "No data read from file to beta vector" << endl; }
}
