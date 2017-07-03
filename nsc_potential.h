//
//	nsc_potential.h (funcitons taken from nsc_potential.cc)
//
//	written by Troy Joseph Raen on 1/20/2015
//
// nuclear star cluster potential from... Anil, Schodel, Trujillo
//-------- BT oblate spheroidal shells ---------//
//

#ifndef _nscpotential_h
#define _nscpotential_h
//#include <cmath>
#include "coordTransform.h"




// BT spheroidal shell potential
double phiBT(const double ratio, const double a, const double e, const double ref, const double deltaM){
	// ratio = ref/u0
	if (ratio < 1.)	{	return -G*deltaM/a/e*asin(e);	}
	else			{	return -G*deltaM/a/e*asin(1/cosh(ref));	}
}






// calculate NSC potential
double calcNSCpot(const vector<double> &coords){
	// CONSTANTS
	const double halfR = 4.2*pc, a = halfR, q = 0.71, ecc = sqrt(1-q*q), u0 = acosh(1/ecc); double phi = 0.; vector<double> coordsObl(3);
	extern vector<double> betaShells; extern int bsize; extern double mShell;
	for (int n = 0; n<bsize; n++) {
		double ai = betaShells[n]*a;
		transCtoObl(coords, coordsObl, ai*ecc);
		double phiTmp = phiBT(coordsObl[0]/u0, ai, ecc, coordsObl[0], mShell/betaShells[n]);
		phi += phiTmp;
	}
	return phi;
}


// calculate NCS acceleration
// see potGradNSC.mw for derivation from phi
void calcNSCacc(const vector<double> &coords, vector<double> &acc){
	vector<double> coordsObl(3), Unit(3); extern vector<double> betaShells; extern int bsize; double ecc = sqrt(1-0.71*0.71), a = 4.2*pc;
	// WHAT is a here?????????
	transCtoObl(coords, coordsObl, a*ecc);
	double coeff = -G*Msun*1./cosh(coordsObl[0])*tanh(coordsObl[0])/sqrt(cosh(coordsObl[0])*cosh(coordsObl[0])-1+cos(coordsObl[1])*cos(coordsObl[1]))/a/ecc/sqrt(1-1./cosh(coordsObl[0])*1./cosh(coordsObl[0]));
	transOblUnit(Unit, coordsObl);
	for (int i = 0; i<3; i++){ acc[i] = coeff*Unit[i]; }
}



#endif
