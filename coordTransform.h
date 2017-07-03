//
//	coordTransform.h
//
//	written by Troy Joseph Raen on 12/10/2014
//

#ifndef _coordTransform_
#define _coordTransform_

//#include <iostream>
#include <vector>
#include <complex>
//#include <cmath>
//#include <tgmath.h>
//#include <complex>
//#include "/Users/troyraen/dropbox/research/headers/ctemplates.h"
//#include "/Users/troyraen/dropbox/research/headers/constants.h"
//using namespace std;


// cartesian to oblate spheroidal
// short axis is z
void transCtoObl(const vector<double> &qCart, vector<double> &qObl, const double ae){
	double r = sqrt(qCart[0]*qCart[0]+qCart[1]*qCart[1]);
	complex<double> arg(r/ae, qCart[2]/ae);
	qObl[0] = real(acosh(arg)); // u
	qObl[1] = imag(acosh(arg)); // v
	qObl[2] = atan(qCart[1]/qCart[0]); // phi
}

void transObltoC(vector<double> &qCart, const vector<double> &qObl, const double ae){
	qCart[0] = ae*cosh(qObl[0])*cos(qObl[1])*cos(qObl[2]);
	qCart[1] = ae*cosh(qObl[0])*cos(qObl[1])*sin(qObl[2]);
	qCart[2] = ae*sinh(qObl[0])*sin(qObl[1]);
}

void transOblUnit(vector<double> &qCart, const vector<double> &qObl){
	double scale = 1./sqrt(sinh(qObl[0])*sinh(qObl[0])+sin(qObl[1])*sin(qObl[1]));
	qCart[0] = scale*sinh(qObl[0])*cos(qObl[1])*cos(qObl[2]);
	qCart[1] = scale*sinh(qObl[0])*cos(qObl[1])*sin(qObl[2]);
	qCart[2] = scale*cosh(qObl[0])*sin(qObl[1]);
}


#endif