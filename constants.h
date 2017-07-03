//  Constants.h
//  Troy Raen
//  6/14/2014
//
//  Useful Astronomical Constants and Conversions
//


#ifndef _constants_h
#define _constants_h


#include <cmath>
#include <math.h>
#include <iostream>
using namespace std;


//---------------Constants--------------//
double pi = 4.*atan(1.),
	Msun = 1.99e33,  // grams
	AU = 1.496e13,  // cm
	pc = 3.086e18,  // cm
	G = 6.67259e-8; // cm^3/g*sec^2



//---------------Conversions-------------//
double radeg =	pi/180.,  // radians per degree
	secyr = pow(75.,4.);  // seconds per year


double asdeg = 3600.,  // arcseconds per degree
	masdeg = 1e3*3600.;  // milli.arcseconds per degree







#endif