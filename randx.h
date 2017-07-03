#ifndef _randx_h
#define _randx_h

#include <math.h>

#ifndef DRANDX
#include <stdlib.h>
#define _DRANDX_ drand48
#else
#define _DRANDX_ DRANDX
#endif

/* gaussian dev gen. (Num. Rec.+minor change) */
inline double gaussrand()
{
    static int iset = 0;
    static double gset;
    double fac,r,v1,v2;
	
    if (iset == 0) {
		do {
			v1 = 2.*_DRANDX_()-1.;
			v2 = 2.*_DRANDX_()-1.;
			r = v1*v1+v2*v2;
		} while (r >= 1. || r == 0.);
		fac = sqrt(-2.*log(r)/r);
		gset = v1*fac;
        iset = 1;
		return(v2*fac);
    } else {
        iset = 0;
    	return(gset);
    }
}
#undef _DRANDX_

#endif