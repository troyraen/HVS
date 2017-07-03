//	run.h
//
//  ORBITAL PROPERTIES OF HIGH VELOCITY STARS NEAR SGR A*
//
//	Created by Troy Raen on 1/28/15.
//
//	header file for main.cpp
//	run Class
//




#ifndef _run_h
#define _run_h




class Run {
private:
	
public:
	std::vector<double>q, vel, acc;
	double r, R, startE, currentE;
	int potential;
	double calcR(const bool ind);
	void Run::setR(void)
	void getStartE();

	
	//				.............. CONSTRUCTOR .............
	Run(Star &star){
	}
	
};








void Run::setR(void){
	r = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
	R = sqrt(q[0]*q[0] + q[2]*q[2]);
}



void Run::getStartE(void){	startE = calcE();	}












#endif
