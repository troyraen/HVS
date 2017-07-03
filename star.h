//
//  star.h
//  ORBITAL PROPERTIES OF HIGH VELOCITY STARS NEAR SGR A*
//
//  Created by Troy Joseph Raen on 8/3/2014
//
//  STAR CLASS
//  in observed data, y is short axis. when
//
//

#ifndef _star_h
#define _star_h




class Star {
	
private:
	vector<double> Data;		//	[x, dx, y, dy, z, dz, Vx, dvx, Vy, dvy, Vz, dvz]
	vector<double> orbit;

public:
	string ID;
	double rmax,	rmin,	tinpc;
	int passes;
	vector<double> IC;	//	[x, y, z, vx, vy, vz]
	

	// constructor - enter data as in observed coord system >x' ^y' .z'
	Star(string IDtmp="IRS9", double x=6.64659e17, double dx=0., double y=-7.68314e17, double dy=0., double z=0., double dz=1.543e18, double Vx=1.09173e7, double dvx=8.27070e5, double Vy=1.14619e7, double dvy=9.45705e5, double Vz=-3.403e7, double dvz=59.5e5) {
		
		ID = IDtmp;
//		cout << "Constructor called for object: " << ID << ":\t Data" << endl;
		// star constructed so that z is the short axis or (x=x', y=-z', z=y') or phi=0, theta=pi/2
		Data.push_back(x);
		Data.push_back(dx);
		Data.push_back(-z);
		Data.push_back(-dz);
		Data.push_back(y);
		Data.push_back(dy);
		
		Data.push_back(Vx);
		Data.push_back(dvx);
		Data.push_back(-Vz);
		Data.push_back(-dvz);
		Data.push_back(Vy);
		Data.push_back(dvy);
	}

	
	void randomize(void);
	void printProperties(const bool err);
	void storePosition(const vector<double>& q);
	void storeProperties(const int numPasses, const double rmaxtmp, const double rmintmp, const double tinpctmp);
	void printOrbit(const bool errInd);

};


void Star::printOrbit(const bool errInd){	printObj(orbit,orbit.size(),3,errInd);	}


void Star::storePosition(const vector<double>& q){
	for (int ii=0; ii<3; ii++) { orbit.push_back(q[ii]/pc);		}
}


void Star::storeProperties(const int numPasses, const double rmaxtmp, const double rmintmp, const double tinpctmp){
	passes = numPasses;
	rmax = rmaxtmp;
	rmin = rmintmp;
	tinpc = tinpctmp;
}



// print properties
/*void Star::printProperties(const bool err){
	if (err==0){
		cout.precision(9);
		cout << rmax/pc << "\t" << rmin/pc << "\t" << passes << "\t" << tinpc << endl;
	}
	else {
		cerr.precision(9);
		cerr << rmax/pc << "\t" << rmin/pc << "\t" << passes << "\t" << tinpc << endl;
	}

}
*/
// randomize
void Star::randomize(void){
	IC.clear();
	
	IC.push_back(Data[0] + Data[1]*gaussrand());
	IC.push_back(Data[2] + Data[3]*gaussrand());
	IC.push_back(Data[4] + Data[5]*gaussrand());
	
	IC.push_back(Data[6] + Data[7]*gaussrand());
	IC.push_back(Data[8] + Data[9]*gaussrand());
	IC.push_back(Data[10] + Data[11]*gaussrand());
	
}




#endif