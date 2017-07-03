//  ctemplates.h
//
//  Created by Troy Raen on 5/28/14.
//
// << std::setprecision(std::numeric_limits<long double>::digits10 + 1) 


#ifndef _ctemplates_h
#define _ctemplates_h
//#include <iostream>


//using namespace std;


template<class TYPE>
inline void printObj(TYPE &obj, int row, int col, bool err) {
	
	if (err == 0) {
		for ( int i = 0; i < row; i++ ) {
			if (i!=0 && i%col == 0) {	cout << "\n";	}
//			for ( int j=0; j < col; j++) {
			std::cout << obj[i] << "\t";
//			}
//			cout << endl;
		}
//		cout << endl;
	}
	else {
		for ( int i = 0; i < row; i++ ) {
//			for (int j=0; j < col; j++){
//				if (j!=0 && j%col == 0) {	cerr << "\n";	}
			std::cerr << obj[i] << "\t";
//			}
			std::cerr << endl;
		}
		std::cerr << endl;
	}
}



#endif
