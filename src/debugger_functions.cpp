#include "debugger_functions.h"
#include "helper_functions.h"

using std::cout;
using std::endl;

void LOGLandmarkObsIDs(const std::vector<LandmarkObs> myvector){
  	cout<< "Printing LandMarksObs Ids" <<endl;
	for(unsigned int i=0; i<myvector.size();++i){
    	cout<< i <<": " <<myvector[i].id << "\t";
    }
  	cout<< endl;
}

