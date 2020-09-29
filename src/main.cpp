#include <iostream>
#include "Burden.h"
#include "Person.h"
#include "Parameters.h"
#include "Population.h"
#include "Simulation.h"

#include "scenarios.h"

using namespace std;

int main(int argc, const char * argv[]) {

    
    if(argc==1){
      cout<<"No scenario selected";
      return 0;
    }
    if(atoi(argv[1])==0){
      scenario_square();
    }
    if(atoi(argv[1])==1){
      scenario_swz();
    }
    return 0;

    //scenario_swz();
    //return 0;
}
