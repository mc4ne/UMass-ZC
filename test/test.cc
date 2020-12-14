// test program to utilize UMass-ZC 

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

#include "zcInputParameters.hh"
#include "UMZeroCrossing.hh"
#include "CSVManager.hh"

int main(){

   char inpath[200]; 
   sprintf(inpath,"./input/PNMR_He_0.15s_0001_20200830_192539.dat"); 

   CSVManager *myCSV = new CSVManager(0,"tsv");
   myCSV->ReadFile(inpath,true); // header in file   

   std::vector<double> v; 
   myCSV->GetColumn_byName<double>("X-chan(V)",v);  

   UMass::zcInputParameters par;
   par.sampleFrequency       = 200E+3; 
   par.expectedFrequency     = 330;
   par.minTime               = 82E-3;   // 82 ms  
   par.maxTime               = 110E-3;  // 110 ms  
   par.timeThreshold         = 1E-3;    // 1 ms for noise baseline determination 
   par.useTimeRange          = true; 
   par.useMidpoint           = true;  
   par.useLinearInterp       = true;  
   par.useLeastSquares       = true;  
   par.useT2Time             = false;  
   par.useBaselineCorrection = false;  

   UMass::ZeroCrossing *myZC = new UMass::ZeroCrossing(par); 

   // analyze data 
   std::vector<double> res; 
   myZC->Analyze(v,res); 

   // print to screen 
   std::vector<std::string> label; 
   label.push_back("midpoint"); 
   label.push_back("linear-interp"); 
   label.push_back("least-squares"); 
   label.push_back("midpoint (phase)"); 
   label.push_back("linear-interp (phase)"); 
   label.push_back("least-squares (phase)"); 

   std::cout << "Zero Crossing Analysis Result: " << std::endl;
   char msg[200]; 
   const int N = res.size();
   for(int i=0;i<N;i++){
      sprintf(msg,"%20s: %.3lf Hz",label[i].c_str(),res[i]);
      std::cout << msg << std::endl;
   }

   delete myZC;  

   return 0;
}
