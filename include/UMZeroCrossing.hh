#ifndef UMASS_ZERO_CROSSING_H
#define UMASS_ZERO_CROSSING_H

// a class that does zero crossing analysis to extract the frequency of a signal 

#include <cstdlib> 
#include <iostream>

#include "zcInputParameters.hh"
#include "NMRPulse.hh"
#include "UMUtility.hh"
#include "UMBaselineCorrection.hh"

#define MAX_SIZE 1E+6 

namespace UMass { 
   class ZCTest{
      private:
	 int a;
      public:
	 ZCTest();
   };

   class ZeroCrossing {

      private:
	 zcInputParameters fParameters;                        // a container for all input parameters 

	 int fNPTS;                                              // number of points to use in least squares fitting 
	 int fNPTSUseable;                                       // number of points to use in least squares fitting (useable; we may have less -- see StoreData method) 
	 int fStep;                                              // number of points to skip in counting zero crossings
	 int *fZC;                                               // number of zero crossings 
	 int fNZC;                                               // number of zero crossings (single pulse) 

	 double *fX,*fY,*fEY;                                    // "analysis arrays": store data for fitting here  
	 double *fFREQ;                                          // final frequency results 
	 double *fFREQ_ph;                                       // final frequency results (phase fit of t_zc)  
	 double *fNC;                                            // number of cycles

	 int *fNCrossing,*fCrossingIndex;

	 double *fTcross,*fVcross;
	 double *fFreqAtCrossing,*fNumCycles;

	 NMRPulse *fPulse; 

	 BaselineCorrection *fBaseline;  

	 void Init(); 
	 void Reset(); 
	 void ClearAnaArrays(); 
	 void ClearVectors(); 
	 void InitAnaArrays();
	 void PrintVectorData(int Type,int PulseNumber); 
	 void CountZeroCrossings(int method,NMRPulse *aPulse);   // uses the UMMath::CountZeroCrossings function 
	 void SetNumPointsForFits(int n); 
	 void SetStepSize(int i)                         {fStep = i;} 
	 void UpdateParameters();                                // update fNPTS and fStep 

	 int CalculateFrequencies(int &TrueNumCrossings,double &FreqFullRange);
	 int Calculate(NMRPulse *aPulse);                        // runs the calculation based on choices 

	 int GetNumAnaBins()                       const {return fNZC;}  
	 int GetCrossingNumber(int i)              const {return fNCrossing[i];}  
	 int GetCrossingIndex(int i)               const {return fCrossingIndex[i];}  
	 int GetNumZeroCrossingsMidpoint()         const {return fZC[0];} 
	 int GetNumZeroCrossingsLinearInterp()     const {return fZC[1];} 
	 int GetNumZeroCrossingsLeastSquares()     const {return fZC[2];} 

	 double GetFrequencyFromPhaseFit();                        // compute the frequency fitting the zero crossings 
	 double GetFrequencyFromPhaseFit_nllsq(const int fitFunc); // compute the frequency fitting the zero crossings (NL-LSQ) 

	 double GetNumCyclesMidpoint()             const {return fNC[0];} 
	 double GetNumCyclesLinearInterp()         const {return fNC[1];} 
	 double GetNumCyclesLeastSquares()         const {return fNC[2];} 
	 double GetTimeAtCrossing(int i)           const {return fTcross[i];} 
	 double GetVoltageAtCrossing(int i)        const {return fVcross[i];} 
	 double GetFrequencyAtCrossing(int i)      const {return fFreqAtCrossing[i];} 
	 double GetNumberOfCylces(int i)           const {return fNumCycles[i];}

      public:
	 ZeroCrossing( zcInputParameters par = UMass::zcInputParameters() );
	 ~ZeroCrossing();

	 void UseAllMethods(){
	    UseMidpoint(); 
	    UseLinearInterpolation();
	    UseLeastSquares();
	 }

	 void UseTimeRange(bool val=true){
	    if(val) std::cout << "[ZeroCrossing]: Will use time range in counting zero crossings." << std::endl;
	    fParameters.useTimeRange = val;
	 }  
	 void UseMidpoint(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use midpoint method." << std::endl;
	    fParameters.useMidpoint = v;
	 }
	 void UseLinearInterpolation(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use linear interpolation method." << std::endl;
	    fParameters.useLinearInterp = v;
	 }
	 void UseLeastSquares(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use least squares method." << std::endl;
	    fParameters.useLeastSquares = v;
	 }
	 void UseIntegerCycles(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use integer number of cycles ONLY." << std::endl;
	    fParameters.useIntegerCycles = v;
	 }
	 void UseT2Time(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use T2 time in evaluating the frequency" << std::endl;
	    fParameters.useT2Time = v;
	 }
	 void UseBaselineCorrection(bool v=true){
	    if(v) std::cout << "[ZeroCrossing]: Will use the baseline correction" << std::endl;
	    fParameters.useBaselineCorrection = v;
	 }

	 void SetParameters(zcInputParameters par); 

	 void SetVerbosity(int v)                       {fParameters.verbosity         = v; } 
	 void SetSampleFrequency(double f)              {fParameters.sampleFrequency   = f; } 
	 void SetExpectedFrequency(double f)            {fParameters.expectedFrequency = f; } 
	 void SetMinTime(double t)                      {fParameters.minTime           = t; } 
	 void SetMaxTime(double t)                      {fParameters.maxTime           = t; }
	 void SetTimeRange(double tmin,double tmax){
	    fParameters.minTime = tmin;
	    fParameters.maxTime = tmax;
	 }     

	 // -- main function a user calls --  
	 // ampl is a vector of voltage readings;
	 // result is a vector of frequency results
	 // - order is: midpoint, linear interp, least squares, 
	 //             midpoint (phase), linear interp (phase), least squares (phase)
	 int Analyze(std::vector<double> ampl,std::vector<double> &result);        

	 // obtain frequency results individually 
	 // use the enumerated type anaType in the Utility namespace as the input type  
	 double GetFrequency(int type) const; 

   };

} //::UMass

#endif 
