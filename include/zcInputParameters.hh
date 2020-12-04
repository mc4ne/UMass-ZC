#ifndef UMASS_ZERO_CROSSING_INPUT_PARAMETERS_H
#define UMASS_ZERO_CROSSING_INPUT_PARAMETERS_H

// a data struct to hold input parameters

namespace UMass {  

   struct zcInputParameters {
      double sampleFrequency;             // sampling frequency of the ADC (Hz) 
      double expectedFrequency;           // anticipated mixdown frequency (Hz) 
      double minTime;                     // starting time of zero crossing (sec)
      double maxTime;                     // end time of zero crossing (sec) 
      double timeThreshold;               // length of time at end of FID to determine noise level (sec)  
      int device;                         // set the device type (Utility::kPlungingProbe, Utility::kTrolley, Utility::kFixedProbe) 
      int offsetOrder;                    // what order offset correction to apply? (0 up to 5)  
      int verbosity;                      // how much info to print to screen (higher number = more text) 
      bool debug;                         // turn on debug mode
      bool useTimeRange;                  // use the time range 
      bool useMidpoint;                   // use the midpoint method to find the crossing  
      bool useLinearInterp;               // use a linear interpolation to find the crossing 
      bool useLeastSquares;               // use a least squares fit to find the crossing 
      bool useIntegerCycles;              // use an integer number of cycles when computing the frequency 
      bool useT2Time;                     // use T2 time in getting the frequency 
      bool useBaselineCorrection;         // use baseline correction on the FID  

      zcInputParameters(double sf=0,double ef=0,double t_min=0,double t_max=0,double t_thr=0,
	                int dev=0,int off_ord=0,int v=0,
	                bool db=false,bool tr=false,bool mid=false,bool lin=false,bool lsq=false,bool ic=false,bool t2=false,bool bl=false):
	 sampleFrequency(sf),expectedFrequency(ef),minTime(t_min),maxTime(t_max),
	 timeThreshold(t_thr),device(dev),offsetOrder(off_ord),verbosity(v),debug(db),useTimeRange(tr),useMidpoint(mid),useLinearInterp(lin),
	 useLeastSquares(lsq),useIntegerCycles(ic),useT2Time(t2),useBaselineCorrection(bl) {}

   } ;

} //::UMass

#endif  
