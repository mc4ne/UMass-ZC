#ifndef UMASS_UTIL_HH
#define UMASS_UTIL_HH

// miscellaneous functions needed by various classes 

#include <cstdlib> 
#include <iostream>
#include <fstream> 
#include <string.h>

// #if USEGPUFID == 1
// #define DSPINC "dsp_gpu.h"
// #else
// #define DSPINC "dsp_cpu.h"
// #endif
// #include DSPINC

#include "NMRPulse.hh" 

namespace UMass { 
   namespace Utility{ 

      enum anaType{
	 kMidpoint                 = 1,
	 kLinearInterpolation      = 2,
	 kLeastSquares             = 3,
	 kMidpointPhase            = 4,
	 kLinearInterpolationPhase = 5,
	 kLeastSquaresPhase        = 6
      };
     
      enum adcType{
	 k3316 = 1,
	 k3302 = 2
      };

      enum deviceType{
	 kPlungingProbe = 1,
	 kTrolley       = 2,
	 kFixedProbe    = 3
      }; 
 
      void ClearAnaArrays(int N,double *X,double *Y,double *EY); 

      double ConvertToVoltage3316(double adc_reading); 
      double GetT2Time(NMRPulse *aPulse); 
      double LinearInterpolation(double x,double x0,double y0,double x1,double y1);
      double GetTimeOfCrossing(int verbosity,int method,int NPTSUseable,double X[],double Y[],double EY[],
	    double t_current,double v_current,double v_current_err,
	    double t_next   ,double v_next   ,double v_next_err);

      int ConvertToVoltage(int type,std::vector<double> &v);
      int PrintSignalToFile(const char *outpath,NMRPulse *aPulse);  
      int CutPlungingProbeData(std::vector<double> &t,std::vector<double> &v);
      int AreEquivStrings(const char *s1,const char *s2);
      int StoreData(int verbosity,int i,int NPTS,NMRPulse *aPulse,double *X,double *Y,double *EY);

      int CountZeroCrossings(int verbosity,int method,int NPTS,int step,
	    bool UseTimeRange,double tMin,double tMax,
	    NMRPulse *aPulse,
	    double *X,double *Y,double *EY,
	    int *NCrossing,int *CrossingIndex,double *tCross,double *vCross); 

      int LeastSquaresFitting(int N,double *x,double *y,double &a,double &b,double &r); 
      int LeastSquaresFitting(std::vector<double> x,std::vector<double> y,double &a,double &b,double &r); 

   } //::Utility 
} //::UMass

#endif 
