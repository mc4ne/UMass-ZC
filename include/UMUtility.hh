#ifndef UMASS_UTIL_HH
#define UMASS_UTIL_HH

// miscellaneous functions needed by various classes 

#include <cstdlib> 
#include <iostream>
#include <fstream> 
#include <cmath>
#include <cfloat>
#include <vector> 
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

// #if USEGPUFID == 1
// #define DSPINC "dsp_gpu.h"
// #else
// #define DSPINC "dsp_cpu.h"
// #endif
// #include DSPINC

#include "NMRPulse.hh" 

namespace UMass { 
   namespace Utility{

      // data struct for nonlinear least squares fitting
      typedef struct data {
	 size_t n;
	 double *x;
	 double *y;
      } data_t; 

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

      double GetMean(std::vector<double> v);
      double GetVariance(std::vector<double> v);
      double GetStandardDeviation(std::vector<double> v);
      double GetCovariance(std::vector<double> x,std::vector<double> y); 

      double ConvertToVoltage3316(double adc_reading); 

      double LinearInterpolationForX(double y,double x0,double y0,double x1,double y1); 
      double LinearInterpolationForY(double x,double x0,double y0,double x1,double y1);
      // double LinearInterpolation(double x,double x0,double y0,double x1,double y1);

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

      double GetT2Time_old(NMRPulse *aPulse); 
      double GetT2Time(int startIndex,NMRPulse *aPulse);

      int FindLocalMaxima(int startIndex,NMRPulse *aPulse,std::vector<double> &T,std::vector<double> &V); 
      int RebinData(int stepSize,std::vector<double> x,std::vector<double> y,
	    std::vector<double> &X,std::vector<double> &Y);

      // non-linear least squares fit 
      int NonLinearLeastSquaresFitting(std::vector<double> x,std::vector<double> y,std::vector<double> dy,
	    int (*F)(const gsl_vector *x,void *data,gsl_vector *f),int (*DF)(const gsl_vector *x,void *data,gsl_matrix *J),
	    std::vector<double> &par,std::vector<double> &parErr,const int NPAR,const int verbosity); 
      void callbackFunction(const size_t iter, void *params,const gsl_multifit_nlinear_workspace *w);      

      // phase fit functions 
      int poly3(const gsl_vector * x, void *data,gsl_vector * f);
      int poly3_df(const gsl_vector *x,void *data,gsl_matrix *J);
      int poly5(const gsl_vector * x, void *data,gsl_vector * f);
      int poly5_df(const gsl_vector *x,void *data,gsl_matrix *J);
      int poly7(const gsl_vector * x, void *data,gsl_vector * f);
      int poly7_df(const gsl_vector *x,void *data,gsl_matrix *J);

      int AdjustTimeWindow(NMRPulse *aPulse,double &tStart,double &tStop);  

   } //::Utility 
} //::UMass

#endif 
