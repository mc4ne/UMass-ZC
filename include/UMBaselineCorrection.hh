#ifndef UMASS_BASELINE_CORRECTION_H
#define UMASS_BASELINE_CORRECTION_H

// a class to handle the baseline corrections 

//#include "dsp_cpu.h"
//#include "NMRPulse.h"
#include "zcInputParameters.hh"
#include "UMUtility.hh"

namespace UMass { 

   class BaselineCorrection {

      private:
	zcInputParameters fParameters; // for input parameters
 
	 bool fOffsetFail; 

	 int fSIZE; 
	 int *fNCrossing,*fCrossingIndex;

	 double *fX,*fY,*fEY,*fTcross,*fVcross; 

	 void Init();
         void ClearDataArrays(); 
         void ClearNZCArrays();

 	 void ApplyOffset(double offset,NMRPulse *aPulse);
	 void ApplyOffsetLinear(double *par,NMRPulse *aPulse);
	 void GetOffsetLinear(double input_offset,NMRPulse *aPulse,double *offset); 

	 int CheckOffset(double offset_old,double offset_new,double t_diff_old,double t_diff_new,double slope); 

	 double GetOffset(double t_thr,NMRPulse *aPulse); 
	 double GetOffsetZC(double input_offset,int stepSize,NMRPulse *aPulse); 
	 double GetTDiff(int nzc,double *tCross,double &delta_t_even_nc,double &delta_t_odd_nc);

      public:
	 BaselineCorrection( zcInputParameters par = UMass::zcInputParameters()); 
	 ~BaselineCorrection();

	 void SetVerbosity(int v) {fParameters.verbosity = v;} 
         void SetParameters(zcInputParameters par); 
 
         // -- function the user calls -- 
	 int ApplyBaselineCorrection(NMRPulse *aPulse); 

   }; 

} //::UMass 

#endif 
