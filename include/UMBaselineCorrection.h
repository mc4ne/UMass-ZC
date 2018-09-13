#ifndef UMASS_BASELINE_CORRECTION_H
#define UMASS_BASELINE_CORRECTION_H

// a class to handle the baseline corrections 

#include "UMMath.h"

class UMBaselineCorrection{

   private:
      double fTMin,fTMax;       // window in which to determine baseline  
      
      int CheckOffset(double offset_old,double offset_new,double t_diff_old,double t_diff_new,double slope); 

      double GetOffset(std::vector<double> v);
      double GetOffsetZC(double input_offset,std::vector<double> &v); 
      double ApplyOffset(double offset,std::vector<double> &v);
      double ApplyOffsetLinear(vector<double> par,std::vector<double> &v);
      double GetTDiff(vector<double> tCross,double &delta_t_even_nc,double &delta_t_odd_nc);

   public:
      UMBaselineCorrection(); 
      ~UMBaselineCorrection();

      void Init();
      void ApplyBaselineCorrection(std::vector<double> &v); 
  
      void SetTimeRange(double tMin,double tMax)             { fTMin = tMin; fTMax = tMax; } 

}; 

#endif 
