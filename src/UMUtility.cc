#include "UMUtility.h"
//______________________________________________________________________________
namespace UMass { 
   namespace Utility {
      //______________________________________________________________________________
      int PrintSignalToFile(const char *outpath,NMRPulse *aPulse){
	 const int N = aPulse->GetNumPoints();
	 char outStr[500];
	 std::ofstream outfile;
         outfile.open(outpath);
	 if( outfile.fail() ){
	    std::cout << "[Utility::PrintSignalToFile]: Cannot open the file: " << outpath << std::endl;
	    return 1;
	 }else{
	    for(int i=0;i<N;i++){
	       sprintf(outStr,"%.3lf",aPulse->GetVoltage(i));
	       outfile << outStr << std::endl;
	    }
	    outfile.close();
	    std::cout << "[Utility::PrintSignalToFile]: The data has been written to the file: " << outpath << std::endl;
	 }
	 return 0;
      }
      //______________________________________________________________________________
      int CutPlungingProbeData(std::vector<double> &t,std::vector<double> &v){
	 // a data-quality cut on PP data
         // - tMin: We know that we have to cut out the first 1 us of data to avoid 
         //   an odd "ground level synchronization" that we see between the DAQ and the ADC
         // - tMax: We need to avoid the last 0.5 ms of data because of the Radant MEMS 
         //   which causes a level change 
         int N = t.size();
         double tMin = 1E-6; 
         double tMax = t[N-1] - 0.5E-3;
	 std::vector<double> tt,vv;
         for(int i=0;i<N;i++){
	    if(t[i]>tMin && t[i]<tMax) {
               tt.push_back(t[i]); 
	       vv.push_back(v[i]); 
	    }
         }
         // clear original vectors and set with cut values  
	 t.clear();
         v.clear();
         N = tt.size();
         for(int i=0;i<N;i++){
	    t.push_back(tt[i]); 
	    v.push_back(vv[i]); 
         } 
         return 0;
      }
      //______________________________________________________________________________
      int ConvertToVoltage(int type,std::vector<double> &v){
	 const int N = v.size();
         for(int i=0;i<N;i++){
	    if(type==k3316) v[i] = ConvertToVoltage3316(v[i]); 
         } 
	 return 0;
      }
      //______________________________________________________________________________
      double ConvertToVoltage3316(double adc_reading){
	 // unsigned short int = 16-bit-wide unsigned integer 
	 // convert the 16-bit integer to a voltage
	 // ADC details 
	 // Manufacturer: Struck 
	 // Model:        SIS3316 125 MHz 16-bit (or 250 MHz 14-bit) ADC
	 // Input Span:   5 V (-2.5,..,2.5 V)  
	 // Impedance:    50 Ohms 
	 // Input Type:   ?? 
	 // use a linear fit to ADC counts vs voltage 
	 double adc       = (double)adc_reading;
	 // partial range
	 // double p0        = 32582.3;  
	 // double p1        = 12630.5;  
	 // full range
	 double p0        = 32558.5;
	 double p1        = 12629.5;
	 double v         = (adc - p0)/p1;
	 return v;
      } 
      //______________________________________________________________________________
      double GetT2Time(NMRPulse *aPulse){
	 // find the T2 time of the signal
	 // find max amplitude 
	 double vmax=-300,v=0;
	 const int N = aPulse->GetNumPoints();
	 for(int i=0;i<N;i++){
	    v = aPulse->GetVoltage(i);
	    if(vmax<v) vmax = v;
	 }
	 // find T2 time 
	 double t2_time=0;
	 double e_const = exp(1);
	 double v_lo   = (vmax/e_const)*(1. - 0.1);
	 double v_hi   = (vmax/e_const)*(1. + 0.1);
	 // std::cout << "VMAX = " << vmax << " VLO = " << v_lo << " VHI = " << v_hi << std::endl; 
	 for(int i=0;i<N;i++){
	    v = aPulse->GetVoltage(i);
	    if( fabs(v)>v_lo && fabs(v)<v_hi ){
	       t2_time = aPulse->GetTime(i);
	    }
	 }
	 // std::cout << "T2 TIME IS " << t2_time << std::endl; 
	 return t2_time;
      }
      //______________________________________________________________________________
      int StoreData(int verbosity,int i,int NPTS,NMRPulse *aPulse,double *X,double *Y,double *EY){

	 int N     = aPulse->GetNumPoints(); 
	 int start = i - NPTS/2;
	 int end   = i + NPTS/2;

	 double v_current = aPulse->GetVoltage(i); 

	 if(v_current!=0){
	    // do nothing
	 }else{
	    // voltage of the zero crossing is zero! 
	    start = i - 3; 
	    end   = i + 3; 
	    if(verbosity>=3){  
	       std::cout << "[UMUtility::StoreData]: Voltage at zero crossing is zero!" << std::endl;
	       std::cout << "                              start = " << start << std::endl;
	       std::cout << "                              end   = " << end   << std::endl;
	    }
	 }

	 // prevent unrealistic bounds: use 6 data points for the fit 
	 if(start < 0){
	    start = i; 
	    end   = i + 6; 
	    if(verbosity>=3){ 
	       std::cout << "[UMUtility::StoreData]: Invalid start point!  Setting to index = " << i << std::endl; 
	       std::cout << "                              start = " << start << std::endl;
	       std::cout << "                              end   = " << end   << std::endl;
	    }
	 }else if(start==0){
	    start = 0; 
	    end   = 7; 
	    if(verbosity>=3){
	       std::cout << "[UMUtility::StoreData]: Starting at index = " << start << std::endl; 
	       std::cout << "                              start = " << start << std::endl;
	       std::cout << "                              end   = " << end   << std::endl;
	    }
	 }

	 if(end > N){
	    start = N - NPTS; 
	    end   = N;
	 } 

	 int k=0;
	 for(int j=start;j<end;j++){  
	    X[k]  = aPulse->GetTime(j); 
	    Y[k]  = aPulse->GetVoltage(j);  
	    EY[k] = aPulse->GetVoltageErr(j);
	    k++; 
	 }

	 int NPTSUseable= NPTS; 

	 if(k!=NPTS){
	    NPTSUseable = k; 
	    if(verbosity>=3){
	       std::cout << "[UMUtility::StoreData]: WARNING!  Do not have the expected number of data points! " << std::endl;
	       std::cout << "                              k    = " << k    << std::endl;
	       std::cout << "                              NPTS = " << NPTS << std::endl;
	       std::cout << "                              Using k data points in fit..." << std::endl;
	    }
	 }

	 return NPTSUseable; 

      }
      //______________________________________________________________________________
      void ClearAnaArrays(int N,double *X,double *Y,double *EY){
	 for(int i=0;i<N;i++){
	    X[i]  = 0; 
	    Y[i]  = 0; 
	    EY[i] = 0; 
	 }
      }
      // //______________________________________________________________________________
      // double LinearInterpolationForY(double x,double x0,double y0,double x1,double y1){
      //    double b = (x-x0)/(x1-x0);
      //    double y = y0 + b*(y1-y0);
      //    return y;
      // }
      // //______________________________________________________________________________
      // double LinearInterpolationForX(double y,double x0,double y0,double x1,double y1){
      //    double b = (y-y0)/(y1-y0);
      //    double x = x0 + b*(x1-x0);
      //    return x;
      // }
      //______________________________________________________________________________
      double LinearInterpolation(double x,double x0,double y0,double x1,double y1){
	 double b = (x-x0)/(x1-x0);
	 double y = y0 + b*(y1-y0);
	 return y;
      }
      //______________________________________________________________________________
      int LeastSquaresFitting(int N,double *x,double *y,double &a,double &b,double &r){

	 // linear regression to find slope b and y-intercept a of 
	 // f(x) = a + bx 

	 int ret_val=0;
	 double num=0,rsq=0;

	 // FIXME: make the function call above to vectors... 
	 std::vector<double> X(N),Y(N);
	 for(int i=0;i<N;i++){
	    X[i] = x[i]; 
	    Y[i] = y[i]; 
	 } 

	 double mean_x = dsp::mean(X);
	 double mean_y = dsp::mean(Y);
	 double var_x  = dsp::variance(X);
	 double var_y  = dsp::variance(Y);
	 double cov_xy = dsp::covariance(X,Y);

	 double ss_xx = ( (double)N )*var_x;
	 double ss_yy = ( (double)N )*var_y;
	 double ss_xy = ( (double)N )*cov_xy;

	 double den = ss_xx*ss_yy;
	 if(den==0){
	    // singular matrix. can't solve the problem.
	    a       = 0;
	    b       = 0;
	    r       = 0;
	    ret_val = 1;
	 }else{
	    b   = cov_xy/var_x;
	    a   = mean_y - b*mean_x;
	    num = ss_xy*ss_xy;
	    rsq = num/den;
	    r   = sqrt(rsq);
	 }

	 return ret_val;
      }
      //______________________________________________________________________________
      int LeastSquaresFitting(std::vector<double> x,std::vector<double> y,double &a,double &b,double &r){

	 // linear regression to find slope b and y-intercept a of 
	 // f(x) = a + bx 

	 int ret_val=0;
	 double num=0,rsq=0;

	 int N         = x.size();
	 double mean_x = dsp::mean(x);
	 double mean_y = dsp::mean(y);
	 double var_x  = dsp::variance(x);
	 double var_y  = dsp::variance(y);
	 double cov_xy = dsp::covariance(x,y);

	 double ss_xx = ( (double)N )*var_x;
	 double ss_yy = ( (double)N )*var_y;
	 double ss_xy = ( (double)N )*cov_xy;

	 double den = ss_xx*ss_yy;
	 if(den==0){
	    // singular matrix. can't solve the problem.
	    a       = 0;
	    b       = 0;
	    r       = 0;
	    ret_val = 1;
	 }else{
	    b   = cov_xy/var_x;
	    a   = mean_y - b*mean_x;
	    num = ss_xy*ss_xy;
	    rsq = num/den;
	    r   = sqrt(rsq);
	 }

	 return ret_val;
      }
      //______________________________________________________________________________
      double GetTimeOfCrossing(int verbosity,int method,int NPTSUseable,double X[],double Y[],double EY[], 
	    double t_current,double v_current,double v_current_err,
	    double t_next   ,double v_next   ,double v_next_err){

	 const int SIZE = NPTSUseable;

	 int ret_val=0; 
	 double v0=0,t0=0,a=0,b=0,r=0;

	 if(method==kMidpoint){
	    // method 1: take midpoint between t_current and t_next 
	    t0 = (t_current + t_next)/2.;
	 }else if(method==kLinearInterpolation){
	    // method 2: get time at V = 0, linear interpolation  
	    // t0 = LinearInterpolationForX(v0,t_current,v_current,t_next,v_next);
	    t0 = LinearInterpolation(v0,v_current,t_current,v_next,t_next);
	    // std::cout << "linear interpolation: t_current = " << t_current << "\t" 
	    //           << "t_next = " << t_next << "\t" << "t0 = " << t0 << std::endl;
	 }else if(method==kLeastSquares){
	    // method 3: least squares fit to neighboring points
	    // to find fit parameters a and b in f(x) = a + bx 
	    ret_val = LeastSquaresFitting(SIZE,X,Y,a,b,r);
	    if(b!=0){
	       t0 = -a/b;
	    }else{
	       t0 = 0;
	    }
	    // make sure t0 is bound properly 
	    if(t0<X[0] || t0>X[SIZE-1]){
	       if(verbosity>=3){
		  std::cout << "[UMUtility::GetTimeOfCrossing]: ERROR!  t0 is out of bounds!  "; 
		  std::cout << "Using linear interpolation instead..." << std::endl;
		  std::cout << "                                t_min = " << X[0] << "\t" 
		     << "t_max = " << X[SIZE-1] << "\t" << "t0 = " << t0 << std::endl;
	       }
	       // t0 = LinearInterpolationForX(v0,t_current,v_current,t_next,v_next);
	       t0 = LinearInterpolation(v0,v_current,t_current,v_next,t_next);
	       if(verbosity>=3){
		  std::cout << "[UMUtility::GetTimeOfCrossing]: linear interpolation: t_current = " << t_current << "\t" 
		     << " t_next = " << t_next << "\t" << "t0 = " << t0 << std::endl;
	       }
	    }
	 }else{
	    // invalid method, set to -1
	    t0 = -1; 
	 }

	 if(verbosity>=3){  
	    if(t0<0){
	       std::cout << "[UMUtility::GetTimeOfCrossing]: BAD CROSSING TIME!" << std::endl;
	       std::cout << "                                      t0        = " << t0        << std::endl;
	       std::cout << "                                      method    = " << method    << std::endl;
	       std::cout << "                                      t_current = " << t_current << std::endl;
	       std::cout << "                                      t_next    = " << t_next    << std::endl;
	       if(method==3){
		  std::cout << "                                      offset = " << a << std::endl;
		  std::cout << "                                      slope  = " << b << std::endl;
		  for(int i=0;i<SIZE;i++){
		     std::cout << "                                      t = " << X[i] << "\t" << "v = " << Y[i] << std::endl;
		  }
	       }
	    }
	 }

	 // to get rid of warnings 
	 ret_val += 0;

	 return t0; 

      }
      //______________________________________________________________________________
      int CountZeroCrossings(int verbosity,int method,int NPTS,int step,
	    bool UseTimeRange,double tMin,double tMax,
	    NMRPulse *aPulse,
	    double *X,double *Y,double *EY, 
	    int *NCrossing,int *CrossingIndex,double *tCross,double *vCross){

	 if(verbosity>=3) std::cout << "[UMUtility]: Counting zero crossings..." << std::endl;

	 // NPTS = number of points for linear fit
	 // Step = how many points to skip ahead in counting zero crossings 

	 int NPTSUseable = 0;

	 const int N  = aPulse->GetNumPoints();

	 int cntr         = 0;
	 // int cntr_prev    = 0;

	 double v0           = 0;
	 double target       = 0;
	 double t0           = 0;
	 // double elapsed_time = 0;  

	 double v_prod=0;
	 double delta_v=0; 
	 double v_current=0,v_next=0;
	 double t_current=0,t_next=0;
	 double v_current_err=0,v_next_err=0;

	 int i=0;  // index for NMRPulse data 
	 do{
	    t_current     = aPulse->GetTime(i);
	    t_next        = aPulse->GetTime(i+1);
	    v_current     = aPulse->GetVoltage(i);
	    v_next        = aPulse->GetVoltage(i+1);
	    v_current_err = aPulse->GetVoltageErr(i);
	    v_next_err    = aPulse->GetVoltageErr(i+1);
	    v_prod        = v_current*v_next;
	    if(v_prod > target){
	       // positive number, no crossing
	       // increment i by 1 
	       i++; 
	    }else if( v_prod <= target ){
	       // negative number or ZERO, we had a crossing
	       if(UseTimeRange){
		  // use the fit range
		  if(t_current > tMin && t_next < tMax){ 
		     // count the crossing 
		     cntr++;
		     delta_v     = fabs(v_current-v_next); 
		     // fill X, Y, EY arrays for fit method 
		     NPTSUseable = StoreData(verbosity,i,NPTS,aPulse,X,Y,EY); 
		     // get time of crossing  
		     t0 = GetTimeOfCrossing(verbosity,method,NPTSUseable,X,Y,EY,t_current,v_current,v_current_err,t_next,v_next,v_next_err);
		     // fill arrays
		     NCrossing[cntr-1]     = cntr;
		     CrossingIndex[cntr-1] = i;
		     tCross[cntr-1]        = t0;  
		     vCross[cntr-1]        = v0;  
		  }
	       }else{
		  // don't use the fit range 
		  // count the crossing 
		  cntr++;
		  delta_v = fabs(v_current-v_next); 
		  NPTSUseable = StoreData(verbosity,i,NPTS,aPulse,X,Y,EY); 
		  // get time of crossing  
		  t0 = GetTimeOfCrossing(verbosity,method,NPTSUseable,X,Y,EY,t_current,v_current,v_current_err,t_next,v_next,v_next_err); 
		  // fill arrays
		  NCrossing[cntr-1]     = cntr;
		  CrossingIndex[cntr-1] = i;
		  tCross[cntr-1]        = t0;  
		  vCross[cntr-1]        = v0;  
	       }
	       // check t0
	       if(t0<0 || delta_v>0.100){   // we shouldn't see a 100 mV jump during the zero crossing 
		  if(verbosity>=4){
		     std::cout << "[UMUtility::CountZeroCrossings]: bad crossing for Zc = " 
			<< cntr << "!  Trying next crossing..." << std::endl;
		     std::cout << "                               t0      = " << t0      << std::endl;
		     std::cout << "                               delta_v = " << delta_v << std::endl;
		  }
		  cntr--;   // don't count the crossing, decrement cntr 
		  i += step;   // move to next bin 
		  // clear analysis arrays 
		  ClearAnaArrays(NPTSUseable,X,Y,EY);           // clears X, Y, EY  (sets to zero) 
		  if(verbosity>=4){ 
		     std::cout << "[UMUtility::CountZeroCrossings]: zero crossing counter reset to: " << cntr << std::endl;
		     std::cout << "[UMUtility::CountZeroCrossings]: moving to index:                " << i       << std::endl;
		  }
		  continue; 
	       }
	       // passed t0 check 
	       i += step;  // move to next bin  
	       // set up for next data point 
	       ClearAnaArrays(NPTSUseable,X,Y,EY);             // clears X, Y, EY  (sets to zero) 
	    }
	 }while( i<(N-1) ); 

	 if(verbosity>=3) std::cout << "[UMUtility::CountZeroCrossings]: Done." << std::endl;
	 ClearAnaArrays(NPTSUseable,X,Y,EY);                   // clears fX, fY, fEY  (sets to zero) 

	 return cntr;   // return number of zero crossings  
      }
   } //::Utility 
} //::UMass
