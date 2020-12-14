#include "UMBaselineCorrection.hh"

namespace UMass {
   //______________________________________________________________________________
   BaselineCorrection::BaselineCorrection(zcInputParameters par){
      Init();
      SetParameters(par);
   }
   //______________________________________________________________________________
   BaselineCorrection::~BaselineCorrection(){
      delete[] fX;
      delete[] fY;
      delete[] fEY;
      delete[] fNCrossing;
      delete[] fCrossingIndex;
      delete[] fTcross;
      delete[] fVcross;
   }
   //______________________________________________________________________________
   void BaselineCorrection::Init(){
      fOffsetFail    = false;
      // fOffsetOrder   = 0; 
      // fParameters.verbosity     = 0;
      fSIZE          = 10E+6;
      const int N    = fSIZE;
      fX             = new double[N];
      fY             = new double[N];
      fEY            = new double[N];
      fNCrossing     = new int[N];
      fCrossingIndex = new int[N];
      fTcross        = new double[N];
      fVcross        = new double[N];
      ClearDataArrays();
   }
   //______________________________________________________________________________
   void BaselineCorrection::ClearNZCArrays(){
      for(int i=0;i<fSIZE;i++){
	 fNCrossing[i]     = 0;
	 fCrossingIndex[i] = 0;
	 fTcross[i]        = 0;
	 fVcross[i]        = 0;
      }
   }
   //______________________________________________________________________________
   void BaselineCorrection::ClearDataArrays(){
      for(int i=0;i<fSIZE;i++){
	 fX[i]         = 0;
	 fY[i]         = 0;
	 fEY[i]        = 0;
      }
   }
   //______________________________________________________________________________
   void BaselineCorrection::SetParameters(zcInputParameters par){
      fParameters.offsetOrder       = par.offsetOrder; 
      fParameters.verbosity         = par.verbosity; 
      fParameters.timeThreshold     = par.timeThreshold;
      fParameters.expectedFrequency = par.expectedFrequency; 
      fParameters.sampleFrequency   = par.sampleFrequency; 
      fParameters.minTime           = par.minTime; 
      fParameters.maxTime           = par.maxTime; 
   }
   //______________________________________________________________________________
   int BaselineCorrection::ApplyBaselineCorrection(NMRPulse *aPulse){

      int order = fParameters.offsetOrder;

      if(order<1){
	 std::cout << "[BaselineCorrection]: No offset correction applied." << std::endl;
	 return 0;
      }
      
      if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Doing offset correction..." << std::endl;

      int N          = aPulse->GetNumPoints();
      // double t_start = aPulse->GetTime(0); 
      // double t_end   = aPulse->GetTime(N-1);

      double offset1=0;
      double offset2=0;
      double offset4=0;

      const int SIZE = 2;
      double *offset3 = new double[SIZE];
      double *offset5 = new double[SIZE];

      for(int i=0;i<SIZE;i++){
	 offset3[i] = 0;
	 offset5[i] = 0;
      }

      // test NMRPulse object 
      NMRPulse *myPulse  = new NMRPulse();

      // double t_thr = par.t_thr; 

      // order 1 is ALWAYS applied to "raw" data (aPulse) 
      if(order>=1){
	 offset1 = GetOffset(fParameters.timeThreshold,aPulse);
	 ApplyOffset(offset1,aPulse);
	 if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Offset1 = " << offset1 << std::endl;
      }
      
      // get RMS noise here. 
      // RMSNoise = GetRMSNoise(fParameters.timeThreshold,aPulse);

      // set data in test pulse object: use chunk of data desired in ZC analysis 

      int j=0;
      double time=0,voltage=0,vErr=0;
      for(int i=0;i<N;i++){
	 time = aPulse->GetTime(i);
	 if( time>=fParameters.minTime && time<=fParameters.maxTime){
	    j++;
	 }
      }

      int NN = j;
      j = 0;
 
      myPulse->SetNumPoints(NN);
      for(int i=0;i<N;i++){
	 time    = aPulse->GetTime(i);
	 voltage = aPulse->GetVoltage(i);
	 vErr    = aPulse->GetVoltageErr(i);
	 if( time>=fParameters.minTime && time<=fParameters.maxTime){
	    myPulse->SetDataPoint(j,time,voltage,vErr);
	    j++;
	 }
      }

      // int NCycles = (int)( fParameters.expectedFrequency*(fParameters.maxTime - fParameters.minTime)/10 );

      double T_exp  = 1./fParameters.expectedFrequency;
      double N_exp  = T_exp*fParameters.sampleFrequency;       // number of points for one period 
      int step      = (int)( (1./8.)*N_exp );                  // skip 1/8 of a period 

      if(order>=2){
	 fOffsetFail = false;
	 offset2 = GetOffsetZC(offset1,step,myPulse);
	 ApplyOffset(offset2,myPulse);
	 ApplyOffset(offset2,aPulse);
	 if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Offset2 = " << offset2 << std::endl;
      }
      if(order>=3){
	 fOffsetFail = false;
	 GetOffsetLinear(offset2,myPulse,offset3);
	 ApplyOffsetLinear(offset3,myPulse);
	 ApplyOffsetLinear(offset3,aPulse);
	 if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Offset3: constant = " << offset3[0] << " slope = " << offset3[1] << std::endl;
      }
      if(order>=4){
	 fOffsetFail = false;
	 offset4 = GetOffsetZC(offset3[0],step,myPulse);
	 ApplyOffset(offset4,myPulse);
	 ApplyOffset(offset4,aPulse);
	 if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Offset4 = " << offset4 << std::endl;
      }

      if(order>=5){
	 fOffsetFail = false;
	 GetOffsetLinear(offset4,myPulse,offset5);
	 ApplyOffsetLinear(offset5,myPulse);
	 ApplyOffsetLinear(offset5,aPulse);
	 if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Offset5: constant = " << offset5[0] << " slope = " << offset5[1] << std::endl;
      }

      // VMax = GetVMax(aPulse);

      if(fParameters.verbosity>=2) std::cout << "[BaselineCorrection]: Done." << std::endl;

      delete[] offset3;
      delete[] offset5;
      delete myPulse;

      return 0;
   }
   //______________________________________________________________________________
   double BaselineCorrection::GetOffset(double t_thr,NMRPulse *aPulse){
      std::vector<double> Noise;

      const int N   = aPulse->GetNumPoints();
      double t_max  = aPulse->GetTime(N-1);
      double tstart = t_max - t_thr;
      double tend   = t_max;

      // std::cout << "tmax   = " << t_max  << std::endl; 
      // std::cout << "tstart = " << tstart << std::endl;
      // std::cout << "tend   = " << tend   << std::endl;

      double t=0,v=0;
      for(int i=0;i<N;i++){
	 t = aPulse->GetTime(i);
	 v = aPulse->GetVoltage(i);
	 if( (t>tstart)&&(t<tend) ){
	    Noise.push_back(v);
	 }
      }

      // double rms    = dsp::rms(Noise); 
      double mean   = Utility::GetMean(Noise);
      double offset = mean;
      Noise.clear();

      // std::cout << "NOISE LEVEL: " << mean << std::endl;

      return offset;
   }
   //______________________________________________________________________________
   void BaselineCorrection::GetOffsetLinear(double input_offset,NMRPulse *aPulse,double *offset){

      // for each complete N cycles, compute the offset of the data by considering 
      // zero crossings  

      int counter=0;

      double t=0,v=0,ev=0;
      double v_mean=0,t_mean=0;

      std::vector<double> x,tOff,vOff;

      double NumCycles = (int)( fParameters.expectedFrequency*(fParameters.maxTime - fParameters.minTime)/10 ) ;   // cut down by a factor of 10   
      int NSamples     = (int)( (NumCycles/fParameters.expectedFrequency)*fParameters.sampleFrequency );  // length of time to average over 

      NMRPulse *myPulse = new NMRPulse(0,NSamples);

      if(fParameters.verbosity>=4) std::cout << "Averaging over " << NSamples << " samples..." << std::endl;

      // int dummy=0;
      // cin >> dummy;

      double T_exp  = 1./fParameters.expectedFrequency;
      double N_exp  = T_exp*fParameters.sampleFrequency;       // number of points for one period 
      int step      = (int)( (1./8.)*N_exp );  // skip 1/8 of a period 

      const int N = aPulse->GetNumPoints();
      for(int i=0;i<N;i++){
	 counter++;
	 if(counter<=NSamples){
	    t  = aPulse->GetTime(i);
	    v  = aPulse->GetVoltage(i);
	    ev = aPulse->GetVoltageErr(i);
	    myPulse->SetDataPoint(counter-1,t,v,ev);
	 }else{
	    // std::cout << "counter = " << counter << "\t" << "NSamples = " << NSamples << std::endl;
	    // achieved the desired number of samples 
	    for(int j=0;j<NSamples;j++) x.push_back( myPulse->GetTime(j) );
	    t_mean = Utility::GetMean(x);
	    v_mean = GetOffsetZC(input_offset,step,myPulse);
	    // std::cout << "t_mean = " << t_mean << "\t" << "v_mean = " << v_mean << std::endl;
	    // cin >> dummy;
	    tOff.push_back(t_mean);
	    vOff.push_back(v_mean);
	    // clean up 
	    counter = 0;
	    myPulse->ClearVectorData();
	    x.clear();
	 }
      }

      delete myPulse;

      // Do least squares fit to f(t) = a + bt  
      double a=0,b=0,r=0,rc=0;

      if(fOffsetFail){
	 offset[0] = 0;
	 offset[1] = 0;
      }else{
	 rc = Utility::LeastSquaresFitting(tOff,vOff,a,b,r);
	 if(rc!=0){
	    std::cout << "[BaselineCorrection::GetOffsetLinear]: Linear fit failed! Setting fit coefficients to zero..." << std::endl;
	    offset[0] = 0;
	    offset[1] = 0;
	 }else{
	    offset[0] = a;
	    offset[1] = b;
	 }
      }

   }
   //______________________________________________________________________________
   void BaselineCorrection::ApplyOffset(double offset,NMRPulse *aPulse){
      // apply an offset to the data based on the mean voltage
      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection]: Applying offset of " << offset << " V to data..." << std::endl;

      double t=0,v=0,ev=0;

      int N = aPulse->GetNumPoints();
      for(int i=0;i<N;i++){
	 t  = aPulse->GetTime(i);
	 v  = aPulse->GetVoltage(i) - offset;
	 ev = aPulse->GetVoltageErr(i);
	 aPulse->SetDataPoint(i,t,v,ev);
      }

      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection]: Done." << std::endl;

   }
   //______________________________________________________________________________
   void BaselineCorrection::ApplyOffsetLinear(double *par,NMRPulse *aPulse){
      // apply an offset to the data based on the mean voltage
      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection]: Applying offset to data..." << std::endl;

      double offset=0,t=0,v=0,ev=0;

      int N = aPulse->GetNumPoints();
      for(int i=0;i<N;i++){
	 t      = aPulse->GetTime(i);
	 offset = par[0] + par[1]*t;
	 v      = aPulse->GetVoltage(i) - offset;
	 ev     = aPulse->GetVoltageErr(i);
	 aPulse->SetDataPoint(i,t,v,ev);
      }

      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection]: Done." << std::endl;

   }
   //______________________________________________________________________________
   double BaselineCorrection::GetOffsetZC(double input_offset,int stepSize,NMRPulse *aPulse){

      // find the offset using zero crossings to determine if time between crossings is constant,
      // as is expected for pure sine waves

      // NMRPulse *myPulse = aPulse->Clone();          // FIXME: this is where the memory leak is occuring! 
      NMRPulse *myPulse = new NMRPulse(aPulse);        // FIXME: this is where the memory leak is occuring! 
      // NMRPulse *myPulse = aPulse;                // object initialization; copy constructor called; doesn't really work...

      int NPTS = stepSize/2; 

      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection::GetOffsetZC]: Finding additional offset..." << std::endl;

      if(fParameters.verbosity>=3){
	 std::cout << "[BaselineCorrection::GetOffsetZC]: Parameters: " << std::endl;
	 std::cout << "                                   Points in fit: "  << NPTS << std::endl;
	 std::cout << "                                   Step size:     "  << stepSize << std::endl;
      }

      // settings for counting zero crossings
      bool UseRange  = false;
      double tMin    = 0;    // in seconds 
      double tMax    = 1;    // in seconds                         
      int type       = Utility::kLeastSquares;

      // vector<int> NCrossing,CrossingIndex;
      // vector<double> tCross,vCross;
      // vector<double> T,V,Slope;

      // int NN=0;
      // int is_nan=0; 
      int rc=0,counter=0;

      int NNa = aPulse->GetNumPoints();
      int NNb = myPulse->GetNumPoints();
      if(NNb!=NNa){
	 std::cout << "[BaselineCorrection::GetOffsetZC]: NMRPulse objects don't match! " << std::endl;
	 std::cout << "                               input pulse size: " << NNa << std::endl;
	 std::cout << "                               test pulse size:  " << NNb << std::endl;
	 return 0;
      }

      double t_even=0,t_odd=0;
      double err   = 1E-16;

      double t_diff_abs=0,t_diff_abs_2=0;
      t_diff_abs_2+=0;

      // first calculation 
      int nzc = Utility::CountZeroCrossings(fParameters.verbosity,type,NPTS,stepSize,UseRange,tMin,tMax,fParameters.useT2Time,
	    myPulse,fX,fY,fEY,fNCrossing,fCrossingIndex,fTcross,fVcross);
      double t_diff_old = GetTDiff(nzc,fTcross,t_even,t_odd);
      ClearNZCArrays();

      if(fOffsetFail){
	 delete myPulse;
	 return 0;
      }

      double offset_old  = input_offset;
      double offset_new  = offset_old*(1. + 0.01);  // add 1%  

      if(offset_old==0) offset_new = 1E-6;

      if(fParameters.verbosity>=4){
	 std::cout << "----------------------------------------------------------------" << std::endl;
	 std::cout << "First calculation: " << std::endl;
	 printf("offset_old = %.7E \n",offset_old);
	 printf("t_diff_old = %.7E \n",t_diff_old);
	 printf("offset_new = %.7E \n",offset_new);
	 std::cout << "----------------------------------------------------------------" << std::endl;
      }

      ApplyOffset(offset_new,myPulse);

      nzc = Utility::CountZeroCrossings(fParameters.verbosity,type,NPTS,stepSize,UseRange,tMin,tMax,fParameters.useT2Time,
	    myPulse,fX,fY,fEY,fNCrossing,fCrossingIndex,fTcross,fVcross);
      double t_diff_new = GetTDiff(nzc,fTcross,t_even,t_odd);
      ClearNZCArrays();

      if(fOffsetFail){
	 delete myPulse;
	 return 0;
      }

      // reset the pulse 
      // myPulse = aPulse->Clone();
      delete myPulse;
      myPulse = new NMRPulse(aPulse);

      double slope = (t_diff_new - t_diff_old)/(offset_new - offset_old);

      rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope);
      if(rc>0){
	 delete myPulse;
	 return 0;
      }

      if(fParameters.verbosity>=4){
	 std::cout << "----------------------------------------------------------------" << std::endl;
	 std::cout << "Second calculation: " << std::endl;
	 printf("offset_old = %.7E \n",offset_old);
	 printf("t_diff_old = %.7E \n",t_diff_old);
	 printf("offset_new = %.7E \n",offset_new);
	 printf("t_diff_new = %.7E \n",t_diff_new);
	 printf("slope      = %.7E \n",slope)     ;
	 std::cout << "----------------------------------------------------------------" << std::endl;
      }

      offset_new = offset_old - t_diff_old/slope;
      double root_diff  = fabs(offset_new-offset_old);

      rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope);
      if(rc>0){
	 delete myPulse;
	 return 0;
      }

      if(t_diff_new==0) offset_new = offset_old;          // if new time difference is identically zero, we use old offset.
      if(fParameters.verbosity>4){
	 std::cout << "trial     = 0"                   << std::endl;
	 printf("slope      = %.7E \n",slope);
	 printf("offset_old = %.7E   ",offset_old);
	 printf("t_diff_old = %.7E \n",t_diff_old);
	 printf("offset_new = %.7E \n",offset_new);
	 printf("diff(offset-offset_prev) = %.7E \n",root_diff);
	 std::cout << "----------------------------------------------------------------" << std::endl;
      }

      // we can't have a massive offset -- that's just not true. 
      if( fabs(offset_new)>0.005 ) offset_new = offset_old;

      // update values 
      offset_old = offset_new;
      t_diff_old = t_diff_new;

      // clear arrays before starting
      ClearNZCArrays();

      do{
	 offset_new = offset_old - t_diff_old/slope;
	 // check the new offset  
	 rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope);
	 if(rc>0) break;
	 ApplyOffset(offset_new,myPulse);
	 nzc = Utility::CountZeroCrossings(fParameters.verbosity,type,NPTS,stepSize,UseRange,tMin,tMax,fParameters.useT2Time,
	       myPulse,fX,fY,fEY,fNCrossing,fCrossingIndex,fTcross,fVcross);
	 t_diff_new = GetTDiff(nzc,fTcross,t_even,t_odd);
	 slope      = (t_diff_new - t_diff_old)/(offset_new - offset_old);
	 root_diff  = fabs(offset_new - offset_old);
	 if(fParameters.verbosity>4){
	    std::cout << "trial     = " << counter+1 << std::endl;
	    printf("slope      = %.7E \n",slope);
	    printf("offset_old = %.7E   ",offset_old);
	    printf("t_diff_old = %.7E \n",t_diff_old);
	    printf("offset_new = %.7E   ",offset_new);
	    printf("t_diff_new = %.7E \n",t_diff_new);
	    printf("diff(offset-offset_prev) = %.7E \n",root_diff);
	    std::cout << "----------------------------------------------------------------" << std::endl;
	 }
	 if(fOffsetFail) break;
	 // set up for next calc 
	 // myPulse = aPulse->Clone();       // reset the pulse   
	 delete myPulse;
	 myPulse = new NMRPulse(aPulse);
	 t_diff_abs   = fabs(t_diff_new);
	 t_diff_abs_2 = fabs(t_diff_new-t_diff_old);
	 t_diff_old   = t_diff_new;
	 offset_old   = offset_new;
	 ClearNZCArrays();
	 counter++;
      }while( (t_diff_abs>err)&&(counter<20) );

      if(counter==20){
	 if(fParameters.verbosity>=3){
	    std::cout << "Reached max counter (20).  Giving up." << std::endl;
	    printf("offset_old: %.7E \n",offset_old);
	    printf("offset_new: %.7E \n",offset_new);
	    printf("t_diff_old: %.7E \n",t_diff_old);
	    printf("t_diff_new: %.7E \n",t_diff_new);
	    printf("slope:      %.7E \n",slope     );
	 }
      }

      if(fParameters.verbosity>=3) std::cout << "[BaselineCorrection::GetOffsetZC]: Done." << std::endl;

      delete myPulse;
      ClearNZCArrays();

      if(rc>0 || fOffsetFail==true){
	 return 0;
      }else{
	 return offset_new;
      }

   }
   //______________________________________________________________________________
   double BaselineCorrection::GetTDiff(int nzc,double *tCross,double &delta_t_even_nc,double &delta_t_odd_nc){
      delta_t_odd_nc=0;
      delta_t_even_nc=0;
      int counter_odd=0,counter_even=0;
      if(nzc<2){
	 if(fParameters.verbosity>=4) std::cout << "[BaselineCorrection::GetTDiff]: NOT ENOUGH DATA POINTS!" << std::endl;
	 if(fParameters.verbosity>=4) std::cout << "[BaselineCorrection::GetTDiff]: Number of data points: " << nzc << std::endl;
	 fOffsetFail = true;
	 return -1;
      }else{
	 for(int i=1;i<nzc;i++){
	    if(i%2==0){
	       delta_t_odd_nc  += (tCross[i]-tCross[i-1]);   // we're going by even/odd Zc -- not Nc, so we interchange odd/even wrt Nc. 
	       counter_odd++;
	    }else if(i%2!=0){
	       delta_t_even_nc += (tCross[i]-tCross[i-1]);
	       counter_even++;
	    }
	 }
	 delta_t_odd_nc  /= ( (double)counter_odd );
	 delta_t_even_nc /= ( (double)counter_even);
      }
      double t_diff = delta_t_odd_nc-delta_t_even_nc;
      return t_diff;
   }
   //______________________________________________________________________________
   int BaselineCorrection::CheckOffset(double offset_old,double offset_new,double t_diff_old,double t_diff_new,double slope){

      int is_nan_t_diff_old = std::isnan(t_diff_old);
      int is_nan_t_diff_new = std::isnan(t_diff_new);
      int is_nan_offset_old = std::isnan(offset_old);
      int is_nan_offset_new = std::isnan(offset_new);
      int is_nan_slope      = std::isnan(slope);

      int rc = 0,rc_tot=0;
      rc+=0;

      if(is_nan_t_diff_old){
	 rc = 1; 
	 rc_tot++;
      }
      if(is_nan_t_diff_new){
	 rc = 2; 
	 rc_tot++;
      }
      if(is_nan_offset_old){
	 rc = 3;
	 rc_tot++;
      }
      if(is_nan_offset_new){
	 rc = 4; 
	 rc_tot++;
      }
      if(is_nan_slope){
	 rc = 5; 
	 rc_tot++;
      }

      // int dummy=0; 
      if(rc_tot>0){
	 printf("[BaselineCorrection]: WARNING: One of the offset values below is NAN! \n");
	 printf("                  offset_old: %.7E \n",offset_old);
	 printf("                  offset_new: %.7E \n",offset_new);
	 printf("                  t_diff_old: %.7E \n",t_diff_old);
	 printf("                  t_diff_new: %.7E \n",t_diff_new);
	 printf("                  slope:      %.7E \n",slope     );
	 // printf("Enter any number to continue: ");
	 // cin  >> dummy; 
      }

      return rc_tot;

   }

}
