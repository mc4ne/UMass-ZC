#include "UMZeroCrossing.hh"

namespace UMass { 
   ZCTest::ZCTest(){
     ;
   }
   //______________________________________________________________________________
   ZeroCrossing::ZeroCrossing(zcInputParameters par){
      Init();
      SetParameters(par);
   }
   //______________________________________________________________________________
   ZeroCrossing::~ZeroCrossing(){
      delete[] fX; 
      delete[] fY; 
      delete[] fEY;
      delete[] fFREQ;
      delete[] fFREQ_ph;
      delete[] fZC; 
      delete[] fNC;
      delete[] fNCrossing;
      delete[] fCrossingIndex;
      delete[] fTcross;
      delete[] fVcross;
      delete[] fFreqAtCrossing; 
      delete[] fNumCycles;
      delete fPulse; 
      delete fBaseline; 
   }
   //______________________________________________________________________________
   void ZeroCrossing::Init(){
      fNPTS             = 1; 
      fNPTSUseable      = 0; 
      fStep             = 0; 
      fNZC              = 0;
      const int N = 3;
      fZC      = new int[N];
      fFREQ    = new double[N];
      fFREQ_ph = new double[N];
      fNC      = new double[N];
      for(int i=0;i<N;i++){
	 fFREQ[i]    = 0;  
	 fFREQ_ph[i] = 0;  
	 fZC[i]      = 0;  
	 fNC[i]      = 0; 
      } 
      const int NN = fNPTS;
      fX  = new double[NN]; 
      fY  = new double[NN]; 
      fEY = new double[NN];
      ClearAnaArrays();
      const int MAX = MAX_SIZE; 
      fPulse          = new NMRPulse(0,MAX); 
      fNCrossing      = new int[MAX]; 
      fCrossingIndex  = new int[MAX]; 
      fTcross         = new double[MAX]; 
      fVcross         = new double[MAX];
      fFreqAtCrossing = new double[MAX];  
      fNumCycles      = new double[MAX];
      ClearVectors();
      fBaseline       = new BaselineCorrection();
   }
   //______________________________________________________________________________
   void ZeroCrossing::SetParameters(zcInputParameters par){
      // fParameters.expectedFrequency     = par.expectedFrequency;
      // fParameters.sampleFrequency       = par.sampleFrequency;
      // fParameters.minTime               = par.minTime;
      // fParameters.maxTime               = par.maxTime; 
      // fParameters.device                = par.device; 
      // fParameters.verbosity             = par.verbosity; 
      // fParameters.debug                 = par.debug;  
      // fParameters.useMidpoint           = par.useMidpoint; 
      // fParameters.useLinearInterp       = par.useLinearInterp;
      // fParameters.useLeastSquares       = par.useLeastSquares; 
      // fParameters.useTimeRange          = par.useTimeRange;   
      // fParameters.useIntegerCycles      = par.useIntegerCycles; 
      // fParameters.useT2Time             = par.useT2Time; 
      // fParameters.useBaselineCorrection = par.useBaselineCorrection; 
      fParameters = par; 
      UpdateParameters(); 
      // update the baseline values as well
      fBaseline->SetParameters(par);  
   }
   //______________________________________________________________________________
   void ZeroCrossing::UpdateParameters(){
      // Freq_exp   = expected frequency
      // SampleFreq = Frequency of ADC clock
      // step_size  = number of points to skip after finding a zero crossing
      // NPTS       = number of points to use in fitting for zero crossing 

      double T_exp  = 1./fParameters.expectedFrequency;
      double N_exp  = T_exp*fParameters.sampleFrequency;  // number of points for one period 
      int step_size = (int)( (1./16.)*N_exp );            // skip 1/16 of a period 
      int NPTS      = step_size/2;                        // use step_size/2 
      //Make sure number of points in fit is at least 2
      if (NPTS<2){
	NPTS = 2;
      }

      SetStepSize(step_size); 
      SetNumPointsForFits(NPTS); 

      if(fParameters.verbosity>=3){ 
	 printf("[ZeroCrossing::UpdateParameters]: Expected Frequency: %.7lf \n",fParameters.expectedFrequency   );
	 printf("[ZeroCrossing::UpdateParameters]: Sampling Frequency: %.7lf \n",fParameters.sampleFrequency);
	 printf("[ZeroCrossing::UpdateParameters]: Step size:          %d    \n",step_size  );
	 printf("[ZeroCrossing::UpdateParameters]: Number of points:   %d    \n",NPTS       );
      }

   }
   //______________________________________________________________________________
   double ZeroCrossing::GetFrequency(int type) const{
      double freq=0;
      switch(type){ 
	 case Utility::kMidpoint:
	    freq = fFREQ[0];
	    break;
	 case Utility::kLinearInterpolation:
	    freq = fFREQ[1];
	    break;
	 case Utility::kLeastSquares:
	    freq = fFREQ[2];
	    break;
	 case Utility::kMidpointPhase:
	    freq = fFREQ_ph[0];
	    break;
	 case Utility::kLinearInterpolationPhase:
	    freq = fFREQ_ph[1];
	    break;
	 case Utility::kLeastSquaresPhase:
	    freq = fFREQ_ph[2];
	    break;
	 default: 
	    std::cout << "[ZeroCrossing::GetFrequency]: Invalid frequency type" << std::endl;
	    freq = -1;
      }
      return freq;
   }
   //______________________________________________________________________________
   int ZeroCrossing::Analyze(std::vector<double> ampl,std::vector<double> &result){

      // ampl   = vector of voltage samples 
      // result = frequency results.  order is: mid, lin, lsq, mid (ph), lin (ph), lsq (ph)   

      int rc=0;

      if(fParameters.sampleFrequency==0){
	 std::cout << "[ZeroCrossing::Analyze]: Invalid sampling frequency = " << fParameters.sampleFrequency << std::endl;
	 return 1;
      }

      if(fParameters.expectedFrequency==0){
	 std::cout << "[ZeroCrossing::Analyze]: Invalid expected frequency = " << fParameters.expectedFrequency << std::endl;
	 return 1;
      }

      // build the time vector 
      int N = ampl.size();
      double t_arg=0;
      std::vector<double> time; 
      for(int i=0;i<N;i++){
         t_arg = ( (double)i )/fParameters.sampleFrequency; 
	 time.push_back(t_arg);
      }

      // convert to voltage if necessary 
      if(fParameters.device==Utility::kPlungingProbe){
	 Utility::CutPlungingProbeData(time,ampl); 
         Utility::ConvertToVoltage(Utility::k3316,ampl);  
      }

      if(fParameters.device==Utility::kFixedProbe){
	 Utility::ConvertToVoltage(Utility::k3316,ampl);  
      }

      // add the data to our NMRPulse object
      const int NN = time.size();
      fPulse->SetNumPoints(NN); 
      for(int i=0;i<NN;i++){
	 fPulse->SetDataPoint(i,time[i],ampl[i],0.);
      }

      // print to file 
      char outpath[500];
      if(fParameters.debug){
	 sprintf(outpath,"./umass-zc_trace.csv"); 
	 Utility::PrintSignalToFile(outpath,fPulse);
      }

      // do the baseline correction if necessary 
      if(fParameters.useBaselineCorrection){
	 fBaseline->ApplyBaselineCorrection(fPulse); 
      }

      // if we're using the T2 time as the endpoint, calculate it here
      // double t2Time = Utility::GetT2Time_old(fPulse);
      double t2Time = Utility::GetT2Time_v3a(0,fPulse,1); // start from index zero, last parameter is verbosity
      fPulse->SetT2Time(t2Time);

      char msg[200];  

      if(fParameters.useT2Time){
	 fParameters.maxTime = t2Time; 
	 if(fParameters.verbosity>=2 || fParameters.debug==true){
	    sprintf(msg,"[ZeroCrossing::Analyze]: Using the T2 time = %.3lf ms",fParameters.maxTime/1E-3); 
	    std::cout << msg << std::endl;
	 }
      }

      // print to file 
      if(fParameters.debug){
	 sprintf(outpath,"./umass-zc_trace_blc.csv"); 
	 Utility::PrintSignalToFile(outpath,fPulse);
      } 

      // run the calculation and store results 
      rc = Calculate(fPulse); 

      for(int i=0;i<3;i++) result.push_back(fFREQ[i]   );   
      for(int i=0;i<3;i++) result.push_back(fFREQ_ph[i]); 

      // clear NMRPulse for next iteration 
      fPulse->ClearData();  

      return rc; 
   }
   //______________________________________________________________________________
   void ZeroCrossing::SetNumPointsForFits(int n){
      fNPTS = n;
      InitAnaArrays();
   }
   //______________________________________________________________________________
   void ZeroCrossing::InitAnaArrays(){
      // analysis arrays: data points to use in finding a zero crossing
      delete[] fX;  
      delete[] fY;  
      delete[] fEY;  
      const int N = fNPTS;
      if(N>0){
	 fX  = new double[N]; 
	 fY  = new double[N]; 
	 fEY = new double[N];
	 ClearAnaArrays();  
      }else{
	 std::cout << "[ZeroCrossing::InitAnaArrays]: Cannot initialize arrays!  Is NPTS set?" << std::endl;
	 exit(1);
      }
   }
   //______________________________________________________________________________
   void ZeroCrossing::Reset(){
      // resets all data members (useful for multiple zero crossing calculations) 
      ClearAnaArrays();
      ClearVectors();
   }
   //______________________________________________________________________________
   void ZeroCrossing::ClearAnaArrays(){
      for(int i=0;i<fNPTS;i++){
	 fX[i]  = 0; 
	 fY[i]  = 0; 
	 fEY[i] = 0; 
      }
   }
   //______________________________________________________________________________
   void ZeroCrossing::ClearVectors(){
      for(int i=0;i<MAX_SIZE;i++){
	 fNCrossing[i]      = 0;
	 fCrossingIndex[i]  = 0; 
	 fTcross[i]         = 0; 
	 fVcross[i]         = 0;
	 fFreqAtCrossing[i] = 0;  
	 fNumCycles[i]      = 0; 
      }
   }
   //______________________________________________________________________________
   int ZeroCrossing::Calculate(NMRPulse *aPulse){

      InitAnaArrays(); 

      bool useNLLSQ = fParameters.useNonLinearLSQ; 
      int fitFunc   = fParameters.phaseFitFunc; 

      // int PulseNumber = aPulse->GetPulseNumber(); 
      int rc_fr=0,rc=0; 
      int zc_mid=0,zc_lin=0,zc_lsq=0;
      double nc_mid=0,nc_lin=0,nc_lsq=0; 
      double freq_mid=0,freq_lin=0,freq_lsq=0;
      double freq_mid_ph=0,freq_lin_ph=0,freq_lsq_ph=0;  // fit num cycles vs t_zc to a line 

      if(fParameters.useMidpoint){
	 CountZeroCrossings(UMass::Utility::kMidpoint,aPulse);
	 rc_fr       = CalculateFrequencies(zc_mid,freq_mid);
         if(useNLLSQ){
	    freq_mid_ph = GetFrequencyFromPhaseFit_nllsq(fitFunc);
         }else{
	    freq_mid_ph = GetFrequencyFromPhaseFit();
         } 
	 nc_mid      = ( (double)zc_mid - 1.)/2.; 
	 // PrintVectorData(1,PulseNumber); 
	 Reset(); 
      }
      rc += rc_fr; 
      if(fParameters.useLinearInterp){
	 CountZeroCrossings(UMass::Utility::kLinearInterpolation,aPulse);
	 rc_fr       = CalculateFrequencies(zc_lin,freq_lin);
	 freq_lin_ph = GetFrequencyFromPhaseFit(); 
         if(useNLLSQ){
	    freq_lin_ph = GetFrequencyFromPhaseFit_nllsq(fitFunc);
         }else{
	    freq_lin_ph = GetFrequencyFromPhaseFit();
         } 
	 nc_lin      = ( (double)zc_lin - 1.)/2.; 
	 // PrintVectorData(2,PulseNumber); 
	 Reset(); 
      }
      rc += rc_fr; 
      if(fParameters.useLeastSquares){
	 CountZeroCrossings(UMass::Utility::kLeastSquares,aPulse);
	 rc_fr       = CalculateFrequencies(zc_lsq,freq_lsq);
         if(useNLLSQ){
	    freq_lsq_ph = GetFrequencyFromPhaseFit_nllsq(fitFunc);
         }else{
	    freq_lsq_ph = GetFrequencyFromPhaseFit();
         } 
	 nc_lsq      = ( (double)zc_lsq - 1.)/2.; 
	 // PrintVectorData(3,PulseNumber); 
	 Reset(); 
      }
      rc += rc_fr; 

      // store results

      // check if not a number 
      if( !std::isnan(freq_mid) )    fFREQ[0] = freq_mid; 
      if( !std::isnan(freq_lin) )    fFREQ[1] = freq_lin; 
      if( !std::isnan(freq_lsq) )    fFREQ[2] = freq_lsq; 

      if( !std::isnan(freq_mid_ph) ) fFREQ_ph[0] = freq_mid_ph; 
      if( !std::isnan(freq_lin_ph) ) fFREQ_ph[1] = freq_lin_ph;
      if( !std::isnan(freq_lsq_ph) ) fFREQ_ph[2] = freq_lsq_ph;

      fZC[0]   = zc_mid;    
      fZC[1]   = zc_lin;    
      fZC[2]   = zc_lsq;    

      fNC[0]   = nc_mid;    
      fNC[1]   = nc_lin;    
      fNC[2]   = nc_lsq;    

      Reset(); 

      return rc; 

   }
   //______________________________________________________________________________
   void ZeroCrossing::CountZeroCrossings(int method,NMRPulse *aPulse){
      fNZC = Utility::CountZeroCrossings(fParameters.verbosity,method,
                                         fNPTS,fStep,fParameters.useTimeRange,fParameters.minTime,fParameters.maxTime,
					 fParameters.useT2Time,
	                                 aPulse,fX,fY,fEY,fNCrossing,fCrossingIndex,fTcross,fVcross);

      if(fParameters.debug){
	 Utility::PrintArraysToFile("zero-crossings.csv",fNZC,fTcross,fVcross,fNCrossing); 
      }

   }
   //______________________________________________________________________________
   int ZeroCrossing::CalculateFrequencies(int &TrueNumCrossings,double &FreqFullRange){

      if(fParameters.verbosity>=3) std::cout << "[ZeroCrossing]: Calculating frequencies..." << std::endl;

      int ret_code = 0; 

      int SIZE = fNZC; // fNCrossing.size();  
      int NumCrossings = SIZE;

      if(NumCrossings==0) std::cout << "[ZeroCrossing::CalculateFrequencies]: NumCrossings is zero!" << std::endl;

      if(fParameters.useIntegerCycles){
	 /// if Zc is even, then number of cycles is odd; we want even cycles so we step back 1 bin
	 if(NumCrossings%2==0){
	    if(fParameters.verbosity>=3) std::cout << "[ZeroCrossing::CalculateFrequencies]: NumCrossings = " << NumCrossings << std::endl;
	    NumCrossings -= 1;
	    if(fParameters.verbosity>=3) std::cout << "[ZeroCrossing::CalculateFrequencies]: Adjusted to NumCrossings = " << NumCrossings << std::endl;
	 }
      }

      TrueNumCrossings = NumCrossings;

      double Zc = (double)NumCrossings;  
      double NC = (Zc - 1.)/2.;  

      if(NumCrossings<5){
	 std::cout << "[ZeroCrossings::CalculateFrequencies]: ERROR!  Number of zero crossings is less than 5! " << std::endl;
	 FreqFullRange = 0; 
	 return 1;  
      }

      // find frequency using all zero crossings 
      // find the index of last crossing  
      int LastCrossingIndex  = NumCrossings - 1;  // the -1 is to account for arrays starting at 0 and not 1

      // compute the time 
      double delta_t  = fTcross[LastCrossingIndex] - fTcross[0];      
      // compute the frequency  
      FreqFullRange = NC/delta_t;   

      if(fParameters.verbosity>=4){
	 printf("[ZeroCrossing::CalculateFrequencies]: time at first crossing: %.7f us \n",1E+6*fTcross[0]                );
	 printf("[ZeroCrossing::CalculateFrequencies]: time at last crossing:  %.7f us \n",1E+6*fTcross[LastCrossingIndex]);
      }

      double iNc=0,ifreq=0; 

      // now calculate the frequency using each crossing  
      for(int i=0;i<NumCrossings;i++){
	 iNc        = ( (double)fNCrossing[i] - 1. )/2.;  
	 if(i==0){
	    delta_t = 0; 
	    ifreq   = 0;
	 }else{
	    delta_t = fTcross[i] - fTcross[0]; 
	    ifreq   = iNc/delta_t; 
	 }
	 fNumCycles[i]      = iNc;
	 fFreqAtCrossing[i] = ifreq; 
      }

      if(fParameters.verbosity>=4) printf("[ZeroCrossing::CalculateFrequencies]: NumCycles =   %.3f \t delta_t = %.7f us \t f = %.7f Hz \n",NC,1E+6*delta_t,FreqFullRange); 

      return ret_code; 

   }
   //______________________________________________________________________________
   double ZeroCrossing::GetFrequencyFromPhaseFit(){
      // linear least squares 
      double freq=0,intercept=0,r=0;
      Utility::LeastSquaresFitting(fNZC,fTcross,fNumCycles,intercept,freq,r); 
      return freq;
   }
   //______________________________________________________________________________
   double ZeroCrossing::GetFrequencyFromPhaseFit_nllsq(const int fitFunc){
      // nonlinear least squares
      int np=0;
      if(fitFunc==1) np = 2; 
      if(fitFunc==3) np = 3; 
      if(fitFunc==5) np = 4; 
      if(fitFunc==7) np = 5; 
      if(fitFunc==9) np = 6; 

      const int NPAR = np; 
 
      std::vector<double> par,parErr;
      for(int i=0;i<NPAR;i++){
	 par.push_back(1); 
	 parErr.push_back(0); 
      }
      const int NPTS = fNZC;
      std::vector<double> x,y,dy;
      for(int i=0;i<NPTS;i++){
	 x.push_back( fTcross[i] );
	 y.push_back( fNumCycles[i] );  
	 dy.push_back(0);               // FIXME: Accurate error estimate?  
      }
      double freq=0; // the frequency is the p1 term 
      int rc=-1; 
      if(fitFunc==1) rc = Utility::NonLinearLeastSquaresFitting(x,y,dy,Utility::poly1,Utility::poly1_df,par,parErr,NPAR,0);
      if(fitFunc==3) rc = Utility::NonLinearLeastSquaresFitting(x,y,dy,Utility::poly3,Utility::poly3_df,par,parErr,NPAR,0);
      if(fitFunc==5) rc = Utility::NonLinearLeastSquaresFitting(x,y,dy,Utility::poly5,Utility::poly5_df,par,parErr,NPAR,0);
      if(fitFunc==7) rc = Utility::NonLinearLeastSquaresFitting(x,y,dy,Utility::poly7,Utility::poly7_df,par,parErr,NPAR,0);
      if(fitFunc==9) rc = Utility::NonLinearLeastSquaresFitting(x,y,dy,Utility::poly9,Utility::poly9_df,par,parErr,NPAR,0);
      if(rc!=0) freq = -1; 
      freq = par[1]; // the frequency is the p1 term 

      // print to file 
      char fileName[200];
      if(fParameters.debug){
	 sprintf(fileName,"fit-func-%d_phase-fit-parameters.csv",fitFunc); 
	 rc = Utility::PrintVectorsToFile(fileName,par,parErr);
      } 

      return freq;
   }
} //::UMass
