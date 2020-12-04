#ifndef NMRPULSE_HH
#define NMRPULSE_HH 

// a generic NMR pulse 

#include <cstdlib> 
#include <cstdio> 
#include <iostream>

namespace UMass {  

   class NMRPulse{

      protected:                    // can access these if one inherits from NMRPulse 
	 int fVerbosity; 
	 int fNumPoints; 
	 int fPulseNumber;
	 double fAmpl,fNoiseRMS,fSNR;
	 double fTemperature,fXCoordinate,fYCoordinate,fZCoordinate;  
	 double *fTime,*fVoltage,*fVoltageErr;
	 unsigned long long fTimeStamp;

      public:
	 // constructor and destructor 
	 NMRPulse(int PulseNumber=0,int NPTS=1);
	 ~NMRPulse();

	 // copy constructors 
	 NMRPulse(const NMRPulse &aPulse);         
	 NMRPulse(const NMRPulse *aPulse);         

	 // overload assignment operator  
	 NMRPulse& operator=(const NMRPulse &aPulse);  
	 NMRPulse* operator=(const NMRPulse *aPulse); 

	 void Print();  
	 void ClearData();
	 void ClearVectorData();
	 void SetNumPoints(int p);   
	 void SetData(int N,double t[],double v[],double ev[]); 
	 void SetDataPoint(int i,double t,double v,double ev); 
	 void SetPulseNumber(int p)                             {fPulseNumber = p;     } 
	 void SetTimeStamp(unsigned long long t)                {fTimeStamp   = t;     } 
	 void SetTemperature(double temp)                       {fTemperature = temp;  } 
	 void SetXCoordinate(double x)                          {fXCoordinate = x;     } 
	 void SetYCoordinate(double y)                          {fYCoordinate = y;     } 
	 void SetZCoordinate(double z)                          {fZCoordinate = z;     } 
	 void SetAmplitude(double a)                            {fAmpl        = a;     }
	 void SetNoiseRMS(double n)                             {fNoiseRMS    = n;     }
	 void SetSignalToNoiseRatio(double snr)                 {fSNR         = snr;   }
	 void SetVerbosity(int v)                               {fVerbosity   = v;     } 

	 int GetNumPoints()                               const {return fNumPoints;    } 
	 int GetPulseNumber()                             const {return fPulseNumber;  }  

	 double GetTemperature()                          const {return fTemperature;  } 
	 double GetXCoordinate()                          const {return fXCoordinate;  } 
	 double GetYCoordinate()                          const {return fYCoordinate;  } 
	 double GetZCoordinate()                          const {return fZCoordinate;  } 

	 double GetAmplitude()                            const {return fAmpl;         }
	 double GetNoiseRMS()                             const {return fNoiseRMS;     }
	 double GetSignalToNoiseRatio()                   const {return fSNR;          }
	 double GetTime(int i)                            const {return fTime[i];      } 
	 double GetVoltage(int i)                         const {return fVoltage[i];   } 
	 double GetVoltageErr(int i)                      const {return fVoltageErr[i];} 

	 unsigned long long GetTimeStamp()                const {return fTimeStamp;    } 

   };

} //::UMass  

#endif 