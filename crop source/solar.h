// Class CSolar
//PFD Photon Flux Density
//tau atmospheric transmissivity, HGJones 1991 considered a constant, it can be calculated from
// measured radiation/potential radiation at atmospheric surface
//SolRad Solar Radiation
// PARfr
//UseObs UseTau
//PAR Photosynthetically Active Radiation (what is PFD)
// Altitude
//time - time of day
//CosTheta Zenith angle see \eqn 11.1 Campbell and Norman (1998)
// consider using equation from GLYCIM elevation angle = 90 - zenith angle since the are measured
// differently
#pragma once
#define cPFD 4.6f  // conversion factor from PAR (W/m2) to PFD (umol m-2 s) for solar radiation, see Campbell and Norman (1994) p 149
// 4.6 is a conversion factor from W to photons for visible solar radiation.
// Amthor 1994, McCree 1981, Challa 1995, Campbell and Norman 1998.
// Some use 4.55, see Goudriaan and van Laar (1994)
#define     SolarConst 1367.0f   // solar constant, Iqbal (1983) Watts m-2
#define     FDIV_GUARD 1.0e-8f // floating point divison guard (tolerance)
#define     PI 3.1415926f

class CSolar
{
private:
	int JDay;
// Environmental Components
	float Time, Latitude, Longitude, altitude, tau;
	float DayLength, Declination, SolarNoon, SinElevation, Azimuth;
	float Sunrise, Sunset,HalfDay, Elevation, CosElevation;
	float CosTheta;
// solar components
	float PAR, SolarRadiation, PotentialSolarTotal, PotentialSolarDirect, PotentialSolarDiffuse;
	float NIRTotal, NIRDiffuse, NIRDirect, PARTotal, PARDiffuse, PARDirect;
	float PFD, NIR, PARFraction, NIRFraction;
	float PotentialPARDiffuse, PotentialPARDirect, PotentialPARTotal;				//PAR: photosynthesis Active Rad
	float PotentialNIRDiffuse, PotentialNIRDirect, PotentialNIRTotal;				//NIR: Near-Infrared Light Rad
	float SolarDirect, SolarDiffuse, FracDiffuse, PFDDirect, PFDDiffuse;			//PFD: photosynthetic flux density (umol photons m-2 s-1)
	float FracPARDirect, FracPARDiffuse, FracNIRDirect, FracNIRDiffuse, FracSolarDiffuse, FracSolarDirect;
	float FracPFDDirect, FracPFDDiffuse;
	float RowAzimuth;																//Angle of row orientation measured from 180?

	bool useObs, useTau;
	static const float SDERP[9]; 
	//used for declination calcs
	// To do: incorporate the computation of global irradiance including NIR
	// PAR is in Wm-1, PFD is in umol m-2 s-1 for visible wavebands

	void SetDayLength();
	void SetDeclination();// Solar declination as a f(Jday)
	void SetSolarNoon();
	void SetSolarElevation() ;    // Solar height (=elevation)
	void SetAzimuth() ; // Solar azimuth measured from ?
	float press();
	float m();

    void SetPotentialSolar();
	void SetPotentialPAR();
	void SetPotentialNIR();

	void SetNIRTotal();
    void SetNIRDirect();
	void SetNIRDiffuse();

    void SetPARTotal();
    void SetPARDirect();
	void SetPARDiffuse();

	void SetFracPARDirect();
	void SetFracNIRDirect();
	void SetPARFraction(float Fraction);	// when measured PAR and Radiation are both available
	void SetPARFraction();			//When only radiation or PAR are available


public:
	CSolar(void);
	~CSolar(void);
	void SetVal(int Day, float Time, float Lat, float Longi, float Alti, float SolRad0);

	float GetDayLength() { return DayLength; }
	float GetDeclination() { return Declination; }
	float GetSolarNoon() { return SolarNoon; }
	float GetSunrise() { return Sunrise; }
	float GetSunset() { return Sunset; }
	float GetSinElevation() { return SinElevation; }
	float GetSolarElevation() { return Elevation; }
	float GetAzimuth() { return Azimuth; }

	float GetPotentialSolarDiffuse() { return PotentialSolarDiffuse; }
	float GetPotentialSolarTotal() { return PotentialSolarTotal; }
	float GetPotentialSolarDirect() { return PotentialSolarDirect; }
	float GetSolarRadiation() { return SolarRadiation; }

	float GetPAR() { return PAR; }
	float GetPARFraction() { return PARFraction; }
	float GetFracPARDirect() { return FracPARDirect; }
	//float GetFracPARDiffuse() { return 1.0 - FracPARDirect; }
	float GetFracNIRDirect() { return FracNIRDirect; }
	//float GetFracNIRDiffuse() { return 1.0 - FracNIRDirect; }
	//float GetNIRFraction() { return 1.0 - PARFraction; }
	float GetPFD() { return PAR * cPFD; }
	float GetPFDDiffuse() { return PAR * cPFD * (1.0f - FracPARDirect); }
	float GetPFDDirect() { return PAR * cPFD * FracPARDirect; }
	float GetPFDTotal() { return PAR * cPFD; } // sane as GetPFD
	float GetNIR() { return NIR; }
	float GetNIRTotal() { return NIRTotal; }
	float GetNIRDiffuse() { return NIRDiffuse; }
	float GetNIRDirect() { return NIRDirect; }

};
