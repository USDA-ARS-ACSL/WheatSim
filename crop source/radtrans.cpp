/*unit CanopyRadTrans;
{Basic canopy architecture parameters, 10/10/00 S.Kim
 modified to represent heterogeneous canopies
 Uniform continuouse canopy: Width2 = 0
 Hedgerow canopy : Width1 = row width, Width2 = interrow, Height1 = hedge height, Height2 = 0
 Intercropping canopy: Height1, Width1, LA1 for Crop1, and so on
 Rose bent canopy: Height1=Upright canopy, Height2 = bent portion height, 10/16/02 S.Kim}

 This is the interface into the Solar Class but calculates transmission coefficients. To use it you must
 initialize the Solar object. See Solar.cpp for more information.

 This routine adds calculations that require use of LAI and LAF (Leaf Angle Factor)
 You must pass parameters to the Solar Class to initialize that class and use the methods therin.
 The LeafAngleFactor can be numeric in which case it is used directly in one equation to calculate
 Kd as a function of zenith angle (eqn 15.4 in Campbell and Norman). In this case use 'IsLeafAngleFactorUsed' as true
 Otherwise you can use shapes as an enumeration.
*/
// dateutils ;

#include <cmath>  // need for math functions
#include <algorithm>  // need for max min functions
#include "radtrans.h"

using namespace std;

inline static float cot(float a) { return 1 / tan(a); }
inline static float sqr(float a) { return (a * a); }


// include Math unit for Real mathmatical doubles

CRadTrans::CRadTrans()
{
//Z LeafAngle needs to be updated, while the "IsLeafAngleFactorUsed" indicates it is not used here
//Z still recommend to be updated, since is updated later.
//  directly use "LeafAngleFactor", see Table 15.1 in Campbell and Normann 1998, see "wheat"
	LeafAngle = Spherical;
	LeafAngleFactor = 0.96f; //Z average leaf angle factor for rye grass
	IsLeafAngleFactorUsed = false;
	absorp = 0.85f;   //leaf absorptivity for PAR
	clump = 1.0f;
	rho_soil = 0.10f; // soil reflectivity for PAR band

//Z get other variables initial values, such that
	IrradianceDirect = IrradianceDiffuse = LAI = Elev = Zenith = LeafAngleFactor = KbVal = KdVal = 0.0f;
}

CRadTrans::~CRadTrans() {}

 void CRadTrans::SetVal(CSolar Irradiance, float SLAI, float LeafAngleFactorIn)
 {
	 IrradianceDirect = Irradiance.GetPFDDirect();		//Z PFD (umol m-2 s) direct from solar irradiance
	 IrradianceDiffuse = Irradiance.GetPFDDiffuse();	//Z PFD (umol m-2 s) diffuse from solar irradiance
// the transmittance values obtained from Day and Bailey (1999), chap 3, 
// ecosystems of the world 20: greenhouse ecosystem, page 76
     LAI = SLAI;
     Elev =  Irradiance.GetSolarElevation();			//Z solar elevation angle in rad
	 Zenith = fmin(abs(PI / 2.0f - Elev), 1.56f);
     LeafAngleFactor = LeafAngleFactorIn;				//Z input leaf angle factor
     IsLeafAngleFactorUsed = true;
	 Kb(Zenith);										//Z KbVal: light interception by canopies, beam radiation extinction coefficient 
	 Kd(LAI);											//Z KdVal: extinction coefficient for black leaf in diffusive radiation
}

float CRadTrans::Reflect()
 {
	//Z for private variable, no need to use get function, can use this->
	// return (1-sqrt(absorp))/(1+sqrt(absorp))*(2*GetKb()/(GetKb()+GetKd()));
	return (1.0f - sqrt(absorp)) / (1.0f + sqrt(absorp)) * (2.0f * KbVal / (KbVal + KdVal));
 }


//Z Kb (a ratio for radiation)
//  Campbell, p 251, Eq. 15.4, Ratio of projected area to hemi-surface area for an ellisoid
//  x is a leaf angle distribution parameter
void CRadTrans::Kb(float theta) 
{
	float x, tmp;
	tmp = 0.5f;
	if (IsLeafAngleFactorUsed == true) { x = LeafAngleFactor; }
	else
	{
		switch (LeafAngle)
		{
		case Spherical:
			x = 1.0f;
			break;
		case Horizontal:
			x = 10.0f;
			break;
		case Vertical:
			x = 0.0f;
			break;
		case Corn:
			x = 1.37f;
			break;
		default:
			x = 1.0f;
		}
	}
	tmp = sqrt(sqr(x) + sqr(tan(theta))) / (x + 1.774f * pow(x + 1.182f, -0.733f));
	KbVal = tmp * clump;
}

//Z Kd, extinction coefficient 
//  for black leaf in diffusive radiation
void CRadTrans::Kd(float LA)
{
	const float gauss3[3] = { -0.774597f,0.000000f,0.774597f };
	const float weight3[3] = { 0.555556f,0.888889f,0.555556f };

	float K, FDiffuse, tmp, angle, x;
	if (IsLeafAngleFactorUsed == true) { x = LeafAngleFactor; }
	else
	{
		switch (LeafAngle)
		{
		case Spherical:
			x = 1.0f;
			break;
		case Horizontal:
			x = 10.0f;
			break;
		case Vertical:
			x = 0.0f;
			break;
		case Corn:
			x = 1.37f;
			break;
		default:
			x = 1.0f;
		}
	}
	
	FDiffuse = 0.0f;
	for (int ii = 0; ii < 3; ii++)  //diffused light ratio to ambient, itegrated over all incident angles from -90 to 90
	{
		angle = (PI / 2.0f) / 2.0f * (gauss3[ii]) + (PI / 2.0f) / 2.0f;
		tmp = sqrt(sqr(x) + sqr(tan(angle))) / (x + 1.774f * pow(x + 1.182f, -0.733f));
		FDiffuse = FDiffuse + (PI / 2.0f) / 2.0f * (2.0f * exp(-tmp * LA) * sin(angle) * cos(angle)) * weight3[ii];
	}
	if (LA <= 0.0f) { K = 0.0f; }
	else { K = -log(FDiffuse) / LA; }
	KdVal = K * clump;
  }

//Z: total irradiance at the top of the canopy, passed over from either observed PAR or TSolar or TIrradiance
float CRadTrans::Irradiancetot() 
{
	//IrradianceDirect: beam radiation at top of canopy, IrradianceDiffuse: diffuse radiation at top.
	// 
	//Z To be exact:
	//Z IrradianceDirect = Irradiance.GetPFDDirect()	 PFD (umol m-2 s) direct from solar irradiance
	//Z IrradianceDiffuse = Irradiance.GetPFDDiffuse() 	 PFD (umol m-2 s) diffuse from solar irradiance
	
	return (IrradianceDirect + IrradianceDiffuse);  
}

//Z: total irradiance (Direct and/or Diffuse) at depth L, simple Beer–Lambert law with beam/diffusive coefficients
float CRadTrans::Qtot(float L) 
{
	//return Irradiancetot() * exp(-sqrt(absorp) * (KdVal + KbVal) / 2.0 * L);

	//Z separate the direct and diffusive radiation and then sum together
	//  refer to Campbell and Norman, 1998 P258 Eqs 15.15-17
	return IrradianceDirect * exp(-sqrt(absorp) * KbVal * L) + IrradianceDiffuse * exp(-sqrt(absorp) * KdVal * L);
}
float CRadTrans::Qbt(float L) // total beam radiation at depth L
{
	return IrradianceDirect * exp(-sqrt(absorp) * KbVal * L);
}
float CRadTrans::Qb(float L) // unintercepted beam (direct beam) flux at depth of L within canopy
{
	return IrradianceDirect * exp(-KbVal * L);
}
float CRadTrans::Qd(float L) // net diffuse flux at depth of L within canopy
{
	return IrradianceDiffuse * exp(-sqrt(absorp) * KdVal * L);
}

//Z weighted average absorved diffuse flux over depth of L within canopy accounting for exponential decay
//  the return value is the averaged photon flux density on an "one-LAI" representative surface

//Z For example, suppose your LAI = 10, and you think that it is a 10-layer structure (each layer has LAI=1). 
//  Choose 1 layer (LAI=1) to represent that 10-layer canopy, 
//  then those equations are calculating the flux density received on that one representative layer with LAI=1.
//  so there is no difference between leaf and ground area numerically.

//Z the computed value can be used as leaf-area based computations, by multiplying correct LAI fractions
//  Campbel and Normal called it "leaf-hemi-surface area". 
float CRadTrans::Qdm()
{
	if (LAI <= 0.0f) { return 0.0f; }
	else {
		//  Integral Qd / Integral L
		//Z see Campbell and Norman section 15.9 P261 for an example
		//  but based on the "(1-exp)/exp" form, it looks like a result from integration
		return IrradianceDiffuse * (1.0f - exp(-sqrt(absorp) * KdVal * LAI)) / (sqrt(absorp) * KdVal * LAI);
	}
 }

//Z mean flux density on sunlit leaves: 
//  two components, one is direct PFD with extinction factor, one is diffusive PFD
//  because "sunlit", use GetKb() or KbVal is sufficient

//Z add absorp=0.85 as the leaf absorped photon flux density
//  Qsh() should already have absorp=0.85 in its own function
float CRadTrans::Qsl()
{
	return absorp * KbVal * IrradianceDirect + Qsh();
}

//Z add absorp=0.85 as the leaf absorped photon flux density
//  Qsh(L) should already have absorp=0.85 in its own function
float CRadTrans::Qsl(float L) // flux density on sunlit leaves at depth L (depth is marked via LAI incremental)
{
	return absorp * KbVal * IrradianceDirect + Qsh(L);
}

//Z mean flux density on shaded leaves over LAI, 
//  just one diffusive PFD component (but has three parts)

//Z add absorp=0.85 as the leaf absorped photon flux density
float CRadTrans::Qsh()
{
	return absorp * (Qdm() + Qsc() + Qsoilm());  // include soil reflection
}

float CRadTrans::Qsh(float L) // diffuse flux density on shaded leaves at depth L
{
	return absorp * (Qd(L) + Qsc(L) + Qsoilm());   // include soil reflection
}

// weighted average of Soil reflectance over canopy accounting for exponential decay

//Z the computed value can be used as leaf-area based computations, by multiplying correct LAI fractions
//  Campbel and Normal called it "leaf-hemi-surface area". 
//  That means it is the average PFD value on a representative, LAI=1, intersection surface
//  similar to Qdm()
float CRadTrans::Qsoilm() 
{
	if (LAI <= 0.0f) { return 0.0f; }
	else {
		//  Integral Qd / Integral L
		//Z to understand this equation
		//   "Qsoil()" is the radiation reach to soil surface
		//   "Qsoil() * rho_soil" is the radiation reflected by soil surface
		//   "whole equation" is the absorption of the reflected equation took by the canopy when the radiation is transporting upwards
		return Qsoil() * rho_soil * (1.0f - exp(-sqrt(absorp) * KdVal * LAI)) / (sqrt(absorp) * KdVal * LAI); 
	}
 }

// weighted average scattered radiation within canopy
float CRadTrans::Qsc() 
{
	float totBeam, nonscatt;
	if (LAI == 0.0f) { return 0.0f; }
	else
	{
		//total beam including scattered absorbed by canopy
		totBeam = IrradianceDirect * (1.0f - exp(-sqrt(absorp) * KbVal * LAI)) / (sqrt(absorp) * KbVal);
		//non scattered beam absorbed by canopy
		nonscatt = IrradianceDirect * (1.0f - exp(-KbVal * LAI)) / (KbVal);
		//mean scattered flux density
		return (totBeam - nonscatt) / LAI; 
	}
}

// scattered radiation at depth L in the canopy
float CRadTrans::Qsc(float L)
{
	return Qbt(L) - Qb(L); // total beam - nonscattered beam at depth L
}

// total PFD at the soil sufrace under the canopy
float CRadTrans::Qsoil() 
{
	return Qtot(LAI);
}

// sunlit LAI assuming closed canopy; thus not accurate for row or isolated canopy
float CRadTrans::LAIsl()
{
	if (Elev <= 0.01f) { return 0.0f; }
	else { return (1.0f - exp(-KbVal * LAI)) / KbVal; }
}

// shaded LAI assuming closed canopy
float CRadTrans::LAIsh()
{
    return LAI - LAIsl();
}

// sunlit fraction of current layer
float CRadTrans::Fsl(float L) 
{
	if (Elev <= 0.01f) { return 0.0f; }
	else { return exp(-KbVal * L); }
}

// shaded fraction of current layer
float CRadTrans::Fsh(float L)
{
	return 1.0f - Fsl(L);
}

//Z this is for testing, total PAR absorption (umol m^-2 ground sec^-1), but could provide useful information
//  total radiation absorbed by the plant canopy
float CRadTrans::get_RadAbsorbTot()
{
	if (LAI <= 0.0f) { return 0.0f; }
	else { return (Qsl() * LAIsl() + Qsh() * LAIsh()); }
}