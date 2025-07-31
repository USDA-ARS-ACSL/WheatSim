#include "solar.h"
class CRadTrans
{
private:

	enum CLeafAngle {Spherical, Horizontal, Vertical, Diaheliotropic, Empirical, Ellipsoidal, Corn};
	//type TCover = (Glass, Acrylic, polyethyl, doublepoly, whitewashed, NoCover);
	float absorp;			//Z leaf absorptivity for PAR
	float clump;
	float rho_soil;			//Z soil reflectivity for PAR band
	float IrradianceDirect, IrradianceDiffuse, LAI, Elev, Zenith, LeafAngleFactor, KbVal, KdVal;
	CLeafAngle LeafAngle;	//Z Solar elevation of a beam and cumulative LAI at the layer, diffused fraction (fdf)
	bool IsLeafAngleFactorUsed;

	float Qbt(float L);		//Z total beam radiation at depth L
	float Qb(float L);		//Z unintercepted beam (direct beam) flux at depth of L within canopy
	float Qd(float L);		//Z net diffuse flux at depth of L within canopy
	float Qsoil();			//Z total PFD at the soil surface under the canopy
	void Kb(float theta);	//Z Ratio of projected area to hemi-surface area for an ellisoid 
	void Kd(float LA);		//Z diffused light ratio to ambient, integrated over all incident angles from - 90 to 90

public:
	CRadTrans(void);		//sdf: diffused fraction of solar radiation
	~CRadTrans(void);

	void SetVal(CSolar Irradiance, float LAI, float leafAngleFactor);
	float Qsc();			//Z weighted average scattered radiation within canopy
	float Qsc(float L);		//Z scattered radiation at depth L in the canopy
	float Qtot(float L);	//Z total irradiance (dir + dif) at depth L, simple empirical approach
	float Irradiancetot();	//Z total PAR at top of the canopy
	float Qsl();			//Z mean flux density on sunlit leaves
	float Qsh();			//Z mean flux density on shaded leaves over LAI
	float Qsl(float L);		//Z flux density on sunlit leaves at delpth L
	float Qsh(float L);		//Z diffuse flux density on shaded leaves at depth L
	float Qdm() ;			//Z weighted average absorved diffuse flux over depth of L within canopy accounting for exponential decay
	float Qsoilm();			//Z weighted average of Soil reflectance over canopy accounting for exponential decay
	float Reflect();
	float LAIsl();			//Z sunlit LAI assuming closed canopy; thus not accurate for row or isolated canopy
	float LAIsh();			//Z shaded LAI assuming closed canopy
	float Fsl(float L);	//Z sunlit fraction of current layer
	float Fsh(float L);	//Z shaded fraction of current layer, "1.0 - sunlit"


	float GetKb() { return KbVal; }    // extiction coefficient assuming spherical leaf dist
	float GetKd() { return KdVal; }

	float get_RadAbsorbTot();	//Z this is for testing, total PAR absorption (umol m^-2 ground sec^-1), but could provide useful information
};
