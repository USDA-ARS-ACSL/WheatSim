#pragma once
#include "gas_ex_species_param.h"
#include "initinfo.h"

/*! \class CGasExchange 
* \brief Class for gas exchange calculations\n
* This class simulates gas exchange in plant leaves for C3 and C4 plants. \n
* \par Usage
	- Use <b>SetParams</b> to initialize the GasExchange object with parameters for a specific variety of plant
	- Use <b>SetVal</b> to pass environmental variables for a simulation and return a structure with output. \n
* See \ref Interface for details on input and output	
*/

 
class CGasExchange  
{
public:
	CGasExchange(const TGasExSpeciesParam& photoparam);
	~CGasExchange(void);

	//Z sets input values for calculations for a particular set of environmental variables
	void SetVal(float PhotoFluxDensity, const TInitInfo info, float Tair, float CO2,
		float RH, float wind, float Press, float width, float psil, float lfNContent, float BiomassLeafOver);
	
	//Z These functions return results of calculations 
	//  both variables and io functions coexist, it is not necessary, 
	//  e.g., if you have "get_ANet()", the "A_net" should be private

	float get_ANet() { return A_net; }						//Z return net photosynthesis (umol CO2 m-2 s-1) 
	float get_AGross() { return A_gross; }					//Z return gross photosynthesis  (umol CO2 m-2 s-1)
	float get_Transpiration() { return ET; }				//Z return transpiration rate (umol H2O m-2 s-1)
	float get_LeafTemperature() { return Tleaf; }			//Z return leaf temperature (C)
	float get_Ci() { return Ci; }							//Z return internal CO2 concentration (umol mol-1)
	float get_StomatalConductance() { return gs; }			//Z return stomatal conductance to water vapor (mol m-2 s-1)
	float get_BoundaryLayerConductance() { return gb; }		//Z return boundary layer conductance (mol m-2 s-1)
	float get_VPD() { return VPD; }							//Z return vapor pressure deficit (kpa)
	float get_leafpEffect() { return this->psileaf_stress; }

	float A_gross;											//Z gross photosynthesis (umol CO2 m-2 s-1) 
	float A_net;											//Z net photosynthesis (umol CO2 m-2 s-1) 
	float ET;												//Z transpiration rate (umol H2O m-2 s-1) 
	float Tleaf;											//Z leaf temperature (C) 
	float Ci;												//Z internal CO2 concentration (umol mol-1) 
	float gs;												//Z stomatal conductance to water vapor (mol m-2 s-1) 
	float gb;												//Z boundary layer conductance (mol m-2 s-1) 
	float Rdc;												//Z Plant dark respiration umol m-2 s-1. 
	float VPD;												//Z vapor pressure deficit (kpa) 
	float temp;

private:
	// variables passed as arguments to the constructor
	TGasExSpeciesParam sParms;
	enum gsModel { BBW, L, H }; //Z the model for stomata conductance

	float PhotoFluxDensity;		// Photosynthetic Flux Density (umol photons m-2 s-1) 
	float R_abs;				// Absorbed incident radiation (watts m-2)        
	float Tair;					// Air temperature at 2m, (C) 
	float CO2;					// CO2 concentration (umol mol-1 air) 
	float RH;					// Relative Humidity (%, i.e., 80) 
	float wind;					// Windspeed at 2 meters (m s-1) 
	float width;				// Leaf width (cm), when using leaf width, can find width/100.0, converting cm to m 
	float Press;				// Air pressure (kPa) 
	float Theta;				// Initial slope of CO2 response(umol m2 s - 1) - de Pury(1997)
	float psileaf;				// leaf water potential, MPa
	float psileaf_stress;		// 0 to 1 factor for stomatal closure effect of leaf water potential, i.e., effect of leaf water pressure (Mpa) on stomatal conductance.
	float lfNContent;			//YY lfNContent is leaf nitrogen content in unit of g m-2(leaf)
	float leaf_age;
	float BiomassLeafOver;

	float Ci_Ca;				//Z Ratio of internal to external CO2, unitless
	float errTolerance;			//Z error tolerance for iterations
	float eqlTolerance;			//Z equality tolerance
	int    iter_total;			//Z holds total number of iterations
	int    iter1, iter2;		//Z holds iteration counters
	int    iter_Ci;				//Z iteration value for Ci umol mol-1, internal CO2 concentration
	bool   isCiConverged;		//Z true if Ci (internal CO2 concentration) iterations have converged

	// Main module to calculate gas exchange rates
	// produce photosynthesis and leaf temperature
	// Incorporate water stress effect
	void GasEx_psil(float psileaf, const TInitInfo info);
	//Z C3 photosynthesis calculation
	void Photosynthesis(float Ci, float psileaf, const TInitInfo info);  
	//Z calculates leaf temperature and transpiration
	void EnergyBalance();     

	//Z called iterively to find optimal internal CO2 concentration returns optimal internal CO2 concentration (CO2i) 
	//  secant search to find optimal internal CO2 concentration
	float SearchCi(float CO2i, float psileaf, const TInitInfo info); 

	//Z Calls photosynthesis modules to evaluate Ci dependent calculations returns difference between old Ci and new Ci
	//  Calculates a new value for Ci for the current values of photosynthesis and stomatal conductance
	//  used in "SearchCi" to fulfull the Ci iteration
	float EvalCi(float Ci, float psileaf, const TInitInfo info);   

	//Z Stomatal conductance (mol m-2 s-1)
	//  stomatal conductance for water vapor in mol m-2 s-1
	//  need a parameter to convert to CO2
	float gsw(float pressure, const TInitInfo info);   
	
	//Z Conductance for turbulant vapor transfer in air - forced convection (mol m-2 s-1)
	//  boundary layer conductance to vapor
	float gbw();           

	//Z Saturated vapor pressure at given temperature. kPa
	float Es(float Temperature);      
	//Z Slope of the vapor pressure curve, first order derivative of Es with respect to T
	//  sometimes present as "gamma" 
	float Slope(float Temperature);  

	//Z mitochondrial respiration
	float Rd();

	//Z Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
	float set_PSIleafeffect(float pressure, const TInitInfo info);

	float QuadSolnUpper(float a, float b, float c); //Z Upper part of quadratic equation solution
	float QuadSolnLower(float a, float b, float c); //Z Lower part of quadratic equation solution
};


 
