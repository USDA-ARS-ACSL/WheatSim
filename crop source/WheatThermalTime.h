#pragma once
#ifndef _WHEAT_THERMALTIME_H_
#define _WHEAT_THERMALTIME_H_

#include <iostream>
#include <algorithm>
#include "stdafx.h"

#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.00069444444f

/*Z two ways to compute thermal time
1. Dssat/apsim thermal time: 
	will be used in thermal time driven organ growth, 
    unit in thermal time gdd
2. WheatFspm: 
	refer to equivalent time under reference temperature, 
	used in the BETA function governed growth, 
	unit in ordinary time scale or thermal time scale,
*/

class WheatThermalTime
{
public:
	WheatThermalTime(void);
	~WheatThermalTime(void);
	void initialize(float timestep, float gddbase, float gddopt, float gddmax, float T_ref, float T_trans, float T_EaR, float T_DsR, float T_DhR);
	//Z: IO for thermal time
	//   dssat based thermal time calculator
	float get_TmprCur() { return airTmprCur; }						//Z get current temperature
	float get_TTdSum_plt_dssat() { return TTdSum_plt_dssat; }		//Z get cumulated thermal time "observed by plant"
	float get_TTdSum_env_dssat() { return TTdSum_env_dssat; }		//Z get cumulated thermal time purely based on weather
	float get_TTdCur_dssat() { return TTd_Cur; }					//Z get current thermal time for the current hour
	//    wheat fspm based thermal time calculator
	float get_TimeTeqCur_fspm() { return TimeTeqCur; }				//Z get time equivalent to a reference temperature i.e. temperature-compensated time (Parent, 2010).
	float get_TTdSum_fspm() { return TTdSum_fspm; }					//Z get cumulated fspm thermal time, adjust based on ref tmpr and arrhenius function

	//    airtmpr moving average
	float get_movingAveTmpr() { return Tave; }

	//Z: gdd computation
	void add_TTd_dssat(float AirTmpr, float SoilTmpr, bool Elong);
	void add_TTd_fspm(float AirTmpr, float SoilTmpr, bool Elong);
	void update(float airTmpr, float airTmax, float airTmin, float Daylength, float soilTmpr, bool elongation);
	void FactorPhotoperiod(float Daylength);
	void FactorVernalisation(float airTmax, float airTmin);

	float ArriheniusFspm(float tmpr);
	float modifyArriheniusFspm(float tmpr);

	void TmprMovingAve(float airTmpr);

private:

	// ************
	//Z first for dssat gdd computation

	float airTmprCur;		//Z current temperature
	float soilTmprCur;
	bool  elongation;

	float Tbase;			//Z based temperature for gdd computation
	float Topt;				//Z optimal temperature for gdd computation, plant growth
	float Tmax;				//Z max temperature for gdd computation, beyond that, plant growth stops and gdd=0
	float T_diff_bo;		//Z temperature difference, opt-base
	float T_diff_mo;		//Z temperature difference, max-opt

	float TTd_Cur;			//Z current gdd at this hour period

	float dTmpr;			//Z "apparent tmpr" based on tmpr, opt tmpr, ...,
	float timeStep;			//Z timestep, should be 1/24 day
	float decreasingBeta;	//Z Topt<Tcur<Tmax, the slope for gdd decreasing

	//Z sum of the gdd based on our "peak function" for thermal time
	float TTdSum_plt_dssat;		//Z the sum "observed" by the plants, including photoperiod and vernalisation factor
	float TTdSum_env_dssat;		//Z the gdd based on weather

	//Z add daily max and min temperature for Vernalisation
	//Z add daylength for photoperiod
	float TDmax, TDmin, DayL;
	//Z add photoperiod and vernalisation factors
	float fD, fV, VernalTotal;

	// ************
	//Z second for fspm wheat model
	//i.e., the temerpature response function
	float Temp_Tref;		//Z =12 Arbitrary reference temperature (°C)
	float Temp_Ttransition;	//Z =9 Below this temperature f = linear function of temperature instead of Arrhenius - like(°C)
	// parameters for the Arrhenius and modified Arrhenius models
	float Temp_Ea_R;		//Z =8900 Parameter Ea / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (K)
	float Temp_DS_R;		//Z =68.432 Parameter deltaS / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (dimensionless)
	float Temp_DH_R;		//Z =20735.5 Parameter deltaH / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (K)
	
	float TimeTeqCur;			//Z time equivalent to a reference temperature i.e. temperature-compensated time (Parent, 2010).
	float TTdSum_fspm;			//Z get cumulated fspm thermal time, adjust based on ref tmpr and arrhenius function
	float ArrheniusTrans_fspm;			//Z auxillary function that compute the arrhenius for the transition tmpr
	float modifiedArrheniusTref_fspm;	//Z auxillary function that compute the modified arrhenius for the reference tmpr
	
	// ************
	//Z thrid a 5-day moving average of airtmpr
	float Tave;
	float TRecord[120];
	int Tcurrentposition;
	float TRecordNumber;

};
#endif