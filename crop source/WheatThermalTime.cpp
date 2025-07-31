#include "WheatThermalTime.h"

using namespace std;

WheatThermalTime::WheatThermalTime(void)
{
	//Z default initialization for dssat thermal time calculator
	airTmprCur = 25.0f;
	soilTmprCur = 25.0f;
	elongation = false;

	Tbase = 0.0f;
	Topt = 26.0f;
	Tmax = 34.0f;

	T_diff_bo = Topt - Tbase;
	T_diff_mo = Tmax - Topt;

	TDmax = 0.0f;
	TDmin = 0.0f;
	DayL = 0.0f;
	TTdSum_plt_dssat = 0.0f;
	TTdSum_env_dssat = 0.0f;
	dTmpr = 0.0f;
	timeStep = 60.0f;
	VernalTotal = 0.0f;
	fD = fV = 0.0f;
	TTd_Cur = 0.0f;
	decreasingBeta = Topt / (Tmax - Topt);

	//Z default initialization for the what fspm model
	Temp_Tref = 12.0f;
	Temp_Ttransition = 9.0f;
	Temp_Ea_R = 8900.0f;
	Temp_DS_R = 68.432f;
	Temp_DH_R = 20735.5f;

	TimeTeqCur = 0.0f;
	TTdSum_fspm = 0.0f;
	ArrheniusTrans_fspm = 0.0f;
	modifiedArrheniusTref_fspm = 0.0f;

	//Z reset a 5-day moving average of air tmpr
	Tave = 0.0f;
	for (int ii = 0; ii < 120; ii++) { TRecord[ii] = 0.0f; }
	Tcurrentposition = 0;
	TRecordNumber = 0.0f;
}

WheatThermalTime::~WheatThermalTime(void) {}

void WheatThermalTime::initialize(float dt, float gddbase, float gddopt, float gddmax, 
	float T_ref, float T_trans, float T_EaR, float T_DsR, float T_DhR)
{
	airTmprCur = 0.0f;
	soilTmprCur = 0.0f;
	elongation = false;

	TDmax = 0.0f;
	TDmin = 0.0f;
	DayL = 0.0f;
	timeStep = dt;
	TTdSum_plt_dssat = 0.0f;
	TTdSum_env_dssat = 0.0f;

	//Z initialize the dssat thermal time calculator
	Tbase = gddbase;
	Topt = gddopt;
	Tmax = gddmax;
	decreasingBeta = Topt / (Tmax - Topt);

	//Z initialize the wheat fspm thermal time caluclator
	Temp_Tref = T_ref;
	Temp_Ttransition = T_trans;
	Temp_Ea_R = T_EaR;
	Temp_DS_R = T_DsR;
	Temp_DH_R = T_DhR;

	TimeTeqCur = 0.0f;
	TTdSum_fspm = 0.0f;
	ArrheniusTrans_fspm = ArriheniusFspm(Temp_Ttransition);
	modifiedArrheniusTref_fspm = modifyArriheniusFspm(Temp_Tref);

	//Z operate the thermal time calulator
	add_TTd_dssat(airTmprCur, soilTmprCur, elongation);
}

void WheatThermalTime::add_TTd_dssat(float AirTmpr, float SoilTmpr, bool Elong)
{
	airTmprCur = AirTmpr;
	soilTmprCur = SoilTmpr;
	elongation = Elong;

	//Z timestep in develop in minutes, convert it to days
	float dD = timeStep * DAYPERMINUTES;

	//Z TTd cumulated based on air tmpr, use for comparing with weather station data
	if (airTmprCur <= Tbase)
	{
		dTmpr = 0.0f;
	}
	else if (airTmprCur > Tbase && airTmprCur <= Topt)
	{
		dTmpr = airTmprCur - Tbase;
	}
	else if (airTmprCur > Topt && airTmprCur <= Tmax)
	{
		dTmpr = airTmprCur - Tbase;
	}
	else
	{
		dTmpr = 0.0f;
	}
	TTdSum_env_dssat += dTmpr * dD;

	//Z TTd cumulated based on air/soil tmpr, observed by plant
	float Tcur = airTmprCur;
	if (!elongation) Tcur = std::max((0.75f * airTmprCur + 0.25f * soilTmprCur), 0.0f);

	if (Tcur <= Tbase)
	{
		dTmpr = 2.0f;
	}
	else if (Tcur > Tbase && Tcur <= Topt)
	{
		dTmpr = std::max(Tcur - Tbase, 2.0f);
	}
	else if (Tcur > Topt && Tcur <= Tmax)
	{
		dTmpr = decreasingBeta * (Tmax - Tcur);
	}
	else
	{
		dTmpr = 0.0f;
	}

	TTd_Cur = (dTmpr * dD) * std::min(fD, fV);
	TTdSum_plt_dssat += TTd_Cur;
}

void WheatThermalTime::FactorPhotoperiod(float Daylength)
{
	DayL = Daylength;
	fD = 1.0f - 0.002f * 1.5f * (20.0f - DayL) * (20.0f - DayL);
	fD = std::max(fD, 0.0f);
	fD = std::min(fD, 1.0f);
}

void WheatThermalTime::FactorVernalisation(float airTmax, float airTmin)
{
	TDmax = airTmax;
	TDmin = airTmin;
	if (TDmax <= 20.0f && TDmin <= 10.0f)
	{
		float DeltaT = TDmax - TDmin;
		float aa = 1.4f - 0.0778f * airTmprCur;
		float bb = 0.5f + 13.44f * airTmprCur / (DeltaT + 3.0f) / (DeltaT + 3.0f);
		VernalTotal += (std::min(aa, bb) / 24.0f);
	}

	if (TDmax > 20.0f && VernalTotal <= 10.0f)
	{
		float VT = VernalTotal;
		float aa = 0.5f * (TDmax - 20.0f);
		VernalTotal -= (std::min(aa, VT) / 24.0f);
	}

	fV = 1.0f - (0.0054545f * 2.0f + 0.0003f) * (50.0f - VernalTotal);
	fV = std::max(fV, 0.0f);
	fV = std::min(fV, 1.0f);

	//	cout << "VernalTotal: " << VernalTotal << " fV: " << fV << '\n';
}

//Z the Arrihenious function used in Fspm to adjust the tmpr response shape
//  tmpr unit in oC but internal unit is in K
float WheatThermalTime::ArriheniusFspm(float t)
{
	float tmpr = t + 273.15f;
	float arr = tmpr * expf(-Temp_Ea_R / tmpr) / (1.0f + expf(Temp_DS_R - Temp_DH_R / tmpr));
	return arr;
}

//Z the modified Arrihenious function used in Fspm to adjust the tmpr response shape
//  Return value of equation from Johnson and Lewin (1946) for temperature. The equation is modified to return zero below zero degree.
//  tmpr unit in oC but internal unit is in K
float WheatThermalTime::modifyArriheniusFspm(float tmpr)
{
	if (tmpr < 0.0f) { return 0.0f; }
	else if (tmpr < Temp_Ttransition) { return tmpr * ArrheniusTrans_fspm / Temp_Ttransition; }
	else { return ArriheniusFspm(tmpr); }
}

//Z thermal time calculator for wheatfpsm
void WheatThermalTime::add_TTd_fspm(float AirTmpr, float SoilTmpr, bool Elong)
{
	airTmprCur = AirTmpr;
	soilTmprCur = SoilTmpr;
	elongation = Elong;

	//Z timestep in develop in minutes, convert it to days
	float dD = timeStep * DAYPERMINUTES;
	//Z TTd cumulated based on air/soil tmpr, observed by plant
	float Tcur = airTmprCur;
	if (!elongation) { Tcur = std::max((0.75f * airTmprCur + 0.25f * soilTmprCur), 0.0f); }

	//Z tmperature equivalent time step 
	TimeTeqCur = dD * modifyArriheniusFspm(Tcur) / modifiedArrheniusTref_fspm;
	//std::cout << TimeTeqCur << '\n';

	//Z cumulate the thermal time
	if (Tcur < 0.0f) { TTdSum_fspm += 0.0f; }
	else { TTdSum_fspm += (TimeTeqCur * Temp_Tref); }
}


void WheatThermalTime::TmprMovingAve(float airTmpr)
{
	//Z do a 5-day moving average of air tmpr
	if (TRecordNumber <= 0.0)
	{
		TRecord[0] = airTmpr;
		Tave = airTmpr;
		Tcurrentposition = 1;
		TRecordNumber = 1.0f;
	}
	else
	{
		float tmpr_kick_out = TRecord[Tcurrentposition];
		float T_sum = Tave * TRecordNumber;
		T_sum = T_sum - tmpr_kick_out + airTmpr;
		TRecord[Tcurrentposition] = airTmpr;

		TRecordNumber = std::min(TRecordNumber + 1.0f, 120.0f);
		Tave = T_sum / TRecordNumber;

		Tcurrentposition = Tcurrentposition + 1;
		if (Tcurrentposition == 120)
		{
			Tcurrentposition = 0;
		}
	}
}

void WheatThermalTime::update(float airTmpr, float airTmax, float airTmin, float Daylength, float soilTmpr, bool elongation)
{
	FactorPhotoperiod(Daylength);
	FactorVernalisation(airTmax, airTmin);
	add_TTd_dssat(airTmpr, soilTmpr, elongation);
	add_TTd_fspm(airTmpr, soilTmpr, elongation);
	TmprMovingAve(airTmpr);
}


