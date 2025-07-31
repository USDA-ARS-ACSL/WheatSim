#include "stdafx.h"
#include <cmath>
#include "WheatDevelopment.h"
#include <iostream>
#include <string>
#include <algorithm>
//#using <mscorlib.dll>
#define PHYLLOCHRON 106.0f
#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.00069444444f
#define DAYPERHOUR 0.041666666667f

using namespace std;

WheatDevelopment::WheatDevelopment(const TInitInfo& info)
{
	if (info.latitude >= 0) { NorthHemisphere = true; }
	else { NorthHemisphere = false; }

	//Z GDD numbers and thermal time
	Cur_TTd = 0.0f;
	TimeTeqCur = 0.0f;
	//Z time step
	TimeStep = info.timeStep;
	DayTime = 0.0f;
	DaysAfterAntss = 0.0f;
	WheatTTd.initialize(static_cast<float>(TimeStep), info.gdd_base, info.gdd_opt, info.gdd_max, 
		info.Temp_Tref, info. Temp_Ttransition, info.Temp_Ea_R, info.Temp_DS_R, info.Temp_DH_R);

	// initialization
	//****** Germination & Emergence Group **************

	temcon = 0;
	condit = -1;			//Z -1 marks before the simulation
	oldcon = -1;			//Z -1 marks before the simulation
	doy = -1;
	doyRecord = -1;			//Z just need to make sure this marks the starting of the code, should be some nonsense number here

	SeedGerminationRate = 0.0f;
	GerminFrac = 0.0f, FracGermin_Real = 0.0f;
	EmergFrac = 0.0f, FracEmerg_Real = 0.0f;
	SeedNum = static_cast<int>(info.plantDensity), GerminateSeedNum = 0, PlantNum = 0;
	TTd_seed = 0.0f, Cur_Germin_TTd = 0.0f, Cur_Emerg_TTd = 0.0f;

	//Z plant ratio
	//  wheat (similar to cereal rye) is different from maize, that the plant number will change over time
	//  possible reason: 1. germination/emergence have a certian fraction
	//                   2. plant may die in winter or due to ambient stress
	//Z "PlantLivingFraction" = plant number / initialized number
	//  can be used to adjust any variables like plant density, popslab ......
	PlantLivingFraction = 1.0f;

	//emercv : Coefficient of variation for emergence curves.
	emercv = 0.20f;
	//mgerm(1:3) : Mean degree-days needed for germination for each of the three soil water condition values.
	mgerm[0] = 80.0f, mgerm[1] = 90.0f, mgerm[2] = 110.0f;
	//melrat(1:3) : Mean elongation rate (for seedlings) for each of the three soil water condition values (growing degree-days/cm).
	melrat[0] = 20.0f, melrat[1] = 25.0f, melrat[2] = 30.0f;
	//Z seedDepth (cm) in the "info" structure, emergence gdd depends on seed depth
	for (int ii = 0; ii < 3; ii++)
	{
		memerg[ii] = melrat[ii] * info.seedDepth + mgerm[ii];
		sgerm[ii] = emercv * mgerm[ii];
		semerg[ii] = emercv * memerg[ii];
		upcutg[ii] = mgerm[ii] + 3.0f * sgerm[ii];
		locutg[ii] = mgerm[ii] - 3.0f * sgerm[ii];
		upcute[ii] = memerg[ii] + 3.0f * semerg[ii];
		locute[ii] = memerg[ii] - 3.0f * semerg[ii];
	}

	//****** Growth Stage Till Jointing **************
	//Z gdd from doy==1 to till Jointing
	TTd_joint = 0.0f;
	TTd_plant = 0.0f;
	TTd_2_singleRidge = PHYLLOCHRON * info.srpa;
	TTd_since_singleRidge = 0.0f;
	TTd_since_jointing = 0.0f;
	TTd_2_doubleRidge = PHYLLOCHRON * (info.srpa + info.drpa);
	TTd_2_flowerPrimInit = PHYLLOCHRON * (info.srpa + info.drpa + info.fpibpa);
	TTd_2_terSpikeGrow = PHYLLOCHRON * (info.srpa + info.drpa + info.tspa);
	TTd_2_elongation = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa);
	TTd_2_jointStart = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa);
	TTd_2_bootStart = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa);
	TTd_2_headStart = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + info.headpa);
	TTd_2_anthesisStart = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + info.headpa + info.antspa);
	TTd_2_anthesisEnd= PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + info.headpa + info.antspa+info.antepa);
	TTd_2_maturity = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + info.headpa + info.antspa) + info.matpa;
	

	TTd_FlagLf_min = 0.0f;
	TTd_FlagLf_min_aux = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + 0.5f);

	//Z chaff growth
	TTd_2_chaffStart = PHYLLOCHRON * (info.srpa + info.drpa + 0.5f);
	TTd_2_chaffEnd = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + info.headpa + info.antspa) + info.matpa * 0.5f;
	DurationChaffGrowth = PHYLLOCHRON * (info.jtpa - 0.5f + info.bootpa + info.headpa + info.antspa) + info.matpa * 0.5f;

	//Z rachis growth
	DurationRachisGrowth = PHYLLOCHRON * info.jtpa;
	
	//****** Plant Death Factors **************
	//rye death should enable a cumulation computation
	LowKill_Frac = 0.0f, HighKill_Frac = 0.0f, PresentKill_Frac = 0.0f, TotalKill_Frac = 0.0f;

	T_cur = 0.0f;
	T_air = 0.0f;
	T_ave = 0.0f;

	Cold_Time = 0.0f;
	Cold_Time_Ratio_Plant = 1.0f;
	Cold_Time_Ratio_Joint = 1.0f;

	//Z this predawn leaf water potential (bar) is stored but not processed in this class, save and offer it to leaf class
	PredawnLWP = -0.05f;

	//Z shaded effects, compute in plant class, assigned to development and stored, used in leaf class
	ShadeEffect = 1.0f;

	//Z 
	pltLivingBiomass = 0.0f;
	pltLivingNitrogenMass = 0.0f;
	pltTillerNum = 0.0f;

	//****** parameters for grain and kernel *******
	this->KernelFitCurv(0, fwtc);
	this->KernelFitCurv(0, durc);
}

WheatDevelopment::~WheatDevelopment() {}

void WheatDevelopment::WheatDelpUpdate(const TWeather& wthr)
{
	//Z get necessary climate data for this module.
	T_cur = std::max(0.0f, wthr.airT);
	T_air = wthr.airT;
	doy = wthr.doy;
	PredawnLWP = wthr.PredawnLWP;
	DayTime=wthr.daytime;

	//Z first update the thermal time computation for the whole plant
	// WheatTTd.update(T_cur, wthr.airTDmax, wthr.airTDmin, wthr.dayLength, wthr.soilT, elongationStart.done);
	WheatTTd.update(T_cur, wthr.airTDmax, wthr.airTDmin, wthr.dayLength, wthr.soilT, singleRidge.done);

	Cur_TTd = WheatTTd.get_TTdCur_dssat();
	TimeTeqCur = WheatTTd.get_TimeTeqCur_fspm();
	T_ave = WheatTTd.get_movingAveTmpr();

	//Z Germination and Emergence
	//  The criterion is if the final emergence fraction = 1.0, i.e., the whole emergence is complete
	if (FracEmerg_Real < 1.0f)
	{
		//Z water saturation, need 3 cut-off values
		if (wthr.wfps >= 0.35f) { temcon = 0; }			//Z wfps > 0.35 (averaged degree of saturation)
		else if (wthr.wfps >= 0.25f) { temcon = 1; }	//Z wfps = 0.25 - 0.35
		else if (wthr.wfps >= 0.15f) { temcon = 2; }	//Z wfps = 0.15 - 0.25
		else { temcon = 3; }							//Z wfps < 0.15

		//Z first call to set up a germination rate
		if (GerminateSeedNum == 0 && temcon != 3)
		{
			//Z make germination as a continous function w.r.t. wfps
			if (wthr.wfps >= 0.35f) { SeedGerminationRate = 1.0f; }
			else if (wthr.wfps >= 0.25f) { SeedGerminationRate = 0.9f + (wthr.wfps - 0.25f); }
			else if (wthr.wfps >= 0.15f) { SeedGerminationRate = 0.8f + (wthr.wfps - 0.15f); }
			else if (wthr.wfps >= 0.05f) { SeedGerminationRate = 0.7f + (wthr.wfps - 0.05f); }
			else { SeedGerminationRate = 0.7f; }

			GerminateSeedNum = static_cast<int>(SeedGerminationRate * static_cast<float>(SeedNum) + 0.5);
			PlantLivingFraction = static_cast<float>(GerminateSeedNum) / static_cast<float>(SeedNum);
		}

		//Z start the germination and emergence procedure
		/*
		Use today and previous soil conditions to determine whether development occurs today and the algorithm to calculate it.
		Four general situations must be handled, in shootgro, temcon = 1,2,3,4 and condit = 0, 
		while in c++, we use 0,1,2,3 and -1 (always minus 1)
			Condit = yesterday's soil condition index
			temcon = today's soil condition index

		(1) If soil conditions are too dry for any development (temcon = 3),
			and this is either the day of planting (condit = -1)
			or not the first day of the dry period (condit = 3),
			do nothing except reset condit to today's value of temcon in preparation for tomorrow.

		(2) If this is the first dry day after a period of conditions wet enough to allow growth (temcon = 3, condit = 0,1,2)
			kill any plants that have germinated but not yet emerged,
			then set oldcon to yesterday's soil condition value (condit) for comparison when conditions improve.

		(3) If soil conditions are adequate for seed activity (temcon = 0,1,2)
			and there has been a change in the soil condition value since the previous day on which development occurred,
			either with (condit = 3, oldcon != temcon) or without an intervening dry period (condit = 1,2, or 3)
			recalculate progress using a new growth curve, then perform the calculations associated with growth.

		(4) If soil conditions are adequate for seed activity (temcon = 0,1,2)
			and this is the first day of activity (oldcon = -1) either with (condit = 3) or without (condit = -1) an initial dry spell,
			or if soil conditions are the same as on the previous day on which any activity occurred,
			either with (oldcon = temcon) or without (condit = temcon) an intervening dry spell,
			perform the calculations associated with growth (there is NO need to update the growth curve,
			OTHERWISE will be no normal germination or emergence).
		*/

		//Z  If soil water conditions are not too dry for growth to occur today, case (3) and (4)
		//   no need to perform curve recalculation for (4)
		if (temcon != 3)
		{
			/*
			If soil conditions have changed since the previous day on which seed activity occurred
			(either with case (1) or without case (2) an intervening severely dry spell)
			recalculate the germination and combined emergence curves.

			Gfrac and efrac are the current cumulative fractions of plants that have germinated or emerged, respectively.
			The iterative function RECURV calculates the GDD values corresponding to gfrac and efrac on the curve associated with the current soil condition (temcon).

			The germination and emergence curves that result from changes in soil conditions will lie between the curves for condition 0 and 2.
			Can be shown by plotting gfrac and efrac  against gdds, with gdds on the same scale as bgddg and bgdde.
			(Gdds is the cumulative GDD since planting, not including days at condition 3.)
			The first parameter passed to RECURV indicates whether the germination curves (1) or the emergence curves (2) are to be used.
			*/

			//Z General case 3, intervening dry spell:
			if (condit != 3 && oldcon != temcon && oldcon != -1)
			{
				if (GerminFrac < 1.0f) { Cur_Germin_TTd = ReCurv(1, GerminFrac); }	//first argument = 1 germination
				if (EmergFrac < 1.0f) { Cur_Emerg_TTd = ReCurv(2, EmergFrac); }		//first argument = 2 emergence
			}
			//Z General case 3, no intervening dry spell :
			if (condit != temcon && condit != -1 && condit != 3)
			{
				if (GerminFrac < 1.0f) { Cur_Germin_TTd = ReCurv(1, GerminFrac); }	//first argument = 1 germination
				if (EmergFrac < 1.0f) { Cur_Emerg_TTd = ReCurv(2, EmergFrac); }		//first argument = 2 emergence
			}

			TTd_seed += Cur_TTd;
			Cur_Germin_TTd += Cur_TTd;
			Cur_Emerg_TTd += Cur_TTd;

			if (Cur_Germin_TTd > upcutg[temcon])
			{
				GerminFrac = 1.0f;
				FracGermin_Real = 1.0f;
			}
			else
			{
				GerminFrac = EmCurv(Cur_Germin_TTd, 1);
				FracGermin_Real = GerminFrac;
				if (Cur_Germin_TTd < locutg[temcon])
				{
					FracGermin_Real = 0.0f;
				}
			}

			if (Cur_Emerg_TTd > upcute[temcon])
			{
				EmergFrac = 1.0f;
				FracEmerg_Real = 1.0f;
			}
			else
			{
				EmergFrac = EmCurv(Cur_Emerg_TTd, 2);
				FracEmerg_Real = EmergFrac;
				if (Cur_Emerg_TTd < locute[temcon])
				{
					FracEmerg_Real = 0.0f;
				}
			}

			//Z set up the germination and emergence marker
			//Z first germination initiation/start, use a small number here, but same as the one for emergence
			//  germination initiation must be earlier than emergence
			if (!germinInit.done && FracGermin_Real >= 0.01f)
			{
				germinInit.daytime = DayTime;
				germinInit.done = true;

			}
			if (!germination.done && FracGermin_Real >= 0.5f)
			{
				germination.daytime = DayTime;
				germination.done = true;

			}
			if (!emergence.done && FracEmerg_Real >= 0.01f)
			{
				emergence.daytime = DayTime;
				emergence.done = true;
			}

			if (TTd_seed >= 350.0f)
			{
				// exceed the target emergence day, kill all the plants
				PlantDie(1.0f);
			}
		}
		// If it is the first day of the dry period, but not the day of planting execute the general case (2) :
		else if (condit != -1 && condit != 3)
		{
			//Z the ambient water condition is very pool
			//  all germinated plants want to die
			PlantDie(FracGermin_Real);
			if (doy > doyRecord)
			{
				//Z doyRecord mark day + 1, where the water condition should be updated
				oldcon = condit;
				doyRecord = doy;
			}
		}

		//Z daily update the previous day water condition
		if (doy > doyRecord) condit = temcon;

		//Z based on the emergence rate
		//  a certain number of plant will grow
		if (FracEmerg_Real < 0.01f)
		{
			PlantNum = 0;
			PlantLivingFraction = 0.001f;
		}
		else if (FracEmerg_Real < 1.0f)
		{
			PlantNum = static_cast<int>((FracEmerg_Real - TotalKill_Frac) * static_cast<float>(GerminateSeedNum) + 0.5f);
			PlantNum = std::max(PlantNum, 1);
			PlantLivingFraction = static_cast<float>(PlantNum) / static_cast<float>(SeedNum);
		}
		else
		{
			PlantNum = static_cast<int>((1.0f - TotalKill_Frac) * static_cast<float>(GerminateSeedNum) + 0.5f);
			PlantNum = std::max(PlantNum, 1);
			PlantLivingFraction = static_cast<float>(PlantNum) / static_cast<float>(SeedNum);
		}

	}

	//Z use the growth gdd (not joint gdd) to determine leaf emergence
	//Z use the gdd to joint to determine start of elongation
	if (emergence.done)
	{
		//Z gdd since emergence
		TTd_plant = TTd_plant + Cur_TTd;
		//Z in shootgro, TTd_joint start after vernalization, aka, it starts on Jan. 1 for the following year in north hemisphere
		//  or Jul 1 in south hemisphere
		if (NorthHemisphere) 
		{
			if (doy < 213)
			{
				TTd_joint = TTd_joint + Cur_TTd;
			}
		}
		else 
		{
			if (doy > 183 || doy < 30)
			{
				TTd_joint = TTd_joint + Cur_TTd;
			}
		}
		
		//Z judge if the single ridge stage is reached
		if (!singleRidge.done)
		{
			if (TTd_joint >= TTd_2_singleRidge)
			{
				singleRidge.daytime = DayTime;
				singleRidge.done = true;
			}
		}

		//Z judge if growth accerlation should occur under the condition that the single ridge is reached
		if (!acceleration.done)
		{
			if (TTd_joint >= TTd_2_singleRidge)
			{
				acceleration.daytime = DayTime;
				acceleration.done = true;
			}
		}

		//Z judge if the double ridge stage is reached
		if (!doubleRidge.done)
		{
			if (TTd_joint >= TTd_2_doubleRidge)
			{
				doubleRidge.daytime = DayTime;
				doubleRidge.done = true;
			}
		}

		//Z flower primordium initiation begins (0.333 PHYLLOCHRON after double ridge)
		if (!flowerPrimInit.done)
		{
			if (TTd_joint >= TTd_2_flowerPrimInit)
			{
				flowerPrimInit.daytime = DayTime;
				flowerPrimInit.done = true;
			}
		}

		//Z judge if the elongation starts
		if (!elongationStart.done)
		{
			if (TTd_joint >= TTd_2_elongation)
			{
				elongationStart.daytime = DayTime;
				elongationStart.done = true;
				TTd_FlagLf_min = TTd_plant - TTd_joint + TTd_FlagLf_min_aux;
			}
		}

		//Z judge terminal spikelet growth stages
		if (!terSpikeGrow.done)
		{
			if (TTd_joint >= TTd_2_terSpikeGrow)
			{
				terSpikeGrow.daytime = DayTime;
				terSpikeGrow.done = true;
			}
		}

		//Z judge jointing stages
		if (!jointStart.done)
		{
			if (TTd_joint >= TTd_2_jointStart)
			{
				jointStart.daytime = DayTime;
				jointStart.done = true;
			}
		}

		//Z judge heading stages
		if (!headStart.done)
		{
			if (TTd_joint >= TTd_2_headStart)
			{
				headStart.daytime = DayTime;
				headStart.done = true;
			}
		}

		//Z judge anthesis starting point
		if (!anthesisStart.done)
		{
			if (TTd_joint >= TTd_2_anthesisStart)
			{
				anthesisStart.daytime = DayTime;
				anthesisStart.done = true;
			}
			DaysAfterAntss += TimeStep * DAYPERMINUTES;
		}

		//Z judge anthesis ending point
		if (!anthesisEnd.done)
		{
			if (TTd_joint >= TTd_2_anthesisEnd)
			{
				anthesisEnd.daytime = DayTime;
				anthesisEnd.done = true;
			}
		}

		//Z judge maturity
		if (!maturity.done)
		{
			if (TTd_joint >= TTd_2_maturity)
			{
				maturity.daytime = DayTime;
				maturity.done = true;

				death.daytime = DayTime;
				death.done = true;
			}
		}

		//Z chaff growth stages
		if (!chaffStart.done)
		{
			if (TTd_joint >= TTd_2_chaffStart)
			{
				chaffStart.daytime = DayTime;
				chaffStart.done = true;
			}
		}
		if (!chaffEnd.done)
		{
			if (TTd_joint >= TTd_2_chaffEnd)
			{
				chaffEnd.daytime = DayTime;
				chaffEnd.done = true;
			}
		}

		//Z low temperature conditions for the tiller
		float base_Tmpr = 5.0f;
		if (doy < 50 || doy > 300) {
			if (T_ave < base_Tmpr) { Cold_Time += DAYPERHOUR; }
		}
		if (singleRidge.done) { TTd_since_singleRidge += Cur_TTd; }
		if (jointStart.done) { TTd_since_jointing += Cur_TTd; }
		else
		{
			Cold_Time_Ratio_Plant = std::min(Cold_Time / TTd_plant, 1.0f);
			Cold_Time_Ratio_Joint = std::min(Cold_Time / TTd_joint, 1.0f);
		}
	}
}

//Z this compute plant dead fraction, from 0.0 to 1.0
void WheatDevelopment::PlantDie(float top)
{
	float lobnd = 0.0f, upbnd = 1.0f;
	LowKill_Frac = std::max(FracEmerg_Real, HighKill_Frac);
	HighKill_Frac = top;
	if (HighKill_Frac < upbnd)
	{
		PresentKill_Frac = std::min(HighKill_Frac - LowKill_Frac, HighKill_Frac - lobnd);
	}
	else
	{
		PresentKill_Frac = std::min(upbnd - LowKill_Frac, upbnd - lobnd);
	}
	PresentKill_Frac = std::max(PresentKill_Frac, 0.0f);
	TotalKill_Frac += PresentKill_Frac;
	TotalKill_Frac = std::min(TotalKill_Frac, 1.0f);

	if ((1.0f - TotalKill_Frac) <= 0.000001f)
	{
		//Z plant is totally killed, need to stop this program.
		TotalKill_Frac = 1.0f;
	}
}

//Z germinaiton and emergence gdd curves
//  input gdd and obtain a fraction
float WheatDevelopment::EmCurv(float xval, int phase)
{
	float b1 = 0.319381530f, b2 = -0.356563782f, b3 = 1.781477937f, b4 = -1.821255978f, 
		b5 = 1.330274429f, p = 0.2316419f, pi = 3.1415927f;
	float dev = 0.0f;
	if (phase == 1) { dev = (xval - mgerm[temcon]) / sgerm[temcon]; }	//Z germination
	else { dev = (xval - memerg[temcon]) / semerg[temcon]; }// emergence
	float abdev = fabs(dev);
	float t1 = 1.0f / (1.0f + p * abdev);
	float yval = 1.0f - 1.0f / sqrt(2.0f * pi) * exp(-pow(abdev, 2.0f) / 2.0f) * t1 * (b1 + t1 * (b2 + t1 * (b3 + t1 * (b4 + t1 * b5))));

	float EmcurvRes = 0.0f;
	if (dev >= 0.0f) { EmcurvRes = yval; }
	else { EmcurvRes = 1.0f - yval; }
	return EmcurvRes;
}

//Z germinaiton and emergence gdd INVERSE curves
//  input a fraction and output gdd
//  with weather condition changes, given a germination/emergence fraction
//  may need to adjust the gdd needed to reach germination/emergence
float WheatDevelopment::ReCurv(int phase, float yval)
{
	int iterno = 0;
	float highx = 0.0f, lowx = 0.0f, medx = 0.0f, tempy = 0.0f, tol = 0.00001f;
	highx = upcute[temcon] + 10.0f;
	float RecurvRes = 0.0f;
	while (true)
	{
		iterno = iterno + 1;
		medx = (highx + lowx) / 2.0f;
		tempy = EmCurv(medx, phase);
		if (abs(yval - tempy) < tol)
		{
			RecurvRes = medx;
			break;
		}
		else if (yval < tempy)
		{
			highx = medx;
		}
		else
		{
			lowx = medx;
		}
	}
	return RecurvRes;
}

//Z This subroutine estimates a (log)linear regression procedure to estimate the duration of
//     kernel growth(DURC)
//     potential finalweight(FWTC)
//     as a function of average temperature during the grain filling period(TAVGR).
//     The equations describing the curves are output to the log file.
//  opt = 0 (option = 0): kernel weight
//  opt = 1 (option = 1): kernel growth duration
void WheatDevelopment::KernelFitCurv(int opt, float* arr)
{
	arr[0] = 0.0f;
	arr[1] = 0.0f;
	arr[2] = 0.0f;
	arr[3] = 0.0f;

	int nbond = -1;

	for (int ii = 0; ii < kerpts; ii++)
	{
		if (kerpot[ii][0] > kertem[opt][0]) { nbond = ii; }
	}

	//Z first line
	float sx = 0.0f, sy = 0.0f, sxy = 0.0f, sx2 = 0.0f;
	for (int ii = 0; ii < nbond; ii++)
	{
		sx += kerpot[ii][0];
		sx2 += (kerpot[ii][0] * kerpot[ii][0]);
		sy += logf(kerpot[ii][opt + 1]);
		sxy += (kerpot[ii][0] * logf(kerpot[ii][opt + 1]));
	}

	arr[1] = (static_cast<float>(nbond) * sxy - sx * sy) / (static_cast<float>(nbond) * sx2 - sx * sx);
	arr[0] = (sy - arr[1] * sx) / static_cast<float>(nbond);

	//Z second line
	sx = 0.0f, sy = 0.0f, sxy = 0.0f, sx2 = 0.0f;
	for (int ii = nbond-1; ii < kerpts; ii++)
	{
		sx += kerpot[ii][0];
		sx2 += (kerpot[ii][0] * kerpot[ii][0]);
		sy += logf(kerpot[ii][opt + 1]);
		sxy += (kerpot[ii][0] * logf(kerpot[ii][opt + 1]));
	}
	
	nbond = kerpts - nbond;
	arr[3] = (static_cast<float>(nbond) * sxy - sx * sy) / (static_cast<float>(nbond) * sx2 - sx * sx);
	arr[2] = (sy - arr[3] * sx) / static_cast<float>(nbond);


}