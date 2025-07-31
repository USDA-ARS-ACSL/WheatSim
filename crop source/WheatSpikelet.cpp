#include "WheatSpikelet.h"

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Spikelet Data Structure

WheatSpikelet::WheatSpikelet()
{
	rank = 1;
	livingFrac = 1.0f;
	livingFrac_old = 1.0f;
	force_to_death_current_step = false;
	dead = false;

	cur_TTd = 0.0f;
	phyAge = 0.0f;
	phyAge_flowerIni = 0.0f;
	phyAge_flowerAbort = 0.0f;
	phyAge_flowerFert = 0.0f;

	SpikSinkStrength = 1.0f;

	//Z this marks the behavior of flowering on that spikelet
	flowerIni = false;
	flowerAbort = false;
	flowerFert = false;

	//Z Wheat Flower matrix
	//	Wheat flower mainly include 0-1 for status identifers
	//Z some notation: Basel Flower: the first flower for each spikelet, technically grows on the tiller at the bottom of that spikelet
	//                 Non-Basel Flower: the flower growing on spikelet, growing on the spikelet.
	//                 Flower is counted from 0
	FlwrNum = 0;
	FlwrInitNum = 0;
	FlwrAbortNum = 0;
	FlwrFertNum = 0;
	KrNum = 0;
	//Z the following three should be bool but for easy calculation, we use int
	//  1: flower exists and good
	//  2: flower not exists or not fertilized
	//     therefore, FWINIT and FWFERT are initialized to 0 and FWABORT is initialized to 1 (and then changed to 0 for flower abortion)
	for (int ii = 1; ii < MAXFLWRNUM; ii++)
	{
		FwInit[ii] = false;			//Z if the flower primordium is initiated
		FwAbort[ii] = true;			//Z if the flower primordium is abortion
		FwFert[ii] = false;			//Z if the flower is fertilized
		KrInitDay[ii] = 0;			//Z kernel intiation hour since the first day of kernel initiation/flower fertilization, relative to the time for the first flower fertilization
		KrGroComplete[ii] = false;	//Z kernel growth complete, 0 not complete, 1 complete
	}

	//Z Wheat Kernel mass and N mass data and matrix
	//  each flower is one kernel
	KrMassSum = 0.0f;
	KrNitrogenMassSum = 0.0f;
	ptnKrMassIncreaseSum = 0.0f;
	ptnKrNitrogenMassIncreaseSum = 0.0f;
	for (int ii = 1; ii < MAXFLWRNUM; ii++)
	{
		KrMass[ii] = 0.0f;
		KrNitrogenMass[ii] = 0.0f;
		maxKrMass[ii] = 0.0f;
		maxKrDuration[ii] = 0.0f;
		KrDuration[ii] = 0.0f;
		KrMassIncrease[ii] = 0.0f;
		KrNitrogenMassIncrease[ii] = 0.0f;

		ptnKrMassIncrease[ii] = 0.0f;
		ptnKrNitrogenMassIncrease[ii] = 0.0f;
	}

	//Z representative kernel mass and nitrogen mass values
	KrNum_Rep = 0.0f;
	KrMassSum_Rep = 0.0f;
	KrNitrogenMassSum_Rep = 0.0f;
	ptnKrMassIncreaseSum_Rep = 0.0f;
	ptnKrNitrogenMassIncreaseSum_Rep = 0.0f;

	//Z old mass values use to compute incremental
	KrNum_old = 0;
	KrMassSum_old = 0.0f;
	KrNitrogenMassSum_old = 0.0f;
}

WheatSpikelet::WheatSpikelet(int rank, float livingFrac) :
	rank(rank), livingFrac(livingFrac), livingFrac_old(livingFrac)
{
	force_to_death_current_step = false;
	dead = false;

	cur_TTd = 0.0f;
	phyAge = 0.0f;
	phyAge_flowerIni = 0.0f;
	phyAge_flowerAbort = 0.0f;
	phyAge_flowerFert = 0.0f;

	SpikSinkStrength = 1.0f;

	//Z this marks the behavior of flowering on that spikelet
	flowerIni = false;
	flowerAbort = false;
	flowerFert = false;

	//Z Wheat Flower matrix
	//	Wheat flower mainly include 0-1 for status identifers
	//Z some notation: Basel Flower: the first flower for each spikelet, technically grows on the tiller at the bottom of that spikelet
	//                 Non-Basel Flower: the flower growing on spikelet, growing on the spikelet.
	//                 Flower is counted from 0
	FlwrNum = 0;
	FlwrInitNum = 0;
	FlwrAbortNum = 0;
	FlwrFertNum = 0;
	KrNum = 0;
	//Z the following three should be bool but for easy calculation, we use int
	//  1: flower exists and good
	//  2: flower not exists or not fertilized
	//     therefore, FWINIT and FWFERT are initialized to 0 and FWABORT is initialized to 1 (and then changed to 0 for flower abortion)
	for (int ii = 1; ii < MAXFLWRNUM; ii++)
	{
		FwInit[ii] = false;			//Z if the flower primordium is initiated
		FwAbort[ii] = true;			//Z if the flower primordium is abortion
		FwFert[ii] = false;			//Z if the flower is fertilized
		KrInitDay[ii] = 0;			//Z kernel intiation hour since the first day of kernel initiation/flower fertilization, relative to the time for the first flower fertilization
		KrGroComplete[ii] = false;	//Z kernel growth complete, 0 not complete, 1 complete
	}

	//Z Wheat Kernel mass and N mass data and matrix
	//  each flower is one kernel
	KrMassSum = 0.0f;
	KrNitrogenMassSum = 0.0f;
	ptnKrMassIncreaseSum = 0.0f;
	ptnKrNitrogenMassIncreaseSum = 0.0f;
	for (int ii = 1; ii < MAXFLWRNUM; ii++)
	{
		KrMass[ii] = 0.0f;
		KrNitrogenMass[ii] = 0.0f;
		maxKrMass[ii] = 0.0f;
		maxKrDuration[ii] = 0.0f;
		KrDuration[ii] = 0.0f;
		KrMassIncrease[ii] = 0.0f;
		KrNitrogenMassIncrease[ii] = 0.0f;

		ptnKrMassIncrease[ii] = 0.0f;
		ptnKrNitrogenMassIncrease[ii] = 0.0f;
	}

	//Z representative kernel mass and nitrogen mass values
	KrNum_Rep = 0.0f;
	KrMassSum_Rep = 0.0f;
	KrNitrogenMassSum_Rep = 0.0f;
	ptnKrMassIncreaseSum_Rep = 0.0f;
	ptnKrNitrogenMassIncreaseSum_Rep = 0.0f;

	//Z old mass values use to compute incremental
	KrNum_old = 0;
	KrMassSum_old = 0.0f;
	KrNitrogenMassSum_old = 0.0f;
}


// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Spikelet Update Operator Structure
//  the constructor is define as follows
WheatSpikeletUpdate::WheatSpikeletUpdate()
{
	current_tmpr = 0.0f;
	current_gdd = 0.0f;
	flowergro_gdd = PHYLLOCHRON / FLWRGROPERPHYLLOCHRON;
	flowerdie_gdd = PHYLLOCHRON / FLWRDIEPERPHYLLOCHRON;
	baselFlwrFert_gdd = TTDBASELFLORET;
	nonbaselFlwrFert_gdd = TTDNONBASELFLORET;

	//Z data from development module
	DayTime = 0.0f;
	DayTimekernelFertStart = 0.0f;
	dayAfterAnthesisStart = 0.0f;

	PltBiomass = 0.0f;
	PltNitrogenMass = 0.0f;
	Pltpsi = -0.5f;

	N_effect = 1.0f;
	psi_effect = 1.0f;
	rai = 1.0f;

	for (int ii = 0; ii < MAXFLWRNUM; ii++) { KrSinkStrengthAdjustCoef[ii] = sinkf[0][ii]; }

	Tave30 = 0.0f;
	TaveDaily = 0.0f;
	for (int ii = 0; ii < 30; ii++) { TRecord30[ii] = 0.0f; TaveRecord30[ii] = 0.0f; }
	Tcurrentposition30 = 0;
	TRecordNumber30 = 0.0f;
	TRecordNumberDay = 0.0f;

	this->kwt00 = 0.0f;
	this->kwt01 = 0.0f;
	this->kwt10 = 0.0f;
	this->kwt11 = 0.0f;
	this->kdur00 = 0.0f;
	this->kdur01 = 0.0f;
	this->kdur10 = 0.0f;
	this->kdur11 = 0.0f;

}

//Z update leaf growth gdd and conditions for this current hourly step
//  leaf stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatSpikeletUpdate::SetSpikeletUpdateCondition(float gdd, float tmpr,
	float psi_predawn, float pltbiomass, float pltnitrogenmass,
	float daytime, float dayfirstKrStart, float dayafterAssStar,
	float kwt00, float kwt01, float kwt10, float kwt11,
	float kdur00, float kdur01, float kdur10, float kdur11)
{
	current_tmpr = tmpr;
	current_gdd = gdd;

	Pltpsi = psi_predawn;
	PltBiomass = pltbiomass;
	PltNitrogenMass = pltnitrogenmass;

	DayTime = daytime;
	DayTimekernelFertStart = dayfirstKrStart;
	dayAfterAnthesisStart = dayafterAssStar;

	this->kwt00 = kwt00;
	this->kwt01 = kwt01;
	this->kwt10 = kwt10;
	this->kwt11 = kwt11;
	this->kdur00 = kdur00;
	this->kdur01 = kdur01;
	this->kdur10 = kdur10;
	this->kdur11 = kdur11;

	//Z first compute the plant level N stress
	//Z N in mass fraction (mg/g biomass)
	N_effect = 2.0f / (1.0f + expf(-2.9f * (max(MIN_N_PCT, 10.0f * PltNitrogenMass / PltBiomass) - MAX_N_PCT))) - 1.0f;
	N_effect = max(min(N_effect, 1.0f), 0.1f);

	//Z water effect
	//Z some constant first, use the "leaf expansion parameters"
	const float psi_f = -1.4251f;	// -1.0, was -1.4251   later changed to -2.3;
	const float s_f = 0.4258f;		// was 0.4258 0.5;
	//Z expanding
	const float psi_threshold_bars = -0.8657f;	//Mpa
	psi_effect = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (Pltpsi - psi_threshold_bars))));
	psi_effect = min(max(psi_effect, 0.1f), 1.0f);

	//Z compute the rai value
	//  raiwc is a constant for the "worst condition"
	rai = max(min(N_effect, psi_effect), raiwc);
	float rai_slope = (expf(rai) - expf(raiwc)) / (expf(raiopt) - expf(raiwc));

	for (int ii = 0; ii < MAXFLWRNUM; ii++) { KrSinkStrengthAdjustCoef[ii] = rai_slope * (sinkf[0][ii] - sinkf[1][ii]) + sinkf[1][ii]; }

	//Z update moving temperature
	this->TmprMovingAve30(current_tmpr);

}

//Z this is the main function for plant spikelet growth
//  every time updating, hourly gdd, tmpr, water potential, shade are universal
//                       N stress are local
//                       living fractor is the input variable for each spikelet object
//  variable list: 1 Spikelet 
//                 2 livingfraction
//                 3 force_to_death by the plant;
//				   4-6 flower initiation, abortion, and fertilization stage for that spikelet	
void WheatSpikeletUpdate::SpikUpdateRun(WheatSpikelet& wSpik, float lf, bool f2d, bool flwrInit, bool flwrAbort, bool flwrFert)
const {

	//Z Extract the flower status
	//WheatSpikelet& wSpik = std::get<0>(t);
	//wSpik.livingFrac = std::get<1>(t);
	//wSpik.force_to_death_current_step = std::get<2>(t);
	//wSpik.flowerIni = std::get<3>(t);
	//wSpik.flowerAbort = std::get<4>(t);
	//wSpik.flowerFert = std::get<5>(t);

	wSpik.livingFrac = lf;
	wSpik.force_to_death_current_step = f2d;
	wSpik.flowerIni = flwrInit;
	wSpik.flowerAbort = flwrAbort;
	wSpik.flowerFert = flwrFert;

	wSpik.cur_TTd = current_gdd;
	wSpik.phyAge += current_gdd;

	//Z reset some parameters
	wSpik.ptnKrMassIncreaseSum = 0.0f;
	wSpik.ptnKrNitrogenMassIncreaseSum = 0.0f;
	wSpik.ptnKrMassIncreaseSum_Rep = 0.0f;
	wSpik.ptnKrNitrogenMassIncreaseSum_Rep = 0.0f;

	//Step Zero Force To Death This Step
	//Z this condition is imposed by the parent tiller, stronger than any other conditions
	//Z we will assume the kernel and flower will stop growing, but the organs already exist will not drop, grains there is still there
	if ((!wSpik.dead) && wSpik.force_to_death_current_step) {
		wSpik.flowerIni = false;
		wSpik.flowerAbort = false;
		wSpik.flowerFert = false;
		for (int ii = 0; ii < MAXFLWRNUM; ii++) { wSpik.KrGroComplete[ii] = true; }
	}

	//Step One Flower Initiation And Abortion
	if ((!wSpik.dead) && (wSpik.flowerIni && (!wSpik.flowerAbort)))
	{
		wSpik.phyAge_flowerIni += current_gdd;
		if (wSpik.phyAge_flowerIni >= flowergro_gdd && wSpik.FlwrNum < MAXFLWRNUM)
		{
			wSpik.FwInit[wSpik.FlwrInitNum] = true;
			wSpik.FlwrInitNum += 1;
			wSpik.FlwrNum = wSpik.FlwrInitNum;
			wSpik.phyAge_flowerIni -= flowergro_gdd;
		}
	}
	if ((!wSpik.dead) && (wSpik.flowerAbort))
	{
		wSpik.phyAge_flowerAbort += current_gdd;
		if (wSpik.phyAge_flowerAbort >= flowerdie_gdd && wSpik.FlwrNum >= 0)
		{
			wSpik.FwAbort[wSpik.FlwrNum - 1] = false;
			wSpik.FlwrAbortNum += 1;
			wSpik.FlwrNum -= 1;
			wSpik.phyAge_flowerAbort -= flowerdie_gdd;
		}
	}

	//Step Two Flower Fertilizaiton
	//Z the first two spikelet will not fertilize any flowers
	if ((!wSpik.dead) && (wSpik.flowerFert && wSpik.rank > 2))
	{
		//Z the basel flower is fertilzerd as the spikelet is ready to fertilize all its flowers
		if (!wSpik.FwFert[0])
		{
			wSpik.FwFert[0] = true;
			wSpik.SpikSinkStrength = max(1.0f - dayAfterAnthesisStart * sinkst_basel, 0.0f);
			//Z compute the fertilization time (hour) relative to the first flower fertilization time
			wSpik.KrInitDay[0] = max(static_cast<int>(floor((DayTime - DayTimekernelFertStart))), 0);
			wSpik.FlwrFertNum = 1;
			wSpik.KrNum = 1;
		}
		wSpik.phyAge_flowerFert += current_gdd;
		if (wSpik.phyAge_flowerFert >= nonbaselFlwrFert_gdd && wSpik.FlwrFertNum < wSpik.FlwrNum)
		{
			wSpik.FwFert[wSpik.FlwrFertNum] = true;
			//Z compute the fertilization time (hour) relative to the first flower fertilization time
			wSpik.KrInitDay[wSpik.FlwrFertNum] = max(static_cast<int>(floor((DayTime - DayTimekernelFertStart))), 0);
			wSpik.FlwrFertNum += 1;
			wSpik.KrNum += 1;
			wSpik.phyAge_flowerFert -= nonbaselFlwrFert_gdd;
		}
	}

	//Step Three Kernal Growth
		//Z first update the target mass and potential growth duration
	if (!wSpik.dead)
	{
		wSpik.KrMassSum = 0.0f;
		wSpik.KrNitrogenMassSum = 0.0f;
		for (int ii = 0; ii < wSpik.FlwrFertNum; ii++)
		{
			//Z when the "Kernek" growth is stopped, we do not insert biomass and N to the kernel anymore
			//  when the "Leaf" expansion is stopped, we still allow biomass and N filling
			//  should pay attention to this difference
			wSpik.ptnKrMassIncrease[ii] = 0.0f;
			wSpik.ptnKrNitrogenMassIncrease[ii] = 0.0f;
			if (!wSpik.KrGroComplete[ii])
			{
				int n = wSpik.KrInitDay[ii];
				float kernelTave = Tave30;
				if (n < 30) { kernelTave = TaveRecord30[n]; }
				float fwt = (kwt00 + kwt01 * kernelTave) * static_cast<float>(kernelTave < kertem[0][0])
					+ (kwt10 + kwt11 * min(kernelTave, kertem[0][1])) * static_cast<float>(kernelTave >= kertem[0][0]);
				fwt = exp(fwt) * 0.001f;	//Z g kernal biomass
				float dur = (kdur00 + kdur01 * kernelTave) * static_cast<float>(kernelTave < kertem[1][0])
					+ (kdur10 + kdur11 * min(kernelTave, kertem[1][1])) * static_cast<float>(kernelTave >= kertem[1][0]);
				dur = exp(dur);				//Z gdd growing duration in thermal time

				//Z Determine the ratio of the number of growing degree days available
				//  for growth to the number of growing degree days for today's growth.
				//  TDUR will be 1.0 except during the last day of kernel growth.
				float tdur = min(current_gdd, dur - wSpik.KrDuration[ii]) / current_gdd;
				//Z the new max (target) kernel mass
				wSpik.maxKrMass[ii] = wSpik.KrMass[ii] + (fwt / dur) * tdur * current_gdd * wSpik.SpikSinkStrength * KrSinkStrengthAdjustCoef[ii];
				//Z update the kernel growth duration
				wSpik.KrDuration[ii] += tdur * current_gdd;

				//Z determine when the kernal growth is complete
				if (wSpik.KrDuration[ii] >= dur) { wSpik.KrGroComplete[ii] = true; }

				//Z biomass demand (g) and N demand (mg)
				wSpik.ptnKrMassIncrease[ii] = max(wSpik.maxKrMass[ii] - wSpik.KrMass[ii], 0.0f);
				wSpik.ptnKrNitrogenMassIncrease[ii] = max(wSpik.maxKrMass[ii] * 10.0f * KERNELNITROGENCONTENTMAX - wSpik.KrNitrogenMass[ii], 0.0f);

				//Z sum for each spikelet
				wSpik.ptnKrMassIncreaseSum += wSpik.ptnKrMassIncrease[ii];
				wSpik.ptnKrNitrogenMassIncreaseSum += wSpik.ptnKrNitrogenMassIncrease[ii];
			}
			//Z no matter growing or not, mass values should always be summed 
			wSpik.KrMassSum += wSpik.KrMass[ii];
			wSpik.KrNitrogenMassSum += wSpik.KrNitrogenMass[ii];
		}
	}

	//Step Four Representative Mapping -----------------------------------------------------------
	//Z Simulation of representative internode, adjusted by the livingFrac

		//Z representative kernel biomass and kernel N mass
		//  assume kernel will not fall, and mater filled there will not release
		//  that is because kernel filling is almost the end of the plant growth
	wSpik.KrNum_Rep += static_cast<float>(max(wSpik.KrNum - wSpik.KrNum_old, 0)) * wSpik.livingFrac;
	wSpik.KrMassSum_Rep += max(wSpik.KrMassSum - wSpik.KrMassSum_old, 0.0f) * wSpik.livingFrac;
	wSpik.KrNitrogenMassSum_Rep += max(wSpik.KrNitrogenMassSum - wSpik.KrNitrogenMassSum_old, 0.0f) * wSpik.livingFrac;

	//Z representative kernel demand biomass and N
	wSpik.ptnKrMassIncreaseSum_Rep = wSpik.ptnKrMassIncreaseSum * wSpik.livingFrac;
	wSpik.ptnKrNitrogenMassIncreaseSum_Rep = wSpik.ptnKrNitrogenMassIncreaseSum * wSpik.livingFrac;

	//Z update old mass and N mass from previous step
	wSpik.KrMassSum_old = wSpik.KrMassSum;
	wSpik.KrNitrogenMassSum_old = wSpik.KrNitrogenMassSum;
	wSpik.KrNum_old = wSpik.KrNum;

	//Z update the living fraction
	wSpik.livingFrac_old = wSpik.livingFrac;

	if (wSpik.dead || (wSpik.livingFrac == 0.0f))
	{
		wSpik.livingFrac = 0.0f;
		wSpik.livingFrac_old = 0.0f;

		wSpik.ptnKrMassIncreaseSum_Rep = 0.0f;
		wSpik.ptnKrNitrogenMassIncreaseSum_Rep = 0.0f;
	}
}

// Z hourly average tmpr function
//   should be removed because it generates an float array that too large for gpu computation
void WheatSpikeletUpdate::TmprMovingAve30(float airTmpr)
{
	//Z do a 30-day moving average of air tmpr

	//First Upscale hourly temperature to daily temperature
	if (TRecordNumberDay <= 0.0f)
	{
		TaveDaily = airTmpr;
		TRecordNumberDay += 1.0f;
	}
	else
	{
		TaveDaily = (TaveDaily * TRecordNumberDay + airTmpr) / (TRecordNumberDay + 1.0f);
		TRecordNumberDay += 1.0f;
	}
	//Second add daily temperature to a 30 moving average list
	if (TRecordNumber30 <= 0.0f)
	{
		TRecord30[0] = TaveDaily;
		Tave30 = TaveDaily;
		TaveRecord30[0] = TaveDaily;
		Tcurrentposition30 = 1;
		TRecordNumber30 = 1.0f;
	}
	else
	{
		float tmpr_kick_out = TRecord30[Tcurrentposition30];
		float T_sum = Tave30 * TRecordNumber30;
		T_sum = T_sum - tmpr_kick_out + TaveDaily * (TRecordNumberDay + 1.0f) / 24.0f;
		TRecord30[Tcurrentposition30] = TaveDaily;

		if (TRecordNumber30 < 30.0f)
		{
			TRecordNumber30 = std::min(TRecordNumber30 + 1.0f / 24.0f, 30.0f);
		}
		else
		{
			TRecordNumber30 = std::min(TRecordNumber30 - 1.0f + 1.0f / 24.0f, 30.0f);
		}
		Tave30 = T_sum / TRecordNumber30;
		TaveRecord30[Tcurrentposition30] = Tave30;

		Tcurrentposition30 = Tcurrentposition30 + static_cast<int>(std::floor(TRecordNumberDay / 23.9f));
	}
	if (Tcurrentposition30 >= 30) { Tcurrentposition30 = 0; }
	if (TRecordNumberDay > 23.9f) { TRecordNumberDay = 0.0f; }

}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Spikelet (Kernel) Mass Distribution Operator Structure

//Z update spikelet mass and 
//  spikelet stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatSpikeletMassAssignment::SetSpikeletMassAssignmentCondition(float biomassRate, float nitrogenRate)
{
	//Z receive biomass and nitrogen distribution rates for spikelet over the entire plant based plant-level biomass allocation
	biomassIncomeRate = biomassRate;
	nitrogenIncomeRate = nitrogenRate;
}

void WheatSpikeletMassAssignment::SpikMassRun(WheatSpikelet& wSpik)
const {
	//Z reset for each time step
	//  store biomass (g) and N (mg) allocated from plant
	for (int ii = 1; ii < MAXFLWRNUM; ii++)
	{
		wSpik.KrMassIncrease[ii] = 0.0f;
		wSpik.KrNitrogenMassIncrease[ii] = 0.0f;
	}

	//Z separate cases explicitly, good for GPU stream control.
	// BIOMASS
	// biomass incremental
	if (wSpik.ptnKrMassIncreaseSum_Rep > 0.0f)
	{
		//Z this "adjustment" variable is necessary and important,
		//  essentially, " / wSpik.livingFrac" is the important term
		//  this convert the mass allocation for a "representative kernel" to "individual kernel"
		float biomassIncreaseIndividual = biomassIncomeRate / wSpik.livingFrac;
		for (int ii = 1; ii < MAXFLWRNUM; ii++)
		{
			wSpik.KrMassIncrease[ii] = wSpik.ptnKrMassIncrease[ii] * biomassIncreaseIndividual;
			wSpik.KrMass[ii] += wSpik.KrMassIncrease[ii];
		}
	}
	// NITROGEM
	// N incremental
	if (wSpik.ptnKrNitrogenMassIncreaseSum_Rep > 0.0f)
	{
		float nitrogenMassIncreaseIndividual = nitrogenIncomeRate / wSpik.livingFrac;
		for (int ii = 1; ii < MAXFLWRNUM; ii++)
		{
			wSpik.KrNitrogenMassIncrease[ii] = wSpik.ptnKrNitrogenMassIncrease[ii] * nitrogenMassIncreaseIndividual;
			wSpik.KrNitrogenMass[ii] += wSpik.KrNitrogenMassIncrease[ii];
		}
	}
	return;
}