#pragma once

#ifndef _WHEAT_SPIKELET_H_
#define _WHEAT_SPIKELET_H_

#define MAXFLWRNUM 10
#define PHYLLOCHRON 106.0f
#define FLWRGROPERPHYLLOCHRON 2.5f	//Z flower primordia initiated per phyllochron per spikelet
#define FLWRDIEPERPHYLLOCHRON 1.34f	//Z flower primordia aborted per phyllochron per spikelet
#define TTDBASELFLORET 9.1f			//Z Gdd required for basal floret pair fertilization
#define TTDNONBASELFLORET 27.3f		//Z Gdd required for non-basal floret fertilization
#define KERNELNITROGENCONTENTMAX 3.5f //Z usually 2.5

#define MAX_N_PCT 3.5f				//Z plant level N mass fraction (%), max, min
#define MIN_N_PCT 0.5f

#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;

//Z spikelet is the "suborgan" from the tiller but "mainstem" for flowers and grains, this class is for each spikelet
//Z New spikelet initiation is controlled by the tiller
//Z for this class, it controls the existence of the spikelet and flowers
//      function: 1. flower primordium initiation and abortion position
//                2. flower fertilization with kernel initiated
//                3. initiate kernel that will grow ultramately

//Z all the mass is carried by the kernel
//  the flower is more or less as a phenology event, rather than something really carries biomass and N mass

struct WheatSpikelet
{
	//Z some property about the spikelet itself
	int rank;							//Z rank of this spikelet on the parent tiller, start from 1 (in the revised version, only mainstem tiller has number "0")
	float livingFrac;
	float livingFrac_old;
	bool force_to_death_current_step;
	bool dead;

	float cur_TTd;
	float phyAge;
	float phyAge_flowerIni;
	float phyAge_flowerAbort;
	float phyAge_flowerFert;

	float SpikSinkStrength;

	//Z this marks the behavior of flowering on that spikelet
	bool flowerIni;
	bool flowerAbort;
	bool flowerFert;

	//Z Wheat Flower matrix
	//	Wheat flower mainly include 0-1 for status identifers
	//Z some notation: Basel Flower: the first flower for each spikelet, technically grows on the tiller at the bottom of that spikelet
	//                 Non-Basel Flower: the flower growing on spikelet, growing on the spikelet.
	//                 Flower is counted from 0
	int FlwrNum;
	int FlwrInitNum;
	int FlwrAbortNum;
	int FlwrFertNum;
	int KrNum;

	//Z the following three should be bool but for easy calculation, we use int
	//  1: flower exists and good
	//  2: flower not exists or not fertilized
	//     therefore, FWINIT and FWFERT are initialized to 0 and FWABORT is initialized to 1 (and then changed to 0 for flower abortion)
	bool FwInit[MAXFLWRNUM];		//Z if the flower primordium is initiated
	bool FwAbort[MAXFLWRNUM];		//Z if the flower primordium is abortion
	bool FwFert[MAXFLWRNUM];		//Z if the flower is fertilized
	int  KrInitDay[MAXFLWRNUM];		//Z kernel intiation hour since the first day of kernel initiation/flower fertilization, relative to the time for the first flower fertilization
	bool KrGroComplete[MAXFLWRNUM]; //Z kernel growth complete, 0 not complete, 1 complete

	//Z Wheat Kernel mass and N mass data and matrix
	//  each flower is one kernel
	float KrMassSum;
	float KrNitrogenMassSum;
	float KrMass[MAXFLWRNUM];
	float KrNitrogenMass[MAXFLWRNUM];
	float maxKrMass[MAXFLWRNUM];
	float maxKrDuration[MAXFLWRNUM];
	float KrDuration[MAXFLWRNUM];
	float KrMassIncrease[MAXFLWRNUM];
	float KrNitrogenMassIncrease[MAXFLWRNUM];

	float ptnKrMassIncrease[MAXFLWRNUM];
	float ptnKrNitrogenMassIncrease[MAXFLWRNUM];
	float ptnKrMassIncreaseSum;
	float ptnKrNitrogenMassIncreaseSum;

	//Z representative kernel mass and nitrogen mass values
	float KrMassSum_Rep;
	float KrNitrogenMassSum_Rep;
	float ptnKrMassIncreaseSum_Rep;
	float ptnKrNitrogenMassIncreaseSum_Rep;

	//Z representative flower and kernel number
	//Z counting flower and kernel number is different from leaf and internode
	//  because leaf can die and disappear, but flower/kernel, especially kernel, when exists, always exists (may be stop growing and become low quality)
	//  so flower and kernel numbers should always be a cumulative value.
	float KrNum_Rep;

	//Z old mass values use to compute incremental
	int KrNum_old;
	float KrMassSum_old;
	float KrNitrogenMassSum_old;

	//Z ----------- WHEAT SPIKELET CONSTRUCTOR ---------------------------------------	
	WheatSpikelet();
	WheatSpikelet(int rank, float livingFrac);
};

//Z Wheat Leaf Operator

//Z Wheat Leaf Update Operator
struct WheatSpikeletUpdate
{
	//Z ambient parameters
		//Z data from development module
	float current_tmpr;
	float current_gdd;
	float flowergro_gdd;
	float flowerdie_gdd;
	float baselFlwrFert_gdd;
	float nonbaselFlwrFert_gdd;

	//Z data from development module
	float DayTime;
	float DayTimekernelFertStart;
	float dayAfterAnthesisStart;

	const float sinkst_basel = 0.1f;
	const float sinkf[2][MAXFLWRNUM] = {
		{1.0f, 1.0f, 0.7f, 0.5f, 0.3f, 0.2f, 0.1f, 0.1f, 0.1f, 0.1f},
		{1.0f, 0.8f, 0.5f, 0.4f, 0.2f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f} };
	float PltBiomass;
	float PltNitrogenMass;
	float Pltpsi;
	float N_effect;
	float psi_effect;
	float rai;		//Z N_effect*psi_effect
	const float raiopt = 1.0f;
	const float raiwc = 0.2f;
	float KrSinkStrengthAdjustCoef[MAXFLWRNUM];

	//Z kernel weight and growth duration coefficents
	//  we use log-linear two-piece functionals
	float kwt00, kwt01, kwt10, kwt11;
	float kdur00, kdur01, kdur10, kdur11;
	const float kertem[2][2] = { {25.0f,36.0f},{25.0f,36.0f} };		//Z Temperature breakpoints and cutoffs for kernal growth curves

	//Z do a 30-day moving average of tmpr (for kernel)
	//  this should definitely be the air tmpr because it is after elongation and all the wheat should grow up very well
	//  this is daily averaged temperature, each index means one day
	float TRecordNumberDay;
	float TaveDaily;
	float Tave30;
	float TRecord30[30];
	float TaveRecord30[30];
	int Tcurrentposition30;
	float TRecordNumber30;

	//Z funtional operators
	WheatSpikeletUpdate();
	void SetSpikeletUpdateCondition(float gdd, float tmpr, float psi_predawn,
		float pltbiomass, float pltnitrogenmass,
		float daytime, float dayfirstKrStart, float dayafterAssStar,
		float kwt00, float kwt01, float kwt10, float kwt11,
		float kdur00, float kdur01, float kdur10, float kdur11);

	void TmprMovingAve30(float current_tmpr);
	//Z this is the main function for plant spikelet growth
	//  every time updating, hourly gdd, tmpr, water potential, shade are universal
	//                       N stress are local
	//                       living fractor is the input variable for each spikelet object
	/*Z parameters : 1. spikelet structure;
	*                2. current living fraction (from tiller);
	*                3. force_to_death by the plant;
	*				 4-6. flower initiation, abortion, and fertilization stage for that spikelet
	*/
	//void operator() (std::tuple<WheatSpikelet&, float, bool, bool, bool, bool> t)
	void SpikUpdateRun(WheatSpikelet& wSpik, float lf, bool f2d, bool flwrInit, bool flwrAbort, bool flwrFert) const;
};


//Z Wheat Spikelet (Kernel) Mass Distribution Operator
struct WheatSpikeletMassAssignment {

	//Z input for this operator should be listed as a member variable
	// 	the variable listed here should be universal among spikelet, under the whole plant scale.
	float biomassIncomeRate;	//Z biomass and N increase rates from the plant level mass budget
	float nitrogenIncomeRate;
	//Z funtional operators
	WheatSpikeletMassAssignment() : biomassIncomeRate(0.0f), nitrogenIncomeRate(0.0f) {}
	void SetSpikeletMassAssignmentCondition(float biomassRate, float nitrogenRate);
	void SpikMassRun(WheatSpikelet& wSpik) const;
};

//**********************************************
//Z Lambda function for transformation and reduced sum
//  but we encapsulate the lambda function within structures
struct SpExtRep_KrNum { float operator() (const WheatSpikelet& ws) const { return ws.KrNum_Rep; } };
struct SpExtRep_KrMassSum { float operator() (const WheatSpikelet& ws) const { return ws.KrMassSum_Rep; } };
struct SpExtRep_KrNitrogenMassSum { float operator() (const WheatSpikelet& ws) const { return ws.KrNitrogenMassSum_Rep; } };
struct SpExtRep_ptnKrMassIncreaseSum { float operator() (const WheatSpikelet& ws) const { return ws.ptnKrMassIncreaseSum_Rep; } };
struct SpExtRep_ptnKrNitrogenMassIncreaseSum { float operator() (const WheatSpikelet& ws) const { return ws.ptnKrNitrogenMassIncreaseSum_Rep; } };

#endif