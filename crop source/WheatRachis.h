#pragma once
#ifndef _WHEAT_RACHIS_
#define _WHEAT_RACHIS_

#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;

#define MAXRACHISLENGTH 10.0f		//Z max rachis length in cm, (100 mm = 10 cm), named maxral
#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.000694f
#define PHYLLOCHRON 106.0f
#define HALFPHYLLOCHRON 53.0f
#define MAX_N_PCT 4.0f
#define MIN_N_PCT 0.5f

//Z: Wheat Rachis Data Structure
struct WheatRachis {

	//Z ----------- INDIVIDUAL Rachis ---------------------------------------

	//Z Growing time/gdd measure
	bool mature;                //Z mark the max rachis lengh is reached
	bool dead;
	bool force_to_death_current_step;

	float cur_TTd;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	float physAge;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	float cumu_growthTdd;		//Z cumu gdd growth age, stressed
	float elongAge;				//Z unstressed growth phyAge

	//Z Rachis Length (cm)
	float RaLength;				//Z rachis length in cm
	float ptnRaLengthIncrease;  //Z potential rachis length increase
	float RaLengthIncrease;		//Z actual Rachis length increase

	//Z Rachis Mass g
	float srl;					//Z rachis specific weight g cm^-1
	float RaMass;
	float deadRaMass;
	float RaMassIncrease;
	float ptnRaMassIncrease;

	//Z Rachis Nitrogen Mass mg
	float RaNitrogenMass;
	float RaNitrogenIncrease;
	float ptnRaNitrogenMassIncrease;

	//Z Rachis Nitrogne conent and release fraction
	float RaNitrogenContent;
	float RaNitrogenReleaseLowBdd;		//Z the min rachis N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	float RaNitrogenReleaseMaxPtge;		//Z the rachis N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	float RaNitrogenRelease;
	float RaNitrogenDead;				//Z stepwise dead portion of rachis nitrogen, will be in the dead plant tissue and not removable
	float RaNitrogenReturn;				//Z stepwise returned portion of rachis nitrogen, will supply plant future usage 
	float RaNitrogenMassDeadCumu;		//Z cumulative dead portion of rachis nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	float N_effect;

	//Z ----------- REPRESENTATIVE Rachis ---------------------------------------

	//Z Rachis length
	float RaLength_Rep;	//Z this should be a cumulative value, as the living fraction changes, 
	//  it cumulated different amount (based on living fraction)

//Z living fraction
	float livingFrac;
	float livingFrac_old;
	float livingFrac_ini;

	//Z representative mass values
	float RaMass_Rep;
	float deadRaMass_Rep;
	float RaNitrogenMass_Rep;

	//Z for Rachis, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	float RaNitrogenReturn_Rep;
	float deadRaNitrogenMass_Rep;

	float ptnRaMassIncrease_Rep;
	float ptnRaNitrogenMassIncrease_Rep;

	//Z ----------- WHEAT rachis CONSTRUCTOR ---------------------------------------
	WheatRachis();
	WheatRachis(float livingFrac);
};

//Z Wheat Rachis Operator

//Z Wheat Rachis Update Operator
struct WheatRachisUpdate {

	//Z input for this operator should be listed as a member variable
	//Z "WheatPlant" should call "SetRachisCondition" function and update the values 
	float current_gdd;			//Z hourly gdd incremental
	float predawn_psi;			//Z predawn water potential in Mpa
	float growthDuration;		//Z jtpa*pchron is the growth period, rachis growth between double ridge and jointing
	float growthRate;			//Z use simple gdd function for rachis growth, should be a negative value as the initial value
	//Z only compute once at double ridge
//  suggest to compute the growth during in development class for just once.

	bool rachisEnd;				//Z this marks the end of the rachis growth,
	//Z at the time of jointing

//Z "WheatRachis" should compute the following effects before executing the update pipeline
//  N effect / stl effect on rachis will vary among rachis, so will be local variables
	float psi_effect_elong;

	//Z ----------- WHEAT rachis operator CONSTRUCTOR ---------------------------------------
	explicit WheatRachisUpdate();
	void SetRachisUpdateCondition(float gdd, float psi_predawn, float growthDuration, bool rachisEnd);
	/*Z parameters : 1. Internode structure;
	*                2. current living fraction (from tiller);
	*		         3. force_to_death by the plant;
	*/
	//void operator() (std::tuple<WheatRachis&, float, bool> t)
	void RchsUpdateRun(WheatRachis& wRchs, float lf, bool f2d) const;
};

//Z Wheat Rachis Mass Distribution Operator
struct WheatRachisMassAssignment {

	//Z input for this operator should be listed as a member variable
	// 	the variable listed here should be universal among rachis, under the whole plant scale.

	//Z "WheatPlant" should call "SetRachisMassAssignmentCondition" function and update the values 
	float biomassIncomeRate;	//Z biomass and N increase rates from the plant level mass budget
	float nitrogenIncomeRate;

	//Z funtional operators
	WheatRachisMassAssignment() : biomassIncomeRate(0.0f), nitrogenIncomeRate(0.0f) {}
	void SetRachisMassAssignmentCondition(float biomassRate, float nitrogenRate);
	void RchsMassRun(WheatRachis& wRchs) const;
};

//**********************************************
//Z Lambda function for transformation and reduced sum
//  but we encapsulate the lambda function within structures
struct RaExtRep_RaLength { float operator() (const WheatRachis& wr) const { return wr.RaLength_Rep; } };
struct RaExtRep_RaMass { float operator() (const WheatRachis& wr) const { return wr.RaMass_Rep; } };
struct RaExtRep_deadRaMass { float operator() (const WheatRachis& wr) const { return wr.deadRaMass_Rep; } };
struct RaExtRep_RaNitrogenMass { float operator() (const WheatRachis& wr) const { return wr.RaNitrogenMass_Rep; } };
struct RaExtRep_RaNitrogenReturn { float operator() (const WheatRachis& wr) const { return wr.RaNitrogenReturn_Rep; } };
struct RaExtRep_deadRaNitrogenMass { float operator() (const WheatRachis& wr) const { return wr.deadRaNitrogenMass_Rep; } };
struct RaExtRep_ptnRaMassIncrease { float operator() (const WheatRachis& wr) const { return wr.ptnRaMassIncrease_Rep; } };
struct RaExtRep_ptnRaNitrogenMassIncrease { float operator() (const WheatRachis& wr) const { return wr.ptnRaNitrogenMassIncrease_Rep; } };

#endif