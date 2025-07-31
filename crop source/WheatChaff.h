#pragma once
#ifndef _WHEAT_CHAFF_
#define _WHEAT_CHAFF_

#include <cmath>
#include <algorithm>
#include <tuple>

#define MAXCHAFFWEIGHT 0.45f  //Z max chaff weight in g
#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.000694f
#define PHYLLOCHRON 106.0f
#define HALFPHYLLOCHRON 53.0f
#define MAX_N_PCT 4.0f
#define MIN_N_PCT 0.5f

using namespace std;

//Z: Wheat Chaff Data Structure
struct WheatChaff {

	//Z ----------- INDIVIDUAL chaff ---------------------------------------

		//Z Growing time/gdd measure
	bool mature;                //Z mark the max chaff lengh is reached
	bool dead;
	bool force_to_death_current_step;

	float cur_TTd;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	float physAge;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	float cumu_growthTdd;		//Z cumu gdd growth age, stressed
	float elongAge;				//						unstressed (this is NOT the elongation age for internodes)

	//Z Chaff Mass g
	float CfMass;
	float deadCfMass;
	float CfMassIncrease;
	float ptnCfMassIncrease;

	//Z Chaff Nitrogen Mass mg
	float CfNitrogenMass;
	float CfNitrogenIncrease;
	float ptnCfNitrogenMassIncrease;

	//Z Chaff Nitrogne conent and release fraction
	float CfNitrogenContent;
	float CfNitrogenReleaseLowBdd;		//Z the min chaff N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	float CfNitrogenReleaseMaxPtge;		//Z the chaff N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	float CfNitrogenRelease;
	float CfNitrogenDead;				//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	float CfNitrogenReturn;				//Z stepwise returned portion of leaf nitrogen, will supply plant future usage 
	float CfNitrogenMassDeadCumu;		//Z cumulative dead portion of leaf nitrogen, will be in the dead plant tissue and not removable


	//Z environemtal data
	float N_effect;

	//Z ----------- REPRESENTATIVE Chaff ---------------------------------------

		//Z living fraction
	float livingFrac;
	float livingFrac_old;
	float livingFrac_ini;

	//Z representative mass values
	float CfMass_Rep;
	float deadCfMass_Rep;
	float CfNitrogenMass_Rep;

	//Z for chaff, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	float CfNitrogenReturn_Rep;
	float deadCfNitrogenMass_Rep;

	float ptnCfMassIncrease_Rep;
	float ptnCfNitrogenMassIncrease_Rep;

	//Z ----------- WHEAT chaff CONSTRUCTOR ---------------------------------------
	WheatChaff();
	WheatChaff(float livingFrac);
};

//Z Wheat Chaff Operator

//Z Wheat Chaff Update Operator
struct WheatChaffUpdate {

	//Z input for this operator should be listed as a member variable
	//Z "WheatPlant" should call "SetChaffCondition" function and update the values 
	float current_gdd;			//Z hourly gdd incremental
	float predawn_psi;			//Z predawn water potential in Mpa
	float growthDuration;       //Z use simple gdd function for chaff growth, 
	//  suggest to compute the growth during in development class for just once.

	bool chaffEnd;				//Z this marks the end of the chaff growth, halfway to the final maturity

	//Z "WheatChaff" should compute the following effects before executing the update pipeline
	//  N effect / stl effect on chaff will vary among chaffs, so will be local variables
	float psi_effect_elong;

	//Z ----------- WHEAT chaff operator CONSTRUCTOR ---------------------------------------
	explicit WheatChaffUpdate();
	void SetChaffUpdateCondition(float gdd, float psi_predawn, float growthDuration, bool chaffEnd);
	/*Z parameters : 1. Internode structure;
	*                2. current living fraction (from tiller);
	*		         3. force_to_death by the plant;
	*/
	//void operator() (std::tuple<WheatChaff&, float, bool> t)
	void ChafUpdateRun(WheatChaff& wChaf, float lf, bool f2d) const;
};

//Z Wheat Chaff Mass Distribution Operator
struct WheatChaffMassAssignment {

	//Z input for this operator should be listed as a member variable
	// 	the variable listed here should be universal among chaffs, under the whole plant scale.

	//Z "WheatPlant" should call "SetChaffMassAssignmentCondition" function and update the values 
	float biomassIncomeRate;	//Z biomass and N increase rates from the plant level mass budget
	float nitrogenIncomeRate;

	//Z funtional operators
	WheatChaffMassAssignment() : biomassIncomeRate(0.0f), nitrogenIncomeRate(0.0f) {}
	void SetChaffMassAssignmentCondition(float biomassRate, float nitrogenRate);
	void ChafMassRun(WheatChaff& wChaf) const;
};

//**********************************************
//Z Lambda function for transformation and reduced sum
//  but we encapsulate the lambda function within structures
struct CfExtRep_CfMass { float operator() (const WheatChaff& wc) const { return wc.CfMass_Rep; } };
struct CfExtRep_deadCfMass { float operator() (const WheatChaff& wc) const { return wc.deadCfMass_Rep; } };
struct CfExtRep_CfNitrogenMass { float operator() (const WheatChaff& wc) const { return wc.CfNitrogenMass_Rep; } };
struct CfExtRep_CfNitrogenReturn { float operator() (const WheatChaff& wc) const { return wc.CfNitrogenReturn_Rep; } };
struct CfExtRep_deadCfNitrogenMass { float operator() (const WheatChaff& wc) const { return wc.deadCfNitrogenMass_Rep; } };
struct CfExtRep_ptnCfMassIncrease { float operator() (const WheatChaff& wc) const { return wc.ptnCfMassIncrease_Rep; } };
struct CfExtRep_ptnCfNitrogenMassIncrease { float operator() (const WheatChaff& wc) const { return wc.ptnCfNitrogenMassIncrease_Rep; } };

#endif