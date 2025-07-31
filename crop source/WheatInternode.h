#pragma once
#ifndef _WHEAT_INTERNODE_
#define _WHEAT_INTERNODE_

#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.000694f
#define PHYLLOCHRON 106.0f
#define SENESDURATION 400.0f
#define HALFPHYLLOCHRON 53.0f
#define MAX_N_PCT 3.5f
#define MIN_N_PCT 0.5f

//Z internode growth parameters from fspm
//Z time (day) difference between leaf and internode from the same position
#define TIMEDIFF_LF_INTR 25.7f//31.7083333f	//Z day (=2739600 sec) of equivalent time between leaf and internode initiation
//Z the initial internode length
#define INTR_LENGTH_INI	0.005f	//Z the initial internode length
//Z  parameter for automate elongation, presented in "s" but cast it to "day" here
#define INTR_TE 27.64615f		//Z end of leaf elongation in automate growth (sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define INTR_TM 21.064983f		//Z time at which leaf elongation rate is maximal in automate growth(sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define INTR_TB -31.190983f		//Z beginning of leaf elongation in automate growth(sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define RER_KC	100.0f			//Z affinity coefficient of RER to C (µmol g-1)
#define RER_KN	15.0f			//Z affinity coefficient of RER to N (µmol g-1)


//Z leaf mass growth parameter
#define IN_ALPHA 0.106f					//Z Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
#define IN_BETA 1.28f					//Z Parameter of the relation between leaf mass and leaf length (dimensionless)
#define IN_DRYMASS_2_STRCTMASS 0.35f	//Z converter between structure mass to dry mass, others can be functional materials
#define IN_STRCTMASS_2_DRYMASS 2.85f

#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>

using namespace std;


//Z: Wheat Internode Data Structure
struct WheatInternode {

	//Z ----------- INDIVIDUAL INTERNODE ---------------------------------------

	//Z Position
	int InRank;			//Z [input] number of the internode at that tiller (mainstem) counted from root to top, start with 1
	int TlOrder;		//Z [input] the order of tiller that this internode grows on, 0 means mainstem, 1 means the tillers from mainstem
	int TlRank;			//Z [input] the cumuRank of the tiller, sum of leaf rank (from 1) where tiller emerge, but mainstem tiller rank is 0
	bool mainstem;		//Z [input] if this internode is from the mainstem

	//Z Growing Stage
	bool wait_for_initiation;	//Z this internode is waiting for initiation
	bool initiation;	
	bool emerge;
	bool ligulation_preL;		//Z roughly say, that means previous leaf is matured or fully expanded, different from leaf (previous leaf emergence)
	bool ligulation_preL_done;
	bool stop_elongation;		//Z this internode is stop elongation because the internode already achieve the max length 
	bool mature;
	bool senesce;
	bool dead;
	bool force_to_death_current_step;
	bool this_internode_elongation;

	//Z Growing time/gdd measure
	float physAge;			//Z gdd time (in reference to endGrowth and lifeSpan, days)
	float InPseudoAge;		//Z day not gdd, time elasped at equivalent temperature, or say "internode pseudo age"
	float InPseudoAge_prev;	//Z InPseudoAge @ previous time step, i.e., "InPseudoAge" @ the previous time step
	float emergeAge;		//Z gdd internode emerge age
	float emergePseudoAge;	//Z internode expansion time based age (use to compute the Beta function)
	float elongAge;			//Z internode elongation age
	float seneAge;			//Z gdd since the end of activeAge, period for senescence
	float senePseudoAge;	//Z leaf senescene pseudo age scale
	float seneDuration;		//Z total time (not gdd) for senesces

	//Z environemtal data
	float N_effect;
	float slw_effect;

	//Z Hidden Part Growth
	float InWholeEmergeDist;    //Z [input] a tricky parameter, "the ligulation height of the latest emerged leaf (for that tiller) - the cumu internode heights below this internode"
								//	that is the whole length for the leaf to emerge, computer in the tiller class
	float InLength_Hidden;		//Z hidden zone length
	float Dist2Emerge;			//Z distance for the current leaf to emerge cm
	float RERmax_I_fit[6] = { 0.20736f,0.186624f,0.15552f,0.165024f,0.16416f,0.152064f };	//Z convert to day-1, Optimal RERmax (s-1 at 12°C) allowing to simulate internode dimensions
																							//  RERmax_dict_IN = {3: 2.4E-06, 4: 2.4E-06, 5: 2.4E-06, 6: 2.4E-06, 7: 2.16E-06, 8: 1.8E-06, 9: 1.91E-06, 10: 1.9E-06, 11: 1.76E-06}  #: s-1 at 12°C FIT jan 20
	float RERmax_I;																//Z RERmax_I_fit for that internode
	float Inconce_Sucrose_default = 5900.0f;									//Z the default leaf sucrose concentration "sucrose/mstruct", based on elong wheat example in wheat fspm github
	float Inconce_Aminoacids_default[4] = { 99000.0f, 193000.0f, 377000.0f };	//Z the default leaf N concentration "amino acids/Nstruct"
	float Inconce_Sucrose;
	float Inconce_Sucrose_mean;													//Z moving average of the C concentration
	float Inconce_N;

	//Z Internode Length (cm)
	float InLength;					//Z internode length cm
	float maxInLength_lig;			//Z max internode length computed at the time of ligulation (previous leaf)
	float maxInLength_ref;			//Z reference max internode, from ryesim that determine internode growth after elongation
	float maxInLength;				//Z max internode length
	float maxElongationLength;		//Z the length for elongation, the max internode length - the internode length already exists
	float maxElongationRate;		//Z max elongation rate during elongation stage
	float ptnInLengthIncrease;		//Z potnetial internode enlongation, cm
	float ptnInLengthDecrease;		//Z potential internode decrease due to death, cm
	float greenInLength;            //Z green or living internode length, after elongation, cm
	float matureInLength;			//Z matured internode length, after elongation, cm
	float seneInLength;				//Z senescent internode length, cm

	//Z Internode biomass (g)
	float lsiwfit[6] = { 0.028f,0.023f,0.017f,0.016f,0.014f,0.007f };
									//Z internode specific weight g cm^-1, first 8 internodes use the first value
	float ratio_LSIW_LSSW;			//Z ratio lineic structural internode mass / lineic structural sheath mass  of the specific structural dry masses(from data of J.Bertheloot, 2004)
	float LSIW;						//Z lineic structural weight of internode
	float InMass;					//Z internode biomass in g
	float InMassIncrease;			//Z biomass assigned to this internode in g
	float ptnInMassIncrease;		//Z potential internode biomass in g
	float greenInMass;				//Z green internode mass in g
	float deadInMass;				//Z dead internode mass in g

	//Z Internode Nitrogen (mg)
	float InNitrogenContent;		//Z mass fraction in % (e.g., 3.0%)
	float InNitrogenMass;			//Z internode nitrogen mass mg
	float InNitrogenIncrease;		//Z nitrogen mass assigned to this internode in mg
	float ptnInNitrogenMassIncrease;//Z potential internode nitrogen mass in mg

	//Z Internode Nitrogen (mg): output after organ senecent or dropped
	float InNitrogenReleaseLowBdd;	//Z the min N mass content, below that, no N will be released even that organ is dead, assume to be 0.5%
	float InNitrogenReleaseMaxPtge;	//Z during plant organ dying, max percent of N releasing

	//Z Individual internode
	//   Release: single internode N output due to senecense and dropping, need to divide into remobile (return) portion and dead (with internode yellow) portion
	//   Dead: at that step, reduced N goes to residue and immobilized
	//   Return: at that step, remobilized and can reused for new plant organs
	//   DeadCumu: cumulated N mass in the senecense/dead plant organs
	float InNitrogenRelease;		//Z nitrogen release due to internode decreasing or dying in mg
	float InNitrogenDead;			//Z stepwise dead portion of intr nitrogen, will be in the dead plant tissue and not removable
	float InNitrogenReturn;			//Z stepwise returned portion of intr nitrogen, will supply plant future usage 
	float InNitrogenMassDeadCumu;	//Z cumulated dead portion of nitrogen, in the dead plant tissue and not removable

	//Z ----------- REPRESENTATIVE INTERNODE ---------------------------------------

	//Z living fraction
	//	trace the living fraction changes for one internode, it is living or dead;
	//	for a representative plant, need to consider statistical effects -- therefore, a internode can be partially living.
	//	when living fraction changes, caused by tiller living fraction, those portion should be considered as dead and dropped
	//	therefore, computation of droping should be done and N should be recycled
	float livingFrac;
	float livingFrac_old;
	float livingFrac_ini;

	float greenInLength_Rep;
	float seneInLength_Rep;
	float InLength_Rep;			//Z internode length include everything, = green + sene, and internode will not drop as internode

	float InMass_Rep;			//Z representative internode mass, including everything
	float greenInMass_Rep;		//Z green internode mass
	float deadInMass_Rep;		//Z dead internode mass, internode will not drop, so we end at internode death

	float InNitrogenMass_Rep;
	//Z Recall:
	//  Release (single internode): single internode N output due to senecense and dropping, need to divide
	//  Dead (representative internode, later): at that step, reduced N goes to residue and immobilized
	//  Return (representative internode, later): at that step, remobilized and can reused for new plant organs
	//  DeadCumu (representative internode, later): cumulated N mass in the senecense/dead plant organs
	//Z only return and dead N have correlation with plant level mass budget
	//  so we only need to report return (to plant) and deat (fall with dropped internode) N
	//  leaf needs "sene+drop=dead"
	//  internode does not need, since internode will not drop
	float InNitrogenReturn_Rep;
	float deadInNitrogenMass_Rep;

	float ptnInMassIncrease_Rep;
	float ptnInNitrogenMassIncrease_Rep;

	//Z ----------- WHEAT internode CONSTRUCTOR ---------------------------------------
	WheatInternode();
	WheatInternode(int rank, int order, int tlRank, bool mainstem, float livingFrac);

};

//Z Wheat Internode Operator

//Z Wheat Internode Update Operator
struct WheatInternodeUpdate {

	//Z input for this operator should be listed as a member variable
	//Z "WheatPlant" should call "SetInternodeCondition" function and update the values 
	float current_gdd;			//Z hourly gdd incremental
	float current_tmpr;			//Z hourly internode temperature, just air tempearture
	float current_teq;			//Z equivalent time (day) at reference temperature (12oC)
	float predawn_psi;			//Z predawn water potential in Mpa
	float shade_effect;			//Z hourly shade effect
	bool plant_elongation;		//Z the whole plant enters elongation, so internode can grow up, this is plant level growth stage
	bool plant_antepa;			//Z the end of anthesis, where assume the stems start to be yellow, this is plant level growth stage

	//Z "WheatInternode" should compute the following effects before executing the update pipeline
	//  N effect / stl effect on internodes will vary among internodes, so will be local variables
	float psi_effect_elong;
	float psi_effect_sene;
	float tmpr_effect;

	//Z Internode charateristic temperatures
	const float Tmax_Intr = 30.0f;
	const float Topt_Intr = 15.0f;
	const float Tmin_Intr = 0.0f;

	//Z funtional operators
	WheatInternodeUpdate() : current_gdd(0.0f), current_tmpr(0.0f), current_teq(0.0f), predawn_psi(0.0f), shade_effect(1.0f),
		psi_effect_elong(1.0f), psi_effect_sene(1.0f), tmpr_effect(1.0f), plant_elongation(false), plant_antepa(false) {}
	void SetInternodeUpdateCondition(float gdd, float tmpr, float teq, float psi_predawn, float shaded, bool elong, bool antepa);
	//Z Input tuple for this function includes
	//		1. Internode structure;
	//      2. current living fraction (from tiller);
	//      3. force_to_death by the plant;
	//		4. if previous leaf is emerged
	//		5. a tricky parameter, "the ligulation height of the latest emerged leaf (for that tiller) - the internode heights below this leaf"
	//						that is the whole length for the leaf to emerge, computer in the tiller class
	//      6. this is the elongating internode;
	//void operator() (std::tuple<WheatInternode&, float, bool, bool> t)
	void IntrUpdateRun(WheatInternode& wIntr, float lf, bool f2d, bool ligulation_preL, float ligulationDist, bool elong) const;
};

//Z Wheat Internode Mass Distribution Operator
struct WheatInternodeMassAssignment {

	//Z input for this operator should be listed as a member variable
	// 	the variable listed here should be universal among internodes, under the whole plant scale.

	//Z "WheatPlant" should call "SetInternodeMassAssignmentCondition" function and update the values 
	float biomassIncomeRate;	//Z biomass and N increase rates from the plant level mass budget
	float nitrogenIncomeRate;

	//Z funtional operators
	WheatInternodeMassAssignment() : biomassIncomeRate(0.0f), nitrogenIncomeRate(0.0f) {}
	void SetInternodeMassAssignmentCondition(float biomassRate, float nitrogenRate);
	void IntrMassRun(WheatInternode& wIntr) const;
};

//**********************************************
//Z Lambda function for transformation and reduced sum
//  but we encapsulate the lambda function within structures
struct InExtRep_InLength { float operator() (const WheatInternode& wi) const { return wi.InLength_Rep; } };
struct InExtRep_InMass { float operator() (const WheatInternode& wi) const { return wi.InMass_Rep; } };
struct InExtRep_InNitrogenMass { float operator() (const WheatInternode& wi) const { return wi.InNitrogenMass_Rep; } };
struct InExtRep_InNitrogenReturn { float operator() (const WheatInternode& wi) const { return wi.InNitrogenReturn_Rep; } };
struct InExtRep_deadInNitrogenMass { float operator() (const WheatInternode& wi) const { return wi.deadInNitrogenMass_Rep; } };
struct InExtRep_ptnInMassIncrease { float operator() (const WheatInternode& wi) const { return wi.ptnInMassIncrease_Rep; } };
struct InExtRep_ptnInNitrogenMassIncrease { float operator() (const WheatInternode& wi) const { return wi.ptnInNitrogenMassIncrease_Rep; } };


#endif