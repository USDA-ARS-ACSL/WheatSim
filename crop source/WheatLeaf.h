#pragma once
#ifndef _WHEAT_LEAF_
#define _WHEAT_LEAF_

#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.000694f
#define PHYLLOCHRON 106.0f
#define HALFPHYLLOCHRON 53.0f
#define MAX_N_PCT 3.5f
#define MIN_N_PCT 0.5f

//Z leaf shape parameters from fspm
//Z time (day) difference for the next primordia
#define PLASTOCHRONE  10.34167f	//Z day (= 76.1 / 12 * 24 * 3600 sec) Leaf plastochron(sec at 12°c but casted to "day" in this study) calculated from Ljutovac 2002 with primordia of 5E-5 m(76 dd); Malvoisin 35dd associated with init 3E-5 m
//Z max leaf length
#define LEAF_LMAX_ABS 70.0f
//Z  parameter for automate elongation, presented in "s" but cast it to "day" here
#define LEAF_TE 35.0f//25.0f		//Z end of leaf elongation in automate growth (sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define LEAF_TM 17.05f		//Z time at which leaf elongation rate is maximal in automate growth(sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define LEAF_TB -9.525f		//Z beginning of leaf elongation in automate growth(sec at 12°c but casted to "day" in this study); fitted from adapted data from Fournier 2005
#define RER_KC	100.0f		//Z affinity coefficient of RER to C (µmol g-1)
#define RER_KN	15.0f		//Z affinity coefficient of RER to N (µmol g-1)
//Z sheat leaf ratio paramters
#define SL_RATIO_A	-0.0021f
#define SL_RATIO_B  0.037f
#define SL_RATIO_C  -0.1527f
#define SL_RATIO_D	0.4962f
#define SL_RATIO_MIN  0.35f
//Z leaf length width ration parameter #: m (Ljutovac 2002)
#define LEAF_WL_BASE	0.05f
#define LEAF_WL_REGUL_MIN 0.5f
#define LEAF_WL_REGUL_MAX 2.0f
#define LEAF_WL_INT_MIN	200.0f
#define LEAF_WL_INT_MAX 5000.0f
//Z leaf specific mass parameter
#define LF_SSLW_MIN	0.0005f				//Z g/cm^2
#define LF_SSLW_MAX 0.0045f
#define LF_SSLW_INT_MIN	400.0f			//Z these two are just coefficients
#define LF_SSLW_INT_MAX 5200.0f
//Z sheath specific mass parameter
#define SH_LSSW_MIN	0.0005f				//Z g/cm
#define SH_LSSW_MAX 0.0080f
#define SH_LSSW_A	0.00005f			//Z these two are just coefficients
#define SH_LSSW_INT_MIN 1700.0f
#define SH_LSSW_NOMINAL_A 0.0403f
#define SH_LSSW_NOMINAL_B -0.0099f
//Z leaf mass growth parameter
#define LF_ALPHA 0.106f					//Z Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
#define LF_BETA 1.28f					//Z Parameter of the relation between leaf mass and leaf length (dimensionless)
#define LF_DRYMASS_2_STRCTMASS 0.35f	//Z converter between structure mass to dry mass, others can be functional materials
#define LF_STRCTMASS_2_DRYMASS 2.85f

#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>

using namespace std;

//Z: Wheat Leaf Data Structure
struct WheatLeaf {

	//Z ----------- INDIVIDUAL LEAF ---------------------------------------

	//Z Position, leaf rank starts from "1" such that the higher order index can be automatically cumulated
	//            if start from "0", then the TlRank (cumulative value) will be constantly 0 or 1 
	int LfRank;					//Z [input] number of the leaf at that tiller (mainstem) counted from root to top, start with 0
	int TlOrder;				//Z [input] the order of tiller that this leaf grows on, 0 means mainstem, 1 means the tillers from mainstem
	int TlRank;					//Z [input] the cumuRank of the tiller, sum of leaf rank (from 1) where tiller emerge, but mainstem tiller rank is 0
	bool mainstem;				//Z [input] if this leaf is from the mainstem

	//Z Growing Stage
	bool firstupdate;
	bool initiation;			//Z new leaf initiate
	bool emerge;				//Z new leaf enters "expand" stage automatically (after emergence)
	bool emerge_preL;			//Z previous leaf emerge, some operation on the current leaf, e.g., compute "leaf length max em"
	bool emerge_preL_done;		//Z previous leaf emerge_done
	bool photosyn_Start;		//Z suppose to the time when leaf length >= 1/2 max leaf length,  
	bool mature;				//Z automate elongation stops
	bool senesce;				//Z sene starts
	bool dead;					//Z sene finished
	bool force_to_death_current_step;	//Z the parent tiller dies, then all the leaf will be forced to die

	//Z Growing time/gdd measure
	float physAge;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	float LfPseudoAge;          //Z day not gdd, time elasped at equivalent temperature, or say "leaf pseudo age"
	float LfPseudoAge_prev;		//Z LfPseudoAge @ previous time step
	float emergeAge;			//Z gdd leaf expansion age
	float emergePseudoAge;		//Z leaf expansion time based age (use to compute the Beta function), equivalent time until leaf mature
	float mature_TTd;			//Z gdd at mature, after fully expansion
	float activeAge;			//Z gdd since mature, the total functioning time period
	float seneAge;				//Z gdd since the end of activeAge, period for senescence
	float senePseudoAge;		//Z leaf senescene pseudo age scale
	float seneDuration; 		//Z total equivalent time for senesces

	//Z environemtal data
	float N_effect;
	float slw_effect;

	//Z Hidden Part Growth
	float LfWholeEmergeDist;    //Z [input] a tricky parameter, "the ligulation height of the latest matured leaf (for that tiller) - the internode heights below this current leaf"
								//	that is the whole length for the leaf to emerge, computer in the tiller class
	float LfPseudoStem;			//Z part of the sheath, equal to "LfWholeEmergeDist" before emergence, but will decrease after emergence because the internode will grow.
								//  before emergence, all the leaf components are "pseudoStem"
	float Dist2Emerge;			//Z distance for the current leaf to emerge cm
	float RERmax_L_fit[7] = { 0.2592f,0.1512f,0.141696f,0.133056f, 0.130464f,0.115776f,0.111456f };	//Z convert to day-1, Optimal RERmax (s-1 at 12°C) allowing to simulate leaf dimensions of Ljutovac (2002) for leaf 5-11, 
																									//  { 0.000003, 0.00000175, 0.00000164, 0.00000154, 0.00000151, 0.00000134, 0.00000129 };
	float RERmax_L;																//Z RERmax_L_fit for that leaf
	float Lfconce_Sucrose_default = 5900.0f;									//Z the default leaf sucrose concentration "sucrose/mstruct", based on elong wheat example in wheat fspm github
	float Lfconce_Aminoacids_default[4] = { 99000.0f, 193000.0f, 377000.0f };	//Z the default leaf N concentration "amino acids/Nstruct"
	float Lfconce_Sucrose;
	float Lfconce_Sucrose_mean;													//Z moving average of the C concentration
	float Lfconce_N;
	float SL_ratio;				//Z sheath:lamina ratio
	
	//Z The previous leaf is emerged, now compute some leaf geometry
	float LfLength;				//Z leaf length in cm
	float maxLfLength;			//Z max leaf length based on leaf rank and if it is at the mainstem
	float maxLfLength_em;		//Z max leaf length computed @ prev leaf emergence
	float maxLfLength_ref;		//Z max leaf length reference @ using shootgro
	float maxLfWidth;			//Z max leaf width in cm
	float maxLfWidth_fit[9] = { 0.4f,0.45f,0.56f,0.75f,1.0f,1.2f,1.3f,1.4f,1.8f }; //Z cm max LaWidth
	float maxLfWidth_ref;		//Z reference leaf max width for that leaf position	

	float LaLength;				//Z Lamina length in cm, only happens after leaf is emerged
	float LaWidth;				//Z Lamina width in cm
	float maxLaLength;			//Z max lamina length in cm
	
	float shapeFrac;			//Z shape factor "area = shape factor * length * width", dimensionless
	float LaArea;				//Z total lamina area changing over time, cm^2
	float matureLaArea;			//Z lamina area at the end of expansion, i.e., the actual max lamina area at one time (different from LaminaArea)
	float seneLaArea;			//Z aged lamina area changing over time, cm^2
	float dropLaArea;			//Z drop lamina area, the lamina area at drop (dead) time point
	float greenLaArea;			//Z green (active) lamina area, i.e., lamina area before aging, and part of the lamina area after aging
	float maxLaArea;			//Z max lamina area based on leaf rank and if it is at the mainstem
	float ptnLaAreaIncrease;	//Z potential lamina area increase, the values after one step of potential growth
	float ptnLaAreaDecrease;	//Z potential lamina area decrease, the values after one step of aging

	float ShLength;				//Z sheath length
	float maxShLength;			//Z max sheath length in cm
	float ptnShLengthIncrease;	//Z sheath length increase at each time step

	//Z leaf biomass + sheath biomass (g)
	//	dropped mass is still mass, so even if the leaf is dropped, it should be counted for leafmass
	//	therefore, if needed, LeafMass-DropLfMass = attached leafmass
	//	(for single leaf, either attached leafmass 0 or = leafmass, since leaf is either attached or dropped)

	float sslwfit[11] = { 0.0015f, 0.0023f, 0.0025f, 0.0024f, 0.0021f, 0.0018f, 0.0016f, 0.0018f, 0.0021f, 0.0026f, 0.0033f };  //Z g/cm^-2 # Manip NEMA 05 / 06 traitments N + (from data of J.Bertheloot, 2004) sauf pour F7 / F8
	float lsswfit[11] = { 0.0008f, 0.0009f, 0.0011f, 0.0018f, 0.0017f, 0.0021f, 0.0024f, 0.0040f, 0.0050f, 0.0055f, 0.0065f };	//Z g/cm # Manip NEMA 05/06 Soissons N+ (from data of J. Bertheloot, 2004)
	float SSLW;					//Z structural specific lamina weight 
	float LSSW;					//Z lineic structural sheath weight
	float LaMass;				//Z lamina biomass no matter it is green or not, dropped or not
	float ShMass;				//Z sheath biomass no matter it is green or not, dropped or not
	float greenLaMass;			//Z green (living) lamina mass, exclude senescent (and dropped) mass
	float greenShMass;			//Z green (living) sheath mass, exclude senescent (and dropped) mass
	float dropLaMass;			//Z dropped lamina biomass, 0 or = lamina mass depending on lamina drop or not
	float dropShMass;			//Z dropped sheath biomass, 0 or = sheathmass depending on lamina drop or not
	float LaBiomassIncrease;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	float ShBiomassIncrease;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	float ptnLaMassIncrease;	//Z cumulated value, potential lamina mass increase due to lamina area increase at this time step and biomass is not fully assigned in previous steps
	float ptnShMassIncrease;	//Z cumulated value, potential sheath mass increase due to lamina area increase at this time step and biomass is not fully assigned in previous steps


//Z Leaf Sheath Nitrogen (mg): income part
	float LaNitrogenMass;				//Z lamina nitorgen mass
	float ShNitrogenMass;				//Z sheath nitrogen mass
	float LaNitrogenContent;			//Z lamina N% fraction in % (e.g., 3.5% as the max value)
	float ShNitrogenContent;			//Z sheath N% fraction in % (e.g., 3.5% as the max value)
	float LaNitrogenIncrease;			//Z lamina nitrogen mass increase after nitrogen assignment 
	float ShNitrogenIncrease;			//Z sheath nitrogen mass increase after nitrogen assignment 
	float ptnLaNitrogenMassIncrease;	//Z cumulated value, potential lamina N mass increase due to lamina area increase at this time step, and N is not fully assigned in previous steps
	float ptnShNitrogenMassIncrease;	//Z cumulated value, potential sheath N mass increase due to lamina area increase at this time step, and N is not fully assigned in previous steps


//Z Lamina Sheath Nitrogen (mg): output after organ senecent or dropped
	float NitrogenReleaseLowBdd;		//Z the min lamina N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	float NitrogenReleaseMaxPtge;		//Z the lamina N release percentage, default to be 80%, that means 80% of lamina N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

//Z Individual lamina
//  Release: single lamina N output due to senecense and dropping, need to divide into remobile (released) portion and dead (with lamina drop) portion
//  Dead: at that step, reduced N goes to residue and immobilized
//  Return: at that step, remobilized and can reused for new plant organs
//  DeadCumu: cumulated N mass in the senecense/dead plant organs

	float LaNitrogenRelease;
	float ShNitrogenRelease;
	float LaNitrogenDead;				//Z stepwise dead portion of lamina nitrogen, will be in the dead plant tissue and not removable
	float ShNitrogenDead;				//Z stepwise dead portion of sheath nitrogen, will be in the dead plant tissue and not removable
	float LaNitrogenReturn;				//Z stepwise returned portion of lamina nitrogen, will supply plant future usage 
	float ShNitrogenReturn;				//Z stepwise returned portion of sheath nitrogen, will supply plant future usage 
	float LaNitrogenMassDeadCumu;		//Z cumulative dead portion of lamina nitrogen, will be in the dead plant tissue and not removable
	float ShNitrogenMassDeadCumu;		//Z cumulative dead portion of sheath nitrogen, will be in the dead plant tissue and not removable


	//Z ----------- REPRESENTATIVE LEAF ---------------------------------------

		//Z Leaf geometry and mass properties after considering "livingFrac",
		//	this will make this plant and leaf become "representative plant" or "representative leaf" from the statistical prespective
		//	individual leaf growth and shape is computated based on "one single leaf"
		//	while we use "livingFrac" to adjust the "one single leaf" to be a "representative leaf"
		//	we use "_Rep" to indicate "Representative".

		//Z why we adjust "livingfrac" based on leaf, a very deep / elementry level class ?
		//	that is because our model design that leaf handles its own shape and growth STRICTLY within the leaf class


		//Z living fraction
		//	trace the living fraction changes for one leaf, it is living or dead;
		//	for a representative plant, need to consider statistical effects -- therefore, a leaf can be partially living.
		//	when living fraction changes, caused by tiller living fraction, those portion should be considered as dead and dropped
		//	therefore, computation of droping should be done and N should be recycled
	float livingFrac;
	float livingFrac_old;
	float livingFrac_ini;

	float greenLaArea_Rep;
	float seneLaArea_Rep;
	float dropLaArea_Rep;		//Z dropped area = 0 and changed to sene area when the lamina is totally dead
	float LaArea_Rep;			//Z lamina area include everything, lamina area = green + sene, 
	

	float LaMass_Rep;
	float ShMass_Rep;
	float greenLaMass_Rep;
	float greenShMass_Rep;
	float dropLaMass_Rep;
	float dropShMass_Rep;

	float LaNitrogenMass_Rep;
	float ShNitrogenMass_Rep;

	//Z  Recall:
	//   Release (single leaf): single leaf N output due to senecense and dropping, need to divide
	//   Dead (representative leaf, later): at that step, reduced N goes to residue and immobilized
	//   Return (representative leaf, later): at that step, remobilized and can reused for new plant organs
	//   DeadCumu (representative leaf, later): cumulated N mass in the senecense/dead plant organs
	//Z  only return and dead N have correlation with plant level mass budget
	//   so we only need to report return (to plant) and deat (fall with dropped leaf) N
	float LaNitrogenReturn_Rep;
	float ShNitrogenReturn_Rep;
	float deadLaNitrogenMass_Rep;
	float deadShNitrogenMass_Rep;

	//Z for dead leaf/sheath nitrogen mass rep, we further partition that into
	//  "sene": N in the senescent portion but not dropped from the stem
	//  "drop": N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     Dead = sene(senescent) + drop

	float seneLaNitrogenMass_Rep;
	float seneShNitrogenMass_Rep;
	float dropLaNitrogenMass_Rep;
	float dropShNitrogenMass_Rep;

	float ptnLaMassIncrease_Rep;
	float ptnShMassIncrease_Rep;
	float ptnLaNitrogenMassIncrease_Rep;
	float ptnShNitrogenMassIncrease_Rep;

	//Z ----------- WHEAT LEAF CONSTRUCTOR ---------------------------------------
	WheatLeaf();
	WheatLeaf(int rank, int order, int tlRank, bool mainstem, float livingFrac);

};

//Z Wheat Leaf Operator
//Z Wheat Leaf Update Operator
struct WheatLeafUpdate {

	//Z input for this operator should be listed as a member variable
	//Z "WheatPlant" should call "SetLeafCondition" function and update the values 
	float current_gdd;			//Z hourly gdd incremental
	float current_tmpr;			//Z hourly leaf temperature
	float current_teq;			//Z equivalent time (day) at reference temperature (12oC)
	float predawn_psi;			//Z predawn leaf water potential in Mpa
	float shade_effect;			//Z hourly shade effect
	bool acceleration;			//Z the whole plant enters acceleration sheath growth or not (single ridge)

	//Z "WheatLeaf" should compute the following effects before executing the update pipeline
	//  N effect / slw effect on leaves will vary among leaves, so will be local variables
	float psi_effect_expand;
	float psi_effect_sene;
	float tmpr_effect;

	//Z Leaf charateristic temperatures
	const float Tmax_Leaf = 30.0f;
	const float Topt_Leaf = 15.0f;
	const float Tmin_Leaf = 0.0f;

	//Z funtional operators
	WheatLeafUpdate() : current_gdd(0.0f), current_tmpr(0.0f), current_teq(0.0f), predawn_psi(0.0f), shade_effect(1.0f),
		psi_effect_expand(1.0f), psi_effect_sene(1.0f), tmpr_effect(1.0f), acceleration(false) {}
	void SetLeafUpdateCondition(float gdd, float tmpr, float teq, float psi_predawn, float shaded, bool accel);
	//Z parameters : 1. Leaf structure;
	//               2. current living fraction (from tiller);
	//               3. force_to_death by the plant;
	//				 4. if previous leaf is emerged
	//				 5. a tricky parameter, "the ligulation height of the latest emerged leaf (for that tiller) - the internode heights below this leaf"
	//					that is the whole length for the leaf to emerge, computer in the tiller class
	//void operator() (std::tuple<WheatLeaf&, float, bool> t)
	void LeafUpdateRun(WheatLeaf& wLeaf, float lf, bool f2d, bool emerge_preL, float ligulationDist) const;
};

//Z Wheat Leaf Mass Distribution Operator
struct WheatLeafMassAssignment {

	//Z input for this operator should be listed as a member variable
	// 	the variable listed here should be universal among leaves, under the whole plant scale.

	//Z "WheatPlant" should call "SetLeafMassAssignmentCondition" function and update the values 
	bool acceleration;			//Z the whole plant enters acceleration sheath growth or not (single ridge)
	float biomassIncomeRate;	//Z biomass and N increase rates from the plant level mass budget
	float nitrogenIncomeRate;

	//Z funtional operators
	WheatLeafMassAssignment() : acceleration(false), biomassIncomeRate(0.0f), nitrogenIncomeRate(0.0f) {}
	void SetLeafMassAssignmentCondition(bool accel, float biomassRate, float nitrogenRate);
	void LeafMassRun(WheatLeaf& wLeaf) const;
};

//**********************************************
//Z Lambda function for transformation and reduced sum
//  but we encapsulate the lambda function within structures
struct LfExtRep_LfNum { float operator() (const WheatLeaf& wl) const { return wl.livingFrac_ini; } };
struct LfExtRep_greenLfNum { float operator() (const WheatLeaf& wl) const { return wl.livingFrac; } };
struct LfExtRep_LaArea { float operator() (const WheatLeaf& wl) const { return wl.LaArea_Rep; } };
struct LfExtRep_greenLaArea { float operator() (const WheatLeaf& wl) const { return wl.greenLaArea_Rep; } };
struct LfExtRep_seneLaArea { float operator() (const WheatLeaf& wl) const { return wl.seneLaArea_Rep; } };
struct LfExtRep_dropLaArea { float operator() (const WheatLeaf& wl) const { return wl.dropLaArea_Rep; } };
struct LfExtRep_LaMass { float operator() (const WheatLeaf& wl) const { return wl.LaMass_Rep; } };
struct LfExtRep_ShMass { float operator() (const WheatLeaf& wl) const { return wl.ShMass_Rep; } };
struct LfExtRep_dropLaMass { float operator() (const WheatLeaf& wl) const { return wl.dropLaMass_Rep; } };
struct LfExtRep_dropShMass { float operator() (const WheatLeaf& wl) const { return wl.dropShMass_Rep; } };
struct LfExtRep_LaNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.LaNitrogenMass_Rep; } };
struct LfExtRep_ShNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.ShNitrogenMass_Rep; } };
struct LfExtRep_LaNitrogenReturn { float operator() (const WheatLeaf& wl) const { return wl.LaNitrogenReturn_Rep; } };
struct LfExtRep_ShNitrogenReturn { float operator() (const WheatLeaf& wl) const { return wl.ShNitrogenReturn_Rep; } };
struct LfExtRep_deadLaNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.deadLaNitrogenMass_Rep; } };
struct LfExtRep_deadShNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.deadShNitrogenMass_Rep; } };
struct LfExtRep_seneLaNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.seneLaNitrogenMass_Rep; } };
struct LfExtRep_seneShNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.seneShNitrogenMass_Rep; } };
struct LfExtRep_dropLaNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.dropLaNitrogenMass_Rep; } };
struct LfExtRep_dropShNitrogenMass { float operator() (const WheatLeaf& wl) const { return wl.dropShNitrogenMass_Rep; } };
struct LfExtRep_ptnLaMassIncrease { float operator() (const WheatLeaf& wl) const { return wl.ptnLaMassIncrease_Rep; } };
struct LfExtRep_ptnShMassIncrease { float operator() (const WheatLeaf& wl) const { return wl.ptnShMassIncrease_Rep; } };
struct LfExtRep_ptnLaNitrogenMassIncrease { float operator() (const WheatLeaf& wl) const { return wl.ptnLaNitrogenMassIncrease_Rep; } };
struct LfExtRep_ptnShNitrogenMassIncrease { float operator() (const WheatLeaf& wl) const { return wl.ptnShNitrogenMassIncrease_Rep; } };

#endif