#include "WheatLeaf.h"

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Leaf Data Structure

WheatLeaf::WheatLeaf() :
	LfRank(1), TlOrder(0), TlRank(0), mainstem(false), livingFrac(1.0f), livingFrac_old(1.0f), livingFrac_ini(1.0f)
{
	//Z default initialization
	firstupdate = true;
	initiation = true;
	emerge = false;
	emerge_preL = false;
	emerge_preL_done = false;
	photosyn_Start = false;
	mature = false;
	senesce = false;
	dead = false;
	force_to_death_current_step = false;

	//Z Growing time/gdd measure
	physAge = 0.0f;				//Z gdd time
	LfPseudoAge = 0.0f;			//Z day not gdd, time elasped at equivalent temperature, or say "leaf pseudo age"
	LfPseudoAge_prev = 0.0f;	//Z LfPseudoAge @ previous time step, i.e., "LfPseudoAge" @ the previous time step
	emergeAge = 0.0f;			//Z gdd leaf expansion age
	emergePseudoAge = 0.0f;		//Z leaf expansion time based age (use to compute the Beta function)
	mature_TTd = 0.0f;			//Z gdd at mature, after fully expansion
	activeAge = 0.0f;			//Z gdd since mature, the total functioning time period
	seneAge = 0.0f;				//Z gdd since the end of activeAge, period for senescence
	senePseudoAge = 0.0f;		//Z leaf senescene pseudo age scale
	seneDuration = 0.0f;		//Z total equivalent time for senesces

	//Z environemtal data
	N_effect = 1.0f;
	slw_effect = 1.0f;

	//Z Hidden Part Growth
	LfWholeEmergeDist = 2.94f;
	LfPseudoStem = 0.005f;
	Dist2Emerge = 2.89f;
	RERmax_L = static_cast<float>(LfRank <= 5) * RERmax_L_fit[0]
		+ static_cast<float>(LfRank == 6) * RERmax_L_fit[1]
		+ static_cast<float>(LfRank == 7) * RERmax_L_fit[2]
		+ static_cast<float>(LfRank == 8) * RERmax_L_fit[3]
		+ static_cast<float>(LfRank == 9) * RERmax_L_fit[4]
		+ static_cast<float>(LfRank == 10) * RERmax_L_fit[5]
		+ static_cast<float>(LfRank >= 11) * RERmax_L_fit[6];
	Lfconce_Sucrose = Lfconce_Sucrose_default;
	Lfconce_N = static_cast<float>(LfRank == 1 && TlOrder == 0) * Lfconce_Aminoacids_default[0]
		+ static_cast<float>((LfRank + TlRank) == 2) * Lfconce_Aminoacids_default[1]
		+ static_cast<float>((LfRank + TlRank) >= 2) * Lfconce_Aminoacids_default[2];
	Lfconce_Sucrose_mean = 2100.0f;

	SL_ratio = max(SL_RATIO_A * powf(static_cast<float>(LfRank), 3.0f)
		+ SL_RATIO_B * powf(static_cast<float>(LfRank), 2.0f)
		+ SL_RATIO_C * static_cast<float>(LfRank)
		+ SL_RATIO_D, SL_RATIO_MIN);

	//Z The previous leaf is emerged, now compute some leaf geometry
	LfLength = 0.005f;

	maxLfLength_ref= (0.9f * static_cast<float>(LfRank) + 8.0f) * static_cast<float>(LfRank < 5)
		+ (4.5f * static_cast<float>(LfRank - 5) + 12.5f) * static_cast<float>(LfRank >= 5 && LfRank < 12)
		+ (30.5f - 5.0f * static_cast<float>(LfRank - 9)) * static_cast<float>(LfRank >= 12);
	maxLfLength_ref = min(min(maxLfLength_ref, 30.5f) * (1.0f + SL_ratio), LEAF_LMAX_ABS);

	maxLfLength = maxLfLength_ref;
	maxLfLength_em = maxLfLength_ref;

	maxLfWidth_ref = static_cast<float>(LfRank <= 3) * maxLfWidth_fit[0]
		+ static_cast<float>(LfRank == 4) * maxLfWidth_fit[1]
		+ static_cast<float>(LfRank == 5) * maxLfWidth_fit[2]
		+ static_cast<float>(LfRank == 6) * maxLfWidth_fit[3]
		+ static_cast<float>(LfRank == 7) * maxLfWidth_fit[4]
		+ static_cast<float>(LfRank == 8) * maxLfWidth_fit[5]
		+ static_cast<float>(LfRank == 9) * maxLfWidth_fit[6]
		+ static_cast<float>(LfRank == 10) * maxLfWidth_fit[7]
		+ static_cast<float>(LfRank >= 11) * maxLfWidth_fit[8];

	//Z lamina length and width based on the whole leaf geometry
	LaLength = 0.0f;
	maxLaLength = 0.0f;
	LaWidth = 0.0f;
	
	//Z compute max (potential) leaf area
	//  leaf shape factor 0.74 (most of the leaves) or 0.83 (first leaf on the main tiller)
	shapeFrac = 0.74f + 0.09f * static_cast<float>(LfRank == 1 && TlOrder == 0);
	LaArea = 0.0f;				//Z total leaf area changing over time, cm^2
	matureLaArea = 0.0f;		//Z leaf area at the end of expansion, i.e., the actual max leaf area at one time (different from LeafArea)
	seneLaArea = 0.0f;			//Z aged leaf area changing over time, cm^2
	dropLaArea = 0.0f;			//Z drop leaf area, the leaf area at drop (dead) time point
	greenLaArea = 0.0f;			//Z green (active) leaf area, i.e., leaf area before aging, and part of the leaf area after aging
	ptnLaAreaIncrease = 0.0f;	//Z potential leaf area increase, the values after one step of potential growth
	ptnLaAreaDecrease = 0.0f;	//Z potential leaf area decrease, the values after one step of aging

	ShLength = 0.005f;
	maxShLength = 0.005f;
	ptnShLengthIncrease = 0.0f;

	// leaf biomass + sheath biomass (g)
	if (LfRank <= 11)
	{
		SSLW = sslwfit[LfRank - 1];
		LSSW = lsswfit[LfRank - 1];
	}
	else
	{
		SSLW = sslwfit[10];
		LSSW = lsswfit[10];
	}

	LaMass = 0.0f;				//Z leaf biomass no matter it is green or not, dropped or not
	ShMass = 0.0f;				//Z sheath biomass no matter it is green or not, dropped or not
	greenLaMass = 0.0f;			//Z green (living) leaf mass, exclude senescent (and dropped) mass
	greenShMass = 0.0f;			//Z green (living) sheath mass, exclude senescent (and dropped) mass
	dropLaMass = 0.0f;			//Z dropped leaf biomass, 0 or = leafmass depending on leaf drop or not
	dropShMass = 0.0f;			//Z dropped sheath biomass, 0 or = sheathmass depending on leaf drop or not
	LaBiomassIncrease = 0.0f;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	ShBiomassIncrease = 0.0f;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	//Z cumulated mass demands, new leaf area/sheath mass will be added,
	//  if certian amount of mass is added, then reduce the ptn value
	ptnLaMassIncrease = 0.0f;
	ptnShMassIncrease = 0.0f;

	//Z Leaf Sheath Nitrogen (mg): income part
	LaNitrogenMass = 0.0f;				//Z leaf nitorgen mass
	ShNitrogenMass = 0.0f;				//Z sheath nitrogen mass
	LaNitrogenContent = MAX_N_PCT;		//Z leaf N% fraction in % (e.g., 3.5% as the max value)
	ShNitrogenContent = MAX_N_PCT;		//Z sheath N% fraction in % (e.g., 3.5% as the max value)
	LaNitrogenIncrease = 0.0f;			//Z leaf nitrogen mass increase after nitrogen assignment 
	ShNitrogenIncrease = 0.0f;			//Z sheath nitrogen mass increase after nitrogen assignment 
	//Z cumulated mass demands, new leaf area/sheath mass will be added,
	//  if certian amount of mass is added, then reduce the ptn value
	ptnLaNitrogenMassIncrease = 0.0f;
	ptnShNitrogenMassIncrease = 0.0f;

	//Z Leaf Sheath Nitrogen (mg): output after organ senecent or dropped
	NitrogenReleaseLowBdd = 0.5f;		//Z the min leaf N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	NitrogenReleaseMaxPtge = 0.8f;		//Z the leaf N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	//Z Individual leaf
	//   Release: single leaf N output due to senecense and dropping, need to divide into remobile (released) portion and dead (with leaf drop) portion
	//   Dead: at that step, reduced N goes to residue and immobilized
	//   Return: at that step, remobilized and can reused for new plant organs
	//   DeadCumu: cumulated N mass in the senecense/dead plant organs
	LaNitrogenRelease = 0.0f;
	ShNitrogenRelease = 0.0f;
	LaNitrogenDead = 0.0f;				//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	ShNitrogenDead = 0.0f;				//Z stepwise dead portion of sheath nitrogen, will be in the dead plant tissue and not removable
	LaNitrogenReturn = 0.0f;			//Z stepwise returned portion of leaf nitrogen, will supply plant future usage 
	ShNitrogenReturn = 0.0f;			//Z stepwise returned portion of sheath nitrogen, will supply plant future usage 
	LaNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	ShNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of sheath nitrogen, will be in the dead plant tissue and not removable


	//Z ----------- REPRESENTATIVE LEAF ---------------------------------------
	//greenLfLength_Rep = 0.0f;
	//greenLfWidth_Rep = 0.0f;
	greenLaArea_Rep = 0.0f;
	seneLaArea_Rep = 0.0f;
	dropLaArea_Rep = 0.0f;
	LaArea_Rep = 0.0f;

	LaMass_Rep = 0.0f;
	ShMass_Rep = 0.0f;
	greenLaMass_Rep = 0.0f;
	greenShMass_Rep = 0.0f;
	dropLaMass_Rep = 0.0f;
	dropShMass_Rep = 0.0f;

	LaNitrogenMass_Rep = 0.0f;
	ShNitrogenMass_Rep = 0.0f;

	//Z Recall:
	//   Release (single leaf): single leaf N output due to senecense and dropping, need to divide
	//   Dead (representative leaf, later): at that step, reduced N goes to residue and immobilized
	//   Return (representative leaf, later): at that step, remobilized and can reused for new plant organs
	//   DeadCumu (representative leaf, later): cumulated N mass in the senecense/dead plant organs
	//Z only return and dead N have correlation with plant level mass budget
	//   so we only need to report return (to plant) and deat (fall with dropped leaf) N
	LaNitrogenReturn_Rep = 0.0f;
	ShNitrogenReturn_Rep = 0.0f;
	deadLaNitrogenMass_Rep = 0.0f;
	deadShNitrogenMass_Rep = 0.0f;

	//Z for dead leaf/sheath nitrogen mass rep, we further partition that into
	//  "sene": N in the senescent portion but not dropped from the stem
	//  "drop": N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     Dead = sene(senescent) + drop
	seneLaNitrogenMass_Rep = 0.0f;
	seneShNitrogenMass_Rep = 0.0f;
	dropLaNitrogenMass_Rep = 0.0f;
	dropShNitrogenMass_Rep = 0.0f;

	ptnLaMassIncrease_Rep = 0.0f;
	ptnShMassIncrease_Rep = 0.0f;
	ptnLaNitrogenMassIncrease_Rep = 0.0f;
	ptnShNitrogenMassIncrease_Rep = 0.0f;
}

WheatLeaf::WheatLeaf(int rank, int order, int tlrank, bool mainstem, float livingFrac) :
	LfRank(rank), TlOrder(order), TlRank(tlrank), mainstem(mainstem), livingFrac(livingFrac), livingFrac_old(livingFrac), livingFrac_ini(livingFrac)
{
	//Z default initialization
	firstupdate = true;
	initiation = true;
	emerge = false;
	emerge_preL = false;
	emerge_preL_done = false;
	photosyn_Start = false;
	mature = false;
	senesce = false;
	dead = false;
	force_to_death_current_step = false;

	//Z Growing time/gdd measure
	physAge = 0.0f;				//Z gdd time
	LfPseudoAge = 0.0f;			//Z day not gdd, time elasped at equivalent temperature, or say "leaf pseudo age"
	LfPseudoAge_prev = 0.0f;	//Z LfPseudoAge @ previous time step, i.e., "LfPseudoAge" @ the previous time step
	emergeAge = 0.0f;			//Z gdd leaf expansion age
	emergePseudoAge = 0.0f;		//Z leaf expansion time based age (use to compute the Beta function)
	mature_TTd = 0.0f;			//Z gdd at mature, after fully expansion
	activeAge = 0.0f;			//Z gdd since mature, the total functioning time period
	seneAge = 0.0f;				//Z gdd since the end of activeAge, period for senescence
	senePseudoAge = 0.0f;		//Z leaf senescene pseudo age scale
	seneDuration = 0.0f;		//Z total equivalent time for senesces

	//Z environemtal data
	N_effect = 1.0f;
	slw_effect = 1.0f;

	//Z Hidden Part Growth
	LfWholeEmergeDist = 2.94f;
	LfPseudoStem = 0.005f;
	Dist2Emerge = 2.89f;
	RERmax_L = static_cast<float>(LfRank <= 5) * RERmax_L_fit[0]
		+ static_cast<float>(LfRank == 6) * RERmax_L_fit[1]
		+ static_cast<float>(LfRank == 7) * RERmax_L_fit[2]
		+ static_cast<float>(LfRank == 8) * RERmax_L_fit[3]
		+ static_cast<float>(LfRank == 9) * RERmax_L_fit[4]
		+ static_cast<float>(LfRank == 10) * RERmax_L_fit[5]
		+ static_cast<float>(LfRank >= 11) * RERmax_L_fit[6];
	Lfconce_Sucrose = Lfconce_Sucrose_default;
	Lfconce_N = static_cast<float>(LfRank == 1 && TlOrder == 0) * Lfconce_Aminoacids_default[0]
		+ static_cast<float>((LfRank + TlRank) == 2) * Lfconce_Aminoacids_default[1]
		+ static_cast<float>((LfRank + TlRank) >= 2) * Lfconce_Aminoacids_default[2];
	Lfconce_Sucrose_mean = 2100.0f;

	SL_ratio = max(SL_RATIO_A * powf(static_cast<float>(LfRank), 3.0f) 
		+ SL_RATIO_B * powf(static_cast<float>(LfRank), 2.0f) 
		+ SL_RATIO_C * static_cast<float>(LfRank) 
		+ SL_RATIO_D, SL_RATIO_MIN);

	//Z The previous leaf is emerged, now compute some leaf geometry
	LfLength = 0.005f;
	maxLfLength_ref = (0.9f * static_cast<float>(LfRank) + 8.0f) * static_cast<float>(LfRank < 5)
		+ (4.5f * static_cast<float>(LfRank - 5) + 12.5f) * static_cast<float>(LfRank >= 5 && LfRank < 12)
		+ (30.5f - 5.0f * static_cast<float>(LfRank - 9)) * static_cast<float>(LfRank >= 12);
	maxLfLength_ref = min(min(maxLfLength_ref, 30.5f) * (1.0f + SL_ratio), LEAF_LMAX_ABS);

	maxLfLength = maxLfLength_ref;
	maxLfLength_em = maxLfLength_ref;

	maxLfWidth_ref = static_cast<float>(LfRank <= 3) * maxLfWidth_fit[0]
		+ static_cast<float>(LfRank == 4) * maxLfWidth_fit[1]
		+ static_cast<float>(LfRank == 5) * maxLfWidth_fit[2]
		+ static_cast<float>(LfRank == 6) * maxLfWidth_fit[3]
		+ static_cast<float>(LfRank == 7) * maxLfWidth_fit[4]
		+ static_cast<float>(LfRank == 8) * maxLfWidth_fit[5]
		+ static_cast<float>(LfRank == 9) * maxLfWidth_fit[6]
		+ static_cast<float>(LfRank == 10) * maxLfWidth_fit[7]
		+ static_cast<float>(LfRank >= 11) * maxLfWidth_fit[8];
	maxLfWidth = maxLfWidth_ref;

	//Z lamina length and width based on the whole leaf geometry
	LaLength = 0.0f;
	maxLaLength = 0.0f;
	LaWidth = 0.0f;
	
	//Z compute max (potential) leaf area
	//  leaf shape factor 0.74 (most of the leaves) or 0.83 (first leaf on the main tiller)
	shapeFrac = 0.74f + 0.09f * static_cast<float>(LfRank == 1 && TlOrder == 0);
	LaArea = 0.0f;				//Z total leaf area changing over time, cm^2
	matureLaArea = 0.0f;		//Z leaf area at the end of expansion, i.e., the actual max leaf area at one time (different from LeafArea)
	seneLaArea = 0.0f;			//Z aged leaf area changing over time, cm^2
	dropLaArea = 0.0f;			//Z drop leaf area, the leaf area at drop (dead) time point
	greenLaArea = 0.0f;			//Z green (active) leaf area, i.e., leaf area before aging, and part of the leaf area after aging
	ptnLaAreaIncrease = 0.0f;	//Z potential leaf area increase, the values after one step of potential growth
	ptnLaAreaDecrease = 0.0f;	//Z potential leaf area decrease, the values after one step of aging

	ShLength = 0.005f;
	maxShLength = 0.005f;
	ptnShLengthIncrease = 0.0f;

	// leaf biomass + sheath biomass (g)
	if (LfRank <= 11) 
	{ 
		SSLW = sslwfit[LfRank - 1]; 
		LSSW = lsswfit[LfRank - 1]; 
	}
	else 
	{ 
		SSLW = sslwfit[10]; 
		LSSW = lsswfit[10]; 
	}

	LaMass = 0.0f;				//Z leaf biomass no matter it is green or not, dropped or not
	ShMass = 0.0f;				//Z sheath biomass no matter it is green or not, dropped or not
	greenLaMass = 0.0f;			//Z green (living) leaf mass, exclude senescent (and dropped) mass
	greenShMass = 0.0f;			//Z green (living) sheath mass, exclude senescent (and dropped) mass
	dropLaMass = 0.0f;			//Z dropped leaf biomass, 0 or = leafmass depending on leaf drop or not
	dropShMass = 0.0f;			//Z dropped sheath biomass, 0 or = sheathmass depending on leaf drop or not
	LaBiomassIncrease = 0.0f;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	ShBiomassIncrease = 0.0f;	//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	//Z cumulated mass demands, new leaf area/sheath mass will be added,
	//  if certian amount of mass is added, then reduce the ptn value
	ptnLaMassIncrease = 0.0f;	
	ptnShMassIncrease = 0.0f;

//Z Leaf Sheath Nitrogen (mg): income part
	LaNitrogenMass = 0.0f;				//Z leaf nitorgen mass
	ShNitrogenMass = 0.0f;				//Z sheath nitrogen mass
	LaNitrogenContent = MAX_N_PCT;		//Z leaf N% fraction in % (e.g., 3.5% as the max value)
	ShNitrogenContent = MAX_N_PCT;		//Z sheath N% fraction in % (e.g., 3.5% as the max value)
	LaNitrogenIncrease = 0.0f;			//Z leaf nitrogen mass increase after nitrogen assignment 
	ShNitrogenIncrease = 0.0f;			//Z sheath nitrogen mass increase after nitrogen assignment 
	//Z cumulated mass demands, new leaf area/sheath mass will be added,
	//  if certian amount of mass is added, then reduce the ptn value
	ptnLaNitrogenMassIncrease = 0.0f;
	ptnShNitrogenMassIncrease = 0.0f;

//Z Leaf Sheath Nitrogen (mg): output after organ senecent or dropped
	NitrogenReleaseLowBdd = 0.5f;		//Z the min leaf N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	NitrogenReleaseMaxPtge = 0.8f;		//Z the leaf N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	//Z Individual leaf
	//   Release: single leaf N output due to senecense and dropping, need to divide into remobile (released) portion and dead (with leaf drop) portion
	//   Dead: at that step, reduced N goes to residue and immobilized
	//   Return: at that step, remobilized and can reused for new plant organs
	//   DeadCumu: cumulated N mass in the senecense/dead plant organs
	LaNitrogenRelease = 0.0f;
	ShNitrogenRelease = 0.0f;
	LaNitrogenDead = 0.0f;				//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	ShNitrogenDead = 0.0f;				//Z stepwise dead portion of sheath nitrogen, will be in the dead plant tissue and not removable
	LaNitrogenReturn = 0.0f;			//Z stepwise returned portion of leaf nitrogen, will supply plant future usage 
	ShNitrogenReturn = 0.0f;			//Z stepwise returned portion of sheath nitrogen, will supply plant future usage 
	LaNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	ShNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of sheath nitrogen, will be in the dead plant tissue and not removable


	//Z ----------- REPRESENTATIVE LEAF ---------------------------------------
	//greenLfLength_Rep = 0.0f;
	//greenLfWidth_Rep = 0.0f;
	greenLaArea_Rep = 0.0f;
	seneLaArea_Rep = 0.0f;
	dropLaArea_Rep = 0.0f;
	LaArea_Rep = 0.0f;

	LaMass_Rep = 0.0f;
	ShMass_Rep = 0.0f;
	greenLaMass_Rep = 0.0f;
	greenShMass_Rep = 0.0f;
	dropLaMass_Rep = 0.0f;
	dropShMass_Rep = 0.0f;

	LaNitrogenMass_Rep = 0.0f;
	ShNitrogenMass_Rep = 0.0f;

	//Z Recall:
	//   Release (single leaf): single leaf N output due to senecense and dropping, need to divide
	//   Dead (representative leaf, later): at that step, reduced N goes to residue and immobilized
	//   Return (representative leaf, later): at that step, remobilized and can reused for new plant organs
	//   DeadCumu (representative leaf, later): cumulated N mass in the senecense/dead plant organs
	//Z only return and dead N have correlation with plant level mass budget
	//   so we only need to report return (to plant) and deat (fall with dropped leaf) N
	LaNitrogenReturn_Rep = 0.0f;
	ShNitrogenReturn_Rep = 0.0f;
	deadLaNitrogenMass_Rep = 0.0f;
	deadShNitrogenMass_Rep = 0.0f;

	//Z for dead leaf/sheath nitrogen mass rep, we further partition that into
	//  "sene": N in the senescent portion but not dropped from the stem
	//  "drop": N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     Dead = sene(senescent) + drop
	seneLaNitrogenMass_Rep = 0.0f;
	seneShNitrogenMass_Rep = 0.0f;
	dropLaNitrogenMass_Rep = 0.0f;
	dropShNitrogenMass_Rep = 0.0f;

	ptnLaMassIncrease_Rep = 0.0f;
	ptnShMassIncrease_Rep = 0.0f;
	ptnLaNitrogenMassIncrease_Rep = 0.0f;
	ptnShNitrogenMassIncrease_Rep = 0.0f;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Leaf Update Operator Structure

//Z update leaf growth gdd and conditions for this current hourly step
//  leaf stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatLeafUpdate::SetLeafUpdateCondition(float gdd, float tmpr, float teq, float psi_predawn, float shaded, bool accel)
{
	current_gdd = gdd;
	current_tmpr = tmpr;
	current_teq = teq;
	predawn_psi = psi_predawn;
	shade_effect = shaded;
	acceleration = accel;

	//Z compute temperature effect
	//  tmpr effects on growth from K.Paff
	//float AlphaGrowth = logf(2.0f) / (logf((Tmax_Leaf - Tmin_Leaf) / (Topt_Leaf - Tmin_Leaf)));
	//if (current_tmpr<Tmin_Leaf || current_tmpr>Tmax_Leaf) { tmpr_effect = 0.1f; }
	//else {
	//	float temp1 = powf((current_tmpr - Tmin_Leaf), AlphaGrowth);
	//	float temp2 = powf((Topt_Leaf - Tmin_Leaf), AlphaGrowth);
	//	tmpr_effect = (2.0f * temp1 * temp2 - powf(temp1, 2.0f)) / powf(temp2, 2.0f);
	//}
	//tmpr_effect = min(max(tmpr_effect, 0.1f), 1.0f);

	//Z compute water potential effect for expanding and senescence
	//Z some constant first
	const float psi_f = -1.4251f; // -1.0, was -1.4251   later changed to -2.3;
	const float s_f = 0.4258f; // was 0.4258 0.5;
	//Z expanding
	const float psi_threshold_bars_expand = -0.8657f;	//Mpa
	psi_effect_expand = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_expand))));
	psi_effect_expand = min(max(psi_effect_expand, 0.1f), 1.0f);
	//z senescence
	const float psi_threshold_bars_senesc = -4.0f;		//Mpa
	psi_effect_sene = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_senesc))));
	psi_effect_sene = min(max(psi_effect_sene, 0.1f), 1.0f);
}


//Z this is the main function for plant leaf growth
//  every time updating, hourly gdd, tmpr, water potential, shade are universal
//                       N stress are local
//                       living fractor is the input variable for each leaf object
//Z parameters :  1. Leaf structure;
//                2. current living fraction (from tiller);
//                3. force_to_death by the plant;
//				  4. if previous leaf is emerged
//				  5. a tricky parameter, "the ligulation height of the latest emerged leaf (for that tiller) - the internode heights below this leaf"
//						that is the whole length for the leaf to emerge, computer in the tiller class
//void operator() (std::tuple<WheatLeaf&, float, bool> t)
void WheatLeafUpdate::LeafUpdateRun(WheatLeaf& wLeaf, float lf, bool f2d, bool emerge_preL, float ligulationDist)
const {
	//Z Simulation of inidividual leaf

	wLeaf.livingFrac = lf;
	wLeaf.force_to_death_current_step = f2d;
	wLeaf.emerge_preL = emerge_preL;
	wLeaf.LfWholeEmergeDist = ligulationDist;	//Z need to input every step, since the ligulation and internode grow with respect to time
	wLeaf.Dist2Emerge = wLeaf.LfWholeEmergeDist - wLeaf.LfPseudoStem;
	wLeaf.physAge += current_gdd;

	//Z reset some critical parameters
	//  do not reset "wLeaf.ptnLaMassIncrease; wLeaf.ptnShMassIncrease; wLeaf.ptnLaNitrogenMassIncrease; wLeaf.ptnShNitrogenMassIncrease;"
	//  that will be cumulative demand that carrys over the unrealized demand from the previous steps.
	wLeaf.ptnLaAreaIncrease = 0.0f;
	wLeaf.ptnShLengthIncrease = 0.0f;
	wLeaf.ptnLaAreaDecrease = 0.0f;

	wLeaf.LaNitrogenRelease = 0.0f;
	wLeaf.ShNitrogenRelease = 0.0f;
	wLeaf.LaNitrogenReturn = 0.0f;
	wLeaf.ShNitrogenReturn = 0.0f;
	wLeaf.LaNitrogenDead = 0.0f;
	wLeaf.ShNitrogenDead = 0.0f;

	//Z compute the N stress for each inidividual leaves, N in mass fraction (mg/g biomass)
	wLeaf.N_effect = 2.0f / (1.0f + expf(-2.9f * (max(MIN_N_PCT, wLeaf.LaNitrogenContent) - MIN_N_PCT))) - 1.0f;
	wLeaf.N_effect = max(min(wLeaf.N_effect, 1.0f), 0.1f);

	float q10fn = powf(2.0f, 0.1f * (current_tmpr - Topt_Leaf));

	//Step ZERO Force to death ---------------------------------------------
	//     this condition is imposed by the parent tiller, stronger than any other conditions
	if ((!wLeaf.dead) && wLeaf.force_to_death_current_step)		//Z prevent kill the leaf twice
	{
		//Z make the leaf senecence to the end, 100%
		if (wLeaf.senesce) { 
			wLeaf.ptnLaAreaDecrease = max(0.0f, (wLeaf.matureLaArea - wLeaf.seneLaArea)); 
		}
		else {
			//Z suppose green leaf area is non-zero quantity after maturity
			//  in leaf expansion, just use Leaf area
			wLeaf.ptnLaAreaDecrease = max(0.0f, wLeaf.LaArea);
		}
		wLeaf.seneLaArea += wLeaf.ptnLaAreaDecrease;
		wLeaf.greenLaArea = 0.0f;
		wLeaf.greenLaMass = 0.0f;
		wLeaf.greenShMass = 0.0f;

		wLeaf.ptnLaMassIncrease = 0.0f;
		wLeaf.ptnShMassIncrease = 0.0f;
		wLeaf.ptnLaNitrogenMassIncrease = 0.0f;
	    wLeaf.ptnShNitrogenMassIncrease = 0.0f;

		//Z DO NOT change "livingFrac" here
		//  DO NOT set "GreenLf/Width/Length, LeafNitrogenMass, SheathNitrogenMass" to 0
		//  That is because we need those values to calculate N return and convert greenlf to droplf in the "LivingFractionAdjustment" function
		//  The correctness can be viewed by computing N mass fraction from the output.
		//  The N mass% ranges 0.25% to 3.5%
		
		wLeaf.LaNitrogenRelease = wLeaf.LaNitrogenMass;
		wLeaf.ShNitrogenRelease = wLeaf.ShNitrogenMass;
		wLeaf.LaNitrogenMass = 0.0f;
		wLeaf.ShNitrogenMass = 0.0f;
		//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
		float LaNitrogenReturnFrac = min(max(wLeaf.LaNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.LaNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
		float LaNitrogenDeadFrac = 1.0f - LaNitrogenReturnFrac;
		float ShNitrogenReturnFrac = min(max(wLeaf.ShNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.ShNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
		float ShNitrogenDeadFrac = 1.0f - ShNitrogenReturnFrac;
		wLeaf.LaNitrogenReturn = wLeaf.LaNitrogenRelease * LaNitrogenReturnFrac;
		wLeaf.ShNitrogenReturn = wLeaf.ShNitrogenRelease * ShNitrogenReturnFrac;
		wLeaf.LaNitrogenDead = wLeaf.LaNitrogenRelease * LaNitrogenDeadFrac;
		wLeaf.ShNitrogenDead = wLeaf.ShNitrogenRelease * ShNitrogenDeadFrac;
		wLeaf.LaNitrogenMassDeadCumu += wLeaf.LaNitrogenDead;
		wLeaf.ShNitrogenMassDeadCumu += wLeaf.ShNitrogenDead;

		//Z make the leaf drop
		wLeaf.dropLaMass = wLeaf.LaMass;
		wLeaf.dropShMass = wLeaf.ShMass;
		wLeaf.dropLaArea = wLeaf.seneLaArea;

		//Z set every growing stage, because we do not need this leaf anymore
		//  the last one mark the leaf death is forced by tiller death

		wLeaf.emerge = false;
		wLeaf.mature = false;
		wLeaf.senesce = true;
		wLeaf.dead = true;

		goto Leaf_LivingFracAdjustment;
	}

	//Step ONE pre-emerge -----------------------------------------------------------
	//     hidden part growth, growing pattern changes when the previous leaf is emerged
	//     after leaf is long enough to be presented from previous pseudostem(sheath), the current leaf is emerged
	if ((!wLeaf.dead) && wLeaf.initiation && (!wLeaf.emerge))
	{
		//Z compute the leaf growth rate
		//  if previous leaf not emerge, exponential-like growth function
		//  if previous leaf emerge, use beta function, compute SSLW and LSSW

		//Z first calculate some beta function parameters for preparation
		float beta_func_0 = 0.0f;
		float beta_func_t0 = 0.0f;
		float beta_func_t1 = 0.0f;
		beta_func_0 = abs((1.0f + (max(0.0f, (LEAF_TE - 0.0f)) / (LEAF_TE - LEAF_TM))) 
			* powf(min(1.0f, (0.0f - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));

		float delta_L = 0.0f;		//Z the length incremental cm 
		float delta_MSH = 0.0f;		//Z the sheath biomass incremental g 

		if ((!wLeaf.emerge_preL) && (!wLeaf.emerge_preL_done)) {

			//Z update the leaf "time based" life length, after previous leaf emergence
			//  there the "Beta" function growth starts
			wLeaf.LfPseudoAge += this->current_teq;

			if (wLeaf.firstupdate) {
				wLeaf.firstupdate = false;
				delta_L = 0.1f;
			}
			else {
				delta_L = wLeaf.LfLength * wLeaf.RERmax_L * this->current_teq
					 / (1.0f + RER_KC / wLeaf.Lfconce_Sucrose) / (1.0f + RER_KN / wLeaf.Lfconce_N);	//Z cm day^-1
			}

			wLeaf.LfLength += delta_L;
			wLeaf.LfPseudoStem = wLeaf.LfLength;
			wLeaf.ShLength = wLeaf.LfLength;

			//Z when the leaf is hidden, leaf does not have lamina and sheath yet, so just assume biomass demand goes to sheath
			//  need to use 0.01f to adjust "cm" to "m" (to use fspm functions)
			//  biomass in g and N mass in mg
			delta_MSH = LF_ALPHA * LF_BETA * powf(wLeaf.LfLength * 0.01f, LF_BETA - 1.0f) * delta_L * 0.01f;
			wLeaf.ptnShMassIncrease += delta_MSH;
			wLeaf.ptnShNitrogenMassIncrease += (delta_MSH * MAX_N_PCT * 10.0f);
		}
		if ((wLeaf.emerge_preL || ((wLeaf.LfPseudoAge >= 0.2f * PLASTOCHRONE) && (wLeaf.LfRank == 1))) && (!wLeaf.emerge_preL_done)) {
		//Z just run once, has a function to force the switch of the growth pattern
		//  for the forced condition, may need to restrict to the first two leaves
			wLeaf.emerge_preL_done = true;
			//Z based on fspm supplement doc, Table S1.1 Eq. 5-6, this should be the first time to compute
			//  1. the max leaf length determined at the previous leaf emergence (time point)
			wLeaf.maxLfLength_em = min(wLeaf.LfLength / beta_func_0, wLeaf.maxLfLength_ref);
			wLeaf.maxLfLength_em = 0.5f * (wLeaf.maxLfLength_em + wLeaf.maxLfLength_ref);  //Z this is a tricky term for making the total leaf length longer
			
			wLeaf.maxLfLength = wLeaf.maxLfLength_em;
			//  2. compute max leaf lamina and sheath scale
			wLeaf.maxLaLength = wLeaf.maxLfLength_em / (1.0f + wLeaf.SL_ratio);
			wLeaf.maxShLength = wLeaf.maxLaLength * wLeaf.SL_ratio;

			//Z reset the PseudoAge for future beta function growth
			wLeaf.LfPseudoAge = 0.0f;
		}
		if (wLeaf.emerge_preL_done) {
			//Z update the leaf "time based" life length, after previous leaf emergence
			//  there the "Beta" function growth starts
			wLeaf.LfPseudoAge += this->current_teq;

			if (wLeaf.LfPseudoAge <= LEAF_TB) {
				delta_L = wLeaf.LfLength - beta_func_0 * wLeaf.maxLfLength;
			}
			if (wLeaf.LfPseudoAge > LEAF_TB && wLeaf.LfPseudoAge <= LEAF_TE) {
				//Z the difference between to beta function is the amount of growth
				beta_func_t1= abs((1.0f + (max(0.0f, (LEAF_TE - wLeaf.LfPseudoAge)) / (LEAF_TE - LEAF_TM))) *
					powf(min(1.0f, (wLeaf.LfPseudoAge - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));
				beta_func_t0= abs((1.0f + (max(0.0f, (LEAF_TE - wLeaf.LfPseudoAge_prev)) / (LEAF_TE - LEAF_TM))) *
					powf(min(1.0f, (wLeaf.LfPseudoAge_prev - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));
				delta_L = min(wLeaf.maxLfLength, wLeaf.maxLfLength * (beta_func_t1 - beta_func_t0));
				delta_L = delta_L / (1.0f + RER_KC / wLeaf.Lfconce_Sucrose) / (1.0f + RER_KN / wLeaf.Lfconce_N);
			}
			if (wLeaf.LfPseudoAge > LEAF_TE) {
				delta_L = 0.0f;
			}

			delta_L = max(0.0f, delta_L);

			wLeaf.LfLength += delta_L;
			wLeaf.LfPseudoStem = wLeaf.LfLength;
			wLeaf.ShLength = wLeaf.LfLength;
			wLeaf.LfPseudoAge_prev = wLeaf.LfPseudoAge;

			//Z update the width length ratio
			float regul_WL_ratio = (LEAF_WL_REGUL_MAX - LEAF_WL_REGUL_MIN) / (LEAF_WL_INT_MAX - LEAF_WL_INT_MIN) * wLeaf.Lfconce_Sucrose_mean +
				(LEAF_WL_REGUL_MIN * LEAF_WL_INT_MAX - LEAF_WL_REGUL_MAX * LEAF_WL_INT_MIN) / (LEAF_WL_INT_MAX - LEAF_WL_INT_MIN);
			regul_WL_ratio = max(regul_WL_ratio, LEAF_WL_REGUL_MIN);
			regul_WL_ratio = min(regul_WL_ratio, LEAF_WL_REGUL_MAX);
			regul_WL_ratio *= LEAF_WL_BASE;

			//Z update the sheath specific linear weigth
	//		wLeaf.LSSW = (SH_LSSW_NOMINAL_A * static_cast<float>(wLeaf.LfRank) + SH_LSSW_NOMINAL_B)
	//					+ SH_LSSW_A * (wLeaf.Lfconce_Sucrose_mean - SH_LSSW_INT_MIN);
	//		wLeaf.LSSW = max(wLeaf.LSSW, SH_LSSW_MIN);
	//		wLeaf.LSSW = min(wLeaf.LSSW, SH_LSSW_MAX);
	//		wLeaf.LSSW *= LF_STRCTMASS_2_DRYMASS;

			//Z when the leaf is hidden, leaf does not have lamina and sheath yet, so just assume biomass demand goes to sheath
			delta_MSH = wLeaf.LSSW * delta_L;
			wLeaf.ptnShMassIncrease += delta_MSH;
			wLeaf.ptnShNitrogenMassIncrease += (delta_MSH * MAX_N_PCT * 10.0f);

			//Z update the leaf max length/width
			wLeaf.maxLfLength = wLeaf.LfLength + wLeaf.maxLfLength_em * (1.0f - beta_func_t1);
			wLeaf.maxLfWidth = min(wLeaf.maxLfLength * regul_WL_ratio, wLeaf.maxLfWidth_ref);
		}

		//Z the hidden part, before emergence will be sheath
		wLeaf.ShLength = wLeaf.LfPseudoStem;

		//Z decide if the leaf can emerge
		if (delta_L > wLeaf.Dist2Emerge) 
		{ 
			wLeaf.emerge = true; 
		}
		
		goto Leaf_LivingFracAdjustment;
	}

	//Step TWO -----------------------------------------------------------
	//	   leaf automate expansion
	//	   now leaf has "lamina" and "sheath",
	//	   while temporarily, we do not consider the mass of division zone
	//     and leaf pseudostem is just one portion of the 
	if ((!wLeaf.dead) && wLeaf.emerge && (!wLeaf.mature))
	{
		//Z leaf is emerged, part of the sheath will be (visible) sheath, while part of the sheath will be wLeaf.LfPseudoStem;
		wLeaf.LfPseudoStem = max(0.0f, ligulationDist);

		//Z compute the leaf growth rate
		//Z first calculate some beta function parameters for preparation
		float beta_func_0 = 0.0f;
		float beta_func_t0 = 0.0f;
		float beta_func_t1 = 0.0f;
		beta_func_0 = abs((1.0f + (max(0.0f, (LEAF_TE - 0.0f)) / (LEAF_TE - LEAF_TM))) *
			powf(min(1.0f, (0.0f - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));

		float delta_L;			//Z the length incremental cm 
		float delta_MLA;		//Z the lamina biomass incremental g
		float delta_MSH;		//Z the sheath biomass incremental g 

		//Z update the leaf "time based" life length, after this leaf emergence
		wLeaf.LfPseudoAge += this->current_teq;

		if (wLeaf.LfPseudoAge <= LEAF_TB) {
			delta_L = wLeaf.LfLength - beta_func_0 * wLeaf.maxLfLength;
		}
		if (wLeaf.LfPseudoAge > LEAF_TB && wLeaf.LfPseudoAge <= LEAF_TE) {
			beta_func_t1 = abs((1.0f + (max(0.0f, (LEAF_TE - wLeaf.LfPseudoAge)) / (LEAF_TE - LEAF_TM))) *
				powf(min(1.0f, (wLeaf.LfPseudoAge - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));
			beta_func_t0 = abs((1.0f + (max(0.0f, (LEAF_TE - wLeaf.LfPseudoAge_prev)) / (LEAF_TE - LEAF_TM))) *
				powf(min(1.0f, (wLeaf.LfPseudoAge_prev - LEAF_TB) / (LEAF_TE - LEAF_TB)), (LEAF_TE - LEAF_TB) / (LEAF_TE - LEAF_TM)));
			delta_L = min(wLeaf.maxLfLength, wLeaf.maxLfLength * (beta_func_t1 - beta_func_t0));
			delta_L = delta_L / (1.0f + RER_KC / wLeaf.Lfconce_Sucrose) / (1.0f + RER_KN / wLeaf.Lfconce_N);
		}
		if (wLeaf.LfPseudoAge > LEAF_TE) {
			delta_L = 0.0f;
		}

		delta_L = max(0.0f, delta_L);

		wLeaf.LfLength += delta_L;
		wLeaf.LfPseudoAge_prev = wLeaf.LfPseudoAge;
		
		//Z update the width length ratio
		float regul_WL_ratio = (LEAF_WL_REGUL_MAX - LEAF_WL_REGUL_MIN) / (LEAF_WL_INT_MAX - LEAF_WL_INT_MIN) * wLeaf.Lfconce_Sucrose_mean +
			(LEAF_WL_REGUL_MIN * LEAF_WL_INT_MAX - LEAF_WL_REGUL_MAX * LEAF_WL_INT_MIN) / (LEAF_WL_INT_MAX - LEAF_WL_INT_MIN);
		regul_WL_ratio = max(regul_WL_ratio, LEAF_WL_REGUL_MIN);
		regul_WL_ratio = min(regul_WL_ratio, LEAF_WL_REGUL_MAX);
		regul_WL_ratio *= LEAF_WL_BASE;

		//Z compute current leaf lamina length, sheath length, lamina width
		float LaLengthFraction = wLeaf.LaLength / wLeaf.LfLength;
		float LaLengthFraction_ref = 1.0f / (1.0f + wLeaf.SL_ratio);
		if (LaLengthFraction < LaLengthFraction_ref)
		{
			wLeaf.LaLength += delta_L;
			wLeaf.LaWidth = min(wLeaf.LaLength * regul_WL_ratio, wLeaf.maxLfWidth);

			float LaArea_new = wLeaf.LaLength * wLeaf.LaWidth * wLeaf.shapeFrac;
			wLeaf.ptnLaAreaIncrease = max(LaArea_new - wLeaf.LaArea, 0.0f);
			wLeaf.LaArea = LaArea_new;
			wLeaf.ptnShLengthIncrease = 0.0f;
		}
		else
		{
			wLeaf.LaLength += delta_L / (1.0f + wLeaf.SL_ratio);
			wLeaf.LaWidth = min(wLeaf.LaLength * regul_WL_ratio, wLeaf.maxLfWidth);
			wLeaf.ShLength += delta_L * wLeaf.SL_ratio / (1.0f + wLeaf.SL_ratio);

			float LaArea_new = wLeaf.LaLength * wLeaf.LaWidth * wLeaf.shapeFrac;
			wLeaf.ptnLaAreaIncrease = max(LaArea_new - wLeaf.LaArea, 0.0f);
			wLeaf.LaArea = LaArea_new;
			wLeaf.ptnShLengthIncrease = delta_L * wLeaf.SL_ratio / (1.0f + wLeaf.SL_ratio);
		}

		//Z update leaf SSLW (for leaf) and LSSW (for sheath)
	//	wLeaf.SSLW = (LF_SSLW_MAX - LF_SSLW_MIN) / (LF_SSLW_INT_MAX - LF_SSLW_INT_MIN) * wLeaf.Lfconce_Sucrose_mean
	//				+ (LF_SSLW_MIN * LF_SSLW_INT_MAX - LF_SSLW_MAX * LF_SSLW_INT_MIN) / (LF_SSLW_INT_MAX - LF_SSLW_INT_MIN);
	//	wLeaf.SSLW = max(wLeaf.SSLW, LF_SSLW_MIN);
	//	wLeaf.SSLW = min(wLeaf.SSLW, LF_SSLW_MAX);
	//	wLeaf.SSLW *= LF_STRCTMASS_2_DRYMASS;

	//	wLeaf.LSSW = (SH_LSSW_NOMINAL_A * static_cast<float>(wLeaf.LfRank + 1) + SH_LSSW_NOMINAL_B)
	//				+ SH_LSSW_A * (wLeaf.Lfconce_Sucrose_mean - SH_LSSW_INT_MIN);
	//	wLeaf.LSSW = max(wLeaf.LSSW, SH_LSSW_MIN);
	//	wLeaf.LSSW = min(wLeaf.LSSW, SH_LSSW_MAX);
	//	wLeaf.LSSW *= LF_STRCTMASS_2_DRYMASS;

		//Z leaf is emerge now, so compute both lamina and sheath
		delta_MLA = wLeaf.SSLW * wLeaf.ptnLaAreaIncrease;
		delta_MSH = wLeaf.LSSW * wLeaf.ptnShLengthIncrease;
		wLeaf.ptnLaMassIncrease += delta_MLA;
		wLeaf.ptnLaNitrogenMassIncrease += (delta_MLA * MAX_N_PCT * 10.0f);
		wLeaf.ptnShMassIncrease += delta_MSH;
		wLeaf.ptnShNitrogenMassIncrease += (delta_MSH * MAX_N_PCT * 10.0f);

		//Z update the leaf max length/width
		wLeaf.maxLfLength = wLeaf.LfLength + wLeaf.maxLfLength_em * (1.0f - beta_func_t1);
		wLeaf.maxLfWidth = max(wLeaf.maxLfLength * regul_WL_ratio, wLeaf.maxLfWidth_ref);

		if (abs(wLeaf.LfLength - wLeaf.maxLfLength) < 0.01f) {
			wLeaf.mature = true;
			wLeaf.emergePseudoAge = wLeaf.LfPseudoAge;
			//cout << wLeaf.emergePseudoAge << '\n';
			wLeaf.seneDuration = wLeaf.emergePseudoAge;
			wLeaf.mature_TTd = wLeaf.physAge;
			wLeaf.matureLaArea = wLeaf.LaArea;
			wLeaf.greenLaArea = wLeaf.LaArea;
		}
		goto Leaf_LivingFracAdjustment;
	}

	//Step Three -----------------------------------------------------------
	//     second leaf matured
	//     Leaf is matured and actively produce photosynthesis
	//     Always update the Leaf stay green period, in case a bad condition happens and leaf must die immediately
	if ((!wLeaf.dead) && wLeaf.mature && (!wLeaf.senesce))
	{
		wLeaf.activeAge += current_gdd;

		//Z leaf stay green period
		float staygreen = PHYLLOCHRON * max(6.5f * min(wLeaf.N_effect, psi_effect_sene), 2.5f);

		//Z update leaf mass, preserve green leaf area since leaf is active
		//  only green leaf area support photosynthesis
		wLeaf.greenLaArea = wLeaf.LaArea;
		wLeaf.greenLaMass = wLeaf.LaMass;
		wLeaf.greenShMass = wLeaf.ShMass;

		if (wLeaf.activeAge >= staygreen)
		{
			wLeaf.senesce = true;
			wLeaf.seneDuration = wLeaf.emergePseudoAge;
		}
		goto Leaf_LivingFracAdjustment;
	}

	//Step Four -----------------------------------------------------------
	//Z	   third leaf senesce
	//Z    the senesce time period is the same as the growth period
	if ((!wLeaf.dead) && wLeaf.senesce)
	{

		//Z during sene, the leaf will not accept new biomass and N
		wLeaf.ptnLaMassIncrease = 0.0f;
		wLeaf.ptnShMassIncrease = 0.0f;
		wLeaf.ptnLaNitrogenMassIncrease = 0.0f;
		wLeaf.ptnShNitrogenMassIncrease = 0.0f;

		//Z senesce start
		wLeaf.senePseudoAge += this->current_teq;
		wLeaf.seneDuration = max(0.0f, wLeaf.seneDuration - 0.5f * this->current_teq * (1.0f - min(psi_effect_sene, wLeaf.N_effect)));
		wLeaf.senePseudoAge = min(wLeaf.senePseudoAge, wLeaf.seneDuration);
		wLeaf.seneLaArea = max(0.0f, wLeaf.matureLaArea * (1.0f + 2.0f * (1.0f - wLeaf.senePseudoAge / wLeaf.seneDuration)) * powf(wLeaf.senePseudoAge / wLeaf.seneDuration, 2.0f));

		if ((wLeaf.seneLaArea > wLeaf.matureLaArea) || (wLeaf.senePseudoAge >= wLeaf.seneDuration))
		{
			wLeaf.seneLaArea = wLeaf.matureLaArea;
			wLeaf.dropLaArea = wLeaf.matureLaArea;
			wLeaf.dropLaMass = wLeaf.LaMass;
			wLeaf.dropShMass = wLeaf.ShMass;
			wLeaf.dead = true;
			wLeaf.livingFrac = 0.0f;
		}

		//Z leaf area decreasing
		wLeaf.ptnLaAreaDecrease = max(0.0f, wLeaf.seneLaArea + wLeaf.greenLaArea - wLeaf.matureLaArea);
		//Z this is the relative sene fraction relative to the current living portion of the plant organ, max to one
		float relative_sene_fraction = min(wLeaf.ptnLaAreaDecrease / wLeaf.greenLaArea, 1.0f);
		//Z then update the green/living leaf portion
		wLeaf.greenLaArea = max(0.0f, wLeaf.matureLaArea - wLeaf.seneLaArea);

		//Z determine N release in mg
		// Release: single leaf N output due to senecense and dropping, need to divide into remobile(released) portion and dead(with leaf drop) portion
		// Dead : at that step, reduced N goes to residue and immobilized
		// Return : at that step, remobilized and can reused for new plant organs
		// DeadCumu : cumulated N mass in the senecense / dead plant organs
		
		wLeaf.LaNitrogenRelease = relative_sene_fraction * wLeaf.LaNitrogenMass;
		wLeaf.ShNitrogenRelease = relative_sene_fraction * wLeaf.ShNitrogenMass;
		wLeaf.LaNitrogenMass = wLeaf.LaNitrogenMass - wLeaf.LaNitrogenRelease;
		wLeaf.ShNitrogenMass = wLeaf.ShNitrogenMass - wLeaf.ShNitrogenRelease;
		//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
		float LaNitrogenReturnFrac = min(max(wLeaf.LaNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.LaNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
		float LaNitrogenDeadFrac = 1.0f - LaNitrogenReturnFrac;
		float ShNitrogenReturnFrac = min(max(wLeaf.ShNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.ShNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
		float ShNitrogenDeadFrac = 1.0f - ShNitrogenReturnFrac;
		wLeaf.LaNitrogenReturn = wLeaf.LaNitrogenRelease * LaNitrogenReturnFrac;
		wLeaf.ShNitrogenReturn = wLeaf.ShNitrogenRelease * ShNitrogenReturnFrac;
		wLeaf.LaNitrogenDead = wLeaf.LaNitrogenRelease * LaNitrogenDeadFrac;
		wLeaf.ShNitrogenDead = wLeaf.ShNitrogenRelease * ShNitrogenDeadFrac;
		wLeaf.LaNitrogenMassDeadCumu += wLeaf.LaNitrogenDead;
		wLeaf.ShNitrogenMassDeadCumu += wLeaf.ShNitrogenDead;

		//Z determine Leaf biomass in g
		wLeaf.greenLaMass = (1.0f - relative_sene_fraction) * wLeaf.greenLaMass;
		wLeaf.greenShMass = (1.0f - relative_sene_fraction) * wLeaf.greenShMass;

		goto Leaf_LivingFracAdjustment;
	}

	//From Single Plant to Representative Plant, Mapping -----------------------------------------------------------
	//Z Simulation of representative leaf, adjusted by the livingFrac
Leaf_LivingFracAdjustment:

	wLeaf.ptnLaMassIncrease_Rep = 0.0f;
	wLeaf.ptnShMassIncrease_Rep = 0.0f;
	wLeaf.ptnLaNitrogenMassIncrease_Rep = 0.0f;
	wLeaf.ptnShNitrogenMassIncrease_Rep = 0.0f;
	wLeaf.LaNitrogenReturn_Rep = 0.0f;
	wLeaf.ShNitrogenReturn_Rep = 0.0f;

	//Z livingfrac must be smaller than or equal to livingfrac old
	float livingDiff = max(wLeaf.livingFrac_old - wLeaf.livingFrac, 0.0f);

	//-------------------------
	//Z adjust living leaf geometry, mass and potential mass required for leaf growth 
	//  green part
	wLeaf.greenLaArea_Rep = wLeaf.greenLaArea * wLeaf.livingFrac;
	wLeaf.greenLaMass_Rep = wLeaf.greenLaMass * wLeaf.livingFrac;
	wLeaf.greenShMass_Rep = wLeaf.greenShMass * wLeaf.livingFrac;

	//  senecent part
	//  senecent leaf portion has to be partitioned
	//	1. senecent of living leaf
	//	2. dead/dropped whole leaf
	//	This is a cumulative value so must be initialized to 0
	wLeaf.seneLaArea_Rep += (wLeaf.ptnLaAreaDecrease * wLeaf.livingFrac + (wLeaf.greenLaArea + wLeaf.ptnLaAreaDecrease) * livingDiff);
	
	//  drop part, for the equation below,
	//  1. if tiller dead, leaf will be dropped, so "DropLfArea_Rep" is a cumulative value
	//	2. current living leaf can drop, shown in the first term
	//	3. leafArea is everything, i.e., green, senecent and dead/dropped (exclusive), so no matter how those leaves shown, once tiller dead, it goes to the second term
	//	4. "exclusive" in 3 means if DropLfArea!=0, it must = LeafArea = SenescentLfArea, and livingFrac = livingDiff = 0
	wLeaf.dropLaArea_Rep += (wLeaf.dropLaArea * wLeaf.livingFrac + wLeaf.matureLaArea * livingDiff);

	//  total leaf area, SenescentLfArea includes but not limit to drop leaf area
	wLeaf.LaArea_Rep = wLeaf.greenLaArea_Rep + wLeaf.seneLaArea_Rep;

	//-------------------------
	//Z adjust leaf/sheath biomass
	//  1. drop part, similar to the "DropLfArea", this should be a cumulative value
	//  2. note that usually (until the last step), DropLfMass should be always 0
	wLeaf.dropLaMass_Rep += (wLeaf.dropLaMass * wLeaf.livingFrac + wLeaf.LaMass * livingDiff);
	wLeaf.dropShMass_Rep += (wLeaf.dropShMass * wLeaf.livingFrac + wLeaf.ShMass * livingDiff);

	//  total leaf/sheath mass
	//	1. if leaf is dropped, then rely on the second term
	//	2. if leaf is living, "DropLfMass" without caring the "livingfrac" will be 0, then both term shows
	//	3. similar to the sheath
	wLeaf.LaMass_Rep = (wLeaf.LaMass - wLeaf.dropLaMass) * wLeaf.livingFrac + wLeaf.dropLaMass_Rep;
	wLeaf.ShMass_Rep = (wLeaf.ShMass - wLeaf.dropShMass) * wLeaf.livingFrac + wLeaf.dropShMass_Rep;

	//-------------------------
	//Z adjust leaf/sheath nitrogen mass and their contribution to nitrogen releasing
	//  mass computation is simple, since we assume only dead portion totally recycle nitrogen
	wLeaf.LaNitrogenMass_Rep = wLeaf.LaNitrogenMass * wLeaf.livingFrac;
	wLeaf.ShNitrogenMass_Rep = wLeaf.ShNitrogenMass * wLeaf.livingFrac;

	//  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
	//	2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
	//	3. NitrogenDecline is just a stepwise variable, and should eventually put into the dead leaf or dead sheath N portion
	wLeaf.deadLaNitrogenMass_Rep += wLeaf.LaNitrogenDead * wLeaf.livingFrac;
	wLeaf.deadShNitrogenMass_Rep += wLeaf.ShNitrogenDead * wLeaf.livingFrac;
	wLeaf.LaNitrogenReturn_Rep = wLeaf.LaNitrogenReturn * wLeaf.livingFrac;
	wLeaf.ShNitrogenReturn_Rep = wLeaf.ShNitrogenReturn * wLeaf.livingFrac;
	//Z should be (LeafNitrogenMass+LfNitrogenDecline+LeafNitrogenRelease), 
	//  and those additional terms "LfNitrogenDecline+LeafNitrogenRelease=LeafNitrogenReduce_Single" can be higher order infinitesimal
	float LaNitrogenFree_temp = (wLeaf.LaNitrogenMass + wLeaf.LaNitrogenRelease) * livingDiff;
	float ShNitrogenFree_temp = (wLeaf.ShNitrogenMass + wLeaf.ShNitrogenRelease) * livingDiff;

	//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
	float LaNitrogenReturnFrac = min(max(wLeaf.LaNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.LaNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
	float LaNitrogenDeadFrac = 1.0f - LaNitrogenReturnFrac;
	float ShNitrogenReturnFrac = min(max(wLeaf.ShNitrogenContent - wLeaf.NitrogenReleaseLowBdd, 0.0f) / wLeaf.ShNitrogenContent, wLeaf.NitrogenReleaseMaxPtge);
	float ShNitrogenDeadFrac = 1.0f - ShNitrogenReturnFrac;
	wLeaf.LaNitrogenReturn_Rep += LaNitrogenFree_temp * LaNitrogenReturnFrac;
	wLeaf.ShNitrogenReturn_Rep += ShNitrogenFree_temp * ShNitrogenReturnFrac;

	// This part of decline N will be the dead,
	// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller 
	// so leaf on dead tiller means "dropped leaf", while "dead leaf" can be dropped or still on the tiller (yellow or sene leaf)
	float aaaa = LaNitrogenFree_temp * LaNitrogenDeadFrac;
	float bbbb = ShNitrogenFree_temp * ShNitrogenDeadFrac;
	wLeaf.deadLaNitrogenMass_Rep += aaaa;
	wLeaf.deadShNitrogenMass_Rep += bbbb;
	wLeaf.dropLaNitrogenMass_Rep += aaaa;
	wLeaf.dropShNitrogenMass_Rep += bbbb;

	// Nitorgen in Sene(Yellow) Leaf, dead (senescent and non-mobile anymore) but not dropped, should be the difference
	// sene(yellow) = Dead - Dropped
	wLeaf.seneLaNitrogenMass_Rep = wLeaf.deadLaNitrogenMass_Rep - wLeaf.dropLaNitrogenMass_Rep;
	wLeaf.seneShNitrogenMass_Rep = wLeaf.deadShNitrogenMass_Rep - wLeaf.dropShNitrogenMass_Rep;

	//-------------------------
	//Z adjust potential leaf growth
	//  only living portion can grow
	wLeaf.ptnLaMassIncrease_Rep = wLeaf.ptnLaMassIncrease * wLeaf.livingFrac;
	wLeaf.ptnLaNitrogenMassIncrease_Rep = wLeaf.ptnLaNitrogenMassIncrease * wLeaf.livingFrac;
	wLeaf.ptnShMassIncrease_Rep = wLeaf.ptnShMassIncrease * wLeaf.livingFrac;
	wLeaf.ptnShNitrogenMassIncrease_Rep = wLeaf.ptnShNitrogenMassIncrease * wLeaf.livingFrac;

	//Z finally update the living fraction numbers for this leaf
	wLeaf.livingFrac_old = wLeaf.livingFrac;

	if (wLeaf.dead)
	{
		wLeaf.livingFrac_old = 0.0f;
		wLeaf.livingFrac = 0.0f;

		wLeaf.greenLaArea = 0.0f;
		wLeaf.greenLaMass = 0.0f;
		wLeaf.greenShMass = 0.0f;
		wLeaf.LaNitrogenMass = 0.0f;
		wLeaf.ShNitrogenMass = 0.0f;

		wLeaf.ptnLaMassIncrease_Rep = 0.0f;
		wLeaf.ptnShMassIncrease_Rep = 0.0f;
		wLeaf.ptnLaNitrogenMassIncrease_Rep = 0.0f;
		wLeaf.ptnShNitrogenMassIncrease_Rep = 0.0f;

		wLeaf.seneLaNitrogenMass_Rep = 0.0f;
		wLeaf.seneShNitrogenMass_Rep = 0.0f;
		wLeaf.dropLaNitrogenMass_Rep = wLeaf.deadLaNitrogenMass_Rep;
		wLeaf.dropShNitrogenMass_Rep = wLeaf.deadShNitrogenMass_Rep;
	}

	if (wLeaf.livingFrac == 0.0f)
	{
		wLeaf.dead = true;
		wLeaf.livingFrac_old = 0.0f;

		wLeaf.greenLaArea = 0.0f;
		wLeaf.greenLaMass = 0.0f;
		wLeaf.greenShMass = 0.0f;
		wLeaf.LaNitrogenMass = 0.0f;
		wLeaf.ShNitrogenMass = 0.0f;

		wLeaf.ptnLaMassIncrease_Rep = 0.0f;
		wLeaf.ptnShMassIncrease_Rep = 0.0f;
		wLeaf.ptnLaNitrogenMassIncrease_Rep = 0.0f;
		wLeaf.ptnShNitrogenMassIncrease_Rep = 0.0f;

		wLeaf.seneLaNitrogenMass_Rep = 0.0f;
		wLeaf.seneShNitrogenMass_Rep = 0.0f;
		wLeaf.dropLaNitrogenMass_Rep = wLeaf.deadLaNitrogenMass_Rep;
		wLeaf.dropShNitrogenMass_Rep = wLeaf.deadShNitrogenMass_Rep;
	}

}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Leaf Mass Distribution Operator Structure

//Z update leaf mass and 
//  leaf stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatLeafMassAssignment::SetLeafMassAssignmentCondition(bool accel, float biomassRate, float nitrogenRate)
{
	//Z receive biomass and nitrogen distribution rates for leaves over the entire plant based plant-level biomass allocation
	acceleration = accel;
	biomassIncomeRate = biomassRate;
	nitrogenIncomeRate = nitrogenRate;
}

//Z leaf mass partition
//  should include leaf and sheath
//  the input parameters are based on the "representative" values of the leaf,
//  but leaf growth and mass should be baesd on "one single leaf"
//  Thus, we need to make a conversion between "Representative" and "one single leaf" using livingfrac
void WheatLeafMassAssignment::LeafMassRun(WheatLeaf& wLeaf)
const {
	//Z reset for each time step
	//  store biomass (g) and N (mg) allocated from plant
	wLeaf.LaBiomassIncrease = 0.0f;
	wLeaf.ShBiomassIncrease = 0.0f;
	wLeaf.LaNitrogenIncrease = 0.0f;
	wLeaf.ShNitrogenIncrease = 0.0f;

	//Z separate cases explicitly, good for GPU stream control.
	// BIOMASS in g
	if (wLeaf.ptnLaMassIncrease_Rep > 0.0f)
	{
		//Z this "adjustment" variable is necessary and important,
		//  essentially, " / wLeaf.livingFrac" is the important term
		//  this convert the mass allocation for a "representative leaf" to "individual"
		wLeaf.LaBiomassIncrease = wLeaf.ptnLaMassIncrease_Rep * biomassIncomeRate / wLeaf.livingFrac;
	}
	if (wLeaf.ptnShMassIncrease_Rep > 0.0f)
	{
		wLeaf.ShBiomassIncrease = wLeaf.ptnShMassIncrease_Rep * biomassIncomeRate / wLeaf.livingFrac;
	}

	// NITROGEM in mg
	if (wLeaf.ptnLaNitrogenMassIncrease_Rep > 0.0f)
	{
		wLeaf.LaNitrogenIncrease = wLeaf.ptnLaNitrogenMassIncrease_Rep * nitrogenIncomeRate / wLeaf.livingFrac;
	}
	if (wLeaf.ptnShNitrogenMassIncrease_Rep > 0.0f)
	{
		wLeaf.ShNitrogenIncrease = wLeaf.ptnShNitrogenMassIncrease_Rep * nitrogenIncomeRate / wLeaf.livingFrac;
	}

	//Z update leaf and sheath bio mass, N mass
	wLeaf.LaMass += wLeaf.LaBiomassIncrease;
	wLeaf.ShMass += wLeaf.ShBiomassIncrease;
	wLeaf.greenLaMass += wLeaf.LaBiomassIncrease;
	wLeaf.greenShMass += wLeaf.ShBiomassIncrease;
	wLeaf.LaNitrogenMass += wLeaf.LaNitrogenIncrease;
	wLeaf.ShNitrogenMass += wLeaf.ShNitrogenIncrease;

	//Z update leaf and sheath ptn mass and N mass, i.e., the biomass and N demand
	wLeaf.ptnLaMassIncrease -= wLeaf.LaBiomassIncrease;
	wLeaf.ptnShMassIncrease -= wLeaf.ShBiomassIncrease;
	wLeaf.ptnLaNitrogenMassIncrease -= wLeaf.LaNitrogenIncrease;
	wLeaf.ptnShNitrogenMassIncrease -= wLeaf.ShNitrogenIncrease;


	//Z update leaf and sheath N content
	//Z OK to over the total leafmass since for a single plant leaf, at this time, leaf is still growing in mass
	//  0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
	if (wLeaf.LaMass > 0) { wLeaf.LaNitrogenContent = max(min(0.10f * wLeaf.LaNitrogenMass / wLeaf.LaMass, MAX_N_PCT), 0.10f); }
	else { wLeaf.LaNitrogenContent = MAX_N_PCT; }
	if (wLeaf.ShMass > 0) { wLeaf.ShNitrogenContent = max(min(0.10f * wLeaf.ShNitrogenMass / wLeaf.ShMass, MAX_N_PCT), 0.10f); }
	else { wLeaf.ShNitrogenContent = MAX_N_PCT; }
}
