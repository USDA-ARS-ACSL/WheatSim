#include "WheatChaff.h"

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Chaff Data Structure

WheatChaff::WheatChaff() :
	livingFrac(0.0f), livingFrac_old(0.0f), livingFrac_ini(0.0f)
{
	mature = false;             //Z mark the max chaff lengh is reached
	dead = false;
	force_to_death_current_step = false;

	cur_TTd = 0.0f;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	physAge = 0.0f;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	cumu_growthTdd = 0.0f;		//Z cumu gdd growth age, stressed
	elongAge = 0.0f;			//						unstressed	

	//Z Chaff Mass g
	CfMass = 0.0f;
	deadCfMass = 0.0f;
	CfMassIncrease = 0.0f;
	ptnCfMassIncrease = 0.0f;

	//Z Chaff Nitrogen Mass mg
	CfNitrogenMass = 0.0f;
	CfNitrogenIncrease = 0.0f;
	ptnCfNitrogenMassIncrease = 0.0f;

	//Z Chaff Nitrogne conent and release fraction
	CfNitrogenContent = MAX_N_PCT;
	CfNitrogenReleaseLowBdd = 0.5f;			//Z the min chaff N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	CfNitrogenReleaseMaxPtge = 0.8f;		//Z the chaff N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	CfNitrogenRelease = 0.0f;
	CfNitrogenDead = 0.0f;				//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	CfNitrogenReturn = 0.0f;			//Z stepwise returned portion of leaf nitrogen, will supply plant future usage 
	CfNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of leaf nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	N_effect = 1.0;

	//Z ----------- REPRESENTATIVE Chaff ---------------------------------------

	//Z representative mass values
	CfMass_Rep = 0.0f;
	deadCfMass_Rep = 0.0f;
	CfNitrogenMass_Rep = 0.0f;

	//Z for chaff, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	CfNitrogenReturn_Rep = 0.0f;
	deadCfNitrogenMass_Rep = 0.0f;

	ptnCfMassIncrease_Rep = 0.0f;
	ptnCfNitrogenMassIncrease_Rep = 0.0f;
}

WheatChaff::WheatChaff(float livingFrac) :
	livingFrac(livingFrac), livingFrac_old(livingFrac), livingFrac_ini(livingFrac)
{
	mature = false;             //Z mark the max chaff lengh is reached
	dead = false;
	force_to_death_current_step = false;

	cur_TTd = 0.0f;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	physAge = 0.0f;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	cumu_growthTdd = 0.0f;		//Z cumu gdd growth age, stressed
	elongAge = 0.0f;			//						unstressed	

	//Z Chaff Mass g
	CfMass = 0.0f;
	deadCfMass = 0.0f;
	CfMassIncrease = 0.0f;
	ptnCfMassIncrease = 0.0f;

	//Z Chaff Nitrogen Mass mg
	CfNitrogenMass = 0.0f;
	CfNitrogenIncrease = 0.0f;
	ptnCfNitrogenMassIncrease = 0.0f;

	//Z Chaff Nitrogne conent and release fraction
	CfNitrogenContent = MAX_N_PCT;
	CfNitrogenReleaseLowBdd = 0.5f;			//Z the min chaff N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	CfNitrogenReleaseMaxPtge = 0.8f;		//Z the chaff N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	CfNitrogenRelease = 0.0f;
	CfNitrogenDead = 0.0f;				//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	CfNitrogenReturn = 0.0f;			//Z stepwise returned portion of leaf nitrogen, will supply plant future usage 
	CfNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of leaf nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	N_effect = 1.0;

	//Z ----------- REPRESENTATIVE Chaff ---------------------------------------

	//Z representative mass values
	CfMass_Rep = 0.0f;
	deadCfMass_Rep = 0.0f;
	CfNitrogenMass_Rep = 0.0f;

	//Z for chaff, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	CfNitrogenReturn_Rep = 0.0f;
	deadCfNitrogenMass_Rep = 0.0f;

	ptnCfMassIncrease_Rep = 0.0f;
	ptnCfNitrogenMassIncrease_Rep = 0.0f;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Chaff Update Operator Structure

WheatChaffUpdate::WheatChaffUpdate() :
	current_gdd(0.0f), predawn_psi(0.0f), psi_effect_elong(1.0f), growthDuration(0.0f), chaffEnd(false) {};

//Z update chaff growth gdd and conditions for this current hourly step
//  stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatChaffUpdate::SetChaffUpdateCondition(float gdd, float psi_predawn, float growthDur, bool chaffEnd)
{
	current_gdd = gdd;
	predawn_psi = psi_predawn;
	//Z should be gdd based. "g/gdd"
	growthDuration = growthDur;

	this->chaffEnd = chaffEnd;

	//Z compute water potential effect for expanding and senescence
	//Z some constant first
	const float psi_f = -1.4251f; // -1.0, was -1.4251   later changed to -2.3;
	const float s_f = 0.4258f; // was 0.4258 0.5;
	//Z expanding
	const float psi_threshold_bars_expand = -0.8657f;	//Mpa
	psi_effect_elong = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_expand))));
	psi_effect_elong = min(max(psi_effect_elong, 0.1f), 1.0f);

}

//Z this is the main function for plant chaff growth
//  every time updating, hourly gdd, tmpr, water potential, shade are universal
//                       N stress are local
//                       living fractor is the input variable for each chafff object
/*Z Input tuple for this function includes
*			1. Chaff structure;
*           2. current living fraction (from tiller);
*           3. force_to_death by the plant;
*/
void WheatChaffUpdate::ChafUpdateRun(WheatChaff& wChaf, float lf, bool f2d)
const {
	//Z Simulation of chaff
	//WheatChaff& wChaf = std::get<0>(t);
	//wChaf.livingFrac = std::get<1>(t);
	//wChaf.force_to_death_current_step = std::get<2>(t);

	wChaf.livingFrac = lf;
	wChaf.force_to_death_current_step = f2d;

	//Z reset some critical parameters
	wChaf.ptnCfMassIncrease = 0.0f;
	wChaf.ptnCfNitrogenMassIncrease = 0.0f;

	//Z chaff finished it growth, no need to update
	if (chaffEnd) return;

	//Z compute the N stress for each inidividual leaves, N in mass fraction (mg/g biomass)
	wChaf.N_effect = 2.0f / (1.0f + expf(-2.9f * (max(MIN_N_PCT, wChaf.CfNitrogenContent) - MIN_N_PCT))) - 1.0f;
	wChaf.N_effect = max(min(wChaf.N_effect, 1.0f), 0.1f);

	//Step ZERO Force to death ---------------------------------------------
		//Z this condition is imposed by the parent tiller, stronger than any other conditions
			//Z prevent kill the leaf twice
	if (!wChaf.dead && wChaf.force_to_death_current_step)
	{
		//Z DO NOT change "livingFrac" here
		//  DO NOT set "GreenInLength, InNitrogenMass" to 0
		//  That is because we need those values to calculate N return and convert green Chaff length to dead length in the "LivingFractionAdjustment" function
		//  The correctness can be viewed by computing N mass fraction from the output.
		//  The N mass% ranges 0.25% to 4.0%

		/*Z determine N release in mg
		* Release: single Chaff N output due to senecense and dropping, need to divide into remobile(released) portion and dead(within yellow stem) portion
		* Dead : at that step, reduced N goes to residue and immobilized
		* Return : at that step, remobilized and can reused for new plant organs
		* DeadCumu : cumulated N mass in the senecense / dead plant organs
		*/

		wChaf.CfNitrogenRelease = wChaf.CfNitrogenMass;
		wChaf.CfNitrogenMass = 0.0f;
		float NitrogenReturnFrac = min(max(wChaf.CfNitrogenContent - wChaf.CfNitrogenReleaseLowBdd, 0.0f) / wChaf.CfNitrogenContent, wChaf.CfNitrogenReleaseMaxPtge);
		float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
		wChaf.CfNitrogenReturn = wChaf.CfNitrogenRelease * NitrogenReturnFrac;
		wChaf.CfNitrogenDead = wChaf.CfNitrogenRelease * NitrogenDeadFrac;
		wChaf.CfNitrogenMassDeadCumu += wChaf.CfNitrogenDead;

		//Z make the internode dead
		wChaf.deadCfMass = wChaf.CfMass;

		//Z set every growing stage, because we do not need to operate this chaff
		wChaf.mature = false;
		wChaf.dead = true;

		goto Chaf_LivingFracAdjustment;
	}

	//Step First -----------------------------------------------------------
		//Z first chaff growth
	if ((!wChaf.dead) && (!wChaf.mature))
	{
		wChaf.cumu_growthTdd += min(wChaf.N_effect, psi_effect_elong) * current_gdd;
		wChaf.elongAge += current_gdd;
		wChaf.ptnCfMassIncrease = max(min(1.0f, wChaf.cumu_growthTdd / growthDuration) * MAXCHAFFWEIGHT - wChaf.CfMass, 0.0f);
		wChaf.ptnCfNitrogenMassIncrease = max(min(1.0f, wChaf.cumu_growthTdd / growthDuration) * MAXCHAFFWEIGHT * MAX_N_PCT * 10.0f - wChaf.CfNitrogenMass, 0.0f);

		if (wChaf.elongAge >= growthDuration)
		{
			wChaf.mature = true;
		}

		goto Chaf_LivingFracAdjustment;
	}

	//Step Second -----------------------------------------------------------
		//Z second chaff
		//  if chaff is matured, no need to update mass like leaf and internode

	//Step Third -----------------------------------------------------------
		//Z third chaff senesce
		//  chaff never dies itself, unless the whole plant dead

	//Step Four Representative Mapping -----------------------------------------------------------
Chaf_LivingFracAdjustment:

	wChaf.ptnCfMassIncrease_Rep = 0.0f;
	wChaf.ptnCfNitrogenMassIncrease_Rep = 0.0f;
	wChaf.CfNitrogenReturn_Rep = 0.0f;

	//Z livingfrac must be smaller than or equal to livingfrac old
	float livingDiff = max(wChaf.livingFrac_old - wChaf.livingFrac, 0.0f);

	//-------------------------
	//Z adjust chaff biomass, only dead and living parts
	/*  1. dead part, should be a cumulative value
	*   2. note that usually (until the last step), deadCfMass for chaff should be always 0
	*/
	wChaf.deadCfMass_Rep += (wChaf.deadCfMass * wChaf.livingFrac + wChaf.CfMass * livingDiff);

	//  total chaff mass
	wChaf.CfMass_Rep = (wChaf.CfMass - wChaf.deadCfMass) * wChaf.livingFrac + wChaf.deadCfMass_Rep;

	//-------------------------
	/*Z adjust chaff nitrogen mass and their contribution to nitrogen releasing
		1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
		*/
	wChaf.CfNitrogenMass_Rep = wChaf.CfNitrogenMass * wChaf.livingFrac;

	/*  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
		2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
		3. NitrogenDecline is just a stepwise variable, and should eventually put into the dead leaf or dead sheath N portion
	*/
	wChaf.deadCfNitrogenMass_Rep += wChaf.CfNitrogenDead * wChaf.livingFrac;
	wChaf.CfNitrogenReturn_Rep = wChaf.CfNitrogenReturn * wChaf.livingFrac;
	//Z N release due to the living fraction changes
	float CfNitrogenFree_temp = (wChaf.CfNitrogenMass + wChaf.CfNitrogenRelease) * livingDiff;

	//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
	float NitrogenReturnFrac = min(max(wChaf.CfNitrogenContent - wChaf.CfNitrogenReleaseLowBdd, 0.0f) / wChaf.CfNitrogenContent, wChaf.CfNitrogenReleaseMaxPtge);
	float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
	wChaf.CfNitrogenReturn_Rep += CfNitrogenFree_temp * NitrogenReturnFrac;

	// This part of decline N will be the dead,
	// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller 
	// so chaff on dead tiller means "dead chaff" (may be not dropped)
	wChaf.deadCfNitrogenMass_Rep += (CfNitrogenFree_temp * NitrogenDeadFrac);

	//-------------------------
	//Z adjust potential chaff
	//  only living portion can grow

	wChaf.ptnCfMassIncrease_Rep = wChaf.ptnCfMassIncrease * wChaf.livingFrac;
	wChaf.ptnCfNitrogenMassIncrease_Rep = wChaf.ptnCfNitrogenMassIncrease * wChaf.livingFrac;

	//Z finally update the living fraction numbers for this leaf
	wChaf.livingFrac_old = wChaf.livingFrac;

	if (wChaf.dead)
	{
		wChaf.livingFrac_old = 0.0f;
		wChaf.livingFrac = 0.0f;
		wChaf.CfNitrogenMass = 0.0f;
		wChaf.ptnCfMassIncrease_Rep = 0.0f;
		wChaf.ptnCfNitrogenMassIncrease_Rep = 0.0f;
	}
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Chaff Mass Distribution Operator Structure

//Z update chaff mass and 
//  chaff stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatChaffMassAssignment::SetChaffMassAssignmentCondition(float biomassRate, float nitrogenRate)
{
	//Z receive biomass and nitrogen distribution rates for internode over the entire plant based plant-level biomass allocation
	biomassIncomeRate = biomassRate;
	nitrogenIncomeRate = nitrogenRate;
}

//Z chaff mass partition
//  the input parameters are based on the "representative" values of the chaff,
//  but chaff growth and mass should be baesd on "one individual chaff"
//  Thus, we need to make a conversion between "Representative" and "one single chaff" using livingfrac
void WheatChaffMassAssignment::ChafMassRun(WheatChaff& wChaf)
const {
	//Z reset for each time step
//  store biomass (g) and N (mg) allocated from plant
	wChaf.CfMassIncrease = 0.0f;
	wChaf.CfNitrogenIncrease = 0.0f;

	//Z separate cases explicitly, good for GPU stream control.
	// BIOMASS
	// biomass incremental
	if (wChaf.ptnCfMassIncrease_Rep > 0.0f)
	{
		//Z this "adjustment" variable is necessary and important,
		//  essentially, " / wChaf.livingFrac" is the important term
		//  this convert the mass allocation for a "representative chaff" to "individual chaff"
		wChaf.CfMassIncrease = wChaf.ptnCfMassIncrease_Rep * biomassIncomeRate / wChaf.livingFrac;	// mass in g
	}
	// NITROGEM
	// N incremental
	if (wChaf.ptnCfNitrogenMassIncrease_Rep > 0.0f)
	{
		wChaf.CfNitrogenIncrease = wChaf.ptnCfNitrogenMassIncrease_Rep * nitrogenIncomeRate / wChaf.livingFrac;	// mass in mg
	}

	//Z update chaff bio mass, N mass
	wChaf.CfMass += wChaf.CfMassIncrease;
	wChaf.CfNitrogenMass += wChaf.CfNitrogenIncrease;

	//Z update internnode N content
	//Z OK to over the total internode mass since for a single plant internode, at this time, inernode is still growing in mass
	//  0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
	if (wChaf.CfMass) { wChaf.CfNitrogenContent = max(min(0.10f * wChaf.CfNitrogenMass / wChaf.CfMass, MAX_N_PCT), 0.1f); }
	else { wChaf.CfNitrogenContent = MAX_N_PCT; }
	return;
}
