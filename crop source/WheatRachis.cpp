#include "WheatRachis.h"

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Rachis Data Structure

WheatRachis::WheatRachis() :
	livingFrac(0.0f), livingFrac_old(0.0f), livingFrac_ini(0.0f)
{
	mature = false;             //Z mark the max chaff lengh is reached
	dead = false;
	force_to_death_current_step = false;

	cur_TTd = 0.0f;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	physAge = 0.0f;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	elongAge = 0.0f;
	cumu_growthTdd = 0.0f;		//Z cumu gdd growth age, stressed

	//Z Rachis Length (cm)
	RaLength = 0.0f;				//Z rachis length in cm
	ptnRaLengthIncrease = 0.0f;		//Z potential rachis length increase
	RaLengthIncrease = 0.0f;		//Z actual Rachis length increase

	//Z Rachis Mass g
	srl = 0.01f;
	RaMass = 0.0f;
	deadRaMass = 0.0f;
	RaMassIncrease = 0.0f;
	ptnRaMassIncrease = 0.0f;

	//Z Rachis Nitrogen Mass mg
	RaNitrogenMass = 0.0f;
	RaNitrogenIncrease = 0.0f;
	ptnRaNitrogenMassIncrease = 0.0f;

	//Z Rachis Nitrogne conent and release fraction
	RaNitrogenContent = MAX_N_PCT;
	RaNitrogenReleaseLowBdd = 0.5f;			//Z the min rachis N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	RaNitrogenReleaseMaxPtge = 0.8f;		//Z the rachis N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	RaNitrogenRelease = 0.0f;
	RaNitrogenDead = 0.0f;				//Z stepwise dead portion of rachis nitrogen, will be in the dead plant tissue and not removable
	RaNitrogenReturn = 0.0f;			//Z stepwise returned portion of rachis nitrogen, will supply plant future usage 
	RaNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of rachis nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	N_effect = 1.0;

	//Z ----------- REPRESENTATIVE Rachis ---------------------------------------

	//Z Rachis length
	RaLength_Rep = 0.0f;	//Z this should be a cumulative value, as the living fraction changes, 
	//  it cumulated different amount (based on living fraction)

//Z representative mass values
	RaMass_Rep = 0.0f;
	deadRaMass_Rep = 0.0f;
	RaNitrogenMass_Rep = 0.0f;

	//Z for rachis, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	RaNitrogenReturn_Rep = 0.0f;
	deadRaNitrogenMass_Rep = 0.0f;

	ptnRaMassIncrease_Rep = 0.0f;
	ptnRaNitrogenMassIncrease_Rep = 0.0f;
}

WheatRachis::WheatRachis(float livingFrac) :
	livingFrac(livingFrac), livingFrac_old(livingFrac), livingFrac_ini(livingFrac)
{
	mature = false;             //Z mark the max chaff lengh is reached
	dead = false;
	force_to_death_current_step = false;

	cur_TTd = 0.0f;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	physAge = 0.0f;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	elongAge = 0.0f;
	cumu_growthTdd = 0.0f;		//Z cumu gdd growth age, stressed

	//Z Rachis Length (cm)
	RaLength = 0.0f;				//Z rachis length in cm
	ptnRaLengthIncrease = 0.0f;		//Z potential rachis length increase
	RaLengthIncrease = 0.0f;		//Z actual Rachis length increase

	//Z Rachis Mass g
	srl = 0.01f;
	RaMass = 0.0f;
	deadRaMass = 0.0f;
	RaMassIncrease = 0.0f;
	ptnRaMassIncrease = 0.0f;

	//Z Rachis Nitrogen Mass mg
	RaNitrogenMass = 0.0f;
	RaNitrogenIncrease = 0.0f;
	ptnRaNitrogenMassIncrease = 0.0f;

	//Z Rachis Nitrogne conent and release fraction
	RaNitrogenContent = MAX_N_PCT;
	RaNitrogenReleaseLowBdd = 0.5f;			//Z the min rachis N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	RaNitrogenReleaseMaxPtge = 0.8f;		//Z the rachis N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	//  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined

	RaNitrogenRelease = 0.0f;
	RaNitrogenDead = 0.0f;				//Z stepwise dead portion of rachis nitrogen, will be in the dead plant tissue and not removable
	RaNitrogenReturn = 0.0f;			//Z stepwise returned portion of rachis nitrogen, will supply plant future usage 
	RaNitrogenMassDeadCumu = 0.0f;		//Z cumulative dead portion of rachis nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	N_effect = 1.0;

	//Z ----------- REPRESENTATIVE Rachis ---------------------------------------

	//Z Rachis length
	RaLength_Rep = 0.0f;	//Z this should be a cumulative value, as the living fraction changes, 
	//  it cumulated different amount (based on living fraction)

//Z representative mass values
	RaMass_Rep = 0.0f;
	deadRaMass_Rep = 0.0f;
	RaNitrogenMass_Rep = 0.0f;

	//Z for rachis, there should be no death,
	//  but the living fraction may change (decrease), so there may be still some N release and return 
	RaNitrogenReturn_Rep = 0.0f;
	deadRaNitrogenMass_Rep = 0.0f;

	ptnRaMassIncrease_Rep = 0.0f;
	ptnRaNitrogenMassIncrease_Rep = 0.0f;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Rachis Update Operator Structure

WheatRachisUpdate::WheatRachisUpdate() :
	current_gdd(0.0f), predawn_psi(0.0f), psi_effect_elong(1.0f), growthDuration(0.0f), growthRate(-100.0f), rachisEnd(false) {};

//Z update rachis growth gdd and conditions for this current hourly step
//  stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatRachisUpdate::SetRachisUpdateCondition(float gdd, float psi_predawn, float growthDur, bool rachisEnd)
{
	current_gdd = gdd;
	predawn_psi = psi_predawn;
	//Z should be gdd based. "g/gdd"
	growthDuration = growthDur;

	this->rachisEnd = rachisEnd;

	//Z for the first time calling the function, compute the rachis growing speed
	if (growthRate <= 0.0f) { growthRate = MAXRACHISLENGTH / growthDur; }

	//Z compute water potential effect for expanding and senescence
	//Z some constant first
	const float psi_f = -1.4251f; // -1.0, was -1.4251   later changed to -2.3;
	const float s_f = 0.4258f; // was 0.4258 0.5;
	//Z expanding
	const float psi_threshold_bars_expand = -0.8657f;	//Mpa
	psi_effect_elong = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_expand))));
	psi_effect_elong = min(max(psi_effect_elong, 0.1f), 1.0f);

}

//Z this is the main function for plant rachis growth
//  every time updating, hourly gdd, tmpr, water potential, shade are universal
//                       N stress are local
//                       living fractor is the input variable for each rachis object
/*Z Input tuple for this function includes
*			1. Rachis structure;
*           2. current living fraction (from tiller);
*           3. force_to_death by the plant;
*/
void WheatRachisUpdate::RchsUpdateRun(WheatRachis& wRchs, float lf, bool f2d)
const {
	//Z Simulation of rachis
	//WheatRachis& wRchs = std::get<0>(t);
	//wRchs.livingFrac = std::get<1>(t);
	//wRchs.force_to_death_current_step = std::get<2>(t);

	wRchs.livingFrac = lf;
	wRchs.force_to_death_current_step = f2d;

	//Z reset some critical parameters
	wRchs.ptnRaMassIncrease = 0.0f;
	wRchs.ptnRaNitrogenMassIncrease = 0.0f;

	if (rachisEnd) return;

	//Z compute the N stress for each inidividual leaves, N in mass fraction (mg/g biomass)
	wRchs.N_effect = 2.0f / (1.0f + expf(-2.9f * (max(MIN_N_PCT, wRchs.RaNitrogenContent) - MIN_N_PCT))) - 1.0f;
	wRchs.N_effect = max(min(wRchs.N_effect, 1.0f), 0.1f);

	//Step ZERO Force to death ---------------------------------------------
		//Z this condition is imposed by the parent tiller, stronger than any other conditions
		//Z prevent kill the leaf twice
	if (!wRchs.dead && wRchs.force_to_death_current_step)
	{
		//Z DO NOT change "livingFrac" here
		//  DO NOT reset "Length, NitrogenMass" to 0
		//  That is because we need those values to calculate N return and convert green Chaff length to dead length in the "LivingFractionAdjustment" function
		//  The correctness can be viewed by computing N mass fraction from the output.
		//  The N mass% ranges 0.25% to 4.0%

		//Z determine N release in mg
		//  Release: single Chaff N output due to senecense and dropping, need to divide into remobile(released) portion and dead(within yellow stem) portion
		//  Dead : at that step, reduced N goes to residue and immobilized
		//  Return : at that step, remobilized and can reused for new plant organs
		//  DeadCumu : cumulated N mass in the senecense / dead plant organs


		wRchs.RaNitrogenRelease = wRchs.RaNitrogenMass;
		wRchs.RaNitrogenMass = 0.0f;
		float NitrogenReturnFrac = min(max(wRchs.RaNitrogenContent - wRchs.RaNitrogenReleaseLowBdd, 0.0f) / wRchs.RaNitrogenContent, wRchs.RaNitrogenReleaseMaxPtge);
		float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
		wRchs.RaNitrogenReturn = wRchs.RaNitrogenRelease * NitrogenReturnFrac;
		wRchs.RaNitrogenDead = wRchs.RaNitrogenRelease * NitrogenDeadFrac;
		wRchs.RaNitrogenMassDeadCumu += wRchs.RaNitrogenDead;

		//Z make the internode dead
		wRchs.deadRaMass = wRchs.RaMass;

		//Z set every growing stage, because we do not need to operate this chaff
		wRchs.mature = false;
		wRchs.dead = true;

		goto Rchs_LivingFracAdjustment;
	}

	//Step First -----------------------------------------------------------
		//Z first chaff growth
	if ((!wRchs.dead) && (!wRchs.mature))
	{
		//Z growth is based on length
		float cur_TTd_stressed = min(wRchs.N_effect, psi_effect_elong) * current_gdd;
		wRchs.cumu_growthTdd += cur_TTd_stressed;
		wRchs.elongAge += current_gdd;

		wRchs.ptnRaLengthIncrease = growthRate * cur_TTd_stressed;
		wRchs.RaLengthIncrease = min(wRchs.ptnRaLengthIncrease + wRchs.RaLength, MAXRACHISLENGTH) - wRchs.RaLength;
		wRchs.RaLength += wRchs.ptnRaLengthIncrease;

		//Z convert length growth to biomass and N mass
		wRchs.ptnRaMassIncrease = max(wRchs.RaLength * wRchs.srl - wRchs.RaMass, 0.0f);
		wRchs.ptnRaNitrogenMassIncrease = max(wRchs.RaLength * wRchs.srl * MAX_N_PCT * 10.0f - wRchs.RaNitrogenMass, 0.0f);

		if (wRchs.elongAge >= growthDuration)
		{
			wRchs.mature = true;
		}

		goto Rchs_LivingFracAdjustment;
	}

	//Step Second -----------------------------------------------------------
		//Z second rachis activation
		//  if rachis is matured, no need to update mass like leaf and internode

	//Step Third -----------------------------------------------------------
		//Z third rachis senesce
		//  rachis never dies itself, unless the whole plant dead

	//Step Four Representative Mapping -----------------------------------------------------------
Rchs_LivingFracAdjustment:

	wRchs.ptnRaMassIncrease_Rep = 0.0f;
	wRchs.ptnRaNitrogenMassIncrease_Rep = 0.0f;
	wRchs.RaNitrogenReturn_Rep = 0.0f;

	//Z livingfrac must be smaller than or equal to livingfrac old
	float livingDiff = max(wRchs.livingFrac_old - wRchs.livingFrac, 0.0f);

	//Z rachis representative length
	wRchs.RaLength_Rep += wRchs.RaLengthIncrease * wRchs.livingFrac;

	//-------------------------
	//Z adjust rachis biomass, only dead and living parts
	//  1. dead part, should be a cumulative value
	//  2. note that usually (until the last step), deadRaMass for rachis should be always 0
	wRchs.deadRaMass_Rep += (wRchs.deadRaMass * wRchs.livingFrac + wRchs.RaMass * livingDiff);

	//  total rachis mass
	wRchs.RaMass_Rep = (wRchs.RaMass - wRchs.deadRaMass) * wRchs.livingFrac + wRchs.deadRaMass_Rep;

	//-------------------------
	//Z adjust rachis nitrogen mass and their contribution to nitrogen releasing
	//	1. mass computation is simple, since we assume only dead portion totally recycle nitrogen

	wRchs.RaNitrogenMass_Rep = wRchs.RaNitrogenMass * wRchs.livingFrac;

	//  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
	//	2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
	//	3. NitrogenDecline is just a stepwise variable, and should eventually put into the dead leaf or dead sheath N portion

	wRchs.deadRaNitrogenMass_Rep += wRchs.RaNitrogenDead * wRchs.livingFrac;
	wRchs.RaNitrogenReturn_Rep = wRchs.RaNitrogenReturn * wRchs.livingFrac;
	//Z N release due to the living fraction changes
	float RaNitrogenFree_temp = (wRchs.RaNitrogenMass + wRchs.RaNitrogenRelease) * livingDiff;

	//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
	float NitrogenReturnFrac = min(max(wRchs.RaNitrogenContent - wRchs.RaNitrogenReleaseLowBdd, 0.0f) / wRchs.RaNitrogenContent, wRchs.RaNitrogenReleaseMaxPtge);
	float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
	wRchs.RaNitrogenReturn_Rep += RaNitrogenFree_temp * NitrogenReturnFrac;

	// This part of decline N will be the dead,
	// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller
	// so chaff on dead tiller means "dead chaff" (may be not dropped)
	wRchs.deadRaNitrogenMass_Rep += (RaNitrogenFree_temp * NitrogenDeadFrac);

	//-------------------------
	//Z adjust potential rachis
	//  only living portion can grow

	wRchs.ptnRaMassIncrease_Rep = wRchs.ptnRaMassIncrease * wRchs.livingFrac;
	wRchs.ptnRaNitrogenMassIncrease_Rep = wRchs.ptnRaNitrogenMassIncrease * wRchs.livingFrac;

	//Z finally update the living fraction numbers for this leaf
	wRchs.livingFrac_old = wRchs.livingFrac;

	if (wRchs.dead)
	{
		wRchs.livingFrac_old = 0.0f;
		wRchs.livingFrac = 0.0f;
		wRchs.RaNitrogenMass = 0.0f;
		wRchs.ptnRaMassIncrease_Rep = 0.0f;
		wRchs.ptnRaNitrogenMassIncrease_Rep = 0.0f;
	}
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Rachis Mass Distribution Operator Structure

//Z update rachis mass and 
//  rachis stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatRachisMassAssignment::SetRachisMassAssignmentCondition(float biomassRate, float nitrogenRate)
{
	//Z receive biomass and nitrogen distribution rates for internode over the entire plant based plant-level biomass allocation
	biomassIncomeRate = biomassRate;
	nitrogenIncomeRate = nitrogenRate;
}

//Z rachis mass partition
//  the input parameters are based on the "representative" values of the rachis,
//  but rachis growth and mass should be baesd on "one individual rachis"
//  Thus, we need to make a conversion between "Representative" and "one single rachis" using livingfrac
void WheatRachisMassAssignment::RchsMassRun(WheatRachis& wRchs)
const {
	//Z reset for each time step
	//  store biomass (g) and N (mg) allocated from plant
	wRchs.RaMassIncrease = 0.0f;
	wRchs.RaNitrogenIncrease = 0.0f;

	//Z separate cases explicitly, good for GPU stream control.
	// BIOMASS
	// biomass incremental
	if (wRchs.ptnRaMassIncrease_Rep > 0.0f)
	{
		//Z this "adjustment" variable is necessary and important,
		//  essentially, " / wRchs.livingFrac" is the important term
		//  this convert the mass allocation for a "representative chaff" to "individual chaff"
		wRchs.RaMassIncrease = wRchs.ptnRaMassIncrease_Rep * biomassIncomeRate / wRchs.livingFrac;	// mass in g
	}
	// NITROGEM
	// N incremental
	if (wRchs.ptnRaNitrogenMassIncrease_Rep > 0.0f)
	{
		wRchs.RaNitrogenIncrease = wRchs.ptnRaNitrogenMassIncrease_Rep * nitrogenIncomeRate / wRchs.livingFrac;	// mass in mg
	}

	//Z update chaff bio mass, N mass
	wRchs.RaMass += wRchs.RaMassIncrease;
	wRchs.RaNitrogenMass += wRchs.RaNitrogenIncrease;

	//Z update internnode N content
	//Z OK to over the total internode mass since for a single plant internode, at this time, inernode is still growing in mass
	//  0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
	if (wRchs.RaMass > 0.0f) { wRchs.RaNitrogenContent = max(min(0.10f * wRchs.RaNitrogenMass / wRchs.RaMass, MAX_N_PCT), 0.1f); }
	else { wRchs.RaNitrogenContent = MAX_N_PCT; }
	return;
}