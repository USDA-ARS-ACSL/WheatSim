#include "WheatInternode.h"

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Internode Data Structure

WheatInternode::WheatInternode() :
	InRank(1), TlOrder(0), TlRank(0), mainstem(false), livingFrac(0.0f), livingFrac_old(0.0f), livingFrac_ini(0.0f)
{
	//Z default initialization
	wait_for_initiation = true;
	initiation = false;
	emerge = false;
	ligulation_preL = false;
	ligulation_preL_done = false;
	stop_elongation = false;
	mature = false;
	senesce = false;
	dead = false;
	force_to_death_current_step = false;
	this_internode_elongation = false;

	//Z Growing time/gdd measure
	physAge = 0.0f;
	InPseudoAge = 0.0f;
	InPseudoAge_prev = 0.0f;
	emergeAge = 0.0f;
	emergePseudoAge = 0.0f;
	elongAge = 0.0f;
	seneAge = 0.0f;
	senePseudoAge = 0.0f;
	seneDuration = 0.0f;

	//Z environemtal data
	N_effect = 1.0f;
	slw_effect = 1.0f;

	//Z Hidden Part Growth
	InWholeEmergeDist = INTR_LENGTH_INI;
	InLength_Hidden = INTR_LENGTH_INI;
	Dist2Emerge = 2.0f * INTR_LENGTH_INI;
	RERmax_I = static_cast<float>(InRank <= 6) * RERmax_I_fit[0]
		+ static_cast<float>(InRank == 7) * RERmax_I_fit[1]
		+ static_cast<float>(InRank == 8) * RERmax_I_fit[2]
		+ static_cast<float>(InRank == 9) * RERmax_I_fit[3]
		+ static_cast<float>(InRank == 10) * RERmax_I_fit[4]
		+ static_cast<float>(InRank >= 11) * RERmax_I_fit[5];
	Inconce_Sucrose = Inconce_Sucrose_default;
	Inconce_N = static_cast<float>(InRank == 1 && TlOrder == 0) * Inconce_Aminoacids_default[0]
		+ static_cast<float>((InRank + TlRank) == 2) * Inconce_Aminoacids_default[1]
		+ static_cast<float>((InRank + TlRank) >= 2) * Inconce_Aminoacids_default[2];
	Inconce_Sucrose_mean = 2100.0f;

	//Z Internode Length (cm)
	//Z max internode length based on the internode position, need to compute that at the initialization
	maxInLength_ref = 0.5f + 1.09f * powf(static_cast<float>(InRank - 1), 1.73f) * static_cast<float>(InRank >= 3 && mainstem)
		+ 1.09f * powf(static_cast<float>(InRank), 1.73f) * static_cast<float>(InRank >= 2 && TlRank == 1)
		+ 1.09f * powf(static_cast<float>(InRank), 1.73f) * static_cast<float>(TlRank > 1);
	maxInLength_ref = min(maxInLength_ref, 30.5f);

	InLength = INTR_LENGTH_INI;			//Z internode length cm
	maxInLength_lig = maxInLength_ref;
	maxInLength = maxInLength_ref;
	maxElongationLength = 0.0f;
	maxElongationRate = 0.0f;
	ptnInLengthIncrease = 0.0f;			//Z potnetial internode enlongation, cm
	ptnInLengthDecrease = 0.0f;			//Z potential internode decrease due to death, cm
	matureInLength = 0.1f;				//Z matured internode length, after elongation, cm
	greenInLength = 0.1f;				//Z green or living internode length, after elongation, cm
	seneInLength = 0.0f;				//Z senescent internode length, cm, 

	//Z Internode biomass (g)
	LSIW = static_cast<float>(InRank <= 8) * lsiwfit[0]
		+ static_cast<float>(InRank == 9) * lsiwfit[1]
		+ static_cast<float>(InRank == 10) * lsiwfit[2]
		+ static_cast<float>(InRank == 11) * lsiwfit[3]
		+ static_cast<float>(InRank == 12) * lsiwfit[4]
		+ static_cast<float>(InRank >= 13) * lsiwfit[5];//Z internode specific weight g cm^-1
	ratio_LSIW_LSSW = 2.5f;
	InMass = 0.0f;						//Z internode biomass in g
	InMassIncrease = 0.0f;				//Z biomass assigned to this internode in g
	ptnInMassIncrease = 0.0f;
	deadInMass = 0.0f;					//Z dead internode mass in g
	greenInMass = 0.0f;

	//Z Internode Nitrogen (mg)
	InNitrogenContent = MAX_N_PCT;   	//Z mass fraction in % (e.g., 3.0%)
	InNitrogenMass = 0.0f;				//Z internode nitrogen mass mg
	InNitrogenIncrease = 0.0f;			//Z nitrogen mass assigned to this internode in mg
	ptnInNitrogenMassIncrease = 0.0f;

	//Z Internode Nitrogen (mg): output after organ senecent or dropped
	InNitrogenReleaseLowBdd = 0.0f;		//Z the min N mass content, below that, no N will be released even that organ is dead, assume to be 0.5%
	InNitrogenReleaseMaxPtge = 0.0f;	//Z during plant organ dying, max percent of N releasing

	//Z Individual internode
	//   Release: single internode N output due to senecense and dropping, need to divide into remobile (return) portion and dead (with internode yellow) portion
	//   Dead: at that step, reduced N goes to residue and immobilized
	//   Return: at that step, remobilized and can reused for new plant organs
	//   DeadCumu: cumulated N mass in the senecense/dead plant organs
	InNitrogenRelease = 0.0f;			//Z nitrogen release due to internode decreasing or dying in mg
	InNitrogenDead = 0.0f;				//Z stepwise dead portion of intr nitrogen, will be in the dead plant tissue and not removable
	InNitrogenReturn = 0.0f;			//Z stepwise returned portion of intr nitrogen, will supply plant future usage 
	InNitrogenMassDeadCumu = 0.0f;		//Z cumulated dead portion of nitrogen, in the dead plant tissue and not removable

	//Z ----------- REPRESENTATIVE INTERNODE ---------------------------------------

	greenInLength_Rep = 0.0f;
	seneInLength_Rep = 0.0f;
	InLength_Rep = 0.0f;				//Z internode length include everything, = green + sene, 
	//  and internode will not drop as internode

	InMass_Rep = 0.0f;					//Z representative internode mass, including everything
	greenInMass_Rep = 0.0f;				//Z green internode mass
	deadInMass_Rep = 0.0f;				//Z dead internode mass, internode will not drop, so we end at internode death

	InNitrogenMass_Rep = 0.0f;
	//Z Recall:
	//   Release (single internode): single internode N output due to senecense and dropping, need to divide
	//   Dead (representative internode, later): at that step, reduced N goes to residue and immobilized
	//   Return (representative internode, later): at that step, remobilized and can reused for new plant organs
	//   DeadCumu (representative internode, later): cumulated N mass in the senecense/dead plant organs
	//Z only return and dead N have correlation with plant level mass budget
	//  so we only need to report return (to plant) and deat (fall with dropped internode) N
	//  leaf needs "sene(yellow)+drop=dead"
	//  internode does not need, since internode will not drop
	InNitrogenReturn_Rep = 0.0f;
	deadInNitrogenMass_Rep = 0.0f;

	ptnInMassIncrease_Rep = 0.0f;
	ptnInNitrogenMassIncrease_Rep = 0.0f;
}

WheatInternode::WheatInternode(int rank, int order, int tlRank, bool mainstem, float livingFrac) :
	InRank(rank), TlOrder(order), TlRank(tlRank), mainstem(mainstem), livingFrac(livingFrac), livingFrac_old(livingFrac), livingFrac_ini(livingFrac)
{
	//Z default initialization
	wait_for_initiation = true;
	initiation = false;
	emerge = false;
	ligulation_preL = false;
	ligulation_preL_done = false;
	stop_elongation = false;
	mature = false;
	senesce = false;
	dead = false;
	force_to_death_current_step = false;
	this_internode_elongation = false;

	//Z Growing time/gdd measure
	physAge = 0.0f;
	InPseudoAge = 0.0f;
	InPseudoAge_prev = 0.0f;
	emergeAge = 0.0f;
	emergePseudoAge = 0.0f;
	elongAge = 0.0f;
	seneAge = 0.0f;
	senePseudoAge = 0.0f;
	seneDuration = 0.0f;

	//Z environemtal data
	N_effect = 1.0f;
	slw_effect = 1.0f;

	//Z Hidden Part Growth
	InWholeEmergeDist = INTR_LENGTH_INI;
	InLength_Hidden = INTR_LENGTH_INI;
	Dist2Emerge = 2.0f * INTR_LENGTH_INI;
	RERmax_I = static_cast<float>(InRank <= 6) * RERmax_I_fit[0]
		+ static_cast<float>(InRank == 7) * RERmax_I_fit[1]
		+ static_cast<float>(InRank == 8) * RERmax_I_fit[2]
		+ static_cast<float>(InRank == 9) * RERmax_I_fit[3]
		+ static_cast<float>(InRank == 10) * RERmax_I_fit[4]
		+ static_cast<float>(InRank >= 11) * RERmax_I_fit[5];
	Inconce_Sucrose = Inconce_Sucrose_default;
	Inconce_N = static_cast<float>(InRank == 1 && TlOrder == 0) * Inconce_Aminoacids_default[0]
		+ static_cast<float>((InRank + TlRank) == 2) * Inconce_Aminoacids_default[1]
		+ static_cast<float>((InRank + TlRank) >= 2) * Inconce_Aminoacids_default[2];
	Inconce_Sucrose_mean = 2100.0f;

	//Z Internode Length (cm)
	//Z max internode length based on the internode position, need to compute that at the initialization
	maxInLength_ref = 0.5f + 1.09f * powf(static_cast<float>(InRank - 1), 1.73f) * static_cast<float>(InRank >= 3 && mainstem)
		+ 1.09f * powf(static_cast<float>(InRank), 1.73f) * static_cast<float>(InRank >= 2 && TlRank == 1)
		+ 1.09f * powf(static_cast<float>(InRank), 1.73f) * static_cast<float>(TlRank > 1);
	maxInLength_ref = min(maxInLength_ref, 30.5f);

	InLength = INTR_LENGTH_INI;			//Z internode length cm
	maxInLength_lig = maxInLength_ref;
	maxInLength = maxInLength_ref;
	maxElongationLength = 0.0f;
	maxElongationRate = 0.0f;
	ptnInLengthIncrease = 0.0f;			//Z potnetial internode enlongation, cm
	ptnInLengthDecrease = 0.0f;			//Z potential internode decrease due to death, cm
	matureInLength = 0.1f;				//Z matured internode length, after elongation, cm
	greenInLength = 0.1f;				//Z green or living internode length, after elongation, cm
	seneInLength = 0.0f;				//Z senescent internode length, cm, 

	//Z Internode biomass (g)
	LSIW = static_cast<float>(InRank <= 8) * lsiwfit[0]
		+ static_cast<float>(InRank == 9) * lsiwfit[1]
		+ static_cast<float>(InRank == 10) * lsiwfit[2]
		+ static_cast<float>(InRank == 11) * lsiwfit[3]
		+ static_cast<float>(InRank == 12) * lsiwfit[4]
		+ static_cast<float>(InRank >= 13) * lsiwfit[5];//Z internode specific weight g cm^-1
	ratio_LSIW_LSSW = 2.5f;
	InMass = 0.0f;						//Z internode biomass in g
	InMassIncrease = 0.0f;				//Z biomass assigned to this internode in g
	ptnInMassIncrease = 0.0f;
	deadInMass = 0.0f;					//Z dead internode mass in g
	greenInMass = 0.0f;

	//Z Internode Nitrogen (mg)
	InNitrogenContent = MAX_N_PCT;   	//Z mass fraction in % (e.g., 3.0%)
	InNitrogenMass = 0.0f;				//Z internode nitrogen mass mg
	InNitrogenIncrease = 0.0f;			//Z nitrogen mass assigned to this internode in mg
	ptnInNitrogenMassIncrease = 0.0f;

	//Z Internode Nitrogen (mg): output after organ senecent or dropped
	InNitrogenReleaseLowBdd = 0.0f;		//Z the min N mass content, below that, no N will be released even that organ is dead, assume to be 0.5%
	InNitrogenReleaseMaxPtge = 0.0f;	//Z during plant organ dying, max percent of N releasing

	//Z Individual internode
	//   Release: single internode N output due to senecense and dropping, need to divide into remobile (return) portion and dead (with internode yellow) portion
	//   Dead: at that step, reduced N goes to residue and immobilized
	//   Return: at that step, remobilized and can reused for new plant organs
	//   DeadCumu: cumulated N mass in the senecense/dead plant organs
	InNitrogenRelease = 0.0f;			//Z nitrogen release due to internode decreasing or dying in mg
	InNitrogenDead = 0.0f;				//Z stepwise dead portion of intr nitrogen, will be in the dead plant tissue and not removable
	InNitrogenReturn = 0.0f;			//Z stepwise returned portion of intr nitrogen, will supply plant future usage 
	InNitrogenMassDeadCumu = 0.0f;		//Z cumulated dead portion of nitrogen, in the dead plant tissue and not removable

	//Z ----------- REPRESENTATIVE INTERNODE ---------------------------------------

	greenInLength_Rep = 0.0f;
	seneInLength_Rep = 0.0f;
	InLength_Rep = 0.0f;				//Z internode length include everything, = green + sene, 
	//  and internode will not drop as internode

	InMass_Rep = 0.0f;					//Z representative internode mass, including everything
	greenInMass_Rep = 0.0f;				//Z green internode mass
	deadInMass_Rep = 0.0f;				//Z dead internode mass, internode will not drop, so we end at internode death

	InNitrogenMass_Rep = 0.0f;
	//Z Recall:
	//   Release (single internode): single internode N output due to senecense and dropping, need to divide
	//   Dead (representative internode, later): at that step, reduced N goes to residue and immobilized
	//   Return (representative internode, later): at that step, remobilized and can reused for new plant organs
	//   DeadCumu (representative internode, later): cumulated N mass in the senecense/dead plant organs
	//Z only return and dead N have correlation with plant level mass budget
	//  so we only need to report return (to plant) and deat (fall with dropped internode) N
	//  leaf needs "sene(yellow)+drop=dead"
	//  internode does not need, since internode will not drop
	InNitrogenReturn_Rep = 0.0f;
	deadInNitrogenMass_Rep = 0.0f;

	ptnInMassIncrease_Rep = 0.0f;
	ptnInNitrogenMassIncrease_Rep = 0.0f;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Internode Update Operator Structure

//Z update internode growth gdd and conditions for this current hourly step
//  internode stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatInternodeUpdate::SetInternodeUpdateCondition(float gdd, float tmpr, float teq, float psi_predawn, float shaded, bool elong, bool antepa)
{
	current_gdd = gdd;
	current_tmpr = tmpr;
	current_teq = teq;
	predawn_psi = psi_predawn;
	shade_effect = shaded;
	plant_elongation = elong;		//Z elongation starts
	plant_antepa = antepa;			//Z anthesis, internode should start dying

	//Z compute temperature effect
	//  tmpr effects on growth from K.Paff
	float AlphaGrowth = logf(2.0f) / (logf((Tmax_Intr - Tmin_Intr) / (Topt_Intr - Tmin_Intr)));
	if (current_tmpr<Tmin_Intr || current_tmpr>Tmax_Intr) { tmpr_effect = 0.1f; }
	else {
		float temp1 = powf((current_tmpr - Tmin_Intr), AlphaGrowth);
		float temp2 = powf((Topt_Intr - Tmin_Intr), AlphaGrowth);
		tmpr_effect = (2.0f * temp1 * temp2 - powf(temp1, 2.0f)) / powf(temp2, 2.0f);
	}
	tmpr_effect = min(max(tmpr_effect, 0.1f), 1.0f);

	//Z compute water potential effect for expanding and senescence
	//Z some constant first
	const float psi_f = -1.4251f; // -1.0, was -1.4251   later changed to -2.3;
	const float s_f = 0.4258f; // was 0.4258 0.5;
	//Z expanding
	const float psi_threshold_bars_expand = -0.8657f;	//Mpa
	psi_effect_elong = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_expand))));
	psi_effect_elong = min(max(psi_effect_elong, 0.1f), 1.0f);
	//z senescence
	const float psi_threshold_bars_senesc = -4.0f;		//Mpa
	psi_effect_sene = (1.0f + expf(psi_f * s_f)) / (1.0f + expf(s_f * (psi_f - (predawn_psi - psi_threshold_bars_senesc))));
	psi_effect_sene = min(max(psi_effect_sene, 0.1f), 1.0f);
}

//Z this is the main function for plant internode growth
//  every time updating, hourly gdd, tmpr, water potential, shade are universal
//                       N stress are local
//                       living fractor is the input variable for each leaf object
//Z Input tuple for this function includes
//		1. Internode structure;
//      2. current living fraction (from tiller);
//      3. force_to_death by the plant;
//		4. if previous leaf is emerged
//		5. a tricky parameter, "the ligulation height of the latest emerged leaf (for that tiller) - the internode heights below this leaf"
//						that is the whole length for the leaf to emerge, computer in the tiller class
//      6. this is the elongating internode;

void WheatInternodeUpdate::IntrUpdateRun(WheatInternode& wIntr, float lf, bool f2d, bool ligulation_preL, float ligulationDist, bool elong)
const {
	//Z Simulation of inidividual internode

	//WheatInternode& wIntr = std::get<0>(t);
	//wIntr.livingFrac = std::get<1>(t);
	//wIntr.force_to_death_current_step = std::get<2>(t);
	//wIntr.this_internode_elongation = std::get<3>(t);

	wIntr.livingFrac = lf;
	wIntr.force_to_death_current_step = f2d;
	wIntr.ligulation_preL = ligulation_preL;
	wIntr.InWholeEmergeDist = ligulationDist;
	wIntr.Dist2Emerge = wIntr.InWholeEmergeDist - wIntr.InLength_Hidden;
	wIntr.this_internode_elongation = elong;
	wIntr.physAge += current_gdd;

	//Z reset some critical parameters
	//  do not reset "wIntr.ptnInMassIncrease; wIntr.ptnInNitrogenMassIncrease"
	//  that will be cumulative demand that carrys over the unrealized demand from the previous steps.
	wIntr.ptnInLengthIncrease = 0.0f;
	wIntr.ptnInLengthDecrease = 0.0f;

	wIntr.InNitrogenRelease = 0.0f;
	wIntr.InNitrogenReturn = 0.0f;
	wIntr.InNitrogenDead = 0.0f;

	//Z compute the N stress for each inidividual leaves, N in mass fraction (mg/g biomass)
	wIntr.N_effect = 2.0f / (1.0f + expf(-2.9f * (max(MIN_N_PCT, wIntr.InNitrogenContent) - MIN_N_PCT))) - 1.0f;
	wIntr.N_effect = max(min(wIntr.N_effect, 1.0f), 0.1f);

	float q10fn = powf(2.0f, 0.1f * (current_tmpr - Topt_Intr));

	//Step ZERO Force to death ---------------------------------------------
	//          this condition is imposed by the parent tiller, stronger than any other conditions
	if (!wIntr.dead && wIntr.force_to_death_current_step)	//Z prevent kill the leaf twice
	{
		//Z make the internode senecence to the end, 100%
		if (wIntr.senesce) {
			wIntr.ptnInLengthDecrease = max(0.0f, (wIntr.matureInLength - wIntr.seneInLength));
		}
		else {
			//Z suppose green internode length is non-zero quantity after the elongation, or say senescence is not reached, 
			//  then just use the total length
			wIntr.ptnInLengthDecrease = max(0.0f, wIntr.InLength);
		}
		wIntr.seneInLength += wIntr.ptnInLengthDecrease;
		wIntr.greenInLength = 0.0f;
		wIntr.greenInMass = 0.0f;

		//Z DO NOT change "livingFrac" here
		//  DO NOT set "GreenInLength, InNitrogenMass" to 0
		//  That is because we need those values to calculate N return and convert green internode length to dead length in the "LivingFractionAdjustment" function
		//  The correctness can be viewed by computing N mass fraction from the output.
		//  The N mass% ranges 0.25% to 3.5%

		wIntr.InNitrogenRelease = wIntr.InNitrogenMass;
		wIntr.InNitrogenMass = 0.0f;
		//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
		float NitrogenReturnFrac = min(max(wIntr.InNitrogenContent - wIntr.InNitrogenReleaseLowBdd, 0.0f) / wIntr.InNitrogenContent, wIntr.InNitrogenReleaseMaxPtge);
		float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
		wIntr.InNitrogenReturn = wIntr.InNitrogenRelease * NitrogenReturnFrac;
		wIntr.InNitrogenDead = wIntr.InNitrogenRelease * NitrogenDeadFrac;
		wIntr.InNitrogenMassDeadCumu += wIntr.InNitrogenDead;

		//Z make the internode dead
		wIntr.deadInMass = wIntr.InMass;

		//Z set every growing stage, because we do not need to operate this internode
		wIntr.emerge = false;
		wIntr.mature = false;
		wIntr.senesce = false;
		wIntr.dead = true;
		wIntr.this_internode_elongation = false;

		goto Intr_LivingFracAdjustment;
	}

	//Step One Wait for Initiation ---------------------------------------------
	//     this condition is imposed by the parent tiller, stronger than any other conditions
	if ((!wIntr.dead) && wIntr.wait_for_initiation && (!wIntr.initiation)) 
	{
		wIntr.InPseudoAge += this->current_teq;
		if (wIntr.InPseudoAge >= TIMEDIFF_LF_INTR) {
			wIntr.wait_for_initiation = true;
			wIntr.initiation = true;
			wIntr.InPseudoAge = 0.0f;
		}
		if (plant_antepa) {
			wIntr.wait_for_initiation = false;
			wIntr.initiation = false;
			wIntr.emerge = false;
			wIntr.mature = false;
			wIntr.matureInLength = wIntr.InLength;
			wIntr.greenInMass = wIntr.InMass;

			//Z internode does not grow yet, so it directly dead
			//  no mass or geometry features assigned to this internode
			wIntr.dead = true;
		}
	//	std::cout << wIntr.InPseudoAge << '\n';
		goto Intr_LivingFracAdjustment;
	}

	//Step Two Pre-emerge -----------------------------------------------------------
	//     hidden part growth, growing pattern changes when the previous leaf is emerged
	//     after internode is long enough to be presented from previous pseudostem(sheath), the current internode is emerged
	if ((!wIntr.dead) && wIntr.initiation && (!wIntr.emerge) && (!wIntr.mature)) 
	{
		//Z compute the internode growth rate
		float beta_func_0 = 0.0f;
		float beta_func_t0 = 0.0f;
		float beta_func_t1 = 0.0f;
		beta_func_0 = abs((1.0f + (max(0.0f, (INTR_TE - 0.0f)) / (INTR_TE - INTR_TM)))
			* powf(min(1.0f, (0.0f - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));

		float delta_L = 0.0f;		//Z the length incremental cm 
		float delta_M_IN = 0.0f;	//Z the internode biomass incremental g 
		if (!wIntr.ligulation_preL) {
			delta_L = wIntr.InLength * wIntr.RERmax_I * this->current_teq
				/ (1.0f + RER_KC / wIntr.Inconce_Sucrose) / (1.0f + RER_KN / wIntr.Inconce_N);
			wIntr.InLength += delta_L;
			wIntr.InLength_Hidden = wIntr.InLength;

			//Z compute the internode mass
			//  need to use 0.01f to adjust "cm" to "m" (to use fspm functions)
			//  biomass in g and N mass in mg
			delta_M_IN = IN_ALPHA * IN_BETA * powf(wIntr.InLength * 0.01f, IN_BETA - 1.0f) * delta_L * 0.01f;
			wIntr.ptnInMassIncrease += delta_M_IN;
			wIntr.ptnInNitrogenMassIncrease += (delta_M_IN * MAX_N_PCT * 10.0f);
		}
		if (wIntr.ligulation_preL && (!wIntr.ligulation_preL_done)) {//Z just run once
			wIntr.ligulation_preL_done = true;
			//Z compute the potential internode length
			wIntr.maxInLength_lig = min(wIntr.InLength / beta_func_0, wIntr.maxInLength_ref);
			wIntr.maxInLength = wIntr.maxInLength_lig;
		}
		if (wIntr.ligulation_preL_done) {
			//Z update the internode "time based" life length, after previous leaf ligulation/matured
			//  there the "Beta" function growth starts
			wIntr.InPseudoAge += this->current_teq;
			if (wIntr.InPseudoAge <= INTR_TB) {
				delta_L = wIntr.InLength - beta_func_0 * wIntr.maxInLength;
			}
			if (wIntr.InPseudoAge > INTR_TB && wIntr.InPseudoAge <= INTR_TE) {
				//Z the difference between beta function is the amount of growth
				beta_func_t1 = abs((1.0f + (max(0.0f, (INTR_TE - wIntr.InPseudoAge)) / (INTR_TE - INTR_TM))) *
					powf(min(1.0f, (wIntr.InPseudoAge - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));
				beta_func_t0 = abs((1.0f + (max(0.0f, (INTR_TE - wIntr.InPseudoAge_prev)) / (INTR_TE - INTR_TM))) *
					powf(min(1.0f, (wIntr.InPseudoAge_prev - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));
				delta_L = min(wIntr.maxInLength, wIntr.maxInLength * (beta_func_t1 - beta_func_t0));
				delta_L = delta_L / (1.0f + RER_KC / wIntr.Inconce_Sucrose) / (1.0f + RER_KN / wIntr.Inconce_N);
			}
			if (wIntr.InPseudoAge > INTR_TE) {
				delta_L = 0.0f;
			}

			delta_L = max(0.0f, delta_L);
			wIntr.InLength += delta_L;
			wIntr.InLength_Hidden = wIntr.InLength;
			wIntr.greenInLength = wIntr.InLength;
			wIntr.greenInMass = wIntr.InMass;
			wIntr.InPseudoAge_prev = wIntr.InPseudoAge;

			//Z update the internode specific linear weight
			//  use the default value for now

			//Z update the bioamss demand
			delta_M_IN = wIntr.LSIW * delta_L;
			wIntr.ptnInMassIncrease += delta_M_IN;
			wIntr.ptnInNitrogenMassIncrease += (delta_M_IN * MAX_N_PCT * 10.0f);

			//Z update the internode max length
			wIntr.maxInLength = wIntr.InLength + wIntr.maxInLength_lig * (1.0f - beta_func_t1);
			wIntr.maxInLength = min(wIntr.maxInLength, wIntr.maxInLength_ref);
		}

		//Z decide if the internode can emerge
		if (delta_L > wIntr.Dist2Emerge) {
			wIntr.emerge = true;
		}

		//Z some initial internode can be really small and neven emerge, we still need to encompass maturity condition
		if (abs(wIntr.InLength - wIntr.maxInLength) < 0.01f) {
			wIntr.emerge = true;
			wIntr.mature = true;
			wIntr.matureInLength = wIntr.InLength;
			wIntr.emergePseudoAge = wIntr.InPseudoAge;
			wIntr.seneDuration = wIntr.emergePseudoAge;
		}

		//Z this is a plant level command for internode senesence start
		if (plant_antepa) {
			wIntr.emerge = false;
			wIntr.mature = true;
			wIntr.matureInLength = wIntr.InLength;

			wIntr.emergePseudoAge = wIntr.InPseudoAge;
			wIntr.seneDuration = wIntr.emergePseudoAge;
		}
		goto Intr_LivingFracAdjustment;
	}

	//Step Three -----------------------------------------------------------
	//		internode growth after emergence
	//      not sure if this portion is necessary, but this part should be very slow
	//      becuase most of internode growht should be in the elongation period.
	if ((!wIntr.dead) && (wIntr.emerge) && (!wIntr.this_internode_elongation) && (!wIntr.mature))
	{
		//Z compute the internode growth rate
		float beta_func_0 = 0.0f;
		float beta_func_t0 = 0.0f;
		float beta_func_t1 = 0.0f;
		beta_func_0 = abs((1.0f + (max(0.0f, (INTR_TE - 0.0f)) / (INTR_TE - INTR_TM)))
			* powf(min(1.0f, (0.0f - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));

		float delta_L = 0.0f;		//Z the length incremental cm 
		float delta_M_IN = 0.0f;	//Z the internode biomass incremental g 
	
		//Z update the internode "time based" life length, after this internode is emerged
		wIntr.InPseudoAge += this->current_teq;
	
		if (wIntr.InPseudoAge <= INTR_TB) {
			delta_L = wIntr.InLength - beta_func_0 * wIntr.maxInLength;
		}
		if (wIntr.InPseudoAge > INTR_TB && wIntr.InPseudoAge <= INTR_TE) {
			//Z the difference between beta function is the amount of growth
			beta_func_t1 = abs((1.0f + (max(0.0f, (INTR_TE - wIntr.InPseudoAge)) / (INTR_TE - INTR_TM))) *
				powf(min(1.0f, (wIntr.InPseudoAge - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));
			beta_func_t0 = abs((1.0f + (max(0.0f, (INTR_TE - wIntr.InPseudoAge_prev)) / (INTR_TE - INTR_TM))) *
				powf(min(1.0f, (wIntr.InPseudoAge_prev - INTR_TB) / (INTR_TE - INTR_TB)), (INTR_TE - INTR_TB) / (INTR_TE - INTR_TM)));
			delta_L = min(wIntr.maxInLength, wIntr.maxInLength * (beta_func_t1 - beta_func_t0));
			delta_L = delta_L / (1.0f + RER_KC / wIntr.Inconce_Sucrose) / (1.0f + RER_KN / wIntr.Inconce_N);
		}
		if (wIntr.InPseudoAge > INTR_TE) {
			delta_L = 0.0f;
		}

		delta_L = max(0.0f, delta_L);
		wIntr.InLength += delta_L;
		wIntr.greenInLength = wIntr.InLength;
		wIntr.greenInMass = wIntr.InMass;
		wIntr.InPseudoAge_prev = wIntr.InPseudoAge;

		//Z update the internode specific linear weight
		//  use the default value for now
		
		//Z update the bioamss demand
		delta_M_IN = wIntr.LSIW * delta_L;
		wIntr.ptnInMassIncrease += delta_M_IN;
		wIntr.ptnInNitrogenMassIncrease += (delta_M_IN * MAX_N_PCT * 10.0f);

		//Z update the internode max length
		wIntr.maxInLength = wIntr.InLength + wIntr.maxInLength_lig * (1.0f - beta_func_t1);
		wIntr.maxInLength = min(wIntr.maxInLength, wIntr.maxInLength_ref);

		//Z update max elongtation scale and rate
		//  the one grow during elongation
		wIntr.maxElongationLength = max((wIntr.maxInLength_ref - wIntr.InLength), 0.0f);
		wIntr.maxElongationRate = wIntr.maxElongationLength / PHYLLOCHRON;

		//Z decide if the internode can mature
		if (abs(wIntr.InLength - wIntr.maxInLength) < 0.01f) {
			wIntr.mature = true;
			wIntr.matureInLength = wIntr.InLength;
			wIntr.emergePseudoAge = wIntr.InPseudoAge;
			wIntr.seneDuration = wIntr.emergePseudoAge;
		}

		//Z this is a plant level command for internode senesence start
		if (plant_antepa) {
			wIntr.mature = true;
			wIntr.matureInLength = wIntr.InLength;
			wIntr.emergePseudoAge = wIntr.InPseudoAge;
			wIntr.seneDuration = wIntr.emergePseudoAge;
		}
		goto Intr_LivingFracAdjustment;
	}

	//Step FOUR --------------------------------------------------
	//		internode elongation expansion
	//		Note that in maizsim, elongation is different from growth, 
	//      growth: cell number and size grow, and biomass infiltrate,
	//		elongation: cell number will not grow, but size will elongate rapidly, and biomass will also inserted.
	if ((!wIntr.dead) && wIntr.this_internode_elongation && (!wIntr.stop_elongation))
	{
		if (wIntr.emergePseudoAge <= 0.0f) 
		{
			wIntr.emergePseudoAge = wIntr.InPseudoAge;
			wIntr.seneDuration = wIntr.emergePseudoAge;
		}
		//Z beta function for current "plant growth time"
		wIntr.elongAge += current_gdd;

		//Z potential incremental for the internode length
		float maxGrowRate_intr = tmpr_effect * wIntr.maxElongationRate;
		wIntr.ptnInLengthIncrease = maxGrowRate_intr * current_gdd;

		//Z stressed incremental for the internode length
		wIntr.ptnInLengthIncrease *= (min(wIntr.N_effect, psi_effect_elong) * shade_effect);
		wIntr.InLength += wIntr.ptnInLengthIncrease;
		wIntr.greenInLength = wIntr.InLength;
		wIntr.greenInMass = wIntr.InMass;

		//Z potential biomass increase internode (g), nitrogen increase internode (mg), demand
		//  internode may not necessarily grow at the stage
		//  but can still fill biomass until the slw number is reached
		wIntr.ptnInMassIncrease += wIntr.ptnInLengthIncrease * wIntr.LSIW;
		wIntr.ptnInNitrogenMassIncrease += wIntr.ptnInLengthIncrease * wIntr.LSIW * MAX_N_PCT * 10.0f;

		if ((wIntr.InLength >= wIntr.maxElongationLength) || (wIntr.elongAge > PHYLLOCHRON) || (plant_antepa))
		{
			wIntr.stop_elongation = true;
			wIntr.mature = true;
			wIntr.matureInLength = wIntr.InLength;
		}
		goto Intr_LivingFracAdjustment;
	}

	//Step FIVE -----------------------------------------------------------
	//      the internode just stay, do nothing
	//		Internode not in elongation should also update its mass
	//		plant_antepa is the time that all internode must mature and start sene
	if ((!wIntr.dead) && wIntr.mature && (!plant_antepa))
	{
		//Z ensure the internode length will not change
		wIntr.greenInLength = wIntr.InLength;
		wIntr.greenInMass = wIntr.InMass;

		goto Intr_LivingFracAdjustment;
	}

	//Step SIX -----------------------------------------------------------
	//		internode senesce
	//		internode senesce should be at the end of the wheat growing season, after the flowers are completed but the grain is still milking
	//		all the intenrode should die at the same time, but we appled a very slow dying speed, i.e., 500 gdd
	//		otherwise bottom (old) internodes die first, and then cut the connection to the top region of that stem, which is silly even for a plant
	if ((!wIntr.dead) && plant_antepa)
	{
		//Z during sene, the internode will not accept new biomass and N
		wIntr.ptnInMassIncrease = 0.0f;
		wIntr.ptnInNitrogenMassIncrease = 0.0f;

		//Z senesce start
		wIntr.senesce = true;
		wIntr.senePseudoAge += this->current_teq;
		wIntr.seneDuration = max(0.0f, wIntr.seneDuration - 0.5f * this->current_teq * (1.0f - min(psi_effect_sene, wIntr.N_effect)));
		wIntr.senePseudoAge = min(wIntr.senePseudoAge, wIntr.seneDuration);

		wIntr.seneInLength = max(0.0f, wIntr.matureInLength * (1.0f + 2.0f * (1.0f - wIntr.senePseudoAge / wIntr.seneDuration)) * powf(wIntr.senePseudoAge / wIntr.seneDuration, 2.0f));

		if (wIntr.seneInLength > wIntr.matureInLength || wIntr.senePseudoAge >= wIntr.seneDuration)
		{
			wIntr.seneInLength = wIntr.matureInLength;
			wIntr.deadInMass = wIntr.InMass;
			wIntr.dead = true;
			wIntr.livingFrac = 0.0f;
		}

		//Z manage internode mass and N mass based on internode decreasing (actually internode is "yellowing", since internode never drop)
		wIntr.ptnInLengthDecrease = max(0.0f, wIntr.seneInLength + wIntr.greenInLength - wIntr.matureInLength);
		//Z this is the relative sene fraction relative to the current living portion of the plant organ, max to one
		float relative_sene_fraction = min(wIntr.ptnInLengthDecrease / wIntr.greenInLength, 1.0f);
		//Z then update the green/living internode portion
		wIntr.greenInLength = max(0.0f, wIntr.matureInLength - wIntr.seneInLength);

		//Z determine N release in mg
		//  Release: single internode N output due to senecense and dropping, need to divide into remobile(released) portion and dead portion (stay in the yellow part of the intenrode)
		//  Dead : at that step, reduced N goes to residue and immobilized
		//  Return : at that step, remobilized and can reused for new plant organs
		//  DeadCumu : cumulated N mass in the senecense / dead plant organs

		wIntr.InNitrogenRelease = relative_sene_fraction * wIntr.InNitrogenMass;
		wIntr.InNitrogenMass = wIntr.InNitrogenMass - wIntr.InNitrogenRelease;
		//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
		float NitrogenReturnFrac = min(max(wIntr.InNitrogenContent - wIntr.InNitrogenReleaseLowBdd, 0.0f) / wIntr.InNitrogenContent, wIntr.InNitrogenReleaseMaxPtge);
		float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
		wIntr.InNitrogenReturn = wIntr.InNitrogenRelease * NitrogenReturnFrac;
		wIntr.InNitrogenDead = wIntr.InNitrogenRelease * NitrogenDeadFrac;
		wIntr.InNitrogenMassDeadCumu += wIntr.InNitrogenDead;

		//Z determine Internode biomass in g
		wIntr.greenInMass = (1.0f - relative_sene_fraction) * wIntr.greenInMass;

		goto Intr_LivingFracAdjustment;
	}

	//Step Four Representative Mapping -----------------------------------------------------------
	//Z Simulation of representative internode, adjusted by the livingFrac
Intr_LivingFracAdjustment:

	wIntr.ptnInMassIncrease_Rep = 0.0f;
	wIntr.ptnInNitrogenMassIncrease_Rep = 0.0f;
	wIntr.InNitrogenReturn_Rep = 0.0f;

	//Z livingfrac must be smaller than or equal to livingfrac old
	float livingDiff = max(wIntr.livingFrac_old - wIntr.livingFrac, 0.0f);

	//Z adjust living internode geometry, mass and potential mass required for internode growth
	//green or living part
	wIntr.greenInLength_Rep = wIntr.greenInLength * wIntr.livingFrac;
	wIntr.greenInMass_Rep = wIntr.greenInMass * wIntr.livingFrac;

	//  senecent part
	//  senecent internode portion has to parts
	//	1. senecent of living internode
	//	2. dead of the whole internode
	//	This is a cumulative value so must be initialized to 0
	wIntr.seneInLength_Rep += (wIntr.ptnInLengthDecrease * wIntr.livingFrac + (wIntr.greenInLength + wIntr.ptnInLengthDecrease) * livingDiff);

	//  comparing to leaf, no drop part for internode since internode will never drop
	//  total internode length, 
	wIntr.InLength_Rep = wIntr.greenInLength_Rep + wIntr.seneInLength_Rep;

	//-------------------------
	//Z adjust internode biomass
	//  1. dead part, should be a cumulative value
	//  2. note that usually (until the last step), deadInMass for individual internode should be always 0
	wIntr.deadInMass_Rep += (wIntr.deadInMass * wIntr.livingFrac + wIntr.InMass * livingDiff);

	//  total internode mass
	wIntr.InMass_Rep = (wIntr.InMass - wIntr.deadInMass) * wIntr.livingFrac + wIntr.deadInMass_Rep;

	//-------------------------
	//Z adjust internode nitrogen mass and their contribution to nitrogen releasing
	//	1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
	wIntr.InNitrogenMass_Rep = wIntr.InNitrogenMass * wIntr.livingFrac;

	//  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
	//	2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
	//	3. NitrogenDecline is just a stepwise variable, and should eventually put into the dead leaf or dead sheath N portion
	wIntr.deadInNitrogenMass_Rep += wIntr.InNitrogenDead * wIntr.livingFrac;
	wIntr.InNitrogenReturn_Rep = wIntr.InNitrogenReturn * wIntr.livingFrac;
	//Z N release due to the living fraction changes
	float InNitrogenFree_temp = (wIntr.InNitrogenMass + wIntr.InNitrogenRelease) * livingDiff;

	//Z redistribute (divide the nitrogen released) into dead, return, then cumulate the dead portion
	float NitrogenReturnFrac = min(max(wIntr.InNitrogenContent - wIntr.InNitrogenReleaseLowBdd, 0.0f) / wIntr.InNitrogenContent, wIntr.InNitrogenReleaseMaxPtge);
	float NitrogenDeadFrac = 1.0f - NitrogenReturnFrac;
	wIntr.InNitrogenReturn_Rep += InNitrogenFree_temp * NitrogenReturnFrac;

	// This part of decline N will be the dead,
	// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller 
	// so internode on dead tiller means "dead internode" (may be not dropped)
	//Z in another word: "leaf N should be divided by dead=sene+dropped", while internode N is just "dropped"
	wIntr.deadInNitrogenMass_Rep += (InNitrogenFree_temp * NitrogenDeadFrac);

	//-------------------------
	//Z adjust potential internode
	//  only living portion can grow

	wIntr.ptnInMassIncrease_Rep = wIntr.ptnInMassIncrease * wIntr.livingFrac;
	wIntr.ptnInNitrogenMassIncrease_Rep = wIntr.ptnInNitrogenMassIncrease * wIntr.livingFrac;

	//Z finally update the living fraction numbers for this leaf
	wIntr.livingFrac_old = wIntr.livingFrac;

	if (wIntr.dead)
	{
		wIntr.livingFrac_old = 0.0f;
		wIntr.livingFrac = 0.0f;

		wIntr.greenInLength = 0.0f;
		wIntr.greenInMass = 0.0f;
		wIntr.InNitrogenMass = 0.0f;

		wIntr.ptnInMassIncrease_Rep = 0.0f;
		wIntr.ptnInNitrogenMassIncrease_Rep = 0.0f;
	}

	if (wIntr.livingFrac == 0.0f)
	{
		wIntr.dead = true;
		wIntr.livingFrac_old = 0.0f;

		wIntr.greenInLength = 0.0f;
		wIntr.greenInMass = 0.0f;
		wIntr.InNitrogenMass = 0.0f;

		wIntr.ptnInMassIncrease_Rep = 0.0f;
		wIntr.ptnInNitrogenMassIncrease_Rep = 0.0f;
	}
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//Z Internode Mass Distribution Operator Structure

//Z update internode mass and 
//  internode stress ranges 0-1.0,
//  but numerically, we assume 0.1-1.0
void WheatInternodeMassAssignment::SetInternodeMassAssignmentCondition(float biomassRate, float nitrogenRate)
{
	//Z receive biomass and nitrogen distribution rates for internode over the entire plant based plant-level biomass allocation
	biomassIncomeRate = biomassRate;
	nitrogenIncomeRate = nitrogenRate;
}

//Z internode mass partition
//  the input parameters are based on the "representative" values of the internode,
//  but internode growth and mass should be baesd on "one individual internode"
//  Thus, we need to make a conversion between "Representative" and "one single internode" using livingfrac
void WheatInternodeMassAssignment::IntrMassRun(WheatInternode& wIntr)
const {
	//Z reset for each time step
	//  store biomass (g) and N (mg) allocated from plant
	wIntr.InMassIncrease = 0.0f;
	wIntr.InNitrogenIncrease = 0.0f;

	//Z separate cases explicitly, good for GPU stream control.
	// BIOMASS
	// biomass incremental
	if (wIntr.ptnInMassIncrease_Rep > 0.0f)
	{
		//Z this "adjustment" variable is necessary and important,
		//  essentially, " / wIntr.livingFrac" is the important term
		//  this convert the mass allocation for a "representative internode" to "individual internode"
		wIntr.InMassIncrease = wIntr.ptnInMassIncrease_Rep * biomassIncomeRate / wIntr.livingFrac;	// mass in g
	}
	// NITROGEM
	// N incremental
	if (wIntr.ptnInNitrogenMassIncrease_Rep > 0.0f)
	{
		wIntr.InNitrogenIncrease = wIntr.ptnInNitrogenMassIncrease_Rep * nitrogenIncomeRate / wIntr.livingFrac;	// mass in mg
	}

	//Z update internnode bio mass, N mass
	wIntr.InMass += wIntr.InMassIncrease;
	wIntr.greenInMass += wIntr.InMassIncrease;
	wIntr.InNitrogenMass += wIntr.InNitrogenIncrease;

	//Z update leaf and sheath ptn mass and N mass, i.e., the biomass and N demand
	wIntr.ptnInMassIncrease -= wIntr.InMassIncrease;
	wIntr.ptnInNitrogenMassIncrease -= wIntr.InNitrogenIncrease;

	//Z update internnode N content
	//Z OK to over the total internode mass since for a single plant internode, at this time, inernode is still growing in mass
	//  0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
	if (wIntr.InMass > 0.0f) { wIntr.InNitrogenContent = max(min(0.10f * wIntr.InNitrogenMass / wIntr.InMass, MAX_N_PCT), 0.1f); }
	else { wIntr.InNitrogenContent = MAX_N_PCT; }
	return;
}