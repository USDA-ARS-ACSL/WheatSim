#include "stdafx.h"
#include "WheatTiller.h"


//Z constructor, need to identify if it is the mainstem
WheatTiller::WheatTiller(int n, int o, int cr, bool m, WheatDevelopment* dv, WheatOrganDataFrame* ddf, float tillerLivingFrac)
{
	//-----------------------------------------------------
	//Z Tiller Position Initiation
		//Z initialize the growing stage
	rank = n;
	order = o;
	cumurank = cr;
	mainstem = m;
	develop = dv;
	organDataFrame = ddf;
	mainstemInitializaiton = false;	//Z a special case for mainstem, immediately grow two leaves
	death_2_finalize = false;

	gddpBranch = develop->get_TTd_Plant();
	physAge = 0.0f;
	TlPseudoAge = 0.0f;
	living = true;						//Z the tiller is living, reserve for counting tiller death, living=1, death=0
	livingFrac = tillerLivingFrac;		//Z DIFFERENT FROM LEAF,INTERNODE, leaf and internode livingfrac follows the tiller, while tiller and subtiller livingfrac management independently
	//  The initial subtiller livingfrac follows the tiller livingfrac, after initialization, subtiller will be on its own, EVEN the parent tiller is dead.
	//  This setting is meant to mimic "perenial cases" where in the next year, all tillers are new subtillers 
	//  For some crops, tiller can also has its own roots, so parent tiller death does not imply subtiller death

	livingFrac_ini = tillerLivingFrac;
	force_to_death_current_step = false;    //Z marker for tiller death due to external/ambient reasons, force to die

	for (int ii = 0; ii < MAXLEAFNUM; ii++) 
	{ 
		TillerLfPhyllchron[ii] = static_cast<float>(ii + 1) * PHYLLOCHRON;
		TillerLfPlastoChrone[ii] = static_cast<float>(ii + 1) * PLASTOCHRONE;
		cumuIntrLength[ii] = 0.0f;
		cumuLigulationHeight[ii] = 0.0f;
	}
	FlagLfPhyllchron_Init = false;
	FlagLfPhyllchron_Cplt = TillerLfPhyllchron[MAXLEAFNUM - 1];

	//-----------------------------------------------------
	//Z Growing Stage
	//  jointing stage is a killer for the tillers
	//  for the whole plant, such that each tiller follow that stage
	singleRidge = develop->is_singleRidge();	//Z the singleRidge stage is reached for the whole plant
	born_after_singleRidge = singleRidge;

	terSpikeInit = develop->is_terSpikeGrow();
	flowerPrimInit = develop->is_flowerPrimInit();
	anthesisStart = develop->is_startAnthesis();
	anthesisEnd = develop->is_endAnthesis();
	flowerPrimInit_already = false;
	anthesisStart_already = false;

	jointing = develop->is_startJoint();
	jointing_already = jointing;			//Z combine with "jointing", s.t. tiller death operation on jointing dates
	jointing_kill = false;

	//  elongation ridge stage is 
	elongationStart = false;					//Z the enlongation is initiated
	elongationStart_already = false;			//Z if elongation already occurs, combine with "elongationStart" to trace the first entry of elongation
	//Z between these two, we determine the elongation index

	//  jointing stage
	jointing = develop->is_startJoint();

	//Z marker for the elongation internode ranks
	elongation_first = -1;
	elongation_last = -1;
	rank_2_elongation = -1;

	//Z marker and parameters for spikelet/flower growth
	phyAge_SpikInit = 0.0f;
	phyAge_FlwFert = 0.0f;
	gddSpike = PHYLLOCHRON / SPIKEPERPHYLLOCHRON;
	SpikAdd = 0;	//Z this step add one spikelet, should only be 0 or 1
	FlwFertAdd = 0;	//Z this step fertilize one spikelet's flower, based on basel flower number
	spiint = 0;		//Z initial rank of spikelet that can have flower initiation
	spifrt = 0;		//Z first fertilization spikelet
	spitop = 0;		//Z top rank of the spikelet that can have flower initiation
	spibot = 0;		//Z bottom rank of the spikelet that can have flower initiation
	ferint = 0;		//Z initial fertilization spikelet
	ferbot = 0;		//Z bottom fertilization spikelet
	fertop = 0;		//Z top fertilization spikelet

	//-----------------------------------------------------
	//Z Environmental conditions	
		//Z initialize tiller temperature
	Cur_TTd = 0.0f;
	Cur_TEq = 0.0f;
	StressTiller = 1.0f;
	StressTillerMin = 1.0f;
	TTd_since_singleRidge = 0.0f;
	TTd_since_jointing = 0.0f;
	N_effect = 1.0f;
	water_effect = 1.0f;
	shade_effect = 1.0f;
	tmpr_effect = 1.0f;
	tmpr_effect_terminal = 1.0f;
	Cold_Time = 0.0f;

	//Z integer counts for the tillers
	//  tiller is just a "box", leaf and internode are the real growing parts
	ChaffNum = 0;
	LeafNum = 0;					//Z Leaf number = green leaf number + dropped leaf number
	EmergeLfNum = 0;
	GreenLfNum = 0;
	DropLfNum = 0;
	InterNodeNum = 0;
	SpikNum = 0;
	SubTillerNum = 0;
	RachisNum = 0;

	for (int ii = 0; ii < MAXLEAFNUM; ii++) { SubLeafIdx[ii] = -1; }
	for (int ii = 0; ii < MAXINTRNUM; ii++) { SubInterNodeIdx[ii] = -1; }
	for (int ii = 0; ii < MAXSPIKNUM; ii++) { SubSpikIdx[ii] = -1; }
	SubChaffIdx = -1;
	SubRachisIdx = -1;

	cout << "Tiller Emerging: Order = " << order << ", Rank = " << rank << "\n";
}

//Z destructor,
//  in all the models, death does not mean destruction,
//    dead organ still holds informations until the entire program ends

//  Therefore, this destructor will not involve mass transfer issues
//   Recall for any recursive process, only call it from the root element, aka the mainstem

//   In WheatSim, since we use smart pointer, the destructor can be exetremely simple
WheatTiller::~WheatTiller() {}

//Z special case of mainstem initializaiton between germinationand emergence
//	apsim initialized 2 leaves, in the function
//	a. initial "FIRST" and "SECOND" leaves and internodes
//	b. BUT DO NOT revise the TillerLfPhyllchron order,
//	   such that the plant still wait for 2 Phyllchron to get the first tiller
//
//Z why this is an isolated function separate from the constructor?
//	because constructor for mainstem is in the constructor of plant
//	and this function is called between germination and emergence
//
//Z because this is the mainstem, "livingFrac=1.0" for all its leaves and internodes initially

void WheatTiller::MainstemInitialize()
{
	//Z this function only called once for mainstem
	//  but we still put if-statment for protection
	if (mainstem && (!mainstemInitializaiton))
	{
		mainstemInitializaiton = true;

		//Z add new leaves, parameter: 
		//	index of leaf on the tiller: start from 1; 
		//  order of the tiller: start from 0 (mainstem)
		//  rank of the tiller: start from 0 (mainstem)
		//  bool if tiller is mainstem; 
		//  initial livingFrac
		SubLeafIdx[LeafNum] = organDataFrame->h_gLeafNum;
		LeafNum += 1;
		GreenLfNum += 1;
		organDataFrame->h_gLeafNum += 1;
		WheatLeaf wlf(LeafNum, 0, 0, true, 1.0f);
		organDataFrame->h_gLeaf.push_back(wlf);
		organDataFrame->h_gLeafLivingFrac.push_back(1.0f);
		organDataFrame->h_gLeafForce2Death.push_back(false);
		organDataFrame->h_gLeafPreLemerge.push_back(false);
		organDataFrame->h_gLeafLigulationDist.push_back(1.0f);
		
		//Z add new internode, parameters: 
		//  index of internode on the tiller: start from 1; 
		//  order of the tiller: start from 0 (mainstem);
		//  rank of the tiller: start from 0 (mainstem)
		//  bool if tiller is mainstem; 
		//  initial livingFrac
		SubInterNodeIdx[InterNodeNum] = organDataFrame->h_gIntrNum;
		InterNodeNum += 1;
		organDataFrame->h_gIntrNum += 1;
		WheatInternode wit(InterNodeNum, 0, 0, true, 1.0f);
		organDataFrame->h_gIntr.push_back(wit);
		organDataFrame->h_gIntrLivingFrac.push_back(1.0f);
		organDataFrame->h_gIntrForce2Death.push_back(false);
		organDataFrame->h_gIntrElongation.push_back(false);
		organDataFrame->h_gIntrPreLigulation.push_back(false);
		organDataFrame->h_gIntrLigulationDist.push_back(1.0f);

		for (int ii = 0; ii < MAXLEAFNUM; ii++)
		{
			TillerLfPhyllchron[ii] = static_cast<float>(ii) * PHYLLOCHRON;
			TillerLfPlastoChrone[ii] = static_cast<float>(ii) * PLASTOCHRONE;
		}
		FlagLfPhyllchron_Init = false;
		FlagLfPhyllchron_Cplt = TillerLfPhyllchron[MAXLEAFNUM - 1];
	}
}

// *************** GROWTH AND DEVELOP *******************************************
// As object-oriented programming, one tiller should only care about itself, and set up update for leaves, internodes, flowers, and kernels
// As data-oriented programming, actual update of leaf and internode will be on GPU 
// multiple tiller recursive operation should be done in a upper class -> WheatPlant

//Z NEW leaf and internode only on this tiller
//     include new leaf, internode, and subtillers
//     because new organs will emerge, the tiller has to be living
// 
//Z this is the first function to call when tiller is operatored
void WheatTiller::WheatTillerSingleMorph(void)
{
	if (living)
	{
		//Z get the TTd for this time step and adding to the age of this tiller
		//  TTd is the thermal time for this 1 hour period
		Cur_TTd = develop->get_TTdCur_dssat();
		Cur_TEq = develop->get_TimeTeqCur_fspm();
		elongationStart = develop->is_startEnlongation();
		physAge += Cur_TTd;
		TlPseudoAge += Cur_TEq;

		//Z spikelet grows between single ridge stage and the spike growth termination stage
		singleRidge = develop->is_singleRidge();
		terSpikeInit = develop->is_terSpikeGrow();

		//Z these variables are used to judge where to initialize, abort, and fertilizer flowers
		//Z it will be only 0 or 1;
		SpikAdd = 0;
		FlwFertAdd = 0;

		//Z new leaf, tiller, and internode should be initiated if the following condition is satisfied
		if (TlPseudoAge >= TillerLfPlastoChrone[LeafNum])
		{
			// emerge leaf
			SubLeafIdx[LeafNum] = organDataFrame->h_gLeafNum;
			LeafNum += 1;
			GreenLfNum += 1;
			organDataFrame->h_gLeafNum += 1;
			WheatLeaf wlf(LeafNum, this->order, this->cumurank, this->mainstem, this->livingFrac);
			organDataFrame->h_gLeaf.push_back(wlf);
			organDataFrame->h_gLeafLivingFrac.push_back(this->livingFrac);
			organDataFrame->h_gLeafForce2Death.push_back(false);
			organDataFrame->h_gLeafPreLemerge.push_back(false);
			organDataFrame->h_gLeafLigulationDist.push_back(1.0f);
			
			// second emerge internode
			SubInterNodeIdx[InterNodeNum] = organDataFrame->h_gIntrNum;
			organDataFrame->h_gIntrNum += 1;
			InterNodeNum += 1;
			WheatInternode witr(this->InterNodeNum, this->order, this->cumurank, this->mainstem, this->livingFrac);
			organDataFrame->h_gIntr.push_back(witr);
			organDataFrame->h_gIntrLivingFrac.push_back(this->livingFrac);
			organDataFrame->h_gIntrForce2Death.push_back(false);
			organDataFrame->h_gIntrElongation.push_back(false);
			organDataFrame->h_gIntrPreLigulation.push_back(false);
			organDataFrame->h_gIntrLigulationDist.push_back(1.0f);
		}

		//Z adding tillers, need emerged leaf number
		EmergeLfNum = 0;
		for (int ii = 0; ii < LeafNum; ii++)
		{
			//Z leaf area update for existing leaves
			int IdxLeaf = SubLeafIdx[ii];
			if (IdxLeaf != -1)
			{
				EmergeLfNum += static_cast<int>(organDataFrame->h_gLeaf[IdxLeaf].emerge);
			}
		}

		//Z second new sub-tiller
		// if leaf appears, then tiller at pre-specified position will also appear
		// first argument means "mainstem=false", that is correct because subtiller cannot be the mainstem
		//		1. if elongation starts, no new tillers
		//		2. if tiller killed during singleRidge, no subtiller from it but existing subtillers will continue
		if ((!elongationStart) && (!jointing_kill))
		{
			if (mainstem && (EmergeLfNum >= 3))
			{
				if (SubTillerNum <= EmergeLfNum - 3)
				{
					//Z recall the input parameters, bool if sub - tiller is the mainstem(of course not)
					//								 develop obj
					//								 initial sub-tiller livingfrac follow the tiller livingfrac
					float tillerStress = this->livingFrac;
					tillerStress = std::max(std::min(tillerStress, N_effect), 0.4f);
					SubTiller.push_back(make_unique<WheatTiller>(this->EmergeLfNum - 2, this->order + 1, this->EmergeLfNum - 2, false, develop, organDataFrame, tillerStress));
					SubTillerNum += 1;
				}
			}
			if (!mainstem && (EmergeLfNum >= 2))
			{
				if (SubTillerNum <= EmergeLfNum - 2)
				{
					//Z recall the input parameters, bool if sub - tiller is the mainstem(of course not)
					//								 develop obj
					//								 initial sub-tiller livingfrac follow the tiller livingfrac
					float tillerStress = this->livingFrac;
					tillerStress = std::max(std::min(tillerStress, N_effect), 0.4f);
					SubTiller.push_back(make_unique<WheatTiller>(this->EmergeLfNum - 2, this->order + 1, this->cumurank + this->EmergeLfNum - 2, false, develop, organDataFrame, tillerStress));
					SubTillerNum += 1;
				}
			}
		}

		//Z fourth spikelet primordia, flower and kernel iniation
		//  because we do hourly simulation, we can add spikelet one by one even the TTd interval is just ~10 gdd
		if ((singleRidge) && (!terSpikeInit))
		{
			phyAge_SpikInit += Cur_TTd;
			//Z at this step, need to add one "addSpike" spike
			if ((phyAge_SpikInit > gddSpike) && (SpikNum < MAXSPIKNUM))
			{
				SubSpikIdx[SpikNum] = organDataFrame->h_gSpikNum;
				SpikNum += 1;
				SpikAdd = 1;
				organDataFrame->h_gSpikNum += 1;
				WheatSpikelet wspk(this->SpikNum, this->livingFrac);
				organDataFrame->h_gSpik.push_back(wspk);
				organDataFrame->h_gSpikLivingFrac.push_back(this->livingFrac);
				organDataFrame->h_gSpikForce2Death.push_back(false);
				organDataFrame->h_gSpikFlowerInit.push_back(false);
				organDataFrame->h_gSpikFlowerAbort.push_back(false);
				organDataFrame->h_gSpikFlowerFert.push_back(false);
				
				//Z reduce the phyAge for spikelet initiation and reset-accumulation for the next spikelet
				phyAge_SpikInit -= gddSpike;
			}
		}

		//Z fifth chaff
		//  start from chaffStart and End with chaffEnd
		if (develop->is_chaffStart() && ChaffNum == 0)
		{
			SubChaffIdx = organDataFrame->h_gChafNum;
			organDataFrame->h_gChafNum += 1;
			ChaffNum = 1;
			WheatChaff wchf(this->livingFrac);
			organDataFrame->h_gChaf.push_back(wchf);
			organDataFrame->h_gChafLivingFrac.push_back(this->livingFrac);
			organDataFrame->h_gChafForce2Death.push_back(false);
		}

		//Z six rachis
		if (develop->is_doubleRidge() && RachisNum == 0)
		{
			SubRachisIdx = organDataFrame->h_gRchsNum;
			organDataFrame->h_gRchsNum += 1;
			RachisNum = 1;
			WheatRachis wrhs(this->livingFrac);
			organDataFrame->h_gRchs.push_back(wrhs);
			organDataFrame->h_gRchsLivingFrac.push_back(this->livingFrac);
			organDataFrame->h_gRchsForce2Death.push_back(false);
		}
	}
	//Z TODO growth stage, HAUN, FEEKERS, ZODAK
}

//Z DEATH of leaf and internode only for this tiller
//     not a destruction of the tiller, internode or leaf at this step
//     because those classes should exists to hold drop leaf and dead internode mass
void WheatTiller::WheatTillerSingleDeath(void)
{
	//Z prevent kill the same tiller twice
	if (living)
	{
		//Z tiller (and its leaf and internode) dies within the "tiller update" function
		//	therefore, although we set living = false and livingFrac = 0.0 here
		//	it will still call leaf and internode update to finish the death
		
		//Z even the tiller is dead,
		//	single tiller summary function will be always called to account living AND dead tissues
		
		living = false;
		livingFrac = 0.0f;
		force_to_death_current_step = true;
		
		//Z call leaf and internode to die
		for (int ii = 0; ii < LeafNum; ii++)
		{
			//Z leaf death, which will also set leaf's "force_to_death_current_step=true",
			//  s.t., later leaf area update and leaf senescence can be called normally.
			//  do not reset the "livingFrac" to zero for the leaf, in "force_to_death", it will be taken cared by the leaf module
			organDataFrame->h_gLeafForce2Death[SubLeafIdx[ii]] = true;
		}
		for (int ii = 0; ii < InterNodeNum; ii++)
		{
			//Z internode death, which will also set internode's "force_to_death_current_step=true",
			//  s.t., later internode update (and possibly future senescence) can be called normally.
			//  do not reset the "livingFrac" to zero for the intenrode, in "force_to_death", it will be taken cared by the internode module
			organDataFrame->h_gIntrForce2Death[SubInterNodeIdx[ii]] = true;
		}

		for (int ii = 0; ii < SpikNum; ii++)
		{
			//Z spike death, which will also set spike's "force_to_death_current_step=true",
			//  s.t., later spike update can be called normally.
			//  do not reset the "livingFrac" to zero for those spikes, in "force_to_death", it will be taken cared by the internode module
			organDataFrame->h_gSpikForce2Death[SubSpikIdx[ii]] = true;
		}

		//Z chaff death
		organDataFrame->h_gChafForce2Death[SubChaffIdx] = true;
		//Z rachis death
		organDataFrame->h_gRchsForce2Death[SubRachisIdx] = true;
	}
}

//Z GROWTH of leaf, internode, etc. only for this tiller
//     include leaf, internode, etc. but DO NOT sub-tiller, 
//     sub-tiller manages its own growth 
void WheatTiller::WheatTillerSingleUpdate(void)
{

	//Z this function call the update function for each leaf and internode in this tiller
	//  shape growth only
	//  no mass growth, mass growth will be based on the assignment later
	if (living)
	{
		//------------------ compute the morphology of the current tiller ---------------------
		//Z the highest ligulation, until the previous leaf
		//  the cumulated internode heights
		float cumuIntrLength_temp = 0.0f;
		for (int ii = 0; ii < InterNodeNum; ii++) {
			cumuIntrLength_temp += organDataFrame->h_gIntr[SubInterNodeIdx[ii]].InLength;
			cumuIntrLength[ii] = cumuIntrLength_temp;
		}
		cumuLigulationHeight[0] = 3.0f;
		if (InterNodeNum > 0) {
			cumuLigulationHeight[1] = organDataFrame->h_gIntr[SubInterNodeIdx[0]].InLength + 3.0f;
		}
		for (int ii = 2; ii < LeafNum; ii++) {
			cumuLigulationHeight[ii] = cumuIntrLength[ii-1] + 3.0f;
			for (int jj = ii - 1; jj > 1; jj--) {
				if (organDataFrame->h_gLeaf[SubLeafIdx[jj]].mature) {
					cumuLigulationHeight[ii] = cumuIntrLength[jj] + organDataFrame->h_gLeaf[SubLeafIdx[jj]].ShLength;
					break;
				}
			}
		}
		
		//------------------ Computing Env Stress and Tiller living Fraction ------------------
		// 
		//Z adjust tiller living fraction based on external stress
		//	external stress are based on nitrogen and water
		//	in SHOOTGRO, this process is done in "Evalst.for"
		singleRidge = develop->is_singleRidge();
		flowerPrimInit = develop->is_flowerPrimInit();
		anthesisStart = develop->is_startAnthesis();
		anthesisEnd = develop->is_endAnthesis();

		//Z this is just compariable but it is simplified, Z does not know how to correctly use it, and hence use the livingfrac functions
		//	1. first compute nitrogen, water stress, and temperature stresses for the tiller
		//	2. then determine if a new minimal stress value (new highest stress occurs)
		//	3. use the new stress value to reduce the livingfraction
		//	4. if livingfraction < 0.01 and it is not the mainstem, then call tiller death function

		//Z use plant level N content to count for the tiller growth stress
		//Z note that N is in mg and biomass is in g, and the content is based on %
		float pltbiomass = develop->get_pltLivingBiomass();
		float pltnitrogenmass = develop->get_pltLivingNitrogenMass();
		float pltNitrogenCond = 0.10f * pltnitrogenmass / pltbiomass;
		pltNitrogenCond = std::max(MIN_N_PCT, pltNitrogenCond);
		//Z small plant biomass implies young plant does not have N stress
		if (pltbiomass < 0.01f) { pltNitrogenCond = MAX_N_PCT; }

		//Z to give more Discrimination at high N end
		N_effect = std::max(0.0f, (2.0f / (1.0f + exp(-1.0f * (pltbiomass - MIN_N_PCT))) - 0.6f));
		N_effect = std::min(1.0f, N_effect);
		N_effect = std::max(0.1f, N_effect);

		water_effect = 1.0f;
		shade_effect = 1.0f;

		//Z shaded effects on tiller, marked based on R:FR ratio, only affect plant with more than 4 tillers.
		if (develop->get_pltTillerNum() >= 4.0f) { shade_effect = develop->get_shadeEffect(); }

		Cold_Time = develop->get_ColdTime();
		if (singleRidge) {
			TTd_since_singleRidge = develop->get_TTd_sinceSingleRidge();
			tmpr_effect = std::min(std::max(1.0f - (1.0f - tmpr_effect_terminal) / PHYLLOCHRON * (TTd_since_singleRidge - PHYLLOCHRON), tmpr_effect_terminal), 1.0f);
		}
		tmpr_effect = std::min(tmpr_effect, 1.0f);
		tmpr_effect = std::max(tmpr_effect, 0.01f);

		StressTiller = 1.0f;
		//StressTiller = tmpr_effect * water_effect * N_effect * shade_effect;

		if (StressTiller < StressTillerMin)
		{
			livingFrac = livingFrac * (StressTiller / StressTillerMin);
			StressTillerMin = StressTiller;
			if (mainstem) {
				livingFrac = std::max(livingFrac, 0.90f);
			}
			else {
				if (livingFrac < 0.01f) { this->WheatTillerSingleDeath(); }
			}
		}

		//------------------------ Programmed Tiller Abortion ------------------

		//Z call tiller death on jointing stage
		//  when jointing occurs, tiller with leaf number < 4 or born after single ridge will be killed

		if (jointing && (!jointing_already))
		{
			int jointingKill_LeafNum = 4;
			/*
			float gdddiff = develop->get_TTd_Plant() - develop->get_TTd_Joint();
			float cold_time_frac = develop->get_ColdTimeRatioJoint();
			if (Cold_Time <= 30.0f)
			{
				jointingKill_LeafNum = 0;
				tmpr_effect_terminal = 1.0f;
			}
			else
			{
				jointingKill_LeafNum = 1;
				tmpr_effect_terminal = 2.0f / (1.0f + exp(4.0f * std::max(cold_time_frac - 0.1f, 0.0f)));
			}
			*/
			jointing_already = jointing;
			TTd_since_jointing = 0.0f;
			if (EmergeLfNum < jointingKill_LeafNum && (!mainstem)) { jointing_kill = true; }
			if (born_after_singleRidge) { jointing_kill = true; }
		}

		//Z this statement means, if the tiller should be aborted, let it not happen immediately, but in a period of, say, 1 PHYLLOCHRON
		if (jointing_kill && livingFrac >= 0.01f)
		{
			TTd_since_jointing = develop->get_TTd_sinceJointing();
			livingFrac = livingFrac * (PHYLLOCHRON - TTd_since_jointing - Cur_TTd) / (PHYLLOCHRON - TTd_since_jointing);
			if (livingFrac < 0.01f) { this->WheatTillerSingleDeath(); }
		}

		//------------------------ Elongation and Flag Leaf Calculation ------------------

		//Z internode elongation or not
		//   only run this part once to determine the fist and the last internode to elongate for this tiller
		//   internode status is for the whole plant, but the internode that can grow is based on each tiller's situation

		//Z no need to worry about tiller after that
		//   by definition, after elongation, there will be no new tillers

		if (elongationStart && (!elongationStart_already))
		{
			if (Cold_Time >= 50.0f)
			{
				tmpr_effect_terminal = 0.1f;
			}

			//Z first get gdd parameters from development
			float mingdf = develop->get_TTd_FlagLf_min();
			float gddj = develop->get_TTd_Joint();
			float gddp = develop->get_TTd_Plant();
			float gdde = develop->get_TTd_Elong();

			//Z how many leaves can grow on this tiller in future
			//Z higher order tiller has less leaf and less leaf to grow
			int lnum = std::max(static_cast<int>((mingdf - TillerLfPhyllchron[EmergeLfNum - 1] - gddpBranch) / PHYLLOCHRON), 0);
			//Z flag leaf may not be the MAXLEAF, it depends on ambient and plant structure
			//  see "gddf" in the shootgro
			FlagLfPhyllchron_Cplt = TillerLfPhyllchron[EmergeLfNum - 1] + static_cast<float>(lnum + 1) * PHYLLOCHRON;

			for (int ii = EmergeLfNum + lnum + 1; ii < MAXLEAFNUM; ii++)
			{
				//Z no more leaf can grow beyond current leaf num + future leaf number "lnum", set Phyllchron for leaf emergence to arbitrary large
				TillerLfPhyllchron[ii] = 9999.0f;
			}

			//Z first (lowest) internode that can elongate
			//  the idea is if the current leaf is complete, then use "LeafNum - 1"
			//  if the current leaf is still expanding, the use "LeafNum - 2"
			//  finally, the first two internode should never elongate

			//Z however, the "rank" variable should always start from 0, 
			//  so addtional "-1" should be take out from LeafNum
			//  same for the max(), where "3" should be "2" since the first two internodes are 0 and 1

			//Z as an simple comment, the following code means the elongated internode is "2 nodes" prior to the current growing leaf
			if ((TillerLfPhyllchron[EmergeLfNum - 1] + gddpBranch + PHYLLOCHRON) <= (gddp - gddj + gdde))
			{
				//elongation_first = LeafNum - 1;
				elongation_first = EmergeLfNum - 2;
			}
			else
			{
				//elongation_first = LeafNum - 2;
				elongation_first = EmergeLfNum - 3;
			}
			//elongation_first = max(elongation_first, 3);
			elongation_first = std::max(elongation_first, 2);

			//Z last (highest) internode that can elongate
			elongation_last = elongation_first + lnum;
			elongation_last = std::min(elongation_last, MAXINTRNUM);

			//Z make sure this block only run once.
			elongationStart_already = elongationStart;
		}

		//Z under elongation status, control which internode should grow
		if (elongationStart)
		{
			//Z first collect plant gdd number need to use
			float gddp = develop->get_TTd_Plant();

			//case 1: If the last internode on the culm (= peduncle) has completed elongation, after initiation of the last leaf (= flag leaf).
			//		  It is 3*PCHRON because internodes grow for 1 phyllochron after a lag of 2 phyllochrons from when the associated leaf appears.
			if (elongation_last > 0)
			{
				//Z the current leaf is already the last one, so the even the last internode loses its chance 
				//  the "FlagLfPhyllchron_Cplt" is computed based on "this tiller", so no need to add "gddpBranch"
				if ((TillerLfPhyllchron[EmergeLfNum - 1] + PHYLLOCHRON) >= FlagLfPhyllchron_Cplt)
				{
					//Z here the comparison is based on "gddp", so need to add "gddpBranch"
					if (gddp >= (TillerLfPhyllchron[elongation_last - 1] + 3.0f * PHYLLOCHRON + gddpBranch))
					{
						//Z all the internodes that can elongate should already finish their elongation growth
						rank_2_elongation = -1;
					}
				}
			}

			//case 2: Grow the peduncle (or last internode to elongate on the culm) if appropriate.
			if (elongation_last > 0)
			{
				//Z the current leaf is already the last one, the last one can still grow
				//  the "FlagLfPhyllchron_Cplt" is computed based on "this tiller", so no need to add "gddpBranch"
				if ((TillerLfPhyllchron[EmergeLfNum - 1] + PHYLLOCHRON) >= FlagLfPhyllchron_Cplt)
				{
					//Z here the comparison is based on "gddp", so need to add "gddpBranch"
					if ((gddp >= (TillerLfPhyllchron[elongation_last - 1] + 2.0f * PHYLLOCHRON + gddpBranch)) && (EmergeLfNum > 1))
					{
						//Z at this time, the flagleaf is growing and the last internode should grow
						//  again, because the "rank" starts from 0, should make out "1" to match the order.
						//rank_2_elongation = LeafNum;
						rank_2_elongation = EmergeLfNum - 1;
					}
				}
			}

			//case 3: Grow the penultimate internode on the culm if appropriate.
			if ((elongation_last - 1 > 0) && (EmergeLfNum > 1))
			{
				//Z gdd is sufficient for plant growth and the corresponding leaf exists
				//Z here the comparison is based on "gddp", so need to add "gddpBranch"
				if ((gddp >= (TillerLfPhyllchron[elongation_last - 2] + 2.0f * PHYLLOCHRON + gddpBranch)) && SubLeafIdx[elongation_last - 2] != -1)
				{
					//Z at this time, the last-but-one leaf could grow and the last-but-two internode can grow 
					//  again, because the "rank" starts from 0, should make out "1" to match the order.
					//rank_2_elongation = LeafNum - 1;
					rank_2_elongation = EmergeLfNum - 2;
				}
			}

			//case 4: If only two leaves have appeared on the culm at the beginning of today, no internode will grow today
			//        That tiller is not ready to elongate
			if (EmergeLfNum <= 2)
			{
				//Z here the comparison is based on "gddp", so need to add "gddpBranch"
				if (gddp >= (TillerLfPhyllchron[EmergeLfNum - 1] + PHYLLOCHRON + gddpBranch))
				{
					rank_2_elongation = -1;
				}
			}

			//case 5: The currently growing internode is the one corresponding to the currently growing leaf - 2 and is not the peduncle or penultimate internode.
			if (EmergeLfNum >= 3)
			{
				//Z note that "rank" starts from 0
				rank_2_elongation = EmergeLfNum - 3;
			}
		}

		//------------------------ Flag Leaf Initiation And Flower Position ------------------
		//Z the initialiation of the flag leaf marks the initiation of the flower abortion
		if (physAge >= (FlagLfPhyllchron_Cplt - PHYLLOCHRON)) { FlagLfPhyllchron_Init = true; }
		//Z flower initiation
		if (flowerPrimInit && (!FlagLfPhyllchron_Init))
		{
			if (!flowerPrimInit_already)
			{
				flowerPrimInit_already = flowerPrimInit;
				spiint = SpikNum - SpikAdd;
				spifrt = spiint;
				spitop = spiint;
				spibot = spiint;
			}
			//Z every time, this will be the range of spikelet that will grow flowers
			spitop = min(spitop + SpikAdd, SpikNum);
			spibot = max(spibot - SpikAdd, 0);
		}
		//Z flower abortion
		//  between "FlagLfPhyllchron_Init" and "antess" all spikelet starts flower primabortion
		//Z flower fertilization
		if (anthesisStart && (!anthesisEnd))
		{
			phyAge_FlwFert += Cur_TTd;
			if (!anthesisStart_already)
			{
				anthesisStart_already = anthesisStart;
				ferint = spifrt;
				ferbot = ferint;
				fertop = ferint;
				this->develop->set_kernelFertStart();

			}
			if (phyAge_FlwFert > TTDBASELFLORET)
			{
				FlwFertAdd = 1;
				//Z reduce the phyAge for basel flower fertilization and wait for the next flower
				phyAge_FlwFert -= TTDBASELFLORET;
			}
			fertop = min(fertop + FlwFertAdd, spitop);
			ferbot = max(ferbot - FlwFertAdd, 0);
		}

		//Z -------------------------------------------------------------------------------------------
		//  - finish the phenology computation and start setting for the input tuple for organ update -
		//  -------------------------------------------------------------------------------------------
		//Z set up parameters for organ operation
				
		//Z set up leaf growth data
		for (int ii = 0; ii < LeafNum; ii++)
		{
			//Z leaf area update for existing leaves
			int IdxLeaf = SubLeafIdx[ii];
			if (IdxLeaf != -1) 
			{ 
				organDataFrame->h_gLeafLivingFrac[IdxLeaf] = this->livingFrac; 
				organDataFrame->h_gLeafForce2Death[IdxLeaf] = this->force_to_death_current_step;
				if (ii == 0) { 
					organDataFrame->h_gLeafPreLemerge[IdxLeaf] = false; 
				}
				else {
					organDataFrame->h_gLeafPreLemerge[IdxLeaf] = organDataFrame->h_gLeaf[SubLeafIdx[ii - 1]].emerge;
				}
				organDataFrame->h_gLeafLigulationDist[IdxLeaf] = max(cumuLigulationHeight[ii] - cumuIntrLength[ii], 0.0f);
			}	
		}

		//Z set up internode growth data
		for (int ii = 0; ii < InterNodeNum; ii++)
		{
			//Z internode update for existing internodes
			int IdxIntr = SubInterNodeIdx[ii];
			if (IdxIntr != -1)
			{
				//Z tiller controls which internode (or rank) should elongate
				organDataFrame->h_gIntrElongation[IdxIntr] = static_cast<bool>(ii == rank_2_elongation);
				organDataFrame->h_gIntrLivingFrac[IdxIntr] = this->livingFrac;
				organDataFrame->h_gIntrForce2Death[IdxIntr] = this->force_to_death_current_step;
				if (ii == 0) {
					organDataFrame->h_gIntrPreLigulation[IdxIntr] = organDataFrame->h_gLeaf[IdxIntr].mature;
				}
				else {
					organDataFrame->h_gIntrPreLigulation[IdxIntr] = organDataFrame->h_gLeaf[SubLeafIdx[ii - 1]].mature;
				}
				organDataFrame->h_gIntrLigulationDist[IdxIntr] = max(cumuLigulationHeight[ii] - cumuIntrLength[ii], 0.0f);
			}
		}

		//Z set up chaff growth data
		if (SubChaffIdx != -1)
		{
			organDataFrame->h_gChafLivingFrac[SubChaffIdx] = this->livingFrac;
			organDataFrame->h_gChafForce2Death[SubChaffIdx] = this->force_to_death_current_step;
		}

		//Z set up rachis growth data
		if (SubRachisIdx != -1)
		{
			organDataFrame->h_gRchsLivingFrac[SubRachisIdx] = this->livingFrac;
			organDataFrame->h_gRchsForce2Death[SubRachisIdx] = this->force_to_death_current_step;
		}

		//Z set up spike growing data
		bool aaaa = (flowerPrimInit && (!FlagLfPhyllchron_Init));	//for flower initiation
		bool bbbb = (FlagLfPhyllchron_Init && (!anthesisStart));	//for flower abortion
		bool cccc = (anthesisStart && (!anthesisEnd));				//for flower fertilization
		for (int ii = 0; ii < SpikNum; ii++)
		{
			//Z spike update for exisitng spikes
			int IdxSpik = SubSpikIdx[ii];
			if (IdxSpik != -1) 
			{
				organDataFrame->h_gSpikLivingFrac[IdxSpik] = this->livingFrac; 
				organDataFrame->h_gSpikForce2Death[IdxSpik] = this->force_to_death_current_step;
			}
			//Z for each spike, control its flower initiation, abortion and fertilization
			if (aaaa)		
			{
				organDataFrame->h_gSpikFlowerInit[IdxSpik] = static_cast<bool>((ii >= spibot) && (ii < spitop));
			}
			organDataFrame->h_gSpikFlowerAbort[IdxSpik] = bbbb;
			if (cccc)
			{
				organDataFrame->h_gSpikFlowerFert[IdxSpik] = static_cast<bool>((ii >= ferbot) && (ii < fertop));
			}
		}
	}
	else
	{
		//Z not living anymore
		//  still need to call once to finalize suborgan, leaf and internode, calculations for once
		
		if (death_2_finalize == false)
		{
			//Z set up leaf growth data
			for (int ii = 0; ii < LeafNum; ii++)
			{
				int IdxLeaf = SubLeafIdx[ii];
				if (IdxLeaf != -1) 
				{ 
					organDataFrame->h_gLeafLivingFrac[IdxLeaf] = 0.0f;
					organDataFrame->h_gLeafForce2Death[IdxLeaf] = this->force_to_death_current_step;
				}
			}

			//Z set up internode growth data
			for (int ii = 0; ii < InterNodeNum; ii++)
			{
				int IdxIntr = SubInterNodeIdx[ii];
				if (IdxIntr != -1)
				{
					organDataFrame->h_gIntrElongation[IdxIntr] = false;
					organDataFrame->h_gIntrLivingFrac[IdxIntr] = 0.0f;
					organDataFrame->h_gIntrForce2Death[IdxIntr] = this->force_to_death_current_step;
				}
			}

			//Z set up chaff growth data
			if (SubChaffIdx != -1)
			{
				organDataFrame->h_gChafLivingFrac[SubChaffIdx] = 0.0f;
				organDataFrame->h_gChafForce2Death[SubChaffIdx] = this->force_to_death_current_step;
			}

			//Z set up rachis growth data
			if (SubRachisIdx != -1)
			{
				organDataFrame->h_gRchsLivingFrac[SubRachisIdx] = 0.0f;
				organDataFrame->h_gRchsForce2Death[SubRachisIdx] = this->force_to_death_current_step;
			}

			//Z set up spike growth data
			for (int ii = 0; ii < SpikNum; ii++)
			{
				int IdxSpik = SubSpikIdx[ii];
				if (IdxSpik != -1) 
				{ 
					organDataFrame->h_gSpikLivingFrac[IdxSpik] = 0.0f;
					organDataFrame->h_gSpikForce2Death[IdxSpik] = this->force_to_death_current_step;
				}
			}

			death_2_finalize = true;
		}
	}
}

