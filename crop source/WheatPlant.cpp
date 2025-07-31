#include "WheatPlant.h"


WheatPlant::WheatPlant(const TInitInfo& info, TGasExSpeciesParam& photoparam)
{
	initInfo = info;
	gasExparam = photoparam;

	//Z initialize the plant body
	develop = make_unique<WheatDevelopment>(initInfo);
	//Z allocate an element for the public data frrame
	organDataFrame = make_unique<WheatOrganDataFrame>();
	//Z we start the plant by initialize the first tiller, aka the mainstem
	//  the initial living fraction of the mainstem MUST be 1.0
	//Z the mainstem does not have "rank" and its order is 0;
	mainstem = make_unique<WheatTiller>(0, 0, 0, true, develop.get(), organDataFrame.get(), 1.0f);
	//Z set pointer to NULL before allocate memory
	sunlit = NULL;
	shaded = NULL;

	//Z set some plant update stages regarding to flower
	flowerPrimInit = false;
	anthesisStart = false;
	anthesisEnd = false;

	//Z initialize plant part sizes
	SeedMass = 0.04f; // g seed^-1
	SeedNitrogenMass = SeedMass * MAX_N_PCT * 10.0f; // total weight of N per seed (mg N seed^-1, aka, mg N plant^-1)
	WheatMass = SeedMass;
	ShootMass = SeedMass * 0.6f;
	RootMass = SeedMass * 0.4f;
	PlantLivingFraction = 1.0f;

	//Z initialize leaf number, area, mass
	//  by default, this will be representative numbers, rather than "one single plant" number
	//  length in cm, area in cm^2, biomass in g, nitrogen mass in mg
	LfNumPlt = 0.0f;
	greenLfNumPlt = 0.0f;
	dropLfNumPlt = 0.0f;

	LaAreaPlt = 0.0f;
	greenLaAreaPlt = 0.0f;
	seneLaAreaPlt = 0.0f;
	dropLaAreaPlt = 0.0f;

	LaMassPlt = ShootMass * 0.9f;	// 90% for lamina, while 10% for internode
	ShMassPlt = 0.0f;
	dropLaMassPlt = 0.0f;
	dropShMassPlt = 0.0f;
	LaNitrogenMassPlt = 0.0f;
	ShNitrogenMassPlt = 0.0f;
	LaNitrogenReturnPlt = 0.0f;
	ShNitrogenReturnPlt = 0.0f;

	deadLaNitrogenMassPlt = 0.0f;
	deadShNitrogenMassPlt = 0.0f;
	seneLaNitrogenMassPlt = 0.0f;
	seneShNitrogenMassPlt = 0.0f;
	dropLaNitrogenMassPlt = 0.0f;
	dropShNitrogenMassPlt = 0.0f;

	ptnLaMassIncreasePlt = 0.0f;
	ptnShMassIncreasePlt = 0.0f;
	ptnLfMassIncreasePlt = 0.0f;
	ptnLaNitrogenMassIncreasePlt = 0.0f;
	ptnShNitrogenMassIncreasePlt = 0.0f;
	ptnLfNitrogenMassIncreasePlt = 0.0f;

	//Z initialize tiller and internode number, mass, length
	//  length in cm, area in cm^2, biomass in g, nitrogen mass in mg
	TlNumSingle = 0;
	TlNumPlt = 0.0f;
	InLengthPlt = 0.0f;
	InMassPlt = ShootMass * 0.1f;	// 10% for internode, while 90% for lamina
	InNitrogenMassPlt = 0.0f;
	InNitrogenReturnPlt = 0.0f;
	deadInNitrogenMassPlt = 0.0f;

	ptnInMassIncreasePlt = 0.0f;
	ptnInNitrogenMassIncreasePlt = 0.0f;

	//***** Plant Chaff Group (cm or g) *********
	CfMassPlt = 0.0f;
	deadCfMassPlt = 0.0f;
	CfNitrogenMassPlt = 0.0f;
	CfNitrogenReturnPlt = 0.0f;
	deadCfNitrogenMassPlt = 0.0f;

	ptnCfMassIncreasePlt = 0.0f;
	ptnCfNitrogenMassIncreasePlt = 0.0f;

	//***** Plant Rachis Group (cm or g) *********
	RaLengthPlt = 0.0f;
	RaMassPlt = 0.0f;
	deadRaMassPlt = 0.0f;
	RaNitrogenMassPlt = 0.0f;
	RaNitrogenReturnPlt = 0.0f;
	deadRaNitrogenMassPlt = 0.0f;

	ptnRaMassIncreasePlt = 0.0f;
	ptnRaNitrogenMassIncreasePlt = 0.0f;

	//***** Plant Spikelet ( Flower + Kernel ) Group (cm or g) *********
	KrNumPlt = 0.0f;
	KrMassSumPlt = 0.0f;
	KrNitrogenMassSumPlt = 0.0f;

	ptnKrMassIncreaseSumPlt = 0.0f;
	ptnKrNitrogenMassIncreaseSumPlt = 0.0f;

	//Z initialize root mass, actual root growth and hence the N demand
	//  root biomass allocation has its own rule, 
	//  while root N is affliated to root biomass
	//  therefore, we only assume potential root N mass
	RtNitrogenMassPlt = 0.0f;
	ptnRtNitrogenGrowth = 0.0f;

	//Z initialize plant based mass and N variables, % number
	NitrogenMassPlt = SeedNitrogenMass; // total weight of N per seed (mg N seed^-1, aka, mg N plant^-1)
	NitrogenRatio = 1.0f;
	LaNitrogenContentPlt = MAX_N_PCT;

	//Z photosynthesis and gas exchange model
	sunlit_LAI = 0.0f;		shaded_LAI = 0.0f;
	sunlit_PFD = 0.0f;		shaded_PFD = 0.0f;
	sunlit_A_gross = 0.0f;	shaded_A_gross = 0.0f;
	sunlit_A_net = 0.0f;	shaded_A_net = 0.0f;
	assimilate = 0.0f;
	photosynthesis_gross = 0.0f;
	photosynthesis_net = 0.0f;
	transpiration = 0.0f;
	transpirationOld = 0.0f;
	temperature = develop->get_TmprCur();
	VPD = 0.0f;
	sunlit_gs = 0.0f; shaded_gs = 0.0f;
	conductance = 0.0f;

	//Z initialize the biomass partitioning factors
	//  pools and storage factors
	BiomassReserve = 0.00001f;
	BiomassPool = SeedMass;
	BiomassSupply = 0.0f;
	BiomassDemand = 0.0f;
	BiomassLeftover = 0.0f;
	BiomassRoot = 0.00001f;

	BiomassLeafGrowth_ptn = 0.0f;
	BiomassIntrNodeGrowth_ptn = 0.0f;
	BiomassChaffGrowth_ptn = 0.0f;
	BiomassRachisGrowth_ptn = 0.0f;
	BiomassKernelGrowth_ptn = 0.0f;
	BiomassShootGrowth_ptn = 0.0f;
	BiomassPltGrowth_ptn = 0.0f;

	ActuralRootBiomassAssignment_PCRL = 0.0f;
	ActuralShootBiomassAssignment = 0.0f;

	NitrogenPool = SeedNitrogenMass;

	//  mass partitioning to plant organs
	shootPart = 0.0f;
	rootPart = 0.0f;
	shootPart_old = 0.0f;
	rootPart_old = 0.0f;
	leafPart = 0.0f;
	internodePart = 0.0f;
	chaffPart = 0.0f;
	rachisPart = 0.0f;
	spikeletPart = 0.0f;

	leafPartNitrogen = 0.0f;
	internodePartNitrogen = 0.0f;
	chaffPartNitrogen = 0.0f;
	rachisPartNitrogen = 0.0f;
	spikeletPartNitrogen = 0.0f;

	LfBiomassAssignmentRate = 0.0f;
	InBiomassAssignmentRate = 0.0f;
	CfBiomassAssignmentRate = 0.0f;
	RaBiomassAssignmentRate = 0.0f;
	SpBiomassAssignmentRate = 0.0f;

	maintRespiration = 0.0f;

	//Z all the following nitrogen is based on mg plant^-1
	HourlyNitrogenDemand = 0.0f;
	CumulativeNitrogenDemand = 0.0f;
	HourlyNitrogenSoilUptake = 0.0f;
	CumulativeNitrogenSoilUptake = 0.0f;

	NitrogenLeafGrowth_ptn = 0.0f;
	NitrogenIntrNodeGrowth_ptn = 0.0f;
	NitrogenChaffGrowth_ptn = 0.0f;
	NitrogenRachisGrowth_ptn = 0.0f;
	NitrogenKernelGrowth_ptn = 0.0f;

	NitrogenDemand = 0.0f;
	NitrogenShootGrowth_ptn = 0.0f;
	ShootNitrogenAvailiableAllocation = 0.0f;
	RootNitrogenAvailiableAllocation = 0.0f;

	LfNitrogenAssignmentRate = 0.0f;
	InNitrogenAssignmentRate = 0.0f;
	CfNitrogenAssignmentRate = 0.0f;
	RaNitrogenAssignmentRate = 0.0f;
	SpNitrogenAssignmentRate = 0.0f;

	//***** Plant Information *********
	sowingDay = 1.0f;
	age = 0.0f;
	emerge_gdd = 0.0f;
	SunlitRatio = 0.0f;
	C2_effect = 0.0f;

}

//Z even plant died, this does not mean plant objective will be destructed.
//  destruction only occurs when simulation is finished
//  smart pointer will be destoryed automatically
WheatPlant::~WheatPlant()
{
}

void WheatPlant::WheatPlantUpdate(const TWeather& weather)
{
	//Z update the development stage
	develop->WheatDelpUpdate(weather);
	//Z this living fraction is NOT the same as the tiller, leaf, internode livingfrac
	//  this is for emergence percentage
	PlantLivingFraction = develop->get_plantLivingFraction();

	//Z update the plant morphology
	//Z first germination starts (not complete), just mark the time
	//  germination start, namely "germination initiation, germinInit", should be always earlier than anything
	if (!develop->is_germinInit())
	{
		//Z seed is not germinated yet, nothing need to be done
		temperature = develop->get_TmprCur();
	}

	//Z germination marks 1/2 of the "active" seeds germinate, based on a germination ratio
	//  may complete after emergence
	//  so germination seems not that useful for plant.cpp, because in the module, we still tend to simulate one representative plant
	//  but maybe important for output
	if ((develop->is_germinInit()) && (!develop->is_germination()))
	{
		temperature = develop->get_TmprCur();
	}

	//Z between "germination start" and "emergence start"
	//  representative plant, start to trace
	//  seed is the source of biomass and nitrogen
	if ((develop->is_germinInit()) && (!develop->is_emerge()))
	{
		temperature = develop->get_TmprCur();
		//Z mainstem initialization, 
		//  apsim assumes two leaves grow on the mainstem
		//  therefore, initialize the two leaves before emergence
		if (!mainstem->is_mainstemInitilized())
		{
			//Z mainstem initialization is just a two step thing. 
			//  so we just push data to both host and device vectors, with small communication cost
			mainstem->MainstemInitialize();
		}

		//Z update morph features, like leaf area, stem length, mass and N etc.
		//  this need to be done after we initialize 2 leaves and 2 internodes

		this->calsSetMorphology();

		//Z plant mass (above ground part): including biomass and nitrogen mass
		ShootMass = LaMassPlt + ShMassPlt + InMassPlt + CfMassPlt + RaMassPlt + KrMassSumPlt;
		WheatMass = SeedMass + ShootMass + RootMass;
		NitrogenMassPlt = SeedNitrogenMass + LaNitrogenMassPlt + ShNitrogenMassPlt + InNitrogenMassPlt + CfNitrogenMassPlt + RaNitrogenMassPlt + KrNitrogenMassSumPlt
			+ deadLaNitrogenMassPlt + deadShNitrogenMassPlt + deadInNitrogenMassPlt + deadCfNitrogenMassPlt + deadRaNitrogenMassPlt
			+ RtNitrogenMassPlt;

		//Z plant based leaf N content (aged leaves release N)
		//  in photosynthesis gas_exchange model, unit mg N per m^2 leaf area
		if (greenLaAreaPlt <= 0.1f) {
			// default N mass (g) per leaf area m^2, with N mass fraction 3.5% and specific leaf weight 0.000045 g/mm^2
			LaNitrogenContentPlt = 1.5750f / 3.5f * MAX_N_PCT;
		}
		else {
			// convert leafN mass from mg to g
			// convert green leaf area to m^2
			LaNitrogenContentPlt = (LaNitrogenMassPlt / 1000.0f) / (greenLaAreaPlt / 10000.0f);
		}

		//Z sync host and device position and start mass distribution
		//  then sum all the mass components
		WheatUpdateSummary();

		//Z potential growth can be consider as the demand
		//  biomass g; nitrogen mass mg
		BiomassLeafGrowth_ptn = ptnLfMassIncreasePlt;
		BiomassIntrNodeGrowth_ptn = ptnInMassIncreasePlt;
		BiomassChaffGrowth_ptn = ptnCfMassIncreasePlt;
		BiomassRachisGrowth_ptn = ptnRaMassIncreasePlt;
		BiomassKernelGrowth_ptn = ptnKrMassIncreaseSumPlt;

		BiomassDemand = BiomassKernelGrowth_ptn;
		BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn + BiomassChaffGrowth_ptn + BiomassRachisGrowth_ptn;

		NitrogenLeafGrowth_ptn = ptnLfNitrogenMassIncreasePlt;
		NitrogenIntrNodeGrowth_ptn = ptnInNitrogenMassIncreasePlt;
		NitrogenChaffGrowth_ptn = ptnCfNitrogenMassIncreasePlt;
		NitrogenRachisGrowth_ptn = ptnRaNitrogenMassIncreasePlt;
		NitrogenKernelGrowth_ptn = ptnKrNitrogenMassIncreaseSumPlt;

		NitrogenDemand = NitrogenKernelGrowth_ptn;
		NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn + NitrogenChaffGrowth_ptn + NitrogenRachisGrowth_ptn;

		RootMass += weather.pcrs;
		ptnRtNitrogenGrowth = std::max((RootMass * MAX_N_PCT_ROOT * 10.0f - RtNitrogenMassPlt), 0.0f);

		//Z maintanence respiration rate
		calcMaintRespiration(weather);

		//Z biomass allocation
		//  in this case, biomass pools are all from seeds, so need to decrease SeedMass step by step
		//  then set mass to the plant organs
		calsBiomassAllocation(weather);

		//Z nitrogen allocation
		//  in this case, nitrogen are all from seeds, with seed N fraction (3.4%)
		//  then set nitrogen mass to the plant organs
		//Z Nitrogen release from seed should be propotional to the BiomassSupply based on (3.4%)
		calsNitrogenAllocation();

		SeedMass = std::max(0.0f, SeedMass - BiomassSupply);
		SeedNitrogenMass = std::max(0.0f, SeedNitrogenMass - RootNitrogenAvailiableAllocation - ShootNitrogenAvailiableAllocation);

		calsSetMass();

		//Z sum the mass component at this time
		//  including biomass and nitrogen mass
		//  note seed mass is changing, but the nitrogen fraction in the seed mass is a constant
		ShootMass = LaMassPlt + ShMassPlt + InMassPlt + CfMassPlt + RaMassPlt + KrMassSumPlt;
		WheatMass = SeedMass + ShootMass + RootMass;
		NitrogenMassPlt = SeedNitrogenMass + LaNitrogenMassPlt + ShNitrogenMassPlt + InNitrogenMassPlt + CfNitrogenMassPlt + RaNitrogenMassPlt + KrNitrogenMassSumPlt
			+ deadLaNitrogenMassPlt + deadShNitrogenMassPlt + deadInNitrogenMassPlt + deadCfNitrogenMassPlt + deadRaNitrogenMassPlt
			+ RtNitrogenMassPlt;
	}
	//Z after emergence
	//Z need to obtain hourly root nitrogen uptake from 2DSOIL as nitrogen input
	if ((develop->is_emerge()) && (!develop->is_dead()))
	{
		//Z from this point, plant will grow in a "normal way" or say, on its "path"
		//  update the temperature and then compute the R:FR ratio, which will furtherly used as a limiter for tiller number
		temperature = develop->get_TmprCur();
		calcRed_FRedRatio(weather);

		//Z update morph features, like leaf area, stem length, mass and N etc.

		this->calsSetMorphology();

		//Z plant mass (above ground part)
		//  including biomass and nitrogen mass
		ShootMass = LaMassPlt + ShMassPlt + InMassPlt + CfMassPlt + RaMassPlt + KrMassSumPlt;
		WheatMass = SeedMass + ShootMass + RootMass;
		NitrogenMassPlt = SeedNitrogenMass + LaNitrogenMassPlt + ShNitrogenMassPlt + InNitrogenMassPlt + CfNitrogenMassPlt + RaNitrogenMassPlt + KrNitrogenMassSumPlt
			+ deadLaNitrogenMassPlt + deadShNitrogenMassPlt + deadInNitrogenMassPlt + deadCfNitrogenMassPlt + deadRaNitrogenMassPlt
			+ RtNitrogenMassPlt;

		//Z plant based leaf N content (aged leaves release N)
		if (greenLaAreaPlt <= 1.0f)
		{
			// default N mass (g) per leaf area m^2, with N mass fraction 3.5% and specific leaf weight 0.000045 g/mm^2
			LaNitrogenContentPlt = 1.5750f / 3.5f * MAX_N_PCT;
		}
		else
		{
			// convert leafN mass from mg to g
			// convert green leaf area cm^2 to m^2
			LaNitrogenContentPlt = (LaNitrogenMassPlt / 1000.0f) / (greenLaAreaPlt / 10000.0f);
		}

		//Z call the photosynthesis method
		//  the criterion for maize is "GreenLfAreaPlt>10 cm^2", but that could be too large for rye
		if (greenLaAreaPlt > 1.0f)
		{
			calcGasExchange(weather, gasExparam);
		}

		//Z sync host and device position and start mass distribution
		//  then sum all the mass components
		WheatUpdateSummary();

		//Z potential growth can be consider as the demand
		//  biomass g; nitrogen mass mg
		BiomassLeafGrowth_ptn = ptnLfMassIncreasePlt;
		BiomassIntrNodeGrowth_ptn = ptnInMassIncreasePlt;
		BiomassChaffGrowth_ptn = ptnCfMassIncreasePlt;
		BiomassRachisGrowth_ptn = ptnRaMassIncreasePlt;
		BiomassKernelGrowth_ptn = ptnKrMassIncreaseSumPlt;

		BiomassDemand = BiomassKernelGrowth_ptn;
		BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn + BiomassChaffGrowth_ptn + BiomassRachisGrowth_ptn;

		NitrogenLeafGrowth_ptn = ptnLfNitrogenMassIncreasePlt;
		NitrogenIntrNodeGrowth_ptn = ptnInNitrogenMassIncreasePlt;
		NitrogenChaffGrowth_ptn = ptnCfNitrogenMassIncreasePlt;
		NitrogenRachisGrowth_ptn = ptnRaNitrogenMassIncreasePlt;
		NitrogenKernelGrowth_ptn = ptnKrNitrogenMassIncreaseSumPlt;

		NitrogenDemand = NitrogenKernelGrowth_ptn;
		NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn + NitrogenChaffGrowth_ptn + NitrogenRachisGrowth_ptn;

		RootMass += weather.pcrs;
		ptnRtNitrogenGrowth = std::max((RootMass * MAX_N_PCT_ROOT * 10.0f - RtNitrogenMassPlt), 0.0f);

		//Z call maintainance respiration
		calcMaintRespiration(weather);

		//Z biomass partitioning between roots and shoots
		//  then partitioning shoots part to leaf (with sheath) and internode
		calsBiomassAllocation(weather);

		//Z nitrogen allocation
		//  in this case, nitrogen are from plant nitrogen release and root N uptake
		//  then set nitrogen mass to the plant organs
		//Z at this time, biomass and nitrogen may not be allocated at a fixed ratio, 
		//  and hence, plant organs will have 
		calsNitrogenAllocation();

		//Z assign photosynthesis to the biomass pools
		//  this comes after the biomass allocation, that means, the photosynthesis biomass will be used in the next time step
		if (abs(weather.time) < 0.0001f)
		{
			//Z by the end of a calender day, zero short term biomass pool
			//  add short term biomass pool to the long term biomass reservior
			BiomassReserve += std::max(0.0f, BiomassPool);
			BiomassPool = 0.0f; //reset shorterm C_pool to zero at midnight, needs to be more mechanistic
		}
		else
		{
			//Z if not at the end of a calender day, take photosynthesis biomass output
			//  should be assigned to the fast/short term biomass pool
			BiomassPool += std::max(photosynthesis_gross, 0.0f);
		}

		calsSetMass();

		//Z sum the mass component at this time
		//  including biomass and nitrogen mass
		//  note seed mass is changing, but the nitrogen fraction in the seed mass is a constant
		ShootMass = LaMassPlt + ShMassPlt + InMassPlt + CfMassPlt + RaMassPlt + KrMassSumPlt;
		WheatMass = SeedMass + ShootMass + RootMass;
		NitrogenMassPlt = SeedNitrogenMass + LaNitrogenMassPlt + ShNitrogenMassPlt + InNitrogenMassPlt + CfNitrogenMassPlt + RaNitrogenMassPlt + KrNitrogenMassSumPlt
			+ deadLaNitrogenMassPlt + deadShNitrogenMassPlt + deadInNitrogenMassPlt + deadCfNitrogenMassPlt + deadRaNitrogenMassPlt
			+ RtNitrogenMassPlt;

		//Z judge if plant is nearly dead
		//  currently, we use if the elongation start as the temporary criterion
		if ((greenLaAreaPlt <= 0.00005f * LaAreaPlt) && develop->is_startEnlongation())
		{
			develop->set_maturity(true, weather.daytime);
			develop->set_death(true, weather.daytime);
		}
	}
}

/* BLOCK
   ******** PLANT ORGANS MORPHOLOGY ********
*/
void WheatPlant::calsSetMorphology()
{
	//Z first reset tiller number for value update
	TlNumSingle = 0;
	TlNumPlt = 0.0f;

	//Z set some plant update stages regarding to flower
	flowerPrimInit = develop->is_flowerPrimInit();
	anthesisStart = develop->is_startAnthesis();
	anthesisEnd = develop->is_endAnthesis();

	//Z call the two-step functions for tiller number update
	//  step 1: new tiller and new organs
	//  step 2: update livingfractor, growing stage
	//  recursively loop all the tillers
	this->WheatPlantMorphoUpdate(mainstem.get());
	this->WheatPlantSetupUpdate();
	this->WheatPlantUpdateRun();
}

//Z Wheat plant morphological update
//  set up the tiller resursive computation, such that the preparation for the GPU computation can be ready
void WheatPlant::WheatPlantMorphoUpdate(WheatTiller* wt)
{
	//Z morph updates
	//  new leaves, new internodes, ...
	//  the tiller will say it grow new organs (add record), the plant will management the data structure, i.e., add data to the GPU data structure
	wt->WheatTillerSingleMorph();
	
	//Z set up the growth of existing organs,
	//  i.e., living factors, which internode should elongate ...
	wt->WheatTillerSingleUpdate();
	
	//Z add the current tiller number/living fraction to the cumulative values
	TlNumSingle += wt->is_tillerLiving();
	TlNumPlt += wt->get_tillerLivingFrac();

	//Z recursive process, starting with the mainstem
	int StNum = wt->get_StNum();
	for (int ii = 0; ii < StNum; ii++)
	{
		//Z first name this subtiller "subrt"
		WheatTiller* subrt = wt->get_subtiller(ii);
		if (subrt != nullptr) { WheatPlantMorphoUpdate(subrt); }
	}
}
//Z update the growth operators first,
//  this prepare the operator on GPU
void WheatPlant::WheatPlantSetupUpdate()
{
	//Z setup parameters
	float curTdd = develop->get_TTdCur_dssat();
	float curTmpr = develop->get_TmprCur();
	float curTeq = develop->get_TimeTeqCur_fspm();
	float psiPredawn = develop->get_PredawnLWP();
	float curShadeEffect = develop->get_shadeEffect();
	bool curAccel = develop->is_accel();
	bool curElong = develop->is_startEnlongation();
	bool curAntepa = develop->is_endAnthesis();

	//Z update leaf growth operator
	organDataFrame->wLeafUpdate.SetLeafUpdateCondition(curTdd, curTmpr, curTeq, psiPredawn, curShadeEffect, curAccel);
	//Z update internode growth operator
	organDataFrame->wIntrUpdate.SetInternodeUpdateCondition(curTdd, curTmpr, curTeq, psiPredawn, curShadeEffect, curElong, curAntepa);
	//Z update chaff growth opeartor
	organDataFrame->wChafUpdate.SetChaffUpdateCondition(curTdd, psiPredawn,
		develop->get_DurationChaffGrowth(), develop->is_chaffEnd());
	//Z update rachis growth operator
	organDataFrame->wRchsUpdate.SetRachisUpdateCondition(curTdd, psiPredawn,
		develop->get_DurationRachisGrowth(), develop->is_startJoint());
	//Z update spikelet (flower+kernel) operator
	organDataFrame->wSpikUpdate.SetSpikeletUpdateCondition(curTdd, curTmpr, psiPredawn,
		develop->get_pltLivingBiomass(), develop->get_pltLivingNitrogenMass(),
		develop->get_DayTime(), develop->get_DayTime_kernelFertStart(), develop->get_DayTime_afterAntss(),
		develop->get_kwt00(), develop->get_kwt01(), develop->get_kwt10(), develop->get_kwt11(),
		develop->get_kdur00(), develop->get_kdur01(), develop->get_kdur10(), develop->get_kdur11());
}
//Z invoke the GPU computation
//  first call "WheatPlantMorphoUpdate"+"WheatPlantSetupUpdate"

void WheatPlant::WheatPlantUpdateRun()
{
	//Z invoke the GPU async for_each method
	//  step 1: construct tuple
	//  step 2: run the code
	/*
	auto leafOperator = zip(organDataFrame->h_gLeaf, organDataFrame->h_gLeafLivingFrac, organDataFrame->h_gLeafForce2Death);
	auto intrOperator = zip(organDataFrame->h_gIntr, organDataFrame->h_gIntrLivingFrac, organDataFrame->h_gIntrForce2Death, organDataFrame->h_gIntrElongation);
	auto chafOperator = zip(organDataFrame->h_gChaf, organDataFrame->h_gChafLivingFrac, organDataFrame->h_gChafForce2Death);
	auto rchsOperator = zip(organDataFrame->h_gRchs, organDataFrame->h_gRchsLivingFrac, organDataFrame->h_gRchsForce2Death);
	auto spikOperator = zip(organDataFrame->h_gSpik, organDataFrame->h_gSpikLivingFrac, organDataFrame->h_gSpikForce2Death, organDataFrame->h_gSpikFlowerInit, organDataFrame->h_gSpikFlowerAbort, organDataFrame->h_gSpikFlowerFert);

	std::for_each(leafOperator.begin(), leafOperator.end(), organDataFrame->wLeafUpdate);
	std::for_each(intrOperator.begin(), intrOperator.end(), organDataFrame->wIntrUpdate);
	std::for_each(chafOperator.begin(), chafOperator.end(), organDataFrame->wChafUpdate);
	std::for_each(rchsOperator.begin(), rchsOperator.end(), organDataFrame->wRchsUpdate);
	std::for_each(spikOperator.begin(), spikOperator.end(), organDataFrame->wSpikUpdate);
	*/
	
	for (int ii = 0; ii < organDataFrame->h_gLeafNum; ii++) { organDataFrame->wLeafUpdate.LeafUpdateRun(organDataFrame->h_gLeaf[ii], organDataFrame->h_gLeafLivingFrac[ii], organDataFrame->h_gLeafForce2Death[ii], organDataFrame->h_gLeafPreLemerge[ii], organDataFrame->h_gLeafLigulationDist[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gIntrNum; ii++) { organDataFrame->wIntrUpdate.IntrUpdateRun(organDataFrame->h_gIntr[ii], organDataFrame->h_gIntrLivingFrac[ii], organDataFrame->h_gIntrForce2Death[ii], organDataFrame->h_gIntrPreLigulation[ii], organDataFrame->h_gIntrLigulationDist[ii], organDataFrame->h_gIntrElongation[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gChafNum; ii++) { organDataFrame->wChafUpdate.ChafUpdateRun(organDataFrame->h_gChaf[ii], organDataFrame->h_gChafLivingFrac[ii], organDataFrame->h_gChafForce2Death[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gRchsNum; ii++) { organDataFrame->wRchsUpdate.RchsUpdateRun(organDataFrame->h_gRchs[ii], organDataFrame->h_gRchsLivingFrac[ii], organDataFrame->h_gRchsForce2Death[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gSpikNum; ii++) { organDataFrame->wSpikUpdate.SpikUpdateRun(organDataFrame->h_gSpik[ii], organDataFrame->h_gSpikLivingFrac[ii], organDataFrame->h_gSpikForce2Death[ii], organDataFrame->h_gSpikFlowerInit[ii], organDataFrame->h_gSpikFlowerAbort[ii], organDataFrame->h_gSpikFlowerFert[ii]); }

}

/* END BLOCK
   ******** PLANT ORGANS MORPHOLOGY ********
*/

/* BLOCK
   ******** PLANT ORGANS MORPHOLOGY SUMMARY AND MASS DEMAND ********
*/
//Z during GPU update the plant morphology, the CPU portion can perform photosynthesis computation at the same time
//  then, CPU and GPU should be sync, such that the updated plant morphology can be summed (based on representative sense) and transfer back to CPU
//  then, CPU can get potential biomass and N demand and start the mass distribution
//Z-Caution!!! always sync the CPU and GPU before this function
void WheatPlant::WheatUpdateSummary()
{
	//Z need to reset the plant parameters before sum
	//  recall that tiller is a box that hold leaves and internodes
//	std::chrono::time_point<std::chrono::system_clock> t1, t2, t3, t4;
//	std::chrono::duration<double> elapsed_1, elapsed_2, elapsed_3;

	LfNumPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.livingFrac_ini; });
	greenLfNumPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.livingFrac; });
	dropLfNumPlt = LfNumPlt - greenLfNumPlt;
	LaAreaPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.LaArea_Rep; });
	greenLaAreaPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.greenLaArea_Rep; });
	seneLaAreaPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.seneLaArea_Rep; });
	dropLaAreaPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.dropLaArea_Rep; });
	LaMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.LaMass_Rep; });
	ShMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ShMass_Rep; });
	dropLaMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.dropLaMass_Rep; });
	dropShMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.dropShMass_Rep; });
	LaNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.LaNitrogenMass_Rep; });
	ShNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ShNitrogenMass_Rep; });
	LaNitrogenReturnPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.LaNitrogenReturn_Rep; });
	ShNitrogenReturnPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ShNitrogenReturn_Rep; });
	deadLaNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.deadLaNitrogenMass_Rep; });
	deadShNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.deadShNitrogenMass_Rep; });
	seneLaNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.seneLaNitrogenMass_Rep; });
	seneShNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.seneShNitrogenMass_Rep; });
	dropLaNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.dropLaNitrogenMass_Rep; });
	dropShNitrogenMassPlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.dropShNitrogenMass_Rep; });
	ptnLaMassIncreasePlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ptnLaMassIncrease_Rep; });
	ptnShMassIncreasePlt= std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ptnShMassIncrease_Rep; });
	ptnLfMassIncreasePlt = ptnLaMassIncreasePlt + ptnShMassIncreasePlt;
	ptnLaNitrogenMassIncreasePlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ptnLaNitrogenMassIncrease_Rep; });
	ptnShNitrogenMassIncreasePlt = std::accumulate(organDataFrame->h_gLeaf.begin(), organDataFrame->h_gLeaf.end(), 0.0f, [](float sum, const WheatLeaf& curr) {return sum + curr.ptnShNitrogenMassIncrease_Rep; });
	ptnLfNitrogenMassIncreasePlt = ptnLaNitrogenMassIncreasePlt + ptnShNitrogenMassIncreasePlt;

	InLengthPlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.InLength_Rep; });
	InMassPlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.InMass_Rep; });
	InNitrogenMassPlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.InNitrogenMass_Rep; });
	InNitrogenReturnPlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.InNitrogenReturn_Rep; });
	deadInNitrogenMassPlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.deadInNitrogenMass_Rep; });
	ptnInMassIncreasePlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.ptnInMassIncrease_Rep; });
	ptnInNitrogenMassIncreasePlt = std::accumulate(organDataFrame->h_gIntr.begin(), organDataFrame->h_gIntr.end(), 0.0f, [](float sum, const WheatInternode& curr) {return sum + curr.ptnInNitrogenMassIncrease_Rep; });

	CfMassPlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.CfMass_Rep; });
	deadCfMassPlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.deadCfMass_Rep; });
	CfNitrogenMassPlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.CfNitrogenMass_Rep; });
	CfNitrogenReturnPlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.CfNitrogenReturn_Rep; });
	deadCfNitrogenMassPlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.deadCfNitrogenMass_Rep; });
	ptnCfMassIncreasePlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.ptnCfMassIncrease_Rep; });
	ptnCfNitrogenMassIncreasePlt = std::accumulate(organDataFrame->h_gChaf.begin(), organDataFrame->h_gChaf.end(), 0.0f, [](float sum, const WheatChaff& curr) {return sum + curr.ptnCfNitrogenMassIncrease_Rep; });

	KrNumPlt = std::accumulate(organDataFrame->h_gSpik.begin(), organDataFrame->h_gSpik.end(), 0.0f, [](float sum, const WheatSpikelet& curr) {return sum + curr.KrNum_Rep; });
	KrMassSumPlt = std::accumulate(organDataFrame->h_gSpik.begin(), organDataFrame->h_gSpik.end(), 0.0f, [](float sum, const WheatSpikelet& curr) {return sum + curr.KrMassSum_Rep; });
	KrNitrogenMassSumPlt = std::accumulate(organDataFrame->h_gSpik.begin(), organDataFrame->h_gSpik.end(), 0.0f, [](float sum, const WheatSpikelet& curr) {return sum + curr.KrNitrogenMassSum_Rep; });
	ptnKrMassIncreaseSumPlt = std::accumulate(organDataFrame->h_gSpik.begin(), organDataFrame->h_gSpik.end(), 0.0f, [](float sum, const WheatSpikelet& curr) {return sum + curr.ptnKrMassIncreaseSum_Rep; });
	ptnKrNitrogenMassIncreaseSumPlt = std::accumulate(organDataFrame->h_gSpik.begin(), organDataFrame->h_gSpik.end(), 0.0f, [](float sum, const WheatSpikelet& curr) {return sum + curr.ptnKrNitrogenMassIncreaseSum_Rep; });

	RaLengthPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.RaLength_Rep; });
	RaMassPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.RaMass_Rep; });
	deadRaMassPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.deadRaMass_Rep; });
	RaNitrogenMassPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.RaNitrogenMass_Rep; });
	RaNitrogenReturnPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.RaNitrogenReturn_Rep; });
	deadRaNitrogenMassPlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.deadRaNitrogenMass_Rep; });
	ptnRaMassIncreasePlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.ptnRaMassIncrease_Rep; });
	ptnRaNitrogenMassIncreasePlt = std::accumulate(organDataFrame->h_gRchs.begin(), organDataFrame->h_gRchs.end(), 0.0f, [](float sum, const WheatRachis& curr) {return sum + curr.ptnRaNitrogenMassIncrease_Rep; });


}
/* END BLOCK
   ******** PLANT ORGANS MORPHOLOGY SUMMARY AND MASS DEMAND ********
*/


/* BLOCK
   ******** PLANT GAS EXCHANGE AND PHOTOSYNTHESIS ********
*/
//Z Coupled Photosynthesis and Stomatal Conductance Model
//  Farquhar-von Caemmerer-Berry (FvCB), Ball-Woodrow-Berry (BWB), leaf-level energy balance model 
void WheatPlant::calcGasExchange(const TWeather& weather, const TGasExSpeciesParam& photoparam)
{
	//const float tau = 0.50f;					//Z atmospheric transmittance, to be implemented as a variable => done
	const float LAF = 1.035f;					//Z leaf angle factor for rye grass leaf, Campbell and Norman (1998), Table 15.1, for rye (NOT ryegrass)
	const float leafwidth = 1.2f;				//Z greenleaf width compute during plant growth every hour
	const float atmPressure = 101.3f;			//Z kPa, to be predicted using altitude

	//Z GreenLfAreaPlt increases and decreases as plant growth and getting old
	//  LeafAreaPlt increases as plant growth but will not decrease, even leaf drops, leaf area counts
	float activeLeafRatio = this->greenLaAreaPlt / this->LaAreaPlt;

	//Z parameter to define single plant 2 field scale, i.e., how many plant in each field
	//  this is not "representative plant" but the whole field
	float plant2field = std::max(1.0f, initInfo.plantDensity * PlantLivingFraction);
	float field2plant = 1.0f / plant2field;

	//Z not all seeds germinate, not all germinated seed emerge, ......
	//  any accident can happen, so put "initInfo.plantDensity * PlantLivingFraction" as the real plant density
	//  10000.0 cm^2=1 m^2, which convert the plant leaf area unit from cm^2 to m^2
	float LAI = this->greenLaAreaPlt * plant2field / 10000.0f;

	//Z testing code
	if ((weather.doy <= 30) || (weather.doy >= 200)) { BiomassLeftover = 0.0f; }
	else { BiomassLeftover = (BiomassPool + BiomassReserve) / LAI; }

	unique_ptr<CGasExchange> sunlit = make_unique<CGasExchange>(photoparam);
	unique_ptr<CGasExchange> shaded = make_unique<CGasExchange>(photoparam);

	unique_ptr<CSolar> sun = make_unique<CSolar>();
	unique_ptr<CRadTrans> light = make_unique<CRadTrans>();

	//Z determine DOY and set it into the sun(CSolar) and light(CRadTrans) category
	Timer timer;
	int mm, dd, yy;
	timer.caldat(weather.jday, mm, dd, yy);
	int jday = timer.julday(1, 1, yy);
	sun->SetVal(weather.jday - jday + 1, weather.time, initInfo.latitude, initInfo.longitude, initInfo.altitude, weather.solRad);
	light->SetVal(*sun, LAI, LAF);

	sunlit_PFD = light->Qsl();		//Z CRadTrans object => mean photon flux density on sunlit leaves (umol m-2 s, unit converted from PAR w/m^2)
	shaded_PFD = light->Qsh();		//Z CRadTrans object => mean photon flux density on sunshade leaves 
	sunlit_LAI = light->LAIsl();	//Z CRadTrans object => LAI portion under sun shine
	shaded_LAI = light->LAIsh();	//Z CRadTrans object => LAI portion under sun shade

	float leaf_psi = weather.LeafWP;
	//??? cauwzj test 
	//leaf_psi = -0.05;
	//this->LeafNitrogenContentPlt = 1.5;
	//Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
	sunlit->SetVal(sunlit_PFD, initInfo, weather.airT, weather.CO2,
		weather.RH, weather.wind, atmPressure, leafwidth,
		leaf_psi, this->LaNitrogenContentPlt, BiomassLeftover);
	shaded->SetVal(shaded_PFD, initInfo, weather.airT, weather.CO2,
		weather.RH, weather.wind, atmPressure, leafwidth,
		leaf_psi, this->LaNitrogenContentPlt, BiomassLeftover);

	//Z Get immediate values from photosynthesis
	sunlit_A_net = sunlit->A_net;
	shaded_A_net = shaded->A_net;
	sunlit_A_gross = sunlit->A_gross;
	shaded_A_gross = shaded->A_gross;
	sunlit_gs = sunlit->gs;
	shaded_gs = shaded->gs;

	//plantsPerMeterSquare units are umol CO2 m-2 ground s-1 ;
	photosynthesis_gross = sunlit_A_gross * sunlit_LAI + shaded_A_gross * shaded_LAI;
	photosynthesis_net = sunlit_A_net * sunlit_LAI + shaded_A_net * shaded_LAI;

	//when outputting the previous step transpiration is compared to the current step's water uptake
	transpirationOld = transpiration;
	transpiration = 0;
	if (sunlit_LAI > 0) transpiration = sunlit->ET * sunlit_LAI;
	if (shaded_LAI > 0) transpiration += shaded->ET * shaded_LAI;
	// plantsPerMeterSquare units are grams per plant per hour;
	// Transpiration g (Water) plant^-1 hour^-1
	// Units of Transpiration from sunlit->ET/shaded->ET are mmol m-2 (leaf area) s-1
	// Calculation of transpiration from ET involves the conversion to gr per plant per hour 
	transpiration = transpiration * (60.0f * initInfo.timeStep) * 18.01f / 1000.0f * field2plant;

	//psi_l = (sunlit->get_psi()*sunlitLAI + shaded->get_psi()*shadedLAI)/LAI;
	this->VPD = sunlit->VPD;

	// photosynthesis_gross is umol CO2 m-2 leaf s-1
	// in the following we convert to g C plant-1 per hour
	// assimilate: grams CO2 per plant per hour
	assimilate = (photosynthesis_gross * CO2_MW * 1.0e-6f) * (60.0f * initInfo.timeStep) * field2plant;
	// photosynthesis_gross: grams biomass per plant per hour
	photosynthesis_gross = photosynthesis_gross * CH2O_MW * 1.0e-6f * (60.0f * initInfo.timeStep) * field2plant;
	// photosynthesis_net: grams biomass per plant per hour
	photosynthesis_net = photosynthesis_net * CH2O_MW * 1.0e-6f * (60.0f * initInfo.timeStep) * field2plant;

	if (LAI != 0.0f)
	{
		//Z leaf temperature
		temperature = (sunlit->Tleaf * sunlit_LAI + shaded->Tleaf * shaded_LAI) / LAI;
		//YY average stomatal conductance
		this->conductance = std::max(0.0f, ((sunlit_gs * sunlit_LAI + shaded_gs * shaded_LAI) / LAI));
	}
	else
	{
		temperature = sunlit->Tleaf;
		this->conductance = 0.0f;
	}
}
/* END BLOCK
   ******** PLANT GAS EXCHANGE AND PHOTOSYNTHESIS ********
*/

/* BLOCK
   ******** MASS OUTLETS AND DISTRIBUTIONS ********
*/
//Z maintainance respiration, output "g biomass plant^-1"
// based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
// units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
void WheatPlant::calcMaintRespiration(const TWeather& wthr)
{
	const float Q10 = 2.0f; // typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596 - could try 1.8
	float dt = initInfo.timeStep * DAYPERMINUTES;
	//	const double maintCoeff = 0.015; // gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61
	const float maintCoeff = 0.018f;
	float agefn = (this->greenLaAreaPlt + 1.0f) / (this->LaAreaPlt + 1.0f); // as more leaves senesce maint cost should go down, added 1 to both denom and numer to avoid division by zero. 
	//no maint cost for dead materials but needs to be more mechanistic, SK
	//agefn=1.0;
	float q10fn = pow(Q10, (wthr.airT - 20.0f) / 10.0f); // should be soil temperature or leaf or combination of use as --> (-stemMass*stem_coef) to reduce
	// total mass. Implement later after testing

//Z stem_coef is not used 
//double stem_coef = __min(1.0, this->DropLeafMassPlt / this->LeafMassPlt);

	maintRespiration = q10fn * maintCoeff * agefn * (this->ShootMass - this->dropLaMassPlt - this->dropShMassPlt) * dt;// gCH2O dt-1, agefn effect removed. 11/17/14. SK.
}

//Z biomass allocation
//  before doing any computation, C from photosynthesis must be assgined to pools
//  then two steps, 
//		first determine how much biomass is assigned to leaf (including sheath), internode (stem) and root
//		second determine the biomass rates for one unit leaf area, stem length
void WheatPlant::calsBiomassAllocation(const TWeather& wthr)
{
	float BiomassReserver2PoolFrac = 0.2f;	//Z during BiomassPool send mass out, BiomassReserve also wants to send 20% (refer to apsim wheat)
	float BiomassRootReleaseFrac = 0.8f;	//Z root biomass pool want to add 20% storage to root partition if biomass is sufficient

	// ********* step 1 ************************
	//Z organize the budge between reserviors and pools
	//  we will mainly argue the biomass split among the "BiomassSupply", "BiomassDemand", "BiomassPool" and "BiomassReserve".

	float b1 = 2.325152587f; // Normalized (0 to 1) temperature response fn parameters, Pasian and Lieth (1990)
	// Lieth and Pasian Scientifica Hortuculturae 46:109-128 1991
	// parameters were fit using rose data -
	float b2 = 0.185418876f; // I'm using this because it can have broad optimal region unlike beta fn or Arrhenius eqn
	float b3 = 0.203535650f;
	const float Td = 48.6f; //High temperature compensation point

	float g1 = 1.0f + exp(b1 - b2 * wthr.airT);
	float g2 = 0.0f;
	if (wthr.airT < Td) g2 = 1.0f - exp(-b3 * (Td - wthr.airT));

	float tmprEffect = g2 / g1; //Z must <1 since g2<1<g1
	float grofac = 1.0f / (5.0f * 60.0f / initInfo.timeStep); // translocation limitation and lag, assume it takes 1 hours to complete, 0.2=5 hrs

	BiomassSupply = 0.0f;

	//Z start to rebalance crop biomass budget
	//  always remember these
	//  BiomassPool ---- Checking Account
	//  BiomassReserve ---- Saving Account
	//  BiomassSupply ---- Withdraw

	//Z the first "if statement" handles maintrespiration
	//  the current maintrespiration is computed in "calcMaintRespiration"
	//  which is placed one step higher than the "calsBiomassAllocation", therefore, we can always have the current maintrespiration

	//Z feed some biomass from reservior to rapid pool
	BiomassPool += BiomassReserver2PoolFrac * std::max(BiomassReserve, 0.0f);
	BiomassReserve -= BiomassReserver2PoolFrac * std::max(BiomassReserve, 0.0f);

	//Z first "pay for" the maint respiration
	if (maintRespiration > 0.0f)
	{
		//Z short term pool can pay for the maintrespiration
		if (BiomassPool > maintRespiration)
		{
			BiomassSupply = maintRespiration;
			BiomassPool -= BiomassSupply;
		}
		//Z short term pool cannot pay for the maintrespiration
		//  but long term pool can pay for the maintrespiration
		else if (BiomassReserve > maintRespiration)
		{
			//Z in this block, we already assume that "BiomassPool <= maintRespiration"
			//  Thus, BiomassReserve is used to satisfy the maintainance respiration demand.
			//  Then BiomassPool is used to replenish BiomassReserve and then BiomassPool is set to 0. 
			//  In other words, BiomassPool is depleted first and then BiomassReserve is used to supplement whatever demand that's left
			BiomassSupply = maintRespiration;
			BiomassReserve -= BiomassSupply;
			BiomassReserve += BiomassPool;	// send remaining C (not enough to meet the demand) in shorterm pool to reserve
			BiomassPool = 0.0f;				// empty the BiomassPool
		}
		//Z both short term and long term pools cannot pay for the maintrespiration
		//  then combine them and pay as much as possible
		//  the worst case is that "BiomassPool+BiomassReserve" are exhausted
		else
		{
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0f;
			BiomassSupply = std::min(BiomassReserve, maintRespiration);
			BiomassReserve -= BiomassSupply;
		}
	}
	//Z no maintRespiration, then there is no BiomassSupply based on it
	else
	{
		BiomassSupply = 0.0f;
	}

	//Z the second "if statement" handles the biomass supply to heading, grain filling or relatied issues will also be added to BiomassSupply
	//  Since BiomassDemand==0 based on the current RYESIM design (never enter reproduction stages), 
	//  just reserve this function here
	//Z BiomassDemand does not enter into equations until grain fill
	//  upto this point BiomassSupply only has the possible maintrespiration
	//  or say, maintRepsiration is taken out from the short or/and long term pools
	if (BiomassDemand > 0.0f)
	{
		//Z pay for heading, grain filling or relatied issues using short term pool first and then use long term reserve
		//  but based on the payment, it is proportional to "BiomassPool" rather than "BiomassDemand"
		//  why use "__max(BiomassPool * tmprEffect * grofac, 0)"
		//    rather than "__max(BiomassDemand * tmprEffect * grofac, 0)"
		//  that is because the plant is RICH in biomass and generous
		if (BiomassPool > BiomassDemand)
		{
			float biomassSupplyGrain = std::max(BiomassPool * tmprEffect * grofac, 0.0f); //CADD from Grant
			BiomassSupply += biomassSupplyGrain;
			BiomassPool -= biomassSupplyGrain;
		}
		//Z if the short term pool BiomassPool not able to pay for the demand, 
		//  then we have to use the long term Reserve pool
		else if (BiomassReserve > BiomassDemand)
		{
			//Z first add BiomassPool to BiomassReserve
			//  BiomassPool always pay first as the short term pool, so it will be zeroed definitely
			//  so by the end, only BiomassReserve will have positive storage
			//  must be a positive storage? yes, because "BiomassReserve > BiomassDemand"
			BiomassReserve += BiomassPool; // BiomassPool negative here, add instead of subtract
			BiomassPool = 0.0f;

			//Z this supply is based on "BiomassDemand"
			//  that means the plant can give sufficient biomass but not that generous since it needs to use the saving account
			float biomassSupplyGrain = std::max(BiomassDemand * tmprEffect * grofac, 0.0f);
			BiomassSupply += biomassSupplyGrain;
			BiomassReserve -= biomassSupplyGrain;
		}
		//Z now neither short term and long term pool cannot satisfy the "BiomassDemand"
		else
		{
			//Z first add BiomassPool to BiomassReserve
			//  BiomassPool always pay first as the short term pool, so it will be zeroed definitely
			//  so by the end, only BiomassReserve will have positive storage
			//  must be a positive storage? yes, because "BiomassReserve > BiomassDemand"
			BiomassReserve += BiomassPool; // BiomassPool negative here, add instead of subtract
			BiomassPool = 0.0f;

			//Z this supply is based on "BiomassReserve"
			//  that means the plant may not be even able to supply based on "BiomassDemand"
			//  the plant should really careful about its biomass savings
			float biomassSupplyGrain = std::max(std::min(BiomassReserve, BiomassDemand) * tmprEffect * grofac, 0.0f); //CADD from Grant
			BiomassSupply += biomassSupplyGrain;
			BiomassReserve -= biomassSupplyGrain;
		}
	}
	//Z no reproductive growth and BiomassDemand, then there is no BiomassSupply based on it
	else
	{
		BiomassSupply += 0.0f;
	}

	//Z the third "if statement" handles the biomass supply to plant part, since "BiomassDemand" are only for grains
	//  we just need to determine if the plant is willing to give

	//Z may also need to put rootMassHere
	//  then total ptn mass will be root + shoot 
	float RootShootFractionController = std::min(0.90f, 0.5f + 0.5f * develop->get_TTd_Joint() / 1000.0f);
	float Yg = 1.0f; // synthesis efficiency

	//Z this is only for the vegetative part
	BiomassPltGrowth_ptn = std::max(BiomassShootGrowth_ptn, 0.0f) / RootShootFractionController / Yg;
	{
		//Z ideally after supplying the BiomassPool, 
		//  biomasspool is sufficient to pay for the shoot growth
		if (BiomassPool > BiomassPltGrowth_ptn)
		{
			//Z use BiomassPool to pay for max biomass demand 
			BiomassPool -= BiomassPltGrowth_ptn;
			BiomassSupply += BiomassPltGrowth_ptn;
		}
		//  biomasspool is not sufficient, but the overall biomass is good to pay for the shoot growth
		else if ((BiomassPool + BiomassReserve) > BiomassPltGrowth_ptn)
		{
			//Z plant almost has no short term storage to spend
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0f;
			BiomassSupply += BiomassPltGrowth_ptn;
			BiomassReserve -= BiomassPltGrowth_ptn;
		}
		//  not sufficient to pay for the demand
		//  shoot and leaf can still grow due to water pressure
		else
		{
			//Z plant almost has no storage to spend later
			//  if BiomassReserve<= 0 (possibly true), there may be even no C supply to plant growth
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0f;
			float BiomassSupply_PltGrowth = std::max(std::min(BiomassPltGrowth_ptn, BiomassReserve), 0.0f);
			BiomassSupply += BiomassSupply_PltGrowth;
			BiomassReserve -= BiomassSupply_PltGrowth;
		}
	}

	// ********* step 2 ************************
	//Z after get BiomassSupply, and shoot and root parts, 
	//  will need to distribute the values among plant organs
	//  we do not have a lot of growing stage, such as tassel..., so we just run until maturaty

	//Z first determine the supply, and assign 0;
	float shootPart_real = 0.0f;
	float rootPart_real = 0.0f;

	if (!develop->is_germinInit())
	{
		//Z no seed germinate yet (actually < 1% seeds germinate), pass and wait for the first seed germinating
		return;
	}
	else
	{
		//Z first take away the biomass demand part
		spikeletPart = max(min(BiomassDemand, BiomassSupply - maintRespiration), 0.0f);

		//Z biomass (g) partitioned to shoot and root after the first seed germinates
		//  why not emergence? because when seed germinates, there is already shoot/root growth
		//  start with 0.5,0.5 and then 
		shootPart = std::max(0.0f, Yg * (RootShootFractionController * max(BiomassSupply - maintRespiration - spikeletPart, 0.0f)));
		rootPart = std::max(0.0f, Yg * ((1.0f - RootShootFractionController) * max(BiomassSupply - maintRespiration - spikeletPart, 0.0f)));
		shootPart_real = shootPart;
		rootPart_real = rootPart;

		if (!develop->is_mature())
		{
			//Z in time step t-1, the value of pcrs is higher than that of pcrl, i.e., use more than expectation, or say already use some from pcrq
			if (wthr.pcrs > rootPart_old)
			{
				//Z the root biomass seems to be overcharged, however, there may exists some root storage,
				//  we add a test to see if the root growth can be paid by the root biomass storage before using the current shoot part
				//
				float rootMassOvercharge = wthr.pcrs - rootPart_old;
				//Z root biomass assignment is sufficient (at least combine with the BiomassRoot Pool), so return leftover to the biomassRoot
				//  then the biomassRoot will release a certion fraction (may be 20%) of the biomass storage to root growth
				if (BiomassRoot >= rootMassOvercharge)
				{
					BiomassRoot -= rootMassOvercharge;
					BiomassRoot = std::max(BiomassRoot, 0.0f);
					shootPart_real = shootPart;
					rootPart_real = rootPart + BiomassRootReleaseFrac * BiomassRoot;
					BiomassRoot = (1.0f - BiomassRootReleaseFrac) * BiomassRoot;
					rootPart_old = rootPart_real;
				}
				// root biomass storage is not sufficient
				// take some shoot part to pay for the root overcharge
				//Z always make the "rootPart_real" as the root growth values assigned at this time step. 
				else
				{
					// first let the shoot part pay for the overcharge
					if (rootPart > (rootMassOvercharge - BiomassRoot))
					{
						rootPart_real = std::max(rootPart - (rootMassOvercharge - BiomassRoot), 0.0f);
						shootPart_real = shootPart;
						BiomassRoot = 0.0f;
						rootPart_old = rootPart_real;
					}
					// second let the shoot part and root part pay for the overcharge
					else if (shootPart > (rootMassOvercharge - BiomassRoot - rootPart))
					{
						rootPart_real = 0.0f;
						shootPart_real = std::max(shootPart - (rootMassOvercharge - BiomassRoot - rootPart), 0.0f);
						BiomassRoot = 0.0f;
						rootPart_old = rootPart_real;
					}
					// eventually previous day overcharge is huge, stop growth
					else
					{
						//Z rootPart_old will leave a negative number
						//  need to be fed back in the following step
						//  during the fed back period, technically no root can grow
						//  the apparent root assignment will be (rootPart_real+BiomassRoot)=0
						rootPart_old = -rootMassOvercharge + BiomassRoot + shootPart + rootPart;
						shootPart_real = 0.0f;
						rootPart_real = 0.0f;
						BiomassRoot = 0.0f;
					}
				}
			}
			//Z root biomass assignment is sufficient, so return leftover to the biomassRoot
			//  then the biomassRoot will release a certion fraction (may be 20%) of the biomass storage to root growth
			else
			{
				float rootMassLeftOver = rootPart_old - wthr.pcrs;
				BiomassRoot += rootMassLeftOver;
				rootPart_real = rootPart + BiomassRootReleaseFrac * BiomassRoot;
				BiomassRoot = (1.0f - BiomassRootReleaseFrac) * BiomassRoot;
				rootPart_old = rootPart_real;
			}

			//Z fill the shoot part again
			if (shootPart_real >= shootPart)
			{
				//Z shoot growth is satisfied, do nothing
			}
			else if ((shootPart_real + BiomassPool) >= shootPart)
			{
				BiomassPool -= (shootPart - shootPart_real);
				shootPart_real = shootPart;
			}
			else if ((shootPart_real + BiomassPool + BiomassReserve) >= shootPart)
			{
				BiomassReserve -= (shootPart - shootPart_real - BiomassPool);
				shootPart_real = shootPart;
				BiomassPool = 0.0f;
			}
			else
			{
				shootPart_real += (BiomassReserve + BiomassPool);
				BiomassPool = 0.0f;
				BiomassReserve = 0.0f;
			}

			//Z partition carbon to plant organs
			leafPart = shootPart_real * max(BiomassLeafGrowth_ptn / BiomassShootGrowth_ptn, 0.0f);
			internodePart = shootPart_real * max(BiomassIntrNodeGrowth_ptn / BiomassShootGrowth_ptn, 0.0f);
			chaffPart = shootPart_real * max(BiomassChaffGrowth_ptn / BiomassShootGrowth_ptn, 0.0f);
			rachisPart = shootPart_real * max(BiomassRachisGrowth_ptn / BiomassShootGrowth_ptn, 0.0f);

			//Z judge if plant give leaf or internode too much
			if (leafPart <= BiomassLeafGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((leafPart - BiomassLeafGrowth_ptn) / Yg);
				leafPart = BiomassLeafGrowth_ptn;
			}

			if (internodePart <= BiomassIntrNodeGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((internodePart - BiomassIntrNodeGrowth_ptn) / Yg);
				internodePart = BiomassIntrNodeGrowth_ptn;
			}

			if (chaffPart <= BiomassChaffGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((chaffPart - BiomassChaffGrowth_ptn) / Yg);
				chaffPart = BiomassChaffGrowth_ptn;
			}

			if (rachisPart <= BiomassRachisGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((rachisPart - BiomassRachisGrowth_ptn) / Yg);
				rachisPart = BiomassRachisGrowth_ptn;
			}
		}
	}
	//??? cauwzj test if root biomass is correctly assigned
	//rootPart_real = 0.05;

	ActuralRootBiomassAssignment_PCRL = rootPart_real;
	ActuralShootBiomassAssignment = shootPart_real;
}

//Z Nitrogen allocation
//		first get the availiable nitrogen either from seeds or from plant root uptake
//		second assign the nitrogen to plant organs
void WheatPlant::calsNitrogenAllocation()
{
	//Z first determine the nitrogen quantity that can be allocated
	//  similar to the biomass allocation
	ShootNitrogenAvailiableAllocation = 0.0f;
	RootNitrogenAvailiableAllocation = 0.0f;
	leafPartNitrogen = 0.0f;
	internodePartNitrogen = 0.0f;
	chaffPartNitrogen = 0.0f;
	rachisPartNitrogen = 0.0f;
	spikeletPartNitrogen = 0.0f;

	if (!develop->is_germinInit())
	{
		return;
	}
	else
	{
		if (!develop->is_emerge())
		{
			//Z nitrogen solely comes from seed at 3.4%
			//  measured in mg (1000.0 converts g to mg)
			//  NitrogenPool += leafPart * MAX_N_PCT * 10.0;
			NitrogenPool = SeedNitrogenMass;

			//Z first root
			if (NitrogenPool >= ptnRtNitrogenGrowth)
			{
				RootNitrogenAvailiableAllocation = ptnRtNitrogenGrowth;
				RtNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool -= RootNitrogenAvailiableAllocation;
			}
			else
			{
				RootNitrogenAvailiableAllocation = NitrogenPool;
				RtNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool = 0.0f;
			}

			//Z second shoot
			if (NitrogenPool >= NitrogenShootGrowth_ptn)
			{
				ShootNitrogenAvailiableAllocation = NitrogenShootGrowth_ptn;
				NitrogenPool -= ShootNitrogenAvailiableAllocation;
			}
			else
			{
				ShootNitrogenAvailiableAllocation = NitrogenPool;
				NitrogenPool = 0.0f;
			}
		}
		else
		{
			//Z nitrogen comes from root uptake during this hour and the N released from plant aging
			NitrogenPool += (HourlyNitrogenSoilUptake + LaNitrogenReturnPlt + ShNitrogenReturnPlt
				+ InNitrogenReturnPlt + CfNitrogenReturnPlt + RaNitrogenReturnPlt);

			//Z first root
			if (NitrogenPool >= ptnRtNitrogenGrowth)
			{
				RootNitrogenAvailiableAllocation = ptnRtNitrogenGrowth;
				RtNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool -= RootNitrogenAvailiableAllocation;
			}
			else
			{
				RootNitrogenAvailiableAllocation = NitrogenPool;
				RtNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool = 0.0f;
			}

			//Z second kernel
			if (NitrogenPool >= NitrogenKernelGrowth_ptn)
			{
				spikeletPartNitrogen = NitrogenKernelGrowth_ptn;
				NitrogenPool -= NitrogenKernelGrowth_ptn;
			}
			else
			{
				spikeletPartNitrogen = NitrogenPool;
				NitrogenPool = 0.0f;
			}

			//Z second shoot
			if (NitrogenPool >= NitrogenShootGrowth_ptn)
			{
				ShootNitrogenAvailiableAllocation = NitrogenShootGrowth_ptn;
				NitrogenPool -= ShootNitrogenAvailiableAllocation;
			}
			else
			{
				ShootNitrogenAvailiableAllocation = NitrogenPool;
				NitrogenPool = 0.0f;
			}
		}
	}

	//Z second assign nitrogen to shoot organs
	if (!develop->is_germinInit())
	{
		return;
	}
	else
	{
		if (!develop->is_mature())
		{
			//Z partition carbon to plant organs

			//Z internode is growing, both leaves and internode receive nitrogen
			//  acceleration, leaf and sheath grow at the same speed
			//  pre-accel, sheath grow at 0.5 speed of the leaf grow
			if (NitrogenShootGrowth_ptn > 0.0f) 
			{ 
				leafPartNitrogen = ShootNitrogenAvailiableAllocation * std::max(NitrogenLeafGrowth_ptn / NitrogenShootGrowth_ptn, 0.0f); 
				internodePartNitrogen = ShootNitrogenAvailiableAllocation * std::max(NitrogenIntrNodeGrowth_ptn / NitrogenShootGrowth_ptn, 0.0f);
				chaffPartNitrogen = ShootNitrogenAvailiableAllocation * std::max(NitrogenChaffGrowth_ptn / NitrogenShootGrowth_ptn, 0.0f);
				rachisPartNitrogen = ShootNitrogenAvailiableAllocation * std::max(NitrogenRachisGrowth_ptn / NitrogenShootGrowth_ptn, 0.0f);
			}
			else 
			{ 
				leafPartNitrogen = 0.0f; 
				internodePartNitrogen = 0.0f;
				chaffPartNitrogen = 0.0f;
				rachisPartNitrogen = 0.0f;
			}
		}
	}
}

/* BLOCK
   ******** MASS-ASSIGNMENT TO ORGANS ********
*/
//Z after allocating the biomass, set the allocated biomass "AND N"
//  back to plant organs for the growth
void WheatPlant::calsSetMass()
{
	//Z first compute the assignment rate
	//  for example, we have the potential leaf area incremental and the biomass assigned for leaves,
	//  then, it will be easy to compute the leaf biomass per leaf area

	//Z note that we do not need to distinct the "acceleration" growth or not, 
	//  just use "leafPart/(leaf area incremental)" although some leafPart will go to sheath
	//  the leaf/sheath ratio will be adjusted in the RyeLeaf class
	if (ptnLfMassIncreasePlt > 0.0f) { LfBiomassAssignmentRate = max(leafPart / ptnLfMassIncreasePlt, 0.0f); }
	else { LfBiomassAssignmentRate = 0.0f; }
	if (ptnLfNitrogenMassIncreasePlt > 0.0f) { LfNitrogenAssignmentRate = max(leafPartNitrogen / ptnLfNitrogenMassIncreasePlt, 0.0f); }
	else { LfNitrogenAssignmentRate = 0.0f; }

	if (ptnInMassIncreasePlt > 0.0f) { InBiomassAssignmentRate = max(internodePart / ptnInMassIncreasePlt, 0.0f); }
	else { InBiomassAssignmentRate = 0.0f;}
	if (ptnInNitrogenMassIncreasePlt > 0.0f) { InNitrogenAssignmentRate = max(internodePartNitrogen / ptnInNitrogenMassIncreasePlt, 0.0f); }
	else { InNitrogenAssignmentRate = 0.0f; }

	if (ptnCfMassIncreasePlt > 0.0f) { CfBiomassAssignmentRate = max(chaffPart / ptnCfMassIncreasePlt, 0.0f); }
	else { CfBiomassAssignmentRate = 0.0f; }
	if (ptnCfNitrogenMassIncreasePlt > 0.0f) { CfNitrogenAssignmentRate = max(chaffPartNitrogen / ptnCfNitrogenMassIncreasePlt, 0.0f); }
	else { CfNitrogenAssignmentRate = 0.0f; }

	if (ptnRaMassIncreasePlt > 0.0f) { RaBiomassAssignmentRate = max(rachisPart / ptnRaMassIncreasePlt, 0.0f); }
	else { RaBiomassAssignmentRate = 0.0f; }
	if (ptnRaNitrogenMassIncreasePlt > 0.0f) { RaNitrogenAssignmentRate = max(rachisPartNitrogen / ptnRaNitrogenMassIncreasePlt, 0.0f); }
	else { RaNitrogenAssignmentRate = 0.0f; }

	if (ptnKrMassIncreaseSumPlt > 0.0f) { SpBiomassAssignmentRate = max(spikeletPart / ptnKrMassIncreaseSumPlt, 0.0f); }
	else { SpBiomassAssignmentRate = 0.0f; }
	if (ptnKrNitrogenMassIncreaseSumPlt > 0.0f) { SpNitrogenAssignmentRate = max(spikeletPartNitrogen / ptnKrNitrogenMassIncreaseSumPlt, 0.0f); }
	else { SpNitrogenAssignmentRate = 0.0f; }

	//Z assign biomass and N mass
	this->WheatMassDistributionSetup();
	this->WheatMassDistributionRun();

}
//Z Two steps to call wheat model to distribute mass and N among the organs
//  first update the mass allocation operators
//  second invoke the GPU operator to finish mass update
void WheatPlant::WheatMassDistributionSetup()
{
	//Z update leaf mass distribution operator
	organDataFrame->wLeafMassAssign.SetLeafMassAssignmentCondition(develop->is_accel(), this->LfBiomassAssignmentRate, this->LfNitrogenAssignmentRate);
	//Z update internode mass distribution operator
	organDataFrame->wIntrMassAssign.SetInternodeMassAssignmentCondition(this->InBiomassAssignmentRate, this->InNitrogenAssignmentRate);
	//Z update chaff mass distribution operator
	organDataFrame->wChafMassAssign.SetChaffMassAssignmentCondition(this->CfBiomassAssignmentRate, this->CfNitrogenAssignmentRate);
	//Z update rachis mass distribution operator
	organDataFrame->wRchsMassAssign.SetRachisMassAssignmentCondition(this->RaBiomassAssignmentRate, this->RaNitrogenAssignmentRate);
	//Z update spikelet mass distribution operator
	organDataFrame->wSpikMassAssign.SetSpikeletMassAssignmentCondition(this->SpBiomassAssignmentRate, this->SpNitrogenAssignmentRate);
}

void WheatPlant::WheatMassDistributionRun()
{
	//Z invoke the GPU for_each method
	//Z use for each method with zip iterator
	
	for (int ii = 0; ii < organDataFrame->h_gLeafNum; ii++) { organDataFrame->wLeafMassAssign.LeafMassRun(organDataFrame->h_gLeaf[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gIntrNum; ii++) { organDataFrame->wIntrMassAssign.IntrMassRun(organDataFrame->h_gIntr[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gChafNum; ii++) { organDataFrame->wChafMassAssign.ChafMassRun(organDataFrame->h_gChaf[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gRchsNum; ii++) { organDataFrame->wRchsMassAssign.RchsMassRun(organDataFrame->h_gRchs[ii]); }
	for (int ii = 0; ii < organDataFrame->h_gSpikNum; ii++) { organDataFrame->wSpikMassAssign.SpikMassRun(organDataFrame->h_gSpik[ii]); }
}
/* END BLOCK
   ******** MASS-ASSIGNMENT TO ORGANS ********
*/


void WheatPlant::calcRed_FRedRatio(const TWeather& weather)
// this function calculates an estimate of the Red to Far red light ratio from sunlit and shaded ET. This 
// ration is used to estimate the effects of plant density on leaf expansion and LAI. 
// A daily mean ratio is calculated. We use a 3 parameter sigmoid function to model the effect
{
	//double Xo=0.43, B=0.05, A=1.2;
	//double Xo=0.6, B=0.13, A=2.0; original
	//double Xo=0.9, B=0.43, A=2.0;
	float Xo = 0.85f, B = 0.65f, A = 2.0f;
	float dt = initInfo.timeStep * DAYPERMINUTES;
	float C2_effectTemp;
	//First set counter to 0 if it is the beginning of the day. 
	if (abs(weather.time) < 0.0001f)
	{// have to rename C2_effect to Light_effect
		//Zhu et al. Journal of Experimental Botany, Vol. 65, No. 2, pp. 641–653, 2014
		C2_effectTemp = exp(-(SunlitRatio - Xo) / B);
		C2_effect = std::min(1.0f, A / (1.0f + C2_effectTemp));
		develop->set_shadeEffect(C2_effect);
		SunlitRatio = 0.0f;
	}
	else
	{
		// calculate from emergence
		if (sunlit_LAI / (sunlit_LAI + shaded_LAI) > 0.05f) { SunlitRatio += sunlit_LAI / (sunlit_LAI + shaded_LAI) * dt; }
		else { SunlitRatio += 1.0f * dt; }
	}
}