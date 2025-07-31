#pragma once
#ifndef _WHEAT_PLANT_H_
#define _WHEAT_PLANT_H_

#include "WheatDevelopment.h"
#include "WheatTiller.h"
#include "WheatLeaf.h"
#include "WheatInternode.h"
#include "WheatChaff.h"
#include "WheatRachis.h"
#include "WheatSpikelet.h"
#include "WheatOrganDataFrame.h"
#include "gas_exchange.h"
#include "initinfo.h"
#include "weather.h"

#include "stdafx.h"
#include "radtrans.h"
#include "timer.h"

#include <iostream>
#include <string>
#include <memory>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <functional>

#include <chrono>

#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.00069444444f
#define CO2_MW 44.0098f
#define C_MW 12.011f
#define CH2O_MW 30.03f
#define MAX_N_PCT 4.0f
#define MIN_N_PCT 0.5f
#define MAX_N_PCT_ROOT 1.5f
#define cPFD 4.60f

#define MAXLEAFNUM 20
#define MAXINTRNUM 20
#define MAXTILLNUM 20
#define MAXCHAFNUM 1

using namespace std;

// Copy from lvalue containers, move from rvalue containers.
template<typename ...Cs>
auto zip(Cs... vecs) {
	std::vector<std::tuple<typename std::decay_t<Cs>::value_type...>> vec;

	auto len = std::min({ vecs.size()... });
	vec.reserve(len);
	for (std::size_t i = 0; i < len; ++i) {
		vec.emplace_back(std::move(vecs[i])...);
	}
	return vec;
};

class WheatPlant
{
public:
	WheatPlant(const TInitInfo&, TGasExSpeciesParam&);
	~WheatPlant();

	//Z control the morphological setting of the wheat plant
	void WheatPlantUpdate(const TWeather& weather);

	//Z for each wheat plant, this only return the mainstem
	WheatTiller* get_mainstem() { return mainstem.get(); }
	//Z get rye development
	WheatDevelopment* get_develop() { return develop.get(); }
	//Z get data frame pointer
	WheatOrganDataFrame* get_organDataFrame() { return organDataFrame.get(); }
	//Z return two pointers that point to a sunlit and a shaded leaf
	CGasExchange* get_sunlit() { return this->sunlit; }
	CGasExchange* get_shaded() { return this->shaded; }

	//********** Z Plant manupulate tiller update and summary **********************
	//Z tiller (including leaf/internode) will not do anything until plant call them to do it
	//  Recursive UPDATE and SUMMARY of leaf mass and internode mass, call by plant (upper class) and act on tiller (lower class)
	//  Plant class has authority on tiller class
	void WheatPlantMorphoUpdate(WheatTiller* wheatTiller);
	//Z  update the organ growth operators
	void WheatPlantSetupUpdate();
	//Z  summary the previous two steps "WheatPlantMorphoUpdate"+"WheatPlantSetupUpdate" and invoke GPU computaiton
	void WheatPlantUpdateRun();
	//Z sum the morphological update, data transfer from GPU to CPU
	void WheatUpdateSummary();
	//Z mass distribution among the tillers
	//  plant class (upper class) will call and authorize such action
	//  plant level mass and N budget will set up the mass assignment quantity
	//  then plant will call an mass assignment operator to complete that via GPU
	void WheatMassDistributionSetup();
	void WheatMassDistributionRun();

	//Z IO function for GET
//********** Z reorganize the code ****************
	float get_LfNum() { return LfNumPlt; }
	float get_greenLfNum() { return greenLfNumPlt; }
	float get_TlNum() { return TlNumPlt; }							//Z living tiller number
	float get_plantLivingFrac() { return PlantLivingFraction; }		//Z not all seeds are germinated or emerged, there is a fraction
	float get_LaArea() { return LaAreaPlt; }						//Z plant scale leaf lamina area in cm^2
	float get_greenLaArea() { return greenLaAreaPlt; }				//Z return green leaf lamina area, leaf area = green + aged (not dropped)

	float get_PltMass() { return WheatMass; }						//Z rye plant mass,   g/plant
	float get_PltTotalNitrogen() { return NitrogenMassPlt; }		//Z Total N in plant, mg/plant
	float get_grossPhotosynthesis() { return photosynthesis_gross; }//Z Gross Photosynthesis, grams biomass per plant per hour
	float get_netPhotosynthesis() { return photosynthesis_net; }	//Z Net Photosynthesis, grams biomass per plant per hour
	float get_ET() { return transpiration; }						//Z transpiration, gr per plant per hour 
	float get_ET_Old() { return transpirationOld; }					//Z previous transpiration, gr per plant per hour 
	float get_PltTmpr() { return temperature; }						//Z plant temperature
	float get_MaintenanceRespiration() { return maintRespiration; }	//Z maintenance respiration g biomass/plant/hour
	float get_stemMass() { return InMassPlt; }						//Z stem mass is the internode mass, since there may be multiple tillers, the existing part at this time step, not include potential growth
	float get_leafMass() { return LaMassPlt + ShMassPlt; }			//Z leaf mass per plant, lamina + sheath
	float get_shootMass() { return ShootMass; }						//Z leaf mass + sheath mass + internode mass
	float get_rootMass() { return RootMass; }						//Z root mass, need to be assigned from 2DSOIL

	float get_leafPart() { return leafPart; }						//Z leaf biomass assignment g/plant
	float get_actualRootBiomassAssignment_PCRL() { return ActuralRootBiomassAssignment_PCRL; }	//Z root biomass assignment g/plant, already adjusted with the BiomassRoot pool, different from maizsim (which done in crop)
	float get_actualShootBiomassAssignment() { return ActuralShootBiomassAssignment; }			//Z shoot biomass assignment g/plant

	float get_conductance() { return conductance; }										//Z leaf surface conductance, output from photosynthesis
	float get_VPD() { return VPD; }														//Z leaf surface VPD, output from photosynthesis
	float get_leafNitrogenMass() { return LaNitrogenMassPlt + ShNitrogenMassPlt; }		//Z plant scale leaf nitrogen mass, lamina + sheath mg/plant

	float get_cumulativeNitrogenDemand() { return CumulativeNitrogenDemand; }			//Z cumulated N demand to be assigned externally, "gN"/plant, 
	float get_cumulativeNitrogenSoilUptake() { return CumulativeNitrogenSoilUptake; }	//Z cumulated N uptake to be assigned externally, "gN"/plant, 
	float get_dropLaArea() { return dropLaAreaPlt; }									//Z plant scale drop leaf area, cm^2/plant

	//Z photosynthesis variables to be output
	float get_sunlit_LAI() { return sunlit_LAI; }
	float get_shaded_LAI() { return shaded_LAI; }
	float get_sunlit_PFD() { return sunlit_PFD; }
	float get_shaded_PFD() { return shaded_PFD; }
	float get_sunlit_A_net() { return sunlit_A_net; }
	float get_shaded_A_net() { return shaded_A_net; }
	float get_sunlit_A_gross() { return sunlit_A_gross; }
	float get_shaded_A_gross() { return shaded_A_gross; }
	float get_sunlit_gs() { return sunlit_gs; }
	float get_shaded_gs() { return shaded_gs; }
	float get_C2Effect() { return C2_effect; }
	float get_SunlitRatio() { return SunlitRatio; }
	float get_MaintRespiration() { return maintRespiration; }

	//Z IO function for biomass pools
	//  these function should be used for model evaluation purpose, maybe not that useful for model application
	//  biomass pools are more like internal variables
	float get_biomassRootPool() { return BiomassRoot; }				//Z root biomass storage, partially release to root biomass at current time g/plant
	float get_biomassReserve() { return BiomassReserve; }			//Z long term biomass reservior before allocation. 
	float get_nitrogenPool() { return NitrogenPool; }
	float get_HourlyNitrogenSoilUptake() { return HourlyNitrogenSoilUptake; }
	float get_LaNitrogenReturn() { return LaNitrogenReturnPlt; }
	float get_ShNitrogenReturn() { return ShNitrogenReturnPlt; }
	float get_InNitrogenReturn() { return InNitrogenReturnPlt; }
	float get_PltShootNitrogenAssignment() { return ShootNitrogenAvailiableAllocation; }
	float get_PltRootNitrogenAssignment() { return RootNitrogenAvailiableAllocation; }


	//cZ IO function for set
	void set_HourlyNitrogenDemand(float x) { HourlyNitrogenDemand = x; }
	void set_CumulativeNitrogenDemand(float x) { CumulativeNitrogenDemand = x; }
	void set_HourlyNitrogenSoilUptake(float x) { HourlyNitrogenSoilUptake = x; }
	void set_CumulativeNitrogenSoilUptake(float x) { CumulativeNitrogenSoilUptake = x; }
	void set_NitrogenRatio(float x) { NitrogenRatio = x; }

	void calcGasExchange(const TWeather& weather, const TGasExSpeciesParam& photoparam);
	void calcMaintRespiration(const TWeather&);
	void calsBiomassAllocation(const TWeather&);
	void calsNitrogenAllocation();
	void calsSetMass();
	void calsSetMorphology();
	void calcRed_FRedRatio(const TWeather&);

	float get_averagedBiomassLeftover() { return BiomassLeftover; }

private:
	TInitInfo initInfo;
	TGasExSpeciesParam gasExparam;
	//Z for Wheat plant, it only trace the mainstem,
	//  the plant will "own" this mainstem pointer, so use unique ptr
	unique_ptr<WheatTiller> mainstem;
	//Z wheat development object
	//  the plant will "own" this development pointer, so use unique ptr
	unique_ptr<WheatDevelopment> develop;
	//Z wheat organ data frame
	//  the data frame is owned by the plant module, so use unique ptr here
	unique_ptr<WheatOrganDataFrame> organDataFrame;

	//Z declare two pointers that point to a sunlit and a shaded leaf in the photosynthesis model
	CGasExchange* sunlit;
	CGasExchange* shaded;

	//Z all the following parts are for representative plants
	//Z some labels:
	// Lf = leaf
	// La = lamina
	// Sh = sheath
	// In = internode
	// Cf = chaff
	// Ra = rachis
	// Tl = tiller
	// Rt = root
	// Sp = spikelet
	// ptn = potential
	// Plt = plant
	//***** Plant Mass Group (g biomass) *********
	float WheatMass;
	float SeedMass, SeedNitrogenMass;
	float ShootMass;
	float RootMass;
	float PlantLivingFraction;		//Z  NOT ALL seeds become living plant
	//Z  used to adjust plant density, popslab ... 

//***** Plant Leaf Area Group (cm^2 or g) *********
//Z Dropped (dead) leaves must be Senescent leaves, so SenescentLfAreaPlt include dropped parts
//  Senescent parts may still stay on the plant, so Dropped leaf will not be all the senescent leaf
//  Thus, DropXXX is part of the total XXX
//Z These three leaf numbers are actually a statistical results, and it is actually sums of living fraction variable
//  e.g., leaf number is the total initial living fraction
//        green leaf number is the total current living fraction
//        drop leaf number is the difference between the two
	float LfNumPlt;
	float greenLfNumPlt;
	float dropLfNumPlt;

	float LaAreaPlt;
	float greenLaAreaPlt;
	float seneLaAreaPlt;
	float dropLaAreaPlt;	//Z leaf area dropped must be senecent leaf, 
	//  but senecent (part) of the leaf may not drop	
	//  so, drop leaf area is a subset of senescent leaf area

	float LaMassPlt;
	float ShMassPlt;
	float dropLaMassPlt;
	float dropShMassPlt;
	float LaNitrogenMassPlt;
	float ShNitrogenMassPlt;
	float LaNitrogenReturnPlt;
	float ShNitrogenReturnPlt;

	//Z for dead plant nitrogen mass, in dead organ and cannot be recycled, we further partition that into
	//  "sene": N in the senecent portion but not dropped from the stem
	//  "drop" N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     dead = sene + drop
	float deadLaNitrogenMassPlt;
	float deadShNitrogenMassPlt;
	float seneLaNitrogenMassPlt;
	float seneShNitrogenMassPlt;
	float dropLaNitrogenMassPlt;
	float dropShNitrogenMassPlt;

	float ptnLaMassIncreasePlt;
	float ptnShMassIncreasePlt;
	float ptnLfMassIncreasePlt;
	float ptnLaNitrogenMassIncreasePlt;
	float ptnShNitrogenMassIncreasePlt;
	float ptnLfNitrogenMassIncreasePlt;

	//***** Plant Internode Group (cm or g) *********
	int TlNumSingle;
	float TlNumPlt;
	float InLengthPlt;
	float InMassPlt;
	float InNitrogenMassPlt;
	float InNitrogenReturnPlt;
	float deadInNitrogenMassPlt;

	float ptnInMassIncreasePlt;
	float ptnInNitrogenMassIncreasePlt;

	//***** Plant Chaff Group (cm or g) *********
	float CfMassPlt;
	float deadCfMassPlt;
	float CfNitrogenMassPlt;
	float CfNitrogenReturnPlt;
	float deadCfNitrogenMassPlt;

	float ptnCfMassIncreasePlt;
	float ptnCfNitrogenMassIncreasePlt;

	//***** Plant Rachis Group (cm or g) *********
	float RaLengthPlt;
	float RaMassPlt;
	float deadRaMassPlt;
	float RaNitrogenMassPlt;
	float RaNitrogenReturnPlt;
	float deadRaNitrogenMassPlt;

	float ptnRaMassIncreasePlt;
	float ptnRaNitrogenMassIncreasePlt;

	//***** Plant Spikelet ( Flower + Kernel ) Group (cm or g) *********
	float KrNumPlt;
	float KrMassSumPlt;
	float KrNitrogenMassSumPlt;

	float ptnKrMassIncreaseSumPlt;
	float ptnKrNitrogenMassIncreaseSumPlt;

	//**** Plant Root N
	float RtNitrogenMassPlt;
	float ptnRtNitrogenGrowth;

	//***** Plant Based N group *********
	float NitrogenMassPlt;					//Z TotalNitrogen per plant, mg
	float NitrogenRatio;					//Z optimal N ratio according to N Dilution ratio	
	float LaNitrogenContentPlt;				//Z averaged leaf N content (only lamina because this is for photosynthesis)

	//***** PhotoSynthesis Group *********
	float sunlit_LAI, shaded_LAI;			//Z sunlit and shaded LAI values
	float sunlit_PFD, shaded_PFD;			//Z sunlit and shaded PFD (umol m-2 s) PFD Photon Flux Density
	float sunlit_A_gross, shaded_A_gross;	//Z leaf area based gross assmilation for one unit sunlit and shaded LAI
	float sunlit_A_net, shaded_A_net;		//Z leaf area based net assmilation for one unit sunlit and shaded LAI
	float assimilate;						//Z assimilation flux, gCO2/plant/timstep (hour)
	float photosynthesis_gross;				//Z gross photosynthesis (biomass version of the assimilation), unit converted from "umol CO2 m-2 s-1" to "g Biomass plant^-1 hour^-1"
	float photosynthesis_net;				//Z net photosynthesis, unit converted from "umol CO2 m-2 s-1" to "g Biomass plant^-1 hour^-1"
	float transpiration;					//Z current and previous values of transpiration, g/plant/hr
	float transpirationOld;
	float temperature;						//Z leaf temperature in the photosynthesis model
	float VPD;								//Z vapor pressure deficit
	float sunlit_gs, shaded_gs;				//Z sunlit and shaded stomatal conductance
	float conductance;						//Z averaged stomatal conductance

	//***** Carbon/Nitorgen Partitioning and Stroage (C in g; N in mg) *********

	float BiomassReserve;					//Z long term biomass pool in g per plant
	float BiomassPool;						//Z short term biomass pool in g per plant
	float BiomassSupply;					//Z biomass supply to plant growth in g per plant
	float BiomassDemand;					//Z biomass demand for REPRODUCTIVE GROWTH in g per plant
	float BiomassLeftover;					//Z BiomassPool+BiomassReserve and averaged over the green leaf area
	float BiomassRoot;						//Z biomass root pool in g per plant
	float BiomassLeafGrowth_ptn;			//Z potential (ideal) leaf (lamina + sheath) mass increases g per plant
	float BiomassIntrNodeGrowth_ptn;		//Z potential (ideal) internode mass increases g per plant
	float BiomassChaffGrowth_ptn;			//Z potential (ideal) chaff mass increases g per plant
	float BiomassRachisGrowth_ptn;			//Z potential (ideal) rachis mass increases g per plant
	float BiomassKernelGrowth_ptn;			//Z potential (ideal) kernel mass increases g per plant

	float BiomassShootGrowth_ptn;			//Z potential (ideal) leaf and internode mass increases g per plant
	float BiomassPltGrowth_ptn;				//Z potential shoot + root mass increases g per plant
	float ActuralRootBiomassAssignment_PCRL;
	float ActuralShootBiomassAssignment;

	float NitrogenPool;						//Z nitrogen pool in mg per plant

	float shootPart;						//g per plant Carbohydrate partitioined to shoot
	float rootPart;							//g per plant Carbohydrate partitioned to root
	float shootPart_old;					//shoot/root in the previous time step
	float rootPart_old;

	float leafPart;							//g per plant biomass partitioned to leaf (in ryesim, this include sheath)
	float internodePart;
	float chaffPart;
	float rachisPart;
	float spikeletPart;

	float leafPartNitrogen;
	float internodePartNitrogen;
	float chaffPartNitrogen;
	float rachisPartNitrogen;
	float spikeletPartNitrogen;

	float LfBiomassAssignmentRate;
	float InBiomassAssignmentRate;
	float CfBiomassAssignmentRate;
	float RaBiomassAssignmentRate;
	float SpBiomassAssignmentRate;

	float maintRespiration;

	float HourlyNitrogenDemand;				// Nitrogen demand in mg N plant-1
	float CumulativeNitrogenDemand;			// cumulativeNitrogen demand in mg N plant-1
	float HourlyNitrogenSoilUptake;			// Nitrogen uptake from the soil mg N plant-1
	float CumulativeNitrogenSoilUptake;		// Nitrogen uptake from the soil mg N plant-1

	float NitrogenLeafGrowth_ptn;
	float NitrogenIntrNodeGrowth_ptn;
	float NitrogenChaffGrowth_ptn;
	float NitrogenRachisGrowth_ptn;
	float NitrogenKernelGrowth_ptn;

	float NitrogenDemand;					// Kernel part, for reproductive growth
	float NitrogenShootGrowth_ptn;			// Vegetative part Based on the potential increase of the plant mass (leaf and internode), there should be a potential N requirement in mg N plant^-1
	float ShootNitrogenAvailiableAllocation;
	float RootNitrogenAvailiableAllocation;

	float LfNitrogenAssignmentRate;
	float InNitrogenAssignmentRate;
	float CfNitrogenAssignmentRate;
	float RaNitrogenAssignmentRate;
	float SpNitrogenAssignmentRate;

	//***** Plant Information *********
	float sowingDay;
	float age;
	float emerge_gdd;						//records thermal time needed for plant to emergy
	float SunlitRatio;						// daily ratio of sunlit leaf area to use for scaling leaf expansion due to carbon stress
	float C2_effect;

	//Z set some plant update stages regarding to flower
	bool flowerPrimInit;
	bool anthesisStart;
	bool anthesisEnd;


};

#endif