#pragma once

#ifndef _WHEAT_TILLER_H_
#define _WHEAT_TILLER_H_

#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>

#include "weather.h"
#include "initinfo.h"
#include "WheatThermalTime.h"
#include "WheatDevelopment.h"
#include "WheatInterNode.h"
#include "WheatLeaf.h"
#include "WheatSpikelet.h"
#include "WheatRachis.h"

#include "WheatOrganDataFrame.h"

#define MINUTESPERDAY 1440.0f
#define DAYPERMINUTES 0.00069444444f
#define DAYPERHOUR 0.041666666667f
#define PHYLLOCHRON 106.0f
#define SPIKEPERPHYLLOCHRON 9.3f	//Z number of spikelet primordia per phyllochron
#define TTDBASELFLORET 9.1f			//Z Gdd required for basal floret pair fertilization
#define TTDNONBASELFLORET 27.3f		//Z Gdd required for non-basal floret fertilization
#define CO2_MW 44.0098f
#define C_MW 12.011f
#define CH2O_MW 30.03f
#define MAX_N_PCT 4.0f
#define MIN_N_PCT 0.5f

#define MAXLEAFNUM 20
#define MAXINTRNUM 20
#define MAXTILLNUM 20
#define MAXCHAFNUM 1
#define MAXSPIKNUM 30	//Z max spikelet per tiller: each spikelet will have max flower number

//Z leaf shape parameters from fspm
//Z time (day) difference for the next primordia
#define PLASTOCHRONE  10.34167f	//Z day (= 76.1 / 12 * 24 * 3600 sec) Leaf plastochron(sec at 12°c but casted to "day" in this study) calculated from Ljutovac 2002 with primordia of 5E-5 m(76 dd); Malvoisin 35dd associated with init 3E-5 m
//Z time (day) difference between leaf and internode from the same position
#define TIMEDIFF_LF_INTR 31.7083333f	//Z day (=2739600 sec) of equivalent time between leaf and internode initiation

using namespace std;

class WheatTiller
{
public:
	WheatTiller(int rank, int order, int cumurank, bool mainstem, WheatDevelopment* dv, WheatOrganDataFrame* ddf, float livingFrac);
	~WheatTiller();
	//Z update the growth of the current tiller, 
	// based on single tiller structure
	// and then upscaling to representative values 

	//Z tiller will only care about its own morphology
	//  that is to say, tiller prepare information for its sub-organ updates,
	//  while the real update is in the wheat plant module, where the gpu function is invoked.
	void WheatTillerSingleMorph(void);
	void WheatTillerSingleDeath(void);
	void WheatTillerSingleUpdate(void);

	//Z global update based on a data-oriented programming
	//  embed in the wheatplant class, because wheatplant class is a "upper class" and wheattiller class is a "lower class"
	//  plant class has the authority to control the tiller class
	//  while tiller class can control (get_information from) even lower classes, such as leaf, internode, etc.

	//***** Ordinary IO function for tiller-leaf-internode *********
	float get_physAge() { return physAge; }
	float get_FlagLfPhyllchron() { return FlagLfPhyllchron_Cplt; }
	WheatTiller* get_subtiller(int ii) { return SubTiller[ii].get(); }

	//***** mainstem development *********
	//Z get order, start from 0 (mainstem)
	int get_TlOrder() { return order; }
	bool is_mainstemInitilized() { return (mainstemInitializaiton > 1); }
	bool is_living() { return living; }
	bool is_force_to_death_current_step() { return force_to_death_current_step; }
	bool is_FlagLfPhyllchron_Init() { return FlagLfPhyllchron_Init; }
	bool is_death2fnalize() { return death_2_finalize; }
	void set_death2fnalize_true() { death_2_finalize = true; }
	void MainstemInitialize();

	//Z get organ number and suborgan index
	int get_StNum() { return SubTillerNum; }

	//Z tiller living fraction
		//  can be used as "representative" tiller number if living, not necessarily integer because it report a statistical value
	float get_tillerLivingFrac() { return livingFrac; }
	//  should be a boolean function, but want to use it as an integer number of living tiller in a plant
	int is_tillerLiving()
	{
		if (living) { return 1; }
		else { return 0; }
	}
	void set_tillerLivingFrac(float x)
	{
		if (!living) { livingFrac = 0.0f; }
		else { livingFrac = x; }
	}

private:

	//Z trace the main plant development status
	//  "SHALLOW COPY" so should not be the smart pointer.
	WheatDevelopment* develop;

	//Z trace the main plant organ data frame pointer
	//  "SHALLOW COPY" so should not be the smart pointer.
	WheatOrganDataFrame* organDataFrame;

	//-----------------------------------------------------------------------------------------------------
	//Sub-organ Group
		//Z LINKED-LIST for tiller hierarchy
		//  tiller structure is build in CPU, so do not put subtiller into the global space.
		//  this tiller is only for itself, and index to the suborgans, e.g., leaves, internodes, ...
	int SubTillerNum;
	vector<unique_ptr<WheatTiller>> SubTiller;

	//Z Internodes are derived from the tiller
	//  internode data in device
	int InterNodeNum;
	int SubInterNodeIdx[MAXINTRNUM];

	//Z Leaves are derived from the tiller
	//  leaf data in device
	int LeafNum;
	int EmergeLfNum;
	int GreenLfNum;
	int DropLfNum;
	int SubLeafIdx[MAXLEAFNUM];

	//Z Chaff is for each tiller (this is only a model simplification)
	//  chaff should be surround the seeds, but we only assume 1 chaff per tiller
	//  and only chaff mass (and N) is considered, no geometrical shape
	int ChaffNum;  //Z 0 or 1
	int SubChaffIdx;

	//Z Rachis is for each tiller (this is only a model simplification)
	//  Rachis should be surround the seeds, but we only assume 1 Rachis per tiller
	//  and only Rachis mass (and N) is considered, no geometrical shape
	int RachisNum;  //Z 0 or 1
	int SubRachisIdx;

	//Z Spikelet for the tiller
	//  each tiller can have upto 30 spikelets, and each spikelet has multiple flowers, and flowers has grains
	int SpikNum;
	int SubSpikIdx[MAXSPIKNUM];

	//-----------------------------------------------------------------------------------------------------
	//Sub-organ Group

	//Z position of this tiller
	int rank;						//Z number of the tiller at its parent tiller (maybe mainstem) counted from bottom, start from 1
	int order;						//Z the order of tiller, 0 for mainstem, 1 for the tiller branching from the mainstem, 2 for the tiller branching from the order 1 tiller...
	int cumurank;					//Z cumulative tiller rank and order, mainstem is 0, first tiller of the mainstem is 0+1, first subtiller of the first tiller is 0+1+1
	bool mainstem;					//Z if this tiller is the mainstem or not
	bool mainstemInitializaiton;	//Z a special case for mainstem, immediately grow two leaves, so two steps

	//Z single ridge killer
	bool singleRidge;				// the singleRidge stage is reached for the whole plant
	bool born_after_singleRidge;	// this tiller is born after singleRidge, will not have this single ridge killer issue
	bool terSpikeInit;				// the termination point for spikelet emerging
	bool flowerPrimInit;			// the initiation of flower primordium formation
	bool flowerPrimInit_already;	// already flower primordium formation, initialized as false
	bool jointing;					// the jointing stage is reached
	bool jointing_already;			// combine with "jointing", s.t. tiller death operation on that date only called/determined once
	bool jointing_kill;				// this tiller should be killed, tiller at jointing with leaf number < 4 or tiller branched after single ridge will be killed
	bool anthesisStart;				// beginning of anthesis, flower initiation/abortion stops, real flower should grow up
	bool anthesisStart_already;     // already anthesis
	bool anthesisEnd;				// end of anthesis, which flower will be the actual flower will be determined at this point

	//Z internode growth calculator
	bool elongationStart;			// receive start elongation from the plant and development class
	bool elongationStart_already;	// if elongation already occurs, combine with "elongationStart" to trace the first entry of elongation
	int elongation_first;			// the first internode for elongation
	int elongation_last;			// the last internode for elongation
	int rank_2_elongation;			// mark the rank of internode, which will elong at the current time step


	//Z living fraction and force to die variables
	bool living;						// reserve for counting tiller death, living=1, death=0
	float livingFrac;					// livingFrac of the tiller, 
	float livingFrac_ini;				// livingfrac when this tiller is initialized
	bool force_to_death_current_step;	// marker for tiller death due to external/ambient reasons, force to die, can only be true on one step when tiller death function is called
	bool death_2_finalize;				// a computational marker, when tiller is killed or dead, use this to update its leaf and internode for the last time

	//Z time order to make growth progress
	//  gdd/equivalent time based
	float gddpBranch;						//Z gdd number (from planting) when the tiller is branched
	float physAge;							//Z gdd age (in reference to endGrowth and lifeSpan, days)
	float Cur_TTd;							//Z Thermal time for this step (1 hour)
	float TlPseudoAge;						//Z the equivalent time age of the tiller
	float Cur_TEq;							//Z the equivalent time for each time step
	float TillerLfPhyllchron[MAXLEAFNUM];	//Z leaf initiation array based on gdd
	float TillerLfPlastoChrone[MAXLEAFNUM];	//Z leaf initiation array based on equivalent time
	float cumuIntrLength[MAXINTRNUM];		//Z these two (the difference) are used to compute the emergence
	float cumuLigulationHeight[MAXLEAFNUM];

	bool FlagLfPhyllchron_Init;				//Z flag leaf initiate or not
	float FlagLfPhyllchron_Cplt;			//Z flag leaf complete

	//Z spike initiation and spike growth parameters
	float phyAge_SpikInit;					//Z cumulated gdd since single ridge to determine spikelet initiation
	float phyAge_FlwFert;					//Z cumulated gdd for that tiller determining flower fertilization positions
	float gddSpike;							//Z parameter, gdd need for each 
	int SpikAdd;							//Z this step add one spikelet, should only be 0 or 1
	int FlwFertAdd;							//Z this step fertilize one spikelet's flower, based on basel flower number
	//Z spikelet that can initiate flower prim
	int spiint;								//Z initial rank of spikelet that can have flower initiation
	int spifrt;								//Z first fertilization spikelet
	int spitop;								//Z top rank of the spikelet that can have flower initiation
	int spibot;								//Z bottom rank of the spikelet that can have flower initiation
	int ferint;								//Z initial fertilization spikelet
	int ferbot;								//Z bottom fertilization spikelet
	int fertop;								//Z top fertilization spikelet

	//***** Environmental Factors *********
	float StressTiller;						//Z Tiller water and nitrogen stress, 1 means no stress and 0 means full stress
	float StressTillerMin;					//Z the minimal stress value ever recorded
	float TTd_since_singleRidge;
	float TTd_since_jointing;
	float N_effect;
	float water_effect;
	float shade_effect;
	float tmpr_effect;
	float tmpr_effect_terminal;
	float Cold_Time;

};

#endif