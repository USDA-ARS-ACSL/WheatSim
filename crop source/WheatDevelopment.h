#pragma once
#ifndef _WHEAT_DEVELOPMENT_H_
#define _WHEAT_DEVELOPMENT_H_
#include "weather.h"
#include "initinfo.h"
#include "WheatThermalTime.h"
#include <iostream>
#include <string>
using namespace std;

struct TEvent
{
public:
	TEvent() { daytime = 0.0f;  done = false; }
	float daytime;
	bool done;
};

class WheatDevelopment
{
public:
	WheatDevelopment(const TInitInfo&);
	~WheatDevelopment();
	void WheatDelpUpdate(const TWeather& wthr);
	void PlantDie(float frac);
	float EmCurv(float xval, int phase);
	float ReCurv(int phase, float yval);
	void KernelFitCurv(int opt, float* arr);

	//Z IO functions
	WheatThermalTime get_WheatTTd() { return WheatTTd; }
	float get_TmprCur() { return T_cur; }		//Z get current temperature (maybe airT or soilT)
	float get_Tair() { return T_air; }		//Z get air temperature
	float get_Tave() { return T_ave; }
	float get_TimeStep() { return TimeStep; }
	float get_DayTime() { return DayTime; }
	float get_ColdTime() { return Cold_Time; }
	float get_ColdTimeRatioPlant() { return Cold_Time_Ratio_Plant; }
	float get_ColdTimeRatioJoint() { return Cold_Time_Ratio_Joint; }

	float get_TTdCur_dssat() { return Cur_TTd; }	//Z incremental of Gdd for this time step
	float get_TTd_Plant() { return TTd_plant; }
	float get_TTd_Joint() { return TTd_joint; }
	float get_TTd_Elong() { return TTd_2_elongation; }
	float get_TTd_FlagLf_min() { return TTd_FlagLf_min; }  //Z the min flag leaf complete gdd, to compute elongation internodes
	float get_TTd_sinceSingleRidge() { return TTd_since_singleRidge; }
	float get_TTd_sinceJointing() { return TTd_since_jointing; }
	float get_DayTime_kernelFertStart() { return kernelFertStart.daytime; }
	float get_DayTime_afterAntss() { return DaysAfterAntss; }

	float get_TimeTeqCur_fspm() { return TimeTeqCur; }

	//Z IO for effects group
	float get_PredawnLWP() { return PredawnLWP; }
	float get_shadeEffect() { return ShadeEffect; }
	float get_plantLivingFraction() { return PlantLivingFraction; }
	float get_pltLivingBiomass() { return pltLivingBiomass; }
	float get_pltLivingNitrogenMass() { return pltLivingNitrogenMass; }
	float get_pltTillerNum() { return pltTillerNum; }

	int get_doy() { return doy; }

	//Z get kernel weight and duration fitting results
	float get_kwt00() { return fwtc[0]; }
	float get_kwt01() { return fwtc[1]; }
	float get_kwt10() { return fwtc[2]; }
	float get_kwt11() { return fwtc[3]; }
	float get_kdur00() { return durc[0]; }
	float get_kdur01() { return durc[1]; }
	float get_kdur10() { return durc[2]; }
	float get_kdur11() { return durc[3]; }

	//Yes-No stage functions
	bool is_germinInit() { return germinInit.done; }			//Z if germination starts, emergence may come before germination completes.
	bool is_germination() { return germination.done; }			//Z if germinated, yes then start root initialization
	bool is_emerge() { return emergence.done; }					//Z if emerged, yes then start plant growth
	bool is_singleRidge() { return singleRidge.done; }			//Z if single-ridge stage is reached, yes then there will be some tiller forced to die (leaf number <4)
	bool is_doubleRidge() { return doubleRidge.done; }			//Z if double-ridge stage is reached, yes then "spikelet" will start to be initiated
	bool is_flowerPrimInit() { return flowerPrimInit.done; }	//Z if the flower primordium initiation begins
	bool is_terSpikeGrow() { return terSpikeGrow.done; }		//Z if terminal spikelet is reached,
	bool is_startEnlongation() { return elongationStart.done; }	//Z if the elongation starts
	bool is_accel() { return acceleration.done; }				//Z if single ridge is reached and growth acceleration occurs
	bool is_startJoint() { return jointStart.done; }			//Z if jointing stage is triggered 
	bool is_startHead() { return headStart.done; }				//Z if heading stage is triggered
	bool is_startAnthesis() { return anthesisStart.done; }		//Z flower start to form
	bool is_endAnthesis() { return anthesisEnd.done; }			//Z flower formation is ended
	bool is_grainFillBegan() { return grainFillBegan.done; }	//Z start the grainfilling, i.e., kernel should grow now
	bool is_mature() { return maturity.done; }					//Z kernel grow end, everything should go away
	bool is_dead() { return death.done; }						//Z trigger the whole program to end

	//Z chaff growth stages
	bool is_chaffStart() { return chaffStart.done; }
	bool is_chaffEnd() { return chaffEnd.done; }
	float get_DurationChaffGrowth() { return DurationChaffGrowth; }

	//Z rachis growth duration
	float get_DurationRachisGrowth() { return DurationRachisGrowth; }

	//Z Set stage functions
	void set_maturity(bool d, float daytime) { maturity.done = d; maturity.daytime = daytime; }
	void set_death(bool d, float daytime) { death.done = d; death.daytime = daytime; }
	void set_shadeEffect(float x) { ShadeEffect = x; }
	void set_pltLivingBiomass(float x) { pltLivingBiomass = x; }
	void set_pltLivingNitrogenMass(float x) { pltLivingNitrogenMass = x; }
	void set_pltTillerNum(float x) { pltTillerNum = x; }
	//Z this function will be called by some tiller whose first do the kernel fertilization
	void set_kernelFertStart(void) {
		if (!kernelFertStart.done)
		{
			kernelFertStart.daytime = DayTime;
			kernelFertStart.done = true;
		}
	}

private:

	WheatThermalTime WheatTTd;

	//****** GDD numbers and growing stage **************
	float Cur_TTd;												//Z gdd at this time step
	float TimeTeqCur;
	//Z germination start > 0.01 (small number)
	TEvent germinInit;
	//Z germination rate > 0.5, the whole germination will be 14 days and overlap the emergence time
	TEvent germination;
	//Z emergence rate > 0.01 (small number)
	TEvent emergence;
	//Z single ridge
	TEvent singleRidge;
	//Z double ridge
	TEvent doubleRidge;
	//Z flower primordium initiation begins
	TEvent flowerPrimInit;
	//Z terminal spikelet growth stages, tss
	TEvent terSpikeGrow;
	//Z elongation stage start, for the whole plant
	TEvent elongationStart;
	//Z jointing
	TEvent jointStart;
	//Z booting: flag leaf fully grown
	TEvent bootStart;
	//Z heading
	TEvent headStart;
	//Z beginning of anthesis
	TEvent anthesisStart;
	//Z end of anthesis
	TEvent anthesisEnd;
	//Z grainFillingBegin
	TEvent grainFillBegan;
	//Z maturity
	TEvent maturity;
	//Z death
	TEvent death;
	//Z acceleration is not a stage, but we treat it as one since it triggers faster growth
	TEvent acceleration;
	//Z flower fertilization start, also means the kernel will initiate
	//  the first day of kernel fertilization (kernel starts to grow)
	//  should be set externally
	TEvent kernelFertStart;

	//Z chaff start and end
	TEvent chaffStart;
	TEvent chaffEnd;

	//****** Environment Conditions **************
	bool NorthHemisphere;	//Z true if north hemisphere, false if south hemisphere
	float T_cur, T_air, T_ave;
	//Z this predawn leaf water potential is stored but not processed in this class, save and offer it to leaf class
	float PredawnLWP;
	//Z shaded effects, compute in plant class, assigned to development and stored, used in leaf class
	float ShadeEffect;
	//Z time step
	float TimeStep;
	float DayTime;
	float DaysAfterAntss;

	//****** Germination & Emergence Group **************

	//soil water condition
	int temcon;				//Z current water content status;
	int condit;				//Z soil condition class
	int	oldcon;				//Z old soil condition class
	int doy;				//Z record doy for future usage;
	int doyRecord;          //Z record the DOY s.t. ready to update daily soil water condition

	//Z germination rate;
	float SeedGerminationRate;
	//Z emergence fraction and its real (truncated) fraction to be used outside;
	float EmergFrac, FracEmerg_Real;
	//Z germination fraction and its real (truncated) fraction to be used outside;
	float GerminFrac, FracGermin_Real;
	//Z seed number should be input, since there exist emergence fraction, so actual plant number changes;
	int SeedNum, GerminateSeedNum, PlantNum;
	float PlantLivingFraction;

	//Z GDD for germination and emergence, TTd_seed will be from sowing.
	float TTd_seed, Cur_Germin_TTd, Cur_Emerg_TTd;

	//Z germination and emergence truncation limits
	float upcutg[3], locutg[3], upcute[3], locute[3];
	//Z mean elong rate(gdd / cm) * planting depth(cm) + mean germination rate(gdd) = memerg
	//     sgerm, semerg, std of germination and emergence gdd
	//     emercv, the covariance coefficients
	float melrat[3], mgerm[3], memerg[3], sgerm[3], semerg[3], emercv;

	//****** Growth Stage Till Jointing **************
	//Z gdd from doy==1 to till Jointing
	float TTd_joint;
	//Z gdd for plant growth, from emergence
	float TTd_plant;
	//Z gdd min for the earliest time that the flag leaf will complete its growth (MINGDF in shootgro)
	//  This is used in determining exactly when the flag leaf will complete growth (GDDF)
	//  and which will be the last internode to undergo accelerated elongation(ILAST, the peduncle).
	float TTd_FlagLf_min;
	float TTd_FlagLf_min_aux;	//Z auxiliary variable, the last leaf will appear within .5 phyllochron of boot, and will thus complete growth no earlier than .5 phyllochron after boot.
	float TTd_2_singleRidge;	//Z gdd to single ridge based on TTd_joint
	float TTd_since_singleRidge;//Z gdd since single ridge start
	float TTd_2_doubleRidge;	//Z gdd to double ridge based on TTd_joint
	float TTd_2_flowerPrimInit; //Z gdd to flower primordium initiation begins
	float TTd_2_elongation;		//Z gdd to elongation based on TTd_joint
	float TTd_2_terSpikeGrow;	//Z gdd to terminal spikelet growth stages based on TTd_joint
	float TTd_2_jointStart;		//Z gdd to jointing stages based on TTd_joint
	float TTd_since_jointing;	//Z good since jointing stage is reached
	float TTd_2_bootStart;		//Z gdd to booting stages based on TTd_joint
	float TTd_2_headStart;		//Z gdd to heading stages based on TTd_joint
	float TTd_2_anthesisStart;	//Z gdd to the beginning of anthesis based on TTd_joint
	float TTd_2_anthesisEnd;	//Z gdd to the end of anthesis based on TTd_joint
	float TTd_2_maturity;		//Z gdd to maturity
	
	//****** Plant Death Factors **************
	//     4 fractions to mark the kill percentage range
	//     Cumulate Present value to avoid 
	float LowKill_Frac, HighKill_Frac, PresentKill_Frac, TotalKill_Frac;

	//****** cumulate the cold temperature time that the plant experienced
	float Cold_Time;
	float Cold_Time_Ratio_Plant;
	float Cold_Time_Ratio_Joint;

	// chaff growth duration
	float TTd_2_chaffStart;			//Z gdd half after the double ridge
	float TTd_2_chaffEnd;			//Z gdd to half maturity
	float DurationChaffGrowth;		//Z chaff growth duration

	//Z rachis growth
	float DurationRachisGrowth;     //Z rachis growth duration

	//****** some plant level morphology and mass data for tiller to call *******
	float pltLivingBiomass;
	float pltLivingNitrogenMass;
	float pltTillerNum;

	//****** parameters for grain and kernel *******
	const float kertem[2][2] = { {25.0f,36.0f},{25.0f,36.0f} };		//Z Temperature breakpoints and cutoffs for kernal growth curves
	const int kerpts = 6;
	const float kerpot[6][3] = {
		{10.0f,60.0f,790.0f},
		{15.0f,55.0f,745.0f},
		{19.6f,51.0f,720.0f},
		{25.0f,43.2f,600.0f},
		{30.0f,27.0f,500.0f},
		{35.0f,13.5f,400.0f} };
	float fwtc[4] = { 0.0f,0.0f,0.0f,0.0f };		//Z kernel weight computation coefficients, 0,1 => first line, 2,3 => second line
	float durc[4] = { 0.0f,0.0f,0.0f,0.0f };		//Z kernel duration computation coeffcients

};


#endif