
//Class Controller for rye
//
// RyeController.h
//
//Based on CPM simulation_controller
#pragma once
#ifndef _RYE_CONTROLLER_H_
#define _RYE_CONTROLLER_H_
#include "timer.h"
#include "WheatDevelopment.h"
#include "WheatPlant.h"
#include "weather.h"
#include "initinfo.h"
#include "gas_ex_species_param.h"

#include <memory>


#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define MINUTESPERDAY (24*60);
#endif

class WheatController
{
public:
	WheatController(const char*, const char*, const char*, TInitInfo);
	~WheatController();

	//Z read files and prepare the initilization
	void initialize();

	//Z read files from 2dsoil module
	void readWeatherFrom2DSOIL(const TWeather& wthr);

	//Z output to the crop results g01
	void outputToCropFile();

	int run(const TWeather& wthr);

	//Z IO functions
	WheatPlant* get_plant(void) { return wheatPlant.get(); }
	TInitInfo* get_initInfo(void) { return &initInfo; } //Z alhtough initInfo is defined as a structure, we return a pointer to avoid copy the whole structure

	//Z other IO functions
	int get_sowingDay() { return SowingDay; }

	//Z reserve for error catching
	int get_errStatus() { return 0; }

private:
	TInitInfo initInfo;
	TGasExSpeciesParam  GasExParam;
	unique_ptr<WheatPlant> wheatPlant;
	unique_ptr<Timer> time;
	TWeather* weather;

	//Z INPUT, name of the crop variety file
	char varietyFile[133];
	//Z OUTPUT, name of the crop output file as g01
	char cropFile[133];
	//Z OUTPUT, name of the debug file
	char DebugFile[133];
	//Z OUTPUT, name of the leaf output file
	char LeafFile[133];
	char outputFile[133];
	char logFile[133];

	//Z current record number, used to index weather files
	int iCur;

	//Z simulation time recorder
	int firstDayOfSim;
	int lastDayOfSim;
	int SowingDay;

	//Z error flag, maybe not that 
	int errorFlag;

	//Z actual supply of water to plant, mol (water) m-2 (leaf) s-1
	float ET_supply;

	float RootWeightFrom2DSOIL;
	float MaxRootDepth;
	float AvailableWater;

};
#endif