#include "stdafx.h"
#include "WheatController.h"
#include "initinfo.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#ifndef FLOAT_EQ
#define EPSILON 0.001f   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif
#define comma ","
#define MINUTESPERDAY 1440.0f

// const a = 17.27; b = 237.7; //constant in deg C
// inline static float E_sat(float T) { return 0.6105f * exp(17.27f * T / (237.7f + T)); }
//#using <mscorlib.dll>
using namespace std;


WheatController::WheatController(const char* filename, const char* outfile, const char* LFile, TInitInfo iniInfo)
{
	//Z we use smart pointer, no need to initialize
	//time = NULL;
	//wheatPlant = NULL;

	//Z but we do not own "weather" instance, that is owned by the "Crop" class
	weather = NULL;

	char* pch = (char*)calloc(133, sizeof(char));
	char* next = (char*)calloc(133, sizeof(char));
	char* ext = (char*)".dbg";
	char* temp = (char*)calloc(133, sizeof(char));

	//Z copy varitey input filename to the "varietyFile" buffer
	//  strcpy_s is the thread-safe version of strcpy
	strcpy_s(varietyFile, filename);

	//Z copy g01 output filename to the "cropFile" buffer
	strcpy_s(cropFile, outfile);

	//Z copy g01 output filename to the "temp" buffer
	strcpy_s(temp, 133, cropFile);

	//Z split "temp" using deliminator "." and return the first token, the rest tokens will be assigned "next" as 2D char array, so it is a pointer to pointer to char
	//  strtok_s is the thread-safe version of strtok
	pch = strtok_s(temp, ".", &next);

	//Z concatenate two strings, will produce "outputfilename.dbg"
	//  TODO, update to "strcat_s()"
	temp = strcat(pch, ext);

	//Z copy debug output filename to the "DebugFile" buffer
	strcpy_s(DebugFile, temp);

	//Z copy leaf output file to the "LeafFile" buffer
	strcpy_s(LeafFile, LFile);

	//Z a local version of the initial information file
	initInfo = iniInfo;

	//Z iterator for the weather
	iCur = 0;

	//Z simulation time recorder
	firstDayOfSim = 0;
	lastDayOfSim = 365;

	//Z initialize and read files
	this->initialize();

	//Z initialize the error flag
	errorFlag = 0;

}

WheatController::~WheatController()
{
	if (weather != NULL)
	{
		delete[] weather;
		weather = NULL;
	}
}

//Z initialize the simulation
//  mainly initialized the parameters and read fromt the input files
void WheatController::initialize()
{
	//Z define a timer to convert dates
	//  present days based on mm/dd/yy
	Timer dConvert;
	int mm, dd, yy;
	cout << "Initializing Controller object...." << endl << endl << endl;
	//Z setiosflags is an io function to control the alignment
	cout << setiosflags(ios::left) << endl
		<< " ***********************************************************" << endl
		<< " *      WheatSim: C++ wheat growth simulator               *" << endl
		<< " *                 VERSION  0.2 2023                       *" << endl
		<< " *   USDA-ARS, Adaptive Cropping Sysems Laboratory         *" << endl
		<< " ***********************************************************" << endl
		<< endl << endl;

	//Z read parameters first
	//  read variety file here and fill in data structures
		//Z to hold duumy strings from variety file
		//  unimportant thing goes to BUFFER
	char* Buffer = (char*)calloc(256, sizeof(char));
	try
	{
		ifstream cfs(varietyFile, ios::in);
		if (!cfs)
		{
			throw "Variety File not found.";
		}
		cfs.getline(initInfo.description, sizeof(initInfo.description), '\n');
		string cult = initInfo.description;
		//Z bypass two description lines, 
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> initInfo.gdd_base >> initInfo.gdd_opt >> initInfo.gdd_max
			>> initInfo.Temp_Tref >> initInfo.Temp_Ttransition
			>> initInfo.Temp_Ea_R >> initInfo.Temp_DS_R >> initInfo.Temp_DH_R;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> initInfo.fpibpa >> initInfo.srpa >> initInfo.drpa >> initInfo.tspa
			>> initInfo.iepa >> initInfo.jtpa >> initInfo.bootpa >> initInfo.headpa
			>> initInfo.antspa >> initInfo.antepa >> initInfo.matpa;;
		initInfo.GDD_rating = 1900;
		// end reading cultivar specific data from variety file
		// now read species specific data at the end of the file
		// loop until we find the location in the file
		string location = "[Gas_Exchange Species Parameters]";
		int result = -1;
		char strTest[255];
		string strTest2;
		do
		{
			//Z loop until finding the "[Gas_Exchange Species Parameters]"
			cfs.getline(Buffer, 256, '\n');
			strcpy_s(strTest, Buffer);
			strTest2.assign(strTest);
			result = strTest2.find(location);

		} while (result == -1);

		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.EaVp >> GasExParam.EaVc >> GasExParam.Eaj >> GasExParam.Hj >>
			GasExParam.Sj >> GasExParam.TPU25 >> GasExParam.Vcm25 >> GasExParam.Jm25 >>
			GasExParam.Rd25 >> GasExParam.Ear >> GasExParam.g0 >> GasExParam.g1;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.f >> GasExParam.scatt >> GasExParam.Kc25 >> GasExParam.Ko25 >>
			GasExParam.Eac >> GasExParam.Eao;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.sf >> GasExParam.phyf >> GasExParam.stomaRatio;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.internalCO2Ratio >> GasExParam.SC_param >> GasExParam.BLC_param;

		//Z finish reading, close the variety file
		if (cfs.eof()) cfs.close();
	}
	catch (const char* message)
	{
		cerr << message << "\n";
		exit(1);
	}

	firstDayOfSim = initInfo.beginDay;
	lastDayOfSim = initInfo.endDay;
	SowingDay = initInfo.sowingDay;

	//Z Timer class gets stepsize in hours
	dConvert.caldat(firstDayOfSim, mm, dd, yy);
	time = make_unique<Timer>(dd, mm, yy, initInfo.timeStep / 60.0f);

	//Z Weather class for the whole
	//  counting total records of weather data
	int dim = static_cast<int>(static_cast<float>((lastDayOfSim + 1) - firstDayOfSim) * (MINUTESPERDAY / initInfo.timeStep));
	weather = new TWeather[dim];

	//Z initialize the wheat plant here
	wheatPlant = make_unique<WheatPlant>(initInfo, GasExParam);

	//Z set heading for output (see function "output to crop file")
	//  TODO may need remove the "setiosflags(ios::left)"
	{
		ofstream cropOut(cropFile, ios::out);
		cropOut << setiosflags(ios::left)
			<< setiosflags(ios::fixed)
			<< setw(12) << "date,"
			<< setw(12) << "jday,"
			<< setw(10) << "time,"
			<< setw(12) << "LeaveNum,"			// number per representative plant
			<< setw(12) << "LeaveArea,"			// cm2 per representative plant
			<< setw(12) << "GreenLfNum,"		// number per representative plant
			<< setw(12) << "GreenLfArea,"		// cm2 per representative plant
			<< setw(12) << "TillerNum,"			// number per representative plant
			<< setw(13) << "ShootMass,"			// grams per representative plant
			<< setw(14) << "RootMass,"			// grams per representative plant
			<< setw(13) << "NitroMass,"			// mg per representative plant
			<< setw(13) << "GrossPhotosyn,"		// grams per representative plant
			<< setw(13) << "NetPhotosyn,"		// grams per representative plant
			<< setw(15) << "CumuTTd(plant),"
			<< setw(15) << "CumuTTd(climate),"
			<< setw(15) << "PltLivingFrac,"
			<< setw(15) << "BiomassReserve,"
			<< setw(15) << "MaintRespiration,"
			<< setw(15) << "NitrogenAssignment,"
			<< setw(15) << "LeafNRelease,"
			<< setw(15) << "SheathNRelease,"
			<< setw(15) << "InternodeNRelease,"
			<< setw(15) << "HourlyNitrUp,"
			<< setw(15) << "NPool,"
			<< setw(15) << "Sunlit_PFD,"
			<< setw(15) << "Shaded_PFD,"
			<< setw(15) << "Sunlit_A_gross,"
			<< setw(15) << "Shaded_A_gross,"
			<< setw(15) << "RUE"
			<< endl;
	}
}

//Z read files from 2dsoil module
//  weather data format
void WheatController::readWeatherFrom2DSOIL(const TWeather& wthr)
{
	weather[iCur] = wthr;
	//Z actual supply of water to plant, mol (water) m-2 (leaf) s-1
	ET_supply = wthr.ET_supply;
	weather[iCur].daytime = wthr.jday + wthr.time;
}

//Z main engine of the rye model
//  weather is from soil, and gas exchange parameters are read in this controller
//
int WheatController::run(const TWeather& wthr)
{
	//Z first read the weather file
	readWeatherFrom2DSOIL(wthr);

	//Z use the weather data to call the plant variable once,
	//  then process the rye growth computation.
	if (weather[iCur].jday >= initInfo.sowingDay && weather[iCur].jday <= lastDayOfSim)
	{
		RootWeightFrom2DSOIL = wthr.TotalRootWeight;
		wheatPlant->WheatPlantUpdate(weather[iCur]);
		MaxRootDepth = wthr.MaxRootDepth;
		AvailableWater = wthr.ThetaAvail;
		outputToCropFile();

		iCur++;
		time->step();
	}
	return 0;
}

//Z output to the crop results g01
void WheatController::outputToCropFile()
{
	int mm, id, iyyy;
	string DateForOutput;
	float av_gs = wheatPlant->get_conductance();
	if (!wheatPlant->get_develop()->is_emerge())
	{
		av_gs = 0;
	}
	float vpd = wheatPlant->get_VPD();
	if (vpd < 0)
	{
		vpd = 0;
	}
	time->caldat(weather[iCur].jday, mm, id, iyyy);
#if 0
	DateForOutput.Format("%.2d/%.2d/%4i", mm, id, iyyy);
#else
	char DateForOutputBuff[16];
	sprintf(DateForOutputBuff, "%.2d/%.2d/%4i", mm, id, iyyy);
	DateForOutput = DateForOutputBuff;
#endif
	ofstream ostr(cropFile, ios::app);
	ostr << setiosflags(ios::right)
		<< setiosflags(ios::fixed)
		<< setw(9) << DateForOutput << comma
		<< setw(6) << weather[iCur].jday << comma
		<< setw(8) << setprecision(0) << weather[iCur].time * 24.0 << comma
		<< setw(11) << setprecision(2) << wheatPlant->get_LfNum() << comma
		<< setw(11) << setprecision(2) << wheatPlant->get_LaArea() << comma
		<< setw(11) << setprecision(2) << wheatPlant->get_greenLfNum() << comma
		<< setw(12) << setprecision(2) << wheatPlant->get_greenLaArea() << comma
		<< setw(12) << setprecision(2) << wheatPlant->get_TlNum() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_shootMass() << comma
		<< setw(13) << setprecision(2) << RootWeightFrom2DSOIL << comma //Z it is interesting that root weight per plant from 2DSOIL is set as a component of weather
		<< setw(13) << setprecision(6) << wheatPlant->get_PltTotalNitrogen() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_grossPhotosynthesis() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_netPhotosynthesis() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_develop()->get_WheatTTd().get_TTdSum_plt_dssat() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_develop()->get_WheatTTd().get_TTdSum_env_dssat() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_develop()->get_plantLivingFraction() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_biomassReserve() << comma
		<< setw(13) << setprecision(2) << wheatPlant->get_MaintenanceRespiration() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_PltShootNitrogenAssignment() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_LaNitrogenReturn() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_ShNitrogenReturn() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_InNitrogenReturn() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_HourlyNitrogenSoilUptake() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_nitrogenPool() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_sunlit_PFD() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_shaded_PFD() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_sunlit_A_gross() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_shaded_A_gross() << comma
		<< setw(13) << setprecision(6) << wheatPlant->get_averagedBiomassLeftover()
		<< endl;
	ostr.close();
}