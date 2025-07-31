#pragma once
#ifndef _WEATHER_H_
#define _WEATHER_H_
/* I use this structure to hold variables transfered to and from 2DSOIL */
struct TWeather
{
public:
	TWeather()
	{
		year = 2004;
		jday = 1;
		time = 0.0f;
		daytime = static_cast<float>(jday) + time;
		CO2 = 370.0f;
		airT = 20.0f, soilT = airT;
		PFD = 0.0f;
		solRad = 0.0f;
		RH = 50.0f;
		wind = 1.5f;
		rain = 0.0f;
		dayLength = 12.0f;
		SoilMP_med = -100.0f;
		LeafWP = -0.5f;
		PredawnLWP = -0.05f;

		DailyOutput = HourlyOutput = 0;

		//Z initializaiton parameters
		ET_supply = 0.0f;
		MaxRootDepth = 0.0f;
		TotalRootWeight = 0.0f;
		ThetaAvail = 0.0f;
		pcrl = 0.0f;
		pcrq = 0.0f;
		pcrs = 0.0f;


		//Z additional variables
		airTDmax = 0.0f;
		airTDmin = 0.0f;
		wfps = 0.0f;
		doy = 0;

	}
	int year;
	int jday;

	float time, daytime;
	//dt added pcrl and pcrq to check balances
	float CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT, ET_supply,
		LeafWP, PredawnLWP,
		pcrs, pcrl,
		pcrq, TotalRootWeight, SoilMP_med;
	float MaxRootDepth, ThetaAvail;
	//double CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT;
	// DT added output flags from 2DSOIL
	int DailyOutput, HourlyOutput;

	//************ additional variable for ryesim *************
	//Z daily max and min temperature for Vernalisation
	float airTDmax, airTDmin;
	//Z soil water satration degree, namely water content / ThFull (max availiable water)
	float wfps;
	//Z day of year
	int doy;
};
#endif