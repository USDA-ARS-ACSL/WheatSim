#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY 1440.0f
#ifndef FLOAT_EQ
#define EPSILON 0.001f   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define MINUTESPERDAY 1440.0f
#endif

struct TInitInfo
{
public:
	TInitInfo()
	{
		GDD_rating = 1331;
		genericLeafNo = 15;
		latitude = 38.0; longitude = 0.0; altitude = 50.0;
		sowingDay = 150;
		beginDay = 1; endDay = 365;
		year = 2004;
		timeStep = 5.0f;
		plantDensity = 8.0f;
		CO2 = 370.0f;
		Rmax_LIR = 0.0978f;
		Rmax_LTAR = 0.53f;
		DayLengthSensitive = true;
		PhyllochronsToSilk = 8.0f;
		PhyllochronsToTassel = 1.0f;
		stayGreen = 4.5f;
		LM_min = 125.0f;

		//Z make the following variables in this structure for ryesim
		fpibpa = srpa = drpa 
			= tspa = iepa = jtpa 
			= bootpa = headpa = antspa 
			= antepa = matpa = 0.0f;
		seedDepth = 0.0f;

		//Z two ways to compute thermal time
		//first the dssat way
		gdd_base = 0.0f;
		gdd_opt = 26.0f;
		gdd_max = 34.0f;
		//second the fpsm wheat way using BETA functions
		Temp_Tref = 12.0f;
		Temp_Ttransition = 9.0f;
		Temp_Ea_R = 8900.0f;
		Temp_DS_R = 68.432f;
		Temp_DH_R = 20735.5f;

	};
	char description[255] = "\0";
	char cultivar[255] = "\0";
	int GDD_rating; // GDD or GTI rating of the cv, see Stewart 1999 for conversion between MRMR and other ratings
	int genericLeafNo; // leaf number at the end of juvenile phase independent of environmental ques of leaf initiation
	float plantDensity;
	float latitude, longitude, altitude;
	int sowingDay, beginDay, endDay;
	float CO2;
	int year;
	float timeStep;
	bool DayLengthSensitive; //1 if daylength sensitive
	float Rmax_LIR, Rmax_LTAR; //  Maximum Leaf tip initiation and appearance rates
	float stayGreen;  // staygreen trait of the hybrid (originally 4.5)
	float LM_min; //Length of largest leaf
	// stay green for this value times growth period after peaking before senescence begins
	// An analogy for this is that with no other stresses involved, it takes 15 years to grow up, stays active for 60 years, and age the last 15 year if it were for a 90 year life span creature.
	//Once fully grown, the clock works differently so that the hotter it is quicker it ages
	float PhyllochronsToSilk; //number of phyllochrons from tassel initiation for 75% silking.
	float PhyllochronsToTassel; // number of phyllochrons past tassel initiation when tassels are fully emerged. (not input yet)
	//todo these above 2 variables are also in development - need to remove them from there.
	//check units

	// ***************************************************
	//Z make the following variables in this structure for WHEATSIM

	//Z phyllochrons numbers for cereal phenology, the last one will be based on gdd
	float fpibpa;	//number of phyllochrons after double ridge that flower primordium initiation begins.
	float srpa;	//number of phyllochrons from vernalization to single ridge.
	float drpa;	//number of phyllochrons between singleridge and double ridge.
	float tspa;	//number of phyllochrons between double ridge and terminal spikelet growth stages.
	float iepa;	//number of phyllochrons between double ridge and start of internode elongation growth stages.
	float jtpa;	//number of phyllochrons between beginning of internode elongation and jointing.
	float bootpa;	//number of phyllochrons between the beginning of jointing and booting (i.e., flag leaf fully grown).
	float headpa;	//number of phyllochrons between the beginning of booting and heading.
	float antspa;	//number of phyllochrons between the beginning of heading and beginning of anthesis.
	float antepa;	//number of phyllochrons between the beginning of anthesis and end of anthesis (duration of anthesis).
	float matpa;	//number of growing degree-days between the beginning of anthesis and physiological maturity.

	float seedDepth;  //seeding depth

	//Z two ways to compute thermal time
	//first the dssat way
	float gdd_base; //base tmpr for computing gdd
	float gdd_opt;  //optimal tmpr for computing gdd
	float gdd_max;  //max tmpr for computing gdd, beyond that, no gdd and plant stop morph growth
	//second the fpsm wheat way using BETA functions
	//i.e., the temerpature response function
	float Temp_Tref;		//Z =12 Arbitrary reference temperature (°C)
	float Temp_Ttransition;	//Z =9 Below this temperature f = linear function of temperature instead of Arrhenius - like(°C)
	// parameters for the Arrhenius and modified Arrhenius models
	float Temp_Ea_R;		//Z =8900 Parameter Ea / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (K)
	float Temp_DS_R;		//Z =68.432 Parameter deltaS / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (dimensionless)
	float Temp_DH_R;		//Z =20735.5 Parameter deltaH / R in Eyring equation from Johnson and Lewin(1946) - Parameter value fitted from Kemp and Blacklow(1982) (K)
};
#endif