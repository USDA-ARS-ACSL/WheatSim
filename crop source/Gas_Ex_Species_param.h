#pragma once
#ifndef _GAS_EX_SPECIES_PARAM_H_
#define _GAS_EX_SPECIES_PARAM_H_

/*!    \struct tparms gas_exchange.h
			   \brief Structure to hold plant parameters for the photosynthesis model. \n
			   Some parameters are specific for C3 or C4 type Plants \n
		*/

		/*!
		//                parameter        description
		@param 		Vcm25	4		Photosynthetic Rubisco Capacity at 25C (umol m-2 s-1)
		@param 		Jm25	5		Potential Rate of electron transport at 25C  (umol m-2 s-1)
		@param 		Vpm25	6		C4 Carboxylation rate at 25C (C4, umol m-2 s-1)
		@param 		TPU25	7		Rate if Triose Phosphate Utilization at 25C (C3, umol m-2 s-1)
		@param 		RD25	8		Mitochondrial respiration in the light at 25C (umol m-2 s-1)
		@param 		Theta   9		Initial slope of CO2 response (umol m2 s-1) - de Pury (1997)
		@param 		EaVc   10		Activation energy for Arrhenius function used to calculate temperature dependence for Vcmax (kJ mol-1)
		@param 		Eaj    11		Activation energy for Arrhenius function used to calculate temperature dependence for J (kJ mol-1)
		@param 		Hj     12		Curvature parameter of the temperature dpendence of Jmax (kJ mol-1)
		@param		Sj	   13		Electron transport temperature response parameter for Jmax (J mole-1 K-1)
		@param      Hv     14       Curvature parameter of the temperature dependence of Vcmax (J mole-1)
		@param 	    EaVp   15		Activation energy for Arrhenius function used to calculate temperature dependence for Vpmax (kJ mol-1)
		@param 	    Sv	   16		Electron transport temperature response parameter for Vcmax (J mole-1 K-1)
		@param 		EAP	   17		Activation energy for Arrhenius function used to calculate temperature dependence for TPU (kJ mol-1)
		@param 		EAR	   18	 	Activation energy for Arrhenius function used to calculate temperature dependence for respiration (kJ mol-1)
		@param 		g0	   19		Minimum stomatal conductance to water vapor at the light compensation point in the BWB model (mol m-2 s-1)
		@param 	    g1	   20   	Empirical coefficient for the sensitivity of Stomatal Conductance to A, Cs and hs in BWB model (no units)
		@param 		StomRatio	21	Stomatal Ratio (fraction)
		@param  	LfWidth     22	Leaf Width (m)
		@param 		LfAngFact	23	Leaf Angle Factor
		@param 		Remark		24	Text
		*/

struct TGasExSpeciesParam
{
public:
	/*********************/
	/* Assigns parameters to the 3 components of Farquar C3 model
	/* Includes T-dependencies
	/*********************/

	TGasExSpeciesParam()
	{
		/*Z Modified Arrhenius function Parameters for
			maximum carboxylation rate of Rubisco, or reaction speed Vcmax "Photosynthetic Rubisco Capacity"
			refer to Liang Fang, 2022 Neglecting acclimation of photosynthesis under drought can cause significant errors in predicting leaf photosynthesis in wheat
			Table S5 well watered treatment
			*/

		//Z CO2-saturated maximum carboxylation rate of Rubisco, umol m-2 s-1
		//  also see Zhao 2020 "Quantifying key model parameters for wheat leaf gas exchange under different environmental conditions"
		Vcm25 = 158.82f;
		//Z Activation energy for temperature dependence, J mol-1 Ex
		EaVc = 76701.0f;
		//Z Deactivation energy, J mol-1 Dx
		Hv = 203791.0f;
		//Z Entropy factor, J K-1 mol-1 Sx
		Sv = 679.0f;

		/*Z Modified Arrhenius function Parameters for
			maximum potential electron transport rate under saturating irradiance, or J "Potential Rate of electron transport"
			refer to Liang Fang, 2022 Neglecting acclimation of photosynthesis under drought can cause significant errors in predicting leaf photosynthesis in wheat
			Table S5 well watered treatment
		*/

		//Z Light - saturated potential electron transport rate, umol m-2 s-1
		//  also see Zhao 2020 "Quantifying key model parameters for wheat leaf gas exchange under different environmental conditions"
		Jm25 = 205.26f;
		//Z Activation energy for temperature dependence, J mol-1
		Eaj = 51049.0f;
		//cccz: Deactivation energy, J mol-1
		Hj = 205069.0f;
		//cccz: Entropy factor, J K-1 mol-1
		Sj = 680.0f;

		//Z  Arrhenius function used to calculate temperature dependence for TPU
		/*Z  Rate of Triose Phosphate Utilization, TPU, at 25C(C3, umol m - 2 s - 1)
			 refer to Liang Fang, 2022 Neglecting acclimation of photosynthesis under drought can cause significant errors in predicting leaf photosynthesis in wheat
			 Figure 5 ef
		*/
		TPU25 = 9.09f;
		//Z Activation energy (J mol - 1)
		/*  refer to a barley value at Braune, 2009 Integrating eects of leaf nitrogen, age,
			rank, and growth temperature into the photosynthesis-stomatal conductance model LEAFC3-N
			parameterised for barley
		*/
		EaVp = 47100.0f;


		/*Z Mitochondrial respiration in the light at 25C umol m-2 s-1
			refer to Liang Fang, 2022 Neglecting acclimation of photosynthesis under drought can cause significant errors in predicting leaf photosynthesis in wheat, Figure S5 ed
			Bernacchi, 2009 Modeling the Temperature Dependence of C3 Photosynthesis Fig 10.5
			Yu, 2004 Simulation of the Stomatal Conductance of Winter Wheat in Response to Light, Temperature and CO2 Changes Fig. 1
		*/
		Rd25 = 1.0f;
		//Z Activation energy (J mol - 1)
		/*  exponential rate of arrhenious function for mitochondrial respiration(J mol)
		*/
		Ear = 66400.0f;

		//Z Min stomatal conductance to water vapor at the light compensation point in the BWB model (mol m-2 s-1)
		//  refer to Zhao 2020 Quantifying key model parameters for wheat leaf gas exchange under different environmental conditions
		g0 = 0.023f;

		//Z Empirical coefficient for the sensitivity of Stomatal Conductance to A, Cs and hs in BWB model (coefficient)
		//  refer to Zhao 2020 Quantifying key model parameters for wheat leaf gas exchange under different environmental conditions
		g1 = 12.172f;

		//Z seems not used or defined in the current code, make it 1.0 to facilitate the computation.
		g2 = 1.0f;

		//Z parameter determine the relation between intercellular CO2 partial pressure c_i and CO2 partial pressure at a leaf surface C_s
		//  relative to the air CO2 pressure Ca
		//  partial pressure in unit ubar

		/*Z C_i = C_a - A[1.6 / gs + 1.37 / gb]Pa
			gs - critical values are the vapor conductance through the stomata
			gb - vapor conductance near the leaf surface
			SC_param=1.57, SLC_param=1.36

			gs, gb are for vapor need to convert to CO2
			conversion factors such as 1.6 and 1.37 are the approximation to ratios between the vapor conductance and CO2 conductance
			refer to Table 7.6 in “An Introduction to Environmental Biophysics”, Campbell and Norman, 1998
		*/
		SC_param = 1.57f;
		BLC_param = 1.36f;

		//Z parameters for set_PSIleafeffect()
		//  Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
		//  sf - sensitive factor
		//  phyf - Reference_Potential_(phyf, bars)
		sf = 4.2f;
		phyf = -0.5f;

		//Z parameter in boundary layer conductance to vapor, function gbw()
		stomaRatio = 0.8f;

		//Z holds factor to multiply CO2 by to get first estimate of Ci
		//  a default value in the GasEx_psil
		internalCO2Ratio = 0.7f;

		//Z in the photosynthesis function
		//  scatt - change irradiance to absorbed irradiance
		//  f     - change irradiance to radiation goes to PSII
		scatt = 0.15f;	//Z leaf reflectance + transmittance
		f = 0.15f;		//Z spectral correction

		//Z RuBisCO speed constants
		//Z The co2 compensation point which is also used for defining gs is derived from carboxylation (vc) as well as oxygenation (vo) velocities
		//  using an empirical estimate of O2 concentrations in the plant tissue, 
		//  kinetic values that are derived from temperature-corrected species-specific parameteres (KO25, KC25) and enzyme activities (AEKC, AEKO).

		//Z kc25, ko25: temperature corrected Michaelis-Menten parameters for carboxylation and oxygenation @25
		//  Eac, Eao:   activation energy for temperature relations
		Kc25 = 291.0f;
		Ko25 = 194.0f;
		Eac = 80990.0f;
		Eao = 23720.0f;
	}
public:
	//Z Modified Arrhenius function Parameters for max carboxylation rate of Rubisco, or reaction speed Vcmax "Photosynthetic Rubisco Capacity"
	float Vcm25, EaVc, Hv, Sv;

	//Z Modified Arrhenius function Parameters for max potential electron transport rate under saturating irradiance, or J "Potential Rate of electron transport"
	float Jm25, Eaj, Hj, Sj;

	//Z Arrhenius function used to calculate temperature dependence for TPU
	float TPU25, EaVp;

	//Z Mitochondrial respiration in the light at 25C umol m-2 s-1
	float Rd25, Ear;

	//Z Min stomatal conductance and Empirical coefficient for the sensitivity of Stomatal Conductance to water vapor
	float g0, g1, g2;

	//Z conductance conversion factor for CO2 and vapor transfer through stomata
	float SC_param, BLC_param;

	//Z parameters for set_PSIleafeffect(), sensitivity factor and reference pressure
	float sf, phyf;

	//Z parameter in boundary layer conductance to vapor, function gbw()
	float stomaRatio;

	//Z a default value in the GasEx_psil, 0.7 air CO2 will be intercellular CO2
	float internalCO2Ratio;

	//Z in the photosynthesis function
	//  scatt - change irradiance to absorbed irradiance
	//  f     - change irradiance to radiation goes to PSII
	float scatt;
	float f;

	//Z RuBisCO speed constants
	float Kc25, Ko25;
	float Eac, Eao;

};
#endif