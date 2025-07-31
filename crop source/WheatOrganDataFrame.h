#pragma once
#ifndef _WHEAT_OrganDataFrame_H_
#define _WHEAT_OrganDataFrame_H_

#include "WheatLeaf.h"
#include "WheatInternode.h"
#include "WheatChaff.h"
#include "WheatRachis.h"
#include "WheatSpikelet.h"

#include <memory>
#include <vector>

struct WheatOrganDataFrame
{
	//***** vectors for GPU facility *************
//Z d means device, h means host, g means global
//Z data frame for all the leaves
	vector<WheatLeaf> h_gLeaf;
	vector<WheatInternode> h_gIntr;
	vector<WheatChaff> h_gChaf;
	vector<WheatRachis> h_gRchs;
	vector<WheatSpikelet> h_gSpik;

	//Z vector for all the living fractors
	vector<float> h_gLeafLivingFrac;
	vector<float> h_gIntrLivingFrac;
	vector<float> h_gChafLivingFrac;
	vector<float> h_gRchsLivingFrac;
	vector<float> h_gSpikLivingFrac;

	//Z vector for all the force_to_death factors
	vector<bool> h_gLeafForce2Death;
	vector<bool> h_gIntrForce2Death;
	vector<bool> h_gChafForce2Death;
	vector<bool> h_gRchsForce2Death;
	vector<bool> h_gSpikForce2Death;

	//Z vector for emerge/mature of previous leaves
	//  and emergence distance
	vector<bool> h_gLeafPreLemerge;
	vector<bool> h_gIntrPreLigulation;
	vector<float> h_gLeafLigulationDist;
	vector<float> h_gIntrLigulationDist;


	//Z other growth controller vectors
	vector<bool> h_gIntrElongation;	//Z internode elongation
	vector<bool> h_gSpikFlowerInit;	//Z spikelet flower initiaiton, abortion and fertilizaiton
	vector<bool> h_gSpikFlowerAbort;
	vector<bool> h_gSpikFlowerFert;

	//Z organ growth operator
	WheatLeafUpdate wLeafUpdate;
	WheatInternodeUpdate wIntrUpdate;
	WheatChaffUpdate wChafUpdate;
	WheatRachisUpdate wRchsUpdate;
	WheatSpikeletUpdate wSpikUpdate;

	//Z organ mass assignment operator
	WheatLeafMassAssignment wLeafMassAssign;
	WheatInternodeMassAssignment wIntrMassAssign;
	WheatChaffMassAssignment wChafMassAssign;
	WheatRachisMassAssignment wRchsMassAssign;
	WheatSpikeletMassAssignment wSpikMassAssign;

	//Z record the plant organ numbers 
	int h_gLeafNum;
	int h_gIntrNum;
	int h_gChafNum;
	int h_gRchsNum;
	int h_gSpikNum;

	vector<float> dvec_LfNum;
	vector<float> dvec_greenLfNum;
	vector<float> dvec_LfArea;
	vector<float> dvec_greenLfArea;
	vector<float> dvec_seneLfArea;
	vector<float> dvec_dropLfArea;
	vector<float> dvec_LfMass;
	vector<float> dvec_ShMass;
	vector<float> dvec_dropLfMass;
	vector<float> dvec_dropShMass;
	vector<float> dvec_LfNitrogenMass;
	vector<float> dvec_ShNitrogenMass;
	vector<float> dvec_LfNitrogenReturn;
	vector<float> dvec_ShNitrogenReturn;
	vector<float> dvec_deadLfNitrogenMass;
	vector<float> dvec_deadShNitrogenMass;
	vector<float> dvec_seneLfNitrogenMass;
	vector<float> dvec_seneShNitrogenMass;
	vector<float> dvec_dropLfNitrogenMass;
	vector<float> dvec_dropShNitrogenMass;
	vector<float> dvec_ptnLfMassIncrease;
	vector<float> dvec_ptnLfNitrogenMassIncrease;

	vector<float> dvec_InLength;
	vector<float> dvec_InMass;
	vector<float> dvec_InNitrogenMass;
	vector<float> dvec_InNitrogenReturn;
	vector<float> dvec_deadInNitrogenMass;
	vector<float> dvec_ptnInMassIncrease;
	vector<float> dvec_ptnInNitrogenMassIncrease;

	vector<float> dvec_CfMass;
	vector<float> dvec_deadCfMass;
	vector<float> dvec_CfNitrogenMass;
	vector<float> dvec_CfNitrogenReturn;
	vector<float> dvec_deadCfNitrogenMass;
	vector<float> dvec_ptnCfMassIncrease;
	vector<float> dvec_ptnCfNitrogenMassIncrease;

	vector<float> dvec_KrNum;
	vector<float> dvec_KrMassSum;
	vector<float> dvec_KrNitrogenMassSum;
	vector<float> dvec_ptnKrMassIncreaseSum;
	vector<float> dvec_ptnKrNitrogenMassIncreaseSum;

	vector<float> dvec_RaLength;
	vector<float> dvec_RaMass;
	vector<float> dvec_deadRaMass;
	vector<float> dvec_RaNitrogenMass;
	vector<float> dvec_RaNitrogenReturn;
	vector<float> dvec_deadRaNitrogenMass;
	vector<float> dvec_ptnRaMassIncrease;
	vector<float> dvec_ptnRaNitrogenMassIncrease;

	//Z ----------- WHEAT LEAF CONSTRUCTOR ---------------------------------------
	WheatOrganDataFrame();

};

#endif