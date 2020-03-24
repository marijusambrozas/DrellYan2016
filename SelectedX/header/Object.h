///////////////////////////////////////////////////////////////
// -- 2015.03.21: Replace MuonID/ElectronID with NtupleHandle
// -- 2015.04.02: Fix some minor bugs(missing penenthesis)
// -- 2015.04.22: Correct Muon TightID
// -- 2016.01.02: Change the class name as "Objects" and include GenOthers, photon, jet and MET
// -- 2017.04.27: Adding "Electron_passMediumID" (by Dalmin Pai)
// -- 2017.07.19: Updating Electron variables with EGM corrections (by Dalmin Pai)
// -- 2017.07.19: Finishing "isMediumElectron_2016dataFor80X()" (by Dalmin Pai)
// -- 2017.07.26: Varialbe type of "Electron_passConvVeto" was changed to "Bool_t" (by Dalmin Pai)
// -- 2017.10.31: Changing Electron variables (by Dalmin Pai)
// -- 2017.12.01: Updating GenOthers with quark information
// -- 2018.01.24: Add "RelPFIso_dBeta" in muon
///////////////////////////////////////////////////////////////

#pragma once

#include <TLorentzVector.h>

//customized header files
#include "NtupleHandle.h"
#include "SelectedX.h"

#define M_Mu 0.1056583715 // -- GeV -- //
#define M_Elec 0.000510998 // -- GeV -- //
#define M_Tau 1.77682 // -- GeV -- //
#define M_Top 172.5 // -- GeV -- //

class Object
{
public:
	Double_t Pt;
	Double_t Et;
	Double_t eta;
	Double_t phi;
	TLorentzVector Momentum;
};

class GenLepton : public Object
{
public:
	Double_t Px;
	Double_t Py;
	Double_t Pz;
	Double_t Mother;
	Int_t ID;
	Double_t charge;
	Double_t Mass;
	Int_t Status;

	//GenFlags(after 7_4_X)
	Int_t isPrompt;
	Int_t isPromptFinalState;
	Int_t isTauDecayProduct;
	Int_t isPromptTauDecayProduct;
	Int_t isDirectPromptTauDecayProductFinalState;
	Int_t isHardProcess;
	Int_t isLastCopy;
	Int_t isLastCopyBeforeFSR;
	Int_t isPromptDecayed;
	Int_t isDecayedLeptonHadron;
	Int_t fromHardProcessBeforeFSR;
	Int_t fromHardProcessDecayed;
	Int_t fromHardProcessFinalState;

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Pt 		= ntuple->GenLepton_pT[index];
		eta 	= ntuple->GenLepton_eta[index];
		phi 	= ntuple->GenLepton_phi[index];
		Px 		= ntuple->GenLepton_px[index];
		Py 		= ntuple->GenLepton_py[index];
		Pz 		= ntuple->GenLepton_pz[index];
		Mother 	= ntuple->GenLepton_mother[index];
		ID 		= ntuple->GenLepton_ID[index];
		charge 	= ntuple->GenLepton_charge[index];
		Status 	= ntuple->GenLepton_status[index];
		
                if(ID == -11 || ID == 11)
		    Mass = M_Elec;
                else if(ID == -13 || ID == 13)
		    Mass = M_Mu;
                else if(ID == -15 || ID == 15)
		    Mass = M_Tau;

		Double_t E = sqrt(Px*Px + Py*Py + Pz*Pz + Mass*Mass);
		Momentum.SetPxPyPzE(Px, Py, Pz, E);

		isPrompt 									= ntuple->GenLepton_isPrompt[index];
		isPromptFinalState 							= ntuple->GenLepton_isPromptFinalState[index];
		isTauDecayProduct 							= ntuple->GenLepton_isTauDecayProduct[index];
		isPromptTauDecayProduct 					= ntuple->GenLepton_isPromptTauDecayProduct[index];
		isDirectPromptTauDecayProductFinalState 	= ntuple->GenLepton_isDirectPromptTauDecayProductFinalState[index];
		isHardProcess 								= ntuple->GenLepton_isHardProcess[index];
		isLastCopy 									= ntuple->GenLepton_isLastCopy[index];
		isLastCopyBeforeFSR 						= ntuple->GenLepton_isLastCopyBeforeFSR[index];
		isPromptDecayed 							= ntuple->GenLepton_isPromptDecayed[index];
		isDecayedLeptonHadron 						= ntuple->GenLepton_isDecayedLeptonHadron[index];
		fromHardProcessBeforeFSR 					= ntuple->GenLepton_fromHardProcessBeforeFSR[index];
		fromHardProcessDecayed 						= ntuple->GenLepton_fromHardProcessDecayed[index];
		fromHardProcessFinalState 					= ntuple->GenLepton_fromHardProcessFinalState[index];
	}

	Bool_t isElectron()
	{
		Bool_t isElec = kFALSE;
                if(abs(ID) == 11)
			isElec = kTRUE;

		return isElec;
	}

	Bool_t isMuon()
	{
		Bool_t isMu = kFALSE;
                if(abs(ID) == 13)
			isMu = kTRUE;

		return isMu;
	}

	Bool_t isMotherZ()
	{
		Bool_t isZ = kFALSE;
                if(Mother == 23)
			isZ = kTRUE;

		return isZ;
	}

};

class GenOthers : public Object
{
public:
	Double_t Px;
	Double_t Py;
	Double_t Pz;
	Double_t Mother;
	Int_t ID;
	Double_t charge;
	Double_t Mass;
	Int_t Status;

	//GenFlags(after 7_4_X)
	Int_t isPrompt;
	Int_t isPromptFinalState;
	Int_t isTauDecayProduct;
	Int_t isPromptTauDecayProduct;
	Int_t isDirectPromptTauDecayProductFinalState;
	Int_t isHardProcess;
	Int_t isLastCopy;
	Int_t isLastCopyBeforeFSR;
	Int_t isPromptDecayed;
	Int_t isDecayedLeptonHadron;
	Int_t fromHardProcessBeforeFSR;
	Int_t fromHardProcessDecayed;
	Int_t fromHardProcessFinalState;

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Pt      = ntuple->GenOthers_pT[index];
		eta     = ntuple->GenOthers_eta[index];
		phi     = ntuple->GenOthers_phi[index];
		Px      = ntuple->GenOthers_px[index];
		Py      = ntuple->GenOthers_py[index];
		Pz      = ntuple->GenOthers_pz[index];
		Mother  = ntuple->GenOthers_mother[index];
		ID      = ntuple->GenOthers_ID[index];
		charge  = ntuple->GenOthers_charge[index];
		Status  = ntuple->GenOthers_status[index];
		
                if(abs(ID) == 22) // -- Photon -- //
		    Mass = 0;
                else if(ID == -6 || ID == 6)
		    Mass = M_Top;

		Double_t E = sqrt(Px*Px + Py*Py + Pz*Pz + Mass*Mass);
		Momentum.SetPxPyPzE(Px, Py, Pz, E);

		isPrompt                                    = ntuple->GenOthers_isPrompt[index];
		isPromptFinalState                          = ntuple->GenOthers_isPromptFinalState[index];
		isTauDecayProduct                           = ntuple->GenOthers_isTauDecayProduct[index];
		isPromptTauDecayProduct                     = ntuple->GenOthers_isPromptTauDecayProduct[index];
		isDirectPromptTauDecayProductFinalState     = ntuple->GenOthers_isDirectPromptTauDecayProductFinalState[index];
		isHardProcess                               = ntuple->GenOthers_isHardProcess[index];
		isLastCopy                                  = ntuple->GenOthers_isLastCopy[index];
		isLastCopyBeforeFSR                         = ntuple->GenOthers_isLastCopyBeforeFSR[index];
		isPromptDecayed                             = ntuple->GenOthers_isPromptDecayed[index];
		isDecayedLeptonHadron                       = ntuple->GenOthers_isDecayedLeptonHadron[index];
		fromHardProcessBeforeFSR                    = ntuple->GenOthers_fromHardProcessBeforeFSR[index];
		fromHardProcessDecayed                      = ntuple->GenOthers_fromHardProcessDecayed[index];
		fromHardProcessFinalState                   = ntuple->GenOthers_fromHardProcessFinalState[index];
	}

};

class Electron : public Object
{
public:
	Double_t Energy;
        Double_t Energy_uncorr;
	Int_t charge;
	Double_t gsfpT;
	Double_t gsfPx;
	Double_t gsfPy;
	Double_t gsfPz;
	Double_t gsfEta;
	Double_t gsfPhi;
	Double_t gsfCharge;
	Double_t etaSC;
	Double_t phiSC;
	Double_t etaWidth;
	Double_t phiWidth;
	Double_t dEtaIn;
	Double_t dEtaInSeed; // updated at 19 Jul. 2017 by Dalmin
	Double_t dPhiIn;
	Double_t sigmaIEtaIEta;
	Double_t Full5x5_SigmaIEtaIEta; // updated at 19 Jul. 2017 by Dalmin
	Double_t HoverE;
	Double_t fbrem;
	Double_t eOverP;
	Double_t InvEminusInvP;
	Double_t dxyVTX;
	Double_t dzVTX;
	Double_t dxy;
	Double_t dz;
	Double_t dxyBS;
	Double_t dzBS;
	Double_t chIso03;
	Double_t nhIso03;
	Double_t phIso03;
	Double_t ChIso03FromPU;
	Int_t mHits;
	Double_t EnergySC;
	Double_t preEnergySC;
	Double_t rawEnergySC;
	Double_t etSC;
	Double_t E15;
	Double_t E25;
	Double_t E55;
	Double_t RelPFIso_dBeta;
	Double_t RelPFIso_Rho;
	Double_t r9;
	Int_t ecalDriven;
//	Int_t passConvVeto;
	Bool_t passConvVeto; // modified at 26 Jul. 2017 by Dalmin
	Bool_t passMediumID; // updated at 27 Apr. 2017 by Dalmin

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Energy = ntuple->Electron_Energy[index];
                Energy_uncorr = ntuple->Electron_EnergyUnCorr[index];
		Pt = ntuple->Electron_pT[index];
		eta = ntuple->Electron_eta[index];
		phi = ntuple->Electron_phi[index];
		charge = ntuple->Electron_charge[index];
		gsfpT = ntuple->Electron_gsfpT[index];
		gsfPx = ntuple->Electron_gsfPx[index];
		gsfPy = ntuple->Electron_gsfPy[index];
		gsfPz = ntuple->Electron_gsfPz[index];
		gsfEta = ntuple->Electron_gsfEta[index];
		gsfPhi = ntuple->Electron_gsfPhi[index];
		gsfCharge = ntuple->Electron_gsfCharge[index];
		etaSC = ntuple->Electron_etaSC[index];
		phiSC = ntuple->Electron_phiSC[index];
		etaWidth = ntuple->Electron_etaWidth[index];
		phiWidth = ntuple->Electron_phiWidth[index];
		dEtaIn = ntuple->Electron_dEtaIn[index];
		dEtaInSeed = ntuple->Electron_dEtaInSeed[index]; // updated at 19 Jul. 2017 by Dalmin
		dPhiIn = ntuple->Electron_dPhiIn[index];
		sigmaIEtaIEta = ntuple->Electron_sigmaIEtaIEta[index];
		Full5x5_SigmaIEtaIEta = ntuple->Electron_Full5x5_SigmaIEtaIEta[index]; // updated at 19 Jul. 2017 by Dalmin
		HoverE = ntuple->Electron_HoverE[index];
		fbrem = ntuple->Electron_fbrem[index];
		eOverP = ntuple->Electron_eOverP[index];
		InvEminusInvP = ntuple->Electron_InvEminusInvP[index];
		dxyVTX = ntuple->Electron_dxyVTX[index];
		dzVTX = ntuple->Electron_dzVTX[index];
		dxy = ntuple->Electron_dxy[index];
		dz = ntuple->Electron_dz[index];
		dxyBS = ntuple->Electron_dxyBS[index];
		dzBS = ntuple->Electron_dzBS[index];
		chIso03 = ntuple->Electron_chIso03[index];
		nhIso03 = ntuple->Electron_nhIso03[index];
		phIso03 = ntuple->Electron_phIso03[index];
		ChIso03FromPU = ntuple->Electron_ChIso03FromPU[index];
		mHits = ntuple->Electron_mHits[index];
		EnergySC = ntuple->Electron_EnergySC[index];
		preEnergySC = ntuple->Electron_preEnergySC[index];
		rawEnergySC = ntuple->Electron_rawEnergySC[index];
		etSC = ntuple->Electron_etSC[index];
		E15 = ntuple->Electron_E15[index];
		E25 = ntuple->Electron_E25[index];
		E55 = ntuple->Electron_E55[index];
		RelPFIso_dBeta = ntuple->Electron_RelPFIso_dBeta[index];
		RelPFIso_Rho = ntuple->Electron_RelPFIso_Rho[index];
		r9 = ntuple->Electron_r9[index];
		ecalDriven = ntuple->Electron_ecalDriven[index];
		passConvVeto = ntuple->Electron_passConvVeto[index];
		passMediumID = ntuple->Electron_passMediumID[index]; // updated at 27 Apr. 2017 by Dalmin

		TLorentzVector v;
		v.SetPtEtaPhiM(Pt, eta, phi, M_Elec);
		Momentum = v;
	}

        void FillFromSelectedX(LongSelectedEE_t *ntuple, Int_t index)   // added at 2018.08.06 by Marijus Ambrozas
        {
                Energy = ntuple->Electron_Energy->at(index);
                Energy_uncorr = ntuple->Electron_Energy_uncorr->at(index);
                Pt = ntuple->Electron_pT->at(index);
                eta = ntuple->Electron_eta->at(index);
                phi = ntuple->Electron_phi->at(index);
                charge = ntuple->Electron_charge->at(index);
                gsfpT = ntuple->Electron_gsfpT->at(index);
                gsfPx = ntuple->Electron_gsfPx->at(index);
                gsfPy = ntuple->Electron_gsfPy->at(index);
                gsfPz = ntuple->Electron_gsfPz->at(index);
                gsfEta = ntuple->Electron_gsfEta->at(index);
                gsfPhi = ntuple->Electron_gsfPhi->at(index);
                gsfCharge = ntuple->Electron_gsfCharge->at(index);
                etaSC = ntuple->Electron_etaSC->at(index);
                phiSC = ntuple->Electron_phiSC->at(index);
                etaWidth = ntuple->Electron_etaWidth->at(index);
                phiWidth = ntuple->Electron_phiWidth->at(index);
                dEtaIn = ntuple->Electron_dEtaIn->at(index);
                dEtaInSeed = ntuple->Electron_dEtaInSeed->at(index);
                dPhiIn = ntuple->Electron_dPhiIn->at(index);
                sigmaIEtaIEta = ntuple->Electron_sigmaIEtaIEta->at(index);
                Full5x5_SigmaIEtaIEta = ntuple->Electron_Full5x5_SigmaIEtaIEta->at(index);
                HoverE = ntuple->Electron_HoverE->at(index);
                fbrem = ntuple->Electron_fbrem->at(index);
                eOverP = ntuple->Electron_eOverP->at(index);
                InvEminusInvP = ntuple->Electron_InvEminusInvP->at(index);
                dxyVTX = ntuple->Electron_dxyVTX->at(index);
                dzVTX = ntuple->Electron_dzVTX->at(index);
                dxy = ntuple->Electron_dxy->at(index);
                dz = ntuple->Electron_dz->at(index);
                dxyBS = ntuple->Electron_dxyBS->at(index);
                dzBS = ntuple->Electron_dzBS->at(index);
                chIso03 = ntuple->Electron_chIso03->at(index);
                nhIso03 = ntuple->Electron_nhIso03->at(index);
                phIso03 = ntuple->Electron_phIso03->at(index);
                ChIso03FromPU = ntuple->Electron_ChIso03FromPU->at(index);
                mHits = ntuple->Electron_mHits->at(index);
                EnergySC = ntuple->Electron_EnergySC->at(index);
                preEnergySC = ntuple->Electron_preEnergySC->at(index);
                rawEnergySC = ntuple->Electron_rawEnergySC->at(index);
                etSC = ntuple->Electron_etSC->at(index);
                E15 = ntuple->Electron_E15->at(index);
                E25 = ntuple->Electron_E25->at(index);
                E55 = ntuple->Electron_E55->at(index);
                RelPFIso_dBeta = ntuple->Electron_RelPFIso_dBeta->at(index);
                RelPFIso_Rho = ntuple->Electron_RelPFIso_Rho->at(index);
                r9 = ntuple->Electron_r9->at(index);
                ecalDriven = ntuple->Electron_ecalDriven->at(index);
                passConvVeto = ntuple->Electron_passConvVeto->at(index);
                passMediumID = ntuple->Electron_passMediumID->at(index);

                TLorentzVector v;
                v.SetPtEtaPhiM(Pt, eta, phi, M_Elec);
                Momentum = v;
        }

        void FillFromSelectedX(LongSelectedEMu_t *ntuple)   // added at 2018.08.07 by Marijus Ambrozas
        {
                Energy = ntuple->Electron_Energy;
                Energy_uncorr = ntuple->Electron_Energy_uncorr;
                Pt = ntuple->Electron_pT;
                eta = ntuple->Electron_eta;
                phi = ntuple->Electron_phi;
                charge = ntuple->Electron_charge;
                gsfpT = ntuple->Electron_gsfpT;
                gsfPx = ntuple->Electron_gsfPx;
                gsfPy = ntuple->Electron_gsfPy;
                gsfPz = ntuple->Electron_gsfPz;
                gsfEta = ntuple->Electron_gsfEta;
                gsfPhi = ntuple->Electron_gsfPhi;
                gsfCharge = ntuple->Electron_gsfCharge;
                etaSC = ntuple->Electron_etaSC;
                phiSC = ntuple->Electron_phiSC;
                etaWidth = ntuple->Electron_etaWidth;
                phiWidth = ntuple->Electron_phiWidth;
                dEtaIn = ntuple->Electron_dEtaIn;
                dEtaInSeed = ntuple->Electron_dEtaInSeed;
                dPhiIn = ntuple->Electron_dPhiIn;
                sigmaIEtaIEta = ntuple->Electron_sigmaIEtaIEta;
                Full5x5_SigmaIEtaIEta = ntuple->Electron_Full5x5_SigmaIEtaIEta;
                HoverE = ntuple->Electron_HoverE;
                fbrem = ntuple->Electron_fbrem;
                eOverP = ntuple->Electron_eOverP;
                InvEminusInvP = ntuple->Electron_InvEminusInvP;
                dxyVTX = ntuple->Electron_dxyVTX;
                dzVTX = ntuple->Electron_dzVTX;
                dxy = ntuple->Electron_dxy;
                dz = ntuple->Electron_dz;
                dxyBS = ntuple->Electron_dxyBS;
                dzBS = ntuple->Electron_dzBS;
                chIso03 = ntuple->Electron_chIso03;
                nhIso03 = ntuple->Electron_nhIso03;
                phIso03 = ntuple->Electron_phIso03;
                ChIso03FromPU = ntuple->Electron_ChIso03FromPU;
                mHits = ntuple->Electron_mHits;
                EnergySC = ntuple->Electron_EnergySC;
                preEnergySC = ntuple->Electron_preEnergySC;
                rawEnergySC = ntuple->Electron_rawEnergySC;
                etSC = ntuple->Electron_etSC;
                E15 = ntuple->Electron_E15;
                E25 = ntuple->Electron_E25;
                E55 = ntuple->Electron_E55;
                RelPFIso_dBeta = ntuple->Electron_RelPFIso_dBeta;
                RelPFIso_Rho = ntuple->Electron_RelPFIso_Rho;
                r9 = ntuple->Electron_r9;
                ecalDriven = ntuple->Electron_ecalDriven;
                passConvVeto = ntuple->Electron_passConvVeto;
                passMediumID = ntuple->Electron_passMediumID;

                TLorentzVector v;
                v.SetPtEtaPhiM(Pt, eta, phi, M_Elec);
                Momentum = v;
        }

	// -- Ref: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 -- // by Dalmin Pai
	Bool_t isMediumElectron_2016dataFor80X()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_Full5x5_SigmaIEtaIEta()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
//			if(Full5x5_SigmaIEtaIEta < 0.00998
                                if(fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
//			if(Full5x5_SigmaIEtaIEta < 0.0298
                                if(fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_dEtaInSeed()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
//				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
//				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_dPhiIn()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
//				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
//				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_HoverE()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
//				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
//				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_RelPFIso_Rho()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
//				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
//				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_InvEminusInvP()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
//				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
//				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_mHits()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
//				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
//				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_2016dataFor80X_minus_passConvVeto()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.00998
				&& fabs(dEtaInSeed) < 0.00311
				&& fabs(dPhiIn) < 0.103
				&& HoverE < 0.253
				&& RelPFIso_Rho < 0.0695
				&& InvEminusInvP < 0.134
//				&& fabs(dxyVTX) < 0.05
//				&& fabs(dzVTX) < 0.1
				&& mHits <= 1
//				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(Full5x5_SigmaIEtaIEta < 0.0298
				&& fabs(dEtaInSeed) < 0.00609
				&& fabs(dPhiIn) < 0.045
				&& HoverE < 0.0878
				&& RelPFIso_Rho < 0.0821
				&& InvEminusInvP < 0.13
//				&& fabs(dxyVTX) < 0.1
//				&& fabs(dzVTX) < 0.2
				&& mHits <= 1
//				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	// -- Ref: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 -- //
	Bool_t isMediumElectron_Spring25ns()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(sigmaIEtaIEta < 0.0101
				&& fabs(dEtaIn) < 0.0103
				&& fabs(dPhiIn) < 0.0336
				&& HoverE < 0.0876
				&& RelPFIso_Rho < 0.0766
				&& InvEminusInvP < 0.0174
				&& fabs(dxyVTX) < 0.0118
				&& fabs(dzVTX) < 0.373
				&& mHits <= 2 
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(sigmaIEtaIEta < 0.0283
				&& fabs(dEtaIn) < 0.00733
				&& fabs(dPhiIn) < 0.114
				&& HoverE < 0.0678
				&& RelPFIso_Rho < 0.0678
				&& InvEminusInvP < 0.0898
				&& fabs(dxyVTX) < 0.0739
				&& fabs(dzVTX) < 0.602
				&& mHits <= 1
				&& passConvVeto == kTRUE
				)
				isPass = kTRUE;
		}

		return isPass;
	}

	Bool_t isMediumElectron_Spring25ns_minus_PFIso()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(sigmaIEtaIEta < 0.0101
				&& fabs(dEtaIn) < 0.0103
				&& fabs(dPhiIn) < 0.0336
				&& HoverE < 0.0876
				// && RelPFIso_Rho < 0.0766
				&& InvEminusInvP < 0.0174
				&& fabs(dxyVTX) < 0.0118
				&& fabs(dzVTX) < 0.373
				&& mHits <= 2 
				// && passConvVeto == 1 
				)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(sigmaIEtaIEta < 0.0283
				&& fabs(dEtaIn) < 0.00733
				&& fabs(dPhiIn) < 0.114
				&& HoverE < 0.0678
				// && RelPFIso_Rho < 0.0678
				&& InvEminusInvP < 0.0898
				&& fabs(dxyVTX) < 0.0739
				&& fabs(dzVTX) < 0.602
				&& mHits <= 1
				// && passConvVeto == 1 
				)
				isPass = kTRUE;
		}

		return isPass;
	}

        Bool_t isTrigMatched(NtupleHandle *nh, TString HLT, Double_t *HLT_pT=NULL, Double_t *HLT_PS=NULL)
	{
		vector<string> *hlt_trigName = nh->HLT_trigName;
		Int_t hlt_ntrig = nh->HLT_ntrig;

		Bool_t isTrigMatch = false;
                for(Int_t k = 0; k < hlt_ntrig; k++)
		{
                        if((hlt_trigName->at((unsigned int)k)) == HLT)
			{
                            Double_t Lepton_pT = this->Pt;
                            Double_t Lepton_eta = this->eta;
                            Double_t Lepton_phi = this->phi;
                            Double_t Trig_pT = nh->HLT_trigPt[k];
                            Double_t Trig_eta = nh->HLT_trigEta[k];
                            Double_t Trig_phi = nh->HLT_trigPhi[k];

                            Double_t dR = sqrt((Lepton_eta - Trig_eta)*(Lepton_eta - Trig_eta) + (Lepton_phi - Trig_phi)*(Lepton_phi - Trig_phi));
                            Double_t dpT = fabs(Lepton_pT - Trig_pT) / Trig_pT;
                            if(dR < 0.3 && fabs(Lepton_eta) < 2.5)
                            {
//                                cout << "HLTname: " << hlt_trigName->at((unsigned int)k) <<"    pT: " << Lepton_pT << "   HLT pT: " << Trig_pT << endl;
                                if (dpT < 0.2)
                                {
                                    isTrigMatch = true;
                                    if (HLT_pT) *HLT_pT = Trig_pT;
                                    if (HLT_PS) *HLT_PS = nh->HLT_trigPS->at(k)/*[k]*/;
                                    break;
                                }
                            }
			}
		}
		return isTrigMatch;
	}

};

class Muon : public Object
{
public:
	// -- Muon ID Variables -- //
	Int_t isGLB;
	Int_t isPF;
	Int_t isTRK;
	Int_t isSTA;
	Int_t charge;
	Double_t chi2dof;
	Int_t muonHits;
	Int_t nSegments;
	Int_t nMatches;
	Int_t trackerLayers;
	Int_t pixelHits;
	Double_t dxyVTX;
	Double_t dzVTX;

	// -- Isolations -- //
	Double_t trkiso;
	Double_t relPFiso;
	Double_t RelPFIso_dBeta;

	// -- Various Track Information -- //
	Double_t Best_pT;
	Double_t Best_pTError;
	Double_t Best_Px;
	Double_t Best_Py;
	Double_t Best_Pz;
	Double_t Best_eta;
	Double_t Best_phi;

	Double_t Inner_pT;
	Double_t Inner_pTError;
	Double_t Inner_Px;
	Double_t Inner_Py;
	Double_t Inner_Pz;
	Double_t Inner_eta;
	Double_t Inner_phi;
	TLorentzVector Momentum_Inner;

	Double_t Outer_pT;
	Double_t Outer_pTError;
	Double_t Outer_Px;
	Double_t Outer_Py;
	Double_t Outer_Pz;
	Double_t Outer_eta;
	Double_t Outer_phi;

	Double_t GLB_pT;
	Double_t GLB_pTError;
	Double_t GLB_Px;
	Double_t GLB_Py;
	Double_t GLB_Pz;
	Double_t GLB_eta;
	Double_t GLB_phi;

	Double_t TuneP_pT;
	Double_t TuneP_pTError;
	Double_t TuneP_Px;
	Double_t TuneP_Py;
	Double_t TuneP_Pz;
	Double_t TuneP_eta;
	Double_t TuneP_phi;

        // Muon ID flags
        Bool_t passLooseID;
        Bool_t passMediumID;
        Bool_t passTightID;
        Bool_t passHighPtID;

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Pt = ntuple->Muon_pT[index];
		eta = ntuple->Muon_eta[index];
		phi = ntuple->Muon_phi[index];
		isGLB = ntuple->isGLBmuon[index];
		isPF = ntuple->isPFmuon[index];
		isTRK = ntuple->isTRKmuon[index];
		charge = ntuple->Muon_charge[index];

		chi2dof = ntuple->Muon_chi2dof[index];
		muonHits = ntuple->Muon_muonHits[index];
		nSegments = ntuple->Muon_nSegments[index];
		nMatches = ntuple->Muon_nMatches[index];
		trackerLayers = ntuple->Muon_trackerLayers[index];
		pixelHits = ntuple->Muon_pixelHits[index];
		dxyVTX = ntuple->Muon_dxyVTX[index];
		dzVTX = ntuple->Muon_dzVTX[index];
		trkiso = ntuple->Muon_trkiso[index] / ntuple->Muon_pT[index];

		Double_t pfChargedIso = ntuple->Muon_PfChargedHadronIsoR04[index];
		Double_t pfNeutralIso = ntuple->Muon_PfNeutralHadronIsoR04[index];
		Double_t pfGammaIso = ntuple->Muon_PfGammaIsoR04[index];
		relPFiso = (pfChargedIso + pfNeutralIso + pfGammaIso) / ntuple->Muon_pT[index];

		Double_t pfSumPUPt = ntuple->Muon_PFSumPUIsoR04[index];
                this->RelPFIso_dBeta = (pfChargedIso + max(0.0, pfNeutralIso + pfGammaIso - 0.5*pfSumPUPt)) / this->Pt;

		Double_t px = ntuple->Muon_Px[index];
		Double_t py = ntuple->Muon_Py[index];
		Double_t pz = ntuple->Muon_Pz[index];
                Double_t Mu_E = sqrt(px*px + py*py + pz*pz + M_Mu*M_Mu);
		TLorentzVector v;
		v.SetPxPyPzE(px, py, pz, Mu_E);
		Momentum = v;

		Best_pT = ntuple->Muon_Best_pT[index];
		Best_pTError = ntuple->Muon_Best_pTError[index];
		Best_Px = ntuple->Muon_Best_Px[index];
		Best_Py = ntuple->Muon_Best_Py[index];
		Best_Pz = ntuple->Muon_Best_Pz[index];
		Best_eta = ntuple->Muon_Best_eta[index];
		Best_phi = ntuple->Muon_Best_phi[index];

		Inner_pT = ntuple->Muon_Inner_pT[index];
		Inner_pTError = ntuple->Muon_Inner_pTError[index];
		Inner_Px = ntuple->Muon_Inner_Px[index];
		Inner_Py = ntuple->Muon_Inner_Py[index];
		Inner_Pz = ntuple->Muon_Inner_Pz[index];
		Inner_eta = ntuple->Muon_Inner_eta[index];
		Inner_phi = ntuple->Muon_Inner_phi[index];

		TLorentzVector v_inner;
                Double_t Mu_Inner_E = sqrt(Inner_Px*Inner_Px + Inner_Py*Inner_Py + Inner_Pz*Inner_Pz + M_Mu*M_Mu);
		v_inner.SetPxPyPzE(Inner_Px, Inner_Py, Inner_Pz, Mu_Inner_E);
		Momentum_Inner =  v_inner;

		Outer_pT = ntuple->Muon_Outer_pT[index];
		Outer_pTError = ntuple->Muon_Outer_pTError[index];
		Outer_Px = ntuple->Muon_Outer_Px[index];
		Outer_Py = ntuple->Muon_Outer_Py[index];
		Outer_Pz = ntuple->Muon_Outer_Pz[index];
		Outer_eta = ntuple->Muon_Outer_eta[index];
		Outer_phi = ntuple->Muon_Outer_phi[index];

		GLB_pT = ntuple->Muon_GLB_pT[index];
		GLB_pTError = ntuple->Muon_GLB_pTError[index];
		GLB_Px = ntuple->Muon_GLB_Px[index];
		GLB_Py = ntuple->Muon_GLB_Py[index];
		GLB_Pz = ntuple->Muon_GLB_Pz[index];
		GLB_eta = ntuple->Muon_GLB_eta[index];
		GLB_phi = ntuple->Muon_GLB_phi[index];

		TuneP_pT = ntuple->Muon_TuneP_pT[index];
		TuneP_pTError = ntuple->Muon_TuneP_pTError[index];
		TuneP_Px = ntuple->Muon_TuneP_Px[index];
		TuneP_Py = ntuple->Muon_TuneP_Py[index];
		TuneP_Pz = ntuple->Muon_TuneP_Pz[index];
		TuneP_eta = ntuple->Muon_TuneP_eta[index];
		TuneP_phi = ntuple->Muon_TuneP_phi[index];

                // Muon ID flags
                passLooseID = ntuple->Muon_passLooseID[index];
                passMediumID = ntuple->Muon_passMediumID[index];
                passTightID = ntuple->Muon_passTightID[index];
                passHighPtID = ntuple->Muon_passHighPtID[index];
	}

        void FillFromSelectedX(LongSelectedMuMu_t *ntuple, Int_t index)
        {
                Pt = ntuple->Muon_pT->at(index);
                eta = ntuple->Muon_eta->at(index);
                phi = ntuple->Muon_phi->at(index);
                isGLB = ntuple->isGLBmuon->at(index);
                isPF = ntuple->isPFmuon->at(index);
                isTRK = ntuple->isTRKmuon->at(index);
                charge = ntuple->Muon_charge->at(index);

                chi2dof = ntuple->Muon_chi2dof->at(index);
                muonHits = ntuple->Muon_muonHits->at(index);
                nSegments = ntuple->Muon_nSegments->at(index);
                nMatches = ntuple->Muon_nMatches->at(index);
                trackerLayers = ntuple->Muon_trackerLayers->at(index);
                pixelHits = ntuple->Muon_pixelHits->at(index);
                dxyVTX = ntuple->Muon_dxyVTX->at(index);
                dzVTX = ntuple->Muon_dzVTX->at(index);
                trkiso = ntuple->Muon_trkiso->at(index) / ntuple->Muon_pT->at(index);

                Double_t pfChargedIso = ntuple->Muon_PfChargedHadronIsoR04->at(index);
                Double_t pfNeutralIso = ntuple->Muon_PfNeutralHadronIsoR04->at(index);
                Double_t pfGammaIso = ntuple->Muon_PfGammaIsoR04->at(index);
                relPFiso = (pfChargedIso + pfNeutralIso + pfGammaIso) / ntuple->Muon_pT->at(index);

                Double_t pfSumPUPt = ntuple->Muon_PFSumPUIsoR04->at(index);
                this->RelPFIso_dBeta = (pfChargedIso + max(0.0, pfNeutralIso + pfGammaIso - 0.5*pfSumPUPt)) / this->Pt;

                Double_t px = ntuple->Muon_Px->at(index);
                Double_t py = ntuple->Muon_Py->at(index);
                Double_t pz = ntuple->Muon_Pz->at(index);
                Double_t Mu_E = sqrt(px*px + py*py + pz*pz + M_Mu*M_Mu);
                TLorentzVector v;
                v.SetPxPyPzE(px, py, pz, Mu_E);
                Momentum = v;

                Best_pT = ntuple->Muon_Best_pT->at(index);
                Best_pTError = ntuple->Muon_Best_pTError->at(index);
                Best_Px = ntuple->Muon_Best_Px->at(index);
                Best_Py = ntuple->Muon_Best_Py->at(index);
                Best_Pz = ntuple->Muon_Best_Pz->at(index);
                Best_eta = ntuple->Muon_Best_eta->at(index);
                Best_phi = ntuple->Muon_Best_phi->at(index);

                Inner_pT = ntuple->Muon_Inner_pT->at(index);
                Inner_pTError = ntuple->Muon_Inner_pTError->at(index);
                Inner_Px = ntuple->Muon_Inner_Px->at(index);
                Inner_Py = ntuple->Muon_Inner_Py->at(index);
                Inner_Pz = ntuple->Muon_Inner_Pz->at(index);
                Inner_eta = ntuple->Muon_Inner_eta->at(index);
                Inner_phi = ntuple->Muon_Inner_phi->at(index);

                TLorentzVector v_inner;
                Double_t Mu_Inner_E = sqrt(Inner_Px*Inner_Px + Inner_Py*Inner_Py + Inner_Pz*Inner_Pz + M_Mu*M_Mu);
                v_inner.SetPxPyPzE(Inner_Px, Inner_Py, Inner_Pz, Mu_Inner_E);
                Momentum_Inner =  v_inner;

                Outer_pT = ntuple->Muon_Outer_pT->at(index);
                Outer_pTError = ntuple->Muon_Outer_pTError->at(index);
                Outer_Px = ntuple->Muon_Outer_Px->at(index);
                Outer_Py = ntuple->Muon_Outer_Py->at(index);
                Outer_Pz = ntuple->Muon_Outer_Pz->at(index);
                Outer_eta = ntuple->Muon_Outer_eta->at(index);
                Outer_phi = ntuple->Muon_Outer_phi->at(index);

                GLB_pT = ntuple->Muon_GLB_pT->at(index);
                GLB_pTError = ntuple->Muon_GLB_pTError->at(index);
                GLB_Px = ntuple->Muon_GLB_Px->at(index);
                GLB_Py = ntuple->Muon_GLB_Py->at(index);
                GLB_Pz = ntuple->Muon_GLB_Pz->at(index);
                GLB_eta = ntuple->Muon_GLB_eta->at(index);
                GLB_phi = ntuple->Muon_GLB_phi->at(index);

                TuneP_pT = ntuple->Muon_TuneP_pT->at(index);
                TuneP_pTError = ntuple->Muon_TuneP_pTError->at(index);
                TuneP_Px = ntuple->Muon_TuneP_Px->at(index);
                TuneP_Py = ntuple->Muon_TuneP_Py->at(index);
                TuneP_Pz = ntuple->Muon_TuneP_Pz->at(index);
                TuneP_eta = ntuple->Muon_TuneP_eta->at(index);
                TuneP_phi = ntuple->Muon_TuneP_phi->at(index);

                // Muon ID flags
                passLooseID = ntuple->Muon_passLooseID->at(index);
                passMediumID = ntuple->Muon_passMediumID->at(index);
                passTightID = ntuple->Muon_passTightID->at(index);
                passHighPtID = ntuple->Muon_passHighPtID->at(index);
        }

        void FillFromSelectedX(LongSelectedEMu_t *ntuple)
        {
                Pt = ntuple->Muon_pT;
                eta = ntuple->Muon_eta;
                phi = ntuple->Muon_phi;
                isGLB = ntuple->isGLBmuon;
                isPF = ntuple->isPFmuon;
                isTRK = ntuple->isTRKmuon;
                charge = ntuple->Muon_charge;

                chi2dof = ntuple->Muon_chi2dof;
                muonHits = ntuple->Muon_muonHits;
                nSegments = ntuple->Muon_nSegments;
                nMatches = ntuple->Muon_nMatches;
                trackerLayers = ntuple->Muon_trackerLayers;
                pixelHits = ntuple->Muon_pixelHits;
                dxyVTX = ntuple->Muon_dxyVTX;
                dzVTX = ntuple->Muon_dzVTX;
                trkiso = ntuple->Muon_trkiso / ntuple->Muon_pT;

                Double_t pfChargedIso = ntuple->Muon_PfChargedHadronIsoR04;
                Double_t pfNeutralIso = ntuple->Muon_PfNeutralHadronIsoR04;
                Double_t pfGammaIso = ntuple->Muon_PfGammaIsoR04;
                relPFiso = (pfChargedIso + pfNeutralIso + pfGammaIso) / ntuple->Muon_pT;

                Double_t pfSumPUPt = ntuple->Muon_PFSumPUIsoR04;
                this->RelPFIso_dBeta = (pfChargedIso + max(0.0, pfNeutralIso + pfGammaIso - 0.5*pfSumPUPt)) / this->Pt;

                Double_t px = ntuple->Muon_Px;
                Double_t py = ntuple->Muon_Py;
                Double_t pz = ntuple->Muon_Pz;
                Double_t Mu_E = sqrt(px*px + py*py + pz*pz + M_Mu*M_Mu);
                TLorentzVector v;
                v.SetPxPyPzE(px, py, pz, Mu_E);
                Momentum = v;

                Best_pT = ntuple->Muon_Best_pT;
                Best_pTError = ntuple->Muon_Best_pTError;
                Best_Px = ntuple->Muon_Best_Px;
                Best_Py = ntuple->Muon_Best_Py;
                Best_Pz = ntuple->Muon_Best_Pz;
                Best_eta = ntuple->Muon_Best_eta;
                Best_phi = ntuple->Muon_Best_phi;

                Inner_pT = ntuple->Muon_Inner_pT;
                Inner_pTError = ntuple->Muon_Inner_pTError;
                Inner_Px = ntuple->Muon_Inner_Px;
                Inner_Py = ntuple->Muon_Inner_Py;
                Inner_Pz = ntuple->Muon_Inner_Pz;
                Inner_eta = ntuple->Muon_Inner_eta;
                Inner_phi = ntuple->Muon_Inner_phi;

                TLorentzVector v_inner;
                Double_t Mu_Inner_E = sqrt(Inner_Px*Inner_Px + Inner_Py*Inner_Py + Inner_Pz*Inner_Pz + M_Mu*M_Mu);
                v_inner.SetPxPyPzE(Inner_Px, Inner_Py, Inner_Pz, Mu_Inner_E);
                Momentum_Inner =  v_inner;

                Outer_pT = ntuple->Muon_Outer_pT;
                Outer_pTError = ntuple->Muon_Outer_pTError;
                Outer_Px = ntuple->Muon_Outer_Px;
                Outer_Py = ntuple->Muon_Outer_Py;
                Outer_Pz = ntuple->Muon_Outer_Pz;
                Outer_eta = ntuple->Muon_Outer_eta;
                Outer_phi = ntuple->Muon_Outer_phi;

                GLB_pT = ntuple->Muon_GLB_pT;
                GLB_pTError = ntuple->Muon_GLB_pTError;
                GLB_Px = ntuple->Muon_GLB_Px;
                GLB_Py = ntuple->Muon_GLB_Py;
                GLB_Pz = ntuple->Muon_GLB_Pz;
                GLB_eta = ntuple->Muon_GLB_eta;
                GLB_phi = ntuple->Muon_GLB_phi;

                TuneP_pT = ntuple->Muon_TuneP_pT;
                TuneP_pTError = ntuple->Muon_TuneP_pTError;
                TuneP_Px = ntuple->Muon_TuneP_Px;
                TuneP_Py = ntuple->Muon_TuneP_Py;
                TuneP_Pz = ntuple->Muon_TuneP_Pz;
                TuneP_eta = ntuple->Muon_TuneP_eta;
                TuneP_phi = ntuple->Muon_TuneP_phi;

                // Muon ID flags
                passLooseID = ntuple->Muon_passLooseID;
                passMediumID = ntuple->Muon_passMediumID;
                passTightID = ntuple->Muon_passTightID;
                passHighPtID = ntuple->Muon_passHighPtID;
        }


	Bool_t isTrigMatched(NtupleHandle *nh, TString HLT)
	{
		vector<string> *hlt_trigName = nh->HLT_trigName;
		Int_t hlt_ntrig = nh->HLT_ntrig;

		Bool_t isTrigMatch = false;
                for(Int_t k = 0; k < hlt_ntrig; k++)
		{
                        if((hlt_trigName->at((unsigned int)k)) == HLT)
			{
				Double_t Lepton_eta = this->eta;
				Double_t Lepton_phi = this->phi;
				Double_t Trig_eta = nh->HLT_trigEta[k];
				Double_t Trig_phi = nh->HLT_trigPhi[k];

                                Double_t dR = sqrt((Lepton_eta - Trig_eta)*(Lepton_eta - Trig_eta) + (Lepton_phi - Trig_phi)*(Lepton_phi - Trig_phi));
                                if(dR < 0.3 && fabs(Lepton_eta) < 2.4)
				{
					isTrigMatch = true;
					break;
				}
			}
		}
		return isTrigMatch;
	}

        Bool_t isTrigMatched(LongSelectedMuMu_t *nh, TString HLT)
        {
                vector<string> *hlt_trigName = nh->HLT_trigName;
                Int_t hlt_ntrig = nh->HLT_ntrig;

                Bool_t isTrigMatch = false;
                for(Int_t k = 0; k < hlt_ntrig; k++)
                {
                        if((hlt_trigName->at((unsigned int)k)) == HLT)
                        {
                                Double_t Lepton_eta = this->eta;
                                Double_t Lepton_phi = this->phi;
                                Double_t Trig_eta = nh->HLT_trigEta->at(k);
                                Double_t Trig_phi = nh->HLT_trigPhi->at(k);

                                Double_t dR = sqrt((Lepton_eta - Trig_eta)*(Lepton_eta - Trig_eta) + (Lepton_phi - Trig_phi)*(Lepton_phi - Trig_phi));
                                if(dR < 0.3 && fabs(Lepton_eta) < 2.4)
                                {
                                        isTrigMatch = true;
                                        break;
                                }
                        }
                }
                return isTrigMatch;
        }

        Bool_t isTrigMatched(LongSelectedEMu_t *nh, TString HLT)
        {
                vector<string> *hlt_trigName = nh->HLT_trigName;
                Int_t hlt_ntrig = nh->HLT_ntrig;

                Bool_t isTrigMatch = false;
                for(Int_t k = 0; k < hlt_ntrig; k++)
                {
                        if((hlt_trigName->at((unsigned int)k)) == HLT)
                        {
                                Double_t Lepton_eta = this->eta;
                                Double_t Lepton_phi = this->phi;
                                Double_t Trig_eta = nh->HLT_trigEta->at(k);
                                Double_t Trig_phi = nh->HLT_trigPhi->at(k);

                                Double_t dR = sqrt((Lepton_eta - Trig_eta)*(Lepton_eta - Trig_eta) + (Lepton_phi - Trig_phi)*(Lepton_phi - Trig_phi));
                                if(dR < 0.3 && fabs(Lepton_eta) < 2.4)
                                {
                                        isTrigMatch = true;
                                        break;
                                }
                        }
                }
                return isTrigMatch;
        }

	Bool_t isHighPtMuon()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_isGLB()
	{
                if(  //this->isGLB == 1
			this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_muonHits()
	{
                if(  this->isGLB == 1
			//&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_nMatches()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			//&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_dpT_over_pT()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        //&& (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_dxyVTX()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
//			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_dzVTX()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
//			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_pixelHits()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
//			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isHighPtMuon_minus_trackerLayers()
	{
                if(  this->isGLB == 1
			&& this->muonHits > 0
			&& this->nMatches > 1
                        && (this->Best_pTError / this->Best_pT) < 0.3
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
//			&& this->trackerLayers > 5
			)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_isGLB()
	{
                if(  // this->isGLB == 1
			this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_isPF()
	{
                if(  this->isGLB == 1
			// && this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_chi2dof()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			// && this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_muonHits()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			// && this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_nMatches()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			// && this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_dxyVTX()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			// && fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_dzVTX()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			// && fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_pixelHits()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			// && this->pixelHits > 0
                        && this->trackerLayers > 5)
		{
			return 1;
		}

		return 0;
	}

	Bool_t isTightMuon_minus_trackerLayers()
	{
                if(  this->isGLB == 1
			&& this->isPF == 1
			&& this->chi2dof < 10
			&& this->muonHits > 0
			&& this->nMatches > 1
			&& fabs(this->dxyVTX) < 0.2
			&& fabs(this->dzVTX) < 0.5
			&& this->pixelHits > 0
			// && this->trackerLayers > 5 
			)
		{
			return 1;
		}

		return 0;
	}

};

class Photon : public Object
{
public:
	Int_t hasPixelSeed;
	Double_t etaSC;
	Double_t phiSC;
	Double_t HoverE;
	Double_t Full5x5_SigmaIEtaIEta;
	Double_t ChIso;
	Double_t NhIso;
	Double_t PhIso;
	Double_t ChIsoWithEA;
	Double_t NhIsoWithEA;
	Double_t PhIsoWithEA;

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Pt = ntuple->Photon_pT[index];
		eta = ntuple->Photon_eta[index];
		phi = ntuple->Photon_phi[index];

		hasPixelSeed = ntuple->Photon_hasPixelSeed[index];
		etaSC = ntuple->Photon_etaSC[index];
		phiSC = ntuple->Photon_phiSC[index];
		HoverE = ntuple->Photon_HoverE[index];
		Full5x5_SigmaIEtaIEta = ntuple->Photon_Full5x5_SigmaIEtaIEta[index];
		ChIso = ntuple->Photon_ChIso[index];
		NhIso = ntuple->Photon_NhIso[index];
		PhIso = ntuple->Photon_PhIso[index];
		ChIsoWithEA = ntuple->Photon_ChIsoWithEA[index];
		NhIsoWithEA = ntuple->Photon_NhIsoWithEA[index];
		PhIsoWithEA = ntuple->Photon_PhIsoWithEA[index];
	}

	Bool_t isMediumPhoton()
	{
		Bool_t isPass = kFALSE;

		Bool_t isBarrel = kFALSE;
		Bool_t isEndcap = kFALSE;
                if(fabs(etaSC) <= 1.479)
			isBarrel = kTRUE;
                else if(fabs(etaSC) > 1.479 && fabs(etaSC) < 2.5)
			isEndcap = kTRUE;

                if(isBarrel == kTRUE)
		{
                        if(HoverE < 0.05
				&& Full5x5_SigmaIEtaIEta < 0.0102
				&& ChIso < 1.37
				&& NhIsoWithEA < 1.06 + 0.014*Pt + 0.000019*Pt*Pt
                                && PhIsoWithEA < 0.28 + 0.0053*Pt)
				isPass = kTRUE;
		}
                else if(isEndcap == kTRUE)
		{
                        if(HoverE < 0.05
				&& Full5x5_SigmaIEtaIEta < 0.0268
				&& ChIso < 1.10
				&& NhIsoWithEA < 2.69 + 0.0139*Pt + 0.000025*Pt*Pt
                                && PhIsoWithEA < 0.39 + 0.0034*Pt )
				isPass = kTRUE;
		}

		return isPass;		
	}
	
};

class Jet : public Object
{
public:
	Double_t Charge;
	Int_t Flavor;

	Double_t bTag;
	Double_t CHfrac;
	Double_t NHfrac;
	Double_t NHEMfrac;
	Double_t CHEMfrac;
	Int_t CHmulti;
	Int_t NHmulti;

	void FillFromNtuple(NtupleHandle *ntuple, Int_t index)
	{
		Pt = ntuple->Jet_pT[index];
		eta = ntuple->Jet_eta[index];
		phi = ntuple->Jet_phi[index];
		Charge = ntuple->Jet_Charge[index];
		Flavor = ntuple->Jet_Flavor[index];
		bTag = ntuple->Jet_bTag[index];
		CHfrac = ntuple->Jet_CHfrac[index];
		NHfrac = ntuple->Jet_NHfrac[index];
		NHEMfrac = ntuple->Jet_NHEMfrac[index];
		CHEMfrac = ntuple->Jet_CHEMfrac[index];
		CHmulti = ntuple->Jet_CHmulti[index];
		NHmulti = ntuple->Jet_NHmulti[index];
	}

};

