// Header for file management
// Created by Marijus Ambrozas 2018.08.14

#pragma once

#include <TGraphErrors.h>
#include <TH2.h>
#include <iostream>
#include <sstream>
#include <vector>

#define Lumi 35867 // -- from Run2016B to Run2016H, JSON. unit: /pb, Updated at 2017.07.30 -- //
#define Lumi_HLTv4p2 865.919 // -- integrated luminosity before Run 257933 -- //
#define nMassBin 43

using namespace std;

enum Process_t
{
    _None=0,
    // "Normal" processes - related to the divisions of files
    _DY_10to50, _DY_50to100, _DY_100to200, _DY_200to400, _DY_400to500, _DY_500to700, _DY_700to800, _DY_800to1000, _DY_1000to1500, _DY_1500to2000, _DY_2000to3000,
    _EndOf_DY_Normal,
    _DYMuMu_10to50, _DYMuMu_50to100, _DYMuMu_100to200, _DYMuMu_200to400, _DYMuMu_400to500, _DYMuMu_500to700, _DYMuMu_700to800, _DYMuMu_800to1000,
    _DYMuMu_1000to1500, _DYMuMu_1500to2000, _DYMuMu_2000to3000, _EndOf_DYMuMu_Normal,
    _DYEE_10to50, _DYEE_50to100, _DYEE_100to200, _DYEE_200to400, _DYEE_400to500, _DYEE_500to700, _DYEE_700to800, _DYEE_800to1000, _DYEE_1000to1500,
    _DYEE_1500to2000, _DYEE_2000to3000, _EndOf_DYEE_Normal,
    _EndOf_MCsignal_Normal,
    _DYTauTau_10to50, _DYTauTau_50toInf, _EndOf_DYTauTau_Normal,
    _ttbar, _ttbar_700to1000, _ttbar_1000toInf, _EndOf_ttbar_Normal,
    _tW, _tbarW, _ZZ, _WZ, _WW, _EndOf_VVnST_Normal,
    _WJets, _WJets_ext2v5, _EndOf_WJets_Normal,
    _QCDMuEnriched_15to20, _QCDMuEnriched_20to30, _QCDMuEnriched_30to50, _QCDMuEnriched_50to80, _QCDMuEnriched_80to120, _QCDMuEnriched_120to170,
    _QCDMuEnriched_170to300, _QCDMuEnriched_300to470, _QCDMuEnriched_470to600, _QCDMuEnriched_600to800, _QCDMuEnriched_800to1000, _QCDMuEnriched_1000toInf,
    _EndOf_QCDMuEnriched_Normal,
    _QCDEMEnriched_20to30, _QCDEMEnriched_30to50, _QCDEMEnriched_50to80, _QCDEMEnriched_80to120, _QCDEMEnriched_120to170, _QCDEMEnriched_170to300,
    _QCDEMEnriched_300toInf, _EndOf_QCDEMEnriched_Normal,
    _GJets_20to100, _GJets_100to200, _GJets_200to500, _GJets_500to1000, _GJets_1000to2000, _GJets_2000to5000, _EndOf_GJets_Normal,
    _EndOf_MCbkg_Normal,
    _DoubleEG_B, _DoubleEG_C, _DoubleEG_D, _DoubleEG_E, _DoubleEG_F, _DoubleEG_G, _DoubleEG_H,
    _EndOf_DoubleEG_Normal,
    _SingleMuon_B, _SingleMuon_C, _SingleMuon_D, _SingleMuon_E, _SingleMuon_F, _SingleMuon_G, _SingleMuon_H,
    _EndOf_SingleMuon_Normal,
    _SingleElectron_B, _SingleElectron_C, _SingleElectron_D, _SingleElectron_E, _SingleElectron_F, _SingleElectron_G, _SingleElectron_H,
    _EndOf_SingleElectron_Normal,
    _SinglePhoton_B, _SinglePhoton_C, _SinglePhoton_D, _SinglePhoton_E, _SinglePhoton_F, _SinglePhoton_G, _SinglePhoton_H,
    _EndOf_SinglePhoton_Normal,
    _EndOf_Data_Normal,
    // "Special" processes - similar processes are combined
    _DY_Full, _DYMuMu_Full, _DYEE_Full,
    _EndOf_MCsignal_Special,
    _DYTauTau_Full, _ttbar_Full, _VVnST, _WJets_Full, _QCDMuEnriched_Full, _QCDEMEnriched_Full, _GJets_Full, _bkg_Full,
    _EndOf_MCbkg_Special, // there is no WJets in bkgSpecial
    _DoubleEG_Full, _SingleMuon_Full, _SingleElectron_Full, _SinglePhoton_Full,
    _EndOf_Data_Special,
    // Processes for testing at local pc
    _Test_MuMu, _Test_EE, _Test_EMu,
    _EndOf_Test,
    // Alternatively generated MC processes for comparing
    _A_DY_50to100, _A_DY_100to250, _A_DY_250to400, _A_DY_400to650, _A_DY_650toInf, _EndOf_A_DY_Normal,
    _A_WJets, _A_ZZ, _A_WZ, _A_WW, _EndOf_A_MCbkg_Normal,
    _A_DY_Full,
    _EndOf_Alternatives
};

inline
Process_t previous (Process_t pr)    // Processes that begin with "EndOf" will be skipped by this
{
  if (pr == _None)
      return pr;
  else if (pr == _DYTauTau_10to50 || pr == _DoubleEG_B || pr == _DY_Full)
      return Process_t(int(pr)-3);
  else if (pr == _DYMuMu_10to50 || pr == _DYEE_10to50 || pr == _EndOf_MCsignal_Normal || pr == _ttbar || pr == _tW ||
           pr == _WJets || pr == _QCDMuEnriched_15to20 || pr == _QCDEMEnriched_20to30 || pr == _GJets_20to100 || pr == _EndOf_MCbkg_Normal ||
           pr == _SingleMuon_B || pr == _SingleElectron_B || pr == _SinglePhoton_B || pr == _EndOf_Data_Normal || pr == _DYTauTau_Full ||
           pr == _DoubleEG_Full || pr == _Test_MuMu || pr == _A_DY_50to100 || pr == _A_WJets || pr == _A_DY_Full)
      return Process_t(int(pr)-2);
  else
      return Process_t(int(pr)-1);
}


inline
Process_t next (Process_t pr)    // Processes that begin with "EndOf" will be skipped by this
{
  if (pr == _EndOf_Alternatives)
      return pr;
  else if (pr == _DYEE_2000to3000 || pr == _GJets_2000to5000 || pr == _SinglePhoton_H)
      return Process_t(int(pr)+3);
  else if (pr == _DY_2000to3000 || pr == _DYMuMu_2000to3000 || pr == _EndOf_DYEE_Normal || pr == _DYTauTau_50toInf || pr == _ttbar_1000toInf ||
           pr == _WW || pr == _WJets_ext2v5 || pr == _QCDMuEnriched_1000toInf || pr == _QCDEMEnriched_300toInf || pr == _EndOf_GJets_Normal ||
           pr == _DoubleEG_H || pr == _SingleMuon_H || pr == _SingleElectron_H || pr == _EndOf_SinglePhoton_Normal || pr == _DYEE_Full ||
           pr == _bkg_Full || pr == _SinglePhoton_Full || pr == _Test_EMu || pr == _A_DY_650toInf || pr == _A_WW)
      return Process_t(int(pr)+2);
  else
      return Process_t(int(pr)+1);
}


class FileMgr
{
public:

        Process_t CurrentProc;
        vector<TString> Tag;
        vector<TString> FullLocation;
        vector<TString> FileLocation;
        vector<TString> NtupleName;
        vector<TString> TreeName;
        vector<Double_t> Xsec;
        vector<Double_t> Wsum;
        vector<Double_t> nEvents;
        vector<Int_t> nNtuples;

        TString BaseLocation;
        TString Type;       
        Bool_t isMC;       

        map<Process_t, TString> Procname;

        // -- Constructor -- //
        FileMgr (Process_t pr = _None);

        vector<Process_t> FindProc (TString search, Bool_t notify = kTRUE, Bool_t instaGet = kFALSE);
        void NextProc ();
        void SetProc (Process_t pr = _None, Bool_t ClearOld = kTRUE);
        void ClearProc ();
        void SetupChain(Int_t i_tuple, TChain *chain);

private:
        Bool_t namesSet = kFALSE;
        Bool_t processesChecked = kFALSE;

        void PrepareProcNames ();
        void CheckProcesses ();

};// end of class definition


// ---------- Constructor ---------- //

FileMgr::FileMgr (Process_t pr)
{
    if (namesSet == kFALSE) { this->PrepareProcNames(); namesSet = kTRUE; }
    if (processesChecked == kFALSE) { this->CheckProcesses(); processesChecked = kTRUE; }
    CurrentProc = pr;
    this->SetProc(CurrentProc, kTRUE);
}


// ----------- Functions ----------- //

void FileMgr::NextProc()
{
    CurrentProc = next(CurrentProc);
    this->SetProc(CurrentProc, kTRUE);
}


void FileMgr::ClearProc()
{
    if (CurrentProc != _None)
    {
        CurrentProc = _None;
        BaseLocation = "";
        Type = "";
        isMC = kFALSE;
        this->SetProc(CurrentProc, kTRUE);
    }
}


void FileMgr::SetProc (Process_t pr, Bool_t ClearOld)
{
    if (ClearOld == kTRUE)
    {
        Tag.clear();
        FullLocation.clear();
        FileLocation.clear();
        TreeName.clear();
        NtupleName.clear();
        Xsec.clear();
        Wsum.clear();
        nEvents.clear();
        nNtuples.clear();
    }
    TString Location;
    CurrentProc = pr;

    if (pr != _None)
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/kplee/DYntuple/v2.8/";

    if (pr == _DY_10to50)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M10to50_v1"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(30920476);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237); // One ntuple missing in v2.7 (236)
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // v2.8 name
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_M10to50_v2"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_M10to50_ext1v1"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);
    }
    else if(pr == _DY_50to100)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M50to100"); Xsec.push_back(1873.52 * 3.0); Wsum.push_back(78026576.0/*81780984.0*/); nEvents.push_back(121265767/*122055296*/);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_100to200)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M100to200"); Xsec.push_back(76.2401 * 3.0); Wsum.push_back(9536931.0); nEvents.push_back(14705939/*38422582*/);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        // Not found in v2.6
//        Tag.push_back("DY_M100to200_ext"); Xsec.push_back(76.2401 * 3.0); Wsum.push_back(703034 + 9607589); nEvents.push_back(15120677);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/";
//        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_200to400)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M200to400"); Xsec.push_back(2.67606 * 3.0); Wsum.push_back(1679621.0); nEvents.push_back(2925885);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_400to500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M400to500"); Xsec.push_back(0.139728 * 3.0); Wsum.push_back(151190.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_500to700)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M500to700"); Xsec.push_back(0.0792496 * 3.0); Wsum.push_back(144096.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_700to800)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M700to800"); Xsec.push_back(0.0123176 * 3.0); Wsum.push_back(136892.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_800to1000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M800to1000"); Xsec.push_back(0.01042 * 3.0); Wsum.push_back(131586.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_1000to1500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M1000to1500"); Xsec.push_back(0.00552772 * 3.0); Wsum.push_back(120010.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_1500to2000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M1500to2000"); Xsec.push_back(0.000741613 * 3.0); Wsum.push_back(111709.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DY_2000to3000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M2000to3000"); Xsec.push_back(0.000178737 * 3.0); Wsum.push_back(101298.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DY_Full)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_M10to50_v1"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(30920476);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_M10to50_v2"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_M10to50_ext1v1"); Xsec.push_back(6016.88 * 3.0); Wsum.push_back(22498434.0 + 47946333.0 + 29386420.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_M50to100"); Xsec.push_back(1873.52 * 3.0); Wsum.push_back(78026576.0); nEvents.push_back(121265767);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M100to200"); Xsec.push_back(76.2401 * 3.0); Wsum.push_back(9536931.0); nEvents.push_back(14705939);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M200to400"); Xsec.push_back(2.67606 * 3.0); Wsum.push_back(1679621.0); nEvents.push_back(295242);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M400to500"); Xsec.push_back(0.139728 * 3.0); Wsum.push_back(151190.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M500to700"); Xsec.push_back(0.0792496 * 3.0); Wsum.push_back(144096.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M700to800"); Xsec.push_back(0.0123176 * 3.0); Wsum.push_back(136892.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M800to1000"); Xsec.push_back(0.01042 * 3.0); Wsum.push_back(131586.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M1000to1500"); Xsec.push_back(0.00552772 * 3.0); Wsum.push_back(120010.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M1500to2000"); Xsec.push_back(0.000741613 * 3.0); Wsum.push_back(111709.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_M2000to3000"); Xsec.push_back(0.000178737 * 3.0); Wsum.push_back(101298.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYMuMu_10to50) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322.0); nEvents.push_back(30920476); // UPDATE WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DYMuMu_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DYMuMu_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);
    }
    else if(pr == _DYMuMu_50to100) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M50to100"); Xsec.push_back(1873.52); Wsum.push_back(26175605.0); nEvents.push_back(122055296);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_100to200)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M100to200"); Xsec.push_back(76.2401); Wsum.push_back(3176762.0); nEvents.push_back(14705939); // CHECK
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_200to400)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M200to400"); Xsec.push_back(2.67606); Wsum.push_back(560322.0); nEvents.push_back(295242);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_400to500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M400to500"); Xsec.push_back(0.139728); Wsum.push_back(50136.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_500to700)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M500to700"); Xsec.push_back(0.0792496); Wsum.push_back(48188.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_700to800)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M700to800"); Xsec.push_back(0.0123176); Wsum.push_back(44984.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_800to1000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M800to1000"); Xsec.push_back(0.01042); Wsum.push_back(43496.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_1000to1500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M1000to1500"); Xsec.push_back(0.00552772); Wsum.push_back(40110.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_1500to2000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M1500to2000"); Xsec.push_back(0.000741613); Wsum.push_back(37176.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYMuMu_2000to3000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M2000to3000"); Xsec.push_back(0.000178737); Wsum.push_back(33360.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYMuMu_Full)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYMuMu_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322); nEvents.push_back(30920476); // FIX WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DYMuMu_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DYMuMu_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7512584.0 + 16016651.0 + 9815322.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DYMuMu_M50to100"); Xsec.push_back(1873.52); Wsum.push_back(26175605.0); nEvents.push_back(122055296);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M100to200"); Xsec.push_back(76.2401); Wsum.push_back(3176762.0); nEvents.push_back(14705939);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M200to400"); Xsec.push_back(2.67606); Wsum.push_back(560322.0); nEvents.push_back(295242);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M400to500"); Xsec.push_back(0.139728); Wsum.push_back(50136.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M500to700"); Xsec.push_back(0.0792496); Wsum.push_back(48188.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M700to800"); Xsec.push_back(0.0123176); Wsum.push_back(44984.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M800to1000"); Xsec.push_back(0.01042); Wsum.push_back(43496.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M1000to1500"); Xsec.push_back(0.00552772); Wsum.push_back(40110.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M1500to2000"); Xsec.push_back(0.000741613); Wsum.push_back(37176.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYMuMu_M2000to3000"); Xsec.push_back(0.000178737); Wsum.push_back(33360.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYEE_10to50)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(30920476); // FIX WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_50to100) // Only EE evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M50to100"); Xsec.push_back(1873.52); Wsum.push_back(26166194.0); nEvents.push_back(122055296);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_100to200)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M100to200"); Xsec.push_back(76.2401); Wsum.push_back(3179506.0); nEvents.push_back(14705939);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_200to400)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M200to400"); Xsec.push_back(2.67606); Wsum.push_back(560818.0); nEvents.push_back(295242);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_400to500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M400to500"); Xsec.push_back(0.139728); Wsum.push_back(50420.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_500to700)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M500to700"); Xsec.push_back(0.0792496); Wsum.push_back(48039.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_700to800)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M700to800"); Xsec.push_back(0.0123176); Wsum.push_back(46114.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_800to1000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M800to1000"); Xsec.push_back(0.01042); Wsum.push_back(44256.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _DYEE_1000to1500)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M1000to1500"); Xsec.push_back(0.00552772); Wsum.push_back(39712.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYEE_1500to2000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M1500to2000"); Xsec.push_back(0.000741613); Wsum.push_back(37287.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYEE_2000to3000)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M2000to3000"); Xsec.push_back(0.000178737); Wsum.push_back(34031.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYEE_Full)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYEE_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(30920476); // FIX WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7512656.0 + 16016761.0 + 9811434.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M50to100"); Xsec.push_back(1873.52); Wsum.push_back(26166194.0); nEvents.push_back(122055296);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M100to200"); Xsec.push_back(76.2401); Wsum.push_back(3179506.0); nEvents.push_back(14705939);
//        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/";
        Location = "DYLL_M100to200_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(213);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M200to400"); Xsec.push_back(2.67606); Wsum.push_back(560818.0); nEvents.push_back(295242);
//        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/";
        Location = "DYLL_M200to400_fixed/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(75);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M400to500"); Xsec.push_back(0.139728); Wsum.push_back(50420.0); nEvents.push_back(287262);
//        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/";
        Location = "DYLL_M400to500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(4);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M500to700"); Xsec.push_back(0.0792496); Wsum.push_back(48039.0); nEvents.push_back(280940);
//        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/";
        Location = "DYLL_M500to700/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(16);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M700to800"); Xsec.push_back(0.0123176); Wsum.push_back(46114.0); nEvents.push_back(276234);
//        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/";
        Location = "DYLL_M700to800/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M800to1000"); Xsec.push_back(0.01042); Wsum.push_back(44256.0); nEvents.push_back(271768);
//        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/";
        Location = "DYLL_M800to1000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(10);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M1000to1500"); Xsec.push_back(0.00552772); Wsum.push_back(39712.0); nEvents.push_back(258620);
//        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/";
        Location = "DYLL_M1000to1500/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(11);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M1500to2000"); Xsec.push_back(0.000741613); Wsum.push_back(37287.0); nEvents.push_back(258625);
//        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/";
        Location = "DYLL_M1500to2000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(20);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYEE_M2000to3000"); Xsec.push_back(0.000178737); Wsum.push_back(34031.0); nEvents.push_back(255342);
//        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/";
        Location = "DYLL_M2000to3000/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(8);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYTauTau_10to50) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYTauTau_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(30920476); // FIX WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYTauTau_50toInf) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYTauTau_M50toInf"); Xsec.push_back(1952.68432327); Wsum.push_back(27277866.0); nEvents.push_back(122055296); //  NNLO Xsec
//        Xsec->push_back(6104.0/3.0);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DYTauTau_Full)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYTauTau_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(30920476); // FIX WSUM
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
//        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(237);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(7473194.0 + 15912921.0 + 9759664.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M50toInf"); Xsec.push_back(1952.68432327); Wsum.push_back(27277866.0); nEvents.push_back(122055296); //  NNLO Xsec
//        Xsec->push_back(6104.0/3.0);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _ttbar)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ttbar"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77081149); //M(ttbar) < 700GeV, ttbar+ttbarBackup
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/";
        Location = "ttbar/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(982);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbarBackup"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77867729);
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/";
        Location = "ttbarBackup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(1190);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _ttbar_700to1000)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ttbar_M700to1000"); Xsec.push_back(76.605); Wsum.push_back(38578324.0); nEvents.push_back(38578324);
//        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/";
        Location = "ttbar_M700to1000/";
//        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(746); // One ntuple missing in v2.7 (745)
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(746); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _ttbar_1000toInf)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ttbar_M1000toInf"); Xsec.push_back(20.578); Wsum.push_back(24561630.0); nEvents.push_back(24561630);
//        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/";
        Location = "ttbar_M1000toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(518);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _ttbar_Full)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ttbar"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77081149); //M(ttbar) < 700GeV, ttbar+ttbarBackup
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/";
        Location = "ttbar/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(982);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbarBackup"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77867729);
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/";
        Location = "ttbarBackup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(1190);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbar_M700to1000"); Xsec.push_back(76.605); Wsum.push_back(38578324.0); nEvents.push_back(38578324);
//        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/";
        Location = "ttbar_M700to1000/";
//        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(746);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(746); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbar_M1000toInf"); Xsec.push_back(20.578); Wsum.push_back(24561630.0); nEvents.push_back(24561630);
//        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/";
        Location = "ttbar_M1000toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(518);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _tW)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("tW"); Xsec.push_back(35.85); Wsum.push_back(6952830.0); nEvents.push_back(6952830);
//        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/";
        Location = "ST_tW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(161);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _tbarW)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("tbarW"); Xsec.push_back(35.85); Wsum.push_back(6933093.0); nEvents.push_back(6933093);
//        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/";
        Location = "ST_tbarW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(144);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _ZZ)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ZZ"); Xsec.push_back(16.523); Wsum.push_back(998034.0); nEvents.push_back(998034);
//        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/";
        Location = "ZZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(66);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _WZ)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WZ"); Xsec.push_back(47.13); Wsum.push_back(2995828.0); nEvents.push_back(2995828);
//        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/";
        Location = "WZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _WW)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WW"); Xsec.push_back(118.7); Wsum.push_back(6987123.0); nEvents.push_back(6987123);
//        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/";
        Location = "WW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _VVnST)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("tW"); Xsec.push_back(35.85); Wsum.push_back(6952830.0); nEvents.push_back(6952830);
//        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/";
        Location = "ST_tW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(161);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("tbarW"); Xsec.push_back(35.85); Wsum.push_back(6933093.0); nEvents.push_back(6933093);
//        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/";
        Location = "ST_tbarW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(144);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ZZ"); Xsec.push_back(16.523); Wsum.push_back(998034.0); nEvents.push_back(998034);
//        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/";
        Location = "ZZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(66);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WZ"); Xsec.push_back(47.13); Wsum.push_back(2995828.0); nEvents.push_back(2995828);
//        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/";
        Location = "WZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WW"); Xsec.push_back(118.7); Wsum.push_back(6987123.0); nEvents.push_back(6987123);
//        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/";
        Location = "WW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _WJets) // amcatnlo
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WJetsToLNu"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0/*253933112.0*/); nEvents.push_back(24120297);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/";
        Location = "WJetsToLNu_amcatnlo/";
//        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(300); // One ntuple missing in v2.7 (299)
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(300); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJetsToLNu_ext"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0/*253933112.0*/); nEvents.push_back(235575647);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/";
        Location = "WJetsToLNu_amcatnlo_ext/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(2672);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _WJets_ext2v5)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("WJetsToLNu_ext2v5"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0/*253933112.0*/); nEvents.push_back(371209392);
        Location = "WJetsToLNu_amcatnlo_ext2v5/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(3568);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _WJets_Full) // amcatnlo
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WJetsToLNu"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(24120297);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/";
        Location = "WJetsToLNu_amcatnlo/";
//        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(300);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(300); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJetsToLNu_ext"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(235575647);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/";
        Location = "WJetsToLNu_amcatnlo_ext/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(2672);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);;

        Tag.push_back("WJetsToLNu_ext2v5"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(371209392);
        Location = "WJetsToLNu_amcatnlo_ext2v5/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(3568);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_15to20)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt15to20"); Xsec.push_back(1273190000.0 * 0.003); Wsum.push_back(4141251.0); nEvents.push_back(4141251);
        Location = "QCDMuEnriched_Pt15to20/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(22);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_20to30)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt20to30"); Xsec.push_back(558528000.0 * 0.0053); Wsum.push_back(31475154.0); nEvents.push_back(31475154);
        Location = "QCDMuEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(152);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_30to50)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt30to50"); Xsec.push_back(139803000.0 * 0.01182); Wsum.push_back(37147024.0); nEvents.push_back(37147024);
        Location = "QCDMuEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(302);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_50to80)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt50to80"); Xsec.push_back(19222500.0 * 0.02276); Wsum.push_back(19806914.0); nEvents.push_back(19806914);
        Location = "QCDMuEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(96);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_80to120)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt80to120"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(13786969);
        Location = "QCDMuEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(146); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt80to120_ext1"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(9797243);
        Location = "QCDMuEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(113);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_120to170)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt120to170"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(8042720);
        Location = "QCDMuEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt120to170_backup"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(11938137);
        Location = "QCDMuEnriched_Pt120to170_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(137);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_170to300)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt170to300"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(7947158);
        Location = "QCDMuEnriched_Pt170to300/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(53);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_ext1"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(9403070);
        Location = "QCDMuEnriched_Pt170to300_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_backup"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(19607775);
        Location = "QCDMuEnriched_Pt170to300_backup/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(383); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_300to470)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt300to470"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(7937587);
        Location = "QCDMuEnriched_Pt300to470/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(67);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext1"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(16452587);
        Location = "QCDMuEnriched_Pt300to470_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext2"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(24605502);
        Location = "QCDMuEnriched_Pt300to470_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(169);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_470to600)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt470to600"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(3851523);
        Location = "QCDMuEnriched_Pt470to600/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext1"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(5663755);
        Location = "QCDMuEnriched_Pt470to600_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext2"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(9847662);
        Location = "QCDMuEnriched_Pt470to600_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_600to800)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt600to800"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(4010135);
        Location = "QCDMuEnriched_Pt600to800/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_ext1"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(5971173);
        Location = "QCDMuEnriched_Pt600to800_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_backup"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(9756852);
        Location = "QCDMuEnriched_Pt600to800_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(189);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_800to1000)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt800to1000"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(3962747);
        Location = "QCDMuEnriched_Pt800to1000/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(35);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext1"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(5838539);
        Location = "QCDMuEnriched_Pt800to1000_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext2"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(9966146);
        Location = "QCDMuEnriched_Pt800to1000_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_1000toInf)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt1000toInf"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(3990113);
        Location = "QCDMuEnriched_Pt1000toInf/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(67); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt1000toInf_ext1"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(9609820);
        Location = "QCDMuEnriched_Pt1000toInf_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(166);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _QCDMuEnriched_Full)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDMuEnriched_Pt15to20"); Xsec.push_back(1273190000.0 * 0.003); Wsum.push_back(4141251.0); nEvents.push_back(4141251);
        Location = "QCDMuEnriched_Pt15to20/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(22);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt20to30"); Xsec.push_back(558528000.0 * 0.0053); Wsum.push_back(31475154.0); nEvents.push_back(31475154);
        Location = "QCDMuEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(152);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt30to50"); Xsec.push_back(139803000.0 * 0.01182); Wsum.push_back(37147024.0); nEvents.push_back(37147024);
        Location = "QCDMuEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(302);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt50to80"); Xsec.push_back(19222500.0 * 0.02276); Wsum.push_back(19806914.0); nEvents.push_back(19806914);
        Location = "QCDMuEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(96);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt80to120"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(13786969);
        Location = "QCDMuEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(146);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt80to120_ext1"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(9797243);
        Location = "QCDMuEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(113);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt120to170"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(8042720);
        Location = "QCDMuEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt120to170_backup"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(11938137);
        Location = "QCDMuEnriched_Pt120to170_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(137);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(7947158);
        Location = "QCDMuEnriched_Pt170to300/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(53);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_ext1"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(9403070);
        Location = "QCDMuEnriched_Pt170to300_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_backup"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(19607775);
        Location = "QCDMuEnriched_Pt170to300_backup/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(383);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(7937587);
        Location = "QCDMuEnriched_Pt300to470/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(67);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext1"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(16452587);
        Location = "QCDMuEnriched_Pt300to470_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext2"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(24605502);
        Location = "QCDMuEnriched_Pt300to470_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(169);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(3851523);
        Location = "QCDMuEnriched_Pt470to600/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext1"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(5663755);
        Location = "QCDMuEnriched_Pt470to600_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext2"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(9847662);
        Location = "QCDMuEnriched_Pt470to600_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(4010135);
        Location = "QCDMuEnriched_Pt600to800/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_ext1"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(5971173);
        Location = "QCDMuEnriched_Pt600to800_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_backup"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(9756852);
        Location = "QCDMuEnriched_Pt600to800_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(189);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(3962747);
        Location = "QCDMuEnriched_Pt800to1000/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(35);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext1"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(5838539);
        Location = "QCDMuEnriched_Pt800to1000_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext2"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(9966146);
        Location = "QCDMuEnriched_Pt800to1000_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt1000toInf"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(3990113);
        Location = "QCDMuEnriched_Pt1000toInf/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(67); // name in v2.8
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt1000toInf_ext1"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(9609820);
        Location = "QCDMuEnriched_Pt1000toInf_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(166);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_20to30)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt20to30"); Xsec.push_back(557600000.0 * 0.0096); Wsum.push_back(9218952.0); nEvents.push_back(9218952);
//        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/";
        Location = "QCDEMEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(112);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_30to50)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt30to50"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(4730195);
//        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/";
        Location = "QCDEMEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(34);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt30to50_ext1"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(6768384);
//        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/";
        Location = "QCDEMEnriched_Pt30to50_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_50to80)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt50to80"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(22337068);
//        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/";
        Location = "QCDEMEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(186);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt50to80_ext1"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(23474168);
//        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/";
        Location = "QCDEMEnriched_Pt50to80_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(148);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_80to120)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt80to120"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(35841780);
//        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/";
        Location = "QCDEMEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(271);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt80to120_ext1"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(41853502); // UPDATE
//        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/";
        Location = "QCDEMEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(279);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_120to170)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt120to170"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(35817276);
//        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/";
        Location = "QCDEMEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(316); // One ntuple missing in v2.7 (315)
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt120to170_ext1"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(41954033);
//        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/";
        Location = "QCDEMEnriched_Pt120to170_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(306);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_170to300)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt170to300"); Xsec.push_back(114000.0 * 0.165); Wsum.push_back(11540162.0); nEvents.push_back(11540162);
//        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/";
        Location = "QCDEMEnriched_Pt170to300/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_300toInf)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt300toInf"); Xsec.push_back(9000.0 * 0.15); Wsum.push_back(7373633.0); nEvents.push_back(7373633);
//        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/";
        Location = "QCDEMEnriched_Pt300toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(93);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _QCDEMEnriched_Full)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("QCDEMEnriched_Pt20to30"); Xsec.push_back(557600000.0 * 0.0096); Wsum.push_back(9218952.0); nEvents.push_back(9218952);
//        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/";
        Location = "QCDEMEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(112);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt30to50"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(4730195);
//        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/";
        Location = "QCDEMEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(34);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt30to50_ext1"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(6768384);
//        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/";
        Location = "QCDEMEnriched_Pt30to50_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt50to80"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(22337068);
//        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/";
        Location = "QCDEMEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(186);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt50to80_ext1"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(23474168);
//        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/";
        Location = "QCDEMEnriched_Pt50to80_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(148);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt80to120"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(35841780);
//        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/";
        Location = "QCDEMEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(271);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt80to120_ext1"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(41853502);
//        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/";
        Location = "QCDEMEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(279);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt120to170"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(35817276);
//        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/";
        Location = "QCDEMEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(316);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt120to170_ext1"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(41954033);
//        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/";
        Location = "QCDEMEnriched_Pt120to170_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(306);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt170to300"); Xsec.push_back(114000.0 * 0.165); Wsum.push_back(11540162.0); nEvents.push_back(11540162);
//        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/";
        Location = "QCDEMEnriched_Pt170to300/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt300toInf"); Xsec.push_back(9000.0 * 0.15); Wsum.push_back(7373633.0); nEvents.push_back(7373633);
//        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/";
        Location = "QCDEMEnriched_Pt300toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(93);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_20to100)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt20to100"); Xsec.push_back(137800.0); Wsum.push_back(1.73524e+11); nEvents.push_back(498342); // Using actual wsum (weights not changed to +-1)
        Location = "GJets_Pt20to100_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(60);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_100to200)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt100to200"); Xsec.push_back(1024.0); Wsum.push_back(1.30426e+09); nEvents.push_back(496553);
        Location = "GJets_Pt100to200_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_200to500)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt200to500"); Xsec.push_back(68.66); Wsum.push_back(8.2977e+07); nEvents.push_back(476054);
        Location = "GJets_Pt200to500_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(19);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_500to1000)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt500to1000"); Xsec.push_back(1.014); Wsum.push_back(1.42083e+06); nEvents.push_back(484865);
        Location = "GJets_Pt500to1000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(38);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_1000to2000)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt1000to2000"); Xsec.push_back(0.02092); Wsum.push_back(24997.5); nEvents.push_back(395409);
        Location = "GJets_Pt1000to2000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(56);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_2000to5000)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt2000to5000"); Xsec.push_back(7.476e-05); Wsum.push_back(89.7406); nEvents.push_back(299739);
        Location = "GJets_Pt2000to5000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(30);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _GJets_Full)
    {
        isMC = kTRUE;
        Type = "BKG";

        Tag.push_back("GammaJets_Pt20to100"); Xsec.push_back(137800.0); Wsum.push_back(479734.0); nEvents.push_back(498342);
        Location = "GJets_Pt20to100_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(60);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt100to200"); Xsec.push_back(1024.0); Wsum.push_back(485683.0); nEvents.push_back(496553);
        Location = "GJets_Pt100to200_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt200to500"); Xsec.push_back(68.66); Wsum.push_back(465906.0); nEvents.push_back(476054);
        Location = "GJets_Pt200to500_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(19);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt500to1000"); Xsec.push_back(1.014); Wsum.push_back(473831.0); nEvents.push_back(484865);
        Location = "GJets_Pt500to1000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(38);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt1000to2000"); Xsec.push_back(0.02092); Wsum.push_back(383501.0); nEvents.push_back(395409);
        Location = "GJets_Pt1000to2000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(56);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt2000to5000"); Xsec.push_back(7.476e-05); Wsum.push_back(225925.0); nEvents.push_back(299739);
        Location = "GJets_Pt2000to5000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(30);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _bkg_Full)
    {
        isMC = kTRUE;
        Type = "BKG";
//        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DYTauTau_M10to50_v1"); Xsec.push_back(6016.88); Wsum.push_back(33080379.0); nEvents.push_back(30920476);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v1/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(237);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_v2"); Xsec.push_back(6016.88); Wsum.push_back(33080379.0); nEvents.push_back(65887977);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_v2/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(765);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M10to50_ext1v1"); Xsec.push_back(6016.88); Wsum.push_back(33080379.0); nEvents.push_back(40381246);
//        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/";
        Location = "DYLL_M10to50_ext1v1/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(467);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DYTauTau_M50toInf"); Xsec.push_back(1952.68432327); Wsum.push_back(27277866.0); nEvents.push_back(122055296); //  NNLO Xsec
//        Xsec->push_back(6104.0/3.0);
//        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/";
        Location = "DYLL_M50toInf/";
        NtupleName.push_back("ntuple_skim_sig_"); nNtuples.push_back(1451);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbar"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77081149); //M(ttbar) < 700GeV, ttbar+ttbarBackup
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/";
        Location = "ttbar/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(982);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbarBackup"); Xsec.push_back(734.577); Wsum.push_back(67632273.0 + 68317507.0); nEvents.push_back(77867729);
//        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/";
        Location = "ttbarBackup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(1190);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbar_M700to1000"); Xsec.push_back(76.605); Wsum.push_back(38578324.0); nEvents.push_back(38578324);
//        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/";
        Location = "ttbar_M700to1000/";
//        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(746);
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(746);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ttbar_M1000toInf"); Xsec.push_back(20.578); Wsum.push_back(24561630.0); nEvents.push_back(24561630);
//        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/";
        Location = "ttbar_M1000toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(518);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("tW"); Xsec.push_back(35.85); Wsum.push_back(6952830.0); nEvents.push_back(6952830);
//        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/";
        Location = "ST_tW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(161);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("tbarW"); Xsec.push_back(35.85); Wsum.push_back(6933093.0); nEvents.push_back(6933093);
//        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/";
        Location = "ST_tbarW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(144);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("ZZ"); Xsec.push_back(16.523); Wsum.push_back(998034.0); nEvents.push_back(998034);
//        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/";
        Location = "ZZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(66);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WZ"); Xsec.push_back(47.13); Wsum.push_back(2995828.0); nEvents.push_back(2995828);
//        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/";
        Location = "WZ/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WW"); Xsec.push_back(118.7); Wsum.push_back(6987123.0); nEvents.push_back(6987123);
//        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/";
        Location = "WW/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJetsToLNu"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(24120297);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/";
        Location = "WJetsToLNu_amcatnlo/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(300);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJetsToLNu_ext"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(235575647);
//        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/";
        Location = "WJetsToLNu_amcatnlo_ext/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(2672);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJetsToLNu_ext2v5"); Xsec.push_back(61526.7); Wsum.push_back(16497009.0 + 161144203.0 + 253966720.0); nEvents.push_back(177139200);
        Location = "WJetsToLNu_amcatnlo_ext2v5/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(3568);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt15to20"); Xsec.push_back(1273190000.0 * 0.003); Wsum.push_back(4141251.0); nEvents.push_back(4141251);
        Location = "QCDMuEnriched_Pt15to20/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(22);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt20to30"); Xsec.push_back(558528000.0 * 0.0053); Wsum.push_back(31475154.0); nEvents.push_back(31475154);
        Location = "QCDMuEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(152);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt30to50"); Xsec.push_back(139803000.0 * 0.01182); Wsum.push_back(37147024.0); nEvents.push_back(37147024);
        Location = "QCDMuEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(302);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt50to80"); Xsec.push_back(19222500.0 * 0.02276); Wsum.push_back(19806914.0); nEvents.push_back(19806914);
        Location = "QCDMuEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(96);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt80to120"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(13786969);
        Location = "QCDMuEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(146);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt80to120_ext1"); Xsec.push_back(2758420.0 * 0.03844); Wsum.push_back(13786969.0 + 9797243.0); nEvents.push_back(9797243);
        Location = "QCDMuEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(113);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt120to170"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(8042720);
        Location = "QCDMuEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt120to170_backup"); Xsec.push_back(469797.0 * 0.05362); Wsum.push_back(8042720.0 + 11938137.0); nEvents.push_back(11938137);
        Location = "QCDMuEnriched_Pt120to170_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(137);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(7947158);
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(53);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_ext1"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(9403070);
        Location = "QCDMuEnriched_Pt170to300_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt170to300_backup"); Xsec.push_back(117989.0 * 0.07335); Wsum.push_back(7947158.0 + 9403070.0 + 19607775.0); nEvents.push_back(19607775);
        Location = "QCDMuEnriched_Pt170to300_backup/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(383);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(7937587);
        Location = "QCDMuEnriched_Pt300to470/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(67);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext1"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(16452587);
        Location = "QCDMuEnriched_Pt300to470_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(128);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt300to470_ext2"); Xsec.push_back(7820.25 * 0.10196); Wsum.push_back(7937587.0 + 16452587.0 + 24605502.0); nEvents.push_back(24605502);
        Location = "QCDMuEnriched_Pt300to470_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(169);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(3851523);
        Location = "QCDMuEnriched_Pt470to600/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext1"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(5663755);
        Location = "QCDMuEnriched_Pt470to600_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(45);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt470to600_ext2"); Xsec.push_back(645.528 * 0.12242); Wsum.push_back(3851523.0 + 5663755.0 + 9847662.0); nEvents.push_back(9847662);
        Location = "QCDMuEnriched_Pt470to600_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(4010135);
        Location = "QCDMuEnriched_Pt600to800/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(32);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_ext1"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(5971173);
        Location = "QCDMuEnriched_Pt600to800_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt600to800_backup"); Xsec.push_back(187.109 * 0.13412); Wsum.push_back(4010135.0 + 5971173.0 + 9756852.0); nEvents.push_back(9756852);
        Location = "QCDMuEnriched_Pt600to800_backup/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(189);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(3962747);
        Location = "QCDMuEnriched_Pt800to1000/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(35);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext1"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(5838539);
        Location = "QCDMuEnriched_Pt800to1000_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(50);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt800to1000_ext2"); Xsec.push_back(32.3486 * 0.14552); Wsum.push_back(3962747.0 + 5838539.0 + 9966146.0); nEvents.push_back(9966146);
        Location = "QCDMuEnriched_Pt800to1000_ext2/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(76);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt1000toInf"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(3990113);
        Location = "QCDMuEnriched_Pt1000toInf/";
        NtupleName.push_back("ntuple_mc_"); nNtuples.push_back(67);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDMuEnriched_Pt1000toInf_ext1"); Xsec.push_back(10.4305 * 0.15544); Wsum.push_back(3990113.0 + 9609820.0); nEvents.push_back(9609820);
        Location = "QCDMuEnriched_Pt1000toInf_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(166);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt20to30"); Xsec.push_back(557600000.0 * 0.0096); Wsum.push_back(9218952.0); nEvents.push_back(9218952);
        Location = "QCDEMEnriched_Pt20to30/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(112);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt30to50"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(4730195);
        Location = "QCDEMEnriched_Pt30to50/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(34);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt30to50_ext1"); Xsec.push_back(136000000.0 * 0.073); Wsum.push_back(4730195.0 + 6768384.0); nEvents.push_back(6768384);
        Location = "QCDEMEnriched_Pt30to50_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(58);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt50to80"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(22337068);
        Location = "QCDEMEnriched_Pt50to80/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(186);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt50to80_ext1"); Xsec.push_back(19800000.0 * 0.146); Wsum.push_back(22337068.0 + 23474168.0); nEvents.push_back(23474168);
        Location = "QCDEMEnriched_Pt50to80_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(148);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt80to120"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(35841780);
        Location = "QCDEMEnriched_Pt80to120/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(271);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt80to120_ext1"); Xsec.push_back(2800000.0 * 0.125); Wsum.push_back(35841780.0 + 41853502.0); nEvents.push_back(41853502);
        Location = "QCDEMEnriched_Pt80to120_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(279);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt120to170"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(35817276);
        Location = "QCDEMEnriched_Pt120to170/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(316);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt120to170_ext1"); Xsec.push_back(477000.0 * 0.132); Wsum.push_back(35817276.0 + 41954033.0); nEvents.push_back(41954033);
        Location = "QCDEMEnriched_Pt120to170_ext1/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(306);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt170to300"); Xsec.push_back(114000.0 * 0.165); Wsum.push_back(11540162.0); nEvents.push_back(11540162);
        Location = "QCDEMEnriched_Pt170to300/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(116);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("QCDEMEnriched_Pt300toInf"); Xsec.push_back(9000.0 * 0.15); Wsum.push_back(7373633.0); nEvents.push_back(7373633);
        Location = "QCDEMEnriched_Pt300toInf/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(93);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt20to100"); Xsec.push_back(137800.0); Wsum.push_back(479734.0); nEvents.push_back(498342);
        Location = "GJets_Pt20to100_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(60);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt100to200"); Xsec.push_back(1024.0); Wsum.push_back(485683.0); nEvents.push_back(496553);
        Location = "GJets_Pt100to200_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(17);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt200to500"); Xsec.push_back(68.66); Wsum.push_back(465906.0); nEvents.push_back(476054);
        Location = "GJets_Pt200to500_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(19);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt500to1000"); Xsec.push_back(1.014); Wsum.push_back(473831.0); nEvents.push_back(484865);
        Location = "GJets_Pt500to1000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(38);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt1000to2000"); Xsec.push_back(0.02092); Wsum.push_back(383501.0); nEvents.push_back(395409);
        Location = "GJets_Pt1000to2000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(56);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("GammaJets_Pt2000to5000"); Xsec.push_back(7.476e-05); Wsum.push_back(225925.0); nEvents.push_back(299739);
        Location = "GJets_Pt2000to5000_Sherpa/";
        NtupleName.push_back("ntuple_skim_mc_"); nNtuples.push_back(30);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_B)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_B"); nEvents.push_back(103625724);
        Location = "DoubleEG_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1318);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_C)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_C"); nEvents.push_back(45521797);
        Location = "DoubleEG_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_D)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_D"); nEvents.push_back(52422569);
        Location = "DoubleEG_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_E)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_E"); nEvents.push_back(47326656);
        Location = "DoubleEG_Run2016E/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_F)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_F"); nEvents.push_back(33943052);
        Location = "DoubleEG_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_G)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_G"); nEvents.push_back(71864512);
        Location = "DoubleEG_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_H)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_Hver2"); nEvents.push_back(68821231);
        Location = "DoubleEG_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_Hver3"); nEvents.push_back(2021309);
        Location = "DoubleEG_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _DoubleEG_Full)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("DoubleEG_B"); nEvents.push_back(103625724);
        Location = "DoubleEG_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1318);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_C"); nEvents.push_back(45521797);
        Location = "DoubleEG_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_D"); nEvents.push_back(52422569);
        Location = "DoubleEG_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_E"); nEvents.push_back(47326656);
        Location = "DoubleEG_Run2016E/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_F"); nEvents.push_back(33943052);
        Location = "DoubleEG_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_G"); nEvents.push_back(71864512);
        Location = "DoubleEG_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_Hver2"); nEvents.push_back(68821231);
        Location = "DoubleEG_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DoubleEG_Hver3"); nEvents.push_back(2021309);
        Location = "DoubleEG_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_B)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_B"); nEvents.push_back(108561074);
        Location = "SingleMuon_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1317);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_C)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_C"); nEvents.push_back(64715287);
        Location = "SingleMuon_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_D)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_D"); nEvents.push_back(96652779);
        Location = "SingleMuon_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_E)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_E"); nEvents.push_back(87358348);
        Location = "SingleMuon_Run2016E/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_F)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_F"); nEvents.push_back(64986568);
        Location = "SingleMuon_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_G)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_G"); nEvents.push_back(138710659);
        Location = "SingleMuon_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_H)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_Hver2"); nEvents.push_back(141936183);
        Location = "SingleMuon_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_Hver3"); nEvents.push_back(4386928);
        Location = "SingleMuon_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleMuon_Full)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cluster142.knu.ac.kr:1094//store/user/dpai/_prime_/DYntuple/v2.6/";

        Tag.push_back("SingleMuon_B"); nEvents.push_back(108561074);
        Location = "SingleMuon_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1317);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_C"); nEvents.push_back(64715287);
        Location = "SingleMuon_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_D"); nEvents.push_back(96652779);
        Location = "SingleMuon_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_E"); nEvents.push_back(87358348);
        Location = "SingleMuon_Run2016E/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_F"); nEvents.push_back(64986568);
        Location = "SingleMuon_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_G"); nEvents.push_back(138710659);
        Location = "SingleMuon_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_Hver2"); nEvents.push_back(141936183);
        Location = "SingleMuon_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleMuon_Hver3"); nEvents.push_back(4386928);
        Location = "SingleMuon_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_B)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_B_0000"); nEvents.push_back(174105617);
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_B_0001"); nEvents.push_back(174105617);
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0001/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(319);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_C)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_C"); nEvents.push_back(93325367);
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_D)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_D"); nEvents.push_back(146493465);
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_E)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_E"); nEvents.push_back(113168502);
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_F)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_F"); nEvents.push_back(70085191);
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_G)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_G_0000"); nEvents.push_back(143169219);
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_G_0001"); nEvents.push_back(143169219);
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0001/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(69);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_H)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_Hver2_0000"); nEvents.push_back(106262454);
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_Hver2_0001"); nEvents.push_back(106262454);
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0001/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(157);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_Hver3"); nEvents.push_back(3187483);
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SingleElectron_Full)
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("SingleElectron_B_0000"); nEvents.push_back(174105617);
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_B_0001"); nEvents.push_back(174105617);
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0001/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(319);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_C"); nEvents.push_back(93325367);
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_D"); nEvents.push_back(146493465);
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_E"); nEvents.push_back(113168502);
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_F"); nEvents.push_back(70085191);
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_G_0000"); nEvents.push_back(143169219);
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_G_0001"); nEvents.push_back(143169219);
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0001/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(69);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_Hver2_0000"); nEvents.push_back(106262454);
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(999);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_Hver2_0001"); nEvents.push_back(106262454);
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(157);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SingleElectron_Hver3"); nEvents.push_back(3187483);
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_B)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_B"); nEvents.push_back(53774727);
        Location = "SinglePhoton_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1318);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_C)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_C"); nEvents.push_back(22037173);
        Location = "SinglePhoton_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_D)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_D"); nEvents.push_back(29000395);
        Location = "SinglePhoton_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_E)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_E"); nEvents.push_back(21053677);
        Location = "SinglePhoton_Run2016E/";
        NtupleName.push_back("ntuple_data_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_F)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_F"); nEvents.push_back(14291439);
        Location = "SinglePhoton_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_G)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_G"); nEvents.push_back(31866122);
        Location = "SinglePhoton_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _SinglePhoton_H)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_Hver2"); nEvents.push_back(32767792);
        Location = "SinglePhoton_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_Hver3"); nEvents.push_back(866051);
        Location = "SinglePhoton_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }

    else if (pr == _SinglePhoton_Full)
    {
        isMC = kFALSE;
        Type = "DATA";

        Tag.push_back("SinglePhoton_B"); nEvents.push_back(53774727);
        Location = "SinglePhoton_Run2016B/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1318);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_C"); nEvents.push_back(22037173);
        Location = "SinglePhoton_Run2016C/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(435);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_D"); nEvents.push_back(29000395);
        Location = "SinglePhoton_Run2016D/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(729);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_E"); nEvents.push_back(21053677);
        Location = "SinglePhoton_Run2016E/";
        NtupleName.push_back("ntuple_data_"); nNtuples.push_back(620);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_F"); nEvents.push_back(14291439);
        Location = "SinglePhoton_Run2016F/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(452);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_G"); nEvents.push_back(31866122);
        Location = "SinglePhoton_Run2016G/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1068);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_Hver2"); nEvents.push_back(32767792);
        Location = "SinglePhoton_Run2016Hver2/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(1156);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("SinglePhoton_Hver3"); nEvents.push_back(866051);
        Location = "SinglePhoton_Run2016Hver3/";
        NtupleName.push_back("ntuple_skim_"); nNtuples.push_back(31);
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _Test_MuMu)
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back("ZToMuMu_M4500to6000_4"); Xsec.push_back(1.0); Wsum.push_back(10200.0); nEvents.push_back(10200);
        Location = "ZToMuMu_M4500to6000_4.root";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _Test_EE)
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back("ZToEE_M4500to6000_2"); Xsec.push_back(1.0); Wsum.push_back(39200.0); nEvents.push_back(39200);
        Location = "ZToEE_M4500to6000_2.root";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if (pr == _Test_EMu)
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back("WW_34"); Xsec.push_back(1.0); Wsum.push_back(10788.0); nEvents.push_back(10788);
//        Location = "WW_34.root";
        Location = "TEST_v26.root";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    if (pr == _A_DY_50to100)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt50to100"); Xsec.push_back(363.81428); Wsum.push_back(1); nEvents.push_back(21890432);
        Location = "DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt50to100/180326_142950/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_Pt50to100_ext3"); Xsec.push_back(363.81428); Wsum.push_back(1); nEvents.push_back(108692157);
        Location = "DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt50to100_ext3/180326_143053/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);
    }
    else if(pr == _A_DY_100to250)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt100to250"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2040596);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250/180326_143142/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext1"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2950812);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext1/180326_143238/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext2"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2991815);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext2/180326_143323/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext5"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(75702951);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext5/180326_143408/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_DY_250to400)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt250to400"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(423976);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400/180326_143530/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext1"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(590806);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext1/180326_143611/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext2"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(594317);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext2/180326_143654/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext5"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(19575946);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext5/180326_143749/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_DY_400to650)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt400to650"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(432056);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650/180326_143837/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt400to650_ext1"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(589842);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650_ext1/180326_143921/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt400to650_ext2"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(604038);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650_ext2/180326_144003/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_DY_650toInf)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt650toInf"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(430691);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf/180326_144113/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt650toInf_ext1"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(599665);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf_ext1/180326_144200/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt650toInf_ext2"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(597526);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf_ext2/180326_144249/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_DY_Full)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("DY_Pt50to100"); Xsec.push_back(363.81428); Wsum.push_back(1); nEvents.push_back(21890432);
        Location = "DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt50to100/180326_142950/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_Pt50to100_ext3"); Xsec.push_back(363.81428); Wsum.push_back(1); nEvents.push_back(108692157);
        Location = "DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt50to100_ext3/180326_143053/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation + Location);

        Tag.push_back("DY_Pt100to250"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2040596);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250/180326_143142/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext1"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2950812);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext1/180326_143238/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext2"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(2991815);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext2/180326_143323/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt100to250_ext5"); Xsec.push_back(84.014804); Wsum.push_back(1); nEvents.push_back(75702951);
        Location = "DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt100to250_ext5/180326_143408/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(423976);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400/180326_143530/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext1"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(590806);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext1/180326_143611/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext2"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(594317);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext2/180326_143654/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt250to400_ext5"); Xsec.push_back(3.228256512); Wsum.push_back(1); nEvents.push_back(19575946);
        Location = "DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt250to400_ext5/180326_143749/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt400to650"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(432056);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650/180326_143837/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt400to650_ext1"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(589842);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650_ext1/180326_143921/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt400to650_ext2"); Xsec.push_back(0.436041144); Wsum.push_back(1); nEvents.push_back(604038);
        Location = "DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt400to650_ext2/180326_144003/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt650toInf"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(430691);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf/180326_144113/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt650toInf_ext1"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(599665);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf_ext1/180326_144200/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("DY_Pt650toInf_ext2"); Xsec.push_back(0.040981055); Wsum.push_back(1); nEvents.push_back(597526);
        Location = "DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_Pt650toInf_ext2/180326_144249/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_WJets)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WJets_madgraph"); Xsec.push_back(61526.7); Wsum.push_back(1); nEvents.push_back(29705748);
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu/180326_143021/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);

        Tag.push_back("WJets_madgraph_ext"); Xsec.push_back(61526.7); Wsum.push_back(1); nEvents.push_back(57026058);
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_ext/180326_143105/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_ZZ)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("ZZ_powheg"); Xsec.push_back(1.256); Wsum.push_back(1); nEvents.push_back(6669988);
        Location = "ZZTo4L_13TeV_powheg_pythia8/crab_ZZto4L/180326_143705/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_WZ)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WZ_powheg"); Xsec.push_back(4.4297); Wsum.push_back(1); nEvents.push_back(17990100);
        Location = "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu/180326_143554/0000/ntuple_skim_120.root";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }
    else if(pr == _A_WW)
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/dpai/_v2p3_/"; // old location

        Tag.push_back("WW_powheg"); Xsec.push_back(12.178); Wsum.push_back(1); nEvents.push_back(1999000);
        Location = "WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu/180326_143324/0000/";
        TreeName.push_back("recoTree/DYTree"); FileLocation.push_back(Location); FullLocation.push_back(BaseLocation+Location);
    }

}// end of SetProc()


vector<Process_t> FileMgr::FindProc (TString search, Bool_t notify, Bool_t instaGet)
{
    TString srch = search;
    srch.ToUpper();
    srch.ReplaceAll("-", "TO");
    srch.ReplaceAll(" TO", "TO"); srch.ReplaceAll("_TO", "TO");
    srch.ReplaceAll("TO ", "TO"); srch.ReplaceAll("TO_", "TO");
    srch.ReplaceAll("DRELLYAN", "DY");
    srch.ReplaceAll("DRELL-YAN", "DY");
    srch.ReplaceAll("ZW", "WZ");
    srch.ReplaceAll("SINGLETOP", "TW");
    srch.ReplaceAll("SINGLEANTITOP", "TBARW") ;
    srch.ReplaceAll("TOPANTITOP", "TTBAR");
    srch.ReplaceAll("DIMUON", "MUMU");
    srch.ReplaceAll("DIELECTRON", "EE");
    srch.ReplaceAll("DITAU", "TAUTAU");
    srch.ReplaceAll("BACKGROUND", "BKG");
    srch.ReplaceAll("POWHEG", "ALTERNATIVE");
    srch.ReplaceAll("MADGRAPH", "ALTERNATIVE");
    srch.ReplaceAll("ALT", "ALTERNATIVE");

    if (notify == kTRUE) cout << "Searched for: " << search << "\nFound: ";
    vector<Process_t> Result;
    Process_t first = _None, last = _None;
    if (srch.Contains("DY"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            if (srch.Contains("FULL"))
            {
                Result.push_back(_A_DY_Full);
                if (notify == kTRUE) cout << Procname[_A_DY_Full] << "." << endl;
            }
            // Checking for various intervals
            if (srch.Contains("INFTO"))
                first = _EndOf_A_DY_Normal;
            else if (srch.Contains("650TO"))
                first = _A_DY_650toInf;
            else if (srch.Contains("400TO"))
                first = _A_DY_400to650;
            else if (srch.Contains("250TO"))
                first = _A_DY_250to400;
            else if (srch.Contains("100TO"))
                first = _A_DY_100to250;
            else if (srch.Contains("50TO"))
                first = _A_DY_50to100;
            else first = _None;

            if (srch.Contains("TOINF"))
                last = _A_DY_650toInf;
            else if (srch.Contains("TO650"))
                last = _A_DY_400to650;
            else if (srch.Contains("TO400"))
                last = _A_DY_250to400;
            else if (srch.Contains("TO250"))
                last = _A_DY_100to250;
            else if (srch.Contains("TO100"))
                last = _A_DY_50to100;
            else if (srch.Contains("TO50"))
                last = _EndOf_Test;
            else last = _None;

            //Swapping first with last if necessary
            if (int(first) > int(last))
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if (first == _A_DY_50to100 && last == _A_DY_650toInf)
            {
                Result.push_back(_A_DY_Full);
                if (notify == kTRUE) cout << Procname[_A_DY_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_A_DY_Full);
                if (notify == kTRUE) cout << Procname[_A_DY_Full] << "." << endl;
            }

        }// End of if(ALTERNATIVE)

        else
        {
            if (srch.Contains("MUMU"))
            {
                // Checking for various intervals
                if (srch.Contains("INFTO"))
                    first = _EndOf_DYMuMu_Normal;
                else if (srch.Contains("3000TO"))
                    first = _EndOf_DYMuMu_Normal;
                else if (srch.Contains("2000TO"))
                    first = _DYMuMu_2000to3000;
                else if (srch.Contains("1500TO"))
                    first = _DYMuMu_1500to2000;
                else if (srch.Contains("1000TO"))
                    first = _DYMuMu_1000to1500;
                else if (srch.Contains("800TO"))
                    first = _DYMuMu_800to1000;
                else if (srch.Contains("700TO"))
                    first = _DYMuMu_700to800;
                else if (srch.Contains("500TO"))
                    first = _DYMuMu_500to700;
                else if (srch.Contains("400TO"))
                    first = _DYMuMu_400to500;
                else if (srch.Contains("200TO"))
                    first = _DYMuMu_200to400;
                else if (srch.Contains("100TO"))
                    first = _DYMuMu_100to200;
                else if (srch.Contains("50TO"))
                    first = _DYMuMu_50to100;
                else if (srch.Contains("10TO"))
                    first = _DYMuMu_10to50;

                if (srch.Contains("TOINF"))
                    last = _DYMuMu_2000to3000;
                else if (srch.Contains("TO3000"))
                    last = _DYMuMu_2000to3000;
                else if (srch.Contains("TO2000"))
                    last = _DYMuMu_1500to2000;
                else if (srch.Contains("TO1500"))
                    last = _DYMuMu_1000to1500;
                else if (srch.Contains("TO1000"))
                    last = _DYMuMu_800to1000;
                else if (srch.Contains("TO800"))
                    last = _DYMuMu_700to800;
                else if (srch.Contains("TO700"))
                    last = _DYMuMu_500to700;
                else if (srch.Contains("TO500"))
                    last = _DYMuMu_400to500;
                else if (srch.Contains("TO400"))
                    last = _DYMuMu_200to400;
                else if (srch.Contains("TO200"))
                    last = _DYMuMu_100to200;
                else if (srch.Contains("TO100"))
                    last = _DYMuMu_50to100;
                else if (srch.Contains("TO50"))
                    last = _DYMuMu_10to50;
                else if (srch.Contains("TO10"))
                    last = _EndOf_DY_Normal;

                // Swapping first with last if necessary
                if (int(first)>int(last) && last!=_None)
                {
                    Process_t NewLast = Process_t(int(first)-1);
                    first = Process_t(int(last)+1);
                    last = NewLast;
                }
                if (first == _DYMuMu_10to50 && last == _DYMuMu_2000to3000)
                {
                    Result.push_back(_DYMuMu_Full);
                    if (notify == kTRUE) cout << Procname[_DYMuMu_Full] << "." << endl;
                }
                else if (first != _None && last != _None)
                {
                    for (Process_t pr=first; pr<=last; pr=next(pr))
                    {
                        Result.push_back(pr);
                        if (notify == kTRUE)
                        {
                            if (pr != last) cout << Procname[pr] << ", ";
                            else cout << Procname[pr] << "." << endl;
                        }
                    }
                }
                else
                {
                    Result.push_back(_DYMuMu_Full);
                    if (notify == kTRUE) cout << Procname[_DYMuMu_Full] << "." << endl;
                }

            }// end of if(DYMuMu)

            else if (srch.Contains("EE"))
            {
                // Checking for various intervals
                if (srch.Contains("INFTO"))
                    first = _EndOf_DYEE_Normal;
                else if (srch.Contains("3000TO"))
                    first = _EndOf_DYEE_Normal;
                else if (srch.Contains("2000TO"))
                    first = _DYEE_2000to3000;
                else if (srch.Contains("1500TO"))
                    first = _DYEE_1500to2000;
                else if (srch.Contains("1000TO"))
                    first = _DYEE_1000to1500;
                else if (srch.Contains("800TO"))
                    first = _DYEE_800to1000;
                else if (srch.Contains("700TO"))
                    first = _DYEE_700to800;
                else if (srch.Contains("500TO"))
                    first = _DYEE_500to700;
                else if (srch.Contains("400TO"))
                    first = _DYEE_400to500;
                else if (srch.Contains("200TO"))
                    first = _DYEE_200to400;
                else if (srch.Contains("100TO"))
                    first = _DYEE_100to200;
                else if (srch.Contains("50TO"))
                    first = _DYEE_50to100;
                else if (srch.Contains("10TO"))
                    first = _DYEE_10to50;

                if (srch.Contains("TOINF"))
                    last = _DYEE_2000to3000;
                else if (srch.Contains("TO3000"))
                    last = _DYEE_2000to3000;
                else if (srch.Contains("TO2000"))
                    last = _DYEE_1500to2000;
                else if (srch.Contains("TO1500"))
                    last = _DYEE_1000to1500;
                else if (srch.Contains("TO1000"))
                    last = _DYEE_800to1000;
                else if (srch.Contains("TO800"))
                    last = _DYEE_700to800;
                else if (srch.Contains("TO700"))
                    last = _DYEE_500to700;
                else if (srch.Contains("TO500"))
                    last = _DYEE_400to500;
                else if (srch.Contains("TO400"))
                    last = _DYEE_200to400;
                else if (srch.Contains("TO200"))
                    last = _DYEE_100to200;
                else if (srch.Contains("TO100"))
                    last = _DYEE_50to100;
                else if (srch.Contains("TO50"))
                    last = _DYEE_10to50;
                else if (srch.Contains("TO10"))
                    last = _EndOf_DYMuMu_Normal;

                // Swapping first with last if necessary
                if (int(first)>int(last) && last!=_None)
                {
                    Process_t NewLast = Process_t(int(first)-1);
                    first = Process_t(int(last)+1);
                    last = NewLast;
                }
                if (first == _DYEE_10to50 && last == _DYEE_2000to3000)
                {
                    Result.push_back(_DYEE_Full);
                    if (notify == kTRUE) cout << Procname[_DYEE_Full] << "." << endl;
                }
                else if (first != _None && last != _None)
                {
                    for (Process_t pr=first; pr<=last; pr=next(pr))
                    {
                        Result.push_back(pr);
                        if (notify == kTRUE)
                        {
                            if (pr != last) cout << Procname[pr] << ", ";
                            else cout << Procname[pr] << "." << endl;
                        }
                    }
                }
                else
                {
                    Result.push_back(_DYEE_Full);
                    if (notify == kTRUE) cout << Procname[_DYEE_Full] << "." << endl;
                }

            }// end of if(DYEE)

            else if (srch.Contains("TAUTAU"))
            {
                if (srch.Contains("FULL"))
                {
                    Result.push_back(_DYTauTau_Full);
                    if (notify == kTRUE) cout << Procname[_DYTauTau_Full] << "." << endl;
                }
                else
                {
                    // Checking for various intervals
                    if (srch.Contains("INFTO"))
                        first = _EndOf_DYTauTau_Normal;
                    else if (srch.Contains("3000TO"))
                        first = _EndOf_DYTauTau_Normal;
                    else if (srch.Contains("2000TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("1500TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("1000TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("800TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("700TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("500TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("400TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("200TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("100TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("50TO"))
                        first = _DYTauTau_50toInf;
                    else if (srch.Contains("10TO"))
                        first = _DYTauTau_10to50;

                    if (srch.Contains("TOINF"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO3000"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO2000"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO1500"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO1000"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO800"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO700"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO500"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO400"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO200"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO100"))
                        last = _DYTauTau_50toInf;
                    else if (srch.Contains("TO50"))
                        last = _DYTauTau_10to50;
                    else if (srch.Contains("TO10"))
                        last = _EndOf_MCsignal_Normal;

                    // Swapping first with last if necessary
                    if (int(first)>int(last) && last!=_None)
                    {
                        Process_t NewLast = Process_t(int(first)-1);
                        first = Process_t(int(last)+1);
                        last = NewLast;
                    }
                    if (first == _DYTauTau_10to50 && last == _DYTauTau_50toInf)
                    {
                        Result.push_back(_DYTauTau_Full);
                        if (notify == kTRUE) cout << Procname[_DYTauTau_Full] << "." << endl;
                    }
                    else if (first != _None && last != _None)
                    {
                        for (Process_t pr=first; pr<=last; pr=next(pr))
                        {
                            Result.push_back(pr);
                            if (notify == kTRUE)
                            {
                                if (pr != last) cout << Procname[pr] << ", ";
                                else cout << Procname[pr] << "." << endl;
                            }
                        }
                    }
                    else
                    {
                        Result.push_back(_DYTauTau_Full);
                        if (notify == kTRUE) cout << Procname[_DYTauTau_Full] << "." << endl;
                    }
                }

            }// end of if(DYTauTau)
            else
            {
                // Checking for various intervals
                if (srch.Contains("INFTO"))
                    first = _EndOf_DY_Normal;
                else if (srch.Contains("3000TO"))
                    first = _EndOf_DY_Normal;
                else if (srch.Contains("2000TO"))
                    first = _DY_2000to3000;
                else if (srch.Contains("1500TO"))
                    first = _DY_1500to2000;
                else if (srch.Contains("1000TO"))
                    first = _DY_1000to1500;
                else if (srch.Contains("800TO"))
                    first = _DY_800to1000;
                else if (srch.Contains("700TO"))
                    first = _DY_700to800;
                else if (srch.Contains("500TO"))
                    first = _DY_500to700;
                else if (srch.Contains("400TO"))
                    first = _DY_400to500;
                else if (srch.Contains("200TO"))
                    first = _DY_200to400;
                else if (srch.Contains("100TO"))
                    first = _DY_100to200;
                else if (srch.Contains("50TO"))
                    first = _DY_50to100;
                else if (srch.Contains("10TO"))
                    first = _DY_10to50;

                if (srch.Contains("TOINF"))
                    last = _DY_2000to3000;
                else if (srch.Contains("TO3000"))
                    last = _DY_2000to3000;
                else if (srch.Contains("TO2000"))
                    last = _DY_1500to2000;
                else if (srch.Contains("TO1500"))
                    last = _DY_1000to1500;
                else if (srch.Contains("TO1000"))
                    last = _DY_800to1000;
                else if (srch.Contains("TO800"))
                    last = _DY_700to800;
                else if (srch.Contains("TO700"))
                    last = _DY_500to700;
                else if (srch.Contains("TO500"))
                    last = _DY_400to500;
                else if (srch.Contains("TO400"))
                    last = _DY_200to400;
                else if (srch.Contains("TO200"))
                    last = _DY_100to200;
                else if (srch.Contains("TO100"))
                    last = _DY_50to100;
                else if (srch.Contains("TO50"))
                    last = _DY_10to50;
                else if (srch.Contains("TO10"))
                    last = _None;

                //Swapping first with last if necessary
                if (int(first) > int(last))
                {
                    Process_t NewLast = Process_t(int(first)-1);
                    first = Process_t(int(last)+1);
                    last = NewLast;
                }
                if (first == _DY_10to50 && last == _DY_2000to3000)
                {
                    Result.push_back(_DY_Full);
                    if (notify == kTRUE) cout << Procname[_DY_Full] << "." << endl;
                }
                else if (first != _None && last != _None)
                {
                    for (Process_t pr=first; pr<=last; pr=next(pr))
                    {
                        Result.push_back(pr);
                        if (notify == kTRUE)
                        {
                            if (pr != last) cout << Procname[pr] << ", ";
                            else cout << Procname[pr] << "." << endl;
                        }
                    }
                }
                else
                {
                    Result.push_back(_DY_Full);
                    if (notify == kTRUE) cout << Procname[_DY_Full] << "." << endl;
                }
            }// end of else (DY)

        }// end of else (not alternative)

    }// end of if(DrellYan)

    else if (srch.Contains("TT") || srch.Contains("TTBAR"))
    {
        if (srch.Contains("FULL"))
        {
            Result.push_back(_ttbar_Full);
            if (notify == kTRUE) cout << Procname[_ttbar_Full] << "." << endl;
        }
        // Checking for various intervals
        else
        {
            if (srch.Contains("700TO"))
                first = _ttbar_700to1000;
            else if (srch.Contains("1000TO"))
                first = _ttbar_1000toInf;
            else if (srch.Contains("INFTO"))
                first = _EndOf_ttbar_Normal;
            else first = _ttbar;

            if (srch.Contains("TO700"))
                last = _ttbar;
            else if (srch.Contains("TO1000"))
                last = _ttbar_700to1000;
            else if (srch.Contains("TO1500"))
                last = _ttbar_1000toInf;
            else if (srch.Contains("TO2000"))
                last = _ttbar_1000toInf;
            else if (srch.Contains("TO3000"))
                last = _ttbar_1000toInf;
            else if (srch.Contains("TOINF"))
                last = _ttbar_1000toInf;
            else last = _ttbar;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if (first == _ttbar && last == _ttbar_1000toInf)
            {
                Result.push_back(_ttbar_Full);
                if (notify == kTRUE) cout << Procname[_ttbar_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_ttbar_Full);
                if (notify == kTRUE) cout << Procname[_ttbar_Full] << "." << endl;
            }
        }// end of else (not "Full")

    }// end of if(ttbar)

    else if (srch.Contains("TW"))
    {
        Result.push_back(_tW);
        if (notify == kTRUE) cout << Procname[_tW] << "." << endl;
    }

    else if (srch.Contains("TBARW"))
    {
        Result.push_back(_tbarW);
        if (notify == kTRUE) cout << Procname[_tbarW] << "." << endl;
    }

    else if (srch.Contains("ZZ"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            Result.push_back(_A_ZZ);
            if (notify == kTRUE) cout << Procname[_A_ZZ] << "." << endl;
        }
        else
        {
            Result.push_back(_ZZ);
            if (notify == kTRUE) cout << Procname[_ZZ] << "." << endl;
        }
    }

    else if (srch.Contains("WZ"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            Result.push_back(_A_WZ);
            if (notify == kTRUE) cout << Procname[_A_WZ] << "." << endl;
        }
        else
        {
            Result.push_back(_WZ);
            if (notify == kTRUE) cout << Procname[_WZ] << "." << endl;
        }
    }

    else if (srch.Contains("WW"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            Result.push_back(_A_WW);
            if (notify == kTRUE) cout << Procname[_A_WW] << "." << endl;
        }
        else
        {
            Result.push_back(_WW);
            if (notify == kTRUE) cout << Procname[_WW] << "." << endl;
        }
    }

    else if (srch.Contains("DIBOSON"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            Result.push_back(_A_ZZ);
            Result.push_back(_A_WZ);
            Result.push_back(_A_WW);
            if (notify == kTRUE) cout << Procname[_A_ZZ] << ", " << Procname[_A_WZ] << ", " << Procname[_A_WW] << "." << endl;
        }
        else
        {
            Result.push_back(_ZZ);
            Result.push_back(_WZ);
            Result.push_back(_WW);
            if (notify == kTRUE) cout << Procname[_ZZ] << ", " << Procname[_WZ] << ", " << Procname[_WW] << "." << endl;
        }
    }

    else if (srch.Contains("VVNST"))
    {
        Result.push_back(_VVnST);
        if (notify == kTRUE) cout << Procname[_VVnST] << "." << endl;
    }

    else if (srch.Contains("W") && srch.Contains("JET"))
    {
        if (srch.Contains("ALTERNATIVE"))
        {
            Result.push_back(_A_WJets);
            if (notify == kTRUE) cout << Procname[_A_WJets] << "." << endl;
        }
        else
        {
            if (srch.Contains("FULL"))
            {
                Result.push_back(_WJets_Full);
                if (notify == kTRUE) cout << Procname[_WJets_Full] << "." << endl;
            }
            else if (srch.Contains("EXT"))
            {
                Result.push_back(_WJets_ext2v5);
                if (notify == kTRUE) cout << Procname[_WJets_ext2v5] << "." << endl;
            }
            else
            {
                Result.push_back(_WJets);
                if (notify == kTRUE) cout << Procname[_WJets] << "." << endl;
            }
        }
    }

    else if (srch.Contains("QCD"))
    {
        if (srch.Contains("MU"))
        {
            // Checking for various intervals
            if (srch.Contains("INFTO") || srch.Contains("INF_TO")  || srch.Contains("INF-") ||
                 srch.Contains("INF_-"))
                first = _EndOf_QCDMuEnriched_Normal;
            else if (srch.Contains("1000TO"))
                first = _QCDMuEnriched_1000toInf;
            else if (srch.Contains("800TO"))
                first = _QCDMuEnriched_800to1000;
            else if (srch.Contains("600TO"))
                first = _QCDMuEnriched_600to800;
            else if (srch.Contains("470TO"))
                first = _QCDMuEnriched_470to600;
            else if (srch.Contains("300TO"))
                first = _QCDMuEnriched_300to470;
            else if (srch.Contains("170TO"))
                first = _QCDMuEnriched_170to300;
            else if (srch.Contains("120TO"))
                first = _QCDMuEnriched_120to170;
            else if (srch.Contains("80TO"))
                first = _QCDMuEnriched_80to120;
            else if (srch.Contains("50TO"))
                first = _QCDMuEnriched_50to80;
            else if (srch.Contains("30TO"))
                first = _QCDMuEnriched_30to50;
            else if (srch.Contains("20TO"))
                first = _QCDMuEnriched_20to30;
            else if (srch.Contains("15TO"))
                first = _QCDMuEnriched_15to20;

            if (srch.Contains("TOINF"))
                last = _QCDMuEnriched_1000toInf;
            else if (srch.Contains("TO1000"))
                last = _QCDMuEnriched_800to1000;
            else if (srch.Contains("TO800"))
                last = _QCDMuEnriched_600to800;
            else if (srch.Contains("TO600"))
                last = _QCDMuEnriched_470to600;
            else if (srch.Contains("TO470"))
                last = _QCDMuEnriched_300to470;
            else if (srch.Contains("TO300"))
                last = _QCDMuEnriched_170to300;
            else if (srch.Contains("TO170"))
                last = _QCDMuEnriched_120to170;
            else if (srch.Contains("TO120"))
                last = _QCDMuEnriched_80to120;
            else if (srch.Contains("TO80"))
                last = _QCDMuEnriched_50to80;
            else if (srch.Contains("TO50"))
                last = _QCDMuEnriched_30to50;
            else if (srch.Contains("TO30"))
                last = _QCDMuEnriched_20to30;
            else if (srch.Contains("TO20"))
                last = _QCDMuEnriched_15to20;
            else if (srch.Contains("TO15"))
                last = _EndOf_WJets_Normal;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if (first == _QCDMuEnriched_15to20 && last == _QCDMuEnriched_1000toInf)
            {
                Result.push_back(_QCDMuEnriched_Full);
                if (notify == kTRUE) cout << Procname[_QCDMuEnriched_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_QCDMuEnriched_Full);
                if (notify == kTRUE) cout << Procname[_QCDMuEnriched_Full] << "." << endl;
            }

        }// end of if(MuEnriched)
        else if (srch.Contains("EM"))
        {
            // Checking for various intervals
            if (srch.Contains("INFTO"))
                first = _EndOf_QCDEMEnriched_Normal;
            else if (srch.Contains("300TO"))
                first = _QCDEMEnriched_300toInf;
            else if (srch.Contains("170TO"))
                first = _QCDEMEnriched_170to300;
            else if (srch.Contains("120TO"))
                first = _QCDEMEnriched_120to170;
            else if (srch.Contains("80TO"))
                first = _QCDEMEnriched_80to120;
            else if (srch.Contains("50TO"))
                first = _QCDEMEnriched_50to80;
            else if (srch.Contains("30TO"))
                first = _QCDEMEnriched_30to50;
            else if (srch.Contains("20TO"))
                first = _QCDEMEnriched_20to30;

            if (srch.Contains("TOINF"))
                last = _QCDEMEnriched_300toInf;
            else if (srch.Contains("TO1000"))
                last = _QCDEMEnriched_300toInf;
            else if (srch.Contains("TO800"))
                last = _QCDEMEnriched_300toInf;
            else if (srch.Contains("TO600"))
                last = _QCDEMEnriched_300toInf;
            else if (srch.Contains("TO470"))
                last = _QCDEMEnriched_300toInf;
            else if (srch.Contains("TO300"))
                last = _QCDEMEnriched_170to300;
            else if (srch.Contains("TO170"))
                last = _QCDEMEnriched_120to170;
            else if (srch.Contains("TO120"))
                last = _QCDEMEnriched_80to120;
            else if (srch.Contains("TO80"))
                last = _QCDEMEnriched_50to80;
            else if (srch.Contains("TO50"))
                last = _QCDEMEnriched_30to50;
            else if (srch.Contains("TO30"))
                last = _QCDEMEnriched_20to30;
            else if (srch.Contains("TO20"))
                last = _EndOf_QCDMuEnriched_Normal;

            // Swapping first with last if necessary
            if (int(first)>int(last)  && last!=_None)
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if (first == _QCDEMEnriched_20to30 && last == _QCDEMEnriched_300toInf)
            {
                Result.push_back(_QCDEMEnriched_Full);
                if (notify == kTRUE) cout << Procname[_QCDEMEnriched_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_QCDEMEnriched_Full);
                if (notify == kTRUE) cout << Procname[_QCDEMEnriched_Full] << "." << endl;
            }

        }// end of if(EMEnriched)

    }// end of if(QCD)

    else if (srch.Contains("G") && srch.Contains("JET"))
    {
        // Checking for various intervals
        if (srch.Contains("INFTO"))
            first = _EndOf_GJets_Normal;
        else if (srch.Contains("5000TO"))
            first = _EndOf_GJets_Normal;
        else if (srch.Contains("2000TO"))
            first = _GJets_2000to5000;
        else if (srch.Contains("1000TO"))
            first = _GJets_1000to2000;
        else if (srch.Contains("500TO"))
            first = _GJets_500to1000;
        else if (srch.Contains("200TO"))
            first = _GJets_200to500;
        else if (srch.Contains("100TO"))
            first = _GJets_100to200;
        else if (srch.Contains("20TO"))
            first = _GJets_20to100;

        if (srch.Contains("TOINF"))
            last = _GJets_2000to5000;
        else if (srch.Contains("TO5000"))
            last = _GJets_2000to5000;
        else if (srch.Contains("TO2000"))
            last = _GJets_1000to2000;
        else if (srch.Contains("TO1000"))
            last = _GJets_500to1000;
        else if (srch.Contains("TO500"))
            last = _GJets_200to500;
        else if (srch.Contains("TO200"))
            last = _GJets_100to200;
        else if (srch.Contains("TO100"))
            last = _GJets_20to100;
        else if (srch.Contains("TO20"))
            last = _EndOf_QCDEMEnriched_Normal;

        // Swapping first with last if necessary
        if (int(first)>int(last) && last!=_None)
        {
            Process_t NewLast = Process_t(int(first)-1);
            first = Process_t(int(last)+1);
            last = NewLast;
        }
        if (first == _GJets_20to100 && last == _GJets_2000to5000)
        {
            Result.push_back(_GJets_Full);
            if (notify == kTRUE) cout << Procname[_GJets_Full] << "." << endl;
        }
        else if (first != _None && last != _None)
        {
            for (Process_t pr=first; pr<=last; pr=next(pr))
            {
                Result.push_back(pr);
                if (notify == kTRUE)
                {
                    if (pr != last) cout << Procname[pr] << ", ";
                    else cout << Procname[pr] << "." << endl;
                }
            }
        }
        else
        {
            Result.push_back(_GJets_Full);
            if (notify == kTRUE) cout << Procname[_GJets_Full] << "." << endl;
        }

    }// end of if(GammaJets)

    else if (srch.Contains("DOUBLE"))
    {
        if (srch.Contains("EG"))
        {
            // Checking if search contains intervals
            if (srch.Contains("BTO"))
                first = _DoubleEG_B;
            else if (srch.Contains("CTO"))
                first = _DoubleEG_C;
            else if (srch.Contains("DTO"))
                first = _DoubleEG_D;
            else if (srch.Contains("ETO"))
                first = _DoubleEG_E;
            else if (srch.Contains("FTO"))
                first = _DoubleEG_F;
            else if (srch.Contains("GTO"))
                first = _DoubleEG_G;
            else if (srch.Contains("HTO"))
                first = _DoubleEG_H;

            if (srch.Contains("TOB"))
                last = _DoubleEG_B;
            else if (srch.Contains("TOC"))
                last = _DoubleEG_C;
            else if (srch.Contains("TOD"))
                last = _DoubleEG_D;
            else if (srch.Contains("TOE"))
                last = _DoubleEG_E;
            else if (srch.Contains("TOF"))
                last = _DoubleEG_F;
            else if (srch.Contains("TOG"))
                last = _DoubleEG_G;
            else if (srch.Contains("TOH"))
                last = _DoubleEG_H;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if (first == _DoubleEG_B && last == _DoubleEG_H)
            {
                Result.push_back(_DoubleEG_Full);
                if (notify == kTRUE) cout << Procname[_DoubleEG_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if (srch.Contains("FULL"))
            {
                Result.push_back(_DoubleEG_Full);
                if (notify == kTRUE) cout << Procname[_DoubleEG_Full] << "." << endl;
            }
            else if (srch.Contains("_B") || srch.Contains("RUNB"))
            {
                Result.push_back(_DoubleEG_B);
                if (notify == kTRUE) cout << Procname[_DoubleEG_B] << "." << endl;
            }
            else if (srch.Contains("_C") || srch.Contains("RUNC"))
            {
                Result.push_back(_DoubleEG_C);
                if (notify == kTRUE) cout << Procname[_DoubleEG_C] << "." << endl;
            }
            else if (srch.Contains("_D") || srch.Contains("RUND"))
            {
                Result.push_back(_DoubleEG_D);
                if (notify == kTRUE) cout << Procname[_DoubleEG_D] << "." << endl;
            }
            else if (srch.Contains("_E") || srch.Contains("RUNE"))
            {
                Result.push_back(_DoubleEG_E);
                if (notify == kTRUE) cout << Procname[_DoubleEG_E] << "." << endl;
            }
            else if (srch.Contains("_F") || srch.Contains("RUNF"))
            {
                Result.push_back(_DoubleEG_F);
                if (notify == kTRUE) cout << Procname[_DoubleEG_F] << "." << endl;
            }
            else if (srch.Contains("_G") || srch.Contains("RUNG"))
            {
                Result.push_back(_DoubleEG_G);
                if (notify == kTRUE) cout << Procname[_DoubleEG_G] << "." << endl;
            }
            else if (srch.Contains("_H") || srch.Contains("RUNH"))
            {
                Result.push_back(_DoubleEG_H);
                if (notify == kTRUE) cout << Procname[_DoubleEG_H] << "." << endl;
            }
            else
            {
                Result.push_back(_DoubleEG_Full);
                if (notify == kTRUE) cout << Procname[_DoubleEG_Full] << "." << endl;
            }

        }// end of if (EG)

    }// end of if(Double)

    else if (srch.Contains("SINGLE"))
    {
        if (srch.Contains("MU"))
        {
            // Checking if search contains intervals
            if (srch.Contains("BTO"))
                first = _SingleMuon_B;
            else if (srch.Contains("CTO"))
                first = _SingleMuon_C;
            else if (srch.Contains("DTO"))
                first = _SingleMuon_D;
            else if (srch.Contains("ETO"))
                first = _SingleMuon_E;
            else if (srch.Contains("FTO"))
                first = _SingleMuon_F;
            else if (srch.Contains("GTO"))
                first = _SingleMuon_G;
            else if (srch.Contains("HTO"))
                first = _SingleMuon_H;

            if (srch.Contains("TOB"))
                last = _SingleMuon_B;
            else if (srch.Contains("TOC"))
                last = _SingleMuon_C;
            else if (srch.Contains("TOD"))
                last = _SingleMuon_D;
            else if (srch.Contains("TOE"))
                last = _SingleMuon_E;
            else if (srch.Contains("TOF"))
                last = _SingleMuon_F;
            else if (srch.Contains("TOG"))
                last = _SingleMuon_G;
            else if (srch.Contains("TOH"))
                last = _SingleMuon_H;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if (first == _SingleMuon_B && last == _SingleMuon_H)
            {
                Result.push_back(_SingleMuon_Full);
                if (notify == kTRUE) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if (srch.Contains("FULL"))
            {
                Result.push_back(_SingleMuon_Full);
                if (notify == kTRUE) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
            else if (srch.Contains("_B")  || srch.Contains("RUNB"))
            {
                Result.push_back(_SingleMuon_B);
                if (notify == kTRUE) cout << Procname[_SingleMuon_B] << "." << endl;
            }
            else if (srch.Contains("_C") || srch.Contains("RUNC"))
            {
                Result.push_back(_SingleMuon_C);
                if (notify == kTRUE) cout << Procname[_SingleMuon_C] << "." << endl;
            }
            else if (srch.Contains("_D") || srch.Contains("RUND"))
            {
                Result.push_back(_SingleMuon_D);
                if (notify == kTRUE) cout << Procname[_SingleMuon_D] << "." << endl;
            }
            else if (srch.Contains("_E") || srch.Contains("RUNE"))
            {
                Result.push_back(_SingleMuon_E);
                if (notify == kTRUE) cout << Procname[_SingleMuon_E] << "." << endl;
            }
            else if (srch.Contains("_F") || srch.Contains("RUNF"))
            {
                Result.push_back(_SingleMuon_F);
                if (notify == kTRUE) cout << Procname[_SingleMuon_F] << "." << endl;
            }
            else if (srch.Contains("_G") || srch.Contains("RUNG"))
            {
                Result.push_back(_SingleMuon_G);
                if (notify == kTRUE) cout << Procname[_SingleMuon_G] << "." << endl;
            }
            else if (srch.Contains("_H") || srch.Contains("RUNH"))
            {
                Result.push_back(_SingleMuon_H);
                if (notify == kTRUE) cout << Procname[_SingleMuon_H] << "." << endl;
            }
            else
            {
                Result.push_back(_SingleMuon_Full);
                if (notify == kTRUE) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
        }// end of if(SingleMuon)

        if (srch.Contains("ELEC"))
        {
            // Checking if search contains intervals
            if (srch.Contains("BTO"))
                first = _SingleElectron_B;
            else if (srch.Contains("CTO"))
                first = _SingleElectron_C;
            else if (srch.Contains("DTO"))
                first = _SingleElectron_D;
            else if (srch.Contains("ETO"))
                first = _SingleElectron_E;
            else if (srch.Contains("FTO"))
                first = _SingleElectron_F;
            else if (srch.Contains("GTO"))
                first = _SingleElectron_G;
            else if (srch.Contains("HTO"))
                first = _SingleElectron_H;

            if (srch.Contains("TOB"))
                last = _SingleElectron_B;
            else if (srch.Contains("TOC"))
                last = _SingleElectron_C;
            else if (srch.Contains("TOD"))
                last = _SingleElectron_D;
            else if (srch.Contains("TOE"))
                last = _SingleElectron_E;
            else if (srch.Contains("TOF"))
                last = _SingleElectron_F;
            else if (srch.Contains("TOG"))
                last = _SingleElectron_G;
            else if (srch.Contains("TOH"))
                last = _SingleElectron_H;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if (first == _SingleElectron_B && last == _SingleElectron_H)
            {
                Result.push_back(_SingleElectron_Full);
                if (notify == kTRUE) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if (srch.Contains("FULL"))
            {
                Result.push_back(_SingleElectron_Full);
                if (notify == kTRUE) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
            else if (srch.Contains("_B")  || srch.Contains("RUNB"))
            {
                Result.push_back(_SingleElectron_B);
                if (notify == kTRUE) cout << Procname[_SingleElectron_B] << "." << endl;
            }
            else if (srch.Contains("_C") || srch.Contains("RUNC"))
            {
                Result.push_back(_SingleElectron_C);
                if (notify == kTRUE) cout << Procname[_SingleElectron_C] << "." << endl;
            }
            else if (srch.Contains("_D") || srch.Contains("RUND"))
            {
                Result.push_back(_SingleElectron_D);
                if (notify == kTRUE) cout << Procname[_SingleElectron_D] << "." << endl;
            }
            else if (srch.Contains("_E") || srch.Contains("RUNE"))
            {
                Result.push_back(_SingleElectron_E);
                if (notify == kTRUE) cout << Procname[_SingleElectron_E] << "." << endl;
            }
            else if (srch.Contains("_F") || srch.Contains("RUNF"))
            {
                Result.push_back(_SingleElectron_F);
                if (notify == kTRUE) cout << Procname[_SingleElectron_F] << "." << endl;
            }
            else if (srch.Contains("_G")|| srch.Contains("RUNG"))
            {
                Result.push_back(_SingleElectron_G);
                if (notify == kTRUE) cout << Procname[_SingleElectron_G] << "." << endl;
            }
            else if (srch.Contains("_H")|| srch.Contains("RUNH"))
            {
                Result.push_back(_SingleElectron_H);
                if (notify == kTRUE) cout << Procname[_SingleElectron_H] << "." << endl;
            }
            else
            {
                Result.push_back(_SingleElectron_Full);
                if (notify == kTRUE) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
        }// end of if(SingleElectron)

        if (srch.Contains("PHOT") || srch.Contains("GAMMA"))
        {
            // Checking if search contains intervals
            if (srch.Contains("BTO"))
                first = _SinglePhoton_B;
            else if (srch.Contains("CTO"))
                first = _SinglePhoton_C;
            else if (srch.Contains("DTO"))
                first = _SinglePhoton_D;
            else if (srch.Contains("ETO"))
                first = _SinglePhoton_E;
            else if (srch.Contains("FTO"))
                first = _SinglePhoton_F;
            else if (srch.Contains("GTO"))
                first = _SinglePhoton_G;
            else if (srch.Contains("HTO"))
                first = _SinglePhoton_H;

            if (srch.Contains("TOB"))
                last = _SinglePhoton_B;
            else if (srch.Contains("TOC"))
                last = _SinglePhoton_C;
            else if (srch.Contains("TOD"))
                last = _SinglePhoton_D;
            else if (srch.Contains("TOE"))
                last = _SinglePhoton_E;
            else if (srch.Contains("TOF"))
                last = _SinglePhoton_F;
            else if (srch.Contains("TOG"))
                last = _SinglePhoton_G;
            else if (srch.Contains("TOH"))
                last = _SinglePhoton_H;

            // Swapping first with last if necessary
            if (int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if (first == _SinglePhoton_B && last == _SinglePhoton_H)
            {
                Result.push_back(_SinglePhoton_Full);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_Full] << "." << endl;
            }
            else if (first != _None && last != _None)
            {
                for (Process_t pr=first; pr<=last; pr=next(pr))
                {
                    Result.push_back(pr);
                    if (notify == kTRUE)
                    {
                        if (pr != last) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if (srch.Contains("FULL"))
            {
                Result.push_back(_SinglePhoton_Full);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_Full] << "." << endl;
            }
            else if (srch.Contains("_B")  || srch.Contains("RUNB"))
            {
                Result.push_back(_SinglePhoton_B);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_B] << "." << endl;
            }
            else if (srch.Contains("_C") || srch.Contains("RUNC"))
            {
                Result.push_back(_SinglePhoton_C);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_C] << "." << endl;
            }
            else if (srch.Contains("_D") || srch.Contains("RUND"))
            {
                Result.push_back(_SinglePhoton_D);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_D] << "." << endl;
            }
            else if (srch.Contains("_E") || srch.Contains("RUNE"))
            {
                Result.push_back(_SinglePhoton_E);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_E] << "." << endl;
            }
            else if (srch.Contains("_F") || srch.Contains("RUNF"))
            {
                Result.push_back(_SinglePhoton_F);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_F] << "." << endl;
            }
            else if (srch.Contains("_G")|| srch.Contains("RUNG"))
            {
                Result.push_back(_SinglePhoton_G);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_G] << "." << endl;
            }
            else if (srch.Contains("_H")|| srch.Contains("RUNH"))
            {
                Result.push_back(_SinglePhoton_H);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_H] << "." << endl;
            }
            else
            {
                Result.push_back(_SinglePhoton_Full);
                if (notify == kTRUE) cout << Procname[_SinglePhoton_Full] << "." << endl;
            }
        }// end of if(SinglePhoton)

    }// end of if(Single)

    else if (srch.Contains("SIGNAL"))
    {
        Result.push_back(_DY_Full);
        if (notify == kTRUE) cout << Procname[_DY_Full] << "." << endl;
    }

    else if (srch.Contains("BKG"))
    {
        Result.push_back(_bkg_Full);
        if (notify == kTRUE)  cout << Procname[_bkg_Full] << "." << endl;
    }

    else if (srch.Contains("TEST"))
    {
        if (srch.Contains("EE"))
        {
            Result.push_back(_Test_EE);
            if (notify == kTRUE)  cout << Procname[_Test_EE] << "." << endl;
        }
        else if (srch.Contains("MUMU"))
        {
            Result.push_back(_Test_MuMu);
            if (notify == kTRUE)  cout << Procname[_Test_MuMu] << "." << endl;
        }
        else if (srch.Contains("EMU"))
        {
            Result.push_back(_Test_EMu);
            if (notify == kTRUE)  cout << Procname[_Test_EMu] << "." << endl;
        }
    }

    else
    {
        Result.push_back(_None);
        if (notify == kTRUE) cout << Procname[_None] << "." << endl;
    }

    if (instaGet == kTRUE)    // Applying SetProc if asked
    {
        cout << "Getting information for returned processes..";
        for (Int_t i=0; i<(int(Result.size())); i++)
        {
            if (Result[i] != _None) this->SetProc(Result[i], kFALSE);
        }
        cout << " Finished." << endl;
    }

    return Result;

}// end of FindProc()


void FileMgr::PrepareProcNames ()
{
    Procname[_None] = "None";
    Procname[_DY_10to50] = "DY_10to50";
    Procname[_DY_50to100] = "DY_50to100";
    Procname[_DY_100to200] = "DY_100to200";
    Procname[_DY_200to400] = "DY_200to400";
    Procname[_DY_400to500] = "DY_400to500";
    Procname[_DY_500to700] = "DY_500to700";
    Procname[_DY_700to800] = "DY_700to800";
    Procname[_DY_800to1000] = "DY_800to1000";
    Procname[_DY_1000to1500] = "DY_1000to1500";
    Procname[_DY_1500to2000] = "DY_1500to2000";
    Procname[_DY_2000to3000] = "DY_2000to3000";
    Procname[_EndOf_DY_Normal] = "EndOf_DY_Normal";
    Procname[_DYMuMu_10to50] = "DYMuMu_10to50";
    Procname[_DYMuMu_50to100] = "DYMuMu_50to100";
    Procname[_DYMuMu_100to200] = "DYMuMu_100to200";
    Procname[_DYMuMu_200to400] = "DYMuMu_200to400";
    Procname[_DYMuMu_400to500] = "DYMuMu_400to500";
    Procname[_DYMuMu_500to700] = "DYMuMu_500to700";
    Procname[_DYMuMu_700to800] = "DYMuMu_700to800";
    Procname[_DYMuMu_800to1000] = "DYMuMu_800to1000";
    Procname[_DYMuMu_1000to1500] = "DYMuMu_1000to1500";
    Procname[_DYMuMu_1500to2000] = "DYMuMu_1500to2000";
    Procname[_DYMuMu_2000to3000] = "DYMuMu_2000to3000";
    Procname[_EndOf_DYMuMu_Normal] = "EndOf_DYMuMu_Normal";
    Procname[_DYEE_10to50] = "DYEE_10to50";
    Procname[_DYEE_50to100] = "DYEE_50to100";
    Procname[_DYEE_100to200] = "DYEE_100to200";
    Procname[_DYEE_200to400] = "DYEE_200to400";
    Procname[_DYEE_400to500] = "DYEE_400to500";
    Procname[_DYEE_500to700] = "DYEE_500to700";
    Procname[_DYEE_700to800] = "DYEE_700to800";
    Procname[_DYEE_800to1000] = "DYEE_800to1000";
    Procname[_DYEE_1000to1500] = "DYEE_1000to1500";
    Procname[_DYEE_1500to2000] = "DYEE_1500to2000";
    Procname[_DYEE_2000to3000] = "DYEE_2000to3000";
    Procname[_EndOf_DYEE_Normal] = "EndOf_DYEE_Normal";
    Procname[_EndOf_MCsignal_Normal] = "EndOf_MCsignal_Normal";
    Procname[_DYTauTau_10to50] = "DYTauTau_10to50";
    Procname[_DYTauTau_50toInf] = "DYTauTau_50toInf";
    Procname[_EndOf_DYTauTau_Normal] = "EndOf_DYTauTau_Normal";
    Procname[_ttbar] = "ttbar";
    Procname[_ttbar_700to1000] = "ttbar_700to1000";
    Procname[_ttbar_1000toInf] = "ttbar_1000toInf";
    Procname[_EndOf_ttbar_Normal] = "EndOf_ttbar_Normal";
    Procname[_tW] = "tW";
    Procname[_tbarW] = "tbarW";
    Procname[_ZZ] = "ZZ";
    Procname[_WZ] = "WZ";
    Procname[_WW] = "WW";
    Procname[_EndOf_VVnST_Normal] = "EndOf_VVnST_Normal";
    Procname[_WJets] = "WJets";
    Procname[_WJets_ext2v5] = "WJets_ext2v5";
    Procname[_EndOf_WJets_Normal] = "EndOf_WJets";
    Procname[_QCDMuEnriched_15to20] = "QCDMuEnriched_15to20";
    Procname[_QCDMuEnriched_20to30] = "QCDMuEnriched_20to30";
    Procname[_QCDMuEnriched_30to50] = "QCDMuEnriched_30to50";
    Procname[_QCDMuEnriched_50to80] = "QCDMuEnriched_50to80";
    Procname[_QCDMuEnriched_80to120] = "QCDMuEnriched_80to120";
    Procname[_QCDMuEnriched_120to170] = "QCDMuEnriched_120to170";
    Procname[_QCDMuEnriched_170to300] = "QCDMuEnriched_170to300";
    Procname[_QCDMuEnriched_300to470] = "QCDMuEnriched_300to470";
    Procname[_QCDMuEnriched_470to600] = "QCDMuEnriched_470to600";
    Procname[_QCDMuEnriched_600to800] = "QCDMuEnriched_600to800";
    Procname[_QCDMuEnriched_800to1000] = "QCDMuEnriched_800to1000";
    Procname[_QCDMuEnriched_1000toInf] = "QCDMuEnriched_1000toInf";
    Procname[_EndOf_QCDMuEnriched_Normal] = "EndOf_QCDMuEnriched_Normal";
    Procname[_QCDEMEnriched_20to30] = "QCDEMEnriched_20to30";
    Procname[_QCDEMEnriched_30to50] = "QCDEMEnriched_30to50";
    Procname[_QCDEMEnriched_50to80] = "QCDEMEnriched_50to80";
    Procname[_QCDEMEnriched_80to120] = "QCDEMEnriched_80to120";
    Procname[_QCDEMEnriched_120to170] = "QCDEMEnriched_120to170";
    Procname[_QCDEMEnriched_170to300] = "QCDEMEnriched_170to300";
    Procname[_QCDEMEnriched_300toInf] = "QCDEMEnriched_300toInf";
    Procname[_EndOf_QCDEMEnriched_Normal] = "EndOf_QCDEMEnriched_Normal";
    Procname[_GJets_20to100] = "GJets_20to100";
    Procname[_GJets_100to200] = "GJets_100to200";
    Procname[_GJets_200to500] = "GJets_200to500";
    Procname[_GJets_500to1000] = "GJets_500to1000";
    Procname[_GJets_1000to2000] = "GJets_1000to2000";
    Procname[_GJets_2000to5000] = "GJets_2000to5000";
    Procname[_EndOf_GJets_Normal] = "EndOf_GJets_Normal";
    Procname[_EndOf_MCbkg_Normal] = "EndOf_MCbkg_Normal";
    Procname[_DoubleEG_B] = "DoubleEG_B";
    Procname[_DoubleEG_C] = "DoubleEG_C";
    Procname[_DoubleEG_D] = "DoubleEG_D";
    Procname[_DoubleEG_E] = "DoubleEG_E";
    Procname[_DoubleEG_F] = "DoubleEG_F";
    Procname[_DoubleEG_G] = "DoubleEG_G";
    Procname[_DoubleEG_H] = "DoubleEG_H";
    Procname[_EndOf_DoubleEG_Normal] = "EndOf_DoubleEG_Normal";
    Procname[_SingleMuon_B] = "SingleMuon_B";
    Procname[_SingleMuon_C] = "SingleMuon_C";
    Procname[_SingleMuon_D] = "SingleMuon_D";
    Procname[_SingleMuon_E] = "SingleMuon_E";
    Procname[_SingleMuon_F] = "SingleMuon_F";
    Procname[_SingleMuon_G] = "SingleMuon_G";
    Procname[_SingleMuon_H] = "SingleMuon_H";
    Procname[_EndOf_SingleMuon_Normal] = "EndOf_SinglMuon_Normal";
    Procname[_SingleElectron_B] = "SingleElectron_B";
    Procname[_SingleElectron_C] = "SingleElectron_C";
    Procname[_SingleElectron_D] = "SingleElectron_D";
    Procname[_SingleElectron_E] = "SingleElectron_E";
    Procname[_SingleElectron_F] = "SingleElectron_F";
    Procname[_SingleElectron_G] = "SingleElectron_G";
    Procname[_SingleElectron_H] = "SingleElectron_H";
    Procname[_SinglePhoton_B] = "SinglePhoton_B";
    Procname[_SinglePhoton_C] = "SinglePhoton_C";
    Procname[_SinglePhoton_D] = "SinglePhoton_D";
    Procname[_SinglePhoton_E] = "SinglePhoton_E";
    Procname[_SinglePhoton_F] = "SinglePhoton_F";
    Procname[_SinglePhoton_G] = "SinglePhoton_G";
    Procname[_SinglePhoton_H] = "SinglePhoton_H";
    Procname[_EndOf_SingleElectron_Normal] = "EndOf_SingleElectron_Normal";
    Procname[_EndOf_Data_Normal] = "EndOf_Data_Normal";
    Procname[_DY_Full] = "DY_Full";
    Procname[_DYMuMu_Full] = "DYMuMu_Full";
    Procname[_DYEE_Full] = "DYEE_Full";
    Procname[_EndOf_MCsignal_Special] = "EndOf_MCsignal_Special";
    Procname[_DYTauTau_Full] = "DYTauTau_Full";
    Procname[_ttbar_Full] = "ttbar_Full";
    Procname[_VVnST] = "VVnST";
    Procname[_WJets_Full] = "WJets_Full";
    Procname[_QCDMuEnriched_Full] = "QCDMuEnriched_Full";
    Procname[_QCDEMEnriched_Full] = "QCDEMEnriched_Full";
    Procname[_GJets_Full] = "GJets_Full";
    Procname[_bkg_Full] = "bkg_Full";
    Procname[_EndOf_MCbkg_Special] = "EndOf_MCbkg_Special";
    Procname[_DoubleEG_Full] = "DoubleEG_Full";
    Procname[_SingleMuon_Full] = "SingleMuon_Full";
    Procname[_SingleElectron_Full] = "SingleElectron_Full";
    Procname[_SinglePhoton_Full] = "SinglePhoton_Full";
    Procname[_EndOf_Data_Special] = "EndOf_Data_Special";
    Procname[_Test_MuMu] = "Test_MuMu";
    Procname[_Test_EE] = "Test_EE";
    Procname[_Test_EMu] = "Test_EMu";
    Procname[_EndOf_Test] = "EndOf_Test";
    Procname[_A_DY_50to100] = "Alt_DY50to100";
    Procname[_A_DY_100to250] = "Alt_DY100to250";
    Procname[_A_DY_250to400] = "Alt_DY250to400";
    Procname[_A_DY_400to650] = "Alt_DY400to650";
    Procname[_A_DY_650toInf] = "Alt_DY650toInf";
    Procname[_EndOf_A_DY_Normal] = "EndOf_Alt_DY_Normal";
    Procname[_A_WJets] = "Alt_WJets";
    Procname[_A_ZZ] = "Alt_ZZ";
    Procname[_A_WZ] = "Alt_WZ";
    Procname[_A_WW] = "Alt_WW";
    Procname[_EndOf_A_MCbkg_Normal] = "EndOf_Alt_MCbkg_Normal";
    Procname[_A_DY_Full] = "Alt_DY_Full";
    Procname[_EndOf_Alternatives] = "EndOf_Alternatives";
    return;

}// end of PrepareProcNames()


void FileMgr::CheckProcesses()
{
    Bool_t allOk = kTRUE;
    cout << "Checking processes: " << endl;
    for (Process_t pr=_DY_10to50; pr<_EndOf_Alternatives; pr=next(pr))
    {
        if (pr > _EndOf_Data_Special && pr <_EndOf_Test) continue;
        this->SetProc(pr, kTRUE);
        if (!Procname[pr].Length())
        {
            cout << "Process " << pr << ": no Procname found (providing its int expression)." << endl;
            allOk = kFALSE;
        }
        if (!CurrentProc)
        {
            cout << "Process " << Procname[pr] << ": no CurrentProc found." << endl;
            allOk = kFALSE;
        }
        if (!Tag.size())
        {
            cout << "Process " << Procname[pr] << ": no Tags found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(Tag.size())); i++)
            {
                if (!Tag[i].Length())
                {
                    cout << "Process " << Procname[pr] << ": no Tag[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if (!FileLocation.size())
        {
            cout << "Process " << Procname[pr] << ": no FileLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(FileLocation.size())); i++)
            {
                if (!FileLocation[i].Length())
                {
                    cout << "Process " << Procname[pr] << ": no FileLocation[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
                else if (pr<_EndOf_Data_Special && !FileLocation[i].Contains("/"))
                {
                    cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not contain '/'." << endl;
                    allOk = kFALSE;
                }
                /*else if (pr>_EndOf_Data_Special && !FileLocation[i].Contains(".root"))
                {
                    cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not contain '.root'." << endl;
                    allOk = kFALSE;
                }*/
            }

        }
        if (!FullLocation.size())
        {
            cout << "Process " << Procname[pr] << ": no FullLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(FullLocation.size())); i++)
            {
                if (!FullLocation[i].Length())
                {
                    cout << "Process " << Procname[pr] << ": no FullLocation[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
                else
                {
                    if (pr<_EndOf_Data_Special && !FullLocation[i].Contains("/"))
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not contain '/'" << endl;
                        allOk = kFALSE;
                    }
                    /*if (pr>_EndOf_Data_Special && !FullLocation[i].Contains(".root"))
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not contain '.root'" << endl;
                        allOk = kFALSE;
                    }
                    if (FullLocation[i][0] != '/')
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not begin with '/'" << endl;
                        allOk = kFALSE;
                    }*/
                }
            }

        }
        if (!TreeName.size())
        {
            cout << "Process " << Procname[pr] << ": no TreeName found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(TreeName.size())); i++)
            {
                if (!TreeName[i].Length())
                {
                    cout << "Process " << Procname[pr] << ": no TreeName[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if (!BaseLocation.Length())
        {
            cout << "Process " << Procname[pr] << ": no BaseLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if (BaseLocation[BaseLocation.Length()-1] != '/') {
                cout << "Process " << Procname[pr] << ": BaseLocation does not end with '/'" << endl;
                allOk = kFALSE;
            }
            /*if (BaseLocation[0] != '/') {
                cout << "Process " << Procname[pr] << ": BaseLocation does not begin with '/'" << endl;
                allOk = kFALSE;
            }*/

        }
        if (!nEvents.size())
        {
            cout << "Process " << Procname[pr] << ": no nEvents found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(nEvents.size())); i++)
            {
                if (!nEvents[i])
                {
                    cout << "Process " << Procname[pr] << ": no nEvents[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if (pr < _DoubleEG_B || (pr > _EndOf_Data_Normal && pr < _EndOf_MCbkg_Special) || pr > _EndOf_Data_Special)
        {
            if (isMC == kFALSE)
            {
                cout << "Process " << Procname[pr] << ": is said to be NOT MC." << endl;
                allOk = kFALSE;
            }
            if (!Xsec.size())
            {
                cout << "Process " << Procname[pr] << ": no Xsec found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<(int(Xsec.size())); i++)
                {
                    if (!Xsec[i])
                    {
                        cout << "Process " << Procname[pr] << ": no Xsec[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
            if (!Wsum.size())
            {
                cout << "Process " << Procname[pr] << ": no Wsum found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<(int(Wsum.size())); i++)
                {
                    if (!Wsum[i])
                    {
                        cout << "Process " << Procname[pr] << ": no Wsum[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
            if (Tag.size() != FileLocation.size() || FileLocation.size() != FullLocation.size() || FullLocation.size() != TreeName.size() ||
                 TreeName.size() != Xsec.size() || Xsec.size() != Wsum.size() || Wsum.size() != nEvents.size())
            {
                cout << "Process " << Procname[pr] << " Vector sizes do not match." << endl;
                allOk = kFALSE;
            }
        }// end of if(MC)
        else
        {
            if (isMC == kTRUE)
            {
                cout << "Process " << Procname[pr] << ": is said to be MC." << endl;
                allOk = kFALSE;
            }
            if (Tag.size() != FileLocation.size() || FileLocation.size() != FullLocation.size() || FullLocation.size() != TreeName.size() ||
                 TreeName.size() != nEvents.size())
            {
                cout << "Process " << Procname[pr] << " Vector sizes do not match." << endl;
                allOk = kFALSE;
            }
        }// end of else()
        if (!Type.Length())
        {
            cout << "Process " << Procname[pr] << ": no Type found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if (pr < _EndOf_MCsignal_Normal || (pr > _EndOf_Data_Normal && pr < _EndOf_MCsignal_Special))
            {
                if (Type != "SIGNAL")
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT SIGNAL." << endl;
                    allOk = kFALSE;
                }
            }
            else if ((pr > _EndOf_MCsignal_Normal && pr < _EndOf_MCbkg_Normal) || (pr > _EndOf_MCsignal_Special && pr < _EndOf_MCbkg_Special))
            {
                if (Type != "BKG")
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT BKG." << endl;
                    allOk = kFALSE;
                }
            }
            else if ((pr > _EndOf_MCbkg_Normal && pr < _EndOf_Data_Normal) || (pr > _EndOf_MCbkg_Special && pr < _EndOf_Data_Special))
            {
                if (Type != "DATA")
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT DATA." << endl;
                    allOk = kFALSE;
                }
            }
            else if (pr > _EndOf_Data_Special && pr < _EndOf_Test)
            {
                if (Type != "TEST")
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT TEST." << endl;
                    allOk = kFALSE;
                }
            }
        }       
        vector<Process_t> forChecking = FindProc(Procname[pr], kFALSE, kFALSE);
        if (!forChecking.size())
        {
            cout << "Process " << Procname[pr] << " No output recorded from FindProc()." << endl;
            allOk = kFALSE;
        }
        else if (forChecking[0] != pr)
        {
            cout << "Process " << Procname[pr] << ": FindProc() did not find this process." << endl;
            allOk = kFALSE;
        }
        this->ClearProc();
        if (CurrentProc != _None)
        {
            cout << "Process " << Procname[pr] << ": Current proc is not _None after ClearProc()." << endl;
            allOk = kFALSE;
        }
        if (Tag.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Tag." << endl;
            allOk = kFALSE;
        }
        if (Xsec.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Xsec." << endl;
            allOk = kFALSE;
        }
        if (Wsum.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Wsum." << endl;
            allOk = kFALSE;
        }
        if (nEvents.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the nEvents." << endl;
            allOk = kFALSE;
        }
        if (BaseLocation.Length())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the BaseLocation." << endl;
            allOk = kFALSE;
        }
        if (FileLocation.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Tag." << endl;
            allOk = kFALSE;
        }
        if (FullLocation.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the FullLocation." << endl;
            allOk = kFALSE;
        }
        if (TreeName.size())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the TreeName." << endl;
            allOk = kFALSE;
        }
        if (Type.Length())
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Type." << endl;
            allOk = kFALSE;
        }
    }// end of for()
    if (allOk == kTRUE) cout << "All OK." << endl;
    else cout << "Problems were detected." << endl;

}// end of CheckProcesses()


void FileMgr::SetupChain(Int_t i_tuple, TChain *chain)
{
    if (!CurrentProc)
    {
        cout << "No process set!" << endl;
        return;
    }
    for (Int_t i_tup=0; i_tup<((int)(nNtuples.size())); i_tup++)
    {
        if (i_tup != i_tuple && i_tuple != -1) continue;
        for (Int_t i=1; i<=nNtuples[i_tup]; i++)
        {
            if ((CurrentProc == _DY_50to100 || CurrentProc == _DYEE_50to100 || CurrentProc == _DYMuMu_50to100 || CurrentProc == _DYTauTau_50toInf ||
                ((CurrentProc == _DY_Full || CurrentProc == _DYEE_Full || CurrentProc == _DYMuMu_Full || CurrentProc == _DYTauTau_Full) && i_tuple == 3)) &&
                (i == 6 || i == 270 || i == 278 || i == 307 || i == 308 || i == 318 || i == 663 || i == 968)) continue; // Skipping several ntuple_skim_sig_X.root because they have some errors

            stringstream ss;
            ss << i;
            chain->Add(FullLocation[i_tup]+NtupleName[i_tup]+ss.str()+".root");
        }
    }
} // end of SetupChain()
