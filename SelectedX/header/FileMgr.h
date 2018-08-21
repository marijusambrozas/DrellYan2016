// Header for file management
// Created by Marijus Ambrozas 2018.08.14
// 2018.08.16: Changed FindProc to return vector of processes (so you can type a larger invM interval and get more processes), changed searching mechanisms.
// 2018.08.17: Added TString Type  ("DATA", "SIGNAL", or "BKG") and vector<TString> TreeName.
// 2018.08.20: Added Test processes (for testing in local PC). Moved CheckProcesses and PrepareProcNames to private functions. CheckProcesses now checks if ClearProc works.

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
    _WJets, _EndOf_WJets,
    _QCDMuEnriched_15to20, _QCDMuEnriched_20to30, _QCDMuEnriched_30to50, _QCDMuEnriched_50to80, _QCDMuEnriched_80to120, _QCDMuEnriched_120to170,
    _QCDMuEnriched_170to300, _QCDMuEnriched_300to470, _QCDMuEnriched_470to600, _QCDMuEnriched_600to800, _QCDMuEnriched_800to1000, _QCDMuEnriched_1000toInf,
    _EndOf_QCDMuEnriched_Normal,
    _QCDEMEnriched_20to30, _QCDEMEnriched_30to50, _QCDEMEnriched_50to80, _QCDEMEnriched_80to120, _QCDEMEnriched_120to170, _QCDEMEnriched_170to300,
    _QCDEMEnriched_300toInf, _EndOf_QCDEMEnriched_Normal,
    _EndOf_MCbkg_Normal,
    _DoubleEG_B, _DoubleEG_C, _DoubleEG_D, _DoubleEG_E, _DoubleEG_F, _DoubleEG_G, _DoubleEG_H,
    _EndOf_DoubleEG_Normal,
    _SingleMuon_B, _SingleMuon_C, _SingleMuon_D, _SingleMuon_E, _SingleMuon_F, _SingleMuon_G, _SingleMuon_H,
    _EndOf_SinglMuon_Normal,
    _SingleElectron_B, _SingleElectron_C, _SingleElectron_D, _SingleElectron_E, _SingleElectron_F, _SingleElectron_G, _SingleElectron_H,
    _EndOf_SingleElectron_Normal,
    _EndOf_Data_Normal,
    // "Special" processes - similar processes are combined
    _DY_Full, _DYMuMu_Full, _DYEE_Full,
    _EndOf_MCsignal_Special,
    _DYTauTau_Full, _ttbar_Full, _VVnST, _QCDMuEnriched_Full, _QCDEMEnriched_Full, _bkg_Full,
    _EndOf_MCbkg_Special, // there is no WJets in bkgSpecial
    _DoubleEG_Full, _SingleMuon_Full, _SingleElectron_Full,
    _EndOf_Data_Special,
    // Processes for testing at local pc
    _Test_MuMu, _Test_EE, _Test_EMu,
    _EndOf_Test
};

inline
Process_t next ( Process_t pr )    // Processes that begin with "EndOf" will be skipped by this
{
  if ( pr == _EndOf_Test )
      return pr;
  else if ( pr == _DYEE_2000to3000 || pr == _QCDEMEnriched_300toInf || pr == _SingleElectron_H )
      return Process_t(int(pr)+3);
  else if ( pr == _DY_2000to3000 || pr == _DYMuMu_2000to3000 || pr == _EndOf_DYEE_Normal || pr == _DYTauTau_50toInf || pr == _ttbar_1000toInf ||
            pr == _WW || pr == _WJets || pr == _QCDMuEnriched_1000toInf || pr == _EndOf_QCDEMEnriched_Normal || pr == _DoubleEG_H ||
            pr == _SingleMuon_H || pr == _EndOf_SingleElectron_Normal || pr == _DYEE_Full || pr == _bkg_Full || pr == _SingleElectron_Full )
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
        vector<TString> TreeName;
        vector<Double_t> Xsec;
        vector<Double_t> Wsum;
        vector<Double_t> nEvents;

        TString BaseLocation;
        TString Type;       
        Bool_t isMC;       

        map<Process_t, TString> Procname;

        // -- Constructor -- //
        FileMgr ( Process_t pr = _None );

        vector<Process_t> FindProc ( TString search, Bool_t notify = kTRUE, Bool_t instaGet = kFALSE );
        void NextProc ();
        void GetProc ( Process_t pr = _None, Bool_t ClearOld = kTRUE );
        void ClearProc ();        

private:
        Bool_t namesSet = kFALSE;
        Bool_t processesChecked = kFALSE;

        void PrepareProcNames ();
        void CheckProcesses ();

};// end of class definition


// ---------- Constructor ---------- //

FileMgr::FileMgr ( Process_t pr )
{
    if ( namesSet == kFALSE ) { this->PrepareProcNames(); namesSet = kTRUE; }
    if ( processesChecked == kFALSE ) { this->CheckProcesses(); processesChecked = kTRUE; }
    CurrentProc = pr;
    this->GetProc(CurrentProc, kTRUE);
}


// ----------- Functions ----------- //

void FileMgr::NextProc()
{
    CurrentProc = next(CurrentProc);
    this->GetProc(CurrentProc, kTRUE);
}


void FileMgr::ClearProc()
{
    if ( CurrentProc != _None )
    {
        CurrentProc = _None;
        BaseLocation = "";
        Type = "";
        isMC = kFALSE;
        this->GetProc(CurrentProc, kTRUE);
    }
}


void FileMgr::GetProc ( Process_t pr, Bool_t ClearOld )
{
    if ( ClearOld == kTRUE )
    {
        Tag.clear();
        FullLocation.clear();
        FileLocation.clear();
        TreeName.clear();
        Xsec.clear();
        Wsum.clear();
        nEvents.clear();
    }
    TString Location;
    CurrentProc = pr;

    if ( pr == _DY_10to50 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M10to50_v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_v2" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_ext1v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr == _DY_50to100 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M50to100" ); Xsec.push_back( 5869.58346 ); Wsum.push_back( 81780984 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M100to200" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200_ext" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M200to400" ); Xsec.push_back( 7.67 ); Wsum.push_back( 169676 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M400to500" ); Xsec.push_back( 0.423 ); Wsum.push_back( 151190 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M500to700" ); Xsec.push_back( 0.24 ); Wsum.push_back( 144096 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M700to800" ); Xsec.push_back( 0.035 ); Wsum.push_back( 136892 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M800to1000" ); Xsec.push_back( 0.03 ); Wsum.push_back( 131586 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M1000to1500" ); Xsec.push_back( 0.016 ); Wsum.push_back( 120010 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M1500to2000" ); Xsec.push_back( 0.002 ); Wsum.push_back( 111709 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M2000to3000" ); Xsec.push_back( 0.00054 ); Wsum.push_back( 101298 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DY_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M10to50_v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_v2" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_ext1v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M50to100" ); Xsec.push_back( 5869.58346 ); Wsum.push_back( 81780984 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200_ext" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M200to400" ); Xsec.push_back( 7.67 ); Wsum.push_back( 169676 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M400to500" ); Xsec.push_back( 0.423 ); Wsum.push_back( 151190 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M500to700" ); Xsec.push_back( 0.24 ); Wsum.push_back( 144096 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M700to800" ); Xsec.push_back( 0.035 ); Wsum.push_back( 136892 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M800to1000" ); Xsec.push_back( 0.03 ); Wsum.push_back( 131586 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M1000to1500" ); Xsec.push_back( 0.016 ); Wsum.push_back( 120010 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M1500to2000" ); Xsec.push_back( 0.002 ); Wsum.push_back( 111709 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M2000to3000" ); Xsec.push_back( 0.00054 ); Wsum.push_back( 101298 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYMuMu_10to50 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr == _DYMuMu_50to100 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37176.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYMuMu_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37176.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_10to50 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 306508623 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_50to100 ) // Only EE evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 306508623 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _tW )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 23944342 ); // I get Wsum=137540054
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 177139200 );
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/*.root";        // There also is madgraph version
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_15to20 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 ); nEvents.push_back( 4141251 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_20to30 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 ); nEvents.push_back( 31302080 ); // One file did not have any trees in it
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_30to50 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 ); nEvents.push_back( 29717171 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_50to80 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 ); nEvents.push_back( 19806914 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_80to120 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 13555323 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 9797243 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_120to170 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 8042720 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 11938137 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_170to300 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 7947158 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 9403070 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 19607775 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_300to470 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 7937587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 16452587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 24605502 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_470to600 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 3851523 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 5663755 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_600to800 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 4010135 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 5971173 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 9756852 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_800to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 3962747 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 5838539 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 9966146 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 3861436 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 9609820 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 ); nEvents.push_back( 4141251 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 ); nEvents.push_back( 31302080 );
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 ); nEvents.push_back( 29717171 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 ); nEvents.push_back( 19806914 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 13555323 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 9797243 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 8042720 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 11938137 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 7947158 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 9403070 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 19607775 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 7937587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 16452587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 24605502 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 3851523 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 5663755 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 4010135 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 5971173 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 9756852 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 3962747 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 5838539 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 9966146 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 3861436 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 9609820  );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_20to30 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 9218952 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_30to50 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 4730195 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 6768384 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_50to80 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 22337068 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 23474168 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_80to120 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 35841780 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 41853502 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_120to170 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 35817276 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 41954033 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_170to300 )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162.0 ); nEvents.push_back( 11540162 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_300toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633.0 ); nEvents.push_back( 7373633 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 9218952 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 4730195 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 6768384 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 22337068 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 23474168 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 35841780 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 41853502 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 35817276 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 41954033 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162.0 ); nEvents.push_back( 11540162 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633.0 ); nEvents.push_back( 7373633 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 23944342 ); // I get Wsum=137540054
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 177139200 );
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/*.root";        // There also is madgraph version
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 ); nEvents.push_back( 4141251 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 ); nEvents.push_back( 31302080 );
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 ); nEvents.push_back( 29717171 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 ); nEvents.push_back( 19806914 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 13555323 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 ); nEvents.push_back( 9797243 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 8042720 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 ); nEvents.push_back( 11938137 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 7947158 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 9403070 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 ); nEvents.push_back( 19607775 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 7937587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 16452587 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 ); nEvents.push_back( 24605502 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 3851523 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 ); nEvents.push_back( 5663755 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 4010135 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 5971173 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 ); nEvents.push_back( 9756852 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 3962747 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 5838539 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 ); nEvents.push_back( 9966146 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 3861436 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 ); nEvents.push_back( 9609820 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 9218952 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 4730195 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 6768384 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 22337068 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 23474168 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 35841780 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 41853502 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 35817276 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 41954033 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162 ); nEvents.push_back( 11540162 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633 ); nEvents.push_back( 7373633 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_B_0000" ); nEvents.push_back( 103625724 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_B_0001" ); nEvents.push_back( 33031246 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_C" ); nEvents.push_back( 45521797 );
        Location = "DoubleEG/crab_DoubleEG_RunC/180326_143612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_D" ); nEvents.push_back( 52422569 );
        Location = "DoubleEG/crab_DoubleEG_RunD/180326_143654/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_E" ); nEvents.push_back( 47326656 );
        Location = "DoubleEG/crab_DoubleEG_RunE/180326_143750/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_F" ); nEvents.push_back( 33943052 );
        Location = "DoubleEG/crab_DoubleEG_RunF/180326_143846/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_G_0000" ); nEvents.push_back( 71864512 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0001" ); nEvents.push_back( 4669958 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_Hver2_0000" ); nEvents.push_back( 68821231 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0001" ); nEvents.push_back( 11645108 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver3" ); nEvents.push_back( 2021309 );
        Location = "DoubleEG/crab_DoubleEG_RunHver3/180326_144719/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_B_0000" ); nEvents.push_back( 103625724 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_B_0001" ); nEvents.push_back( 33031246 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_C" ); nEvents.push_back( 45521797 );
        Location = "DoubleEG/crab_DoubleEG_RunC/180326_143612/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_D" ); nEvents.push_back( 52422569 );
        Location = "DoubleEG/crab_DoubleEG_RunD/180326_143654/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_E" ); nEvents.push_back( 47326656 );
        Location = "DoubleEG/crab_DoubleEG_RunE/180326_143750/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_F" ); nEvents.push_back( 33943052 );
        Location = "DoubleEG/crab_DoubleEG_RunF/180326_143846/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0000" ); nEvents.push_back( 71864512 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0001" ); nEvents.push_back( 4669958 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0000" ); nEvents.push_back( 68821231 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0001" ); nEvents.push_back( 11645108 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0001/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver3" ); nEvents.push_back( 2021309 );
        Location = "DoubleEG/crab_DoubleEG_RunHver3/180326_144719/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_B" ); nEvents.push_back( 108561074 );
        Location = "SingleMuon/crab_SingleMuon_RunB/180326_143105/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_C" ); nEvents.push_back( 64715287 );
        Location = "SingleMuon/crab_SingleMuon_RunC/180326_143152/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_D" ); nEvents.push_back( 96652779 );
        Location = "SingleMuon/crab_SingleMuon_RunD/180326_143257/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_E" ); nEvents.push_back( 87358348 );
        Location = "SingleMuon/crab_SingleMuon_RunE/180326_143338/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_F" ); nEvents.push_back( 64986568 );
        Location = "SingleMuon/crab_SingleMuon_RunF/180326_143419/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_G" ); nEvents.push_back( 138710659 );
        Location = "SingleMuon/crab_SingleMuon_RunG/180326_144335/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_Hver2" ); nEvents.push_back( 141936183 );
        Location = "SingleMuon/crab_SingleMuon_RunHver2/180326_144412/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver3" ); nEvents.push_back( 4386928 );
        Location = "SingleMuon/crab_SingleMuon_RunHver3/180326_144454/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_B" ); nEvents.push_back( 108561074 );
        Location = "SingleMuon/crab_SingleMuon_RunB/180326_143105/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_C" ); nEvents.push_back( 64715287 );
        Location = "SingleMuon/crab_SingleMuon_RunC/180326_143152/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_D" ); nEvents.push_back( 96652779 );
        Location = "SingleMuon/crab_SingleMuon_RunD/180326_143257/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_E" ); nEvents.push_back( 87358348 );
        Location = "SingleMuon/crab_SingleMuon_RunE/180326_143338/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_F" ); nEvents.push_back( 64986568 );
        Location = "SingleMuon/crab_SingleMuon_RunF/180326_143419/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_G" ); nEvents.push_back( 138710659 );
        Location = "SingleMuon/crab_SingleMuon_RunG/180326_144335/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver2" ); nEvents.push_back( 141936183 );
        Location = "SingleMuon/crab_SingleMuon_RunHver2/180326_144412/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver3" ); nEvents.push_back( 4386928 );
        Location = "SingleMuon/crab_SingleMuon_RunHver3/180326_144454/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_B" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_C" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_D" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_E" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_G" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        BaseLocation = "/xrootd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_B" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_C" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_D" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_E" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_G" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/*.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _Test_MuMu )
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back( "ZToMuMu_M4500to6000_4"); Xsec.push_back( 1.0 ); Wsum.push_back( 10200.0 ); nEvents.push_back( 10200 );
        Location = "ZToMuMu_M4500to6000_4.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _Test_EE )
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back( "ZToEE_M4500to6000_2"); Xsec.push_back( 1.0 ); Wsum.push_back( 39200.0 ); nEvents.push_back( 39200 );
        Location = "ZToEE_M4500to6000_2.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _Test_EMu )
    {
        isMC = kTRUE;
        Type = "TEST";
        BaseLocation = "/media/sf_DATA/test/";

        Tag.push_back( "WW_34"); Xsec.push_back( 1.0 ); Wsum.push_back( 10788.0 ); nEvents.push_back( 10788 );
        Location = "WW_34.root";
        TreeName.push_back( "recoTree/DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
}// end of Get()


vector<Process_t> FileMgr::FindProc ( TString search, Bool_t notify, Bool_t instaGet )
{
    TString srch = search;
    if ( notify == kTRUE ) cout << "Searched for: " << search << "\nFound: ";
    vector<Process_t> Result;
    Process_t first = _None, last = _None;
    if ( srch.Contains("DY") || srch.Contains("dy") || srch.Contains("dY") || srch.Contains("Dy") || srch.Contains("DrellYan") || srch.Contains("DRELLYAN") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            // Checking for various intervals
            if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                 srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                 srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                 srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                 srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                 srch.Contains("Inf -") )
                first = _EndOf_DYMuMu_Normal;
            else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                      srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                      srch.Contains("3000_TO") )
                first = _EndOf_DYMuMu_Normal;
            else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                      srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                      srch.Contains("2000_TO") )
                first = _DYMuMu_2000to3000;
            else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                      srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                      srch.Contains("1500_TO") )
                first = _DYMuMu_1500to2000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _DYMuMu_1000to1500;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _DYMuMu_800to1000;
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                      srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                      srch.Contains("700_TO") )
                first = _DYMuMu_700to800;
            else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                      srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                      srch.Contains("500_TO") )
                first = _DYMuMu_500to700;
            else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                      srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                      srch.Contains("400_TO") )
                first = _DYMuMu_400to500;
            else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                      srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                      srch.Contains("200_TO") )
                first = _DYMuMu_200to400;
            else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                      srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                      srch.Contains("100_TO") )
                first = _DYMuMu_100to200;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                first = _DYMuMu_50to100;
            else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                      srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                first = _DYMuMu_10to50;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _DYMuMu_2000to3000;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _DYMuMu_2000to3000;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _DYMuMu_1500to2000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _DYMuMu_1000to1500;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _DYMuMu_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _DYMuMu_700to800;
            else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                      srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                      srch.Contains("TO_700") )
                last = _DYMuMu_500to700;
            else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                      srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                      srch.Contains("TO_500") )
                last = _DYMuMu_400to500;
            else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                      srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                      srch.Contains("TO_400") )
                last = _DYMuMu_200to400;
            else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                      srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                      srch.Contains("TO_200") )
                last = _DYMuMu_100to200;
            else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                      srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                      srch.Contains("TO_100") )
                last = _DYMuMu_50to100;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _DYMuMu_10to50;
            else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                      srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                      srch.Contains("TO_10") )
                last = _EndOf_DY_Normal;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None)
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _DYMuMu_10to50 && last == _DYMuMu_2000to3000 )
            {
                Result.push_back(_DYMuMu_Full);
                if ( notify == kTRUE ) cout << Procname[_DYMuMu_Full] << "." << endl;
            }
            else if ( first != _None && last != _None )
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_DYMuMu_Full);
                if ( notify == kTRUE ) cout << Procname[_DYMuMu_Full] << "." << endl;
            }

        }// end of if(DYMuMu)

        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
             srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            // Checking for various intervals
            if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                 srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                 srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                 srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                 srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                 srch.Contains("Inf -") )
                first = _EndOf_DYEE_Normal;
            else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                      srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                      srch.Contains("3000_TO") )
                first = _EndOf_DYEE_Normal;
            else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                      srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                      srch.Contains("2000_TO") )
                first = _DYEE_2000to3000;
            else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                      srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                      srch.Contains("1500_TO") )
                first = _DYEE_1500to2000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _DYEE_1000to1500;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _DYEE_800to1000;
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                      srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                      srch.Contains("700_TO") )
                first = _DYEE_700to800;
            else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                      srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                      srch.Contains("500_TO") )
                first = _DYEE_500to700;
            else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                      srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                      srch.Contains("400_TO") )
                first = _DYEE_400to500;
            else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                      srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                      srch.Contains("200_TO") )
                first = _DYEE_200to400;
            else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                      srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                      srch.Contains("100_TO") )
                first = _DYEE_100to200;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                first = _DYEE_50to100;
            else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                      srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                first = _DYEE_10to50;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _DYEE_2000to3000;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _DYEE_2000to3000;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _DYEE_1500to2000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _DYEE_1000to1500;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _DYEE_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _DYEE_700to800;
            else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                      srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                      srch.Contains("TO_700") )
                last = _DYEE_500to700;
            else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                      srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                      srch.Contains("TO_500") )
                last = _DYEE_400to500;
            else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                      srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                      srch.Contains("TO_400") )
                last = _DYEE_200to400;
            else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                      srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                      srch.Contains("TO_200") )
                last = _DYEE_100to200;
            else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                      srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                      srch.Contains("TO_100") )
                last = _DYEE_50to100;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _DYEE_10to50;
            else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                      srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                      srch.Contains("TO_10") )
                last = _EndOf_DYMuMu_Normal;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _DYEE_10to50 && last == _DYEE_2000to3000 )
            {
                Result.push_back(_DYEE_Full);
                if ( notify == kTRUE ) cout << Procname[_DYEE_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_DYEE_Full);
                if ( notify == kTRUE ) cout << Procname[_DYEE_Full] << "." << endl;
            }

        }// end of if(DYEE)

        else if ( srch.Contains("TauTau") || srch.Contains("tautau") || srch.Contains("TAUTAU") || srch.Contains("Ditau") || srch.Contains("DiTau") ||
             srch.Contains("ditau") || srch.Contains("diTau") || srch.Contains("DITAU") )
        {
            if ( srch.Contains("full") || srch.Contains("Full") || srch.Contains("FULL") )
            {
                Result.push_back(_DYTauTau_Full);
                if ( notify == kTRUE ) cout << Procname[_DYTauTau_Full] << "." << endl;
            }
            else
            {
                // Checking for various intervals
                if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                     srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                     srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                     srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                     srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                     srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                     srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                     srch.Contains("Inf -") )
                    first = _EndOf_DYTauTau_Normal;
                else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                          srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                          srch.Contains("3000_TO") )
                    first = _EndOf_DYTauTau_Normal;
                else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                          srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                          srch.Contains("2000_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                          srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                          srch.Contains("1500_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                          srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                          srch.Contains("1000_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                          srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                          srch.Contains("800_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                          srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                          srch.Contains("700_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                          srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                          srch.Contains("500_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                          srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                          srch.Contains("400_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                          srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                          srch.Contains("200_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                          srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                          srch.Contains("100_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                          srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                    first = _DYTauTau_50toInf;
                else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                          srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                    first = _DYTauTau_10to50;

                if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                     srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                     srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                     srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                     srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                     srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                     srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                     srch.Contains("- Inf") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                          srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                          srch.Contains("TO_3000") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                          srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                          srch.Contains("TO_2000") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                          srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                          srch.Contains("TO_1500") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                          srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                          srch.Contains("TO_1000") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                          srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                          srch.Contains("TO_800") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                          srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                          srch.Contains("TO_700") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                          srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                          srch.Contains("TO_500") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                          srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                          srch.Contains("TO_400") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                          srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                          srch.Contains("TO_200") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                          srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                          srch.Contains("TO_100") )
                    last = _DYTauTau_50toInf;
                else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                          srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                          srch.Contains("TO_50") )
                    last = _DYTauTau_10to50;
                else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                          srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                          srch.Contains("TO_10") )
                    last = _EndOf_MCsignal_Normal;

                // Swapping first with last if necessary
                if ( int(first)>int(last) && last!=_None )
                {
                    Process_t NewLast = Process_t(int(first)-1);
                    first = Process_t(int(last)+1);
                    last = NewLast;
                }
                if ( first == _DYTauTau_10to50 && last == _DYTauTau_50toInf )
                {
                    Result.push_back(_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_DYTauTau_Full] << "." << endl;
                }
                else if ( first != _None && last != _None )
                {
                    for ( Process_t pr=first; pr<=last; pr=next(pr) )
                    {
                        Result.push_back(pr);
                        if ( notify == kTRUE )
                        {
                            if ( pr != last ) cout << Procname[pr] << ", ";
                            else cout << Procname[pr] << "." << endl;
                        }
                    }
                }
                else
                {
                    Result.push_back(_DYEE_Full);
                    if ( notify == kTRUE ) cout << Procname[_DYEE_Full] << "." << endl;
                }
            }

        }// end of if(DYTauTau)
        else
        {
            // Checking for various intervals
            if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                 srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                 srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                 srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                 srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                 srch.Contains("Inf -") )
                first = _EndOf_DY_Normal;
            else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                      srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                      srch.Contains("3000_TO") )
                first = _EndOf_DY_Normal;
            else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                      srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                      srch.Contains("2000_TO") )
                first = _DY_2000to3000;
            else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                      srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                      srch.Contains("1500_TO") )
                first = _DY_1500to2000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _DY_1000to1500;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _DY_800to1000;
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                      srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                      srch.Contains("700_TO") )
                first = _DY_700to800;
            else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                      srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                      srch.Contains("500_TO") )
                first = _DY_500to700;
            else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                      srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                      srch.Contains("400_TO") )
                first = _DY_400to500;
            else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                      srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                      srch.Contains("200_TO") )
                first = _DY_200to400;
            else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                      srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                      srch.Contains("100_TO") )
                first = _DY_100to200;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                first = _DY_50to100;
            else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                      srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                first = _DY_10to50;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _DY_2000to3000;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _DY_2000to3000;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _DY_1500to2000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _DY_1000to1500;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _DY_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _DY_700to800;
            else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                      srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                      srch.Contains("TO_700") )
                last = _DY_500to700;
            else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                      srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                      srch.Contains("TO_500") )
                last = _DY_400to500;
            else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                      srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                      srch.Contains("TO_400") )
                last = _DY_200to400;
            else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                      srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                      srch.Contains("TO_200") )
                last = _DY_100to200;
            else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                      srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                      srch.Contains("TO_100") )
                last = _DY_50to100;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") || srch.Contains("TO 50") ||
                      srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") || srch.Contains("TO_50") )
                last = _DY_10to50;
            else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") || srch.Contains("TO 10") ||
                      srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") || srch.Contains("TO_10") )
                last = _None;

            //Swapping first with last if necessary
            if ( int(first) > int(last) )
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _DY_10to50 && last == _DY_2000to3000 )
            {
                Result.push_back(_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_DY_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_DY_Full] << "." << endl;
            }

        }// end of else (DY)

    }// end of if(DrellYan)

    else if ( srch.Contains("ttbar") || srch.Contains("tt") || srch.Contains("TTbar") || srch.Contains("TT") || srch.Contains("ttBAR") ||
         srch.Contains("TTBAR") || srch.Contains("TantiT") || srch.Contains("tantit") || srch.Contains("TANTIT") || srch.Contains("tANTIt") ||
         srch.Contains("TopAntiTop") )
    {
        if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
        {
            Result.push_back(_ttbar_Full);
            if ( notify == kTRUE ) cout << Procname[_ttbar_Full] << "." << endl;
        }
        // Checking for various intervals
        else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") || srch.Contains("700 TO") ||
             srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") || srch.Contains("700_TO") )
            first = _ttbar_700to1000;
        else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                  srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                  srch.Contains("1000_TO") )
            first = _ttbar_1000toInf;
        else if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                  srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                  srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                  srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                  srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                  srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                  srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                  srch.Contains("Inf -") )
            first = _EndOf_ttbar_Normal;
        else first = _ttbar;

        if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") || srch.Contains("TO 700") ||
             srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") || srch.Contains("TO_700") )
            last = _ttbar;
        else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                  srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                  srch.Contains("TO_1000") )
            last = _ttbar_700to1000;
        else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                  srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                  srch.Contains("TO_1500") )
            last = _ttbar_1000toInf;
        else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                  srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                  srch.Contains("TO_2000") )
            last = _ttbar_1000toInf;
        else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                  srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                  srch.Contains("TO_3000") )
            last = _ttbar_1000toInf;
        else if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                  srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                  srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                  srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                  srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                  srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                  srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                  srch.Contains("- Inf") )
            last = _ttbar_1000toInf;
        else last = _ttbar;

        // Swapping first with last if necessary
        if ( int(first)>int(last) && last!=_None )
        {
            Process_t NewLast = Process_t(int(first)-1);
            first = Process_t(int(last)+1);
            last = NewLast;
        }
        if ( first == _ttbar && last == _ttbar_1000toInf )
        {
            Result.push_back(_ttbar_Full);
            if ( notify == kTRUE ) cout << Procname[_ttbar_Full] << "." << endl;
        }
        else if ( first != _None && last != _None )
        {
            for ( Process_t pr=first; pr<=last; pr=next(pr) )
            {
                Result.push_back(pr);
                if ( notify == kTRUE )
                {
                    if ( pr != last ) cout << Procname[pr] << ", ";
                    else cout << Procname[pr] << "." << endl;
                }
            }
        }
        else
        {
            Result.push_back(_ttbar_Full);
            if ( notify == kTRUE ) cout << Procname[_ttbar_Full] << "." << endl;
        }

    }// end of if(ttbar)

    else if ( srch.Contains("tW") || srch.Contains("TW") || srch.Contains("SingleTop") || srch.Contains("singleTop") || srch.Contains("singletop") ||
         srch.Contains("SingleTOP") || srch.Contains("singleTOP") || srch.Contains("SINGLETOP") )
    {
        Result.push_back(_tW);
        if ( notify == kTRUE ) cout << Procname[_tW] << "." << endl;
    }

    else if ( srch.Contains("tbarW") || srch.Contains("TbarW") || srch.Contains("TBarW") || srch.Contains("tBarW") || srch.Contains("TBARW") ||
              srch.Contains("SingleAntiTop") || srch.Contains("SingleAntitop") || srch.Contains("singleAntiTop") || srch.Contains("singleAntitop") ||
              srch.Contains("SingleAntiTOP") || srch.Contains("singleAntiTOP") || srch.Contains("SingleANTITOP") || srch.Contains("singleANTITOP") ||
              srch.Contains("SINGLEANTITOP") )
    {
        Result.push_back(_tbarW);
        if ( notify == kTRUE ) cout << Procname[_tbarW] << "." << endl;
    }

    else if ( srch.Contains("ZZ") || srch.Contains("zz") )
    {
        Result.push_back(_ZZ);
        if ( notify == kTRUE ) cout << Procname[_ZZ] << "." << endl;
    }

    else if ( srch.Contains("WZ") || srch.Contains("wz") || srch.Contains("ZW") || srch.Contains("zw") )
    {
        Result.push_back(_WZ);
        if ( notify == kTRUE ) cout << Procname[_WZ] << "." << endl;
    }

    else if ( srch.Contains("WW") || srch.Contains("ww") )
    {
        Result.push_back(_WW);
        if ( notify == kTRUE ) cout << Procname[_WW] << "." << endl;
    }

    else if ( srch.Contains("DIBOSON") || srch.Contains("diboson") || srch.Contains("Diboson") || srch.Contains("diBoson"))
    {
        Result.push_back(_ZZ);
        Result.push_back(_WZ);
        Result.push_back(_WW);
        if ( notify == kTRUE ) cout << Procname[_ZZ] << ", " << Procname[_WZ] << ", " << Procname[_WW] << "." << endl;
    }

    else if ( srch.Contains("VVnST") || srch.Contains("VVNST") || srch.Contains("vvnst") || srch.Contains("VVnst") || srch.Contains("vvnST") ||
              srch.Contains("VV") || srch.Contains("vv") )
    {
        Result.push_back(_VVnST);
        if ( notify == kTRUE ) cout << Procname[_VVnST] << "." << endl;
    }

    else if ( srch.Contains("WJets") || srch.Contains("Wjets") || srch.Contains("wjets") || srch.Contains("WJETS") || srch.Contains("W+jets") ||
              srch.Contains("W+Jets") || srch.Contains("W+JETS") )
    {
        Result.push_back(_WJets);
        if ( notify == kTRUE ) cout << Procname[_WJets] << "." << endl;
    }

    else if ( srch.Contains("QCD") || srch.Contains("qcd") )
    {
        if ( srch.Contains("Mu") || srch.Contains("mu") || srch.Contains("MU") )
        {
            // Checking for various intervals
            if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                 srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                 srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                 srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                 srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                 srch.Contains("Inf -") )
                first = _EndOf_QCDMuEnriched_Normal;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _QCDMuEnriched_1000toInf;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _QCDMuEnriched_800to1000;
            else if ( srch.Contains("600to") || srch.Contains("600To") || srch.Contains("600TO") || srch.Contains("600 to") || srch.Contains("600 To") ||
                      srch.Contains("600 TO") || srch.Contains("600-") || srch.Contains("600 -") || srch.Contains("600_to") || srch.Contains("600_To") ||
                      srch.Contains("600_TO") )
                first = _QCDMuEnriched_600to800;
            else if ( srch.Contains("470to") || srch.Contains("470To") || srch.Contains("470TO") || srch.Contains("470 to") || srch.Contains("470 To") ||
                      srch.Contains("470 TO") || srch.Contains("470-") || srch.Contains("470 -") || srch.Contains("470_to") || srch.Contains("470_To") ||
                      srch.Contains("470_TO") )
                first = _QCDMuEnriched_470to600;
            else if ( srch.Contains("300to") || srch.Contains("300To") || srch.Contains("300TO") || srch.Contains("300 to") || srch.Contains("300 To") ||
                      srch.Contains("300 TO") || srch.Contains("300-") || srch.Contains("300 -") || srch.Contains("300_to") || srch.Contains("300_To") ||
                      srch.Contains("300_TO") )
                first = _QCDMuEnriched_300to470;
            else if ( srch.Contains("170to") || srch.Contains("170To") || srch.Contains("170TO") || srch.Contains("170 to") || srch.Contains("170 To") ||
                      srch.Contains("170 TO") || srch.Contains("170-") || srch.Contains("170 -") || srch.Contains("170_to") || srch.Contains("170_To") ||
                      srch.Contains("170_TO") )
                first = _QCDMuEnriched_170to300;
            else if ( srch.Contains("120to") || srch.Contains("120To") || srch.Contains("120TO") || srch.Contains("120 to") || srch.Contains("120 To") ||
                      srch.Contains("120 TO") || srch.Contains("120-") || srch.Contains("120 -") || srch.Contains("120_to") || srch.Contains("120_To") ||
                      srch.Contains("120_TO") )
                first = _QCDMuEnriched_120to170;
            else if ( srch.Contains("80to") || srch.Contains("80To") || srch.Contains("80TO") || srch.Contains("80 to") || srch.Contains("80 To") ||
                      srch.Contains("80 TO") || srch.Contains("80-") || srch.Contains("80 -") || srch.Contains("80_to") || srch.Contains("80_To") ||
                      srch.Contains("80_TO") )
                first = _QCDMuEnriched_80to120;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") ||
                      srch.Contains("50_TO") )
                first = _QCDMuEnriched_50to80;
            else if ( srch.Contains("30to") || srch.Contains("30To") || srch.Contains("30TO") || srch.Contains("30 to") || srch.Contains("30 To") ||
                      srch.Contains("30 TO") || srch.Contains("30-") || srch.Contains("30 -") || srch.Contains("30_to") || srch.Contains("30_To") ||
                      srch.Contains("30_TO") )
                first = _QCDMuEnriched_30to50;
            else if ( srch.Contains("20to") || srch.Contains("20To") || srch.Contains("20TO") || srch.Contains("20 to") || srch.Contains("20 To") ||
                      srch.Contains("20 TO") || srch.Contains("20-") || srch.Contains("20 -") || srch.Contains("20_to") || srch.Contains("20_To") ||
                      srch.Contains("20_TO") )
                first = _QCDMuEnriched_20to30;
            else if ( srch.Contains("15to") || srch.Contains("15To") || srch.Contains("15TO") || srch.Contains("15 to") || srch.Contains("15 To") ||
                      srch.Contains("15 TO") || srch.Contains("15-") || srch.Contains("15 -") || srch.Contains("15_to") || srch.Contains("15_To") ||
                      srch.Contains("15_TO") )
                first = _QCDMuEnriched_15to20;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _QCDMuEnriched_1000toInf;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _QCDMuEnriched_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _QCDMuEnriched_600to800;
            else if ( srch.Contains("to600") || srch.Contains("To600") || srch.Contains("TO600") || srch.Contains("to 600") || srch.Contains("To 600") ||
                      srch.Contains("TO 600") || srch.Contains("-600") || srch.Contains("- 600") || srch.Contains("to_600") || srch.Contains("To_600") ||
                      srch.Contains("TO_600") )
                last = _QCDMuEnriched_470to600;
            else if ( srch.Contains("to470") || srch.Contains("To470") || srch.Contains("TO470") || srch.Contains("to 470") || srch.Contains("To 470") ||
                      srch.Contains("TO 470") || srch.Contains("-470") || srch.Contains("- 470") || srch.Contains("to_470") || srch.Contains("To_470") ||
                      srch.Contains("TO_470") )
                last = _QCDMuEnriched_300to470;
            else if ( srch.Contains("to300") || srch.Contains("To300") || srch.Contains("TO300") || srch.Contains("to 300") || srch.Contains("To 300") ||
                      srch.Contains("TO 300") || srch.Contains("-300") || srch.Contains("- 300") || srch.Contains("to_300") || srch.Contains("To_300") ||
                      srch.Contains("TO_300") )
                last = _QCDMuEnriched_170to300;
            else if ( srch.Contains("to170") || srch.Contains("To170") || srch.Contains("TO170") || srch.Contains("to 170") || srch.Contains("To 170") ||
                      srch.Contains("TO 170") || srch.Contains("-170") || srch.Contains("- 170") || srch.Contains("to_170") || srch.Contains("To_170") ||
                      srch.Contains("TO_170") )
                last = _QCDMuEnriched_120to170;
            else if ( srch.Contains("to120") || srch.Contains("To120") || srch.Contains("TO120") || srch.Contains("to 120") || srch.Contains("To 120") ||
                      srch.Contains("TO 120") || srch.Contains("-120") || srch.Contains("- 120") || srch.Contains("to_120") || srch.Contains("To_120") ||
                      srch.Contains("TO_120") )
                last = _QCDMuEnriched_80to120;
            else if ( srch.Contains("to80") || srch.Contains("To80") || srch.Contains("TO80") || srch.Contains("to 80") || srch.Contains("To 80") ||
                      srch.Contains("TO 80") || srch.Contains("-80") || srch.Contains("- 80") || srch.Contains("to_80") || srch.Contains("To_80") ||
                      srch.Contains("TO_80") )
                last = _QCDMuEnriched_50to80;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _QCDMuEnriched_30to50;
            else if ( srch.Contains("to30") || srch.Contains("To30") || srch.Contains("TO30") || srch.Contains("to 30") || srch.Contains("To 30") ||
                      srch.Contains("TO 30") || srch.Contains("-30") || srch.Contains("- 30") || srch.Contains("to_30") || srch.Contains("To_30") ||
                      srch.Contains("TO_30") )
                last = _QCDMuEnriched_20to30;
            else if ( srch.Contains("to20") || srch.Contains("To20") || srch.Contains("TO20") || srch.Contains("to 20") || srch.Contains("To 20") ||
                      srch.Contains("TO 20") || srch.Contains("-20") || srch.Contains("- 20") || srch.Contains("to_20") || srch.Contains("To_20") ||
                      srch.Contains("TO_20") )
                last = _QCDMuEnriched_15to20;
            else if ( srch.Contains("to15") || srch.Contains("To15") || srch.Contains("TO15") || srch.Contains("to 15") || srch.Contains("To 15") ||
                      srch.Contains("TO 15") || srch.Contains("-15") || srch.Contains("- 15") || srch.Contains("to_15") || srch.Contains("To_15") ||
                      srch.Contains("TO_15") )
                last = _EndOf_WJets;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _QCDMuEnriched_15to20 && last == _QCDMuEnriched_1000toInf )
            {
                Result.push_back(_QCDMuEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_QCDMuEnriched_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_QCDMuEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_QCDMuEnriched_Full] << "." << endl;
            }

        }// end of if(MuEnriched)
        else if ( srch.Contains("EM") || srch.Contains("em") || srch.Contains("Em") || srch.Contains("eM") )
        {
            // Checking for various intervals
            if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                 srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                 srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                 srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                 srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                 srch.Contains("Inf -") )
                first = _EndOf_QCDEMEnriched_Normal;
            else if ( srch.Contains("300to") || srch.Contains("300To") || srch.Contains("300TO") || srch.Contains("300 to") || srch.Contains("300 To") ||
                      srch.Contains("300 TO") || srch.Contains("300-") || srch.Contains("300 -") || srch.Contains("300_to") || srch.Contains("300_To") ||
                      srch.Contains("300_TO") )
                first = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("170to") || srch.Contains("170To") || srch.Contains("170TO") || srch.Contains("170 to") || srch.Contains("170 To") ||
                      srch.Contains("170 TO") || srch.Contains("170-") || srch.Contains("170 -") || srch.Contains("170_to") || srch.Contains("170_To") ||
                      srch.Contains("170_TO") )
                first = _QCDEMEnriched_170to300;
            else if ( srch.Contains("120to") || srch.Contains("120To") || srch.Contains("120TO") || srch.Contains("120 to") || srch.Contains("120 To") ||
                      srch.Contains("120 TO") || srch.Contains("120-") || srch.Contains("120 -") || srch.Contains("120_to") || srch.Contains("120_To") ||
                      srch.Contains("120_TO") )
                first = _QCDEMEnriched_120to170;
            else if ( srch.Contains("80to") || srch.Contains("80To") || srch.Contains("80TO") || srch.Contains("80 to") || srch.Contains("80 To") ||
                      srch.Contains("80 TO") || srch.Contains("80-") || srch.Contains("80 -") || srch.Contains("80_to") || srch.Contains("80_To") ||
                      srch.Contains("80_TO") )
                first = _QCDEMEnriched_80to120;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") ||
                      srch.Contains("50_TO") )
                first = _QCDEMEnriched_50to80;
            else if ( srch.Contains("30to") || srch.Contains("30To") || srch.Contains("30TO") || srch.Contains("30 to") || srch.Contains("30 To") ||
                      srch.Contains("30 TO") || srch.Contains("30-") || srch.Contains("30 -") || srch.Contains("30_to") || srch.Contains("30_To") ||
                      srch.Contains("30_TO") )
                first = _QCDEMEnriched_30to50;
            else if ( srch.Contains("20to") || srch.Contains("20To") || srch.Contains("20TO") || srch.Contains("20 to") || srch.Contains("20 To") ||
                      srch.Contains("20 TO") || srch.Contains("20-") || srch.Contains("20 -") || srch.Contains("20_to") || srch.Contains("20_To") ||
                      srch.Contains("20_TO") )
                first = _QCDEMEnriched_20to30;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("to600") || srch.Contains("To600") || srch.Contains("TO600") || srch.Contains("to 600") || srch.Contains("To 600") ||
                      srch.Contains("TO 600") || srch.Contains("-600") || srch.Contains("- 600") || srch.Contains("to_600") || srch.Contains("To_600") ||
                      srch.Contains("TO_600") )
                last = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("to470") || srch.Contains("To470") || srch.Contains("TO470") || srch.Contains("to 470") || srch.Contains("To 470") ||
                      srch.Contains("TO 470") || srch.Contains("-470") || srch.Contains("- 470") || srch.Contains("to_470") || srch.Contains("To_470") ||
                      srch.Contains("TO_470") )
                last = _QCDEMEnriched_300toInf;
            else if ( srch.Contains("to300") || srch.Contains("To300") || srch.Contains("TO300") || srch.Contains("to 300") || srch.Contains("To 300") ||
                      srch.Contains("TO 300") || srch.Contains("-300") || srch.Contains("- 300") || srch.Contains("to_300") || srch.Contains("To_300") ||
                      srch.Contains("TO_300") )
                last = _QCDEMEnriched_170to300;
            else if ( srch.Contains("to170") || srch.Contains("To170") || srch.Contains("TO170") || srch.Contains("to 170") || srch.Contains("To 170") ||
                      srch.Contains("TO 170") || srch.Contains("-170") || srch.Contains("- 170") || srch.Contains("to_170") || srch.Contains("To_170") ||
                      srch.Contains("TO_170") )
                last = _QCDEMEnriched_120to170;
            else if ( srch.Contains("to120") || srch.Contains("To120") || srch.Contains("TO120") || srch.Contains("to 120") || srch.Contains("To 120") ||
                      srch.Contains("TO 120") || srch.Contains("-120") || srch.Contains("- 120") || srch.Contains("to_120") || srch.Contains("To_120") ||
                      srch.Contains("TO_120") )
                last = _QCDEMEnriched_80to120;
            else if ( srch.Contains("to80") || srch.Contains("To80") || srch.Contains("TO80") || srch.Contains("to 80") || srch.Contains("To 80") ||
                      srch.Contains("TO 80") || srch.Contains("-80") || srch.Contains("- 80") || srch.Contains("to_80") || srch.Contains("To_80") ||
                      srch.Contains("TO_80") )
                last = _QCDEMEnriched_50to80;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _QCDEMEnriched_30to50;
            else if ( srch.Contains("to30") || srch.Contains("To30") || srch.Contains("TO30") || srch.Contains("to 30") || srch.Contains("To 30") ||
                      srch.Contains("TO 30") || srch.Contains("-30") || srch.Contains("- 30") || srch.Contains("to_30") || srch.Contains("To_30") ||
                      srch.Contains("TO_30") )
                last = _QCDEMEnriched_20to30;
            else if ( srch.Contains("to20") || srch.Contains("To20") || srch.Contains("TO20") || srch.Contains("to 20") || srch.Contains("To 20") ||
                      srch.Contains("TO 20") || srch.Contains("-20") || srch.Contains("- 20") || srch.Contains("to_20") || srch.Contains("To_20") ||
                      srch.Contains("TO_20") )
                last = _EndOf_QCDMuEnriched_Normal;

            // Swapping first with last if necessary
            if ( int(first)>int(last)  && last!=_None)
            {
                Process_t NewLast = Process_t(int(first)-1);
                first = Process_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _QCDEMEnriched_20to30 && last == _QCDEMEnriched_300toInf )
            {
                Result.push_back(_QCDEMEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_QCDEMEnriched_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            else
            {
                Result.push_back(_QCDEMEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_QCDEMEnriched_Full] << "." << endl;
            }

        }// end of if(EMEnriched)

    }// end of if(QCD)

    else if (  srch.Contains("Double") || srch.Contains("double") || srch.Contains("DOUBLE") )
    {
        if ( srch.Contains("Eg") || srch.Contains("eg") || srch.Contains("EG") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                 srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                 srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                 srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                first = _DoubleEG_B;
            else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                      srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                      srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                      srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                      srch.Contains("c_To") || srch.Contains("c_TO") )
                first = _DoubleEG_C;
            else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                      srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                      srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                      srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                      srch.Contains("d_To") || srch.Contains("d_TO") )
                first = _DoubleEG_D;
            else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                      srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                      srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                      srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                      srch.Contains("e_To") || srch.Contains("e_TO") )
                first = _DoubleEG_E;
            else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                      srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                      srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                      srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                      srch.Contains("f_To") || srch.Contains("f_TO") )
                first = _DoubleEG_F;
            else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                      srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                      srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                      srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                      srch.Contains("g_To") || srch.Contains("g_TO") )
                first = _DoubleEG_G;
            else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                      srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                      srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                      srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                      srch.Contains("h_To") || srch.Contains("h_TO") )
                first = _DoubleEG_H;

            if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                 srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                 srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                 srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                 srch.Contains("To_b") || srch.Contains("TO_b") )
                last = _DoubleEG_B;
            else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                      srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                      srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                      srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                      srch.Contains("To_c") || srch.Contains("TO_c") )
                last = _DoubleEG_C;
            else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                      srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                      srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                      srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                      srch.Contains("To_d") || srch.Contains("TO_d") )
                last = _DoubleEG_D;
            else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                      srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                      srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                      srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                      srch.Contains("To_e") || srch.Contains("TO_e") )
                last = _DoubleEG_E;
            else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                      srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                      srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                      srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                      srch.Contains("To_f") || srch.Contains("TO_f") )
                last = _DoubleEG_F;
            else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                      srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                      srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                      srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                      srch.Contains("To_g") || srch.Contains("TO_g") )
                last = _DoubleEG_G;
            else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                      srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                      srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                      srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                      srch.Contains("To_h") || srch.Contains("TO_h") )
                last = _DoubleEG_H;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if ( first == _DoubleEG_B && last == _DoubleEG_H )
            {
                Result.push_back(_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                Result.push_back(_DoubleEG_B);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                Result.push_back(_DoubleEG_C);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                Result.push_back(_DoubleEG_D);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_D] << "." << endl;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                Result.push_back(_DoubleEG_E);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_E] << "." << endl;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                Result.push_back(_DoubleEG_F);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_F] << "." << endl;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                Result.push_back(_DoubleEG_G);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_G] << "." << endl;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                Result.push_back(_DoubleEG_H);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_H] << "." << endl;
            }
            else
            {
                Result.push_back(_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_DoubleEG_Full] << "." << endl;
            }

        }// end of if (EG)

    }// end of if(Double)

    else if ( srch.Contains("Single") || srch.Contains("single") || srch.Contains("SINGLE") )
    {
        if ( srch.Contains("Mu") || srch.Contains("mu") || srch.Contains("MU") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                 srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                 srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                 srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                first = _SingleMuon_B;
            else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                      srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                      srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                      srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                      srch.Contains("c_To") || srch.Contains("c_TO") )
                first = _SingleMuon_C;
            else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                      srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                      srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                      srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                      srch.Contains("d_To") || srch.Contains("d_TO") )
                first = _SingleMuon_D;
            else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                      srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                      srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                      srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                      srch.Contains("e_To") || srch.Contains("e_TO") )
                first = _SingleMuon_E;
            else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                      srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                      srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                      srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                      srch.Contains("f_To") || srch.Contains("f_TO") )
                first = _SingleMuon_F;
            else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                      srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                      srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                      srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                      srch.Contains("g_To") || srch.Contains("g_TO") )
                first = _SingleMuon_G;
            else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                      srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                      srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                      srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                      srch.Contains("h_To") || srch.Contains("h_TO") )
                first = _SingleMuon_H;

            if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                 srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                 srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                 srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                 srch.Contains("To_b") || srch.Contains("TO_b") )
                last = _SingleMuon_B;
            else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                      srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                      srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                      srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                      srch.Contains("To_c") || srch.Contains("TO_c") )
                last = _SingleMuon_C;
            else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                      srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                      srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                      srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                      srch.Contains("To_d") || srch.Contains("TO_d") )
                last = _SingleMuon_D;
            else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                      srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                      srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                      srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                      srch.Contains("To_e") || srch.Contains("TO_e") )
                last = _SingleMuon_E;
            else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                      srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                      srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                      srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                      srch.Contains("To_f") || srch.Contains("TO_f") )
                last = _SingleMuon_F;
            else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                      srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                      srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                      srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                      srch.Contains("To_g") || srch.Contains("TO_g") )
                last = _SingleMuon_G;
            else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                      srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                      srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                      srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                      srch.Contains("To_h") || srch.Contains("TO_h") )
                last = _SingleMuon_H;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if ( first == _SingleMuon_B && last == _SingleMuon_H )
            {
                Result.push_back(_SingleMuon_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_SingleMuon_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                Result.push_back(_SingleMuon_B);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                Result.push_back(_SingleMuon_C);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                Result.push_back(_SingleMuon_D);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_D] << "." << endl;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                Result.push_back(_SingleMuon_E);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_E] << "." << endl;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                Result.push_back(_SingleMuon_F);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_F] << "." << endl;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                Result.push_back(_SingleMuon_G);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_G] << "." << endl;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                Result.push_back(_SingleMuon_H);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_H] << "." << endl;
            }
            else
            {
                Result.push_back(_SingleMuon_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleMuon_Full] << "." << endl;
            }
        }// end of if(SingleMuon)

        if ( srch.Contains("Elec") || srch.Contains("elec") || srch.Contains("ELEC") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                 srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                 srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                 srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                first = _SingleElectron_B;
            else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                      srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                      srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                      srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                      srch.Contains("c_To") || srch.Contains("c_TO") )
                first = _SingleElectron_C;
            else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                      srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                      srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                      srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                      srch.Contains("d_To") || srch.Contains("d_TO") )
                first = _SingleElectron_D;
            else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                      srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                      srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                      srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                      srch.Contains("e_To") || srch.Contains("e_TO") )
                first = _SingleElectron_E;
            else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                      srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                      srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                      srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                      srch.Contains("f_To") || srch.Contains("f_TO") )
                first = _SingleElectron_F;
            else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                      srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                      srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                      srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                      srch.Contains("g_To") || srch.Contains("g_TO") )
                first = _SingleElectron_G;
            else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                      srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                      srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                      srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                      srch.Contains("h_To") || srch.Contains("h_TO") )
                first = _SingleElectron_H;

            if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                 srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                 srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                 srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                 srch.Contains("To_b") || srch.Contains("TO_b") )
                last = _SingleElectron_B;
            else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                      srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                      srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                      srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                      srch.Contains("To_c") || srch.Contains("TO_c") )
                last = _SingleElectron_C;
            else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                      srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                      srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                      srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                      srch.Contains("To_d") || srch.Contains("TO_d") )
                last = _SingleElectron_D;
            else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                      srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                      srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                      srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                      srch.Contains("To_e") || srch.Contains("TO_e") )
                last = _SingleElectron_E;
            else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                      srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                      srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                      srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                      srch.Contains("To_f") || srch.Contains("TO_f") )
                last = _SingleElectron_F;
            else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                      srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                      srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                      srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                      srch.Contains("To_g") || srch.Contains("TO_g") )
                last = _SingleElectron_G;
            else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                      srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                      srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                      srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                      srch.Contains("To_h") || srch.Contains("TO_h") )
                last = _SingleElectron_H;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                Process_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if ( first == _SingleElectron_B && last == _SingleElectron_H )
            {
                Result.push_back(_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( Process_t pr=first; pr<=last; pr=next(pr) )
                {
                    Result.push_back(pr);
                    if ( notify == kTRUE )
                    {
                        if ( pr != last ) cout << Procname[pr] << ", ";
                        else cout << Procname[pr] << "." << endl;
                    }
                }
            }
            // Goes on if there are no intervals
            else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                Result.push_back(_SingleElectron_B);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                Result.push_back(_SingleElectron_C);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                Result.push_back(_SingleElectron_D);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_D] << "." << endl;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                Result.push_back(_SingleElectron_E);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_E] << "." << endl;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                Result.push_back(_SingleElectron_F);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_F] << "." << endl;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                Result.push_back(_SingleElectron_G);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_G] << "." << endl;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                Result.push_back(_SingleElectron_H);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_H] << "." << endl;
            }
            else
            {
                Result.push_back(_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_SingleElectron_Full] << "." << endl;
            }
        }// end of if(SingleElectron)

    }// end of if(Single)

    else if ( srch.Contains("Signal") || srch.Contains("signal") || srch.Contains("SIGNAL") )
    {
        Result.push_back(_DY_Full);
        if ( notify == kTRUE ) cout << Procname[_DY_Full] << "." << endl;
    }

    else if ( srch.Contains("Bkg") || srch.Contains("bkg") || srch.Contains("BKG") || srch.Contains("Background") || srch.Contains("background") ||
              srch.Contains("BACKGROUND") )
    {
        Result.push_back(_bkg_Full);
        if ( notify == kTRUE )  cout << Procname[_bkg_Full] << "." << endl;
    }

    else if ( srch.Contains("Test") || srch.Contains("test") || srch.Contains("TEST") )
    {
        if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
             srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_Test_EE);
            if ( notify == kTRUE )  cout << Procname[_Test_EE] << "." << endl;
        }
        else if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
                  srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_Test_MuMu);
            if ( notify == kTRUE )  cout << Procname[_Test_MuMu] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_Test_EMu);
            if ( notify == kTRUE )  cout << Procname[_Test_EMu] << "." << endl;
        }
    }

    else
    {
        Result.push_back(_None);
        if ( notify == kTRUE ) cout << Procname[_None] << "." << endl;
    }

    if ( instaGet == kTRUE )    // Applying GetProc if asked
    {
        cout << "Getting information for returned processes..";
        for ( Int_t i=0; i<(int(Result.size())); i++ )
        {
            if ( Result[i] != _None ) this->GetProc( Result[i], kFALSE );
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
    Procname[_EndOf_WJets] = "EndOf_WJets";
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
    Procname[_EndOf_SinglMuon_Normal] = "EndOf_SinglMuon_Normal";
    Procname[_SingleElectron_B] = "SingleElectron_B";
    Procname[_SingleElectron_C] = "SingleElectron_C";
    Procname[_SingleElectron_D] = "SingleElectron_D";
    Procname[_SingleElectron_E] = "SingleElectron_E";
    Procname[_SingleElectron_F] = "SingleElectron_F";
    Procname[_SingleElectron_G] = "SingleElectron_G";
    Procname[_SingleElectron_H] = "SingleElectron_H";
    Procname[_EndOf_SingleElectron_Normal] = "EndOf_SingleElectron_Normal";
    Procname[_EndOf_Data_Normal] = "EndOf_Data_Normal";
    Procname[_DY_Full] = "DY_Full";
    Procname[_DYMuMu_Full] = "DYMuMu_Full";
    Procname[_DYEE_Full] = "DYEE_Full";
    Procname[_EndOf_MCsignal_Special] = "EndOf_MCsignal_Special";
    Procname[_DYTauTau_Full] = "DYTauTau_Full";
    Procname[_ttbar_Full] = "ttbar_Full";
    Procname[_VVnST] = "VVnST";
    Procname[_QCDMuEnriched_Full] = "QCDMuEnriched_Full";
    Procname[_QCDEMEnriched_Full] = "QCDEMEnriched_Full";
    Procname[_bkg_Full] = "bkg_Full";
    Procname[_EndOf_MCbkg_Special] = "EndOf_MCbkg_Special";
    Procname[_DoubleEG_Full] = "DoubleEG_Full";
    Procname[_SingleMuon_Full] = "SingleMuon_Full";
    Procname[_SingleElectron_Full] = "SingleElectron_Full";
    Procname[_EndOf_Data_Special] = "EndOf_Data_Special";
    Procname[_Test_MuMu] = "Test_MuMu";
    Procname[_Test_EE] = "Test_EE";
    Procname[_Test_EMu] = "Test_EMu";
    Procname[_EndOf_Test] = "EndOf_Test";
    return;

}// end of PrepareProcNames()


void FileMgr::CheckProcesses()
{
    Bool_t allOk = kTRUE;
    cout << "Checking processes: " << endl;
    for ( Process_t pr=_DY_10to50; pr<_EndOf_Data_Special; pr=next(pr) )
    {
        this->GetProc(pr, kTRUE);
        if ( !Procname[pr].Length() )
        {
            cout << "Process " << pr << ": no Procname found (providing its int expression)." << endl;
            allOk = kFALSE;
        }
        if ( !CurrentProc )
        {
            cout << "Process " << Procname[pr] << ": no CurrentProc found." << endl;
            allOk = kFALSE;
        }
        if ( !Tag.size() )
        {
            cout << "Process " << Procname[pr] << ": no Tags found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(Tag.size())); i++ )
            {
                if ( !Tag[i].Length() )
                {
                    cout << "Process " << Procname[pr] << ": no Tag[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if ( !FileLocation.size() )
        {
            cout << "Process " << Procname[pr] << ": no FileLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for ( Int_t i=0; i<(int(FileLocation.size())); i++ )
            {
                if ( !FileLocation[i].Length() )
                {
                    cout << "Process " << Procname[pr] << ": no FileLocation[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
                else if ( pr<_EndOf_Data_Special && !FileLocation[i].Contains("/*.root") )
                {
                    cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not end with '/*.root'." << endl;
                    allOk = kFALSE;
                }
                else if ( pr>_EndOf_Data_Special && !FileLocation[i].Contains(".root") )
                {
                    cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not end with '.root'." << endl;
                    allOk = kFALSE;
                }
            }

        }
        if ( !FullLocation.size() )
        {
            cout << "Process " << Procname[pr] << ": no FullLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for ( Int_t i=0; i<(int(FullLocation.size())); i++ )
            {
                if ( !FullLocation[i].Length() )
                {
                    cout << "Process " << Procname[pr] << ": no FullLocation[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
                else
                {
                    if ( pr<_EndOf_Data_Special && !FullLocation[i].Contains("/*.root") )
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not end with '/*.root'" << endl;
                        allOk = kFALSE;
                    }
                    if ( pr>_EndOf_Data_Special && !FullLocation[i].Contains(".root") )
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not end with '.root'" << endl;
                        allOk = kFALSE;
                    }
                    if ( FullLocation[i][0] != '/')
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not begin with '/'" << endl;
                        allOk = kFALSE;
                    }
                }
            }

        }
        if ( !TreeName.size() )
        {
            cout << "Process " << Procname[pr] << ": no TreeName found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for ( Int_t i=0; i<(int(TreeName.size())); i++ )
            {
                if ( !TreeName[i].Length() )
                {
                    cout << "Process " << Procname[pr] << ": no TreeName[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if ( !BaseLocation.Length() )
        {
            cout << "Process " << Procname[pr] << ": no BaseLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if ( BaseLocation[BaseLocation.Length()-1] != '/') {
                cout << "Process " << Procname[pr] << ": BaseLocation does not end with '/'" << endl;
                allOk = kFALSE;
            }
            if ( BaseLocation[0] != '/') {
                cout << "Process " << Procname[pr] << ": BaseLocation does not begin with '/'" << endl;
                allOk = kFALSE;
            }

        }
        if ( !nEvents.size() )
        {
            cout << "Process " << Procname[pr] << ": no nEvents found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(nEvents.size())); i++)
            {
                if ( !nEvents[i] )
                {
                    cout << "Process " << Procname[pr] << ": no nEvents[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if ( pr < _DoubleEG_B || ( pr > _EndOf_Data_Normal && pr < _EndOf_MCbkg_Special ) || ( pr > _EndOf_Data_Special && pr < _EndOf_Test ) )
        {
            if ( isMC == kFALSE )
            {
                cout << "Process " << Procname[pr] << ": is said to be NOT MC." << endl;
                allOk = kFALSE;
            }
            if ( !Xsec.size() )
            {
                cout << "Process " << Procname[pr] << ": no Xsec found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<(int(Xsec.size())); i++)
                {
                    if ( !Xsec[i] )
                    {
                        cout << "Process " << Procname[pr] << ": no Xsec[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
            if ( !Wsum.size() )
            {
                cout << "Process " << Procname[pr] << ": no Wsum found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<(int(Wsum.size())); i++)
                {
                    if ( !Wsum[i] )
                    {
                        cout << "Process " << Procname[pr] << ": no Wsum[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
            if ( Tag.size() != FileLocation.size() || FileLocation.size() != FullLocation.size() || FullLocation.size() != TreeName.size() ||
                 TreeName.size() != Xsec.size() || Xsec.size() != Wsum.size() || Wsum.size() != nEvents.size() )
            {
                cout << "Process " << Procname[pr] << " Vector sizes do not match." << endl;
                allOk = kFALSE;
            }
        }// end of if(MC)
        else
        {
            if ( isMC == kTRUE )
            {
                cout << "Process " << Procname[pr] << ": is said to be MC." << endl;
                allOk = kFALSE;
            }
            if ( Tag.size() != FileLocation.size() || FileLocation.size() != FullLocation.size() || FullLocation.size() != TreeName.size() ||
                 TreeName.size() != nEvents.size() )
            {
                cout << "Process " << Procname[pr] << " Vector sizes do not match." << endl;
                allOk = kFALSE;
            }
        }// end of else()
        if ( !Type.Length() )
        {
            cout << "Process " << Procname[pr] << ": no Type found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if ( pr < _EndOf_MCsignal_Normal || ( pr > _EndOf_Data_Normal && pr < _EndOf_MCsignal_Special ) )
            {
                if ( Type != "SIGNAL" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT SIGNAL." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( ( pr > _EndOf_MCsignal_Normal && pr < _EndOf_MCbkg_Normal ) || ( pr > _EndOf_MCsignal_Special && pr < _EndOf_MCbkg_Special ) )
            {
                if ( Type != "BKG" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT BKG." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( ( pr > _EndOf_MCbkg_Normal && pr < _EndOf_Data_Normal ) || ( pr > _EndOf_MCbkg_Special && pr < _EndOf_Data_Special ) )
            {
                if ( Type != "DATA" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT DATA." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( pr > _EndOf_Data_Special && pr < _EndOf_Test )
            {
                if ( Type != "TEST" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT TEST." << endl;
                    allOk = kFALSE;
                }
            }
        }       
        vector<Process_t> forChecking = FindProc(Procname[pr], kFALSE, kFALSE);
        if ( !forChecking.size() )
        {
            cout << "Process " << Procname[pr] << " No output recorded from FindProc()." << endl;
            allOk = kFALSE;
        }
        else if ( forChecking[0] != pr )
        {
            cout << "Process " << Procname[pr] << ": FindProc() did not find this process." << endl;
            allOk = kFALSE;
        }
        this->ClearProc();
        if ( CurrentProc != _None )
        {
            cout << "Process " << Procname[pr] << ": Current proc is not _None after ClearProc()." << endl;
            allOk = kFALSE;
        }
        if ( Tag.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Tag." << endl;
            allOk = kFALSE;
        }
        if ( Xsec.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Xsec." << endl;
            allOk = kFALSE;
        }
        if ( Wsum.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Wsum." << endl;
            allOk = kFALSE;
        }
        if ( nEvents.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the nEvents." << endl;
            allOk = kFALSE;
        }
        if ( BaseLocation.Length() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the BaseLocation." << endl;
            allOk = kFALSE;
        }
        if ( FileLocation.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Tag." << endl;
            allOk = kFALSE;
        }
        if ( FullLocation.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the FullLocation." << endl;
            allOk = kFALSE;
        }
        if ( TreeName.size() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the TreeName." << endl;
            allOk = kFALSE;
        }
        if ( Type.Length() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the Type." << endl;
            allOk = kFALSE;
        }
    }// end of for()
    if ( allOk == kTRUE ) cout << "All OK." << endl;
    else cout << "Problems were detected." << endl;

}// end of CheckProcesses()
