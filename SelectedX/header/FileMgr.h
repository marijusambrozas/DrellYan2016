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
    _EndOf_Data_Special
};

inline
Process_t next(Process_t pr)    // Processes that begin with "EndOf" will be skipped by this
{
  if ( pr == _EndOf_Data_Special )
      return pr;
  else if ( pr == _DYEE_2000to3000 || pr == _QCDEMEnriched_300toInf || pr == _SingleElectron_H )
      return Process_t(int(pr)+3);
  else if ( pr == _DY_2000to3000 || pr == _DYMuMu_2000to3000 || pr == _EndOf_DYEE_Normal || pr == _DYTauTau_50toInf || pr == _ttbar_1000toInf ||
            pr == _WW || pr == _WJets || pr == _QCDMuEnriched_1000toInf || pr == _EndOf_QCDEMEnriched_Normal || pr == _DoubleEG_H ||
            pr == _SingleMuon_H || pr == _EndOf_SingleElectron_Normal || pr == _DYEE_Full || pr == _bkg_Full )
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
        vector<Double_t> Xsec;
        vector<Double_t> Wsum;
        vector<Double_t> nEvents;
        TString BaseLocation;
        Bool_t isMC;

        // -- Constructors -- //
        FileMgr ( Process_t pr = _None );
        FileMgr ( TString search = "" );

        Process_t FindProc ( TString search );
        void NextProc ();
        void GetProc ( Process_t pr = _None, Bool_t ClearOld = kTRUE );
        void ClearProc ();
        void CheckProcesses ();
};

FileMgr::FileMgr ( Process_t pr )
{
        CurrentProc = pr;
        this->GetProc(CurrentProc, kTRUE);
}

FileMgr::FileMgr ( TString search )
{
        CurrentProc = this->FindProc(search);
        this->GetProc(CurrentProc, kTRUE);
}

void FileMgr::GetProc ( Process_t pr, Bool_t ClearOld )
{
    if ( ClearOld == kTRUE )
    {
        Tag.clear();
        FullLocation.clear();
        FileLocation.clear();
        Xsec.clear();
        Wsum.clear();
        nEvents.clear();
    }
    TString Location;
    CurrentProc = pr;

    if ( pr == _DY_10to50 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M10to50_v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_v2" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_ext1v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr == _DY_50to100 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M50to100" ); Xsec.push_back( 5869.58346 ); Wsum.push_back( 81780984 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_100to200 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M100to200" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200_ext" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_200to400 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M200to400" ); Xsec.push_back( 7.67 ); Wsum.push_back( 169676 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_400to500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M400to500" ); Xsec.push_back( 0.423 ); Wsum.push_back( 151190 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_500to700 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M500to700" ); Xsec.push_back( 0.24 ); Wsum.push_back( 144096 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_700to800 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M700to800" ); Xsec.push_back( 0.035 ); Wsum.push_back( 136892 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_800to1000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M800to1000" ); Xsec.push_back( 0.03 ); Wsum.push_back( 131586 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_1000to1500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M1000to1500" ); Xsec.push_back( 0.016 ); Wsum.push_back( 120010 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_1500to2000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M1500to2000" ); Xsec.push_back( 0.002 ); Wsum.push_back( 111709 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DY_2000to3000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M2000to3000" ); Xsec.push_back( 0.00054 ); Wsum.push_back( 101298 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DY_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DY_M10to50_v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_v2" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M10to50_ext1v1" ); Xsec.push_back( 18610.0 ); Wsum.push_back( 22301710 + 47946333 + 29386420 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DY_M50to100" ); Xsec.push_back( 5869.58346 ); Wsum.push_back( 81780984 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M100to200_ext" ); Xsec.push_back( 226 ); Wsum.push_back( 703034 + 9607589 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M200to400" ); Xsec.push_back( 7.67 ); Wsum.push_back( 169676 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M400to500" ); Xsec.push_back( 0.423 ); Wsum.push_back( 151190 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M500to700" ); Xsec.push_back( 0.24 ); Wsum.push_back( 144096 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M700to800" ); Xsec.push_back( 0.035 ); Wsum.push_back( 136892 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M800to1000" ); Xsec.push_back( 0.03 ); Wsum.push_back( 131586 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M1000to1500" ); Xsec.push_back( 0.016 ); Wsum.push_back( 120010 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M1500to2000" ); Xsec.push_back( 0.002 ); Wsum.push_back( 111709 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DY_M2000to3000" ); Xsec.push_back( 0.00054 ); Wsum.push_back( 101298 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYMuMu_10to50 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr == _DYMuMu_50to100 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_100to200 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3433295.0 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3433295.0 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_200to400 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_400to500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_500to700 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_700to800 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_800to1000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_1000to1500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_1500to2000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYMuMu_2000to3000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYMuMu_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33278866.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3433295.0 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3433295.0 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_10to50 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 306508623 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_50to100 ) // Only EE evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_100to200 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3437885.0 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3437885.0 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_200to400 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_400to500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_500to700 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_700to800 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_800to1000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _DYEE_1000to1500 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_1500to2000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_2000to3000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYEE_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 306508623 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33275218.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 122055296 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3437885.0 ); nEvents.push_back( 38422582 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 3437885.0 ); nEvents.push_back( 15120677 );
        Location = "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 295242 );
        Location = "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 287262 );
        Location = "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 280940 );
        Location = "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 276234 );
        Location = "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 271768 );
        Location = "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 258620 );
        Location = "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 258625 );
        Location = "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 255342 );
        Location = "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DYTauTau_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar_700to1000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _ttbar_1000toInf )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ttbar_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _tW )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _tbarW )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _ZZ )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _WZ )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _WW )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _VVnST )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _WJets )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 23944342 ); // I get Wsum=137540054
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 177139200 );
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/*.root";        // There also is madgraph version
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_15to20 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_20to30 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_30to50 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_50to80 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_80to120 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_120to170 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_170to300 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_300to470 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_470to600 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_600to800 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_800to1000 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_1000toInf )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _QCDMuEnriched_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_20to30 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_30to50 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_50to80 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_80to120 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_120to170 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_170to300 )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_300toInf )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _QCDEMEnriched_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _bkg_Full )
    {
        isMC = kTRUE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 30650862 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 65887977 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 ); nEvents.push_back( 40381246 );
        Location = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 122055296 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77081149 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); nEvents.push_back( 77867729 );
        Location = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 38422582 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 24561630 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 6952830 );
        Location = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 6933093 );
        Location = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 998034 );
        Location = "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000/*.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 2995828 );
        Location = "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000/*.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 6987123 );
        Location = "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000/*.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 23944342 ); // I get Wsum=137540054
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 177139200 );
        Location = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000/*.root";        // There also is madgraph version
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
        Location = "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_B )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_B_0000" ); nEvents.push_back( 103625724 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_B_0001" ); nEvents.push_back( 33031246 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_C )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_C" ); nEvents.push_back( 45521797 );
        Location = "DoubleEG/crab_DoubleEG_RunC/180326_143612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_D )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_D" ); nEvents.push_back( 52422569 );
        Location = "DoubleEG/crab_DoubleEG_RunD/180326_143654/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_E )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_E" ); nEvents.push_back( 47326656 );
        Location = "DoubleEG/crab_DoubleEG_RunE/180326_143750/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_F )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_F" ); nEvents.push_back( 33943052 );
        Location = "DoubleEG/crab_DoubleEG_RunF/180326_143846/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_G )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_G_0000" ); nEvents.push_back( 71864512 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0001" ); nEvents.push_back( 4669958 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_H )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_Hver2_0000" ); nEvents.push_back( 68821231 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0001" ); nEvents.push_back( 11645108 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver3" ); nEvents.push_back( 2021309 );
        Location = "DoubleEG/crab_DoubleEG_RunHver3/180326_144719/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _DoubleEG_Full )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "DoubleEG_B_0000" ); nEvents.push_back( 103625724 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_B_0001" ); nEvents.push_back( 33031246 );
        Location = "DoubleEG/crab_DoubleEG_RunB/180326_143532/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_C" ); nEvents.push_back( 45521797 );
        Location = "DoubleEG/crab_DoubleEG_RunC/180326_143612/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_D" ); nEvents.push_back( 52422569 );
        Location = "DoubleEG/crab_DoubleEG_RunD/180326_143654/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_E" ); nEvents.push_back( 47326656 );
        Location = "DoubleEG/crab_DoubleEG_RunE/180326_143750/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_F" ); nEvents.push_back( 33943052 );
        Location = "DoubleEG/crab_DoubleEG_RunF/180326_143846/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0000" ); nEvents.push_back( 71864512 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_G_0001" ); nEvents.push_back( 4669958 );
        Location = "DoubleEG/crab_DoubleEG_RunG/180326_144559/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0000" ); nEvents.push_back( 68821231 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver2_0001" ); nEvents.push_back( 11645108 );
        Location = "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0001/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "DoubleEG_Hver3" ); nEvents.push_back( 2021309 );
        Location = "DoubleEG/crab_DoubleEG_RunHver3/180326_144719/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_B )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_B" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunB/180326_143105/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_C )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_C" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunC/180326_143152/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_D )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_D" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunD/180326_143257/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_E )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_E" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunE/180326_143338/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_F )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_F" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunF/180326_143419/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_G )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_G" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunG/180326_144335/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_H )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunHver2/180326_144412/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunHver3/180326_144454/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleMuon_Full )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleMuon_B" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunB/180326_143105/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_C" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunC/180326_143152/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_D" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunD/180326_143257/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_E" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunE/180326_143338/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_F" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunF/180326_143419/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_G" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunG/180326_144335/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunHver2/180326_144412/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleMuon_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleMuon/crab_SingleMuon_RunHver3/180326_144454/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_B )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_B" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_C )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_C" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_D )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_D" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_E )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_E" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_F )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_G )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_H )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _SingleElectron_Full )
    {
        isMC = kFALSE;
        BaseLocation = "/xrd/store/user/dpai/_v2p3_/";

        Tag.push_back( "SingleElectron_B" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunB/180326_143935/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_C" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunC/180326_144015/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_D" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunD/180326_144117/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_E" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunE/180326_144202/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunF/180326_144247/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_F" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunG/180326_144755/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver2" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver2/180326_144832/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SingleElectron_Hver3" ); nEvents.push_back( 0 );
        Location = "SingleElectron/crab_SingleElectron_RunHver3/180326_144908/0000/*.root";
        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
}// end of Get()

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
        isMC = kFALSE;
        this->GetProc(CurrentProc, kTRUE);
    }
}

Process_t FileMgr::FindProc ( TString search )
{
    TString srch = search;
    if ( srch.Contains("DY") || srch.Contains("dy") || srch.Contains("dY") || srch.Contains("Dy") || srch.Contains("DrellYan") || srch.Contains("DRELLYAN") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            if ( srch.Contains("10to50") || srch.Contains("10To50") || srch.Contains("10-50") )
            {
                return _DYMuMu_10to50;
            }
            else if ( srch.Contains("50to100") || srch.Contains("50To100") || srch.Contains("50-100") )
            {
                return _DYMuMu_50to100;
            }
            else if ( srch.Contains("100to200") || srch.Contains("100To200") || srch.Contains("100-200") )
            {
                return _DYMuMu_100to200;
            }
            else if ( srch.Contains("200to400") || srch.Contains("200To400") || srch.Contains("200-400") )
            {
                return _DYMuMu_200to400;
            }
            else if ( srch.Contains("400to500") || srch.Contains("400To500") || srch.Contains("400-500") )
            {
                return _DYMuMu_400to500;
            }
            else if ( srch.Contains("500to700") || srch.Contains("500To700") || srch.Contains("500-700") )
            {
                return _DYMuMu_500to700;
            }
            else if ( srch.Contains("700to800") || srch.Contains("700To800") || srch.Contains("700-800") )
            {
                return _DYMuMu_700to800;
            }
            else if ( srch.Contains("800to1000") || srch.Contains("800To1000") || srch.Contains("800-1000") )
            {
                return _DYMuMu_800to1000;
            }
            else if ( srch.Contains("1000to1500") || srch.Contains("1000To1500") || srch.Contains("1000-1500") )
            {
                return _DYMuMu_1000to1500;
            }
            else if ( srch.Contains("1500to2000") || srch.Contains("1500To2000") || srch.Contains("1500-2000") )
            {
                return _DYEE_1500to2000;
            }
            else if ( srch.Contains("2000to3000") || srch.Contains("2000To3000") || srch.Contains("2000-3000") )
            {
                return _DYMuMu_2000to3000;
            }
            else
            {
                return _DYMuMu_Full;
            }
        }// end of if(DYMuMu)

        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
             srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            if ( srch.Contains("10to50") || srch.Contains("10To50") || srch.Contains("10-50") )
            {
                return _DYEE_10to50;
            }
            else if ( srch.Contains("50to100") || srch.Contains("50To100") || srch.Contains("50-100") )
            {
                return _DYEE_50to100;
            }
            else if ( srch.Contains("100to200") || srch.Contains("100To200") || srch.Contains("100-200") )
            {
                return _DYEE_100to200;
            }
            else if ( srch.Contains("200to400") || srch.Contains("200To400") || srch.Contains("200-400") )
            {
                return _DYEE_200to400;
            }
            else if ( srch.Contains("400to500") || srch.Contains("400To500") || srch.Contains("400-500") )
            {
                return _DYEE_400to500;
            }
            else if ( srch.Contains("500to700") || srch.Contains("500To700") || srch.Contains("500-700") )
            {
                return _DYEE_500to700;
            }
            else if ( srch.Contains("700to800") || srch.Contains("700To800") || srch.Contains("700-800") )
            {
                return _DYEE_700to800;
            }
            else if ( srch.Contains("800to1000") || srch.Contains("800To1000") || srch.Contains("800-1000") )
            {
                return _DYEE_800to1000;
            }
            else if ( srch.Contains("1000to1500") || srch.Contains("1000To1500") || srch.Contains("1000-1500") )
            {
                return _DYEE_1000to1500;
            }
            else if ( srch.Contains("1500to2000") || srch.Contains("1500To2000") || srch.Contains("1500-2000") )
            {
                return _DYEE_1500to2000;
            }
            else if ( srch.Contains("2000to3000") || srch.Contains("2000To3000") || srch.Contains("2000-3000") )
            {
                return _DYEE_2000to3000;
            }
            else
            {
                return _DYEE_Full;
            }
        }// end of if(DYEE)

        else if ( srch.Contains("TauTau") || srch.Contains("tautau") || srch.Contains("TAUTAU") || srch.Contains("Ditau") || srch.Contains("DiTau") ||
             srch.Contains("ditau") || srch.Contains("diTau") || srch.Contains("DITAU") )
        {
            if ( srch.Contains("10to50") || srch.Contains("10To50") || srch.Contains("10-50") )
            {
                return _DYTauTau_10to50;
            }
            else if ( srch.Contains("50to100") || srch.Contains("50To100") || srch.Contains("50-100") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("100to200") || srch.Contains("100To200") || srch.Contains("100-200") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("200to400") || srch.Contains("200To400") || srch.Contains("200-400") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("400to500") || srch.Contains("400To500") || srch.Contains("400-500") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("500to700") || srch.Contains("500To700") || srch.Contains("500-700") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("700to800") || srch.Contains("700To800") || srch.Contains("700-800") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("800to1000") || srch.Contains("800To1000") || srch.Contains("800-1000") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("1000to1500") || srch.Contains("1000To1500") || srch.Contains("1000-1500") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("1500to2000") || srch.Contains("1500To2000") || srch.Contains("1500-2000") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("2000to3000") || srch.Contains("2000To3000") || srch.Contains("2000-3000") )
            {
                return _DYTauTau_50toInf;
            }
            else if ( srch.Contains("50toInf") || srch.Contains("50ToInf") || srch.Contains("50-Inf") || srch.Contains("50-inf") || srch.Contains("50-INF") ||
                      srch.Contains("50toINF") || srch.Contains("50ToINF") )
            {
                return _DYTauTau_50toInf;
            }
            else
            {
                return _DYTauTau_Full;
            }
        }// end of if(DYTauTau)
        else
        {
            if ( srch.Contains("10to50") || srch.Contains("10To50") || srch.Contains("10-50") )
            {
                return _DY_10to50;
            }
            else if ( srch.Contains("50to100") || srch.Contains("50To100") || srch.Contains("50-100") )
            {
                return _DY_50to100;
            }
            else if ( srch.Contains("100to200") || srch.Contains("100To200") || srch.Contains("100-200") )
            {
                return _DY_100to200;
            }
            else if ( srch.Contains("200to400") || srch.Contains("200To400") || srch.Contains("200-400") )
            {
                return _DY_200to400;
            }
            else if ( srch.Contains("400to500") || srch.Contains("400To500") || srch.Contains("400-500") )
            {
                return _DY_400to500;
            }
            else if ( srch.Contains("500to700") || srch.Contains("500To700") || srch.Contains("500-700") )
            {
                return _DY_500to700;
            }
            else if ( srch.Contains("700to800") || srch.Contains("700To800") || srch.Contains("700-800") )
            {
                return _DY_700to800;
            }
            else if ( srch.Contains("800to1000") || srch.Contains("800To1000") || srch.Contains("800-1000") )
            {
                return _DY_800to1000;
            }
            else if ( srch.Contains("1000to1500") || srch.Contains("1000To1500") || srch.Contains("1000-1500") )
            {
                return _DY_1000to1500;
            }
            else if ( srch.Contains("1500to2000") || srch.Contains("1500To2000") || srch.Contains("1500-2000") )
            {
                return _DY_1500to2000;
            }
            else if ( srch.Contains("2000to3000") || srch.Contains("2000To3000") || srch.Contains("2000-3000") )
            {
                return _DY_2000to3000;
            }
            else
            {
                return _DY_Full;
            }
        }// end of else (DY)

    }// end of if(DrellYan)

    else if ( srch.Contains("ttbar") || srch.Contains("tt") || srch.Contains("TTbar") || srch.Contains("TT") || srch.Contains("ttBAR") ||
         srch.Contains("TTBAR") || srch.Contains("TantiT") || srch.Contains("tantit") || srch.Contains("TANTIT") || srch.Contains("tANTIt") ||
         srch.Contains("TopAntiTop") )
    {
        if ( srch.Contains("700to1000") || srch.Contains("700To1000") || srch.Contains("700-1000") )
        {
            return _ttbar_700to1000;
        }
        else if ( srch.Contains("1000toInf") || srch.Contains("1000ToInf") || srch.Contains("1000-Inf") || srch.Contains("1000-inf") ||
                  srch.Contains("1000-INF") ||
                  srch.Contains("1000toINF") || srch.Contains("1000ToINF") )
        {
            return _ttbar_1000toInf;
        }
        else if ( srch.Contains("full") || srch.Contains("Full") || srch.Contains("FULL") )
        {
            return _ttbar_Full;
        }
        else
        {
            return _ttbar;
        }
    }// end of if(ttbar)

    else if ( srch.Contains("tW") || srch.Contains("TW") || srch.Contains("SingleTop") || srch.Contains("singleTop") || srch.Contains("singletop") ||
         srch.Contains("SingleTOP") || srch.Contains("singleTOP") || srch.Contains("SINGLETOP") )
    {
        return _tW;
    }

    else if ( srch.Contains("tbarW") || srch.Contains("TbarW") || srch.Contains("TBarW") || srch.Contains("tBarW") || srch.Contains("TBARW") ||
              srch.Contains("SingleAntiTop") || srch.Contains("SingleAntitop") || srch.Contains("singleAntiTop") || srch.Contains("singleAntitop") ||
              srch.Contains("SingleAntiTOP") || srch.Contains("singleAntiTOP") || srch.Contains("SingleANTITOP") || srch.Contains("singleANTITOP") ||
              srch.Contains("SINGLEANTITOP") )
    {
        return _tbarW;
    }

    else if ( srch.Contains("ZZ") || srch.Contains("zz") )
    {
        return _ZZ;
    }

    else if ( srch.Contains("WZ") || srch.Contains("wz") || srch.Contains("ZW") || srch.Contains("zw") )
    {
        return _WZ;
    }

    else if ( srch.Contains("WW") || srch.Contains("ww") )
    {
        return _WW;
    }

    else if ( srch.Contains("DIBOSON") || srch.Contains("diboson") || srch.Contains("Diboson") || srch.Contains("diBoson"))
    {
        return _VVnST;
    }

    else if ( srch.Contains("VVnST") || srch.Contains("VVNST") || srch.Contains("vvnst") || srch.Contains("VVnst") || srch.Contains("vvnST") ||
              srch.Contains("VV") || srch.Contains("vv") )
    {
        return _VVnST;
    }

    else if ( srch.Contains("WJets") || srch.Contains("Wjets") || srch.Contains("wjets") || srch.Contains("WJETS") || srch.Contains("W+jets") ||
              srch.Contains("W+Jets") || srch.Contains("W+JETS") )
    {
        return _WJets;
    }

    else if ( srch.Contains("QCD") || srch.Contains("qcd") )
    {
        if ( srch.Contains("Mu") || srch.Contains("mu") || srch.Contains("MU") )
        {
            if ( srch.Contains("15to20") || srch.Contains("15To20") || srch.Contains("15-20") )
            {
                return _QCDMuEnriched_15to20;
            }
            else if ( srch.Contains("20to30") || srch.Contains("20To30") || srch.Contains("20-30") )
            {
                return _QCDMuEnriched_20to30;
            }
            else if ( srch.Contains("30to50") || srch.Contains("30To50") || srch.Contains("30-50") )
            {
                return _QCDMuEnriched_30to50;
            }
            else if ( srch.Contains("50to80") || srch.Contains("50To80") || srch.Contains("50-80") )
            {
                return _QCDMuEnriched_50to80;
            }
            else if ( srch.Contains("80to120") || srch.Contains("80To120") || srch.Contains("80-120") )
            {
                return _QCDMuEnriched_80to120;
            }
            else if ( srch.Contains("120to170") || srch.Contains("120To170") || srch.Contains("120-170") )
            {
                return _QCDMuEnriched_120to170;
            }
            else if ( srch.Contains("170to300") || srch.Contains("170To300") || srch.Contains("170-300") )
            {
                return _QCDMuEnriched_170to300;
            }
            else if ( srch.Contains("300to470") || srch.Contains("300To470") || srch.Contains("300-470") )
            {
                return _QCDMuEnriched_300to470;
            }
            else if ( srch.Contains("470to600") || srch.Contains("470To600") || srch.Contains("470-600") )
            {
                return _QCDMuEnriched_470to600;
            }
            else if ( srch.Contains("600to800") || srch.Contains("600To800") || srch.Contains("600-800") )
            {
                return _QCDMuEnriched_600to800;
            }
            else if ( srch.Contains("800to1000") || srch.Contains("800To1000") || srch.Contains("800-1000") )
            {
                return _QCDMuEnriched_800to1000;
            }
            else if ( srch.Contains("1000toInf") || srch.Contains("1000ToInf") || srch.Contains("1000-Inf") || srch.Contains("1000-inf") ||
                      srch.Contains("1000-INF") || srch.Contains("1000toINF") || srch.Contains("1000ToINF") )
            {
                return _QCDMuEnriched_1000toInf;
            }
            else
            {
                return _QCDMuEnriched_Full;

            }
        }// end of if(MuEnriched)
        else if ( srch.Contains("EM") || srch.Contains("em") || srch.Contains("Em") || srch.Contains("eM") )
        {
            if ( srch.Contains("20to30") || srch.Contains("20To30") || srch.Contains("20-30") )
            {
                return _QCDEMEnriched_20to30;
            }
            else if ( srch.Contains("30to50") || srch.Contains("30To50") || srch.Contains("30-50") )
            {
                return _QCDEMEnriched_30to50;
            }
            else if ( srch.Contains("50to80") || srch.Contains("50To80") || srch.Contains("50-80") )
            {
                return _QCDEMEnriched_50to80;
            }
            else if ( srch.Contains("80to120") || srch.Contains("80To120") || srch.Contains("80-120") )
            {
                return _QCDEMEnriched_80to120;
            }
            else if ( srch.Contains("120to170") || srch.Contains("120To170") || srch.Contains("120-170") )
            {
                return _QCDEMEnriched_120to170;
            }
            else if ( srch.Contains("170to300") || srch.Contains("170To300") || srch.Contains("170-300") )
            {
                return _QCDEMEnriched_170to300;
            }
            else if ( srch.Contains("300toInf") || srch.Contains("300ToInf") || srch.Contains("300-Inf") || srch.Contains("300-inf") ||
                      srch.Contains("300-INF") || srch.Contains("300toINF") || srch.Contains("300ToINF") )
            {
                return _QCDEMEnriched_300toInf;
            }
            else
            {
                return _QCDEMEnriched_Full;

            }
        }// end of if(EMEnriched)

    }// end of if(QCD)

    else if (  srch.Contains("Double") || srch.Contains("double") || srch.Contains("DOUBLE") )
    {
        if ( srch.Contains("Eg") || srch.Contains("eg") || srch.Contains("EG") )
        {
            if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                return _DoubleEG_B;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                return _DoubleEG_C;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                return _DoubleEG_D;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                return _DoubleEG_E;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                return _DoubleEG_F;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                return _DoubleEG_G;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                return _DoubleEG_H;
            }
            else
            {
                return _DoubleEG_Full;
            }
        }// end of if (EG)

    }// end of if(Double)

    else if ( srch.Contains("Single") || srch.Contains("single") || srch.Contains("SINGLE") )
    {
        if ( srch.Contains("Mu") || srch.Contains("mu") || srch.Contains("MU") )
        {
            if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                return _SingleMuon_B;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                return _SingleMuon_C;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                return _SingleMuon_D;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                return _SingleMuon_E;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                return _SingleMuon_F;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                return _SingleMuon_G;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                return _SingleMuon_H;
            }
            else
            {
                return _SingleMuon_Full;
            }
        }// end of if(SingleMuon)

        if ( srch.Contains("Elec") || srch.Contains("elec") || srch.Contains("ELEC") )
        {
            if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                return _SingleElectron_B;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                return _SingleElectron_C;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                return _SingleElectron_D;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                return _SingleElectron_E;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                return _SingleElectron_F;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                return _SingleElectron_G;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
            {
                return _SingleElectron_H;
            }
            else
            {
                return _SingleElectron_Full;
            }
        }// end of if(SingleElectron)

    }// end of if(Single)

    else if ( srch.Contains("Signal") || srch.Contains("signal") || srch.Contains("SIGNAL") )
    {
        return _DY_Full;
    }

    else if ( srch.Contains("Bkg") || srch.Contains("bkg") || srch.Contains("BKG") || srch.Contains("Background") || srch.Contains("background") ||
              srch.Contains("BACKGROUND") )
    {
        return _bkg_Full;
    }

    else
    {
        return _None;
    }
    return _None;
}// end of FindProc()

void FileMgr::CheckProcesses()
{
    Bool_t allOk = kTRUE;
    cout << "Checking processes: " << endl;
    for ( Process_t pr=_DYEE_10to50; pr<_EndOf_Data_Special; pr=next(pr) )
    {
        this->GetProc(pr, kTRUE);
        if ( !CurrentProc )
        {
            cout << "Process " << pr << ": no CurrentProc found." << endl;
            allOk = kFALSE;
        }
        if ( !Tag.size() )
        {
            cout << "Process " << pr << ": no Tags found." << endl;
            allOk = kFALSE;
        }
        if ( !FileLocation.size() )
        {
            cout << "Process " << pr << ": no FileLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for ( Int_t i=0; i<FileLocation.size(); i++ )
            {
                if ( !FileLocation[i].Contains("/*.root") ) {
                    cout << "Process " << pr << ": FileLocation[" << i << "] does not have '/*.root'" << endl;
                    allOk = kFALSE;
                }
            }

        }
        if ( !FullLocation.size() )
        {
            cout << "Process " << pr << ": no FullLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for ( Int_t i=0; i<FullLocation.size(); i++ )
            {
                if ( !FullLocation[i].Contains("/*.root") ) {
                    cout << "Process " << pr << ": FullLocation[" << i << "] does not have '/*.root'" << endl;
                    allOk = kFALSE;
                }
                if ( FullLocation[i][0] != '/') {
                    cout << "Process " << pr << ": FullLocation[" << i << "] does not begin with '/'" << endl;
                    allOk = kFALSE;
                }
            }

        }
        if ( !BaseLocation.Length() )
        {
            cout << "Process " << pr << ": no BaseLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if ( BaseLocation[BaseLocation.Length()-1] != '/') {
                cout << "Process " << pr << ": BaseLocation does not end with '/'" << endl;
                allOk = kFALSE;
            }
            if ( BaseLocation[0] != '/') {
                cout << "Process " << pr << ": BaseLocation does not begin with '/'" << endl;
                allOk = kFALSE;
            }

        }
        if ( !nEvents.size() )
        {
            cout << "Process " << pr << ": no nEvents found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<nEvents.size(); i++)
            {
                if ( nEvents[i] == 0 )
                {
                    cout << "Process " << pr << ": no nEvents[" << i << "] found." << endl;
                    allOk = kFALSE;
                }
            }
        }
        if ( pr < _DoubleEG_B || ( pr < _DoubleEG_Full && pr > _EndOf_Data_Normal ))
        {
            if ( isMC == kFALSE )
            {
                cout << "Process " << pr << ": is said to be NOT MC." << endl;
                allOk = kFALSE;
            }
            if ( !Xsec.size() )
            {
                cout << "Process " << pr << ": no Xsec found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<Xsec.size(); i++)
                {
                    if ( Xsec[i] == 0 )
                    {
                        cout << "Process " << pr << ": no Xsec[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
            if ( !Wsum.size() )
            {
                cout << "Process " << pr << ": no Wsum found." << endl;
                allOk = kFALSE;
            }
            else
            {
                for (Int_t i=0; i<Wsum.size(); i++)
                {
                    if ( Wsum[i] == 0 )
                    {
                        cout << "Process " << pr << ": no Wsum[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
        }// end of if(MC)
        else
        {
            if ( isMC == kTRUE )
            {
                cout << "Process " << pr << ": is said to be MC." << endl;
                allOk = kFALSE;
            }
        }// end of else()
    }
}
