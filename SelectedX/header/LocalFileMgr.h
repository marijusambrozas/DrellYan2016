// Header for file management
// First working version created by Marijus Ambrozas 2018.08.23
// 2018.08.14: Added HistLocation for every process
// 2018.08.31: Added an ability to change between files that contain events selected before and after the Rochester correction (you need to change ROCCORR variable to switch)
// 2019.02.25: Updated to v2.6 data.

#pragma once

#include <TGraphErrors.h>
#include <TH2.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

#define Lumi 35867 // -- from Run2016B to Run2016H, JSON. unit: /pb, Updated at 2017.07.30 -- //
#define Lumi_HLTv4p2 865.919 // -- integrated luminosity before Run 257933 -- //
#define nMassBin 43

using namespace std;

enum SelProc_t
{
    _None=0,
    // "Normal" processes - related to the divisions of files
    _MuMu_DY_10to50, _MuMu_DY_50to100, _MuMu_DY_100to200, _MuMu_DY_200to400, _MuMu_DY_400to500, _MuMu_DY_500to700, _MuMu_DY_700to800, _MuMu_DY_800to1000,
    _MuMu_DY_1000to1500, _MuMu_DY_1500to2000, _MuMu_DY_2000to3000,
    _EndOf_MuMu_MCsignal_Normal,
    _MuMu_DYTauTau_10to50, _MuMu_DYTauTau_50toInf, _EndOf_MuMu_DYTauTau_Normal,
    _MuMu_ttbar, _MuMu_ttbar_700to1000, _MuMu_ttbar_1000toInf, _EndOf_MuMu_ttbar_Normal,
    _MuMu_tW, _MuMu_tbarW, _MuMu_ZZ, _MuMu_WZ, _MuMu_WW, _EndOf_MuMu_VVnST_Normal,
    _MuMu_WJets, _MuMu_WJets_ext2v5, _EndOf_MuMu_WJets_Normal,
    _MuMu_QCDMuEnriched_15to20, _MuMu_QCDMuEnriched_20to30, _MuMu_QCDMuEnriched_30to50, _MuMu_QCDMuEnriched_50to80, _MuMu_QCDMuEnriched_80to120,
    _MuMu_QCDMuEnriched_120to170, _MuMu_QCDMuEnriched_170to300, _MuMu_QCDMuEnriched_300to470, _MuMu_QCDMuEnriched_470to600, _MuMu_QCDMuEnriched_600to800,
    _MuMu_QCDMuEnriched_800to1000, _MuMu_QCDMuEnriched_1000toInf,
    _EndOf_MuMu_QCDMuEnriched_Normal,
    _EndOf_MuMu_MCbkg_Normal,
    _MuMu_SingleMuon_B, _MuMu_SingleMuon_C, _MuMu_SingleMuon_D, _MuMu_SingleMuon_E, _MuMu_SingleMuon_F, _MuMu_SingleMuon_G, _MuMu_SingleMuon_H,
    _EndOf_MuMu_Data_Normal,
    _MuMu_DY_Full,
    _MuMu_DYTauTau_Full, _MuMu_ttbar_Full, _MuMu_VVnST, _MuMu_WJets_Full, _MuMu_QCDMuEnriched_Full, _MuMu_Bkg_Full,
    _MuMu_SingleMuon_Full,
    _EndOf_MuMu_Special,
    _Test_MuMu,
    _EndOf_MuMu,

    _EE_DY_10to50, _EE_DY_50to100, _EE_DY_100to200, _EE_DY_200to400, _EE_DY_400to500, _EE_DY_500to700, _EE_DY_700to800, _EE_DY_800to1000, _EE_DY_1000to1500,
    _EE_DY_1500to2000, _EE_DY_2000to3000,
    _EndOf_EE_MCsignal_Normal,
    _EE_DYTauTau_10to50, _EE_DYTauTau_50toInf, _EndOf_EE_DYTauTau_Normal,
    _EE_ttbar, _EE_ttbar_700to1000, _EE_ttbar_1000toInf, _EndOf_EE_ttbar_Normal,
    _EE_tW, _EE_tbarW, _EE_ZZ, _EE_WZ, _EE_WW, _EndOf_EE_VVnST_Normal,
    _EE_WJets, _EE_WJets_ext2v5, _EndOf_EE_WJets_Normal,
    _EE_QCDEMEnriched_20to30, _EE_QCDEMEnriched_30to50, _EE_QCDEMEnriched_50to80, _EE_QCDEMEnriched_80to120, _EE_QCDEMEnriched_120to170, _EE_QCDEMEnriched_170to300,
    _EE_QCDEMEnriched_300toInf,
    _EndOf_EE_QCDEMEnriched_Normal,
    _EndOf_EE_MCbkg_Normal,
    _EE_DoubleEG_B, _EE_DoubleEG_C, _EE_DoubleEG_D, _EE_DoubleEG_E, _EE_DoubleEG_F, _EE_DoubleEG_G, _EE_DoubleEG_H,
    _EndOf_EE_DoubleEG_Normal,
    _EE_SingleElectron_B, _EE_SingleElectron_C, _EE_SingleElectron_D, _EE_SingleElectron_E, _EE_SingleElectron_F, _EE_SingleElectron_G, _EE_SingleElectron_H,
    _EndOf_EE_SingleElectron_Normal,
    _EndOf_EE_Data_Normal,
    _EE_DY_Full,
    _EE_DYTauTau_Full, _EE_ttbar_Full, _EE_VVnST, _EE_WJets_Full, _EE_QCDEMEnriched_Full, _EE_Bkg_Full,
    _EE_DoubleEG_Full,  _EE_SingleElectron_Full,
    _EndOf_EE_Special,
    _Test_EE,
    _EndOf_EE,

    _EMu_DYTauTau_10to50, _EMu_DYTauTau_50toInf, _EndOf_EMu_DYTauTau_Normal,
    _EMu_ttbar, _EMu_ttbar_700to1000, _EMu_ttbar_1000toInf, _EndOf_EMu_ttbar_Normal,
    _EMu_tW, _EMu_tbarW, _EMu_ZZ, _EMu_WZ, _EMu_WW, _EndOf_EMu_VVnST_Normal,
    _EMu_WJets, _EMu_WJets_ext2v5, _EndOf_EMu_WJets_Normal,
    _EndOf_EMu_MCbkg_Normal,
    _EMu_SingleMuon_B, _EMu_SingleMuon_C, _EMu_SingleMuon_D, _EMu_SingleMuon_E, _EMu_SingleMuon_F, _EMu_SingleMuon_G, _EMu_SingleMuon_H,
    _EndOf_EMu_Data_Normal,
    _EMu_DYTauTau_Full, _EMu_ttbar_Full, _EMu_VVnST, _EMu_WJets_Full, _EMu_Bkg_Full,
    _EMu_SingleMuon_Full,
    _EndOf_EMu_Special,
    _Test_EMu,
    _EndOf_EMu
};

inline
SelProc_t next ( SelProc_t pr )    // Processes that begin with "EndOf" will be skipped by this
{
  if ( pr == _EndOf_EMu )
      return pr;
  else if ( pr == _MuMu_QCDMuEnriched_1000toInf || pr == _EE_QCDEMEnriched_300toInf || pr == _EE_SingleElectron_H || pr == _EMu_WJets_ext2v5 )
      return SelProc_t(int(pr)+3);
  else if ( pr == _MuMu_DY_2000to3000 || pr == _MuMu_DYTauTau_50toInf || pr == _MuMu_ttbar_1000toInf || pr == _MuMu_WW || pr == _MuMu_WJets_ext2v5 ||
            pr == _EndOf_MuMu_QCDMuEnriched_Normal || pr == _MuMu_SingleMuon_H || pr == _MuMu_SingleMuon_Full || pr == _Test_MuMu ||
            pr == _EE_DY_2000to3000 || pr == _EE_DYTauTau_50toInf || pr == _EE_ttbar_1000toInf || pr == _EE_WW || pr == _EE_WJets_ext2v5 ||
            pr == _EndOf_EE_QCDEMEnriched_Normal || pr == _EE_DoubleEG_H || pr == _EndOf_EE_SingleElectron_Normal ||  pr == _EE_SingleElectron_Full || pr == _Test_EE ||
            pr == _EMu_DYTauTau_50toInf || pr == _EMu_ttbar_1000toInf || pr == _EMu_WW || pr == _EMu_SingleMuon_H || pr == _EMu_SingleMuon_Full )
      return SelProc_t(int(pr)+2);
  else
      return SelProc_t(int(pr)+1);
}


class LocalFileMgr
{
public:

        SelProc_t CurrentProc;
        vector<TString> Tag;
        vector<TString> FullLocation;
        vector<TString> FileLocation;
        vector<TString> TreeName;
        vector<Double_t> Xsec;
        vector<Double_t> Wsum;
        vector<Double_t> nEvents;

        TString BaseLocation;
        TString HistLocation;
        TString Type;       
        Bool_t isMC;       

        map<SelProc_t, TString> Procname;

        // -- Constructor -- //
        LocalFileMgr ( SelProc_t pr = _None );

        vector<SelProc_t> FindProc ( TString search, Bool_t notify = kTRUE, Bool_t instaGet = kFALSE );
        void NextProc ();
        void SetProc ( SelProc_t pr = _None, Bool_t ClearOld = kTRUE );
        void ClearProc ();

        void SwitchROCCORR();

private:
        Bool_t namesSet = kFALSE;
        Bool_t processesChecked = kFALSE;
        Bool_t ROCCORR = kFALSE; // DEPRECATED

        void PrepareProcNames ();
        void CheckProcesses ();

};// end of class definition


// ---------- Constructor ---------- //

LocalFileMgr::LocalFileMgr ( SelProc_t pr )
{
    if ( namesSet == kFALSE ) { this->PrepareProcNames(); namesSet = kTRUE; }
    if ( processesChecked == kFALSE ) { this->CheckProcesses(); processesChecked = kTRUE; }
    CurrentProc = pr;
    this->SetProc(CurrentProc, kTRUE);
}


// ----------- Functions ----------- //

void LocalFileMgr::NextProc()
{
    CurrentProc = next(CurrentProc);
    this->SetProc(CurrentProc, kTRUE);
}


void LocalFileMgr::ClearProc()
{
    if ( CurrentProc != _None )
    {
        CurrentProc = _None;
        BaseLocation = "";
        HistLocation = "";
        Type = "";
        isMC = kFALSE;
        this->SetProc(CurrentProc, kTRUE);
    }
}


void LocalFileMgr::SetProc ( SelProc_t pr, Bool_t ClearOld )
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
    TString RocCorr = "";
    CurrentProc = pr;

    if ( ROCCORR == kTRUE ) RocCorr = "_roccor";

    if ( pr != _None )BaseLocation = "/media/sf_DATA/";

    if ( pr == _MuMu_DY_10to50 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7471260+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 20366 );
        else nEvents.push_back( 20366 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v2" ); Xsec.push_back( 6016.88  ); Wsum.push_back( 7471260+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43076 );
        else nEvents.push_back( 43076 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 6016.88  ); Wsum.push_back( 7471260+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26401 );
        else nEvents.push_back( 26401 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr ==_MuMu_DY_50to100 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";
        Tag.push_back( "SelectedMuMu_DYMuMu_M50to100" ); Xsec.push_back( 1873.52 ); Wsum.push_back( 26175605.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 12650938 );
        else nEvents.push_back( 12650938 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M50to100"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200" ); Xsec.push_back( 76.2401 ); Wsum.push_back( 3176762.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1849390 );
        else nEvents.push_back( 1849390 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M200to400" ); Xsec.push_back( 2.67606 ); Wsum.push_back( 560322.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 479764 );
        else nEvents.push_back( 479764 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M200to400"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M400to500" ); Xsec.push_back( 0.139728 ); Wsum.push_back( 50136.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 55166 );
        else nEvents.push_back( 55166 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M400to500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M500to700" ); Xsec.push_back( 0.0792496 ); Wsum.push_back( 48188.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 58374 );
        else nEvents.push_back( 58374 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M500to700"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M700to800" ); Xsec.push_back( 0.0123176 ); Wsum.push_back( 44984.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 60781 );
        else nEvents.push_back( 60781 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M700to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M800to1000" ); Xsec.push_back( 0.01042 ); Wsum.push_back( 43496.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 61650 );
        else nEvents.push_back( 61650 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M1000to1500" ); Xsec.push_back( 0.00552772 ); Wsum.push_back( 40110.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 61732 );
        else nEvents.push_back( 61732 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1000to1500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M1500to2000" ); Xsec.push_back( 0.000741613 ); Wsum.push_back( 37176.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 63924 );
        else nEvents.push_back( 63924 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1500to2000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M2000to3000" ); Xsec.push_back( 0.000178737 ); Wsum.push_back( 33360.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 64224 );
        else nEvents.push_back( 64224 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M2000to3000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr ==_MuMu_DY_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v1" ); Wsum.push_back( 7471260+16016651+9815322 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 20366 );
        else nEvents.push_back( 20366 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v2" ); Wsum.push_back( 7471260+16016651+9815322 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43076 );
        else nEvents.push_back( 43076 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_ext1v1" ); Wsum.push_back( 7471260+16016651+9815322 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26401 );
        else nEvents.push_back( 26401 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M50to100" ); Wsum.push_back( 26175605.0 );
        Xsec.push_back( 1873.52 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 12650938 );
        else nEvents.push_back( 12650938 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M50to100"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200" ); Wsum.push_back( 3176762.0 );
        Xsec.push_back( 76.2401 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1849390 );
        else nEvents.push_back( 1849390 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M200to400" ); Wsum.push_back( 560322.0 );
        Xsec.push_back( 2.67606 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 479764 );
        else nEvents.push_back( 479764 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M200to400"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M400to500" ); Wsum.push_back( 50136.0 );
        Xsec.push_back( 0.139728 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 55166 );
        else nEvents.push_back( 55166 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M400to500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M500to700" ); Wsum.push_back( 48188.0 );
        Xsec.push_back( 0.0792496 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 58374 );
        else nEvents.push_back( 58374 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M500to700"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M700to800" ); Wsum.push_back( 44984.0 );
        Xsec.push_back( 0.0123176 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 60781 );
        else nEvents.push_back( 60781 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M700to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M800to1000" ); Wsum.push_back( 43496.0 );
        Xsec.push_back( 0.01042 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 61650 );
        else nEvents.push_back( 61650 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M1000to1500" ); Wsum.push_back( 40110.0 );
        Xsec.push_back( 0.00552772 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 61732 );
        else nEvents.push_back( 61732 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1000to1500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M1500to2000" ); Wsum.push_back( 37176.0 );
        Xsec.push_back( 0.000741613 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 63924 );
        else nEvents.push_back( 63924 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1500to2000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M2000to3000" ); Wsum.push_back( 33360.0 );
        Xsec.push_back( 0.000178737 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 64224 );
        else nEvents.push_back( 64224 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M2000to3000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28 );
        else nEvents.push_back( 28 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28629 );
        else nEvents.push_back( 28629 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Wsum.push_back( 7432050+15912921+9759664 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Wsum.push_back( 7432050+15912921+9759664 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Wsum.push_back( 7432050+15912921+9759664 );
        Xsec.push_back( 6016.88 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28 );
        else nEvents.push_back( 28 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Wsum.push_back( 27277866.0 );  //  NNLO Xsec
        Xsec.push_back( 1952.68432327 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28629 );
        else nEvents.push_back( 28629 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );  //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 414023 );
        else nEvents.push_back( 414023 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 419022 );
        else nEvents.push_back( 419022 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 231295 );
        else nEvents.push_back( 231295 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 123358 );
        else nEvents.push_back( 123358 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ttbar_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 414023 );
        else nEvents.push_back( 414023 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 419022 );
        else nEvents.push_back( 419022 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 231295 );
        else nEvents.push_back( 231295 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 123358 );
        else nEvents.push_back( 123358 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_tW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47435 );
        else nEvents.push_back( 47435 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47291 );
        else nEvents.push_back( 47291 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 27989 );
        else nEvents.push_back( 27989 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43781 );
        else nEvents.push_back( 43781 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32935 );
        else nEvents.push_back( 32935 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47435 );
        else nEvents.push_back( 47435 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47291 );
        else nEvents.push_back( 47291 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 27989 );
        else nEvents.push_back( 27989 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43781 );
        else nEvents.push_back( 43781 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32935 );
        else nEvents.push_back( 32935 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47 );
        else nEvents.push_back( 47 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 485 );
        else nEvents.push_back( 485 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_WJets_ext2v5 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 748 );
        else nEvents.push_back( 748 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_WJets_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47 );
        else nEvents.push_back( 47 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 485 );
        else nEvents.push_back( 485 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 748 );
        else nEvents.push_back( 748 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_15to20 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt15to20"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_20to30 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt20to30"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_30to50 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt30to50"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_50to80 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_80to120 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_120to170 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_170to300 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_300to470 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_470to600 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_600to800 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_800to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt15to20"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt20to30"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt30to50"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_Bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28 );
        else nEvents.push_back( 28 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 );  //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 28629 );
        else nEvents.push_back( 28629 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 414023 );
        else nEvents.push_back( 414023 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 419022 );
        else nEvents.push_back( 419022 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 231295 );
        else nEvents.push_back( 231295 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 123358 );
        else nEvents.push_back( 123358 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47435 );
        else nEvents.push_back( 47435 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47291 );
        else nEvents.push_back( 47291 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 27989 );
        else nEvents.push_back( 27989 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43781 );
        else nEvents.push_back( 43781 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32935 );
        else nEvents.push_back( 32935 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47 );
        else nEvents.push_back( 47 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 485 );
        else nEvents.push_back( 485 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 748 );
        else nEvents.push_back( 748 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 4141251.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt15to20"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 31302080.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt20to30"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 29717171.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt30to50"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 19806914.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 0 );
        else nEvents.push_back( 0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_B ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_B_0000" ); Wsum.push_back( 154033697 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3711485 );
        else nEvents.push_back( 3711485 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_B_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_C ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1599013 );
        else nEvents.push_back( 1599013 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_D ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2654466 );
        else nEvents.push_back( 2654466 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_E ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2434555 );
        else nEvents.push_back( 2434555 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_F ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1877436 );
        else nEvents.push_back( 1877436 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_G ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_G_0000" ); Wsum.push_back( 147937361 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4968832 );
        else nEvents.push_back( 4968832 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_G_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_H ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver2_0000" ); Wsum.push_back( 166695051 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5597930 );
        else nEvents.push_back( 5597930 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver2_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 144964 );
        else nEvents.push_back( 144964 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_Full ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_B_0000" ); Wsum.push_back( 154033697 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3711485 );
        else nEvents.push_back( 3711485 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_B_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1599013 );
        else nEvents.push_back( 1599013 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2654466 );
        else nEvents.push_back( 2654466 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2434555 );
        else nEvents.push_back( 2434555 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1877436 );
        else nEvents.push_back( 1877436 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_G_0000" ); Wsum.push_back( 147937361 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4968832 );
        else nEvents.push_back( 4968832 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_G_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver2_0000" ); Wsum.push_back( 166695051 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5597930 );
        else nEvents.push_back( 5597930 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver2_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 144964 );
        else nEvents.push_back( 144964 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _Test_MuMu )
    {
        isMC = kTRUE;
        Type = "TEST";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ZToMuMu_M4500to6000_4"); Xsec.push_back( 1.0 ); Wsum.push_back( 10200.0 ); nEvents.push_back( 10200 );
        Location = "test/SelectedMuMu_ZToMuMu_M4500to6000_4"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_10to50 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 10438 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 22623 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 13853 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_50to100 ) // Only EE evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M50to100" ); Xsec.push_back( 1873.52 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 7706465 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M50to100.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M100to200" ); Xsec.push_back( 76.2401 ); Wsum.push_back( 3179506 ); nEvents.push_back( 1174892 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M200to400" ); Xsec.push_back( 2.67606 ); Wsum.push_back( 560818.0 ); nEvents.push_back( 337026 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M200to400.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M400to500" ); Xsec.push_back( 0.139728 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 40028 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M400to500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M500to700" ); Xsec.push_back( 0.0792496 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 42760 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M500to700.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M700to800" ); Xsec.push_back( 0.0123176 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 45355 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M700to800.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M800to1000" ); Xsec.push_back( 0.01042 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 46883 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M800to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M1000to1500" ); Xsec.push_back( 0.00552772 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 47242 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1000to1500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M1500to2000" ); Xsec.push_back( 0.000741613 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 50391 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1500to2000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M2000to3000" ); Xsec.push_back( 0.000178737 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 50990 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M2000to3000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M10to50_v1" ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 10438 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_v2" ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 22623 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_ext1v1" ); Wsum.push_back( 7471588+16016761+9811434 ); nEvents.push_back( 13853 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M50to100" ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 7706465 );
        Xsec.push_back( 1873.52 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M50to100.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M100to200" ); Wsum.push_back( 3179506 ); nEvents.push_back( 1174892 );
        Xsec.push_back( 76.2401 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M200to400" ); Wsum.push_back( 560818.0 ); nEvents.push_back( 337026 );
        Xsec.push_back( 2.67606 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M200to400.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M400to500" ); Wsum.push_back( 50420.0 ); nEvents.push_back( 40028 );
        Xsec.push_back( 0.139728 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M400to500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M500to700" ); Wsum.push_back( 48039.0 ); nEvents.push_back( 42760 );
        Xsec.push_back( 0.0792496 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M500to700.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M700to800" ); Wsum.push_back( 46114.0 ); nEvents.push_back( 45355 );
        Xsec.push_back( 0.0123176 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M700to800.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M800to1000" ); Wsum.push_back( 44256.0 ); nEvents.push_back( 46883 );
        Xsec.push_back( 0.01042 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M800to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M1000to1500" ); Wsum.push_back( 39712.0 ); nEvents.push_back( 47242 );
        Xsec.push_back( 0.00552772 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1000to1500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M1500to2000" ); Wsum.push_back( 37287.0 ); nEvents.push_back( 50391 );
        Xsec.push_back( 0.000741613 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1500to2000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M2000to3000" ); Wsum.push_back( 34031.0 ); nEvents.push_back( 50990 );
        Xsec.push_back( 0.000178737 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M2000to3000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }   
    else if ( pr == _EE_DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 7 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 26 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 18 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20255 ); //  NNLO Xsec
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 7 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 26 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 18 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20255 ); //  NNLO Xsec
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277412 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279765 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162430 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89588 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277412 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279765 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162430 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89588 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_tW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32190 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32117 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17817 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28718 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21057 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32190 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32117 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17817 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28718 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21057 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 198 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 1774 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_WJets_ext2v5 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 3002 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext2v5.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_WJets_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 198 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 1774 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 3002 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext2v5.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_20to30 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt20to30.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_30to50 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_50to80 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 2 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_80to120 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_120to170 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back(35817276+41954033 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_170to300 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162.0 ); nEvents.push_back( 3 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt170to300.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_300toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633.0 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt300toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_QCDEMEnriched_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt20to30.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 2 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162.0 ); nEvents.push_back( 3 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt170to300.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633.0 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt300toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_Bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 7 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 26 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Wsum.push_back( 7432050+15912921+9759664 ); nEvents.push_back( 18 );
        Xsec.push_back( 6016.88 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20255 ); //  NNLO Xsec
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277412 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279765 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162430 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89588 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32190 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32117 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17817 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28718 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21057 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 198 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 1774 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253691042 ); nEvents.push_back( 3002 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext2v5.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt20to30" ); Xsec.push_back( 557600000*0.0096 ); Wsum.push_back( 9218952.0 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt20to30.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt30to50_ext1" ); Xsec.push_back( 136000000*0.073 ); Wsum.push_back( 4730195+6768384 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt30to50_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 2 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt50to80_ext1" ); Xsec.push_back( 19800000*0.146 ); Wsum.push_back( 22337068+23474168 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt50to80_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 1 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt80to120_ext1" ); Xsec.push_back( 2800000*0.125 ); Wsum.push_back( 35841780+41853502 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt80to120_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 6 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt120to170_ext1" ); Xsec.push_back( 477000*0.132 ); Wsum.push_back( 35817276+41954033 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt120to170_ext1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt170to300" ); Xsec.push_back( 114000*0.165 ); Wsum.push_back( 11540162.0 ); nEvents.push_back( 3 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt170to300.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_QCDEMEnriched_Pt300toInf" ); Xsec.push_back( 9000*0.15 ); Wsum.push_back( 7373633.0 ); nEvents.push_back( 8 );
        Location = "SelectedEE/MC_bkg/SelectedEE_QCDEMEnriched_Pt300toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_B_0000" ); Wsum.push_back( 103625724 ); nEvents.push_back( 1477919 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_B_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_B_0001" ); Wsum.push_back( 33031246 ); nEvents.push_back( 662642 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_B_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_C" ); Wsum.push_back( 45521797 ); nEvents.push_back( 912737 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_C.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_D" ); Wsum.push_back( 52422569 ); nEvents.push_back( 1596627 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_D.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_E" ); Wsum.push_back( 47326656 ); nEvents.push_back( 1388519 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_E.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_F" ); Wsum.push_back( 33943052 ); nEvents.push_back( 1074053 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_F.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_G_0000" ); Wsum.push_back( 71864512 ); nEvents.push_back( 2708188 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_G_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_G_0001" ); Wsum.push_back( 4669958 ); nEvents.push_back( 176543 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_G_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_Hver2_0000" ); Wsum.push_back( 68821231 ); nEvents.push_back( 2624157 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver2_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_Hver2_0001" ); Wsum.push_back( 11645108 ); nEvents.push_back( 432957 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver2_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_Hver3" ); Wsum.push_back( 2021309 ); nEvents.push_back( 78568 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver3.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DoubleEG_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DoubleEG_B_0000" ); Wsum.push_back( 103625724 ); nEvents.push_back( 1477919 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_B_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_B_0001" ); Wsum.push_back( 33031246 ); nEvents.push_back( 662642 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_B_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_C" ); Wsum.push_back( 45521797 ); nEvents.push_back( 912737 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_C.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_D" ); Wsum.push_back( 52422569 ); nEvents.push_back( 1596627 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_D.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_E" ); Wsum.push_back( 47326656 ); nEvents.push_back( 1388519 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_E.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_F" ); Wsum.push_back( 33943052 ); nEvents.push_back( 1074053 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_F.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_G_0000" ); Wsum.push_back( 71864512 ); nEvents.push_back( 2708188 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_G_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_G_0001" ); Wsum.push_back( 4669958 ); nEvents.push_back( 176543 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_G_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_Hver2_0000" ); Wsum.push_back( 68821231 ); nEvents.push_back( 2624157 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver2_0000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_Hver2_0001" ); Wsum.push_back( 11645108 ); nEvents.push_back( 432957 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver2_0001.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DoubleEG_Hver3" ); Wsum.push_back( 2021309 ); nEvents.push_back( 78568 );
        Location = "SelectedEE/Data/SelectedEE_DoubleEG_Hver3.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_B" ); Wsum.push_back( 174105617 ); nEvents.push_back( 1489369 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_B.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_C" ); Wsum.push_back( 93325367 ); nEvents.push_back( 883477 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_C.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_D" ); Wsum.push_back( 146493465 ); nEvents.push_back( 1545296 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_D.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_E" ); Wsum.push_back( 113168502 ); nEvents.push_back( 1338042 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_E.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_F" ); Wsum.push_back( 70085191 ); nEvents.push_back( 1034009 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_F.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_G" ); Wsum.push_back( 143169219 ); nEvents.push_back( 2627807 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_G.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_Hver2" ); Wsum.push_back( 106262454 ); nEvents.push_back( 2512996 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_Hver2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_Hver3" ); Wsum.push_back( 3187483 ); nEvents.push_back( 75757 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_Hver3.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_SingleElectron_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_SingleElectron_B" ); Wsum.push_back( 174105617 ); nEvents.push_back( 1489369 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_B.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_C" ); Wsum.push_back( 93325367 ); nEvents.push_back( 883477 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_C.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_D" ); Wsum.push_back( 146493465 ); nEvents.push_back( 1545296 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_D.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_E" ); Wsum.push_back( 113168502 ); nEvents.push_back( 1338042 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_E.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_F" ); Wsum.push_back( 70085191 ); nEvents.push_back( 1034009 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_F.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_G" ); Wsum.push_back( 143169219 ); nEvents.push_back( 2627807 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_G.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_Hver2" ); Wsum.push_back( 106262454 ); nEvents.push_back( 2512996 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_Hver2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_SingleElectron_Hver3" ); Wsum.push_back( 3187483 ); nEvents.push_back( 75757 );
        Location = "SelectedEE/Data/SelectedEE_SingleElectron_Hver3.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }    
    else if ( pr == _Test_EE )
    {
        isMC = kTRUE;
        Type = "TEST";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ZToEE_M4500to6000_2"); Xsec.push_back( 1.0 ); Wsum.push_back( 39200.0 ); nEvents.push_back( 39200 );
        Location = "test/SelectedEE_ZToEE_M4500to6000_2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 18 );
        else nEvents.push_back( 18 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35 );
        else nEvents.push_back( 35 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 22 );
        else nEvents.push_back( 22 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32195 );
        else nEvents.push_back( 32195 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 18 );
        else nEvents.push_back( 18 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35 );
        else nEvents.push_back( 35 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 22 );
        else nEvents.push_back( 22 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 ); //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32195 );
        else nEvents.push_back( 32195 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );   //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 553157 );
        else nEvents.push_back( 553157 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 559696 );
        else nEvents.push_back( 559696 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 326098 );
        else nEvents.push_back( 326098 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 176131 );
        else nEvents.push_back( 176131 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ttbar_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );   //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 553157 );
        else nEvents.push_back( 553157 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 559696 );
        else nEvents.push_back( 559696 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 326098 );
        else nEvents.push_back( 326098 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 176131 );
        else nEvents.push_back( 176131 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_tW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65085 );
        else nEvents.push_back( 65085 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65376 );
        else nEvents.push_back( 65376 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1754 );
        else nEvents.push_back( 1754 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8817 );
        else nEvents.push_back( 8817 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43508 );
        else nEvents.push_back( 43508 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65085 );
        else nEvents.push_back( 65085 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65376 );
        else nEvents.push_back( 65376 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1754 );
        else nEvents.push_back( 1754 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8817 );
        else nEvents.push_back( 8817 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43508 );
        else nEvents.push_back( 43508 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 274 );
        else nEvents.push_back( 274 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2602 );
        else nEvents.push_back( 2602 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_WJets_ext2v5 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4103 );
        else nEvents.push_back( 4103 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_WJets_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 274 );
        else nEvents.push_back( 274 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2602 );
        else nEvents.push_back( 2602 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4103 );
        else nEvents.push_back( 4103 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_Bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 18 );
        else nEvents.push_back( 18 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35 );
        else nEvents.push_back( 35 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 6016.88 ); Wsum.push_back( 7432050+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 22 );
        else nEvents.push_back( 22 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1952.68432327 ); Wsum.push_back( 27277866.0 ); //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 32195 );
        else nEvents.push_back( 32195 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 553157 );
        else nEvents.push_back( 553157 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 559696 );
        else nEvents.push_back( 559696 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 326098 );
        else nEvents.push_back( 326098 );
        Wsum.push_back( 38487178.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 176131 );
        else nEvents.push_back( 176131 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65085 );
        else nEvents.push_back( 65085 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65376 );
        else nEvents.push_back( 65376 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1754 );
        else nEvents.push_back( 1754 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8817 );
        else nEvents.push_back( 8817 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 43508 );
        else nEvents.push_back( 43508 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 274 );
        else nEvents.push_back( 274 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2602 );
        else nEvents.push_back( 2602 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext2v5" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 16433848+161144203+253933112 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4103 );
        else nEvents.push_back( 4103 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext2v5"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt15to20" ); Xsec.push_back( 720648000*0.00042 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt15to20"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt20to30" ); Xsec.push_back( 1273190000*0.003 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt20to30"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt30to50" ); Xsec.push_back( 139803000*0.01182 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt30to50"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt50to80" ); Xsec.push_back( 19222500*0.02276 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//// DID NOT FIND THIS ONE
////        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
////        Location = "SelectedEMu/MC_bkg/QCDMuEnriched_Pt470to600_ext2";
////        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "SelectedEMu/MC_bkg/SelectedEMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_B )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_B_0000" ); Wsum.push_back( 154033697 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 60094 );
        else nEvents.push_back( 60094 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_B_0001" ); Wsum.push_back( 108561074 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 49756 );
//        else nEvents.push_back( 49749 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 25104 );
        else nEvents.push_back( 25104 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 41770 );
        else nEvents.push_back( 41770 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 37155 );
        else nEvents.push_back( 37155 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 29019 );
        else nEvents.push_back( 29019 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_G_0000" ); Wsum.push_back( 147937361 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 76472 );
        else nEvents.push_back( 76472 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_G_0001" ); Wsum.push_back( 138710659 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 84902 );
//        else nEvents.push_back( 84896 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_Hver2_0000" ); Wsum.push_back( 166695051 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 84092 );
        else nEvents.push_back( 84092 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_Hver2_0001" ); Wsum.push_back( 141936183 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 85368 );
//        else nEvents.push_back( 85351 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2226 );
        else nEvents.push_back( 2226 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_B_0000" ); Wsum.push_back( 154033697 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 60094 );
        else nEvents.push_back( 60094 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_B_0001" ); Wsum.push_back( 108561074 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 49756 );
//        else nEvents.push_back( 49749 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 25104 );
        else nEvents.push_back( 25104 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 41770 );
        else nEvents.push_back( 41770 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 37155 );
        else nEvents.push_back( 37155 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 29019 );
        else nEvents.push_back( 29019 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_G_0000" ); Wsum.push_back( 147937361 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 76472 );
        else nEvents.push_back( 76472 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_G_0001" ); Wsum.push_back( 138710659 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 84902 );
//        else nEvents.push_back( 84896 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver2_0000" ); Wsum.push_back( 166695051 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 84092 );
        else nEvents.push_back( 84092 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2_0000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

//        Tag.push_back( "SelectedEMu_SingleMuon_Hver2_0001" ); Wsum.push_back( 141936183 );
//        if ( ROCCORR == kTRUE ) nEvents.push_back( 85368 );
//        else nEvents.push_back( 85351 );
//        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2_0001"+RocCorr+".root";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2226 );
        else nEvents.push_back( 2226 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _Test_EMu )
    {
        isMC = kTRUE;
        Type = "TEST";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WW_34"); Xsec.push_back( 1.0 ); Wsum.push_back( 10788.0 ); nEvents.push_back( 10788 );
        Location = "test/SelectedEMu_WW_34"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
}// end of Get()


vector<SelProc_t> LocalFileMgr::FindProc ( TString search, Bool_t notify, Bool_t instaGet )
{
    TString srch = search;
    srch.ToUpper();
    srch.ReplaceAll("-", "TO");
    srch.ReplaceAll(" TO", "TO"); srch.ReplaceAll("_TO", "TO");
    srch.ReplaceAll("TO ", "TO"); srch.ReplaceAll("TO_", "TO");
    srch.ReplaceAll("DRELLYAN", "DY");
    srch.ReplaceAll("ZW", "WZ");
    srch.ReplaceAll("SINGLETOP", "TW");
    srch.ReplaceAll("SINGLEANTITOP", "TBARW");
    srch.ReplaceAll("TOPANTITOP", "TTBAR");
    srch.ReplaceAll("TANTIT", "TTBAR");
    srch.ReplaceAll("DIMUON", "MUMU");
    srch.ReplaceAll("DIELECTRON", "EE");
    srch.ReplaceAll("DITAU", "TAUTAU");
    srch.ReplaceAll("BACKGROUND", "BKG");

    if ( notify == kTRUE ) cout << "Searched for: " << search << "\nFound: ";
    vector<SelProc_t> Result;
    SelProc_t first = _None, last = _None;
    if ( srch.Contains("DY") )
    {
        if ( srch.Contains("TAUTAU") )
        {
            if ( srch.Contains("MUMU") )
            {
                if ( srch.Contains("FULL") )
                {
                    Result.push_back(_MuMu_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_DYTauTau_Full] << "." << endl;
                }
                else
                {
                    // Checking for various intervals
                    if ( srch.Contains("INFTO") )
                        first = _EndOf_MuMu_DYTauTau_Normal;
                    else if ( srch.Contains("3000TO") )
                        first = _EndOf_MuMu_DYTauTau_Normal;
                    else if ( srch.Contains("2000TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1500TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1000TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("800TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("700TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("500TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("400TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("200TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("100TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("50TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("10TO") )
                        first = _MuMu_DYTauTau_10to50;

                    if ( srch.Contains("TOINF") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO3000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO2000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1500") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO800") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO700") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO500") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO400") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO200") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO100") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO50") )
                        last = _MuMu_DYTauTau_10to50;
                    else if ( srch.Contains("TO10") )
                        last = _EndOf_MuMu_MCsignal_Normal;

                    // Swapping first with last if necessary
                    if ( int(first)>int(last) && last!=_None )
                    {
                        SelProc_t NewLast = SelProc_t(int(first)-1);
                        first = SelProc_t(int(last)+1);
                        last = NewLast;
                    }
                    if ( first == _MuMu_DYTauTau_10to50 && last == _MuMu_DYTauTau_50toInf )
                    {
                        Result.push_back(_MuMu_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_MuMu_DYTauTau_Full] << "." << endl;
                    }
                    else if ( first != _None && last != _None )
                    {
                        for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                        Result.push_back(_MuMu_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_MuMu_DYTauTau_Full] << "." << endl;
                    }
                }
            }// end of if(MuMu)

            else if ( srch.Contains("EE") )
            {
                if ( srch.Contains("FULL") )
                {
                    Result.push_back(_EE_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_EE_DYTauTau_Full] << "." << endl;
                }
                else
                {
                    // Checking for various intervals
                    if ( srch.Contains("INFTO") )
                        first = _EndOf_EE_DYTauTau_Normal;
                    else if ( srch.Contains("3000TO") )
                        first = _EndOf_EE_DYTauTau_Normal;
                    else if ( srch.Contains("2000TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("1500TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("1000TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("800TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("700TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("500TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("400TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("200TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("100TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("50TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("10TO") )
                        first = _EE_DYTauTau_10to50;

                    if ( srch.Contains("TOINF") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO3000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO2000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1500") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO800") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO700") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO500") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO400") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO200") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO100") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("TO50") )
                        last = _EE_DYTauTau_10to50;
                    else if ( srch.Contains("TO10") )
                        last = _EndOf_EE_MCsignal_Normal;

                    // Swapping first with last if necessary
                    if ( int(first)>int(last) && last!=_None )
                    {
                        SelProc_t NewLast = SelProc_t(int(first)-1);
                        first = SelProc_t(int(last)+1);
                        last = NewLast;
                    }
                    if ( first == _EE_DYTauTau_10to50 && last == _EE_DYTauTau_50toInf )
                    {
                        Result.push_back(_EE_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_EE_DYTauTau_Full] << "." << endl;
                    }
                    else if ( first != _None && last != _None )
                    {
                        for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                        Result.push_back(_EE_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_EE_DYTauTau_Full] << "." << endl;
                    }
                }
            }// end of if(EE)

            else if ( srch.Contains("EMU") )
            {
                if ( srch.Contains("FULL") )
                {
                    Result.push_back(_EMu_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_DYTauTau_Full] << "." << endl;
                }
                else
                {
                    // Checking for various intervals
                    if ( srch.Contains("INFTO") )
                        first = _EndOf_EMu_DYTauTau_Normal;
                    else if ( srch.Contains("3000TO") )
                        first = _EndOf_EMu_DYTauTau_Normal;
                    else if ( srch.Contains("2000TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1500TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1000TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("800TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("700TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("500TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("400TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("200TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("100TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("50TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("10TO") )
                        first = _EMu_DYTauTau_10to50;

                    if ( srch.Contains("TOINF") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO3000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO2000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1500") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO1000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO800") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO700") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO500") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO400") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO200") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO100") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("TO50") )
                        last = _EMu_DYTauTau_10to50;
                    else if ( srch.Contains("TO10") )
                        last = _EndOf_EE;

                    // Swapping first with last if necessary
                    if ( int(first)>int(last) && last!=_None )
                    {
                        SelProc_t NewLast = SelProc_t(int(first)-1);
                        first = SelProc_t(int(last)+1);
                        last = NewLast;
                    }
                    if ( first == _EMu_DYTauTau_10to50 && last == _EMu_DYTauTau_50toInf )
                    {
                        Result.push_back(_EMu_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_EMu_DYTauTau_Full] << "." << endl;
                    }
                    else if ( first != _None && last != _None )
                    {
                        for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                        Result.push_back(_EMu_DYTauTau_Full);
                        if ( notify == kTRUE ) cout << Procname[_EMu_DYTauTau_Full] << "." << endl;
                    }
                }
            }// end of if(EMu)

        }// end of if(DYTauTau)

        else if ( srch.Contains("MUMU") )
        {
            // Checking for various intervals
            if ( srch.Contains("INFTO") )
                first = _EndOf_MuMu_MCsignal_Normal;
            else if ( srch.Contains("3000TO") )
                first = _EndOf_MuMu_MCsignal_Normal;
            else if ( srch.Contains("2000TO") )
                first =_MuMu_DY_2000to3000;
            else if ( srch.Contains("1500TO") )
                first =_MuMu_DY_1500to2000;
            else if ( srch.Contains("1000TO") )
                first =_MuMu_DY_1000to1500;
            else if ( srch.Contains("800TO") )
                first =_MuMu_DY_800to1000;
            else if ( srch.Contains("700TO") )
                first =_MuMu_DY_700to800;
            else if ( srch.Contains("500TO") )
                first =_MuMu_DY_500to700;
            else if ( srch.Contains("400TO") )
                first =_MuMu_DY_400to500;
            else if ( srch.Contains("200TO") )
                first =_MuMu_DY_200to400;
            else if (  srch.Contains("100TO") )
                first =_MuMu_DY_100to200;
            else if ( srch.Contains("50TO") )
                first =_MuMu_DY_50to100;
            else if ( srch.Contains("10TO") )
                first = _MuMu_DY_10to50;

            if ( srch.Contains("TOINF") )
                last =_MuMu_DY_2000to3000;
            else if ( srch.Contains("TO3000") )
                last =_MuMu_DY_2000to3000;
            else if ( srch.Contains("TO2000") )
                last =_MuMu_DY_1500to2000;
            else if ( srch.Contains("TO1500") )
                last =_MuMu_DY_1000to1500;
            else if ( srch.Contains("TO1000") )
                last =_MuMu_DY_800to1000;
            else if ( srch.Contains("TO800") )
                last =_MuMu_DY_700to800;
            else if ( srch.Contains("TO700") )
                last =_MuMu_DY_500to700;
            else if ( srch.Contains("TO500") )
                last =_MuMu_DY_400to500;
            else if ( srch.Contains("TO400") )
                last =_MuMu_DY_200to400;
            else if ( srch.Contains("TO200") )
                last =_MuMu_DY_100to200;
            else if ( srch.Contains("TO100") )
                last =_MuMu_DY_50to100;
            else if ( srch.Contains("TO50") )
                last = _MuMu_DY_10to50;
            else if ( srch.Contains("TO10") )
                last = _None;

            // Swapping first with last if necessary
            if ( int(first)>int(last) )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _MuMu_DY_10to50 && last ==_MuMu_DY_2000to3000 )
            {
                Result.push_back(_MuMu_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_DY_Full] << "." << endl;
            }
            else if ( first != _None && last != _None )
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_MuMu_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_DY_Full] << "." << endl;
            }

        }// end of if(DYMuMu)

        else if ( srch.Contains("EE") )
        {
            // Checking for various intervals
            if ( srch.Contains("INFTO") )
                first = _EndOf_EE_MCsignal_Normal;
            else if ( srch.Contains("3000TO") )
                first = _EndOf_EE_MCsignal_Normal;
            else if ( srch.Contains("2000TO") )
                first = _EE_DY_2000to3000;
            else if ( srch.Contains("1500TO") )
                first = _EE_DY_1500to2000;
            else if ( srch.Contains("1000TO") )
                first = _EE_DY_1000to1500;
            else if ( srch.Contains("800TO") )
                first = _EE_DY_800to1000;
            else if ( srch.Contains("700TO") )
                first = _EE_DY_700to800;
            else if ( srch.Contains("500TO") )
                first = _EE_DY_500to700;
            else if ( srch.Contains("400TO") )
                first = _EE_DY_400to500;
            else if ( srch.Contains("200TO") )
                first = _EE_DY_200to400;
            else if ( srch.Contains("100TO") )
                first = _EE_DY_100to200;
            else if ( srch.Contains("50TO") )
                first = _EE_DY_50to100;
            else if ( srch.Contains("10TO") )
                first = _EE_DY_10to50;

            if ( srch.Contains("TOINF") )
                last = _EE_DY_2000to3000;
            else if ( srch.Contains("TO3000") )
                last = _EE_DY_2000to3000;
            else if ( srch.Contains("TO2000") )
                last = _EE_DY_1500to2000;
            else if ( srch.Contains("TO1500") )
                last = _EE_DY_1000to1500;
            else if ( srch.Contains("TO1000") )
                last = _EE_DY_800to1000;
            else if ( srch.Contains("TO800") )
                last = _EE_DY_700to800;
            else if ( srch.Contains("TO700") )
                last = _EE_DY_500to700;
            else if ( srch.Contains("TO500") )
                last = _EE_DY_400to500;
            else if ( srch.Contains("TO400") )
                last = _EE_DY_200to400;
            else if ( srch.Contains("TO200") )
                last = _EE_DY_100to200;
            else if ( srch.Contains("TO100") )
                last = _EE_DY_50to100;
            else if ( srch.Contains("TO50") )
                last = _EE_DY_10to50;
            else if ( srch.Contains("TO10") )
                last = _EndOf_MuMu;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _EE_DY_10to50 && last == _EE_DY_2000to3000 )
            {
                Result.push_back(_EE_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DY_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_EE_DY_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DY_Full] << "." << endl;
            }

        }// end of if(DYEE)       

    }// end of if(DrellYan)

    else if ( srch.Contains("TTBAR") )
    {
        if ( srch.Contains("MUMU") )
        {
            if ( srch.Contains("FULL") )
            {
                Result.push_back(_MuMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700TO") )
                first = _MuMu_ttbar_700to1000;
            else if ( srch.Contains("1000TO") )
                first = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("INFTO") )
                first = _EndOf_MuMu_ttbar_Normal;
            else first = _MuMu_ttbar;

            if ( srch.Contains("TO700") )
                last = _MuMu_ttbar;
            else if ( srch.Contains("TO1000") )
                last = _MuMu_ttbar_700to1000;
            else if ( srch.Contains("TO1500") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("TO2000") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("TO3000") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("TOINF") )
                last = _MuMu_ttbar_1000toInf;
            else last = _MuMu_ttbar;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _MuMu_ttbar && last == _MuMu_ttbar_1000toInf )
            {
                Result.push_back(_MuMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_ttbar_Full] << "." << endl;
            }
            else if ( first != _None && last != _None )
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_MuMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_ttbar_Full] << "." << endl;
            }
        }// end of if(MuMu)

        else if ( srch.Contains("EE") )
        {
            if ( srch.Contains("FULL") )
            {
                Result.push_back(_EE_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700TO") )
                first = _EE_ttbar_700to1000;
            else if ( srch.Contains("1000TO") )
                first = _EE_ttbar_1000toInf;
            else if ( srch.Contains("INFTO") )
                first = _EndOf_EE_ttbar_Normal;
            else first = _EE_ttbar;

            if ( srch.Contains("TO700"))
                last = _EE_ttbar;
            else if ( srch.Contains("TO1000") )
                last = _EE_ttbar_700to1000;
            else if ( srch.Contains("TO1500") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("TO2000") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("TO3000") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("TOINF") )
                last = _EE_ttbar_1000toInf;
            else last = _EE_ttbar;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _EE_ttbar && last == _EE_ttbar_1000toInf )
            {
                Result.push_back(_EE_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_ttbar_Full] << "." << endl;
            }
            else if ( first != _None && last != _None )
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_EE_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_ttbar_Full] << "." << endl;
            }
        }// end of if(EE)

        else if ( srch.Contains("EMU") )
        {
            if ( srch.Contains("FULL") )
            {
                Result.push_back(_EMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EMu_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700TO") )
                first = _EMu_ttbar_700to1000;
            else if ( srch.Contains("1000TO") )
                first = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("INFTO") )
                first = _EndOf_EMu_ttbar_Normal;
            else first = _EMu_ttbar;

            if ( srch.Contains("TO700") )
                last = _EMu_ttbar;
            else if ( srch.Contains("TO1000") )
                last = _EMu_ttbar_700to1000;
            else if ( srch.Contains("TO1500") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("TO2000") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("TO3000") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("TOINF") )
                last = _EMu_ttbar_1000toInf;
            else last = _EMu_ttbar;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _EMu_ttbar && last == _EMu_ttbar_1000toInf )
            {
                Result.push_back(_EMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EMu_ttbar_Full] << "." << endl;
            }
            else if ( first != _None && last != _None )
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_EMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EMu_ttbar_Full] << "." << endl;
            }
        }// end of if(EMu)

    }// end of if(ttbar)

    else if ( srch.Contains("TW") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_tW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_tW] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_tW);
            if ( notify == kTRUE ) cout << Procname[_EE_tW] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_tW);
            if ( notify == kTRUE ) cout << Procname[_EMu_tW] << "." << endl;
        }

    }

    else if ( srch.Contains("TBARW") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_tbarW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_tbarW] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_tbarW);
            if ( notify == kTRUE ) cout << Procname[_EE_tbarW] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_tbarW);
            if ( notify == kTRUE ) cout << Procname[_EMu_tbarW] << "." << endl;
        }
    }

    else if ( srch.Contains("ZZ") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_ZZ);
            if ( notify == kTRUE ) cout << Procname[_MuMu_ZZ] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_ZZ);
            if ( notify == kTRUE ) cout << Procname[_EE_ZZ] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_ZZ);
            if ( notify == kTRUE ) cout << Procname[_EMu_ZZ] << "." << endl;
        }
    }

    else if ( srch.Contains("WZ") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_WZ);
            if ( notify == kTRUE ) cout << Procname[_MuMu_WZ] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_WZ);
            if ( notify == kTRUE ) cout << Procname[_EE_WZ] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_WZ);
            if ( notify == kTRUE ) cout << Procname[_EMu_WZ] << "." << endl;
        }
    }

    else if ( srch.Contains("WW") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_WW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_WW] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_WW);
            if ( notify == kTRUE ) cout << Procname[_EE_WW] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_WW);
            if ( notify == kTRUE ) cout << Procname[_EMu_WW] << "." << endl;
        }
    }

    else if ( srch.Contains("DIBOSON") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_ZZ);
            Result.push_back(_MuMu_WZ);
            Result.push_back(_MuMu_WW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_ZZ] << ", " << Procname[_MuMu_WZ] << ", " << Procname[_MuMu_WW] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_ZZ);
            Result.push_back(_EE_WZ);
            Result.push_back(_EE_WW);
            if ( notify == kTRUE ) cout << Procname[_EE_ZZ] << ", " << Procname[_EE_WZ] << ", " << Procname[_EE_WW] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_ZZ);
            Result.push_back(_EMu_WZ);
            Result.push_back(_EMu_WW);
            if ( notify == kTRUE ) cout << Procname[_EMu_ZZ] << ", " << Procname[_EMu_WZ] << ", " << Procname[_EMu_WW] << "." << endl;
        }
    }

    else if ( srch.Contains("VVNST") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_VVnST);
            if ( notify == kTRUE ) cout << Procname[_MuMu_VVnST] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_VVnST);
            if ( notify == kTRUE ) cout << Procname[_EE_VVnST] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_VVnST);
            if ( notify == kTRUE ) cout << Procname[_EMu_VVnST] << "." << endl;
        }
    }

    else if ( srch.Contains("WJETS") || srch.Contains("W+JETS") )
    {
        if ( srch.Contains("MUMU") )
        {
            if (srch.Contains("FULL"))
            {
                Result.push_back(_MuMu_WJets_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_WJets_Full] << "." << endl;
            }
            if (srch.Contains("EXT"))
            {
                Result.push_back(_MuMu_WJets_ext2v5);
                if ( notify == kTRUE ) cout << Procname[_MuMu_WJets_ext2v5] << "." << endl;
            }
            else
            {
                Result.push_back(_MuMu_WJets);
                if ( notify == kTRUE ) cout << Procname[_MuMu_WJets] << "." << endl;
            }
        }
        else if ( srch.Contains("EE") )
        {
            if (srch.Contains("FULL"))
            {
                Result.push_back(_EE_WJets_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_WJets_Full] << "." << endl;
            }
            if (srch.Contains("EXT"))
            {
                Result.push_back(_EE_WJets_ext2v5);
                if ( notify == kTRUE ) cout << Procname[_EE_WJets_ext2v5] << "." << endl;
            }
            else
            {
                Result.push_back(_EE_WJets);
                if ( notify == kTRUE ) cout << Procname[_EE_WJets] << "." << endl;
            }
        }
        else if ( srch.Contains("EMU") )
        {
            if (srch.Contains("FULL"))
            {
                Result.push_back(_EMu_WJets_Full);
                if ( notify == kTRUE ) cout << Procname[_EMu_WJets_Full] << "." << endl;
            }
            if (srch.Contains("EXT"))
            {
                Result.push_back(_EMu_WJets_ext2v5);
                if ( notify == kTRUE ) cout << Procname[_EMu_WJets_ext2v5] << "." << endl;
            }
            else
            {
                Result.push_back(_EMu_WJets);
                if ( notify == kTRUE ) cout << Procname[_EMu_WJets] << "." << endl;
            }
        }
    }

    else if ( srch.Contains("QCD") )
    {
        if ( srch.Contains("Mu") || srch.Contains("mu") || srch.Contains("MU") )
        {
            // Checking for various intervals
            if ( srch.Contains("INFTO") )
                first = _EndOf_MuMu_QCDMuEnriched_Normal;
            else if ( srch.Contains("1000TO") )
                first = _MuMu_QCDMuEnriched_1000toInf;
            else if ( srch.Contains("800TO") )
                first = _MuMu_QCDMuEnriched_800to1000;
            else if ( srch.Contains("600TO") )
                first = _MuMu_QCDMuEnriched_600to800;
            else if ( srch.Contains("470TO") )
                first = _MuMu_QCDMuEnriched_470to600;
            else if ( srch.Contains("300TO") )
                first = _MuMu_QCDMuEnriched_300to470;
            else if ( srch.Contains("170TO") )
                first = _MuMu_QCDMuEnriched_170to300;
            else if ( srch.Contains("120TO") )
                first = _MuMu_QCDMuEnriched_120to170;
            else if ( srch.Contains("80TO") )
                first = _MuMu_QCDMuEnriched_80to120;
            else if ( srch.Contains("50TO") )
                first = _MuMu_QCDMuEnriched_50to80;
            else if ( srch.Contains("30TO") )
                first = _MuMu_QCDMuEnriched_30to50;
            else if ( srch.Contains("20TO") )
                first = _MuMu_QCDMuEnriched_20to30;
            else if ( srch.Contains("15TO") )
                first = _MuMu_QCDMuEnriched_15to20;

            if ( srch.Contains("TOINF") )
                last = _MuMu_QCDMuEnriched_1000toInf;
            else if ( srch.Contains("TO1000") )
                last = _MuMu_QCDMuEnriched_800to1000;
            else if ( srch.Contains("TO800") )
                last = _MuMu_QCDMuEnriched_600to800;
            else if ( srch.Contains("TO600") )
                last = _MuMu_QCDMuEnriched_470to600;
            else if ( srch.Contains("TO470") )
                last = _MuMu_QCDMuEnriched_300to470;
            else if ( srch.Contains("TO300") )
                last = _MuMu_QCDMuEnriched_170to300;
            else if ( srch.Contains("TO170") )
                last = _MuMu_QCDMuEnriched_120to170;
            else if ( srch.Contains("TO120") )
                last = _MuMu_QCDMuEnriched_80to120;
            else if ( srch.Contains("TO80") )
                last = _MuMu_QCDMuEnriched_50to80;
            else if ( srch.Contains("TO50") )
                last = _MuMu_QCDMuEnriched_30to50;
            else if ( srch.Contains("TO30") )
                last = _MuMu_QCDMuEnriched_20to30;
            else if ( srch.Contains("TO20") )
                last = _MuMu_QCDMuEnriched_15to20;
            else if ( srch.Contains("TO15") )
                last = _EndOf_MuMu_WJets_Normal;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _MuMu_QCDMuEnriched_15to20 && last == _MuMu_QCDMuEnriched_1000toInf )
            {
                Result.push_back(_MuMu_QCDMuEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_QCDMuEnriched_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_MuMu_QCDMuEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_QCDMuEnriched_Full] << "." << endl;
            }

        }// end of if(MuEnriched)
        else if ( srch.Contains("EM") || srch.Contains("EE") )
        {
            // Checking for various intervals
            if ( srch.Contains("INFTO") )
                first = _EndOf_EE_QCDEMEnriched_Normal;
            else if ( srch.Contains("300TO") )
                first = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("170TO") )
                first = _EE_QCDEMEnriched_170to300;
            else if ( srch.Contains("120TO") )
                first = _EE_QCDEMEnriched_120to170;
            else if ( srch.Contains("80TO") )
                first = _EE_QCDEMEnriched_80to120;
            else if ( srch.Contains("50TO") )
                first = _EE_QCDEMEnriched_50to80;
            else if ( srch.Contains("30TO") )
                first = _EE_QCDEMEnriched_30to50;
            else if ( srch.Contains("20TO") )
                first = _EE_QCDEMEnriched_20to30;

            if ( srch.Contains("TOINF") )
                last = _EE_QCDEMEnriched_300toInf;
            else if (srch.Contains("TO1000") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("TO800") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("TO600") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("TO470") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("TO300") )
                last = _EE_QCDEMEnriched_170to300;
            else if ( srch.Contains("TO170") )
                last = _EE_QCDEMEnriched_120to170;
            else if ( srch.Contains("TO120") )
                last = _EE_QCDEMEnriched_80to120;
            else if ( srch.Contains("TO80") )
                last = _EE_QCDEMEnriched_50to80;
            else if ( srch.Contains("TO50") )
                last = _EE_QCDEMEnriched_30to50;
            else if ( srch.Contains("TO30") )
                last = _EE_QCDEMEnriched_20to30;
            else if ( srch.Contains("TO20") )
                last = _EndOf_EE_WJets_Normal;

            // Swapping first with last if necessary
            if ( int(first)>int(last)  && last!=_None)
            {
                SelProc_t NewLast = SelProc_t(int(first)-1);
                first = SelProc_t(int(last)+1);
                last = NewLast;
            }
            if ( first == _EE_QCDEMEnriched_20to30 && last == _EE_QCDEMEnriched_300toInf )
            {
                Result.push_back(_EE_QCDEMEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_QCDEMEnriched_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                Result.push_back(_EE_QCDEMEnriched_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_QCDEMEnriched_Full] << "." << endl;
            }

        }// end of if(EMEnriched)

    }// end of if(QCD)

    else if ( srch.Contains("DOUBLE") )
    {
        if ( srch.Contains("EG") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("BTO") )
                first = _EE_DoubleEG_B;
            else if ( srch.Contains("CTO") )
                first = _EE_DoubleEG_C;
            else if ( srch.Contains("DTO") )
                first = _EE_DoubleEG_D;
            else if ( srch.Contains("ETO") )
                first = _EE_DoubleEG_E;
            else if ( srch.Contains("FTO") )
                first = _EE_DoubleEG_F;
            else if ( srch.Contains("GTO") )
                first = _EE_DoubleEG_G;
            else if ( srch.Contains("HTO") )
                first = _EE_DoubleEG_H;

            if ( srch.Contains("TOB") )
                last = _EE_DoubleEG_B;
            else if ( srch.Contains("TOC") )
                last = _EE_DoubleEG_C;
            else if ( srch.Contains("TOD") )
                last = _EE_DoubleEG_D;
            else if ( srch.Contains("TOE") )
                last = _EE_DoubleEG_E;
            else if ( srch.Contains("TOF") )
                last = _EE_DoubleEG_F;
            else if ( srch.Contains("TOG") )
                last = _EE_DoubleEG_G;
            else if ( srch.Contains("TOH") )
                last = _EE_DoubleEG_H;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if ( first == _EE_DoubleEG_B && last == _EE_DoubleEG_H )
            {
                Result.push_back(_EE_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
            else if ( srch.Contains("FULL") )
            {
                Result.push_back(_EE_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("RUNB") )
            {
                Result.push_back(_EE_DoubleEG_B);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("RUNC") )
            {
                Result.push_back(_EE_DoubleEG_C);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("RUND") )
            {
                if ( srch.Contains("_E") || srch.Contains("RUNE") )
                {
                    Result.push_back(_EE_DoubleEG_E);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("RUNF") )
                {
                    Result.push_back(_EE_DoubleEG_F);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("RUNG") )
                {
                    Result.push_back(_EE_DoubleEG_G);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("RUNH") )
                {
                    Result.push_back(_EE_DoubleEG_H);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_H] << "." << endl;
                }
                else
                {
                    Result.push_back(_EE_DoubleEG_D);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_D] << "." << endl;
                }
            }            
            else
            {
                Result.push_back(_EE_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_Full] << "." << endl;
            }

        }// end of if (EG)

    }// end of if(Double)

    else if ( srch.Contains("SINGLE") )
    {
        if ( srch.Contains("MUON") )
        {
            if ( srch.Contains("MUMU") )
            {
                // Checking if search contains intervals
                if ( srch.Contains("BTO") )
                    first = _MuMu_SingleMuon_B;
                else if ( srch.Contains("CTO") )
                    first = _MuMu_SingleMuon_C;
                else if ( srch.Contains("DTO") )
                    first = _MuMu_SingleMuon_D;
                else if ( srch.Contains("ETO") )
                    first = _MuMu_SingleMuon_E;
                else if ( srch.Contains("FTO") )
                    first = _MuMu_SingleMuon_F;
                else if ( srch.Contains("GTO") )
                    first = _MuMu_SingleMuon_G;
                else if ( srch.Contains("HTO") )
                    first = _MuMu_SingleMuon_H;

                if ( srch.Contains("TOB") )
                    last = _MuMu_SingleMuon_B;
                else if ( srch.Contains("TOC") )
                    last = _MuMu_SingleMuon_C;
                else if ( srch.Contains("TOD") )
                    last = _MuMu_SingleMuon_D;
                else if ( srch.Contains("TOE") )
                    last = _MuMu_SingleMuon_E;
                else if ( srch.Contains("TOF") )
                    last = _MuMu_SingleMuon_F;
                else if ( srch.Contains("TOG") )
                    last = _MuMu_SingleMuon_G;
                else if ( srch.Contains("TOH") )
                    last = _MuMu_SingleMuon_H;

                // Swapping first with last if necessary
                if ( int(first)>int(last) && last!=_None )
                {
                    SelProc_t NewLast = first;
                    first = last;
                    last = NewLast;
                }
                if ( first == _MuMu_SingleMuon_B && last == _MuMu_SingleMuon_H )
                {
                    Result.push_back(_MuMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_Full] << "." << endl;
                }
                else if ( first != _None && last != _None)
                {
                    for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                else if ( srch.Contains("FULL") )
                {
                    Result.push_back(_MuMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_Full] << "." << endl;
                }
                else if ( srch.Contains("_B") || srch.Contains("RUNB") )
                {
                    Result.push_back(_MuMu_SingleMuon_B);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_B] << "." << endl;
                }
                else if ( srch.Contains("_C") || srch.Contains("RUNC") )
                {
                    Result.push_back(_MuMu_SingleMuon_C);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_C] << "." << endl;
                }
                else if ( srch.Contains("_D") || srch.Contains("RUND") )
                {
                    Result.push_back(_MuMu_SingleMuon_D);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_D] << "." << endl;
                }
                else if ( srch.Contains("_E") || srch.Contains("RUNE") )
                {
                    Result.push_back(_MuMu_SingleMuon_E);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("RUNF") )
                {
                    Result.push_back(_MuMu_SingleMuon_F);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("RUNG") )
                {
                    Result.push_back(_MuMu_SingleMuon_G);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("RUNH") )
                {
                    Result.push_back(_MuMu_SingleMuon_H);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_H] << "." << endl;
                }
                else
                {
                    Result.push_back(_MuMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_Full] << "." << endl;
                }
            }// end of if(MuMu)

            else if ( srch.Contains("EMU") )
            {
                // Checking if search contains intervals
                if ( srch.Contains("BTO") )
                    first = _EMu_SingleMuon_B;
                else if ( srch.Contains("CTO") )
                    first = _EMu_SingleMuon_C;
                else if ( srch.Contains("DTO") )
                    first = _EMu_SingleMuon_D;
                else if ( srch.Contains("ETO") )
                    first = _EMu_SingleMuon_E;
                else if ( srch.Contains("FTO") )
                    first = _EMu_SingleMuon_F;
                else if ( srch.Contains("GTO") )
                    first = _EMu_SingleMuon_G;
                else if ( srch.Contains("HTO") )
                    first = _EMu_SingleMuon_H;

                if ( srch.Contains("TOB") )
                    last = _EMu_SingleMuon_B;
                else if ( srch.Contains("TOC") )
                    last = _EMu_SingleMuon_C;
                else if ( srch.Contains("TOD") )
                    last = _EMu_SingleMuon_D;
                else if ( srch.Contains("TOE") )
                    last = _EMu_SingleMuon_E;
                else if ( srch.Contains("TOF") )
                    last = _EMu_SingleMuon_F;
                else if ( srch.Contains("TOG") )
                    last = _EMu_SingleMuon_G;
                else if ( srch.Contains("TOH") )
                    last = _EMu_SingleMuon_H;

                // Swapping first with last if necessary
                if ( int(first)>int(last) && last!=_None )
                {
                    SelProc_t NewLast = first;
                    first = last;
                    last = NewLast;
                }
                if ( first == _EMu_SingleMuon_B && last == _EMu_SingleMuon_H )
                {
                    Result.push_back(_EMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_Full] << "." << endl;
                }
                else if ( first != _None && last != _None)
                {
                    for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
                else if ( srch.Contains("FULL") )
                {
                    Result.push_back(_EMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_Full] << "." << endl;
                }
                else if ( srch.Contains("_B") || srch.Contains("RUNB") )
                {
                    Result.push_back(_EMu_SingleMuon_B);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_B] << "." << endl;
                }
                else if ( srch.Contains("_C") || srch.Contains("RUNC") )
                {
                    Result.push_back(_EMu_SingleMuon_C);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_C] << "." << endl;
                }
                else if ( srch.Contains("_D") || srch.Contains("RUND") )
                {
                    Result.push_back(_EMu_SingleMuon_D);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_D] << "." << endl;
                }
                else if ( srch.Contains("_E") || srch.Contains("RUNE") )
                {
                    Result.push_back(_EMu_SingleMuon_E);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("RUNF") )
                {
                    Result.push_back(_EMu_SingleMuon_F);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("RUNG") )
                {
                    Result.push_back(_EMu_SingleMuon_G);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("RUNH") )
                {
                    Result.push_back(_EMu_SingleMuon_H);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_H] << "." << endl;
                }
                else
                {
                    Result.push_back(_EMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_Full] << "." << endl;
                }
            }// end of if(EMu)

        }// end of if(SingleMuon)

        if ( srch.Contains("ELEC") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("BTO") )
                first = _EE_SingleElectron_B;
            else if ( srch.Contains("CTO") )
                first = _EE_SingleElectron_C;
            else if ( srch.Contains("DTO") )
                first = _EE_SingleElectron_D;
            else if ( srch.Contains("ETO") )
                first = _EE_SingleElectron_E;
            else if ( srch.Contains("FTO") )
                first = _EE_SingleElectron_F;
            else if ( srch.Contains("GTO") )
                first = _EE_SingleElectron_G;
            else if ( srch.Contains("HTO") )
                first = _EE_SingleElectron_H;

            if ( srch.Contains("TOB") )
                last = _EE_SingleElectron_B;
            else if ( srch.Contains("TOC") )
                last = _EE_SingleElectron_C;
            else if ( srch.Contains("TOD") )
                last = _EE_SingleElectron_D;
            else if ( srch.Contains("TOE") )
                last = _EE_SingleElectron_E;
            else if ( srch.Contains("TOF"))
                last = _EE_SingleElectron_F;
            else if ( srch.Contains("TOG") )
                last = _EE_SingleElectron_G;
            else if ( srch.Contains("TOH") )
                last = _EE_SingleElectron_H;

            // Swapping first with last if necessary
            if ( int(first)>int(last) && last!=_None )
            {
                SelProc_t NewLast = first;
                first = last;
                last = NewLast;
            }
            if ( first == _EE_SingleElectron_B && last == _EE_SingleElectron_H )
            {
                Result.push_back(_EE_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_Full] << "." << endl;
            }
            else if ( first != _None && last != _None)
            {
                for ( SelProc_t pr=first; pr<=last; pr=next(pr) )
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
            else if ( srch.Contains("FULL") )
            {
                Result.push_back(_EE_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("RUNB") )
            {
                Result.push_back(_EE_SingleElectron_B);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("RUNC") )
            {
                Result.push_back(_EE_SingleElectron_C);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("RUND") )
            {
                Result.push_back(_EE_SingleElectron_D);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_D] << "." << endl;
            }
            else if ( srch.Contains("_E") || srch.Contains("RUNE") )
            {
                Result.push_back(_EE_SingleElectron_E);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_E] << "." << endl;
            }
            else if ( srch.Contains("_F") || srch.Contains("RUNF") )
            {
                Result.push_back(_EE_SingleElectron_F);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_F] << "." << endl;
            }
            else if ( srch.Contains("_G") || srch.Contains("RUNG") )
            {
                Result.push_back(_EE_SingleElectron_G);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_G] << "." << endl;
            }
            else if ( srch.Contains("_H") || srch.Contains("RUNH") )
            {
                Result.push_back(_EE_SingleElectron_H);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_H] << "." << endl;
            }
            else
            {
                Result.push_back(_EE_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_Full] << "." << endl;
            }
        }// end of if(SingleElectron)

    }// end of if(Single)

    else if ( srch.Contains("SIGNAL") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_DY_Full);
            if ( notify == kTRUE ) cout << Procname[_MuMu_DY_Full] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_DY_Full);
            if ( notify == kTRUE ) cout << Procname[_EE_DY_Full] << "." << endl;
        }
    }

    else if ( srch.Contains("BKG") )
    {
        if ( srch.Contains("MUMU") )
        {
            Result.push_back(_MuMu_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_MuMu_Bkg_Full] << "." << endl;
        }
        else if ( srch.Contains("EE") )
        {
            Result.push_back(_EE_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_EE_Bkg_Full] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
        {
            Result.push_back(_EMu_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_EMu_Bkg_Full] << "." << endl;
        }
    }

    else if ( srch.Contains("TEST") )
    {
        if ( srch.Contains("EE") )
        {
            Result.push_back(_Test_EE);
            if ( notify == kTRUE )  cout << Procname[_Test_EE] << "." << endl;
        }
        else if ( srch.Contains("MUMU") )
        {
            Result.push_back(_Test_MuMu);
            if ( notify == kTRUE )  cout << Procname[_Test_MuMu] << "." << endl;
        }
        else if ( srch.Contains("EMU") )
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

    if ( instaGet == kTRUE )    // Applying SetProc if asked
    {
        cout << "Getting information for returned processes..";
        for ( Int_t i=0; i<(int(Result.size())); i++ )
        {
            if ( Result[i] != _None ) this->SetProc( Result[i], kFALSE );
        }
        cout << " Finished." << endl;
    }

    return Result;

}// end of FindProc()


void LocalFileMgr::SwitchROCCORR()
{
    if ( ROCCORR == kTRUE )
    {
        ROCCORR = kFALSE;
        cout << "ROCCOR has been changed from kTRUE to kFALSE." << endl;
        return;
    }
    else
    {
        ROCCORR = kTRUE;
        cout << "ROCCOR has been changed from kFALSE to kTRUE." << endl;
        return;
    }
}


void LocalFileMgr::PrepareProcNames ()
{
    Procname[_None] = "None";
    Procname[_MuMu_DY_10to50] = "SelectedMuMu_DY_10to50";
    Procname[_MuMu_DY_50to100] = "SelectedMuMu_DY_50to100";
    Procname[_MuMu_DY_100to200] = "SelectedMuMu_DY_100to200";
    Procname[_MuMu_DY_200to400] = "SelectedMuMu_DY_200to400";
    Procname[_MuMu_DY_400to500] = "SelectedMuMu_DY_400to500";
    Procname[_MuMu_DY_500to700] = "SelectedMuMu_DY_500to700";
    Procname[_MuMu_DY_700to800] = "SelectedMuMu_DY_700to800";
    Procname[_MuMu_DY_800to1000] = "SelectedMuMu_DY_800to1000";
    Procname[_MuMu_DY_1000to1500] = "SelectedMuMu_DY_1000to1500";
    Procname[_MuMu_DY_1500to2000] = "SelectedMuMu_DY_1500to2000";
    Procname[_MuMu_DY_2000to3000] = "SelectedMuMu_DY_2000to3000";
    Procname[_EndOf_MuMu_MCsignal_Normal] = "EndOf_SelectedMuMu_MCsignal_Normal";
    Procname[_MuMu_DYTauTau_10to50] = "SelectedMuMu_DYTauTau_10to50";
    Procname[_MuMu_DYTauTau_50toInf] = "SelectedMuMu_DYTauTau_50toInf";
    Procname[_EndOf_MuMu_DYTauTau_Normal] = "EndOf_SelectedMuMu_DYTauTau_Normal";
    Procname[_MuMu_ttbar] = "SelectedMuMu_ttbar";
    Procname[_MuMu_ttbar_700to1000] = "SelectedMuMu_ttbar_700to1000";
    Procname[_MuMu_ttbar_1000toInf] = "SelectedMuMu_ttbar_1000toInf";
    Procname[_EndOf_MuMu_ttbar_Normal] = "EndOf_SelectedMuMu_ttbar_Normal";
    Procname[_MuMu_tW] = "SelectedMuMu_tW";
    Procname[_MuMu_tbarW] = "SelectedMuMu_tbarW";
    Procname[_MuMu_ZZ] = "SelectedMuMu_ZZ";
    Procname[_MuMu_WZ] = "SelectedMuMu_WZ";
    Procname[_MuMu_WW] = "SelectedMuMu_WW";
    Procname[_EndOf_MuMu_VVnST_Normal] = "EndOf_SelectedMuMu_VVnST_Normal";
    Procname[_MuMu_WJets] = "SelectedMuMu_WJets";
    Procname[_MuMu_WJets_ext2v5] = "SelectedMuMu_WJets_ext2v5";
    Procname[_EndOf_MuMu_WJets_Normal] = "EndOf_SelectedMuMu_WJets";
    Procname[_MuMu_QCDMuEnriched_15to20] = "SelectedMuMu_QCDMuEnriched_15to20";
    Procname[_MuMu_QCDMuEnriched_20to30] = "SelectedMuMu_QCDMuEnriched_20to30";
    Procname[_MuMu_QCDMuEnriched_30to50] = "SelectedMuMu_QCDMuEnriched_30to50";
    Procname[_MuMu_QCDMuEnriched_50to80] = "SelectedMuMu_QCDMuEnriched_50to80";
    Procname[_MuMu_QCDMuEnriched_80to120] = "SelectedMuMu_QCDMuEnriched_80to120";
    Procname[_MuMu_QCDMuEnriched_120to170] = "SelectedMuMu_QCDMuEnriched_120to170";
    Procname[_MuMu_QCDMuEnriched_170to300] = "SelectedMuMu_QCDMuEnriched_170to300";
    Procname[_MuMu_QCDMuEnriched_300to470] = "SelectedMuMu_QCDMuEnriched_300to470";
    Procname[_MuMu_QCDMuEnriched_470to600] = "SelectedMuMu_QCDMuEnriched_470to600";
    Procname[_MuMu_QCDMuEnriched_600to800] = "SelectedMuMu_QCDMuEnriched_600to800";
    Procname[_MuMu_QCDMuEnriched_800to1000] = "SelectedMuMu_QCDMuEnriched_800to1000";
    Procname[_MuMu_QCDMuEnriched_1000toInf] = "SelectedMuMu_QCDMuEnriched_1000toInf";
    Procname[_EndOf_MuMu_QCDMuEnriched_Normal] = "EndOf_SelectedMuMu_QCDMuEnriched_Normal";
    Procname[_EndOf_MuMu_MCbkg_Normal] = "EndOf_SelectedMuMu_MCbkg_Normal";
    Procname[_MuMu_SingleMuon_B] = "SelectedMuMu_SingleMuon_B";
    Procname[_MuMu_SingleMuon_C] = "SelectedMuMu_SingleMuon_C";
    Procname[_MuMu_SingleMuon_D] = "SelectedMuMu_SingleMuon_D";
    Procname[_MuMu_SingleMuon_E] = "SelectedMuMu_SingleMuon_E";
    Procname[_MuMu_SingleMuon_F] = "SelectedMuMu_SingleMuon_F";
    Procname[_MuMu_SingleMuon_G] = "SelectedMuMu_SingleMuon_G";
    Procname[_MuMu_SingleMuon_H] = "SelectedMuMu_SingleMuon_H";
    Procname[_EndOf_MuMu_Data_Normal] = "EndOf_SelectedMuMu_Data_Normal";
    Procname[_MuMu_DY_Full] = "SelectedMuMu_DY_Full";
    Procname[_MuMu_DYTauTau_Full] = "SelectedMuMu_DYTauTau_Full";
    Procname[_MuMu_ttbar_Full] = "SelectedMuMu_ttbar_Full";
    Procname[_MuMu_VVnST] = "SelectedMuMu_VVnST";
    Procname[_MuMu_WJets_Full] = "SelectedMuMu_WJets_Full";
    Procname[_MuMu_QCDMuEnriched_Full] = "SelectedMuMu_QCDMuEnriched_Full";
    Procname[_MuMu_Bkg_Full] = "SelectedMuMu_Bkg_Full";
    Procname[_MuMu_SingleMuon_Full] = "SelectedMuMu_SingleMuon_Full";
    Procname[_EndOf_MuMu_Special] = "EndOf_SelectedMuMu_Special";
    Procname[_Test_MuMu] = "Test_SelectedMuMu";
    Procname[_EndOf_MuMu] = "EndOf_SelectedMuMu";

    Procname[_EE_DY_10to50] = "SelectedEE_DY_10to50";
    Procname[_EE_DY_50to100] = "SelectedEE_DY_50to100";
    Procname[_EE_DY_100to200] = "SelectedEE_DY_100to200";
    Procname[_EE_DY_200to400] = "SelectedEE_DY_200to400";
    Procname[_EE_DY_400to500] = "SelectedEE_DY_400to500";
    Procname[_EE_DY_500to700] = "SelectedEE_DY_500to700";
    Procname[_EE_DY_700to800] = "SelectedEE_DY_700to800";
    Procname[_EE_DY_800to1000] = "SelectedEE_DY_800to1000";
    Procname[_EE_DY_1000to1500] = "SelectedEE_DY_1000to1500";
    Procname[_EE_DY_1500to2000] = "SelectedEE_DY_1500to2000";
    Procname[_EE_DY_2000to3000] = "SelectedEE_DY_2000to3000";
    Procname[_EndOf_EE_MCsignal_Normal] = "EndOf_SelectedEE_MCsignal_Normal";
    Procname[_EE_DYTauTau_10to50] = "SelectedEE_DYTauTau_10to50";
    Procname[_EE_DYTauTau_50toInf] = "SelectedEE_DYTauTau_50toInf";
    Procname[_EndOf_EE_DYTauTau_Normal] = "EndOf_SelectedEE_DYTauTau_Normal";
    Procname[_EE_ttbar] = "SelectedEE_ttbar";
    Procname[_EE_ttbar_700to1000] = "SelectedEE_ttbar_700to1000";
    Procname[_EE_ttbar_1000toInf] = "SelectedEE_ttbar_1000toInf";
    Procname[_EndOf_EE_ttbar_Normal] = "EndOf_SelectedEE_ttbar_Normal";
    Procname[_EE_tW] = "SelectedEE_tW";
    Procname[_EE_tbarW] = "SelectedEE_tbarW";
    Procname[_EE_ZZ] = "SelectedEE_ZZ";
    Procname[_EE_WZ] = "SelectedEE_WZ";
    Procname[_EE_WW] = "SelectedEE_WW";
    Procname[_EndOf_EE_VVnST_Normal] = "EndOf_SelectedEE_VVnST_Normal";
    Procname[_EE_WJets] = "SelectedEE_WJets";
    Procname[_EE_WJets_ext2v5] = "SelectedEE_WJets_ext2v5";
    Procname[_EndOf_EE_WJets_Normal] = "EndOf_SelectedEE_WJets";
    Procname[_EE_QCDEMEnriched_20to30] = "SelectedEE_QCDEMEnriched_20to30";
    Procname[_EE_QCDEMEnriched_30to50] = "SelectedEE_QCDEMEnriched_30to50";
    Procname[_EE_QCDEMEnriched_50to80] = "SelectedEE_QCDEMEnriched_50to80";
    Procname[_EE_QCDEMEnriched_80to120] = "SelectedEE_QCDEMEnriched_80to120";
    Procname[_EE_QCDEMEnriched_120to170] = "SelectedEE_QCDEMEnriched_120to170";
    Procname[_EE_QCDEMEnriched_170to300] = "SelectedEE_QCDEMEnriched_170to300";
    Procname[_EE_QCDEMEnriched_300toInf] = "SelectedEE_QCDEMEnriched_300toInf";
    Procname[_EndOf_EE_QCDEMEnriched_Normal] = "EndOf_SelectedEE_QCDEMEnriched_Normal";
    Procname[_EndOf_EE_MCbkg_Normal] = "EndOf_SelectedEE_MCbkg_Normal";
    Procname[_EE_DoubleEG_B] = "SelectedEE_DoubleEG_B";
    Procname[_EE_DoubleEG_C] = "SelectedEE_DoubleEG_C";
    Procname[_EE_DoubleEG_D] = "SelectedEE_DoubleEG_D";
    Procname[_EE_DoubleEG_E] = "SelectedEE_DoubleEG_E";
    Procname[_EE_DoubleEG_F] = "SelectedEE_DoubleEG_F";
    Procname[_EE_DoubleEG_G] = "SelectedEE_DoubleEG_G";
    Procname[_EE_DoubleEG_H] = "SelectedEE_DoubleEG_H";
    Procname[_EndOf_EE_DoubleEG_Normal] = "EndOf_SelectedEE_DoubleEG_Normal";
    Procname[_EE_SingleElectron_B] = "SelectedEE_SingleElectron_B";
    Procname[_EE_SingleElectron_C] = "SelectedEE_SingleElectron_C";
    Procname[_EE_SingleElectron_D] = "SelectedEE_SingleElectron_D";
    Procname[_EE_SingleElectron_E] = "SelectedEE_SingleElectron_E";
    Procname[_EE_SingleElectron_F] = "SelectedEE_SingleElectron_F";
    Procname[_EE_SingleElectron_G] = "SelectedEE_SingleElectron_G";
    Procname[_EE_SingleElectron_H] = "SelectedEE_SingleElectron_H";
    Procname[_EndOf_EE_SingleElectron_Normal] = "EndOf_SelectedEE_SingleElectron_Normal";
    Procname[_EndOf_EE_Data_Normal] = "EndOf_SelectedEE_Data_Normal";
    Procname[_EE_DY_Full] = "SelectedEE_DY_Full";
    Procname[_EE_DYTauTau_Full] = "SelectedEE_DYTauTau_Full";
    Procname[_EE_ttbar_Full] = "SelectedEE_ttbar_Full";
    Procname[_EE_VVnST] = "SelectedEE_VVnST";
    Procname[_EE_WJets_Full] = "SelectedEE_WJets_Full";
    Procname[_EE_QCDEMEnriched_Full] = "SelectedEE_QCDEMEnriched_Full";
    Procname[_EE_Bkg_Full] = "SelectedEE_Bkg_Full";
    Procname[_EE_DoubleEG_Full] = "SelectedEE_DoubleEG_Full";
    Procname[_EE_SingleElectron_Full] = "SelectedEE_SingleElectron_Full";
    Procname[_EndOf_EE_Special] = "EndOf_SelectedEE_Special";
    Procname[_Test_EE] = "Test_SelectedEE";
    Procname[_EndOf_EE] = "EndOf_SelectedEE";

    Procname[_EMu_DYTauTau_10to50] = "SelectedEMu_DYTauTau_10to50";
    Procname[_EMu_DYTauTau_50toInf] = "SelectedEMu_DYTauTau_50toInf";
    Procname[_EndOf_EMu_DYTauTau_Normal] = "EndOf_SelectedEMu_DYTauTau_Normal";
    Procname[_EMu_ttbar] = "SelectedEMu_ttbar";
    Procname[_EMu_ttbar_700to1000] = "SelectedEMu_ttbar_700to1000";
    Procname[_EMu_ttbar_1000toInf] = "SelectedEMu_ttbar_1000toInf";
    Procname[_EndOf_EMu_ttbar_Normal] = "EndOf_SelectedEMu_ttbar_Normal";
    Procname[_EMu_tW] = "SelectedEMu_tW";
    Procname[_EMu_tbarW] = "SelectedEMu_tbarW";
    Procname[_EMu_ZZ] = "SelectedEMu_ZZ";
    Procname[_EMu_WZ] = "SelectedEMu_WZ";
    Procname[_EMu_WW] = "SelectedEMu_WW";
    Procname[_EndOf_EMu_VVnST_Normal] = "EndOf_SelectedEMu_VVnST_Normal";
    Procname[_EMu_WJets] = "SelectedEMu_WJets";
    Procname[_EMu_WJets_ext2v5] = "SelectedEMu_WJets_ext2v5";
    Procname[_EndOf_EMu_WJets_Normal] = "EndOf_SelectedEMu_WJets";
    Procname[_EndOf_EMu_MCbkg_Normal] = "EndOf_SelectedEMu_MCbkg_Normal";
    Procname[_EMu_SingleMuon_B] = "SelectedEMu_SingleMuon_B";
    Procname[_EMu_SingleMuon_C] = "SelectedEMu_SingleMuon_C";
    Procname[_EMu_SingleMuon_D] = "SelectedEMu_SingleMuon_D";
    Procname[_EMu_SingleMuon_E] = "SelectedEMu_SingleMuon_E";
    Procname[_EMu_SingleMuon_F] = "SelectedEMu_SingleMuon_F";
    Procname[_EMu_SingleMuon_G] = "SelectedEMu_SingleMuon_G";
    Procname[_EMu_SingleMuon_H] = "SelectedEMu_SingleMuon_H";
    Procname[_EndOf_EMu_Data_Normal] = "EndOf_SelectedEMu_Data_Normal";
    Procname[_EMu_DYTauTau_Full] = "SelectedEMu_DYTauTau_Full";
    Procname[_EMu_ttbar_Full] = "SelectedEMu_ttbar_Full";
    Procname[_EMu_VVnST] = "SelectedEMu_VVnST";
    Procname[_EMu_WJets_Full] = "SelectedEMu_WJets_Full";
    Procname[_EMu_Bkg_Full] = "SelectedEMu_Bkg_Full";
    Procname[_EMu_SingleMuon_Full] = "SelectedEMu_SingleMuon_Full";
    Procname[_EndOf_EMu_Special] = "EndOf_SelectedEMu_Special";
    Procname[_Test_EMu] = "Test_SelectedEMu";
    Procname[_EndOf_EMu] = "EndOf_SelectedEMu";
    return;

}// end of PrepareProcNames()


void LocalFileMgr::CheckProcesses()
{
    Bool_t allOk = kTRUE;
    cout << "Checking processes: " << endl;
    for ( SelProc_t pr=_MuMu_DY_10to50; pr<_EndOf_EMu; pr=next(pr) )
    {
        this->SetProc(pr, kTRUE);
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
                else
                {
                    if ( pr < _EndOf_MuMu && !Tag[i].Contains("SelectedMuMu_") )
                    {
                        cout << "Process " << Procname[pr] << ": Tag[" << i << "] does not have 'SelectedMuMu_'." << endl;
                        allOk = kFALSE;
                    }
                    if ( pr > _EndOf_MuMu && pr < _EndOf_EE && !Tag[i].Contains("SelectedEE_") )
                    {
                        cout << "Process " << Procname[pr] << ": Tag[" << i << "] does not have 'SelectedEE_'." << endl;
                        allOk = kFALSE;
                    }
                    if ( pr > _EndOf_EE && pr < _EndOf_EMu && !Tag[i].Contains("SelectedEMu_") )
                    {
                        cout << "Process " << Procname[pr] << ": Tag[" << i << "] does not have 'SelectedEMu_'." << endl;
                        allOk = kFALSE;
                    }
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
                else
                {
                    if ( FileLocation[i].Contains("/*.root") )
                    {
                        cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] has '/*.root' instead of '.root'." << endl;
                        allOk = kFALSE;
                    }
                    else if ( !FileLocation[i].Contains(".root") )
                    {
                        cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not have '.root'." << endl;
                        allOk = kFALSE;
                    }
                    if ( ( pr > _None && pr < _EndOf_MuMu ) || ( pr > _EndOf_EE && pr < _EndOf_EMu ) )
                    {
                        if ( ROCCORR == kTRUE && !FileLocation[i].Contains("_roccor.root") )
                        {
                            cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not have '_roccor.root'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( pr < _EndOf_MuMu && !FileLocation[i].Contains("SelectedMuMu_") )
                    {
                        cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not have 'SelectedMuMu_'." << endl;
                        allOk = kFALSE;
                    }
                    if ( pr > _EndOf_MuMu && pr < _EndOf_EE && !FileLocation[i].Contains("SelectedEE_") )
                    {
                        cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not have 'SelectedEE_'." << endl;
                        allOk = kFALSE;
                    }
                    if ( pr > _EndOf_EE && pr < _EndOf_EMu && !FileLocation[i].Contains("SelectedEMu_") )
                    {
                        cout << "Process " << Procname[pr] << ": FileLocation[" << i << "] does not have 'SelectedEMu_'." << endl;
                        allOk = kFALSE;
                    }
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
                    if ( FullLocation[i].Contains("/*.root") )
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] has '/*.root' instead of '.root'." << endl;
                        allOk = kFALSE;
                    }
                    else if ( !FullLocation[i].Contains(".root") )
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '.root'." << endl;
                        allOk = kFALSE;
                    }
                    if ( ( pr > _None && pr < _EndOf_MuMu ) || ( pr > _EndOf_EE && pr < _EndOf_EMu ) )
                    {
                        if ( ROCCORR == kTRUE && !FullLocation[i].Contains("_roccor.root") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '_roccor.root'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _None && pr < _EndOf_MuMu_MCsignal_Normal) || pr == _MuMu_DY_Full )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedMuMu/MC_signal/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedMuMu/MC_signal/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FileLocation[i].Contains("SelectedMuMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedMuMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_MuMu_MCsignal_Normal && pr < _EndOf_MuMu_QCDMuEnriched_Normal) || ( pr > _MuMu_DY_Full && pr < _MuMu_SingleMuon_Full ) )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedMuMu/MC_bkg/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedMuMu/MC_bkg/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FileLocation[i].Contains("SelectedMuMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedMuMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_MuMu_MCbkg_Normal && pr < _EndOf_MuMu_Data_Normal) || pr == _MuMu_SingleMuon_Full )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedMuMu/Data/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedMuMu/Data/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedMuMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedMuMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( pr == _Test_MuMu )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/test/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/test/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedMuMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedMuMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_MuMu && pr < _EndOf_EE_MCsignal_Normal) || pr == _EE_DY_Full )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedEE/MC_signal/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedEE/MC_signal/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEE_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEE_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_EE_MCsignal_Normal && pr < _EndOf_EE_QCDEMEnriched_Normal) || ( pr > _EE_DY_Full && pr < _EE_DoubleEG_Full ) )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedEE/MC_bkg/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedEE/MC_bkg/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEE_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEE_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_EE_MCbkg_Normal && pr < _EndOf_EE_SingleElectron_Normal) || pr == _EE_DoubleEG_Full || pr == _EE_SingleElectron_Full )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedEE/Data/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedEE/Data/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEE_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEE_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( pr == _Test_EE )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/test/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/test/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEE_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEE_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_EE && pr < _EndOf_EMu_MCbkg_Normal) || ( pr > _EndOf_EMu_Data_Normal && pr < _EMu_SingleMuon_Full ) )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedEMu/MC_bkg/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedEMu/MC_bkg/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( ( pr > _EndOf_EMu_MCbkg_Normal && pr < _EndOf_EMu_Data_Normal) || pr == _EMu_SingleMuon_Full )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/SelectedEMu/Data/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/SelectedEMu/Data/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( pr == _Test_EMu )
                    {
                        if ( !FullLocation[i].Contains("/media/sf_DATA/test/") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have '/media/sf_DATA/test/'." << endl;
                            allOk = kFALSE;
                        }
                        if ( !FullLocation[i].Contains("SelectedEMu_") )
                        {
                            cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not have 'SelectedEMu_'." << endl;
                            allOk = kFALSE;
                        }
                    }
                    if ( FullLocation[i][0] != '/')
                    {
                        cout << "Process " << Procname[pr] << ": FullLocation[" << i << "] does not begin with '/'." << endl;
                        allOk = kFALSE;
                    }
                }
            }// End of for(FullLocation[i])

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
        if ( !HistLocation.Length() )
        {
            cout << "Process " << Procname[pr] << ": no HistLocation found." << endl;
            allOk = kFALSE;
        }
        else
        {
            if ( HistLocation[HistLocation.Length()-1] != '/') {
                cout << "Process " << Procname[pr] << ": HistLocation does not end with '/'" << endl;
                allOk = kFALSE;
            }
            if ( HistLocation[0] != '/') {
                cout << "Process " << Procname[pr] << ": HistLocation does not begin with '/'" << endl;
                allOk = kFALSE;
            }
            if ( pr < _EE_DY_10to50 && HistLocation != "/media/sf_DATA/SelectedMuMu/Histos/" )
            {
                cout << "Process " << Procname[pr] << ": HistLocation does is not '/media/sf_DATA/SelectedMuMu/Histos/'" << endl;
                allOk = kFALSE;
            }
            if ( pr > _EndOf_MuMu && pr < _EMu_DYTauTau_10to50 && HistLocation != "/media/sf_DATA/SelectedEE/Histos/" )
            {
                cout << "Process " << Procname[pr] << ": HistLocation does is not '/media/sf_DATA/SelectedEE/Histos/'" << endl;
                allOk = kFALSE;
            }
            if ( pr > _EndOf_EE && HistLocation != "/media/sf_DATA/SelectedEMu/Histos/" )
            {
                cout << "Process " << Procname[pr] << ": HistLocation does is not '/media/sf_DATA/SelectedEMu/Histos/'" << endl;
                allOk = kFALSE;
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
        if ( !nEvents.size() )
        {
            cout << "Process " << Procname[pr] << ": no nEvents found." << endl;
            allOk = kFALSE;
        }
        else
        {
            for (Int_t i=0; i<(int(nEvents.size())); i++)
            {
                Bool_t ok = kTRUE; // Not all processes are supposed to have events. This will help to prevent unnecessary messages.
                if ( !nEvents[i] && pr != _MuMu_QCDMuEnriched_15to20 && pr != _MuMu_QCDMuEnriched_20to30 && pr != _EE_QCDEMEnriched_20to30 )
                {
                    if ( ( pr != _MuMu_QCDMuEnriched_Full || i > 1 ) && ( pr != _MuMu_Bkg_Full || i < 15 || i > 16 ) &&
                         ( pr != _EE_QCDEMEnriched_30to50 || i > 0 ) && ( pr != _EE_QCDEMEnriched_Full || i > 1 ) &&
                         ( pr != _EE_Bkg_Full || i < 15 || i > 16 ) )
                        ok = kFALSE;

                    if ( ok == kFALSE )
                    {
                        cout << "Process " << Procname[pr] << ": no nEvents[" << i << "] found." << endl;
                        allOk = kFALSE;
                    }
                }
            }
        }
        if ( pr < _EndOf_MuMu_QCDMuEnriched_Normal || ( pr > _EndOf_MuMu_Data_Normal && pr < _MuMu_SingleMuon_Full ) |
             ( pr > _EndOf_MuMu && pr < _EndOf_EE_QCDEMEnriched_Normal ) || ( pr > _EndOf_EE_Data_Normal && pr < _EE_DoubleEG_Full ) ||
             ( pr > _EndOf_EE && pr < _EndOf_EMu_MCbkg_Normal ) || ( pr > _EndOf_EMu_Data_Normal && pr < _EMu_SingleMuon_Full ) ||
             pr == _Test_MuMu || pr == _Test_EE || pr == _Test_EMu )
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
                 TreeName.size() != nEvents.size() || nEvents.size() != Wsum.size() )
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
            if ( pr < _EndOf_MuMu_MCsignal_Normal || pr == _MuMu_DY_Full || ( pr > _EndOf_MuMu && pr < _EndOf_EE_MCsignal_Normal ) || pr == _EE_DY_Full )
            {
                if ( Type != "SIGNAL" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT SIGNAL." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( ( pr > _EndOf_MuMu_MCsignal_Normal && pr < _EndOf_MuMu_QCDMuEnriched_Normal ) || ( pr > _MuMu_DY_Full && pr < _MuMu_SingleMuon_Full ) ||
                      ( pr > _EndOf_EE_MCsignal_Normal && pr < _EndOf_EE_QCDEMEnriched_Normal ) || ( pr > _EE_DY_Full && pr < _EE_DoubleEG_Full ) ||
                      ( pr > _EndOf_EE && pr < _EndOf_EMu_MCbkg_Normal ) || ( pr > _EndOf_EMu_Data_Normal && pr < _EMu_SingleMuon_Full ) )
            {
                if ( Type != "BKG" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT BKG." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( ( pr > _EndOf_MuMu_MCbkg_Normal && pr < _EndOf_MuMu_Data_Normal ) || pr == _MuMu_SingleMuon_Full ||
                      ( pr > _EndOf_EE_MCbkg_Normal && pr < _EndOf_EE_SingleElectron_Normal ) || pr == _EE_DoubleEG_Full || pr == _EE_SingleElectron_Full ||
                      ( pr > _EndOf_EMu_MCbkg_Normal && pr < _EndOf_EMu_Data_Normal ) || pr == _EMu_SingleMuon_Full )
            {
                if ( Type != "DATA" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT DATA." << endl;
                    allOk = kFALSE;
                }
            }
            else if ( pr == _Test_MuMu || pr == _Test_EE || pr == _Test_EMu )
            {
                if ( Type != "TEST" )
                {
                    cout << "Process " << Procname[pr] << ": is said to be NOT TEST." << endl;
                    allOk = kFALSE;
                }
            }
        }       
        vector<SelProc_t> forChecking = FindProc(Procname[pr], kFALSE, kFALSE);
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
        if ( HistLocation.Length() )
        {
            cout << "Process " << Procname[pr] << ": ClearProc() did not clear the HistLocation." << endl;
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
