// Header for file management
// First working version created by Marijus Ambrozas 2018.08.23
// 2018.08.14: Added HistLocation for every process
// 2018.08.31: Added an ability to change between files that contain events selected before and after the Rochester correction (you need to change ROCCORR variable to switch)

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
    _MuMu_WJets, _EndOf_MuMu_WJets,
    _MuMu_QCDMuEnriched_15to20, _MuMu_QCDMuEnriched_20to30, _MuMu_QCDMuEnriched_30to50, _MuMu_QCDMuEnriched_50to80, _MuMu_QCDMuEnriched_80to120,
    _MuMu_QCDMuEnriched_120to170, _MuMu_QCDMuEnriched_170to300, _MuMu_QCDMuEnriched_300to470, _MuMu_QCDMuEnriched_470to600, _MuMu_QCDMuEnriched_600to800,
    _MuMu_QCDMuEnriched_800to1000, _MuMu_QCDMuEnriched_1000toInf,
    _EndOf_MuMu_QCDMuEnriched_Normal,
    _EndOf_MuMu_MCbkg_Normal,
    _MuMu_SingleMuon_B, _MuMu_SingleMuon_C, _MuMu_SingleMuon_D, _MuMu_SingleMuon_E, _MuMu_SingleMuon_F, _MuMu_SingleMuon_G, _MuMu_SingleMuon_H,
    _EndOf_MuMu_Data_Normal,
    _MuMu_DY_Full,
    _MuMu_DYTauTau_Full, _MuMu_ttbar_Full, _MuMu_VVnST, _MuMu_QCDMuEnriched_Full, _MuMu_Bkg_Full,
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
    _EE_WJets, _EndOf_EE_WJets,
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
    _EE_DYTauTau_Full, _EE_ttbar_Full, _EE_VVnST, _EE_QCDEMEnriched_Full, _EE_Bkg_Full,
    _EE_DoubleEG_Full,  _EE_SingleElectron_Full,
    _EndOf_EE_Special,
    _Test_EE,
    _EndOf_EE,

    _EMu_DYTauTau_10to50, _EMu_DYTauTau_50toInf, _EndOf_EMu_DYTauTau_Normal,
    _EMu_ttbar, _EMu_ttbar_700to1000, _EMu_ttbar_1000toInf, _EndOf_EMu_ttbar_Normal,
    _EMu_tW, _EMu_tbarW, _EMu_ZZ, _EMu_WZ, _EMu_WW, _EndOf_EMu_VVnST_Normal,
    _EMu_WJets, _EndOf_EMu_WJets,
    _EndOf_EMu_MCbkg_Normal,
    _EMu_SingleMuon_B, _EMu_SingleMuon_C, _EMu_SingleMuon_D, _EMu_SingleMuon_E, _EMu_SingleMuon_F, _EMu_SingleMuon_G, _EMu_SingleMuon_H,
    _EndOf_EMu_Data_Normal,
    _EMu_DYTauTau_Full, _EMu_ttbar_Full, _EMu_VVnST, _EMu_Bkg_Full,
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
  else if ( pr == _MuMu_QCDMuEnriched_1000toInf || pr == _EE_QCDEMEnriched_300toInf || pr == _EE_SingleElectron_H || pr == _EMu_WJets )
      return SelProc_t(int(pr)+3);
  else if ( pr == _MuMu_DY_2000to3000 || pr == _MuMu_DYTauTau_50toInf || pr == _MuMu_ttbar_1000toInf || pr == _MuMu_WW || pr == _MuMu_WJets ||
            pr == _EndOf_MuMu_QCDMuEnriched_Normal || pr == _MuMu_SingleMuon_H || pr == _MuMu_SingleMuon_Full || pr == _Test_MuMu ||
            pr == _EE_DY_2000to3000 || pr == _EE_DYTauTau_50toInf || pr == _EE_ttbar_1000toInf || pr == _EE_WW || pr == _EE_WJets ||
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
        void GetProc ( SelProc_t pr = _None, Bool_t ClearOld = kTRUE );
        void ClearProc ();

        void SwitchROCCORR();

private:
        Bool_t namesSet = kFALSE;
        Bool_t processesChecked = kFALSE;
        Bool_t ROCCORR = kTRUE; // Change this to kFALSE if you don't want to work with events that were selected after the Rochester correction was applied

        void PrepareProcNames ();
        void CheckProcesses ();

};// end of class definition


// ---------- Constructor ---------- //

LocalFileMgr::LocalFileMgr ( SelProc_t pr )
{
    if ( namesSet == kFALSE ) { this->PrepareProcNames(); namesSet = kTRUE; }
    if ( processesChecked == kFALSE ) { this->CheckProcesses(); processesChecked = kTRUE; }
    CurrentProc = pr;
    this->GetProc(CurrentProc, kTRUE);
}


// ----------- Functions ----------- //

void LocalFileMgr::NextProc()
{
    CurrentProc = next(CurrentProc);
    this->GetProc(CurrentProc, kTRUE);
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
        this->GetProc(CurrentProc, kTRUE);
    }
}


void LocalFileMgr::GetProc ( SelProc_t pr, Bool_t ClearOld )
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

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 20822 );
        else nEvents.push_back( 20814 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 44224 );
        else nEvents.push_back( 44186 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 27278 );
        else nEvents.push_back( 27207 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );
    }
    else if( pr ==_MuMu_DY_50to100 ) // Only MuMu evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";
        Tag.push_back( "SelectedMuMu_DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 13702787 );
        else nEvents.push_back( 13699599 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M50to100"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 145625 );
        else nEvents.push_back( 145608 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2060461 );
        else nEvents.push_back( 2060233 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 51625 );
        else nEvents.push_back( 51624 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M200to400"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 58616 );
        else nEvents.push_back( 58615 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M400to500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 62122 );
        else nEvents.push_back( 62123 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M500to700"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 64895 );
        else nEvents.push_back( 64893 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M700to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65828 );
        else nEvents.push_back( 65829 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 66149 );
        else nEvents.push_back( 66149 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1000to1500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37176.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 68769 );
        else nEvents.push_back( 68770 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1500to2000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr ==_MuMu_DY_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 69286 );
        else nEvents.push_back( 69286 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M2000to3000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr ==_MuMu_DY_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 20822 );
        else nEvents.push_back( 20814 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 44224 );
        else nEvents.push_back( 44186 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7446893+16016651+9815322 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 27278 );
        else nEvents.push_back( 27207 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation + Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26175605.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 13702787 );
        else nEvents.push_back( 13699599 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M50to100"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 145625 );
        else nEvents.push_back( 145608 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 233822+3199473 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2060461 );
        else nEvents.push_back( 2060233 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M100to200_ext"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56340.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 51625 );
        else nEvents.push_back( 51624 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M200to400"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50136.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 58616 );
        else nEvents.push_back( 58615 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M400to500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48188.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 62122 );
        else nEvents.push_back( 62123 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M500to700"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 44984.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 64895 );
        else nEvents.push_back( 64893 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M700to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 43496.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 65828 );
        else nEvents.push_back( 65829 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 40110.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 66149 );
        else nEvents.push_back( 66149 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1000to1500"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37176.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 68769 );
        else nEvents.push_back( 68770 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M1500to2000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYMuMu_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 33360.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 69286 );
        else nEvents.push_back( 69286 );
        Location = "SelectedMuMu/MC_signal/SelectedMuMu_DYMuMu_M2000to3000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 25 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 );  //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 31777 );
        else nEvents.push_back( 31733 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 25 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 );  //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 31777 );
        else nEvents.push_back( 31733 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );  //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 473907 );
        else nEvents.push_back( 473806 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 479397 );
        else nEvents.push_back( 479287 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 269661 );
        else nEvents.push_back( 269627 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 152053 );
        else nEvents.push_back( 152026 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 473907 );
        else nEvents.push_back( 473806 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 479397 );
        else nEvents.push_back( 479287 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 269661 );
        else nEvents.push_back( 269627 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 152053 );
        else nEvents.push_back( 152026 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52500 );
        else nEvents.push_back( 52490 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52489 );
        else nEvents.push_back( 52492 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 30567 );
        else nEvents.push_back( 30557 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47584 );
        else nEvents.push_back( 47568 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35594 );
        else nEvents.push_back( 35585 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52500 );
        else nEvents.push_back( 52490 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52489 );
        else nEvents.push_back( 52492 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 30567 );
        else nEvents.push_back( 30557 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47584 );
        else nEvents.push_back( 47568 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35594 );
        else nEvents.push_back( 35585 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 ); // I get Wsum=137540054
        if ( ROCCORR == kTRUE ) nEvents.push_back( 92 );
        else nEvents.push_back( 93 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 823 );
        else nEvents.push_back( 823 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext"+RocCorr+".root";        // There also is madgraph version
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_80to120 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 9 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_120to170 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5 );
        else nEvents.push_back( 5 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_300to470 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_800to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 6 );
        else nEvents.push_back( 6 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _MuMu_QCDMuEnriched_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8 );
        else nEvents.push_back( 8 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 2 );
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 9 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5 );
        else nEvents.push_back( 5 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 6 );
        else nEvents.push_back( 6 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8 );
        else nEvents.push_back( 8 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_Bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 11 );
        else nEvents.push_back( 11 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 26 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 26 );
        else nEvents.push_back( 25 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 );  //  NNLO Xsec
        if ( ROCCORR == kTRUE ) nEvents.push_back( 31777 );
        else nEvents.push_back( 31733 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE ) nEvents.push_back( 473907 );
        else nEvents.push_back( 473806 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 479397 );
        else nEvents.push_back( 479287 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 269661 );
        else nEvents.push_back( 269627 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 152053 );
        else nEvents.push_back( 152026 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52500 );
        else nEvents.push_back( 52490 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 52489 );
        else nEvents.push_back( 52492 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 30567 );
        else nEvents.push_back( 30557 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_ZZ"+RocCorr+".root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 47584 );
        else nEvents.push_back( 47568 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WZ"+RocCorr+".root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 35594 );
        else nEvents.push_back( 35585 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WW"+RocCorr+".root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); // Wsum.push_back( 86731698.0 ); // I get Wsum=137540054
        if ( ROCCORR == kTRUE ) nEvents.push_back( 92 );
        else nEvents.push_back( 93 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); // Wsum.push_back( 86731698.0 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 823 );
        else nEvents.push_back( 823 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_WJetsToLNu_ext"+RocCorr+".root";        // There also is madgraph version
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
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt50to80"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120" ); Xsec.push_back( 2758420*0.03844 ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 9 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt80to120_ext1" ); Xsec.push_back( 2758420*0.03844  ); Wsum.push_back( 13555323+9797243 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt80to120_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170" ); Xsec.push_back( 469797*0.05362 ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5 );
        else nEvents.push_back( 5 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt120to170_backup" ); Xsec.push_back( 469797*0.05362  ); Wsum.push_back( 8042720+11938137 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt120to170_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_ext1" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt170to300_backup" ); Xsec.push_back( 117989*0.07335 ); Wsum.push_back( 7947158+9403070+19607775 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt170to300_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext1" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt300to470_ext2" ); Xsec.push_back( 7820.25*0.10196 ); Wsum.push_back( 7937587+16452587+24605502 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 10 );
        else nEvents.push_back( 10 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt300to470_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt470to600_ext1" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 3851523+5663755 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1 );
        else nEvents.push_back( 1 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt470to600_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

// DID NOT FIND THIS ONE
//        Tag.push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec.push_back( 645.528*0.12242 ); Wsum.push_back( 1.0 ); nEvents.push_back( 0 );
//        Location = "QCDMuEnriched_Pt470to600_ext2";
//        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_ext1" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 4 );
        else nEvents.push_back( 4 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt600to800_backup" ); Xsec.push_back( 187.109*0.13412 ); Wsum.push_back( 4010135+5971173+9756852 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt600to800_backup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 6 );
        else nEvents.push_back( 6 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2" ); Xsec.push_back( 32.3486*0.14552 ); Wsum.push_back( 3962747+5838539+9966146 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 3 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt800to1000_ext2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 8 );
        else nEvents.push_back( 8 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1" ); Xsec.push_back( 10.4305*0.15544 ); Wsum.push_back( 3861436+9609820 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 3 );
        else nEvents.push_back( 2 );
        Location = "SelectedMuMu/MC_bkg/SelectedMuMu_QCDMuEnriched_Pt1000toInf_ext1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_B ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_B" ); Wsum.push_back( 108561074 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2862839 );
        else nEvents.push_back( 2862076 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_B"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_C ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1735877 );
        else nEvents.push_back( 1735466 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_D ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2880898 );
        else nEvents.push_back( 2880299 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_E ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2652397 );
        else nEvents.push_back( 2651780 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_F ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2044728 );
        else nEvents.push_back( 2044183 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_G ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_G" ); Wsum.push_back( 138710659 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5055758 );
        else nEvents.push_back( 5054446 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_G"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_H ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver2" ); Wsum.push_back( 141936183 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5170278 );
        else nEvents.push_back( 5168884 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 157691 );
        else nEvents.push_back( 157648 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _MuMu_SingleMuon_Full ) // Number of events in the original sample (before selection) are written into Wsum
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedMuMu/Histos/";

        Tag.push_back( "SelectedMuMu_SingleMuon_B" ); Wsum.push_back( 108561074 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2862839 );
        else nEvents.push_back( 2862076 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_B"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 1735877 );
        else nEvents.push_back( 1735466 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2880898 );
        else nEvents.push_back( 2880299 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2652397 );
        else nEvents.push_back( 2651780 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2044728 );
        else nEvents.push_back( 2044183 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_G" ); Wsum.push_back( 138710659 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5055758 );
        else nEvents.push_back( 5054446 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_G"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver2" ); Wsum.push_back( 141936183 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 5170278 );
        else nEvents.push_back( 5168884 );
        Location = "SelectedMuMu/Data/SelectedMuMu_SingleMuon_Hver2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedMuMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 157691 );
        else nEvents.push_back( 157648 );
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

        Tag.push_back( "SelectedEE_DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 10426 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 22623 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 13836 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_50to100 ) // Only EE evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 7706926 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M50to100.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_100to200 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 86786 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 1224185 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200_ext.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_200to400 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 33801 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M200to400.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_400to500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 40021 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M400to500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_500to700 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 42770 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M500to700.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_700to800 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 45351 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M700to800.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_800to1000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 46888 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M800to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_DY_1000to1500 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 47223 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1000to1500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_1500to2000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 50387 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1500to2000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_2000to3000 )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 50986 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M2000to3000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DY_Full )
    {
        isMC = kTRUE;
        Type = "SIGNAL";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYEE_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 10426 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 22623 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7447023+16016761+9811434 ); nEvents.push_back( 13836 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M50to100" ); Xsec.push_back( 5869.58346/3.0 ); Wsum.push_back( 26166194.0 ); nEvents.push_back( 7706926 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M50to100.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M100to200" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 86786 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M100to200_ext" ); Xsec.push_back( 226/3.0 ); Wsum.push_back( 234322+3203563 ); nEvents.push_back( 1224185 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M100to200_ext.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M200to400" ); Xsec.push_back( 7.67/3.0 ); Wsum.push_back( 56144.0 ); nEvents.push_back( 33801 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M200to400.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M400to500" ); Xsec.push_back( 0.423/3.0 ); Wsum.push_back( 50420.0 ); nEvents.push_back( 40021 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M400to500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M500to700" ); Xsec.push_back( 0.24/3.0 ); Wsum.push_back( 48039.0 ); nEvents.push_back( 42770 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M500to700.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M700to800" ); Xsec.push_back( 0.035/3.0 ); Wsum.push_back( 46114.0 ); nEvents.push_back( 45351 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M700to800.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M800to1000" ); Xsec.push_back( 0.03/3.0 ); Wsum.push_back( 44256.0 ); nEvents.push_back( 46888 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M800to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M1000to1500" ); Xsec.push_back( 0.016/3.0 ); Wsum.push_back( 39712.0 ); nEvents.push_back( 47223 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1000to1500.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M1500to2000" ); Xsec.push_back( 0.002/3.0 ); Wsum.push_back( 37287.0 ); nEvents.push_back( 50387 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M1500to2000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYEE_M2000to3000" ); Xsec.push_back( 0.00054/3.0 ); Wsum.push_back( 34031.0 ); nEvents.push_back( 50986 );
        Location = "SelectedEE/MC_signal/SelectedEE_DYEE_M2000to3000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }   
    else if ( pr == _EE_DYTauTau_10to50 ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 7 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 26 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 17 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20257 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 7 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 26 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 17 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20257 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277386 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279721 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162178 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89582 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ttbar_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277386 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279721 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162178 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89582 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_tW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32181 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32096 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17813 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28695 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21068 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EE_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32181 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32096 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17813 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28695 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21068 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EE_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEE/Histos/";

        Tag.push_back( "SelectedEE_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 200 ); // I get Wsum=137540054
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 1338 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext.root";        // There also is madgraph version
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

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 7 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 26 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_v2.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 ); nEvents.push_back( 17 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M10to50_ext1v1.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); nEvents.push_back( 20257 ); //  NNLO Xsec
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEE/MC_bkg/SelectedEE_DYTauTau_M50toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 277386 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 ); nEvents.push_back( 279721 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbarBackup.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M700to1000" ); Xsec.push_back( 76.605 ); nEvents.push_back( 162178 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M700to1000.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ttbar_M1000toInf" ); Xsec.push_back( 20.578 ); nEvents.push_back( 89582 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEE/MC_bkg/SelectedEE_ttbar_M1000toInf.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 ); nEvents.push_back( 32181 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 ); nEvents.push_back( 32096 );
        Location = "SelectedEE/MC_bkg/SelectedEE_tbarW.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 ); nEvents.push_back( 17813 );
        Location = "SelectedEE/MC_bkg/SelectedEE_ZZ.root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 ); nEvents.push_back( 28695 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WZ.root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 ); nEvents.push_back( 21068 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WW.root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 200 ); // I get Wsum=137540054
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu.root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEE_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 86731698.0 ); nEvents.push_back( 1338 );
        Location = "SelectedEE/MC_bkg/SelectedEE_WJetsToLNu_ext.root";        // There also is madgraph version
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

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 17 );
        else nEvents.push_back( 17 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 43 );
        else nEvents.push_back( 43 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 32 );
        else nEvents.push_back( 32 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_DYTauTau_50toInf ) // Only TauTau evens are counted in Wsum and nEvents
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); //  NNLO Xsec
        if ( ROCCORR == kTRUE) nEvents.push_back( 39143 );
        else nEvents.push_back( 39126 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_DYTauTau_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 17 );
        else nEvents.push_back( 17 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 43 );
        else nEvents.push_back( 43 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 7407794+15912921+9759664 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 32 );
        else nEvents.push_back( 32 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); //  NNLO Xsec
        if ( ROCCORR == kTRUE) nEvents.push_back( 39143 );
        else nEvents.push_back( 39126 );
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
        if ( ROCCORR == kTRUE) nEvents.push_back( 630864 );
        else nEvents.push_back( 630854 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 638032 );
        else nEvents.push_back( 638011 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ttbar_700to1000 )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 373448 );
        else nEvents.push_back( 373417 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_SelectedMuMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_ttbar_1000toInf )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 208898 );
        else nEvents.push_back( 208900 );
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
        if ( ROCCORR == kTRUE) nEvents.push_back( 630864 );
        else nEvents.push_back( 630854 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 67632273+68317507 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 638032 );
        else nEvents.push_back( 638011 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 373448 );
        else nEvents.push_back( 373417 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 208898 );
        else nEvents.push_back( 208900 );
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
        if ( ROCCORR == kTRUE) nEvents.push_back( 72663 );
        else nEvents.push_back( 72654 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_tbarW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 72852 );
        else nEvents.push_back( 72851 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_ZZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 1911 );
        else nEvents.push_back( 1911 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_WZ )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 9595 );
        else nEvents.push_back( 9598 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_WW )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 47772 );
        else nEvents.push_back( 47774 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_VVnST )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 72663 );
        else nEvents.push_back( 72654 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 72852 );
        else nEvents.push_back( 72851 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 1911 );
        else nEvents.push_back( 1911 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 9595 );
        else nEvents.push_back( 9598 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 47772 );
        else nEvents.push_back( 47774 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if( pr == _EMu_WJets )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 ); // I get Wsum=137540054
        if ( ROCCORR == kTRUE) nEvents.push_back( 469 );
        else nEvents.push_back( 469 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 3191 );
        else nEvents.push_back( 3186 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext"+RocCorr+".root";        // There also is madgraph version
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_Bkg_Full )
    {
        isMC = kTRUE;
        Type = "BKG";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 17 );
        else nEvents.push_back( 17 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_v2" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 43 );
        else nEvents.push_back( 43 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_v2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M10to50_ext1v1" ); Xsec.push_back( 18610.0/3.0 ); Wsum.push_back( 33080379.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 32 );
        else nEvents.push_back( 32 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M10to50_ext1v1"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_DYTauTau_M50toInf" ); Xsec.push_back( 1921.8 ); Wsum.push_back( 27277866.0 ); //  NNLO Xsec
        if ( ROCCORR == kTRUE) nEvents.push_back( 39143 );
        else nEvents.push_back( 39126 );
//        Xsec->push_back( 6104.0/3.0 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_DYTauTau_M50toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
        if ( ROCCORR == kTRUE) nEvents.push_back( 630864 );
        else nEvents.push_back( 630854 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbarBackup" ); Xsec.push_back( 734.577 ); Wsum.push_back( 135949780.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 638032 );
        else nEvents.push_back( 638011 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbarBackup"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M700to1000" ); Xsec.push_back( 76.605 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 373448 );
        else nEvents.push_back( 373417 );
        Wsum.push_back( 38422582.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M700to1000"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ttbar_M1000toInf" ); Xsec.push_back( 20.578 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 208898 );
        else nEvents.push_back( 208900 );
        Wsum.push_back( 24561630.0 );                                       //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ttbar_M1000toInf"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6952830.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 72663 );
        else nEvents.push_back( 72654 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_tbarW" ); Xsec.push_back( 35.85 ); Wsum.push_back( 6933093.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 72852 );
        else nEvents.push_back( 72851 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_tbarW"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_ZZ" ); Xsec.push_back( 16.523 ); Wsum.push_back( 998034.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 1911 );
        else nEvents.push_back( 1911 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_ZZ"+RocCorr+".root";                  // NOT SURE (there also is ZZTo4L), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WZ" ); Xsec.push_back( 47.13 ); Wsum.push_back( 2995828.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 9595 );
        else nEvents.push_back( 9598 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WZ"+RocCorr+".root";                  // NOT SURE (there also is WZTo3LNu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WW" ); Xsec.push_back( 118.7 ); Wsum.push_back( 6987123.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 47772 );
        else nEvents.push_back( 47774 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WW"+RocCorr+".root";                  // NOT SURE (there also is WWTo2L2Nu), but probably ok
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 ); // I get Wsum=137540054
        if ( ROCCORR == kTRUE) nEvents.push_back( 469 );
        else nEvents.push_back( 469 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_WJetsToLNu_ext" ); Xsec.push_back( 61526.7 ); Wsum.push_back( 137540054.0 ); //Wsum.push_back( 86731698.0 );
        if ( ROCCORR == kTRUE) nEvents.push_back( 3191 );
        else nEvents.push_back( 3186 );
        Location = "SelectedEMu/MC_bkg/SelectedEMu_WJetsToLNu_ext"+RocCorr+".root";        // There also is madgraph version
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

        Tag.push_back( "SelectedEMu_SingleMuon_B" ); Wsum.push_back( 108561074 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 49756 );
        else nEvents.push_back( 49749 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_C )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 29763 );
        else nEvents.push_back( 29758 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_D )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 49397 );
        else nEvents.push_back( 49393 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_E )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 44041 );
        else nEvents.push_back( 44030 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_F )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 34504 );
        else nEvents.push_back( 34498 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_G )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_G" ); Wsum.push_back( 138710659 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 84902 );
        else nEvents.push_back( 84896 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_H )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_Hver2" ); Wsum.push_back( 141936183 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 85368 );
        else nEvents.push_back( 85351 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2603 );
        else nEvents.push_back( 2604 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver3"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );
    }
    else if ( pr == _EMu_SingleMuon_Full )
    {
        isMC = kFALSE;
        Type = "DATA";
        HistLocation = "/media/sf_DATA/SelectedEMu/Histos/";

        Tag.push_back( "SelectedEMu_SingleMuon_B" ); Wsum.push_back( 108561074 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 49756 );
        else nEvents.push_back( 49749 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_B"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_C" ); Wsum.push_back( 64715287 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 29763 );
        else nEvents.push_back( 29758 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_C"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_D" ); Wsum.push_back( 96652779 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 49397 );
        else nEvents.push_back( 49393 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_D"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_E" ); Wsum.push_back( 87358348 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 44041 );
        else nEvents.push_back( 44030 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_E"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_F" ); Wsum.push_back( 64986568 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 34504 );
        else nEvents.push_back( 34498 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_F"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_G" ); Wsum.push_back( 138710659 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 84902 );
        else nEvents.push_back( 84896 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_G"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver2" ); Wsum.push_back( 141936183 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 85368 );
        else nEvents.push_back( 85351 );
        Location = "SelectedEMu/Data/SelectedEMu_SingleMuon_Hver2"+RocCorr+".root";
        TreeName.push_back( "DYTree" ); FileLocation.push_back( Location ); FullLocation.push_back( BaseLocation+Location );

        Tag.push_back( "SelectedEMu_SingleMuon_Hver3" ); Wsum.push_back( 4386928 );
        if ( ROCCORR == kTRUE ) nEvents.push_back( 2603 );
        else nEvents.push_back( 2604 );
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
    if ( notify == kTRUE ) cout << "Searched for: " << search << "\nFound: ";
    vector<SelProc_t> Result;
    SelProc_t first = _None, last = _None;
    if ( srch.Contains("DY") || srch.Contains("dy") || srch.Contains("dY") || srch.Contains("Dy") || srch.Contains("DrellYan") || srch.Contains("DRELLYAN") )
    {
        if ( srch.Contains("TauTau") || srch.Contains("tautau") || srch.Contains("TAUTAU") || srch.Contains("Ditau") || srch.Contains("DiTau") ||
             srch.Contains("ditau") || srch.Contains("diTau") || srch.Contains("DITAU") )
        {
            if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
                 srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
            {
                if ( srch.Contains("full") || srch.Contains("Full") || srch.Contains("FULL") )
                {
                    Result.push_back(_MuMu_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_DYTauTau_Full] << "." << endl;
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
                        first = _EndOf_MuMu_DYTauTau_Normal;
                    else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                              srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                              srch.Contains("3000_TO") )
                        first = _EndOf_MuMu_DYTauTau_Normal;
                    else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                              srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                              srch.Contains("2000_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                              srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                              srch.Contains("1500_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                              srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                              srch.Contains("1000_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                              srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                              srch.Contains("800_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                              srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                              srch.Contains("700_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                              srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                              srch.Contains("500_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                              srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                              srch.Contains("400_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                              srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                              srch.Contains("200_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                              srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                              srch.Contains("100_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                              srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                        first = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                              srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                        first = _MuMu_DYTauTau_10to50;

                    if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                         srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                         srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                         srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                         srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                         srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                         srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                         srch.Contains("- Inf") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                              srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                              srch.Contains("TO_3000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                              srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                              srch.Contains("TO_2000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                              srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                              srch.Contains("TO_1500") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                              srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                              srch.Contains("TO_1000") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                              srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                              srch.Contains("TO_800") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                              srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                              srch.Contains("TO_700") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                              srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                              srch.Contains("TO_500") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                              srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                              srch.Contains("TO_400") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                              srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                              srch.Contains("TO_200") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                              srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                              srch.Contains("TO_100") )
                        last = _MuMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                              srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                              srch.Contains("TO_50") )
                        last = _MuMu_DYTauTau_10to50;
                    else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                              srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                              srch.Contains("TO_10") )
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

            else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                      srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
            {
                if ( srch.Contains("full") || srch.Contains("Full") || srch.Contains("FULL") )
                {
                    Result.push_back(_EE_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_EE_DYTauTau_Full] << "." << endl;
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
                        first = _EndOf_EE_DYTauTau_Normal;
                    else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                              srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                              srch.Contains("3000_TO") )
                        first = _EndOf_EE_DYTauTau_Normal;
                    else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                              srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                              srch.Contains("2000_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                              srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                              srch.Contains("1500_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                              srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                              srch.Contains("1000_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                              srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                              srch.Contains("800_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                              srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                              srch.Contains("700_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                              srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                              srch.Contains("500_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                              srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                              srch.Contains("400_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                              srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                              srch.Contains("200_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                              srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                              srch.Contains("100_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                              srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                        first = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                              srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                        first = _EE_DYTauTau_10to50;

                    if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                         srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                         srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                         srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                         srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                         srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                         srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                         srch.Contains("- Inf") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                              srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                              srch.Contains("TO_3000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                              srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                              srch.Contains("TO_2000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                              srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                              srch.Contains("TO_1500") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                              srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                              srch.Contains("TO_1000") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                              srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                              srch.Contains("TO_800") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                              srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                              srch.Contains("TO_700") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                              srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                              srch.Contains("TO_500") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                              srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                              srch.Contains("TO_400") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                              srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                              srch.Contains("TO_200") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                              srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                              srch.Contains("TO_100") )
                        last = _EE_DYTauTau_50toInf;
                    else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                              srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                              srch.Contains("TO_50") )
                        last = _EE_DYTauTau_10to50;
                    else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                              srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                              srch.Contains("TO_10") )
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

            else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                      srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
            {
                if ( srch.Contains("full") || srch.Contains("Full") || srch.Contains("FULL") )
                {
                    Result.push_back(_EMu_DYTauTau_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_DYTauTau_Full] << "." << endl;
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
                        first = _EndOf_EMu_DYTauTau_Normal;
                    else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                              srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                              srch.Contains("3000_TO") )
                        first = _EndOf_EMu_DYTauTau_Normal;
                    else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                              srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                              srch.Contains("2000_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                              srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                              srch.Contains("1500_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                              srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                              srch.Contains("1000_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                              srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                              srch.Contains("800_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                              srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                              srch.Contains("700_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                              srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                              srch.Contains("500_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                              srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                              srch.Contains("400_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                              srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                              srch.Contains("200_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                              srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                              srch.Contains("100_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                              srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                        first = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                              srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                        first = _EMu_DYTauTau_10to50;

                    if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                         srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                         srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                         srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                         srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                         srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                         srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                         srch.Contains("- Inf") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                              srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                              srch.Contains("TO_3000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                              srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                              srch.Contains("TO_2000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                              srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                              srch.Contains("TO_1500") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                              srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                              srch.Contains("TO_1000") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                              srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                              srch.Contains("TO_800") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                              srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                              srch.Contains("TO_700") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                              srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                              srch.Contains("TO_500") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                              srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                              srch.Contains("TO_400") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                              srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                              srch.Contains("TO_200") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                              srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                              srch.Contains("TO_100") )
                        last = _EMu_DYTauTau_50toInf;
                    else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                              srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                              srch.Contains("TO_50") )
                        last = _EMu_DYTauTau_10to50;
                    else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                              srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                              srch.Contains("TO_10") )
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

        else if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
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
                first = _EndOf_MuMu_MCsignal_Normal;
            else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                      srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                      srch.Contains("3000_TO") )
                first = _EndOf_MuMu_MCsignal_Normal;
            else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                      srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                      srch.Contains("2000_TO") )
                first =_MuMu_DY_2000to3000;
            else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                      srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                      srch.Contains("1500_TO") )
                first =_MuMu_DY_1500to2000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first =_MuMu_DY_1000to1500;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first =_MuMu_DY_800to1000;
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                      srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                      srch.Contains("700_TO") )
                first =_MuMu_DY_700to800;
            else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                      srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                      srch.Contains("500_TO") )
                first =_MuMu_DY_500to700;
            else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                      srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                      srch.Contains("400_TO") )
                first =_MuMu_DY_400to500;
            else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                      srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                      srch.Contains("200_TO") )
                first =_MuMu_DY_200to400;
            else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                      srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                      srch.Contains("100_TO") )
                first =_MuMu_DY_100to200;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                first =_MuMu_DY_50to100;
            else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                      srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                first = _MuMu_DY_10to50;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last =_MuMu_DY_2000to3000;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last =_MuMu_DY_2000to3000;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last =_MuMu_DY_1500to2000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last =_MuMu_DY_1000to1500;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last =_MuMu_DY_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last =_MuMu_DY_700to800;
            else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                      srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                      srch.Contains("TO_700") )
                last =_MuMu_DY_500to700;
            else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                      srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                      srch.Contains("TO_500") )
                last =_MuMu_DY_400to500;
            else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                      srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                      srch.Contains("TO_400") )
                last =_MuMu_DY_200to400;
            else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                      srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                      srch.Contains("TO_200") )
                last =_MuMu_DY_100to200;
            else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                      srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                      srch.Contains("TO_100") )
                last =_MuMu_DY_50to100;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _MuMu_DY_10to50;
            else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                      srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                      srch.Contains("TO_10") )
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
                first = _EndOf_EE_MCsignal_Normal;
            else if ( srch.Contains("3000to") || srch.Contains("3000To") || srch.Contains("3000TO") || srch.Contains("3000 to") || srch.Contains("3000 To") ||
                      srch.Contains("3000 TO") || srch.Contains("3000-") || srch.Contains("3000 -") || srch.Contains("3000_to") || srch.Contains("3000_To") ||
                      srch.Contains("3000_TO") )
                first = _EndOf_EE_MCsignal_Normal;
            else if ( srch.Contains("2000to") || srch.Contains("2000To") || srch.Contains("2000TO") || srch.Contains("2000 to") || srch.Contains("2000 To") ||
                      srch.Contains("2000 TO") || srch.Contains("2000-") || srch.Contains("2000 -") || srch.Contains("2000_to") || srch.Contains("2000_To") ||
                      srch.Contains("2000_TO") )
                first = _EE_DY_2000to3000;
            else if ( srch.Contains("1500to") || srch.Contains("1500To") || srch.Contains("1500TO") || srch.Contains("1500 to") || srch.Contains("1500 To") ||
                      srch.Contains("1500 TO") || srch.Contains("1500-") || srch.Contains("1500 -") || srch.Contains("1500_to") || srch.Contains("1500_To") ||
                      srch.Contains("1500_TO") )
                first = _EE_DY_1500to2000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _EE_DY_1000to1500;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _EE_DY_800to1000;
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") ||
                      srch.Contains("700 TO") || srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") ||
                      srch.Contains("700_TO") )
                first = _EE_DY_700to800;
            else if ( srch.Contains("500to") || srch.Contains("500To") || srch.Contains("500TO") || srch.Contains("500 to") || srch.Contains("500 To") ||
                      srch.Contains("500 TO") || srch.Contains("500-") || srch.Contains("500 -") || srch.Contains("500_to") || srch.Contains("500_To") ||
                      srch.Contains("500_TO") )
                first = _EE_DY_500to700;
            else if ( srch.Contains("400to") || srch.Contains("400To") || srch.Contains("400TO") || srch.Contains("400 to") || srch.Contains("400 To") ||
                      srch.Contains("400 TO") || srch.Contains("400-") || srch.Contains("400 -") || srch.Contains("400_to") || srch.Contains("400_To") ||
                      srch.Contains("400_TO") )
                first = _EE_DY_400to500;
            else if ( srch.Contains("200to") || srch.Contains("200To") || srch.Contains("200TO") || srch.Contains("200 to") || srch.Contains("200 To") ||
                      srch.Contains("200 TO") || srch.Contains("200-") || srch.Contains("200 -") || srch.Contains("200_to") || srch.Contains("200_To") ||
                      srch.Contains("200_TO") )
                first = _EE_DY_200to400;
            else if ( srch.Contains("100to") || srch.Contains("100To") || srch.Contains("100TO") || srch.Contains("100 to") || srch.Contains("100 To") ||
                      srch.Contains("100 TO") || srch.Contains("100-") || srch.Contains("100 -") || srch.Contains("100_to") || srch.Contains("100_To") ||
                      srch.Contains("100_TO") )
                first = _EE_DY_100to200;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") || srch.Contains("50_TO") )
                first = _EE_DY_50to100;
            else if ( srch.Contains("10to") || srch.Contains("10To") || srch.Contains("10TO") || srch.Contains("10 to") || srch.Contains("10 To") || srch.Contains("10 TO") ||
                      srch.Contains("10-") || srch.Contains("10 -") || srch.Contains("10_to") || srch.Contains("10_To") || srch.Contains("10_TO") )
                first = _EE_DY_10to50;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _EE_DY_2000to3000;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _EE_DY_2000to3000;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _EE_DY_1500to2000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _EE_DY_1000to1500;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _EE_DY_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _EE_DY_700to800;
            else if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") ||
                      srch.Contains("TO 700") || srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") ||
                      srch.Contains("TO_700") )
                last = _EE_DY_500to700;
            else if ( srch.Contains("to500") || srch.Contains("To500") || srch.Contains("TO500") || srch.Contains("to 500") || srch.Contains("To 500") ||
                      srch.Contains("TO 500") || srch.Contains("-500") || srch.Contains("- 500") || srch.Contains("to_500") || srch.Contains("To_500") ||
                      srch.Contains("TO_500") )
                last = _EE_DY_400to500;
            else if ( srch.Contains("to400") || srch.Contains("To400") || srch.Contains("TO400") || srch.Contains("to 400") || srch.Contains("To 400") ||
                      srch.Contains("TO 400") || srch.Contains("-400") || srch.Contains("- 400") || srch.Contains("to_400") || srch.Contains("To_400") ||
                      srch.Contains("TO_400") )
                last = _EE_DY_200to400;
            else if ( srch.Contains("to200") || srch.Contains("To200") || srch.Contains("TO200") || srch.Contains("to 200") || srch.Contains("To 200") ||
                      srch.Contains("TO 200") || srch.Contains("-200") || srch.Contains("- 200") || srch.Contains("to_200") || srch.Contains("To_200") ||
                      srch.Contains("TO_200") )
                last = _EE_DY_100to200;
            else if ( srch.Contains("to100") || srch.Contains("To100") || srch.Contains("TO100") || srch.Contains("to 100") || srch.Contains("To 100") ||
                      srch.Contains("TO 100") || srch.Contains("-100") || srch.Contains("- 100") || srch.Contains("to_100") || srch.Contains("To_100") ||
                      srch.Contains("TO_100") )
                last = _EE_DY_50to100;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _EE_DY_10to50;
            else if ( srch.Contains("to10") || srch.Contains("To10") || srch.Contains("TO10") || srch.Contains("to 10") || srch.Contains("To 10") ||
                      srch.Contains("TO 10") || srch.Contains("-10") || srch.Contains("- 10") || srch.Contains("to_10") || srch.Contains("To_10") ||
                      srch.Contains("TO_10") )
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

    else if ( srch.Contains("ttbar") || srch.Contains("tt") || srch.Contains("TTbar") || srch.Contains("TT") || srch.Contains("ttBAR") ||
         srch.Contains("TTBAR") || srch.Contains("TantiT") || srch.Contains("tantit") || srch.Contains("TANTIT") || srch.Contains("tANTIt") ||
         srch.Contains("TopAntiTop") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_MuMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_MuMu_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") || srch.Contains("700 TO") ||
                 srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") || srch.Contains("700_TO") )
                first = _MuMu_ttbar_700to1000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                      srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                      srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                      srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                      srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                      srch.Contains("Inf -") )
                first = _EndOf_MuMu_ttbar_Normal;
            else first = _MuMu_ttbar;

            if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") || srch.Contains("TO 700") ||
                 srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") || srch.Contains("TO_700") )
                last = _MuMu_ttbar;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _MuMu_ttbar_700to1000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _MuMu_ttbar_1000toInf;
            else if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                      srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                      srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                      srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                      srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                      srch.Contains("- Inf") )
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

        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_EE_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") || srch.Contains("700 TO") ||
                 srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") || srch.Contains("700_TO") )
                first = _EE_ttbar_700to1000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _EE_ttbar_1000toInf;
            else if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                      srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                      srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                      srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                      srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                      srch.Contains("Inf -") )
                first = _EndOf_EE_ttbar_Normal;
            else first = _EE_ttbar;

            if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") || srch.Contains("TO 700") ||
                 srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") || srch.Contains("TO_700") )
                last = _EE_ttbar;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _EE_ttbar_700to1000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _EE_ttbar_1000toInf;
            else if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                      srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                      srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                      srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                      srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                      srch.Contains("- Inf") )
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

        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_EMu_ttbar_Full);
                if ( notify == kTRUE ) cout << Procname[_EMu_ttbar_Full] << "." << endl;
            }
            // Checking for various intervals
            else if ( srch.Contains("700to") || srch.Contains("700To") || srch.Contains("700TO") || srch.Contains("700 to") || srch.Contains("700 To") || srch.Contains("700 TO") ||
                 srch.Contains("700-") || srch.Contains("700 -") || srch.Contains("700_to") || srch.Contains("700_To") || srch.Contains("700_TO") )
                first = _EMu_ttbar_700to1000;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("Infto") || srch.Contains("InfTo") || srch.Contains("InfTO") || srch.Contains("Inf to") || srch.Contains("Inf To") ||
                      srch.Contains("Inf TO") ||  srch.Contains("Inf_to") || srch.Contains("Inf_To") || srch.Contains("Inf_TO") ||
                      srch.Contains("infto") || srch.Contains("infTo") || srch.Contains("infTO") || srch.Contains("inf to") || srch.Contains("inf To") ||
                      srch.Contains("inf TO") ||  srch.Contains("inf_to") || srch.Contains("inf_To") || srch.Contains("inf_TO") ||
                      srch.Contains("INFto") || srch.Contains("INFTo") || srch.Contains("INFTO") || srch.Contains("INF to") || srch.Contains("INF To") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("inf-") || srch.Contains("INF-") || srch.Contains("In-f") || srch.Contains("inf -") || srch.Contains("INF -") ||
                      srch.Contains("Inf -") )
                first = _EndOf_EMu_ttbar_Normal;
            else first = _EMu_ttbar;

            if ( srch.Contains("to700") || srch.Contains("To700") || srch.Contains("TO700") || srch.Contains("to 700") || srch.Contains("To 700") || srch.Contains("TO 700") ||
                 srch.Contains("-700") || srch.Contains("- 700") || srch.Contains("to_700") || srch.Contains("To_700") || srch.Contains("TO_700") )
                last = _EMu_ttbar;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _EMu_ttbar_700to1000;
            else if ( srch.Contains("to1500") || srch.Contains("To1500") || srch.Contains("TO1500") || srch.Contains("to 1500") || srch.Contains("To 1500") ||
                      srch.Contains("TO 1500") || srch.Contains("-1500") || srch.Contains("- 1500") || srch.Contains("to_1500") || srch.Contains("To_1500") ||
                      srch.Contains("TO_1500") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("to2000") || srch.Contains("To2000") || srch.Contains("TO2000") || srch.Contains("to 2000") || srch.Contains("To 2000") ||
                      srch.Contains("TO 2000") || srch.Contains("-2000") || srch.Contains("- 2000") || srch.Contains("to_2000") || srch.Contains("To_2000") ||
                      srch.Contains("TO_2000") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("to3000") || srch.Contains("To3000") || srch.Contains("TO3000") || srch.Contains("to 3000") || srch.Contains("To 3000") ||
                      srch.Contains("TO 3000") || srch.Contains("-3000") || srch.Contains("- 3000") || srch.Contains("to_3000") || srch.Contains("To_3000") ||
                      srch.Contains("TO_3000") )
                last = _EMu_ttbar_1000toInf;
            else if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                      srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                      srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                      srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                      srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                      srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                      srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                      srch.Contains("- Inf") )
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

    else if ( srch.Contains("tW") || srch.Contains("TW") || srch.Contains("SingleTop") || srch.Contains("singleTop") || srch.Contains("singletop") ||
         srch.Contains("SingleTOP") || srch.Contains("singleTOP") || srch.Contains("SINGLETOP") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_tW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_tW] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_tW);
            if ( notify == kTRUE ) cout << Procname[_EE_tW] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_tW);
            if ( notify == kTRUE ) cout << Procname[_EMu_tW] << "." << endl;
        }

    }

    else if ( srch.Contains("tbarW") || srch.Contains("TbarW") || srch.Contains("TBarW") || srch.Contains("tBarW") || srch.Contains("TBARW") ||
              srch.Contains("SingleAntiTop") || srch.Contains("SingleAntitop") || srch.Contains("singleAntiTop") || srch.Contains("singleAntitop") ||
              srch.Contains("SingleAntiTOP") || srch.Contains("singleAntiTOP") || srch.Contains("SingleANTITOP") || srch.Contains("singleANTITOP") ||
              srch.Contains("SINGLEANTITOP") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_tbarW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_tbarW] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_tbarW);
            if ( notify == kTRUE ) cout << Procname[_EE_tbarW] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_tbarW);
            if ( notify == kTRUE ) cout << Procname[_EMu_tbarW] << "." << endl;
        }
    }

    else if ( srch.Contains("ZZ") || srch.Contains("zz") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_ZZ);
            if ( notify == kTRUE ) cout << Procname[_MuMu_ZZ] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_ZZ);
            if ( notify == kTRUE ) cout << Procname[_EE_ZZ] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_ZZ);
            if ( notify == kTRUE ) cout << Procname[_EMu_ZZ] << "." << endl;
        }
    }

    else if ( srch.Contains("WZ") || srch.Contains("wz") || srch.Contains("ZW") || srch.Contains("zw") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_WZ);
            if ( notify == kTRUE ) cout << Procname[_MuMu_WZ] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_WZ);
            if ( notify == kTRUE ) cout << Procname[_EE_WZ] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_WZ);
            if ( notify == kTRUE ) cout << Procname[_EMu_WZ] << "." << endl;
        }
    }

    else if ( srch.Contains("WW") || srch.Contains("ww") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_WW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_WW] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_WW);
            if ( notify == kTRUE ) cout << Procname[_EE_WW] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_WW);
            if ( notify == kTRUE ) cout << Procname[_EMu_WW] << "." << endl;
        }
    }

    else if ( srch.Contains("DIBOSON") || srch.Contains("diboson") || srch.Contains("Diboson") || srch.Contains("diBoson"))
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_ZZ);
            Result.push_back(_MuMu_WZ);
            Result.push_back(_MuMu_WW);
            if ( notify == kTRUE ) cout << Procname[_MuMu_ZZ] << ", " << Procname[_MuMu_WZ] << ", " << Procname[_MuMu_WW] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_ZZ);
            Result.push_back(_EE_WZ);
            Result.push_back(_EE_WW);
            if ( notify == kTRUE ) cout << Procname[_EE_ZZ] << ", " << Procname[_EE_WZ] << ", " << Procname[_EE_WW] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_ZZ);
            Result.push_back(_EMu_WZ);
            Result.push_back(_EMu_WW);
            if ( notify == kTRUE ) cout << Procname[_EMu_ZZ] << ", " << Procname[_EMu_WZ] << ", " << Procname[_EMu_WW] << "." << endl;
        }
    }

    else if ( srch.Contains("VVnST") || srch.Contains("VVNST") || srch.Contains("vvnst") || srch.Contains("VVnst") || srch.Contains("vvnST") ||
              srch.Contains("VV") || srch.Contains("vv") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_VVnST);
            if ( notify == kTRUE ) cout << Procname[_MuMu_VVnST] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_VVnST);
            if ( notify == kTRUE ) cout << Procname[_EE_VVnST] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_VVnST);
            if ( notify == kTRUE ) cout << Procname[_EMu_VVnST] << "." << endl;
        }
    }

    else if ( srch.Contains("WJets") || srch.Contains("Wjets") || srch.Contains("wjets") || srch.Contains("WJETS") || srch.Contains("W+jets") ||
              srch.Contains("W+Jets") || srch.Contains("W+JETS") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_WJets);
            if ( notify == kTRUE ) cout << Procname[_MuMu_WJets] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_WJets);
            if ( notify == kTRUE ) cout << Procname[_EE_WJets] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_WJets);
            if ( notify == kTRUE ) cout << Procname[_EMu_WJets] << "." << endl;
        }
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
                first = _EndOf_MuMu_QCDMuEnriched_Normal;
            else if ( srch.Contains("1000to") || srch.Contains("1000To") || srch.Contains("1000TO") || srch.Contains("1000 to") || srch.Contains("1000 To") ||
                      srch.Contains("1000 TO") || srch.Contains("1000-") || srch.Contains("1000 -") || srch.Contains("1000_to") || srch.Contains("1000_To") ||
                      srch.Contains("1000_TO") )
                first = _MuMu_QCDMuEnriched_1000toInf;
            else if ( srch.Contains("800to") || srch.Contains("800To") || srch.Contains("800TO") || srch.Contains("800 to") || srch.Contains("800 To") ||
                      srch.Contains("800 TO") || srch.Contains("800-") || srch.Contains("800 -") || srch.Contains("800_to") || srch.Contains("800_To") ||
                      srch.Contains("800_TO") )
                first = _MuMu_QCDMuEnriched_800to1000;
            else if ( srch.Contains("600to") || srch.Contains("600To") || srch.Contains("600TO") || srch.Contains("600 to") || srch.Contains("600 To") ||
                      srch.Contains("600 TO") || srch.Contains("600-") || srch.Contains("600 -") || srch.Contains("600_to") || srch.Contains("600_To") ||
                      srch.Contains("600_TO") )
                first = _MuMu_QCDMuEnriched_600to800;
            else if ( srch.Contains("470to") || srch.Contains("470To") || srch.Contains("470TO") || srch.Contains("470 to") || srch.Contains("470 To") ||
                      srch.Contains("470 TO") || srch.Contains("470-") || srch.Contains("470 -") || srch.Contains("470_to") || srch.Contains("470_To") ||
                      srch.Contains("470_TO") )
                first = _MuMu_QCDMuEnriched_470to600;
            else if ( srch.Contains("300to") || srch.Contains("300To") || srch.Contains("300TO") || srch.Contains("300 to") || srch.Contains("300 To") ||
                      srch.Contains("300 TO") || srch.Contains("300-") || srch.Contains("300 -") || srch.Contains("300_to") || srch.Contains("300_To") ||
                      srch.Contains("300_TO") )
                first = _MuMu_QCDMuEnriched_300to470;
            else if ( srch.Contains("170to") || srch.Contains("170To") || srch.Contains("170TO") || srch.Contains("170 to") || srch.Contains("170 To") ||
                      srch.Contains("170 TO") || srch.Contains("170-") || srch.Contains("170 -") || srch.Contains("170_to") || srch.Contains("170_To") ||
                      srch.Contains("170_TO") )
                first = _MuMu_QCDMuEnriched_170to300;
            else if ( srch.Contains("120to") || srch.Contains("120To") || srch.Contains("120TO") || srch.Contains("120 to") || srch.Contains("120 To") ||
                      srch.Contains("120 TO") || srch.Contains("120-") || srch.Contains("120 -") || srch.Contains("120_to") || srch.Contains("120_To") ||
                      srch.Contains("120_TO") )
                first = _MuMu_QCDMuEnriched_120to170;
            else if ( srch.Contains("80to") || srch.Contains("80To") || srch.Contains("80TO") || srch.Contains("80 to") || srch.Contains("80 To") ||
                      srch.Contains("80 TO") || srch.Contains("80-") || srch.Contains("80 -") || srch.Contains("80_to") || srch.Contains("80_To") ||
                      srch.Contains("80_TO") )
                first = _MuMu_QCDMuEnriched_80to120;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") ||
                      srch.Contains("50_TO") )
                first = _MuMu_QCDMuEnriched_50to80;
            else if ( srch.Contains("30to") || srch.Contains("30To") || srch.Contains("30TO") || srch.Contains("30 to") || srch.Contains("30 To") ||
                      srch.Contains("30 TO") || srch.Contains("30-") || srch.Contains("30 -") || srch.Contains("30_to") || srch.Contains("30_To") ||
                      srch.Contains("30_TO") )
                first = _MuMu_QCDMuEnriched_30to50;
            else if ( srch.Contains("20to") || srch.Contains("20To") || srch.Contains("20TO") || srch.Contains("20 to") || srch.Contains("20 To") ||
                      srch.Contains("20 TO") || srch.Contains("20-") || srch.Contains("20 -") || srch.Contains("20_to") || srch.Contains("20_To") ||
                      srch.Contains("20_TO") )
                first = _MuMu_QCDMuEnriched_20to30;
            else if ( srch.Contains("15to") || srch.Contains("15To") || srch.Contains("15TO") || srch.Contains("15 to") || srch.Contains("15 To") ||
                      srch.Contains("15 TO") || srch.Contains("15-") || srch.Contains("15 -") || srch.Contains("15_to") || srch.Contains("15_To") ||
                      srch.Contains("15_TO") )
                first = _MuMu_QCDMuEnriched_15to20;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _MuMu_QCDMuEnriched_1000toInf;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _MuMu_QCDMuEnriched_800to1000;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _MuMu_QCDMuEnriched_600to800;
            else if ( srch.Contains("to600") || srch.Contains("To600") || srch.Contains("TO600") || srch.Contains("to 600") || srch.Contains("To 600") ||
                      srch.Contains("TO 600") || srch.Contains("-600") || srch.Contains("- 600") || srch.Contains("to_600") || srch.Contains("To_600") ||
                      srch.Contains("TO_600") )
                last = _MuMu_QCDMuEnriched_470to600;
            else if ( srch.Contains("to470") || srch.Contains("To470") || srch.Contains("TO470") || srch.Contains("to 470") || srch.Contains("To 470") ||
                      srch.Contains("TO 470") || srch.Contains("-470") || srch.Contains("- 470") || srch.Contains("to_470") || srch.Contains("To_470") ||
                      srch.Contains("TO_470") )
                last = _MuMu_QCDMuEnriched_300to470;
            else if ( srch.Contains("to300") || srch.Contains("To300") || srch.Contains("TO300") || srch.Contains("to 300") || srch.Contains("To 300") ||
                      srch.Contains("TO 300") || srch.Contains("-300") || srch.Contains("- 300") || srch.Contains("to_300") || srch.Contains("To_300") ||
                      srch.Contains("TO_300") )
                last = _MuMu_QCDMuEnriched_170to300;
            else if ( srch.Contains("to170") || srch.Contains("To170") || srch.Contains("TO170") || srch.Contains("to 170") || srch.Contains("To 170") ||
                      srch.Contains("TO 170") || srch.Contains("-170") || srch.Contains("- 170") || srch.Contains("to_170") || srch.Contains("To_170") ||
                      srch.Contains("TO_170") )
                last = _MuMu_QCDMuEnriched_120to170;
            else if ( srch.Contains("to120") || srch.Contains("To120") || srch.Contains("TO120") || srch.Contains("to 120") || srch.Contains("To 120") ||
                      srch.Contains("TO 120") || srch.Contains("-120") || srch.Contains("- 120") || srch.Contains("to_120") || srch.Contains("To_120") ||
                      srch.Contains("TO_120") )
                last = _MuMu_QCDMuEnriched_80to120;
            else if ( srch.Contains("to80") || srch.Contains("To80") || srch.Contains("TO80") || srch.Contains("to 80") || srch.Contains("To 80") ||
                      srch.Contains("TO 80") || srch.Contains("-80") || srch.Contains("- 80") || srch.Contains("to_80") || srch.Contains("To_80") ||
                      srch.Contains("TO_80") )
                last = _MuMu_QCDMuEnriched_50to80;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _MuMu_QCDMuEnriched_30to50;
            else if ( srch.Contains("to30") || srch.Contains("To30") || srch.Contains("TO30") || srch.Contains("to 30") || srch.Contains("To 30") ||
                      srch.Contains("TO 30") || srch.Contains("-30") || srch.Contains("- 30") || srch.Contains("to_30") || srch.Contains("To_30") ||
                      srch.Contains("TO_30") )
                last = _MuMu_QCDMuEnriched_20to30;
            else if ( srch.Contains("to20") || srch.Contains("To20") || srch.Contains("TO20") || srch.Contains("to 20") || srch.Contains("To 20") ||
                      srch.Contains("TO 20") || srch.Contains("-20") || srch.Contains("- 20") || srch.Contains("to_20") || srch.Contains("To_20") ||
                      srch.Contains("TO_20") )
                last = _MuMu_QCDMuEnriched_15to20;
            else if ( srch.Contains("to15") || srch.Contains("To15") || srch.Contains("TO15") || srch.Contains("to 15") || srch.Contains("To 15") ||
                      srch.Contains("TO 15") || srch.Contains("-15") || srch.Contains("- 15") || srch.Contains("to_15") || srch.Contains("To_15") ||
                      srch.Contains("TO_15") )
                last = _EndOf_MuMu_WJets;

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
        else if ( srch.Contains("EM") || srch.Contains("em") || srch.Contains("Em") || srch.Contains("eM") || srch.Contains("EE") || srch.Contains("ee")
                  || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") || srch.Contains("Dielectron") || srch.Contains("diElectron")
                  || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
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
                first = _EndOf_EE_QCDEMEnriched_Normal;
            else if ( srch.Contains("300to") || srch.Contains("300To") || srch.Contains("300TO") || srch.Contains("300 to") || srch.Contains("300 To") ||
                      srch.Contains("300 TO") || srch.Contains("300-") || srch.Contains("300 -") || srch.Contains("300_to") || srch.Contains("300_To") ||
                      srch.Contains("300_TO") )
                first = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("170to") || srch.Contains("170To") || srch.Contains("170TO") || srch.Contains("170 to") || srch.Contains("170 To") ||
                      srch.Contains("170 TO") || srch.Contains("170-") || srch.Contains("170 -") || srch.Contains("170_to") || srch.Contains("170_To") ||
                      srch.Contains("170_TO") )
                first = _EE_QCDEMEnriched_170to300;
            else if ( srch.Contains("120to") || srch.Contains("120To") || srch.Contains("120TO") || srch.Contains("120 to") || srch.Contains("120 To") ||
                      srch.Contains("120 TO") || srch.Contains("120-") || srch.Contains("120 -") || srch.Contains("120_to") || srch.Contains("120_To") ||
                      srch.Contains("120_TO") )
                first = _EE_QCDEMEnriched_120to170;
            else if ( srch.Contains("80to") || srch.Contains("80To") || srch.Contains("80TO") || srch.Contains("80 to") || srch.Contains("80 To") ||
                      srch.Contains("80 TO") || srch.Contains("80-") || srch.Contains("80 -") || srch.Contains("80_to") || srch.Contains("80_To") ||
                      srch.Contains("80_TO") )
                first = _EE_QCDEMEnriched_80to120;
            else if ( srch.Contains("50to") || srch.Contains("50To") || srch.Contains("50TO") || srch.Contains("50 to") || srch.Contains("50 To") ||
                      srch.Contains("50 TO") || srch.Contains("50-") || srch.Contains("50 -") || srch.Contains("50_to") || srch.Contains("50_To") ||
                      srch.Contains("50_TO") )
                first = _EE_QCDEMEnriched_50to80;
            else if ( srch.Contains("30to") || srch.Contains("30To") || srch.Contains("30TO") || srch.Contains("30 to") || srch.Contains("30 To") ||
                      srch.Contains("30 TO") || srch.Contains("30-") || srch.Contains("30 -") || srch.Contains("30_to") || srch.Contains("30_To") ||
                      srch.Contains("30_TO") )
                first = _EE_QCDEMEnriched_30to50;
            else if ( srch.Contains("20to") || srch.Contains("20To") || srch.Contains("20TO") || srch.Contains("20 to") || srch.Contains("20 To") ||
                      srch.Contains("20 TO") || srch.Contains("20-") || srch.Contains("20 -") || srch.Contains("20_to") || srch.Contains("20_To") ||
                      srch.Contains("20_TO") )
                first = _EE_QCDEMEnriched_20to30;

            if ( srch.Contains("toInf") || srch.Contains("ToInf") || srch.Contains("TOInf") || srch.Contains("to Inf") || srch.Contains("To Inf") ||
                 srch.Contains("TO Inf") ||  srch.Contains("to_Inf") || srch.Contains("To_Inf") || srch.Contains("TO_Inf") ||
                 srch.Contains("toinf") || srch.Contains("Toinf") || srch.Contains("TOinf") || srch.Contains("to inf") || srch.Contains("To inf") ||
                 srch.Contains("TO inf") ||  srch.Contains("to_inf") || srch.Contains("To_inf") || srch.Contains("TO_inf") ||
                 srch.Contains("toINF") || srch.Contains("ToINF") || srch.Contains("TOINF") || srch.Contains("to INF") || srch.Contains("To INF") ||
                 srch.Contains("TO INF") ||  srch.Contains("to_INF") || srch.Contains("To_INF") || srch.Contains("TO_INF") ||
                 srch.Contains("-inf") || srch.Contains("-INF") || srch.Contains("-Inf") || srch.Contains("- inf") || srch.Contains("- INF") ||
                 srch.Contains("- Inf") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("to1000") || srch.Contains("To1000") || srch.Contains("TO1000") || srch.Contains("to 1000") || srch.Contains("To 1000") ||
                      srch.Contains("TO 1000") || srch.Contains("-1000") || srch.Contains("- 1000") || srch.Contains("to_1000") || srch.Contains("To_1000") ||
                      srch.Contains("TO_1000") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("to800") || srch.Contains("To800") || srch.Contains("TO800") || srch.Contains("to 800") || srch.Contains("To 800") ||
                      srch.Contains("TO 800") || srch.Contains("-800") || srch.Contains("- 800") || srch.Contains("to_800") || srch.Contains("To_800") ||
                      srch.Contains("TO_800") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("to600") || srch.Contains("To600") || srch.Contains("TO600") || srch.Contains("to 600") || srch.Contains("To 600") ||
                      srch.Contains("TO 600") || srch.Contains("-600") || srch.Contains("- 600") || srch.Contains("to_600") || srch.Contains("To_600") ||
                      srch.Contains("TO_600") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("to470") || srch.Contains("To470") || srch.Contains("TO470") || srch.Contains("to 470") || srch.Contains("To 470") ||
                      srch.Contains("TO 470") || srch.Contains("-470") || srch.Contains("- 470") || srch.Contains("to_470") || srch.Contains("To_470") ||
                      srch.Contains("TO_470") )
                last = _EE_QCDEMEnriched_300toInf;
            else if ( srch.Contains("to300") || srch.Contains("To300") || srch.Contains("TO300") || srch.Contains("to 300") || srch.Contains("To 300") ||
                      srch.Contains("TO 300") || srch.Contains("-300") || srch.Contains("- 300") || srch.Contains("to_300") || srch.Contains("To_300") ||
                      srch.Contains("TO_300") )
                last = _EE_QCDEMEnriched_170to300;
            else if ( srch.Contains("to170") || srch.Contains("To170") || srch.Contains("TO170") || srch.Contains("to 170") || srch.Contains("To 170") ||
                      srch.Contains("TO 170") || srch.Contains("-170") || srch.Contains("- 170") || srch.Contains("to_170") || srch.Contains("To_170") ||
                      srch.Contains("TO_170") )
                last = _EE_QCDEMEnriched_120to170;
            else if ( srch.Contains("to120") || srch.Contains("To120") || srch.Contains("TO120") || srch.Contains("to 120") || srch.Contains("To 120") ||
                      srch.Contains("TO 120") || srch.Contains("-120") || srch.Contains("- 120") || srch.Contains("to_120") || srch.Contains("To_120") ||
                      srch.Contains("TO_120") )
                last = _EE_QCDEMEnriched_80to120;
            else if ( srch.Contains("to80") || srch.Contains("To80") || srch.Contains("TO80") || srch.Contains("to 80") || srch.Contains("To 80") ||
                      srch.Contains("TO 80") || srch.Contains("-80") || srch.Contains("- 80") || srch.Contains("to_80") || srch.Contains("To_80") ||
                      srch.Contains("TO_80") )
                last = _EE_QCDEMEnriched_50to80;
            else if ( srch.Contains("to50") || srch.Contains("To50") || srch.Contains("TO50") || srch.Contains("to 50") || srch.Contains("To 50") ||
                      srch.Contains("TO 50") || srch.Contains("-50") || srch.Contains("- 50") || srch.Contains("to_50") || srch.Contains("To_50") ||
                      srch.Contains("TO_50") )
                last = _EE_QCDEMEnriched_30to50;
            else if ( srch.Contains("to30") || srch.Contains("To30") || srch.Contains("TO30") || srch.Contains("to 30") || srch.Contains("To 30") ||
                      srch.Contains("TO 30") || srch.Contains("-30") || srch.Contains("- 30") || srch.Contains("to_30") || srch.Contains("To_30") ||
                      srch.Contains("TO_30") )
                last = _EE_QCDEMEnriched_20to30;
            else if ( srch.Contains("to20") || srch.Contains("To20") || srch.Contains("TO20") || srch.Contains("to 20") || srch.Contains("To 20") ||
                      srch.Contains("TO 20") || srch.Contains("-20") || srch.Contains("- 20") || srch.Contains("to_20") || srch.Contains("To_20") ||
                      srch.Contains("TO_20") )
                last = _EndOf_EE_WJets;

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

    else if (  srch.Contains("Double") || srch.Contains("double") || srch.Contains("DOUBLE") )
    {
        if ( srch.Contains("Eg") || srch.Contains("eg") || srch.Contains("EG") || srch.Contains("eG") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                 srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                 srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                 srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                first = _EE_DoubleEG_B;
            else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                      srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                      srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                      srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                      srch.Contains("c_To") || srch.Contains("c_TO") )
                first = _EE_DoubleEG_C;
            else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                      srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                      srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                      srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                      srch.Contains("d_To") || srch.Contains("d_TO") )
                first = _EE_DoubleEG_D;
            else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                      srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                      srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                      srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                      srch.Contains("e_To") || srch.Contains("e_TO") )
                first = _EE_DoubleEG_E;
            else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                      srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                      srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                      srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                      srch.Contains("f_To") || srch.Contains("f_TO") )
                first = _EE_DoubleEG_F;
            else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                      srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                      srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                      srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                      srch.Contains("g_To") || srch.Contains("g_TO") )
                first = _EE_DoubleEG_G;
            else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                      srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                      srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                      srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                      srch.Contains("h_To") || srch.Contains("h_TO") )
                first = _EE_DoubleEG_H;

            if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                 srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                 srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                 srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                 srch.Contains("To_b") || srch.Contains("TO_b") )
                last = _EE_DoubleEG_B;
            else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                      srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                      srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                      srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                      srch.Contains("To_c") || srch.Contains("TO_c") )
                last = _EE_DoubleEG_C;
            else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                      srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                      srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                      srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                      srch.Contains("To_d") || srch.Contains("TO_d") )
                last = _EE_DoubleEG_D;
            else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                      srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                      srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                      srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                      srch.Contains("To_e") || srch.Contains("TO_e") )
                last = _EE_DoubleEG_E;
            else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                      srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                      srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                      srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                      srch.Contains("To_f") || srch.Contains("TO_f") )
                last = _EE_DoubleEG_F;
            else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                      srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                      srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                      srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                      srch.Contains("To_g") || srch.Contains("TO_g") )
                last = _EE_DoubleEG_G;
            else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                      srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                      srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                      srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                      srch.Contains("To_h") || srch.Contains("TO_h") )
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
            else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_EE_DoubleEG_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                Result.push_back(_EE_DoubleEG_B);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                Result.push_back(_EE_DoubleEG_C);
                if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
                {
                    Result.push_back(_EE_DoubleEG_E);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
                {
                    Result.push_back(_EE_DoubleEG_F);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
                {
                    Result.push_back(_EE_DoubleEG_G);
                    if ( notify == kTRUE ) cout << Procname[_EE_DoubleEG_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
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

    else if ( srch.Contains("Single") || srch.Contains("single") || srch.Contains("SINGLE") )
    {
        if ( srch.Contains("Muon") || srch.Contains("muon") || srch.Contains("MUON") )
        {
            if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
                 srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
            {
                // Checking if search contains intervals
                if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                     srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                     srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                     srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                    first = _MuMu_SingleMuon_B;
                else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                          srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                          srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                          srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                          srch.Contains("c_To") || srch.Contains("c_TO") )
                    first = _MuMu_SingleMuon_C;
                else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                          srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                          srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                          srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                          srch.Contains("d_To") || srch.Contains("d_TO") )
                    first = _MuMu_SingleMuon_D;
                else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                          srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                          srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                          srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                          srch.Contains("e_To") || srch.Contains("e_TO") )
                    first = _MuMu_SingleMuon_E;
                else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                          srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                          srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                          srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                          srch.Contains("f_To") || srch.Contains("f_TO") )
                    first = _MuMu_SingleMuon_F;
                else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                          srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                          srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                          srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                          srch.Contains("g_To") || srch.Contains("g_TO") )
                    first = _MuMu_SingleMuon_G;
                else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                          srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                          srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                          srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                          srch.Contains("h_To") || srch.Contains("h_TO") )
                    first = _MuMu_SingleMuon_H;

                if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                     srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                     srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                     srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                     srch.Contains("To_b") || srch.Contains("TO_b") )
                    last = _MuMu_SingleMuon_B;
                else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                          srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                          srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                          srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                          srch.Contains("To_c") || srch.Contains("TO_c") )
                    last = _MuMu_SingleMuon_C;
                else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                          srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                          srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                          srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                          srch.Contains("To_d") || srch.Contains("TO_d") )
                    last = _MuMu_SingleMuon_D;
                else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                          srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                          srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                          srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                          srch.Contains("To_e") || srch.Contains("TO_e") )
                    last = _MuMu_SingleMuon_E;
                else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                          srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                          srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                          srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                          srch.Contains("To_f") || srch.Contains("TO_f") )
                    last = _MuMu_SingleMuon_F;
                else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                          srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                          srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                          srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                          srch.Contains("To_g") || srch.Contains("TO_g") )
                    last = _MuMu_SingleMuon_G;
                else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                          srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                          srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                          srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                          srch.Contains("To_h") || srch.Contains("TO_h") )
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
                else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
                {
                    Result.push_back(_MuMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_Full] << "." << endl;
                }
                else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
                {
                    Result.push_back(_MuMu_SingleMuon_B);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_B] << "." << endl;
                }
                else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
                {
                    Result.push_back(_MuMu_SingleMuon_C);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_C] << "." << endl;
                }
                else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
                {
                    Result.push_back(_MuMu_SingleMuon_D);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_D] << "." << endl;
                }
                else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
                {
                    Result.push_back(_MuMu_SingleMuon_E);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
                {
                    Result.push_back(_MuMu_SingleMuon_F);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
                {
                    Result.push_back(_MuMu_SingleMuon_G);
                    if ( notify == kTRUE ) cout << Procname[_MuMu_SingleMuon_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
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

            else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                      srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
            {
                // Checking if search contains intervals
                if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                     srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                     srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                     srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                    first = _EMu_SingleMuon_B;
                else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                          srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                          srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                          srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                          srch.Contains("c_To") || srch.Contains("c_TO") )
                    first = _EMu_SingleMuon_C;
                else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                          srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                          srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                          srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                          srch.Contains("d_To") || srch.Contains("d_TO") )
                    first = _EMu_SingleMuon_D;
                else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                          srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                          srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                          srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                          srch.Contains("e_To") || srch.Contains("e_TO") )
                    first = _EMu_SingleMuon_E;
                else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                          srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                          srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                          srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                          srch.Contains("f_To") || srch.Contains("f_TO") )
                    first = _EMu_SingleMuon_F;
                else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                          srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                          srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                          srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                          srch.Contains("g_To") || srch.Contains("g_TO") )
                    first = _EMu_SingleMuon_G;
                else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                          srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                          srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                          srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                          srch.Contains("h_To") || srch.Contains("h_TO") )
                    first = _EMu_SingleMuon_H;

                if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                     srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                     srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                     srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                     srch.Contains("To_b") || srch.Contains("TO_b") )
                    last = _EMu_SingleMuon_B;
                else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                          srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                          srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                          srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                          srch.Contains("To_c") || srch.Contains("TO_c") )
                    last = _EMu_SingleMuon_C;
                else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                          srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                          srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                          srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                          srch.Contains("To_d") || srch.Contains("TO_d") )
                    last = _EMu_SingleMuon_D;
                else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                          srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                          srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                          srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                          srch.Contains("To_e") || srch.Contains("TO_e") )
                    last = _EMu_SingleMuon_E;
                else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                          srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                          srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                          srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                          srch.Contains("To_f") || srch.Contains("TO_f") )
                    last = _EMu_SingleMuon_F;
                else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                          srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                          srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                          srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                          srch.Contains("To_g") || srch.Contains("TO_g") )
                    last = _EMu_SingleMuon_G;
                else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                          srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                          srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                          srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                          srch.Contains("To_h") || srch.Contains("TO_h") )
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
                else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
                {
                    Result.push_back(_EMu_SingleMuon_Full);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_Full] << "." << endl;
                }
                else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
                {
                    Result.push_back(_EMu_SingleMuon_B);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_B] << "." << endl;
                }
                else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
                {
                    Result.push_back(_EMu_SingleMuon_C);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_C] << "." << endl;
                }
                else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
                {
                    Result.push_back(_EMu_SingleMuon_D);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_D] << "." << endl;
                }
                else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
                {
                    Result.push_back(_EMu_SingleMuon_E);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_E] << "." << endl;
                }
                else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
                {
                    Result.push_back(_EMu_SingleMuon_F);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_F] << "." << endl;
                }
                else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
                {
                    Result.push_back(_EMu_SingleMuon_G);
                    if ( notify == kTRUE ) cout << Procname[_EMu_SingleMuon_G] << "." << endl;
                }
                else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
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

        if ( srch.Contains("Elec") || srch.Contains("elec") || srch.Contains("ELEC") )
        {
            // Checking if search contains intervals
            if ( srch.Contains("Bto") || srch.Contains("BTo") || srch.Contains("BTO") || srch.Contains("B to") || srch.Contains("B To") ||
                 srch.Contains("B TO") || srch.Contains("B-") || srch.Contains("B -") || srch.Contains("B_to") || srch.Contains("B_To") || srch.Contains("B_TO") ||
                 srch.Contains("bto") || srch.Contains("bTo") || srch.Contains("bTO") || srch.Contains("b to") || srch.Contains("b To") || srch.Contains("c TO") ||
                 srch.Contains("b-") || srch.Contains("b -") || srch.Contains("b_to") || srch.Contains("b_To") || srch.Contains("b_TO"))
                first = _EE_SingleElectron_B;
            else if ( srch.Contains("Cto") || srch.Contains("CTo") || srch.Contains("CTO") || srch.Contains("C to") || srch.Contains("C To") ||
                      srch.Contains("C TO") || srch.Contains("C-") || srch.Contains("C -") || srch.Contains("C_to") || srch.Contains("C_To") ||
                      srch.Contains("C_TO") || srch.Contains("cto") || srch.Contains("cTo") || srch.Contains("cTO") || srch.Contains("c to") ||
                      srch.Contains("c To") || srch.Contains("c TO") || srch.Contains("c-") || srch.Contains("c -") || srch.Contains("c_to") ||
                      srch.Contains("c_To") || srch.Contains("c_TO") )
                first = _EE_SingleElectron_C;
            else if ( srch.Contains("Dto") || srch.Contains("DTo") || srch.Contains("DTO") || srch.Contains("D to") || srch.Contains("D To") ||
                      srch.Contains("D TO") || srch.Contains("D-") || srch.Contains("D -") || srch.Contains("D_to") || srch.Contains("D_To") ||
                      srch.Contains("D_TO") || srch.Contains("dto") || srch.Contains("dTo") || srch.Contains("dTO") || srch.Contains("d to") ||
                      srch.Contains("d To") || srch.Contains("d TO") || srch.Contains("d-") || srch.Contains("d -") || srch.Contains("d_to") ||
                      srch.Contains("d_To") || srch.Contains("d_TO") )
                first = _EE_SingleElectron_D;
            else if ( srch.Contains("Eto") || srch.Contains("ETo") || srch.Contains("ETO") || srch.Contains("E to") || srch.Contains("E To") ||
                      srch.Contains("E TO") || srch.Contains("E-") || srch.Contains("E -") || srch.Contains("E_to") || srch.Contains("E_To") ||
                      srch.Contains("E_TO") || srch.Contains("eto") || srch.Contains("eTo") || srch.Contains("eTO") || srch.Contains("e to") ||
                      srch.Contains("e To") || srch.Contains("e TO") || srch.Contains("e-") || srch.Contains("e -") || srch.Contains("e_to") ||
                      srch.Contains("e_To") || srch.Contains("e_TO") )
                first = _EE_SingleElectron_E;
            else if ( srch.Contains("Fto") || srch.Contains("FTo") || srch.Contains("FTO") || srch.Contains("F to") || srch.Contains("F To") ||
                      srch.Contains("F TO") || srch.Contains("F-") || srch.Contains("F -") || srch.Contains("F_to") || srch.Contains("F_To") ||
                      srch.Contains("F_TO") || srch.Contains("fto") || srch.Contains("fTo") || srch.Contains("fTO") || srch.Contains("f to") ||
                      srch.Contains("f To") || srch.Contains("f TO") || srch.Contains("f-") || srch.Contains("f -") || srch.Contains("f_to") ||
                      srch.Contains("f_To") || srch.Contains("f_TO") )
                first = _EE_SingleElectron_F;
            else if ( srch.Contains("Gto") || srch.Contains("GTo") || srch.Contains("GTO") || srch.Contains("G to") || srch.Contains("G To") ||
                      srch.Contains("G TO") || srch.Contains("G-") || srch.Contains("G -") || srch.Contains("G_to") || srch.Contains("G_To") ||
                      srch.Contains("G_TO") || srch.Contains("gto") || srch.Contains("gTo") || srch.Contains("gTO") || srch.Contains("g to") ||
                      srch.Contains("g To") || srch.Contains("g TO") || srch.Contains("g-") || srch.Contains("g -") || srch.Contains("g_to") ||
                      srch.Contains("g_To") || srch.Contains("g_TO") )
                first = _EE_SingleElectron_G;
            else if ( srch.Contains("Hto") || srch.Contains("HTo") || srch.Contains("HTO") || srch.Contains("H to") || srch.Contains("H To") ||
                      srch.Contains("H TO") || srch.Contains("H-") || srch.Contains("H -") || srch.Contains("H_to") || srch.Contains("H_To") ||
                      srch.Contains("H_TO") || srch.Contains("hto") || srch.Contains("hTo") || srch.Contains("hTO") || srch.Contains("h to") ||
                      srch.Contains("h To") || srch.Contains("h TO") || srch.Contains("h-") || srch.Contains("h -") || srch.Contains("h_to") ||
                      srch.Contains("h_To") || srch.Contains("h_TO") )
                first = _EE_SingleElectron_H;

            if ( srch.Contains("toB") || srch.Contains("ToB") || srch.Contains("TOB") || srch.Contains("to B") || srch.Contains("To B") ||
                 srch.Contains("TO B") || srch.Contains("-B") || srch.Contains("- B") || srch.Contains("to_B") || srch.Contains("To_B") ||
                 srch.Contains("TO_B") || srch.Contains("tob") || srch.Contains("Tob") || srch.Contains("TOb") || srch.Contains("to b") ||
                 srch.Contains("To b") || srch.Contains("TO b") || srch.Contains("-b") || srch.Contains("- b") || srch.Contains("to_b") ||
                 srch.Contains("To_b") || srch.Contains("TO_b") )
                last = _EE_SingleElectron_B;
            else if ( srch.Contains("toC") || srch.Contains("ToC") || srch.Contains("TOC") || srch.Contains("to C") || srch.Contains("To C") ||
                      srch.Contains("TO C") || srch.Contains("-C") || srch.Contains("- C") || srch.Contains("to_C") || srch.Contains("To_C") ||
                      srch.Contains("TO_C") || srch.Contains("toc") || srch.Contains("Toc") || srch.Contains("TOc") || srch.Contains("to c") ||
                      srch.Contains("To c") || srch.Contains("TO c") || srch.Contains("-c") || srch.Contains("- c") || srch.Contains("to_c") ||
                      srch.Contains("To_c") || srch.Contains("TO_c") )
                last = _EE_SingleElectron_C;
            else if ( srch.Contains("toD") || srch.Contains("ToD") || srch.Contains("TOD") || srch.Contains("to D") || srch.Contains("To D") ||
                      srch.Contains("TO D") || srch.Contains("-D") || srch.Contains("- D") || srch.Contains("to_D") || srch.Contains("To_D") ||
                      srch.Contains("TO_D") || srch.Contains("tod") || srch.Contains("Tod") || srch.Contains("TOd") || srch.Contains("to d") ||
                      srch.Contains("To d") || srch.Contains("TO d") || srch.Contains("-d") || srch.Contains("- d") || srch.Contains("to_d") ||
                      srch.Contains("To_d") || srch.Contains("TO_d") )
                last = _EE_SingleElectron_D;
            else if ( srch.Contains("toE") || srch.Contains("ToE") || srch.Contains("TOE") || srch.Contains("to E") || srch.Contains("To E") ||
                      srch.Contains("TO E") || srch.Contains("-E") || srch.Contains("- E") || srch.Contains("to_E") || srch.Contains("To_E") ||
                      srch.Contains("TO_E") || srch.Contains("toe") || srch.Contains("Toe") || srch.Contains("TOe") || srch.Contains("to e") ||
                      srch.Contains("To e") || srch.Contains("TO e") || srch.Contains("-e") || srch.Contains("- e") || srch.Contains("to_e") ||
                      srch.Contains("To_e") || srch.Contains("TO_e") )
                last = _EE_SingleElectron_E;
            else if ( srch.Contains("toF") || srch.Contains("ToF") || srch.Contains("TOF") || srch.Contains("to F") || srch.Contains("To F") ||
                      srch.Contains("TO F") || srch.Contains("-F") || srch.Contains("- F") || srch.Contains("to_F") || srch.Contains("To_F") ||
                      srch.Contains("TO_F") || srch.Contains("tof") || srch.Contains("Tof") || srch.Contains("TOf") || srch.Contains("to f") ||
                      srch.Contains("To f") || srch.Contains("TO f") || srch.Contains("-f") || srch.Contains("- f") || srch.Contains("to_f") ||
                      srch.Contains("To_f") || srch.Contains("TO_f") )
                last = _EE_SingleElectron_F;
            else if ( srch.Contains("toG") || srch.Contains("ToG") || srch.Contains("TOG") || srch.Contains("to G") || srch.Contains("To G") ||
                      srch.Contains("TO G") || srch.Contains("-G") || srch.Contains("- G") || srch.Contains("to_G") || srch.Contains("To_G") ||
                      srch.Contains("TO_G") || srch.Contains("tog") || srch.Contains("Tog") || srch.Contains("TOg") || srch.Contains("to g") ||
                      srch.Contains("To g") || srch.Contains("TO g") || srch.Contains("-g") || srch.Contains("- g") || srch.Contains("to_g") ||
                      srch.Contains("To_g") || srch.Contains("TO_g") )
                last = _EE_SingleElectron_G;
            else if ( srch.Contains("toH") || srch.Contains("ToH") || srch.Contains("TOH") || srch.Contains("to H") || srch.Contains("To H") ||
                      srch.Contains("TO H") || srch.Contains("-H") || srch.Contains("- H") || srch.Contains("to_H") || srch.Contains("To_H") ||
                      srch.Contains("TO_H") || srch.Contains("toh") || srch.Contains("Toh") || srch.Contains("TOh") || srch.Contains("to h") ||
                      srch.Contains("To h") || srch.Contains("TO h") || srch.Contains("-h") || srch.Contains("- h") || srch.Contains("to_h") ||
                      srch.Contains("To_h") || srch.Contains("TO_h") )
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
            else if ( srch.Contains("Full") || srch.Contains("full") || srch.Contains("FULL") )
            {
                Result.push_back(_EE_SingleElectron_Full);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_Full] << "." << endl;
            }
            else if ( srch.Contains("_B") || srch.Contains("_b") || srch.Contains("runB") || srch.Contains("RUNB") || srch.Contains("RunB") || srch.Contains("runb") )
            {
                Result.push_back(_EE_SingleElectron_B);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_B] << "." << endl;
            }
            else if ( srch.Contains("_C") || srch.Contains("_c") || srch.Contains("runC") || srch.Contains("RUNC") || srch.Contains("RunC") || srch.Contains("runc") )
            {
                Result.push_back(_EE_SingleElectron_C);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_C] << "." << endl;
            }
            else if ( srch.Contains("_D") || srch.Contains("_d") || srch.Contains("runD") || srch.Contains("RUND") || srch.Contains("RunD") || srch.Contains("runD") )
            {
                Result.push_back(_EE_SingleElectron_D);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_D] << "." << endl;
            }
            else if ( srch.Contains("_E") || srch.Contains("_e") || srch.Contains("runE") || srch.Contains("RUNE") || srch.Contains("RunE") || srch.Contains("runE") )
            {
                Result.push_back(_EE_SingleElectron_E);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_E] << "." << endl;
            }
            else if ( srch.Contains("_F") || srch.Contains("_f") || srch.Contains("runF") || srch.Contains("RUNF") || srch.Contains("RunF") || srch.Contains("runF") )
            {
                Result.push_back(_EE_SingleElectron_F);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_F] << "." << endl;
            }
            else if ( srch.Contains("_G") || srch.Contains("_g") || srch.Contains("runG") || srch.Contains("RUNG") || srch.Contains("RunG") || srch.Contains("runG") )
            {
                Result.push_back(_EE_SingleElectron_G);
                if ( notify == kTRUE ) cout << Procname[_EE_SingleElectron_G] << "." << endl;
            }
            else if ( srch.Contains("_H") || srch.Contains("_h") || srch.Contains("runH") || srch.Contains("RUNH") || srch.Contains("RunH") || srch.Contains("runH") )
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

    else if ( srch.Contains("Signal") || srch.Contains("signal") || srch.Contains("SIGNAL") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_DY_Full);
            if ( notify == kTRUE ) cout << Procname[_MuMu_DY_Full] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_DY_Full);
            if ( notify == kTRUE ) cout << Procname[_EE_DY_Full] << "." << endl;
        }
    }

    else if ( srch.Contains("Bkg") || srch.Contains("bkg") || srch.Contains("BKG") || srch.Contains("Background") || srch.Contains("background") ||
              srch.Contains("BACKGROUND") )
    {
        if ( srch.Contains("MuMu") || srch.Contains("mumu") || srch.Contains("MUMU") || srch.Contains("Dimuon") || srch.Contains("diMuon") ||
             srch.Contains("DiMuon") || srch.Contains("dimuon") || srch.Contains("DIMUON") )
        {
            Result.push_back(_MuMu_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_MuMu_Bkg_Full] << "." << endl;
        }
        else if ( srch.Contains("EE") || srch.Contains("ee") || srch.Contains("Ee") || srch.Contains("eE") || srch.Contains("dielectron") ||
                  srch.Contains("Dielectron") || srch.Contains("diElectron") || srch.Contains("DiElectron") || srch.Contains("DIELECTRON") )
        {
            Result.push_back(_EE_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_EE_Bkg_Full] << "." << endl;
        }
        else if ( srch.Contains("EMu") || srch.Contains("emu") || srch.Contains("EMU") || srch.Contains("eMu") || srch.Contains("eMU") ||
                  srch.Contains("Emu") || ( srch.Contains("Ele") && srch.Contains("Mu") ) )
        {
            Result.push_back(_EMu_Bkg_Full);
            if ( notify == kTRUE )  cout << Procname[_EMu_Bkg_Full] << "." << endl;
        }
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
    Procname[_EndOf_MuMu_WJets] = "EndOf_SelectedMuMu_WJets";
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
    Procname[_EndOf_EE_WJets] = "EndOf_SelectedEE_WJets";
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
    Procname[_EndOf_EMu_WJets] = "EndOf_SelectedEMu_WJets";
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
