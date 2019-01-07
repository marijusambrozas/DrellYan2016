#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TLine.h>
#include <TColor.h>
#include <TF1.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

void EE_HistDrawer (  TString whichGraphs, TString type );
void MuMu_HistDrawer ( TString whichGraphs, TString type );
void EMu_HistDrawer ( TString whichGraphs, TString type );

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


void HistDrawer ( TString WhichX = "", TString WhichGraphs = "ALL", TString type = "" )
{
    TString whichX = WhichX;
    whichX.ToUpper();
    TString whichGraphs = WhichGraphs;
    whichGraphs.ToUpper();
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") )
    {
        Xselected++;
        cout << "\n*******      EE_HistDrawer ( " << whichGraphs << " " << type << " )      *******" << endl;
        EE_HistDrawer( whichGraphs, type );
    }
    if ( whichX.Contains("MUMU") )
    {
        Xselected++;
        cout << "\n*****  MuMu_HistDrawer ( " << whichGraphs << " " << type << " )  *****" << endl;
        MuMu_HistDrawer( whichGraphs, type );
    }
    if ( whichX.Contains("EMU") )
    {
        Xselected++;
        cout << "\n*****   EMu_HistDrawer ( " << whichGraphs << " " << type << " )  *****" << endl;
        EMu_HistDrawer( whichGraphs, type );
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ----------------------------- Electron Channel ------------------------------ ///
void EE_HistDrawer ( TString whichGraphs, TString type )
{
    if ( !whichGraphs.Length() )
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc( _EE_DY_Full );
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _EE_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _EE_DoubleEG_Full );
//    Mgr.SetProc( _EE_SingleElectron_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=1.3107e+07, MCintegral;

//################################# INVARIANT MASS #################################################

    if( whichGraphs=="ALL" || whichGraphs=="INVMASS" )
    {
        count_drawn++;

        THStack *s_mass_fine_before_PUCorr, *s_mass_fine_before_EffCorr, *s_mass_fine, *s_mass;
        s_mass_fine_before_PUCorr = new THStack("s_mass_fine_before_PUCorr", "");
        s_mass_fine_before_EffCorr = new THStack("s_mass_fine_before_EffCorr", "");
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine_before_PUCorr[9], *h_bkg_mass_fine_before_EffCorr[9],
             *h_bkg_mass_fine[9], *h_bkg_mass[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            Color_t color = kBlack;
            if ( pr == _EE_QCDEMEnriched_Full ) color = kRed + 3;
            if ( pr == _EE_WJets ) color = kRed - 2;
            if ( pr == _EE_WW ) color = kMagenta - 5;
            if ( pr == _EE_WZ ) color = kMagenta - 2;
            if ( pr == _EE_ZZ ) color = kMagenta - 6;
            if ( pr == _EE_tbarW ) color = kGreen - 2;
            if ( pr == _EE_tW ) color = kGreen + 2;
            if ( pr == _EE_ttbar_Full ) color = kCyan + 2;
            if ( pr == _EE_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_fine[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetFillColor(color);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_fine[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetLineColor(color);

            h_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
            s_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );

            iter++;

            if ( pr == _EE_QCDEMEnriched_Full ) // next - WJets
                pr = _EndOf_EE_WJets;
            if ( pr == _EE_WJets ) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if ( pr == _EE_tW ) // next -- ttbar
                pr = _EE_VVnST;
            if ( pr == _EE_DYTauTau_Full ) // last
                break;

        } // End of for(bkg)       

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_fine_before_PUCorr, *h_DY_mass_fine_before_EffCorr, *h_DY_mass_fine, *h_DY_mass;

        f_DY->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine_before_PUCorr );
        f_DY->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine_before_EffCorr );
        f_DY->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine );
        f_DY->GetObject( "h_mass_"+Mgr.Procname[_EE_DY_Full], h_DY_mass );

        h_DY_mass_fine_before_PUCorr->SetFillColor(kOrange);
        h_DY_mass_fine_before_EffCorr->SetFillColor(kOrange);
        h_DY_mass_fine->SetFillColor(kOrange);
        h_DY_mass->SetFillColor(kOrange);

        h_DY_mass_fine_before_PUCorr->SetLineColor(kOrange);
        h_DY_mass_fine_before_EffCorr->SetLineColor(kOrange);
        h_DY_mass_fine->SetLineColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);

        h_DY_mass_fine_before_PUCorr->SetDirectory(0);
        h_DY_mass_fine_before_EffCorr->SetDirectory(0);
        h_DY_mass_fine->SetDirectory(0);
        h_DY_mass->SetDirectory(0);

        s_mass_fine_before_PUCorr->Add( h_DY_mass_fine_before_PUCorr );
        s_mass_fine_before_EffCorr->Add( h_DY_mass_fine_before_EffCorr );
        s_mass_fine->Add( h_DY_mass_fine );
        s_mass->Add( h_DY_mass );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine_before_PUCorr, *h_data_mass_fine_before_EffCorr, *h_data_mass_fine, *h_data_mass;

        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass );

        h_data_mass_fine_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine_before_PUCorr->SetMarkerColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetMarkerColor(kBlack);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetLineColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetLineColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetDirectory(0);
        h_data_mass_fine_before_EffCorr->SetDirectory(0);
        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_fine_before_PUCorr, *RP_mass_fine_before_EffCorr, *RP_mass_fine, *RP_mass;
        RP_mass_fine_before_PUCorr = new myRatioPlot_t( "RP_mass_fine_before_PUCorr", s_mass_fine_before_PUCorr, h_data_mass_fine_before_PUCorr );
        RP_mass_fine_before_EffCorr = new myRatioPlot_t( "RP_mass_fine_before_EffCorr", s_mass_fine_before_EffCorr, h_data_mass_fine_before_EffCorr );
        RP_mass_fine = new myRatioPlot_t( "RP_mass_fine", s_mass_fine, h_data_mass_fine );
        RP_mass = new myRatioPlot_t( "RP_mass", s_mass, h_data_mass );

        RP_mass_fine_before_PUCorr->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_mass_fine_before_EffCorr->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before Efficiency SF)", 15, 3000);
        RP_mass_fine->SetPlots("Dielectron invariant mass [GeV/c^{2}]", 15, 3000);
        RP_mass->SetPlots("Dielectron invariant mass [GeV/c^{2}]", 15, 3000);

        RP_mass_fine_before_PUCorr->SetLegend(0.75, 0.4);
        RP_mass_fine_before_EffCorr->SetLegend(0.75, 0.4);
        RP_mass_fine->SetLegend(0.75, 0.4);
        RP_mass->SetLegend(0.75, 0.4);

        // Legend data
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_data_mass_fine_before_PUCorr, "Data", "lp");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_data_mass_fine_before_EffCorr, "Data", "lp");
        RP_mass_fine->AddLegendEntry(h_data_mass_fine, "Data", "lp");
        RP_mass->AddLegendEntry(h_data_mass, "Data", "lp");

        // Legend MC signal
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_DY_mass_fine_before_PUCorr, "DY#rightarrow ee", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_DY_mass_fine_before_EffCorr, "DY#rightarrow ee", "f");
        RP_mass_fine->AddLegendEntry(h_DY_mass_fine, "DY#rightarrow ee", "f");
        RP_mass->AddLegendEntry(h_DY_mass, "DY#rightarrow ee", "f");

        // Legend MC BKG
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[8], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[7], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[7], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[7], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[7], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[6], "tW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[6], "tW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[6], "tW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[6], "tW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[5], "#bar{t}W", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[5], "#bar{t}W", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[5], "#bar{t}W", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[5], "#bar{t}W", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4], "ZZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4], "ZZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4], "ZZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4], "ZZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "WZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "WZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "WZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "WZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2], "WW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2], "WW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2], "WW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2], "WW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "W+Jets", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "W+Jets", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "W+Jets", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "W+Jets", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "QCD", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "QCD", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "QCD", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "QCD", "f");

        RP_mass_fine_before_PUCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_EffCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine->Draw(0.5, 3e6, 1);
        RP_mass->Draw(0.5, 3e6, 1);

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################

    if ( whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI") )
    {
        count_drawn++;

        THStack *s_pT, *s_rapi;
        s_rapi = new THStack( "s_rapi", "" );
        s_pT = new THStack( "s_pT", "" );

        THStack *s_pT_lead_before_PUCorr = new THStack( "s_pT_lead_before_PUCorr", "" );
        THStack *s_pT_sublead_before_PUCorr = new THStack( "s_pT_sublead_before_PUCorr", "" );
        THStack *s_eta_lead_before_PUCorr = new THStack( "s_eta_lead_before_PUCorr", "" );
        THStack *s_eta_sublead_before_PUCorr = new THStack( "s_eta_sublead_before_PUCorr", "" );
        THStack *s_phi_lead_before_PUCorr = new THStack( "s_phi_lead_before_PUCorr", "" );
        THStack *s_phi_sublead_before_PUCorr = new THStack( "s_phi_sublead_before_PUCorr", "" );

        THStack *s_pT_lead_before_EffCorr = new THStack( "s_pT_lead_before_EffCorr", "" );
        THStack *s_pT_sublead_before_EffCorr = new THStack( "s_pT_sublead_before_EffCorr", "" );
        THStack *s_eta_lead_before_EffCorr = new THStack( "s_eta_lead_before_EffCorr", "" );
        THStack *s_eta_sublead_before_EffCorr = new THStack( "s_eta_sublead_before_EffCorr", "" );
        THStack *s_phi_lead_before_EffCorr = new THStack( "s_phi_lead_before_EffCorr", "" );
        THStack *s_phi_sublead_before_EffCorr = new THStack( "s_phi_sublead_before_EffCorr", "" );

        THStack *s_pT_lead = new THStack( "s_pT_lead", "" );
        THStack *s_pT_sublead = new THStack( "s_pT_sublead", "" );
        THStack *s_eta_lead = new THStack( "s_eta_lead", "" );
        THStack *s_eta_sublead = new THStack( "s_eta_sublead", "" );
        THStack *s_phi_lead = new THStack( "s_phi_lead", "" );
        THStack *s_phi_sublead = new THStack( "s_phi_sublead", "" );


//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_pT[9], *h_bkg_rapi[9], *h_bkg_pT_lead_before_PUCorr[9], *h_bkg_pT_lead_before_EffCorr[9], *h_bkg_pT_lead[9],
             *h_bkg_pT_sublead_before_PUCorr[9], *h_bkg_pT_sublead_before_EffCorr[9], *h_bkg_pT_sublead[9],
             *h_bkg_eta_lead_before_PUCorr[9], *h_bkg_eta_lead_before_EffCorr[9], *h_bkg_eta_lead[9],
             *h_bkg_eta_sublead_before_PUCorr[9], *h_bkg_eta_sublead_before_EffCorr[9], *h_bkg_eta_sublead[9],
             *h_bkg_phi_lead_before_PUCorr[9], *h_bkg_phi_lead_before_EffCorr[9], *h_bkg_phi_lead[9],
             *h_bkg_phi_sublead_before_PUCorr[9], *h_bkg_phi_sublead_before_EffCorr[9], *h_bkg_phi_sublead[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg->GetObject( "h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter] );
            f_bkg->GetObject( "h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter] );
            f_bkg->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_pT_lead_"+Mgr.Procname[pr], h_bkg_pT_lead[iter] );
            f_bkg->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_pT_sublead_"+Mgr.Procname[pr], h_bkg_pT_sublead[iter] );
            f_bkg->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_eta_lead_"+Mgr.Procname[pr], h_bkg_eta_lead[iter] );
            f_bkg->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_eta_sublead_"+Mgr.Procname[pr], h_bkg_eta_sublead[iter] );
            f_bkg->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_phi_lead_"+Mgr.Procname[pr], h_bkg_phi_lead[iter] );
            f_bkg->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PUCorr[iter] );
            f_bkg->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_phi_sublead_"+Mgr.Procname[pr], h_bkg_phi_sublead[iter] );

            Color_t color = kBlack;
            if ( pr == _EE_QCDEMEnriched_Full ) color = kRed + 3;
            if ( pr == _EE_WJets ) color = kRed - 2;
            if ( pr == _EE_WW ) color = kMagenta - 5;
            if ( pr == _EE_WZ ) color = kMagenta - 2;
            if ( pr == _EE_ZZ ) color = kMagenta - 6;
            if ( pr == _EE_tbarW ) color = kGreen - 2;
            if ( pr == _EE_tW ) color = kGreen + 2;
            if ( pr == _EE_ttbar_Full ) color = kCyan + 2;
            if ( pr == _EE_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_pT[iter]->SetFillColor(color);
            h_bkg_rapi[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead[iter]->SetFillColor(color);

            h_bkg_pT[iter]->SetLineColor(color);
            h_bkg_rapi[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead[iter]->SetLineColor(color);


            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead[iter]->SetDirectory(0);


            s_pT->Add( h_bkg_pT[iter] );
            s_rapi->Add( h_bkg_rapi[iter] );
            s_pT_lead_before_PUCorr->Add( h_bkg_pT_lead_before_PUCorr[iter] );
            s_pT_lead_before_EffCorr->Add( h_bkg_pT_lead_before_EffCorr[iter] );
            s_pT_lead->Add( h_bkg_pT_lead[iter] );
            s_pT_sublead_before_PUCorr->Add( h_bkg_pT_sublead_before_PUCorr[iter] );
            s_pT_sublead_before_EffCorr->Add( h_bkg_pT_sublead_before_EffCorr[iter] );
            s_pT_sublead->Add( h_bkg_pT_sublead[iter] );
            s_eta_lead_before_PUCorr->Add( h_bkg_eta_lead_before_PUCorr[iter] );
            s_eta_lead_before_EffCorr->Add( h_bkg_eta_lead_before_EffCorr[iter] );
            s_eta_lead->Add( h_bkg_eta_lead[iter] );
            s_eta_sublead_before_PUCorr->Add( h_bkg_eta_sublead_before_PUCorr[iter] );
            s_eta_sublead_before_EffCorr->Add( h_bkg_eta_sublead_before_EffCorr[iter] );
            s_eta_sublead->Add( h_bkg_eta_sublead[iter] );
            s_phi_lead_before_PUCorr->Add( h_bkg_phi_lead_before_PUCorr[iter] );
            s_phi_lead_before_EffCorr->Add( h_bkg_phi_lead_before_EffCorr[iter] );
            s_phi_lead->Add( h_bkg_phi_lead[iter] );
            s_phi_sublead_before_PUCorr->Add( h_bkg_phi_sublead_before_PUCorr[iter] );
            s_phi_sublead_before_EffCorr->Add( h_bkg_phi_sublead_before_EffCorr[iter] );
            s_phi_sublead->Add( h_bkg_phi_sublead[iter] );

            iter++;

            if ( pr == _EE_QCDEMEnriched_Full ) // next - WJets
                pr = _EndOf_EE_WJets;
            if ( pr == _EE_WJets ) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if ( pr == _EE_tW ) // next -- ttbar
                pr = _EE_VVnST;
            if ( pr == _EE_DYTauTau_Full ) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_pT, *h_DY_rapi, *h_DY_pT_lead_before_PUCorr, *h_DY_pT_lead_before_EffCorr, *h_DY_pT_lead,
             *h_DY_pT_sublead_before_PUCorr, *h_DY_pT_sublead_before_EffCorr, *h_DY_pT_sublead,
             *h_DY_eta_lead_before_PUCorr, *h_DY_eta_lead_before_EffCorr, *h_DY_eta_lead,
             *h_DY_eta_sublead_before_PUCorr, *h_DY_eta_sublead_before_EffCorr, *h_DY_eta_sublead,
             *h_DY_phi_lead_before_PUCorr, *h_DY_phi_lead_before_EffCorr, *h_DY_phi_lead,
             *h_DY_phi_sublead_before_PUCorr, *h_DY_phi_sublead_before_EffCorr, *h_DY_phi_sublead;

        f_DY->GetObject( "h_pT_"+Mgr.Procname[_EE_DY_Full], h_DY_pT );
        f_DY->GetObject( "h_rapi_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi );
        f_DY->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_PUCorr );
        f_DY->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_EffCorr );
        f_DY->GetObject( "h_pT_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead );
        f_DY->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_PUCorr );
        f_DY->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_EffCorr );
        f_DY->GetObject( "h_pT_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead );
        f_DY->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_PUCorr );
        f_DY->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_EffCorr );
        f_DY->GetObject( "h_eta_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead );
        f_DY->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_PUCorr );
        f_DY->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_EffCorr );
        f_DY->GetObject( "h_eta_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead );
        f_DY->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_PUCorr );
        f_DY->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_EffCorr );
        f_DY->GetObject( "h_phi_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead );
        f_DY->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_PUCorr );
        f_DY->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_EffCorr );
        f_DY->GetObject( "h_phi_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead );


        h_DY_pT->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_lead->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_sublead->SetFillColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_lead->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_sublead->SetFillColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_lead->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_sublead->SetFillColor(kOrange);


        h_DY_pT->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_lead->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_sublead->SetLineColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_lead->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_sublead->SetLineColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_lead->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_sublead->SetLineColor(kOrange);


        h_DY_pT->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT_lead_before_PUCorr->SetDirectory(0);
        h_DY_pT_lead_before_EffCorr->SetDirectory(0);
        h_DY_pT_lead->SetDirectory(0);
        h_DY_pT_sublead_before_PUCorr->SetDirectory(0);
        h_DY_pT_sublead_before_EffCorr->SetDirectory(0);
        h_DY_pT_sublead->SetDirectory(0);
        h_DY_eta_lead_before_PUCorr->SetDirectory(0);
        h_DY_eta_lead_before_EffCorr->SetDirectory(0);
        h_DY_eta_lead->SetDirectory(0);
        h_DY_eta_sublead_before_PUCorr->SetDirectory(0);
        h_DY_eta_sublead_before_EffCorr->SetDirectory(0);
        h_DY_eta_sublead->SetDirectory(0);
        h_DY_phi_lead_before_PUCorr->SetDirectory(0);
        h_DY_phi_lead_before_EffCorr->SetDirectory(0);
        h_DY_phi_lead->SetDirectory(0);
        h_DY_phi_sublead_before_PUCorr->SetDirectory(0);
        h_DY_phi_sublead_before_EffCorr->SetDirectory(0);
        h_DY_phi_sublead->SetDirectory(0);


        s_pT->Add( h_DY_pT );
        s_rapi->Add( h_DY_rapi );
        s_pT_lead_before_PUCorr->Add( h_DY_pT_lead_before_PUCorr );
        s_pT_lead_before_EffCorr->Add( h_DY_pT_lead_before_EffCorr );
        s_pT_lead->Add( h_DY_pT_lead );
        s_pT_sublead_before_PUCorr->Add( h_DY_pT_sublead_before_PUCorr );
        s_pT_sublead_before_EffCorr->Add( h_DY_pT_sublead_before_EffCorr );
        s_pT_sublead->Add( h_DY_pT_sublead );
        s_eta_lead_before_PUCorr->Add( h_DY_eta_lead_before_PUCorr );
        s_eta_lead_before_EffCorr->Add( h_DY_eta_lead_before_EffCorr );
        s_eta_lead->Add( h_DY_eta_lead );
        s_eta_sublead_before_PUCorr->Add( h_DY_eta_sublead_before_PUCorr );
        s_eta_sublead_before_EffCorr->Add( h_DY_eta_sublead_before_EffCorr );
        s_eta_sublead->Add( h_DY_eta_sublead );
        s_phi_lead_before_PUCorr->Add( h_DY_phi_lead_before_PUCorr );
        s_phi_lead_before_EffCorr->Add( h_DY_phi_lead_before_EffCorr );
        s_phi_lead->Add( h_DY_phi_lead );
        s_phi_sublead_before_PUCorr->Add( h_DY_phi_sublead_before_PUCorr );
        s_phi_sublead_before_EffCorr->Add( h_DY_phi_sublead_before_EffCorr );
        s_phi_sublead->Add( h_DY_phi_sublead );


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT, *h_data_rapi, *h_data_pT_lead_before_PUCorr, *h_data_pT_lead_before_EffCorr, *h_data_pT_lead,
             *h_data_pT_sublead_before_PUCorr, *h_data_pT_sublead_before_EffCorr, *h_data_pT_sublead,
             *h_data_eta_lead_before_PUCorr, *h_data_eta_lead_before_EffCorr, *h_data_eta_lead,
             *h_data_eta_sublead_before_PUCorr, *h_data_eta_sublead_before_EffCorr, *h_data_eta_sublead,
             *h_data_phi_lead_before_PUCorr, *h_data_phi_lead_before_EffCorr, *h_data_phi_lead,
             *h_data_phi_sublead_before_PUCorr, *h_data_phi_sublead_before_EffCorr, *h_data_phi_sublead;

        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject( "h_pT_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT );
        f_data->GetObject( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], h_data_rapi );
        f_data->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_lead_before_PUCorr );
        f_data->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_lead_before_EffCorr );
        f_data->GetObject( "h_pT_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_lead );
        f_data->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_sublead_before_PUCorr );
        f_data->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_sublead_before_EffCorr );
        f_data->GetObject( "h_pT_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_sublead );
        f_data->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_lead_before_PUCorr );
        f_data->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_lead_before_EffCorr );
        f_data->GetObject( "h_eta_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_lead );
        f_data->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_sublead_before_PUCorr );
        f_data->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_sublead_before_EffCorr );
        f_data->GetObject( "h_eta_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_sublead );
        f_data->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_lead_before_PUCorr );
        f_data->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_lead_before_EffCorr );
        f_data->GetObject( "h_phi_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_lead );
        f_data->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_sublead_before_PUCorr );
        f_data->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_sublead_before_EffCorr );
        f_data->GetObject( "h_phi_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_sublead );


        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead->SetMarkerStyle(kFullDotLarge);


        h_data_pT->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);
        h_data_pT_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_pT_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_pT_lead->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_pT_sublead->SetMarkerColor(kBlack);
        h_data_eta_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_eta_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_eta_lead->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_eta_sublead->SetMarkerColor(kBlack);
        h_data_phi_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_phi_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_phi_lead->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_phi_sublead->SetMarkerColor(kBlack);


        h_data_pT->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);
        h_data_pT_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_pT_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_pT_lead->SetLineColor(kBlack);
        h_data_pT_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_pT_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_pT_sublead->SetLineColor(kBlack);
        h_data_eta_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_eta_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_eta_lead->SetLineColor(kBlack);
        h_data_eta_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_eta_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_eta_sublead->SetLineColor(kBlack);
        h_data_phi_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_phi_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_phi_lead->SetLineColor(kBlack);
        h_data_phi_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_phi_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_phi_sublead->SetLineColor(kBlack);


        h_data_pT->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT_lead_before_PUCorr->SetDirectory(0);
        h_data_pT_lead_before_EffCorr->SetDirectory(0);
        h_data_pT_lead->SetDirectory(0);
        h_data_pT_sublead_before_PUCorr->SetDirectory(0);
        h_data_pT_sublead_before_EffCorr->SetDirectory(0);
        h_data_pT_sublead->SetDirectory(0);
        h_data_eta_lead_before_PUCorr->SetDirectory(0);
        h_data_eta_lead_before_EffCorr->SetDirectory(0);
        h_data_eta_lead->SetDirectory(0);
        h_data_eta_sublead_before_PUCorr->SetDirectory(0);
        h_data_eta_sublead_before_EffCorr->SetDirectory(0);
        h_data_eta_sublead->SetDirectory(0);
        h_data_phi_lead_before_PUCorr->SetDirectory(0);
        h_data_phi_lead_before_EffCorr->SetDirectory(0);
        h_data_phi_lead->SetDirectory(0);
        h_data_phi_sublead_before_PUCorr->SetDirectory(0);
        h_data_phi_sublead_before_EffCorr->SetDirectory(0);
        h_data_phi_sublead->SetDirectory(0);


//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_pT, *RP_rapi, *RP_pT_lead_before_PUCorr, *RP_pT_lead_before_EffCorr, *RP_pT_lead,
                      *RP_pT_sublead_before_PUCorr, *RP_pT_sublead_before_EffCorr, *RP_pT_sublead,
                      *RP_eta_lead_before_PUCorr, *RP_eta_lead_before_EffCorr, *RP_eta_lead,
                      *RP_eta_sublead_before_PUCorr, *RP_eta_sublead_before_EffCorr, *RP_eta_sublead,
                      *RP_phi_lead_before_PUCorr, *RP_phi_lead_before_EffCorr, *RP_phi_lead,
                      *RP_phi_sublead_before_PUCorr, *RP_phi_sublead_before_EffCorr, *RP_phi_sublead;
        RP_pT = new myRatioPlot_t( "RP_pT", s_pT, h_data_pT );
        RP_rapi = new myRatioPlot_t( "RP_rapi", s_rapi, h_data_rapi );
        RP_pT_lead_before_PUCorr = new myRatioPlot_t( "RP_pT_lead_before_PUCorr", s_pT_lead_before_PUCorr, h_data_pT_lead_before_PUCorr );
        RP_pT_lead_before_EffCorr = new myRatioPlot_t( "RP_pT_lead_before_EffCorr", s_pT_lead_before_EffCorr, h_data_pT_lead_before_EffCorr );
        RP_pT_lead = new myRatioPlot_t( "RP_pT_lead", s_pT_lead, h_data_pT_lead );
        RP_pT_sublead_before_PUCorr = new myRatioPlot_t( "RP_pT_sublead_before_PUCorr", s_pT_sublead_before_PUCorr, h_data_pT_sublead_before_PUCorr );
        RP_pT_sublead_before_EffCorr = new myRatioPlot_t( "RP_pT_sublead_before_EffCorr", s_pT_sublead_before_EffCorr, h_data_pT_sublead_before_EffCorr );
        RP_pT_sublead = new myRatioPlot_t( "RP_pT_sublead", s_pT_sublead, h_data_pT_sublead );
        RP_eta_lead_before_PUCorr = new myRatioPlot_t( "RP_eta_lead_before_PUCorr", s_eta_lead_before_PUCorr, h_data_eta_lead_before_PUCorr );
        RP_eta_lead_before_EffCorr = new myRatioPlot_t( "RP_eta_lead_before_EffCorr", s_eta_lead_before_EffCorr, h_data_eta_lead_before_EffCorr );
        RP_eta_lead = new myRatioPlot_t( "RP_eta_lead", s_eta_lead, h_data_eta_lead );
        RP_eta_sublead_before_PUCorr = new myRatioPlot_t( "RP_eta_sublead_before_PUCorr", s_eta_sublead_before_PUCorr, h_data_eta_sublead_before_PUCorr );
        RP_eta_sublead_before_EffCorr = new myRatioPlot_t( "RP_eta_sublead_before_EffCorr", s_eta_sublead_before_EffCorr, h_data_eta_sublead_before_EffCorr );
        RP_eta_sublead = new myRatioPlot_t( "RP_eta_sublead", s_eta_sublead, h_data_eta_sublead );
        RP_phi_lead_before_PUCorr = new myRatioPlot_t( "RP_phi_lead_before_PUCorr", s_phi_lead_before_PUCorr, h_data_phi_lead_before_PUCorr );
        RP_phi_lead_before_EffCorr = new myRatioPlot_t( "RP_phi_lead_before_EffCorr", s_phi_lead_before_EffCorr, h_data_phi_lead_before_EffCorr );
        RP_phi_lead = new myRatioPlot_t( "RP_phi_lead", s_phi_lead, h_data_phi_lead );
        RP_phi_sublead_before_PUCorr = new myRatioPlot_t( "RP_phi_sublead_before_PUCorr", s_phi_sublead_before_PUCorr, h_data_phi_sublead_before_PUCorr );
        RP_phi_sublead_before_EffCorr = new myRatioPlot_t( "RP_phi_sublead_before_EffCorr", s_phi_sublead_before_EffCorr, h_data_phi_sublead_before_EffCorr );
        RP_phi_sublead = new myRatioPlot_t( "RP_phi_sublead", s_phi_sublead, h_data_phi_sublead );


        RP_pT->SetPlots("Dielectron p_{T} [GeV/c]", 0, 1000);
        RP_rapi->SetPlots("Dielectron rapidity", -4, 4);
        RP_pT_lead_before_PUCorr->SetPlots("Leading electron p_{T} [GeV/c] (before PU correction)", 0, 1000);
        RP_pT_lead_before_EffCorr->SetPlots("Leading electron p_{T} [GeV/c] (before Efficiency SF)", 0, 1000);
        RP_pT_lead->SetPlots("Leading electron p_{T} [GeV/c]", 0, 1000);
        RP_pT_sublead_before_PUCorr->SetPlots("Subleading electron p_{T} [GeV/c] (before PU correction)", 0, 1000);
        RP_pT_sublead_before_EffCorr->SetPlots("Subleading electron p_{T} [GeV/c] (before Efficiency SF)", 0, 1000);
        RP_pT_sublead->SetPlots("Subleading electron p_{T} [GeV/c]", 0, 1000);
        RP_eta_lead_before_PUCorr->SetPlots("Leading electron #eta_{SC} (before PU correction)", -3.5, 3.5);
        RP_eta_lead_before_EffCorr->SetPlots("Leading electron #eta_{SC} (before Efficiency SF)", -3.5, 3.5);
        RP_eta_lead->SetPlots("Leading electron #eta_{SC}", -3.5, 3.5);
        RP_eta_sublead_before_PUCorr->SetPlots("Subleading electron #eta_{SC} (before PU correction)", -3.5, 3.5);
        RP_eta_sublead_before_EffCorr->SetPlots("Subleading electron #eta_{SC} (before Efficiency SF)", -3.5, 3.5);
        RP_eta_sublead->SetPlots("Subleading electron #eta_{SC}", -3.5, 3.5);
        RP_phi_lead_before_PUCorr->SetPlots("Leading electron #phi (before PU correction)", -4, 4);
        RP_phi_lead_before_EffCorr->SetPlots("Leading electron #phi (before Efficiency SF)", -4, 4);
        RP_phi_lead->SetPlots("Leading electron #phi", -4, 4);
        RP_phi_sublead_before_PUCorr->SetPlots("Subleading electron #phi (before PU correction)", -4, 4);
        RP_phi_sublead_before_EffCorr->SetPlots("Subleading electron #phi (before Efficiency SF)", -4, 4);
        RP_phi_sublead->SetPlots("Subleading electron #phi", -4, 4);


        RP_pT->SetLegend(0.75, 0.4);
        RP_rapi->SetLegend(0.75, 0.4);
        RP_pT_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_pT_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_pT_lead->SetLegend(0.75, 0.4);
        RP_pT_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_pT_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_pT_sublead->SetLegend(0.75, 0.4);
        RP_eta_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_eta_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_eta_lead->SetLegend(0.75, 0.4);
        RP_eta_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_eta_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_eta_sublead->SetLegend(0.75, 0.4);
        RP_phi_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_phi_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_phi_lead->SetLegend(0.75, 0.4);
        RP_phi_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_phi_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_phi_sublead->SetLegend(0.75, 0.4);


        // Legend data
        RP_pT->AddLegendEntry(h_data_pT, "Data", "lp");
        RP_rapi->AddLegendEntry(h_data_rapi, "Data", "lp");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_data_pT_lead_before_PUCorr, "Data", "lp");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_data_pT_lead_before_EffCorr, "Data", "lp");
        RP_pT_lead->AddLegendEntry(h_data_pT_lead, "Data", "lp");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_data_pT_sublead_before_PUCorr, "Data", "lp");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_data_pT_sublead_before_EffCorr, "Data", "lp");
        RP_pT_sublead->AddLegendEntry(h_data_pT_sublead, "Data", "lp");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_data_eta_lead_before_PUCorr, "Data", "lp");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_data_eta_lead_before_EffCorr, "Data", "lp");
        RP_eta_lead->AddLegendEntry(h_data_eta_lead, "Data", "lp");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_data_eta_sublead_before_PUCorr, "Data", "lp");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_data_eta_sublead_before_EffCorr, "Data", "lp");
        RP_eta_sublead->AddLegendEntry(h_data_eta_sublead, "Data", "lp");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_data_phi_lead_before_PUCorr, "Data", "lp");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_data_phi_lead_before_EffCorr, "Data", "lp");
        RP_phi_lead->AddLegendEntry(h_data_phi_lead, "Data", "lp");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_data_phi_sublead_before_PUCorr, "Data", "lp");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_data_phi_sublead_before_EffCorr, "Data", "lp");
        RP_phi_sublead->AddLegendEntry(h_data_phi_sublead, "Data", "lp");


        // Legend MC signal
        RP_pT->AddLegendEntry(h_DY_pT, "DY#rightarrow ee", "f");
        RP_rapi->AddLegendEntry(h_DY_rapi, "DY#rightarrow ee", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_DY_pT_lead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_DY_pT_lead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_pT_lead->AddLegendEntry(h_DY_pT_lead, "DY#rightarrow ee", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_DY_pT_sublead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_DY_pT_sublead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_pT_sublead->AddLegendEntry(h_DY_pT_sublead, "DY#rightarrow ee", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_DY_eta_lead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_DY_eta_lead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_eta_lead->AddLegendEntry(h_DY_eta_lead, "DY#rightarrow ee", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_DY_eta_sublead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_DY_eta_sublead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_eta_sublead->AddLegendEntry(h_DY_eta_sublead, "DY#rightarrow ee", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_DY_phi_lead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_DY_phi_lead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_phi_lead->AddLegendEntry(h_DY_phi_lead, "DY#rightarrow ee", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_DY_phi_sublead_before_PUCorr, "DY#rightarrow ee", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_DY_phi_sublead_before_EffCorr, "DY#rightarrow ee", "f");
        RP_phi_sublead->AddLegendEntry(h_DY_phi_sublead, "DY#rightarrow ee", "f");


        // Legend MC BKG        
        RP_pT->AddLegendEntry(h_bkg_pT[8], "DY#rightarrow #tau#tau", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[8], "DY#rightarrow #tau#tau", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[7], "t#bar{t}", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[7], "t#bar{t}", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[7], "t#bar{t}", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[7], "t#bar{t}", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[7], "t#bar{t}", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[7], "t#bar{t}", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[7], "t#bar{t}", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[7], "t#bar{t}", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[6], "tW", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[6], "tW", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[6], "tW", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[6], "tW", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[6], "tW", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[6], "tW", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[6], "tW", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[6], "tW", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[6], "tW", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[6], "tW", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[6], "tW", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[6], "tW", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[6], "tW", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[6], "tW", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[6], "tW", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[6], "tW", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[6], "tW", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[6], "tW", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[6], "tW", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[6], "tW", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[5], "#bar{t}W", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[5], "#bar{t}W", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[5], "#bar{t}W", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[5], "#bar{t}W", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[5], "#bar{t}W", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[5], "#bar{t}W", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[5], "#bar{t}W", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[5], "#bar{t}W", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[4], "ZZ", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[4], "ZZ", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[4], "ZZ", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[4], "ZZ", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[4], "ZZ", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[4], "ZZ", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[4], "ZZ", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[4], "ZZ", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[4], "ZZ", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[4], "ZZ", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[4], "ZZ", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[1], "ZZ", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[4], "ZZ", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[4], "ZZ", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[4], "ZZ", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[4], "ZZ", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[4], "ZZ", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[4], "ZZ", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[4], "ZZ", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[4], "ZZ", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[3], "WZ", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[3], "WZ", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[3], "WZ", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[3], "WZ", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[3], "WZ", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[3], "WZ", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[3], "WZ", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[3], "WZ", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[3], "WZ", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[3], "WZ", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[3], "WZ", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[3], "WZ", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[3], "WZ", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[3], "WZ", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[3], "WZ", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[3], "WZ", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[3], "WZ", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[3], "WZ", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[3], "WZ", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[3], "WZ", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[2], "WW", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[2], "WW", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[2], "WW", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[2], "WW", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[2], "WW", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[2], "WW", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[2], "WW", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[2], "WW", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[2], "WW", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[2], "WW", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[2], "WW", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[2], "WW", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[2], "WW", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[2], "WW", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[2], "WW", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[2], "WW", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[2], "WW", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[2], "WW", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[2], "WW", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[2], "WW", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[1], "W+Jets", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[1], "W+Jets", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[1], "W+Jets", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[1], "W+Jets", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[1], "W+Jets", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[1], "W+Jets", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[1], "W+Jets", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[1], "W+Jets", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[1], "W+Jets", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[1], "W+Jets", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[1], "W+Jets", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[1], "W+Jets", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[1], "W+Jets", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[1], "W+Jets", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[0], "QCD", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[0], "QCD", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[0], "QCD", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[0], "QCD", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[0], "QCD", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[0], "QCD", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[0], "QCD", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[0], "QCD", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[0], "QCD", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[0], "QCD", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[0], "QCD", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[0], "QCD", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[0], "QCD", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[0], "QCD", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[0], "QCD", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[0], "QCD", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[0], "QCD", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[0], "QCD", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[0], "QCD", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[0], "QCD", "f");

        RP_pT->Draw(0.5, 8e6, 0);
        RP_rapi->Draw(0.5, 8e6, 0);
        RP_pT_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_pT_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_pT_lead->Draw(0.5, 8e6, 0);
        RP_pT_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_pT_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_pT_sublead->Draw(0.5, 8e6, 0);
        RP_eta_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_eta_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_eta_lead->Draw(0.5, 8e6, 0);
        RP_eta_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_eta_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_eta_sublead->Draw(0.5, 8e6, 0);
        RP_phi_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_phi_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_phi_lead->Draw(0.5, 8e6, 0);
        RP_phi_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_phi_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_phi_sublead->Draw(0.5, 8e6, 0);

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nVTX #################################################

    if( whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[9], *h_bkg_nVTX_before_EffCorr[9], *h_bkg_nVTX[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter] );

            Color_t color = kBlack;
            if ( pr == _EE_QCDEMEnriched_Full ) color = kRed + 3;
            if ( pr == _EE_WJets ) color = kRed - 2;
            if ( pr == _EE_WW ) color = kMagenta - 5;
            if ( pr == _EE_WZ ) color = kMagenta - 2;
            if ( pr == _EE_ZZ ) color = kMagenta - 6;
            if ( pr == _EE_tbarW ) color = kGreen - 2;
            if ( pr == _EE_tW ) color = kGreen + 2;
            if ( pr == _EE_ttbar_Full ) color = kCyan + 2;
            if ( pr == _EE_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_nVTX[iter]->SetFillColor(color);

            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_nVTX[iter]->SetLineColor(color);

            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nVTX[iter]->SetDirectory(0);

            s_nVTX_before_PUCorr->Add( h_bkg_nVTX_before_PUCorr[iter] );
            s_nVTX_before_EffCorr->Add( h_bkg_nVTX_before_EffCorr[iter] );
            s_nVTX->Add( h_bkg_nVTX[iter] );

            iter++;

            if ( pr == _EE_QCDEMEnriched_Full ) // next - WJets
                pr = _EndOf_EE_WJets;
            if ( pr == _EE_WJets ) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if ( pr == _EE_tW ) // next -- ttbar
                pr = _EE_VVnST;
            if ( pr == _EE_DYTauTau_Full ) // last
                break;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nVTX_before_PUCorr, *h_DY_nVTX_before_EffCorr, *h_DY_nVTX;

        f_DY->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX_before_PUCorr );
        f_DY->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX_before_EffCorr );
        f_DY->GetObject( "h_nVTX_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX );

        h_DY_nVTX_before_PUCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetFillColor(kOrange);
        h_DY_nVTX->SetFillColor(kOrange);

        h_DY_nVTX_before_PUCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetLineColor(kOrange);
        h_DY_nVTX->SetLineColor(kOrange);

        h_DY_nVTX_before_PUCorr->SetDirectory(0);
        h_DY_nVTX_before_EffCorr->SetDirectory(0);
        h_DY_nVTX->SetDirectory(0);

        s_nVTX_before_PUCorr->Add( h_DY_nVTX_before_PUCorr );
        s_nVTX_before_EffCorr->Add( h_DY_nVTX_before_EffCorr );
        s_nVTX->Add( h_DY_nVTX );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;
        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], h_data_nVTX );

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t( "RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX );
        RP_nVTX_before_EffCorr = new myRatioPlot_t( "RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX );
        RP_nVTX = new myRatioPlot_t( "RP_nVTX", s_nVTX, h_data_nVTX );

        RP_nVTX_before_PUCorr->SetPlots("nVTX (before PU correction)", 0, 50);
        RP_nVTX_before_EffCorr->SetPlots("nVTX (before Efficiency correction)", 0, 50);
        RP_nVTX->SetPlots("nVTX", 0, 50);

        RP_nVTX_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nVTX_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nVTX->SetLegend(0.75, 0.4);

        // Legend data
        RP_nVTX_before_PUCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX->AddLegendEntry(h_data_nVTX, "Data", "lp");

        // Legend MC signal
        RP_nVTX_before_PUCorr->AddLegendEntry(h_DY_nVTX_before_PUCorr, "DY#rightarrow ee", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_DY_nVTX_before_EffCorr, "DY#rightarrow ee", "f");
        RP_nVTX->AddLegendEntry(h_DY_nVTX, "DY#rightarrow ee", "f");

        // Legend MC BKG
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[7], "t#bar{t}", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[7], "t#bar{t}", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[7], "t#bar{t}", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[6], "tW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[6], "tW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[6], "tW", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[5], "#bar{t}W", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[5], "#bar{t}W", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[5], "#bar{t}W", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[4], "ZZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[4], "ZZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[4], "ZZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[3], "WZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[3], "WZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[3], "WZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[2], "WW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[2], "WW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[2], "WW", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[1], "W+Jets", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[1], "W+Jets", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[1], "W+Jets", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[0], "QCD", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[0], "QCD", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[0], "QCD", "f");

        RP_nVTX_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nVTX->Draw(0.5, 3e6, 0);

    } // End of if(nVTX)

} // End of EE_HistDrawer()


/// -------------------------------- Muon Channel ------------------------------------ ///
void MuMu_HistDrawer ( TString whichGraphs , TString type)
{
    if ( !whichGraphs.Length() )
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc( _MuMu_DY_Full );
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+"_roccor.root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+"_roccor.root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+"_roccor.root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+"_roccor.root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+"_roccor.root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+"_roccor.root" << " opened successfully" << endl;

    // Files without Rochester correction
    Mgr.SetProc( _MuMu_DY_Full );
    TString name_DY_noRC = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY_noRC = new TFile( name_DY_noRC, "READ" );
    if (f_DY_noRC->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_Bkg_Full );
    TString name_bkg_noRC = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg_noRC = new TFile( name_bkg_noRC, "READ" );
    if (f_bkg_noRC->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_SingleMuon_Full );
    TString name_data_noRC = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data_noRC = new TFile( name_data_noRC, "READ" );
    if (f_data_noRC->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=2.25081e+07, MCintegral;

//################################# INVARIANT MASS #################################################

    if( whichGraphs=="ALL" || whichGraphs=="INVM" || whichGraphs=="INVMASS" )
    {
        count_drawn++;

        THStack *s_mass_fine_before_PUCorr, *s_mass_fine_before_RoccoR, *s_mass_fine_before_EffCorr, *s_mass_fine, *s_mass;
        s_mass_fine_before_PUCorr = new THStack("s_mass_fine_before_PUCorr", "");
        s_mass_fine_before_RoccoR = new THStack("s_mass_fine_before_RoccoR", "");
        s_mass_fine_before_EffCorr = new THStack("s_mass_fine_before_EffCorr", "");
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine_before_PUCorr[9], *h_bkg_mass_fine_before_RoccoR[9], *h_bkg_mass_fine_before_EffCorr[9],
             *h_bkg_mass_fine[9], *h_bkg_mass[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg_noRC->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            Color_t color = kBlack;
            if ( pr == _MuMu_QCDMuEnriched_Full ) color = kRed + 3;
            if ( pr == _MuMu_WJets ) color = kRed - 2;
            if ( pr == _MuMu_WW ) color = kMagenta - 5;
            if ( pr == _MuMu_WZ ) color = kMagenta - 2;
            if ( pr == _MuMu_ZZ ) color = kMagenta - 6;
            if ( pr == _MuMu_tbarW ) color = kGreen - 2;
            if ( pr == _MuMu_tW ) color = kGreen + 2;
            if ( pr == _MuMu_ttbar_Full ) color = kCyan + 2;
            if ( pr == _MuMu_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_fine[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetFillColor(color);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_fine[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetLineColor(color);

            h_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
            s_mass_fine_before_RoccoR->Add( h_bkg_mass_fine_before_RoccoR[iter] );
            s_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );

            iter++;

            if ( pr == _MuMu_QCDMuEnriched_Full )
                pr = _EndOf_MuMu_WJets; // next - WJets
            if ( pr == _MuMu_WJets )
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if ( pr == _MuMu_tW )
                pr = _MuMu_VVnST; // next - ttbar
            if ( pr == _MuMu_DYTauTau_Full ) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_fine_before_PUCorr, *h_DY_mass_fine_before_RoccoR, *h_DY_mass_fine_before_EffCorr, *h_DY_mass_fine, *h_DY_mass;

        f_DY_noRC->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine_before_PUCorr );
        f_DY_noRC->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine_before_RoccoR );
        f_DY->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine_before_EffCorr );
        f_DY->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine );
        f_DY->GetObject( "h_mass_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass );

        h_DY_mass_fine_before_PUCorr->SetFillColor(kOrange);
        h_DY_mass_fine_before_RoccoR->SetFillColor(kOrange);
        h_DY_mass_fine_before_EffCorr->SetFillColor(kOrange);
        h_DY_mass_fine->SetFillColor(kOrange);
        h_DY_mass->SetFillColor(kOrange);

        h_DY_mass_fine_before_PUCorr->SetLineColor(kOrange);
        h_DY_mass_fine_before_RoccoR->SetLineColor(kOrange);
        h_DY_mass_fine_before_EffCorr->SetLineColor(kOrange);
        h_DY_mass_fine->SetLineColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);

        h_DY_mass_fine_before_PUCorr->SetDirectory(0);
        h_DY_mass_fine_before_RoccoR->SetDirectory(0);
        h_DY_mass_fine_before_EffCorr->SetDirectory(0);
        h_DY_mass_fine->SetDirectory(0);
        h_DY_mass->SetDirectory(0);

        s_mass_fine_before_PUCorr->Add( h_DY_mass_fine_before_PUCorr );
        s_mass_fine_before_RoccoR->Add( h_DY_mass_fine_before_RoccoR );
        s_mass_fine_before_EffCorr->Add( h_DY_mass_fine_before_EffCorr );
        s_mass_fine->Add( h_DY_mass_fine );
        s_mass->Add( h_DY_mass );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine_before_PUCorr, *h_data_mass_fine_before_RoccoR, *h_data_mass_fine_before_EffCorr, *h_data_mass_fine, *h_data_mass;

        f_data_noRC->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine_before_PUCorr );
        f_data_noRC->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine_before_RoccoR );
        f_data->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass );

        h_data_mass_fine_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine_before_PUCorr->SetMarkerColor(kBlack);
        h_data_mass_fine_before_RoccoR->SetMarkerColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetMarkerColor(kBlack);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetLineColor(kBlack);
        h_data_mass_fine_before_RoccoR->SetLineColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetLineColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetDirectory(0);
        h_data_mass_fine_before_RoccoR->SetDirectory(0);
        h_data_mass_fine_before_EffCorr->SetDirectory(0);
        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_fine_before_PUCorr, *RP_mass_fine_before_RoccoR, *RP_mass_fine_before_EffCorr, *RP_mass_fine, *RP_mass;
        RP_mass_fine_before_PUCorr = new myRatioPlot_t( "RP_mass_fine_before_PUCorr", s_mass_fine_before_PUCorr, h_data_mass_fine_before_PUCorr );
        RP_mass_fine_before_RoccoR = new myRatioPlot_t( "RP_mass_fine_before_RoccoR", s_mass_fine_before_RoccoR, h_data_mass_fine_before_RoccoR );
        RP_mass_fine_before_EffCorr = new myRatioPlot_t( "RP_mass_fine_before_EffCorr", s_mass_fine_before_EffCorr, h_data_mass_fine_before_EffCorr );
        RP_mass_fine = new myRatioPlot_t( "RP_mass_fine", s_mass_fine, h_data_mass_fine );
        RP_mass = new myRatioPlot_t( "RP_mass", s_mass, h_data_mass );

        RP_mass_fine_before_PUCorr->SetPlots("Dimuon invariant mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_mass_fine_before_RoccoR->SetPlots("Dimuon invariant mass [GeV/c^{2}] (before Rochester correction)", 15, 3000);
        RP_mass_fine_before_EffCorr->SetPlots("Dimuon invariant mass [GeV/c^{2}] (before Efficiency SF)", 15, 3000);
        RP_mass_fine->SetPlots("Dimuon invariant mass [GeV/c^{2}]", 15, 3000);
        RP_mass->SetPlots("Dimuon invariant mass [GeV/c^{2}]", 15, 3000);

        RP_mass_fine_before_PUCorr->SetLegend(0.75, 0.4);
        RP_mass_fine_before_RoccoR->SetLegend(0.75, 0.4);
        RP_mass_fine_before_EffCorr->SetLegend(0.75, 0.4);
        RP_mass_fine->SetLegend(0.75, 0.4);
        RP_mass->SetLegend(0.75, 0.4);

        // Legend data
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_data_mass_fine_before_PUCorr, "Data", "lp");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_data_mass_fine_before_RoccoR, "Data", "lp");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_data_mass_fine_before_EffCorr, "Data", "lp");
        RP_mass_fine->AddLegendEntry(h_data_mass_fine, "Data", "lp");
        RP_mass->AddLegendEntry(h_data_mass, "Data", "lp");

        // Legend MC signal
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_DY_mass_fine_before_PUCorr, "DY#rightarrow#mu#mu", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_DY_mass_fine_before_RoccoR, "DY#rightarrow#mu#mu", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_DY_mass_fine_before_EffCorr, "DY#rightarrow#mu#mu", "f");
        RP_mass_fine->AddLegendEntry(h_DY_mass_fine, "DY#rightarrow#mu#mu", "f");
        RP_mass->AddLegendEntry(h_DY_mass, "DY#rightarrow#mu#mu", "f");

        // Legend MC BKG
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[8], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[8], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[7], "t#bar{t}", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[7], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[7], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[7], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[7], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[6], "tW", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[6], "tW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[6], "tW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[6], "tW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[6], "tW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[5], "#bar{t}W", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[5], "#bar{t}W", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[5], "#bar{t}W", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[5], "#bar{t}W", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[5], "#bar{t}W", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4], "ZZ", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[4], "ZZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4], "ZZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4], "ZZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4], "ZZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "WZ", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[3], "WZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "WZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "WZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "WZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2], "WW", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[2], "WW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2], "WW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2], "WW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2], "WW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "W+Jets", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[1], "W+Jets", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "W+Jets", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "W+Jets", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "W+Jets", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "QCD", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[0], "QCD", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "QCD", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "QCD", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "QCD", "f");

        RP_mass_fine_before_PUCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_RoccoR->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_EffCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine->Draw(0.5, 3e6, 1);
        RP_mass->Draw(0.5, 1e7, 1);

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################

    if ( whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI") )
    {
        count_drawn++;

        THStack *s_pT, *s_rapi;
        s_rapi = new THStack( "s_rapi", "" );
        s_pT = new THStack( "s_pT", "" );

        THStack *s_pT_lead_before_PUCorr = new THStack( "s_pT_lead_before_PUCorr", "" );
        THStack *s_pT_sublead_before_PUCorr = new THStack( "s_pT_sublead_before_PUCorr", "" );
        THStack *s_eta_lead_before_PUCorr = new THStack( "s_eta_lead_before_PUCorr", "" );
        THStack *s_eta_sublead_before_PUCorr = new THStack( "s_eta_sublead_before_PUCorr", "" );
        THStack *s_phi_lead_before_PUCorr = new THStack( "s_phi_lead_before_PUCorr", "" );
        THStack *s_phi_sublead_before_PUCorr = new THStack( "s_phi_sublead_before_PUCorr", "" );

        THStack *s_pT_lead_before_RoccoR = new THStack( "s_pT_lead_before_RoccoR", "" );
        THStack *s_pT_sublead_before_RoccoR = new THStack( "s_pT_sublead_before_RoccoR", "" );
        THStack *s_eta_lead_before_RoccoR = new THStack( "s_eta_lead_before_RoccoR", "" );
        THStack *s_eta_sublead_before_RoccoR = new THStack( "s_eta_sublead_before_RoccoR", "" );
        THStack *s_phi_lead_before_RoccoR = new THStack( "s_phi_lead_before_RoccoR", "" );
        THStack *s_phi_sublead_before_RoccoR = new THStack( "s_phi_sublead_before_RoccoR", "" );

        THStack *s_pT_lead_before_EffCorr = new THStack( "s_pT_lead_before_EffCorr", "" );
        THStack *s_pT_sublead_before_EffCorr = new THStack( "s_pT_sublead_before_EffCorr", "" );
        THStack *s_eta_lead_before_EffCorr = new THStack( "s_eta_lead_before_EffCorr", "" );
        THStack *s_eta_sublead_before_EffCorr = new THStack( "s_eta_sublead_before_EffCorr", "" );
        THStack *s_phi_lead_before_EffCorr = new THStack( "s_phi_lead_before_EffCorr", "" );
        THStack *s_phi_sublead_before_EffCorr = new THStack( "s_phi_sublead_before_EffCorr", "" );

        THStack *s_pT_lead = new THStack( "s_pT_lead", "" );
        THStack *s_pT_sublead = new THStack( "s_pT_sublead", "" );
        THStack *s_eta_lead = new THStack( "s_eta_lead", "" );
        THStack *s_eta_sublead = new THStack( "s_eta_sublead", "" );
        THStack *s_phi_lead = new THStack( "s_phi_lead", "" );
        THStack *s_phi_sublead = new THStack( "s_phi_sublead", "" );

//----------------------------------- MC bkg -------------------------------------------------------

        TH1D *h_bkg_pT[9], *h_bkg_rapi[9], *h_bkg_pT_lead_before_PUCorr[9], *h_bkg_pT_lead_before_RoccoR[9], *h_bkg_pT_lead_before_EffCorr[9], *h_bkg_pT_lead[9],
             *h_bkg_pT_sublead_before_PUCorr[9],*h_bkg_pT_sublead_before_RoccoR[9],  *h_bkg_pT_sublead_before_EffCorr[9], *h_bkg_pT_sublead[9],
             *h_bkg_eta_lead_before_PUCorr[9], *h_bkg_eta_lead_before_RoccoR[9], *h_bkg_eta_lead_before_EffCorr[9], *h_bkg_eta_lead[9],
             *h_bkg_eta_sublead_before_PUCorr[9], *h_bkg_eta_sublead_before_RoccoR[9], *h_bkg_eta_sublead_before_EffCorr[9], *h_bkg_eta_sublead[9],
             *h_bkg_phi_lead_before_PUCorr[9], *h_bkg_phi_lead_before_RoccoR[9], *h_bkg_phi_lead_before_EffCorr[9], *h_bkg_phi_lead[9],
             *h_bkg_phi_sublead_before_PUCorr[9], *h_bkg_phi_sublead_before_RoccoR[9], *h_bkg_phi_sublead_before_EffCorr[9], *h_bkg_phi_sublead[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg->GetObject( "h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter] );
            f_bkg->GetObject( "h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter] );
            f_bkg_noRC->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_pT_lead_"+Mgr.Procname[pr], h_bkg_pT_lead[iter] );
            f_bkg_noRC->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_pT_sublead_"+Mgr.Procname[pr], h_bkg_pT_sublead[iter] );
            f_bkg_noRC->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_eta_lead_"+Mgr.Procname[pr], h_bkg_eta_lead[iter] );
            f_bkg_noRC->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_eta_sublead_"+Mgr.Procname[pr], h_bkg_eta_sublead[iter] );
            f_bkg_noRC->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_phi_lead_"+Mgr.Procname[pr], h_bkg_phi_lead[iter] );
            f_bkg_noRC->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_RoccoR[iter] );
            f_bkg->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_EffCorr[iter] );
            f_bkg->GetObject( "h_phi_sublead_"+Mgr.Procname[pr], h_bkg_phi_sublead[iter] );

            Color_t color = kBlack;
            if ( pr == _MuMu_QCDMuEnriched_Full ) color = kRed + 3;
            if ( pr == _MuMu_WJets ) color = kRed - 2;
            if ( pr == _MuMu_WW ) color = kMagenta - 5;
            if ( pr == _MuMu_WZ ) color = kMagenta - 2;
            if ( pr == _MuMu_ZZ ) color = kMagenta - 6;
            if ( pr == _MuMu_tbarW ) color = kGreen - 2;
            if ( pr == _MuMu_tW ) color = kGreen + 2;
            if ( pr == _MuMu_ttbar_Full ) color = kCyan + 2;
            if ( pr == _MuMu_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_pT[iter]->SetFillColor(color);
            h_bkg_rapi[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead[iter]->SetFillColor(color);

            h_bkg_pT[iter]->SetLineColor(color);
            h_bkg_rapi[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead[iter]->SetLineColor(color);

            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead[iter]->SetDirectory(0);

            s_pT->Add( h_bkg_pT[iter] );
            s_rapi->Add( h_bkg_rapi[iter] );
            s_pT_lead_before_PUCorr->Add( h_bkg_pT_lead_before_PUCorr[iter] );
            s_pT_lead_before_RoccoR->Add( h_bkg_pT_lead_before_RoccoR[iter] );
            s_pT_lead_before_EffCorr->Add( h_bkg_pT_lead_before_EffCorr[iter] );
            s_pT_lead->Add( h_bkg_pT_lead[iter] );
            s_pT_sublead_before_PUCorr->Add( h_bkg_pT_sublead_before_PUCorr[iter] );
            s_pT_sublead_before_RoccoR->Add( h_bkg_pT_sublead_before_RoccoR[iter] );
            s_pT_sublead_before_EffCorr->Add( h_bkg_pT_sublead_before_EffCorr[iter] );
            s_pT_sublead->Add( h_bkg_pT_sublead[iter] );
            s_eta_lead_before_PUCorr->Add( h_bkg_eta_lead_before_PUCorr[iter] );
            s_eta_lead_before_RoccoR->Add( h_bkg_eta_lead_before_RoccoR[iter] );
            s_eta_lead_before_EffCorr->Add( h_bkg_eta_lead_before_EffCorr[iter] );
            s_eta_lead->Add( h_bkg_eta_lead[iter] );
            s_eta_sublead_before_PUCorr->Add( h_bkg_eta_sublead_before_PUCorr[iter] );
            s_eta_sublead_before_RoccoR->Add( h_bkg_eta_sublead_before_RoccoR[iter] );
            s_eta_sublead_before_EffCorr->Add( h_bkg_eta_sublead_before_EffCorr[iter] );
            s_eta_sublead->Add( h_bkg_eta_sublead[iter] );
            s_phi_lead_before_PUCorr->Add( h_bkg_phi_lead_before_PUCorr[iter] );
            s_phi_lead_before_RoccoR->Add( h_bkg_phi_lead_before_RoccoR[iter] );
            s_phi_lead_before_EffCorr->Add( h_bkg_phi_lead_before_EffCorr[iter] );
            s_phi_lead->Add( h_bkg_phi_lead[iter] );
            s_phi_sublead_before_PUCorr->Add( h_bkg_phi_sublead_before_PUCorr[iter] );
            s_phi_sublead_before_RoccoR->Add( h_bkg_phi_sublead_before_RoccoR[iter] );
            s_phi_sublead_before_EffCorr->Add( h_bkg_phi_sublead_before_EffCorr[iter] );
            s_phi_sublead->Add( h_bkg_phi_sublead[iter] );

            iter++;

            if ( pr == _MuMu_QCDMuEnriched_Full )
                pr = _EndOf_MuMu_WJets; // next - WJets
            if ( pr == _MuMu_WJets )
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if ( pr == _MuMu_tW )
                pr = _MuMu_VVnST; // next - ttbar
            if ( pr == _MuMu_DYTauTau_Full ) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_pT, *h_DY_rapi, *h_DY_pT_lead_before_PUCorr, *h_DY_pT_lead_before_RoccoR, *h_DY_pT_lead_before_EffCorr, *h_DY_pT_lead,
             *h_DY_pT_sublead_before_PUCorr, *h_DY_pT_sublead_before_RoccoR, *h_DY_pT_sublead_before_EffCorr, *h_DY_pT_sublead,
             *h_DY_eta_lead_before_PUCorr, *h_DY_eta_lead_before_RoccoR, *h_DY_eta_lead_before_EffCorr, *h_DY_eta_lead,
             *h_DY_eta_sublead_before_PUCorr, *h_DY_eta_sublead_before_RoccoR,*h_DY_eta_sublead_before_EffCorr, *h_DY_eta_sublead,
             *h_DY_phi_lead_before_PUCorr, *h_DY_phi_lead_before_RoccoR, *h_DY_phi_lead_before_EffCorr, *h_DY_phi_lead,
             *h_DY_phi_sublead_before_PUCorr, *h_DY_phi_sublead_before_RoccoR, *h_DY_phi_sublead_before_EffCorr, *h_DY_phi_sublead;

        f_DY->GetObject( "h_pT_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT );
        f_DY->GetObject( "h_rapi_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi );
        f_DY_noRC->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_PUCorr );
        f_DY_noRC->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_RoccoR );
        f_DY->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_EffCorr );
        f_DY->GetObject( "h_pT_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead );
        f_DY_noRC->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_PUCorr );
        f_DY_noRC->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_RoccoR );
        f_DY->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_EffCorr );
        f_DY->GetObject( "h_pT_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead );
        f_DY_noRC->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_PUCorr );
        f_DY_noRC->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_RoccoR );
        f_DY->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_EffCorr );
        f_DY->GetObject( "h_eta_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead );
        f_DY_noRC->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_PUCorr );
        f_DY_noRC->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_RoccoR );
        f_DY->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_EffCorr );
        f_DY->GetObject( "h_eta_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead );
        f_DY_noRC->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_PUCorr );
        f_DY_noRC->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_RoccoR );
        f_DY->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_EffCorr );
        f_DY->GetObject( "h_phi_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead );
        f_DY_noRC->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_PUCorr );
        f_DY_noRC->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_RoccoR );
        f_DY->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_EffCorr );
        f_DY->GetObject( "h_phi_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead );

        h_DY_pT->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_lead->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_RoccoR->SetFillColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_sublead->SetFillColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_lead->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_RoccoR->SetFillColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_sublead->SetFillColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_lead->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_RoccoR->SetFillColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_sublead->SetFillColor(kOrange);

        h_DY_pT->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_lead->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_RoccoR->SetLineColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_sublead->SetLineColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_lead->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_RoccoR->SetLineColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_sublead->SetLineColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_lead->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_RoccoR->SetLineColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_sublead->SetLineColor(kOrange);

        h_DY_pT->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT_lead_before_PUCorr->SetDirectory(0);
        h_DY_pT_lead_before_RoccoR->SetDirectory(0);
        h_DY_pT_lead_before_EffCorr->SetDirectory(0);
        h_DY_pT_lead->SetDirectory(0);
        h_DY_pT_sublead_before_PUCorr->SetDirectory(0);
        h_DY_pT_sublead_before_RoccoR->SetDirectory(0);
        h_DY_pT_sublead_before_EffCorr->SetDirectory(0);
        h_DY_pT_sublead->SetDirectory(0);
        h_DY_eta_lead_before_PUCorr->SetDirectory(0);
        h_DY_eta_lead_before_RoccoR->SetDirectory(0);
        h_DY_eta_lead_before_EffCorr->SetDirectory(0);
        h_DY_eta_lead->SetDirectory(0);
        h_DY_eta_sublead_before_PUCorr->SetDirectory(0);
        h_DY_eta_sublead_before_RoccoR->SetDirectory(0);
        h_DY_eta_sublead_before_EffCorr->SetDirectory(0);
        h_DY_eta_sublead->SetDirectory(0);
        h_DY_phi_lead_before_PUCorr->SetDirectory(0);
        h_DY_phi_lead_before_RoccoR->SetDirectory(0);
        h_DY_phi_lead_before_EffCorr->SetDirectory(0);
        h_DY_phi_lead->SetDirectory(0);
        h_DY_phi_sublead_before_PUCorr->SetDirectory(0);
        h_DY_phi_sublead_before_RoccoR->SetDirectory(0);
        h_DY_phi_sublead_before_EffCorr->SetDirectory(0);
        h_DY_phi_sublead->SetDirectory(0);

        s_pT->Add( h_DY_pT );
        s_rapi->Add( h_DY_rapi );
        s_pT_lead_before_PUCorr->Add( h_DY_pT_lead_before_PUCorr );
        s_pT_lead_before_RoccoR->Add( h_DY_pT_lead_before_RoccoR );
        s_pT_lead_before_EffCorr->Add( h_DY_pT_lead_before_EffCorr );
        s_pT_lead->Add( h_DY_pT_lead );
        s_pT_sublead_before_PUCorr->Add( h_DY_pT_sublead_before_PUCorr );
        s_pT_sublead_before_RoccoR->Add( h_DY_pT_sublead_before_RoccoR );
        s_pT_sublead_before_EffCorr->Add( h_DY_pT_sublead_before_EffCorr );
        s_pT_sublead->Add( h_DY_pT_sublead );
        s_eta_lead_before_PUCorr->Add( h_DY_eta_lead_before_PUCorr );
        s_eta_lead_before_RoccoR->Add( h_DY_eta_lead_before_RoccoR );
        s_eta_lead_before_EffCorr->Add( h_DY_eta_lead_before_EffCorr );
        s_eta_lead->Add( h_DY_eta_lead );
        s_eta_sublead_before_PUCorr->Add( h_DY_eta_sublead_before_PUCorr );
        s_eta_sublead_before_RoccoR->Add( h_DY_eta_sublead_before_RoccoR );
        s_eta_sublead_before_EffCorr->Add( h_DY_eta_sublead_before_EffCorr );
        s_eta_sublead->Add( h_DY_eta_sublead );
        s_phi_lead_before_PUCorr->Add( h_DY_phi_lead_before_PUCorr );
        s_phi_lead_before_RoccoR->Add( h_DY_phi_lead_before_RoccoR );
        s_phi_lead_before_EffCorr->Add( h_DY_phi_lead_before_EffCorr );
        s_phi_lead->Add( h_DY_phi_lead );
        s_phi_sublead_before_PUCorr->Add( h_DY_phi_sublead_before_PUCorr );
        s_phi_sublead_before_RoccoR->Add( h_DY_phi_sublead_before_RoccoR );
        s_phi_sublead_before_EffCorr->Add( h_DY_phi_sublead_before_EffCorr );
        s_phi_sublead->Add( h_DY_phi_sublead );


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT, *h_data_rapi, *h_data_pT_lead_before_PUCorr, *h_data_pT_lead_before_RoccoR, *h_data_pT_lead_before_EffCorr, *h_data_pT_lead,
             *h_data_pT_sublead_before_PUCorr, *h_data_pT_sublead_before_RoccoR, *h_data_pT_sublead_before_EffCorr, *h_data_pT_sublead,
             *h_data_eta_lead_before_PUCorr, *h_data_eta_lead_before_RoccoR, *h_data_eta_lead_before_EffCorr, *h_data_eta_lead,
             *h_data_eta_sublead_before_PUCorr, *h_data_eta_sublead_before_RoccoR, *h_data_eta_sublead_before_EffCorr, *h_data_eta_sublead,
             *h_data_phi_lead_before_PUCorr, *h_data_phi_lead_before_RoccoR, *h_data_phi_lead_before_EffCorr, *h_data_phi_lead,
             *h_data_phi_sublead_before_PUCorr, *h_data_phi_sublead_before_RoccoR, *h_data_phi_sublead_before_EffCorr, *h_data_phi_sublead;

        f_data->GetObject( "h_pT_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT );
        f_data->GetObject( "h_rapi_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_rapi );
        f_data_noRC->GetObject( "h_pT_lead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead_before_PUCorr );
        f_data_noRC->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead_before_RoccoR );
        f_data->GetObject( "h_pT_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead_before_EffCorr );
        f_data->GetObject( "h_pT_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead );
        f_data_noRC->GetObject( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead_before_PUCorr );
        f_data_noRC->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead_before_RoccoR );
        f_data->GetObject( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead_before_EffCorr );
        f_data->GetObject( "h_pT_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead );
        f_data_noRC->GetObject( "h_eta_lead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead_before_PUCorr );
        f_data_noRC->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead_before_RoccoR );
        f_data->GetObject( "h_eta_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead_before_EffCorr );
        f_data->GetObject( "h_eta_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead );
        f_data_noRC->GetObject( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead_before_PUCorr );
        f_data_noRC->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead_before_RoccoR );
        f_data->GetObject( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead_before_EffCorr );
        f_data->GetObject( "h_eta_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead );
        f_data_noRC->GetObject( "h_phi_lead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead_before_PUCorr );
        f_data_noRC->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead_before_RoccoR );
        f_data->GetObject( "h_phi_lead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead_before_EffCorr );
        f_data->GetObject( "h_phi_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead );
        f_data_noRC->GetObject( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead_before_PUCorr );
        f_data_noRC->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead_before_RoccoR );
        f_data->GetObject( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead_before_EffCorr );
        f_data->GetObject( "h_phi_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead );

        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead->SetMarkerStyle(kFullDotLarge);

        h_data_pT->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);
        h_data_pT_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_pT_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_pT_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_pT_lead->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_pT_sublead->SetMarkerColor(kBlack);
        h_data_eta_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_eta_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_eta_lead->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_eta_sublead->SetMarkerColor(kBlack);
        h_data_phi_lead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_phi_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_lead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_phi_lead->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_PUCorr->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_EffCorr->SetMarkerColor(kBlack);
        h_data_phi_sublead->SetMarkerColor(kBlack);

        h_data_pT->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);
        h_data_pT_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_pT_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_pT_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_pT_lead->SetLineColor(kBlack);
        h_data_pT_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_pT_sublead_before_RoccoR->SetLineColor(kBlack);
        h_data_pT_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_pT_sublead->SetLineColor(kBlack);
        h_data_eta_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_eta_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_eta_lead->SetLineColor(kBlack);
        h_data_eta_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_eta_sublead_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_eta_sublead->SetLineColor(kBlack);
        h_data_phi_lead_before_PUCorr->SetLineColor(kBlack);
        h_data_phi_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_lead_before_EffCorr->SetLineColor(kBlack);
        h_data_phi_lead->SetLineColor(kBlack);
        h_data_phi_sublead_before_PUCorr->SetLineColor(kBlack);
        h_data_phi_sublead_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_sublead_before_EffCorr->SetLineColor(kBlack);
        h_data_phi_sublead->SetLineColor(kBlack);

        h_data_pT->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT_lead_before_PUCorr->SetDirectory(0);
        h_data_pT_lead_before_RoccoR->SetDirectory(0);
        h_data_pT_lead_before_EffCorr->SetDirectory(0);
        h_data_pT_lead->SetDirectory(0);
        h_data_pT_sublead_before_PUCorr->SetDirectory(0);
        h_data_pT_sublead_before_RoccoR->SetDirectory(0);
        h_data_pT_sublead_before_EffCorr->SetDirectory(0);
        h_data_pT_sublead->SetDirectory(0);
        h_data_eta_lead_before_PUCorr->SetDirectory(0);
        h_data_eta_lead_before_RoccoR->SetDirectory(0);
        h_data_eta_lead_before_EffCorr->SetDirectory(0);
        h_data_eta_lead->SetDirectory(0);
        h_data_eta_sublead_before_PUCorr->SetDirectory(0);
        h_data_eta_sublead_before_RoccoR->SetDirectory(0);
        h_data_eta_sublead_before_EffCorr->SetDirectory(0);
        h_data_eta_sublead->SetDirectory(0);
        h_data_phi_lead_before_PUCorr->SetDirectory(0);
        h_data_phi_lead_before_RoccoR->SetDirectory(0);
        h_data_phi_lead_before_EffCorr->SetDirectory(0);
        h_data_phi_lead->SetDirectory(0);
        h_data_phi_sublead_before_PUCorr->SetDirectory(0);
        h_data_phi_sublead_before_RoccoR->SetDirectory(0);
        h_data_phi_sublead_before_EffCorr->SetDirectory(0);
        h_data_phi_sublead->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_pT, *RP_rapi, *RP_pT_lead_before_PUCorr, *RP_pT_lead_before_RoccoR, *RP_pT_lead_before_EffCorr, *RP_pT_lead,
                      *RP_pT_sublead_before_PUCorr, *RP_pT_sublead_before_RoccoR, *RP_pT_sublead_before_EffCorr, *RP_pT_sublead,
                      *RP_eta_lead_before_PUCorr, *RP_eta_lead_before_RoccoR, *RP_eta_lead_before_EffCorr, *RP_eta_lead,
                      *RP_eta_sublead_before_PUCorr, *RP_eta_sublead_before_RoccoR, *RP_eta_sublead_before_EffCorr, *RP_eta_sublead,
                      *RP_phi_lead_before_PUCorr, *RP_phi_lead_before_RoccoR, *RP_phi_lead_before_EffCorr, *RP_phi_lead,
                      *RP_phi_sublead_before_PUCorr, *RP_phi_sublead_before_RoccoR, *RP_phi_sublead_before_EffCorr, *RP_phi_sublead;
        RP_pT = new myRatioPlot_t( "RP_pT", s_pT, h_data_pT );
        RP_rapi = new myRatioPlot_t( "RP_rapi", s_rapi, h_data_rapi );
        RP_pT_lead_before_PUCorr = new myRatioPlot_t( "RP_pT_lead_before_PUCorr", s_pT_lead_before_PUCorr, h_data_pT_lead_before_PUCorr );
        RP_pT_lead_before_RoccoR = new myRatioPlot_t( "RP_pT_lead_before_RoccoR", s_pT_lead_before_RoccoR, h_data_pT_lead_before_RoccoR );
        RP_pT_lead_before_EffCorr = new myRatioPlot_t( "RP_pT_lead_before_EffCorr", s_pT_lead_before_EffCorr, h_data_pT_lead_before_EffCorr );
        RP_pT_lead = new myRatioPlot_t( "RP_pT_lead", s_pT_lead, h_data_pT_lead );
        RP_pT_sublead_before_PUCorr = new myRatioPlot_t( "RP_pT_sublead_before_PUCorr", s_pT_sublead_before_PUCorr, h_data_pT_sublead_before_PUCorr );
        RP_pT_sublead_before_RoccoR = new myRatioPlot_t( "RP_pT_sublead_before_RoccoR", s_pT_sublead_before_RoccoR, h_data_pT_sublead_before_RoccoR );
        RP_pT_sublead_before_EffCorr = new myRatioPlot_t( "RP_pT_sublead_before_EffCorr", s_pT_sublead_before_EffCorr, h_data_pT_sublead_before_EffCorr );
        RP_pT_sublead = new myRatioPlot_t( "RP_pT_sublead", s_pT_sublead, h_data_pT_sublead );
        RP_eta_lead_before_PUCorr = new myRatioPlot_t( "RP_eta_lead_before_PUCorr", s_eta_lead_before_PUCorr, h_data_eta_lead_before_PUCorr );
        RP_eta_lead_before_RoccoR = new myRatioPlot_t( "RP_eta_lead_before_RoccoR", s_eta_lead_before_RoccoR, h_data_eta_lead_before_RoccoR );
        RP_eta_lead_before_EffCorr = new myRatioPlot_t( "RP_eta_lead_before_EffCorr", s_eta_lead_before_EffCorr, h_data_eta_lead_before_EffCorr );
        RP_eta_lead = new myRatioPlot_t( "RP_eta_lead", s_eta_lead, h_data_eta_lead );
        RP_eta_sublead_before_PUCorr = new myRatioPlot_t( "RP_eta_sublead_before_PUCorr", s_eta_sublead_before_PUCorr, h_data_eta_sublead_before_PUCorr );
        RP_eta_sublead_before_RoccoR = new myRatioPlot_t( "RP_eta_sublead_before_RoccoR", s_eta_sublead_before_RoccoR, h_data_eta_sublead_before_RoccoR );
        RP_eta_sublead_before_EffCorr = new myRatioPlot_t( "RP_eta_sublead_before_EffCorr", s_eta_sublead_before_EffCorr, h_data_eta_sublead_before_EffCorr );
        RP_eta_sublead = new myRatioPlot_t( "RP_eta_sublead", s_eta_sublead, h_data_eta_sublead );
        RP_phi_lead_before_PUCorr = new myRatioPlot_t( "RP_phi_lead_before_PUCorr", s_phi_lead_before_PUCorr, h_data_phi_lead_before_PUCorr );
        RP_phi_lead_before_RoccoR = new myRatioPlot_t( "RP_phi_lead_before_RoccoR", s_phi_lead_before_RoccoR, h_data_phi_lead_before_RoccoR );
        RP_phi_lead_before_EffCorr = new myRatioPlot_t( "RP_phi_lead_before_EffCorr", s_phi_lead_before_EffCorr, h_data_phi_lead_before_EffCorr );
        RP_phi_lead = new myRatioPlot_t( "RP_phi_lead", s_phi_lead, h_data_phi_lead );
        RP_phi_sublead_before_PUCorr = new myRatioPlot_t( "RP_phi_sublead_before_PUCorr", s_phi_sublead_before_PUCorr, h_data_phi_sublead_before_PUCorr );
        RP_phi_sublead_before_RoccoR = new myRatioPlot_t( "RP_phi_sublead_before_RoccoR", s_phi_sublead_before_RoccoR, h_data_phi_sublead_before_RoccoR );
        RP_phi_sublead_before_EffCorr = new myRatioPlot_t( "RP_phi_sublead_before_EffCorr", s_phi_sublead_before_EffCorr, h_data_phi_sublead_before_EffCorr );
        RP_phi_sublead = new myRatioPlot_t( "RP_phi_sublead", s_phi_sublead, h_data_phi_sublead );

        RP_pT->SetPlots("Dimuon p_{T} [GeV/c]", 0, 1000);
        RP_rapi->SetPlots("Dimuon rapidity", -3, 3);
        RP_pT_lead_before_PUCorr->SetPlots("Leading muon p_{T} [GeV/c] (before PU correction)", 0, 1000);
        RP_pT_lead_before_RoccoR->SetPlots("Leading muon p_{T} [GeV/c] (before Rochester correction)", 0, 1000);
        RP_pT_lead_before_EffCorr->SetPlots("Leading muon p_{T} [GeV/c] (before Efficiency SF)", 0, 1000);
        RP_pT_lead->SetPlots("Leading muon p_{T} [GeV/c]", 0, 1000);
        RP_pT_sublead_before_PUCorr->SetPlots("Subleading muon p_{T} [GeV/c] (before PU correction)", 0, 1000);
        RP_pT_sublead_before_RoccoR->SetPlots("Subleading muon p_{T} [GeV/c] (before Rochester correction)", 0, 1000);
        RP_pT_sublead_before_EffCorr->SetPlots("Subleading muon p_{T} [GeV/c] (before Efficiency SF)", 0, 1000);
        RP_pT_sublead->SetPlots("Subleading muon p_{T} [GeV/c]", 0, 1000);
        RP_eta_lead_before_PUCorr->SetPlots("Leading muon #eta (before PU correction)", -4, 4);
        RP_eta_lead_before_RoccoR->SetPlots("Leading muon #eta (before Rochester correction)", -4, 4);
        RP_eta_lead_before_EffCorr->SetPlots("Leading muon #eta (before Efficiency SF)", -4, 4);
        RP_eta_lead->SetPlots("Leading muon #eta", -4, 4);
        RP_eta_sublead_before_PUCorr->SetPlots("Subleading muon #eta (before PU correction)", -4, 4);
        RP_eta_sublead_before_RoccoR->SetPlots("Subleading muon #eta (before Rochester correction)", -4, 4);
        RP_eta_sublead_before_EffCorr->SetPlots("Subleading muon #eta (before Efficiency SF)", -4, 4);
        RP_eta_sublead->SetPlots("Subleading muon #eta", -4, 4);
        RP_phi_lead_before_PUCorr->SetPlots("Leading muon #phi (before PU correction)", -4, 4);
        RP_phi_lead_before_RoccoR->SetPlots("Leading muon #phi (before Rochester correction)", -4, 4);
        RP_phi_lead_before_EffCorr->SetPlots("Leading muon #phi (before Efficiency SF)", -4, 4);
        RP_phi_lead->SetPlots("Leading muon #phi", -4, 4);
        RP_phi_sublead_before_PUCorr->SetPlots("Subleading muon #phi (before PU correction)", -4, 4);
        RP_phi_sublead_before_RoccoR->SetPlots("Subleading muon #phi (before Rochester correction)", -4, 4);
        RP_phi_sublead_before_EffCorr->SetPlots("Subleading muon #phi (before Efficiency SF)", -4, 4);
        RP_phi_sublead->SetPlots("Subleading muon #phi", -4, 4);

        RP_pT->SetLegend(0.75, 0.4);
        RP_rapi->SetLegend(0.75, 0.4);
        RP_pT_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_pT_lead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_pT_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_pT_lead->SetLegend(0.75, 0.4);
        RP_pT_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_pT_sublead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_pT_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_pT_sublead->SetLegend(0.75, 0.4);
        RP_eta_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_eta_lead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_eta_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_eta_lead->SetLegend(0.75, 0.4);
        RP_eta_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_eta_sublead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_eta_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_eta_sublead->SetLegend(0.75, 0.4);
        RP_phi_lead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_phi_lead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_phi_lead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_phi_lead->SetLegend(0.75, 0.4);
        RP_phi_sublead_before_PUCorr->SetLegend(0.75, 0.4);
        RP_phi_sublead_before_RoccoR->SetLegend(0.75, 0.4);
        RP_phi_sublead_before_EffCorr->SetLegend(0.75, 0.4);
        RP_phi_sublead->SetLegend(0.75, 0.4);


        // Legend data
        RP_pT->AddLegendEntry(h_data_pT, "Data", "lp");
        RP_rapi->AddLegendEntry(h_data_rapi, "Data", "lp");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_data_pT_lead_before_PUCorr, "Data", "lp");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_data_pT_lead_before_RoccoR, "Data", "lp");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_data_pT_lead_before_EffCorr, "Data", "lp");
        RP_pT_lead->AddLegendEntry(h_data_pT_lead, "Data", "lp");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_data_pT_sublead_before_PUCorr, "Data", "lp");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_data_pT_sublead_before_RoccoR, "Data", "lp");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_data_pT_sublead_before_EffCorr, "Data", "lp");
        RP_pT_sublead->AddLegendEntry(h_data_pT_sublead, "Data", "lp");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_data_eta_lead_before_PUCorr, "Data", "lp");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_data_eta_lead_before_RoccoR, "Data", "lp");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_data_eta_lead_before_EffCorr, "Data", "lp");
        RP_eta_lead->AddLegendEntry(h_data_eta_lead, "Data", "lp");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_data_eta_sublead_before_PUCorr, "Data", "lp");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_data_eta_sublead_before_RoccoR, "Data", "lp");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_data_eta_sublead_before_EffCorr, "Data", "lp");
        RP_eta_sublead->AddLegendEntry(h_data_eta_sublead, "Data", "lp");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_data_phi_lead_before_PUCorr, "Data", "lp");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_data_phi_lead_before_RoccoR, "Data", "lp");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_data_phi_lead_before_EffCorr, "Data", "lp");
        RP_phi_lead->AddLegendEntry(h_data_phi_lead, "Data", "lp");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_data_phi_sublead_before_PUCorr, "Data", "lp");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_data_phi_sublead_before_RoccoR, "Data", "lp");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_data_phi_sublead_before_EffCorr, "Data", "lp");
        RP_phi_sublead->AddLegendEntry(h_data_phi_sublead, "Data", "lp");


        // Legend MC signal
        RP_pT->AddLegendEntry(h_DY_pT, "DY#rightarrow #mu#mu", "f");
        RP_rapi->AddLegendEntry(h_DY_rapi, "DY#rightarrow #mu#mu", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_DY_pT_lead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_DY_pT_lead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_DY_pT_lead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_pT_lead->AddLegendEntry(h_DY_pT_lead, "DY#rightarrow #mu#mu", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_DY_pT_sublead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_DY_pT_sublead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_DY_pT_sublead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_pT_sublead->AddLegendEntry(h_DY_pT_sublead, "DY#rightarrow #mu#mu", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_DY_eta_lead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_DY_eta_lead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_DY_eta_lead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_eta_lead->AddLegendEntry(h_DY_eta_lead, "DY#rightarrow #mu#mu", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_DY_eta_sublead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_DY_eta_sublead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_DY_eta_sublead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_eta_sublead->AddLegendEntry(h_DY_eta_sublead, "DY#rightarrow #mu#mu", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_DY_phi_lead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_DY_phi_lead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_DY_phi_lead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_phi_lead->AddLegendEntry(h_DY_phi_lead, "DY#rightarrow #mu#mu", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_DY_phi_sublead_before_PUCorr, "DY#rightarrow #mu#mu", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_DY_phi_sublead_before_RoccoR, "DY#rightarrow #mu#mu", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_DY_phi_sublead_before_EffCorr, "DY#rightarrow #mu#mu", "f");
        RP_phi_sublead->AddLegendEntry(h_DY_phi_sublead, "DY#rightarrow #mu#mu", "f");


        // Legend MC BKG
        RP_pT->AddLegendEntry(h_bkg_pT[8], "DY#rightarrow #tau#tau", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[8], "DY#rightarrow #tau#tau", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[7], "t#bar{t}", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[7], "t#bar{t}", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[7], "t#bar{t}", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[7], "t#bar{t}", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[7], "t#bar{t}", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[7], "t#bar{t}", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[7], "t#bar{t}", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[7], "t#bar{t}", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[7], "t#bar{t}", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[7], "t#bar{t}", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[7], "t#bar{t}", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[7], "t#bar{t}", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[7], "t#bar{t}", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[7], "t#bar{t}", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[7], "t#bar{t}", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[7], "t#bar{t}", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[7], "t#bar{t}", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[7], "t#bar{t}", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[6], "tW", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[6], "tW", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[6], "tW", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[6], "tW", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[6], "tW", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[6], "tW", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[6], "tW", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[6], "tW", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[6], "tW", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[6], "tW", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[6], "tW", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[6], "tW", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[6], "tW", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[6], "tW", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[6], "tW", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[6], "tW", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[6], "tW", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[6], "tW", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[6], "tW", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[6], "tW", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[6], "tW", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[6], "tW", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[6], "tW", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[6], "tW", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[6], "tW", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[6], "tW", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[5], "#bar{t}W", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[5], "#bar{t}W", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[5], "#bar{t}W", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[5], "#bar{t}W", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[5], "#bar{t}W", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[5], "#bar{t}W", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[5], "#bar{t}W", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[5], "#bar{t}W", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[5], "#bar{t}W", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[5], "#bar{t}W", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[5], "#bar{t}W", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[5], "#bar{t}W", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[5], "#bar{t}W", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[5], "#bar{t}W", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[5], "#bar{t}W", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[5], "#bar{t}W", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[5], "#bar{t}W", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[5], "#bar{t}W", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[4], "ZZ", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[4], "ZZ", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[4], "ZZ", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[4], "ZZ", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[4], "ZZ", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[4], "ZZ", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[4], "ZZ", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[4], "ZZ", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[4], "ZZ", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[4], "ZZ", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[4], "ZZ", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[4], "ZZ", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[4], "ZZ", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[4], "ZZ", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[4], "ZZ", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[4], "ZZ", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[4], "ZZ", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[4], "ZZ", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[4], "ZZ", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[4], "ZZ", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[4], "ZZ", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[4], "ZZ", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[4], "ZZ", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[4], "ZZ", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[4], "ZZ", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[4], "ZZ", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[3], "WZ", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[3], "WZ", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[3], "WZ", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[3], "WZ", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[3], "WZ", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[3], "WZ", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[3], "WZ", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[3], "WZ", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[3], "WZ", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[3], "WZ", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[3], "WZ", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[3], "WZ", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[3], "WZ", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[3], "WZ", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[3], "WZ", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[3], "WZ", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[3], "WZ", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[3], "WZ", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[3], "WZ", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[3], "WZ", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[3], "WZ", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[3], "WZ", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[3], "WZ", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[3], "WZ", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[3], "WZ", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[3], "WZ", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[2], "WW", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[2], "WW", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[2], "WW", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[2], "WW", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[2], "WW", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[2], "WW", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[2], "WW", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[2], "WW", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[2], "WW", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[2], "WW", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[2], "WW", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[2], "WW", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[2], "WW", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[2], "WW", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[2], "WW", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[2], "WW", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[2], "WW", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[2], "WW", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[2], "WW", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[2], "WW", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[2], "WW", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[2], "WW", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[2], "WW", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[2], "WW", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[2], "WW", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[2], "WW", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[1], "W+Jets", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[1], "W+Jets", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[1], "W+Jets", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[1], "W+Jets", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[1], "W+Jets", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[1], "W+Jets", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[1], "W+Jets", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[1], "W+Jets", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[1], "W+Jets", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[1], "W+Jets", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[1], "W+Jets", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[1], "W+Jets", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[1], "W+Jets", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[1], "W+Jets", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[1], "W+Jets", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[1], "W+Jets", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[1], "W+Jets", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[1], "W+Jets", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[1], "W+Jets", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[1], "W+Jets", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[1], "W+Jets", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[1], "W+Jets", "f");

        RP_pT->AddLegendEntry(h_bkg_pT[0], "QCD", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[0], "QCD", "f");
        RP_pT_lead_before_PUCorr->AddLegendEntry(h_bkg_pT_lead_before_PUCorr[0], "QCD", "f");
        RP_pT_lead_before_RoccoR->AddLegendEntry(h_bkg_pT_lead_before_RoccoR[0], "QCD", "f");
        RP_pT_lead_before_EffCorr->AddLegendEntry(h_bkg_pT_lead_before_EffCorr[0], "QCD", "f");
        RP_pT_lead->AddLegendEntry(h_bkg_pT_lead[0], "QCD", "f");
        RP_pT_sublead_before_PUCorr->AddLegendEntry(h_bkg_pT_sublead_before_PUCorr[0], "QCD", "f");
        RP_pT_sublead_before_RoccoR->AddLegendEntry(h_bkg_pT_sublead_before_RoccoR[0], "QCD", "f");
        RP_pT_sublead_before_EffCorr->AddLegendEntry(h_bkg_pT_sublead_before_EffCorr[0], "QCD", "f");
        RP_pT_sublead->AddLegendEntry(h_bkg_pT_sublead[0], "QCD", "f");
        RP_eta_lead_before_PUCorr->AddLegendEntry(h_bkg_eta_lead_before_PUCorr[0], "QCD", "f");
        RP_eta_lead_before_RoccoR->AddLegendEntry(h_bkg_eta_lead_before_RoccoR[0], "QCD", "f");
        RP_eta_lead_before_EffCorr->AddLegendEntry(h_bkg_eta_lead_before_EffCorr[0], "QCD", "f");
        RP_eta_lead->AddLegendEntry(h_bkg_eta_lead[0], "QCD", "f");
        RP_eta_sublead_before_PUCorr->AddLegendEntry(h_bkg_eta_sublead_before_PUCorr[0], "QCD", "f");
        RP_eta_sublead_before_RoccoR->AddLegendEntry(h_bkg_eta_sublead_before_RoccoR[0], "QCD", "f");
        RP_eta_sublead_before_EffCorr->AddLegendEntry(h_bkg_eta_sublead_before_EffCorr[0], "QCD", "f");
        RP_eta_sublead->AddLegendEntry(h_bkg_eta_sublead[0], "QCD", "f");
        RP_phi_lead_before_PUCorr->AddLegendEntry(h_bkg_phi_lead_before_PUCorr[0], "QCD", "f");
        RP_phi_lead_before_RoccoR->AddLegendEntry(h_bkg_phi_lead_before_RoccoR[0], "QCD", "f");
        RP_phi_lead_before_EffCorr->AddLegendEntry(h_bkg_phi_lead_before_EffCorr[0], "QCD", "f");
        RP_phi_lead->AddLegendEntry(h_bkg_phi_lead[0], "QCD", "f");
        RP_phi_sublead_before_PUCorr->AddLegendEntry(h_bkg_phi_sublead_before_PUCorr[0], "QCD", "f");
        RP_phi_sublead_before_RoccoR->AddLegendEntry(h_bkg_phi_sublead_before_RoccoR[0], "QCD", "f");
        RP_phi_sublead_before_EffCorr->AddLegendEntry(h_bkg_phi_sublead_before_EffCorr[0], "QCD", "f");
        RP_phi_sublead->AddLegendEntry(h_bkg_phi_sublead[0], "QCD", "f");

        RP_pT->Draw(0.5, 8e6, 0);
        RP_rapi->Draw(0.5, 8e6, 0);
        RP_pT_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_pT_lead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_pT_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_pT_lead->Draw(0.5, 8e6, 0);
        RP_pT_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_pT_sublead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_pT_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_pT_sublead->Draw(0.5, 8e6, 0);
        RP_eta_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_eta_lead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_eta_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_eta_lead->Draw(0.5, 8e6, 0);
        RP_eta_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_eta_sublead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_eta_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_eta_sublead->Draw(0.5, 8e6, 0);
        RP_phi_lead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_phi_lead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_phi_lead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_phi_lead->Draw(0.5, 8e6, 0);
        RP_phi_sublead_before_PUCorr->Draw(0.5, 8e6, 0);
        RP_phi_sublead_before_RoccoR->Draw(0.5, 8e6, 0);
        RP_phi_sublead_before_EffCorr->Draw(0.5, 8e6, 0);
        RP_phi_sublead->Draw(0.5, 8e6, 0);

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nVTX #################################################

    if( whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[9], *h_bkg_nVTX_before_EffCorr[9], *h_bkg_nVTX[9];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            f_bkg->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter] );

            Color_t color = kBlack;
            if ( pr == _MuMu_QCDMuEnriched_Full ) color = kRed + 3;
            if ( pr == _MuMu_WJets ) color = kRed - 2;
            if ( pr == _MuMu_WW ) color = kMagenta - 5;
            if ( pr == _MuMu_WZ ) color = kMagenta - 2;
            if ( pr == _MuMu_ZZ ) color = kMagenta - 6;
            if ( pr == _MuMu_tbarW ) color = kGreen - 2;
            if ( pr == _MuMu_tW ) color = kGreen + 2;
            if ( pr == _MuMu_ttbar_Full ) color = kCyan + 2;
            if ( pr == _MuMu_DYTauTau_Full ) color = kOrange - 5;

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_nVTX[iter]->SetFillColor(color);

            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_nVTX[iter]->SetLineColor(color);

            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nVTX[iter]->SetDirectory(0);

            s_nVTX_before_PUCorr->Add( h_bkg_nVTX_before_PUCorr[iter] );
            s_nVTX_before_EffCorr->Add( h_bkg_nVTX_before_EffCorr[iter] );
            s_nVTX->Add( h_bkg_nVTX[iter] );           

            iter++;

            if ( pr == _MuMu_QCDMuEnriched_Full )
                pr = _EndOf_MuMu_WJets; // next - WJets
            if ( pr == _MuMu_WJets )
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if ( pr == _MuMu_tW )
                pr = _MuMu_VVnST; // next - ttbar
            if ( pr == _MuMu_DYTauTau_Full ) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nVTX_before_PUCorr, *h_DY_nVTX_before_EffCorr, *h_DY_nVTX;

        f_DY->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX_before_PUCorr );
        f_DY->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX_before_EffCorr );
        f_DY->GetObject( "h_nVTX_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX );

        h_DY_nVTX_before_PUCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetFillColor(kOrange);
        h_DY_nVTX->SetFillColor(kOrange);

        h_DY_nVTX_before_PUCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetLineColor(kOrange);
        h_DY_nVTX->SetLineColor(kOrange);

        h_DY_nVTX_before_PUCorr->SetDirectory(0);
        h_DY_nVTX_before_EffCorr->SetDirectory(0);
        h_DY_nVTX->SetDirectory(0);

        s_nVTX_before_PUCorr->Add( h_DY_nVTX_before_PUCorr );
        s_nVTX_before_EffCorr->Add( h_DY_nVTX_before_EffCorr );
        s_nVTX->Add( h_DY_nVTX );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;

        f_data->GetObject( "h_nVTX_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_nVTX );

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t( "RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX );
        RP_nVTX_before_EffCorr = new myRatioPlot_t( "RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX );
        RP_nVTX = new myRatioPlot_t( "RP_nVTX", s_nVTX, h_data_nVTX );

        RP_nVTX_before_PUCorr->SetPlots("# Primary Vertices (before PU correction)", 0, 50);
        RP_nVTX_before_EffCorr->SetPlots("# Primary Vertices (before Efficiency SF)", 0, 50);
        RP_nVTX->SetPlots("# Primary Vertices", 0, 50);

        RP_nVTX_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nVTX_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nVTX->SetLegend(0.75, 0.4);

        // Legend data
        RP_nVTX_before_PUCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX->AddLegendEntry(h_data_nVTX, "Data", "lp");

        // Legend MC signal
        RP_nVTX_before_PUCorr->AddLegendEntry(h_DY_nVTX_before_PUCorr, "DY#rightarrow#mu#mu", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_DY_nVTX_before_EffCorr, "DY#rightarrow#mu#mu", "f");
        RP_nVTX->AddLegendEntry(h_DY_nVTX, "DY#rightarrow#mu#mu", "f");

        // Legend MC BKG
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[8], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[7], "t#bar{t}", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[7], "t#bar{t}", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[7], "t#bar{t}", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[6], "tW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[6], "tW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[6], "tW", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[5], "#bar{t}W", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[5], "#bar{t}W", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[5], "#bar{t}W", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[4], "ZZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[4], "ZZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[4], "ZZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[3], "WZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[3], "WZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[3], "WZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[2], "WW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[2], "WW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[2], "WW", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[1], "W+Jets", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[1], "W+Jets", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[1], "W+Jets", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[0], "QCD", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[0], "QCD", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[0], "QCD", "f");

        RP_nVTX_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nVTX->Draw(0.5, 3e6, 0);

    } // End of if(nVTX)

} // End of MuMu_HistDrawer()

/// -------------------------------- EMu events ------------------------------------ ///
void EMu_HistDrawer ( TString whichGraphs , TString type)
{
    if ( !whichGraphs.Length() )
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t isWJ = 0;
    isWJ = 1; // UNCOMMENT THIS IF YOU WANT TO INCLUDE W+JETS INTO HISTOGRAMS
    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc( _EMu_Bkg_Full );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _EMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    // Files without Rochester correction
    Mgr.SetProc( _EMu_Bkg_Full );
    TString name_bkg_noRC = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root";
    TFile* f_bkg_noRC = new TFile( name_bkg_noRC, "READ" );
    if (f_bkg_noRC->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _EMu_SingleMuon_Full );
    TString name_data_noRC = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root";
    TFile* f_data_noRC = new TFile( name_data_noRC, "READ" );
    if (f_data_noRC->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=348650, MCintegral;

//################################# INVARIANT MASS #################################################

    if( whichGraphs=="ALL" || whichGraphs=="INVM" || whichGraphs=="INVMASS" )
    {
        count_drawn++;

        THStack *s_mass_fine_before_PUCorr, *s_mass_fine_before_RoccoR, *s_mass_fine_before_EffCorr, *s_mass_fine, *s_mass,
                *s_SS_mass_fine_before_PUCorr, *s_SS_mass_fine_before_RoccoR, *s_SS_mass_fine_before_EffCorr, *s_SS_mass_fine, *s_SS_mass;
        s_mass_fine_before_PUCorr = new THStack("s_mass_fine_before_PUCorr", "");
        s_mass_fine_before_RoccoR = new THStack("s_mass_fine_before_RoccoR", "");
        s_mass_fine_before_EffCorr = new THStack("s_mass_fine_before_EffCorr", "");
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");
        s_SS_mass_fine_before_PUCorr = new THStack("s_SS_mass_fine_before_PUCorr", "");
        s_SS_mass_fine_before_RoccoR = new THStack("s_SS_mass_fine_before_RoccoR", "");
        s_SS_mass_fine_before_EffCorr = new THStack("s_SS_mass_fine_before_EffCorr", "");
        s_SS_mass_fine = new THStack("s_SS_mass_fine", "");
        s_SS_mass = new THStack("s_SS_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine_before_PUCorr[8], *h_bkg_mass_fine_before_EffCorr[8], *h_bkg_mass_fine_before_RoccoR[8],
             *h_bkg_mass_fine[8], *h_bkg_mass[8], *h_SS_bkg_mass_fine_before_PUCorr[8], *h_SS_bkg_mass_fine_before_EffCorr[8],
             *h_SS_bkg_mass_fine_before_RoccoR[8], *h_SS_bkg_mass_fine[8], *h_SS_bkg_mass[8];
        Int_t iter = 0;

        Double_t ratio_WJets_SSvsOS;
        Double_t avg_ratio_WJets_SSvsOS = 0;

        for ( SelProc_t pr = _EMu_WJets; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( !isWJ && pr == _EMu_WJets ) continue;
            if ( pr == _EndOf_EMu_VVnST_Normal ) continue;

            f_bkg_noRC->GetObject( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_emu_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_emu_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );
            f_bkg_noRC->GetObject( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg_noRC->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_"+Mgr.Procname[pr], h_SS_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_emuSS_mass_"+Mgr.Procname[pr], h_SS_bkg_mass[iter] );

            if ( pr == _EMu_WJets )
            {
                ratio_WJets_SSvsOS = h_bkg_mass[iter]->Integral( 1, h_bkg_mass[iter]->GetSize()-2 );
                ratio_WJets_SSvsOS /= h_SS_bkg_mass[iter]->Integral( 1, h_SS_bkg_mass[iter]->GetSize()-2 );

                for ( Int_t i=1; i<h_bkg_mass[iter]->GetSize()-1; i++ )
                {
                    if ( h_SS_bkg_mass[iter]->GetBinContent(i) != 0 )
                        avg_ratio_WJets_SSvsOS += h_bkg_mass[iter]->GetBinContent(i) / h_SS_bkg_mass[iter]->GetBinContent(i);
                }
                avg_ratio_WJets_SSvsOS /= h_SS_bkg_mass[iter]->GetSize()-2;
            }

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_bkg_mass[iter]->SetFillColor(iter+2);
            h_SS_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
            h_SS_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter+2);
            h_SS_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
            h_SS_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_SS_bkg_mass[iter]->SetFillColor(iter+2);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_bkg_mass[iter]->SetLineColor(iter+2);
            h_SS_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
            h_SS_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter+2);
            h_SS_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
            h_SS_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_SS_bkg_mass[iter]->SetLineColor(iter+2);

            h_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);
            h_SS_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
            h_SS_bkg_mass_fine_before_RoccoR[iter]->SetDirectory(0);
            h_SS_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
            h_SS_bkg_mass_fine[iter]->SetDirectory(0);
            h_SS_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
            s_mass_fine_before_RoccoR->Add( h_bkg_mass_fine_before_RoccoR[iter] );
            s_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );
            s_SS_mass_fine_before_PUCorr->Add( h_SS_bkg_mass_fine_before_PUCorr[iter] );
            s_SS_mass_fine_before_RoccoR->Add( h_SS_bkg_mass_fine_before_RoccoR[iter] );
            s_SS_mass_fine_before_EffCorr->Add( h_SS_bkg_mass_fine_before_EffCorr[iter] );
            s_SS_mass_fine->Add( h_SS_bkg_mass_fine[iter] );
            s_SS_mass->Add( h_SS_bkg_mass[iter] );

            iter++;

            if ( pr == _EMu_tW ) // next - ttbar
                pr = _EMu_VVnST;
            if ( pr == _EMu_DYTauTau_Full ) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine_before_PUCorr, *h_data_mass_fine_before_RoccoR, *h_data_mass_fine_before_EffCorr, *h_data_mass_fine, *h_data_mass,
             *h_SS_data_mass_fine_before_PUCorr, *h_SS_data_mass_fine_before_RoccoR, *h_SS_data_mass_fine_before_EffCorr, *h_SS_data_mass_fine, *h_SS_data_mass;

        f_data_noRC->GetObject( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_PUCorr );
        f_data_noRC->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_RoccoR );
        f_data->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_emu_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine );
        f_data->GetObject( "h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass );
        f_data_noRC->GetObject( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine_before_PUCorr );
        f_data_noRC->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine_before_RoccoR );
        f_data->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_emuSS_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine );
        f_data->GetObject( "h_emuSS_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass );

        h_data_mass_fine_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_fine_before_PUCorr->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_fine_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_fine_before_EffCorr->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine_before_PUCorr->SetMarkerColor(kBlack);
        h_data_mass_fine_before_RoccoR->SetMarkerColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetMarkerColor(kBlack);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);
        h_SS_data_mass_fine_before_PUCorr->SetMarkerColor(kBlack);
        h_SS_data_mass_fine_before_RoccoR->SetMarkerColor(kBlack);
        h_SS_data_mass_fine_before_EffCorr->SetMarkerColor(kBlack);
        h_SS_data_mass_fine->SetMarkerColor(kBlack);
        h_SS_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetLineColor(kBlack);
        h_data_mass_fine_before_RoccoR->SetLineColor(kBlack);
        h_data_mass_fine_before_EffCorr->SetLineColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);
        h_SS_data_mass_fine_before_PUCorr->SetLineColor(kBlack);
        h_SS_data_mass_fine_before_RoccoR->SetLineColor(kBlack);
        h_SS_data_mass_fine_before_EffCorr->SetLineColor(kBlack);
        h_SS_data_mass_fine->SetLineColor(kBlack);
        h_SS_data_mass->SetLineColor(kBlack);

        h_data_mass_fine_before_PUCorr->SetDirectory(0);
        h_data_mass_fine_before_RoccoR->SetDirectory(0);
        h_data_mass_fine_before_EffCorr->SetDirectory(0);
        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);
        h_SS_data_mass_fine_before_PUCorr->SetDirectory(0);
        h_SS_data_mass_fine_before_RoccoR->SetDirectory(0);
        h_SS_data_mass_fine_before_EffCorr->SetDirectory(0);
        h_SS_data_mass_fine->SetDirectory(0);
        h_SS_data_mass->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_fine_before_PUCorr, *RP_mass_fine_before_RoccoR, *RP_mass_fine_before_EffCorr, *RP_mass_fine, *RP_mass,
                      *RP_SS_mass_fine_before_PUCorr, *RP_SS_mass_fine_before_RoccoR, *RP_SS_mass_fine_before_EffCorr, *RP_SS_mass_fine, *RP_SS_mass;

        RP_mass_fine_before_PUCorr = new myRatioPlot_t( "RP_mass_fine_before_PUCorr", s_mass_fine_before_PUCorr, h_data_mass_fine_before_PUCorr );
        RP_mass_fine_before_RoccoR = new myRatioPlot_t( "RP_mass_fine_before_RoccoR", s_mass_fine_before_RoccoR, h_data_mass_fine_before_RoccoR );
        RP_mass_fine_before_EffCorr = new myRatioPlot_t( "RP_mass_fine_before_EffCorr", s_mass_fine_before_EffCorr, h_data_mass_fine_before_EffCorr );
        RP_mass_fine = new myRatioPlot_t( "RP_mass_fine", s_mass_fine, h_data_mass_fine );
        RP_mass = new myRatioPlot_t( "RP_mass", s_mass, h_data_mass );
        RP_SS_mass_fine_before_PUCorr = new myRatioPlot_t( "RP_SS_mass_fine_before_PUCorr", s_SS_mass_fine_before_PUCorr, h_SS_data_mass_fine_before_PUCorr );
        RP_SS_mass_fine_before_RoccoR = new myRatioPlot_t( "RP_SS_mass_fine_before_RoccoR", s_SS_mass_fine_before_RoccoR, h_SS_data_mass_fine_before_RoccoR );
        RP_SS_mass_fine_before_EffCorr = new myRatioPlot_t( "RP_SS_mass_fine_before_EffCorr", s_SS_mass_fine_before_EffCorr, h_SS_data_mass_fine_before_EffCorr );
        RP_SS_mass_fine = new myRatioPlot_t( "RP_SS_mass_fine", s_SS_mass_fine, h_SS_data_mass_fine );
        RP_SS_mass = new myRatioPlot_t( "RP_SS_mass", s_SS_mass, h_SS_data_mass );

        RP_mass_fine_before_PUCorr->SetPlots("e#mu mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_mass_fine_before_RoccoR->SetPlots("e#mu mass [GeV/c^{2}] (before Rochester correction)", 15, 3000);
        RP_mass_fine_before_EffCorr->SetPlots("e#mu mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
        RP_mass_fine->SetPlots("e#mu mass [GeV/c^{2}]", 15, 3000);
        RP_mass->SetPlots("e#mu mass [GeV/c^{2}]", 15, 3000);
        RP_SS_mass_fine_before_PUCorr->SetPlots("e#mu (same-signed) mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_SS_mass_fine_before_RoccoR->SetPlots("e#mu (same-signed) mass [GeV/c^{2}] (before Rochester correction)", 15, 3000);
        RP_SS_mass_fine_before_EffCorr->SetPlots("e#mu (same-signed) mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
        RP_SS_mass_fine->SetPlots("e#mu (same-signed) mass [GeV/c^{2}]", 15, 3000);
        RP_SS_mass->SetPlots("e#mu (same-signed) mass [GeV/c^{2}]", 15, 3000);

        RP_mass_fine_before_PUCorr->SetLegend(0.75, 0.4);
        RP_mass_fine_before_RoccoR->SetLegend(0.75, 0.4);
        RP_mass_fine_before_EffCorr->SetLegend(0.75, 0.4);
        RP_mass_fine->SetLegend(0.75, 0.4);
        RP_mass->SetLegend(0.75, 0.4);
        RP_SS_mass_fine_before_PUCorr->SetLegend(0.75, 0.4);
        RP_SS_mass_fine_before_RoccoR->SetLegend(0.75, 0.4);
        RP_SS_mass_fine_before_EffCorr->SetLegend(0.75, 0.4);
        RP_SS_mass_fine->SetLegend(0.75, 0.4);
        RP_SS_mass->SetLegend(0.75, 0.4);

        // Legend data
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_data_mass_fine_before_PUCorr, "Data", "lp");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_data_mass_fine_before_RoccoR, "Data", "lp");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_data_mass_fine_before_EffCorr, "Data", "lp");
        RP_mass_fine->AddLegendEntry(h_data_mass_fine, "Data", "lp");
        RP_mass->AddLegendEntry(h_data_mass, "Data", "lp");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_data_mass_fine_before_PUCorr, "Data", "lp");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_data_mass_fine_before_RoccoR, "Data", "lp");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_data_mass_fine_before_EffCorr, "Data", "lp");
        RP_SS_mass_fine->AddLegendEntry(h_SS_data_mass_fine, "Data", "lp");
        RP_SS_mass->AddLegendEntry(h_SS_data_mass, "Data", "lp");

        // Legend MC BKG        
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[5+isWJ], "t#bar{t}", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[5+isWJ], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[5+isWJ], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[5+isWJ], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[5+isWJ], "t#bar{t}", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[5+isWJ], "t#bar{t}", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[5+isWJ], "t#bar{t}", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[5+isWJ], "t#bar{t}", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[5+isWJ], "t#bar{t}", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[5+isWJ], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4+isWJ], "tW", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[4+isWJ], "tW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4+isWJ], "tW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4+isWJ], "tW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4+isWJ], "tW", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[4+isWJ], "tW", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[4+isWJ], "tW", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[4+isWJ], "tW", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[4+isWJ], "tW", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[4+isWJ], "tW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3+isWJ], "#bar{t}W", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[3+isWJ], "#bar{t}W", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3+isWJ], "#bar{t}W", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3+isWJ], "#bar{t}W", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3+isWJ], "#bar{t}W", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[3+isWJ], "#bar{t}W", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[3+isWJ], "#bar{t}W", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[3+isWJ], "#bar{t}W", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[3+isWJ], "#bar{t}W", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[3+isWJ], "#bar{t}WW", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2+isWJ], "ZZ", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[2+isWJ], "ZZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2+isWJ], "ZZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2+isWJ], "ZZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2+isWJ], "ZZ", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[2+isWJ], "ZZ", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[2+isWJ], "ZZ", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[2+isWJ], "ZZ", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[2+isWJ], "ZZ", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[2+isWJ], "ZZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1+isWJ], "WZ", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[1+isWJ], "WZ", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1+isWJ], "WZ", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1+isWJ], "WZ", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1+isWJ], "WZ", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[1+isWJ], "WZ", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[1+isWJ], "WZ", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[1+isWJ], "WZ", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[1+isWJ], "WZ", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[1+isWJ], "WZ", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0+isWJ], "WW", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[0+isWJ], "WW", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0+isWJ], "WW", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0+isWJ], "WW", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0+isWJ], "WW", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[0+isWJ], "WW", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[0+isWJ], "WW", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[0+isWJ], "WW", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[0+isWJ], "WW", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[0+isWJ], "WW", "f");
        if ( isWJ )
        {
            RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[0], "W+Jets", "f");
            RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[0], "W+Jets", "f");
            RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[0], "W+Jets", "f");
            RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[0], "W+Jets", "f");
            RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[0], "W+Jets", "f");
        }

        RP_mass_fine_before_PUCorr->Draw(0.5, 1e6, 1);
        RP_mass_fine_before_RoccoR->Draw(0.5, 1e6, 1);
        RP_mass_fine_before_EffCorr->Draw(0.5, 1e6, 1);
        RP_mass_fine->Draw(0.5, 1e6, 1);
        RP_mass->Draw(0.5, 1e6, 1);
        RP_SS_mass_fine_before_PUCorr->Draw(0.5, 1e6, 1);
        RP_SS_mass_fine_before_RoccoR->Draw(0.5, 1e6, 1);
        RP_SS_mass_fine_before_EffCorr->Draw(0.5, 1e6, 1);
        RP_SS_mass_fine->Draw(0.5, 1e6, 1);
        RP_SS_mass->Draw(0.5, 1e6, 1);

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;
        if ( isWJ )
        {
            std::cout << "OS/SS ratio of WJets events: " << ratio_WJets_SSvsOS << endl;
            std::cout << "Average OS/SS ratio of WJets events per bin: " << avg_ratio_WJets_SSvsOS << endl;
        }

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################

    if ( whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI") )
    {
        count_drawn++;

        THStack *s_pT_ele, *s_eta_ele, *s_phi_ele, *s_pT_mu, *s_eta_mu, *s_phi_mu,
                *s_pT_eleSS, *s_eta_eleSS, *s_phi_eleSS, *s_pT_muSS, *s_eta_muSS, *s_phi_muSS;

        s_pT_ele = new THStack("s_pT_ele", "");
        s_eta_ele = new THStack("s_eta_ele", "");
        s_phi_ele = new THStack("s_phi_ele", "");
        s_pT_mu = new THStack("s_pT_mu", "");
        s_eta_mu = new THStack("s_eta_mu", "");
        s_phi_mu = new THStack("s_phi_mu", "");
        s_pT_eleSS = new THStack("s_pT_eleSS", "");
        s_eta_eleSS = new THStack("s_eta_eleSS", "");
        s_phi_eleSS = new THStack("s_phi_eleSS", "");
        s_pT_muSS = new THStack("s_pT_muSS", "");
        s_eta_muSS = new THStack("s_eta_muSS", "");
        s_phi_muSS = new THStack("s_phi_muSS", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_pT_ele[8], *h_bkg_eta_ele[8], *h_bkg_phi_ele[8], *h_bkg_pT_mu[8], *h_bkg_eta_mu[8], *h_bkg_phi_mu[8],
             *h_bkg_pT_eleSS[8], *h_bkg_eta_eleSS[8], *h_bkg_phi_eleSS[8], *h_bkg_pT_muSS[8], *h_bkg_eta_muSS[8], *h_bkg_phi_muSS[8];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_WJets; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( !isWJ && pr == _EMu_WJets ) continue;
            if ( pr == _EndOf_EMu_VVnST_Normal ) continue;

            f_bkg->GetObject( "h_ele_pT_"+Mgr.Procname[pr], h_bkg_pT_ele[iter] );
            f_bkg->GetObject( "h_ele_eta_"+Mgr.Procname[pr], h_bkg_eta_ele[iter] );
            f_bkg->GetObject( "h_ele_phi_"+Mgr.Procname[pr], h_bkg_phi_ele[iter] );
            f_bkg->GetObject( "h_mu_pT_"+Mgr.Procname[pr], h_bkg_pT_mu[iter] );
            f_bkg->GetObject( "h_mu_eta_"+Mgr.Procname[pr], h_bkg_eta_mu[iter] );
            f_bkg->GetObject( "h_mu_phi_"+Mgr.Procname[pr], h_bkg_phi_mu[iter] );
            f_bkg->GetObject( "h_eleSS_pT_"+Mgr.Procname[pr], h_bkg_pT_eleSS[iter] );
            f_bkg->GetObject( "h_eleSS_eta_"+Mgr.Procname[pr], h_bkg_eta_eleSS[iter] );
            f_bkg->GetObject( "h_eleSS_phi_"+Mgr.Procname[pr], h_bkg_phi_eleSS[iter] );
            f_bkg->GetObject( "h_muSS_pT_"+Mgr.Procname[pr], h_bkg_pT_muSS[iter] );
            f_bkg->GetObject( "h_muSS_eta_"+Mgr.Procname[pr], h_bkg_eta_muSS[iter] );
            f_bkg->GetObject( "h_muSS_phi_"+Mgr.Procname[pr], h_bkg_phi_muSS[iter] );

            h_bkg_pT_ele[iter]->SetFillColor(iter+2);
            h_bkg_eta_ele[iter]->SetFillColor(iter+2);
            h_bkg_phi_ele[iter]->SetFillColor(iter+2);
            h_bkg_pT_mu[iter]->SetFillColor(iter+2);
            h_bkg_eta_mu[iter]->SetFillColor(iter+2);
            h_bkg_phi_mu[iter]->SetFillColor(iter+2);
            h_bkg_pT_eleSS[iter]->SetFillColor(iter+2);
            h_bkg_eta_eleSS[iter]->SetFillColor(iter+2);
            h_bkg_phi_eleSS[iter]->SetFillColor(iter+2);
            h_bkg_pT_muSS[iter]->SetFillColor(iter+2);
            h_bkg_eta_muSS[iter]->SetFillColor(iter+2);
            h_bkg_phi_muSS[iter]->SetFillColor(iter+2);

            h_bkg_pT_ele[iter]->SetLineColor(iter+2);
            h_bkg_eta_ele[iter]->SetLineColor(iter+2);
            h_bkg_phi_ele[iter]->SetLineColor(iter+2);
            h_bkg_pT_mu[iter]->SetLineColor(iter+2);
            h_bkg_eta_mu[iter]->SetLineColor(iter+2);
            h_bkg_phi_mu[iter]->SetLineColor(iter+2);
            h_bkg_pT_eleSS[iter]->SetLineColor(iter+2);
            h_bkg_eta_eleSS[iter]->SetLineColor(iter+2);
            h_bkg_phi_eleSS[iter]->SetLineColor(iter+2);
            h_bkg_pT_muSS[iter]->SetLineColor(iter+2);
            h_bkg_eta_muSS[iter]->SetLineColor(iter+2);
            h_bkg_phi_muSS[iter]->SetLineColor(iter+2);

            h_bkg_pT_ele[iter]->SetDirectory(0);
            h_bkg_eta_ele[iter]->SetDirectory(0);
            h_bkg_phi_ele[iter]->SetDirectory(0);
            h_bkg_pT_mu[iter]->SetDirectory(0);
            h_bkg_eta_mu[iter]->SetDirectory(0);
            h_bkg_phi_mu[iter]->SetDirectory(0);
            h_bkg_pT_eleSS[iter]->SetDirectory(0);
            h_bkg_eta_eleSS[iter]->SetDirectory(0);
            h_bkg_phi_eleSS[iter]->SetDirectory(0);
            h_bkg_pT_muSS[iter]->SetDirectory(0);
            h_bkg_eta_muSS[iter]->SetDirectory(0);
            h_bkg_phi_muSS[iter]->SetDirectory(0);

            s_pT_ele->Add( h_bkg_pT_ele[iter] );
            s_eta_ele->Add( h_bkg_eta_ele[iter] );
            s_phi_ele->Add( h_bkg_phi_ele[iter] );
            s_pT_mu->Add( h_bkg_pT_mu[iter] );
            s_eta_mu->Add( h_bkg_eta_mu[iter] );
            s_phi_mu->Add( h_bkg_phi_mu[iter] );
            s_pT_eleSS->Add( h_bkg_pT_eleSS[iter] );
            s_eta_eleSS->Add( h_bkg_eta_eleSS[iter] );
            s_phi_eleSS->Add( h_bkg_phi_eleSS[iter] );
            s_pT_muSS->Add( h_bkg_pT_muSS[iter] );
            s_eta_muSS->Add( h_bkg_eta_muSS[iter] );
            s_phi_muSS->Add( h_bkg_phi_muSS[iter] );

            iter++;

            if ( pr == _EMu_tW ) // next - ttbar
                pr = _EMu_VVnST;
            if ( pr == _EMu_DYTauTau_Full ) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT_ele, *h_data_eta_ele, *h_data_phi_ele, *h_data_pT_mu, *h_data_eta_mu, *h_data_phi_mu,
             *h_data_pT_eleSS, *h_data_eta_eleSS, *h_data_phi_eleSS, *h_data_pT_muSS, *h_data_eta_muSS, *h_data_phi_muSS;

        f_data->GetObject( "h_ele_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_ele );
        f_data->GetObject( "h_ele_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_ele );
        f_data->GetObject( "h_ele_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_ele );
        f_data->GetObject( "h_mu_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_mu );
        f_data->GetObject( "h_mu_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_mu );
        f_data->GetObject( "h_mu_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_mu );
        f_data->GetObject( "h_eleSS_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_eleSS );
        f_data->GetObject( "h_eleSS_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_eleSS );
        f_data->GetObject( "h_eleSS_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_eleSS );
        f_data->GetObject( "h_muSS_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_muSS );
        f_data->GetObject( "h_muSS_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_muSS );
        f_data->GetObject( "h_muSS_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_muSS );

        h_data_pT_ele->SetMarkerStyle(kFullDotLarge);
        h_data_eta_ele->SetMarkerStyle(kFullDotLarge);
        h_data_phi_ele->SetMarkerStyle(kFullDotLarge);
        h_data_pT_mu->SetMarkerStyle(kFullDotLarge);
        h_data_eta_mu->SetMarkerStyle(kFullDotLarge);
        h_data_phi_mu->SetMarkerStyle(kFullDotLarge);
        h_data_pT_eleSS->SetMarkerStyle(kFullDotLarge);
        h_data_eta_eleSS->SetMarkerStyle(kFullDotLarge);
        h_data_phi_eleSS->SetMarkerStyle(kFullDotLarge);
        h_data_pT_muSS->SetMarkerStyle(kFullDotLarge);
        h_data_eta_muSS->SetMarkerStyle(kFullDotLarge);
        h_data_phi_muSS->SetMarkerStyle(kFullDotLarge);

        h_data_pT_ele->SetMarkerColor(kBlack);
        h_data_eta_ele->SetMarkerColor(kBlack);
        h_data_phi_ele->SetMarkerColor(kBlack);
        h_data_pT_mu->SetMarkerColor(kBlack);
        h_data_eta_mu->SetMarkerColor(kBlack);
        h_data_phi_mu->SetMarkerColor(kBlack);
        h_data_pT_eleSS->SetMarkerColor(kBlack);
        h_data_eta_eleSS->SetMarkerColor(kBlack);
        h_data_phi_eleSS->SetMarkerColor(kBlack);
        h_data_pT_muSS->SetMarkerColor(kBlack);
        h_data_eta_muSS->SetMarkerColor(kBlack);
        h_data_phi_muSS->SetMarkerColor(kBlack);

        h_data_pT_ele->SetLineColor(kBlack);
        h_data_eta_ele->SetLineColor(kBlack);
        h_data_phi_ele->SetLineColor(kBlack);
        h_data_pT_mu->SetLineColor(kBlack);
        h_data_eta_mu->SetLineColor(kBlack);
        h_data_phi_mu->SetLineColor(kBlack);
        h_data_pT_eleSS->SetLineColor(kBlack);
        h_data_eta_eleSS->SetLineColor(kBlack);
        h_data_phi_eleSS->SetLineColor(kBlack);
        h_data_pT_muSS->SetLineColor(kBlack);
        h_data_eta_muSS->SetLineColor(kBlack);
        h_data_phi_muSS->SetLineColor(kBlack);

        h_data_pT_ele->SetDirectory(0);
        h_data_eta_ele->SetDirectory(0);
        h_data_phi_ele->SetDirectory(0);
        h_data_pT_mu->SetDirectory(0);
        h_data_eta_mu->SetDirectory(0);
        h_data_phi_mu->SetDirectory(0);
        h_data_pT_eleSS->SetDirectory(0);
        h_data_eta_eleSS->SetDirectory(0);
        h_data_phi_eleSS->SetDirectory(0);
        h_data_pT_muSS->SetDirectory(0);
        h_data_eta_muSS->SetDirectory(0);
        h_data_phi_muSS->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_pT_ele, *RP_eta_ele, *RP_phi_ele, *RP_pT_mu, *RP_eta_mu, *RP_phi_mu,
                      *RP_pT_eleSS, *RP_eta_eleSS, *RP_phi_eleSS, *RP_pT_muSS, *RP_eta_muSS, *RP_phi_muSS;

        RP_pT_ele = new myRatioPlot_t( "RP_pT_ele", s_pT_ele, h_data_pT_ele );
        RP_eta_ele = new myRatioPlot_t( "RP_eta_ele", s_eta_ele, h_data_eta_ele );
        RP_phi_ele = new myRatioPlot_t( "RP_phi_ele", s_phi_ele, h_data_phi_ele );
        RP_pT_mu = new myRatioPlot_t( "RP_pT_mu", s_pT_mu, h_data_pT_mu );
        RP_eta_mu = new myRatioPlot_t( "RP_eta_mu", s_eta_mu, h_data_eta_mu );
        RP_phi_mu = new myRatioPlot_t( "RP_phi_mu", s_phi_mu, h_data_phi_mu );
        RP_pT_eleSS = new myRatioPlot_t( "RP_pT_eleSS", s_pT_eleSS, h_data_pT_eleSS );
        RP_eta_eleSS = new myRatioPlot_t( "RP_eta_eleSS", s_eta_eleSS, h_data_eta_eleSS );
        RP_phi_eleSS = new myRatioPlot_t( "RP_phi_eleSS", s_phi_eleSS, h_data_phi_eleSS );
        RP_pT_muSS = new myRatioPlot_t( "RP_pT_muSS", s_pT_muSS, h_data_pT_muSS );
        RP_eta_muSS = new myRatioPlot_t( "RP_eta_muSS", s_eta_muSS, h_data_eta_muSS );
        RP_phi_muSS = new myRatioPlot_t( "RP_phi_muSS", s_phi_muSS, h_data_phi_muSS );

        RP_pT_ele->SetPlots("Electron p_{T} [GeV/c]", 0, 600);
        RP_eta_ele->SetPlots("Electron #eta", -4, 4);
        RP_phi_ele->SetPlots("Electron #phi", -4, 4);
        RP_pT_mu->SetPlots("Muon p_{T} [GeV/c]", 0, 600);
        RP_eta_mu->SetPlots("Muon #eta", -4, 4);
        RP_phi_mu->SetPlots("Muon #phi", -4, 4);
        RP_pT_eleSS->SetPlots("Electron p_{T} (same-sign e#mu) [GeV/c]", 0, 600);
        RP_eta_eleSS->SetPlots("Electron #eta (same-sign e#mu)", -4, 4);
        RP_phi_eleSS->SetPlots("Electron #phi (same-sign e#mu)", -4, 4);
        RP_pT_muSS->SetPlots("Muon p_{T} (same-sign e#mu) [GeV/c]", 0, 600);
        RP_eta_muSS->SetPlots("Muon #eta (same-sign e#mu)", -4, 4);
        RP_phi_muSS->SetPlots("Muon #phi (same-sign e#mu)", -4, 4);

        RP_pT_ele->SetLegend(0.75, 0.4);
        RP_eta_ele->SetLegend(0.75, 0.4);
        RP_phi_ele->SetLegend(0.75, 0.4);
        RP_pT_mu->SetLegend(0.75, 0.4);
        RP_eta_mu->SetLegend(0.75, 0.4);
        RP_phi_mu->SetLegend(0.75, 0.4);
        RP_pT_eleSS->SetLegend(0.75, 0.4);
        RP_eta_eleSS->SetLegend(0.75, 0.4);
        RP_phi_eleSS->SetLegend(0.75, 0.4);
        RP_pT_muSS->SetLegend(0.75, 0.4);
        RP_eta_muSS->SetLegend(0.75, 0.4);
        RP_phi_muSS->SetLegend(0.75, 0.4);

        // Legend data
        RP_pT_ele->AddLegendEntry(h_data_pT_ele, "Data", "lp");
        RP_eta_ele->AddLegendEntry(h_data_eta_ele, "Data", "lp");
        RP_phi_ele->AddLegendEntry(h_data_phi_ele, "Data", "lp");
        RP_pT_mu->AddLegendEntry(h_data_pT_mu, "Data", "lp");
        RP_eta_mu->AddLegendEntry(h_data_eta_mu, "Data", "lp");
        RP_phi_mu->AddLegendEntry(h_data_phi_mu, "Data", "lp");
        RP_pT_eleSS->AddLegendEntry(h_data_pT_eleSS, "Data", "lp");
        RP_eta_eleSS->AddLegendEntry(h_data_eta_eleSS, "Data", "lp");
        RP_phi_eleSS->AddLegendEntry(h_data_phi_eleSS, "Data", "lp");
        RP_pT_muSS->AddLegendEntry(h_data_pT_muSS, "Data", "lp");
        RP_eta_muSS->AddLegendEntry(h_data_eta_muSS, "Data", "lp");
        RP_phi_muSS->AddLegendEntry(h_data_phi_muSS, "Data", "lp");

        // Legend MC BKG
        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[6+isWJ], "DY#rightarrow #tau#tau", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[5+isWJ], "t#bar{t}", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[5+isWJ], "t#bar{t}", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[5+isWJ], "t#bar{t}", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[5+isWJ], "t#bar{t}", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[5+isWJ], "t#bar{t}", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[5+isWJ], "t#bar{t}", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[5+isWJ], "t#bar{t}", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[5+isWJ], "t#bar{t}", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[5+isWJ], "t#bar{t}", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[5+isWJ], "t#bar{t}", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[5+isWJ], "t#bar{t}", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[5+isWJ], "t#bar{t}", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[4+isWJ], "tW", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[4+isWJ], "tW", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[4+isWJ], "tW", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[4+isWJ], "tW", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[4+isWJ], "tW", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[4+isWJ], "tW", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[4+isWJ], "tW", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[4+isWJ], "tW", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[4+isWJ], "tW", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[4+isWJ], "tW", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[4+isWJ], "tW", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[4+isWJ], "tW", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[3+isWJ], "#bar{t}W", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[3+isWJ], "#bar{t}W", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[3+isWJ], "#bar{t}W", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[3+isWJ], "#bar{t}W", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[3+isWJ], "#bar{t}W", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[3+isWJ], "#bar{t}W", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[3+isWJ], "#bar{t}W", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[3+isWJ], "#bar{t}W", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[3+isWJ], "#bar{t}W", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[3+isWJ], "#bar{t}W", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[3+isWJ], "#bar{t}W", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[3+isWJ], "#bar{t}W", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[2+isWJ], "ZZ", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[2+isWJ], "ZZ", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[2+isWJ], "ZZ", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[2+isWJ], "ZZ", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[2+isWJ], "ZZ", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[2+isWJ], "ZZ", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[2+isWJ], "ZZ", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[2+isWJ], "ZZ", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[2+isWJ], "ZZ", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[2+isWJ], "ZZ", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[2+isWJ], "ZZ", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[2+isWJ], "ZZ", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[1+isWJ], "WZ", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[1+isWJ], "WZ", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[1+isWJ], "WZ", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[1+isWJ], "WZ", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[1+isWJ], "WZ", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[1+isWJ], "WZ", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[1+isWJ], "WZ", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[1+isWJ], "WZ", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[1+isWJ], "WZ", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[1+isWJ], "WZ", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[1+isWJ], "WZ", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[1+isWJ], "WZ", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[0+isWJ], "WW", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[0+isWJ], "WW", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[0+isWJ], "WW", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[0+isWJ], "WW", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[0+isWJ], "WW", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[0+isWJ], "WW", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[0+isWJ], "WW", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[0+isWJ], "WW", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[0+isWJ], "WW", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[0+isWJ], "WW", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[0+isWJ], "WW", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[0+isWJ], "WW", "f");

        if ( isWJ )
        {
            RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[0], "W+Jets", "f");
            RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[0], "W+Jets", "f");
            RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[0], "W+Jets", "f");
            RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[0], "W+Jets", "f");
            RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[0], "W+Jets", "f");
            RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[0], "W+Jets", "f");
            RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[0], "W+Jets", "f");
            RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[0], "W+Jets", "f");
            RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[0], "W+Jets", "f");
            RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[0], "W+Jets", "f");
            RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[0], "W+Jets", "f");
            RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[0], "W+Jets", "f");
        }

        RP_pT_ele->Draw(0.5, 1e6, 0);
        RP_eta_ele->Draw(0.5, 1e6, 0);
        RP_phi_ele->Draw(0.5, 1e6, 0);
        RP_pT_mu->Draw(0.5, 1e6, 0);
        RP_eta_mu->Draw(0.5, 1e6, 0);
        RP_phi_mu->Draw(0.5, 1e6, 0);
        RP_pT_eleSS->Draw(0.5, 1e6, 0);
        RP_eta_eleSS->Draw(0.5, 1e6, 0);
        RP_phi_eleSS->Draw(0.5, 1e6, 0);
        RP_pT_muSS->Draw(0.5, 1e6, 0);
        RP_eta_muSS->Draw(0.5, 1e6, 0);
        RP_phi_muSS->Draw(0.5, 1e6, 0);

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nVTX #################################################

    if( whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[8], *h_bkg_nVTX_before_EffCorr[8], *h_bkg_nVTX[8];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_WJets; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( !isWJ && pr == _EMu_WJets ) continue;
            if ( pr == _EndOf_EMu_VVnST_Normal ) continue;

            f_bkg->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter] );

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_nVTX[iter]->SetFillColor(iter+2);

            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_nVTX[iter]->SetLineColor(iter+2);

            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nVTX[iter]->SetDirectory(0);

            s_nVTX_before_PUCorr->Add( h_bkg_nVTX_before_PUCorr[iter] );
            s_nVTX_before_EffCorr->Add( h_bkg_nVTX_before_EffCorr[iter] );
            s_nVTX->Add( h_bkg_nVTX[iter] );

            iter++;

            if ( pr == _EMu_tW ) // next - ttbar
                pr = _EMu_VVnST;
            if ( pr == _EMu_DYTauTau_Full ) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;

        f_data->GetObject( "h_nVTX_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_nVTX );

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t( "RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX );
        RP_nVTX_before_EffCorr = new myRatioPlot_t( "RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX );
        RP_nVTX = new myRatioPlot_t( "RP_nVTX", s_nVTX, h_data_nVTX );

        RP_nVTX_before_PUCorr->SetPlots("# Primary Vertices (before PU correction)", 0, 50);
        RP_nVTX_before_EffCorr->SetPlots("# Primary Vertices (before Efficiency correction)", 0, 50);
        RP_nVTX->SetPlots("# Primary Vertices", 0, 50);

        RP_nVTX_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nVTX_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nVTX->SetLegend(0.75, 0.4);

        // Legend data
        RP_nVTX_before_PUCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_data_nVTX, "Data", "lp");
        RP_nVTX->AddLegendEntry(h_data_nVTX, "Data", "lp");

        // Legend MC BKG
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[6+isWJ], "DY#rightarrow #tau#tau", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[5+isWJ], "t#bar{t}", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[5+isWJ], "t#bar{t}", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[5+isWJ], "t#bar{t}", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[4+isWJ], "tW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[4+isWJ], "tW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[4+isWJ], "tW", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[3+isWJ], "#bar{t}W", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[3+isWJ], "#bar{t}W", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[3+isWJ], "#bar{t}W", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[2+isWJ], "ZZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[2+isWJ], "ZZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[2+isWJ], "ZZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[1+isWJ], "WZ", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[1+isWJ], "WZ", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[1+isWJ], "WZ", "f");
        RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[0+isWJ], "WW", "f");
        RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[0+isWJ], "WW", "f");
        RP_nVTX->AddLegendEntry(h_bkg_nVTX[0+isWJ], "WW", "f");
        if ( isWJ )
        {
            RP_nVTX_before_PUCorr->AddLegendEntry(h_bkg_nVTX_before_PUCorr[0], "W+Jets", "f");
            RP_nVTX_before_EffCorr->AddLegendEntry(h_bkg_nVTX_before_EffCorr[0], "W+Jets", "f");
            RP_nVTX->AddLegendEntry(h_bkg_nVTX[0], "W+Jets", "f");
        }



        RP_nVTX_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nVTX->Draw(0.5, 3e6, 0);

    } // End of if(nVTX)

} // End of EMu_HistDrawer()
