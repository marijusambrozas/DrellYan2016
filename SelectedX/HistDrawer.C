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


void HistDrawer ( TString whichX = "", TString whichGraphs = "ALL", TString type = "" )
{
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        cout << "\n*******      EE_HistDrawer ( " << whichGraphs << " " << type << " )      *******" << endl;
        EE_HistDrawer( whichGraphs, type );
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        cout << "\n*****  MuMu_HistDrawer ( " << whichGraphs << " " << type << " )  *****" << endl;
        MuMu_HistDrawer( whichGraphs, type );
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
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

    Mgr.GetProc( _EE_DY_Full );
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.GetProc( _EE_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.GetProc( _EE_DoubleEG_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DoubleEG_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_DoubleEG_Full]+".root" << " opened successfully" << endl;

//################################# INVARIANT MASS #################################################

    if(whichGraphs=="all" || whichGraphs=="All" || whichGraphs=="ALL" || whichGraphs=="invm" || whichGraphs=="InvM" ||
       whichGraphs=="Invm" || whichGraphs=="invM" || whichGraphs=="INVM" || whichGraphs=="InvMass"||
       whichGraphs=="invmass" || whichGraphs=="invMass" || whichGraphs=="INVMASS")
    {
        count_drawn++;

        THStack *s_mass_fine_before_PUCorr, *s_mass_fine_before_EffCorr, *s_mass_fine, *s_mass;
        s_mass_fine_before_PUCorr = new THStack("s_mass_fine_before_PUCorr", "");
        s_mass_fine_before_EffCorr = new THStack("s_mass_fine_before_EffCorr", "");
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine_before_PUCorr[5], *h_bkg_mass_fine_before_EffCorr[5],
             *h_bkg_mass_fine[5], *h_bkg_mass[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EE_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_bkg_mass[iter]->SetFillColor(iter+2);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_bkg_mass[iter]->SetLineColor(iter+2);

            h_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
            s_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine_before_PUCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine_before_EffCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_mass_"+Mgr.Procname[_EE_WJets], h_bkg_mass[iter] );

                h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
                h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
                h_bkg_mass_fine[iter]->SetFillColor(iter+2);
                h_bkg_mass[iter]->SetFillColor(iter+2);

                h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
                h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
                h_bkg_mass_fine[iter]->SetLineColor(iter+2);
                h_bkg_mass[iter]->SetLineColor(iter+2);

                h_bkg_mass_fine_before_PUCorr[iter]->SetDirectory(0);
                h_bkg_mass_fine_before_EffCorr[iter]->SetDirectory(0);
                h_bkg_mass_fine[iter]->SetDirectory(0);
                h_bkg_mass[iter]->SetDirectory(0);

                s_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
                s_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
                s_mass_fine->Add( h_bkg_mass_fine[iter] );
                s_mass->Add( h_bkg_mass[iter] );

            } // End of WJets

            iter++;            
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

        f_data->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass );

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
        RP_mass_fine_before_EffCorr->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
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
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "QCD", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "QCD", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "QCD", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "QCD", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "VVnST", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "VVnST", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "VVnST", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "VVnST", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2], "W+Jets", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2], "W+Jets", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2], "W+Jets", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2], "W+Jets", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4], "DY#rightarrow #tau#tau", "f");

        RP_mass_fine_before_PUCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_EffCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine->Draw(0.5, 3e6, 1);
        RP_mass->Draw(0.5, 3e6, 1);
        Double_t dataerror, MCerror;
        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(invm)

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

    Mgr.GetProc( _EE_DY_Full );
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    if (f_DY->IsOpen()) std::cout << "File " << name_DY << " opened successfully" << endl;
    Mgr.GetProc( _EE_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << name_bkg << " opened successfully" << endl;
    Mgr.GetProc( _EE_DoubleEG_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DoubleEG_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << name_data << " opened successfully" << endl;

//################################# INVARIANT MASS #################################################

    if(whichGraphs=="all" || whichGraphs=="All" || whichGraphs=="ALL" || whichGraphs=="invm" || whichGraphs=="InvM" ||
       whichGraphs=="Invm" || whichGraphs=="invM" || whichGraphs=="INVM" || whichGraphs=="InvMass"||
       whichGraphs=="invmass" || whichGraphs=="invMass" || whichGraphs=="INVMASS")
    {
        count_drawn++;

        THStack *s_mass_fine_before_PUCorr, *s_mass_fine_before_RoccoR, *s_mass_fine_before_EffCorr, *s_mass_fine, *s_mass;
        s_mass_fine_before_PUCorr = new THStack("s_mass_fine_before_PUCorr", "");
        s_mass_fine_before_RoccoR = new THStack("s_mass_fine_before_RoccoR", "");
        s_mass_fine_before_EffCorr = new THStack("s_mass_fine_before_EffCorr", "");
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine_before_PUCorr[5], *h_bkg_mass_fine_before_RoccoR[5], *h_bkg_mass_fine_before_EffCorr[5],
             *h_bkg_mass_fine[5], *h_bkg_mass[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EE_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[pr], h_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter);
            h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter);
            h_bkg_mass_fine[iter]->SetFillColor(iter);
            h_bkg_mass[iter]->SetFillColor(iter);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter);
            h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter);
            h_bkg_mass_fine[iter]->SetLineColor(iter);
            h_bkg_mass[iter]->SetLineColor(iter);

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

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine_before_PUCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine_before_RoccoR[iter] );
                f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine_before_EffCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_mass_"+Mgr.Procname[_EE_WJets], h_bkg_mass[iter] );

                h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter);
                h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter);
                h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter);
                h_bkg_mass_fine[iter]->SetFillColor(iter);
                h_bkg_mass[iter]->SetFillColor(iter);

                h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter);
                h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter);
                h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter);
                h_bkg_mass_fine[iter]->SetLineColor(iter);
                h_bkg_mass[iter]->SetLineColor(iter);

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
            } // End of WJets

            iter++;
        } // End of for(bkg)



//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_fine_before_PUCorr, *h_DY_mass_fine_before_RoccoR, *h_DY_mass_fine_before_EffCorr, *h_DY_mass_fine, *h_DY_mass;

        f_DY->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine_before_PUCorr );
        f_DY->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine_before_RoccoR );
        f_DY->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine_before_EffCorr );
        f_DY->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine );
        f_DY->GetObject( "h_mass_"+Mgr.Procname[_EE_DY_Full], h_DY_mass );

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

        f_data->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine_before_RoccoR );
        f_data->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[_EE_DoubleEG_Full], h_data_mass );

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

        RP_mass_fine_before_PUCorr->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_mass_fine_before_RoccoR->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before Rochester correction)", 15, 3000);
        RP_mass_fine_before_EffCorr->SetPlots("Dielectron invariant mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
        RP_mass_fine->SetPlots("Dielectron invariant mass [GeV/c^{2}]", 15, 3000);
        RP_mass->SetPlots("Dielectron invariant mass [GeV/c^{2}]", 15, 3000);

        RP_mass_fine_before_PUCorr->SetLegend(0.75, 0.4);
        RP_mass_fine_before_RoccoR->SetLegend(0.75, 0.4);
        RP_mass_fine_before_EffCorr->SetLegend(0.75, 0.4);
        RP_mass_fine->SetLegend(0.75, 0.4);
        RP_mass->SetLegend(0.75, 0.4);

        // Legend data
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_data_mass_fine_before_PUCorr, "Data", "lp");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_data_mass_fine_before_RoccoR, "Data","lp");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_data_mass_fine_before_EffCorr, "Data", "lp");
        RP_mass_fine->AddLegendEntry(h_data_mass_fine, "Data", "lp");
        RP_mass->AddLegendEntry(h_data_mass, "Data", "lp");

        // Legend MC signal
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_DY_mass_fine_before_PUCorr, "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_DY_mass_fine_before_RoccoR, "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_DY_mass_fine_before_EffCorr, "Data", "f");
        RP_mass_fine->AddLegendEntry(h_DY_mass_fine, "Data", "f");
        RP_mass->AddLegendEntry(h_DY_mass, "Data", "f");

        // Legend MC BKG
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[0], "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "Data", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "Data", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "Data", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[1], "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "Data", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "Data", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "Data", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[iter], "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[iter], "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[iter], "Data", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[iter], "Data", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[iter], "Data", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[3], "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "Data", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "Data", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "Data", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4], "Data", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[4], "Data","f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4], "Data", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4], "Data", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4], "Data", "f");

        RP_mass_fine_before_PUCorr->Draw(10e-2, 4e5, 1);
        RP_mass_fine_before_RoccoR->Draw(10e-2, 4e5, 1);
        RP_mass_fine_before_EffCorr->Draw(10e-2, 4e5, 1);
        RP_mass_fine->Draw(10e-2, 4e5, 1);
        RP_mass->Draw(10e-2, 4e5, 1);
        Double_t dataerror, MCerror;
        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(invm)

} // End of MuMu_HistDrawer()


/// --------------------------------- EMu events --------------------------------- ///
void EMu_HistDrawer ( TString whichGraphs, TString type )
{

    return;
} // End of EMu_HistDrawer()
