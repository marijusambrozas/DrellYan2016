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
    TString name_PileUp = "./etc/PileUp/80X/ROOTFile_PUReWeight_80X_v20170817_64mb.root";
    TFile *f_PileUp = new TFile( name_PileUp, "READ" );
    if (f_PileUp->IsOpen()) std::cout << "File " << "ROOTFile_PUReWeight_80X_v20170817_64mb.root" << " opened successfully" << endl;


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

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################

    if ( whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI") )
    {
        count_drawn++;

        THStack *s_Pt, *s_rapi, *s_pT, *s_eta, *s_phi;
        s_Pt = new THStack("s_Pt", "");
        s_rapi = new THStack("s_rapi", "");
        s_pT = new THStack("s_pT", "");
        s_eta = new THStack("s_eta", "");
        s_phi = new THStack("s_phi", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_Pt[5], *h_bkg_rapi[5], *h_bkg_pT[5], *h_bkg_eta[5], *h_bkg_phi[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EE_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_Pt_"+Mgr.Procname[pr], h_bkg_Pt[iter] );
            f_bkg->GetObject( "h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter] );
            f_bkg->GetObject( "h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter] );
            f_bkg->GetObject( "h_eta_"+Mgr.Procname[pr], h_bkg_eta[iter] );
            f_bkg->GetObject( "h_phi_"+Mgr.Procname[pr], h_bkg_phi[iter] );

            h_bkg_Pt[iter]->SetFillColor(iter+2);
            h_bkg_rapi[iter]->SetFillColor(iter+2);
            h_bkg_pT[iter]->SetFillColor(iter+2);
            h_bkg_eta[iter]->SetFillColor(iter+2);
            h_bkg_phi[iter]->SetFillColor(iter+2);

            h_bkg_Pt[iter]->SetLineColor(iter+2);
            h_bkg_rapi[iter]->SetLineColor(iter+2);
            h_bkg_pT[iter]->SetLineColor(iter+2);
            h_bkg_eta[iter]->SetLineColor(iter+2);
            h_bkg_phi[iter]->SetLineColor(iter+2);


            h_bkg_Pt[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_eta[iter]->SetDirectory(0);
            h_bkg_phi[iter]->SetDirectory(0);


            s_Pt->Add( h_bkg_Pt[iter] );
            s_rapi->Add( h_bkg_rapi[iter] );
            s_pT->Add( h_bkg_pT[iter] );
            s_eta->Add( h_bkg_eta[iter] );
            s_phi->Add( h_bkg_phi[iter] );


            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_Pt_"+Mgr.Procname[_EE_WJets], h_bkg_Pt[iter] );
                f_bkg->GetObject( "h_rapi_"+Mgr.Procname[_EE_WJets], h_bkg_rapi[iter] );
                f_bkg->GetObject( "h_pT_"+Mgr.Procname[_EE_WJets], h_bkg_pT[iter] );
                f_bkg->GetObject( "h_eta_"+Mgr.Procname[_EE_WJets], h_bkg_eta[iter] );
                f_bkg->GetObject( "h_phi_"+Mgr.Procname[_EE_WJets], h_bkg_phi[iter] );


                h_bkg_Pt[iter]->SetFillColor(iter+2);
                h_bkg_rapi[iter]->SetFillColor(iter+2);
                h_bkg_pT[iter]->SetFillColor(iter+2);
                h_bkg_eta[iter]->SetFillColor(iter+2);
                h_bkg_phi[iter]->SetFillColor(iter+2);


                h_bkg_Pt[iter]->SetLineColor(iter+2);
                h_bkg_rapi[iter]->SetLineColor(iter+2);
                h_bkg_pT[iter]->SetLineColor(iter+2);
                h_bkg_eta[iter]->SetLineColor(iter+2);
                h_bkg_phi[iter]->SetLineColor(iter+2);


                h_bkg_Pt[iter]->SetDirectory(0);
                h_bkg_rapi[iter]->SetDirectory(0);
                h_bkg_pT[iter]->SetDirectory(0);
                h_bkg_eta[iter]->SetDirectory(0);
                h_bkg_phi[iter]->SetDirectory(0);


                s_Pt->Add( h_bkg_Pt[iter] );
                s_rapi->Add( h_bkg_rapi[iter] );
                s_pT->Add( h_bkg_pT[iter] );
                s_eta->Add( h_bkg_eta[iter] );
                s_phi->Add( h_bkg_phi[iter] );
            } // End of WJets

            iter++;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_Pt, *h_DY_rapi, *h_DY_pT, *h_DY_eta, *h_DY_phi;

        f_DY->GetObject( "h_Pt_"+Mgr.Procname[_EE_DY_Full], h_DY_Pt );
        f_DY->GetObject( "h_rapi_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi );
        f_DY->GetObject( "h_pT_"+Mgr.Procname[_EE_DY_Full], h_DY_pT );
        f_DY->GetObject( "h_eta_"+Mgr.Procname[_EE_DY_Full], h_DY_eta );
        f_DY->GetObject( "h_phi_"+Mgr.Procname[_EE_DY_Full], h_DY_phi );


        h_DY_Pt->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT->SetFillColor(kOrange);
        h_DY_eta->SetFillColor(kOrange);
        h_DY_phi->SetFillColor(kOrange);


        h_DY_Pt->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT->SetLineColor(kOrange);
        h_DY_eta->SetLineColor(kOrange);
        h_DY_phi->SetLineColor(kOrange);


        h_DY_Pt->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT->SetDirectory(0);
        h_DY_eta->SetDirectory(0);
        h_DY_phi->SetDirectory(0);


        s_Pt->Add( h_DY_Pt );
        s_rapi->Add( h_DY_rapi );
        s_pT->Add( h_DY_pT );
        s_eta->Add( h_DY_eta );
        s_phi->Add( h_DY_phi );


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_Pt, *h_data_rapi, *h_data_pT, *h_data_eta, *h_data_phi;

        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc], h_data_Pt );
        f_data->GetObject( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], h_data_rapi );
        f_data->GetObject( "h_pT_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT );
        f_data->GetObject( "h_eta_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta );
        f_data->GetObject( "h_phi_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi );


        h_data_Pt->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_eta->SetMarkerStyle(kFullDotLarge);
        h_data_phi->SetMarkerStyle(kFullDotLarge);


        h_data_Pt->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);
        h_data_pT->SetMarkerColor(kBlack);
        h_data_eta->SetMarkerColor(kBlack);
        h_data_phi->SetMarkerColor(kBlack);


        h_data_Pt->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);
        h_data_pT->SetLineColor(kBlack);
        h_data_eta->SetLineColor(kBlack);
        h_data_phi->SetLineColor(kBlack);


        h_data_Pt->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT->SetDirectory(0);
        h_data_eta->SetDirectory(0);
        h_data_phi->SetDirectory(0);


//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_Pt, *RP_rapi, *RP_pT, *RP_eta, *RP_phi;
        RP_Pt = new myRatioPlot_t( "RP_Pt", s_Pt, h_data_Pt );
        RP_rapi = new myRatioPlot_t( "RP_rapi", s_rapi, h_data_rapi );
        RP_pT = new myRatioPlot_t( "RP_pT", s_pT, h_data_pT );
        RP_eta = new myRatioPlot_t( "RP_eta", s_eta, h_data_eta );
        RP_phi = new myRatioPlot_t( "RP_phi", s_phi, h_data_phi );


        RP_Pt->SetPlots("Dielectron p_{T} [GeV/c]", 0, 600);
        RP_rapi->SetPlots("Dielectron rapidity", -3, 3);
        RP_pT->SetPlots("Electron p_{T} (leading+subleading) [GeV/c]", 0, 600);
        RP_eta->SetPlots("Electron #eta (leading+subleading)", -3, 3);
        RP_phi->SetPlots("Electron #phi (leading+subleading)", -3.3, 3.3);


        RP_Pt->SetLegend(0.75, 0.4);
        RP_rapi->SetLegend(0.75, 0.4);
        RP_pT->SetLegend(0.75, 0.4);
        RP_eta->SetLegend(0.75, 0.4);
        RP_phi->SetLegend(0.75, 0.4);


        // Legend data
        RP_Pt->AddLegendEntry(h_data_Pt, "Data", "lp");
        RP_rapi->AddLegendEntry(h_data_rapi, "Data", "lp");
        RP_pT->AddLegendEntry(h_data_pT, "Data", "lp");
        RP_eta->AddLegendEntry(h_data_eta, "Data", "lp");
        RP_phi->AddLegendEntry(h_data_phi, "Data", "lp");


        // Legend MC signal
        RP_Pt->AddLegendEntry(h_DY_Pt, "DY#rightarrow ee", "f");
        RP_rapi->AddLegendEntry(h_DY_rapi, "DY#rightarrow ee", "f");
        RP_pT->AddLegendEntry(h_DY_pT, "DY#rightarrow ee", "f");
        RP_eta->AddLegendEntry(h_DY_eta, "DY#rightarrow ee", "f");
        RP_phi->AddLegendEntry(h_DY_phi, "DY#rightarrow ee", "f");


        // Legend MC BKG
        RP_Pt->AddLegendEntry(h_bkg_Pt[0], "QCD", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[0], "QCD", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[0], "QCD", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[0], "QCD", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[0], "QCD", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[1], "VVnST", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[1], "VVnST", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[1], "VVnST", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[1], "VVnST", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[1], "VVnST", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[2], "W+Jets", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[2], "W+Jets", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[2], "W+Jets", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[2], "W+Jets", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[2], "W+Jets", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[3], "t#bar{t}", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[3], "t#bar{t}", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[3], "t#bar{t}", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[3], "t#bar{t}", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[3], "t#bar{t}", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[4], "DY#rightarrow #tau#tau", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[4], "DY#rightarrow #tau#tau", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[4], "DY#rightarrow #tau#tau", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[4], "DY#rightarrow #tau#tau", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[4], "DY#rightarrow #tau#tau", "f");

        RP_Pt->Draw(0.5, 3e6, 0);
        RP_rapi->Draw(0.5, 3e6, 0);
        RP_pT->Draw(0.5, 3e6, 0);
        RP_eta->Draw(0.5, 3e6, 0);
        RP_phi->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nPU #################################################

    if( whichGraphs=="ALL" || whichGraphs=="NPU" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nPU_before_PUCorr, *s_nPU_before_EffCorr, *s_nPU;
        s_nPU_before_PUCorr = new THStack("s_nPU_before_PUCorr", "");
        s_nPU_before_EffCorr = new THStack("s_nPU_before_EffCorr", "");
        s_nPU = new THStack("s_nPU", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nPU_before_PUCorr[5], *h_bkg_nPU_before_EffCorr[5], *h_bkg_nPU[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EE_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nPU_"+Mgr.Procname[pr], h_bkg_nPU[iter] );

            h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU[iter]->SetFillColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU[iter]->SetLineColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nPU[iter]->SetDirectory(0);

            s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
            s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
            s_nPU->Add( h_bkg_nPU[iter] );

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[_EE_WJets], h_bkg_nPU_before_PUCorr[iter] );
                f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[_EE_WJets], h_bkg_nPU_before_EffCorr[iter] );
                f_bkg->GetObject( "h_nPU_"+Mgr.Procname[_EE_WJets], h_bkg_nPU[iter] );

                h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU[iter]->SetFillColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU[iter]->SetLineColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
                h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
                h_bkg_nPU[iter]->SetDirectory(0);

                s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
                s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
                s_nPU->Add( h_bkg_nPU[iter] );

            } // End of WJets

            iter++;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nPU_before_PUCorr, *h_DY_nPU_before_EffCorr, *h_DY_nPU;

        f_DY->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nPU_before_PUCorr );
        f_DY->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nPU_before_EffCorr );
        f_DY->GetObject( "h_nPU_"+Mgr.Procname[_EE_DY_Full], h_DY_nPU );

        h_DY_nPU_before_PUCorr->SetFillColor(kOrange);
        h_DY_nPU_before_EffCorr->SetFillColor(kOrange);
        h_DY_nPU->SetFillColor(kOrange);

        h_DY_nPU_before_PUCorr->SetLineColor(kOrange);
        h_DY_nPU_before_EffCorr->SetLineColor(kOrange);
        h_DY_nPU->SetLineColor(kOrange);

        h_DY_nPU_before_PUCorr->SetDirectory(0);
        h_DY_nPU_before_EffCorr->SetDirectory(0);
        h_DY_nPU->SetDirectory(0);

        s_nPU_before_PUCorr->Add( h_DY_nPU_before_PUCorr );
        s_nPU_before_EffCorr->Add( h_DY_nPU_before_EffCorr );
        s_nPU->Add( h_DY_nPU );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nPU;
        TH1D *h_ones = new TH1D( "ones", "", 75, 0, 75 );
        for ( int i=0; i<=75; i++ )
        {
            h_ones->Fill(i);
        }

        f_PileUp->GetObject( "pileup", h_data_nPU );

        h_data_nPU->Multiply( h_ones, h_data_nPU, 1, dataintegral );

        h_data_nPU->SetMarkerStyle(kFullDotLarge);
        h_data_nPU->SetMarkerColor(kBlack);
        h_data_nPU->SetLineColor(kBlack);

        h_data_nPU->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nPU_before_PUCorr, *RP_nPU_before_EffCorr, *RP_nPU;
        RP_nPU_before_PUCorr = new myRatioPlot_t( "RP_nPU_before_PUCorr", s_nPU_before_PUCorr, h_data_nPU );
        RP_nPU_before_EffCorr = new myRatioPlot_t( "RP_nPU_before_EffCorr", s_nPU_before_EffCorr, h_data_nPU );
        RP_nPU = new myRatioPlot_t( "RP_nPU", s_nPU, h_data_nPU );

        RP_nPU_before_PUCorr->SetPlots("# Primary Vertices (before PU correction)", 0, 50);
        RP_nPU_before_EffCorr->SetPlots("# Primary Vertices (before Efficiency correction)", 0, 50);
        RP_nPU->SetPlots("# Primary Vertices", 0, 50);

        RP_nPU_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nPU_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nPU->SetLegend(0.75, 0.4);

        // Legend data
        RP_nPU_before_PUCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU_before_EffCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU->AddLegendEntry(h_data_nPU, "Data", "lp");

        // Legend MC signal
        RP_nPU_before_PUCorr->AddLegendEntry(h_DY_nPU_before_PUCorr, "DY#rightarrow ee", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_DY_nPU_before_EffCorr, "DY#rightarrow ee", "f");
        RP_nPU->AddLegendEntry(h_DY_nPU, "DY#rightarrow ee", "f");

        // Legend MC BKG
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[0], "QCD", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[0], "QCD", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[0], "QCD", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[1], "VVnST", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[1], "VVnST", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[1], "VVnST", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[2], "W+Jets", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[2], "W+Jets", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[2], "W+Jets", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[3], "t#bar{t}", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[3], "t#bar{t}", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[3], "t#bar{t}", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[4], "DY#rightarrow #tau#tau", "f");

        RP_nPU_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nPU_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nPU->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(nPU)

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
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _MuMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;
    TString name_PileUp = "./etc/PileUp/80X/ROOTFile_PUReWeight_80X_v20170817_64mb.root";
    TFile *f_PileUp = new TFile( name_PileUp, "READ" );
    if (f_PileUp->IsOpen()) std::cout << "File " << "ROOTFile_PUReWeight_80X_v20170817_64mb.root" << " opened successfully" << endl;


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
        TH1D *h_bkg_mass_fine_before_PUCorr[5], *h_bkg_mass_fine_before_EffCorr[5], *h_bkg_mass_fine_before_RoccoR[5],
             *h_bkg_mass_fine[5], *h_bkg_mass[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _MuMu_DY_Full; pr=SelProc_t((int)(pr-1)) )
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

            h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_bkg_mass[iter]->SetFillColor(iter+2);

            h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_bkg_mass[iter]->SetLineColor(iter+2);

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
                f_bkg->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass_fine_before_PUCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass_fine_before_RoccoR[iter] );
                f_bkg->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass_fine_before_EffCorr[iter] );
                f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_mass_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass[iter] );

                h_bkg_mass_fine_before_PUCorr[iter]->SetFillColor(iter+2);
                h_bkg_mass_fine_before_RoccoR[iter]->SetFillColor(iter+2);
                h_bkg_mass_fine_before_EffCorr[iter]->SetFillColor(iter+2);
                h_bkg_mass_fine[iter]->SetFillColor(iter+2);
                h_bkg_mass[iter]->SetFillColor(iter+2);

                h_bkg_mass_fine_before_PUCorr[iter]->SetLineColor(iter+2);
                h_bkg_mass_fine_before_RoccoR[iter]->SetLineColor(iter+2);
                h_bkg_mass_fine_before_EffCorr[iter]->SetLineColor(iter+2);
                h_bkg_mass_fine[iter]->SetLineColor(iter+2);
                h_bkg_mass[iter]->SetLineColor(iter+2);

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

        f_DY->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine_before_PUCorr );
        f_DY->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine_before_RoccoR );
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

        f_data->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine_before_RoccoR );
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
        RP_mass_fine_before_EffCorr->SetPlots("Dimuon invariant mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
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
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "QCD", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[0], "QCD", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "QCD", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "QCD", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "QCD", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "VVnST", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[1], "VVnST", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "VVnST", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "VVnST", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "VVnST", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2], "W+Jets", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[2], "W+Jets", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2], "W+Jets", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2], "W+Jets", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2], "W+Jets", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "t#bar{t}", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[3], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[4], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[4], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[4], "DY#rightarrow #tau#tau", "f");

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

        THStack *s_Pt, *s_rapi, *s_pT, *s_eta, *s_phi;
        s_Pt = new THStack("s_Pt", "");
        s_rapi = new THStack("s_rapi", "");
        s_pT = new THStack("s_pT", "");
        s_eta = new THStack("s_eta", "");
        s_phi = new THStack("s_phi", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_Pt[5], *h_bkg_rapi[5], *h_bkg_pT[5], *h_bkg_eta[5], *h_bkg_phi[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _MuMu_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_Pt_"+Mgr.Procname[pr], h_bkg_Pt[iter] );
            f_bkg->GetObject( "h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter] );
            f_bkg->GetObject( "h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter] );
            f_bkg->GetObject( "h_eta_"+Mgr.Procname[pr], h_bkg_eta[iter] );
            f_bkg->GetObject( "h_phi_"+Mgr.Procname[pr], h_bkg_phi[iter] );

            h_bkg_Pt[iter]->SetFillColor(iter+2);
            h_bkg_rapi[iter]->SetFillColor(iter+2);
            h_bkg_pT[iter]->SetFillColor(iter+2);
            h_bkg_eta[iter]->SetFillColor(iter+2);
            h_bkg_phi[iter]->SetFillColor(iter+2);

            h_bkg_Pt[iter]->SetLineColor(iter+2);
            h_bkg_rapi[iter]->SetLineColor(iter+2);
            h_bkg_pT[iter]->SetLineColor(iter+2);
            h_bkg_eta[iter]->SetLineColor(iter+2);
            h_bkg_phi[iter]->SetLineColor(iter+2);


            h_bkg_Pt[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_eta[iter]->SetDirectory(0);
            h_bkg_phi[iter]->SetDirectory(0);


            s_Pt->Add( h_bkg_Pt[iter] );
            s_rapi->Add( h_bkg_rapi[iter] );
            s_pT->Add( h_bkg_pT[iter] );
            s_eta->Add( h_bkg_eta[iter] );
            s_phi->Add( h_bkg_phi[iter] );


            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_Pt_"+Mgr.Procname[_MuMu_WJets], h_bkg_Pt[iter] );
                f_bkg->GetObject( "h_rapi_"+Mgr.Procname[_MuMu_WJets], h_bkg_rapi[iter] );
                f_bkg->GetObject( "h_pT_"+Mgr.Procname[_MuMu_WJets], h_bkg_pT[iter] );
                f_bkg->GetObject( "h_eta_"+Mgr.Procname[_MuMu_WJets], h_bkg_eta[iter] );
                f_bkg->GetObject( "h_phi_"+Mgr.Procname[_MuMu_WJets], h_bkg_phi[iter] );


                h_bkg_Pt[iter]->SetFillColor(iter+2);
                h_bkg_rapi[iter]->SetFillColor(iter+2);
                h_bkg_pT[iter]->SetFillColor(iter+2);
                h_bkg_eta[iter]->SetFillColor(iter+2);
                h_bkg_phi[iter]->SetFillColor(iter+2);


                h_bkg_Pt[iter]->SetLineColor(iter+2);
                h_bkg_rapi[iter]->SetLineColor(iter+2);
                h_bkg_pT[iter]->SetLineColor(iter+2);
                h_bkg_eta[iter]->SetLineColor(iter+2);
                h_bkg_phi[iter]->SetLineColor(iter+2);


                h_bkg_Pt[iter]->SetDirectory(0);
                h_bkg_rapi[iter]->SetDirectory(0);
                h_bkg_pT[iter]->SetDirectory(0);
                h_bkg_eta[iter]->SetDirectory(0);
                h_bkg_phi[iter]->SetDirectory(0);


                s_Pt->Add( h_bkg_Pt[iter] );
                s_rapi->Add( h_bkg_rapi[iter] );
                s_pT->Add( h_bkg_pT[iter] );
                s_eta->Add( h_bkg_eta[iter] );
                s_phi->Add( h_bkg_phi[iter] );
            } // End of WJets

            iter++;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_Pt, *h_DY_rapi, *h_DY_pT, *h_DY_eta, *h_DY_phi;

        f_DY->GetObject( "h_Pt_"+Mgr.Procname[_MuMu_DY_Full], h_DY_Pt );
        f_DY->GetObject( "h_rapi_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi );
        f_DY->GetObject( "h_pT_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT );
        f_DY->GetObject( "h_eta_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta );
        f_DY->GetObject( "h_phi_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi );


        h_DY_Pt->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT->SetFillColor(kOrange);
        h_DY_eta->SetFillColor(kOrange);
        h_DY_phi->SetFillColor(kOrange);


        h_DY_Pt->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT->SetLineColor(kOrange);
        h_DY_eta->SetLineColor(kOrange);
        h_DY_phi->SetLineColor(kOrange);


        h_DY_Pt->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT->SetDirectory(0);
        h_DY_eta->SetDirectory(0);
        h_DY_phi->SetDirectory(0);


        s_Pt->Add( h_DY_Pt );
        s_rapi->Add( h_DY_rapi );
        s_pT->Add( h_DY_pT );
        s_eta->Add( h_DY_eta );
        s_phi->Add( h_DY_phi );


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_Pt, *h_data_rapi, *h_data_pT, *h_data_eta, *h_data_phi;

        f_data->GetObject( "h_Pt_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_Pt );
        f_data->GetObject( "h_rapi_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_rapi );
        f_data->GetObject( "h_pT_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT );
        f_data->GetObject( "h_eta_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta );
        f_data->GetObject( "h_phi_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi );


        h_data_Pt->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_eta->SetMarkerStyle(kFullDotLarge);
        h_data_phi->SetMarkerStyle(kFullDotLarge);


        h_data_Pt->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);
        h_data_pT->SetMarkerColor(kBlack);
        h_data_eta->SetMarkerColor(kBlack);
        h_data_phi->SetMarkerColor(kBlack);


        h_data_Pt->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);
        h_data_pT->SetLineColor(kBlack);
        h_data_eta->SetLineColor(kBlack);
        h_data_phi->SetLineColor(kBlack);


        h_data_Pt->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT->SetDirectory(0);
        h_data_eta->SetDirectory(0);
        h_data_phi->SetDirectory(0);


//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_Pt, *RP_rapi, *RP_pT, *RP_eta, *RP_phi;
        RP_Pt = new myRatioPlot_t( "RP_Pt", s_Pt, h_data_Pt );
        RP_rapi = new myRatioPlot_t( "RP_rapi", s_rapi, h_data_rapi );
        RP_pT = new myRatioPlot_t( "RP_pT", s_pT, h_data_pT );
        RP_eta = new myRatioPlot_t( "RP_eta", s_eta, h_data_eta );
        RP_phi = new myRatioPlot_t( "RP_phi", s_phi, h_data_phi );


        RP_Pt->SetPlots("Dimuon p_{T} [GeV/c]", 0, 600);
        RP_rapi->SetPlots("Dimuon rapidity", -3, 3);
        RP_pT->SetPlots("Muon p_{T} (leading+subleading) [GeV/c]", 0, 600);
        RP_eta->SetPlots("Muon #eta (leading+subleading)", -3, 3);
        RP_phi->SetPlots("Muon #phi (leading+subleading)", -3.3, 3.3);


        RP_Pt->SetLegend(0.75, 0.4);
        RP_rapi->SetLegend(0.75, 0.4);
        RP_pT->SetLegend(0.75, 0.4);
        RP_eta->SetLegend(0.75, 0.4);
        RP_phi->SetLegend(0.75, 0.4);


        // Legend data
        RP_Pt->AddLegendEntry(h_data_Pt, "Data", "lp");
        RP_rapi->AddLegendEntry(h_data_rapi, "Data", "lp");
        RP_pT->AddLegendEntry(h_data_pT, "Data", "lp");
        RP_eta->AddLegendEntry(h_data_eta, "Data", "lp");
        RP_phi->AddLegendEntry(h_data_phi, "Data", "lp");


        // Legend MC signal
        RP_Pt->AddLegendEntry(h_DY_Pt, "DY#rightarrow#mu#mu", "f");
        RP_rapi->AddLegendEntry(h_DY_rapi, "DY#rightarrow#mu#mu", "f");
        RP_pT->AddLegendEntry(h_DY_pT, "DY#rightarrow#mu#mu", "f");
        RP_eta->AddLegendEntry(h_DY_eta, "DY#rightarrow#mu#mu", "f");
        RP_phi->AddLegendEntry(h_DY_phi, "DY#rightarrow#mu#mu", "f");


        // Legend MC BKG
        RP_Pt->AddLegendEntry(h_bkg_Pt[0], "QCD", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[0], "QCD", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[0], "QCD", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[0], "QCD", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[0], "QCD", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[1], "VVnST", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[1], "VVnST", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[1], "VVnST", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[1], "VVnST", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[1], "VVnST", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[2], "W+Jets", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[2], "W+Jets", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[2], "W+Jets", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[2], "W+Jets", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[2], "W+Jets", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[3], "t#bar{t}", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[3], "t#bar{t}", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[3], "t#bar{t}", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[3], "t#bar{t}", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[3], "t#bar{t}", "f");
        RP_Pt->AddLegendEntry(h_bkg_Pt[4], "DY#rightarrow #tau#tau", "f");
        RP_rapi->AddLegendEntry(h_bkg_rapi[4], "DY#rightarrow #tau#tau", "f");
        RP_pT->AddLegendEntry(h_bkg_pT[4], "DY#rightarrow #tau#tau", "f");
        RP_eta->AddLegendEntry(h_bkg_eta[4], "DY#rightarrow #tau#tau", "f");
        RP_phi->AddLegendEntry(h_bkg_phi[4], "DY#rightarrow #tau#tau", "f");

        RP_Pt->Draw(0.5, 3e6, 0);
        RP_rapi->Draw(0.5, 3e6, 0);
        RP_pT->Draw(0.5, 3e6, 0);
        RP_eta->Draw(0.5, 3e6, 0);
        RP_phi->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nPU #################################################

    if( whichGraphs=="ALL" || whichGraphs=="NPU" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nPU_before_PUCorr, *s_nPU_before_EffCorr, *s_nPU;
        s_nPU_before_PUCorr = new THStack("s_nPU_before_PUCorr", "");
        s_nPU_before_EffCorr = new THStack("s_nPU_before_EffCorr", "");
        s_nPU = new THStack("s_nPU", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nPU_before_PUCorr[5], *h_bkg_nPU_before_EffCorr[5], *h_bkg_nPU[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _MuMu_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nPU_"+Mgr.Procname[pr], h_bkg_nPU[iter] );

            h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU[iter]->SetFillColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU[iter]->SetLineColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nPU[iter]->SetDirectory(0);

            s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
            s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
            s_nPU->Add( h_bkg_nPU[iter] );

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[_MuMu_WJets], h_bkg_nPU_before_PUCorr[iter] );
                f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[_MuMu_WJets], h_bkg_nPU_before_EffCorr[iter] );
                f_bkg->GetObject( "h_nPU_"+Mgr.Procname[_MuMu_WJets], h_bkg_nPU[iter] );

                h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU[iter]->SetFillColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU[iter]->SetLineColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
                h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
                h_bkg_nPU[iter]->SetDirectory(0);

                s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
                s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
                s_nPU->Add( h_bkg_nPU[iter] );

            } // End of WJets

            iter++;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nPU_before_PUCorr, *h_DY_nPU_before_EffCorr, *h_DY_nPU;

        f_DY->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nPU_before_PUCorr );
        f_DY->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nPU_before_EffCorr );
        f_DY->GetObject( "h_nPU_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nPU );

        h_DY_nPU_before_PUCorr->SetFillColor(kOrange);
        h_DY_nPU_before_EffCorr->SetFillColor(kOrange);
        h_DY_nPU->SetFillColor(kOrange);

        h_DY_nPU_before_PUCorr->SetLineColor(kOrange);
        h_DY_nPU_before_EffCorr->SetLineColor(kOrange);
        h_DY_nPU->SetLineColor(kOrange);

        h_DY_nPU_before_PUCorr->SetDirectory(0);
        h_DY_nPU_before_EffCorr->SetDirectory(0);
        h_DY_nPU->SetDirectory(0);

        s_nPU_before_PUCorr->Add( h_DY_nPU_before_PUCorr );
        s_nPU_before_EffCorr->Add( h_DY_nPU_before_EffCorr );
        s_nPU->Add( h_DY_nPU );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nPU;
        TH1D *h_ones = new TH1D( "ones", "", 75, 0, 75 );
        for ( int i=0; i<=75; i++ )
        {
            h_ones->Fill(i);
        }

        f_PileUp->GetObject( "pileup", h_data_nPU );

        h_data_nPU->Multiply( h_ones, h_data_nPU, 1, dataintegral );

        h_data_nPU->SetMarkerStyle(kFullDotLarge);
        h_data_nPU->SetMarkerColor(kBlack);
        h_data_nPU->SetLineColor(kBlack);

        h_data_nPU->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nPU_before_PUCorr, *RP_nPU_before_EffCorr, *RP_nPU;
        RP_nPU_before_PUCorr = new myRatioPlot_t( "RP_nPU_before_PUCorr", s_nPU_before_PUCorr, h_data_nPU );
        RP_nPU_before_EffCorr = new myRatioPlot_t( "RP_nPU_before_EffCorr", s_nPU_before_EffCorr, h_data_nPU );
        RP_nPU = new myRatioPlot_t( "RP_nPU", s_nPU, h_data_nPU );

        RP_nPU_before_PUCorr->SetPlots("# Primary Vertices (before PU correction)", 0, 50);
        RP_nPU_before_EffCorr->SetPlots("# Primary Vertices (before Efficiency correction)", 0, 50);
        RP_nPU->SetPlots("# Primary Vertices", 0, 50);

        RP_nPU_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nPU_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nPU->SetLegend(0.75, 0.4);

        // Legend data
        RP_nPU_before_PUCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU_before_EffCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU->AddLegendEntry(h_data_nPU, "Data", "lp");

        // Legend MC signal
        RP_nPU_before_PUCorr->AddLegendEntry(h_DY_nPU_before_PUCorr, "DY#rightarrow#mu#mu", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_DY_nPU_before_EffCorr, "DY#rightarrow#mu#mu", "f");
        RP_nPU->AddLegendEntry(h_DY_nPU, "DY#rightarrow#mu#mu", "f");

        // Legend MC BKG
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[0], "QCD", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[0], "QCD", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[0], "QCD", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[1], "VVnST", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[1], "VVnST", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[1], "VVnST", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[2], "W+Jets", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[2], "W+Jets", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[2], "W+Jets", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[3], "t#bar{t}", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[3], "t#bar{t}", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[3], "t#bar{t}", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[4], "DY#rightarrow #tau#tau", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[4], "DY#rightarrow #tau#tau", "f");

        RP_nPU_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nPU_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nPU->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(nPU)

} // End of MuMu_HistDrawer()

/// -------------------------------- EMu events ------------------------------------ ///
void EMu_HistDrawer ( TString whichGraphs , TString type)
{
    if ( !whichGraphs.Length() )
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc( _EMu_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc( _EMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root" << " opened successfully" << endl;
    TString name_PileUp = "./etc/PileUp/80X/ROOTFile_PUReWeight_80X_v20170817_64mb.root";
    TFile *f_PileUp = new TFile( name_PileUp, "READ" );
    if (f_PileUp->IsOpen()) std::cout << "File " << "ROOTFile_PUReWeight_80X_v20170817_64mb.root" << " opened successfully" << endl;


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
        TH1D *h_bkg_mass_fine_before_PUCorr[4], *h_bkg_mass_fine_before_EffCorr[4], *h_bkg_mass_fine_before_RoccoR[4],
             *h_bkg_mass_fine[5], *h_bkg_mass[4], *h_SS_bkg_mass_fine_before_PUCorr[4], *h_SS_bkg_mass_fine_before_EffCorr[4],
             *h_SS_bkg_mass_fine_before_RoccoR[4], *h_SS_bkg_mass_fine[4], *h_SS_bkg_mass[4];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_VVnST; pr > _EndOf_EMu_Data_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 3 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg->GetObject( "h_emu_mass_fine_before_RoccoR_"+Mgr.Procname[pr], h_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_emu_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_emu_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_PUCorr[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_before_RoccoR_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_RoccoR[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_fine_before_EffCorr[iter] );
            f_bkg->GetObject( "h_emuSS_mass_fine_"+Mgr.Procname[pr], h_SS_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_emuSS_mass_"+Mgr.Procname[pr], h_SS_bkg_mass[iter] );

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

            if ( iter == 0 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_WJets], h_bkg_mass_fine_before_PUCorr[iter] );
                f_bkg->GetObject( "h_emu_mass_fine_before_RoccoR_"+Mgr.Procname[_EMu_WJets], h_bkg_mass_fine_before_RoccoR[iter] );
                f_bkg->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_WJets], h_bkg_mass_fine_before_EffCorr[iter] );
                f_bkg->GetObject( "h_emu_mass_fine_"+Mgr.Procname[_EMu_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_emu_mass_"+Mgr.Procname[_EMu_WJets], h_bkg_mass[iter] );
                f_bkg->GetObject( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_WJets], h_SS_bkg_mass_fine_before_PUCorr[iter] );
                f_bkg->GetObject( "h_emuSS_mass_fine_before_RoccoR_"+Mgr.Procname[_EMu_WJets], h_SS_bkg_mass_fine_before_RoccoR[iter] );
                f_bkg->GetObject( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_WJets], h_SS_bkg_mass_fine_before_EffCorr[iter] );
                f_bkg->GetObject( "h_emuSS_mass_fine_"+Mgr.Procname[_EMu_WJets], h_SS_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_emuSS_mass_"+Mgr.Procname[_EMu_WJets], h_SS_bkg_mass[iter] );

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
                s_SS_mass_fine_before_PUCorr->Add( h_bkg_mass_fine_before_PUCorr[iter] );
                s_SS_mass_fine_before_RoccoR->Add( h_bkg_mass_fine_before_RoccoR[iter] );
                s_SS_mass_fine_before_EffCorr->Add( h_bkg_mass_fine_before_EffCorr[iter] );
                s_SS_mass_fine->Add( h_bkg_mass_fine[iter] );
                s_SS_mass->Add( h_bkg_mass[iter] );

            } // End of WJets

            iter++;
        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine_before_PUCorr, *h_data_mass_fine_before_RoccoR, *h_data_mass_fine_before_EffCorr, *h_data_mass_fine, *h_data_mass,
             *h_SS_data_mass_fine_before_PUCorr, *h_SS_data_mass_fine_before_RoccoR, *h_SS_data_mass_fine_before_EffCorr, *h_SS_data_mass_fine, *h_SS_data_mass;

        f_data->GetObject( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_emu_mass_fine_before_RoccoR_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_RoccoR );
        f_data->GetObject( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine_before_EffCorr );
        f_data->GetObject( "h_emu_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine );
        f_data->GetObject( "h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass );
        f_data->GetObject( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine_before_PUCorr );
        f_data->GetObject( "h_emuSS_mass_fine_before_RoccoR_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine_before_RoccoR );
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
        RP_SS_mass_fine_before_PUCorr->SetPlots("e#mu (single-signed) mass [GeV/c^{2}] (before PU correction)", 15, 3000);
        RP_SS_mass_fine_before_RoccoR->SetPlots("e#mu (single-signed) mass [GeV/c^{2}] (before Rochester correction)", 15, 3000);
        RP_SS_mass_fine_before_EffCorr->SetPlots("e#mu (single-signed) mass [GeV/c^{2}] (before Efficiency correction)", 15, 3000);
        RP_SS_mass_fine->SetPlots("e#mu (single-signed) mass [GeV/c^{2}]", 15, 3000);
        RP_SS_mass->SetPlots("e#mu (single-signed) mass [GeV/c^{2}]", 15, 3000);

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
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[0], "VVnST", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[0], "VVnST", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[0], "VVnST", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[0], "VVnST", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[0], "VVnST", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[0], "VVnST", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[0], "VVnST", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[0], "VVnST", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[0], "VVnST", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[0], "VVnST", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[1], "W+Jets", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[1], "W+Jets", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[1], "W+Jets", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[1], "W+Jets", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[1], "W+Jets", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[1], "W+Jets", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[1], "W+Jets", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[1], "W+Jets", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[1], "W+Jets", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[1], "W+Jets", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[2], "t#bar{t}", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[2], "t#bar{t}", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[2], "t#bar{t}", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[2], "t#bar{t}", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[2], "t#bar{t}", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[2], "t#bar{t}", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[2], "t#bar{t}", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[2], "t#bar{t}", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[2], "t#bar{t}", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[2], "t#bar{t}", "f");
        RP_mass_fine_before_PUCorr->AddLegendEntry(h_bkg_mass_fine_before_PUCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_RoccoR->AddLegendEntry(h_bkg_mass_fine_before_RoccoR[3], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine_before_EffCorr->AddLegendEntry(h_bkg_mass_fine_before_EffCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_mass_fine->AddLegendEntry(h_bkg_mass_fine[3], "DY#rightarrow #tau#tau", "f");
        RP_mass->AddLegendEntry(h_bkg_mass[3], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_PUCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_PUCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_RoccoR->AddLegendEntry(h_SS_bkg_mass_fine_before_RoccoR[3], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine_before_EffCorr->AddLegendEntry(h_SS_bkg_mass_fine_before_EffCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass_fine->AddLegendEntry(h_SS_bkg_mass_fine[3], "DY#rightarrow #tau#tau", "f");
        RP_SS_mass->AddLegendEntry(h_SS_bkg_mass[3], "DY#rightarrow #tau#tau", "f");

        RP_mass_fine_before_PUCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_RoccoR->Draw(0.5, 3e6, 1);
        RP_mass_fine_before_EffCorr->Draw(0.5, 3e6, 1);
        RP_mass_fine->Draw(0.5, 3e6, 1);
        RP_mass->Draw(0.5, 1e7, 1);
        RP_SS_mass_fine_before_PUCorr->Draw(0.5, 3e6, 1);
        RP_SS_mass_fine_before_RoccoR->Draw(0.5, 3e6, 1);
        RP_SS_mass_fine_before_EffCorr->Draw(0.5, 3e6, 1);
        RP_SS_mass_fine->Draw(0.5, 3e6, 1);
        RP_SS_mass->Draw(0.5, 1e7, 1);

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

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
        TH1D *h_bkg_pT_ele[4], *h_bkg_eta_ele[4], *h_bkg_phi_ele[4], *h_bkg_pT_mu[4], *h_bkg_eta_mu[4], *h_bkg_phi_mu[4],
             *h_bkg_pT_eleSS[4], *h_bkg_eta_eleSS[4], *h_bkg_phi_eleSS[4], *h_bkg_pT_muSS[4], *h_bkg_eta_muSS[4], *h_bkg_phi_muSS[4];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_VVnST; pr > _EndOf_EMu_Data_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 3 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

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

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_ele_pT_"+Mgr.Procname[_EMu_WJets], h_bkg_pT_ele[iter] );
                f_bkg->GetObject( "h_ele_eta_"+Mgr.Procname[_EMu_WJets], h_bkg_eta_ele[iter] );
                f_bkg->GetObject( "h_ele_phi_"+Mgr.Procname[_EMu_WJets], h_bkg_phi_ele[iter] );
                f_bkg->GetObject( "h_mu_pT_"+Mgr.Procname[_EMu_WJets], h_bkg_pT_mu[iter] );
                f_bkg->GetObject( "h_mu_eta_"+Mgr.Procname[_EMu_WJets], h_bkg_eta_mu[iter] );
                f_bkg->GetObject( "h_mu_phi_"+Mgr.Procname[_EMu_WJets], h_bkg_phi_mu[iter] );
                f_bkg->GetObject( "h_eleSS_pT_"+Mgr.Procname[_EMu_WJets], h_bkg_pT_eleSS[iter] );
                f_bkg->GetObject( "h_eleSS_eta_"+Mgr.Procname[_EMu_WJets], h_bkg_eta_eleSS[iter] );
                f_bkg->GetObject( "h_eleSS_phi_"+Mgr.Procname[_EMu_WJets], h_bkg_phi_eleSS[iter] );
                f_bkg->GetObject( "h_muSS_pT_"+Mgr.Procname[_EMu_WJets], h_bkg_pT_muSS[iter] );
                f_bkg->GetObject( "h_muSS_eta_"+Mgr.Procname[_EMu_WJets], h_bkg_eta_muSS[iter] );
                f_bkg->GetObject( "h_muSS_phi_"+Mgr.Procname[_EMu_WJets], h_bkg_phi_muSS[iter] );

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
            } // End of WJets

            iter++;
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
        RP_eta_ele->SetPlots("Electron #eta", -3, 3);
        RP_phi_ele->SetPlots("Electron #phi", -3.3, 3.3);
        RP_pT_mu->SetPlots("Muon p_{T} [GeV/c]", 0, 600);
        RP_eta_mu->SetPlots("Muon #eta", -3, 3);
        RP_phi_mu->SetPlots("Muon #phi", -3.3, 3.3);
        RP_pT_eleSS->SetPlots("Electron p_{T} (single-sign e#mu) [GeV/c]", 0, 600);
        RP_eta_eleSS->SetPlots("Electron #eta (single-sign e#mu)", -3, 3);
        RP_phi_eleSS->SetPlots("Electron #phi (single-sign e#mu)", -3.3, 3.3);
        RP_pT_muSS->SetPlots("Muon p_{T} (single-sign e#mu) [GeV/c]", 0, 600);
        RP_eta_muSS->SetPlots("Muon #eta (single-sign e#mu)", -3, 3);
        RP_phi_muSS->SetPlots("Muon #phi (single-sign e#mu)", -3.3, 3.3);

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
        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[0], "VVnST", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[0], "VVnST", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[0], "VVnST", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[0], "VVnST", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[0], "VVnST", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[0], "VVnST", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[0], "VVnST", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[0], "VVnST", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[0], "VVnST", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[0], "VVnST", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[0], "VVnST", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[0], "VVnST", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[1], "W+Jets", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[1], "W+Jets", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[1], "W+Jets", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[1], "W+Jets", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[1], "W+Jets", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[1], "W+Jets", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[1], "W+Jets", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[1], "W+Jets", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[1], "W+Jets", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[1], "W+Jets", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[1], "W+Jets", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[1], "W+Jets", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[2], "t#bar{t}", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[2], "t#bar{t}", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[2], "t#bar{t}", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[2], "t#bar{t}", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[2], "t#bar{t}", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[2], "t#bar{t}", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[2], "t#bar{t}", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[2], "t#bar{t}", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[2], "t#bar{t}", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[2], "t#bar{t}", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[2], "t#bar{t}", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[2], "t#bar{t}", "f");

        RP_pT_ele->AddLegendEntry(h_bkg_pT_ele[3], "DY#rightarrow #tau#tau", "f");
        RP_eta_ele->AddLegendEntry(h_bkg_eta_ele[3], "DY#rightarrow #tau#tau", "f");
        RP_phi_ele->AddLegendEntry(h_bkg_phi_ele[3], "DY#rightarrow #tau#tau", "f");
        RP_pT_mu->AddLegendEntry(h_bkg_pT_mu[3], "DY#rightarrow #tau#tau", "f");
        RP_eta_mu->AddLegendEntry(h_bkg_eta_mu[3], "DY#rightarrow #tau#tau", "f");
        RP_phi_mu->AddLegendEntry(h_bkg_phi_mu[3], "DY#rightarrow #tau#tau", "f");
        RP_pT_eleSS->AddLegendEntry(h_bkg_pT_eleSS[3], "DY#rightarrow #tau#tau", "f");
        RP_eta_eleSS->AddLegendEntry(h_bkg_eta_eleSS[3], "DY#rightarrow #tau#tau", "f");
        RP_phi_eleSS->AddLegendEntry(h_bkg_phi_eleSS[3], "DY#rightarrow #tau#tau", "f");
        RP_pT_muSS->AddLegendEntry(h_bkg_pT_muSS[3], "DY#rightarrow #tau#tau", "f");
        RP_eta_muSS->AddLegendEntry(h_bkg_eta_muSS[3], "DY#rightarrow #tau#tau", "f");
        RP_phi_muSS->AddLegendEntry(h_bkg_phi_muSS[3], "DY#rightarrow #tau#tau", "f");

        RP_pT_ele->Draw(0.5, 3e6, 0);
        RP_eta_ele->Draw(0.5, 3e6, 0);
        RP_phi_ele->Draw(0.5, 3e6, 0);
        RP_pT_mu->Draw(0.5, 3e6, 0);
        RP_eta_mu->Draw(0.5, 3e6, 0);
        RP_phi_mu->Draw(0.5, 3e6, 0);
        RP_pT_eleSS->Draw(0.5, 3e6, 0);
        RP_eta_eleSS->Draw(0.5, 3e6, 0);
        RP_phi_eleSS->Draw(0.5, 3e6, 0);
        RP_pT_muSS->Draw(0.5, 3e6, 0);
        RP_eta_muSS->Draw(0.5, 3e6, 0);
        RP_phi_muSS->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nPU #################################################

    if( whichGraphs=="ALL" || whichGraphs=="NPU" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP") )
    {
        count_drawn++;

        THStack *s_nPU_before_PUCorr, *s_nPU_before_EffCorr, *s_nPU;
        s_nPU_before_PUCorr = new THStack("s_nPU_before_PUCorr", "");
        s_nPU_before_EffCorr = new THStack("s_nPU_before_EffCorr", "");
        s_nPU = new THStack("s_nPU", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nPU_before_PUCorr[4], *h_bkg_nPU_before_EffCorr[4], *h_bkg_nPU[4];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_VVnST; pr > _EndOf_EMu_Data_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 3 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_PUCorr[iter] );
            f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nPU_before_EffCorr[iter] );
            f_bkg->GetObject( "h_nPU_"+Mgr.Procname[pr], h_bkg_nPU[iter] );

            h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
            h_bkg_nPU[iter]->SetFillColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
            h_bkg_nPU[iter]->SetLineColor(iter+2);

            h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_nPU[iter]->SetDirectory(0);

            s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
            s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
            s_nPU->Add( h_bkg_nPU[iter] );

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[_EMu_WJets], h_bkg_nPU_before_PUCorr[iter] );
                f_bkg->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[_EMu_WJets], h_bkg_nPU_before_EffCorr[iter] );
                f_bkg->GetObject( "h_nPU_"+Mgr.Procname[_EMu_WJets], h_bkg_nPU[iter] );

                h_bkg_nPU_before_PUCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetFillColor(iter+2);
                h_bkg_nPU[iter]->SetFillColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU_before_EffCorr[iter]->SetLineColor(iter+2);
                h_bkg_nPU[iter]->SetLineColor(iter+2);

                h_bkg_nPU_before_PUCorr[iter]->SetDirectory(0);
                h_bkg_nPU_before_EffCorr[iter]->SetDirectory(0);
                h_bkg_nPU[iter]->SetDirectory(0);

                s_nPU_before_PUCorr->Add( h_bkg_nPU_before_PUCorr[iter] );
                s_nPU_before_EffCorr->Add( h_bkg_nPU_before_EffCorr[iter] );
                s_nPU->Add( h_bkg_nPU[iter] );

            } // End of WJets

            iter++;
        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nPU;
        TH1D *h_ones = new TH1D( "ones", "", 75, 0, 75 );
        for ( int i=0; i<=75; i++ )
        {
            h_ones->Fill(i);
        }

        f_PileUp->GetObject( "pileup", h_data_nPU );

        h_data_nPU->Multiply( h_ones, h_data_nPU, 1, dataintegral );

        h_data_nPU->SetMarkerStyle(kFullDotLarge);
        h_data_nPU->SetMarkerColor(kBlack);
        h_data_nPU->SetLineColor(kBlack);

        h_data_nPU->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nPU_before_PUCorr, *RP_nPU_before_EffCorr, *RP_nPU;
        RP_nPU_before_PUCorr = new myRatioPlot_t( "RP_nPU_before_PUCorr", s_nPU_before_PUCorr, h_data_nPU );
        RP_nPU_before_EffCorr = new myRatioPlot_t( "RP_nPU_before_EffCorr", s_nPU_before_EffCorr, h_data_nPU );
        RP_nPU = new myRatioPlot_t( "RP_nPU", s_nPU, h_data_nPU );

        RP_nPU_before_PUCorr->SetPlots("# Primary Vertices (before PU correction)", 0, 50);
        RP_nPU_before_EffCorr->SetPlots("# Primary Vertices (before Efficiency correction)", 0, 50);
        RP_nPU->SetPlots("# Primary Vertices", 0, 50);

        RP_nPU_before_PUCorr->SetLegend(0.75, 0.4);
        RP_nPU_before_EffCorr->SetLegend(0.75, 0.4);
        RP_nPU->SetLegend(0.75, 0.4);

        // Legend data
        RP_nPU_before_PUCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU_before_EffCorr->AddLegendEntry(h_data_nPU, "Data", "lp");
        RP_nPU->AddLegendEntry(h_data_nPU, "Data", "lp");

        // Legend MC BKG
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[0], "VVnST", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[0], "VVnST", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[0], "VVnST", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[1], "W+Jets", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[1], "W+Jets", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[1], "W+Jets", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[2], "t#bar{t}", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[2], "t#bar{t}", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[2], "t#bar{t}", "f");
        RP_nPU_before_PUCorr->AddLegendEntry(h_bkg_nPU_before_PUCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_nPU_before_EffCorr->AddLegendEntry(h_bkg_nPU_before_EffCorr[3], "DY#rightarrow #tau#tau", "f");
        RP_nPU->AddLegendEntry(h_bkg_nPU[3], "DY#rightarrow #tau#tau", "f");

        RP_nPU_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nPU_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nPU->Draw(0.5, 3e6, 0);
//        Double_t dataerror, MCerror;
//        Double_t dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
//        Double_t MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

//        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
//        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

    } // End of if(nPU)

} // End of EMu_HistDrawer()
