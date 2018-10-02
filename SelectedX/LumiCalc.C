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

void EE_LumiCalc ();
void MuMu_LumiCalc ();
void EMu_LumiCalc ();

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

const Double_t LumiDefault = 35867.0;


void LumiCalc ( TString whichX = "" )
{
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        cout << "\n**********       EE_HistDrawer       **********" << endl;
        EE_LumiCalc();
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        cout << "\n**********      MuMu_HistDrawer      **********" << endl;
        MuMu_LumiCalc();
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        cout << "\n**********      EMu_HistDrawer     ************" << endl;
        EMu_LumiCalc();
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ----------------------------- Electron Channel ------------------------------ ///
void EE_LumiCalc ()
{
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
//    Mgr.GetProc( _EE_SingleElectron_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=1.3107e+07, MCintegral;

//################################# INVARIANT MASS #################################################

        THStack *s_mass_fine, *s_mass;
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine[5], *h_bkg_mass[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EE_DY_Full; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 4 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }
            f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            h_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_bkg_mass[iter]->SetFillColor(iter+2);

            h_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_bkg_mass[iter]->SetLineColor(iter+2);

            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );

            if ( iter == 1 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_mass_"+Mgr.Procname[_EE_WJets], h_bkg_mass[iter] );

                h_bkg_mass_fine[iter]->SetFillColor(iter+2);
                h_bkg_mass[iter]->SetFillColor(iter+2);

                h_bkg_mass_fine[iter]->SetLineColor(iter+2);
                h_bkg_mass[iter]->SetLineColor(iter+2);

                h_bkg_mass_fine[iter]->SetDirectory(0);
                h_bkg_mass[iter]->SetDirectory(0);

                s_mass_fine->Add( h_bkg_mass_fine[iter] );
                s_mass->Add( h_bkg_mass[iter] );

            } // End of WJets

            iter++;            
        } // End of for(bkg)       

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_fine, *h_DY_mass;

        f_DY->GetObject( "h_mass_fine_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine );
        f_DY->GetObject( "h_mass_"+Mgr.Procname[_EE_DY_Full], h_DY_mass );

        h_DY_mass_fine->SetFillColor(kOrange);
        h_DY_mass->SetFillColor(kOrange);

        h_DY_mass_fine->SetLineColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);

        h_DY_mass_fine->SetDirectory(0);
        h_DY_mass->SetDirectory(0);

        s_mass_fine->Add( h_DY_mass_fine );
        s_mass->Add( h_DY_mass );

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine, *h_data_mass;

        Mgr.GetProc(_EE_DoubleEG_Full);
//        Mgr.GetProc(_EE_SingleElectron_Full);

        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass );

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);

        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);

// ----------------------------- Calculation ---------------------------------------- //

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events (default): " << MCintegral << "+-" << MCerror << endl;

        TH1D* h_MCstack_mass_fine = ((TH1D*)(s_mass_fine->GetStack()->Last()));
        TH1D* h_MCstack_mass = ((TH1D*)(s_mass->GetStack()->Last()));

        TH1D* h_MCstack_mass_fine_MULTIPLIED;
        TH1D* h_MCstack_mass_MULTIPLIED;

        Int_t size = 200;
        Double_t ChiSq_fine[size];
        Double_t ChiSq[size];
        Double_t factor[size];
        Double_t factor_start = 0.9;
        Double_t factor_end = 1.1;
        Double_t step = (factor_end - factor_start) / size;

        for ( Int_t i=0; i<size; i++ )
        {
            factor[i] = factor_start + ( i * step );

            h_MCstack_mass_fine_MULTIPLIED = ( (TH1D*) (h_MCstack_mass_fine->Clone( "h_MCstack_mass_fine_MULTIPLIED" ) ) );
            h_MCstack_mass_MULTIPLIED = ( (TH1D*) ( h_MCstack_mass->Clone( "h_MCstack_mass_MULTIPLIED" ) ) );

            for ( Int_t k=1; k<h_MCstack_mass_fine->GetSize()-1; k++ )
            {
                Double_t cont = h_MCstack_mass_fine_MULTIPLIED->GetBinContent(k);
                h_MCstack_mass_fine_MULTIPLIED->SetBinContent( k, cont * factor[i] );
                if ( k < h_MCstack_mass_MULTIPLIED->GetSize()-1 )
                {
                    cont = h_MCstack_mass_MULTIPLIED->GetBinContent(k);
                    h_MCstack_mass_MULTIPLIED->SetBinContent( k, cont * factor[i] );
                }
            }

            ChiSq[i] = 0;
            ChiSq_fine[i] = 0;
            for (Int_t j=1; j<h_MCstack_mass_fine_MULTIPLIED->GetSize()-1; j++)
            {
                if ( h_data_mass_fine->GetBinContent(j) != 0 )
                {
                    ChiSq_fine[i] += (h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j)) *
                                     (h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j)) / h_data_mass_fine->GetBinContent(j);
                }
                else if ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) != 0 )
                {
                    ChiSq_fine[i] += (h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j)) *
                                     (h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j)) /
                                      h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j);
                }

                if ( j < h_MCstack_mass_MULTIPLIED->GetSize()-1 )
                {
                    if ( h_data_mass_fine->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += (h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j)) *
                                         (h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j)) / h_data_mass->GetBinContent(j);
                    }
                    else if ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += (h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j)) *
                                         (h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j)) / h_MCstack_mass_MULTIPLIED->GetBinContent(j);
                    }

                }
            }
        }

        TGraph* ChiSq_graph_fine = new TGraph(size, factor, ChiSq_fine);
        TGraph* ChiSq_graph = new TGraph(size, factor, ChiSq);

        TCanvas* c_fine = new TCanvas("ChiSq_fineBin", "ChiSq_fineBin", 1000, 1000);
        ChiSq_graph_fine->Draw();
        c_fine->Update();

        TCanvas* c = new TCanvas("ChiSq", "ChiSq", 1000, 1000);
        c->cd();
        ChiSq_graph->Draw();
        c->Update();

} // End of EE_LumiCalc()


/// -------------------------------- Muon Channel ------------------------------------ ///
void MuMu_LumiCalc()
{
    LocalFileMgr Mgr;

    Mgr.GetProc( _MuMu_DY_Full );
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY = new TFile( name_DY, "READ" );
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.GetProc( _MuMu_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.GetProc( _MuMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=2.25081e+07, MCintegral;

//################################# INVARIANT MASS #################################################

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

} // End of MuMu_HistDrawer()

/// -------------------------------- EMu events ------------------------------------ ///
void EMu_LumiCalc ()
{
    LocalFileMgr Mgr;

    Mgr.GetProc( _EMu_Bkg_Full );
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile( name_bkg, "READ" );
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.GetProc( _EMu_SingleMuon_Full );
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile( name_data, "READ" );
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    Double_t dataerror, MCerror, dataintegral=348650, MCintegral;

//################################# INVARIANT MASS #################################################

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

} // End of EMu_HistDrawer()
