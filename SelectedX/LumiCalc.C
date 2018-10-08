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
        Double_t factor_start = 0.85;
        Double_t factor_end = 1.05;
        Double_t step = (factor_end - factor_start) / size;
        Double_t smallest_chisq_fine = 1e6;
        Double_t smallest_factor_fine = -1;
        Double_t smallest_chisq = 1e6;
        Double_t smallest_factor = -1;

        myProgressBar_t bar(size);

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
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_data_mass_fine->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }
                else if ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) != 0 )
                {
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }

                if ( j < h_MCstack_mass_MULTIPLIED->GetSize()-1 )
                {
                    if ( h_data_mass_fine->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_data_mass->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }
                    else if ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }

                }
            }
            if ( ChiSq_fine[i] < smallest_chisq_fine )
            {
                smallest_chisq_fine = ChiSq_fine[i];
                smallest_factor_fine = factor[i];
            }
            if ( ChiSq[i] < smallest_chisq )
            {
                smallest_chisq = ChiSq[i];
                smallest_factor = factor[i];
            }
            bar.Draw(i);
        }

        cout << "The luminosity value that gives the smallest chi^2 for fine-binned ee mass histogram is " << LumiDefault * smallest_factor_fine;
        cout << "   (factor: " << smallest_factor_fine << "; chi^2 value: " << smallest_chisq_fine << ")" << endl;
        cout << "The luminosity value that gives the smallest chi^2 for ee mass histogram is " << LumiDefault * smallest_factor;
        cout << "   (factor: " << smallest_factor << "; chi^2 value: " << smallest_chisq << ")" << endl;

// ----------------------------- Drawing ------------------------------------ //

        TCanvas* c_fine = new TCanvas( "ChiSq_fineBin", "ChiSq_fineBin", 1000, 1000 );
        TGraph* ChiSq_graph_fine = new TGraph( size, factor, ChiSq_fine );
        ChiSq_graph_fine->SetTitle( "Fine-binned ee invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph_fine->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph_fine->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph_fine->SetLineWidth(2);
        ChiSq_graph_fine->GetYaxis()->SetRangeUser(0, 26);

        TLine* l_vertical_fine = new TLine( smallest_factor_fine, 0, smallest_factor_fine, 26 );
        l_vertical_fine->SetLineColor(kRed);
        l_vertical_fine->SetLineWidth(2);

        ChiSq_graph_fine->Draw();
        l_vertical_fine->Draw("same");
        c_fine->SetGridx();
        c_fine->SetTickx();
        c_fine->SetGridy();
        c_fine->SetTicky();
        c_fine->Update();

        TCanvas* c = new TCanvas( "ChiSq", "ChiSq", 1000, 1000 );
        TGraph* ChiSq_graph = new TGraph(size, factor, ChiSq);
        ChiSq_graph->SetTitle( "ee invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph->SetLineWidth(2);
        ChiSq_graph->GetYaxis()->SetRangeUser(0, 5700);

        TLine* l_vertical = new TLine( smallest_factor, 0, smallest_factor, 5700 );
        l_vertical->SetLineColor(kRed);
        l_vertical->SetLineWidth(2);

        c->cd();
        ChiSq_graph->Draw();
        l_vertical->Draw("same");
        c->SetGridx();
        c->SetTickx();
        c->SetGridy();
        c->SetTicky();
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

        THStack *s_mass_fine, *s_mass;
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine[5], *h_bkg_mass[5];
        Int_t iter = 0;

        for ( SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _MuMu_DY_Full; pr=SelProc_t((int)(pr-1)) )
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
                f_bkg->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_mass_"+Mgr.Procname[_MuMu_WJets], h_bkg_mass[iter] );

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

        f_DY->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine );
        f_DY->GetObject( "h_mass_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass );

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

        f_data->GetObject( "h_mass_fine_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine );
        f_data->GetObject( "h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass );

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);

        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);

// --------------------------------- Calculation -------------------------------------- //

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

        TH1D* h_MCstack_mass_fine = ((TH1D*)(s_mass_fine->GetStack()->Last()));
        TH1D* h_MCstack_mass = ((TH1D*)(s_mass->GetStack()->Last()));

        TH1D* h_MCstack_mass_fine_MULTIPLIED;
        TH1D* h_MCstack_mass_MULTIPLIED;

        Int_t size = 200;
        Double_t ChiSq_fine[size];
        Double_t ChiSq[size];
        Double_t factor[size];
        Double_t factor_start = 0.8;
        Double_t factor_end = 1.0;
        Double_t step = (factor_end - factor_start) / size;
        Double_t smallest_chisq_fine = 1e6;
        Double_t smallest_factor_fine = -1;
        Double_t smallest_chisq = 1e6;
        Double_t smallest_factor = -1;

        myProgressBar_t bar(size);

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
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_data_mass_fine->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }
                else if ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) != 0 )
                {
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }

                if ( j < h_MCstack_mass_MULTIPLIED->GetSize()-1 )
                {
                    if ( h_data_mass_fine->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_data_mass->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }
                    else if ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }

                }
            }
            if ( ChiSq_fine[i] < smallest_chisq_fine )
            {
                smallest_chisq_fine = ChiSq_fine[i];
                smallest_factor_fine = factor[i];
            }
            if ( ChiSq[i] < smallest_chisq )
            {
                smallest_chisq = ChiSq[i];
                smallest_factor = factor[i];
            }
            bar.Draw(i);
        }

        cout << "The luminosity value that gives the smallest chi^2 for fine-binned MuMu mass histogram is " << LumiDefault * smallest_factor_fine;
        cout << "   (factor: " << smallest_factor_fine << "; chi^2 value: " << smallest_chisq_fine << ")" << endl;
        cout << "The luminosity value that gives the smallest chi^2 for MuMu mass histogram is " << LumiDefault * smallest_factor;
        cout << "   (factor: " << smallest_factor << "; chi^2 value: " << smallest_chisq << ")" << endl;

// ----------------------------- Drawing ------------------------------------ //

        TCanvas* c_fine = new TCanvas( "ChiSq_fineBin", "ChiSq_fineBin", 1000, 1000 );
        TGraph* ChiSq_graph_fine = new TGraph( size, factor, ChiSq_fine );
        ChiSq_graph_fine->SetTitle( "Fine-binned #mu#mu invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph_fine->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph_fine->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph_fine->SetLineWidth(2);
        ChiSq_graph_fine->GetYaxis()->SetRangeUser(0, 50);

        TLine* l_vertical_fine = new TLine( smallest_factor_fine, 0, smallest_factor_fine, 50 );
        l_vertical_fine->SetLineColor(kRed);
        l_vertical_fine->SetLineWidth(2);

        ChiSq_graph_fine->Draw();
        l_vertical_fine->Draw("same");
        c_fine->SetGridx();
        c_fine->SetTickx();
        c_fine->SetGridy();
        c_fine->SetTicky();
        c_fine->Update();

        TCanvas* c = new TCanvas( "ChiSq", "ChiSq", 1000, 1000 );
        TGraph* ChiSq_graph = new TGraph(size, factor, ChiSq);
        ChiSq_graph->SetTitle( "#mu#mu invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph->SetLineWidth(2);
        ChiSq_graph->GetYaxis()->SetRangeUser(0, 9000);

        TLine* l_vertical = new TLine( smallest_factor, 0, smallest_factor, 9000 );
        l_vertical->SetLineColor(kRed);
        l_vertical->SetLineWidth(2);

        c->cd();
        ChiSq_graph->Draw();
        l_vertical->Draw("same");
        c->SetGridx();
        c->SetTickx();
        c->SetGridy();
        c->SetTicky();
        c->Update();

} // End of MuMu_LumiCalc()

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

        THStack *s_mass_fine, *s_mass, *s_SS_mass_fine, *s_SS_mass;
        s_mass_fine = new THStack("s_mass_fine", "");
        s_mass = new THStack("s_mass", "");
        s_SS_mass_fine = new THStack("s_SS_mass_fine", "");
        s_SS_mass = new THStack("s_SS_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_fine[5], *h_bkg_mass[4];
        Int_t iter = 0;

        for ( SelProc_t pr = _EMu_VVnST; pr > _EndOf_EMu_Data_Normal; pr=SelProc_t((int)(pr-1)) )
        {
            if ( iter > 3 )
            {
                cout << "Error: Iteration exceeds histogram limits!!!" << endl;
                break;
            }

            f_bkg->GetObject( "h_emu_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter] );
            f_bkg->GetObject( "h_emu_mass_"+Mgr.Procname[pr], h_bkg_mass[iter] );

            h_bkg_mass_fine[iter]->SetFillColor(iter+2);
            h_bkg_mass[iter]->SetFillColor(iter+2);

            h_bkg_mass_fine[iter]->SetLineColor(iter+2);
            h_bkg_mass[iter]->SetLineColor(iter+2);

            h_bkg_mass_fine[iter]->SetDirectory(0);
            h_bkg_mass[iter]->SetDirectory(0);

            s_mass_fine->Add( h_bkg_mass_fine[iter] );
            s_mass->Add( h_bkg_mass[iter] );

            if ( iter == 0 )
            {   // ---- WJETS ---- //
                iter++;
                f_bkg->GetObject( "h_emu_mass_fine_"+Mgr.Procname[_EMu_WJets], h_bkg_mass_fine[iter] );
                f_bkg->GetObject( "h_emu_mass_"+Mgr.Procname[_EMu_WJets], h_bkg_mass[iter] );

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

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_fine, *h_data_mass;

        f_data->GetObject( "h_emu_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine );
        f_data->GetObject( "h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass );

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerStyle(kFullDotLarge);

        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass->SetMarkerColor(kBlack);

        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass->SetLineColor(kBlack);

        h_data_mass_fine->SetDirectory(0);
        h_data_mass->SetDirectory(0);

// --------------------------------- Calculation -------------------------------------- //

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);

        std::cout << "data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;

        TH1D* h_MCstack_mass_fine = ((TH1D*)(s_mass_fine->GetStack()->Last()));
        TH1D* h_MCstack_mass = ((TH1D*)(s_mass->GetStack()->Last()));

        TH1D* h_MCstack_mass_fine_MULTIPLIED;
        TH1D* h_MCstack_mass_MULTIPLIED;

        Int_t size = 200;
        Double_t ChiSq_fine[size];
        Double_t ChiSq[size];
        Double_t factor[size];
        Double_t factor_start = 0.8;
        Double_t factor_end = 1.0;
        Double_t step = (factor_end - factor_start) / size;
        Double_t smallest_chisq_fine = 1e6;
        Double_t smallest_factor_fine = -1;
        Double_t smallest_chisq = 1e6;
        Double_t smallest_factor = -1;

        myProgressBar_t bar(size);

        for ( Int_t i=0; i<size; i++ )
        {
            factor[i] = factor_start + ( i * step );

            h_MCstack_mass_fine_MULTIPLIED = ( (TH1D*) (h_MCstack_mass_fine->Clone( "h_MCstack_mass_fine_MULTIPLIED" ) ) );
            h_MCstack_mass_MULTIPLIED = ( (TH1D*) ( h_MCstack_mass->Clone( "h_MCstack_mass_MULTIPLIED" ) ) );

            for ( Int_t k=1; k<h_MCstack_mass_fine->GetSize()-1; k++ ) // Multiplying each MC bin by a certain factor (changing luminosity)
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
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_data_mass_fine->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }
                else if ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) != 0 )
                {
                    ChiSq_fine[i] += ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) *
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) - h_data_mass_fine->GetBinContent(j) ) /
                                     ( h_MCstack_mass_fine_MULTIPLIED->GetBinContent(j) * (h_data_mass_fine->GetSize() - 2) );
                }

                if ( j < h_MCstack_mass_MULTIPLIED->GetSize()-1 )
                {
                    if ( h_data_mass_fine->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_data_mass->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }
                    else if ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) != 0 )
                    {
                        ChiSq[i] += ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) *
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) - h_data_mass->GetBinContent(j) ) /
                                    ( h_MCstack_mass_MULTIPLIED->GetBinContent(j) * (h_data_mass->GetSize() - 2) );
                    }

                }
            }
            if ( ChiSq_fine[i] < smallest_chisq_fine )
            {
                smallest_chisq_fine = ChiSq_fine[i];
                smallest_factor_fine = factor[i];
            }
            if ( ChiSq[i] < smallest_chisq )
            {
                smallest_chisq = ChiSq[i];
                smallest_factor = factor[i];
            }
            bar.Draw(i);
        }

        cout << "The luminosity value that gives the smallest chi^2 for fine-binned EMu mass histogram is " << LumiDefault * smallest_factor_fine;
        cout << "   (factor: " << smallest_factor_fine << "; chi^2 value: " << smallest_chisq_fine << ")" << endl;
        cout << "The luminosity value that gives the smallest chi^2 for EMu mass histogram is " << LumiDefault * smallest_factor;
        cout << "   (factor: " << smallest_factor << "; chi^2 value: " << smallest_chisq << ")" << endl;

// ----------------------------- Drawing ------------------------------------ //

        TCanvas* c_fine = new TCanvas( "ChiSq_fineBin", "ChiSq_fineBin", 1000, 1000 );
        TGraph* ChiSq_graph_fine = new TGraph( size, factor, ChiSq_fine );
        ChiSq_graph_fine->SetTitle( "Fine-binned e#mu invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph_fine->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph_fine->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph_fine->SetLineWidth(2);
        ChiSq_graph_fine->GetYaxis()->SetRangeUser(0, 2);

        TLine* l_vertical_fine = new TLine( smallest_factor_fine, 0, smallest_factor_fine, 2 );
        l_vertical_fine->SetLineColor(kRed);
        l_vertical_fine->SetLineWidth(2);

        ChiSq_graph_fine->Draw();
        l_vertical_fine->Draw("same");
        c_fine->SetGridx();
        c_fine->SetTickx();
        c_fine->SetGridy();
        c_fine->SetTicky();
        c_fine->Update();

        TCanvas* c = new TCanvas( "ChiSq", "ChiSq", 1000, 1000 );
        TGraph* ChiSq_graph = new TGraph(size, factor, ChiSq);
        ChiSq_graph->SetTitle( "e#mu invariant mass histogram's #chi^{2} dependancy on luminosity change" );
        ChiSq_graph->GetXaxis()->SetTitle( "Luminosity factor" );
        ChiSq_graph->GetYaxis()->SetTitle( "#frac{#chi^{2}}{ndof}" );
        ChiSq_graph->SetLineWidth(2);
        ChiSq_graph->GetYaxis()->SetRangeUser(0, 200);

        TLine* l_vertical = new TLine( smallest_factor, 0, smallest_factor, 200 );
        l_vertical->SetLineColor(kRed);
        l_vertical->SetLineWidth(2);

        c->cd();
        ChiSq_graph->Draw();
        l_vertical->Draw("same");
        c->SetGridx();
        c->SetTickx();
        c->SetGridy();
        c->SetTicky();
        c->Update();

} // End of EMu_HistDrawer()
