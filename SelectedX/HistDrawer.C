#include <TROOT.h>
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
#include <TVectorT.h>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

void EE_HistDrawer (TString whichGraphs, TString type);
void MuMu_HistDrawer (TString whichGraphs, TString type);
void EMu_HistDrawer (TString whichGraphs, TString type);
void Est_HistDrawer (Int_t FR_systErr);
void TEST_HistDrawer (TString whichGraphs , TString type);
void Test_RocCorr ();
Double_t CompChiSquared (TH1D *h_data, THStack *s_MC);
Double_t CompAvgDataMCDifference (TH1D *h_data, THStack *s_MC);
void removeNegativeBins(TH1D *h);

// -- Drell-Yan mass bins -- //
const Int_t binnum = 43;
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

// -- Mass bins for background estimation -- //
const Int_t binnum2 = 86;
const Double_t massbins2[87] = {15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64,
                                66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5,
                                115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200,
                                210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830,
                                915, 1000, 1250, 1500, 2250, 3000};


void HistDrawer (TString WhichX = "", Int_t FR_systErr = 1, TString WhichGraphs = "ALL", TString type = "")
{
    TString whichX = WhichX;
    whichX.ToUpper();
    TString whichGraphs = WhichGraphs;
    whichGraphs.ToUpper();
    Int_t Xselected = 0;
    if (whichX.Contains("EE"))
    {
        Xselected++;
        cout << "\n*******      EE_HistDrawer (" << whichGraphs << " " << type << ")      *******" << endl;
        EE_HistDrawer(whichGraphs, type);
    }
    if (whichX.Contains("MUMU"))
    {
        Xselected++;
        cout << "\n*****  MuMu_HistDrawer (" << whichGraphs << " " << type << ")  *****" << endl;
        MuMu_HistDrawer(whichGraphs, type);
    }
    if (whichX.Contains("EMU"))
    {
        Xselected++;
        cout << "\n*****   EMu_HistDrawer (" << whichGraphs << " " << type << ")  *****" << endl;
        EMu_HistDrawer(whichGraphs, type);
    }
    if (whichX.Contains("EST") && !whichX.Contains("TEST"))
    {
        Xselected++;
        cout << "\n*****   Est_HistDrawer (" << FR_systErr << ")  *****" << endl;
        Est_HistDrawer(FR_systErr);
    }
    if (whichX.Contains("TEST"))
    {
        Xselected++;
        cout << "\n*****   Test_RocCorr ()  *****" << endl;
        Test_RocCorr();
    }
    if (whichX.Contains("DENSITY"))
    {
        Xselected++;
        cout << "\n*****   TEST_HistDrawer ()  *****" << endl;
        TEST_HistDrawer(whichGraphs,type);
    }
    if (Xselected == 0) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ############################################################################# ///
/// ----------------------------- Electron Channel ------------------------------ ///
/// ############################################################################# ///
void EE_HistDrawer (TString whichGraphs, TString type)
{
    if (!whichGraphs.Length())
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc(_EE_DY_Full);
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root";
    TFile* f_DY = new TFile(name_DY, "READ");
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_EE_Bkg_Full);
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root";
    TFile* f_bkg = new TFile(name_bkg, "READ");
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_EE_DoubleEG_Full);
//    Mgr.SetProc(_EE_SingleElectron_Full);
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root";
    TFile* f_data = new TFile(name_data, "READ");
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root" << " opened successfully" << endl;

//################################# INVARIANT MASS #################################################

    if(whichGraphs=="ALL" || whichGraphs=="INVMASS")
    {
        count_drawn++;

        THStack *s_mass_before_PUCorr = new THStack("s_mass_before_PUCorr", "");
        THStack *s_mass_before_EffCorr = new THStack("s_mass_before_EffCorr", "");
        THStack *s_mass_before_PVzCorr = new THStack("s_mass_before_PVzCorr", "");
        THStack *s_mass_before_L1Corr = new THStack("s_mass_before_L1Corr", "");
        THStack *s_mass_before_TopPtCorr = new THStack("s_mass_before_TopPtCorr", "");
        THStack *s_mass_fine = new THStack("s_mass_fine", "");
        THStack *s_mass = new THStack("s_mass", "");
        THStack *s_mass2 = new THStack("s_mass2", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_before_PUCorr[9], *h_bkg_mass_before_EffCorr[9], *h_bkg_mass_before_PVzCorr[9], *h_bkg_mass_before_L1Corr[9],
             *h_bkg_mass_before_TopPtCorr[9], *h_bkg_mass_fine[9], *h_bkg_mass[9], *h_bkg_mass2[9];
        Int_t iter = 0;

        for (SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _EE_QCDEMEnriched_Full)
            {
                iter++;
                continue;
            }

            f_bkg->GetObject("h_mass_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PUCorr[iter]);
            f_bkg->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_before_EffCorr[iter]);
            f_bkg->GetObject("h_mass_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PVzCorr[iter]);
            f_bkg->GetObject("h_mass_before_L1Corr_"+Mgr.Procname[pr], h_bkg_mass_before_L1Corr[iter]);
            f_bkg->GetObject("h_mass_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_mass_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter]);
            f_bkg->GetObject("h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter]);
            f_bkg->GetObject("h_mass2_"+Mgr.Procname[pr], h_bkg_mass2[iter]);
            removeNegativeBins(h_bkg_mass_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_mass_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_mass_fine[iter]);
            removeNegativeBins(h_bkg_mass[iter]);
            removeNegativeBins(h_bkg_mass2[iter]);

            Color_t color = kBlack;
            if (pr == _EE_QCDEMEnriched_Full) color = kRed + 3;
            if (pr == _EE_WJets_Full) color = kRed - 2;
            if (pr == _EE_WW) color = kMagenta - 5;
            if (pr == _EE_WZ) color = kMagenta - 2;
            if (pr == _EE_ZZ) color = kMagenta - 6;
            if (pr == _EE_tbarW) color = kGreen - 2;
            if (pr == _EE_tW) color = kGreen + 2;
            if (pr == _EE_ttbar_Full) color = kCyan + 2;
            if (pr == _EE_DYTauTau_Full) color = kOrange - 5;

            h_bkg_mass_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetDirectory(0);
            s_mass_before_PUCorr->Add(h_bkg_mass_before_PUCorr[iter]);

            h_bkg_mass_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
            s_mass_before_EffCorr->Add(h_bkg_mass_before_EffCorr[iter]);

            h_bkg_mass_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetDirectory(0);
            s_mass_before_PVzCorr->Add(h_bkg_mass_before_PVzCorr[iter]);

            h_bkg_mass_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetDirectory(0);
            s_mass_before_L1Corr->Add(h_bkg_mass_before_L1Corr[iter]);

            h_bkg_mass_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetDirectory(0);
            s_mass_before_TopPtCorr->Add(h_bkg_mass_before_TopPtCorr[iter]);

            h_bkg_mass_fine[iter]->SetFillColor(color);
            h_bkg_mass_fine[iter]->SetLineColor(color);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            s_mass_fine->Add(h_bkg_mass_fine[iter]);

            h_bkg_mass[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetDirectory(0);            
            s_mass->Add(h_bkg_mass[iter]);

            h_bkg_mass2[iter]->SetFillColor(color);
            h_bkg_mass2[iter]->SetLineColor(color);
            h_bkg_mass2[iter]->SetDirectory(0);
            s_mass2->Add(h_bkg_mass2[iter]);

            iter++;

            if (pr == _EE_WJets_Full) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if (pr == _EE_tW) // next -- ttbar
                pr = _EE_VVnST;
            if (pr == _EE_DYTauTau_Full) // last
                break;

        } // End of for(bkg)       

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_before_PUCorr, *h_DY_mass_before_EffCorr, *h_DY_mass_before_PVzCorr, *h_DY_mass_before_L1Corr,
             *h_DY_mass_before_TopPtCorr, *h_DY_mass_fine, *h_DY_mass, *h_DY_mass2;

        f_DY->GetObject("h_mass_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_before_PUCorr);
        f_DY->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_before_EffCorr);
        f_DY->GetObject("h_mass_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_before_PVzCorr);
        f_DY->GetObject("h_mass_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_before_L1Corr);
        f_DY->GetObject("h_mass_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_before_TopPtCorr);
        f_DY->GetObject("h_mass_fine_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_fine);
        f_DY->GetObject("h_mass_"+Mgr.Procname[_EE_DY_Full], h_DY_mass);
        f_DY->GetObject("h_mass2_"+Mgr.Procname[_EE_DY_Full], h_DY_mass2);
        removeNegativeBins(h_DY_mass_before_PUCorr);
        removeNegativeBins(h_DY_mass_before_EffCorr);
        removeNegativeBins(h_DY_mass_before_PVzCorr);
        removeNegativeBins(h_DY_mass_before_L1Corr);
        removeNegativeBins(h_DY_mass_before_TopPtCorr);
        removeNegativeBins(h_DY_mass_fine);
        removeNegativeBins(h_DY_mass);
        removeNegativeBins(h_DY_mass2);

        h_DY_mass_before_PUCorr->SetFillColor(kOrange);
        h_DY_mass_before_PUCorr->SetLineColor(kOrange);
        h_DY_mass_before_PUCorr->SetDirectory(0);
        s_mass_before_PUCorr->Add(h_DY_mass_before_PUCorr);

        h_DY_mass_before_EffCorr->SetFillColor(kOrange);
        h_DY_mass_before_EffCorr->SetLineColor(kOrange);
        h_DY_mass_before_EffCorr->SetDirectory(0);
        s_mass_before_EffCorr->Add(h_DY_mass_before_EffCorr);

        h_DY_mass_before_PVzCorr->SetFillColor(kOrange);
        h_DY_mass_before_PVzCorr->SetLineColor(kOrange);
        h_DY_mass_before_PVzCorr->SetDirectory(0);
        s_mass_before_PVzCorr->Add(h_DY_mass_before_PVzCorr);

        h_DY_mass_before_L1Corr->SetFillColor(kOrange);
        h_DY_mass_before_L1Corr->SetLineColor(kOrange);
        h_DY_mass_before_L1Corr->SetDirectory(0);
        s_mass_before_L1Corr->Add(h_DY_mass_before_L1Corr);

        h_DY_mass_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_mass_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_mass_before_TopPtCorr->SetDirectory(0);
        s_mass_before_TopPtCorr->Add(h_DY_mass_before_TopPtCorr);

        h_DY_mass_fine->SetFillColor(kOrange);
        h_DY_mass_fine->SetLineColor(kOrange);
        h_DY_mass_fine->SetDirectory(0);
        s_mass_fine->Add(h_DY_mass_fine);

        h_DY_mass->SetFillColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);
        h_DY_mass->SetDirectory(0);
        s_mass->Add(h_DY_mass);

        h_DY_mass2->SetFillColor(kOrange);
        h_DY_mass2->SetLineColor(kOrange);
        h_DY_mass2->SetDirectory(0);
        s_mass2->Add(h_DY_mass2);

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_before_PUCorr, *h_data_mass_before_EffCorr, *h_data_mass_before_PVzCorr, *h_data_mass_before_L1Corr,
             *h_data_mass_before_TopPtCorr, *h_data_mass_fine, *h_data_mass, *h_data_mass2;

        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject("h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_fine);
        f_data->GetObject("h_mass_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass);
        f_data->GetObject("h_mass2_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass2);

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass_fine->SetDirectory(0);

        h_data_mass->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerColor(kBlack);
        h_data_mass->SetLineColor(kBlack);
        h_data_mass->SetDirectory(0);

        h_data_mass2->SetMarkerStyle(kFullDotLarge);
        h_data_mass2->SetMarkerColor(kBlack);
        h_data_mass2->SetLineColor(kBlack);
        h_data_mass2->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_before_PUCorr = new myRatioPlot_t("RP_mass_before_PUCorr", s_mass_before_PUCorr, h_data_mass);
        myRatioPlot_t *RP_mass_before_EffCorr = new myRatioPlot_t("RP_mass_before_EffCorr", s_mass_before_EffCorr, h_data_mass);
        myRatioPlot_t *RP_mass_before_PVzCorr = new myRatioPlot_t("RP_mass_before_PVzCorr", s_mass_before_PVzCorr, h_data_mass);
        myRatioPlot_t *RP_mass_before_L1Corr = new myRatioPlot_t("RP_mass_before_L1Corr", s_mass_before_L1Corr, h_data_mass);
        myRatioPlot_t *RP_mass_before_TopPtCorr = new myRatioPlot_t("RP_mass_before_TopPtCorr", s_mass_before_TopPtCorr, h_data_mass);
        myRatioPlot_t *RP_mass_fine = new myRatioPlot_t("RP_mass_fine", s_mass_fine, h_data_mass_fine);
        myRatioPlot_t *RP_mass = new myRatioPlot_t("RP_mass", s_mass, h_data_mass);
        myRatioPlot_t *RP_mass2 = new myRatioPlot_t("RP_mass2", s_mass2, h_data_mass2);

        RP_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] before PU correction", 15, 3000);
//        RP_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] prie#check{s} visas pataisas", 15, 3000, "Eksp./MC");
        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] before Efficiency SF", 15, 3000);
//        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] prie#check{s} efektyvumo pataisas", 15, 3000, "Eksp./MC");
        RP_mass_before_PVzCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] before PVz correction", 15, 3000);
        RP_mass_before_L1Corr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] before L1 prefiring correction", 15, 3000);
        RP_mass_before_TopPtCorr->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] before Top quark p_{T} reweighting", 15, 3000);
        RP_mass_fine->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 60, 120);
        RP_mass->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
//        RP_mass->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}] pritaikius pataisas", 15, 3000, "Eksp./MC");
        RP_mass2->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);

        // Legend
        TLegend *legend = new TLegend(0.8, 0.45, 0.95, 0.95);

        legend->AddEntry(h_data_mass, "Data", "lp");
//        legend->AddEntry(h_data_mass, "Matavimas", "lp");
        legend->AddEntry(h_DY_mass, "DY#rightarrow #font[12]{ee}", "f");
        legend->AddEntry(h_bkg_mass[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_mass[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_mass[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_mass[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_mass[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_mass[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_mass[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_mass[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_mass[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_mass_before_PUCorr->ImportLegend(legend);
        RP_mass_before_EffCorr->ImportLegend(legend);
        RP_mass_before_PVzCorr->ImportLegend(legend);
        RP_mass_before_L1Corr->ImportLegend(legend);
        RP_mass_before_TopPtCorr->ImportLegend(legend);
        RP_mass_fine->ImportLegend(legend);
        RP_mass->ImportLegend(legend);
        RP_mass2->ImportLegend(legend);

        // Drawing
        RP_mass_before_PUCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_EffCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_PVzCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_L1Corr->Draw(0.5, 1e7, 1);
        RP_mass_before_TopPtCorr->Draw(0.5, 1e7, 1);
        RP_mass_fine->Draw(0.5, 1e7, 0);
        RP_mass->Draw(0.5, 1e7, 1);
        RP_mass2->Draw(0.5, 1e7, 1);

        Double_t dataerror, MCerror, MCerror_noSF, DYerror, dataintegral=1.3107e+07, MCintegral, MCintegral_noSF, DYintegral;
        Double_t dataerrorZ, MCerrorZ, DYerrorZ, dataintegralZ=1.3107e+07, MCintegralZ, DYintegralZ;
        Double_t dataerror_noZ=0, MCerror_noZ=0, DYerror_noZ=0, dataintegral_noZ=1.3107e+07, MCintegral_noZ, DYintegral_noZ, temp_noZ;

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);
        MCintegral_noSF = ((TH1D*)(s_mass_before_EffCorr->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror_noSF);
        DYintegral = h_DY_mass->IntegralAndError(1, h_DY_mass->GetSize()-2, DYerror);

        dataintegralZ = h_data_mass->IntegralAndError(10, 22, dataerrorZ);
        MCintegralZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ);
        DYintegralZ = h_DY_mass->IntegralAndError(10, 22, DYerrorZ);

        dataintegral_noZ = h_data_mass->IntegralAndError(1, 9, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataintegral_noZ += h_data_mass->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataerror_noZ = sqrt(dataerror_noZ);

        MCintegral_noZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCintegral_noZ += ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCerror_noZ = sqrt(MCerror_noZ);

        DYintegral_noZ = h_DY_mass->IntegralAndError(1, 9, temp_noZ);
        DYerror_noZ += temp_noZ * temp_noZ;
        DYintegral_noZ += h_DY_mass->IntegralAndError(23, h_DY_mass->GetSize()-2, temp_noZ);
        DYerror_noZ += temp_noZ * temp_noZ;
        DYerror_noZ = sqrt(DYerror_noZ);

        std::cout << "Data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;
        std::cout << "MC/Obs: " << MCintegral / dataintegral << "+-" <<
                     sqrt((dataerror / dataintegral) * (dataerror / dataintegral) +
                           (MCerror / MCintegral) * (MCerror / MCintegral)) << endl;
        std::cout << "MC events before corrections: " << MCintegral_noSF << "+-" << MCerror_noSF << endl;
        std::cout << "MC DY events: " << DYintegral << "+-" << DYerror << endl;
        std::cout << "Avg. Data and MC relative difference: " << CompAvgDataMCDifference(h_data_mass, s_mass) << endl;
        std::cout << "Chi^2: " << CompChiSquared(h_data_mass, s_mass) << endl << endl;

        std::cout << "Data events around Z: " << dataintegralZ << "+-" << dataerrorZ << endl;
        std::cout << "MC events around Z: " << MCintegralZ << "+-" << MCerrorZ << endl;
        std::cout << "DY events around Z: " << DYintegralZ << "+-" << DYerrorZ << endl << endl;
        std::cout << "Data events outside Z: " << dataintegral_noZ << "+-" << dataerror_noZ << endl;
        std::cout << "MC events outside Z: " << MCintegral_noZ << "+-" << MCerror_noZ << endl;
        std::cout << "DY events outside Z: " << DYintegral_noZ << "+-" << DYerror_noZ << endl << endl;

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################

    if (whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI"))
    {
        count_drawn++;

        THStack *s_rapi_before_PUCorr = new THStack("s_rapi_before_PUCorr", "");
        THStack *s_pT_before_PUCorr = new THStack("s_pT_before_PUCorr", "");
        THStack *s_pT_lead_before_PUCorr = new THStack("s_pT_lead_before_PUCorr", "");
        THStack *s_pT_sublead_before_PUCorr = new THStack("s_pT_sublead_before_PUCorr", "");
        THStack *s_eta_lead_before_PUCorr = new THStack("s_eta_lead_before_PUCorr", "");
        THStack *s_eta_sublead_before_PUCorr = new THStack("s_eta_sublead_before_PUCorr", "");
        THStack *s_phi_lead_before_PUCorr = new THStack("s_phi_lead_before_PUCorr", "");
        THStack *s_phi_sublead_before_PUCorr = new THStack("s_phi_sublead_before_PUCorr", "");

        THStack *s_rapi_before_EffCorr = new THStack("s_rapi_before_EffCorr", "");
        THStack *s_pT_before_EffCorr = new THStack("s_pT_before_EffCorr", "");
        THStack *s_pT_lead_before_EffCorr = new THStack("s_pT_lead_before_EffCorr", "");
        THStack *s_pT_sublead_before_EffCorr = new THStack("s_pT_sublead_before_EffCorr", "");
        THStack *s_eta_lead_before_EffCorr = new THStack("s_eta_lead_before_EffCorr", "");
        THStack *s_eta_sublead_before_EffCorr = new THStack("s_eta_sublead_before_EffCorr", "");
        THStack *s_phi_lead_before_EffCorr = new THStack("s_phi_lead_before_EffCorr", "");
        THStack *s_phi_sublead_before_EffCorr = new THStack("s_phi_sublead_before_EffCorr", "");

        THStack *s_rapi_before_PVzCorr = new THStack("s_rapi_before_PVzCorr", "");
        THStack *s_pT_before_PVzCorr = new THStack("s_pT_before_PVzCorr", "");
        THStack *s_pT_lead_before_PVzCorr = new THStack("s_pT_lead_before_PVzCorr", "");
        THStack *s_pT_sublead_before_PVzCorr = new THStack("s_pT_sublead_before_PVzCorr", "");
        THStack *s_eta_lead_before_PVzCorr = new THStack("s_eta_lead_before_PVzCorr", "");
        THStack *s_eta_sublead_before_PVzCorr = new THStack("s_eta_sublead_before_PVzCorr", "");
        THStack *s_phi_lead_before_PVzCorr = new THStack("s_phi_lead_before_PVzCorr", "");
        THStack *s_phi_sublead_before_PVzCorr = new THStack("s_phi_sublead_before_PVzCorr", "");

        THStack *s_rapi_before_L1Corr = new THStack("s_rapi_before_L1Corr", "");
        THStack *s_pT_before_L1Corr = new THStack("s_pT_before_L1Corr", "");
        THStack *s_pT_lead_before_L1Corr = new THStack("s_pT_lead_before_L1Corr", "");
        THStack *s_pT_sublead_before_L1Corr = new THStack("s_pT_sublead_before_L1Corr", "");
        THStack *s_eta_lead_before_L1Corr = new THStack("s_eta_lead_before_L1Corr", "");
        THStack *s_eta_sublead_before_L1Corr = new THStack("s_eta_sublead_before_L1Corr", "");
        THStack *s_phi_lead_before_L1Corr = new THStack("s_phi_lead_before_L1Corr", "");
        THStack *s_phi_sublead_before_L1Corr = new THStack("s_phi_sublead_before_L1Corr", "");

        THStack *s_rapi_before_TopPtCorr = new THStack("s_rapi_before_TopPtCorr", "");
        THStack *s_pT_before_TopPtCorr = new THStack("s_pT_before_TopPtCorr", "");
        THStack *s_pT_lead_before_TopPtCorr = new THStack("s_pT_lead_before_TopPtCorr", "");
        THStack *s_pT_sublead_before_TopPtCorr = new THStack("s_pT_sublead_before_TopPtCorr", "");
        THStack *s_eta_lead_before_TopPtCorr = new THStack("s_eta_lead_before_TopPtCorr", "");
        THStack *s_eta_sublead_before_TopPtCorr = new THStack("s_eta_sublead_before_TopPtCorr", "");
        THStack *s_phi_lead_before_TopPtCorr = new THStack("s_phi_lead_before_TopPtCorr", "");
        THStack *s_phi_sublead_before_TopPtCorr = new THStack("s_phi_sublead_before_TopPtCorr", "");

        THStack *s_rapi = new THStack("s_rapi", "");
        THStack *s_pT = new THStack("s_pT", "");
        THStack *s_pT_lead = new THStack("s_pT_lead", "");
        THStack *s_pT_sublead = new THStack("s_pT_sublead", "");
        THStack *s_eta_lead = new THStack("s_eta_lead", "");
        THStack *s_eta_sublead = new THStack("s_eta_sublead", "");
        THStack *s_phi_lead = new THStack("s_phi_lead", "");
        THStack *s_phi_sublead = new THStack("s_phi_sublead", "");


//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_pT_before_PUCorr[9], *h_bkg_pT_before_EffCorr[9], *h_bkg_pT_before_PVzCorr[9],
             *h_bkg_pT_before_L1Corr[9], *h_bkg_pT_before_TopPtCorr[9], *h_bkg_pT[9],
             *h_bkg_rapi_before_PUCorr[9], *h_bkg_rapi_before_EffCorr[9], *h_bkg_rapi_before_PVzCorr[9],
             *h_bkg_rapi_before_L1Corr[9], *h_bkg_rapi_before_TopPtCorr[9], *h_bkg_rapi[9],
             *h_bkg_pT_lead_before_PUCorr[9], *h_bkg_pT_lead_before_EffCorr[9], *h_bkg_pT_lead_before_PVzCorr[9],
             *h_bkg_pT_lead_before_L1Corr[9], *h_bkg_pT_lead_before_TopPtCorr[9], *h_bkg_pT_lead[9],
             *h_bkg_pT_sublead_before_PUCorr[9], *h_bkg_pT_sublead_before_EffCorr[9], *h_bkg_pT_sublead_before_PVzCorr[9],
             *h_bkg_pT_sublead_before_L1Corr[9], *h_bkg_pT_sublead_before_TopPtCorr[9], *h_bkg_pT_sublead[9],
             *h_bkg_eta_lead_before_PUCorr[9], *h_bkg_eta_lead_before_EffCorr[9], *h_bkg_eta_lead_before_PVzCorr[9],
             *h_bkg_eta_lead_before_L1Corr[9], *h_bkg_eta_lead_before_TopPtCorr[9], *h_bkg_eta_lead[9],
             *h_bkg_eta_sublead_before_PUCorr[9], *h_bkg_eta_sublead_before_EffCorr[9], *h_bkg_eta_sublead_before_PVzCorr[9],
             *h_bkg_eta_sublead_before_L1Corr[9], *h_bkg_eta_sublead_before_TopPtCorr[9], *h_bkg_eta_sublead[9],
             *h_bkg_phi_lead_before_PUCorr[9], *h_bkg_phi_lead_before_EffCorr[9], *h_bkg_phi_lead_before_PVzCorr[9],
             *h_bkg_phi_lead_before_L1Corr[9], *h_bkg_phi_lead_before_TopPtCorr[9], *h_bkg_phi_lead[9],
             *h_bkg_phi_sublead_before_PUCorr[9], *h_bkg_phi_sublead_before_EffCorr[9], *h_bkg_phi_sublead_before_PVzCorr[9],
             *h_bkg_phi_sublead_before_L1Corr[9], *h_bkg_phi_sublead_before_TopPtCorr[9], *h_bkg_phi_sublead[9];
        Int_t iter = 0;

        for (SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _EE_QCDEMEnriched_Full)
            {
                iter++;
                continue;
            }
            f_bkg->GetObject("h_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_before_PUCorr[iter]);
            f_bkg->GetObject("h_rapi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_PUCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PUCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PUCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_PUCorr[iter]);

            f_bkg->GetObject("h_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_before_EffCorr[iter]);
            f_bkg->GetObject("h_rapi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_EffCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_EffCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_EffCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_EffCorr[iter]);

            f_bkg->GetObject("h_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_before_PVzCorr[iter]);
            f_bkg->GetObject("h_rapi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_PVzCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_PVzCorr[iter]);

            f_bkg->GetObject("h_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_before_L1Corr[iter]);
            f_bkg->GetObject("h_rapi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_rapi_before_L1Corr[iter]);
            f_bkg->GetObject("h_pT_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_L1Corr[iter]);
            f_bkg->GetObject("h_eta_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_L1Corr[iter]);
            f_bkg->GetObject("h_phi_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_rapi_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_L1Corr[iter]);

            f_bkg->GetObject("h_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_rapi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_TopPtCorr[iter]);

            f_bkg->GetObject("h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter]);
            f_bkg->GetObject("h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter]);
            f_bkg->GetObject("h_pT_lead_"+Mgr.Procname[pr], h_bkg_pT_lead[iter]);
            f_bkg->GetObject("h_pT_sublead_"+Mgr.Procname[pr], h_bkg_pT_sublead[iter]);
            f_bkg->GetObject("h_eta_lead_"+Mgr.Procname[pr], h_bkg_eta_lead[iter]);
            f_bkg->GetObject("h_eta_sublead_"+Mgr.Procname[pr], h_bkg_eta_sublead[iter]);
            f_bkg->GetObject("h_phi_lead_"+Mgr.Procname[pr], h_bkg_phi_lead[iter]);
            f_bkg->GetObject("h_phi_sublead_"+Mgr.Procname[pr], h_bkg_phi_sublead[iter]);
            removeNegativeBins(h_bkg_pT[iter]);
            removeNegativeBins(h_bkg_rapi[iter]);
            removeNegativeBins(h_bkg_pT_lead[iter]);
            removeNegativeBins(h_bkg_pT_sublead[iter]);
            removeNegativeBins(h_bkg_eta_lead[iter]);
            removeNegativeBins(h_bkg_eta_sublead[iter]);
            removeNegativeBins(h_bkg_phi_lead[iter]);
            removeNegativeBins(h_bkg_phi_sublead[iter]);

            Color_t color = kBlack;
            if (pr == _EE_QCDEMEnriched_Full) color = kRed + 3;
            if (pr == _EE_WJets_Full) color = kRed - 2;
            if (pr == _EE_WW) color = kMagenta - 5;
            if (pr == _EE_WZ) color = kMagenta - 2;
            if (pr == _EE_ZZ) color = kMagenta - 6;
            if (pr == _EE_tbarW) color = kGreen - 2;
            if (pr == _EE_tW) color = kGreen + 2;
            if (pr == _EE_ttbar_Full) color = kCyan + 2;
            if (pr == _EE_DYTauTau_Full) color = kOrange - 5;

            // Before PU
            h_bkg_pT_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetDirectory(0);
            //
            s_pT_before_PUCorr->Add(h_bkg_pT_before_PUCorr[iter]);
            s_rapi_before_PUCorr->Add(h_bkg_rapi_before_PUCorr[iter]);
            s_pT_lead_before_PUCorr->Add(h_bkg_pT_lead_before_PUCorr[iter]);
            s_pT_sublead_before_PUCorr->Add(h_bkg_pT_sublead_before_PUCorr[iter]);
            s_eta_lead_before_PUCorr->Add(h_bkg_eta_lead_before_PUCorr[iter]);
            s_eta_sublead_before_PUCorr->Add(h_bkg_eta_sublead_before_PUCorr[iter]);
            s_phi_lead_before_PUCorr->Add(h_bkg_phi_lead_before_PUCorr[iter]);
            s_phi_sublead_before_PUCorr->Add(h_bkg_phi_sublead_before_PUCorr[iter]);

            // Before Eff SF
            h_bkg_pT_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetDirectory(0);
            //
            s_pT_before_EffCorr->Add(h_bkg_pT_before_EffCorr[iter]);
            s_rapi_before_EffCorr->Add(h_bkg_rapi_before_EffCorr[iter]);
            s_pT_lead_before_EffCorr->Add(h_bkg_pT_lead_before_EffCorr[iter]);
            s_pT_sublead_before_EffCorr->Add(h_bkg_pT_sublead_before_EffCorr[iter]);
            s_eta_lead_before_EffCorr->Add(h_bkg_eta_lead_before_EffCorr[iter]);
            s_eta_sublead_before_EffCorr->Add(h_bkg_eta_sublead_before_EffCorr[iter]);
            s_phi_lead_before_EffCorr->Add(h_bkg_phi_lead_before_EffCorr[iter]);
            s_phi_sublead_before_EffCorr->Add(h_bkg_phi_sublead_before_EffCorr[iter]);

            // Before PVz
            h_bkg_pT_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetDirectory(0);
            //
            s_pT_before_PVzCorr->Add(h_bkg_pT_before_PVzCorr[iter]);
            s_rapi_before_PVzCorr->Add(h_bkg_rapi_before_PVzCorr[iter]);
            s_pT_lead_before_PVzCorr->Add(h_bkg_pT_lead_before_PVzCorr[iter]);
            s_pT_sublead_before_PVzCorr->Add(h_bkg_pT_sublead_before_PVzCorr[iter]);
            s_eta_lead_before_PVzCorr->Add(h_bkg_eta_lead_before_PVzCorr[iter]);
            s_eta_sublead_before_PVzCorr->Add(h_bkg_eta_sublead_before_PVzCorr[iter]);
            s_phi_lead_before_PVzCorr->Add(h_bkg_phi_lead_before_PVzCorr[iter]);
            s_phi_sublead_before_PVzCorr->Add(h_bkg_phi_sublead_before_PVzCorr[iter]);

            // Before L1 prefiring
            h_bkg_pT_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_rapi_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_rapi_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_rapi_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetDirectory(0);
            //
            s_pT_before_L1Corr->Add(h_bkg_pT_before_L1Corr[iter]);
            s_rapi_before_L1Corr->Add(h_bkg_rapi_before_L1Corr[iter]);
            s_pT_lead_before_L1Corr->Add(h_bkg_pT_lead_before_L1Corr[iter]);
            s_pT_sublead_before_L1Corr->Add(h_bkg_pT_sublead_before_L1Corr[iter]);
            s_eta_lead_before_L1Corr->Add(h_bkg_eta_lead_before_L1Corr[iter]);
            s_eta_sublead_before_L1Corr->Add(h_bkg_eta_sublead_before_L1Corr[iter]);
            s_phi_lead_before_L1Corr->Add(h_bkg_phi_lead_before_L1Corr[iter]);
            s_phi_sublead_before_L1Corr->Add(h_bkg_phi_sublead_before_L1Corr[iter]);

            // Before top pT prefiring
            h_bkg_pT_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            //
            s_pT_before_TopPtCorr->Add(h_bkg_pT_before_TopPtCorr[iter]);
            s_rapi_before_TopPtCorr->Add(h_bkg_rapi_before_TopPtCorr[iter]);
            s_pT_lead_before_TopPtCorr->Add(h_bkg_pT_lead_before_TopPtCorr[iter]);
            s_pT_sublead_before_TopPtCorr->Add(h_bkg_pT_sublead_before_TopPtCorr[iter]);
            s_eta_lead_before_TopPtCorr->Add(h_bkg_eta_lead_before_TopPtCorr[iter]);
            s_eta_sublead_before_TopPtCorr->Add(h_bkg_eta_sublead_before_TopPtCorr[iter]);
            s_phi_lead_before_TopPtCorr->Add(h_bkg_phi_lead_before_TopPtCorr[iter]);
            s_phi_sublead_before_TopPtCorr->Add(h_bkg_phi_sublead_before_TopPtCorr[iter]);

            // After all corrections
            h_bkg_pT[iter]->SetFillColor(color);
            h_bkg_rapi[iter]->SetFillColor(color);
            h_bkg_pT_lead[iter]->SetFillColor(color);
            h_bkg_pT_sublead[iter]->SetFillColor(color);
            h_bkg_eta_lead[iter]->SetFillColor(color);
            h_bkg_eta_sublead[iter]->SetFillColor(color);
            h_bkg_phi_lead[iter]->SetFillColor(color);
            h_bkg_phi_sublead[iter]->SetFillColor(color);
            //
            h_bkg_pT[iter]->SetLineColor(color);
            h_bkg_rapi[iter]->SetLineColor(color);
            h_bkg_pT_lead[iter]->SetLineColor(color);
            h_bkg_pT_sublead[iter]->SetLineColor(color);
            h_bkg_eta_lead[iter]->SetLineColor(color);
            h_bkg_eta_sublead[iter]->SetLineColor(color);
            h_bkg_phi_lead[iter]->SetLineColor(color);
            h_bkg_phi_sublead[iter]->SetLineColor(color);
            //
            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT_lead[iter]->SetDirectory(0);
            h_bkg_pT_sublead[iter]->SetDirectory(0);
            h_bkg_eta_lead[iter]->SetDirectory(0);
            h_bkg_eta_sublead[iter]->SetDirectory(0);
            h_bkg_phi_lead[iter]->SetDirectory(0);
            h_bkg_phi_sublead[iter]->SetDirectory(0);
            //
            s_pT->Add(h_bkg_pT[iter]);
            s_rapi->Add(h_bkg_rapi[iter]);
            s_pT_lead->Add(h_bkg_pT_lead[iter]);
            s_pT_sublead->Add(h_bkg_pT_sublead[iter]);
            s_eta_lead->Add(h_bkg_eta_lead[iter]);
            s_eta_sublead->Add(h_bkg_eta_sublead[iter]);
            s_phi_lead->Add(h_bkg_phi_lead[iter]);
            s_phi_sublead->Add(h_bkg_phi_sublead[iter]);

            iter++;

            if (pr == _EE_WJets_Full) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if (pr == _EE_tW) // next -- ttbar
                pr = _EE_VVnST;
            if (pr == _EE_DYTauTau_Full) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_pT_before_PUCorr, *h_DY_pT_before_EffCorr, *h_DY_pT_before_PVzCorr,
             *h_DY_pT_before_L1Corr, *h_DY_pT_before_TopPtCorr, *h_DY_pT,
             *h_DY_rapi_before_PUCorr, *h_DY_rapi_before_EffCorr, *h_DY_rapi_before_PVzCorr,
             *h_DY_rapi_before_L1Corr, *h_DY_rapi_before_TopPtCorr, *h_DY_rapi,
             *h_DY_pT_lead_before_PUCorr, *h_DY_pT_lead_before_EffCorr, *h_DY_pT_lead_before_PVzCorr,
             *h_DY_pT_lead_before_L1Corr, *h_DY_pT_lead_before_TopPtCorr, *h_DY_pT_lead,
             *h_DY_pT_sublead_before_PUCorr, *h_DY_pT_sublead_before_EffCorr, *h_DY_pT_sublead_before_PVzCorr,
             *h_DY_pT_sublead_before_L1Corr, *h_DY_pT_sublead_before_TopPtCorr, *h_DY_pT_sublead,
             *h_DY_eta_lead_before_PUCorr, *h_DY_eta_lead_before_EffCorr, *h_DY_eta_lead_before_PVzCorr,
             *h_DY_eta_lead_before_L1Corr, *h_DY_eta_lead_before_TopPtCorr, *h_DY_eta_lead,
             *h_DY_eta_sublead_before_PUCorr, *h_DY_eta_sublead_before_EffCorr, *h_DY_eta_sublead_before_PVzCorr,
             *h_DY_eta_sublead_before_L1Corr, *h_DY_eta_sublead_before_TopPtCorr, *h_DY_eta_sublead,
             *h_DY_phi_lead_before_PUCorr, *h_DY_phi_lead_before_EffCorr, *h_DY_phi_lead_before_PVzCorr,
             *h_DY_phi_lead_before_L1Corr, *h_DY_phi_lead_before_TopPtCorr, *h_DY_phi_lead,
             *h_DY_phi_sublead_before_PUCorr, *h_DY_phi_sublead_before_EffCorr, *h_DY_phi_sublead_before_PVzCorr,
             *h_DY_phi_sublead_before_L1Corr, *h_DY_phi_sublead_before_TopPtCorr, *h_DY_phi_sublead;

        f_DY->GetObject("h_pT_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_before_PUCorr);
        f_DY->GetObject("h_rapi_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi_before_PUCorr);
        f_DY->GetObject("h_pT_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_PUCorr);
        f_DY->GetObject("h_pT_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_PUCorr);
        f_DY->GetObject("h_eta_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_PUCorr);
        f_DY->GetObject("h_eta_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_PUCorr);
        f_DY->GetObject("h_phi_lead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_PUCorr);
        f_DY->GetObject("h_phi_sublead_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_PUCorr);
        removeNegativeBins(h_DY_pT_before_PUCorr);
        removeNegativeBins(h_DY_rapi_before_PUCorr);
        removeNegativeBins(h_DY_pT_lead_before_PUCorr);
        removeNegativeBins(h_DY_pT_sublead_before_PUCorr);
        removeNegativeBins(h_DY_eta_lead_before_PUCorr);
        removeNegativeBins(h_DY_eta_sublead_before_PUCorr);
        removeNegativeBins(h_DY_phi_lead_before_PUCorr);
        removeNegativeBins(h_DY_phi_sublead_before_PUCorr);
        //
        f_DY->GetObject("h_pT_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_before_EffCorr);
        f_DY->GetObject("h_rapi_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi_before_EffCorr);
        f_DY->GetObject("h_pT_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_EffCorr);
        f_DY->GetObject("h_pT_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_EffCorr);
        f_DY->GetObject("h_eta_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_EffCorr);
        f_DY->GetObject("h_eta_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_EffCorr);
        f_DY->GetObject("h_phi_lead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_EffCorr);
        f_DY->GetObject("h_phi_sublead_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_EffCorr);
        removeNegativeBins(h_DY_pT_before_EffCorr);
        removeNegativeBins(h_DY_rapi_before_EffCorr);
        removeNegativeBins(h_DY_pT_lead_before_EffCorr);
        removeNegativeBins(h_DY_pT_sublead_before_EffCorr);
        removeNegativeBins(h_DY_eta_lead_before_EffCorr);
        removeNegativeBins(h_DY_eta_sublead_before_EffCorr);
        removeNegativeBins(h_DY_phi_lead_before_EffCorr);
        removeNegativeBins(h_DY_phi_sublead_before_EffCorr);
        //
        f_DY->GetObject("h_pT_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_before_PVzCorr);
        f_DY->GetObject("h_rapi_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi_before_PVzCorr);
        f_DY->GetObject("h_pT_lead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_PVzCorr);
        f_DY->GetObject("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_PVzCorr);
        f_DY->GetObject("h_eta_lead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_PVzCorr);
        f_DY->GetObject("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_PVzCorr);
        f_DY->GetObject("h_phi_lead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_PVzCorr);
        f_DY->GetObject("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_pT_before_PVzCorr);
        removeNegativeBins(h_DY_rapi_before_PVzCorr);
        removeNegativeBins(h_DY_pT_lead_before_PVzCorr);
        removeNegativeBins(h_DY_pT_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_eta_lead_before_PVzCorr);
        removeNegativeBins(h_DY_eta_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_phi_lead_before_PVzCorr);
        removeNegativeBins(h_DY_phi_sublead_before_PVzCorr);
        //
        f_DY->GetObject("h_pT_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_before_L1Corr);
        f_DY->GetObject("h_rapi_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi_before_L1Corr);
        f_DY->GetObject("h_pT_lead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_L1Corr);
        f_DY->GetObject("h_pT_sublead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_L1Corr);
        f_DY->GetObject("h_eta_lead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_L1Corr);
        f_DY->GetObject("h_eta_sublead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_L1Corr);
        f_DY->GetObject("h_phi_lead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_L1Corr);
        f_DY->GetObject("h_phi_sublead_before_L1Corr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_L1Corr);
        removeNegativeBins(h_DY_pT_before_L1Corr);
        removeNegativeBins(h_DY_rapi_before_L1Corr);
        removeNegativeBins(h_DY_pT_lead_before_L1Corr);
        removeNegativeBins(h_DY_pT_sublead_before_L1Corr);
        removeNegativeBins(h_DY_eta_lead_before_L1Corr);
        removeNegativeBins(h_DY_eta_sublead_before_L1Corr);
        removeNegativeBins(h_DY_phi_lead_before_L1Corr);
        removeNegativeBins(h_DY_phi_sublead_before_L1Corr);
        //
        f_DY->GetObject("h_pT_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_before_TopPtCorr);
        f_DY->GetObject("h_rapi_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi_before_TopPtCorr);
        f_DY->GetObject("h_pT_lead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead_before_TopPtCorr);
        f_DY->GetObject("h_pT_sublead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead_before_TopPtCorr);
        f_DY->GetObject("h_eta_lead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead_before_TopPtCorr);
        f_DY->GetObject("h_eta_sublead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead_before_TopPtCorr);
        f_DY->GetObject("h_phi_lead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead_before_TopPtCorr);
        f_DY->GetObject("h_phi_sublead_before_TopPtCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_before_TopPtCorr);
        removeNegativeBins(h_DY_rapi_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_eta_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_eta_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_phi_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_phi_sublead_before_TopPtCorr);
        //
        f_DY->GetObject("h_pT_"+Mgr.Procname[_EE_DY_Full], h_DY_pT);
        f_DY->GetObject("h_rapi_"+Mgr.Procname[_EE_DY_Full], h_DY_rapi);
        f_DY->GetObject("h_pT_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_lead);
        f_DY->GetObject("h_pT_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_pT_sublead);
        f_DY->GetObject("h_eta_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_lead);
        f_DY->GetObject("h_eta_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_eta_sublead);
        f_DY->GetObject("h_phi_lead_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_lead);
        f_DY->GetObject("h_phi_sublead_"+Mgr.Procname[_EE_DY_Full], h_DY_phi_sublead);
        removeNegativeBins(h_DY_pT);
        removeNegativeBins(h_DY_rapi);
        removeNegativeBins(h_DY_pT_lead);
        removeNegativeBins(h_DY_pT_sublead);
        removeNegativeBins(h_DY_eta_lead);
        removeNegativeBins(h_DY_eta_sublead);
        removeNegativeBins(h_DY_phi_lead);
        removeNegativeBins(h_DY_phi_sublead);

        h_DY_pT_before_PUCorr->SetFillColor(kOrange);
        h_DY_rapi_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_PUCorr->SetLineColor(kOrange);
        h_DY_rapi_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_PUCorr->SetDirectory(0);
        h_DY_rapi_before_PUCorr->SetDirectory(0);
        h_DY_pT_lead_before_PUCorr->SetDirectory(0);
        h_DY_pT_sublead_before_PUCorr->SetDirectory(0);
        h_DY_eta_lead_before_PUCorr->SetDirectory(0);
        h_DY_eta_sublead_before_PUCorr->SetDirectory(0);
        h_DY_phi_lead_before_PUCorr->SetDirectory(0);
        h_DY_phi_sublead_before_PUCorr->SetDirectory(0);
        //
        s_pT_before_PUCorr->Add(h_DY_pT_before_PUCorr);
        s_rapi_before_PUCorr->Add(h_DY_rapi_before_PUCorr);
        s_pT_lead_before_PUCorr->Add(h_DY_pT_lead_before_PUCorr);
        s_pT_sublead_before_PUCorr->Add(h_DY_pT_sublead_before_PUCorr);
        s_eta_lead_before_PUCorr->Add(h_DY_eta_lead_before_PUCorr);
        s_eta_sublead_before_PUCorr->Add(h_DY_eta_sublead_before_PUCorr);
        s_phi_lead_before_PUCorr->Add(h_DY_phi_lead_before_PUCorr);
        s_phi_sublead_before_PUCorr->Add(h_DY_phi_sublead_before_PUCorr);

        h_DY_pT_before_EffCorr->SetFillColor(kOrange);
        h_DY_rapi_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_EffCorr->SetLineColor(kOrange);
        h_DY_rapi_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_EffCorr->SetDirectory(0);
        h_DY_rapi_before_EffCorr->SetDirectory(0);
        h_DY_pT_lead_before_EffCorr->SetDirectory(0);
        h_DY_pT_sublead_before_EffCorr->SetDirectory(0);
        h_DY_eta_lead_before_EffCorr->SetDirectory(0);
        h_DY_eta_sublead_before_EffCorr->SetDirectory(0);
        h_DY_phi_lead_before_EffCorr->SetDirectory(0);
        h_DY_phi_sublead_before_EffCorr->SetDirectory(0);
        //
        s_pT_before_EffCorr->Add(h_DY_pT_before_EffCorr);
        s_rapi_before_EffCorr->Add(h_DY_rapi_before_EffCorr);
        s_pT_lead_before_EffCorr->Add(h_DY_pT_lead_before_EffCorr);
        s_pT_sublead_before_EffCorr->Add(h_DY_pT_sublead_before_EffCorr);
        s_eta_lead_before_EffCorr->Add(h_DY_eta_lead_before_EffCorr);
        s_eta_sublead_before_EffCorr->Add(h_DY_eta_sublead_before_EffCorr);
        s_phi_sublead_before_EffCorr->Add(h_DY_phi_sublead_before_EffCorr);
        s_phi_lead_before_EffCorr->Add(h_DY_phi_lead_before_EffCorr);

        h_DY_pT_before_PVzCorr->SetFillColor(kOrange);
        h_DY_rapi_before_PVzCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PVzCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_PVzCorr->SetLineColor(kOrange);
        h_DY_rapi_before_PVzCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PVzCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_PVzCorr->SetDirectory(0);
        h_DY_rapi_before_PVzCorr->SetDirectory(0);
        h_DY_pT_lead_before_PVzCorr->SetDirectory(0);
        h_DY_pT_sublead_before_PVzCorr->SetDirectory(0);
        h_DY_eta_lead_before_PVzCorr->SetDirectory(0);
        h_DY_eta_sublead_before_PVzCorr->SetDirectory(0);
        h_DY_phi_lead_before_PVzCorr->SetDirectory(0);
        h_DY_phi_sublead_before_PVzCorr->SetDirectory(0);
        //
        s_pT_before_PVzCorr->Add(h_DY_pT_before_PVzCorr);
        s_rapi_before_PVzCorr->Add(h_DY_rapi_before_PVzCorr);
        s_pT_lead_before_PVzCorr->Add(h_DY_pT_lead_before_PVzCorr);
        s_pT_sublead_before_PVzCorr->Add(h_DY_pT_sublead_before_PVzCorr);
        s_eta_lead_before_PVzCorr->Add(h_DY_eta_lead_before_PVzCorr);
        s_eta_sublead_before_PVzCorr->Add(h_DY_eta_sublead_before_PVzCorr);
        s_phi_lead_before_PVzCorr->Add(h_DY_phi_lead_before_PVzCorr);
        s_phi_sublead_before_PVzCorr->Add(h_DY_phi_sublead_before_PVzCorr);

        h_DY_pT_before_L1Corr->SetFillColor(kOrange);
        h_DY_rapi_before_L1Corr->SetFillColor(kOrange);
        h_DY_pT_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_L1Corr->SetFillColor(kOrange);
        h_DY_eta_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_L1Corr->SetFillColor(kOrange);
        h_DY_phi_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_L1Corr->SetFillColor(kOrange);
        //
        h_DY_pT_before_L1Corr->SetLineColor(kOrange);
        h_DY_rapi_before_L1Corr->SetLineColor(kOrange);
        h_DY_pT_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_L1Corr->SetLineColor(kOrange);
        h_DY_eta_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_L1Corr->SetLineColor(kOrange);
        h_DY_phi_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_L1Corr->SetLineColor(kOrange);
        //
        h_DY_pT_before_L1Corr->SetDirectory(0);
        h_DY_rapi_before_L1Corr->SetDirectory(0);
        h_DY_pT_lead_before_L1Corr->SetDirectory(0);
        h_DY_pT_sublead_before_L1Corr->SetDirectory(0);
        h_DY_eta_lead_before_L1Corr->SetDirectory(0);
        h_DY_eta_sublead_before_L1Corr->SetDirectory(0);
        h_DY_phi_lead_before_L1Corr->SetDirectory(0);
        h_DY_phi_sublead_before_L1Corr->SetDirectory(0);
        //
        s_pT_before_L1Corr->Add(h_DY_pT_before_L1Corr);
        s_rapi_before_L1Corr->Add(h_DY_rapi_before_L1Corr);
        s_pT_lead_before_L1Corr->Add(h_DY_pT_lead_before_L1Corr);
        s_pT_sublead_before_L1Corr->Add(h_DY_pT_sublead_before_L1Corr);
        s_eta_lead_before_L1Corr->Add(h_DY_eta_lead_before_L1Corr);
        s_eta_sublead_before_L1Corr->Add(h_DY_eta_sublead_before_L1Corr);
        s_phi_lead_before_L1Corr->Add(h_DY_phi_lead_before_L1Corr);
        s_phi_sublead_before_L1Corr->Add(h_DY_phi_sublead_before_L1Corr);

        h_DY_pT_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_rapi_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_TopPtCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_rapi_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_TopPtCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_TopPtCorr->SetDirectory(0);
        h_DY_rapi_before_TopPtCorr->SetDirectory(0);
        h_DY_pT_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_pT_sublead_before_TopPtCorr->SetDirectory(0);
        h_DY_eta_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_eta_sublead_before_TopPtCorr->SetDirectory(0);
        h_DY_phi_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_phi_sublead_before_TopPtCorr->SetDirectory(0);
        //
        s_pT_before_TopPtCorr->Add(h_DY_pT_before_TopPtCorr);
        s_rapi_before_TopPtCorr->Add(h_DY_rapi_before_TopPtCorr);
        s_pT_lead_before_TopPtCorr->Add(h_DY_pT_lead_before_TopPtCorr);
        s_pT_sublead_before_TopPtCorr->Add(h_DY_pT_sublead_before_TopPtCorr);
        s_eta_lead_before_TopPtCorr->Add(h_DY_eta_lead_before_TopPtCorr);
        s_eta_sublead_before_TopPtCorr->Add(h_DY_eta_sublead_before_TopPtCorr);
        s_phi_lead_before_TopPtCorr->Add(h_DY_phi_lead_before_TopPtCorr);
        s_phi_sublead_before_TopPtCorr->Add(h_DY_phi_sublead_before_TopPtCorr);

        h_DY_pT->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT_lead->SetFillColor(kOrange);
        h_DY_pT_sublead->SetFillColor(kOrange);
        h_DY_eta_lead->SetFillColor(kOrange);
        h_DY_eta_sublead->SetFillColor(kOrange);
        h_DY_phi_lead->SetFillColor(kOrange);
        h_DY_phi_sublead->SetFillColor(kOrange);
        //
        h_DY_pT->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT_lead->SetLineColor(kOrange);
        h_DY_pT_sublead->SetLineColor(kOrange);
        h_DY_eta_lead->SetLineColor(kOrange);
        h_DY_eta_sublead->SetLineColor(kOrange);
        h_DY_phi_lead->SetLineColor(kOrange);
        h_DY_phi_sublead->SetLineColor(kOrange);
        //
        h_DY_pT->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT_lead->SetDirectory(0);
        h_DY_pT_sublead->SetDirectory(0);
        h_DY_eta_lead->SetDirectory(0);
        h_DY_eta_sublead->SetDirectory(0);
        h_DY_phi_lead->SetDirectory(0);
        h_DY_phi_sublead->SetDirectory(0);
        //
        s_pT->Add(h_DY_pT);
        s_rapi->Add(h_DY_rapi);
        s_pT_lead->Add(h_DY_pT_lead);
        s_pT_sublead->Add(h_DY_pT_sublead);
        s_eta_lead->Add(h_DY_eta_lead);
        s_eta_sublead->Add(h_DY_eta_sublead);
        s_phi_lead->Add(h_DY_phi_lead);
        s_phi_sublead->Add(h_DY_phi_sublead);


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT, *h_data_rapi, *h_data_pT_lead, *h_data_pT_sublead, *h_data_eta_lead, *h_data_eta_sublead,
             *h_data_phi_lead, *h_data_phi_sublead;

        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject("h_pT_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT);
        f_data->GetObject("h_rapi_"+Mgr.Procname[Mgr.CurrentProc], h_data_rapi);
        f_data->GetObject("h_pT_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_lead);
        f_data->GetObject("h_pT_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_pT_sublead);
        f_data->GetObject("h_eta_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_lead);
        f_data->GetObject("h_eta_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_eta_sublead);
        f_data->GetObject("h_phi_lead_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_lead);
        f_data->GetObject("h_phi_sublead_"+Mgr.Procname[Mgr.CurrentProc], h_data_phi_sublead);


        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead->SetMarkerStyle(kFullDotLarge);
        //
        h_data_pT->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);       
        h_data_pT_lead->SetMarkerColor(kBlack);       
        h_data_pT_sublead->SetMarkerColor(kBlack);       
        h_data_eta_lead->SetMarkerColor(kBlack);       
        h_data_eta_sublead->SetMarkerColor(kBlack);        
        h_data_phi_lead->SetMarkerColor(kBlack);        
        h_data_phi_sublead->SetMarkerColor(kBlack);
        //
        h_data_pT->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);       
        h_data_pT_lead->SetLineColor(kBlack);
        h_data_pT_sublead->SetLineColor(kBlack);
        h_data_eta_lead->SetLineColor(kBlack);
        h_data_eta_sublead->SetLineColor(kBlack);
        h_data_phi_lead->SetLineColor(kBlack);
        h_data_phi_sublead->SetLineColor(kBlack);


        h_data_pT->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT_lead->SetDirectory(0);
        h_data_pT_sublead->SetDirectory(0);
        h_data_eta_lead->SetDirectory(0);
        h_data_eta_sublead->SetDirectory(0);
        h_data_phi_lead->SetDirectory(0);
        h_data_phi_sublead->SetDirectory(0);


//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_pT_before_PUCorr, *RP_pT_before_EffCorr, *RP_pT_before_PVzCorr,
                      *RP_pT_before_L1Corr, *RP_pT_before_TopPtCorr, *RP_pT,
                      *RP_rapi_before_PUCorr, *RP_rapi_before_EffCorr, *RP_rapi_before_PVzCorr,
                      *RP_rapi_before_L1Corr, *RP_rapi_before_TopPtCorr, *RP_rapi,
                      *RP_pT_lead_before_PUCorr, *RP_pT_lead_before_EffCorr, *RP_pT_lead_before_PVzCorr,
                      *RP_pT_lead_before_L1Corr, *RP_pT_lead_before_TopPtCorr, *RP_pT_lead,
                      *RP_pT_sublead_before_PUCorr, *RP_pT_sublead_before_EffCorr, *RP_pT_sublead_before_PVzCorr,
                      *RP_pT_sublead_before_L1Corr, *RP_pT_sublead_before_TopPtCorr, *RP_pT_sublead,
                      *RP_eta_lead_before_PUCorr, *RP_eta_lead_before_EffCorr, *RP_eta_lead_before_PVzCorr,
                      *RP_eta_lead_before_L1Corr, *RP_eta_lead_before_TopPtCorr, *RP_eta_lead,
                      *RP_eta_sublead_before_PUCorr, *RP_eta_sublead_before_EffCorr, *RP_eta_sublead_before_PVzCorr,
                      *RP_eta_sublead_before_L1Corr, *RP_eta_sublead_before_TopPtCorr, *RP_eta_sublead,
                      *RP_phi_lead_before_PUCorr, *RP_phi_lead_before_EffCorr, *RP_phi_lead_before_PVzCorr,
                      *RP_phi_lead_before_L1Corr, *RP_phi_lead_before_TopPtCorr, *RP_phi_lead,
                      *RP_phi_sublead_before_PUCorr, *RP_phi_sublead_before_EffCorr, *RP_phi_sublead_before_PVzCorr,
                      *RP_phi_sublead_before_L1Corr, *RP_phi_sublead_before_TopPtCorr, *RP_phi_sublead;
        RP_pT_before_PUCorr = new myRatioPlot_t("RP_pT_before_PUCorr", s_pT_before_PUCorr, h_data_pT);
        RP_rapi_before_PUCorr = new myRatioPlot_t("RP_rapi_before_PUCorr", s_rapi_before_PUCorr, h_data_rapi);
        RP_pT_lead_before_PUCorr = new myRatioPlot_t("RP_pT_lead_before_PUCorr", s_pT_lead_before_PUCorr, h_data_pT_lead);
        RP_pT_sublead_before_PUCorr = new myRatioPlot_t("RP_pT_sublead_before_PUCorr", s_pT_sublead_before_PUCorr, h_data_pT_sublead);
        RP_eta_lead_before_PUCorr = new myRatioPlot_t("RP_eta_lead_before_PUCorr", s_eta_lead_before_PUCorr, h_data_eta_lead);
        RP_eta_sublead_before_PUCorr = new myRatioPlot_t("RP_eta_sublead_before_PUCorr", s_eta_sublead_before_PUCorr, h_data_eta_sublead);
        RP_phi_lead_before_PUCorr = new myRatioPlot_t("RP_phi_lead_before_PUCorr", s_phi_lead_before_PUCorr, h_data_phi_lead);
        RP_phi_sublead_before_PUCorr = new myRatioPlot_t("RP_phi_sublead_before_PUCorr", s_phi_sublead_before_PUCorr, h_data_phi_sublead);
        RP_pT_before_EffCorr = new myRatioPlot_t("RP_pT_before_EffCorr", s_pT_before_EffCorr, h_data_pT);
        RP_rapi_before_EffCorr = new myRatioPlot_t("RP_rapi_before_EffCorr", s_rapi_before_EffCorr, h_data_rapi);
        RP_pT_lead_before_EffCorr = new myRatioPlot_t("RP_pT_lead_before_EffCorr", s_pT_lead_before_EffCorr, h_data_pT_lead);
        RP_pT_sublead_before_EffCorr = new myRatioPlot_t("RP_pT_sublead_before_EffCorr", s_pT_sublead_before_EffCorr, h_data_pT_sublead);
        RP_eta_lead_before_EffCorr = new myRatioPlot_t("RP_eta_lead_before_EffCorr", s_eta_lead_before_EffCorr, h_data_eta_lead);
        RP_eta_sublead_before_EffCorr = new myRatioPlot_t("RP_eta_sublead_before_EffCorr", s_eta_sublead_before_EffCorr, h_data_eta_sublead);
        RP_phi_lead_before_EffCorr = new myRatioPlot_t("RP_phi_lead_before_EffCorr", s_phi_lead_before_EffCorr, h_data_phi_lead);
        RP_phi_sublead_before_EffCorr = new myRatioPlot_t("RP_phi_sublead_before_EffCorr", s_phi_sublead_before_EffCorr, h_data_phi_sublead);
        RP_pT_before_PVzCorr = new myRatioPlot_t("RP_pT_before_PVzCorr", s_pT_before_PVzCorr, h_data_pT);
        RP_rapi_before_PVzCorr = new myRatioPlot_t("RP_rapi_before_PVzCorr", s_rapi_before_PVzCorr, h_data_rapi);
        RP_pT_lead_before_PVzCorr = new myRatioPlot_t("RP_pT_lead_before_PVzCorr", s_pT_lead_before_PVzCorr, h_data_pT_lead);
        RP_pT_sublead_before_PVzCorr = new myRatioPlot_t("RP_pT_sublead_before_PVzCorr", s_pT_sublead_before_PVzCorr, h_data_pT_sublead);
        RP_eta_lead_before_PVzCorr = new myRatioPlot_t("RP_eta_lead_before_PVzCorr", s_eta_lead_before_PVzCorr, h_data_eta_lead);
        RP_eta_sublead_before_PVzCorr = new myRatioPlot_t("RP_eta_sublead_before_PVzCorr", s_eta_sublead_before_PVzCorr, h_data_eta_sublead);
        RP_phi_lead_before_PVzCorr = new myRatioPlot_t("RP_phi_lead_before_PVzCorr", s_phi_lead_before_PVzCorr, h_data_phi_lead);
        RP_phi_sublead_before_PVzCorr = new myRatioPlot_t("RP_phi_sublead_before_PVzCorr", s_phi_sublead_before_PVzCorr, h_data_phi_sublead);
        RP_pT_before_L1Corr = new myRatioPlot_t("RP_pT_before_L1Corr", s_pT_before_L1Corr, h_data_pT);
        RP_rapi_before_L1Corr = new myRatioPlot_t("RP_rapi_before_L1Corr", s_rapi_before_L1Corr, h_data_rapi);
        RP_pT_lead_before_L1Corr = new myRatioPlot_t("RP_pT_lead_before_L1Corr", s_pT_lead_before_L1Corr, h_data_pT_lead);
        RP_pT_sublead_before_L1Corr = new myRatioPlot_t("RP_pT_sublead_before_L1Corr", s_pT_sublead_before_L1Corr, h_data_pT_sublead);
        RP_eta_lead_before_L1Corr = new myRatioPlot_t("RP_eta_lead_before_L1Corr", s_eta_lead_before_L1Corr, h_data_eta_lead);
        RP_eta_sublead_before_L1Corr = new myRatioPlot_t("RP_eta_sublead_before_L1Corr", s_eta_sublead_before_L1Corr, h_data_eta_sublead);
        RP_phi_lead_before_L1Corr = new myRatioPlot_t("RP_phi_lead_before_L1Corr", s_phi_lead_before_L1Corr, h_data_phi_lead);
        RP_phi_sublead_before_L1Corr = new myRatioPlot_t("RP_phi_sublead_before_L1Corr", s_phi_sublead_before_L1Corr, h_data_phi_sublead);
        RP_pT_before_TopPtCorr = new myRatioPlot_t("RP_pT_before_TopPtCorr", s_pT_before_TopPtCorr, h_data_pT);
        RP_rapi_before_TopPtCorr = new myRatioPlot_t("RP_rapi_before_TopPtCorr", s_rapi_before_TopPtCorr, h_data_rapi);
        RP_pT_lead_before_TopPtCorr = new myRatioPlot_t("RP_pT_lead_before_TopPtCorr", s_pT_lead_before_TopPtCorr, h_data_pT_lead);
        RP_pT_sublead_before_TopPtCorr = new myRatioPlot_t("RP_pT_sublead_before_TopPtCorr", s_pT_sublead_before_TopPtCorr, h_data_pT_sublead);
        RP_eta_lead_before_TopPtCorr = new myRatioPlot_t("RP_eta_lead_before_TopPtCorr", s_eta_lead_before_TopPtCorr, h_data_eta_lead);
        RP_eta_sublead_before_TopPtCorr = new myRatioPlot_t("RP_eta_sublead_before_TopPtCorr", s_eta_sublead_before_TopPtCorr, h_data_eta_sublead);
        RP_phi_lead_before_TopPtCorr = new myRatioPlot_t("RP_phi_lead_before_TopPtCorr", s_phi_lead_before_TopPtCorr, h_data_phi_lead);
        RP_phi_sublead_before_TopPtCorr = new myRatioPlot_t("RP_phi_sublead_before_TopPtCorr", s_phi_sublead_before_TopPtCorr, h_data_phi_sublead);
        RP_pT = new myRatioPlot_t("RP_pT", s_pT, h_data_pT);
        RP_rapi = new myRatioPlot_t("RP_rapi", s_rapi, h_data_rapi);
        RP_pT_lead = new myRatioPlot_t("RP_pT_lead", s_pT_lead, h_data_pT_lead);
        RP_pT_sublead = new myRatioPlot_t("RP_pT_sublead", s_pT_sublead, h_data_pT_sublead);
        RP_eta_lead = new myRatioPlot_t("RP_eta_lead", s_eta_lead, h_data_eta_lead);
        RP_eta_sublead = new myRatioPlot_t("RP_eta_sublead", s_eta_sublead, h_data_eta_sublead);
        RP_phi_lead = new myRatioPlot_t("RP_phi_lead", s_phi_lead, h_data_phi_lead);
        RP_phi_sublead = new myRatioPlot_t("RP_phi_sublead", s_phi_sublead, h_data_phi_sublead);

        RP_pT_before_PUCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c] before PU correction", 0, 1000);
        RP_rapi_before_PUCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} before PU correction", -3, 3);
//        RP_rapi_before_PUCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} prie#check{s} visas pataisas", -3, 3, "Eksp./MC");
        RP_pT_lead_before_PUCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c] before PU correction", 0, 1000);
        RP_pT_sublead_before_PUCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c] before PU correction", 0, 1000);
        RP_eta_lead_before_PUCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}}) before PU correction", -3.5, 3.5);
//        RP_eta_lead_before_PUCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{1}}) prie#check{s} visas pataisas", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_PUCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}}) before PU correction", -3.5, 3.5);
//        RP_eta_sublead_before_PUCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{2}}) prie#check{s} visas pataisas", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_PUCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}}) before PU correction", -4, 4);
        RP_phi_sublead_before_PUCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}}) before PU correction", -4, 4);

        RP_pT_before_EffCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c] before Efficiency SF", 0, 1000);
        RP_rapi_before_EffCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} before Efficiency SF", -3, 3);
        RP_pT_lead_before_EffCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c] before Efficiency SF", 0, 1000);
        RP_pT_sublead_before_EffCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c] before Efficiency SF", 0, 1000);
        RP_eta_lead_before_EffCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}}) before Efficiency SF", -3.5, 3.5);
        RP_eta_sublead_before_EffCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}}) before Efficiency SF", -3.5, 3.5);
        RP_phi_lead_before_EffCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}}) before Efficiency SF", -4, 4);
        RP_phi_sublead_before_EffCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}}) before Efficiency SF", -4, 4);

        RP_pT_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c] before PVz correction", 0, 1000);
        RP_rapi_before_PVzCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} before PVz correction", -3, 3);
//        RP_rapi_before_PVzCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3, 3, "Eksp./MC");
        RP_pT_lead_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c] before PVz correction", 0, 1000);
        RP_pT_sublead_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c] before PVz correction", 0, 1000);
        RP_eta_lead_before_PVzCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}}) before PVz correction", -3.5, 3.5);
//        RP_eta_lead_before_PVzCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{1}}) prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_PVzCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}}) before PVz correction", -3.5, 3.5);
//        RP_eta_sublead_before_PVzCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{2}}) prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_PVzCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}}) before PVz correction", -4, 4);
        RP_phi_sublead_before_PVzCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}}) before PVz correction", -4, 4);

        RP_pT_before_L1Corr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c] before L1 prefiring correction", 0, 1000);
        RP_rapi_before_L1Corr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} before L1 prefiring correction", -3, 3);
//        RP_rapi_before_L1Corr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3, 3, "Eksp./MC");
        RP_pT_lead_before_L1Corr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c] before L1 prefiring correction", 0, 1000);
        RP_pT_sublead_before_L1Corr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c] before L1 prefiring correction", 0, 1000);
        RP_eta_lead_before_L1Corr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}}) before L1 prefiring correction", -3.5, 3.5);
//        RP_eta_lead_before_L1Corr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{1}}) prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_L1Corr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}}) before L1 prefiring correction", -3.5, 3.5);
//        RP_eta_sublead_before_L1Corr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{2}}) prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_L1Corr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}}) before L1 prefiring correction", -4, 4);
        RP_phi_sublead_before_L1Corr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}}) before L1 prefiring correction", -4, 4);

        RP_pT_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 1000);
        RP_rapi_before_TopPtCorr->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} before top p_{#lower[-0.2]{T}} reweighting", -3, 3);
        RP_pT_lead_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 1000);
        RP_pT_sublead_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 1000);
        RP_eta_lead_before_TopPtCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}}) before top p_{#lower[-0.2]{T}} reweighting", -3.5, 3.5);
        RP_eta_sublead_before_TopPtCorr->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}}) before top p_{#lower[-0.2]{T}} reweighting", -3.5, 3.5);
        RP_phi_lead_before_TopPtCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}}) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_phi_sublead_before_TopPtCorr->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}}) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);

        RP_pT->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#font[12]{#scale[1.3]{ee}}}} [GeV/c]", 0, 1000);
        RP_rapi->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}}", -3, 3);
//        RP_rapi->SetPlots("y_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} pritaikius pataisas", -3, 3, "Eksp./MC");
        RP_pT_lead->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{lead}}) [GeV/c]", 0, 1000);
        RP_pT_sublead->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}_{#lower[-0.2]{sublead}}) [GeV/c]", 0, 1000);
        RP_eta_lead->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{lead}})", -3.5, 3.5);
//        RP_eta_lead->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{1}}) pritaikius pataisas", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{sublead}})", -3.5, 3.5);
//        RP_eta_sublead->SetPlots("#eta_{#lower[-0.3]{SC}} (#font[12]{e}_{#lower[-0.2]{2}}) pritaikius pataisas", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{lead}})", -4, 4);
        RP_phi_sublead->SetPlots("#phi (#font[12]{e}_{#lower[-0.4]{sublead}})", -4, 4);

        TLegend *legend = new TLegend(0.5, 0.7, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_pT, "Data", "lp");
//        legend->AddEntry(h_data_pT, "Matavimas", "lp");
        legend->AddEntry(h_DY_pT, "DY#rightarrow #font[12]{ee}", "f");
        legend->AddEntry(h_bkg_pT[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_pT[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_pT[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_pT[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_pT[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_pT[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_pT[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_pT[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_pT[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_pT_before_PUCorr->ImportLegend(legend);
        RP_rapi_before_PUCorr->ImportLegend(legend);
        RP_pT_lead_before_PUCorr->ImportLegend(legend);
        RP_pT_sublead_before_PUCorr->ImportLegend(legend);
        RP_eta_lead_before_PUCorr->ImportLegend(legend);
        RP_eta_sublead_before_PUCorr->ImportLegend(legend);
        RP_phi_lead_before_PUCorr->ImportLegend(legend);
        RP_phi_sublead_before_PUCorr->ImportLegend(legend);
        RP_pT_before_EffCorr->ImportLegend(legend);
        RP_rapi_before_EffCorr->ImportLegend(legend);
        RP_pT_lead_before_EffCorr->ImportLegend(legend);
        RP_pT_sublead_before_EffCorr->ImportLegend(legend);
        RP_eta_lead_before_EffCorr->ImportLegend(legend);
        RP_eta_sublead_before_EffCorr->ImportLegend(legend);
        RP_phi_lead_before_EffCorr->ImportLegend(legend);
        RP_phi_sublead_before_EffCorr->ImportLegend(legend);
        RP_pT_before_PVzCorr->ImportLegend(legend);
        RP_rapi_before_PVzCorr->ImportLegend(legend);
        RP_pT_lead_before_PVzCorr->ImportLegend(legend);
        RP_pT_sublead_before_PVzCorr->ImportLegend(legend);
        RP_eta_lead_before_PVzCorr->ImportLegend(legend);
        RP_eta_sublead_before_PVzCorr->ImportLegend(legend);
        RP_phi_lead_before_PVzCorr->ImportLegend(legend);
        RP_phi_sublead_before_PVzCorr->ImportLegend(legend);
        RP_pT_before_L1Corr->ImportLegend(legend);
        RP_rapi_before_L1Corr->ImportLegend(legend);
        RP_pT_lead_before_L1Corr->ImportLegend(legend);
        RP_pT_sublead_before_L1Corr->ImportLegend(legend);
        RP_eta_lead_before_L1Corr->ImportLegend(legend);
        RP_eta_sublead_before_L1Corr->ImportLegend(legend);
        RP_phi_lead_before_L1Corr->ImportLegend(legend);
        RP_phi_sublead_before_L1Corr->ImportLegend(legend);
        RP_pT_before_TopPtCorr->ImportLegend(legend);
        RP_rapi_before_TopPtCorr->ImportLegend(legend);
        RP_pT_lead_before_TopPtCorr->ImportLegend(legend);
        RP_pT_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_eta_lead_before_TopPtCorr->ImportLegend(legend);
        RP_eta_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_phi_lead_before_TopPtCorr->ImportLegend(legend);
        RP_phi_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_pT->ImportLegend(legend);
        RP_rapi->ImportLegend(legend);
        RP_pT_lead->ImportLegend(legend);
        RP_pT_sublead->ImportLegend(legend);
        RP_eta_lead->ImportLegend(legend);
        RP_eta_sublead->ImportLegend(legend);
        RP_phi_lead->ImportLegend(legend);
        RP_phi_sublead->ImportLegend(legend);

        RP_pT_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_rapi_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_L1Corr->Draw(0.8, 4e7, 0);
        RP_phi_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_pT->Draw(0.8, 6e7, 0);
        RP_rapi->Draw(0.8, 6e7, 0);
        RP_pT_lead->Draw(0.8, 6e7, 0);
        RP_pT_sublead->Draw(0.8, 6e7, 0);
        RP_eta_lead->Draw(0.8, 6e7, 0);
        RP_eta_sublead->Draw(0.8, 6e7, 0);
        RP_phi_lead->Draw(0.8, 6e7, 0);
        RP_phi_sublead->Draw(0.8, 6e7, 0);

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nVTX #################################################

    if(whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP"))
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[9], *h_bkg_nVTX_before_EffCorr[9], *h_bkg_nVTX[9];
        Int_t iter = 0;

        for (SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _EE_QCDEMEnriched_Full)
            {
                iter++;
                continue;
            }
            f_bkg->GetObject("h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter]);
            f_bkg->GetObject("h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter]);
            f_bkg->GetObject("h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter]);
            removeNegativeBins(h_bkg_nVTX_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_nVTX_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_nVTX[iter]);

            Color_t color = kBlack;
            if (pr == _EE_QCDEMEnriched_Full) color = kRed + 3;
            if (pr == _EE_WJets_Full) color = kRed - 2;
            if (pr == _EE_WW) color = kMagenta - 5;
            if (pr == _EE_WZ) color = kMagenta - 2;
            if (pr == _EE_ZZ) color = kMagenta - 6;
            if (pr == _EE_tbarW) color = kGreen - 2;
            if (pr == _EE_tW) color = kGreen + 2;
            if (pr == _EE_ttbar_Full) color = kCyan + 2;
            if (pr == _EE_DYTauTau_Full) color = kOrange - 5;

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            s_nVTX_before_PUCorr->Add(h_bkg_nVTX_before_PUCorr[iter]);

            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            s_nVTX_before_EffCorr->Add(h_bkg_nVTX_before_EffCorr[iter]);

            h_bkg_nVTX[iter]->SetFillColor(color);
            h_bkg_nVTX[iter]->SetLineColor(color);
            h_bkg_nVTX[iter]->SetDirectory(0);
            s_nVTX->Add(h_bkg_nVTX[iter]);

            iter++;

            if (pr == _EE_WJets_Full) // next - WW
                pr = _EndOf_EE_VVnST_Normal;
            if (pr == _EE_tW) // next -- ttbar
                pr = _EE_VVnST;
            if (pr == _EE_DYTauTau_Full) // last
                break;
        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nVTX_before_PUCorr, *h_DY_nVTX_before_EffCorr, *h_DY_nVTX;

        f_DY->GetObject("h_nVTX_before_PUCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX_before_PUCorr);
        f_DY->GetObject("h_nVTX_before_EffCorr_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX_before_EffCorr);
        f_DY->GetObject("h_nVTX_"+Mgr.Procname[_EE_DY_Full], h_DY_nVTX);
        removeNegativeBins(h_DY_nVTX_before_PUCorr);
        removeNegativeBins(h_DY_nVTX_before_EffCorr);
        removeNegativeBins(h_DY_nVTX);

        h_DY_nVTX_before_PUCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_PUCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_PUCorr->SetDirectory(0);
        s_nVTX_before_PUCorr->Add(h_DY_nVTX_before_PUCorr);

        h_DY_nVTX_before_EffCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetDirectory(0);
        s_nVTX_before_EffCorr->Add(h_DY_nVTX_before_EffCorr);

        h_DY_nVTX->SetFillColor(kOrange);
        h_DY_nVTX->SetLineColor(kOrange);
        h_DY_nVTX->SetDirectory(0);
        s_nVTX->Add(h_DY_nVTX);

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;
        Mgr.SetProc(_EE_DoubleEG_Full);
//        Mgr.SetProc(_EE_SingleElectron_Full);

        f_data->GetObject("h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], h_data_nVTX);

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t("RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX);
        RP_nVTX_before_EffCorr = new myRatioPlot_t("RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX);
        RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_data_nVTX);

        RP_nVTX_before_PUCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.3]{#font[12]{#scale[1.2]{ee}}}} before PU correction", 0, 50);
//        RP_nVTX_before_PUCorr->SetPlots("N_{#lower[-0.2]{PV} #lower[-0.3]{#font[12]{#scale[1.2]{ee}}}} prie#check{s} pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", 0, 50, "Eksp./MC");
        RP_nVTX_before_EffCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.3]{#font[12]{#scale[1.2]{ee}}}} before Efficiency correction", 0, 50);
        RP_nVTX->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.3]{#font[12]{#scale[1.2]{ee}}}}", 0, 50);
//        RP_nVTX->SetPlots("N_{#lower[-0.2]{PV} #lower[-0.3]{#font[12]{#scale[1.2]{ee}}}}  pritaikius pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", 0, 50, "Eksp./MC");

        TLegend *legend = new TLegend(0.6, 0.75, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_nVTX, "Data", "lp");;
//        legend->AddEntry(h_data_nVTX, "Matavimas", "lp");;
        legend->AddEntry(h_DY_nVTX, "DY#rightarrow #font[12]{ee}", "f");
        legend->AddEntry(h_bkg_nVTX[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_nVTX[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_nVTX[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_nVTX[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_nVTX[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_nVTX[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_nVTX[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_nVTX[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_nVTX[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_nVTX_before_PUCorr->ImportLegend(legend);
        RP_nVTX_before_EffCorr->ImportLegend(legend);
        RP_nVTX->ImportLegend(legend);

        RP_nVTX_before_PUCorr->Draw(0.5, 3e7, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e7, 0);
        RP_nVTX->Draw(0.5, 3e7, 0);

        cout << "nVTX Chi^2 before PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX_before_PUCorr) << endl;
        cout << "nVTX Chi^2 after PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX) << endl;

        Double_t EvtPercentage = 0;
        for (Int_t i=11; i<32; i++)
        {
            EvtPercentage += h_data_nVTX->GetBinContent(i);
        }
        EvtPercentage /= h_data_nVTX->Integral(1, h_data_nVTX->GetSize()-2) * 0.01;
        cout << "There are " << EvtPercentage << "% of events in 10-30 nVTX range" << endl;

    } // End of if(nVTX)

    f_DY->Close();
    f_bkg->Close();
    f_data->Close();

} // End of EE_HistDrawer()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void MuMu_HistDrawer (TString whichGraphs , TString type)
{
    if (!whichGraphs.Length())
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc(_MuMu_DY_Full);
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY = new TFile(name_DY, "READ");
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_Bkg_Full);
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile(name_bkg, "READ");
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_SingleMuon_Full);
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile(name_data, "READ");
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

//################################# INVARIANT MASS #################################################

    if(whichGraphs=="ALL" || whichGraphs=="INVM" || whichGraphs=="INVMASS")
    {
        count_drawn++;

        THStack *s_mass_before_PUCorr = new THStack("s_mass_before_PUCorr", "");
        THStack *s_mass_before_RoccoR = new THStack("s_mass_before_RoccoR", "");
        THStack *s_mass_before_EffCorr = new THStack("s_mass_before_EffCorr", "");
        THStack *s_mass_before_PVzCorr = new THStack("s_mass_before_PVzCorr", "");
        THStack *s_mass_before_L1Corr = new THStack("s_mass_before_L1Corr", "");
        THStack *s_mass_before_TopPtCorr = new THStack("s_mass_before_TopPtCorr", "");
        THStack *s_mass_fine = new THStack("s_mass_fine", "");
        THStack *s_mass = new THStack("s_mass", "");
        THStack *s_mass2 = new THStack("s_mass2", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_before_PUCorr[9], *h_bkg_mass_before_RoccoR[9], *h_bkg_mass_before_EffCorr[9], *h_bkg_mass_before_PVzCorr[9],
             *h_bkg_mass_before_L1Corr[9], *h_bkg_mass_before_TopPtCorr[9], *h_bkg_mass_fine[9], *h_bkg_mass[9], *h_bkg_mass2[9];
        Int_t iter = 0;

        for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
//            if (pr == _MuMu_QCDMuEnriched_Full)
//            {
//                iter++;
//                continue;
//            }
            f_bkg->GetObject("h_mass_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PUCorr[iter]);
            f_bkg->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[pr], h_bkg_mass_before_RoccoR[iter]);
            f_bkg->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_before_EffCorr[iter]);
            f_bkg->GetObject("h_mass_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PVzCorr[iter]);
            f_bkg->GetObject("h_mass_before_L1Corr_"+Mgr.Procname[pr], h_bkg_mass_before_L1Corr[iter]);
            f_bkg->GetObject("h_mass_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_mass_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter]);
            f_bkg->GetObject("h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter]);
            f_bkg->GetObject("h_mass2_"+Mgr.Procname[pr], h_bkg_mass2[iter]);
            removeNegativeBins(h_bkg_mass_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_mass_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_mass_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_mass_fine[iter]);
            removeNegativeBins(h_bkg_mass[iter]);
            removeNegativeBins(h_bkg_mass2[iter]);


            Color_t color = kBlack;
            if (pr == _MuMu_QCDMuEnriched_Full) color = kRed + 3;
            if (pr == _MuMu_WJets_Full) color = kRed - 2;
            if (pr == _MuMu_WW) color = kMagenta - 5;
            if (pr == _MuMu_WZ) color = kMagenta - 2;
            if (pr == _MuMu_ZZ) color = kMagenta - 6;
            if (pr == _MuMu_tbarW) color = kGreen - 2;
            if (pr == _MuMu_tW) color = kGreen + 2;
            if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
            if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_mass_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetDirectory(0);
            s_mass_before_PUCorr->Add(h_bkg_mass_before_PUCorr[iter]);

            h_bkg_mass_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_mass_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_mass_before_RoccoR[iter]->SetDirectory(0);
            s_mass_before_RoccoR->Add(h_bkg_mass_before_RoccoR[iter]);

            h_bkg_mass_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
            s_mass_before_EffCorr->Add(h_bkg_mass_before_EffCorr[iter]);

            h_bkg_mass_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetDirectory(0);
            s_mass_before_PVzCorr->Add(h_bkg_mass_before_PVzCorr[iter]);

            h_bkg_mass_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetDirectory(0);
            s_mass_before_L1Corr->Add(h_bkg_mass_before_L1Corr[iter]);

            h_bkg_mass_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetDirectory(0);
            s_mass_before_TopPtCorr->Add(h_bkg_mass_before_TopPtCorr[iter]);

            h_bkg_mass_fine[iter]->SetFillColor(color);
            h_bkg_mass_fine[iter]->SetLineColor(color);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            s_mass_fine->Add(h_bkg_mass_fine[iter]);

            h_bkg_mass[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetDirectory(0);
            s_mass->Add(h_bkg_mass[iter]);

            h_bkg_mass2[iter]->SetFillColor(color);
            h_bkg_mass2[iter]->SetLineColor(color);
            h_bkg_mass2[iter]->SetDirectory(0);
            s_mass2->Add(h_bkg_mass2[iter]);

            iter++;

            if (pr == _MuMu_WJets_Full)
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if (pr == _MuMu_tW)
                pr = _MuMu_VVnST; // next - ttbar
            if (pr == _MuMu_DYTauTau_Full) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_before_PUCorr, *h_DY_mass_before_RoccoR, *h_DY_mass_before_EffCorr, *h_DY_mass_before_PVzCorr,
             *h_DY_mass_before_L1Corr, *h_DY_mass_before_TopPtCorr, *h_DY_mass_fine, *h_DY_mass, *h_DY_mass2;

        f_DY->GetObject("h_mass_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_PUCorr);
        f_DY->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_RoccoR);
        f_DY->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_EffCorr);
        f_DY->GetObject("h_mass_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_PVzCorr);
        f_DY->GetObject("h_mass_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_L1Corr);
        f_DY->GetObject("h_mass_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_TopPtCorr);
        f_DY->GetObject("h_mass_fine_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_fine);
        f_DY->GetObject("h_mass_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass);
        f_DY->GetObject("h_mass2_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass2);
        removeNegativeBins(h_DY_mass_before_PUCorr);
        removeNegativeBins(h_DY_mass_before_RoccoR);
        removeNegativeBins(h_DY_mass_before_EffCorr);
        removeNegativeBins(h_DY_mass_before_PVzCorr);
        removeNegativeBins(h_DY_mass_before_L1Corr);
        removeNegativeBins(h_DY_mass_before_TopPtCorr);
        removeNegativeBins(h_DY_mass_fine);
        removeNegativeBins(h_DY_mass);
        removeNegativeBins(h_DY_mass2);


        h_DY_mass_before_PUCorr->SetFillColor(kOrange);
        h_DY_mass_before_PUCorr->SetLineColor(kOrange);
        h_DY_mass_before_PUCorr->SetDirectory(0);
        s_mass_before_PUCorr->Add(h_DY_mass_before_PUCorr);

        h_DY_mass_before_RoccoR->SetFillColor(kOrange);
        h_DY_mass_before_RoccoR->SetLineColor(kOrange);
        h_DY_mass_before_RoccoR->SetDirectory(0);
        s_mass_before_RoccoR->Add(h_DY_mass_before_RoccoR);

        h_DY_mass_before_EffCorr->SetFillColor(kOrange);
        h_DY_mass_before_EffCorr->SetLineColor(kOrange);
        h_DY_mass_before_EffCorr->SetDirectory(0);
        s_mass_before_EffCorr->Add(h_DY_mass_before_EffCorr);

        h_DY_mass_before_PVzCorr->SetFillColor(kOrange);
        h_DY_mass_before_PVzCorr->SetLineColor(kOrange);
        h_DY_mass_before_PVzCorr->SetDirectory(0);
        s_mass_before_PVzCorr->Add(h_DY_mass_before_PVzCorr);

        h_DY_mass_before_L1Corr->SetFillColor(kOrange);
        h_DY_mass_before_L1Corr->SetLineColor(kOrange);
        h_DY_mass_before_L1Corr->SetDirectory(0);
        s_mass_before_L1Corr->Add(h_DY_mass_before_L1Corr);

        h_DY_mass_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_mass_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_mass_before_TopPtCorr->SetDirectory(0);
        s_mass_before_TopPtCorr->Add(h_DY_mass_before_TopPtCorr);

        h_DY_mass_fine->SetFillColor(kOrange);
        h_DY_mass_fine->SetLineColor(kOrange);
        h_DY_mass_fine->SetDirectory(0);
        s_mass_fine->Add(h_DY_mass_fine);

        h_DY_mass->SetFillColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);
        h_DY_mass->SetDirectory(0);
        s_mass->Add(h_DY_mass);

        h_DY_mass2->SetFillColor(kOrange);
        h_DY_mass2->SetLineColor(kOrange);
        h_DY_mass2->SetDirectory(0);
        s_mass2->Add(h_DY_mass2);

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_before_RoccoR, *h_data_mass_fine, *h_data_mass, *h_data_mass2;

        f_data->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_before_RoccoR);
        f_data->GetObject("h_mass_fine_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_fine);
        f_data->GetObject("h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass);
        f_data->GetObject("h_mass2_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass2);

        h_data_mass_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_mass_before_RoccoR->SetMarkerColor(kBlack);
        h_data_mass_before_RoccoR->SetLineColor(kBlack);
        h_data_mass_before_RoccoR->SetDirectory(0);

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass_fine->SetDirectory(0);

        h_data_mass->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerColor(kBlack);
        h_data_mass->SetLineColor(kBlack);
        h_data_mass->SetDirectory(0);

        h_data_mass2->SetMarkerStyle(kFullDotLarge);
        h_data_mass2->SetMarkerColor(kBlack);
        h_data_mass2->SetLineColor(kBlack);
        h_data_mass2->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_before_PUCorr, *RP_mass_before_RoccoR, *RP_mass_before_EffCorr, *RP_mass_before_PVzCorr,
                      *RP_mass_before_L1Corr, *RP_mass_before_TopPtCorr, *RP_mass_fine, *RP_mass, *RP_mass2;
        RP_mass_before_PUCorr = new myRatioPlot_t("RP_mass_before_PUCorr", s_mass_before_PUCorr, h_data_mass_before_RoccoR);
        RP_mass_before_RoccoR = new myRatioPlot_t("RP_mass_before_RoccoR", s_mass_before_RoccoR, h_data_mass_before_RoccoR);
        RP_mass_before_EffCorr = new myRatioPlot_t("RP_mass_before_EffCorr", s_mass_before_EffCorr, h_data_mass);
        RP_mass_before_PVzCorr = new myRatioPlot_t("RP_mass_before_PVzCorr", s_mass_before_PVzCorr, h_data_mass);
        RP_mass_before_L1Corr = new myRatioPlot_t("RP_mass_before_L1Corr", s_mass_before_L1Corr, h_data_mass);
        RP_mass_before_TopPtCorr = new myRatioPlot_t("RP_mass_before_TopPtCorr", s_mass_before_TopPtCorr, h_data_mass);
        RP_mass_fine = new myRatioPlot_t("RP_mass_fine", s_mass_fine, h_data_mass_fine);
        RP_mass = new myRatioPlot_t("RP_mass", s_mass, h_data_mass);
        RP_mass2 = new myRatioPlot_t("RP_mass2", s_mass2, h_data_mass2);

        RP_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before PU correction", 15, 3000);
//        RP_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] prie#check{s} visas pataisas", 15, 3000, "Eksp./MC");
        RP_mass_before_RoccoR->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before Rochester correction", 15, 3000);
        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before Efficiency SF", 15, 3000);
//        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] prie#check{s} efektyvumo pataisas", 15, 3000, "Eksp./MC");
        RP_mass_before_PVzCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before PVz correction", 15, 3000);
        RP_mass_before_L1Corr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before L1 prefiring correction", 15, 3000);
        RP_mass_before_TopPtCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before top p_{#lower[-0.2]{T}} reweighting", 15, 3000);
        RP_mass_fine->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 60, 120);
        RP_mass->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
//        RP_mass->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] pritaikius pataisas", 15, 3000, "Eksp./MC");
        RP_mass2->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

        TLegend *legend = new TLegend(0.8, 0.45, 0.95, 0.95);

        legend->AddEntry(h_data_mass, "Data", "lp");
//        legend->AddEntry(h_data_mass, "Matavimas", "lp");
        legend->AddEntry(h_DY_mass, "DY#rightarrow#mu#mu", "f");
        legend->AddEntry(h_bkg_mass[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_mass[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_mass[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_mass[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_mass[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_mass[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_mass[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_mass[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
        legend->AddEntry(h_bkg_mass[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_mass_before_PUCorr->ImportLegend(legend);
        RP_mass_before_RoccoR->ImportLegend(legend);
        RP_mass_before_EffCorr->ImportLegend(legend);
        RP_mass_before_PVzCorr->ImportLegend(legend);
        RP_mass_before_L1Corr->ImportLegend(legend);
        RP_mass_before_TopPtCorr->ImportLegend(legend);
        RP_mass_fine->ImportLegend(legend);
        RP_mass->ImportLegend(legend);
        RP_mass2->ImportLegend(legend);

        RP_mass_before_PUCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_RoccoR->Draw(0.5, 1e7, 1);
        RP_mass_before_EffCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_PVzCorr->Draw(0.5, 1e7, 1);
        RP_mass_before_L1Corr->Draw(0.5, 1e7, 1);
        RP_mass_before_TopPtCorr->Draw(0.5, 1e7, 1);
        RP_mass_fine->Draw(0.5, 1e7, 0);
        RP_mass->Draw(0.5, 1e7, 1);
        RP_mass2->Draw(0.5, 1e7, 1);

        Double_t dataerror, MCerror, MCerror_noSF, DYerror, dataintegral=2.25081e+07, MCintegral, MCintegral_noSF, DYintegral;
        Double_t dataerrorZ, MCerrorZ, DYerrorZ, dataintegralZ=2.25081e+07, MCintegralZ, DYintegralZ;
        Double_t dataerror_noZ=0, MCerror_noZ=0, DYerror_noZ=0, dataintegral_noZ=2.25081e+07, MCintegral_noZ, DYintegral_noZ, temp_noZ;

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);
        MCintegral_noSF = ((TH1D*)(s_mass_before_RoccoR->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror_noSF);
        DYintegral = h_DY_mass->IntegralAndError(1, h_DY_mass->GetSize()-2, DYerror);

        dataintegralZ = h_data_mass->IntegralAndError(10, 22, dataerrorZ);
        MCintegralZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ);
        DYintegralZ = h_DY_mass->IntegralAndError(10, 22, DYerrorZ);

        dataintegral_noZ = h_data_mass->IntegralAndError(1, 9, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataintegral_noZ += h_data_mass->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataerror_noZ = sqrt(dataerror_noZ);

        MCintegral_noZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCintegral_noZ += ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCerror_noZ = sqrt(MCerror_noZ);

        DYintegral_noZ = h_DY_mass->IntegralAndError(1, 9, temp_noZ);
        DYerror_noZ += temp_noZ * temp_noZ;
        DYintegral_noZ += h_DY_mass->IntegralAndError(23, h_DY_mass->GetSize()-2, temp_noZ);
        DYerror_noZ += temp_noZ * temp_noZ;
        DYerror_noZ = sqrt(DYerror_noZ);

        std::cout << "Data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;
        std::cout << "MC/Obs: " << MCintegral / dataintegral << "+-" <<
                     sqrt((dataerror / dataintegral) * (dataerror / dataintegral) +
                           (MCerror / MCintegral) * (MCerror / MCintegral)) << endl;
        std::cout << "MC events before corrections: " << MCintegral_noSF << "+-" << MCerror_noSF << endl;
        std::cout << "DY events: " << DYintegral << "+-" << DYerror << endl;
        std::cout << "Avg. Data and MC relative difference: " << CompAvgDataMCDifference(h_data_mass, s_mass) << endl;
        std::cout << "Chi^2: " << CompChiSquared(h_data_mass, s_mass) << endl << endl;

        std::cout << "Data events around Z: " << dataintegralZ << "+-" << dataerrorZ << endl;
        std::cout << "MC events around Z: " << MCintegralZ << "+-" << MCerrorZ << endl;
        std::cout << "DY events around Z: " << DYintegralZ << "+-" << DYerrorZ << endl << endl;
        std::cout << "Data events outside Z: " << dataintegral_noZ << "+-" << dataerror_noZ << endl;
        std::cout << "MC events outside Z: " << MCintegral_noZ << "+-" << MCerror_noZ << endl;
        std::cout << "DY events outside Z: " << DYintegral_noZ << "+-" << DYerror_noZ << endl << endl;


    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################
/*
    if (whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI"))
    {
        count_drawn++;

        THStack *s_rapi_before_PUCorr = new THStack("s_rapi_before_PUCorr", "");
        THStack *s_pT_before_PUCorr = new THStack("s_pT_before_PUCorr", "");
        THStack *s_pT_lead_before_PUCorr = new THStack("s_pT_lead_before_PUCorr", "");
        THStack *s_pT_sublead_before_PUCorr = new THStack("s_pT_sublead_before_PUCorr", "");
        THStack *s_eta_lead_before_PUCorr = new THStack("s_eta_lead_before_PUCorr", "");
        THStack *s_eta_sublead_before_PUCorr = new THStack("s_eta_sublead_before_PUCorr", "");
        THStack *s_phi_lead_before_PUCorr = new THStack("s_phi_lead_before_PUCorr", "");
        THStack *s_phi_sublead_before_PUCorr = new THStack("s_phi_sublead_before_PUCorr", "");

        THStack *s_rapi_before_RoccoR = new THStack("s_rapi_before_RoccoR", "");
        THStack *s_pT_before_RoccoR = new THStack("s_pT_before_RoccoR", "");
        THStack *s_pT_lead_before_RoccoR = new THStack("s_pT_lead_before_RoccoR", "");
        THStack *s_pT_sublead_before_RoccoR = new THStack("s_pT_sublead_before_RoccoR", "");
        THStack *s_eta_lead_before_RoccoR = new THStack("s_eta_lead_before_RoccoR", "");
        THStack *s_eta_sublead_before_RoccoR = new THStack("s_eta_sublead_before_RoccoR", "");
        THStack *s_phi_lead_before_RoccoR = new THStack("s_phi_lead_before_RoccoR", "");
        THStack *s_phi_sublead_before_RoccoR = new THStack("s_phi_sublead_before_RoccoR", "");

        THStack *s_rapi_before_EffCorr = new THStack("s_rapi_before_EffCorr", "");
        THStack *s_pT_before_EffCorr = new THStack("s_pT_before_EffCorr", "");
        THStack *s_pT_lead_before_EffCorr = new THStack("s_pT_lead_before_EffCorr", "");
        THStack *s_pT_sublead_before_EffCorr = new THStack("s_pT_sublead_before_EffCorr", "");
        THStack *s_eta_lead_before_EffCorr = new THStack("s_eta_lead_before_EffCorr", "");
        THStack *s_eta_sublead_before_EffCorr = new THStack("s_eta_sublead_before_EffCorr", "");
        THStack *s_phi_lead_before_EffCorr = new THStack("s_phi_lead_before_EffCorr", "");
        THStack *s_phi_sublead_before_EffCorr = new THStack("s_phi_sublead_before_EffCorr", "");

        THStack *s_rapi_before_PVzCorr = new THStack("s_rapi_before_PVzCorr", "");
        THStack *s_pT_before_PVzCorr = new THStack("s_pT_before_PVzCorr", "");
        THStack *s_pT_lead_before_PVzCorr = new THStack("s_pT_lead_before_PVzCorr", "");
        THStack *s_pT_sublead_before_PVzCorr = new THStack("s_pT_sublead_before_PVzCorr", "");
        THStack *s_eta_lead_before_PVzCorr = new THStack("s_eta_lead_before_PVzCorr", "");
        THStack *s_eta_sublead_before_PVzCorr = new THStack("s_eta_sublead_before_PVzCorr", "");
        THStack *s_phi_lead_before_PVzCorr = new THStack("s_phi_lead_before_PVzCorr", "");
        THStack *s_phi_sublead_before_PVzCorr = new THStack("s_phi_sublead_before_PVzCorr", "");

        THStack *s_rapi_before_L1Corr = new THStack("s_rapi_before_L1Corr", "");
        THStack *s_pT_before_L1Corr = new THStack("s_pT_before_L1Corr", "");
        THStack *s_pT_lead_before_L1Corr = new THStack("s_pT_lead_before_L1Corr", "");
        THStack *s_pT_sublead_before_L1Corr = new THStack("s_pT_sublead_before_L1Corr", "");
        THStack *s_eta_lead_before_L1Corr = new THStack("s_eta_lead_before_L1Corr", "");
        THStack *s_eta_sublead_before_L1Corr = new THStack("s_eta_sublead_before_L1Corr", "");
        THStack *s_phi_lead_before_L1Corr = new THStack("s_phi_lead_before_L1Corr", "");
        THStack *s_phi_sublead_before_L1Corr = new THStack("s_phi_sublead_before_L1Corr", "");

        THStack *s_rapi_before_TopPtCorr = new THStack("s_rapi_before_TopPtCorr", "");
        THStack *s_pT_before_TopPtCorr = new THStack("s_pT_before_TopPtCorr", "");
        THStack *s_pT_lead_before_TopPtCorr = new THStack("s_pT_lead_before_TopPtCorr", "");
        THStack *s_pT_sublead_before_TopPtCorr = new THStack("s_pT_sublead_before_TopPtCorr", "");
        THStack *s_eta_lead_before_TopPtCorr = new THStack("s_eta_lead_before_TopPtCorr", "");
        THStack *s_eta_sublead_before_TopPtCorr = new THStack("s_eta_sublead_before_TopPtCorr", "");
        THStack *s_phi_lead_before_TopPtCorr = new THStack("s_phi_lead_before_TopPtCorr", "");
        THStack *s_phi_sublead_before_TopPtCorr = new THStack("s_phi_sublead_before_TopPtCorr", "");

        THStack *s_rapi = new THStack("s_rapi", "");
        THStack *s_pT = new THStack("s_pT", "");
        THStack *s_pT_lead = new THStack("s_pT_lead", "");
        THStack *s_pT_sublead = new THStack("s_pT_sublead", "");
        THStack *s_eta_lead = new THStack("s_eta_lead", "");
        THStack *s_eta_sublead = new THStack("s_eta_sublead", "");
        THStack *s_phi_lead = new THStack("s_phi_lead", "");
        THStack *s_phi_sublead = new THStack("s_phi_sublead", "");

//----------------------------------- MC bkg -------------------------------------------------------

        TH1D *h_bkg_pT_before_PUCorr[9], *h_bkg_pT_before_RoccoR[9], *h_bkg_pT_before_EffCorr[9],
             *h_bkg_pT_before_PVzCorr[9], *h_bkg_pT_before_L1Corr[9], *h_bkg_pT_before_TopPtCorr[9], *h_bkg_pT[9],
             *h_bkg_rapi_before_PUCorr[9], *h_bkg_rapi_before_RoccoR[9], *h_bkg_rapi_before_EffCorr[9],
             *h_bkg_rapi_before_PVzCorr[9], *h_bkg_rapi_before_L1Corr[9], *h_bkg_rapi_before_TopPtCorr[9], *h_bkg_rapi[9],
             *h_bkg_pT_lead_before_PUCorr[9], *h_bkg_pT_lead_before_RoccoR[9], *h_bkg_pT_lead_before_EffCorr[9],
             *h_bkg_pT_lead_before_PVzCorr[9], *h_bkg_pT_lead_before_L1Corr[9], *h_bkg_pT_lead_before_TopPtCorr[9], *h_bkg_pT_lead[9],
             *h_bkg_pT_sublead_before_PUCorr[9],*h_bkg_pT_sublead_before_RoccoR[9],  *h_bkg_pT_sublead_before_EffCorr[9],
             *h_bkg_pT_sublead_before_PVzCorr[9], *h_bkg_pT_sublead_before_L1Corr[9], *h_bkg_pT_sublead_before_TopPtCorr[9], *h_bkg_pT_sublead[9],
             *h_bkg_eta_lead_before_PUCorr[9], *h_bkg_eta_lead_before_RoccoR[9], *h_bkg_eta_lead_before_EffCorr[9],
             *h_bkg_eta_lead_before_PVzCorr[9], *h_bkg_eta_lead_before_L1Corr[9], *h_bkg_eta_lead_before_TopPtCorr[9], *h_bkg_eta_lead[9],
             *h_bkg_eta_sublead_before_PUCorr[9], *h_bkg_eta_sublead_before_RoccoR[9], *h_bkg_eta_sublead_before_EffCorr[9],
             *h_bkg_eta_sublead_before_PVzCorr[9], *h_bkg_eta_sublead_before_L1Corr[9], *h_bkg_eta_sublead_before_TopPtCorr[9], *h_bkg_eta_sublead[9],
             *h_bkg_phi_lead_before_PUCorr[9], *h_bkg_phi_lead_before_RoccoR[9], *h_bkg_phi_lead_before_EffCorr[9],
             *h_bkg_phi_lead_before_PVzCorr[9], *h_bkg_phi_lead_before_L1Corr[9], *h_bkg_phi_lead_before_TopPtCorr[9], *h_bkg_phi_lead[9],
             *h_bkg_phi_sublead_before_PUCorr[9], *h_bkg_phi_sublead_before_RoccoR[9], *h_bkg_phi_sublead_before_EffCorr[9],
             *h_bkg_phi_sublead_before_PVzCorr[9], *h_bkg_phi_sublead_before_L1Corr[9], *h_bkg_phi_sublead_before_TopPtCorr[9], *h_bkg_phi_sublead[9];
        Int_t iter = 0;

        for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _MuMu_QCDMuEnriched_Full)
            {
                iter++;
                continue;
            }
            f_bkg->GetObject("h_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_before_PUCorr[iter]);
            f_bkg->GetObject("h_rapi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_PUCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PUCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PUCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PUCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_PUCorr[iter]);

            f_bkg->GetObject("h_pT_before_RocCorr_"+Mgr.Procname[pr], h_bkg_pT_before_RoccoR[iter]);
            f_bkg->GetObject("h_rapi_before_RocCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_RoccoR[iter]);
            f_bkg->GetObject("h_pT_lead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_RoccoR[iter]);
            f_bkg->GetObject("h_pT_sublead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_RoccoR[iter]);
            f_bkg->GetObject("h_eta_lead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_RoccoR[iter]);
            f_bkg->GetObject("h_eta_sublead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_RoccoR[iter]);
            f_bkg->GetObject("h_phi_lead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_RoccoR[iter]);
            f_bkg->GetObject("h_phi_sublead_before_RocCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_pT_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_rapi_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_RoccoR[iter]);

            f_bkg->GetObject("h_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_before_EffCorr[iter]);
            f_bkg->GetObject("h_rapi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_EffCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_EffCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_EffCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_EffCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_EffCorr[iter]);

            f_bkg->GetObject("h_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_before_PVzCorr[iter]);
            f_bkg->GetObject("h_rapi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_PVzCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_PVzCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_PVzCorr[iter]);

            f_bkg->GetObject("h_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_before_L1Corr[iter]);
            f_bkg->GetObject("h_rapi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_rapi_before_L1Corr[iter]);
            f_bkg->GetObject("h_pT_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_L1Corr[iter]);
            f_bkg->GetObject("h_eta_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_L1Corr[iter]);
            f_bkg->GetObject("h_phi_lead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_L1Corr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_rapi_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_L1Corr[iter]);

            f_bkg->GetObject("h_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_rapi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_rapi_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_pT_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_pT_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_sublead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eta_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eta_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_sublead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_phi_lead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_lead_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_phi_sublead_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_rapi_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_sublead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_lead_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_sublead_before_TopPtCorr[iter]);

            f_bkg->GetObject("h_pT_"+Mgr.Procname[pr], h_bkg_pT[iter]);
            f_bkg->GetObject("h_rapi_"+Mgr.Procname[pr], h_bkg_rapi[iter]);
            f_bkg->GetObject("h_pT_lead_"+Mgr.Procname[pr], h_bkg_pT_lead[iter]);
            f_bkg->GetObject("h_pT_sublead_"+Mgr.Procname[pr], h_bkg_pT_sublead[iter]);
            f_bkg->GetObject("h_eta_lead_"+Mgr.Procname[pr], h_bkg_eta_lead[iter]);
            f_bkg->GetObject("h_eta_sublead_"+Mgr.Procname[pr], h_bkg_eta_sublead[iter]);
            f_bkg->GetObject("h_phi_lead_"+Mgr.Procname[pr], h_bkg_phi_lead[iter]);
            f_bkg->GetObject("h_phi_sublead_"+Mgr.Procname[pr], h_bkg_phi_sublead[iter]);
            removeNegativeBins(h_bkg_pT[iter]);
            removeNegativeBins(h_bkg_rapi[iter]);
            removeNegativeBins(h_bkg_pT_lead[iter]);
            removeNegativeBins(h_bkg_pT_sublead[iter]);
            removeNegativeBins(h_bkg_eta_lead[iter]);
            removeNegativeBins(h_bkg_eta_sublead[iter]);
            removeNegativeBins(h_bkg_phi_lead[iter]);
            removeNegativeBins(h_bkg_phi_sublead[iter]);

            Color_t color = kBlack;
            if (pr == _MuMu_QCDMuEnriched_Full) color = kRed + 3;
            if (pr == _MuMu_WJets_Full) color = kRed - 2;
            if (pr == _MuMu_WW) color = kMagenta - 5;
            if (pr == _MuMu_WZ) color = kMagenta - 2;
            if (pr == _MuMu_ZZ) color = kMagenta - 6;
            if (pr == _MuMu_tbarW) color = kGreen - 2;
            if (pr == _MuMu_tW) color = kGreen + 2;
            if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
            if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_pT_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PUCorr[iter]->SetDirectory(0);
            //
            s_pT_before_PUCorr->Add(h_bkg_pT_before_PUCorr[iter]);
            s_rapi_before_PUCorr->Add(h_bkg_rapi_before_PUCorr[iter]);
            s_pT_lead_before_PUCorr->Add(h_bkg_pT_lead_before_PUCorr[iter]);
            s_pT_sublead_before_PUCorr->Add(h_bkg_pT_sublead_before_PUCorr[iter]);
            s_eta_lead_before_PUCorr->Add(h_bkg_eta_lead_before_PUCorr[iter]);
            s_eta_sublead_before_PUCorr->Add(h_bkg_eta_sublead_before_PUCorr[iter]);
            s_phi_lead_before_PUCorr->Add(h_bkg_phi_lead_before_PUCorr[iter]);
            s_phi_sublead_before_PUCorr->Add(h_bkg_phi_sublead_before_PUCorr[iter]);

            h_bkg_pT_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_rapi_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_rapi_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_rapi_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_RoccoR[iter]->SetDirectory(0);
            //
            s_pT_before_RoccoR->Add(h_bkg_pT_before_RoccoR[iter]);
            s_rapi_before_RoccoR->Add(h_bkg_rapi_before_RoccoR[iter]);
            s_pT_lead_before_RoccoR->Add(h_bkg_pT_lead_before_RoccoR[iter]);
            s_pT_sublead_before_RoccoR->Add(h_bkg_pT_sublead_before_RoccoR[iter]);
            s_eta_lead_before_RoccoR->Add(h_bkg_eta_lead_before_RoccoR[iter]);
            s_eta_sublead_before_RoccoR->Add(h_bkg_eta_sublead_before_RoccoR[iter]);
            s_phi_lead_before_RoccoR->Add(h_bkg_phi_lead_before_RoccoR[iter]);
            s_phi_sublead_before_RoccoR->Add(h_bkg_phi_sublead_before_RoccoR[iter]);

            h_bkg_pT_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_EffCorr[iter]->SetDirectory(0);
            //
            s_pT_before_EffCorr->Add(h_bkg_pT_before_EffCorr[iter]);
            s_rapi_before_EffCorr->Add(h_bkg_rapi_before_EffCorr[iter]);
            s_pT_lead_before_EffCorr->Add(h_bkg_pT_lead_before_EffCorr[iter]);
            s_pT_sublead_before_EffCorr->Add(h_bkg_pT_sublead_before_EffCorr[iter]);
            s_eta_lead_before_EffCorr->Add(h_bkg_eta_lead_before_EffCorr[iter]);
            s_eta_sublead_before_EffCorr->Add(h_bkg_eta_sublead_before_EffCorr[iter]);
            s_phi_lead_before_EffCorr->Add(h_bkg_phi_lead_before_EffCorr[iter]);
            s_phi_sublead_before_EffCorr->Add(h_bkg_phi_sublead_before_EffCorr[iter]);

            h_bkg_pT_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_PVzCorr[iter]->SetDirectory(0);
            //
            s_pT_before_PVzCorr->Add(h_bkg_pT_before_PVzCorr[iter]);
            s_rapi_before_PVzCorr->Add(h_bkg_rapi_before_PVzCorr[iter]);
            s_pT_lead_before_PVzCorr->Add(h_bkg_pT_lead_before_PVzCorr[iter]);
            s_pT_sublead_before_PVzCorr->Add(h_bkg_pT_sublead_before_PVzCorr[iter]);
            s_eta_lead_before_PVzCorr->Add(h_bkg_eta_lead_before_PVzCorr[iter]);
            s_eta_sublead_before_PVzCorr->Add(h_bkg_eta_sublead_before_PVzCorr[iter]);
            s_phi_lead_before_PVzCorr->Add(h_bkg_phi_lead_before_PVzCorr[iter]);
            s_phi_sublead_before_PVzCorr->Add(h_bkg_phi_sublead_before_PVzCorr[iter]);

            h_bkg_pT_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_rapi_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_rapi_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_rapi_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_L1Corr[iter]->SetDirectory(0);
            //
            s_pT_before_L1Corr->Add(h_bkg_pT_before_L1Corr[iter]);
            s_rapi_before_L1Corr->Add(h_bkg_rapi_before_L1Corr[iter]);
            s_pT_lead_before_L1Corr->Add(h_bkg_pT_lead_before_L1Corr[iter]);
            s_pT_sublead_before_L1Corr->Add(h_bkg_pT_sublead_before_L1Corr[iter]);
            s_eta_lead_before_L1Corr->Add(h_bkg_eta_lead_before_L1Corr[iter]);
            s_eta_sublead_before_L1Corr->Add(h_bkg_eta_sublead_before_L1Corr[iter]);
            s_phi_lead_before_L1Corr->Add(h_bkg_phi_lead_before_L1Corr[iter]);
            s_phi_sublead_before_L1Corr->Add(h_bkg_phi_sublead_before_L1Corr[iter]);

            h_bkg_pT_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_rapi_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_rapi_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_rapi_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_lead_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_sublead_before_TopPtCorr[iter]->SetDirectory(0);
            //
            s_pT_before_TopPtCorr->Add(h_bkg_pT_before_TopPtCorr[iter]);
            s_rapi_before_TopPtCorr->Add(h_bkg_rapi_before_TopPtCorr[iter]);
            s_pT_lead_before_TopPtCorr->Add(h_bkg_pT_lead_before_TopPtCorr[iter]);
            s_pT_sublead_before_TopPtCorr->Add(h_bkg_pT_sublead_before_TopPtCorr[iter]);
            s_eta_lead_before_TopPtCorr->Add(h_bkg_eta_lead_before_TopPtCorr[iter]);
            s_eta_sublead_before_TopPtCorr->Add(h_bkg_eta_sublead_before_TopPtCorr[iter]);
            s_phi_lead_before_TopPtCorr->Add(h_bkg_phi_lead_before_TopPtCorr[iter]);
            s_phi_sublead_before_TopPtCorr->Add(h_bkg_phi_sublead_before_TopPtCorr[iter]);

            h_bkg_pT[iter]->SetFillColor(color);
            h_bkg_rapi[iter]->SetFillColor(color);
            h_bkg_pT_lead[iter]->SetFillColor(color);
            h_bkg_pT_sublead[iter]->SetFillColor(color);
            h_bkg_eta_lead[iter]->SetFillColor(color);
            h_bkg_eta_sublead[iter]->SetFillColor(color);
            h_bkg_phi_lead[iter]->SetFillColor(color);
            h_bkg_phi_sublead[iter]->SetFillColor(color);
            //
            h_bkg_pT[iter]->SetLineColor(color);
            h_bkg_rapi[iter]->SetLineColor(color);
            h_bkg_pT_lead[iter]->SetLineColor(color);
            h_bkg_pT_sublead[iter]->SetLineColor(color);
            h_bkg_eta_lead[iter]->SetLineColor(color);
            h_bkg_eta_sublead[iter]->SetLineColor(color);
            h_bkg_phi_lead[iter]->SetLineColor(color);
            h_bkg_phi_sublead[iter]->SetLineColor(color);
            //
            h_bkg_pT[iter]->SetDirectory(0);
            h_bkg_rapi[iter]->SetDirectory(0);
            h_bkg_pT_lead[iter]->SetDirectory(0);
            h_bkg_pT_sublead[iter]->SetDirectory(0);
            h_bkg_eta_lead[iter]->SetDirectory(0);
            h_bkg_eta_sublead[iter]->SetDirectory(0);
            h_bkg_phi_lead[iter]->SetDirectory(0);
            h_bkg_phi_sublead[iter]->SetDirectory(0);
            //
            s_pT->Add(h_bkg_pT[iter]);
            s_rapi->Add(h_bkg_rapi[iter]);
            s_pT_lead->Add(h_bkg_pT_lead[iter]);
            s_pT_sublead->Add(h_bkg_pT_sublead[iter]);
            s_eta_lead->Add(h_bkg_eta_lead[iter]);
            s_eta_sublead->Add(h_bkg_eta_sublead[iter]);
            s_phi_lead_before_L1Corr->Add(h_bkg_phi_lead_before_L1Corr[iter]);
            s_phi_lead->Add(h_bkg_phi_lead[iter]);
            s_phi_sublead->Add(h_bkg_phi_sublead[iter]);

            iter++;

            if (pr == _MuMu_WJets_Full)
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if (pr == _MuMu_tW)
                pr = _MuMu_VVnST; // next - ttbar
            if (pr == _MuMu_DYTauTau_Full) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_pT_before_PUCorr, *h_DY_pT_before_RoccoR, *h_DY_pT_before_EffCorr,
             *h_DY_pT_before_PVzCorr, *h_DY_pT_before_L1Corr, *h_DY_pT_before_TopPtCorr, *h_DY_pT,
             *h_DY_rapi_before_PUCorr, *h_DY_rapi_before_RoccoR, *h_DY_rapi_before_EffCorr,
             *h_DY_rapi_before_PVzCorr, *h_DY_rapi_before_L1Corr, *h_DY_rapi_before_TopPtCorr, *h_DY_rapi,
             *h_DY_pT_lead_before_PUCorr, *h_DY_pT_lead_before_RoccoR, *h_DY_pT_lead_before_EffCorr,
             *h_DY_pT_lead_before_PVzCorr, *h_DY_pT_lead_before_L1Corr, *h_DY_pT_lead_before_TopPtCorr, *h_DY_pT_lead,
             *h_DY_pT_sublead_before_PUCorr, *h_DY_pT_sublead_before_RoccoR, *h_DY_pT_sublead_before_EffCorr,
             *h_DY_pT_sublead_before_PVzCorr, *h_DY_pT_sublead_before_L1Corr, *h_DY_pT_sublead_before_TopPtCorr, *h_DY_pT_sublead,
             *h_DY_eta_lead_before_PUCorr, *h_DY_eta_lead_before_RoccoR, *h_DY_eta_lead_before_EffCorr,
             *h_DY_eta_lead_before_PVzCorr, *h_DY_eta_lead_before_L1Corr, *h_DY_eta_lead_before_TopPtCorr, *h_DY_eta_lead,
             *h_DY_eta_sublead_before_PUCorr, *h_DY_eta_sublead_before_RoccoR, *h_DY_eta_sublead_before_EffCorr,
             *h_DY_eta_sublead_before_PVzCorr, *h_DY_eta_sublead_before_L1Corr, *h_DY_eta_sublead_before_TopPtCorr, *h_DY_eta_sublead,
             *h_DY_phi_lead_before_PUCorr, *h_DY_phi_lead_before_RoccoR, *h_DY_phi_lead_before_EffCorr,
             *h_DY_phi_lead_before_PVzCorr, *h_DY_phi_lead_before_L1Corr, *h_DY_phi_lead_before_TopPtCorr, *h_DY_phi_lead,
             *h_DY_phi_sublead_before_PUCorr, *h_DY_phi_sublead_before_RoccoR, *h_DY_phi_sublead_before_EffCorr,
             *h_DY_phi_sublead_before_PVzCorr, *h_DY_phi_sublead_before_L1Corr, *h_DY_phi_sublead_before_TopPtCorr, *h_DY_phi_sublead;

        f_DY->GetObject("h_pT_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_PUCorr);
        f_DY->GetObject("h_rapi_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_PUCorr);
        f_DY->GetObject("h_pT_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_PUCorr);
        f_DY->GetObject("h_pT_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_PUCorr);
        f_DY->GetObject("h_eta_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_PUCorr);
        f_DY->GetObject("h_eta_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_PUCorr);
        f_DY->GetObject("h_phi_lead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_PUCorr);
        f_DY->GetObject("h_phi_sublead_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_PUCorr);
        removeNegativeBins(h_DY_pT_before_PUCorr);
        removeNegativeBins(h_DY_rapi_before_PUCorr);
        removeNegativeBins(h_DY_pT_lead_before_PUCorr);
        removeNegativeBins(h_DY_pT_sublead_before_PUCorr);
        removeNegativeBins(h_DY_eta_lead_before_PUCorr);
        removeNegativeBins(h_DY_eta_sublead_before_PUCorr);
        removeNegativeBins(h_DY_phi_lead_before_PUCorr);
        removeNegativeBins(h_DY_phi_sublead_before_PUCorr);

        f_DY->GetObject("h_pT_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_RoccoR);
        f_DY->GetObject("h_rapi_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_RoccoR);
        f_DY->GetObject("h_pT_lead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_RoccoR);
        f_DY->GetObject("h_pT_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_RoccoR);
        f_DY->GetObject("h_eta_lead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_RoccoR);
        f_DY->GetObject("h_eta_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_RoccoR);
        f_DY->GetObject("h_phi_lead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_RoccoR);
        f_DY->GetObject("h_phi_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_RoccoR);
        removeNegativeBins(h_DY_pT_before_RoccoR);
        removeNegativeBins(h_DY_rapi_before_RoccoR);
        removeNegativeBins(h_DY_pT_lead_before_RoccoR);
        removeNegativeBins(h_DY_pT_sublead_before_RoccoR);
        removeNegativeBins(h_DY_eta_lead_before_RoccoR);
        removeNegativeBins(h_DY_eta_sublead_before_RoccoR);
        removeNegativeBins(h_DY_phi_lead_before_RoccoR);
        removeNegativeBins(h_DY_phi_sublead_before_RoccoR);

        f_DY->GetObject("h_pT_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_EffCorr);
        f_DY->GetObject("h_rapi_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_EffCorr);
        f_DY->GetObject("h_pT_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_EffCorr);
        f_DY->GetObject("h_pT_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_EffCorr);
        f_DY->GetObject("h_eta_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_EffCorr);
        f_DY->GetObject("h_eta_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_EffCorr);
        f_DY->GetObject("h_phi_lead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_EffCorr);
        f_DY->GetObject("h_phi_sublead_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_EffCorr);
        removeNegativeBins(h_DY_pT_before_EffCorr);
        removeNegativeBins(h_DY_rapi_before_EffCorr);
        removeNegativeBins(h_DY_pT_lead_before_EffCorr);
        removeNegativeBins(h_DY_pT_sublead_before_EffCorr);
        removeNegativeBins(h_DY_eta_lead_before_EffCorr);
        removeNegativeBins(h_DY_eta_sublead_before_EffCorr);
        removeNegativeBins(h_DY_phi_lead_before_EffCorr);
        removeNegativeBins(h_DY_phi_sublead_before_EffCorr);

        f_DY->GetObject("h_pT_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_PVzCorr);
        f_DY->GetObject("h_rapi_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_PVzCorr);
        f_DY->GetObject("h_pT_lead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_PVzCorr);
        f_DY->GetObject("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_PVzCorr);
        f_DY->GetObject("h_eta_lead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_PVzCorr);
        f_DY->GetObject("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_PVzCorr);
        f_DY->GetObject("h_phi_lead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_PVzCorr);
        f_DY->GetObject("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_pT_before_PVzCorr);
        removeNegativeBins(h_DY_rapi_before_PVzCorr);
        removeNegativeBins(h_DY_pT_lead_before_PVzCorr);
        removeNegativeBins(h_DY_pT_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_eta_lead_before_PVzCorr);
        removeNegativeBins(h_DY_eta_sublead_before_PVzCorr);
        removeNegativeBins(h_DY_phi_lead_before_PVzCorr);
        removeNegativeBins(h_DY_phi_sublead_before_PVzCorr);

        f_DY->GetObject("h_pT_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_L1Corr);
        f_DY->GetObject("h_rapi_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_L1Corr);
        f_DY->GetObject("h_pT_lead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_L1Corr);
        f_DY->GetObject("h_pT_sublead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_L1Corr);
        f_DY->GetObject("h_eta_lead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_L1Corr);
        f_DY->GetObject("h_eta_sublead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_L1Corr);
        f_DY->GetObject("h_phi_lead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_L1Corr);
        f_DY->GetObject("h_phi_sublead_before_L1Corr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_L1Corr);
        removeNegativeBins(h_DY_pT_before_L1Corr);
        removeNegativeBins(h_DY_rapi_before_L1Corr);
        removeNegativeBins(h_DY_pT_lead_before_L1Corr);
        removeNegativeBins(h_DY_pT_sublead_before_L1Corr);
        removeNegativeBins(h_DY_eta_lead_before_L1Corr);
        removeNegativeBins(h_DY_eta_sublead_before_L1Corr);
        removeNegativeBins(h_DY_phi_lead_before_L1Corr);
        removeNegativeBins(h_DY_phi_sublead_before_L1Corr);

        f_DY->GetObject("h_pT_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_before_TopPtCorr);
        f_DY->GetObject("h_rapi_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi_before_TopPtCorr);
        f_DY->GetObject("h_pT_lead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead_before_TopPtCorr);
        f_DY->GetObject("h_pT_sublead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead_before_TopPtCorr);
        f_DY->GetObject("h_eta_lead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead_before_TopPtCorr);
        f_DY->GetObject("h_eta_sublead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead_before_TopPtCorr);
        f_DY->GetObject("h_phi_lead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead_before_TopPtCorr);
        f_DY->GetObject("h_phi_sublead_before_TopPtCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_before_TopPtCorr);
        removeNegativeBins(h_DY_rapi_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_pT_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_eta_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_eta_sublead_before_TopPtCorr);
        removeNegativeBins(h_DY_phi_lead_before_TopPtCorr);
        removeNegativeBins(h_DY_phi_sublead_before_TopPtCorr);

        f_DY->GetObject("h_pT_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT);
        f_DY->GetObject("h_rapi_"+Mgr.Procname[_MuMu_DY_Full], h_DY_rapi);
        f_DY->GetObject("h_pT_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_lead);
        f_DY->GetObject("h_pT_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_pT_sublead);
        f_DY->GetObject("h_eta_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_lead);
        f_DY->GetObject("h_eta_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_eta_sublead);
        f_DY->GetObject("h_phi_lead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_lead);
        f_DY->GetObject("h_phi_sublead_"+Mgr.Procname[_MuMu_DY_Full], h_DY_phi_sublead);
        removeNegativeBins(h_DY_pT);
        removeNegativeBins(h_DY_rapi);
        removeNegativeBins(h_DY_pT_lead);
        removeNegativeBins(h_DY_pT_sublead);
        removeNegativeBins(h_DY_eta_lead);
        removeNegativeBins(h_DY_eta_sublead);
        removeNegativeBins(h_DY_phi_lead);
        removeNegativeBins(h_DY_phi_sublead);

        h_DY_pT_before_PUCorr->SetFillColor(kOrange);
        h_DY_rapi_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_PUCorr->SetLineColor(kOrange);
        h_DY_rapi_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_PUCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PUCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_PUCorr->SetDirectory(0);
        h_DY_rapi_before_PUCorr->SetDirectory(0);
        h_DY_pT_lead_before_PUCorr->SetDirectory(0);
        h_DY_pT_sublead_before_PUCorr->SetDirectory(0);
        h_DY_eta_lead_before_PUCorr->SetDirectory(0);
        h_DY_eta_sublead_before_PUCorr->SetDirectory(0);
        h_DY_phi_lead_before_PUCorr->SetDirectory(0);
        h_DY_phi_sublead_before_PUCorr->SetDirectory(0);
        //
        s_pT_before_PUCorr->Add(h_DY_pT_before_PUCorr);
        s_rapi_before_PUCorr->Add(h_DY_rapi_before_PUCorr);
        s_pT_lead_before_PUCorr->Add(h_DY_pT_lead_before_PUCorr);
        s_pT_sublead_before_PUCorr->Add(h_DY_pT_sublead_before_PUCorr);
        s_eta_lead_before_PUCorr->Add(h_DY_eta_lead_before_PUCorr);
        s_eta_sublead_before_PUCorr->Add(h_DY_eta_sublead_before_PUCorr);
        s_phi_lead_before_PUCorr->Add(h_DY_phi_lead_before_PUCorr);
        s_phi_sublead_before_PUCorr->Add(h_DY_phi_sublead_before_PUCorr);

        h_DY_pT_before_RoccoR->SetFillColor(kOrange);
        h_DY_rapi_before_RoccoR->SetFillColor(kOrange);
        h_DY_pT_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_pT_sublead_before_RoccoR->SetFillColor(kOrange);
        h_DY_eta_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_eta_sublead_before_RoccoR->SetFillColor(kOrange);
        h_DY_phi_lead_before_RoccoR->SetFillColor(kOrange);
        h_DY_phi_sublead_before_RoccoR->SetFillColor(kOrange);
        //
        h_DY_pT_before_RoccoR->SetLineColor(kOrange);
        h_DY_rapi_before_RoccoR->SetLineColor(kOrange);
        h_DY_pT_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_pT_sublead_before_RoccoR->SetLineColor(kOrange);
        h_DY_eta_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_eta_sublead_before_RoccoR->SetLineColor(kOrange);
        h_DY_phi_lead_before_RoccoR->SetLineColor(kOrange);
        h_DY_phi_sublead_before_RoccoR->SetLineColor(kOrange);
        //
        h_DY_pT_before_RoccoR->SetDirectory(0);
        h_DY_rapi_before_RoccoR->SetDirectory(0);
        h_DY_pT_lead_before_RoccoR->SetDirectory(0);
        h_DY_pT_sublead_before_RoccoR->SetDirectory(0);
        h_DY_eta_lead_before_RoccoR->SetDirectory(0);
        h_DY_eta_sublead_before_RoccoR->SetDirectory(0);
        h_DY_phi_lead_before_RoccoR->SetDirectory(0);
        h_DY_phi_sublead_before_RoccoR->SetDirectory(0);
        //
        s_pT_before_RoccoR->Add(h_DY_pT_before_RoccoR);
        s_rapi_before_RoccoR->Add(h_DY_rapi_before_RoccoR);
        s_pT_lead_before_RoccoR->Add(h_DY_pT_lead_before_RoccoR);
        s_pT_sublead_before_RoccoR->Add(h_DY_pT_sublead_before_RoccoR);
        s_eta_lead_before_RoccoR->Add(h_DY_eta_lead_before_RoccoR);
        s_eta_sublead_before_RoccoR->Add(h_DY_eta_sublead_before_RoccoR);
        s_phi_lead_before_RoccoR->Add(h_DY_phi_lead_before_RoccoR);
        s_phi_sublead_before_RoccoR->Add(h_DY_phi_sublead_before_RoccoR);

        h_DY_pT_before_EffCorr->SetFillColor(kOrange);
        h_DY_rapi_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_EffCorr->SetLineColor(kOrange);
        h_DY_rapi_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_EffCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_EffCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_EffCorr->SetDirectory(0);
        h_DY_rapi_before_EffCorr->SetDirectory(0);
        h_DY_pT_lead_before_EffCorr->SetDirectory(0);
        h_DY_pT_sublead_before_EffCorr->SetDirectory(0);
        h_DY_eta_lead_before_EffCorr->SetDirectory(0);
        h_DY_eta_sublead_before_EffCorr->SetDirectory(0);
        h_DY_phi_lead_before_EffCorr->SetDirectory(0);
        h_DY_phi_sublead_before_EffCorr->SetDirectory(0);
        //
        s_pT_before_EffCorr->Add(h_DY_pT_before_EffCorr);
        s_rapi_before_EffCorr->Add(h_DY_rapi_before_EffCorr);
        s_pT_lead_before_EffCorr->Add(h_DY_pT_lead_before_EffCorr);
        s_pT_sublead_before_EffCorr->Add(h_DY_pT_sublead_before_EffCorr);
        s_eta_lead_before_EffCorr->Add(h_DY_eta_lead_before_EffCorr);
        s_eta_sublead_before_EffCorr->Add(h_DY_eta_sublead_before_EffCorr);
        s_phi_lead_before_EffCorr->Add(h_DY_phi_lead_before_EffCorr);
        s_phi_sublead_before_EffCorr->Add(h_DY_phi_sublead_before_EffCorr);

        h_DY_pT_before_PVzCorr->SetFillColor(kOrange);
        h_DY_rapi_before_PVzCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_PVzCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_PVzCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_PVzCorr->SetLineColor(kOrange);
        h_DY_rapi_before_PVzCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_PVzCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_PVzCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_PVzCorr->SetDirectory(0);
        h_DY_rapi_before_PVzCorr->SetDirectory(0);
        h_DY_pT_lead_before_PVzCorr->SetDirectory(0);
        h_DY_pT_sublead_before_PVzCorr->SetDirectory(0);
        h_DY_eta_lead_before_PVzCorr->SetDirectory(0);
        h_DY_eta_sublead_before_PVzCorr->SetDirectory(0);
        h_DY_phi_lead_before_PVzCorr->SetDirectory(0);
        h_DY_phi_sublead_before_PVzCorr->SetDirectory(0);
        //
        s_pT_before_PVzCorr->Add(h_DY_pT_before_PVzCorr);
        s_rapi_before_PVzCorr->Add(h_DY_rapi_before_PVzCorr);
        s_pT_lead_before_PVzCorr->Add(h_DY_pT_lead_before_PVzCorr);
        s_pT_sublead_before_PVzCorr->Add(h_DY_pT_sublead_before_PVzCorr);
        s_eta_lead_before_PVzCorr->Add(h_DY_eta_lead_before_PVzCorr);
        s_eta_sublead_before_PVzCorr->Add(h_DY_eta_sublead_before_PVzCorr);
        s_phi_lead_before_PVzCorr->Add(h_DY_phi_lead_before_PVzCorr);
        s_phi_sublead_before_PVzCorr->Add(h_DY_phi_sublead_before_PVzCorr);

        h_DY_pT_before_L1Corr->SetFillColor(kOrange);
        h_DY_rapi_before_L1Corr->SetFillColor(kOrange);
        h_DY_pT_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_L1Corr->SetFillColor(kOrange);
        h_DY_eta_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_L1Corr->SetFillColor(kOrange);
        h_DY_phi_lead_before_L1Corr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_L1Corr->SetFillColor(kOrange);
        //
        h_DY_pT_before_L1Corr->SetLineColor(kOrange);
        h_DY_rapi_before_L1Corr->SetLineColor(kOrange);
        h_DY_pT_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_L1Corr->SetLineColor(kOrange);
        h_DY_eta_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_L1Corr->SetLineColor(kOrange);
        h_DY_phi_lead_before_L1Corr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_L1Corr->SetLineColor(kOrange);
        //
        h_DY_pT_before_L1Corr->SetDirectory(0);
        h_DY_rapi_before_L1Corr->SetDirectory(0);
        h_DY_pT_lead_before_L1Corr->SetDirectory(0);
        h_DY_pT_sublead_before_L1Corr->SetDirectory(0);
        h_DY_eta_lead_before_L1Corr->SetDirectory(0);
        h_DY_eta_sublead_before_L1Corr->SetDirectory(0);
        h_DY_phi_lead_before_L1Corr->SetDirectory(0);
        h_DY_phi_sublead_before_L1Corr->SetDirectory(0);
        //
        s_pT_before_L1Corr->Add(h_DY_pT_before_L1Corr);
        s_rapi_before_L1Corr->Add(h_DY_rapi_before_L1Corr);
        s_pT_lead_before_L1Corr->Add(h_DY_pT_lead_before_L1Corr);
        s_pT_sublead_before_L1Corr->Add(h_DY_pT_sublead_before_L1Corr);
        s_eta_lead_before_L1Corr->Add(h_DY_eta_lead_before_L1Corr);
        s_eta_sublead_before_L1Corr->Add(h_DY_eta_sublead_before_L1Corr);
        s_phi_lead_before_L1Corr->Add(h_DY_phi_lead_before_L1Corr);
        s_phi_sublead_before_L1Corr->Add(h_DY_phi_sublead_before_L1Corr);

        h_DY_pT_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_rapi_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_pT_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_pT_sublead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_eta_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_eta_sublead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_phi_lead_before_TopPtCorr->SetFillColor(kOrange);
        h_DY_phi_sublead_before_TopPtCorr->SetFillColor(kOrange);
        //
        h_DY_pT_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_rapi_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_pT_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_pT_sublead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_eta_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_eta_sublead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_phi_lead_before_TopPtCorr->SetLineColor(kOrange);
        h_DY_phi_sublead_before_TopPtCorr->SetLineColor(kOrange);
        //
        h_DY_pT_before_TopPtCorr->SetDirectory(0);
        h_DY_rapi_before_TopPtCorr->SetDirectory(0);
        h_DY_pT_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_pT_sublead_before_TopPtCorr->SetDirectory(0);
        h_DY_eta_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_eta_sublead_before_TopPtCorr->SetDirectory(0);
        h_DY_phi_lead_before_TopPtCorr->SetDirectory(0);
        h_DY_phi_sublead_before_TopPtCorr->SetDirectory(0);
        //
        s_pT_before_TopPtCorr->Add(h_DY_pT_before_TopPtCorr);
        s_rapi_before_TopPtCorr->Add(h_DY_rapi_before_TopPtCorr);
        s_pT_lead_before_TopPtCorr->Add(h_DY_pT_lead_before_TopPtCorr);
        s_pT_sublead_before_TopPtCorr->Add(h_DY_pT_sublead_before_TopPtCorr);
        s_eta_lead_before_TopPtCorr->Add(h_DY_eta_lead_before_TopPtCorr);
        s_eta_sublead_before_TopPtCorr->Add(h_DY_eta_sublead_before_TopPtCorr);
        s_phi_lead_before_TopPtCorr->Add(h_DY_phi_lead_before_TopPtCorr);
        s_phi_sublead_before_TopPtCorr->Add(h_DY_phi_sublead_before_TopPtCorr);

        h_DY_pT->SetFillColor(kOrange);
        h_DY_rapi->SetFillColor(kOrange);
        h_DY_pT_lead->SetFillColor(kOrange);
        h_DY_pT_sublead->SetFillColor(kOrange);
        h_DY_eta_lead->SetFillColor(kOrange);
        h_DY_eta_sublead->SetFillColor(kOrange);
        h_DY_phi_lead->SetFillColor(kOrange);
        h_DY_phi_sublead->SetFillColor(kOrange);
        //
        h_DY_pT->SetLineColor(kOrange);
        h_DY_rapi->SetLineColor(kOrange);
        h_DY_pT_lead->SetLineColor(kOrange);
        h_DY_pT_sublead->SetLineColor(kOrange);
        h_DY_eta_lead->SetLineColor(kOrange);
        h_DY_eta_sublead->SetLineColor(kOrange);
        h_DY_phi_lead->SetLineColor(kOrange);
        h_DY_phi_sublead->SetLineColor(kOrange);
        //
        h_DY_pT->SetDirectory(0);
        h_DY_rapi->SetDirectory(0);
        h_DY_pT_lead->SetDirectory(0);
        h_DY_pT_sublead->SetDirectory(0);
        h_DY_eta_lead->SetDirectory(0);
        h_DY_eta_sublead->SetDirectory(0);
        h_DY_phi_lead->SetDirectory(0);
        h_DY_phi_sublead->SetDirectory(0);
        //
        s_pT->Add(h_DY_pT);
        s_rapi->Add(h_DY_rapi);
        s_pT_lead->Add(h_DY_pT_lead);
        s_pT_sublead->Add(h_DY_pT_sublead);
        s_eta_lead->Add(h_DY_eta_lead);
        s_eta_sublead->Add(h_DY_eta_sublead);
        s_phi_lead->Add(h_DY_phi_lead);
        s_phi_sublead->Add(h_DY_phi_sublead);


//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT_before_RoccoR, *h_data_pT,
             *h_data_rapi_before_RoccoR, *h_data_rapi,
             *h_data_pT_lead_before_RoccoR, *h_data_pT_lead,
             *h_data_pT_sublead_before_RoccoR, *h_data_pT_sublead,
             *h_data_eta_lead_before_RoccoR, *h_data_eta_lead,
             *h_data_eta_sublead_before_RoccoR, *h_data_eta_sublead,
             *h_data_phi_lead_before_RoccoR, *h_data_phi_lead,
             *h_data_phi_sublead_before_RoccoR, *h_data_phi_sublead;

        f_data->GetObject("h_pT_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_before_RoccoR);
        f_data->GetObject("h_rapi_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_rapi_before_RoccoR);
        f_data->GetObject("h_pT_lead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead_before_RoccoR);
        f_data->GetObject("h_pT_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead_before_RoccoR);
        f_data->GetObject("h_eta_lead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead_before_RoccoR);
        f_data->GetObject("h_eta_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead_before_RoccoR);
        f_data->GetObject("h_phi_lead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead_before_RoccoR);
        f_data->GetObject("h_phi_sublead_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead_before_RoccoR);
        f_data->GetObject("h_pT_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT);
        f_data->GetObject("h_rapi_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_rapi);
        f_data->GetObject("h_pT_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_lead);
        f_data->GetObject("h_pT_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_pT_sublead);
        f_data->GetObject("h_eta_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_lead);
        f_data->GetObject("h_eta_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_eta_sublead);
        f_data->GetObject("h_phi_lead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_lead);
        f_data->GetObject("h_phi_sublead_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_phi_sublead);

        h_data_pT_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_rapi_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        //
        h_data_pT_before_RoccoR->SetMarkerColor(kBlack);
        h_data_rapi_before_RoccoR->SetMarkerColor(kBlack);
        h_data_pT_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_pT_sublead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_sublead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_lead_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_sublead_before_RoccoR->SetMarkerColor(kBlack);
        //
        h_data_pT_before_RoccoR->SetLineColor(kBlack);
        h_data_rapi_before_RoccoR->SetLineColor(kBlack);
        h_data_pT_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_pT_sublead_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_sublead_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_lead_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_sublead_before_RoccoR->SetLineColor(kBlack);
        //
        h_data_pT_before_RoccoR->SetDirectory(0);
        h_data_rapi_before_RoccoR->SetDirectory(0);
        h_data_pT_lead_before_RoccoR->SetDirectory(0);
        h_data_pT_sublead_before_RoccoR->SetDirectory(0);
        h_data_eta_lead_before_RoccoR->SetDirectory(0);
        h_data_eta_sublead_before_RoccoR->SetDirectory(0);
        h_data_phi_lead_before_RoccoR->SetDirectory(0);
        h_data_phi_sublead_before_RoccoR->SetDirectory(0);

        h_data_pT->SetMarkerStyle(kFullDotLarge);
        h_data_rapi->SetMarkerStyle(kFullDotLarge);
        h_data_pT_lead->SetMarkerStyle(kFullDotLarge);
        h_data_pT_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_lead->SetMarkerStyle(kFullDotLarge);
        h_data_eta_sublead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_lead->SetMarkerStyle(kFullDotLarge);
        h_data_phi_sublead->SetMarkerStyle(kFullDotLarge);
        //
        h_data_pT->SetMarkerColor(kBlack);
        h_data_rapi->SetMarkerColor(kBlack);
        h_data_pT_lead->SetMarkerColor(kBlack);
        h_data_pT_sublead->SetMarkerColor(kBlack);
        h_data_eta_lead->SetMarkerColor(kBlack);
        h_data_eta_sublead->SetMarkerColor(kBlack);
        h_data_phi_lead->SetMarkerColor(kBlack);
        h_data_phi_sublead->SetMarkerColor(kBlack);
        //
        h_data_pT->SetLineColor(kBlack);
        h_data_rapi->SetLineColor(kBlack);
        h_data_pT_lead->SetLineColor(kBlack);
        h_data_pT_sublead->SetLineColor(kBlack);
        h_data_eta_lead->SetLineColor(kBlack);
        h_data_eta_sublead->SetLineColor(kBlack);
        h_data_phi_lead->SetLineColor(kBlack);
        h_data_phi_sublead->SetLineColor(kBlack);
        //
        h_data_pT->SetDirectory(0);
        h_data_rapi->SetDirectory(0);
        h_data_pT_lead->SetDirectory(0);
        h_data_pT_sublead->SetDirectory(0);
        h_data_eta_lead->SetDirectory(0);
        h_data_eta_sublead->SetDirectory(0);
        h_data_phi_lead->SetDirectory(0);
        h_data_phi_sublead->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_pT_before_PUCorr, *RP_pT_before_RoccoR, *RP_pT_before_EffCorr,
                      *RP_pT_before_PVzCorr, *RP_pT_before_L1Corr, *RP_pT_before_TopPtCorr, *RP_pT,
                      *RP_rapi_before_PUCorr, *RP_rapi_before_RoccoR, *RP_rapi_before_EffCorr,
                      *RP_rapi_before_PVzCorr, *RP_rapi_before_L1Corr, *RP_rapi_before_TopPtCorr, *RP_rapi,
                      *RP_pT_lead_before_PUCorr, *RP_pT_lead_before_RoccoR, *RP_pT_lead_before_EffCorr,
                      *RP_pT_lead_before_PVzCorr, *RP_pT_lead_before_L1Corr, *RP_pT_lead_before_TopPtCorr, *RP_pT_lead,
                      *RP_pT_sublead_before_PUCorr, *RP_pT_sublead_before_RoccoR, *RP_pT_sublead_before_EffCorr,
                      *RP_pT_sublead_before_PVzCorr, *RP_pT_sublead_before_L1Corr, *RP_pT_sublead_before_TopPtCorr, *RP_pT_sublead,
                      *RP_eta_lead_before_PUCorr, *RP_eta_lead_before_RoccoR, *RP_eta_lead_before_EffCorr,
                      *RP_eta_lead_before_PVzCorr, *RP_eta_lead_before_L1Corr, *RP_eta_lead_before_TopPtCorr, *RP_eta_lead,
                      *RP_eta_sublead_before_PUCorr, *RP_eta_sublead_before_RoccoR, *RP_eta_sublead_before_EffCorr,
                      *RP_eta_sublead_before_PVzCorr, *RP_eta_sublead_before_L1Corr, *RP_eta_sublead_before_TopPtCorr, *RP_eta_sublead,
                      *RP_phi_lead_before_PUCorr, *RP_phi_lead_before_RoccoR, *RP_phi_lead_before_EffCorr,
                      *RP_phi_lead_before_PVzCorr, *RP_phi_lead_before_L1Corr, *RP_phi_lead_before_TopPtCorr, *RP_phi_lead,
                      *RP_phi_sublead_before_PUCorr, *RP_phi_sublead_before_RoccoR, *RP_phi_sublead_before_EffCorr,
                      *RP_phi_sublead_before_PVzCorr, *RP_phi_sublead_before_L1Corr, *RP_phi_sublead_before_TopPtCorr, *RP_phi_sublead;

        RP_pT_before_PUCorr = new myRatioPlot_t("RP_pT_before_PUCorr", s_pT_before_PUCorr, h_data_pT_before_RoccoR);
        RP_rapi_before_PUCorr = new myRatioPlot_t("RP_rapi_before_PUCorr", s_rapi_before_PUCorr, h_data_rapi_before_RoccoR);
        RP_pT_lead_before_PUCorr = new myRatioPlot_t("RP_pT_lead_before_PUCorr", s_pT_lead_before_PUCorr, h_data_pT_lead_before_RoccoR);
        RP_pT_sublead_before_PUCorr = new myRatioPlot_t("RP_pT_sublead_before_PUCorr", s_pT_sublead_before_PUCorr, h_data_pT_sublead_before_RoccoR);
        RP_eta_lead_before_PUCorr = new myRatioPlot_t("RP_eta_lead_before_PUCorr", s_eta_lead_before_PUCorr, h_data_eta_lead_before_RoccoR);
        RP_eta_sublead_before_PUCorr = new myRatioPlot_t("RP_eta_sublead_before_PUCorr", s_eta_sublead_before_PUCorr, h_data_eta_sublead_before_RoccoR);
        RP_phi_lead_before_PUCorr = new myRatioPlot_t("RP_phi_lead_before_PUCorr", s_phi_lead_before_PUCorr, h_data_phi_lead_before_RoccoR);
        RP_phi_sublead_before_PUCorr = new myRatioPlot_t("RP_phi_sublead_before_PUCorr", s_phi_sublead_before_PUCorr, h_data_phi_sublead_before_RoccoR);

        RP_pT_before_RoccoR = new myRatioPlot_t("RP_pT_before_RoccoR", s_pT_before_RoccoR, h_data_pT_before_RoccoR);
        RP_rapi_before_RoccoR = new myRatioPlot_t("RP_rapi_before_RoccoR", s_rapi_before_RoccoR, h_data_rapi_before_RoccoR);
        RP_pT_lead_before_RoccoR = new myRatioPlot_t("RP_pT_lead_before_RoccoR", s_pT_lead_before_RoccoR, h_data_pT_lead_before_RoccoR);
        RP_pT_sublead_before_RoccoR = new myRatioPlot_t("RP_pT_sublead_before_RoccoR", s_pT_sublead_before_RoccoR, h_data_pT_sublead_before_RoccoR);
        RP_eta_lead_before_RoccoR = new myRatioPlot_t("RP_eta_lead_before_RoccoR", s_eta_lead_before_RoccoR, h_data_eta_lead_before_RoccoR);
        RP_eta_sublead_before_RoccoR = new myRatioPlot_t("RP_eta_sublead_before_RoccoR", s_eta_sublead_before_RoccoR, h_data_eta_sublead_before_RoccoR);
        RP_phi_lead_before_RoccoR = new myRatioPlot_t("RP_phi_lead_before_RoccoR", s_phi_lead_before_RoccoR, h_data_phi_lead_before_RoccoR);
        RP_phi_sublead_before_RoccoR = new myRatioPlot_t("RP_phi_sublead_before_RoccoR", s_phi_sublead_before_RoccoR, h_data_phi_sublead_before_RoccoR);

        RP_pT_before_EffCorr = new myRatioPlot_t("RP_pT_before_EffCorr", s_pT_before_EffCorr, h_data_pT);
        RP_rapi_before_EffCorr = new myRatioPlot_t("RP_rapi_before_EffCorr", s_rapi_before_EffCorr, h_data_rapi);
        RP_pT_lead_before_EffCorr = new myRatioPlot_t("RP_pT_lead_before_EffCorr", s_pT_lead_before_EffCorr, h_data_pT_lead);
        RP_pT_sublead_before_EffCorr = new myRatioPlot_t("RP_pT_sublead_before_EffCorr", s_pT_sublead_before_EffCorr, h_data_pT_sublead);
        RP_eta_lead_before_EffCorr = new myRatioPlot_t("RP_eta_lead_before_EffCorr", s_eta_lead_before_EffCorr, h_data_eta_lead);
        RP_eta_sublead_before_EffCorr = new myRatioPlot_t("RP_eta_sublead_before_EffCorr", s_eta_sublead_before_EffCorr, h_data_eta_sublead);
        RP_phi_lead_before_EffCorr = new myRatioPlot_t("RP_phi_lead_before_EffCorr", s_phi_lead_before_EffCorr, h_data_phi_lead);
        RP_phi_sublead_before_EffCorr = new myRatioPlot_t("RP_phi_sublead_before_EffCorr", s_phi_sublead_before_EffCorr, h_data_phi_sublead);

        RP_pT_before_PVzCorr = new myRatioPlot_t("RP_pT_before_PVzCorr", s_pT_before_PVzCorr, h_data_pT);
        RP_rapi_before_PVzCorr = new myRatioPlot_t("RP_rapi_before_PVzCorr", s_rapi_before_PVzCorr, h_data_rapi);
        RP_pT_lead_before_PVzCorr = new myRatioPlot_t("RP_pT_lead_before_PVzCorr", s_pT_lead_before_PVzCorr, h_data_pT_lead);
        RP_pT_sublead_before_PVzCorr = new myRatioPlot_t("RP_pT_sublead_before_PVzCorr", s_pT_sublead_before_PVzCorr, h_data_pT_sublead);
        RP_eta_lead_before_PVzCorr = new myRatioPlot_t("RP_eta_lead_before_PVzCorr", s_eta_lead_before_PVzCorr, h_data_eta_lead);
        RP_eta_sublead_before_PVzCorr = new myRatioPlot_t("RP_eta_sublead_before_PVzCorr", s_eta_sublead_before_PVzCorr, h_data_eta_sublead);
        RP_phi_lead_before_PVzCorr = new myRatioPlot_t("RP_phi_lead_before_PVzCorr", s_phi_lead_before_PVzCorr, h_data_phi_lead);
        RP_phi_sublead_before_PVzCorr = new myRatioPlot_t("RP_phi_sublead_before_PVzCorr", s_phi_sublead_before_PVzCorr, h_data_phi_sublead);

        RP_pT_before_L1Corr = new myRatioPlot_t("RP_pT_before_L1Corr", s_pT_before_L1Corr, h_data_pT);
        RP_rapi_before_L1Corr = new myRatioPlot_t("RP_rapi_before_L1Corr", s_rapi_before_L1Corr, h_data_rapi);
        RP_pT_lead_before_L1Corr = new myRatioPlot_t("RP_pT_lead_before_L1Corr", s_pT_lead_before_L1Corr, h_data_pT_lead);
        RP_pT_sublead_before_L1Corr = new myRatioPlot_t("RP_pT_sublead_before_L1Corr", s_pT_sublead_before_L1Corr, h_data_pT_sublead);
        RP_eta_lead_before_L1Corr = new myRatioPlot_t("RP_eta_lead_before_L1Corr", s_eta_lead_before_L1Corr, h_data_eta_lead);
        RP_eta_sublead_before_L1Corr = new myRatioPlot_t("RP_eta_sublead_before_L1Corr", s_eta_sublead_before_L1Corr, h_data_eta_sublead);
        RP_phi_lead_before_L1Corr = new myRatioPlot_t("RP_phi_lead_before_L1Corr", s_phi_lead_before_L1Corr, h_data_phi_lead);
        RP_phi_sublead_before_L1Corr = new myRatioPlot_t("RP_phi_sublead_before_L1Corr", s_phi_sublead_before_L1Corr, h_data_phi_sublead);

        RP_pT_before_TopPtCorr = new myRatioPlot_t("RP_pT_before_TopPtCorr", s_pT_before_TopPtCorr, h_data_pT);
        RP_rapi_before_TopPtCorr = new myRatioPlot_t("RP_rapi_before_TopPtCorr", s_rapi_before_TopPtCorr, h_data_rapi);
        RP_pT_lead_before_TopPtCorr = new myRatioPlot_t("RP_pT_lead_before_TopPtCorr", s_pT_lead_before_TopPtCorr, h_data_pT_lead);
        RP_pT_sublead_before_TopPtCorr = new myRatioPlot_t("RP_pT_sublead_before_TopPtCorr", s_pT_sublead_before_TopPtCorr, h_data_pT_sublead);
        RP_eta_lead_before_TopPtCorr = new myRatioPlot_t("RP_eta_lead_before_TopPtCorr", s_eta_lead_before_TopPtCorr, h_data_eta_lead);
        RP_eta_sublead_before_TopPtCorr = new myRatioPlot_t("RP_eta_sublead_before_TopPtCorr", s_eta_sublead_before_TopPtCorr, h_data_eta_sublead);
        RP_phi_lead_before_TopPtCorr = new myRatioPlot_t("RP_phi_lead_before_TopPtCorr", s_phi_lead_before_TopPtCorr, h_data_phi_lead);
        RP_phi_sublead_before_TopPtCorr = new myRatioPlot_t("RP_phi_sublead_before_TopPtCorr", s_phi_sublead_before_TopPtCorr, h_data_phi_sublead);

        RP_pT = new myRatioPlot_t("RP_pT", s_pT, h_data_pT);
        RP_rapi = new myRatioPlot_t("RP_rapi", s_rapi, h_data_rapi);
        RP_pT_lead = new myRatioPlot_t("RP_pT_lead", s_pT_lead, h_data_pT_lead);
        RP_pT_sublead = new myRatioPlot_t("RP_pT_sublead", s_pT_sublead, h_data_pT_sublead);
        RP_eta_lead = new myRatioPlot_t("RP_eta_lead", s_eta_lead, h_data_eta_lead);
        RP_eta_sublead = new myRatioPlot_t("RP_eta_sublead", s_eta_sublead, h_data_eta_sublead);
        RP_phi_lead = new myRatioPlot_t("RP_phi_lead", s_phi_lead, h_data_phi_lead);
        RP_phi_sublead = new myRatioPlot_t("RP_phi_sublead", s_phi_sublead, h_data_phi_sublead);

        RP_pT_before_PUCorr->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before PU correction", 0, 1000);
        RP_rapi_before_PUCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before PU correction", -3, 3);
//        RP_rapi_before_PUCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} prie#check{s} visas pataisas", -3, 3, "Eksp./MC");
        RP_pT_lead_before_PUCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before PU correction", 0, 1000);
        RP_pT_sublead_before_PUCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before PU correction", 0, 1000);
        RP_eta_lead_before_PUCorr->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before PU correction", -3.5, 3.5);
//        RP_eta_lead_before_PUCorr->SetPlots("#eta (#mu_{#lower[-0.4]{1}}) prie#check{s} visas pataisas", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_PUCorr->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before PU correction", -3.5, 3.5);
//        RP_eta_sublead_before_PUCorr->SetPlots("#eta (#mu_{#lower[-0.4]{2}}) prie#check{s} visas pataisas", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_PUCorr->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before PU correction", -4, 4);
        RP_phi_sublead_before_PUCorr->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before PU correction", -4, 4);

        RP_pT_before_RoccoR->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before Rochester correction", 0, 1000);
        RP_rapi_before_RoccoR->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before Rochester correction", -3, 3);
        RP_pT_lead_before_RoccoR->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before Rochester correction", 0, 1000);
        RP_pT_sublead_before_RoccoR->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before Rochester correction", 0, 1000);
        RP_eta_lead_before_RoccoR->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before Rochester correction", -3.5, 3.5);
        RP_eta_sublead_before_RoccoR->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before Rochester correction", -3.5, 3.5);
        RP_phi_lead_before_RoccoR->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before Rochester correction", -4, 4);
        RP_phi_sublead_before_RoccoR->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before Rochester correction", -4, 4);

        RP_pT_before_EffCorr->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before Efficiency SF", 0, 1000);
        RP_rapi_before_EffCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before Efficiency SF", -3, 3);
        RP_pT_lead_before_EffCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before Efficiency SF", 0, 1000);
        RP_pT_sublead_before_EffCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before Efficiency SF", 0, 1000);
        RP_eta_lead_before_EffCorr->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before Efficiency SF", -3.5, 3.5);
        RP_eta_sublead_before_EffCorr->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before Efficiency SF", -3.5, 3.5);
        RP_phi_lead_before_EffCorr->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before Efficiency SF", -4, 4);
        RP_phi_sublead_before_EffCorr->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before Efficiency SF", -4, 4);

        RP_pT_before_PVzCorr->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before PVz correction", 0, 1000);
        RP_rapi_before_PVzCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before PVz correction", -3, 3);
//        RP_rapi_before_PVzCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3, 3, "Eksp./MC");
        RP_pT_lead_before_PVzCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before PVz correction", 0, 1000);
        RP_pT_sublead_before_PVzCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before PVz correction", 0, 1000);
        RP_eta_lead_before_PVzCorr->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before PVz correction", -3.5, 3.5);
//        RP_eta_lead_before_PVzCorr->SetPlots("#eta (#mu_{#lower[-0.4]{1}}) prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_PVzCorr->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before PVz correction", -3.5, 3.5);
//        RP_eta_sublead_before_PVzCorr->SetPlots("#eta (#mu_{#lower[-0.4]{2}}) prie#check{s} PVz pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_PVzCorr->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before PVz correction", -4, 4);
        RP_phi_sublead_before_PVzCorr->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before PVz correction", -4, 4);

        RP_pT_before_L1Corr->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before L1 prefiring correction", 0, 1000);
        RP_rapi_before_L1Corr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before L1 prefiring correction", -3, 3);
//        RP_rapi_before_L1Corr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3, 3, "Eksp./MC");
        RP_pT_lead_before_L1Corr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before  L1 prefiring correction", 0, 1000);
        RP_pT_sublead_before_L1Corr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before  L1 prefiring correction", 0, 1000);
        RP_eta_lead_before_L1Corr->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before L1 prefiring correction", -3.5, 3.5);
//        RP_eta_lead_before_L1Corr->SetPlots("#eta (#mu_{#lower[-0.4]{1}}) prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_eta_sublead_before_L1Corr->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before  L1 prefiring correction", -3.5, 3.5);
//        RP_eta_sublead_before_L1Corr->SetPlots("#eta (#mu_{#lower[-0.4]{2}}) prie#check{s} trigerio pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead_before_L1Corr->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before  L1 prefiring correction", -4, 4);
        RP_phi_sublead_before_L1Corr->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before L1 prefiring correction", -4, 4);

        RP_pT_before_TopPtCorr->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c] before top p_{#lower[-0.25]{T}} reweighting", 0, 1000);
        RP_rapi_before_TopPtCorr->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} before top p_{#lower[-0.25]{T}} reweighting", -3, 3);
        RP_pT_lead_before_TopPtCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c] before top p_{#lower[-0.25]{T}} reweighting", 0, 1000);
        RP_pT_sublead_before_TopPtCorr->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c] before top p_{#lower[-0.25]{T}} reweighting", 0, 1000);
        RP_eta_lead_before_TopPtCorr->SetPlots("#eta (#mu_{#lower[-0.4]{lead}}) before top p_{#lower[-0.25]{T}} reweighting", -3.5, 3.5);
        RP_eta_sublead_before_TopPtCorr->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}}) before top p_{#lower[-0.25]{T}} reweighting", -3.5, 3.5);
        RP_phi_lead_before_TopPtCorr->SetPlots("#phi (#mu_{#lower[-0.4]{lead}}) before top p_{#lower[-0.25]{T}} reweighting", -4, 4);
        RP_phi_sublead_before_TopPtCorr->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}}) before top p_{#lower[-0.25]{T}} reweighting", -4, 4);

        RP_pT->SetPlots("p_{#lower[-0.25]{T #scale[1.2]{#mu#mu}}} [GeV/c]", 0, 1000);
        RP_rapi->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}}", -3, 3);
//        RP_rapi->SetPlots("y_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} pritaikius pataisas", -3, 3, "Eksp./MC");
        RP_pT_lead->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{lead}}) [GeV/c]", 0, 1000);
        RP_pT_sublead->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{sublead}}) [GeV/c]", 0, 1000);
        RP_eta_lead->SetPlots("#eta (#mu_{#lower[-0.4]{lead}})", -3.5, 3.5);
//        RP_eta_lead->SetPlots("#eta (#mu_{#lower[-0.4]{1}}) pritaikius pataisas", -3.5, 3.5, "Eksp./MC");
//        RP_eta_sublead->SetPlots("#eta (#mu_{#lower[-0.4]{sublead}})", -3.5, 3.5);
        RP_eta_sublead->SetPlots("#eta (#mu_{#lower[-0.4]{2}}) pritaikius pataisas", -3.5, 3.5, "Eksp./MC");
        RP_phi_lead->SetPlots("#phi (#mu_{#lower[-0.4]{lead}})", -4, 4);
        RP_phi_sublead->SetPlots("#phi (#mu_{#lower[-0.4]{sublead}})", -4, 4);

        TLegend *legend = new TLegend(0.55, 0.72, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_pT, "Data", "lp");
//        legend->AddEntry(h_data_pT, "Matavimas", "lp");
        legend->AddEntry(h_DY_pT, "DY#rightarrow #mu#mu", "f");
        legend->AddEntry(h_bkg_pT[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_pT[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_pT[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_pT[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_pT[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_pT[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_pT[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_pT[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_pT[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_pT_before_PUCorr->ImportLegend(legend);
        RP_rapi_before_PUCorr->ImportLegend(legend);
        RP_pT_lead_before_PUCorr->ImportLegend(legend);
        RP_pT_sublead_before_PUCorr->ImportLegend(legend);
        RP_eta_lead_before_PUCorr->ImportLegend(legend);
        RP_eta_sublead_before_PUCorr->ImportLegend(legend);
        RP_phi_lead_before_PUCorr->ImportLegend(legend);
        RP_phi_sublead_before_PUCorr->ImportLegend(legend);
        RP_pT_before_RoccoR->ImportLegend(legend);
        RP_rapi_before_RoccoR->ImportLegend(legend);
        RP_pT_lead_before_RoccoR->ImportLegend(legend);;
        RP_pT_sublead_before_RoccoR->ImportLegend(legend);
        RP_eta_lead_before_RoccoR->ImportLegend(legend);
        RP_eta_sublead_before_RoccoR->ImportLegend(legend);
        RP_phi_lead_before_RoccoR->ImportLegend(legend);
        RP_phi_sublead_before_RoccoR->ImportLegend(legend);
        RP_pT_before_EffCorr->ImportLegend(legend);
        RP_rapi_before_EffCorr->ImportLegend(legend);
        RP_pT_lead_before_EffCorr->ImportLegend(legend);
        RP_pT_sublead_before_EffCorr->ImportLegend(legend);
        RP_eta_lead_before_EffCorr->ImportLegend(legend);
        RP_eta_sublead_before_EffCorr->ImportLegend(legend);
        RP_phi_lead_before_EffCorr->ImportLegend(legend);
        RP_phi_sublead_before_EffCorr->ImportLegend(legend);
        RP_pT_before_PVzCorr->ImportLegend(legend);
        RP_rapi_before_PVzCorr->ImportLegend(legend);
        RP_pT_lead_before_PVzCorr->ImportLegend(legend);
        RP_pT_sublead_before_PVzCorr->ImportLegend(legend);
        RP_eta_lead_before_PVzCorr->ImportLegend(legend);
        RP_eta_sublead_before_PVzCorr->ImportLegend(legend);
        RP_phi_lead_before_PVzCorr->ImportLegend(legend);
        RP_phi_sublead_before_PVzCorr->ImportLegend(legend);
        RP_pT_before_L1Corr->ImportLegend(legend);
        RP_rapi_before_L1Corr->ImportLegend(legend);
        RP_pT_lead_before_L1Corr->ImportLegend(legend);
        RP_pT_sublead_before_L1Corr->ImportLegend(legend);
        RP_eta_lead_before_L1Corr->ImportLegend(legend);
        RP_eta_sublead_before_L1Corr->ImportLegend(legend);
        RP_phi_lead_before_L1Corr->ImportLegend(legend);
        RP_phi_sublead_before_L1Corr->ImportLegend(legend);
        RP_pT_before_TopPtCorr->ImportLegend(legend);
        RP_rapi_before_TopPtCorr->ImportLegend(legend);
        RP_pT_lead_before_TopPtCorr->ImportLegend(legend);
        RP_pT_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_eta_lead_before_TopPtCorr->ImportLegend(legend);
        RP_eta_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_phi_lead_before_TopPtCorr->ImportLegend(legend);
        RP_phi_sublead_before_TopPtCorr->ImportLegend(legend);
        RP_pT->ImportLegend(legend);
        RP_rapi->ImportLegend(legend);
        RP_pT_lead->ImportLegend(legend);
        RP_pT_sublead->ImportLegend(legend);
        RP_eta_lead->ImportLegend(legend);
        RP_eta_sublead->ImportLegend(legend);
        RP_phi_lead->ImportLegend(legend);
        RP_phi_sublead->ImportLegend(legend);

        RP_pT_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_PUCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_rapi_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_RoccoR->Draw(0.8, 6e7, 0);
        RP_pT_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_EffCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_PVzCorr->Draw(0.8, 6e7, 0);
        RP_pT_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_rapi_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_L1Corr->Draw(0.8, 6e7, 0);
        RP_pT_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_rapi_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_pT_lead_before_TopPtCorr->Draw(0.8, 1e7, 0);
        RP_pT_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_eta_lead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_eta_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_phi_lead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_phi_sublead_before_TopPtCorr->Draw(0.8, 6e7, 0);
        RP_pT->Draw(0.8, 6e7, 0);
        RP_rapi->Draw(0.8, 6e7, 0);
        RP_pT_lead->Draw(0.8, 6e7, 0);
        RP_pT_sublead->Draw(0.8, 6e7, 0);
        RP_eta_lead->Draw(0.8, 6e7, 0);
        RP_eta_sublead->Draw(0.8, 6e7, 0);
        RP_phi_lead->Draw(0.8, 6e7, 0);
        RP_phi_sublead->Draw(0.8, 6e7, 0);

    } // End of if(Pt, rapi, pT, eta, phi)

//################################# nVTX #################################################

    if(whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP"))
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[9], *h_bkg_nVTX_before_EffCorr[9], *h_bkg_nVTX[9];
        Int_t iter = 0;

        for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _MuMu_QCDMuEnriched_Full)
            {
                iter++;
                continue;
            }
            f_bkg->GetObject("h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter]);
            f_bkg->GetObject("h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter]);
            f_bkg->GetObject("h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter]);
            removeNegativeBins(h_bkg_nVTX_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_nVTX_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_nVTX[iter]);

            Color_t color = kBlack;
            if (pr == _MuMu_QCDMuEnriched_Full) color = kRed + 3;
            if (pr == _MuMu_WJets_Full) color = kRed - 2;
            if (pr == _MuMu_WW) color = kMagenta - 5;
            if (pr == _MuMu_WZ) color = kMagenta - 2;
            if (pr == _MuMu_ZZ) color = kMagenta - 6;
            if (pr == _MuMu_tbarW) color = kGreen - 2;
            if (pr == _MuMu_tW) color = kGreen + 2;
            if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
            if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            s_nVTX_before_PUCorr->Add(h_bkg_nVTX_before_PUCorr[iter]);

            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            s_nVTX_before_EffCorr->Add(h_bkg_nVTX_before_EffCorr[iter]);

            h_bkg_nVTX[iter]->SetFillColor(color);
            h_bkg_nVTX[iter]->SetLineColor(color);
            h_bkg_nVTX[iter]->SetDirectory(0);
            s_nVTX->Add(h_bkg_nVTX[iter]);

            iter++;

            if (pr == _MuMu_WJets_Full)
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if (pr == _MuMu_tW)
                pr = _MuMu_VVnST; // next - ttbar
            if (pr == _MuMu_DYTauTau_Full) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_nVTX_before_PUCorr, *h_DY_nVTX_before_EffCorr, *h_DY_nVTX;

        f_DY->GetObject("h_nVTX_before_PUCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX_before_PUCorr);
        f_DY->GetObject("h_nVTX_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX_before_EffCorr);
        f_DY->GetObject("h_nVTX_"+Mgr.Procname[_MuMu_DY_Full], h_DY_nVTX);
        removeNegativeBins(h_DY_nVTX_before_PUCorr);
        removeNegativeBins(h_DY_nVTX_before_EffCorr);
        removeNegativeBins(h_DY_nVTX);

        h_DY_nVTX_before_PUCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_PUCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_PUCorr->SetDirectory(0);
        s_nVTX_before_PUCorr->Add(h_DY_nVTX_before_PUCorr);

        h_DY_nVTX_before_EffCorr->SetFillColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetLineColor(kOrange);
        h_DY_nVTX_before_EffCorr->SetDirectory(0);
        s_nVTX_before_EffCorr->Add(h_DY_nVTX_before_EffCorr);

        h_DY_nVTX->SetFillColor(kOrange);
        h_DY_nVTX->SetLineColor(kOrange);
        h_DY_nVTX->SetDirectory(0);
        s_nVTX->Add(h_DY_nVTX);

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;

        f_data->GetObject("h_nVTX_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_nVTX);

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t("RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX);
        RP_nVTX_before_EffCorr = new myRatioPlot_t("RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX);
        RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_data_nVTX);

        RP_nVTX_before_PUCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.15]{#mu#mu}}} before PU correction", 0, 50);
//        RP_nVTX_before_PUCorr->SetPlots("N_{#lower[-0.2]{PV} #lower[-0.25]{#scale[1.15]{#mu#mu}}} prie#check{s} pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", 0, 50, "Eksp./MC");
        RP_nVTX_before_EffCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.15]{#mu#mu}}} before Efficiency SF", 0, 50);
        RP_nVTX->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.15]{#mu#mu}}}", 0, 50);
//        RP_nVTX->SetPlots("N_{#lower[-0.2]{PV} #lower[-0.25]{#scale[1.15]{#mu#mu}}} pritaikius pataisa_{#kern[-0.95]{#lower[-0.4]{#scale[0.8]{c}}}}", 0, 50, "Eksp./MC");

        TLegend *legend = new TLegend(0.6, 0.75, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_nVTX, "Data", "lp");
//        legend->AddEntry(h_data_nVTX, "Matavimas", "lp");
        legend->AddEntry(h_DY_nVTX_before_PUCorr, "DY#rightarrow#mu#mu", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_nVTX_before_PUCorr[0], "#font[12]{#scale[1.1]{QCD}}", "f");

        RP_nVTX_before_PUCorr->ImportLegend(legend);
        RP_nVTX_before_EffCorr->ImportLegend(legend);
        RP_nVTX->ImportLegend(legend);

        RP_nVTX_before_PUCorr->Draw(0.5, 3e7, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e7, 0);
        RP_nVTX->Draw(0.5, 3e7, 0);

        cout << "nVTX Chi^2 before PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX_before_PUCorr) << endl;
        cout << "nVTX Chi^2 after PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX) << endl;

        Double_t EvtPercentage = 0;
        for (Int_t i=11; i<32; i++)
        {
            EvtPercentage += h_data_nVTX->GetBinContent(i);
        }
        EvtPercentage /= h_data_nVTX->Integral(1, h_data_nVTX->GetSize()-2) * 0.01;
        cout << "There are " << EvtPercentage << "% of events in 10-30 nVTX range" << endl;

    } // End of if(nVTX)
*/

    f_DY->Close();
    f_bkg->Close();
    f_data->Close();

} // End of MuMu_HistDrawer()


/// ################################################################################ ///
/// -------------------------------- EMu events ------------------------------------ ///
/// ################################################################################ ///
void EMu_HistDrawer (TString whichGraphs , TString type)
{
    if (!whichGraphs.Length())
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t isWJ = 0;
    isWJ = 1; // UNCOMMENT THIS IF YOU WANT TO INCLUDE W+JETS INTO HISTOGRAMS
    Int_t UseFR = 0;
    UseFR = 1; // UNCOMMENT THIS IF YOU WANT TO WANT TO INCLUDE QCD AND WJETS ESTIMATIONS FROM FR
    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc(_EMu_Bkg_Full);
    cout << "Hists location: " << Mgr.HistLocation << endl;
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile(name_bkg, "READ");
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_EMu_SingleMuon_Full);
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile(name_data, "READ");
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

//################################# INVARIANT MASS #################################################

    if(whichGraphs=="ALL" || whichGraphs=="INVM" || whichGraphs=="INVMASS")
    {
        count_drawn++;

        THStack *s_mass_before_PUCorr = new THStack("s_mass_before_PUCorr", "");
        THStack *s_mass_before_RoccoR = new THStack("s_mass_before_RoccoR", "");
        THStack *s_mass_before_EffCorr = new THStack("s_mass_before_EffCorr", "");
        THStack *s_mass_before_PVzCorr = new THStack("s_mass_before_PVzCorr", "");
        THStack *s_mass_before_L1Corr = new THStack("s_mass_before_L1Corr", "");
        THStack *s_mass_before_TopPtCorr = new THStack("s_mass_before_TopPtCorr", "");
        THStack *s_mass_fine = new THStack("s_mass_fine", "");
        THStack *s_mass = new THStack("s_mass", "");
        THStack *s_mass2 = new THStack("s_mass2", "");
        THStack *s_SS_mass_before_PUCorr = new THStack("s_SS_mass_before_PUCorr", "");
        THStack *s_SS_mass_before_RoccoR = new THStack("s_SS_mass_before_RoccoR", "");
        THStack *s_SS_mass_before_EffCorr = new THStack("s_SS_mass_before_EffCorr", "");
        THStack *s_SS_mass_before_PVzCorr = new THStack("s_SS_mass_before_PVzCorr", "");
        THStack *s_SS_mass_before_L1Corr = new THStack("s_SS_mass_before_L1Corr", "");
        THStack *s_SS_mass_before_TopPtCorr = new THStack("s_SS_mass_before_TopPtCorr", "");
        THStack *s_SS_mass_fine = new THStack("s_SS_mass_fine", "");
        THStack *s_SS_mass = new THStack("s_SS_mass", "");
        THStack *s_SS_mass2 = new THStack("s_SS_mass2", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_before_PUCorr[8], *h_bkg_mass_before_RoccoR[8], *h_bkg_mass_before_EffCorr[8],
             *h_bkg_mass_before_PVzCorr[8], *h_bkg_mass_before_L1Corr[8], *h_bkg_mass_before_TopPtCorr[8],
             *h_bkg_mass_fine[8], *h_bkg_mass[8], *h_bkg_mass2[8],
             *h_SS_bkg_mass_before_PUCorr[8], *h_SS_bkg_mass_before_RoccoR[8], *h_SS_bkg_mass_before_EffCorr[8],
             *h_SS_bkg_mass_before_PVzCorr[8], *h_SS_bkg_mass_before_L1Corr[8], *h_SS_bkg_mass_before_TopPtCorr[8],
             *h_SS_bkg_mass_fine[8], *h_SS_bkg_mass[8], *h_SS_bkg_mass2[8];
        Int_t iter = 0;

        Double_t ratio_WJets_SSvsOS;
        Double_t avg_ratio_WJets_SSvsOS = 0;

        for (SelProc_t pr = _EMu_WJets_Full; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (!isWJ && pr == _EMu_WJets_Full) continue;
            if (pr == _EndOf_EMu_VVnST_Normal) continue;

            f_bkg->GetObject("h_emu_mass_before_PUCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PUCorr[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_PUCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_PUCorr[iter]);
            f_bkg->GetObject("h_emu_mass_before_RocCorr_"+Mgr.Procname[pr], h_bkg_mass_before_RoccoR[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_RocCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_RoccoR[iter]);
            f_bkg->GetObject("h_emu_mass_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_before_EffCorr[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_EffCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_EffCorr[iter]);
            f_bkg->GetObject("h_emu_mass_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_mass_before_PVzCorr[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_PVzCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_PVzCorr[iter]);
            f_bkg->GetObject("h_emu_mass_before_L1Corr_"+Mgr.Procname[pr], h_bkg_mass_before_L1Corr[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_L1Corr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_L1Corr[iter]);
            f_bkg->GetObject("h_emu_mass_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_mass_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_emuSS_mass_before_TopPtCorr_"+Mgr.Procname[pr], h_SS_bkg_mass_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_emu_mass_fine_"+Mgr.Procname[pr], h_bkg_mass_fine[iter]);
            f_bkg->GetObject("h_emuSS_mass_fine_"+Mgr.Procname[pr], h_SS_bkg_mass_fine[iter]);
            f_bkg->GetObject("h_emu_mass_"+Mgr.Procname[pr], h_bkg_mass[iter]);
            f_bkg->GetObject("h_emuSS_mass_"+Mgr.Procname[pr], h_SS_bkg_mass[iter]);
            f_bkg->GetObject("h_emu_mass2_"+Mgr.Procname[pr], h_bkg_mass2[iter]);
            f_bkg->GetObject("h_emuSS_mass2_"+Mgr.Procname[pr], h_SS_bkg_mass2[iter]);
            removeNegativeBins(h_bkg_mass_before_PUCorr[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_RoccoR[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_mass_before_EffCorr[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_PVzCorr[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_mass_before_L1Corr[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_mass_before_TopPtCorr[iter]);
            removeNegativeBins(h_SS_bkg_mass_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_mass_fine[iter]);
            removeNegativeBins(h_SS_bkg_mass_fine[iter]);
            removeNegativeBins(h_bkg_mass[iter]);
            removeNegativeBins(h_SS_bkg_mass[iter]);
            removeNegativeBins(h_bkg_mass2[iter]);
            removeNegativeBins(h_SS_bkg_mass2[iter]);

            if (pr == _EMu_WJets_Full)
            {
                ratio_WJets_SSvsOS = h_bkg_mass[iter]->Integral(1, h_bkg_mass[iter]->GetSize()-2);
                ratio_WJets_SSvsOS /= h_SS_bkg_mass[iter]->Integral(1, h_SS_bkg_mass[iter]->GetSize()-2);

                for (Int_t i=1; i<h_bkg_mass[iter]->GetSize()-1; i++)
                {
                    if (h_SS_bkg_mass[iter]->GetBinContent(i) != 0)
                        avg_ratio_WJets_SSvsOS += h_bkg_mass[iter]->GetBinContent(i) / h_SS_bkg_mass[iter]->GetBinContent(i);
                }
                avg_ratio_WJets_SSvsOS /= h_SS_bkg_mass[iter]->GetSize()-2;
            }

            Color_t color = kBlack;
            if (pr == _EMu_WJets_Full) color = kRed - 2;
            if (pr == _EMu_WW) color = kMagenta - 5;
            if (pr == _EMu_WZ) color = kMagenta - 2;
            if (pr == _EMu_ZZ) color = kMagenta - 6;
            if (pr == _EMu_tbarW) color = kGreen - 2;
            if (pr == _EMu_tW) color = kGreen + 2;
            if (pr == _EMu_ttbar_Full) color = kCyan + 2;
            if (pr == _EMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_mass_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PUCorr[iter]->SetDirectory(0);
            s_mass_before_PUCorr->Add(h_bkg_mass_before_PUCorr[iter]);
            h_SS_bkg_mass_before_PUCorr[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_PUCorr[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_PUCorr[iter]->SetDirectory(0);
            s_SS_mass_before_PUCorr->Add(h_SS_bkg_mass_before_PUCorr[iter]);

            h_bkg_mass_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_mass_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_mass_before_RoccoR[iter]->SetDirectory(0);
            s_mass_before_RoccoR->Add(h_bkg_mass_before_RoccoR[iter]);
            h_SS_bkg_mass_before_RoccoR[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_RoccoR[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_RoccoR[iter]->SetDirectory(0);
            s_SS_mass_before_RoccoR->Add(h_SS_bkg_mass_before_RoccoR[iter]);

            h_bkg_mass_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
            s_mass_before_EffCorr->Add(h_bkg_mass_before_EffCorr[iter]);
            h_SS_bkg_mass_before_EffCorr[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_EffCorr[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
            s_SS_mass_before_EffCorr->Add(h_SS_bkg_mass_before_EffCorr[iter]);

            h_bkg_mass_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_PVzCorr[iter]->SetDirectory(0);
            s_mass_before_PVzCorr->Add(h_bkg_mass_before_PVzCorr[iter]);
            h_SS_bkg_mass_before_PVzCorr[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_PVzCorr[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_PVzCorr[iter]->SetDirectory(0);
            s_SS_mass_before_PVzCorr->Add(h_SS_bkg_mass_before_PVzCorr[iter]);

            h_bkg_mass_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_mass_before_L1Corr[iter]->SetDirectory(0);
            s_mass_before_L1Corr->Add(h_bkg_mass_before_L1Corr[iter]);
            h_SS_bkg_mass_before_L1Corr[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_L1Corr[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_L1Corr[iter]->SetDirectory(0);
            s_SS_mass_before_L1Corr->Add(h_SS_bkg_mass_before_L1Corr[iter]);

            h_bkg_mass_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_TopPtCorr[iter]->SetDirectory(0);
            s_mass_before_TopPtCorr->Add(h_bkg_mass_before_TopPtCorr[iter]);
            h_SS_bkg_mass_before_TopPtCorr[iter]->SetFillColor(color);
            h_SS_bkg_mass_before_TopPtCorr[iter]->SetLineColor(color);
            h_SS_bkg_mass_before_TopPtCorr[iter]->SetDirectory(0);
            s_SS_mass_before_TopPtCorr->Add(h_SS_bkg_mass_before_TopPtCorr[iter]);

            h_bkg_mass_fine[iter]->SetFillColor(color);
            h_bkg_mass_fine[iter]->SetLineColor(color);
            h_bkg_mass_fine[iter]->SetDirectory(0);
            s_mass_fine->Add(h_bkg_mass_fine[iter]);
            h_SS_bkg_mass_fine[iter]->SetFillColor(color);
            h_SS_bkg_mass_fine[iter]->SetLineColor(color);
            h_SS_bkg_mass_fine[iter]->SetDirectory(0);
            s_SS_mass_fine->Add(h_SS_bkg_mass_fine[iter]);

            h_bkg_mass[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetDirectory(0);
            s_mass->Add(h_bkg_mass[iter]);
            h_SS_bkg_mass[iter]->SetFillColor(color);
            h_SS_bkg_mass[iter]->SetLineColor(color);
            h_SS_bkg_mass[iter]->SetDirectory(0);
            s_SS_mass->Add(h_SS_bkg_mass[iter]);

            h_bkg_mass2[iter]->SetFillColor(color);
            h_bkg_mass2[iter]->SetLineColor(color);
            h_bkg_mass2[iter]->SetDirectory(0);
            s_mass2->Add(h_bkg_mass2[iter]);
            h_SS_bkg_mass2[iter]->SetFillColor(color);
            h_SS_bkg_mass2[iter]->SetLineColor(color);
            h_SS_bkg_mass2[iter]->SetDirectory(0);
            s_SS_mass2->Add(h_SS_bkg_mass2[iter]);

            iter++;

            if (pr == _EMu_WJets_Full) // next -- WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) // next - ttbar
                pr = _EMu_VVnST;
            if (pr == _EMu_DYTauTau_Full) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass_before_RoccoR, *h_data_mass_fine, *h_data_mass, *h_data_mass2,
             *h_SS_data_mass_before_RoccoR, *h_SS_data_mass_fine, *h_SS_data_mass, *h_SS_data_mass2;

        f_data->GetObject("h_emu_mass_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_before_RoccoR);
        f_data->GetObject("h_emuSS_mass_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_before_RoccoR);
        f_data->GetObject("h_emu_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass_fine);
        f_data->GetObject("h_emuSS_mass_fine_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass_fine);
        f_data->GetObject("h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass);
        f_data->GetObject("h_emuSS_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass);
        f_data->GetObject("h_emu_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_mass2);
        f_data->GetObject("h_emuSS_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_SS_data_mass2);

        h_data_mass_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_mass_before_RoccoR->SetMarkerColor(kBlack);
        h_data_mass_before_RoccoR->SetLineColor(kBlack);
        h_data_mass_before_RoccoR->SetDirectory(0);
        h_SS_data_mass_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_before_RoccoR->SetMarkerColor(kBlack);
        h_SS_data_mass_before_RoccoR->SetLineColor(kBlack);
        h_SS_data_mass_before_RoccoR->SetDirectory(0);

        h_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_data_mass_fine->SetMarkerColor(kBlack);
        h_data_mass_fine->SetLineColor(kBlack);
        h_data_mass_fine->SetDirectory(0);
        h_SS_data_mass_fine->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass_fine->SetMarkerColor(kBlack);
        h_SS_data_mass_fine->SetLineColor(kBlack);
        h_SS_data_mass_fine->SetDirectory(0);

        h_data_mass->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerColor(kBlack);
        h_data_mass->SetLineColor(kBlack);
        h_data_mass->SetDirectory(0);
        h_SS_data_mass->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass->SetMarkerColor(kBlack);
        h_SS_data_mass->SetLineColor(kBlack);
        h_SS_data_mass->SetDirectory(0);

        h_data_mass2->SetMarkerStyle(kFullDotLarge);
        h_data_mass2->SetMarkerColor(kBlack);
        h_data_mass2->SetLineColor(kBlack);
        h_data_mass2->SetDirectory(0);
        h_SS_data_mass2->SetMarkerStyle(kFullDotLarge);
        h_SS_data_mass2->SetMarkerColor(kBlack);
        h_SS_data_mass2->SetLineColor(kBlack);
        h_SS_data_mass2->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_before_PUCorr, *RP_mass_before_RoccoR, *RP_mass_before_EffCorr,
                      *RP_mass_before_PVzCorr, *RP_mass_before_L1Corr, *RP_mass_before_TopPtCorr,
                      *RP_mass_fine, *RP_mass, *RP_mass2,
                      *RP_SS_mass_before_PUCorr, *RP_SS_mass_before_RoccoR, *RP_SS_mass_before_EffCorr,
                      *RP_SS_mass_before_PVzCorr, *RP_SS_mass_before_L1Corr,*RP_SS_mass_before_TopPtCorr,
                      *RP_SS_mass_fine, *RP_SS_mass, *RP_SS_mass2;

        RP_mass_before_PUCorr = new myRatioPlot_t("RP_mass_before_PUCorr", s_mass_before_PUCorr, h_data_mass_before_RoccoR);
        RP_SS_mass_before_PUCorr = new myRatioPlot_t("RP_SS_mass_before_PUCorr", s_SS_mass_before_PUCorr, h_SS_data_mass_before_RoccoR);
        RP_mass_before_RoccoR = new myRatioPlot_t("RP_mass_before_RoccoR", s_mass_before_RoccoR, h_data_mass_before_RoccoR);
        RP_SS_mass_before_RoccoR = new myRatioPlot_t("RP_SS_mass_before_RoccoR", s_SS_mass_before_RoccoR, h_SS_data_mass_before_RoccoR);
        RP_mass_before_EffCorr = new myRatioPlot_t("RP_mass_before_EffCorr", s_mass_before_EffCorr, h_data_mass);
        RP_SS_mass_before_EffCorr = new myRatioPlot_t("RP_SS_mass_before_EffCorr", s_SS_mass_before_EffCorr, h_SS_data_mass);
        RP_mass_before_PVzCorr = new myRatioPlot_t("RP_mass_before_PVzCorr", s_mass_before_PVzCorr, h_data_mass);
        RP_SS_mass_before_PVzCorr = new myRatioPlot_t("RP_SS_mass_before_PVzCorr", s_SS_mass_before_PVzCorr, h_SS_data_mass);
        RP_mass_before_L1Corr = new myRatioPlot_t("RP_mass_before_L1Corr", s_mass_before_L1Corr, h_data_mass);
        RP_SS_mass_before_L1Corr = new myRatioPlot_t("RP_SS_mass_before_L1Corr", s_SS_mass_before_L1Corr, h_SS_data_mass);
        RP_mass_before_TopPtCorr = new myRatioPlot_t("RP_mass_before_TopPtCorr", s_mass_before_TopPtCorr, h_data_mass);
        RP_SS_mass_before_TopPtCorr = new myRatioPlot_t("RP_SS_mass_before_TopPtCorr", s_SS_mass_before_TopPtCorr, h_SS_data_mass);
        RP_mass_fine = new myRatioPlot_t("RP_mass_fine", s_mass_fine, h_data_mass_fine);
        RP_SS_mass_fine = new myRatioPlot_t("RP_SS_mass_fine", s_SS_mass_fine, h_SS_data_mass_fine);
        RP_mass = new myRatioPlot_t("RP_mass", s_mass, h_data_mass);
        RP_SS_mass = new myRatioPlot_t("RP_SS_mass", s_SS_mass, h_SS_data_mass);
        RP_mass2 = new myRatioPlot_t("RP_mass2", s_mass2, h_data_mass2);
        RP_SS_mass2 = new myRatioPlot_t("RP_SS_mass2", s_SS_mass2, h_SS_data_mass2);

        RP_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}] before PU correction", 15, 3000);
        RP_SS_mass_before_PUCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before PU correction", 15, 3000);
        RP_mass_before_RoccoR->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}] before Rochester correction", 15, 3000);
        RP_SS_mass_before_RoccoR->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before Rochester correction", 15, 3000);
        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}] before Efficiency SF", 15, 3000);
        RP_SS_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before Efficiency SF", 15, 3000);
        RP_mass_before_PVzCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}] before PVz correction", 15, 3000);
        RP_SS_mass_before_PVzCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before PVz correction", 15, 3000);
        RP_mass_before_L1Corr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}] before L1 prefiring correction", 15, 3000);
        RP_SS_mass_before_L1Corr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before L1 prefiring correction", 15, 3000);
        RP_mass_before_TopPtCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}] before top p_{#lower[-0.2]{T}} reweighting", 15, 3000);
        RP_SS_mass_before_TopPtCorr->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}] before top p_{#lower[-0.2]{T}} reweighting", 15, 3000);
        RP_mass_fine->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}]", 60, 120);
        RP_SS_mass_fine->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}]", 60, 120);
        RP_mass->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
        RP_SS_mass->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}]", 15, 3000);
        RP_mass2->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
        RP_SS_mass2->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e}#mu}}} (same-sign) [GeV/c^{2}]", 15, 3000);

        TLegend *legend = new TLegend(0.8, 0.45, 0.95, 0.95);

        legend->AddEntry(h_SS_data_mass, "Data", "lp");
        legend->AddEntry(h_bkg_mass[6+isWJ], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_SS_bkg_mass[5+isWJ], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_SS_bkg_mass[4+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_SS_bkg_mass[3+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_SS_bkg_mass[2+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_SS_bkg_mass[1+isWJ], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_SS_bkg_mass[0+isWJ], "#font[12]{#scale[1.1]{WW}}", "f");
        if (isWJ)
        {
            legend->AddEntry(h_SS_bkg_mass[0], "#font[12]{#scale[1.1]{W}}+Jets", "f");
        }

        RP_mass_before_PUCorr->ImportLegend(legend);
        RP_SS_mass_before_PUCorr->ImportLegend(legend);
        RP_mass_before_RoccoR->ImportLegend(legend);
        RP_SS_mass_before_RoccoR->ImportLegend(legend);
        RP_mass_before_EffCorr->ImportLegend(legend);
        RP_SS_mass_before_EffCorr->ImportLegend(legend);
        RP_mass_before_PVzCorr->ImportLegend(legend);
        RP_SS_mass_before_PVzCorr->ImportLegend(legend);
        RP_mass_before_L1Corr->ImportLegend(legend);
        RP_SS_mass_before_L1Corr->ImportLegend(legend);
        RP_mass_before_TopPtCorr->ImportLegend(legend);
        RP_SS_mass_before_TopPtCorr->ImportLegend(legend);
        RP_mass_fine->ImportLegend(legend);
        RP_SS_mass_fine->ImportLegend(legend);
        RP_mass->ImportLegend(legend);
        RP_SS_mass->ImportLegend(legend);
        RP_mass2->ImportLegend(legend);
        RP_SS_mass2->ImportLegend(legend);

        RP_mass_before_PUCorr->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_PUCorr->Draw(4e-1, 7e4, 1);
        RP_mass_before_RoccoR->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_RoccoR->Draw(4e-1, 7e4, 1);
        RP_mass_before_EffCorr->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_EffCorr->Draw(4e-1, 7e4, 1);
        RP_mass_before_PVzCorr->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_PVzCorr->Draw(4e-1, 7e4, 1);
        RP_mass_before_L1Corr->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_L1Corr->Draw(4e-1, 7e4, 1);
        RP_mass_before_TopPtCorr->Draw(4e-1, 7e4, 1);
        RP_SS_mass_before_TopPtCorr->Draw(4e-1,7e4, 1);
        RP_mass_fine->Draw(4e-1, 7e4, 0);
        RP_SS_mass_fine->Draw(4e-1, 7e4, 0);
        RP_mass->Draw(4e-1, 7e4, 1);
        RP_SS_mass->Draw(4e-1, 7e4, 1);
        RP_mass2->Draw(4e-1, 7e4, 1);
        RP_SS_mass2->Draw(4e-1, 7e4, 1);

//------------------------------ Get FR estimates -----------------------------------

        TH1D *h_QCDest_mass, *h_QCDest_SS_mass, *h_WJETest_mass, *h_WJETest_SS_mass;
        TFile *f_QCDFR = new TFile(Mgr.HistLocation+"EstQCD_EMu.root", "READ");
        f_QCDFR->GetObject("h_QCD_est", h_QCDest_mass);
        f_QCDFR->GetObject("h_QCD_est_SS", h_QCDest_SS_mass);
        h_QCDest_mass->SetDirectory(0);
        h_QCDest_mass->SetFillColor(kRed+3);
        h_QCDest_mass->SetLineColor(kRed+3);
        h_QCDest_SS_mass->SetDirectory(0);
        h_QCDest_SS_mass->SetFillColor(kRed+3);
        h_QCDest_SS_mass->SetLineColor(kRed+3);
        f_QCDFR->Close();

        TFile *f_WJETFR = new TFile(Mgr.HistLocation+"EstWJets_EMu.root", "READ");
        f_WJETFR->GetObject("h_WJET_est", h_WJETest_mass);
        f_WJETFR->GetObject("h_WJET_est_SS", h_WJETest_SS_mass);
        h_WJETest_mass->SetDirectory(0);
        h_WJETest_mass->SetFillColor(kRed-2);
        h_WJETest_mass->SetLineColor(kRed-2);
        h_WJETest_SS_mass->SetDirectory(0);
        h_WJETest_SS_mass->SetFillColor(kRed-2);
        h_WJETest_SS_mass->SetLineColor(kRed-2);
        f_WJETFR->Close();

        THStack *s_mass_wFR = new THStack("s_mass_wFR", "");
        s_mass_wFR->Add(h_QCDest_mass);
        s_mass_wFR->Add(h_WJETest_mass);
        s_mass_wFR->Add(h_bkg_mass[0+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[1+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[2+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[3+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[4+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[5+isWJ]);
        s_mass_wFR->Add(h_bkg_mass[6+isWJ]);
        THStack *s_mass_SS_wFR = new THStack("s_mass_SS_wFR", "");
        s_mass_SS_wFR->Add(h_QCDest_SS_mass);
        s_mass_SS_wFR->Add(h_WJETest_SS_mass);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[0+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[1+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[2+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[3+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[4+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[5+isWJ]);
        s_mass_SS_wFR->Add(h_SS_bkg_mass[6+isWJ]);

        TLegend *legend_wFR = new TLegend(0.8, 0.45, 0.95, 0.95);
        legend_wFR->AddEntry(h_data_mass, "Data", "lp");
        legend_wFR->AddEntry(h_bkg_mass[6+isWJ], "DY#rightarrow #tau#tau (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[5+isWJ], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[4+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[3+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[2+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[1+isWJ], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
        legend_wFR->AddEntry(h_bkg_mass[0+isWJ], "#font[12]{#scale[1.1]{WW}} (MC)", "f");
        legend_wFR->AddEntry(h_WJETest_mass, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
        legend_wFR->AddEntry(h_QCDest_mass, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");

        myRatioPlot_t *RP_mass_wFR = new myRatioPlot_t("RP_mass_wFR", s_mass_wFR, h_data_mass);
        myRatioPlot_t *RP_mass_SS_wFR = new myRatioPlot_t("RP_mass_SS_wFR", s_mass_SS_wFR, h_SS_data_mass);
        RP_mass_wFR->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
        RP_mass_SS_wFR->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
        RP_mass_wFR->ImportLegend(legend_wFR);
        RP_mass_SS_wFR->ImportLegend(legend_wFR);
        RP_mass_wFR->Draw(4e-1, 7e4, 1);
        RP_mass_SS_wFR->Draw(4e-1, 7e4, 1);

// --------------------------- Printing some numbers --------------------------------

        Double_t dataerror, MCerror, MCerror_noSF, dataintegral=348650, MCintegral, MCintegral_noSF;
        Double_t dataerrorSS, MCerrorSS, MCerrorSS_noSF, dataintegralSS=348650, MCintegralSS, MCintegralSS_noSF;
        Double_t dataerrorZ, MCerrorZ, dataintegralZ=2.25081e+07, MCintegralZ;
        Double_t dataerror_noZ=0, MCerror_noZ=0, dataintegral_noZ=2.25081e+07, MCintegral_noZ, temp_noZ;

        dataintegral = h_data_mass->IntegralAndError(1, h_data_mass->GetSize()-2, dataerror);
        MCintegral = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror);
        MCintegral_noSF = ((TH1D*)(s_mass_before_RoccoR->GetStack()->Last()))->IntegralAndError(1, h_data_mass->GetSize()-2, MCerror_noSF);
        dataintegralSS = h_SS_data_mass->IntegralAndError(1, h_SS_data_mass->GetSize()-2, dataerrorSS);
        MCintegralSS = ((TH1D*)(s_SS_mass->GetStack()->Last()))->IntegralAndError(1, h_SS_data_mass->GetSize()-2, MCerrorSS);
        MCintegralSS_noSF = ((TH1D*)(s_SS_mass_before_RoccoR->GetStack()->Last()))->IntegralAndError(1, h_SS_data_mass->GetSize()-2, MCerrorSS_noSF);

        dataintegralZ = h_data_mass->IntegralAndError(10, 22, dataerrorZ);
        MCintegralZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ);

        dataintegral_noZ = h_data_mass->IntegralAndError(1, 9, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataintegral_noZ += h_data_mass->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        dataerror_noZ += temp_noZ * temp_noZ;
        dataerror_noZ = sqrt(dataerror_noZ);

        MCintegral_noZ = ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCintegral_noZ += ((TH1D*)(s_mass->GetStack()->Last()))->IntegralAndError(23, h_data_mass->GetSize()-2, temp_noZ);
        MCerror_noZ += temp_noZ * temp_noZ;
        MCerror_noZ = sqrt(MCerror_noZ);

        std::cout << "Data events: " << dataintegral << "+-" << dataerror << endl;
        std::cout << "MC events: " << MCintegral << "+-" << MCerror << endl;
        std::cout << "MC/Obs: " << MCintegral / dataintegral << "+-" <<
                     sqrt((dataerror / dataintegral) * (dataerror / dataintegral) +
                           (MCerror / MCintegral) * (MCerror / MCintegral)) << endl;
        std::cout << "MC events before corrections: " << MCintegral_noSF << "+-" << MCerror_noSF << endl;
        std::cout << "Avg. Data and MC relative difference: " << CompAvgDataMCDifference(h_data_mass, s_mass) << endl;

        std::cout << "Same-sign data events: " << dataintegralSS << "+-" << dataerrorSS << endl;
        std::cout << "Same-sign MC events: " << MCintegralSS << "+-" << MCerrorSS << endl;
        std::cout << "Same-sign MC/Obs: " << MCintegralSS / dataintegralSS << "+-" <<
                     sqrt((dataerrorSS / dataintegralSS) * (dataerrorSS / dataintegralSS) +
                           (MCerrorSS / MCintegralSS) * (MCerrorSS / MCintegralSS)) << endl;
        std::cout << "Same-sign MC events before corrections: " << MCintegralSS_noSF << "+-" << MCerrorSS_noSF << endl << endl;

        std::cout << "Data events around Z: " << dataintegralZ << "+-" << dataerrorZ << endl;
        std::cout << "MC events around Z: " << MCintegralZ << "+-" << MCerrorZ << endl;
        std::cout << "Data events outside Z: " << dataintegral_noZ << "+-" << dataerror_noZ << endl;
        std::cout << "MC events outside Z: " << MCintegral_noZ << "+-" << MCerror_noZ << endl;

        if (isWJ)
        {
            std::cout << "OS WJets events: " << h_bkg_mass[0]->Integral(1, h_bkg_mass[0]->GetSize()-2) << endl;
            std::cout << "SS WJets events: " << h_SS_bkg_mass[0]->Integral(1, h_SS_bkg_mass[0]->GetSize()-2) << endl;
            std::cout << "OS/SS ratio of WJets events: " << ratio_WJets_SSvsOS << endl;
            std::cout << "Average OS/SS ratio of WJets events per bin: " << avg_ratio_WJets_SSvsOS << endl << endl;
        }

    } // End of if(invm)

//################################# Pt, RAPI, pT, ETA, PHI #################################################
/*
    if (whichGraphs=="ALL" || whichGraphs.Contains("PT") || whichGraphs.Contains("ETA") || whichGraphs.Contains("PHI") || whichGraphs.Contains("RAPI"))
    {
        count_drawn++;

        THStack *s_pT_ele_before_PUCorr = new THStack("s_pT_ele_before_PUCorr", "");
        THStack *s_pT_mu_before_PUCorr = new THStack("s_pT_mu_before_PUCorr", "");
        THStack *s_pT_eleSS_before_PUCorr = new THStack("s_pT_eleSS_before_PUCorr", "");
        THStack *s_pT_muSS_before_PUCorr = new THStack("s_pT_muSS_before_PUCorr", "");
        THStack *s_eta_ele_before_PUCorr = new THStack("s_eta_ele_before_PUCorr", "");
        THStack *s_eta_mu_before_PUCorr = new THStack("s_eta_mu_before_PUCorr", "");
        THStack *s_eta_eleSS_before_PUCorr = new THStack("s_eta_eleSS_before_PUCorr", "");
        THStack *s_eta_muSS_before_PUCorr = new THStack("s_eta_muSS_before_PUCorr", "");
        THStack *s_phi_ele_before_PUCorr = new THStack("s_phi_ele_before_PUCorr", "");
        THStack *s_phi_mu_before_PUCorr = new THStack("s_phi_mu_before_PUCorr", "");
        THStack *s_phi_eleSS_before_PUCorr = new THStack("s_phi_eleSS_before_PUCorr", "");
        THStack *s_phi_muSS_before_PUCorr = new THStack("s_phi_muSS_before_PUCorr", "");

        THStack *s_pT_ele_before_RoccoR = new THStack("s_pT_ele_before_RoccoR", "");
        THStack *s_pT_mu_before_RoccoR = new THStack("s_pT_mu_before_RoccoR", "");
        THStack *s_pT_eleSS_before_RoccoR = new THStack("s_pT_eleSS_before_RoccoR", "");
        THStack *s_pT_muSS_before_RoccoR = new THStack("s_pT_muSS_before_RoccoR", "");
        THStack *s_eta_ele_before_RoccoR = new THStack("s_eta_ele_before_RoccoR", "");
        THStack *s_eta_mu_before_RoccoR = new THStack("s_eta_mu_before_RoccoR", "");
        THStack *s_eta_eleSS_before_RoccoR = new THStack("s_eta_eleSS_before_RoccoR", "");
        THStack *s_eta_muSS_before_RoccoR = new THStack("s_eta_muSS_before_RoccoR", "");
        THStack *s_phi_ele_before_RoccoR = new THStack("s_phi_ele_before_RoccoR", "");
        THStack *s_phi_mu_before_RoccoR = new THStack("s_phi_mu_before_RoccoR", "");
        THStack *s_phi_eleSS_before_RoccoR = new THStack("s_phi_eleSS_before_RoccoR", "");
        THStack *s_phi_muSS_before_RoccoR = new THStack("s_phi_muSS_before_RoccoR", "");

        THStack *s_pT_ele_before_EffCorr = new THStack("s_pT_ele_before_EffCorr", "");
        THStack *s_pT_mu_before_EffCorr = new THStack("s_pT_mu_before_EffCorr", "");
        THStack *s_pT_eleSS_before_EffCorr = new THStack("s_pT_eleSS_before_EffCorr", "");
        THStack *s_pT_muSS_before_EffCorr = new THStack("s_pT_muSS_before_EffCorr", "");
        THStack *s_eta_ele_before_EffCorr = new THStack("s_eta_ele_before_EffCorr", "");
        THStack *s_eta_mu_before_EffCorr = new THStack("s_eta_mu_before_EffCorr", "");
        THStack *s_eta_eleSS_before_EffCorr = new THStack("s_eta_eleSS_before_EffCorr", "");
        THStack *s_eta_muSS_before_EffCorr = new THStack("s_eta_muSS_before_EffCorr", "");
        THStack *s_phi_ele_before_EffCorr = new THStack("s_phi_ele_before_EffCorr", "");
        THStack *s_phi_mu_before_EffCorr = new THStack("s_phi_mu_before_EffCorr", "");
        THStack *s_phi_eleSS_before_EffCorr = new THStack("s_phi_eleSS_before_EffCorr", "");
        THStack *s_phi_muSS_before_EffCorr = new THStack("s_phi_muSS_before_EffCorr", "");

        THStack *s_pT_ele_before_PVzCorr = new THStack("s_pT_ele_before_PVzCorr", "");
        THStack *s_pT_mu_before_PVzCorr = new THStack("s_pT_mu_before_PVzCorr", "");
        THStack *s_pT_eleSS_before_PVzCorr = new THStack("s_pT_eleSS_before_PVzCorr", "");
        THStack *s_pT_muSS_before_PVzCorr = new THStack("s_pT_muSS_before_PVzCorr", "");
        THStack *s_eta_ele_before_PVzCorr = new THStack("s_eta_ele_before_PVzCorr", "");
        THStack *s_eta_mu_before_PVzCorr = new THStack("s_eta_mu_before_PVzCorr", "");
        THStack *s_eta_eleSS_before_PVzCorr = new THStack("s_eta_eleSS_before_PVzCorr", "");
        THStack *s_eta_muSS_before_PVzCorr = new THStack("s_eta_muSS_before_PVzCorr", "");
        THStack *s_phi_ele_before_PVzCorr = new THStack("s_phi_ele_before_PVzCorr", "");
        THStack *s_phi_mu_before_PVzCorr = new THStack("s_phi_mu_before_PVzCorr", "");
        THStack *s_phi_eleSS_before_PVzCorr = new THStack("s_phi_eleSS_before_PVzCorr", "");
        THStack *s_phi_muSS_before_PVzCorr = new THStack("s_phi_muSS_before_PVzCorr", "");

        THStack *s_pT_ele_before_L1Corr = new THStack("s_pT_ele_before_L1Corr", "");
        THStack *s_pT_mu_before_L1Corr = new THStack("s_pT_mu_before_L1Corr", "");
        THStack *s_pT_eleSS_before_L1Corr = new THStack("s_pT_eleSS_before_L1Corr", "");
        THStack *s_pT_muSS_before_L1Corr = new THStack("s_pT_muSS_before_L1Corr", "");
        THStack *s_eta_ele_before_L1Corr = new THStack("s_eta_ele_before_L1Corr", "");
        THStack *s_eta_mu_before_L1Corr = new THStack("s_eta_mu_before_L1Corr", "");
        THStack *s_eta_eleSS_before_L1Corr = new THStack("s_eta_eleSS_before_L1Corr", "");
        THStack *s_eta_muSS_before_L1Corr = new THStack("s_eta_muSS_before_L1Corr", "");
        THStack *s_phi_ele_before_L1Corr = new THStack("s_phi_ele_before_L1Corr", "");
        THStack *s_phi_mu_before_L1Corr = new THStack("s_phi_mu_before_L1Corr", "");
        THStack *s_phi_eleSS_before_L1Corr = new THStack("s_phi_eleSS_before_L1Corr", "");
        THStack *s_phi_muSS_before_L1Corr = new THStack("s_phi_muSS_before_L1Corr", "");

        THStack *s_pT_ele_before_TopPtCorr = new THStack("s_pT_ele_before_TopPtCorr", "");
        THStack *s_pT_mu_before_TopPtCorr = new THStack("s_pT_mu_before_TopPtCorr", "");
        THStack *s_pT_eleSS_before_TopPtCorr = new THStack("s_pT_eleSS_before_TopPtCorr", "");
        THStack *s_pT_muSS_before_TopPtCorr = new THStack("s_pT_muSS_before_TopPtCorr", "");
        THStack *s_eta_ele_before_TopPtCorr = new THStack("s_eta_ele_before_TopPtCorr", "");
        THStack *s_eta_mu_before_TopPtCorr = new THStack("s_eta_mu_before_TopPtCorr", "");
        THStack *s_eta_eleSS_before_TopPtCorr = new THStack("s_eta_eleSS_before_TopPtCorr", "");
        THStack *s_eta_muSS_before_TopPtCorr = new THStack("s_eta_muSS_before_TopPtCorr", "");
        THStack *s_phi_ele_before_TopPtCorr = new THStack("s_phi_ele_before_TopPtCorr", "");
        THStack *s_phi_mu_before_TopPtCorr = new THStack("s_phi_mu_before_TopPtCorr", "");
        THStack *s_phi_eleSS_before_TopPtCorr = new THStack("s_phi_eleSS_before_TopPtCorr", "");
        THStack *s_phi_muSS_before_TopPtCorr = new THStack("s_phi_muSS_before_TopPtCorr", "");

        THStack *s_pT_ele = new THStack("s_pT_ele", "");
        THStack *s_pT_mu = new THStack("s_pT_mu", "");
        THStack *s_pT_eleSS = new THStack("s_pT_eleSS", "");
        THStack *s_pT_muSS = new THStack("s_pT_muSS", "");
        THStack *s_eta_ele = new THStack("s_eta_ele", "");
        THStack *s_eta_mu = new THStack("s_eta_mu", "");
        THStack *s_eta_eleSS = new THStack("s_eta_eleSS", "");
        THStack *s_eta_muSS = new THStack("s_eta_muSS", "");
        THStack *s_phi_ele = new THStack("s_phi_ele", "");
        THStack *s_phi_mu = new THStack("s_phi_mu", "");
        THStack *s_phi_eleSS = new THStack("s_phi_eleSS", "");
        THStack *s_phi_muSS = new THStack("s_phi_muSS", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_pT_ele_before_PUCorr[8], *h_bkg_eta_ele_before_PUCorr[8], *h_bkg_phi_ele_before_PUCorr[8],
             *h_bkg_pT_mu_before_PUCorr[8],*h_bkg_eta_mu_before_PUCorr[8], *h_bkg_phi_mu_before_PUCorr[8],
             *h_bkg_pT_eleSS_before_PUCorr[8], *h_bkg_eta_eleSS_before_PUCorr[8], *h_bkg_phi_eleSS_before_PUCorr[8],
             *h_bkg_pT_muSS_before_PUCorr[8], *h_bkg_eta_muSS_before_PUCorr[8], *h_bkg_phi_muSS_before_PUCorr[8],
             *h_bkg_pT_ele_before_RoccoR[8], *h_bkg_eta_ele_before_RoccoR[8], *h_bkg_phi_ele_before_RoccoR[8],
             *h_bkg_pT_mu_before_RoccoR[8], *h_bkg_eta_mu_before_RoccoR[8], *h_bkg_phi_mu_before_RoccoR[8],
             *h_bkg_pT_eleSS_before_RoccoR[8], *h_bkg_eta_eleSS_before_RoccoR[8], *h_bkg_phi_eleSS_before_RoccoR[8],
             *h_bkg_pT_muSS_before_RoccoR[8], *h_bkg_eta_muSS_before_RoccoR[8], *h_bkg_phi_muSS_before_RoccoR[8],
             *h_bkg_pT_ele_before_EffCorr[8], *h_bkg_eta_ele_before_EffCorr[8], *h_bkg_phi_ele_before_EffCorr[8],
             *h_bkg_pT_mu_before_EffCorr[8], *h_bkg_eta_mu_before_EffCorr[8], *h_bkg_phi_mu_before_EffCorr[8],
             *h_bkg_pT_eleSS_before_EffCorr[8], *h_bkg_eta_eleSS_before_EffCorr[8], *h_bkg_phi_eleSS_before_EffCorr[8],
             *h_bkg_pT_muSS_before_EffCorr[8], *h_bkg_eta_muSS_before_EffCorr[8], *h_bkg_phi_muSS_before_EffCorr[8],
             *h_bkg_pT_ele_before_PVzCorr[8], *h_bkg_eta_ele_before_PVzCorr[8], *h_bkg_phi_ele_before_PVzCorr[8],
             *h_bkg_pT_mu_before_PVzCorr[8], *h_bkg_eta_mu_before_PVzCorr[8], *h_bkg_phi_mu_before_PVzCorr[8],
             *h_bkg_pT_eleSS_before_PVzCorr[8], *h_bkg_eta_eleSS_before_PVzCorr[8], *h_bkg_phi_eleSS_before_PVzCorr[8],
             *h_bkg_pT_muSS_before_PVzCorr[8], *h_bkg_eta_muSS_before_PVzCorr[8], *h_bkg_phi_muSS_before_PVzCorr[8],
             *h_bkg_pT_ele_before_L1Corr[8], *h_bkg_eta_ele_before_L1Corr[8], *h_bkg_phi_ele_before_L1Corr[8],
             *h_bkg_pT_mu_before_L1Corr[8], *h_bkg_eta_mu_before_L1Corr[8], *h_bkg_phi_mu_before_L1Corr[8],
             *h_bkg_pT_eleSS_before_L1Corr[8], *h_bkg_eta_eleSS_before_L1Corr[8], *h_bkg_phi_eleSS_before_L1Corr[8],
             *h_bkg_pT_muSS_before_L1Corr[8], *h_bkg_eta_muSS_before_L1Corr[8], *h_bkg_phi_muSS_before_L1Corr[8],
             *h_bkg_pT_ele_before_TopPtCorr[8], *h_bkg_eta_ele_before_TopPtCorr[8], *h_bkg_phi_ele_before_TopPtCorr[8],
             *h_bkg_pT_mu_before_TopPtCorr[8], *h_bkg_eta_mu_before_TopPtCorr[8], *h_bkg_phi_mu_before_TopPtCorr[8],
             *h_bkg_pT_eleSS_before_TopPtCorr[8], *h_bkg_eta_eleSS_before_TopPtCorr[8], *h_bkg_phi_eleSS_before_TopPtCorr[8],
             *h_bkg_pT_muSS_before_TopPtCorr[8], *h_bkg_eta_muSS_before_TopPtCorr[8], *h_bkg_phi_muSS_before_TopPtCorr[8],
             *h_bkg_pT_ele[8], *h_bkg_eta_ele[8], *h_bkg_phi_ele[8], *h_bkg_pT_mu[8], *h_bkg_eta_mu[8], *h_bkg_phi_mu[8],
             *h_bkg_pT_eleSS[8], *h_bkg_eta_eleSS[8], *h_bkg_phi_eleSS[8], *h_bkg_pT_muSS[8], *h_bkg_eta_muSS[8], *h_bkg_phi_muSS[8];
        Int_t iter = 0;

        for (SelProc_t pr = _EMu_WJets_Full; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (!isWJ && pr == _EMu_WJets_Full) continue;
            if (pr == _EndOf_EMu_VVnST_Normal) continue;

            f_bkg->GetObject("h_ele_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_ele_before_PUCorr[iter]);
            f_bkg->GetObject("h_ele_eta_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_ele_before_PUCorr[iter]);
            f_bkg->GetObject("h_ele_phi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_ele_before_PUCorr[iter]);
            f_bkg->GetObject("h_mu_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_PUCorr[iter]);
            f_bkg->GetObject("h_mu_eta_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_PUCorr[iter]);
            f_bkg->GetObject("h_mu_phi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_PUCorr[iter]);
            f_bkg->GetObject("h_eleSS_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_eleSS_before_PUCorr[iter]);
            f_bkg->GetObject("h_eleSS_eta_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_eleSS_before_PUCorr[iter]);
            f_bkg->GetObject("h_eleSS_phi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_eleSS_before_PUCorr[iter]);
            f_bkg->GetObject("h_muSS_pT_before_PUCorr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_PUCorr[iter]);
            f_bkg->GetObject("h_muSS_eta_before_PUCorr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_PUCorr[iter]);
            f_bkg->GetObject("h_muSS_phi_before_PUCorr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_ele_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_ele_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_ele_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_eleSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_eleSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_eleSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_PUCorr[iter]);

            f_bkg->GetObject("h_mu_pT_before_RocCorr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_RoccoR[iter]);
            f_bkg->GetObject("h_mu_eta_before_RocCorr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_RoccoR[iter]);
            f_bkg->GetObject("h_mu_phi_before_RocCorr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_RoccoR[iter]);
            f_bkg->GetObject("h_muSS_pT_before_RocCorr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_RoccoR[iter]);
            f_bkg->GetObject("h_muSS_eta_before_RocCorr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_RoccoR[iter]);
            f_bkg->GetObject("h_muSS_phi_before_RocCorr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_RoccoR[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_RoccoR[iter]);

            f_bkg->GetObject("h_ele_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_ele_before_EffCorr[iter]);
            f_bkg->GetObject("h_ele_eta_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_ele_before_EffCorr[iter]);
            f_bkg->GetObject("h_ele_phi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_ele_before_EffCorr[iter]);
            f_bkg->GetObject("h_mu_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_EffCorr[iter]);
            f_bkg->GetObject("h_mu_eta_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_EffCorr[iter]);
            f_bkg->GetObject("h_mu_phi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_EffCorr[iter]);
            f_bkg->GetObject("h_eleSS_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_eleSS_before_EffCorr[iter]);
            f_bkg->GetObject("h_eleSS_eta_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_eleSS_before_EffCorr[iter]);
            f_bkg->GetObject("h_eleSS_phi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_eleSS_before_EffCorr[iter]);
            f_bkg->GetObject("h_muSS_pT_before_EffCorr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_EffCorr[iter]);
            f_bkg->GetObject("h_muSS_eta_before_EffCorr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_EffCorr[iter]);
            f_bkg->GetObject("h_muSS_phi_before_EffCorr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_ele_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_ele_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_ele_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_eleSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_eleSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_eleSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_EffCorr[iter]);

            f_bkg->GetObject("h_ele_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_ele_before_PVzCorr[iter]);
            f_bkg->GetObject("h_ele_eta_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_ele_before_PVzCorr[iter]);
            f_bkg->GetObject("h_ele_phi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_ele_before_PVzCorr[iter]);
            f_bkg->GetObject("h_mu_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_PVzCorr[iter]);
            f_bkg->GetObject("h_mu_eta_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_PVzCorr[iter]);
            f_bkg->GetObject("h_mu_phi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eleSS_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_eleSS_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eleSS_eta_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_eleSS_before_PVzCorr[iter]);
            f_bkg->GetObject("h_eleSS_phi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_eleSS_before_PVzCorr[iter]);
            f_bkg->GetObject("h_muSS_pT_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_PVzCorr[iter]);
            f_bkg->GetObject("h_muSS_eta_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_PVzCorr[iter]);
            f_bkg->GetObject("h_muSS_phi_before_PVzCorr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_ele_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_ele_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_ele_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_eleSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_eleSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_eleSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_PVzCorr[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_PVzCorr[iter]);

            f_bkg->GetObject("h_ele_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_ele_before_L1Corr[iter]);
            f_bkg->GetObject("h_ele_eta_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_ele_before_L1Corr[iter]);
            f_bkg->GetObject("h_ele_phi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_ele_before_L1Corr[iter]);
            f_bkg->GetObject("h_mu_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_L1Corr[iter]);
            f_bkg->GetObject("h_mu_eta_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_L1Corr[iter]);
            f_bkg->GetObject("h_mu_phi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_L1Corr[iter]);
            f_bkg->GetObject("h_eleSS_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_eleSS_before_L1Corr[iter]);
            f_bkg->GetObject("h_eleSS_eta_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_eleSS_before_L1Corr[iter]);
            f_bkg->GetObject("h_eleSS_phi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_eleSS_before_L1Corr[iter]);
            f_bkg->GetObject("h_muSS_pT_before_L1Corr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_L1Corr[iter]);
            f_bkg->GetObject("h_muSS_eta_before_L1Corr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_L1Corr[iter]);
            f_bkg->GetObject("h_muSS_phi_before_L1Corr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_ele_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_ele_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_ele_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_eleSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_eleSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_eleSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_L1Corr[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_L1Corr[iter]);

            f_bkg->GetObject("h_ele_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_ele_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_ele_eta_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_ele_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_ele_phi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_ele_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_mu_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_mu_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_mu_eta_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_mu_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_mu_phi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_mu_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eleSS_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_eleSS_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eleSS_eta_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_eleSS_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_eleSS_phi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_eleSS_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_muSS_pT_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_pT_muSS_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_muSS_eta_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_eta_muSS_before_TopPtCorr[iter]);
            f_bkg->GetObject("h_muSS_phi_before_TopPtCorr_"+Mgr.Procname[pr], h_bkg_phi_muSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_ele_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_ele_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_ele_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_mu_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_mu_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_mu_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_eleSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_eleSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_eleSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_pT_muSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_eta_muSS_before_TopPtCorr[iter]);
            removeNegativeBins(h_bkg_phi_muSS_before_TopPtCorr[iter]);

            f_bkg->GetObject("h_ele_pT_"+Mgr.Procname[pr], h_bkg_pT_ele[iter]);
            f_bkg->GetObject("h_ele_eta_"+Mgr.Procname[pr], h_bkg_eta_ele[iter]);
            f_bkg->GetObject("h_ele_phi_"+Mgr.Procname[pr], h_bkg_phi_ele[iter]);
            f_bkg->GetObject("h_mu_pT_"+Mgr.Procname[pr], h_bkg_pT_mu[iter]);
            f_bkg->GetObject("h_mu_eta_"+Mgr.Procname[pr], h_bkg_eta_mu[iter]);
            f_bkg->GetObject("h_mu_phi_"+Mgr.Procname[pr], h_bkg_phi_mu[iter]);
            f_bkg->GetObject("h_eleSS_pT_"+Mgr.Procname[pr], h_bkg_pT_eleSS[iter]);
            f_bkg->GetObject("h_eleSS_eta_"+Mgr.Procname[pr], h_bkg_eta_eleSS[iter]);
            f_bkg->GetObject("h_eleSS_phi_"+Mgr.Procname[pr], h_bkg_phi_eleSS[iter]);
            f_bkg->GetObject("h_muSS_pT_"+Mgr.Procname[pr], h_bkg_pT_muSS[iter]);
            f_bkg->GetObject("h_muSS_eta_"+Mgr.Procname[pr], h_bkg_eta_muSS[iter]);
            f_bkg->GetObject("h_muSS_phi_"+Mgr.Procname[pr], h_bkg_phi_muSS[iter]);
            removeNegativeBins(h_bkg_pT_ele[iter]);
            removeNegativeBins(h_bkg_eta_el[iter]);
            removeNegativeBins(h_bkg_phi_ele[iter]);
            removeNegativeBins(h_bkg_pT_mu[iter]);
            removeNegativeBins(h_bkg_eta_mu[iter]);
            removeNegativeBins(h_bkg_phi_mu[iter]);
            removeNegativeBins(h_bkg_pT_eleSS[iter]);
            removeNegativeBins(h_bkg_eta_eleSS[iter]);
            removeNegativeBins(h_bkg_phi_eleSS[iter]);
            removeNegativeBins(h_bkg_pT_muSS[iter]);
            removeNegativeBins(h_bkg_eta_muSS[iter]);
            removeNegativeBins(h_bkg_phi_muSS[iter]);

            Color_t color = kBlack;
            if (pr == _EMu_WJets_Full) color = kRed - 2;
            if (pr == _EMu_WW) color = kMagenta - 5;
            if (pr == _EMu_WZ) color = kMagenta - 2;
            if (pr == _EMu_ZZ) color = kMagenta - 6;
            if (pr == _EMu_tbarW) color = kGreen - 2;
            if (pr == _EMu_tW) color = kGreen + 2;
            if (pr == _EMu_ttbar_Full) color = kCyan + 2;
            if (pr == _EMu_DYTauTau_Full) color = kOrange - 5;

            // Before PU correction
            h_bkg_pT_ele_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_ele_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_ele_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_mu_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_eleSS_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_eleSS_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_eleSS_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_PUCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_ele_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_ele_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_mu_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_eleSS_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_eleSS_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_eleSS_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_PUCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_ele_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_ele_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_ele_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_mu_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_eleSS_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_eleSS_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_eleSS_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_PUCorr[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_PUCorr[iter]->SetDirectory(0);
            //
            s_pT_ele_before_PUCorr->Add(h_bkg_pT_ele_before_PUCorr[iter]);
            s_eta_ele_before_PUCorr->Add(h_bkg_eta_ele_before_PUCorr[iter]);
            s_phi_ele_before_PUCorr->Add(h_bkg_phi_ele_before_PUCorr[iter]);
            s_pT_mu_before_PUCorr->Add(h_bkg_pT_mu_before_PUCorr[iter]);
            s_eta_mu_before_PUCorr->Add(h_bkg_eta_mu_before_PUCorr[iter]);
            s_phi_mu_before_PUCorr->Add(h_bkg_phi_mu_before_PUCorr[iter]);
            s_pT_eleSS_before_PUCorr->Add(h_bkg_pT_eleSS_before_PUCorr[iter]);
            s_eta_eleSS_before_PUCorr->Add(h_bkg_eta_eleSS_before_PUCorr[iter]);
            s_phi_eleSS_before_PUCorr->Add(h_bkg_phi_eleSS_before_PUCorr[iter]);
            s_pT_muSS_before_PUCorr->Add(h_bkg_pT_muSS_before_PUCorr[iter]);
            s_eta_muSS_before_PUCorr->Add(h_bkg_eta_muSS_before_PUCorr[iter]);
            s_phi_muSS_before_PUCorr->Add(h_bkg_phi_muSS_before_PUCorr[iter]);

            // Before Rochester correction
            h_bkg_pT_mu_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_RoccoR[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_RoccoR[iter]->SetFillColor(color);
            //
            h_bkg_pT_mu_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_RoccoR[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_RoccoR[iter]->SetLineColor(color);
            //
            h_bkg_pT_mu_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_RoccoR[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_RoccoR[iter]->SetDirectory(0);
            //
            s_pT_mu_before_RoccoR->Add(h_bkg_pT_mu_before_RoccoR[iter]);
            s_eta_mu_before_RoccoR->Add(h_bkg_eta_mu_before_RoccoR[iter]);
            s_phi_mu_before_RoccoR->Add(h_bkg_phi_mu_before_RoccoR[iter]);
            s_pT_muSS_before_RoccoR->Add(h_bkg_pT_muSS_before_RoccoR[iter]);
            s_eta_muSS_before_RoccoR->Add(h_bkg_eta_muSS_before_RoccoR[iter]);
            s_phi_muSS_before_RoccoR->Add(h_bkg_phi_muSS_before_RoccoR[iter]);

            // Before efficiency SF
            h_bkg_pT_ele_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_ele_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_ele_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_mu_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_eleSS_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_eleSS_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_eleSS_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_EffCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_ele_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_ele_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_mu_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_eleSS_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_eleSS_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_eleSS_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_EffCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_ele_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_ele_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_ele_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_mu_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_eleSS_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_eleSS_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_eleSS_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_EffCorr[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_EffCorr[iter]->SetDirectory(0);
            //
            s_pT_ele_before_EffCorr->Add(h_bkg_pT_ele_before_EffCorr[iter]);
            s_eta_ele_before_EffCorr->Add(h_bkg_eta_ele_before_EffCorr[iter]);
            s_phi_ele_before_EffCorr->Add(h_bkg_phi_ele_before_EffCorr[iter]);
            s_pT_mu_before_EffCorr->Add(h_bkg_pT_mu_before_EffCorr[iter]);
            s_eta_mu_before_EffCorr->Add(h_bkg_eta_mu_before_EffCorr[iter]);
            s_phi_mu_before_EffCorr->Add(h_bkg_phi_mu_before_EffCorr[iter]);
            s_pT_eleSS_before_EffCorr->Add(h_bkg_pT_eleSS_before_EffCorr[iter]);
            s_eta_eleSS_before_EffCorr->Add(h_bkg_eta_eleSS_before_EffCorr[iter]);
            s_phi_eleSS_before_EffCorr->Add(h_bkg_phi_eleSS_before_EffCorr[iter]);
            s_pT_muSS_before_EffCorr->Add(h_bkg_pT_muSS_before_EffCorr[iter]);
            s_eta_muSS_before_EffCorr->Add(h_bkg_eta_muSS_before_EffCorr[iter]);
            s_phi_muSS_before_EffCorr->Add(h_bkg_phi_muSS_before_EffCorr[iter]);

            // Before PVz reweighting
            h_bkg_pT_ele_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_ele_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_ele_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_mu_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_eleSS_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_eleSS_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_eleSS_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_PVzCorr[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_PVzCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_ele_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_ele_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_mu_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_eleSS_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_eleSS_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_eleSS_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_PVzCorr[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_PVzCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_ele_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_ele_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_ele_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_mu_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_eleSS_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_eleSS_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_eleSS_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_PVzCorr[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_PVzCorr[iter]->SetDirectory(0);
            //
            s_pT_ele_before_PVzCorr->Add(h_bkg_pT_ele_before_PVzCorr[iter]);
            s_eta_ele_before_PVzCorr->Add(h_bkg_eta_ele_before_PVzCorr[iter]);
            s_phi_ele_before_PVzCorr->Add(h_bkg_phi_ele_before_PVzCorr[iter]);
            s_pT_mu_before_PVzCorr->Add(h_bkg_pT_mu_before_PVzCorr[iter]);
            s_eta_mu_before_PVzCorr->Add(h_bkg_eta_mu_before_PVzCorr[iter]);
            s_phi_mu_before_PVzCorr->Add(h_bkg_phi_mu_before_PVzCorr[iter]);
            s_pT_eleSS_before_PVzCorr->Add(h_bkg_pT_eleSS_before_PVzCorr[iter]);
            s_eta_eleSS_before_PVzCorr->Add(h_bkg_eta_eleSS_before_PVzCorr[iter]);
            s_phi_eleSS_before_PVzCorr->Add(h_bkg_phi_eleSS_before_PVzCorr[iter]);
            s_pT_muSS_before_PVzCorr->Add(h_bkg_pT_muSS_before_PVzCorr[iter]);
            s_eta_muSS_before_PVzCorr->Add(h_bkg_eta_muSS_before_PVzCorr[iter]);
            s_phi_muSS_before_PVzCorr->Add(h_bkg_phi_muSS_before_PVzCorr[iter]);

            // Before L1 prefiring reweighting
            h_bkg_pT_ele_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_ele_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_ele_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_mu_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_eleSS_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_eleSS_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_eleSS_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_L1Corr[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_L1Corr[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_ele_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_ele_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_mu_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_eleSS_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_eleSS_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_eleSS_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_L1Corr[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_L1Corr[iter]->SetLineColor(color);
            //
            h_bkg_pT_ele_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_ele_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_ele_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_mu_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_eleSS_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_eleSS_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_eleSS_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_L1Corr[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_L1Corr[iter]->SetDirectory(0);
            //
            s_pT_ele_before_L1Corr->Add(h_bkg_pT_ele_before_L1Corr[iter]);
            s_eta_ele_before_L1Corr->Add(h_bkg_eta_ele_before_L1Corr[iter]);
            s_phi_ele_before_L1Corr->Add(h_bkg_phi_ele_before_L1Corr[iter]);
            s_pT_mu_before_L1Corr->Add(h_bkg_pT_mu_before_L1Corr[iter]);
            s_eta_mu_before_L1Corr->Add(h_bkg_eta_mu_before_L1Corr[iter]);
            s_phi_mu_before_L1Corr->Add(h_bkg_phi_mu_before_L1Corr[iter]);
            s_pT_eleSS_before_L1Corr->Add(h_bkg_pT_eleSS_before_L1Corr[iter]);
            s_eta_eleSS_before_L1Corr->Add(h_bkg_eta_eleSS_before_L1Corr[iter]);
            s_phi_eleSS_before_L1Corr->Add(h_bkg_phi_eleSS_before_L1Corr[iter]);
            s_pT_muSS_before_L1Corr->Add(h_bkg_pT_muSS_before_L1Corr[iter]);
            s_eta_muSS_before_L1Corr->Add(h_bkg_eta_muSS_before_L1Corr[iter]);
            s_phi_muSS_before_L1Corr->Add(h_bkg_phi_muSS_before_L1Corr[iter]);

            // Before top pT reweighting
            h_bkg_pT_ele_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_ele_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_ele_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_mu_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_mu_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_mu_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_eleSS_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_eleSS_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_eleSS_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_pT_muSS_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_eta_muSS_before_TopPtCorr[iter]->SetFillColor(color);
            h_bkg_phi_muSS_before_TopPtCorr[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_ele_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_ele_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_mu_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_mu_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_mu_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_eleSS_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_eleSS_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_eleSS_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_pT_muSS_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_eta_muSS_before_TopPtCorr[iter]->SetLineColor(color);
            h_bkg_phi_muSS_before_TopPtCorr[iter]->SetLineColor(color);
            //
            h_bkg_pT_ele_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_ele_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_ele_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_mu_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_mu_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_mu_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_eleSS_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_eleSS_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_eleSS_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_pT_muSS_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_eta_muSS_before_TopPtCorr[iter]->SetDirectory(0);
            h_bkg_phi_muSS_before_TopPtCorr[iter]->SetDirectory(0);
            //
            s_pT_ele_before_TopPtCorr->Add(h_bkg_pT_ele_before_TopPtCorr[iter]);
            s_eta_ele_before_TopPtCorr->Add(h_bkg_eta_ele_before_TopPtCorr[iter]);
            s_phi_ele_before_TopPtCorr->Add(h_bkg_phi_ele_before_TopPtCorr[iter]);
            s_pT_mu_before_TopPtCorr->Add(h_bkg_pT_mu_before_TopPtCorr[iter]);
            s_eta_mu_before_TopPtCorr->Add(h_bkg_eta_mu_before_TopPtCorr[iter]);
            s_phi_mu_before_TopPtCorr->Add(h_bkg_phi_mu_before_TopPtCorr[iter]);
            s_pT_eleSS_before_TopPtCorr->Add(h_bkg_pT_eleSS_before_TopPtCorr[iter]);
            s_eta_eleSS_before_TopPtCorr->Add(h_bkg_eta_eleSS_before_TopPtCorr[iter]);
            s_phi_eleSS_before_TopPtCorr->Add(h_bkg_phi_eleSS_before_TopPtCorr[iter]);
            s_pT_muSS_before_TopPtCorr->Add(h_bkg_pT_muSS_before_TopPtCorr[iter]);
            s_eta_muSS_before_TopPtCorr->Add(h_bkg_eta_muSS_before_TopPtCorr[iter]);
            s_phi_muSS_before_TopPtCorr->Add(h_bkg_phi_muSS_before_TopPtCorr[iter]);

            // After all corrections
            h_bkg_pT_ele[iter]->SetFillColor(color);
            h_bkg_eta_ele[iter]->SetFillColor(color);
            h_bkg_phi_ele[iter]->SetFillColor(color);
            h_bkg_pT_mu[iter]->SetFillColor(color);
            h_bkg_eta_mu[iter]->SetFillColor(color);
            h_bkg_phi_mu[iter]->SetFillColor(color);
            h_bkg_pT_eleSS[iter]->SetFillColor(color);
            h_bkg_eta_eleSS[iter]->SetFillColor(color);
            h_bkg_phi_eleSS[iter]->SetFillColor(color);
            h_bkg_pT_muSS[iter]->SetFillColor(color);
            h_bkg_eta_muSS[iter]->SetFillColor(color);
            h_bkg_phi_muSS[iter]->SetFillColor(color);
            //
            h_bkg_pT_ele[iter]->SetLineColor(color);
            h_bkg_eta_ele[iter]->SetLineColor(color);
            h_bkg_phi_ele[iter]->SetLineColor(color);
            h_bkg_pT_mu[iter]->SetLineColor(color);
            h_bkg_eta_mu[iter]->SetLineColor(color);
            h_bkg_phi_mu[iter]->SetLineColor(color);
            h_bkg_pT_eleSS[iter]->SetLineColor(color);
            h_bkg_eta_eleSS[iter]->SetLineColor(color);
            h_bkg_phi_eleSS[iter]->SetLineColor(color);
            h_bkg_pT_muSS[iter]->SetLineColor(color);
            h_bkg_eta_muSS[iter]->SetLineColor(color);
            h_bkg_phi_muSS[iter]->SetLineColor(color);
            //
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
            //
            s_pT_ele->Add(h_bkg_pT_ele[iter]);
            s_eta_ele->Add(h_bkg_eta_ele[iter]);
            s_phi_ele->Add(h_bkg_phi_ele[iter]);
            s_pT_mu->Add(h_bkg_pT_mu[iter]);
            s_eta_mu->Add(h_bkg_eta_mu[iter]);
            s_phi_mu->Add(h_bkg_phi_mu[iter]);
            s_pT_eleSS->Add(h_bkg_pT_eleSS[iter]);
            s_eta_eleSS->Add(h_bkg_eta_eleSS[iter]);
            s_phi_eleSS->Add(h_bkg_phi_eleSS[iter]);
            s_pT_muSS->Add(h_bkg_pT_muSS[iter]);
            s_eta_muSS->Add(h_bkg_eta_muSS[iter]);
            s_phi_muSS->Add(h_bkg_phi_muSS[iter]);

            iter++;

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) // next - ttbar
                pr = _EMu_VVnST;
            if (pr == _EMu_DYTauTau_Full) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_pT_ele_before_RoccoR, *h_data_eta_ele_before_RoccoR, *h_data_phi_ele_before_RoccoR,
             *h_data_pT_mu_before_RoccoR, *h_data_eta_mu_before_RoccoR, *h_data_phi_mu_before_RoccoR,
             *h_data_pT_eleSS_before_RoccoR, *h_data_eta_eleSS_before_RoccoR, *h_data_phi_eleSS_before_RoccoR,
             *h_data_pT_muSS_before_RoccoR, *h_data_eta_muSS_before_RoccoR, *h_data_phi_muSS_before_RoccoR,
             *h_data_pT_ele, *h_data_eta_ele, *h_data_phi_ele, *h_data_pT_mu, *h_data_eta_mu, *h_data_phi_mu,
             *h_data_pT_eleSS, *h_data_eta_eleSS, *h_data_phi_eleSS, *h_data_pT_muSS, *h_data_eta_muSS, *h_data_phi_muSS;

        f_data->GetObject("h_mu_pT_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_mu_before_RoccoR);
        f_data->GetObject("h_mu_eta_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_mu_before_RoccoR);
        f_data->GetObject("h_mu_phi_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_mu_before_RoccoR);
        f_data->GetObject("h_muSS_pT_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_muSS_before_RoccoR);
        f_data->GetObject("h_muSS_eta_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_muSS_before_RoccoR);
        f_data->GetObject("h_muSS_phi_before_RocCorr_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_muSS_before_RoccoR);

        f_data->GetObject("h_ele_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_ele);
        f_data->GetObject("h_ele_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_ele);
        f_data->GetObject("h_ele_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_ele);
        f_data->GetObject("h_mu_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_mu);
        f_data->GetObject("h_mu_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_mu);
        f_data->GetObject("h_mu_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_mu);
        f_data->GetObject("h_eleSS_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_eleSS);
        f_data->GetObject("h_eleSS_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_eleSS);
        f_data->GetObject("h_eleSS_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_eleSS);
        f_data->GetObject("h_muSS_pT_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_pT_muSS);
        f_data->GetObject("h_muSS_eta_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_eta_muSS);
        f_data->GetObject("h_muSS_phi_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_phi_muSS);

        // Before Rochester correction
        h_data_pT_mu_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_mu_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_mu_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_pT_muSS_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_eta_muSS_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        h_data_phi_muSS_before_RoccoR->SetMarkerStyle(kFullDotLarge);
        //
        h_data_pT_mu_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_mu_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_mu_before_RoccoR->SetMarkerColor(kBlack);
        h_data_pT_muSS_before_RoccoR->SetMarkerColor(kBlack);
        h_data_eta_muSS_before_RoccoR->SetMarkerColor(kBlack);
        h_data_phi_muSS_before_RoccoR->SetMarkerColor(kBlack);
        //
        h_data_pT_mu_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_mu_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_mu_before_RoccoR->SetLineColor(kBlack);
        h_data_pT_muSS_before_RoccoR->SetLineColor(kBlack);
        h_data_eta_muSS_before_RoccoR->SetLineColor(kBlack);
        h_data_phi_muSS_before_RoccoR->SetLineColor(kBlack);
        //
        h_data_pT_mu_before_RoccoR->SetDirectory(0);
        h_data_eta_mu_before_RoccoR->SetDirectory(0);
        h_data_phi_mu_before_RoccoR->SetDirectory(0);
        h_data_pT_muSS_before_RoccoR->SetDirectory(0);
        h_data_eta_muSS_before_RoccoR->SetDirectory(0);
        h_data_phi_muSS_before_RoccoR->SetDirectory(0);

        // After all corrections
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
        //
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
        //
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
        //
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

        myRatioPlot_t *RP_pT_ele_before_PUCorr, *RP_eta_ele_before_PUCorr, *RP_phi_ele_before_PUCorr,
                      *RP_pT_mu_before_PUCorr, *RP_eta_mu_before_PUCorr, *RP_phi_mu_before_PUCorr,
                      *RP_pT_eleSS_before_PUCorr, *RP_eta_eleSS_before_PUCorr, *RP_phi_eleSS_before_PUCorr,
                      *RP_pT_muSS_before_PUCorr, *RP_eta_muSS_before_PUCorr, *RP_phi_muSS_before_PUCorr,
                      *RP_pT_ele_before_RoccoR, *RP_eta_ele_before_RoccoR, *RP_phi_ele_before_RoccoR,
                      *RP_pT_mu_before_RoccoR, *RP_eta_mu_before_RoccoR, *RP_phi_mu_before_RoccoR,
                      *RP_pT_eleSS_before_RoccoR, *RP_eta_eleSS_before_RoccoR, *RP_phi_eleSS_before_RoccoR,
                      *RP_pT_muSS_before_RoccoR, *RP_eta_muSS_before_RoccoR, *RP_phi_muSS_before_RoccoR,
                      *RP_pT_ele_before_EffCorr, *RP_eta_ele_before_EffCorr, *RP_phi_ele_before_EffCorr,
                      *RP_pT_mu_before_EffCorr, *RP_eta_mu_before_EffCorr, *RP_phi_mu_before_EffCorr,
                      *RP_pT_eleSS_before_EffCorr, *RP_eta_eleSS_before_EffCorr, *RP_phi_eleSS_before_EffCorr,
                      *RP_pT_muSS_before_EffCorr, *RP_eta_muSS_before_EffCorr, *RP_phi_muSS_before_EffCorr,
                      *RP_pT_ele_before_PVzCorr, *RP_eta_ele_before_PVzCorr, *RP_phi_ele_before_PVzCorr,
                      *RP_pT_mu_before_PVzCorr, *RP_eta_mu_before_PVzCorr, *RP_phi_mu_before_PVzCorr,
                      *RP_pT_eleSS_before_PVzCorr, *RP_eta_eleSS_before_PVzCorr, *RP_phi_eleSS_before_PVzCorr,
                      *RP_pT_muSS_before_PVzCorr, *RP_eta_muSS_before_PVzCorr, *RP_phi_muSS_before_PVzCorr,
                      *RP_pT_ele_before_L1Corr, *RP_eta_ele_before_L1Corr, *RP_phi_ele_before_L1Corr,
                      *RP_pT_mu_before_L1Corr, *RP_eta_mu_before_L1Corr, *RP_phi_mu_before_L1Corr,
                      *RP_pT_eleSS_before_L1Corr, *RP_eta_eleSS_before_L1Corr, *RP_phi_eleSS_before_L1Corr,
                      *RP_pT_muSS_before_L1Corr, *RP_eta_muSS_before_L1Corr, *RP_phi_muSS_before_L1Corr,
                      *RP_pT_ele_before_TopPtCorr, *RP_eta_ele_before_TopPtCorr, *RP_phi_ele_before_TopPtCorr,
                      *RP_pT_mu_before_TopPtCorr, *RP_eta_mu_before_TopPtCorr, *RP_phi_mu_before_TopPtCorr,
                      *RP_pT_eleSS_before_TopPtCorr, *RP_eta_eleSS_before_TopPtCorr, *RP_phi_eleSS_before_TopPtCorr,
                      *RP_pT_muSS_before_TopPtCorr, *RP_eta_muSS_before_TopPtCorr, *RP_phi_muSS_before_TopPtCorr,
                      *RP_pT_ele, *RP_eta_ele, *RP_phi_ele, *RP_pT_mu, *RP_eta_mu, *RP_phi_mu,
                      *RP_pT_eleSS, *RP_eta_eleSS, *RP_phi_eleSS, *RP_pT_muSS, *RP_eta_muSS, *RP_phi_muSS;

        RP_pT_ele_before_PUCorr = new myRatioPlot_t("RP_pT_ele_before_PUCorr", s_pT_ele_before_PUCorr, h_data_pT_ele);
        RP_eta_ele_before_PUCorr = new myRatioPlot_t("RP_eta_ele_before_PUCorr", s_eta_ele_before_PUCorr, h_data_eta_ele);
        RP_phi_ele_before_PUCorr = new myRatioPlot_t("RP_phi_ele_before_PUCorr", s_phi_ele_before_PUCorr, h_data_phi_ele);
        RP_pT_mu_before_PUCorr = new myRatioPlot_t("RP_pT_mu_before_PUCorr", s_pT_mu_before_PUCorr, h_data_pT_mu_before_RoccoR);
        RP_eta_mu_before_PUCorr = new myRatioPlot_t("RP_eta_mu_before_PUCorr", s_eta_mu_before_PUCorr, h_data_eta_mu_before_RoccoR);
        RP_phi_mu_before_PUCorr = new myRatioPlot_t("RP_phi_mu_before_PUCorr", s_phi_mu_before_PUCorr, h_data_phi_mu_before_RoccoR);
        RP_pT_eleSS_before_PUCorr = new myRatioPlot_t("RP_pT_eleSS_before_PUCorr", s_pT_eleSS_before_PUCorr, h_data_pT_eleSS);
        RP_eta_eleSS_before_PUCorr = new myRatioPlot_t("RP_eta_eleSS_before_PUCorr", s_eta_eleSS_before_PUCorr, h_data_eta_eleSS);
        RP_phi_eleSS_before_PUCorr = new myRatioPlot_t("RP_phi_eleSS_before_PUCorr", s_phi_eleSS_before_PUCorr, h_data_phi_eleSS);
        RP_pT_muSS_before_PUCorr = new myRatioPlot_t("RP_pT_muSS_before_PUCorr", s_pT_muSS_before_PUCorr, h_data_pT_muSS_before_RoccoR);
        RP_eta_muSS_before_PUCorr = new myRatioPlot_t("RP_eta_muSS_before_PUCorr", s_eta_muSS_before_PUCorr, h_data_eta_muSS_before_RoccoR);
        RP_phi_muSS_before_PUCorr = new myRatioPlot_t("RP_phi_muSS_before_PUCorr", s_phi_muSS_before_PUCorr, h_data_phi_muSS_before_RoccoR);

        RP_pT_mu_before_RoccoR = new myRatioPlot_t("RP_pT_mu_before_RoccoR", s_pT_mu_before_RoccoR, h_data_pT_mu_before_RoccoR);
        RP_eta_mu_before_RoccoR = new myRatioPlot_t("RP_eta_mu_before_RoccoR", s_eta_mu_before_RoccoR, h_data_eta_mu_before_RoccoR);
        RP_phi_mu_before_RoccoR = new myRatioPlot_t("RP_phi_mu_before_RoccoR", s_phi_mu_before_RoccoR, h_data_phi_mu_before_RoccoR);
        RP_pT_muSS_before_RoccoR = new myRatioPlot_t("RP_pT_muSS_before_RoccoR", s_pT_muSS_before_RoccoR, h_data_pT_muSS_before_RoccoR);
        RP_eta_muSS_before_RoccoR = new myRatioPlot_t("RP_eta_muSS_before_RoccoR", s_eta_muSS_before_RoccoR, h_data_eta_muSS_before_RoccoR);
        RP_phi_muSS_before_RoccoR = new myRatioPlot_t("RP_phi_muSS_before_RoccoR", s_phi_muSS_before_RoccoR, h_data_phi_muSS_before_RoccoR);

        RP_pT_ele_before_EffCorr = new myRatioPlot_t("RP_pT_ele_before_EffCorr", s_pT_ele_before_EffCorr, h_data_pT_ele);
        RP_eta_ele_before_EffCorr = new myRatioPlot_t("RP_eta_ele_before_EffCorr", s_eta_ele_before_EffCorr, h_data_eta_ele);
        RP_phi_ele_before_EffCorr = new myRatioPlot_t("RP_phi_ele_before_EffCorr", s_phi_ele_before_EffCorr, h_data_phi_ele);
        RP_pT_mu_before_EffCorr = new myRatioPlot_t("RP_pT_mu_before_EffCorr", s_pT_mu_before_EffCorr, h_data_pT_mu);
        RP_eta_mu_before_EffCorr = new myRatioPlot_t("RP_eta_mu_before_EffCorr", s_eta_mu_before_EffCorr, h_data_eta_mu);
        RP_phi_mu_before_EffCorr = new myRatioPlot_t("RP_phi_mu_before_EffCorr", s_phi_mu_before_EffCorr, h_data_phi_mu);
        RP_pT_eleSS_before_EffCorr = new myRatioPlot_t("RP_pT_eleSS_before_EffCorr", s_pT_eleSS_before_EffCorr, h_data_pT_eleSS);
        RP_eta_eleSS_before_EffCorr = new myRatioPlot_t("RP_eta_eleSS_before_EffCorr", s_eta_eleSS_before_EffCorr, h_data_eta_eleSS);
        RP_phi_eleSS_before_EffCorr = new myRatioPlot_t("RP_phi_eleSS_before_EffCorr", s_phi_eleSS_before_EffCorr, h_data_phi_eleSS);
        RP_pT_muSS_before_EffCorr = new myRatioPlot_t("RP_pT_muSS_before_EffCorr", s_pT_muSS_before_EffCorr, h_data_pT_muSS);
        RP_eta_muSS_before_EffCorr = new myRatioPlot_t("RP_eta_muSS_before_EffCorr", s_eta_muSS_before_EffCorr, h_data_eta_muSS);
        RP_phi_muSS_before_EffCorr = new myRatioPlot_t("RP_phi_muSS_before_EffCorr", s_phi_muSS_before_EffCorr, h_data_phi_muSS);

        RP_pT_ele_before_PVzCorr = new myRatioPlot_t("RP_pT_ele_before_PVzCorr", s_pT_ele_before_PVzCorr, h_data_pT_ele);
        RP_eta_ele_before_PVzCorr = new myRatioPlot_t("RP_eta_ele_before_PVzCorr", s_eta_ele_before_PVzCorr, h_data_eta_ele);
        RP_phi_ele_before_PVzCorr = new myRatioPlot_t("RP_phi_ele_before_PVzCorr", s_phi_ele_before_PVzCorr, h_data_phi_ele);
        RP_pT_mu_before_PVzCorr = new myRatioPlot_t("RP_pT_mu_before_PVzCorr", s_pT_mu_before_PVzCorr, h_data_pT_mu);
        RP_eta_mu_before_PVzCorr = new myRatioPlot_t("RP_eta_mu_before_PVzCorr", s_eta_mu_before_PVzCorr, h_data_eta_mu);
        RP_phi_mu_before_PVzCorr = new myRatioPlot_t("RP_phi_mu_before_PVzCorr", s_phi_mu_before_PVzCorr, h_data_phi_mu);
        RP_pT_eleSS_before_PVzCorr = new myRatioPlot_t("RP_pT_eleSS_before_PVzCorr", s_pT_eleSS_before_PVzCorr, h_data_pT_eleSS);
        RP_eta_eleSS_before_PVzCorr = new myRatioPlot_t("RP_eta_eleSS_before_PVzCorr", s_eta_eleSS_before_PVzCorr, h_data_eta_eleSS);
        RP_phi_eleSS_before_PVzCorr = new myRatioPlot_t("RP_phi_eleSS_before_PVzCorr", s_phi_eleSS_before_PVzCorr, h_data_phi_eleSS);
        RP_pT_muSS_before_PVzCorr = new myRatioPlot_t("RP_pT_muSS_before_PVzCorr", s_pT_muSS_before_PVzCorr, h_data_pT_muSS);
        RP_eta_muSS_before_PVzCorr = new myRatioPlot_t("RP_eta_muSS_before_PVzCorr", s_eta_muSS_before_PVzCorr, h_data_eta_muSS);
        RP_phi_muSS_before_PVzCorr = new myRatioPlot_t("RP_phi_muSS_before_PVzCorr", s_phi_muSS_before_PVzCorr, h_data_phi_muSS);

        RP_pT_ele_before_L1Corr = new myRatioPlot_t("RP_pT_ele_before_L1Corr", s_pT_ele_before_L1Corr, h_data_pT_ele);
        RP_eta_ele_before_L1Corr = new myRatioPlot_t("RP_eta_ele_before_L1Corr", s_eta_ele_before_L1Corr, h_data_eta_ele);
        RP_phi_ele_before_L1Corr = new myRatioPlot_t("RP_phi_ele_before_L1Corr", s_phi_ele_before_L1Corr, h_data_phi_ele);
        RP_pT_mu_before_L1Corr = new myRatioPlot_t("RP_pT_mu_before_L1Corr", s_pT_mu_before_L1Corr, h_data_pT_mu);
        RP_eta_mu_before_L1Corr = new myRatioPlot_t("RP_eta_mu_before_L1Corr", s_eta_mu_before_L1Corr, h_data_eta_mu);
        RP_phi_mu_before_L1Corr = new myRatioPlot_t("RP_phi_mu_before_L1Corr", s_phi_mu_before_L1Corr, h_data_phi_mu);
        RP_pT_eleSS_before_L1Corr = new myRatioPlot_t("RP_pT_eleSS_before_L1Corr", s_pT_eleSS_before_L1Corr, h_data_pT_eleSS);
        RP_eta_eleSS_before_L1Corr = new myRatioPlot_t("RP_eta_eleSS_before_L1Corr", s_eta_eleSS_before_L1Corr, h_data_eta_eleSS);
        RP_phi_eleSS_before_L1Corr = new myRatioPlot_t("RP_phi_eleSS_before_L1Corr", s_phi_eleSS_before_L1Corr, h_data_phi_eleSS);
        RP_pT_muSS_before_L1Corr = new myRatioPlot_t("RP_pT_muSS_before_L1Corr", s_pT_muSS_before_L1Corr, h_data_pT_muSS);
        RP_eta_muSS_before_L1Corr = new myRatioPlot_t("RP_eta_muSS_before_L1Corr", s_eta_muSS_before_L1Corr, h_data_eta_muSS);
        RP_phi_muSS_before_L1Corr = new myRatioPlot_t("RP_phi_muSS_before_L1Corr", s_phi_muSS_before_L1Corr, h_data_phi_muSS);

        RP_pT_ele_before_TopPtCorr = new myRatioPlot_t("RP_pT_ele_before_TopPtCorr", s_pT_ele_before_TopPtCorr, h_data_pT_ele);
        RP_eta_ele_before_TopPtCorr = new myRatioPlot_t("RP_eta_ele_before_TopPtCorr", s_eta_ele_before_TopPtCorr, h_data_eta_ele);
        RP_phi_ele_before_TopPtCorr = new myRatioPlot_t("RP_phi_ele_before_TopPtCorr", s_phi_ele_before_TopPtCorr, h_data_phi_ele);
        RP_pT_mu_before_TopPtCorr = new myRatioPlot_t("RP_pT_mu_before_TopPtCorr", s_pT_mu_before_TopPtCorr, h_data_pT_mu);
        RP_eta_mu_before_TopPtCorr = new myRatioPlot_t("RP_eta_mu_before_TopPtCorr", s_eta_mu_before_TopPtCorr, h_data_eta_mu);
        RP_phi_mu_before_TopPtCorr = new myRatioPlot_t("RP_phi_mu_before_TopPtCorr", s_phi_mu_before_TopPtCorr, h_data_phi_mu);
        RP_pT_eleSS_before_TopPtCorr = new myRatioPlot_t("RP_pT_eleSS_before_TopPtCorr", s_pT_eleSS_before_TopPtCorr, h_data_pT_eleSS);
        RP_eta_eleSS_before_TopPtCorr = new myRatioPlot_t("RP_eta_eleSS_before_TopPtCorr", s_eta_eleSS_before_TopPtCorr, h_data_eta_eleSS);
        RP_phi_eleSS_before_TopPtCorr = new myRatioPlot_t("RP_phi_eleSS_before_TopPtCorr", s_phi_eleSS_before_TopPtCorr, h_data_phi_eleSS);
        RP_pT_muSS_before_TopPtCorr = new myRatioPlot_t("RP_pT_muSS_before_TopPtCorr", s_pT_muSS_before_TopPtCorr, h_data_pT_muSS);
        RP_eta_muSS_before_TopPtCorr = new myRatioPlot_t("RP_eta_muSS_before_TopPtCorr", s_eta_muSS_before_TopPtCorr, h_data_eta_muSS);
        RP_phi_muSS_before_TopPtCorr = new myRatioPlot_t("RP_phi_muSS_before_TopPtCorr", s_phi_muSS_before_TopPtCorr, h_data_phi_muSS);

        RP_pT_ele = new myRatioPlot_t("RP_pT_ele", s_pT_ele, h_data_pT_ele);
        RP_eta_ele = new myRatioPlot_t("RP_eta_ele", s_eta_ele, h_data_eta_ele);
        RP_phi_ele = new myRatioPlot_t("RP_phi_ele", s_phi_ele, h_data_phi_ele);
        RP_pT_mu = new myRatioPlot_t("RP_pT_mu", s_pT_mu, h_data_pT_mu);
        RP_eta_mu = new myRatioPlot_t("RP_eta_mu", s_eta_mu, h_data_eta_mu);
        RP_phi_mu = new myRatioPlot_t("RP_phi_mu", s_phi_mu, h_data_phi_mu);
        RP_pT_eleSS = new myRatioPlot_t("RP_pT_eleSS", s_pT_eleSS, h_data_pT_eleSS);
        RP_eta_eleSS = new myRatioPlot_t("RP_eta_eleSS", s_eta_eleSS, h_data_eta_eleSS);
        RP_phi_eleSS = new myRatioPlot_t("RP_phi_eleSS", s_phi_eleSS, h_data_phi_eleSS);
        RP_pT_muSS = new myRatioPlot_t("RP_pT_muSS", s_pT_muSS, h_data_pT_muSS);
        RP_eta_muSS = new myRatioPlot_t("RP_eta_muSS", s_eta_muSS, h_data_eta_muSS);
        RP_phi_muSS = new myRatioPlot_t("RP_phi_muSS", s_phi_muSS, h_data_phi_muSS);

        RP_pT_ele_before_PUCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c]  before PU correction", 0, 600);
        RP_eta_ele_before_PUCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}}  before PU correction", -4, 4);
        RP_phi_ele_before_PUCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}}  before PU correction", -4, 4);
        RP_pT_mu_before_PUCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c]  before PU correction", 0, 600);
        RP_eta_mu_before_PUCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}}  before PU correction", -4, 4);
        RP_phi_mu_before_PUCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}}  before PU correction", -4, 4);
        RP_pT_eleSS_before_PUCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c  before PU correction]", 0, 600);
        RP_eta_eleSS_before_PUCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign)  before PU correction", -4, 4);
        RP_phi_eleSS_before_PUCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign)  before PU correction", -4, 4);
        RP_pT_muSS_before_PUCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c]  before PU correction", 0, 600);
        RP_eta_muSS_before_PUCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign)  before PU correction", -4, 4);
        RP_phi_muSS_before_PUCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign)  before PU correction", -4, 4);

        RP_pT_mu_before_RoccoR->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c] before Rochester correction", 0, 600);
        RP_eta_mu_before_RoccoR->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} before Rochester correction", -4, 4);
        RP_phi_mu_before_RoccoR->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} before Rochester correction", -4, 4);
        RP_pT_muSS_before_RoccoR->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c] before Rochester correction", 0, 600);
        RP_eta_muSS_before_RoccoR->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before Rochester correction", -4, 4);
        RP_phi_muSS_before_RoccoR->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before Rochester correction", -4, 4);

        RP_pT_ele_before_EffCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c] before Efficiency SF", 0, 600);
        RP_eta_ele_before_EffCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} before Efficiency SF", -4, 4);
        RP_phi_ele_before_EffCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} before Efficiency SF", -4, 4);
        RP_pT_mu_before_EffCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c] before Efficiency SF", 0, 600);
        RP_eta_mu_before_EffCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} before Efficiency SF", -4, 4);
        RP_phi_mu_before_EffCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} before Efficiency SF", -4, 4);
        RP_pT_eleSS_before_EffCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c] before Efficiency SF", 0, 600);
        RP_eta_eleSS_before_EffCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign) before Efficiency SF", -4, 4);
        RP_phi_eleSS_before_EffCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign) before Efficiency SF", -4, 4);
        RP_pT_muSS_before_EffCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c] before Efficiency SF", 0, 600);
        RP_eta_muSS_before_EffCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before Efficiency SF", -4, 4);
        RP_phi_muSS_before_EffCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before Efficiency SF", -4, 4);

        RP_pT_ele_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c] before PVz correction", 0, 600);
        RP_eta_ele_before_PVzCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} before PVz correction", -4, 4);
        RP_phi_ele_before_PVzCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} before PVz correction", -4, 4);
        RP_pT_mu_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c] before PVz correction", 0, 600);
        RP_eta_mu_before_PVzCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} before PVz correction", -4, 4);
        RP_phi_mu_before_PVzCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} before PVz correction", -4, 4);
        RP_pT_eleSS_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c] before PVz correction", 0, 600);
        RP_eta_eleSS_before_PVzCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign) before PVz correction", -4, 4);
        RP_phi_eleSS_before_PVzCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign) before PVz correction", -4, 4);
        RP_pT_muSS_before_PVzCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c] before PVz correction", 0, 600);
        RP_eta_muSS_before_PVzCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before PVz correction", -4, 4);
        RP_phi_muSS_before_PVzCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before PVz correction", -4, 4);

        RP_pT_ele_before_L1Corr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c] before L1 prefiring correction", 0, 600);
        RP_eta_ele_before_L1Corr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} before L1 prefiring correction", -4, 4);
        RP_phi_ele_before_L1Corr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} before L1 prefiring correction", -4, 4);
        RP_pT_mu_before_L1Corr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c] before L1 prefiring correction", 0, 600);
        RP_eta_mu_before_L1Corr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} before L1 prefiring correction", -4, 4);
        RP_phi_mu_before_L1Corr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} before L1 prefiring correction", -4, 4);
        RP_pT_eleSS_before_L1Corr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c] before L1 prefiring correction", 0, 600);
        RP_eta_eleSS_before_L1Corr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign) before L1 prefiring correction", -4, 4);
        RP_phi_eleSS_before_L1Corr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign) before L1 prefiring correction", -4, 4);
        RP_pT_muSS_before_L1Corr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c] before L1 prefiring correction", 0, 600);
        RP_eta_muSS_before_L1Corr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before L1 prefiring correction", -4, 4);
        RP_phi_muSS_before_L1Corr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before L1 prefiring correction", -4, 4);

        RP_pT_ele_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 600);
        RP_eta_ele_before_TopPtCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_phi_ele_before_TopPtCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_pT_mu_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 600);
        RP_eta_mu_before_TopPtCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_phi_mu_before_TopPtCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_pT_eleSS_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 600);
        RP_eta_eleSS_before_TopPtCorr->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_phi_eleSS_before_TopPtCorr->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_pT_muSS_before_TopPtCorr->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c] before top p_{#lower[-0.2]{T}} reweighting", 0, 600);
        RP_eta_muSS_before_TopPtCorr->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);
        RP_phi_muSS_before_TopPtCorr->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign) before top p_{#lower[-0.2]{T}} reweighting", -4, 4);

        RP_pT_ele->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} [GeV/c]", 0, 600);
        RP_eta_ele->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}}", -4, 4);
        RP_phi_ele->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}}", -4, 4);
        RP_pT_mu->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} [GeV/c]", 0, 600);
        RP_eta_mu->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}}", -4, 4);
        RP_phi_mu->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}}", -4, 4);
        RP_pT_eleSS->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.3]{#font[12]{e}}}} (same-sign) [GeV/c]", 0, 600);
        RP_eta_eleSS->SetPlots("#eta_{#font[12]{#lower[-0.4]{#scale[1.2]{e}}}} (same-sign)", -4, 4);
        RP_phi_eleSS->SetPlots("#phi_{#font[12]{#lower[-0.3]{#scale[1.2]{e}}}} (same-sign)", -4, 4);
        RP_pT_muSS->SetPlots("p_{#lower[-0.2]{T} #lower[-0.3]{#scale[1.15]{#mu}}} (same-sign) [GeV/c]", 0, 600);
        RP_eta_muSS->SetPlots("#eta_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign)", -4, 4);
        RP_phi_muSS->SetPlots("#phi_{#lower[-0.25]{#scale[1.2]{#mu}}} (same-sign)", -4, 4);

        TLegend *legend = new TLegend(0.62, 0.77, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_pT_ele, "Data", "lp");
        legend->AddEntry(h_bkg_pT_ele[6+isWJ], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_pT_ele[5+isWJ], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_pT_ele[4+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_pT_ele[3+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_pT_ele[2+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_pT_ele[1+isWJ], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_pT_ele[0+isWJ], "#font[12]{#scale[1.1]{WW}}", "f");
        if (isWJ)
        {
            legend->AddEntry(h_bkg_pT_ele[0], "#font[12]{#scale[1.1]{W}}+Jets", "f");
        }

        RP_pT_ele_before_PUCorr->ImportLegend(legend);
        RP_eta_ele_before_PUCorr->ImportLegend(legend);
        RP_phi_ele_before_PUCorr->ImportLegend(legend);
        RP_pT_mu_before_PUCorr->ImportLegend(legend);
        RP_eta_mu_before_PUCorr->ImportLegend(legend);
        RP_phi_mu_before_PUCorr->ImportLegend(legend);
        RP_pT_eleSS_before_PUCorr->ImportLegend(legend);
        RP_eta_eleSS_before_PUCorr->ImportLegend(legend);
        RP_phi_eleSS_before_PUCorr->ImportLegend(legend);
        RP_pT_muSS_before_PUCorr->ImportLegend(legend);
        RP_eta_muSS_before_PUCorr->ImportLegend(legend);
        RP_phi_muSS_before_PUCorr->ImportLegend(legend);

        RP_pT_mu_before_RoccoR->ImportLegend(legend);
        RP_eta_mu_before_RoccoR->ImportLegend(legend);
        RP_phi_mu_before_RoccoR->ImportLegend(legend);
        RP_pT_muSS_before_RoccoR->ImportLegend(legend);
        RP_eta_muSS_before_RoccoR->ImportLegend(legend);
        RP_phi_muSS_before_RoccoR->ImportLegend(legend);

        RP_pT_ele_before_EffCorr->ImportLegend(legend);
        RP_eta_ele_before_EffCorr->ImportLegend(legend);
        RP_phi_ele_before_EffCorr->ImportLegend(legend);
        RP_pT_mu_before_EffCorr->ImportLegend(legend);
        RP_eta_mu_before_EffCorr->ImportLegend(legend);
        RP_phi_mu_before_EffCorr->ImportLegend(legend);
        RP_pT_eleSS_before_EffCorr->ImportLegend(legend);
        RP_eta_eleSS_before_EffCorr->ImportLegend(legend);
        RP_phi_eleSS_before_EffCorr->ImportLegend(legend);
        RP_pT_muSS_before_EffCorr->ImportLegend(legend);
        RP_eta_muSS_before_EffCorr->ImportLegend(legend);
        RP_phi_muSS_before_EffCorr->ImportLegend(legend);

        RP_pT_ele_before_PVzCorr->ImportLegend(legend);
        RP_eta_ele_before_PVzCorr->ImportLegend(legend);
        RP_phi_ele_before_PVzCorr->ImportLegend(legend);
        RP_pT_mu_before_PVzCorr->ImportLegend(legend);
        RP_eta_mu_before_PVzCorr->ImportLegend(legend);
        RP_phi_mu_before_PVzCorr->ImportLegend(legend);
        RP_pT_eleSS_before_PVzCorr->ImportLegend(legend);
        RP_eta_eleSS_before_PVzCorr->ImportLegend(legend);
        RP_phi_eleSS_before_PVzCorr->ImportLegend(legend);
        RP_pT_muSS_before_PVzCorr->ImportLegend(legend);
        RP_eta_muSS_before_PVzCorr->ImportLegend(legend);
        RP_phi_muSS_before_PVzCorr->ImportLegend(legend);

        RP_pT_ele_before_L1Corr->ImportLegend(legend);
        RP_eta_ele_before_L1Corr->ImportLegend(legend);
        RP_phi_ele_before_L1Corr->ImportLegend(legend);
        RP_pT_mu_before_L1Corr->ImportLegend(legend);
        RP_eta_mu_before_L1Corr->ImportLegend(legend);
        RP_phi_mu_before_L1Corr->ImportLegend(legend);
        RP_pT_eleSS_before_L1Corr->ImportLegend(legend);
        RP_eta_eleSS_before_L1Corr->ImportLegend(legend);
        RP_phi_eleSS_before_L1Corr->ImportLegend(legend);
        RP_pT_muSS_before_L1Corr->ImportLegend(legend);
        RP_eta_muSS_before_L1Corr->ImportLegend(legend);
        RP_phi_muSS_before_L1Corr->ImportLegend(legend);

        RP_pT_ele_before_TopPtCorr->ImportLegend(legend);
        RP_eta_ele_before_TopPtCorr->ImportLegend(legend);
        RP_phi_ele_before_TopPtCorr->ImportLegend(legend);
        RP_pT_mu_before_TopPtCorr->ImportLegend(legend);
        RP_eta_mu_before_TopPtCorr->ImportLegend(legend);
        RP_phi_mu_before_TopPtCorr->ImportLegend(legend);
        RP_pT_eleSS_before_TopPtCorr->ImportLegend(legend);
        RP_eta_eleSS_before_TopPtCorr->ImportLegend(legend);
        RP_phi_eleSS_before_TopPtCorr->ImportLegend(legend);
        RP_pT_muSS_before_TopPtCorr->ImportLegend(legend);
        RP_eta_muSS_before_TopPtCorr->ImportLegend(legend);
        RP_phi_muSS_before_TopPtCorr->ImportLegend(legend);

        RP_pT_ele->ImportLegend(legend);
        RP_eta_ele->ImportLegend(legend);
        RP_phi_ele->ImportLegend(legend);
        RP_pT_mu->ImportLegend(legend);
        RP_eta_mu->ImportLegend(legend);
        RP_phi_mu->ImportLegend(legend);
        RP_pT_eleSS->ImportLegend(legend);
        RP_eta_eleSS->ImportLegend(legend);
        RP_phi_eleSS->ImportLegend(legend);
        RP_pT_muSS->ImportLegend(legend);
        RP_eta_muSS->ImportLegend(legend);
        RP_phi_muSS->ImportLegend(legend);

        RP_pT_ele_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_eta_ele_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_phi_ele_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_pT_mu_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_pT_eleSS_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_eta_eleSS_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_phi_eleSS_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_PUCorr->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_PUCorr->Draw(0.5, 1e6, 0);

        RP_pT_mu_before_RoccoR->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_RoccoR->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_RoccoR->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_RoccoR->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_RoccoR->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_RoccoR->Draw(0.5, 1e6, 0);

        RP_pT_ele_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_eta_ele_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_phi_ele_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_pT_mu_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_pT_eleSS_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_eta_eleSS_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_phi_eleSS_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_EffCorr->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_EffCorr->Draw(0.5, 1e6, 0);

        RP_pT_ele_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_eta_ele_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_phi_ele_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_pT_mu_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_pT_eleSS_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_eta_eleSS_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_phi_eleSS_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_PVzCorr->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_PVzCorr->Draw(0.5, 1e6, 0);

        RP_pT_ele_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_eta_ele_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_phi_ele_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_pT_mu_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_pT_eleSS_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_eta_eleSS_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_phi_eleSS_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_L1Corr->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_L1Corr->Draw(0.5, 1e6, 0);

        RP_pT_ele_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_eta_ele_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_phi_ele_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_pT_mu_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_eta_mu_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_phi_mu_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_pT_eleSS_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_eta_eleSS_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_phi_eleSS_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_pT_muSS_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_eta_muSS_before_TopPtCorr->Draw(0.5, 1e6, 0);
        RP_phi_muSS_before_TopPtCorr->Draw(0.5, 1e6, 0);

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

    if(whichGraphs=="ALL" || whichGraphs=="nVTX" || whichGraphs=="NPV" || whichGraphs.Contains("PILEUP"))
    {
        count_drawn++;

        THStack *s_nVTX_before_PUCorr, *s_nVTX_before_EffCorr, *s_nVTX;
        s_nVTX_before_PUCorr = new THStack("s_nVTX_before_PUCorr", "");
        s_nVTX_before_EffCorr = new THStack("s_nVTX_before_EffCorr", "");
        s_nVTX = new THStack("s_nVTX", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_nVTX_before_PUCorr[8], *h_bkg_nVTX_before_EffCorr[8], *h_bkg_nVTX[8];
        Int_t iter = 0;

        for (SelProc_t pr = _EMu_WJets_Full; pr > _EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (!isWJ && pr == _EMu_WJets) continue;
            if (pr == _EndOf_EMu_VVnST_Normal) continue;

            f_bkg->GetObject("h_nVTX_before_PUCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_PUCorr[iter]);
            f_bkg->GetObject("h_nVTX_before_EffCorr_"+Mgr.Procname[pr], h_bkg_nVTX_before_EffCorr[iter]);
            f_bkg->GetObject("h_nVTX_"+Mgr.Procname[pr], h_bkg_nVTX[iter]);
            removeNegativeBins(h_bkg_nVTX_before_PUCorr[iter]);
            removeNegativeBins(h_bkg_nVTX_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_nVTX[iter]);


            Color_t color = kBlack;
            if (pr == _EMu_WJets_Full) color = kRed - 2;
            if (pr == _EMu_WW) color = kMagenta - 5;
            if (pr == _EMu_WZ) color = kMagenta - 2;
            if (pr == _EMu_ZZ) color = kMagenta - 6;
            if (pr == _EMu_tbarW) color = kGreen - 2;
            if (pr == _EMu_tW) color = kGreen + 2;
            if (pr == _EMu_ttbar_Full) color = kCyan + 2;
            if (pr == _EMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_nVTX_before_PUCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_PUCorr[iter]->SetDirectory(0);
            s_nVTX_before_PUCorr->Add(h_bkg_nVTX_before_PUCorr[iter]);

            h_bkg_nVTX_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_nVTX_before_EffCorr[iter]->SetDirectory(0);
            s_nVTX_before_EffCorr->Add(h_bkg_nVTX_before_EffCorr[iter]);

            h_bkg_nVTX[iter]->SetFillColor(color);
            h_bkg_nVTX[iter]->SetLineColor(color);
            h_bkg_nVTX[iter]->SetDirectory(0);
            s_nVTX->Add(h_bkg_nVTX[iter]);

            iter++;

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) // next - ttbar
                pr = _EMu_VVnST;
            if (pr == _EMu_DYTauTau_Full) // last process
                break;

        } // End of for(bkg)

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_nVTX;

        f_data->GetObject("h_nVTX_"+Mgr.Procname[_EMu_SingleMuon_Full], h_data_nVTX);

        h_data_nVTX->SetMarkerStyle(kFullDotLarge);
        h_data_nVTX->SetMarkerColor(kBlack);
        h_data_nVTX->SetLineColor(kBlack);

        h_data_nVTX->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_nVTX_before_PUCorr, *RP_nVTX_before_EffCorr, *RP_nVTX;
        RP_nVTX_before_PUCorr = new myRatioPlot_t("RP_nVTX_before_PUCorr", s_nVTX_before_PUCorr, h_data_nVTX);
        RP_nVTX_before_EffCorr = new myRatioPlot_t("RP_nVTX_before_EffCorr", s_nVTX_before_EffCorr, h_data_nVTX);
        RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_data_nVTX);

        RP_nVTX_before_PUCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.2]{#font[12]{e}#mu}}} before PU correction", 0, 50);
        RP_nVTX_before_EffCorr->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.2]{#font[12]{e}#mu}}} before Efficiency correction", 0, 50);
        RP_nVTX->SetPlots("N_{#lower[-0.2]{VTX} #lower[-0.25]{#scale[1.2]{#font[12]{e}#mu}}}", 0, 50);

        TLegend *legend = new TLegend(0.6, 0.75, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(h_data_nVTX, "Data", "lp");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[6+isWJ], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[5+isWJ], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[4+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[3+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[2+isWJ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[1+isWJ], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_nVTX_before_PUCorr[0+isWJ], "#font[12]{#scale[1.1]{WW}}", "f");
        if (isWJ)
        {
            legend->AddEntry(h_bkg_nVTX_before_PUCorr[0], "#font[12]{#scale[1.1]{W}}+Jets", "f");
        }

        RP_nVTX_before_PUCorr->ImportLegend(legend);
        RP_nVTX_before_EffCorr->ImportLegend(legend);
        RP_nVTX->ImportLegend(legend);

        RP_nVTX_before_PUCorr->Draw(0.5, 3e6, 0);
        RP_nVTX_before_EffCorr->Draw(0.5, 3e6, 0);
        RP_nVTX->Draw(0.5, 3e6, 0);

        cout << "nVTX Chi^2 before PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX_before_PUCorr) << endl;
        cout << "nVTX Chi^2 after PU reweight: " << CompChiSquared(h_data_nVTX, s_nVTX) << endl;

    } // End of if(nVTX)
    */

} // End of EMu_HistDrawer()


/// ############################################################################### ///
/// ----------------------------- BKG ESTIMATIONS --------------------------------- ///
/// ############################################################################### ///
void Est_HistDrawer(Int_t FR_systErr)
{
    LocalFileMgr Mgr;

//############################ ELECTRON CHANNEL #####################################

    Mgr.SetProc(_EE_DY_Full);
    TString name_DY_ee = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root";
    TFile* f_DY_ee = new TFile(name_DY_ee, "READ");
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY_ee->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_EE_Bkg_Full);
    TString name_bkg_ee = Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root";
    TFile* f_bkg_ee = new TFile(name_bkg_ee, "READ");
    if (f_bkg_ee->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root" << " opened successfully" << endl;
    TString name_bkg_est_ee = Mgr.HistLocation+"EstBkg_EE.root";
    TFile* f_bkg_est_ee = new TFile(name_bkg_est_ee, "READ");
    if (f_bkg_est_ee->IsOpen()) std::cout << "File " << "EstBkg_EE.root" << " opened successfully" << endl;
    TString name_WJets_est_ee = Mgr.HistLocation+"EstWJets_EE.root";
    TFile* f_WJets_est_ee = new TFile(name_WJets_est_ee, "READ");
    if (f_WJets_est_ee->IsOpen()) std::cout << "File " << "EstWJets_EE.root" << " opened successfully" << endl;
    TString name_QCD_est_ee = Mgr.HistLocation+"EstQCD_EE.root";
    TFile* f_QCD_est_ee = new TFile(name_QCD_est_ee, "READ");
    if (f_QCD_est_ee->IsOpen()) std::cout << "File " << "EstQCD_EE.root" << " opened successfully" << endl;
    Mgr.SetProc(_EE_DoubleEG_Full);
    TString name_data_ee = Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root";
    TFile* f_data_ee = new TFile(name_data_ee, "READ");
    if (f_data_ee->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root" << " opened successfully\n" << endl;

    THStack *s_mass_ee = new THStack("s_mass_ee", "");
    THStack *s_mass_ee_wFR = new THStack("s_mass_ee_wFR", "");
    THStack *s_mass_ee2 = new THStack("s_mass_ee2", "");

//--------------------------------- FakeEle bkg -----------------------------------------------------
    TH1D *h_fakes_mass_ee[2];
    f_QCD_est_ee->GetObject("h_QCD_est_fit", h_fakes_mass_ee[0]);
    removeNegativeBins(h_fakes_mass_ee[0]);
    h_fakes_mass_ee[0]->SetFillColor(kRed + 3);
    h_fakes_mass_ee[0]->SetLineColor(kRed + 3);
    h_fakes_mass_ee[0]->SetDirectory(0);
    s_mass_ee_wFR->Add(h_fakes_mass_ee[0]);

    f_WJets_est_ee->GetObject("h_WJET_est_fit", h_fakes_mass_ee[1]);
    removeNegativeBins(h_fakes_mass_ee[1]);
    h_fakes_mass_ee[1]->SetFillColor(kRed - 2);
    h_fakes_mass_ee[1]->SetLineColor(kRed - 2);
    h_fakes_mass_ee[1]->SetDirectory(0);
    h_fakes_mass_ee[1]->GetYaxis()->SetRangeUser(0.01,1e7);
    s_mass_ee_wFR->Add(h_fakes_mass_ee[1]);

//----------------------------------- MC bkg -------------------------------------------------------
    TH1D *h_bkg_mass_ee[9], *h_bkg_mass_ee2[9];
    Int_t iter = 0;

    for (SelProc_t pr = _EE_QCDEMEnriched_Full; pr > _EndOf_EE_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
//        if (pr == _EE_QCDEMEnriched_Full)
//        {
//            iter++;
//            continue;
//        }

//        if (pr == _EE_WJets_Full || pr == _EE_QCDEMEnriched_Full)
        if (pr == _EE_WJets_Full || pr == _EE_QCDEMEnriched_Full || pr == _EE_WZ || pr == _EE_ZZ)
        {
            f_bkg_ee->GetObject("h_mass_"+Mgr.Procname[pr], h_bkg_mass_ee[iter]);
            f_bkg_ee->GetObject("h_mass2_"+Mgr.Procname[pr], h_bkg_mass_ee2[iter]);
        }
        else
        {
            f_bkg_est_ee->GetObject("h_ee_mass_Est_"+Mgr.Procname[pr], h_bkg_mass_ee[iter]);
            f_bkg_est_ee->GetObject("h_ee_mass_Est2_"+Mgr.Procname[pr], h_bkg_mass_ee2[iter]);
        }
        removeNegativeBins(h_bkg_mass_ee[iter]);
        removeNegativeBins(h_bkg_mass_ee2[iter]);

        Color_t color = kBlack;
        if (pr == _EE_QCDEMEnriched_Full) color = kRed + 3;
        if (pr == _EE_WJets_Full) color = kRed - 2;
        if (pr == _EE_WW) color = kMagenta - 5;
        if (pr == _EE_WZ) color = kMagenta - 2;
        if (pr == _EE_ZZ) color = kMagenta - 6;
        if (pr == _EE_tbarW) color = kGreen - 2;
        if (pr == _EE_tW) color = kGreen + 2;
        if (pr == _EE_ttbar_Full) color = kCyan + 2;
        if (pr == _EE_DYTauTau_Full) color = kOrange - 5;

        h_bkg_mass_ee[iter]->SetFillColor(color);
        h_bkg_mass_ee[iter]->SetLineColor(color);
        h_bkg_mass_ee[iter]->SetDirectory(0);
        s_mass_ee->Add(h_bkg_mass_ee[iter]);
        h_bkg_mass_ee2[iter]->SetFillColor(color);
        h_bkg_mass_ee2[iter]->SetLineColor(color);
        h_bkg_mass_ee2[iter]->SetDirectory(0);
        s_mass_ee2->Add(h_bkg_mass_ee2[iter]);

        if ((pr < _EndOf_EE_VVnST_Normal || pr > _EndOf_EE_MCbkg_Normal) && pr != _EE_QCDEMEnriched_Full && pr != _EE_WJets_Full)
            s_mass_ee_wFR->Add(h_bkg_mass_ee[iter]);

        iter++;

        if (pr == _EE_WJets_Full) // next - WW
            pr = _EndOf_EE_VVnST_Normal;
        if (pr == _EE_tW) // next -- ttbar
            pr = _EE_VVnST;
        if (pr == _EE_DYTauTau_Full) // last
            break;

    } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

    TH1D *h_DY_mass_ee, *h_DY_mass_ee2;
    f_DY_ee->GetObject("h_mass_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_ee);
    f_DY_ee->GetObject("h_mass2_"+Mgr.Procname[_EE_DY_Full], h_DY_mass_ee2);
    removeNegativeBins(h_DY_mass_ee);
    removeNegativeBins(h_DY_mass_ee2);
    h_DY_mass_ee->SetFillColor(kOrange);
    h_DY_mass_ee->SetLineColor(kOrange);
    h_DY_mass_ee->SetDirectory(0);
    s_mass_ee->Add(h_DY_mass_ee);
    s_mass_ee_wFR->Add(h_DY_mass_ee);
    h_DY_mass_ee2->SetFillColor(kOrange);
    h_DY_mass_ee2->SetLineColor(kOrange);
    h_DY_mass_ee2->SetDirectory(0);
    s_mass_ee2->Add(h_DY_mass_ee2);

//--------------------------------------- DATA -----------------------------------------------------

    TH1D *h_data_mass_ee, *h_data_mass_ee2;

    Mgr.SetProc(_EE_DoubleEG_Full);
    f_data_ee->GetObject("h_mass_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_ee);
    f_data_ee->GetObject("h_mass2_"+Mgr.Procname[Mgr.CurrentProc], h_data_mass_ee2);
    h_data_mass_ee->SetMarkerStyle(kFullDotLarge);
    h_data_mass_ee->SetMarkerColor(kBlack);
    h_data_mass_ee->SetLineColor(kBlack);
    h_data_mass_ee->SetDirectory(0);
    h_data_mass_ee2->SetMarkerStyle(kFullDotLarge);
    h_data_mass_ee2->SetMarkerColor(kBlack);
    h_data_mass_ee2->SetLineColor(kBlack);
    h_data_mass_ee2->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

// SYSTEMATICS
    TVectorD estSyst_ee = *((TVectorD*)(f_bkg_est_ee->Get("estSystematics")));
    TVectorD estSyst_ee2 = *((TVectorD*)(f_bkg_est_ee->Get("estSystematics2")));
    Double_t systematics_ee[43], systematics_ee2[86], estSystematics_ee[43], estSystematics_ee2[86];
    for (int i=0; i<86; i++)
    {
        if (i < 43)
        {
            estSystematics_ee[i] = estSyst_ee[i];
            systematics_ee[i] = estSyst_ee[i] / ((TH1D*)(s_mass_ee->GetStack()->Last()))->GetBinContent(i+1);
        }
        systematics_ee2[i] = estSyst_ee2[i] / ((TH1D*)(s_mass_ee2->GetStack()->Last()))->GetBinContent(i+1);
        estSystematics_ee2[i] = estSyst_ee2[i];
    }


    myRatioPlot_t *RP_mass_ee = new myRatioPlot_t("RP_mass_ee", s_mass_ee, h_data_mass_ee);
    myRatioPlot_t *RP_mass_ee2 = new myRatioPlot_t("RP_mass_ee2", s_mass_ee2, h_data_mass_ee2);
    myRatioPlot_t *RP_mass_ee_wFR = new myRatioPlot_t("RP_mass_ee_wFR", s_mass_ee_wFR, h_data_mass_ee);
    RP_mass_ee->SetPlots("m_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)   ");
//    RP_mass_ee->SetPlots("m_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Eksp./(MC+i_{#kern[-0.65]{#lower[-0.3]{#scale[0.7]{c}}}}v.)      ");
    RP_mass_ee2->SetPlots("m_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)   ");
//    RP_mass_ee2->SetPlots("m_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Eksp./(MC+i_{#kern[-0.8]{#lower[-0.3]{#scale[0.7]{c}}}}v.)     ");
    RP_mass_ee_wFR->SetPlots("m_{#font[12]{#lower[-0.2]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)   ");
    RP_mass_ee->SetSystematics(NULL, estSystematics_ee, systematics_ee);
    RP_mass_ee2->SetSystematics(NULL, estSystematics_ee2, systematics_ee2);
    if (FR_systErr > 0)
    {
        TH1D *h_QCD_systErr, *h_WJets_systErr;
        f_QCD_est_ee->GetObject("h_QCD_fullsysterr", h_QCD_systErr);
        f_WJets_est_ee->GetObject("h_WJET_fullsysterr", h_WJets_systErr);
        for (Int_t j=0; j<43; j++)
        {
            estSystematics_ee[j] = sqrt(h_QCD_systErr->GetBinContent(j+1)*h_QCD_systErr->GetBinContent(j+1)+
                                        h_WJets_systErr->GetBinContent(j+1)*h_WJets_systErr->GetBinContent(j+1)+
                                        estSystematics_ee[j]*estSystematics_ee[j]);
            systematics_ee[j] = estSystematics_ee[j] / ((TH1D*)(s_mass_ee->GetStack()->Last()))->GetBinContent(j+1);
        }
        RP_mass_ee_wFR->SetSystematics(NULL, estSystematics_ee, systematics_ee);
    }

    TLegend *legend_ee = new TLegend(0.8, 0.45, 0.95, 0.95);
    legend_ee->AddEntry(h_data_mass_ee, "Data", "lp");
    legend_ee->AddEntry(h_DY_mass_ee, "DY#rightarrow #font[12]{ee} (MC)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[8], "DY#rightarrow #tau#tau (est.)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (est.)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (est.)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (est.)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[3], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[2], "#font[12]{#scale[1.1]{WW}} (est.)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[1], "#font[12]{#scale[1.1]{W}}+Jets (MC)", "f");
    legend_ee->AddEntry(h_bkg_mass_ee[0], "#font[12]{#scale[1.1]{QCD}} (MC)", "f");

    TLegend *legend_ee_wFR = new TLegend(0.8, 0.45, 0.95, 0.95);
    legend_ee_wFR->AddEntry(h_data_mass_ee, "Data", "lp");
    legend_ee_wFR->AddEntry(h_DY_mass_ee, "DY#rightarrow#mu#mu (MC)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[8], "DY#rightarrow #tau#tau (est.)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (est.)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (est.)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (est.)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[3], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_ee_wFR->AddEntry(h_bkg_mass_ee[2], "#font[12]{#scale[1.1]{WW}} (est.)", "f");
    legend_ee_wFR->AddEntry(h_fakes_mass_ee[1], "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
    legend_ee_wFR->AddEntry(h_fakes_mass_ee[0], "#font[12]{#scale[1.1]{QCD}} (est.)", "f");

    RP_mass_ee->ImportLegend(legend_ee);
    RP_mass_ee2->ImportLegend(legend_ee);
    RP_mass_ee_wFR->ImportLegend(legend_ee_wFR);
    RP_mass_ee->Draw(0.5, 1e7, 1);
    RP_mass_ee2->Draw(0.5, 1e7, 1);
    RP_mass_ee_wFR->Draw(0.5, 1e7, 1);

    Double_t dataerror_ee, MCerror_ee, MCerror_ee_wFR, dataintegral_ee=1.3107e+07, MCintegral_ee, MCintegral_ee_wFR;
    Double_t dataerrorZ_ee, MCerrorZ_ee, DYerrorZ_ee, dataintegralZ_ee, MCintegralZ_ee, DYintegralZ_ee;
    Double_t dataerror_noZ_ee=0, MCerror_noZ_ee=0, dataintegral_noZ_ee=0, MCintegral_noZ_ee, temp_noZ_ee;
    Double_t qcdMCerror_ee=0, wjetsMCerror_ee=0, qcdMCintegral_ee=0, wjetsMCintegral_ee=0;
    Double_t qcdDDerror_ee=0, wjetsDDerror_ee=0, qcdDDintegral_ee=0, wjetsDDintegral_ee=0;

    dataintegral_ee = h_data_mass_ee->IntegralAndError(1, h_data_mass_ee->GetSize()-2, dataerror_ee);
    MCintegral_ee = ((TH1D*)(s_mass_ee->GetStack()->Last()))->IntegralAndError(1, h_data_mass_ee->GetSize()-2, MCerror_ee);
    MCintegral_ee_wFR = ((TH1D*)(s_mass_ee_wFR->GetStack()->Last()))->IntegralAndError(1, h_data_mass_ee->GetSize()-2, MCerror_ee_wFR);

    dataintegralZ_ee = h_data_mass_ee->IntegralAndError(10, 22, dataerrorZ_ee);
    MCintegralZ_ee = ((TH1D*)(s_mass_ee->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ_ee);

    dataintegral_noZ_ee = h_data_mass_ee->IntegralAndError(1, 9, temp_noZ_ee);
    dataerror_noZ_ee += temp_noZ_ee * temp_noZ_ee;
    dataintegral_noZ_ee += h_data_mass_ee->IntegralAndError(23, h_data_mass_ee->GetSize()-2, temp_noZ_ee);
    dataerror_noZ_ee += temp_noZ_ee * temp_noZ_ee;
    dataerror_noZ_ee = sqrt(dataerror_noZ_ee);

    MCintegral_noZ_ee = ((TH1D*)(s_mass_ee->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ_ee);
    MCerror_noZ_ee += temp_noZ_ee * temp_noZ_ee;
    MCintegral_noZ_ee += ((TH1D*)(s_mass_ee->GetStack()->Last()))->IntegralAndError(23, h_data_mass_ee->GetSize()-2, temp_noZ_ee);
    MCerror_noZ_ee += temp_noZ_ee * temp_noZ_ee;
    MCerror_noZ_ee = sqrt(MCerror_noZ_ee);

    qcdMCintegral_ee = h_bkg_mass_ee[0]->IntegralAndError(1, h_bkg_mass_ee[0]->GetSize()-2, qcdMCerror_ee);
    wjetsMCintegral_ee = h_bkg_mass_ee[1]->IntegralAndError(1, h_bkg_mass_ee[1]->GetSize()-2, wjetsMCerror_ee);
    qcdDDintegral_ee = h_fakes_mass_ee[0]->IntegralAndError(1, h_fakes_mass_ee[0]->GetSize()-2, qcdDDerror_ee);
    wjetsDDintegral_ee = h_fakes_mass_ee[1]->IntegralAndError(1, h_fakes_mass_ee[1]->GetSize()-2, wjetsDDerror_ee);

    std::cout << "ee Data events: " << dataintegral_ee << "+-" << dataerror_ee << endl;
    std::cout << "ee MC+DD events (only EMu): " << MCintegral_ee << "+-" << MCerror_ee << endl;
    std::cout << "ee MC/Obs (only EMu): " << MCintegral_ee/dataintegral_ee << "+-" <<
                 sqrt((dataerror_ee / dataintegral_ee) * (dataerror_ee / dataintegral_ee) +
                       (MCerror_ee / MCintegral_ee) * (MCerror_ee / MCintegral_ee)) << endl;
    std::cout << "ee Avg. Data and Est relative difference (only EMu): " << CompAvgDataMCDifference(h_data_mass_ee, s_mass_ee) << endl;
    std::cout << "ee Chi^2 (only EMu): " << CompChiSquared(h_data_mass_ee, s_mass_ee) << endl << endl;

    std::cout << "ee Data events around Z (only EMu): " << dataintegralZ_ee << "+-" << dataerrorZ_ee << endl;
    std::cout << "ee MC events around Z (only EMu): " << MCintegralZ_ee << "+-" << MCerrorZ_ee << endl;
    std::cout << "ee Data events outside Z (only EMu): " << dataintegral_noZ_ee << "+-" << dataerror_noZ_ee << endl;
    std::cout << "ee MC events outside Z (only EMu): " << MCintegral_noZ_ee << "+-" << MCerror_noZ_ee << endl << endl;

    std::cout << "ee MC+DD events (EMu+FR): " << MCintegral_ee_wFR << "+-" << MCerror_ee_wFR << endl;
    std::cout << "ee MC/Obs (EMu+FR): " << MCintegral_ee_wFR/dataintegral_ee << "+-" <<
                 sqrt((dataerror_ee / dataintegral_ee) * (dataerror_ee / dataintegral_ee) +
                       (MCerror_ee_wFR / MCintegral_ee_wFR) * (MCerror_ee_wFR / MCintegral_ee_wFR)) << endl;

    std::cout << "ee QCD events (MC): " << qcdMCintegral_ee << "+-" << qcdMCerror_ee << endl;
    std::cout << "ee W+Jets events (MC): " << wjetsMCintegral_ee << "+-" << wjetsMCerror_ee << endl;
    std::cout << "ee QCD events (FR): " << qcdDDintegral_ee << "+-" << qcdDDerror_ee << endl;
    std::cout << "ee W+Jets events (FR): " << wjetsDDintegral_ee << "+-" << wjetsDDerror_ee << endl << endl;

//############################## MUON CHANNEL ###############################################

    Mgr.SetProc(_MuMu_DY_Full);
    TString name_DY_mumu = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY_mumu = new TFile(name_DY_mumu, "READ");
    cout << "Hists location: " << Mgr.HistLocation << endl;
    if (f_DY_mumu->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_Bkg_Full);
    TString name_bkg_mumu = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg_mumu = new TFile(name_bkg_mumu, "READ");
    if (f_bkg_mumu->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    TString name_bkg_est_mumu = Mgr.HistLocation+"EstBkg_MuMu.root";
    TFile* f_bkg_est_mumu = new TFile(name_bkg_est_mumu, "READ");
    if (f_bkg_est_mumu->IsOpen()) std::cout << "File " << "EstBkg_MuMu.root" << " opened successfully" << endl;
    TString name_WJets_est_mumu = Mgr.HistLocation+"EstWJets_MuMu.root";
    TFile* f_WJets_est_mumu = new TFile(name_WJets_est_mumu, "READ");
    if (f_WJets_est_mumu->IsOpen()) std::cout << "File " << "EstWJets_MuMu.root" << " opened successfully" << endl;
    TString name_QCD_est_mumu = Mgr.HistLocation+"EstQCD_MuMu.root";
    TFile* f_QCD_est_mumu = new TFile(name_QCD_est_mumu, "READ");
    if (f_QCD_est_mumu->IsOpen()) std::cout << "File " << "EstQCD_MuMu.root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_SingleMuon_Full);
    TString name_data_mumu = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data_mumu = new TFile(name_data_mumu, "READ");
    if (f_data_mumu->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully\n" << endl;

    THStack *s_mass_mumu = new THStack("s_mass_mumu", "");
    THStack *s_mass_mumu_wFR = new THStack("s_mass_mumu_wFR", "");
    THStack *s_mass_mumu2 = new THStack("s_mass_mumu2", "");

    //--------------------------------- FakeMu bkg -----------------------------------------------------
    TH1D *h_fakes_mass_mumu[2];
    f_QCD_est_mumu->GetObject("h_QCD_est", h_fakes_mass_mumu[0]);
    removeNegativeBins(h_fakes_mass_mumu[0]);
    h_fakes_mass_mumu[0]->SetFillColor(kRed + 3);
    h_fakes_mass_mumu[0]->SetLineColor(kRed + 3);
    h_fakes_mass_mumu[0]->SetDirectory(0);
    s_mass_mumu_wFR->Add(h_fakes_mass_mumu[0]);

    f_WJets_est_mumu->GetObject("h_WJET_fit_SS", h_fakes_mass_mumu[1]);
    removeNegativeBins(h_fakes_mass_mumu[1]);
    h_fakes_mass_mumu[1]->SetFillColor(kRed - 2);
    h_fakes_mass_mumu[1]->SetLineColor(kRed - 2);
    h_fakes_mass_mumu[1]->SetDirectory(0);
    h_fakes_mass_mumu[1]->GetYaxis()->SetRangeUser(0.01,1e7);
    s_mass_mumu_wFR->Add(h_fakes_mass_mumu[1]);

    //----------------------------------- MC bkg -------------------------------------------------------
    TH1D *h_bkg_mass_mumu[9], *h_bkg_mass_mumu2[9];
    iter = 0;

    for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
//        if (pr == _MuMu_QCDMuEnriched_Full)
//        {
//            iter++;
//            continue;
//        }

//        if (pr == _MuMu_WJets || pr == _MuMu_QCDMuEnriched_Full)
        if (pr == _MuMu_WJets_Full || pr == _MuMu_QCDMuEnriched_Full || pr == _MuMu_WZ || pr == _MuMu_ZZ)
        {
            f_bkg_mumu->GetObject("h_mass_"+Mgr.Procname[pr], h_bkg_mass_mumu[iter]);
            f_bkg_mumu->GetObject("h_mass2_"+Mgr.Procname[pr], h_bkg_mass_mumu2[iter]);
        }
        else
        {
            f_bkg_est_mumu->GetObject("h_MuMu_mass_Est_"+Mgr.Procname[pr], h_bkg_mass_mumu[iter]);
            f_bkg_est_mumu->GetObject("h_MuMu_mass_Est2_"+Mgr.Procname[pr], h_bkg_mass_mumu2[iter]);
        }
        removeNegativeBins(h_bkg_mass_mumu[iter]);
        removeNegativeBins(h_bkg_mass_mumu2[iter]);

        Color_t color = kBlack;
        if (pr == _MuMu_QCDMuEnriched_Full) color = kRed + 3;
        if (pr == _MuMu_WJets_Full) color = kRed - 2;
        if (pr == _MuMu_WW) color = kMagenta - 5;
        if (pr == _MuMu_WZ) color = kMagenta - 2;
        if (pr == _MuMu_ZZ) color = kMagenta - 6;
        if (pr == _MuMu_tbarW) color = kGreen - 2;
        if (pr == _MuMu_tW) color = kGreen + 2;
        if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
        if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

        h_bkg_mass_mumu[iter]->SetFillColor(color);
        h_bkg_mass_mumu[iter]->SetLineColor(color);
        h_bkg_mass_mumu[iter]->SetDirectory(0);
        s_mass_mumu->Add(h_bkg_mass_mumu[iter]);
        h_bkg_mass_mumu2[iter]->SetFillColor(color);
        h_bkg_mass_mumu2[iter]->SetLineColor(color);
        h_bkg_mass_mumu2[iter]->SetDirectory(0);
        s_mass_mumu2->Add(h_bkg_mass_mumu2[iter]);
        
        if ((pr < _EndOf_MuMu_VVnST_Normal || pr > _EndOf_MuMu_MCbkg_Normal) && pr != _MuMu_QCDMuEnriched_Full && pr != _MuMu_WJets_Full)
            s_mass_mumu_wFR->Add(h_bkg_mass_mumu[iter]);

        iter++;

        if (pr == _MuMu_WJets_Full)
            pr = _EndOf_MuMu_VVnST_Normal; // next - WW
        if (pr == _MuMu_tW)
            pr = _MuMu_VVnST; // next - ttbar
        if (pr == _MuMu_DYTauTau_Full) // last
            break;

    } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

    TH1D *h_DY_mass_mumu, *h_DY_mass_mumu2;
    f_DY_mumu->GetObject("h_mass_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_mumu);
    f_DY_mumu->GetObject("h_mass2_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_mumu2);
    removeNegativeBins(h_DY_mass_mumu);
    removeNegativeBins(h_DY_mass_mumu2);
    h_DY_mass_mumu->SetFillColor(kOrange);
    h_DY_mass_mumu->SetLineColor(kOrange);
    h_DY_mass_mumu->SetDirectory(0);
    s_mass_mumu->Add(h_DY_mass_mumu);
    h_DY_mass_mumu2->SetFillColor(kOrange);
    h_DY_mass_mumu2->SetLineColor(kOrange);
    h_DY_mass_mumu2->SetDirectory(0);
    s_mass_mumu2->Add(h_DY_mass_mumu2);

    s_mass_mumu_wFR->Add(h_DY_mass_mumu);

//--------------------------------------- DATA -----------------------------------------------------

    TH1D *h_data_mass_mumu, *h_data_mass_mumu2;
    f_data_mumu->GetObject("h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_mumu);
    f_data_mumu->GetObject("h_mass2_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_mumu2);

    h_data_mass_mumu->SetMarkerStyle(kFullDotLarge);
    h_data_mass_mumu->SetMarkerColor(kBlack);
    h_data_mass_mumu->SetLineColor(kBlack);
    h_data_mass_mumu->SetDirectory(0);
    h_data_mass_mumu2->SetMarkerStyle(kFullDotLarge);
    h_data_mass_mumu2->SetMarkerColor(kBlack);
    h_data_mass_mumu2->SetLineColor(kBlack);
    h_data_mass_mumu2->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

// SYSTEMATICS
    TVectorD estSyst_mumu = *((TVectorD*)(f_bkg_est_mumu->Get("estSystematics")));
    TVectorD estSyst_mumu2 = *((TVectorD*)(f_bkg_est_mumu->Get("estSystematics2")));
    Double_t systematics_mumu[43], systematics_mumu2[86], estSystematics_mumu[43], estSystematics_mumu2[86];
    for (int i=0; i<86; i++)
    {
        if (i < 43)
        {
            estSystematics_mumu[i] = estSyst_mumu[i];
            systematics_mumu[i] = estSyst_mumu[i] / ((TH1D*)(s_mass_mumu->GetStack()->Last()))->GetBinContent(i+1);
        }
        systematics_mumu2[i] = estSyst_mumu2[i] / ((TH1D*)(s_mass_mumu2->GetStack()->Last()))->GetBinContent(i+1);
        estSystematics_mumu2[i] = estSyst_mumu2[i];
    }

    myRatioPlot_t *RP_mass_mumu = new myRatioPlot_t("RP_mass_mumu", s_mass_mumu, h_data_mass_mumu);
    myRatioPlot_t *RP_mass_mumu2 = new myRatioPlot_t("RP_mass_mumu2", s_mass_mumu2, h_data_mass_mumu2);
    myRatioPlot_t *RP_mass_mumu_wFR = new myRatioPlot_t("RP_mass_mumu_wFR", s_mass_mumu_wFR, h_data_mass_mumu);

    RP_mass_mumu->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)   ");
    RP_mass_mumu2->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)    ");
    RP_mass_mumu_wFR->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Data/(MC+est)");
    RP_mass_mumu->SetSystematics(NULL, estSystematics_mumu, systematics_mumu);
    RP_mass_mumu2->SetSystematics(NULL, estSystematics_mumu2, systematics_mumu2);
    if (FR_systErr > 0)
    {
        TH1D *h_QCD_systErr, *h_WJets_systErr;
        f_QCD_est_mumu->GetObject("h_QCD_fullsysterr", h_QCD_systErr);
        f_WJets_est_mumu->GetObject("h_WJET_fullsysterr", h_WJets_systErr);
        for (Int_t j=0; j<43; j++)
        {
            estSystematics_mumu[j] = sqrt(h_QCD_systErr->GetBinContent(j+1)*h_QCD_systErr->GetBinContent(j+1)+
                                          h_WJets_systErr->GetBinContent(j+1)*h_WJets_systErr->GetBinContent(j+1)+
                                          estSystematics_mumu[j]*estSystematics_mumu[j]);
            systematics_mumu[j] = estSystematics_mumu[j] / ((TH1D*)(s_mass_mumu->GetStack()->Last()))->GetBinContent(j+1);
        }
        RP_mass_mumu_wFR->SetSystematics(NULL, estSystematics_mumu, systematics_mumu);
    }

    TLegend *legend_mumu = new TLegend(0.8, 0.45, 0.95, 0.95);
    TLegend *legend_mumu_wFR = new TLegend(0.8, 0.45, 0.95, 0.95);

    // Legend
    legend_mumu->AddEntry(h_data_mass_mumu, "Data", "lp");
    legend_mumu->AddEntry(h_DY_mass_mumu, "DY#rightarrow#mu#mu (MC)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[8], "DY#rightarrow #tau#tau (est.)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (est.)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (est.)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (est.)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[3], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[2], "#font[12]{#scale[1.1]{WW}} (est.)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[1], "#font[12]{#scale[1.1]{W}}+Jets (MC)", "f");
    legend_mumu->AddEntry(h_bkg_mass_mumu[0], "#font[12]{#scale[1.1]{QCD}} (MC)", "f");

    legend_mumu_wFR->AddEntry(h_data_mass_mumu, "Data", "lp");
    legend_mumu_wFR->AddEntry(h_DY_mass_mumu, "DY#rightarrow#mu#mu (MC)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[8], "DY#rightarrow #tau#tau (est.)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (est.)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (est.)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (est.)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[3], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_mumu_wFR->AddEntry(h_bkg_mass_mumu[2], "#font[12]{#scale[1.1]{WW}} (est.)", "f");
    legend_mumu_wFR->AddEntry(h_fakes_mass_mumu[1], "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
    legend_mumu_wFR->AddEntry(h_fakes_mass_mumu[0], "#font[12]{#scale[1.1]{QCD}} (est.)", "f");

    RP_mass_mumu->ImportLegend(legend_mumu);
    RP_mass_mumu->Draw(0.5, 1e7, 1);
    RP_mass_mumu2->ImportLegend(legend_mumu);
    RP_mass_mumu2->Draw(0.5, 1e7, 1);
    RP_mass_mumu_wFR->ImportLegend(legend_mumu_wFR);
    RP_mass_mumu_wFR->Draw(0.5, 1e7, 1);

    Double_t dataerror_mumu, MCerror_mumu, MCerror_mumu_wFR, dataintegral_mumu=2.25081e+07, MCintegral_mumu, MCintegral_mumu_wFR;
    Double_t dataerrorZ_mumu, MCerrorZ_mumu, DYerrorZ_mumu, dataintegralZ_mumu, MCintegralZ_mumu, DYintegralZ_mumu;
    Double_t dataerror_noZ_mumu=0, MCerror_noZ_mumu=0, dataintegral_noZ_mumu=0, MCintegral_noZ_mumu, temp_noZ_mumu;
    Double_t qcdMCerror_mumu=0, wjetsMCerror_mumu=0, qcdMCintegral_mumu=0, wjetsMCintegral_mumu=0;
    Double_t qcdDDerror_mumu=0, wjetsDDerror_mumu=0, qcdDDintegral_mumu=0, wjetsDDintegral_mumu=0;

    dataintegral_mumu = h_data_mass_mumu->IntegralAndError(1, h_data_mass_mumu->GetSize()-2, dataerror_mumu);
    MCintegral_mumu = ((TH1D*)(s_mass_mumu->GetStack()->Last()))->IntegralAndError(1, h_data_mass_mumu->GetSize()-2, MCerror_mumu);
    MCintegral_mumu_wFR = ((TH1D*)(s_mass_mumu_wFR->GetStack()->Last()))->IntegralAndError(1, h_data_mass_mumu->GetSize()-2, MCerror_mumu_wFR);

    dataintegralZ_mumu = h_data_mass_mumu->IntegralAndError(10, 22, dataerrorZ_mumu);
    MCintegralZ_mumu = ((TH1D*)(s_mass_mumu->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ_mumu);

    dataintegral_noZ_mumu = h_data_mass_mumu->IntegralAndError(1, 9, temp_noZ_mumu);
    dataerror_noZ_mumu += temp_noZ_mumu * temp_noZ_mumu;
    dataintegral_noZ_mumu += h_data_mass_mumu->IntegralAndError(23, h_data_mass_mumu->GetSize()-2, temp_noZ_mumu);
    dataerror_noZ_mumu += temp_noZ_mumu * temp_noZ_mumu;
    dataerror_noZ_mumu = sqrt(dataerror_noZ_mumu);

    MCintegral_noZ_mumu = ((TH1D*)(s_mass_mumu->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ_mumu);
    MCerror_noZ_mumu += temp_noZ_mumu * temp_noZ_mumu;
    MCintegral_noZ_mumu += ((TH1D*)(s_mass_mumu->GetStack()->Last()))->IntegralAndError(23, h_data_mass_mumu->GetSize()-2, temp_noZ_mumu);
    MCerror_noZ_mumu += temp_noZ_mumu * temp_noZ_mumu;
    MCerror_noZ_mumu = sqrt(MCerror_noZ_mumu);

    qcdMCintegral_mumu = h_bkg_mass_mumu[0]->IntegralAndError(1, h_bkg_mass_mumu[0]->GetSize()-2, qcdMCerror_mumu);
    wjetsMCintegral_mumu = h_bkg_mass_mumu[1]->IntegralAndError(1, h_bkg_mass_mumu[1]->GetSize()-2, wjetsMCerror_mumu);
    qcdDDintegral_mumu = h_fakes_mass_mumu[0]->IntegralAndError(1, h_fakes_mass_mumu[0]->GetSize()-2, qcdDDerror_mumu);
    wjetsDDintegral_mumu = h_fakes_mass_mumu[1]->IntegralAndError(1, h_fakes_mass_mumu[1]->GetSize()-2, wjetsDDerror_mumu);

    std::cout << "MuMu Data events: " << dataintegral_mumu << "+-" << dataerror_mumu << endl;
    std::cout << "MuMu MC+DD events (only EMu): " << MCintegral_mumu << "+-" << MCerror_mumu << endl;
    std::cout << "MuMu MC/Obs (only EMu): " << MCintegral_mumu/dataintegral_mumu << "+-" <<
                 sqrt((dataerror_mumu / dataintegral_mumu) * (dataerror_mumu / dataintegral_mumu) +
                       (MCerror_mumu / MCintegral_mumu) * (MCerror_mumu / MCintegral_mumu)) << endl;
    std::cout << "MuMu Avg. Data and Est relative difference (only EMu): " << CompAvgDataMCDifference(h_data_mass_mumu, s_mass_mumu) << endl;
    std::cout << "MuMu Chi^2 (only EMu): " << CompChiSquared(h_data_mass_mumu, s_mass_mumu) << endl << endl;

    std::cout << "MuMu Data events around Z (only EMu): " << dataintegralZ_mumu << "+-" << dataerrorZ_mumu << endl;
    std::cout << "MuMu MC events around Z (only EMu): " << MCintegralZ_mumu << "+-" << MCerrorZ_mumu << endl;
    std::cout << "MuMu Data events outside Z (only EMu): " << dataintegral_noZ_mumu << "+-" << dataerror_noZ_mumu << endl;
    std::cout << "MuMu MC events outside Z (only EMu): " << MCintegral_noZ_mumu << "+-" << MCerror_noZ_mumu << endl << endl;

    std::cout << "MuMu MC+DD events (EMu+FR): " << MCintegral_mumu_wFR << "+-" << MCerror_mumu_wFR << endl;
    std::cout << "MuMu MC/Obs (EMu+FR): " << MCintegral_mumu_wFR/dataintegral_mumu << "+-" <<
                 sqrt((dataerror_mumu / dataintegral_mumu) * (dataerror_mumu / dataintegral_mumu) +
                       (MCerror_mumu_wFR / MCintegral_mumu_wFR) * (MCerror_mumu_wFR / MCintegral_mumu_wFR)) << endl;

    std::cout << "MuMu QCD events (MC): " << qcdMCintegral_mumu << "+-" << qcdMCerror_mumu << endl;
    std::cout << "MuMu W+Jets events (MC): " << wjetsMCintegral_mumu << "+-" << wjetsMCerror_mumu << endl;
    std::cout << "MuMu QCD events (FR): " << qcdDDintegral_mumu << "+-" << qcdDDerror_mumu << endl;
    std::cout << "MuMu W+Jets events (FR): " << wjetsDDintegral_mumu << "+-" << wjetsDDerror_mumu << endl;

    f_DY_ee->Close();
    f_DY_mumu->Close();
    f_bkg_ee->Close();
    f_bkg_mumu->Close();
    f_bkg_est_ee->Close();
    f_bkg_est_mumu->Close();
    f_data_ee->Close();
    f_data_mumu->Close();

} // End of Est_HistDrawer()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Test_RocCorr()
{
    LocalFileMgr Mgr;

    Mgr.SetProc(_MuMu_DY_Full);
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY = new TFile(name_DY, "READ");
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_Bkg_Full);
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile(name_bkg, "READ");
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_SingleMuon_Full);
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile(name_data, "READ");
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

    THStack *s_mass_before_RoccoR = new THStack("s_mass_before_RoccoR", "");
    THStack *s_mass_before_EffCorr = new THStack("s_mass_before_EffCorr", "");

//----------------------------------- MC bkg -------------------------------------------------------
    TH1D *h_bkg_mass_before_RoccoR[9], *h_bkg_mass_before_EffCorr[9];
    Int_t iter = 0;

    for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
        if (pr == _MuMu_QCDMuEnriched_Full)
        {
            iter++;
            continue;
        }
        f_bkg->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[pr], h_bkg_mass_before_RoccoR[iter]);
        f_bkg->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_before_EffCorr[iter]);

        removeNegativeBins(h_bkg_mass_before_RoccoR[iter]);
        removeNegativeBins(h_bkg_mass_before_EffCorr[iter]);

        h_bkg_mass_before_RoccoR[iter]->SetFillColor(kRed);
        h_bkg_mass_before_RoccoR[iter]->SetLineColor(kRed);
        h_bkg_mass_before_RoccoR[iter]->SetDirectory(0);
        s_mass_before_RoccoR->Add(h_bkg_mass_before_RoccoR[iter]);

        h_bkg_mass_before_EffCorr[iter]->SetMarkerStyle(kFullDotLarge);
        h_bkg_mass_before_EffCorr[iter]->SetMarkerColor(kBlack);
        h_bkg_mass_before_EffCorr[iter]->SetLineColor(kBlack);
        h_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
        s_mass_before_EffCorr->Add(h_bkg_mass_before_EffCorr[iter]);

        iter++;

        if (pr == _MuMu_WJets_Full)
            pr = _EndOf_MuMu_VVnST_Normal; // next - WW
        if (pr == _MuMu_tW)
            pr = _MuMu_VVnST; // next - ttbar
        if (pr == _MuMu_DYTauTau_Full) // last
            break;

    } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

    TH1D *h_DY_mass_before_RoccoR, *h_DY_mass_before_EffCorr;

    f_DY->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_RoccoR);
    f_DY->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_EffCorr);

    removeNegativeBins(h_DY_mass_before_RoccoR);
    removeNegativeBins(h_DY_mass_before_EffCorr);

    h_DY_mass_before_RoccoR->SetFillColor(kRed);
    h_DY_mass_before_RoccoR->SetLineColor(kRed);
    h_DY_mass_before_RoccoR->SetDirectory(0);
    s_mass_before_RoccoR->Add(h_DY_mass_before_RoccoR);

    h_DY_mass_before_EffCorr->SetMarkerStyle(kFullDotLarge);
    h_DY_mass_before_EffCorr->SetMarkerColor(kBlack);
    h_DY_mass_before_EffCorr->SetLineColor(kBlack);
    h_DY_mass_before_EffCorr->SetDirectory(0);
    s_mass_before_EffCorr->Add(h_DY_mass_before_EffCorr);

//--------------------------------------- DATA -----------------------------------------------------

    TH1D *h_data_mass_before_RoccoR, *h_data_mass;

    f_data->GetObject("h_mass_before_RocCorr_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass_before_RoccoR);
    f_data->GetObject("h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass);

    h_data_mass_before_RoccoR->SetFillColor(kRed);
    h_data_mass_before_RoccoR->SetLineColor(kRed);
    h_data_mass_before_RoccoR->SetDirectory(0);

    h_data_mass->SetMarkerStyle(kFullDotLarge);
    h_data_mass->SetMarkerColor(kBlack);
    h_data_mass->SetLineColor(kBlack);
    h_data_mass->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

    myRatioPlot_t *RP_mass_test_MC, *RP_mass_test_Data, *RP_mass_test_uncorrMC_corrData, *RP_mass_test_corrMC_uncorrData;
    RP_mass_test_MC = new myRatioPlot_t("RP_mass_test_MC",
                                        ((TH1D*)(s_mass_before_RoccoR->GetStack()->Last())),
                                        ((TH1D*)(s_mass_before_EffCorr->GetStack()->Last()))
                                       );
    RP_mass_test_Data = new myRatioPlot_t("RP_mass_test_Data", h_data_mass_before_RoccoR, h_data_mass);
    RP_mass_test_uncorrMC_corrData = new myRatioPlot_t("RP_mass_test_uncorrMC_corrData",
                                                       ((TH1D*)(s_mass_before_RoccoR->GetStack()->Last())),
                                                       h_data_mass);
    RP_mass_test_corrMC_uncorrData = new myRatioPlot_t("RP_mass_test_corrMC_uncorrData",
                                                       ((TH1D*)(s_mass_before_EffCorr->GetStack()->Last())),
                                                       h_data_mass_before_RoccoR);

    RP_mass_test_MC->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] MC", 15, 3000, "Corr/Uncorr");
    RP_mass_test_Data->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] Data", 15, 3000, "Corr/Uncorr");
    RP_mass_test_uncorrMC_corrData->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] MC=uncorr, Data=Corr", 15, 3000, "Corr/Uncorr");
    RP_mass_test_corrMC_uncorrData->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] MC=Corr, Data=uncorr", 15, 3000, "Corr/Uncorr");

    TLegend *legend = new TLegend(0.8, 0.45, 0.95, 0.95);
    legend->AddEntry(h_data_mass, "After RoccoR", "lp");
    legend->AddEntry(h_DY_mass_before_RoccoR, "Before RoccoR", "f");

    TLegend *legend_UC = new TLegend(0.8, 0.45, 0.95, 0.95);
    legend_UC->AddEntry(h_data_mass, "Corrected Data", "lp");
    legend_UC->AddEntry(h_DY_mass_before_RoccoR, "Uncorrected MC", "f");

    TLegend *legend_CU = new TLegend(0.8, 0.45, 0.95, 0.95);
    legend_CU->AddEntry(h_data_mass_before_RoccoR, "UnCorrected Data", "lp");
    legend_CU->AddEntry(h_DY_mass_before_EffCorr, "Corrected MC", "f");

    RP_mass_test_MC->ImportLegend(legend);
    RP_mass_test_Data->ImportLegend(legend);
    RP_mass_test_uncorrMC_corrData->ImportLegend(legend_UC);
    RP_mass_test_corrMC_uncorrData->ImportLegend(legend_CU);

    RP_mass_test_MC->Draw(0.5, 1e7, 1);
    RP_mass_test_Data->Draw(0.5, 1e7, 1);
    RP_mass_test_uncorrMC_corrData->Draw(0.5, 1e7, 1);
    RP_mass_test_corrMC_uncorrData->Draw(0.5, 1e7, 1);

    f_DY->Close();
    f_bkg->Close();
    f_data->Close();

} // End of Test_RocCorr()


/// ------------------------------- COMP CHI^2 ---------------------------------- ///
Double_t CompChiSquared (TH1D *h_data, THStack *s_MC)
{
    Double_t ChiSquared = 0;
    Int_t size_data = h_data->GetSize();
    Int_t size_MC = ((TH1D*)(s_MC->GetStack()->Last()))->GetSize();
    Int_t temp_data, temp_MC;
    if (size_data != size_MC)
    {
        cout << "CompChiSquared: Sizes do not match!" << endl;
        return -999;
    }
    else for (Int_t i=1; i<size_data-1; i++)
    {
        temp_data = h_data->GetBinContent(i);
        temp_MC = ((TH1D*)(s_MC->GetStack()->Last()))->GetBinContent(i);
        if (temp_data != 0)
            ChiSquared += ((temp_MC - temp_data) * (temp_MC - temp_data)) / temp_data;
    }
    return ChiSquared / (size_data - 2);
} // End of CompChiSquared()


/// ------------------------------- COMP AVG |DATA/MC-1| ---------------------------------- ///
Double_t CompAvgDataMCDifference (TH1D *h_data, THStack *s_MC)
{
    TH1D* h_MC = ((TH1D*)(s_MC->GetStack()->Last()));
    Double_t AvgDataMCRatio = 0;
    Int_t size_data = h_data->GetSize();
    Int_t size_MC = h_MC->GetSize();
    if (size_data != size_MC)
    {
        cout << "CompAvgDataMCDifference: Sizes do not match!" << endl;
        return -999;
    }
    else for (Int_t i=1; i<size_data-1; i++)
    {
        if (h_MC->GetBinContent(i) != 0)
            AvgDataMCRatio += fabs(h_data->GetBinContent(i) / h_MC->GetBinContent(i) - 1);
        else size_MC--;
    }
    return AvgDataMCRatio / (size_MC - 2);
} // End of CompChiSquared()


void removeNegativeBins(TH1D *h)
{
    for (int i=0; i<h->GetSize(); i++)
    {
        if (h->GetBinContent(i) < 0)
        {
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }
}


void TEST_HistDrawer (TString whichGraphs , TString type)
{
    if (!whichGraphs.Length())
    {
        cout << "Error: no whichGraphs specified!" << endl;
        return;
    }

    Int_t count_drawn = 0;
    LocalFileMgr Mgr;

    Mgr.SetProc(_MuMu_DY_Full);
    TString name_DY = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root";
    TFile* f_DY = new TFile(name_DY, "READ");
    if (f_DY->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_Bkg_Full);
    TString name_bkg = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root";
    TFile* f_bkg = new TFile(name_bkg, "READ");
    if (f_bkg->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root" << " opened successfully" << endl;
    Mgr.SetProc(_MuMu_SingleMuon_Full);
    TString name_data = Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root";
    TFile* f_data = new TFile(name_data, "READ");
    if (f_data->IsOpen()) std::cout << "File " << "Hist_"+Mgr.Procname[_MuMu_SingleMuon_Full]+".root" << " opened successfully" << endl;

        THStack *s_mass_before_EffCorr = new THStack("s_mass_before_EffCorr", "");
        THStack *s_mass = new THStack("s_mass", "");

//----------------------------------- MC bkg -------------------------------------------------------
        TH1D *h_bkg_mass_before_EffCorr[9], *h_bkg_mass[9];
        Int_t iter = 0;

        for (SelProc_t pr = _MuMu_QCDMuEnriched_Full; pr > _EndOf_MuMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _MuMu_QCDMuEnriched_Full)
            {
                iter++;
                continue;
            }
            f_bkg->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[pr], h_bkg_mass_before_EffCorr[iter]);
            f_bkg->GetObject("h_mass_"+Mgr.Procname[pr], h_bkg_mass[iter]);
            removeNegativeBins(h_bkg_mass_before_EffCorr[iter]);
            removeNegativeBins(h_bkg_mass[iter]);

            // Converting to event density (dividing by bin width)
            for (Int_t i=1; i<=binnum; i++)
            {
                h_bkg_mass_before_EffCorr[iter]->SetBinContent(i, h_bkg_mass_before_EffCorr[iter]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
                h_bkg_mass_before_EffCorr[iter]->SetBinError(i, h_bkg_mass_before_EffCorr[iter]->GetBinError(i)/(massbins[i]-massbins[i-1]));
                h_bkg_mass[iter]->SetBinContent(i, h_bkg_mass[iter]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
                h_bkg_mass[iter]->SetBinError(i, h_bkg_mass[iter]->GetBinError(i)/(massbins[i]-massbins[i-1]));
            }

            Color_t color = kBlack;
            if (pr == _MuMu_QCDMuEnriched_Full) color = kRed + 3;
            if (pr == _MuMu_WJets_Full) color = kRed - 2;
            if (pr == _MuMu_WW) color = kMagenta - 5;
            if (pr == _MuMu_WZ) color = kMagenta - 2;
            if (pr == _MuMu_ZZ) color = kMagenta - 6;
            if (pr == _MuMu_tbarW) color = kGreen - 2;
            if (pr == _MuMu_tW) color = kGreen + 2;
            if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
            if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

            h_bkg_mass_before_EffCorr[iter]->SetFillColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetLineColor(color);
            h_bkg_mass_before_EffCorr[iter]->SetDirectory(0);
            s_mass_before_EffCorr->Add(h_bkg_mass_before_EffCorr[iter]);

            h_bkg_mass[iter]->SetFillColor(color);
            h_bkg_mass[iter]->SetLineColor(color);
            h_bkg_mass[iter]->SetDirectory(0);
            s_mass->Add(h_bkg_mass[iter]);

            iter++;

            if (pr == _MuMu_WJets_Full)
                pr = _EndOf_MuMu_VVnST_Normal; // next - WW
            if (pr == _MuMu_tW)
                pr = _MuMu_VVnST; // next - ttbar
            if (pr == _MuMu_DYTauTau_Full) // last
                break;

        } // End of for(bkg)

//---------------------------------- MC signal -----------------------------------------------------

        TH1D *h_DY_mass_before_EffCorr, *h_DY_mass;

        f_DY->GetObject("h_mass_before_EffCorr_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass_before_EffCorr);
        f_DY->GetObject("h_mass_"+Mgr.Procname[_MuMu_DY_Full], h_DY_mass);
        removeNegativeBins(h_DY_mass_before_EffCorr);
        removeNegativeBins(h_DY_mass);

        // Converting to event density (dividing by bin width)
        for (Int_t i=1; i<=binnum; i++)
        {
            h_DY_mass_before_EffCorr->SetBinContent(i, h_DY_mass_before_EffCorr->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_DY_mass_before_EffCorr->SetBinError(i, h_DY_mass_before_EffCorr->GetBinError(i)/(massbins[i]-massbins[i-1]));
            h_DY_mass->SetBinContent(i, h_DY_mass->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_DY_mass->SetBinError(i, h_DY_mass->GetBinError(i)/(massbins[i]-massbins[i-1]));
        }

        h_DY_mass_before_EffCorr->SetFillColor(kOrange);
        h_DY_mass_before_EffCorr->SetLineColor(kOrange);
        h_DY_mass_before_EffCorr->SetDirectory(0);
        s_mass_before_EffCorr->Add(h_DY_mass_before_EffCorr);

        h_DY_mass->SetFillColor(kOrange);
        h_DY_mass->SetLineColor(kOrange);
        h_DY_mass->SetDirectory(0);
        s_mass->Add(h_DY_mass);

//--------------------------------------- DATA -----------------------------------------------------

        TH1D *h_data_mass;

        f_data->GetObject("h_mass_"+Mgr.Procname[_MuMu_SingleMuon_Full], h_data_mass);

        // Converting to event density (dividing by bin width)
        for (Int_t i=1; i<=binnum; i++)
        {
            h_data_mass->SetBinContent(i, h_data_mass->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_data_mass->SetBinError(i, h_data_mass->GetBinError(i)/(massbins[i]-massbins[i-1]));
        }


        h_data_mass->SetMarkerStyle(kFullDotLarge);
        h_data_mass->SetMarkerColor(kBlack);
        h_data_mass->SetLineColor(kBlack);
        h_data_mass->SetDirectory(0);

//--------------------------------- Ratio Plot --------------------------------------

        myRatioPlot_t *RP_mass_before_EffCorr, *RP_mass;
        RP_mass_before_EffCorr = new myRatioPlot_t("RP_mass_before_EffCorr", s_mass_before_EffCorr, h_data_mass);
        RP_mass = new myRatioPlot_t("RP_mass", s_mass, h_data_mass);

        RP_mass_before_EffCorr->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}] before Efficiency SF", 15, 3000);
        RP_mass->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

        TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

        legend->AddEntry(h_data_mass, "Data", "lp");
//        legend->AddEntry(h_data_mass, "Matavimas", "lp");
        legend->AddEntry(h_DY_mass, "DY#rightarrow#mu#mu", "f");
        legend->AddEntry(h_bkg_mass[8], "DY#rightarrow #tau#tau", "f");
        legend->AddEntry(h_bkg_mass[7], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
        legend->AddEntry(h_bkg_mass[6], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
        legend->AddEntry(h_bkg_mass[5], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
        legend->AddEntry(h_bkg_mass[4], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
        legend->AddEntry(h_bkg_mass[3], "#font[12]{#scale[1.1]{WZ}}", "f");
        legend->AddEntry(h_bkg_mass[2], "#font[12]{#scale[1.1]{WW}}", "f");
        legend->AddEntry(h_bkg_mass[1], "#font[12]{#scale[1.1]{W}}+Jets", "f");
//        legend->AddEntry(h_bkg_mass[0], "#font[12]{#scale[1.1]{QCD}}", "f");
        legend->SetNColumns(2);

        RP_mass_before_EffCorr->ImportLegend(legend);
        RP_mass->ImportLegend(legend);

        RP_mass_before_EffCorr->Draw(0.5, 1e6, 1);
        RP_mass_before_EffCorr->s_stackedProcesses->GetYaxis()->SetTitle("N_{#lower[-0.25]{EVT}}/GeV");
        RP_mass_before_EffCorr->pad1->cd();
        RP_mass_before_EffCorr->s_stackedProcesses->Draw("sameaxis");
        RP_mass_before_EffCorr->legend->Draw();
        RP_mass_before_EffCorr->pad1->Update();
        RP_mass_before_EffCorr->canvas->Update();

        RP_mass->Draw(0.5, 1e6, 1);
        RP_mass->s_stackedProcesses->GetYaxis()->SetTitle("N_{#lower[-0.25]{EVT}}/GeV");
        RP_mass->pad1->cd();
        RP_mass->s_stackedProcesses->Draw("sameaxis");
        RP_mass->legend->Draw();
        RP_mass->pad1->Update();
        RP_mass->canvas->Update();

} // End of TEST_HistDrawer()
