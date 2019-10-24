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
#include "./header/FileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

void E_HistDrawer(Int_t type);
void Mu_HistDrawer(Int_t type);
void Test_HistDrawer(Int_t type);
void Fit_HistDrawer();
Double_t CompChiSquared(TH1D *h_data, THStack *s_MC);
Double_t CompAvgDataMCDifference(TH1D *h_data, THStack *s_MC);
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
const Double_t L_B2F = 19721.0;
const Double_t L_G2H = 16146.0;
const Double_t L_B2H = 35867.0;

void TestVariablesForFR (TString WhichX = "", Int_t type = 2)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    if (type < 1 || type > 2) // 1 -- draw histograms obtained from MakeSelectionForFR.C; 2 -- draw histograms obtained from FRgraphMaker.C
    {
        cout << "Wrong type!" << endl;
        return;
    }
    if (whichX.Contains("E"))
    {
        Xselected++;
        cout << "\n*******      EE_HistDrawer(" << type << ")      *******" << endl;
        E_HistDrawer(type);
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*******     MuMu_HistDrawer(" << type << ")     *******" << endl;
        Mu_HistDrawer(type);
    }
    if (whichX.Contains("TEST"))
    {
        Xselected++;
        cout << "\n*******     Test_HistDrawer(" << type << ")     *******" << endl;
        Test_HistDrawer(type);
    }
    if (whichX.Contains("FIT"))
    {
        Xselected++;
        cout << "\n*******     Fit_HistDrawer()     *******" << endl;
        Fit_HistDrawer();
    }
    if (Xselected == 0) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ############################################################################# ///
/// ----------------------------- Electron Channel ------------------------------ ///
/// ############################################################################# ///
void E_HistDrawer(Int_t type)
{
   return; // NOT READY YET
} // End of EE_HistDrawer()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_HistDrawer(Int_t type)
{
    FileMgr fm;
    THStack *s_PFiso_barrel_deno = new THStack("s_PFiso_barrel_deno", "");
    THStack *s_PFiso_endcap_deno = new THStack("s_PFiso_endcap_deno", "");
    THStack *s_TRKiso_barrel_deno = new THStack("s_TRKiso_barrel_deno", "");
    THStack *s_TRKiso_endcap_deno = new THStack("s_TRKiso_endcap_deno", "");
    THStack *s_PFiso_barrel_nume = new THStack("s_PFiso_barrel_nume", "");
    THStack *s_PFiso_endcap_nume = new THStack("s_PFiso_endcap_nume", "");
    THStack *s_TRKiso_barrel_nume = new THStack("s_TRKiso_barrel_nume", "");
    THStack *s_TRKiso_endcap_nume = new THStack("s_TRKiso_endcap_nume", "");
    THStack *s_pT_barrel = new THStack("s_pT_barrel", "");
    THStack *s_pT_endcap = new THStack("s_pT_endcap", "");
    THStack *s_eta = new THStack("s_eta", "");
    THStack *s_nVTX = new THStack("s_nVTX", "");

    TH1D *h_PFiso_barrel_MC_deno[_EndOf_Data_Special], *h_PFiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_TRKiso_barrel_MC_deno[_EndOf_Data_Special], *h_TRKiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_PFiso_barrel_MC_nume[_EndOf_Data_Special], *h_PFiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_TRKiso_barrel_MC_nume[_EndOf_Data_Special], *h_TRKiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC[_EndOf_Data_Special], *h_pT_endcap_MC[_EndOf_Data_Special],
         *h_eta_MC[_EndOf_Data_Special], *h_nVTX_MC[_EndOf_Data_Special],
         *h_PFiso_barrel_data_deno, *h_PFiso_endcap_data_deno, *h_TRKiso_barrel_data_deno, *h_TRKiso_endcap_data_deno,
         *h_PFiso_barrel_data_nume, *h_PFiso_endcap_data_nume, *h_TRKiso_barrel_data_nume, *h_TRKiso_endcap_data_nume,
         *h_pT_barrel_data, *h_pT_endcap_data, *h_eta_data, *h_nVTX_data;

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr1]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr1]);
        file->GetObject("h_TRKiso_barrel_deno", h_TRKiso_barrel_MC_deno[pr1]);
        file->GetObject("h_TRKiso_endcap_deno", h_TRKiso_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr1]);
        file->GetObject("h_TRKiso_barrel_nume", h_TRKiso_barrel_MC_nume[pr1]);
        file->GetObject("h_TRKiso_endcap_nume", h_TRKiso_endcap_MC_nume[pr1]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC[pr1]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC[pr1]);
        file->GetObject("h_eta_deno", h_eta_MC[pr1]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr1]);

        removeNegativeBins(h_PFiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_TRKiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_TRKiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_TRKiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_TRKiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_pT_barrel_MC[pr1]);
        removeNegativeBins(h_pT_endcap_MC[pr1]);
        removeNegativeBins(h_eta_MC[pr1]);
        removeNegativeBins(h_nVTX_MC[pr1]);

        Color_t color = kBlack;
        if (pr1 == _WJets || pr1 == _WJets_ext2v5) color = kRed - 2;
        if (pr1 == _VVnST) color = kMagenta - 5;
        if (pr1 == _WW) color = kMagenta - 5;
        if (pr1 == _WZ) color = kMagenta - 2;
        if (pr1 == _ZZ) color = kMagenta - 6;
        if (pr1 == _tbarW) color = kGreen - 2;
        if (pr1 == _tW) color = kGreen + 2;
        if (pr1 == _ttbar || pr1 == _ttbar_700to1000 || pr1 == _ttbar_1000toInf) color = kCyan + 2;

        h_PFiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_TRKiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_TRKiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_TRKiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_TRKiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_pT_barrel_MC[pr1]->SetFillColor(color);
        h_pT_endcap_MC[pr1]->SetFillColor(color);
        h_eta_MC[pr1]->SetFillColor(color);
        h_nVTX_MC[pr1]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_TRKiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_TRKiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_TRKiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_TRKiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_pT_barrel_MC[pr1]->SetLineColor(color);
        h_pT_endcap_MC[pr1]->SetLineColor(color);
        h_eta_MC[pr1]->SetLineColor(color);
        h_nVTX_MC[pr1]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_TRKiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_TRKiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_TRKiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_TRKiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_pT_barrel_MC[pr1]->SetDirectory(0);
        h_pT_endcap_MC[pr1]->SetDirectory(0);
        h_eta_MC[pr1]->SetDirectory(0);
        h_nVTX_MC[pr1]->SetDirectory(0);

//        if (pr1==_WJets)
//        {
//            h_PFiso_barrel_MC_deno[pr1]->Scale(0.916710636); // 162787688/177578051
//            h_PFiso_endcap_MC_deno[pr1]->Scale(0.916710636);
//            h_TRKiso_barrel_MC_deno[pr1]->Scale(0.916710636);
//            h_TRKiso_endcap_MC_deno[pr1]->Scale(0.916710636);
//            h_PFiso_barrel_MC_nume[pr1]->Scale(0.916710636);
//            h_PFiso_endcap_MC_nume[pr1]->Scale(0.916710636);
//            h_TRKiso_barrel_MC_nume[pr1]->Scale(0.916710636);
//            h_TRKiso_endcap_MC_nume[pr1]->Scale(0.916710636);
//            h_pT_barrel_MC[pr1]->Scale(0.916710636);
//            h_pT_endcap_MC[pr1]->Scale(0.916710636);
//            h_eta_MC[pr1]->Scale(0.916710636);
//            h_nVTX_MC[pr1]->Scale(0.916710636);
//        }

        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr1]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr1]);
        s_TRKiso_barrel_deno->Add(h_TRKiso_barrel_MC_deno[pr1]);
        s_TRKiso_endcap_deno->Add(h_TRKiso_endcap_MC_deno[pr1]);
        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr1]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr1]);
        s_TRKiso_barrel_nume->Add(h_TRKiso_barrel_MC_nume[pr1]);
        s_TRKiso_endcap_nume->Add(h_TRKiso_endcap_MC_nume[pr1]);
        s_pT_barrel->Add(h_pT_barrel_MC[pr1]);
        s_pT_endcap->Add(h_pT_endcap_MC[pr1]);
        s_eta->Add(h_eta_MC[pr1]);
        s_nVTX->Add(h_nVTX_MC[pr1]);

        file->Close();

        if (pr1 == _WW) {pr1 = _WZ; continue;}
        if (pr1 == _WZ) {pr1 = _ZZ; continue;}
        if (pr1 == _ZZ) {pr1 = _tbarW; continue;}
        if (pr1 == _tbarW) {pr1 = _tW; continue;}
        if (pr1 == _tW) {pr1 = _ttbar; continue;}
        if (pr1 == _ttbar) {pr1 = _ttbar_700to1000; continue;}
        if (pr1 == _ttbar_700to1000) {pr1 = _ttbar_1000toInf; continue;}
        if (pr1 == _ttbar_1000toInf) {pr1 = _WJets; continue;}
        if (pr1 == _WJets) {pr1 = _WJets_ext2v5; continue;}
        if (pr1 == _WJets_ext2v5) {stop = 1;}
    }

    // Drell-Yan
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr]);
        file->GetObject("h_TRKiso_barrel_deno", h_TRKiso_barrel_MC_deno[pr]);
        file->GetObject("h_TRKiso_endcap_deno", h_TRKiso_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr]);
        file->GetObject("h_TRKiso_barrel_nume", h_TRKiso_barrel_MC_nume[pr]);
        file->GetObject("h_TRKiso_endcap_nume", h_TRKiso_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
//        removeNegativeBins(h_PFiso_barrel_MC_deno[pr]);
//        removeNegativeBins(h_PFiso_endcap_MC_deno[pr]);
//        removeNegativeBins(h_TRKiso_barrel_MC_deno[pr]);
//        removeNegativeBins(h_TRKiso_endcap_MC_deno[pr]);
//        removeNegativeBins(h_PFiso_barrel_MC_nume[pr]);
//        removeNegativeBins(h_PFiso_endcap_MC_nume[pr]);
//        removeNegativeBins(h_TRKiso_barrel_MC_nume[pr]);
//        removeNegativeBins(h_TRKiso_endcap_MC_nume[pr]);
//        removeNegativeBins(h_pT_barrel_MC[pr]);
//        removeNegativeBins(h_pT_endcap_MC[pr]);
//        removeNegativeBins(h_eta_MC[pr]);
//        removeNegativeBins(h_nVTX_MC[pr]);

        if (pr == _DY_10to50)
        {
            h_PFiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_barrel_MC_deno[pr]->Clone("h_PFiso_barrel_deno_DY")));
            h_PFiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_endcap_MC_deno[pr]->Clone("h_PFiso_endcap_deno_DY")));
            h_TRKiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_TRKiso_barrel_MC_deno[pr]->Clone("h_TRKiso_barrel_deno_DY")));
            h_TRKiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_TRKiso_endcap_MC_deno[pr]->Clone("h_TRKiso_endcap_deno_DY")));
            h_PFiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_barrel_MC_nume[pr]->Clone("h_PFiso_barrel_nume_DY")));
            h_PFiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_endcap_MC_nume[pr]->Clone("h_PFiso_endcap_nume_DY")));
            h_TRKiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_TRKiso_barrel_MC_nume[pr]->Clone("h_TRKiso_barrel_nume_DY")));
            h_TRKiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_TRKiso_endcap_MC_nume[pr]->Clone("h_TRKiso_endcap_nume_DY")));
            h_pT_barrel_MC[_DY_Full] = ((TH1D*)(h_pT_barrel_MC[pr]->Clone("h_pT_barrel_deno_DY")));
            h_pT_endcap_MC[_DY_Full] = ((TH1D*)(h_pT_endcap_MC[pr]->Clone("h_pT_endcap_deno_DY")));
            h_eta_MC[_DY_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_DY")));
            h_nVTX_MC[_DY_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_DY")));
            h_PFiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_TRKiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_TRKiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_TRKiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_TRKiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC[_DY_Full]->SetDirectory(0);
            h_eta_MC[_DY_Full]->SetDirectory(0);
            h_nVTX_MC[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_barrel_MC_deno[_DY_Full]->Add(h_PFiso_barrel_MC_deno[pr]);
            h_PFiso_endcap_MC_deno[_DY_Full]->Add(h_PFiso_endcap_MC_deno[pr]);
            h_TRKiso_barrel_MC_deno[_DY_Full]->Add(h_TRKiso_barrel_MC_deno[pr]);
            h_TRKiso_endcap_MC_deno[_DY_Full]->Add(h_TRKiso_endcap_MC_deno[pr]);
            h_PFiso_barrel_MC_nume[_DY_Full]->Add(h_PFiso_barrel_MC_nume[pr]);
            h_PFiso_endcap_MC_nume[_DY_Full]->Add(h_PFiso_endcap_MC_nume[pr]);
            h_TRKiso_barrel_MC_nume[_DY_Full]->Add(h_TRKiso_barrel_MC_nume[pr]);
            h_TRKiso_endcap_MC_nume[_DY_Full]->Add(h_TRKiso_endcap_MC_nume[pr]);
            h_pT_barrel_MC[_DY_Full]->Add(h_pT_barrel_MC[pr]);
            h_pT_endcap_MC[_DY_Full]->Add(h_pT_endcap_MC[pr]);
            h_eta_MC[_DY_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_DY_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kOrange - 5;
        h_PFiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_TRKiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_TRKiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_TRKiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_TRKiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC[pr]->SetFillColor(color);
        h_pT_endcap_MC[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_TRKiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_TRKiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_TRKiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_TRKiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC[pr]->SetLineColor(color);
        h_pT_endcap_MC[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_TRKiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_TRKiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_TRKiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_TRKiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC[pr]->SetDirectory(0);
        h_pT_endcap_MC[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr]);
        s_TRKiso_barrel_deno->Add(h_TRKiso_barrel_MC_deno[pr]);
        s_TRKiso_endcap_deno->Add(h_TRKiso_endcap_MC_deno[pr]);
        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr]);
        s_TRKiso_barrel_nume->Add(h_TRKiso_barrel_MC_nume[pr]);
        s_TRKiso_endcap_nume->Add(h_TRKiso_endcap_MC_nume[pr]);
        s_pT_barrel->Add(h_pT_barrel_MC[pr]);
        s_pT_endcap->Add(h_pT_endcap_MC[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);

        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr]);
        file->GetObject("h_TRKiso_barrel_deno", h_TRKiso_barrel_MC_deno[pr]);
        file->GetObject("h_TRKiso_endcap_deno", h_TRKiso_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr]);
        file->GetObject("h_TRKiso_barrel_nume", h_TRKiso_barrel_MC_nume[pr]);
        file->GetObject("h_TRKiso_endcap_nume", h_TRKiso_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_TRKiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_TRKiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_TRKiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_TRKiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC[pr]);
        removeNegativeBins(h_pT_endcap_MC[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_barrel_MC_deno[pr]->Clone("h_PFiso_barrel_deno_QCD")));
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_endcap_MC_deno[pr]->Clone("h_PFiso_endcap_deno_QCD")));
            h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_TRKiso_barrel_MC_deno[pr]->Clone("h_TRKiso_barrel_deno_QCD")));
            h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_TRKiso_endcap_MC_deno[pr]->Clone("h_TRKiso_endcap_deno_QCD")));
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_barrel_MC_nume[pr]->Clone("h_PFiso_barrel_nume_QCD")));
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_endcap_MC_nume[pr]->Clone("h_PFiso_endcap_nume_QCD")));
            h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_TRKiso_barrel_MC_nume[pr]->Clone("h_TRKiso_barrel_nume_QCD")));
            h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_TRKiso_endcap_MC_nume[pr]->Clone("h_TRKiso_endcap_nume_QCD")));
            h_pT_barrel_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC[pr]->Clone("h_pT_barrel_deno_QCD")));
            h_pT_endcap_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC[pr]->Clone("h_pT_endcap_deno_QCD")));
            h_eta_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_QCD")));
            h_nVTX_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_QCD")));
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_eta_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_nVTX_MC[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_PFiso_barrel_MC_deno[pr]);
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_PFiso_endcap_MC_deno[pr]);
            h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_TRKiso_barrel_MC_deno[pr]);
            h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_TRKiso_endcap_MC_deno[pr]);
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_PFiso_barrel_MC_nume[pr]);
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_PFiso_endcap_MC_nume[pr]);
            h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_TRKiso_barrel_MC_nume[pr]);
            h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_TRKiso_endcap_MC_nume[pr]);
            h_pT_barrel_MC[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC[pr]);
            h_pT_endcap_MC[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC[pr]);
            h_eta_MC[_QCDMuEnriched_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_QCDMuEnriched_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kRed + 3;
        h_PFiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_TRKiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_TRKiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_TRKiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_TRKiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC[pr]->SetFillColor(color);
        h_pT_endcap_MC[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_TRKiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_TRKiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_TRKiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_TRKiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC[pr]->SetLineColor(color);
        h_pT_endcap_MC[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_TRKiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_TRKiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_TRKiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_TRKiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC[pr]->SetDirectory(0);
        h_pT_endcap_MC[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr]);
        s_TRKiso_barrel_deno->Add(h_TRKiso_barrel_MC_deno[pr]);
        s_TRKiso_endcap_deno->Add(h_TRKiso_endcap_MC_deno[pr]);
        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr]);
        s_TRKiso_barrel_nume->Add(h_TRKiso_barrel_MC_nume[pr]);
        s_TRKiso_endcap_nume->Add(h_TRKiso_endcap_MC_nume[pr]);
        s_pT_barrel->Add(h_pT_barrel_MC[pr]);
        s_pT_endcap->Add(h_pT_endcap_MC[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);

        file->Close();
    }

     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
             h_pT_barrel_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
             h_pT_endcap_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
                   h_eta_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
             h_pT_barrel_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
             h_pT_endcap_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
                   h_eta_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
    h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
    h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
    h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
    h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
             h_pT_barrel_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
             h_pT_endcap_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
                   h_eta_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        TH1D *h_temp[12];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_data_deno);
            file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_data_deno);
            file->GetObject("h_TRKiso_barrel_deno", h_TRKiso_barrel_data_deno);
            file->GetObject("h_TRKiso_endcap_deno", h_TRKiso_endcap_data_deno);
            file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_data_nume);
            file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_data_nume);
            file->GetObject("h_TRKiso_barrel_nume", h_TRKiso_barrel_data_nume);
            file->GetObject("h_TRKiso_endcap_nume", h_TRKiso_endcap_data_nume);
            file->GetObject("h_pT_barrel_deno", h_pT_barrel_data);
            file->GetObject("h_pT_endcap_deno", h_pT_endcap_data);
            file->GetObject("h_eta_deno", h_eta_data);
            file->GetObject("h_nVTX", h_nVTX_data);
            removeNegativeBins(h_PFiso_barrel_data_deno);
            removeNegativeBins(h_PFiso_endcap_data_deno);
            removeNegativeBins(h_TRKiso_barrel_data_deno);
            removeNegativeBins(h_TRKiso_endcap_data_deno);
            removeNegativeBins(h_PFiso_barrel_data_nume);
            removeNegativeBins(h_PFiso_endcap_data_nume);
            removeNegativeBins(h_TRKiso_barrel_data_nume);
            removeNegativeBins(h_TRKiso_endcap_data_nume);
            removeNegativeBins(h_pT_barrel_data);
            removeNegativeBins(h_pT_endcap_data);
            removeNegativeBins(h_eta_data);
            removeNegativeBins(h_nVTX_data);
        }
        else
        {
            file->GetObject("h_PFiso_barrel_deno", h_temp[0]);
            file->GetObject("h_PFiso_endcap_deno", h_temp[1]);
            file->GetObject("h_TRKiso_barrel_deno", h_temp[2]);
            file->GetObject("h_TRKiso_endcap_deno", h_temp[3]);
            file->GetObject("h_PFiso_barrel_nume", h_temp[4]);
            file->GetObject("h_PFiso_endcap_nume", h_temp[5]);
            file->GetObject("h_TRKiso_barrel_nume", h_temp[6]);
            file->GetObject("h_TRKiso_endcap_nume", h_temp[7]);
            file->GetObject("h_pT_barrel_deno", h_temp[8]);
            file->GetObject("h_pT_endcap_deno", h_temp[9]);
            file->GetObject("h_eta_deno", h_temp[10]);
            file->GetObject("h_nVTX", h_temp[11]);
            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            removeNegativeBins(h_temp[4]);
            removeNegativeBins(h_temp[5]);
            removeNegativeBins(h_temp[6]);
            removeNegativeBins(h_temp[7]);
            removeNegativeBins(h_temp[8]);
            removeNegativeBins(h_temp[9]);
            removeNegativeBins(h_temp[10]);
            removeNegativeBins(h_temp[11]);
            h_PFiso_barrel_data_deno->Add(h_temp[0]);
            h_PFiso_endcap_data_deno->Add(h_temp[1]);
            h_TRKiso_barrel_data_deno->Add(h_temp[2]);
            h_TRKiso_endcap_data_deno->Add(h_temp[3]);
            h_PFiso_barrel_data_nume->Add(h_temp[4]);
            h_PFiso_endcap_data_nume->Add(h_temp[5]);
            h_TRKiso_barrel_data_nume->Add(h_temp[6]);
            h_TRKiso_endcap_data_nume->Add(h_temp[7]);
            h_pT_barrel_data->Add(h_temp[8]);
            h_pT_endcap_data->Add(h_temp[9]);
            h_eta_data->Add(h_temp[10]);
            h_nVTX_data->Add(h_temp[11]);
        }
    }

    h_PFiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_TRKiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_TRKiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_TRKiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_TRKiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data->SetMarkerStyle(kFullDotLarge);
    h_eta_data->SetMarkerStyle(kFullDotLarge);
    h_nVTX_data->SetMarkerStyle(kFullDotLarge);
    h_PFiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_PFiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_TRKiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_TRKiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_PFiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_PFiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_TRKiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_TRKiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_pT_barrel_data->SetMarkerColor(kBlack);
    h_pT_endcap_data->SetMarkerColor(kBlack);
    h_eta_data->SetMarkerColor(kBlack);
    h_nVTX_data->SetMarkerColor(kBlack);
    h_PFiso_barrel_data_deno->SetLineColor(kBlack);
    h_PFiso_endcap_data_deno->SetLineColor(kBlack);
    h_TRKiso_barrel_data_deno->SetLineColor(kBlack);
    h_TRKiso_endcap_data_deno->SetLineColor(kBlack);
    h_PFiso_barrel_data_nume->SetLineColor(kBlack);
    h_PFiso_endcap_data_nume->SetLineColor(kBlack);
    h_TRKiso_barrel_data_nume->SetLineColor(kBlack);
    h_TRKiso_endcap_data_nume->SetLineColor(kBlack);
    h_pT_barrel_data->SetLineColor(kBlack);
    h_pT_endcap_data->SetLineColor(kBlack);
    h_eta_data->SetLineColor(kBlack);
    h_nVTX_data->SetLineColor(kBlack);
    h_PFiso_barrel_data_deno->SetDirectory(0);
    h_PFiso_endcap_data_deno->SetDirectory(0);
    h_TRKiso_barrel_data_deno->SetDirectory(0);
    h_TRKiso_endcap_data_deno->SetDirectory(0);
    h_PFiso_barrel_data_nume->SetDirectory(0);
    h_PFiso_endcap_data_nume->SetDirectory(0);
    h_TRKiso_barrel_data_nume->SetDirectory(0);
    h_TRKiso_endcap_data_nume->SetDirectory(0);
    h_pT_barrel_data->SetDirectory(0);
    h_pT_endcap_data->SetDirectory(0);
    h_eta_data->SetDirectory(0);
    h_nVTX_data->SetDirectory(0);

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_PFiso_barrel_deno = new myRatioPlot_t("RP_PFiso_barrel_deno", s_PFiso_barrel_deno, h_PFiso_barrel_data_deno);
    myRatioPlot_t *RP_PFiso_endcap_deno = new myRatioPlot_t("RP_PFiso_endcap_deno", s_PFiso_endcap_deno, h_PFiso_endcap_data_deno);
    myRatioPlot_t *RP_TRKiso_barrel_deno = new myRatioPlot_t("RP_TRKiso_barrel_deno", s_TRKiso_barrel_deno, h_TRKiso_barrel_data_deno);
    myRatioPlot_t *RP_TRKiso_endcap_deno = new myRatioPlot_t("RP_TRKiso_endcap_deno", s_TRKiso_endcap_deno, h_TRKiso_endcap_data_deno);
    myRatioPlot_t *RP_PFiso_barrel_nume = new myRatioPlot_t("RP_PFiso_barrel_nume", s_PFiso_barrel_nume, h_PFiso_barrel_data_nume);
    myRatioPlot_t *RP_PFiso_endcap_nume = new myRatioPlot_t("RP_PFiso_endcap_nume", s_PFiso_endcap_nume, h_PFiso_endcap_data_nume);
    myRatioPlot_t *RP_TRKiso_barrel_nume = new myRatioPlot_t("RP_TRKiso_barrel_nume", s_TRKiso_barrel_nume, h_TRKiso_barrel_data_nume);
    myRatioPlot_t *RP_TRKiso_endcap_nume = new myRatioPlot_t("RP_TRKiso_endcap_nume", s_TRKiso_endcap_nume, h_TRKiso_endcap_data_nume);
    myRatioPlot_t *RP_pT_barrel = new myRatioPlot_t("RP_pT_barrel", s_pT_barrel, h_pT_barrel_data);
    myRatioPlot_t *RP_pT_endcap = new myRatioPlot_t("RP_pT_endcap", s_pT_endcap, h_pT_endcap_data);
    myRatioPlot_t *RP_eta = new myRatioPlot_t("RP_eta", s_eta, h_eta_data);
    myRatioPlot_t *RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_nVTX_data);

    RP_PFiso_barrel_deno->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_PFiso_endcap_deno->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_TRKiso_barrel_deno->SetPlots("I_{TRK}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_TRKiso_endcap_deno->SetPlots("I_{TRK}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_PFiso_barrel_nume->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})", 0, 0.15);
    RP_PFiso_endcap_nume->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})", 0, 0.15);
    RP_TRKiso_barrel_nume->SetPlots("I_{TRK}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})", 0, 1);
    RP_TRKiso_endcap_nume->SetPlots("I_{TRK}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})", 0, 1);
    RP_pT_barrel->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}) [GeV/c]", 50, 500);
    RP_pT_endcap->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}})", 50, 500);
    RP_eta->SetPlots("#eta (#mu)", -3, 3);
    RP_nVTX->SetPlots("N_{#lower[-0.25]{VTX}} (#mu)", 0, 50);


    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_eta_data, "Data", "lp");
//    legend->AddEntry(h_eta_data, "Matavimas", "lp");
    legend->AddEntry(h_eta_MC[_DY_50to100], "DY", "f");
    legend->AddEntry(h_eta_MC[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_eta_MC[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_eta_MC[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_eta_MC[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_eta_MC[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_eta_MC[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
//    legend->AddEntry(h_eta_MC[_VVnST], "Diboson+#font[12]{#scale[1.1]{tW}}+#font[12]{#scale[1.1]{#bar{t}W}}", "f");
    legend->AddEntry(h_eta_MC[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_eta_MC[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->SetNColumns(2);

    RP_PFiso_barrel_deno->ImportLegend(legend);
    RP_PFiso_endcap_deno->ImportLegend(legend);
    RP_TRKiso_barrel_deno->ImportLegend(legend);
    RP_TRKiso_endcap_deno->ImportLegend(legend);
    RP_PFiso_barrel_nume->ImportLegend(legend);
    RP_PFiso_endcap_nume->ImportLegend(legend);
    RP_TRKiso_barrel_nume->ImportLegend(legend);
    RP_TRKiso_endcap_nume->ImportLegend(legend);
    RP_pT_barrel->ImportLegend(legend);
    RP_pT_endcap->ImportLegend(legend);
    RP_eta->ImportLegend(legend);
    RP_nVTX->ImportLegend(legend);

    RP_PFiso_barrel_deno->Draw(1, 1e8, 0);
    RP_PFiso_barrel_deno->DrawOnTop(h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]);

    RP_PFiso_endcap_deno->Draw(1, 1e8, 0);
    RP_PFiso_endcap_deno->DrawOnTop(h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]);

    RP_TRKiso_barrel_deno->Draw(1, 1e8, 0);
    RP_TRKiso_barrel_deno->DrawOnTop(h_TRKiso_barrel_MC_deno[_QCDMuEnriched_Full]);

    RP_TRKiso_endcap_deno->Draw(1, 1e8, 0);
    RP_TRKiso_endcap_deno->DrawOnTop(h_TRKiso_endcap_MC_deno[_QCDMuEnriched_Full]);

    RP_PFiso_barrel_nume->Draw(1, 1e8, 0);
    RP_PFiso_barrel_nume->DrawOnTop(h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]);

    RP_PFiso_endcap_nume->Draw(1, 1e8, 0);
    RP_PFiso_endcap_nume->DrawOnTop(h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]);

    RP_TRKiso_barrel_nume->Draw(1, 1e8, 0);
    RP_TRKiso_barrel_nume->DrawOnTop(h_TRKiso_barrel_MC_nume[_QCDMuEnriched_Full]);

    RP_TRKiso_endcap_nume->Draw(1, 1e8, 0);
    RP_TRKiso_endcap_nume->DrawOnTop(h_TRKiso_endcap_MC_nume[_QCDMuEnriched_Full]);

    RP_pT_barrel->Draw(1, 1e8, 0);
    RP_pT_barrel->DrawOnTop(h_pT_barrel_MC[_QCDMuEnriched_Full]);

    RP_pT_endcap->Draw(1, 1e8, 0);
    RP_pT_endcap->DrawOnTop(h_pT_endcap_MC[_QCDMuEnriched_Full]);

    RP_eta->Draw(1, 1e8, 0);
    RP_nVTX->Draw(1, 1e8, 0);

    cout << "MC PFiso integral: " << ((TH1D*)(s_PFiso_barrel_deno->GetStack()->Last()))->Integral() +
                                     ((TH1D*)(s_PFiso_endcap_deno->GetStack()->Last()))->Integral() << endl;
    cout << "Data PFiso integral: " << h_PFiso_barrel_data_deno->Integral() + h_PFiso_endcap_data_deno->Integral() << endl;
    cout << "--------\nQCD PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "--------\nWJets PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_WJets]->Integral() << endl;
    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_WJets]->Integral() << endl;
    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_WJets]->Integral() << endl;
    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_WJets]->Integral() << endl;

} // End of Mu_HistDrawer


/// ################################################################################## ///
/// ------------------------------------- TEST --------------------------------------- ///
/// ################################################################################## ///
void Test_HistDrawer(Int_t type)
{
    FileMgr fm;
    THStack *s_mass = new THStack("s_mass", "");
    THStack *s_nVTX = new THStack("s_nVTX", "");

    TH1D *h_mass_MC[_EndOf_Data_Special], *h_nVTX_MC[_EndOf_Data_Special],
         *h_mass_data, *h_nVTX_data;

//----------------------------------- MC bkg -------------------------------------------------------

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else return;
        file->GetObject("h_mass", h_mass_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        removeNegativeBins(h_mass_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);

        // Converting to event density (dividing by bin width)
        for (Int_t i = 1; i<=binnum; i++)
        {
            h_mass_MC[pr]->SetBinContent(i, h_mass_MC[pr]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_mass_MC[pr]->SetBinError(i, h_mass_MC[pr]->GetBinError(i)/(massbins[i]-massbins[i-1]));
        }

        if (pr == _QCDMuEnriched_15to20)
        {
            h_mass_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_mass_MC[pr]->Clone("h_mass_QCD")));
            h_nVTX_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_QCD")));
            h_mass_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_nVTX_MC[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_MC[_QCDMuEnriched_Full]->Add(h_mass_MC[pr]);
            h_nVTX_MC[_QCDMuEnriched_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kRed + 3;
        h_mass_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_mass_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_mass_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_mass->Add(h_mass_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);

        file->Close();
    }

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr1]+"_TEST.root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+"_TEST.root", "READ");
        else return;
        file->GetObject("h_mass", h_mass_MC[pr1]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr1]);

        removeNegativeBins(h_mass_MC[pr1]);
        removeNegativeBins(h_nVTX_MC[pr1]);

        // Converting to event density (dividing by bin width)
        for (Int_t i = 1; i<=binnum; i++)
        {
            h_mass_MC[pr1]->SetBinContent(i, h_mass_MC[pr1]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_mass_MC[pr1]->SetBinError(i, h_mass_MC[pr1]->GetBinError(i)/(massbins[i]-massbins[i-1]));
        }

        Color_t color = kBlack;
        if (pr1 == _WJets || pr1 == _WJets_ext2v5) color = kRed - 2;
        if (pr1 == _VVnST) color = kMagenta - 5;
        if (pr1 == _WW) color = kMagenta - 5;
        if (pr1 == _WZ) color = kMagenta - 2;
        if (pr1 == _ZZ) color = kMagenta - 6;
        if (pr1 == _tbarW) color = kGreen - 2;
        if (pr1 == _tW) color = kGreen + 2;
        if (pr1 == _ttbar || pr1 == _ttbar_700to1000 || pr1 == _ttbar_1000toInf) color = kCyan + 2;

        h_mass_MC[pr1]->SetFillColor(color);
        h_nVTX_MC[pr1]->SetFillColor(color);
        h_mass_MC[pr1]->SetLineColor(color);
        h_nVTX_MC[pr1]->SetLineColor(color);
        h_mass_MC[pr1]->SetDirectory(0);
        h_nVTX_MC[pr1]->SetDirectory(0);

        s_mass->Add(h_mass_MC[pr1]);
        s_nVTX->Add(h_nVTX_MC[pr1]);

        file->Close();

        if (pr1 == _WW) {pr1 = _WZ; continue;}
        if (pr1 == _WZ) {pr1 = _ZZ; continue;}
        if (pr1 == _ZZ) {pr1 = _tbarW; continue;}
        if (pr1 == _tbarW) {pr1 = _tW; continue;}
        if (pr1 == _tW) {pr1 = _ttbar; continue;}
        if (pr1 == _ttbar) {pr1 = _ttbar_700to1000; continue;}
        if (pr1 == _ttbar_700to1000) {pr1 = _ttbar_1000toInf; continue;}
        if (pr1 == _ttbar_1000toInf) {pr1 = _WJets; continue;}
        if (pr1 == _WJets) {pr1 = _WJets_ext2v5; continue;}
        if (pr1 == _WJets_ext2v5) {stop = 1;}
    }

    // Drell-Yan
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else return;
        file->GetObject("h_mass", h_mass_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        removeNegativeBins(h_mass_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);

        // Converting to event density (dividing by bin width)
        for (Int_t i = 1; i<=binnum; i++)
        {
            h_mass_MC[pr]->SetBinContent(i, h_mass_MC[pr]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
            h_mass_MC[pr]->SetBinError(i, h_mass_MC[pr]->GetBinError(i)/(massbins[i]-massbins[i-1]));
        }

        if (pr == _DY_10to50)
        {
            h_mass_MC[_DY_Full] = ((TH1D*)(h_mass_MC[pr]->Clone("h_mass_DY")));
            h_nVTX_MC[_DY_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_DY")));
            h_mass_MC[_DY_Full]->SetDirectory(0);
            h_nVTX_MC[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_MC[_DY_Full]->Add(h_mass_MC[pr]);
            h_nVTX_MC[_DY_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kOrange - 5;
        h_mass_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_mass_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_mass_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_mass->Add(h_mass_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);

        file->Close();
    }

    h_mass_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_nVTX_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
    h_mass_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_nVTX_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
    h_mass_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
    h_nVTX_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+"_TEST.root", "READ");
        else return;
        TH1D *h_temp[2];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_mass", h_mass_data);
            file->GetObject("h_nVTX", h_nVTX_data);
            removeNegativeBins(h_mass_data);
            removeNegativeBins(h_nVTX_data);
            // Converting to event density (dividing by bin width)
            for (Int_t i = 1; i<=binnum; i++)
            {
                h_mass_data->SetBinContent(i, h_mass_data->GetBinContent(i)/(massbins[i]-massbins[i-1]));
                h_mass_data->SetBinError(i, h_mass_data->GetBinError(i)/(massbins[i]-massbins[i-1]));
            }
        }
        else
        {
            file->GetObject("h_mass", h_temp[0]);
            file->GetObject("h_nVTX", h_temp[1]);
            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            // Converting to event density (dividing by bin width)
            for (Int_t i = 1; i<=binnum; i++)
            {
                h_temp[0]->SetBinContent(i, h_temp[0]->GetBinContent(i)/(massbins[i]-massbins[i-1]));
                h_temp[0]->SetBinError(i, h_temp[0]->GetBinError(i)/(massbins[i]-massbins[i-1]));
            }
            h_mass_data->Add(h_temp[0]);
            h_nVTX_data->Add(h_temp[1]);
        }
    }

    h_mass_data->SetMarkerStyle(kFullDotLarge);
    h_nVTX_data->SetMarkerStyle(kFullDotLarge);
    h_mass_data->SetMarkerColor(kBlack);
    h_nVTX_data->SetMarkerColor(kBlack);
    h_mass_data->SetLineColor(kBlack);
    h_nVTX_data->SetLineColor(kBlack);
    h_mass_data->SetDirectory(0);
    h_nVTX_data->SetDirectory(0);

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_mass = new myRatioPlot_t("RP_mass", s_mass, h_mass_data);
    myRatioPlot_t *RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_nVTX_data);

    RP_mass->SetPlots("m_{#mu#mu} [GeV/c^{2}]", -3, 3);
    RP_nVTX->SetPlots("N_{#lower[-0.25]{VTX}} (#mu)", 0, 50);


    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_mass_data, "Data", "lp");
//    legend->AddEntry(h_mass_data, "Matavimas", "lp");
    legend->AddEntry(h_mass_MC[_DY_50to100], "DY", "f");
    legend->AddEntry(h_mass_MC[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_mass_MC[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_mass_MC[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_mass_MC[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_mass_MC[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_mass_MC[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
//    legend->AddEntry(h_mass_MC[_VVnST], "Diboson+#font[12]{#scale[1.1]{tW}}+#font[12]{#scale[1.1]{#bar{t}W}}", "f");
    legend->AddEntry(h_mass_MC[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_mass_MC[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->SetNColumns(2);

    RP_mass->ImportLegend(legend);
    RP_nVTX->ImportLegend(legend);

    RP_mass->Draw(1, 1e6, 1);
    RP_mass->s_stackedProcesses->GetYaxis()->SetTitle("N_{#lower[-0.25]{EVT}}/GeV");
    RP_mass->pad1->cd();
    RP_mass->s_stackedProcesses->Draw("sameaxis");
    RP_mass->legend->Draw();
    RP_mass->pad1->Update();
    RP_mass->canvas->Update();
    RP_nVTX->Draw(1, 1e8, 0);

    cout << "MC mass integral: " << ((TH1D*)(s_mass->GetStack()->Last()))->Integral() << endl;
    cout << "Data mass integral: " << h_mass_data->Integral() << endl;

} // End of Test_HistDrawer


/// ############################################################################# ///
/// ---------------------------- Template fit test ------------------------------ ///
/// ############################################################################# ///
void Fit_HistDrawer()
{
    FileMgr fm;

    TH1D *h_barrel_MC_deno_50to70  [_EndOf_Data_Special],
         *h_barrel_MC_nume_50to70  [_EndOf_Data_Special],
         *h_endcap_MC_deno_50to70  [_EndOf_Data_Special],
         *h_endcap_MC_nume_50to70  [_EndOf_Data_Special],
         *h_barrel_MC_deno_70to100 [_EndOf_Data_Special],
         *h_barrel_MC_nume_70to100 [_EndOf_Data_Special],
         *h_endcap_MC_deno_70to100 [_EndOf_Data_Special],
         *h_endcap_MC_nume_70to100 [_EndOf_Data_Special],
         *h_barrel_MC_deno_100to200[_EndOf_Data_Special],
         *h_barrel_MC_nume_100to200[_EndOf_Data_Special],
         *h_endcap_MC_deno_100to200[_EndOf_Data_Special],
         *h_endcap_MC_nume_100to200[_EndOf_Data_Special],
         *h_barrel_MC_deno_200to500[_EndOf_Data_Special],
         *h_barrel_MC_nume_200to500[_EndOf_Data_Special],
         *h_endcap_MC_deno_200to500[_EndOf_Data_Special],
         *h_endcap_MC_nume_200to500[_EndOf_Data_Special],
         *h_barrel_data_deno_50to70,
         *h_barrel_data_nume_50to70,
         *h_endcap_data_deno_50to70,
         *h_endcap_data_nume_50to70,
         *h_barrel_data_deno_70to100,
         *h_barrel_data_nume_70to100,
         *h_endcap_data_deno_70to100,
         *h_endcap_data_nume_70to100,
         *h_barrel_data_deno_100to200,
         *h_barrel_data_nume_100to200,
         *h_endcap_data_deno_100to200,
         *h_endcap_data_nume_100to200,
         *h_barrel_data_deno_200to500,
         *h_barrel_data_nume_200to500,
         *h_endcap_data_deno_200to500,
         *h_endcap_data_nume_200to500;

    THStack *s_barrel_deno_50to70,
            *s_barrel_nume_50to70,
            *s_endcap_deno_50to70,
            *s_endcap_nume_50to70,
            *s_barrel_deno_70to100,
            *s_barrel_nume_70to100,
            *s_endcap_deno_70to100,
            *s_endcap_nume_70to100,
            *s_barrel_deno_100to200,
            *s_barrel_nume_100to200,
            *s_endcap_deno_100to200,
            *s_endcap_nume_100to200,
            *s_barrel_deno_200to500,
            *s_barrel_nume_200to500,
            *s_endcap_deno_200to500,
            *s_endcap_nume_200to500;

// ############################# SETUP ################################# //
//----------------------------- MC bkg ------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+".root", "READ");
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr1]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr1]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr1]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr1]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr1]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr1]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr1]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr1]);
        file->GetObject("h_PFiso_barrel_deno_100to200", h_barrel_MC_deno_100to200[pr1]);
        file->GetObject("h_PFiso_endcap_deno_100to200", h_endcap_MC_deno_100to200[pr1]);
        file->GetObject("h_PFiso_barrel_nume_100to200", h_barrel_MC_nume_100to200[pr1]);
        file->GetObject("h_PFiso_endcap_nume_100to200", h_endcap_MC_nume_100to200[pr1]);
        file->GetObject("h_PFiso_barrel_deno_200to500", h_barrel_MC_deno_200to500[pr1]);
        file->GetObject("h_PFiso_endcap_deno_200to500", h_endcap_MC_deno_200to500[pr1]);
        file->GetObject("h_PFiso_barrel_nume_200to500", h_barrel_MC_nume_200to500[pr1]);
        file->GetObject("h_PFiso_endcap_nume_200to500", h_endcap_MC_nume_200to500[pr1]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_deno_100to200[pr1]);
        removeNegativeBins(h_endcap_MC_deno_100to200[pr1]);
        removeNegativeBins(h_barrel_MC_nume_100to200[pr1]);
        removeNegativeBins(h_endcap_MC_nume_100to200[pr1]);
        removeNegativeBins(h_barrel_MC_deno_200to500[pr1]);
        removeNegativeBins(h_endcap_MC_deno_200to500[pr1]);
        removeNegativeBins(h_barrel_MC_nume_200to500[pr1]);
        removeNegativeBins(h_endcap_MC_nume_200to500[pr1]);

        h_barrel_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_deno_100to200[pr1]->SetDirectory(0);
        h_endcap_MC_deno_100to200[pr1]->SetDirectory(0);
        h_barrel_MC_nume_100to200[pr1]->SetDirectory(0);
        h_endcap_MC_nume_100to200[pr1]->SetDirectory(0);
        h_barrel_MC_deno_200to500[pr1]->SetDirectory(0);
        h_endcap_MC_deno_200to500[pr1]->SetDirectory(0);
        h_barrel_MC_nume_200to500[pr1]->SetDirectory(0);
        h_endcap_MC_nume_200to500[pr1]->SetDirectory(0);

        Color_t color = kBlack;
        if (pr1 == _WJets || pr1 == _WJets_ext2v5) color = kRed - 2;
        if (pr1 == _VVnST) color = kMagenta - 5;
        if (pr1 == _WW) color = kMagenta - 5;
        if (pr1 == _WZ) color = kMagenta - 2;
        if (pr1 == _ZZ) color = kMagenta - 6;
        if (pr1 == _tbarW) color = kGreen - 2;
        if (pr1 == _tW) color = kGreen + 2;
        if (pr1 == _ttbar || pr1 == _ttbar_700to1000 || pr1 == _ttbar_1000toInf) color = kCyan + 2;

        h_barrel_MC_deno_50to70  [pr1]->SetFillColor(color);
        h_endcap_MC_deno_50to70  [pr1]->SetFillColor(color);
        h_barrel_MC_nume_50to70  [pr1]->SetFillColor(color);
        h_endcap_MC_nume_50to70  [pr1]->SetFillColor(color);
        h_barrel_MC_deno_70to100 [pr1]->SetFillColor(color);
        h_endcap_MC_deno_70to100 [pr1]->SetFillColor(color);
        h_barrel_MC_nume_70to100 [pr1]->SetFillColor(color);
        h_endcap_MC_nume_70to100 [pr1]->SetFillColor(color);
        h_barrel_MC_deno_100to200[pr1]->SetFillColor(color);
        h_endcap_MC_deno_100to200[pr1]->SetFillColor(color);
        h_barrel_MC_nume_100to200[pr1]->SetFillColor(color);
        h_endcap_MC_nume_100to200[pr1]->SetFillColor(color);
        h_barrel_MC_deno_200to500[pr1]->SetFillColor(color);
        h_endcap_MC_deno_200to500[pr1]->SetFillColor(color);
        h_barrel_MC_nume_200to500[pr1]->SetFillColor(color);
        h_endcap_MC_nume_200to500[pr1]->SetFillColor(color);

        h_barrel_MC_deno_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_deno_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_deno_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_deno_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_deno_100to200[pr1]->SetLineColor(color);
        h_endcap_MC_deno_100to200[pr1]->SetLineColor(color);
        h_barrel_MC_nume_100to200[pr1]->SetLineColor(color);
        h_endcap_MC_nume_100to200[pr1]->SetLineColor(color);
        h_barrel_MC_deno_200to500[pr1]->SetLineColor(color);
        h_endcap_MC_deno_200to500[pr1]->SetLineColor(color);
        h_barrel_MC_nume_200to500[pr1]->SetLineColor(color);
        h_endcap_MC_nume_200to500[pr1]->SetLineColor(color);

        file->Close();

        if (pr1 == _WW) {pr1 = _WZ; continue;}
        if (pr1 == _WZ) {pr1 = _ZZ; continue;}
        if (pr1 == _ZZ) {pr1 = _tbarW; continue;}
        if (pr1 == _tbarW) {pr1 = _tW; continue;}
        if (pr1 == _tW) {pr1 = _ttbar; continue;}
        if (pr1 == _ttbar) {pr1 = _ttbar_700to1000; continue;}
        if (pr1 == _ttbar_700to1000) {pr1 = _ttbar_1000toInf; continue;}
        if (pr1 == _ttbar_1000toInf) {pr1 = _WJets; continue;}
        if (pr1 == _WJets) {pr1 = _WJets_ext2v5; continue;}
        if (pr1 == _WJets_ext2v5) {stop = 1;}
    }
    h_barrel_MC_deno_50to70  [_ttbar]->Add(h_barrel_MC_deno_50to70  [_ttbar_700to1000]);
    h_endcap_MC_deno_50to70  [_ttbar]->Add(h_endcap_MC_deno_50to70  [_ttbar_700to1000]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_700to1000]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_700to1000]);
    h_barrel_MC_deno_70to100 [_ttbar]->Add(h_barrel_MC_deno_70to100 [_ttbar_700to1000]);
    h_endcap_MC_deno_70to100 [_ttbar]->Add(h_endcap_MC_deno_70to100 [_ttbar_700to1000]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_700to1000]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_700to1000]);
    h_barrel_MC_deno_100to200[_ttbar]->Add(h_barrel_MC_deno_100to200[_ttbar_700to1000]);
    h_endcap_MC_deno_100to200[_ttbar]->Add(h_endcap_MC_deno_100to200[_ttbar_700to1000]);
    h_barrel_MC_nume_100to200[_ttbar]->Add(h_barrel_MC_nume_100to200[_ttbar_700to1000]);
    h_endcap_MC_nume_100to200[_ttbar]->Add(h_endcap_MC_nume_100to200[_ttbar_700to1000]);
    h_barrel_MC_deno_200to500[_ttbar]->Add(h_barrel_MC_deno_200to500[_ttbar_700to1000]);
    h_endcap_MC_deno_200to500[_ttbar]->Add(h_endcap_MC_deno_200to500[_ttbar_700to1000]);
    h_barrel_MC_nume_200to500[_ttbar]->Add(h_barrel_MC_nume_200to500[_ttbar_700to1000]);
    h_endcap_MC_nume_200to500[_ttbar]->Add(h_endcap_MC_nume_200to500[_ttbar_700to1000]);

    h_barrel_MC_deno_50to70  [_ttbar]->Add(h_barrel_MC_deno_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_deno_50to70  [_ttbar]->Add(h_endcap_MC_deno_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_deno_70to100 [_ttbar]->Add(h_barrel_MC_deno_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_deno_70to100 [_ttbar]->Add(h_endcap_MC_deno_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_deno_100to200[_ttbar]->Add(h_barrel_MC_deno_100to200[_ttbar_1000toInf]);
    h_endcap_MC_deno_100to200[_ttbar]->Add(h_endcap_MC_deno_100to200[_ttbar_1000toInf]);
    h_barrel_MC_nume_100to200[_ttbar]->Add(h_barrel_MC_nume_100to200[_ttbar_1000toInf]);
    h_endcap_MC_nume_100to200[_ttbar]->Add(h_endcap_MC_nume_100to200[_ttbar_1000toInf]);
    h_barrel_MC_deno_200to500[_ttbar]->Add(h_barrel_MC_deno_200to500[_ttbar_1000toInf]);
    h_endcap_MC_deno_200to500[_ttbar]->Add(h_endcap_MC_deno_200to500[_ttbar_1000toInf]);
    h_barrel_MC_nume_200to500[_ttbar]->Add(h_barrel_MC_nume_200to500[_ttbar_1000toInf]);
    h_endcap_MC_nume_200to500[_ttbar]->Add(h_endcap_MC_nume_200to500[_ttbar_1000toInf]);

    h_barrel_MC_deno_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_deno_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_deno_100to200[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_100to200[_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_100to200[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_100to200[_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_deno_200to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_200to500[_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_200to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_200to500[_ttbar]->SetFillColor(kCyan+2);

    h_barrel_MC_deno_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_deno_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_deno_100to200[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_100to200[_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_100to200[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_100to200[_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_deno_200to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_200to500[_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_200to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_200to500[_ttbar]->SetLineColor(kCyan+2);


    h_barrel_MC_deno_50to70  [_WJets]->Add(h_barrel_MC_deno_50to70  [_WJets_ext2v5]);
    h_endcap_MC_deno_50to70  [_WJets]->Add(h_endcap_MC_deno_50to70  [_WJets_ext2v5]);
    h_barrel_MC_nume_50to70  [_WJets]->Add(h_barrel_MC_nume_50to70  [_WJets_ext2v5]);
    h_endcap_MC_nume_50to70  [_WJets]->Add(h_endcap_MC_nume_50to70  [_WJets_ext2v5]);
    h_barrel_MC_deno_70to100 [_WJets]->Add(h_barrel_MC_deno_70to100 [_WJets_ext2v5]);
    h_endcap_MC_deno_70to100 [_WJets]->Add(h_endcap_MC_deno_70to100 [_WJets_ext2v5]);
    h_barrel_MC_nume_70to100 [_WJets]->Add(h_barrel_MC_nume_70to100 [_WJets_ext2v5]);
    h_endcap_MC_nume_70to100 [_WJets]->Add(h_endcap_MC_nume_70to100 [_WJets_ext2v5]);
    h_barrel_MC_deno_100to200[_WJets]->Add(h_barrel_MC_deno_100to200[_WJets_ext2v5]);
    h_endcap_MC_deno_100to200[_WJets]->Add(h_endcap_MC_deno_100to200[_WJets_ext2v5]);
    h_barrel_MC_nume_100to200[_WJets]->Add(h_barrel_MC_nume_100to200[_WJets_ext2v5]);
    h_endcap_MC_nume_100to200[_WJets]->Add(h_endcap_MC_nume_100to200[_WJets_ext2v5]);
    h_barrel_MC_deno_200to500[_WJets]->Add(h_barrel_MC_deno_200to500[_WJets_ext2v5]);
    h_endcap_MC_deno_200to500[_WJets]->Add(h_endcap_MC_deno_200to500[_WJets_ext2v5]);
    h_barrel_MC_nume_200to500[_WJets]->Add(h_barrel_MC_nume_200to500[_WJets_ext2v5]);
    h_endcap_MC_nume_200to500[_WJets]->Add(h_endcap_MC_nume_200to500[_WJets_ext2v5]);

    h_barrel_MC_deno_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_deno_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_deno_100to200[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_100to200[_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_100to200[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_100to200[_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_deno_200to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_200to500[_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_200to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_200to500[_WJets]->SetFillColor(kRed-2);

    h_barrel_MC_deno_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_deno_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_deno_100to200[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_100to200[_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_100to200[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_100to200[_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_deno_200to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_200to500[_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_200to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_200to500[_WJets]->SetLineColor(kRed-2);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_deno_100to200", h_barrel_MC_deno_100to200[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to200", h_endcap_MC_deno_100to200[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to200", h_barrel_MC_nume_100to200[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to200", h_endcap_MC_nume_100to200[pr]);
        file->GetObject("h_PFiso_barrel_deno_200to500", h_barrel_MC_deno_200to500[pr]);
        file->GetObject("h_PFiso_endcap_deno_200to500", h_endcap_MC_deno_200to500[pr]);
        file->GetObject("h_PFiso_barrel_nume_200to500", h_barrel_MC_nume_200to500[pr]);
        file->GetObject("h_PFiso_endcap_nume_200to500", h_endcap_MC_nume_200to500[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno_100to200[pr]);
        removeNegativeBins(h_endcap_MC_deno_100to200[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to200[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to200[pr]);
        removeNegativeBins(h_barrel_MC_deno_200to500[pr]);
        removeNegativeBins(h_endcap_MC_deno_200to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_200to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_200to500[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_deno_100to200[pr]->SetDirectory(0);
        h_endcap_MC_deno_100to200[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to200[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to200[pr]->SetDirectory(0);
        h_barrel_MC_deno_200to500[pr]->SetDirectory(0);
        h_endcap_MC_deno_200to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_200to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_200to500[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_barrel_MC_deno_50to70  [_DY_Full] = ((TH1D*)(h_barrel_MC_deno_50to70  [pr]->Clone("h_barrel_MC_deno_DY_50to70")));
            h_endcap_MC_deno_50to70  [_DY_Full] = ((TH1D*)(h_endcap_MC_deno_50to70  [pr]->Clone("h_endcap_MC_deno_DY_50to70")));
            h_barrel_MC_nume_50to70  [_DY_Full] = ((TH1D*)(h_barrel_MC_nume_50to70  [pr]->Clone("h_barrel_MC_nume_DY_50to70")));
            h_endcap_MC_nume_50to70  [_DY_Full] = ((TH1D*)(h_endcap_MC_nume_50to70  [pr]->Clone("h_endcap_MC_nume_DY_50to70")));
            h_barrel_MC_deno_70to100 [_DY_Full] = ((TH1D*)(h_barrel_MC_deno_70to100 [pr]->Clone("h_barrel_MC_deno_DY_70to100")));
            h_endcap_MC_deno_70to100 [_DY_Full] = ((TH1D*)(h_endcap_MC_deno_70to100 [pr]->Clone("h_endcap_MC_deno_DY_70to100")));
            h_barrel_MC_nume_70to100 [_DY_Full] = ((TH1D*)(h_barrel_MC_nume_70to100 [pr]->Clone("h_barrel_MC_nume_DY_70to100")));
            h_endcap_MC_nume_70to100 [_DY_Full] = ((TH1D*)(h_endcap_MC_nume_70to100 [pr]->Clone("h_endcap_MC_nume_DY_70to100")));
            h_barrel_MC_deno_100to200[_DY_Full] = ((TH1D*)(h_barrel_MC_deno_100to200[pr]->Clone("h_barrel_MC_deno_DY_100to200")));
            h_endcap_MC_deno_100to200[_DY_Full] = ((TH1D*)(h_endcap_MC_deno_100to200[pr]->Clone("h_endcap_MC_deno_DY_100to200")));
            h_barrel_MC_nume_100to200[_DY_Full] = ((TH1D*)(h_barrel_MC_nume_100to200[pr]->Clone("h_barrel_MC_nume_DY_100to200")));
            h_endcap_MC_nume_100to200[_DY_Full] = ((TH1D*)(h_endcap_MC_nume_100to200[pr]->Clone("h_endcap_MC_nume_DY_100to200")));
            h_barrel_MC_deno_200to500[_DY_Full] = ((TH1D*)(h_barrel_MC_deno_200to500[pr]->Clone("h_barrel_MC_deno_DY_200to500")));
            h_endcap_MC_deno_200to500[_DY_Full] = ((TH1D*)(h_endcap_MC_deno_200to500[pr]->Clone("h_endcap_MC_deno_DY_200to500")));
            h_barrel_MC_nume_200to500[_DY_Full] = ((TH1D*)(h_barrel_MC_nume_200to500[pr]->Clone("h_barrel_MC_nume_DY_200to500")));
            h_endcap_MC_nume_200to500[_DY_Full] = ((TH1D*)(h_endcap_MC_nume_200to500[pr]->Clone("h_endcap_MC_nume_DY_200to500")));

            h_barrel_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_100to200[_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_100to200[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to200[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to200[_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_200to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_200to500[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_200to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_200to500[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_deno_50to70  [_DY_Full]->Add(h_barrel_MC_deno_50to70  [pr]);
            h_endcap_MC_deno_50to70  [_DY_Full]->Add(h_endcap_MC_deno_50to70  [pr]);
            h_barrel_MC_nume_50to70  [_DY_Full]->Add(h_barrel_MC_nume_50to70  [pr]);
            h_endcap_MC_nume_50to70  [_DY_Full]->Add(h_endcap_MC_nume_50to70  [pr]);
            h_barrel_MC_deno_70to100 [_DY_Full]->Add(h_barrel_MC_deno_70to100 [pr]);
            h_endcap_MC_deno_70to100 [_DY_Full]->Add(h_endcap_MC_deno_70to100 [pr]);
            h_barrel_MC_nume_70to100 [_DY_Full]->Add(h_barrel_MC_nume_70to100 [pr]);
            h_endcap_MC_nume_70to100 [_DY_Full]->Add(h_endcap_MC_nume_70to100 [pr]);
            h_barrel_MC_deno_100to200[_DY_Full]->Add(h_barrel_MC_deno_100to200[pr]);
            h_endcap_MC_deno_100to200[_DY_Full]->Add(h_endcap_MC_deno_100to200[pr]);
            h_barrel_MC_nume_100to200[_DY_Full]->Add(h_barrel_MC_nume_100to200[pr]);
            h_endcap_MC_nume_100to200[_DY_Full]->Add(h_endcap_MC_nume_100to200[pr]);
            h_barrel_MC_deno_200to500[_DY_Full]->Add(h_barrel_MC_deno_200to500[pr]);
            h_endcap_MC_deno_200to500[_DY_Full]->Add(h_endcap_MC_deno_200to500[pr]);
            h_barrel_MC_nume_200to500[_DY_Full]->Add(h_barrel_MC_nume_200to500[pr]);
            h_endcap_MC_nume_200to500[_DY_Full]->Add(h_endcap_MC_nume_200to500[pr]);
        }
        file->Close();
    }

    h_barrel_MC_deno_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_deno_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_deno_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_deno_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_deno_100to200[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_deno_100to200[_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_100to200[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_100to200[_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_deno_200to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_deno_200to500[_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_200to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_200to500[_DY_Full]->SetFillColor(kOrange-5);

    h_barrel_MC_deno_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_deno_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_deno_100to200[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_100to200[_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_100to200[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_100to200[_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_deno_200to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_200to500[_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_200to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_200to500[_DY_Full]->SetLineColor(kOrange-5);

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_barrel_deno_100to200", h_barrel_MC_deno_100to200[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to200", h_endcap_MC_deno_100to200[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to200", h_barrel_MC_nume_100to200[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to200", h_endcap_MC_nume_100to200[pr]);
        file->GetObject("h_PFiso_barrel_deno_200to500", h_barrel_MC_deno_200to500[pr]);
        file->GetObject("h_PFiso_endcap_deno_200to500", h_endcap_MC_deno_200to500[pr]);
        file->GetObject("h_PFiso_barrel_nume_200to500", h_barrel_MC_nume_200to500[pr]);
        file->GetObject("h_PFiso_endcap_nume_200to500", h_endcap_MC_nume_200to500[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno_100to200[pr]);
        removeNegativeBins(h_endcap_MC_deno_100to200[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to200[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to200[pr]);
        removeNegativeBins(h_barrel_MC_deno_200to500[pr]);
        removeNegativeBins(h_endcap_MC_deno_200to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_200to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_200to500[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_deno_100to200[pr]->SetDirectory(0);
        h_endcap_MC_deno_100to200[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to200[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to200[pr]->SetDirectory(0);
        h_barrel_MC_deno_200to500[pr]->SetDirectory(0);
        h_endcap_MC_deno_200to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_200to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_200to500[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno_50to70  [pr]->Clone("h_barrel_MC_deno_QCD_50to70")));
            h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno_50to70  [pr]->Clone("h_endcap_MC_deno_QCD_50to70")));
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_50to70  [pr]->Clone("h_barrel_MC_nume_QCD_50to70")));
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_50to70  [pr]->Clone("h_endcap_MC_nume_QCD_50to70")));
            h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno_70to100 [pr]->Clone("h_barrel_MC_deno_QCD_70to100")));
            h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno_70to100 [pr]->Clone("h_endcap_MC_deno_QCD_70to100")));
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_70to100 [pr]->Clone("h_barrel_MC_nume_QCD_70to100")));
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_70to100 [pr]->Clone("h_endcap_MC_nume_QCD_70to100")));
            h_barrel_MC_deno_100to200[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno_100to200[pr]->Clone("h_barrel_MC_deno_QCD_100to200")));
            h_endcap_MC_deno_100to200[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno_100to200[pr]->Clone("h_endcap_MC_deno_QCD_100to200")));
            h_barrel_MC_nume_100to200[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_100to200[pr]->Clone("h_barrel_MC_nume_QCD_100to200")));
            h_endcap_MC_nume_100to200[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_100to200[pr]->Clone("h_endcap_MC_nume_QCD_100to200")));
            h_barrel_MC_deno_200to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno_200to500[pr]->Clone("h_barrel_MC_deno_QCD_200to500")));
            h_endcap_MC_deno_200to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno_200to500[pr]->Clone("h_endcap_MC_deno_QCD_200to500")));
            h_barrel_MC_nume_200to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_200to500[pr]->Clone("h_barrel_MC_nume_QCD_200to500")));
            h_endcap_MC_nume_200to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_200to500[pr]->Clone("h_endcap_MC_nume_QCD_200to500")));

            h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_100to200[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_100to200[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to200[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to200[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_200to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_200to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_200to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_200to500[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->Add(h_barrel_MC_deno_50to70 [pr]);
            h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->Add(h_endcap_MC_deno_50to70 [pr]);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_50to70 [pr]);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_50to70 [pr]);
            h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->Add(h_barrel_MC_deno_70to100[pr]);
            h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->Add(h_endcap_MC_deno_70to100[pr]);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_70to100[pr]);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_70to100[pr]);
            h_barrel_MC_deno_100to200[_QCDMuEnriched_Full]->Add(h_barrel_MC_deno_100to200[pr]);
            h_endcap_MC_deno_100to200[_QCDMuEnriched_Full]->Add(h_endcap_MC_deno_100to200[pr]);
            h_barrel_MC_nume_100to200[_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_100to200[pr]);
            h_endcap_MC_nume_100to200[_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_100to200[pr]);
            h_barrel_MC_deno_200to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_deno_200to500[pr]);
            h_endcap_MC_deno_200to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_deno_200to500[pr]);
            h_barrel_MC_nume_200to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_200to500[pr]);
            h_endcap_MC_nume_200to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_200to500[pr]);
        }
        file->Close();
    }

    h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_deno_100to200[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_deno_100to200[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_100to200[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_100to200[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_deno_200to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_deno_200to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_200to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_200to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);

    h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_deno_100to200[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_100to200[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_100to200[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_100to200[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_deno_200to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_200to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_200to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_200to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        TH1D *h_temp[16];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_data_deno_50to70 );
            file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_data_deno_50to70 );
            file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_data_nume_50to70 );
            file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_data_nume_50to70 );
            file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_data_deno_70to100);
            file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_data_deno_70to100);
            file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_data_nume_70to100);
            file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_data_nume_70to100);
            file->GetObject("h_PFiso_barrel_deno_100to200", h_barrel_data_deno_100to200);
            file->GetObject("h_PFiso_endcap_deno_100to200", h_endcap_data_deno_100to200);
            file->GetObject("h_PFiso_barrel_nume_100to200", h_barrel_data_nume_100to200);
            file->GetObject("h_PFiso_endcap_nume_100to200", h_endcap_data_nume_100to200);
            file->GetObject("h_PFiso_barrel_deno_200to500", h_barrel_data_deno_200to500);
            file->GetObject("h_PFiso_endcap_deno_200to500", h_endcap_data_deno_200to500);
            file->GetObject("h_PFiso_barrel_nume_200to500", h_barrel_data_nume_200to500);
            file->GetObject("h_PFiso_endcap_nume_200to500", h_endcap_data_nume_200to500);

            removeNegativeBins(h_barrel_data_deno_50to70);
            removeNegativeBins(h_endcap_data_deno_50to70);
            removeNegativeBins(h_barrel_data_nume_50to70);
            removeNegativeBins(h_endcap_data_nume_50to70);
            removeNegativeBins(h_barrel_data_deno_70to100);
            removeNegativeBins(h_endcap_data_deno_70to100);
            removeNegativeBins(h_barrel_data_nume_70to100);
            removeNegativeBins(h_endcap_data_nume_70to100);
            removeNegativeBins(h_barrel_data_deno_100to200);
            removeNegativeBins(h_endcap_data_deno_100to200);
            removeNegativeBins(h_barrel_data_nume_100to200);
            removeNegativeBins(h_endcap_data_nume_100to200);
            removeNegativeBins(h_barrel_data_deno_200to500);
            removeNegativeBins(h_endcap_data_deno_200to500);
            removeNegativeBins(h_barrel_data_nume_200to500);
            removeNegativeBins(h_endcap_data_nume_200to500);
        }
        else
        {
            file->GetObject("h_PFiso_barrel_deno_50to70",   h_temp[0]);
            file->GetObject("h_PFiso_endcap_deno_50to70",   h_temp[1]);
            file->GetObject("h_PFiso_barrel_nume_50to70",   h_temp[2]);
            file->GetObject("h_PFiso_endcap_nume_50to70",   h_temp[3]);
            file->GetObject("h_PFiso_barrel_deno_70to100",  h_temp[4]);
            file->GetObject("h_PFiso_endcap_deno_70to100",  h_temp[5]);
            file->GetObject("h_PFiso_barrel_nume_70to100",  h_temp[6]);
            file->GetObject("h_PFiso_endcap_nume_70to100",  h_temp[7]);
            file->GetObject("h_PFiso_barrel_deno_100to200", h_temp[8]);
            file->GetObject("h_PFiso_endcap_deno_100to200", h_temp[9]);
            file->GetObject("h_PFiso_barrel_nume_100to200", h_temp[10]);
            file->GetObject("h_PFiso_endcap_nume_100to200", h_temp[11]);
            file->GetObject("h_PFiso_barrel_deno_200to500", h_temp[12]);
            file->GetObject("h_PFiso_endcap_deno_200to500", h_temp[13]);
            file->GetObject("h_PFiso_barrel_nume_200to500", h_temp[14]);
            file->GetObject("h_PFiso_endcap_nume_200to500", h_temp[15]);

            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            removeNegativeBins(h_temp[4]);
            removeNegativeBins(h_temp[5]);
            removeNegativeBins(h_temp[6]);
            removeNegativeBins(h_temp[7]);
            removeNegativeBins(h_temp[8]);
            removeNegativeBins(h_temp[9]);
            removeNegativeBins(h_temp[10]);
            removeNegativeBins(h_temp[11]);
            removeNegativeBins(h_temp[12]);
            removeNegativeBins(h_temp[13]);
            removeNegativeBins(h_temp[14]);
            removeNegativeBins(h_temp[15]);

            h_barrel_data_deno_50to70  ->Add(h_temp[0]);
            h_endcap_data_deno_50to70  ->Add(h_temp[1]);
            h_barrel_data_nume_50to70  ->Add(h_temp[2]);
            h_endcap_data_nume_50to70  ->Add(h_temp[3]);
            h_barrel_data_deno_70to100 ->Add(h_temp[4]);
            h_endcap_data_deno_70to100 ->Add(h_temp[5]);
            h_barrel_data_nume_70to100 ->Add(h_temp[6]);
            h_endcap_data_nume_70to100 ->Add(h_temp[7]);
            h_barrel_data_deno_100to200->Add(h_temp[8]);
            h_endcap_data_deno_100to200->Add(h_temp[9]);
            h_barrel_data_nume_100to200->Add(h_temp[10]);
            h_endcap_data_nume_100to200->Add(h_temp[11]);
            h_barrel_data_deno_200to500->Add(h_temp[12]);
            h_endcap_data_deno_200to500->Add(h_temp[13]);
            h_barrel_data_nume_200to500->Add(h_temp[14]);
            h_endcap_data_nume_200to500->Add(h_temp[15]);
        }
    }

    h_barrel_data_deno_50to70  ->SetDirectory(0);
    h_endcap_data_deno_50to70  ->SetDirectory(0);
    h_barrel_data_nume_50to70  ->SetDirectory(0);
    h_endcap_data_nume_50to70  ->SetDirectory(0);
    h_barrel_data_deno_70to100 ->SetDirectory(0);
    h_endcap_data_deno_70to100 ->SetDirectory(0);
    h_barrel_data_nume_70to100 ->SetDirectory(0);
    h_endcap_data_nume_70to100 ->SetDirectory(0);
    h_barrel_data_deno_100to200->SetDirectory(0);
    h_endcap_data_deno_100to200->SetDirectory(0);
    h_barrel_data_nume_100to200->SetDirectory(0);
    h_endcap_data_nume_100to200->SetDirectory(0);
    h_barrel_data_deno_200to500->SetDirectory(0);
    h_endcap_data_deno_200to500->SetDirectory(0);
    h_barrel_data_nume_200to500->SetDirectory(0);
    h_endcap_data_nume_200to500->SetDirectory(0);

    h_barrel_data_deno_50to70  ->SetLineColor(kBlack);
    h_endcap_data_deno_50to70  ->SetLineColor(kBlack);
    h_barrel_data_nume_50to70  ->SetLineColor(kBlack);
    h_endcap_data_nume_50to70  ->SetLineColor(kBlack);
    h_barrel_data_deno_70to100 ->SetLineColor(kBlack);
    h_endcap_data_deno_70to100 ->SetLineColor(kBlack);
    h_barrel_data_nume_70to100 ->SetLineColor(kBlack);
    h_endcap_data_nume_70to100 ->SetLineColor(kBlack);
    h_barrel_data_deno_100to200->SetLineColor(kBlack);
    h_endcap_data_deno_100to200->SetLineColor(kBlack);
    h_barrel_data_nume_100to200->SetLineColor(kBlack);
    h_endcap_data_nume_100to200->SetLineColor(kBlack);
    h_barrel_data_deno_200to500->SetLineColor(kBlack);
    h_endcap_data_deno_200to500->SetLineColor(kBlack);
    h_barrel_data_nume_200to500->SetLineColor(kBlack);
    h_endcap_data_nume_200to500->SetLineColor(kBlack);

    h_barrel_data_deno_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_deno_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_deno_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_deno_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_deno_100to200->SetMarkerColor(kBlack);
    h_endcap_data_deno_100to200->SetMarkerColor(kBlack);
    h_barrel_data_nume_100to200->SetMarkerColor(kBlack);
    h_endcap_data_nume_100to200->SetMarkerColor(kBlack);
    h_barrel_data_deno_200to500->SetMarkerColor(kBlack);
    h_endcap_data_deno_200to500->SetMarkerColor(kBlack);
    h_barrel_data_nume_200to500->SetMarkerColor(kBlack);
    h_endcap_data_nume_200to500->SetMarkerColor(kBlack);

    h_barrel_data_deno_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_deno_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_deno_100to200->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_100to200->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_100to200->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_100to200->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_deno_200to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_200to500->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_200to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_200to500->SetMarkerStyle(kFullDotLarge);

//-------------------------------------- SCALING -----------------------------------------------------

    h_barrel_MC_deno_50to70[_DY_Full]           ->Scale(2.0067e+06 / h_barrel_MC_deno_50to70[_DY_Full]           ->Integral());
//    h_barrel_MC_deno_50to70[_DY_Full]           ->Scale(3.0767e+06 / h_barrel_MC_deno_50to70[_DY_Full]           ->Integral());
    h_barrel_MC_deno_50to70[_WW]                ->Scale(4.5226e+04 / h_barrel_MC_deno_50to70[_WW]                ->Integral());
    h_barrel_MC_deno_50to70[_WZ]                ->Scale(1.6319e+04 / h_barrel_MC_deno_50to70[_WZ]                ->Integral());
    h_barrel_MC_deno_50to70[_ZZ]                ->Scale(5.0282e+03 / h_barrel_MC_deno_50to70[_ZZ]                ->Integral());
    h_barrel_MC_deno_50to70[_tW]                ->Scale(2.4659e+04 / h_barrel_MC_deno_50to70[_tW]                ->Integral());
    h_barrel_MC_deno_50to70[_tbarW]             ->Scale(2.4717e+04 / h_barrel_MC_deno_50to70[_tbarW]             ->Integral());
    h_barrel_MC_deno_50to70[_ttbar]             ->Scale(10.9162e+05 / h_barrel_MC_deno_50to70[_ttbar]             ->Integral());
//    h_barrel_MC_deno_50to70[_ttbar]             ->Scale(5.9162e+05 / h_barrel_MC_deno_50to70[_ttbar]             ->Integral());
    h_barrel_MC_deno_50to70[_WJets]             ->Scale(11.1813e+06 / h_barrel_MC_deno_50to70[_WJets]             ->Integral());
//    h_barrel_MC_deno_50to70[_WJets]             ->Scale(9.1813e+06 / h_barrel_MC_deno_50to70[_WJets]             ->Integral());
    h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]->Scale(1.5578e+07 / h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]->Integral());
//    h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]->Scale(1.5578e+07 / h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]->Integral());

    s_barrel_deno_50to70 = new THStack("s_barrel_deno_50to70", "");
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_DY_Full]           );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_WW]                );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_WZ]                );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_ZZ]                );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_tW]                );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_tbarW]             );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_ttbar]             );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_WJets]             );
    s_barrel_deno_50to70->Add(h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]);


// ------------------------------------- DRAWING -----------------------------------------------------
    myRatioPlot_t *RP_barrel_deno_50to70 = new myRatioPlot_t("RP_barrel_deno_50to70", s_barrel_deno_50to70, h_barrel_data_deno_50to70);
    RP_barrel_deno_50to70->SetPlots("I_{PF}^{rel.}", 0, 5);

    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_barrel_data_deno_50to70, "Data", "lp");
    legend->AddEntry(h_barrel_MC_deno_50to70[_DY_Full], "DY", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_barrel_MC_deno_50to70[_QCDMuEnriched_Full], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->SetNColumns(2);

    RP_barrel_deno_50to70->ImportLegend(legend);
    RP_barrel_deno_50to70->Draw(1, 1e7, 0);

} // End of Fit_HistDrawer()


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
