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

void Mu_QCDest_HistDrawer(Int_t remNegBins);
void Mu_WJETest_HistDrawer(Int_t remNegBins);

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

void FR_HistDrawer (TString WhichX = "", Int_t type = 2)
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
        cout << "\n*******      E_HistDrawer(" << type << ")      *******" << endl;
        E_HistDrawer(type);
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*******     Mu_QCDest_HistDrawer(" << type << ")     *******" << endl;
            Mu_QCDest_HistDrawer(type);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*******     Mu_WJETest_HistDrawer(" << type << ")     *******" << endl;
            Mu_WJETest_HistDrawer(type);
        }
        else
        {
            cout << "\n*******     Mu_HistDrawer(" << type << ")     *******" << endl;
            Mu_HistDrawer(type);
        }
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
         *h_barrel_MC_deno_100to500[_EndOf_Data_Special],
         *h_barrel_MC_nume_100to500[_EndOf_Data_Special],
         *h_endcap_MC_deno_100to500[_EndOf_Data_Special],
         *h_endcap_MC_nume_100to500[_EndOf_Data_Special],
         *h_barrel_data_deno_50to70,
         *h_barrel_data_nume_50to70,
         *h_endcap_data_deno_50to70,
         *h_endcap_data_nume_50to70,
         *h_barrel_data_deno_70to100,
         *h_barrel_data_nume_70to100,
         *h_endcap_data_deno_70to100,
         *h_endcap_data_nume_70to100,
         *h_barrel_data_deno_100to500,
         *h_barrel_data_nume_100to500,
         *h_endcap_data_deno_100to500,
         *h_endcap_data_nume_100to500;

// ############################# SETUP ################################# //
//----------------------------- MC bkg ------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+".root", "READ");
        file->GetObject("h_pT_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr1]);
        file->GetObject("h_pT_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr1]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr1]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr1]);
        file->GetObject("h_pT_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr1]);
        file->GetObject("h_pT_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr1]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr1]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr1]);
        file->GetObject("h_pT_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr1]);
        file->GetObject("h_pT_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr1]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr1]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr1]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_deno_100to500[pr1]);
        removeNegativeBins(h_endcap_MC_deno_100to500[pr1]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr1]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr1]);

        h_barrel_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_deno_100to500[pr1]->SetDirectory(0);
        h_endcap_MC_deno_100to500[pr1]->SetDirectory(0);
        h_barrel_MC_nume_100to500[pr1]->SetDirectory(0);
        h_endcap_MC_nume_100to500[pr1]->SetDirectory(0);

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
        h_barrel_MC_deno_100to500[pr1]->SetFillColor(color);
        h_endcap_MC_deno_100to500[pr1]->SetFillColor(color);
        h_barrel_MC_nume_100to500[pr1]->SetFillColor(color);
        h_endcap_MC_nume_100to500[pr1]->SetFillColor(color);

        h_barrel_MC_deno_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_deno_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_deno_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_deno_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_deno_100to500[pr1]->SetLineColor(color);
        h_endcap_MC_deno_100to500[pr1]->SetLineColor(color);
        h_barrel_MC_nume_100to500[pr1]->SetLineColor(color);
        h_endcap_MC_nume_100to500[pr1]->SetLineColor(color);

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
    h_barrel_MC_deno_100to500[_ttbar]->Add(h_barrel_MC_deno_100to500[_ttbar_700to1000]);
    h_endcap_MC_deno_100to500[_ttbar]->Add(h_endcap_MC_deno_100to500[_ttbar_700to1000]);
    h_barrel_MC_nume_100to500[_ttbar]->Add(h_barrel_MC_nume_100to500[_ttbar_700to1000]);
    h_endcap_MC_nume_100to500[_ttbar]->Add(h_endcap_MC_nume_100to500[_ttbar_700to1000]);

    h_barrel_MC_deno_50to70  [_ttbar]->Add(h_barrel_MC_deno_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_deno_50to70  [_ttbar]->Add(h_endcap_MC_deno_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_deno_70to100 [_ttbar]->Add(h_barrel_MC_deno_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_deno_70to100 [_ttbar]->Add(h_endcap_MC_deno_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_deno_100to500[_ttbar]->Add(h_barrel_MC_deno_100to500[_ttbar_1000toInf]);
    h_endcap_MC_deno_100to500[_ttbar]->Add(h_endcap_MC_deno_100to500[_ttbar_1000toInf]);
    h_barrel_MC_nume_100to500[_ttbar]->Add(h_barrel_MC_nume_100to500[_ttbar_1000toInf]);
    h_endcap_MC_nume_100to500[_ttbar]->Add(h_endcap_MC_nume_100to500[_ttbar_1000toInf]);

    h_barrel_MC_deno_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_deno_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_deno_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_deno_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_100to500[_ttbar]->SetFillColor(kCyan+2);

    h_barrel_MC_deno_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_deno_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_deno_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_deno_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_100to500[_ttbar]->SetLineColor(kCyan+2);


    h_barrel_MC_deno_50to70  [_WJets]->Add(h_barrel_MC_deno_50to70  [_WJets_ext2v5]);
    h_endcap_MC_deno_50to70  [_WJets]->Add(h_endcap_MC_deno_50to70  [_WJets_ext2v5]);
    h_barrel_MC_nume_50to70  [_WJets]->Add(h_barrel_MC_nume_50to70  [_WJets_ext2v5]);
    h_endcap_MC_nume_50to70  [_WJets]->Add(h_endcap_MC_nume_50to70  [_WJets_ext2v5]);
    h_barrel_MC_deno_70to100 [_WJets]->Add(h_barrel_MC_deno_70to100 [_WJets_ext2v5]);
    h_endcap_MC_deno_70to100 [_WJets]->Add(h_endcap_MC_deno_70to100 [_WJets_ext2v5]);
    h_barrel_MC_nume_70to100 [_WJets]->Add(h_barrel_MC_nume_70to100 [_WJets_ext2v5]);
    h_endcap_MC_nume_70to100 [_WJets]->Add(h_endcap_MC_nume_70to100 [_WJets_ext2v5]);
    h_barrel_MC_deno_100to500[_WJets]->Add(h_barrel_MC_deno_100to500[_WJets_ext2v5]);
    h_endcap_MC_deno_100to500[_WJets]->Add(h_endcap_MC_deno_100to500[_WJets_ext2v5]);
    h_barrel_MC_nume_100to500[_WJets]->Add(h_barrel_MC_nume_100to500[_WJets_ext2v5]);
    h_endcap_MC_nume_100to500[_WJets]->Add(h_endcap_MC_nume_100to500[_WJets_ext2v5]);

    h_barrel_MC_deno_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_deno_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_deno_100to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_deno_100to500[_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_100to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_100to500[_WJets]->SetFillColor(kRed-2);

    h_barrel_MC_deno_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_deno_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_deno_100to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_deno_100to500[_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_100to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_100to500[_WJets]->SetLineColor(kRed-2);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_pT_barrel_deno_50to70",   h_barrel_MC_deno_50to70  [pr]);
        file->GetObject("h_pT_endcap_deno_50to70",   h_endcap_MC_deno_50to70  [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70  [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70  [pr]);
        file->GetObject("h_pT_barrel_deno_70to100",  h_barrel_MC_deno_70to100 [pr]);
        file->GetObject("h_pT_endcap_deno_70to100",  h_endcap_MC_deno_70to100 [pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100 [pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100 [pr]);
        file->GetObject("h_pT_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr]);
        file->GetObject("h_pT_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno_100to500[pr]);
        removeNegativeBins(h_endcap_MC_deno_100to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);;
        h_barrel_MC_deno_100to500[pr]->SetDirectory(0);
        h_endcap_MC_deno_100to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to500[pr]->SetDirectory(0);

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
            h_barrel_MC_deno_100to500[_DY_Full] = ((TH1D*)(h_barrel_MC_deno_100to500[pr]->Clone("h_barrel_MC_deno_DY_100to500")));
            h_endcap_MC_deno_100to500[_DY_Full] = ((TH1D*)(h_endcap_MC_deno_100to500[pr]->Clone("h_endcap_MC_deno_DY_100to500")));
            h_barrel_MC_nume_100to500[_DY_Full] = ((TH1D*)(h_barrel_MC_nume_100to500[pr]->Clone("h_barrel_MC_nume_DY_100to500")));
            h_endcap_MC_nume_100to500[_DY_Full] = ((TH1D*)(h_endcap_MC_nume_100to500[pr]->Clone("h_endcap_MC_nume_DY_100to500")));

            h_barrel_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_100to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_100to500[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to500[_DY_Full]->SetDirectory(0);
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
            h_barrel_MC_deno_100to500[_DY_Full]->Add(h_barrel_MC_deno_100to500[pr]);
            h_endcap_MC_deno_100to500[_DY_Full]->Add(h_endcap_MC_deno_100to500[pr]);
            h_barrel_MC_nume_100to500[_DY_Full]->Add(h_barrel_MC_nume_100to500[pr]);
            h_endcap_MC_nume_100to500[_DY_Full]->Add(h_endcap_MC_nume_100to500[pr]);
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
    h_barrel_MC_deno_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_deno_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_100to500[_DY_Full]->SetFillColor(kOrange-5);

    h_barrel_MC_deno_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_deno_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_deno_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_deno_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_100to500[_DY_Full]->SetLineColor(kOrange-5);

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_pT_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr]);
        file->GetObject("h_pT_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr]);
        file->GetObject("h_pT_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr]);
        file->GetObject("h_pT_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr]);
        file->GetObject("h_pT_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr]);
        file->GetObject("h_pT_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno_100to500[pr]);
        removeNegativeBins(h_endcap_MC_deno_100to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_deno_100to500[pr]->SetDirectory(0);
        h_endcap_MC_deno_100to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to500[pr]->SetDirectory(0);

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
            h_barrel_MC_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno_100to500[pr]->Clone("h_barrel_MC_deno_QCD_100to500")));
            h_endcap_MC_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno_100to500[pr]->Clone("h_endcap_MC_deno_QCD_100to500")));
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_100to500[pr]->Clone("h_barrel_MC_nume_QCD_100to500")));
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_100to500[pr]->Clone("h_endcap_MC_nume_QCD_100to500")));

            h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
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
            h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_deno_100to500[pr]);
            h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_deno_100to500[pr]);
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_100to500[pr]);
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_100to500[pr]);
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
    h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);

    h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        TH1D *h_temp[12];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_pT_barrel_deno_50to70",   h_barrel_data_deno_50to70 );
            file->GetObject("h_pT_endcap_deno_50to70",   h_endcap_data_deno_50to70 );
            file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_data_nume_50to70 );
            file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_data_nume_50to70 );
            file->GetObject("h_pT_barrel_deno_70to100",  h_barrel_data_deno_70to100);
            file->GetObject("h_pT_endcap_deno_70to100",  h_endcap_data_deno_70to100);
            file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_data_nume_70to100);
            file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_data_nume_70to100);
            file->GetObject("h_pT_barrel_deno_100to500", h_barrel_data_deno_100to500);
            file->GetObject("h_pT_endcap_deno_100to500", h_endcap_data_deno_100to500);
            file->GetObject("h_pT_barrel_nume_100to500", h_barrel_data_nume_100to500);
            file->GetObject("h_pT_endcap_nume_100to500", h_endcap_data_nume_100to500);

            removeNegativeBins(h_barrel_data_deno_50to70);
            removeNegativeBins(h_endcap_data_deno_50to70);
            removeNegativeBins(h_barrel_data_nume_50to70);
            removeNegativeBins(h_endcap_data_nume_50to70);
            removeNegativeBins(h_barrel_data_deno_70to100);
            removeNegativeBins(h_endcap_data_deno_70to100);
            removeNegativeBins(h_barrel_data_nume_70to100);
            removeNegativeBins(h_endcap_data_nume_70to100);
            removeNegativeBins(h_barrel_data_deno_100to500);
            removeNegativeBins(h_endcap_data_deno_100to500);
            removeNegativeBins(h_barrel_data_nume_100to500);
            removeNegativeBins(h_endcap_data_nume_100to500);
        }
        else
        {
            file->GetObject("h_pT_barrel_deno_50to70",   h_temp[0]);
            file->GetObject("h_pT_endcap_deno_50to70",   h_temp[1]);
            file->GetObject("h_pT_barrel_nume_50to70",   h_temp[2]);
            file->GetObject("h_pT_endcap_nume_50to70",   h_temp[3]);
            file->GetObject("h_pT_barrel_deno_70to100",  h_temp[4]);
            file->GetObject("h_pT_endcap_deno_70to100",  h_temp[5]);
            file->GetObject("h_pT_barrel_nume_70to100",  h_temp[6]);
            file->GetObject("h_pT_endcap_nume_70to100",  h_temp[7]);
            file->GetObject("h_pT_barrel_deno_100to500", h_temp[8]);
            file->GetObject("h_pT_endcap_deno_100to500", h_temp[9]);
            file->GetObject("h_pT_barrel_nume_100to500", h_temp[10]);
            file->GetObject("h_pT_endcap_nume_100to500", h_temp[11]);

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

            h_barrel_data_deno_50to70  ->Add(h_temp[0]);
            h_endcap_data_deno_50to70  ->Add(h_temp[1]);
            h_barrel_data_nume_50to70  ->Add(h_temp[2]);
            h_endcap_data_nume_50to70  ->Add(h_temp[3]);
            h_barrel_data_deno_70to100 ->Add(h_temp[4]);
            h_endcap_data_deno_70to100 ->Add(h_temp[5]);
            h_barrel_data_nume_70to100 ->Add(h_temp[6]);
            h_endcap_data_nume_70to100 ->Add(h_temp[7]);
            h_barrel_data_deno_100to500->Add(h_temp[8]);
            h_endcap_data_deno_100to500->Add(h_temp[9]);
            h_barrel_data_nume_100to500->Add(h_temp[10]);
            h_endcap_data_nume_100to500->Add(h_temp[11]);
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
    h_barrel_data_deno_100to500->SetDirectory(0);
    h_endcap_data_deno_100to500->SetDirectory(0);
    h_barrel_data_nume_100to500->SetDirectory(0);
    h_endcap_data_nume_100to500->SetDirectory(0);

    h_barrel_data_deno_50to70  ->SetLineColor(kBlack);
    h_endcap_data_deno_50to70  ->SetLineColor(kBlack);
    h_barrel_data_nume_50to70  ->SetLineColor(kBlack);
    h_endcap_data_nume_50to70  ->SetLineColor(kBlack);
    h_barrel_data_deno_70to100 ->SetLineColor(kBlack);
    h_endcap_data_deno_70to100 ->SetLineColor(kBlack);
    h_barrel_data_nume_70to100 ->SetLineColor(kBlack);
    h_endcap_data_nume_70to100 ->SetLineColor(kBlack);
    h_barrel_data_deno_100to500->SetLineColor(kBlack);
    h_endcap_data_deno_100to500->SetLineColor(kBlack);
    h_barrel_data_nume_100to500->SetLineColor(kBlack);
    h_endcap_data_nume_100to500->SetLineColor(kBlack);

    h_barrel_data_deno_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_deno_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_deno_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_deno_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_deno_100to500->SetMarkerColor(kBlack);
    h_endcap_data_deno_100to500->SetMarkerColor(kBlack);
    h_barrel_data_nume_100to500->SetMarkerColor(kBlack);
    h_endcap_data_nume_100to500->SetMarkerColor(kBlack);

    h_barrel_data_deno_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_deno_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_deno_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_deno_100to500->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_100to500->SetMarkerStyle(kFullDotLarge);

//-------------------------------------- SCALING -----------------------------------------------------

    h_barrel_MC_deno_50to70[_DY_Full]             ->Scale(6.0283e+05 / h_barrel_MC_deno_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_deno_50to70[_WW]                  ->Scale(2.9979e+04 / h_barrel_MC_deno_50to70[_WW]                 ->Integral());
    h_barrel_MC_deno_50to70[_WZ]                  ->Scale(7.5397e+03 / h_barrel_MC_deno_50to70[_WZ]                 ->Integral());
    h_barrel_MC_deno_50to70[_ZZ]                  ->Scale(1.4975e+03 / h_barrel_MC_deno_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_deno_50to70[_tW]                  ->Scale(1.3394e+04 / h_barrel_MC_deno_50to70[_tW]                 ->Integral());
    h_barrel_MC_deno_50to70[_tbarW]               ->Scale(1.3531e+04 / h_barrel_MC_deno_50to70[_tbarW]              ->Integral());
    h_barrel_MC_deno_50to70[_ttbar]               ->Scale(2.7678e+05 / h_barrel_MC_deno_50to70[_ttbar]              ->Integral());
    h_barrel_MC_deno_50to70[_WJets]               ->Scale(9.2318e+06 / h_barrel_MC_deno_50to70[_WJets]              ->Integral());
    h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]  ->Scale(1.3530e+07 / h_barrel_MC_deno_50to70[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_deno_70to100[_DY_Full]            ->Scale(2.8544e+05 / h_barrel_MC_deno_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_deno_70to100[_WW]                 ->Scale(4.5226e+04 / h_barrel_MC_deno_50to70[_WW]                 ->Integral());
    h_barrel_MC_deno_70to100[_WZ]                 ->Scale(5.0573e+03 / h_barrel_MC_deno_50to70[_WZ]                 ->Integral());
    h_barrel_MC_deno_70to100[_ZZ]                 ->Scale(7.2227e+02 / h_barrel_MC_deno_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_deno_70to100[_tW]                 ->Scale(1.2115e+04 / h_barrel_MC_deno_50to70[_tW]                 ->Integral());
    h_barrel_MC_deno_70to100[_tbarW]              ->Scale(1.2028e+04 / h_barrel_MC_deno_50to70[_tbarW]              ->Integral());
    h_barrel_MC_deno_70to100[_ttbar]              ->Scale(2.3233e+05 / h_barrel_MC_deno_50to70[_ttbar]              ->Integral());
    h_barrel_MC_deno_70to100[_WJets]              ->Scale(3.4985e+06 / h_barrel_MC_deno_50to70[_WJets]              ->Integral());
    h_barrel_MC_deno_70to100[_QCDMuEnriched_Full] ->Scale(3.7908e+06 / h_barrel_MC_deno_50to70[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_deno_100to500[_DY_Full]           ->Scale(1.0562e+05 / h_barrel_MC_deno_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_deno_100to500[_WW]                ->Scale(1.4394e+04 / h_barrel_MC_deno_50to70[_WW]                 ->Integral());
    h_barrel_MC_deno_100to500[_WZ]                ->Scale(5.1568e+03 / h_barrel_MC_deno_50to70[_WZ]                 ->Integral());
    h_barrel_MC_deno_100to500[_ZZ]                ->Scale(8.3444e+02 / h_barrel_MC_deno_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_deno_100to500[_tW]                ->Scale(1.2316e+04 / h_barrel_MC_deno_50to70[_tW]                 ->Integral());
    h_barrel_MC_deno_100to500[_tbarW]             ->Scale(1.2003e+04 / h_barrel_MC_deno_50to70[_tbarW]              ->Integral());
    h_barrel_MC_deno_100to500[_ttbar]             ->Scale(2.7302e+05 / h_barrel_MC_deno_50to70[_ttbar]              ->Integral());
    h_barrel_MC_deno_100to500[_WJets]             ->Scale(1.3415e+06 / h_barrel_MC_deno_50to70[_WJets]              ->Integral());
    h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->Scale(7.9078e+05 / h_barrel_MC_deno_50to70[_QCDMuEnriched_Full] ->Integral());

    THStack * s_barrel_deno = new THStack("s_barrel_deno", "");
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_ZZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_tW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_tbarW]             );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_ttbar]             );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_DY_Full]           );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WJets]             );
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_ZZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_tW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_tbarW]             );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_ttbar]             );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_DY_Full]           );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WJets]             );
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_ZZ]                );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_tW]                );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_tbarW]             );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_ttbar]             );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_DY_Full]           );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WJets]             );
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]);

    h_barrel_data_deno_50to70->Add(h_barrel_data_deno_70to100);
    h_barrel_data_deno_50to70->Add(h_barrel_data_deno_100to500);


// ------------------------------------- DRAWING -----------------------------------------------------
    myRatioPlot_t *RP_barrel_deno = new myRatioPlot_t("RP_barrel_deno_50to70", s_barrel_deno, h_barrel_data_deno_70to100);
    RP_barrel_deno->SetPlots("p_{T} (#mu_{barrel}^{deno})", 0, 5);

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

    RP_barrel_deno->ImportLegend(legend);
    RP_barrel_deno->Draw(1, 1e7, 0);

} // End of Fit_HistDrawer()


void Mu_QCDest_HistDrawer(Int_t remNegBins)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";

    f = new TFile(Dir+"QCDest_Mu.root", "READ");

    TH1D * h_mass[_EndOf_Data_Special];
    TH1D * h_QCD_est;
    THStack * s_mass_wQCD = new THStack("s_mass_wQCD", "");
    THStack * s_mass_woQCD = new THStack("s_mass_woQCD", "");
    Color_t color = kBlack;

    for (Process_t pr=_DY_10to50; pr<_EndOf_SinglMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        h_mass[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1) removeNegativeBins(h_mass[pr]);

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_QCDMuEnriched_Normal) color = kRed + 3;

        if (pr < _SingleMuon_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
        }

        if (pr < _SingleMuon_B) s_mass_wQCD->Add(h_mass[pr]);
        if (pr < _QCDMuEnriched_15to20) s_mass_woQCD->Add(h_mass[pr]);

        // Adding up for convenience
        if (pr == _DY_10to50)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal) h_mass[_DY_Full]->Add(h_mass[pr]);
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal) h_mass[_ttbar_Full]->Add(h_mass[pr]);
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal) h_mass[_VVnST]->Add(h_mass[pr]);
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal) h_mass[_WJets_Full]->Add(h_mass[pr]);
        else if (pr == _QCDMuEnriched_15to20)
        {
            h_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDMuEnriched_Full")));
            h_mass[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_QCDMuEnriched_Normal) h_mass[_QCDMuEnriched_Full]->Add(h_mass[pr]);
        else if (pr == _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_SinglMuon_Normal) h_mass[_SingleMuon_Full]->Add(h_mass[pr]);


        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

    } // End of pr iteration

    // QCD estimation
    h_QCD_est = ((TH1D*)(h_mass[_SingleMuon_Full]->Clone("h_QCD_est")));
    h_QCD_est->SetTitle("");
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Add(h_mass[_DY_Full], -1);
    h_QCD_est->Add(h_mass[_ttbar_Full], -1);
//    h_QCD_est->Add(h_mass[_VVnST], -1);
//    h_QCD_est->Add(h_mass[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kRed + 3);

    myRatioPlot_t *RP_mass_wQCD = new myRatioPlot_t("c_mass_wQCD", s_mass_wQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_SingleMuon_Full]);

    RP_mass_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend * legend_wQCD = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend_wQCD->AddEntry(h_mass[_SingleMuon_B], "Data", "pl");
    legend_wQCD->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_wQCD->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_wQCD->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend_wQCD->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend_wQCD->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend_wQCD->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend_wQCD->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend_wQCD->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    TLegend * legend_woQCD = ((TLegend*)(legend_wQCD->Clone()));
    legend_wQCD->AddEntry(h_mass[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wQCD->ImportLegend(legend_wQCD);
    RP_mass_woQCD->ImportLegend(legend_woQCD);

    RP_mass_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD->Draw(1e-3, 1e3, 1);

    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    TCanvas * c_QCD_est = new TCanvas("c_QCD_est", "c_QCD_est", 750, 850);
    c_QCD_est->SetTopMargin(0.05);
    c_QCD_est->SetRightMargin(0.05);
    c_QCD_est->SetBottomMargin(0.15);
    c_QCD_est->SetLeftMargin(0.15);
    h_QCD_est->Draw("hist");
    h_QCD_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_QCD_est->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est->GetXaxis()->SetMoreLogLabels();
    h_QCD_est->GetXaxis()->SetNoExponent();
    h_QCD_est->GetYaxis()->SetTitle("Number of events");
    h_QCD_est->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est->GetYaxis()->SetTitleOffset(1.2);
    h_QCD_est->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est->GetYaxis()->SetMoreLogLabels();
    h_QCD_est->GetYaxis()->SetNoExponent();
    l_QCD_est->Draw();
    c_QCD_est->SetLogx();
    c_QCD_est->SetLogy();
    c_QCD_est->SetGridx();
    c_QCD_est->SetGridy();
    c_QCD_est->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_Mu.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_Mu.root" << " COULD NOT BE CLOSED!\n" << endl;

} // End of Mu_QCDest_HistDrawer()


void Mu_WJETest_HistDrawer(Int_t remNegBins)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";

    f = new TFile(Dir+"WJETest_Mu.root", "READ");

    TH1D * h_mass[_EndOf_Data_Special];
    TH1D * h_WJET_est;
    THStack * s_mass_wWJET = new THStack("s_mass_wWJET", "");
    THStack * s_mass_woWJET = new THStack("s_mass_woWJET", "");
    Color_t color = kBlack;

    for (Process_t pr=_SingleMuon_H; pr>=_DY_10to50; pr=previous(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        h_mass[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1) removeNegativeBins(h_mass[pr]);

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_QCDMuEnriched_Normal) color = kRed + 3;

        if (pr < _SingleMuon_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
        }

        if (pr < _SingleMuon_B)
        {
            s_mass_wWJET->Add(h_mass[pr]);
            if (pr < _WJets || pr > _EndOf_WJets_Normal) s_mass_woWJET->Add(h_mass[pr]);
        }

        // Adding up for convenience
        if (pr == _SingleMuon_H)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr >= _SingleMuon_B) h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
        else if (pr == _QCDMuEnriched_1000toInf)
        {
            h_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDMuEnriched_Full")));
            h_mass[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else if (pr >= _QCDMuEnriched_15to20) h_mass[_QCDMuEnriched_Full]->Add(h_mass[pr]);
        else if (pr == _WJets_ext2v5)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
        }
        else if (pr >= _WJets) h_mass[_WJets_Full]->Add(h_mass[pr]);
        else if (pr == _WW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
        }
        else if (pr >= _tW) h_mass[_VVnST]->Add(h_mass[pr]);
        else if (pr == _ttbar_1000toInf)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr >= _ttbar) h_mass[_ttbar_Full]->Add(h_mass[pr]);
        else if (pr == _DY_2000to3000)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
        }
        else if (pr >= _DY_10to50) h_mass[_DY_Full]->Add(h_mass[pr]);

        if (pr == _SingleMuon_B) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCD_Mu_1000toInf
        if (pr == _ttbar) pr = _EndOf_DY_Normal; // next -- DY_2000to3000

    } // End of pr iteration

    // For drawing MC W+Jets on top
    h_mass[_WJets_Full]->SetFillColor(kGray);
    h_mass[_WJets_Full]->SetFillStyle(3002);
    h_mass[_WJets_Full]->SetLineColor(kRed);

    // W+Jets estimation
    h_WJET_est = ((TH1D*)(h_mass[_SingleMuon_Full]->Clone("h_WJET_est")));
    h_WJET_est->SetTitle("");
    h_WJET_est->SetDirectory(0);
    h_WJET_est->Add(h_mass[_DY_Full], -1);
    h_WJET_est->Add(h_mass[_ttbar_Full], -1);
    h_WJET_est->Add(h_mass[_VVnST], -1);
    h_WJET_est->Add(h_mass[_QCDMuEnriched_Full], -1);
    removeNegativeBins(h_WJET_est);
    h_WJET_est->SetFillColor(kRed - 2);
    h_WJET_est->SetLineColor(kRed - 2);

    myRatioPlot_t *RP_mass_wWJET = new myRatioPlot_t("c_mass_wWJET", s_mass_wWJET, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woWJET = new myRatioPlot_t("c_mass_woWJET", s_mass_woWJET, h_mass[_SingleMuon_Full]);

    RP_mass_wWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend * legend_wWJET = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend_wWJET->AddEntry(h_mass[_SingleMuon_B], "Data", "pl");
    legend_wWJET->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_wWJET->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend_wWJET->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend_wWJET->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend_wWJET->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend_wWJET->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    TLegend * legend_woWJET = ((TLegend*)(legend_wWJET->Clone()));
    legend_wWJET->AddEntry(h_mass[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);

    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
//    RP_mass_wWJET->DrawOnTop(h_mass[_WJets_Full]);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
    TCanvas * c_WJET_est = new TCanvas("c_WJET_est", "c_WJET_est", 750, 850);
    c_WJET_est->SetTopMargin(0.05);
    c_WJET_est->SetRightMargin(0.05);
    c_WJET_est->SetBottomMargin(0.15);
    c_WJET_est->SetLeftMargin(0.15);
    h_WJET_est->Draw("hist");
    h_WJET_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_WJET_est->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est->GetXaxis()->SetMoreLogLabels();
    h_WJET_est->GetXaxis()->SetNoExponent();
    h_WJET_est->GetYaxis()->SetTitle("Number of events");
    h_WJET_est->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est->GetYaxis()->SetTitleOffset(1.2);
    h_WJET_est->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est->GetYaxis()->SetMoreLogLabels();
    h_WJET_est->GetYaxis()->SetNoExponent();
    l_WJET_est->Draw();
    c_WJET_est->SetLogx();
    c_WJET_est->SetLogy();
    c_WJET_est->SetGridx();
    c_WJET_est->SetGridy();
    c_WJET_est->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_Mu.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_Mu.root" << " COULD NOT BE CLOSED!\n" << endl;

} // End of Mu_WJETest_HistDrawer()


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
