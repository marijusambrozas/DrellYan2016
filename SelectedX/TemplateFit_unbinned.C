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
#include <RooFit.h>
#include <RooDataHist.h>
//#include <RooHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>


// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

using namespace RooFit;

void E_Tfit(Int_t type);
void Mu_Tfit(Int_t type);
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

void TemplateFit_unbinned (TString WhichX = "", Int_t type = 2)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    if (type < 1 || type > 2)
    {
        cout << "Wrong type!" << endl;
        return;
    }
    if (whichX.Contains("E"))
    {
        Xselected++;
        cout << "\n*******      E_Tfit(" << type << ")      *******" << endl;
        E_Tfit(type);
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*******     Mu_Tfit(" << type << ")     *******" << endl;
        Mu_Tfit(type);
    }
    if (Xselected == 0) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ############################################################################# ///
/// ----------------------------- Electron Channel ------------------------------ ///
/// ############################################################################# ///
void E_Tfit(Int_t type)
{
   return; // NOT READY YET
} // End of E_Tfit()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_Tfit(Int_t type)
{
    FileMgr fm;

    TH1D *h_barrel_MC_deno[_EndOf_Data_Special], *h_barrel_MC_nume[_EndOf_Data_Special],
         *h_endcap_MC_deno[_EndOf_Data_Special], *h_endcap_MC_nume[_EndOf_Data_Special],
         *h_barrel_data_deno, *h_barrel_data_nume, *h_endcap_data_deno, *h_endcap_data_nume;

// ############################# SETUP ################################# //
//----------------------------- MC bkg ------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr1]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_endcap_deno", h_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_barrel_nume", h_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_endcap_nume", h_endcap_MC_nume[pr1]);

        removeNegativeBins(h_barrel_MC_deno[pr1]);
        removeNegativeBins(h_endcap_MC_deno[pr1]);
        removeNegativeBins(h_barrel_MC_nume[pr1]);
        removeNegativeBins(h_endcap_MC_nume[pr1]);

        h_barrel_MC_deno[pr1]->SetDirectory(0);
        h_endcap_MC_deno[pr1]->SetDirectory(0);
        h_barrel_MC_nume[pr1]->SetDirectory(0);
        h_endcap_MC_nume[pr1]->SetDirectory(0);

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
    h_barrel_MC_deno[_ttbar]->Add(h_barrel_MC_deno[_ttbar_700to1000]);
    h_endcap_MC_deno[_ttbar]->Add(h_endcap_MC_deno[_ttbar_700to1000]);
    h_barrel_MC_nume[_ttbar]->Add(h_barrel_MC_nume[_ttbar_700to1000]);
    h_endcap_MC_nume[_ttbar]->Add(h_endcap_MC_nume[_ttbar_700to1000]);
    h_barrel_MC_deno[_ttbar]->Add(h_barrel_MC_deno[_ttbar_1000toInf]);
    h_endcap_MC_deno[_ttbar]->Add(h_endcap_MC_deno[_ttbar_1000toInf]);
    h_barrel_MC_nume[_ttbar]->Add(h_barrel_MC_nume[_ttbar_1000toInf]);
    h_endcap_MC_nume[_ttbar]->Add(h_endcap_MC_nume[_ttbar_1000toInf]);
    h_barrel_MC_deno[_WJets]->Add(h_barrel_MC_deno[_WJets_ext2v5]);
    h_endcap_MC_deno[_WJets]->Add(h_endcap_MC_deno[_WJets_ext2v5]);
    h_barrel_MC_nume[_WJets]->Add(h_barrel_MC_nume[_WJets_ext2v5]);
    h_endcap_MC_nume[_WJets]->Add(h_endcap_MC_nume[_WJets_ext2v5]);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume", h_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_endcap_MC_nume[pr]);
        removeNegativeBins(h_barrel_MC_deno[pr]);
        removeNegativeBins(h_endcap_MC_deno[pr]);
        removeNegativeBins(h_barrel_MC_nume[pr]);
        removeNegativeBins(h_endcap_MC_nume[pr]);

        h_barrel_MC_deno[pr]->SetDirectory(0);
        h_endcap_MC_deno[pr]->SetDirectory(0);
        h_barrel_MC_nume[pr]->SetDirectory(0);
        h_endcap_MC_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_barrel_MC_deno[pr]->Clone("h_barrel_MC_deno_DY")));
            h_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_endcap_MC_deno[pr]->Clone("h_endcap_MC_deno_DY")));
            h_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_barrel_MC_nume[pr]->Clone("h_barrel_MC_nume_DY")));
            h_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_endcap_MC_nume[pr]->Clone("h_endcap_MC_nume_DY")));
            h_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_deno[_DY_Full]->Add(h_barrel_MC_deno[pr]);
            h_endcap_MC_deno[_DY_Full]->Add(h_endcap_MC_deno[pr]);
            h_barrel_MC_nume[_DY_Full]->Add(h_barrel_MC_nume[pr]);
            h_endcap_MC_nume[_DY_Full]->Add(h_endcap_MC_nume[pr]);
        }
        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno", h_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume", h_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_endcap_MC_nume[pr]);
        removeNegativeBins(h_barrel_MC_deno[pr]);
        removeNegativeBins(h_endcap_MC_deno[pr]);
        removeNegativeBins(h_barrel_MC_nume[pr]);
        removeNegativeBins(h_endcap_MC_nume[pr]);

        h_barrel_MC_deno[pr]->SetDirectory(0);
        h_endcap_MC_deno[pr]->SetDirectory(0);
        h_barrel_MC_nume[pr]->SetDirectory(0);
        h_endcap_MC_nume[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno[pr]->Clone("h_barrel_MC_deno_QCD")));
            h_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno[pr]->Clone("h_endcap_MC_deno_QCD")));
            h_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume[pr]->Clone("h_barrel_MC_nume_QCD")));
            h_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume[pr]->Clone("h_endcap_MC_nume_QCD")));
            h_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_barrel_MC_deno[pr]);
            h_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_endcap_MC_deno[pr]);
            h_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_barrel_MC_nume[pr]);
            h_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_endcap_MC_nume[pr]);
        }
        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        TH1D *h_temp[4];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_PFiso_barrel_deno", h_barrel_data_deno);
            file->GetObject("h_PFiso_endcap_deno", h_endcap_data_deno);
            file->GetObject("h_PFiso_barrel_nume", h_barrel_data_nume);
            file->GetObject("h_PFiso_endcap_nume", h_endcap_data_nume);
            removeNegativeBins(h_barrel_data_deno);
            removeNegativeBins(h_endcap_data_deno);
            removeNegativeBins(h_barrel_data_nume);
            removeNegativeBins(h_endcap_data_nume);
        }
        else
        {
            file->GetObject("h_PFiso_barrel_deno", h_temp[0]);
            file->GetObject("h_PFiso_endcap_deno", h_temp[1]);
            file->GetObject("h_PFiso_barrel_nume", h_temp[2]);
            file->GetObject("h_PFiso_endcap_nume", h_temp[3]);
            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            h_barrel_data_deno->Add(h_temp[0]);
            h_endcap_data_deno->Add(h_temp[1]);
            h_barrel_data_nume->Add(h_temp[2]);
            h_endcap_data_nume->Add(h_temp[3]);
        }
    }

    h_barrel_data_deno->SetDirectory(0);
    h_endcap_data_deno->SetDirectory(0);
    h_barrel_data_nume->SetDirectory(0);
    h_endcap_data_nume->SetDirectory(0);

// ######################## MODEL BUILDING ##########################

    // Making RooDataHist
    RooRealVar iso_nume("iso", "PFiso/p_{T}", 0, 0.15);
    RooRealVar iso_deno("iso", "PFiso/p_{T}", 0, 5);

    RooDataHist *rh_barrel_nume_QCD = new RooDataHist("rh_barrel_nume_QCD", "RooHist_barrel_nume_QCD", iso_nume, h_barrel_MC_nume[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_nume_WJets = new RooDataHist("rh_barrel_nume_WJets", "RooHist_barrel_nume_WJets", iso_nume, h_barrel_MC_nume[_WJets]);
    RooDataHist *rh_barrel_nume_DY = new RooDataHist("rh_barrel_nume_DY", "RooHist_barrel_nume_DY", iso_nume, h_barrel_MC_nume[_DY_Full]);
    RooDataHist *rh_barrel_nume_ttbar = new RooDataHist("rh_barrel_nume_ttbar", "RooHist_barrel_nume_ttbar", iso_nume, h_barrel_MC_nume[_ttbar]);
    RooDataHist *rh_barrel_nume_tW = new RooDataHist("rh_barrel_nume_tW", "RooHist_barrel_nume_tW", iso_nume, h_barrel_MC_nume[_tW]);
    RooDataHist *rh_barrel_nume_tbarW = new RooDataHist("rh_barrel_nume_tbarW", "RooHist_barrel_nume_tbarW", iso_nume, h_barrel_MC_nume[_tbarW]);
    RooDataHist *rh_barrel_nume_WW = new RooDataHist("rh_barrel_nume_WW", "RooHist_barrel_nume_WW", iso_nume, h_barrel_MC_nume[_WW]);
    RooDataHist *rh_barrel_nume_WZ = new RooDataHist("rh_barrel_nume_WZ", "RooHist_barrel_nume_WZ", iso_nume, h_barrel_MC_nume[_WZ]);
    RooDataHist *rh_barrel_nume_ZZ = new RooDataHist("rh_barrel_nume_ZZ", "RooHist_barrel_nume_ZZ", iso_nume, h_barrel_MC_nume[_ZZ]);
    RooDataHist *rh_barrel_nume_data = new RooDataHist("rh_barrel_nume_data", "RooHist_barrel_nume_data", iso_nume, h_barrel_data_nume);

    RooDataHist *rh_endcap_nume_QCD = new RooDataHist("rh_endcap_nume_QCD", "RooHist_endcap_nume_QCD", iso_nume, h_endcap_MC_nume[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_nume_WJets = new RooDataHist("rh_endcap_nume_WJets", "RooHist_endcap_nume_WJets", iso_nume, h_endcap_MC_nume[_WJets]);
    RooDataHist *rh_endcap_nume_DY = new RooDataHist("rh_endcap_nume_DY", "RooHist_endcap_nume_DY", iso_nume, h_endcap_MC_nume[_DY_Full]);
    RooDataHist *rh_endcap_nume_ttbar = new RooDataHist("rh_endcap_nume_ttbar", "RooHist_endcap_nume_ttbar", iso_nume, h_endcap_MC_nume[_ttbar]);
    RooDataHist *rh_endcap_nume_tW = new RooDataHist("rh_endcap_nume_tW", "RooHist_endcap_nume_tW", iso_nume, h_endcap_MC_nume[_tW]);
    RooDataHist *rh_endcap_nume_tbarW = new RooDataHist("rh_endcap_nume_tbarW", "RooHist_endcap_nume_tbarW", iso_nume, h_endcap_MC_nume[_tbarW]);
    RooDataHist *rh_endcap_nume_WW = new RooDataHist("rh_endcap_nume_WW", "RooHist_endcap_nume_WW", iso_nume, h_endcap_MC_nume[_WW]);
    RooDataHist *rh_endcap_nume_WZ = new RooDataHist("rh_endcap_nume_WZ", "RooHist_endcap_nume_WZ", iso_nume, h_endcap_MC_nume[_WZ]);
    RooDataHist *rh_endcap_nume_ZZ = new RooDataHist("rh_endcap_nume_ZZ", "RooHist_endcap_nume_ZZ", iso_nume, h_endcap_MC_nume[_ZZ]);
    RooDataHist *rh_endcap_nume_data = new RooDataHist("rh_endcap_nume_data", "RooHist_endcap_nume_data", iso_nume, h_endcap_data_nume);

    RooDataHist *rh_barrel_deno_QCD = new RooDataHist("rh_barrel_deno_QCD", "RooHist_barrel_deno_QCD", iso_deno, h_barrel_MC_deno[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_deno_WJets = new RooDataHist("rh_barrel_deno_WJets", "RooHist_barrel_deno_WJets", iso_deno, h_barrel_MC_deno[_WJets]);
    RooDataHist *rh_barrel_deno_DY = new RooDataHist("rh_barrel_deno_DY", "RooHist_barrel_deno_DY", iso_deno, h_barrel_MC_deno[_DY_Full]);
    RooDataHist *rh_barrel_deno_ttbar = new RooDataHist("rh_barrel_deno_ttbar", "RooHist_barrel_deno_ttbar", iso_deno, h_barrel_MC_deno[_ttbar]);
    RooDataHist *rh_barrel_deno_tW = new RooDataHist("rh_barrel_deno_tW", "RooHist_barrel_deno_tW", iso_deno, h_barrel_MC_deno[_tW]);
    RooDataHist *rh_barrel_deno_tbarW = new RooDataHist("rh_barrel_deno_tbarW", "RooHist_barrel_deno_tbarW", iso_deno, h_barrel_MC_deno[_tbarW]);
    RooDataHist *rh_barrel_deno_WW = new RooDataHist("rh_barrel_deno_WW", "RooHist_barrel_deno_WW", iso_deno, h_barrel_MC_deno[_WW]);
    RooDataHist *rh_barrel_deno_WZ = new RooDataHist("rh_barrel_deno_WZ", "RooHist_barrel_deno_WZ", iso_deno, h_barrel_MC_deno[_WZ]);
    RooDataHist *rh_barrel_deno_ZZ = new RooDataHist("rh_barrel_deno_ZZ", "RooHist_barrel_deno_ZZ", iso_deno, h_barrel_MC_deno[_ZZ]);
    RooDataHist *rh_barrel_deno_data = new RooDataHist("rh_barrel_deno_data", "RooHist_barrel_deno_data", iso_deno, h_barrel_data_deno);

    RooDataHist *rh_endcap_deno_QCD = new RooDataHist("rh_endcap_deno_QCD", "RooHist_endcap_deno_QCD", iso_deno, h_endcap_MC_deno[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_deno_WJets = new RooDataHist("rh_endcap_deno_WJets", "RooHist_endcap_deno_WJets", iso_deno, h_endcap_MC_deno[_WJets]);
    RooDataHist *rh_endcap_deno_DY = new RooDataHist("rh_endcap_deno_DY", "RooHist_endcap_deno_DY", iso_deno, h_endcap_MC_deno[_DY_Full]);
    RooDataHist *rh_endcap_deno_ttbar = new RooDataHist("rh_endcap_deno_ttbar", "RooHist_endcap_deno_ttbar", iso_deno, h_endcap_MC_deno[_ttbar]);
    RooDataHist *rh_endcap_deno_tW = new RooDataHist("rh_endcap_deno_tW", "RooHist_endcap_deno_tW", iso_deno, h_endcap_MC_deno[_tW]);
    RooDataHist *rh_endcap_deno_tbarW = new RooDataHist("rh_endcap_deno_tbarW", "RooHist_endcap_deno_tbarW", iso_deno, h_endcap_MC_deno[_tbarW]);
    RooDataHist *rh_endcap_deno_WW = new RooDataHist("rh_endcap_deno_WW", "RooHist_endcap_deno_WW", iso_deno, h_endcap_MC_deno[_WW]);
    RooDataHist *rh_endcap_deno_WZ = new RooDataHist("rh_endcap_deno_WZ", "RooHist_endcap_deno_WZ", iso_deno, h_endcap_MC_deno[_WZ]);
    RooDataHist *rh_endcap_deno_ZZ = new RooDataHist("rh_endcap_deno_ZZ", "RooHist_endcap_deno_ZZ", iso_deno, h_endcap_MC_deno[_ZZ]);
    RooDataHist *rh_endcap_deno_data = new RooDataHist("rh_endcap_deno_data", "RooHist_endcap_deno_data", iso_deno, h_endcap_data_deno);

    // Making RooHistPdf
    RooHistPdf *pdf_barrel_nume_QCD   = new RooHistPdf("pdf_barrel_nume_QCD",   "Numerator barrel MC QCD template",    iso_nume, *rh_barrel_nume_QCD, 0);
    RooHistPdf *pdf_barrel_nume_WJets = new RooHistPdf("pdf_barrel_nume_WJets", "Numerator barrel MC W+Jets template", iso_nume, *rh_barrel_nume_WJets, 0);
    RooHistPdf *pdf_barrel_nume_DY    = new RooHistPdf("pdf_barrel_nume_DY",    "Numerator barrel MC DY template",     iso_nume, *rh_barrel_nume_DY, 0);
    RooHistPdf *pdf_barrel_nume_ttbar = new RooHistPdf("pdf_barrel_nume_ttbar", "Numerator barrel MC ttbar template",  iso_nume, *rh_barrel_nume_ttbar, 0);
    RooHistPdf *pdf_barrel_nume_tW    = new RooHistPdf("pdf_barrel_nume_tW",    "Numerator barrel MC tW template",     iso_nume, *rh_barrel_nume_tW, 0);
    RooHistPdf *pdf_barrel_nume_tbarW = new RooHistPdf("pdf_barrel_nume_tbarW", "Numerator barrel MC tbarW template",  iso_nume, *rh_barrel_nume_tbarW, 0);
    RooHistPdf *pdf_barrel_nume_WW    = new RooHistPdf("pdf_barrel_nume_WW",    "Numerator barrel MC WW template",     iso_nume, *rh_barrel_nume_WW, 0);
    RooHistPdf *pdf_barrel_nume_WZ    = new RooHistPdf("pdf_barrel_nume_WZ",    "Numerator barrel MC WZ template",     iso_nume, *rh_barrel_nume_WZ, 0);
    RooHistPdf *pdf_barrel_nume_ZZ    = new RooHistPdf("pdf_barrel_nume_ZZ",    "Numerator barrel MC ZZ template",     iso_nume, *rh_barrel_nume_ZZ, 0);

    RooHistPdf *pdf_endcap_nume_QCD   = new RooHistPdf("pdf_endcap_nume_QCD",   "Numerator endcap MC QCD template",    iso_nume, *rh_endcap_nume_QCD, 0);
    RooHistPdf *pdf_endcap_nume_WJets = new RooHistPdf("pdf_endcap_nume_WJets", "Numerator endcap MC W+Jets template", iso_nume, *rh_endcap_nume_WJets, 0);
    RooHistPdf *pdf_endcap_nume_DY    = new RooHistPdf("pdf_endcap_nume_DY",    "Numerator endcap MC DY template",     iso_nume, *rh_endcap_nume_DY, 0);
    RooHistPdf *pdf_endcap_nume_ttbar = new RooHistPdf("pdf_endcap_nume_ttbar", "Numerator endcap MC ttbar template",  iso_nume, *rh_endcap_nume_ttbar, 0);
    RooHistPdf *pdf_endcap_nume_tW    = new RooHistPdf("pdf_endcap_nume_tW",    "Numerator endcap MC tW template",     iso_nume, *rh_endcap_nume_tW, 0);
    RooHistPdf *pdf_endcap_nume_tbarW = new RooHistPdf("pdf_endcap_nume_tbarW", "Numerator endcap MC tbarW template",  iso_nume, *rh_endcap_nume_tbarW, 0);
    RooHistPdf *pdf_endcap_nume_WW    = new RooHistPdf("pdf_endcap_nume_WW",    "Numerator endcap MC WW template",     iso_nume, *rh_endcap_nume_WW, 0);
    RooHistPdf *pdf_endcap_nume_WZ    = new RooHistPdf("pdf_endcap_nume_WZ",    "Numerator endcap MC WZ template",     iso_nume, *rh_endcap_nume_WZ, 0);
    RooHistPdf *pdf_endcap_nume_ZZ    = new RooHistPdf("pdf_endcap_nume_ZZ",    "Numerator endcap MC ZZ template",     iso_nume, *rh_endcap_nume_ZZ, 0);

    RooHistPdf *pdf_barrel_deno_QCD   = new RooHistPdf("pdf_barrel_deno_QCD",   "Denominator barrel MC QCD template",    iso_deno, *rh_barrel_deno_QCD, 0);
    RooHistPdf *pdf_barrel_deno_WJets = new RooHistPdf("pdf_barrel_deno_WJets", "Denominator barrel MC W+Jets template", iso_deno, *rh_barrel_deno_WJets, 0);
    RooHistPdf *pdf_barrel_deno_DY    = new RooHistPdf("pdf_barrel_deno_DY",    "Denominator barrel MC DY template",     iso_deno, *rh_barrel_deno_DY, 0);
    RooHistPdf *pdf_barrel_deno_ttbar = new RooHistPdf("pdf_barrel_deno_ttbar", "Denominator barrel MC ttbar template",  iso_deno, *rh_barrel_deno_ttbar, 0);
    RooHistPdf *pdf_barrel_deno_tW    = new RooHistPdf("pdf_barrel_deno_tW",    "Denominator barrel MC tW template",     iso_deno, *rh_barrel_deno_tW, 0);
    RooHistPdf *pdf_barrel_deno_tbarW = new RooHistPdf("pdf_barrel_deno_tbarW", "Denominator barrel MC tbarW template",  iso_deno, *rh_barrel_deno_tbarW, 0);
    RooHistPdf *pdf_barrel_deno_WW    = new RooHistPdf("pdf_barrel_deno_WW",    "Denominator barrel MC WW template",     iso_deno, *rh_barrel_deno_WW, 0);
    RooHistPdf *pdf_barrel_deno_WZ    = new RooHistPdf("pdf_barrel_deno_WZ",    "Denominator barrel MC WZ template",     iso_deno, *rh_barrel_deno_WZ, 0);
    RooHistPdf *pdf_barrel_deno_ZZ    = new RooHistPdf("pdf_barrel_deno_ZZ",    "Denominator barrel MC ZZ template",     iso_deno, *rh_barrel_deno_ZZ, 0);

    RooHistPdf *pdf_endcap_deno_QCD   = new RooHistPdf("pdf_endcap_deno_QCD",   "Denominator endcap MC QCD template",    iso_deno, *rh_endcap_deno_QCD, 0);
    RooHistPdf *pdf_endcap_deno_WJets = new RooHistPdf("pdf_endcap_deno_WJets", "Denominator endcap MC W+Jets template", iso_deno, *rh_endcap_deno_WJets, 0);
    RooHistPdf *pdf_endcap_deno_DY    = new RooHistPdf("pdf_endcap_deno_DY",    "Denominator endcap MC DY template",     iso_deno, *rh_endcap_deno_DY, 0);
    RooHistPdf *pdf_endcap_deno_ttbar = new RooHistPdf("pdf_endcap_deno_ttbar", "Denominator endcap MC ttbar template",  iso_deno, *rh_endcap_deno_ttbar, 0);
    RooHistPdf *pdf_endcap_deno_tW    = new RooHistPdf("pdf_endcap_deno_tW",    "Denominator endcap MC tW template",     iso_deno, *rh_endcap_deno_tW, 0);
    RooHistPdf *pdf_endcap_deno_tbarW = new RooHistPdf("pdf_endcap_deno_tbarW", "Denominator endcap MC tbarW template",  iso_deno, *rh_endcap_deno_tbarW, 0);
    RooHistPdf *pdf_endcap_deno_WW    = new RooHistPdf("pdf_endcap_deno_WW",    "Denominator endcap MC WW template",     iso_deno, *rh_endcap_deno_WW, 0);
    RooHistPdf *pdf_endcap_deno_WZ    = new RooHistPdf("pdf_endcap_deno_WZ",    "Denominator endcap MC WZ template",     iso_deno, *rh_endcap_deno_WZ, 0);
    RooHistPdf *pdf_endcap_deno_ZZ    = new RooHistPdf("pdf_endcap_deno_ZZ",    "Denominator endcap MC ZZ template",     iso_deno, *rh_endcap_deno_ZZ, 0);

    // Constraints for integrals
//    Double_t N_barrel_nume_ttbar = (L_B2H * 831.76 / 154948878) * h_barrel_MC_nume[_ttbar]->Integral();
//    Double_t N_endcap_nume_ttbar = (L_B2H * 831.76 / 154948878) * h_endcap_MC_nume[_ttbar]->Integral();
//    Double_t N_barrel_deno_ttbar = (L_B2H * 831.76 / 154948878) * h_barrel_MC_deno[_ttbar]->Integral();
//    Double_t N_endcap_deno_ttbar = (L_B2H * 831.76 / 154948878) * h_endcap_MC_deno[_ttbar]->Integral();
    Double_t N_barrel_nume_ttbar = h_barrel_MC_nume[_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar = h_endcap_MC_nume[_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar = h_barrel_MC_deno[_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar = h_endcap_MC_deno[_ttbar]->Integral();

//    Double_t n_barrel_nume_WJets = (L_B2H * 61526 / 162787688) * h_barrel_MC_nume[_WJets]->Integral();
//    Double_t n_endcap_nume_WJets = (L_B2H * 61526 / 162787688) * h_endcap_MC_nume[_WJets]->Integral();
//    Double_t n_barrel_deno_WJets = (L_B2H * 61526 / 162787688) * h_barrel_MC_deno[_WJets]->Integral();
//    Double_t n_endcap_deno_WJets = (L_B2H * 61526 / 162787688) * h_endcap_MC_deno[_WJets]->Integral();
    Double_t N_barrel_nume_WJets = h_barrel_MC_nume[_WJets]->Integral();
    Double_t N_endcap_nume_WJets = h_endcap_MC_nume[_WJets]->Integral();
    Double_t N_barrel_deno_WJets = h_barrel_MC_deno[_WJets]->Integral();
    Double_t N_endcap_deno_WJets = h_endcap_MC_deno[_WJets]->Integral();

//    Double_t N_barrel_nume_DY = (L_B2H * 3 * 1952.68432327 / 81780984) * h_barrel_MC_nume[_DY_50to100]->Integral(); // UPDATE THE NUMBER!!
//    Double_t N_endcap_nume_DY = (L_B2H * 3 * 1952.68432327 / 81780984) * h_endcap_MC_nume[_DY_50to100]->Integral(); // THIS IS ONLY 50toINF
//    Double_t N_barrel_deno_DY = (L_B2H * 3 * 1952.68432327 / 81780984) * h_barrel_MC_deno[_DY_50to100]->Integral();
//    Double_t N_endcap_deno_DY = (L_B2H * 3 * 1952.68432327 / 81780984) * h_endcap_MC_deno[_DY_50to100]->Integral();
    Double_t N_barrel_nume_DY = h_barrel_MC_nume[_DY_50to100]->Integral();
    Double_t N_endcap_nume_DY = h_endcap_MC_nume[_DY_50to100]->Integral();
    Double_t N_barrel_deno_DY = h_barrel_MC_deno[_DY_50to100]->Integral();
    Double_t N_endcap_deno_DY = h_endcap_MC_deno[_DY_50to100]->Integral();

//    FileMgr Mgr;
//    Double_t xsec_QCD = 0, wsum_QCD = 0;
//    for (Process_t pr=_QCDMuEnriched_15to20; pr<=_QCDMuEnriched_1000toInf; pr=next(pr))
//    {
//        Mgr.SetProc(pr);
//        xsec_QCD += Mgr.Xsec[0];
//        wsum_QCD += Mgr.Wsum[0];
//    }
//    Double_t N_barrel_nume_QCD = (L_B2H * xsec_QCD / wsum_QCD) * h_barrel_MC_nume[_QCDMuEnriched_Full]->Integral();
//    Double_t N_endcap_nume_QCD = (L_B2H * xsec_QCD / wsum_QCD) * h_endcap_MC_nume[_QCDMuEnriched_Full]->Integral();
//    Double_t N_barrel_deno_QCD = (L_B2H * xsec_QCD / wsum_QCD) * h_barrel_MC_deno[_QCDMuEnriched_Full]->Integral();
//    Double_t N_endcap_deno_QCD = (L_B2H * xsec_QCD / wsum_QCD) * h_endcap_MC_deno[_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_nume_QCD = h_barrel_MC_nume[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD = h_endcap_MC_nume[_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD = h_barrel_MC_deno[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD = h_endcap_MC_deno[_QCDMuEnriched_Full]->Integral();

//    Mgr.SetProc(_tW);
//    Double_t N_barrel_nume_tW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_nume[_tW]->Integral();
//    Double_t N_endcap_nume_tW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_nume[_tW]->Integral();
//    Double_t N_barrel_deno_tW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_deno[_tW]->Integral();
//    Double_t N_endcap_deno_tW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_deno[_tW]->Integral();
    Double_t N_barrel_nume_tW = h_barrel_MC_nume[_tW]->Integral();
    Double_t N_endcap_nume_tW = h_endcap_MC_nume[_tW]->Integral();
    Double_t N_barrel_deno_tW = h_barrel_MC_deno[_tW]->Integral();
    Double_t N_endcap_deno_tW = h_endcap_MC_deno[_tW]->Integral();

//    Mgr.SetProc(_tbarW);
//    Double_t N_barrel_nume_tbarW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_nume[_tbarW]->Integral();
//    Double_t N_endcap_nume_tbarW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_nume[_tbarW]->Integral();
//    Double_t N_barrel_deno_tbarW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_deno[_tbarW]->Integral();
//    Double_t N_endcap_deno_tbarW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_deno[_tbarW]->Integral();
    Double_t N_barrel_nume_tbarW = h_barrel_MC_nume[_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW = h_endcap_MC_nume[_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW = h_barrel_MC_deno[_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW = h_endcap_MC_deno[_tbarW]->Integral();

//    Mgr.SetProc(_WW);
//    Double_t N_barrel_nume_WW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_nume[_WW]->Integral();
//    Double_t N_endcap_nume_WW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_nume[_WW]->Integral();
//    Double_t N_barrel_deno_WW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_deno[_WW]->Integral();
//    Double_t N_endcap_deno_WW = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_deno[_WW]->Integral();
    Double_t N_barrel_nume_WW = h_barrel_MC_nume[_WW]->Integral();
    Double_t N_endcap_nume_WW = h_endcap_MC_nume[_WW]->Integral();
    Double_t N_barrel_deno_WW = h_barrel_MC_deno[_WW]->Integral();
    Double_t N_endcap_deno_WW = h_endcap_MC_deno[_WW]->Integral();

//    Mgr.SetProc(_WZ);
//    Double_t N_barrel_nume_WZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_nume[_WZ]->Integral();
//    Double_t N_endcap_nume_WZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_nume[_WZ]->Integral();
//    Double_t N_barrel_deno_WZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_deno[_WZ]->Integral();
//    Double_t N_endcap_deno_WZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_deno[_WZ]->Integral();
    Double_t N_barrel_nume_WZ = h_barrel_MC_nume[_WZ]->Integral();
    Double_t N_endcap_nume_WZ = h_endcap_MC_nume[_WZ]->Integral();
    Double_t N_barrel_deno_WZ = h_barrel_MC_deno[_WZ]->Integral();
    Double_t N_endcap_deno_WZ = h_endcap_MC_deno[_WZ]->Integral();

//    Mgr.SetProc(_ZZ);
//    Double_t N_barrel_nume_ZZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_nume[_ZZ]->Integral();
//    Double_t N_endcap_nume_ZZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_nume[_ZZ]->Integral();
//    Double_t N_barrel_deno_ZZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_barrel_MC_deno[_ZZ]->Integral();
//    Double_t N_endcap_deno_ZZ = (L_B2H * Mgr.Xsec[0] / Mgr.Wsum[0]) * h_endcap_MC_deno[_ZZ]->Integral();
    Double_t N_barrel_nume_ZZ = h_barrel_MC_nume[_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ = h_endcap_MC_nume[_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ = h_barrel_MC_deno[_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ = h_endcap_MC_deno[_ZZ]->Integral();

    Double_t N_barrel_nume_total = N_barrel_nume_ttbar + N_barrel_nume_WJets + N_barrel_nume_DY + N_barrel_nume_QCD +
                                   N_barrel_nume_tW    + N_barrel_nume_tbarW + N_barrel_nume_WW + N_barrel_nume_WZ  + N_barrel_nume_ZZ;
    Double_t N_endcap_nume_total = N_endcap_nume_ttbar + N_endcap_nume_WJets + N_endcap_nume_DY + N_endcap_nume_QCD +
                                   N_endcap_nume_tW    + N_endcap_nume_tbarW + N_endcap_nume_WW + N_endcap_nume_WZ  + N_endcap_nume_ZZ;
    Double_t N_barrel_deno_total = N_barrel_deno_ttbar + N_barrel_deno_WJets + N_barrel_deno_DY + N_barrel_deno_QCD +
                                   N_barrel_deno_tW    + N_barrel_deno_tbarW + N_barrel_deno_WW + N_barrel_deno_WZ  + N_barrel_deno_ZZ;
    Double_t N_endcap_deno_total = N_endcap_deno_ttbar + N_endcap_deno_WJets + N_endcap_deno_DY + N_endcap_deno_QCD +
                                   N_endcap_deno_tW    + N_endcap_deno_tbarW + N_endcap_deno_WW + N_endcap_deno_WZ  + N_endcap_deno_ZZ;

    Double_t Nnorm_barrel_nume_ttbar = N_barrel_nume_ttbar * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_ttbar = N_endcap_nume_ttbar * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_ttbar = N_barrel_deno_ttbar * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_ttbar = N_endcap_deno_ttbar * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WJets = N_barrel_nume_WJets * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WJets = N_endcap_nume_WJets * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WJets = N_barrel_deno_WJets * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WJets = N_endcap_deno_WJets * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_DY = N_barrel_nume_DY * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_DY = N_endcap_nume_DY * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_DY = N_barrel_deno_DY * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_DY = N_endcap_deno_DY * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_tW = N_barrel_nume_tW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_tW = N_endcap_nume_tW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_tW = N_barrel_deno_tW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_tW = N_endcap_deno_tW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_tbarW = N_barrel_nume_tbarW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_tbarW = N_endcap_nume_tbarW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_tbarW = N_barrel_deno_tbarW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_tbarW = N_endcap_deno_tbarW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WW = N_barrel_nume_WW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WW = N_endcap_nume_WW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WW = N_barrel_deno_WW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WW = N_endcap_deno_WW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WZ = N_barrel_nume_WZ * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WZ = N_endcap_nume_WZ * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WZ = N_barrel_deno_WZ * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WZ = N_endcap_deno_WZ * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_ZZ = N_barrel_nume_ZZ * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_ZZ = N_endcap_nume_ZZ * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_ZZ = N_barrel_deno_ZZ * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_ZZ = N_endcap_deno_ZZ * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_QCD = N_barrel_nume_QCD * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_QCD = N_endcap_nume_QCD * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_QCD = N_barrel_deno_QCD * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_QCD = N_endcap_deno_QCD * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    // Fit constraints
//    RooRealVar n_barrel_nume_ttbar("n_barrel_nume_ttbar", "n_barrel_nume_ttbar", Nnorm_barrel_nume_ttbar, Nnorm_barrel_nume_ttbar*0.95, Nnorm_barrel_nume_ttbar*1.05);
//    RooRealVar n_endcap_nume_ttbar("n_endcap_nume_ttbar", "n_endcap_nume_ttbar", Nnorm_endcap_nume_ttbar, Nnorm_endcap_nume_ttbar*0.95, Nnorm_endcap_nume_ttbar*1.05);
//    RooRealVar n_barrel_deno_ttbar("n_barrel_deno_ttbar", "n_barrel_deno_ttbar", Nnorm_barrel_deno_ttbar, Nnorm_barrel_deno_ttbar*0.95, Nnorm_barrel_deno_ttbar*1.05);
//    RooRealVar n_endcap_deno_ttbar("n_endcap_deno_ttbar", "n_endcap_deno_ttbar", Nnorm_endcap_deno_ttbar, Nnorm_endcap_deno_ttbar*0.95, Nnorm_endcap_deno_ttbar*1.05);

//    RooRealVar n_barrel_nume_WJets("n_barrel_nume_WJets", "n_barrel_nume_WJets", Nnorm_barrel_nume_WJets, Nnorm_barrel_nume_WJets*0.95, Nnorm_barrel_nume_WJets*1.05);
//    RooRealVar n_endcap_nume_WJets("n_endcap_nume_WJets", "n_endcap_nume_WJets", Nnorm_endcap_nume_WJets, Nnorm_endcap_nume_WJets*0.95, Nnorm_endcap_nume_WJets*1.05);
//    RooRealVar n_barrel_deno_WJets("n_barrel_deno_WJets", "n_barrel_deno_WJets", Nnorm_barrel_deno_WJets, Nnorm_barrel_deno_WJets*0.95, Nnorm_barrel_deno_WJets*1.05);
//    RooRealVar n_endcap_deno_WJets("n_endcap_deno_WJets", "n_endcap_deno_WJets", Nnorm_endcap_deno_WJets, Nnorm_endcap_deno_WJets*0.95, Nnorm_endcap_deno_WJets*1.05);

//    RooRealVar n_barrel_nume_DY("n_barrel_nume_DY", "n_barrel_nume_DY", Nnorm_barrel_nume_DY, Nnorm_barrel_nume_DY*0.95, Nnorm_barrel_nume_DY*1.05);
//    RooRealVar n_endcap_nume_DY("n_endcap_nume_DY", "n_endcap_nume_DY", Nnorm_endcap_nume_DY, Nnorm_endcap_nume_DY*0.95, Nnorm_endcap_nume_DY*1.05);
//    RooRealVar n_barrel_deno_DY("n_barrel_deno_DY", "n_barrel_deno_DY", Nnorm_barrel_deno_DY, Nnorm_barrel_deno_DY*0.95, Nnorm_barrel_deno_DY*1.05);
//    RooRealVar n_endcap_deno_DY("n_endcap_deno_DY", "n_endcap_deno_DY", Nnorm_endcap_deno_DY, Nnorm_endcap_deno_DY*0.95, Nnorm_endcap_deno_DY*1.05);

//    RooRealVar n_barrel_nume_tW("n_barrel_nume_tW", "n_barrel_nume_tW", Nnorm_barrel_nume_tW, Nnorm_barrel_nume_tW*0.95, Nnorm_barrel_nume_tW*1.05);
//    RooRealVar n_endcap_nume_tW("n_endcap_nume_tW", "n_endcap_nume_tW", Nnorm_endcap_nume_tW, Nnorm_endcap_nume_tW*0.95, Nnorm_endcap_nume_tW*1.05);
//    RooRealVar n_barrel_deno_tW("n_barrel_deno_tW", "n_barrel_deno_tW", Nnorm_barrel_deno_tW, Nnorm_barrel_deno_tW*0.95, Nnorm_barrel_deno_tW*1.05);
//    RooRealVar n_endcap_deno_tW("n_endcap_deno_tW", "n_endcap_deno_tW", Nnorm_endcap_deno_tW, Nnorm_endcap_deno_tW*0.95, Nnorm_endcap_deno_tW*1.05);

//    RooRealVar n_barrel_nume_tbarW("n_barrel_nume_tbarW", "n_barrel_nume_tbarW", Nnorm_barrel_nume_tbarW, Nnorm_barrel_nume_tbarW*0.95, Nnorm_barrel_nume_tbarW*1.05);
//    RooRealVar n_endcap_nume_tbarW("n_endcap_nume_tbarW", "n_endcap_nume_tbarW", Nnorm_endcap_nume_tbarW, Nnorm_endcap_nume_tbarW*0.95, Nnorm_endcap_nume_tbarW*1.05);
//    RooRealVar n_barrel_deno_tbarW("n_barrel_deno_tbarW", "n_barrel_deno_tbarW", Nnorm_barrel_deno_tbarW, Nnorm_barrel_deno_tbarW*0.95, Nnorm_barrel_deno_tbarW*1.05);
//    RooRealVar n_endcap_deno_tbarW("n_endcap_deno_tbarW", "n_endcap_deno_tbarW", Nnorm_endcap_deno_tbarW, Nnorm_endcap_deno_tbarW*0.95, Nnorm_endcap_deno_tbarW*1.05);

//    RooRealVar n_barrel_nume_WW("n_barrel_nume_WW", "n_barrel_nume_WW", Nnorm_barrel_nume_WW, Nnorm_barrel_nume_WW*0.95, Nnorm_barrel_nume_WW*1.05);
//    RooRealVar n_endcap_nume_WW("n_endcap_nume_WW", "n_endcap_nume_WW", Nnorm_endcap_nume_WW, Nnorm_endcap_nume_WW*0.95, Nnorm_endcap_nume_WW*1.05);
//    RooRealVar n_barrel_deno_WW("n_barrel_deno_WW", "n_barrel_deno_WW", Nnorm_barrel_deno_WW, Nnorm_barrel_deno_WW*0.95, Nnorm_barrel_deno_WW*1.05);
//    RooRealVar n_endcap_deno_WW("n_endcap_deno_WW", "n_endcap_deno_WW", Nnorm_endcap_deno_WW, Nnorm_endcap_deno_WW*0.95, Nnorm_endcap_deno_WW*1.05);

//    RooRealVar n_barrel_nume_WZ("n_barrel_nume_WZ", "n_barrel_nume_WZ", Nnorm_barrel_nume_WZ, Nnorm_barrel_nume_WZ*0.95, Nnorm_barrel_nume_WZ*1.05);
//    RooRealVar n_endcap_nume_WZ("n_endcap_nume_WZ", "n_endcap_nume_WZ", Nnorm_endcap_nume_WZ, Nnorm_endcap_nume_WZ*0.95, Nnorm_endcap_nume_WZ*1.05);
//    RooRealVar n_barrel_deno_WZ("n_barrel_deno_WZ", "n_barrel_deno_WZ", Nnorm_barrel_deno_WZ, Nnorm_barrel_deno_WZ*0.95, Nnorm_barrel_deno_WZ*1.05);
//    RooRealVar n_endcap_deno_WZ("n_endcap_deno_WZ", "n_endcap_deno_WZ", Nnorm_endcap_deno_WZ, Nnorm_endcap_deno_WZ*0.95, Nnorm_endcap_deno_WZ*1.05);

//    RooRealVar n_barrel_nume_ZZ("n_barrel_nume_ZZ", "n_barrel_nume_ZZ", Nnorm_barrel_nume_ZZ, Nnorm_barrel_nume_ZZ*0.95, Nnorm_barrel_nume_ZZ*1.05);
//    RooRealVar n_endcap_nume_ZZ("n_endcap_nume_ZZ", "n_endcap_nume_ZZ", Nnorm_endcap_nume_ZZ, Nnorm_endcap_nume_ZZ*0.95, Nnorm_endcap_nume_ZZ*1.05);
//    RooRealVar n_barrel_deno_ZZ("n_barrel_deno_ZZ", "n_barrel_deno_ZZ", Nnorm_barrel_deno_ZZ, Nnorm_barrel_deno_ZZ*0.95, Nnorm_barrel_deno_ZZ*1.05);
//    RooRealVar n_endcap_deno_ZZ("n_endcap_deno_ZZ", "n_endcap_deno_ZZ", Nnorm_endcap_deno_ZZ, Nnorm_endcap_deno_ZZ*0.95, Nnorm_endcap_deno_ZZ*1.05);

    RooRealVar n_barrel_nume_ttbar("n_barrel_nume_ttbar", "n_barrel_nume_ttbar", Nnorm_barrel_nume_ttbar, Nnorm_barrel_nume_ttbar*0.9, Nnorm_barrel_nume_ttbar*1.1);
    RooRealVar n_endcap_nume_ttbar("n_endcap_nume_ttbar", "n_endcap_nume_ttbar", Nnorm_endcap_nume_ttbar, Nnorm_endcap_nume_ttbar*0.9, Nnorm_endcap_nume_ttbar*1.1);
    RooRealVar n_barrel_deno_ttbar("n_barrel_deno_ttbar", "n_barrel_deno_ttbar", Nnorm_barrel_deno_ttbar, Nnorm_barrel_deno_ttbar*0.9, Nnorm_barrel_deno_ttbar*1.1);
    RooRealVar n_endcap_deno_ttbar("n_endcap_deno_ttbar", "n_endcap_deno_ttbar", Nnorm_endcap_deno_ttbar, Nnorm_endcap_deno_ttbar*0.9, Nnorm_endcap_deno_ttbar*1.1);

    RooRealVar n_barrel_nume_WJets("n_barrel_nume_WJets", "n_barrel_nume_WJets", Nnorm_barrel_nume_WJets, Nnorm_barrel_nume_WJets*0.9, Nnorm_barrel_nume_WJets*1.1);
    RooRealVar n_endcap_nume_WJets("n_endcap_nume_WJets", "n_endcap_nume_WJets", Nnorm_endcap_nume_WJets, Nnorm_endcap_nume_WJets*0.9, Nnorm_endcap_nume_WJets*1.1);
    RooRealVar n_barrel_deno_WJets("n_barrel_deno_WJets", "n_barrel_deno_WJets", Nnorm_barrel_deno_WJets, Nnorm_barrel_deno_WJets*0.9, Nnorm_barrel_deno_WJets*1.1);
    RooRealVar n_endcap_deno_WJets("n_endcap_deno_WJets", "n_endcap_deno_WJets", Nnorm_endcap_deno_WJets, Nnorm_endcap_deno_WJets*0.9, Nnorm_endcap_deno_WJets*1.1);

    RooRealVar n_barrel_nume_DY("n_barrel_nume_DY", "n_barrel_nume_DY", Nnorm_barrel_nume_DY, Nnorm_barrel_nume_DY*0.9, Nnorm_barrel_nume_DY*1.1);
    RooRealVar n_endcap_nume_DY("n_endcap_nume_DY", "n_endcap_nume_DY", Nnorm_endcap_nume_DY, Nnorm_endcap_nume_DY*0.9, Nnorm_endcap_nume_DY*1.1);
    RooRealVar n_barrel_deno_DY("n_barrel_deno_DY", "n_barrel_deno_DY", Nnorm_barrel_deno_DY, Nnorm_barrel_deno_DY*0.9, Nnorm_barrel_deno_DY*1.1);
    RooRealVar n_endcap_deno_DY("n_endcap_deno_DY", "n_endcap_deno_DY", Nnorm_endcap_deno_DY, Nnorm_endcap_deno_DY*0.9, Nnorm_endcap_deno_DY*1.1);

    RooRealVar n_barrel_nume_tW("n_barrel_nume_tW", "n_barrel_nume_tW", Nnorm_barrel_nume_tW, Nnorm_barrel_nume_tW*0.9, Nnorm_barrel_nume_tW*1.1);
    RooRealVar n_endcap_nume_tW("n_endcap_nume_tW", "n_endcap_nume_tW", Nnorm_endcap_nume_tW, Nnorm_endcap_nume_tW*0.9, Nnorm_endcap_nume_tW*1.1);
    RooRealVar n_barrel_deno_tW("n_barrel_deno_tW", "n_barrel_deno_tW", Nnorm_barrel_deno_tW, Nnorm_barrel_deno_tW*0.9, Nnorm_barrel_deno_tW*1.1);
    RooRealVar n_endcap_deno_tW("n_endcap_deno_tW", "n_endcap_deno_tW", Nnorm_endcap_deno_tW, Nnorm_endcap_deno_tW*0.9, Nnorm_endcap_deno_tW*1.1);

    RooRealVar n_barrel_nume_tbarW("n_barrel_nume_tbarW", "n_barrel_nume_tbarW", Nnorm_barrel_nume_tbarW, Nnorm_barrel_nume_tbarW*0.9, Nnorm_barrel_nume_tbarW*1.1);
    RooRealVar n_endcap_nume_tbarW("n_endcap_nume_tbarW", "n_endcap_nume_tbarW", Nnorm_endcap_nume_tbarW, Nnorm_endcap_nume_tbarW*0.9, Nnorm_endcap_nume_tbarW*1.1);
    RooRealVar n_barrel_deno_tbarW("n_barrel_deno_tbarW", "n_barrel_deno_tbarW", Nnorm_barrel_deno_tbarW, Nnorm_barrel_deno_tbarW*0.9, Nnorm_barrel_deno_tbarW*1.1);
    RooRealVar n_endcap_deno_tbarW("n_endcap_deno_tbarW", "n_endcap_deno_tbarW", Nnorm_endcap_deno_tbarW, Nnorm_endcap_deno_tbarW*0.9, Nnorm_endcap_deno_tbarW*1.1);

    RooRealVar n_barrel_nume_WW("n_barrel_nume_WW", "n_barrel_nume_WW", Nnorm_barrel_nume_WW, Nnorm_barrel_nume_WW*0.9, Nnorm_barrel_nume_WW*1.1);
    RooRealVar n_endcap_nume_WW("n_endcap_nume_WW", "n_endcap_nume_WW", Nnorm_endcap_nume_WW, Nnorm_endcap_nume_WW*0.9, Nnorm_endcap_nume_WW*1.1);
    RooRealVar n_barrel_deno_WW("n_barrel_deno_WW", "n_barrel_deno_WW", Nnorm_barrel_deno_WW, Nnorm_barrel_deno_WW*0.9, Nnorm_barrel_deno_WW*1.1);
    RooRealVar n_endcap_deno_WW("n_endcap_deno_WW", "n_endcap_deno_WW", Nnorm_endcap_deno_WW, Nnorm_endcap_deno_WW*0.9, Nnorm_endcap_deno_WW*1.1);

    RooRealVar n_barrel_nume_WZ("n_barrel_nume_WZ", "n_barrel_nume_WZ", Nnorm_barrel_nume_WZ, Nnorm_barrel_nume_WZ*0.9, Nnorm_barrel_nume_WZ*1.1);
    RooRealVar n_endcap_nume_WZ("n_endcap_nume_WZ", "n_endcap_nume_WZ", Nnorm_endcap_nume_WZ, Nnorm_endcap_nume_WZ*0.9, Nnorm_endcap_nume_WZ*1.1);
    RooRealVar n_barrel_deno_WZ("n_barrel_deno_WZ", "n_barrel_deno_WZ", Nnorm_barrel_deno_WZ, Nnorm_barrel_deno_WZ*0.9, Nnorm_barrel_deno_WZ*1.1);
    RooRealVar n_endcap_deno_WZ("n_endcap_deno_WZ", "n_endcap_deno_WZ", Nnorm_endcap_deno_WZ, Nnorm_endcap_deno_WZ*0.9, Nnorm_endcap_deno_WZ*1.1);

    RooRealVar n_barrel_nume_ZZ("n_barrel_nume_ZZ", "n_barrel_nume_ZZ", Nnorm_barrel_nume_ZZ, Nnorm_barrel_nume_ZZ*0.9, Nnorm_barrel_nume_ZZ*1.1);
    RooRealVar n_endcap_nume_ZZ("n_endcap_nume_ZZ", "n_endcap_nume_ZZ", Nnorm_endcap_nume_ZZ, Nnorm_endcap_nume_ZZ*0.9, Nnorm_endcap_nume_ZZ*1.1);
    RooRealVar n_barrel_deno_ZZ("n_barrel_deno_ZZ", "n_barrel_deno_ZZ", Nnorm_barrel_deno_ZZ, Nnorm_barrel_deno_ZZ*0.9, Nnorm_barrel_deno_ZZ*1.1);
    RooRealVar n_endcap_deno_ZZ("n_endcap_deno_ZZ", "n_endcap_deno_ZZ", Nnorm_endcap_deno_ZZ, Nnorm_endcap_deno_ZZ*0.9, Nnorm_endcap_deno_ZZ*1.1);

    RooRealVar n_barrel_nume_QCD("n_barrel_nume_QCD", "n_barrel_nume_QCD", Nnorm_barrel_nume_QCD, Nnorm_barrel_nume_QCD*0.5, Nnorm_barrel_nume_QCD*1.5);
    RooRealVar n_endcap_nume_QCD("n_endcap_nume_QCD", "n_endcap_nume_QCD", Nnorm_endcap_nume_QCD, Nnorm_endcap_nume_QCD*0.5, Nnorm_endcap_nume_QCD*1.5);
    RooRealVar n_barrel_deno_QCD("n_barrel_deno_QCD", "n_barrel_deno_QCD", Nnorm_barrel_deno_QCD, Nnorm_barrel_deno_QCD*0.5, Nnorm_barrel_deno_QCD*1.5);
    RooRealVar n_endcap_deno_QCD("n_endcap_deno_QCD", "n_endcap_deno_QCD", Nnorm_endcap_deno_QCD, Nnorm_endcap_deno_QCD*0.5, Nnorm_endcap_deno_QCD*1.5);

    // Models
    RooAddPdf model_barrel_nume("model_barrel_nume", "model_barrel_nume", RooArgList(*pdf_barrel_nume_QCD,   *pdf_barrel_nume_WJets, *pdf_barrel_nume_DY,
                                                                                     *pdf_barrel_nume_ttbar, *pdf_barrel_nume_tW,    *pdf_barrel_nume_tbarW,
                                                                                     *pdf_barrel_nume_WW,    *pdf_barrel_nume_WZ,    *pdf_barrel_nume_ZZ),
                                RooArgList(n_barrel_nume_QCD,   n_barrel_nume_WJets, n_barrel_nume_DY, n_barrel_nume_ttbar, n_barrel_nume_tW,
                                           n_barrel_nume_tbarW, n_barrel_nume_WW,    n_barrel_nume_WZ, n_barrel_nume_ZZ));

    RooAddPdf model_endcap_nume("model_endcap_nume", "model_endcap_nume", RooArgList(*pdf_endcap_nume_QCD,   *pdf_endcap_nume_WJets, *pdf_endcap_nume_DY,
                                                                                     *pdf_endcap_nume_ttbar, *pdf_endcap_nume_tW,    *pdf_endcap_nume_tbarW,
                                                                                     *pdf_endcap_nume_WW,    *pdf_endcap_nume_WZ,    *pdf_endcap_nume_ZZ),
                                RooArgList(n_endcap_nume_QCD,   n_endcap_nume_WJets, n_endcap_nume_DY, n_endcap_nume_ttbar, n_endcap_nume_tW,
                                           n_endcap_nume_tbarW, n_endcap_nume_WW,    n_endcap_nume_WZ, n_endcap_nume_ZZ));

    RooAddPdf model_barrel_deno("model_barrel_deno", "model_barrel_deno", RooArgList(*pdf_barrel_deno_QCD,   *pdf_barrel_deno_WJets, *pdf_barrel_deno_DY,
                                                                                     *pdf_barrel_deno_ttbar, *pdf_barrel_deno_tW,    *pdf_barrel_deno_tbarW,
                                                                                     *pdf_barrel_deno_WW,    *pdf_barrel_deno_WZ,    *pdf_barrel_deno_ZZ),
                                RooArgList(n_barrel_deno_QCD,   n_barrel_deno_WJets, n_barrel_deno_DY, n_barrel_deno_ttbar, n_barrel_deno_tW,
                                           n_barrel_deno_tbarW, n_barrel_deno_WW,    n_barrel_deno_WZ, n_barrel_deno_ZZ));

    RooAddPdf model_endcap_deno("model_endcap_deno", "model_endcap_deno", RooArgList(*pdf_endcap_deno_QCD,   *pdf_endcap_deno_WJets, *pdf_endcap_deno_DY,
                                                                                     *pdf_endcap_deno_ttbar, *pdf_endcap_deno_tW,    *pdf_endcap_deno_tbarW,
                                                                                     *pdf_endcap_deno_WW,    *pdf_endcap_deno_WZ,    *pdf_endcap_deno_ZZ),
                                RooArgList(n_endcap_deno_QCD,   n_endcap_deno_WJets, n_endcap_deno_DY, n_endcap_deno_ttbar, n_endcap_deno_tW,
                                           n_endcap_deno_tbarW, n_endcap_deno_WW,    n_endcap_deno_WZ, n_endcap_deno_ZZ));

    // Fitting
    RooFitResult* fit_barrel_nume = model_barrel_nume.fitTo(*rh_barrel_nume_data, Save());
    RooFitResult* fit_endcap_nume = model_endcap_nume.fitTo(*rh_endcap_nume_data, Save());
    RooFitResult* fit_barrel_deno = model_barrel_deno.fitTo(*rh_barrel_deno_data, Save());
    RooFitResult* fit_endcap_deno = model_endcap_deno.fitTo(*rh_endcap_deno_data, Save());

    // DRAWING NUMERATOR BARREL
    TCanvas *c_fit_barrel_nume = new TCanvas("c_fit_barrel_nume", "c_fit_barrel_nume", 800, 800);
    c_fit_barrel_nume->cd();

    //Top Pad
    TPad *c1_barrel_nume = new TPad("padc1_barrel_nume","padc1_barrel_nume",0.01,0.01,0.99,0.99);
    c1_barrel_nume->Draw();
    c1_barrel_nume->cd();
    c1_barrel_nume->SetTopMargin(0.01);
    c1_barrel_nume->SetBottomMargin(0.35);
    c1_barrel_nume->SetRightMargin(0.03);
    c1_barrel_nume->SetLeftMargin(0.13);
    c1_barrel_nume->SetFillStyle(1);
    c1_barrel_nume->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_nume = iso_nume.frame(Title(" "));
    rh_barrel_nume_data->plotOn(frame_barrel_nume, DataError(RooAbsData::SumW2));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW,"
                                                           "pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,pdf_barrel_nume_DY,"
                                                           "pdf_barrel_nume_WJets,pdf_barrel_nume_QCD"),
                             LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW,"
                                                           "pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,pdf_barrel_nume_DY,"
                                                           "pdf_barrel_nume_WJets"),
                             LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW,"
                                                           "pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,pdf_barrel_nume_DY"),
                             LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW,"
                                                           "pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar"),
                             LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW,"
                                                           "pdf_barrel_nume_tbarW"),
                             LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,pdf_barrel_nume_tW"),
                             LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW"),
                             LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ"),
                             LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ"),
                             LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_nume_data->plotOn(frame_barrel_nume, DataError(RooAbsData::SumW2));
    frame_barrel_nume->Draw();
    fit_barrel_nume->Print();

    // Legend
    TLegend *legend = new TLegend(0.65, 0.7, 0.95, 0.97);
    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->AddEntry(frame_barrel_nume->nameOf(0), "Data", "LP");
    legend->AddEntry(frame_barrel_nume->nameOf(1), "#font[12]{#scale[1.1]{QCD}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(2), "#font[12]{#scale[1.1]{W}}+Jets", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(3), "DY", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(4), "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(5), "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(6), "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(7), "#font[12]{#scale[1.1]{WW}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(8), "#font[12]{#scale[1.1]{WZ}}", "F");
    legend->AddEntry(frame_barrel_nume->nameOf(9), "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "F");
    legend->SetNColumns(2);

    legend->Draw();

    frame_barrel_nume->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_nume->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_nume->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    //Top Pad
    TPad *c2_barrel_nume = new TPad("padc2_barrel_nume","padc2_barrel_nume",0.01,0.01,0.99,0.35);
    c2_barrel_nume->Draw();
    c2_barrel_nume->cd();
    c2_barrel_nume->SetTopMargin(0.05);
    c2_barrel_nume->SetBottomMargin(0.33);
    c2_barrel_nume->SetRightMargin(0.02);
    c2_barrel_nume->SetLeftMargin(0.12);
    c2_barrel_nume->SetFillStyle(0);
    c2_barrel_nume->SetGrid();

    // Ratio plot
    TH1D *h_barrel_nume_MC_fit = ((TH1D*)(model_barrel_nume.createHistogram("h_barrel_nume_MC_fit", iso_nume)));
    Double_t N_barrel_nume_data = h_barrel_data_nume->Integral();
    Double_t N_barrel_nume_MC   = h_barrel_nume_MC_fit->Integral();
    h_barrel_nume_MC_fit->Scale(N_barrel_nume_data/N_barrel_nume_MC); // Why would I wanna do that???
    cout << "----- NUMERATOR BARREL -----" << endl;
    cout << "Data integral: "   << N_barrel_nume_data << endl;
    cout << "MC integral: "     << h_barrel_nume_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_nume->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_nume_MC_fit->GetBinContent(1) << endl;

    TH1D *h_barrel_nume_ratio = ((TH1D*)(h_barrel_data_nume->Clone("h_barrel_nume_ratio")));
    h_barrel_data_nume->Sumw2(); h_barrel_nume_MC_fit->Sumw2();
    h_barrel_nume_ratio->Divide(h_barrel_data_nume, h_barrel_nume_MC_fit);
    h_barrel_nume_ratio->SetTitle("");
    h_barrel_nume_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_nume_ratio->GetXaxis()->SetNoExponent(1);
    h_barrel_nume_ratio->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})");
    h_barrel_nume_ratio->GetXaxis()->SetTitleSize(0.17);
    h_barrel_nume_ratio->GetXaxis()->SetLabelSize(0.125);
    h_barrel_nume_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_nume_ratio->GetYaxis()->SetTitle("Data/MC");
    h_barrel_nume_ratio->GetYaxis()->SetTitleSize(0.114);
    h_barrel_nume_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_nume_ratio->GetYaxis()->SetLabelSize(0.11);
    h_barrel_nume_ratio->GetYaxis()->SetTickLength(0.01);
    h_barrel_nume_ratio->GetYaxis()->SetDecimals(1);
    h_barrel_nume_ratio->SetMaximum(1.25);
    h_barrel_nume_ratio->SetMinimum(0.75);
    h_barrel_nume_ratio->GetYaxis()->SetNdivisions(5);
    h_barrel_nume_ratio->SetLineWidth(1);
    h_barrel_nume_ratio->SetLineColor(kBlack);
    h_barrel_nume_ratio->SetMarkerStyle(kFullDotLarge);
    h_barrel_nume_ratio->SetMarkerColor(kBlack);
    h_barrel_nume_ratio->SetStats(kFALSE);

    h_barrel_nume_ratio->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line = ((TH1D*)(h_barrel_nume_ratio->Clone("h_line")));
    h_line->Reset("ICES");
    for (Int_t i=1; i<=h_line->GetNbinsX(); i++)
        h_line->SetBinContent(i, 1);
    h_line->SetLineColor(kRed);
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_nume = model_barrel_nume.createChi2(*rh_barrel_nume_data);
    cout << "chi2: " << chi2_barrel_nume->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_nume->getVal() / ((Double_t)h_barrel_data_nume->GetNbinsX()) << endl;

    // DRAWING NUMERATOR ENDCAP
    TCanvas *c_fit_endcap_nume = new TCanvas("c_fit_endcap_nume", "c_fit_endcap_nume", 800, 800);
    c_fit_endcap_nume->cd();

    //Top Pad
    TPad *c1_endcap_nume = new TPad("padc1_endcap_nume","padc1_endcap_nume",0.01,0.01,0.99,0.99);
    c1_endcap_nume->Draw();
    c1_endcap_nume->cd();
    c1_endcap_nume->SetTopMargin(0.01);
    c1_endcap_nume->SetBottomMargin(0.35);
    c1_endcap_nume->SetRightMargin(0.03);
    c1_endcap_nume->SetLeftMargin(0.13);
    c1_endcap_nume->SetFillStyle(1);
    c1_endcap_nume->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_nume = iso_nume.frame(Title(" "));
    rh_endcap_nume_data->plotOn(frame_endcap_nume, DataError(RooAbsData::SumW2));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW,"
                                                           "pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,pdf_endcap_nume_DY,"
                                                           "pdf_endcap_nume_WJets,pdf_endcap_nume_QCD"),
                             LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW,"
                                                           "pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,pdf_endcap_nume_DY,"
                                                           "pdf_endcap_nume_WJets"),
                             LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW,"
                                                           "pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,pdf_endcap_nume_DY"),
                             LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW,"
                                                           "pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar"),
                             LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW,"
                                                           "pdf_endcap_nume_tbarW"),
                             LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,pdf_endcap_nume_tW"),
                             LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW"),
                             LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ"),
                             LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ"),
                             LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_nume_data->plotOn(frame_endcap_nume, DataError(RooAbsData::SumW2));
    frame_endcap_nume->Draw();
    fit_endcap_nume->Print();

    // Legend
    legend->Draw();

    frame_endcap_nume->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_nume->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_nume->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    //Top Pad
    TPad *c2_endcap_nume = new TPad("padc2_endcap_nume","padc2_endcap_nume",0.01,0.01,0.99,0.35);
    c2_endcap_nume->Draw();
    c2_endcap_nume->cd();
    c2_endcap_nume->SetTopMargin(0.05);
    c2_endcap_nume->SetBottomMargin(0.33);
    c2_endcap_nume->SetRightMargin(0.02);
    c2_endcap_nume->SetLeftMargin(0.12);
    c2_endcap_nume->SetFillStyle(0);
    c2_endcap_nume->SetGrid();

    // Ratio plot
    TH1D *h_endcap_nume_MC_fit = ((TH1D*)(model_endcap_nume.createHistogram("h_endcap_nume_MC_fit", iso_nume)));
    Double_t N_endcap_nume_data = h_endcap_data_nume->Integral();
    Double_t N_endcap_nume_MC   = h_endcap_nume_MC_fit->Integral();
    h_endcap_nume_MC_fit->Scale(N_endcap_nume_data/N_endcap_nume_MC); // Why would I wanna do that???
    cout << "----- NUMERATOR ENDCAP -----" << endl;
    cout << "Data integral: "   << N_endcap_nume_data << endl;
    cout << "MC integral: "     << h_endcap_nume_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_nume->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_nume_MC_fit->GetBinContent(1) << endl;

    TH1D *h_endcap_nume_ratio = ((TH1D*)(h_endcap_data_nume->Clone("h_endcap_nume_ratio")));
    h_endcap_data_nume->Sumw2(); h_endcap_nume_MC_fit->Sumw2();
    h_endcap_nume_ratio->Divide(h_endcap_data_nume, h_endcap_nume_MC_fit);
    h_endcap_nume_ratio->SetTitle("");
    h_endcap_nume_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_nume_ratio->GetXaxis()->SetNoExponent(1);
    h_endcap_nume_ratio->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})");
    h_endcap_nume_ratio->GetXaxis()->SetTitleSize(0.17);
    h_endcap_nume_ratio->GetXaxis()->SetLabelSize(0.125);
    h_endcap_nume_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_nume_ratio->GetYaxis()->SetTitle("Data/MC");
    h_endcap_nume_ratio->GetYaxis()->SetTitleSize(0.114);
    h_endcap_nume_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_nume_ratio->GetYaxis()->SetLabelSize(0.11);
    h_endcap_nume_ratio->GetYaxis()->SetTickLength(0.01);
    h_endcap_nume_ratio->GetYaxis()->SetDecimals(1);
    h_endcap_nume_ratio->SetMaximum(1.25);
    h_endcap_nume_ratio->SetMinimum(0.75);
    h_endcap_nume_ratio->GetYaxis()->SetNdivisions(5);
    h_endcap_nume_ratio->SetLineWidth(1);
    h_endcap_nume_ratio->SetLineColor(kBlack);
    h_endcap_nume_ratio->SetMarkerStyle(kFullDotLarge);
    h_endcap_nume_ratio->SetMarkerColor(kBlack);
    h_endcap_nume_ratio->SetStats(kFALSE);

    h_endcap_nume_ratio->Draw("E1P");

    // Red line at Data/MC=1
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_nume = model_endcap_nume.createChi2(*rh_endcap_nume_data);
    cout << "chi2: " << chi2_endcap_nume->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_nume->getVal() / ((Double_t)h_endcap_data_nume->GetNbinsX()) << endl;

    // DRAWING DENOMINATOR BARREL
    TCanvas *c_fit_barrel_deno = new TCanvas("c_fit_barrel_deno", "c_fit_barrel_deno", 800, 800);
    c_fit_barrel_deno->cd();

    //Top Pad
    TPad *c1_barrel_deno = new TPad("padc1_barrel_deno","padc1_barrel_deno",0.01,0.01,0.99,0.99);
    c1_barrel_deno->Draw();
    c1_barrel_deno->cd();
    c1_barrel_deno->SetTopMargin(0.01);
    c1_barrel_deno->SetBottomMargin(0.35);
    c1_barrel_deno->SetRightMargin(0.03);
    c1_barrel_deno->SetLeftMargin(0.13);
    c1_barrel_deno->SetFillStyle(1);
    c1_barrel_deno->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_deno = iso_deno.frame(Title(" "));
    rh_barrel_deno_data->plotOn(frame_barrel_deno, DataError(RooAbsData::SumW2));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW,"
                                                           "pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,pdf_barrel_deno_DY,"
                                                           "pdf_barrel_deno_WJets,pdf_barrel_deno_QCD"),
                             LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW,"
                                                           "pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,pdf_barrel_deno_DY,"
                                                           "pdf_barrel_deno_WJets"),
                             LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW,"
                                                           "pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,pdf_barrel_deno_DY"),
                             LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW,"
                                                           "pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar"),
                             LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW,"
                                                           "pdf_barrel_deno_tbarW"),
                             LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,pdf_barrel_deno_tW"),
                             LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW"),
                             LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ"),
                             LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ"),
                             LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_deno_data->plotOn(frame_barrel_deno, DataError(RooAbsData::SumW2));
    frame_barrel_deno->Draw();
    fit_barrel_deno->Print();

    // Legend
    legend->Draw();

    frame_barrel_deno->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_deno->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_deno->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    //Top Pad
    TPad *c2_barrel_deno = new TPad("padc2_barrel_deno","padc2_barrel_deno",0.01,0.01,0.99,0.35);
    c2_barrel_deno->Draw();
    c2_barrel_deno->cd();
    c2_barrel_deno->SetTopMargin(0.05);
    c2_barrel_deno->SetBottomMargin(0.33);
    c2_barrel_deno->SetRightMargin(0.02);
    c2_barrel_deno->SetLeftMargin(0.12);
    c2_barrel_deno->SetFillStyle(0);
    c2_barrel_deno->SetGrid();

    // Ratio plot
    TH1D *h_barrel_deno_MC_fit = ((TH1D*)(model_barrel_deno.createHistogram("h_barrel_deno_MC_fit", iso_deno)));
    Double_t N_barrel_deno_data = h_barrel_data_deno->Integral();
    Double_t N_barrel_deno_MC   = h_barrel_deno_MC_fit->Integral();
    h_barrel_deno_MC_fit->Scale(N_barrel_deno_data/N_barrel_deno_MC); // Why is this necessary???
    cout << "----- DENOMINATOR BARREL -----" << endl;
    cout << "Data integral: "   << N_barrel_deno_data << endl;
    cout << "MC integral: "     << h_barrel_deno_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_deno->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_deno_MC_fit->GetBinContent(1) << endl;

    TH1D *h_barrel_deno_ratio = ((TH1D*)(h_barrel_data_deno->Clone("h_barrel_deno_ratio")));
    h_barrel_data_deno->Sumw2(); h_barrel_deno_MC_fit->Sumw2();
    h_barrel_deno_ratio->Divide(h_barrel_data_deno, h_barrel_deno_MC_fit);
    h_barrel_deno_ratio->SetTitle("");
    h_barrel_deno_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_deno_ratio->GetXaxis()->SetNoExponent(1);
    h_barrel_deno_ratio->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_barrel_deno_ratio->GetXaxis()->SetTitleSize(0.17);
    h_barrel_deno_ratio->GetXaxis()->SetLabelSize(0.125);
    h_barrel_deno_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_deno_ratio->GetYaxis()->SetTitle("Data/MC");
    h_barrel_deno_ratio->GetYaxis()->SetTitleSize(0.114);
    h_barrel_deno_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_deno_ratio->GetYaxis()->SetLabelSize(0.11);
    h_barrel_deno_ratio->GetYaxis()->SetTickLength(0.01);
    h_barrel_deno_ratio->GetYaxis()->SetDecimals(1);
    h_barrel_deno_ratio->SetMaximum(1.25);
    h_barrel_deno_ratio->SetMinimum(0.75);
    h_barrel_deno_ratio->GetYaxis()->SetNdivisions(5);
    h_barrel_deno_ratio->SetLineWidth(1);
    h_barrel_deno_ratio->SetLineColor(kBlack);
    h_barrel_deno_ratio->SetMarkerStyle(kFullDotLarge);
    h_barrel_deno_ratio->SetMarkerColor(kBlack);
    h_barrel_deno_ratio->SetStats(kFALSE);

    h_barrel_deno_ratio->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line_deno = ((TH1D*)(h_barrel_deno_ratio->Clone("h_line_deno")));
    h_line_deno->Reset("ICES");
    for (Int_t i=1; i<=h_line_deno->GetNbinsX(); i++)
        h_line_deno->SetBinContent(i, 1);
    h_line_deno->SetLineColor(kRed);
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_deno = model_barrel_deno.createChi2(*rh_barrel_deno_data);
    cout << "chi2: " << chi2_barrel_deno->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_deno->getVal() / ((Double_t)h_barrel_data_deno->GetNbinsX()) << endl;

    // DRAWING DENOMINATOR ENDCAP
    TCanvas *c_fit_endcap_deno = new TCanvas("c_fit_endcap_deno", "c_fit_endcap_deno", 800, 800);
    c_fit_endcap_deno->cd();

    //Top Pad
    TPad *c1_endcap_deno = new TPad("padc1_endcap_deno","padc1_endcap_deno",0.01,0.01,0.99,0.99);
    c1_endcap_deno->Draw();
    c1_endcap_deno->cd();
    c1_endcap_deno->SetTopMargin(0.01);
    c1_endcap_deno->SetBottomMargin(0.35);
    c1_endcap_deno->SetRightMargin(0.03);
    c1_endcap_deno->SetLeftMargin(0.13);
    c1_endcap_deno->SetFillStyle(1);
    c1_endcap_deno->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_deno = iso_deno.frame(Title(" "));
    rh_endcap_deno_data->plotOn(frame_endcap_deno, DataError(RooAbsData::SumW2));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW,"
                                                           "pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,pdf_endcap_deno_DY,"
                                                           "pdf_endcap_deno_WJets,pdf_endcap_deno_QCD"),
                             LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW,"
                                                           "pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,pdf_endcap_deno_DY,"
                                                           "pdf_endcap_deno_WJets"),
                             LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW,"
                                                           "pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,pdf_endcap_deno_DY"),
                             LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW,"
                                                           "pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar"),
                             LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW,"
                                                           "pdf_endcap_deno_tbarW"),
                             LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,pdf_endcap_deno_tW"),
                             LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW"),
                             LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ"),
                             LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ"),
                             LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_deno_data->plotOn(frame_endcap_deno, DataError(RooAbsData::SumW2));
    frame_endcap_deno->Draw();
    fit_endcap_deno->Print();

    // Legend
    legend->Draw();

    frame_endcap_deno->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_deno->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_deno->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    //Top Pad
    TPad *c2_endcap_deno = new TPad("padc2_endcap_deno","padc2_endcap_deno",0.01,0.01,0.99,0.35);
    c2_endcap_deno->Draw();
    c2_endcap_deno->cd();
    c2_endcap_deno->SetTopMargin(0.05);
    c2_endcap_deno->SetBottomMargin(0.33);
    c2_endcap_deno->SetRightMargin(0.02);
    c2_endcap_deno->SetLeftMargin(0.12);
    c2_endcap_deno->SetFillStyle(0);
    c2_endcap_deno->SetGrid();

    // Ratio plot
    TH1D *h_endcap_deno_MC_fit = ((TH1D*)(model_endcap_deno.createHistogram("h_endcap_deno_MC_fit", iso_deno)));
    Double_t N_endcap_deno_data = h_endcap_data_deno->Integral();
    Double_t N_endcap_deno_MC   = h_endcap_deno_MC_fit->Integral();
    h_endcap_deno_MC_fit->Scale(N_endcap_deno_data/N_endcap_deno_MC); // Why is this necessary???
    cout << "----- DENOMINATOR ENDCAP -----" << endl;
    cout << "Data integral: "   << N_endcap_deno_data << endl;
    cout << "MC integral: "     << h_endcap_deno_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_deno->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_deno_MC_fit->GetBinContent(1) << endl;

    TH1D *h_endcap_deno_ratio = ((TH1D*)(h_endcap_data_deno->Clone("h_endcap_deno_ratio")));
    h_endcap_data_deno->Sumw2(); h_endcap_deno_MC_fit->Sumw2();
    h_endcap_deno_ratio->Divide(h_endcap_data_deno, h_endcap_deno_MC_fit);
    h_endcap_deno_ratio->SetTitle("");
    h_endcap_deno_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_deno_ratio->GetXaxis()->SetNoExponent(1);
    h_endcap_deno_ratio->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_endcap_deno_ratio->GetXaxis()->SetTitleSize(0.17);
    h_endcap_deno_ratio->GetXaxis()->SetLabelSize(0.125);
    h_endcap_deno_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_deno_ratio->GetYaxis()->SetTitle("Data/MC");
    h_endcap_deno_ratio->GetYaxis()->SetTitleSize(0.114);
    h_endcap_deno_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_deno_ratio->GetYaxis()->SetLabelSize(0.11);
    h_endcap_deno_ratio->GetYaxis()->SetTickLength(0.01);
    h_endcap_deno_ratio->GetYaxis()->SetDecimals(1);
    h_endcap_deno_ratio->SetMaximum(1.25);
    h_endcap_deno_ratio->SetMinimum(0.75);
    h_endcap_deno_ratio->GetYaxis()->SetNdivisions(5);
    h_endcap_deno_ratio->SetLineWidth(1);
    h_endcap_deno_ratio->SetLineColor(kBlack);
    h_endcap_deno_ratio->SetMarkerStyle(kFullDotLarge);
    h_endcap_deno_ratio->SetMarkerColor(kBlack);
    h_endcap_deno_ratio->SetStats(kFALSE);

    h_endcap_deno_ratio->Draw("E1P");

    // Red line at Data/MC=1
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_deno = model_endcap_deno.createChi2(*rh_endcap_deno_data);
    cout << "chi2: " << chi2_endcap_deno->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_deno->getVal() / ((Double_t)h_endcap_data_deno->GetNbinsX()) << endl;


//    // Writing
//    TFile *file_FR = new TFile("/media/sf_DATA/FR/Muon/FakeRate_muon.root", "RECREATE");
//    if (file_FR->IsOpen()) cout << "File '/media/sf_DATA/FR/Muon/FakeRate_muon.root' has been created. Writing histograms.." << endl;
//    file_FR->cd();
//    h_FRratio_barrel->Write();
//    h_FRratio_endcap->Write();
//    h_FRtemplate_barrel->Write();
//    h_FRtemplate_endcap->Write();
//    cout << "Finished. Closing the file.." << endl;
//    file_FR->Close();
//    if (!file_FR->IsOpen()) cout << "File '/media/sf_DATA/FR/Muon/FakeRate_muon.root' has been closed successfully." << endl;
//    else cout << "File did not close!" << endl;

//    // Drawing
//    TCanvas *c_FR_barrel = new TCanvas("c_FR_barrel", "c_FR_barrel", 800, 800);
//    c_FR_barrel->cd();
//    h_FRratio_barrel->Draw();
//    h_FRtemplate_barrel->Draw("same");
//    c_FR_barrel->Update();
//    TCanvas *c_FR_endcap = new TCanvas("c_FR_endcap", "c_FR_endcap", 800, 800);
//    c_FR_endcap->cd();
//    h_FRratio_endcap->Draw();
//    h_FRtemplate_endcap->Draw("same");
//    c_FR_endcap->Update();

} // End of Mu_EstFR()


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
