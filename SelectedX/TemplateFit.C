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

void TemplateFit (TString WhichX = "", Int_t type = 2)
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
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr1]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr1]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr1]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr1]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr1]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr1]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr1]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr1]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr1]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr1]);
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr1]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr1]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr1]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr1]);

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

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

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

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno_100to500[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno_100to500[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

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
            file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_data_deno_50to70 );
            file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_data_deno_50to70 );
            file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_data_nume_50to70 );
            file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_data_nume_50to70 );
            file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_data_deno_70to100);
            file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_data_deno_70to100);
            file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_data_nume_70to100);
            file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_data_nume_70to100);
            file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_data_deno_100to500);
            file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_data_deno_100to500);
            file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_data_nume_100to500);
            file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_data_nume_100to500);

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
            file->GetObject("h_PFiso_barrel_deno_50to70",   h_temp[0]);
            file->GetObject("h_PFiso_endcap_deno_50to70",   h_temp[1]);
            file->GetObject("h_PFiso_barrel_nume_50to70",   h_temp[2]);
            file->GetObject("h_PFiso_endcap_nume_50to70",   h_temp[3]);
            file->GetObject("h_PFiso_barrel_deno_70to100",  h_temp[4]);
            file->GetObject("h_PFiso_endcap_deno_70to100",  h_temp[5]);
            file->GetObject("h_PFiso_barrel_nume_70to100",  h_temp[6]);
            file->GetObject("h_PFiso_endcap_nume_70to100",  h_temp[7]);
            file->GetObject("h_PFiso_barrel_deno_100to500", h_temp[8]);
            file->GetObject("h_PFiso_endcap_deno_100to500", h_temp[9]);
            file->GetObject("h_PFiso_barrel_nume_100to500", h_temp[10]);
            file->GetObject("h_PFiso_endcap_nume_100to500", h_temp[11]);

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

// ######################## MODEL BUILDING ##########################

    // Making RooDataHist
    RooRealVar iso_nume("iso", "PFiso/p_{T}", 0, 0.15);
    RooRealVar iso_deno("iso", "PFiso/p_{T}", 0, 5);

    RooDataHist *rh_barrel_nume_QCD_50to70   = new RooDataHist("rh_barrel_nume_QCD_50to70",   "RooHist_barrel_nume_QCD_50to70",   iso_nume, h_barrel_MC_nume_50to70 [_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_nume_QCD_50to70   = new RooDataHist("rh_endcap_nume_QCD_50to70",   "RooHist_endcap_nume_QCD_50to70",   iso_nume, h_endcap_MC_nume_50to70 [_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_deno_QCD_50to70   = new RooDataHist("rh_barrel_deno_QCD_50to70",   "RooHist_barrel_deno_QCD_50to70",   iso_deno, h_barrel_MC_deno_50to70 [_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_deno_QCD_50to70   = new RooDataHist("rh_endcap_deno_QCD_50to70",   "RooHist_endcap_deno_QCD_50to70",   iso_deno, h_endcap_MC_deno_50to70 [_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_nume_QCD_70to100  = new RooDataHist("rh_barrel_nume_QCD_70to100",  "RooHist_barrel_nume_QCD_70to100",  iso_nume, h_barrel_MC_nume_70to100[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_nume_QCD_70to100  = new RooDataHist("rh_endcap_nume_QCD_70to100",  "RooHist_endcap_nume_QCD_70to100",  iso_nume, h_endcap_MC_nume_70to100[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_deno_QCD_70to100  = new RooDataHist("rh_barrel_deno_QCD_70to100",  "RooHist_barrel_deno_QCD_70to100",  iso_deno, h_barrel_MC_deno_70to100[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_deno_QCD_70to100  = new RooDataHist("rh_endcap_deno_QCD_70to100",  "RooHist_endcap_deno_QCD_70to100",  iso_deno, h_endcap_MC_deno_70to100[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_nume_QCD_100to500 = new RooDataHist("rh_barrel_nume_QCD_100to500", "RooHist_barrel_nume_QCD_100to500", iso_nume, h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_nume_QCD_100to500 = new RooDataHist("rh_endcap_nume_QCD_100to500", "RooHist_endcap_nume_QCD_100to500", iso_nume, h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_deno_QCD_100to500 = new RooDataHist("rh_barrel_deno_QCD_100to500", "RooHist_barrel_deno_QCD_100to500", iso_deno, h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_deno_QCD_100to500 = new RooDataHist("rh_endcap_deno_QCD_100to500", "RooHist_endcap_deno_QCD_100to500", iso_deno, h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]);

    RooDataHist *rh_barrel_nume_WJets_50to70   = new RooDataHist("rh_barrel_nume_WJets_50to70",   "RooHist_barrel_nume_WJets_50to70",   iso_nume, h_barrel_MC_nume_50to70 [_WJets]);
    RooDataHist *rh_endcap_nume_WJets_50to70   = new RooDataHist("rh_endcap_nume_WJets_50to70",   "RooHist_endcap_nume_WJets_50to70",   iso_nume, h_endcap_MC_nume_50to70 [_WJets]);
    RooDataHist *rh_barrel_deno_WJets_50to70   = new RooDataHist("rh_barrel_deno_WJets_50to70",   "RooHist_barrel_deno_WJets_50to70",   iso_deno, h_barrel_MC_deno_50to70 [_WJets]);
    RooDataHist *rh_endcap_deno_WJets_50to70   = new RooDataHist("rh_endcap_deno_WJets_50to70",   "RooHist_endcap_deno_WJets_50to70",   iso_deno, h_endcap_MC_deno_50to70 [_WJets]);
    RooDataHist *rh_barrel_nume_WJets_70to100  = new RooDataHist("rh_barrel_nume_WJets_70to100",  "RooHist_barrel_nume_WJets_70to100",  iso_nume, h_barrel_MC_nume_70to100[_WJets]);
    RooDataHist *rh_endcap_nume_WJets_70to100  = new RooDataHist("rh_endcap_nume_WJets_70to100",  "RooHist_endcap_nume_WJets_70to100",  iso_nume, h_endcap_MC_nume_70to100[_WJets]);
    RooDataHist *rh_barrel_deno_WJets_70to100  = new RooDataHist("rh_barrel_deno_WJets_70to100",  "RooHist_barrel_deno_WJets_70to100",  iso_deno, h_barrel_MC_deno_70to100[_WJets]);
    RooDataHist *rh_endcap_deno_WJets_70to100  = new RooDataHist("rh_endcap_deno_WJets_70to100",  "RooHist_endcap_deno_WJets_70to100",  iso_deno, h_endcap_MC_deno_70to100[_WJets]);
    RooDataHist *rh_barrel_nume_WJets_100to500 = new RooDataHist("rh_barrel_nume_WJets_100to500", "RooHist_barrel_nume_WJets_100to500", iso_nume, h_barrel_MC_nume_100to500[_WJets]);
    RooDataHist *rh_endcap_nume_WJets_100to500 = new RooDataHist("rh_endcap_nume_WJets_100to500", "RooHist_endcap_nume_WJets_100to500", iso_nume, h_endcap_MC_nume_100to500[_WJets]);
    RooDataHist *rh_barrel_deno_WJets_100to500 = new RooDataHist("rh_barrel_deno_WJets_100to500", "RooHist_barrel_deno_WJets_100to500", iso_deno, h_barrel_MC_deno_100to500[_WJets]);
    RooDataHist *rh_endcap_deno_WJets_100to500 = new RooDataHist("rh_endcap_deno_WJets_100to500", "RooHist_endcap_deno_WJets_100to500", iso_deno, h_endcap_MC_deno_100to500[_WJets]);

    RooDataHist *rh_barrel_nume_DY_50to70   = new RooDataHist("rh_barrel_nume_DY_50to70",   "RooHist_barrel_nume_DY_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_DY_Full]);
    RooDataHist *rh_endcap_nume_DY_50to70   = new RooDataHist("rh_endcap_nume_DY_50to70",   "RooHist_endcap_nume_DY_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_DY_Full]);
    RooDataHist *rh_barrel_deno_DY_50to70   = new RooDataHist("rh_barrel_deno_DY_50to70",   "RooHist_barrel_deno_DY_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_DY_Full]);
    RooDataHist *rh_endcap_deno_DY_50to70   = new RooDataHist("rh_endcap_deno_DY_50to70",   "RooHist_endcap_deno_DY_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_DY_Full]);
    RooDataHist *rh_barrel_nume_DY_70to100  = new RooDataHist("rh_barrel_nume_DY_70to100",  "RooHist_barrel_nume_DY_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_DY_Full]);
    RooDataHist *rh_endcap_nume_DY_70to100  = new RooDataHist("rh_endcap_nume_DY_70to100",  "RooHist_endcap_nume_DY_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_DY_Full]);
    RooDataHist *rh_barrel_deno_DY_70to100  = new RooDataHist("rh_barrel_deno_DY_70to100",  "RooHist_barrel_deno_DY_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_DY_Full]);
    RooDataHist *rh_endcap_deno_DY_70to100  = new RooDataHist("rh_endcap_deno_DY_70to100",  "RooHist_endcap_deno_DY_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_DY_Full]);
    RooDataHist *rh_barrel_nume_DY_100to500 = new RooDataHist("rh_barrel_nume_DY_100to500", "RooHist_barrel_nume_DY_100to500", iso_nume, h_barrel_MC_nume_100to500[_DY_Full]);
    RooDataHist *rh_endcap_nume_DY_100to500 = new RooDataHist("rh_endcap_nume_DY_100to500", "RooHist_endcap_nume_DY_100to500", iso_nume, h_endcap_MC_nume_100to500[_DY_Full]);
    RooDataHist *rh_barrel_deno_DY_100to500 = new RooDataHist("rh_barrel_deno_DY_100to500", "RooHist_barrel_deno_DY_100to500", iso_deno, h_barrel_MC_deno_100to500[_DY_Full]);
    RooDataHist *rh_endcap_deno_DY_100to500 = new RooDataHist("rh_endcap_deno_DY_100to500", "RooHist_endcap_deno_DY_100to500", iso_deno, h_endcap_MC_deno_100to500[_DY_Full]);

    RooDataHist *rh_barrel_nume_ttbar_50to70   = new RooDataHist("rh_barrel_nume_ttbar_50to70",   "RooHist_barrel_nume_ttbar_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar_50to70   = new RooDataHist("rh_endcap_nume_ttbar_50to70",   "RooHist_endcap_nume_ttbar_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar_50to70   = new RooDataHist("rh_barrel_deno_ttbar_50to70",   "RooHist_barrel_deno_ttbar_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar_50to70   = new RooDataHist("rh_endcap_deno_ttbar_50to70",   "RooHist_endcap_deno_ttbar_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_ttbar]);
    RooDataHist *rh_barrel_nume_ttbar_70to100  = new RooDataHist("rh_barrel_nume_ttbar_70to100",  "RooHist_barrel_nume_ttbar_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar_70to100  = new RooDataHist("rh_endcap_nume_ttbar_70to100",  "RooHist_endcap_nume_ttbar_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar_70to100  = new RooDataHist("rh_barrel_deno_ttbar_70to100",  "RooHist_barrel_deno_ttbar_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar_70to100  = new RooDataHist("rh_endcap_deno_ttbar_70to100",  "RooHist_endcap_deno_ttbar_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_ttbar]);
    RooDataHist *rh_barrel_nume_ttbar_100to500 = new RooDataHist("rh_barrel_nume_ttbar_100to500", "RooHist_barrel_nume_ttbar_100to500", iso_nume, h_barrel_MC_nume_100to500[_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar_100to500 = new RooDataHist("rh_endcap_nume_ttbar_100to500", "RooHist_endcap_nume_ttbar_100to500", iso_nume, h_endcap_MC_nume_100to500[_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar_100to500 = new RooDataHist("rh_barrel_deno_ttbar_100to500", "RooHist_barrel_deno_ttbar_100to500", iso_deno, h_barrel_MC_deno_100to500[_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar_100to500 = new RooDataHist("rh_endcap_deno_ttbar_100to500", "RooHist_endcap_deno_ttbar_100to500", iso_deno, h_endcap_MC_deno_100to500[_ttbar]);

    RooDataHist *rh_barrel_nume_tW_50to70   = new RooDataHist("rh_barrel_nume_tW_50to70",   "RooHist_barrel_nume_tW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_tW]);
    RooDataHist *rh_endcap_nume_tW_50to70   = new RooDataHist("rh_endcap_nume_tW_50to70",   "RooHist_endcap_nume_tW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_tW]);
    RooDataHist *rh_barrel_deno_tW_50to70   = new RooDataHist("rh_barrel_deno_tW_50to70",   "RooHist_barrel_deno_tW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_tW]);
    RooDataHist *rh_endcap_deno_tW_50to70   = new RooDataHist("rh_endcap_deno_tW_50to70",   "RooHist_endcap_deno_tW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_tW]);
    RooDataHist *rh_barrel_nume_tW_70to100  = new RooDataHist("rh_barrel_nume_tW_70to100",  "RooHist_barrel_nume_tW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_tW]);
    RooDataHist *rh_endcap_nume_tW_70to100  = new RooDataHist("rh_endcap_nume_tW_70to100",  "RooHist_endcap_nume_tW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_tW]);
    RooDataHist *rh_barrel_deno_tW_70to100  = new RooDataHist("rh_barrel_deno_tW_70to100",  "RooHist_barrel_deno_tW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_tW]);
    RooDataHist *rh_endcap_deno_tW_70to100  = new RooDataHist("rh_endcap_deno_tW_70to100",  "RooHist_endcap_deno_tW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_tW]);
    RooDataHist *rh_barrel_nume_tW_100to500 = new RooDataHist("rh_barrel_nume_tW_100to500", "RooHist_barrel_nume_tW_100to500", iso_nume, h_barrel_MC_nume_100to500[_tW]);
    RooDataHist *rh_endcap_nume_tW_100to500 = new RooDataHist("rh_endcap_nume_tW_100to500", "RooHist_endcap_nume_tW_100to500", iso_nume, h_endcap_MC_nume_100to500[_tW]);
    RooDataHist *rh_barrel_deno_tW_100to500 = new RooDataHist("rh_barrel_deno_tW_100to500", "RooHist_barrel_deno_tW_100to500", iso_deno, h_barrel_MC_deno_100to500[_tW]);
    RooDataHist *rh_endcap_deno_tW_100to500 = new RooDataHist("rh_endcap_deno_tW_100to500", "RooHist_endcap_deno_tW_100to500", iso_deno, h_endcap_MC_deno_100to500[_tW]);

    RooDataHist *rh_barrel_nume_tbarW_50to70   = new RooDataHist("rh_barrel_nume_tbarW_50to70",   "RooHist_barrel_nume_tbarW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW_50to70   = new RooDataHist("rh_endcap_nume_tbarW_50to70",   "RooHist_endcap_nume_tbarW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW_50to70   = new RooDataHist("rh_barrel_deno_tbarW_50to70",   "RooHist_barrel_deno_tbarW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW_50to70   = new RooDataHist("rh_endcap_deno_tbarW_50to70",   "RooHist_endcap_deno_tbarW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_tbarW]);
    RooDataHist *rh_barrel_nume_tbarW_70to100  = new RooDataHist("rh_barrel_nume_tbarW_70to100",  "RooHist_barrel_nume_tbarW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW_70to100  = new RooDataHist("rh_endcap_nume_tbarW_70to100",  "RooHist_endcap_nume_tbarW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW_70to100  = new RooDataHist("rh_barrel_deno_tbarW_70to100",  "RooHist_barrel_deno_tbarW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW_70to100  = new RooDataHist("rh_endcap_deno_tbarW_70to100",  "RooHist_endcap_deno_tbarW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_tbarW]);
    RooDataHist *rh_barrel_nume_tbarW_100to500 = new RooDataHist("rh_barrel_nume_tbarW_100to500", "RooHist_barrel_nume_tbarW_100to500", iso_nume, h_barrel_MC_nume_100to500[_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW_100to500 = new RooDataHist("rh_endcap_nume_tbarW_100to500", "RooHist_endcap_nume_tbarW_100to500", iso_nume, h_endcap_MC_nume_100to500[_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW_100to500 = new RooDataHist("rh_barrel_deno_tbarW_100to500", "RooHist_barrel_deno_tbarW_100to500", iso_deno, h_barrel_MC_deno_100to500[_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW_100to500 = new RooDataHist("rh_endcap_deno_tbarW_100to500", "RooHist_endcap_deno_tbarW_100to500", iso_deno, h_endcap_MC_deno_100to500[_tbarW]);

    RooDataHist *rh_barrel_nume_WW_50to70  =  new RooDataHist("rh_barrel_nume_WW_50to70",   "RooHist_barrel_nume_WW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_WW]);
    RooDataHist *rh_endcap_nume_WW_50to70  =  new RooDataHist("rh_endcap_nume_WW_50to70",   "RooHist_endcap_nume_WW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_WW]);
    RooDataHist *rh_barrel_deno_WW_50to70  =  new RooDataHist("rh_barrel_deno_WW_50to70",   "RooHist_barrel_deno_WW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_WW]);
    RooDataHist *rh_endcap_deno_WW_50to70  =  new RooDataHist("rh_endcap_deno_WW_50to70",   "RooHist_endcap_deno_WW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_WW]);
    RooDataHist *rh_barrel_nume_WW_70to100 =  new RooDataHist("rh_barrel_nume_WW_70to100",  "RooHist_barrel_nume_WW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_WW]);
    RooDataHist *rh_endcap_nume_WW_70to100 =  new RooDataHist("rh_endcap_nume_WW_70to100",  "RooHist_endcap_nume_WW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_WW]);
    RooDataHist *rh_barrel_deno_WW_70to100 =  new RooDataHist("rh_barrel_deno_WW_70to100",  "RooHist_barrel_deno_WW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_WW]);
    RooDataHist *rh_endcap_deno_WW_70to100 =  new RooDataHist("rh_endcap_deno_WW_70to100",  "RooHist_endcap_deno_WW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_WW]);
    RooDataHist *rh_barrel_nume_WW_100to500 = new RooDataHist("rh_barrel_nume_WW_100to500", "RooHist_barrel_nume_WW_100to500", iso_nume, h_barrel_MC_nume_100to500[_WW]);
    RooDataHist *rh_endcap_nume_WW_100to500 = new RooDataHist("rh_endcap_nume_WW_100to500", "RooHist_endcap_nume_WW_100to500", iso_nume, h_endcap_MC_nume_100to500[_WW]);
    RooDataHist *rh_barrel_deno_WW_100to500 = new RooDataHist("rh_barrel_deno_WW_100to500", "RooHist_barrel_deno_WW_100to500", iso_deno, h_barrel_MC_deno_100to500[_WW]);
    RooDataHist *rh_endcap_deno_WW_100to500 = new RooDataHist("rh_endcap_deno_WW_100to500", "RooHist_endcap_deno_WW_100to500", iso_deno, h_endcap_MC_deno_100to500[_WW]);

    RooDataHist *rh_barrel_nume_WZ_50to70   = new RooDataHist("rh_barrel_nume_WZ_50to70",   "RooHist_barrel_nume_WZ_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_WZ]);
    RooDataHist *rh_endcap_nume_WZ_50to70   = new RooDataHist("rh_endcap_nume_WZ_50to70",   "RooHist_endcap_nume_WZ_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_WZ]);
    RooDataHist *rh_barrel_deno_WZ_50to70   = new RooDataHist("rh_barrel_deno_WZ_50to70",   "RooHist_barrel_deno_WZ_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_WZ]);
    RooDataHist *rh_endcap_deno_WZ_50to70   = new RooDataHist("rh_endcap_deno_WZ_50to70",   "RooHist_endcap_deno_WZ_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_WZ]);
    RooDataHist *rh_barrel_nume_WZ_70to100  = new RooDataHist("rh_barrel_nume_WZ_70to100",  "RooHist_barrel_nume_WZ_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_WZ]);
    RooDataHist *rh_endcap_nume_WZ_70to100  = new RooDataHist("rh_endcap_nume_WZ_70to100",  "RooHist_endcap_nume_WZ_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_WZ]);
    RooDataHist *rh_barrel_deno_WZ_70to100  = new RooDataHist("rh_barrel_deno_WZ_70to100",  "RooHist_barrel_deno_WZ_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_WZ]);
    RooDataHist *rh_endcap_deno_WZ_70to100  = new RooDataHist("rh_endcap_deno_WZ_70to100",  "RooHist_endcap_deno_WZ_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_WZ]);
    RooDataHist *rh_barrel_nume_WZ_100to500 = new RooDataHist("rh_barrel_nume_WZ_100to500", "RooHist_barrel_nume_WZ_100to500", iso_nume, h_barrel_MC_nume_100to500[_WZ]);
    RooDataHist *rh_endcap_nume_WZ_100to500 = new RooDataHist("rh_endcap_nume_WZ_100to500", "RooHist_endcap_nume_WZ_100to500", iso_nume, h_endcap_MC_nume_100to500[_WZ]);
    RooDataHist *rh_barrel_deno_WZ_100to500 = new RooDataHist("rh_barrel_deno_WZ_100to500", "RooHist_barrel_deno_WZ_100to500", iso_deno, h_barrel_MC_deno_100to500[_WZ]);
    RooDataHist *rh_endcap_deno_WZ_100to500 = new RooDataHist("rh_endcap_deno_WZ_100to500", "RooHist_endcap_deno_WZ_100to500", iso_deno, h_endcap_MC_deno_100to500[_WZ]);

    RooDataHist *rh_barrel_nume_ZZ_50to70   = new RooDataHist("rh_barrel_nume_ZZ_50to70",   "RooHist_barrel_nume_ZZ_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ_50to70   = new RooDataHist("rh_endcap_nume_ZZ_50to70",   "RooHist_endcap_nume_ZZ_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ_50to70   = new RooDataHist("rh_barrel_deno_ZZ_50to70",   "RooHist_barrel_deno_ZZ_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ_50to70   = new RooDataHist("rh_endcap_deno_ZZ_50to70",   "RooHist_endcap_deno_ZZ_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_ZZ]);
    RooDataHist *rh_barrel_nume_ZZ_70to100  = new RooDataHist("rh_barrel_nume_ZZ_70to100",  "RooHist_barrel_nume_ZZ_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ_70to100  = new RooDataHist("rh_endcap_nume_ZZ_70to100",  "RooHist_endcap_nume_ZZ_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ_70to100  = new RooDataHist("rh_barrel_deno_ZZ_70to100",  "RooHist_barrel_deno_ZZ_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ_70to100  = new RooDataHist("rh_endcap_deno_ZZ_70to100",  "RooHist_endcap_deno_ZZ_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_ZZ]);
    RooDataHist *rh_barrel_nume_ZZ_100to500 = new RooDataHist("rh_barrel_nume_ZZ_100to500", "RooHist_barrel_nume_ZZ_100to500", iso_nume, h_barrel_MC_nume_100to500[_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ_100to500 = new RooDataHist("rh_endcap_nume_ZZ_100to500", "RooHist_endcap_nume_ZZ_100to500", iso_nume, h_endcap_MC_nume_100to500[_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ_100to500 = new RooDataHist("rh_barrel_deno_ZZ_100to500", "RooHist_barrel_deno_ZZ_100to500", iso_deno, h_barrel_MC_deno_100to500[_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ_100to500 = new RooDataHist("rh_endcap_deno_ZZ_100to500", "RooHist_endcap_deno_ZZ_100to500", iso_deno, h_endcap_MC_deno_100to500[_ZZ]);

    RooDataHist *rh_barrel_nume_data_50to70   = new RooDataHist("rh_barrel_nume_data_50to70",   "RooHist_barrel_nume_data_50to70",   iso_nume, h_barrel_data_nume_50to70 );
    RooDataHist *rh_endcap_nume_data_50to70   = new RooDataHist("rh_endcap_nume_data_50to70",   "RooHist_endcap_nume_data_50to70",   iso_nume, h_endcap_data_nume_50to70 );
    RooDataHist *rh_barrel_deno_data_50to70   = new RooDataHist("rh_barrel_deno_data_50to70",   "RooHist_barrel_deno_data_50to70",   iso_deno, h_barrel_data_deno_50to70 );
    RooDataHist *rh_endcap_deno_data_50to70   = new RooDataHist("rh_endcap_deno_data_50to70",   "RooHist_endcap_deno_data_50to70",   iso_deno, h_endcap_data_deno_50to70 );
    RooDataHist *rh_barrel_nume_data_70to100  = new RooDataHist("rh_barrel_nume_data_70to100",  "RooHist_barrel_nume_data_70to100",  iso_nume, h_barrel_data_nume_70to100);
    RooDataHist *rh_endcap_nume_data_70to100  = new RooDataHist("rh_endcap_nume_data_70to100",  "RooHist_endcap_nume_data_70to100",  iso_nume, h_endcap_data_nume_70to100);
    RooDataHist *rh_barrel_deno_data_70to100  = new RooDataHist("rh_barrel_deno_data_70to100",  "RooHist_barrel_deno_data_70to100",  iso_deno, h_barrel_data_deno_70to100);
    RooDataHist *rh_endcap_deno_data_70to100  = new RooDataHist("rh_endcap_deno_data_70to100",  "RooHist_endcap_deno_data_70to100",  iso_deno, h_endcap_data_deno_70to100);
    RooDataHist *rh_barrel_nume_data_100to500 = new RooDataHist("rh_barrel_nume_data_100to500", "RooHist_barrel_nume_data_100to500", iso_nume, h_barrel_data_nume_100to500);
    RooDataHist *rh_endcap_nume_data_100to500 = new RooDataHist("rh_endcap_nume_data_100to500", "RooHist_endcap_nume_data_100to500", iso_nume, h_endcap_data_nume_100to500);
    RooDataHist *rh_barrel_deno_data_100to500 = new RooDataHist("rh_barrel_deno_data_100to500", "RooHist_barrel_deno_data_100to500", iso_deno, h_barrel_data_deno_100to500);
    RooDataHist *rh_endcap_deno_data_100to500 = new RooDataHist("rh_endcap_deno_data_100to500", "RooHist_endcap_deno_data_100to500", iso_deno, h_endcap_data_deno_100to500);

    // Making RooHistPdf
    RooHistPdf *pdf_barrel_nume_QCD_50to70   = new RooHistPdf("pdf_barrel_nume_QCD_50to70",   "Numerator barrel MC QCD template 50 to 70",     iso_nume, *rh_barrel_nume_QCD_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_QCD_50to70   = new RooHistPdf("pdf_endcap_nume_QCD_50to70",   "Numerator endcap MC QCD template 50 to 70",     iso_nume, *rh_endcap_nume_QCD_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_QCD_50to70   = new RooHistPdf("pdf_barrel_deno_QCD_50to70",   "Denominator barrel MC QCD template 50 to 70",   iso_deno, *rh_barrel_deno_QCD_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_QCD_50to70   = new RooHistPdf("pdf_endcap_deno_QCD_50to70",   "Denominator endcap MC QCD template 50 to 70",   iso_deno, *rh_endcap_deno_QCD_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_QCD_70to100  = new RooHistPdf("pdf_barrel_nume_QCD_70to100",  "Numerator barrel MC QCD template 70 to 100",    iso_nume, *rh_barrel_nume_QCD_70to100,  0);
    RooHistPdf *pdf_endcap_nume_QCD_70to100  = new RooHistPdf("pdf_endcap_nume_QCD_70to100",  "Numerator endcap MC QCD template 70 to 100",    iso_nume, *rh_endcap_nume_QCD_70to100,  0);
    RooHistPdf *pdf_barrel_deno_QCD_70to100  = new RooHistPdf("pdf_barrel_deno_QCD_70to100",  "Denominator barrel MC QCD template 70 to 100",  iso_deno, *rh_barrel_deno_QCD_70to100,  0);
    RooHistPdf *pdf_endcap_deno_QCD_70to100  = new RooHistPdf("pdf_endcap_deno_QCD_70to100",  "Denominator endcap MC QCD template 70 to 100",  iso_deno, *rh_endcap_deno_QCD_70to100,  0);
    RooHistPdf *pdf_barrel_nume_QCD_100to500 = new RooHistPdf("pdf_barrel_nume_QCD_100to500", "Numerator barrel MC QCD template 200 to 500",   iso_nume, *rh_barrel_nume_QCD_100to500, 0);
    RooHistPdf *pdf_endcap_nume_QCD_100to500 = new RooHistPdf("pdf_endcap_nume_QCD_100to500", "Numerator endcap MC QCD template 200 to 500",   iso_nume, *rh_endcap_nume_QCD_100to500, 0);
    RooHistPdf *pdf_barrel_deno_QCD_100to500 = new RooHistPdf("pdf_barrel_deno_QCD_100to500", "Denominator barrel MC QCD template 200 to 500", iso_deno, *rh_barrel_deno_QCD_100to500, 0);
    RooHistPdf *pdf_endcap_deno_QCD_100to500 = new RooHistPdf("pdf_endcap_deno_QCD_100to500", "Denominator endcap MC QCD template 200 to 500", iso_deno, *rh_endcap_deno_QCD_100to500, 0);

    RooHistPdf *pdf_barrel_nume_WJets_50to70   = new RooHistPdf("pdf_barrel_nume_WJets_50to70",   "Numerator barrel MC W+Jets template 50 to 70",     iso_nume, *rh_barrel_nume_WJets_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WJets_50to70   = new RooHistPdf("pdf_endcap_nume_WJets_50to70",   "Numerator endcap MC W+Jets template 50 to 70",     iso_nume, *rh_endcap_nume_WJets_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WJets_50to70   = new RooHistPdf("pdf_barrel_deno_WJets_50to70",   "Denominator barrel MC W+Jets template 50 to 70",   iso_deno, *rh_barrel_deno_WJets_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WJets_50to70   = new RooHistPdf("pdf_endcap_deno_WJets_50to70",   "Denominator endcap MC W+Jets template 50 to 70",   iso_deno, *rh_endcap_deno_WJets_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WJets_70to100  = new RooHistPdf("pdf_barrel_nume_WJets_70to100",  "Numerator barrel MC W+Jets template 70 to 100",    iso_nume, *rh_barrel_nume_WJets_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WJets_70to100  = new RooHistPdf("pdf_endcap_nume_WJets_70to100",  "Numerator endcap MC W+Jets template 70 to 100",    iso_nume, *rh_endcap_nume_WJets_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WJets_70to100  = new RooHistPdf("pdf_barrel_deno_WJets_70to100",  "Denominator barrel MC W+Jets template 70 to 100",  iso_deno, *rh_barrel_deno_WJets_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WJets_70to100  = new RooHistPdf("pdf_endcap_deno_WJets_70to100",  "Denominator endcap MC W+Jets template 70 to 100",  iso_deno, *rh_endcap_deno_WJets_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WJets_100to500 = new RooHistPdf("pdf_barrel_nume_WJets_100to500", "Numerator barrel MC W+Jets template 200 to 500",   iso_nume, *rh_barrel_nume_WJets_100to500, 0);
    RooHistPdf *pdf_endcap_nume_WJets_100to500 = new RooHistPdf("pdf_endcap_nume_WJets_100to500", "Numerator endcap MC W+Jets template 200 to 500",   iso_nume, *rh_endcap_nume_WJets_100to500, 0);
    RooHistPdf *pdf_barrel_deno_WJets_100to500 = new RooHistPdf("pdf_barrel_deno_WJets_100to500", "Denominator barrel MC W+Jets template 200 to 500", iso_deno, *rh_barrel_deno_WJets_100to500, 0);
    RooHistPdf *pdf_endcap_deno_WJets_100to500 = new RooHistPdf("pdf_endcap_deno_WJets_100to500", "Denominator endcap MC W+Jets template 200 to 500", iso_deno, *rh_endcap_deno_WJets_100to500, 0);

    RooHistPdf *pdf_barrel_nume_DY_50to70   = new RooHistPdf("pdf_barrel_nume_DY_50to70",   "Numerator barrel MC DY template 50 to 70",     iso_nume, *rh_barrel_nume_DY_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_DY_50to70   = new RooHistPdf("pdf_endcap_nume_DY_50to70",   "Numerator endcap MC DY template 50 to 70",     iso_nume, *rh_endcap_nume_DY_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_DY_50to70   = new RooHistPdf("pdf_barrel_deno_DY_50to70",   "Denominator barrel MC DY template 50 to 70",   iso_deno, *rh_barrel_deno_DY_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_DY_50to70   = new RooHistPdf("pdf_endcap_deno_DY_50to70",   "Denominator endcap MC DY template 50 to 70",   iso_deno, *rh_endcap_deno_DY_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_DY_70to100  = new RooHistPdf("pdf_barrel_nume_DY_70to100",  "Numerator barrel MC DY template 70 to 100",    iso_nume, *rh_barrel_nume_DY_70to100,  0);
    RooHistPdf *pdf_endcap_nume_DY_70to100  = new RooHistPdf("pdf_endcap_nume_DY_70to100",  "Numerator endcap MC DY template 70 to 100",    iso_nume, *rh_endcap_nume_DY_70to100,  0);
    RooHistPdf *pdf_barrel_deno_DY_70to100  = new RooHistPdf("pdf_barrel_deno_DY_70to100",  "Denominator barrel MC DY template 70 to 100",  iso_deno, *rh_barrel_deno_DY_70to100,  0);
    RooHistPdf *pdf_endcap_deno_DY_70to100  = new RooHistPdf("pdf_endcap_deno_DY_70to100",  "Denominator endcap MC DY template 70 to 100",  iso_deno, *rh_endcap_deno_DY_70to100,  0);
    RooHistPdf *pdf_barrel_nume_DY_100to500 = new RooHistPdf("pdf_barrel_nume_DY_100to500", "Numerator barrel MC DY template 200 to 500",   iso_nume, *rh_barrel_nume_DY_100to500, 0);
    RooHistPdf *pdf_endcap_nume_DY_100to500 = new RooHistPdf("pdf_endcap_nume_DY_100to500", "Numerator endcap MC DY template 200 to 500",   iso_nume, *rh_endcap_nume_DY_100to500, 0);
    RooHistPdf *pdf_barrel_deno_DY_100to500 = new RooHistPdf("pdf_barrel_deno_DY_100to500", "Denominator barrel MC DY template 200 to 500", iso_deno, *rh_barrel_deno_DY_100to500, 0);
    RooHistPdf *pdf_endcap_deno_DY_100to500 = new RooHistPdf("pdf_endcap_deno_DY_100to500", "Denominator endcap MC DY template 200 to 500", iso_deno, *rh_endcap_deno_DY_100to500, 0);

    RooHistPdf *pdf_barrel_nume_ttbar_50to70   = new RooHistPdf("pdf_barrel_nume_ttbar_50to70",   "Numerator barrel MC ttbar template 50 to 70",     iso_nume, *rh_barrel_nume_ttbar_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_ttbar_50to70   = new RooHistPdf("pdf_endcap_nume_ttbar_50to70",   "Numerator endcap MC ttbar template 50 to 70",     iso_nume, *rh_endcap_nume_ttbar_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_ttbar_50to70   = new RooHistPdf("pdf_barrel_deno_ttbar_50to70",   "Denominator barrel MC ttbar template 50 to 70",   iso_deno, *rh_barrel_deno_ttbar_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_ttbar_50to70   = new RooHistPdf("pdf_endcap_deno_ttbar_50to70",   "Denominator endcap MC ttbar template 50 to 70",   iso_deno, *rh_endcap_deno_ttbar_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_ttbar_70to100  = new RooHistPdf("pdf_barrel_nume_ttbar_70to100",  "Numerator barrel MC ttbar template 70 to 100",    iso_nume, *rh_barrel_nume_ttbar_70to100,  0);
    RooHistPdf *pdf_endcap_nume_ttbar_70to100  = new RooHistPdf("pdf_endcap_nume_ttbar_70to100",  "Numerator endcap MC ttbar template 70 to 100",    iso_nume, *rh_endcap_nume_ttbar_70to100,  0);
    RooHistPdf *pdf_barrel_deno_ttbar_70to100  = new RooHistPdf("pdf_barrel_deno_ttbar_70to100",  "Denominator barrel MC ttbar template 70 to 100",  iso_deno, *rh_barrel_deno_ttbar_70to100,  0);
    RooHistPdf *pdf_endcap_deno_ttbar_70to100  = new RooHistPdf("pdf_endcap_deno_ttbar_70to100",  "Denominator endcap MC ttbar template 70 to 100",  iso_deno, *rh_endcap_deno_ttbar_70to100,  0);
    RooHistPdf *pdf_barrel_nume_ttbar_100to500 = new RooHistPdf("pdf_barrel_nume_ttbar_100to500", "Numerator barrel MC ttbar template 200 to 500",   iso_nume, *rh_barrel_nume_ttbar_100to500, 0);
    RooHistPdf *pdf_endcap_nume_ttbar_100to500 = new RooHistPdf("pdf_endcap_nume_ttbar_100to500", "Numerator endcap MC ttbar template 200 to 500",   iso_nume, *rh_endcap_nume_ttbar_100to500, 0);
    RooHistPdf *pdf_barrel_deno_ttbar_100to500 = new RooHistPdf("pdf_barrel_deno_ttbar_100to500", "Denominator barrel MC ttbar template 200 to 500", iso_deno, *rh_barrel_deno_ttbar_100to500, 0);
    RooHistPdf *pdf_endcap_deno_ttbar_100to500 = new RooHistPdf("pdf_endcap_deno_ttbar_100to500", "Denominator endcap MC ttbar template 200 to 500", iso_deno, *rh_endcap_deno_ttbar_100to500, 0);

    RooHistPdf *pdf_barrel_nume_tW_50to70   = new RooHistPdf("pdf_barrel_nume_tW_50to70",   "Numerator barrel MC tW template 50 to 70",     iso_nume, *rh_barrel_nume_tW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_tW_50to70   = new RooHistPdf("pdf_endcap_nume_tW_50to70",   "Numerator endcap MC tW template 50 to 70",     iso_nume, *rh_endcap_nume_tW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_tW_50to70   = new RooHistPdf("pdf_barrel_deno_tW_50to70",   "Denominator barrel MC tW template 50 to 70",   iso_deno, *rh_barrel_deno_tW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_tW_50to70   = new RooHistPdf("pdf_endcap_deno_tW_50to70",   "Denominator endcap MC tW template 50 to 70",   iso_deno, *rh_endcap_deno_tW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_tW_70to100  = new RooHistPdf("pdf_barrel_nume_tW_70to100",  "Numerator barrel MC tW template 70 to 100",    iso_nume, *rh_barrel_nume_tW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_tW_70to100  = new RooHistPdf("pdf_endcap_nume_tW_70to100",  "Numerator endcap MC tW template 70 to 100",    iso_nume, *rh_endcap_nume_tW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_tW_70to100  = new RooHistPdf("pdf_barrel_deno_tW_70to100",  "Denominator barrel MC tW template 70 to 100",  iso_deno, *rh_barrel_deno_tW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_tW_70to100  = new RooHistPdf("pdf_endcap_deno_tW_70to100",  "Denominator endcap MC tW template 70 to 100",  iso_deno, *rh_endcap_deno_tW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_tW_100to500 = new RooHistPdf("pdf_barrel_nume_tW_100to500", "Numerator barrel MC tW template 200 to 500",   iso_nume, *rh_barrel_nume_tW_100to500, 0);
    RooHistPdf *pdf_endcap_nume_tW_100to500 = new RooHistPdf("pdf_endcap_nume_tW_100to500", "Numerator endcap MC tW template 200 to 500",   iso_nume, *rh_endcap_nume_tW_100to500, 0);
    RooHistPdf *pdf_barrel_deno_tW_100to500 = new RooHistPdf("pdf_barrel_deno_tW_100to500", "Denominator barrel MC tW template 200 to 500", iso_deno, *rh_barrel_deno_tW_100to500, 0);
    RooHistPdf *pdf_endcap_deno_tW_100to500 = new RooHistPdf("pdf_endcap_deno_tW_100to500", "Denominator endcap MC tW template 200 to 500", iso_deno, *rh_endcap_deno_tW_100to500, 0);

    RooHistPdf *pdf_barrel_nume_tbarW_50to70   = new RooHistPdf("pdf_barrel_nume_tbarW_50to70",   "Numerator barrel MC tbarW template 50 to 70",     iso_nume, *rh_barrel_nume_tbarW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_tbarW_50to70   = new RooHistPdf("pdf_endcap_nume_tbarW_50to70",   "Numerator endcap MC tbarW template 50 to 70",     iso_nume, *rh_endcap_nume_tbarW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_tbarW_50to70   = new RooHistPdf("pdf_barrel_deno_tbarW_50to70",   "Denominator barrel MC tbarW template 50 to 70",   iso_deno, *rh_barrel_deno_tbarW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_tbarW_50to70   = new RooHistPdf("pdf_endcap_deno_tbarW_50to70",   "Denominator endcap MC tbarW template 50 to 70",   iso_deno, *rh_endcap_deno_tbarW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_tbarW_70to100  = new RooHistPdf("pdf_barrel_nume_tbarW_70to100",  "Numerator barrel MC tbarW template 70 to 100",    iso_nume, *rh_barrel_nume_tbarW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_tbarW_70to100  = new RooHistPdf("pdf_endcap_nume_tbarW_70to100",  "Numerator endcap MC tbarW template 70 to 100",    iso_nume, *rh_endcap_nume_tbarW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_tbarW_70to100  = new RooHistPdf("pdf_barrel_deno_tbarW_70to100",  "Denominator barrel MC tbarW template 70 to 100",  iso_deno, *rh_barrel_deno_tbarW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_tbarW_70to100  = new RooHistPdf("pdf_endcap_deno_tbarW_70to100",  "Denominator endcap MC tbarW template 70 to 100",  iso_deno, *rh_endcap_deno_tbarW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_tbarW_100to500 = new RooHistPdf("pdf_barrel_nume_tbarW_100to500", "Numerator barrel MC tbarW template 200 to 500",   iso_nume, *rh_barrel_nume_tbarW_100to500, 0);
    RooHistPdf *pdf_endcap_nume_tbarW_100to500 = new RooHistPdf("pdf_endcap_nume_tbarW_100to500", "Numerator endcap MC tbarW template 200 to 500",   iso_nume, *rh_endcap_nume_tbarW_100to500, 0);
    RooHistPdf *pdf_barrel_deno_tbarW_100to500 = new RooHistPdf("pdf_barrel_deno_tbarW_100to500", "Denominator barrel MC tbarW template 200 to 500", iso_deno, *rh_barrel_deno_tbarW_100to500, 0);
    RooHistPdf *pdf_endcap_deno_tbarW_100to500 = new RooHistPdf("pdf_endcap_deno_tbarW_100to500", "Denominator endcap MC tbarW template 200 to 500", iso_deno, *rh_endcap_deno_tbarW_100to500, 0);

    RooHistPdf *pdf_barrel_nume_WW_50to70   = new RooHistPdf("pdf_barrel_nume_WW_50to70",   "Numerator barrel MC WW template 50 to 70",     iso_nume, *rh_barrel_nume_WW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WW_50to70   = new RooHistPdf("pdf_endcap_nume_WW_50to70",   "Numerator endcap MC WW template 50 to 70",     iso_nume, *rh_endcap_nume_WW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WW_50to70   = new RooHistPdf("pdf_barrel_deno_WW_50to70",   "Denominator barrel MC WW template 50 to 70",   iso_deno, *rh_barrel_deno_WW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WW_50to70   = new RooHistPdf("pdf_endcap_deno_WW_50to70",   "Denominator endcap MC WW template 50 to 70",   iso_deno, *rh_endcap_deno_WW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WW_70to100  = new RooHistPdf("pdf_barrel_nume_WW_70to100",  "Numerator barrel MC WW template 70 to 100",    iso_nume, *rh_barrel_nume_WW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WW_70to100  = new RooHistPdf("pdf_endcap_nume_WW_70to100",  "Numerator endcap MC WW template 70 to 100",    iso_nume, *rh_endcap_nume_WW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WW_70to100  = new RooHistPdf("pdf_barrel_deno_WW_70to100",  "Denominator barrel MC WW template 70 to 100",  iso_deno, *rh_barrel_deno_WW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WW_70to100  = new RooHistPdf("pdf_endcap_deno_WW_70to100",  "Denominator endcap MC WW template 70 to 100",  iso_deno, *rh_endcap_deno_WW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WW_100to500 = new RooHistPdf("pdf_barrel_nume_WW_100to500", "Numerator barrel MC WW template 200 to 500",   iso_nume, *rh_barrel_nume_WW_100to500, 0);
    RooHistPdf *pdf_endcap_nume_WW_100to500 = new RooHistPdf("pdf_endcap_nume_WW_100to500", "Numerator endcap MC WW template 200 to 500",   iso_nume, *rh_endcap_nume_WW_100to500, 0);
    RooHistPdf *pdf_barrel_deno_WW_100to500 = new RooHistPdf("pdf_barrel_deno_WW_100to500", "Denominator barrel MC WW template 200 to 500", iso_deno, *rh_barrel_deno_WW_100to500, 0);
    RooHistPdf *pdf_endcap_deno_WW_100to500 = new RooHistPdf("pdf_endcap_deno_WW_100to500", "Denominator endcap MC WW template 200 to 500", iso_deno, *rh_endcap_deno_WW_100to500, 0);

    RooHistPdf *pdf_barrel_nume_WZ_50to70   = new RooHistPdf("pdf_barrel_nume_WZ_50to70",   "Numerator barrel MC WZ template 50 to 70",     iso_nume, *rh_barrel_nume_WZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WZ_50to70   = new RooHistPdf("pdf_endcap_nume_WZ_50to70",   "Numerator endcap MC WZ template 50 to 70",     iso_nume, *rh_endcap_nume_WZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WZ_50to70   = new RooHistPdf("pdf_barrel_deno_WZ_50to70",   "Denominator barrel MC WZ template 50 to 70",   iso_deno, *rh_barrel_deno_WZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WZ_50to70   = new RooHistPdf("pdf_endcap_deno_WZ_50to70",   "Denominator endcap MC WZ template 50 to 70",   iso_deno, *rh_endcap_deno_WZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WZ_70to100  = new RooHistPdf("pdf_barrel_nume_WZ_70to100",  "Numerator barrel MC WZ template 70 to 100",    iso_nume, *rh_barrel_nume_WZ_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WZ_70to100  = new RooHistPdf("pdf_endcap_nume_WZ_70to100",  "Numerator endcap MC WZ template 70 to 100",    iso_nume, *rh_endcap_nume_WZ_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WZ_70to100  = new RooHistPdf("pdf_barrel_deno_WZ_70to100",  "Denominator barrel MC WZ template 70 to 100",  iso_deno, *rh_barrel_deno_WZ_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WZ_70to100  = new RooHistPdf("pdf_endcap_deno_WZ_70to100",  "Denominator endcap MC WZ template 70 to 100",  iso_deno, *rh_endcap_deno_WZ_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WZ_100to500 = new RooHistPdf("pdf_barrel_nume_WZ_100to500", "Numerator barrel MC WZ template 200 to 500",   iso_nume, *rh_barrel_nume_WZ_100to500, 0);
    RooHistPdf *pdf_endcap_nume_WZ_100to500 = new RooHistPdf("pdf_endcap_nume_WZ_100to500", "Numerator endcap MC WZ template 200 to 500",   iso_nume, *rh_endcap_nume_WZ_100to500, 0);
    RooHistPdf *pdf_barrel_deno_WZ_100to500 = new RooHistPdf("pdf_barrel_deno_WZ_100to500", "Denominator barrel MC WZ template 200 to 500", iso_deno, *rh_barrel_deno_WZ_100to500, 0);
    RooHistPdf *pdf_endcap_deno_WZ_100to500 = new RooHistPdf("pdf_endcap_deno_WZ_100to500", "Denominator endcap MC WZ template 200 to 500", iso_deno, *rh_endcap_deno_WZ_100to500, 0);

    RooHistPdf *pdf_barrel_nume_ZZ_50to70   = new RooHistPdf("pdf_barrel_nume_ZZ_50to70",   "Numerator barrel MC ZZ template 50 to 70",     iso_nume, *rh_barrel_nume_ZZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_ZZ_50to70   = new RooHistPdf("pdf_endcap_nume_ZZ_50to70",   "Numerator endcap MC ZZ template 50 to 70",     iso_nume, *rh_endcap_nume_ZZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_ZZ_50to70   = new RooHistPdf("pdf_barrel_deno_ZZ_50to70",   "Denominator barrel MC ZZ template 50 to 70",   iso_deno, *rh_barrel_deno_ZZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_ZZ_50to70   = new RooHistPdf("pdf_endcap_deno_ZZ_50to70",   "Denominator endcap MC ZZ template 50 to 70",   iso_deno, *rh_endcap_deno_ZZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_ZZ_70to100  = new RooHistPdf("pdf_barrel_nume_ZZ_70to100",  "Numerator barrel MC ZZ template 70 to 100",    iso_nume, *rh_barrel_nume_ZZ_70to100,  0);
    RooHistPdf *pdf_endcap_nume_ZZ_70to100  = new RooHistPdf("pdf_endcap_nume_ZZ_70to100",  "Numerator endcap MC ZZ template 70 to 100",    iso_nume, *rh_endcap_nume_ZZ_70to100,  0);
    RooHistPdf *pdf_barrel_deno_ZZ_70to100  = new RooHistPdf("pdf_barrel_deno_ZZ_70to100",  "Denominator barrel MC ZZ template 70 to 100",  iso_deno, *rh_barrel_deno_ZZ_70to100,  0);
    RooHistPdf *pdf_endcap_deno_ZZ_70to100  = new RooHistPdf("pdf_endcap_deno_ZZ_70to100",  "Denominator endcap MC ZZ template 70 to 100",  iso_deno, *rh_endcap_deno_ZZ_70to100,  0);
    RooHistPdf *pdf_barrel_nume_ZZ_100to500 = new RooHistPdf("pdf_barrel_nume_ZZ_100to500", "Numerator barrel MC ZZ template 200 to 500",   iso_nume, *rh_barrel_nume_ZZ_100to500, 0);
    RooHistPdf *pdf_endcap_nume_ZZ_100to500 = new RooHistPdf("pdf_endcap_nume_ZZ_100to500", "Numerator endcap MC ZZ template 200 to 500",   iso_nume, *rh_endcap_nume_ZZ_100to500, 0);
    RooHistPdf *pdf_barrel_deno_ZZ_100to500 = new RooHistPdf("pdf_barrel_deno_ZZ_100to500", "Denominator barrel MC ZZ template 200 to 500", iso_deno, *rh_barrel_deno_ZZ_100to500, 0);
    RooHistPdf *pdf_endcap_deno_ZZ_100to500 = new RooHistPdf("pdf_endcap_deno_ZZ_100to500", "Denominator endcap MC ZZ template 200 to 500", iso_deno, *rh_endcap_deno_ZZ_100to500, 0);

    // Constraints for integrals
    Double_t N_barrel_nume_ttbar_50to70   = h_barrel_MC_nume_50to70  [_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar_50to70   = h_endcap_MC_nume_50to70  [_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar_50to70   = h_barrel_MC_deno_50to70  [_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar_50to70   = h_endcap_MC_deno_50to70  [_ttbar]->Integral();
    Double_t N_barrel_nume_ttbar_70to100  = h_barrel_MC_nume_70to100 [_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar_70to100  = h_endcap_MC_nume_70to100 [_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar_70to100  = h_barrel_MC_deno_70to100 [_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar_70to100  = h_endcap_MC_deno_70to100 [_ttbar]->Integral();
    Double_t N_barrel_nume_ttbar_100to500 = h_barrel_MC_nume_100to500[_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar_100to500 = h_endcap_MC_nume_100to500[_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar_100to500 = h_barrel_MC_deno_100to500[_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar_100to500 = h_endcap_MC_deno_100to500[_ttbar]->Integral();

    Double_t N_barrel_nume_WJets_50to70   = h_barrel_MC_nume_50to70  [_WJets]->Integral();
    Double_t N_endcap_nume_WJets_50to70   = h_endcap_MC_nume_50to70  [_WJets]->Integral();
    Double_t N_barrel_deno_WJets_50to70   = h_barrel_MC_deno_50to70  [_WJets]->Integral();
    Double_t N_endcap_deno_WJets_50to70   = h_endcap_MC_deno_50to70  [_WJets]->Integral();
    Double_t N_barrel_nume_WJets_70to100  = h_barrel_MC_nume_70to100 [_WJets]->Integral();
    Double_t N_endcap_nume_WJets_70to100  = h_endcap_MC_nume_70to100 [_WJets]->Integral();
    Double_t N_barrel_deno_WJets_70to100  = h_barrel_MC_deno_70to100 [_WJets]->Integral();
    Double_t N_endcap_deno_WJets_70to100  = h_endcap_MC_deno_70to100 [_WJets]->Integral();
    Double_t N_barrel_nume_WJets_100to500 = h_barrel_MC_nume_100to500[_WJets]->Integral();
    Double_t N_endcap_nume_WJets_100to500 = h_endcap_MC_nume_100to500[_WJets]->Integral();
    Double_t N_barrel_deno_WJets_100to500 = h_barrel_MC_deno_100to500[_WJets]->Integral();
    Double_t N_endcap_deno_WJets_100to500 = h_endcap_MC_deno_100to500[_WJets]->Integral();

    Double_t N_barrel_nume_DY_50to70   = h_barrel_MC_nume_50to70  [_DY_Full]->Integral();
    Double_t N_endcap_nume_DY_50to70   = h_endcap_MC_nume_50to70  [_DY_Full]->Integral();
    Double_t N_barrel_deno_DY_50to70   = h_barrel_MC_deno_50to70  [_DY_Full]->Integral();
    Double_t N_endcap_deno_DY_50to70   = h_endcap_MC_deno_50to70  [_DY_Full]->Integral();
    Double_t N_barrel_nume_DY_70to100  = h_barrel_MC_nume_70to100 [_DY_Full]->Integral();
    Double_t N_endcap_nume_DY_70to100  = h_endcap_MC_nume_70to100 [_DY_Full]->Integral();
    Double_t N_barrel_deno_DY_70to100  = h_barrel_MC_deno_70to100 [_DY_Full]->Integral();
    Double_t N_endcap_deno_DY_70to100  = h_endcap_MC_deno_70to100 [_DY_Full]->Integral();
    Double_t N_barrel_nume_DY_100to500 = h_barrel_MC_nume_100to500[_DY_Full]->Integral();
    Double_t N_endcap_nume_DY_100to500 = h_endcap_MC_nume_100to500[_DY_Full]->Integral();
    Double_t N_barrel_deno_DY_100to500 = h_barrel_MC_deno_100to500[_DY_Full]->Integral();
    Double_t N_endcap_deno_DY_100to500 = h_endcap_MC_deno_100to500[_DY_Full]->Integral();

    Double_t N_barrel_nume_QCD_50to70   = h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD_50to70   = h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD_50to70   = h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD_50to70   = h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_nume_QCD_70to100  = h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD_70to100  = h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD_70to100  = h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD_70to100  = h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_nume_QCD_100to500 = h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD_100to500 = h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD_100to500 = h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD_100to500 = h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->Integral();

    Double_t N_barrel_nume_tW_50to70   = h_barrel_MC_nume_50to70  [_tW]->Integral();
    Double_t N_endcap_nume_tW_50to70   = h_endcap_MC_nume_50to70  [_tW]->Integral();
    Double_t N_barrel_deno_tW_50to70   = h_barrel_MC_deno_50to70  [_tW]->Integral();
    Double_t N_endcap_deno_tW_50to70   = h_endcap_MC_deno_50to70  [_tW]->Integral();
    Double_t N_barrel_nume_tW_70to100  = h_barrel_MC_nume_70to100 [_tW]->Integral();
    Double_t N_endcap_nume_tW_70to100  = h_endcap_MC_nume_70to100 [_tW]->Integral();
    Double_t N_barrel_deno_tW_70to100  = h_barrel_MC_deno_70to100 [_tW]->Integral();
    Double_t N_endcap_deno_tW_70to100  = h_endcap_MC_deno_70to100 [_tW]->Integral();
    Double_t N_barrel_nume_tW_100to500 = h_barrel_MC_nume_100to500[_tW]->Integral();
    Double_t N_endcap_nume_tW_100to500 = h_endcap_MC_nume_100to500[_tW]->Integral();
    Double_t N_barrel_deno_tW_100to500 = h_barrel_MC_deno_100to500[_tW]->Integral();
    Double_t N_endcap_deno_tW_100to500 = h_endcap_MC_deno_100to500[_tW]->Integral();

    Double_t N_barrel_nume_tbarW_50to70   = h_barrel_MC_nume_50to70  [_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW_50to70   = h_endcap_MC_nume_50to70  [_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW_50to70   = h_barrel_MC_deno_50to70  [_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW_50to70   = h_endcap_MC_deno_50to70  [_tbarW]->Integral();
    Double_t N_barrel_nume_tbarW_70to100  = h_barrel_MC_nume_70to100 [_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW_70to100  = h_endcap_MC_nume_70to100 [_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW_70to100  = h_barrel_MC_deno_70to100 [_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW_70to100  = h_endcap_MC_deno_70to100 [_tbarW]->Integral();
    Double_t N_barrel_nume_tbarW_100to500 = h_barrel_MC_nume_100to500[_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW_100to500 = h_endcap_MC_nume_100to500[_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW_100to500 = h_barrel_MC_deno_100to500[_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW_100to500 = h_endcap_MC_deno_100to500[_tbarW]->Integral();

    Double_t N_barrel_nume_WW_50to70   = h_barrel_MC_nume_50to70  [_WW]->Integral();
    Double_t N_endcap_nume_WW_50to70   = h_endcap_MC_nume_50to70  [_WW]->Integral();
    Double_t N_barrel_deno_WW_50to70   = h_barrel_MC_deno_50to70  [_WW]->Integral();
    Double_t N_endcap_deno_WW_50to70   = h_endcap_MC_deno_50to70  [_WW]->Integral();
    Double_t N_barrel_nume_WW_70to100  = h_barrel_MC_nume_70to100 [_WW]->Integral();
    Double_t N_endcap_nume_WW_70to100  = h_endcap_MC_nume_70to100 [_WW]->Integral();
    Double_t N_barrel_deno_WW_70to100  = h_barrel_MC_deno_70to100 [_WW]->Integral();
    Double_t N_endcap_deno_WW_70to100  = h_endcap_MC_deno_70to100 [_WW]->Integral();
    Double_t N_barrel_nume_WW_100to500 = h_barrel_MC_nume_100to500[_WW]->Integral();
    Double_t N_endcap_nume_WW_100to500 = h_endcap_MC_nume_100to500[_WW]->Integral();
    Double_t N_barrel_deno_WW_100to500 = h_barrel_MC_deno_100to500[_WW]->Integral();
    Double_t N_endcap_deno_WW_100to500 = h_endcap_MC_deno_100to500[_WW]->Integral();

    Double_t N_barrel_nume_WZ_50to70   = h_barrel_MC_nume_50to70  [_WZ]->Integral();
    Double_t N_endcap_nume_WZ_50to70   = h_endcap_MC_nume_50to70  [_WZ]->Integral();
    Double_t N_barrel_deno_WZ_50to70   = h_barrel_MC_deno_50to70  [_WZ]->Integral();
    Double_t N_endcap_deno_WZ_50to70   = h_endcap_MC_deno_50to70  [_WZ]->Integral();
    Double_t N_barrel_nume_WZ_70to100  = h_barrel_MC_nume_70to100 [_WZ]->Integral();
    Double_t N_endcap_nume_WZ_70to100  = h_endcap_MC_nume_70to100 [_WZ]->Integral();
    Double_t N_barrel_deno_WZ_70to100  = h_barrel_MC_deno_70to100 [_WZ]->Integral();
    Double_t N_endcap_deno_WZ_70to100  = h_endcap_MC_deno_70to100 [_WZ]->Integral();
    Double_t N_barrel_nume_WZ_100to500 = h_barrel_MC_nume_100to500[_WZ]->Integral();
    Double_t N_endcap_nume_WZ_100to500 = h_endcap_MC_nume_100to500[_WZ]->Integral();
    Double_t N_barrel_deno_WZ_100to500 = h_barrel_MC_deno_100to500[_WZ]->Integral();
    Double_t N_endcap_deno_WZ_100to500 = h_endcap_MC_deno_100to500[_WZ]->Integral();

    Double_t N_barrel_nume_ZZ_50to70   = h_barrel_MC_nume_50to70  [_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ_50to70   = h_endcap_MC_nume_50to70  [_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ_50to70   = h_barrel_MC_deno_50to70  [_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ_50to70   = h_endcap_MC_deno_50to70  [_ZZ]->Integral();
    Double_t N_barrel_nume_ZZ_70to100  = h_barrel_MC_nume_70to100 [_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ_70to100  = h_endcap_MC_nume_70to100 [_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ_70to100  = h_barrel_MC_deno_70to100 [_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ_70to100  = h_endcap_MC_deno_70to100 [_ZZ]->Integral();
    Double_t N_barrel_nume_ZZ_100to500 = h_barrel_MC_nume_100to500[_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ_100to500 = h_endcap_MC_nume_100to500[_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ_100to500 = h_barrel_MC_deno_100to500[_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ_100to500 = h_endcap_MC_deno_100to500[_ZZ]->Integral();

    Double_t N_barrel_nume_total_50to70 = N_barrel_nume_ttbar_50to70 + N_barrel_nume_WJets_50to70 + N_barrel_nume_DY_50to70 + N_barrel_nume_QCD_50to70 +
                                          N_barrel_nume_tW_50to70    + N_barrel_nume_tbarW_50to70 + N_barrel_nume_WW_50to70 + N_barrel_nume_WZ_50to70  +
                                          N_barrel_nume_ZZ_50to70;
    Double_t N_endcap_nume_total_50to70 = N_endcap_nume_ttbar_50to70 + N_endcap_nume_WJets_50to70 + N_endcap_nume_DY_50to70 + N_endcap_nume_QCD_50to70 +
                                          N_endcap_nume_tW_50to70    + N_endcap_nume_tbarW_50to70 + N_endcap_nume_WW_50to70 + N_endcap_nume_WZ_50to70  +
                                          N_endcap_nume_ZZ_50to70;
    Double_t N_barrel_deno_total_50to70 = N_barrel_deno_ttbar_50to70 + N_barrel_deno_WJets_50to70 + N_barrel_deno_DY_50to70 + N_barrel_deno_QCD_50to70 +
                                          N_barrel_deno_tW_50to70    + N_barrel_deno_tbarW_50to70 + N_barrel_deno_WW_50to70 + N_barrel_deno_WZ_50to70  +
                                          N_barrel_deno_ZZ_50to70;
    Double_t N_endcap_deno_total_50to70 = N_endcap_deno_ttbar_50to70 + N_endcap_deno_WJets_50to70 + N_endcap_deno_DY_50to70 + N_endcap_deno_QCD_50to70 +
                                          N_endcap_deno_tW_50to70    + N_endcap_deno_tbarW_50to70 + N_endcap_deno_WW_50to70 + N_endcap_deno_WZ_50to70  +
                                          N_endcap_deno_ZZ_50to70;
    Double_t N_barrel_nume_total_70to100 = N_barrel_nume_ttbar_70to100 + N_barrel_nume_WJets_70to100 + N_barrel_nume_DY_70to100 + N_barrel_nume_QCD_70to100 +
                                           N_barrel_nume_tW_70to100    + N_barrel_nume_tbarW_70to100 + N_barrel_nume_WW_70to100 + N_barrel_nume_WZ_70to100  +
                                           N_barrel_nume_ZZ_70to100;
    Double_t N_endcap_nume_total_70to100 = N_endcap_nume_ttbar_70to100 + N_endcap_nume_WJets_70to100 + N_endcap_nume_DY_70to100 + N_endcap_nume_QCD_70to100 +
                                           N_endcap_nume_tW_70to100    + N_endcap_nume_tbarW_70to100 + N_endcap_nume_WW_70to100 + N_endcap_nume_WZ_70to100  +
                                           N_endcap_nume_ZZ_70to100;
    Double_t N_barrel_deno_total_70to100 = N_barrel_deno_ttbar_70to100 + N_barrel_deno_WJets_70to100 + N_barrel_deno_DY_70to100 + N_barrel_deno_QCD_70to100 +
                                           N_barrel_deno_tW_70to100    + N_barrel_deno_tbarW_70to100 + N_barrel_deno_WW_70to100 + N_barrel_deno_WZ_70to100  +
                                           N_barrel_deno_ZZ_70to100;
    Double_t N_endcap_deno_total_70to100 = N_endcap_deno_ttbar_70to100 + N_endcap_deno_WJets_70to100 + N_endcap_deno_DY_70to100 + N_endcap_deno_QCD_70to100 +
                                           N_endcap_deno_tW_70to100    + N_endcap_deno_tbarW_70to100 + N_endcap_deno_WW_70to100 + N_endcap_deno_WZ_70to100  +
                                           N_endcap_deno_ZZ_70to100;
    Double_t N_barrel_nume_total_100to500 = N_barrel_nume_ttbar_100to500 + N_barrel_nume_WJets_100to500 + N_barrel_nume_DY_100to500 + N_barrel_nume_QCD_100to500 +
                                            N_barrel_nume_tW_100to500    + N_barrel_nume_tbarW_100to500 + N_barrel_nume_WW_100to500 + N_barrel_nume_WZ_100to500  +
                                            N_barrel_nume_ZZ_100to500;
    Double_t N_endcap_nume_total_100to500 = N_endcap_nume_ttbar_100to500 + N_endcap_nume_WJets_100to500 + N_endcap_nume_DY_100to500 + N_endcap_nume_QCD_100to500 +
                                            N_endcap_nume_tW_100to500    + N_endcap_nume_tbarW_100to500 + N_endcap_nume_WW_100to500 + N_endcap_nume_WZ_100to500  +
                                            N_endcap_nume_ZZ_100to500;
    Double_t N_barrel_deno_total_100to500 = N_barrel_deno_ttbar_100to500 + N_barrel_deno_WJets_100to500 + N_barrel_deno_DY_100to500 + N_barrel_deno_QCD_100to500 +
                                            N_barrel_deno_tW_100to500    + N_barrel_deno_tbarW_100to500 + N_barrel_deno_WW_100to500 + N_barrel_deno_WZ_100to500  +
                                            N_barrel_deno_ZZ_100to500;
    Double_t N_endcap_deno_total_100to500 = N_endcap_deno_ttbar_100to500 + N_endcap_deno_WJets_100to500 + N_endcap_deno_DY_100to500 + N_endcap_deno_QCD_100to500 +
                                            N_endcap_deno_tW_100to500    + N_endcap_deno_tbarW_100to500 + N_endcap_deno_WW_100to500 + N_endcap_deno_WZ_100to500  +
                                            N_endcap_deno_ZZ_100to500;

    Double_t Nnorm_barrel_nume_ttbar_50to70   = N_barrel_nume_ttbar_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_ttbar_50to70   = N_endcap_nume_ttbar_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_ttbar_50to70   = N_barrel_deno_ttbar_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_ttbar_50to70   = N_endcap_deno_ttbar_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_ttbar_70to100  = N_barrel_nume_ttbar_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_ttbar_70to100  = N_endcap_nume_ttbar_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_ttbar_70to100  = N_barrel_deno_ttbar_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_ttbar_70to100  = N_endcap_deno_ttbar_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_ttbar_100to500 = N_barrel_nume_ttbar_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_ttbar_100to500 = N_endcap_nume_ttbar_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_ttbar_100to500 = N_barrel_deno_ttbar_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_ttbar_100to500 = N_endcap_deno_ttbar_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_WJets_50to70   = N_barrel_nume_WJets_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WJets_50to70   = N_endcap_nume_WJets_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WJets_50to70   = N_barrel_deno_WJets_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WJets_50to70   = N_endcap_deno_WJets_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WJets_70to100  = N_barrel_nume_WJets_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WJets_70to100  = N_endcap_nume_WJets_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WJets_70to100  = N_barrel_deno_WJets_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WJets_70to100  = N_endcap_deno_WJets_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WJets_100to500 = N_barrel_nume_WJets_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_WJets_100to500 = N_endcap_nume_WJets_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_WJets_100to500 = N_barrel_deno_WJets_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_WJets_100to500 = N_endcap_deno_WJets_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_DY_50to70   = N_barrel_nume_DY_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_DY_50to70   = N_endcap_nume_DY_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_DY_50to70   = N_barrel_deno_DY_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_DY_50to70   = N_endcap_deno_DY_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_DY_70to100  = N_barrel_nume_DY_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_DY_70to100  = N_endcap_nume_DY_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_DY_70to100  = N_barrel_deno_DY_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_DY_70to100  = N_endcap_deno_DY_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_DY_100to500 = N_barrel_nume_DY_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_DY_100to500 = N_endcap_nume_DY_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_DY_100to500 = N_barrel_deno_DY_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_DY_100to500 = N_endcap_deno_DY_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_tW_50to70   = N_barrel_nume_tW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_tW_50to70   = N_endcap_nume_tW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_tW_50to70   = N_barrel_deno_tW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_tW_50to70   = N_endcap_deno_tW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_tW_70to100  = N_barrel_nume_tW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_tW_70to100  = N_endcap_nume_tW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_tW_70to100  = N_barrel_deno_tW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_tW_70to100  = N_endcap_deno_tW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_tW_100to500 = N_barrel_nume_tW_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_tW_100to500 = N_endcap_nume_tW_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_tW_100to500 = N_barrel_deno_tW_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_tW_100to500 = N_endcap_deno_tW_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_tbarW_50to70   = N_barrel_nume_tbarW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_tbarW_50to70   = N_endcap_nume_tbarW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_tbarW_50to70   = N_barrel_deno_tbarW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_tbarW_50to70   = N_endcap_deno_tbarW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_tbarW_70to100  = N_barrel_nume_tbarW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_tbarW_70to100  = N_endcap_nume_tbarW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_tbarW_70to100  = N_barrel_deno_tbarW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_tbarW_70to100  = N_endcap_deno_tbarW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_tbarW_100to500 = N_barrel_nume_tbarW_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_tbarW_100to500 = N_endcap_nume_tbarW_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_tbarW_100to500 = N_barrel_deno_tbarW_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_tbarW_100to500 = N_endcap_deno_tbarW_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_WW_50to70   = N_barrel_nume_WW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WW_50to70   = N_endcap_nume_WW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WW_50to70   = N_barrel_deno_WW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WW_50to70   = N_endcap_deno_WW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WW_70to100  = N_barrel_nume_WW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WW_70to100  = N_endcap_nume_WW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WW_70to100  = N_barrel_deno_WW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WW_70to100  = N_endcap_deno_WW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WW_100to500 = N_barrel_nume_WW_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_WW_100to500 = N_endcap_nume_WW_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_WW_100to500 = N_barrel_deno_WW_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_WW_100to500 = N_endcap_deno_WW_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_WZ_50to70   = N_barrel_nume_WZ_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WZ_50to70   = N_endcap_nume_WZ_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WZ_50to70   = N_barrel_deno_WZ_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WZ_50to70   = N_endcap_deno_WZ_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WZ_70to100  = N_barrel_nume_WZ_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WZ_70to100  = N_endcap_nume_WZ_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WZ_70to100  = N_barrel_deno_WZ_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WZ_70to100  = N_endcap_deno_WZ_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WZ_100to500 = N_barrel_nume_WZ_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_WZ_100to500 = N_endcap_nume_WZ_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_WZ_100to500 = N_barrel_deno_WZ_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_WZ_100to500 = N_endcap_deno_WZ_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_ZZ_50to70   = N_barrel_nume_ZZ_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_ZZ_50to70   = N_endcap_nume_ZZ_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_ZZ_50to70   = N_barrel_deno_ZZ_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_ZZ_50to70   = N_endcap_deno_ZZ_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_ZZ_70to100  = N_barrel_nume_ZZ_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_ZZ_70to100  = N_endcap_nume_ZZ_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_ZZ_70to100  = N_barrel_deno_ZZ_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_ZZ_70to100  = N_endcap_deno_ZZ_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_ZZ_100to500 = N_barrel_nume_ZZ_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_ZZ_100to500 = N_endcap_nume_ZZ_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_ZZ_100to500 = N_barrel_deno_ZZ_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_ZZ_100to500 = N_endcap_deno_ZZ_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    Double_t Nnorm_barrel_nume_QCD_50to70   = N_barrel_nume_QCD_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_QCD_50to70   = N_endcap_nume_QCD_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_QCD_50to70   = N_barrel_deno_QCD_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_QCD_50to70   = N_endcap_deno_QCD_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_QCD_70to100  = N_barrel_nume_QCD_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_QCD_70to100  = N_endcap_nume_QCD_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_QCD_70to100  = N_barrel_deno_QCD_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_QCD_70to100  = N_endcap_deno_QCD_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_QCD_100to500 = N_barrel_nume_QCD_100to500 * h_barrel_data_nume_100to500->Integral() / N_barrel_nume_total_100to500;
    Double_t Nnorm_endcap_nume_QCD_100to500 = N_endcap_nume_QCD_100to500 * h_endcap_data_nume_100to500->Integral() / N_endcap_nume_total_100to500;
    Double_t Nnorm_barrel_deno_QCD_100to500 = N_barrel_deno_QCD_100to500 * h_barrel_data_deno_100to500->Integral() / N_barrel_deno_total_100to500;
    Double_t Nnorm_endcap_deno_QCD_100to500 = N_endcap_deno_QCD_100to500 * h_endcap_data_deno_100to500->Integral() / N_endcap_deno_total_100to500;

    // Fit constraints
    RooRealVar n_barrel_nume_ttbar_50to70  ("n_barrel_nume_ttbar_50to70",   "n_barrel_nume_ttbar_50to70",   Nnorm_barrel_nume_ttbar_50to70 ,  Nnorm_barrel_nume_ttbar_50to70  *0.75, Nnorm_barrel_nume_ttbar_50to70  *1.25);
    RooRealVar n_endcap_nume_ttbar_50to70  ("n_endcap_nume_ttbar_50to70",   "n_endcap_nume_ttbar_50to70",   Nnorm_endcap_nume_ttbar_50to70 ,  Nnorm_endcap_nume_ttbar_50to70  *0.75, Nnorm_endcap_nume_ttbar_50to70  *1.25);
    RooRealVar n_barrel_deno_ttbar_50to70  ("n_barrel_deno_ttbar_50to70",   "n_barrel_deno_ttbar_50to70",   Nnorm_barrel_deno_ttbar_50to70 ,  Nnorm_barrel_deno_ttbar_50to70  *0.75, Nnorm_barrel_deno_ttbar_50to70  *1.25);
    RooRealVar n_endcap_deno_ttbar_50to70  ("n_endcap_deno_ttbar_50to70",   "n_endcap_deno_ttbar_50to70",   Nnorm_endcap_deno_ttbar_50to70 ,  Nnorm_endcap_deno_ttbar_50to70  *0.75, Nnorm_endcap_deno_ttbar_50to70  *1.25);
    RooRealVar n_barrel_nume_ttbar_70to100 ("n_barrel_nume_ttbar_70to100",  "n_barrel_nume_ttbar_70to100",  Nnorm_barrel_nume_ttbar_70to100,  Nnorm_barrel_nume_ttbar_70to100 *0.75, Nnorm_barrel_nume_ttbar_70to100 *1.25);
    RooRealVar n_endcap_nume_ttbar_70to100 ("n_endcap_nume_ttbar_70to100",  "n_endcap_nume_ttbar_70to100",  Nnorm_endcap_nume_ttbar_70to100,  Nnorm_endcap_nume_ttbar_70to100 *0.75, Nnorm_endcap_nume_ttbar_70to100 *1.25);
    RooRealVar n_barrel_deno_ttbar_70to100 ("n_barrel_deno_ttbar_70to100",  "n_barrel_deno_ttbar_70to100",  Nnorm_barrel_deno_ttbar_70to100,  Nnorm_barrel_deno_ttbar_70to100 *0.75, Nnorm_barrel_deno_ttbar_70to100 *1.25);
    RooRealVar n_endcap_deno_ttbar_70to100 ("n_endcap_deno_ttbar_70to100",  "n_endcap_deno_ttbar_70to100",  Nnorm_endcap_deno_ttbar_70to100,  Nnorm_endcap_deno_ttbar_70to100 *0.75, Nnorm_endcap_deno_ttbar_70to100 *1.25);
    RooRealVar n_barrel_nume_ttbar_100to500("n_barrel_nume_ttbar_100to500", "n_barrel_nume_ttbar_100to500", Nnorm_barrel_nume_ttbar_100to500, Nnorm_barrel_nume_ttbar_100to500*0.75, Nnorm_barrel_nume_ttbar_100to500*1.25);
    RooRealVar n_endcap_nume_ttbar_100to500("n_endcap_nume_ttbar_100to500", "n_endcap_nume_ttbar_100to500", Nnorm_endcap_nume_ttbar_100to500, Nnorm_endcap_nume_ttbar_100to500*0.75, Nnorm_endcap_nume_ttbar_100to500*1.25);
    RooRealVar n_barrel_deno_ttbar_100to500("n_barrel_deno_ttbar_100to500", "n_barrel_deno_ttbar_100to500", Nnorm_barrel_deno_ttbar_100to500, Nnorm_barrel_deno_ttbar_100to500*0.75, Nnorm_barrel_deno_ttbar_100to500*1.25);
    RooRealVar n_endcap_deno_ttbar_100to500("n_endcap_deno_ttbar_100to500", "n_endcap_deno_ttbar_100to500", Nnorm_endcap_deno_ttbar_100to500, Nnorm_endcap_deno_ttbar_100to500*0.75, Nnorm_endcap_deno_ttbar_100to500*1.25);

    RooRealVar n_barrel_nume_WJets_50to70  ("n_barrel_nume_WJets_50to70",   "n_barrel_nume_WJets_50to70",   Nnorm_barrel_nume_WJets_50to70 ,  Nnorm_barrel_nume_WJets_50to70  *0.5, Nnorm_barrel_nume_WJets_50to70  *1.5);
    RooRealVar n_endcap_nume_WJets_50to70  ("n_endcap_nume_WJets_50to70",   "n_endcap_nume_WJets_50to70",   Nnorm_endcap_nume_WJets_50to70 ,  Nnorm_endcap_nume_WJets_50to70  *0.5, Nnorm_endcap_nume_WJets_50to70  *1.5);
    RooRealVar n_barrel_deno_WJets_50to70  ("n_barrel_deno_WJets_50to70",   "n_barrel_deno_WJets_50to70",   Nnorm_barrel_deno_WJets_50to70 ,  Nnorm_barrel_deno_WJets_50to70  *0.5, Nnorm_barrel_deno_WJets_50to70  *1.5);
    RooRealVar n_endcap_deno_WJets_50to70  ("n_endcap_deno_WJets_50to70",   "n_endcap_deno_WJets_50to70",   Nnorm_endcap_deno_WJets_50to70 ,  Nnorm_endcap_deno_WJets_50to70  *0.5, Nnorm_endcap_deno_WJets_50to70  *1.5);
    RooRealVar n_barrel_nume_WJets_70to100 ("n_barrel_nume_WJets_70to100",  "n_barrel_nume_WJets_70to100",  Nnorm_barrel_nume_WJets_70to100,  Nnorm_barrel_nume_WJets_70to100 *0.5, Nnorm_barrel_nume_WJets_70to100 *1.5);
    RooRealVar n_endcap_nume_WJets_70to100 ("n_endcap_nume_WJets_70to100",  "n_endcap_nume_WJets_70to100",  Nnorm_endcap_nume_WJets_70to100,  Nnorm_endcap_nume_WJets_70to100 *0.5, Nnorm_endcap_nume_WJets_70to100 *1.5);
    RooRealVar n_barrel_deno_WJets_70to100 ("n_barrel_deno_WJets_70to100",  "n_barrel_deno_WJets_70to100",  Nnorm_barrel_deno_WJets_70to100,  Nnorm_barrel_deno_WJets_70to100 *0.5, Nnorm_barrel_deno_WJets_70to100 *1.5);
    RooRealVar n_endcap_deno_WJets_70to100 ("n_endcap_deno_WJets_70to100",  "n_endcap_deno_WJets_70to100",  Nnorm_endcap_deno_WJets_70to100,  Nnorm_endcap_deno_WJets_70to100 *0.5, Nnorm_endcap_deno_WJets_70to100 *1.5);
    RooRealVar n_barrel_nume_WJets_100to500("n_barrel_nume_WJets_100to500", "n_barrel_nume_WJets_100to500", Nnorm_barrel_nume_WJets_100to500, Nnorm_barrel_nume_WJets_100to500*0.5, Nnorm_barrel_nume_WJets_100to500*1.5);
    RooRealVar n_endcap_nume_WJets_100to500("n_endcap_nume_WJets_100to500", "n_endcap_nume_WJets_100to500", Nnorm_endcap_nume_WJets_100to500, Nnorm_endcap_nume_WJets_100to500*0.5, Nnorm_endcap_nume_WJets_100to500*1.5);
    RooRealVar n_barrel_deno_WJets_100to500("n_barrel_deno_WJets_100to500", "n_barrel_deno_WJets_100to500", Nnorm_barrel_deno_WJets_100to500, Nnorm_barrel_deno_WJets_100to500*0.5, Nnorm_barrel_deno_WJets_100to500*1.5);
    RooRealVar n_endcap_deno_WJets_100to500("n_endcap_deno_WJets_100to500", "n_endcap_deno_WJets_100to500", Nnorm_endcap_deno_WJets_100to500, Nnorm_endcap_deno_WJets_100to500*0.5, Nnorm_endcap_deno_WJets_100to500*1.5);

    RooRealVar n_barrel_nume_DY_50to70  ("n_barrel_nume_DY_50to70",   "n_barrel_nume_DY_50to70",   Nnorm_barrel_nume_DY_50to70 ,  Nnorm_barrel_nume_DY_50to70  *0.8, Nnorm_barrel_nume_DY_50to70  *1.2);
    RooRealVar n_endcap_nume_DY_50to70  ("n_endcap_nume_DY_50to70",   "n_endcap_nume_DY_50to70",   Nnorm_endcap_nume_DY_50to70 ,  Nnorm_endcap_nume_DY_50to70  *0.8, Nnorm_endcap_nume_DY_50to70  *1.2);
    RooRealVar n_barrel_deno_DY_50to70  ("n_barrel_deno_DY_50to70",   "n_barrel_deno_DY_50to70",   Nnorm_barrel_deno_DY_50to70 ,  Nnorm_barrel_deno_DY_50to70  *0.8, Nnorm_barrel_deno_DY_50to70  *1.2);
    RooRealVar n_endcap_deno_DY_50to70  ("n_endcap_deno_DY_50to70",   "n_endcap_deno_DY_50to70",   Nnorm_endcap_deno_DY_50to70 ,  Nnorm_endcap_deno_DY_50to70  *0.8, Nnorm_endcap_deno_DY_50to70  *1.2);
    RooRealVar n_barrel_nume_DY_70to100 ("n_barrel_nume_DY_70to100",  "n_barrel_nume_DY_70to100",  Nnorm_barrel_nume_DY_70to100,  Nnorm_barrel_nume_DY_70to100 *0.8, Nnorm_barrel_nume_DY_70to100 *1.2);
    RooRealVar n_endcap_nume_DY_70to100 ("n_endcap_nume_DY_70to100",  "n_endcap_nume_DY_70to100",  Nnorm_endcap_nume_DY_70to100,  Nnorm_endcap_nume_DY_70to100 *0.8, Nnorm_endcap_nume_DY_70to100 *1.2);
    RooRealVar n_barrel_deno_DY_70to100 ("n_barrel_deno_DY_70to100",  "n_barrel_deno_DY_70to100",  Nnorm_barrel_deno_DY_70to100,  Nnorm_barrel_deno_DY_70to100 *0.8, Nnorm_barrel_deno_DY_70to100 *1.2);
    RooRealVar n_endcap_deno_DY_70to100 ("n_endcap_deno_DY_70to100",  "n_endcap_deno_DY_70to100",  Nnorm_endcap_deno_DY_70to100,  Nnorm_endcap_deno_DY_70to100 *0.8, Nnorm_endcap_deno_DY_70to100 *1.2);
    RooRealVar n_barrel_nume_DY_100to500("n_barrel_nume_DY_100to500", "n_barrel_nume_DY_100to500", Nnorm_barrel_nume_DY_100to500, Nnorm_barrel_nume_DY_100to500*0.8, Nnorm_barrel_nume_DY_100to500*1.2);
    RooRealVar n_endcap_nume_DY_100to500("n_endcap_nume_DY_100to500", "n_endcap_nume_DY_100to500", Nnorm_endcap_nume_DY_100to500, Nnorm_endcap_nume_DY_100to500*0.8, Nnorm_endcap_nume_DY_100to500*1.2);
    RooRealVar n_barrel_deno_DY_100to500("n_barrel_deno_DY_100to500", "n_barrel_deno_DY_100to500", Nnorm_barrel_deno_DY_100to500, Nnorm_barrel_deno_DY_100to500*0.8, Nnorm_barrel_deno_DY_100to500*1.2);
    RooRealVar n_endcap_deno_DY_100to500("n_endcap_deno_DY_100to500", "n_endcap_deno_DY_100to500", Nnorm_endcap_deno_DY_100to500, Nnorm_endcap_deno_DY_100to500*0.8, Nnorm_endcap_deno_DY_100to500*1.2);

    RooRealVar n_barrel_nume_tW_50to70  ("n_barrel_nume_tW_50to70",   "n_barrel_nume_tW_50to70",   Nnorm_barrel_nume_tW_50to70 ,  Nnorm_barrel_nume_tW_50to70  *0.8, Nnorm_barrel_nume_tW_50to70  *1.2);
    RooRealVar n_endcap_nume_tW_50to70  ("n_endcap_nume_tW_50to70",   "n_endcap_nume_tW_50to70",   Nnorm_endcap_nume_tW_50to70 ,  Nnorm_endcap_nume_tW_50to70  *0.8, Nnorm_endcap_nume_tW_50to70  *1.2);
    RooRealVar n_barrel_deno_tW_50to70  ("n_barrel_deno_tW_50to70",   "n_barrel_deno_tW_50to70",   Nnorm_barrel_deno_tW_50to70 ,  Nnorm_barrel_deno_tW_50to70  *0.8, Nnorm_barrel_deno_tW_50to70  *1.2);
    RooRealVar n_endcap_deno_tW_50to70  ("n_endcap_deno_tW_50to70",   "n_endcap_deno_tW_50to70",   Nnorm_endcap_deno_tW_50to70 ,  Nnorm_endcap_deno_tW_50to70  *0.8, Nnorm_endcap_deno_tW_50to70  *1.2);
    RooRealVar n_barrel_nume_tW_70to100 ("n_barrel_nume_tW_70to100",  "n_barrel_nume_tW_70to100",  Nnorm_barrel_nume_tW_70to100,  Nnorm_barrel_nume_tW_70to100 *0.8, Nnorm_barrel_nume_tW_70to100 *1.2);
    RooRealVar n_endcap_nume_tW_70to100 ("n_endcap_nume_tW_70to100",  "n_endcap_nume_tW_70to100",  Nnorm_endcap_nume_tW_70to100,  Nnorm_endcap_nume_tW_70to100 *0.8, Nnorm_endcap_nume_tW_70to100 *1.2);
    RooRealVar n_barrel_deno_tW_70to100 ("n_barrel_deno_tW_70to100",  "n_barrel_deno_tW_70to100",  Nnorm_barrel_deno_tW_70to100,  Nnorm_barrel_deno_tW_70to100 *0.8, Nnorm_barrel_deno_tW_70to100 *1.2);
    RooRealVar n_endcap_deno_tW_70to100 ("n_endcap_deno_tW_70to100",  "n_endcap_deno_tW_70to100",  Nnorm_endcap_deno_tW_70to100,  Nnorm_endcap_deno_tW_70to100 *0.8, Nnorm_endcap_deno_tW_70to100 *1.2);
    RooRealVar n_barrel_nume_tW_100to500("n_barrel_nume_tW_100to500", "n_barrel_nume_tW_100to500", Nnorm_barrel_nume_tW_100to500, Nnorm_barrel_nume_tW_100to500*0.8, Nnorm_barrel_nume_tW_100to500*1.2);
    RooRealVar n_endcap_nume_tW_100to500("n_endcap_nume_tW_100to500", "n_endcap_nume_tW_100to500", Nnorm_endcap_nume_tW_100to500, Nnorm_endcap_nume_tW_100to500*0.8, Nnorm_endcap_nume_tW_100to500*1.2);
    RooRealVar n_barrel_deno_tW_100to500("n_barrel_deno_tW_100to500", "n_barrel_deno_tW_100to500", Nnorm_barrel_deno_tW_100to500, Nnorm_barrel_deno_tW_100to500*0.8, Nnorm_barrel_deno_tW_100to500*1.2);
    RooRealVar n_endcap_deno_tW_100to500("n_endcap_deno_tW_100to500", "n_endcap_deno_tW_100to500", Nnorm_endcap_deno_tW_100to500, Nnorm_endcap_deno_tW_100to500*0.8, Nnorm_endcap_deno_tW_100to500*1.2);

    RooRealVar n_barrel_nume_tbarW_50to70  ("n_barrel_nume_tbarW_50to70",   "n_barrel_nume_tbarW_50to70",   Nnorm_barrel_nume_tbarW_50to70 ,  Nnorm_barrel_nume_tbarW_50to70  *0.8, Nnorm_barrel_nume_tbarW_50to70  *1.2);
    RooRealVar n_endcap_nume_tbarW_50to70  ("n_endcap_nume_tbarW_50to70",   "n_endcap_nume_tbarW_50to70",   Nnorm_endcap_nume_tbarW_50to70 ,  Nnorm_endcap_nume_tbarW_50to70  *0.8, Nnorm_endcap_nume_tbarW_50to70  *1.2);
    RooRealVar n_barrel_deno_tbarW_50to70  ("n_barrel_deno_tbarW_50to70",   "n_barrel_deno_tbarW_50to70",   Nnorm_barrel_deno_tbarW_50to70 ,  Nnorm_barrel_deno_tbarW_50to70  *0.8, Nnorm_barrel_deno_tbarW_50to70  *1.2);
    RooRealVar n_endcap_deno_tbarW_50to70  ("n_endcap_deno_tbarW_50to70",   "n_endcap_deno_tbarW_50to70",   Nnorm_endcap_deno_tbarW_50to70 ,  Nnorm_endcap_deno_tbarW_50to70  *0.8, Nnorm_endcap_deno_tbarW_50to70  *1.2);
    RooRealVar n_barrel_nume_tbarW_70to100 ("n_barrel_nume_tbarW_70to100",  "n_barrel_nume_tbarW_70to100",  Nnorm_barrel_nume_tbarW_70to100,  Nnorm_barrel_nume_tbarW_70to100 *0.8, Nnorm_barrel_nume_tbarW_70to100 *1.2);
    RooRealVar n_endcap_nume_tbarW_70to100 ("n_endcap_nume_tbarW_70to100",  "n_endcap_nume_tbarW_70to100",  Nnorm_endcap_nume_tbarW_70to100,  Nnorm_endcap_nume_tbarW_70to100 *0.8, Nnorm_endcap_nume_tbarW_70to100 *1.2);
    RooRealVar n_barrel_deno_tbarW_70to100 ("n_barrel_deno_tbarW_70to100",  "n_barrel_deno_tbarW_70to100",  Nnorm_barrel_deno_tbarW_70to100,  Nnorm_barrel_deno_tbarW_70to100 *0.8, Nnorm_barrel_deno_tbarW_70to100 *1.2);
    RooRealVar n_endcap_deno_tbarW_70to100 ("n_endcap_deno_tbarW_70to100",  "n_endcap_deno_tbarW_70to100",  Nnorm_endcap_deno_tbarW_70to100,  Nnorm_endcap_deno_tbarW_70to100 *0.8, Nnorm_endcap_deno_tbarW_70to100 *1.2);
    RooRealVar n_barrel_nume_tbarW_100to500("n_barrel_nume_tbarW_100to500", "n_barrel_nume_tbarW_100to500", Nnorm_barrel_nume_tbarW_100to500, Nnorm_barrel_nume_tbarW_100to500*0.8, Nnorm_barrel_nume_tbarW_100to500*1.2);
    RooRealVar n_endcap_nume_tbarW_100to500("n_endcap_nume_tbarW_100to500", "n_endcap_nume_tbarW_100to500", Nnorm_endcap_nume_tbarW_100to500, Nnorm_endcap_nume_tbarW_100to500*0.8, Nnorm_endcap_nume_tbarW_100to500*1.2);
    RooRealVar n_barrel_deno_tbarW_100to500("n_barrel_deno_tbarW_100to500", "n_barrel_deno_tbarW_100to500", Nnorm_barrel_deno_tbarW_100to500, Nnorm_barrel_deno_tbarW_100to500*0.8, Nnorm_barrel_deno_tbarW_100to500*1.2);
    RooRealVar n_endcap_deno_tbarW_100to500("n_endcap_deno_tbarW_100to500", "n_endcap_deno_tbarW_100to500", Nnorm_endcap_deno_tbarW_100to500, Nnorm_endcap_deno_tbarW_100to500*0.8, Nnorm_endcap_deno_tbarW_100to500*1.2);

    RooRealVar n_barrel_nume_WW_50to70  ("n_barrel_nume_WW_50to70",   "n_barrel_nume_WW_50to70",   Nnorm_barrel_nume_WW_50to70 ,  Nnorm_barrel_nume_WW_50to70  *0.8, Nnorm_barrel_nume_WW_50to70  *1.2);
    RooRealVar n_endcap_nume_WW_50to70  ("n_endcap_nume_WW_50to70",   "n_endcap_nume_WW_50to70",   Nnorm_endcap_nume_WW_50to70 ,  Nnorm_endcap_nume_WW_50to70  *0.8, Nnorm_endcap_nume_WW_50to70  *1.2);
    RooRealVar n_barrel_deno_WW_50to70  ("n_barrel_deno_WW_50to70",   "n_barrel_deno_WW_50to70",   Nnorm_barrel_deno_WW_50to70 ,  Nnorm_barrel_deno_WW_50to70  *0.8, Nnorm_barrel_deno_WW_50to70  *1.2);
    RooRealVar n_endcap_deno_WW_50to70  ("n_endcap_deno_WW_50to70",   "n_endcap_deno_WW_50to70",   Nnorm_endcap_deno_WW_50to70 ,  Nnorm_endcap_deno_WW_50to70  *0.8, Nnorm_endcap_deno_WW_50to70  *1.2);
    RooRealVar n_barrel_nume_WW_70to100 ("n_barrel_nume_WW_70to100",  "n_barrel_nume_WW_70to100",  Nnorm_barrel_nume_WW_70to100,  Nnorm_barrel_nume_WW_70to100 *0.8, Nnorm_barrel_nume_WW_70to100 *1.2);
    RooRealVar n_endcap_nume_WW_70to100 ("n_endcap_nume_WW_70to100",  "n_endcap_nume_WW_70to100",  Nnorm_endcap_nume_WW_70to100,  Nnorm_endcap_nume_WW_70to100 *0.8, Nnorm_endcap_nume_WW_70to100 *1.2);
    RooRealVar n_barrel_deno_WW_70to100 ("n_barrel_deno_WW_70to100",  "n_barrel_deno_WW_70to100",  Nnorm_barrel_deno_WW_70to100,  Nnorm_barrel_deno_WW_70to100 *0.8, Nnorm_barrel_deno_WW_70to100 *1.2);
    RooRealVar n_endcap_deno_WW_70to100 ("n_endcap_deno_WW_70to100",  "n_endcap_deno_WW_70to100",  Nnorm_endcap_deno_WW_70to100,  Nnorm_endcap_deno_WW_70to100 *0.8, Nnorm_endcap_deno_WW_70to100 *1.2);
    RooRealVar n_barrel_nume_WW_100to500("n_barrel_nume_WW_100to500", "n_barrel_nume_WW_100to500", Nnorm_barrel_nume_WW_100to500, Nnorm_barrel_nume_WW_100to500*0.8, Nnorm_barrel_nume_WW_100to500*1.2);
    RooRealVar n_endcap_nume_WW_100to500("n_endcap_nume_WW_100to500", "n_endcap_nume_WW_100to500", Nnorm_endcap_nume_WW_100to500, Nnorm_endcap_nume_WW_100to500*0.8, Nnorm_endcap_nume_WW_100to500*1.2);
    RooRealVar n_barrel_deno_WW_100to500("n_barrel_deno_WW_100to500", "n_barrel_deno_WW_100to500", Nnorm_barrel_deno_WW_100to500, Nnorm_barrel_deno_WW_100to500*0.8, Nnorm_barrel_deno_WW_100to500*1.2);
    RooRealVar n_endcap_deno_WW_100to500("n_endcap_deno_WW_100to500", "n_endcap_deno_WW_100to500", Nnorm_endcap_deno_WW_100to500, Nnorm_endcap_deno_WW_100to500*0.8, Nnorm_endcap_deno_WW_100to500*1.2);

    RooRealVar n_barrel_nume_WZ_50to70  ("n_barrel_nume_WZ_50to70",   "n_barrel_nume_WZ_50to70",   Nnorm_barrel_nume_WZ_50to70 ,  Nnorm_barrel_nume_WZ_50to70  *0.8, Nnorm_barrel_nume_WZ_50to70  *1.2);
    RooRealVar n_endcap_nume_WZ_50to70  ("n_endcap_nume_WZ_50to70",   "n_endcap_nume_WZ_50to70",   Nnorm_endcap_nume_WZ_50to70 ,  Nnorm_endcap_nume_WZ_50to70  *0.8, Nnorm_endcap_nume_WZ_50to70  *1.2);
    RooRealVar n_barrel_deno_WZ_50to70  ("n_barrel_deno_WZ_50to70",   "n_barrel_deno_WZ_50to70",   Nnorm_barrel_deno_WZ_50to70 ,  Nnorm_barrel_deno_WZ_50to70  *0.8, Nnorm_barrel_deno_WZ_50to70  *1.2);
    RooRealVar n_endcap_deno_WZ_50to70  ("n_endcap_deno_WZ_50to70",   "n_endcap_deno_WZ_50to70",   Nnorm_endcap_deno_WZ_50to70 ,  Nnorm_endcap_deno_WZ_50to70  *0.8, Nnorm_endcap_deno_WZ_50to70  *1.2);
    RooRealVar n_barrel_nume_WZ_70to100 ("n_barrel_nume_WZ_70to100",  "n_barrel_nume_WZ_70to100",  Nnorm_barrel_nume_WZ_70to100,  Nnorm_barrel_nume_WZ_70to100 *0.8, Nnorm_barrel_nume_WZ_70to100 *1.2);
    RooRealVar n_endcap_nume_WZ_70to100 ("n_endcap_nume_WZ_70to100",  "n_endcap_nume_WZ_70to100",  Nnorm_endcap_nume_WZ_70to100,  Nnorm_endcap_nume_WZ_70to100 *0.8, Nnorm_endcap_nume_WZ_70to100 *1.2);
    RooRealVar n_barrel_deno_WZ_70to100 ("n_barrel_deno_WZ_70to100",  "n_barrel_deno_WZ_70to100",  Nnorm_barrel_deno_WZ_70to100,  Nnorm_barrel_deno_WZ_70to100 *0.8, Nnorm_barrel_deno_WZ_70to100 *1.2);
    RooRealVar n_endcap_deno_WZ_70to100 ("n_endcap_deno_WZ_70to100",  "n_endcap_deno_WZ_70to100",  Nnorm_endcap_deno_WZ_70to100,  Nnorm_endcap_deno_WZ_70to100 *0.8, Nnorm_endcap_deno_WZ_70to100 *1.2);
    RooRealVar n_barrel_nume_WZ_100to500("n_barrel_nume_WZ_100to500", "n_barrel_nume_WZ_100to500", Nnorm_barrel_nume_WZ_100to500, Nnorm_barrel_nume_WZ_100to500*0.8, Nnorm_barrel_nume_WZ_100to500*1.2);
    RooRealVar n_endcap_nume_WZ_100to500("n_endcap_nume_WZ_100to500", "n_endcap_nume_WZ_100to500", Nnorm_endcap_nume_WZ_100to500, Nnorm_endcap_nume_WZ_100to500*0.8, Nnorm_endcap_nume_WZ_100to500*1.2);
    RooRealVar n_barrel_deno_WZ_100to500("n_barrel_deno_WZ_100to500", "n_barrel_deno_WZ_100to500", Nnorm_barrel_deno_WZ_100to500, Nnorm_barrel_deno_WZ_100to500*0.8, Nnorm_barrel_deno_WZ_100to500*1.2);
    RooRealVar n_endcap_deno_WZ_100to500("n_endcap_deno_WZ_100to500", "n_endcap_deno_WZ_100to500", Nnorm_endcap_deno_WZ_100to500, Nnorm_endcap_deno_WZ_100to500*0.8, Nnorm_endcap_deno_WZ_100to500*1.2);

    RooRealVar n_barrel_nume_ZZ_50to70  ("n_barrel_nume_ZZ_50to70",   "n_barrel_nume_ZZ_50to70",   Nnorm_barrel_nume_ZZ_50to70 ,  Nnorm_barrel_nume_ZZ_50to70  *0.8, Nnorm_barrel_nume_ZZ_50to70  *1.2);
    RooRealVar n_endcap_nume_ZZ_50to70  ("n_endcap_nume_ZZ_50to70",   "n_endcap_nume_ZZ_50to70",   Nnorm_endcap_nume_ZZ_50to70 ,  Nnorm_endcap_nume_ZZ_50to70  *0.8, Nnorm_endcap_nume_ZZ_50to70  *1.2);
    RooRealVar n_barrel_deno_ZZ_50to70  ("n_barrel_deno_ZZ_50to70",   "n_barrel_deno_ZZ_50to70",   Nnorm_barrel_deno_ZZ_50to70 ,  Nnorm_barrel_deno_ZZ_50to70  *0.8, Nnorm_barrel_deno_ZZ_50to70  *1.2);
    RooRealVar n_endcap_deno_ZZ_50to70  ("n_endcap_deno_ZZ_50to70",   "n_endcap_deno_ZZ_50to70",   Nnorm_endcap_deno_ZZ_50to70 ,  Nnorm_endcap_deno_ZZ_50to70  *0.8, Nnorm_endcap_deno_ZZ_50to70  *1.2);
    RooRealVar n_barrel_nume_ZZ_70to100 ("n_barrel_nume_ZZ_70to100",  "n_barrel_nume_ZZ_70to100",  Nnorm_barrel_nume_ZZ_70to100,  Nnorm_barrel_nume_ZZ_70to100 *0.8, Nnorm_barrel_nume_ZZ_70to100 *1.2);
    RooRealVar n_endcap_nume_ZZ_70to100 ("n_endcap_nume_ZZ_70to100",  "n_endcap_nume_ZZ_70to100",  Nnorm_endcap_nume_ZZ_70to100,  Nnorm_endcap_nume_ZZ_70to100 *0.8, Nnorm_endcap_nume_ZZ_70to100 *1.2);
    RooRealVar n_barrel_deno_ZZ_70to100 ("n_barrel_deno_ZZ_70to100",  "n_barrel_deno_ZZ_70to100",  Nnorm_barrel_deno_ZZ_70to100,  Nnorm_barrel_deno_ZZ_70to100 *0.8, Nnorm_barrel_deno_ZZ_70to100 *1.2);
    RooRealVar n_endcap_deno_ZZ_70to100 ("n_endcap_deno_ZZ_70to100",  "n_endcap_deno_ZZ_70to100",  Nnorm_endcap_deno_ZZ_70to100,  Nnorm_endcap_deno_ZZ_70to100 *0.8, Nnorm_endcap_deno_ZZ_70to100 *1.2);
    RooRealVar n_barrel_nume_ZZ_100to500("n_barrel_nume_ZZ_100to500", "n_barrel_nume_ZZ_100to500", Nnorm_barrel_nume_ZZ_100to500, Nnorm_barrel_nume_ZZ_100to500*0.8, Nnorm_barrel_nume_ZZ_100to500*1.2);
    RooRealVar n_endcap_nume_ZZ_100to500("n_endcap_nume_ZZ_100to500", "n_endcap_nume_ZZ_100to500", Nnorm_endcap_nume_ZZ_100to500, Nnorm_endcap_nume_ZZ_100to500*0.8, Nnorm_endcap_nume_ZZ_100to500*1.2);
    RooRealVar n_barrel_deno_ZZ_100to500("n_barrel_deno_ZZ_100to500", "n_barrel_deno_ZZ_100to500", Nnorm_barrel_deno_ZZ_100to500, Nnorm_barrel_deno_ZZ_100to500*0.8, Nnorm_barrel_deno_ZZ_100to500*1.2);
    RooRealVar n_endcap_deno_ZZ_100to500("n_endcap_deno_ZZ_100to500", "n_endcap_deno_ZZ_100to500", Nnorm_endcap_deno_ZZ_100to500, Nnorm_endcap_deno_ZZ_100to500*0.8, Nnorm_endcap_deno_ZZ_100to500*1.2);

    RooRealVar n_barrel_nume_QCD_50to70  ("n_barrel_nume_QCD_50to70",   "n_barrel_nume_QCD_50to70",   Nnorm_barrel_nume_QCD_50to70 ,  Nnorm_barrel_nume_QCD_50to70  *0.5, Nnorm_barrel_nume_QCD_50to70  *1.5);
    RooRealVar n_endcap_nume_QCD_50to70  ("n_endcap_nume_QCD_50to70",   "n_endcap_nume_QCD_50to70",   Nnorm_endcap_nume_QCD_50to70 ,  Nnorm_endcap_nume_QCD_50to70  *0.5, Nnorm_endcap_nume_QCD_50to70  *1.5);
    RooRealVar n_barrel_deno_QCD_50to70  ("n_barrel_deno_QCD_50to70",   "n_barrel_deno_QCD_50to70",   Nnorm_barrel_deno_QCD_50to70 ,  Nnorm_barrel_deno_QCD_50to70  *0.5, Nnorm_barrel_deno_QCD_50to70  *1.5);
    RooRealVar n_endcap_deno_QCD_50to70  ("n_endcap_deno_QCD_50to70",   "n_endcap_deno_QCD_50to70",   Nnorm_endcap_deno_QCD_50to70 ,  Nnorm_endcap_deno_QCD_50to70  *0.5, Nnorm_endcap_deno_QCD_50to70  *1.5);
    RooRealVar n_barrel_nume_QCD_70to100 ("n_barrel_nume_QCD_70to100",  "n_barrel_nume_QCD_70to100",  Nnorm_barrel_nume_QCD_70to100,  Nnorm_barrel_nume_QCD_70to100 *0.5, Nnorm_barrel_nume_QCD_70to100 *1.5);
    RooRealVar n_endcap_nume_QCD_70to100 ("n_endcap_nume_QCD_70to100",  "n_endcap_nume_QCD_70to100",  Nnorm_endcap_nume_QCD_70to100,  Nnorm_endcap_nume_QCD_70to100 *0.5, Nnorm_endcap_nume_QCD_70to100 *1.5);
    RooRealVar n_barrel_deno_QCD_70to100 ("n_barrel_deno_QCD_70to100",  "n_barrel_deno_QCD_70to100",  Nnorm_barrel_deno_QCD_70to100,  Nnorm_barrel_deno_QCD_70to100 *0.5, Nnorm_barrel_deno_QCD_70to100 *1.5);
    RooRealVar n_endcap_deno_QCD_70to100 ("n_endcap_deno_QCD_70to100",  "n_endcap_deno_QCD_70to100",  Nnorm_endcap_deno_QCD_70to100,  Nnorm_endcap_deno_QCD_70to100 *0.5, Nnorm_endcap_deno_QCD_70to100 *1.5);
    RooRealVar n_barrel_nume_QCD_100to500("n_barrel_nume_QCD_100to500", "n_barrel_nume_QCD_100to500", Nnorm_barrel_nume_QCD_100to500, Nnorm_barrel_nume_QCD_100to500*0.5, Nnorm_barrel_nume_QCD_100to500*1.5);
    RooRealVar n_endcap_nume_QCD_100to500("n_endcap_nume_QCD_100to500", "n_endcap_nume_QCD_100to500", Nnorm_endcap_nume_QCD_100to500, Nnorm_endcap_nume_QCD_100to500*0.5, Nnorm_endcap_nume_QCD_100to500*1.5);
    RooRealVar n_barrel_deno_QCD_100to500("n_barrel_deno_QCD_100to500", "n_barrel_deno_QCD_100to500", Nnorm_barrel_deno_QCD_100to500, Nnorm_barrel_deno_QCD_100to500*0.5, Nnorm_barrel_deno_QCD_100to500*1.5);
    RooRealVar n_endcap_deno_QCD_100to500("n_endcap_deno_QCD_100to500", "n_endcap_deno_QCD_100to500", Nnorm_endcap_deno_QCD_100to500, Nnorm_endcap_deno_QCD_100to500*0.5, Nnorm_endcap_deno_QCD_100to500*1.5);

    // Models
    RooAddPdf model_barrel_nume_50to70("model_barrel_nume_50to70", "model_barrel_nume_50to70",
                                        RooArgList(*pdf_barrel_nume_QCD_50to70,   *pdf_barrel_nume_WJets_50to70, *pdf_barrel_nume_DY_50to70,
                                                   *pdf_barrel_nume_ttbar_50to70, *pdf_barrel_nume_tW_50to70,    *pdf_barrel_nume_tbarW_50to70,
                                                   *pdf_barrel_nume_WW_50to70,    *pdf_barrel_nume_WZ_50to70,    *pdf_barrel_nume_ZZ_50to70),
                                        RooArgList(n_barrel_nume_QCD_50to70,   n_barrel_nume_WJets_50to70, n_barrel_nume_DY_50to70,
                                                   n_barrel_nume_ttbar_50to70, n_barrel_nume_tW_50to70,    n_barrel_nume_tbarW_50to70,
                                                   n_barrel_nume_WW_50to70,    n_barrel_nume_WZ_50to70,    n_barrel_nume_ZZ_50to70));

    RooAddPdf model_endcap_nume_50to70("model_endcap_nume_50to70", "model_endcap_nume_50to70",
                                        RooArgList(*pdf_endcap_nume_QCD_50to70,   *pdf_endcap_nume_WJets_50to70, *pdf_endcap_nume_DY_50to70,
                                                   *pdf_endcap_nume_ttbar_50to70, *pdf_endcap_nume_tW_50to70,    *pdf_endcap_nume_tbarW_50to70,
                                                   *pdf_endcap_nume_WW_50to70,    *pdf_endcap_nume_WZ_50to70,    *pdf_endcap_nume_ZZ_50to70),
                                        RooArgList(n_endcap_nume_QCD_50to70,   n_endcap_nume_WJets_50to70, n_endcap_nume_DY_50to70,
                                                   n_endcap_nume_ttbar_50to70, n_endcap_nume_tW_50to70,    n_endcap_nume_tbarW_50to70,
                                                   n_endcap_nume_WW_50to70,    n_endcap_nume_WZ_50to70,    n_endcap_nume_ZZ_50to70));

    RooAddPdf model_barrel_deno_50to70("model_barrel_deno_50to70", "model_barrel_deno_50to70",
                                        RooArgList(*pdf_barrel_deno_QCD_50to70,   *pdf_barrel_deno_WJets_50to70, *pdf_barrel_deno_DY_50to70,
                                                   *pdf_barrel_deno_ttbar_50to70, *pdf_barrel_deno_tW_50to70,    *pdf_barrel_deno_tbarW_50to70,
                                                   *pdf_barrel_deno_WW_50to70,    *pdf_barrel_deno_WZ_50to70,    *pdf_barrel_deno_ZZ_50to70),
                                        RooArgList(n_barrel_deno_QCD_50to70,   n_barrel_deno_WJets_50to70, n_barrel_deno_DY_50to70,
                                                   n_barrel_deno_ttbar_50to70, n_barrel_deno_tW_50to70,    n_barrel_deno_tbarW_50to70,
                                                   n_barrel_deno_WW_50to70,    n_barrel_deno_WZ_50to70,    n_barrel_deno_ZZ_50to70));

    RooAddPdf model_endcap_deno_50to70("model_endcap_deno_50to70", "model_endcap_deno_50to70",
                                        RooArgList(*pdf_endcap_deno_QCD_50to70,   *pdf_endcap_deno_WJets_50to70, *pdf_endcap_deno_DY_50to70,
                                                   *pdf_endcap_deno_ttbar_50to70, *pdf_endcap_deno_tW_50to70,    *pdf_endcap_deno_tbarW_50to70,
                                                   *pdf_endcap_deno_WW_50to70,    *pdf_endcap_deno_WZ_50to70,    *pdf_endcap_deno_ZZ_50to70),
                                        RooArgList(n_endcap_deno_QCD_50to70,   n_endcap_deno_WJets_50to70, n_endcap_deno_DY_50to70,
                                                   n_endcap_deno_ttbar_50to70, n_endcap_deno_tW_50to70,    n_endcap_deno_tbarW_50to70,
                                                   n_endcap_deno_WW_50to70,    n_endcap_deno_WZ_50to70,    n_endcap_deno_ZZ_50to70));


    RooAddPdf model_barrel_nume_70to100("model_barrel_nume_70to100", "model_barrel_nume_70to100",
                                         RooArgList(*pdf_barrel_nume_QCD_70to100,   *pdf_barrel_nume_WJets_70to100, *pdf_barrel_nume_DY_70to100,
                                                    *pdf_barrel_nume_ttbar_70to100, *pdf_barrel_nume_tW_70to100,    *pdf_barrel_nume_tbarW_70to100,
                                                    *pdf_barrel_nume_WW_70to100,    *pdf_barrel_nume_WZ_70to100,    *pdf_barrel_nume_ZZ_70to100),
                                         RooArgList(n_barrel_nume_QCD_70to100,   n_barrel_nume_WJets_70to100, n_barrel_nume_DY_70to100,
                                                    n_barrel_nume_ttbar_70to100, n_barrel_nume_tW_70to100,    n_barrel_nume_tbarW_70to100,
                                                    n_barrel_nume_WW_70to100,    n_barrel_nume_WZ_70to100, n_barrel_nume_ZZ_70to100));

    RooAddPdf model_endcap_nume_70to100("model_endcap_nume_70to100", "model_endcap_nume_70to100",
                                         RooArgList(*pdf_endcap_nume_QCD_70to100,   *pdf_endcap_nume_WJets_70to100, *pdf_endcap_nume_DY_70to100,
                                                    *pdf_endcap_nume_ttbar_70to100, *pdf_endcap_nume_tW_70to100,    *pdf_endcap_nume_tbarW_70to100,
                                                    *pdf_endcap_nume_WW_70to100,    *pdf_endcap_nume_WZ_70to100,    *pdf_endcap_nume_ZZ_70to100),
                                         RooArgList(n_endcap_nume_QCD_70to100,   n_endcap_nume_WJets_70to100, n_endcap_nume_DY_70to100,
                                                    n_endcap_nume_ttbar_70to100, n_endcap_nume_tW_70to100,    n_endcap_nume_tbarW_70to100,
                                                    n_endcap_nume_WW_70to100,    n_endcap_nume_WZ_70to100,    n_endcap_nume_ZZ_70to100));

    RooAddPdf model_barrel_deno_70to100("model_barrel_deno_70to100", "model_barrel_deno_70to100",
                                         RooArgList(*pdf_barrel_deno_QCD_70to100,   *pdf_barrel_deno_WJets_70to100, *pdf_barrel_deno_DY_70to100,
                                                    *pdf_barrel_deno_ttbar_70to100, *pdf_barrel_deno_tW_70to100,    *pdf_barrel_deno_tbarW_70to100,
                                                    *pdf_barrel_deno_WW_70to100,    *pdf_barrel_deno_WZ_70to100,    *pdf_barrel_deno_ZZ_70to100),
                                         RooArgList(n_barrel_deno_QCD_70to100,   n_barrel_deno_WJets_70to100, n_barrel_deno_DY_70to100,
                                                    n_barrel_deno_ttbar_70to100, n_barrel_deno_tW_70to100,    n_barrel_deno_tbarW_70to100,
                                                    n_barrel_deno_WW_70to100,    n_barrel_deno_WZ_70to100,    n_barrel_deno_ZZ_70to100));

    RooAddPdf model_endcap_deno_70to100("model_endcap_deno_70to100", "model_endcap_deno_70to100",
                                         RooArgList(*pdf_endcap_deno_QCD_70to100,   *pdf_endcap_deno_WJets_70to100, *pdf_endcap_deno_DY_70to100,
                                                    *pdf_endcap_deno_ttbar_70to100, *pdf_endcap_deno_tW_70to100,    *pdf_endcap_deno_tbarW_70to100,
                                                    *pdf_endcap_deno_WW_70to100,    *pdf_endcap_deno_WZ_70to100,    *pdf_endcap_deno_ZZ_70to100),
                                         RooArgList(n_endcap_deno_QCD_70to100,   n_endcap_deno_WJets_70to100, n_endcap_deno_DY_70to100,
                                                    n_endcap_deno_ttbar_70to100, n_endcap_deno_tW_70to100,    n_endcap_deno_tbarW_70to100,
                                                    n_endcap_deno_WW_70to100,    n_endcap_deno_WZ_70to100,    n_endcap_deno_ZZ_70to100));

    RooAddPdf model_barrel_nume_100to500("model_barrel_nume_100to500", "model_barrel_nume_100to500",
                                         RooArgList(*pdf_barrel_nume_QCD_100to500,   *pdf_barrel_nume_WJets_100to500, *pdf_barrel_nume_DY_100to500,
                                                    *pdf_barrel_nume_ttbar_100to500, *pdf_barrel_nume_tW_100to500,    *pdf_barrel_nume_tbarW_100to500,
                                                    *pdf_barrel_nume_WW_100to500,    *pdf_barrel_nume_WZ_100to500,    *pdf_barrel_nume_ZZ_100to500),
                                         RooArgList(n_barrel_nume_QCD_100to500,   n_barrel_nume_WJets_100to500, n_barrel_nume_DY_100to500,
                                                    n_barrel_nume_ttbar_100to500, n_barrel_nume_tW_100to500, n_barrel_nume_tbarW_100to500,
                                                    n_barrel_nume_WW_100to500,    n_barrel_nume_WZ_100to500, n_barrel_nume_ZZ_100to500));

    RooAddPdf model_endcap_nume_100to500("model_endcap_nume_100to500", "model_endcap_nume_100to500",
                                         RooArgList(*pdf_endcap_nume_QCD_100to500,   *pdf_endcap_nume_WJets_100to500, *pdf_endcap_nume_DY_100to500,
                                                    *pdf_endcap_nume_ttbar_100to500, *pdf_endcap_nume_tW_100to500,    *pdf_endcap_nume_tbarW_100to500,
                                                    *pdf_endcap_nume_WW_100to500,    *pdf_endcap_nume_WZ_100to500,    *pdf_endcap_nume_ZZ_100to500),
                                         RooArgList(n_endcap_nume_QCD_100to500,   n_endcap_nume_WJets_100to500, n_endcap_nume_DY_100to500,
                                                    n_endcap_nume_ttbar_100to500, n_endcap_nume_tW_100to500,    n_endcap_nume_tbarW_100to500,
                                                    n_endcap_nume_WW_100to500,    n_endcap_nume_WZ_100to500,    n_endcap_nume_ZZ_100to500));

    RooAddPdf model_barrel_deno_100to500("model_barrel_deno_100to500", "model_barrel_deno_100to500",
                                         RooArgList(*pdf_barrel_deno_QCD_100to500,   *pdf_barrel_deno_WJets_100to500, *pdf_barrel_deno_DY_100to500,
                                                    *pdf_barrel_deno_ttbar_100to500, *pdf_barrel_deno_tW_100to500,    *pdf_barrel_deno_tbarW_100to500,
                                                    *pdf_barrel_deno_WW_100to500,    *pdf_barrel_deno_WZ_100to500,    *pdf_barrel_deno_ZZ_100to500),
                                         RooArgList(n_barrel_deno_QCD_100to500,   n_barrel_deno_WJets_100to500, n_barrel_deno_DY_100to500,
                                                    n_barrel_deno_ttbar_100to500, n_barrel_deno_tW_100to500,    n_barrel_deno_tbarW_100to500,
                                                    n_barrel_deno_WW_100to500,    n_barrel_deno_WZ_100to500,    n_barrel_deno_ZZ_100to500));

    RooAddPdf model_endcap_deno_100to500("model_endcap_deno_100to500", "model_endcap_deno_100to500",
                                         RooArgList(*pdf_endcap_deno_QCD_100to500,   *pdf_endcap_deno_WJets_100to500, *pdf_endcap_deno_DY_100to500,
                                                    *pdf_endcap_deno_ttbar_100to500, *pdf_endcap_deno_tW_100to500,    *pdf_endcap_deno_tbarW_100to500,
                                                    *pdf_endcap_deno_WW_100to500,    *pdf_endcap_deno_WZ_100to500,    *pdf_endcap_deno_ZZ_100to500),
                                         RooArgList(n_endcap_deno_QCD_100to500,   n_endcap_deno_WJets_100to500, n_endcap_deno_DY_100to500,
                                                    n_endcap_deno_ttbar_100to500, n_endcap_deno_tW_100to500,    n_endcap_deno_tbarW_100to500,
                                                    n_endcap_deno_WW_100to500,    n_endcap_deno_WZ_100to500,    n_endcap_deno_ZZ_100to500));

    // Fitting
    RooFitResult* fit_barrel_nume_50to70   = model_barrel_nume_50to70  .fitTo(*rh_barrel_nume_data_50to70  , Save());
    RooFitResult* fit_endcap_nume_50to70   = model_endcap_nume_50to70  .fitTo(*rh_endcap_nume_data_50to70  , Save());
    RooFitResult* fit_barrel_deno_50to70   = model_barrel_deno_50to70  .fitTo(*rh_barrel_deno_data_50to70  , Save());
    RooFitResult* fit_endcap_deno_50to70   = model_endcap_deno_50to70  .fitTo(*rh_endcap_deno_data_50to70  , Save());
    RooFitResult* fit_barrel_nume_70to100  = model_barrel_nume_70to100 .fitTo(*rh_barrel_nume_data_70to100 , Save());
    RooFitResult* fit_endcap_nume_70to100  = model_endcap_nume_70to100 .fitTo(*rh_endcap_nume_data_70to100 , Save());
    RooFitResult* fit_barrel_deno_70to100  = model_barrel_deno_70to100 .fitTo(*rh_barrel_deno_data_70to100 , Save());
    RooFitResult* fit_endcap_deno_70to100  = model_endcap_deno_70to100 .fitTo(*rh_endcap_deno_data_70to100 , Save());
    RooFitResult* fit_barrel_nume_100to500 = model_barrel_nume_100to500.fitTo(*rh_barrel_nume_data_100to500, Save());
    RooFitResult* fit_endcap_nume_100to500 = model_endcap_nume_100to500.fitTo(*rh_endcap_nume_data_100to500, Save());
    RooFitResult* fit_barrel_deno_100to500 = model_barrel_deno_100to500.fitTo(*rh_barrel_deno_data_100to500, Save());
    RooFitResult* fit_endcap_deno_100to500 = model_endcap_deno_100to500.fitTo(*rh_endcap_deno_data_100to500, Save());

    /// DRAWING NUMERATOR BARREL 50 to 70
    cout << "\n----- NUMERATOR BARREL pT 50to70 -----" << endl;
    TCanvas *c_fit_barrel_nume_50to70 = new TCanvas("c_fit_barrel_nume_50to70", "c_fit_barrel_nume_50to70", 800, 800);
    c_fit_barrel_nume_50to70->cd();

    //Top Pad
    TPad *c1_barrel_nume_50to70 = new TPad("padc1_barrel_nume_50to70","padc1_barrel_nume_50to70",0.01,0.01,0.99,0.99);
    c1_barrel_nume_50to70->Draw();
    c1_barrel_nume_50to70->cd();
    c1_barrel_nume_50to70->SetTopMargin(0.01);
    c1_barrel_nume_50to70->SetBottomMargin(0.35);
    c1_barrel_nume_50to70->SetRightMargin(0.03);
    c1_barrel_nume_50to70->SetLeftMargin(0.13);
    c1_barrel_nume_50to70->SetFillStyle(1);
    c1_barrel_nume_50to70->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_nume_50to70 = iso_nume.frame(Title(" "));
    rh_barrel_nume_data_50to70->plotOn(frame_barrel_nume_50to70, DataError(RooAbsData::SumW2));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70,pdf_barrel_nume_tbarW_50to70,pdf_barrel_nume_ttbar_50to70,"
                                                                         "pdf_barrel_nume_DY_50to70,pdf_barrel_nume_WJets_50to70,pdf_barrel_nume_QCD_50to70"),
                                    LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70,pdf_barrel_nume_tbarW_50to70,pdf_barrel_nume_ttbar_50to70,"
                                                                         "pdf_barrel_nume_DY_50to70,pdf_barrel_nume_WJets_50to70"),
                                    LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70,pdf_barrel_nume_tbarW_50to70,pdf_barrel_nume_ttbar_50to70,"
                                                                         "pdf_barrel_nume_DY_50to70"),
                                    LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70,pdf_barrel_nume_tbarW_50to70,pdf_barrel_nume_ttbar_50to70"),
                                    LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70,pdf_barrel_nume_tbarW_50to70"),
                                    LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70,"
                                                                         "pdf_barrel_nume_tW_50to70"),
                                    LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70,pdf_barrel_nume_WW_50to70"),
                                    LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70,pdf_barrel_nume_WZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_nume_50to70.plotOn(frame_barrel_nume_50to70, Components("pdf_barrel_nume_ZZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_nume_data_50to70->plotOn(frame_barrel_nume_50to70, DataError(RooAbsData::SumW2));
    frame_barrel_nume_50to70->Draw();
    fit_barrel_nume_50to70->Print();

    // Legend
    TLegend *legend = new TLegend(0.65, 0.7, 0.95, 0.97);
    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(0), "Data", "LP");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(1), "#font[12]{#scale[1.1]{QCD}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(2), "#font[12]{#scale[1.1]{W}}+Jets", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(3), "DY", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(4), "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(5), "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(6), "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(7), "#font[12]{#scale[1.1]{WW}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(8), "#font[12]{#scale[1.1]{WZ}}", "F");
    legend->AddEntry(frame_barrel_nume_50to70->nameOf(9), "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "F");
    legend->SetNColumns(2);

    legend->Draw();

    frame_barrel_nume_50to70->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_nume_50to70->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_nume_50to70->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_barrel_nume_50to70 = new TPad("padc2_barrel_nume_50to70","padc2_barrel_nume_50to70",0.01,0.01,0.99,0.35);
    c2_barrel_nume_50to70->Draw();
    c2_barrel_nume_50to70->cd();
    c2_barrel_nume_50to70->SetTopMargin(0.05);
    c2_barrel_nume_50to70->SetBottomMargin(0.33);
    c2_barrel_nume_50to70->SetRightMargin(0.02);
    c2_barrel_nume_50to70->SetLeftMargin(0.12);
    c2_barrel_nume_50to70->SetFillStyle(0);
    c2_barrel_nume_50to70->SetGrid();

    // Ratio plot
    TH1D *h_barrel_nume_MC_fit_50to70 = ((TH1D*)(model_barrel_nume_50to70.createHistogram("h_barrel_nume_MC_fit_50to70", iso_nume)));
    Double_t N_barrel_nume_data_50to70 = h_barrel_data_nume_50to70->Integral();
    Double_t N_barrel_nume_MC_50to70   = h_barrel_nume_MC_fit_50to70->Integral();
    h_barrel_nume_MC_fit_50to70->Scale(N_barrel_nume_data_50to70/N_barrel_nume_MC_50to70); // Why would I wanna do that???
    cout << "\nData integral: " << N_barrel_nume_data_50to70 << endl;
    cout << "MC integral: "     << h_barrel_nume_MC_fit_50to70->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_nume_50to70->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_nume_MC_fit_50to70->GetBinContent(1) << endl;

    TH1D *h_barrel_nume_ratio_50to70 = ((TH1D*)(h_barrel_data_nume_50to70->Clone("h_barrel_nume_ratio_50to70")));
    h_barrel_data_nume_50to70->Sumw2(); h_barrel_nume_MC_fit_50to70->Sumw2();
    h_barrel_nume_ratio_50to70->Divide(h_barrel_data_nume_50to70, h_barrel_nume_MC_fit_50to70);
    h_barrel_nume_ratio_50to70->SetTitle("");
    h_barrel_nume_ratio_50to70->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_nume_ratio_50to70->GetXaxis()->SetNoExponent(1);
    h_barrel_nume_ratio_50to70->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})");
    h_barrel_nume_ratio_50to70->GetXaxis()->SetTitleSize(0.17);
    h_barrel_nume_ratio_50to70->GetXaxis()->SetLabelSize(0.125);
    h_barrel_nume_ratio_50to70->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetTitle("Data/MC");
    h_barrel_nume_ratio_50to70->GetYaxis()->SetTitleSize(0.114);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetLabelSize(0.11);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetTickLength(0.01);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetDecimals(1);
    h_barrel_nume_ratio_50to70->SetMaximum(1.25);
    h_barrel_nume_ratio_50to70->SetMinimum(0.75);
    h_barrel_nume_ratio_50to70->GetYaxis()->SetNdivisions(5);
    h_barrel_nume_ratio_50to70->SetLineWidth(1);
    h_barrel_nume_ratio_50to70->SetLineColor(kBlack);
    h_barrel_nume_ratio_50to70->SetMarkerStyle(kFullDotLarge);
    h_barrel_nume_ratio_50to70->SetMarkerColor(kBlack);
    h_barrel_nume_ratio_50to70->SetStats(kFALSE);

    h_barrel_nume_ratio_50to70->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line = ((TH1D*)(h_barrel_nume_ratio_50to70->Clone("h_line")));
    h_line->Reset("ICES");
    for (Int_t i=1; i<=h_line->GetNbinsX(); i++)
        h_line->SetBinContent(i, 1);
    h_line->SetLineColor(kRed);
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_nume_50to70 = model_barrel_nume_50to70.createChi2(*rh_barrel_nume_data_50to70);
    cout << "chi2: " << chi2_barrel_nume_50to70->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_nume_50to70->getVal() / ((Double_t)h_barrel_data_nume_50to70->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING NUMERATOR ENDCAP 50to70
    cout << "\n----- NUMERATOR ENDCAP 50to70 -----" << endl;
    TCanvas *c_fit_endcap_nume_50to70 = new TCanvas("c_fit_endcap_nume_50to70", "c_fit_endcap_nume_50to70", 800, 800);
    c_fit_endcap_nume_50to70->cd();

    //Top Pad
    TPad *c1_endcap_nume_50to70 = new TPad("padc1_endcap_nume_50to70","padc1_endcap_nume_50to70",0.01,0.01,0.99,0.99);
    c1_endcap_nume_50to70->Draw();
    c1_endcap_nume_50to70->cd();
    c1_endcap_nume_50to70->SetTopMargin(0.01);
    c1_endcap_nume_50to70->SetBottomMargin(0.35);
    c1_endcap_nume_50to70->SetRightMargin(0.03);
    c1_endcap_nume_50to70->SetLeftMargin(0.13);
    c1_endcap_nume_50to70->SetFillStyle(1);
    c1_endcap_nume_50to70->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_nume_50to70 = iso_nume.frame(Title(" "));
    rh_endcap_nume_data_50to70->plotOn(frame_endcap_nume_50to70, DataError(RooAbsData::SumW2));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70,pdf_endcap_nume_tbarW_50to70,pdf_endcap_nume_ttbar_50to70,"
                                                                         "pdf_endcap_nume_DY_50to70,pdf_endcap_nume_WJets_50to70,pdf_endcap_nume_QCD_50to70"),
                                    LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70,pdf_endcap_nume_tbarW_50to70,pdf_endcap_nume_ttbar_50to70,"
                                                                         "pdf_endcap_nume_DY_50to70,pdf_endcap_nume_WJets_50to70"),
                                    LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70,pdf_endcap_nume_tbarW_50to70,pdf_endcap_nume_ttbar_50to70,"
                                                                         "pdf_endcap_nume_DY_50to70"),
                                    LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70,pdf_endcap_nume_tbarW_50to70,pdf_endcap_nume_ttbar_50to70"),
                                    LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70,pdf_endcap_nume_tbarW_50to70"),
                                    LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70,"
                                                                         "pdf_endcap_nume_tW_50to70"),
                                    LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70,pdf_endcap_nume_WW_50to70"),
                                    LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70,pdf_endcap_nume_WZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_nume_50to70.plotOn(frame_endcap_nume_50to70, Components("pdf_endcap_nume_ZZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_nume_data_50to70->plotOn(frame_endcap_nume_50to70, DataError(RooAbsData::SumW2));
    frame_endcap_nume_50to70->Draw();
    fit_endcap_nume_50to70->Print();

    // Legend
    legend->Draw();

    frame_endcap_nume_50to70->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_nume_50to70->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_nume_50to70->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_nume_50to70 = new TPad("padc2_endcap_nume_50to70","padc2_endcap_nume_50to70",0.01,0.01,0.99,0.35);
    c2_endcap_nume_50to70->Draw();
    c2_endcap_nume_50to70->cd();
    c2_endcap_nume_50to70->SetTopMargin(0.05);
    c2_endcap_nume_50to70->SetBottomMargin(0.33);
    c2_endcap_nume_50to70->SetRightMargin(0.02);
    c2_endcap_nume_50to70->SetLeftMargin(0.12);
    c2_endcap_nume_50to70->SetFillStyle(0);
    c2_endcap_nume_50to70->SetGrid();

    // Ratio plot
    TH1D *h_endcap_nume_MC_fit_50to70 = ((TH1D*)(model_endcap_nume_50to70.createHistogram("h_endcap_nume_MC_fit_50to70", iso_nume)));
    Double_t N_endcap_nume_data_50to70 = h_endcap_data_nume_50to70->Integral();
    Double_t N_endcap_nume_MC_50to70   = h_endcap_nume_MC_fit_50to70->Integral();
    h_endcap_nume_MC_fit_50to70->Scale(N_endcap_nume_data_50to70/N_endcap_nume_MC_50to70); // Why would I wanna do that???
    cout << "\nData integral: " << N_endcap_nume_data_50to70 << endl;
    cout << "MC integral: "     << h_endcap_nume_MC_fit_50to70->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_nume_50to70->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_nume_MC_fit_50to70->GetBinContent(1) << endl;

    TH1D *h_endcap_nume_ratio_50to70 = ((TH1D*)(h_endcap_data_nume_50to70->Clone("h_endcap_nume_ratio_50to70")));
    h_endcap_data_nume_50to70->Sumw2(); h_endcap_nume_MC_fit_50to70->Sumw2();
    h_endcap_nume_ratio_50to70->Divide(h_endcap_data_nume_50to70, h_endcap_nume_MC_fit_50to70);
    h_endcap_nume_ratio_50to70->SetTitle("");
    h_endcap_nume_ratio_50to70->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_nume_ratio_50to70->GetXaxis()->SetNoExponent(1);
    h_endcap_nume_ratio_50to70->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})");
    h_endcap_nume_ratio_50to70->GetXaxis()->SetTitleSize(0.17);
    h_endcap_nume_ratio_50to70->GetXaxis()->SetLabelSize(0.125);
    h_endcap_nume_ratio_50to70->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetTitle("Data/MC");
    h_endcap_nume_ratio_50to70->GetYaxis()->SetTitleSize(0.114);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetLabelSize(0.11);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetTickLength(0.01);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetDecimals(1);
    h_endcap_nume_ratio_50to70->SetMaximum(1.25);
    h_endcap_nume_ratio_50to70->SetMinimum(0.75);
    h_endcap_nume_ratio_50to70->GetYaxis()->SetNdivisions(5);
    h_endcap_nume_ratio_50to70->SetLineWidth(1);
    h_endcap_nume_ratio_50to70->SetLineColor(kBlack);
    h_endcap_nume_ratio_50to70->SetMarkerStyle(kFullDotLarge);
    h_endcap_nume_ratio_50to70->SetMarkerColor(kBlack);
    h_endcap_nume_ratio_50to70->SetStats(kFALSE);

    h_endcap_nume_ratio_50to70->Draw("E1P");

    // Red line at Data/MC=1
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_nume_50to70 = model_endcap_nume_50to70.createChi2(*rh_endcap_nume_data_50to70);
    cout << "chi2: " << chi2_endcap_nume_50to70->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_nume_50to70->getVal() / ((Double_t)h_endcap_data_nume_50to70->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR BARREL 50to70
    cout << "\n----- DENOMINATOR BARREL 50to70 -----" << endl;
    TCanvas *c_fit_barrel_deno_50to70 = new TCanvas("c_fit_barrel_deno_50to70", "c_fit_barrel_deno_50to70", 800, 800);
    c_fit_barrel_deno_50to70->cd();

    //Top Pad
    TPad *c1_barrel_deno_50to70 = new TPad("padc1_barrel_deno_50to70","padc1_barrel_deno_50to70",0.01,0.01,0.99,0.99);
    c1_barrel_deno_50to70->Draw();
    c1_barrel_deno_50to70->cd();
    c1_barrel_deno_50to70->SetTopMargin(0.01);
    c1_barrel_deno_50to70->SetBottomMargin(0.35);
    c1_barrel_deno_50to70->SetRightMargin(0.03);
    c1_barrel_deno_50to70->SetLeftMargin(0.13);
    c1_barrel_deno_50to70->SetFillStyle(1);
    c1_barrel_deno_50to70->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_deno_50to70 = iso_deno.frame(Title(" "));
    rh_barrel_deno_data_50to70->plotOn(frame_barrel_deno_50to70, DataError(RooAbsData::SumW2));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70,pdf_barrel_deno_tbarW_50to70,pdf_barrel_deno_ttbar_50to70,"
                                                                         "pdf_barrel_deno_DY_50to70,pdf_barrel_deno_WJets_50to70,pdf_barrel_deno_QCD_50to70"),
                                    LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70,pdf_barrel_deno_tbarW_50to70,pdf_barrel_deno_ttbar_50to70,"
                                                                         "pdf_barrel_deno_DY_50to70,pdf_barrel_deno_WJets_50to70"),
                                    LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70,pdf_barrel_deno_tbarW_50to70,pdf_barrel_deno_ttbar_50to70,"
                                                                         "pdf_barrel_deno_DY_50to70"),
                                    LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70,pdf_barrel_deno_tbarW_50to70,pdf_barrel_deno_ttbar_50to70"),
                                    LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70,pdf_barrel_deno_tbarW_50to70"),
                                    LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70,"
                                                                         "pdf_barrel_deno_tW_50to70"),
                                    LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70,pdf_barrel_deno_WW_50to70"),
                                    LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70,pdf_barrel_deno_WZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_deno_50to70.plotOn(frame_barrel_deno_50to70, Components("pdf_barrel_deno_ZZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_deno_data_50to70->plotOn(frame_barrel_deno_50to70, DataError(RooAbsData::SumW2));
    frame_barrel_deno_50to70->Draw();
    fit_barrel_deno_50to70->Print();

    // Legend
    legend->Draw();

    frame_barrel_deno_50to70->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_deno_50to70->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_deno_50to70->GetXaxis()->SetLabelSize(0);

    // Bottom Pad
    TPad *c2_barrel_deno_50to70 = new TPad("padc2_barrel_deno_50to70","padc2_barrel_deno_50to70",0.01,0.01,0.99,0.35);
    c2_barrel_deno_50to70->Draw();
    c2_barrel_deno_50to70->cd();
    c2_barrel_deno_50to70->SetTopMargin(0.05);
    c2_barrel_deno_50to70->SetBottomMargin(0.33);
    c2_barrel_deno_50to70->SetRightMargin(0.02);
    c2_barrel_deno_50to70->SetLeftMargin(0.12);
    c2_barrel_deno_50to70->SetFillStyle(0);
    c2_barrel_deno_50to70->SetGrid();

    // Ratio plot
    TH1D *h_barrel_deno_MC_fit_50to70 = ((TH1D*)(model_barrel_deno_50to70.createHistogram("h_barrel_deno_MC_fit_50to70", iso_deno)));
    Double_t N_barrel_deno_data_50to70 = h_barrel_data_deno_50to70->Integral();
    Double_t N_barrel_deno_MC_50to70   = h_barrel_deno_MC_fit_50to70->Integral();
    h_barrel_deno_MC_fit_50to70->Scale(N_barrel_deno_data_50to70/N_barrel_deno_MC_50to70); // Why is this necessary???
    cout << "\nData integral: " << N_barrel_deno_data_50to70 << endl;
    cout << "MC integral: "     << h_barrel_deno_MC_fit_50to70->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_deno_50to70->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_deno_MC_fit_50to70->GetBinContent(1) << endl;

    TH1D *h_barrel_deno_ratio_50to70 = ((TH1D*)(h_barrel_data_deno_50to70->Clone("h_barrel_deno_ratio_50to70")));
    h_barrel_data_deno_50to70->Sumw2(); h_barrel_deno_MC_fit_50to70->Sumw2();
    h_barrel_deno_ratio_50to70->Divide(h_barrel_data_deno_50to70, h_barrel_deno_MC_fit_50to70);
    h_barrel_deno_ratio_50to70->SetTitle("");
    h_barrel_deno_ratio_50to70->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_deno_ratio_50to70->GetXaxis()->SetNoExponent(1);
    h_barrel_deno_ratio_50to70->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_barrel_deno_ratio_50to70->GetXaxis()->SetTitleSize(0.17);
    h_barrel_deno_ratio_50to70->GetXaxis()->SetLabelSize(0.125);
    h_barrel_deno_ratio_50to70->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetTitle("Data/MC");
    h_barrel_deno_ratio_50to70->GetYaxis()->SetTitleSize(0.114);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetLabelSize(0.11);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetTickLength(0.01);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetDecimals(1);
    h_barrel_deno_ratio_50to70->SetMaximum(1.25);
    h_barrel_deno_ratio_50to70->SetMinimum(0.75);
    h_barrel_deno_ratio_50to70->GetYaxis()->SetNdivisions(5);
    h_barrel_deno_ratio_50to70->SetLineWidth(1);
    h_barrel_deno_ratio_50to70->SetLineColor(kBlack);
    h_barrel_deno_ratio_50to70->SetMarkerStyle(kFullDotLarge);
    h_barrel_deno_ratio_50to70->SetMarkerColor(kBlack);
    h_barrel_deno_ratio_50to70->SetStats(kFALSE);

    h_barrel_deno_ratio_50to70->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line_deno = ((TH1D*)(h_barrel_deno_ratio_50to70->Clone("h_line_deno")));
    h_line_deno->Reset("ICES");
    for (Int_t i=1; i<=h_line_deno->GetNbinsX(); i++)
        h_line_deno->SetBinContent(i, 1);
    h_line_deno->SetLineColor(kRed);
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_deno_50to70 = model_barrel_deno_50to70.createChi2(*rh_barrel_deno_data_50to70);
    cout << "chi2: " << chi2_barrel_deno_50to70->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_deno_50to70->getVal() / ((Double_t)h_barrel_data_deno_50to70->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR ENDCAP
    cout << "\n----- DENOMINATOR ENDCAP 50to70 -----" << endl;
    TCanvas *c_fit_endcap_deno_50to70 = new TCanvas("c_fit_endcap_deno_50to70", "c_fit_endcap_deno_50to70", 800, 800);
    c_fit_endcap_deno_50to70->cd();

    //Top Pad
    TPad *c1_endcap_deno_50to70 = new TPad("padc1_endcap_deno_50to70","padc1_endcap_deno_50to70",0.01,0.01,0.99,0.99);
    c1_endcap_deno_50to70->Draw();
    c1_endcap_deno_50to70->cd();
    c1_endcap_deno_50to70->SetTopMargin(0.01);
    c1_endcap_deno_50to70->SetBottomMargin(0.35);
    c1_endcap_deno_50to70->SetRightMargin(0.03);
    c1_endcap_deno_50to70->SetLeftMargin(0.13);
    c1_endcap_deno_50to70->SetFillStyle(1);
    c1_endcap_deno_50to70->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_deno_50to70 = iso_deno.frame(Title(" "));
    rh_endcap_deno_data_50to70->plotOn(frame_endcap_deno_50to70, DataError(RooAbsData::SumW2));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70,pdf_endcap_deno_tbarW_50to70,pdf_endcap_deno_ttbar_50to70,"
                                                                         "pdf_endcap_deno_DY_50to70,pdf_endcap_deno_WJets_50to70,pdf_endcap_deno_QCD_50to70"),
                                    LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70,pdf_endcap_deno_tbarW_50to70,pdf_endcap_deno_ttbar_50to70,"
                                                                         "pdf_endcap_deno_DY_50to70,pdf_endcap_deno_WJets_50to70"),
                                    LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70,pdf_endcap_deno_tbarW_50to70,pdf_endcap_deno_ttbar_50to70,"
                                                                         "pdf_endcap_deno_DY_50to70"),
                                    LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70,pdf_endcap_deno_tbarW_50to70,pdf_endcap_deno_ttbar_50to70"),
                                    LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70,pdf_endcap_deno_tbarW_50to70"),
                                    LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70,"
                                                                         "pdf_endcap_deno_tW_50to70"),
                                    LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70,pdf_endcap_deno_WW_50to70"),
                                    LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70,pdf_endcap_deno_WZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_deno_50to70.plotOn(frame_endcap_deno_50to70, Components("pdf_endcap_deno_ZZ_50to70"),
                                    LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_deno_data_50to70->plotOn(frame_endcap_deno_50to70, DataError(RooAbsData::SumW2));
    frame_endcap_deno_50to70->Draw();
    fit_endcap_deno_50to70->Print();

    // Legend
    legend->Draw();

    frame_endcap_deno_50to70->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_deno_50to70->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_deno_50to70->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_deno_50to70 = new TPad("padc2_endcap_deno_50to70","padc2_endcap_deno_50to70",0.01,0.01,0.99,0.35);
    c2_endcap_deno_50to70->Draw();
    c2_endcap_deno_50to70->cd();
    c2_endcap_deno_50to70->SetTopMargin(0.05);
    c2_endcap_deno_50to70->SetBottomMargin(0.33);
    c2_endcap_deno_50to70->SetRightMargin(0.02);
    c2_endcap_deno_50to70->SetLeftMargin(0.12);
    c2_endcap_deno_50to70->SetFillStyle(0);
    c2_endcap_deno_50to70->SetGrid();

    // Ratio plot
    TH1D *h_endcap_deno_MC_fit_50to70 = ((TH1D*)(model_endcap_deno_50to70.createHistogram("h_endcap_deno_MC_fit_50to70", iso_deno)));
    Double_t N_endcap_deno_data_50to70 = h_endcap_data_deno_50to70->Integral();
    Double_t N_endcap_deno_MC_50to70   = h_endcap_deno_MC_fit_50to70->Integral();
    h_endcap_deno_MC_fit_50to70->Scale(N_endcap_deno_data_50to70/N_endcap_deno_MC_50to70); // Why is this necessary???
    cout << "\nData integral: " << N_endcap_deno_data_50to70 << endl;
    cout << "MC integral: "     << h_endcap_deno_MC_fit_50to70->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_deno_50to70->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_deno_MC_fit_50to70->GetBinContent(1) << endl;

    TH1D *h_endcap_deno_ratio_50to70 = ((TH1D*)(h_endcap_data_deno_50to70->Clone("h_endcap_deno_ratio_50to70")));
    h_endcap_data_deno_50to70->Sumw2(); h_endcap_deno_MC_fit_50to70->Sumw2();
    h_endcap_deno_ratio_50to70->Divide(h_endcap_data_deno_50to70, h_endcap_deno_MC_fit_50to70);
    h_endcap_deno_ratio_50to70->SetTitle("");
    h_endcap_deno_ratio_50to70->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_deno_ratio_50to70->GetXaxis()->SetNoExponent(1);
    h_endcap_deno_ratio_50to70->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_endcap_deno_ratio_50to70->GetXaxis()->SetTitleSize(0.17);
    h_endcap_deno_ratio_50to70->GetXaxis()->SetLabelSize(0.125);
    h_endcap_deno_ratio_50to70->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetTitle("Data/MC");
    h_endcap_deno_ratio_50to70->GetYaxis()->SetTitleSize(0.114);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetLabelSize(0.11);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetTickLength(0.01);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetDecimals(1);
    h_endcap_deno_ratio_50to70->SetMaximum(1.25);
    h_endcap_deno_ratio_50to70->SetMinimum(0.75);
    h_endcap_deno_ratio_50to70->GetYaxis()->SetNdivisions(5);
    h_endcap_deno_ratio_50to70->SetLineWidth(1);
    h_endcap_deno_ratio_50to70->SetLineColor(kBlack);
    h_endcap_deno_ratio_50to70->SetMarkerStyle(kFullDotLarge);
    h_endcap_deno_ratio_50to70->SetMarkerColor(kBlack);
    h_endcap_deno_ratio_50to70->SetStats(kFALSE);

    h_endcap_deno_ratio_50to70->Draw("E1P");

    // Red line at Data/MC=1
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_deno_50to70 = model_endcap_deno_50to70.createChi2(*rh_endcap_deno_data_50to70);
    cout << "chi2: " << chi2_endcap_deno_50to70->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_deno_50to70->getVal() / ((Double_t)h_endcap_data_deno_50to70->GetNbinsX()) << "\n\n" << endl;


    /// DRAWING NUMERATOR BARREL 70 to 100
    cout << "\n----- NUMERATOR BARREL pT 70to100 -----" << endl;
    TCanvas *c_fit_barrel_nume_70to100 = new TCanvas("c_fit_barrel_nume_70to100", "c_fit_barrel_nume_70to100", 800, 800);
    c_fit_barrel_nume_70to100->cd();

    //Top Pad
    TPad *c1_barrel_nume_70to100 = new TPad("padc1_barrel_nume_70to100","padc1_barrel_nume_70to100",0.01,0.01,0.99,0.99);
    c1_barrel_nume_70to100->Draw();
    c1_barrel_nume_70to100->cd();
    c1_barrel_nume_70to100->SetTopMargin(0.01);
    c1_barrel_nume_70to100->SetBottomMargin(0.35);
    c1_barrel_nume_70to100->SetRightMargin(0.03);
    c1_barrel_nume_70to100->SetLeftMargin(0.13);
    c1_barrel_nume_70to100->SetFillStyle(1);
    c1_barrel_nume_70to100->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_nume_70to100 = iso_nume.frame(Title(" "));
    rh_barrel_nume_data_70to100->plotOn(frame_barrel_nume_70to100, DataError(RooAbsData::SumW2));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100,pdf_barrel_nume_tbarW_70to100,pdf_barrel_nume_ttbar_70to100,"
                                                                           "pdf_barrel_nume_DY_70to100,pdf_barrel_nume_WJets_70to100,pdf_barrel_nume_QCD_70to100"),
                                     LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100,pdf_barrel_nume_tbarW_70to100,pdf_barrel_nume_ttbar_70to100,"
                                                                           "pdf_barrel_nume_DY_70to100,pdf_barrel_nume_WJets_70to100"),
                                     LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100,pdf_barrel_nume_tbarW_70to100,pdf_barrel_nume_ttbar_70to100,"
                                                                           "pdf_barrel_nume_DY_70to100"),
                                     LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100,pdf_barrel_nume_tbarW_70to100,pdf_barrel_nume_ttbar_70to100"),
                                     LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100,pdf_barrel_nume_tbarW_70to100"),
                                     LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100,"
                                                                           "pdf_barrel_nume_tW_70to100"),
                                     LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100,pdf_barrel_nume_WW_70to100"),
                                     LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100,pdf_barrel_nume_WZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_nume_70to100.plotOn(frame_barrel_nume_70to100, Components("pdf_barrel_nume_ZZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_nume_data_70to100->plotOn(frame_barrel_nume_70to100, DataError(RooAbsData::SumW2));
    frame_barrel_nume_70to100->Draw();
    fit_barrel_nume_70to100->Print();

    legend->Draw();

    frame_barrel_nume_70to100->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_nume_70to100->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_nume_70to100->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_barrel_nume_70to100 = new TPad("padc2_barrel_nume_70to100","padc2_barrel_nume_70to100",0.01,0.01,0.99,0.35);
    c2_barrel_nume_70to100->Draw();
    c2_barrel_nume_70to100->cd();
    c2_barrel_nume_70to100->SetTopMargin(0.05);
    c2_barrel_nume_70to100->SetBottomMargin(0.33);
    c2_barrel_nume_70to100->SetRightMargin(0.02);
    c2_barrel_nume_70to100->SetLeftMargin(0.12);
    c2_barrel_nume_70to100->SetFillStyle(0);
    c2_barrel_nume_70to100->SetGrid();

    // Ratio plot
    TH1D *h_barrel_nume_MC_fit_70to100 = ((TH1D*)(model_barrel_nume_70to100.createHistogram("h_barrel_nume_MC_fit_70to100", iso_nume)));
    Double_t N_barrel_nume_data_70to100 = h_barrel_data_nume_70to100->Integral();
    Double_t N_barrel_nume_MC_70to100 = h_barrel_nume_MC_fit_70to100->Integral();
    h_barrel_nume_MC_fit_70to100->Scale(N_barrel_nume_data_70to100/N_barrel_nume_MC_70to100); // Why would I wanna do that???
    cout << "\nData integral: " << N_barrel_nume_data_70to100 << endl;
    cout << "MC integral: "     << h_barrel_nume_MC_fit_70to100->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_nume_70to100->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_nume_MC_fit_70to100->GetBinContent(1) << endl;

    TH1D *h_barrel_nume_ratio_70to100 = ((TH1D*)(h_barrel_data_nume_70to100->Clone("h_barrel_nume_ratio_70to100")));
    h_barrel_data_nume_70to100->Sumw2(); h_barrel_nume_MC_fit_70to100->Sumw2();
    h_barrel_nume_ratio_70to100->Divide(h_barrel_data_nume_70to100, h_barrel_nume_MC_fit_70to100);
    h_barrel_nume_ratio_70to100->SetTitle("");
    h_barrel_nume_ratio_70to100->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_nume_ratio_70to100->GetXaxis()->SetNoExponent(1);
    h_barrel_nume_ratio_70to100->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})");
    h_barrel_nume_ratio_70to100->GetXaxis()->SetTitleSize(0.17);
    h_barrel_nume_ratio_70to100->GetXaxis()->SetLabelSize(0.125);
    h_barrel_nume_ratio_70to100->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetTitle("Data/MC");
    h_barrel_nume_ratio_70to100->GetYaxis()->SetTitleSize(0.114);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetLabelSize(0.11);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetTickLength(0.01);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetDecimals(1);
    h_barrel_nume_ratio_70to100->SetMaximum(1.25);
    h_barrel_nume_ratio_70to100->SetMinimum(0.75);
    h_barrel_nume_ratio_70to100->GetYaxis()->SetNdivisions(5);
    h_barrel_nume_ratio_70to100->SetLineWidth(1);
    h_barrel_nume_ratio_70to100->SetLineColor(kBlack);
    h_barrel_nume_ratio_70to100->SetMarkerStyle(kFullDotLarge);
    h_barrel_nume_ratio_70to100->SetMarkerColor(kBlack);
    h_barrel_nume_ratio_70to100->SetStats(kFALSE);

    h_barrel_nume_ratio_70to100->Draw("E1P");
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_nume_70to100 = model_barrel_nume_70to100.createChi2(*rh_barrel_nume_data_70to100);
    cout << "chi2: " << chi2_barrel_nume_70to100->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_nume_70to100->getVal() / ((Double_t)h_barrel_data_nume_70to100->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING NUMERATOR ENDCAP 70to100
    cout << "\n----- NUMERATOR ENDCAP 70to100 -----" << endl;
    TCanvas *c_fit_endcap_nume_70to100 = new TCanvas("c_fit_endcap_nume_70to100", "c_fit_endcap_nume_70to100", 800, 800);
    c_fit_endcap_nume_70to100->cd();

    //Top Pad
    TPad *c1_endcap_nume_70to100 = new TPad("padc1_endcap_nume_70to100","padc1_endcap_nume_70to100",0.01,0.01,0.99,0.99);
    c1_endcap_nume_70to100->Draw();
    c1_endcap_nume_70to100->cd();
    c1_endcap_nume_70to100->SetTopMargin(0.01);
    c1_endcap_nume_70to100->SetBottomMargin(0.35);
    c1_endcap_nume_70to100->SetRightMargin(0.03);
    c1_endcap_nume_70to100->SetLeftMargin(0.13);
    c1_endcap_nume_70to100->SetFillStyle(1);
    c1_endcap_nume_70to100->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_nume_70to100 = iso_nume.frame(Title(" "));
    rh_endcap_nume_data_70to100->plotOn(frame_endcap_nume_70to100, DataError(RooAbsData::SumW2));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100,pdf_endcap_nume_tbarW_70to100,pdf_endcap_nume_ttbar_70to100,"
                                                                           "pdf_endcap_nume_DY_70to100,pdf_endcap_nume_WJets_70to100,pdf_endcap_nume_QCD_70to100"),
                                     LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100,pdf_endcap_nume_tbarW_70to100,pdf_endcap_nume_ttbar_70to100,"
                                                                           "pdf_endcap_nume_DY_70to100,pdf_endcap_nume_WJets_70to100"),
                                     LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100,pdf_endcap_nume_tbarW_70to100,pdf_endcap_nume_ttbar_70to100,"
                                                                           "pdf_endcap_nume_DY_70to100"),
                                     LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100,pdf_endcap_nume_tbarW_70to100,pdf_endcap_nume_ttbar_70to100"),
                                     LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100,pdf_endcap_nume_tbarW_70to100"),
                                     LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100,"
                                                                           "pdf_endcap_nume_tW_70to100"),
                                     LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100,pdf_endcap_nume_WW_70to100"),
                                     LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100,pdf_endcap_nume_WZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_nume_70to100.plotOn(frame_endcap_nume_70to100, Components("pdf_endcap_nume_ZZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_nume_data_70to100->plotOn(frame_endcap_nume_70to100, DataError(RooAbsData::SumW2));
    frame_endcap_nume_70to100->Draw();
    fit_endcap_nume_70to100->Print();

    // Legend
    legend->Draw();

    frame_endcap_nume_70to100->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_nume_70to100->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_nume_70to100->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_nume_70to100 = new TPad("padc2_endcap_nume_70to100","padc2_endcap_nume_70to100",0.01,0.01,0.99,0.35);
    c2_endcap_nume_70to100->Draw();
    c2_endcap_nume_70to100->cd();
    c2_endcap_nume_70to100->SetTopMargin(0.05);
    c2_endcap_nume_70to100->SetBottomMargin(0.33);
    c2_endcap_nume_70to100->SetRightMargin(0.02);
    c2_endcap_nume_70to100->SetLeftMargin(0.12);
    c2_endcap_nume_70to100->SetFillStyle(0);
    c2_endcap_nume_70to100->SetGrid();

    // Ratio plot
    TH1D *h_endcap_nume_MC_fit_70to100 = ((TH1D*)(model_endcap_nume_70to100.createHistogram("h_endcap_nume_MC_fit_70to100", iso_nume)));
    Double_t N_endcap_nume_data_70to100 = h_endcap_data_nume_70to100->Integral();
    Double_t N_endcap_nume_MC_70to100   = h_endcap_nume_MC_fit_70to100->Integral();
    h_endcap_nume_MC_fit_70to100->Scale(N_endcap_nume_data_70to100/N_endcap_nume_MC_70to100); // Why would I wanna do that???
    cout << "\nData integral: " << N_endcap_nume_data_70to100 << endl;
    cout << "MC integral: "     << h_endcap_nume_MC_fit_70to100->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_nume_70to100->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_nume_MC_fit_70to100->GetBinContent(1) << endl;

    TH1D *h_endcap_nume_ratio_70to100 = ((TH1D*)(h_endcap_data_nume_70to100->Clone("h_endcap_nume_ratio_70to100")));
    h_endcap_data_nume_70to100->Sumw2(); h_endcap_nume_MC_fit_70to100->Sumw2();
    h_endcap_nume_ratio_70to100->Divide(h_endcap_data_nume_70to100, h_endcap_nume_MC_fit_70to100);
    h_endcap_nume_ratio_70to100->SetTitle("");
    h_endcap_nume_ratio_70to100->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_nume_ratio_70to100->GetXaxis()->SetNoExponent(1);
    h_endcap_nume_ratio_70to100->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})");
    h_endcap_nume_ratio_70to100->GetXaxis()->SetTitleSize(0.17);
    h_endcap_nume_ratio_70to100->GetXaxis()->SetLabelSize(0.125);
    h_endcap_nume_ratio_70to100->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetTitle("Data/MC");
    h_endcap_nume_ratio_70to100->GetYaxis()->SetTitleSize(0.114);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetLabelSize(0.11);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetTickLength(0.01);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetDecimals(1);
    h_endcap_nume_ratio_70to100->SetMaximum(1.25);
    h_endcap_nume_ratio_70to100->SetMinimum(0.75);
    h_endcap_nume_ratio_70to100->GetYaxis()->SetNdivisions(5);
    h_endcap_nume_ratio_70to100->SetLineWidth(1);
    h_endcap_nume_ratio_70to100->SetLineColor(kBlack);
    h_endcap_nume_ratio_70to100->SetMarkerStyle(kFullDotLarge);
    h_endcap_nume_ratio_70to100->SetMarkerColor(kBlack);
    h_endcap_nume_ratio_70to100->SetStats(kFALSE);

    h_endcap_nume_ratio_70to100->Draw("E1P");

    // Red line at Data/MC=1
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_nume_70to100 = model_endcap_nume_70to100.createChi2(*rh_endcap_nume_data_70to100);
    cout << "chi2: " << chi2_endcap_nume_70to100->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_nume_70to100->getVal() / ((Double_t)h_endcap_data_nume_70to100->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR BARREL 70to100
    cout << "\n----- DENOMINATOR BARREL 70to100 -----" << endl;
    TCanvas *c_fit_barrel_deno_70to100 = new TCanvas("c_fit_barrel_deno_70to100", "c_fit_barrel_deno_70to100", 800, 800);
    c_fit_barrel_deno_70to100->cd();

    //Top Pad
    TPad *c1_barrel_deno_70to100 = new TPad("padc1_barrel_deno_70to100","padc1_barrel_deno_70to100",0.01,0.01,0.99,0.99);
    c1_barrel_deno_70to100->Draw();
    c1_barrel_deno_70to100->cd();
    c1_barrel_deno_70to100->SetTopMargin(0.01);
    c1_barrel_deno_70to100->SetBottomMargin(0.35);
    c1_barrel_deno_70to100->SetRightMargin(0.03);
    c1_barrel_deno_70to100->SetLeftMargin(0.13);
    c1_barrel_deno_70to100->SetFillStyle(1);
    c1_barrel_deno_70to100->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_deno_70to100 = iso_deno.frame(Title(" "));
    rh_barrel_deno_data_70to100->plotOn(frame_barrel_deno_70to100, DataError(RooAbsData::SumW2));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100,pdf_barrel_deno_tbarW_70to100,pdf_barrel_deno_ttbar_70to100,"
                                                                           "pdf_barrel_deno_DY_70to100,pdf_barrel_deno_WJets_70to100,pdf_barrel_deno_QCD_70to100"),
                                     LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100,pdf_barrel_deno_tbarW_70to100,pdf_barrel_deno_ttbar_70to100,"
                                                                           "pdf_barrel_deno_DY_70to100,pdf_barrel_deno_WJets_70to100"),
                                     LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100,pdf_barrel_deno_tbarW_70to100,pdf_barrel_deno_ttbar_70to100,"
                                                                           "pdf_barrel_deno_DY_70to100"),
                                     LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100,pdf_barrel_deno_tbarW_70to100,pdf_barrel_deno_ttbar_70to100"),
                                     LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100,pdf_barrel_deno_tbarW_70to100"),
                                     LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100,"
                                                                           "pdf_barrel_deno_tW_70to100"),
                                     LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100,pdf_barrel_deno_WW_70to100"),
                                     LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100,pdf_barrel_deno_WZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_deno_70to100.plotOn(frame_barrel_deno_70to100, Components("pdf_barrel_deno_ZZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_deno_data_70to100->plotOn(frame_barrel_deno_70to100, DataError(RooAbsData::SumW2));
    frame_barrel_deno_70to100->Draw();
    fit_barrel_deno_70to100->Print();

    // Legend
    legend->Draw();

    frame_barrel_deno_70to100->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_deno_70to100->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_deno_70to100->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_barrel_deno_70to100 = new TPad("padc2_barrel_deno_70to100","padc2_barrel_deno_70to100",0.01,0.01,0.99,0.35);
    c2_barrel_deno_70to100->Draw();
    c2_barrel_deno_70to100->cd();
    c2_barrel_deno_70to100->SetTopMargin(0.05);
    c2_barrel_deno_70to100->SetBottomMargin(0.33);
    c2_barrel_deno_70to100->SetRightMargin(0.02);
    c2_barrel_deno_70to100->SetLeftMargin(0.12);
    c2_barrel_deno_70to100->SetFillStyle(0);
    c2_barrel_deno_70to100->SetGrid();

    // Ratio plot
    TH1D *h_barrel_deno_MC_fit_70to100 = ((TH1D*)(model_barrel_deno_70to100.createHistogram("h_barrel_deno_MC_fit_70to100", iso_deno)));
    Double_t N_barrel_deno_data_70to100 = h_barrel_data_deno_70to100->Integral();
    Double_t N_barrel_deno_MC_70to100   = h_barrel_deno_MC_fit_70to100->Integral();
    h_barrel_deno_MC_fit_70to100->Scale(N_barrel_deno_data_70to100/N_barrel_deno_MC_70to100); // Why is this necessary???
    cout << "\nData integral: " << N_barrel_deno_data_70to100 << endl;
    cout << "MC integral: "     << h_barrel_deno_MC_fit_70to100->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_deno_70to100->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_deno_MC_fit_70to100->GetBinContent(1) << endl;

    TH1D *h_barrel_deno_ratio_70to100 = ((TH1D*)(h_barrel_data_deno_70to100->Clone("h_barrel_deno_ratio_70to100")));
    h_barrel_data_deno_70to100->Sumw2(); h_barrel_deno_MC_fit_70to100->Sumw2();
    h_barrel_deno_ratio_70to100->Divide(h_barrel_data_deno_70to100, h_barrel_deno_MC_fit_70to100);
    h_barrel_deno_ratio_70to100->SetTitle("");
    h_barrel_deno_ratio_70to100->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_deno_ratio_70to100->GetXaxis()->SetNoExponent(1);
    h_barrel_deno_ratio_70to100->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_barrel_deno_ratio_70to100->GetXaxis()->SetTitleSize(0.17);
    h_barrel_deno_ratio_70to100->GetXaxis()->SetLabelSize(0.125);
    h_barrel_deno_ratio_70to100->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetTitle("Data/MC");
    h_barrel_deno_ratio_70to100->GetYaxis()->SetTitleSize(0.114);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetLabelSize(0.11);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetTickLength(0.01);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetDecimals(1);
    h_barrel_deno_ratio_70to100->SetMaximum(1.25);
    h_barrel_deno_ratio_70to100->SetMinimum(0.75);
    h_barrel_deno_ratio_70to100->GetYaxis()->SetNdivisions(5);
    h_barrel_deno_ratio_70to100->SetLineWidth(1);
    h_barrel_deno_ratio_70to100->SetLineColor(kBlack);
    h_barrel_deno_ratio_70to100->SetMarkerStyle(kFullDotLarge);
    h_barrel_deno_ratio_70to100->SetMarkerColor(kBlack);
    h_barrel_deno_ratio_70to100->SetStats(kFALSE);

    h_barrel_deno_ratio_70to100->Draw("E1P");
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_deno_70to100 = model_barrel_deno_70to100.createChi2(*rh_barrel_deno_data_70to100);
    cout << "chi2: " << chi2_barrel_deno_70to100->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_deno_70to100->getVal() / ((Double_t)h_barrel_data_deno_70to100->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR ENDCAP 70to100
    cout << "\n----- DENOMINATOR ENDCAP 70to100 -----" << endl;
    TCanvas *c_fit_endcap_deno_70to100 = new TCanvas("c_fit_endcap_deno_70to100", "c_fit_endcap_deno_70to100", 800, 800);
    c_fit_endcap_deno_70to100->cd();

    //Top Pad
    TPad *c1_endcap_deno_70to100 = new TPad("padc1_endcap_deno_70to100","padc1_endcap_deno_70to100",0.01,0.01,0.99,0.99);
    c1_endcap_deno_70to100->Draw();
    c1_endcap_deno_70to100->cd();
    c1_endcap_deno_70to100->SetTopMargin(0.01);
    c1_endcap_deno_70to100->SetBottomMargin(0.35);
    c1_endcap_deno_70to100->SetRightMargin(0.03);
    c1_endcap_deno_70to100->SetLeftMargin(0.13);
    c1_endcap_deno_70to100->SetFillStyle(1);
    c1_endcap_deno_70to100->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_deno_70to100 = iso_deno.frame(Title(" "));
    rh_endcap_deno_data_70to100->plotOn(frame_endcap_deno_70to100, DataError(RooAbsData::SumW2));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100,pdf_endcap_deno_tbarW_70to100,pdf_endcap_deno_ttbar_70to100,"
                                                                           "pdf_endcap_deno_DY_70to100,pdf_endcap_deno_WJets_70to100,pdf_endcap_deno_QCD_70to100"),
                                     LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100,pdf_endcap_deno_tbarW_70to100,pdf_endcap_deno_ttbar_70to100,"
                                                                           "pdf_endcap_deno_DY_70to100,pdf_endcap_deno_WJets_70to100"),
                                     LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100,pdf_endcap_deno_tbarW_70to100,pdf_endcap_deno_ttbar_70to100,"
                                                                           "pdf_endcap_deno_DY_70to100"),
                                     LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100,pdf_endcap_deno_tbarW_70to100,pdf_endcap_deno_ttbar_70to100"),
                                     LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100,pdf_endcap_deno_tbarW_70to100"),
                                     LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100,"
                                                                           "pdf_endcap_deno_tW_70to100"),
                                     LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100,pdf_endcap_deno_WW_70to100"),
                                     LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100,pdf_endcap_deno_WZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_deno_70to100.plotOn(frame_endcap_deno_70to100, Components("pdf_endcap_deno_ZZ_70to100"),
                                     LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_deno_data_70to100->plotOn(frame_endcap_deno_70to100, DataError(RooAbsData::SumW2));
    frame_endcap_deno_70to100->Draw();
    fit_endcap_deno_70to100->Print();

    // Legend
    legend->Draw();

    frame_endcap_deno_70to100->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_deno_70to100->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_deno_70to100->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_deno_70to100 = new TPad("padc2_endcap_deno_70to100","padc2_endcap_deno_70to100",0.01,0.01,0.99,0.35);
    c2_endcap_deno_70to100->Draw();
    c2_endcap_deno_70to100->cd();
    c2_endcap_deno_70to100->SetTopMargin(0.05);
    c2_endcap_deno_70to100->SetBottomMargin(0.33);
    c2_endcap_deno_70to100->SetRightMargin(0.02);
    c2_endcap_deno_70to100->SetLeftMargin(0.12);
    c2_endcap_deno_70to100->SetFillStyle(0);
    c2_endcap_deno_70to100->SetGrid();

    // Ratio plot
    TH1D *h_endcap_deno_MC_fit_70to100 = ((TH1D*)(model_endcap_deno_70to100.createHistogram("h_endcap_deno_MC_fit_70to100", iso_deno)));
    Double_t N_endcap_deno_data_70to100 = h_endcap_data_deno_70to100->Integral();
    Double_t N_endcap_deno_MC_70to100   = h_endcap_deno_MC_fit_70to100->Integral();
    h_endcap_deno_MC_fit_70to100->Scale(N_endcap_deno_data_70to100/N_endcap_deno_MC_70to100); // Why is this necessary???
    cout << "\nData integral: " << N_endcap_deno_data_70to100 << endl;
    cout << "MC integral: "     << h_endcap_deno_MC_fit_70to100->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_deno_70to100->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_deno_MC_fit_70to100->GetBinContent(1) << endl;

    TH1D *h_endcap_deno_ratio_70to100 = ((TH1D*)(h_endcap_data_deno_70to100->Clone("h_endcap_deno_ratio_70to100")));
    h_endcap_data_deno_70to100->Sumw2(); h_endcap_deno_MC_fit_70to100->Sumw2();
    h_endcap_deno_ratio_70to100->Divide(h_endcap_data_deno_70to100, h_endcap_deno_MC_fit_70to100);
    h_endcap_deno_ratio_70to100->SetTitle("");
    h_endcap_deno_ratio_70to100->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_deno_ratio_70to100->GetXaxis()->SetNoExponent(1);
    h_endcap_deno_ratio_70to100->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_endcap_deno_ratio_70to100->GetXaxis()->SetTitleSize(0.17);
    h_endcap_deno_ratio_70to100->GetXaxis()->SetLabelSize(0.125);
    h_endcap_deno_ratio_70to100->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetTitle("Data/MC");
    h_endcap_deno_ratio_70to100->GetYaxis()->SetTitleSize(0.114);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetLabelSize(0.11);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetTickLength(0.01);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetDecimals(1);
    h_endcap_deno_ratio_70to100->SetMaximum(1.25);
    h_endcap_deno_ratio_70to100->SetMinimum(0.75);
    h_endcap_deno_ratio_70to100->GetYaxis()->SetNdivisions(5);
    h_endcap_deno_ratio_70to100->SetLineWidth(1);
    h_endcap_deno_ratio_70to100->SetLineColor(kBlack);
    h_endcap_deno_ratio_70to100->SetMarkerStyle(kFullDotLarge);
    h_endcap_deno_ratio_70to100->SetMarkerColor(kBlack);
    h_endcap_deno_ratio_70to100->SetStats(kFALSE);

    h_endcap_deno_ratio_70to100->Draw("E1P");

    // Red line at Data/MC=1
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_deno_70to100 = model_endcap_deno_70to100.createChi2(*rh_endcap_deno_data_70to100);
    cout << "chi2: " << chi2_endcap_deno_70to100->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_deno_70to100->getVal() / ((Double_t)h_endcap_data_deno_70to100->GetNbinsX()) << "\n\n" << endl;


    /// DRAWING NUMERATOR BARREL 200 to 500
    cout << "\n----- NUMERATOR BARREL pT 100to500 -----" << endl;
    TCanvas *c_fit_barrel_nume_100to500 = new TCanvas("c_fit_barrel_nume_100to500", "c_fit_barrel_nume_100to500", 800, 800);
    c_fit_barrel_nume_100to500->cd();

    //Top Pad
    TPad *c1_barrel_nume_100to500 = new TPad("padc1_barrel_nume_100to500","padc1_barrel_nume_100to500",0.01,0.01,0.99,0.99);
    c1_barrel_nume_100to500->Draw();
    c1_barrel_nume_100to500->cd();
    c1_barrel_nume_100to500->SetTopMargin(0.01);
    c1_barrel_nume_100to500->SetBottomMargin(0.35);
    c1_barrel_nume_100to500->SetRightMargin(0.03);
    c1_barrel_nume_100to500->SetLeftMargin(0.13);
    c1_barrel_nume_100to500->SetFillStyle(1);
    c1_barrel_nume_100to500->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_nume_100to500 = iso_nume.frame(Title(" "));
    rh_barrel_nume_data_100to500->plotOn(frame_barrel_nume_100to500, DataError(RooAbsData::SumW2));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500,pdf_barrel_nume_tbarW_100to500,pdf_barrel_nume_ttbar_100to500,"
                                                                             "pdf_barrel_nume_DY_100to500,pdf_barrel_nume_WJets_100to500,pdf_barrel_nume_QCD_100to500"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500,pdf_barrel_nume_tbarW_100to500,pdf_barrel_nume_ttbar_100to500,"
                                                                             "pdf_barrel_nume_DY_100to500,pdf_barrel_nume_WJets_100to500"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500,pdf_barrel_nume_tbarW_100to500,pdf_barrel_nume_ttbar_100to500,"
                                                                             "pdf_barrel_nume_DY_100to500"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500,pdf_barrel_nume_tbarW_100to500,pdf_barrel_nume_ttbar_100to500"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500,pdf_barrel_nume_tbarW_100to500"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500,"
                                                                             "pdf_barrel_nume_tW_100to500"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500,pdf_barrel_nume_WW_100to500"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500,pdf_barrel_nume_WZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_nume_100to500.plotOn(frame_barrel_nume_100to500, Components("pdf_barrel_nume_ZZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_nume_data_100to500->plotOn(frame_barrel_nume_100to500, DataError(RooAbsData::SumW2));
    frame_barrel_nume_100to500->Draw();
    fit_barrel_nume_100to500->Print();

    legend->Draw();

    frame_barrel_nume_100to500->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_nume_100to500->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_nume_100to500->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_barrel_nume_100to500 = new TPad("padc2_barrel_nume_100to500","padc2_barrel_nume_100to500",0.01,0.01,0.99,0.35);
    c2_barrel_nume_100to500->Draw();
    c2_barrel_nume_100to500->cd();
    c2_barrel_nume_100to500->SetTopMargin(0.05);
    c2_barrel_nume_100to500->SetBottomMargin(0.33);
    c2_barrel_nume_100to500->SetRightMargin(0.02);
    c2_barrel_nume_100to500->SetLeftMargin(0.12);
    c2_barrel_nume_100to500->SetFillStyle(0);
    c2_barrel_nume_100to500->SetGrid();

    // Ratio plot
    TH1D *h_barrel_nume_MC_fit_100to500 = ((TH1D*)(model_barrel_nume_100to500.createHistogram("h_barrel_nume_MC_fit_100to500", iso_nume)));
    Double_t N_barrel_nume_data_100to500 = h_barrel_data_nume_100to500->Integral();
    Double_t N_barrel_nume_MC_100to500 = h_barrel_nume_MC_fit_100to500->Integral();
    h_barrel_nume_MC_fit_100to500->Scale(N_barrel_nume_data_100to500/N_barrel_nume_MC_100to500); // Why would I wanna do that???
    cout << "\nData integral: " << N_barrel_nume_data_100to500 << endl;
    cout << "MC integral: "     << h_barrel_nume_MC_fit_100to500->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_nume_100to500->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_nume_MC_fit_100to500->GetBinContent(1) << endl;

    TH1D *h_barrel_nume_ratio_100to500 = ((TH1D*)(h_barrel_data_nume_100to500->Clone("h_barrel_nume_ratio_100to500")));
    h_barrel_data_nume_100to500->Sumw2(); h_barrel_nume_MC_fit_100to500->Sumw2();
    h_barrel_nume_ratio_100to500->Divide(h_barrel_data_nume_100to500, h_barrel_nume_MC_fit_100to500);
    h_barrel_nume_ratio_100to500->SetTitle("");
    h_barrel_nume_ratio_100to500->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_nume_ratio_100to500->GetXaxis()->SetNoExponent(1);
    h_barrel_nume_ratio_100to500->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})");
    h_barrel_nume_ratio_100to500->GetXaxis()->SetTitleSize(0.17);
    h_barrel_nume_ratio_100to500->GetXaxis()->SetLabelSize(0.125);
    h_barrel_nume_ratio_100to500->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetTitle("Data/MC");
    h_barrel_nume_ratio_100to500->GetYaxis()->SetTitleSize(0.114);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetLabelSize(0.11);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetTickLength(0.01);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetDecimals(1);
    h_barrel_nume_ratio_100to500->SetMaximum(1.25);
    h_barrel_nume_ratio_100to500->SetMinimum(0.75);
    h_barrel_nume_ratio_100to500->GetYaxis()->SetNdivisions(5);
    h_barrel_nume_ratio_100to500->SetLineWidth(1);
    h_barrel_nume_ratio_100to500->SetLineColor(kBlack);
    h_barrel_nume_ratio_100to500->SetMarkerStyle(kFullDotLarge);
    h_barrel_nume_ratio_100to500->SetMarkerColor(kBlack);
    h_barrel_nume_ratio_100to500->SetStats(kFALSE);

    h_barrel_nume_ratio_100to500->Draw("E1P");
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_nume_100to500 = model_barrel_nume_100to500.createChi2(*rh_barrel_nume_data_100to500);
    cout << "chi2: " << chi2_barrel_nume_100to500->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_nume_100to500->getVal() / ((Double_t)h_barrel_data_nume_100to500->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING NUMERATOR ENDCAP 100to500
    cout << "\n----- NUMERATOR ENDCAP 100to500 -----" << endl;
    TCanvas *c_fit_endcap_nume_100to500 = new TCanvas("c_fit_endcap_nume_100to500", "c_fit_endcap_nume_100to500", 800, 800);
    c_fit_endcap_nume_100to500->cd();

    //Top Pad
    TPad *c1_endcap_nume_100to500 = new TPad("padc1_endcap_nume_100to500","padc1_endcap_nume_100to500",0.01,0.01,0.99,0.99);
    c1_endcap_nume_100to500->Draw();
    c1_endcap_nume_100to500->cd();
    c1_endcap_nume_100to500->SetTopMargin(0.01);
    c1_endcap_nume_100to500->SetBottomMargin(0.35);
    c1_endcap_nume_100to500->SetRightMargin(0.03);
    c1_endcap_nume_100to500->SetLeftMargin(0.13);
    c1_endcap_nume_100to500->SetFillStyle(1);
    c1_endcap_nume_100to500->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_nume_100to500 = iso_nume.frame(Title(" "));
    rh_endcap_nume_data_100to500->plotOn(frame_endcap_nume_100to500, DataError(RooAbsData::SumW2));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                             "pdf_endcap_nume_tW_100to500,pdf_endcap_nume_tbarW_100to500,pdf_endcap_nume_ttbar_100to500,"
                                                                             "pdf_endcap_nume_DY_100to500,pdf_endcap_nume_WJets_100to500,pdf_endcap_nume_QCD_100to500"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                             "pdf_endcap_nume_tW_100to500,pdf_endcap_nume_tbarW_100to500,pdf_endcap_nume_ttbar_100to500,"
                                                                             "pdf_endcap_nume_DY_100to500,pdf_endcap_nume_WJets_100to500"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                             "pdf_endcap_nume_tW_100to500,pdf_endcap_nume_tbarW_100to500,pdf_endcap_nume_ttbar_100to500,"
                                                                             "pdf_endcap_nume_DY_100to500"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                             "pdf_endcap_nume_tW_100to500,pdf_endcap_nume_tbarW_100to500,pdf_endcap_nume_ttbar_100to500"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                             "pdf_endcap_nume_tW_100to500,pdf_endcap_nume_tbarW_100to500"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500,"
                                                                            "pdf_endcap_nume_tW_100to500"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500,pdf_endcap_nume_WW_100to500"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500,pdf_endcap_nume_WZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_nume_100to500.plotOn(frame_endcap_nume_100to500, Components("pdf_endcap_nume_ZZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_nume_data_100to500->plotOn(frame_endcap_nume_100to500, DataError(RooAbsData::SumW2));
    frame_endcap_nume_100to500->Draw();
    fit_endcap_nume_100to500->Print();

    // Legend
    legend->Draw();

    frame_endcap_nume_100to500->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_nume_100to500->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_nume_100to500->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_nume_100to500 = new TPad("padc2_endcap_nume_100to500","padc2_endcap_nume_100to500",0.01,0.01,0.99,0.35);
    c2_endcap_nume_100to500->Draw();
    c2_endcap_nume_100to500->cd();
    c2_endcap_nume_100to500->SetTopMargin(0.05);
    c2_endcap_nume_100to500->SetBottomMargin(0.33);
    c2_endcap_nume_100to500->SetRightMargin(0.02);
    c2_endcap_nume_100to500->SetLeftMargin(0.12);
    c2_endcap_nume_100to500->SetFillStyle(0);
    c2_endcap_nume_100to500->SetGrid();

    // Ratio plot
    TH1D *h_endcap_nume_MC_fit_100to500 = ((TH1D*)(model_endcap_nume_100to500.createHistogram("h_endcap_nume_MC_fit_100to500", iso_nume)));
    Double_t N_endcap_nume_data_100to500 = h_endcap_data_nume_100to500->Integral();
    Double_t N_endcap_nume_MC_100to500   = h_endcap_nume_MC_fit_100to500->Integral();
    h_endcap_nume_MC_fit_100to500->Scale(N_endcap_nume_data_100to500/N_endcap_nume_MC_100to500); // Why would I wanna do that???
    cout << "\nData integral: " << N_endcap_nume_data_100to500 << endl;
    cout << "MC integral: "     << h_endcap_nume_MC_fit_100to500->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_nume_100to500->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_nume_MC_fit_100to500->GetBinContent(1) << endl;

    TH1D *h_endcap_nume_ratio_100to500 = ((TH1D*)(h_endcap_data_nume_100to500->Clone("h_endcap_nume_ratio_100to500")));
    h_endcap_data_nume_100to500->Sumw2(); h_endcap_nume_MC_fit_100to500->Sumw2();
    h_endcap_nume_ratio_100to500->Divide(h_endcap_data_nume_100to500, h_endcap_nume_MC_fit_100to500);
    h_endcap_nume_ratio_100to500->SetTitle("");
    h_endcap_nume_ratio_100to500->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_nume_ratio_100to500->GetXaxis()->SetNoExponent(1);
    h_endcap_nume_ratio_100to500->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})");
    h_endcap_nume_ratio_100to500->GetXaxis()->SetTitleSize(0.17);
    h_endcap_nume_ratio_100to500->GetXaxis()->SetLabelSize(0.125);
    h_endcap_nume_ratio_100to500->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetTitle("Data/MC");
    h_endcap_nume_ratio_100to500->GetYaxis()->SetTitleSize(0.114);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetLabelSize(0.11);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetTickLength(0.01);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetDecimals(1);
    h_endcap_nume_ratio_100to500->SetMaximum(1.25);
    h_endcap_nume_ratio_100to500->SetMinimum(0.75);
    h_endcap_nume_ratio_100to500->GetYaxis()->SetNdivisions(5);
    h_endcap_nume_ratio_100to500->SetLineWidth(1);
    h_endcap_nume_ratio_100to500->SetLineColor(kBlack);
    h_endcap_nume_ratio_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_nume_ratio_100to500->SetMarkerColor(kBlack);
    h_endcap_nume_ratio_100to500->SetStats(kFALSE);

    h_endcap_nume_ratio_100to500->Draw("E1P");

    // Red line at Data/MC=1
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_nume_100to500 = model_endcap_nume_100to500.createChi2(*rh_endcap_nume_data_100to500);
    cout << "chi2: " << chi2_endcap_nume_100to500->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_nume_100to500->getVal() / ((Double_t)h_endcap_data_nume_100to500->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR BARREL 100to500
    cout << "\n----- DENOMINATOR BARREL 100to500 -----" << endl;
    TCanvas *c_fit_barrel_deno_100to500 = new TCanvas("c_fit_barrel_deno_100to500", "c_fit_barrel_deno_100to500", 800, 800);
    c_fit_barrel_deno_100to500->cd();

    //Top Pad
    TPad *c1_barrel_deno_100to500 = new TPad("padc1_barrel_deno_100to500","padc1_barrel_deno_100to500",0.01,0.01,0.99,0.99);
    c1_barrel_deno_100to500->Draw();
    c1_barrel_deno_100to500->cd();
    c1_barrel_deno_100to500->SetTopMargin(0.01);
    c1_barrel_deno_100to500->SetBottomMargin(0.35);
    c1_barrel_deno_100to500->SetRightMargin(0.03);
    c1_barrel_deno_100to500->SetLeftMargin(0.13);
    c1_barrel_deno_100to500->SetFillStyle(1);
    c1_barrel_deno_100to500->SetLogy();

    // Main stack histogram
    RooPlot *frame_barrel_deno_100to500 = iso_deno.frame(Title(" "));
    rh_barrel_deno_data_100to500->plotOn(frame_barrel_deno_100to500, DataError(RooAbsData::SumW2));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500,pdf_barrel_deno_tbarW_100to500,pdf_barrel_deno_ttbar_100to500,"
                                                                             "pdf_barrel_deno_DY_100to500,pdf_barrel_deno_WJets_100to500,pdf_barrel_deno_QCD_100to500"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500,pdf_barrel_deno_tbarW_100to500,pdf_barrel_deno_ttbar_100to500,"
                                                                             "pdf_barrel_deno_DY_100to500,pdf_barrel_deno_WJets_100to500"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500,pdf_barrel_deno_tbarW_100to500,pdf_barrel_deno_ttbar_100to500,"
                                                                             "pdf_barrel_deno_DY_100to500"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500,pdf_barrel_deno_tbarW_100to500,pdf_barrel_deno_ttbar_100to500"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500,pdf_barrel_deno_tbarW_100to500"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500,"
                                                                             "pdf_barrel_deno_tW_100to500"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500,pdf_barrel_deno_WW_100to500"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500,pdf_barrel_deno_WZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_deno_100to500.plotOn(frame_barrel_deno_100to500, Components("pdf_barrel_deno_ZZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_deno_data_100to500->plotOn(frame_barrel_deno_100to500, DataError(RooAbsData::SumW2));
    frame_barrel_deno_100to500->Draw();
    fit_barrel_deno_100to500->Print();

    // Legend
    legend->Draw();

    frame_barrel_deno_100to500->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_deno_100to500->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_deno_100to500->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_barrel_deno_100to500 = new TPad("padc2_barrel_deno_100to500","padc2_barrel_deno_100to500",0.01,0.01,0.99,0.35);
    c2_barrel_deno_100to500->Draw();
    c2_barrel_deno_100to500->cd();
    c2_barrel_deno_100to500->SetTopMargin(0.05);
    c2_barrel_deno_100to500->SetBottomMargin(0.33);
    c2_barrel_deno_100to500->SetRightMargin(0.02);
    c2_barrel_deno_100to500->SetLeftMargin(0.12);
    c2_barrel_deno_100to500->SetFillStyle(0);
    c2_barrel_deno_100to500->SetGrid();

    // Ratio plot
    TH1D *h_barrel_deno_MC_fit_100to500 = ((TH1D*)(model_barrel_deno_100to500.createHistogram("h_barrel_deno_MC_fit_100to500", iso_deno)));
    Double_t N_barrel_deno_data_100to500 = h_barrel_data_deno_100to500->Integral();
    Double_t N_barrel_deno_MC_100to500   = h_barrel_deno_MC_fit_100to500->Integral();
    h_barrel_deno_MC_fit_100to500->Scale(N_barrel_deno_data_100to500/N_barrel_deno_MC_100to500); // Why is this necessary???
    cout << "\nData integral: " << N_barrel_deno_data_100to500 << endl;
    cout << "MC integral: "     << h_barrel_deno_MC_fit_100to500->Integral() << endl;
    cout << "Data in 1st bin: " << h_barrel_data_deno_100to500->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_barrel_deno_MC_fit_100to500->GetBinContent(1) << endl;

    TH1D *h_barrel_deno_ratio_100to500 = ((TH1D*)(h_barrel_data_deno_100to500->Clone("h_barrel_deno_ratio_100to500")));
    h_barrel_data_deno_100to500->Sumw2(); h_barrel_deno_MC_fit_100to500->Sumw2();
    h_barrel_deno_ratio_100to500->Divide(h_barrel_data_deno_100to500, h_barrel_deno_MC_fit_100to500);
    h_barrel_deno_ratio_100to500->SetTitle("");
    h_barrel_deno_ratio_100to500->GetXaxis()->SetMoreLogLabels(1);
    h_barrel_deno_ratio_100to500->GetXaxis()->SetNoExponent(1);
    h_barrel_deno_ratio_100to500->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_barrel_deno_ratio_100to500->GetXaxis()->SetTitleSize(0.17);
    h_barrel_deno_ratio_100to500->GetXaxis()->SetLabelSize(0.125);
    h_barrel_deno_ratio_100to500->GetXaxis()->SetTitleOffset(0.8);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetTitle("Data/MC");
    h_barrel_deno_ratio_100to500->GetYaxis()->SetTitleSize(0.114);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetTitleOffset(0.48);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetLabelSize(0.11);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetTickLength(0.01);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetDecimals(1);
    h_barrel_deno_ratio_100to500->SetMaximum(1.25);
    h_barrel_deno_ratio_100to500->SetMinimum(0.75);
    h_barrel_deno_ratio_100to500->GetYaxis()->SetNdivisions(5);
    h_barrel_deno_ratio_100to500->SetLineWidth(1);
    h_barrel_deno_ratio_100to500->SetLineColor(kBlack);
    h_barrel_deno_ratio_100to500->SetMarkerStyle(kFullDotLarge);
    h_barrel_deno_ratio_100to500->SetMarkerColor(kBlack);
    h_barrel_deno_ratio_100to500->SetStats(kFALSE);

    h_barrel_deno_ratio_100to500->Draw("E1P");
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_deno_100to500 = model_barrel_deno_100to500.createChi2(*rh_barrel_deno_data_100to500);
    cout << "chi2: " << chi2_barrel_deno_100to500->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_deno_100to500->getVal() / ((Double_t)h_barrel_data_deno_100to500->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR ENDCAP 100to500
    cout << "\n----- DENOMINATOR ENDCAP 100to500 -----" << endl;
    TCanvas *c_fit_endcap_deno_100to500 = new TCanvas("c_fit_endcap_deno_100to500", "c_fit_endcap_deno_100to500", 800, 800);
    c_fit_endcap_deno_100to500->cd();

    //Top Pad
    TPad *c1_endcap_deno_100to500 = new TPad("padc1_endcap_deno_100to500","padc1_endcap_deno_100to500",0.01,0.01,0.99,0.99);
    c1_endcap_deno_100to500->Draw();
    c1_endcap_deno_100to500->cd();
    c1_endcap_deno_100to500->SetTopMargin(0.01);
    c1_endcap_deno_100to500->SetBottomMargin(0.35);
    c1_endcap_deno_100to500->SetRightMargin(0.03);
    c1_endcap_deno_100to500->SetLeftMargin(0.13);
    c1_endcap_deno_100to500->SetFillStyle(1);
    c1_endcap_deno_100to500->SetLogy();

    // Main stack histogram
    RooPlot *frame_endcap_deno_100to500 = iso_deno.frame(Title(" "));
    rh_endcap_deno_data_100to500->plotOn(frame_endcap_deno_100to500, DataError(RooAbsData::SumW2));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500,pdf_endcap_deno_tbarW_100to500,pdf_endcap_deno_ttbar_100to500,"
                                                                             "pdf_endcap_deno_DY_100to500,pdf_endcap_deno_WJets_100to500,pdf_endcap_deno_QCD_100to500"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500,pdf_endcap_deno_tbarW_100to500,pdf_endcap_deno_ttbar_100to500,"
                                                                             "pdf_endcap_deno_DY_100to500,pdf_endcap_deno_WJets_100to500"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500,pdf_endcap_deno_tbarW_100to500,pdf_endcap_deno_ttbar_100to500,"
                                                                             "pdf_endcap_deno_DY_100to500"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500,pdf_endcap_deno_tbarW_100to500,pdf_endcap_deno_ttbar_100to500"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500,pdf_endcap_deno_tbarW_100to500"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500,"
                                                                             "pdf_endcap_deno_tW_100to500"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500,pdf_endcap_deno_WW_100to500"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500,pdf_endcap_deno_WZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_deno_100to500.plotOn(frame_endcap_deno_100to500, Components("pdf_endcap_deno_ZZ_100to500"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_deno_data_100to500->plotOn(frame_endcap_deno_100to500, DataError(RooAbsData::SumW2));
    frame_endcap_deno_100to500->Draw();
    fit_endcap_deno_100to500->Print();

    // Legend
    legend->Draw();

    frame_endcap_deno_100to500->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_deno_100to500->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_deno_100to500->GetXaxis()->SetLabelSize(0);

    // Bottom pad
    TPad *c2_endcap_deno_100to500 = new TPad("padc2_endcap_deno_100to500","padc2_endcap_deno_100to500",0.01,0.01,0.99,0.35);
    c2_endcap_deno_100to500->Draw();
    c2_endcap_deno_100to500->cd();
    c2_endcap_deno_100to500->SetTopMargin(0.05);
    c2_endcap_deno_100to500->SetBottomMargin(0.33);
    c2_endcap_deno_100to500->SetRightMargin(0.02);
    c2_endcap_deno_100to500->SetLeftMargin(0.12);
    c2_endcap_deno_100to500->SetFillStyle(0);
    c2_endcap_deno_100to500->SetGrid();

    // Ratio plot
    TH1D *h_endcap_deno_MC_fit_100to500 = ((TH1D*)(model_endcap_deno_100to500.createHistogram("h_endcap_deno_MC_fit_100to500", iso_deno)));
    Double_t N_endcap_deno_data_100to500 = h_endcap_data_deno_100to500->Integral();
    Double_t N_endcap_deno_MC_100to500   = h_endcap_deno_MC_fit_100to500->Integral();
    h_endcap_deno_MC_fit_100to500->Scale(N_endcap_deno_data_100to500/N_endcap_deno_MC_100to500); // Why is this necessary???
    cout << "\nData integral: " << N_endcap_deno_data_100to500 << endl;
    cout << "MC integral: "     << h_endcap_deno_MC_fit_100to500->Integral() << endl;
    cout << "Data in 1st bin: " << h_endcap_data_deno_100to500->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_endcap_deno_MC_fit_100to500->GetBinContent(1) << endl;

    TH1D *h_endcap_deno_ratio_100to500 = ((TH1D*)(h_endcap_data_deno_100to500->Clone("h_endcap_deno_ratio_100to500")));
    h_endcap_data_deno_100to500->Sumw2(); h_endcap_deno_MC_fit_100to500->Sumw2();
    h_endcap_deno_ratio_100to500->Divide(h_endcap_data_deno_100to500, h_endcap_deno_MC_fit_100to500);
    h_endcap_deno_ratio_100to500->SetTitle("");
    h_endcap_deno_ratio_100to500->GetXaxis()->SetMoreLogLabels(1);
    h_endcap_deno_ratio_100to500->GetXaxis()->SetNoExponent(1);
    h_endcap_deno_ratio_100to500->GetXaxis()->SetTitle("I_{#lower[-0.2]{PF}}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})");
    h_endcap_deno_ratio_100to500->GetXaxis()->SetTitleSize(0.17);
    h_endcap_deno_ratio_100to500->GetXaxis()->SetLabelSize(0.125);
    h_endcap_deno_ratio_100to500->GetXaxis()->SetTitleOffset(0.8);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetTitle("Data/MC");
    h_endcap_deno_ratio_100to500->GetYaxis()->SetTitleSize(0.114);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetTitleOffset(0.48);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetLabelSize(0.11);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetTickLength(0.01);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetDecimals(1);
    h_endcap_deno_ratio_100to500->SetMaximum(1.25);
    h_endcap_deno_ratio_100to500->SetMinimum(0.75);
    h_endcap_deno_ratio_100to500->GetYaxis()->SetNdivisions(5);
    h_endcap_deno_ratio_100to500->SetLineWidth(1);
    h_endcap_deno_ratio_100to500->SetLineColor(kBlack);
    h_endcap_deno_ratio_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_deno_ratio_100to500->SetMarkerColor(kBlack);
    h_endcap_deno_ratio_100to500->SetStats(kFALSE);

    h_endcap_deno_ratio_100to500->Draw("E1P");

    // Red line at Data/MC=1
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_endcap_deno_100to500 = model_endcap_deno_100to500.createChi2(*rh_endcap_deno_data_100to500);
    cout << "chi2: " << chi2_endcap_deno_100to500->getVal() << endl;
    cout << "Normalized chi2: " << chi2_endcap_deno_100to500->getVal() / ((Double_t)h_endcap_data_deno_100to500->GetNbinsX()) << endl;


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
