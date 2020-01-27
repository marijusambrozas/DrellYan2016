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
void Mu_WJETSest_Tfit(Int_t type);
void Mu_WJETSest_Tfit_legacy();

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
        if (whichX.Contains("EST"))
        {
            if (whichX.Contains("W") && whichX.Contains("JET"))
            {
                cout << "\n*******     Mu_WJETSest_Tfit(" << type << ")     *******" << endl;
                Mu_WJETSest_Tfit(type);
            }
        }
        else
        {
            cout << "\n*******     Mu_Tfit(" << type << ")     *******" << endl;
            Mu_Tfit(type);
        }
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
         *h_barrel_MC_deno[_EndOf_Data_Special],
         *h_barrel_MC_nume[_EndOf_Data_Special],
         *h_endcap_MC_deno[_EndOf_Data_Special],
         *h_endcap_MC_nume[_EndOf_Data_Special],
         *h_barrel_data_deno_50to70,
         *h_barrel_data_nume_50to70,
         *h_endcap_data_deno_50to70,
         *h_endcap_data_nume_50to70,
         *h_barrel_data_deno_70to100,
         *h_barrel_data_nume_70to100,
         *h_endcap_data_deno_70to100,
         *h_endcap_data_nume_70to100,
         *h_barrel_data_deno,
         *h_barrel_data_nume,
         *h_endcap_data_deno,
         *h_endcap_data_nume;

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
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume[pr1]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_deno[pr1]);
        removeNegativeBins(h_endcap_MC_deno[pr1]);
        removeNegativeBins(h_barrel_MC_nume[pr1]);
        removeNegativeBins(h_endcap_MC_nume[pr1]);

        h_barrel_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr1]->SetDirectory(0);
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
    h_barrel_MC_deno_50to70  [_ttbar]->Add(h_barrel_MC_deno_50to70  [_ttbar_700to1000]);
    h_endcap_MC_deno_50to70  [_ttbar]->Add(h_endcap_MC_deno_50to70  [_ttbar_700to1000]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_700to1000]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_700to1000]);
    h_barrel_MC_deno_70to100 [_ttbar]->Add(h_barrel_MC_deno_70to100 [_ttbar_700to1000]);
    h_endcap_MC_deno_70to100 [_ttbar]->Add(h_endcap_MC_deno_70to100 [_ttbar_700to1000]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_700to1000]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_700to1000]);
    h_barrel_MC_deno[_ttbar]->Add(h_barrel_MC_deno[_ttbar_700to1000]);
    h_endcap_MC_deno[_ttbar]->Add(h_endcap_MC_deno[_ttbar_700to1000]);
    h_barrel_MC_nume[_ttbar]->Add(h_barrel_MC_nume[_ttbar_700to1000]);
    h_endcap_MC_nume[_ttbar]->Add(h_endcap_MC_nume[_ttbar_700to1000]);

    h_barrel_MC_deno_50to70  [_ttbar]->Add(h_barrel_MC_deno_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_deno_50to70  [_ttbar]->Add(h_endcap_MC_deno_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_deno_70to100 [_ttbar]->Add(h_barrel_MC_deno_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_deno_70to100 [_ttbar]->Add(h_endcap_MC_deno_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_deno[_ttbar]->Add(h_barrel_MC_deno[_ttbar_1000toInf]);
    h_endcap_MC_deno[_ttbar]->Add(h_endcap_MC_deno[_ttbar_1000toInf]);
    h_barrel_MC_nume[_ttbar]->Add(h_barrel_MC_nume[_ttbar_1000toInf]);
    h_endcap_MC_nume[_ttbar]->Add(h_endcap_MC_nume[_ttbar_1000toInf]);

    h_barrel_MC_deno_50to70  [_WJets]->Add(h_barrel_MC_deno_50to70  [_WJets_ext2v5]);
    h_endcap_MC_deno_50to70  [_WJets]->Add(h_endcap_MC_deno_50to70  [_WJets_ext2v5]);
    h_barrel_MC_nume_50to70  [_WJets]->Add(h_barrel_MC_nume_50to70  [_WJets_ext2v5]);
    h_endcap_MC_nume_50to70  [_WJets]->Add(h_endcap_MC_nume_50to70  [_WJets_ext2v5]);
    h_barrel_MC_deno_70to100 [_WJets]->Add(h_barrel_MC_deno_70to100 [_WJets_ext2v5]);
    h_endcap_MC_deno_70to100 [_WJets]->Add(h_endcap_MC_deno_70to100 [_WJets_ext2v5]);
    h_barrel_MC_nume_70to100 [_WJets]->Add(h_barrel_MC_nume_70to100 [_WJets_ext2v5]);
    h_endcap_MC_nume_70to100 [_WJets]->Add(h_endcap_MC_nume_70to100 [_WJets_ext2v5]);
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
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70  [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100 [pr]);
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno[pr]);
        removeNegativeBins(h_endcap_MC_deno[pr]);
        removeNegativeBins(h_barrel_MC_nume[pr]);
        removeNegativeBins(h_endcap_MC_nume[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_deno[pr]->SetDirectory(0);
        h_endcap_MC_deno[pr]->SetDirectory(0);
        h_barrel_MC_nume[pr]->SetDirectory(0);
        h_endcap_MC_nume[pr]->SetDirectory(0);

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
            h_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_barrel_MC_deno[pr]->Clone("h_barrel_MC_deno_DY")));
            h_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_endcap_MC_deno[pr]->Clone("h_endcap_MC_deno_DY")));
            h_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_barrel_MC_nume[pr]->Clone("h_barrel_MC_nume_DY")));
            h_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_endcap_MC_nume[pr]->Clone("h_endcap_MC_nume_DY")));

            h_barrel_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume[_DY_Full]->SetDirectory(0);
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
        file->GetObject("h_PFiso_barrel_deno_50to70",   h_barrel_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_deno_50to70",   h_endcap_MC_deno_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr]);
        file->GetObject("h_PFiso_barrel_deno_70to100",  h_barrel_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_endcap_deno_70to100",  h_endcap_MC_deno_70to100[pr]);
        file->GetObject("h_PFiso_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr]);
        file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_MC_nume[pr]);

        removeNegativeBins(h_barrel_MC_deno_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_deno_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_deno_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_deno_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_deno[pr]);
        removeNegativeBins(h_endcap_MC_deno[pr]);
        removeNegativeBins(h_barrel_MC_nume[pr]);
        removeNegativeBins(h_endcap_MC_nume[pr]);

        h_barrel_MC_deno_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_deno_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_deno_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_deno_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_deno[pr]->SetDirectory(0);
        h_endcap_MC_deno[pr]->SetDirectory(0);
        h_barrel_MC_nume[pr]->SetDirectory(0);
        h_endcap_MC_nume[pr]->SetDirectory(0);

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
            h_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_deno[pr]->Clone("h_barrel_MC_deno_QCD")));
            h_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_deno[pr]->Clone("h_endcap_MC_deno_QCD")));
            h_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume[pr]->Clone("h_barrel_MC_nume_QCD")));
            h_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume[pr]->Clone("h_endcap_MC_nume_QCD")));

            h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
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
            file->GetObject("h_PFiso_barrel_deno_100to500", h_barrel_data_deno);
            file->GetObject("h_PFiso_endcap_deno_100to500", h_endcap_data_deno);
            file->GetObject("h_PFiso_barrel_nume_100to500", h_barrel_data_nume);
            file->GetObject("h_PFiso_endcap_nume_100to500", h_endcap_data_nume);

            removeNegativeBins(h_barrel_data_deno_50to70);
            removeNegativeBins(h_endcap_data_deno_50to70);
            removeNegativeBins(h_barrel_data_nume_50to70);
            removeNegativeBins(h_endcap_data_nume_50to70);
            removeNegativeBins(h_barrel_data_deno_70to100);
            removeNegativeBins(h_endcap_data_deno_70to100);
            removeNegativeBins(h_barrel_data_nume_70to100);
            removeNegativeBins(h_endcap_data_nume_70to100);
            removeNegativeBins(h_barrel_data_deno);
            removeNegativeBins(h_endcap_data_deno);
            removeNegativeBins(h_barrel_data_nume);
            removeNegativeBins(h_endcap_data_nume);
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
            h_barrel_data_deno->Add(h_temp[8]);
            h_endcap_data_deno->Add(h_temp[9]);
            h_barrel_data_nume->Add(h_temp[10]);
            h_endcap_data_nume->Add(h_temp[11]);
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
    h_barrel_data_deno->SetDirectory(0);
    h_endcap_data_deno->SetDirectory(0);
    h_barrel_data_nume->SetDirectory(0);
    h_endcap_data_nume->SetDirectory(0);

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
    RooDataHist *rh_barrel_nume_QCD = new RooDataHist("rh_barrel_nume_QCD", "RooHist_barrel_nume_QCD", iso_nume, h_barrel_MC_nume[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_nume_QCD = new RooDataHist("rh_endcap_nume_QCD", "RooHist_endcap_nume_QCD", iso_nume, h_endcap_MC_nume[_QCDMuEnriched_Full]);
    RooDataHist *rh_barrel_deno_QCD = new RooDataHist("rh_barrel_deno_QCD", "RooHist_barrel_deno_QCD", iso_deno, h_barrel_MC_deno[_QCDMuEnriched_Full]);
    RooDataHist *rh_endcap_deno_QCD = new RooDataHist("rh_endcap_deno_QCD", "RooHist_endcap_deno_QCD", iso_deno, h_endcap_MC_deno[_QCDMuEnriched_Full]);

    RooDataHist *rh_barrel_nume_WJets_50to70   = new RooDataHist("rh_barrel_nume_WJets_50to70",   "RooHist_barrel_nume_WJets_50to70",   iso_nume, h_barrel_MC_nume_50to70 [_WJets]);
    RooDataHist *rh_endcap_nume_WJets_50to70   = new RooDataHist("rh_endcap_nume_WJets_50to70",   "RooHist_endcap_nume_WJets_50to70",   iso_nume, h_endcap_MC_nume_50to70 [_WJets]);
    RooDataHist *rh_barrel_deno_WJets_50to70   = new RooDataHist("rh_barrel_deno_WJets_50to70",   "RooHist_barrel_deno_WJets_50to70",   iso_deno, h_barrel_MC_deno_50to70 [_WJets]);
    RooDataHist *rh_endcap_deno_WJets_50to70   = new RooDataHist("rh_endcap_deno_WJets_50to70",   "RooHist_endcap_deno_WJets_50to70",   iso_deno, h_endcap_MC_deno_50to70 [_WJets]);
    RooDataHist *rh_barrel_nume_WJets_70to100  = new RooDataHist("rh_barrel_nume_WJets_70to100",  "RooHist_barrel_nume_WJets_70to100",  iso_nume, h_barrel_MC_nume_70to100[_WJets]);
    RooDataHist *rh_endcap_nume_WJets_70to100  = new RooDataHist("rh_endcap_nume_WJets_70to100",  "RooHist_endcap_nume_WJets_70to100",  iso_nume, h_endcap_MC_nume_70to100[_WJets]);
    RooDataHist *rh_barrel_deno_WJets_70to100  = new RooDataHist("rh_barrel_deno_WJets_70to100",  "RooHist_barrel_deno_WJets_70to100",  iso_deno, h_barrel_MC_deno_70to100[_WJets]);
    RooDataHist *rh_endcap_deno_WJets_70to100  = new RooDataHist("rh_endcap_deno_WJets_70to100",  "RooHist_endcap_deno_WJets_70to100",  iso_deno, h_endcap_MC_deno_70to100[_WJets]);
    RooDataHist *rh_barrel_nume_WJets = new RooDataHist("rh_barrel_nume_WJets", "RooHist_barrel_nume_WJets", iso_nume, h_barrel_MC_nume[_WJets]);
    RooDataHist *rh_endcap_nume_WJets = new RooDataHist("rh_endcap_nume_WJets", "RooHist_endcap_nume_WJets", iso_nume, h_endcap_MC_nume[_WJets]);
    RooDataHist *rh_barrel_deno_WJets = new RooDataHist("rh_barrel_deno_WJets", "RooHist_barrel_deno_WJets", iso_deno, h_barrel_MC_deno[_WJets]);
    RooDataHist *rh_endcap_deno_WJets = new RooDataHist("rh_endcap_deno_WJets", "RooHist_endcap_deno_WJets", iso_deno, h_endcap_MC_deno[_WJets]);

    RooDataHist *rh_barrel_nume_DY_50to70   = new RooDataHist("rh_barrel_nume_DY_50to70",   "RooHist_barrel_nume_DY_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_DY_Full]);
    RooDataHist *rh_endcap_nume_DY_50to70   = new RooDataHist("rh_endcap_nume_DY_50to70",   "RooHist_endcap_nume_DY_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_DY_Full]);
    RooDataHist *rh_barrel_deno_DY_50to70   = new RooDataHist("rh_barrel_deno_DY_50to70",   "RooHist_barrel_deno_DY_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_DY_Full]);
    RooDataHist *rh_endcap_deno_DY_50to70   = new RooDataHist("rh_endcap_deno_DY_50to70",   "RooHist_endcap_deno_DY_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_DY_Full]);
    RooDataHist *rh_barrel_nume_DY_70to100  = new RooDataHist("rh_barrel_nume_DY_70to100",  "RooHist_barrel_nume_DY_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_DY_Full]);
    RooDataHist *rh_endcap_nume_DY_70to100  = new RooDataHist("rh_endcap_nume_DY_70to100",  "RooHist_endcap_nume_DY_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_DY_Full]);
    RooDataHist *rh_barrel_deno_DY_70to100  = new RooDataHist("rh_barrel_deno_DY_70to100",  "RooHist_barrel_deno_DY_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_DY_Full]);
    RooDataHist *rh_endcap_deno_DY_70to100  = new RooDataHist("rh_endcap_deno_DY_70to100",  "RooHist_endcap_deno_DY_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_DY_Full]);
    RooDataHist *rh_barrel_nume_DY = new RooDataHist("rh_barrel_nume_DY", "RooHist_barrel_nume_DY", iso_nume, h_barrel_MC_nume[_DY_Full]);
    RooDataHist *rh_endcap_nume_DY = new RooDataHist("rh_endcap_nume_DY", "RooHist_endcap_nume_DY", iso_nume, h_endcap_MC_nume[_DY_Full]);
    RooDataHist *rh_barrel_deno_DY = new RooDataHist("rh_barrel_deno_DY", "RooHist_barrel_deno_DY", iso_deno, h_barrel_MC_deno[_DY_Full]);
    RooDataHist *rh_endcap_deno_DY = new RooDataHist("rh_endcap_deno_DY", "RooHist_endcap_deno_DY", iso_deno, h_endcap_MC_deno[_DY_Full]);

    RooDataHist *rh_barrel_nume_ttbar_50to70   = new RooDataHist("rh_barrel_nume_ttbar_50to70",   "RooHist_barrel_nume_ttbar_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar_50to70   = new RooDataHist("rh_endcap_nume_ttbar_50to70",   "RooHist_endcap_nume_ttbar_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar_50to70   = new RooDataHist("rh_barrel_deno_ttbar_50to70",   "RooHist_barrel_deno_ttbar_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar_50to70   = new RooDataHist("rh_endcap_deno_ttbar_50to70",   "RooHist_endcap_deno_ttbar_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_ttbar]);
    RooDataHist *rh_barrel_nume_ttbar_70to100  = new RooDataHist("rh_barrel_nume_ttbar_70to100",  "RooHist_barrel_nume_ttbar_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar_70to100  = new RooDataHist("rh_endcap_nume_ttbar_70to100",  "RooHist_endcap_nume_ttbar_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar_70to100  = new RooDataHist("rh_barrel_deno_ttbar_70to100",  "RooHist_barrel_deno_ttbar_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar_70to100  = new RooDataHist("rh_endcap_deno_ttbar_70to100",  "RooHist_endcap_deno_ttbar_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_ttbar]);
    RooDataHist *rh_barrel_nume_ttbar = new RooDataHist("rh_barrel_nume_ttbar", "RooHist_barrel_nume_ttbar", iso_nume, h_barrel_MC_nume[_ttbar]);
    RooDataHist *rh_endcap_nume_ttbar = new RooDataHist("rh_endcap_nume_ttbar", "RooHist_endcap_nume_ttbar", iso_nume, h_endcap_MC_nume[_ttbar]);
    RooDataHist *rh_barrel_deno_ttbar = new RooDataHist("rh_barrel_deno_ttbar", "RooHist_barrel_deno_ttbar", iso_deno, h_barrel_MC_deno[_ttbar]);
    RooDataHist *rh_endcap_deno_ttbar = new RooDataHist("rh_endcap_deno_ttbar", "RooHist_endcap_deno_ttbar", iso_deno, h_endcap_MC_deno[_ttbar]);

    RooDataHist *rh_barrel_nume_tW_50to70   = new RooDataHist("rh_barrel_nume_tW_50to70",   "RooHist_barrel_nume_tW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_tW]);
    RooDataHist *rh_endcap_nume_tW_50to70   = new RooDataHist("rh_endcap_nume_tW_50to70",   "RooHist_endcap_nume_tW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_tW]);
    RooDataHist *rh_barrel_deno_tW_50to70   = new RooDataHist("rh_barrel_deno_tW_50to70",   "RooHist_barrel_deno_tW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_tW]);
    RooDataHist *rh_endcap_deno_tW_50to70   = new RooDataHist("rh_endcap_deno_tW_50to70",   "RooHist_endcap_deno_tW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_tW]);
    RooDataHist *rh_barrel_nume_tW_70to100  = new RooDataHist("rh_barrel_nume_tW_70to100",  "RooHist_barrel_nume_tW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_tW]);
    RooDataHist *rh_endcap_nume_tW_70to100  = new RooDataHist("rh_endcap_nume_tW_70to100",  "RooHist_endcap_nume_tW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_tW]);
    RooDataHist *rh_barrel_deno_tW_70to100  = new RooDataHist("rh_barrel_deno_tW_70to100",  "RooHist_barrel_deno_tW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_tW]);
    RooDataHist *rh_endcap_deno_tW_70to100  = new RooDataHist("rh_endcap_deno_tW_70to100",  "RooHist_endcap_deno_tW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_tW]);
    RooDataHist *rh_barrel_nume_tW = new RooDataHist("rh_barrel_nume_tW", "RooHist_barrel_nume_tW", iso_nume, h_barrel_MC_nume[_tW]);
    RooDataHist *rh_endcap_nume_tW = new RooDataHist("rh_endcap_nume_tW", "RooHist_endcap_nume_tW", iso_nume, h_endcap_MC_nume[_tW]);
    RooDataHist *rh_barrel_deno_tW = new RooDataHist("rh_barrel_deno_tW", "RooHist_barrel_deno_tW", iso_deno, h_barrel_MC_deno[_tW]);
    RooDataHist *rh_endcap_deno_tW = new RooDataHist("rh_endcap_deno_tW", "RooHist_endcap_deno_tW", iso_deno, h_endcap_MC_deno[_tW]);

    RooDataHist *rh_barrel_nume_tbarW_50to70   = new RooDataHist("rh_barrel_nume_tbarW_50to70",   "RooHist_barrel_nume_tbarW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW_50to70   = new RooDataHist("rh_endcap_nume_tbarW_50to70",   "RooHist_endcap_nume_tbarW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW_50to70   = new RooDataHist("rh_barrel_deno_tbarW_50to70",   "RooHist_barrel_deno_tbarW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW_50to70   = new RooDataHist("rh_endcap_deno_tbarW_50to70",   "RooHist_endcap_deno_tbarW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_tbarW]);
    RooDataHist *rh_barrel_nume_tbarW_70to100  = new RooDataHist("rh_barrel_nume_tbarW_70to100",  "RooHist_barrel_nume_tbarW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW_70to100  = new RooDataHist("rh_endcap_nume_tbarW_70to100",  "RooHist_endcap_nume_tbarW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW_70to100  = new RooDataHist("rh_barrel_deno_tbarW_70to100",  "RooHist_barrel_deno_tbarW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW_70to100  = new RooDataHist("rh_endcap_deno_tbarW_70to100",  "RooHist_endcap_deno_tbarW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_tbarW]);
    RooDataHist *rh_barrel_nume_tbarW = new RooDataHist("rh_barrel_nume_tbarW", "RooHist_barrel_nume_tbarW", iso_nume, h_barrel_MC_nume[_tbarW]);
    RooDataHist *rh_endcap_nume_tbarW = new RooDataHist("rh_endcap_nume_tbarW", "RooHist_endcap_nume_tbarW", iso_nume, h_endcap_MC_nume[_tbarW]);
    RooDataHist *rh_barrel_deno_tbarW = new RooDataHist("rh_barrel_deno_tbarW", "RooHist_barrel_deno_tbarW", iso_deno, h_barrel_MC_deno[_tbarW]);
    RooDataHist *rh_endcap_deno_tbarW = new RooDataHist("rh_endcap_deno_tbarW", "RooHist_endcap_deno_tbarW", iso_deno, h_endcap_MC_deno[_tbarW]);

    RooDataHist *rh_barrel_nume_WW_50to70  =  new RooDataHist("rh_barrel_nume_WW_50to70",   "RooHist_barrel_nume_WW_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_WW]);
    RooDataHist *rh_endcap_nume_WW_50to70  =  new RooDataHist("rh_endcap_nume_WW_50to70",   "RooHist_endcap_nume_WW_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_WW]);
    RooDataHist *rh_barrel_deno_WW_50to70  =  new RooDataHist("rh_barrel_deno_WW_50to70",   "RooHist_barrel_deno_WW_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_WW]);
    RooDataHist *rh_endcap_deno_WW_50to70  =  new RooDataHist("rh_endcap_deno_WW_50to70",   "RooHist_endcap_deno_WW_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_WW]);
    RooDataHist *rh_barrel_nume_WW_70to100 =  new RooDataHist("rh_barrel_nume_WW_70to100",  "RooHist_barrel_nume_WW_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_WW]);
    RooDataHist *rh_endcap_nume_WW_70to100 =  new RooDataHist("rh_endcap_nume_WW_70to100",  "RooHist_endcap_nume_WW_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_WW]);
    RooDataHist *rh_barrel_deno_WW_70to100 =  new RooDataHist("rh_barrel_deno_WW_70to100",  "RooHist_barrel_deno_WW_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_WW]);
    RooDataHist *rh_endcap_deno_WW_70to100 =  new RooDataHist("rh_endcap_deno_WW_70to100",  "RooHist_endcap_deno_WW_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_WW]);
    RooDataHist *rh_barrel_nume_WW = new RooDataHist("rh_barrel_nume_WW", "RooHist_barrel_nume_WW", iso_nume, h_barrel_MC_nume[_WW]);
    RooDataHist *rh_endcap_nume_WW = new RooDataHist("rh_endcap_nume_WW", "RooHist_endcap_nume_WW", iso_nume, h_endcap_MC_nume[_WW]);
    RooDataHist *rh_barrel_deno_WW = new RooDataHist("rh_barrel_deno_WW", "RooHist_barrel_deno_WW", iso_deno, h_barrel_MC_deno[_WW]);
    RooDataHist *rh_endcap_deno_WW = new RooDataHist("rh_endcap_deno_WW", "RooHist_endcap_deno_WW", iso_deno, h_endcap_MC_deno[_WW]);

    RooDataHist *rh_barrel_nume_WZ_50to70   = new RooDataHist("rh_barrel_nume_WZ_50to70",   "RooHist_barrel_nume_WZ_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_WZ]);
    RooDataHist *rh_endcap_nume_WZ_50to70   = new RooDataHist("rh_endcap_nume_WZ_50to70",   "RooHist_endcap_nume_WZ_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_WZ]);
    RooDataHist *rh_barrel_deno_WZ_50to70   = new RooDataHist("rh_barrel_deno_WZ_50to70",   "RooHist_barrel_deno_WZ_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_WZ]);
    RooDataHist *rh_endcap_deno_WZ_50to70   = new RooDataHist("rh_endcap_deno_WZ_50to70",   "RooHist_endcap_deno_WZ_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_WZ]);
    RooDataHist *rh_barrel_nume_WZ_70to100  = new RooDataHist("rh_barrel_nume_WZ_70to100",  "RooHist_barrel_nume_WZ_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_WZ]);
    RooDataHist *rh_endcap_nume_WZ_70to100  = new RooDataHist("rh_endcap_nume_WZ_70to100",  "RooHist_endcap_nume_WZ_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_WZ]);
    RooDataHist *rh_barrel_deno_WZ_70to100  = new RooDataHist("rh_barrel_deno_WZ_70to100",  "RooHist_barrel_deno_WZ_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_WZ]);
    RooDataHist *rh_endcap_deno_WZ_70to100  = new RooDataHist("rh_endcap_deno_WZ_70to100",  "RooHist_endcap_deno_WZ_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_WZ]);
    RooDataHist *rh_barrel_nume_WZ = new RooDataHist("rh_barrel_nume_WZ", "RooHist_barrel_nume_WZ", iso_nume, h_barrel_MC_nume[_WZ]);
    RooDataHist *rh_endcap_nume_WZ = new RooDataHist("rh_endcap_nume_WZ", "RooHist_endcap_nume_WZ", iso_nume, h_endcap_MC_nume[_WZ]);
    RooDataHist *rh_barrel_deno_WZ = new RooDataHist("rh_barrel_deno_WZ", "RooHist_barrel_deno_WZ", iso_deno, h_barrel_MC_deno[_WZ]);
    RooDataHist *rh_endcap_deno_WZ = new RooDataHist("rh_endcap_deno_WZ", "RooHist_endcap_deno_WZ", iso_deno, h_endcap_MC_deno[_WZ]);

    RooDataHist *rh_barrel_nume_ZZ_50to70   = new RooDataHist("rh_barrel_nume_ZZ_50to70",   "RooHist_barrel_nume_ZZ_50to70",   iso_nume, h_barrel_MC_nume_50to70  [_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ_50to70   = new RooDataHist("rh_endcap_nume_ZZ_50to70",   "RooHist_endcap_nume_ZZ_50to70",   iso_nume, h_endcap_MC_nume_50to70  [_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ_50to70   = new RooDataHist("rh_barrel_deno_ZZ_50to70",   "RooHist_barrel_deno_ZZ_50to70",   iso_deno, h_barrel_MC_deno_50to70  [_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ_50to70   = new RooDataHist("rh_endcap_deno_ZZ_50to70",   "RooHist_endcap_deno_ZZ_50to70",   iso_deno, h_endcap_MC_deno_50to70  [_ZZ]);
    RooDataHist *rh_barrel_nume_ZZ_70to100  = new RooDataHist("rh_barrel_nume_ZZ_70to100",  "RooHist_barrel_nume_ZZ_70to100",  iso_nume, h_barrel_MC_nume_70to100 [_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ_70to100  = new RooDataHist("rh_endcap_nume_ZZ_70to100",  "RooHist_endcap_nume_ZZ_70to100",  iso_nume, h_endcap_MC_nume_70to100 [_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ_70to100  = new RooDataHist("rh_barrel_deno_ZZ_70to100",  "RooHist_barrel_deno_ZZ_70to100",  iso_deno, h_barrel_MC_deno_70to100 [_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ_70to100  = new RooDataHist("rh_endcap_deno_ZZ_70to100",  "RooHist_endcap_deno_ZZ_70to100",  iso_deno, h_endcap_MC_deno_70to100 [_ZZ]);
    RooDataHist *rh_barrel_nume_ZZ = new RooDataHist("rh_barrel_nume_ZZ", "RooHist_barrel_nume_ZZ", iso_nume, h_barrel_MC_nume[_ZZ]);
    RooDataHist *rh_endcap_nume_ZZ = new RooDataHist("rh_endcap_nume_ZZ", "RooHist_endcap_nume_ZZ", iso_nume, h_endcap_MC_nume[_ZZ]);
    RooDataHist *rh_barrel_deno_ZZ = new RooDataHist("rh_barrel_deno_ZZ", "RooHist_barrel_deno_ZZ", iso_deno, h_barrel_MC_deno[_ZZ]);
    RooDataHist *rh_endcap_deno_ZZ = new RooDataHist("rh_endcap_deno_ZZ", "RooHist_endcap_deno_ZZ", iso_deno, h_endcap_MC_deno[_ZZ]);

    RooDataHist *rh_barrel_nume_data_50to70   = new RooDataHist("rh_barrel_nume_data_50to70",   "RooHist_barrel_nume_data_50to70",   iso_nume, h_barrel_data_nume_50to70 );
    RooDataHist *rh_endcap_nume_data_50to70   = new RooDataHist("rh_endcap_nume_data_50to70",   "RooHist_endcap_nume_data_50to70",   iso_nume, h_endcap_data_nume_50to70 );
    RooDataHist *rh_barrel_deno_data_50to70   = new RooDataHist("rh_barrel_deno_data_50to70",   "RooHist_barrel_deno_data_50to70",   iso_deno, h_barrel_data_deno_50to70 );
    RooDataHist *rh_endcap_deno_data_50to70   = new RooDataHist("rh_endcap_deno_data_50to70",   "RooHist_endcap_deno_data_50to70",   iso_deno, h_endcap_data_deno_50to70 );
    RooDataHist *rh_barrel_nume_data_70to100  = new RooDataHist("rh_barrel_nume_data_70to100",  "RooHist_barrel_nume_data_70to100",  iso_nume, h_barrel_data_nume_70to100);
    RooDataHist *rh_endcap_nume_data_70to100  = new RooDataHist("rh_endcap_nume_data_70to100",  "RooHist_endcap_nume_data_70to100",  iso_nume, h_endcap_data_nume_70to100);
    RooDataHist *rh_barrel_deno_data_70to100  = new RooDataHist("rh_barrel_deno_data_70to100",  "RooHist_barrel_deno_data_70to100",  iso_deno, h_barrel_data_deno_70to100);
    RooDataHist *rh_endcap_deno_data_70to100  = new RooDataHist("rh_endcap_deno_data_70to100",  "RooHist_endcap_deno_data_70to100",  iso_deno, h_endcap_data_deno_70to100);
    RooDataHist *rh_barrel_nume_data = new RooDataHist("rh_barrel_nume_data", "RooHist_barrel_nume_data", iso_nume, h_barrel_data_nume);
    RooDataHist *rh_endcap_nume_data = new RooDataHist("rh_endcap_nume_data", "RooHist_endcap_nume_data", iso_nume, h_endcap_data_nume);
    RooDataHist *rh_barrel_deno_data = new RooDataHist("rh_barrel_deno_data", "RooHist_barrel_deno_data", iso_deno, h_barrel_data_deno);
    RooDataHist *rh_endcap_deno_data = new RooDataHist("rh_endcap_deno_data", "RooHist_endcap_deno_data", iso_deno, h_endcap_data_deno);

    // Making RooHistPdf
    RooHistPdf *pdf_barrel_nume_QCD_50to70   = new RooHistPdf("pdf_barrel_nume_QCD_50to70",   "Numerator barrel MC QCD template 50 to 70",     iso_nume, *rh_barrel_nume_QCD_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_QCD_50to70   = new RooHistPdf("pdf_endcap_nume_QCD_50to70",   "Numerator endcap MC QCD template 50 to 70",     iso_nume, *rh_endcap_nume_QCD_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_QCD_50to70   = new RooHistPdf("pdf_barrel_deno_QCD_50to70",   "Denominator barrel MC QCD template 50 to 70",   iso_deno, *rh_barrel_deno_QCD_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_QCD_50to70   = new RooHistPdf("pdf_endcap_deno_QCD_50to70",   "Denominator endcap MC QCD template 50 to 70",   iso_deno, *rh_endcap_deno_QCD_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_QCD_70to100  = new RooHistPdf("pdf_barrel_nume_QCD_70to100",  "Numerator barrel MC QCD template 70 to 100",    iso_nume, *rh_barrel_nume_QCD_70to100,  0);
    RooHistPdf *pdf_endcap_nume_QCD_70to100  = new RooHistPdf("pdf_endcap_nume_QCD_70to100",  "Numerator endcap MC QCD template 70 to 100",    iso_nume, *rh_endcap_nume_QCD_70to100,  0);
    RooHistPdf *pdf_barrel_deno_QCD_70to100  = new RooHistPdf("pdf_barrel_deno_QCD_70to100",  "Denominator barrel MC QCD template 70 to 100",  iso_deno, *rh_barrel_deno_QCD_70to100,  0);
    RooHistPdf *pdf_endcap_deno_QCD_70to100  = new RooHistPdf("pdf_endcap_deno_QCD_70to100",  "Denominator endcap MC QCD template 70 to 100",  iso_deno, *rh_endcap_deno_QCD_70to100,  0);
    RooHistPdf *pdf_barrel_nume_QCD = new RooHistPdf("pdf_barrel_nume_QCD", "Numerator barrel MC QCD template 200 to 500",   iso_nume, *rh_barrel_nume_QCD, 0);
    RooHistPdf *pdf_endcap_nume_QCD = new RooHistPdf("pdf_endcap_nume_QCD", "Numerator endcap MC QCD template 200 to 500",   iso_nume, *rh_endcap_nume_QCD, 0);
    RooHistPdf *pdf_barrel_deno_QCD = new RooHistPdf("pdf_barrel_deno_QCD", "Denominator barrel MC QCD template 200 to 500", iso_deno, *rh_barrel_deno_QCD, 0);
    RooHistPdf *pdf_endcap_deno_QCD = new RooHistPdf("pdf_endcap_deno_QCD", "Denominator endcap MC QCD template 200 to 500", iso_deno, *rh_endcap_deno_QCD, 0);

    RooHistPdf *pdf_barrel_nume_WJets_50to70   = new RooHistPdf("pdf_barrel_nume_WJets_50to70",   "Numerator barrel MC W+Jets template 50 to 70",     iso_nume, *rh_barrel_nume_WJets_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WJets_50to70   = new RooHistPdf("pdf_endcap_nume_WJets_50to70",   "Numerator endcap MC W+Jets template 50 to 70",     iso_nume, *rh_endcap_nume_WJets_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WJets_50to70   = new RooHistPdf("pdf_barrel_deno_WJets_50to70",   "Denominator barrel MC W+Jets template 50 to 70",   iso_deno, *rh_barrel_deno_WJets_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WJets_50to70   = new RooHistPdf("pdf_endcap_deno_WJets_50to70",   "Denominator endcap MC W+Jets template 50 to 70",   iso_deno, *rh_endcap_deno_WJets_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WJets_70to100  = new RooHistPdf("pdf_barrel_nume_WJets_70to100",  "Numerator barrel MC W+Jets template 70 to 100",    iso_nume, *rh_barrel_nume_WJets_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WJets_70to100  = new RooHistPdf("pdf_endcap_nume_WJets_70to100",  "Numerator endcap MC W+Jets template 70 to 100",    iso_nume, *rh_endcap_nume_WJets_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WJets_70to100  = new RooHistPdf("pdf_barrel_deno_WJets_70to100",  "Denominator barrel MC W+Jets template 70 to 100",  iso_deno, *rh_barrel_deno_WJets_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WJets_70to100  = new RooHistPdf("pdf_endcap_deno_WJets_70to100",  "Denominator endcap MC W+Jets template 70 to 100",  iso_deno, *rh_endcap_deno_WJets_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WJets = new RooHistPdf("pdf_barrel_nume_WJets", "Numerator barrel MC W+Jets template 200 to 500",   iso_nume, *rh_barrel_nume_WJets, 0);
    RooHistPdf *pdf_endcap_nume_WJets = new RooHistPdf("pdf_endcap_nume_WJets", "Numerator endcap MC W+Jets template 200 to 500",   iso_nume, *rh_endcap_nume_WJets, 0);
    RooHistPdf *pdf_barrel_deno_WJets = new RooHistPdf("pdf_barrel_deno_WJets", "Denominator barrel MC W+Jets template 200 to 500", iso_deno, *rh_barrel_deno_WJets, 0);
    RooHistPdf *pdf_endcap_deno_WJets = new RooHistPdf("pdf_endcap_deno_WJets", "Denominator endcap MC W+Jets template 200 to 500", iso_deno, *rh_endcap_deno_WJets, 0);

    RooHistPdf *pdf_barrel_nume_DY_50to70   = new RooHistPdf("pdf_barrel_nume_DY_50to70",   "Numerator barrel MC DY template 50 to 70",     iso_nume, *rh_barrel_nume_DY_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_DY_50to70   = new RooHistPdf("pdf_endcap_nume_DY_50to70",   "Numerator endcap MC DY template 50 to 70",     iso_nume, *rh_endcap_nume_DY_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_DY_50to70   = new RooHistPdf("pdf_barrel_deno_DY_50to70",   "Denominator barrel MC DY template 50 to 70",   iso_deno, *rh_barrel_deno_DY_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_DY_50to70   = new RooHistPdf("pdf_endcap_deno_DY_50to70",   "Denominator endcap MC DY template 50 to 70",   iso_deno, *rh_endcap_deno_DY_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_DY_70to100  = new RooHistPdf("pdf_barrel_nume_DY_70to100",  "Numerator barrel MC DY template 70 to 100",    iso_nume, *rh_barrel_nume_DY_70to100,  0);
    RooHistPdf *pdf_endcap_nume_DY_70to100  = new RooHistPdf("pdf_endcap_nume_DY_70to100",  "Numerator endcap MC DY template 70 to 100",    iso_nume, *rh_endcap_nume_DY_70to100,  0);
    RooHistPdf *pdf_barrel_deno_DY_70to100  = new RooHistPdf("pdf_barrel_deno_DY_70to100",  "Denominator barrel MC DY template 70 to 100",  iso_deno, *rh_barrel_deno_DY_70to100,  0);
    RooHistPdf *pdf_endcap_deno_DY_70to100  = new RooHistPdf("pdf_endcap_deno_DY_70to100",  "Denominator endcap MC DY template 70 to 100",  iso_deno, *rh_endcap_deno_DY_70to100,  0);
    RooHistPdf *pdf_barrel_nume_DY = new RooHistPdf("pdf_barrel_nume_DY", "Numerator barrel MC DY template 200 to 500",   iso_nume, *rh_barrel_nume_DY, 0);
    RooHistPdf *pdf_endcap_nume_DY = new RooHistPdf("pdf_endcap_nume_DY", "Numerator endcap MC DY template 200 to 500",   iso_nume, *rh_endcap_nume_DY, 0);
    RooHistPdf *pdf_barrel_deno_DY = new RooHistPdf("pdf_barrel_deno_DY", "Denominator barrel MC DY template 200 to 500", iso_deno, *rh_barrel_deno_DY, 0);
    RooHistPdf *pdf_endcap_deno_DY = new RooHistPdf("pdf_endcap_deno_DY", "Denominator endcap MC DY template 200 to 500", iso_deno, *rh_endcap_deno_DY, 0);

    RooHistPdf *pdf_barrel_nume_ttbar_50to70   = new RooHistPdf("pdf_barrel_nume_ttbar_50to70",   "Numerator barrel MC ttbar template 50 to 70",     iso_nume, *rh_barrel_nume_ttbar_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_ttbar_50to70   = new RooHistPdf("pdf_endcap_nume_ttbar_50to70",   "Numerator endcap MC ttbar template 50 to 70",     iso_nume, *rh_endcap_nume_ttbar_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_ttbar_50to70   = new RooHistPdf("pdf_barrel_deno_ttbar_50to70",   "Denominator barrel MC ttbar template 50 to 70",   iso_deno, *rh_barrel_deno_ttbar_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_ttbar_50to70   = new RooHistPdf("pdf_endcap_deno_ttbar_50to70",   "Denominator endcap MC ttbar template 50 to 70",   iso_deno, *rh_endcap_deno_ttbar_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_ttbar_70to100  = new RooHistPdf("pdf_barrel_nume_ttbar_70to100",  "Numerator barrel MC ttbar template 70 to 100",    iso_nume, *rh_barrel_nume_ttbar_70to100,  0);
    RooHistPdf *pdf_endcap_nume_ttbar_70to100  = new RooHistPdf("pdf_endcap_nume_ttbar_70to100",  "Numerator endcap MC ttbar template 70 to 100",    iso_nume, *rh_endcap_nume_ttbar_70to100,  0);
    RooHistPdf *pdf_barrel_deno_ttbar_70to100  = new RooHistPdf("pdf_barrel_deno_ttbar_70to100",  "Denominator barrel MC ttbar template 70 to 100",  iso_deno, *rh_barrel_deno_ttbar_70to100,  0);
    RooHistPdf *pdf_endcap_deno_ttbar_70to100  = new RooHistPdf("pdf_endcap_deno_ttbar_70to100",  "Denominator endcap MC ttbar template 70 to 100",  iso_deno, *rh_endcap_deno_ttbar_70to100,  0);
    RooHistPdf *pdf_barrel_nume_ttbar = new RooHistPdf("pdf_barrel_nume_ttbar", "Numerator barrel MC ttbar template 200 to 500",   iso_nume, *rh_barrel_nume_ttbar, 0);
    RooHistPdf *pdf_endcap_nume_ttbar = new RooHistPdf("pdf_endcap_nume_ttbar", "Numerator endcap MC ttbar template 200 to 500",   iso_nume, *rh_endcap_nume_ttbar, 0);
    RooHistPdf *pdf_barrel_deno_ttbar = new RooHistPdf("pdf_barrel_deno_ttbar", "Denominator barrel MC ttbar template 200 to 500", iso_deno, *rh_barrel_deno_ttbar, 0);
    RooHistPdf *pdf_endcap_deno_ttbar = new RooHistPdf("pdf_endcap_deno_ttbar", "Denominator endcap MC ttbar template 200 to 500", iso_deno, *rh_endcap_deno_ttbar, 0);

    RooHistPdf *pdf_barrel_nume_tW_50to70   = new RooHistPdf("pdf_barrel_nume_tW_50to70",   "Numerator barrel MC tW template 50 to 70",     iso_nume, *rh_barrel_nume_tW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_tW_50to70   = new RooHistPdf("pdf_endcap_nume_tW_50to70",   "Numerator endcap MC tW template 50 to 70",     iso_nume, *rh_endcap_nume_tW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_tW_50to70   = new RooHistPdf("pdf_barrel_deno_tW_50to70",   "Denominator barrel MC tW template 50 to 70",   iso_deno, *rh_barrel_deno_tW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_tW_50to70   = new RooHistPdf("pdf_endcap_deno_tW_50to70",   "Denominator endcap MC tW template 50 to 70",   iso_deno, *rh_endcap_deno_tW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_tW_70to100  = new RooHistPdf("pdf_barrel_nume_tW_70to100",  "Numerator barrel MC tW template 70 to 100",    iso_nume, *rh_barrel_nume_tW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_tW_70to100  = new RooHistPdf("pdf_endcap_nume_tW_70to100",  "Numerator endcap MC tW template 70 to 100",    iso_nume, *rh_endcap_nume_tW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_tW_70to100  = new RooHistPdf("pdf_barrel_deno_tW_70to100",  "Denominator barrel MC tW template 70 to 100",  iso_deno, *rh_barrel_deno_tW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_tW_70to100  = new RooHistPdf("pdf_endcap_deno_tW_70to100",  "Denominator endcap MC tW template 70 to 100",  iso_deno, *rh_endcap_deno_tW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_tW = new RooHistPdf("pdf_barrel_nume_tW", "Numerator barrel MC tW template 200 to 500",   iso_nume, *rh_barrel_nume_tW, 0);
    RooHistPdf *pdf_endcap_nume_tW = new RooHistPdf("pdf_endcap_nume_tW", "Numerator endcap MC tW template 200 to 500",   iso_nume, *rh_endcap_nume_tW, 0);
    RooHistPdf *pdf_barrel_deno_tW = new RooHistPdf("pdf_barrel_deno_tW", "Denominator barrel MC tW template 200 to 500", iso_deno, *rh_barrel_deno_tW, 0);
    RooHistPdf *pdf_endcap_deno_tW = new RooHistPdf("pdf_endcap_deno_tW", "Denominator endcap MC tW template 200 to 500", iso_deno, *rh_endcap_deno_tW, 0);

    RooHistPdf *pdf_barrel_nume_tbarW_50to70   = new RooHistPdf("pdf_barrel_nume_tbarW_50to70",   "Numerator barrel MC tbarW template 50 to 70",     iso_nume, *rh_barrel_nume_tbarW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_tbarW_50to70   = new RooHistPdf("pdf_endcap_nume_tbarW_50to70",   "Numerator endcap MC tbarW template 50 to 70",     iso_nume, *rh_endcap_nume_tbarW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_tbarW_50to70   = new RooHistPdf("pdf_barrel_deno_tbarW_50to70",   "Denominator barrel MC tbarW template 50 to 70",   iso_deno, *rh_barrel_deno_tbarW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_tbarW_50to70   = new RooHistPdf("pdf_endcap_deno_tbarW_50to70",   "Denominator endcap MC tbarW template 50 to 70",   iso_deno, *rh_endcap_deno_tbarW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_tbarW_70to100  = new RooHistPdf("pdf_barrel_nume_tbarW_70to100",  "Numerator barrel MC tbarW template 70 to 100",    iso_nume, *rh_barrel_nume_tbarW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_tbarW_70to100  = new RooHistPdf("pdf_endcap_nume_tbarW_70to100",  "Numerator endcap MC tbarW template 70 to 100",    iso_nume, *rh_endcap_nume_tbarW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_tbarW_70to100  = new RooHistPdf("pdf_barrel_deno_tbarW_70to100",  "Denominator barrel MC tbarW template 70 to 100",  iso_deno, *rh_barrel_deno_tbarW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_tbarW_70to100  = new RooHistPdf("pdf_endcap_deno_tbarW_70to100",  "Denominator endcap MC tbarW template 70 to 100",  iso_deno, *rh_endcap_deno_tbarW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_tbarW = new RooHistPdf("pdf_barrel_nume_tbarW", "Numerator barrel MC tbarW template 200 to 500",   iso_nume, *rh_barrel_nume_tbarW, 0);
    RooHistPdf *pdf_endcap_nume_tbarW = new RooHistPdf("pdf_endcap_nume_tbarW", "Numerator endcap MC tbarW template 200 to 500",   iso_nume, *rh_endcap_nume_tbarW, 0);
    RooHistPdf *pdf_barrel_deno_tbarW = new RooHistPdf("pdf_barrel_deno_tbarW", "Denominator barrel MC tbarW template 200 to 500", iso_deno, *rh_barrel_deno_tbarW, 0);
    RooHistPdf *pdf_endcap_deno_tbarW = new RooHistPdf("pdf_endcap_deno_tbarW", "Denominator endcap MC tbarW template 200 to 500", iso_deno, *rh_endcap_deno_tbarW, 0);

    RooHistPdf *pdf_barrel_nume_WW_50to70   = new RooHistPdf("pdf_barrel_nume_WW_50to70",   "Numerator barrel MC WW template 50 to 70",     iso_nume, *rh_barrel_nume_WW_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WW_50to70   = new RooHistPdf("pdf_endcap_nume_WW_50to70",   "Numerator endcap MC WW template 50 to 70",     iso_nume, *rh_endcap_nume_WW_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WW_50to70   = new RooHistPdf("pdf_barrel_deno_WW_50to70",   "Denominator barrel MC WW template 50 to 70",   iso_deno, *rh_barrel_deno_WW_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WW_50to70   = new RooHistPdf("pdf_endcap_deno_WW_50to70",   "Denominator endcap MC WW template 50 to 70",   iso_deno, *rh_endcap_deno_WW_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WW_70to100  = new RooHistPdf("pdf_barrel_nume_WW_70to100",  "Numerator barrel MC WW template 70 to 100",    iso_nume, *rh_barrel_nume_WW_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WW_70to100  = new RooHistPdf("pdf_endcap_nume_WW_70to100",  "Numerator endcap MC WW template 70 to 100",    iso_nume, *rh_endcap_nume_WW_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WW_70to100  = new RooHistPdf("pdf_barrel_deno_WW_70to100",  "Denominator barrel MC WW template 70 to 100",  iso_deno, *rh_barrel_deno_WW_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WW_70to100  = new RooHistPdf("pdf_endcap_deno_WW_70to100",  "Denominator endcap MC WW template 70 to 100",  iso_deno, *rh_endcap_deno_WW_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WW = new RooHistPdf("pdf_barrel_nume_WW", "Numerator barrel MC WW template 200 to 500",   iso_nume, *rh_barrel_nume_WW, 0);
    RooHistPdf *pdf_endcap_nume_WW = new RooHistPdf("pdf_endcap_nume_WW", "Numerator endcap MC WW template 200 to 500",   iso_nume, *rh_endcap_nume_WW, 0);
    RooHistPdf *pdf_barrel_deno_WW = new RooHistPdf("pdf_barrel_deno_WW", "Denominator barrel MC WW template 200 to 500", iso_deno, *rh_barrel_deno_WW, 0);
    RooHistPdf *pdf_endcap_deno_WW = new RooHistPdf("pdf_endcap_deno_WW", "Denominator endcap MC WW template 200 to 500", iso_deno, *rh_endcap_deno_WW, 0);

    RooHistPdf *pdf_barrel_nume_WZ_50to70   = new RooHistPdf("pdf_barrel_nume_WZ_50to70",   "Numerator barrel MC WZ template 50 to 70",     iso_nume, *rh_barrel_nume_WZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_WZ_50to70   = new RooHistPdf("pdf_endcap_nume_WZ_50to70",   "Numerator endcap MC WZ template 50 to 70",     iso_nume, *rh_endcap_nume_WZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_WZ_50to70   = new RooHistPdf("pdf_barrel_deno_WZ_50to70",   "Denominator barrel MC WZ template 50 to 70",   iso_deno, *rh_barrel_deno_WZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_WZ_50to70   = new RooHistPdf("pdf_endcap_deno_WZ_50to70",   "Denominator endcap MC WZ template 50 to 70",   iso_deno, *rh_endcap_deno_WZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_WZ_70to100  = new RooHistPdf("pdf_barrel_nume_WZ_70to100",  "Numerator barrel MC WZ template 70 to 100",    iso_nume, *rh_barrel_nume_WZ_70to100,  0);
    RooHistPdf *pdf_endcap_nume_WZ_70to100  = new RooHistPdf("pdf_endcap_nume_WZ_70to100",  "Numerator endcap MC WZ template 70 to 100",    iso_nume, *rh_endcap_nume_WZ_70to100,  0);
    RooHistPdf *pdf_barrel_deno_WZ_70to100  = new RooHistPdf("pdf_barrel_deno_WZ_70to100",  "Denominator barrel MC WZ template 70 to 100",  iso_deno, *rh_barrel_deno_WZ_70to100,  0);
    RooHistPdf *pdf_endcap_deno_WZ_70to100  = new RooHistPdf("pdf_endcap_deno_WZ_70to100",  "Denominator endcap MC WZ template 70 to 100",  iso_deno, *rh_endcap_deno_WZ_70to100,  0);
    RooHistPdf *pdf_barrel_nume_WZ = new RooHistPdf("pdf_barrel_nume_WZ", "Numerator barrel MC WZ template 200 to 500",   iso_nume, *rh_barrel_nume_WZ, 0);
    RooHistPdf *pdf_endcap_nume_WZ = new RooHistPdf("pdf_endcap_nume_WZ", "Numerator endcap MC WZ template 200 to 500",   iso_nume, *rh_endcap_nume_WZ, 0);
    RooHistPdf *pdf_barrel_deno_WZ = new RooHistPdf("pdf_barrel_deno_WZ", "Denominator barrel MC WZ template 200 to 500", iso_deno, *rh_barrel_deno_WZ, 0);
    RooHistPdf *pdf_endcap_deno_WZ = new RooHistPdf("pdf_endcap_deno_WZ", "Denominator endcap MC WZ template 200 to 500", iso_deno, *rh_endcap_deno_WZ, 0);

    RooHistPdf *pdf_barrel_nume_ZZ_50to70   = new RooHistPdf("pdf_barrel_nume_ZZ_50to70",   "Numerator barrel MC ZZ template 50 to 70",     iso_nume, *rh_barrel_nume_ZZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_nume_ZZ_50to70   = new RooHistPdf("pdf_endcap_nume_ZZ_50to70",   "Numerator endcap MC ZZ template 50 to 70",     iso_nume, *rh_endcap_nume_ZZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_deno_ZZ_50to70   = new RooHistPdf("pdf_barrel_deno_ZZ_50to70",   "Denominator barrel MC ZZ template 50 to 70",   iso_deno, *rh_barrel_deno_ZZ_50to70 ,  0);
    RooHistPdf *pdf_endcap_deno_ZZ_50to70   = new RooHistPdf("pdf_endcap_deno_ZZ_50to70",   "Denominator endcap MC ZZ template 50 to 70",   iso_deno, *rh_endcap_deno_ZZ_50to70 ,  0);
    RooHistPdf *pdf_barrel_nume_ZZ_70to100  = new RooHistPdf("pdf_barrel_nume_ZZ_70to100",  "Numerator barrel MC ZZ template 70 to 100",    iso_nume, *rh_barrel_nume_ZZ_70to100,  0);
    RooHistPdf *pdf_endcap_nume_ZZ_70to100  = new RooHistPdf("pdf_endcap_nume_ZZ_70to100",  "Numerator endcap MC ZZ template 70 to 100",    iso_nume, *rh_endcap_nume_ZZ_70to100,  0);
    RooHistPdf *pdf_barrel_deno_ZZ_70to100  = new RooHistPdf("pdf_barrel_deno_ZZ_70to100",  "Denominator barrel MC ZZ template 70 to 100",  iso_deno, *rh_barrel_deno_ZZ_70to100,  0);
    RooHistPdf *pdf_endcap_deno_ZZ_70to100  = new RooHistPdf("pdf_endcap_deno_ZZ_70to100",  "Denominator endcap MC ZZ template 70 to 100",  iso_deno, *rh_endcap_deno_ZZ_70to100,  0);
    RooHistPdf *pdf_barrel_nume_ZZ = new RooHistPdf("pdf_barrel_nume_ZZ", "Numerator barrel MC ZZ template 200 to 500",   iso_nume, *rh_barrel_nume_ZZ, 0);
    RooHistPdf *pdf_endcap_nume_ZZ = new RooHistPdf("pdf_endcap_nume_ZZ", "Numerator endcap MC ZZ template 200 to 500",   iso_nume, *rh_endcap_nume_ZZ, 0);
    RooHistPdf *pdf_barrel_deno_ZZ = new RooHistPdf("pdf_barrel_deno_ZZ", "Denominator barrel MC ZZ template 200 to 500", iso_deno, *rh_barrel_deno_ZZ, 0);
    RooHistPdf *pdf_endcap_deno_ZZ = new RooHistPdf("pdf_endcap_deno_ZZ", "Denominator endcap MC ZZ template 200 to 500", iso_deno, *rh_endcap_deno_ZZ, 0);

    // Constraints for integrals
    Double_t N_barrel_nume_ttbar_50to70   = h_barrel_MC_nume_50to70  [_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar_50to70   = h_endcap_MC_nume_50to70  [_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar_50to70   = h_barrel_MC_deno_50to70  [_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar_50to70   = h_endcap_MC_deno_50to70  [_ttbar]->Integral();
    Double_t N_barrel_nume_ttbar_70to100  = h_barrel_MC_nume_70to100 [_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar_70to100  = h_endcap_MC_nume_70to100 [_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar_70to100  = h_barrel_MC_deno_70to100 [_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar_70to100  = h_endcap_MC_deno_70to100 [_ttbar]->Integral();
    Double_t N_barrel_nume_ttbar = h_barrel_MC_nume[_ttbar]->Integral();
    Double_t N_endcap_nume_ttbar = h_endcap_MC_nume[_ttbar]->Integral();
    Double_t N_barrel_deno_ttbar = h_barrel_MC_deno[_ttbar]->Integral();
    Double_t N_endcap_deno_ttbar = h_endcap_MC_deno[_ttbar]->Integral();

    Double_t N_barrel_nume_WJets_50to70   = h_barrel_MC_nume_50to70  [_WJets]->Integral();
    Double_t N_endcap_nume_WJets_50to70   = h_endcap_MC_nume_50to70  [_WJets]->Integral();
    Double_t N_barrel_deno_WJets_50to70   = h_barrel_MC_deno_50to70  [_WJets]->Integral();
    Double_t N_endcap_deno_WJets_50to70   = h_endcap_MC_deno_50to70  [_WJets]->Integral();
    Double_t N_barrel_nume_WJets_70to100  = h_barrel_MC_nume_70to100 [_WJets]->Integral();
    Double_t N_endcap_nume_WJets_70to100  = h_endcap_MC_nume_70to100 [_WJets]->Integral();
    Double_t N_barrel_deno_WJets_70to100  = h_barrel_MC_deno_70to100 [_WJets]->Integral();
    Double_t N_endcap_deno_WJets_70to100  = h_endcap_MC_deno_70to100 [_WJets]->Integral();
    Double_t N_barrel_nume_WJets = h_barrel_MC_nume[_WJets]->Integral();
    Double_t N_endcap_nume_WJets = h_endcap_MC_nume[_WJets]->Integral();
    Double_t N_barrel_deno_WJets = h_barrel_MC_deno[_WJets]->Integral();
    Double_t N_endcap_deno_WJets = h_endcap_MC_deno[_WJets]->Integral();

    Double_t N_barrel_nume_DY_50to70   = h_barrel_MC_nume_50to70  [_DY_Full]->Integral();
    Double_t N_endcap_nume_DY_50to70   = h_endcap_MC_nume_50to70  [_DY_Full]->Integral();
    Double_t N_barrel_deno_DY_50to70   = h_barrel_MC_deno_50to70  [_DY_Full]->Integral();
    Double_t N_endcap_deno_DY_50to70   = h_endcap_MC_deno_50to70  [_DY_Full]->Integral();
    Double_t N_barrel_nume_DY_70to100  = h_barrel_MC_nume_70to100 [_DY_Full]->Integral();
    Double_t N_endcap_nume_DY_70to100  = h_endcap_MC_nume_70to100 [_DY_Full]->Integral();
    Double_t N_barrel_deno_DY_70to100  = h_barrel_MC_deno_70to100 [_DY_Full]->Integral();
    Double_t N_endcap_deno_DY_70to100  = h_endcap_MC_deno_70to100 [_DY_Full]->Integral();
    Double_t N_barrel_nume_DY = h_barrel_MC_nume[_DY_Full]->Integral();
    Double_t N_endcap_nume_DY = h_endcap_MC_nume[_DY_Full]->Integral();
    Double_t N_barrel_deno_DY = h_barrel_MC_deno[_DY_Full]->Integral();
    Double_t N_endcap_deno_DY = h_endcap_MC_deno[_DY_Full]->Integral();

    Double_t N_barrel_nume_QCD_50to70   = h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD_50to70   = h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD_50to70   = h_barrel_MC_deno_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD_50to70   = h_endcap_MC_deno_50to70  [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_nume_QCD_70to100  = h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD_70to100  = h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD_70to100  = h_barrel_MC_deno_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD_70to100  = h_endcap_MC_deno_70to100 [_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_nume_QCD = h_barrel_MC_nume[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_nume_QCD = h_endcap_MC_nume[_QCDMuEnriched_Full]->Integral();
    Double_t N_barrel_deno_QCD = h_barrel_MC_deno[_QCDMuEnriched_Full]->Integral();
    Double_t N_endcap_deno_QCD = h_endcap_MC_deno[_QCDMuEnriched_Full]->Integral();

    Double_t N_barrel_nume_tW_50to70   = h_barrel_MC_nume_50to70  [_tW]->Integral();
    Double_t N_endcap_nume_tW_50to70   = h_endcap_MC_nume_50to70  [_tW]->Integral();
    Double_t N_barrel_deno_tW_50to70   = h_barrel_MC_deno_50to70  [_tW]->Integral();
    Double_t N_endcap_deno_tW_50to70   = h_endcap_MC_deno_50to70  [_tW]->Integral();
    Double_t N_barrel_nume_tW_70to100  = h_barrel_MC_nume_70to100 [_tW]->Integral();
    Double_t N_endcap_nume_tW_70to100  = h_endcap_MC_nume_70to100 [_tW]->Integral();
    Double_t N_barrel_deno_tW_70to100  = h_barrel_MC_deno_70to100 [_tW]->Integral();
    Double_t N_endcap_deno_tW_70to100  = h_endcap_MC_deno_70to100 [_tW]->Integral();
    Double_t N_barrel_nume_tW = h_barrel_MC_nume[_tW]->Integral();
    Double_t N_endcap_nume_tW = h_endcap_MC_nume[_tW]->Integral();
    Double_t N_barrel_deno_tW = h_barrel_MC_deno[_tW]->Integral();
    Double_t N_endcap_deno_tW = h_endcap_MC_deno[_tW]->Integral();

    Double_t N_barrel_nume_tbarW_50to70   = h_barrel_MC_nume_50to70  [_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW_50to70   = h_endcap_MC_nume_50to70  [_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW_50to70   = h_barrel_MC_deno_50to70  [_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW_50to70   = h_endcap_MC_deno_50to70  [_tbarW]->Integral();
    Double_t N_barrel_nume_tbarW_70to100  = h_barrel_MC_nume_70to100 [_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW_70to100  = h_endcap_MC_nume_70to100 [_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW_70to100  = h_barrel_MC_deno_70to100 [_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW_70to100  = h_endcap_MC_deno_70to100 [_tbarW]->Integral();
    Double_t N_barrel_nume_tbarW = h_barrel_MC_nume[_tbarW]->Integral();
    Double_t N_endcap_nume_tbarW = h_endcap_MC_nume[_tbarW]->Integral();
    Double_t N_barrel_deno_tbarW = h_barrel_MC_deno[_tbarW]->Integral();
    Double_t N_endcap_deno_tbarW = h_endcap_MC_deno[_tbarW]->Integral();

    Double_t N_barrel_nume_WW_50to70   = h_barrel_MC_nume_50to70  [_WW]->Integral();
    Double_t N_endcap_nume_WW_50to70   = h_endcap_MC_nume_50to70  [_WW]->Integral();
    Double_t N_barrel_deno_WW_50to70   = h_barrel_MC_deno_50to70  [_WW]->Integral();
    Double_t N_endcap_deno_WW_50to70   = h_endcap_MC_deno_50to70  [_WW]->Integral();
    Double_t N_barrel_nume_WW_70to100  = h_barrel_MC_nume_70to100 [_WW]->Integral();
    Double_t N_endcap_nume_WW_70to100  = h_endcap_MC_nume_70to100 [_WW]->Integral();
    Double_t N_barrel_deno_WW_70to100  = h_barrel_MC_deno_70to100 [_WW]->Integral();
    Double_t N_endcap_deno_WW_70to100  = h_endcap_MC_deno_70to100 [_WW]->Integral();
    Double_t N_barrel_nume_WW = h_barrel_MC_nume[_WW]->Integral();
    Double_t N_endcap_nume_WW = h_endcap_MC_nume[_WW]->Integral();
    Double_t N_barrel_deno_WW = h_barrel_MC_deno[_WW]->Integral();
    Double_t N_endcap_deno_WW = h_endcap_MC_deno[_WW]->Integral();

    Double_t N_barrel_nume_WZ_50to70   = h_barrel_MC_nume_50to70  [_WZ]->Integral();
    Double_t N_endcap_nume_WZ_50to70   = h_endcap_MC_nume_50to70  [_WZ]->Integral();
    Double_t N_barrel_deno_WZ_50to70   = h_barrel_MC_deno_50to70  [_WZ]->Integral();
    Double_t N_endcap_deno_WZ_50to70   = h_endcap_MC_deno_50to70  [_WZ]->Integral();
    Double_t N_barrel_nume_WZ_70to100  = h_barrel_MC_nume_70to100 [_WZ]->Integral();
    Double_t N_endcap_nume_WZ_70to100  = h_endcap_MC_nume_70to100 [_WZ]->Integral();
    Double_t N_barrel_deno_WZ_70to100  = h_barrel_MC_deno_70to100 [_WZ]->Integral();
    Double_t N_endcap_deno_WZ_70to100  = h_endcap_MC_deno_70to100 [_WZ]->Integral();
    Double_t N_barrel_nume_WZ = h_barrel_MC_nume[_WZ]->Integral();
    Double_t N_endcap_nume_WZ = h_endcap_MC_nume[_WZ]->Integral();
    Double_t N_barrel_deno_WZ = h_barrel_MC_deno[_WZ]->Integral();
    Double_t N_endcap_deno_WZ = h_endcap_MC_deno[_WZ]->Integral();

    Double_t N_barrel_nume_ZZ_50to70   = h_barrel_MC_nume_50to70  [_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ_50to70   = h_endcap_MC_nume_50to70  [_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ_50to70   = h_barrel_MC_deno_50to70  [_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ_50to70   = h_endcap_MC_deno_50to70  [_ZZ]->Integral();
    Double_t N_barrel_nume_ZZ_70to100  = h_barrel_MC_nume_70to100 [_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ_70to100  = h_endcap_MC_nume_70to100 [_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ_70to100  = h_barrel_MC_deno_70to100 [_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ_70to100  = h_endcap_MC_deno_70to100 [_ZZ]->Integral();
    Double_t N_barrel_nume_ZZ = h_barrel_MC_nume[_ZZ]->Integral();
    Double_t N_endcap_nume_ZZ = h_endcap_MC_nume[_ZZ]->Integral();
    Double_t N_barrel_deno_ZZ = h_barrel_MC_deno[_ZZ]->Integral();
    Double_t N_endcap_deno_ZZ = h_endcap_MC_deno[_ZZ]->Integral();

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
    Double_t N_barrel_nume_total = N_barrel_nume_ttbar + N_barrel_nume_WJets + N_barrel_nume_DY + N_barrel_nume_QCD +
                                            N_barrel_nume_tW    + N_barrel_nume_tbarW + N_barrel_nume_WW + N_barrel_nume_WZ  +
                                            N_barrel_nume_ZZ;
    Double_t N_endcap_nume_total = N_endcap_nume_ttbar + N_endcap_nume_WJets + N_endcap_nume_DY + N_endcap_nume_QCD +
                                            N_endcap_nume_tW    + N_endcap_nume_tbarW + N_endcap_nume_WW + N_endcap_nume_WZ  +
                                            N_endcap_nume_ZZ;
    Double_t N_barrel_deno_total = N_barrel_deno_ttbar + N_barrel_deno_WJets + N_barrel_deno_DY + N_barrel_deno_QCD +
                                            N_barrel_deno_tW    + N_barrel_deno_tbarW + N_barrel_deno_WW + N_barrel_deno_WZ  +
                                            N_barrel_deno_ZZ;
    Double_t N_endcap_deno_total = N_endcap_deno_ttbar + N_endcap_deno_WJets + N_endcap_deno_DY + N_endcap_deno_QCD +
                                            N_endcap_deno_tW    + N_endcap_deno_tbarW + N_endcap_deno_WW + N_endcap_deno_WZ  +
                                            N_endcap_deno_ZZ;

    Double_t Nnorm_barrel_nume_ttbar_50to70   = N_barrel_nume_ttbar_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_ttbar_50to70   = N_endcap_nume_ttbar_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_ttbar_50to70   = N_barrel_deno_ttbar_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_ttbar_50to70   = N_endcap_deno_ttbar_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_ttbar_70to100  = N_barrel_nume_ttbar_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_ttbar_70to100  = N_endcap_nume_ttbar_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_ttbar_70to100  = N_barrel_deno_ttbar_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_ttbar_70to100  = N_endcap_deno_ttbar_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_ttbar = N_barrel_nume_ttbar * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_ttbar = N_endcap_nume_ttbar * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_ttbar = N_barrel_deno_ttbar * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_ttbar = N_endcap_deno_ttbar * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WJets_50to70   = N_barrel_nume_WJets_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WJets_50to70   = N_endcap_nume_WJets_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WJets_50to70   = N_barrel_deno_WJets_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WJets_50to70   = N_endcap_deno_WJets_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WJets_70to100  = N_barrel_nume_WJets_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WJets_70to100  = N_endcap_nume_WJets_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WJets_70to100  = N_barrel_deno_WJets_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WJets_70to100  = N_endcap_deno_WJets_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WJets = N_barrel_nume_WJets * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WJets = N_endcap_nume_WJets * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WJets = N_barrel_deno_WJets * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WJets = N_endcap_deno_WJets * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_DY_50to70   = N_barrel_nume_DY_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_DY_50to70   = N_endcap_nume_DY_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_DY_50to70   = N_barrel_deno_DY_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_DY_50to70   = N_endcap_deno_DY_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_DY_70to100  = N_barrel_nume_DY_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_DY_70to100  = N_endcap_nume_DY_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_DY_70to100  = N_barrel_deno_DY_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_DY_70to100  = N_endcap_deno_DY_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_DY = N_barrel_nume_DY * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_DY = N_endcap_nume_DY * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_DY = N_barrel_deno_DY * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_DY = N_endcap_deno_DY * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_tW_50to70   = N_barrel_nume_tW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_tW_50to70   = N_endcap_nume_tW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_tW_50to70   = N_barrel_deno_tW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_tW_50to70   = N_endcap_deno_tW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_tW_70to100  = N_barrel_nume_tW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_tW_70to100  = N_endcap_nume_tW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_tW_70to100  = N_barrel_deno_tW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_tW_70to100  = N_endcap_deno_tW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_tW = N_barrel_nume_tW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_tW = N_endcap_nume_tW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_tW = N_barrel_deno_tW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_tW = N_endcap_deno_tW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_tbarW_50to70   = N_barrel_nume_tbarW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_tbarW_50to70   = N_endcap_nume_tbarW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_tbarW_50to70   = N_barrel_deno_tbarW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_tbarW_50to70   = N_endcap_deno_tbarW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_tbarW_70to100  = N_barrel_nume_tbarW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_tbarW_70to100  = N_endcap_nume_tbarW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_tbarW_70to100  = N_barrel_deno_tbarW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_tbarW_70to100  = N_endcap_deno_tbarW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_tbarW = N_barrel_nume_tbarW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_tbarW = N_endcap_nume_tbarW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_tbarW = N_barrel_deno_tbarW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_tbarW = N_endcap_deno_tbarW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WW_50to70   = N_barrel_nume_WW_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WW_50to70   = N_endcap_nume_WW_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WW_50to70   = N_barrel_deno_WW_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WW_50to70   = N_endcap_deno_WW_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WW_70to100  = N_barrel_nume_WW_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WW_70to100  = N_endcap_nume_WW_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WW_70to100  = N_barrel_deno_WW_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WW_70to100  = N_endcap_deno_WW_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WW = N_barrel_nume_WW * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WW = N_endcap_nume_WW * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WW = N_barrel_deno_WW * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WW = N_endcap_deno_WW * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_WZ_50to70   = N_barrel_nume_WZ_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_WZ_50to70   = N_endcap_nume_WZ_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_WZ_50to70   = N_barrel_deno_WZ_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_WZ_50to70   = N_endcap_deno_WZ_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_WZ_70to100  = N_barrel_nume_WZ_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_WZ_70to100  = N_endcap_nume_WZ_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_WZ_70to100  = N_barrel_deno_WZ_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_WZ_70to100  = N_endcap_deno_WZ_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_WZ = N_barrel_nume_WZ * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_WZ = N_endcap_nume_WZ * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_WZ = N_barrel_deno_WZ * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_WZ = N_endcap_deno_WZ * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_ZZ_50to70   = N_barrel_nume_ZZ_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_ZZ_50to70   = N_endcap_nume_ZZ_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_ZZ_50to70   = N_barrel_deno_ZZ_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_ZZ_50to70   = N_endcap_deno_ZZ_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_ZZ_70to100  = N_barrel_nume_ZZ_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_ZZ_70to100  = N_endcap_nume_ZZ_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_ZZ_70to100  = N_barrel_deno_ZZ_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_ZZ_70to100  = N_endcap_deno_ZZ_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_ZZ = N_barrel_nume_ZZ * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_ZZ = N_endcap_nume_ZZ * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_ZZ = N_barrel_deno_ZZ * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_ZZ = N_endcap_deno_ZZ * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    Double_t Nnorm_barrel_nume_QCD_50to70   = N_barrel_nume_QCD_50to70   * h_barrel_data_nume_50to70  ->Integral() / N_barrel_nume_total_50to70 ;
    Double_t Nnorm_endcap_nume_QCD_50to70   = N_endcap_nume_QCD_50to70   * h_endcap_data_nume_50to70  ->Integral() / N_endcap_nume_total_50to70 ;
    Double_t Nnorm_barrel_deno_QCD_50to70   = N_barrel_deno_QCD_50to70   * h_barrel_data_deno_50to70  ->Integral() / N_barrel_deno_total_50to70 ;
    Double_t Nnorm_endcap_deno_QCD_50to70   = N_endcap_deno_QCD_50to70   * h_endcap_data_deno_50to70  ->Integral() / N_endcap_deno_total_50to70 ;
    Double_t Nnorm_barrel_nume_QCD_70to100  = N_barrel_nume_QCD_70to100  * h_barrel_data_nume_70to100 ->Integral() / N_barrel_nume_total_70to100;
    Double_t Nnorm_endcap_nume_QCD_70to100  = N_endcap_nume_QCD_70to100  * h_endcap_data_nume_70to100 ->Integral() / N_endcap_nume_total_70to100;
    Double_t Nnorm_barrel_deno_QCD_70to100  = N_barrel_deno_QCD_70to100  * h_barrel_data_deno_70to100 ->Integral() / N_barrel_deno_total_70to100;
    Double_t Nnorm_endcap_deno_QCD_70to100  = N_endcap_deno_QCD_70to100  * h_endcap_data_deno_70to100 ->Integral() / N_endcap_deno_total_70to100;
    Double_t Nnorm_barrel_nume_QCD = N_barrel_nume_QCD * h_barrel_data_nume->Integral() / N_barrel_nume_total;
    Double_t Nnorm_endcap_nume_QCD = N_endcap_nume_QCD * h_endcap_data_nume->Integral() / N_endcap_nume_total;
    Double_t Nnorm_barrel_deno_QCD = N_barrel_deno_QCD * h_barrel_data_deno->Integral() / N_barrel_deno_total;
    Double_t Nnorm_endcap_deno_QCD = N_endcap_deno_QCD * h_endcap_data_deno->Integral() / N_endcap_deno_total;

    // Fit constraints
    RooRealVar n_barrel_nume_ttbar_50to70  ("n_barrel_nume_ttbar_50to70",   "n_barrel_nume_ttbar_50to70",   Nnorm_barrel_nume_ttbar_50to70 ,  Nnorm_barrel_nume_ttbar_50to70  *0.75, Nnorm_barrel_nume_ttbar_50to70  *1.25);
    RooRealVar n_endcap_nume_ttbar_50to70  ("n_endcap_nume_ttbar_50to70",   "n_endcap_nume_ttbar_50to70",   Nnorm_endcap_nume_ttbar_50to70 ,  Nnorm_endcap_nume_ttbar_50to70  *0.75, Nnorm_endcap_nume_ttbar_50to70  *1.25);
    RooRealVar n_barrel_deno_ttbar_50to70  ("n_barrel_deno_ttbar_50to70",   "n_barrel_deno_ttbar_50to70",   Nnorm_barrel_deno_ttbar_50to70 ,  Nnorm_barrel_deno_ttbar_50to70  *0.75, Nnorm_barrel_deno_ttbar_50to70  *1.25);
    RooRealVar n_endcap_deno_ttbar_50to70  ("n_endcap_deno_ttbar_50to70",   "n_endcap_deno_ttbar_50to70",   Nnorm_endcap_deno_ttbar_50to70 ,  Nnorm_endcap_deno_ttbar_50to70  *0.75, Nnorm_endcap_deno_ttbar_50to70  *1.25);
    RooRealVar n_barrel_nume_ttbar_70to100 ("n_barrel_nume_ttbar_70to100",  "n_barrel_nume_ttbar_70to100",  Nnorm_barrel_nume_ttbar_70to100,  Nnorm_barrel_nume_ttbar_70to100 *0.75, Nnorm_barrel_nume_ttbar_70to100 *1.25);
    RooRealVar n_endcap_nume_ttbar_70to100 ("n_endcap_nume_ttbar_70to100",  "n_endcap_nume_ttbar_70to100",  Nnorm_endcap_nume_ttbar_70to100,  Nnorm_endcap_nume_ttbar_70to100 *0.75, Nnorm_endcap_nume_ttbar_70to100 *1.25);
    RooRealVar n_barrel_deno_ttbar_70to100 ("n_barrel_deno_ttbar_70to100",  "n_barrel_deno_ttbar_70to100",  Nnorm_barrel_deno_ttbar_70to100,  Nnorm_barrel_deno_ttbar_70to100 *0.75, Nnorm_barrel_deno_ttbar_70to100 *1.25);
    RooRealVar n_endcap_deno_ttbar_70to100 ("n_endcap_deno_ttbar_70to100",  "n_endcap_deno_ttbar_70to100",  Nnorm_endcap_deno_ttbar_70to100,  Nnorm_endcap_deno_ttbar_70to100 *0.75, Nnorm_endcap_deno_ttbar_70to100 *1.25);
    RooRealVar n_barrel_nume_ttbar("n_barrel_nume_ttbar", "n_barrel_nume_ttbar", Nnorm_barrel_nume_ttbar, Nnorm_barrel_nume_ttbar*0.75, Nnorm_barrel_nume_ttbar*1.25);
    RooRealVar n_endcap_nume_ttbar("n_endcap_nume_ttbar", "n_endcap_nume_ttbar", Nnorm_endcap_nume_ttbar, Nnorm_endcap_nume_ttbar*0.75, Nnorm_endcap_nume_ttbar*1.25);
    RooRealVar n_barrel_deno_ttbar("n_barrel_deno_ttbar", "n_barrel_deno_ttbar", Nnorm_barrel_deno_ttbar, Nnorm_barrel_deno_ttbar*0.75, Nnorm_barrel_deno_ttbar*1.25);
    RooRealVar n_endcap_deno_ttbar("n_endcap_deno_ttbar", "n_endcap_deno_ttbar", Nnorm_endcap_deno_ttbar, Nnorm_endcap_deno_ttbar*0.75, Nnorm_endcap_deno_ttbar*1.25);

    RooRealVar n_barrel_nume_WJets_50to70  ("n_barrel_nume_WJets_50to70",   "n_barrel_nume_WJets_50to70",   Nnorm_barrel_nume_WJets_50to70 ,  Nnorm_barrel_nume_WJets_50to70  *0.5, Nnorm_barrel_nume_WJets_50to70  *1.5);
    RooRealVar n_endcap_nume_WJets_50to70  ("n_endcap_nume_WJets_50to70",   "n_endcap_nume_WJets_50to70",   Nnorm_endcap_nume_WJets_50to70 ,  Nnorm_endcap_nume_WJets_50to70  *0.5, Nnorm_endcap_nume_WJets_50to70  *1.5);
    RooRealVar n_barrel_deno_WJets_50to70  ("n_barrel_deno_WJets_50to70",   "n_barrel_deno_WJets_50to70",   Nnorm_barrel_deno_WJets_50to70 ,  Nnorm_barrel_deno_WJets_50to70  *0.5, Nnorm_barrel_deno_WJets_50to70  *1.5);
    RooRealVar n_endcap_deno_WJets_50to70  ("n_endcap_deno_WJets_50to70",   "n_endcap_deno_WJets_50to70",   Nnorm_endcap_deno_WJets_50to70 ,  Nnorm_endcap_deno_WJets_50to70  *0.5, Nnorm_endcap_deno_WJets_50to70  *1.5);
    RooRealVar n_barrel_nume_WJets_70to100 ("n_barrel_nume_WJets_70to100",  "n_barrel_nume_WJets_70to100",  Nnorm_barrel_nume_WJets_70to100,  Nnorm_barrel_nume_WJets_70to100 *0.5, Nnorm_barrel_nume_WJets_70to100 *1.5);
    RooRealVar n_endcap_nume_WJets_70to100 ("n_endcap_nume_WJets_70to100",  "n_endcap_nume_WJets_70to100",  Nnorm_endcap_nume_WJets_70to100,  Nnorm_endcap_nume_WJets_70to100 *0.5, Nnorm_endcap_nume_WJets_70to100 *1.5);
    RooRealVar n_barrel_deno_WJets_70to100 ("n_barrel_deno_WJets_70to100",  "n_barrel_deno_WJets_70to100",  Nnorm_barrel_deno_WJets_70to100,  Nnorm_barrel_deno_WJets_70to100 *0.5, Nnorm_barrel_deno_WJets_70to100 *1.5);
    RooRealVar n_endcap_deno_WJets_70to100 ("n_endcap_deno_WJets_70to100",  "n_endcap_deno_WJets_70to100",  Nnorm_endcap_deno_WJets_70to100,  Nnorm_endcap_deno_WJets_70to100 *0.5, Nnorm_endcap_deno_WJets_70to100 *1.5);
    RooRealVar n_barrel_nume_WJets("n_barrel_nume_WJets", "n_barrel_nume_WJets", Nnorm_barrel_nume_WJets, Nnorm_barrel_nume_WJets*0.5, Nnorm_barrel_nume_WJets*1.5);
    RooRealVar n_endcap_nume_WJets("n_endcap_nume_WJets", "n_endcap_nume_WJets", Nnorm_endcap_nume_WJets, Nnorm_endcap_nume_WJets*0.5, Nnorm_endcap_nume_WJets*1.5);
    RooRealVar n_barrel_deno_WJets("n_barrel_deno_WJets", "n_barrel_deno_WJets", Nnorm_barrel_deno_WJets, Nnorm_barrel_deno_WJets*0.5, Nnorm_barrel_deno_WJets*1.5);
    RooRealVar n_endcap_deno_WJets("n_endcap_deno_WJets", "n_endcap_deno_WJets", Nnorm_endcap_deno_WJets, Nnorm_endcap_deno_WJets*0.5, Nnorm_endcap_deno_WJets*1.5);

    RooRealVar n_barrel_nume_DY_50to70  ("n_barrel_nume_DY_50to70",   "n_barrel_nume_DY_50to70",   Nnorm_barrel_nume_DY_50to70 ,  Nnorm_barrel_nume_DY_50to70  *0.8, Nnorm_barrel_nume_DY_50to70  *1.2);
    RooRealVar n_endcap_nume_DY_50to70  ("n_endcap_nume_DY_50to70",   "n_endcap_nume_DY_50to70",   Nnorm_endcap_nume_DY_50to70 ,  Nnorm_endcap_nume_DY_50to70  *0.8, Nnorm_endcap_nume_DY_50to70  *1.2);
    RooRealVar n_barrel_deno_DY_50to70  ("n_barrel_deno_DY_50to70",   "n_barrel_deno_DY_50to70",   Nnorm_barrel_deno_DY_50to70 ,  Nnorm_barrel_deno_DY_50to70  *0.8, Nnorm_barrel_deno_DY_50to70  *1.2);
    RooRealVar n_endcap_deno_DY_50to70  ("n_endcap_deno_DY_50to70",   "n_endcap_deno_DY_50to70",   Nnorm_endcap_deno_DY_50to70 ,  Nnorm_endcap_deno_DY_50to70  *0.8, Nnorm_endcap_deno_DY_50to70  *1.2);
    RooRealVar n_barrel_nume_DY_70to100 ("n_barrel_nume_DY_70to100",  "n_barrel_nume_DY_70to100",  Nnorm_barrel_nume_DY_70to100,  Nnorm_barrel_nume_DY_70to100 *0.8, Nnorm_barrel_nume_DY_70to100 *1.2);
    RooRealVar n_endcap_nume_DY_70to100 ("n_endcap_nume_DY_70to100",  "n_endcap_nume_DY_70to100",  Nnorm_endcap_nume_DY_70to100,  Nnorm_endcap_nume_DY_70to100 *0.8, Nnorm_endcap_nume_DY_70to100 *1.2);
    RooRealVar n_barrel_deno_DY_70to100 ("n_barrel_deno_DY_70to100",  "n_barrel_deno_DY_70to100",  Nnorm_barrel_deno_DY_70to100,  Nnorm_barrel_deno_DY_70to100 *0.8, Nnorm_barrel_deno_DY_70to100 *1.2);
    RooRealVar n_endcap_deno_DY_70to100 ("n_endcap_deno_DY_70to100",  "n_endcap_deno_DY_70to100",  Nnorm_endcap_deno_DY_70to100,  Nnorm_endcap_deno_DY_70to100 *0.8, Nnorm_endcap_deno_DY_70to100 *1.2);
    RooRealVar n_barrel_nume_DY("n_barrel_nume_DY", "n_barrel_nume_DY", Nnorm_barrel_nume_DY, Nnorm_barrel_nume_DY*0.8, Nnorm_barrel_nume_DY*1.2);
    RooRealVar n_endcap_nume_DY("n_endcap_nume_DY", "n_endcap_nume_DY", Nnorm_endcap_nume_DY, Nnorm_endcap_nume_DY*0.8, Nnorm_endcap_nume_DY*1.2);
    RooRealVar n_barrel_deno_DY("n_barrel_deno_DY", "n_barrel_deno_DY", Nnorm_barrel_deno_DY, Nnorm_barrel_deno_DY*0.8, Nnorm_barrel_deno_DY*1.2);
    RooRealVar n_endcap_deno_DY("n_endcap_deno_DY", "n_endcap_deno_DY", Nnorm_endcap_deno_DY, Nnorm_endcap_deno_DY*0.8, Nnorm_endcap_deno_DY*1.2);

    RooRealVar n_barrel_nume_tW_50to70  ("n_barrel_nume_tW_50to70",   "n_barrel_nume_tW_50to70",   Nnorm_barrel_nume_tW_50to70 ,  Nnorm_barrel_nume_tW_50to70  *0.8, Nnorm_barrel_nume_tW_50to70  *1.2);
    RooRealVar n_endcap_nume_tW_50to70  ("n_endcap_nume_tW_50to70",   "n_endcap_nume_tW_50to70",   Nnorm_endcap_nume_tW_50to70 ,  Nnorm_endcap_nume_tW_50to70  *0.8, Nnorm_endcap_nume_tW_50to70  *1.2);
    RooRealVar n_barrel_deno_tW_50to70  ("n_barrel_deno_tW_50to70",   "n_barrel_deno_tW_50to70",   Nnorm_barrel_deno_tW_50to70 ,  Nnorm_barrel_deno_tW_50to70  *0.8, Nnorm_barrel_deno_tW_50to70  *1.2);
    RooRealVar n_endcap_deno_tW_50to70  ("n_endcap_deno_tW_50to70",   "n_endcap_deno_tW_50to70",   Nnorm_endcap_deno_tW_50to70 ,  Nnorm_endcap_deno_tW_50to70  *0.8, Nnorm_endcap_deno_tW_50to70  *1.2);
    RooRealVar n_barrel_nume_tW_70to100 ("n_barrel_nume_tW_70to100",  "n_barrel_nume_tW_70to100",  Nnorm_barrel_nume_tW_70to100,  Nnorm_barrel_nume_tW_70to100 *0.8, Nnorm_barrel_nume_tW_70to100 *1.2);
    RooRealVar n_endcap_nume_tW_70to100 ("n_endcap_nume_tW_70to100",  "n_endcap_nume_tW_70to100",  Nnorm_endcap_nume_tW_70to100,  Nnorm_endcap_nume_tW_70to100 *0.8, Nnorm_endcap_nume_tW_70to100 *1.2);
    RooRealVar n_barrel_deno_tW_70to100 ("n_barrel_deno_tW_70to100",  "n_barrel_deno_tW_70to100",  Nnorm_barrel_deno_tW_70to100,  Nnorm_barrel_deno_tW_70to100 *0.8, Nnorm_barrel_deno_tW_70to100 *1.2);
    RooRealVar n_endcap_deno_tW_70to100 ("n_endcap_deno_tW_70to100",  "n_endcap_deno_tW_70to100",  Nnorm_endcap_deno_tW_70to100,  Nnorm_endcap_deno_tW_70to100 *0.8, Nnorm_endcap_deno_tW_70to100 *1.2);
    RooRealVar n_barrel_nume_tW("n_barrel_nume_tW", "n_barrel_nume_tW", Nnorm_barrel_nume_tW, Nnorm_barrel_nume_tW*0.8, Nnorm_barrel_nume_tW*1.2);
    RooRealVar n_endcap_nume_tW("n_endcap_nume_tW", "n_endcap_nume_tW", Nnorm_endcap_nume_tW, Nnorm_endcap_nume_tW*0.8, Nnorm_endcap_nume_tW*1.2);
    RooRealVar n_barrel_deno_tW("n_barrel_deno_tW", "n_barrel_deno_tW", Nnorm_barrel_deno_tW, Nnorm_barrel_deno_tW*0.8, Nnorm_barrel_deno_tW*1.2);
    RooRealVar n_endcap_deno_tW("n_endcap_deno_tW", "n_endcap_deno_tW", Nnorm_endcap_deno_tW, Nnorm_endcap_deno_tW*0.8, Nnorm_endcap_deno_tW*1.2);

    RooRealVar n_barrel_nume_tbarW_50to70  ("n_barrel_nume_tbarW_50to70",   "n_barrel_nume_tbarW_50to70",   Nnorm_barrel_nume_tbarW_50to70 ,  Nnorm_barrel_nume_tbarW_50to70  *0.8, Nnorm_barrel_nume_tbarW_50to70  *1.2);
    RooRealVar n_endcap_nume_tbarW_50to70  ("n_endcap_nume_tbarW_50to70",   "n_endcap_nume_tbarW_50to70",   Nnorm_endcap_nume_tbarW_50to70 ,  Nnorm_endcap_nume_tbarW_50to70  *0.8, Nnorm_endcap_nume_tbarW_50to70  *1.2);
    RooRealVar n_barrel_deno_tbarW_50to70  ("n_barrel_deno_tbarW_50to70",   "n_barrel_deno_tbarW_50to70",   Nnorm_barrel_deno_tbarW_50to70 ,  Nnorm_barrel_deno_tbarW_50to70  *0.8, Nnorm_barrel_deno_tbarW_50to70  *1.2);
    RooRealVar n_endcap_deno_tbarW_50to70  ("n_endcap_deno_tbarW_50to70",   "n_endcap_deno_tbarW_50to70",   Nnorm_endcap_deno_tbarW_50to70 ,  Nnorm_endcap_deno_tbarW_50to70  *0.8, Nnorm_endcap_deno_tbarW_50to70  *1.2);
    RooRealVar n_barrel_nume_tbarW_70to100 ("n_barrel_nume_tbarW_70to100",  "n_barrel_nume_tbarW_70to100",  Nnorm_barrel_nume_tbarW_70to100,  Nnorm_barrel_nume_tbarW_70to100 *0.8, Nnorm_barrel_nume_tbarW_70to100 *1.2);
    RooRealVar n_endcap_nume_tbarW_70to100 ("n_endcap_nume_tbarW_70to100",  "n_endcap_nume_tbarW_70to100",  Nnorm_endcap_nume_tbarW_70to100,  Nnorm_endcap_nume_tbarW_70to100 *0.8, Nnorm_endcap_nume_tbarW_70to100 *1.2);
    RooRealVar n_barrel_deno_tbarW_70to100 ("n_barrel_deno_tbarW_70to100",  "n_barrel_deno_tbarW_70to100",  Nnorm_barrel_deno_tbarW_70to100,  Nnorm_barrel_deno_tbarW_70to100 *0.8, Nnorm_barrel_deno_tbarW_70to100 *1.2);
    RooRealVar n_endcap_deno_tbarW_70to100 ("n_endcap_deno_tbarW_70to100",  "n_endcap_deno_tbarW_70to100",  Nnorm_endcap_deno_tbarW_70to100,  Nnorm_endcap_deno_tbarW_70to100 *0.8, Nnorm_endcap_deno_tbarW_70to100 *1.2);
    RooRealVar n_barrel_nume_tbarW("n_barrel_nume_tbarW", "n_barrel_nume_tbarW", Nnorm_barrel_nume_tbarW, Nnorm_barrel_nume_tbarW*0.8, Nnorm_barrel_nume_tbarW*1.2);
    RooRealVar n_endcap_nume_tbarW("n_endcap_nume_tbarW", "n_endcap_nume_tbarW", Nnorm_endcap_nume_tbarW, Nnorm_endcap_nume_tbarW*0.8, Nnorm_endcap_nume_tbarW*1.2);
    RooRealVar n_barrel_deno_tbarW("n_barrel_deno_tbarW", "n_barrel_deno_tbarW", Nnorm_barrel_deno_tbarW, Nnorm_barrel_deno_tbarW*0.8, Nnorm_barrel_deno_tbarW*1.2);
    RooRealVar n_endcap_deno_tbarW("n_endcap_deno_tbarW", "n_endcap_deno_tbarW", Nnorm_endcap_deno_tbarW, Nnorm_endcap_deno_tbarW*0.8, Nnorm_endcap_deno_tbarW*1.2);

    RooRealVar n_barrel_nume_WW_50to70  ("n_barrel_nume_WW_50to70",   "n_barrel_nume_WW_50to70",   Nnorm_barrel_nume_WW_50to70 ,  Nnorm_barrel_nume_WW_50to70  *0.8, Nnorm_barrel_nume_WW_50to70  *1.2);
    RooRealVar n_endcap_nume_WW_50to70  ("n_endcap_nume_WW_50to70",   "n_endcap_nume_WW_50to70",   Nnorm_endcap_nume_WW_50to70 ,  Nnorm_endcap_nume_WW_50to70  *0.8, Nnorm_endcap_nume_WW_50to70  *1.2);
    RooRealVar n_barrel_deno_WW_50to70  ("n_barrel_deno_WW_50to70",   "n_barrel_deno_WW_50to70",   Nnorm_barrel_deno_WW_50to70 ,  Nnorm_barrel_deno_WW_50to70  *0.8, Nnorm_barrel_deno_WW_50to70  *1.2);
    RooRealVar n_endcap_deno_WW_50to70  ("n_endcap_deno_WW_50to70",   "n_endcap_deno_WW_50to70",   Nnorm_endcap_deno_WW_50to70 ,  Nnorm_endcap_deno_WW_50to70  *0.8, Nnorm_endcap_deno_WW_50to70  *1.2);
    RooRealVar n_barrel_nume_WW_70to100 ("n_barrel_nume_WW_70to100",  "n_barrel_nume_WW_70to100",  Nnorm_barrel_nume_WW_70to100,  Nnorm_barrel_nume_WW_70to100 *0.8, Nnorm_barrel_nume_WW_70to100 *1.2);
    RooRealVar n_endcap_nume_WW_70to100 ("n_endcap_nume_WW_70to100",  "n_endcap_nume_WW_70to100",  Nnorm_endcap_nume_WW_70to100,  Nnorm_endcap_nume_WW_70to100 *0.8, Nnorm_endcap_nume_WW_70to100 *1.2);
    RooRealVar n_barrel_deno_WW_70to100 ("n_barrel_deno_WW_70to100",  "n_barrel_deno_WW_70to100",  Nnorm_barrel_deno_WW_70to100,  Nnorm_barrel_deno_WW_70to100 *0.8, Nnorm_barrel_deno_WW_70to100 *1.2);
    RooRealVar n_endcap_deno_WW_70to100 ("n_endcap_deno_WW_70to100",  "n_endcap_deno_WW_70to100",  Nnorm_endcap_deno_WW_70to100,  Nnorm_endcap_deno_WW_70to100 *0.8, Nnorm_endcap_deno_WW_70to100 *1.2);
    RooRealVar n_barrel_nume_WW("n_barrel_nume_WW", "n_barrel_nume_WW", Nnorm_barrel_nume_WW, Nnorm_barrel_nume_WW*0.8, Nnorm_barrel_nume_WW*1.2);
    RooRealVar n_endcap_nume_WW("n_endcap_nume_WW", "n_endcap_nume_WW", Nnorm_endcap_nume_WW, Nnorm_endcap_nume_WW*0.8, Nnorm_endcap_nume_WW*1.2);
    RooRealVar n_barrel_deno_WW("n_barrel_deno_WW", "n_barrel_deno_WW", Nnorm_barrel_deno_WW, Nnorm_barrel_deno_WW*0.8, Nnorm_barrel_deno_WW*1.2);
    RooRealVar n_endcap_deno_WW("n_endcap_deno_WW", "n_endcap_deno_WW", Nnorm_endcap_deno_WW, Nnorm_endcap_deno_WW*0.8, Nnorm_endcap_deno_WW*1.2);

    RooRealVar n_barrel_nume_WZ_50to70  ("n_barrel_nume_WZ_50to70",   "n_barrel_nume_WZ_50to70",   Nnorm_barrel_nume_WZ_50to70 ,  Nnorm_barrel_nume_WZ_50to70  *0.8, Nnorm_barrel_nume_WZ_50to70  *1.2);
    RooRealVar n_endcap_nume_WZ_50to70  ("n_endcap_nume_WZ_50to70",   "n_endcap_nume_WZ_50to70",   Nnorm_endcap_nume_WZ_50to70 ,  Nnorm_endcap_nume_WZ_50to70  *0.8, Nnorm_endcap_nume_WZ_50to70  *1.2);
    RooRealVar n_barrel_deno_WZ_50to70  ("n_barrel_deno_WZ_50to70",   "n_barrel_deno_WZ_50to70",   Nnorm_barrel_deno_WZ_50to70 ,  Nnorm_barrel_deno_WZ_50to70  *0.8, Nnorm_barrel_deno_WZ_50to70  *1.2);
    RooRealVar n_endcap_deno_WZ_50to70  ("n_endcap_deno_WZ_50to70",   "n_endcap_deno_WZ_50to70",   Nnorm_endcap_deno_WZ_50to70 ,  Nnorm_endcap_deno_WZ_50to70  *0.8, Nnorm_endcap_deno_WZ_50to70  *1.2);
    RooRealVar n_barrel_nume_WZ_70to100 ("n_barrel_nume_WZ_70to100",  "n_barrel_nume_WZ_70to100",  Nnorm_barrel_nume_WZ_70to100,  Nnorm_barrel_nume_WZ_70to100 *0.8, Nnorm_barrel_nume_WZ_70to100 *1.2);
    RooRealVar n_endcap_nume_WZ_70to100 ("n_endcap_nume_WZ_70to100",  "n_endcap_nume_WZ_70to100",  Nnorm_endcap_nume_WZ_70to100,  Nnorm_endcap_nume_WZ_70to100 *0.8, Nnorm_endcap_nume_WZ_70to100 *1.2);
    RooRealVar n_barrel_deno_WZ_70to100 ("n_barrel_deno_WZ_70to100",  "n_barrel_deno_WZ_70to100",  Nnorm_barrel_deno_WZ_70to100,  Nnorm_barrel_deno_WZ_70to100 *0.8, Nnorm_barrel_deno_WZ_70to100 *1.2);
    RooRealVar n_endcap_deno_WZ_70to100 ("n_endcap_deno_WZ_70to100",  "n_endcap_deno_WZ_70to100",  Nnorm_endcap_deno_WZ_70to100,  Nnorm_endcap_deno_WZ_70to100 *0.8, Nnorm_endcap_deno_WZ_70to100 *1.2);
    RooRealVar n_barrel_nume_WZ("n_barrel_nume_WZ", "n_barrel_nume_WZ", Nnorm_barrel_nume_WZ, Nnorm_barrel_nume_WZ*0.8, Nnorm_barrel_nume_WZ*1.2);
    RooRealVar n_endcap_nume_WZ("n_endcap_nume_WZ", "n_endcap_nume_WZ", Nnorm_endcap_nume_WZ, Nnorm_endcap_nume_WZ*0.8, Nnorm_endcap_nume_WZ*1.2);
    RooRealVar n_barrel_deno_WZ("n_barrel_deno_WZ", "n_barrel_deno_WZ", Nnorm_barrel_deno_WZ, Nnorm_barrel_deno_WZ*0.8, Nnorm_barrel_deno_WZ*1.2);
    RooRealVar n_endcap_deno_WZ("n_endcap_deno_WZ", "n_endcap_deno_WZ", Nnorm_endcap_deno_WZ, Nnorm_endcap_deno_WZ*0.8, Nnorm_endcap_deno_WZ*1.2);

    RooRealVar n_barrel_nume_ZZ_50to70  ("n_barrel_nume_ZZ_50to70",   "n_barrel_nume_ZZ_50to70",   Nnorm_barrel_nume_ZZ_50to70 ,  Nnorm_barrel_nume_ZZ_50to70  *0.8, Nnorm_barrel_nume_ZZ_50to70  *1.2);
    RooRealVar n_endcap_nume_ZZ_50to70  ("n_endcap_nume_ZZ_50to70",   "n_endcap_nume_ZZ_50to70",   Nnorm_endcap_nume_ZZ_50to70 ,  Nnorm_endcap_nume_ZZ_50to70  *0.8, Nnorm_endcap_nume_ZZ_50to70  *1.2);
    RooRealVar n_barrel_deno_ZZ_50to70  ("n_barrel_deno_ZZ_50to70",   "n_barrel_deno_ZZ_50to70",   Nnorm_barrel_deno_ZZ_50to70 ,  Nnorm_barrel_deno_ZZ_50to70  *0.8, Nnorm_barrel_deno_ZZ_50to70  *1.2);
    RooRealVar n_endcap_deno_ZZ_50to70  ("n_endcap_deno_ZZ_50to70",   "n_endcap_deno_ZZ_50to70",   Nnorm_endcap_deno_ZZ_50to70 ,  Nnorm_endcap_deno_ZZ_50to70  *0.8, Nnorm_endcap_deno_ZZ_50to70  *1.2);
    RooRealVar n_barrel_nume_ZZ_70to100 ("n_barrel_nume_ZZ_70to100",  "n_barrel_nume_ZZ_70to100",  Nnorm_barrel_nume_ZZ_70to100,  Nnorm_barrel_nume_ZZ_70to100 *0.8, Nnorm_barrel_nume_ZZ_70to100 *1.2);
    RooRealVar n_endcap_nume_ZZ_70to100 ("n_endcap_nume_ZZ_70to100",  "n_endcap_nume_ZZ_70to100",  Nnorm_endcap_nume_ZZ_70to100,  Nnorm_endcap_nume_ZZ_70to100 *0.8, Nnorm_endcap_nume_ZZ_70to100 *1.2);
    RooRealVar n_barrel_deno_ZZ_70to100 ("n_barrel_deno_ZZ_70to100",  "n_barrel_deno_ZZ_70to100",  Nnorm_barrel_deno_ZZ_70to100,  Nnorm_barrel_deno_ZZ_70to100 *0.8, Nnorm_barrel_deno_ZZ_70to100 *1.2);
    RooRealVar n_endcap_deno_ZZ_70to100 ("n_endcap_deno_ZZ_70to100",  "n_endcap_deno_ZZ_70to100",  Nnorm_endcap_deno_ZZ_70to100,  Nnorm_endcap_deno_ZZ_70to100 *0.8, Nnorm_endcap_deno_ZZ_70to100 *1.2);
    RooRealVar n_barrel_nume_ZZ("n_barrel_nume_ZZ", "n_barrel_nume_ZZ", Nnorm_barrel_nume_ZZ, Nnorm_barrel_nume_ZZ*0.8, Nnorm_barrel_nume_ZZ*1.2);
    RooRealVar n_endcap_nume_ZZ("n_endcap_nume_ZZ", "n_endcap_nume_ZZ", Nnorm_endcap_nume_ZZ, Nnorm_endcap_nume_ZZ*0.8, Nnorm_endcap_nume_ZZ*1.2);
    RooRealVar n_barrel_deno_ZZ("n_barrel_deno_ZZ", "n_barrel_deno_ZZ", Nnorm_barrel_deno_ZZ, Nnorm_barrel_deno_ZZ*0.8, Nnorm_barrel_deno_ZZ*1.2);
    RooRealVar n_endcap_deno_ZZ("n_endcap_deno_ZZ", "n_endcap_deno_ZZ", Nnorm_endcap_deno_ZZ, Nnorm_endcap_deno_ZZ*0.8, Nnorm_endcap_deno_ZZ*1.2);

    RooRealVar n_barrel_nume_QCD_50to70  ("n_barrel_nume_QCD_50to70",   "n_barrel_nume_QCD_50to70",   Nnorm_barrel_nume_QCD_50to70 ,  Nnorm_barrel_nume_QCD_50to70  *0.5, Nnorm_barrel_nume_QCD_50to70  *1.5);
    RooRealVar n_endcap_nume_QCD_50to70  ("n_endcap_nume_QCD_50to70",   "n_endcap_nume_QCD_50to70",   Nnorm_endcap_nume_QCD_50to70 ,  Nnorm_endcap_nume_QCD_50to70  *0.5, Nnorm_endcap_nume_QCD_50to70  *1.5);
    RooRealVar n_barrel_deno_QCD_50to70  ("n_barrel_deno_QCD_50to70",   "n_barrel_deno_QCD_50to70",   Nnorm_barrel_deno_QCD_50to70 ,  Nnorm_barrel_deno_QCD_50to70  *0.5, Nnorm_barrel_deno_QCD_50to70  *1.5);
    RooRealVar n_endcap_deno_QCD_50to70  ("n_endcap_deno_QCD_50to70",   "n_endcap_deno_QCD_50to70",   Nnorm_endcap_deno_QCD_50to70 ,  Nnorm_endcap_deno_QCD_50to70  *0.5, Nnorm_endcap_deno_QCD_50to70  *1.5);
    RooRealVar n_barrel_nume_QCD_70to100 ("n_barrel_nume_QCD_70to100",  "n_barrel_nume_QCD_70to100",  Nnorm_barrel_nume_QCD_70to100,  Nnorm_barrel_nume_QCD_70to100 *0.5, Nnorm_barrel_nume_QCD_70to100 *1.5);
    RooRealVar n_endcap_nume_QCD_70to100 ("n_endcap_nume_QCD_70to100",  "n_endcap_nume_QCD_70to100",  Nnorm_endcap_nume_QCD_70to100,  Nnorm_endcap_nume_QCD_70to100 *0.5, Nnorm_endcap_nume_QCD_70to100 *1.5);
    RooRealVar n_barrel_deno_QCD_70to100 ("n_barrel_deno_QCD_70to100",  "n_barrel_deno_QCD_70to100",  Nnorm_barrel_deno_QCD_70to100,  Nnorm_barrel_deno_QCD_70to100 *0.5, Nnorm_barrel_deno_QCD_70to100 *1.5);
    RooRealVar n_endcap_deno_QCD_70to100 ("n_endcap_deno_QCD_70to100",  "n_endcap_deno_QCD_70to100",  Nnorm_endcap_deno_QCD_70to100,  Nnorm_endcap_deno_QCD_70to100 *0.5, Nnorm_endcap_deno_QCD_70to100 *1.5);
    RooRealVar n_barrel_nume_QCD("n_barrel_nume_QCD", "n_barrel_nume_QCD", Nnorm_barrel_nume_QCD, Nnorm_barrel_nume_QCD*0.5, Nnorm_barrel_nume_QCD*1.5);
    RooRealVar n_endcap_nume_QCD("n_endcap_nume_QCD", "n_endcap_nume_QCD", Nnorm_endcap_nume_QCD, Nnorm_endcap_nume_QCD*0.5, Nnorm_endcap_nume_QCD*1.5);
    RooRealVar n_barrel_deno_QCD("n_barrel_deno_QCD", "n_barrel_deno_QCD", Nnorm_barrel_deno_QCD, Nnorm_barrel_deno_QCD*0.5, Nnorm_barrel_deno_QCD*1.5);
    RooRealVar n_endcap_deno_QCD("n_endcap_deno_QCD", "n_endcap_deno_QCD", Nnorm_endcap_deno_QCD, Nnorm_endcap_deno_QCD*0.5, Nnorm_endcap_deno_QCD*1.5);

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

    RooAddPdf model_barrel_nume("model_barrel_nume", "model_barrel_nume",
                                         RooArgList(*pdf_barrel_nume_QCD,   *pdf_barrel_nume_WJets, *pdf_barrel_nume_DY,
                                                    *pdf_barrel_nume_ttbar, *pdf_barrel_nume_tW,    *pdf_barrel_nume_tbarW,
                                                    *pdf_barrel_nume_WW,    *pdf_barrel_nume_WZ,    *pdf_barrel_nume_ZZ),
                                         RooArgList(n_barrel_nume_QCD,   n_barrel_nume_WJets, n_barrel_nume_DY,
                                                    n_barrel_nume_ttbar, n_barrel_nume_tW, n_barrel_nume_tbarW,
                                                    n_barrel_nume_WW,    n_barrel_nume_WZ, n_barrel_nume_ZZ));

    RooAddPdf model_endcap_nume("model_endcap_nume", "model_endcap_nume",
                                         RooArgList(*pdf_endcap_nume_QCD,   *pdf_endcap_nume_WJets, *pdf_endcap_nume_DY,
                                                    *pdf_endcap_nume_ttbar, *pdf_endcap_nume_tW,    *pdf_endcap_nume_tbarW,
                                                    *pdf_endcap_nume_WW,    *pdf_endcap_nume_WZ,    *pdf_endcap_nume_ZZ),
                                         RooArgList(n_endcap_nume_QCD,   n_endcap_nume_WJets, n_endcap_nume_DY,
                                                    n_endcap_nume_ttbar, n_endcap_nume_tW,    n_endcap_nume_tbarW,
                                                    n_endcap_nume_WW,    n_endcap_nume_WZ,    n_endcap_nume_ZZ));

    RooAddPdf model_barrel_deno("model_barrel_deno", "model_barrel_deno",
                                         RooArgList(*pdf_barrel_deno_QCD,   *pdf_barrel_deno_WJets, *pdf_barrel_deno_DY,
                                                    *pdf_barrel_deno_ttbar, *pdf_barrel_deno_tW,    *pdf_barrel_deno_tbarW,
                                                    *pdf_barrel_deno_WW,    *pdf_barrel_deno_WZ,    *pdf_barrel_deno_ZZ),
                                         RooArgList(n_barrel_deno_QCD,   n_barrel_deno_WJets, n_barrel_deno_DY,
                                                    n_barrel_deno_ttbar, n_barrel_deno_tW,    n_barrel_deno_tbarW,
                                                    n_barrel_deno_WW,    n_barrel_deno_WZ,    n_barrel_deno_ZZ));

    RooAddPdf model_endcap_deno("model_endcap_deno", "model_endcap_deno",
                                         RooArgList(*pdf_endcap_deno_QCD,   *pdf_endcap_deno_WJets, *pdf_endcap_deno_DY,
                                                    *pdf_endcap_deno_ttbar, *pdf_endcap_deno_tW,    *pdf_endcap_deno_tbarW,
                                                    *pdf_endcap_deno_WW,    *pdf_endcap_deno_WZ,    *pdf_endcap_deno_ZZ),
                                         RooArgList(n_endcap_deno_QCD,   n_endcap_deno_WJets, n_endcap_deno_DY,
                                                    n_endcap_deno_ttbar, n_endcap_deno_tW,    n_endcap_deno_tbarW,
                                                    n_endcap_deno_WW,    n_endcap_deno_WZ,    n_endcap_deno_ZZ));

    // Fitting
    RooFitResult* fit_barrel_nume_50to70   = model_barrel_nume_50to70  .fitTo(*rh_barrel_nume_data_50to70  , Save());
    RooFitResult* fit_endcap_nume_50to70   = model_endcap_nume_50to70  .fitTo(*rh_endcap_nume_data_50to70  , Save());
    RooFitResult* fit_barrel_deno_50to70   = model_barrel_deno_50to70  .fitTo(*rh_barrel_deno_data_50to70  , Save());
    RooFitResult* fit_endcap_deno_50to70   = model_endcap_deno_50to70  .fitTo(*rh_endcap_deno_data_50to70  , Save());
    RooFitResult* fit_barrel_nume_70to100  = model_barrel_nume_70to100 .fitTo(*rh_barrel_nume_data_70to100 , Save());
    RooFitResult* fit_endcap_nume_70to100  = model_endcap_nume_70to100 .fitTo(*rh_endcap_nume_data_70to100 , Save());
    RooFitResult* fit_barrel_deno_70to100  = model_barrel_deno_70to100 .fitTo(*rh_barrel_deno_data_70to100 , Save());
    RooFitResult* fit_endcap_deno_70to100  = model_endcap_deno_70to100 .fitTo(*rh_endcap_deno_data_70to100 , Save());
    RooFitResult* fit_barrel_nume = model_barrel_nume.fitTo(*rh_barrel_nume_data, Save());
    RooFitResult* fit_endcap_nume = model_endcap_nume.fitTo(*rh_endcap_nume_data, Save());
    RooFitResult* fit_barrel_deno = model_barrel_deno.fitTo(*rh_barrel_deno_data, Save());
    RooFitResult* fit_endcap_deno = model_endcap_deno.fitTo(*rh_endcap_deno_data, Save());

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
    frame_barrel_nume_50to70->GetYaxis()->SetRangeUser(1e0, 2e8);
    frame_barrel_nume_50to70->Draw();
    fit_barrel_nume_50to70->Print();

    // Legend
    TLegend *legend = new TLegend(0.65, 0.8, 0.95, 0.97);
    legend->SetFillColor(kWhite);
//    legend->SetLineColor(kWhite);
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
    frame_endcap_nume_50to70->GetYaxis()->SetRangeUser(1e0, 2e8);
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
    frame_barrel_deno_50to70->GetYaxis()->SetRangeUser(1e0, 2e8);
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
    frame_endcap_deno_50to70->GetYaxis()->SetRangeUser(1e0, 2e7);
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
    frame_barrel_nume_70to100->GetYaxis()->SetRangeUser(1e0, 2e8);
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
    frame_endcap_nume_70to100->GetYaxis()->SetRangeUser(1e0, 2e8);
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
    frame_barrel_deno_70to100->GetYaxis()->SetRangeUser(1e0, 2e7);
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
    frame_endcap_deno_70to100->GetYaxis()->SetRangeUser(1e0, 2e7);
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
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW,pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,"
                                                                             "pdf_barrel_nume_DY,pdf_barrel_nume_WJets,pdf_barrel_nume_QCD"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW,pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,"
                                                                             "pdf_barrel_nume_DY,pdf_barrel_nume_WJets"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW,pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar,"
                                                                             "pdf_barrel_nume_DY"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW,pdf_barrel_nume_tbarW,pdf_barrel_nume_ttbar"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW,pdf_barrel_nume_tbarW"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW,"
                                                                             "pdf_barrel_nume_tW"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ,pdf_barrel_nume_WW"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ,pdf_barrel_nume_WZ"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_nume.plotOn(frame_barrel_nume, Components("pdf_barrel_nume_ZZ"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_nume_data->plotOn(frame_barrel_nume, DataError(RooAbsData::SumW2));
    frame_barrel_nume->GetYaxis()->SetRangeUser(1e0, 5e7);
    frame_barrel_nume->Draw();
    fit_barrel_nume->Print();

    legend->Draw();

    frame_barrel_nume->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_nume->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_nume->GetXaxis()->SetLabelSize(0);

    // Bottom pad
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
    Double_t N_barrel_nume_MC = h_barrel_nume_MC_fit->Integral();
    h_barrel_nume_MC_fit->Scale(N_barrel_nume_data/N_barrel_nume_MC); // Why would I wanna do that???
    cout << "\nData integral: " << N_barrel_nume_data << endl;
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
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_nume = model_barrel_nume.createChi2(*rh_barrel_nume_data);
    cout << "chi2: " << chi2_barrel_nume->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_nume->getVal() / ((Double_t)h_barrel_data_nume->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING NUMERATOR ENDCAP 100to500
    cout << "\n----- NUMERATOR ENDCAP 100to500 -----" << endl;
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
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                             "pdf_endcap_nume_tW,pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,"
                                                                             "pdf_endcap_nume_DY,pdf_endcap_nume_WJets,pdf_endcap_nume_QCD"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                             "pdf_endcap_nume_tW,pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,"
                                                                             "pdf_endcap_nume_DY,pdf_endcap_nume_WJets"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                             "pdf_endcap_nume_tW,pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar,"
                                                                             "pdf_endcap_nume_DY"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                             "pdf_endcap_nume_tW,pdf_endcap_nume_tbarW,pdf_endcap_nume_ttbar"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                             "pdf_endcap_nume_tW,pdf_endcap_nume_tbarW"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW,"
                                                                            "pdf_endcap_nume_tW"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ,pdf_endcap_nume_WW"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ,pdf_endcap_nume_WZ"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_nume.plotOn(frame_endcap_nume, Components("pdf_endcap_nume_ZZ"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_nume_data->plotOn(frame_endcap_nume, DataError(RooAbsData::SumW2));
    frame_endcap_nume->GetYaxis()->SetRangeUser(1e0, 5e7);
    frame_endcap_nume->Draw();
    fit_endcap_nume->Print();

    // Legend
    legend->Draw();

    frame_endcap_nume->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_nume->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_nume->GetXaxis()->SetLabelSize(0);

    // Bottom pad
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
    cout << "\nData integral: " << N_endcap_nume_data << endl;
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
    cout << "Normalized chi2: " << chi2_endcap_nume->getVal() / ((Double_t)h_endcap_data_nume->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR BARREL 100to500
    cout << "\n----- DENOMINATOR BARREL 100to500 -----" << endl;
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
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW,pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,"
                                                                             "pdf_barrel_deno_DY,pdf_barrel_deno_WJets,pdf_barrel_deno_QCD"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW,pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,"
                                                                             "pdf_barrel_deno_DY,pdf_barrel_deno_WJets"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW,pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar,"
                                                                             "pdf_barrel_deno_DY"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW,pdf_barrel_deno_tbarW,pdf_barrel_deno_ttbar"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW,pdf_barrel_deno_tbarW"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW,"
                                                                             "pdf_barrel_deno_tW"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ,pdf_barrel_deno_WW"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ,pdf_barrel_deno_WZ"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_barrel_deno.plotOn(frame_barrel_deno, Components("pdf_barrel_deno_ZZ"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_barrel_deno_data->plotOn(frame_barrel_deno, DataError(RooAbsData::SumW2));
    frame_barrel_deno->GetYaxis()->SetRangeUser(1e0, 2e7);
    frame_barrel_deno->Draw();
    fit_barrel_deno->Print();

    // Legend
    legend->Draw();

    frame_barrel_deno->GetYaxis()->SetTitle("Number of entries");
    frame_barrel_deno->GetYaxis()->SetTitleOffset(1.5);
    frame_barrel_deno->GetXaxis()->SetLabelSize(0);

    // Bottom pad
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
    cout << "\nData integral: " << N_barrel_deno_data << endl;
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
    h_line_deno->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_barrel_deno = model_barrel_deno.createChi2(*rh_barrel_deno_data);
    cout << "chi2: " << chi2_barrel_deno->getVal() << endl;
    cout << "Normalized chi2: " << chi2_barrel_deno->getVal() / ((Double_t)h_barrel_data_deno->GetNbinsX()) << "\n\n" << endl;

    /// DRAWING DENOMINATOR ENDCAP 100to500
    cout << "\n----- DENOMINATOR ENDCAP 100to500 -----" << endl;
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
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW,pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,"
                                                                             "pdf_endcap_deno_DY,pdf_endcap_deno_WJets,pdf_endcap_deno_QCD"),
                                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW,pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,"
                                                                             "pdf_endcap_deno_DY,pdf_endcap_deno_WJets"),
                                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW,pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar,"
                                                                             "pdf_endcap_deno_DY"),
                                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW,pdf_endcap_deno_tbarW,pdf_endcap_deno_ttbar"),
                                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW,pdf_endcap_deno_tbarW"),
                                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW,"
                                                                             "pdf_endcap_deno_tW"),
                                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ,pdf_endcap_deno_WW"),
                                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ,pdf_endcap_deno_WZ"),
                                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_endcap_deno.plotOn(frame_endcap_deno, Components("pdf_endcap_deno_ZZ"),
                                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_endcap_deno_data->plotOn(frame_endcap_deno, DataError(RooAbsData::SumW2));
    frame_endcap_deno->GetYaxis()->SetRangeUser(1e0, 2e7);
    frame_endcap_deno->Draw();
    fit_endcap_deno->Print();

    // Legend
    legend->Draw();

    frame_endcap_deno->GetYaxis()->SetTitle("Number of entries");
    frame_endcap_deno->GetYaxis()->SetTitleOffset(1.5);
    frame_endcap_deno->GetXaxis()->SetLabelSize(0);

    // Bottom pad
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
    cout << "\nData integral: " << N_endcap_deno_data << endl;
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

    cout << "Starting integral values for QCD (check for final values in fit results):" << endl;
    cout << "Numerator barrel pT 50-70  : " << N_barrel_nume_QCD_50to70  << endl;
    cout << "Numerator endcap pT 50-70  : " << N_endcap_nume_QCD_50to70  << endl;
    cout << "Numerator barrel pT 70-100 : " << N_barrel_nume_QCD_70to100 << endl;
    cout << "Numerator endcap pT 70-100 : " << N_endcap_nume_QCD_70to100 << endl;
    cout << "Numerator barrel pT 100-500: " << N_barrel_nume_QCD         << endl;
    cout << "Numerator endcap pT 100-500: " << N_endcap_nume_QCD         << endl;
    cout << "Denominator barrel pT 50-70  : " << N_barrel_deno_QCD_50to70  << endl;
    cout << "Denominator endcap pT 50-70  : " << N_endcap_deno_QCD_50to70  << endl;
    cout << "Denominator barrel pT 70-100 : " << N_barrel_deno_QCD_70to100 << endl;
    cout << "Denominator endcap pT 70-100 : " << N_endcap_deno_QCD_70to100 << endl;
    cout << "Denominator barrel pT 100-500: " << N_barrel_deno_QCD         << endl;
    cout << "Denominator endcap pT 100-500: " << N_endcap_deno_QCD         << endl;

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

} // End of Mu_TFit()


/// ################################################################################## ///
/// ---------------------------- Muon W+Jets estimation ------------------------------ ///
/// ################################################################################## ///
void Mu_WJETSest_Tfit(Int_t type)
{
    FileMgr Mgr;

    TH1D *h_MC_mass[_EndOf_Data_Special], *h_MC_mass_SS[_EndOf_Data_Special], *h_MC_mass_dijet[_EndOf_Data_Special], *h_MC_mass_SS_dijet[_EndOf_Data_Special],
         *h_data_mass, *h_data_mass_SS, *h_data_mass_dijet, *h_data_mass_SS_dijet;

    TFile *file = new TFile("/media/sf_DATA/FR/Muon/WJETest_Mu.root", "READ");
    TFile *file2 = new TFile("/media/sf_DATA/FR/Muon/QCDest_Mu.root", "READ");

// ############################# SETUP ################################# //
//----------------------------- MC bkg ------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        if (type == 1)
        {
            file->GetObject("h_mass_forFit_"+Mgr.Procname[pr1],   h_MC_mass[pr1]);
            file->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr1],   h_MC_mass_SS[pr1]);
            file2->GetObject("h_mass_forFit_"+Mgr.Procname[pr1],   h_MC_mass_dijet[pr1]);
            file2->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr1],   h_MC_mass_SS_dijet[pr1]);
        }
        else if (type == 2) // Fitting with regular binning
        {
            file->GetObject("h_mass_"+Mgr.Procname[pr1],   h_MC_mass[pr1]);
            file->GetObject("h_mass_SS_"+Mgr.Procname[pr1],   h_MC_mass_SS[pr1]);
            file2->GetObject("h_mass_"+Mgr.Procname[pr1],   h_MC_mass_dijet[pr1]);
            file2->GetObject("h_mass_SS_"+Mgr.Procname[pr1],   h_MC_mass_SS_dijet[pr1]);
        }
        else
        {
            cout << "Wrong 'type'!\nSelect '1' to fit histograms with 5 GeV bins or '2' to fit with regular binning" << endl;\
            return;
        }

        removeNegativeBins(h_MC_mass[pr1]);
        removeNegativeBins(h_MC_mass_dijet[pr1]);
        removeNegativeBins(h_MC_mass_SS[pr1]);
        removeNegativeBins(h_MC_mass_SS_dijet[pr1]);
        h_MC_mass[pr1]->SetDirectory(0);
        h_MC_mass_dijet[pr1]->SetDirectory(0);
        h_MC_mass_SS[pr1]->SetDirectory(0);
        h_MC_mass_SS_dijet[pr1]->SetDirectory(0);

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
    h_MC_mass[_ttbar]->Add(h_MC_mass[_ttbar_700to1000]);
    h_MC_mass[_ttbar]->Add(h_MC_mass[_ttbar_1000toInf]);
    h_MC_mass[_WJets]->Add(h_MC_mass[_WJets_ext2v5]);
    h_MC_mass_SS[_ttbar]->Add(h_MC_mass_SS[_ttbar_700to1000]);
    h_MC_mass_SS[_ttbar]->Add(h_MC_mass_SS[_ttbar_1000toInf]);
    h_MC_mass_SS[_WJets]->Add(h_MC_mass_SS[_WJets_ext2v5]);
    h_MC_mass_dijet[_ttbar]->Add(h_MC_mass_dijet[_ttbar_700to1000]);
    h_MC_mass_dijet[_ttbar]->Add(h_MC_mass_dijet[_ttbar_1000toInf]);
    h_MC_mass_dijet[_WJets]->Add(h_MC_mass_dijet[_WJets_ext2v5]);
    h_MC_mass_SS_dijet[_ttbar]->Add(h_MC_mass_SS_dijet[_ttbar_700to1000]);
    h_MC_mass_SS_dijet[_ttbar]->Add(h_MC_mass_SS_dijet[_ttbar_1000toInf]);
    h_MC_mass_SS_dijet[_WJets]->Add(h_MC_mass_SS_dijet[_WJets_ext2v5]);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        if (type == 1)
        {
            file->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_MC_mass[pr]);
            file->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_MC_mass_SS[pr]);
            file2->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_MC_mass_dijet[pr]);
            file2->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_MC_mass_SS_dijet[pr]);
        }
        else if (type == 2)
        {
            file->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass[pr]);
            file->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_MC_mass_SS[pr]);
            file2->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass_dijet[pr]);
            file2->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_MC_mass_SS_dijet[pr]);
        }

        removeNegativeBins(h_MC_mass[pr]);
        removeNegativeBins(h_MC_mass_SS[pr]);
        removeNegativeBins(h_MC_mass_dijet[pr]);
        removeNegativeBins(h_MC_mass_SS_dijet[pr]);
        h_MC_mass[pr]->SetDirectory(0);
        h_MC_mass_SS[pr]->SetDirectory(0);
        h_MC_mass_dijet[pr]->SetDirectory(0);
        h_MC_mass_SS_dijet[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_MC_mass[_DY_Full] = ((TH1D*)(h_MC_mass[pr]->Clone("h_MC_mass_DY")));
            h_MC_mass_SS[_DY_Full] = ((TH1D*)(h_MC_mass_SS[pr]->Clone("h_MC_mass_SS_DY")));
            h_MC_mass_dijet[_DY_Full] = ((TH1D*)(h_MC_mass_dijet[pr]->Clone("h_MC_mass_dijet_DY")));
            h_MC_mass_SS_dijet[_DY_Full] = ((TH1D*)(h_MC_mass_SS_dijet[pr]->Clone("h_MC_mass_SS_dijet_DY")));
            h_MC_mass[_DY_Full]->SetDirectory(0);
            h_MC_mass_SS[_DY_Full]->SetDirectory(0);
            h_MC_mass_dijet[_DY_Full]->SetDirectory(0);
            h_MC_mass_SS_dijet[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_MC_mass[_DY_Full]->Add(h_MC_mass[pr]);
            h_MC_mass_SS[_DY_Full]->Add(h_MC_mass_SS[pr]);
            h_MC_mass_dijet[_DY_Full]->Add(h_MC_mass_dijet[pr]);
            h_MC_mass_SS_dijet[_DY_Full]->Add(h_MC_mass_SS_dijet[pr]);
        }
    }

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        if (type == 1)
        {
            file->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_MC_mass[pr]);
            file->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_MC_mass_SS[pr]);
            file2->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_MC_mass_dijet[pr]);
            file2->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_MC_mass_SS_dijet[pr]);
        }
        else if (type == 2)
        {
            file->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass[pr]);
            file->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_MC_mass_SS[pr]);
            file2->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass_dijet[pr]);
            file2->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_MC_mass_SS_dijet[pr]);
        }

        removeNegativeBins(h_MC_mass[pr]);
        removeNegativeBins(h_MC_mass_SS[pr]);
        removeNegativeBins(h_MC_mass_dijet[pr]);
        removeNegativeBins(h_MC_mass_SS_dijet[pr]);
        h_MC_mass[pr]->SetDirectory(0);
        h_MC_mass_SS[pr]->SetDirectory(0);
        h_MC_mass_SS_dijet[pr]->SetDirectory(0);
        h_MC_mass_dijet[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_MC_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_MC_mass[pr]->Clone("h_MC_mass_QCD")));
            h_MC_mass_SS[_QCDMuEnriched_Full] = ((TH1D*)(h_MC_mass_SS[pr]->Clone("h_MC_mass_SS_QCD")));
            h_MC_mass_dijet[_QCDMuEnriched_Full] = ((TH1D*)(h_MC_mass_dijet[pr]->Clone("h_MC_mass_dijet_QCD")));
            h_MC_mass_SS_dijet[_QCDMuEnriched_Full] = ((TH1D*)(h_MC_mass_SS_dijet[pr]->Clone("h_MC_mass_SS_dijet_QCD")));
            h_MC_mass[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MC_mass_SS[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MC_mass_dijet[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MC_mass_SS_dijet[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_MC_mass[_QCDMuEnriched_Full]->Add(h_MC_mass[pr]);
            h_MC_mass_SS[_QCDMuEnriched_Full]->Add(h_MC_mass_SS[pr]);
            h_MC_mass_dijet[_QCDMuEnriched_Full]->Add(h_MC_mass_dijet[pr]);
            h_MC_mass_SS_dijet[_QCDMuEnriched_Full]->Add(h_MC_mass_SS_dijet[pr]);
        }
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TH1D *h_temp[4];
        if (pr == _SingleMuon_B)
        {
            if (type == 1)
            {
                file->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_data_mass);
                file->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_data_mass_SS);
                file2->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_data_mass_dijet);
                file2->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_data_mass_SS_dijet);
            }
            else if (type == 2)
            {
                file->GetObject("h_mass_"+Mgr.Procname[pr], h_data_mass);
                file->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_data_mass_SS);
                file2->GetObject("h_mass_"+Mgr.Procname[pr], h_data_mass_dijet);
                file2->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_data_mass_SS_dijet);
            }
            removeNegativeBins(h_data_mass);
            removeNegativeBins(h_data_mass_SS);
            removeNegativeBins(h_data_mass_dijet);
            removeNegativeBins(h_data_mass_SS_dijet);
        }
        else
        {
            if (type == 1)
            {
                file->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_temp[0]);
                file->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_temp[1]);
                file2->GetObject("h_mass_forFit_"+Mgr.Procname[pr], h_temp[2]);
                file2->GetObject("h_mass_SS_forFit_"+Mgr.Procname[pr], h_temp[3]);
            }
            else if (type == 2)
            {
                file->GetObject("h_mass_"+Mgr.Procname[pr], h_temp[0]);
                file->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_temp[1]);
                file2->GetObject("h_mass_"+Mgr.Procname[pr], h_temp[2]);
                file2->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_temp[3]);
            }
            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            h_data_mass->Add(h_temp[0]);
            h_data_mass_SS->Add(h_temp[1]);
            h_data_mass_dijet->Add(h_temp[2]);
            h_data_mass_SS_dijet->Add(h_temp[3]);
        }
    }
    h_data_mass->SetDirectory(0);
    h_data_mass_SS->SetDirectory(0);
    h_data_mass_dijet->SetDirectory(0);
    h_data_mass_SS_dijet->SetDirectory(0);

// ######################## MODEL BUILDING ##########################

    // Making data-driven QCD templates
    TH1D *h_QCD_est = ((TH1D*)(h_data_mass_dijet->Clone("h_QCD_est")));
    TH1D *h_QCD_est_SS = ((TH1D*)(h_data_mass_SS_dijet->Clone("h_QCD_est_SS")));
    h_QCD_est->Add(h_MC_mass_dijet[_DY_Full], -1);
    h_QCD_est->Add(h_MC_mass_dijet[_ttbar], -1);
    h_QCD_est_SS->Add(h_MC_mass_SS_dijet[_DY_Full], -1);
    h_QCD_est_SS->Add(h_MC_mass_SS_dijet[_ttbar], -1);
    removeNegativeBins(h_QCD_est);
    removeNegativeBins(h_QCD_est_SS);

    // Making data-driven W+Jets template from same-sign distributions
    TH1D *h_WJets_est_SS = ((TH1D*)(h_data_mass_SS->Clone("h_WJets_est_SS")));
    h_WJets_est_SS->Add(h_MC_mass_SS[_ttbar], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_DY_Full], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_tW], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_tbarW], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_WW], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_WZ], -1);
    h_WJets_est_SS->Add(h_MC_mass_SS[_ZZ], -1);
    h_WJets_est_SS->Add(h_QCD_est_SS, -2);
    removeNegativeBins(h_WJets_est_SS);
    THStack *s = new THStack("a","A");
    s->Add(h_MC_mass_SS[_WJets]);
    s->Add(h_MC_mass_SS[_ttbar]);
    s->Add(h_QCD_est_SS);
    s->Add(h_QCD_est_SS);
    myRatioPlot_t *p = new myRatioPlot_t("asd", s, h_data_mass_SS);
    p->SetPlots("mass", 0, 3000);
//    p->Draw(0.1,1e4, 1);
    TCanvas *c_test = new TCanvas("test", "WJets template", 800, 800);
    h_WJets_est_SS->Draw("hist");
    c_test->SetLogx();
    c_test->Update();

    // Making RooDataHist
    Double_t upper_limit = 15;
    if (type == 1) upper_limit = 200;
    else if (type == 2) upper_limit = 3000;
    RooRealVar mass("mass", "m_{#mu#mu} [GeV]", 15, upper_limit);

    RooDataHist *rh_mass_QCD   = new RooDataHist("rh_mass_QCD",   "RooHist_mass_QCD",          mass, h_QCD_est);
    RooDataHist *rh_mass_WJets = new RooDataHist("rh_mass_WJets", "RooHist_mass_WJets",        mass, h_WJets_est_SS);
    RooDataHist *rh_mass_DY    = new RooDataHist("rh_mass_DY",    "RooHist_endcap_deno_DY",    mass, h_MC_mass[_DY_Full]);
    RooDataHist *rh_mass_ttbar = new RooDataHist("rh_mass_ttbar", "RooHist_endcap_deno_ttbar", mass, h_MC_mass[_ttbar]);
    RooDataHist *rh_mass_tW    = new RooDataHist("rh_mass_tW",    "RooHist_endcap_deno_tW",    mass, h_MC_mass[_tW]);
    RooDataHist *rh_mass_tbarW = new RooDataHist("rh_mass_tbarW", "RooHist_endcap_deno_tbarW", mass, h_MC_mass[_tbarW]);
    RooDataHist *rh_mass_WW    = new RooDataHist("rh_mass_WW",    "RooHist_endcap_deno_WW",    mass, h_MC_mass[_WW]);
    RooDataHist *rh_mass_WZ    = new RooDataHist("rh_mass_WZ",    "RooHist_endcap_deno_WZ",    mass, h_MC_mass[_WZ]);
    RooDataHist *rh_mass_ZZ    = new RooDataHist("rh_mass_ZZ",    "RooHist_endcap_deno_ZZ",    mass, h_MC_mass[_ZZ]);
    RooDataHist *rh_mass_data  = new RooDataHist("rh_mass_data",  "RooHist_endcap_deno_data",  mass, h_data_mass);

    // Making RooHistPdf
    RooHistPdf *pdf_mass_QCD   = new RooHistPdf("pdf_mass_QCD",   "MC QCD mass template",    mass, *rh_mass_QCD,   0);
    RooHistPdf *pdf_mass_WJets = new RooHistPdf("pdf_mass_WJets", "MC W+Jets mass template", mass, *rh_mass_WJets, 0);
    RooHistPdf *pdf_mass_DY    = new RooHistPdf("pdf_mass_DY",    "MC DY mass template",     mass, *rh_mass_DY,    0);
    RooHistPdf *pdf_mass_ttbar = new RooHistPdf("pdf_mass_ttbar", "MC ttbar mass template",  mass, *rh_mass_ttbar, 0);
    RooHistPdf *pdf_mass_tW    = new RooHistPdf("pdf_mass_tW",    "MC tW mass template",     mass, *rh_mass_tW,    0);
    RooHistPdf *pdf_mass_tbarW = new RooHistPdf("pdf_mass_tbarW", "MC tbarW template",       mass, *rh_mass_tbarW, 0);
    RooHistPdf *pdf_mass_WW    = new RooHistPdf("pdf_mass_WW",    "MC WW mass template",     mass, *rh_mass_WW,    0);
    RooHistPdf *pdf_mass_WZ    = new RooHistPdf("pdf_mass_WZ",    "MC WZ mass template",     mass, *rh_mass_WZ,    0);
    RooHistPdf *pdf_mass_ZZ    = new RooHistPdf("pdf_mass_ZZ",    "MC ZZ mass template",     mass, *rh_mass_ZZ,    0);

    // Constraints for integrals
    Double_t N_mass_ttbar = h_MC_mass[_ttbar]->Integral();
    Double_t N_mass_WJets = h_WJets_est_SS->Integral() * 3;
    Double_t N_mass_DY    = h_MC_mass[_DY_Full]->Integral();
    Double_t N_mass_QCD   = h_QCD_est->Integral() * 2;
    Double_t N_mass_tW    = h_MC_mass[_tW]->Integral();
    Double_t N_mass_tbarW = h_MC_mass[_tbarW]->Integral();
    Double_t N_mass_WW    = h_MC_mass[_WW]->Integral();
    Double_t N_mass_WZ    = h_MC_mass[_WZ]->Integral();
    Double_t N_mass_ZZ    = h_MC_mass[_ZZ]->Integral();
    cout << "# of same sign W+Jets events: " << h_WJets_est_SS->Integral() << endl;
    cout << "Expected # of opposite sign W+Jets events: ~" << N_mass_WJets << endl;

    Double_t N_mass_total = N_mass_ttbar + N_mass_WJets + N_mass_DY +
                            N_mass_QCD/*   + N_mass_tW    + N_mass_tbarW +
                            N_mass_WW    + N_mass_WZ    + N_mass_ZZ*/;

    Double_t Nnorm_mass_ttbar = N_mass_ttbar * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WJets = N_mass_WJets * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_DY    = N_mass_DY    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_tW    = N_mass_tW    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_tbarW = N_mass_tbarW * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WW    = N_mass_WW    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WZ    = N_mass_WZ    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_ZZ    = N_mass_ZZ    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_QCD   = N_mass_QCD   * h_data_mass->Integral() / N_mass_total;

    // Fit constraints
    RooRealVar n_mass_ttbar("n_mass_ttbar", "n_mass_ttbar", Nnorm_mass_ttbar, Nnorm_mass_ttbar *0.8, Nnorm_mass_ttbar *1.2);
    RooRealVar n_mass_WJets("n_mass_WJets", "n_mass_WJets", Nnorm_mass_WJets, Nnorm_mass_WJets *0.5,  Nnorm_mass_WJets *1.5 );
    RooRealVar n_mass_DY   ("n_mass_DY",    "n_mass_DY",    Nnorm_mass_DY,    Nnorm_mass_DY    *0.8,  Nnorm_mass_DY    *1.2 );
    RooRealVar n_mass_tW   ("n_mass_tW",    "n_mass_tW",    Nnorm_mass_tW,    Nnorm_mass_tW    *0.8,  Nnorm_mass_tW    *1.2 );
    RooRealVar n_mass_tbarW("n_mass_tbarW", "n_mass_tbarW", Nnorm_mass_tbarW, Nnorm_mass_tbarW *0.8,  Nnorm_mass_tbarW *1.2 );
    RooRealVar n_mass_WW   ("n_mass_WW",    "n_mass_WW",    Nnorm_mass_WW,    Nnorm_mass_WW    *0.8,  Nnorm_mass_WW    *1.2 );
    RooRealVar n_mass_WZ   ("n_mass_WZ",    "n_mass_WZ",    Nnorm_mass_WZ,    Nnorm_mass_WZ    *0.8,  Nnorm_mass_WZ    *1.2 );
    RooRealVar n_mass_ZZ   ("n_mass_ZZ",    "n_mass_ZZ",    Nnorm_mass_ZZ,    Nnorm_mass_ZZ    *0.8,  Nnorm_mass_ZZ    *1.2 );
    RooRealVar n_mass_QCD  ("n_mass_QCD",   "n_mass_QCD",   Nnorm_mass_QCD,   Nnorm_mass_QCD   *0.5,  Nnorm_mass_QCD   *1.5 );

    // Models
    RooAddPdf model_mass("model_mass", "model_mass",
                         RooArgList(*pdf_mass_QCD,   *pdf_mass_WJets, *pdf_mass_DY,
                                    *pdf_mass_ttbar, *pdf_mass_tW,    *pdf_mass_tbarW,
                                    *pdf_mass_WW,    *pdf_mass_WZ,    *pdf_mass_ZZ),
                         RooArgList(n_mass_QCD,   n_mass_WJets, n_mass_DY,
                                    n_mass_ttbar, n_mass_tW,    n_mass_tbarW,
                                    n_mass_WW,    n_mass_WZ,    n_mass_ZZ));

    // Fitting
    RooFitResult* fit_mass = model_mass.fitTo(*rh_mass_data, Save());


    /// DRAWING
    cout << "\n----- FIT RESULT -----" << endl;
    TCanvas *c_fit = new TCanvas("c_fit", "c_fit", 800, 800);
    c_fit->cd();

    //Top Pad
    TPad *c1 = new TPad("padc1","padc1",0.01,0.01,0.99,0.99);
    c1->Draw();
    c1->cd();
    c1->SetTopMargin(0.01);
    c1->SetBottomMargin(0.35);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.13);
    c1->SetFillStyle(1);
    c1->SetLogx();
    c1->SetLogy();

    // Main stack histogram
    RooPlot *frame = mass.frame(Title(" "));
    rh_mass_data->plotOn(frame, DataError(RooAbsData::SumW2));

    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar,pdf_mass_DY"),
                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJetspdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar"),
                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW"),
                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW"),
                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW"),
                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ,pdf_mass_WZ"),
                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ZZ"),
                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets"),
                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD"),
                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
/*
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ttbar,pdf_mass_DY"),
                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets,pdf_mass_ttbar"),
                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD,pdf_mass_WJets"),
                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_QCD"),
                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
*/

    rh_mass_data->plotOn(frame, DataError(RooAbsData::SumW2));
    frame->GetYaxis()->SetRangeUser(1e0, 2e6);
    frame->Draw();
    fit_mass->Print();

    // Legend
    TLegend *legend = new TLegend(0.65, 0.8, 0.95, 0.97);
    legend->SetFillColor(kWhite);
    legend->SetNColumns(2);
//    legend->SetLineColor(kWhite);
    legend->AddEntry(frame->nameOf(0), "Data", "LP");
    legend->AddEntry(frame->nameOf(1), "DY", "F");
    legend->AddEntry(frame->nameOf(2), "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "F");
    legend->AddEntry(frame->nameOf(3), "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "F");
    legend->AddEntry(frame->nameOf(4), "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "F");
    legend->AddEntry(frame->nameOf(5), "#font[12]{#scale[1.1]{WW}}", "F");
    legend->AddEntry(frame->nameOf(6), "#font[12]{#scale[1.1]{WZ}}", "F");
    legend->AddEntry(frame->nameOf(7), "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "F");
    legend->AddEntry(frame->nameOf(8), "#font[12]{#scale[1.1]{W}}+Jets", "F");
    legend->AddEntry(frame->nameOf(9), "#font[12]{#scale[1.1]{QCD}}", "F");

//    legend->SetNColumns(2);

    legend->Draw();

    frame->GetYaxis()->SetTitle("Number of entries");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);

    // Legend
    legend->Draw();

    // Bottom pad
    TPad *c2 = new TPad("padc2","padc2",0.01,0.01,0.99,0.35);
    c2->Draw();
    c2->cd();
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.33);
    c2->SetRightMargin(0.02);
    c2->SetLeftMargin(0.12);
    c2->SetFillStyle(0);
    c2->SetLogx();
    c2->SetGrid();

    // Ratio plot
    TH1D *h_mass_MC_fit = ((TH1D*)(model_mass.createHistogram("h_mass_MC_fit", mass)));
    Double_t N_mass_data = h_data_mass  ->Integral();
    Double_t N_mass_MC   = h_mass_MC_fit->Integral();
    h_mass_MC_fit->Scale(N_mass_data/N_mass_MC); // Why is this necessary???
    cout << "\nData integral: " << N_mass_data << endl;
    cout << "MC integral: "     << h_mass_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_data_mass->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_mass_MC_fit->GetBinContent(1) << endl;

    TH1D *h_mass_ratio = ((TH1D*)(h_data_mass->Clone("h_mass_ratio")));
    h_data_mass->Sumw2(); h_mass_MC_fit->Sumw2();
    h_mass_ratio->Divide(h_data_mass, h_mass_MC_fit);
    h_mass_ratio->SetTitle("");
    h_mass_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_mass_ratio->GetXaxis()->SetNoExponent(1);
    h_mass_ratio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
    h_mass_ratio->GetXaxis()->SetTitleSize(0.17);
    h_mass_ratio->GetXaxis()->SetLabelSize(0.125);
    h_mass_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_mass_ratio->GetYaxis()->SetTitle("Data/MC");
    h_mass_ratio->GetYaxis()->SetTitleSize(0.114);
    h_mass_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_mass_ratio->GetYaxis()->SetLabelSize(0.11);
    h_mass_ratio->GetYaxis()->SetTickLength(0.01);
    h_mass_ratio->GetYaxis()->SetDecimals(1);
    h_mass_ratio->SetMaximum(1.25);
    h_mass_ratio->SetMinimum(0.75);
    h_mass_ratio->GetYaxis()->SetNdivisions(5);
    h_mass_ratio->SetLineWidth(1);
    h_mass_ratio->SetLineColor(kBlack);
    h_mass_ratio->SetMarkerStyle(kFullDotLarge);
    h_mass_ratio->SetMarkerColor(kBlack);
    h_mass_ratio->SetStats(kFALSE);

    h_mass_ratio->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line = ((TH1D*)(h_data_mass->Clone("h_line")));
    h_line->Reset("ICES");
    for (Int_t i=1; i<=h_line->GetNbinsX(); i++)
        h_line->SetBinContent(i, 1);
    h_line->SetLineColor(kRed);
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_mass = model_mass.createChi2(*rh_mass_data);
    cout << "chi2: " << chi2_mass->getVal() << endl;
    cout << "Normalized chi2: " << chi2_mass->getVal() / ((Double_t)h_data_mass->GetNbinsX()) << endl;

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

} // End of Mu_WJETSest_TFit()


void Mu_WJETSest_Tfit_legacy()
{
    FileMgr Mgr;

    TH1D *h_MC_mass  [_EndOf_Data_Special],
         *h_data_mass;

    TFile *file = new TFile("/media/sf_DATA/FR/Muon/WJETest_Mu.root", "READ");

// ############################# SETUP ################################# //
//----------------------------- MC bkg ------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        file->GetObject("h_mass_"+Mgr.Procname[pr1],   h_MC_mass[pr1]);

        removeNegativeBins(h_MC_mass[pr1]);
        h_MC_mass[pr1]->SetDirectory(0);

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
    h_MC_mass[_ttbar]->Add(h_MC_mass[_ttbar_700to1000]);
    h_MC_mass[_ttbar]->Add(h_MC_mass[_ttbar_1000toInf]);
    h_MC_mass[_WJets]->Add(h_MC_mass[_WJets_ext2v5]);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        file->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass[pr]);

        removeNegativeBins(h_MC_mass[pr]);
        h_MC_mass[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_MC_mass[_DY_Full] = ((TH1D*)(h_MC_mass[pr]->Clone("h_MC_mass_DY")));
            h_MC_mass[_DY_Full]->SetDirectory(0);
        }
        else
            h_MC_mass[_DY_Full]->Add(h_MC_mass[pr]);
    }

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        file->GetObject("h_mass_"+Mgr.Procname[pr], h_MC_mass[pr]);

        removeNegativeBins(h_MC_mass[pr]);
        h_MC_mass[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_MC_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_MC_mass[pr]->Clone("h_MC_mass_QCD")));
            h_MC_mass[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
            h_MC_mass[_QCDMuEnriched_Full]->Add(h_MC_mass[pr]);
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TH1D *h_temp;
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_mass_"+Mgr.Procname[pr], h_data_mass);
            removeNegativeBins(h_data_mass);
        }
        else
        {
            file->GetObject("h_mass_"+Mgr.Procname[pr], h_temp);
            removeNegativeBins(h_temp);
            h_data_mass->Add(h_temp);
        }
    }
    h_data_mass->SetDirectory(0);

// ######################## MODEL BUILDING ##########################

    // Making RooDataHist
    RooRealVar mass("mass", "m_{#mu#mu} [GeV]", 15, 3000);

    RooDataHist *rh_mass_QCD   = new RooDataHist("rh_mass_QCD",   "RooHist_mass_QCD",          mass, h_MC_mass[_QCDMuEnriched_Full]);
    RooDataHist *rh_mass_WJets = new RooDataHist("rh_mass_WJets", "RooHist_mass_WJets",        mass, h_MC_mass[_WJets]);
    RooDataHist *rh_mass_DY    = new RooDataHist("rh_mass_DY",    "RooHist_endcap_deno_DY",    mass, h_MC_mass[_DY_Full]);
    RooDataHist *rh_mass_ttbar = new RooDataHist("rh_mass_ttbar", "RooHist_endcap_deno_ttbar", mass, h_MC_mass[_ttbar]);
    RooDataHist *rh_mass_tW    = new RooDataHist("rh_mass_tW",    "RooHist_endcap_deno_tW",    mass, h_MC_mass[_tW]);
    RooDataHist *rh_mass_tbarW = new RooDataHist("rh_mass_tbarW", "RooHist_endcap_deno_tbarW", mass, h_MC_mass[_tbarW]);
    RooDataHist *rh_mass_WW    = new RooDataHist("rh_mass_WW",    "RooHist_endcap_deno_WW",    mass, h_MC_mass[_WW]);
    RooDataHist *rh_mass_WZ    = new RooDataHist("rh_mass_WZ",    "RooHist_endcap_deno_WZ",    mass, h_MC_mass[_WZ]);
    RooDataHist *rh_mass_ZZ    = new RooDataHist("rh_mass_ZZ",    "RooHist_endcap_deno_ZZ",    mass, h_MC_mass[_ZZ]);
    RooDataHist *rh_mass_data  = new RooDataHist("rh_mass_data",  "RooHist_endcap_deno_data",  mass, h_data_mass);

    // Making RooHistPdf
    RooHistPdf *pdf_mass_QCD   = new RooHistPdf("pdf_mass_QCD",   "MC QCD mass template",    mass, *rh_mass_QCD,   0);
    RooHistPdf *pdf_mass_WJets = new RooHistPdf("pdf_mass_WJets", "MC W+Jets mass template", mass, *rh_mass_WJets, 0);
    RooHistPdf *pdf_mass_DY    = new RooHistPdf("pdf_mass_DY",    "MC DY mass template",     mass, *rh_mass_DY,    0);
    RooHistPdf *pdf_mass_ttbar = new RooHistPdf("pdf_mass_ttbar", "MC ttbar mass template",  mass, *rh_mass_ttbar, 0);
    RooHistPdf *pdf_mass_tW    = new RooHistPdf("pdf_mass_tW",    "MC tW mass template",     mass, *rh_mass_tW,    0);
    RooHistPdf *pdf_mass_tbarW = new RooHistPdf("pdf_mass_tbarW", "MC tbarW template",       mass, *rh_mass_tbarW, 0);
    RooHistPdf *pdf_mass_WW    = new RooHistPdf("pdf_mass_WW",    "MC WW mass template",     mass, *rh_mass_WW,    0);
    RooHistPdf *pdf_mass_WZ    = new RooHistPdf("pdf_mass_WZ",    "MC WZ mass template",     mass, *rh_mass_WZ,    0);
    RooHistPdf *pdf_mass_ZZ    = new RooHistPdf("pdf_mass_ZZ",    "MC ZZ mass template",     mass, *rh_mass_ZZ,    0);

    // Constraints for integrals
    Double_t N_mass_ttbar = h_MC_mass[_ttbar]->Integral();
    Double_t N_mass_WJets = h_MC_mass[_WJets]->Integral();
    Double_t N_mass_DY    = h_MC_mass[_DY_Full]->Integral();
    Double_t N_mass_QCD   = h_MC_mass[_QCDMuEnriched_Full]->Integral();
    Double_t N_mass_tW    = h_MC_mass[_tW]->Integral();
    Double_t N_mass_tbarW = h_MC_mass[_tbarW]->Integral();
    Double_t N_mass_WW    = h_MC_mass[_WW]->Integral();
    Double_t N_mass_WZ    = h_MC_mass[_WZ]->Integral();
    Double_t N_mass_ZZ    = h_MC_mass[_ZZ]->Integral();

    Double_t N_mass_total = N_mass_ttbar + N_mass_WJets + N_mass_DY +
                            N_mass_QCD   + N_mass_tW    + N_mass_tbarW +
                            N_mass_WW    + N_mass_WZ    + N_mass_ZZ;

    Double_t Nnorm_mass_ttbar = N_mass_ttbar * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WJets = N_mass_WJets * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_DY    = N_mass_DY    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_tW    = N_mass_tW    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_tbarW = N_mass_tbarW * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WW    = N_mass_WW    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_WZ    = N_mass_WZ    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_ZZ    = N_mass_ZZ    * h_data_mass->Integral() / N_mass_total;
    Double_t Nnorm_mass_QCD   = N_mass_QCD   * h_data_mass->Integral() / N_mass_total;

    // Fit constraints
    RooRealVar n_mass_ttbar("n_mass_ttbar", "n_mass_ttbar", Nnorm_mass_ttbar, Nnorm_mass_ttbar *0.75, Nnorm_mass_ttbar *1.25);
    RooRealVar n_mass_WJets("n_mass_WJets", "n_mass_WJets", Nnorm_mass_WJets, Nnorm_mass_WJets *0.5,  Nnorm_mass_WJets *1.5 );
    RooRealVar n_mass_DY   ("n_mass_DY",    "n_mass_DY",    Nnorm_mass_DY,    Nnorm_mass_DY    *0.8,  Nnorm_mass_DY    *1.2 );
    RooRealVar n_mass_tW   ("n_mass_tW",    "n_mass_tW",    Nnorm_mass_tW,    Nnorm_mass_tW    *0.8,  Nnorm_mass_tW    *1.2 );
    RooRealVar n_mass_tbarW("n_mass_tbarW", "n_mass_tbarW", Nnorm_mass_tbarW, Nnorm_mass_tbarW *0.8,  Nnorm_mass_tbarW *1.2 );
    RooRealVar n_mass_WW   ("n_mass_WW",    "n_mass_WW",    Nnorm_mass_WW,    Nnorm_mass_WW    *0.8,  Nnorm_mass_WW    *1.2 );
    RooRealVar n_mass_WZ   ("n_mass_WZ",    "n_mass_WZ",    Nnorm_mass_WZ,    Nnorm_mass_WZ    *0.8,  Nnorm_mass_WZ    *1.2 );
    RooRealVar n_mass_ZZ   ("n_mass_ZZ",    "n_mass_ZZ",    Nnorm_mass_ZZ,    Nnorm_mass_ZZ    *0.8,  Nnorm_mass_ZZ    *1.2 );
    RooRealVar n_mass_QCD  ("n_mass_QCD",   "n_mass_QCD",   Nnorm_mass_QCD,   Nnorm_mass_QCD   *0.5,  Nnorm_mass_QCD   *1.5 );

    // Models
    RooAddPdf model_mass("model_mass", "model_mass",
                         RooArgList(*pdf_mass_QCD,   *pdf_mass_WJets, *pdf_mass_DY,
                                    *pdf_mass_ttbar, *pdf_mass_tW,    *pdf_mass_tbarW,
                                    *pdf_mass_WW,    *pdf_mass_WZ,    *pdf_mass_ZZ),
                         RooArgList(n_mass_QCD,   n_mass_WJets, n_mass_DY,
                                    n_mass_ttbar, n_mass_tW,    n_mass_tbarW,
                                    n_mass_WW,    n_mass_WZ,    n_mass_ZZ));

    // Fitting
    RooFitResult* fit_mass = model_mass.fitTo(*rh_mass_data, Save());


    /// DRAWING
    cout << "\n----- FIT RESULT -----" << endl;
    TCanvas *c_fit = new TCanvas("c_fit", "c_fit", 800, 800);
    c_fit->cd();

    //Top Pad
    TPad *c1 = new TPad("padc1","padc1",0.01,0.01,0.99,0.99);
    c1->Draw();
    c1->cd();
    c1->SetTopMargin(0.01);
    c1->SetBottomMargin(0.35);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.13);
    c1->SetFillStyle(1);
    c1->SetLogx();
    c1->SetLogy();

    // Main stack histogram
    RooPlot *frame = mass.frame(Title(" "));
    rh_mass_data->plotOn(frame, DataError(RooAbsData::SumW2));

    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar,"
                                        "pdf_mass_DY,pdf_mass_WJets,pdf_mass_QCD"),
                      LineColor(0), FillColor(kRed+3), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar,"
                                        "pdf_mass_DY,pdf_mass_WJets"),
                      LineColor(0), FillColor(kRed-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar,"
                                        "pdf_mass_DY"),
                      LineColor(0), FillColor(kOrange), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW,pdf_mass_ttbar"),
                      LineColor(0), FillColor(kCyan+2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW,pdf_mass_tbarW"),
                      LineColor(0), FillColor(kGreen-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW,"
                                        "pdf_mass_tW"),
                      LineColor(0), FillColor(kGreen+2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ,pdf_mass_WW"),
                      LineColor(0), FillColor(kMagenta-5), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ,pdf_mass_WZ"),
                      LineColor(0), FillColor(kMagenta-2), DrawOption("F"));
    model_mass.plotOn(frame, Components("pdf_mass_ZZ"),
                      LineColor(0), FillColor(kMagenta-6), DrawOption("F"));

    rh_mass_data->plotOn(frame, DataError(RooAbsData::SumW2));
    frame->GetYaxis()->SetRangeUser(1e0, 2e6);
    frame->Draw();
    fit_mass->Print();

    // Legend
    TLegend *legend = new TLegend(0.65, 0.8, 0.95, 0.97);
    legend->SetFillColor(kWhite);
//    legend->SetLineColor(kWhite);
    legend->AddEntry(frame->nameOf(0), "Data", "LP");
    legend->AddEntry(frame->nameOf(1), "#font[12]{#scale[1.1]{QCD}}", "F");
    legend->AddEntry(frame->nameOf(2), "#font[12]{#scale[1.1]{W}}+Jets", "F");
    legend->AddEntry(frame->nameOf(3), "DY", "F");
    legend->AddEntry(frame->nameOf(4), "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "F");
    legend->AddEntry(frame->nameOf(5), "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "F");
    legend->AddEntry(frame->nameOf(6), "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "F");
    legend->AddEntry(frame->nameOf(7), "#font[12]{#scale[1.1]{WW}}", "F");
    legend->AddEntry(frame->nameOf(8), "#font[12]{#scale[1.1]{WZ}}", "F");
    legend->AddEntry(frame->nameOf(9), "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "F");
    legend->SetNColumns(2);

    legend->Draw();

    frame->GetYaxis()->SetTitle("Number of entries");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);

    // Legend
    legend->Draw();

    // Bottom pad
    TPad *c2 = new TPad("padc2","padc2",0.01,0.01,0.99,0.35);
    c2->Draw();
    c2->cd();
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.33);
    c2->SetRightMargin(0.02);
    c2->SetLeftMargin(0.12);
    c2->SetFillStyle(0);
    c2->SetLogx();
    c2->SetGrid();

    // Ratio plot
    TH1D *h_mass_MC_fit = ((TH1D*)(model_mass.createHistogram("h_mass_MC_fit", mass)));
    Double_t N_mass_data = h_data_mass  ->Integral();
    Double_t N_mass_MC   = h_mass_MC_fit->Integral();
    h_mass_MC_fit->Scale(N_mass_data/N_mass_MC); // Why is this necessary???
    cout << "\nData integral: " << N_mass_data << endl;
    cout << "MC integral: "     << h_mass_MC_fit->Integral() << endl;
    cout << "Data in 1st bin: " << h_data_mass->GetBinContent(1) << endl;
    cout << "MC in 1st bin: "   << h_mass_MC_fit->GetBinContent(1) << endl;

    TH1D *h_mass_ratio = ((TH1D*)(h_data_mass->Clone("h_mass_ratio")));
    h_data_mass->Sumw2(); h_mass_MC_fit->Sumw2();
    h_mass_ratio->Divide(h_data_mass, h_mass_MC_fit);
    h_mass_ratio->SetTitle("");
    h_mass_ratio->GetXaxis()->SetMoreLogLabels(1);
    h_mass_ratio->GetXaxis()->SetNoExponent(1);
    h_mass_ratio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
    h_mass_ratio->GetXaxis()->SetTitleSize(0.17);
    h_mass_ratio->GetXaxis()->SetLabelSize(0.125);
    h_mass_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_mass_ratio->GetYaxis()->SetTitle("Data/MC");
    h_mass_ratio->GetYaxis()->SetTitleSize(0.114);
    h_mass_ratio->GetYaxis()->SetTitleOffset(0.48);
    h_mass_ratio->GetYaxis()->SetLabelSize(0.11);
    h_mass_ratio->GetYaxis()->SetTickLength(0.01);
    h_mass_ratio->GetYaxis()->SetDecimals(1);
    h_mass_ratio->SetMaximum(1.25);
    h_mass_ratio->SetMinimum(0.75);
    h_mass_ratio->GetYaxis()->SetNdivisions(5);
    h_mass_ratio->SetLineWidth(1);
    h_mass_ratio->SetLineColor(kBlack);
    h_mass_ratio->SetMarkerStyle(kFullDotLarge);
    h_mass_ratio->SetMarkerColor(kBlack);
    h_mass_ratio->SetStats(kFALSE);

    h_mass_ratio->Draw("E1P");

    // Red line at Data/MC=1
    TH1D *h_line = ((TH1D*)(h_data_mass->Clone("h_line")));
    h_line->Reset("ICES");
    for (Int_t i=1; i<=h_line->GetNbinsX(); i++)
        h_line->SetBinContent(i, 1);
    h_line->SetLineColor(kRed);
    h_line->Draw("LSAME");

    //Chi^2
    RooAbsReal *chi2_mass = model_mass.createChi2(*rh_mass_data);
    cout << "chi2: " << chi2_mass->getVal() << endl;
    cout << "Normalized chi2: " << chi2_mass->getVal() / ((Double_t)h_data_mass->GetNbinsX()) << endl;


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

} // End of Mu_WJETSest_TFit_legacy()


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
