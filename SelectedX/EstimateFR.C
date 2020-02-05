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
#include <TText.h>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

void E_EstFR (Int_t type);
void Mu_EstFR (Int_t type);
Double_t CompChiSquared (TH1D *h_data, THStack *s_MC);
Double_t CompAvgDataMCDifference (TH1D *h_data, THStack *s_MC);
void removeNegativeBins(TH1D *h);

const Double_t L_B2F = 19721.0;
const Double_t L_G2H = 16146.0;
const Double_t L_B2H = 35867.0;

void EstimateFR (TString WhichX = "", Int_t type = 2)
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
        cout << "\n*******      E_EstFR(" << type << ")      *******" << endl;
        E_EstFR(type);
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*******     Mu_EstFR(" << type << ")     *******" << endl;
        Mu_EstFR(type);
    }
    if (Xselected == 0) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ############################################################################# ///
/// ----------------------------- Electron Channel ------------------------------ ///
/// ############################################################################# ///
void E_EstFR(Int_t type)
{
   return; // NOT READY YET
} // End of E_EstFR()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_EstFR(Int_t type)
{
    FileMgr fm;

    TH1D *h_pT_barrel_deno, *h_pT_barrel_nume, *h_pT_endcap_deno, *h_pT_endcap_nume, *h_pT_barrel_ctrl, *h_pT_endcap_ctrl,
         *h_FRratio_barrel, *h_FRratio_endcap, *h_FRtemplate_barrel, *h_FRtemplate_endcap,
         *h_FRdalmin_barrel, *h_FRdalmin_endcap;

    TH1D *h_pT_barrel_MC_deno[_EndOf_Data_Special], *h_pT_barrel_MC_nume[_EndOf_Data_Special],
         *h_pT_endcap_MC_deno[_EndOf_Data_Special], *h_pT_endcap_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_data_deno, *h_pT_barrel_data_nume, *h_pT_endcap_data_deno, *h_pT_endcap_data_nume,
         *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl;

    TH1D *h_pT_barrel_template_deno_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_deno_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_deno_100to500[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_100to500[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_100to500[_EndOf_Data_Special];


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
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr1]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr1]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);

        removeNegativeBins(h_pT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);

        h_pT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);

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
    h_pT_barrel_MC_deno[_ttbar]->Add(h_pT_barrel_MC_deno[_ttbar_700to1000]);
    h_pT_endcap_MC_deno[_ttbar]->Add(h_pT_endcap_MC_deno[_ttbar_700to1000]);
    h_pT_barrel_MC_ctrl[_ttbar]->Add(h_pT_barrel_MC_ctrl[_ttbar_700to1000]);
    h_pT_endcap_MC_ctrl[_ttbar]->Add(h_pT_endcap_MC_ctrl[_ttbar_700to1000]);
    h_pT_barrel_MC_nume[_ttbar]->Add(h_pT_barrel_MC_nume[_ttbar_700to1000]);
    h_pT_endcap_MC_nume[_ttbar]->Add(h_pT_endcap_MC_nume[_ttbar_700to1000]);
    h_pT_barrel_MC_deno[_ttbar]->Add(h_pT_barrel_MC_deno[_ttbar_1000toInf]);
    h_pT_endcap_MC_deno[_ttbar]->Add(h_pT_endcap_MC_deno[_ttbar_1000toInf]);
    h_pT_barrel_MC_ctrl[_ttbar]->Add(h_pT_barrel_MC_ctrl[_ttbar_1000toInf]);
    h_pT_endcap_MC_ctrl[_ttbar]->Add(h_pT_endcap_MC_ctrl[_ttbar_1000toInf]);
    h_pT_barrel_MC_nume[_ttbar]->Add(h_pT_barrel_MC_nume[_ttbar_1000toInf]);
    h_pT_endcap_MC_nume[_ttbar]->Add(h_pT_endcap_MC_nume[_ttbar_1000toInf]);
    h_pT_barrel_MC_deno[_WJets]->Add(h_pT_barrel_MC_deno[_WJets_ext2v5]);
    h_pT_endcap_MC_deno[_WJets]->Add(h_pT_endcap_MC_deno[_WJets_ext2v5]);
    h_pT_barrel_MC_ctrl[_WJets]->Add(h_pT_barrel_MC_ctrl[_WJets_ext2v5]);
    h_pT_endcap_MC_ctrl[_WJets]->Add(h_pT_endcap_MC_ctrl[_WJets_ext2v5]);
    h_pT_barrel_MC_nume[_WJets]->Add(h_pT_barrel_MC_nume[_WJets_ext2v5]);
    h_pT_endcap_MC_nume[_WJets]->Add(h_pT_endcap_MC_nume[_WJets_ext2v5]);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);

        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_DY")));
            h_pT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_DY")));
            h_pT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_deno[_DY_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_DY_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
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
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno_50to70",   h_pT_barrel_template_deno_50to70  [pr]);
        file->GetObject("h_pT_endcap_deno_50to70",   h_pT_endcap_template_deno_50to70  [pr]);
        file->GetObject("h_pT_barrel_ctrl_50to70",   h_pT_barrel_template_ctrl_50to70  [pr]);
        file->GetObject("h_pT_endcap_ctrl_50to70",   h_pT_endcap_template_ctrl_50to70  [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_pT_barrel_template_nume_50to70  [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_pT_endcap_template_nume_50to70  [pr]);
        file->GetObject("h_pT_barrel_deno_70to100",  h_pT_barrel_template_deno_70to100 [pr]);
        file->GetObject("h_pT_endcap_deno_70to100",  h_pT_endcap_template_deno_70to100 [pr]);
        file->GetObject("h_pT_barrel_ctrl_70to100",  h_pT_barrel_template_ctrl_70to100 [pr]);
        file->GetObject("h_pT_endcap_ctrl_70to100",  h_pT_endcap_template_ctrl_70to100 [pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_pT_barrel_template_nume_70to100 [pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_pT_endcap_template_nume_70to100 [pr]);
        file->GetObject("h_pT_barrel_deno_100to500", h_pT_barrel_template_deno_100to500[pr]);
        file->GetObject("h_pT_endcap_deno_100to500", h_pT_endcap_template_deno_100to500[pr]);
        file->GetObject("h_pT_barrel_ctrl_100to500", h_pT_barrel_template_ctrl_100to500[pr]);
        file->GetObject("h_pT_endcap_ctrl_100to500", h_pT_endcap_template_ctrl_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_pT_barrel_template_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_pT_endcap_template_nume_100to500[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_template_deno_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_deno_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_nume_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_nume_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_deno_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_deno_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_nume_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_nume_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_deno_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_deno_100to500[pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_100to500[pr]);
        removeNegativeBins(h_pT_barrel_template_nume_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_nume_100to500[pr]);

        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_template_deno_50to70  [pr]->SetDirectory(0);
        h_pT_endcap_template_deno_50to70  [pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_50to70  [pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_50to70  [pr]->SetDirectory(0);
        h_pT_barrel_template_nume_50to70  [pr]->SetDirectory(0);
        h_pT_endcap_template_nume_50to70  [pr]->SetDirectory(0);
        h_pT_barrel_template_deno_70to100 [pr]->SetDirectory(0);
        h_pT_endcap_template_deno_70to100 [pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_70to100 [pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_70to100 [pr]->SetDirectory(0);
        h_pT_barrel_template_nume_70to100 [pr]->SetDirectory(0);
        h_pT_endcap_template_nume_70to100 [pr]->SetDirectory(0);
        h_pT_barrel_template_deno_100to500[pr]->SetDirectory(0);
        h_pT_endcap_template_deno_100to500[pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_100to500[pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_100to500[pr]->SetDirectory(0);
        h_pT_barrel_template_nume_100to500[pr]->SetDirectory(0);
        h_pT_endcap_template_nume_100to500[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_QCD")));
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_QCD")));
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_QCD")));
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_QCD")));
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_QCD")));
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_50to70  [pr]->Clone("h_pT_barrel_deno_QCD_50to70"  )));
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_50to70  [pr]->Clone("h_pT_endcap_deno_QCD_50to70"  )));
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_50to70  [pr]->Clone("h_pT_barrel_ctrl_QCD_50to70"  )));
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_50to70  [pr]->Clone("h_pT_endcap_ctrl_QCD_50to70"  )));
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_50to70  [pr]->Clone("h_pT_barrel_nume_QCD_50to70"  )));
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_50to70  [pr]->Clone("h_pT_endcap_nume_QCD_50to70"  )));
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_70to100 [pr]->Clone("h_pT_barrel_deno_QCD_70to100" )));
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_70to100 [pr]->Clone("h_pT_endcap_deno_QCD_70to100" )));
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_70to100 [pr]->Clone("h_pT_barrel_ctrl_QCD_70to100" )));
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_70to100 [pr]->Clone("h_pT_endcap_ctrl_QCD_70to100" )));
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_70to100 [pr]->Clone("h_pT_barrel_nume_QCD_70to100" )));
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_70to100 [pr]->Clone("h_pT_endcap_nume_QCD_70to100" )));
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_100to500[pr]->Clone("h_pT_barrel_deno_QCD_100to500")));
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_100to500[pr]->Clone("h_pT_endcap_deno_QCD_100to500")));
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_100to500[pr]->Clone("h_pT_barrel_ctrl_QCD_100to500")));
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_100to500[pr]->Clone("h_pT_endcap_ctrl_QCD_100to500")));
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_100to500[pr]->Clone("h_pT_barrel_nume_QCD_100to500")));
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_100to500[pr]->Clone("h_pT_endcap_nume_QCD_100to500")));
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_50to70[pr]);
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_50to70[pr]);
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_50to70[pr]);
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_50to70[pr]);
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_50to70[pr]);
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_50to70[pr]);
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_70to100[pr]);
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_70to100[pr]);
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_70to100[pr]);
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_70to100[pr]);
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_70to100[pr]);
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_70to100[pr]);
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_100to500[pr]);
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_100to500[pr]);
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_100to500[pr]);
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_100to500[pr]);
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_100to500[pr]);
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_100to500[pr]);
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
        TH1D *h_temp[6];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_pT_barrel_deno", h_pT_barrel_data_deno);
            file->GetObject("h_pT_endcap_deno", h_pT_endcap_data_deno);
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            removeNegativeBins(h_pT_barrel_data_deno);
            removeNegativeBins(h_pT_endcap_data_deno);
            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_nume);
        }
        else
        {
            file->GetObject("h_pT_barrel_deno", h_temp[0]);
            file->GetObject("h_pT_endcap_deno", h_temp[1]);
            file->GetObject("h_pT_barrel_ctrl", h_temp[2]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[3]);
            file->GetObject("h_pT_barrel_nume", h_temp[4]);
            file->GetObject("h_pT_endcap_nume", h_temp[5]);
            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            removeNegativeBins(h_temp[4]);
            removeNegativeBins(h_temp[5]);
            h_pT_barrel_data_deno->Add(h_temp[0]);
            h_pT_endcap_data_deno->Add(h_temp[1]);
            h_pT_barrel_data_ctrl->Add(h_temp[2]);
            h_pT_endcap_data_ctrl->Add(h_temp[3]);
            h_pT_barrel_data_nume->Add(h_temp[4]);
            h_pT_endcap_data_nume->Add(h_temp[5]);
        }
    }

    h_pT_barrel_data_deno->SetDirectory(0);
    h_pT_endcap_data_deno->SetDirectory(0);
    h_pT_barrel_data_ctrl->SetDirectory(0);
    h_pT_endcap_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume->SetDirectory(0);
    h_pT_endcap_data_nume->SetDirectory(0);

//--------------------------------- FR by ratio --------------------------------------

    // ####### Numerator ####### //
    // Barrel
    h_pT_barrel_nume = ((TH1D*)(h_pT_barrel_MC_deno[_DY_Full]->Clone("h_pT_barrel_nume")));
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_ttbar]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_tW]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_tbarW]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_WW]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_WZ]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_ZZ]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_WJets]);
    h_pT_barrel_nume->Add(h_pT_barrel_MC_deno[_QCDMuEnriched_Full]);
    h_pT_barrel_nume->Multiply(h_pT_barrel_MC_nume[_QCDMuEnriched_Full]);
    h_pT_barrel_nume->Multiply(h_pT_barrel_data_nume);
    // Endcap
    h_pT_endcap_nume = ((TH1D*)(h_pT_endcap_MC_deno[_DY_Full]->Clone("h_pT_endcap_nume")));
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_ttbar]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_tW]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_tbarW]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_WW]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_WZ]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_ZZ]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_WJets]);
    h_pT_endcap_nume->Add(h_pT_endcap_MC_deno[_QCDMuEnriched_Full]);
    h_pT_endcap_nume->Multiply(h_pT_endcap_MC_nume[_QCDMuEnriched_Full]);
    h_pT_endcap_nume->Multiply(h_pT_endcap_data_nume);

    // ####### Denominator ####### //
    // Barrel
    h_pT_barrel_deno = ((TH1D*)(h_pT_barrel_MC_nume[_DY_Full]->Clone("h_pT_barrel_deno")));
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_ttbar]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_tW]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_tbarW]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_WW]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_WZ]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_ZZ]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_WJets]);
    h_pT_barrel_deno->Add(h_pT_barrel_MC_nume[_QCDMuEnriched_Full]);
    h_pT_barrel_deno->Multiply(h_pT_barrel_MC_deno[_QCDMuEnriched_Full]);
    h_pT_barrel_deno->Multiply(h_pT_barrel_data_deno);
    // Endcap
    h_pT_endcap_deno = ((TH1D*)(h_pT_endcap_MC_nume[_DY_Full]->Clone("h_pT_endcap_deno")));
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_ttbar]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_tW]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_tbarW]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_WW]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_WZ]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_ZZ]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_WJets]);
    h_pT_endcap_deno->Add(h_pT_endcap_MC_nume[_QCDMuEnriched_Full]);
    h_pT_endcap_deno->Multiply(h_pT_endcap_MC_deno[_QCDMuEnriched_Full]);
    h_pT_endcap_deno->Multiply(h_pT_endcap_data_deno);

    // ######## FR ######## //
    // Barrel
    h_FRratio_barrel = ((TH1D*)(h_pT_barrel_nume->Clone("h_FRratio_barrel")));
    h_FRratio_barrel->Divide(h_pT_barrel_deno);
    h_FRratio_barrel->SetDirectory(0);
    // Endcap
    h_FRratio_endcap = ((TH1D*)(h_pT_endcap_nume->Clone("h_FRratio_endcap")));
    h_FRratio_endcap->Divide(h_pT_endcap_deno);
    h_FRratio_endcap->SetDirectory(0);



//--------------------------------- FR by template --------------------------------------  

    // Barrel
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.0393e+06/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // OLD All muons' pT>52GeV
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.5567e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.7123e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(4.4126e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(6.0141e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(9.3704e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.0467e+06/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // All muons' pT>52GeV
    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.5551e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.7505e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(4.4041e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(5.9729e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(9.3823e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.1587e+06/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // Only one muon's pT>52GeV
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.7258e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.9382e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(4.4139e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(5.9937e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(9.3962e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());


//    Double_t nevt=0, err=0;
    TH1D *h_FRtemplate_barrel_deno = ((TH1D*)(h_pT_barrel_template_deno_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_barrel_deno")));
    h_FRtemplate_barrel_deno->Add(h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_barrel_deno->Add(h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]);
//    nevt = h_FRtemplate_barrel_deno->IntegralAndError(1, h_FRtemplate_barrel_deno->GetSize()-2, err);
//    cout << "Deno barrel events: " << nevt << "+-" << err << endl;

    h_FRtemplate_barrel = ((TH1D*)(h_pT_barrel_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_barrel")));
    h_FRtemplate_barrel->Add(h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_barrel->Add(h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]);

    h_FRtemplate_barrel->Divide(h_FRtemplate_barrel_deno);
    h_FRtemplate_barrel->SetDirectory(0);

    // For systematic errors
    TH1D *h_FRtemplate_barrel_up   = ((TH1D*)(h_FRtemplate_barrel->Clone("h_FRtemplate_barrel_up")));
    TH1D *h_FRtemplate_barrel_down = ((TH1D*)(h_FRtemplate_barrel->Clone("h_FRtemplate_barrel_down")));
    for (Int_t i=1; i<h_FRtemplate_barrel->GetSize()-1; i++)
    {
        h_FRtemplate_barrel_up  ->SetBinContent(i, h_FRtemplate_barrel->GetBinContent(i)+h_FRtemplate_barrel->GetBinError(i));
        h_FRtemplate_barrel_down->SetBinContent(i, h_FRtemplate_barrel->GetBinContent(i)-h_FRtemplate_barrel->GetBinError(i));
    }

    // Endcap
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.1732e+06/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // OLD All muons' pT>52GeV
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(7.7612e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.9061e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.9792e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(6.0998e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(3.4973e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.1835e+06/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // All muons' pT>52GeV
    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(7.7520e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.9318e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.9852e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(6.0587e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(3.4908e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.4732e+06/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // Only one muon's pT>52GeV
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(9.3080e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.9061e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.9804e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(6.0817e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(3.5033e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());

    TH1D *h_FRtemplate_endcap_deno = ((TH1D*)(h_pT_endcap_template_deno_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_endcap_deno")));
    h_FRtemplate_endcap_deno->Add(h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_endcap_deno->Add(h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]);

    h_FRtemplate_endcap = ((TH1D*)(h_pT_endcap_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_endcap")));
    h_FRtemplate_endcap->Add(h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_endcap->Add(h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]);

    h_FRtemplate_endcap->Divide(h_FRtemplate_endcap_deno);
    h_FRtemplate_endcap->SetDirectory(0);

    // For systematic errors
    TH1D *h_FRtemplate_endcap_up   = ((TH1D*)(h_FRtemplate_endcap->Clone("h_FRtemplate_endcap_up")));
    TH1D *h_FRtemplate_endcap_down = ((TH1D*)(h_FRtemplate_endcap->Clone("h_FRtemplate_endcap_down")));
    for (Int_t i=1; i<h_FRtemplate_endcap->GetSize()-1; i++)
    {
        h_FRtemplate_endcap_up  ->SetBinContent(i, h_FRtemplate_endcap->GetBinContent(i)+h_FRtemplate_endcap->GetBinError(i));
        h_FRtemplate_endcap_down->SetBinContent(i, h_FRtemplate_endcap->GetBinContent(i)-h_FRtemplate_endcap->GetBinError(i));
    }


//--------------------------------- Mixed FR --------------------------------------
    // Numerator -- from ratio method, denominator -- from template fitting

    // ----- Numerator ----- //
    // Barrel
    TH1D * h_FRmixed_barrel = ((TH1D*)(h_pT_barrel_data_nume->Clone("h_FRmixed_barrel"))); // so far it is just a numerator of numerator
    h_FRmixed_barrel->SetDirectory(0);
    h_FRmixed_barrel->Multiply(h_pT_barrel_MC_nume[_QCDMuEnriched_Full]);

    TH1D * h_pT_barrel_mixed_nume_2 = ((TH1D*)(h_pT_barrel_MC_nume[_DY_Full]->Clone("h_pT_barrel_deno"))); // denominator of numerator
    h_pT_barrel_mixed_nume_2->SetDirectory(0);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_ttbar]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_tW]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_tbarW]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_WW]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_WZ]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_ZZ]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_WJets]);
    h_pT_barrel_mixed_nume_2->Add(h_pT_barrel_MC_nume[_QCDMuEnriched_Full]);

    h_FRmixed_barrel->Divide(h_pT_barrel_mixed_nume_2);

    // Endcap
    // Barrel
    TH1D * h_FRmixed_endcap = ((TH1D*)(h_pT_endcap_data_nume->Clone("h_FRmixed_endcap"))); // so far it is just a numerator of numerator
    h_FRmixed_endcap->SetDirectory(0);
    h_FRmixed_endcap->Multiply(h_pT_endcap_MC_nume[_QCDMuEnriched_Full]);

    TH1D * h_pT_endcap_mixed_nume_2 = ((TH1D*)(h_pT_endcap_MC_nume[_DY_Full]->Clone("h_pT_endcap_deno"))); // denominator of numerator
    h_pT_endcap_mixed_nume_2->SetDirectory(0);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_ttbar]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_tW]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_tbarW]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_WW]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_WZ]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_ZZ]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_WJets]);
    h_pT_endcap_mixed_nume_2->Add(h_pT_endcap_MC_nume[_QCDMuEnriched_Full]);

    h_FRmixed_endcap->Divide(h_pT_endcap_mixed_nume_2);

    // ----- Denominator ----- //
    TH1D* h_pT_barrel_mixed_deno = ((TH1D*)(h_FRtemplate_barrel_deno->Clone("h_pT_barrel_mixed_deno")));
    TH1D* h_pT_endcap_mixed_deno = ((TH1D*)(h_FRtemplate_endcap_deno->Clone("h_pT_endcap_mixed_deno")));
    h_pT_barrel_mixed_deno->SetDirectory(0);
    h_pT_endcap_mixed_deno->SetDirectory(0);

    // ----- Fake Rate ----- //
    h_FRmixed_barrel->Divide(h_pT_barrel_mixed_deno);
    h_FRmixed_endcap->Divide(h_pT_endcap_mixed_deno);


//-------------------------- FR from signal+control fit ------------------------------
    // Signal and control regions -- from template fit
    // Numerator = Signal, Denominator = Signal + Control
    Double_t err_signal_barrel[18], err_signal_endcap[9], err_control_barrel[18], err_control_endcap[9];
    Double_t val_signal_barrel[18], val_signal_endcap[9], val_control_barrel[18], val_control_endcap[9];

    // ----- Numerator ----- //
    TH1D *h_FRsigCtrl_template_barrel = ((TH1D*)(h_pT_barrel_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRsigCtrl_template_barrel")));
    h_FRsigCtrl_template_barrel->Add(h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRsigCtrl_template_barrel->Add(h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]);
    TH1D *h_FRsigCtrl_template_endcap = ((TH1D*)(h_pT_endcap_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRsigCtrl_template_endcap")));
    h_FRsigCtrl_template_endcap->Add(h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRsigCtrl_template_endcap->Add(h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]);

    // ---- Denominator ---- //
    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.4763e+07/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(4.1953e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(8.8544e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(6.8025e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.7406e+06/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(2.9938e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());

    TH1D *h_DenoSigCtrl_template_barrel = ((TH1D*)(h_pT_barrel_template_ctrl_50to70[_QCDMuEnriched_Full]->Clone("h_DenoSigCtrl_template_barrel")));
    h_DenoSigCtrl_template_barrel->Add(h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]);
    h_DenoSigCtrl_template_barrel->Add(h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]);
    TH1D *h_DenoSigCtrl_template_endcap = ((TH1D*)(h_pT_endcap_template_ctrl_50to70[_QCDMuEnriched_Full]->Clone("h_DenoSigCtrl_template_endcap")));
    h_DenoSigCtrl_template_endcap->Add(h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]);
    h_DenoSigCtrl_template_endcap->Add(h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]);

    // Getting the errrors
    for (Int_t i_bin=0; i_bin<18; i_bin++)
    {
        err_signal_barrel[i_bin] = sqrt(h_FRsigCtrl_template_barrel->GetBinContent(i_bin+1));
        err_control_barrel[i_bin] = sqrt(h_DenoSigCtrl_template_barrel->GetBinContent(i_bin+1));
        val_signal_barrel[i_bin] = h_FRsigCtrl_template_barrel->GetBinContent(i_bin+1);
        val_control_barrel[i_bin] = h_DenoSigCtrl_template_barrel->GetBinContent(i_bin+1);
        if (i_bin < 9)
        {
            err_signal_endcap[i_bin] = sqrt(h_FRsigCtrl_template_endcap->GetBinContent(i_bin+1));
            err_control_endcap[i_bin] = sqrt(h_DenoSigCtrl_template_endcap->GetBinContent(i_bin+1));
            val_signal_endcap[i_bin] = h_FRsigCtrl_template_endcap->GetBinContent(i_bin+1);
            val_control_endcap[i_bin] = h_DenoSigCtrl_template_endcap->GetBinContent(i_bin+1);
        }
    }

    h_DenoSigCtrl_template_barrel->Add(h_FRsigCtrl_template_barrel); // Deno = Signal+Control
    h_DenoSigCtrl_template_endcap->Add(h_FRsigCtrl_template_endcap);

    // ----- Fake rate ----- //
    h_FRsigCtrl_template_barrel->Divide(h_DenoSigCtrl_template_barrel);
    h_FRsigCtrl_template_endcap->Divide(h_DenoSigCtrl_template_endcap);

    // ------- Errors ------- //
    // Delta(A/A+B) = AB/(A+B)^2 * SQRT[(DeltaA/A)^2+(DeltaB/B)^2]
    // Where A -- signal, B -- control, and their Deltas are their square roots
    for (Int_t i_bin=0; i_bin<18; i_bin++)
    {
        Double_t err_barrel;
        err_barrel = (val_signal_barrel[i_bin] * val_control_barrel[i_bin]) /
                    ((val_signal_barrel[i_bin] + val_control_barrel[i_bin]) * (val_signal_barrel[i_bin] + val_control_barrel[i_bin])) *
                     sqrt(1/val_signal_barrel[i_bin] + 1/val_control_barrel[i_bin]);
        h_FRsigCtrl_template_barrel->SetBinError(i_bin+1, err_barrel);
        if (i_bin < 9)
        {
            Double_t err_endcap;
            err_endcap = (val_signal_endcap[i_bin] * val_control_endcap[i_bin]) /
                        ((val_signal_endcap[i_bin] + val_control_endcap[i_bin]) * (val_signal_endcap[i_bin] + val_control_endcap[i_bin])) *
                         sqrt(1/val_signal_endcap[i_bin] + 1/val_control_endcap[i_bin]);
            h_FRsigCtrl_template_endcap->SetBinError(i_bin+1, err_endcap);
        }
    }

//--------------------------------- Dalmin's FR --------------------------------------
    h_FRdalmin_barrel = ((TH1D*)(h_FRmixed_barrel->Clone("h_FRdalmin_barrel")));
    h_FRdalmin_barrel->SetDirectory(0);
    h_FRdalmin_endcap = ((TH1D*)(h_FRmixed_endcap->Clone("h_FRdalmin_endcap")));
    h_FRdalmin_endcap->SetDirectory(0);

    h_FRdalmin_barrel->SetBinContent(1, 0.073598);
    h_FRdalmin_barrel->SetBinContent(2, 0.0748619);
    h_FRdalmin_barrel->SetBinContent(3, 0.0735089);
    h_FRdalmin_barrel->SetBinContent(4, 0.0835437);
    h_FRdalmin_barrel->SetBinContent(5, 0.0837622);
    h_FRdalmin_barrel->SetBinContent(6, 0.0810606);
    h_FRdalmin_barrel->SetBinContent(7, 0.0819881);
    h_FRdalmin_barrel->SetBinContent(8, 0.0986144);
    h_FRdalmin_barrel->SetBinContent(9, 0.0899716);
    h_FRdalmin_barrel->SetBinContent(10, 0.114714);
    h_FRdalmin_barrel->SetBinContent(11, 0.114298);
    h_FRdalmin_barrel->SetBinContent(12, 0.138463);
    h_FRdalmin_barrel->SetBinContent(13, 0.170341);
    h_FRdalmin_barrel->SetBinContent(14, 0.188865);
    h_FRdalmin_barrel->SetBinContent(15, 0.244211);
    h_FRdalmin_barrel->SetBinContent(16, 0.263503);

    h_FRdalmin_endcap->SetBinContent(1, 0.180487);
    h_FRdalmin_endcap->SetBinContent(2, 0.187814);
    h_FRdalmin_endcap->SetBinContent(3, 0.205887);
    h_FRdalmin_endcap->SetBinContent(4, 0.193403);
    h_FRdalmin_endcap->SetBinContent(5, 0.209049);
    h_FRdalmin_endcap->SetBinContent(6, 0.220315);
    h_FRdalmin_endcap->SetBinContent(7, 0.307998);
    h_FRdalmin_endcap->SetBinContent(8, 0.395422);


    // Writing
    TFile *file_FR = new TFile("/media/sf_DATA/FR/Muon/FakeRate_muon.root", "RECREATE");
    if (file_FR->IsOpen()) cout << "File '/media/sf_DATA/FR/Muon/FakeRate_muon.root' has been created. Writing histograms.." << endl;
    file_FR->cd();
    h_FRratio_barrel->Write();
    h_FRratio_endcap->Write();
    h_FRtemplate_barrel->Write();
    h_FRtemplate_endcap->Write();
    h_FRtemplate_barrel_up->Write();
    h_FRtemplate_endcap_up->Write();
    h_FRtemplate_barrel_down->Write();
    h_FRtemplate_endcap_down->Write();
    h_FRmixed_barrel->Write();
    h_FRmixed_endcap->Write();
    h_FRsigCtrl_template_barrel->Write();
    h_FRsigCtrl_template_endcap->Write();
    h_FRdalmin_barrel->Write();
    h_FRdalmin_endcap->Write();
    cout << "Finished. Closing the file.." << endl;
    file_FR->Close();
    if (!file_FR->IsOpen()) cout << "File '/media/sf_DATA/FR/Muon/FakeRate_muon.root' has been closed successfully." << endl;
    else cout << "File did not close!" << endl;

    // Drawing
    TCanvas *c_FR_barrel = new TCanvas("c_FR_barrel", "c_FR_barrel", 800, 800);
    c_FR_barrel->cd();
    c_FR_barrel->SetGrid(1);
    c_FR_barrel->SetLogx(1);
    c_FR_barrel->SetRightMargin(0.05);
    c_FR_barrel->SetTopMargin(0.05);
    c_FR_barrel->SetBottomMargin(0.12);
    c_FR_barrel->SetLeftMargin(0.13);
    h_FRratio_barrel->SetMarkerStyle(kFullSquare);
    h_FRratio_barrel->SetMarkerColor(kRed);
    h_FRratio_barrel->SetLineColor(kRed);
    h_FRratio_barrel->SetStats(kFALSE);
    h_FRratio_barrel->SetTitle("");
    h_FRratio_barrel->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FRratio_barrel->GetXaxis()->SetTitleOffset(1);
    h_FRratio_barrel->GetXaxis()->SetTitleSize(0.05);
    h_FRratio_barrel->GetXaxis()->SetLabelSize(0.04);
    h_FRratio_barrel->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_barrel->GetYaxis()->SetTitleSize(0.05);
    h_FRratio_barrel->GetYaxis()->SetTitleOffset(1.25);
    h_FRratio_barrel->GetYaxis()->SetLabelSize(0.04);
    h_FRratio_barrel->GetXaxis()->SetNoExponent(1);
    h_FRratio_barrel->GetXaxis()->SetMoreLogLabels(1);
    h_FRratio_barrel->GetXaxis()->SetRangeUser(52, 1000);
    h_FRratio_barrel->GetYaxis()->SetRangeUser(0, 0.55);
    h_FRratio_barrel->Draw();
    h_FRmixed_barrel->SetMarkerStyle(22);
    h_FRmixed_barrel->SetMarkerColor(kBlue);
    h_FRmixed_barrel->SetLineColor(kBlue);
    h_FRmixed_barrel->SetStats(kFALSE);
    h_FRmixed_barrel->Draw("same");
    h_FRtemplate_barrel->SetMarkerStyle(33);
    h_FRtemplate_barrel->SetMarkerColor(kGreen+2);
    h_FRtemplate_barrel->SetLineColor(kGreen+2);
    h_FRtemplate_barrel->SetStats(kFALSE);
    h_FRtemplate_barrel->Draw("same");
    h_FRsigCtrl_template_barrel->SetMarkerStyle(kFullDotLarge);
    h_FRsigCtrl_template_barrel->SetMarkerColor(kBlack);
    h_FRsigCtrl_template_barrel->SetLineColor(kBlack);
    h_FRsigCtrl_template_barrel->SetStats(kFALSE);
    h_FRsigCtrl_template_barrel->Draw("same");

    TLegend *legend = new TLegend(0.13, 0.77, 0.6, 0.95);
    legend->AddEntry(h_FRratio_barrel, "Ratio", "LP");
    legend->AddEntry(h_FRmixed_barrel, "Template (deno)", "LP");
    legend->AddEntry(h_FRtemplate_barrel, "Template (nume, deno)", "LP");
    legend->AddEntry(h_FRsigCtrl_template_barrel, "Template (signal, non-signal)", "LP");
    legend->Draw();
    legend->Draw();
    TText *textb = new TText (0.45, 0.6, "Barrel");
    textb->SetTextAlign(11);
    textb->SetTextSize(0.05);
    textb->SetNDC(true);
    textb->Draw();
    c_FR_barrel->Update();

    TCanvas *c_FR_endcap = new TCanvas("c_FR_endcap", "c_FR_endcap", 800, 800);
    c_FR_endcap->cd();
    c_FR_endcap->cd();
    c_FR_endcap->SetGrid(1);
    c_FR_endcap->SetLogx(1);
    c_FR_endcap->SetRightMargin(0.05);
    c_FR_endcap->SetTopMargin(0.05);
    c_FR_endcap->SetBottomMargin(0.12);
    c_FR_endcap->SetLeftMargin(0.13);
    h_FRratio_endcap->SetMarkerStyle(kFullSquare);
    h_FRratio_endcap->SetMarkerColor(kRed);
    h_FRratio_endcap->SetLineColor(kRed);
    h_FRratio_endcap->SetStats(kFALSE);
    h_FRratio_endcap->SetTitle("");
    h_FRratio_endcap->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FRratio_endcap->GetXaxis()->SetTitleOffset(1);
    h_FRratio_endcap->GetXaxis()->SetTitleSize(0.05);
    h_FRratio_endcap->GetXaxis()->SetLabelSize(0.04);
    h_FRratio_endcap->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_endcap->GetYaxis()->SetTitleSize(0.05);
    h_FRratio_endcap->GetYaxis()->SetTitleOffset(1.25);
    h_FRratio_endcap->GetYaxis()->SetLabelSize(0.04);
    h_FRratio_endcap->GetXaxis()->SetNoExponent(1);
    h_FRratio_endcap->GetXaxis()->SetMoreLogLabels(1);
    h_FRratio_endcap->GetXaxis()->SetRangeUser(52, 1000);
    h_FRratio_endcap->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_endcap->Draw();
    h_FRmixed_endcap->SetMarkerStyle(22);
    h_FRmixed_endcap->SetMarkerColor(kBlue);
    h_FRmixed_endcap->SetLineColor(kBlue);
    h_FRmixed_endcap->SetStats(kFALSE);
    h_FRmixed_endcap->Draw("same");
    h_FRtemplate_endcap->SetMarkerStyle(33);
    h_FRtemplate_endcap->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap->SetLineColor(kGreen+2);
    h_FRtemplate_endcap->SetStats(kFALSE);
    h_FRtemplate_endcap->Draw("same");
    h_FRsigCtrl_template_endcap->SetMarkerStyle(kFullDotLarge);
    h_FRsigCtrl_template_endcap->SetMarkerColor(kBlack);
    h_FRsigCtrl_template_endcap->SetLineColor(kBlack);
    h_FRsigCtrl_template_endcap->SetStats(kFALSE);
    h_FRsigCtrl_template_endcap->Draw("same");
    legend->Draw();
    TText *texte = new TText (0.45, 0.6, "Endcap");
    texte->SetTextAlign(11);
    texte->SetTextSize(0.05);
    texte->SetNDC(true);
    texte->Draw();
    c_FR_endcap->Update();

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
