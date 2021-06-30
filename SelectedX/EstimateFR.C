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
#include <TText.h>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"
#include "./header/myRatioPlot_t.h"
#include "./etc/RoccoR/RoccoR.cc"

void E_EstFR (Int_t type);
void E_EstFR_alt();
void E_EstFR_alt2();
void E_EstPR(); // PROMPT RATE
void E_EstPR_alt(); // PROMPT RATE
void E_EstFRandPR_MC();
void Mu_EstFR (Int_t type);
void Mu_EstPR();
void Mu_EstPR_alt();
void Mu_EstFRandPR_MC();

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
        if (whichX.Contains("MC"))
        {
            cout << "\n*******      E_EstFRandPR_MC()      *******" << endl;
            E_EstFRandPR_MC();
        }
        else if (whichX.Contains("PR"))
        {
            if (whichX.Contains("ALT"))
            {
                cout << "\n*******      E_EstPR_alt()      *******" << endl;
                E_EstPR_alt();
            }
            else
            {
                cout << "\n*******      E_EstPR()      *******" << endl;
                E_EstPR();
            }
        }
        else if (whichX.Contains("ALT") && whichX.Contains("2"))
        {
            cout << "\n*******      E_EstFR_alt2()      *******" << endl;
            E_EstFR_alt2();
        }
        else if (whichX.Contains("ALT"))
        {
            cout << "\n*******      E_EstFR_alt()      *******" << endl;
            E_EstFR_alt();
        }
        else
        {
            cout << "\n*******      E_EstFR(" << type << ")      *******" << endl;
            E_EstFR(type);
        }
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("MC"))
        {
            cout << "\n*******      Mu_EstFRandPR_MC()      *******" << endl;
            Mu_EstFRandPR_MC();
        }
        else if (whichX.Contains("PR"))
        {
            if (whichX.Contains("ALT"))
            {
                cout << "\n*******      Mu_EstPR_alt()      *******" << endl;
                Mu_EstPR_alt();
            }
            else
            {
                cout << "\n*******      Mu_EstPR()      *******" << endl;
                Mu_EstPR();
            }
        }
        else
        {
            cout << "\n*******     Mu_EstFR(" << type << ")     *******" << endl;
            Mu_EstFR(type);
        }
    }
    if (Xselected == 0) cout << "Wrong arument!" << endl;

} // End of HistDrawer()


/// ############################################################################# ///
/// ----------------------------- Electron Channel ------------------------------ ///
/// ############################################################################# ///
void E_EstFR(Int_t type)
{
    FileMgr fm;

//    TH1D *h_pT_barrel_nume,      *h_pT_endcap_nume,      *h_pT_barrel_deno,      *h_pT_endcap_deno,
//         *h_pT_barrel_nume_sub,  *h_pT_endcap_nume_sub,  *h_pT_barrel_deno_sub,  *h_pT_endcap_deno_sub,
//         *h_pT_barrel_nume_abcd, *h_pT_endcap_nume_abcd, *h_pT_barrel_deno_abcd, *h_pT_endcap_deno_abcd,
//         *h_FRratio_barrel,    *h_FRratio_endcap,
//         *h_FRsubtract_barrel, *h_FRsubtract_endcap,
//         *h_FRtemplate_barrel, *h_FRtemplate_endcap;

//    TH1D *h_eta_nume,     *h_eta_deno,
//         *h_eta_nume_sub, *h_eta_deno_sub,
//         *h_FRsubtract_eta, *h_FRratio_eta;

    TH1D *h_pT_barrel_MC_nume[_EndOf_Data_Special],     *h_pT_endcap_MC_nume[_EndOf_Data_Special], *h_pT_endcap2_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special],     *h_pT_endcap_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap2_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_data_nume,                        *h_pT_endcap_data_nume,                    *h_pT_endcap2_data_nume,
         *h_pT_barrel_data_ctrl,                        *h_pT_endcap_data_ctrl,                    *h_pT_endcap2_data_ctrl,
         *h_eta_MC_nume[_EndOf_Data_Special], *h_eta_MC_ctrl[_EndOf_Data_Special], *h_eta_data_nume, *h_eta_data_ctrl;

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr1]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr1]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr1]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr1]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr1]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr1]);
        removeNegativeBins(h_eta_MC_ctrl[pr1]);
        removeNegativeBins(h_eta_MC_nume[pr1]);

        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr1]->SetDirectory(0);
        h_eta_MC_ctrl[pr1]->SetDirectory(0);
        h_eta_MC_nume[pr1]->SetDirectory(0);

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

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_DY")));
            h_pT_endcap2_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap2_MC_ctrl[pr]->Clone("h_pT_endcap2_MC_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_DY")));
            h_pT_endcap2_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap2_MC_nume[pr]->Clone("h_pT_endcap2_MC_nume_DY")));
            h_eta_MC_ctrl[_DY_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_DY")));
            h_eta_MC_nume[_DY_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap2_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap2_MC_nume[_DY_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_MC_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_endcap2_MC_ctrl[_DY_Full]->Add(h_pT_endcap2_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_endcap2_MC_nume[_DY_Full]->Add(h_pT_endcap2_MC_nume[pr]);
            h_eta_MC_ctrl[_DY_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_DY_Full]->Add(h_eta_MC_nume[pr]);
        }
        file->Close();
    }

    // Gamma+Jets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _GJets_20to100)
        {
            h_pT_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_GJets")));
            h_pT_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_GJets")));
            h_pT_endcap2_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap2_MC_ctrl[pr]->Clone("h_pT_endcap2_MC_ctrl_GJets")));
            h_pT_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_GJets")));
            h_pT_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_GJets")));
            h_pT_endcap2_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap2_MC_nume[pr]->Clone("h_pT_endcap2_MC_nume_GJets")));
            h_eta_MC_ctrl[_GJets_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_GJets")));
            h_eta_MC_nume[_GJets_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_GJets")));

            h_pT_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap2_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap2_MC_nume[_GJets_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_eta_MC_nume[_GJets_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_GJets_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_GJets_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_endcap2_MC_ctrl[_GJets_Full]->Add(h_pT_endcap2_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_GJets_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_GJets_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_endcap2_MC_nume[_GJets_Full]->Add(h_pT_endcap2_MC_nume[pr]);
            h_eta_MC_ctrl[_GJets_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_GJets_Full]->Add(h_eta_MC_nume[pr]);
        }
        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_QCD")));
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_QCD")));
            h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap2_MC_ctrl[pr]->Clone("h_pT_endcap2_MC_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_QCD")));
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_QCD")));
            h_pT_endcap2_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap2_MC_nume[pr]->Clone("h_pT_endcap2_MC_nume_QCD")));
            h_eta_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_QCD")));
            h_eta_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_QCD")));

            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap2_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap2_MC_nume[pr]);
            h_eta_MC_ctrl[_QCDEMEnriched_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_QCDEMEnriched_Full]->Add(h_eta_MC_nume[pr]);
        }
        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        TH1D *h_temp[8];
        if (pr == _SinglePhoton_B)
        {
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_data_nume);
            file->GetObject("h_eta_ctrl", h_eta_data_ctrl);
            file->GetObject("h_eta_nume", h_eta_data_nume);

            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_pT_endcap2_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_nume);
            removeNegativeBins(h_pT_endcap2_data_nume);
            removeNegativeBins(h_eta_data_ctrl);
            removeNegativeBins(h_eta_data_nume);
        }
        else
        {
            file->GetObject("h_pT_barrel_ctrl",  h_temp[0]);
            file->GetObject("h_pT_endcap_ctrl",  h_temp[1]);
            file->GetObject("h_pT_endcap2_ctrl", h_temp[2]);
            file->GetObject("h_pT_barrel_nume",  h_temp[3]);
            file->GetObject("h_pT_endcap_nume",  h_temp[4]);
            file->GetObject("h_pT_endcap2_nume", h_temp[5]);
            file->GetObject("h_eta_ctrl",        h_temp[6]);
            file->GetObject("h_eta_nume",        h_temp[7]);

            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            removeNegativeBins(h_temp[4]);
            removeNegativeBins(h_temp[5]);
            removeNegativeBins(h_temp[6]);
            removeNegativeBins(h_temp[7]);

            h_pT_barrel_data_ctrl ->Add(h_temp[0]);
            h_pT_endcap_data_ctrl ->Add(h_temp[1]);
            h_pT_endcap2_data_ctrl->Add(h_temp[2]);
            h_pT_barrel_data_nume ->Add(h_temp[3]);
            h_pT_endcap_data_nume ->Add(h_temp[4]);
            h_pT_endcap2_data_nume->Add(h_temp[5]);
            h_eta_data_ctrl       ->Add(h_temp[6]);
            h_eta_data_nume       ->Add(h_temp[7]);
        }
    }

    h_pT_barrel_data_ctrl ->SetDirectory(0);
    h_pT_endcap_data_ctrl ->SetDirectory(0);
    h_pT_endcap2_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume ->SetDirectory(0);
    h_pT_endcap_data_nume ->SetDirectory(0);
    h_pT_endcap2_data_nume->SetDirectory(0);
    h_eta_data_ctrl       ->SetDirectory(0);
    h_eta_data_nume       ->SetDirectory(0);


//--------------------------------- FR from QCD MC -------------------------------------- (deno = nume + ctrl)

    //             QCD_MC_nume
    // FR = -------------------------
    //      QCD_MC_nume + QCD_MC_ctrl

    // ------ Numerator ------ //
    TH1D *h_pT_barrel_nume_fMC = ((TH1D*)(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Clone("h_pT_barrel_nume_fMC")));
    TH1D *h_pT_endcap_nume_fMC = ((TH1D*)(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Clone("h_pT_endcap_nume_fMC")));
    TH1D *h_pT_endcap2_nume_fMC = ((TH1D*)(h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]->Clone("h_pT_endcap2_nume_fMC")));
    TH1D *h_eta_nume_fMC = ((TH1D*)(h_eta_MC_nume[_QCDEMEnriched_Full]->Clone("h_eta_nume_fMC")));

    // ------ Denominator ------ //
    TH1D *h_pT_barrel_deno_fMC = ((TH1D*)(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_barrel_deno_fMC")));
    h_pT_barrel_deno_fMC->Add(h_pT_barrel_nume_fMC); // deno = sig+ctrl
    TH1D *h_pT_endcap_deno_fMC = ((TH1D*)(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_endcap_deno_fMC")));
    h_pT_endcap_deno_fMC->Add(h_pT_endcap_nume_fMC); // deno = sig+ctrl
    TH1D *h_pT_endcap2_deno_fMC = ((TH1D*)(h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_endcap2_deno_fMC")));
    h_pT_endcap2_deno_fMC->Add(h_pT_endcap2_nume_fMC); // deno = sig+ctrl
    TH1D *h_eta_deno_fMC = ((TH1D*)(h_eta_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_eta_deno_fMC")));
    h_eta_deno_fMC->Add(h_eta_nume_fMC); // deno = sig+ctrl

    // ------ FR ------ //
    // Barrel
    TH1D *h_FRMC_barrel = ((TH1D*)(h_pT_barrel_nume_fMC->Clone("h_FRMC_barrel")));
    h_FRMC_barrel->Divide(h_pT_barrel_deno_fMC);
    h_FRMC_barrel->SetDirectory(0);
    cout << "Numerator barrel (QCD MC): " << h_pT_barrel_nume_fMC->Integral() << endl;
    cout << "Denominator barrel (QCD MC): " << h_pT_barrel_deno_fMC->Integral() << endl;
    // Endcap
    TH1D *h_FRMC_endcap = ((TH1D*)(h_pT_endcap_nume_fMC->Clone("h_FRMC_endcap")));
    h_FRMC_endcap->Divide(h_pT_endcap_deno_fMC);
    h_FRMC_endcap->SetDirectory(0);
    cout << "Numerator endcap (QCD MC): " << h_pT_endcap_nume_fMC->Integral() << endl;
    cout << "Denominator endcap (QCD MC): " << h_pT_endcap_deno_fMC->Integral() << endl;
    // Far endcap
    TH1D *h_FRMC_endcap2 = ((TH1D*)(h_pT_endcap2_nume_fMC->Clone("h_FRMC_endcap2")));
    h_FRMC_endcap2->Divide(h_pT_endcap2_deno_fMC);
    h_FRMC_endcap2->SetDirectory(0);
    cout << "Numerator far endcap (QCD MC): " << h_pT_endcap2_nume_fMC->Integral() << endl;
    cout << "Denominator far endcap (QCD MC): " << h_pT_endcap2_deno_fMC->Integral() << endl;
    // Eta
    TH1D *h_FRMC_eta = ((TH1D*)(h_eta_nume_fMC->Clone("h_FRMC_eta")));
    h_FRMC_eta->Divide(h_eta_deno_fMC);
    h_FRMC_eta->SetDirectory(0);

//--------------------------------- FR by ratio -------------------------------------- (deno = nume + ctrl)

    //            DATA_nume * QCD_nume * sum(allMC_nume + allMC_ctrl)
    // FR = -----------------------------------------------------------------
    //      (DATA_nume + DATA_ctrl) * (QCD_nume + QCD_ctrl) * sum(allMC_nume)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_div = ((TH1D*)(h_pT_barrel_MC_nume[_DY_Full]->Clone("h_pT_barrel_nume_div")));
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_ttbar]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_tW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_tbarW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_WW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_WZ]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_ZZ]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_WJets]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_GJets_Full]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_nume = ((TH1D*)(h_pT_barrel_data_nume->Clone("h_pT_barrel_nume")));
    h_pT_barrel_nume->Multiply(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]);
    h_pT_barrel_nume->Divide(h_pT_barrel_nume_div);

    // Endcap
    TH1D *h_pT_endcap_nume_div = ((TH1D*)(h_pT_endcap_MC_nume[_DY_Full]->Clone("h_pT_endcap_nume_div")));
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_ttbar]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_tW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_tbarW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_WW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_WZ]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_ZZ]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_WJets]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_GJets_Full]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_nume = ((TH1D*)(h_pT_endcap_data_nume->Clone("h_pT_endcap_nume")));
    h_pT_endcap_nume->Multiply(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]);
    h_pT_endcap_nume->Divide(h_pT_endcap_nume_div);

    // Far endcap
    TH1D *h_pT_endcap2_nume_div = ((TH1D*)(h_pT_endcap2_MC_nume[_DY_Full]->Clone("h_pT_endcap2_nume_div")));
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_ttbar]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_tW]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_tbarW]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_WW]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_WZ]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_ZZ]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_WJets]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_GJets_Full]);
    h_pT_endcap2_nume_div->Add(h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap2_nume = ((TH1D*)(h_pT_endcap2_data_nume->Clone("h_pT_endcap2_nume")));
    h_pT_endcap2_nume->Multiply(h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]);
    h_pT_endcap2_nume->Divide(h_pT_endcap2_nume_div);

    // Eta
    TH1D *h_eta_nume_div = ((TH1D*)(h_eta_MC_nume[_DY_Full]->Clone("h_eta_nume_div")));
    h_eta_nume_div->Add(h_eta_MC_nume[_ttbar]);
    h_eta_nume_div->Add(h_eta_MC_nume[_tW]);
    h_eta_nume_div->Add(h_eta_MC_nume[_tbarW]);
    h_eta_nume_div->Add(h_eta_MC_nume[_WW]);
    h_eta_nume_div->Add(h_eta_MC_nume[_WZ]);
    h_eta_nume_div->Add(h_eta_MC_nume[_ZZ]);
    h_eta_nume_div->Add(h_eta_MC_nume[_WJets]);
    h_eta_nume_div->Add(h_eta_MC_nume[_GJets_Full]);
    h_eta_nume_div->Add(h_eta_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_eta_nume = ((TH1D*)(h_eta_data_nume->Clone("h_eta_nume")));
    h_eta_nume->Multiply(h_eta_MC_nume[_QCDEMEnriched_Full]);
    h_eta_nume->Divide(h_eta_nume_div);

    // ####### Control ####### //
    // Barrel
    TH1D *h_pT_barrel_ctrl_div = ((TH1D*)(h_pT_barrel_MC_ctrl[_DY_Full]->Clone("h_pT_barrel_ctrl_div")));
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_ttbar]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_tW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_tbarW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_WW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_WZ]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_ZZ]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_WJets]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_GJets_Full]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_ctrl = ((TH1D*)(h_pT_barrel_data_ctrl->Clone("h_pT_barrel_ctrl")));
    h_pT_barrel_ctrl->Multiply(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_barrel_ctrl->Divide(h_pT_barrel_ctrl_div);

    // Endcap
    TH1D *h_pT_endcap_ctrl_div = ((TH1D*)(h_pT_endcap_MC_ctrl[_DY_Full]->Clone("h_pT_endcap_ctrl_div")));
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_ttbar]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_tW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_tbarW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_WW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_WZ]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_ZZ]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_WJets]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_GJets_Full]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_ctrl = ((TH1D*)(h_pT_endcap_data_ctrl->Clone("h_pT_endcap_ctrl")));
    h_pT_endcap_ctrl->Multiply(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_endcap_ctrl->Divide(h_pT_endcap_ctrl_div);

    // Endcap
    TH1D *h_pT_endcap2_ctrl_div = ((TH1D*)(h_pT_endcap2_MC_ctrl[_DY_Full]->Clone("h_pT_endcap2_ctrl_div")));
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_ttbar]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_tW]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_tbarW]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_WW]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_WZ]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_ZZ]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_WJets]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_GJets_Full]);
    h_pT_endcap2_ctrl_div->Add(h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap2_ctrl = ((TH1D*)(h_pT_endcap2_data_ctrl->Clone("h_pT_endcap2_ctrl")));
    h_pT_endcap2_ctrl->Multiply(h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_endcap2_ctrl->Divide(h_pT_endcap2_ctrl_div);

    // Eta
    TH1D *h_eta_ctrl_div = ((TH1D*)(h_eta_MC_ctrl[_DY_Full]->Clone("h_eta_ctrl_div")));
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_ttbar]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_tW]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_tbarW]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_WW]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_WZ]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_ZZ]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_WJets]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_GJets_Full]);
    h_eta_ctrl_div->Add(h_eta_MC_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_eta_ctrl = ((TH1D*)(h_eta_data_ctrl->Clone("h_eta_ctrl")));
    h_eta_ctrl->Multiply(h_eta_MC_ctrl[_QCDEMEnriched_Full]);
    h_eta_ctrl->Divide(h_eta_ctrl_div);

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_deno_div = ((TH1D*)(h_pT_barrel_MC_ctrl[_DY_Full]->Clone("h_pT_barrel_deno_div")));
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_ttbar]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_tW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_tbarW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_WW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_WZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_ZZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_WJets]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_GJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_DY_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_ttbar]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_tW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_tbarW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_WW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_WZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_ZZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_WJets]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_GJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_deno_mult = ((TH1D*)(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_barrel_deno_mult")));
    h_pT_barrel_deno_mult->Add(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_deno = ((TH1D*)(h_pT_barrel_data_ctrl->Clone("h_pT_barrel_deno")));
    h_pT_barrel_deno->Add(h_pT_barrel_data_nume);
    h_pT_barrel_deno->Multiply(h_pT_barrel_deno_mult);
    h_pT_barrel_deno->Divide(h_pT_barrel_deno_div);

    // Endcap
    TH1D *h_pT_endcap_deno_div = ((TH1D*)(h_pT_endcap_MC_ctrl[_DY_Full]->Clone("h_pT_endcap_deno_div")));
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_ttbar]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_tW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_tbarW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_WW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_WZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_ZZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_WJets]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_GJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_DY_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_ttbar]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_tW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_tbarW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_WW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_WZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_ZZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_WJets]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_GJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_deno_mult = ((TH1D*)(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_endcap_deno_mult")));
    h_pT_endcap_deno_mult->Add(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_deno = ((TH1D*)(h_pT_endcap_data_ctrl->Clone("h_pT_endcap_deno")));
    h_pT_endcap_deno->Add(h_pT_endcap_data_nume);
    h_pT_endcap_deno->Multiply(h_pT_endcap_deno_mult);
    h_pT_endcap_deno->Divide(h_pT_endcap_deno_div);

    // Far endcap
    TH1D *h_pT_endcap2_deno_div = ((TH1D*)(h_pT_endcap2_MC_ctrl[_DY_Full]->Clone("h_pT_endcap2_deno_div")));
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_ttbar]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_tW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_tbarW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_WW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_WZ]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_ZZ]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_WJets]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_GJets_Full]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_DY_Full]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_ttbar]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_tW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_tbarW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_WW]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_WZ]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_ZZ]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_WJets]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_GJets_Full]);
    h_pT_endcap2_deno_div->Add(h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap2_deno_mult = ((TH1D*)(h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_endcap2_deno_mult")));
    h_pT_endcap2_deno_mult->Add(h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap2_deno = ((TH1D*)(h_pT_endcap2_data_ctrl->Clone("h_pT_endcap2_deno")));
    h_pT_endcap2_deno->Add(h_pT_endcap2_data_nume);
    h_pT_endcap2_deno->Multiply(h_pT_endcap2_deno_mult);
    h_pT_endcap2_deno->Divide(h_pT_endcap2_deno_div);

    // Eta
    TH1D *h_eta_deno_div = ((TH1D*)(h_eta_MC_ctrl[_DY_Full]->Clone("h_eta_deno_div")));
    h_eta_deno_div->Add(h_eta_MC_ctrl[_ttbar]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_tW]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_tbarW]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_WW]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_WZ]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_ZZ]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_WJets]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_GJets_Full]);
    h_eta_deno_div->Add(h_eta_MC_ctrl[_QCDEMEnriched_Full]);
    h_eta_deno_div->Add(h_eta_MC_nume[_DY_Full]);
    h_eta_deno_div->Add(h_eta_MC_nume[_ttbar]);
    h_eta_deno_div->Add(h_eta_MC_nume[_tW]);
    h_eta_deno_div->Add(h_eta_MC_nume[_tbarW]);
    h_eta_deno_div->Add(h_eta_MC_nume[_WW]);
    h_eta_deno_div->Add(h_eta_MC_nume[_WZ]);
    h_eta_deno_div->Add(h_eta_MC_nume[_ZZ]);
    h_eta_deno_div->Add(h_eta_MC_nume[_WJets]);
    h_eta_deno_div->Add(h_eta_MC_nume[_GJets_Full]);
    h_eta_deno_div->Add(h_eta_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_eta_deno_mult = ((TH1D*)(h_eta_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_eta_deno_mult")));
    h_eta_deno_mult->Add(h_eta_MC_nume[_QCDEMEnriched_Full]);
    TH1D *h_eta_deno = ((TH1D*)(h_eta_data_ctrl->Clone("h_eta_deno")));
    h_eta_deno->Add(h_eta_data_nume);
    h_eta_deno->Multiply(h_eta_deno_mult);
    h_eta_deno->Divide(h_eta_deno_div);

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRratio_barrel = ((TH1D*)(h_pT_barrel_nume->Clone("h_FRratio_barrel")));
    h_FRratio_barrel->Divide(h_pT_barrel_deno);
    h_FRratio_barrel->SetDirectory(0);
    cout << "Numerator barrel (ratio): " << h_pT_barrel_nume->Integral() << endl;
    cout << "Denominator barrel (ratio): " << h_pT_barrel_deno->Integral() << endl;
    // Endcap
    TH1D *h_FRratio_endcap = ((TH1D*)(h_pT_endcap_nume->Clone("h_FRratio_endcap")));
    h_FRratio_endcap->Divide(h_pT_endcap_deno);
    h_FRratio_endcap->SetDirectory(0);
    cout << "Numerator endcap (ratio): " << h_pT_endcap_nume->Integral() << endl;
    cout << "Denominator endcap (ratio): " << h_pT_endcap_deno->Integral() << endl;
    // Far endcap
    TH1D *h_FRratio_endcap2 = ((TH1D*)(h_pT_endcap2_nume->Clone("h_FRratio_endcap2")));
    h_FRratio_endcap2->Divide(h_pT_endcap2_deno);
    h_FRratio_endcap2->SetDirectory(0);
    cout << "Numerator far endcap (ratio): " << h_pT_endcap2_nume->Integral() << endl;
    cout << "Denominator far endcap (ratio): " << h_pT_endcap2_deno->Integral() << endl;
    // Eta
    TH1D *h_FRratio_eta = ((TH1D*)(h_eta_nume->Clone("h_FRratio_eta")));
    h_FRratio_eta->Divide(h_eta_deno);
    h_FRratio_eta->SetDirectory(0);


//--------------------------------- FR by subtraction -------------------------------------- (deno = nume + ctrl)

    //                     DATA_nume - sum(nonQCD_MC_nume)
    // FR = --------------------------------------------------------------
    //      (DATA_nume + DATA_ctrl) - sum(nonQCD_MC_nume + nonQCD_MC_ctrl)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_sub = ((TH1D*)(h_pT_barrel_data_nume->Clone("h_pT_barrel_nume_sub")));
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_ttbar], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_tW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_tbarW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_WW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_WZ], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_ZZ], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_WJets], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_GJets_Full], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_MC_nume[_DY_Full], -1);
    // Endcap
    TH1D *h_pT_endcap_nume_sub = ((TH1D*)(h_pT_endcap_data_nume->Clone("h_pT_endcap_nume_sub")));
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_ttbar], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_tW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_tbarW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_WW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_WZ], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_ZZ], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_WJets], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_GJets_Full], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_MC_nume[_DY_Full], -1);
    // Far endcap
    TH1D *h_pT_endcap2_nume_sub = ((TH1D*)(h_pT_endcap2_data_nume->Clone("h_pT_endcap2_nume_sub")));
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_ttbar], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_tW], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_tbarW], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_WW], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_WZ], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_ZZ], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_WJets], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_GJets_Full], -1);
    h_pT_endcap2_nume_sub->Add(h_pT_endcap2_MC_nume[_DY_Full], -1);
    // Eta
    TH1D *h_eta_nume_sub = ((TH1D*)(h_eta_data_nume->Clone("h_eta_nume_sub")));
    h_eta_nume_sub->Add(h_eta_MC_nume[_ttbar], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_tW], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_tbarW], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_WW], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_WZ], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_ZZ], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_WJets], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_GJets_Full], -1);
    h_eta_nume_sub->Add(h_eta_MC_nume[_DY_Full], -1);

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_ctrl_sub = ((TH1D*)(h_pT_barrel_data_ctrl->Clone("h_pT_barrel_ctrl_sub")));
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_ttbar], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_tW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_tbarW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_WW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_WZ], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_ZZ], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_WJets], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_GJets_Full], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_MC_ctrl[_DY_Full], -1);
    TH1D *h_pT_barrel_deno_sub = ((TH1D*)(h_pT_barrel_ctrl_sub->Clone("h_pT_barrel_deno_sub")));
    h_pT_barrel_deno_sub->Add(h_pT_barrel_nume_sub); // deno = sig+ctrl
    // Endcap
    TH1D *h_pT_endcap_ctrl_sub = ((TH1D*)(h_pT_endcap_data_ctrl->Clone("h_pT_endcap_ctrl_sub")));
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_ttbar], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_tW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_tbarW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_WW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_WZ], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_ZZ], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_WJets], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_GJets_Full], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_MC_ctrl[_DY_Full], -1);
    TH1D *h_pT_endcap_deno_sub = ((TH1D*)(h_pT_endcap_ctrl_sub->Clone("h_pT_endcap_deno_sub")));
    h_pT_endcap_deno_sub->Add(h_pT_endcap_nume_sub); // deno = sig+ctrl
    // Far endcap
    TH1D *h_pT_endcap2_ctrl_sub = ((TH1D*)(h_pT_endcap2_data_ctrl->Clone("h_pT_endcap2_ctrl_sub")));
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_ttbar], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_tW], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_tbarW], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_WW], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_WZ], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_ZZ], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_WJets], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_GJets_Full], -1);
    h_pT_endcap2_ctrl_sub->Add(h_pT_endcap2_MC_ctrl[_DY_Full], -1);
    TH1D *h_pT_endcap2_deno_sub = ((TH1D*)(h_pT_endcap2_ctrl_sub->Clone("h_pT_endcap2_deno_sub")));
    h_pT_endcap2_deno_sub->Add(h_pT_endcap2_nume_sub); // deno = sig+ctrl
    // Eta
    TH1D *h_eta_ctrl_sub = ((TH1D*)(h_eta_data_ctrl->Clone("h_eta_ctrl_sub")));
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_ttbar], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_tW], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_tbarW], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_WW], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_WZ], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_ZZ], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_WJets], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_GJets_Full], -1);
    h_eta_ctrl_sub->Add(h_eta_MC_ctrl[_DY_Full], -1);
    TH1D *h_eta_deno_sub = ((TH1D*)(h_eta_ctrl_sub->Clone("h_eta_deno_sub")));
    h_eta_deno_sub->Add(h_eta_nume_sub); // deno = sig+ctrl

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRsubtract_barrel = ((TH1D*)(h_pT_barrel_nume_sub->Clone("h_FRsubtract_barrel")));
    h_FRsubtract_barrel->Divide(h_pT_barrel_deno_sub);
    h_FRsubtract_barrel->SetDirectory(0);
    cout << "Numerator barrel (subtraction): " << h_pT_barrel_nume_sub->Integral() << endl;
    cout << "Denominator barrel (subtraction): " << h_pT_barrel_deno_sub->Integral() << endl;
    // Endcap
    TH1D *h_FRsubtract_endcap = ((TH1D*)(h_pT_endcap_nume_sub->Clone("h_FRsubtract_endcap")));
    h_FRsubtract_endcap->Divide(h_pT_endcap_deno_sub);
    h_FRsubtract_endcap->SetDirectory(0);
    cout << "Numerator endcap (subtraction): " << h_pT_endcap_nume_sub->Integral() << endl;
    cout << "Denominator endcap (subtraction): " << h_pT_endcap_deno_sub->Integral() << endl;
    // Far endcap
    TH1D *h_FRsubtract_endcap2 = ((TH1D*)(h_pT_endcap2_nume_sub->Clone("h_FRsubtract_endcap2")));
    h_FRsubtract_endcap2->Divide(h_pT_endcap2_deno_sub);
    h_FRsubtract_endcap2->SetDirectory(0);
    cout << "Numerator far endcap (subtraction): " << h_pT_endcap2_nume_sub->Integral() << endl;
    cout << "Denominator far endcap (subtraction): " << h_pT_endcap2_deno_sub->Integral() << endl;
    // Eta
    TH1D *h_FRsubtract_eta = ((TH1D*)(h_eta_nume_sub->Clone("h_FRsubtract_eta")));
    h_FRsubtract_eta->Divide(h_eta_deno_sub);
    h_FRsubtract_eta->SetDirectory(0);

    // NOT PREPARED
//--------------------------------- FR by template -------------------------------------- (deno = nume + ctrl)
/*
    //                  QCD_nume(from ABCD)
    // FR = --------------------------------------------
    //      (QCD_nume + DATA_ctrl) - sum(nonQCD_MC_ctrl)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_abcd, *h_pT_endcap_nume_abcd;
    TFile *f_abcd = new TFile("/media/sf_DATA/FR/Electron/ABCD_hists.root", "READ");
    f_abcd->GetObject("h_pT_QCD_nume_barrel", h_pT_barrel_nume_abcd);
    f_abcd->GetObject("h_pT_QCD_nume_endcap", h_pT_endcap_nume_abcd);
    h_pT_barrel_nume_abcd->SetDirectory(0);
    h_pT_endcap_nume_abcd->SetDirectory(0);
    f_abcd->Close();

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_deno_abcd = ((TH1D*)(h_pT_barrel_data_ctrl->Clone("h_pT_barrel_deno_abcd")));
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_ttbar], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_tW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_tbarW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WZ], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_ZZ], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WJets], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_GJets_Full], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_DY_Full], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_nume_abcd); // deno = sig+ctrl
    // Endcap
    TH1D *h_pT_endcap_deno_abcd = ((TH1D*)(h_pT_endcap_data_ctrl->Clone("h_pT_endcap_deno_abcd")));
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_ttbar], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_tW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_tbarW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WZ], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_ZZ], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WJets], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_GJets_Full], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_DY_Full], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_nume_abcd); // deno = sig+ctrl
    // Far endcap
    TH1D *h_pT_endcap2_deno_abcd = ((TH1D*)(h_pT_endcap2_data_ctrl->Clone("h_pT_endcap2_deno_abcd")));
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_ttbar], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_tW], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_tbarW], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_WW], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_WZ], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_ZZ], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_WJets], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_GJets_Full], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_MC_ctrl[_DY_Full], -1);
    h_pT_endcap2_deno_abcd->Add(h_pT_endcap2_nume_abcd); // deno = sig+ctrl

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRtemplate_barrel = ((TH1D*)(h_pT_barrel_nume_abcd->Clone("h_FRtemplate_barrel")));
    h_FRtemplate_barrel->Divide(h_pT_barrel_deno_abcd);
    h_FRtemplate_barrel->SetDirectory(0);
    cout << "Numerator barrel (abcd): " << h_pT_barrel_nume_abcd->Integral() << endl;
    cout << "Denominator barrel (abcd): " << h_pT_barrel_deno_abcd->Integral() << endl;
    // Endcap
    TH1D *h_FRtemplate_endcap = ((TH1D*)(h_pT_endcap_nume_abcd->Clone("h_FRtemplate_endcap")));
    h_FRtemplate_endcap->Divide(h_pT_endcap_deno_abcd);
    h_FRtemplate_endcap->SetDirectory(0);
    cout << "Numerator endcap (abcd): " << h_pT_endcap_nume_abcd->Integral() << endl;
    cout << "Denominator endcap (abcd): " << h_pT_endcap_deno_abcd->Integral() << endl;
    // Far endcap
    TH1D *h_FRtemplate_endcap2 = ((TH1D*)(h_pT_endcap2_nume_abcd->Clone("h_FRtemplate_endcap2")));
    h_FRtemplate_endcap2->Divide(h_pT_endcap2_deno_abcd);
    h_FRtemplate_endcap2->SetDirectory(0);
    cout << "Numerator far endcap (abcd): " << h_pT_endcap2_nume_abcd->Integral() << endl;
    cout << "Denominator far endcap (abcd): " << h_pT_endcap2_deno_abcd->Integral() << endl;
    */

// --------------------- Geting the right errors -------------------- //
//                           ___________________________________
//       /  A  \      A*B    |  / Delta(A) \^2   / Delta(B) \^2 |
// Delta( ----- ) = -------  | ( ---------- ) + ( ---------- )    ;    Here A = Signal region, B = Control region
//       \ A+B /    (A+B)^2 \/  \    A     /     \    B     /          (numerator=signal, denominator=signal+control)

    for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
    {
        Double_t sig_barrel_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_barrel_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_barrel_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_barrel_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);

        Double_t sig_endcap_MC = h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_endcap_MC = h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_endcap_MC = h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_endcap_MC = h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);

        Double_t sig_endcap2_MC = h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_endcap2_MC = h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_endcap2_MC = h_pT_endcap2_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_endcap2_MC = h_pT_endcap2_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);


        Double_t sig_barrel_ratio = h_pT_barrel_nume->GetBinContent(i_bin);
        Double_t ctrl_barrel_ratio = h_pT_barrel_ctrl->GetBinContent(i_bin);
        Double_t err_sig_barrel_ratio = h_pT_barrel_nume->GetBinError(i_bin);
        Double_t err_ctrl_barrel_ratio = h_pT_barrel_ctrl->GetBinError(i_bin);

        Double_t sig_endcap_ratio = h_pT_endcap_nume->GetBinContent(i_bin);
        Double_t ctrl_endcap_ratio = h_pT_endcap_ctrl->GetBinContent(i_bin);
        Double_t err_sig_endcap_ratio = h_pT_endcap_nume->GetBinError(i_bin);
        Double_t err_ctrl_endcap_ratio = h_pT_endcap_ctrl->GetBinError(i_bin);

        Double_t sig_endcap2_ratio = h_pT_endcap2_nume->GetBinContent(i_bin);
        Double_t ctrl_endcap2_ratio = h_pT_endcap2_ctrl->GetBinContent(i_bin);
        Double_t err_sig_endcap2_ratio = h_pT_endcap2_nume->GetBinError(i_bin);
        Double_t err_ctrl_endcap2_ratio = h_pT_endcap2_ctrl->GetBinError(i_bin);

        Double_t sig_barrel_sub = h_pT_barrel_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_barrel_sub = h_pT_barrel_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_barrel_sub = h_pT_barrel_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_barrel_sub = h_pT_barrel_ctrl_sub->GetBinError(i_bin);

        Double_t sig_endcap_sub = h_pT_endcap_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_endcap_sub = h_pT_endcap_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_endcap_sub = h_pT_endcap_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_endcap_sub = h_pT_endcap_ctrl_sub->GetBinError(i_bin);

        Double_t sig_endcap2_sub = h_pT_endcap2_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_endcap2_sub = h_pT_endcap2_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_endcap2_sub = h_pT_endcap2_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_endcap2_sub = h_pT_endcap2_ctrl_sub->GetBinError(i_bin);

        Double_t err_FR_barrel_MC = sqrt((err_sig_barrel_MC / sig_barrel_MC) * (err_sig_barrel_MC / sig_barrel_MC) +
                                         (err_ctrl_barrel_MC / ctrl_barrel_MC) * (err_ctrl_barrel_MC / ctrl_barrel_MC));
        err_FR_barrel_MC *= sig_barrel_MC * ctrl_barrel_MC / ((sig_barrel_MC + ctrl_barrel_MC) * (sig_barrel_MC + ctrl_barrel_MC));

        Double_t err_FR_endcap_MC = sqrt((err_sig_endcap_MC / sig_endcap_MC) * (err_sig_endcap_MC / sig_endcap_MC) +
                                         (err_ctrl_endcap_MC / ctrl_endcap_MC) * (err_ctrl_endcap_MC / ctrl_endcap_MC));
        err_FR_endcap_MC *= sig_endcap_MC * ctrl_endcap_MC / ((sig_endcap_MC + ctrl_endcap_MC) * (sig_endcap_MC + ctrl_endcap_MC));

        Double_t err_FR_endcap2_MC = sqrt((err_sig_endcap2_MC / sig_endcap2_MC) * (err_sig_endcap2_MC / sig_endcap2_MC) +
                                          (err_ctrl_endcap2_MC / ctrl_endcap2_MC) * (err_ctrl_endcap2_MC / ctrl_endcap2_MC));
        err_FR_endcap2_MC *= sig_endcap2_MC * ctrl_endcap2_MC / ((sig_endcap2_MC + ctrl_endcap2_MC) * (sig_endcap2_MC + ctrl_endcap2_MC));

        Double_t err_FR_barrel_ratio = sqrt((err_sig_barrel_ratio / sig_barrel_ratio) * (err_sig_barrel_ratio / sig_barrel_ratio) +
                                            (err_ctrl_barrel_ratio / ctrl_barrel_ratio) * (err_ctrl_barrel_ratio / ctrl_barrel_ratio));
        err_FR_barrel_ratio *= sig_barrel_ratio * ctrl_barrel_ratio / ((sig_barrel_ratio + ctrl_barrel_ratio) * (sig_barrel_ratio + ctrl_barrel_ratio));

        Double_t err_FR_endcap_ratio = sqrt((err_sig_endcap_ratio / sig_endcap_ratio) * (err_sig_endcap_ratio / sig_endcap_ratio) +
                                            (err_ctrl_endcap_ratio / ctrl_endcap_ratio) * (err_ctrl_endcap_ratio / ctrl_endcap_ratio));
        err_FR_endcap_ratio *= sig_endcap_ratio * ctrl_endcap_ratio / ((sig_endcap_ratio + ctrl_endcap_ratio) * (sig_endcap_ratio + ctrl_endcap_ratio));

        Double_t err_FR_endcap2_ratio = sqrt((err_sig_endcap2_ratio / sig_endcap2_ratio) * (err_sig_endcap2_ratio / sig_endcap2_ratio) +
                                             (err_ctrl_endcap2_ratio / ctrl_endcap2_ratio) * (err_ctrl_endcap2_ratio / ctrl_endcap2_ratio));
        err_FR_endcap2_ratio *= sig_endcap2_ratio * ctrl_endcap2_ratio / ((sig_endcap2_ratio + ctrl_endcap2_ratio) * (sig_endcap2_ratio + ctrl_endcap2_ratio));

        Double_t err_FR_barrel_sub = sqrt((err_sig_barrel_sub / sig_barrel_sub) * (err_sig_barrel_sub / sig_barrel_sub) +
                                          (err_ctrl_barrel_sub / ctrl_barrel_sub) * (err_ctrl_barrel_sub / ctrl_barrel_sub));
        err_FR_barrel_sub *= sig_barrel_sub * ctrl_barrel_sub / ((sig_barrel_sub + ctrl_barrel_sub) * (sig_barrel_sub + ctrl_barrel_sub));

        Double_t err_FR_endcap_sub = sqrt((err_sig_endcap_sub / sig_endcap_sub) * (err_sig_endcap_sub / sig_endcap_sub) +
                                          (err_ctrl_endcap_sub / ctrl_endcap_sub) * (err_ctrl_endcap_sub / ctrl_endcap_sub));
        err_FR_endcap_sub *= sig_endcap_sub * ctrl_endcap_sub / ((sig_endcap_sub + ctrl_endcap_sub) * (sig_endcap_sub + ctrl_endcap_sub));

        Double_t err_FR_endcap2_sub = sqrt((err_sig_endcap2_sub / sig_endcap2_sub) * (err_sig_endcap2_sub / sig_endcap2_sub) +
                                           (err_ctrl_endcap2_sub / ctrl_endcap2_sub) * (err_ctrl_endcap2_sub / ctrl_endcap2_sub));
        err_FR_endcap2_sub *= sig_endcap2_sub * ctrl_endcap2_sub / ((sig_endcap2_sub + ctrl_endcap2_sub) * (sig_endcap2_sub + ctrl_endcap2_sub));

        h_FRMC_barrel->SetBinError(i_bin, err_FR_barrel_MC);
        h_FRMC_endcap->SetBinError(i_bin, err_FR_endcap_MC);
        h_FRMC_endcap2->SetBinError(i_bin, err_FR_endcap2_MC);
        h_FRratio_barrel->SetBinError(i_bin, err_FR_barrel_ratio);
        h_FRratio_endcap->SetBinError(i_bin, err_FR_endcap_ratio);
        h_FRratio_endcap2->SetBinError(i_bin, err_FR_endcap2_ratio);
        h_FRsubtract_barrel->SetBinError(i_bin, err_FR_barrel_sub);
        h_FRsubtract_endcap->SetBinError(i_bin, err_FR_endcap_sub);
        h_FRsubtract_endcap2->SetBinError(i_bin, err_FR_endcap2_sub);
    }

    for (Int_t i_bin=1; i_bin<=50; i_bin++)
    {
        Double_t sig_eta_MC = h_eta_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_eta_MC = h_eta_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_eta_MC = h_eta_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_eta_MC = h_eta_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);

        Double_t sig_eta_ratio = h_eta_nume->GetBinContent(i_bin);
        Double_t ctrl_eta_ratio = h_eta_ctrl->GetBinContent(i_bin);
        Double_t err_sig_eta_ratio = h_eta_nume->GetBinError(i_bin);
        Double_t err_ctrl_eta_ratio = h_eta_ctrl->GetBinError(i_bin);

        Double_t sig_eta_sub = h_eta_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_eta_sub = h_eta_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_eta_sub = h_eta_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_eta_sub = h_eta_ctrl_sub->GetBinError(i_bin);

        Double_t err_FR_eta_MC = sqrt((err_sig_eta_MC / sig_eta_MC) * (err_sig_eta_MC / sig_eta_MC) + (err_ctrl_eta_MC / ctrl_eta_MC) * (err_ctrl_eta_MC / ctrl_eta_MC));
        err_FR_eta_MC *= sig_eta_MC * ctrl_eta_MC / ((sig_eta_MC + ctrl_eta_MC) * (sig_eta_MC + ctrl_eta_MC));

        Double_t err_FR_eta_ratio = sqrt((err_sig_eta_ratio / sig_eta_ratio) * (err_sig_eta_ratio / sig_eta_ratio) + (err_ctrl_eta_ratio / ctrl_eta_ratio) * (err_ctrl_eta_ratio / ctrl_eta_ratio));
        err_FR_eta_ratio *= sig_eta_ratio * ctrl_eta_ratio / ((sig_eta_ratio + ctrl_eta_ratio) * (sig_eta_ratio + ctrl_eta_ratio));

        Double_t err_FR_eta_sub = sqrt((err_sig_eta_sub / sig_eta_sub) * (err_sig_eta_sub / sig_eta_sub) + (err_ctrl_eta_sub / ctrl_eta_sub) * (err_ctrl_eta_sub / ctrl_eta_sub));
        err_FR_eta_sub *= sig_eta_sub * ctrl_eta_sub / ((sig_eta_sub + ctrl_eta_sub) * (sig_eta_sub + ctrl_eta_sub));

        if (h_FRMC_eta->GetBinContent(i_bin) > 0)
            h_FRMC_eta->SetBinError(i_bin, err_FR_eta_MC);
        else h_FRMC_eta->SetBinError(i_bin, 0);
        if (h_FRratio_eta->GetBinContent(i_bin) > 0)
            h_FRratio_eta->SetBinError(i_bin, err_FR_eta_ratio);
        else h_FRratio_eta->SetBinError(i_bin, 0);
        if (h_FRsubtract_eta->GetBinContent(i_bin) > 0)
            h_FRsubtract_eta->SetBinError(i_bin, err_FR_eta_sub);
        else h_FRsubtract_eta->SetBinError(i_bin, 0);
    }


    // Drawing
    TCanvas *c_FR_barrel = new TCanvas("c_FR_barrel", "c_FR_barrel", 800, 800);
    c_FR_barrel->cd();
    c_FR_barrel->SetGrid(1);
    c_FR_barrel->SetLogx(1);
    c_FR_barrel->SetRightMargin(0.05);
    c_FR_barrel->SetTopMargin(0.05);
    c_FR_barrel->SetBottomMargin(0.12);
    c_FR_barrel->SetLeftMargin(0.13);
    h_FRMC_barrel->SetMarkerStyle(23);
    h_FRMC_barrel->SetMarkerColor(kYellow);
    h_FRMC_barrel->SetLineColor(kYellow);
    h_FRMC_barrel->SetStats(kFALSE);
    h_FRsubtract_barrel->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_barrel->SetMarkerColor(kBlack);
    h_FRsubtract_barrel->SetLineColor(kBlack);
    h_FRsubtract_barrel->SetStats(kFALSE);
    h_FRratio_barrel->SetMarkerStyle(kFullSquare);
    h_FRratio_barrel->SetMarkerColor(kRed);
    h_FRratio_barrel->SetLineColor(kRed);
    h_FRratio_barrel->SetStats(kFALSE);
    /*h_FRtemplate_barrel->SetMarkerStyle(33);
    h_FRtemplate_barrel->SetMarkerSize(1.5);
    h_FRtemplate_barrel->SetMarkerColor(kGreen+2);
    h_FRtemplate_barrel->SetLineColor(kGreen+2);
    h_FRtemplate_barrel->SetStats(kFALSE);*/
    h_FRsubtract_barrel->SetTitle("");
    h_FRsubtract_barrel->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_FRsubtract_barrel->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_barrel->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_barrel->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_barrel->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_barrel->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_barrel->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_barrel->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_barrel->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_barrel->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_barrel->GetXaxis()->SetRangeUser(25, 3000);
    h_FRsubtract_barrel->GetYaxis()->SetRangeUser(0, 0.3);
    h_FRsubtract_barrel->Draw();
//    h_FRratio_barrel->Draw("same");
//    h_FRMC_barrel->Draw("same");
//    h_FRtemplate_barrel->Draw("same");
    // Fit
    TF1 *f_barrel_25to38 = new TF1("f_barrel_25to38","gaus",25,50);
    TF1 *f_barrel_38to180 = new TF1("f_barrel_38to180","[0]+[1]*x+[2]*x^2",37,180);
    TF1 *f_barrel_180to350 = new TF1("f_barrel_180to350","[0]+[1]*x+[2]*x^2+[3]*x^3",130,470);
    TF1 *f_barrel_350to5000 = new TF1("f_barrel_350to5000","[0]+[1]*x+[2]*x^2+[3]*exp([4]*x)",300,5000);
    f_barrel_25to38->SetLineColor(kMagenta-5);
    f_barrel_38to180->SetLineColor(kMagenta-5);
    f_barrel_180to350->SetLineColor(kMagenta-5);
    f_barrel_350to5000->SetLineColor(kMagenta-5);
    f_barrel_350to5000->SetParameter(4, -0.1);
    h_FRsubtract_barrel->Fit("f_barrel_25to38", "R");
    h_FRsubtract_barrel->Fit("f_barrel_38to180", "R");
    h_FRsubtract_barrel->Fit("f_barrel_180to350", "R");
    h_FRsubtract_barrel->Fit("f_barrel_350to5000", "R");
//    f_barrel_25to38->Draw("same");
//    f_barrel_38to180->Draw("same");
//    f_barrel_180to350->Draw("same");
//    f_barrel_350to5000->Draw("same");
    Double_t chi2_barrel_50 = f_barrel_25to38->GetChisquare();
    Double_t chi2_barrel_180 = f_barrel_38to180->GetChisquare();
    Double_t chi2_barrel_500 = f_barrel_180to350->GetChisquare();
    Double_t chi2_barrel_5000 = f_barrel_350to5000->GetChisquare();
    cout << "Barrel fit chi squares: " << chi2_barrel_50 << "  " << chi2_barrel_180 << "  " << chi2_barrel_500 << "  " << chi2_barrel_5000 << endl;
    cout << "\n\n";

    TLegend *legend = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend->AddEntry(h_FRMC_barrel, "QCD MC", "LP");
//    legend->AddEntry(h_FRratio_barrel, "Ratio", "LP");
    legend->AddEntry(h_FRsubtract_barrel, "Subtraction", "LP");
    TLegend *legend_noABCD = ((TLegend*)legend->Clone());
//    legend->AddEntry(h_FRtemplate_barrel, "ABCD", "LP");
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
    h_FRMC_endcap->SetMarkerStyle(23);
    h_FRMC_endcap->SetMarkerColor(kYellow);
    h_FRMC_endcap->SetLineColor(kYellow);
    h_FRMC_endcap->SetStats(kFALSE);
    h_FRsubtract_endcap->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_endcap->SetMarkerColor(kBlack);
    h_FRsubtract_endcap->SetLineColor(kBlack);
    h_FRsubtract_endcap->SetStats(kFALSE);
    h_FRratio_endcap->SetMarkerStyle(kFullSquare);
    h_FRratio_endcap->SetMarkerColor(kRed);
    h_FRratio_endcap->SetLineColor(kRed);
    h_FRratio_endcap->SetStats(kFALSE);
    /*h_FRtemplate_endcap->SetMarkerStyle(33);
    h_FRtemplate_endcap->SetMarkerSize(1.5);
    h_FRtemplate_endcap->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap->SetLineColor(kGreen+2);
    h_FRtemplate_endcap->SetStats(kFALSE);*/
    h_FRsubtract_endcap->SetTitle("");
    h_FRsubtract_endcap->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_FRsubtract_endcap->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_endcap->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_endcap->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_endcap->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_endcap->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_endcap->GetXaxis()->SetRangeUser(25, 5000);
    h_FRsubtract_endcap->GetYaxis()->SetRangeUser(0, 0.3);
    h_FRsubtract_endcap->Draw();
    // Fit
    TF1 *f_endcap_25to108 = new TF1("f_endcap_25to108","[0]+[1]*x+[2]*x^2",25,110);
    f_endcap_25to108->SetParameter(0,1);
    f_endcap_25to108->SetParameter(1,-0.01);
    f_endcap_25to108->SetParameter(2,0.01);
    f_endcap_25to108->SetLineColor(kAzure+1);
    TF1 *f_endcap_108to350 = new TF1("f_endcap_108to350","[0]+[1]*x+[2]*x^2",95,400);
    f_endcap_108to350->SetParameter(0,0);
    f_endcap_108to350->SetParameter(1,0.01);
    f_endcap_108to350->SetParameter(2,0.01);
    f_endcap_108to350->SetLineColor(kAzure+1);
    TF1 *f_endcap_350to5000 = new TF1("f_endcap_350to5000","[0]+[1]*exp([2]*x)",300,5000);
    f_endcap_350to5000->SetParameter(0,0.2);
    f_endcap_350to5000->SetParameter(1,-0.01);
    f_endcap_350to5000->SetParameter(1,0.01);
    f_endcap_350to5000->SetLineColor(kAzure+1);
    h_FRsubtract_endcap->Fit("f_endcap_25to108", "R");
    h_FRsubtract_endcap->Fit("f_endcap_108to350", "R");
    h_FRsubtract_endcap->Fit("f_endcap_350to5000", "R");
//    f_endcap_25to108->Draw("same");
//    f_endcap_108to350->Draw("same");
//    f_endcap_350to5000->Draw("same");
    Double_t chi2_endcap_100 = f_endcap_25to108->GetChisquare();
    Double_t chi2_endcap_400 = f_endcap_108to350->GetChisquare();
    Double_t chi2_endcap_1000 = f_endcap_350to5000->GetChisquare();
    cout << "Endcap fit chi squares: " << chi2_endcap_100 << "  " << chi2_endcap_400 << "  " << chi2_endcap_1000 << endl;
    cout << "\n\n";


//    h_FRratio_endcap->Draw("same");
//    h_FRMC_endcap->Draw("same");
//    h_FRtemplate_endcap->Draw("same");
    legend->Draw();
    TLatex *texte = new TLatex (0.45, 0.6, "Endcap #eta < 2.2");
    texte->SetTextAlign(11);
    texte->SetTextSize(0.05);
    texte->SetNDC(true);
    texte->Draw();
    c_FR_endcap->Update();

    TCanvas *c_FR_endcap2 = new TCanvas("c_FR_endcap2", "c_FR_endcap2", 800, 800);
    c_FR_endcap2->cd();
    c_FR_endcap2->cd();
    c_FR_endcap2->SetGrid(1);
    c_FR_endcap2->SetLogx(1);
    c_FR_endcap2->SetRightMargin(0.05);
    c_FR_endcap2->SetTopMargin(0.05);
    c_FR_endcap2->SetBottomMargin(0.12);
    c_FR_endcap2->SetLeftMargin(0.13);
    h_FRMC_endcap2->SetMarkerStyle(23);
    h_FRMC_endcap2->SetMarkerColor(kYellow);
    h_FRMC_endcap2->SetLineColor(kYellow);
    h_FRMC_endcap2->SetStats(kFALSE);
    h_FRsubtract_endcap2->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_endcap2->SetMarkerColor(kBlack);
    h_FRsubtract_endcap2->SetLineColor(kBlack);
    h_FRsubtract_endcap2->SetStats(kFALSE);
    h_FRratio_endcap2->SetMarkerStyle(kFullSquare);
    h_FRratio_endcap2->SetMarkerColor(kRed);
    h_FRratio_endcap2->SetLineColor(kRed);
    h_FRratio_endcap2->SetStats(kFALSE);
    /*h_FRtemplate_endcap2->SetMarkerStyle(33);
    h_FRtemplate_endcap2->SetMarkerSize(1.5);
    h_FRtemplate_endcap2->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap2->SetLineColor(kGreen+2);
    h_FRtemplate_endcap2->SetStats(kFALSE);*/
    h_FRsubtract_endcap2->SetTitle("");
    h_FRsubtract_endcap2->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_FRsubtract_endcap2->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_endcap2->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap2->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap2->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_endcap2->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap2->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_endcap2->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap2->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_endcap2->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_endcap2->GetXaxis()->SetRangeUser(25, 3000);
    h_FRsubtract_endcap2->GetYaxis()->SetRangeUser(0, 0.3);
    h_FRsubtract_endcap2->Draw();
//    h_FRratio_endcap2->Draw("same");
//    h_FRMC_endcap2->Draw("same");
//    h_FRtemplate_endcap2->Draw("same");
    TF1 *f_endcap2_25to200 = new TF1("f_endcap2_25to200","[0]+[1]*x+[2]*x^2+[3]*x^3",25,200);
    f_endcap2_25to200->SetLineColor(kOrange+5);
    TF1 *f_endcap2_200to350 = new TF1("f_endcap2_200to350","gaus+[3]",160, 400);
    f_endcap2_200to350->SetParameter(0, 0.11);
    f_endcap2_200to350->SetParLimits(0, 0.11, 0.14);
    f_endcap2_200to350->SetParameter(1, 350);
    f_endcap2_200to350->SetParLimits(1, 330, 360);
    f_endcap2_200to350->SetParameter(2, 140);
    f_endcap2_200to350->SetParLimits(2, 130,150);
    f_endcap2_200to350->SetParameter(3, 0.04);
    f_endcap2_200to350->SetParLimits(3, 0.025, 0.1);
    f_endcap2_200to350->SetLineColor(kOrange+5);
    TF1 *f_endcap2_350to5000 = new TF1("f_endcap2_350to5000","landau+[3]",300,5000);
    f_endcap2_350to5000->SetParameter(0, 0.1);
    f_endcap2_350to5000->SetParameter(1, 350);
    f_endcap2_350to5000->SetParLimits(2, 100, 200);
    f_endcap2_350to5000->SetParLimits(3, 0.11, 0.12);
    f_endcap2_350to5000->SetLineColor(kOrange+5);
    h_FRsubtract_endcap2->Fit("f_endcap2_25to200", "R");
    h_FRsubtract_endcap2->Fit("f_endcap2_200to350", "R");
    h_FRsubtract_endcap2->Fit("f_endcap2_350to5000", "R");
//    f_endcap2_25to200->Draw("same");
//    f_endcap2_200to350->Draw("same");
//    f_endcap2_350to5000->Draw("same");
    Double_t chi2_endcap2_200 = f_endcap2_25to200->GetChisquare();
    Double_t chi2_endcap2_400 = f_endcap2_200to350->GetChisquare();
    Double_t chi2_endcap2_5000 = f_endcap2_350to5000->GetChisquare();
    cout << "Far endcap fit chi squares: " << chi2_endcap2_200 << "  " << chi2_endcap2_400 << "  " << chi2_endcap2_5000 << endl;

    legend->Draw();
    TLatex *texte2 = new TLatex (0.45, 0.6, "Endcap #eta #geq 2.2");
    texte2->SetTextAlign(11);
    texte2->SetTextSize(0.05);
    texte2->SetNDC(true);
    texte2->Draw();
    c_FR_endcap2->Update();

    TCanvas *c_FR_allin1 = new TCanvas("c_FR_allin1", "c_FR_allin1", 800, 800);
    c_FR_allin1->cd();
    c_FR_allin1->SetGrid(1);
    c_FR_allin1->SetRightMargin(0.05);
    c_FR_allin1->SetTopMargin(0.05);
    c_FR_allin1->SetBottomMargin(0.12);
    c_FR_allin1->SetLeftMargin(0.13);
    TH1D *h_FRsubtract_endcap_plot = ((TH1D*)(h_FRsubtract_endcap->Clone("h_FRsubtract_endcap_plot")));
    h_FRsubtract_endcap_plot->SetMarkerStyle(kFullSquare);
    h_FRsubtract_endcap_plot->SetMarkerColor(kBlue);
    h_FRsubtract_endcap_plot->SetLineColor(kBlue);
    h_FRsubtract_endcap_plot->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_endcap_plot->GetYaxis()->SetTitleSize(0.045);
    h_FRsubtract_endcap_plot->GetYaxis()->SetTitleOffset(1.4);
    h_FRsubtract_endcap_plot->Draw();
    TH1D *h_FRsubtract_endcap2_plot = ((TH1D*)(h_FRsubtract_endcap2->Clone("h_FRsubtract_endcap2_plot")));
    h_FRsubtract_endcap2_plot->SetMarkerStyle(33);
    h_FRsubtract_endcap2_plot->SetMarkerSize(1.5);
    h_FRsubtract_endcap2_plot->SetMarkerColor(kOrange-3);
    h_FRsubtract_endcap2_plot->SetLineColor(kOrange-3);
    h_FRsubtract_endcap2_plot->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_endcap2_plot->GetYaxis()->SetTitleSize(0.045);
    h_FRsubtract_endcap2_plot->GetYaxis()->SetTitleOffset(1.4);
    h_FRsubtract_endcap2_plot->Draw("same");
    h_FRsubtract_barrel->Draw("same");
    TLegend *legend1 = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend1->AddEntry(h_FRsubtract_barrel, "Subtraction (barrel)", "LP");
//    legend1->AddEntry(h_FRsubtract_endcap_plot, "Subtraction (endcap |#eta| < 2.2)", "LP");
//    legend1->AddEntry(h_FRsubtract_endcap2_plot, "Subtraction (endcap |#eta| #geq 2.2)", "LP");
    legend1->AddEntry(h_FRsubtract_barrel, "|#eta| < 1.4442", "LP");
    legend1->AddEntry(h_FRsubtract_endcap_plot, "1.566 < |#eta| < 2.2", "LP");
    legend1->AddEntry(h_FRsubtract_endcap2_plot, "|#eta| #geq 2.2", "LP");
    legend1->Draw();
//    f_barrel_25to38->Draw("same");
//    f_barrel_38to180->Draw("same");
//    f_barrel_180to350->Draw("same");
//    f_barrel_350to5000->Draw("same");
//    f_endcap_25to108->Draw("same");
//    f_endcap_108to350->Draw("same");
//    f_endcap_350to5000->Draw("same");
//    f_endcap2_25to200->Draw("same");
//    f_endcap2_200to350->Draw("same");
//    f_endcap2_350to5000->Draw("same");
    c_FR_allin1->SetLogx();
    c_FR_allin1->Update();

    TCanvas *c_FR_eta = new TCanvas("c_FR_eta", "c_FR_eta", 800, 800);
    c_FR_eta->cd();
    c_FR_eta->SetGrid(1);
    c_FR_eta->SetRightMargin(0.05);
    c_FR_eta->SetTopMargin(0.05);
    c_FR_eta->SetBottomMargin(0.12);
    c_FR_eta->SetLeftMargin(0.13);
    h_FRMC_eta->SetMarkerStyle(23);
    h_FRMC_eta->SetMarkerColor(kOrange-3);
    h_FRMC_eta->SetLineColor(kOrange-3);
    h_FRMC_eta->SetStats(kFALSE);
    h_FRsubtract_eta->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_eta->SetMarkerColor(kBlack);
    h_FRsubtract_eta->SetLineColor(kBlack);
    h_FRsubtract_eta->SetStats(kFALSE);
    h_FRratio_eta->SetMarkerStyle(kFullSquare);
    h_FRratio_eta->SetMarkerColor(kRed);
    h_FRratio_eta->SetLineColor(kRed);
    h_FRratio_eta->SetStats(kFALSE);
    h_FRsubtract_eta->SetTitle("");
    h_FRsubtract_eta->GetXaxis()->SetTitle("#eta (#mu)");
    h_FRsubtract_eta->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_eta->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_eta->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_eta->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_eta->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_eta->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_eta->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_eta->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_eta->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_eta->GetXaxis()->SetRangeUser(25, 1000);
    h_FRsubtract_eta->GetYaxis()->SetRangeUser(0, 1);
    h_FRsubtract_eta->Draw();
//    h_FRratio_eta->Draw("same");
//    h_FRMC_eta->Draw("same");
    legend_noABCD->Draw();
    c_FR_eta->Update();


// ---------------------------- Writing------------------------------ //
    TString filename = "/media/sf_DATA/FR/Electron/FakeRate_electron.root";
    TFile *file_FR = new TFile(filename, "RECREATE");
    if (file_FR->IsOpen()) cout << "File '" << filename << "' has been created. Writing histograms.." << endl;
    file_FR->cd();
    h_FRMC_barrel->Write();
    h_FRMC_endcap->Write();
    h_FRMC_endcap2->Write();
    h_FRratio_barrel->Write();
    h_FRratio_endcap->Write();
    h_FRratio_endcap2->Write();
    h_FRsubtract_barrel->Write();
    h_FRsubtract_endcap->Write();
    h_FRsubtract_endcap2->Write();
//    h_FRtemplate_barrel->Write();
//    h_FRtemplate_endcap->Write();
//    h_FRtemplate_endcap2->Write();
    f_barrel_25to38->Write();
    f_barrel_38to180->Write();
    f_barrel_180to350->Write();
    f_barrel_350to5000->Write();
    f_endcap_25to108->Write();
    f_endcap_108to350->Write();
    f_endcap_350to5000->Write();
    f_endcap2_25to200->Write();
    f_endcap2_200to350->Write();
    f_endcap2_350to5000->Write();
    cout << "Finished. Closing the file.." << endl;
    file_FR->Close();
    if (!file_FR->IsOpen()) cout << "File '" << filename << "' has been closed successfully." << endl;
    else cout << "File did not close!" << endl;

} // End of E_EstFR()


void E_EstFR_alt()
{
    TFile *f = new TFile("/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root", "READ");

    TH1D *h_pT_barrel_lead_pass[4],
         *h_pT_endcap_lead_pass[4],
         *h_pT_barrel_lead_fail[4],
         *h_pT_endcap_lead_fail[4],
         *h_pT_barrel_sub_pass[4],
         *h_pT_endcap_sub_pass[4],
         *h_pT_barrel_sub_fail[4],
         *h_pT_endcap_sub_fail[4];
    TString type[4] = {"data", "DY", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_lead_pass_"+type[i], h_pT_barrel_lead_pass[i]);
        f->GetObject("h_pT_endcap_lead_pass_"+type[i], h_pT_endcap_lead_pass[i]);
        f->GetObject("h_pT_barrel_lead_fail_"+type[i], h_pT_barrel_lead_fail[i]);
        f->GetObject("h_pT_endcap_lead_fail_"+type[i], h_pT_endcap_lead_fail[i]);
        f->GetObject("h_pT_barrel_sub_pass_"+type[i], h_pT_barrel_sub_pass[i]);
        f->GetObject("h_pT_endcap_sub_pass_"+type[i], h_pT_endcap_sub_pass[i]);
        f->GetObject("h_pT_barrel_sub_fail_"+type[i], h_pT_barrel_sub_fail[i]);
        f->GetObject("h_pT_endcap_sub_fail_"+type[i], h_pT_endcap_sub_fail[i]);
        h_pT_barrel_lead_pass[i]->SetDirectory(0);
        h_pT_endcap_lead_pass[i]->SetDirectory(0);
        h_pT_barrel_lead_fail[i]->SetDirectory(0);
        h_pT_endcap_lead_fail[i]->SetDirectory(0);
        h_pT_barrel_sub_pass[i]->SetDirectory(0);
        h_pT_endcap_sub_pass[i]->SetDirectory(0);
        h_pT_barrel_sub_fail[i]->SetDirectory(0);
        h_pT_endcap_sub_fail[i]->SetDirectory(0);
    }

    TH1D *h_FR_subtract_barrel_lead = ((TH1D*)(h_pT_barrel_lead_pass[0]->Clone("h_FR_subtract_lead_barrel")));
    TH1D *h_FR_subtract_endcap_lead = ((TH1D*)(h_pT_endcap_lead_pass[0]->Clone("h_FR_subtract_lead_endcap")));
    TH1D *h_FR_MC_barrel_lead = ((TH1D*)(h_pT_barrel_lead_pass[3]->Clone("h_FR_MC_lead_barrel")));
    TH1D *h_FR_MC_endcap_lead = ((TH1D*)(h_pT_endcap_lead_pass[3]->Clone("h_FR_MC_lead_endcap")));
    TH1D *h_FR_subtract_barrel_sub = ((TH1D*)(h_pT_barrel_sub_pass[0]->Clone("h_FR_subtract_sub_barrel")));
    TH1D *h_FR_subtract_endcap_sub = ((TH1D*)(h_pT_endcap_sub_pass[0]->Clone("h_FR_subtract_sub_endcap")));
    TH1D *h_FR_MC_barrel_sub = ((TH1D*)(h_pT_barrel_sub_pass[3]->Clone("h_FR_MC_sub_barrel")));
    TH1D *h_FR_MC_endcap_sub = ((TH1D*)(h_pT_endcap_sub_pass[3]->Clone("h_FR_MC_sub_endcap")));
    h_FR_subtract_barrel_lead->Add(h_pT_barrel_lead_pass[1], -1);
    h_FR_subtract_endcap_lead->Add(h_pT_endcap_lead_pass[1], -1);
    h_FR_subtract_barrel_lead->Add(h_pT_barrel_lead_pass[2], -1);
    h_FR_subtract_endcap_lead->Add(h_pT_endcap_lead_pass[2], -1);
    h_FR_subtract_barrel_sub->Add(h_pT_barrel_sub_pass[1], -1);
    h_FR_subtract_endcap_sub->Add(h_pT_endcap_sub_pass[1], -1);
    h_FR_subtract_barrel_sub->Add(h_pT_barrel_sub_pass[2], -1);
    h_FR_subtract_endcap_sub->Add(h_pT_endcap_sub_pass[2], -1);


    TH1D *h_FR_subtract_barrel_lead_deno = ((TH1D*)(h_pT_barrel_lead_fail[0]->Clone("h_FR_subtract_lead_barrel_deno")));
    TH1D *h_FR_subtract_endcap_lead_deno = ((TH1D*)(h_pT_barrel_lead_fail[0]->Clone("h_FR_subtract_lead_endcap_deno")));
    TH1D *h_FR_MC_barrel_lead_deno = ((TH1D*)(h_pT_barrel_lead_fail[3]->Clone("h_FR_MC_lead_barrel_deno")));
    TH1D *h_FR_MC_endcap_lead_deno = ((TH1D*)(h_pT_barrel_lead_fail[3]->Clone("h_FR_MC_lead_endcap_deno")));
    TH1D *h_FR_subtract_barrel_sub_deno = ((TH1D*)(h_pT_barrel_sub_fail[0]->Clone("h_FR_subtract_sub_barrel_deno")));
    TH1D *h_FR_subtract_endcap_sub_deno = ((TH1D*)(h_pT_barrel_sub_fail[0]->Clone("h_FR_subtract_sub_endcap_deno")));
    TH1D *h_FR_MC_barrel_sub_deno = ((TH1D*)(h_pT_barrel_sub_fail[3]->Clone("h_FR_MC_sub_barrel_deno")));
    TH1D *h_FR_MC_endcap_sub_deno = ((TH1D*)(h_pT_barrel_sub_fail[3]->Clone("h_FR_MC_sub_endcap_deno")));
    h_FR_subtract_barrel_lead_deno->Add(h_pT_barrel_lead_pass[0]);
    h_FR_subtract_endcap_lead_deno->Add(h_pT_endcap_lead_pass[0]);
    h_FR_MC_barrel_lead_deno->Add(h_pT_barrel_lead_pass[3]);
    h_FR_MC_endcap_lead_deno->Add(h_pT_endcap_lead_pass[3]);
    h_FR_subtract_barrel_lead_deno->Add(h_pT_barrel_lead_fail[1], -1);
    h_FR_subtract_endcap_lead_deno->Add(h_pT_endcap_lead_fail[1], -1);
    h_FR_subtract_barrel_lead_deno->Add(h_pT_barrel_lead_fail[2], -1);
    h_FR_subtract_endcap_lead_deno->Add(h_pT_endcap_lead_fail[2], -1);
    h_FR_subtract_barrel_lead_deno->Add(h_pT_barrel_lead_pass[1], -1);
    h_FR_subtract_endcap_lead_deno->Add(h_pT_endcap_lead_pass[1], -1);
    h_FR_subtract_barrel_lead_deno->Add(h_pT_barrel_lead_pass[2], -1);
    h_FR_subtract_endcap_lead_deno->Add(h_pT_endcap_lead_pass[2], -1);
    h_FR_subtract_barrel_sub_deno->Add(h_pT_barrel_sub_pass[0]);
    h_FR_subtract_endcap_sub_deno->Add(h_pT_endcap_sub_pass[0]);
    h_FR_MC_barrel_sub_deno->Add(h_pT_barrel_sub_pass[3]);
    h_FR_MC_endcap_sub_deno->Add(h_pT_endcap_sub_pass[3]);
    h_FR_subtract_barrel_sub_deno->Add(h_pT_barrel_sub_fail[1], -1);
    h_FR_subtract_endcap_sub_deno->Add(h_pT_endcap_sub_fail[1], -1);
    h_FR_subtract_barrel_sub_deno->Add(h_pT_barrel_sub_fail[2], -1);
    h_FR_subtract_endcap_sub_deno->Add(h_pT_endcap_sub_fail[2], -1);
    h_FR_subtract_barrel_sub_deno->Add(h_pT_barrel_sub_pass[1], -1);
    h_FR_subtract_endcap_sub_deno->Add(h_pT_endcap_sub_pass[1], -1);
    h_FR_subtract_barrel_sub_deno->Add(h_pT_barrel_sub_pass[2], -1);
    h_FR_subtract_endcap_sub_deno->Add(h_pT_endcap_sub_pass[2], -1);

    h_FR_subtract_barrel_lead->Divide(h_FR_subtract_barrel_lead_deno);
    h_FR_subtract_endcap_lead->Divide(h_FR_subtract_endcap_lead_deno);
    h_FR_MC_barrel_lead->Divide(h_FR_MC_barrel_lead_deno);
    h_FR_MC_endcap_lead->Divide(h_FR_MC_endcap_lead_deno);
    h_FR_subtract_barrel_sub->Divide(h_FR_subtract_barrel_sub_deno);
    h_FR_subtract_endcap_sub->Divide(h_FR_subtract_endcap_sub_deno);
    h_FR_MC_barrel_sub->Divide(h_FR_MC_barrel_sub_deno);
    h_FR_MC_endcap_sub->Divide(h_FR_MC_endcap_sub_deno);

    h_FR_subtract_barrel_lead->SetDirectory(0);
    h_FR_subtract_barrel_lead->SetStats(0);
    h_FR_subtract_barrel_lead->SetMarkerStyle(20);
    h_FR_subtract_barrel_lead->SetMarkerColor(kBlack);
    h_FR_subtract_barrel_lead->SetLineColor(kBlack);
    h_FR_subtract_endcap_lead->SetDirectory(0);
    h_FR_subtract_endcap_lead->SetStats(0);
    h_FR_subtract_endcap_lead->SetMarkerStyle(20);
    h_FR_subtract_endcap_lead->SetMarkerColor(kBlack);
    h_FR_subtract_endcap_lead->SetLineColor(kBlack);
    h_FR_MC_barrel_lead->SetDirectory(0);
    h_FR_MC_barrel_lead->SetStats(0);
    h_FR_MC_barrel_lead->SetMarkerStyle(20);
    h_FR_MC_barrel_lead->SetMarkerColor(kRed);
    h_FR_MC_barrel_lead->SetLineColor(kRed);
    h_FR_MC_endcap_lead->SetDirectory(0);
    h_FR_MC_endcap_lead->SetStats(0);
    h_FR_MC_endcap_lead->SetMarkerStyle(20);
    h_FR_MC_endcap_lead->SetMarkerColor(kRed);
    h_FR_MC_endcap_lead->SetLineColor(kRed);
    h_FR_subtract_barrel_sub->SetDirectory(0);
    h_FR_subtract_barrel_sub->SetStats(0);
    h_FR_subtract_barrel_sub->SetMarkerStyle(21);
    h_FR_subtract_barrel_sub->SetMarkerColor(kBlack);
    h_FR_subtract_barrel_sub->SetLineColor(kBlack);
    h_FR_subtract_endcap_sub->SetDirectory(0);
    h_FR_subtract_endcap_sub->SetStats(0);
    h_FR_subtract_endcap_sub->SetMarkerStyle(21);
    h_FR_subtract_endcap_sub->SetMarkerColor(kBlack);
    h_FR_subtract_endcap_sub->SetLineColor(kBlack);
    h_FR_MC_barrel_sub->SetDirectory(0);
    h_FR_MC_barrel_sub->SetStats(0);
    h_FR_MC_barrel_sub->SetMarkerStyle(21);
    h_FR_MC_barrel_sub->SetMarkerColor(kRed);
    h_FR_MC_barrel_sub->SetLineColor(kRed);
    h_FR_MC_endcap_sub->SetDirectory(0);
    h_FR_MC_endcap_sub->SetStats(0);
    h_FR_MC_endcap_sub->SetMarkerStyle(21);
    h_FR_MC_endcap_sub->SetMarkerColor(kRed);
    h_FR_MC_endcap_sub->SetLineColor(kRed);


    // Creating and drawing FR plots
    TLegend * legend = new TLegend(0.6, 0.8, 0.95, 0.95);
    legend->AddEntry(h_FR_subtract_barrel_lead, "Subtraction (lead)", "pl");
    legend->AddEntry(h_FR_subtract_barrel_sub, "Subtraction (sub)", "pl");
    legend->AddEntry(h_FR_MC_barrel_lead, "MC (lead)", "pl");
    legend->AddEntry(h_FR_MC_barrel_sub, "MC (sub)", "pl");

    TCanvas *c_FR_barrel = new TCanvas("c_FR_barrel", "c_FR_barrel", 800, 800);
    c_FR_barrel->SetTopMargin(0.05);
    c_FR_barrel->SetRightMargin(0.05);
    c_FR_barrel->SetBottomMargin(0.15);
    c_FR_barrel->SetLeftMargin(0.15);
    h_FR_subtract_barrel_lead->SetTitle("");
    h_FR_subtract_barrel_lead->GetXaxis()->SetTitle("p_{#lower[-0.2]{T}} [GeV/c]");
    h_FR_subtract_barrel_lead->GetXaxis()->SetTitleSize(0.062);
    h_FR_subtract_barrel_lead->GetXaxis()->SetTitleOffset(0.9);
    h_FR_subtract_barrel_lead->GetXaxis()->SetLabelSize(0.048);
    h_FR_subtract_barrel_lead->GetXaxis()->SetMoreLogLabels();
    h_FR_subtract_barrel_lead->GetXaxis()->SetNoExponent();
    h_FR_subtract_barrel_lead->GetYaxis()->SetTitle("Fake Rate");
    h_FR_subtract_barrel_lead->GetYaxis()->SetTitleSize(0.05);
    h_FR_subtract_barrel_lead->GetYaxis()->SetTitleOffset(1.25);
    h_FR_subtract_barrel_lead->GetYaxis()->SetLabelSize(0.043);
    h_FR_subtract_barrel_lead->GetYaxis()->SetMoreLogLabels();
    h_FR_subtract_barrel_lead->GetYaxis()->SetNoExponent();
    h_FR_subtract_barrel_lead->Draw();
    h_FR_subtract_barrel_lead->GetYaxis()->SetRangeUser(0, 1);
    h_FR_subtract_barrel_sub->Draw("same");
    h_FR_MC_barrel_lead->Draw("same");
    h_FR_MC_barrel_sub->Draw("same");
    legend->Draw();
    c_FR_barrel->SetLogx();
    c_FR_barrel->SetGridx();
    c_FR_barrel->SetGridy();
    c_FR_barrel->Update();

    TCanvas *c_FR_endcap = new TCanvas("c_FR_endcap", "c_FR_endcap", 800, 800);
    c_FR_endcap->SetTopMargin(0.05);
    c_FR_endcap->SetRightMargin(0.05);
    c_FR_endcap->SetBottomMargin(0.15);
    c_FR_endcap->SetLeftMargin(0.15);
    h_FR_subtract_endcap_lead->SetTitle("");
    h_FR_subtract_endcap_lead->GetXaxis()->SetTitle("p_{#lower[-0.2]{T}} [GeV/c]");
    h_FR_subtract_endcap_lead->GetXaxis()->SetTitleSize(0.062);
    h_FR_subtract_endcap_lead->GetXaxis()->SetTitleOffset(0.9);
    h_FR_subtract_endcap_lead->GetXaxis()->SetLabelSize(0.048);
    h_FR_subtract_endcap_lead->GetXaxis()->SetMoreLogLabels();
    h_FR_subtract_endcap_lead->GetXaxis()->SetNoExponent();
    h_FR_subtract_endcap_lead->GetYaxis()->SetTitle("Fake Rate");
    h_FR_subtract_endcap_lead->GetYaxis()->SetTitleSize(0.05);
    h_FR_subtract_endcap_lead->GetYaxis()->SetTitleOffset(1.25);
    h_FR_subtract_endcap_lead->GetYaxis()->SetLabelSize(0.043);
    h_FR_subtract_endcap_lead->GetYaxis()->SetMoreLogLabels();
    h_FR_subtract_endcap_lead->GetYaxis()->SetNoExponent();
    h_FR_subtract_endcap_lead->Draw();
    h_FR_subtract_endcap_lead->GetYaxis()->SetRangeUser(0, 1);
    h_FR_subtract_endcap_sub->Draw("same");
    h_FR_MC_endcap_lead->Draw("same");
    h_FR_MC_endcap_sub->Draw("same");
    legend->Draw();
    c_FR_endcap->SetLogx();
    c_FR_endcap->SetGridx();
    c_FR_endcap->SetGridy();
    c_FR_endcap->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << "/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " <<"/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root" << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/FR/Electron/FakeRate_electron_alt.root", "RECREATE");
    f_out->cd();
    h_FR_subtract_barrel_lead->Write();
    h_FR_subtract_endcap_lead->Write();
    h_FR_MC_barrel_lead->Write();
    h_FR_MC_endcap_lead->Write();
    h_FR_subtract_barrel_sub->Write();
    h_FR_subtract_endcap_sub->Write();
    h_FR_MC_barrel_sub->Write();
    h_FR_MC_endcap_sub->Write();
    f_out->Close();

} // End of E_EstFR_alt


void E_EstFR_alt2()
{
    FileMgr fm;

    TH1D *h_pT_barrel_nume[_EndOf_Data_Special], *h_pT_endcap_nume[_EndOf_Data_Special],
         *h_pT_barrel_ctrl[_EndOf_Data_Special], *h_pT_endcap_ctrl[_EndOf_Data_Special],
         *h_eta_nume[_EndOf_Data_Special], *h_eta_ctrl[_EndOf_Data_Special];

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr = _WW;
    while (!stop)
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);

        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);

        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_eta_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]->SetDirectory(0);

        if (pr == _ttbar)
        {
            h_pT_barrel_ctrl[_ttbar_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_ttbar")));
            h_pT_endcap_ctrl[_ttbar_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_ttbar")));
            h_pT_barrel_nume[_ttbar_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_ttbar")));
            h_pT_endcap_nume[_ttbar_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_ttbar")));
            h_eta_ctrl[_ttbar_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_ttbar")));
            h_eta_nume[_ttbar_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_ttbar")));

            h_pT_barrel_ctrl[_ttbar_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_ttbar_Full]->SetDirectory(0);
            h_pT_barrel_nume[_ttbar_Full]->SetDirectory(0);
            h_pT_endcap_nume[_ttbar_Full]->SetDirectory(0);
            h_eta_ctrl[_ttbar_Full]->SetDirectory(0);
            h_eta_nume[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr == _ttbar_700to1000 || pr == _ttbar_1000toInf)
        {
            h_pT_barrel_ctrl[_ttbar_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_ttbar_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_ttbar_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_ttbar_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_ttbar_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_ttbar_Full]->Add(h_eta_nume[pr]);
        }
        if (pr == _WJets)
        {
            h_pT_barrel_ctrl[_WJets_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_WJets")));
            h_pT_endcap_ctrl[_WJets_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_WJets")));
            h_pT_barrel_nume[_WJets_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_WJets")));
            h_pT_endcap_nume[_WJets_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_WJets")));
            h_eta_ctrl[_WJets_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_WJets")));
            h_eta_nume[_WJets_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_WJets")));

            h_pT_barrel_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_nume[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_nume[_WJets_Full]->SetDirectory(0);
            h_eta_ctrl[_WJets_Full]->SetDirectory(0);
            h_eta_nume[_WJets_Full]->SetDirectory(0);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_pT_barrel_ctrl[_WJets_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_WJets_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_WJets_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_WJets_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_WJets_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_WJets_Full]->Add(h_eta_nume[pr]);
        }


        file->Close();

        if (pr == _WW) {pr = _WZ; continue;}
        if (pr == _WZ) {pr = _ZZ; continue;}
        if (pr == _ZZ) {pr = _tbarW; continue;}
        if (pr == _tbarW) {pr = _tW; continue;}
        if (pr == _tW) {pr = _ttbar; continue;}
        if (pr == _ttbar) {pr = _ttbar_700to1000; continue;}
        if (pr == _ttbar_700to1000) {pr = _ttbar_1000toInf; continue;}
        if (pr == _ttbar_1000toInf) {pr = _WJets; continue;}
        if (pr == _WJets) {pr = _WJets_ext2v5; continue;}
        if (pr == _WJets_ext2v5) {stop = 1;}
    }

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);

        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);

        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_eta_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_DY")));
            h_pT_endcap_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_DY")));
            h_pT_barrel_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_DY")));
            h_pT_endcap_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_DY")));
            h_eta_ctrl[_DY_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_DY")));
            h_eta_nume[_DY_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_DY")));
            h_pT_barrel_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_nume[_DY_Full]->SetDirectory(0);
            h_eta_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_ctrl[_DY_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_DY_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_DY_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_DY_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_DY_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_DY_Full]->Add(h_eta_nume[pr]);
        }
        file->Close();
    }

    // Gamma+Jets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);

        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);

        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_eta_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]->SetDirectory(0);

        if (pr == _GJets_20to100)
        {
            h_pT_barrel_ctrl[_GJets_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_GJets")));
            h_pT_endcap_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_GJets")));
            h_pT_barrel_nume[_GJets_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_GJets")));
            h_pT_endcap_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_GJets")));
            h_eta_ctrl[_GJets_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_GJets")));
            h_eta_nume[_GJets_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_GJets")));

            h_pT_barrel_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_nume[_GJets_Full]->SetDirectory(0);
            h_eta_ctrl[_GJets_Full]->SetDirectory(0);
            h_eta_nume[_GJets_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_ctrl[_GJets_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_GJets_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_GJets_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_GJets_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_GJets_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_GJets_Full]->Add(h_eta_nume[pr]);
        }
        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);

        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);

        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_eta_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]->SetDirectory(0);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_pT_barrel_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_QCD")));
            h_pT_endcap_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_QCD")));
            h_pT_barrel_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_QCD")));
            h_pT_endcap_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_QCD")));
            h_eta_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_QCD")));
            h_eta_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_QCD")));

            h_pT_barrel_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_nume[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_ctrl[_QCDEMEnriched_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_QCDEMEnriched_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_QCDEMEnriched_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_QCDEMEnriched_Full]->Add(h_eta_nume[pr]);
        }
        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_DoubleEG_B; pr<=_DoubleEG_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);

        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);

        if (pr == _DoubleEG_B)
        {
            h_pT_barrel_ctrl[_DoubleEG_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_data")));
            h_pT_endcap_ctrl[_DoubleEG_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_data")));
            h_pT_barrel_nume[_DoubleEG_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_data")));
            h_pT_endcap_nume[_DoubleEG_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_data")));
            h_eta_ctrl[_DoubleEG_Full] = ((TH1D*)(h_eta_ctrl[pr]->Clone("h_eta_ctrl_data")));
            h_eta_nume[_DoubleEG_Full] = ((TH1D*)(h_eta_nume[pr]->Clone("h_eta_nume_data")));

            h_pT_barrel_ctrl[_DoubleEG_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_DoubleEG_Full]->SetDirectory(0);
            h_pT_barrel_nume[_DoubleEG_Full]->SetDirectory(0);
            h_pT_endcap_nume[_DoubleEG_Full]->SetDirectory(0);
            h_eta_ctrl[_DoubleEG_Full]->SetDirectory(0);
            h_eta_nume[_DoubleEG_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_ctrl[_DoubleEG_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_DoubleEG_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_pT_barrel_nume[_DoubleEG_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_DoubleEG_Full]->Add(h_pT_endcap_nume[pr]);
            h_eta_ctrl[_DoubleEG_Full]->Add(h_eta_ctrl[pr]);
            h_eta_nume[_DoubleEG_Full]->Add(h_eta_nume[pr]);
        }
    }


//--------------------------------- FR from DY MC -------------------------------------- (deno = nume + ctrl)

    //            DY_MC_nume
    // FR = -----------------------
    //      DY_MC_nume + DY_MC_ctrl

    // ------ Numerator ------ //
    TH1D *h_pT_barrel_nume_fMC = ((TH1D*)(h_pT_barrel_nume[_DY_Full]->Clone("h_pT_barrel_nume_fMC")));
    TH1D *h_pT_endcap_nume_fMC = ((TH1D*)(h_pT_endcap_nume[_DY_Full]->Clone("h_pT_endcap_nume_fMC")));
    TH1D *h_eta_nume_fMC = ((TH1D*)(h_eta_nume[_DY_Full]->Clone("h_eta_nume_fMC")));

    // ------ Denominator ------ //
    TH1D *h_pT_barrel_deno_fMC = ((TH1D*)(h_pT_barrel_ctrl[_DY_Full]->Clone("h_pT_barrel_deno_fMC")));
    h_pT_barrel_deno_fMC->Add(h_pT_barrel_nume_fMC); // deno = sig+ctrl
    TH1D *h_pT_endcap_deno_fMC = ((TH1D*)(h_pT_endcap_ctrl[_DY_Full]->Clone("h_pT_endcap_deno_fMC")));
    h_pT_endcap_deno_fMC->Add(h_pT_endcap_nume_fMC); // deno = sig+ctrl
    TH1D *h_eta_deno_fMC = ((TH1D*)(h_eta_ctrl[_DY_Full]->Clone("h_eta_deno_fMC")));
    h_eta_deno_fMC->Add(h_eta_nume_fMC); // deno = sig+ctrl

    // ------ FR ------ //
    // Barrel
    TH1D *h_FRMC_barrel = ((TH1D*)(h_pT_barrel_nume_fMC->Clone("h_FRMC_barrel")));
    h_FRMC_barrel->Divide(h_pT_barrel_deno_fMC);
    h_FRMC_barrel->SetDirectory(0);
    cout << "Numerator barrel (DY MC): " << h_pT_barrel_nume_fMC->Integral() << endl;
    cout << "Denominator barrel (DY MC): " << h_pT_barrel_deno_fMC->Integral() << endl;
    // Endcap
    TH1D *h_FRMC_endcap = ((TH1D*)(h_pT_endcap_nume_fMC->Clone("h_FRMC_endcap")));
    h_FRMC_endcap->Divide(h_pT_endcap_deno_fMC);
    h_FRMC_endcap->SetDirectory(0);
    cout << "Numerator endcap (DY MC): " << h_pT_endcap_nume_fMC->Integral() << endl;
    cout << "Denominator endcap (DY MC): " << h_pT_endcap_deno_fMC->Integral() << endl;
    // Eta
    TH1D *h_FRMC_eta = ((TH1D*)(h_eta_nume_fMC->Clone("h_FRMC_eta")));
    h_FRMC_eta->Divide(h_eta_deno_fMC);
    h_FRMC_eta->SetDirectory(0);

//--------------------------------- FR by ratio -------------------------------------- (deno = nume + ctrl)

    //            DATA_nume * DY_nume * sum(allMC_nume + allMC_ctrl)
    // FR = -----------------------------------------------------------------
    //      (DATA_nume + DATA_ctrl) * (DY_nume + DY_ctrl) * sum(allMC_nume)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_div = ((TH1D*)(h_pT_barrel_nume[_DY_Full]->Clone("h_pT_barrel_nume_div")));
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_ttbar_Full]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_tW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_tbarW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_WW]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_WZ]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_ZZ]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_WJets_Full]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_GJets_Full]);
    h_pT_barrel_nume_div->Add(h_pT_barrel_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_nume_full = ((TH1D*)(h_pT_barrel_nume[_DoubleEG_Full]->Clone("h_pT_barrel_nume")));
    h_pT_barrel_nume_full->Multiply(h_pT_barrel_nume[_DY_Full]);
    h_pT_barrel_nume_full->Divide(h_pT_barrel_nume_div);

    // Endcap
    TH1D *h_pT_endcap_nume_div = ((TH1D*)(h_pT_endcap_nume[_DY_Full]->Clone("h_pT_endcap_nume_div")));
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_ttbar_Full]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_tW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_tbarW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_WW]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_WZ]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_ZZ]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_WJets_Full]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_GJets_Full]);
    h_pT_endcap_nume_div->Add(h_pT_endcap_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_nume_full = ((TH1D*)(h_pT_endcap_nume[_DoubleEG_Full]->Clone("h_pT_endcap_nume")));
    h_pT_endcap_nume_full->Multiply(h_pT_endcap_nume[_DY_Full]);
    h_pT_endcap_nume_full->Divide(h_pT_endcap_nume_div);

    // Eta
    TH1D *h_eta_nume_div = ((TH1D*)(h_eta_nume[_DY_Full]->Clone("h_eta_nume_div")));
    h_eta_nume_div->Add(h_eta_nume[_ttbar_Full]);
    h_eta_nume_div->Add(h_eta_nume[_tW]);
    h_eta_nume_div->Add(h_eta_nume[_tbarW]);
    h_eta_nume_div->Add(h_eta_nume[_WW]);
    h_eta_nume_div->Add(h_eta_nume[_WZ]);
    h_eta_nume_div->Add(h_eta_nume[_ZZ]);
    h_eta_nume_div->Add(h_eta_nume[_WJets_Full]);
    h_eta_nume_div->Add(h_eta_nume[_GJets_Full]);
    h_eta_nume_div->Add(h_eta_nume[_QCDEMEnriched_Full]);
    TH1D *h_eta_nume_full = ((TH1D*)(h_eta_nume[_DoubleEG_Full]->Clone("h_eta_nume")));
    h_eta_nume_full->Multiply(h_eta_nume[_DY_Full]);
    h_eta_nume_full->Divide(h_eta_nume_div);

    // ####### Control ####### //
    // Barrel
    TH1D *h_pT_barrel_ctrl_div = ((TH1D*)(h_pT_barrel_ctrl[_DY_Full]->Clone("h_pT_barrel_ctrl_div")));
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_ttbar_Full]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_tW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_tbarW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_WW]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_WZ]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_ZZ]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_WJets_Full]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_GJets_Full]);
    h_pT_barrel_ctrl_div->Add(h_pT_barrel_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_ctrl_full = ((TH1D*)(h_pT_barrel_ctrl[_DoubleEG_Full]->Clone("h_pT_barrel_ctrl")));
    h_pT_barrel_ctrl_full->Multiply(h_pT_barrel_ctrl[_DY_Full]);
    h_pT_barrel_ctrl_full->Divide(h_pT_barrel_ctrl_div);

    // Endcap
    TH1D *h_pT_endcap_ctrl_div = ((TH1D*)(h_pT_endcap_ctrl[_DY_Full]->Clone("h_pT_endcap_ctrl_div")));
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_ttbar_Full]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_tW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_tbarW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_WW]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_WZ]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_ZZ]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_WJets_Full]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_GJets_Full]);
    h_pT_endcap_ctrl_div->Add(h_pT_endcap_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_ctrl_full = ((TH1D*)(h_pT_endcap_ctrl[_DoubleEG_Full]->Clone("h_pT_endcap_ctrl")));
    h_pT_endcap_ctrl_full->Multiply(h_pT_endcap_ctrl[_DY_Full]);
    h_pT_endcap_ctrl_full->Divide(h_pT_endcap_ctrl_div);

    // Eta
    TH1D *h_eta_ctrl_div = ((TH1D*)(h_eta_ctrl[_DY_Full]->Clone("h_eta_ctrl_div")));
    h_eta_ctrl_div->Add(h_eta_ctrl[_ttbar_Full]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_tW]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_tbarW]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_WW]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_WZ]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_ZZ]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_WJets_Full]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_GJets_Full]);
    h_eta_ctrl_div->Add(h_eta_ctrl[_QCDEMEnriched_Full]);
    TH1D *h_eta_ctrl_full = ((TH1D*)(h_eta_ctrl[_DoubleEG_Full]->Clone("h_eta_ctrl")));
    h_eta_ctrl_full->Multiply(h_eta_ctrl[_DY_Full]);
    h_eta_ctrl_full->Divide(h_eta_ctrl_div);

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_deno_div = ((TH1D*)(h_pT_barrel_ctrl[_DY_Full]->Clone("h_pT_barrel_deno_div")));
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_ttbar_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_tW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_tbarW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_WW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_WZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_ZZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_WJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_GJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_ctrl[_QCDEMEnriched_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_DY_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_ttbar_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_tW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_tbarW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_WW]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_WZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_ZZ]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_WJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_GJets_Full]);
    h_pT_barrel_deno_div->Add(h_pT_barrel_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_barrel_deno_mult = ((TH1D*)(h_pT_barrel_ctrl[_DY_Full]->Clone("h_pT_barrel_deno_mult")));
    h_pT_barrel_deno_mult->Add(h_pT_barrel_nume[_DY_Full]);
    TH1D *h_pT_barrel_deno_full = ((TH1D*)(h_pT_barrel_ctrl[_DoubleEG_Full]->Clone("h_pT_barrel_deno")));
    h_pT_barrel_deno_full->Add(h_pT_barrel_nume[_DoubleEG_Full]);
    h_pT_barrel_deno_full->Multiply(h_pT_barrel_deno_mult);
    h_pT_barrel_deno_full->Divide(h_pT_barrel_deno_div);

    // Endcap
    TH1D *h_pT_endcap_deno_div = ((TH1D*)(h_pT_endcap_ctrl[_DY_Full]->Clone("h_pT_endcap_deno_div")));
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_ttbar_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_tW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_tbarW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_WW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_WZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_ZZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_WJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_GJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_ctrl[_QCDEMEnriched_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_DY_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_ttbar_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_tW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_tbarW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_WW]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_WZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_ZZ]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_WJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_GJets_Full]);
    h_pT_endcap_deno_div->Add(h_pT_endcap_nume[_QCDEMEnriched_Full]);
    TH1D *h_pT_endcap_deno_mult = ((TH1D*)(h_pT_endcap_ctrl[_DY_Full]->Clone("h_pT_endcap_deno_mult")));
    h_pT_endcap_deno_mult->Add(h_pT_endcap_nume[_DY_Full]);
    TH1D *h_pT_endcap_deno_full = ((TH1D*)(h_pT_endcap_ctrl[_DoubleEG_Full]->Clone("h_pT_endcap_deno")));
    h_pT_endcap_deno_full->Add(h_pT_endcap_nume[_DoubleEG_Full]);
    h_pT_endcap_deno_full->Multiply(h_pT_endcap_deno_mult);
    h_pT_endcap_deno_full->Divide(h_pT_endcap_deno_div);

    // Eta
    TH1D *h_eta_deno_div = ((TH1D*)(h_eta_ctrl[_DY_Full]->Clone("h_eta_deno_div")));
    h_eta_deno_div->Add(h_eta_ctrl[_ttbar]);
    h_eta_deno_div->Add(h_eta_ctrl[_tW]);
    h_eta_deno_div->Add(h_eta_ctrl[_tbarW]);
    h_eta_deno_div->Add(h_eta_ctrl[_WW]);
    h_eta_deno_div->Add(h_eta_ctrl[_WZ]);
    h_eta_deno_div->Add(h_eta_ctrl[_ZZ]);
    h_eta_deno_div->Add(h_eta_ctrl[_WJets]);
    h_eta_deno_div->Add(h_eta_ctrl[_GJets_Full]);
    h_eta_deno_div->Add(h_eta_ctrl[_QCDEMEnriched_Full]);
    h_eta_deno_div->Add(h_eta_nume[_DY_Full]);
    h_eta_deno_div->Add(h_eta_nume[_ttbar]);
    h_eta_deno_div->Add(h_eta_nume[_tW]);
    h_eta_deno_div->Add(h_eta_nume[_tbarW]);
    h_eta_deno_div->Add(h_eta_nume[_WW]);
    h_eta_deno_div->Add(h_eta_nume[_WZ]);
    h_eta_deno_div->Add(h_eta_nume[_ZZ]);
    h_eta_deno_div->Add(h_eta_nume[_WJets]);
    h_eta_deno_div->Add(h_eta_nume[_GJets_Full]);
    h_eta_deno_div->Add(h_eta_nume[_QCDEMEnriched_Full]);
    TH1D *h_eta_deno_mult = ((TH1D*)(h_eta_ctrl[_DY_Full]->Clone("h_eta_deno_mult")));
    h_eta_deno_mult->Add(h_eta_nume[_DY_Full]);
    TH1D *h_eta_deno_full = ((TH1D*)(h_eta_ctrl[_DoubleEG_Full]->Clone("h_eta_deno")));
    h_eta_deno_full->Add(h_eta_nume[_DoubleEG_Full]);
    h_eta_deno_full->Multiply(h_eta_deno_mult);
    h_eta_deno_full->Divide(h_eta_deno_div);

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRratio_barrel = ((TH1D*)(h_pT_barrel_nume_full->Clone("h_FRratio_barrel")));
    h_FRratio_barrel->Divide(h_pT_barrel_deno_full);
    h_FRratio_barrel->SetDirectory(0);
    cout << "Numerator barrel (ratio): " << h_pT_barrel_nume_full->Integral() << endl;
    cout << "Denominator barrel (ratio): " << h_pT_barrel_deno_full->Integral() << endl;
    // Endcap
    TH1D *h_FRratio_endcap = ((TH1D*)(h_pT_endcap_nume_full->Clone("h_FRratio_endcap")));
    h_FRratio_endcap->Divide(h_pT_endcap_deno_full);
    h_FRratio_endcap->SetDirectory(0);
    cout << "Numerator endcap (ratio): " << h_pT_endcap_nume_full->Integral() << endl;
    cout << "Denominator endcap (ratio): " << h_pT_endcap_deno_full->Integral() << endl;
    // Eta
    TH1D *h_FRratio_eta = ((TH1D*)(h_eta_nume_full->Clone("h_FRratio_eta")));
    h_FRratio_eta->Divide(h_eta_deno_full);
    h_FRratio_eta->SetDirectory(0);


//--------------------------------- FR by subtraction -------------------------------------- (deno = nume + ctrl)

    //                     DATA_nume - sum(prompt_MC_nume)
    // FR = --------------------------------------------------------------
    //      (DATA_nume + DATA_ctrl) - sum(prompt_MC_nume + prompt_MC_ctrl)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_sub = ((TH1D*)(h_pT_barrel_nume[_DoubleEG_Full]->Clone("h_pT_barrel_nume_sub")));
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_ttbar_Full], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_tW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_tbarW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_WW], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_WZ], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_ZZ], -1);
    h_pT_barrel_nume_sub->Add(h_pT_barrel_nume[_GJets_Full], -1);
    // Endcap
    TH1D *h_pT_endcap_nume_sub = ((TH1D*)(h_pT_endcap_nume[_DoubleEG_Full]->Clone("h_pT_endcap_nume_sub")));
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_ttbar_Full], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_tW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_tbarW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_WW], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_WZ], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_ZZ], -1);
    h_pT_endcap_nume_sub->Add(h_pT_endcap_nume[_GJets_Full], -1);
    // Eta
    TH1D *h_eta_nume_sub = ((TH1D*)(h_eta_nume[_DoubleEG_Full]->Clone("h_eta_nume_sub")));
    h_eta_nume_sub->Add(h_eta_nume[_ttbar_Full], -1);
    h_eta_nume_sub->Add(h_eta_nume[_tW], -1);
    h_eta_nume_sub->Add(h_eta_nume[_tbarW], -1);
    h_eta_nume_sub->Add(h_eta_nume[_WW], -1);
    h_eta_nume_sub->Add(h_eta_nume[_WZ], -1);
    h_eta_nume_sub->Add(h_eta_nume[_ZZ], -1);
    h_eta_nume_sub->Add(h_eta_nume[_GJets_Full], -1);

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_ctrl_sub = ((TH1D*)(h_pT_barrel_ctrl[_DoubleEG_Full]->Clone("h_pT_barrel_ctrl_sub")));
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_ttbar_Full], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_tW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_tbarW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_WW], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_WZ], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_ZZ], -1);
    h_pT_barrel_ctrl_sub->Add(h_pT_barrel_ctrl[_GJets_Full], -1);
    TH1D *h_pT_barrel_deno_sub = ((TH1D*)(h_pT_barrel_ctrl_sub->Clone("h_pT_barrel_deno_sub")));
    h_pT_barrel_deno_sub->Add(h_pT_barrel_nume_sub); // deno = sig+ctrl
    // Endcap
    TH1D *h_pT_endcap_ctrl_sub = ((TH1D*)(h_pT_endcap_ctrl[_DoubleEG_Full]->Clone("h_pT_endcap_ctrl_sub")));
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_ttbar_Full], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_tW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_tbarW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_WW], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_WZ], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_ZZ], -1);
    h_pT_endcap_ctrl_sub->Add(h_pT_endcap_ctrl[_GJets_Full], -1);
    TH1D *h_pT_endcap_deno_sub = ((TH1D*)(h_pT_endcap_ctrl_sub->Clone("h_pT_endcap_deno_sub")));
    h_pT_endcap_deno_sub->Add(h_pT_endcap_nume_sub); // deno = sig+ctrl
    // Eta
    TH1D *h_eta_ctrl_sub = ((TH1D*)(h_eta_ctrl[_DoubleEG_Full]->Clone("h_eta_ctrl_sub")));
    h_eta_ctrl_sub->Add(h_eta_ctrl[_ttbar_Full], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_tW], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_tbarW], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_WW], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_WZ], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_ZZ], -1);
    h_eta_ctrl_sub->Add(h_eta_ctrl[_GJets_Full], -1);
    TH1D *h_eta_deno_sub = ((TH1D*)(h_eta_ctrl_sub->Clone("h_eta_deno_sub")));
    h_eta_deno_sub->Add(h_eta_nume_sub); // deno = sig+ctrl

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRsubtract_barrel = ((TH1D*)(h_pT_barrel_nume_sub->Clone("h_FRsubtract_barrel")));
    h_FRsubtract_barrel->Divide(h_pT_barrel_deno_sub);
    h_FRsubtract_barrel->SetDirectory(0);
    cout << "Numerator barrel (subtraction): " << h_pT_barrel_nume_sub->Integral() << endl;
    cout << "Denominator barrel (subtraction): " << h_pT_barrel_deno_sub->Integral() << endl;
    // Endcap
    TH1D *h_FRsubtract_endcap = ((TH1D*)(h_pT_endcap_nume_sub->Clone("h_FRsubtract_endcap")));
    h_FRsubtract_endcap->Divide(h_pT_endcap_deno_sub);
    h_FRsubtract_endcap->SetDirectory(0);
    cout << "Numerator endcap (subtraction): " << h_pT_endcap_nume_sub->Integral() << endl;
    cout << "Denominator endcap (subtraction): " << h_pT_endcap_deno_sub->Integral() << endl;
    // Eta
    TH1D *h_FRsubtract_eta = ((TH1D*)(h_eta_nume_sub->Clone("h_FRsubtract_eta")));
    h_FRsubtract_eta->Divide(h_eta_deno_sub);
    h_FRsubtract_eta->SetDirectory(0);

//--------------------------------- FR by template -------------------------------------- (deno = nume + ctrl)
/*
    //                  QCD_nume(from ABCD)
    // FR = --------------------------------------------
    //      (QCD_nume + DATA_ctrl) - sum(nonQCD_MC_ctrl)

    // ####### Numerator ####### //
    // Barrel
    TH1D *h_pT_barrel_nume_abcd, *h_pT_endcap_nume_abcd;
    TFile *f_abcd = new TFile("/media/sf_DATA/FR/Electron/ABCD_hists.root", "READ");
    f_abcd->GetObject("h_pT_QCD_nume_barrel", h_pT_barrel_nume_abcd);
    f_abcd->GetObject("h_pT_QCD_nume_endcap", h_pT_endcap_nume_abcd);
    h_pT_barrel_nume_abcd->SetDirectory(0);
    h_pT_endcap_nume_abcd->SetDirectory(0);
    f_abcd->Close();

    // ####### Denominator ####### //
    // Barrel
    TH1D *h_pT_barrel_deno_abcd = ((TH1D*)(h_pT_barrel_data_ctrl->Clone("h_pT_barrel_deno_abcd")));
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_ttbar], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_tW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_tbarW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WW], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WZ], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_ZZ], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_WJets], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_GJets_Full], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_MC_ctrl[_DY_Full], -1);
    h_pT_barrel_deno_abcd->Add(h_pT_barrel_nume_abcd); // deno = sig+ctrl
    // Endcap
    TH1D *h_pT_endcap_deno_abcd = ((TH1D*)(h_pT_endcap_data_ctrl->Clone("h_pT_endcap_deno_abcd")));
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_ttbar], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_tW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_tbarW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WW], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WZ], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_ZZ], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_WJets], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_GJets_Full], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_MC_ctrl[_DY_Full], -1);
    h_pT_endcap_deno_abcd->Add(h_pT_endcap_nume_abcd); // deno = sig+ctrl

    // ######## FR ######## //
    // Barrel
    TH1D *h_FRtemplate_barrel = ((TH1D*)(h_pT_barrel_nume_abcd->Clone("h_FRtemplate_barrel")));
    h_FRtemplate_barrel->Divide(h_pT_barrel_deno_abcd);
    h_FRtemplate_barrel->SetDirectory(0);
    cout << "Numerator barrel (abcd): " << h_pT_barrel_nume_abcd->Integral() << endl;
    cout << "Denominator barrel (abcd): " << h_pT_barrel_deno_abcd->Integral() << endl;
    // Endcap
    TH1D *h_FRtemplate_endcap = ((TH1D*)(h_pT_endcap_nume_abcd->Clone("h_FRtemplate_endcap")));
    h_FRtemplate_endcap->Divide(h_pT_endcap_deno_abcd);
    h_FRtemplate_endcap->SetDirectory(0);
    cout << "Numerator endcap (abcd): " << h_pT_endcap_nume_abcd->Integral() << endl;
    cout << "Denominator endcap (abcd): " << h_pT_endcap_deno_abcd->Integral() << endl;
    */

// --------------------- Geting the right errors -------------------- //
//                           ___________________________________
//       /  A  \      A*B    |  / Delta(A) \^2   / Delta(B) \^2 |
// Delta( ----- ) = -------  | ( ---------- ) + ( ---------- )    ;    Here A = Signal region, B = Control region
//       \ A+B /    (A+B)^2 \/  \    A     /     \    B     /          (numerator=signal, denominator=signal+control)

    for (Int_t i_bin=1; i_bin<=nPtBin_ele_alt; i_bin++)
    {
        Double_t sig_barrel_MC = h_pT_barrel_nume[_DY_Full]->GetBinContent(i_bin);
        Double_t ctrl_barrel_MC = h_pT_barrel_ctrl[_DY_Full]->GetBinContent(i_bin);
        Double_t err_sig_barrel_MC = h_pT_barrel_nume[_DY_Full]->GetBinError(i_bin);
        Double_t err_ctrl_barrel_MC = h_pT_barrel_ctrl[_DY_Full]->GetBinError(i_bin);

        Double_t sig_endcap_MC = h_pT_barrel_nume[_DY_Full]->GetBinContent(i_bin);
        Double_t ctrl_endcap_MC = h_pT_barrel_ctrl[_DY_Full]->GetBinContent(i_bin);
        Double_t err_sig_endcap_MC = h_pT_barrel_nume[_DY_Full]->GetBinError(i_bin);
        Double_t err_ctrl_endcap_MC = h_pT_barrel_ctrl[_DY_Full]->GetBinError(i_bin);

        Double_t sig_barrel_ratio = h_pT_barrel_nume_full->GetBinContent(i_bin);
        Double_t ctrl_barrel_ratio = h_pT_barrel_ctrl_full->GetBinContent(i_bin);
        Double_t err_sig_barrel_ratio = h_pT_barrel_nume_full->GetBinError(i_bin);
        Double_t err_ctrl_barrel_ratio = h_pT_barrel_ctrl_full->GetBinError(i_bin);

        Double_t sig_endcap_ratio = h_pT_endcap_nume_full->GetBinContent(i_bin);
        Double_t ctrl_endcap_ratio = h_pT_endcap_ctrl_full->GetBinContent(i_bin);
        Double_t err_sig_endcap_ratio = h_pT_endcap_nume_full->GetBinError(i_bin);
        Double_t err_ctrl_endcap_ratio = h_pT_endcap_ctrl_full->GetBinError(i_bin);

        Double_t sig_barrel_sub = h_pT_barrel_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_barrel_sub = h_pT_barrel_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_barrel_sub = h_pT_barrel_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_barrel_sub = h_pT_barrel_ctrl_sub->GetBinError(i_bin);

        Double_t sig_endcap_sub = h_pT_endcap_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_endcap_sub = h_pT_endcap_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_endcap_sub = h_pT_endcap_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_endcap_sub = h_pT_endcap_ctrl_sub->GetBinError(i_bin);

        Double_t err_FR_barrel_MC = sqrt((err_sig_barrel_MC / sig_barrel_MC) * (err_sig_barrel_MC / sig_barrel_MC) + (err_ctrl_barrel_MC / ctrl_barrel_MC) * (err_ctrl_barrel_MC / ctrl_barrel_MC));
        err_FR_barrel_MC *= sig_barrel_MC * ctrl_barrel_MC / ((sig_barrel_MC + ctrl_barrel_MC) * (sig_barrel_MC + ctrl_barrel_MC));

        Double_t err_FR_endcap_MC = sqrt((err_sig_endcap_MC / sig_endcap_MC) * (err_sig_endcap_MC / sig_endcap_MC) + (err_ctrl_endcap_MC / ctrl_endcap_MC) * (err_ctrl_endcap_MC / ctrl_endcap_MC));
        err_FR_endcap_MC *= sig_endcap_MC * ctrl_endcap_MC / ((sig_endcap_MC + ctrl_endcap_MC) * (sig_endcap_MC + ctrl_endcap_MC));

        Double_t err_FR_barrel_ratio = sqrt((err_sig_barrel_ratio / sig_barrel_ratio) * (err_sig_barrel_ratio / sig_barrel_ratio) + (err_ctrl_barrel_ratio / ctrl_barrel_ratio) * (err_ctrl_barrel_ratio / ctrl_barrel_ratio));
        err_FR_barrel_ratio *= sig_barrel_ratio * ctrl_barrel_ratio / ((sig_barrel_ratio + ctrl_barrel_ratio) * (sig_barrel_ratio + ctrl_barrel_ratio));

        Double_t err_FR_endcap_ratio = sqrt((err_sig_endcap_ratio / sig_endcap_ratio) * (err_sig_endcap_ratio / sig_endcap_ratio) + (err_ctrl_endcap_ratio / ctrl_endcap_ratio) * (err_ctrl_endcap_ratio / ctrl_endcap_ratio));
        err_FR_endcap_ratio *= sig_endcap_ratio * ctrl_endcap_ratio / ((sig_endcap_ratio + ctrl_endcap_ratio) * (sig_endcap_ratio + ctrl_endcap_ratio));

        Double_t err_FR_barrel_sub = sqrt((err_sig_barrel_sub / sig_barrel_sub) * (err_sig_barrel_sub / sig_barrel_sub) + (err_ctrl_barrel_sub / ctrl_barrel_sub) * (err_ctrl_barrel_sub / ctrl_barrel_sub));
        err_FR_barrel_sub *= sig_barrel_sub * ctrl_barrel_sub / ((sig_barrel_sub + ctrl_barrel_sub) * (sig_barrel_sub + ctrl_barrel_sub));

        Double_t err_FR_endcap_sub = sqrt((err_sig_endcap_sub / sig_endcap_sub) * (err_sig_endcap_sub / sig_endcap_sub) + (err_ctrl_endcap_sub / ctrl_endcap_sub) * (err_ctrl_endcap_sub / ctrl_endcap_sub));
        err_FR_endcap_sub *= sig_endcap_sub * ctrl_endcap_sub / ((sig_endcap_sub + ctrl_endcap_sub) * (sig_endcap_sub + ctrl_endcap_sub));

        h_FRMC_barrel->SetBinError(i_bin, err_FR_barrel_MC);
        h_FRMC_endcap->SetBinError(i_bin, err_FR_endcap_MC);
        h_FRratio_barrel->SetBinError(i_bin, err_FR_barrel_ratio);
        h_FRratio_endcap->SetBinError(i_bin, err_FR_endcap_ratio);
        h_FRsubtract_barrel->SetBinError(i_bin, err_FR_barrel_sub);
        h_FRsubtract_endcap->SetBinError(i_bin, err_FR_endcap_sub);
    }

    for (Int_t i_bin=1; i_bin<=50; i_bin++)
    {
        Double_t sig_eta_MC = h_eta_nume[_DY_Full]->GetBinContent(i_bin);
        Double_t ctrl_eta_MC = h_eta_ctrl[_DY_Full]->GetBinContent(i_bin);
        Double_t err_sig_eta_MC = h_eta_nume[_DY_Full]->GetBinError(i_bin);
        Double_t err_ctrl_eta_MC = h_eta_ctrl[_DY_Full]->GetBinError(i_bin);

        Double_t sig_eta_ratio = h_eta_nume_full->GetBinContent(i_bin);
        Double_t ctrl_eta_ratio = h_eta_ctrl_full->GetBinContent(i_bin);
        Double_t err_sig_eta_ratio = h_eta_nume_full->GetBinError(i_bin);
        Double_t err_ctrl_eta_ratio = h_eta_ctrl_full->GetBinError(i_bin);

        Double_t sig_eta_sub = h_eta_nume_sub->GetBinContent(i_bin);
        Double_t ctrl_eta_sub = h_eta_ctrl_sub->GetBinContent(i_bin);
        Double_t err_sig_eta_sub = h_eta_nume_sub->GetBinError(i_bin);
        Double_t err_ctrl_eta_sub = h_eta_ctrl_sub->GetBinError(i_bin);

        Double_t err_FR_eta_MC = sqrt((err_sig_eta_MC / sig_eta_MC) * (err_sig_eta_MC / sig_eta_MC) + (err_ctrl_eta_MC / ctrl_eta_MC) * (err_ctrl_eta_MC / ctrl_eta_MC));
        err_FR_eta_MC *= sig_eta_MC * ctrl_eta_MC / ((sig_eta_MC + ctrl_eta_MC) * (sig_eta_MC + ctrl_eta_MC));

        Double_t err_FR_eta_ratio = sqrt((err_sig_eta_ratio / sig_eta_ratio) * (err_sig_eta_ratio / sig_eta_ratio) + (err_ctrl_eta_ratio / ctrl_eta_ratio) * (err_ctrl_eta_ratio / ctrl_eta_ratio));
        err_FR_eta_ratio *= sig_eta_ratio * ctrl_eta_ratio / ((sig_eta_ratio + ctrl_eta_ratio) * (sig_eta_ratio + ctrl_eta_ratio));

        Double_t err_FR_eta_sub = sqrt((err_sig_eta_sub / sig_eta_sub) * (err_sig_eta_sub / sig_eta_sub) + (err_ctrl_eta_sub / ctrl_eta_sub) * (err_ctrl_eta_sub / ctrl_eta_sub));
        err_FR_eta_sub *= sig_eta_sub * ctrl_eta_sub / ((sig_eta_sub + ctrl_eta_sub) * (sig_eta_sub + ctrl_eta_sub));

        if (h_FRMC_eta->GetBinContent(i_bin) > 0)
            h_FRMC_eta->SetBinError(i_bin, err_FR_eta_MC);
        else h_FRMC_eta->SetBinError(i_bin, 0);
        if (h_FRratio_eta->GetBinContent(i_bin) > 0)
            h_FRratio_eta->SetBinError(i_bin, err_FR_eta_ratio);
        else h_FRratio_eta->SetBinError(i_bin, 0);
        if (h_FRsubtract_eta->GetBinContent(i_bin) > 0)
            h_FRsubtract_eta->SetBinError(i_bin, err_FR_eta_sub);
        else h_FRsubtract_eta->SetBinError(i_bin, 0);
    }

// ---------------------------- Writing------------------------------ //
    TString filename = "/media/sf_DATA/FR/Electron/FakeRate_electron_alt.root";
    TFile *file_FR = new TFile(filename, "RECREATE");
    if (file_FR->IsOpen()) cout << "File '" << filename << "' has been created. Writing histograms.." << endl;
    file_FR->cd();
    h_FRMC_barrel->Write();
    h_FRMC_endcap->Write();
    h_FRratio_barrel->Write();
    h_FRratio_endcap->Write();
    h_FRsubtract_barrel->Write();
    h_FRsubtract_endcap->Write();
//    h_FRtemplate_barrel->Write();
//    h_FRtemplate_endcap->Write();
    cout << "Finished. Closing the file.." << endl;
    file_FR->Close();
    if (!file_FR->IsOpen()) cout << "File '" << filename << "' has been closed successfully." << endl;
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
    h_FRMC_barrel->SetMarkerStyle(23);
    h_FRMC_barrel->SetMarkerColor(kYellow);
    h_FRMC_barrel->SetLineColor(kYellow);
    h_FRMC_barrel->SetStats(kFALSE);
    h_FRsubtract_barrel->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_barrel->SetMarkerColor(kBlack);
    h_FRsubtract_barrel->SetLineColor(kBlack);
    h_FRsubtract_barrel->SetStats(kFALSE);
    h_FRratio_barrel->SetMarkerStyle(kFullSquare);
    h_FRratio_barrel->SetMarkerColor(kRed);
    h_FRratio_barrel->SetLineColor(kRed);
    h_FRratio_barrel->SetStats(kFALSE);
    /*h_FRtemplate_barrel->SetMarkerStyle(33);
    h_FRtemplate_barrel->SetMarkerSize(1.5);
    h_FRtemplate_barrel->SetMarkerColor(kGreen+2);
    h_FRtemplate_barrel->SetLineColor(kGreen+2);
    h_FRtemplate_barrel->SetStats(kFALSE);*/
    h_FRsubtract_barrel->SetTitle("");
    h_FRsubtract_barrel->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_FRsubtract_barrel->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_barrel->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_barrel->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_barrel->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_barrel->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_barrel->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_barrel->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_barrel->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_barrel->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_barrel->GetXaxis()->SetRangeUser(25, 3000);
    h_FRsubtract_barrel->GetYaxis()->SetRangeUser(0, 0.3);
    h_FRsubtract_barrel->Draw();
//    h_FRratio_barrel->Draw("same");
//    h_FRMC_barrel->Draw("same");
//    h_FRtemplate_barrel->Draw("same");

    TLegend *legend = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend->AddEntry(h_FRMC_barrel, "QCD MC", "LP");
//    legend->AddEntry(h_FRratio_barrel, "Ratio", "LP");
    legend->AddEntry(h_FRsubtract_barrel, "Subtraction", "LP");
    TLegend *legend_noABCD = ((TLegend*)legend->Clone());
//    legend->AddEntry(h_FRtemplate_barrel, "ABCD", "LP");
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
    h_FRMC_endcap->SetMarkerStyle(23);
    h_FRMC_endcap->SetMarkerColor(kYellow);
    h_FRMC_endcap->SetLineColor(kYellow);
    h_FRMC_endcap->SetStats(kFALSE);
    h_FRsubtract_endcap->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_endcap->SetMarkerColor(kBlack);
    h_FRsubtract_endcap->SetLineColor(kBlack);
    h_FRsubtract_endcap->SetStats(kFALSE);
    h_FRratio_endcap->SetMarkerStyle(kFullSquare);
    h_FRratio_endcap->SetMarkerColor(kRed);
    h_FRratio_endcap->SetLineColor(kRed);
    h_FRratio_endcap->SetStats(kFALSE);
    /*h_FRtemplate_endcap->SetMarkerStyle(33);
    h_FRtemplate_endcap->SetMarkerSize(1.5);
    h_FRtemplate_endcap->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap->SetLineColor(kGreen+2);
    h_FRtemplate_endcap->SetStats(kFALSE);*/
    h_FRsubtract_endcap->SetTitle("");
    h_FRsubtract_endcap->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_FRsubtract_endcap->GetXaxis()->SetTitleOffset(1);
    h_FRsubtract_endcap->GetXaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap->GetXaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap->GetYaxis()->SetTitle("Fake rate");
    h_FRsubtract_endcap->GetYaxis()->SetTitleSize(0.05);
    h_FRsubtract_endcap->GetYaxis()->SetTitleOffset(1.25);
    h_FRsubtract_endcap->GetYaxis()->SetLabelSize(0.04);
    h_FRsubtract_endcap->GetXaxis()->SetNoExponent(1);
    h_FRsubtract_endcap->GetXaxis()->SetMoreLogLabels(1);
    h_FRsubtract_endcap->GetXaxis()->SetRangeUser(25, 3000);
    h_FRsubtract_endcap->GetYaxis()->SetRangeUser(0, 0.3);
    h_FRsubtract_endcap->Draw();
//    h_FRratio_endcap->Draw("same");
//    h_FRMC_endcap->Draw("same");
//    h_FRtemplate_endcap->Draw("same");
    legend->Draw();
    TText *texte = new TText (0.45, 0.6, "Endcap");
    texte->SetTextAlign(11);
    texte->SetTextSize(0.05);
    texte->SetNDC(true);
    texte->Draw();
    c_FR_endcap->Update();

    TCanvas *c_FR_allin1 = new TCanvas("c_FR_allin1", "c_FR_allin1", 800, 800);
    c_FR_allin1->cd();
    c_FR_allin1->SetGrid(1);
    c_FR_allin1->SetRightMargin(0.05);
    c_FR_allin1->SetTopMargin(0.05);
    c_FR_allin1->SetBottomMargin(0.12);
    c_FR_allin1->SetLeftMargin(0.13);
    TH1D *h_FRsubtract_endcap1 = ((TH1D*)(h_FRsubtract_endcap->Clone("h_FRsubtract_endcap1")));
    h_FRsubtract_endcap1->SetMarkerStyle(kFullSquare);
    h_FRsubtract_endcap1->SetMarkerColor(kBlue);
    h_FRsubtract_endcap1->GetYaxis()->SetTitle("Signalo srities tikimyb#dot{e}");
    h_FRsubtract_endcap1->GetYaxis()->SetTitleSize(0.045);
    h_FRsubtract_endcap1->GetYaxis()->SetTitleOffset(1.4);
    h_FRsubtract_endcap1->Draw();
    h_FRsubtract_barrel->Draw("same");
    TLegend *legend1 = new TLegend(0.13, 0.77, 0.6, 0.95);
    legend1->AddEntry(h_FRsubtract_barrel, "Atimties metodas (|#eta_{#lower[-0.2]{SC}}| < 1.4442)", "LP");
    legend1->AddEntry(h_FRsubtract_endcap1, "Atimties metodas (|#eta_{#lower[-0.2]{SC}}| > 1.566)", "LP");
    legend1->Draw();
    c_FR_allin1->SetLogx();
    c_FR_allin1->Update();

    TCanvas *c_FR_eta = new TCanvas("c_FR_eta", "c_FR_eta", 800, 800);
    c_FR_eta->cd();
    c_FR_eta->SetGrid(1);
    c_FR_eta->SetRightMargin(0.05);
    c_FR_eta->SetTopMargin(0.05);
    c_FR_eta->SetBottomMargin(0.12);
    c_FR_eta->SetLeftMargin(0.13);
    h_FRMC_eta->SetMarkerStyle(23);
    h_FRMC_eta->SetMarkerColor(kYellow);
    h_FRMC_eta->SetLineColor(kYellow);
    h_FRMC_eta->SetStats(kFALSE);
    h_FRsubtract_eta->SetMarkerStyle(kFullDotLarge);
    h_FRsubtract_eta->SetMarkerColor(kBlack);
    h_FRsubtract_eta->SetLineColor(kBlack);
    h_FRsubtract_eta->SetStats(kFALSE);
    h_FRratio_eta->SetMarkerStyle(kFullSquare);
    h_FRratio_eta->SetMarkerColor(kRed);
    h_FRratio_eta->SetLineColor(kRed);
    h_FRratio_eta->SetStats(kFALSE);
    h_FRratio_eta->SetTitle("");
    h_FRratio_eta->GetXaxis()->SetTitle("#eta (#mu)");
    h_FRratio_eta->GetXaxis()->SetTitleOffset(1);
    h_FRratio_eta->GetXaxis()->SetTitleSize(0.05);
    h_FRratio_eta->GetXaxis()->SetLabelSize(0.04);
    h_FRratio_eta->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_eta->GetYaxis()->SetTitleSize(0.05);
    h_FRratio_eta->GetYaxis()->SetTitleOffset(1.25);
    h_FRratio_eta->GetYaxis()->SetLabelSize(0.04);
    h_FRratio_eta->GetXaxis()->SetNoExponent(1);
    h_FRratio_eta->GetXaxis()->SetMoreLogLabels(1);
    h_FRratio_eta->GetXaxis()->SetRangeUser(25, 1000);
    h_FRratio_eta->GetYaxis()->SetRangeUser(0, 0.1);
    h_FRratio_eta->Draw();
    h_FRMC_eta->Draw("same");
    h_FRsubtract_eta->Draw("same");
    legend_noABCD->Draw();
    c_FR_eta->Update();

} // End of E_EstFR_alt2()


void E_EstFRandPR_MC()
{
    TString filename = "/media/sf_DATA/FR/Electron/MCFR_Hists_E.root";
    TFile *f = new TFile(filename, "READ");
    cout << "Reading from file " << filename << " ... ";

    TH1D *h_FR_pT_barrel_nume[4], *h_FR_pT_endcap_nume[4], *h_FR_pT_endcap2_nume[4];
    TH1D *h_FR_pT_barrel_deno[4], *h_FR_pT_endcap_deno[4], *h_FR_pT_endcap2_deno[4];
    TH1D *h_FR_eta_nume[4],       *h_FR_eta_deno[4];
    TH1D *h_PR_pT_barrel_nume[4], *h_PR_pT_endcap_nume[4], *h_PR_pT_endcap2_nume[4];
    TH1D *h_PR_pT_barrel_deno[4], *h_PR_pT_endcap_deno[4], *h_PR_pT_endcap2_deno[4];
    TH1D *h_PR_eta_nume[4],       *h_PR_eta_deno[4];
    TH1D *h_FR_nume_sum[4],       *h_FR_deno_sum[4],        *h_PR_nume_sum[4],       *h_PR_deno_sum[4];

    TString type[4] = {"DY", "ttbar", "WJets", "QCD"};
    Color_t colors[4] = {kOrange, kCyan+2, kRed-2, kRed+3};
    Style_t markers[4] = {kFullDotLarge, kFullSquare, kFullDiamond, kFullTriangleDown};
    Float_t sizes[4] = {1, 1, 1.5, 1};

    for (Int_t i=0; i<4; i++)
    {
        f->GetObject("h_FR_pT_barrel_nume_"+type[i],  h_FR_pT_barrel_nume[i]);
        f->GetObject("h_FR_pT_endcap_nume_"+type[i],  h_FR_pT_endcap_nume[i]);
        f->GetObject("h_FR_pT_endcap2_nume_"+type[i], h_FR_pT_endcap2_nume[i]);
        f->GetObject("h_FR_eta_nume_"+type[i],        h_FR_eta_nume[i]);
        f->GetObject("h_FR_pT_barrel_deno_"+type[i],  h_FR_pT_barrel_deno[i]);
        f->GetObject("h_FR_pT_endcap_deno_"+type[i],  h_FR_pT_endcap_deno[i]);
        f->GetObject("h_FR_pT_endcap2_deno_"+type[i], h_FR_pT_endcap2_deno[i]);
        f->GetObject("h_FR_eta_deno_"+type[i],        h_FR_eta_deno[i]);
        f->GetObject("h_PR_pT_barrel_nume_"+type[i],  h_PR_pT_barrel_nume[i]);
        f->GetObject("h_PR_pT_endcap_nume_"+type[i],  h_PR_pT_endcap_nume[i]);
        f->GetObject("h_PR_pT_endcap2_nume_"+type[i], h_PR_pT_endcap2_nume[i]);
        f->GetObject("h_PR_eta_nume_"+type[i],        h_PR_eta_nume[i]);
        f->GetObject("h_PR_pT_barrel_deno_"+type[i],  h_PR_pT_barrel_deno[i]);
        f->GetObject("h_PR_pT_endcap_deno_"+type[i],  h_PR_pT_endcap_deno[i]);
        f->GetObject("h_PR_pT_endcap2_deno_"+type[i], h_PR_pT_endcap2_deno[i]);
        f->GetObject("h_PR_eta_deno_"+type[i],        h_PR_eta_deno[i]);
        h_FR_pT_barrel_nume[i]->SetDirectory(0);
        h_FR_pT_endcap_nume[i]->SetDirectory(0);
        h_FR_pT_endcap2_nume[i]->SetDirectory(0);
        h_FR_eta_nume[i]->SetDirectory(0);
        h_FR_pT_barrel_deno[i]->SetDirectory(0);
        h_FR_pT_endcap_deno[i]->SetDirectory(0);
        h_FR_pT_endcap2_deno[i]->SetDirectory(0);
        h_FR_eta_deno[i]->SetDirectory(0);
        h_PR_pT_barrel_nume[i]->SetDirectory(0);
        h_PR_pT_endcap_nume[i]->SetDirectory(0);
        h_PR_pT_endcap2_nume[i]->SetDirectory(0);
        h_PR_eta_nume[i]->SetDirectory(0);
        h_PR_pT_barrel_deno[i]->SetDirectory(0);
        h_PR_pT_endcap_deno[i]->SetDirectory(0);
        h_PR_pT_endcap2_deno[i]->SetDirectory(0);
        h_PR_eta_deno[i]->SetDirectory(0);

        removeNegativeBins(h_FR_pT_barrel_nume[i] );
        removeNegativeBins(h_FR_pT_endcap_nume[i] );
        removeNegativeBins(h_FR_pT_endcap2_nume[i]);
        removeNegativeBins(h_FR_eta_nume[i]       );
        removeNegativeBins(h_PR_pT_barrel_nume[i] );
        removeNegativeBins(h_PR_pT_endcap_nume[i] );
        removeNegativeBins(h_PR_pT_endcap2_nume[i]);
        removeNegativeBins(h_PR_eta_nume[i]       );

        h_FR_pT_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_pT_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_pT_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_eta_nume[i]       ->SetLineColor(colors[i]);
        h_PR_pT_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_eta_nume[i]       ->SetLineColor(colors[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_eta_nume[i]       ->SetMarkerColor(colors[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_eta_nume[i]       ->SetMarkerColor(colors[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_eta_nume[i]       ->SetMarkerStyle(markers[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_eta_nume[i]       ->SetMarkerStyle(markers[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_eta_nume[i]       ->SetMarkerSize(sizes[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_eta_nume[i]       ->SetMarkerSize(sizes[i]);

        h_FR_pT_barrel_nume[i] ->SetStats(0);
        h_FR_pT_endcap_nume[i] ->SetStats(0);
        h_FR_pT_endcap2_nume[i]->SetStats(0);
        h_FR_eta_nume[i]       ->SetStats(0);
        h_PR_pT_barrel_nume[i] ->SetStats(0);
        h_PR_pT_endcap_nume[i] ->SetStats(0);
        h_PR_pT_endcap2_nume[i]->SetStats(0);
        h_PR_eta_nume[i]       ->SetStats(0);

        // Calculating rates
        // Average rates
        h_FR_nume_sum[i] = ((TH1D*)(h_FR_pT_barrel_nume[i]->Clone("h_FR_nume_sum_"+type[i])));
        h_FR_nume_sum[i]->SetDirectory(0);
        h_FR_nume_sum[i]->Add(h_FR_pT_endcap_nume[i]);
        h_FR_nume_sum[i]->Add(h_FR_pT_endcap2_nume[i]);
        h_FR_deno_sum[i] = ((TH1D*)(h_FR_pT_barrel_deno[i]->Clone("h_FR_deno_sum_"+type[i])));
        h_FR_deno_sum[i]->SetDirectory(0);
        h_FR_deno_sum[i]->Add(h_FR_pT_endcap_deno[i]);
        h_FR_deno_sum[i]->Add(h_FR_pT_endcap2_deno[i]);

        h_PR_nume_sum[i] = ((TH1D*)(h_PR_pT_barrel_nume[i]->Clone("h_FR_nume_sum_"+type[i])));
        h_PR_nume_sum[i]->SetDirectory(0);
        h_PR_nume_sum[i]->Add(h_PR_pT_endcap_nume[i]);
        h_PR_nume_sum[i]->Add(h_PR_pT_endcap2_nume[i]);
        h_PR_deno_sum[i] = ((TH1D*)(h_PR_pT_barrel_deno[i]->Clone("h_FR_deno_sum_"+type[i])));
        h_PR_deno_sum[i]->SetDirectory(0);
        h_PR_deno_sum[i]->Add(h_PR_pT_endcap_deno[i]);
        h_PR_deno_sum[i]->Add(h_PR_pT_endcap2_deno[i]);

        // Full rates
        h_FR_pT_barrel_nume[i] ->Divide(h_FR_pT_barrel_deno[i] );
        h_FR_pT_endcap_nume[i] ->Divide(h_FR_pT_endcap_deno[i] );
        h_FR_pT_endcap2_nume[i]->Divide(h_FR_pT_endcap2_deno[i]);
        h_FR_eta_nume[i]       ->Divide(h_FR_eta_deno[i]       );
        h_PR_pT_barrel_nume[i] ->Divide(h_PR_pT_barrel_deno[i] );
        h_PR_pT_endcap_nume[i] ->Divide(h_PR_pT_endcap_deno[i] );
        h_PR_pT_endcap2_nume[i]->Divide(h_PR_pT_endcap2_deno[i]);
        h_PR_eta_nume[i]       ->Divide(h_PR_eta_deno[i]       );
    }

    cout << "finished." << endl;
    f->Close();
    if (!f->IsOpen()) cout << "File has been closed successfully.\n" << endl;
    else cout << "FILE " << filename << " COULD NOT BE CLOSED!\n" << endl;

    // Drawing
    TLegend *legend = new TLegend(0.65, 0.8, 0.95, 0.95);
    legend->AddEntry(h_FR_pT_barrel_nume[0], "DY", "lp");
    legend->AddEntry(h_FR_pT_barrel_nume[1], "t#bar{t}", "lp");
    legend->AddEntry(h_FR_pT_barrel_nume[2], "W+Jets", "lp");
    legend->AddEntry(h_FR_pT_barrel_nume[3], "QCD", "lp");

    TCanvas *c_FR_barrel  = new TCanvas("c_FR_barrel",  "FR barrel", 800, 800);
    c_FR_barrel->cd();
    c_FR_barrel->SetGrid(1);
    c_FR_barrel->SetLogx(1);
    c_FR_barrel->SetRightMargin(0.05);
    c_FR_barrel->SetTopMargin(0.05);
    c_FR_barrel->SetBottomMargin(0.12);
    c_FR_barrel->SetLeftMargin(0.13);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_barrel_nume[0]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_barrel_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_barrel_nume[0]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_pT_barrel_nume[0]->Draw();
    h_FR_pT_barrel_nume[1]->Draw("same");
    h_FR_pT_barrel_nume[2]->Draw("same");
    h_FR_pT_barrel_nume[3]->Draw("same");
    legend->Draw();
    c_FR_barrel->Update();

    TCanvas *c_FR_endcap  = new TCanvas("c_FR_endcap",  "FR endcap", 800, 800);
    c_FR_endcap->cd();
    c_FR_endcap->SetGrid(1);
    c_FR_endcap->SetLogx(1);
    c_FR_endcap->SetRightMargin(0.05);
    c_FR_endcap->SetTopMargin(0.05);
    c_FR_endcap->SetBottomMargin(0.12);
    c_FR_endcap->SetLeftMargin(0.13);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap_nume[0]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_endcap_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_endcap_nume[0]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_pT_endcap_nume[0]->Draw();
    h_FR_pT_endcap_nume[1]->Draw("same");
    h_FR_pT_endcap_nume[2]->Draw("same");
    h_FR_pT_endcap_nume[3]->Draw("same");
    legend->Draw();
    c_FR_endcap->Update();

    TCanvas *c_FR_endcap2 = new TCanvas("c_FR_endcap2", "FR endcap2", 800, 800);
    c_FR_endcap2->cd();
    c_FR_endcap2->SetGrid(1);
    c_FR_endcap2->SetLogx(1);
    c_FR_endcap2->SetRightMargin(0.05);
    c_FR_endcap2->SetTopMargin(0.05);
    c_FR_endcap2->SetBottomMargin(0.12);
    c_FR_endcap2->SetLeftMargin(0.13);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap2_nume[0]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_endcap2_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_endcap2_nume[0]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_pT_endcap2_nume[0]->Draw();
    h_FR_pT_endcap2_nume[1]->Draw("same");
    h_FR_pT_endcap2_nume[2]->Draw("same");
    h_FR_pT_endcap2_nume[3]->Draw("same");
    legend->Draw();
    c_FR_endcap2->Update();

    TCanvas *c_FR_eta = new TCanvas("c_FR_eta", "FR eta", 800, 800);
    c_FR_eta->cd();
    c_FR_eta->SetGrid(1);
    c_FR_eta->SetRightMargin(0.05);
    c_FR_eta->SetTopMargin(0.05);
    c_FR_eta->SetBottomMargin(0.12);
    c_FR_eta->SetLeftMargin(0.13);
    h_FR_eta_nume[0]->GetXaxis()->SetTitle("#eta (#mu)");
    h_FR_eta_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_FR_eta_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_FR_eta_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_FR_eta_nume[0]->GetYaxis()->SetTitle("Fake rate");
    h_FR_eta_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_FR_eta_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_eta_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_FR_eta_nume[0]->GetXaxis()->SetNoExponent(1);
    h_FR_eta_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_eta_nume[0]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_FR_eta_nume[0]->GetYaxis()->SetRangeUser(0, 0.7);
    h_FR_eta_nume[0]->Draw();
    h_FR_eta_nume[1]->Draw("same");
    h_FR_eta_nume[2]->Draw("same");
    h_FR_eta_nume[3]->Draw("same");
    legend->Draw();
    c_FR_eta->Update();

    TCanvas *c_PR_barrel  = new TCanvas("c_PR_barrel",  "PR barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_barrel_nume[0]->Draw();
    h_PR_pT_barrel_nume[1]->Draw("same");
    h_PR_pT_barrel_nume[2]->Draw("same");
    legend->Draw();
    c_PR_barrel->Update();

    TCanvas *c_PR_endcap  = new TCanvas("c_PR_endcap",  "PR endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_endcap_nume[0]->Draw();
    h_PR_pT_endcap_nume[1]->Draw("same");
    h_PR_pT_endcap_nume[2]->Draw("same");
    legend->Draw();
    c_PR_endcap->Update();

    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "PR endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_endcap2_nume[0]->Draw();
    h_PR_pT_endcap2_nume[1]->Draw("same");
    h_PR_pT_endcap2_nume[2]->Draw("same");
    legend->Draw();
    c_PR_endcap2->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "PR eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
    h_PR_eta_nume[0]->GetXaxis()->SetTitle("#eta (#mu)");
    h_PR_eta_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_eta_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_eta_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_eta_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_eta_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_eta_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_eta_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_eta_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_eta_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_eta_nume[0]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_PR_eta_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_eta_nume[0]->Draw();
    h_PR_eta_nume[1]->Draw("same");
    h_PR_eta_nume[2]->Draw("same");
    legend->Draw();
    c_PR_eta->Update();

//    DYAnalyzer a("");
//    a.SaveCanvases("ALL", "~/Desktop/LastCopy_DYQCD", "png");

    // Writing
    TString outname = "/media/sf_DATA/FR/Electron/MC_FakeRates.root";
    cout << "Writing to " << outname << " ... ";
    TFile *f_out = new TFile(outname, "RECREATE");
    for (Int_t i=0; i<4; i++)
    {
        f_out->WriteObjectAny(h_FR_pT_barrel_nume[i],  "TH1D", "h_FR_barrel_" +type[i]);
        f_out->WriteObjectAny(h_FR_pT_endcap_nume[i],  "TH1D", "h_FR_endcap_" +type[i]);
        f_out->WriteObjectAny(h_FR_pT_endcap2_nume[i], "TH1D", "h_FR_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_FR_eta_nume[i],        "TH1D", "h_FR_eta_"    +type[i]);
        f_out->WriteObjectAny(h_PR_pT_barrel_nume[i],  "TH1D", "h_PR_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_endcap_nume[i],  "TH1D", "h_PR_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_endcap2_nume[i], "TH1D", "h_PR_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_eta_nume[i],        "TH1D", "h_PR_eta_"    +type[i]);
    }
    f_out->Close();
    cout << "Finished." << endl;
    if (!f->IsOpen()) cout << "File has been closed successfully.\n" << endl;
    else cout << "FILE COULD NOT BE CLOSED!\n" << endl;

    cout << "Average rates:" << endl;
    for (Int_t i=0; i<4; i++)
    {
        Double_t FRavg = h_FR_nume_sum[i]->Integral() / h_FR_deno_sum[i]->Integral();
        Double_t PRavg = h_PR_nume_sum[i]->Integral() / h_PR_deno_sum[i]->Integral();
        cout << type[i] << ": FR = " << FRavg << "   PR = " << PRavg << endl;
    }
    cout << endl;

} // End of E_EstFRandPR_MC()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_EstFR(Int_t type)
{
    FileMgr fm;

    TH1D *h_pT_barrel_deno, *h_pT_endcap_deno, *h_pT_endcap2_deno, *h_eta_deno,
         *h_pT_barrel_nume, *h_pT_endcap_nume, *h_pT_endcap2_nume, *h_eta_nume,
         *h_pT_barrel_ctrl, *h_pT_endcap_ctrl, *h_pT_endcap2_ctrl, *h_eta_ctrl,
         *h_FRratio_barrel,    *h_FRratio_endcap,    *h_FRratio_endcap2,    *h_FRratio_eta,
         *h_FRtemplate_barrel, *h_FRtemplate_endcap, *h_FRtemplate_endcap2,
         *h_FRdalmin_barrel,   *h_FRdalmin_endcap;

    TH1D *h_pT_barrel_MC_deno[_EndOf_Data_Special], *h_pT_endcap_MC_deno[_EndOf_Data_Special], *h_pT_endcap2_MC_deno[_EndOf_Data_Special],
         *h_pT_barrel_MC_nume[_EndOf_Data_Special], *h_pT_endcap_MC_nume[_EndOf_Data_Special], *h_pT_endcap2_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap2_MC_ctrl[_EndOf_Data_Special],
         *h_eta_MC_deno[_EndOf_Data_Special], *h_eta_MC_nume[_EndOf_Data_Special], *h_eta_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_data_deno, *h_pT_endcap_data_deno, *h_pT_endcap2_data_deno, *h_eta_data_deno,
         *h_pT_barrel_data_nume, *h_pT_endcap_data_nume, *h_pT_endcap2_data_nume, *h_eta_data_nume,
         *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl, *h_pT_endcap2_data_ctrl, *h_eta_data_ctrl;

    TH1D *h_pT_barrel_template_deno_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_50to70[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_50to70[_EndOf_Data_Special],
         *h_pT_endcap2_template_deno_50to70[_EndOf_Data_Special],
         *h_pT_endcap2_template_ctrl_50to70[_EndOf_Data_Special],
         *h_pT_endcap2_template_nume_50to70[_EndOf_Data_Special],
         *h_pT_barrel_template_deno_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_70to100[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_70to100[_EndOf_Data_Special],
         *h_pT_endcap2_template_deno_70to100[_EndOf_Data_Special],
         *h_pT_endcap2_template_ctrl_70to100[_EndOf_Data_Special],
         *h_pT_endcap2_template_nume_70to100[_EndOf_Data_Special],
         *h_pT_barrel_template_deno_100to500[_EndOf_Data_Special],
         *h_pT_barrel_template_ctrl_100to500[_EndOf_Data_Special],
         *h_pT_barrel_template_nume_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_deno_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_ctrl_100to500[_EndOf_Data_Special],
         *h_pT_endcap_template_nume_100to500[_EndOf_Data_Special],
         *h_pT_endcap2_template_deno_100to500[_EndOf_Data_Special],
         *h_pT_endcap2_template_ctrl_100to500[_EndOf_Data_Special],
         *h_pT_endcap2_template_nume_100to500[_EndOf_Data_Special];


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
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr1]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);
        file->GetObject("h_pT_endcap2_deno", h_pT_endcap2_MC_deno[pr1]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr1]);
        file->GetObject("h_eta_deno", h_eta_MC_deno[pr1]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr1]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr1]);

        removeNegativeBins(h_pT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap2_MC_deno[pr1]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr1]);
        removeNegativeBins(h_eta_MC_deno[pr1]);
        removeNegativeBins(h_eta_MC_ctrl[pr1]);
        removeNegativeBins(h_eta_MC_nume[pr1]);

        h_pT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap2_MC_deno[pr1]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr1]->SetDirectory(0);
        h_eta_MC_deno[pr1]->SetDirectory(0);
        h_eta_MC_ctrl[pr1]->SetDirectory(0);
        h_eta_MC_nume[pr1]->SetDirectory(0);

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
    h_pT_barrel_MC_ctrl[_ttbar]->Add(h_pT_barrel_MC_ctrl[_ttbar_700to1000]);
    h_pT_barrel_MC_nume[_ttbar]->Add(h_pT_barrel_MC_nume[_ttbar_700to1000]);
    h_pT_endcap_MC_deno[_ttbar]->Add(h_pT_endcap_MC_deno[_ttbar_700to1000]);
    h_pT_endcap_MC_ctrl[_ttbar]->Add(h_pT_endcap_MC_ctrl[_ttbar_700to1000]);
    h_pT_endcap_MC_nume[_ttbar]->Add(h_pT_endcap_MC_nume[_ttbar_700to1000]);
    h_pT_endcap2_MC_deno[_ttbar]->Add(h_pT_endcap2_MC_deno[_ttbar_700to1000]);
    h_pT_endcap2_MC_ctrl[_ttbar]->Add(h_pT_endcap2_MC_ctrl[_ttbar_700to1000]);
    h_pT_endcap2_MC_nume[_ttbar]->Add(h_pT_endcap2_MC_nume[_ttbar_700to1000]);
    h_eta_MC_deno[_ttbar]->Add(h_eta_MC_deno[_ttbar_700to1000]);
    h_eta_MC_ctrl[_ttbar]->Add(h_eta_MC_ctrl[_ttbar_700to1000]);
    h_eta_MC_nume[_ttbar]->Add(h_eta_MC_nume[_ttbar_700to1000]);
    h_pT_barrel_MC_deno[_ttbar]->Add(h_pT_barrel_MC_deno[_ttbar_1000toInf]);
    h_pT_barrel_MC_ctrl[_ttbar]->Add(h_pT_barrel_MC_ctrl[_ttbar_1000toInf]);
    h_pT_barrel_MC_nume[_ttbar]->Add(h_pT_barrel_MC_nume[_ttbar_1000toInf]);
    h_pT_endcap_MC_deno[_ttbar]->Add(h_pT_endcap_MC_deno[_ttbar_1000toInf]);
    h_pT_endcap_MC_ctrl[_ttbar]->Add(h_pT_endcap_MC_ctrl[_ttbar_1000toInf]);
    h_pT_endcap_MC_nume[_ttbar]->Add(h_pT_endcap_MC_nume[_ttbar_1000toInf]);
    h_pT_endcap2_MC_deno[_ttbar]->Add(h_pT_endcap2_MC_deno[_ttbar_1000toInf]);
    h_pT_endcap2_MC_ctrl[_ttbar]->Add(h_pT_endcap2_MC_ctrl[_ttbar_1000toInf]);
    h_pT_endcap2_MC_nume[_ttbar]->Add(h_pT_endcap2_MC_nume[_ttbar_1000toInf]);
    h_eta_MC_deno[_ttbar]->Add(h_eta_MC_deno[_ttbar_1000toInf]);
    h_eta_MC_ctrl[_ttbar]->Add(h_eta_MC_ctrl[_ttbar_1000toInf]);
    h_eta_MC_nume[_ttbar]->Add(h_eta_MC_nume[_ttbar_1000toInf]);
    h_pT_barrel_MC_deno[_WJets]->Add(h_pT_barrel_MC_deno[_WJets_ext2v5]);
    h_pT_barrel_MC_ctrl[_WJets]->Add(h_pT_barrel_MC_ctrl[_WJets_ext2v5]);
    h_pT_barrel_MC_nume[_WJets]->Add(h_pT_barrel_MC_nume[_WJets_ext2v5]);
    h_pT_endcap_MC_deno[_WJets]->Add(h_pT_endcap_MC_deno[_WJets_ext2v5]);
    h_pT_endcap_MC_ctrl[_WJets]->Add(h_pT_endcap_MC_ctrl[_WJets_ext2v5]);
    h_pT_endcap_MC_nume[_WJets]->Add(h_pT_endcap_MC_nume[_WJets_ext2v5]);
    h_pT_endcap2_MC_deno[_WJets]->Add(h_pT_endcap2_MC_deno[_WJets_ext2v5]);
    h_pT_endcap2_MC_ctrl[_WJets]->Add(h_pT_endcap2_MC_ctrl[_WJets_ext2v5]);
    h_pT_endcap2_MC_nume[_WJets]->Add(h_pT_endcap2_MC_nume[_WJets_ext2v5]);
    h_eta_MC_deno[_WJets]->Add(h_eta_MC_deno[_WJets_ext2v5]);
    h_eta_MC_ctrl[_WJets]->Add(h_eta_MC_ctrl[_WJets_ext2v5]);
    h_eta_MC_nume[_WJets]->Add(h_eta_MC_nume[_WJets_ext2v5]);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_endcap2_deno", h_pT_endcap2_MC_deno[pr]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr]);
        file->GetObject("h_eta_deno", h_eta_MC_deno[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap2_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_deno[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap2_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_deno[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_DY")));
            h_pT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_DY")));
            h_pT_endcap2_MC_deno[_DY_Full] = ((TH1D*)(h_pT_endcap2_MC_deno[pr]->Clone("h_pT_endcap2_MC_deno_DY")));
            h_pT_endcap2_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap2_MC_ctrl[pr]->Clone("h_pT_endcap2_MC_ctrl_DY")));
            h_pT_endcap2_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap2_MC_nume[pr]->Clone("h_pT_endcap2_MC_nume_DY")));
            h_eta_MC_deno[_DY_Full] = ((TH1D*)(h_eta_MC_deno[pr]->Clone("h_eta_MC_deno_DY")));
            h_eta_MC_ctrl[_DY_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_DY")));
            h_eta_MC_nume[_DY_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_DY")));
            h_pT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap2_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_endcap2_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap2_MC_nume[_DY_Full]->SetDirectory(0);
            h_eta_MC_deno[_DY_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_MC_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_deno[_DY_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_deno[_DY_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_endcap2_MC_deno[_DY_Full]->Add(h_pT_endcap2_MC_deno[pr]);
            h_pT_endcap2_MC_ctrl[_DY_Full]->Add(h_pT_endcap2_MC_ctrl[pr]);
            h_pT_endcap2_MC_nume[_DY_Full]->Add(h_pT_endcap2_MC_nume[pr]);
            h_eta_MC_deno[_DY_Full]->Add(h_eta_MC_deno[pr]);
            h_eta_MC_ctrl[_DY_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_DY_Full]->Add(h_eta_MC_nume[pr]);
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
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_endcap2_deno", h_pT_endcap2_MC_deno[pr]);
        file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_MC_nume[pr]);
        file->GetObject("h_eta_deno", h_eta_MC_deno[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno_50to70",   h_pT_barrel_template_deno_50to70  [pr]);
        file->GetObject("h_pT_barrel_ctrl_50to70",   h_pT_barrel_template_ctrl_50to70  [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_pT_barrel_template_nume_50to70  [pr]);
        file->GetObject("h_pT_endcap_deno_50to70",   h_pT_endcap_template_deno_50to70  [pr]);
        file->GetObject("h_pT_endcap_ctrl_50to70",   h_pT_endcap_template_ctrl_50to70  [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_pT_endcap_template_nume_50to70  [pr]);
        file->GetObject("h_pT_endcap2_deno_50to70",   h_pT_endcap2_template_deno_50to70  [pr]);
        file->GetObject("h_pT_endcap2_ctrl_50to70",   h_pT_endcap2_template_ctrl_50to70  [pr]);
        file->GetObject("h_pT_endcap2_nume_50to70",   h_pT_endcap2_template_nume_50to70  [pr]);
        file->GetObject("h_pT_barrel_deno_70to100",  h_pT_barrel_template_deno_70to100 [pr]);
        file->GetObject("h_pT_barrel_ctrl_70to100",  h_pT_barrel_template_ctrl_70to100 [pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_pT_barrel_template_nume_70to100 [pr]);
        file->GetObject("h_pT_endcap_deno_70to100",  h_pT_endcap_template_deno_70to100 [pr]);
        file->GetObject("h_pT_endcap_ctrl_70to100",  h_pT_endcap_template_ctrl_70to100 [pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_pT_endcap_template_nume_70to100 [pr]);
        file->GetObject("h_pT_endcap2_deno_70to100",  h_pT_endcap2_template_deno_70to100 [pr]);
        file->GetObject("h_pT_endcap2_ctrl_70to100",  h_pT_endcap2_template_ctrl_70to100 [pr]);
        file->GetObject("h_pT_endcap2_nume_70to100",  h_pT_endcap2_template_nume_70to100 [pr]);
        file->GetObject("h_pT_barrel_deno_100to500", h_pT_barrel_template_deno_100to500[pr]);
        file->GetObject("h_pT_barrel_ctrl_100to500", h_pT_barrel_template_ctrl_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_pT_barrel_template_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_deno_100to500", h_pT_endcap_template_deno_100to500[pr]);
        file->GetObject("h_pT_endcap_ctrl_100to500", h_pT_endcap_template_ctrl_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_pT_endcap_template_nume_100to500[pr]);
        file->GetObject("h_pT_endcap2_deno_100to500", h_pT_endcap2_template_deno_100to500[pr]);
        file->GetObject("h_pT_endcap2_ctrl_100to500", h_pT_endcap2_template_ctrl_100to500[pr]);
        file->GetObject("h_pT_endcap2_nume_100to500", h_pT_endcap2_template_nume_100to500[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap2_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap2_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap2_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_deno[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_template_deno_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_nume_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_deno_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_50to70  [pr]);
        removeNegativeBins(h_pT_endcap_template_nume_50to70  [pr]);
        removeNegativeBins(h_pT_endcap2_template_deno_50to70  [pr]);
        removeNegativeBins(h_pT_endcap2_template_ctrl_50to70  [pr]);
        removeNegativeBins(h_pT_endcap2_template_nume_50to70  [pr]);
        removeNegativeBins(h_pT_barrel_template_deno_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_nume_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_deno_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_70to100 [pr]);
        removeNegativeBins(h_pT_endcap_template_nume_70to100 [pr]);
        removeNegativeBins(h_pT_endcap2_template_deno_70to100 [pr]);
        removeNegativeBins(h_pT_endcap2_template_ctrl_70to100 [pr]);
        removeNegativeBins(h_pT_endcap2_template_nume_70to100 [pr]);
        removeNegativeBins(h_pT_barrel_template_deno_100to500[pr]);
        removeNegativeBins(h_pT_barrel_template_ctrl_100to500[pr]);
        removeNegativeBins(h_pT_barrel_template_nume_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_deno_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_ctrl_100to500[pr]);
        removeNegativeBins(h_pT_endcap_template_nume_100to500[pr]);
        removeNegativeBins(h_pT_endcap2_template_deno_100to500[pr]);
        removeNegativeBins(h_pT_endcap2_template_ctrl_100to500[pr]);
        removeNegativeBins(h_pT_endcap2_template_nume_100to500[pr]);

        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap2_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap2_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap2_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_deno[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_template_deno_50to70   [pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_50to70   [pr]->SetDirectory(0);
        h_pT_barrel_template_nume_50to70   [pr]->SetDirectory(0);
        h_pT_endcap_template_deno_50to70   [pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_50to70   [pr]->SetDirectory(0);
        h_pT_endcap_template_nume_50to70   [pr]->SetDirectory(0);
        h_pT_endcap2_template_deno_50to70  [pr]->SetDirectory(0);
        h_pT_endcap2_template_ctrl_50to70  [pr]->SetDirectory(0);
        h_pT_endcap2_template_nume_50to70  [pr]->SetDirectory(0);
        h_pT_barrel_template_deno_70to100  [pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_70to100  [pr]->SetDirectory(0);
        h_pT_barrel_template_nume_70to100  [pr]->SetDirectory(0);
        h_pT_endcap_template_deno_70to100  [pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_70to100  [pr]->SetDirectory(0);
        h_pT_endcap_template_nume_70to100  [pr]->SetDirectory(0);
        h_pT_endcap2_template_deno_70to100 [pr]->SetDirectory(0);
        h_pT_endcap2_template_ctrl_70to100 [pr]->SetDirectory(0);
        h_pT_endcap2_template_nume_70to100 [pr]->SetDirectory(0);
        h_pT_barrel_template_deno_100to500 [pr]->SetDirectory(0);
        h_pT_barrel_template_ctrl_100to500 [pr]->SetDirectory(0);
        h_pT_barrel_template_nume_100to500 [pr]->SetDirectory(0);
        h_pT_endcap_template_deno_100to500 [pr]->SetDirectory(0);
        h_pT_endcap_template_ctrl_100to500 [pr]->SetDirectory(0);
        h_pT_endcap_template_nume_100to500 [pr]->SetDirectory(0);
        h_pT_endcap2_template_deno_100to500[pr]->SetDirectory(0);
        h_pT_endcap2_template_ctrl_100to500[pr]->SetDirectory(0);
        h_pT_endcap2_template_nume_100to500[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_QCD")));
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_QCD")));
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_QCD")));
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_QCD")));
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_QCD")));
            h_pT_endcap2_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_MC_deno[pr]->Clone("h_pT_endcap2_MC_deno_QCD")));
            h_pT_endcap2_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_MC_ctrl[pr]->Clone("h_pT_endcap2_MC_ctrl_QCD")));
            h_pT_endcap2_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_MC_nume[pr]->Clone("h_pT_endcap2_MC_nume_QCD")));
            h_eta_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_eta_MC_deno[pr]->Clone("h_eta_MC_deno_QCD")));
            h_eta_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_QCD")));
            h_eta_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_QCD")));
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_50to70  [pr]->Clone("h_pT_barrel_deno_QCD_50to70"  )));
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_50to70  [pr]->Clone("h_pT_barrel_ctrl_QCD_50to70"  )));
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_50to70  [pr]->Clone("h_pT_barrel_nume_QCD_50to70"  )));
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_50to70  [pr]->Clone("h_pT_endcap_deno_QCD_50to70"  )));
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_50to70  [pr]->Clone("h_pT_endcap_ctrl_QCD_50to70"  )));
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_50to70  [pr]->Clone("h_pT_endcap_nume_QCD_50to70"  )));
            h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_deno_50to70  [pr]->Clone("h_pT_endcap2_deno_QCD_50to70"  )));
            h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_ctrl_50to70  [pr]->Clone("h_pT_endcap2_ctrl_QCD_50to70"  )));
            h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_nume_50to70  [pr]->Clone("h_pT_endcap2_nume_QCD_50to70"  )));
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_70to100 [pr]->Clone("h_pT_barrel_deno_QCD_70to100" )));
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_70to100 [pr]->Clone("h_pT_barrel_ctrl_QCD_70to100" )));
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_70to100 [pr]->Clone("h_pT_barrel_nume_QCD_70to100" )));
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_70to100 [pr]->Clone("h_pT_endcap_deno_QCD_70to100" )));
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_70to100 [pr]->Clone("h_pT_endcap_ctrl_QCD_70to100" )));
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_70to100 [pr]->Clone("h_pT_endcap_nume_QCD_70to100" )));
            h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_deno_70to100 [pr]->Clone("h_pT_endcap2_deno_QCD_70to100" )));
            h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_ctrl_70to100 [pr]->Clone("h_pT_endcap2_ctrl_QCD_70to100" )));
            h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_nume_70to100 [pr]->Clone("h_pT_endcap2_nume_QCD_70to100" )));
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_deno_100to500[pr]->Clone("h_pT_barrel_deno_QCD_100to500")));
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_ctrl_100to500[pr]->Clone("h_pT_barrel_ctrl_QCD_100to500")));
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_template_nume_100to500[pr]->Clone("h_pT_barrel_nume_QCD_100to500")));
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_deno_100to500[pr]->Clone("h_pT_endcap_deno_QCD_100to500")));
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_ctrl_100to500[pr]->Clone("h_pT_endcap_ctrl_QCD_100to500")));
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_template_nume_100to500[pr]->Clone("h_pT_endcap_nume_QCD_100to500")));
            h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_deno_100to500[pr]->Clone("h_pT_endcap2_deno_QCD_100to500")));
            h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_ctrl_100to500[pr]->Clone("h_pT_endcap2_ctrl_QCD_100to500")));
            h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap2_template_nume_100to500[pr]->Clone("h_pT_endcap2_nume_QCD_100to500")));
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_eta_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_eta_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_endcap2_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_endcap2_MC_deno[pr]);
            h_pT_endcap2_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_endcap2_MC_ctrl[pr]);
            h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_endcap2_MC_nume[pr]);
            h_eta_MC_deno[_QCDMuEnriched_Full]->Add(h_eta_MC_deno[pr]);
            h_eta_MC_ctrl[_QCDMuEnriched_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_QCDMuEnriched_Full]->Add(h_eta_MC_nume[pr]);
            h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_50to70[pr]);
            h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_50to70[pr]);
            h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_50to70[pr]);
            h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_50to70[pr]);
            h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_50to70[pr]);
            h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_50to70[pr]);
            h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_deno_50to70[pr]);
            h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_ctrl_50to70[pr]);
            h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_nume_50to70[pr]);
            h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_70to100[pr]);
            h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_70to100[pr]);
            h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_70to100[pr]);
            h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_70to100[pr]);
            h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_70to100[pr]);
            h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_70to100[pr]);
            h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_deno_70to100[pr]);
            h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_ctrl_70to100[pr]);
            h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_nume_70to100[pr]);
            h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_deno_100to500[pr]);
            h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_ctrl_100to500[pr]);
            h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Add(h_pT_barrel_template_nume_100to500[pr]);
            h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_deno_100to500[pr]);
            h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_ctrl_100to500[pr]);
            h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap_template_nume_100to500[pr]);
            h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_deno_100to500[pr]);
            h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_ctrl_100to500[pr]);
            h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Add(h_pT_endcap2_template_nume_100to500[pr]);
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
            file->GetObject("h_pT_barrel_deno", h_pT_barrel_data_deno);
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_deno", h_pT_endcap_data_deno);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            file->GetObject("h_pT_endcap2_deno", h_pT_endcap2_data_deno);
            file->GetObject("h_pT_endcap2_ctrl", h_pT_endcap2_data_ctrl);
            file->GetObject("h_pT_endcap2_nume", h_pT_endcap2_data_nume);
            file->GetObject("h_eta_deno", h_eta_data_deno);
            file->GetObject("h_eta_ctrl", h_eta_data_ctrl);
            file->GetObject("h_eta_nume", h_eta_data_nume);
            removeNegativeBins(h_pT_barrel_data_deno);
            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_deno);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_nume);
            removeNegativeBins(h_pT_endcap2_data_deno);
            removeNegativeBins(h_pT_endcap2_data_ctrl);
            removeNegativeBins(h_pT_endcap2_data_nume);
            removeNegativeBins(h_eta_data_deno);
            removeNegativeBins(h_eta_data_ctrl);
            removeNegativeBins(h_eta_data_nume);
        }
        else
        {
            file->GetObject("h_pT_barrel_deno", h_temp[0]);
            file->GetObject("h_pT_barrel_ctrl", h_temp[1]);
            file->GetObject("h_pT_barrel_nume", h_temp[2]);
            file->GetObject("h_pT_endcap_deno", h_temp[3]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[4]);
            file->GetObject("h_pT_endcap_nume", h_temp[5]);
            file->GetObject("h_pT_endcap2_deno", h_temp[6]);
            file->GetObject("h_pT_endcap2_ctrl", h_temp[7]);
            file->GetObject("h_pT_endcap2_nume", h_temp[8]);
            file->GetObject("h_eta_deno", h_temp[9]);
            file->GetObject("h_eta_ctrl", h_temp[10]);
            file->GetObject("h_eta_nume", h_temp[11]);
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
            h_pT_barrel_data_deno->Add(h_temp[0]);
            h_pT_barrel_data_ctrl->Add(h_temp[1]);
            h_pT_barrel_data_nume->Add(h_temp[2]);
            h_pT_endcap_data_deno->Add(h_temp[3]);
            h_pT_endcap_data_ctrl->Add(h_temp[4]);
            h_pT_endcap_data_nume->Add(h_temp[5]);
            h_pT_endcap2_data_deno->Add(h_temp[6]);
            h_pT_endcap2_data_ctrl->Add(h_temp[7]);
            h_pT_endcap2_data_nume->Add(h_temp[8]);
            h_eta_data_deno->Add(h_temp[9]);
            h_eta_data_ctrl->Add(h_temp[10]);
            h_eta_data_nume->Add(h_temp[11]);
        }
    }

    h_pT_barrel_data_deno->SetDirectory(0);
    h_pT_barrel_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume->SetDirectory(0);
    h_pT_endcap_data_deno->SetDirectory(0);
    h_pT_endcap_data_ctrl->SetDirectory(0);
    h_pT_endcap_data_nume->SetDirectory(0);
    h_pT_endcap2_data_deno->SetDirectory(0);
    h_pT_endcap2_data_ctrl->SetDirectory(0);
    h_pT_endcap2_data_nume->SetDirectory(0);
    h_eta_data_deno->SetDirectory(0);
    h_eta_data_ctrl->SetDirectory(0);
    h_eta_data_nume->SetDirectory(0);

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
    // Far endcap
    h_pT_endcap2_nume = ((TH1D*)(h_pT_endcap2_MC_deno[_DY_Full]->Clone("h_pT_endcap2_nume")));
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_ttbar]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_tW]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_tbarW]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_WW]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_WZ]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_ZZ]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_WJets]);
    h_pT_endcap2_nume->Add(h_pT_endcap2_MC_deno[_QCDMuEnriched_Full]);
    h_pT_endcap2_nume->Multiply(h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]);
    h_pT_endcap2_nume->Multiply(h_pT_endcap2_data_nume);
    // Eta
    h_eta_nume = ((TH1D*)(h_eta_MC_deno[_DY_Full]->Clone("h_eta_nume")));
    h_eta_nume->Add(h_eta_MC_deno[_ttbar]);
    h_eta_nume->Add(h_eta_MC_deno[_tW]);
    h_eta_nume->Add(h_eta_MC_deno[_tbarW]);
    h_eta_nume->Add(h_eta_MC_deno[_WW]);
    h_eta_nume->Add(h_eta_MC_deno[_WZ]);
    h_eta_nume->Add(h_eta_MC_deno[_ZZ]);
    h_eta_nume->Add(h_eta_MC_deno[_WJets]);
    h_eta_nume->Add(h_eta_MC_deno[_QCDMuEnriched_Full]);
    h_eta_nume->Multiply(h_eta_MC_nume[_QCDMuEnriched_Full]);
    h_eta_nume->Multiply(h_eta_data_nume);

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
    // Far endcap
    h_pT_endcap2_deno = ((TH1D*)(h_pT_endcap2_MC_nume[_DY_Full]->Clone("h_pT_endcap2_deno")));
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_ttbar]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_tW]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_tbarW]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_WW]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_WZ]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_ZZ]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_WJets]);
    h_pT_endcap2_deno->Add(h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]);
    h_pT_endcap2_deno->Multiply(h_pT_endcap2_MC_deno[_QCDMuEnriched_Full]);
    h_pT_endcap2_deno->Multiply(h_pT_endcap2_data_deno);
    // Eta
    h_eta_deno = ((TH1D*)(h_eta_MC_nume[_DY_Full]->Clone("h_eta_deno")));
    h_eta_deno->Add(h_eta_MC_nume[_ttbar]);
    h_eta_deno->Add(h_eta_MC_nume[_tW]);
    h_eta_deno->Add(h_eta_MC_nume[_tbarW]);
    h_eta_deno->Add(h_eta_MC_nume[_WW]);
    h_eta_deno->Add(h_eta_MC_nume[_WZ]);
    h_eta_deno->Add(h_eta_MC_nume[_ZZ]);
    h_eta_deno->Add(h_eta_MC_nume[_WJets]);
    h_eta_deno->Add(h_eta_MC_nume[_QCDMuEnriched_Full]);
    h_eta_deno->Multiply(h_eta_MC_deno[_QCDMuEnriched_Full]);
    h_eta_deno->Multiply(h_eta_data_deno);

    // ######## FR ######## //
    // Barrel
    h_FRratio_barrel = ((TH1D*)(h_pT_barrel_nume->Clone("h_FRratio_barrel")));
    h_FRratio_barrel->Divide(h_pT_barrel_deno);
    h_FRratio_barrel->SetDirectory(0);
    // Endcap
    h_FRratio_endcap = ((TH1D*)(h_pT_endcap_nume->Clone("h_FRratio_endcap")));
    h_FRratio_endcap->Divide(h_pT_endcap_deno);
    h_FRratio_endcap->SetDirectory(0);
    // Far endcap
    h_FRratio_endcap2 = ((TH1D*)(h_pT_endcap2_nume->Clone("h_FRratio_endcap2")));
    h_FRratio_endcap2->Divide(h_pT_endcap2_deno);
    h_FRratio_endcap2->SetDirectory(0);
    // Eta
    h_FRratio_eta = ((TH1D*)(h_eta_nume->Clone("h_FRratio_eta")));
    h_FRratio_eta->Divide(h_eta_deno);
    h_FRratio_eta->SetDirectory(0);



//--------------------------------- FR by template --------------------------------------  

    // Barrel
    // Normal
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.0265e+06/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.5464e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.6875e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(4.3278e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(5.6599e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(8.9442e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets>7
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(5.5608e+05/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.0083e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.5093e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(2.8957e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(3.4871e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(6.1411e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets<=7
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(2.3906e+05/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(3.3006e+06/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(5.1176e+04/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(8.0012e+05/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(9.1969e+03/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(1.3859e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // Alternative binning
    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(8.1604e+05/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.3397e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.0720e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(3.6988e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(4.4790e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(7.5344e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());


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
    // Normal
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(6.2621e+05/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(5.0765e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.4418e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.2941e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(2.5725e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(2.2921e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets>7
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(3.5338e+05/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(3.2782e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(8.7577e+04/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(8.5841e+05/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(1.5271e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(1.5379e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets<=7
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.5174e+05/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.1162e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.5009e+04/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(2.4386e+05/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(3.6723e+03/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(3.7746e+04/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // Alternative binning
    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(5.1724e+05/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(4.4065e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.1506e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.1024e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(1.9495e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(1.9229e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());

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

    // Far endcap
    // Normal
//    h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(5.1056e+05/h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(2.5062e+06/h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.2484e+05/h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(5.9138e+05/h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Scale(2.3828e+04/h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Scale(8.7950e+04/h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njet>7
//    h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(3.0241e+05/h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.6197e+06/h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(7.4037e+04/h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(3.8883e+05/h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Scale(1.4858e+04/h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Scale(6.0151e+04/h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njet<=7
//    h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.2311e+05/h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(5.8173e+05/h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.3830e+04/h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.2590e+05/h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Scale(4.0893e+03/h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Scale(1.5250e+04/h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    // Alternative binning
    h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(4.3393e+05/h_pT_endcap2_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(2.2046e+06/h_pT_endcap2_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.0097e+05/h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(5.1535e+05/h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Scale(1.9238e+04/h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Scale(7.5564e+04/h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]->Integral());

    TH1D *h_FRtemplate_endcap2_deno = ((TH1D*)(h_pT_endcap2_template_deno_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_endcap2_deno")));
    h_FRtemplate_endcap2_deno->Add(h_pT_endcap2_template_deno_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_endcap2_deno->Add(h_pT_endcap2_template_deno_100to500[_QCDMuEnriched_Full]);

    h_FRtemplate_endcap2 = ((TH1D*)(h_pT_endcap2_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRtemplate_endcap2")));
    h_FRtemplate_endcap2->Add(h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRtemplate_endcap2->Add(h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]);

    h_FRtemplate_endcap2->Divide(h_FRtemplate_endcap2_deno);
    h_FRtemplate_endcap2->SetDirectory(0);

    // For systematic errors
    TH1D *h_FRtemplate_endcap2_up   = ((TH1D*)(h_FRtemplate_endcap2->Clone("h_FRtemplate_endcap2_up")));
    TH1D *h_FRtemplate_endcap2_down = ((TH1D*)(h_FRtemplate_endcap2->Clone("h_FRtemplate_endcap2_down")));
    for (Int_t i=1; i<h_FRtemplate_endcap2->GetSize()-1; i++)
    {
        h_FRtemplate_endcap2_up  ->SetBinContent(i, h_FRtemplate_endcap2->GetBinContent(i)+h_FRtemplate_endcap2->GetBinError(i));
        h_FRtemplate_endcap2_down->SetBinContent(i, h_FRtemplate_endcap2->GetBinContent(i)-h_FRtemplate_endcap2->GetBinError(i));
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

    // Far endcap
    TH1D * h_FRmixed_endcap2 = ((TH1D*)(h_pT_endcap2_data_nume->Clone("h_FRmixed_endcap2"))); // so far it is just a numerator of numerator
    h_FRmixed_endcap2->SetDirectory(0);
    h_FRmixed_endcap2->Multiply(h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]);

    TH1D * h_pT_endcap2_mixed_nume_2 = ((TH1D*)(h_pT_endcap2_MC_nume[_DY_Full]->Clone("h_pT_endcap2_deno"))); // denominator of numerator
    h_pT_endcap2_mixed_nume_2->SetDirectory(0);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_ttbar]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_tW]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_tbarW]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_WW]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_WZ]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_ZZ]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_WJets]);
    h_pT_endcap2_mixed_nume_2->Add(h_pT_endcap2_MC_nume[_QCDMuEnriched_Full]);

    h_FRmixed_endcap2->Divide(h_pT_endcap2_mixed_nume_2);

    // ----- Denominator ----- //
    TH1D* h_pT_barrel_mixed_deno = ((TH1D*)(h_FRtemplate_barrel_deno->Clone("h_pT_barrel_mixed_deno")));
    TH1D* h_pT_endcap_mixed_deno = ((TH1D*)(h_FRtemplate_endcap_deno->Clone("h_pT_endcap_mixed_deno")));
    TH1D* h_pT_endcap2_mixed_deno = ((TH1D*)(h_FRtemplate_endcap2_deno->Clone("h_pT_endcap2_mixed_deno")));
    h_pT_barrel_mixed_deno->SetDirectory(0);
    h_pT_endcap_mixed_deno->SetDirectory(0);
    h_pT_endcap2_mixed_deno->SetDirectory(0);

    // ----- Fake Rate ----- //
    h_FRmixed_barrel->Divide(h_pT_barrel_mixed_deno);
    h_FRmixed_endcap->Divide(h_pT_endcap_mixed_deno);
    h_FRmixed_endcap2->Divide(h_pT_endcap2_mixed_deno);

//-------------------------- FR from signal+control fit ------------------------------
    // Signal and control regions -- from template fit
    // Numerator = Signal, Denominator = Signal + Control
    Double_t err_signal_barrel[18],  err_signal_endcap[9],  err_signal_endcap2[9],
            err_control_barrel[18], err_control_endcap[9], err_control_endcap2[9];
    Double_t val_signal_barrel[18],  val_signal_endcap[9],  val_signal_endcap2[9],
            val_control_barrel[18], val_control_endcap[9], val_control_endcap2[9];

    // ----- Numerator ----- //
    TH1D *h_FRsigCtrl_template_barrel = ((TH1D*)(h_pT_barrel_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRsigCtrl_template_barrel")));
    h_FRsigCtrl_template_barrel->Add(h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRsigCtrl_template_barrel->Add(h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]);
    TH1D *h_FRsigCtrl_template_endcap = ((TH1D*)(h_pT_endcap_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRsigCtrl_template_endcap")));
    h_FRsigCtrl_template_endcap->Add(h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRsigCtrl_template_endcap->Add(h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]);
    TH1D *h_FRsigCtrl_template_endcap2 = ((TH1D*)(h_pT_endcap2_template_nume_50to70[_QCDMuEnriched_Full]->Clone("h_FRsigCtrl_template_endcap2")));
    h_FRsigCtrl_template_endcap2->Add(h_pT_endcap2_template_nume_70to100 [_QCDMuEnriched_Full]);
    h_FRsigCtrl_template_endcap2->Add(h_pT_endcap2_template_nume_100to500[_QCDMuEnriched_Full]);

    // ---- Denominator ---- //
    // Normal
//    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.4666e+07/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(4.1323e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(8.5056e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(4.5911e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.1797e+06/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(2.0436e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(2.0989e+06/h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(4.8647e+05/h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(7.1678e+04/h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets>7
//    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(9.6273e+06/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(2.7780e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(5.9027e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(2.9891e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(7.8610e+05/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(1.3830e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.3618e+06/h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(3.2843e+05/h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(4.9485e+04/h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    // FR njets<=7
//    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(3.1309e+06/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(7.6302e+05/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(1.3297e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.0007e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(2.2238e+05/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(3.4064e+04/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(4.8883e+05/h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.0038e+05/h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(1.2567e+04/h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    // Alternative binning
    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.2773e+07/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(3.5414e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(7.2390e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(3.9930e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.0106e+06/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(1.7285e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.8537e+06/h_pT_endcap2_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(4.2948e+05/h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(6.1643e+04/h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());


    TH1D *h_DenoSigCtrl_template_barrel = ((TH1D*)(h_pT_barrel_template_ctrl_50to70[_QCDMuEnriched_Full]->Clone("h_DenoSigCtrl_template_barrel")));
    h_DenoSigCtrl_template_barrel->Add(h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]);
    h_DenoSigCtrl_template_barrel->Add(h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]);
    TH1D *h_DenoSigCtrl_template_endcap = ((TH1D*)(h_pT_endcap_template_ctrl_50to70[_QCDMuEnriched_Full]->Clone("h_DenoSigCtrl_template_endcap")));
    h_DenoSigCtrl_template_endcap->Add(h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]);
    h_DenoSigCtrl_template_endcap->Add(h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]);
    TH1D *h_DenoSigCtrl_template_endcap2 = ((TH1D*)(h_pT_endcap2_template_ctrl_50to70[_QCDMuEnriched_Full]->Clone("h_DenoSigCtrl_template_endcap2")));
    h_DenoSigCtrl_template_endcap2->Add(h_pT_endcap2_template_ctrl_70to100 [_QCDMuEnriched_Full]);
    h_DenoSigCtrl_template_endcap2->Add(h_pT_endcap2_template_ctrl_100to500[_QCDMuEnriched_Full]);

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

            err_signal_endcap2[i_bin] = sqrt(h_FRsigCtrl_template_endcap2->GetBinContent(i_bin+1));
            err_control_endcap2[i_bin] = sqrt(h_DenoSigCtrl_template_endcap2->GetBinContent(i_bin+1));
            val_signal_endcap2[i_bin] = h_FRsigCtrl_template_endcap2->GetBinContent(i_bin+1);
            val_control_endcap2[i_bin] = h_DenoSigCtrl_template_endcap2->GetBinContent(i_bin+1);
        }
    }

    h_DenoSigCtrl_template_barrel->Add(h_FRsigCtrl_template_barrel); // Deno = Signal+Control
    h_DenoSigCtrl_template_endcap->Add(h_FRsigCtrl_template_endcap);
    h_DenoSigCtrl_template_endcap2->Add(h_FRsigCtrl_template_endcap2);

    // ----- Fake rate ----- //
    h_FRsigCtrl_template_barrel->Divide(h_DenoSigCtrl_template_barrel);
    h_FRsigCtrl_template_endcap->Divide(h_DenoSigCtrl_template_endcap);
    h_FRsigCtrl_template_endcap2->Divide(h_DenoSigCtrl_template_endcap2);

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
            Double_t err_endcap, err_endcap2;
            err_endcap = (val_signal_endcap[i_bin] * val_control_endcap[i_bin]) /
                        ((val_signal_endcap[i_bin] + val_control_endcap[i_bin]) * (val_signal_endcap[i_bin] + val_control_endcap[i_bin])) *
                         sqrt(1/val_signal_endcap[i_bin] + 1/val_control_endcap[i_bin]);
            err_endcap2 = (val_signal_endcap2[i_bin] * val_control_endcap2[i_bin]) /
                         ((val_signal_endcap2[i_bin] + val_control_endcap2[i_bin]) * (val_signal_endcap2[i_bin] + val_control_endcap2[i_bin])) *
                         sqrt(1/val_signal_endcap2[i_bin] + 1/val_control_endcap2[i_bin]);
            h_FRsigCtrl_template_endcap->SetBinError(i_bin+1, err_endcap);
            h_FRsigCtrl_template_endcap2->SetBinError(i_bin+1, err_endcap2);
        }
    }

    TH1D* h_FRsigCtrl_template_barrel_plus   = ((TH1D*)(h_FRsigCtrl_template_barrel ->Clone("h_FRsigCtrl_template_plus_barrel")));
    TH1D* h_FRsigCtrl_template_endcap_plus   = ((TH1D*)(h_FRsigCtrl_template_endcap ->Clone("h_FRsigCtrl_template_plus_endcap")));
    TH1D* h_FRsigCtrl_template_endcap2_plus  = ((TH1D*)(h_FRsigCtrl_template_endcap2->Clone("h_FRsigCtrl_template_plus_endcap2")));
    TH1D* h_FRsigCtrl_template_barrel_minus  = ((TH1D*)(h_FRsigCtrl_template_barrel ->Clone("h_FRsigCtrl_template_minus_barrel")));
    TH1D* h_FRsigCtrl_template_endcap_minus  = ((TH1D*)(h_FRsigCtrl_template_endcap ->Clone("h_FRsigCtrl_template_minus_endcap")));
    TH1D* h_FRsigCtrl_template_endcap2_minus = ((TH1D*)(h_FRsigCtrl_template_endcap2->Clone("h_FRsigCtrl_template_minus_endcap2")));
    h_FRsigCtrl_template_barrel_plus  ->SetDirectory(0);
    h_FRsigCtrl_template_endcap_plus  ->SetDirectory(0);
    h_FRsigCtrl_template_endcap2_plus ->SetDirectory(0);
    h_FRsigCtrl_template_barrel_minus ->SetDirectory(0);
    h_FRsigCtrl_template_endcap_minus ->SetDirectory(0);
    h_FRsigCtrl_template_endcap2_minus->SetDirectory(0);
    for (Int_t i_bin=1; i_bin<=18; i_bin++)
    {
        if (h_FRsigCtrl_template_barrel_plus ->GetBinContent(i_bin) + h_FRsigCtrl_template_barrel_plus ->GetBinError(i_bin) < 1)
            h_FRsigCtrl_template_barrel_plus ->SetBinContent(i_bin, h_FRsigCtrl_template_barrel_plus ->GetBinContent(i_bin) + h_FRsigCtrl_template_barrel_plus ->GetBinError(i_bin));
        else h_FRsigCtrl_template_barrel_plus ->SetBinContent(i_bin, 1);

        if (h_FRsigCtrl_template_barrel_minus ->GetBinContent(i_bin) - h_FRsigCtrl_template_barrel_minus ->GetBinError(i_bin) > 0)
            h_FRsigCtrl_template_barrel_minus ->SetBinContent(i_bin, h_FRsigCtrl_template_barrel_minus ->GetBinContent(i_bin) - h_FRsigCtrl_template_barrel_minus ->GetBinError(i_bin));
        else h_FRsigCtrl_template_barrel_minus ->SetBinContent(i_bin, 0);

        if (i_bin <= 9)
        {
            if (h_FRsigCtrl_template_endcap_plus ->GetBinContent(i_bin) + h_FRsigCtrl_template_endcap_plus ->GetBinError(i_bin) < 1)
                h_FRsigCtrl_template_endcap_plus ->SetBinContent(i_bin, h_FRsigCtrl_template_endcap_plus ->GetBinContent(i_bin) + h_FRsigCtrl_template_endcap_plus ->GetBinError(i_bin));
            else h_FRsigCtrl_template_endcap_plus ->SetBinContent(i_bin, 1);
            if (h_FRsigCtrl_template_endcap2_plus->GetBinContent(i_bin) + h_FRsigCtrl_template_endcap2_plus->GetBinError(i_bin) < 1)
                h_FRsigCtrl_template_endcap2_plus->SetBinContent(i_bin, h_FRsigCtrl_template_endcap2_plus->GetBinContent(i_bin) + h_FRsigCtrl_template_endcap2_plus->GetBinError(i_bin));
            else h_FRsigCtrl_template_endcap2_plus->SetBinContent(i_bin, 1);

            if (h_FRsigCtrl_template_endcap_minus ->GetBinContent(i_bin) - h_FRsigCtrl_template_endcap_minus ->GetBinError(i_bin) > 0)
                h_FRsigCtrl_template_endcap_minus ->SetBinContent(i_bin, h_FRsigCtrl_template_endcap_minus ->GetBinContent(i_bin) - h_FRsigCtrl_template_endcap_minus ->GetBinError(i_bin));
            else h_FRsigCtrl_template_endcap_minus ->SetBinContent(i_bin, 0);
            if (h_FRsigCtrl_template_endcap2_minus->GetBinContent(i_bin) - h_FRsigCtrl_template_endcap2_minus->GetBinError(i_bin) > 0)
                h_FRsigCtrl_template_endcap2_minus->SetBinContent(i_bin, h_FRsigCtrl_template_endcap2_minus->GetBinContent(i_bin) - h_FRsigCtrl_template_endcap2_minus->GetBinError(i_bin));
            else h_FRsigCtrl_template_endcap2_minus->SetBinContent(i_bin,0);
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
    h_FRratio_endcap2->Write();
    h_FRratio_eta->Write();
    h_FRtemplate_barrel->Write();
    h_FRtemplate_endcap->Write();
    h_FRtemplate_endcap2->Write();
    h_FRtemplate_barrel_up->Write();
    h_FRtemplate_endcap_up->Write();
    h_FRtemplate_endcap2_up->Write();
    h_FRtemplate_barrel_down->Write();
    h_FRtemplate_endcap_down->Write();
    h_FRtemplate_endcap2_down->Write();
    h_FRmixed_barrel->Write();
    h_FRmixed_endcap->Write();
    h_FRmixed_endcap2->Write();
    h_FRsigCtrl_template_barrel->Write();
    h_FRsigCtrl_template_endcap->Write();
    h_FRsigCtrl_template_endcap2->Write();
    h_FRsigCtrl_template_barrel_plus ->Write();
    h_FRsigCtrl_template_endcap_plus ->Write();
    h_FRsigCtrl_template_endcap2_plus ->Write();
    h_FRsigCtrl_template_barrel_minus ->Write();
    h_FRsigCtrl_template_endcap_minus ->Write();
    h_FRsigCtrl_template_endcap2_minus->Write();
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
//    h_FRmixed_barrel->Draw("same");
    h_FRtemplate_barrel->SetMarkerStyle(33);
    h_FRtemplate_barrel->SetMarkerColor(kGreen+2);
    h_FRtemplate_barrel->SetLineColor(kGreen+2);
    h_FRtemplate_barrel->SetStats(kFALSE);
//    h_FRtemplate_barrel->Draw("same");
    h_FRsigCtrl_template_barrel->SetMarkerStyle(kFullDotLarge);
    h_FRsigCtrl_template_barrel->SetMarkerColor(kBlack);
    h_FRsigCtrl_template_barrel->SetLineColor(kBlack);
    h_FRsigCtrl_template_barrel->SetStats(kFALSE);
    h_FRsigCtrl_template_barrel->Draw("same");

    TLegend *legend = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend->AddEntry(h_FRratio_barrel, "Ratio", "LP");
//    legend->AddEntry(h_FRmixed_barrel, "Template (deno)", "LP");
//    legend->AddEntry(h_FRtemplate_barrel, "Template (nume, deno)", "LP");
//    legend->AddEntry(h_FRsigCtrl_template_barrel, "Template (signal, non-signal)", "LP");
    legend->AddEntry(h_FRratio_barrel, "MC ratio", "LP");
    legend->AddEntry(h_FRsigCtrl_template_barrel, "Template fit", "LP");
    legend->Draw();
    TText *textb = new TText (0.45, 0.6, "Barrel");
    textb->SetTextAlign(11);
    textb->SetTextSize(0.05);
    textb->SetNDC(true);
    textb->Draw();
    c_FR_barrel->Update();

    TCanvas *c_FR_endcap = new TCanvas("c_FR_endcap", "c_FR_endcap", 800, 800);
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
//    h_FRmixed_endcap->Draw("same");
    h_FRtemplate_endcap->SetMarkerStyle(33);
    h_FRtemplate_endcap->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap->SetLineColor(kGreen+2);
    h_FRtemplate_endcap->SetStats(kFALSE);
//    h_FRtemplate_endcap->Draw("same");
    h_FRsigCtrl_template_endcap->SetMarkerStyle(kFullDotLarge);
    h_FRsigCtrl_template_endcap->SetMarkerColor(kBlack);
    h_FRsigCtrl_template_endcap->SetLineColor(kBlack);
    h_FRsigCtrl_template_endcap->SetStats(kFALSE);
    h_FRsigCtrl_template_endcap->Draw("same");
    legend->Draw();
    TLatex *texte = new TLatex (0.45, 0.6, "Endcap |#eta < 1.8|");
    texte->SetTextAlign(11);
    texte->SetTextSize(0.05);
    texte->SetNDC(true);
    texte->Draw();
    c_FR_endcap->Update();

    TCanvas *c_FR_endcap2 = new TCanvas("c_FR_endcap2", "c_FR_endcap2", 800, 800);
    c_FR_endcap2->cd();
    c_FR_endcap2->SetGrid(1);
    c_FR_endcap2->SetLogx(1);
    c_FR_endcap2->SetRightMargin(0.05);
    c_FR_endcap2->SetTopMargin(0.05);
    c_FR_endcap2->SetBottomMargin(0.12);
    c_FR_endcap2->SetLeftMargin(0.13);
    h_FRratio_endcap2->SetMarkerStyle(kFullSquare);
    h_FRratio_endcap2->SetMarkerColor(kRed);
    h_FRratio_endcap2->SetLineColor(kRed);
    h_FRratio_endcap2->SetStats(kFALSE);
    h_FRratio_endcap2->SetTitle("");
    h_FRratio_endcap2->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FRratio_endcap2->GetXaxis()->SetTitleOffset(1);
    h_FRratio_endcap2->GetXaxis()->SetTitleSize(0.05);
    h_FRratio_endcap2->GetXaxis()->SetLabelSize(0.04);
    h_FRratio_endcap2->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_endcap2->GetYaxis()->SetTitleSize(0.05);
    h_FRratio_endcap2->GetYaxis()->SetTitleOffset(1.25);
    h_FRratio_endcap2->GetYaxis()->SetLabelSize(0.04);
    h_FRratio_endcap2->GetXaxis()->SetNoExponent(1);
    h_FRratio_endcap2->GetXaxis()->SetMoreLogLabels(1);
    h_FRratio_endcap2->GetXaxis()->SetRangeUser(52, 1000);
    h_FRratio_endcap2->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_endcap2->Draw();
    h_FRmixed_endcap2->SetMarkerStyle(22);
    h_FRmixed_endcap2->SetMarkerColor(kBlue);
    h_FRmixed_endcap2->SetLineColor(kBlue);
    h_FRmixed_endcap2->SetStats(kFALSE);
//    h_FRmixed_endcap2->Draw("same");
    h_FRtemplate_endcap2->SetMarkerStyle(33);
    h_FRtemplate_endcap2->SetMarkerColor(kGreen+2);
    h_FRtemplate_endcap2->SetLineColor(kGreen+2);
    h_FRtemplate_endcap2->SetStats(kFALSE);
//    h_FRtemplate_endcap2->Draw("same");
    h_FRsigCtrl_template_endcap2->SetMarkerStyle(kFullDotLarge);
    h_FRsigCtrl_template_endcap2->SetMarkerColor(kBlack);
    h_FRsigCtrl_template_endcap2->SetLineColor(kBlack);
    h_FRsigCtrl_template_endcap2->SetStats(kFALSE);
    h_FRsigCtrl_template_endcap2->Draw("same");
    legend->Draw();
    TLatex *texte2 = new TLatex (0.45, 0.6, "Endcap |1.8 #leq #eta < 2.4|");
    texte2->SetTextAlign(11);
    texte2->SetTextSize(0.05);
    texte2->SetNDC(true);
    texte2->Draw();
    c_FR_endcap2->Update();

    TCanvas *c_FR_eta = new TCanvas("c_FR_eta", "c_FR_eta", 800, 800);
    c_FR_eta->cd();
    c_FR_eta->SetGrid(1);
    c_FR_eta->SetRightMargin(0.05);
    c_FR_eta->SetTopMargin(0.05);
    c_FR_eta->SetBottomMargin(0.12);
    c_FR_eta->SetLeftMargin(0.13);
    h_FRratio_eta->SetMarkerStyle(kFullSquare);
    h_FRratio_eta->SetMarkerColor(kRed);
    h_FRratio_eta->SetLineColor(kRed);
    h_FRratio_eta->SetStats(kFALSE);
    h_FRratio_eta->SetTitle("");
    h_FRratio_eta->GetXaxis()->SetTitle("#eta_{#mu} [GeV/c]");
    h_FRratio_eta->GetXaxis()->SetTitleOffset(1);
    h_FRratio_eta->GetXaxis()->SetTitleSize(0.05);
    h_FRratio_eta->GetXaxis()->SetLabelSize(0.04);
    h_FRratio_eta->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_eta->GetYaxis()->SetTitleSize(0.05);
    h_FRratio_eta->GetYaxis()->SetTitleOffset(1.25);
    h_FRratio_eta->GetYaxis()->SetLabelSize(0.04);
    h_FRratio_eta->GetXaxis()->SetRangeUser(-2.42, 2.42);
    h_FRratio_eta->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_eta->Draw();
    TLegend *legend_eta = new TLegend(0.13, 0.87, 0.45, 0.95);
    legend_eta->AddEntry(h_FRratio_eta, "MC ratio", "LP");
    legend_eta->Draw();
    c_FR_eta->Update();

    TCanvas *c_FR_allin1 = new TCanvas("c_FR_allin1", "c_FR_allin1", 800, 800);
    c_FR_allin1->cd();
    c_FR_allin1->SetGrid(1);
    c_FR_allin1->SetRightMargin(0.05);
    c_FR_allin1->SetTopMargin(0.05);
    c_FR_allin1->SetBottomMargin(0.12);
    c_FR_allin1->SetLeftMargin(0.13);
    h_FRsigCtrl_template_barrel->SetTitle("");
    TH1D *h_FRsigCtrl_template_endcap_draw = ((TH1D*)(h_FRsigCtrl_template_endcap->Clone("h_FRsigCtrl_template_endcap_draw")));
    TH1D *h_FRsigCtrl_template_endcap2_draw = ((TH1D*)(h_FRsigCtrl_template_endcap2->Clone("h_FRsigCtrl_template_endcap2_draw")));
    TH1D *h_FRratio_barrel_draw = ((TH1D*)(h_FRratio_barrel->Clone("h_FRratio_barrel_draw")));
    TH1D *h_FRratio_endcap_draw = ((TH1D*)(h_FRratio_endcap->Clone("h_FRratio_endcap_draw")));
    TH1D *h_FRratio_endcap2_draw = ((TH1D*)(h_FRratio_endcap2->Clone("h_FRratio_endcap2_draw")));
    h_FRsigCtrl_template_endcap2_draw->SetTitle("");
    h_FRsigCtrl_template_endcap2_draw->SetMarkerStyle(kFullSquare);
    h_FRsigCtrl_template_endcap2_draw->SetMarkerColor(kBlue);
    h_FRsigCtrl_template_endcap2_draw->SetLineColor(kBlue);
    h_FRsigCtrl_template_endcap_draw->SetTitle("");
    h_FRsigCtrl_template_endcap_draw->SetMarkerStyle(33);
    h_FRsigCtrl_template_endcap_draw->SetMarkerSize(1.5);
    h_FRsigCtrl_template_endcap_draw->SetMarkerColor(kOrange+8);
    h_FRsigCtrl_template_endcap_draw->SetLineColor(kOrange+8);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetNoExponent(1);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
    h_FRsigCtrl_template_endcap_draw->GetXaxis()->SetRangeUser(52, 1000);
    h_FRsigCtrl_template_endcap_draw->GetYaxis()->SetTitle("Fake rate");
    h_FRsigCtrl_template_endcap_draw->GetYaxis()->SetTitleSize(0.045);
    h_FRsigCtrl_template_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
    h_FRsigCtrl_template_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_FRsigCtrl_template_endcap_draw->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_endcap_draw->SetMarkerColor(kMagenta+3);
    h_FRratio_endcap_draw->SetLineColor(kMagenta+3);
    h_FRratio_endcap2_draw->GetYaxis()->SetTitle("Fake rate");
    h_FRratio_endcap2_draw->GetYaxis()->SetTitleSize(0.045);
    h_FRratio_endcap2_draw->GetYaxis()->SetTitleSize(0.045);
    h_FRratio_endcap2_draw->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_endcap2_draw->SetMarkerStyle(33);
    h_FRratio_endcap2_draw->SetMarkerSize(1.5);
    h_FRratio_endcap2_draw->SetMarkerColor(kOrange+1);
    h_FRratio_endcap2_draw->SetLineColor(kOrange+1);
    h_FRratio_barrel_draw->SetMarkerStyle(kFullDotLarge);
//    h_FRratio_endcap_draw->Draw();
//    h_FRratio_endcap2_draw->Draw("same");
//    h_FRratio_barrel_draw->Draw("same");
    h_FRsigCtrl_template_endcap_draw->Draw();
    h_FRsigCtrl_template_barrel->Draw("same");
    h_FRsigCtrl_template_endcap2_draw->Draw("same");
    TLegend *legend1 = new TLegend(0.13, 0.81, 0.5, 0.95);
//    legend1->AddEntry(h_FRratio_barrel_draw, "Ratio (barrel)", "LP");
//    legend1->AddEntry(h_FRratio_endcap_draw, "Ratio (endcap |#eta<1.8|)", "LP");
//    legend1->AddEntry(h_FRratio_endcap2_draw, "Ratio (endcap |1.8#leq#eta<2.4|)", "LP");
//    legend1->AddEntry(h_FRsigCtrl_template_barrel, "Template (barrel)", "LP");
//    legend1->AddEntry(h_FRsigCtrl_template_endcap_draw, "Template (endcap |#eta<1.8|)", "LP");
//    legend1->AddEntry(h_FRsigCtrl_template_endcap2_draw, "Template (endcap |1.8#leq#eta<2.4|)", "LP");
    legend1->AddEntry(h_FRsigCtrl_template_barrel, "|#eta|<1.2", "LP");
    legend1->AddEntry(h_FRsigCtrl_template_endcap_draw, "1.2#leq|#eta|<1.8", "LP");
    legend1->AddEntry(h_FRsigCtrl_template_endcap2_draw, "1.8#leq|#eta|<2.4", "LP");
    legend1->Draw();
    c_FR_allin1->SetLogx();
    c_FR_allin1->Update();

} // End of Mu_EstFR()


void E_EstPR()
{
    TString inName = "/media/sf_DATA/FR/Electron/PR_Hist_E.root";
    TFile *f = new TFile(inName, "READ");

    TH1D *h_pT_barrel_pass[4],
         *h_pT_endcap_pass[4],
         *h_pT_endcap2_pass[4],
         *h_eta_pass[4],
         *h_pT_barrel_fail[4],
         *h_pT_endcap_fail[4],
         *h_pT_endcap2_fail[4],
         *h_eta_fail[4];
    TString type[4] = {"data", "DY", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_endcap2_pass_"+type[i], h_pT_endcap2_pass[i]);
        f->GetObject("h_eta_pass_"+type[i], h_eta_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        f->GetObject("h_pT_endcap2_fail_"+type[i], h_pT_endcap2_fail[i]);
        f->GetObject("h_eta_fail_"+type[i], h_eta_fail[i]);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_endcap2_pass[i]->SetDirectory(0);
        h_eta_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
        h_pT_endcap2_fail[i]->SetDirectory(0);
        h_eta_fail[i]->SetDirectory(0);
    }

    TH1D *h_eff_data_barrel = ((TH1D*)(h_pT_barrel_pass[0]->Clone("h_PR_subtract_barrel")));
    TH1D *h_eff_data_endcap = ((TH1D*)(h_pT_endcap_pass[0]->Clone("h_PR_subtract_endcap")));
    TH1D *h_eff_data_endcap2 = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_subtract_endcap2")));
    TH1D *h_eff_data_eta = ((TH1D*)(h_eta_pass[0]->Clone("h_PR_subtract_eta")));
    TH1D *h_eff_MC_barrel = ((TH1D*)(h_pT_barrel_pass[1]->Clone("h_PR_MC_barrel")));
    TH1D *h_eff_MC_endcap = ((TH1D*)(h_pT_endcap_pass[1]->Clone("h_PR_MC_endcap")));
    TH1D *h_eff_MC_endcap2 = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_MC_endcap2")));
    TH1D *h_eff_MC_eta = ((TH1D*)(h_eta_pass[1]->Clone("h_PR_MC_eta")));
    TH1D *h_eff_bkgr_barrel = ((TH1D*)(h_pT_barrel_pass[2]->Clone("h_PR_bkgr_barrel")));
    TH1D *h_eff_bkgr_endcap = ((TH1D*)(h_pT_endcap_pass[2]->Clone("h_PR_bkgr_endcap")));
    TH1D *h_eff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_pass[2]->Clone("h_PR_bkgr_endcap2")));
    TH1D *h_eff_bkgr_eta = ((TH1D*)(h_eta_pass[2]->Clone("h_PR_bkgr_eta")));
    TH1D *h_eff_bkgf_barrel = ((TH1D*)(h_pT_barrel_pass[3]->Clone("h_PR_bkgf_barrel")));
    TH1D *h_eff_bkgf_endcap = ((TH1D*)(h_pT_endcap_pass[3]->Clone("h_PR_bkgf_endcap")));
    TH1D *h_eff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_pass[3]->Clone("h_PR_bkgf_endcap2")));
    TH1D *h_eff_bkgf_eta = ((TH1D*)(h_eta_pass[3]->Clone("h_PR_bkgf_eta")));
    h_eff_data_barrel->Add(h_eff_bkgf_barrel, -1);
    h_eff_data_endcap->Add(h_eff_bkgf_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgf_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgf_eta, -1);
    h_eff_data_barrel->Add(h_eff_bkgr_barrel, -1);
    h_eff_data_endcap->Add(h_eff_bkgr_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgr_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgr_eta, -1);
    removeNegativeBins(h_eff_data_barrel);
    removeNegativeBins(h_eff_data_endcap);
    removeNegativeBins(h_eff_data_endcap2);
    removeNegativeBins(h_eff_data_eta);

    TH1D *h_ineff_data_barrel = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_1-PR_subtract_barrel")));
    TH1D *h_ineff_data_endcap = ((TH1D*)(h_pT_endcap_fail[0]->Clone("h_1-PR_subtract_endcap")));
    TH1D *h_ineff_data_endcap2 = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_1-PR_subtract_endcap2")));
    TH1D *h_ineff_data_eta = ((TH1D*)(h_eta_fail[0]->Clone("h_1-PR_subtract_eta")));
    TH1D *h_ineff_MC_barrel = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_1-PR_MC_barrel")));
    TH1D *h_ineff_MC_endcap = ((TH1D*)(h_pT_endcap_fail[1]->Clone("h_1-PR_MC_endcap")));
    TH1D *h_ineff_MC_endcap2 = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_1-PR_MC_endcap2")));
    TH1D *h_ineff_MC_eta = ((TH1D*)(h_eta_fail[1]->Clone("h_1-PR_MC_eta")));
    TH1D *h_ineff_bkgr_barrel = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_1-PR_bkgr_barrel")));
    TH1D *h_ineff_bkgr_endcap = ((TH1D*)(h_pT_endcap_fail[2]->Clone("h_1-PR_bkgr_endcap")));
    TH1D *h_ineff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2")));
    TH1D *h_ineff_bkgr_eta = ((TH1D*)(h_eta_fail[2]->Clone("h_1-PR_bkgr_eta")));
    TH1D *h_ineff_bkgf_barrel = ((TH1D*)(h_pT_barrel_fail[3]->Clone("h_1-PR_bkgf_barrel")));
    TH1D *h_ineff_bkgf_endcap = ((TH1D*)(h_pT_endcap_fail[3]->Clone("h_1-PR_bkgf_endcap")));
    TH1D *h_ineff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2")));
    TH1D *h_ineff_bkgf_eta = ((TH1D*)(h_eta_fail[3]->Clone("h_1-PR_bkgf_eta")));
    h_ineff_data_barrel->Add(h_ineff_bkgf_barrel, -1);
    h_ineff_data_endcap->Add(h_ineff_bkgf_endcap, -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgf_endcap2, -1);
    h_ineff_data_eta->Add(h_ineff_bkgf_eta, -1);
    h_ineff_data_barrel->Add(h_ineff_bkgr_barrel, -1);
    h_ineff_data_endcap->Add(h_ineff_bkgr_endcap, -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgr_endcap2, -1);
    h_ineff_data_eta->Add(h_ineff_bkgr_eta, -1);
    removeNegativeBins(h_ineff_data_barrel);
    removeNegativeBins(h_ineff_data_endcap);
    removeNegativeBins(h_ineff_data_endcap2);
    removeNegativeBins(h_ineff_data_eta);

    TH1D *h_eff_data_barrel_deno = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_subtract_barrel_deno")));
    TH1D *h_eff_data_endcap_deno = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_subtract_endcap_deno")));
    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_endcap2_deno")));
    TH1D *h_eff_data_eta_deno = ((TH1D*)(h_eff_data_eta->Clone("h_PR_subtract_eta_deno")));
    h_eff_data_barrel_deno->Add(h_ineff_data_barrel);
    h_eff_data_endcap_deno->Add(h_ineff_data_endcap);
    h_eff_data_endcap2_deno->Add(h_ineff_data_endcap2);
    h_eff_data_eta_deno->Add(h_ineff_data_eta);

    TH1D *h_eff_MC_barrel_deno = ((TH1D*)(h_eff_MC_barrel->Clone("h_PR_MC_barrel_deno")));
    TH1D *h_eff_MC_endcap_deno = ((TH1D*)(h_eff_MC_endcap->Clone("h_PR_MC_endcap_deno")));
    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_eff_MC_endcap2->Clone("h_PR_MC_endcap2_deno")));
    TH1D *h_eff_MC_eta_deno = ((TH1D*)(h_eff_MC_eta->Clone("h_PR_MC_eta_deno")));
    h_eff_MC_barrel_deno->Add(h_ineff_MC_barrel);
    h_eff_MC_endcap_deno->Add(h_ineff_MC_endcap);
    h_eff_MC_endcap2_deno->Add(h_ineff_MC_endcap2);
    h_eff_MC_eta_deno->Add(h_ineff_MC_eta);

//    TH1D *h_eff_data_barrel_deno = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_PR_subtract_barrel_deno")));
//    TH1D *h_eff_data_endcap_deno = ((TH1D*)(h_pT_endcap_fail[0]->Clone("h_PR_subtract_endcap_deno")));
//    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_PR_subtract_endcap2_deno")));
//    TH1D *h_eff_data_eta_deno = ((TH1D*)(h_eta_fail[0]->Clone("h_PR_subtract_eta_deno")));
//    TH1D *h_eff_MC_barrel_deno = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_PR_MC_barrel_deno")));
//    TH1D *h_eff_MC_endcap_deno = ((TH1D*)(h_pT_endcap_fail[1]->Clone("h_PR_MC_endcap_deno")));
//    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_PR_MC_endcap2_deno")));
//    TH1D *h_eff_MC_eta_deno = ((TH1D*)(h_eta_fail[1]->Clone("h_PR_MC_eta_deno")));
//    TH1D *h_eff_bkgr_barrel_deno = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_1-PR_bkgr_barrel_deno")));
//    TH1D *h_eff_bkgr_endcap_deno = ((TH1D*)(h_pT_endcap_fail[2]->Clone("h_1-PR_bkgr_endcap_deno")));
//    TH1D *h_eff_bkgr_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2_deno")));
//    TH1D *h_eff_bkgr_eta_deno = ((TH1D*)(h_eta_fail[2]->Clone("h_1-PR_bkgr_eta_deno")));
//    TH1D *h_eff_bkgf_barrel_deno = ((TH1D*)(h_pT_barrel_fail[3]->Clone("h_1-PR_bkgf_barrel_deno")));
//    TH1D *h_eff_bkgf_endcap_deno = ((TH1D*)(h_pT_endcap_fail[3]->Clone("h_1-PR_bkgf_endcap_deno")));
//    TH1D *h_eff_bkgf_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2_deno")));
//    TH1D *h_eff_bkgf_eta_deno = ((TH1D*)(h_eta_fail[3]->Clone("h_1-PR_bkgf_eta_deno")));
//    h_eff_data_barrel_deno->Add(h_pT_barrel_pass[0]);
//    h_eff_data_endcap_deno->Add(h_pT_endcap_pass[0]);
//    h_eff_data_endcap2_deno->Add(h_pT_endcap2_pass[0]);
//    h_eff_data_eta_deno->Add(h_eta_pass[0]);
//    h_eff_MC_barrel_deno->Add(h_pT_barrel_pass[1]);
//    h_eff_MC_endcap_deno->Add(h_pT_endcap_pass[1]);
//    h_eff_MC_endcap2_deno->Add(h_pT_endcap2_pass[1]);
//    h_eff_MC_eta_deno->Add(h_eta_pass[1]);
//    h_eff_bkgr_barrel_deno->Add(h_pT_barrel_pass[2]);
//    h_eff_bkgr_endcap_deno->Add(h_pT_endcap_pass[2]);
//    h_eff_bkgr_endcap2_deno->Add(h_pT_endcap2_pass[2]);
//    h_eff_bkgr_eta_deno->Add(h_eta_pass[2]);
//    h_eff_bkgf_barrel_deno->Add(h_pT_barrel_pass[3]);
//    h_eff_bkgf_endcap_deno->Add(h_pT_endcap_pass[3]);
//    h_eff_bkgf_endcap2_deno->Add(h_pT_endcap2_pass[3]);
//    h_eff_bkgf_eta_deno->Add(h_eta_pass[3]);
//    h_eff_data_barrel_deno->Add(h_eff_bkgf_barrel_deno, -1);
//    h_eff_data_endcap_deno->Add(h_eff_bkgf_endcap_deno, -1);
//    h_eff_data_endcap2_deno->Add(h_eff_bkgf_endcap2_deno, -1);
//    h_eff_data_eta_deno->Add(h_eff_bkgf_eta_deno, -1);
//    h_eff_data_barrel_deno->Add(h_eff_bkgr_barrel_deno, -1);
//    h_eff_data_endcap_deno->Add(h_eff_bkgr_endcap_deno, -1);
//    h_eff_data_endcap2_deno->Add(h_eff_bkgr_endcap2_deno, -1);
//    h_eff_data_eta_deno->Add(h_eff_bkgr_eta_deno, -1);

    h_eff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_eff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_eff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_eff_data_eta->Divide(h_eff_data_eta_deno);
    h_eff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_eff_MC_endcap->Divide(h_eff_MC_endcap_deno);
    h_eff_MC_endcap2->Divide(h_eff_MC_endcap2_deno);
    h_eff_MC_eta->Divide(h_eff_MC_eta_deno);
    h_ineff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_ineff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_ineff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_ineff_data_eta->Divide(h_eff_data_eta_deno);
    h_ineff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_ineff_MC_endcap->Divide(h_eff_MC_endcap_deno);
    h_ineff_MC_endcap2->Divide(h_eff_MC_endcap2_deno);
    h_ineff_MC_eta->Divide(h_eff_MC_eta_deno);

    h_eff_data_barrel->SetDirectory(0);
    h_eff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel->SetMarkerColor(kBlack);
    h_eff_data_barrel->SetLineColor(kBlack);
    h_eff_data_endcap->SetDirectory(0);
    h_eff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap->SetMarkerColor(kBlack);
    h_eff_data_endcap->SetLineColor(kBlack);
    h_eff_data_endcap2->SetDirectory(0);
    h_eff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2->SetMarkerColor(kBlack);
    h_eff_data_endcap2->SetLineColor(kBlack);
    h_eff_data_eta->SetDirectory(0);
    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_data_eta->SetMarkerColor(kBlack);
    h_eff_data_eta->SetLineColor(kBlack);
    h_eff_MC_barrel->SetDirectory(0);
    h_eff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel->SetMarkerColor(kRed);
    h_eff_MC_barrel->SetLineColor(kRed);
    h_eff_MC_endcap->SetDirectory(0);
    h_eff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap->SetMarkerColor(kRed);
    h_eff_MC_endcap->SetLineColor(kRed);
    h_eff_MC_endcap2->SetDirectory(0);
    h_eff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap2->SetMarkerColor(kRed);
    h_eff_MC_endcap2->SetLineColor(kRed);
    h_eff_MC_eta->SetDirectory(0);
    h_eff_MC_eta->SetMarkerStyle(kFullSquare);
    h_eff_MC_eta->SetMarkerColor(kRed);
    h_eff_MC_eta->SetLineColor(kRed);

    h_ineff_data_barrel->SetDirectory(0);
    h_ineff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel->SetMarkerColor(kBlack);
    h_ineff_data_barrel->SetLineColor(kBlack);
    h_ineff_data_endcap->SetDirectory(0);
    h_ineff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap->SetMarkerColor(kBlack);
    h_ineff_data_endcap->SetLineColor(kBlack);
    h_ineff_data_endcap2->SetDirectory(0);
    h_ineff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap2->SetMarkerColor(kBlack);
    h_ineff_data_endcap2->SetLineColor(kBlack);
    h_ineff_data_eta->SetDirectory(0);
    h_ineff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_eta->SetMarkerColor(kBlack);
    h_ineff_data_eta->SetLineColor(kBlack);
    h_ineff_MC_barrel->SetDirectory(0);
    h_ineff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel->SetMarkerColor(kRed);
    h_ineff_MC_barrel->SetLineColor(kRed);
    h_ineff_MC_endcap->SetDirectory(0);
    h_ineff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap->SetMarkerColor(kRed);
    h_ineff_MC_endcap->SetLineColor(kRed);
    h_ineff_MC_endcap2->SetDirectory(0);
    h_ineff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap2->SetMarkerColor(kRed);
    h_ineff_MC_endcap2->SetLineColor(kRed);
    h_ineff_MC_eta->SetDirectory(0);
    h_ineff_MC_eta->SetMarkerStyle(kFullSquare);
    h_ineff_MC_eta->SetMarkerColor(kRed);
    h_ineff_MC_eta->SetLineColor(kRed);

    TH1D *h_eff_ratio_barrel = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_ratio_barrel")));
    TH1D *h_eff_ratio_endcap = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_ratio_endcap")));
    TH1D *h_eff_ratio_endcap2 = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_ratio_endcap2")));
    TH1D *h_eff_ratio_eta = ((TH1D*)(h_eff_data_eta->Clone("h_PR_ratio_eta")));
    TH1D *h_ineff_ratio_barrel = ((TH1D*)(h_ineff_data_barrel->Clone("h_1-PR_ratio_barrel")));
    TH1D *h_ineff_ratio_endcap = ((TH1D*)(h_ineff_data_endcap->Clone("h_1-PR_ratio_endcap")));
    TH1D *h_ineff_ratio_endcap2 = ((TH1D*)(h_ineff_data_endcap2->Clone("h_1-PR_ratio_endcap2")));
    TH1D *h_ineff_ratio_eta = ((TH1D*)(h_ineff_data_eta->Clone("h_1-PR_ratio_eta")));
    h_eff_ratio_barrel->Divide(h_eff_MC_barrel);
    h_eff_ratio_endcap->Divide(h_eff_MC_endcap);
    h_eff_ratio_endcap2->Divide(h_eff_MC_endcap2);
    h_eff_ratio_eta->Divide(h_eff_MC_eta);
    h_eff_ratio_barrel->SetDirectory(0);
    h_eff_ratio_endcap->SetDirectory(0);
    h_eff_ratio_endcap->SetMarkerColor(kBlue);
    h_eff_ratio_endcap->SetLineColor(kBlue);
    h_eff_ratio_endcap2->SetDirectory(0);
    h_eff_ratio_endcap2->SetMarkerColor(kBlue);
    h_eff_ratio_endcap2->SetLineColor(kBlue);
    h_eff_ratio_eta->SetDirectory(0);
    h_eff_ratio_eta->SetMarkerColor(kBlue);
    h_eff_ratio_eta->SetLineColor(kBlue);
    h_ineff_ratio_barrel->Divide(h_ineff_MC_barrel);
    h_ineff_ratio_endcap->Divide(h_ineff_MC_endcap);
    h_ineff_ratio_endcap2->Divide(h_ineff_MC_endcap2);
    h_ineff_ratio_eta->Divide(h_ineff_MC_eta);
    h_ineff_ratio_barrel->SetDirectory(0);
    h_ineff_ratio_endcap->SetDirectory(0);
    h_ineff_ratio_endcap->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap->SetLineColor(kBlue);
    h_ineff_ratio_endcap2->SetDirectory(0);
    h_ineff_ratio_endcap2->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap2->SetLineColor(kBlue);
    h_ineff_ratio_eta->SetDirectory(0);
    h_ineff_ratio_eta->SetMarkerColor(kBlue);
    h_ineff_ratio_eta->SetLineColor(kBlue);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_eff_barrel = new myRatioPlot_t("RP_PR_barrel", h_eff_MC_barrel, h_eff_data_barrel);
    myRatioPlot_t *RP_eff_endcap = new myRatioPlot_t("RP_PR_endcap", h_eff_MC_endcap, h_eff_data_endcap);
    myRatioPlot_t *RP_eff_endcap2 = new myRatioPlot_t("RP_PR_endcap2", h_eff_MC_endcap2, h_eff_data_endcap2);
    myRatioPlot_t *RP_eff_eta = new myRatioPlot_t("RP_PR_eta", h_eff_MC_eta, h_eff_data_eta);
    myRatioPlot_t *RP_ineff_barrel = new myRatioPlot_t("RP_1-PR_barrel", h_ineff_MC_barrel, h_ineff_data_barrel);
    myRatioPlot_t *RP_ineff_endcap = new myRatioPlot_t("RP_1-PR_endcap", h_ineff_MC_endcap, h_ineff_data_endcap);
    myRatioPlot_t *RP_ineff_endcap2 = new myRatioPlot_t("RP_1-PR_endcap2", h_ineff_MC_endcap2, h_ineff_data_endcap2);
    myRatioPlot_t *RP_ineff_eta = new myRatioPlot_t("RP_1-PR_eta", h_ineff_MC_eta, h_ineff_data_eta);

    RP_eff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_eta->SetPlots("#eta", -2.4, 2.4);
    RP_ineff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_eta->SetPlots("#eta", -2.4, 2.4);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_eff_data_barrel, "Data", "pl");
    legend->AddEntry(h_eff_MC_barrel, "DY MC", "pl");

    RP_eff_barrel->ImportLegend(legend);
    RP_eff_endcap->ImportLegend(legend);
    RP_eff_endcap2->ImportLegend(legend);
    RP_eff_eta->ImportLegend(legend);
    RP_ineff_barrel->ImportLegend(legend);
    RP_ineff_endcap->ImportLegend(legend);
    RP_ineff_endcap2->ImportLegend(legend);
    RP_ineff_eta->ImportLegend(legend);

//    RP_eff_barrel->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap2->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_eta->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_ineff_barrel->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_eta->Draw(0.01, 1, 1, "", "1-PR");

//    RP_eff_barrel->pad1->SetLogy(0);
//    RP_eff_endcap->pad1->SetLogy(0);
//    RP_eff_endcap2->pad1->SetLogy(0);
//    RP_eff_eta->pad1->SetLogy(0);
//    RP_ineff_barrel->pad1->SetLogy(0);
//    RP_ineff_endcap->pad1->SetLogy(0);
//    RP_ineff_endcap2->pad1->SetLogy(0);
//    RP_ineff_eta->pad1->SetLogy(0);

    TCanvas *c_PR_barrel = new TCanvas("c_PR_barrel", "c_PR_barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
    TH1D *h_eff_data_barrel_draw = ((TH1D*)(h_eff_data_barrel->Clone("h_eff_data_barrel_draw")));
    h_eff_data_barrel_draw->SetDirectory(0);
    h_eff_data_barrel_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel_draw->SetMarkerColor(kBlack);
    h_eff_data_barrel_draw->SetLineColor(kBlack);
    h_eff_data_barrel_draw->SetStats(kFALSE);
    h_eff_data_barrel_draw->SetTitle("");
    h_eff_data_barrel_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_data_barrel_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_barrel_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_barrel_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_barrel_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_barrel_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_barrel_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_barrel_draw->GetYaxis()->SetRangeUser(0.3, 1.1);
    h_eff_data_barrel_draw->Draw();
    // Fit
    TF1 *f_barrel_17to43 = new TF1("f_barrel_17to43","[0]+[1]*x+[2]*sqrt(x)",17,45);
    f_barrel_17to43->SetLineColor(kMagenta-5);
    TF1 *f_barrel_43to56 = new TF1("f_barrel_43to56","[0]+exp([1]*([2]-x))",40,70);
    f_barrel_43to56->SetParameter(0, 0.55);
    f_barrel_43to56->SetParameter(1, 0.55);
    f_barrel_43to56->SetParameter(2, 50);
    f_barrel_43to56->SetLineColor(kMagenta-5);
    TF1 *f_barrel_56to95 = new TF1("f_barrel_56to95","[0]+[1]*x+[2]*x^2",50,100);
    f_barrel_56to95->SetLineColor(kMagenta-5);
    TF1 *f_barrel_95to5000 = new TF1("f_barrel_95to5000","[0]+gaus(1)",90,5000);
    f_barrel_95to5000->SetParameter(0, 0.4);
    f_barrel_95to5000->SetParameter(1, 0.3);
    f_barrel_95to5000->SetParameter(2, 90);
    f_barrel_95to5000->SetParameter(3, 150);
    f_barrel_95to5000->SetLineColor(kMagenta-5);
    h_eff_data_barrel_draw->Fit("f_barrel_17to43", "R");
    h_eff_data_barrel_draw->Fit("f_barrel_43to56", "R");
    h_eff_data_barrel_draw->Fit("f_barrel_56to95", "R");
    h_eff_data_barrel_draw->Fit("f_barrel_95to5000", "R");
    f_barrel_17to43->Draw("same");
    f_barrel_43to56->Draw("same");
    f_barrel_56to95->Draw("same");
    f_barrel_95to5000->Draw("same");
    Double_t chi2_barrel_40 = f_barrel_17to43->GetChisquare();
    Double_t chi2_barrel_70 = f_barrel_43to56->GetChisquare();
    Double_t chi2_barrel_100 = f_barrel_56to95->GetChisquare();
    Double_t chi2_barrel_5000 = f_barrel_95to5000->GetChisquare();
    cout << "Barrel fit chi squares: " << chi2_barrel_40 << "  " << chi2_barrel_70<< "  " << chi2_barrel_100  << "  " << chi2_barrel_5000 << endl;
    cout << "\n\n";
    c_PR_barrel->Update();


    TCanvas *c_PR_endcap = new TCanvas("c_PR_endcap", "c_PR_endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap_draw = ((TH1D*)(h_eff_data_endcap->Clone("h_eff_data_endcap_draw")));
    h_eff_data_endcap_draw->SetDirectory(0);
    h_eff_data_endcap_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap_draw->SetLineColor(kBlack);
    h_eff_data_endcap_draw->SetStats(kFALSE);
    h_eff_data_endcap_draw->SetTitle("");
    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.3, 1.1);
    h_eff_data_endcap_draw->Draw();
    // Fit
    TF1 *f_endcap_17to42p5 = new TF1("f_endcap_17to42p5","[0]+[1]*x+[2]*sqrt(x)",17,43);
    f_endcap_17to42p5->SetLineColor(kAzure+1);
    TF1 *f_endcap_42p5to55 = new TF1("f_endcap_42p5to55","[0]+[1]/(x-[2])+[3]*x",40,60);
    f_endcap_42p5to55->SetLineColor(kAzure+1);
    TF1 *f_endcap_55to5000 = new TF1("f_endcap_55to5000","[0]+[1]/(1+exp([2]*(x-[3])))",50,5000);
    f_endcap_55to5000->SetParameter(0, h_eff_data_endcap_draw->GetBinContent(16));
    f_endcap_55to5000->SetParameter(1, h_eff_data_endcap_draw->GetBinContent(11)-h_eff_data_endcap_draw->GetBinContent(16));
    f_endcap_55to5000->SetParameter(2, 0.05);
    f_endcap_55to5000->SetParLimits(2, 0.04, 0.15);
    f_endcap_55to5000->SetParameter(3, 200);
    f_endcap_55to5000->SetParLimits(3, 200, 400);
    f_endcap_55to5000->SetLineColor(kAzure+1);
    h_eff_data_endcap_draw->Fit("f_endcap_17to42p5", "R");
    h_eff_data_endcap_draw->Fit("f_endcap_42p5to55", "R");
    h_eff_data_endcap_draw->Fit("f_endcap_55to5000", "R");
    f_endcap_17to42p5->Draw("same");
    f_endcap_42p5to55->Draw("same");
    f_endcap_55to5000->Draw("same");
    Double_t chi2_endcap_60 = f_endcap_17to42p5->GetChisquare();
    Double_t chi2_endcap_90 = f_endcap_42p5to55->GetChisquare();
    Double_t chi2_endcap_5000 = f_endcap_55to5000->GetChisquare();
    cout << "Endcap fit chi squares: " << chi2_endcap_60 << "  " << chi2_endcap_90 << "  " << chi2_endcap_5000 << endl;
    cout << "\n\n";
    c_PR_endcap->Update();


    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "c_PR_endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap2_draw = ((TH1D*)(h_eff_data_endcap2->Clone("h_eff_data_endcap2_draw")));
    h_eff_data_endcap2_draw->SetDirectory(0);
    h_eff_data_endcap2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap2_draw->SetLineColor(kBlack);
    h_eff_data_endcap2_draw->SetStats(kFALSE);
    h_eff_data_endcap2_draw->SetTitle("");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap2_draw->GetYaxis()->SetRangeUser(0, 0.8);
    h_eff_data_endcap2_draw->Draw();
    // Fit
    TF1 *f_endcap2_17to42p5 = new TF1("f_endcap2_17to42p5","[0]+[1]*sqrt(x)+[2]*x",17,47);
    f_endcap2_17to42p5->SetLineColor(kOrange+5);
    TF1 *f_endcap2_42p5to56 = new TF1("f_endcap2_42p5to56","[0]+[1]*x+[2]/(x-[3])",42,67);
    f_endcap2_42p5to56->SetParLimits(1, 0.001, 0.05);
    f_endcap2_42p5to56->SetParLimits(2, 50, 70);
    f_endcap2_42p5to56->SetParLimits(3, 10, 50);
    f_endcap2_42p5to56->SetLineColor(kOrange+5);
    TF1 *f_endcap2_56to76 = new TF1("f_endcap2_56to76","[0]+gaus(1)",56,110);
    f_endcap2_56to76->FixParameter(0, h_eff_data_endcap2_draw->GetBinContent(10));
    f_endcap2_56to76->SetParLimits(1, 0.03, 0.07);
    f_endcap2_56to76->SetParLimits(2, 74, 76);
    f_endcap2_56to76->SetParLimits(3, 7, 15);
    f_endcap2_56to76->SetLineColor(kOrange+5);
    TF1 *f_endcap2_76to5000 = new TF1("f_endcap2_76to5000","[0]+exp([1]*(x-[2]))",70,5000);
    f_endcap2_76to5000->FixParameter(0, h_eff_data_endcap2_draw->GetBinContent(16));
    f_endcap2_76to5000->SetParameter(1, -0.3);
    f_endcap2_76to5000->SetParameter(2, 30);
    f_endcap2_76to5000->SetLineColor(kOrange+5);
    h_eff_data_endcap2_draw->Fit("f_endcap2_17to42p5", "R");
    h_eff_data_endcap2_draw->Fit("f_endcap2_42p5to56", "R");
    h_eff_data_endcap2_draw->Fit("f_endcap2_56to76", "R");
    h_eff_data_endcap2_draw->Fit("f_endcap2_76to5000", "R");
    f_endcap2_17to42p5->Draw("same");
    f_endcap2_42p5to56->Draw("same");
    f_endcap2_56to76->Draw("same");
    f_endcap2_76to5000->Draw("same");
    Double_t chi2_endcap2_45 = f_endcap2_17to42p5->GetChisquare();
    Double_t chi2_endcap2_65 = f_endcap2_42p5to56->GetChisquare();
    Double_t chi2_endcap2_90 = f_endcap2_56to76->GetChisquare();
    Double_t chi2_endcap2_5000 = f_endcap2_76to5000->GetChisquare();
    cout << "Endcap fit chi squares: " << chi2_endcap2_45 << "  " << chi2_endcap2_65 << "  " << chi2_endcap2_90 << "  " << chi2_endcap2_5000 << endl;
    cout << "\n\n";
    c_PR_endcap->Update();


    TCanvas *c_PR_allin1 = new TCanvas("c_PR_allin1", "c_PR_allin1", 800, 800);
    c_PR_allin1->cd();
    c_PR_allin1->SetGrid(1);
    c_PR_allin1->SetRightMargin(0.05);
    c_PR_allin1->SetTopMargin(0.05);
    c_PR_allin1->SetBottomMargin(0.12);
    c_PR_allin1->SetLeftMargin(0.13);
    h_eff_data_endcap_draw->SetDirectory(0);
    h_eff_data_endcap_draw->SetTitle("");
    h_eff_data_endcap_draw->SetMarkerStyle(kFullSquare);
    h_eff_data_endcap_draw->SetMarkerColor(kBlue);
    h_eff_data_endcap_draw->SetLineColor(kBlue);
    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.2, 0.9);
    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(0, 5000);
    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent();
    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels();
    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->SetDirectory(0);
    h_eff_data_endcap2_draw->SetMarkerStyle(33);
    h_eff_data_endcap2_draw->SetMarkerSize(1.5);
    h_eff_data_endcap2_draw->SetMarkerColor(kOrange-3);
    h_eff_data_endcap2_draw->SetLineColor(kOrange-3);
    h_eff_data_endcap_draw->Draw();
    h_eff_data_endcap2_draw->Draw("same");
    h_eff_data_barrel->Draw("same");
    TLegend *legend1 = new TLegend(0.5, 0.72, 0.95, 0.95);
    legend1->AddEntry(h_eff_data_barrel, "|#eta| < 1.4442", "LP");
    legend1->AddEntry(h_eff_data_endcap_draw, "1.566 < |#eta| < 2.1", "LP");
    legend1->AddEntry(h_eff_data_endcap2_draw, "|#eta| #geq 2.1", "LP");
    legend1->Draw();
    f_barrel_17to43->Draw("same");
    f_barrel_43to56->Draw("same");
    f_barrel_56to95->Draw("same");
    f_barrel_95to5000->Draw("same");
    f_endcap_17to42p5->Draw("same");
    f_endcap_42p5to55->Draw("same");
    f_endcap_55to5000->Draw("same");
    f_endcap2_17to42p5->Draw("same");
    f_endcap2_42p5to56->Draw("same");
    f_endcap2_56to76->Draw("same");
    f_endcap2_76to5000->Draw("same");
    c_PR_allin1->SetLogx();
    c_PR_allin1->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "c_PR_eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
    h_eff_data_eta->SetDirectory(0);
    h_eff_data_eta->SetTitle("");
    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_data_eta->SetMarkerColor(kBlack);
    h_eff_data_eta->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_eta->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_eta->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_eta->GetYaxis()->SetRangeUser(0.3, 1);
    h_eff_data_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_eff_data_eta->GetXaxis()->SetNoExponent();
    h_eff_data_eta->GetXaxis()->SetMoreLogLabels();
    h_eff_data_eta->GetXaxis()->SetTitle("#eta");
    h_eff_data_eta->GetXaxis()->SetTitleOffset(1);
    h_eff_data_eta->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_eta->Draw();
    TLegend *legend_eta = new TLegend(0.13, 0.77, 0.6, 0.95);
    legend_eta->AddEntry(h_eff_data_eta, "Subtraction", "LP");
    legend_eta->Draw();
    c_PR_eta->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << inName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << inName << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/FR/Electron/PromptRate_electron.root", "RECREATE");
    f_out->cd();
    h_eff_data_barrel->Write();
    h_eff_data_endcap->Write();
    h_eff_data_endcap2->Write();
    h_eff_data_eta->Write();
    h_eff_MC_barrel->Write();
    h_eff_MC_endcap->Write();
    h_eff_MC_endcap2->Write();
    h_eff_MC_eta->Write();
    h_ineff_data_barrel->Write();
    h_ineff_data_endcap->Write();
    h_ineff_data_endcap2->Write();
    h_ineff_data_eta->Write();
    h_ineff_MC_barrel->Write();
    h_ineff_MC_endcap->Write();
    h_ineff_MC_endcap2->Write();
    h_ineff_MC_eta->Write();
    h_eff_ratio_barrel->Write();
    h_eff_ratio_endcap->Write();
    h_eff_ratio_endcap2->Write();
    h_eff_ratio_eta->Write();
    h_ineff_ratio_barrel->Write();
    h_ineff_ratio_endcap->Write();
    h_ineff_ratio_endcap2->Write();
    h_ineff_ratio_eta->Write();
    f_barrel_17to43->Write();
    f_barrel_43to56->Write();
    f_barrel_56to95->Write();
    f_barrel_95to5000->Write();
    f_endcap_17to42p5->Write();
    f_endcap_42p5to55->Write();
    f_endcap_55to5000->Write();
    f_endcap2_17to42p5->Write();
    f_endcap2_42p5to56->Write();
    f_endcap2_56to76->Write();
    f_endcap2_76to5000->Write();
    f_out->Close();

} // End of E_EstimatePR


void E_EstPR_alt()
{
    TString inName = "/media/sf_DATA/FR/Electron/PR_Hist_E_alt.root";
    TFile *f = new TFile(inName, "READ");

    TH1D *h_pT_barrel_pass[4],
         *h_pT_endcap_pass[4],
         *h_pT_endcap2_pass[4],
         *h_eta_pass[4],
         *h_pT_barrel_fail[4],
         *h_pT_endcap_fail[4],
         *h_pT_endcap2_fail[4],
         *h_eta_fail[4];
    TString type[4] = {"data", "wjet", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_endcap2_pass_"+type[i], h_pT_endcap2_pass[i]);
        f->GetObject("h_eta_pass_"+type[i], h_eta_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        f->GetObject("h_pT_endcap2_fail_"+type[i], h_pT_endcap2_fail[i]);
        f->GetObject("h_eta_fail_"+type[i], h_eta_fail[i]);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_endcap2_pass[i]->SetDirectory(0);
        h_eta_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
        h_pT_endcap2_fail[i]->SetDirectory(0);
        h_eta_fail[i]->SetDirectory(0);
    }

    TH1D *h_eff_data_barrel = ((TH1D*)(h_pT_barrel_pass[0]->Clone("h_PR_subtract_barrel")));
    TH1D *h_eff_data_endcap = ((TH1D*)(h_pT_endcap_pass[0]->Clone("h_PR_subtract_endcap")));
    TH1D *h_eff_data_endcap2 = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_subtract_endcap2")));
    TH1D *h_eff_data_eta = ((TH1D*)(h_eta_pass[0]->Clone("h_PR_subtract_eta")));
    TH1D *h_eff_MC_barrel = ((TH1D*)(h_pT_barrel_pass[1]->Clone("h_PR_MC_barrel")));
    TH1D *h_eff_MC_endcap = ((TH1D*)(h_pT_endcap_pass[1]->Clone("h_PR_MC_endcap")));
    TH1D *h_eff_MC_endcap2 = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_MC_endcap2")));
    TH1D *h_eff_MC_eta = ((TH1D*)(h_eta_pass[1]->Clone("h_PR_MC_eta")));
    TH1D *h_eff_bkgr_barrel = ((TH1D*)(h_pT_barrel_pass[2]->Clone("h_PR_bkgr_barrel")));
    TH1D *h_eff_bkgr_endcap = ((TH1D*)(h_pT_endcap_pass[2]->Clone("h_PR_bkgr_endcap")));
    TH1D *h_eff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_pass[2]->Clone("h_PR_bkgr_endcap2")));
    TH1D *h_eff_bkgr_eta = ((TH1D*)(h_eta_pass[2]->Clone("h_PR_bkgr_eta")));
    TH1D *h_eff_bkgf_barrel = ((TH1D*)(h_pT_barrel_pass[3]->Clone("h_PR_bkgf_barrel")));
    TH1D *h_eff_bkgf_endcap = ((TH1D*)(h_pT_endcap_pass[3]->Clone("h_PR_bkgf_endcap")));
    TH1D *h_eff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_pass[3]->Clone("h_PR_bkgf_endcap2")));
    TH1D *h_eff_bkgf_eta = ((TH1D*)(h_eta_pass[3]->Clone("h_PR_bkgf_eta")));
//    h_eff_data_barrel->Add(h_eff_bkgf_barrel, -1);
//    h_eff_data_endcap->Add(h_eff_bkgf_endcap, -1);
//    h_eff_data_endcap2->Add(h_eff_bkgf_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgf_eta, -1);
    h_eff_data_barrel->Add(h_eff_bkgr_barrel, -1);
    h_eff_data_endcap->Add(h_eff_bkgr_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgr_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgr_eta, -1);
    removeNegativeBins(h_eff_data_barrel);
    removeNegativeBins(h_eff_data_endcap);
    removeNegativeBins(h_eff_data_endcap2);
    removeNegativeBins(h_eff_data_eta);

    TH1D *h_ineff_data_barrel = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_1-PR_subtract_barrel")));
    TH1D *h_ineff_data_endcap = ((TH1D*)(h_pT_endcap_fail[0]->Clone("h_1-PR_subtract_endcap")));
    TH1D *h_ineff_data_endcap2 = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_1-PR_subtract_endcap2")));
    TH1D *h_ineff_data_eta = ((TH1D*)(h_eta_fail[0]->Clone("h_1-PR_subtract_eta")));
    TH1D *h_ineff_MC_barrel = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_1-PR_MC_barrel")));
    TH1D *h_ineff_MC_endcap = ((TH1D*)(h_pT_endcap_fail[1]->Clone("h_1-PR_MC_endcap")));
    TH1D *h_ineff_MC_endcap2 = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_1-PR_MC_endcap2")));
    TH1D *h_ineff_MC_eta = ((TH1D*)(h_eta_fail[1]->Clone("h_1-PR_MC_eta")));
    TH1D *h_ineff_bkgr_barrel = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_1-PR_bkgr_barrel")));
    TH1D *h_ineff_bkgr_endcap = ((TH1D*)(h_pT_endcap_fail[2]->Clone("h_1-PR_bkgr_endcap")));
    TH1D *h_ineff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2")));
    TH1D *h_ineff_bkgr_eta = ((TH1D*)(h_eta_fail[2]->Clone("h_1-PR_bkgr_eta")));
    TH1D *h_ineff_bkgf_barrel = ((TH1D*)(h_pT_barrel_fail[3]->Clone("h_1-PR_bkgf_barrel")));
    TH1D *h_ineff_bkgf_endcap = ((TH1D*)(h_pT_endcap_fail[3]->Clone("h_1-PR_bkgf_endcap")));
    TH1D *h_ineff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2")));
    TH1D *h_ineff_bkgf_eta = ((TH1D*)(h_eta_fail[3]->Clone("h_1-PR_bkgf_eta")));
//    h_ineff_data_barrel->Add(h_ineff_bkgf_barrel, -1);
//    h_ineff_data_endcap->Add(h_ineff_bkgf_endcap, -1);
//    h_ineff_data_endcap2->Add(h_ineff_bkgf_endcap2, -1);
    h_ineff_data_eta->Add(h_ineff_bkgf_eta, -1);
    h_ineff_data_barrel->Add(h_ineff_bkgr_barrel, -1);
    h_ineff_data_endcap->Add(h_ineff_bkgr_endcap, -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgr_endcap2, -1);
    h_ineff_data_eta->Add(h_ineff_bkgr_eta, -1);
    removeNegativeBins(h_ineff_data_barrel);
    removeNegativeBins(h_ineff_data_endcap);
    removeNegativeBins(h_ineff_data_endcap2);
    removeNegativeBins(h_ineff_data_eta);

    TH1D *h_eff_data_barrel_deno = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_subtract_barrel_deno")));
    TH1D *h_eff_data_endcap_deno = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_subtract_endcap_deno")));
    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_endcap2_deno")));
    TH1D *h_eff_data_eta_deno = ((TH1D*)(h_eff_data_eta->Clone("h_PR_subtract_eta_deno")));
    h_eff_data_barrel_deno->Add(h_ineff_data_barrel);
    h_eff_data_endcap_deno->Add(h_ineff_data_endcap);
    h_eff_data_endcap2_deno->Add(h_ineff_data_endcap2);
    h_eff_data_eta_deno->Add(h_ineff_data_eta);

    TH1D *h_eff_MC_barrel_deno = ((TH1D*)(h_eff_MC_barrel->Clone("h_PR_MC_barrel_deno")));
    TH1D *h_eff_MC_endcap_deno = ((TH1D*)(h_eff_MC_endcap->Clone("h_PR_MC_endcap_deno")));
    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_eff_MC_endcap2->Clone("h_PR_MC_endcap2_deno")));
    TH1D *h_eff_MC_eta_deno = ((TH1D*)(h_eff_MC_eta->Clone("h_PR_MC_eta_deno")));
    h_eff_MC_barrel_deno->Add(h_ineff_MC_barrel);
    h_eff_MC_endcap_deno->Add(h_ineff_MC_endcap);
    h_eff_MC_endcap2_deno->Add(h_ineff_MC_endcap2);
    h_eff_MC_eta_deno->Add(h_ineff_MC_eta);

//    TH1D *h_eff_data_barrel_deno = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_PR_subtract_barrel_deno")));
//    TH1D *h_eff_data_endcap_deno = ((TH1D*)(h_pT_endcap_fail[0]->Clone("h_PR_subtract_endcap_deno")));
//    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_PR_subtract_endcap2_deno")));
//    TH1D *h_eff_data_eta_deno = ((TH1D*)(h_eta_fail[0]->Clone("h_PR_subtract_eta_deno")));
//    TH1D *h_eff_MC_barrel_deno = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_PR_MC_barrel_deno")));
//    TH1D *h_eff_MC_endcap_deno = ((TH1D*)(h_pT_endcap_fail[1]->Clone("h_PR_MC_endcap_deno")));
//    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_PR_MC_endcap2_deno")));
//    TH1D *h_eff_MC_eta_deno = ((TH1D*)(h_eta_fail[1]->Clone("h_PR_MC_eta_deno")));
//    TH1D *h_eff_bkgr_barrel_deno = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_1-PR_bkgr_barrel_deno")));
//    TH1D *h_eff_bkgr_endcap_deno = ((TH1D*)(h_pT_endcap_fail[2]->Clone("h_1-PR_bkgr_endcap_deno")));
//    TH1D *h_eff_bkgr_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2_deno")));
//    TH1D *h_eff_bkgr_eta_deno = ((TH1D*)(h_eta_fail[2]->Clone("h_1-PR_bkgr_eta_deno")));
//    TH1D *h_eff_bkgf_barrel_deno = ((TH1D*)(h_pT_barrel_fail[3]->Clone("h_1-PR_bkgf_barrel_deno")));
//    TH1D *h_eff_bkgf_endcap_deno = ((TH1D*)(h_pT_endcap_fail[3]->Clone("h_1-PR_bkgf_endcap_deno")));
//    TH1D *h_eff_bkgf_endcap2_deno = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2_deno")));
//    TH1D *h_eff_bkgf_eta_deno = ((TH1D*)(h_eta_fail[3]->Clone("h_1-PR_bkgf_eta_deno")));
//    h_eff_data_barrel_deno->Add(h_pT_barrel_pass[0]);
//    h_eff_data_endcap_deno->Add(h_pT_endcap_pass[0]);
//    h_eff_data_endcap2_deno->Add(h_pT_endcap2_pass[0]);
//    h_eff_data_eta_deno->Add(h_eta_pass[0]);
//    h_eff_MC_barrel_deno->Add(h_pT_barrel_pass[1]);
//    h_eff_MC_endcap_deno->Add(h_pT_endcap_pass[1]);
//    h_eff_MC_endcap2_deno->Add(h_pT_endcap2_pass[1]);
//    h_eff_MC_eta_deno->Add(h_eta_pass[1]);
//    h_eff_bkgr_barrel_deno->Add(h_pT_barrel_pass[2]);
//    h_eff_bkgr_endcap_deno->Add(h_pT_endcap_pass[2]);
//    h_eff_bkgr_endcap2_deno->Add(h_pT_endcap2_pass[2]);
//    h_eff_bkgr_eta_deno->Add(h_eta_pass[2]);
//    h_eff_bkgf_barrel_deno->Add(h_pT_barrel_pass[3]);
//    h_eff_bkgf_endcap_deno->Add(h_pT_endcap_pass[3]);
//    h_eff_bkgf_endcap2_deno->Add(h_pT_endcap2_pass[3]);
//    h_eff_bkgf_eta_deno->Add(h_eta_pass[3]);
//    h_eff_data_barrel_deno->Add(h_eff_bkgf_barrel_deno, -1);
//    h_eff_data_endcap_deno->Add(h_eff_bkgf_endcap_deno, -1);
//    h_eff_data_endcap2_deno->Add(h_eff_bkgf_endcap2_deno, -1);
//    h_eff_data_eta_deno->Add(h_eff_bkgf_eta_deno, -1);
//    h_eff_data_barrel_deno->Add(h_eff_bkgr_barrel_deno, -1);
//    h_eff_data_endcap_deno->Add(h_eff_bkgr_endcap_deno, -1);
//    h_eff_data_endcap2_deno->Add(h_eff_bkgr_endcap2_deno, -1);
//    h_eff_data_eta_deno->Add(h_eff_bkgr_eta_deno, -1);

    h_eff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_eff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_eff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_eff_data_eta->Divide(h_eff_data_eta_deno);
    h_eff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_eff_MC_endcap->Divide(h_eff_MC_endcap_deno);
    h_eff_MC_endcap2->Divide(h_eff_MC_endcap2_deno);
    h_eff_MC_eta->Divide(h_eff_MC_eta_deno);
    h_ineff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_ineff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_ineff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_ineff_data_eta->Divide(h_eff_data_eta_deno);
    h_ineff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_ineff_MC_endcap->Divide(h_eff_MC_endcap_deno);
    h_ineff_MC_endcap2->Divide(h_eff_MC_endcap2_deno);
    h_ineff_MC_eta->Divide(h_eff_MC_eta_deno);

    h_eff_data_barrel->SetDirectory(0);
    h_eff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel->SetMarkerColor(kBlack);
    h_eff_data_barrel->SetLineColor(kBlack);
    h_eff_data_endcap->SetDirectory(0);
    h_eff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap->SetMarkerColor(kBlack);
    h_eff_data_endcap->SetLineColor(kBlack);
    h_eff_data_endcap2->SetDirectory(0);
    h_eff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2->SetMarkerColor(kBlack);
    h_eff_data_endcap2->SetLineColor(kBlack);
    h_eff_data_eta->SetDirectory(0);
    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_data_eta->SetMarkerColor(kBlack);
    h_eff_data_eta->SetLineColor(kBlack);
    h_eff_MC_barrel->SetDirectory(0);
    h_eff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel->SetMarkerColor(kRed);
    h_eff_MC_barrel->SetLineColor(kRed);
    h_eff_MC_endcap->SetDirectory(0);
    h_eff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap->SetMarkerColor(kRed);
    h_eff_MC_endcap->SetLineColor(kRed);
    h_eff_MC_endcap2->SetDirectory(0);
    h_eff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap2->SetMarkerColor(kRed);
    h_eff_MC_endcap2->SetLineColor(kRed);
    h_eff_MC_eta->SetDirectory(0);
    h_eff_MC_eta->SetMarkerStyle(kFullSquare);
    h_eff_MC_eta->SetMarkerColor(kRed);
    h_eff_MC_eta->SetLineColor(kRed);

    h_ineff_data_barrel->SetDirectory(0);
    h_ineff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel->SetMarkerColor(kBlack);
    h_ineff_data_barrel->SetLineColor(kBlack);
    h_ineff_data_endcap->SetDirectory(0);
    h_ineff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap->SetMarkerColor(kBlack);
    h_ineff_data_endcap->SetLineColor(kBlack);
    h_ineff_data_endcap2->SetDirectory(0);
    h_ineff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap2->SetMarkerColor(kBlack);
    h_ineff_data_endcap2->SetLineColor(kBlack);
    h_ineff_data_eta->SetDirectory(0);
    h_ineff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_eta->SetMarkerColor(kBlack);
    h_ineff_data_eta->SetLineColor(kBlack);
    h_ineff_MC_barrel->SetDirectory(0);
    h_ineff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel->SetMarkerColor(kRed);
    h_ineff_MC_barrel->SetLineColor(kRed);
    h_ineff_MC_endcap->SetDirectory(0);
    h_ineff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap->SetMarkerColor(kRed);
    h_ineff_MC_endcap->SetLineColor(kRed);
    h_ineff_MC_endcap2->SetDirectory(0);
    h_ineff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap2->SetMarkerColor(kRed);
    h_ineff_MC_endcap2->SetLineColor(kRed);
    h_ineff_MC_eta->SetDirectory(0);
    h_ineff_MC_eta->SetMarkerStyle(kFullSquare);
    h_ineff_MC_eta->SetMarkerColor(kRed);
    h_ineff_MC_eta->SetLineColor(kRed);

    TH1D *h_eff_ratio_barrel = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_ratio_barrel")));
    TH1D *h_eff_ratio_endcap = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_ratio_endcap")));
    TH1D *h_eff_ratio_endcap2 = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_ratio_endcap2")));
    TH1D *h_eff_ratio_eta = ((TH1D*)(h_eff_data_eta->Clone("h_PR_ratio_eta")));
    TH1D *h_ineff_ratio_barrel = ((TH1D*)(h_ineff_data_barrel->Clone("h_1-PR_ratio_barrel")));
    TH1D *h_ineff_ratio_endcap = ((TH1D*)(h_ineff_data_endcap->Clone("h_1-PR_ratio_endcap")));
    TH1D *h_ineff_ratio_endcap2 = ((TH1D*)(h_ineff_data_endcap2->Clone("h_1-PR_ratio_endcap2")));
    TH1D *h_ineff_ratio_eta = ((TH1D*)(h_ineff_data_eta->Clone("h_1-PR_ratio_eta")));
    h_eff_ratio_barrel->Divide(h_eff_MC_barrel);
    h_eff_ratio_endcap->Divide(h_eff_MC_endcap);
    h_eff_ratio_endcap2->Divide(h_eff_MC_endcap2);
    h_eff_ratio_eta->Divide(h_eff_MC_eta);
    h_eff_ratio_barrel->SetDirectory(0);
    h_eff_ratio_endcap->SetDirectory(0);
    h_eff_ratio_endcap->SetMarkerColor(kBlue);
    h_eff_ratio_endcap->SetLineColor(kBlue);
    h_eff_ratio_endcap2->SetDirectory(0);
    h_eff_ratio_endcap2->SetMarkerColor(kBlue);
    h_eff_ratio_endcap2->SetLineColor(kBlue);
    h_eff_ratio_eta->SetDirectory(0);
    h_eff_ratio_eta->SetMarkerColor(kBlue);
    h_eff_ratio_eta->SetLineColor(kBlue);
    h_ineff_ratio_barrel->Divide(h_ineff_MC_barrel);
    h_ineff_ratio_endcap->Divide(h_ineff_MC_endcap);
    h_ineff_ratio_endcap2->Divide(h_ineff_MC_endcap2);
    h_ineff_ratio_eta->Divide(h_ineff_MC_eta);
    h_ineff_ratio_barrel->SetDirectory(0);
    h_ineff_ratio_endcap->SetDirectory(0);
    h_ineff_ratio_endcap->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap->SetLineColor(kBlue);
    h_ineff_ratio_endcap2->SetDirectory(0);
    h_ineff_ratio_endcap2->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap2->SetLineColor(kBlue);
    h_ineff_ratio_eta->SetDirectory(0);
    h_ineff_ratio_eta->SetMarkerColor(kBlue);
    h_ineff_ratio_eta->SetLineColor(kBlue);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_eff_barrel = new myRatioPlot_t("RP_PR_barrel", h_eff_MC_barrel, h_eff_data_barrel);
    myRatioPlot_t *RP_eff_endcap = new myRatioPlot_t("RP_PR_endcap", h_eff_MC_endcap, h_eff_data_endcap);
    myRatioPlot_t *RP_eff_endcap2 = new myRatioPlot_t("RP_PR_endcap2", h_eff_MC_endcap2, h_eff_data_endcap2);
    myRatioPlot_t *RP_eff_eta = new myRatioPlot_t("RP_PR_eta", h_eff_MC_eta, h_eff_data_eta);
    myRatioPlot_t *RP_ineff_barrel = new myRatioPlot_t("RP_1-PR_barrel", h_ineff_MC_barrel, h_ineff_data_barrel);
    myRatioPlot_t *RP_ineff_endcap = new myRatioPlot_t("RP_1-PR_endcap", h_ineff_MC_endcap, h_ineff_data_endcap);
    myRatioPlot_t *RP_ineff_endcap2 = new myRatioPlot_t("RP_1-PR_endcap2", h_ineff_MC_endcap2, h_ineff_data_endcap2);
    myRatioPlot_t *RP_ineff_eta = new myRatioPlot_t("RP_1-PR_eta", h_ineff_MC_eta, h_ineff_data_eta);

    RP_eff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_eta->SetPlots("#eta", -2.4, 2.4);
    RP_ineff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_eta->SetPlots("#eta", -2.4, 2.4);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_eff_data_barrel, "Data", "pl");
    legend->AddEntry(h_eff_MC_barrel, "DY MC", "pl");

    RP_eff_barrel->ImportLegend(legend);
    RP_eff_endcap->ImportLegend(legend);
    RP_eff_endcap2->ImportLegend(legend);
    RP_eff_eta->ImportLegend(legend);
    RP_ineff_barrel->ImportLegend(legend);
    RP_ineff_endcap->ImportLegend(legend);
    RP_ineff_endcap2->ImportLegend(legend);
    RP_ineff_eta->ImportLegend(legend);

//    RP_eff_barrel->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap2->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_eta->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_ineff_barrel->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_eta->Draw(0.01, 1, 1, "", "1-PR");

//    RP_eff_barrel->pad1->SetLogy(0);
//    RP_eff_endcap->pad1->SetLogy(0);
//    RP_eff_endcap2->pad1->SetLogy(0);
//    RP_eff_eta->pad1->SetLogy(0);
//    RP_ineff_barrel->pad1->SetLogy(0);
//    RP_ineff_endcap->pad1->SetLogy(0);
//    RP_ineff_endcap2->pad1->SetLogy(0);
//    RP_ineff_eta->pad1->SetLogy(0);

    TCanvas *c_PR_barrel = new TCanvas("c_PR_barrel", "c_PR_barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
//    TH1D *h_eff_data_barrel_draw = ((TH1D*)(h_eff_data_barrel->Clone("h_eff_data_barrel_draw")));
//    h_eff_data_barrel_draw->SetDirectory(0);
//    h_eff_data_barrel_draw->SetMarkerStyle(kFullDotLarge);
//    h_eff_data_barrel_draw->SetMarkerColor(kBlack);
//    h_eff_data_barrel_draw->SetLineColor(kBlack);
//    h_eff_data_barrel_draw->SetStats(kFALSE);
//    h_eff_data_barrel_draw->SetTitle("");
//    h_eff_data_barrel_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
//    h_eff_data_barrel_draw->GetXaxis()->SetTitleOffset(1);
//    h_eff_data_barrel_draw->GetXaxis()->SetTitleSize(0.05);
//    h_eff_data_barrel_draw->GetXaxis()->SetLabelSize(0.04);
//    h_eff_data_barrel_draw->GetYaxis()->SetTitle("Prompt rate");
//    h_eff_data_barrel_draw->GetYaxis()->SetTitleSize(0.05);
//    h_eff_data_barrel_draw->GetYaxis()->SetTitleOffset(1.25);
//    h_eff_data_barrel_draw->GetYaxis()->SetLabelSize(0.04);
//    h_eff_data_barrel_draw->GetXaxis()->SetNoExponent(1);
//    h_eff_data_barrel_draw->GetXaxis()->SetMoreLogLabels(1);
//    h_eff_data_barrel_draw->GetXaxis()->SetRangeUser(17, 3000);
//    h_eff_data_barrel_draw->GetYaxis()->SetRangeUser(0.3, 1.1);
//    h_eff_data_barrel_draw->Draw();

    TH1D *h_eff_MC_barrel_draw = ((TH1D*)(h_eff_MC_barrel->Clone("h_eff_MC_barrel_draw")));
    h_eff_MC_barrel_draw->SetDirectory(0);
    h_eff_MC_barrel_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_MC_barrel_draw->SetMarkerColor(kBlack);
    h_eff_MC_barrel_draw->SetLineColor(kBlack);
    h_eff_MC_barrel_draw->SetStats(kFALSE);
    h_eff_MC_barrel_draw->SetTitle("");
    h_eff_MC_barrel_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_MC_barrel_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_MC_barrel_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_MC_barrel_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_MC_barrel_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_MC_barrel_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_MC_barrel_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_MC_barrel_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_MC_barrel_draw->GetXaxis()->SetNoExponent(1);
    h_eff_MC_barrel_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_MC_barrel_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_MC_barrel_draw->GetYaxis()->SetRangeUser(0.5, 1.1);
    h_eff_MC_barrel_draw->Draw();
    c_PR_barrel->Update();


    TCanvas *c_PR_endcap = new TCanvas("c_PR_endcap", "c_PR_endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
//    TH1D *h_eff_data_endcap_draw = ((TH1D*)(h_eff_data_endcap->Clone("h_eff_data_endcap_draw")));
//    h_eff_data_endcap_draw->SetDirectory(0);
//    h_eff_data_endcap_draw->SetMarkerStyle(kFullDotLarge);
//    h_eff_data_endcap_draw->SetMarkerColor(kBlack);
//    h_eff_data_endcap_draw->SetLineColor(kBlack);
//    h_eff_data_endcap_draw->SetStats(kFALSE);
//    h_eff_data_endcap_draw->SetTitle("");
//    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
//    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
//    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
//    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
//    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent(1);
//    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
//    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(17, 3000);
//    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.3, 1.1);
//    h_eff_data_endcap_draw->Draw();

    TH1D *h_eff_MC_endcap_draw = ((TH1D*)(h_eff_MC_endcap->Clone("h_eff_MC_endcap_draw")));
    h_eff_MC_endcap_draw->SetDirectory(0);
    h_eff_MC_endcap_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_MC_endcap_draw->SetMarkerColor(kBlack);
    h_eff_MC_endcap_draw->SetLineColor(kBlack);
    h_eff_MC_endcap_draw->SetStats(kFALSE);
    h_eff_MC_endcap_draw->SetTitle("");
    h_eff_MC_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_MC_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_MC_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_MC_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_MC_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap_draw->GetXaxis()->SetNoExponent(1);
    h_eff_MC_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_MC_endcap_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_MC_endcap_draw->GetYaxis()->SetRangeUser(0.5, 1.1);
    h_eff_MC_endcap_draw->Draw();
    c_PR_endcap->Update();


    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "c_PR_endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
//    TH1D *h_eff_data_endcap2_draw = ((TH1D*)(h_eff_data_endcap2->Clone("h_eff_data_endcap2_draw")));
//    h_eff_data_endcap2_draw->SetDirectory(0);
//    h_eff_data_endcap2_draw->SetMarkerStyle(kFullDotLarge);
//    h_eff_data_endcap2_draw->SetMarkerColor(kBlack);
//    h_eff_data_endcap2_draw->SetLineColor(kBlack);
//    h_eff_data_endcap2_draw->SetStats(kFALSE);
//    h_eff_data_endcap2_draw->SetTitle("");
//    h_eff_data_endcap2_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
//    h_eff_data_endcap2_draw->GetXaxis()->SetTitleOffset(1);
//    h_eff_data_endcap2_draw->GetXaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap2_draw->GetXaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap2_draw->GetYaxis()->SetTitle("Prompt rate");
//    h_eff_data_endcap2_draw->GetYaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap2_draw->GetYaxis()->SetTitleOffset(1.25);
//    h_eff_data_endcap2_draw->GetYaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap2_draw->GetXaxis()->SetNoExponent(1);
//    h_eff_data_endcap2_draw->GetXaxis()->SetMoreLogLabels(1);
//    h_eff_data_endcap2_draw->GetXaxis()->SetRangeUser(17, 3000);
//    h_eff_data_endcap2_draw->GetYaxis()->SetRangeUser(0, 0.8);
//    h_eff_data_endcap2_draw->Draw();

    TH1D *h_eff_MC_endcap2_draw = ((TH1D*)(h_eff_MC_endcap2->Clone("h_eff_MC_endcap2_draw")));
    h_eff_MC_endcap2_draw->SetDirectory(0);
    h_eff_MC_endcap2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_MC_endcap2_draw->SetMarkerColor(kBlack);
    h_eff_MC_endcap2_draw->SetLineColor(kBlack);
    h_eff_MC_endcap2_draw->SetStats(kFALSE);
    h_eff_MC_endcap2_draw->SetTitle("");
    h_eff_MC_endcap2_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_MC_endcap2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_MC_endcap2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_MC_endcap2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_MC_endcap2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_MC_endcap2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_MC_endcap2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_MC_endcap2_draw->GetYaxis()->SetRangeUser(0.5, 1.1);
    h_eff_MC_endcap2_draw->Draw();
    c_PR_endcap->Update();


    TCanvas *c_PR_allin1 = new TCanvas("c_PR_allin1", "c_PR_allin1", 800, 800);
    c_PR_allin1->cd();
    c_PR_allin1->SetGrid(1);
    c_PR_allin1->SetRightMargin(0.05);
    c_PR_allin1->SetTopMargin(0.05);
    c_PR_allin1->SetBottomMargin(0.12);
    c_PR_allin1->SetLeftMargin(0.13);
//    h_eff_data_endcap_draw->SetDirectory(0);
//    h_eff_data_endcap_draw->SetTitle("");
//    h_eff_data_endcap_draw->SetMarkerStyle(kFullSquare);
//    h_eff_data_endcap_draw->SetMarkerColor(kBlue);
//    h_eff_data_endcap_draw->SetLineColor(kBlue);
//    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
//    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.12);
//    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.2, 0.9);
//    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(0, 5000);
//    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent();
//    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels();
//    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
//    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
//    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
//    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
//    h_eff_data_endcap2_draw->SetDirectory(0);
//    h_eff_data_endcap2_draw->SetMarkerStyle(33);
//    h_eff_data_endcap2_draw->SetMarkerSize(1.5);
//    h_eff_data_endcap2_draw->SetMarkerColor(kOrange-3);
//    h_eff_data_endcap2_draw->SetLineColor(kOrange-3);
//    h_eff_data_endcap_draw->Draw();
//    h_eff_data_endcap2_draw->Draw("same");
//    h_eff_data_barrel->Draw("same");
//    TLegend *legend1 = new TLegend(0.5, 0.72, 0.95, 0.95);
//    legend1->AddEntry(h_eff_data_barrel, "|#eta| < 1.4442", "LP");
//    legend1->AddEntry(h_eff_data_endcap_draw, "1.566 < |#eta| < 2.1", "LP");
//    legend1->AddEntry(h_eff_data_endcap2_draw, "|#eta| #geq 2.1", "LP");
//    legend1->Draw();

    h_eff_MC_endcap_draw->SetDirectory(0);
    h_eff_MC_endcap_draw->SetTitle("");
    h_eff_MC_endcap_draw->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap_draw->SetMarkerColor(kBlue);
    h_eff_MC_endcap_draw->SetLineColor(kBlue);
    h_eff_MC_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_MC_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap_draw->GetYaxis()->SetTitleOffset(1.12);
    h_eff_MC_endcap_draw->GetYaxis()->SetRangeUser(0.5, 1.1);
    h_eff_MC_endcap_draw->GetXaxis()->SetRangeUser(0, 5000);
    h_eff_MC_endcap_draw->GetXaxis()->SetNoExponent();
    h_eff_MC_endcap_draw->GetXaxis()->SetMoreLogLabels();
    h_eff_MC_endcap_draw->GetXaxis()->SetTitle("p_{T} (#font[12]{e}) [GeV/c]");
    h_eff_MC_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_MC_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_MC_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_MC_endcap2_draw->SetDirectory(0);
    h_eff_MC_endcap2_draw->SetMarkerStyle(33);
    h_eff_MC_endcap2_draw->SetMarkerSize(1.5);
    h_eff_MC_endcap2_draw->SetMarkerColor(kOrange-3);
    h_eff_MC_endcap2_draw->SetLineColor(kOrange-3);
    h_eff_MC_endcap_draw->Draw();
    h_eff_MC_endcap2_draw->Draw("same");
    h_eff_MC_barrel->Draw("same");
    TLegend *legend1 = new TLegend(0.6, 0.77, 0.95, 0.95);
    legend1->AddEntry(h_eff_MC_barrel, "|#eta| < 1.4442", "LP");
    legend1->AddEntry(h_eff_MC_endcap_draw, "1.566 < |#eta| < 2.1", "LP");
    legend1->AddEntry(h_eff_MC_endcap2_draw, "|#eta| #geq 2.1", "LP");
    legend1->Draw();
    c_PR_allin1->SetLogx();
    c_PR_allin1->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "c_PR_eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
//    h_eff_data_eta->SetDirectory(0);
//    h_eff_data_eta->SetTitle("");
//    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
//    h_eff_data_eta->SetMarkerColor(kBlack);
//    h_eff_data_eta->GetYaxis()->SetTitle("Prompt rate");
//    h_eff_data_eta->GetYaxis()->SetTitleSize(0.05);
//    h_eff_data_eta->GetYaxis()->SetLabelSize(0.04);
//    h_eff_data_eta->GetYaxis()->SetTitleOffset(1.12);
//    h_eff_data_eta->GetYaxis()->SetRangeUser(0.5, 1);
//    h_eff_data_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
//    h_eff_data_eta->GetXaxis()->SetNoExponent();
//    h_eff_data_eta->GetXaxis()->SetMoreLogLabels();
//    h_eff_data_eta->GetXaxis()->SetTitle("#eta");
//    h_eff_data_eta->GetXaxis()->SetTitleOffset(1);
//    h_eff_data_eta->GetXaxis()->SetTitleSize(0.05);
//    h_eff_data_eta->GetXaxis()->SetLabelSize(0.04);
//    h_eff_data_eta->Draw();
//    TLegend *legend_eta = new TLegend(0.13, 0.85, 0.4, 0.95);
//    legend_eta->AddEntry(h_eff_data_eta, "Subtraction", "LP");
//    legend_eta->Draw();
    h_eff_MC_eta->SetDirectory(0);
    h_eff_MC_eta->SetTitle("");
    h_eff_MC_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_MC_eta->SetMarkerColor(kBlack);
    h_eff_MC_eta->GetYaxis()->SetTitle("Prompt rate");
    h_eff_MC_eta->GetYaxis()->SetTitleSize(0.05);
    h_eff_MC_eta->GetYaxis()->SetLabelSize(0.04);
    h_eff_MC_eta->GetYaxis()->SetTitleOffset(1.12);
    h_eff_MC_eta->GetYaxis()->SetRangeUser(0.5, 1.1);
    h_eff_MC_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_eff_MC_eta->GetXaxis()->SetNoExponent();
    h_eff_MC_eta->GetXaxis()->SetMoreLogLabels();
    h_eff_MC_eta->GetXaxis()->SetTitle("#eta");
    h_eff_MC_eta->GetXaxis()->SetTitleOffset(1);
    h_eff_MC_eta->GetXaxis()->SetTitleSize(0.05);
    h_eff_MC_eta->GetXaxis()->SetLabelSize(0.04);
    h_eff_MC_eta->Draw();
    TLegend *legend_eta = new TLegend(0.13, 0.85, 0.4, 0.95);
    legend_eta->AddEntry(h_eff_MC_eta, "W+Jets MC", "LP");
    legend_eta->Draw();
    c_PR_eta->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << inName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << inName << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/FR/Electron/PromptRate_electron_alt.root", "RECREATE");
    f_out->cd();
    h_eff_data_barrel->Write();
    h_eff_data_endcap->Write();
    h_eff_data_endcap2->Write();
    h_eff_data_eta->Write();
    h_eff_MC_barrel->Write();
    h_eff_MC_endcap->Write();
    h_eff_MC_endcap2->Write();
    h_eff_MC_eta->Write();
    h_ineff_data_barrel->Write();
    h_ineff_data_endcap->Write();
    h_ineff_data_endcap2->Write();
    h_ineff_data_eta->Write();
    h_ineff_MC_barrel->Write();
    h_ineff_MC_endcap->Write();
    h_ineff_MC_endcap2->Write();
    h_ineff_MC_eta->Write();
    h_eff_ratio_barrel->Write();
    h_eff_ratio_endcap->Write();
    h_eff_ratio_endcap2->Write();
    h_eff_ratio_eta->Write();
    h_ineff_ratio_barrel->Write();
    h_ineff_ratio_endcap->Write();
    h_ineff_ratio_endcap2->Write();
    h_ineff_ratio_eta->Write();
    f_out->Close();

} // End of E_EstimatePR_alt



void Mu_EstPR()
{
    TString inName = "/media/sf_DATA/FR/Muon/PR_Hist_Mu.root";
    TFile *f = new TFile(inName, "READ");

    TH1D *h_pT_barrel_pass[4],
         *h_pT_barrel2_pass[4],
         *h_pT_endcap_pass[4],
         *h_pT_endcap2_pass[4],
         *h_eta_pass[4],
         *h_MET_pass[4],
         *h_MT_pass[4],
         *h_nJet_pass[4],
         *h_dRjet_pass[4],
         *h_dPhiMET_pass[4],
         *h_pT_barrel_fail[4],
         *h_pT_barrel2_fail[4],
         *h_pT_endcap_fail[4],
         *h_pT_endcap2_fail[4],
         *h_eta_fail[4],
         *h_MET_fail[4],
         *h_MT_fail[4],
         *h_nJet_fail[4],
         *h_dRjet_fail[4],
         *h_dPhiMET_fail[4];

    TString type[4] = {"data", "DY", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_barrel2_pass_"+type[i], h_pT_barrel2_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_endcap2_pass_"+type[i], h_pT_endcap2_pass[i]);
        f->GetObject("h_eta_pass_"+type[i], h_eta_pass[i]);
        f->GetObject("h_MET_pass_"+type[i], h_MET_pass[i]);
        f->GetObject("h_MT_pass_"+type[i], h_MT_pass[i]);
        f->GetObject("h_nJet_pass_"+type[i], h_nJet_pass[i]);
        f->GetObject("h_dRjet_pass_"+type[i], h_dRjet_pass[i]);
        f->GetObject("h_dPhiMET_pass_"+type[i], h_dPhiMET_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_barrel2_fail_"+type[i], h_pT_barrel2_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        f->GetObject("h_pT_endcap2_fail_"+type[i], h_pT_endcap2_fail[i]);
        f->GetObject("h_eta_fail_"+type[i], h_eta_fail[i]);
        f->GetObject("h_MET_fail_"+type[i], h_MET_fail[i]);
        f->GetObject("h_MT_fail_"+type[i], h_MT_fail[i]);
        f->GetObject("h_nJet_fail_"+type[i], h_nJet_fail[i]);
        f->GetObject("h_dRjet_fail_"+type[i], h_dRjet_fail[i]);
        f->GetObject("h_dPhiMET_fail_"+type[i], h_dPhiMET_fail[i]);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_barrel2_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_endcap2_pass[i]->SetDirectory(0);
        h_eta_pass[i]->SetDirectory(0);
        h_MET_pass[i]->SetDirectory(0);
        h_MT_pass[i]->SetDirectory(0);
        h_nJet_pass[i]->SetDirectory(0);
        h_dRjet_pass[i]->SetDirectory(0);
        h_dPhiMET_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_barrel2_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
        h_pT_endcap2_fail[i]->SetDirectory(0);
        h_eta_fail[i]->SetDirectory(0);
        h_MET_fail[i]->SetDirectory(0);
        h_MT_fail[i]->SetDirectory(0);
        h_nJet_fail[i]->SetDirectory(0);
        h_dRjet_fail[i]->SetDirectory(0);
        h_dPhiMET_fail[i]->SetDirectory(0);
    }


    // Subtraction and MC methods
    TH1D *h_eff_data_barrel = ((TH1D*)(h_pT_barrel_pass[0]->Clone("h_PR_subtract_barrel")));
    TH1D *h_eff_data_barrel2 = ((TH1D*)(h_pT_barrel2_pass[0]->Clone("h_PR_subtract_barrel2")));
    TH1D *h_eff_data_endcap = ((TH1D*)(h_pT_endcap_pass[0]->Clone("h_PR_subtract_endcap")));
    TH1D *h_eff_data_endcap2 = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_subtract_endcap2")));
    TH1D *h_eff_data_eta = ((TH1D*)(h_eta_pass[0]->Clone("h_PR_subtract_eta")));
    TH1D *h_eff_data_MET = ((TH1D*)(h_MET_pass[0]->Clone("h_PR_subtract_MET")));
    TH1D *h_eff_data_MT = ((TH1D*)(h_MT_pass[0]->Clone("h_PR_subtract_MT")));
    TH1D *h_eff_data_nJet = ((TH1D*)(h_nJet_pass[0]->Clone("h_PR_subtract_nJet")));
    TH1D *h_eff_data_dRjet = ((TH1D*)(h_dRjet_pass[0]->Clone("h_PR_subtract_dRjet")));
    TH1D *h_eff_data_dPhiMET = ((TH1D*)(h_dPhiMET_pass[0]->Clone("h_PR_subtract_dPhiMET")));
    TH1D *h_eff_MC_barrel = ((TH1D*)(h_pT_barrel_pass[1]->Clone("h_PR_MC_barrel")));
    TH1D *h_eff_MC_barrel2 = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_MC_barrel2")));
    TH1D *h_eff_MC_endcap = ((TH1D*)(h_pT_endcap_pass[1]->Clone("h_PR_MC_endcap")));
    TH1D *h_eff_MC_endcap2 = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_MC_endcap2")));
    TH1D *h_eff_MC_eta = ((TH1D*)(h_eta_pass[1]->Clone("h_PR_MC_eta")));
    TH1D *h_eff_MC_MET = ((TH1D*)(h_MET_pass[1]->Clone("h_PR_MC_MET")));
    TH1D *h_eff_MC_MT = ((TH1D*)(h_MT_pass[1]->Clone("h_PR_MC_MT")));
    TH1D *h_eff_MC_nJet = ((TH1D*)(h_nJet_pass[1]->Clone("h_PR_MC_nJet")));
    TH1D *h_eff_MC_dRjet = ((TH1D*)(h_dRjet_pass[1]->Clone("h_PR_MC_dRjet")));
    TH1D *h_eff_MC_dPhiMET = ((TH1D*)(h_dPhiMET_pass[1]->Clone("h_PR_MC_dPhiMET")));
    TH1D *h_eff_bkgr_barrel = ((TH1D*)(h_pT_barrel_pass[2]->Clone("h_PR_bkgr_barrel")));
    TH1D *h_eff_bkgr_barrel2 = ((TH1D*)(h_pT_barrel2_pass[2]->Clone("h_PR_bkgr_barrel2")));
    TH1D *h_eff_bkgr_endcap = ((TH1D*)(h_pT_endcap_pass[2]->Clone("h_PR_bkgr_endcap")));
    TH1D *h_eff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_pass[2]->Clone("h_PR_bkgr_endcap2")));
    TH1D *h_eff_bkgr_eta = ((TH1D*)(h_eta_pass[2]->Clone("h_PR_bkgr_eta")));
    TH1D *h_eff_bkgr_MET = ((TH1D*)(h_MET_pass[2]->Clone("h_PR_bkgr_MET")));
    TH1D *h_eff_bkgr_MT = ((TH1D*)(h_MT_pass[2]->Clone("h_PR_bkgr_MT")));
    TH1D *h_eff_bkgr_nJet = ((TH1D*)(h_nJet_pass[2]->Clone("h_PR_bkgr_nJet")));
    TH1D *h_eff_bkgr_dRjet = ((TH1D*)(h_dRjet_pass[2]->Clone("h_PR_bkgr_dRjet")));
    TH1D *h_eff_bkgr_dPhiMET = ((TH1D*)(h_dPhiMET_pass[2]->Clone("h_PR_bkgr_dPhiMET")));
    TH1D *h_eff_bkgf_barrel = ((TH1D*)(h_pT_barrel_pass[3]->Clone("h_PR_bkgf_barrel")));
    TH1D *h_eff_bkgf_barrel2 = ((TH1D*)(h_pT_barrel2_pass[3]->Clone("h_PR_bkgf_barrel2")));
    TH1D *h_eff_bkgf_endcap = ((TH1D*)(h_pT_endcap_pass[3]->Clone("h_PR_bkgf_endcap")));
    TH1D *h_eff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_pass[3]->Clone("h_PR_bkgf_endcap2")));
    TH1D *h_eff_bkgf_eta = ((TH1D*)(h_eta_pass[3]->Clone("h_PR_bkgf_eta")));
    TH1D *h_eff_bkgf_MET = ((TH1D*)(h_MET_pass[3]->Clone("h_PR_bkgf_MET")));
    TH1D *h_eff_bkgf_MT = ((TH1D*)(h_MT_pass[3]->Clone("h_PR_bkgf_MT")));
    TH1D *h_eff_bkgf_nJet = ((TH1D*)(h_nJet_pass[3]->Clone("h_PR_bkgf_nJet")));
    TH1D *h_eff_bkgf_dRjet = ((TH1D*)(h_dRjet_pass[3]->Clone("h_PR_bkgf_dRjet")));
    TH1D *h_eff_bkgf_dPhiMET = ((TH1D*)(h_dPhiMET_pass[3]->Clone("h_PR_bkgf_dPhiMET")));
    h_eff_data_barrel->Add(h_eff_bkgf_barrel, -1);
    h_eff_data_barrel2->Add(h_eff_bkgf_barrel2, -1);
    h_eff_data_endcap->Add(h_eff_bkgf_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgf_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgf_eta, -1);
    h_eff_data_MET->Add(h_eff_bkgf_MET, -1);
    h_eff_data_MT->Add(h_eff_bkgf_MT, -1);
    h_eff_data_nJet->Add(h_eff_bkgf_nJet, -1);
    h_eff_data_dRjet->Add(h_eff_bkgf_dRjet, -1);
    h_eff_data_dPhiMET->Add(h_eff_bkgf_dPhiMET, -1);
    h_eff_data_barrel->Add(h_eff_bkgr_barrel, -1);
    h_eff_data_barrel2->Add(h_eff_bkgr_barrel2, -1);
    h_eff_data_endcap->Add(h_eff_bkgr_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgr_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgr_eta, -1);
    h_eff_data_MET->Add(h_eff_bkgr_MET, -1);
    h_eff_data_MT->Add(h_eff_bkgr_MT, -1);
    h_eff_data_nJet->Add(h_eff_bkgr_nJet, -1);
    h_eff_data_dRjet->Add(h_eff_bkgr_dRjet, -1);
    h_eff_data_dPhiMET->Add(h_eff_bkgr_dPhiMET, -1);
    removeNegativeBins(h_eff_data_barrel);
    removeNegativeBins(h_eff_data_barrel2);
    removeNegativeBins(h_eff_data_endcap);
    removeNegativeBins(h_eff_data_endcap2);
    removeNegativeBins(h_eff_data_eta);
    removeNegativeBins(h_eff_data_MET);
    removeNegativeBins(h_eff_data_MT);
    removeNegativeBins(h_eff_data_nJet);
    removeNegativeBins(h_eff_data_dRjet);
    removeNegativeBins(h_eff_data_dPhiMET);

    TH1D *h_ineff_data_barrel  = ((TH1D*)(h_pT_barrel_fail [0]->Clone("h_1-PR_subtract_barrel")));
    TH1D *h_ineff_data_barrel2 = ((TH1D*)(h_pT_barrel2_fail[0]->Clone("h_1-PR_subtract2_barrel")));
    TH1D *h_ineff_data_endcap  = ((TH1D*)(h_pT_endcap_fail [0]->Clone("h_1-PR_subtract_endcap")));
    TH1D *h_ineff_data_endcap2 = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_1-PR_subtract_endcap2")));
    TH1D *h_ineff_data_eta     = ((TH1D*)(h_eta_fail       [0]->Clone("h_1-PR_subtract_eta")));
    TH1D *h_ineff_data_MET     = ((TH1D*)(h_MET_fail       [0]->Clone("h_1-PR_subtract_MET")));
    TH1D *h_ineff_data_MT      = ((TH1D*)(h_MT_fail        [0]->Clone("h_1-PR_subtract_MT")));
    TH1D *h_ineff_data_nJet    = ((TH1D*)(h_nJet_fail      [0]->Clone("h_1-PR_subtract_nJet")));
    TH1D *h_ineff_data_dRjet   = ((TH1D*)(h_dRjet_fail     [0]->Clone("h_1-PR_subtract_dRjet")));
    TH1D *h_ineff_data_dPhiMET = ((TH1D*)(h_dPhiMET_fail   [0]->Clone("h_1-PR_subtract_dPhiMET")));
    TH1D *h_ineff_MC_barrel    = ((TH1D*)(h_pT_barrel_fail [1]->Clone("h_1-PR_MC_barrel")));
    TH1D *h_ineff_MC_barrel2   = ((TH1D*)(h_pT_barrel2_fail[1]->Clone("h_1-PR_MC_barrel2")));
    TH1D *h_ineff_MC_endcap    = ((TH1D*)(h_pT_endcap_fail [1]->Clone("h_1-PR_MC_endcap")));
    TH1D *h_ineff_MC_endcap2   = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_1-PR_MC_endcap2")));
    TH1D *h_ineff_MC_eta       = ((TH1D*)(h_eta_fail       [1]->Clone("h_1-PR_MC_eta")));
    TH1D *h_ineff_MC_MET       = ((TH1D*)(h_MET_fail       [1]->Clone("h_1-PR_MC_MET")));
    TH1D *h_ineff_MC_MT        = ((TH1D*)(h_MT_fail        [1]->Clone("h_1-PR_MC_MT")));
    TH1D *h_ineff_MC_nJet      = ((TH1D*)(h_nJet_fail      [1]->Clone("h_1-PR_MC_nJet")));
    TH1D *h_ineff_MC_dRjet     = ((TH1D*)(h_dRjet_fail     [1]->Clone("h_1-PR_MC_dRjet")));
    TH1D *h_ineff_MC_dPhiMET   = ((TH1D*)(h_dPhiMET_fail   [1]->Clone("h_1-PR_MC_dPhiMET")));
    TH1D *h_ineff_bkgr_barrel  = ((TH1D*)(h_pT_barrel_fail [2]->Clone("h_1-PR_bkgr_barrel")));
    TH1D *h_ineff_bkgr_barrel2 = ((TH1D*)(h_pT_barrel2_fail[2]->Clone("h_1-PR_bkgr_barrel2")));
    TH1D *h_ineff_bkgr_endcap  = ((TH1D*)(h_pT_endcap_fail [2]->Clone("h_1-PR_bkgr_endcap")));
    TH1D *h_ineff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2")));
    TH1D *h_ineff_bkgr_eta     = ((TH1D*)(h_eta_fail       [2]->Clone("h_1-PR_bkgr_eta")));
    TH1D *h_ineff_bkgr_MET     = ((TH1D*)(h_MET_fail       [2]->Clone("h_1-PR_bkgr_MET")));
    TH1D *h_ineff_bkgr_MT      = ((TH1D*)(h_MT_fail        [2]->Clone("h_1-PR_bkgr_MT")));
    TH1D *h_ineff_bkgr_nJet    = ((TH1D*)(h_nJet_fail      [2]->Clone("h_1-PR_bkgr_nJet")));
    TH1D *h_ineff_bkgr_dRjet   = ((TH1D*)(h_dRjet_fail     [2]->Clone("h_1-PR_bkgr_dRjet")));
    TH1D *h_ineff_bkgr_dPhiMET = ((TH1D*)(h_dPhiMET_fail   [2]->Clone("h_1-PR_bkgr_dPhiMET")));
    TH1D *h_ineff_bkgf_barrel  = ((TH1D*)(h_pT_barrel_fail [3]->Clone("h_1-PR_bkgf_barrel")));
    TH1D *h_ineff_bkgf_barrel2 = ((TH1D*)(h_pT_barrel2_fail[3]->Clone("h_1-PR_bkgf_barrel2")));
    TH1D *h_ineff_bkgf_endcap  = ((TH1D*)(h_pT_endcap_fail [3]->Clone("h_1-PR_bkgf_endcap")));
    TH1D *h_ineff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2")));
    TH1D *h_ineff_bkgf_eta     = ((TH1D*)(h_eta_fail       [3]->Clone("h_1-PR_bkgf_eta")));
    TH1D *h_ineff_bkgf_MET     = ((TH1D*)(h_MET_fail       [3]->Clone("h_1-PR_bkgf_MET")));
    TH1D *h_ineff_bkgf_MT      = ((TH1D*)(h_MT_fail        [3]->Clone("h_1-PR_bkgf_MT")));
    TH1D *h_ineff_bkgf_nJet    = ((TH1D*)(h_nJet_fail      [3]->Clone("h_1-PR_bkgf_nJet")));
    TH1D *h_ineff_bkgf_dRjet   = ((TH1D*)(h_dRjet_fail     [3]->Clone("h_1-PR_bkgf_dRjet")));
    TH1D *h_ineff_bkgf_dPhiMET = ((TH1D*)(h_dPhiMET_fail   [3]->Clone("h_1-PR_bkgf_dPhiMET")));
    h_ineff_data_barrel ->Add(h_ineff_bkgf_barrel,  -1);
    h_ineff_data_barrel2->Add(h_ineff_bkgf_barrel2, -1);
    h_ineff_data_endcap ->Add(h_ineff_bkgf_endcap,  -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgf_endcap2, -1);
    h_ineff_data_eta    ->Add(h_ineff_bkgf_eta,     -1);
    h_ineff_data_MET    ->Add(h_ineff_bkgf_MET,     -1);
    h_ineff_data_MT     ->Add(h_ineff_bkgf_MT,      -1);
    h_ineff_data_nJet   ->Add(h_ineff_bkgf_nJet,    -1);
    h_ineff_data_dRjet  ->Add(h_ineff_bkgf_dRjet,   -1);
    h_ineff_data_dPhiMET->Add(h_ineff_bkgf_dPhiMET, -1);
    h_ineff_data_barrel ->Add(h_ineff_bkgr_barrel,  -1);
    h_ineff_data_barrel2->Add(h_ineff_bkgr_barrel2, -1);
    h_ineff_data_endcap ->Add(h_ineff_bkgr_endcap,  -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgr_endcap2, -1);
    h_ineff_data_eta    ->Add(h_ineff_bkgr_eta,     -1);
    h_ineff_data_MET    ->Add(h_ineff_bkgr_MET,     -1);
    h_ineff_data_MT     ->Add(h_ineff_bkgr_MT,      -1);
    h_ineff_data_nJet   ->Add(h_ineff_bkgr_nJet,    -1);
    h_ineff_data_dRjet  ->Add(h_ineff_bkgr_dRjet,   -1);
    h_ineff_data_dPhiMET->Add(h_ineff_bkgr_dPhiMET, -1);
    removeNegativeBins(h_ineff_data_barrel);
    removeNegativeBins(h_ineff_data_barrel2);
    removeNegativeBins(h_ineff_data_endcap);
    removeNegativeBins(h_ineff_data_endcap2);
    removeNegativeBins(h_ineff_data_eta);
    removeNegativeBins(h_ineff_data_MET);
    removeNegativeBins(h_ineff_data_MT);
    removeNegativeBins(h_ineff_data_nJet);
    removeNegativeBins(h_ineff_data_dRjet);
    removeNegativeBins(h_ineff_data_dPhiMET);

    TH1D *h_eff_data_barrel_deno  = ((TH1D*)(h_eff_data_barrel ->Clone("h_PR_subtract_barrel_deno")));
    TH1D *h_eff_data_barrel2_deno = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_subtract_barrel2_deno")));
    TH1D *h_eff_data_endcap_deno  = ((TH1D*)(h_eff_data_endcap ->Clone("h_PR_subtract_endcap_deno")));
    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_endcap2_deno")));
    TH1D *h_eff_data_eta_deno     = ((TH1D*)(h_eff_data_eta    ->Clone("h_PR_subtract_eta_deno")));
    TH1D *h_eff_data_MET_deno     = ((TH1D*)(h_eff_data_MET    ->Clone("h_PR_subtract_MET_deno")));
    TH1D *h_eff_data_MT_deno      = ((TH1D*)(h_eff_data_MT     ->Clone("h_PR_subtract_MT_deno")));
    TH1D *h_eff_data_nJet_deno    = ((TH1D*)(h_eff_data_nJet   ->Clone("h_PR_subtract_nJet_deno")));
    TH1D *h_eff_data_dRjet_deno   = ((TH1D*)(h_eff_data_dRjet  ->Clone("h_PR_subtract_dRjet_deno")));
    TH1D *h_eff_data_dPhiMET_deno = ((TH1D*)(h_eff_data_dPhiMET->Clone("h_PR_subtract_dPhiMET_deno")));
    h_eff_data_barrel_deno ->Add(h_ineff_data_barrel);
    h_eff_data_barrel2_deno->Add(h_ineff_data_barrel2);
    h_eff_data_endcap_deno ->Add(h_ineff_data_endcap);
    h_eff_data_endcap2_deno->Add(h_ineff_data_endcap2);
    h_eff_data_eta_deno    ->Add(h_ineff_data_eta);
    h_eff_data_MET_deno    ->Add(h_ineff_data_MET);
    h_eff_data_MT_deno     ->Add(h_ineff_data_MT);
    h_eff_data_nJet_deno   ->Add(h_ineff_data_nJet);
    h_eff_data_dRjet_deno  ->Add(h_ineff_data_dRjet);
    h_eff_data_dPhiMET_deno->Add(h_ineff_data_dPhiMET);

    TH1D *h_eff_MC_barrel_deno  = ((TH1D*)(h_eff_MC_barrel ->Clone("h_PR_MC_barrel_deno")));
    TH1D *h_eff_MC_barrel2_deno = ((TH1D*)(h_eff_MC_barrel2->Clone("h_PR_MC_barrel2_deno")));
    TH1D *h_eff_MC_endcap_deno  = ((TH1D*)(h_eff_MC_endcap ->Clone("h_PR_MC_endcap_deno")));
    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_eff_MC_endcap2->Clone("h_PR_MC_endcap2_deno")));
    TH1D *h_eff_MC_eta_deno     = ((TH1D*)(h_eff_MC_eta    ->Clone("h_PR_MC_eta_deno")));
    TH1D *h_eff_MC_MET_deno     = ((TH1D*)(h_eff_MC_MET    ->Clone("h_PR_MC_MET_deno")));
    TH1D *h_eff_MC_MT_deno      = ((TH1D*)(h_eff_MC_MT     ->Clone("h_PR_MC_MT_deno")));
    TH1D *h_eff_MC_nJet_deno    = ((TH1D*)(h_eff_MC_nJet   ->Clone("h_PR_MC_nJet_deno")));
    TH1D *h_eff_MC_dRjet_deno   = ((TH1D*)(h_eff_MC_dRjet  ->Clone("h_PR_MC_dRjet_deno")));
    TH1D *h_eff_MC_dPhiMET_deno = ((TH1D*)(h_eff_MC_dPhiMET->Clone("h_PR_MC_dPhiMET_deno")));
    h_eff_MC_barrel_deno ->Add(h_ineff_MC_barrel);
    h_eff_MC_barrel2_deno->Add(h_ineff_MC_barrel2);
    h_eff_MC_endcap_deno ->Add(h_ineff_MC_endcap);
    h_eff_MC_endcap2_deno->Add(h_ineff_MC_endcap2);
    h_eff_MC_eta_deno    ->Add(h_ineff_MC_eta);
    h_eff_MC_MET_deno    ->Add(h_ineff_MC_MET);
    h_eff_MC_MT_deno     ->Add(h_ineff_MC_MT);
    h_eff_MC_nJet_deno   ->Add(h_ineff_MC_nJet);
    h_eff_MC_dRjet_deno  ->Add(h_ineff_MC_dRjet);
    h_eff_MC_dPhiMET_deno->Add(h_ineff_MC_dPhiMET);

    h_eff_data_barrel   ->Divide(h_eff_data_barrel_deno);
    h_eff_data_barrel2  ->Divide(h_eff_data_barrel2_deno);
    h_eff_data_endcap   ->Divide(h_eff_data_endcap_deno);
    h_eff_data_endcap2  ->Divide(h_eff_data_endcap2_deno);
    h_eff_data_eta      ->Divide(h_eff_data_eta_deno);
    h_eff_data_MET      ->Divide(h_eff_data_MET_deno);
    h_eff_data_MT       ->Divide(h_eff_data_MT_deno);
    h_eff_data_nJet     ->Divide(h_eff_data_nJet_deno);
    h_eff_data_dRjet    ->Divide(h_eff_data_dRjet_deno);
    h_eff_data_dPhiMET  ->Divide(h_eff_data_dPhiMET_deno);
    h_eff_MC_barrel     ->Divide(h_eff_MC_barrel_deno);
    h_eff_MC_barrel2    ->Divide(h_eff_MC_barrel2_deno);
    h_eff_MC_endcap     ->Divide(h_eff_MC_endcap_deno);
    h_eff_MC_endcap2    ->Divide(h_eff_MC_endcap2_deno);
    h_eff_MC_eta        ->Divide(h_eff_MC_eta_deno);
    h_eff_MC_MET        ->Divide(h_eff_MC_MET_deno);
    h_eff_MC_MT         ->Divide(h_eff_MC_MT_deno);
    h_eff_MC_nJet       ->Divide(h_eff_MC_nJet_deno);
    h_eff_MC_dRjet      ->Divide(h_eff_MC_dRjet_deno);
    h_eff_MC_dPhiMET    ->Divide(h_eff_MC_dPhiMET_deno);
    h_ineff_data_barrel ->Divide(h_eff_data_barrel_deno);
    h_ineff_data_barrel2->Divide(h_eff_data_barrel2_deno);
    h_ineff_data_endcap ->Divide(h_eff_data_endcap_deno);
    h_ineff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_ineff_data_eta    ->Divide(h_eff_data_eta_deno);
    h_ineff_data_MET    ->Divide(h_eff_data_MET_deno);
    h_ineff_data_MT     ->Divide(h_eff_data_MT_deno);
    h_ineff_data_nJet   ->Divide(h_eff_data_nJet_deno);
    h_ineff_data_dRjet  ->Divide(h_eff_data_dRjet_deno);
    h_ineff_data_dPhiMET->Divide(h_eff_data_dPhiMET_deno);
    h_ineff_MC_barrel   ->Divide(h_eff_MC_barrel_deno);
    h_ineff_MC_barrel2  ->Divide(h_eff_MC_barrel2_deno);
    h_ineff_MC_endcap   ->Divide(h_eff_MC_endcap_deno);
    h_ineff_MC_endcap2  ->Divide(h_eff_MC_endcap2_deno);
    h_ineff_MC_eta      ->Divide(h_eff_MC_eta_deno);
    h_ineff_MC_MET      ->Divide(h_eff_MC_MET_deno);
    h_ineff_MC_MT       ->Divide(h_eff_MC_MT_deno);
    h_ineff_MC_nJet     ->Divide(h_eff_MC_nJet_deno);
    h_ineff_MC_dRjet    ->Divide(h_eff_MC_dRjet_deno);
    h_ineff_MC_dPhiMET  ->Divide(h_eff_MC_dPhiMET_deno);

    // Ratio method
    TH1D *h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [1]->Clone("h_PR_ratio_barrel_deno")));
    TH1D *h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_ratio_barrel2_deno")));
    TH1D *h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [1]->Clone("h_PR_ratio_endcap_deno")));
    TH1D *h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_ratio_endcap2_deno")));
    TH1D *h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [1]->Clone("h_PR_ratio_eta_deno")));
    TH1D *h_PR_ratio_MET_deno     = ((TH1D*)(h_MET_pass       [1]->Clone("h_PR_ratio_MET_deno")));
    TH1D *h_PR_ratio_MT_deno      = ((TH1D*)(h_MT_pass        [1]->Clone("h_PR_ratio_MT_deno")));
    TH1D *h_PR_ratio_nJet_deno    = ((TH1D*)(h_nJet_pass      [1]->Clone("h_PR_ratio_nJet_deno")));
    TH1D *h_PR_ratio_dRjet_deno   = ((TH1D*)(h_dRjet_pass     [1]->Clone("h_PR_ratio_dRjet_deno")));
    TH1D *h_PR_ratio_dPhiMET_deno = ((TH1D*)(h_dPhiMET_pass   [1]->Clone("h_PR_ratio_dPhiMET_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_pass [2]);
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_pass [3]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_pass[2]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_pass[3]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_pass [2]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_pass [3]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_pass[2]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_pass[3]);
    h_PR_ratio_eta_deno    ->Add(h_eta_pass       [2]);
    h_PR_ratio_eta_deno    ->Add(h_eta_pass       [3]);
    h_PR_ratio_MET_deno    ->Add(h_MET_pass       [2]);
    h_PR_ratio_MET_deno    ->Add(h_MET_pass       [3]);
    h_PR_ratio_MT_deno     ->Add(h_MT_pass        [2]);
    h_PR_ratio_MT_deno     ->Add(h_MT_pass        [3]);
    h_PR_ratio_nJet_deno   ->Add(h_nJet_pass      [2]);
    h_PR_ratio_nJet_deno   ->Add(h_nJet_pass      [3]);
    h_PR_ratio_dRjet_deno  ->Add(h_dRjet_pass     [2]);
    h_PR_ratio_dRjet_deno  ->Add(h_dRjet_pass     [3]);
    h_PR_ratio_dPhiMET_deno->Add(h_dPhiMET_pass   [2]);
    h_PR_ratio_dPhiMET_deno->Add(h_dPhiMET_pass   [3]);

    TH1D *h_PR_ratio_barrel  = ((TH1D*)(h_pT_barrel_fail [1]->Clone("h_PR_ratio_barrel")));
    TH1D *h_PR_ratio_barrel2 = ((TH1D*)(h_pT_barrel2_fail[1]->Clone("h_PR_ratio_barrel2")));
    TH1D *h_PR_ratio_endcap  = ((TH1D*)(h_pT_endcap_fail [1]->Clone("h_PR_ratio_endcap")));
    TH1D *h_PR_ratio_endcap2 = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_PR_ratio_endcap2")));
    TH1D *h_PR_ratio_eta     = ((TH1D*)(h_eta_fail       [1]->Clone("h_PR_ratio_eta")));
    TH1D *h_PR_ratio_MET     = ((TH1D*)(h_MET_fail       [1]->Clone("h_PR_ratio_MET")));
    TH1D *h_PR_ratio_MT      = ((TH1D*)(h_MT_fail        [1]->Clone("h_PR_ratio_MT")));
    TH1D *h_PR_ratio_nJet    = ((TH1D*)(h_nJet_fail      [1]->Clone("h_PR_ratio_nJet")));
    TH1D *h_PR_ratio_dRjet   = ((TH1D*)(h_dRjet_fail     [1]->Clone("h_PR_ratio_dRjet")));
    TH1D *h_PR_ratio_dPhiMET = ((TH1D*)(h_dPhiMET_fail   [1]->Clone("h_PR_ratio_dPhiMET")));
    h_PR_ratio_barrel ->Add(h_pT_barrel_fail[2]);
    h_PR_ratio_barrel ->Add(h_pT_barrel_fail[3]);
    h_PR_ratio_barrel ->Add(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Add(h_pT_barrel2_fail[2]);
    h_PR_ratio_barrel2->Add(h_pT_barrel2_fail[3]);
    h_PR_ratio_barrel2->Add(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Add(h_pT_endcap_fail[2]);
    h_PR_ratio_endcap ->Add(h_pT_endcap_fail[3]);
    h_PR_ratio_endcap ->Add(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Add(h_pT_endcap2_fail[2]);
    h_PR_ratio_endcap2->Add(h_pT_endcap2_fail[3]);
    h_PR_ratio_endcap2->Add(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Add(h_eta_fail[2]);
    h_PR_ratio_eta    ->Add(h_eta_fail[3]);
    h_PR_ratio_eta    ->Add(h_PR_ratio_eta_deno);
    h_PR_ratio_MET    ->Add(h_MET_fail[2]);
    h_PR_ratio_MET    ->Add(h_MET_fail[3]);
    h_PR_ratio_MET    ->Add(h_PR_ratio_MET_deno);
    h_PR_ratio_MT     ->Add(h_MT_fail[2]);
    h_PR_ratio_MT     ->Add(h_MT_fail[3]);
    h_PR_ratio_MT     ->Add(h_PR_ratio_MT_deno);
    h_PR_ratio_nJet   ->Add(h_nJet_fail[2]);
    h_PR_ratio_nJet   ->Add(h_nJet_fail[3]);
    h_PR_ratio_nJet   ->Add(h_PR_ratio_nJet_deno);
    h_PR_ratio_dRjet  ->Add(h_dRjet_fail[2]);
    h_PR_ratio_dRjet  ->Add(h_dRjet_fail[3]);
    h_PR_ratio_dRjet  ->Add(h_PR_ratio_dRjet_deno);
    h_PR_ratio_dPhiMET->Add(h_dPhiMET_fail[2]);
    h_PR_ratio_dPhiMET->Add(h_dPhiMET_fail[3]);
    h_PR_ratio_dPhiMET->Add(h_PR_ratio_dPhiMET_deno);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);
    h_PR_ratio_MET    ->Divide(h_PR_ratio_MET_deno);
    h_PR_ratio_MT     ->Divide(h_PR_ratio_MT_deno);
    h_PR_ratio_nJet   ->Divide(h_PR_ratio_nJet_deno);
    h_PR_ratio_dRjet  ->Divide(h_PR_ratio_dRjet_deno);
    h_PR_ratio_dPhiMET->Divide(h_PR_ratio_dPhiMET_deno);

    h_PR_ratio_barrel ->Multiply(h_pT_barrel_pass [0]);
    h_PR_ratio_barrel2->Multiply(h_pT_barrel2_pass[0]);
    h_PR_ratio_endcap ->Multiply(h_pT_endcap_pass [0]);
    h_PR_ratio_endcap2->Multiply(h_pT_endcap2_pass[0]);
    h_PR_ratio_eta    ->Multiply(h_eta_pass       [0]);
    h_PR_ratio_MET    ->Multiply(h_MET_pass       [0]);
    h_PR_ratio_MT     ->Multiply(h_MT_pass        [0]);
    h_PR_ratio_nJet   ->Multiply(h_nJet_pass      [0]);
    h_PR_ratio_dRjet  ->Multiply(h_dRjet_pass     [0]);
    h_PR_ratio_dPhiMET->Multiply(h_dPhiMET_pass   [0]);

    h_PR_ratio_barrel ->Multiply(h_pT_barrel_pass [1]);
    h_PR_ratio_barrel2->Multiply(h_pT_barrel2_pass[1]);
    h_PR_ratio_endcap ->Multiply(h_pT_endcap_pass [1]);
    h_PR_ratio_endcap2->Multiply(h_pT_endcap2_pass[1]);
    h_PR_ratio_eta    ->Multiply(h_eta_pass       [1]);
    h_PR_ratio_MET    ->Multiply(h_MET_pass       [1]);
    h_PR_ratio_MT     ->Multiply(h_MT_pass        [1]);
    h_PR_ratio_nJet   ->Multiply(h_nJet_pass      [1]);
    h_PR_ratio_dRjet  ->Multiply(h_dRjet_pass     [1]);
    h_PR_ratio_dPhiMET->Multiply(h_dPhiMET_pass   [1]);

    h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [0]->Clone("h_PR_ratio_barrel_deno")));
    h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[0]->Clone("h_PR_ratio_barrel2_deno")));
    h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [0]->Clone("h_PR_ratio_endcap_deno")));
    h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_ratio_endcap2_deno")));
    h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [0]->Clone("h_PR_ratio_eta_deno")));
    h_PR_ratio_MET_deno     = ((TH1D*)(h_MET_pass       [0]->Clone("h_PR_ratio_MET_deno")));
    h_PR_ratio_MT_deno      = ((TH1D*)(h_MT_pass        [0]->Clone("h_PR_ratio_MT_deno")));
    h_PR_ratio_nJet_deno    = ((TH1D*)(h_nJet_pass      [0]->Clone("h_PR_ratio_nJet_deno")));
    h_PR_ratio_dRjet_deno   = ((TH1D*)(h_dRjet_pass     [0]->Clone("h_PR_ratio_dRjet_deno")));
    h_PR_ratio_dPhiMET_deno = ((TH1D*)(h_dPhiMET_pass   [0]->Clone("h_PR_ratio_dPhiMET_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_fail [0]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_fail[0]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_fail [0]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_fail[0]);
    h_PR_ratio_eta_deno    ->Add(h_eta_fail       [0]);
    h_PR_ratio_MET_deno    ->Add(h_MET_fail       [0]);
    h_PR_ratio_MT_deno     ->Add(h_MT_fail        [0]);
    h_PR_ratio_nJet_deno   ->Add(h_nJet_fail      [0]);
    h_PR_ratio_dRjet_deno  ->Add(h_dRjet_fail     [0]);
    h_PR_ratio_dPhiMET_deno->Add(h_dPhiMET_fail   [0]);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);
    h_PR_ratio_MET    ->Divide(h_PR_ratio_MET_deno);
    h_PR_ratio_MT     ->Divide(h_PR_ratio_MT_deno);
    h_PR_ratio_nJet   ->Divide(h_PR_ratio_nJet_deno);
    h_PR_ratio_dRjet  ->Divide(h_PR_ratio_dRjet_deno);
    h_PR_ratio_dPhiMET->Divide(h_PR_ratio_dPhiMET_deno);

    h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [1]->Clone("h_PR_ratio_barrel_deno")));
    h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_ratio_barrel2_deno")));
    h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [1]->Clone("h_PR_ratio_endcap_deno")));
    h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_ratio_endcap2_deno")));
    h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [1]->Clone("h_PR_ratio_eta_deno")));
    h_PR_ratio_MET_deno     = ((TH1D*)(h_MET_pass       [1]->Clone("h_PR_ratio_MET_deno")));
    h_PR_ratio_MT_deno      = ((TH1D*)(h_MT_pass        [1]->Clone("h_PR_ratio_MT_deno")));
    h_PR_ratio_nJet_deno    = ((TH1D*)(h_nJet_pass      [1]->Clone("h_PR_ratio_nJet_deno")));
    h_PR_ratio_dRjet_deno   = ((TH1D*)(h_dRjet_pass     [1]->Clone("h_PR_ratio_dRjet_deno")));
    h_PR_ratio_dPhiMET_deno = ((TH1D*)(h_dPhiMET_pass   [1]->Clone("h_PR_ratio_dPhiMET_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_fail [1]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_fail[1]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_fail [1]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_fail[1]);
    h_PR_ratio_eta_deno    ->Add(h_eta_fail       [1]);
    h_PR_ratio_MET_deno    ->Add(h_MET_fail       [1]);
    h_PR_ratio_MT_deno     ->Add(h_MT_fail        [1]);
    h_PR_ratio_nJet_deno   ->Add(h_nJet_fail      [1]);
    h_PR_ratio_dRjet_deno  ->Add(h_dRjet_fail     [1]);
    h_PR_ratio_dPhiMET_deno->Add(h_dPhiMET_fail   [1]);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);
    h_PR_ratio_MET    ->Divide(h_PR_ratio_MET_deno);
    h_PR_ratio_MT     ->Divide(h_PR_ratio_MT_deno);
    h_PR_ratio_nJet   ->Divide(h_PR_ratio_nJet_deno);
    h_PR_ratio_dRjet  ->Divide(h_PR_ratio_dRjet_deno);
    h_PR_ratio_dPhiMET->Divide(h_PR_ratio_dPhiMET_deno);


    // Drawing
    h_eff_data_barrel->SetDirectory(0);
    h_eff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel->SetMarkerColor(kBlack);
    h_eff_data_barrel->SetLineColor(kBlack);
    h_eff_data_barrel2->SetDirectory(0);
    h_eff_data_barrel2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel2->SetMarkerColor(kBlack);
    h_eff_data_barrel2->SetLineColor(kBlack);
    h_eff_data_endcap->SetDirectory(0);
    h_eff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap->SetMarkerColor(kBlack);
    h_eff_data_endcap->SetLineColor(kBlack);
    h_eff_data_endcap2->SetDirectory(0);
    h_eff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2->SetMarkerColor(kBlack);
    h_eff_data_endcap2->SetLineColor(kBlack);
    h_eff_data_eta->SetDirectory(0);
    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_data_eta->SetMarkerColor(kBlack);
    h_eff_data_eta->SetLineColor(kBlack);
    h_eff_data_MET->SetDirectory(0);
    h_eff_data_MET->SetMarkerStyle(kFullDotLarge);
    h_eff_data_MET->SetMarkerColor(kBlack);
    h_eff_data_MET->SetLineColor(kBlack);
    h_eff_data_MT->SetDirectory(0);
    h_eff_data_MT->SetMarkerStyle(kFullDotLarge);
    h_eff_data_MT->SetMarkerColor(kBlack);
    h_eff_data_MT->SetLineColor(kBlack);
    h_eff_data_nJet->SetDirectory(0);
    h_eff_data_nJet->SetMarkerStyle(kFullDotLarge);
    h_eff_data_nJet->SetMarkerColor(kBlack);
    h_eff_data_nJet->SetLineColor(kBlack);
    h_eff_data_dRjet->SetDirectory(0);
    h_eff_data_dRjet->SetMarkerStyle(kFullDotLarge);
    h_eff_data_dRjet->SetMarkerColor(kBlack);
    h_eff_data_dRjet->SetLineColor(kBlack);
    h_eff_data_dPhiMET->SetDirectory(0);
    h_eff_data_dPhiMET->SetMarkerStyle(kFullDotLarge);
    h_eff_data_dPhiMET->SetMarkerColor(kBlack);
    h_eff_data_dPhiMET->SetLineColor(kBlack);
    h_eff_MC_barrel->SetDirectory(0);
    h_eff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel->SetMarkerColor(kRed);
    h_eff_MC_barrel->SetLineColor(kRed);
    h_eff_MC_barrel2->SetDirectory(0);
    h_eff_MC_barrel2->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel2->SetMarkerColor(kRed);
    h_eff_MC_barrel2->SetLineColor(kRed);
    h_eff_MC_endcap->SetDirectory(0);
    h_eff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap->SetMarkerColor(kRed);
    h_eff_MC_endcap->SetLineColor(kRed);
    h_eff_MC_endcap2->SetDirectory(0);
    h_eff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap2->SetMarkerColor(kRed);
    h_eff_MC_endcap2->SetLineColor(kRed);
    h_eff_MC_eta->SetDirectory(0);
    h_eff_MC_eta->SetMarkerStyle(kFullSquare);
    h_eff_MC_eta->SetMarkerColor(kRed);
    h_eff_MC_eta->SetLineColor(kRed);
    h_eff_MC_MET->SetDirectory(0);
    h_eff_MC_MET->SetMarkerStyle(kFullSquare);
    h_eff_MC_MET->SetMarkerColor(kRed);
    h_eff_MC_MET->SetLineColor(kRed);
    h_eff_MC_MT->SetDirectory(0);
    h_eff_MC_MT->SetMarkerStyle(kFullSquare);
    h_eff_MC_MT->SetMarkerColor(kRed);
    h_eff_MC_MT->SetLineColor(kRed);
    h_eff_MC_nJet->SetDirectory(0);
    h_eff_MC_nJet->SetMarkerStyle(kFullSquare);
    h_eff_MC_nJet->SetMarkerColor(kRed);
    h_eff_MC_nJet->SetLineColor(kRed);
    h_eff_MC_dRjet->SetDirectory(0);
    h_eff_MC_dRjet->SetMarkerStyle(kFullSquare);
    h_eff_MC_dRjet->SetMarkerColor(kRed);
    h_eff_MC_dRjet->SetLineColor(kRed);
    h_eff_MC_dPhiMET->SetDirectory(0);
    h_eff_MC_dPhiMET->SetMarkerStyle(kFullSquare);
    h_eff_MC_dPhiMET->SetMarkerColor(kRed);
    h_eff_MC_dPhiMET->SetLineColor(kRed);
    h_PR_ratio_barrel->SetDirectory(0);
    h_PR_ratio_barrel->SetMarkerStyle(33);
    h_PR_ratio_barrel->SetMarkerSize(1.5);
    h_PR_ratio_barrel->SetMarkerColor(kBlue);
    h_PR_ratio_barrel->SetLineColor(kBlue);
    h_PR_ratio_barrel2->SetDirectory(0);
    h_PR_ratio_barrel2->SetMarkerStyle(33);
    h_PR_ratio_barrel2->SetMarkerSize(1.5);
    h_PR_ratio_barrel2->SetMarkerColor(kBlue);
    h_PR_ratio_barrel2->SetLineColor(kBlue);
    h_PR_ratio_endcap->SetDirectory(0);
    h_PR_ratio_endcap->SetMarkerStyle(33);
    h_PR_ratio_endcap->SetMarkerSize(1.5);
    h_PR_ratio_endcap->SetMarkerColor(kBlue);
    h_PR_ratio_endcap->SetLineColor(kBlue);
    h_PR_ratio_endcap2->SetDirectory(0);
    h_PR_ratio_endcap2->SetMarkerStyle(33);
    h_PR_ratio_endcap2->SetMarkerSize(1.5);
    h_PR_ratio_endcap2->SetMarkerColor(kBlue);
    h_PR_ratio_endcap2->SetLineColor(kBlue);
    h_PR_ratio_eta->SetDirectory(0);
    h_PR_ratio_eta->SetMarkerStyle(33);
    h_PR_ratio_eta->SetMarkerSize(1.5);
    h_PR_ratio_eta->SetMarkerColor(kBlue);
    h_PR_ratio_eta->SetLineColor(kBlue);
    h_PR_ratio_MET->SetDirectory(0);
    h_PR_ratio_MET->SetMarkerStyle(33);
    h_PR_ratio_MET->SetMarkerSize(1.5);
    h_PR_ratio_MET->SetMarkerColor(kBlue);
    h_PR_ratio_MET->SetLineColor(kBlue);
    h_PR_ratio_MT->SetDirectory(0);
    h_PR_ratio_MT->SetMarkerStyle(33);
    h_PR_ratio_MT->SetMarkerSize(1.5);
    h_PR_ratio_MT->SetMarkerColor(kBlue);
    h_PR_ratio_MT->SetLineColor(kBlue);
    h_PR_ratio_nJet->SetDirectory(0);
    h_PR_ratio_nJet->SetMarkerStyle(33);
    h_PR_ratio_nJet->SetMarkerSize(1.5);
    h_PR_ratio_nJet->SetMarkerColor(kBlue);
    h_PR_ratio_nJet->SetLineColor(kBlue);
    h_PR_ratio_dRjet->SetDirectory(0);
    h_PR_ratio_dRjet->SetMarkerStyle(33);
    h_PR_ratio_dRjet->SetMarkerSize(1.5);
    h_PR_ratio_dRjet->SetMarkerColor(kBlue);
    h_PR_ratio_dRjet->SetLineColor(kBlue);
    h_PR_ratio_dPhiMET->SetDirectory(0);
    h_PR_ratio_dPhiMET->SetMarkerStyle(33);
    h_PR_ratio_dPhiMET->SetMarkerSize(1.5);
    h_PR_ratio_dPhiMET->SetMarkerColor(kBlue);
    h_PR_ratio_dPhiMET->SetLineColor(kBlue);

    h_ineff_data_barrel->SetDirectory(0);
    h_ineff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel->SetMarkerColor(kBlack);
    h_ineff_data_barrel->SetLineColor(kBlack);
    h_ineff_data_barrel2->SetDirectory(0);
    h_ineff_data_barrel2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel2->SetMarkerColor(kBlack);
    h_ineff_data_barrel2->SetLineColor(kBlack);
    h_ineff_data_endcap->SetDirectory(0);
    h_ineff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap->SetMarkerColor(kBlack);
    h_ineff_data_endcap->SetLineColor(kBlack);
    h_ineff_data_endcap2->SetDirectory(0);
    h_ineff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap2->SetMarkerColor(kBlack);
    h_ineff_data_endcap2->SetLineColor(kBlack);
    h_ineff_data_eta->SetDirectory(0);
    h_ineff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_eta->SetMarkerColor(kBlack);
    h_ineff_data_eta->SetLineColor(kBlack);
    h_ineff_data_MET->SetDirectory(0);
    h_ineff_data_MET->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_MET->SetMarkerColor(kBlack);
    h_ineff_data_MET->SetLineColor(kBlack);
    h_ineff_data_MT->SetDirectory(0);
    h_ineff_data_MT->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_MT->SetMarkerColor(kBlack);
    h_ineff_data_MT->SetLineColor(kBlack);
    h_ineff_data_nJet->SetDirectory(0);
    h_ineff_data_nJet->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_nJet->SetMarkerColor(kBlack);
    h_ineff_data_nJet->SetLineColor(kBlack);
    h_ineff_data_dRjet->SetDirectory(0);
    h_ineff_data_dRjet->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_dRjet->SetMarkerColor(kBlack);
    h_ineff_data_dRjet->SetLineColor(kBlack);
    h_ineff_data_dPhiMET->SetDirectory(0);
    h_ineff_data_dPhiMET->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_dPhiMET->SetMarkerColor(kBlack);
    h_ineff_data_dPhiMET->SetLineColor(kBlack);
    h_ineff_MC_barrel->SetDirectory(0);
    h_ineff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel->SetMarkerColor(kRed);
    h_ineff_MC_barrel->SetLineColor(kRed);
    h_ineff_MC_barrel2->SetDirectory(0);
    h_ineff_MC_barrel2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel2->SetMarkerColor(kRed);
    h_ineff_MC_barrel2->SetLineColor(kRed);
    h_ineff_MC_endcap->SetDirectory(0);
    h_ineff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap->SetMarkerColor(kRed);
    h_ineff_MC_endcap->SetLineColor(kRed);
    h_ineff_MC_endcap2->SetDirectory(0);
    h_ineff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap2->SetMarkerColor(kRed);
    h_ineff_MC_endcap2->SetLineColor(kRed);
    h_ineff_MC_eta->SetDirectory(0);
    h_ineff_MC_eta->SetMarkerStyle(kFullSquare);
    h_ineff_MC_eta->SetMarkerColor(kRed);
    h_ineff_MC_eta->SetLineColor(kRed);
    h_ineff_MC_MET->SetDirectory(0);
    h_ineff_MC_MET->SetMarkerStyle(kFullSquare);
    h_ineff_MC_MET->SetMarkerColor(kRed);
    h_ineff_MC_MET->SetLineColor(kRed);
    h_ineff_MC_MT->SetDirectory(0);
    h_ineff_MC_MT->SetMarkerStyle(kFullSquare);
    h_ineff_MC_MT->SetMarkerColor(kRed);
    h_ineff_MC_MT->SetLineColor(kRed);
    h_ineff_MC_nJet->SetDirectory(0);
    h_ineff_MC_nJet->SetMarkerStyle(kFullSquare);
    h_ineff_MC_nJet->SetMarkerColor(kRed);
    h_ineff_MC_nJet->SetLineColor(kRed);
    h_ineff_MC_dRjet->SetDirectory(0);
    h_ineff_MC_dRjet->SetMarkerStyle(kFullSquare);
    h_ineff_MC_dRjet->SetMarkerColor(kRed);
    h_ineff_MC_dRjet->SetLineColor(kRed);
    h_ineff_MC_dPhiMET->SetDirectory(0);
    h_ineff_MC_dPhiMET->SetMarkerStyle(kFullSquare);
    h_ineff_MC_dPhiMET->SetMarkerColor(kRed);
    h_ineff_MC_dPhiMET->SetLineColor(kRed);

    TH1D *h_eff_ratio_barrel    = ((TH1D*)(h_eff_data_barrel   ->Clone("h_eff_ratio_barrel")));
    TH1D *h_eff_ratio_barrel2   = ((TH1D*)(h_eff_data_barrel2  ->Clone("h_eff_ratio_barrel2")));
    TH1D *h_eff_ratio_endcap    = ((TH1D*)(h_eff_data_endcap   ->Clone("h_eff_ratio_endcap")));
    TH1D *h_eff_ratio_endcap2   = ((TH1D*)(h_eff_data_endcap2  ->Clone("h_eff_ratio_endcap2")));
    TH1D *h_eff_ratio_eta       = ((TH1D*)(h_eff_data_eta      ->Clone("h_eff_ratio_eta")));
    TH1D *h_eff_ratio_MET       = ((TH1D*)(h_eff_data_MET      ->Clone("h_eff_ratio_MET")));
    TH1D *h_eff_ratio_MT        = ((TH1D*)(h_eff_data_MT       ->Clone("h_eff_ratio_MT")));
    TH1D *h_eff_ratio_nJet      = ((TH1D*)(h_eff_data_nJet     ->Clone("h_eff_ratio_nJet")));
    TH1D *h_eff_ratio_dRjet     = ((TH1D*)(h_eff_data_dRjet    ->Clone("h_eff_ratio_dRjet")));
    TH1D *h_eff_ratio_dPhiMET   = ((TH1D*)(h_eff_data_dPhiMET  ->Clone("h_eff_ratio_dPhiMET")));
    TH1D *h_ineff_ratio_barrel  = ((TH1D*)(h_ineff_data_barrel ->Clone("h_1-eff_ratio_barrel")));
    TH1D *h_ineff_ratio_barrel2 = ((TH1D*)(h_ineff_data_barrel2->Clone("h_1-eff_ratio_barrel2")));
    TH1D *h_ineff_ratio_endcap  = ((TH1D*)(h_ineff_data_endcap ->Clone("h_1-eff_ratio_endcap")));
    TH1D *h_ineff_ratio_endcap2 = ((TH1D*)(h_ineff_data_endcap2->Clone("h_1-eff_ratio_endcap2")));
    TH1D *h_ineff_ratio_eta     = ((TH1D*)(h_ineff_data_eta    ->Clone("h_1-eff_ratio_eta")));
    TH1D *h_ineff_ratio_MET     = ((TH1D*)(h_ineff_data_MET    ->Clone("h_1-eff_ratio_MET")));
    TH1D *h_ineff_ratio_MT      = ((TH1D*)(h_ineff_data_MT     ->Clone("h_1-eff_ratio_MT")));
    TH1D *h_ineff_ratio_nJet    = ((TH1D*)(h_ineff_data_nJet   ->Clone("h_1-eff_ratio_nJet")));
    TH1D *h_ineff_ratio_dRjet   = ((TH1D*)(h_ineff_data_dRjet  ->Clone("h_1-eff_ratio_dRjet")));
    TH1D *h_ineff_ratio_dPhiMET = ((TH1D*)(h_ineff_data_dPhiMET->Clone("h_1-eff_ratio_dPhiMET")));
    h_eff_ratio_barrel ->Divide(h_eff_MC_barrel);
    h_eff_ratio_barrel2->Divide(h_eff_MC_barrel);
    h_eff_ratio_endcap ->Divide(h_eff_MC_endcap);
    h_eff_ratio_endcap2->Divide(h_eff_MC_endcap2);
    h_eff_ratio_eta    ->Divide(h_eff_MC_eta);
    h_eff_ratio_MET    ->Divide(h_eff_MC_MET);
    h_eff_ratio_MT     ->Divide(h_eff_MC_MT);
    h_eff_ratio_nJet   ->Divide(h_eff_MC_nJet);
    h_eff_ratio_dRjet  ->Divide(h_eff_MC_dRjet);
    h_eff_ratio_dPhiMET->Divide(h_eff_MC_dPhiMET);
    h_eff_ratio_barrel ->SetDirectory(0);
    h_eff_ratio_barrel2->SetDirectory(0);
    h_eff_ratio_barrel2->SetMarkerColor(kBlue);
    h_eff_ratio_barrel2->SetLineColor(kBlue);
    h_eff_ratio_endcap ->SetDirectory(0);
    h_eff_ratio_endcap ->SetMarkerColor(kBlue);
    h_eff_ratio_endcap ->SetLineColor(kBlue);
    h_eff_ratio_endcap2->SetDirectory(0);
    h_eff_ratio_endcap2->SetMarkerColor(kBlue);
    h_eff_ratio_endcap2->SetLineColor(kBlue);
    h_eff_ratio_eta    ->SetDirectory(0);
    h_eff_ratio_eta    ->SetMarkerColor(kBlue);
    h_eff_ratio_eta    ->SetLineColor(kBlue);
    h_eff_ratio_MET    ->SetDirectory(0);
    h_eff_ratio_MET    ->SetMarkerColor(kBlue);
    h_eff_ratio_MET    ->SetLineColor(kBlue);
    h_eff_ratio_MT     ->SetDirectory(0);
    h_eff_ratio_MT     ->SetMarkerColor(kBlue);
    h_eff_ratio_MT     ->SetLineColor(kBlue);
    h_eff_ratio_nJet   ->SetDirectory(0);
    h_eff_ratio_nJet   ->SetMarkerColor(kBlue);
    h_eff_ratio_nJet   ->SetLineColor(kBlue);
    h_eff_ratio_dRjet  ->SetDirectory(0);
    h_eff_ratio_dRjet  ->SetMarkerColor(kBlue);
    h_eff_ratio_dRjet  ->SetLineColor(kBlue);
    h_eff_ratio_dPhiMET->SetDirectory(0);
    h_eff_ratio_dPhiMET->SetMarkerColor(kBlue);
    h_eff_ratio_dPhiMET->SetLineColor(kBlue);
    h_ineff_ratio_barrel ->Divide(h_ineff_MC_barrel);
    h_ineff_ratio_barrel2->Divide(h_ineff_MC_barrel2);
    h_ineff_ratio_endcap ->Divide(h_ineff_MC_endcap);
    h_ineff_ratio_endcap2->Divide(h_ineff_MC_endcap2);
    h_ineff_ratio_eta    ->Divide(h_ineff_MC_eta);
    h_ineff_ratio_MET    ->Divide(h_ineff_MC_MET);
    h_ineff_ratio_MT     ->Divide(h_ineff_MC_MT);
    h_ineff_ratio_nJet   ->Divide(h_ineff_MC_nJet);
    h_ineff_ratio_dRjet  ->Divide(h_ineff_MC_dRjet);
    h_ineff_ratio_dPhiMET->Divide(h_ineff_MC_dPhiMET);
    h_ineff_ratio_barrel ->SetDirectory(0);
    h_ineff_ratio_barrel2->SetDirectory(0);
    h_ineff_ratio_barrel2->SetMarkerColor(kBlue);
    h_ineff_ratio_barrel2->SetLineColor(kBlue);
    h_ineff_ratio_endcap ->SetDirectory(0);
    h_ineff_ratio_endcap ->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap ->SetLineColor(kBlue);
    h_ineff_ratio_endcap2->SetDirectory(0);
    h_ineff_ratio_endcap2->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap2->SetLineColor(kBlue);
    h_ineff_ratio_eta    ->SetDirectory(0);
    h_ineff_ratio_eta    ->SetMarkerColor(kBlue);
    h_ineff_ratio_eta    ->SetLineColor(kBlue);
    h_ineff_ratio_MET    ->SetDirectory(0);
    h_ineff_ratio_MET    ->SetMarkerColor(kBlue);
    h_ineff_ratio_MET    ->SetLineColor(kBlue);
    h_ineff_ratio_MT     ->SetDirectory(0);
    h_ineff_ratio_MT     ->SetMarkerColor(kBlue);
    h_ineff_ratio_MT     ->SetLineColor(kBlue);
    h_ineff_ratio_nJet   ->SetDirectory(0);
    h_ineff_ratio_nJet   ->SetMarkerColor(kBlue);
    h_ineff_ratio_nJet   ->SetLineColor(kBlue);
    h_ineff_ratio_dRjet  ->SetDirectory(0);
    h_ineff_ratio_dRjet  ->SetMarkerColor(kBlue);
    h_ineff_ratio_dRjet  ->SetLineColor(kBlue);
    h_ineff_ratio_dPhiMET->SetDirectory(0);
    h_ineff_ratio_dPhiMET->SetMarkerColor(kBlue);
    h_ineff_ratio_dPhiMET->SetLineColor(kBlue);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_eff_barrel    = new myRatioPlot_t("RP_PR_barrel",    h_eff_MC_barrel,    h_eff_data_barrel);
    myRatioPlot_t *RP_eff_barrel2   = new myRatioPlot_t("RP_PR_barrel2",   h_eff_MC_barrel2,   h_eff_data_barrel2);
    myRatioPlot_t *RP_eff_endcap    = new myRatioPlot_t("RP_PR_endcap",    h_eff_MC_endcap,    h_eff_data_endcap);
    myRatioPlot_t *RP_eff_endcap2   = new myRatioPlot_t("RP_PR_endcap2",   h_eff_MC_endcap2,   h_eff_data_endcap2);
    myRatioPlot_t *RP_eff_eta       = new myRatioPlot_t("RP_PR_eta",       h_eff_MC_eta,       h_eff_data_eta);
    myRatioPlot_t *RP_eff_MET       = new myRatioPlot_t("RP_PR_MET",       h_eff_MC_MET,       h_eff_data_MET);
    myRatioPlot_t *RP_eff_MT        = new myRatioPlot_t("RP_PR_MT",        h_eff_MC_MT,        h_eff_data_MT);
    myRatioPlot_t *RP_eff_nJet      = new myRatioPlot_t("RP_PR_nJet",      h_eff_MC_nJet,      h_eff_data_nJet);
    myRatioPlot_t *RP_eff_dRjet     = new myRatioPlot_t("RP_PR_dRjet",     h_eff_MC_dRjet,     h_eff_data_dRjet);
    myRatioPlot_t *RP_eff_dPhiMET   = new myRatioPlot_t("RP_PR_dPhiMET",   h_eff_MC_dPhiMET,   h_eff_data_dPhiMET);
    myRatioPlot_t *RP_ineff_barrel  = new myRatioPlot_t("RP_1-PR_barrel",  h_ineff_MC_barrel,  h_ineff_data_barrel);
    myRatioPlot_t *RP_ineff_barrel2 = new myRatioPlot_t("RP_1-PR_barrel2", h_ineff_MC_barrel2, h_ineff_data_barrel2);
    myRatioPlot_t *RP_ineff_endcap  = new myRatioPlot_t("RP_1-PR_endcap",  h_ineff_MC_endcap,  h_ineff_data_endcap);
    myRatioPlot_t *RP_ineff_endcap2 = new myRatioPlot_t("RP_1-PR_endcap2", h_ineff_MC_endcap2, h_ineff_data_endcap2);
    myRatioPlot_t *RP_ineff_eta     = new myRatioPlot_t("RP_1-PR_eta",     h_ineff_MC_eta,     h_ineff_data_eta);
    myRatioPlot_t *RP_ineff_MET     = new myRatioPlot_t("RP_1-PR_MET",     h_ineff_MC_MET,     h_ineff_data_MET);
    myRatioPlot_t *RP_ineff_MT      = new myRatioPlot_t("RP_1-PR_MT",      h_ineff_MC_MT,      h_ineff_data_MT);
    myRatioPlot_t *RP_ineff_nJet    = new myRatioPlot_t("RP_1-PR_nJet",    h_ineff_MC_nJet,    h_ineff_data_nJet);
    myRatioPlot_t *RP_ineff_dRjet   = new myRatioPlot_t("RP_1-PR_dRjet",   h_ineff_MC_dRjet,   h_ineff_data_dRjet);
    myRatioPlot_t *RP_ineff_dPhiMET = new myRatioPlot_t("RP_1-PR_dPhiMET", h_ineff_MC_dPhiMET, h_ineff_data_dPhiMET);

    RP_eff_barrel   ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_barrel2  ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap   ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap2  ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_eta      ->SetPlots("#eta", -2.4, 2.4);
    RP_eff_MET      ->SetPlots("E_{#lower[-0.2]{T}}^{miss} [GeV]", 0, 100);
    RP_eff_MT       ->SetPlots("m_{#lower[-0.2]{T}} [GeV/c^{2}]", 0, 200);
    RP_eff_nJet     ->SetPlots("#jets", 0-0.5, 30-0.5);
    RP_eff_dRjet    ->SetPlots("#Delta R_{l, j})", 0, 5);
    RP_eff_dPhiMET  ->SetPlots("#Delta #phi_{l, MET}", 0, 6.3);
    RP_ineff_barrel ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_barrel2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_eta    ->SetPlots("#eta", -2.4, 2.4);
    RP_ineff_MET    ->SetPlots("E_{#lower[-0.2]{T}}^{miss} [GeV]", 0, 100);
    RP_ineff_MT     ->SetPlots("m_{#lower[-0.2]{T}} [GeV/c^{2}]", 0, 200);
    RP_ineff_nJet   ->SetPlots("#jets", 0-0.5, 30-0.5);
    RP_ineff_dRjet  ->SetPlots("#Delta R_{l, j})", 0, 5);
    RP_ineff_dPhiMET->SetPlots("#Delta #phi_{l, MET}", 0, 6.3);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_eff_data_barrel, "Data", "pl");
    legend->AddEntry(h_eff_MC_barrel, "DY MC", "pl");

    RP_eff_barrel   ->ImportLegend(legend);
    RP_eff_barrel2  ->ImportLegend(legend);
    RP_eff_endcap   ->ImportLegend(legend);
    RP_eff_endcap2  ->ImportLegend(legend);
    RP_eff_eta      ->ImportLegend(legend);
    RP_eff_MET      ->ImportLegend(legend);
    RP_eff_MT       ->ImportLegend(legend);
    RP_eff_nJet     ->ImportLegend(legend);
    RP_eff_dRjet    ->ImportLegend(legend);
    RP_eff_dPhiMET  ->ImportLegend(legend);
    RP_ineff_barrel ->ImportLegend(legend);
    RP_ineff_barrel2->ImportLegend(legend);
    RP_ineff_endcap ->ImportLegend(legend);
    RP_ineff_endcap2->ImportLegend(legend);
    RP_ineff_eta    ->ImportLegend(legend);
    RP_ineff_MET    ->ImportLegend(legend);
    RP_ineff_MT     ->ImportLegend(legend);
    RP_ineff_nJet   ->ImportLegend(legend);
    RP_ineff_dRjet  ->ImportLegend(legend);
    RP_ineff_dPhiMET->ImportLegend(legend);

//    RP_eff_barrel   ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_barrel2  ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap   ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap2  ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_eta      ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_eff_MET      ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_eff_MT       ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_eff_nJet     ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_eff_dRjet    ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_eff_dPhiMET  ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_ineff_barrel ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_barrel2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_eta    ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_MET    ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_MT     ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_nJet   ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_dRjet  ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_dPhiMET->Draw(0.01, 1, 1, "", "1-PR");

//    RP_eff_barrel   ->pad1->SetLogy(0);
//    RP_eff_barrel2  ->pad1->SetLogy(0);
//    RP_eff_endcap   ->pad1->SetLogy(0);
//    RP_eff_endcap2  ->pad1->SetLogy(0);
//    RP_eff_eta      ->pad1->SetLogy(0);
//    RP_eff_MET      ->pad1->SetLogy(0);
//    RP_eff_MT       ->pad1->SetLogy(0);
//    RP_eff_nJet     ->pad1->SetLogy(0);
//    RP_eff_dRjet    ->pad1->SetLogy(0);
//    RP_eff_dPhiMET  ->pad1->SetLogy(0);
//    RP_ineff_barrel ->pad1->SetLogy(0);
//    RP_ineff_barrel2->pad1->SetLogy(0);
//    RP_ineff_endcap ->pad1->SetLogy(0);
//    RP_ineff_endcap2->pad1->SetLogy(0);
//    RP_ineff_eta    ->pad1->SetLogy(0);
//    RP_ineff_MET    ->pad1->SetLogy(0);
//    RP_ineff_MT     ->pad1->SetLogy(0);
//    RP_ineff_nJet   ->pad1->SetLogy(0);
//    RP_ineff_dRjet  ->pad1->SetLogy(0);
//    RP_ineff_dPhiMET->pad1->SetLogy(0);

    TCanvas *c_PR_barrel = new TCanvas("c_PR_barrel", "c_PR_barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
    TH1D *h_eff_data_barrel_draw = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_data_barrel_draw")));
    TH1D *h_eff_MC_barrel_draw   = ((TH1D*)(h_eff_MC_barrel  ->Clone("h_PR_MC_barrel_draw")));
    TH1D *h_PR_ratio_barrel_draw = ((TH1D*)(h_PR_ratio_barrel->Clone("h_PR_ratio_barrel_draw")));
    h_eff_data_barrel_draw->SetDirectory(0);
    h_eff_MC_barrel_draw->SetDirectory(0);
    h_PR_ratio_barrel_draw->SetDirectory(0);
    h_eff_data_barrel_draw->SetStats(kFALSE);
    h_eff_data_barrel_draw->SetTitle("");
    h_eff_data_barrel_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_barrel_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_barrel_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_barrel_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_barrel_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_barrel_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_barrel_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_barrel_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_barrel_draw->Draw();
    h_eff_MC_barrel_draw->Draw("same");
    h_PR_ratio_barrel_draw->Draw("same");
    TLegend *legend1 = new TLegend(0.6, 0.75, 0.95, 0.95);
    legend1->AddEntry(h_eff_data_barrel_draw, "Subtraction", "lp");
    legend1->AddEntry(h_PR_ratio_barrel_draw, "Ratio", "lp");
    legend1->AddEntry(h_eff_MC_barrel_draw, "MC", "lp");
    legend1->Draw();
    c_PR_barrel->Update();


    TCanvas *c_PR_barrel2 = new TCanvas("c_PR_barrel2", "c_PR_barrel2", 800, 800);
    c_PR_barrel2->cd();
    c_PR_barrel2->SetGrid(1);
    c_PR_barrel2->SetLogx(1);
    c_PR_barrel2->SetRightMargin(0.05);
    c_PR_barrel2->SetTopMargin(0.05);
    c_PR_barrel2->SetBottomMargin(0.12);
    c_PR_barrel2->SetLeftMargin(0.13);
    TH1D *h_eff_data_barrel2_draw = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_data_barrel2_draw")));
    TH1D *h_eff_MC_barrel2_draw   = ((TH1D*)(h_eff_MC_barrel2  ->Clone("h_PR_MC_barrel2_draw")));
    TH1D *h_PR_ratio_barrel2_draw = ((TH1D*)(h_PR_ratio_barrel2->Clone("h_PR_ratio_barrel2_draw")));
    h_eff_data_barrel2_draw->SetDirectory(0);
    h_eff_MC_barrel2_draw->SetDirectory(0);
    h_PR_ratio_barrel2_draw->SetDirectory(0);
    h_eff_data_barrel2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel2_draw->SetMarkerColor(kBlack);
    h_eff_data_barrel2_draw->SetLineColor(kBlack);
    h_eff_data_barrel2_draw->SetStats(kFALSE);
    h_eff_data_barrel2_draw->SetTitle("");
    h_eff_data_barrel2_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_barrel2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_barrel2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_barrel2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_barrel2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_barrel2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_barrel2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_barrel2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_barrel2_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_barrel2_draw->Draw();
    h_eff_MC_barrel2_draw->Draw("same");
    h_PR_ratio_barrel2_draw->Draw("same");
    legend1->Draw();
    c_PR_barrel2->Update();


    TCanvas *c_PR_endcap = new TCanvas("c_PR_endcap", "c_PR_endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap_draw = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_data_endcap_draw")));
    TH1D *h_eff_MC_endcap_draw   = ((TH1D*)(h_eff_MC_endcap  ->Clone("h_PR_MC_endcap_draw")));
    TH1D *h_PR_ratio_endcap_draw = ((TH1D*)(h_PR_ratio_endcap->Clone("h_PR_ratio_endcap_draw")));
    h_eff_data_endcap_draw->SetDirectory(0);
    h_eff_MC_endcap_draw->SetDirectory(0);
    h_PR_ratio_endcap_draw->SetDirectory(0);
    h_eff_data_endcap_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap_draw->SetLineColor(kBlack);
    h_eff_data_endcap_draw->SetStats(kFALSE);
    h_eff_data_endcap_draw->SetTitle("");
    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_endcap_draw->Draw();
    h_eff_MC_endcap_draw->Draw("same");
    h_PR_ratio_endcap_draw->Draw("same");
    legend1->Draw();
    c_PR_endcap->Update();


    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "c_PR_endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap2_draw = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_data_endcap2_draw")));
    TH1D *h_eff_MC_endcap2_draw   = ((TH1D*)(h_eff_MC_endcap2  ->Clone("h_PR_MC_endcap2_draw")));
    TH1D *h_PR_ratio_endcap2_draw = ((TH1D*)(h_PR_ratio_endcap2->Clone("h_PR_ratio_endcap2_draw")));
    h_eff_data_endcap2_draw->SetDirectory(0);
    h_eff_MC_endcap2_draw->SetDirectory(0);
    h_PR_ratio_endcap2_draw->SetDirectory(0);
    h_eff_data_endcap2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap2_draw->SetLineColor(kBlack);
    h_eff_data_endcap2_draw->SetStats(kFALSE);
    h_eff_data_endcap2_draw->SetTitle("");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap2_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_endcap2_draw->Draw();
    h_eff_MC_endcap2_draw->Draw("same");
    h_PR_ratio_endcap2_draw->Draw("same");
    legend1->Draw();
    c_PR_endcap2->Update();


    TCanvas *c_PR_allin1 = new TCanvas("c_PR_allin1", "c_PR_allin1", 800, 800);
    c_PR_allin1->cd();
    c_PR_allin1->SetGrid(1);
    c_PR_allin1->SetRightMargin(0.05);
    c_PR_allin1->SetTopMargin(0.05);
    c_PR_allin1->SetBottomMargin(0.12);
    c_PR_allin1->SetLeftMargin(0.13);
    h_PR_ratio_endcap_draw->SetDirectory(0);
    h_PR_ratio_endcap_draw->SetTitle("");
    h_PR_ratio_endcap_draw->SetMarkerStyle(kFullSquare);
    h_PR_ratio_endcap_draw->SetMarkerColor(kBlue);
    h_PR_ratio_endcap_draw->SetLineColor(kBlue);
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitleOffset(1.12);
    h_PR_ratio_endcap_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_PR_ratio_endcap_draw->GetXaxis()->SetRangeUser(0, 5000);
    h_PR_ratio_endcap_draw->GetXaxis()->SetNoExponent();
    h_PR_ratio_endcap_draw->GetXaxis()->SetMoreLogLabels();
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_PR_ratio_endcap2_draw->SetDirectory(0);
    h_PR_ratio_endcap2_draw->SetMarkerStyle(33);
    h_PR_ratio_endcap2_draw->SetMarkerSize(1.5);
    h_PR_ratio_endcap2_draw->SetMarkerColor(kOrange-3);
    h_PR_ratio_endcap2_draw->SetLineColor(kOrange-3);
    h_PR_ratio_barrel2_draw->SetDirectory(0);
    h_PR_ratio_barrel2_draw->SetMarkerStyle(23);
    h_PR_ratio_barrel2_draw->SetMarkerColor(kRed-3);
    h_PR_ratio_barrel2_draw->SetLineColor(kRed-3);
    h_PR_ratio_endcap_draw->Draw();
    h_PR_ratio_endcap2_draw->Draw("same");
    h_PR_ratio_barrel->Draw("same");
    h_PR_ratio_barrel2_draw->Draw("same");
    TLegend *legend_allin1 = new TLegend(0.5, 0.72, 0.95, 0.95);
    legend_allin1->AddEntry(h_PR_ratio_barrel, "Barrel |#eta| < 0.7", "LP");
    legend_allin1->AddEntry(h_PR_ratio_barrel2_draw, "Barrel 0.7 #leq |#eta| < 1.2", "LP");
    legend_allin1->AddEntry(h_PR_ratio_endcap_draw, "Endcap 1.2 #leq |#eta| < 1.8", "LP");
    legend_allin1->AddEntry(h_PR_ratio_endcap2_draw, "Endcap 1.8 #leq |#eta| < 2.4", "LP");
    legend_allin1->Draw();\
    c_PR_allin1->SetLogx();
    c_PR_allin1->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "c_PR_eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
    TH1D *h_eff_MC_eta_draw   = ((TH1D*)(h_eff_MC_eta  ->Clone("h_PR_MC_eta_draw")));
    TH1D *h_PR_ratio_eta_draw = ((TH1D*)(h_PR_ratio_eta->Clone("h_PR_ratio_eta_draw")));
    h_eff_data_eta->SetDirectory(0);
    h_eff_MC_eta_draw->SetDirectory(0);
    h_PR_ratio_eta_draw->SetDirectory(0);
    h_eff_data_eta->SetTitle("");
    h_eff_data_eta->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_eta->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_eta->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_eta->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_eff_data_eta->GetXaxis()->SetNoExponent();
    h_eff_data_eta->GetXaxis()->SetMoreLogLabels();
    h_eff_data_eta->GetXaxis()->SetTitle("#eta");
    h_eff_data_eta->GetXaxis()->SetTitleOffset(1);
    h_eff_data_eta->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_eta->Draw();
    h_eff_MC_eta_draw->Draw("same");
    h_PR_ratio_eta_draw->Draw("same");
//    TLegend *legend_eta = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend_eta->AddEntry(h_eff_data_eta, "Subtraction", "LP");
//    legend_eta->Draw();
    legend1->Draw();
    c_PR_eta->Update();

    TCanvas *c_PR_MET = new TCanvas("c_PR_MET", "c_PR_MET", 800, 800);
    c_PR_MET->cd();
    c_PR_MET->SetGrid(1);
    c_PR_MET->SetRightMargin(0.05);
    c_PR_MET->SetTopMargin(0.05);
    c_PR_MET->SetBottomMargin(0.12);
    c_PR_MET->SetLeftMargin(0.13);
    TH1D *h_eff_MC_MET_draw   = ((TH1D*)(h_eff_MC_MET  ->Clone("h_PR_MC_MET_draw")));
    TH1D *h_PR_ratio_MET_draw = ((TH1D*)(h_PR_ratio_MET->Clone("h_PR_ratio_MET_draw")));
    h_eff_data_MET->SetDirectory(0);
    h_eff_MC_MET_draw->SetDirectory(0);
    h_PR_ratio_MET_draw->SetDirectory(0);
    h_eff_data_MET->SetTitle("");
    h_eff_data_MET->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_MET->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_MET->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_MET->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_MET->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_MET->GetXaxis()->SetRangeUser(0, 100);
    h_eff_data_MET->GetXaxis()->SetNoExponent();
    h_eff_data_MET->GetXaxis()->SetMoreLogLabels();
    h_eff_data_MET->GetXaxis()->SetTitle("MET");
    h_eff_data_MET->GetXaxis()->SetTitleOffset(1);
    h_eff_data_MET->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_MET->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_MET->Draw();
    h_eff_MC_MET_draw->Draw("same");
    h_PR_ratio_MET_draw->Draw("same");
    legend1->Draw();
    c_PR_MET->Update();

    TCanvas *c_PR_MT = new TCanvas("c_PR_MT", "c_PR_MT", 800, 800);
    c_PR_MT->cd();
    c_PR_MT->SetGrid(1);
    c_PR_MT->SetRightMargin(0.05);
    c_PR_MT->SetTopMargin(0.05);
    c_PR_MT->SetBottomMargin(0.12);
    c_PR_MT->SetLeftMargin(0.13);
    TH1D *h_eff_MC_MT_draw   = ((TH1D*)(h_eff_MC_MT  ->Clone("h_PR_MC_MT_draw")));
    TH1D *h_PR_ratio_MT_draw = ((TH1D*)(h_PR_ratio_MT->Clone("h_PR_ratio_MT_draw")));
    h_eff_data_MT->SetDirectory(0);
    h_eff_MC_MT_draw->SetDirectory(0);
    h_PR_ratio_MT_draw->SetDirectory(0);
    h_eff_data_MT->SetTitle("");
    h_eff_data_MT->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_MT->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_MT->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_MT->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_MT->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_MT->GetXaxis()->SetRangeUser(0, 200);
    h_eff_data_MT->GetXaxis()->SetNoExponent();
    h_eff_data_MT->GetXaxis()->SetMoreLogLabels();
    h_eff_data_MT->GetXaxis()->SetTitle("MT");
    h_eff_data_MT->GetXaxis()->SetTitleOffset(1);
    h_eff_data_MT->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_MT->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_MT->Draw();
    h_eff_MC_MT_draw->Draw("same");
    h_PR_ratio_MT_draw->Draw("same");
    legend1->Draw();
    c_PR_MT->Update();

    TCanvas *c_PR_nJet = new TCanvas("c_PR_nJet", "c_PR_nJet", 800, 800);
    c_PR_nJet->cd();
    c_PR_nJet->SetGrid(1);
    c_PR_nJet->SetRightMargin(0.05);
    c_PR_nJet->SetTopMargin(0.05);
    c_PR_nJet->SetBottomMargin(0.12);
    c_PR_nJet->SetLeftMargin(0.13);
    TH1D *h_eff_MC_nJet_draw   = ((TH1D*)(h_eff_MC_nJet  ->Clone("h_PR_MC_nJet_draw")));
    TH1D *h_PR_ratio_nJet_draw = ((TH1D*)(h_PR_ratio_nJet->Clone("h_PR_ratio_nJet_draw")));
    h_eff_data_nJet->SetDirectory(0);
    h_eff_MC_nJet_draw->SetDirectory(0);
    h_PR_ratio_nJet_draw->SetDirectory(0);
    h_eff_data_nJet->SetTitle("");
    h_eff_data_nJet->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_nJet->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_nJet->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_nJet->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_nJet->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_nJet->GetXaxis()->SetRangeUser(0-0.5, 30-0.5);
    h_eff_data_nJet->GetXaxis()->SetNoExponent();
    h_eff_data_nJet->GetXaxis()->SetMoreLogLabels();
    h_eff_data_nJet->GetXaxis()->SetTitle("nJet");
    h_eff_data_nJet->GetXaxis()->SetTitleOffset(1);
    h_eff_data_nJet->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_nJet->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_nJet->Draw();
    h_eff_MC_nJet_draw->Draw("same");
    h_PR_ratio_nJet_draw->Draw("same");
    legend1->Draw();
    c_PR_nJet->Update();

    TCanvas *c_PR_dRjet = new TCanvas("c_PR_dRjet", "c_PR_dRjet", 800, 800);
    c_PR_dRjet->cd();
    c_PR_dRjet->SetGrid(1);
    c_PR_dRjet->SetRightMargin(0.05);
    c_PR_dRjet->SetTopMargin(0.05);
    c_PR_dRjet->SetBottomMargin(0.12);
    c_PR_dRjet->SetLeftMargin(0.13);
    TH1D *h_eff_MC_dRjet_draw   = ((TH1D*)(h_eff_MC_dRjet  ->Clone("h_PR_MC_dRjet_draw")));
    TH1D *h_PR_ratio_dRjet_draw = ((TH1D*)(h_PR_ratio_dRjet->Clone("h_PR_ratio_dRjet_draw")));
    h_eff_data_dRjet->SetDirectory(0);
    h_eff_MC_dRjet_draw->SetDirectory(0);
    h_PR_ratio_dRjet_draw->SetDirectory(0);
    h_eff_data_dRjet->SetTitle("");
    h_eff_data_dRjet->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_dRjet->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_dRjet->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_dRjet->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_dRjet->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_dRjet->GetXaxis()->SetRangeUser(0, 5);
    h_eff_data_dRjet->GetXaxis()->SetNoExponent();
    h_eff_data_dRjet->GetXaxis()->SetMoreLogLabels();
    h_eff_data_dRjet->GetXaxis()->SetTitle("dRjet");
    h_eff_data_dRjet->GetXaxis()->SetTitleOffset(1);
    h_eff_data_dRjet->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_dRjet->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_dRjet->Draw();
    h_eff_MC_dRjet_draw->Draw("same");
    h_PR_ratio_dRjet_draw->Draw("same");
    legend1->Draw();
    c_PR_dRjet->Update();

    TCanvas *c_PR_dPhiMET = new TCanvas("c_PR_dPhiMET", "c_PR_dPhiMET", 800, 800);
    c_PR_dPhiMET->cd();
    c_PR_dPhiMET->SetGrid(1);
    c_PR_dPhiMET->SetRightMargin(0.05);
    c_PR_dPhiMET->SetTopMargin(0.05);
    c_PR_dPhiMET->SetBottomMargin(0.12);
    c_PR_dPhiMET->SetLeftMargin(0.13);
    TH1D *h_eff_MC_dPhiMET_draw   = ((TH1D*)(h_eff_MC_dPhiMET  ->Clone("h_PR_MC_dPhiMET_draw")));
    TH1D *h_PR_ratio_dPhiMET_draw = ((TH1D*)(h_PR_ratio_dPhiMET->Clone("h_PR_ratio_dPhiMET_draw")));
    h_eff_data_dPhiMET->SetDirectory(0);
    h_eff_MC_dPhiMET_draw->SetDirectory(0);
    h_PR_ratio_dPhiMET_draw->SetDirectory(0);
    h_eff_data_dPhiMET->SetTitle("");
    h_eff_data_dPhiMET->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_dPhiMET->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_dPhiMET->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_dPhiMET->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_dPhiMET->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_dPhiMET->GetXaxis()->SetRangeUser(0, 6.3);
    h_eff_data_dPhiMET->GetXaxis()->SetNoExponent();
    h_eff_data_dPhiMET->GetXaxis()->SetMoreLogLabels();
    h_eff_data_dPhiMET->GetXaxis()->SetTitle("dPhiMET");
    h_eff_data_dPhiMET->GetXaxis()->SetTitleOffset(1);
    h_eff_data_dPhiMET->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_dPhiMET->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_dPhiMET->Draw();
    h_eff_MC_dPhiMET_draw->Draw("same");
    h_PR_ratio_dPhiMET_draw->Draw("same");
    legend1->Draw();
    c_PR_dPhiMET->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << inName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << inName << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/FR/Muon/PromptRate_muon.root", "RECREATE");
    f_out->cd();
    h_eff_data_barrel    ->Write();
    h_eff_data_barrel2   ->Write();
    h_eff_data_endcap    ->Write();
    h_eff_data_endcap2   ->Write();
    h_eff_data_eta       ->Write();
    h_eff_MC_barrel      ->Write();
    h_eff_MC_barrel2     ->Write();
    h_eff_MC_endcap      ->Write();
    h_eff_MC_endcap2     ->Write();
    h_eff_MC_eta         ->Write();
//    h_ineff_data_barrel  ->Write();
//    h_ineff_data_barrel2 ->Write();
//    h_ineff_data_endcap  ->Write();
//    h_ineff_data_endcap2 ->Write();
//    h_ineff_data_eta     ->Write();
//    h_ineff_MC_barrel    ->Write();
//    h_ineff_MC_barrel2   ->Write();
//    h_ineff_MC_endcap    ->Write();
//    h_ineff_MC_endcap2   ->Write();
//    h_ineff_MC_eta       ->Write();
//    h_eff_ratio_barrel   ->Write();
//    h_eff_ratio_barrel2  ->Write();
//    h_eff_ratio_endcap   ->Write();
//    h_eff_ratio_endcap2  ->Write();
//    h_eff_ratio_eta      ->Write();
//    h_ineff_ratio_barrel ->Write();
//    h_ineff_ratio_barrel2->Write();
//    h_ineff_ratio_endcap ->Write();
//    h_ineff_ratio_endcap2->Write();
//    h_ineff_ratio_eta    ->Write();
    h_PR_ratio_barrel    ->Write();
    h_PR_ratio_barrel2   ->Write();
    h_PR_ratio_endcap    ->Write();
    h_PR_ratio_endcap2   ->Write();
    h_PR_ratio_eta       ->Write();
//    f_barrel_17to43->Write();
//    f_barrel_43to56->Write();
//    f_barrel_56to95->Write();
//    f_barrel_95to5000->Write();
//    f_barrel2_17to43->Write();
//    f_barrel2_43to56->Write();
//    f_barrel2_56to95->Write();
//    f_barrel2_95to5000->Write();
//    f_endcap_17to42p5->Write();
//    f_endcap_42p5to55->Write();
//    f_endcap_55to5000->Write();
//    f_endcap2_17to42p5->Write();
//    f_endcap2_42p5to56->Write();
//    f_endcap2_56to76->Write();
//    f_endcap2_76to5000->Write();
    f_out->Close();

} // End of Mu_EstimatePR


void Mu_EstPR_alt()
{
    TString inName = "/media/sf_DATA/FR/Muon/PR_Hist_Mu_alt.root";
    TFile *f = new TFile(inName, "READ");

    TH1D *h_pT_barrel_pass[4],
         *h_pT_barrel2_pass[4],
         *h_pT_endcap_pass[4],
         *h_pT_endcap2_pass[4],
         *h_eta_pass[4],
         *h_pT_barrel_fail[4],
         *h_pT_barrel2_fail[4],
         *h_pT_endcap_fail[4],
         *h_pT_endcap2_fail[4],
         *h_eta_fail[4];
    TString type[4] = {"data", "wjet", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_barrel2_pass_"+type[i], h_pT_barrel2_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_endcap2_pass_"+type[i], h_pT_endcap2_pass[i]);
        f->GetObject("h_eta_pass_"+type[i], h_eta_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_barrel2_fail_"+type[i], h_pT_barrel2_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        f->GetObject("h_pT_endcap2_fail_"+type[i], h_pT_endcap2_fail[i]);
        f->GetObject("h_eta_fail_"+type[i], h_eta_fail[i]);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_barrel2_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_endcap2_pass[i]->SetDirectory(0);
        h_eta_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_barrel2_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
        h_pT_endcap2_fail[i]->SetDirectory(0);
        h_eta_fail[i]->SetDirectory(0);


        // FOR TEST: AVERAGING
        h_pT_barrel_pass[i]->Add(h_pT_barrel2_pass[i]);
        h_pT_endcap_pass[i]->Add(h_pT_endcap2_pass[i]);
        h_pT_barrel_fail[i]->Add(h_pT_barrel2_fail[i]);
        h_pT_endcap_fail[i]->Add(h_pT_endcap2_fail[i]);
    }

    // Subtraction and MC methods
    TH1D *h_eff_data_barrel  = ((TH1D*)(h_pT_barrel_pass[0]->Clone("h_PR_subtract_barrel")));
    TH1D *h_eff_data_barrel2 = ((TH1D*)(h_pT_barrel2_pass[0]->Clone("h_PR_subtract_barrel2")));
    TH1D *h_eff_data_endcap  = ((TH1D*)(h_pT_endcap_pass[0]->Clone("h_PR_subtract_endcap")));
    TH1D *h_eff_data_endcap2 = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_subtract_endcap2")));
    TH1D *h_eff_data_eta     = ((TH1D*)(h_eta_pass[0]->Clone("h_PR_subtract_eta")));
    TH1D *h_eff_MC_barrel    = ((TH1D*)(h_pT_barrel_pass[1]->Clone("h_PR_MC_barrel")));
    TH1D *h_eff_MC_barrel2   = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_MC_barrel2")));
    TH1D *h_eff_MC_endcap    = ((TH1D*)(h_pT_endcap_pass[1]->Clone("h_PR_MC_endcap")));
    TH1D *h_eff_MC_endcap2   = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_MC_endcap2")));
    TH1D *h_eff_MC_eta       = ((TH1D*)(h_eta_pass[1]->Clone("h_PR_MC_eta")));
    TH1D *h_eff_bkgr_barrel  = ((TH1D*)(h_pT_barrel_pass[2]->Clone("h_PR_bkgr_barrel")));
    TH1D *h_eff_bkgr_barrel2 = ((TH1D*)(h_pT_barrel2_pass[2]->Clone("h_PR_bkgr_barrel2")));
    TH1D *h_eff_bkgr_endcap  = ((TH1D*)(h_pT_endcap_pass[2]->Clone("h_PR_bkgr_endcap")));
    TH1D *h_eff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_pass[2]->Clone("h_PR_bkgr_endcap2")));
    TH1D *h_eff_bkgr_eta     = ((TH1D*)(h_eta_pass[2]->Clone("h_PR_bkgr_eta")));
    TH1D *h_eff_bkgf_barrel  = ((TH1D*)(h_pT_barrel_pass[3]->Clone("h_PR_bkgf_barrel")));
    TH1D *h_eff_bkgf_barrel2 = ((TH1D*)(h_pT_barrel2_pass[3]->Clone("h_PR_bkgf_barrel2")));
    TH1D *h_eff_bkgf_endcap  = ((TH1D*)(h_pT_endcap_pass[3]->Clone("h_PR_bkgf_endcap")));
    TH1D *h_eff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_pass[3]->Clone("h_PR_bkgf_endcap2")));
    TH1D *h_eff_bkgf_eta     = ((TH1D*)(h_eta_pass[3]->Clone("h_PR_bkgf_eta")));
    h_eff_data_barrel ->SetDirectory(0);
    h_eff_data_barrel2->SetDirectory(0);
    h_eff_data_endcap ->SetDirectory(0);
    h_eff_data_endcap2->SetDirectory(0);
    h_eff_data_eta    ->SetDirectory(0);
    h_eff_MC_barrel   ->SetDirectory(0);
    h_eff_MC_barrel2  ->SetDirectory(0);
    h_eff_MC_endcap   ->SetDirectory(0);
    h_eff_MC_endcap2  ->SetDirectory(0);
    h_eff_MC_eta      ->SetDirectory(0);
    h_eff_bkgr_barrel ->SetDirectory(0);
    h_eff_bkgr_barrel2->SetDirectory(0);
    h_eff_bkgr_endcap ->SetDirectory(0);
    h_eff_bkgr_endcap2->SetDirectory(0);
    h_eff_bkgr_eta    ->SetDirectory(0);
    h_eff_bkgf_barrel ->SetDirectory(0);
    h_eff_bkgf_barrel2->SetDirectory(0);
    h_eff_bkgf_endcap ->SetDirectory(0);
    h_eff_bkgf_endcap2->SetDirectory(0);
    h_eff_bkgf_eta    ->SetDirectory(0);
//    h_eff_data_barrel->Add(h_eff_bkgf_barrel, -1);
//    h_eff_data_barrel2->Add(h_eff_bkgf_barrel2, -1);
//    h_eff_data_endcap->Add(h_eff_bkgf_endcap, -1);
//    h_eff_data_endcap2->Add(h_eff_bkgf_endcap2, -1);
//    h_eff_data_eta->Add(h_eff_bkgf_eta, -1);
    h_eff_data_barrel->Add(h_eff_bkgr_barrel, -1);
    h_eff_data_barrel2->Add(h_eff_bkgr_barrel2, -1);
    h_eff_data_endcap->Add(h_eff_bkgr_endcap, -1);
    h_eff_data_endcap2->Add(h_eff_bkgr_endcap2, -1);
    h_eff_data_eta->Add(h_eff_bkgr_eta, -1);
    // FOR ERRORS
    Double_t val_signal_barrel[nPtBinEndcap], val_signal_barrel2[nPtBinEndcap], val_signal_endcap[nPtBinEndcap], val_signal_endcap2[nPtBinEndcap];
    Double_t err_signal_barrel[nPtBinEndcap], err_signal_barrel2[nPtBinEndcap], err_signal_endcap[nPtBinEndcap], err_signal_endcap2[nPtBinEndcap];
    cout << "SIGNAL" << endl;
    for (Int_t i_bin=1; i_bin<=nPtBinEndcap; i_bin++)
    {
        val_signal_barrel[i_bin-1] = h_eff_data_barrel->GetBinContent(i_bin);
        err_signal_barrel[i_bin-1] = h_eff_data_barrel->GetBinError(i_bin);
        val_signal_barrel2[i_bin-1] = h_eff_data_barrel2->GetBinContent(i_bin);
        err_signal_barrel2[i_bin-1] = h_eff_data_barrel2->GetBinError(i_bin);
        val_signal_endcap[i_bin-1] = h_eff_data_endcap->GetBinContent(i_bin);
        err_signal_endcap[i_bin-1] = h_eff_data_endcap->GetBinError(i_bin);
        val_signal_endcap2[i_bin-1] = h_eff_data_endcap2->GetBinContent(i_bin);
        err_signal_endcap2[i_bin-1] = h_eff_data_endcap2->GetBinError(i_bin);
        cout << "Bin" << i_bin << ":\t" << val_signal_barrel[i_bin-1] << "+-" << err_signal_barrel[i_bin-1] << "(b)\t";
        cout << val_signal_barrel2[i_bin-1] << "+-" << err_signal_barrel2[i_bin-1] << "(b2)\t";
        cout << val_signal_endcap[i_bin-1] << "+-" << err_signal_endcap[i_bin-1] << "(e)\t";
        cout << val_signal_endcap2[i_bin-1] << "+-" << err_signal_endcap2[i_bin-1] << "(e2)\t" << endl;
    }
//    removeNegativeBins(h_eff_data_barrel);
//    removeNegativeBins(h_eff_data_barrel2);
//    removeNegativeBins(h_eff_data_endcap);
//    removeNegativeBins(h_eff_data_endcap2);
//    removeNegativeBins(h_eff_data_eta);

    TH1D *h_ineff_data_barrel  = ((TH1D*)(h_pT_barrel_fail [0]->Clone("h_1-PR_subtract_barrel")));
    TH1D *h_ineff_data_barrel2 = ((TH1D*)(h_pT_barrel2_fail[0]->Clone("h_1-PR_subtract2_barrel")));
    TH1D *h_ineff_data_endcap  = ((TH1D*)(h_pT_endcap_fail [0]->Clone("h_1-PR_subtract_endcap")));
    TH1D *h_ineff_data_endcap2 = ((TH1D*)(h_pT_endcap2_fail[0]->Clone("h_1-PR_subtract_endcap2")));
    TH1D *h_ineff_data_eta     = ((TH1D*)(h_eta_fail       [0]->Clone("h_1-PR_subtract_eta")));
    TH1D *h_ineff_MC_barrel    = ((TH1D*)(h_pT_barrel_fail [1]->Clone("h_1-PR_MC_barrel")));
    TH1D *h_ineff_MC_barrel2   = ((TH1D*)(h_pT_barrel2_fail[1]->Clone("h_1-PR_MC_barrel2")));
    TH1D *h_ineff_MC_endcap    = ((TH1D*)(h_pT_endcap_fail [1]->Clone("h_1-PR_MC_endcap")));
    TH1D *h_ineff_MC_endcap2   = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_1-PR_MC_endcap2")));
    TH1D *h_ineff_MC_eta       = ((TH1D*)(h_eta_fail       [1]->Clone("h_1-PR_MC_eta")));
    TH1D *h_ineff_bkgr_barrel  = ((TH1D*)(h_pT_barrel_fail [2]->Clone("h_1-PR_bkgr_barrel")));
    TH1D *h_ineff_bkgr_barrel2 = ((TH1D*)(h_pT_barrel2_fail[2]->Clone("h_1-PR_bkgr_barrel2")));
    TH1D *h_ineff_bkgr_endcap  = ((TH1D*)(h_pT_endcap_fail [2]->Clone("h_1-PR_bkgr_endcap")));
    TH1D *h_ineff_bkgr_endcap2 = ((TH1D*)(h_pT_endcap2_fail[2]->Clone("h_1-PR_bkgr_endcap2")));
    TH1D *h_ineff_bkgr_eta     = ((TH1D*)(h_eta_fail       [2]->Clone("h_1-PR_bkgr_eta")));
    TH1D *h_ineff_bkgf_barrel  = ((TH1D*)(h_pT_barrel_fail [3]->Clone("h_1-PR_bkgf_barrel")));
    TH1D *h_ineff_bkgf_barrel2 = ((TH1D*)(h_pT_barrel2_fail[3]->Clone("h_1-PR_bkgf_barrel2")));
    TH1D *h_ineff_bkgf_endcap  = ((TH1D*)(h_pT_endcap_fail [3]->Clone("h_1-PR_bkgf_endcap")));
    TH1D *h_ineff_bkgf_endcap2 = ((TH1D*)(h_pT_endcap2_fail[3]->Clone("h_1-PR_bkgf_endcap2")));
    TH1D *h_ineff_bkgf_eta     = ((TH1D*)(h_eta_fail       [3]->Clone("h_1-PR_bkgf_eta")));
    h_ineff_data_barrel ->SetDirectory(0);
    h_ineff_data_barrel2->SetDirectory(0);
    h_ineff_data_endcap ->SetDirectory(0);
    h_ineff_data_endcap2->SetDirectory(0);
    h_ineff_data_eta    ->SetDirectory(0);
    h_ineff_MC_barrel   ->SetDirectory(0);
    h_ineff_MC_barrel2  ->SetDirectory(0);
    h_ineff_MC_endcap   ->SetDirectory(0);
    h_ineff_MC_endcap2  ->SetDirectory(0);
    h_ineff_MC_eta      ->SetDirectory(0);
    h_ineff_bkgr_barrel ->SetDirectory(0);
    h_ineff_bkgr_barrel2->SetDirectory(0);
    h_ineff_bkgr_endcap ->SetDirectory(0);
    h_ineff_bkgr_endcap2->SetDirectory(0);
    h_ineff_bkgr_eta    ->SetDirectory(0);
    h_ineff_bkgf_barrel ->SetDirectory(0);
    h_ineff_bkgf_barrel2->SetDirectory(0);
    h_ineff_bkgf_endcap ->SetDirectory(0);
    h_ineff_bkgf_endcap2->SetDirectory(0);
    h_ineff_bkgf_eta    ->SetDirectory(0);
//    h_ineff_data_barrel ->Add(h_ineff_bkgf_barrel,  -1);
//    h_ineff_data_barrel2->Add(h_ineff_bkgf_barrel2, -1);
//    h_ineff_data_endcap ->Add(h_ineff_bkgf_endcap,  -1);
//    h_ineff_data_endcap2->Add(h_ineff_bkgf_endcap2, -1);
//    h_ineff_data_eta    ->Add(h_ineff_bkgf_eta,     -1);
    h_ineff_data_barrel ->Add(h_ineff_bkgr_barrel,  -1);
    h_ineff_data_barrel2->Add(h_ineff_bkgr_barrel2, -1);
    h_ineff_data_endcap ->Add(h_ineff_bkgr_endcap,  -1);
    h_ineff_data_endcap2->Add(h_ineff_bkgr_endcap2, -1);
    h_ineff_data_eta    ->Add(h_ineff_bkgr_eta,     -1);
    // FOR ERRORS
    Double_t val_control_barrel[nPtBinEndcap], val_control_barrel2[nPtBinEndcap], val_control_endcap[nPtBinEndcap], val_control_endcap2[nPtBinEndcap];
    Double_t err_control_barrel[nPtBinEndcap], err_control_barrel2[nPtBinEndcap], err_control_endcap[nPtBinEndcap], err_control_endcap2[nPtBinEndcap];
    cout << "CONTROL" << endl;
    for (Int_t i_bin=1; i_bin<=nPtBinEndcap; i_bin++)
    {
        val_control_barrel[i_bin-1] = h_ineff_data_barrel->GetBinContent(i_bin);
        err_control_barrel[i_bin-1] = h_ineff_data_barrel->GetBinError(i_bin);
        val_control_barrel2[i_bin-1] = h_ineff_data_barrel2->GetBinContent(i_bin);
        err_control_barrel2[i_bin-1] = h_ineff_data_barrel2->GetBinError(i_bin);
        val_control_endcap[i_bin-1] = h_ineff_data_endcap->GetBinContent(i_bin);
        err_control_endcap[i_bin-1] = h_ineff_data_endcap->GetBinError(i_bin);
        val_control_endcap2[i_bin-1] = h_ineff_data_endcap2->GetBinContent(i_bin);
        err_control_endcap2[i_bin-1] = h_ineff_data_endcap2->GetBinError(i_bin);

        cout << "Bin" << i_bin << ":\t" << val_control_barrel[i_bin-1] << "+-" << err_control_barrel[i_bin-1] << "(b)\t";
        cout << val_control_barrel2[i_bin-1] << "+-" << err_control_barrel2[i_bin-1] << "(b2)\t";
        cout << val_control_endcap[i_bin-1] << "+-" << err_control_endcap[i_bin-1] << "(e)\t";
        cout << val_control_endcap2[i_bin-1] << "+-" << err_control_endcap2[i_bin-1] << "(e2)\t" << endl;
    }
//    removeNegativeBins(h_ineff_data_barrel);
//    removeNegativeBins(h_ineff_data_barrel2);
//    removeNegativeBins(h_ineff_data_endcap);
//    removeNegativeBins(h_ineff_data_endcap2);
//    removeNegativeBins(h_ineff_data_eta);

    TH1D *h_eff_data_barrel_deno  = ((TH1D*)(h_eff_data_barrel ->Clone("h_PR_subtract_barrel_deno")));
    TH1D *h_eff_data_barrel2_deno = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_subtract_barrel2_deno")));
    TH1D *h_eff_data_endcap_deno  = ((TH1D*)(h_eff_data_endcap ->Clone("h_PR_subtract_endcap_deno")));
    TH1D *h_eff_data_endcap2_deno = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_endcap2_deno")));
    TH1D *h_eff_data_eta_deno = ((TH1D*)(h_eff_data_eta        ->Clone("h_PR_subtract_eta_deno")));
    h_eff_data_barrel_deno ->Add(h_ineff_data_barrel);
    h_eff_data_barrel2_deno->Add(h_ineff_data_barrel2);
    h_eff_data_endcap_deno ->Add(h_ineff_data_endcap);
    h_eff_data_endcap2_deno->Add(h_ineff_data_endcap2);
    h_eff_data_eta_deno    ->Add(h_ineff_data_eta);

    TH1D *h_eff_MC_barrel_deno  = ((TH1D*)(h_eff_MC_barrel ->Clone("h_PR_MC_barrel_deno")));
    TH1D *h_eff_MC_barrel2_deno = ((TH1D*)(h_eff_MC_barrel2->Clone("h_PR_MC_barrel2_deno")));
    TH1D *h_eff_MC_endcap_deno  = ((TH1D*)(h_eff_MC_endcap ->Clone("h_PR_MC_endcap_deno")));
    TH1D *h_eff_MC_endcap2_deno = ((TH1D*)(h_eff_MC_endcap2->Clone("h_PR_MC_endcap2_deno")));
    TH1D *h_eff_MC_eta_deno     = ((TH1D*)(h_eff_MC_eta    ->Clone("h_PR_MC_eta_deno")));
    h_eff_MC_barrel_deno ->Add(h_ineff_MC_barrel);
    h_eff_MC_barrel2_deno->Add(h_ineff_MC_barrel2);
    h_eff_MC_endcap_deno ->Add(h_ineff_MC_endcap);
    h_eff_MC_endcap2_deno->Add(h_ineff_MC_endcap2);
    h_eff_MC_eta_deno    ->Add(h_ineff_MC_eta);

    h_eff_data_barrel_deno ->SetDirectory(0);
    h_eff_data_barrel2_deno->SetDirectory(0);
    h_eff_data_endcap_deno ->SetDirectory(0);
    h_eff_data_endcap2_deno->SetDirectory(0);
    h_eff_data_eta_deno    ->SetDirectory(0);
    h_eff_MC_barrel_deno   ->SetDirectory(0);
    h_eff_MC_barrel2_deno  ->SetDirectory(0);
    h_eff_MC_endcap_deno   ->SetDirectory(0);
    h_eff_MC_endcap2_deno  ->SetDirectory(0);
    h_eff_MC_eta_deno      ->SetDirectory(0);

    h_eff_data_barrel   ->Divide(h_eff_data_barrel_deno);
    h_eff_data_barrel2  ->Divide(h_eff_data_barrel2_deno);
    h_eff_data_endcap   ->Divide(h_eff_data_endcap_deno);
    h_eff_data_endcap2  ->Divide(h_eff_data_endcap2_deno);
    h_eff_data_eta      ->Divide(h_eff_data_eta_deno);
    h_eff_MC_barrel     ->Divide(h_eff_MC_barrel_deno);
    h_eff_MC_barrel2    ->Divide(h_eff_MC_barrel2_deno);
    h_eff_MC_endcap     ->Divide(h_eff_MC_endcap_deno);
    h_eff_MC_endcap2    ->Divide(h_eff_MC_endcap2_deno);
    h_eff_MC_eta        ->Divide(h_eff_MC_eta_deno);
    h_ineff_data_barrel ->Divide(h_eff_data_barrel_deno);
    h_ineff_data_barrel2->Divide(h_eff_data_barrel2_deno);
    h_ineff_data_endcap ->Divide(h_eff_data_endcap_deno);
    h_ineff_data_endcap2->Divide(h_eff_data_endcap2_deno);
    h_ineff_data_eta    ->Divide(h_eff_data_eta_deno);
    h_ineff_MC_barrel   ->Divide(h_eff_MC_barrel_deno);
    h_ineff_MC_barrel2  ->Divide(h_eff_MC_barrel2_deno);
    h_ineff_MC_endcap   ->Divide(h_eff_MC_endcap_deno);
    h_ineff_MC_endcap2  ->Divide(h_eff_MC_endcap2_deno);
    h_ineff_MC_eta      ->Divide(h_eff_MC_eta_deno);

    // Getting the errrors
    cout << "PR ERRORS" << endl;
    for (Int_t i_bin=1; i_bin<=nPtBinEndcap; i_bin++)
    {
        Double_t err_barrel = 1;
        if (val_signal_barrel[i_bin-1] > 0 && val_control_barrel[i_bin-1] > 0) // calculate normally
        {
            err_barrel = sqrt((err_signal_barrel[i_bin-1]*err_signal_barrel[i_bin-1])/(val_signal_barrel[i_bin-1]*val_signal_barrel[i_bin-1]) +
                              (err_control_barrel[i_bin-1]*err_control_barrel[i_bin-1])/(val_control_barrel[i_bin-1]*val_control_barrel[i_bin-1]));
            err_barrel *= val_signal_barrel[i_bin-1] * val_control_barrel[i_bin-1] / ((val_signal_barrel[i_bin-1] + val_control_barrel[i_bin-1]) *
                                                                                      (val_signal_barrel[i_bin-1] + val_control_barrel[i_bin-1]));
        }
        else if (val_signal_barrel[i_bin-1] > 0) // remove control part from calculation
            err_barrel = err_signal_barrel[i_bin-1]/val_signal_barrel[i_bin-1];
        else if (val_control_barrel[i_bin-1] > 0) // remove control part from calculation
            err_barrel = err_control_barrel[i_bin-1]/val_control_barrel[i_bin-1];

        if (err_barrel > 1) err_barrel = 1;
        h_eff_data_barrel->SetBinError(i_bin, err_barrel);

        Double_t err_barrel2 = 1;
        if (val_signal_barrel2[i_bin-1] > 0 && val_control_barrel2[i_bin-1] > 0) // calculate normally
        {
            err_barrel2 = sqrt((err_signal_barrel2[i_bin-1]*err_signal_barrel2[i_bin-1])/(val_signal_barrel2[i_bin-1]*val_signal_barrel2[i_bin-1]) +
                              (err_control_barrel2[i_bin-1]*err_control_barrel2[i_bin-1])/(val_control_barrel2[i_bin-1]*val_control_barrel2[i_bin-1]));
            err_barrel2 *= val_signal_barrel2[i_bin-1] * val_control_barrel2[i_bin-1] / ((val_signal_barrel2[i_bin-1] + val_control_barrel2[i_bin-1]) *
                                                                                      (val_signal_barrel2[i_bin-1] + val_control_barrel2[i_bin-1]));
        }
        else if (val_signal_barrel2[i_bin-1] > 0) // remove control part from calculation
            err_barrel2 = err_signal_barrel2[i_bin-1]/val_signal_barrel2[i_bin-1];
        else if (val_control_barrel2[i_bin-1] > 0) // remove control part from calculation
            err_barrel2 = err_control_barrel2[i_bin-1]/val_control_barrel2[i_bin-1];

        if (err_barrel2 > 1) err_barrel2 = 1;
        h_eff_data_barrel2->SetBinError(i_bin, err_barrel2);

        Double_t err_endcap = 1;
        if (val_signal_endcap[i_bin-1] > 0 && val_control_endcap[i_bin-1] > 0) // calculate normally
        {
            err_endcap = sqrt((err_signal_endcap[i_bin-1]*err_signal_endcap[i_bin-1])/(val_signal_endcap[i_bin-1]*val_signal_endcap[i_bin-1]) +
                              (err_control_endcap[i_bin-1]*err_control_endcap[i_bin-1])/(val_control_endcap[i_bin-1]*val_control_endcap[i_bin-1]));
            err_endcap *= val_signal_endcap[i_bin-1] * val_control_endcap[i_bin-1] / ((val_signal_endcap[i_bin-1] + val_control_endcap[i_bin-1]) *
                                                                                      (val_signal_endcap[i_bin-1] + val_control_endcap[i_bin-1]));
        }
        else if (val_signal_endcap[i_bin-1] > 0) // remove control part from calculation
            err_endcap = err_signal_endcap[i_bin-1]/val_signal_endcap[i_bin-1];
        else if (val_control_endcap[i_bin-1] > 0) // remove control part from calculation
            err_endcap = err_control_endcap[i_bin-1]/val_control_endcap[i_bin-1];

        if (err_endcap > 1) err_endcap = 1;
        h_eff_data_endcap->SetBinError(i_bin, err_endcap);

        Double_t err_endcap2 = 1;
        if (val_signal_endcap2[i_bin-1] > 0 && val_control_endcap2[i_bin-1] > 0) // calculate normally
        {
            err_endcap2 = sqrt((err_signal_endcap2[i_bin-1]*err_signal_endcap2[i_bin-1])/(val_signal_endcap2[i_bin-1]*val_signal_endcap2[i_bin-1]) +
                               (err_control_endcap2[i_bin-1]*err_control_endcap2[i_bin-1])/(val_control_endcap2[i_bin-1]*val_control_endcap2[i_bin-1]));
            err_endcap2 *= val_signal_endcap2[i_bin-1] * val_control_endcap2[i_bin-1] / ((val_signal_endcap2[i_bin-1] + val_control_endcap2[i_bin-1]) *
                                                                                         (val_signal_endcap2[i_bin-1] + val_control_endcap2[i_bin-1]));
        }
        else if (val_signal_endcap2[i_bin-1] > 0) // remove control part from calculation
            err_endcap2 = err_signal_endcap2[i_bin-1]/val_signal_endcap2[i_bin-1];
        else if (val_control_endcap2[i_bin-1] > 0) // remove control part from calculation
            err_endcap2 = err_control_endcap2[i_bin-1]/val_control_endcap2[i_bin-1];

        if (err_endcap2 > 1) err_endcap2 = 1;
        h_eff_data_endcap2->SetBinError(i_bin, err_endcap2);

        cout << "Bin" << i_bin << ":\t" << err_barrel << "(b)\t";
        cout << err_barrel2 << "(b2)\t";
        cout << err_endcap << "(e)\t";
        cout << err_endcap2 << "(e2)\t" << endl;
    }

    //  For systematics
    TH1D* h_PR_subtract_barrel_plus = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_subtract_plus_barrel")));
    TH1D* h_PR_subtract_barrel2_plus = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_subtract_plus_barrel2")));
    TH1D* h_PR_subtract_endcap_plus = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_subtract_plus_endcap")));
    TH1D* h_PR_subtract_endcap2_plus = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_plus_endcap2")));
    TH1D* h_PR_subtract_barrel_minus = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_subtract_minus_barrel")));
    TH1D* h_PR_subtract_barrel2_minus = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_subtract_minus_barrel2")));
    TH1D* h_PR_subtract_endcap_minus = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_subtract_minus_endcap")));
    TH1D* h_PR_subtract_endcap2_minus = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_subtract_minus_endcap2")));
    h_PR_subtract_barrel_plus  ->SetDirectory(0);
    h_PR_subtract_barrel2_plus  ->SetDirectory(0);
    h_PR_subtract_endcap_plus  ->SetDirectory(0);
    h_PR_subtract_endcap2_plus ->SetDirectory(0);
    h_PR_subtract_barrel_minus ->SetDirectory(0);
    h_PR_subtract_barrel2_minus ->SetDirectory(0);
    h_PR_subtract_endcap_minus ->SetDirectory(0);
    h_PR_subtract_endcap2_minus->SetDirectory(0);
    for (Int_t i_bin=1; i_bin<=nPtBinEndcap; i_bin++)
    {
        if (h_PR_subtract_barrel_plus ->GetBinContent(i_bin) + h_PR_subtract_barrel_plus ->GetBinError(i_bin) < 1)
            h_PR_subtract_barrel_plus ->SetBinContent(i_bin, h_PR_subtract_barrel_plus ->GetBinContent(i_bin) + h_PR_subtract_barrel_plus ->GetBinError(i_bin));
        else h_PR_subtract_barrel_plus ->SetBinContent(i_bin, 1);
        if (h_PR_subtract_barrel2_plus ->GetBinContent(i_bin) + h_PR_subtract_barrel2_plus ->GetBinError(i_bin) < 1)
            h_PR_subtract_barrel2_plus ->SetBinContent(i_bin, h_PR_subtract_barrel2_plus ->GetBinContent(i_bin) + h_PR_subtract_barrel2_plus ->GetBinError(i_bin));
        else h_PR_subtract_barrel2_plus ->SetBinContent(i_bin, 1);
        if (h_PR_subtract_endcap_plus ->GetBinContent(i_bin) + h_PR_subtract_endcap_plus ->GetBinError(i_bin) < 1)
            h_PR_subtract_endcap_plus ->SetBinContent(i_bin, h_PR_subtract_endcap_plus ->GetBinContent(i_bin) + h_PR_subtract_endcap_plus ->GetBinError(i_bin));
        else h_PR_subtract_endcap_plus ->SetBinContent(i_bin, 1);
        if (h_PR_subtract_endcap2_plus->GetBinContent(i_bin) + h_PR_subtract_endcap2_plus->GetBinError(i_bin) < 1)
            h_PR_subtract_endcap2_plus->SetBinContent(i_bin, h_PR_subtract_endcap2_plus->GetBinContent(i_bin) + h_PR_subtract_endcap2_plus->GetBinError(i_bin));
        else h_PR_subtract_endcap2_plus->SetBinContent(i_bin, 1);

        if (h_PR_subtract_barrel_minus ->GetBinContent(i_bin) - h_PR_subtract_barrel_minus ->GetBinError(i_bin) > 0)
            h_PR_subtract_barrel_minus ->SetBinContent(i_bin, h_PR_subtract_barrel_minus ->GetBinContent(i_bin) - h_PR_subtract_barrel_minus ->GetBinError(i_bin));
        else h_PR_subtract_barrel_minus ->SetBinContent(i_bin, 0);
        if (h_PR_subtract_barrel2_minus ->GetBinContent(i_bin) - h_PR_subtract_barrel2_minus ->GetBinError(i_bin) > 0)
            h_PR_subtract_barrel2_minus ->SetBinContent(i_bin, h_PR_subtract_barrel2_minus ->GetBinContent(i_bin) - h_PR_subtract_barrel2_minus ->GetBinError(i_bin));
        else h_PR_subtract_barrel2_minus ->SetBinContent(i_bin, 0);
        if (h_PR_subtract_endcap_minus ->GetBinContent(i_bin) - h_PR_subtract_endcap_minus ->GetBinError(i_bin) > 0)
            h_PR_subtract_endcap_minus ->SetBinContent(i_bin, h_PR_subtract_endcap_minus ->GetBinContent(i_bin) - h_PR_subtract_endcap_minus ->GetBinError(i_bin));
        else h_PR_subtract_endcap_minus ->SetBinContent(i_bin, 0);
        if (h_PR_subtract_endcap2_minus->GetBinContent(i_bin) - h_PR_subtract_endcap2_minus->GetBinError(i_bin) > 0)
            h_PR_subtract_endcap2_minus->SetBinContent(i_bin, h_PR_subtract_endcap2_minus->GetBinContent(i_bin) - h_PR_subtract_endcap2_minus->GetBinError(i_bin));
        else h_PR_subtract_endcap2_minus->SetBinContent(i_bin,0);
    }

    // Ratio method
    TH1D *h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [1]->Clone("h_PR_ratio_barrel_deno")));
    TH1D *h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_ratio_barrel2_deno")));
    TH1D *h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [1]->Clone("h_PR_ratio_endcap_deno")));
    TH1D *h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_ratio_endcap2_deno")));
    TH1D *h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [1]->Clone("h_PR_ratio_eta_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_pass [2]);
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_pass [3]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_pass[2]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_pass[3]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_pass [2]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_pass [3]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_pass[2]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_pass[3]);
    h_PR_ratio_eta_deno    ->Add(h_eta_pass       [2]);
    h_PR_ratio_eta_deno    ->Add(h_eta_pass       [3]);

    TH1D *h_PR_ratio_barrel  = ((TH1D*)(h_pT_barrel_fail [1]->Clone("h_PR_ratio_barrel")));
    TH1D *h_PR_ratio_barrel2 = ((TH1D*)(h_pT_barrel2_fail[1]->Clone("h_PR_ratio_barrel2")));
    TH1D *h_PR_ratio_endcap  = ((TH1D*)(h_pT_endcap_fail [1]->Clone("h_PR_ratio_endcap")));
    TH1D *h_PR_ratio_endcap2 = ((TH1D*)(h_pT_endcap2_fail[1]->Clone("h_PR_ratio_endcap2")));
    TH1D *h_PR_ratio_eta     = ((TH1D*)(h_eta_fail       [1]->Clone("h_PR_ratio_eta")));
    h_PR_ratio_barrel ->Add(h_pT_barrel_fail[2]);
    h_PR_ratio_barrel ->Add(h_pT_barrel_fail[3]);
    h_PR_ratio_barrel ->Add(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Add(h_pT_barrel2_fail[2]);
    h_PR_ratio_barrel2->Add(h_pT_barrel2_fail[3]);
    h_PR_ratio_barrel2->Add(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Add(h_pT_endcap_fail[2]);
    h_PR_ratio_endcap ->Add(h_pT_endcap_fail[3]);
    h_PR_ratio_endcap ->Add(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Add(h_pT_endcap2_fail[2]);
    h_PR_ratio_endcap2->Add(h_pT_endcap2_fail[3]);
    h_PR_ratio_endcap2->Add(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Add(h_eta_fail[2]);
    h_PR_ratio_eta    ->Add(h_eta_fail[3]);
    h_PR_ratio_eta    ->Add(h_PR_ratio_eta_deno);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);

    h_PR_ratio_barrel ->Multiply(h_pT_barrel_pass [0]);
    h_PR_ratio_barrel2->Multiply(h_pT_barrel2_pass[0]);
    h_PR_ratio_endcap ->Multiply(h_pT_endcap_pass [0]);
    h_PR_ratio_endcap2->Multiply(h_pT_endcap2_pass[0]);
    h_PR_ratio_eta    ->Multiply(h_eta_pass       [0]);

    h_PR_ratio_barrel ->Multiply(h_pT_barrel_pass [1]);
    h_PR_ratio_barrel2->Multiply(h_pT_barrel2_pass[1]);
    h_PR_ratio_endcap ->Multiply(h_pT_endcap_pass [1]);
    h_PR_ratio_endcap2->Multiply(h_pT_endcap2_pass[1]);
    h_PR_ratio_eta    ->Multiply(h_eta_pass       [1]);

    h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [0]->Clone("h_PR_ratio_barrel_deno")));
    h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[0]->Clone("h_PR_ratio_barrel2_deno")));
    h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [0]->Clone("h_PR_ratio_endcap_deno")));
    h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[0]->Clone("h_PR_ratio_endcap2_deno")));
    h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [0]->Clone("h_PR_ratio_eta_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_fail [0]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_fail[0]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_fail [0]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_fail[0]);
    h_PR_ratio_eta_deno    ->Add(h_eta_fail       [0]);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);

    h_PR_ratio_barrel_deno  = ((TH1D*)(h_pT_barrel_pass [1]->Clone("h_PR_ratio_barrel_deno")));
    h_PR_ratio_barrel2_deno = ((TH1D*)(h_pT_barrel2_pass[1]->Clone("h_PR_ratio_barrel2_deno")));
    h_PR_ratio_endcap_deno  = ((TH1D*)(h_pT_endcap_pass [1]->Clone("h_PR_ratio_endcap_deno")));
    h_PR_ratio_endcap2_deno = ((TH1D*)(h_pT_endcap2_pass[1]->Clone("h_PR_ratio_endcap2_deno")));
    h_PR_ratio_eta_deno     = ((TH1D*)(h_eta_pass       [1]->Clone("h_PR_ratio_eta_deno")));
    h_PR_ratio_barrel_deno ->Add(h_pT_barrel_fail [1]);
    h_PR_ratio_barrel2_deno->Add(h_pT_barrel2_fail[1]);
    h_PR_ratio_endcap_deno ->Add(h_pT_endcap_fail [1]);
    h_PR_ratio_endcap2_deno->Add(h_pT_endcap2_fail[1]);
    h_PR_ratio_eta_deno    ->Add(h_eta_fail       [1]);

    h_PR_ratio_barrel ->Divide(h_PR_ratio_barrel_deno);
    h_PR_ratio_barrel2->Divide(h_PR_ratio_barrel2_deno);
    h_PR_ratio_endcap ->Divide(h_PR_ratio_endcap_deno);
    h_PR_ratio_endcap2->Divide(h_PR_ratio_endcap2_deno);
    h_PR_ratio_eta    ->Divide(h_PR_ratio_eta_deno);

    // Drawing
    h_eff_data_barrel->SetDirectory(0);
    h_eff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel->SetMarkerColor(kBlack);
    h_eff_data_barrel->SetLineColor(kBlack);
    h_eff_data_barrel2->SetDirectory(0);
    h_eff_data_barrel2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel2->SetMarkerColor(kBlack);
    h_eff_data_barrel2->SetLineColor(kBlack);
    h_eff_data_endcap->SetDirectory(0);
    h_eff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap->SetMarkerColor(kBlack);
    h_eff_data_endcap->SetLineColor(kBlack);
    h_eff_data_endcap2->SetDirectory(0);
    h_eff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2->SetMarkerColor(kBlack);
    h_eff_data_endcap2->SetLineColor(kBlack);
    h_eff_data_eta->SetDirectory(0);
    h_eff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_eff_data_eta->SetMarkerColor(kBlack);
    h_eff_data_eta->SetLineColor(kBlack);
    h_eff_MC_barrel->SetDirectory(0);
    h_eff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel->SetMarkerColor(kRed);
    h_eff_MC_barrel->SetLineColor(kRed);
    h_eff_MC_barrel2->SetDirectory(0);
    h_eff_MC_barrel2->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel2->SetMarkerColor(kRed);
    h_eff_MC_barrel2->SetLineColor(kRed);
    h_eff_MC_endcap->SetDirectory(0);
    h_eff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap->SetMarkerColor(kRed);
    h_eff_MC_endcap->SetLineColor(kRed);
    h_eff_MC_endcap2->SetDirectory(0);
    h_eff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap2->SetMarkerColor(kRed);
    h_eff_MC_endcap2->SetLineColor(kRed);
    h_eff_MC_eta->SetDirectory(0);
    h_eff_MC_eta->SetMarkerStyle(kFullSquare);
    h_eff_MC_eta->SetMarkerColor(kRed);
    h_eff_MC_eta->SetLineColor(kRed);
    h_PR_ratio_barrel->SetDirectory(0);
    h_PR_ratio_barrel->SetMarkerStyle(33);
    h_PR_ratio_barrel->SetMarkerSize(1.5);
    h_PR_ratio_barrel->SetMarkerColor(kBlue);
    h_PR_ratio_barrel->SetLineColor(kBlue);
    h_PR_ratio_barrel2->SetDirectory(0);
    h_PR_ratio_barrel2->SetMarkerStyle(33);
    h_PR_ratio_barrel2->SetMarkerSize(1.5);
    h_PR_ratio_barrel2->SetMarkerColor(kBlue);
    h_PR_ratio_barrel2->SetLineColor(kBlue);
    h_PR_ratio_endcap->SetDirectory(0);
    h_PR_ratio_endcap->SetMarkerStyle(33);
    h_PR_ratio_endcap->SetMarkerSize(1.5);
    h_PR_ratio_endcap->SetMarkerColor(kBlue);
    h_PR_ratio_endcap->SetLineColor(kBlue);
    h_PR_ratio_endcap2->SetDirectory(0);
    h_PR_ratio_endcap2->SetMarkerStyle(33);
    h_PR_ratio_endcap2->SetMarkerSize(1.5);
    h_PR_ratio_endcap2->SetMarkerColor(kBlue);
    h_PR_ratio_endcap2->SetLineColor(kBlue);
    h_PR_ratio_eta->SetDirectory(0);
    h_PR_ratio_eta->SetMarkerStyle(33);
    h_PR_ratio_eta->SetMarkerSize(1.5);
    h_PR_ratio_eta->SetMarkerColor(kBlue);
    h_PR_ratio_eta->SetLineColor(kBlue);

    h_ineff_data_barrel->SetDirectory(0);
    h_ineff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel->SetMarkerColor(kBlack);
    h_ineff_data_barrel->SetLineColor(kBlack);
    h_ineff_data_barrel2->SetDirectory(0);
    h_ineff_data_barrel2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel2->SetMarkerColor(kBlack);
    h_ineff_data_barrel2->SetLineColor(kBlack);
    h_ineff_data_endcap->SetDirectory(0);
    h_ineff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap->SetMarkerColor(kBlack);
    h_ineff_data_endcap->SetLineColor(kBlack);
    h_ineff_data_endcap2->SetDirectory(0);
    h_ineff_data_endcap2->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap2->SetMarkerColor(kBlack);
    h_ineff_data_endcap2->SetLineColor(kBlack);
    h_ineff_data_eta->SetDirectory(0);
    h_ineff_data_eta->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_eta->SetMarkerColor(kBlack);
    h_ineff_data_eta->SetLineColor(kBlack);
    h_ineff_MC_barrel->SetDirectory(0);
    h_ineff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel->SetMarkerColor(kRed);
    h_ineff_MC_barrel->SetLineColor(kRed);
    h_ineff_MC_barrel2->SetDirectory(0);
    h_ineff_MC_barrel2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel2->SetMarkerColor(kRed);
    h_ineff_MC_barrel2->SetLineColor(kRed);
    h_ineff_MC_endcap->SetDirectory(0);
    h_ineff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap->SetMarkerColor(kRed);
    h_ineff_MC_endcap->SetLineColor(kRed);
    h_ineff_MC_endcap2->SetDirectory(0);
    h_ineff_MC_endcap2->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap2->SetMarkerColor(kRed);
    h_ineff_MC_endcap2->SetLineColor(kRed);
    h_ineff_MC_eta->SetDirectory(0);
    h_ineff_MC_eta->SetMarkerStyle(kFullSquare);
    h_ineff_MC_eta->SetMarkerColor(kRed);
    h_ineff_MC_eta->SetLineColor(kRed);

    TH1D *h_eff_ratio_barrel    = ((TH1D*)(h_eff_data_barrel   ->Clone("h_eff_ratio_barrel")));
    TH1D *h_eff_ratio_barrel2   = ((TH1D*)(h_eff_data_barrel2  ->Clone("h_eff_ratio_barrel2")));
    TH1D *h_eff_ratio_endcap    = ((TH1D*)(h_eff_data_endcap   ->Clone("h_eff_ratio_endcap")));
    TH1D *h_eff_ratio_endcap2   = ((TH1D*)(h_eff_data_endcap2  ->Clone("h_eff_ratio_endcap2")));
    TH1D *h_eff_ratio_eta       = ((TH1D*)(h_eff_data_eta      ->Clone("h_eff_ratio_eta")));
    TH1D *h_ineff_ratio_barrel  = ((TH1D*)(h_ineff_data_barrel ->Clone("h_1-eff_ratio_barrel")));
    TH1D *h_ineff_ratio_barrel2 = ((TH1D*)(h_ineff_data_barrel2->Clone("h_1-eff_ratio_barrel2")));
    TH1D *h_ineff_ratio_endcap  = ((TH1D*)(h_ineff_data_endcap ->Clone("h_1-eff_ratio_endcap")));
    TH1D *h_ineff_ratio_endcap2 = ((TH1D*)(h_ineff_data_endcap2->Clone("h_1-eff_ratio_endcap2")));
    TH1D *h_ineff_ratio_eta     = ((TH1D*)(h_ineff_data_eta    ->Clone("h_1-eff_ratio_eta")));
    h_eff_ratio_barrel ->Divide(h_eff_MC_barrel);
    h_eff_ratio_barrel2->Divide(h_eff_MC_barrel);
    h_eff_ratio_endcap ->Divide(h_eff_MC_endcap);
    h_eff_ratio_endcap2->Divide(h_eff_MC_endcap2);
    h_eff_ratio_eta    ->Divide(h_eff_MC_eta);
    h_eff_ratio_barrel ->SetDirectory(0);
    h_eff_ratio_barrel2->SetDirectory(0);
    h_eff_ratio_barrel2->SetMarkerColor(kBlue);
    h_eff_ratio_barrel2->SetLineColor(kBlue);
    h_eff_ratio_endcap ->SetDirectory(0);
    h_eff_ratio_endcap ->SetMarkerColor(kBlue);
    h_eff_ratio_endcap ->SetLineColor(kBlue);
    h_eff_ratio_endcap2->SetDirectory(0);
    h_eff_ratio_endcap2->SetMarkerColor(kBlue);
    h_eff_ratio_endcap2->SetLineColor(kBlue);
    h_eff_ratio_eta    ->SetDirectory(0);
    h_eff_ratio_eta    ->SetMarkerColor(kBlue);
    h_eff_ratio_eta    ->SetLineColor(kBlue);
    h_ineff_ratio_barrel ->Divide(h_ineff_MC_barrel);
    h_ineff_ratio_barrel2->Divide(h_ineff_MC_barrel2);
    h_ineff_ratio_endcap ->Divide(h_ineff_MC_endcap);
    h_ineff_ratio_endcap2->Divide(h_ineff_MC_endcap2);
    h_ineff_ratio_eta    ->Divide(h_ineff_MC_eta);
    h_ineff_ratio_barrel ->SetDirectory(0);
    h_ineff_ratio_barrel2->SetDirectory(0);
    h_ineff_ratio_barrel2->SetMarkerColor(kBlue);
    h_ineff_ratio_barrel2->SetLineColor(kBlue);
    h_ineff_ratio_endcap ->SetDirectory(0);
    h_ineff_ratio_endcap ->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap ->SetLineColor(kBlue);
    h_ineff_ratio_endcap2->SetDirectory(0);
    h_ineff_ratio_endcap2->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap2->SetLineColor(kBlue);
    h_ineff_ratio_eta    ->SetDirectory(0);
    h_ineff_ratio_eta    ->SetMarkerColor(kBlue);
    h_ineff_ratio_eta    ->SetLineColor(kBlue);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_eff_barrel    = new myRatioPlot_t("RP_PR_barrel",    h_eff_MC_barrel,    h_eff_data_barrel);
    myRatioPlot_t *RP_eff_barrel2   = new myRatioPlot_t("RP_PR_barrel2",   h_eff_MC_barrel2,   h_eff_data_barrel2);
    myRatioPlot_t *RP_eff_endcap    = new myRatioPlot_t("RP_PR_endcap",    h_eff_MC_endcap,    h_eff_data_endcap);
    myRatioPlot_t *RP_eff_endcap2   = new myRatioPlot_t("RP_PR_endcap2",   h_eff_MC_endcap2,   h_eff_data_endcap2);
    myRatioPlot_t *RP_eff_eta       = new myRatioPlot_t("RP_PR_eta",       h_eff_MC_eta,       h_eff_data_eta);
    myRatioPlot_t *RP_ineff_barrel  = new myRatioPlot_t("RP_1-PR_barrel",  h_ineff_MC_barrel,  h_ineff_data_barrel);
    myRatioPlot_t *RP_ineff_barrel2 = new myRatioPlot_t("RP_1-PR_barrel2", h_ineff_MC_barrel2, h_ineff_data_barrel2);
    myRatioPlot_t *RP_ineff_endcap  = new myRatioPlot_t("RP_1-PR_endcap",  h_ineff_MC_endcap,  h_ineff_data_endcap);
    myRatioPlot_t *RP_ineff_endcap2 = new myRatioPlot_t("RP_1-PR_endcap2", h_ineff_MC_endcap2, h_ineff_data_endcap2);
    myRatioPlot_t *RP_ineff_eta     = new myRatioPlot_t("RP_1-PR_eta",     h_ineff_MC_eta,     h_ineff_data_eta);

    RP_eff_barrel   ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_barrel2  ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap   ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap2  ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_eta      ->SetPlots("#eta", -2.4, 2.4);
    RP_ineff_barrel ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_barrel2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap ->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap2->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_eta    ->SetPlots("#eta", -2.4, 2.4);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_eff_data_barrel, "Data", "pl");
    legend->AddEntry(h_eff_MC_barrel, "W+Jets MC", "pl");

    RP_eff_barrel   ->ImportLegend(legend);
    RP_eff_barrel2  ->ImportLegend(legend);
    RP_eff_endcap   ->ImportLegend(legend);
    RP_eff_endcap2  ->ImportLegend(legend);
    RP_eff_eta      ->ImportLegend(legend);
    RP_ineff_barrel ->ImportLegend(legend);
    RP_ineff_barrel2->ImportLegend(legend);
    RP_ineff_endcap ->ImportLegend(legend);
    RP_ineff_endcap2->ImportLegend(legend);
    RP_ineff_eta    ->ImportLegend(legend);

//    RP_eff_barrel   ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_barrel2  ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap   ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_endcap2  ->Draw(0.01, 1, 1, "", "Prompt rate");
//    RP_eff_eta      ->Draw(0.01, 1, 0, "", "Prompt rate");
//    RP_ineff_barrel ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_barrel2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap ->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_endcap2->Draw(0.01, 1, 1, "", "1-PR");
//    RP_ineff_eta    ->Draw(0.01, 1, 1, "", "1-PR");

//    RP_eff_barrel   ->pad1->SetLogy(0);
//    RP_eff_barrel2  ->pad1->SetLogy(0);
//    RP_eff_endcap   ->pad1->SetLogy(0);
//    RP_eff_endcap2  ->pad1->SetLogy(0);
//    RP_eff_eta      ->pad1->SetLogy(0);
//    RP_ineff_barrel ->pad1->SetLogy(0);
//    RP_ineff_barrel2->pad1->SetLogy(0);
//    RP_ineff_endcap ->pad1->SetLogy(0);
//    RP_ineff_endcap2->pad1->SetLogy(0);
//    RP_ineff_eta    ->pad1->SetLogy(0);

    TCanvas *c_PR_barrel = new TCanvas("c_PR_barrel", "c_PR_barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
    TH1D *h_eff_data_barrel_draw = ((TH1D*)(h_eff_data_barrel->Clone("h_PR_data_barrel_draw")));
    TH1D *h_eff_MC_barrel_draw   = ((TH1D*)(h_eff_MC_barrel  ->Clone("h_PR_MC_barrel_draw")));
    TH1D *h_PR_ratio_barrel_draw = ((TH1D*)(h_PR_ratio_barrel->Clone("h_PR_ratio_barrel_draw")));
    h_eff_data_barrel_draw->SetDirectory(0);
    h_eff_MC_barrel_draw->SetDirectory(0);
    h_PR_ratio_barrel_draw->SetDirectory(0);
    h_eff_data_barrel_draw->SetStats(kFALSE);
    h_eff_data_barrel_draw->SetTitle("");
    h_eff_data_barrel_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_barrel_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_barrel_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_barrel_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_barrel_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_barrel_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_barrel_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_barrel_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_barrel_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_barrel_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_barrel_draw->Draw();
    h_eff_MC_barrel_draw->Draw("same");
    h_PR_ratio_barrel_draw->Draw("same");
    TLegend *legend1 = new TLegend(0.6, 0.75, 0.95, 0.95);
    legend1->AddEntry(h_eff_data_barrel_draw, "Subtraction", "lp");
    legend1->AddEntry(h_PR_ratio_barrel_draw, "Ratio", "lp");
    legend1->AddEntry(h_eff_MC_barrel_draw, "MC", "lp");
    legend1->Draw();
    // Fit
//    TF1 *f_barrel_17to43 = new TF1("f_barrel_17to43","[0]+[1]*x+[2]*sqrt(x)",17,45);
//    f_barrel_17to43->SetLineColor(kMagenta-5);
//    TF1 *f_barrel_43to56 = new TF1("f_barrel_43to56","[0]+exp([1]*([2]-x))",40,70);
//    f_barrel_43to56->SetParameter(0, 0.55);
//    f_barrel_43to56->SetParameter(1, 0.55);
//    f_barrel_43to56->SetParameter(2, 50);
//    f_barrel_43to56->SetLineColor(kMagenta-5);
//    TF1 *f_barrel_56to95 = new TF1("f_barrel_56to95","[0]+[1]*x+[2]*x^2",50,100);
//    f_barrel_56to95->SetLineColor(kMagenta-5);
//    TF1 *f_barrel_95to5000 = new TF1("f_barrel_95to5000","[0]+gaus(1)",90,5000);
//    f_barrel_95to5000->SetParameter(0, 0.4);
//    f_barrel_95to5000->SetParameter(1, 0.3);
//    f_barrel_95to5000->SetParameter(2, 90);
//    f_barrel_95to5000->SetParameter(3, 150);
//    f_barrel_95to5000->SetLineColor(kMagenta-5);
//    h_eff_data_barrel_draw->Fit("f_barrel_17to43", "R");
//    h_eff_data_barrel_draw->Fit("f_barrel_43to56", "R");
//    h_eff_data_barrel_draw->Fit("f_barrel_56to95", "R");
//    h_eff_data_barrel_draw->Fit("f_barrel_95to5000", "R");
//    f_barrel_17to43->Draw("same");
//    f_barrel_43to56->Draw("same");
//    f_barrel_56to95->Draw("same");
//    f_barrel_95to5000->Draw("same");
//    Double_t chi2_barrel_40 = f_barrel_17to43->GetChisquare();
//    Double_t chi2_barrel_70 = f_barrel_43to56->GetChisquare();
//    Double_t chi2_barrel_100 = f_barrel_56to95->GetChisquare();
//    Double_t chi2_barrel_5000 = f_barrel_95to5000->GetChisquare();
//    cout << "Barrel fit chi squares: " << chi2_barrel_40 << "  " << chi2_barrel_70<< "  " << chi2_barrel_100  << "  " << chi2_barrel_5000 << endl;
//    cout << "\n\n";
    c_PR_barrel->Update();


    TCanvas *c_PR_barrel2 = new TCanvas("c_PR_barrel2", "c_PR_barrel2", 800, 800);
    c_PR_barrel2->cd();
    c_PR_barrel2->SetGrid(1);
    c_PR_barrel2->SetLogx(1);
    c_PR_barrel2->SetRightMargin(0.05);
    c_PR_barrel2->SetTopMargin(0.05);
    c_PR_barrel2->SetBottomMargin(0.12);
    c_PR_barrel2->SetLeftMargin(0.13);
    TH1D *h_eff_data_barrel2_draw = ((TH1D*)(h_eff_data_barrel2->Clone("h_PR_data_barrel2_draw")));
    TH1D *h_eff_MC_barrel2_draw   = ((TH1D*)(h_eff_MC_barrel2  ->Clone("h_PR_MC_barrel2_draw")));
    TH1D *h_PR_ratio_barrel2_draw = ((TH1D*)(h_PR_ratio_barrel2->Clone("h_PR_ratio_barrel2_draw")));
    h_eff_data_barrel2_draw->SetDirectory(0);
    h_eff_MC_barrel2_draw->SetDirectory(0);
    h_PR_ratio_barrel2_draw->SetDirectory(0);
    h_eff_data_barrel2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel2_draw->SetMarkerColor(kBlack);
    h_eff_data_barrel2_draw->SetLineColor(kBlack);
    h_eff_data_barrel2_draw->SetStats(kFALSE);
    h_eff_data_barrel2_draw->SetTitle("");
    h_eff_data_barrel2_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_barrel2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_barrel2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_barrel2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_barrel2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_barrel2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_barrel2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_barrel2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_barrel2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_barrel2_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_barrel2_draw->Draw();
    h_eff_MC_barrel2_draw->Draw("same");
    h_PR_ratio_barrel2_draw->Draw("same");
    legend1->Draw();
    // Fit
//    TF1 *f_barrel2_17to43 = new TF1("f_barrel2_17to43","[0]+[1]*x+[2]*sqrt(x)",17,45);
//    f_barrel2_17to43->SetLineColor(kMagenta-5);
//    TF1 *f_barrel2_43to56 = new TF1("f_barrel2_43to56","[0]+exp([1]*([2]-x))",40,70);
//    f_barrel2_43to56->SetParameter(0, 0.55);
//    f_barrel2_43to56->SetParameter(1, 0.55);
//    f_barrel2_43to56->SetParameter(2, 50);
//    f_barrel2_43to56->SetLineColor(kMagenta-5);
//    TF1 *f_barrel2_56to95 = new TF1("f_barrel2_56to95","[0]+[1]*x+[2]*x^2",50,100);
//    f_barrel2_56to95->SetLineColor(kMagenta-5);
//    TF1 *f_barrel2_95to5000 = new TF1("f_barrel2_95to5000","[0]+gaus(1)",90,5000);
//    f_barrel2_95to5000->SetParameter(0, 0.4);
//    f_barrel2_95to5000->SetParameter(1, 0.3);
//    f_barrel2_95to5000->SetParameter(2, 90);
//    f_barrel2_95to5000->SetParameter(3, 150);
//    f_barrel2_95to5000->SetLineColor(kMagenta-5);
//    h_eff_data_barrel2_draw->Fit("f_barrel2_17to43", "R");
//    h_eff_data_barrel2_draw->Fit("f_barrel2_43to56", "R");
//    h_eff_data_barrel2_draw->Fit("f_barrel2_56to95", "R");
//    h_eff_data_barrel2_draw->Fit("f_barrel2_95to5000", "R");
//    f_barrel2_17to43->Draw("same");
//    f_barrel2_43to56->Draw("same");
//    f_barrel2_56to95->Draw("same");
//    f_barrel2_95to5000->Draw("same");
//    Double_t chi2_barrel2_40 = f_barrel2_17to43->GetChisquare();
//    Double_t chi2_barrel2_70 = f_barrel2_43to56->GetChisquare();
//    Double_t chi2_barrel2_100 = f_barrel2_56to95->GetChisquare();
//    Double_t chi2_barrel2_5000 = f_barrel2_95to5000->GetChisquare();
//    cout << "Far barrel fit chi squares: " << chi2_barrel2_40 << "  " << chi2_barrel2_70<< "  " << chi2_barrel2_100  << "  " << chi2_barrel2_5000 << endl;
//    cout << "\n\n";
    c_PR_barrel2->Update();


    TCanvas *c_PR_endcap = new TCanvas("c_PR_endcap", "c_PR_endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap_draw = ((TH1D*)(h_eff_data_endcap->Clone("h_PR_data_endcap_draw")));
    TH1D *h_eff_MC_endcap_draw   = ((TH1D*)(h_eff_MC_endcap  ->Clone("h_PR_MC_endcap_draw")));
    TH1D *h_PR_ratio_endcap_draw = ((TH1D*)(h_PR_ratio_endcap->Clone("h_PR_ratio_endcap_draw")));
    h_eff_data_endcap_draw->SetDirectory(0);
    h_eff_MC_endcap_draw->SetDirectory(0);
    h_PR_ratio_endcap_draw->SetDirectory(0);
    h_eff_data_endcap_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap_draw->SetLineColor(kBlack);
    h_eff_data_endcap_draw->SetStats(kFALSE);
    h_eff_data_endcap_draw->SetTitle("");
    h_eff_data_endcap_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_endcap_draw->Draw();
    h_eff_MC_endcap_draw->Draw("same");
    h_PR_ratio_endcap_draw->Draw("same");
    legend1->Draw();
    // Fit
//    TF1 *f_endcap_17to42p5 = new TF1("f_endcap_17to42p5","[0]+[1]*x+[2]*sqrt(x)",17,43);
//    f_endcap_17to42p5->SetLineColor(kAzure+1);
//    TF1 *f_endcap_42p5to55 = new TF1("f_endcap_42p5to55","[0]+[1]/(x-[2])+[3]*x",40,60);
//    f_endcap_42p5to55->SetLineColor(kAzure+1);
//    TF1 *f_endcap_55to5000 = new TF1("f_endcap_55to5000","[0]+[1]/(1+exp([2]*(x-[3])))",50,5000);
//    f_endcap_55to5000->SetParameter(0, h_eff_data_endcap_draw->GetBinContent(16));
//    f_endcap_55to5000->SetParameter(1, h_eff_data_endcap_draw->GetBinContent(11)-h_eff_data_endcap_draw->GetBinContent(16));
//    f_endcap_55to5000->SetParameter(2, 0.05);
//    f_endcap_55to5000->SetParLimits(2, 0.04, 0.15);
//    f_endcap_55to5000->SetParameter(3, 200);
//    f_endcap_55to5000->SetParLimits(3, 200, 400);
//    f_endcap_55to5000->SetLineColor(kAzure+1);
//    h_eff_data_endcap_draw->Fit("f_endcap_17to42p5", "R");
//    h_eff_data_endcap_draw->Fit("f_endcap_42p5to55", "R");
//    h_eff_data_endcap_draw->Fit("f_endcap_55to5000", "R");
//    f_endcap_17to42p5->Draw("same");
//    f_endcap_42p5to55->Draw("same");
//    f_endcap_55to5000->Draw("same");
//    Double_t chi2_endcap_60 = f_endcap_17to42p5->GetChisquare();
//    Double_t chi2_endcap_90 = f_endcap_42p5to55->GetChisquare();
//    Double_t chi2_endcap_5000 = f_endcap_55to5000->GetChisquare();
//    cout << "Endcap fit chi squares: " << chi2_endcap_60 << "  " << chi2_endcap_90 << "  " << chi2_endcap_5000 << endl;
//    cout << "\n\n";
    c_PR_endcap->Update();


    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "c_PR_endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
    TH1D *h_eff_data_endcap2_draw = ((TH1D*)(h_eff_data_endcap2->Clone("h_PR_data_endcap2_draw")));
    TH1D *h_eff_MC_endcap2_draw   = ((TH1D*)(h_eff_MC_endcap2  ->Clone("h_PR_MC_endcap2_draw")));
    TH1D *h_PR_ratio_endcap2_draw = ((TH1D*)(h_PR_ratio_endcap2->Clone("h_PR_ratio_endcap2_draw")));
    h_eff_data_endcap2_draw->SetDirectory(0);
    h_eff_MC_endcap2_draw->SetDirectory(0);
    h_PR_ratio_endcap2_draw->SetDirectory(0);
    h_eff_data_endcap2_draw->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap2_draw->SetMarkerColor(kBlack);
    h_eff_data_endcap2_draw->SetLineColor(kBlack);
    h_eff_data_endcap2_draw->SetStats(kFALSE);
    h_eff_data_endcap2_draw->SetTitle("");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleOffset(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_endcap2_draw->GetYaxis()->SetTitleOffset(1.25);
    h_eff_data_endcap2_draw->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_endcap2_draw->GetXaxis()->SetNoExponent(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetMoreLogLabels(1);
    h_eff_data_endcap2_draw->GetXaxis()->SetRangeUser(17, 3000);
    h_eff_data_endcap2_draw->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_endcap2_draw->Draw();
    h_eff_MC_endcap2_draw->Draw("same");
    h_PR_ratio_endcap2_draw->Draw("same");
    legend1->Draw();
    // Fit
//    TF1 *f_endcap2_17to42p5 = new TF1("f_endcap2_17to42p5","[0]+[1]*sqrt(x)+[2]*x",17,47);
//    f_endcap2_17to42p5->SetLineColor(kOrange+5);
//    TF1 *f_endcap2_42p5to56 = new TF1("f_endcap2_42p5to56","[0]+[1]*x+[2]/(x-[3])",42,67);
//    f_endcap2_42p5to56->SetParLimits(1, 0.001, 0.05);
//    f_endcap2_42p5to56->SetParLimits(2, 50, 70);
//    f_endcap2_42p5to56->SetParLimits(3, 10, 50);
//    f_endcap2_42p5to56->SetLineColor(kOrange+5);
//    TF1 *f_endcap2_56to76 = new TF1("f_endcap2_56to76","[0]+gaus(1)",56,110);
//    f_endcap2_56to76->FixParameter(0, h_eff_data_endcap2_draw->GetBinContent(10));
//    f_endcap2_56to76->SetParLimits(1, 0.03, 0.07);
//    f_endcap2_56to76->SetParLimits(2, 74, 76);
//    f_endcap2_56to76->SetParLimits(3, 7, 15);
//    f_endcap2_56to76->SetLineColor(kOrange+5);
//    TF1 *f_endcap2_76to5000 = new TF1("f_endcap2_76to5000","[0]+exp([1]*(x-[2]))",70,5000);
//    f_endcap2_76to5000->FixParameter(0, h_eff_data_endcap2_draw->GetBinContent(16));
//    f_endcap2_76to5000->SetParameter(1, -0.3);
//    f_endcap2_76to5000->SetParameter(2, 30);
//    f_endcap2_76to5000->SetLineColor(kOrange+5);
//    h_eff_data_endcap2_draw->Fit("f_endcap2_17to42p5", "R");
//    h_eff_data_endcap2_draw->Fit("f_endcap2_42p5to56", "R");
//    h_eff_data_endcap2_draw->Fit("f_endcap2_56to76", "R");
//    h_eff_data_endcap2_draw->Fit("f_endcap2_76to5000", "R");
//    f_endcap2_17to42p5->Draw("same");
//    f_endcap2_42p5to56->Draw("same");
//    f_endcap2_56to76->Draw("same");
//    f_endcap2_76to5000->Draw("same");
//    Double_t chi2_endcap2_45 = f_endcap2_17to42p5->GetChisquare();
//    Double_t chi2_endcap2_65 = f_endcap2_42p5to56->GetChisquare();
//    Double_t chi2_endcap2_90 = f_endcap2_56to76->GetChisquare();
//    Double_t chi2_endcap2_5000 = f_endcap2_76to5000->GetChisquare();
//    cout << "Endcap fit chi squares: " << chi2_endcap2_45 << "  " << chi2_endcap2_65 << "  " << chi2_endcap2_90 << "  " << chi2_endcap2_5000 << endl;
//    cout << "\n\n";
    c_PR_endcap2->Update();


    TCanvas *c_PR_allin1 = new TCanvas("c_PR_allin1", "c_PR_allin1", 800, 800);
    c_PR_allin1->cd();
    c_PR_allin1->SetGrid(1);
    c_PR_allin1->SetRightMargin(0.05);
    c_PR_allin1->SetTopMargin(0.05);
    c_PR_allin1->SetBottomMargin(0.12);
    c_PR_allin1->SetLeftMargin(0.13);
    h_PR_ratio_endcap_draw->SetDirectory(0);
    h_PR_ratio_endcap_draw->SetTitle("");
    h_PR_ratio_endcap_draw->SetMarkerStyle(kFullSquare);
    h_PR_ratio_endcap_draw->SetMarkerColor(kBlue);
    h_PR_ratio_endcap_draw->SetLineColor(kBlue);
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitle("Prompt rate");
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_draw->GetYaxis()->SetLabelSize(0.04);
    h_PR_ratio_endcap_draw->GetYaxis()->SetTitleOffset(1.12);
    h_PR_ratio_endcap_draw->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_PR_ratio_endcap_draw->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_ratio_endcap_draw->GetXaxis()->SetNoExponent();
    h_PR_ratio_endcap_draw->GetXaxis()->SetMoreLogLabels();
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitleOffset(1);
    h_PR_ratio_endcap_draw->GetXaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_draw->GetXaxis()->SetLabelSize(0.04);
    h_PR_ratio_barrel->SetMarkerStyle(kFullDotLarge);
    h_PR_ratio_barrel->SetMarkerColor(kBlack);
    h_PR_ratio_barrel->SetLineColor(kBlack);
    h_PR_ratio_endcap2_draw->SetDirectory(0);
    h_PR_ratio_endcap2_draw->SetMarkerStyle(33);
    h_PR_ratio_endcap2_draw->SetMarkerSize(1.5);
    h_PR_ratio_endcap2_draw->SetMarkerColor(kOrange-3);
    h_PR_ratio_endcap2_draw->SetLineColor(kOrange-3);
    h_PR_ratio_barrel2_draw->SetDirectory(0);
    h_PR_ratio_barrel2_draw->SetMarkerStyle(23);
    h_PR_ratio_barrel2_draw->SetMarkerColor(kRed-3);
    h_PR_ratio_barrel2_draw->SetLineColor(kRed-3);
//    h_PR_ratio_endcap_draw->Draw();
//    h_PR_ratio_endcap2_draw->Draw("same");
    // Making a shifted histogram for easier viewing
    Double_t shifted_bins[10] = {52,63,72,83,93,105,160,213,540,1050};
    TH1D *h_PR_ratio_endcap_shifted = new TH1D("h_PR_ratio_endcap_shifted", "", 9, shifted_bins);
    for (Int_t i=1; i<=9; i++)
    {
        h_PR_ratio_endcap_shifted->SetBinContent(i, h_PR_ratio_endcap_draw->GetBinContent(i));
        h_PR_ratio_endcap_shifted->SetBinError(i, h_PR_ratio_endcap_draw->GetBinError(i));
    }
    h_PR_ratio_endcap_shifted->SetDirectory(0);
    h_PR_ratio_endcap_shifted->SetMarkerStyle(kFullSquare);
    h_PR_ratio_endcap_shifted->SetMarkerSize(1.4);
    h_PR_ratio_endcap_shifted->SetMarkerColor(kBlue);
    h_PR_ratio_endcap_shifted->SetLineColor(kBlue);
    h_PR_ratio_endcap_shifted->GetYaxis()->SetTitle("Prompt selection efficiency");
    h_PR_ratio_endcap_shifted->GetYaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_shifted->GetYaxis()->SetLabelSize(0.04);
    h_PR_ratio_endcap_shifted->GetYaxis()->SetTitleOffset(1.12);
    h_PR_ratio_endcap_shifted->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_PR_ratio_endcap_shifted->GetXaxis()->SetRangeUser(52, 500);
    h_PR_ratio_endcap_shifted->GetXaxis()->SetNoExponent();
    h_PR_ratio_endcap_shifted->GetXaxis()->SetMoreLogLabels();
    h_PR_ratio_endcap_shifted->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_ratio_endcap_shifted->GetXaxis()->SetTitleOffset(1);
    h_PR_ratio_endcap_shifted->GetXaxis()->SetTitleSize(0.05);
    h_PR_ratio_endcap_shifted->GetXaxis()->SetLabelSize(0.04);
    h_PR_ratio_endcap_shifted->Draw();

    h_PR_ratio_barrel->Draw("same");
//    h_PR_ratio_barrel2_draw->Draw("same");
    TLegend *legend_allin1 = new TLegend(0.13, 0.8, 0.4, 0.95);
//    legend_allin1->AddEntry(h_PR_ratio_barrel, "|#eta| < 0.7", "LP");
//    legend_allin1->AddEntry(h_PR_ratio_barrel2_draw, "0.7 #leq |#eta| < 1.2", "LP");
//    legend_allin1->AddEntry(h_PR_ratio_endcap_draw, "1.2 #leq |#eta| < 1.8", "LP");
//    legend_allin1->AddEntry(h_PR_ratio_endcap2_draw, "1.8 #leq |#eta| < 2.4", "LP");
    legend_allin1->AddEntry(h_PR_ratio_barrel, "|#eta| < 1.2", "LP");
    legend_allin1->AddEntry(h_PR_ratio_endcap_draw, "1.2 #leq |#eta| < 2.4", "LP");
    legend_allin1->Draw();
//    f_barrel_17to43->Draw("same");
//    f_barrel_43to56->Draw("same");
//    f_barrel_56to95->Draw("same");
//    f_barrel_95to5000->Draw("same");
//    f_barrel2_17to43->Draw("same");
//    f_barrel2_43to56->Draw("same");
//    f_barrel2_56to95->Draw("same");
//    f_barrel2_95to5000->Draw("same");
//    f_endcap_17to42p5->Draw("same");
//    f_endcap_42p5to55->Draw("same");
//    f_endcap_55to5000->Draw("same");
//    f_endcap2_17to42p5->Draw("same");
//    f_endcap2_42p5to56->Draw("same");
//    f_endcap2_56to76->Draw("same");
//    f_endcap2_76to5000->Draw("same");
    c_PR_allin1->SetLogx();
    c_PR_allin1->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "c_PR_eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
    TH1D *h_eff_MC_eta_draw   = ((TH1D*)(h_eff_MC_eta  ->Clone("h_PR_MC_eta_draw")));
    TH1D *h_PR_ratio_eta_draw = ((TH1D*)(h_PR_ratio_eta->Clone("h_PR_ratio_eta_draw")));
    h_eff_data_eta->SetDirectory(0);
    h_eff_MC_eta_draw->SetDirectory(0);
    h_PR_ratio_eta_draw->SetDirectory(0);
    h_eff_data_eta->SetTitle("");
    h_eff_data_eta->GetYaxis()->SetTitle("Prompt rate");
    h_eff_data_eta->GetYaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetYaxis()->SetLabelSize(0.04);
    h_eff_data_eta->GetYaxis()->SetTitleOffset(1.12);
    h_eff_data_eta->GetYaxis()->SetRangeUser(0.4, 1.1);
    h_eff_data_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_eff_data_eta->GetXaxis()->SetNoExponent();
    h_eff_data_eta->GetXaxis()->SetMoreLogLabels();
    h_eff_data_eta->GetXaxis()->SetTitle("#eta");
    h_eff_data_eta->GetXaxis()->SetTitleOffset(1);
    h_eff_data_eta->GetXaxis()->SetTitleSize(0.05);
    h_eff_data_eta->GetXaxis()->SetLabelSize(0.04);
    h_eff_data_eta->Draw();
    h_eff_MC_eta_draw->Draw("same");
    h_PR_ratio_eta_draw->Draw("same");
//    TLegend *legend_eta = new TLegend(0.13, 0.77, 0.6, 0.95);
//    legend_eta->AddEntry(h_eff_data_eta, "Subtraction", "LP");
//    legend_eta->Draw();
    legend1->Draw();
    c_PR_eta->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << inName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << inName << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/FR/Muon/PromptRate_muon_alt.root", "RECREATE");
    f_out->cd();
    h_eff_data_barrel    ->Write();
    h_eff_data_barrel2   ->Write();
    h_eff_data_endcap    ->Write();
    h_eff_data_endcap2   ->Write();
    h_eff_data_eta       ->Write();
    h_PR_subtract_barrel_plus ->Write();
    h_PR_subtract_barrel2_plus ->Write();
    h_PR_subtract_endcap_plus ->Write();
    h_PR_subtract_endcap2_plus ->Write();
    h_PR_subtract_barrel_minus ->Write();
    h_PR_subtract_barrel2_minus ->Write();
    h_PR_subtract_endcap_minus ->Write();
    h_PR_subtract_endcap2_minus->Write();
    h_eff_MC_barrel      ->Write();
    h_eff_MC_barrel2     ->Write();
    h_eff_MC_endcap      ->Write();
    h_eff_MC_endcap2     ->Write();
    h_eff_MC_eta         ->Write();
//    h_ineff_data_barrel  ->Write();
//    h_ineff_data_barrel2 ->Write();
//    h_ineff_data_endcap  ->Write();
//    h_ineff_data_endcap2 ->Write();
//    h_ineff_data_eta     ->Write();
//    h_ineff_MC_barrel    ->Write();
//    h_ineff_MC_barrel2   ->Write();
//    h_ineff_MC_endcap    ->Write();
//    h_ineff_MC_endcap2   ->Write();
//    h_ineff_MC_eta       ->Write();
//    h_eff_ratio_barrel   ->Write();
//    h_eff_ratio_barrel2  ->Write();
//    h_eff_ratio_endcap   ->Write();
//    h_eff_ratio_endcap2  ->Write();
//    h_eff_ratio_eta      ->Write();
//    h_ineff_ratio_barrel ->Write();
//    h_ineff_ratio_barrel2->Write();
//    h_ineff_ratio_endcap ->Write();
//    h_ineff_ratio_endcap2->Write();
//    h_ineff_ratio_eta    ->Write();
    h_PR_ratio_barrel    ->Write();
    h_PR_ratio_barrel2   ->Write();
    h_PR_ratio_endcap    ->Write();
    h_PR_ratio_endcap2   ->Write();
    h_PR_ratio_eta       ->Write();
//    f_barrel_17to43->Write();
//    f_barrel_43to56->Write();
//    f_barrel_56to95->Write();
//    f_barrel_95to5000->Write();
//    f_barrel2_17to43->Write();
//    f_barrel2_43to56->Write();
//    f_barrel2_56to95->Write();
//    f_barrel2_95to5000->Write();
//    f_endcap_17to42p5->Write();
//    f_endcap_42p5to55->Write();
//    f_endcap_55to5000->Write();
//    f_endcap2_17to42p5->Write();
//    f_endcap2_42p5to56->Write();
//    f_endcap2_56to76->Write();
//    f_endcap2_76to5000->Write();
    f_out->Close();

} // End of Mu_EstimatePR_alt


void Mu_EstFRandPR_MC()
{
    TString filename = "/media/sf_DATA/FR/Muon/MCFR_Hists_Mu.root";
    TFile *f = new TFile(filename, "READ");
    cout << "Reading from file " << filename << " ... ";

    TH1D *h_FR_pT_barrel_nume[6], *h_FR_pT_endcap_nume[6], *h_FR_pT_endcap2_nume[6];
    TH1D *h_FR_pT_barrel_deno[6], *h_FR_pT_endcap_deno[6], *h_FR_pT_endcap2_deno[6];
    TH1D *h_FR_eta_nume[6],       *h_FR_eta_deno[6];
    TH1D *h_FR_nJet_barrel_nume[6], *h_FR_nJet_endcap_nume[6], *h_FR_nJet_endcap2_nume[6];
    TH1D *h_FR_nJet_barrel_deno[6], *h_FR_nJet_endcap_deno[6], *h_FR_nJet_endcap2_deno[6];
    TH1D *h_FR_nBJet_barrel_nume[6], *h_FR_nBJet_endcap_nume[6], *h_FR_nBJet_endcap2_nume[6];
    TH1D *h_FR_nBJet_barrel_deno[6], *h_FR_nBJet_endcap_deno[6], *h_FR_nBJet_endcap2_deno[6];
    TH1D *h_FR_dRjet_barrel_nume[6], *h_FR_dRjet_endcap_nume[6], *h_FR_dRjet_endcap2_nume[6];
    TH1D *h_FR_dRjet_barrel_deno[6], *h_FR_dRjet_endcap_deno[6], *h_FR_dRjet_endcap2_deno[6];
    TH1D *h_FR_dPhiMET_barrel_nume[6], *h_FR_dPhiMET_endcap_nume[6], *h_FR_dPhiMET_endcap2_nume[6];
    TH1D *h_FR_dPhiMET_barrel_deno[6], *h_FR_dPhiMET_endcap_deno[6], *h_FR_dPhiMET_endcap2_deno[6];
    TH1D *h_FR_MET_barrel_nume[6], *h_FR_MET_endcap_nume[6], *h_FR_MET_endcap2_nume[6];
    TH1D *h_FR_MET_barrel_deno[6], *h_FR_MET_endcap_deno[6], *h_FR_MET_endcap2_deno[6];
    TH1D *h_PR_pT_barrel_nume[6], *h_PR_pT_endcap_nume[6], *h_PR_pT_endcap2_nume[6];
    TH1D *h_PR_pT_barrel_deno[6], *h_PR_pT_endcap_deno[6], *h_PR_pT_endcap2_deno[6];
    TH1D *h_PR_eta_nume[6],       *h_PR_eta_deno[6];
    TH1D *h_PR_nJet_barrel_nume[6], *h_PR_nJet_endcap_nume[6], *h_PR_nJet_endcap2_nume[6];
    TH1D *h_PR_nJet_barrel_deno[6], *h_PR_nJet_endcap_deno[6], *h_PR_nJet_endcap2_deno[6];
    TH1D *h_PR_nBJet_barrel_nume[6], *h_PR_nBJet_endcap_nume[6], *h_PR_nBJet_endcap2_nume[6];
    TH1D *h_PR_nBJet_barrel_deno[6], *h_PR_nBJet_endcap_deno[6], *h_PR_nBJet_endcap2_deno[6];
    TH1D *h_PR_dRjet_barrel_nume[6], *h_PR_dRjet_endcap_nume[6], *h_PR_dRjet_endcap2_nume[6];
    TH1D *h_PR_dRjet_barrel_deno[6], *h_PR_dRjet_endcap_deno[6], *h_PR_dRjet_endcap2_deno[6];
    TH1D *h_PR_dPhiMET_barrel_nume[6], *h_PR_dPhiMET_endcap_nume[6], *h_PR_dPhiMET_endcap2_nume[6];
    TH1D *h_PR_dPhiMET_barrel_deno[6], *h_PR_dPhiMET_endcap_deno[6], *h_PR_dPhiMET_endcap2_deno[6];
    TH1D *h_PR_MET_barrel_nume[6], *h_PR_MET_endcap_nume[6], *h_PR_MET_endcap2_nume[6];
    TH1D *h_PR_MET_barrel_deno[6], *h_PR_MET_endcap_deno[6], *h_PR_MET_endcap2_deno[6];
    TH1D *h_PR_pT_fromTop_barrel_nume[6], *h_PR_pT_fromTop_endcap_nume[6], *h_PR_pT_fromTop_endcap2_nume[6];
    TH1D *h_PR_pT_fromTop_barrel_deno[6], *h_PR_pT_fromTop_endcap_deno[6], *h_PR_pT_fromTop_endcap2_deno[6];
    TH1D *h_PR_pT_notFromTop_barrel_nume[6], *h_PR_pT_notFromTop_endcap_nume[6], *h_PR_pT_notFromTop_endcap2_nume[6];
    TH1D *h_PR_pT_notFromTop_barrel_deno[6], *h_PR_pT_notFromTop_endcap_deno[6], *h_PR_pT_notFromTop_endcap2_deno[6];
    TH1D *h_PR_eta_fromTop_nume[6],       *h_PR_eta_fromTop_deno[6];
    TH1D *h_PR_eta_notFromTop_nume[6],       *h_PR_eta_notFromTop_deno[6];
    TH1D *h_PR_nJet_fromTop_barrel_nume[6], *h_PR_nJet_fromTop_endcap_nume[6], *h_PR_nJet_fromTop_endcap2_nume[6];
    TH1D *h_PR_nJet_fromTop_barrel_deno[6], *h_PR_nJet_fromTop_endcap_deno[6], *h_PR_nJet_fromTop_endcap2_deno[6];
    TH1D *h_PR_nJet_notFromTop_barrel_nume[6], *h_PR_nJet_notFromTop_endcap_nume[6], *h_PR_nJet_notFromTop_endcap2_nume[6];
    TH1D *h_PR_nJet_notFromTop_barrel_deno[6], *h_PR_nJet_notFromTop_endcap_deno[6], *h_PR_nJet_notFromTop_endcap2_deno[6];
    TH1D *h_PR_MET_fromTop_barrel_nume[6],    *h_PR_MET_fromTop_endcap_nume[6],    *h_PR_MET_fromTop_endcap2_nume[6];
    TH1D *h_PR_MET_fromTop_barrel_deno[6],    *h_PR_MET_fromTop_endcap_deno[6],    *h_PR_MET_fromTop_endcap2_deno[6];
    TH1D *h_PR_MET_notFromTop_barrel_nume[6], *h_PR_MET_notFromTop_endcap_nume[6], *h_PR_MET_notFromTop_endcap2_nume[6];
    TH1D *h_PR_MET_notFromTop_barrel_deno[6], *h_PR_MET_notFromTop_endcap_deno[6], *h_PR_MET_notFromTop_endcap2_deno[6];
    TH1D *h_FR_nume_sum[6],       *h_FR_deno_sum[6],        *h_PR_nume_sum[6],       *h_PR_deno_sum[6];

    TString type[6] = {"DY", "ttbar", "tW", "diboson", "WJets", "QCD"};
    Color_t colors[6] = {kOrange, kCyan+2, kGreen+2, kMagenta-6, kRed-2, kRed+3};
    Style_t markers[6] = {kFullDotLarge, kFullSquare, kFullDiamond, kFullTriangleDown, kFullDotLarge, kFullTriangleUp};
    Float_t sizes[6] = {1, 1, 1.5, 1, 1, 1, };

    for (Int_t i=0; i<6; i++)
    {
        f->GetObject("h_FR_pT_barrel_nume_"+type[i],  h_FR_pT_barrel_nume[i]);
        f->GetObject("h_FR_pT_endcap_nume_"+type[i],  h_FR_pT_endcap_nume[i]);
        f->GetObject("h_FR_pT_endcap2_nume_"+type[i], h_FR_pT_endcap2_nume[i]);
        f->GetObject("h_FR_eta_nume_"+type[i],        h_FR_eta_nume[i]);
        f->GetObject("h_FR_nJet_barrel_nume_"+type[i],  h_FR_nJet_barrel_nume[i]);
        f->GetObject("h_FR_nJet_endcap_nume_"+type[i],  h_FR_nJet_endcap_nume[i]);
        f->GetObject("h_FR_nJet_endcap2_nume_"+type[i], h_FR_nJet_endcap2_nume[i]);
        f->GetObject("h_FR_nBJet_barrel_nume_"+type[i],  h_FR_nBJet_barrel_nume[i]);
        f->GetObject("h_FR_nBJet_endcap_nume_"+type[i],  h_FR_nBJet_endcap_nume[i]);
        f->GetObject("h_FR_nBJet_endcap2_nume_"+type[i], h_FR_nBJet_endcap2_nume[i]);
        f->GetObject("h_FR_dRjet_barrel_nume_"+type[i],  h_FR_dRjet_barrel_nume[i]);
        f->GetObject("h_FR_dRjet_endcap_nume_"+type[i],  h_FR_dRjet_endcap_nume[i]);
        f->GetObject("h_FR_dRjet_endcap2_nume_"+type[i], h_FR_dRjet_endcap2_nume[i]);
        f->GetObject("h_FR_dPhiMET_barrel_nume_"+type[i],  h_FR_dPhiMET_barrel_nume[i]);
        f->GetObject("h_FR_dPhiMET_endcap_nume_"+type[i],  h_FR_dPhiMET_endcap_nume[i]);
        f->GetObject("h_FR_dPhiMET_endcap2_nume_"+type[i], h_FR_dPhiMET_endcap2_nume[i]);
        f->GetObject("h_FR_MET_barrel_nume_"+type[i],  h_FR_MET_barrel_nume[i]);
        f->GetObject("h_FR_MET_endcap_nume_"+type[i],  h_FR_MET_endcap_nume[i]);
        f->GetObject("h_FR_MET_endcap2_nume_"+type[i], h_FR_MET_endcap2_nume[i]);
        f->GetObject("h_FR_pT_barrel_deno_"+type[i],  h_FR_pT_barrel_deno[i]);
        f->GetObject("h_FR_pT_endcap_deno_"+type[i],  h_FR_pT_endcap_deno[i]);
        f->GetObject("h_FR_pT_endcap2_deno_"+type[i], h_FR_pT_endcap2_deno[i]);
        f->GetObject("h_FR_eta_deno_"+type[i],        h_FR_eta_deno[i]);
        f->GetObject("h_FR_nJet_barrel_deno_"+type[i],  h_FR_nJet_barrel_deno[i]);
        f->GetObject("h_FR_nJet_endcap_deno_"+type[i],  h_FR_nJet_endcap_deno[i]);
        f->GetObject("h_FR_nJet_endcap2_deno_"+type[i], h_FR_nJet_endcap2_deno[i]);
        f->GetObject("h_FR_nBJet_barrel_deno_"+type[i],  h_FR_nBJet_barrel_deno[i]);
        f->GetObject("h_FR_nBJet_endcap_deno_"+type[i],  h_FR_nBJet_endcap_deno[i]);
        f->GetObject("h_FR_nBJet_endcap2_deno_"+type[i], h_FR_nBJet_endcap2_deno[i]);
        f->GetObject("h_FR_dRjet_barrel_deno_"+type[i],  h_FR_dRjet_barrel_deno[i]);
        f->GetObject("h_FR_dRjet_endcap_deno_"+type[i],  h_FR_dRjet_endcap_deno[i]);
        f->GetObject("h_FR_dRjet_endcap2_deno_"+type[i], h_FR_dRjet_endcap2_deno[i]);
        f->GetObject("h_FR_dPhiMET_barrel_deno_"+type[i],  h_FR_dPhiMET_barrel_deno[i]);
        f->GetObject("h_FR_dPhiMET_endcap_deno_"+type[i],  h_FR_dPhiMET_endcap_deno[i]);
        f->GetObject("h_FR_dPhiMET_endcap2_deno_"+type[i], h_FR_dPhiMET_endcap2_deno[i]);
        f->GetObject("h_FR_MET_barrel_deno_"+type[i],  h_FR_MET_barrel_deno[i]);
        f->GetObject("h_FR_MET_endcap_deno_"+type[i],  h_FR_MET_endcap_deno[i]);
        f->GetObject("h_FR_MET_endcap2_deno_"+type[i], h_FR_MET_endcap2_deno[i]);
        f->GetObject("h_PR_pT_barrel_nume_"+type[i],  h_PR_pT_barrel_nume[i]);
        f->GetObject("h_PR_pT_endcap_nume_"+type[i],  h_PR_pT_endcap_nume[i]);
        f->GetObject("h_PR_pT_endcap2_nume_"+type[i], h_PR_pT_endcap2_nume[i]);
        f->GetObject("h_PR_eta_nume_"+type[i],        h_PR_eta_nume[i]);
        f->GetObject("h_PR_nJet_barrel_nume_"+type[i],  h_PR_nJet_barrel_nume[i]);
        f->GetObject("h_PR_nJet_endcap_nume_"+type[i],  h_PR_nJet_endcap_nume[i]);
        f->GetObject("h_PR_nJet_endcap2_nume_"+type[i], h_PR_nJet_endcap2_nume[i]);
        f->GetObject("h_PR_nBJet_barrel_nume_"+type[i],  h_PR_nBJet_barrel_nume[i]);
        f->GetObject("h_PR_nBJet_endcap_nume_"+type[i],  h_PR_nBJet_endcap_nume[i]);
        f->GetObject("h_PR_nBJet_endcap2_nume_"+type[i], h_PR_nBJet_endcap2_nume[i]);
        f->GetObject("h_PR_dRjet_barrel_nume_"+type[i],  h_PR_dRjet_barrel_nume[i]);
        f->GetObject("h_PR_dRjet_endcap_nume_"+type[i],  h_PR_dRjet_endcap_nume[i]);
        f->GetObject("h_PR_dRjet_endcap2_nume_"+type[i], h_PR_dRjet_endcap2_nume[i]);
        f->GetObject("h_PR_dPhiMET_barrel_nume_"+type[i],  h_PR_dPhiMET_barrel_nume[i]);
        f->GetObject("h_PR_dPhiMET_endcap_nume_"+type[i],  h_PR_dPhiMET_endcap_nume[i]);
        f->GetObject("h_PR_dPhiMET_endcap2_nume_"+type[i], h_PR_dPhiMET_endcap2_nume[i]);
        f->GetObject("h_PR_MET_barrel_nume_"+type[i],  h_PR_MET_barrel_nume[i]);
        f->GetObject("h_PR_MET_endcap_nume_"+type[i],  h_PR_MET_endcap_nume[i]);
        f->GetObject("h_PR_MET_endcap2_nume_"+type[i], h_PR_MET_endcap2_nume[i]);
        f->GetObject("h_PR_pT_fromTop_barrel_nume_"+type[i],  h_PR_pT_fromTop_barrel_nume[i]);
        f->GetObject("h_PR_pT_fromTop_endcap_nume_"+type[i],  h_PR_pT_fromTop_endcap_nume[i]);
        f->GetObject("h_PR_pT_fromTop_endcap2_nume_"+type[i], h_PR_pT_fromTop_endcap2_nume[i]);
        f->GetObject("h_PR_pT_notFromTop_barrel_nume_"+type[i],  h_PR_pT_notFromTop_barrel_nume[i]);
        f->GetObject("h_PR_pT_notFromTop_endcap_nume_"+type[i],  h_PR_pT_notFromTop_endcap_nume[i]);
        f->GetObject("h_PR_pT_notFromTop_endcap2_nume_"+type[i], h_PR_pT_notFromTop_endcap2_nume[i]);
        f->GetObject("h_PR_eta_fromTop_nume_"+type[i],        h_PR_eta_fromTop_nume[i]);
        f->GetObject("h_PR_eta_notFromTop_nume_"+type[i],        h_PR_eta_notFromTop_nume[i]);
        f->GetObject("h_PR_nJet_fromTop_barrel_nume_"+type[i],  h_PR_nJet_fromTop_barrel_nume[i]);
        f->GetObject("h_PR_nJet_fromTop_endcap_nume_"+type[i],  h_PR_nJet_fromTop_endcap_nume[i]);
        f->GetObject("h_PR_nJet_fromTop_endcap2_nume_"+type[i], h_PR_nJet_fromTop_endcap2_nume[i]);
        f->GetObject("h_PR_nJet_notFromTop_barrel_nume_"+type[i],  h_PR_nJet_notFromTop_barrel_nume[i]);
        f->GetObject("h_PR_nJet_notFromTop_endcap_nume_"+type[i],  h_PR_nJet_notFromTop_endcap_nume[i]);
        f->GetObject("h_PR_nJet_notFromTop_endcap2_nume_"+type[i], h_PR_nJet_notFromTop_endcap2_nume[i]);
        f->GetObject("h_PR_MET_fromTop_barrel_nume_"+type[i],     h_PR_MET_fromTop_barrel_nume[i]);
        f->GetObject("h_PR_MET_fromTop_endcap_nume_"+type[i],     h_PR_MET_fromTop_endcap_nume[i]);
        f->GetObject("h_PR_MET_fromTop_endcap2_nume_"+type[i],    h_PR_MET_fromTop_endcap2_nume[i]);
        f->GetObject("h_PR_MET_notFromTop_barrel_nume_"+type[i],  h_PR_MET_notFromTop_barrel_nume[i]);
        f->GetObject("h_PR_MET_notFromTop_endcap_nume_"+type[i],  h_PR_MET_notFromTop_endcap_nume[i]);
        f->GetObject("h_PR_MET_notFromTop_endcap2_nume_"+type[i], h_PR_MET_notFromTop_endcap2_nume[i]);
        f->GetObject("h_PR_pT_barrel_deno_"+type[i],  h_PR_pT_barrel_deno[i]);
        f->GetObject("h_PR_pT_endcap_deno_"+type[i],  h_PR_pT_endcap_deno[i]);
        f->GetObject("h_PR_pT_endcap2_deno_"+type[i], h_PR_pT_endcap2_deno[i]);
        f->GetObject("h_PR_eta_deno_"+type[i],        h_PR_eta_deno[i]);
        f->GetObject("h_PR_nJet_barrel_deno_"+type[i],  h_PR_nJet_barrel_deno[i]);
        f->GetObject("h_PR_nJet_endcap_deno_"+type[i],  h_PR_nJet_endcap_deno[i]);
        f->GetObject("h_PR_nJet_endcap2_deno_"+type[i], h_PR_nJet_endcap2_deno[i]);
        f->GetObject("h_PR_nBJet_barrel_deno_"+type[i],  h_PR_nBJet_barrel_deno[i]);
        f->GetObject("h_PR_nBJet_endcap_deno_"+type[i],  h_PR_nBJet_endcap_deno[i]);
        f->GetObject("h_PR_nBJet_endcap2_deno_"+type[i], h_PR_nBJet_endcap2_deno[i]);
        f->GetObject("h_PR_dRjet_barrel_deno_"+type[i],  h_PR_dRjet_barrel_deno[i]);
        f->GetObject("h_PR_dRjet_endcap_deno_"+type[i],  h_PR_dRjet_endcap_deno[i]);
        f->GetObject("h_PR_dRjet_endcap2_deno_"+type[i], h_PR_dRjet_endcap2_deno[i]);
        f->GetObject("h_PR_dPhiMET_barrel_deno_"+type[i],  h_PR_dPhiMET_barrel_deno[i]);
        f->GetObject("h_PR_dPhiMET_endcap_deno_"+type[i],  h_PR_dPhiMET_endcap_deno[i]);
        f->GetObject("h_PR_dPhiMET_endcap2_deno_"+type[i], h_PR_dPhiMET_endcap2_deno[i]);
        f->GetObject("h_PR_MET_barrel_deno_"+type[i],  h_PR_MET_barrel_deno[i]);
        f->GetObject("h_PR_MET_endcap_deno_"+type[i],  h_PR_MET_endcap_deno[i]);
        f->GetObject("h_PR_MET_endcap2_deno_"+type[i], h_PR_MET_endcap2_deno[i]);
        f->GetObject("h_PR_pT_fromTop_barrel_deno_"+type[i],  h_PR_pT_fromTop_barrel_deno[i]);
        f->GetObject("h_PR_pT_fromTop_endcap_deno_"+type[i],  h_PR_pT_fromTop_endcap_deno[i]);
        f->GetObject("h_PR_pT_fromTop_endcap2_deno_"+type[i], h_PR_pT_fromTop_endcap2_deno[i]);
        f->GetObject("h_PR_pT_notFromTop_barrel_deno_"+type[i],  h_PR_pT_notFromTop_barrel_deno[i]);
        f->GetObject("h_PR_pT_notFromTop_endcap_deno_"+type[i],  h_PR_pT_notFromTop_endcap_deno[i]);
        f->GetObject("h_PR_pT_notFromTop_endcap2_deno_"+type[i], h_PR_pT_notFromTop_endcap2_deno[i]);
        f->GetObject("h_PR_eta_fromTop_deno_"+type[i],        h_PR_eta_fromTop_deno[i]);
        f->GetObject("h_PR_eta_notFromTop_deno_"+type[i],        h_PR_eta_notFromTop_deno[i]);
        f->GetObject("h_PR_nJet_fromTop_barrel_deno_"+type[i],  h_PR_nJet_fromTop_barrel_deno[i]);
        f->GetObject("h_PR_nJet_fromTop_endcap_deno_"+type[i],  h_PR_nJet_fromTop_endcap_deno[i]);
        f->GetObject("h_PR_nJet_fromTop_endcap2_deno_"+type[i], h_PR_nJet_fromTop_endcap2_deno[i]);
        f->GetObject("h_PR_nJet_notFromTop_barrel_deno_"+type[i],  h_PR_nJet_notFromTop_barrel_deno[i]);
        f->GetObject("h_PR_nJet_notFromTop_endcap_deno_"+type[i],  h_PR_nJet_notFromTop_endcap_deno[i]);
        f->GetObject("h_PR_nJet_notFromTop_endcap2_deno_"+type[i], h_PR_nJet_notFromTop_endcap2_deno[i]);
        f->GetObject("h_PR_MET_fromTop_barrel_deno_"+type[i],     h_PR_MET_fromTop_barrel_deno[i]);
        f->GetObject("h_PR_MET_fromTop_endcap_deno_"+type[i],     h_PR_MET_fromTop_endcap_deno[i]);
        f->GetObject("h_PR_MET_fromTop_endcap2_deno_"+type[i],    h_PR_MET_fromTop_endcap2_deno[i]);
        f->GetObject("h_PR_MET_notFromTop_barrel_deno_"+type[i],  h_PR_MET_notFromTop_barrel_deno[i]);
        f->GetObject("h_PR_MET_notFromTop_endcap_deno_"+type[i],  h_PR_MET_notFromTop_endcap_deno[i]);
        f->GetObject("h_PR_MET_notFromTop_endcap2_deno_"+type[i], h_PR_MET_notFromTop_endcap2_deno[i]);
        h_FR_pT_barrel_nume[i]->SetDirectory(0);
        h_FR_pT_endcap_nume[i]->SetDirectory(0);
        h_FR_pT_endcap2_nume[i]->SetDirectory(0);
        h_FR_eta_nume[i]->SetDirectory(0);
        h_FR_nJet_barrel_nume[i]->SetDirectory(0);
        h_FR_nJet_endcap_nume[i]->SetDirectory(0);
        h_FR_nJet_endcap2_nume[i]->SetDirectory(0);
        h_FR_nBJet_barrel_nume[i]->SetDirectory(0);
        h_FR_nBJet_endcap_nume[i]->SetDirectory(0);
        h_FR_nBJet_endcap2_nume[i]->SetDirectory(0);
        h_FR_dRjet_barrel_nume[i]->SetDirectory(0);
        h_FR_dRjet_endcap_nume[i]->SetDirectory(0);
        h_FR_dRjet_endcap2_nume[i]->SetDirectory(0);
        h_FR_dPhiMET_barrel_nume[i]->SetDirectory(0);
        h_FR_dPhiMET_endcap_nume[i]->SetDirectory(0);
        h_FR_dPhiMET_endcap2_nume[i]->SetDirectory(0);
        h_FR_MET_barrel_nume[i]->SetDirectory(0);
        h_FR_MET_endcap_nume[i]->SetDirectory(0);
        h_FR_MET_endcap2_nume[i]->SetDirectory(0);
        h_FR_pT_barrel_deno[i]->SetDirectory(0);
        h_FR_pT_endcap_deno[i]->SetDirectory(0);
        h_FR_pT_endcap2_deno[i]->SetDirectory(0);
        h_FR_eta_deno[i]->SetDirectory(0);
        h_FR_nJet_barrel_deno[i]->SetDirectory(0);
        h_FR_nJet_endcap_deno[i]->SetDirectory(0);
        h_FR_nJet_endcap2_deno[i]->SetDirectory(0);
        h_FR_nBJet_barrel_deno[i]->SetDirectory(0);
        h_FR_nBJet_endcap_deno[i]->SetDirectory(0);
        h_FR_nBJet_endcap2_deno[i]->SetDirectory(0);
        h_FR_dRjet_barrel_deno[i]->SetDirectory(0);
        h_FR_dRjet_endcap_deno[i]->SetDirectory(0);
        h_FR_dRjet_endcap2_deno[i]->SetDirectory(0);
        h_FR_dPhiMET_barrel_deno[i]->SetDirectory(0);
        h_FR_dPhiMET_endcap_deno[i]->SetDirectory(0);
        h_FR_dPhiMET_endcap2_deno[i]->SetDirectory(0);
        h_FR_MET_barrel_deno[i]->SetDirectory(0);
        h_FR_MET_endcap_deno[i]->SetDirectory(0);
        h_FR_MET_endcap2_deno[i]->SetDirectory(0);
        h_PR_pT_barrel_nume[i]->SetDirectory(0);
        h_PR_pT_endcap_nume[i]->SetDirectory(0);
        h_PR_pT_endcap2_nume[i]->SetDirectory(0);
        h_PR_eta_nume[i]->SetDirectory(0);
        h_PR_nJet_barrel_nume[i]->SetDirectory(0);
        h_PR_nJet_endcap_nume[i]->SetDirectory(0);
        h_PR_nJet_endcap2_nume[i]->SetDirectory(0);
        h_PR_nBJet_barrel_nume[i]->SetDirectory(0);
        h_PR_nBJet_endcap_nume[i]->SetDirectory(0);
        h_PR_nBJet_endcap2_nume[i]->SetDirectory(0);
        h_PR_dRjet_barrel_nume[i]->SetDirectory(0);
        h_PR_dRjet_endcap_nume[i]->SetDirectory(0);
        h_PR_dRjet_endcap2_nume[i]->SetDirectory(0);
        h_PR_dPhiMET_barrel_nume[i]->SetDirectory(0);
        h_PR_dPhiMET_endcap_nume[i]->SetDirectory(0);
        h_PR_dPhiMET_endcap2_nume[i]->SetDirectory(0);
        h_PR_MET_barrel_nume[i]->SetDirectory(0);
        h_PR_MET_endcap_nume[i]->SetDirectory(0);
        h_PR_MET_endcap2_nume[i]->SetDirectory(0);
        h_PR_pT_fromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_pT_fromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_pT_fromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_pT_notFromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_pT_notFromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_eta_fromTop_nume[i]->SetDirectory(0);
        h_PR_eta_notFromTop_nume[i]->SetDirectory(0);
        h_PR_nJet_fromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_nJet_fromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_MET_fromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_MET_fromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_MET_fromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_MET_notFromTop_barrel_nume[i]->SetDirectory(0);
        h_PR_MET_notFromTop_endcap_nume[i]->SetDirectory(0);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetDirectory(0);
        h_PR_pT_barrel_deno[i]->SetDirectory(0);
        h_PR_pT_endcap_deno[i]->SetDirectory(0);
        h_PR_pT_endcap2_deno[i]->SetDirectory(0);
        h_PR_eta_deno[i]->SetDirectory(0);
        h_PR_nJet_barrel_deno[i]->SetDirectory(0);
        h_PR_nJet_endcap_deno[i]->SetDirectory(0);
        h_PR_nJet_endcap2_deno[i]->SetDirectory(0);
        h_PR_nBJet_barrel_deno[i]->SetDirectory(0);
        h_PR_nBJet_endcap_deno[i]->SetDirectory(0);
        h_PR_nBJet_endcap2_deno[i]->SetDirectory(0);
        h_PR_dRjet_barrel_deno[i]->SetDirectory(0);
        h_PR_dRjet_endcap_deno[i]->SetDirectory(0);
        h_PR_dRjet_endcap2_deno[i]->SetDirectory(0);
        h_PR_dPhiMET_barrel_deno[i]->SetDirectory(0);
        h_PR_dPhiMET_endcap_deno[i]->SetDirectory(0);
        h_PR_dPhiMET_endcap2_deno[i]->SetDirectory(0);
        h_PR_MET_barrel_deno[i]->SetDirectory(0);
        h_PR_MET_endcap_deno[i]->SetDirectory(0);
        h_PR_MET_endcap2_deno[i]->SetDirectory(0);
        h_PR_pT_fromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_pT_fromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_pT_fromTop_endcap2_deno[i]->SetDirectory(0);
        h_PR_pT_notFromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_pT_notFromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_pT_notFromTop_endcap2_deno[i]->SetDirectory(0);
        h_PR_eta_fromTop_deno[i]->SetDirectory(0);
        h_PR_eta_notFromTop_deno[i]->SetDirectory(0);
        h_PR_nJet_fromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_nJet_fromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_nJet_fromTop_endcap2_deno[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_nJet_notFromTop_endcap2_deno[i]->SetDirectory(0);
        h_PR_MET_fromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_MET_fromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_MET_fromTop_endcap2_deno[i]->SetDirectory(0);
        h_PR_MET_notFromTop_barrel_deno[i]->SetDirectory(0);
        h_PR_MET_notFromTop_endcap_deno[i]->SetDirectory(0);
        h_PR_MET_notFromTop_endcap2_deno[i]->SetDirectory(0);

        removeNegativeBins(h_FR_pT_barrel_nume[i] );
        removeNegativeBins(h_FR_pT_endcap_nume[i] );
        removeNegativeBins(h_FR_pT_endcap2_nume[i]);
        removeNegativeBins(h_FR_eta_nume[i]       );
        removeNegativeBins(h_FR_nJet_barrel_nume[i] );
        removeNegativeBins(h_FR_nJet_endcap_nume[i] );
        removeNegativeBins(h_FR_nJet_endcap2_nume[i]);
        removeNegativeBins(h_FR_nBJet_barrel_nume[i] );
        removeNegativeBins(h_FR_nBJet_endcap_nume[i] );
        removeNegativeBins(h_FR_nBJet_endcap2_nume[i]);
        removeNegativeBins(h_FR_dRjet_barrel_nume[i] );
        removeNegativeBins(h_FR_dRjet_endcap_nume[i] );
        removeNegativeBins(h_FR_dRjet_endcap2_nume[i]);
        removeNegativeBins(h_FR_dPhiMET_barrel_nume[i] );
        removeNegativeBins(h_FR_dPhiMET_endcap_nume[i] );
        removeNegativeBins(h_FR_dPhiMET_endcap2_nume[i]);
        removeNegativeBins(h_FR_MET_barrel_nume[i] );
        removeNegativeBins(h_FR_MET_endcap_nume[i] );
        removeNegativeBins(h_FR_MET_endcap2_nume[i]);
        removeNegativeBins(h_PR_pT_barrel_nume[i] );
        removeNegativeBins(h_PR_pT_endcap_nume[i] );
        removeNegativeBins(h_PR_pT_endcap2_nume[i]);
        removeNegativeBins(h_PR_eta_nume[i]       );
        removeNegativeBins(h_PR_nJet_barrel_nume[i] );
        removeNegativeBins(h_PR_nJet_endcap_nume[i] );
        removeNegativeBins(h_PR_nJet_endcap2_nume[i]);
        removeNegativeBins(h_PR_nBJet_barrel_nume[i] );
        removeNegativeBins(h_PR_nBJet_endcap_nume[i] );
        removeNegativeBins(h_PR_nBJet_endcap2_nume[i]);
        removeNegativeBins(h_PR_dRjet_barrel_nume[i] );
        removeNegativeBins(h_PR_dRjet_endcap_nume[i] );
        removeNegativeBins(h_PR_dRjet_endcap2_nume[i]);
        removeNegativeBins(h_PR_dPhiMET_barrel_nume[i] );
        removeNegativeBins(h_PR_dPhiMET_endcap_nume[i] );
        removeNegativeBins(h_PR_dPhiMET_endcap2_nume[i]);
        removeNegativeBins(h_PR_MET_barrel_nume[i] );
        removeNegativeBins(h_PR_MET_endcap_nume[i] );
        removeNegativeBins(h_PR_MET_endcap2_nume[i]);
        removeNegativeBins(h_PR_pT_fromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_pT_fromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_pT_fromTop_endcap2_nume[i]);
        removeNegativeBins(h_PR_pT_notFromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_pT_notFromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_pT_notFromTop_endcap2_nume[i]);
        removeNegativeBins(h_PR_eta_fromTop_nume[i]       );
        removeNegativeBins(h_PR_eta_notFromTop_nume[i]       );
        removeNegativeBins(h_PR_nJet_fromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_nJet_fromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_nJet_fromTop_endcap2_nume[i]);
        removeNegativeBins(h_PR_nJet_notFromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_nJet_notFromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_nJet_notFromTop_endcap2_nume[i]);
        removeNegativeBins(h_PR_MET_fromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_MET_fromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_MET_fromTop_endcap2_nume[i]);
        removeNegativeBins(h_PR_MET_notFromTop_barrel_nume[i] );
        removeNegativeBins(h_PR_MET_notFromTop_endcap_nume[i] );
        removeNegativeBins(h_PR_MET_notFromTop_endcap2_nume[i]);

        h_FR_pT_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_pT_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_pT_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_eta_nume[i]       ->SetLineColor(colors[i]);
        h_FR_nJet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_nJet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_nJet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_nBJet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_nBJet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_nBJet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_dRjet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_dRjet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_dRjet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_dPhiMET_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_dPhiMET_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_dPhiMET_endcap2_nume[i]->SetLineColor(colors[i]);
        h_FR_MET_barrel_nume[i] ->SetLineColor(colors[i]);
        h_FR_MET_endcap_nume[i] ->SetLineColor(colors[i]);
        h_FR_MET_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_pT_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_eta_nume[i]       ->SetLineColor(colors[i]);
        h_PR_nJet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_nBJet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_nBJet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_nBJet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_dRjet_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_dRjet_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_dRjet_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_dPhiMET_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_dPhiMET_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_dPhiMET_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_MET_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_pT_fromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_fromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_fromTop_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_pT_notFromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_notFromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_eta_fromTop_nume[i]       ->SetLineColor(colors[i]);
        h_PR_eta_notFromTop_nume[i]       ->SetLineColor(colors[i]);
        h_PR_nJet_fromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_fromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_nJet_notFromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_notFromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_MET_fromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_fromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_fromTop_endcap2_nume[i]->SetLineColor(colors[i]);
        h_PR_MET_notFromTop_barrel_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_notFromTop_endcap_nume[i] ->SetLineColor(colors[i]);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetLineColor(colors[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_eta_nume[i]       ->SetMarkerColor(colors[i]);
        h_FR_nJet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_nJet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_nJet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_nBJet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_nBJet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_nBJet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_dRjet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_dRjet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_dRjet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_dPhiMET_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_dPhiMET_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_dPhiMET_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_FR_MET_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_MET_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_FR_MET_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_eta_nume[i]       ->SetMarkerColor(colors[i]);
        h_PR_nJet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_nBJet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nBJet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nBJet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_dRjet_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_dRjet_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_dRjet_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_dPhiMET_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_dPhiMET_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_dPhiMET_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_MET_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_pT_fromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_fromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_fromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_pT_notFromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_notFromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_eta_fromTop_nume[i]       ->SetMarkerColor(colors[i]);
        h_PR_eta_notFromTop_nume[i]       ->SetMarkerColor(colors[i]);
        h_PR_nJet_fromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_fromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_nJet_notFromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_notFromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_MET_fromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_fromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_fromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);
        h_PR_MET_notFromTop_barrel_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_notFromTop_endcap_nume[i] ->SetMarkerColor(colors[i]);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetMarkerColor(colors[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_eta_nume[i]       ->SetMarkerStyle(markers[i]);
        h_FR_nJet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_nJet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_nJet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_nBJet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_nBJet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_nBJet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_dRjet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_dRjet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_dRjet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_dPhiMET_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_dPhiMET_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_dPhiMET_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_FR_MET_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_MET_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_FR_MET_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_eta_nume[i]       ->SetMarkerStyle(markers[i]);
        h_PR_nJet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_nBJet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nBJet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nBJet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_dRjet_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_dRjet_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_dRjet_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_dPhiMET_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_dPhiMET_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_dPhiMET_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_MET_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_pT_fromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_fromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_fromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_pT_notFromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_notFromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_eta_fromTop_nume[i]       ->SetMarkerStyle(markers[i]);
        h_PR_eta_notFromTop_nume[i]       ->SetMarkerStyle(markers[i]);
        h_PR_nJet_fromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_fromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_nJet_notFromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_notFromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_MET_fromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_fromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_fromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);
        h_PR_MET_notFromTop_barrel_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_notFromTop_endcap_nume[i] ->SetMarkerStyle(markers[i]);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetMarkerStyle(markers[i]);

        h_FR_pT_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_pT_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_pT_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_eta_nume[i]       ->SetMarkerSize(sizes[i]);
        h_FR_nJet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_nJet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_nJet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_nBJet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_nBJet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_nBJet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_dRjet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_dRjet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_dRjet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_dPhiMET_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_dPhiMET_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_dPhiMET_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_FR_MET_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_MET_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_FR_MET_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_pT_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_eta_nume[i]       ->SetMarkerSize(sizes[i]);
        h_PR_nJet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_nBJet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nBJet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nBJet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_dRjet_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_dRjet_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_dRjet_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_dPhiMET_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_dPhiMET_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_dPhiMET_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_MET_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_pT_fromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_fromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_fromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_pT_notFromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_notFromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_eta_fromTop_nume[i]       ->SetMarkerSize(sizes[i]);
        h_PR_eta_notFromTop_nume[i]       ->SetMarkerSize(sizes[i]);
        h_PR_nJet_fromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_fromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_nJet_notFromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_notFromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_MET_fromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_fromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_fromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);
        h_PR_MET_notFromTop_barrel_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_notFromTop_endcap_nume[i] ->SetMarkerSize(sizes[i]);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetMarkerSize(sizes[i]);

        h_FR_pT_barrel_nume[i] ->SetStats(0);
        h_FR_pT_endcap_nume[i] ->SetStats(0);
        h_FR_pT_endcap2_nume[i]->SetStats(0);
        h_FR_eta_nume[i]       ->SetStats(0);
        h_FR_nJet_barrel_nume[i] ->SetStats(0);
        h_FR_nJet_endcap_nume[i] ->SetStats(0);
        h_FR_nJet_endcap2_nume[i]->SetStats(0);
        h_FR_nBJet_barrel_nume[i] ->SetStats(0);
        h_FR_nBJet_endcap_nume[i] ->SetStats(0);
        h_FR_nBJet_endcap2_nume[i]->SetStats(0);
        h_FR_dRjet_barrel_nume[i] ->SetStats(0);
        h_FR_dRjet_endcap_nume[i] ->SetStats(0);
        h_FR_dRjet_endcap2_nume[i]->SetStats(0);
        h_FR_dPhiMET_barrel_nume[i] ->SetStats(0);
        h_FR_dPhiMET_endcap_nume[i] ->SetStats(0);
        h_FR_dPhiMET_endcap2_nume[i]->SetStats(0);
        h_FR_MET_barrel_nume[i] ->SetStats(0);
        h_FR_MET_endcap_nume[i] ->SetStats(0);
        h_FR_MET_endcap2_nume[i]->SetStats(0);
        h_PR_pT_barrel_nume[i] ->SetStats(0);
        h_PR_pT_endcap_nume[i] ->SetStats(0);
        h_PR_pT_endcap2_nume[i]->SetStats(0);
        h_PR_eta_nume[i]       ->SetStats(0);
        h_PR_nJet_barrel_nume[i] ->SetStats(0);
        h_PR_nJet_endcap_nume[i] ->SetStats(0);
        h_PR_nJet_endcap2_nume[i]->SetStats(0);
        h_PR_nBJet_barrel_nume[i] ->SetStats(0);
        h_PR_nBJet_endcap_nume[i] ->SetStats(0);
        h_PR_nBJet_endcap2_nume[i]->SetStats(0);
        h_PR_dRjet_barrel_nume[i] ->SetStats(0);
        h_PR_dRjet_endcap_nume[i] ->SetStats(0);
        h_PR_dRjet_endcap2_nume[i]->SetStats(0);
        h_PR_dPhiMET_barrel_nume[i] ->SetStats(0);
        h_PR_dPhiMET_endcap_nume[i] ->SetStats(0);
        h_PR_dPhiMET_endcap2_nume[i]->SetStats(0);
        h_PR_MET_barrel_nume[i] ->SetStats(0);
        h_PR_MET_endcap_nume[i] ->SetStats(0);
        h_PR_MET_endcap2_nume[i]->SetStats(0);
        h_PR_pT_fromTop_barrel_nume[i] ->SetStats(0);
        h_PR_pT_fromTop_endcap_nume[i] ->SetStats(0);
        h_PR_pT_fromTop_endcap2_nume[i]->SetStats(0);
        h_PR_pT_notFromTop_barrel_nume[i] ->SetStats(0);
        h_PR_pT_notFromTop_endcap_nume[i] ->SetStats(0);
        h_PR_pT_notFromTop_endcap2_nume[i]->SetStats(0);
        h_PR_eta_fromTop_nume[i]       ->SetStats(0);
        h_PR_eta_fromTop_nume[i]       ->SetStats(0);
        h_PR_nJet_fromTop_barrel_nume[i] ->SetStats(0);
        h_PR_nJet_fromTop_endcap_nume[i] ->SetStats(0);
        h_PR_nJet_fromTop_endcap2_nume[i]->SetStats(0);
        h_PR_nJet_notFromTop_barrel_nume[i] ->SetStats(0);
        h_PR_nJet_notFromTop_endcap_nume[i] ->SetStats(0);
        h_PR_nJet_notFromTop_endcap2_nume[i]->SetStats(0);
        h_PR_MET_fromTop_barrel_nume[i] ->SetStats(0);
        h_PR_MET_fromTop_endcap_nume[i] ->SetStats(0);
        h_PR_MET_fromTop_endcap2_nume[i]->SetStats(0);
        h_PR_MET_notFromTop_barrel_nume[i] ->SetStats(0);
        h_PR_MET_notFromTop_endcap_nume[i] ->SetStats(0);
        h_PR_MET_notFromTop_endcap2_nume[i]->SetStats(0);

        // Calculating rates
        // Average rates
        h_FR_nume_sum[i] = ((TH1D*)(h_FR_pT_barrel_nume[i]->Clone("h_FR_nume_sum_"+type[i])));
        h_FR_nume_sum[i]->SetDirectory(0);
        h_FR_nume_sum[i]->Add(h_FR_pT_endcap_nume[i]);
        h_FR_nume_sum[i]->Add(h_FR_pT_endcap2_nume[i]);
        h_FR_deno_sum[i] = ((TH1D*)(h_FR_pT_barrel_deno[i]->Clone("h_FR_deno_sum_"+type[i])));
        h_FR_deno_sum[i]->SetDirectory(0);
        h_FR_deno_sum[i]->Add(h_FR_pT_endcap_deno[i]);
        h_FR_deno_sum[i]->Add(h_FR_pT_endcap2_deno[i]);

        h_PR_nume_sum[i] = ((TH1D*)(h_PR_pT_barrel_nume[i]->Clone("h_FR_nume_sum_"+type[i])));
        h_PR_nume_sum[i]->SetDirectory(0);
        h_PR_nume_sum[i]->Add(h_PR_pT_endcap_nume[i]);
        h_PR_nume_sum[i]->Add(h_PR_pT_endcap2_nume[i]);
        h_PR_deno_sum[i] = ((TH1D*)(h_PR_pT_barrel_deno[i]->Clone("h_FR_deno_sum_"+type[i])));
        h_PR_deno_sum[i]->SetDirectory(0);
        h_PR_deno_sum[i]->Add(h_PR_pT_endcap_deno[i]);
        h_PR_deno_sum[i]->Add(h_PR_pT_endcap2_deno[i]);

        // Full rates
        h_FR_pT_barrel_nume[i] ->Divide(h_FR_pT_barrel_deno[i] );
        h_FR_pT_endcap_nume[i] ->Divide(h_FR_pT_endcap_deno[i] );
        h_FR_pT_endcap2_nume[i]->Divide(h_FR_pT_endcap2_deno[i]);
        h_FR_eta_nume[i]       ->Divide(h_FR_eta_deno[i]       );
        h_FR_nJet_barrel_nume[i] ->Divide(h_FR_nJet_barrel_deno[i] );
        h_FR_nJet_endcap_nume[i] ->Divide(h_FR_nJet_endcap_deno[i] );
        h_FR_nJet_endcap2_nume[i]->Divide(h_FR_nJet_endcap2_deno[i]);
        h_FR_nBJet_barrel_nume[i] ->Divide(h_FR_nBJet_barrel_deno[i] );
        h_FR_nBJet_endcap_nume[i] ->Divide(h_FR_nBJet_endcap_deno[i] );
        h_FR_nBJet_endcap2_nume[i]->Divide(h_FR_nBJet_endcap2_deno[i]);
        h_FR_dRjet_barrel_nume[i] ->Divide(h_FR_dRjet_barrel_deno[i] );
        h_FR_dRjet_endcap_nume[i] ->Divide(h_FR_dRjet_endcap_deno[i] );
        h_FR_dRjet_endcap2_nume[i]->Divide(h_FR_dRjet_endcap2_deno[i]);
        h_FR_dPhiMET_barrel_nume[i] ->Divide(h_FR_dPhiMET_barrel_deno[i] );
        h_FR_dPhiMET_endcap_nume[i] ->Divide(h_FR_dPhiMET_endcap_deno[i] );
        h_FR_dPhiMET_endcap2_nume[i]->Divide(h_FR_dPhiMET_endcap2_deno[i]);
        h_FR_MET_barrel_nume[i] ->Divide(h_FR_MET_barrel_deno[i] );
        h_FR_MET_endcap_nume[i] ->Divide(h_FR_MET_endcap_deno[i] );
        h_FR_MET_endcap2_nume[i]->Divide(h_FR_MET_endcap2_deno[i]);
        h_PR_pT_barrel_nume[i] ->Divide(h_PR_pT_barrel_deno[i] );
        h_PR_pT_endcap_nume[i] ->Divide(h_PR_pT_endcap_deno[i] );
        h_PR_pT_endcap2_nume[i]->Divide(h_PR_pT_endcap2_deno[i]);
        h_PR_eta_nume[i]       ->Divide(h_PR_eta_deno[i]       );
        h_PR_nJet_barrel_nume[i] ->Divide(h_PR_nJet_barrel_deno[i] );
        h_PR_nJet_endcap_nume[i] ->Divide(h_PR_nJet_endcap_deno[i] );
        h_PR_nJet_endcap2_nume[i]->Divide(h_PR_nJet_endcap2_deno[i]);
        h_PR_nBJet_barrel_nume[i] ->Divide(h_PR_nBJet_barrel_deno[i] );
        h_PR_nBJet_endcap_nume[i] ->Divide(h_PR_nBJet_endcap_deno[i] );
        h_PR_nBJet_endcap2_nume[i]->Divide(h_PR_nBJet_endcap2_deno[i]);
        h_PR_dRjet_barrel_nume[i] ->Divide(h_PR_dRjet_barrel_deno[i] );
        h_PR_dRjet_endcap_nume[i] ->Divide(h_PR_dRjet_endcap_deno[i] );
        h_PR_dRjet_endcap2_nume[i]->Divide(h_PR_dRjet_endcap2_deno[i]);
        h_PR_dPhiMET_barrel_nume[i] ->Divide(h_PR_dPhiMET_barrel_deno[i] );
        h_PR_dPhiMET_endcap_nume[i] ->Divide(h_PR_dPhiMET_endcap_deno[i] );
        h_PR_dPhiMET_endcap2_nume[i]->Divide(h_PR_dPhiMET_endcap2_deno[i]);
        h_PR_MET_barrel_nume[i] ->Divide(h_PR_MET_barrel_deno[i] );
        h_PR_MET_endcap_nume[i] ->Divide(h_PR_MET_endcap_deno[i] );
        h_PR_MET_endcap2_nume[i]->Divide(h_PR_MET_endcap2_deno[i]);
        h_PR_pT_fromTop_barrel_nume[i] ->Divide(h_PR_pT_fromTop_barrel_deno[i] );
        h_PR_pT_fromTop_endcap_nume[i] ->Divide(h_PR_pT_fromTop_endcap_deno[i] );
        h_PR_pT_fromTop_endcap2_nume[i]->Divide(h_PR_pT_fromTop_endcap2_deno[i]);
        h_PR_pT_notFromTop_barrel_nume[i] ->Divide(h_PR_pT_notFromTop_barrel_deno[i] );
        h_PR_pT_notFromTop_endcap_nume[i] ->Divide(h_PR_pT_notFromTop_endcap_deno[i] );
        h_PR_pT_notFromTop_endcap2_nume[i]->Divide(h_PR_pT_notFromTop_endcap2_deno[i]);
        h_PR_eta_fromTop_nume[i]          ->Divide(h_PR_eta_fromTop_deno[i]          );
        h_PR_eta_notFromTop_nume[i]       ->Divide(h_PR_eta_notFromTop_deno[i]       );
        h_PR_nJet_fromTop_barrel_nume[i] ->Divide(h_PR_nJet_fromTop_barrel_deno[i] );
        h_PR_nJet_fromTop_endcap_nume[i] ->Divide(h_PR_nJet_fromTop_endcap_deno[i] );
        h_PR_nJet_fromTop_endcap2_nume[i]->Divide(h_PR_nJet_fromTop_endcap2_deno[i]);
        h_PR_nJet_notFromTop_barrel_nume[i] ->Divide(h_PR_nJet_notFromTop_barrel_deno[i] );
        h_PR_nJet_notFromTop_endcap_nume[i] ->Divide(h_PR_nJet_notFromTop_endcap_deno[i] );
        h_PR_nJet_notFromTop_endcap2_nume[i]->Divide(h_PR_nJet_notFromTop_endcap2_deno[i]);
        h_PR_MET_fromTop_barrel_nume[i]    ->Divide(h_PR_MET_fromTop_barrel_deno[i] );
        h_PR_MET_fromTop_endcap_nume[i]    ->Divide(h_PR_MET_fromTop_endcap_deno[i] );
        h_PR_MET_fromTop_endcap2_nume[i]   ->Divide(h_PR_MET_fromTop_endcap2_deno[i]);
        h_PR_MET_notFromTop_barrel_nume[i] ->Divide(h_PR_MET_notFromTop_barrel_deno[i] );
        h_PR_MET_notFromTop_endcap_nume[i] ->Divide(h_PR_MET_notFromTop_endcap_deno[i] );
        h_PR_MET_notFromTop_endcap2_nume[i]->Divide(h_PR_MET_notFromTop_endcap2_deno[i]);
    }

    cout << "finished." << endl;
    f->Close();
    if (!f->IsOpen()) cout << "File has been closed successfully.\n" << endl;
    else cout << "FILE " << filename << " COULD NOT BE CLOSED!\n" << endl;

    // Drawing
    TLegend *legend_FR = new TLegend(0.7, 0.75, 0.95, 0.95);
    legend_FR->AddEntry(h_FR_pT_barrel_nume[1], "t#bar{t}", "lp");
    legend_FR->AddEntry(h_FR_pT_barrel_nume[2], "tW", "lp");
    legend_FR->AddEntry(h_FR_pT_barrel_nume[5], "QCD", "lp");

    TLegend *legend_PR = new TLegend(0.5, 0.78, 0.95, 0.95);
    legend_PR->SetNColumns(2);
    legend_PR->AddEntry(h_FR_pT_barrel_nume[0], "DY", "lp");
    legend_PR->AddEntry(h_FR_pT_barrel_nume[1], "t#bar{t}", "lp");
    legend_PR->AddEntry(h_FR_pT_barrel_nume[2], "tW", "lp");
    legend_PR->AddEntry(h_FR_pT_barrel_nume[3], "diboson", "lp");
    legend_PR->AddEntry(h_FR_pT_barrel_nume[4], "W+Jets", "lp");

    TLegend *legend_fromTop = new TLegend(0.7, 0.75, 0.95, 0.95);
    legend_fromTop->AddEntry(h_FR_pT_barrel_nume[1], "t#bar{t}", "lp");
    legend_fromTop->AddEntry(h_FR_pT_barrel_nume[2], "tW", "lp");

    TLegend *legend_notFromTop = new TLegend(0.6, 0.85, 0.95, 0.95);
    legend_notFromTop->SetNColumns(2);
    legend_notFromTop->AddEntry(h_FR_pT_barrel_nume[0], "DY", "lp");
    legend_notFromTop->AddEntry(h_FR_pT_barrel_nume[2], "tW", "lp");
    legend_notFromTop->AddEntry(h_FR_pT_barrel_nume[3], "diboson", "lp");
    legend_notFromTop->AddEntry(h_FR_pT_barrel_nume[4], "W+Jets", "lp");

    TCanvas *c_FR_barrel  = new TCanvas("c_FR_barrel",  "FR barrel", 800, 800);
    c_FR_barrel->cd();
    c_FR_barrel->SetGrid(1);
    c_FR_barrel->SetLogx(1);
    c_FR_barrel->SetRightMargin(0.05);
    c_FR_barrel->SetTopMargin(0.05);
    c_FR_barrel->SetBottomMargin(0.12);
    c_FR_barrel->SetLeftMargin(0.13);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_barrel_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_pT_barrel_nume[1]->Draw();
    h_FR_pT_barrel_nume[2]->Draw("same");
    h_FR_pT_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_barrel->Update();

    TCanvas *c_FR_endcap  = new TCanvas("c_FR_endcap",  "FR endcap", 800, 800);
    c_FR_endcap->cd();
    c_FR_endcap->SetGrid(1);
    c_FR_endcap->SetLogx(1);
    c_FR_endcap->SetRightMargin(0.05);
    c_FR_endcap->SetTopMargin(0.05);
    c_FR_endcap->SetBottomMargin(0.12);
    c_FR_endcap->SetLeftMargin(0.13);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_endcap_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_pT_endcap_nume[1]->Draw();
    h_FR_pT_endcap_nume[2]->Draw("same");
    h_FR_pT_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_endcap->Update();

    TCanvas *c_FR_endcap2 = new TCanvas("c_FR_endcap2", "FR endcap2", 800, 800);
    c_FR_endcap2->cd();
    c_FR_endcap2->SetGrid(1);
    c_FR_endcap2->SetLogx(1);
    c_FR_endcap2->SetRightMargin(0.05);
    c_FR_endcap2->SetTopMargin(0.05);
    c_FR_endcap2->SetBottomMargin(0.12);
    c_FR_endcap2->SetLeftMargin(0.13);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_pT_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_pT_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_pT_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_pT_endcap2_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_FR_pT_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_pT_endcap2_nume[1]->Draw();
    h_FR_pT_endcap2_nume[2]->Draw("same");
    h_FR_pT_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_endcap2->Update();

    TCanvas *c_FR_eta = new TCanvas("c_FR_eta", "FR eta", 800, 800);
    c_FR_eta->cd();
    c_FR_eta->SetGrid(1);
    c_FR_eta->SetRightMargin(0.05);
    c_FR_eta->SetTopMargin(0.05);
    c_FR_eta->SetBottomMargin(0.12);
    c_FR_eta->SetLeftMargin(0.13);
    h_FR_eta_nume[1]->GetXaxis()->SetTitle("#eta (#mu)");
    h_FR_eta_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_eta_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_eta_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_eta_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_eta_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_eta_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_eta_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_eta_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_eta_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_eta_nume[1]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_FR_eta_nume[1]->GetYaxis()->SetRangeUser(0, 0.7);
    h_FR_eta_nume[1]->Draw();
    h_FR_eta_nume[2]->Draw("same");
    h_FR_eta_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_eta->Update();

    TCanvas *c_FR_nJet_barrel  = new TCanvas("c_FR_nJet_barrel",  "FR(nJet) barrel", 800, 800);
    c_FR_nJet_barrel->cd();
    c_FR_nJet_barrel->SetGrid(1);
    c_FR_nJet_barrel->SetRightMargin(0.05);
    c_FR_nJet_barrel->SetTopMargin(0.05);
    c_FR_nJet_barrel->SetBottomMargin(0.12);
    c_FR_nJet_barrel->SetLeftMargin(0.13);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nJet_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nJet_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nJet_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nJet_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nJet_barrel_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_FR_nJet_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_nJet_barrel_nume[1]->Draw();
    h_FR_nJet_barrel_nume[2]->Draw("same");
    h_FR_nJet_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nJet_barrel->Update();

    TCanvas *c_FR_nJet_endcap  = new TCanvas("c_FR_nJet_endcap",  "FR(nJet) endcap", 800, 800);
    c_FR_nJet_endcap->cd();
    c_FR_nJet_endcap->SetGrid(1);
    c_FR_nJet_endcap->SetRightMargin(0.05);
    c_FR_nJet_endcap->SetTopMargin(0.05);
    c_FR_nJet_endcap->SetBottomMargin(0.12);
    c_FR_nJet_endcap->SetLeftMargin(0.13);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nJet_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nJet_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nJet_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nJet_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nJet_endcap_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_FR_nJet_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_nJet_endcap_nume[1]->Draw();
    h_FR_nJet_endcap_nume[2]->Draw("same");
    h_FR_nJet_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nJet_endcap->Update();

    TCanvas *c_FR_nJet_endcap2 = new TCanvas("c_FR_nJet_endcap2", "FR(nJet) endcap2", 800, 800);
    c_FR_nJet_endcap2->cd();
    c_FR_nJet_endcap2->SetGrid(1);
    c_FR_nJet_endcap2->SetRightMargin(0.05);
    c_FR_nJet_endcap2->SetTopMargin(0.05);
    c_FR_nJet_endcap2->SetBottomMargin(0.12);
    c_FR_nJet_endcap2->SetLeftMargin(0.13);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nJet_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nJet_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nJet_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nJet_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nJet_endcap2_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_FR_nJet_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_nJet_endcap2_nume[1]->Draw();
    h_FR_nJet_endcap2_nume[2]->Draw("same");
    h_FR_nJet_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nJet_endcap2->Update();       

    TCanvas *c_FR_nBJet_barrel  = new TCanvas("c_FR_nBJet_barrel",  "FR(nBJet) barrel", 800, 800);
    c_FR_nBJet_barrel->cd();
    c_FR_nBJet_barrel->SetGrid(1);
    c_FR_nBJet_barrel->SetRightMargin(0.05);
    c_FR_nBJet_barrel->SetTopMargin(0.05);
    c_FR_nBJet_barrel->SetBottomMargin(0.12);
    c_FR_nBJet_barrel->SetLeftMargin(0.13);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetTitle("N_{b-jets}");
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nBJet_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nBJet_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nBJet_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nBJet_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nBJet_barrel_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_FR_nBJet_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_nBJet_barrel_nume[1]->Draw();
    h_FR_nBJet_barrel_nume[2]->Draw("same");
    h_FR_nBJet_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nBJet_barrel->Update();

    TCanvas *c_FR_nBJet_endcap  = new TCanvas("c_FR_nBJet_endcap",  "FR(nBJet) endcap", 800, 800);
    c_FR_nBJet_endcap->cd();
    c_FR_nBJet_endcap->SetGrid(1);
    c_FR_nBJet_endcap->SetRightMargin(0.05);
    c_FR_nBJet_endcap->SetTopMargin(0.05);
    c_FR_nBJet_endcap->SetBottomMargin(0.12);
    c_FR_nBJet_endcap->SetLeftMargin(0.13);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetTitle("N_{b-jets}");
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nBJet_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nBJet_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nBJet_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nBJet_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nBJet_endcap_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_FR_nBJet_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_nBJet_endcap_nume[1]->Draw();
    h_FR_nBJet_endcap_nume[2]->Draw("same");
    h_FR_nBJet_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nBJet_endcap->Update();

    TCanvas *c_FR_nBJet_endcap2 = new TCanvas("c_FR_nBJet_endcap2", "FR(nBJet) endcap2", 800, 800);
    c_FR_nBJet_endcap2->cd();
    c_FR_nBJet_endcap2->SetGrid(1);
    c_FR_nBJet_endcap2->SetRightMargin(0.05);
    c_FR_nBJet_endcap2->SetTopMargin(0.05);
    c_FR_nBJet_endcap2->SetBottomMargin(0.12);
    c_FR_nBJet_endcap2->SetLeftMargin(0.13);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetTitle("N_{b-jets}");
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_nBJet_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_nBJet_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_nBJet_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_nBJet_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_nBJet_endcap2_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_FR_nBJet_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_nBJet_endcap2_nume[1]->Draw();
    h_FR_nBJet_endcap2_nume[2]->Draw("same");
    h_FR_nBJet_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_nBJet_endcap2->Update();

    TCanvas *c_FR_dRjet_barrel  = new TCanvas("c_FR_dRjet_barrel",  "FR(dRjet) barrel", 800, 800);
    c_FR_dRjet_barrel->cd();
    c_FR_dRjet_barrel->SetGrid(1);
    c_FR_dRjet_barrel->SetRightMargin(0.05);
    c_FR_dRjet_barrel->SetTopMargin(0.05);
    c_FR_dRjet_barrel->SetBottomMargin(0.12);
    c_FR_dRjet_barrel->SetLeftMargin(0.13);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dRjet_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dRjet_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dRjet_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dRjet_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dRjet_barrel_nume[1]->GetXaxis()->SetRangeUser(0, 5);
    h_FR_dRjet_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_dRjet_barrel_nume[1]->Draw();
    h_FR_dRjet_barrel_nume[2]->Draw("same");
    h_FR_dRjet_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dRjet_barrel->Update();

    TCanvas *c_FR_dRjet_endcap  = new TCanvas("c_FR_dRjet_endcap",  "FR(dRjet) endcap", 800, 800);
    c_FR_dRjet_endcap->cd();
    c_FR_dRjet_endcap->SetGrid(1);
    c_FR_dRjet_endcap->SetRightMargin(0.05);
    c_FR_dRjet_endcap->SetTopMargin(0.05);
    c_FR_dRjet_endcap->SetBottomMargin(0.12);
    c_FR_dRjet_endcap->SetLeftMargin(0.13);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dRjet_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dRjet_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dRjet_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dRjet_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dRjet_endcap_nume[1]->GetXaxis()->SetRangeUser(0, 5);
    h_FR_dRjet_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_dRjet_endcap_nume[1]->Draw();
    h_FR_dRjet_endcap_nume[2]->Draw("same");
    h_FR_dRjet_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dRjet_endcap->Update();

    TCanvas *c_FR_dRjet_endcap2 = new TCanvas("c_FR_dRjet_endcap2", "FR(dRjet) endcap2", 800, 800);
    c_FR_dRjet_endcap2->cd();
    c_FR_dRjet_endcap2->SetGrid(1);
    c_FR_dRjet_endcap2->SetRightMargin(0.05);
    c_FR_dRjet_endcap2->SetTopMargin(0.05);
    c_FR_dRjet_endcap2->SetBottomMargin(0.12);
    c_FR_dRjet_endcap2->SetLeftMargin(0.13);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dRjet_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dRjet_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dRjet_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dRjet_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dRjet_endcap2_nume[1]->GetXaxis()->SetRangeUser(0, 5);
    h_FR_dRjet_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_dRjet_endcap2_nume[1]->Draw();
    h_FR_dRjet_endcap2_nume[2]->Draw("same");
    h_FR_dRjet_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dRjet_endcap2->Update();

    TCanvas *c_FR_dPhiMET_barrel  = new TCanvas("c_FR_dPhiMET_barrel",  "FR(dPhiMET) barrel", 800, 800);
    c_FR_dPhiMET_barrel->cd();
    c_FR_dPhiMET_barrel->SetGrid(1);
    c_FR_dPhiMET_barrel->SetRightMargin(0.05);
    c_FR_dPhiMET_barrel->SetTopMargin(0.05);
    c_FR_dPhiMET_barrel->SetBottomMargin(0.12);
    c_FR_dPhiMET_barrel->SetLeftMargin(0.13);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetTitle("#Delta $phi_{#mu, MET}");
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dPhiMET_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dPhiMET_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dPhiMET_barrel_nume[1]->GetXaxis()->SetRangeUser(0, 6.3);
    h_FR_dPhiMET_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_dPhiMET_barrel_nume[1]->Draw();
    h_FR_dPhiMET_barrel_nume[2]->Draw("same");
    h_FR_dPhiMET_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dPhiMET_barrel->Update();

    TCanvas *c_FR_dPhiMET_endcap  = new TCanvas("c_FR_dPhiMET_endcap",  "FR(dPhiMET) endcap", 800, 800);
    c_FR_dPhiMET_endcap->cd();
    c_FR_dPhiMET_endcap->SetGrid(1);
    c_FR_dPhiMET_endcap->SetRightMargin(0.05);
    c_FR_dPhiMET_endcap->SetTopMargin(0.05);
    c_FR_dPhiMET_endcap->SetBottomMargin(0.12);
    c_FR_dPhiMET_endcap->SetLeftMargin(0.13);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetTitle("#Delta $phi_{#mu, MET}");
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dPhiMET_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dPhiMET_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dPhiMET_endcap_nume[1]->GetXaxis()->SetRangeUser(0, 6.3);
    h_FR_dPhiMET_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_dPhiMET_endcap_nume[1]->Draw();
    h_FR_dPhiMET_endcap_nume[2]->Draw("same");
    h_FR_dPhiMET_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dPhiMET_endcap->Update();

    TCanvas *c_FR_dPhiMET_endcap2 = new TCanvas("c_FR_dPhiMET_endcap2", "FR(dPhiMET) endcap2", 800, 800);
    c_FR_dPhiMET_endcap2->cd();
    c_FR_dPhiMET_endcap2->SetGrid(1);
    c_FR_dPhiMET_endcap2->SetRightMargin(0.05);
    c_FR_dPhiMET_endcap2->SetTopMargin(0.05);
    c_FR_dPhiMET_endcap2->SetBottomMargin(0.12);
    c_FR_dPhiMET_endcap2->SetLeftMargin(0.13);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetTitle("#Delta $phi_{#mu, MET}");
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_dPhiMET_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_dPhiMET_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_dPhiMET_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_dPhiMET_endcap2_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_FR_dPhiMET_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_dPhiMET_endcap2_nume[1]->Draw();
    h_FR_dPhiMET_endcap2_nume[2]->Draw("same");
    h_FR_dPhiMET_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_dPhiMET_endcap2->Update();

    TCanvas *c_FR_MET_barrel  = new TCanvas("c_FR_MET_barrel",  "FR(MET) barrel", 800, 800);
    c_FR_MET_barrel->cd();
    c_FR_MET_barrel->SetGrid(1);
    c_FR_MET_barrel->SetRightMargin(0.05);
    c_FR_MET_barrel->SetTopMargin(0.05);
    c_FR_MET_barrel->SetBottomMargin(0.12);
    c_FR_MET_barrel->SetLeftMargin(0.13);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_MET_barrel_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_MET_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_MET_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_MET_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_MET_barrel_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_FR_MET_barrel_nume[1]->GetYaxis()->SetRangeUser(0, 0.24);
    h_FR_MET_barrel_nume[1]->Draw();
    h_FR_MET_barrel_nume[2]->Draw("same");
    h_FR_MET_barrel_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_MET_barrel->Update();

    TCanvas *c_FR_MET_endcap  = new TCanvas("c_FR_MET_endcap",  "FR(MET) endcap", 800, 800);
    c_FR_MET_endcap->cd();
    c_FR_MET_endcap->SetGrid(1);
    c_FR_MET_endcap->SetRightMargin(0.05);
    c_FR_MET_endcap->SetTopMargin(0.05);
    c_FR_MET_endcap->SetBottomMargin(0.12);
    c_FR_MET_endcap->SetLeftMargin(0.13);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_MET_endcap_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_MET_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_MET_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_MET_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_MET_endcap_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_FR_MET_endcap_nume[1]->GetYaxis()->SetRangeUser(0, 0.26);
    h_FR_MET_endcap_nume[1]->Draw();
    h_FR_MET_endcap_nume[2]->Draw("same");
    h_FR_MET_endcap_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_MET_endcap->Update();

    TCanvas *c_FR_MET_endcap2 = new TCanvas("c_FR_MET_endcap2", "FR(MET) endcap2", 800, 800);
    c_FR_MET_endcap2->cd();
    c_FR_MET_endcap2->SetGrid(1);
    c_FR_MET_endcap2->SetRightMargin(0.05);
    c_FR_MET_endcap2->SetTopMargin(0.05);
    c_FR_MET_endcap2->SetBottomMargin(0.12);
    c_FR_MET_endcap2->SetLeftMargin(0.13);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_FR_MET_endcap2_nume[1]->GetYaxis()->SetTitle("Fake rate");
    h_FR_MET_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_FR_MET_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_FR_MET_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_FR_MET_endcap2_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_FR_MET_endcap2_nume[1]->GetYaxis()->SetRangeUser(0, 0.64);
    h_FR_MET_endcap2_nume[1]->Draw();
    h_FR_MET_endcap2_nume[2]->Draw("same");
    h_FR_MET_endcap2_nume[5]->Draw("same");
    legend_FR->Draw();
    c_FR_MET_endcap2->Update();

    TCanvas *c_PR_barrel  = new TCanvas("c_PR_barrel",  "PR barrel", 800, 800);
    c_PR_barrel->cd();
    c_PR_barrel->SetGrid(1);
    c_PR_barrel->SetLogx(1);
    c_PR_barrel->SetRightMargin(0.05);
    c_PR_barrel->SetTopMargin(0.05);
    c_PR_barrel->SetBottomMargin(0.12);
    c_PR_barrel->SetLeftMargin(0.13);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_barrel_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_barrel_nume[0]->Draw();
    h_PR_pT_barrel_nume[1]->Draw("same");
    h_PR_pT_barrel_nume[2]->Draw("same");
    h_PR_pT_barrel_nume[3]->Draw("same");
    h_PR_pT_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_barrel->Update();

    TCanvas *c_PR_endcap  = new TCanvas("c_PR_endcap",  "PR endcap", 800, 800);
    c_PR_endcap->cd();
    c_PR_endcap->SetGrid(1);
    c_PR_endcap->SetLogx(1);
    c_PR_endcap->SetRightMargin(0.05);
    c_PR_endcap->SetTopMargin(0.05);
    c_PR_endcap->SetBottomMargin(0.12);
    c_PR_endcap->SetLeftMargin(0.13);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_endcap_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_endcap_nume[0]->Draw();
    h_PR_pT_endcap_nume[1]->Draw("same");
    h_PR_pT_endcap_nume[2]->Draw("same");
    h_PR_pT_endcap_nume[3]->Draw("same");
    h_PR_pT_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_endcap->Update();

    TCanvas *c_PR_endcap2 = new TCanvas("c_PR_endcap2", "PR endcap2", 800, 800);
    c_PR_endcap2->cd();
    c_PR_endcap2->SetGrid(1);
    c_PR_endcap2->SetLogx(1);
    c_PR_endcap2->SetRightMargin(0.05);
    c_PR_endcap2->SetTopMargin(0.05);
    c_PR_endcap2->SetBottomMargin(0.12);
    c_PR_endcap2->SetLeftMargin(0.13);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_endcap2_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_endcap2_nume[0]->Draw();
    h_PR_pT_endcap2_nume[1]->Draw("same");
    h_PR_pT_endcap2_nume[2]->Draw("same");
    h_PR_pT_endcap2_nume[3]->Draw("same");
    h_PR_pT_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_endcap2->Update();

    TCanvas *c_PR_eta = new TCanvas("c_PR_eta", "PR eta", 800, 800);
    c_PR_eta->cd();
    c_PR_eta->SetGrid(1);
    c_PR_eta->SetRightMargin(0.05);
    c_PR_eta->SetTopMargin(0.05);
    c_PR_eta->SetBottomMargin(0.12);
    c_PR_eta->SetLeftMargin(0.13);
    h_PR_eta_nume[0]->GetXaxis()->SetTitle("#eta (#mu)");
    h_PR_eta_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_eta_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_eta_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_eta_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_eta_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_eta_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_eta_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_eta_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_eta_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_eta_nume[0]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_PR_eta_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_eta_nume[0]->Draw();
    h_PR_eta_nume[1]->Draw("same");
    h_PR_eta_nume[2]->Draw("same");
    h_PR_eta_nume[3]->Draw("same");
    h_PR_eta_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_eta->Update();    

    TCanvas *c_PR_nJet_barrel  = new TCanvas("c_PR_nJet_barrel",  "PR(nJet) barrel", 800, 800);
    c_PR_nJet_barrel->cd();
    c_PR_nJet_barrel->SetGrid(1);
    c_PR_nJet_barrel->SetRightMargin(0.05);
    c_PR_nJet_barrel->SetTopMargin(0.05);
    c_PR_nJet_barrel->SetBottomMargin(0.12);
    c_PR_nJet_barrel->SetLeftMargin(0.13);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_barrel_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_PR_nJet_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_barrel_nume[0]->Draw();
    h_PR_nJet_barrel_nume[1]->Draw("same");
    h_PR_nJet_barrel_nume[2]->Draw("same");
    h_PR_nJet_barrel_nume[3]->Draw("same");
    h_PR_nJet_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nJet_barrel->Update();

    TCanvas *c_PR_nJet_endcap  = new TCanvas("c_PR_nJet_endcap",  "PR(nJet) endcap", 800, 800);
    c_PR_nJet_endcap->cd();
    c_PR_nJet_endcap->SetGrid(1);
    c_PR_nJet_endcap->SetRightMargin(0.05);
    c_PR_nJet_endcap->SetTopMargin(0.05);
    c_PR_nJet_endcap->SetBottomMargin(0.12);
    c_PR_nJet_endcap->SetLeftMargin(0.13);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_endcap_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_PR_nJet_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_endcap_nume[0]->Draw();
    h_PR_nJet_endcap_nume[1]->Draw("same");
    h_PR_nJet_endcap_nume[2]->Draw("same");
    h_PR_nJet_endcap_nume[3]->Draw("same");
    h_PR_nJet_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nJet_endcap->Update();

    TCanvas *c_PR_nJet_endcap2 = new TCanvas("c_PR_nJet_endcap2", "PR(nJet) endcap2", 800, 800);
    c_PR_nJet_endcap2->cd();
    c_PR_nJet_endcap2->SetGrid(1);
    c_PR_nJet_endcap2->SetRightMargin(0.05);
    c_PR_nJet_endcap2->SetTopMargin(0.05);
    c_PR_nJet_endcap2->SetBottomMargin(0.12);
    c_PR_nJet_endcap2->SetLeftMargin(0.13);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_endcap2_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 10-0.5);
    h_PR_nJet_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_endcap2_nume[0]->Draw();
    h_PR_nJet_endcap2_nume[1]->Draw("same");
    h_PR_nJet_endcap2_nume[2]->Draw("same");
    h_PR_nJet_endcap2_nume[3]->Draw("same");
    h_PR_nJet_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nJet_endcap2->Update();

    TCanvas *c_PR_nBJet_barrel  = new TCanvas("c_PR_nBJet_barrel",  "PR(nBJet) barrel", 800, 800);
    c_PR_nBJet_barrel->cd();
    c_PR_nBJet_barrel->SetGrid(1);
    c_PR_nBJet_barrel->SetRightMargin(0.05);
    c_PR_nBJet_barrel->SetTopMargin(0.05);
    c_PR_nBJet_barrel->SetBottomMargin(0.12);
    c_PR_nBJet_barrel->SetLeftMargin(0.13);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetTitle("N_{b-jets}");
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nBJet_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nBJet_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nBJet_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nBJet_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nBJet_barrel_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_PR_nBJet_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nBJet_barrel_nume[0]->Draw();
    h_PR_nBJet_barrel_nume[1]->Draw("same");
    h_PR_nBJet_barrel_nume[2]->Draw("same");
    h_PR_nBJet_barrel_nume[3]->Draw("same");
    h_PR_nBJet_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nBJet_barrel->Update();

    TCanvas *c_PR_nBJet_endcap  = new TCanvas("c_PR_nBJet_endcap",  "PR(nBJet) endcap", 800, 800);
    c_PR_nBJet_endcap->cd();
    c_PR_nBJet_endcap->SetGrid(1);
    c_PR_nBJet_endcap->SetRightMargin(0.05);
    c_PR_nBJet_endcap->SetTopMargin(0.05);
    c_PR_nBJet_endcap->SetBottomMargin(0.12);
    c_PR_nBJet_endcap->SetLeftMargin(0.13);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetTitle("N_{b-jets}");
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nBJet_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nBJet_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nBJet_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nBJet_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nBJet_endcap_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_PR_nBJet_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nBJet_endcap_nume[0]->Draw();
    h_PR_nBJet_endcap_nume[1]->Draw("same");
    h_PR_nBJet_endcap_nume[2]->Draw("same");
    h_PR_nBJet_endcap_nume[3]->Draw("same");
    h_PR_nBJet_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nBJet_endcap->Update();

    TCanvas *c_PR_nBJet_endcap2 = new TCanvas("c_PR_nBJet_endcap2", "PR(nBJet) endcap2", 800, 800);
    c_PR_nBJet_endcap2->cd();
    c_PR_nBJet_endcap2->SetGrid(1);
    c_PR_nBJet_endcap2->SetRightMargin(0.05);
    c_PR_nBJet_endcap2->SetTopMargin(0.05);
    c_PR_nBJet_endcap2->SetBottomMargin(0.12);
    c_PR_nBJet_endcap2->SetLeftMargin(0.13);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetTitle("N_{b-jets}");
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nBJet_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nBJet_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nBJet_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nBJet_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nBJet_endcap2_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 6-0.5);
    h_PR_nBJet_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nBJet_endcap2_nume[0]->Draw();
    h_PR_nBJet_endcap2_nume[1]->Draw("same");
    h_PR_nBJet_endcap2_nume[2]->Draw("same");
    h_PR_nBJet_endcap2_nume[3]->Draw("same");
    h_PR_nBJet_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_nBJet_endcap2->Update();

    TCanvas *c_PR_dRjet_barrel  = new TCanvas("c_PR_dRjet_barrel",  "PR(dRjet) barrel", 800, 800);
    c_PR_dRjet_barrel->cd();
    c_PR_dRjet_barrel->SetGrid(1);
    c_PR_dRjet_barrel->SetRightMargin(0.05);
    c_PR_dRjet_barrel->SetTopMargin(0.05);
    c_PR_dRjet_barrel->SetBottomMargin(0.12);
    c_PR_dRjet_barrel->SetLeftMargin(0.13);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dRjet_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dRjet_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dRjet_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dRjet_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dRjet_barrel_nume[0]->GetXaxis()->SetRangeUser(0, 5);
    h_PR_dRjet_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dRjet_barrel_nume[0]->Draw();
    h_PR_dRjet_barrel_nume[1]->Draw("same");
    h_PR_dRjet_barrel_nume[2]->Draw("same");
    h_PR_dRjet_barrel_nume[3]->Draw("same");
    h_PR_dRjet_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dRjet_barrel->Update();

    TCanvas *c_PR_dRjet_endcap  = new TCanvas("c_PR_dRjet_endcap",  "PR(dRjet) endcap", 800, 800);
    c_PR_dRjet_endcap->cd();
    c_PR_dRjet_endcap->SetGrid(1);
    c_PR_dRjet_endcap->SetRightMargin(0.05);
    c_PR_dRjet_endcap->SetTopMargin(0.05);
    c_PR_dRjet_endcap->SetBottomMargin(0.12);
    c_PR_dRjet_endcap->SetLeftMargin(0.13);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dRjet_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dRjet_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dRjet_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dRjet_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dRjet_endcap_nume[0]->GetXaxis()->SetRangeUser(0, 5);
    h_PR_dRjet_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dRjet_endcap_nume[0]->Draw();
    h_PR_dRjet_endcap_nume[1]->Draw("same");
    h_PR_dRjet_endcap_nume[2]->Draw("same");
    h_PR_dRjet_endcap_nume[3]->Draw("same");
    h_PR_dRjet_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dRjet_endcap->Update();

    TCanvas *c_PR_dRjet_endcap2 = new TCanvas("c_PR_dRjet_endcap2", "PR(dRjet) endcap2", 800, 800);
    c_PR_dRjet_endcap2->cd();
    c_PR_dRjet_endcap2->SetGrid(1);
    c_PR_dRjet_endcap2->SetRightMargin(0.05);
    c_PR_dRjet_endcap2->SetTopMargin(0.05);
    c_PR_dRjet_endcap2->SetBottomMargin(0.12);
    c_PR_dRjet_endcap2->SetLeftMargin(0.13);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetTitle("#Delta R_{#mu, j}");
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dRjet_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dRjet_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dRjet_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dRjet_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dRjet_endcap2_nume[0]->GetXaxis()->SetRangeUser(0, 5);
    h_PR_dRjet_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dRjet_endcap2_nume[0]->Draw();
    h_PR_dRjet_endcap2_nume[1]->Draw("same");
    h_PR_dRjet_endcap2_nume[2]->Draw("same");
    h_PR_dRjet_endcap2_nume[3]->Draw("same");
    h_PR_dRjet_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dRjet_endcap2->Update();

    TCanvas *c_PR_dPhiMET_barrel  = new TCanvas("c_PR_dPhiMET_barrel",  "PR(dPhiMET) barrel", 800, 800);
    c_PR_dPhiMET_barrel->cd();
    c_PR_dPhiMET_barrel->SetGrid(1);
    c_PR_dPhiMET_barrel->SetRightMargin(0.05);
    c_PR_dPhiMET_barrel->SetTopMargin(0.05);
    c_PR_dPhiMET_barrel->SetBottomMargin(0.12);
    c_PR_dPhiMET_barrel->SetLeftMargin(0.13);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetTitle("#Delta #phi_{#mu, MET}");
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dPhiMET_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dPhiMET_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dPhiMET_barrel_nume[0]->GetXaxis()->SetRangeUser(0, 6.3);
    h_PR_dPhiMET_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dPhiMET_barrel_nume[0]->Draw();
    h_PR_dPhiMET_barrel_nume[1]->Draw("same");
    h_PR_dPhiMET_barrel_nume[2]->Draw("same");
    h_PR_dPhiMET_barrel_nume[3]->Draw("same");
    h_PR_dPhiMET_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dPhiMET_barrel->Update();

    TCanvas *c_PR_dPhiMET_endcap  = new TCanvas("c_PR_dPhiMET_endcap",  "PR(dPhiMET) endcap", 800, 800);
    c_PR_dPhiMET_endcap->cd();
    c_PR_dPhiMET_endcap->SetGrid(1);
    c_PR_dPhiMET_endcap->SetRightMargin(0.05);
    c_PR_dPhiMET_endcap->SetTopMargin(0.05);
    c_PR_dPhiMET_endcap->SetBottomMargin(0.12);
    c_PR_dPhiMET_endcap->SetLeftMargin(0.13);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetTitle("#Delta #phi_{#mu, MET}");
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dPhiMET_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dPhiMET_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dPhiMET_endcap_nume[0]->GetXaxis()->SetRangeUser(0, 6.3);
    h_PR_dPhiMET_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dPhiMET_endcap_nume[0]->Draw();
    h_PR_dPhiMET_endcap_nume[1]->Draw("same");
    h_PR_dPhiMET_endcap_nume[2]->Draw("same");
    h_PR_dPhiMET_endcap_nume[3]->Draw("same");
    h_PR_dPhiMET_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dPhiMET_endcap->Update();

    TCanvas *c_PR_dPhiMET_endcap2 = new TCanvas("c_PR_dPhiMET_endcap2", "PR(dPhiMET) endcap2", 800, 800);
    c_PR_dPhiMET_endcap2->cd();
    c_PR_dPhiMET_endcap2->SetGrid(1);
    c_PR_dPhiMET_endcap2->SetRightMargin(0.05);
    c_PR_dPhiMET_endcap2->SetTopMargin(0.05);
    c_PR_dPhiMET_endcap2->SetBottomMargin(0.12);
    c_PR_dPhiMET_endcap2->SetLeftMargin(0.13);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetTitle("#Delta #phi_{#mu, MET}");
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_dPhiMET_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_dPhiMET_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_dPhiMET_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_dPhiMET_endcap2_nume[0]->GetXaxis()->SetRangeUser(0, 6.3);
    h_PR_dPhiMET_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_dPhiMET_endcap2_nume[0]->Draw();
    h_PR_dPhiMET_endcap2_nume[1]->Draw("same");
    h_PR_dPhiMET_endcap2_nume[2]->Draw("same");
    h_PR_dPhiMET_endcap2_nume[3]->Draw("same");
    h_PR_dPhiMET_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_dPhiMET_endcap2->Update();

    TCanvas *c_PR_MET_barrel  = new TCanvas("c_PR_MET_barrel",  "PR(MET) barrel", 800, 800);
    c_PR_MET_barrel->cd();
    c_PR_MET_barrel->SetGrid(1);
    c_PR_MET_barrel->SetRightMargin(0.05);
    c_PR_MET_barrel->SetTopMargin(0.05);
    c_PR_MET_barrel->SetBottomMargin(0.12);
    c_PR_MET_barrel->SetLeftMargin(0.13);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_barrel_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_barrel_nume[0]->Draw();
    h_PR_MET_barrel_nume[1]->Draw("same");
    h_PR_MET_barrel_nume[2]->Draw("same");
    h_PR_MET_barrel_nume[3]->Draw("same");
    h_PR_MET_barrel_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_MET_barrel->Update();

    TCanvas *c_PR_MET_endcap  = new TCanvas("c_PR_MET_endcap",  "PR(MET) endcap", 800, 800);
    c_PR_MET_endcap->cd();
    c_PR_MET_endcap->SetGrid(1);
    c_PR_MET_endcap->SetRightMargin(0.05);
    c_PR_MET_endcap->SetTopMargin(0.05);
    c_PR_MET_endcap->SetBottomMargin(0.12);
    c_PR_MET_endcap->SetLeftMargin(0.13);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_endcap_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_endcap_nume[0]->Draw();
    h_PR_MET_endcap_nume[1]->Draw("same");
    h_PR_MET_endcap_nume[2]->Draw("same");
    h_PR_MET_endcap_nume[3]->Draw("same");
    h_PR_MET_endcap_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_MET_endcap->Update();

    TCanvas *c_PR_MET_endcap2 = new TCanvas("c_PR_MET_endcap2", "PR(MET) endcap2", 800, 800);
    c_PR_MET_endcap2->cd();
    c_PR_MET_endcap2->SetGrid(1);
    c_PR_MET_endcap2->SetRightMargin(0.05);
    c_PR_MET_endcap2->SetTopMargin(0.05);
    c_PR_MET_endcap2->SetBottomMargin(0.12);
    c_PR_MET_endcap2->SetLeftMargin(0.13);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_endcap2_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_endcap2_nume[0]->Draw();
    h_PR_MET_endcap2_nume[1]->Draw("same");
    h_PR_MET_endcap2_nume[2]->Draw("same");
    h_PR_MET_endcap2_nume[3]->Draw("same");
    h_PR_MET_endcap2_nume[4]->Draw("same");
    legend_PR->Draw();
    c_PR_MET_endcap2->Update();

    TCanvas *c_PR_fromTop_barrel  = new TCanvas("c_PR_fromTop_barrel",  "PR (from Top) barrel", 800, 800);
    c_PR_fromTop_barrel->cd();
    c_PR_fromTop_barrel->SetGrid(1);
    c_PR_fromTop_barrel->SetLogx(1);
    c_PR_fromTop_barrel->SetRightMargin(0.05);
    c_PR_fromTop_barrel->SetTopMargin(0.05);
    c_PR_fromTop_barrel->SetBottomMargin(0.12);
    c_PR_fromTop_barrel->SetLeftMargin(0.13);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_barrel_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_fromTop_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_fromTop_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_fromTop_barrel_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_fromTop_barrel_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_fromTop_barrel_nume[1]->Draw();
    h_PR_pT_fromTop_barrel_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_fromTop_barrel->Update();

    TCanvas *c_PR_fromTop_endcap  = new TCanvas("c_PR_fromTop_endcap",  "PR (from top) endcap", 800, 800);
    c_PR_fromTop_endcap->cd();
    c_PR_fromTop_endcap->SetGrid(1);
    c_PR_fromTop_endcap->SetLogx(1);
    c_PR_fromTop_endcap->SetRightMargin(0.05);
    c_PR_fromTop_endcap->SetTopMargin(0.05);
    c_PR_fromTop_endcap->SetBottomMargin(0.12);
    c_PR_fromTop_endcap->SetLeftMargin(0.13);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_endcap_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_fromTop_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_fromTop_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_fromTop_endcap_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_fromTop_endcap_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_fromTop_endcap_nume[1]->Draw();
    h_PR_pT_fromTop_endcap_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_fromTop_endcap->Update();

    TCanvas *c_PR_fromTop_endcap2 = new TCanvas("c_PR_fromTop_endcap2", "PR (from top) endcap2", 800, 800);
    c_PR_fromTop_endcap2->cd();
    c_PR_fromTop_endcap2->SetGrid(1);
    c_PR_fromTop_endcap2->SetLogx(1);
    c_PR_fromTop_endcap2->SetRightMargin(0.05);
    c_PR_fromTop_endcap2->SetTopMargin(0.05);
    c_PR_fromTop_endcap2->SetBottomMargin(0.12);
    c_PR_fromTop_endcap2->SetLeftMargin(0.13);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_endcap2_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_fromTop_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_fromTop_endcap2_nume[1]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_fromTop_endcap2_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_fromTop_endcap2_nume[1]->Draw();
    h_PR_pT_fromTop_endcap2_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_fromTop_endcap2->Update();

    TCanvas *c_PR_notFromTop_barrel  = new TCanvas("c_PR_notFromTop_barrel",  "PR (not from Top) barrel", 800, 800);
    c_PR_notFromTop_barrel->cd();
    c_PR_notFromTop_barrel->SetGrid(1);
    c_PR_notFromTop_barrel->SetLogx(1);
    c_PR_notFromTop_barrel->SetRightMargin(0.05);
    c_PR_notFromTop_barrel->SetTopMargin(0.05);
    c_PR_notFromTop_barrel->SetBottomMargin(0.12);
    c_PR_notFromTop_barrel->SetLeftMargin(0.13);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_notFromTop_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_notFromTop_barrel_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_notFromTop_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_notFromTop_barrel_nume[0]->Draw();
    h_PR_pT_notFromTop_barrel_nume[2]->Draw("same");
    h_PR_pT_notFromTop_barrel_nume[3]->Draw("same");
    h_PR_pT_notFromTop_barrel_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_notFromTop_barrel->Update();

    TCanvas *c_PR_notFromTop_endcap  = new TCanvas("c_PR_notFromTop_endcap",  "PR (not from top) endcap", 800, 800);
    c_PR_notFromTop_endcap->cd();
    c_PR_notFromTop_endcap->SetGrid(1);
    c_PR_notFromTop_endcap->SetLogx(1);
    c_PR_notFromTop_endcap->SetRightMargin(0.05);
    c_PR_notFromTop_endcap->SetTopMargin(0.05);
    c_PR_notFromTop_endcap->SetBottomMargin(0.12);
    c_PR_notFromTop_endcap->SetLeftMargin(0.13);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_notFromTop_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_notFromTop_endcap_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_notFromTop_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_notFromTop_endcap_nume[0]->Draw();
    h_PR_pT_notFromTop_endcap_nume[2]->Draw("same");
    h_PR_pT_notFromTop_endcap_nume[3]->Draw("same");
    h_PR_pT_notFromTop_endcap_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_notFromTop_endcap->Update();

    TCanvas *c_PR_notFromTop_endcap2 = new TCanvas("c_PR_notFromTop_endcap2", "PR (not from top) endcap2", 800, 800);
    c_PR_notFromTop_endcap2->cd();
    c_PR_notFromTop_endcap2->SetGrid(1);
    c_PR_notFromTop_endcap2->SetLogx(1);
    c_PR_notFromTop_endcap2->SetRightMargin(0.05);
    c_PR_notFromTop_endcap2->SetTopMargin(0.05);
    c_PR_notFromTop_endcap2->SetBottomMargin(0.12);
    c_PR_notFromTop_endcap2->SetLeftMargin(0.13);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitle("p_{T} (#mu) [GeV/c]");
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_pT_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetXaxis()->SetRangeUser(52, 1000);
    h_PR_pT_notFromTop_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_pT_notFromTop_endcap2_nume[0]->Draw();
    h_PR_pT_notFromTop_endcap2_nume[2]->Draw("same");
    h_PR_pT_notFromTop_endcap2_nume[3]->Draw("same");
    h_PR_pT_notFromTop_endcap2_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_notFromTop_endcap2->Update();

    TCanvas *c_PR_eta_fromTop = new TCanvas("c_PR_eta_fromTop", "PR eta (from top)", 800, 800);
    c_PR_eta_fromTop->cd();
    c_PR_eta_fromTop->SetGrid(1);
    c_PR_eta_fromTop->SetRightMargin(0.05);
    c_PR_eta_fromTop->SetTopMargin(0.05);
    c_PR_eta_fromTop->SetBottomMargin(0.12);
    c_PR_eta_fromTop->SetLeftMargin(0.13);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetTitle("#eta (#mu)");
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_eta_fromTop_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_eta_fromTop_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_eta_fromTop_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_eta_fromTop_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_eta_fromTop_nume[1]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_PR_eta_fromTop_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_eta_fromTop_nume[1]->Draw();
    h_PR_eta_fromTop_nume[2]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_eta_fromTop->Update();

    TCanvas *c_PR_eta_notFromTop = new TCanvas("c_PR_eta_notFromTop", "PR eta (not from top)", 800, 800);
    c_PR_eta_notFromTop->cd();
    c_PR_eta_notFromTop->SetGrid(1);
    c_PR_eta_notFromTop->SetRightMargin(0.05);
    c_PR_eta_notFromTop->SetTopMargin(0.05);
    c_PR_eta_notFromTop->SetBottomMargin(0.12);
    c_PR_eta_notFromTop->SetLeftMargin(0.13);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetTitle("#eta (#mu)");
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_eta_notFromTop_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_eta_notFromTop_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_eta_notFromTop_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_eta_notFromTop_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_eta_notFromTop_nume[0]->GetXaxis()->SetRangeUser(-2.4, 2.4);
    h_PR_eta_notFromTop_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_eta_notFromTop_nume[0]->Draw();
    h_PR_eta_notFromTop_nume[2]->Draw("same");
    h_PR_eta_notFromTop_nume[3]->Draw("same");
    h_PR_eta_notFromTop_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_eta_fromTop->Update();

    TCanvas *c_PR_nJet_fromTop_barrel  = new TCanvas("c_PR_nJet_fromTop_barrel",  "PR(nJet) barrel (from top)", 800, 800);
    c_PR_nJet_fromTop_barrel->cd();
    c_PR_nJet_fromTop_barrel->SetGrid(1);
    c_PR_nJet_fromTop_barrel->SetRightMargin(0.05);
    c_PR_nJet_fromTop_barrel->SetTopMargin(0.05);
    c_PR_nJet_fromTop_barrel->SetBottomMargin(0.12);
    c_PR_nJet_fromTop_barrel->SetLeftMargin(0.13);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_barrel_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_fromTop_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_fromTop_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_fromTop_barrel_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_fromTop_barrel_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_fromTop_barrel_nume[1]->Draw();
    h_PR_nJet_fromTop_barrel_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_nJet_fromTop_barrel->Update();

    TCanvas *c_PR_nJet_fromTop_endcap  = new TCanvas("c_PR_nJet_fromTop_endcap",  "PR(nJet) endcap (from top)", 800, 800);
    c_PR_nJet_fromTop_endcap->cd();
    c_PR_nJet_fromTop_endcap->SetGrid(1);
    c_PR_nJet_fromTop_endcap->SetRightMargin(0.05);
    c_PR_nJet_fromTop_endcap->SetTopMargin(0.05);
    c_PR_nJet_fromTop_endcap->SetBottomMargin(0.12);
    c_PR_nJet_fromTop_endcap->SetLeftMargin(0.13);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_endcap_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_fromTop_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_fromTop_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_fromTop_endcap_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_fromTop_endcap_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_fromTop_endcap_nume[1]->Draw();
    h_PR_nJet_fromTop_endcap_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_nJet_fromTop_endcap->Update();

    TCanvas *c_PR_nJet_fromTop_endcap2 = new TCanvas("c_PR_nJet_fromTop_endcap2", "PR(nJet) endcap2 (from top)", 800, 800);
    c_PR_nJet_fromTop_endcap2->cd();
    c_PR_nJet_fromTop_endcap2->SetGrid(1);
    c_PR_nJet_fromTop_endcap2->SetRightMargin(0.05);
    c_PR_nJet_fromTop_endcap2->SetTopMargin(0.05);
    c_PR_nJet_fromTop_endcap2->SetBottomMargin(0.12);
    c_PR_nJet_fromTop_endcap2->SetLeftMargin(0.13);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_fromTop_endcap2_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_fromTop_endcap2_nume[1]->Draw();
    h_PR_nJet_fromTop_endcap2_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_nJet_fromTop_endcap2->Update();

    TCanvas *c_PR_nJet_notFromTop_barrel  = new TCanvas("c_PR_nJet_notFromTop_barrel",  "PR(nJet) barrel (not from top)", 800, 800);
    c_PR_nJet_notFromTop_barrel->cd();
    c_PR_nJet_notFromTop_barrel->SetGrid(1);
    c_PR_nJet_notFromTop_barrel->SetRightMargin(0.05);
    c_PR_nJet_notFromTop_barrel->SetTopMargin(0.05);
    c_PR_nJet_notFromTop_barrel->SetBottomMargin(0.12);
    c_PR_nJet_notFromTop_barrel->SetLeftMargin(0.13);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_notFromTop_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_notFromTop_barrel_nume[0]->Draw();
    h_PR_nJet_notFromTop_barrel_nume[2]->Draw("same");
    h_PR_nJet_notFromTop_barrel_nume[3]->Draw("same");
    h_PR_nJet_notFromTop_barrel_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_nJet_notFromTop_barrel->Update();

    TCanvas *c_PR_nJet_notFromTop_endcap  = new TCanvas("c_PR_nJet_notFromTop_endcap",  "PR(nJet) endcap (not from top)", 800, 800);
    c_PR_nJet_notFromTop_endcap->cd();
    c_PR_nJet_notFromTop_endcap->SetGrid(1);
    c_PR_nJet_notFromTop_endcap->SetRightMargin(0.05);
    c_PR_nJet_notFromTop_endcap->SetTopMargin(0.05);
    c_PR_nJet_notFromTop_endcap->SetBottomMargin(0.12);
    c_PR_nJet_notFromTop_endcap->SetLeftMargin(0.13);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_notFromTop_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_notFromTop_endcap_nume[0]->Draw();
    h_PR_nJet_notFromTop_endcap_nume[2]->Draw("same");
    h_PR_nJet_notFromTop_endcap_nume[3]->Draw("same");
    h_PR_nJet_notFromTop_endcap_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_nJet_notFromTop_endcap->Update();

    TCanvas *c_PR_nJet_notFromTop_endcap2 = new TCanvas("c_PR_nJet_notFromTop_endcap2", "PR(nJet) endcap2 (from from top)", 800, 800);
    c_PR_nJet_notFromTop_endcap2->cd();
    c_PR_nJet_notFromTop_endcap2->SetGrid(1);
    c_PR_nJet_notFromTop_endcap2->SetRightMargin(0.05);
    c_PR_nJet_notFromTop_endcap2->SetTopMargin(0.05);
    c_PR_nJet_notFromTop_endcap2->SetBottomMargin(0.12);
    c_PR_nJet_notFromTop_endcap2->SetLeftMargin(0.13);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitle("N_{Jets}");
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetXaxis()->SetRangeUser(0-0.5, 20-0.5);
    h_PR_nJet_notFromTop_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_nJet_notFromTop_endcap2_nume[0]->Draw();
    h_PR_nJet_notFromTop_endcap2_nume[2]->Draw("same");
    h_PR_nJet_notFromTop_endcap2_nume[3]->Draw("same");
    h_PR_nJet_notFromTop_endcap2_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_nJet_notFromTop_endcap2->Update();

    TCanvas *c_PR_MET_fromTop_barrel  = new TCanvas("c_PR_MET_fromTop_barrel",  "PR(MET) barrel (from top)", 800, 800);
    c_PR_MET_fromTop_barrel->cd();
    c_PR_MET_fromTop_barrel->SetGrid(1);
    c_PR_MET_fromTop_barrel->SetRightMargin(0.05);
    c_PR_MET_fromTop_barrel->SetTopMargin(0.05);
    c_PR_MET_fromTop_barrel->SetBottomMargin(0.12);
    c_PR_MET_fromTop_barrel->SetLeftMargin(0.13);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_barrel_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_fromTop_barrel_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_barrel_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_fromTop_barrel_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_fromTop_barrel_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_fromTop_barrel_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_fromTop_barrel_nume[1]->Draw();
    h_PR_MET_fromTop_barrel_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_MET_fromTop_barrel->Update();

    TCanvas *c_PR_MET_fromTop_endcap  = new TCanvas("c_PR_MET_fromTop_endcap",  "PR(MET) endcap (from top)", 800, 800);
    c_PR_MET_fromTop_endcap->cd();
    c_PR_MET_fromTop_endcap->SetGrid(1);
    c_PR_MET_fromTop_endcap->SetRightMargin(0.05);
    c_PR_MET_fromTop_endcap->SetTopMargin(0.05);
    c_PR_MET_fromTop_endcap->SetBottomMargin(0.12);
    c_PR_MET_fromTop_endcap->SetLeftMargin(0.13);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_endcap_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_fromTop_endcap_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_endcap_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_fromTop_endcap_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_fromTop_endcap_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_fromTop_endcap_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_fromTop_endcap_nume[1]->Draw();
    h_PR_MET_fromTop_endcap_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_MET_fromTop_endcap->Update();

    TCanvas *c_PR_MET_fromTop_endcap2 = new TCanvas("c_PR_MET_fromTop_endcap2", "PR(MET) endcap2 (from top)", 800, 800);
    c_PR_MET_fromTop_endcap2->cd();
    c_PR_MET_fromTop_endcap2->SetGrid(1);
    c_PR_MET_fromTop_endcap2->SetRightMargin(0.05);
    c_PR_MET_fromTop_endcap2->SetTopMargin(0.05);
    c_PR_MET_fromTop_endcap2->SetBottomMargin(0.12);
    c_PR_MET_fromTop_endcap2->SetLeftMargin(0.13);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_endcap2_nume[1]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_fromTop_endcap2_nume[1]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_fromTop_endcap2_nume[1]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_fromTop_endcap2_nume[1]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_fromTop_endcap2_nume[1]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_fromTop_endcap2_nume[1]->Draw();
    h_PR_MET_fromTop_endcap2_nume[2]->Draw("same");
    legend_fromTop->Draw();
    c_PR_MET_fromTop_endcap2->Update();

    TCanvas *c_PR_MET_notFromTop_barrel  = new TCanvas("c_PR_MET_notFromTop_barrel",  "PR(MET) barrel (not from top)", 800, 800);
    c_PR_MET_notFromTop_barrel->cd();
    c_PR_MET_notFromTop_barrel->SetGrid(1);
    c_PR_MET_notFromTop_barrel->SetRightMargin(0.05);
    c_PR_MET_notFromTop_barrel->SetTopMargin(0.05);
    c_PR_MET_notFromTop_barrel->SetBottomMargin(0.12);
    c_PR_MET_notFromTop_barrel->SetLeftMargin(0.13);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_barrel_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_barrel_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_notFromTop_barrel_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_notFromTop_barrel_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_notFromTop_barrel_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_notFromTop_barrel_nume[0]->Draw();
    h_PR_MET_notFromTop_barrel_nume[2]->Draw("same");
    h_PR_MET_notFromTop_barrel_nume[3]->Draw("same");
    h_PR_MET_notFromTop_barrel_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_MET_notFromTop_barrel->Update();

    TCanvas *c_PR_MET_notFromTop_endcap  = new TCanvas("c_PR_MET_notFromTop_endcap",  "PR(MET) endcap (not from top)", 800, 800);
    c_PR_MET_notFromTop_endcap->cd();
    c_PR_MET_notFromTop_endcap->SetGrid(1);
    c_PR_MET_notFromTop_endcap->SetRightMargin(0.05);
    c_PR_MET_notFromTop_endcap->SetTopMargin(0.05);
    c_PR_MET_notFromTop_endcap->SetBottomMargin(0.12);
    c_PR_MET_notFromTop_endcap->SetLeftMargin(0.13);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_endcap_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_endcap_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_notFromTop_endcap_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_notFromTop_endcap_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_notFromTop_endcap_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_notFromTop_endcap_nume[0]->Draw();
    h_PR_MET_notFromTop_endcap_nume[2]->Draw("same");
    h_PR_MET_notFromTop_endcap_nume[3]->Draw("same");
    h_PR_MET_notFromTop_endcap_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_MET_notFromTop_endcap->Update();

    TCanvas *c_PR_MET_notFromTop_endcap2 = new TCanvas("c_PR_MET_notFromTop_endcap2", "PR(MET) endcap2 (from from top)", 800, 800);
    c_PR_MET_notFromTop_endcap2->cd();
    c_PR_MET_notFromTop_endcap2->SetGrid(1);
    c_PR_MET_notFromTop_endcap2->SetRightMargin(0.05);
    c_PR_MET_notFromTop_endcap2->SetTopMargin(0.05);
    c_PR_MET_notFromTop_endcap2->SetBottomMargin(0.12);
    c_PR_MET_notFromTop_endcap2->SetLeftMargin(0.13);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitle("MET [GeV]");
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleOffset(1);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitle("Prompt rate");
    h_PR_MET_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleSize(0.05);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetYaxis()->SetTitleOffset(1.25);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetYaxis()->SetLabelSize(0.04);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetNoExponent(1);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetMoreLogLabels(1);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetXaxis()->SetRangeUser(0, 80);
    h_PR_MET_notFromTop_endcap2_nume[0]->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_PR_MET_notFromTop_endcap2_nume[0]->Draw();
    h_PR_MET_notFromTop_endcap2_nume[2]->Draw("same");
    h_PR_MET_notFromTop_endcap2_nume[3]->Draw("same");
    h_PR_MET_notFromTop_endcap2_nume[4]->Draw("same");
    legend_notFromTop->Draw();
    c_PR_MET_notFromTop_endcap2->Update();

//    DYAnalyzer a("");
//    a.SaveCanvases("ALL", "~/Desktop/LastCopy_DYQCD", "png");

    // Writing
    TString outname = "/media/sf_DATA/FR/Muon/MC_FakeRates.root";
    cout << "Writing to " << outname << " ... ";
    TFile *f_out = new TFile(outname, "RECREATE");
    for (Int_t i=0; i<6; i++)
    {
        f_out->WriteObjectAny(h_FR_pT_barrel_nume[i],  "TH1D", "h_FR_barrel_" +type[i]);
        f_out->WriteObjectAny(h_FR_pT_endcap_nume[i],  "TH1D", "h_FR_endcap_" +type[i]);
        f_out->WriteObjectAny(h_FR_pT_endcap2_nume[i], "TH1D", "h_FR_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_FR_eta_nume[i],        "TH1D", "h_FR_eta_"    +type[i]);
        f_out->WriteObjectAny(h_FR_nJet_barrel_nume[i],  "TH1D", "h_FR_nJet_barrel_" +type[i]);
        f_out->WriteObjectAny(h_FR_nJet_endcap_nume[i],  "TH1D", "h_FR_nJet_endcap_" +type[i]);
        f_out->WriteObjectAny(h_FR_nJet_endcap2_nume[i], "TH1D", "h_FR_nJet_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_pT_barrel_nume[i],  "TH1D", "h_PR_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_endcap_nume[i],  "TH1D", "h_PR_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_endcap2_nume[i], "TH1D", "h_PR_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_eta_nume[i],        "TH1D", "h_PR_eta_"    +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_barrel_nume[i],  "TH1D", "h_PR_nJet_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_endcap_nume[i],  "TH1D", "h_PR_nJet_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_endcap2_nume[i], "TH1D", "h_PR_nJet_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_pT_fromTop_barrel_nume[i],  "TH1D", "h_PR_fromTop_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_fromTop_endcap_nume[i],  "TH1D", "h_PR_fromTop_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_fromTop_endcap2_nume[i], "TH1D", "h_PR_fromTop_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_pT_notFromTop_barrel_nume[i],  "TH1D", "h_PR_notFromTop_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_notFromTop_endcap_nume[i],  "TH1D", "h_PR_notFromTop_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_pT_notFromTop_endcap2_nume[i], "TH1D", "h_PR_notFromTop_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_eta_fromTop_nume[i],        "TH1D", "h_PR_fromTop_eta_"    +type[i]);
        f_out->WriteObjectAny(h_PR_eta_notFromTop_nume[i],        "TH1D", "h_PR_notFromTop_eta_"    +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_fromTop_barrel_nume[i],  "TH1D", "h_PR_nJet_fromTop_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_fromTop_endcap_nume[i],  "TH1D", "h_PR_nJet_fromTop_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_fromTop_endcap2_nume[i], "TH1D", "h_PR_nJet_fromTop_endcap2_"+type[i]);
        f_out->WriteObjectAny(h_PR_nJet_notFromTop_barrel_nume[i],  "TH1D", "h_PR_nJet_notFromTop_barrel_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_notFromTop_endcap_nume[i],  "TH1D", "h_PR_nJet_notFromTop_endcap_" +type[i]);
        f_out->WriteObjectAny(h_PR_nJet_notFromTop_endcap2_nume[i], "TH1D", "h_PR_nJet_notFromTop_endcap2_"+type[i]);
    }

    f_out->Close();
    cout << "Finished." << endl;
    if (!f->IsOpen()) cout << "File has been closed successfully.\n" << endl;
    else cout << "FILE COULD NOT BE CLOSED!\n" << endl;

    cout << "Average rates:" << endl;
    for (Int_t i=0; i<6; i++)
    {
        Double_t FRavg = h_FR_nume_sum[i]->Integral() / h_FR_deno_sum[i]->Integral();
        Double_t PRavg = h_PR_nume_sum[i]->Integral() / h_PR_deno_sum[i]->Integral();
        cout << type[i] << ": FR = " << FRavg << "   PR = " << PRavg << endl;
    }
    cout << endl;

} // End of Mu_EstFRandPR_MC()


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
            h->SetBinError(i, 1);
        }
    }
}
