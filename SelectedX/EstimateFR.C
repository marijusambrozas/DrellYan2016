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
void Mu_EstFR (Int_t type);
void E_DYefficiency();

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
        if (whichX.Contains("DY"))
        {
            cout << "\n*******      E_DYefficiency()      *******" << endl;
            E_DYefficiency();
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

    TH1D *h_pT_barrel_MC_nume[_EndOf_Data_Special], *h_pT_endcap_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_data_nume,*h_pT_endcap_data_nume, *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl,
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
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr1]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr1]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_eta_MC_ctrl[pr1]);
        removeNegativeBins(h_eta_MC_nume[pr1]);

        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);
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
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_DY")));
            h_eta_MC_ctrl[_DY_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_DY")));
            h_eta_MC_nume[_DY_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_MC_nume[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
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
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _GJets_20to100)
        {
            h_pT_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_GJets")));
            h_pT_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_GJets")));
            h_pT_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_GJets")));
            h_pT_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_GJets")));
            h_eta_MC_ctrl[_GJets_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_GJets")));
            h_eta_MC_nume[_GJets_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_GJets")));

            h_pT_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_eta_MC_nume[_GJets_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_GJets_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_GJets_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_GJets_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_GJets_Full]->Add(h_pT_endcap_MC_nume[pr]);
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
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_MC_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_MC_nume[pr]);

        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_eta_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC_nume[pr]);

        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_eta_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC_nume[pr]->SetDirectory(0);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_MC_ctrl_QCD")));
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_MC_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_MC_nume_QCD")));
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_MC_nume_QCD")));
            h_eta_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_MC_ctrl[pr]->Clone("h_eta_MC_ctrl_QCD")));
            h_eta_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_MC_nume[pr]->Clone("h_eta_MC_nume_QCD")));

            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_eta_MC_ctrl[_QCDEMEnriched_Full]->Add(h_eta_MC_ctrl[pr]);
            h_eta_MC_nume[_QCDEMEnriched_Full]->Add(h_eta_MC_nume[pr]);
        }
        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        TH1D *h_temp[6];
        if (pr == _SinglePhoton_B)
        {
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            file->GetObject("h_eta_ctrl", h_eta_data_ctrl);
            file->GetObject("h_eta_nume", h_eta_data_nume);

            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_nume);
            removeNegativeBins(h_eta_data_ctrl);
            removeNegativeBins(h_eta_data_nume);
        }
        else
        {
            file->GetObject("h_pT_barrel_ctrl", h_temp[0]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[1]);
            file->GetObject("h_pT_barrel_nume", h_temp[2]);
            file->GetObject("h_pT_endcap_nume", h_temp[3]);
            file->GetObject("h_eta_ctrl", h_temp[4]);
            file->GetObject("h_eta_nume", h_temp[5]);

            removeNegativeBins(h_temp[0]);
            removeNegativeBins(h_temp[1]);
            removeNegativeBins(h_temp[2]);
            removeNegativeBins(h_temp[3]);
            removeNegativeBins(h_temp[4]);
            removeNegativeBins(h_temp[5]);

            h_pT_barrel_data_ctrl->Add(h_temp[0]);
            h_pT_endcap_data_ctrl->Add(h_temp[1]);
            h_pT_barrel_data_nume->Add(h_temp[2]);
            h_pT_endcap_data_nume->Add(h_temp[3]);
            h_eta_data_ctrl->Add(h_temp[4]);
            h_eta_data_nume->Add(h_temp[5]);
        }
    }

    h_pT_barrel_data_ctrl->SetDirectory(0);
    h_pT_endcap_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume->SetDirectory(0);
    h_pT_endcap_data_nume->SetDirectory(0);
    h_eta_data_ctrl->SetDirectory(0);
    h_eta_data_nume->SetDirectory(0);


//--------------------------------- FR from QCD MC -------------------------------------- (deno = nume + ctrl)

    //             QCD_MC_nume
    // FR = -------------------------
    //      QCD_MC_nume + QCD_MC_ctrl

    // ------ Numerator ------ //
    TH1D *h_pT_barrel_nume_fMC = ((TH1D*)(h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Clone("h_pT_barrel_nume_fMC")));
    TH1D *h_pT_endcap_nume_fMC = ((TH1D*)(h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Clone("h_pT_endcap_nume_fMC")));
    TH1D *h_eta_nume_fMC = ((TH1D*)(h_eta_MC_nume[_QCDEMEnriched_Full]->Clone("h_eta_nume_fMC")));

    // ------ Denominator ------ //
    TH1D *h_pT_barrel_deno_fMC = ((TH1D*)(h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_barrel_deno_fMC")));
    h_pT_barrel_deno_fMC->Add(h_pT_barrel_nume_fMC); // deno = sig+ctrl
    TH1D *h_pT_endcap_deno_fMC = ((TH1D*)(h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Clone("h_pT_endcap_deno_fMC")));
    h_pT_endcap_deno_fMC->Add(h_pT_endcap_nume_fMC); // deno = sig+ctrl
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

    for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
    {
        Double_t sig_barrel_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_barrel_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_barrel_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_barrel_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);

        Double_t sig_endcap_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t ctrl_endcap_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinContent(i_bin);
        Double_t err_sig_endcap_MC = h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->GetBinError(i_bin);
        Double_t err_ctrl_endcap_MC = h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->GetBinError(i_bin);

        Double_t sig_barrel_ratio = h_pT_barrel_nume->GetBinContent(i_bin);
        Double_t ctrl_barrel_ratio = h_pT_barrel_ctrl->GetBinContent(i_bin);
        Double_t err_sig_barrel_ratio = h_pT_barrel_nume->GetBinError(i_bin);
        Double_t err_ctrl_barrel_ratio = h_pT_barrel_ctrl->GetBinError(i_bin);

        Double_t sig_endcap_ratio = h_pT_endcap_nume->GetBinContent(i_bin);
        Double_t ctrl_endcap_ratio = h_pT_endcap_ctrl->GetBinContent(i_bin);
        Double_t err_sig_endcap_ratio = h_pT_endcap_nume->GetBinError(i_bin);
        Double_t err_ctrl_endcap_ratio = h_pT_endcap_ctrl->GetBinError(i_bin);

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

// ---------------------------- Writing------------------------------ //
    TString filename = "/media/sf_DATA/FR/Electron/FakeRate_electron.root";
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
//    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(6.2192e+05/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // MT<60GeV
//    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.0304e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.4895e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(2.5956e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(3.2371e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(4.8312e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.0278e+06/h_pT_barrel_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(1.5538e+07/h_pT_barrel_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.7131e+05/h_pT_barrel_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(4.4023e+06/h_pT_barrel_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Scale(5.9007e+04/h_pT_barrel_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Scale(9.3672e+05/h_pT_barrel_template_deno_100to500[_QCDMuEnriched_Full]->Integral());


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
//    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(7.5101e+05/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral()); // MT<60GeV
//    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(5.3721e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(1.7019e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.2104e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(2.6946e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(1.8185e+05 /h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Scale(1.2750e+06/h_pT_endcap_template_nume_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Scale(7.7560e+06/h_pT_endcap_template_deno_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Scale(2.8650e+05/h_pT_endcap_template_nume_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Scale(1.9793e+06/h_pT_endcap_template_deno_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Scale(5.9687e+04/h_pT_endcap_template_nume_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Scale(3.4903e+05/h_pT_endcap_template_deno_100to500[_QCDMuEnriched_Full]->Integral());

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
//    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(9.7350e+06/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral()); // MT<60GeV
//    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(2.4554e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(4.5693e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(4.6752e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.0593e+06/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
//    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(1.5689e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(1.4758e+07/h_pT_barrel_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(4.1938e+06/h_pT_barrel_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(8.8585e+05/h_pT_barrel_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Scale(6.8002e+06/h_pT_endcap_template_ctrl_50to70  [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Scale(1.7355e+06/h_pT_endcap_template_ctrl_70to100 [_QCDMuEnriched_Full]->Integral());
    h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Scale(3.0191e+05/h_pT_endcap_template_ctrl_100to500[_QCDMuEnriched_Full]->Integral());

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
    TH1D *h_FRsigCtrl_template_endcap1 = ((TH1D*)(h_FRsigCtrl_template_endcap->Clone("h_FRsigCtrl_template_endcap1")));
    TH1D *h_FRratio_endcap1 = ((TH1D*)(h_FRratio_endcap->Clone("h_FRratio_endcap1")));
    TH1D *h_FRratio_barrel1 = ((TH1D*)(h_FRratio_barrel->Clone("h_FRratio_barrel1")));
    h_FRsigCtrl_template_endcap1->SetTitle("");
    h_FRsigCtrl_template_endcap1->SetMarkerStyle(kFullSquare);
    h_FRsigCtrl_template_endcap1->SetMarkerColor(kBlue);
    h_FRratio_endcap1->GetYaxis()->SetTitle("Signalo srities tikimyb#dot{e}");
    h_FRratio_endcap1->GetYaxis()->SetTitleSize(0.045);
    h_FRratio_endcap1->GetYaxis()->SetTitleSize(0.045);
    h_FRratio_endcap1->GetYaxis()->SetRangeUser(0, 1);
    h_FRratio_endcap1->SetMarkerColor(kOrange+8);
    h_FRratio_barrel1->SetMarkerStyle(kFullDotLarge);
    h_FRratio_endcap1->Draw();
    h_FRsigCtrl_template_endcap1->Draw("same");
    h_FRsigCtrl_template_barrel->Draw("same");
    h_FRratio_barrel1->Draw("same");
    TLegend *legend1 = new TLegend(0.13, 0.77, 0.6, 0.95);
    legend1->AddEntry(h_FRratio_barrel1, "Santykio metodas (|#eta| < 1.2)", "LP");
    legend1->AddEntry(h_FRratio_endcap1, "Santykio metodas (|#eta| #geq 1.2)", "LP");
    legend1->AddEntry(h_FRsigCtrl_template_barrel, "Sablonu priderinimas (|#eta| < 1.2)", "LP");
    legend1->AddEntry(h_FRsigCtrl_template_endcap1, "Sablonu priderinimas (|#eta| #geq 1.2)", "LP");
    legend1->Draw();
    c_FR_allin1->SetLogx();
    c_FR_allin1->Update();

} // End of Mu_EstFR()


void E_DYefficiency()
{
    TFile *f = new TFile("/media/sf_DATA/SelectedEE/Histos/DYefficiency.root", "READ");

    TH1D *h_pT_barrel_pass[3],
         *h_pT_endcap_pass[3],
         *h_pT_barrel_fail[3],
         *h_pT_endcap_fail[3];
    TString type[3] = {"data", "DY", "bkg"};

    for (Int_t i=2; i>=0; i--)
    {
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
    }

    TH1D *h_eff_data_barrel = ((TH1D*)(h_pT_barrel_pass[0]->Clone("h_eff_data_barrel")));
    TH1D *h_eff_data_endcap = ((TH1D*)(h_pT_endcap_pass[0]->Clone("h_eff_data_endcap")));
    TH1D *h_eff_MC_barrel = ((TH1D*)(h_pT_barrel_pass[1]->Clone("h_eff_MC_barrel")));
    TH1D *h_eff_MC_endcap = ((TH1D*)(h_pT_endcap_pass[1]->Clone("h_eff_MC_endcap")));
    TH1D *h_eff_bkg_barrel = ((TH1D*)(h_pT_barrel_pass[2]->Clone("h_eff_bkg_barrel")));
    TH1D *h_eff_bkg_endcap = ((TH1D*)(h_pT_endcap_pass[2]->Clone("h_eff_bkg_endcap")));
    h_eff_data_barrel->Add(h_eff_bkg_barrel, -1);
    h_eff_data_endcap->Add(h_eff_bkg_endcap, -1);

    TH1D *h_ineff_data_barrel = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_ineff_data_barrel")));
    TH1D *h_ineff_data_endcap = ((TH1D*)(h_pT_endcap_fail[0]->Clone("h_ineff_data_endcap")));
    TH1D *h_ineff_MC_barrel = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_ineff_MC_barrel")));
    TH1D *h_ineff_MC_endcap = ((TH1D*)(h_pT_endcap_fail[1]->Clone("h_ineff_MC_endcap")));
    TH1D *h_ineff_bkg_barrel = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_ineff_bkg_barrel")));
    TH1D *h_ineff_bkg_endcap = ((TH1D*)(h_pT_endcap_fail[2]->Clone("h_ineff_bkg_endcap")));
    h_ineff_data_barrel->Add(h_ineff_bkg_barrel, -1);
    h_ineff_data_endcap->Add(h_ineff_bkg_endcap, -1);

    TH1D *h_eff_data_barrel_deno = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_eff_data_barrel_deno")));
    TH1D *h_eff_data_endcap_deno = ((TH1D*)(h_pT_barrel_fail[0]->Clone("h_eff_data_endcap_deno")));
    TH1D *h_eff_MC_barrel_deno = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_eff_MC_barrel_deno")));
    TH1D *h_eff_MC_endcap_deno = ((TH1D*)(h_pT_barrel_fail[1]->Clone("h_eff_MC_endcap_deno")));
    TH1D *h_eff_bkg_barrel_deno = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_eff_bkg_barrel_deno")));
    TH1D *h_eff_bkg_endcap_deno = ((TH1D*)(h_pT_barrel_fail[2]->Clone("h_eff_bkg_endcap_deno")));
    h_eff_data_barrel_deno->Add(h_pT_barrel_pass[0]);
    h_eff_data_endcap_deno->Add(h_pT_endcap_pass[0]);
    h_eff_MC_barrel_deno->Add(h_pT_barrel_pass[1]);
    h_eff_MC_endcap_deno->Add(h_pT_endcap_pass[1]);
    h_eff_bkg_barrel_deno->Add(h_pT_barrel_pass[2]);
    h_eff_bkg_endcap_deno->Add(h_pT_endcap_pass[2]);
    h_eff_data_barrel_deno->Add(h_eff_bkg_barrel_deno, -1);
    h_eff_data_endcap_deno->Add(h_eff_bkg_endcap_deno, -1);

    h_eff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_eff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_eff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_eff_MC_endcap->Divide(h_eff_MC_endcap_deno);
    h_ineff_data_barrel->Divide(h_eff_data_barrel_deno);
    h_ineff_data_endcap->Divide(h_eff_data_endcap_deno);
    h_ineff_MC_barrel->Divide(h_eff_MC_barrel_deno);
    h_ineff_MC_endcap->Divide(h_eff_MC_endcap_deno);

    h_eff_data_barrel->SetDirectory(0);
    h_eff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_eff_data_barrel->SetMarkerColor(kBlack);
    h_eff_data_barrel->SetLineColor(kBlack);
    h_eff_data_endcap->SetDirectory(0);
    h_eff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_eff_data_endcap->SetMarkerColor(kBlack);
    h_eff_data_endcap->SetLineColor(kBlack);
    h_eff_MC_barrel->SetDirectory(0);
    h_eff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_eff_MC_barrel->SetMarkerColor(kRed);
    h_eff_MC_barrel->SetLineColor(kRed);
    h_eff_MC_endcap->SetDirectory(0);
    h_eff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_eff_MC_endcap->SetMarkerColor(kRed);
    h_eff_MC_endcap->SetLineColor(kRed);

    h_ineff_data_barrel->SetDirectory(0);
    h_ineff_data_barrel->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_barrel->SetMarkerColor(kBlack);
    h_ineff_data_barrel->SetLineColor(kBlack);
    h_ineff_data_endcap->SetDirectory(0);
    h_ineff_data_endcap->SetMarkerStyle(kFullDotLarge);
    h_ineff_data_endcap->SetMarkerColor(kBlack);
    h_ineff_data_endcap->SetLineColor(kBlack);
    h_ineff_MC_barrel->SetDirectory(0);
    h_ineff_MC_barrel->SetMarkerStyle(kFullSquare);
    h_ineff_MC_barrel->SetMarkerColor(kRed);
    h_ineff_MC_barrel->SetLineColor(kRed);
    h_ineff_MC_endcap->SetDirectory(0);
    h_ineff_MC_endcap->SetMarkerStyle(kFullSquare);
    h_ineff_MC_endcap->SetMarkerColor(kRed);
    h_ineff_MC_endcap->SetLineColor(kRed);

    TH1D *h_eff_ratio_barrel = ((TH1D*)(h_eff_data_barrel->Clone("h_eff_ratio_barrel")));
    TH1D *h_eff_ratio_endcap = ((TH1D*)(h_eff_data_endcap->Clone("h_eff_ratio_endcap")));
    TH1D *h_ineff_ratio_barrel = ((TH1D*)(h_ineff_data_barrel->Clone("h_ineff_ratio_barrel")));
    TH1D *h_ineff_ratio_endcap = ((TH1D*)(h_ineff_data_endcap->Clone("h_ineff_ratio_endcap")));
    h_eff_ratio_barrel->Divide(h_eff_MC_barrel);
    h_eff_ratio_endcap->Divide(h_eff_MC_endcap);
    h_eff_ratio_barrel->SetDirectory(0);
    h_eff_ratio_endcap->SetDirectory(0);
    h_eff_ratio_endcap->SetMarkerColor(kBlue);
    h_eff_ratio_endcap->SetLineColor(kBlue);
    h_ineff_ratio_barrel->Divide(h_ineff_MC_barrel);
    h_ineff_ratio_endcap->Divide(h_ineff_MC_endcap);
    h_ineff_ratio_barrel->SetDirectory(0);
    h_ineff_ratio_endcap->SetDirectory(0);
    h_ineff_ratio_endcap->SetMarkerColor(kBlue);
    h_ineff_ratio_endcap->SetLineColor(kBlue);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_eff_barrel = new myRatioPlot_t("c_eff_barrel", h_eff_MC_barrel, h_eff_data_barrel);
    myRatioPlot_t *RP_eff_endcap = new myRatioPlot_t("c_eff_endcap", h_eff_MC_endcap, h_eff_data_endcap);
    myRatioPlot_t *RP_ineff_barrel = new myRatioPlot_t("c_ineff_barrel", h_ineff_MC_barrel, h_ineff_data_barrel);
    myRatioPlot_t *RP_ineff_endcap = new myRatioPlot_t("c_ineff_endcap", h_ineff_MC_endcap, h_ineff_data_endcap);

    RP_eff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_eff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_barrel->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);
    RP_ineff_endcap->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 5000);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_eff_data_barrel, "Data", "pl");
    legend->AddEntry(h_eff_MC_barrel, "DY MC", "pl");

    RP_eff_barrel->ImportLegend(legend);
    RP_eff_endcap->ImportLegend(legend);
    RP_ineff_barrel->ImportLegend(legend);
    RP_ineff_endcap->ImportLegend(legend);    

    RP_eff_barrel->Draw(0.01, 1, 1, "", "Efficiency");
    RP_eff_endcap->Draw(0.01, 1, 1, "", "Efficiency");
    RP_ineff_barrel->Draw(0.01, 1, 1, "", "Inefficiency");
    RP_ineff_endcap->Draw(0.01, 1, 1, "", "Inefficiency");

    RP_eff_barrel->pad1->SetLogy(0);
    RP_eff_endcap->pad1->SetLogy(0);
    RP_ineff_barrel->pad1->SetLogy(0);
    RP_ineff_endcap->pad1->SetLogy(0);

    TLegend *legend_ratio = new TLegend(0.7, 0.8, 0.95, 0.95);
    legend_ratio->AddEntry(h_eff_ratio_barrel, "Barrel", "pl");
    legend_ratio->AddEntry(h_eff_ratio_endcap, "Endcap", "pl");

    TCanvas *c_eff_ratio = new TCanvas("c_eff_ratio", "c_eff_ratio", 800, 800);
    c_eff_ratio->SetTopMargin(0.05);
    c_eff_ratio->SetRightMargin(0.05);
    c_eff_ratio->SetBottomMargin(0.15);
    c_eff_ratio->SetLeftMargin(0.15);
    h_eff_ratio_barrel->SetTitle("");
    h_eff_ratio_barrel->GetXaxis()->SetTitle("p_{#lower[-0.2]{T}} [GeV/c]");
    h_eff_ratio_barrel->GetXaxis()->SetTitleSize(0.062);
    h_eff_ratio_barrel->GetXaxis()->SetTitleOffset(0.9);
    h_eff_ratio_barrel->GetXaxis()->SetLabelSize(0.048);
    h_eff_ratio_barrel->GetXaxis()->SetMoreLogLabels();
    h_eff_ratio_barrel->GetXaxis()->SetNoExponent();
    h_eff_ratio_barrel->GetYaxis()->SetTitle("Eff(Data)/Eff(MC)");
    h_eff_ratio_barrel->GetYaxis()->SetTitleSize(0.05);
    h_eff_ratio_barrel->GetYaxis()->SetTitleOffset(1.25);
    h_eff_ratio_barrel->GetYaxis()->SetLabelSize(0.043);
    h_eff_ratio_barrel->GetYaxis()->SetMoreLogLabels();
    h_eff_ratio_barrel->GetYaxis()->SetNoExponent();
    h_eff_ratio_barrel->Draw();
    h_eff_ratio_barrel->GetYaxis()->SetRangeUser(0.3, 2);
    h_eff_ratio_endcap->Draw("same");
    legend_ratio->Draw();
    c_eff_ratio->SetLogx();
    c_eff_ratio->SetGridx();
    c_eff_ratio->SetGridy();
    c_eff_ratio->Update();

    TCanvas *c_ineff_ratio = new TCanvas("c_ineff_ratio", "c_ineff_ratio", 800, 800);
    c_ineff_ratio->SetTopMargin(0.05);
    c_ineff_ratio->SetRightMargin(0.05);
    c_ineff_ratio->SetBottomMargin(0.15);
    c_ineff_ratio->SetLeftMargin(0.15);
    h_ineff_ratio_barrel->SetTitle("");
    h_ineff_ratio_barrel->GetXaxis()->SetTitle("p_{#lower[-0.2]{T}} [GeV/c]");
    h_ineff_ratio_barrel->GetXaxis()->SetTitleSize(0.062);
    h_ineff_ratio_barrel->GetXaxis()->SetTitleOffset(0.9);
    h_ineff_ratio_barrel->GetXaxis()->SetLabelSize(0.048);
    h_ineff_ratio_barrel->GetXaxis()->SetMoreLogLabels();
    h_ineff_ratio_barrel->GetXaxis()->SetNoExponent();
    h_ineff_ratio_barrel->GetYaxis()->SetTitle("Ineff(Data)/Ineff(MC)");
    h_ineff_ratio_barrel->GetYaxis()->SetTitleSize(0.05);
    h_ineff_ratio_barrel->GetYaxis()->SetTitleOffset(1.25);
    h_ineff_ratio_barrel->GetYaxis()->SetLabelSize(0.043);
    h_ineff_ratio_barrel->GetYaxis()->SetMoreLogLabels();
    h_ineff_ratio_barrel->GetYaxis()->SetNoExponent();
    h_ineff_ratio_barrel->Draw();
    h_ineff_ratio_barrel->GetYaxis()->SetRangeUser(0.5, 2);
    h_ineff_ratio_endcap->Draw("same");
    legend_ratio->Draw();
    c_ineff_ratio->SetLogx();
    c_ineff_ratio->SetGridx();
    c_ineff_ratio->SetGridy();
    c_ineff_ratio->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << "/media/sf_DATA/SelectedEE/Histos/DYefficiency.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " <<"/media/sf_DATA/SelectedEE/Histos/DYefficiency.root" << " COULD NOT BE CLOSED!\n" << endl;

    TFile *f_out = new TFile("/media/sf_DATA/SelectedEE/Histos/DY_effcorr.root", "RECREATE");
    f_out->cd();
    h_eff_data_barrel->Write();
    h_eff_data_endcap->Write();
    h_eff_MC_barrel->Write();
    h_eff_MC_endcap->Write();
    h_ineff_data_barrel->Write();
    h_ineff_data_endcap->Write();
    h_ineff_MC_barrel->Write();
    h_ineff_MC_endcap->Write();
    h_eff_ratio_barrel->Write();
    h_eff_ratio_endcap->Write();
    h_ineff_ratio_barrel->Write();
    h_ineff_ratio_endcap->Write();
    f_out->Close();

} // End of E_DYefficiency


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
