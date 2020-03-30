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
#include <TROOT.h>

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

void E_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr);
void E_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr);
void Mu_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr);
void Mu_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr);

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

void FR_HistDrawer (TString WhichX = "", Int_t type = 2, Int_t systErr = 0)
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
        if (whichX.Contains("QCD"))
        {
            cout << "\n*******     E_QCDest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            E_QCDest_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*******     E_WJETest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            E_WJETest_HistDrawer(type, systErr);
        }
        else
        {
            cout << "\n*******      E_HistDrawer(" << type << ")      *******" << endl;
            E_HistDrawer(type);
        }
    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*******     Mu_QCDest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            Mu_QCDest_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*******     Mu_WJETest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            Mu_WJETest_HistDrawer(type, systErr);
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
    FileMgr fm;
    THStack *s_PFiso_dBeta_barrel_nume = new THStack("s_PFiso_dBeta_barrel_nume", "");
    THStack *s_PFiso_dBeta_endcap_nume = new THStack("s_PFiso_dBeta_endcap_nume", "");
    THStack *s_PFiso_dBeta_barrel_deno = new THStack("s_PFiso_dBeta_barrel_deno", "");
    THStack *s_PFiso_dBeta_endcap_deno = new THStack("s_PFiso_dBeta_endcap_deno", "");
    THStack *s_PFiso_dBeta_barrel_ctrl = new THStack("s_PFiso_dBeta_barrel_ctrl", "");
    THStack *s_PFiso_dBeta_endcap_ctrl = new THStack("s_PFiso_dBeta_endcap_ctrl", "");
    THStack *s_PFiso_Rho_barrel_nume = new THStack("s_PFiso_Rho_barrel_nume", "");
    THStack *s_PFiso_Rho_endcap_nume = new THStack("s_PFiso_Rho_endcap_nume", "");
    THStack *s_PFiso_Rho_barrel_deno = new THStack("s_PFiso_Rho_barrel_deno", "");
    THStack *s_PFiso_Rho_endcap_deno = new THStack("s_PFiso_Rho_endcap_deno", "");
    THStack *s_PFiso_Rho_barrel_ctrl = new THStack("s_PFiso_Rho_barrel_ctrl", "");
    THStack *s_PFiso_Rho_endcap_ctrl = new THStack("s_PFiso_Rho_endcap_ctrl", "");
    THStack *s_pT_barrel_nume = new THStack("s_pT_barrel_nume", "");
    THStack *s_pT_endcap_nume = new THStack("s_pT_endcap_nume", "");
    THStack *s_pT_barrel_deno = new THStack("s_pT_barrel_deno", "");
    THStack *s_pT_endcap_deno = new THStack("s_pT_endcap_deno", "");
    THStack *s_pT_barrel_ctrl = new THStack("s_pT_barrel_ctrl", "");
    THStack *s_pT_endcap_ctrl = new THStack("s_pT_endcap_ctrl", "");
    THStack *s_SigmaIEtaIEta_barrel_nume = new THStack("s_SigmaIEtaIEta_barrel_nume", "");
    THStack *s_SigmaIEtaIEta_endcap_nume = new THStack("s_SigmaIEtaIEta_endcap_nume", "");
    THStack *s_SigmaIEtaIEta_barrel_deno = new THStack("s_SigmaIEtaIEta_barrel_deno", "");
    THStack *s_SigmaIEtaIEta_endcap_deno = new THStack("s_SigmaIEtaIEta_endcap_deno", "");
    THStack *s_SigmaIEtaIEta_barrel_ctrl = new THStack("s_SigmaIEtaIEta_barrel_ctrl", "");
    THStack *s_SigmaIEtaIEta_endcap_ctrl = new THStack("s_SigmaIEtaIEta_endcap_ctrl", "");
    THStack *s_dEtaInSeed_barrel_nume = new THStack("s_dEtaInSeed_barrel_nume", "");
    THStack *s_dEtaInSeed_endcap_nume = new THStack("s_dEtaInSeed_endcap_nume", "");
    THStack *s_dEtaInSeed_barrel_deno = new THStack("s_dEtaInSeed_barrel_deno", "");
    THStack *s_dEtaInSeed_endcap_deno = new THStack("s_dEtaInSeed_endcap_deno", "");
    THStack *s_dEtaInSeed_barrel_ctrl = new THStack("s_dEtaInSeed_barrel_ctrl", "");
    THStack *s_dEtaInSeed_endcap_ctrl = new THStack("s_dEtaInSeed_endcap_ctrl", "");
    THStack *s_dPhiIn_barrel_nume = new THStack("s_dPhiIn_barrel_nume", "");
    THStack *s_dPhiIn_endcap_nume = new THStack("s_dPhiIn_endcap_nume", "");
    THStack *s_dPhiIn_barrel_deno = new THStack("s_dPhiIn_barrel_deno", "");
    THStack *s_dPhiIn_endcap_deno = new THStack("s_dPhiIn_endcap_deno", "");
    THStack *s_dPhiIn_barrel_ctrl = new THStack("s_dPhiIn_barrel_ctrl", "");
    THStack *s_dPhiIn_endcap_ctrl = new THStack("s_dPhiIn_endcap_ctrl", "");
    THStack *s_HoverE_barrel_nume = new THStack("s_HoverE_barrel_nume", "");
    THStack *s_HoverE_endcap_nume = new THStack("s_HoverE_endcap_nume", "");
    THStack *s_HoverE_barrel_deno = new THStack("s_HoverE_barrel_deno", "");
    THStack *s_HoverE_endcap_deno = new THStack("s_HoverE_endcap_deno", "");
    THStack *s_HoverE_barrel_ctrl = new THStack("s_HoverE_barrel_ctrl", "");
    THStack *s_HoverE_endcap_ctrl = new THStack("s_HoverE_endcap_ctrl", "");
    THStack *s_InvEminusInvP_barrel_nume = new THStack("s_InvEminusInvP_barrel_nume", "");
    THStack *s_InvEminusInvP_endcap_nume = new THStack("s_InvEminusInvP_endcap_nume", "");
    THStack *s_InvEminusInvP_barrel_deno = new THStack("s_InvEminusInvP_barrel_deno", "");
    THStack *s_InvEminusInvP_endcap_deno = new THStack("s_InvEminusInvP_endcap_deno", "");
    THStack *s_InvEminusInvP_barrel_ctrl = new THStack("s_InvEminusInvP_barrel_ctrl", "");
    THStack *s_InvEminusInvP_endcap_ctrl = new THStack("s_InvEminusInvP_endcap_ctrl", "");
    THStack *s_chiso_barrel_nume = new THStack("s_chiso_barrel_nume", "");
    THStack *s_chiso_endcap_nume = new THStack("s_chiso_endcap_nume", "");
    THStack *s_chiso_barrel_deno = new THStack("s_chiso_barrel_deno", "");
    THStack *s_chiso_endcap_deno = new THStack("s_chiso_endcap_deno", "");
    THStack *s_chiso_barrel_ctrl = new THStack("s_chiso_barrel_ctrl", "");
    THStack *s_chiso_endcap_ctrl = new THStack("s_chiso_endcap_ctrl", "");
    THStack *s_nhiso_barrel_nume = new THStack("s_nhiso_barrel_nume", "");
    THStack *s_nhiso_endcap_nume = new THStack("s_nhiso_endcap_nume", "");
    THStack *s_nhiso_barrel_deno = new THStack("s_nhiso_barrel_deno", "");
    THStack *s_nhiso_endcap_deno = new THStack("s_nhiso_endcap_deno", "");
    THStack *s_nhiso_barrel_ctrl = new THStack("s_nhiso_barrel_ctrl", "");
    THStack *s_nhiso_endcap_ctrl = new THStack("s_nhiso_endcap_ctrl", "");
    THStack *s_phiso_barrel_nume = new THStack("s_phiso_barrel_nume", "");
    THStack *s_phiso_endcap_nume = new THStack("s_phiso_endcap_nume", "");
    THStack *s_phiso_barrel_deno = new THStack("s_phiso_barrel_deno", "");
    THStack *s_phiso_endcap_deno = new THStack("s_phiso_endcap_deno", "");
    THStack *s_phiso_barrel_ctrl = new THStack("s_phiso_barrel_ctrl", "");
    THStack *s_phiso_endcap_ctrl = new THStack("s_phiso_endcap_ctrl", "");
    THStack *s_chisoPU_barrel_nume = new THStack("s_chisoPU_barrel_nume", "");
    THStack *s_chisoPU_endcap_nume = new THStack("s_chisoPU_endcap_nume", "");
    THStack *s_chisoPU_barrel_deno = new THStack("s_chisoPU_barrel_deno", "");
    THStack *s_chisoPU_endcap_deno = new THStack("s_chisoPU_endcap_deno", "");
    THStack *s_chisoPU_barrel_ctrl = new THStack("s_chisoPU_barrel_ctrl", "");
    THStack *s_chisoPU_endcap_ctrl = new THStack("s_chisoPU_endcap_ctrl", "");
    THStack *s_MET = new THStack("s_MET", "");
    THStack *s_MT_barrel_nume = new THStack("s_MT_barrel_nume", "");
    THStack *s_MT_endcap_nume = new THStack("s_MT_endcap_nume", "");
    THStack *s_MT_barrel_deno = new THStack("s_MT_barrel_deno", "");
    THStack *s_MT_endcap_deno = new THStack("s_MT_endcap_deno", "");
    THStack *s_MT_barrel_ctrl = new THStack("s_MT_barrel_ctrl", "");
    THStack *s_MT_endcap_ctrl = new THStack("s_MT_endcap_ctrl", "");
    THStack *s_eta = new THStack("s_eta", "");
    THStack *s_nVTX = new THStack("s_nVTX", "");
    THStack *s_mass_test = new THStack("s_mass_test", "");

    TH1D *h_PFiso_dBeta_barrel_MC_nume[_EndOf_Data_Special], *h_PFiso_dBeta_endcap_MC_nume[_EndOf_Data_Special],
         *h_PFiso_dBeta_barrel_MC_deno[_EndOf_Data_Special], *h_PFiso_dBeta_endcap_MC_deno[_EndOf_Data_Special],
         *h_PFiso_dBeta_barrel_MC_ctrl[_EndOf_Data_Special], *h_PFiso_dBeta_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_PFiso_Rho_barrel_MC_nume[_EndOf_Data_Special], *h_PFiso_Rho_endcap_MC_nume[_EndOf_Data_Special],
         *h_PFiso_Rho_barrel_MC_deno[_EndOf_Data_Special], *h_PFiso_Rho_endcap_MC_deno[_EndOf_Data_Special],
         *h_PFiso_Rho_barrel_MC_ctrl[_EndOf_Data_Special], *h_PFiso_Rho_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_MC_nume[_EndOf_Data_Special], *h_pT_endcap_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_deno[_EndOf_Data_Special], *h_pT_endcap_MC_deno[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_SigmaIEtaIEta_barrel_MC_nume[_EndOf_Data_Special], *h_SigmaIEtaIEta_endcap_MC_nume[_EndOf_Data_Special],
         *h_SigmaIEtaIEta_barrel_MC_deno[_EndOf_Data_Special], *h_SigmaIEtaIEta_endcap_MC_deno[_EndOf_Data_Special],
         *h_SigmaIEtaIEta_barrel_MC_ctrl[_EndOf_Data_Special], *h_SigmaIEtaIEta_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_dEtaInSeed_barrel_MC_nume[_EndOf_Data_Special], *h_dEtaInSeed_endcap_MC_nume[_EndOf_Data_Special],
         *h_dEtaInSeed_barrel_MC_deno[_EndOf_Data_Special], *h_dEtaInSeed_endcap_MC_deno[_EndOf_Data_Special],
         *h_dEtaInSeed_barrel_MC_ctrl[_EndOf_Data_Special], *h_dEtaInSeed_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_dPhiIn_barrel_MC_nume[_EndOf_Data_Special], *h_dPhiIn_endcap_MC_nume[_EndOf_Data_Special],
         *h_dPhiIn_barrel_MC_deno[_EndOf_Data_Special], *h_dPhiIn_endcap_MC_deno[_EndOf_Data_Special],
         *h_dPhiIn_barrel_MC_ctrl[_EndOf_Data_Special], *h_dPhiIn_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_HoverE_barrel_MC_nume[_EndOf_Data_Special], *h_HoverE_endcap_MC_nume[_EndOf_Data_Special],
         *h_HoverE_barrel_MC_deno[_EndOf_Data_Special], *h_HoverE_endcap_MC_deno[_EndOf_Data_Special],
         *h_HoverE_barrel_MC_ctrl[_EndOf_Data_Special], *h_HoverE_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_InvEminusInvP_barrel_MC_nume[_EndOf_Data_Special], *h_InvEminusInvP_endcap_MC_nume[_EndOf_Data_Special],
         *h_InvEminusInvP_barrel_MC_deno[_EndOf_Data_Special], *h_InvEminusInvP_endcap_MC_deno[_EndOf_Data_Special],
         *h_InvEminusInvP_barrel_MC_ctrl[_EndOf_Data_Special], *h_InvEminusInvP_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_chiso_barrel_MC_nume[_EndOf_Data_Special], *h_chiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_chiso_barrel_MC_deno[_EndOf_Data_Special], *h_chiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_chiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_chiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_nhiso_barrel_MC_nume[_EndOf_Data_Special], *h_nhiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_nhiso_barrel_MC_deno[_EndOf_Data_Special], *h_nhiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_nhiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_nhiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_phiso_barrel_MC_nume[_EndOf_Data_Special], *h_phiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_phiso_barrel_MC_deno[_EndOf_Data_Special], *h_phiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_phiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_phiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_chisoPU_barrel_MC_nume[_EndOf_Data_Special], *h_chisoPU_endcap_MC_nume[_EndOf_Data_Special],
         *h_chisoPU_barrel_MC_deno[_EndOf_Data_Special], *h_chisoPU_endcap_MC_deno[_EndOf_Data_Special],
         *h_chisoPU_barrel_MC_ctrl[_EndOf_Data_Special], *h_chisoPU_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_MET_MC[_EndOf_Data_Special], *h_eta_data, *h_nVTX_data,
         *h_MT_barrel_MC_nume[_EndOf_Data_Special], *h_MT_endcap_MC_nume[_EndOf_Data_Special],
         *h_MT_barrel_MC_deno[_EndOf_Data_Special], *h_MT_endcap_MC_deno[_EndOf_Data_Special],
         *h_MT_barrel_MC_ctrl[_EndOf_Data_Special], *h_MT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_eta_MC[_EndOf_Data_Special], *h_nVTX_MC[_EndOf_Data_Special], *h_mass_test_MC[_EndOf_Data_Special],
         *h_PFiso_dBeta_barrel_data_nume, *h_PFiso_dBeta_endcap_data_nume,
         *h_PFiso_dBeta_barrel_data_deno, *h_PFiso_dBeta_endcap_data_deno,
         *h_PFiso_dBeta_barrel_data_ctrl, *h_PFiso_dBeta_endcap_data_ctrl,
         *h_PFiso_Rho_barrel_data_nume, *h_PFiso_Rho_endcap_data_nume,
         *h_PFiso_Rho_barrel_data_deno, *h_PFiso_Rho_endcap_data_deno,
         *h_PFiso_Rho_barrel_data_ctrl, *h_PFiso_Rho_endcap_data_ctrl,
         *h_pT_barrel_data_nume, *h_pT_endcap_data_nume,
         *h_pT_barrel_data_deno, *h_pT_endcap_data_deno,
         *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl,
         *h_SigmaIEtaIEta_barrel_data_nume, *h_SigmaIEtaIEta_endcap_data_nume,
         *h_SigmaIEtaIEta_barrel_data_deno, *h_SigmaIEtaIEta_endcap_data_deno,
         *h_SigmaIEtaIEta_barrel_data_ctrl, *h_SigmaIEtaIEta_endcap_data_ctrl,
         *h_dEtaInSeed_barrel_data_nume, *h_dEtaInSeed_endcap_data_nume,
         *h_dEtaInSeed_barrel_data_deno, *h_dEtaInSeed_endcap_data_deno,
         *h_dEtaInSeed_barrel_data_ctrl, *h_dEtaInSeed_endcap_data_ctrl,
         *h_dPhiIn_barrel_data_nume, *h_dPhiIn_endcap_data_nume,
         *h_dPhiIn_barrel_data_deno, *h_dPhiIn_endcap_data_deno,
         *h_dPhiIn_barrel_data_ctrl, *h_dPhiIn_endcap_data_ctrl,
         *h_HoverE_barrel_data_nume, *h_HoverE_endcap_data_nume,
         *h_HoverE_barrel_data_deno, *h_HoverE_endcap_data_deno,
         *h_HoverE_barrel_data_ctrl, *h_HoverE_endcap_data_ctrl,
         *h_InvEminusInvP_barrel_data_nume, *h_InvEminusInvP_endcap_data_nume,
         *h_InvEminusInvP_barrel_data_deno, *h_InvEminusInvP_endcap_data_deno,
         *h_InvEminusInvP_barrel_data_ctrl, *h_InvEminusInvP_endcap_data_ctrl,
         *h_chiso_barrel_data_nume, *h_chiso_endcap_data_nume,
         *h_chiso_barrel_data_deno, *h_chiso_endcap_data_deno,
         *h_chiso_barrel_data_ctrl, *h_chiso_endcap_data_ctrl,
         *h_nhiso_barrel_data_nume, *h_nhiso_endcap_data_nume,
         *h_nhiso_barrel_data_deno, *h_nhiso_endcap_data_deno,
         *h_nhiso_barrel_data_ctrl, *h_nhiso_endcap_data_ctrl,
         *h_phiso_barrel_data_nume, *h_phiso_endcap_data_nume,
         *h_phiso_barrel_data_deno, *h_phiso_endcap_data_deno,
         *h_phiso_barrel_data_ctrl, *h_phiso_endcap_data_ctrl,
         *h_chisoPU_barrel_data_nume, *h_chisoPU_endcap_data_nume,
         *h_chisoPU_barrel_data_deno, *h_chisoPU_endcap_data_deno,
         *h_chisoPU_barrel_data_ctrl, *h_chisoPU_endcap_data_ctrl,
         *h_MET_data,
         *h_MT_barrel_data_nume, *h_MT_endcap_data_nume,
         *h_MT_barrel_data_deno, *h_MT_endcap_data_deno,
         *h_MT_barrel_data_ctrl, *h_MT_endcap_data_ctrl,
         *h_mass_test_data;

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Electron/SelectedForFR_E_"+fm.Procname[pr1]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr1]+".root", "READ");
        else return;
        file->GetObject("h_PFiso_dBeta_barrel_nume", h_PFiso_dBeta_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_dBeta_endcap_nume", h_PFiso_dBeta_endcap_MC_nume[pr1]);
        file->GetObject("h_PFiso_dBeta_barrel_deno", h_PFiso_dBeta_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_dBeta_endcap_deno", h_PFiso_dBeta_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_PFiso_dBeta_barrel_MC_ctrl[pr1]);
        file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_PFiso_dBeta_endcap_MC_ctrl[pr1]);
        file->GetObject("h_PFiso_Rho_barrel_nume", h_PFiso_Rho_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_Rho_endcap_nume", h_PFiso_Rho_endcap_MC_nume[pr1]);
        file->GetObject("h_PFiso_Rho_barrel_deno", h_PFiso_Rho_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_Rho_endcap_deno", h_PFiso_Rho_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_Rho_barrel_ctrl", h_PFiso_Rho_barrel_MC_ctrl[pr1]);
        file->GetObject("h_PFiso_Rho_endcap_ctrl", h_PFiso_Rho_endcap_MC_ctrl[pr1]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr1]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr1]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_SigmaIEtaIEta_barrel_MC_nume[pr1]);
        file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_SigmaIEtaIEta_endcap_MC_nume[pr1]);
        file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_SigmaIEtaIEta_barrel_MC_deno[pr1]);
        file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_SigmaIEtaIEta_endcap_MC_deno[pr1]);
        file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]);
        file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]);
        file->GetObject("h_dEtaInSeed_barrel_nume", h_dEtaInSeed_barrel_MC_nume[pr1]);
        file->GetObject("h_dEtaInSeed_endcap_nume", h_dEtaInSeed_endcap_MC_nume[pr1]);
        file->GetObject("h_dEtaInSeed_barrel_deno", h_dEtaInSeed_barrel_MC_deno[pr1]);
        file->GetObject("h_dEtaInSeed_endcap_deno", h_dEtaInSeed_endcap_MC_deno[pr1]);
        file->GetObject("h_dEtaInSeed_barrel_ctrl", h_dEtaInSeed_barrel_MC_ctrl[pr1]);
        file->GetObject("h_dEtaInSeed_endcap_ctrl", h_dEtaInSeed_endcap_MC_ctrl[pr1]);
        file->GetObject("h_dPhiIn_barrel_nume", h_dPhiIn_barrel_MC_nume[pr1]);
        file->GetObject("h_dPhiIn_endcap_nume", h_dPhiIn_endcap_MC_nume[pr1]);
        file->GetObject("h_dPhiIn_barrel_deno", h_dPhiIn_barrel_MC_deno[pr1]);
        file->GetObject("h_dPhiIn_endcap_deno", h_dPhiIn_endcap_MC_deno[pr1]);
        file->GetObject("h_dPhiIn_barrel_ctrl", h_dPhiIn_barrel_MC_ctrl[pr1]);
        file->GetObject("h_dPhiIn_endcap_ctrl", h_dPhiIn_endcap_MC_ctrl[pr1]);
        file->GetObject("h_HoverE_barrel_nume", h_HoverE_barrel_MC_nume[pr1]);
        file->GetObject("h_HoverE_endcap_nume", h_HoverE_endcap_MC_nume[pr1]);
        file->GetObject("h_HoverE_barrel_deno", h_HoverE_barrel_MC_deno[pr1]);
        file->GetObject("h_HoverE_endcap_deno", h_HoverE_endcap_MC_deno[pr1]);
        file->GetObject("h_HoverE_barrel_ctrl", h_HoverE_barrel_MC_ctrl[pr1]);
        file->GetObject("h_HoverE_endcap_ctrl", h_HoverE_endcap_MC_ctrl[pr1]);
        file->GetObject("h_InvEminusInvP_barrel_nume", h_InvEminusInvP_barrel_MC_nume[pr1]);
        file->GetObject("h_InvEminusInvP_endcap_nume", h_InvEminusInvP_endcap_MC_nume[pr1]);
        file->GetObject("h_InvEminusInvP_barrel_deno", h_InvEminusInvP_barrel_MC_deno[pr1]);
        file->GetObject("h_InvEminusInvP_endcap_deno", h_InvEminusInvP_endcap_MC_deno[pr1]);
        file->GetObject("h_InvEminusInvP_barrel_ctrl", h_InvEminusInvP_barrel_MC_ctrl[pr1]);
        file->GetObject("h_InvEminusInvP_endcap_ctrl", h_InvEminusInvP_endcap_MC_ctrl[pr1]);
        file->GetObject("h_chiso_barrel_nume", h_chiso_barrel_MC_nume[pr1]);
        file->GetObject("h_chiso_endcap_nume", h_chiso_endcap_MC_nume[pr1]);
        file->GetObject("h_chiso_barrel_deno", h_chiso_barrel_MC_deno[pr1]);
        file->GetObject("h_chiso_endcap_deno", h_chiso_endcap_MC_deno[pr1]);
        file->GetObject("h_chiso_barrel_ctrl", h_chiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_chiso_endcap_ctrl", h_chiso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_nhiso_barrel_nume", h_nhiso_barrel_MC_nume[pr1]);
        file->GetObject("h_nhiso_endcap_nume", h_nhiso_endcap_MC_nume[pr1]);
        file->GetObject("h_nhiso_barrel_deno", h_nhiso_barrel_MC_deno[pr1]);
        file->GetObject("h_nhiso_endcap_deno", h_nhiso_endcap_MC_deno[pr1]);
        file->GetObject("h_nhiso_barrel_ctrl", h_nhiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_nhiso_endcap_ctrl", h_nhiso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_phiso_barrel_nume", h_phiso_barrel_MC_nume[pr1]);
        file->GetObject("h_phiso_endcap_nume", h_phiso_endcap_MC_nume[pr1]);
        file->GetObject("h_phiso_barrel_deno", h_phiso_barrel_MC_deno[pr1]);
        file->GetObject("h_phiso_endcap_deno", h_phiso_endcap_MC_deno[pr1]);
        file->GetObject("h_phiso_barrel_ctrl", h_phiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_phiso_endcap_ctrl", h_phiso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_chisoPU_barrel_nume", h_chisoPU_barrel_MC_nume[pr1]);
        file->GetObject("h_chisoPU_endcap_nume", h_chisoPU_endcap_MC_nume[pr1]);
        file->GetObject("h_chisoPU_barrel_deno", h_chisoPU_barrel_MC_deno[pr1]);
        file->GetObject("h_chisoPU_endcap_deno", h_chisoPU_endcap_MC_deno[pr1]);
        file->GetObject("h_chisoPU_barrel_ctrl", h_chisoPU_barrel_MC_ctrl[pr1]);
        file->GetObject("h_chisoPU_endcap_ctrl", h_chisoPU_endcap_MC_ctrl[pr1]);
        file->GetObject("h_MET", h_MET_MC[pr1]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr1]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr1]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr1]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr1]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_eta_deno", h_eta_MC[pr1]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr1]);
        file->GetObject("h_mass_test", h_mass_test_MC[pr1]);

        removeNegativeBins(h_PFiso_dBeta_barrel_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_nume[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_nume[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_deno[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_deno[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_nume[pr1]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_nume[pr1]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_deno[pr1]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_deno[pr1]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_dPhiIn_barrel_MC_nume[pr1]);
        removeNegativeBins(h_dPhiIn_endcap_MC_nume[pr1]);
        removeNegativeBins(h_dPhiIn_barrel_MC_deno[pr1]);
        removeNegativeBins(h_dPhiIn_endcap_MC_deno[pr1]);
        removeNegativeBins(h_dPhiIn_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_dPhiIn_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_HoverE_barrel_MC_nume[pr1]);
        removeNegativeBins(h_HoverE_endcap_MC_nume[pr1]);
        removeNegativeBins(h_HoverE_barrel_MC_deno[pr1]);
        removeNegativeBins(h_HoverE_endcap_MC_deno[pr1]);
        removeNegativeBins(h_HoverE_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_HoverE_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_nume[pr1]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_nume[pr1]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_deno[pr1]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_deno[pr1]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_chiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_chiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_chiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_chiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_chiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_chiso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_nhiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_nhiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_nhiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_nhiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_nhiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_nhiso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_phiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_phiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_phiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_phiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_phiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_phiso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_chisoPU_barrel_MC_nume[pr1]);
        removeNegativeBins(h_chisoPU_endcap_MC_nume[pr1]);
        removeNegativeBins(h_chisoPU_barrel_MC_deno[pr1]);
        removeNegativeBins(h_chisoPU_endcap_MC_deno[pr1]);
        removeNegativeBins(h_chisoPU_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_chisoPU_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_MET_MC[pr1]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_eta_MC[pr1]);
        removeNegativeBins(h_nVTX_MC[pr1]);
        removeNegativeBins(h_mass_test_MC[pr1]);

        Color_t color = kBlack;
        if (pr1 == _WJets || pr1 == _WJets_ext2v5) color = kRed - 2;
        if (pr1 == _VVnST) color = kMagenta - 5;
        if (pr1 == _WW) color = kMagenta - 5;
        if (pr1 == _WZ) color = kMagenta - 2;
        if (pr1 == _ZZ) color = kMagenta - 6;
        if (pr1 == _tbarW) color = kGreen - 2;
        if (pr1 == _tW) color = kGreen + 2;
        if (pr1 == _ttbar || pr1 == _ttbar_700to1000 || pr1 == _ttbar_1000toInf) color = kCyan + 2;

        h_PFiso_dBeta_barrel_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr1]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr1]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr1]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr1]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr1]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr1]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr1]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr1]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_dPhiIn_barrel_MC_nume[pr1]->SetFillColor(color);
        h_dPhiIn_endcap_MC_nume[pr1]->SetFillColor(color);
        h_dPhiIn_barrel_MC_deno[pr1]->SetFillColor(color);
        h_dPhiIn_endcap_MC_deno[pr1]->SetFillColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_HoverE_barrel_MC_nume[pr1]->SetFillColor(color);
        h_HoverE_endcap_MC_nume[pr1]->SetFillColor(color);
        h_HoverE_barrel_MC_deno[pr1]->SetFillColor(color);
        h_HoverE_endcap_MC_deno[pr1]->SetFillColor(color);
        h_HoverE_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_HoverE_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr1]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr1]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr1]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr1]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_chiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_chiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_chiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_chiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_chiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_chiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_nhiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_nhiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_nhiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_nhiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_nhiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_nhiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_phiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_phiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_phiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_phiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_phiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_phiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_chisoPU_barrel_MC_nume[pr1]->SetFillColor(color);
        h_chisoPU_endcap_MC_nume[pr1]->SetFillColor(color);
        h_chisoPU_barrel_MC_deno[pr1]->SetFillColor(color);
        h_chisoPU_endcap_MC_deno[pr1]->SetFillColor(color);
        h_chisoPU_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_chisoPU_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_MET_MC[pr1]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr1]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr1]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr1]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr1]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_eta_MC[pr1]->SetFillColor(color);
        h_nVTX_MC[pr1]->SetFillColor(color);
        h_mass_test_MC[pr1]->SetFillColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr1]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr1]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr1]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr1]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr1]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr1]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr1]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr1]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_dPhiIn_barrel_MC_nume[pr1]->SetLineColor(color);
        h_dPhiIn_endcap_MC_nume[pr1]->SetLineColor(color);
        h_dPhiIn_barrel_MC_deno[pr1]->SetLineColor(color);
        h_dPhiIn_endcap_MC_deno[pr1]->SetLineColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_HoverE_barrel_MC_nume[pr1]->SetLineColor(color);
        h_HoverE_endcap_MC_nume[pr1]->SetLineColor(color);
        h_HoverE_barrel_MC_deno[pr1]->SetLineColor(color);
        h_HoverE_endcap_MC_deno[pr1]->SetLineColor(color);
        h_HoverE_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_HoverE_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr1]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr1]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr1]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr1]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_chiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_chiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_chiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_chiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_chiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_chiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_nhiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_nhiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_nhiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_nhiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_nhiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_nhiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_phiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_phiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_phiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_phiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_phiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_phiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_chisoPU_barrel_MC_nume[pr1]->SetLineColor(color);
        h_chisoPU_endcap_MC_nume[pr1]->SetLineColor(color);
        h_chisoPU_barrel_MC_deno[pr1]->SetLineColor(color);
        h_chisoPU_endcap_MC_deno[pr1]->SetLineColor(color);
        h_chisoPU_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_chisoPU_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_MET_MC[pr1]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr1]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr1]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr1]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr1]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_eta_MC[pr1]->SetLineColor(color);
        h_nVTX_MC[pr1]->SetLineColor(color);
        h_mass_test_MC[pr1]->SetLineColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_nume[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_nume[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_deno[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_deno[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_nume[pr1]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_nume[pr1]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_deno[pr1]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_deno[pr1]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_dPhiIn_barrel_MC_nume[pr1]->SetDirectory(0);
        h_dPhiIn_endcap_MC_nume[pr1]->SetDirectory(0);
        h_dPhiIn_barrel_MC_deno[pr1]->SetDirectory(0);
        h_dPhiIn_endcap_MC_deno[pr1]->SetDirectory(0);
        h_dPhiIn_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_dPhiIn_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_HoverE_barrel_MC_nume[pr1]->SetDirectory(0);
        h_HoverE_endcap_MC_nume[pr1]->SetDirectory(0);
        h_HoverE_barrel_MC_deno[pr1]->SetDirectory(0);
        h_HoverE_endcap_MC_deno[pr1]->SetDirectory(0);
        h_HoverE_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_HoverE_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_nume[pr1]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_nume[pr1]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_deno[pr1]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_deno[pr1]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_chiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_chiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_chiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_chiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_chiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_chiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_nhiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_nhiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_nhiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_nhiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_nhiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_nhiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_phiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_phiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_phiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_phiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_phiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_phiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_chisoPU_barrel_MC_nume[pr1]->SetDirectory(0);
        h_chisoPU_endcap_MC_nume[pr1]->SetDirectory(0);
        h_chisoPU_barrel_MC_deno[pr1]->SetDirectory(0);
        h_chisoPU_endcap_MC_deno[pr1]->SetDirectory(0);
        h_chisoPU_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_chisoPU_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_MET_MC[pr1]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_eta_MC[pr1]->SetDirectory(0);
        h_nVTX_MC[pr1]->SetDirectory(0);
        h_mass_test_MC[pr1]->SetDirectory(0);

        if (pr1 == _WJets)
        {
            h_PFiso_dBeta_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_nume[pr1]->Clone("h_PFiso_dBeta_barrel_nume_WJets")));
            h_PFiso_dBeta_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_nume[pr1]->Clone("h_PFiso_dBeta_endcap_nume_WJets")));
            h_PFiso_dBeta_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_deno[pr1]->Clone("h_PFiso_dBeta_barrel_deno_WJets")));
            h_PFiso_dBeta_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_deno[pr1]->Clone("h_PFiso_dBeta_endcap_deno_WJets")));
            h_PFiso_dBeta_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_ctrl[pr1]->Clone("h_PFiso_dBeta_barrel_ctrl_WJets")));
            h_PFiso_dBeta_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_ctrl[pr1]->Clone("h_PFiso_dBeta_endcap_ctrl_WJets")));
            h_PFiso_Rho_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_nume[pr1]->Clone("h_PFiso_Rho_barrel_nume_WJets")));
            h_PFiso_Rho_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_nume[pr1]->Clone("h_PFiso_Rho_endcap_nume_WJets")));
            h_PFiso_Rho_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_deno[pr1]->Clone("h_PFiso_Rho_barrel_deno_WJets")));
            h_PFiso_Rho_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_deno[pr1]->Clone("h_PFiso_Rho_endcap_deno_WJets")));
            h_PFiso_Rho_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_ctrl[pr1]->Clone("h_PFiso_Rho_barrel_ctrl_WJets")));
            h_PFiso_Rho_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_ctrl[pr1]->Clone("h_PFiso_Rho_endcap_ctrl_WJets")));
            h_pT_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr1]->Clone("h_pT_barrel_nume_WJets")));
            h_pT_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr1]->Clone("h_pT_endcap_nume_WJets")));
            h_pT_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr1]->Clone("h_pT_barrel_deno_WJets")));
            h_pT_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr1]->Clone("h_pT_endcap_deno_WJets")));
            h_pT_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr1]->Clone("h_pT_barrel_ctrl_WJets")));
            h_pT_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr1]->Clone("h_pT_endcap_ctrl_WJets")));
            h_SigmaIEtaIEta_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_nume[pr1]->Clone("h_SigmaIEtaIEta_barrel_nume_WJets")));
            h_SigmaIEtaIEta_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_nume[pr1]->Clone("h_SigmaIEtaIEta_endcap_nume_WJets")));
            h_SigmaIEtaIEta_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_deno[pr1]->Clone("h_SigmaIEtaIEta_barrel_deno_WJets")));
            h_SigmaIEtaIEta_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_deno[pr1]->Clone("h_SigmaIEtaIEta_endcap_deno_WJets")));
            h_SigmaIEtaIEta_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]->Clone("h_SigmaIEtaIEta_barrel_ctrl_WJets")));
            h_SigmaIEtaIEta_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]->Clone("h_SigmaIEtaIEta_endcap_ctrl_WJets")));
            h_dEtaInSeed_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_nume[pr1]->Clone("h_dEtaInSeed_barrel_nume_WJets")));
            h_dEtaInSeed_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_nume[pr1]->Clone("h_dEtaInSeed_endcap_nume_WJets")));
            h_dEtaInSeed_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_deno[pr1]->Clone("h_dEtaInSeed_barrel_deno_WJets")));
            h_dEtaInSeed_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_deno[pr1]->Clone("h_dEtaInSeed_endcap_deno_WJets")));
            h_dEtaInSeed_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_ctrl[pr1]->Clone("h_dEtaInSeed_barrel_ctrl_WJets")));
            h_dEtaInSeed_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_ctrl[pr1]->Clone("h_dEtaInSeed_endcap_ctrl_WJets")));
            h_dPhiIn_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_nume[pr1]->Clone("h_dPhiIn_barrel_nume_WJets")));
            h_dPhiIn_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_nume[pr1]->Clone("h_dPhiIn_endcap_nume_WJets")));
            h_dPhiIn_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_deno[pr1]->Clone("h_dPhiIn_barrel_deno_WJets")));
            h_dPhiIn_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_deno[pr1]->Clone("h_dPhiIn_endcap_deno_WJets")));
            h_dPhiIn_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_ctrl[pr1]->Clone("h_dPhiIn_barrel_ctrl_WJets")));
            h_dPhiIn_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_ctrl[pr1]->Clone("h_dPhiIn_endcap_ctrl_WJets")));
            h_HoverE_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_nume[pr1]->Clone("h_HoverE_barrel_nume_WJets")));
            h_HoverE_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_nume[pr1]->Clone("h_HoverE_endcap_nume_WJets")));
            h_HoverE_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_deno[pr1]->Clone("h_HoverE_barrel_deno_WJets")));
            h_HoverE_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_deno[pr1]->Clone("h_HoverE_endcap_deno_WJets")));
            h_HoverE_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_ctrl[pr1]->Clone("h_HoverE_barrel_ctrl_WJets")));
            h_HoverE_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_ctrl[pr1]->Clone("h_HoverE_endcap_ctrl_WJets")));
            h_InvEminusInvP_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_nume[pr1]->Clone("h_InvEminusInvP_barrel_nume_WJets")));
            h_InvEminusInvP_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_nume[pr1]->Clone("h_InvEminusInvP_endcap_nume_WJets")));
            h_InvEminusInvP_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_deno[pr1]->Clone("h_InvEminusInvP_barrel_deno_WJets")));
            h_InvEminusInvP_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_deno[pr1]->Clone("h_InvEminusInvP_endcap_deno_WJets")));
            h_InvEminusInvP_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_ctrl[pr1]->Clone("h_InvEminusInvP_barrel_ctrl_WJets")));
            h_InvEminusInvP_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_ctrl[pr1]->Clone("h_InvEminusInvP_endcap_ctrl_WJets")));
            h_chiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_chiso_barrel_MC_nume[pr1]->Clone("h_chiso_barrel_nume_WJets")));
            h_chiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_chiso_endcap_MC_nume[pr1]->Clone("h_chiso_endcap_nume_WJets")));
            h_chiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_chiso_barrel_MC_deno[pr1]->Clone("h_chiso_barrel_deno_WJets")));
            h_chiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_chiso_endcap_MC_deno[pr1]->Clone("h_chiso_endcap_deno_WJets")));
            h_chiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_chiso_barrel_MC_ctrl[pr1]->Clone("h_chiso_barrel_ctrl_WJets")));
            h_chiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_chiso_endcap_MC_ctrl[pr1]->Clone("h_chiso_endcap_ctrl_WJets")));
            h_nhiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_nume[pr1]->Clone("h_nhiso_barrel_nume_WJets")));
            h_nhiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_nume[pr1]->Clone("h_nhiso_endcap_nume_WJets")));
            h_nhiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_deno[pr1]->Clone("h_nhiso_barrel_deno_WJets")));
            h_nhiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_deno[pr1]->Clone("h_nhiso_endcap_deno_WJets")));
            h_nhiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_ctrl[pr1]->Clone("h_nhiso_barrel_ctrl_WJets")));
            h_nhiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_ctrl[pr1]->Clone("h_nhiso_endcap_ctrl_WJets")));
            h_phiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_phiso_barrel_MC_nume[pr1]->Clone("h_phiso_barrel_nume_WJets")));
            h_phiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_phiso_endcap_MC_nume[pr1]->Clone("h_phiso_endcap_nume_WJets")));
            h_phiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_phiso_barrel_MC_deno[pr1]->Clone("h_phiso_barrel_deno_WJets")));
            h_phiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_phiso_endcap_MC_deno[pr1]->Clone("h_phiso_endcap_deno_WJets")));
            h_phiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_phiso_barrel_MC_ctrl[pr1]->Clone("h_phiso_barrel_ctrl_WJets")));
            h_phiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_phiso_endcap_MC_ctrl[pr1]->Clone("h_phiso_endcap_ctrl_WJets")));
            h_chisoPU_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_nume[pr1]->Clone("h_chisoPU_barrel_nume_WJets")));
            h_chisoPU_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_nume[pr1]->Clone("h_chisoPU_endcap_nume_WJets")));
            h_chisoPU_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_deno[pr1]->Clone("h_chisoPU_barrel_deno_WJets")));
            h_chisoPU_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_deno[pr1]->Clone("h_chisoPU_endcap_deno_WJets")));
            h_chisoPU_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_ctrl[pr1]->Clone("h_chisoPU_barrel_ctrl_WJets")));
            h_chisoPU_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_ctrl[pr1]->Clone("h_chisoPU_endcap_ctrl_WJets")));
            h_MET_MC[_WJets_Full] = ((TH1D*)(h_MET_MC[pr1]->Clone("h_MET_DY")));
            h_MT_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr1]->Clone("h_MT_barrel_nume_WJets")));
            h_MT_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr1]->Clone("h_MT_endcap_nume_WJets")));
            h_MT_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr1]->Clone("h_MT_barrel_deno_WJets")));
            h_MT_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr1]->Clone("h_MT_endcap_deno_WJets")));
            h_MT_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr1]->Clone("h_MT_barrel_ctrl_WJets")));
            h_MT_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr1]->Clone("h_MT_endcap_ctrl_WJets")));
            h_eta_MC[_WJets_Full] = ((TH1D*)(h_eta_MC[pr1]->Clone("h_eta_deno_WJets")));
            h_nVTX_MC[_WJets_Full] = ((TH1D*)(h_nVTX_MC[pr1]->Clone("h_nVTX_WJets")));
            h_mass_test_MC[_WJets_Full] = ((TH1D*)(h_mass_test_MC[pr1]->Clone("h_mass_test_WJets")));

            h_PFiso_dBeta_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_MET_MC[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_eta_MC[_WJets_Full]->SetDirectory(0);
            h_nVTX_MC[_WJets_Full]->SetDirectory(0);
            h_mass_test_MC[_WJets_Full]->SetDirectory(0);
        }
        else if (pr1 == _WJets_ext2v5)
        {
            h_PFiso_dBeta_barrel_MC_nume[_WJets_Full]->Add(h_PFiso_dBeta_barrel_MC_nume[pr1]);
            h_PFiso_dBeta_endcap_MC_nume[_WJets_Full]->Add(h_PFiso_dBeta_endcap_MC_nume[pr1]);
            h_PFiso_dBeta_barrel_MC_deno[_WJets_Full]->Add(h_PFiso_dBeta_barrel_MC_deno[pr1]);
            h_PFiso_dBeta_endcap_MC_deno[_WJets_Full]->Add(h_PFiso_dBeta_endcap_MC_deno[pr1]);
            h_PFiso_dBeta_barrel_MC_ctrl[_WJets_Full]->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr1]);
            h_PFiso_dBeta_endcap_MC_ctrl[_WJets_Full]->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr1]);
            h_PFiso_Rho_barrel_MC_nume[_WJets_Full]->Add(h_PFiso_Rho_barrel_MC_nume[pr1]);
            h_PFiso_Rho_endcap_MC_nume[_WJets_Full]->Add(h_PFiso_Rho_endcap_MC_nume[pr1]);
            h_PFiso_Rho_barrel_MC_deno[_WJets_Full]->Add(h_PFiso_Rho_barrel_MC_deno[pr1]);
            h_PFiso_Rho_endcap_MC_deno[_WJets_Full]->Add(h_PFiso_Rho_endcap_MC_deno[pr1]);
            h_PFiso_Rho_barrel_MC_ctrl[_WJets_Full]->Add(h_PFiso_Rho_barrel_MC_ctrl[pr1]);
            h_PFiso_Rho_endcap_MC_ctrl[_WJets_Full]->Add(h_PFiso_Rho_endcap_MC_ctrl[pr1]);
            h_pT_barrel_MC_nume[_WJets_Full]->Add(h_pT_barrel_MC_nume[pr1]);
            h_pT_endcap_MC_nume[_WJets_Full]->Add(h_pT_endcap_MC_nume[pr1]);
            h_pT_barrel_MC_deno[_WJets_Full]->Add(h_pT_barrel_MC_deno[pr1]);
            h_pT_endcap_MC_deno[_WJets_Full]->Add(h_pT_endcap_MC_deno[pr1]);
            h_pT_barrel_MC_ctrl[_WJets_Full]->Add(h_pT_barrel_MC_ctrl[pr1]);
            h_pT_endcap_MC_ctrl[_WJets_Full]->Add(h_pT_endcap_MC_ctrl[pr1]);
            h_SigmaIEtaIEta_barrel_MC_nume[_WJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr1]);
            h_SigmaIEtaIEta_endcap_MC_nume[_WJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr1]);
            h_SigmaIEtaIEta_barrel_MC_deno[_WJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr1]);
            h_SigmaIEtaIEta_endcap_MC_deno[_WJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr1]);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_WJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_WJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]);
            h_dEtaInSeed_barrel_MC_nume[_WJets_Full]->Add(h_dEtaInSeed_barrel_MC_nume[pr1]);
            h_dEtaInSeed_endcap_MC_nume[_WJets_Full]->Add(h_dEtaInSeed_endcap_MC_nume[pr1]);
            h_dEtaInSeed_barrel_MC_deno[_WJets_Full]->Add(h_dEtaInSeed_barrel_MC_deno[pr1]);
            h_dEtaInSeed_endcap_MC_deno[_WJets_Full]->Add(h_dEtaInSeed_endcap_MC_deno[pr1]);
            h_dEtaInSeed_barrel_MC_ctrl[_WJets_Full]->Add(h_dEtaInSeed_barrel_MC_ctrl[pr1]);
            h_dEtaInSeed_endcap_MC_ctrl[_WJets_Full]->Add(h_dEtaInSeed_endcap_MC_ctrl[pr1]);
            h_dPhiIn_barrel_MC_nume[_WJets_Full]->Add(h_dPhiIn_barrel_MC_nume[pr1]);
            h_dPhiIn_endcap_MC_nume[_WJets_Full]->Add(h_dPhiIn_endcap_MC_nume[pr1]);
            h_dPhiIn_barrel_MC_deno[_WJets_Full]->Add(h_dPhiIn_barrel_MC_deno[pr1]);
            h_dPhiIn_endcap_MC_deno[_WJets_Full]->Add(h_dPhiIn_endcap_MC_deno[pr1]);
            h_dPhiIn_barrel_MC_ctrl[_WJets_Full]->Add(h_dPhiIn_barrel_MC_ctrl[pr1]);
            h_dPhiIn_endcap_MC_ctrl[_WJets_Full]->Add(h_dPhiIn_endcap_MC_ctrl[pr1]);
            h_HoverE_barrel_MC_nume[_WJets_Full]->Add(h_HoverE_barrel_MC_nume[pr1]);
            h_HoverE_endcap_MC_nume[_WJets_Full]->Add(h_HoverE_endcap_MC_nume[pr1]);
            h_HoverE_barrel_MC_deno[_WJets_Full]->Add(h_HoverE_barrel_MC_deno[pr1]);
            h_HoverE_endcap_MC_deno[_WJets_Full]->Add(h_HoverE_endcap_MC_deno[pr1]);
            h_HoverE_barrel_MC_ctrl[_WJets_Full]->Add(h_HoverE_barrel_MC_ctrl[pr1]);
            h_HoverE_endcap_MC_ctrl[_WJets_Full]->Add(h_HoverE_endcap_MC_ctrl[pr1]);
            h_InvEminusInvP_barrel_MC_nume[_WJets_Full]->Add(h_InvEminusInvP_barrel_MC_nume[pr1]);
            h_InvEminusInvP_endcap_MC_nume[_WJets_Full]->Add(h_InvEminusInvP_endcap_MC_nume[pr1]);
            h_InvEminusInvP_barrel_MC_deno[_WJets_Full]->Add(h_InvEminusInvP_barrel_MC_deno[pr1]);
            h_InvEminusInvP_endcap_MC_deno[_WJets_Full]->Add(h_InvEminusInvP_endcap_MC_deno[pr1]);
            h_InvEminusInvP_barrel_MC_ctrl[_WJets_Full]->Add(h_InvEminusInvP_barrel_MC_ctrl[pr1]);
            h_InvEminusInvP_endcap_MC_ctrl[_WJets_Full]->Add(h_InvEminusInvP_endcap_MC_ctrl[pr1]);
            h_chiso_barrel_MC_nume[_WJets_Full]->Add(h_chiso_barrel_MC_nume[pr1]);
            h_chiso_endcap_MC_nume[_WJets_Full]->Add(h_chiso_endcap_MC_nume[pr1]);
            h_chiso_barrel_MC_deno[_WJets_Full]->Add(h_chiso_barrel_MC_deno[pr1]);
            h_chiso_endcap_MC_deno[_WJets_Full]->Add(h_chiso_endcap_MC_deno[pr1]);
            h_chiso_barrel_MC_ctrl[_WJets_Full]->Add(h_chiso_barrel_MC_ctrl[pr1]);
            h_chiso_endcap_MC_ctrl[_WJets_Full]->Add(h_chiso_endcap_MC_ctrl[pr1]);
            h_nhiso_barrel_MC_nume[_WJets_Full]->Add(h_nhiso_barrel_MC_nume[pr1]);
            h_nhiso_endcap_MC_nume[_WJets_Full]->Add(h_nhiso_endcap_MC_nume[pr1]);
            h_nhiso_barrel_MC_deno[_WJets_Full]->Add(h_nhiso_barrel_MC_deno[pr1]);
            h_nhiso_endcap_MC_deno[_WJets_Full]->Add(h_nhiso_endcap_MC_deno[pr1]);
            h_nhiso_barrel_MC_ctrl[_WJets_Full]->Add(h_nhiso_barrel_MC_ctrl[pr1]);
            h_nhiso_endcap_MC_ctrl[_WJets_Full]->Add(h_nhiso_endcap_MC_ctrl[pr1]);
            h_phiso_barrel_MC_nume[_WJets_Full]->Add(h_phiso_barrel_MC_nume[pr1]);
            h_phiso_endcap_MC_nume[_WJets_Full]->Add(h_phiso_endcap_MC_nume[pr1]);
            h_phiso_barrel_MC_deno[_WJets_Full]->Add(h_phiso_barrel_MC_deno[pr1]);
            h_phiso_endcap_MC_deno[_WJets_Full]->Add(h_phiso_endcap_MC_deno[pr1]);
            h_phiso_barrel_MC_ctrl[_WJets_Full]->Add(h_phiso_barrel_MC_ctrl[pr1]);
            h_phiso_endcap_MC_ctrl[_WJets_Full]->Add(h_phiso_endcap_MC_ctrl[pr1]);
            h_chisoPU_barrel_MC_nume[_WJets_Full]->Add(h_chisoPU_barrel_MC_nume[pr1]);
            h_chisoPU_endcap_MC_nume[_WJets_Full]->Add(h_chisoPU_endcap_MC_nume[pr1]);
            h_chisoPU_barrel_MC_deno[_WJets_Full]->Add(h_chisoPU_barrel_MC_deno[pr1]);
            h_chisoPU_endcap_MC_deno[_WJets_Full]->Add(h_chisoPU_endcap_MC_deno[pr1]);
            h_chisoPU_barrel_MC_ctrl[_WJets_Full]->Add(h_chisoPU_barrel_MC_ctrl[pr1]);
            h_chisoPU_endcap_MC_ctrl[_WJets_Full]->Add(h_chisoPU_endcap_MC_ctrl[pr1]);
            h_MET_MC[_WJets_Full]->Add(h_MET_MC[pr1]);
            h_MT_barrel_MC_nume[_WJets_Full]->Add(h_MT_barrel_MC_nume[pr1]);
            h_MT_endcap_MC_nume[_WJets_Full]->Add(h_MT_endcap_MC_nume[pr1]);
            h_MT_barrel_MC_deno[_WJets_Full]->Add(h_MT_barrel_MC_deno[pr1]);
            h_MT_endcap_MC_deno[_WJets_Full]->Add(h_MT_endcap_MC_deno[pr1]);
            h_MT_barrel_MC_ctrl[_WJets_Full]->Add(h_MT_barrel_MC_ctrl[pr1]);
            h_MT_endcap_MC_ctrl[_WJets_Full]->Add(h_MT_endcap_MC_ctrl[pr1]);
            h_eta_MC[_WJets_Full]->Add(h_eta_MC[pr1]);
            h_nVTX_MC[_WJets_Full]->Add(h_nVTX_MC[pr1]);
            h_mass_test_MC[_WJets_Full]->Add(h_mass_test_MC[pr1]);
        }

        s_PFiso_dBeta_barrel_nume->Add(h_PFiso_dBeta_barrel_MC_nume[pr1]);
        s_PFiso_dBeta_endcap_nume->Add(h_PFiso_dBeta_endcap_MC_nume[pr1]);
        s_PFiso_dBeta_barrel_deno->Add(h_PFiso_dBeta_barrel_MC_deno[pr1]);
        s_PFiso_dBeta_endcap_deno->Add(h_PFiso_dBeta_endcap_MC_deno[pr1]);
        s_PFiso_dBeta_barrel_ctrl->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr1]);
        s_PFiso_dBeta_endcap_ctrl->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr1]);
        s_PFiso_Rho_barrel_nume->Add(h_PFiso_Rho_barrel_MC_nume[pr1]);
        s_PFiso_Rho_endcap_nume->Add(h_PFiso_Rho_endcap_MC_nume[pr1]);
        s_PFiso_Rho_barrel_deno->Add(h_PFiso_Rho_barrel_MC_deno[pr1]);
        s_PFiso_Rho_endcap_deno->Add(h_PFiso_Rho_endcap_MC_deno[pr1]);
        s_PFiso_Rho_barrel_ctrl->Add(h_PFiso_Rho_barrel_MC_ctrl[pr1]);
        s_PFiso_Rho_endcap_ctrl->Add(h_PFiso_Rho_endcap_MC_ctrl[pr1]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr1]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr1]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr1]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr1]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr1]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr1]);
        s_SigmaIEtaIEta_barrel_nume->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr1]);
        s_SigmaIEtaIEta_endcap_nume->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr1]);
        s_SigmaIEtaIEta_barrel_deno->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr1]);
        s_SigmaIEtaIEta_endcap_deno->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr1]);
        s_SigmaIEtaIEta_barrel_ctrl->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr1]);
        s_SigmaIEtaIEta_endcap_ctrl->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr1]);
        s_dEtaInSeed_barrel_nume->Add(h_dEtaInSeed_barrel_MC_nume[pr1]);
        s_dEtaInSeed_endcap_nume->Add(h_dEtaInSeed_endcap_MC_nume[pr1]);
        s_dEtaInSeed_barrel_deno->Add(h_dEtaInSeed_barrel_MC_deno[pr1]);
        s_dEtaInSeed_endcap_deno->Add(h_dEtaInSeed_endcap_MC_deno[pr1]);
        s_dEtaInSeed_barrel_ctrl->Add(h_dEtaInSeed_barrel_MC_ctrl[pr1]);
        s_dEtaInSeed_endcap_ctrl->Add(h_dEtaInSeed_endcap_MC_ctrl[pr1]);
        s_dPhiIn_barrel_nume->Add(h_dPhiIn_barrel_MC_nume[pr1]);
        s_dPhiIn_endcap_nume->Add(h_dPhiIn_endcap_MC_nume[pr1]);
        s_dPhiIn_barrel_deno->Add(h_dPhiIn_barrel_MC_deno[pr1]);
        s_dPhiIn_endcap_deno->Add(h_dPhiIn_endcap_MC_deno[pr1]);
        s_dPhiIn_barrel_ctrl->Add(h_dPhiIn_barrel_MC_ctrl[pr1]);
        s_dPhiIn_endcap_ctrl->Add(h_dPhiIn_endcap_MC_ctrl[pr1]);
        s_HoverE_barrel_nume->Add(h_HoverE_barrel_MC_nume[pr1]);
        s_HoverE_endcap_nume->Add(h_HoverE_endcap_MC_nume[pr1]);
        s_HoverE_barrel_deno->Add(h_HoverE_barrel_MC_deno[pr1]);
        s_HoverE_endcap_deno->Add(h_HoverE_endcap_MC_deno[pr1]);
        s_HoverE_barrel_ctrl->Add(h_HoverE_barrel_MC_ctrl[pr1]);
        s_HoverE_endcap_ctrl->Add(h_HoverE_endcap_MC_ctrl[pr1]);
        s_InvEminusInvP_barrel_nume->Add(h_InvEminusInvP_barrel_MC_nume[pr1]);
        s_InvEminusInvP_endcap_nume->Add(h_InvEminusInvP_endcap_MC_nume[pr1]);
        s_InvEminusInvP_barrel_deno->Add(h_InvEminusInvP_barrel_MC_deno[pr1]);
        s_InvEminusInvP_endcap_deno->Add(h_InvEminusInvP_endcap_MC_deno[pr1]);
        s_InvEminusInvP_barrel_ctrl->Add(h_InvEminusInvP_barrel_MC_ctrl[pr1]);
        s_InvEminusInvP_endcap_ctrl->Add(h_InvEminusInvP_endcap_MC_ctrl[pr1]);
        s_chiso_barrel_nume->Add(h_chiso_barrel_MC_nume[pr1]);
        s_chiso_endcap_nume->Add(h_chiso_endcap_MC_nume[pr1]);
        s_chiso_barrel_deno->Add(h_chiso_barrel_MC_deno[pr1]);
        s_chiso_endcap_deno->Add(h_chiso_endcap_MC_deno[pr1]);
        s_chiso_barrel_ctrl->Add(h_chiso_barrel_MC_ctrl[pr1]);
        s_chiso_endcap_ctrl->Add(h_chiso_endcap_MC_ctrl[pr1]);
        s_nhiso_barrel_nume->Add(h_nhiso_barrel_MC_nume[pr1]);
        s_nhiso_endcap_nume->Add(h_nhiso_endcap_MC_nume[pr1]);
        s_nhiso_barrel_deno->Add(h_nhiso_barrel_MC_deno[pr1]);
        s_nhiso_endcap_deno->Add(h_nhiso_endcap_MC_deno[pr1]);
        s_nhiso_barrel_ctrl->Add(h_nhiso_barrel_MC_ctrl[pr1]);
        s_nhiso_endcap_ctrl->Add(h_nhiso_endcap_MC_ctrl[pr1]);
        s_phiso_barrel_nume->Add(h_phiso_barrel_MC_nume[pr1]);
        s_phiso_endcap_nume->Add(h_phiso_endcap_MC_nume[pr1]);
        s_phiso_barrel_deno->Add(h_phiso_barrel_MC_deno[pr1]);
        s_phiso_endcap_deno->Add(h_phiso_endcap_MC_deno[pr1]);
        s_phiso_barrel_ctrl->Add(h_phiso_barrel_MC_ctrl[pr1]);
        s_phiso_endcap_ctrl->Add(h_phiso_endcap_MC_ctrl[pr1]);
        s_chisoPU_barrel_nume->Add(h_chisoPU_barrel_MC_nume[pr1]);
        s_chisoPU_endcap_nume->Add(h_chisoPU_endcap_MC_nume[pr1]);
        s_chisoPU_barrel_deno->Add(h_chisoPU_barrel_MC_deno[pr1]);
        s_chisoPU_endcap_deno->Add(h_chisoPU_endcap_MC_deno[pr1]);
        s_chisoPU_barrel_ctrl->Add(h_chisoPU_barrel_MC_ctrl[pr1]);
        s_chisoPU_endcap_ctrl->Add(h_chisoPU_endcap_MC_ctrl[pr1]);
        s_MET->Add(h_MET_MC[pr1]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr1]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr1]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr1]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr1]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr1]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr1]);
        s_eta->Add(h_eta_MC[pr1]);
        s_nVTX->Add(h_nVTX_MC[pr1]);
        s_mass_test->Add(h_mass_test_MC[pr1]);

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
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_PFiso_dBeta_barrel_nume", h_PFiso_dBeta_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_nume", h_PFiso_dBeta_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_deno", h_PFiso_dBeta_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_deno", h_PFiso_dBeta_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_barrel_nume", h_PFiso_Rho_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_endcap_nume", h_PFiso_Rho_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_barrel_deno", h_PFiso_Rho_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_endcap_deno", h_PFiso_Rho_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_barrel_ctrl", h_PFiso_Rho_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_endcap_ctrl", h_PFiso_Rho_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_barrel_nume", h_dEtaInSeed_barrel_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_endcap_nume", h_dEtaInSeed_endcap_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_barrel_deno", h_dEtaInSeed_barrel_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_endcap_deno", h_dEtaInSeed_endcap_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_barrel_ctrl", h_dEtaInSeed_barrel_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_endcap_ctrl", h_dEtaInSeed_endcap_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_barrel_nume", h_dPhiIn_barrel_MC_nume[pr]);
        file->GetObject("h_dPhiIn_endcap_nume", h_dPhiIn_endcap_MC_nume[pr]);
        file->GetObject("h_dPhiIn_barrel_deno", h_dPhiIn_barrel_MC_deno[pr]);
        file->GetObject("h_dPhiIn_endcap_deno", h_dPhiIn_endcap_MC_deno[pr]);
        file->GetObject("h_dPhiIn_barrel_ctrl", h_dPhiIn_barrel_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_endcap_ctrl", h_dPhiIn_endcap_MC_ctrl[pr]);
        file->GetObject("h_HoverE_barrel_nume", h_HoverE_barrel_MC_nume[pr]);
        file->GetObject("h_HoverE_endcap_nume", h_HoverE_endcap_MC_nume[pr]);
        file->GetObject("h_HoverE_barrel_deno", h_HoverE_barrel_MC_deno[pr]);
        file->GetObject("h_HoverE_endcap_deno", h_HoverE_endcap_MC_deno[pr]);
        file->GetObject("h_HoverE_barrel_ctrl", h_HoverE_barrel_MC_ctrl[pr]);
        file->GetObject("h_HoverE_endcap_ctrl", h_HoverE_endcap_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_barrel_nume", h_InvEminusInvP_barrel_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_endcap_nume", h_InvEminusInvP_endcap_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_barrel_deno", h_InvEminusInvP_barrel_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_endcap_deno", h_InvEminusInvP_endcap_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_barrel_ctrl", h_InvEminusInvP_barrel_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_endcap_ctrl", h_InvEminusInvP_endcap_MC_ctrl[pr]);
        file->GetObject("h_chiso_barrel_nume", h_chiso_barrel_MC_nume[pr]);
        file->GetObject("h_chiso_endcap_nume", h_chiso_endcap_MC_nume[pr]);
        file->GetObject("h_chiso_barrel_deno", h_chiso_barrel_MC_deno[pr]);
        file->GetObject("h_chiso_endcap_deno", h_chiso_endcap_MC_deno[pr]);
        file->GetObject("h_chiso_barrel_ctrl", h_chiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_chiso_endcap_ctrl", h_chiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_nhiso_barrel_nume", h_nhiso_barrel_MC_nume[pr]);
        file->GetObject("h_nhiso_endcap_nume", h_nhiso_endcap_MC_nume[pr]);
        file->GetObject("h_nhiso_barrel_deno", h_nhiso_barrel_MC_deno[pr]);
        file->GetObject("h_nhiso_endcap_deno", h_nhiso_endcap_MC_deno[pr]);
        file->GetObject("h_nhiso_barrel_ctrl", h_nhiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_nhiso_endcap_ctrl", h_nhiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_phiso_barrel_nume", h_phiso_barrel_MC_nume[pr]);
        file->GetObject("h_phiso_endcap_nume", h_phiso_endcap_MC_nume[pr]);
        file->GetObject("h_phiso_barrel_deno", h_phiso_barrel_MC_deno[pr]);
        file->GetObject("h_phiso_endcap_deno", h_phiso_endcap_MC_deno[pr]);
        file->GetObject("h_phiso_barrel_ctrl", h_phiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_phiso_endcap_ctrl", h_phiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_barrel_nume", h_chisoPU_barrel_MC_nume[pr]);
        file->GetObject("h_chisoPU_endcap_nume", h_chisoPU_endcap_MC_nume[pr]);
        file->GetObject("h_chisoPU_barrel_deno", h_chisoPU_barrel_MC_deno[pr]);
        file->GetObject("h_chisoPU_endcap_deno", h_chisoPU_endcap_MC_deno[pr]);
        file->GetObject("h_chisoPU_barrel_ctrl", h_chisoPU_barrel_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_endcap_ctrl", h_chisoPU_endcap_MC_ctrl[pr]);
        file->GetObject("h_MET", h_MET_MC[pr]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        file->GetObject("h_mass_test", h_mass_test_MC[pr]);

        removeNegativeBins(h_PFiso_dBeta_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_nume[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_nume[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_deno[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_deno[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_chiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_chiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_chiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_chiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_phiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_phiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_phiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_phiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_MET_MC[pr]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);
        removeNegativeBins(h_mass_test_MC[pr]);

        if (pr == _DY_10to50)
        {
            h_PFiso_dBeta_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_nume[pr]->Clone("h_PFiso_dBeta_barrel_nume_DY")));
            h_PFiso_dBeta_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_nume[pr]->Clone("h_PFiso_dBeta_endcap_nume_DY")));
            h_PFiso_dBeta_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_deno[pr]->Clone("h_PFiso_dBeta_barrel_deno_DY")));
            h_PFiso_dBeta_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_deno[pr]->Clone("h_PFiso_dBeta_endcap_deno_DY")));
            h_PFiso_dBeta_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_ctrl[pr]->Clone("h_PFiso_dBeta_barrel_ctrl_DY")));
            h_PFiso_dBeta_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_ctrl[pr]->Clone("h_PFiso_dBeta_endcap_ctrl_DY")));
            h_PFiso_Rho_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_nume[pr]->Clone("h_PFiso_Rho_barrel_nume_DY")));
            h_PFiso_Rho_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_nume[pr]->Clone("h_PFiso_Rho_endcap_nume_DY")));
            h_PFiso_Rho_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_deno[pr]->Clone("h_PFiso_Rho_barrel_deno_DY")));
            h_PFiso_Rho_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_deno[pr]->Clone("h_PFiso_Rho_endcap_deno_DY")));
            h_PFiso_Rho_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_ctrl[pr]->Clone("h_PFiso_Rho_barrel_ctrl_DY")));
            h_PFiso_Rho_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_ctrl[pr]->Clone("h_PFiso_Rho_endcap_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_nume_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_nume_DY")));
            h_pT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_deno_DY")));
            h_pT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_deno_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_ctrl_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_ctrl_DY")));
            h_SigmaIEtaIEta_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_nume[pr]->Clone("h_SigmaIEtaIEta_barrel_nume_DY")));
            h_SigmaIEtaIEta_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_nume[pr]->Clone("h_SigmaIEtaIEta_endcap_nume_DY")));
            h_SigmaIEtaIEta_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_deno[pr]->Clone("h_SigmaIEtaIEta_barrel_deno_DY")));
            h_SigmaIEtaIEta_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_deno[pr]->Clone("h_SigmaIEtaIEta_endcap_deno_DY")));
            h_SigmaIEtaIEta_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_barrel_ctrl_DY")));
            h_SigmaIEtaIEta_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_endcap_ctrl_DY")));
            h_dEtaInSeed_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_nume[pr]->Clone("h_dEtaInSeed_barrel_nume_DY")));
            h_dEtaInSeed_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_nume[pr]->Clone("h_dEtaInSeed_endcap_nume_DY")));
            h_dEtaInSeed_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_deno[pr]->Clone("h_dEtaInSeed_barrel_deno_DY")));
            h_dEtaInSeed_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_deno[pr]->Clone("h_dEtaInSeed_endcap_deno_DY")));
            h_dEtaInSeed_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_ctrl[pr]->Clone("h_dEtaInSeed_barrel_ctrl_DY")));
            h_dEtaInSeed_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_ctrl[pr]->Clone("h_dEtaInSeed_endcap_ctrl_DY")));
            h_dPhiIn_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_nume[pr]->Clone("h_dPhiIn_barrel_nume_DY")));
            h_dPhiIn_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_nume[pr]->Clone("h_dPhiIn_endcap_nume_DY")));
            h_dPhiIn_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_deno[pr]->Clone("h_dPhiIn_barrel_deno_DY")));
            h_dPhiIn_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_deno[pr]->Clone("h_dPhiIn_endcap_deno_DY")));
            h_dPhiIn_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_ctrl[pr]->Clone("h_dPhiIn_barrel_ctrl_DY")));
            h_dPhiIn_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_ctrl[pr]->Clone("h_dPhiIn_endcap_ctrl_DY")));
            h_HoverE_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_HoverE_barrel_MC_nume[pr]->Clone("h_HoverE_barrel_nume_DY")));
            h_HoverE_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_HoverE_endcap_MC_nume[pr]->Clone("h_HoverE_endcap_nume_DY")));
            h_HoverE_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_HoverE_barrel_MC_deno[pr]->Clone("h_HoverE_barrel_deno_DY")));
            h_HoverE_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_HoverE_endcap_MC_deno[pr]->Clone("h_HoverE_endcap_deno_DY")));
            h_HoverE_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_HoverE_barrel_MC_ctrl[pr]->Clone("h_HoverE_barrel_ctrl_DY")));
            h_HoverE_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_HoverE_endcap_MC_ctrl[pr]->Clone("h_HoverE_endcap_ctrl_DY")));
            h_InvEminusInvP_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_nume[pr]->Clone("h_InvEminusInvP_barrel_nume_DY")));
            h_InvEminusInvP_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_nume[pr]->Clone("h_InvEminusInvP_endcap_nume_DY")));
            h_InvEminusInvP_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_deno[pr]->Clone("h_InvEminusInvP_barrel_deno_DY")));
            h_InvEminusInvP_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_deno[pr]->Clone("h_InvEminusInvP_endcap_deno_DY")));
            h_InvEminusInvP_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_ctrl[pr]->Clone("h_InvEminusInvP_barrel_ctrl_DY")));
            h_InvEminusInvP_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_ctrl[pr]->Clone("h_InvEminusInvP_endcap_ctrl_DY")));
            h_chiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_chiso_barrel_MC_nume[pr]->Clone("h_chiso_barrel_nume_DY")));
            h_chiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_chiso_endcap_MC_nume[pr]->Clone("h_chiso_endcap_nume_DY")));
            h_chiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_chiso_barrel_MC_deno[pr]->Clone("h_chiso_barrel_deno_DY")));
            h_chiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_chiso_endcap_MC_deno[pr]->Clone("h_chiso_endcap_deno_DY")));
            h_chiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_chiso_barrel_MC_ctrl[pr]->Clone("h_chiso_barrel_ctrl_DY")));
            h_chiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_chiso_endcap_MC_ctrl[pr]->Clone("h_chiso_endcap_ctrl_DY")));
            h_nhiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_nhiso_barrel_MC_nume[pr]->Clone("h_nhiso_barrel_nume_DY")));
            h_nhiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_nhiso_endcap_MC_nume[pr]->Clone("h_nhiso_endcap_nume_DY")));
            h_nhiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_nhiso_barrel_MC_deno[pr]->Clone("h_nhiso_barrel_deno_DY")));
            h_nhiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_nhiso_endcap_MC_deno[pr]->Clone("h_nhiso_endcap_deno_DY")));
            h_nhiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_nhiso_barrel_MC_ctrl[pr]->Clone("h_nhiso_barrel_ctrl_DY")));
            h_nhiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_nhiso_endcap_MC_ctrl[pr]->Clone("h_nhiso_endcap_ctrl_DY")));
            h_phiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_phiso_barrel_MC_nume[pr]->Clone("h_phiso_barrel_nume_DY")));
            h_phiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_phiso_endcap_MC_nume[pr]->Clone("h_phiso_endcap_nume_DY")));
            h_phiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_phiso_barrel_MC_deno[pr]->Clone("h_phiso_barrel_deno_DY")));
            h_phiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_phiso_endcap_MC_deno[pr]->Clone("h_phiso_endcap_deno_DY")));
            h_phiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_phiso_barrel_MC_ctrl[pr]->Clone("h_phiso_barrel_ctrl_DY")));
            h_phiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_phiso_endcap_MC_ctrl[pr]->Clone("h_phiso_endcap_ctrl_DY")));
            h_chisoPU_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_chisoPU_barrel_MC_nume[pr]->Clone("h_chisoPU_barrel_nume_DY")));
            h_chisoPU_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_chisoPU_endcap_MC_nume[pr]->Clone("h_chisoPU_endcap_nume_DY")));
            h_chisoPU_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_chisoPU_barrel_MC_deno[pr]->Clone("h_chisoPU_barrel_deno_DY")));
            h_chisoPU_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_chisoPU_endcap_MC_deno[pr]->Clone("h_chisoPU_endcap_deno_DY")));
            h_chisoPU_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_chisoPU_barrel_MC_ctrl[pr]->Clone("h_chisoPU_barrel_ctrl_DY")));
            h_chisoPU_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_chisoPU_endcap_MC_ctrl[pr]->Clone("h_chisoPU_endcap_ctrl_DY")));
            h_MET_MC[_DY_Full] = ((TH1D*)(h_MET_MC[pr]->Clone("h_MET_DY")));
            h_MT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr]->Clone("h_MT_barrel_nume_DY")));
            h_MT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr]->Clone("h_MT_endcap_nume_DY")));
            h_MT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr]->Clone("h_MT_barrel_deno_DY")));
            h_MT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr]->Clone("h_MT_endcap_deno_DY")));
            h_MT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr]->Clone("h_MT_barrel_ctrl_DY")));
            h_MT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr]->Clone("h_MT_endcap_ctrl_DY")));
            h_eta_MC[_DY_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_DY")));
            h_nVTX_MC[_DY_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_DY")));
            h_mass_test_MC[_DY_Full] = ((TH1D*)(h_mass_test_MC[pr]->Clone("h_mass_test_DY")));

            h_PFiso_dBeta_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_chiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_chiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_chiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_chiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_chiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_chiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_phiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_phiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_phiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_phiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_phiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_phiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_MET_MC[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_MC[_DY_Full]->SetDirectory(0);
            h_nVTX_MC[_DY_Full]->SetDirectory(0);
            h_mass_test_MC[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_dBeta_barrel_MC_nume[_DY_Full]->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
            h_PFiso_dBeta_endcap_MC_nume[_DY_Full]->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
            h_PFiso_dBeta_barrel_MC_deno[_DY_Full]->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
            h_PFiso_dBeta_endcap_MC_deno[_DY_Full]->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
            h_PFiso_dBeta_barrel_MC_ctrl[_DY_Full]->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
            h_PFiso_dBeta_endcap_MC_ctrl[_DY_Full]->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
            h_PFiso_Rho_barrel_MC_nume[_DY_Full]->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
            h_PFiso_Rho_endcap_MC_nume[_DY_Full]->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
            h_PFiso_Rho_barrel_MC_deno[_DY_Full]->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
            h_PFiso_Rho_endcap_MC_deno[_DY_Full]->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
            h_PFiso_Rho_barrel_MC_ctrl[_DY_Full]->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
            h_PFiso_Rho_endcap_MC_ctrl[_DY_Full]->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_MC_deno[_DY_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_DY_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_SigmaIEtaIEta_barrel_MC_nume[_DY_Full]->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
            h_SigmaIEtaIEta_endcap_MC_nume[_DY_Full]->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
            h_SigmaIEtaIEta_barrel_MC_deno[_DY_Full]->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
            h_SigmaIEtaIEta_endcap_MC_deno[_DY_Full]->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_DY_Full]->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_DY_Full]->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
            h_dEtaInSeed_barrel_MC_nume[_DY_Full]->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
            h_dEtaInSeed_endcap_MC_nume[_DY_Full]->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
            h_dEtaInSeed_barrel_MC_deno[_DY_Full]->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
            h_dEtaInSeed_endcap_MC_deno[_DY_Full]->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
            h_dEtaInSeed_barrel_MC_ctrl[_DY_Full]->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
            h_dEtaInSeed_endcap_MC_ctrl[_DY_Full]->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
            h_dPhiIn_barrel_MC_nume[_DY_Full]->Add(h_dPhiIn_barrel_MC_nume[pr]);
            h_dPhiIn_endcap_MC_nume[_DY_Full]->Add(h_dPhiIn_endcap_MC_nume[pr]);
            h_dPhiIn_barrel_MC_deno[_DY_Full]->Add(h_dPhiIn_barrel_MC_deno[pr]);
            h_dPhiIn_endcap_MC_deno[_DY_Full]->Add(h_dPhiIn_endcap_MC_deno[pr]);
            h_dPhiIn_barrel_MC_ctrl[_DY_Full]->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
            h_dPhiIn_endcap_MC_ctrl[_DY_Full]->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
            h_HoverE_barrel_MC_nume[_DY_Full]->Add(h_HoverE_barrel_MC_nume[pr]);
            h_HoverE_endcap_MC_nume[_DY_Full]->Add(h_HoverE_endcap_MC_nume[pr]);
            h_HoverE_barrel_MC_deno[_DY_Full]->Add(h_HoverE_barrel_MC_deno[pr]);
            h_HoverE_endcap_MC_deno[_DY_Full]->Add(h_HoverE_endcap_MC_deno[pr]);
            h_HoverE_barrel_MC_ctrl[_DY_Full]->Add(h_HoverE_barrel_MC_ctrl[pr]);
            h_HoverE_endcap_MC_ctrl[_DY_Full]->Add(h_HoverE_endcap_MC_ctrl[pr]);
            h_InvEminusInvP_barrel_MC_nume[_DY_Full]->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
            h_InvEminusInvP_endcap_MC_nume[_DY_Full]->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
            h_InvEminusInvP_barrel_MC_deno[_DY_Full]->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
            h_InvEminusInvP_endcap_MC_deno[_DY_Full]->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
            h_InvEminusInvP_barrel_MC_ctrl[_DY_Full]->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
            h_InvEminusInvP_endcap_MC_ctrl[_DY_Full]->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
            h_chiso_barrel_MC_nume[_DY_Full]->Add(h_chiso_barrel_MC_nume[pr]);
            h_chiso_endcap_MC_nume[_DY_Full]->Add(h_chiso_endcap_MC_nume[pr]);
            h_chiso_barrel_MC_deno[_DY_Full]->Add(h_chiso_barrel_MC_deno[pr]);
            h_chiso_endcap_MC_deno[_DY_Full]->Add(h_chiso_endcap_MC_deno[pr]);
            h_chiso_barrel_MC_ctrl[_DY_Full]->Add(h_chiso_barrel_MC_ctrl[pr]);
            h_chiso_endcap_MC_ctrl[_DY_Full]->Add(h_chiso_endcap_MC_ctrl[pr]);
            h_nhiso_barrel_MC_nume[_DY_Full]->Add(h_nhiso_barrel_MC_nume[pr]);
            h_nhiso_endcap_MC_nume[_DY_Full]->Add(h_nhiso_endcap_MC_nume[pr]);
            h_nhiso_barrel_MC_deno[_DY_Full]->Add(h_nhiso_barrel_MC_deno[pr]);
            h_nhiso_endcap_MC_deno[_DY_Full]->Add(h_nhiso_endcap_MC_deno[pr]);
            h_nhiso_barrel_MC_ctrl[_DY_Full]->Add(h_nhiso_barrel_MC_ctrl[pr]);
            h_nhiso_endcap_MC_ctrl[_DY_Full]->Add(h_nhiso_endcap_MC_ctrl[pr]);
            h_phiso_barrel_MC_nume[_DY_Full]->Add(h_phiso_barrel_MC_nume[pr]);
            h_phiso_endcap_MC_nume[_DY_Full]->Add(h_phiso_endcap_MC_nume[pr]);
            h_phiso_barrel_MC_deno[_DY_Full]->Add(h_phiso_barrel_MC_deno[pr]);
            h_phiso_endcap_MC_deno[_DY_Full]->Add(h_phiso_endcap_MC_deno[pr]);
            h_phiso_barrel_MC_ctrl[_DY_Full]->Add(h_phiso_barrel_MC_ctrl[pr]);
            h_phiso_endcap_MC_ctrl[_DY_Full]->Add(h_phiso_endcap_MC_ctrl[pr]);
            h_chisoPU_barrel_MC_nume[_DY_Full]->Add(h_chisoPU_barrel_MC_nume[pr]);
            h_chisoPU_endcap_MC_nume[_DY_Full]->Add(h_chisoPU_endcap_MC_nume[pr]);
            h_chisoPU_barrel_MC_deno[_DY_Full]->Add(h_chisoPU_barrel_MC_deno[pr]);
            h_chisoPU_endcap_MC_deno[_DY_Full]->Add(h_chisoPU_endcap_MC_deno[pr]);
            h_chisoPU_barrel_MC_ctrl[_DY_Full]->Add(h_chisoPU_barrel_MC_ctrl[pr]);
            h_chisoPU_endcap_MC_ctrl[_DY_Full]->Add(h_chisoPU_endcap_MC_ctrl[pr]);
            h_MET_MC[_DY_Full]->Add(h_MET_MC[pr]);
            h_MT_barrel_MC_nume[_DY_Full]->Add(h_MT_barrel_MC_nume[pr]);
            h_MT_endcap_MC_nume[_DY_Full]->Add(h_MT_endcap_MC_nume[pr]);
            h_MT_barrel_MC_deno[_DY_Full]->Add(h_MT_barrel_MC_deno[pr]);
            h_MT_endcap_MC_deno[_DY_Full]->Add(h_MT_endcap_MC_deno[pr]);
            h_MT_barrel_MC_ctrl[_DY_Full]->Add(h_MT_barrel_MC_ctrl[pr]);
            h_MT_endcap_MC_ctrl[_DY_Full]->Add(h_MT_endcap_MC_ctrl[pr]);
            h_eta_MC[_DY_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_DY_Full]->Add(h_nVTX_MC[pr]);
            h_mass_test_MC[_DY_Full]->Add(h_mass_test_MC[pr]);
        }

        Color_t color = kOrange - 5;
        h_PFiso_dBeta_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_chiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_chiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_chiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_phiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_phiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_phiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_MET_MC[pr]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_mass_test_MC[pr]->SetFillColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_chiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_chiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_chiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_phiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_phiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_phiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_MET_MC[pr]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_mass_test_MC[pr]->SetLineColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_nume[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_nume[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_deno[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_deno[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_chiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_chiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_chiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_chiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_phiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_phiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_phiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_phiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_MET_MC[pr]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);
        h_mass_test_MC[pr]->SetDirectory(0);

        s_PFiso_dBeta_barrel_nume->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
        s_PFiso_dBeta_endcap_nume->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
        s_PFiso_dBeta_barrel_deno->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
        s_PFiso_dBeta_endcap_deno->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
        s_PFiso_dBeta_barrel_ctrl->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        s_PFiso_dBeta_endcap_ctrl->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        s_PFiso_Rho_barrel_nume->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
        s_PFiso_Rho_endcap_nume->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
        s_PFiso_Rho_barrel_deno->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
        s_PFiso_Rho_endcap_deno->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
        s_PFiso_Rho_barrel_ctrl->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        s_PFiso_Rho_endcap_ctrl->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
        s_SigmaIEtaIEta_barrel_nume->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        s_SigmaIEtaIEta_endcap_nume->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        s_SigmaIEtaIEta_barrel_deno->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        s_SigmaIEtaIEta_endcap_deno->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        s_SigmaIEtaIEta_barrel_ctrl->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        s_SigmaIEtaIEta_endcap_ctrl->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        s_dEtaInSeed_barrel_nume->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
        s_dEtaInSeed_endcap_nume->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
        s_dEtaInSeed_barrel_deno->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
        s_dEtaInSeed_endcap_deno->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
        s_dEtaInSeed_barrel_ctrl->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        s_dEtaInSeed_endcap_ctrl->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        s_dPhiIn_barrel_nume->Add(h_dPhiIn_barrel_MC_nume[pr]);
        s_dPhiIn_endcap_nume->Add(h_dPhiIn_endcap_MC_nume[pr]);
        s_dPhiIn_barrel_deno->Add(h_dPhiIn_barrel_MC_deno[pr]);
        s_dPhiIn_endcap_deno->Add(h_dPhiIn_endcap_MC_deno[pr]);
        s_dPhiIn_barrel_ctrl->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
        s_dPhiIn_endcap_ctrl->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
        s_HoverE_barrel_nume->Add(h_HoverE_barrel_MC_nume[pr]);
        s_HoverE_endcap_nume->Add(h_HoverE_endcap_MC_nume[pr]);
        s_HoverE_barrel_deno->Add(h_HoverE_barrel_MC_deno[pr]);
        s_HoverE_endcap_deno->Add(h_HoverE_endcap_MC_deno[pr]);
        s_HoverE_barrel_ctrl->Add(h_HoverE_barrel_MC_ctrl[pr]);
        s_HoverE_endcap_ctrl->Add(h_HoverE_endcap_MC_ctrl[pr]);
        s_InvEminusInvP_barrel_nume->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
        s_InvEminusInvP_endcap_nume->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
        s_InvEminusInvP_barrel_deno->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
        s_InvEminusInvP_endcap_deno->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
        s_InvEminusInvP_barrel_ctrl->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        s_InvEminusInvP_endcap_ctrl->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        s_chiso_barrel_nume->Add(h_chiso_barrel_MC_nume[pr]);
        s_chiso_endcap_nume->Add(h_chiso_endcap_MC_nume[pr]);
        s_chiso_barrel_deno->Add(h_chiso_barrel_MC_deno[pr]);
        s_chiso_endcap_deno->Add(h_chiso_endcap_MC_deno[pr]);
        s_chiso_barrel_ctrl->Add(h_chiso_barrel_MC_ctrl[pr]);
        s_chiso_endcap_ctrl->Add(h_chiso_endcap_MC_ctrl[pr]);
        s_nhiso_barrel_nume->Add(h_nhiso_barrel_MC_nume[pr]);
        s_nhiso_endcap_nume->Add(h_nhiso_endcap_MC_nume[pr]);
        s_nhiso_barrel_deno->Add(h_nhiso_barrel_MC_deno[pr]);
        s_nhiso_endcap_deno->Add(h_nhiso_endcap_MC_deno[pr]);
        s_nhiso_barrel_ctrl->Add(h_nhiso_barrel_MC_ctrl[pr]);
        s_nhiso_endcap_ctrl->Add(h_nhiso_endcap_MC_ctrl[pr]);
        s_phiso_barrel_nume->Add(h_phiso_barrel_MC_nume[pr]);
        s_phiso_endcap_nume->Add(h_phiso_endcap_MC_nume[pr]);
        s_phiso_barrel_deno->Add(h_phiso_barrel_MC_deno[pr]);
        s_phiso_endcap_deno->Add(h_phiso_endcap_MC_deno[pr]);
        s_phiso_barrel_ctrl->Add(h_phiso_barrel_MC_ctrl[pr]);
        s_phiso_endcap_ctrl->Add(h_phiso_endcap_MC_ctrl[pr]);
        s_chisoPU_barrel_nume->Add(h_chisoPU_barrel_MC_nume[pr]);
        s_chisoPU_endcap_nume->Add(h_chisoPU_endcap_MC_nume[pr]);
        s_chisoPU_barrel_deno->Add(h_chisoPU_barrel_MC_deno[pr]);
        s_chisoPU_endcap_deno->Add(h_chisoPU_endcap_MC_deno[pr]);
        s_chisoPU_barrel_ctrl->Add(h_chisoPU_barrel_MC_ctrl[pr]);
        s_chisoPU_endcap_ctrl->Add(h_chisoPU_endcap_MC_ctrl[pr]);
        s_MET->Add(h_MET_MC[pr]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);
        s_mass_test->Add(h_mass_test_MC[pr]);

        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_PFiso_dBeta_barrel_nume", h_PFiso_dBeta_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_nume", h_PFiso_dBeta_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_deno", h_PFiso_dBeta_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_deno", h_PFiso_dBeta_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_barrel_nume", h_PFiso_Rho_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_endcap_nume", h_PFiso_Rho_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_barrel_deno", h_PFiso_Rho_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_endcap_deno", h_PFiso_Rho_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_barrel_ctrl", h_PFiso_Rho_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_endcap_ctrl", h_PFiso_Rho_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_barrel_nume", h_dEtaInSeed_barrel_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_endcap_nume", h_dEtaInSeed_endcap_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_barrel_deno", h_dEtaInSeed_barrel_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_endcap_deno", h_dEtaInSeed_endcap_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_barrel_ctrl", h_dEtaInSeed_barrel_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_endcap_ctrl", h_dEtaInSeed_endcap_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_barrel_nume", h_dPhiIn_barrel_MC_nume[pr]);
        file->GetObject("h_dPhiIn_endcap_nume", h_dPhiIn_endcap_MC_nume[pr]);
        file->GetObject("h_dPhiIn_barrel_deno", h_dPhiIn_barrel_MC_deno[pr]);
        file->GetObject("h_dPhiIn_endcap_deno", h_dPhiIn_endcap_MC_deno[pr]);
        file->GetObject("h_dPhiIn_barrel_ctrl", h_dPhiIn_barrel_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_endcap_ctrl", h_dPhiIn_endcap_MC_ctrl[pr]);
        file->GetObject("h_HoverE_barrel_nume", h_HoverE_barrel_MC_nume[pr]);
        file->GetObject("h_HoverE_endcap_nume", h_HoverE_endcap_MC_nume[pr]);
        file->GetObject("h_HoverE_barrel_deno", h_HoverE_barrel_MC_deno[pr]);
        file->GetObject("h_HoverE_endcap_deno", h_HoverE_endcap_MC_deno[pr]);
        file->GetObject("h_HoverE_barrel_ctrl", h_HoverE_barrel_MC_ctrl[pr]);
        file->GetObject("h_HoverE_endcap_ctrl", h_HoverE_endcap_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_barrel_nume", h_InvEminusInvP_barrel_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_endcap_nume", h_InvEminusInvP_endcap_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_barrel_deno", h_InvEminusInvP_barrel_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_endcap_deno", h_InvEminusInvP_endcap_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_barrel_ctrl", h_InvEminusInvP_barrel_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_endcap_ctrl", h_InvEminusInvP_endcap_MC_ctrl[pr]);
        file->GetObject("h_chiso_barrel_nume", h_chiso_barrel_MC_nume[pr]);
        file->GetObject("h_chiso_endcap_nume", h_chiso_endcap_MC_nume[pr]);
        file->GetObject("h_chiso_barrel_deno", h_chiso_barrel_MC_deno[pr]);
        file->GetObject("h_chiso_endcap_deno", h_chiso_endcap_MC_deno[pr]);
        file->GetObject("h_chiso_barrel_ctrl", h_chiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_chiso_endcap_ctrl", h_chiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_nhiso_barrel_nume", h_nhiso_barrel_MC_nume[pr]);
        file->GetObject("h_nhiso_endcap_nume", h_nhiso_endcap_MC_nume[pr]);
        file->GetObject("h_nhiso_barrel_deno", h_nhiso_barrel_MC_deno[pr]);
        file->GetObject("h_nhiso_endcap_deno", h_nhiso_endcap_MC_deno[pr]);
        file->GetObject("h_nhiso_barrel_ctrl", h_nhiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_nhiso_endcap_ctrl", h_nhiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_phiso_barrel_nume", h_phiso_barrel_MC_nume[pr]);
        file->GetObject("h_phiso_endcap_nume", h_phiso_endcap_MC_nume[pr]);
        file->GetObject("h_phiso_barrel_deno", h_phiso_barrel_MC_deno[pr]);
        file->GetObject("h_phiso_endcap_deno", h_phiso_endcap_MC_deno[pr]);
        file->GetObject("h_phiso_barrel_ctrl", h_phiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_phiso_endcap_ctrl", h_phiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_barrel_nume", h_chisoPU_barrel_MC_nume[pr]);
        file->GetObject("h_chisoPU_endcap_nume", h_chisoPU_endcap_MC_nume[pr]);
        file->GetObject("h_chisoPU_barrel_deno", h_chisoPU_barrel_MC_deno[pr]);
        file->GetObject("h_chisoPU_endcap_deno", h_chisoPU_endcap_MC_deno[pr]);
        file->GetObject("h_chisoPU_barrel_ctrl", h_chisoPU_barrel_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_endcap_ctrl", h_chisoPU_endcap_MC_ctrl[pr]);
        file->GetObject("h_MET", h_MET_MC[pr]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        file->GetObject("h_mass_test", h_mass_test_MC[pr]);

        removeNegativeBins(h_PFiso_dBeta_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_nume[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_nume[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_deno[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_deno[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_chiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_chiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_chiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_chiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_phiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_phiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_phiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_phiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_MET_MC[pr]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);
        removeNegativeBins(h_mass_test_MC[pr]);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_PFiso_dBeta_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_nume[pr]->Clone("h_PFiso_dBeta_barrel_nume_QCD")));
            h_PFiso_dBeta_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_nume[pr]->Clone("h_PFiso_dBeta_endcap_nume_QCD")));
            h_PFiso_dBeta_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_deno[pr]->Clone("h_PFiso_dBeta_barrel_deno_QCD")));
            h_PFiso_dBeta_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_deno[pr]->Clone("h_PFiso_dBeta_endcap_deno_QCD")));
            h_PFiso_dBeta_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_ctrl[pr]->Clone("h_PFiso_dBeta_barrel_ctrl_QCD")));
            h_PFiso_dBeta_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_ctrl[pr]->Clone("h_PFiso_dBeta_endcap_ctrl_QCD")));
            h_PFiso_Rho_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_nume[pr]->Clone("h_PFiso_Rho_barrel_nume_QCD")));
            h_PFiso_Rho_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_nume[pr]->Clone("h_PFiso_Rho_endcap_nume_QCD")));
            h_PFiso_Rho_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_deno[pr]->Clone("h_PFiso_Rho_barrel_deno_QCD")));
            h_PFiso_Rho_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_deno[pr]->Clone("h_PFiso_Rho_endcap_deno_QCD")));
            h_PFiso_Rho_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_ctrl[pr]->Clone("h_PFiso_Rho_barrel_ctrl_QCD")));
            h_PFiso_Rho_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_ctrl[pr]->Clone("h_PFiso_Rho_endcap_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_nume_QCD")));
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_nume_QCD")));
            h_pT_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_deno_QCD")));
            h_pT_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_deno_QCD")));
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_ctrl_QCD")));
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_ctrl_QCD")));
            h_SigmaIEtaIEta_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_nume[pr]->Clone("h_SigmaIEtaIEta_barrel_nume_QCD")));
            h_SigmaIEtaIEta_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_nume[pr]->Clone("h_SigmaIEtaIEta_endcap_nume_QCD")));
            h_SigmaIEtaIEta_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_deno[pr]->Clone("h_SigmaIEtaIEta_barrel_deno_QCD")));
            h_SigmaIEtaIEta_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_deno[pr]->Clone("h_SigmaIEtaIEta_endcap_deno_QCD")));
            h_SigmaIEtaIEta_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_barrel_ctrl_QCD")));
            h_SigmaIEtaIEta_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_endcap_ctrl_QCD")));
            h_dEtaInSeed_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_nume[pr]->Clone("h_dEtaInSeed_barrel_nume_QCD")));
            h_dEtaInSeed_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_nume[pr]->Clone("h_dEtaInSeed_endcap_nume_QCD")));
            h_dEtaInSeed_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_deno[pr]->Clone("h_dEtaInSeed_barrel_deno_QCD")));
            h_dEtaInSeed_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_deno[pr]->Clone("h_dEtaInSeed_endcap_deno_QCD")));
            h_dEtaInSeed_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_ctrl[pr]->Clone("h_dEtaInSeed_barrel_ctrl_QCD")));
            h_dEtaInSeed_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_ctrl[pr]->Clone("h_dEtaInSeed_endcap_ctrl_QCD")));
            h_dPhiIn_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_nume[pr]->Clone("h_dPhiIn_barrel_nume_QCD")));
            h_dPhiIn_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_nume[pr]->Clone("h_dPhiIn_endcap_nume_QCD")));
            h_dPhiIn_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_deno[pr]->Clone("h_dPhiIn_barrel_deno_QCD")));
            h_dPhiIn_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_deno[pr]->Clone("h_dPhiIn_endcap_deno_QCD")));
            h_dPhiIn_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_ctrl[pr]->Clone("h_dPhiIn_barrel_ctrl_QCD")));
            h_dPhiIn_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_ctrl[pr]->Clone("h_dPhiIn_endcap_ctrl_QCD")));
            h_HoverE_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_MC_nume[pr]->Clone("h_HoverE_barrel_nume_QCD")));
            h_HoverE_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_MC_nume[pr]->Clone("h_HoverE_endcap_nume_QCD")));
            h_HoverE_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_MC_deno[pr]->Clone("h_HoverE_barrel_deno_QCD")));
            h_HoverE_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_MC_deno[pr]->Clone("h_HoverE_endcap_deno_QCD")));
            h_HoverE_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_MC_ctrl[pr]->Clone("h_HoverE_barrel_ctrl_QCD")));
            h_HoverE_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_MC_ctrl[pr]->Clone("h_HoverE_endcap_ctrl_QCD")));
            h_InvEminusInvP_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_nume[pr]->Clone("h_InvEminusInvP_barrel_nume_QCD")));
            h_InvEminusInvP_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_nume[pr]->Clone("h_InvEminusInvP_endcap_nume_QCD")));
            h_InvEminusInvP_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_deno[pr]->Clone("h_InvEminusInvP_barrel_deno_QCD")));
            h_InvEminusInvP_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_deno[pr]->Clone("h_InvEminusInvP_endcap_deno_QCD")));
            h_InvEminusInvP_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_ctrl[pr]->Clone("h_InvEminusInvP_barrel_ctrl_QCD")));
            h_InvEminusInvP_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_ctrl[pr]->Clone("h_InvEminusInvP_endcap_ctrl_QCD")));
            h_chiso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_barrel_MC_nume[pr]->Clone("h_chiso_barrel_nume_QCD")));
            h_chiso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_endcap_MC_nume[pr]->Clone("h_chiso_endcap_nume_QCD")));
            h_chiso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_barrel_MC_deno[pr]->Clone("h_chiso_barrel_deno_QCD")));
            h_chiso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_endcap_MC_deno[pr]->Clone("h_chiso_endcap_deno_QCD")));
            h_chiso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_barrel_MC_ctrl[pr]->Clone("h_chiso_barrel_ctrl_QCD")));
            h_chiso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_chiso_endcap_MC_ctrl[pr]->Clone("h_chiso_endcap_ctrl_QCD")));
            h_nhiso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_barrel_MC_nume[pr]->Clone("h_nhiso_barrel_nume_QCD")));
            h_nhiso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_endcap_MC_nume[pr]->Clone("h_nhiso_endcap_nume_QCD")));
            h_nhiso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_barrel_MC_deno[pr]->Clone("h_nhiso_barrel_deno_QCD")));
            h_nhiso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_endcap_MC_deno[pr]->Clone("h_nhiso_endcap_deno_QCD")));
            h_nhiso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_barrel_MC_ctrl[pr]->Clone("h_nhiso_barrel_ctrl_QCD")));
            h_nhiso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_nhiso_endcap_MC_ctrl[pr]->Clone("h_nhiso_endcap_ctrl_QCD")));
            h_phiso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_barrel_MC_nume[pr]->Clone("h_phiso_barrel_nume_QCD")));
            h_phiso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_endcap_MC_nume[pr]->Clone("h_phiso_endcap_nume_QCD")));
            h_phiso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_barrel_MC_deno[pr]->Clone("h_phiso_barrel_deno_QCD")));
            h_phiso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_endcap_MC_deno[pr]->Clone("h_phiso_endcap_deno_QCD")));
            h_phiso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_barrel_MC_ctrl[pr]->Clone("h_phiso_barrel_ctrl_QCD")));
            h_phiso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_phiso_endcap_MC_ctrl[pr]->Clone("h_phiso_endcap_ctrl_QCD")));
            h_chisoPU_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_barrel_MC_nume[pr]->Clone("h_chisoPU_barrel_nume_QCD")));
            h_chisoPU_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_endcap_MC_nume[pr]->Clone("h_chisoPU_endcap_nume_QCD")));
            h_chisoPU_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_barrel_MC_deno[pr]->Clone("h_chisoPU_barrel_deno_QCD")));
            h_chisoPU_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_endcap_MC_deno[pr]->Clone("h_chisoPU_endcap_deno_QCD")));
            h_chisoPU_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_barrel_MC_ctrl[pr]->Clone("h_chisoPU_barrel_ctrl_QCD")));
            h_chisoPU_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_chisoPU_endcap_MC_ctrl[pr]->Clone("h_chisoPU_endcap_ctrl_QCD")));
            h_MET_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_MET_MC[pr]->Clone("h_MET_QCD")));
            h_MT_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr]->Clone("h_MT_barrel_nume_QCD")));
            h_MT_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr]->Clone("h_MT_endcap_nume_QCD")));
            h_MT_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr]->Clone("h_MT_barrel_deno_QCD")));
            h_MT_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr]->Clone("h_MT_endcap_deno_QCD")));
            h_MT_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr]->Clone("h_MT_barrel_ctrl_QCD")));
            h_MT_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr]->Clone("h_MT_endcap_ctrl_QCD")));
            h_eta_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_QCD")));
            h_nVTX_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_QCD")));
            h_mass_test_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_test_MC[pr]->Clone("h_mass_test_QCD")));

            h_PFiso_dBeta_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_phiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MET_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_nVTX_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_test_MC[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_dBeta_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
            h_PFiso_dBeta_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
            h_PFiso_dBeta_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
            h_PFiso_dBeta_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
            h_PFiso_dBeta_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
            h_PFiso_dBeta_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
            h_PFiso_Rho_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
            h_PFiso_Rho_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
            h_PFiso_Rho_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
            h_PFiso_Rho_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
            h_PFiso_Rho_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
            h_PFiso_Rho_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_SigmaIEtaIEta_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
            h_SigmaIEtaIEta_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
            h_SigmaIEtaIEta_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
            h_SigmaIEtaIEta_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
            h_dEtaInSeed_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
            h_dEtaInSeed_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
            h_dEtaInSeed_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
            h_dEtaInSeed_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
            h_dEtaInSeed_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
            h_dEtaInSeed_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
            h_dPhiIn_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_dPhiIn_barrel_MC_nume[pr]);
            h_dPhiIn_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_dPhiIn_endcap_MC_nume[pr]);
            h_dPhiIn_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_dPhiIn_barrel_MC_deno[pr]);
            h_dPhiIn_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_dPhiIn_endcap_MC_deno[pr]);
            h_dPhiIn_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
            h_dPhiIn_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
            h_HoverE_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_MC_nume[pr]);
            h_HoverE_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_MC_nume[pr]);
            h_HoverE_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_MC_deno[pr]);
            h_HoverE_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_MC_deno[pr]);
            h_HoverE_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_MC_ctrl[pr]);
            h_HoverE_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_MC_ctrl[pr]);
            h_InvEminusInvP_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
            h_InvEminusInvP_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
            h_InvEminusInvP_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
            h_InvEminusInvP_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
            h_InvEminusInvP_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
            h_InvEminusInvP_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
            h_chiso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_chiso_barrel_MC_nume[pr]);
            h_chiso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_chiso_endcap_MC_nume[pr]);
            h_chiso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_chiso_barrel_MC_deno[pr]);
            h_chiso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_chiso_endcap_MC_deno[pr]);
            h_chiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_chiso_barrel_MC_ctrl[pr]);
            h_chiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_chiso_endcap_MC_ctrl[pr]);
            h_nhiso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_nhiso_barrel_MC_nume[pr]);
            h_nhiso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_nhiso_endcap_MC_nume[pr]);
            h_nhiso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_nhiso_barrel_MC_deno[pr]);
            h_nhiso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_nhiso_endcap_MC_deno[pr]);
            h_nhiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_nhiso_barrel_MC_ctrl[pr]);
            h_nhiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_nhiso_endcap_MC_ctrl[pr]);
            h_phiso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_phiso_barrel_MC_nume[pr]);
            h_phiso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_phiso_endcap_MC_nume[pr]);
            h_phiso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_phiso_barrel_MC_deno[pr]);
            h_phiso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_phiso_endcap_MC_deno[pr]);
            h_phiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_phiso_barrel_MC_ctrl[pr]);
            h_phiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_phiso_endcap_MC_ctrl[pr]);
            h_chisoPU_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_chisoPU_barrel_MC_nume[pr]);
            h_chisoPU_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_chisoPU_endcap_MC_nume[pr]);
            h_chisoPU_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_chisoPU_barrel_MC_deno[pr]);
            h_chisoPU_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_chisoPU_endcap_MC_deno[pr]);
            h_chisoPU_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_chisoPU_barrel_MC_ctrl[pr]);
            h_chisoPU_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_chisoPU_endcap_MC_ctrl[pr]);
            h_MET_MC[_QCDEMEnriched_Full]->Add(h_MET_MC[pr]);
            h_MT_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_MT_barrel_MC_nume[pr]);
            h_MT_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_MT_endcap_MC_nume[pr]);
            h_MT_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_MT_barrel_MC_deno[pr]);
            h_MT_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_MT_endcap_MC_deno[pr]);
            h_MT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_MT_barrel_MC_ctrl[pr]);
            h_MT_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_MT_endcap_MC_ctrl[pr]);
            h_eta_MC[_QCDEMEnriched_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_QCDEMEnriched_Full]->Add(h_nVTX_MC[pr]);
            h_mass_test_MC[_QCDEMEnriched_Full]->Add(h_mass_test_MC[pr]);
        }

        Color_t color = kRed + 3;
        h_PFiso_dBeta_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_chiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_chiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_chiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_phiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_phiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_phiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_MET_MC[pr]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_mass_test_MC[pr]->SetFillColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_chiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_chiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_chiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_phiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_phiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_phiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_MET_MC[pr]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_mass_test_MC[pr]->SetLineColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_nume[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_nume[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_deno[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_deno[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_chiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_chiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_chiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_chiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_phiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_phiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_phiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_phiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_MET_MC[pr]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);
        h_mass_test_MC[pr]->SetDirectory(0);

        s_PFiso_dBeta_barrel_nume->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
        s_PFiso_dBeta_endcap_nume->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
        s_PFiso_dBeta_barrel_deno->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
        s_PFiso_dBeta_endcap_deno->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
        s_PFiso_dBeta_barrel_ctrl->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        s_PFiso_dBeta_endcap_ctrl->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        s_PFiso_Rho_barrel_nume->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
        s_PFiso_Rho_endcap_nume->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
        s_PFiso_Rho_barrel_deno->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
        s_PFiso_Rho_endcap_deno->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
        s_PFiso_Rho_barrel_ctrl->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        s_PFiso_Rho_endcap_ctrl->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
        s_SigmaIEtaIEta_barrel_nume->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        s_SigmaIEtaIEta_endcap_nume->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        s_SigmaIEtaIEta_barrel_deno->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        s_SigmaIEtaIEta_endcap_deno->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        s_SigmaIEtaIEta_barrel_ctrl->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        s_SigmaIEtaIEta_endcap_ctrl->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        s_dEtaInSeed_barrel_nume->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
        s_dEtaInSeed_endcap_nume->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
        s_dEtaInSeed_barrel_deno->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
        s_dEtaInSeed_endcap_deno->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
        s_dEtaInSeed_barrel_ctrl->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        s_dEtaInSeed_endcap_ctrl->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        s_dPhiIn_barrel_nume->Add(h_dPhiIn_barrel_MC_nume[pr]);
        s_dPhiIn_endcap_nume->Add(h_dPhiIn_endcap_MC_nume[pr]);
        s_dPhiIn_barrel_deno->Add(h_dPhiIn_barrel_MC_deno[pr]);
        s_dPhiIn_endcap_deno->Add(h_dPhiIn_endcap_MC_deno[pr]);
        s_dPhiIn_barrel_ctrl->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
        s_dPhiIn_endcap_ctrl->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
        s_HoverE_barrel_nume->Add(h_HoverE_barrel_MC_nume[pr]);
        s_HoverE_endcap_nume->Add(h_HoverE_endcap_MC_nume[pr]);
        s_HoverE_barrel_deno->Add(h_HoverE_barrel_MC_deno[pr]);
        s_HoverE_endcap_deno->Add(h_HoverE_endcap_MC_deno[pr]);
        s_HoverE_barrel_ctrl->Add(h_HoverE_barrel_MC_ctrl[pr]);
        s_HoverE_endcap_ctrl->Add(h_HoverE_endcap_MC_ctrl[pr]);
        s_InvEminusInvP_barrel_nume->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
        s_InvEminusInvP_endcap_nume->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
        s_InvEminusInvP_barrel_deno->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
        s_InvEminusInvP_endcap_deno->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
        s_InvEminusInvP_barrel_ctrl->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        s_InvEminusInvP_endcap_ctrl->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        s_chiso_barrel_nume->Add(h_chiso_barrel_MC_nume[pr]);
        s_chiso_endcap_nume->Add(h_chiso_endcap_MC_nume[pr]);
        s_chiso_barrel_deno->Add(h_chiso_barrel_MC_deno[pr]);
        s_chiso_endcap_deno->Add(h_chiso_endcap_MC_deno[pr]);
        s_chiso_barrel_ctrl->Add(h_chiso_barrel_MC_ctrl[pr]);
        s_chiso_endcap_ctrl->Add(h_chiso_endcap_MC_ctrl[pr]);
        s_nhiso_barrel_nume->Add(h_nhiso_barrel_MC_nume[pr]);
        s_nhiso_endcap_nume->Add(h_nhiso_endcap_MC_nume[pr]);
        s_nhiso_barrel_deno->Add(h_nhiso_barrel_MC_deno[pr]);
        s_nhiso_endcap_deno->Add(h_nhiso_endcap_MC_deno[pr]);
        s_nhiso_barrel_ctrl->Add(h_nhiso_barrel_MC_ctrl[pr]);
        s_nhiso_endcap_ctrl->Add(h_nhiso_endcap_MC_ctrl[pr]);
        s_phiso_barrel_nume->Add(h_phiso_barrel_MC_nume[pr]);
        s_phiso_endcap_nume->Add(h_phiso_endcap_MC_nume[pr]);
        s_phiso_barrel_deno->Add(h_phiso_barrel_MC_deno[pr]);
        s_phiso_endcap_deno->Add(h_phiso_endcap_MC_deno[pr]);
        s_phiso_barrel_ctrl->Add(h_phiso_barrel_MC_ctrl[pr]);
        s_phiso_endcap_ctrl->Add(h_phiso_endcap_MC_ctrl[pr]);
        s_chisoPU_barrel_nume->Add(h_chisoPU_barrel_MC_nume[pr]);
        s_chisoPU_endcap_nume->Add(h_chisoPU_endcap_MC_nume[pr]);
        s_chisoPU_barrel_deno->Add(h_chisoPU_barrel_MC_deno[pr]);
        s_chisoPU_endcap_deno->Add(h_chisoPU_endcap_MC_deno[pr]);
        s_chisoPU_barrel_ctrl->Add(h_chisoPU_barrel_MC_ctrl[pr]);
        s_chisoPU_endcap_ctrl->Add(h_chisoPU_endcap_MC_ctrl[pr]);
        s_MET->Add(h_MET_MC[pr]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);
        s_mass_test->Add(h_mass_test_MC[pr]);

        file->Close();
    }

    // GammaJets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_PFiso_dBeta_barrel_nume", h_PFiso_dBeta_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_nume", h_PFiso_dBeta_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_deno", h_PFiso_dBeta_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_deno", h_PFiso_dBeta_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_barrel_nume", h_PFiso_Rho_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_endcap_nume", h_PFiso_Rho_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_Rho_barrel_deno", h_PFiso_Rho_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_endcap_deno", h_PFiso_Rho_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_Rho_barrel_ctrl", h_PFiso_Rho_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_Rho_endcap_ctrl", h_PFiso_Rho_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_barrel_nume", h_dEtaInSeed_barrel_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_endcap_nume", h_dEtaInSeed_endcap_MC_nume[pr]);
        file->GetObject("h_dEtaInSeed_barrel_deno", h_dEtaInSeed_barrel_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_endcap_deno", h_dEtaInSeed_endcap_MC_deno[pr]);
        file->GetObject("h_dEtaInSeed_barrel_ctrl", h_dEtaInSeed_barrel_MC_ctrl[pr]);
        file->GetObject("h_dEtaInSeed_endcap_ctrl", h_dEtaInSeed_endcap_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_barrel_nume", h_dPhiIn_barrel_MC_nume[pr]);
        file->GetObject("h_dPhiIn_endcap_nume", h_dPhiIn_endcap_MC_nume[pr]);
        file->GetObject("h_dPhiIn_barrel_deno", h_dPhiIn_barrel_MC_deno[pr]);
        file->GetObject("h_dPhiIn_endcap_deno", h_dPhiIn_endcap_MC_deno[pr]);
        file->GetObject("h_dPhiIn_barrel_ctrl", h_dPhiIn_barrel_MC_ctrl[pr]);
        file->GetObject("h_dPhiIn_endcap_ctrl", h_dPhiIn_endcap_MC_ctrl[pr]);
        file->GetObject("h_HoverE_barrel_nume", h_HoverE_barrel_MC_nume[pr]);
        file->GetObject("h_HoverE_endcap_nume", h_HoverE_endcap_MC_nume[pr]);
        file->GetObject("h_HoverE_barrel_deno", h_HoverE_barrel_MC_deno[pr]);
        file->GetObject("h_HoverE_endcap_deno", h_HoverE_endcap_MC_deno[pr]);
        file->GetObject("h_HoverE_barrel_ctrl", h_HoverE_barrel_MC_ctrl[pr]);
        file->GetObject("h_HoverE_endcap_ctrl", h_HoverE_endcap_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_barrel_nume", h_InvEminusInvP_barrel_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_endcap_nume", h_InvEminusInvP_endcap_MC_nume[pr]);
        file->GetObject("h_InvEminusInvP_barrel_deno", h_InvEminusInvP_barrel_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_endcap_deno", h_InvEminusInvP_endcap_MC_deno[pr]);
        file->GetObject("h_InvEminusInvP_barrel_ctrl", h_InvEminusInvP_barrel_MC_ctrl[pr]);
        file->GetObject("h_InvEminusInvP_endcap_ctrl", h_InvEminusInvP_endcap_MC_ctrl[pr]);
        file->GetObject("h_chiso_barrel_nume", h_chiso_barrel_MC_nume[pr]);
        file->GetObject("h_chiso_endcap_nume", h_chiso_endcap_MC_nume[pr]);
        file->GetObject("h_chiso_barrel_deno", h_chiso_barrel_MC_deno[pr]);
        file->GetObject("h_chiso_endcap_deno", h_chiso_endcap_MC_deno[pr]);
        file->GetObject("h_chiso_barrel_ctrl", h_chiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_chiso_endcap_ctrl", h_chiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_nhiso_barrel_nume", h_nhiso_barrel_MC_nume[pr]);
        file->GetObject("h_nhiso_endcap_nume", h_nhiso_endcap_MC_nume[pr]);
        file->GetObject("h_nhiso_barrel_deno", h_nhiso_barrel_MC_deno[pr]);
        file->GetObject("h_nhiso_endcap_deno", h_nhiso_endcap_MC_deno[pr]);
        file->GetObject("h_nhiso_barrel_ctrl", h_nhiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_nhiso_endcap_ctrl", h_nhiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_phiso_barrel_nume", h_phiso_barrel_MC_nume[pr]);
        file->GetObject("h_phiso_endcap_nume", h_phiso_endcap_MC_nume[pr]);
        file->GetObject("h_phiso_barrel_deno", h_phiso_barrel_MC_deno[pr]);
        file->GetObject("h_phiso_endcap_deno", h_phiso_endcap_MC_deno[pr]);
        file->GetObject("h_phiso_barrel_ctrl", h_phiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_phiso_endcap_ctrl", h_phiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_barrel_nume", h_chisoPU_barrel_MC_nume[pr]);
        file->GetObject("h_chisoPU_endcap_nume", h_chisoPU_endcap_MC_nume[pr]);
        file->GetObject("h_chisoPU_barrel_deno", h_chisoPU_barrel_MC_deno[pr]);
        file->GetObject("h_chisoPU_endcap_deno", h_chisoPU_endcap_MC_deno[pr]);
        file->GetObject("h_chisoPU_barrel_ctrl", h_chisoPU_barrel_MC_ctrl[pr]);
        file->GetObject("h_chisoPU_endcap_ctrl", h_chisoPU_endcap_MC_ctrl[pr]);
        file->GetObject("h_MET", h_MET_MC[pr]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        file->GetObject("h_mass_test", h_mass_test_MC[pr]);

        removeNegativeBins(h_PFiso_dBeta_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_nume[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_deno[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_nume[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_deno[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_nume[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_nume[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_deno[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_deno[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_nume[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_deno[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_chiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_chiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_chiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_chiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_nhiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_nhiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_phiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_phiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_phiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_phiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_phiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_nume[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_deno[pr]);
        removeNegativeBins(h_chisoPU_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_chisoPU_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_MET_MC[pr]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);
        removeNegativeBins(h_mass_test_MC[pr]);

        if (pr == _GJets_20to100)
        {
            h_PFiso_dBeta_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_nume[pr]->Clone("h_PFiso_dBeta_barrel_nume_GJets")));
            h_PFiso_dBeta_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_nume[pr]->Clone("h_PFiso_dBeta_endcap_nume_GJets")));
            h_PFiso_dBeta_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_deno[pr]->Clone("h_PFiso_dBeta_barrel_deno_GJets")));
            h_PFiso_dBeta_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_deno[pr]->Clone("h_PFiso_dBeta_endcap_deno_GJets")));
            h_PFiso_dBeta_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_barrel_MC_ctrl[pr]->Clone("h_PFiso_dBeta_barrel_ctrl_GJets")));
            h_PFiso_dBeta_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_PFiso_dBeta_endcap_MC_ctrl[pr]->Clone("h_PFiso_dBeta_endcap_ctrl_GJets")));
            h_PFiso_Rho_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_nume[pr]->Clone("h_PFiso_Rho_barrel_nume_GJets")));
            h_PFiso_Rho_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_nume[pr]->Clone("h_PFiso_Rho_endcap_nume_GJets")));
            h_PFiso_Rho_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_deno[pr]->Clone("h_PFiso_Rho_barrel_deno_GJets")));
            h_PFiso_Rho_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_deno[pr]->Clone("h_PFiso_Rho_endcap_deno_GJets")));
            h_PFiso_Rho_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_ctrl[pr]->Clone("h_PFiso_Rho_barrel_ctrl_GJets")));
            h_PFiso_Rho_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_ctrl[pr]->Clone("h_PFiso_Rho_endcap_ctrl_GJets")));
            h_pT_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_nume_GJets")));
            h_pT_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_nume_GJets")));
            h_pT_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_deno_GJets")));
            h_pT_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_deno_GJets")));
            h_pT_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_ctrl_GJets")));
            h_pT_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_ctrl_GJets")));
            h_SigmaIEtaIEta_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_nume[pr]->Clone("h_SigmaIEtaIEta_barrel_nume_GJets")));
            h_SigmaIEtaIEta_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_nume[pr]->Clone("h_SigmaIEtaIEta_endcap_nume_GJets")));
            h_SigmaIEtaIEta_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_deno[pr]->Clone("h_SigmaIEtaIEta_barrel_deno_GJets")));
            h_SigmaIEtaIEta_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_deno[pr]->Clone("h_SigmaIEtaIEta_endcap_deno_GJets")));
            h_SigmaIEtaIEta_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_barrel_ctrl_GJets")));
            h_SigmaIEtaIEta_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->Clone("h_SigmaIEtaIEta_endcap_ctrl_GJets")));
            h_dEtaInSeed_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_nume[pr]->Clone("h_dEtaInSeed_barrel_nume_GJets")));
            h_dEtaInSeed_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_nume[pr]->Clone("h_dEtaInSeed_endcap_nume_GJets")));
            h_dEtaInSeed_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_deno[pr]->Clone("h_dEtaInSeed_barrel_deno_GJets")));
            h_dEtaInSeed_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_deno[pr]->Clone("h_dEtaInSeed_endcap_deno_GJets")));
            h_dEtaInSeed_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_ctrl[pr]->Clone("h_dEtaInSeed_barrel_ctrl_GJets")));
            h_dEtaInSeed_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_ctrl[pr]->Clone("h_dEtaInSeed_endcap_ctrl_GJets")));
            h_dPhiIn_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_nume[pr]->Clone("h_dPhiIn_barrel_nume_GJets")));
            h_dPhiIn_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_nume[pr]->Clone("h_dPhiIn_endcap_nume_GJets")));
            h_dPhiIn_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_deno[pr]->Clone("h_dPhiIn_barrel_deno_GJets")));
            h_dPhiIn_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_deno[pr]->Clone("h_dPhiIn_endcap_deno_GJets")));
            h_dPhiIn_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_ctrl[pr]->Clone("h_dPhiIn_barrel_ctrl_GJets")));
            h_dPhiIn_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_ctrl[pr]->Clone("h_dPhiIn_endcap_ctrl_GJets")));
            h_HoverE_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_nume[pr]->Clone("h_HoverE_barrel_nume_GJets")));
            h_HoverE_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_nume[pr]->Clone("h_HoverE_endcap_nume_GJets")));
            h_HoverE_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_deno[pr]->Clone("h_HoverE_barrel_deno_GJets")));
            h_HoverE_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_deno[pr]->Clone("h_HoverE_endcap_deno_GJets")));
            h_HoverE_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_ctrl[pr]->Clone("h_HoverE_barrel_ctrl_GJets")));
            h_HoverE_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_ctrl[pr]->Clone("h_HoverE_endcap_ctrl_GJets")));
            h_InvEminusInvP_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_nume[pr]->Clone("h_InvEminusInvP_barrel_nume_GJets")));
            h_InvEminusInvP_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_nume[pr]->Clone("h_InvEminusInvP_endcap_nume_GJets")));
            h_InvEminusInvP_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_deno[pr]->Clone("h_InvEminusInvP_barrel_deno_GJets")));
            h_InvEminusInvP_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_deno[pr]->Clone("h_InvEminusInvP_endcap_deno_GJets")));
            h_InvEminusInvP_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_ctrl[pr]->Clone("h_InvEminusInvP_barrel_ctrl_GJets")));
            h_InvEminusInvP_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_ctrl[pr]->Clone("h_InvEminusInvP_endcap_ctrl_GJets")));
            h_chiso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_chiso_barrel_MC_nume[pr]->Clone("h_chiso_barrel_nume_GJets")));
            h_chiso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_chiso_endcap_MC_nume[pr]->Clone("h_chiso_endcap_nume_GJets")));
            h_chiso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_chiso_barrel_MC_deno[pr]->Clone("h_chiso_barrel_deno_GJets")));
            h_chiso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_chiso_endcap_MC_deno[pr]->Clone("h_chiso_endcap_deno_GJets")));
            h_chiso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_chiso_barrel_MC_ctrl[pr]->Clone("h_chiso_barrel_ctrl_GJets")));
            h_chiso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_chiso_endcap_MC_ctrl[pr]->Clone("h_chiso_endcap_ctrl_GJets")));
            h_nhiso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_nume[pr]->Clone("h_nhiso_barrel_nume_GJets")));
            h_nhiso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_nume[pr]->Clone("h_nhiso_endcap_nume_GJets")));
            h_nhiso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_deno[pr]->Clone("h_nhiso_barrel_deno_GJets")));
            h_nhiso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_deno[pr]->Clone("h_nhiso_endcap_deno_GJets")));
            h_nhiso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_nhiso_barrel_MC_ctrl[pr]->Clone("h_nhiso_barrel_ctrl_GJets")));
            h_nhiso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_nhiso_endcap_MC_ctrl[pr]->Clone("h_nhiso_endcap_ctrl_GJets")));
            h_phiso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_phiso_barrel_MC_nume[pr]->Clone("h_phiso_barrel_nume_GJets")));
            h_phiso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_phiso_endcap_MC_nume[pr]->Clone("h_phiso_endcap_nume_GJets")));
            h_phiso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_phiso_barrel_MC_deno[pr]->Clone("h_phiso_barrel_deno_GJets")));
            h_phiso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_phiso_endcap_MC_deno[pr]->Clone("h_phiso_endcap_deno_GJets")));
            h_phiso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_phiso_barrel_MC_ctrl[pr]->Clone("h_phiso_barrel_ctrl_GJets")));
            h_phiso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_phiso_endcap_MC_ctrl[pr]->Clone("h_phiso_endcap_ctrl_GJets")));
            h_chisoPU_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_nume[pr]->Clone("h_chisoPU_barrel_nume_GJets")));
            h_chisoPU_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_nume[pr]->Clone("h_chisoPU_endcap_nume_GJets")));
            h_chisoPU_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_deno[pr]->Clone("h_chisoPU_barrel_deno_GJets")));
            h_chisoPU_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_deno[pr]->Clone("h_chisoPU_endcap_deno_GJets")));
            h_chisoPU_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_chisoPU_barrel_MC_ctrl[pr]->Clone("h_chisoPU_barrel_ctrl_GJets")));
            h_chisoPU_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_chisoPU_endcap_MC_ctrl[pr]->Clone("h_chisoPU_endcap_ctrl_GJets")));
            h_MET_MC[_GJets_Full] = ((TH1D*)(h_MET_MC[pr]->Clone("h_MET_GJets")));
            h_MT_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr]->Clone("h_MT_barrel_nume_GJets")));
            h_MT_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr]->Clone("h_MT_endcap_nume_GJets")));
            h_MT_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr]->Clone("h_MT_barrel_deno_GJets")));
            h_MT_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr]->Clone("h_MT_endcap_deno_GJets")));
            h_MT_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr]->Clone("h_MT_barrel_ctrl_GJets")));
            h_MT_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr]->Clone("h_MT_endcap_ctrl_GJets")));
            h_eta_MC[_GJets_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_GJets")));
            h_nVTX_MC[_GJets_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_GJets")));
            h_mass_test_MC[_GJets_Full] = ((TH1D*)(h_mass_test_MC[pr]->Clone("h_mass_test_GJets")));

            h_PFiso_dBeta_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_PFiso_dBeta_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_chiso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_chiso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_nhiso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_nhiso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_phiso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_phiso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_chisoPU_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_chisoPU_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_MET_MC[_GJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_eta_MC[_GJets_Full]->SetDirectory(0);
            h_nVTX_MC[_GJets_Full]->SetDirectory(0);
            h_mass_test_MC[_GJets_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_dBeta_barrel_MC_nume[_GJets_Full]->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
            h_PFiso_dBeta_endcap_MC_nume[_GJets_Full]->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
            h_PFiso_dBeta_barrel_MC_deno[_GJets_Full]->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
            h_PFiso_dBeta_endcap_MC_deno[_GJets_Full]->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
            h_PFiso_dBeta_barrel_MC_ctrl[_GJets_Full]->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
            h_PFiso_dBeta_endcap_MC_ctrl[_GJets_Full]->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
            h_PFiso_Rho_barrel_MC_nume[_GJets_Full]->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
            h_PFiso_Rho_endcap_MC_nume[_GJets_Full]->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
            h_PFiso_Rho_barrel_MC_deno[_GJets_Full]->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
            h_PFiso_Rho_endcap_MC_deno[_GJets_Full]->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
            h_PFiso_Rho_barrel_MC_ctrl[_GJets_Full]->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
            h_PFiso_Rho_endcap_MC_ctrl[_GJets_Full]->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_GJets_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_GJets_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_MC_deno[_GJets_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_GJets_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_GJets_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_GJets_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_SigmaIEtaIEta_barrel_MC_nume[_GJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
            h_SigmaIEtaIEta_endcap_MC_nume[_GJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
            h_SigmaIEtaIEta_barrel_MC_deno[_GJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
            h_SigmaIEtaIEta_endcap_MC_deno[_GJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
            h_SigmaIEtaIEta_barrel_MC_ctrl[_GJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
            h_SigmaIEtaIEta_endcap_MC_ctrl[_GJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
            h_dEtaInSeed_barrel_MC_nume[_GJets_Full]->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
            h_dEtaInSeed_endcap_MC_nume[_GJets_Full]->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
            h_dEtaInSeed_barrel_MC_deno[_GJets_Full]->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
            h_dEtaInSeed_endcap_MC_deno[_GJets_Full]->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
            h_dEtaInSeed_barrel_MC_ctrl[_GJets_Full]->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
            h_dEtaInSeed_endcap_MC_ctrl[_GJets_Full]->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
            h_dPhiIn_barrel_MC_nume[_GJets_Full]->Add(h_dPhiIn_barrel_MC_nume[pr]);
            h_dPhiIn_endcap_MC_nume[_GJets_Full]->Add(h_dPhiIn_endcap_MC_nume[pr]);
            h_dPhiIn_barrel_MC_deno[_GJets_Full]->Add(h_dPhiIn_barrel_MC_deno[pr]);
            h_dPhiIn_endcap_MC_deno[_GJets_Full]->Add(h_dPhiIn_endcap_MC_deno[pr]);
            h_dPhiIn_barrel_MC_ctrl[_GJets_Full]->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
            h_dPhiIn_endcap_MC_ctrl[_GJets_Full]->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
            h_HoverE_barrel_MC_nume[_GJets_Full]->Add(h_HoverE_barrel_MC_nume[pr]);
            h_HoverE_endcap_MC_nume[_GJets_Full]->Add(h_HoverE_endcap_MC_nume[pr]);
            h_HoverE_barrel_MC_deno[_GJets_Full]->Add(h_HoverE_barrel_MC_deno[pr]);
            h_HoverE_endcap_MC_deno[_GJets_Full]->Add(h_HoverE_endcap_MC_deno[pr]);
            h_HoverE_barrel_MC_ctrl[_GJets_Full]->Add(h_HoverE_barrel_MC_ctrl[pr]);
            h_HoverE_endcap_MC_ctrl[_GJets_Full]->Add(h_HoverE_endcap_MC_ctrl[pr]);
            h_InvEminusInvP_barrel_MC_nume[_GJets_Full]->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
            h_InvEminusInvP_endcap_MC_nume[_GJets_Full]->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
            h_InvEminusInvP_barrel_MC_deno[_GJets_Full]->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
            h_InvEminusInvP_endcap_MC_deno[_GJets_Full]->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
            h_InvEminusInvP_barrel_MC_ctrl[_GJets_Full]->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
            h_InvEminusInvP_endcap_MC_ctrl[_GJets_Full]->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
            h_chiso_barrel_MC_nume[_GJets_Full]->Add(h_chiso_barrel_MC_nume[pr]);
            h_chiso_endcap_MC_nume[_GJets_Full]->Add(h_chiso_endcap_MC_nume[pr]);
            h_chiso_barrel_MC_deno[_GJets_Full]->Add(h_chiso_barrel_MC_deno[pr]);
            h_chiso_endcap_MC_deno[_GJets_Full]->Add(h_chiso_endcap_MC_deno[pr]);
            h_chiso_barrel_MC_ctrl[_GJets_Full]->Add(h_chiso_barrel_MC_ctrl[pr]);
            h_chiso_endcap_MC_ctrl[_GJets_Full]->Add(h_chiso_endcap_MC_ctrl[pr]);
            h_nhiso_barrel_MC_nume[_GJets_Full]->Add(h_nhiso_barrel_MC_nume[pr]);
            h_nhiso_endcap_MC_nume[_GJets_Full]->Add(h_nhiso_endcap_MC_nume[pr]);
            h_nhiso_barrel_MC_deno[_GJets_Full]->Add(h_nhiso_barrel_MC_deno[pr]);
            h_nhiso_endcap_MC_deno[_GJets_Full]->Add(h_nhiso_endcap_MC_deno[pr]);
            h_nhiso_barrel_MC_ctrl[_GJets_Full]->Add(h_nhiso_barrel_MC_ctrl[pr]);
            h_nhiso_endcap_MC_ctrl[_GJets_Full]->Add(h_nhiso_endcap_MC_ctrl[pr]);
            h_phiso_barrel_MC_nume[_GJets_Full]->Add(h_phiso_barrel_MC_nume[pr]);
            h_phiso_endcap_MC_nume[_GJets_Full]->Add(h_phiso_endcap_MC_nume[pr]);
            h_phiso_barrel_MC_deno[_GJets_Full]->Add(h_phiso_barrel_MC_deno[pr]);
            h_phiso_endcap_MC_deno[_GJets_Full]->Add(h_phiso_endcap_MC_deno[pr]);
            h_phiso_barrel_MC_ctrl[_GJets_Full]->Add(h_phiso_barrel_MC_ctrl[pr]);
            h_phiso_endcap_MC_ctrl[_GJets_Full]->Add(h_phiso_endcap_MC_ctrl[pr]);
            h_chisoPU_barrel_MC_nume[_GJets_Full]->Add(h_chisoPU_barrel_MC_nume[pr]);
            h_chisoPU_endcap_MC_nume[_GJets_Full]->Add(h_chisoPU_endcap_MC_nume[pr]);
            h_chisoPU_barrel_MC_deno[_GJets_Full]->Add(h_chisoPU_barrel_MC_deno[pr]);
            h_chisoPU_endcap_MC_deno[_GJets_Full]->Add(h_chisoPU_endcap_MC_deno[pr]);
            h_chisoPU_barrel_MC_ctrl[_GJets_Full]->Add(h_chisoPU_barrel_MC_ctrl[pr]);
            h_chisoPU_endcap_MC_ctrl[_GJets_Full]->Add(h_chisoPU_endcap_MC_ctrl[pr]);
            h_MET_MC[_GJets_Full]->Add(h_MET_MC[pr]);
            h_MT_barrel_MC_nume[_GJets_Full]->Add(h_MT_barrel_MC_nume[pr]);
            h_MT_endcap_MC_nume[_GJets_Full]->Add(h_MT_endcap_MC_nume[pr]);
            h_MT_barrel_MC_deno[_GJets_Full]->Add(h_MT_barrel_MC_deno[pr]);
            h_MT_endcap_MC_deno[_GJets_Full]->Add(h_MT_endcap_MC_deno[pr]);
            h_MT_barrel_MC_ctrl[_GJets_Full]->Add(h_MT_barrel_MC_ctrl[pr]);
            h_MT_endcap_MC_ctrl[_GJets_Full]->Add(h_MT_endcap_MC_ctrl[pr]);
            h_eta_MC[_GJets_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_GJets_Full]->Add(h_nVTX_MC[pr]);
            h_mass_test_MC[_GJets_Full]->Add(h_mass_test_MC[pr]);
        }

        Color_t color = kYellow + 3;
        h_PFiso_dBeta_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_chiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_chiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_chiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_phiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_phiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_phiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetFillColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_MET_MC[pr]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_mass_test_MC[pr]->SetFillColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_nume[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_deno[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_nume[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_nume[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_deno[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_deno[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_chiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_chiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_chiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_chiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_nhiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_nhiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_phiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_phiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_phiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_phiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_phiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_nume[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_deno[pr]->SetLineColor(color);
        h_chisoPU_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_chisoPU_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_MET_MC[pr]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_mass_test_MC[pr]->SetLineColor(color);

        h_PFiso_dBeta_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_dBeta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_dBeta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_nume[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_deno[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_nume[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_deno[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_nume[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_deno[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_nume[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_nume[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_deno[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_deno[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_nume[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_deno[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_chiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_chiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_chiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_chiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_nhiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_nhiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_phiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_phiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_phiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_phiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_phiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_nume[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_deno[pr]->SetDirectory(0);
        h_chisoPU_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_chisoPU_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_MET_MC[pr]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);
        h_mass_test_MC[pr]->SetDirectory(0);

        s_PFiso_dBeta_barrel_nume->Add(h_PFiso_dBeta_barrel_MC_nume[pr]);
        s_PFiso_dBeta_endcap_nume->Add(h_PFiso_dBeta_endcap_MC_nume[pr]);
        s_PFiso_dBeta_barrel_deno->Add(h_PFiso_dBeta_barrel_MC_deno[pr]);
        s_PFiso_dBeta_endcap_deno->Add(h_PFiso_dBeta_endcap_MC_deno[pr]);
        s_PFiso_dBeta_barrel_ctrl->Add(h_PFiso_dBeta_barrel_MC_ctrl[pr]);
        s_PFiso_dBeta_endcap_ctrl->Add(h_PFiso_dBeta_endcap_MC_ctrl[pr]);
        s_PFiso_Rho_barrel_nume->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
        s_PFiso_Rho_endcap_nume->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
        s_PFiso_Rho_barrel_deno->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
        s_PFiso_Rho_endcap_deno->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
        s_PFiso_Rho_barrel_ctrl->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
        s_PFiso_Rho_endcap_ctrl->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
        s_SigmaIEtaIEta_barrel_nume->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
        s_SigmaIEtaIEta_endcap_nume->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
        s_SigmaIEtaIEta_barrel_deno->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
        s_SigmaIEtaIEta_endcap_deno->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
        s_SigmaIEtaIEta_barrel_ctrl->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
        s_SigmaIEtaIEta_endcap_ctrl->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
        s_dEtaInSeed_barrel_nume->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
        s_dEtaInSeed_endcap_nume->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
        s_dEtaInSeed_barrel_deno->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
        s_dEtaInSeed_endcap_deno->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
        s_dEtaInSeed_barrel_ctrl->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
        s_dEtaInSeed_endcap_ctrl->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
        s_dPhiIn_barrel_nume->Add(h_dPhiIn_barrel_MC_nume[pr]);
        s_dPhiIn_endcap_nume->Add(h_dPhiIn_endcap_MC_nume[pr]);
        s_dPhiIn_barrel_deno->Add(h_dPhiIn_barrel_MC_deno[pr]);
        s_dPhiIn_endcap_deno->Add(h_dPhiIn_endcap_MC_deno[pr]);
        s_dPhiIn_barrel_ctrl->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
        s_dPhiIn_endcap_ctrl->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
        s_HoverE_barrel_nume->Add(h_HoverE_barrel_MC_nume[pr]);
        s_HoverE_endcap_nume->Add(h_HoverE_endcap_MC_nume[pr]);
        s_HoverE_barrel_deno->Add(h_HoverE_barrel_MC_deno[pr]);
        s_HoverE_endcap_deno->Add(h_HoverE_endcap_MC_deno[pr]);
        s_HoverE_barrel_ctrl->Add(h_HoverE_barrel_MC_ctrl[pr]);
        s_HoverE_endcap_ctrl->Add(h_HoverE_endcap_MC_ctrl[pr]);
        s_InvEminusInvP_barrel_nume->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
        s_InvEminusInvP_endcap_nume->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
        s_InvEminusInvP_barrel_deno->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
        s_InvEminusInvP_endcap_deno->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
        s_InvEminusInvP_barrel_ctrl->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
        s_InvEminusInvP_endcap_ctrl->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
        s_chiso_barrel_nume->Add(h_chiso_barrel_MC_nume[pr]);
        s_chiso_endcap_nume->Add(h_chiso_endcap_MC_nume[pr]);
        s_chiso_barrel_deno->Add(h_chiso_barrel_MC_deno[pr]);
        s_chiso_endcap_deno->Add(h_chiso_endcap_MC_deno[pr]);
        s_chiso_barrel_ctrl->Add(h_chiso_barrel_MC_ctrl[pr]);
        s_chiso_endcap_ctrl->Add(h_chiso_endcap_MC_ctrl[pr]);
        s_nhiso_barrel_nume->Add(h_nhiso_barrel_MC_nume[pr]);
        s_nhiso_endcap_nume->Add(h_nhiso_endcap_MC_nume[pr]);
        s_nhiso_barrel_deno->Add(h_nhiso_barrel_MC_deno[pr]);
        s_nhiso_endcap_deno->Add(h_nhiso_endcap_MC_deno[pr]);
        s_nhiso_barrel_ctrl->Add(h_nhiso_barrel_MC_ctrl[pr]);
        s_nhiso_endcap_ctrl->Add(h_nhiso_endcap_MC_ctrl[pr]);
        s_phiso_barrel_nume->Add(h_phiso_barrel_MC_nume[pr]);
        s_phiso_endcap_nume->Add(h_phiso_endcap_MC_nume[pr]);
        s_phiso_barrel_deno->Add(h_phiso_barrel_MC_deno[pr]);
        s_phiso_endcap_deno->Add(h_phiso_endcap_MC_deno[pr]);
        s_phiso_barrel_ctrl->Add(h_phiso_barrel_MC_ctrl[pr]);
        s_phiso_endcap_ctrl->Add(h_phiso_endcap_MC_ctrl[pr]);
        s_chisoPU_barrel_nume->Add(h_chisoPU_barrel_MC_nume[pr]);
        s_chisoPU_endcap_nume->Add(h_chisoPU_endcap_MC_nume[pr]);
        s_chisoPU_barrel_deno->Add(h_chisoPU_barrel_MC_deno[pr]);
        s_chisoPU_endcap_deno->Add(h_chisoPU_endcap_MC_deno[pr]);
        s_chisoPU_barrel_ctrl->Add(h_chisoPU_barrel_MC_ctrl[pr]);
        s_chisoPU_endcap_ctrl->Add(h_chisoPU_endcap_MC_ctrl[pr]);
        s_MET->Add(h_MET_MC[pr]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);
        s_mass_test->Add(h_mass_test_MC[pr]);

        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");
        TH1D *h_temp[82];
        if (pr == _SinglePhoton_B)
        {
            file->GetObject("h_PFiso_dBeta_barrel_nume", h_PFiso_dBeta_barrel_data_nume);
            file->GetObject("h_PFiso_dBeta_endcap_nume", h_PFiso_dBeta_endcap_data_nume);
            file->GetObject("h_PFiso_dBeta_barrel_deno", h_PFiso_dBeta_barrel_data_deno);
            file->GetObject("h_PFiso_dBeta_endcap_deno", h_PFiso_dBeta_endcap_data_deno);
            file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_PFiso_dBeta_barrel_data_ctrl);
            file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_PFiso_dBeta_endcap_data_ctrl);
            file->GetObject("h_PFiso_Rho_barrel_nume", h_PFiso_Rho_barrel_data_nume);
            file->GetObject("h_PFiso_Rho_endcap_nume", h_PFiso_Rho_endcap_data_nume);
            file->GetObject("h_PFiso_Rho_barrel_deno", h_PFiso_Rho_barrel_data_deno);
            file->GetObject("h_PFiso_Rho_endcap_deno", h_PFiso_Rho_endcap_data_deno);
            file->GetObject("h_PFiso_Rho_barrel_ctrl", h_PFiso_Rho_barrel_data_ctrl);
            file->GetObject("h_PFiso_Rho_endcap_ctrl", h_PFiso_Rho_endcap_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            file->GetObject("h_pT_barrel_deno", h_pT_barrel_data_deno);
            file->GetObject("h_pT_endcap_deno", h_pT_endcap_data_deno);
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_SigmaIEtaIEta_barrel_data_nume);
            file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_SigmaIEtaIEta_endcap_data_nume);
            file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_SigmaIEtaIEta_barrel_data_deno);
            file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_SigmaIEtaIEta_endcap_data_deno);
            file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_SigmaIEtaIEta_barrel_data_ctrl);
            file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_SigmaIEtaIEta_endcap_data_ctrl);
            file->GetObject("h_dEtaInSeed_barrel_nume", h_dEtaInSeed_barrel_data_nume);
            file->GetObject("h_dEtaInSeed_endcap_nume", h_dEtaInSeed_endcap_data_nume);
            file->GetObject("h_dEtaInSeed_barrel_deno", h_dEtaInSeed_barrel_data_deno);
            file->GetObject("h_dEtaInSeed_endcap_deno", h_dEtaInSeed_endcap_data_deno);
            file->GetObject("h_dEtaInSeed_barrel_ctrl", h_dEtaInSeed_barrel_data_ctrl);
            file->GetObject("h_dEtaInSeed_endcap_ctrl", h_dEtaInSeed_endcap_data_ctrl);
            file->GetObject("h_dPhiIn_barrel_nume", h_dPhiIn_barrel_data_nume);
            file->GetObject("h_dPhiIn_endcap_nume", h_dPhiIn_endcap_data_nume);
            file->GetObject("h_dPhiIn_barrel_deno", h_dPhiIn_barrel_data_deno);
            file->GetObject("h_dPhiIn_endcap_deno", h_dPhiIn_endcap_data_deno);
            file->GetObject("h_dPhiIn_barrel_ctrl", h_dPhiIn_barrel_data_ctrl);
            file->GetObject("h_dPhiIn_endcap_ctrl", h_dPhiIn_endcap_data_ctrl);
            file->GetObject("h_HoverE_barrel_nume", h_HoverE_barrel_data_nume);
            file->GetObject("h_HoverE_endcap_nume", h_HoverE_endcap_data_nume);
            file->GetObject("h_HoverE_barrel_deno", h_HoverE_barrel_data_deno);
            file->GetObject("h_HoverE_endcap_deno", h_HoverE_endcap_data_deno);
            file->GetObject("h_HoverE_barrel_ctrl", h_HoverE_barrel_data_ctrl);
            file->GetObject("h_HoverE_endcap_ctrl", h_HoverE_endcap_data_ctrl);
            file->GetObject("h_InvEminusInvP_barrel_nume", h_InvEminusInvP_barrel_data_nume);
            file->GetObject("h_InvEminusInvP_endcap_nume", h_InvEminusInvP_endcap_data_nume);
            file->GetObject("h_InvEminusInvP_barrel_deno", h_InvEminusInvP_barrel_data_deno);
            file->GetObject("h_InvEminusInvP_endcap_deno", h_InvEminusInvP_endcap_data_deno);
            file->GetObject("h_InvEminusInvP_barrel_ctrl", h_InvEminusInvP_barrel_data_ctrl);
            file->GetObject("h_InvEminusInvP_endcap_ctrl", h_InvEminusInvP_endcap_data_ctrl);
            file->GetObject("h_chiso_barrel_nume", h_chiso_barrel_data_nume);
            file->GetObject("h_chiso_endcap_nume", h_chiso_endcap_data_nume);
            file->GetObject("h_chiso_barrel_deno", h_chiso_barrel_data_deno);
            file->GetObject("h_chiso_endcap_deno", h_chiso_endcap_data_deno);
            file->GetObject("h_chiso_barrel_ctrl", h_chiso_barrel_data_ctrl);
            file->GetObject("h_chiso_endcap_ctrl", h_chiso_endcap_data_ctrl);
            file->GetObject("h_nhiso_barrel_nume", h_nhiso_barrel_data_nume);
            file->GetObject("h_nhiso_endcap_nume", h_nhiso_endcap_data_nume);
            file->GetObject("h_nhiso_barrel_deno", h_nhiso_barrel_data_deno);
            file->GetObject("h_nhiso_endcap_deno", h_nhiso_endcap_data_deno);
            file->GetObject("h_nhiso_barrel_ctrl", h_nhiso_barrel_data_ctrl);
            file->GetObject("h_nhiso_endcap_ctrl", h_nhiso_endcap_data_ctrl);
            file->GetObject("h_phiso_barrel_nume", h_phiso_barrel_data_nume);
            file->GetObject("h_phiso_endcap_nume", h_phiso_endcap_data_nume);
            file->GetObject("h_phiso_barrel_deno", h_phiso_barrel_data_deno);
            file->GetObject("h_phiso_endcap_deno", h_phiso_endcap_data_deno);
            file->GetObject("h_phiso_barrel_ctrl", h_phiso_barrel_data_ctrl);
            file->GetObject("h_phiso_endcap_ctrl", h_phiso_endcap_data_ctrl);
            file->GetObject("h_chisoPU_barrel_nume", h_chisoPU_barrel_data_nume);
            file->GetObject("h_chisoPU_endcap_nume", h_chisoPU_endcap_data_nume);
            file->GetObject("h_chisoPU_barrel_deno", h_chisoPU_barrel_data_deno);
            file->GetObject("h_chisoPU_endcap_deno", h_chisoPU_endcap_data_deno);
            file->GetObject("h_chisoPU_barrel_ctrl", h_chisoPU_barrel_data_ctrl);
            file->GetObject("h_chisoPU_endcap_ctrl", h_chisoPU_endcap_data_ctrl);
            file->GetObject("h_MET", h_MET_data);
            file->GetObject("h_MT_barrel_nume", h_MT_barrel_data_nume);
            file->GetObject("h_MT_endcap_nume", h_MT_endcap_data_nume);
            file->GetObject("h_MT_barrel_deno", h_MT_barrel_data_deno);
            file->GetObject("h_MT_endcap_deno", h_MT_endcap_data_deno);
            file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_data_ctrl);
            file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_data_ctrl);
            file->GetObject("h_eta_deno", h_eta_data);
            file->GetObject("h_nVTX", h_nVTX_data);
            file->GetObject("h_mass_test", h_mass_test_data);

            removeNegativeBins(h_PFiso_dBeta_barrel_data_nume);
            removeNegativeBins(h_PFiso_dBeta_endcap_data_nume);
            removeNegativeBins(h_PFiso_dBeta_barrel_data_deno);
            removeNegativeBins(h_PFiso_dBeta_endcap_data_deno);
            removeNegativeBins(h_PFiso_dBeta_barrel_data_ctrl);
            removeNegativeBins(h_PFiso_dBeta_endcap_data_ctrl);
            removeNegativeBins(h_PFiso_Rho_barrel_data_nume);
            removeNegativeBins(h_PFiso_Rho_endcap_data_nume);
            removeNegativeBins(h_PFiso_Rho_barrel_data_deno);
            removeNegativeBins(h_PFiso_Rho_endcap_data_deno);
            removeNegativeBins(h_PFiso_Rho_barrel_data_ctrl);
            removeNegativeBins(h_PFiso_Rho_endcap_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_nume);
            removeNegativeBins(h_pT_barrel_data_deno);
            removeNegativeBins(h_pT_endcap_data_deno);
            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_SigmaIEtaIEta_barrel_data_nume);
            removeNegativeBins(h_SigmaIEtaIEta_endcap_data_nume);
            removeNegativeBins(h_SigmaIEtaIEta_barrel_data_deno);
            removeNegativeBins(h_SigmaIEtaIEta_endcap_data_deno);
            removeNegativeBins(h_SigmaIEtaIEta_barrel_data_ctrl);
            removeNegativeBins(h_SigmaIEtaIEta_endcap_data_ctrl);
            removeNegativeBins(h_dEtaInSeed_barrel_data_nume);
            removeNegativeBins(h_dEtaInSeed_endcap_data_nume);
            removeNegativeBins(h_dEtaInSeed_barrel_data_deno);
            removeNegativeBins(h_dEtaInSeed_endcap_data_deno);
            removeNegativeBins(h_dEtaInSeed_barrel_data_ctrl);
            removeNegativeBins(h_dEtaInSeed_endcap_data_ctrl);
            removeNegativeBins(h_dPhiIn_barrel_data_nume);
            removeNegativeBins(h_dPhiIn_endcap_data_nume);
            removeNegativeBins(h_dPhiIn_barrel_data_deno);
            removeNegativeBins(h_dPhiIn_endcap_data_deno);
            removeNegativeBins(h_dPhiIn_barrel_data_ctrl);
            removeNegativeBins(h_dPhiIn_endcap_data_ctrl);
            removeNegativeBins(h_HoverE_barrel_data_nume);
            removeNegativeBins(h_HoverE_endcap_data_nume);
            removeNegativeBins(h_HoverE_barrel_data_deno);
            removeNegativeBins(h_HoverE_endcap_data_deno);
            removeNegativeBins(h_HoverE_barrel_data_ctrl);
            removeNegativeBins(h_HoverE_endcap_data_ctrl);
            removeNegativeBins(h_InvEminusInvP_barrel_data_nume);
            removeNegativeBins(h_InvEminusInvP_endcap_data_nume);
            removeNegativeBins(h_InvEminusInvP_barrel_data_deno);
            removeNegativeBins(h_InvEminusInvP_endcap_data_deno);
            removeNegativeBins(h_InvEminusInvP_barrel_data_ctrl);
            removeNegativeBins(h_InvEminusInvP_endcap_data_ctrl);
            removeNegativeBins(h_chiso_barrel_data_nume);
            removeNegativeBins(h_chiso_endcap_data_nume);
            removeNegativeBins(h_chiso_barrel_data_deno);
            removeNegativeBins(h_chiso_endcap_data_deno);
            removeNegativeBins(h_chiso_barrel_data_ctrl);
            removeNegativeBins(h_chiso_endcap_data_ctrl);
            removeNegativeBins(h_nhiso_barrel_data_nume);
            removeNegativeBins(h_nhiso_endcap_data_nume);
            removeNegativeBins(h_nhiso_barrel_data_deno);
            removeNegativeBins(h_nhiso_endcap_data_deno);
            removeNegativeBins(h_nhiso_barrel_data_ctrl);
            removeNegativeBins(h_nhiso_endcap_data_ctrl);
            removeNegativeBins(h_phiso_barrel_data_nume);
            removeNegativeBins(h_phiso_endcap_data_nume);
            removeNegativeBins(h_phiso_barrel_data_deno);
            removeNegativeBins(h_phiso_endcap_data_deno);
            removeNegativeBins(h_phiso_barrel_data_ctrl);
            removeNegativeBins(h_phiso_endcap_data_ctrl);
            removeNegativeBins(h_chisoPU_barrel_data_nume);
            removeNegativeBins(h_chisoPU_endcap_data_nume);
            removeNegativeBins(h_chisoPU_barrel_data_deno);
            removeNegativeBins(h_chisoPU_endcap_data_deno);
            removeNegativeBins(h_chisoPU_barrel_data_ctrl);
            removeNegativeBins(h_chisoPU_endcap_data_ctrl);
            removeNegativeBins(h_MET_data);
            removeNegativeBins(h_MT_barrel_data_nume);
            removeNegativeBins(h_MT_endcap_data_nume);
            removeNegativeBins(h_MT_barrel_data_deno);
            removeNegativeBins(h_MT_endcap_data_deno);
            removeNegativeBins(h_MT_barrel_data_ctrl);
            removeNegativeBins(h_MT_endcap_data_ctrl);
            removeNegativeBins(h_eta_data);
            removeNegativeBins(h_nVTX_data);
            removeNegativeBins(h_mass_test_data);
        }
        else
        {
            file->GetObject("h_PFiso_dBeta_barrel_nume", h_temp[0]);
            file->GetObject("h_PFiso_dBeta_endcap_nume", h_temp[1]);
            file->GetObject("h_PFiso_dBeta_barrel_deno", h_temp[2]);
            file->GetObject("h_PFiso_dBeta_endcap_deno", h_temp[3]);
            file->GetObject("h_PFiso_dBeta_barrel_ctrl", h_temp[4]);
            file->GetObject("h_PFiso_dBeta_endcap_ctrl", h_temp[5]);
            file->GetObject("h_PFiso_Rho_barrel_nume", h_temp[6]);
            file->GetObject("h_PFiso_Rho_endcap_nume", h_temp[7]);
            file->GetObject("h_PFiso_Rho_barrel_deno", h_temp[8]);
            file->GetObject("h_PFiso_Rho_endcap_deno", h_temp[9]);
            file->GetObject("h_PFiso_Rho_barrel_ctrl", h_temp[10]);
            file->GetObject("h_PFiso_Rho_endcap_ctrl", h_temp[11]);
            file->GetObject("h_pT_barrel_nume", h_temp[12]);
            file->GetObject("h_pT_endcap_nume", h_temp[13]);
            file->GetObject("h_pT_barrel_deno", h_temp[14]);
            file->GetObject("h_pT_endcap_deno", h_temp[15]);
            file->GetObject("h_pT_barrel_ctrl", h_temp[16]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[17]);
            file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_temp[18]);
            file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_temp[19]);
            file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_temp[20]);
            file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_temp[21]);
            file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_temp[22]);
            file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_temp[23]);
            file->GetObject("h_dEtaInSeed_barrel_nume", h_temp[24]);
            file->GetObject("h_dEtaInSeed_endcap_nume", h_temp[25]);
            file->GetObject("h_dEtaInSeed_barrel_deno", h_temp[26]);
            file->GetObject("h_dEtaInSeed_endcap_deno", h_temp[27]);
            file->GetObject("h_dEtaInSeed_barrel_ctrl", h_temp[28]);
            file->GetObject("h_dEtaInSeed_endcap_ctrl", h_temp[29]);
            file->GetObject("h_dPhiIn_barrel_nume", h_temp[30]);
            file->GetObject("h_dPhiIn_endcap_nume", h_temp[31]);
            file->GetObject("h_dPhiIn_barrel_deno", h_temp[32]);
            file->GetObject("h_dPhiIn_endcap_deno", h_temp[33]);
            file->GetObject("h_dPhiIn_barrel_ctrl", h_temp[34]);
            file->GetObject("h_dPhiIn_endcap_ctrl", h_temp[35]);
            file->GetObject("h_HoverE_barrel_nume", h_temp[36]);
            file->GetObject("h_HoverE_endcap_nume", h_temp[37]);
            file->GetObject("h_HoverE_barrel_deno", h_temp[38]);
            file->GetObject("h_HoverE_endcap_deno", h_temp[39]);
            file->GetObject("h_HoverE_barrel_ctrl", h_temp[40]);
            file->GetObject("h_HoverE_endcap_ctrl", h_temp[41]);
            file->GetObject("h_InvEminusInvP_barrel_nume", h_temp[42]);
            file->GetObject("h_InvEminusInvP_endcap_nume", h_temp[43]);
            file->GetObject("h_InvEminusInvP_barrel_deno", h_temp[44]);
            file->GetObject("h_InvEminusInvP_endcap_deno", h_temp[45]);
            file->GetObject("h_InvEminusInvP_barrel_ctrl", h_temp[46]);
            file->GetObject("h_InvEminusInvP_endcap_ctrl", h_temp[47]);
            file->GetObject("h_chiso_barrel_nume", h_temp[48]);
            file->GetObject("h_chiso_endcap_nume", h_temp[49]);
            file->GetObject("h_chiso_barrel_deno", h_temp[50]);
            file->GetObject("h_chiso_endcap_deno", h_temp[51]);
            file->GetObject("h_chiso_barrel_ctrl", h_temp[52]);
            file->GetObject("h_chiso_endcap_ctrl", h_temp[53]);
            file->GetObject("h_nhiso_barrel_nume", h_temp[54]);
            file->GetObject("h_nhiso_endcap_nume", h_temp[55]);
            file->GetObject("h_nhiso_barrel_deno", h_temp[56]);
            file->GetObject("h_nhiso_endcap_deno", h_temp[57]);
            file->GetObject("h_nhiso_barrel_ctrl", h_temp[58]);
            file->GetObject("h_nhiso_endcap_ctrl", h_temp[59]);
            file->GetObject("h_phiso_barrel_nume", h_temp[60]);
            file->GetObject("h_phiso_endcap_nume", h_temp[61]);
            file->GetObject("h_phiso_barrel_deno", h_temp[62]);
            file->GetObject("h_phiso_endcap_deno", h_temp[63]);
            file->GetObject("h_phiso_barrel_ctrl", h_temp[64]);
            file->GetObject("h_phiso_endcap_ctrl", h_temp[65]);
            file->GetObject("h_chisoPU_barrel_nume", h_temp[66]);
            file->GetObject("h_chisoPU_endcap_nume", h_temp[67]);
            file->GetObject("h_chisoPU_barrel_deno", h_temp[68]);
            file->GetObject("h_chisoPU_endcap_deno", h_temp[69]);
            file->GetObject("h_chisoPU_barrel_ctrl", h_temp[70]);
            file->GetObject("h_chisoPU_endcap_ctrl", h_temp[71]);
            file->GetObject("h_MET", h_temp[72]);
            file->GetObject("h_MT_barrel_nume", h_temp[73]);
            file->GetObject("h_MT_endcap_nume", h_temp[74]);
            file->GetObject("h_MT_barrel_deno", h_temp[75]);
            file->GetObject("h_MT_endcap_deno", h_temp[76]);
            file->GetObject("h_MT_barrel_ctrl", h_temp[77]);
            file->GetObject("h_MT_endcap_ctrl", h_temp[78]);
            file->GetObject("h_eta_deno", h_temp[79]);
            file->GetObject("h_nVTX", h_temp[80]);
            file->GetObject("h_mass_test", h_temp[81]);

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
            removeNegativeBins(h_temp[16]);
            removeNegativeBins(h_temp[17]);
            removeNegativeBins(h_temp[18]);
            removeNegativeBins(h_temp[19]);
            removeNegativeBins(h_temp[20]);
            removeNegativeBins(h_temp[21]);
            removeNegativeBins(h_temp[22]);
            removeNegativeBins(h_temp[23]);
            removeNegativeBins(h_temp[24]);
            removeNegativeBins(h_temp[25]);
            removeNegativeBins(h_temp[26]);
            removeNegativeBins(h_temp[27]);
            removeNegativeBins(h_temp[28]);
            removeNegativeBins(h_temp[29]);
            removeNegativeBins(h_temp[30]);
            removeNegativeBins(h_temp[31]);
            removeNegativeBins(h_temp[32]);
            removeNegativeBins(h_temp[33]);
            removeNegativeBins(h_temp[34]);
            removeNegativeBins(h_temp[35]);
            removeNegativeBins(h_temp[36]);
            removeNegativeBins(h_temp[37]);
            removeNegativeBins(h_temp[38]);
            removeNegativeBins(h_temp[39]);
            removeNegativeBins(h_temp[40]);
            removeNegativeBins(h_temp[41]);
            removeNegativeBins(h_temp[42]);
            removeNegativeBins(h_temp[43]);
            removeNegativeBins(h_temp[44]);
            removeNegativeBins(h_temp[45]);
            removeNegativeBins(h_temp[46]);
            removeNegativeBins(h_temp[47]);
            removeNegativeBins(h_temp[48]);
            removeNegativeBins(h_temp[49]);
            removeNegativeBins(h_temp[50]);
            removeNegativeBins(h_temp[51]);
            removeNegativeBins(h_temp[52]);
            removeNegativeBins(h_temp[53]);
            removeNegativeBins(h_temp[54]);
            removeNegativeBins(h_temp[55]);
            removeNegativeBins(h_temp[56]);
            removeNegativeBins(h_temp[57]);
            removeNegativeBins(h_temp[58]);
            removeNegativeBins(h_temp[59]);
            removeNegativeBins(h_temp[60]);
            removeNegativeBins(h_temp[61]);
            removeNegativeBins(h_temp[62]);
            removeNegativeBins(h_temp[63]);
            removeNegativeBins(h_temp[64]);
            removeNegativeBins(h_temp[65]);
            removeNegativeBins(h_temp[66]);
            removeNegativeBins(h_temp[67]);
            removeNegativeBins(h_temp[68]);
            removeNegativeBins(h_temp[69]);
            removeNegativeBins(h_temp[70]);
            removeNegativeBins(h_temp[71]);
            removeNegativeBins(h_temp[72]);
            removeNegativeBins(h_temp[73]);
            removeNegativeBins(h_temp[74]);
            removeNegativeBins(h_temp[75]);
            removeNegativeBins(h_temp[76]);
            removeNegativeBins(h_temp[77]);
            removeNegativeBins(h_temp[78]);
            removeNegativeBins(h_temp[79]);
            removeNegativeBins(h_temp[80]);
            removeNegativeBins(h_temp[81]);

            h_PFiso_dBeta_barrel_data_nume->Add(h_temp[0]);
            h_PFiso_dBeta_endcap_data_nume->Add(h_temp[1]);
            h_PFiso_dBeta_barrel_data_deno->Add(h_temp[2]);
            h_PFiso_dBeta_endcap_data_deno->Add(h_temp[3]);
            h_PFiso_dBeta_barrel_data_ctrl->Add(h_temp[4]);
            h_PFiso_dBeta_endcap_data_ctrl->Add(h_temp[5]);
            h_PFiso_Rho_barrel_data_nume->Add(h_temp[6]);
            h_PFiso_Rho_endcap_data_nume->Add(h_temp[7]);
            h_PFiso_Rho_barrel_data_deno->Add(h_temp[8]);
            h_PFiso_Rho_endcap_data_deno->Add(h_temp[9]);
            h_PFiso_Rho_barrel_data_ctrl->Add(h_temp[10]);
            h_PFiso_Rho_endcap_data_ctrl->Add(h_temp[11]);
            h_pT_barrel_data_nume->Add(h_temp[12]);
            h_pT_endcap_data_nume->Add(h_temp[13]);
            h_pT_barrel_data_deno->Add(h_temp[14]);
            h_pT_endcap_data_deno->Add(h_temp[15]);
            h_pT_barrel_data_ctrl->Add(h_temp[16]);
            h_pT_endcap_data_ctrl->Add(h_temp[17]);
            h_SigmaIEtaIEta_barrel_data_nume->Add(h_temp[18]);
            h_SigmaIEtaIEta_endcap_data_nume->Add(h_temp[19]);
            h_SigmaIEtaIEta_barrel_data_deno->Add(h_temp[20]);
            h_SigmaIEtaIEta_endcap_data_deno->Add(h_temp[21]);
            h_SigmaIEtaIEta_barrel_data_ctrl->Add(h_temp[22]);
            h_SigmaIEtaIEta_endcap_data_ctrl->Add(h_temp[23]);
            h_dEtaInSeed_barrel_data_nume->Add(h_temp[24]);
            h_dEtaInSeed_endcap_data_nume->Add(h_temp[25]);
            h_dEtaInSeed_barrel_data_deno->Add(h_temp[26]);
            h_dEtaInSeed_endcap_data_deno->Add(h_temp[27]);
            h_dEtaInSeed_barrel_data_ctrl->Add(h_temp[28]);
            h_dEtaInSeed_endcap_data_ctrl->Add(h_temp[29]);
            h_dPhiIn_barrel_data_nume->Add(h_temp[30]);
            h_dPhiIn_endcap_data_nume->Add(h_temp[31]);
            h_dPhiIn_barrel_data_deno->Add(h_temp[32]);
            h_dPhiIn_endcap_data_deno->Add(h_temp[33]);
            h_dPhiIn_barrel_data_ctrl->Add(h_temp[34]);
            h_dPhiIn_endcap_data_ctrl->Add(h_temp[35]);
            h_HoverE_barrel_data_nume->Add(h_temp[36]);
            h_HoverE_endcap_data_nume->Add(h_temp[37]);
            h_HoverE_barrel_data_deno->Add(h_temp[38]);
            h_HoverE_endcap_data_deno->Add(h_temp[39]);
            h_HoverE_barrel_data_ctrl->Add(h_temp[40]);
            h_HoverE_endcap_data_ctrl->Add(h_temp[41]);
            h_InvEminusInvP_barrel_data_nume->Add(h_temp[42]);
            h_InvEminusInvP_endcap_data_nume->Add(h_temp[43]);
            h_InvEminusInvP_barrel_data_deno->Add(h_temp[44]);
            h_InvEminusInvP_endcap_data_deno->Add(h_temp[45]);
            h_InvEminusInvP_barrel_data_ctrl->Add(h_temp[46]);
            h_InvEminusInvP_endcap_data_ctrl->Add(h_temp[47]);
            h_chiso_barrel_data_nume->Add(h_temp[48]);
            h_chiso_endcap_data_nume->Add(h_temp[49]);
            h_chiso_barrel_data_deno->Add(h_temp[50]);
            h_chiso_endcap_data_deno->Add(h_temp[51]);
            h_chiso_barrel_data_ctrl->Add(h_temp[52]);
            h_chiso_endcap_data_ctrl->Add(h_temp[53]);
            h_nhiso_barrel_data_nume->Add(h_temp[54]);
            h_nhiso_endcap_data_nume->Add(h_temp[55]);
            h_nhiso_barrel_data_deno->Add(h_temp[56]);
            h_nhiso_endcap_data_deno->Add(h_temp[57]);
            h_nhiso_barrel_data_ctrl->Add(h_temp[58]);
            h_nhiso_endcap_data_ctrl->Add(h_temp[59]);
            h_phiso_barrel_data_nume->Add(h_temp[60]);
            h_phiso_endcap_data_nume->Add(h_temp[61]);
            h_phiso_barrel_data_deno->Add(h_temp[62]);
            h_phiso_endcap_data_deno->Add(h_temp[63]);
            h_phiso_barrel_data_ctrl->Add(h_temp[64]);
            h_phiso_endcap_data_ctrl->Add(h_temp[65]);
            h_chisoPU_barrel_data_nume->Add(h_temp[66]);
            h_chisoPU_endcap_data_nume->Add(h_temp[67]);
            h_chisoPU_barrel_data_deno->Add(h_temp[68]);
            h_chisoPU_endcap_data_deno->Add(h_temp[69]);
            h_chisoPU_barrel_data_ctrl->Add(h_temp[70]);
            h_chisoPU_endcap_data_ctrl->Add(h_temp[71]);
            h_MET_data->Add(h_temp[72]);
            h_MT_barrel_data_nume->Add(h_temp[73]);
            h_MT_endcap_data_nume->Add(h_temp[74]);
            h_MT_barrel_data_deno->Add(h_temp[75]);
            h_MT_endcap_data_deno->Add(h_temp[76]);
            h_MT_barrel_data_ctrl->Add(h_temp[77]);
            h_MT_endcap_data_ctrl->Add(h_temp[78]);
            h_eta_data->Add(h_temp[79]);
            h_nVTX_data->Add(h_temp[80]);
            h_mass_test_data->Add(h_temp[81]);
        }
    }

    h_PFiso_dBeta_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_dBeta_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_dBeta_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_dBeta_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_dBeta_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_PFiso_dBeta_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_HoverE_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_HoverE_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_HoverE_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_chiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_chiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_chiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_chiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_chiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_chiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_nhiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_nhiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_nhiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_nhiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_nhiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_nhiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_phiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_phiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_phiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_phiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_phiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_phiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_chisoPU_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_MET_data->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_eta_data->SetMarkerStyle(kFullDotLarge);
    h_nVTX_data->SetMarkerStyle(kFullDotLarge);
    h_mass_test_data->SetMarkerStyle(kFullDotLarge);

    h_PFiso_dBeta_barrel_data_nume->SetMarkerColor(kBlack);
    h_PFiso_dBeta_endcap_data_nume->SetMarkerColor(kBlack);
    h_PFiso_dBeta_barrel_data_deno->SetMarkerColor(kBlack);
    h_PFiso_dBeta_endcap_data_deno->SetMarkerColor(kBlack);
    h_PFiso_dBeta_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_PFiso_dBeta_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_PFiso_Rho_barrel_data_nume->SetMarkerColor(kBlack);
    h_PFiso_Rho_endcap_data_nume->SetMarkerColor(kBlack);
    h_PFiso_Rho_barrel_data_deno->SetMarkerColor(kBlack);
    h_PFiso_Rho_endcap_data_deno->SetMarkerColor(kBlack);
    h_PFiso_Rho_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_PFiso_Rho_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_pT_barrel_data_nume->SetMarkerColor(kBlack);
    h_pT_endcap_data_nume->SetMarkerColor(kBlack);
    h_pT_barrel_data_deno->SetMarkerColor(kBlack);
    h_pT_endcap_data_deno->SetMarkerColor(kBlack);
    h_pT_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_pT_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_nume->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_nume->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_deno->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_deno->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_dEtaInSeed_barrel_data_nume->SetMarkerColor(kBlack);
    h_dEtaInSeed_endcap_data_nume->SetMarkerColor(kBlack);
    h_dEtaInSeed_barrel_data_deno->SetMarkerColor(kBlack);
    h_dEtaInSeed_endcap_data_deno->SetMarkerColor(kBlack);
    h_dEtaInSeed_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_dEtaInSeed_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_dPhiIn_barrel_data_nume->SetMarkerColor(kBlack);
    h_dPhiIn_endcap_data_nume->SetMarkerColor(kBlack);
    h_dPhiIn_barrel_data_deno->SetMarkerColor(kBlack);
    h_dPhiIn_endcap_data_deno->SetMarkerColor(kBlack);
    h_dPhiIn_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_dPhiIn_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_HoverE_barrel_data_nume->SetMarkerColor(kBlack);
    h_HoverE_endcap_data_nume->SetMarkerColor(kBlack);
    h_HoverE_barrel_data_deno->SetMarkerColor(kBlack);
    h_HoverE_endcap_data_deno->SetMarkerColor(kBlack);
    h_HoverE_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_HoverE_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_InvEminusInvP_barrel_data_nume->SetMarkerColor(kBlack);
    h_InvEminusInvP_endcap_data_nume->SetMarkerColor(kBlack);
    h_InvEminusInvP_barrel_data_deno->SetMarkerColor(kBlack);
    h_InvEminusInvP_endcap_data_deno->SetMarkerColor(kBlack);
    h_InvEminusInvP_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_InvEminusInvP_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_chiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_chiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_chiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_chiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_chiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_chiso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_nhiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_nhiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_nhiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_nhiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_nhiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_nhiso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_phiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_phiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_phiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_phiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_phiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_phiso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_chisoPU_barrel_data_nume->SetMarkerColor(kBlack);
    h_chisoPU_endcap_data_nume->SetMarkerColor(kBlack);
    h_chisoPU_barrel_data_deno->SetMarkerColor(kBlack);
    h_chisoPU_endcap_data_deno->SetMarkerColor(kBlack);
    h_chisoPU_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_chisoPU_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_MET_data->SetMarkerColor(kBlack);
    h_MT_barrel_data_nume->SetMarkerColor(kBlack);
    h_MT_endcap_data_nume->SetMarkerColor(kBlack);
    h_MT_barrel_data_deno->SetMarkerColor(kBlack);
    h_MT_endcap_data_deno->SetMarkerColor(kBlack);
    h_MT_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_MT_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_eta_data->SetMarkerColor(kBlack);
    h_nVTX_data->SetMarkerColor(kBlack);
    h_mass_test_data->SetMarkerColor(kBlack);

    h_PFiso_dBeta_barrel_data_nume->SetLineColor(kBlack);
    h_PFiso_dBeta_endcap_data_nume->SetLineColor(kBlack);
    h_PFiso_dBeta_barrel_data_deno->SetLineColor(kBlack);
    h_PFiso_dBeta_endcap_data_deno->SetLineColor(kBlack);
    h_PFiso_dBeta_barrel_data_ctrl->SetLineColor(kBlack);
    h_PFiso_dBeta_endcap_data_ctrl->SetLineColor(kBlack);
    h_PFiso_Rho_barrel_data_nume->SetLineColor(kBlack);
    h_PFiso_Rho_endcap_data_nume->SetLineColor(kBlack);
    h_PFiso_Rho_barrel_data_deno->SetLineColor(kBlack);
    h_PFiso_Rho_endcap_data_deno->SetLineColor(kBlack);
    h_PFiso_Rho_barrel_data_ctrl->SetLineColor(kBlack);
    h_PFiso_Rho_endcap_data_ctrl->SetLineColor(kBlack);
    h_pT_barrel_data_nume->SetLineColor(kBlack);
    h_pT_endcap_data_nume->SetLineColor(kBlack);
    h_pT_barrel_data_deno->SetLineColor(kBlack);
    h_pT_endcap_data_deno->SetLineColor(kBlack);
    h_pT_barrel_data_ctrl->SetLineColor(kBlack);
    h_pT_endcap_data_ctrl->SetLineColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_nume->SetLineColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_nume->SetLineColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_deno->SetLineColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_deno->SetLineColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_ctrl->SetLineColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_ctrl->SetLineColor(kBlack);
    h_dEtaInSeed_barrel_data_nume->SetLineColor(kBlack);
    h_dEtaInSeed_endcap_data_nume->SetLineColor(kBlack);
    h_dEtaInSeed_barrel_data_deno->SetLineColor(kBlack);
    h_dEtaInSeed_endcap_data_deno->SetLineColor(kBlack);
    h_dEtaInSeed_barrel_data_ctrl->SetLineColor(kBlack);
    h_dEtaInSeed_endcap_data_ctrl->SetLineColor(kBlack);
    h_dPhiIn_barrel_data_nume->SetLineColor(kBlack);
    h_dPhiIn_endcap_data_nume->SetLineColor(kBlack);
    h_dPhiIn_barrel_data_deno->SetLineColor(kBlack);
    h_dPhiIn_endcap_data_deno->SetLineColor(kBlack);
    h_dPhiIn_barrel_data_ctrl->SetLineColor(kBlack);
    h_dPhiIn_endcap_data_ctrl->SetLineColor(kBlack);
    h_HoverE_barrel_data_nume->SetLineColor(kBlack);
    h_HoverE_endcap_data_nume->SetLineColor(kBlack);
    h_HoverE_barrel_data_deno->SetLineColor(kBlack);
    h_HoverE_endcap_data_deno->SetLineColor(kBlack);
    h_HoverE_barrel_data_ctrl->SetLineColor(kBlack);
    h_HoverE_endcap_data_ctrl->SetLineColor(kBlack);
    h_InvEminusInvP_barrel_data_nume->SetLineColor(kBlack);
    h_InvEminusInvP_endcap_data_nume->SetLineColor(kBlack);
    h_InvEminusInvP_barrel_data_deno->SetLineColor(kBlack);
    h_InvEminusInvP_endcap_data_deno->SetLineColor(kBlack);
    h_InvEminusInvP_barrel_data_ctrl->SetLineColor(kBlack);
    h_InvEminusInvP_endcap_data_ctrl->SetLineColor(kBlack);
    h_chiso_barrel_data_nume->SetLineColor(kBlack);
    h_chiso_endcap_data_nume->SetLineColor(kBlack);
    h_chiso_barrel_data_deno->SetLineColor(kBlack);
    h_chiso_endcap_data_deno->SetLineColor(kBlack);
    h_chiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_chiso_endcap_data_ctrl->SetLineColor(kBlack);
    h_nhiso_barrel_data_nume->SetLineColor(kBlack);
    h_nhiso_endcap_data_nume->SetLineColor(kBlack);
    h_nhiso_barrel_data_deno->SetLineColor(kBlack);
    h_nhiso_endcap_data_deno->SetLineColor(kBlack);
    h_nhiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_nhiso_endcap_data_ctrl->SetLineColor(kBlack);
    h_phiso_barrel_data_nume->SetLineColor(kBlack);
    h_phiso_endcap_data_nume->SetLineColor(kBlack);
    h_phiso_barrel_data_deno->SetLineColor(kBlack);
    h_phiso_endcap_data_deno->SetLineColor(kBlack);
    h_phiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_phiso_endcap_data_ctrl->SetLineColor(kBlack);
    h_chisoPU_barrel_data_nume->SetLineColor(kBlack);
    h_chisoPU_endcap_data_nume->SetLineColor(kBlack);
    h_chisoPU_barrel_data_deno->SetLineColor(kBlack);
    h_chisoPU_endcap_data_deno->SetLineColor(kBlack);
    h_chisoPU_barrel_data_ctrl->SetLineColor(kBlack);
    h_chisoPU_endcap_data_ctrl->SetLineColor(kBlack);
    h_MET_data->SetLineColor(kBlack);
    h_MT_barrel_data_nume->SetLineColor(kBlack);
    h_MT_endcap_data_nume->SetLineColor(kBlack);
    h_MT_barrel_data_deno->SetLineColor(kBlack);
    h_MT_endcap_data_deno->SetLineColor(kBlack);
    h_MT_barrel_data_ctrl->SetLineColor(kBlack);
    h_MT_endcap_data_ctrl->SetLineColor(kBlack);
    h_eta_data->SetLineColor(kBlack);
    h_nVTX_data->SetLineColor(kBlack);
    h_mass_test_data->SetLineColor(kBlack);

    h_PFiso_dBeta_barrel_data_nume->SetDirectory(0);
    h_PFiso_dBeta_endcap_data_nume->SetDirectory(0);
    h_PFiso_dBeta_barrel_data_deno->SetDirectory(0);
    h_PFiso_dBeta_endcap_data_deno->SetDirectory(0);
    h_PFiso_dBeta_barrel_data_ctrl->SetDirectory(0);
    h_PFiso_dBeta_endcap_data_ctrl->SetDirectory(0);
    h_PFiso_Rho_barrel_data_nume->SetDirectory(0);
    h_PFiso_Rho_endcap_data_nume->SetDirectory(0);
    h_PFiso_Rho_barrel_data_deno->SetDirectory(0);
    h_PFiso_Rho_endcap_data_deno->SetDirectory(0);
    h_PFiso_Rho_barrel_data_ctrl->SetDirectory(0);
    h_PFiso_Rho_endcap_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume->SetDirectory(0);
    h_pT_endcap_data_nume->SetDirectory(0);
    h_pT_barrel_data_deno->SetDirectory(0);
    h_pT_endcap_data_deno->SetDirectory(0);
    h_pT_barrel_data_ctrl->SetDirectory(0);
    h_pT_endcap_data_ctrl->SetDirectory(0);
    h_SigmaIEtaIEta_barrel_data_nume->SetDirectory(0);
    h_SigmaIEtaIEta_endcap_data_nume->SetDirectory(0);
    h_SigmaIEtaIEta_barrel_data_deno->SetDirectory(0);
    h_SigmaIEtaIEta_endcap_data_deno->SetDirectory(0);
    h_SigmaIEtaIEta_barrel_data_ctrl->SetDirectory(0);
    h_SigmaIEtaIEta_endcap_data_ctrl->SetDirectory(0);
    h_dEtaInSeed_barrel_data_nume->SetDirectory(0);
    h_dEtaInSeed_endcap_data_nume->SetDirectory(0);
    h_dEtaInSeed_barrel_data_deno->SetDirectory(0);
    h_dEtaInSeed_endcap_data_deno->SetDirectory(0);
    h_dEtaInSeed_barrel_data_ctrl->SetDirectory(0);
    h_dEtaInSeed_endcap_data_ctrl->SetDirectory(0);
    h_dPhiIn_barrel_data_nume->SetDirectory(0);
    h_dPhiIn_endcap_data_nume->SetDirectory(0);
    h_dPhiIn_barrel_data_deno->SetDirectory(0);
    h_dPhiIn_endcap_data_deno->SetDirectory(0);
    h_dPhiIn_barrel_data_ctrl->SetDirectory(0);
    h_dPhiIn_endcap_data_ctrl->SetDirectory(0);
    h_HoverE_barrel_data_nume->SetDirectory(0);
    h_HoverE_endcap_data_nume->SetDirectory(0);
    h_HoverE_barrel_data_deno->SetDirectory(0);
    h_HoverE_endcap_data_deno->SetDirectory(0);
    h_HoverE_barrel_data_ctrl->SetDirectory(0);
    h_HoverE_endcap_data_ctrl->SetDirectory(0);
    h_InvEminusInvP_barrel_data_nume->SetDirectory(0);
    h_InvEminusInvP_endcap_data_nume->SetDirectory(0);
    h_InvEminusInvP_barrel_data_deno->SetDirectory(0);
    h_InvEminusInvP_endcap_data_deno->SetDirectory(0);
    h_InvEminusInvP_barrel_data_ctrl->SetDirectory(0);
    h_InvEminusInvP_endcap_data_ctrl->SetDirectory(0);
    h_chiso_barrel_data_nume->SetDirectory(0);
    h_chiso_endcap_data_nume->SetDirectory(0);
    h_chiso_barrel_data_deno->SetDirectory(0);
    h_chiso_endcap_data_deno->SetDirectory(0);
    h_chiso_barrel_data_ctrl->SetDirectory(0);
    h_chiso_endcap_data_ctrl->SetDirectory(0);
    h_nhiso_barrel_data_nume->SetDirectory(0);
    h_nhiso_endcap_data_nume->SetDirectory(0);
    h_nhiso_barrel_data_deno->SetDirectory(0);
    h_nhiso_endcap_data_deno->SetDirectory(0);
    h_nhiso_barrel_data_ctrl->SetDirectory(0);
    h_nhiso_endcap_data_ctrl->SetDirectory(0);
    h_phiso_barrel_data_nume->SetDirectory(0);
    h_phiso_endcap_data_nume->SetDirectory(0);
    h_phiso_barrel_data_deno->SetDirectory(0);
    h_phiso_endcap_data_deno->SetDirectory(0);
    h_phiso_barrel_data_ctrl->SetDirectory(0);
    h_phiso_endcap_data_ctrl->SetDirectory(0);
    h_chisoPU_barrel_data_nume->SetDirectory(0);
    h_chisoPU_endcap_data_nume->SetDirectory(0);
    h_chisoPU_barrel_data_deno->SetDirectory(0);
    h_chisoPU_endcap_data_deno->SetDirectory(0);
    h_chisoPU_barrel_data_ctrl->SetDirectory(0);
    h_chisoPU_endcap_data_ctrl->SetDirectory(0);
    h_MET_data->SetDirectory(0);
    h_MT_barrel_data_nume->SetDirectory(0);
    h_MT_endcap_data_nume->SetDirectory(0);
    h_MT_barrel_data_deno->SetDirectory(0);
    h_MT_endcap_data_deno->SetDirectory(0);
    h_MT_barrel_data_ctrl->SetDirectory(0);
    h_MT_endcap_data_ctrl->SetDirectory(0);
    h_eta_data->SetDirectory(0);
    h_nVTX_data->SetDirectory(0);
    h_mass_test_data->SetDirectory(0);

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_PFiso_dBeta_barrel_nume = new myRatioPlot_t("RP_PFiso_dBeta_barrel_nume", s_PFiso_dBeta_barrel_nume, h_PFiso_dBeta_barrel_data_nume);
    myRatioPlot_t *RP_PFiso_dBeta_endcap_nume = new myRatioPlot_t("RP_PFiso_dBeta_endcap_nume", s_PFiso_dBeta_endcap_nume, h_PFiso_dBeta_endcap_data_nume);
    myRatioPlot_t *RP_PFiso_dBeta_barrel_deno = new myRatioPlot_t("RP_PFiso_dBeta_barrel_deno", s_PFiso_dBeta_barrel_deno, h_PFiso_dBeta_barrel_data_deno);
    myRatioPlot_t *RP_PFiso_dBeta_endcap_deno = new myRatioPlot_t("RP_PFiso_dBeta_endcap_deno", s_PFiso_dBeta_endcap_deno, h_PFiso_dBeta_endcap_data_deno);
    myRatioPlot_t *RP_PFiso_dBeta_barrel_ctrl = new myRatioPlot_t("RP_PFiso_dBeta_barrel_ctrl", s_PFiso_dBeta_barrel_ctrl, h_PFiso_dBeta_barrel_data_ctrl);
    myRatioPlot_t *RP_PFiso_dBeta_endcap_ctrl = new myRatioPlot_t("RP_PFiso_dBeta_endcap_ctrl", s_PFiso_dBeta_endcap_ctrl, h_PFiso_dBeta_endcap_data_ctrl);
    myRatioPlot_t *RP_PFiso_Rho_barrel_nume = new myRatioPlot_t("RP_PFiso_Rho_barrel_nume", s_PFiso_Rho_barrel_nume, h_PFiso_Rho_barrel_data_nume);
    myRatioPlot_t *RP_PFiso_Rho_endcap_nume = new myRatioPlot_t("RP_PFiso_Rho_endcap_nume", s_PFiso_Rho_endcap_nume, h_PFiso_Rho_endcap_data_nume);
    myRatioPlot_t *RP_PFiso_Rho_barrel_deno = new myRatioPlot_t("RP_PFiso_Rho_barrel_deno", s_PFiso_Rho_barrel_deno, h_PFiso_Rho_barrel_data_deno);
    myRatioPlot_t *RP_PFiso_Rho_endcap_deno = new myRatioPlot_t("RP_PFiso_Rho_endcap_deno", s_PFiso_Rho_endcap_deno, h_PFiso_Rho_endcap_data_deno);
    myRatioPlot_t *RP_PFiso_Rho_barrel_ctrl = new myRatioPlot_t("RP_PFiso_Rho_barrel_ctrl", s_PFiso_Rho_barrel_ctrl, h_PFiso_Rho_barrel_data_ctrl);
    myRatioPlot_t *RP_PFiso_Rho_endcap_ctrl = new myRatioPlot_t("RP_PFiso_Rho_endcap_ctrl", s_PFiso_Rho_endcap_ctrl, h_PFiso_Rho_endcap_data_ctrl);
    myRatioPlot_t *RP_pT_barrel_nume = new myRatioPlot_t("RP_pT_barrel_nume", s_pT_barrel_nume, h_pT_barrel_data_nume);
    myRatioPlot_t *RP_pT_endcap_nume = new myRatioPlot_t("RP_pT_endcap_nume", s_pT_endcap_nume, h_pT_endcap_data_nume);
    myRatioPlot_t *RP_pT_barrel_deno = new myRatioPlot_t("RP_pT_barrel_deno", s_pT_barrel_deno, h_pT_barrel_data_deno);
    myRatioPlot_t *RP_pT_endcap_deno = new myRatioPlot_t("RP_pT_endcap_deno", s_pT_endcap_deno, h_pT_endcap_data_deno);
    myRatioPlot_t *RP_pT_barrel_ctrl = new myRatioPlot_t("RP_pT_barrel_ctrl", s_pT_barrel_ctrl, h_pT_barrel_data_ctrl);
    myRatioPlot_t *RP_pT_endcap_ctrl = new myRatioPlot_t("RP_pT_endcap_ctrl", s_pT_endcap_ctrl, h_pT_endcap_data_ctrl);
    myRatioPlot_t *RP_SigmaIEtaIEta_barrel_nume = new myRatioPlot_t("RP_SigmaIEtaIEta_barrel_nume", s_SigmaIEtaIEta_barrel_nume, h_SigmaIEtaIEta_barrel_data_nume);
    myRatioPlot_t *RP_SigmaIEtaIEta_endcap_nume = new myRatioPlot_t("RP_SigmaIEtaIEta_endcap_nume", s_SigmaIEtaIEta_endcap_nume, h_SigmaIEtaIEta_endcap_data_nume);
    myRatioPlot_t *RP_SigmaIEtaIEta_barrel_deno = new myRatioPlot_t("RP_SigmaIEtaIEta_barrel_deno", s_SigmaIEtaIEta_barrel_deno, h_SigmaIEtaIEta_barrel_data_deno);
    myRatioPlot_t *RP_SigmaIEtaIEta_endcap_deno = new myRatioPlot_t("RP_SigmaIEtaIEta_endcap_deno", s_SigmaIEtaIEta_endcap_deno, h_SigmaIEtaIEta_endcap_data_deno);
    myRatioPlot_t *RP_SigmaIEtaIEta_barrel_ctrl = new myRatioPlot_t("RP_SigmaIEtaIEta_barrel_ctrl", s_SigmaIEtaIEta_barrel_ctrl, h_SigmaIEtaIEta_barrel_data_ctrl);
    myRatioPlot_t *RP_SigmaIEtaIEta_endcap_ctrl = new myRatioPlot_t("RP_SigmaIEtaIEta_endcap_ctrl", s_SigmaIEtaIEta_endcap_ctrl, h_SigmaIEtaIEta_endcap_data_ctrl);
    myRatioPlot_t *RP_dEtaInSeed_barrel_nume = new myRatioPlot_t("RP_dEtaInSeed_barrel_nume", s_dEtaInSeed_barrel_nume, h_dEtaInSeed_barrel_data_nume);
    myRatioPlot_t *RP_dEtaInSeed_endcap_nume = new myRatioPlot_t("RP_dEtaInSeed_endcap_nume", s_dEtaInSeed_endcap_nume, h_dEtaInSeed_endcap_data_nume);
    myRatioPlot_t *RP_dEtaInSeed_barrel_deno = new myRatioPlot_t("RP_dEtaInSeed_barrel_deno", s_dEtaInSeed_barrel_deno, h_dEtaInSeed_barrel_data_deno);
    myRatioPlot_t *RP_dEtaInSeed_endcap_deno = new myRatioPlot_t("RP_dEtaInSeed_endcap_deno", s_dEtaInSeed_endcap_deno, h_dEtaInSeed_endcap_data_deno);
    myRatioPlot_t *RP_dEtaInSeed_barrel_ctrl = new myRatioPlot_t("RP_dEtaInSeed_barrel_ctrl", s_dEtaInSeed_barrel_ctrl, h_dEtaInSeed_barrel_data_ctrl);
    myRatioPlot_t *RP_dEtaInSeed_endcap_ctrl = new myRatioPlot_t("RP_dEtaInSeed_endcap_ctrl", s_dEtaInSeed_endcap_ctrl, h_dEtaInSeed_endcap_data_ctrl);
    myRatioPlot_t *RP_dPhiIn_barrel_nume = new myRatioPlot_t("RP_dPhiIn_barrel_nume", s_dPhiIn_barrel_nume, h_dPhiIn_barrel_data_nume);
    myRatioPlot_t *RP_dPhiIn_endcap_nume = new myRatioPlot_t("RP_dPhiIn_endcap_nume", s_dPhiIn_endcap_nume, h_dPhiIn_endcap_data_nume);
    myRatioPlot_t *RP_dPhiIn_barrel_deno = new myRatioPlot_t("RP_dPhiIn_barrel_deno", s_dPhiIn_barrel_deno, h_dPhiIn_barrel_data_deno);
    myRatioPlot_t *RP_dPhiIn_endcap_deno = new myRatioPlot_t("RP_dPhiIn_endcap_deno", s_dPhiIn_endcap_deno, h_dPhiIn_endcap_data_deno);
    myRatioPlot_t *RP_dPhiIn_barrel_ctrl = new myRatioPlot_t("RP_dPhiIn_barrel_ctrl", s_dPhiIn_barrel_ctrl, h_dPhiIn_barrel_data_ctrl);
    myRatioPlot_t *RP_dPhiIn_endcap_ctrl = new myRatioPlot_t("RP_dPhiIn_endcap_ctrl", s_dPhiIn_endcap_ctrl, h_dPhiIn_endcap_data_ctrl);
    myRatioPlot_t *RP_HoverE_barrel_nume = new myRatioPlot_t("RP_HoverE_barrel_nume", s_HoverE_barrel_nume, h_HoverE_barrel_data_nume);
    myRatioPlot_t *RP_HoverE_endcap_nume = new myRatioPlot_t("RP_HoverE_endcap_nume", s_HoverE_endcap_nume, h_HoverE_endcap_data_nume);
    myRatioPlot_t *RP_HoverE_barrel_deno = new myRatioPlot_t("RP_HoverE_barrel_deno", s_HoverE_barrel_deno, h_HoverE_barrel_data_deno);
    myRatioPlot_t *RP_HoverE_endcap_deno = new myRatioPlot_t("RP_HoverE_endcap_deno", s_HoverE_endcap_deno, h_HoverE_endcap_data_deno);
    myRatioPlot_t *RP_HoverE_barrel_ctrl = new myRatioPlot_t("RP_HoverE_barrel_ctrl", s_HoverE_barrel_ctrl, h_HoverE_barrel_data_ctrl);
    myRatioPlot_t *RP_HoverE_endcap_ctrl = new myRatioPlot_t("RP_HoverE_endcap_ctrl", s_HoverE_endcap_ctrl, h_HoverE_endcap_data_ctrl);
    myRatioPlot_t *RP_InvEminusInvP_barrel_nume = new myRatioPlot_t("RP_InvEminusInvP_barrel_nume", s_InvEminusInvP_barrel_nume, h_InvEminusInvP_barrel_data_nume);
    myRatioPlot_t *RP_InvEminusInvP_endcap_nume = new myRatioPlot_t("RP_InvEminusInvP_endcap_nume", s_InvEminusInvP_endcap_nume, h_InvEminusInvP_endcap_data_nume);
    myRatioPlot_t *RP_InvEminusInvP_barrel_deno = new myRatioPlot_t("RP_InvEminusInvP_barrel_deno", s_InvEminusInvP_barrel_deno, h_InvEminusInvP_barrel_data_deno);
    myRatioPlot_t *RP_InvEminusInvP_endcap_deno = new myRatioPlot_t("RP_InvEminusInvP_endcap_deno", s_InvEminusInvP_endcap_deno, h_InvEminusInvP_endcap_data_deno);
    myRatioPlot_t *RP_InvEminusInvP_barrel_ctrl = new myRatioPlot_t("RP_InvEminusInvP_barrel_ctrl", s_InvEminusInvP_barrel_ctrl, h_InvEminusInvP_barrel_data_ctrl);
    myRatioPlot_t *RP_InvEminusInvP_endcap_ctrl = new myRatioPlot_t("RP_InvEminusInvP_endcap_ctrl", s_InvEminusInvP_endcap_ctrl, h_InvEminusInvP_endcap_data_ctrl);
    myRatioPlot_t *RP_chiso_barrel_nume = new myRatioPlot_t("RP_chiso_barrel_nume", s_chiso_barrel_nume, h_chiso_barrel_data_nume);
    myRatioPlot_t *RP_chiso_endcap_nume = new myRatioPlot_t("RP_chiso_endcap_nume", s_chiso_endcap_nume, h_chiso_endcap_data_nume);
    myRatioPlot_t *RP_chiso_barrel_deno = new myRatioPlot_t("RP_chiso_barrel_deno", s_chiso_barrel_deno, h_chiso_barrel_data_deno);
    myRatioPlot_t *RP_chiso_endcap_deno = new myRatioPlot_t("RP_chiso_endcap_deno", s_chiso_endcap_deno, h_chiso_endcap_data_deno);
    myRatioPlot_t *RP_chiso_barrel_ctrl = new myRatioPlot_t("RP_chiso_barrel_ctrl", s_chiso_barrel_ctrl, h_chiso_barrel_data_ctrl);
    myRatioPlot_t *RP_chiso_endcap_ctrl = new myRatioPlot_t("RP_chiso_endcap_ctrl", s_chiso_endcap_ctrl, h_chiso_endcap_data_ctrl);
    myRatioPlot_t *RP_nhiso_barrel_nume = new myRatioPlot_t("RP_nhiso_barrel_nume", s_nhiso_barrel_nume, h_nhiso_barrel_data_nume);
    myRatioPlot_t *RP_nhiso_endcap_nume = new myRatioPlot_t("RP_nhiso_endcap_nume", s_nhiso_endcap_nume, h_nhiso_endcap_data_nume);
    myRatioPlot_t *RP_nhiso_barrel_deno = new myRatioPlot_t("RP_nhiso_barrel_deno", s_nhiso_barrel_deno, h_nhiso_barrel_data_deno);
    myRatioPlot_t *RP_nhiso_endcap_deno = new myRatioPlot_t("RP_nhiso_endcap_deno", s_nhiso_endcap_deno, h_nhiso_endcap_data_deno);
    myRatioPlot_t *RP_nhiso_barrel_ctrl = new myRatioPlot_t("RP_nhiso_barrel_ctrl", s_nhiso_barrel_ctrl, h_nhiso_barrel_data_ctrl);
    myRatioPlot_t *RP_nhiso_endcap_ctrl = new myRatioPlot_t("RP_nhiso_endcap_ctrl", s_nhiso_endcap_ctrl, h_nhiso_endcap_data_ctrl);
    myRatioPlot_t *RP_phiso_barrel_nume = new myRatioPlot_t("RP_phiso_barrel_nume", s_phiso_barrel_nume, h_phiso_barrel_data_nume);
    myRatioPlot_t *RP_phiso_endcap_nume = new myRatioPlot_t("RP_phiso_endcap_nume", s_phiso_endcap_nume, h_phiso_endcap_data_nume);
    myRatioPlot_t *RP_phiso_barrel_deno = new myRatioPlot_t("RP_phiso_barrel_deno", s_phiso_barrel_deno, h_phiso_barrel_data_deno);
    myRatioPlot_t *RP_phiso_endcap_deno = new myRatioPlot_t("RP_phiso_endcap_deno", s_phiso_endcap_deno, h_phiso_endcap_data_deno);
    myRatioPlot_t *RP_phiso_barrel_ctrl = new myRatioPlot_t("RP_phiso_barrel_ctrl", s_phiso_barrel_ctrl, h_phiso_barrel_data_ctrl);
    myRatioPlot_t *RP_phiso_endcap_ctrl = new myRatioPlot_t("RP_phiso_endcap_ctrl", s_phiso_endcap_ctrl, h_phiso_endcap_data_ctrl);
    myRatioPlot_t *RP_chisoPU_barrel_nume = new myRatioPlot_t("RP_chisoPU_barrel_nume", s_chisoPU_barrel_nume, h_chisoPU_barrel_data_nume);
    myRatioPlot_t *RP_chisoPU_endcap_nume = new myRatioPlot_t("RP_chisoPU_endcap_nume", s_chisoPU_endcap_nume, h_chisoPU_endcap_data_nume);
    myRatioPlot_t *RP_chisoPU_barrel_deno = new myRatioPlot_t("RP_chisoPU_barrel_deno", s_chisoPU_barrel_deno, h_chisoPU_barrel_data_deno);
    myRatioPlot_t *RP_chisoPU_endcap_deno = new myRatioPlot_t("RP_chisoPU_endcap_deno", s_chisoPU_endcap_deno, h_chisoPU_endcap_data_deno);
    myRatioPlot_t *RP_chisoPU_barrel_ctrl = new myRatioPlot_t("RP_chisoPU_barrel_ctrl", s_chisoPU_barrel_ctrl, h_chisoPU_barrel_data_ctrl);
    myRatioPlot_t *RP_chisoPU_endcap_ctrl = new myRatioPlot_t("RP_chisoPU_endcap_ctrl", s_chisoPU_endcap_ctrl, h_chisoPU_endcap_data_ctrl);
    myRatioPlot_t *RP_MET = new myRatioPlot_t("RP_MET", s_MET, h_MET_data);
    myRatioPlot_t *RP_MT_barrel_nume = new myRatioPlot_t("RP_MT_barrel_nume", s_MT_barrel_nume, h_MT_barrel_data_nume);
    myRatioPlot_t *RP_MT_endcap_nume = new myRatioPlot_t("RP_MT_endcap_nume", s_MT_endcap_nume, h_MT_endcap_data_nume);
    myRatioPlot_t *RP_MT_barrel_deno = new myRatioPlot_t("RP_MT_barrel_deno", s_MT_barrel_deno, h_MT_barrel_data_deno);
    myRatioPlot_t *RP_MT_endcap_deno = new myRatioPlot_t("RP_MT_endcap_deno", s_MT_endcap_deno, h_MT_endcap_data_deno);
    myRatioPlot_t *RP_MT_barrel_ctrl = new myRatioPlot_t("RP_MT_barrel_ctrl", s_MT_barrel_ctrl, h_MT_barrel_data_ctrl);
    myRatioPlot_t *RP_MT_endcap_ctrl = new myRatioPlot_t("RP_MT_endcap_ctrl", s_MT_endcap_ctrl, h_MT_endcap_data_ctrl);
    myRatioPlot_t *RP_eta = new myRatioPlot_t("RP_eta", s_eta, h_eta_data);
    myRatioPlot_t *RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_nVTX_data);
    myRatioPlot_t *RP_mass_test = new myRatioPlot_t("RP_mass_test", s_mass_test, h_mass_test_data);

    RP_PFiso_dBeta_barrel_nume->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.3);
    RP_PFiso_dBeta_endcap_nume->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.3);
    RP_PFiso_dBeta_barrel_deno->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_PFiso_dBeta_endcap_deno->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_PFiso_dBeta_barrel_ctrl->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{barrel}}^{control})", 0, 5);
    RP_PFiso_dBeta_endcap_ctrl->SetPlots("relPFiso_dBeta (e_{#lower[-0.4]{endcap}}^{control})", 0, 5);
    RP_PFiso_Rho_barrel_nume->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.2);
    RP_PFiso_Rho_endcap_nume->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.2);
    RP_PFiso_Rho_barrel_deno->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_PFiso_Rho_endcap_deno->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_PFiso_Rho_barrel_ctrl->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{control})", 0, 5);
    RP_PFiso_Rho_endcap_ctrl->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{control})", 0, 5);
    RP_pT_barrel_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{nume}) [GeV/c]", 25, 1000);
    RP_pT_endcap_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{nume}) [GeV/c]", 25, 1000);
    RP_pT_barrel_deno->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{deno}) [GeV/c]", 25, 1000);
    RP_pT_endcap_deno->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{deno}) [GeV/c]", 25, 1000);
    RP_pT_barrel_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{control}) [GeV/c]", 25, 1000);
    RP_pT_endcap_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{control}) [GeV/c]", 25, 1000);
    RP_SigmaIEtaIEta_barrel_nume->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{barrel}}^{nume})",    0, 0.01);
    RP_SigmaIEtaIEta_endcap_nume->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{endcap}}^{nume})",    0, 0.03);
    RP_SigmaIEtaIEta_barrel_deno->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{barrel}}^{deno})",    0, 0.05);
    RP_SigmaIEtaIEta_endcap_deno->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{endcap}}^{deno})",    0, 0.1);
    RP_SigmaIEtaIEta_barrel_ctrl->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{barrel}}^{control})", 0, 0.05);
    RP_SigmaIEtaIEta_endcap_ctrl->SetPlots("#sigma_{i#etai#eta} (e_{#lower[-0.4]{endcap}}^{control})", 0, 0.1);
    RP_dEtaInSeed_barrel_nume->SetPlots("dEtaInSeed (e_{#lower[-0.4]{barrel}}^{nume})", -0.1, 0.1);
    RP_dEtaInSeed_endcap_nume->SetPlots("dEtaInSeed (e_{#lower[-0.4]{endcap}}^{nume})", -0.1, 0.1);
    RP_dEtaInSeed_barrel_deno->SetPlots("dEtaInSeed (e_{#lower[-0.4]{barrel}}^{deno})", -0.1, 0.1);
    RP_dEtaInSeed_endcap_deno->SetPlots("dEtaInSeed (e_{#lower[-0.4]{endcap}}^{deno})", -1, 1);
    RP_dEtaInSeed_barrel_ctrl->SetPlots("dEtaInSeed (e_{#lower[-0.4]{barrel}}^{control})", -0.1, 0.1);
    RP_dEtaInSeed_endcap_ctrl->SetPlots("dEtaInSeed (e_{#lower[-0.4]{endcap}}^{control})", -1, 1);
    RP_dPhiIn_barrel_nume->SetPlots("dPhiIn (e_{#lower[-0.4]{barrel}}^{nume})", -0.1, 0.1);
    RP_dPhiIn_endcap_nume->SetPlots("dPhiIn (e_{#lower[-0.4]{endcap}}^{nume})", -0.1, 0.1);
    RP_dPhiIn_barrel_deno->SetPlots("dPhiIn (e_{#lower[-0.4]{barrel}}^{deno})", -1, 1);
    RP_dPhiIn_endcap_deno->SetPlots("dPhiIn (e_{#lower[-0.4]{endcap}}^{deno})", -1, 1);
    RP_dPhiIn_barrel_ctrl->SetPlots("dPhiIn (e_{#lower[-0.4]{barrel}}^{control})", -1, 1);
    RP_dPhiIn_endcap_ctrl->SetPlots("dPhiIn (e_{#lower[-0.4]{endcap}}^{control})", -1, 1);
    RP_HoverE_barrel_nume->SetPlots("H/E (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.1);
    RP_HoverE_endcap_nume->SetPlots("H/E (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.15);
    RP_HoverE_barrel_deno->SetPlots("H/E (e_{#lower[-0.4]{barrel}}^{deno})", 0, 0.1);
    RP_HoverE_endcap_deno->SetPlots("H/E (e_{#lower[-0.4]{endcap}}^{deno})", 0, 0.15);
    RP_HoverE_barrel_ctrl->SetPlots("H/E (e_{#lower[-0.4]{barrel}}^{control})", 0, 0.1);
    RP_HoverE_endcap_ctrl->SetPlots("H/E (e_{#lower[-0.4]{endcap}}^{control})", 0, 0.15);
    RP_InvEminusInvP_barrel_nume->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.5);
    RP_InvEminusInvP_endcap_nume->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.5);
    RP_InvEminusInvP_barrel_deno->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{deno})", 0, 6);
    RP_InvEminusInvP_endcap_deno->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{deno})", 0, 6);
    RP_InvEminusInvP_barrel_ctrl->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{control})", 0, 6);
    RP_InvEminusInvP_endcap_ctrl->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{control})", 0, 6);
    RP_chiso_barrel_nume->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.2);
    RP_chiso_endcap_nume->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.2);
    RP_chiso_barrel_deno->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{deno})", 0, 10);
    RP_chiso_endcap_deno->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{deno})", 0, 10);
    RP_chiso_barrel_ctrl->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{control})", 0, 10);
    RP_chiso_endcap_ctrl->SetPlots("I_{ch. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{control})", 0, 10);
    RP_nhiso_barrel_nume->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.2);
    RP_nhiso_endcap_nume->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.2);
    RP_nhiso_barrel_deno->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{deno})", 0, 10);
    RP_nhiso_endcap_deno->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{deno})", 0, 10);
    RP_nhiso_barrel_ctrl->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{barrel}}^{control})", 0, 10);
    RP_nhiso_endcap_ctrl->SetPlots("I_{n. had.}^{rel.} (e_{#lower[-0.4]{endcap}}^{control})", 0, 10);
    RP_phiso_barrel_nume->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.2);
    RP_phiso_endcap_nume->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.2);
    RP_phiso_barrel_deno->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{barrel}}^{deno})", 0, 10);
    RP_phiso_endcap_deno->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{endcap}}^{deno})", 0, 10);
    RP_phiso_barrel_ctrl->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{barrel}}^{control})", 0, 10);
    RP_phiso_endcap_ctrl->SetPlots("I_{ph.}^{rel.} (e_{#lower[-0.4]{endcap}}^{control})", 0, 10);
    RP_chisoPU_barrel_nume->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{barrel}}^{nume})", 0, 1);
    RP_chisoPU_endcap_nume->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{endcap}}^{nume})", 0, 1);
    RP_chisoPU_barrel_deno->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_chisoPU_endcap_deno->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_chisoPU_barrel_ctrl->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{barrel}}^{control})", 0, 5);
    RP_chisoPU_endcap_ctrl->SetPlots("I_{ch. had. from PU}^{rel.} (e_{#lower[-0.4]{endcap}}^{control})", 0, 5);
    RP_MET->SetPlots("E_{#lower[-0.25]{T}}^{miss} [GeV]", 0, 500);
    RP_MT_barrel_nume->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{nume}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_nume->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{nume}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_barrel_deno->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{deno}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_deno->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{deno}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_barrel_ctrl->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{control}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_ctrl->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{control}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_eta->SetPlots("#eta (e)", -3, 3);
    RP_nVTX->SetPlots("N_{#lower[-0.25]{VTX}}", 0, 50);
    RP_mass_test->SetPlots("m_{#lower[-0.25]{ee}} [GeV/c^{2}]", 15, 3000);


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
    legend->AddEntry(h_eta_MC[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->AddEntry(h_eta_MC[_GJets_20to100], "#gamma+Jets", "f");
    legend->SetNColumns(2);

    RP_PFiso_dBeta_barrel_nume->ImportLegend(legend);
    RP_PFiso_dBeta_endcap_nume->ImportLegend(legend);
    RP_PFiso_dBeta_barrel_deno->ImportLegend(legend);
    RP_PFiso_dBeta_endcap_deno->ImportLegend(legend);
    RP_PFiso_dBeta_barrel_ctrl->ImportLegend(legend);
    RP_PFiso_dBeta_endcap_ctrl->ImportLegend(legend);
    RP_PFiso_Rho_barrel_nume->ImportLegend(legend);
    RP_PFiso_Rho_endcap_nume->ImportLegend(legend);
    RP_PFiso_Rho_barrel_deno->ImportLegend(legend);
    RP_PFiso_Rho_endcap_deno->ImportLegend(legend);
    RP_PFiso_Rho_barrel_ctrl->ImportLegend(legend);
    RP_PFiso_Rho_endcap_ctrl->ImportLegend(legend);
    RP_pT_barrel_nume->ImportLegend(legend);
    RP_pT_endcap_nume->ImportLegend(legend);
    RP_pT_barrel_deno->ImportLegend(legend);
    RP_pT_endcap_deno->ImportLegend(legend);
    RP_pT_barrel_ctrl->ImportLegend(legend);
    RP_pT_endcap_ctrl->ImportLegend(legend);
    RP_SigmaIEtaIEta_barrel_nume->ImportLegend(legend);
    RP_SigmaIEtaIEta_endcap_nume->ImportLegend(legend);
    RP_SigmaIEtaIEta_barrel_deno->ImportLegend(legend);
    RP_SigmaIEtaIEta_endcap_deno->ImportLegend(legend);
    RP_SigmaIEtaIEta_barrel_ctrl->ImportLegend(legend);
    RP_SigmaIEtaIEta_endcap_ctrl->ImportLegend(legend);
    RP_dEtaInSeed_barrel_nume->ImportLegend(legend);
    RP_dEtaInSeed_endcap_nume->ImportLegend(legend);
    RP_dEtaInSeed_barrel_deno->ImportLegend(legend);
    RP_dEtaInSeed_endcap_deno->ImportLegend(legend);
    RP_dEtaInSeed_barrel_ctrl->ImportLegend(legend);
    RP_dEtaInSeed_endcap_ctrl->ImportLegend(legend);
    RP_dPhiIn_barrel_nume->ImportLegend(legend);
    RP_dPhiIn_endcap_nume->ImportLegend(legend);
    RP_dPhiIn_barrel_deno->ImportLegend(legend);
    RP_dPhiIn_endcap_deno->ImportLegend(legend);
    RP_dPhiIn_barrel_ctrl->ImportLegend(legend);
    RP_dPhiIn_endcap_ctrl->ImportLegend(legend);
    RP_HoverE_barrel_nume->ImportLegend(legend);
    RP_HoverE_endcap_nume->ImportLegend(legend);
    RP_HoverE_barrel_deno->ImportLegend(legend);
    RP_HoverE_endcap_deno->ImportLegend(legend);
    RP_HoverE_barrel_ctrl->ImportLegend(legend);
    RP_HoverE_endcap_ctrl->ImportLegend(legend);
    RP_InvEminusInvP_barrel_nume->ImportLegend(legend);
    RP_InvEminusInvP_endcap_nume->ImportLegend(legend);
    RP_InvEminusInvP_barrel_deno->ImportLegend(legend);
    RP_InvEminusInvP_endcap_deno->ImportLegend(legend);
    RP_InvEminusInvP_barrel_ctrl->ImportLegend(legend);
    RP_InvEminusInvP_endcap_ctrl->ImportLegend(legend);
    RP_chiso_barrel_nume->ImportLegend(legend);
    RP_chiso_endcap_nume->ImportLegend(legend);
    RP_chiso_barrel_deno->ImportLegend(legend);
    RP_chiso_endcap_deno->ImportLegend(legend);
    RP_chiso_barrel_ctrl->ImportLegend(legend);
    RP_chiso_endcap_ctrl->ImportLegend(legend);
    RP_nhiso_barrel_nume->ImportLegend(legend);
    RP_nhiso_endcap_nume->ImportLegend(legend);
    RP_nhiso_barrel_deno->ImportLegend(legend);
    RP_nhiso_endcap_deno->ImportLegend(legend);
    RP_nhiso_barrel_ctrl->ImportLegend(legend);
    RP_nhiso_endcap_ctrl->ImportLegend(legend);
    RP_phiso_barrel_nume->ImportLegend(legend);
    RP_phiso_endcap_nume->ImportLegend(legend);
    RP_phiso_barrel_deno->ImportLegend(legend);
    RP_phiso_endcap_deno->ImportLegend(legend);
    RP_phiso_barrel_ctrl->ImportLegend(legend);
    RP_phiso_endcap_ctrl->ImportLegend(legend);
    RP_chisoPU_barrel_nume->ImportLegend(legend);
    RP_chisoPU_endcap_nume->ImportLegend(legend);
    RP_chisoPU_barrel_deno->ImportLegend(legend);
    RP_chisoPU_endcap_deno->ImportLegend(legend);
    RP_chisoPU_barrel_ctrl->ImportLegend(legend);
    RP_chisoPU_endcap_ctrl->ImportLegend(legend);
    RP_MET->ImportLegend(legend);
    RP_MT_barrel_nume->ImportLegend(legend);
    RP_MT_endcap_nume->ImportLegend(legend);
    RP_MT_barrel_deno->ImportLegend(legend);
    RP_MT_endcap_deno->ImportLegend(legend);
    RP_MT_barrel_ctrl->ImportLegend(legend);
    RP_MT_endcap_ctrl->ImportLegend(legend);
    RP_eta->ImportLegend(legend);
    RP_nVTX->ImportLegend(legend);
    RP_mass_test->ImportLegend(legend);

    RP_PFiso_dBeta_barrel_nume->Draw(1, 1e10, 0);
    RP_PFiso_dBeta_endcap_nume->Draw(1, 1e10, 0);
//    RP_PFiso_dBeta_barrel_deno->Draw(1, 1e10, 0);
//    RP_PFiso_dBeta_endcap_deno->Draw(1, 1e10, 0);
    RP_PFiso_dBeta_barrel_ctrl->Draw(1, 1e10, 0);
    RP_PFiso_dBeta_endcap_ctrl->Draw(1, 1e10, 0);
    RP_PFiso_Rho_barrel_nume->Draw(1, 1e10, 0);
    RP_PFiso_Rho_endcap_nume->Draw(1, 1e10, 0);
//    RP_PFiso_Rho_barrel_deno->Draw(1, 1e10, 0);
//    RP_PFiso_Rho_endcap_deno->Draw(1, 1e10, 0);
    RP_PFiso_Rho_barrel_ctrl->Draw(1, 1e10, 0);
    RP_PFiso_Rho_endcap_ctrl->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_barrel_nume->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_endcap_nume->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_barrel_deno->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_endcap_deno->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_barrel_ctrl->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_endcap_ctrl->Draw(1, 1e10, 0);
    RP_dEtaInSeed_barrel_nume->Draw(1, 1e10, 0);
    RP_dEtaInSeed_endcap_nume->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_barrel_deno->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_endcap_deno->Draw(1, 1e10, 0);
    RP_dEtaInSeed_barrel_ctrl->Draw(1, 1e10, 0);
    RP_dEtaInSeed_endcap_ctrl->Draw(1, 1e10, 0);
    RP_dPhiIn_barrel_nume->Draw(1, 1e10, 0);
    RP_dPhiIn_endcap_nume->Draw(1, 1e10, 0);
//    RP_dPhiIn_barrel_deno->Draw(1, 1e10, 0);
//    RP_dPhiIn_endcap_deno->Draw(1, 1e10, 0);
    RP_dPhiIn_barrel_ctrl->Draw(1, 1e10, 0);
    RP_dPhiIn_endcap_ctrl->Draw(1, 1e10, 0);
    RP_HoverE_barrel_nume->Draw(1, 1e10, 0);
    RP_HoverE_endcap_nume->Draw(1, 1e10, 0);
//    RP_HoverE_barrel_deno->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_deno->Draw(1, 1e10, 0);
    RP_HoverE_barrel_ctrl->Draw(1, 1e10, 0);
    RP_HoverE_endcap_ctrl->Draw(1, 1e10, 0);
    RP_InvEminusInvP_barrel_nume->Draw(1, 1e10, 0);
    RP_InvEminusInvP_endcap_nume->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_barrel_deno->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_endcap_deno->Draw(1, 1e10, 0);
    RP_InvEminusInvP_barrel_ctrl->Draw(1, 1e10, 0);
    RP_InvEminusInvP_endcap_ctrl->Draw(1, 1e10, 0);
    RP_chiso_barrel_nume->Draw(1, 1e10, 0);
    RP_chiso_endcap_nume->Draw(1, 1e10, 0);
//    RP_chiso_barrel_deno->Draw(1, 1e10, 0);
//    RP_chiso_endcap_deno->Draw(1, 1e10, 0);
    RP_chiso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_chiso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_nhiso_barrel_nume->Draw(1, 1e10, 0);
    RP_nhiso_endcap_nume->Draw(1, 1e10, 0);
//    RP_nhiso_barrel_deno->Draw(1, 1e10, 0);
//    RP_nhiso_endcap_deno->Draw(1, 1e10, 0);
    RP_nhiso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_nhiso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_phiso_barrel_nume->Draw(1, 1e10, 0);
    RP_phiso_endcap_nume->Draw(1, 1e10, 0);
//    RP_phiso_barrel_deno->Draw(1, 1e10, 0);
//    RP_phiso_endcap_deno->Draw(1, 1e10, 0);
    RP_phiso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_phiso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_chisoPU_barrel_nume->Draw(1, 1e10, 0);
    RP_chisoPU_endcap_nume->Draw(1, 1e10, 0);
//    RP_chisoPU_barrel_deno->Draw(1, 1e10, 0);
//    RP_chisoPU_endcap_deno->Draw(1, 1e10, 0);
    RP_chisoPU_barrel_ctrl->Draw(1, 1e10, 0);
    RP_chisoPU_endcap_ctrl->Draw(1, 1e10, 0);
    RP_pT_barrel_nume->Draw(1, 1e10, 0);
    RP_pT_endcap_nume->Draw(1, 1e10, 0);
//    RP_pT_barrel_deno->Draw(1, 1e10, 0);
//    RP_pT_endcap_deno->Draw(1, 1e10, 0);
    RP_pT_barrel_ctrl->Draw(1, 1e10, 0);
    RP_pT_endcap_ctrl->Draw(1, 1e10, 0);
    RP_MET->Draw(1, 1e10, 0);
    RP_MT_barrel_nume->Draw(1, 1e10, 0);
    RP_MT_endcap_nume->Draw(1, 1e10, 0);
//    RP_MT_barrel_deno->Draw(1, 1e10, 0);
//    RP_MT_endcap_deno->Draw(1, 1e10, 0);
    RP_MT_barrel_ctrl->Draw(1, 1e10, 0);
    RP_MT_endcap_ctrl->Draw(1, 1e10, 0);
    RP_eta->Draw(1, 1e12, 0);
    RP_nVTX->Draw(1, 1e10, 0);
    RP_mass_test->Draw(1, 1e9, 1);

//    cout << "MC PFiso integral: " << ((TH1D*)(s_PFiso_barrel_deno->GetStack()->Last()))->Integral() +
//                                     ((TH1D*)(s_PFiso_endcap_deno->GetStack()->Last()))->Integral() << endl;
//    cout << "Data PFiso integral: " << h_PFiso_barrel_data_deno->Integral() + h_PFiso_endcap_data_deno->Integral() << endl;
//    cout << "--------\nQCD PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "--------\nWJets PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_WJets_Full]->Integral() << endl;
//    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_WJets_Full]->Integral() << endl;
//    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_WJets_Full]->Integral() << endl;
//    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_WJets_Full]->Integral() << endl;
//    cout << "--------\nDY PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_DY_Full]->Integral() << endl;
//    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_DY_Full]->Integral() << endl;
//    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_DY_Full]->Integral() << endl;
//    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_DY_Full]->Integral() << endl;

    cout << "MC integral (deno): " << ((TH1D*)(s_pT_barrel_deno->GetStack()->Last()))->Integral() +
                                     ((TH1D*)(s_pT_endcap_deno->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (deno): " << h_pT_barrel_data_deno->Integral() + h_pT_endcap_data_deno->Integral() << endl;
    cout << "QCD integral (deno): " << h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "MC integral(nume): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() +
            ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume): " << h_pT_barrel_data_nume->Integral() + h_pT_endcap_data_nume->Integral() << endl;
    cout << "QCD integral (nume): " << h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;

    // ---- TEST OF MT CUTS ---- //
    Double_t QCD_full_nume, QCD_full_deno, QCD_full_ctrl, WJets_full_nume, WJets_full_deno, WJets_full_ctrl;
    Double_t QCD_red_nume[250], QCD_red_deno[250], QCD_red_ctrl[250], WJets_red_nume[250], WJets_red_deno[250], WJets_red_ctrl[250], cuts[250];
    Double_t SSB_nume[250], SSB_deno[250], SSB_ctrl[250];
    QCD_full_nume = h_PFiso_dBeta_barrel_MC_nume[_QCDEMEnriched_Full]->Integral() + h_PFiso_dBeta_endcap_MC_nume[_QCDEMEnriched_Full]->Integral();
    QCD_full_deno = h_PFiso_dBeta_barrel_MC_deno[_QCDEMEnriched_Full]->Integral() + h_PFiso_dBeta_endcap_MC_deno[_QCDEMEnriched_Full]->Integral();
    QCD_full_ctrl = h_PFiso_dBeta_barrel_MC_ctrl[_QCDEMEnriched_Full]->Integral() + h_PFiso_dBeta_endcap_MC_ctrl[_QCDEMEnriched_Full]->Integral();
    WJets_full_nume = h_MT_barrel_MC_nume[_WJets_Full]->Integral() + h_MT_endcap_MC_nume[_WJets_Full]->Integral();
    WJets_full_deno = h_MT_barrel_MC_deno[_WJets_Full]->Integral() + h_MT_endcap_MC_deno[_WJets_Full]->Integral();
    WJets_full_ctrl = h_MT_barrel_MC_ctrl[_WJets_Full]->Integral() + h_MT_endcap_MC_ctrl[_WJets_Full]->Integral();
    for (Int_t i_bin=0; i_bin<250; i_bin++)
    {
        cuts[i_bin] = (i_bin + 1) * 2;
        QCD_red_nume[i_bin] = (h_MT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_dBeta_endcap_MC_nume[_QCDEMEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_nume;
        QCD_red_deno[i_bin] = (h_MT_barrel_MC_deno[_QCDEMEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_dBeta_endcap_MC_deno[_QCDEMEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_deno;
        QCD_red_ctrl[i_bin] = (h_MT_barrel_MC_ctrl[_QCDEMEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_dBeta_endcap_MC_ctrl[_QCDEMEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_ctrl;
        WJets_red_nume[i_bin] = (h_MT_barrel_MC_nume[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_dBeta_endcap_MC_nume[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_nume;
        WJets_red_deno[i_bin] = (h_MT_barrel_MC_deno[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_dBeta_endcap_MC_deno[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_deno;
        WJets_red_ctrl[i_bin] = (h_MT_barrel_MC_ctrl[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_dBeta_endcap_MC_ctrl[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_ctrl;

        SSB_nume[i_bin] = QCD_red_nume[i_bin] / (QCD_red_nume[i_bin] + WJets_red_nume[i_bin]);
        SSB_deno[i_bin] = QCD_red_deno[i_bin] / (QCD_red_deno[i_bin] + WJets_red_deno[i_bin]);
        SSB_ctrl[i_bin] = QCD_red_ctrl[i_bin] / (QCD_red_ctrl[i_bin] + WJets_red_ctrl[i_bin]);
    }
    TGraph *g_QCD_cuts_nume = new TGraph(250, cuts, QCD_red_nume);
    TGraph *g_QCD_cuts_deno = new TGraph(250, cuts, QCD_red_deno);
    TGraph *g_QCD_cuts_ctrl = new TGraph(250, cuts, QCD_red_ctrl);
    TGraph *g_WJets_cuts_nume = new TGraph(250, cuts, WJets_red_nume);
    TGraph *g_WJets_cuts_deno = new TGraph(250, cuts, WJets_red_deno);
    TGraph *g_WJets_cuts_ctrl = new TGraph(250, cuts, WJets_red_ctrl);
    TGraph *g_SSB_nume = new TGraph(250, cuts, SSB_nume);
    TGraph *g_SSB_deno = new TGraph(250, cuts, SSB_deno);
    TGraph *g_SSB_ctrl = new TGraph(250, cuts, SSB_ctrl);
    g_QCD_cuts_nume->SetLineWidth(3);
    g_QCD_cuts_deno->SetLineWidth(3);
    g_QCD_cuts_ctrl->SetLineWidth(3);
    g_WJets_cuts_nume->SetLineWidth(3);
    g_WJets_cuts_deno->SetLineWidth(3);
    g_WJets_cuts_ctrl->SetLineWidth(3);
    g_SSB_nume->SetLineWidth(3);
    g_SSB_deno->SetLineWidth(3);
    g_SSB_ctrl->SetLineWidth(3);
    g_WJets_cuts_nume->SetLineColor(kRed);
    g_WJets_cuts_deno->SetLineColor(kRed);
    g_WJets_cuts_ctrl->SetLineColor(kRed);
    g_SSB_nume->SetLineColor(kGreen-2);
    g_SSB_deno->SetLineColor(kRed);
    g_SSB_ctrl->SetLineColor(kBlue);

    TLegend *l_cuts = new TLegend(0.7, 0.8, 0.95, 0.9);
    l_cuts->AddEntry(g_QCD_cuts_nume, "QCD", "l");
    l_cuts->AddEntry(g_WJets_cuts_nume, "W+Jets", "l");
    TLegend *l_SSB = new TLegend(0.65, 0.7, 0.95, 0.87);
    l_SSB->AddEntry(g_SSB_nume, "Numerator", "l");
    l_SSB->AddEntry(g_SSB_deno, "Denominator", "l");
    l_SSB->AddEntry(g_SSB_ctrl, "Non-signal", "l");

    TCanvas *c_cuts_nume = new TCanvas("c_cuts_nume", "Numerator MT cuts", 800, 800);
    c_cuts_nume->SetRightMargin(0.05);
    c_cuts_nume->SetTopMargin(0.07);
    c_cuts_nume->SetLeftMargin(0.13);
    c_cuts_nume->SetBottomMargin(0.13);
    g_QCD_cuts_nume->SetTitle("Numerator");
    g_QCD_cuts_nume->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_nume->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_nume->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_nume->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_nume->Draw();
    g_QCD_cuts_nume->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_nume->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_nume->Draw("same");
    l_cuts->Draw();
    c_cuts_nume->SetGridx();
    c_cuts_nume->SetGridy();
    c_cuts_nume->Update();
    TCanvas *c_cuts_deno = new TCanvas("c_cuts_deno", "Denominator MT cuts", 800, 800);
    c_cuts_deno->SetRightMargin(0.05);
    c_cuts_deno->SetTopMargin(0.07);
    c_cuts_deno->SetLeftMargin(0.13);
    c_cuts_deno->SetBottomMargin(0.13);
    g_QCD_cuts_deno->SetTitle("Denominator");
    g_QCD_cuts_deno->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_deno->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_deno->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_deno->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_deno->Draw();
    g_QCD_cuts_deno->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_deno->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_deno->Draw("same");
    l_cuts->Draw();
    c_cuts_deno->SetGridx();
    c_cuts_deno->SetGridy();
    c_cuts_deno->Update();
    TCanvas *c_cuts_ctrl = new TCanvas("c_cuts_ctrl", "Non-signal MT cuts", 800, 800);
    c_cuts_ctrl->SetRightMargin(0.05);
    c_cuts_ctrl->SetTopMargin(0.07);
    c_cuts_ctrl->SetLeftMargin(0.13);
    c_cuts_ctrl->SetBottomMargin(0.13);
    g_QCD_cuts_ctrl->SetTitle("Non-signal region");
    g_QCD_cuts_ctrl->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_ctrl->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_ctrl->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_ctrl->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_ctrl->Draw();
    g_QCD_cuts_ctrl->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_ctrl->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_ctrl->Draw("same");
    l_cuts->Draw();
    c_cuts_ctrl->SetGridx();
    c_cuts_ctrl->SetGridy();
    c_cuts_ctrl->Update();
    TCanvas *c_SSB = new TCanvas("c_SSB", "S/(S+B) MT cuts", 800, 800);
    c_SSB->SetRightMargin(0.05);
    c_SSB->SetTopMargin(0.13);
    c_SSB->SetLeftMargin(0.13);
    c_SSB->SetBottomMargin(0.13);
    g_SSB_nume->SetTitle("#frac{QCD}{QCD+WJets}");
    g_SSB_nume->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_SSB_nume->GetYaxis()->SetTitle("QCD/(QCD+WJets)");
    g_SSB_nume->GetXaxis()->SetTitleSize(0.05);
    g_SSB_nume->GetYaxis()->SetTitleSize(0.05);
    g_SSB_nume->Draw();
    g_SSB_nume->GetXaxis()->SetRangeUser(0, 250);
    g_SSB_nume->GetYaxis()->SetRangeUser(0, 1);
    g_SSB_deno->Draw("same");
    g_SSB_ctrl->Draw("same");
    l_SSB->Draw();
    c_SSB->SetGridx();
    c_SSB->SetGridy();
    c_SSB->Update();
} // End of EE_HistDrawer()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_HistDrawer(Int_t type)
{
    FileMgr fm;
    THStack *s_PFiso_barrel_nume = new THStack("s_PFiso_barrel_nume", "");
    THStack *s_PFiso_endcap_nume = new THStack("s_PFiso_endcap_nume", "");
    THStack *s_PFiso_barrel_deno = new THStack("s_PFiso_barrel_deno", "");
    THStack *s_PFiso_endcap_deno = new THStack("s_PFiso_endcap_deno", "");
    THStack *s_PFiso_barrel_ctrl = new THStack("s_PFiso_barrel_ctrl", "");
    THStack *s_PFiso_endcap_ctrl = new THStack("s_PFiso_endcap_ctrl", "");
    THStack *s_pT_barrel_nume = new THStack("s_pT_barrel_nume", "");
    THStack *s_pT_endcap_nume = new THStack("s_pT_endcap_nume", "");
    THStack *s_pT_barrel_deno = new THStack("s_pT_barrel_deno", "");
    THStack *s_pT_endcap_deno = new THStack("s_pT_endcap_deno", "");
    THStack *s_pT_barrel_ctrl = new THStack("s_pT_barrel_ctrl", "");
    THStack *s_pT_endcap_ctrl = new THStack("s_pT_endcap_ctrl", "");
    THStack *s_MET = new THStack("s_MET", "");
    THStack *s_MT_barrel_nume = new THStack("s_MT_barrel_nume", "");
    THStack *s_MT_endcap_nume = new THStack("s_MT_endcap_nume", "");
    THStack *s_MT_barrel_deno = new THStack("s_MT_barrel_deno", "");
    THStack *s_MT_endcap_deno = new THStack("s_MT_endcap_deno", "");
    THStack *s_MT_barrel_ctrl = new THStack("s_MT_barrel_ctrl", "");
    THStack *s_MT_endcap_ctrl = new THStack("s_MT_endcap_ctrl", "");
    THStack *s_eta = new THStack("s_eta", "");
    THStack *s_nVTX = new THStack("s_nVTX", "");

    TH1D *h_PFiso_barrel_MC_nume[_EndOf_Data_Special], *h_PFiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_PFiso_barrel_MC_deno[_EndOf_Data_Special], *h_PFiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_PFiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_PFiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_pT_barrel_MC_nume[_EndOf_Data_Special], *h_pT_endcap_MC_nume[_EndOf_Data_Special],
         *h_pT_barrel_MC_deno[_EndOf_Data_Special], *h_pT_endcap_MC_deno[_EndOf_Data_Special],
         *h_pT_barrel_MC_ctrl[_EndOf_Data_Special], *h_pT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_MET_MC[_EndOf_Data_Special], *h_eta_data, *h_nVTX_data,
         *h_MT_barrel_MC_nume[_EndOf_Data_Special], *h_MT_endcap_MC_nume[_EndOf_Data_Special],
         *h_MT_barrel_MC_deno[_EndOf_Data_Special], *h_MT_endcap_MC_deno[_EndOf_Data_Special],
         *h_MT_barrel_MC_ctrl[_EndOf_Data_Special], *h_MT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_eta_MC[_EndOf_Data_Special], *h_nVTX_MC[_EndOf_Data_Special],
         *h_PFiso_barrel_data_nume, *h_PFiso_endcap_data_nume, *h_PFiso_barrel_data_deno, *h_PFiso_endcap_data_deno,
         *h_PFiso_barrel_data_ctrl, *h_PFiso_endcap_data_ctrl, *h_pT_barrel_data_nume, *h_pT_endcap_data_nume,
         *h_pT_barrel_data_deno, *h_pT_endcap_data_deno, *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl,
         *h_MET_data, *h_MT_barrel_data_nume, *h_MT_endcap_data_nume,
         *h_MT_barrel_data_deno, *h_MT_endcap_data_deno, *h_MT_barrel_data_ctrl, *h_MT_endcap_data_ctrl;

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
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr1]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr1]);
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr1]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr1]);
        file->GetObject("h_PFiso_barrel_ctrl", h_PFiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_PFiso_endcap_ctrl", h_PFiso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr1]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr1]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr1]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr1]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_MET", h_MET_MC[pr1]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr1]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr1]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr1]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr1]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr1]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr1]);
        file->GetObject("h_eta_deno", h_eta_MC[pr1]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr1]);

        removeNegativeBins(h_PFiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_PFiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_PFiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_PFiso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_MET_MC[pr1]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr1]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr1]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr1]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr1]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr1]);
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

        h_PFiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_PFiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_PFiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr1]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr1]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr1]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr1]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_MET_MC[pr1]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr1]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr1]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr1]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr1]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_eta_MC[pr1]->SetFillColor(color);
        h_nVTX_MC[pr1]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_PFiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_PFiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr1]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr1]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr1]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr1]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_MET_MC[pr1]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr1]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr1]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr1]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr1]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_eta_MC[pr1]->SetLineColor(color);
        h_nVTX_MC[pr1]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_PFiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_PFiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_PFiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_MET_MC[pr1]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr1]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr1]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr1]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr1]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_eta_MC[pr1]->SetDirectory(0);
        h_nVTX_MC[pr1]->SetDirectory(0);

        if (pr1 == _WJets)
        {
            h_PFiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_barrel_MC_nume[pr1]->Clone("h_PFiso_barrel_nume_WJets")));
            h_PFiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_PFiso_endcap_MC_nume[pr1]->Clone("h_PFiso_endcap_nume_WJets")));
            h_PFiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_barrel_MC_deno[pr1]->Clone("h_PFiso_barrel_deno_WJets")));
            h_PFiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_PFiso_endcap_MC_deno[pr1]->Clone("h_PFiso_endcap_deno_WJets")));
            h_PFiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_barrel_MC_ctrl[pr1]->Clone("h_PFiso_barrel_ctrl_WJets")));
            h_PFiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_PFiso_endcap_MC_ctrl[pr1]->Clone("h_PFiso_endcap_ctrl_WJets")));
            h_pT_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr1]->Clone("h_pT_barrel_nume_WJets")));
            h_pT_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr1]->Clone("h_pT_endcap_nume_WJets")));
            h_pT_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr1]->Clone("h_pT_barrel_deno_WJets")));
            h_pT_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr1]->Clone("h_pT_endcap_deno_WJets")));
            h_pT_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr1]->Clone("h_pT_barrel_ctrl_WJets")));
            h_pT_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr1]->Clone("h_pT_endcap_ctrl_WJets")));
            h_MET_MC[_WJets_Full] = ((TH1D*)(h_MET_MC[pr1]->Clone("h_MET_DY")));
            h_MT_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr1]->Clone("h_MT_barrel_nume_WJets")));
            h_MT_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr1]->Clone("h_MT_endcap_nume_WJets")));
            h_MT_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr1]->Clone("h_MT_barrel_deno_WJets")));
            h_MT_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr1]->Clone("h_MT_endcap_deno_WJets")));
            h_MT_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr1]->Clone("h_MT_barrel_ctrl_WJets")));
            h_MT_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr1]->Clone("h_MT_endcap_ctrl_WJets")));
            h_eta_MC[_WJets_Full] = ((TH1D*)(h_eta_MC[pr1]->Clone("h_eta_deno_WJets")));
            h_nVTX_MC[_WJets_Full] = ((TH1D*)(h_nVTX_MC[pr1]->Clone("h_nVTX_WJets")));
            h_PFiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_MET_MC[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_eta_MC[_WJets_Full]->SetDirectory(0);
            h_nVTX_MC[_WJets_Full]->SetDirectory(0);
        }
        else if (pr1 == _WJets_ext2v5)
        {
            h_PFiso_barrel_MC_nume[_WJets_Full]->Add(h_PFiso_barrel_MC_nume[pr1]);
            h_PFiso_endcap_MC_nume[_WJets_Full]->Add(h_PFiso_endcap_MC_nume[pr1]);
            h_PFiso_barrel_MC_deno[_WJets_Full]->Add(h_PFiso_barrel_MC_deno[pr1]);
            h_PFiso_endcap_MC_deno[_WJets_Full]->Add(h_PFiso_endcap_MC_deno[pr1]);
            h_PFiso_barrel_MC_ctrl[_WJets_Full]->Add(h_PFiso_barrel_MC_ctrl[pr1]);
            h_PFiso_endcap_MC_ctrl[_WJets_Full]->Add(h_PFiso_endcap_MC_ctrl[pr1]);
            h_pT_barrel_MC_nume[_WJets_Full]->Add(h_pT_barrel_MC_nume[pr1]);
            h_pT_endcap_MC_nume[_WJets_Full]->Add(h_pT_endcap_MC_nume[pr1]);
            h_pT_barrel_MC_deno[_WJets_Full]->Add(h_pT_barrel_MC_deno[pr1]);
            h_pT_endcap_MC_deno[_WJets_Full]->Add(h_pT_endcap_MC_deno[pr1]);
            h_pT_barrel_MC_ctrl[_WJets_Full]->Add(h_pT_barrel_MC_ctrl[pr1]);
            h_pT_endcap_MC_ctrl[_WJets_Full]->Add(h_pT_endcap_MC_ctrl[pr1]);
            h_MET_MC[_WJets_Full]->Add(h_MET_MC[pr1]);
            h_MT_barrel_MC_nume[_WJets_Full]->Add(h_MT_barrel_MC_nume[pr1]);
            h_MT_endcap_MC_nume[_WJets_Full]->Add(h_MT_endcap_MC_nume[pr1]);
            h_MT_barrel_MC_deno[_WJets_Full]->Add(h_MT_barrel_MC_deno[pr1]);
            h_MT_endcap_MC_deno[_WJets_Full]->Add(h_MT_endcap_MC_deno[pr1]);
            h_MT_barrel_MC_ctrl[_WJets_Full]->Add(h_MT_barrel_MC_ctrl[pr1]);
            h_MT_endcap_MC_ctrl[_WJets_Full]->Add(h_MT_endcap_MC_ctrl[pr1]);
            h_eta_MC[_WJets_Full]->Add(h_eta_MC[pr1]);
            h_nVTX_MC[_WJets_Full]->Add(h_nVTX_MC[pr1]);
        }

        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr1]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr1]);
        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr1]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr1]);
        s_PFiso_barrel_ctrl->Add(h_PFiso_barrel_MC_ctrl[pr1]);
        s_PFiso_endcap_ctrl->Add(h_PFiso_endcap_MC_ctrl[pr1]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr1]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr1]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr1]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr1]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr1]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr1]);
        s_MET->Add(h_MET_MC[pr1]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr1]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr1]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr1]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr1]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr1]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr1]);
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
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_ctrl", h_PFiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_endcap_ctrl", h_PFiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_MET", h_MET_MC[pr]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_MET_MC[pr]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);

        if (pr == _DY_10to50)
        {
            h_PFiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_barrel_MC_nume[pr]->Clone("h_PFiso_barrel_nume_DY")));
            h_PFiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_PFiso_endcap_MC_nume[pr]->Clone("h_PFiso_endcap_nume_DY")));
            h_PFiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_barrel_MC_deno[pr]->Clone("h_PFiso_barrel_deno_DY")));
            h_PFiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_PFiso_endcap_MC_deno[pr]->Clone("h_PFiso_endcap_deno_DY")));
            h_PFiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_barrel_MC_ctrl[pr]->Clone("h_PFiso_barrel_ctrl_DY")));
            h_PFiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_PFiso_endcap_MC_ctrl[pr]->Clone("h_PFiso_endcap_ctrl_DY")));
            h_pT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_nume_DY")));
            h_pT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_nume_DY")));
            h_pT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_deno_DY")));
            h_pT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_deno_DY")));
            h_pT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_ctrl_DY")));
            h_pT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_ctrl_DY")));
            h_MET_MC[_DY_Full] = ((TH1D*)(h_MET_MC[pr]->Clone("h_MET_DY")));
            h_MT_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr]->Clone("h_MT_barrel_nume_DY")));
            h_MT_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr]->Clone("h_MT_endcap_nume_DY")));
            h_MT_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr]->Clone("h_MT_barrel_deno_DY")));
            h_MT_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr]->Clone("h_MT_endcap_deno_DY")));
            h_MT_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr]->Clone("h_MT_barrel_ctrl_DY")));
            h_MT_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr]->Clone("h_MT_endcap_ctrl_DY")));
            h_eta_MC[_DY_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_DY")));
            h_nVTX_MC[_DY_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_DY")));
            h_PFiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_MET_MC[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_MC[_DY_Full]->SetDirectory(0);
            h_nVTX_MC[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_barrel_MC_nume[_DY_Full]->Add(h_PFiso_barrel_MC_nume[pr]);
            h_PFiso_endcap_MC_nume[_DY_Full]->Add(h_PFiso_endcap_MC_nume[pr]);
            h_PFiso_barrel_MC_deno[_DY_Full]->Add(h_PFiso_barrel_MC_deno[pr]);
            h_PFiso_endcap_MC_deno[_DY_Full]->Add(h_PFiso_endcap_MC_deno[pr]);
            h_PFiso_barrel_MC_ctrl[_DY_Full]->Add(h_PFiso_barrel_MC_ctrl[pr]);
            h_PFiso_endcap_MC_ctrl[_DY_Full]->Add(h_PFiso_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_DY_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_DY_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_MC_deno[_DY_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_DY_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_DY_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_DY_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_MET_MC[_DY_Full]->Add(h_MET_MC[pr]);
            h_MT_barrel_MC_nume[_DY_Full]->Add(h_MT_barrel_MC_nume[pr]);
            h_MT_endcap_MC_nume[_DY_Full]->Add(h_MT_endcap_MC_nume[pr]);
            h_MT_barrel_MC_deno[_DY_Full]->Add(h_MT_barrel_MC_deno[pr]);
            h_MT_endcap_MC_deno[_DY_Full]->Add(h_MT_endcap_MC_deno[pr]);
            h_MT_barrel_MC_ctrl[_DY_Full]->Add(h_MT_barrel_MC_ctrl[pr]);
            h_MT_endcap_MC_ctrl[_DY_Full]->Add(h_MT_endcap_MC_ctrl[pr]);
            h_eta_MC[_DY_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_DY_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kOrange - 5;
        h_PFiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_MET_MC[pr]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_MET_MC[pr]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_MET_MC[pr]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr]);
        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr]);
        s_PFiso_barrel_ctrl->Add(h_PFiso_barrel_MC_ctrl[pr]);
        s_PFiso_endcap_ctrl->Add(h_PFiso_endcap_MC_ctrl[pr]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
        s_MET->Add(h_MET_MC[pr]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
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
        file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_MC_nume[pr]);
        file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_MC_nume[pr]);
        file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_MC_deno[pr]);
        file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_MC_deno[pr]);
        file->GetObject("h_PFiso_barrel_ctrl", h_PFiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_PFiso_endcap_ctrl", h_PFiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_pT_barrel_nume", h_pT_barrel_MC_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_MC_nume[pr]);
        file->GetObject("h_pT_barrel_deno", h_pT_barrel_MC_deno[pr]);
        file->GetObject("h_pT_endcap_deno", h_pT_endcap_MC_deno[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_MC_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_MC_ctrl[pr]);
        file->GetObject("h_MET", h_MET_MC[pr]);
        file->GetObject("h_MT_barrel_nume", h_MT_barrel_MC_nume[pr]);
        file->GetObject("h_MT_endcap_nume", h_MT_endcap_MC_nume[pr]);
        file->GetObject("h_MT_barrel_deno", h_MT_barrel_MC_deno[pr]);
        file->GetObject("h_MT_endcap_deno", h_MT_endcap_MC_deno[pr]);
        file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_MC_ctrl[pr]);
        file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_MC_ctrl[pr]);
        file->GetObject("h_eta_deno", h_eta_MC[pr]);
        file->GetObject("h_nVTX", h_nVTX_MC[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_PFiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_PFiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_pT_barrel_MC_nume[pr]);
        removeNegativeBins(h_pT_endcap_MC_nume[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno[pr]);
        removeNegativeBins(h_pT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_MET_MC[pr]);
        removeNegativeBins(h_MT_barrel_MC_nume[pr]);
        removeNegativeBins(h_MT_endcap_MC_nume[pr]);
        removeNegativeBins(h_MT_barrel_MC_deno[pr]);
        removeNegativeBins(h_MT_endcap_MC_deno[pr]);
        removeNegativeBins(h_MT_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_MT_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_eta_MC[pr]);
        removeNegativeBins(h_nVTX_MC[pr]);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_barrel_MC_nume[pr]->Clone("h_PFiso_barrel_nume_QCD")));
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_endcap_MC_nume[pr]->Clone("h_PFiso_endcap_nume_QCD")));
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_barrel_MC_deno[pr]->Clone("h_PFiso_barrel_deno_QCD")));
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_endcap_MC_deno[pr]->Clone("h_PFiso_endcap_deno_QCD")));
            h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_barrel_MC_ctrl[pr]->Clone("h_PFiso_barrel_ctrl_QCD")));
            h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_PFiso_endcap_MC_ctrl[pr]->Clone("h_PFiso_endcap_ctrl_QCD")));
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_nume[pr]->Clone("h_pT_barrel_nume_QCD")));
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_nume[pr]->Clone("h_pT_endcap_nume_QCD")));
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_deno_QCD")));
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_deno_QCD")));
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_ctrl[pr]->Clone("h_pT_barrel_ctrl_QCD")));
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_ctrl[pr]->Clone("h_pT_endcap_ctrl_QCD")));
            h_MET_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_MET_MC[pr]->Clone("h_MET_QCD")));
            h_MT_barrel_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_nume[pr]->Clone("h_MT_barrel_nume_QCD")));
            h_MT_endcap_MC_nume[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_nume[pr]->Clone("h_MT_endcap_nume_QCD")));
            h_MT_barrel_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_deno[pr]->Clone("h_MT_barrel_deno_QCD")));
            h_MT_endcap_MC_deno[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_deno[pr]->Clone("h_MT_endcap_deno_QCD")));
            h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_barrel_MC_ctrl[pr]->Clone("h_MT_barrel_ctrl_QCD")));
            h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full] = ((TH1D*)(h_MT_endcap_MC_ctrl[pr]->Clone("h_MT_endcap_ctrl_QCD")));
            h_eta_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_eta_MC[pr]->Clone("h_eta_deno_QCD")));
            h_nVTX_MC[_QCDMuEnriched_Full] = ((TH1D*)(h_nVTX_MC[pr]->Clone("h_nVTX_QCD")));
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MET_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_nume[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_deno[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetDirectory(0);
            h_eta_MC[_QCDMuEnriched_Full]->SetDirectory(0);
            h_nVTX_MC[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_PFiso_barrel_MC_nume[pr]);
            h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_PFiso_endcap_MC_nume[pr]);
            h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_PFiso_barrel_MC_deno[pr]);
            h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_PFiso_endcap_MC_deno[pr]);
            h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->Add(h_PFiso_barrel_MC_ctrl[pr]);
            h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->Add(h_PFiso_endcap_MC_ctrl[pr]);
            h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_nume[pr]);
            h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_nume[pr]);
            h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_deno[pr]);
            h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_deno[pr]);
            h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_ctrl[pr]);
            h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_ctrl[pr]);
            h_MET_MC[_QCDMuEnriched_Full]->Add(h_MET_MC[pr]);
            h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->Add(h_MT_barrel_MC_nume[pr]);
            h_MT_endcap_MC_nume[_QCDMuEnriched_Full]->Add(h_MT_endcap_MC_nume[pr]);
            h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->Add(h_MT_barrel_MC_deno[pr]);
            h_MT_endcap_MC_deno[_QCDMuEnriched_Full]->Add(h_MT_endcap_MC_deno[pr]);
            h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->Add(h_MT_barrel_MC_ctrl[pr]);
            h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]->Add(h_MT_endcap_MC_ctrl[pr]);
            h_eta_MC[_QCDMuEnriched_Full]->Add(h_eta_MC[pr]);
            h_nVTX_MC[_QCDMuEnriched_Full]->Add(h_nVTX_MC[pr]);
        }

        Color_t color = kRed + 3;
        h_PFiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_PFiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_pT_barrel_MC_nume[pr]->SetFillColor(color);
        h_pT_endcap_MC_nume[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno[pr]->SetFillColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_MET_MC[pr]->SetFillColor(color);
        h_MT_barrel_MC_nume[pr]->SetFillColor(color);
        h_MT_endcap_MC_nume[pr]->SetFillColor(color);
        h_MT_barrel_MC_deno[pr]->SetFillColor(color);
        h_MT_endcap_MC_deno[pr]->SetFillColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_eta_MC[pr]->SetFillColor(color);
        h_nVTX_MC[pr]->SetFillColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_PFiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_pT_barrel_MC_nume[pr]->SetLineColor(color);
        h_pT_endcap_MC_nume[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno[pr]->SetLineColor(color);
        h_pT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_MET_MC[pr]->SetLineColor(color);
        h_MT_barrel_MC_nume[pr]->SetLineColor(color);
        h_MT_endcap_MC_nume[pr]->SetLineColor(color);
        h_MT_barrel_MC_deno[pr]->SetLineColor(color);
        h_MT_endcap_MC_deno[pr]->SetLineColor(color);
        h_MT_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_MT_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_eta_MC[pr]->SetLineColor(color);
        h_nVTX_MC[pr]->SetLineColor(color);
        h_PFiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_PFiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_PFiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_pT_barrel_MC_nume[pr]->SetDirectory(0);
        h_pT_endcap_MC_nume[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno[pr]->SetDirectory(0);
        h_pT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_MET_MC[pr]->SetDirectory(0);
        h_MT_barrel_MC_nume[pr]->SetDirectory(0);
        h_MT_endcap_MC_nume[pr]->SetDirectory(0);
        h_MT_barrel_MC_deno[pr]->SetDirectory(0);
        h_MT_endcap_MC_deno[pr]->SetDirectory(0);
        h_MT_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_MT_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_eta_MC[pr]->SetDirectory(0);
        h_nVTX_MC[pr]->SetDirectory(0);

        s_PFiso_barrel_nume->Add(h_PFiso_barrel_MC_nume[pr]);
        s_PFiso_endcap_nume->Add(h_PFiso_endcap_MC_nume[pr]);
        s_PFiso_barrel_deno->Add(h_PFiso_barrel_MC_deno[pr]);
        s_PFiso_endcap_deno->Add(h_PFiso_endcap_MC_deno[pr]);
        s_PFiso_barrel_ctrl->Add(h_PFiso_barrel_MC_ctrl[pr]);
        s_PFiso_endcap_ctrl->Add(h_PFiso_endcap_MC_ctrl[pr]);
        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
        s_MET->Add(h_MET_MC[pr]);
        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
        s_eta->Add(h_eta_MC[pr]);
        s_nVTX->Add(h_nVTX_MC[pr]);

        file->Close();
    }

     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
                   h_MET_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
        h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillColor(kGray);
                   h_eta_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetFillColor(kGray);
     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
                   h_MET_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_endcap_MC_nume[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_endcap_MC_deno[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
        h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetFillStyle(3002);
                   h_eta_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetFillStyle(3002);
     h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
     h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_barrel_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_endcap_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_barrel_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_endcap_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
                   h_MET_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_endcap_MC_nume[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_endcap_MC_deno[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
        h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]->SetLineColor(kRed);
                   h_eta_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);
                  h_nVTX_MC[_QCDMuEnriched_Full]->SetLineColor(kRed);

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        TFile *file;
        if (type == 1) file = new TFile("/media/sf_DATA/FR/Muon/SelectedForFR_Mu_"+fm.Procname[pr]+".root", "READ");
        else if (type == 2) file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        else return;
        TH1D *h_temp[21];
        if (pr == _SingleMuon_B)
        {
            file->GetObject("h_PFiso_barrel_nume", h_PFiso_barrel_data_nume);
            file->GetObject("h_PFiso_endcap_nume", h_PFiso_endcap_data_nume);
            file->GetObject("h_PFiso_barrel_deno", h_PFiso_barrel_data_deno);
            file->GetObject("h_PFiso_endcap_deno", h_PFiso_endcap_data_deno);
            file->GetObject("h_PFiso_barrel_ctrl", h_PFiso_barrel_data_ctrl);
            file->GetObject("h_PFiso_endcap_ctrl", h_PFiso_endcap_data_ctrl);
            file->GetObject("h_pT_barrel_nume", h_pT_barrel_data_nume);
            file->GetObject("h_pT_endcap_nume", h_pT_endcap_data_nume);
            file->GetObject("h_pT_barrel_deno", h_pT_barrel_data_deno);
            file->GetObject("h_pT_endcap_deno", h_pT_endcap_data_deno);
            file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_data_ctrl);
            file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_data_ctrl);
            file->GetObject("h_MET", h_MET_data);
            file->GetObject("h_MT_barrel_nume", h_MT_barrel_data_nume);
            file->GetObject("h_MT_endcap_nume", h_MT_endcap_data_nume);
            file->GetObject("h_MT_barrel_deno", h_MT_barrel_data_deno);
            file->GetObject("h_MT_endcap_deno", h_MT_endcap_data_deno);
            file->GetObject("h_MT_barrel_ctrl", h_MT_barrel_data_ctrl);
            file->GetObject("h_MT_endcap_ctrl", h_MT_endcap_data_ctrl);
            file->GetObject("h_eta_deno", h_eta_data);
            file->GetObject("h_nVTX", h_nVTX_data);
            removeNegativeBins(h_PFiso_barrel_data_nume);
            removeNegativeBins(h_PFiso_endcap_data_nume);
            removeNegativeBins(h_PFiso_barrel_data_deno);
            removeNegativeBins(h_PFiso_endcap_data_deno);
            removeNegativeBins(h_PFiso_barrel_data_ctrl);
            removeNegativeBins(h_PFiso_endcap_data_ctrl);
            removeNegativeBins(h_pT_barrel_data_nume);
            removeNegativeBins(h_pT_endcap_data_nume);
            removeNegativeBins(h_pT_barrel_data_deno);
            removeNegativeBins(h_pT_endcap_data_deno);
            removeNegativeBins(h_pT_barrel_data_ctrl);
            removeNegativeBins(h_pT_endcap_data_ctrl);
            removeNegativeBins(h_MET_data);
            removeNegativeBins(h_MT_barrel_data_nume);
            removeNegativeBins(h_MT_endcap_data_nume);
            removeNegativeBins(h_MT_barrel_data_deno);
            removeNegativeBins(h_MT_endcap_data_deno);
            removeNegativeBins(h_MT_barrel_data_ctrl);
            removeNegativeBins(h_MT_endcap_data_ctrl);
            removeNegativeBins(h_eta_data);
            removeNegativeBins(h_nVTX_data);
        }
        else
        {
            file->GetObject("h_PFiso_barrel_nume", h_temp[0]);
            file->GetObject("h_PFiso_endcap_nume", h_temp[1]);
            file->GetObject("h_PFiso_barrel_deno", h_temp[2]);
            file->GetObject("h_PFiso_endcap_deno", h_temp[3]);
            file->GetObject("h_PFiso_barrel_ctrl", h_temp[4]);
            file->GetObject("h_PFiso_endcap_ctrl", h_temp[5]);
            file->GetObject("h_pT_barrel_nume", h_temp[6]);
            file->GetObject("h_pT_endcap_nume", h_temp[7]);
            file->GetObject("h_pT_barrel_deno", h_temp[8]);
            file->GetObject("h_pT_endcap_deno", h_temp[9]);
            file->GetObject("h_pT_barrel_ctrl", h_temp[10]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[11]);
            file->GetObject("h_MET", h_temp[12]);
            file->GetObject("h_MT_barrel_nume", h_temp[13]);
            file->GetObject("h_MT_endcap_nume", h_temp[14]);
            file->GetObject("h_MT_barrel_deno", h_temp[15]);
            file->GetObject("h_MT_endcap_deno", h_temp[16]);
            file->GetObject("h_MT_barrel_ctrl", h_temp[17]);
            file->GetObject("h_MT_endcap_ctrl", h_temp[18]);
            file->GetObject("h_eta_deno", h_temp[19]);
            file->GetObject("h_nVTX", h_temp[20]);
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
            removeNegativeBins(h_temp[16]);
            removeNegativeBins(h_temp[17]);
            removeNegativeBins(h_temp[18]);
            removeNegativeBins(h_temp[19]);
            removeNegativeBins(h_temp[20]);
            h_PFiso_barrel_data_nume->Add(h_temp[0]);
            h_PFiso_endcap_data_nume->Add(h_temp[1]);
            h_PFiso_barrel_data_deno->Add(h_temp[2]);
            h_PFiso_endcap_data_deno->Add(h_temp[3]);
            h_PFiso_barrel_data_ctrl->Add(h_temp[4]);
            h_PFiso_endcap_data_ctrl->Add(h_temp[5]);
            h_pT_barrel_data_nume->Add(h_temp[6]);
            h_pT_endcap_data_nume->Add(h_temp[7]);
            h_pT_barrel_data_deno->Add(h_temp[8]);
            h_pT_endcap_data_deno->Add(h_temp[9]);
            h_pT_barrel_data_ctrl->Add(h_temp[10]);
            h_pT_endcap_data_ctrl->Add(h_temp[11]);
            h_MET_data->Add(h_temp[12]);
            h_MT_barrel_data_nume->Add(h_temp[13]);
            h_MT_endcap_data_nume->Add(h_temp[14]);
            h_MT_barrel_data_deno->Add(h_temp[15]);
            h_MT_endcap_data_deno->Add(h_temp[16]);
            h_MT_barrel_data_ctrl->Add(h_temp[17]);
            h_MT_endcap_data_ctrl->Add(h_temp[18]);
            h_eta_data->Add(h_temp[19]);
            h_nVTX_data->Add(h_temp[20]);
        }
    }

    h_PFiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_PFiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_PFiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_PFiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_MET_data->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_MT_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_MT_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_eta_data->SetMarkerStyle(kFullDotLarge);
    h_nVTX_data->SetMarkerStyle(kFullDotLarge);
    h_PFiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_PFiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_PFiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_PFiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_PFiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_PFiso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_pT_barrel_data_nume->SetMarkerColor(kBlack);
    h_pT_endcap_data_nume->SetMarkerColor(kBlack);
    h_pT_barrel_data_deno->SetMarkerColor(kBlack);
    h_pT_endcap_data_deno->SetMarkerColor(kBlack);
    h_pT_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_pT_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_MET_data->SetMarkerColor(kBlack);
    h_MT_barrel_data_nume->SetMarkerColor(kBlack);
    h_MT_endcap_data_nume->SetMarkerColor(kBlack);
    h_MT_barrel_data_deno->SetMarkerColor(kBlack);
    h_MT_endcap_data_deno->SetMarkerColor(kBlack);
    h_MT_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_MT_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_eta_data->SetMarkerColor(kBlack);
    h_nVTX_data->SetMarkerColor(kBlack);
    h_PFiso_barrel_data_nume->SetLineColor(kBlack);
    h_PFiso_endcap_data_nume->SetLineColor(kBlack);
    h_PFiso_barrel_data_deno->SetLineColor(kBlack);
    h_PFiso_endcap_data_deno->SetLineColor(kBlack);
    h_PFiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_PFiso_endcap_data_ctrl->SetLineColor(kBlack);
    h_pT_barrel_data_nume->SetLineColor(kBlack);
    h_pT_endcap_data_nume->SetLineColor(kBlack);
    h_pT_barrel_data_deno->SetLineColor(kBlack);
    h_pT_endcap_data_deno->SetLineColor(kBlack);
    h_pT_barrel_data_ctrl->SetLineColor(kBlack);
    h_pT_endcap_data_ctrl->SetLineColor(kBlack);
    h_MET_data->SetLineColor(kBlack);
    h_MT_barrel_data_nume->SetLineColor(kBlack);
    h_MT_endcap_data_nume->SetLineColor(kBlack);
    h_MT_barrel_data_deno->SetLineColor(kBlack);
    h_MT_endcap_data_deno->SetLineColor(kBlack);
    h_MT_barrel_data_ctrl->SetLineColor(kBlack);
    h_MT_endcap_data_ctrl->SetLineColor(kBlack);
    h_eta_data->SetLineColor(kBlack);
    h_nVTX_data->SetLineColor(kBlack);
    h_PFiso_barrel_data_nume->SetDirectory(0);
    h_PFiso_endcap_data_nume->SetDirectory(0);
    h_PFiso_barrel_data_deno->SetDirectory(0);
    h_PFiso_endcap_data_deno->SetDirectory(0);
    h_PFiso_barrel_data_ctrl->SetDirectory(0);
    h_PFiso_endcap_data_ctrl->SetDirectory(0);
    h_pT_barrel_data_nume->SetDirectory(0);
    h_pT_endcap_data_nume->SetDirectory(0);
    h_pT_barrel_data_deno->SetDirectory(0);
    h_pT_endcap_data_deno->SetDirectory(0);
    h_pT_barrel_data_ctrl->SetDirectory(0);
    h_pT_endcap_data_ctrl->SetDirectory(0);
    h_MET_data->SetDirectory(0);
    h_MT_barrel_data_nume->SetDirectory(0);
    h_MT_endcap_data_nume->SetDirectory(0);
    h_MT_barrel_data_deno->SetDirectory(0);
    h_MT_endcap_data_deno->SetDirectory(0);
    h_MT_barrel_data_ctrl->SetDirectory(0);
    h_MT_endcap_data_ctrl->SetDirectory(0);
    h_eta_data->SetDirectory(0);
    h_nVTX_data->SetDirectory(0);

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_PFiso_barrel_nume = new myRatioPlot_t("RP_PFiso_barrel_nume", s_PFiso_barrel_nume, h_PFiso_barrel_data_nume);
    myRatioPlot_t *RP_PFiso_endcap_nume = new myRatioPlot_t("RP_PFiso_endcap_nume", s_PFiso_endcap_nume, h_PFiso_endcap_data_nume);
    myRatioPlot_t *RP_PFiso_barrel_deno = new myRatioPlot_t("RP_PFiso_barrel_deno", s_PFiso_barrel_deno, h_PFiso_barrel_data_deno);
    myRatioPlot_t *RP_PFiso_endcap_deno = new myRatioPlot_t("RP_PFiso_endcap_deno", s_PFiso_endcap_deno, h_PFiso_endcap_data_deno);
    myRatioPlot_t *RP_PFiso_barrel_ctrl = new myRatioPlot_t("RP_PFiso_barrel_ctrl", s_PFiso_barrel_ctrl, h_PFiso_barrel_data_ctrl);
    myRatioPlot_t *RP_PFiso_endcap_ctrl = new myRatioPlot_t("RP_PFiso_endcap_ctrl", s_PFiso_endcap_ctrl, h_PFiso_endcap_data_ctrl);
    myRatioPlot_t *RP_pT_barrel_nume = new myRatioPlot_t("RP_pT_barrel_nume", s_pT_barrel_nume, h_pT_barrel_data_nume);
    myRatioPlot_t *RP_pT_endcap_nume = new myRatioPlot_t("RP_pT_endcap_nume", s_pT_endcap_nume, h_pT_endcap_data_nume);
    myRatioPlot_t *RP_pT_barrel_deno = new myRatioPlot_t("RP_pT_barrel_deno", s_pT_barrel_deno, h_pT_barrel_data_deno);
    myRatioPlot_t *RP_pT_endcap_deno = new myRatioPlot_t("RP_pT_endcap_deno", s_pT_endcap_deno, h_pT_endcap_data_deno);
    myRatioPlot_t *RP_pT_barrel_ctrl = new myRatioPlot_t("RP_pT_barrel_ctrl", s_pT_barrel_ctrl, h_pT_barrel_data_ctrl);
    myRatioPlot_t *RP_pT_endcap_ctrl = new myRatioPlot_t("RP_pT_endcap_ctrl", s_pT_endcap_ctrl, h_pT_endcap_data_ctrl);
    myRatioPlot_t *RP_MET = new myRatioPlot_t("RP_MET", s_MET, h_MET_data);
    myRatioPlot_t *RP_MT_barrel_nume = new myRatioPlot_t("RP_MT_barrel_nume", s_MT_barrel_nume, h_MT_barrel_data_nume);
    myRatioPlot_t *RP_MT_endcap_nume = new myRatioPlot_t("RP_MT_endcap_nume", s_MT_endcap_nume, h_MT_endcap_data_nume);
    myRatioPlot_t *RP_MT_barrel_deno = new myRatioPlot_t("RP_MT_barrel_deno", s_MT_barrel_deno, h_MT_barrel_data_deno);
    myRatioPlot_t *RP_MT_endcap_deno = new myRatioPlot_t("RP_MT_endcap_deno", s_MT_endcap_deno, h_MT_endcap_data_deno);
    myRatioPlot_t *RP_MT_barrel_ctrl = new myRatioPlot_t("RP_MT_barrel_ctrl", s_MT_barrel_ctrl, h_MT_barrel_data_ctrl);
    myRatioPlot_t *RP_MT_endcap_ctrl = new myRatioPlot_t("RP_MT_endcap_ctrl", s_MT_endcap_ctrl, h_MT_endcap_data_ctrl);
    myRatioPlot_t *RP_eta = new myRatioPlot_t("RP_eta", s_eta, h_eta_data);
    myRatioPlot_t *RP_nVTX = new myRatioPlot_t("RP_nVTX", s_nVTX, h_nVTX_data);

    RP_PFiso_barrel_nume->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{nume})", 0, 0.15);
    RP_PFiso_endcap_nume->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{nume})", 0, 0.15);
    RP_PFiso_barrel_deno->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_PFiso_endcap_deno->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_PFiso_barrel_ctrl->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{barrel}}^{control})", 0.15, 5);
    RP_PFiso_endcap_ctrl->SetPlots("I_{PF}^{rel.} (#mu_{#lower[-0.4]{endcap}}^{control})", 0.15, 5);
    RP_pT_barrel_nume->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{nume}) [GeV/c]", 52, 1000);
    RP_pT_endcap_nume->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{nume}) [GeV/c]", 52, 1000);
    RP_pT_barrel_deno->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{deno}) [GeV/c]", 52, 1000);
    RP_pT_endcap_deno->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{deno}) [GeV/c]", 52, 1000);
    RP_pT_barrel_ctrl->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{control}) [GeV/c]", 52, 1000);
    RP_pT_endcap_ctrl->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{control}) [GeV/c]", 52, 1000);
    RP_MET->SetPlots("E_{#lower[-0.25]{T}}^{miss} [GeV]", 0, 500);
    RP_MT_barrel_nume->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{nume}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_nume->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{nume}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_barrel_deno->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{deno}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_deno->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{deno}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_barrel_ctrl->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{control}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_MT_endcap_ctrl->SetPlots("m_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{control}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
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

    RP_PFiso_barrel_nume->ImportLegend(legend);
    RP_PFiso_endcap_nume->ImportLegend(legend);
    RP_PFiso_barrel_deno->ImportLegend(legend);
    RP_PFiso_endcap_deno->ImportLegend(legend);
    RP_PFiso_barrel_ctrl->ImportLegend(legend);
    RP_PFiso_endcap_ctrl->ImportLegend(legend);
    RP_pT_barrel_nume->ImportLegend(legend);
    RP_pT_endcap_nume->ImportLegend(legend);
    RP_pT_barrel_deno->ImportLegend(legend);
    RP_pT_endcap_deno->ImportLegend(legend);
    RP_pT_barrel_ctrl->ImportLegend(legend);
    RP_pT_endcap_ctrl->ImportLegend(legend);
    RP_MET->ImportLegend(legend);
    RP_MT_barrel_nume->ImportLegend(legend);
    RP_MT_endcap_nume->ImportLegend(legend);
    RP_MT_barrel_deno->ImportLegend(legend);
    RP_MT_endcap_deno->ImportLegend(legend);
    RP_MT_barrel_ctrl->ImportLegend(legend);
    RP_MT_endcap_ctrl->ImportLegend(legend);
    RP_eta->ImportLegend(legend);
    RP_nVTX->ImportLegend(legend);

    RP_PFiso_barrel_nume->Draw(1, 1e8, 0);
//    RP_PFiso_barrel_nume->DrawOnTop(h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]);

    RP_PFiso_endcap_nume->Draw(1, 1e8, 0);
//    RP_PFiso_endcap_nume->DrawOnTop(h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]);

    RP_PFiso_barrel_deno->Draw(1, 1e8, 0);
//    RP_PFiso_barrel_deno->DrawOnTop(h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]);

    RP_PFiso_endcap_deno->Draw(1, 1e8, 0);
//    RP_PFiso_endcap_deno->DrawOnTop(h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]);

    RP_PFiso_barrel_ctrl->Draw(1, 1e8, 0);
//    RP_PFiso_barrel_ctrl->DrawOnTop(h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]);

    RP_PFiso_endcap_ctrl->Draw(1, 1e8, 0);
//    RP_PFiso_endcap_ctrl->DrawOnTop(h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]);

    RP_pT_barrel_nume->Draw(1, 1e8, 0);
//    RP_pT_barrel_nume->DrawOnTop(h_pT_barrel_MC_nume[_QCDMuEnriched_Full]);

    RP_pT_endcap_nume->Draw(1, 1e8, 0);
//    RP_pT_endcap_nume->DrawOnTop(h_pT_endcap_MC_nume[_QCDMuEnriched_Full]);

    RP_pT_barrel_deno->Draw(1, 1e8, 0);
//    RP_pT_barrel_deno->DrawOnTop(h_pT_barrel_MC_deno[_QCDMuEnriched_Full]);

    RP_pT_endcap_deno->Draw(1, 1e8, 0);
//    RP_pT_endcap_deno->DrawOnTop(h_pT_endcap_MC_deno[_QCDMuEnriched_Full]);

    RP_pT_barrel_ctrl->Draw(1, 1e8, 0);
//    RP_pT_barrel_ctrl->DrawOnTop(h_pT_barrel_MC_ctrl[_QCDMuEnriched_Full]);

    RP_pT_endcap_ctrl->Draw(1, 1e8, 0);
//    RP_pT_endcap_ctrl->DrawOnTop(h_pT_endcap_MC_ctrl[_QCDMuEnriched_Full]);

    RP_MET->Draw(1, 1e8, 0);
//    RP_MET->DrawOnTop(h_MET_MC[_QCDMuEnriched_Full]);

    RP_MT_barrel_nume->Draw(1, 1e8, 0);
//    RP_MT_barrel_nume->DrawOnTop(h_MT_barrel_MC_nume[_QCDMuEnriched_Full]);

    RP_MT_endcap_nume->Draw(1, 1e8, 0);
//    RP_MT_endcap_nume->DrawOnTop(h_MT_endcap_MC_nume[_QCDMuEnriched_Full]);

    RP_MT_barrel_deno->Draw(1, 1e8, 0);
//    RP_MT_barrel_deno->DrawOnTop(h_MT_barrel_MC_deno[_QCDMuEnriched_Full]);

    RP_MT_endcap_deno->Draw(1, 1e8, 0);
//    RP_MT_endcap_deno->DrawOnTop(h_MT_endcap_MC_deno[_QCDMuEnriched_Full]);

    RP_MT_barrel_ctrl->Draw(1, 1e8, 0);
//    RP_MT_barrel_ctrl->DrawOnTop(h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]);

    RP_MT_endcap_ctrl->Draw(1, 1e8, 0);
//    RP_MT_endcap_ctrl->DrawOnTop(h_MT_endcap_MC_ctrl[_QCDMuEnriched_Full]);

    RP_eta->Draw(1, 1e8, 0);
    RP_nVTX->Draw(1, 1e8, 0);

    cout << "MC PFiso integral: " << ((TH1D*)(s_PFiso_barrel_deno->GetStack()->Last()))->Integral() +
                                     ((TH1D*)(s_PFiso_endcap_deno->GetStack()->Last()))->Integral() << endl;
    cout << "Data PFiso integral: " << h_PFiso_barrel_data_deno->Integral() + h_PFiso_endcap_data_deno->Integral() << endl;
    cout << "--------\nQCD PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Integral() << endl;
    cout << "--------\nWJets PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_WJets_Full]->Integral() << endl;
    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_WJets_Full]->Integral() << endl;
    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_WJets_Full]->Integral() << endl;
    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_WJets_Full]->Integral() << endl;
    cout << "--------\nDY PFiso integral:\nBarrel nume   " << h_PFiso_barrel_MC_nume[_DY_Full]->Integral() << endl;
    cout << "Barrel deno   " << h_PFiso_barrel_MC_deno[_DY_Full]->Integral() << endl;
    cout << "Endcap nume   " << h_PFiso_endcap_MC_nume[_DY_Full]->Integral() << endl;
    cout << "Endcap deno   " << h_PFiso_endcap_MC_deno[_DY_Full]->Integral() << endl;

    // ---- TEST OF MT CUTS ---- //
    Double_t QCD_full_nume, QCD_full_deno, QCD_full_ctrl, WJets_full_nume, WJets_full_deno, WJets_full_ctrl;
    Double_t QCD_red_nume[250], QCD_red_deno[250], QCD_red_ctrl[250], WJets_red_nume[250], WJets_red_deno[250], WJets_red_ctrl[250], cuts[250];
    Double_t SSB_nume[250], SSB_deno[250], SSB_ctrl[250];
    QCD_full_nume = h_PFiso_barrel_MC_nume[_QCDMuEnriched_Full]->Integral() + h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Integral();
    QCD_full_deno = h_PFiso_barrel_MC_deno[_QCDMuEnriched_Full]->Integral() + h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Integral();
    QCD_full_ctrl = h_PFiso_barrel_MC_ctrl[_QCDMuEnriched_Full]->Integral() + h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->Integral();
    WJets_full_nume = h_MT_barrel_MC_nume[_WJets_Full]->Integral() + h_MT_endcap_MC_nume[_WJets_Full]->Integral();
    WJets_full_deno = h_MT_barrel_MC_deno[_WJets_Full]->Integral() + h_MT_endcap_MC_deno[_WJets_Full]->Integral();
    WJets_full_ctrl = h_MT_barrel_MC_ctrl[_WJets_Full]->Integral() + h_MT_endcap_MC_ctrl[_WJets_Full]->Integral();
    for (Int_t i_bin=0; i_bin<250; i_bin++)
    {
        cuts[i_bin] = (i_bin + 1) * 2;
        QCD_red_nume[i_bin] = (h_MT_barrel_MC_nume[_QCDMuEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_endcap_MC_nume[_QCDMuEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_nume;
        QCD_red_deno[i_bin] = (h_MT_barrel_MC_deno[_QCDMuEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_endcap_MC_deno[_QCDMuEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_deno;
        QCD_red_ctrl[i_bin] = (h_MT_barrel_MC_ctrl[_QCDMuEnriched_Full]->Integral(1, i_bin+1) +
                               h_PFiso_endcap_MC_ctrl[_QCDMuEnriched_Full]->Integral(1, i_bin+1)) / QCD_full_ctrl;
        WJets_red_nume[i_bin] = (h_MT_barrel_MC_nume[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_endcap_MC_nume[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_nume;
        WJets_red_deno[i_bin] = (h_MT_barrel_MC_deno[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_endcap_MC_deno[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_deno;
        WJets_red_ctrl[i_bin] = (h_MT_barrel_MC_ctrl[_WJets_Full]->Integral(1, i_bin+1) +
                                 h_PFiso_endcap_MC_ctrl[_WJets_Full]->Integral(1, i_bin+1)) / WJets_full_ctrl;

        SSB_nume[i_bin] = QCD_red_nume[i_bin] / (QCD_red_nume[i_bin] + WJets_red_nume[i_bin]);
        SSB_deno[i_bin] = QCD_red_deno[i_bin] / (QCD_red_deno[i_bin] + WJets_red_deno[i_bin]);
        SSB_ctrl[i_bin] = QCD_red_ctrl[i_bin] / (QCD_red_ctrl[i_bin] + WJets_red_ctrl[i_bin]);
    }
    TGraph *g_QCD_cuts_nume = new TGraph(250, cuts, QCD_red_nume);
    TGraph *g_QCD_cuts_deno = new TGraph(250, cuts, QCD_red_deno);
    TGraph *g_QCD_cuts_ctrl = new TGraph(250, cuts, QCD_red_ctrl);
    TGraph *g_WJets_cuts_nume = new TGraph(250, cuts, WJets_red_nume);
    TGraph *g_WJets_cuts_deno = new TGraph(250, cuts, WJets_red_deno);
    TGraph *g_WJets_cuts_ctrl = new TGraph(250, cuts, WJets_red_ctrl);
    TGraph *g_SSB_nume = new TGraph(250, cuts, SSB_nume);
    TGraph *g_SSB_deno = new TGraph(250, cuts, SSB_deno);
    TGraph *g_SSB_ctrl = new TGraph(250, cuts, SSB_ctrl);
    g_QCD_cuts_nume->SetLineWidth(3);
    g_QCD_cuts_deno->SetLineWidth(3);
    g_QCD_cuts_ctrl->SetLineWidth(3);
    g_WJets_cuts_nume->SetLineWidth(3);
    g_WJets_cuts_deno->SetLineWidth(3);
    g_WJets_cuts_ctrl->SetLineWidth(3);
    g_SSB_nume->SetLineWidth(3);
    g_SSB_deno->SetLineWidth(3);
    g_SSB_ctrl->SetLineWidth(3);
    g_WJets_cuts_nume->SetLineColor(kRed);
    g_WJets_cuts_deno->SetLineColor(kRed);
    g_WJets_cuts_ctrl->SetLineColor(kRed);
    g_SSB_nume->SetLineColor(kGreen-2);
    g_SSB_deno->SetLineColor(kRed);
    g_SSB_ctrl->SetLineColor(kBlue);

    TLegend *l_cuts = new TLegend(0.7, 0.8, 0.95, 0.9);
    l_cuts->AddEntry(g_QCD_cuts_nume, "QCD", "l");
    l_cuts->AddEntry(g_WJets_cuts_nume, "W+Jets", "l");
    TLegend *l_SSB = new TLegend(0.65, 0.7, 0.95, 0.87);
    l_SSB->AddEntry(g_SSB_nume, "Numerator", "l");
    l_SSB->AddEntry(g_SSB_deno, "Denominator", "l");
    l_SSB->AddEntry(g_SSB_ctrl, "Non-signal", "l");

    TCanvas *c_cuts_nume = new TCanvas("c_cuts_nume", "Numerator MT cuts", 800, 800);
    c_cuts_nume->SetRightMargin(0.05);
    c_cuts_nume->SetTopMargin(0.07);
    c_cuts_nume->SetLeftMargin(0.13);
    c_cuts_nume->SetBottomMargin(0.13);
    g_QCD_cuts_nume->SetTitle("Numerator");
    g_QCD_cuts_nume->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_nume->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_nume->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_nume->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_nume->Draw();
    g_QCD_cuts_nume->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_nume->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_nume->Draw("same");
    l_cuts->Draw();
    c_cuts_nume->SetGridx();
    c_cuts_nume->SetGridy();
    c_cuts_nume->Update();
    TCanvas *c_cuts_deno = new TCanvas("c_cuts_deno", "Denominator MT cuts", 800, 800);
    c_cuts_deno->SetRightMargin(0.05);
    c_cuts_deno->SetTopMargin(0.07);
    c_cuts_deno->SetLeftMargin(0.13);
    c_cuts_deno->SetBottomMargin(0.13);
    g_QCD_cuts_deno->SetTitle("Denominator");
    g_QCD_cuts_deno->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_deno->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_deno->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_deno->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_deno->Draw();
    g_QCD_cuts_deno->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_deno->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_deno->Draw("same");
    l_cuts->Draw();
    c_cuts_deno->SetGridx();
    c_cuts_deno->SetGridy();
    c_cuts_deno->Update();
    TCanvas *c_cuts_ctrl = new TCanvas("c_cuts_ctrl", "Non-signal MT cuts", 800, 800);
    c_cuts_ctrl->SetRightMargin(0.05);
    c_cuts_ctrl->SetTopMargin(0.07);
    c_cuts_ctrl->SetLeftMargin(0.13);
    c_cuts_ctrl->SetBottomMargin(0.13);
    g_QCD_cuts_ctrl->SetTitle("Non-signal region");
    g_QCD_cuts_ctrl->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_QCD_cuts_ctrl->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts_ctrl->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts_ctrl->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts_ctrl->Draw();
    g_QCD_cuts_ctrl->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts_ctrl->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts_ctrl->Draw("same");
    l_cuts->Draw();
    c_cuts_ctrl->SetGridx();
    c_cuts_ctrl->SetGridy();
    c_cuts_ctrl->Update();
    TCanvas *c_SSB = new TCanvas("c_SSB", "S/(S+B) MT cuts", 800, 800);
    c_SSB->SetRightMargin(0.05);
    c_SSB->SetTopMargin(0.13);
    c_SSB->SetLeftMargin(0.13);
    c_SSB->SetBottomMargin(0.13);
    g_SSB_nume->SetTitle("#frac{QCD}{QCD+WJets}");
    g_SSB_nume->GetXaxis()->SetTitle("M_{T} cut [GeV/c^{2}]");
    g_SSB_nume->GetYaxis()->SetTitle("QCD/(QCD+WJets)");
    g_SSB_nume->GetXaxis()->SetTitleSize(0.05);
    g_SSB_nume->GetYaxis()->SetTitleSize(0.05);
    g_SSB_nume->Draw();
    g_SSB_nume->GetXaxis()->SetRangeUser(0, 250);
    g_SSB_nume->GetYaxis()->SetRangeUser(0, 1);
    g_SSB_deno->Draw("same");
    g_SSB_ctrl->Draw("same");
    l_SSB->Draw();
    c_SSB->SetGridx();
    c_SSB->SetGridy();
    c_SSB->Update();

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
    h_barrel_MC_deno_70to100[_DY_Full]            ->Scale(2.8544e+05 / h_barrel_MC_deno_70to100[_DY_Full]            ->Integral());
    h_barrel_MC_deno_70to100[_WW]                 ->Scale(4.5226e+04 / h_barrel_MC_deno_70to100[_WW]                 ->Integral());
    h_barrel_MC_deno_70to100[_WZ]                 ->Scale(5.0573e+03 / h_barrel_MC_deno_70to100[_WZ]                 ->Integral());
    h_barrel_MC_deno_70to100[_ZZ]                 ->Scale(7.2227e+02 / h_barrel_MC_deno_70to100[_ZZ]                 ->Integral());
    h_barrel_MC_deno_70to100[_tW]                 ->Scale(1.2115e+04 / h_barrel_MC_deno_70to100[_tW]                 ->Integral());
    h_barrel_MC_deno_70to100[_tbarW]              ->Scale(1.2028e+04 / h_barrel_MC_deno_70to100[_tbarW]              ->Integral());
    h_barrel_MC_deno_70to100[_ttbar]              ->Scale(2.3233e+05 / h_barrel_MC_deno_70to100[_ttbar]              ->Integral());
    h_barrel_MC_deno_70to100[_WJets]              ->Scale(3.4985e+06 / h_barrel_MC_deno_70to100[_WJets]              ->Integral());
    h_barrel_MC_deno_70to100[_QCDMuEnriched_Full] ->Scale(3.7908e+06 / h_barrel_MC_deno_70to100[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_deno_100to500[_DY_Full]           ->Scale(1.0562e+05 / h_barrel_MC_deno_100to500[_DY_Full]            ->Integral());
    h_barrel_MC_deno_100to500[_WW]                ->Scale(1.4394e+04 / h_barrel_MC_deno_100to500[_WW]                 ->Integral());
    h_barrel_MC_deno_100to500[_WZ]                ->Scale(5.1568e+03 / h_barrel_MC_deno_100to500[_WZ]                 ->Integral());
    h_barrel_MC_deno_100to500[_ZZ]                ->Scale(8.3444e+02 / h_barrel_MC_deno_100to500[_ZZ]                 ->Integral());
    h_barrel_MC_deno_100to500[_tW]                ->Scale(1.2316e+04 / h_barrel_MC_deno_100to500[_tW]                 ->Integral());
    h_barrel_MC_deno_100to500[_tbarW]             ->Scale(1.2003e+04 / h_barrel_MC_deno_100to500[_tbarW]              ->Integral());
    h_barrel_MC_deno_100to500[_ttbar]             ->Scale(2.7302e+05 / h_barrel_MC_deno_100to500[_ttbar]              ->Integral());
    h_barrel_MC_deno_100to500[_WJets]             ->Scale(1.3415e+06 / h_barrel_MC_deno_100to500[_WJets]              ->Integral());
    h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]->Scale(7.9078e+05 / h_barrel_MC_deno_100to500[_QCDMuEnriched_Full] ->Integral());

    h_endcap_MC_deno_50to70[_DY_Full]             ->Scale(8.4551e+05 / h_endcap_MC_deno_50to70[_DY_Full]            ->Integral());
    h_endcap_MC_deno_50to70[_WW]                  ->Scale(2.4238e+04 / h_endcap_MC_deno_50to70[_WW]                 ->Integral());
    h_endcap_MC_deno_50to70[_WZ]                  ->Scale(8.1288e+03 / h_endcap_MC_deno_50to70[_WZ]                 ->Integral());
    h_endcap_MC_deno_50to70[_ZZ]                  ->Scale(1.5537e+03 / h_endcap_MC_deno_50to70[_ZZ]                 ->Integral());
    h_endcap_MC_deno_50to70[_tW]                  ->Scale(7.6181e+03 / h_endcap_MC_deno_50to70[_tW]                 ->Integral());
    h_endcap_MC_deno_50to70[_tbarW]               ->Scale(6.0976e+03 / h_endcap_MC_deno_50to70[_tbarW]              ->Integral());
    h_endcap_MC_deno_50to70[_ttbar]               ->Scale(1.4028e+05 / h_endcap_MC_deno_50to70[_ttbar]              ->Integral());
    h_endcap_MC_deno_50to70[_WJets]               ->Scale(7.4277e+06 / h_endcap_MC_deno_50to70[_WJets]              ->Integral());
    h_endcap_MC_deno_50to70[_QCDMuEnriched_Full]  ->Scale(6.7609e+06 / h_endcap_MC_deno_50to70[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_deno_70to100[_DY_Full]            ->Scale(2.8134e+05 / h_endcap_MC_deno_70to100[_DY_Full]            ->Integral());
    h_endcap_MC_deno_70to100[_WW]                 ->Scale(1.9348e+04 / h_endcap_MC_deno_70to100[_WW]                 ->Integral());
    h_endcap_MC_deno_70to100[_WZ]                 ->Scale(3.8243e+03 / h_endcap_MC_deno_70to100[_WZ]                 ->Integral());
    h_endcap_MC_deno_70to100[_ZZ]                 ->Scale(8.4462e+02 / h_endcap_MC_deno_70to100[_ZZ]                 ->Integral());
    h_endcap_MC_deno_70to100[_tW]                 ->Scale(5.0896e+03 / h_endcap_MC_deno_70to100[_tW]                 ->Integral());
    h_endcap_MC_deno_70to100[_tbarW]              ->Scale(5.6085e+03 / h_endcap_MC_deno_70to100[_tbarW]              ->Integral());
    h_endcap_MC_deno_70to100[_ttbar]              ->Scale(1.6028e+05 / h_endcap_MC_deno_70to100[_ttbar]              ->Integral());
    h_endcap_MC_deno_70to100[_WJets]              ->Scale(2.6868e+06 / h_endcap_MC_deno_70to100[_WJets]              ->Integral());
    h_endcap_MC_deno_70to100[_QCDMuEnriched_Full] ->Scale(1.7141e+06 / h_endcap_MC_deno_70to100[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_deno_100to500[_DY_Full]           ->Scale(1.5491e+05 / h_endcap_MC_deno_100to500[_DY_Full]            ->Integral());
    h_endcap_MC_deno_100to500[_WW]                ->Scale(1.2298e+04 / h_endcap_MC_deno_100to500[_WW]                 ->Integral());
    h_endcap_MC_deno_100to500[_WZ]                ->Scale(2.4067e+03 / h_endcap_MC_deno_100to500[_WZ]                 ->Integral());
    h_endcap_MC_deno_100to500[_ZZ]                ->Scale(5.7851e+02 / h_endcap_MC_deno_100to500[_ZZ]                 ->Integral());
    h_endcap_MC_deno_100to500[_tW]                ->Scale(5.8994e+03 / h_endcap_MC_deno_100to500[_tW]                 ->Integral());
    h_endcap_MC_deno_100to500[_tbarW]             ->Scale(5.8587e+03 / h_endcap_MC_deno_100to500[_tbarW]              ->Integral());
    h_endcap_MC_deno_100to500[_ttbar]             ->Scale(1.0494e+05 / h_endcap_MC_deno_100to500[_ttbar]              ->Integral());
    h_endcap_MC_deno_100to500[_WJets]             ->Scale(9.1044e+05 / h_endcap_MC_deno_100to500[_WJets]              ->Integral());
    h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]->Scale(2.9986e+05 / h_endcap_MC_deno_100to500[_QCDMuEnriched_Full] ->Integral());

    h_barrel_MC_nume_50to70[_DY_Full]             ->Scale(6.6645e+05 / h_barrel_MC_nume_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_nume_50to70[_WW]                  ->Scale(4.7497e+04 / h_barrel_MC_nume_50to70[_WW]                 ->Integral());
    h_barrel_MC_nume_50to70[_WZ]                  ->Scale(1.1871e+04 / h_barrel_MC_nume_50to70[_WZ]                 ->Integral());
    h_barrel_MC_nume_50to70[_ZZ]                  ->Scale(1.5384e+03 / h_barrel_MC_nume_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_nume_50to70[_tW]                  ->Scale(1.9176e+04 / h_barrel_MC_nume_50to70[_tW]                 ->Integral());
    h_barrel_MC_nume_50to70[_tbarW]               ->Scale(1.9135e+04 / h_barrel_MC_nume_50to70[_tbarW]              ->Integral());
    h_barrel_MC_nume_50to70[_ttbar]               ->Scale(4.0713e+05 / h_barrel_MC_nume_50to70[_ttbar]              ->Integral());
    h_barrel_MC_nume_50to70[_WJets]               ->Scale(8.6321e+06 / h_barrel_MC_nume_50to70[_WJets]              ->Integral());
    h_barrel_MC_nume_50to70[_QCDMuEnriched_Full]  ->Scale(8.2803e+05 / h_barrel_MC_nume_50to70[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_nume_70to100[_DY_Full]            ->Scale(3.0764e+05 / h_barrel_MC_nume_70to100[_DY_Full]            ->Integral());
    h_barrel_MC_nume_70to100[_WW]                 ->Scale(3.1821e+04 / h_barrel_MC_nume_70to100[_WW]                 ->Integral());
    h_barrel_MC_nume_70to100[_WZ]                 ->Scale(8.1105e+03 / h_barrel_MC_nume_70to100[_WZ]                 ->Integral());
    h_barrel_MC_nume_70to100[_ZZ]                 ->Scale(1.0908e+03 / h_barrel_MC_nume_70to100[_ZZ]                 ->Integral());
    h_barrel_MC_nume_70to100[_tW]                 ->Scale(1.8311e+04 / h_barrel_MC_nume_70to100[_tW]                 ->Integral());
    h_barrel_MC_nume_70to100[_tbarW]              ->Scale(1.8198e+04 / h_barrel_MC_nume_70to100[_tbarW]              ->Integral());
    h_barrel_MC_nume_70to100[_ttbar]              ->Scale(3.7420e+05 / h_barrel_MC_nume_70to100[_ttbar]              ->Integral());
    h_barrel_MC_nume_70to100[_WJets]              ->Scale(3.1926e+06 / h_barrel_MC_nume_70to100[_WJets]              ->Integral());
    h_barrel_MC_nume_70to100[_QCDMuEnriched_Full] ->Scale(2.1451e+05 / h_barrel_MC_nume_70to100[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_nume_100to500[_DY_Full]           ->Scale(1.0814e+05 / h_barrel_MC_nume_100to500[_DY_Full]            ->Integral());
    h_barrel_MC_nume_100to500[_WW]                ->Scale(1.3634e+04 / h_barrel_MC_nume_100to500[_WW]                 ->Integral());
    h_barrel_MC_nume_100to500[_WZ]                ->Scale(3.6344e+03 / h_barrel_MC_nume_100to500[_WZ]                 ->Integral());
    h_barrel_MC_nume_100to500[_ZZ]                ->Scale(6.1481e+02 / h_barrel_MC_nume_100to500[_ZZ]                 ->Integral());
    h_barrel_MC_nume_100to500[_tW]                ->Scale(1.7462e+04 / h_barrel_MC_nume_100to500[_tW]                 ->Integral());
    h_barrel_MC_nume_100to500[_tbarW]             ->Scale(1.7216e+04 / h_barrel_MC_nume_100to500[_tbarW]              ->Integral());
    h_barrel_MC_nume_100to500[_ttbar]             ->Scale(2.7310e+05 / h_barrel_MC_nume_100to500[_ttbar]              ->Integral());
    h_barrel_MC_nume_100to500[_WJets]             ->Scale(1.2825e+06 / h_barrel_MC_nume_100to500[_WJets]              ->Integral());
    h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->Scale(4.5879e+04 / h_barrel_MC_nume_100to500[_QCDMuEnriched_Full] ->Integral());

    h_endcap_MC_nume_50to70[_DY_Full]             ->Scale(6.5670e+05 / h_endcap_MC_nume_50to70[_DY_Full]            ->Integral());
    h_endcap_MC_nume_50to70[_WW]                  ->Scale(3.6323e+04 / h_endcap_MC_nume_50to70[_WW]                 ->Integral());
    h_endcap_MC_nume_50to70[_WZ]                  ->Scale(9.3875e+03 / h_endcap_MC_nume_50to70[_WZ]                 ->Integral());
    h_endcap_MC_nume_50to70[_ZZ]                  ->Scale(1.6381e+03 / h_endcap_MC_nume_50to70[_ZZ]                 ->Integral());
    h_endcap_MC_nume_50to70[_tW]                  ->Scale(8.7482e+03 / h_endcap_MC_nume_50to70[_tW]                 ->Integral());
    h_endcap_MC_nume_50to70[_tbarW]               ->Scale(8.7211e+03 / h_endcap_MC_nume_50to70[_tbarW]              ->Integral());
    h_endcap_MC_nume_50to70[_ttbar]               ->Scale(1.9208e+05 / h_endcap_MC_nume_50to70[_ttbar]              ->Integral());
    h_endcap_MC_nume_50to70[_WJets]               ->Scale(7.2978e+06 / h_endcap_MC_nume_50to70[_WJets]              ->Integral());
    h_endcap_MC_nume_50to70[_QCDMuEnriched_Full]  ->Scale(9.5464e+05 / h_endcap_MC_nume_50to70[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_nume_70to100[_DY_Full]            ->Scale(3.3187e+05 / h_endcap_MC_nume_70to100[_DY_Full]            ->Integral());
    h_endcap_MC_nume_70to100[_WW]                 ->Scale(2.3422e+04 / h_endcap_MC_nume_70to100[_WW]                 ->Integral());
    h_endcap_MC_nume_70to100[_WZ]                 ->Scale(6.0969e+03 / h_endcap_MC_nume_70to100[_WZ]                 ->Integral());
    h_endcap_MC_nume_70to100[_ZZ]                 ->Scale(1.2274e+03 / h_endcap_MC_nume_70to100[_ZZ]                 ->Integral());
    h_endcap_MC_nume_70to100[_tW]                 ->Scale(7.8944e+03 / h_endcap_MC_nume_70to100[_tW]                 ->Integral());
    h_endcap_MC_nume_70to100[_tbarW]              ->Scale(7.7716e+03 / h_endcap_MC_nume_70to100[_tbarW]              ->Integral());
    h_endcap_MC_nume_70to100[_ttbar]              ->Scale(1.6741e+05 / h_endcap_MC_nume_70to100[_ttbar]              ->Integral());
    h_endcap_MC_nume_70to100[_WJets]              ->Scale(2.5483e+06 / h_endcap_MC_nume_70to100[_WJets]              ->Integral());
    h_endcap_MC_nume_70to100[_QCDMuEnriched_Full] ->Scale(2.3836e+05 / h_endcap_MC_nume_70to100[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_nume_100to500[_DY_Full]           ->Scale(1.1091e+05 / h_endcap_MC_nume_100to500[_DY_Full]            ->Integral());
    h_endcap_MC_nume_100to500[_WW]                ->Scale(8.7089e+03 / h_endcap_MC_nume_100to500[_WW]                 ->Integral());
    h_endcap_MC_nume_100to500[_WZ]                ->Scale(2.3481e+03 / h_endcap_MC_nume_100to500[_WZ]                 ->Integral());
    h_endcap_MC_nume_100to500[_ZZ]                ->Scale(5.1436e+02 / h_endcap_MC_nume_100to500[_ZZ]                 ->Integral());
    h_endcap_MC_nume_100to500[_tW]                ->Scale(6.1123e+03 / h_endcap_MC_nume_100to500[_tW]                 ->Integral());
    h_endcap_MC_nume_100to500[_tbarW]             ->Scale(6.1855e+03 / h_endcap_MC_nume_100to500[_tbarW]              ->Integral());
    h_endcap_MC_nume_100to500[_ttbar]             ->Scale(1.0701e+05 / h_endcap_MC_nume_100to500[_ttbar]              ->Integral());
    h_endcap_MC_nume_100to500[_WJets]             ->Scale(9.3238e+05 / h_endcap_MC_nume_100to500[_WJets]              ->Integral());
    h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->Scale(4.8543e+04 / h_endcap_MC_nume_100to500[_QCDMuEnriched_Full] ->Integral());

    // barrel denominator
    THStack * s_barrel_deno = new THStack("s_barrel_deno", "");
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WW]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WW]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WW]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_tW]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_tW]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_tW]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_deno_50to70[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_70to100[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_deno_100to500[_QCDMuEnriched_Full]);

    h_barrel_data_deno_50to70->Add(h_barrel_data_deno_70to100);
    h_barrel_data_deno_50to70->Add(h_barrel_data_deno_100to500);

    // endcap denominator
    THStack * s_endcap_deno = new THStack("s_endcap_deno", "");
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_WW]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_WW]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_WW]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_tW]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_tW]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_tW]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_deno_50to70[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_deno_70to100[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_deno_100to500[_QCDMuEnriched_Full]);

    h_endcap_data_deno_50to70->Add(h_endcap_data_deno_70to100);
    h_endcap_data_deno_50to70->Add(h_endcap_data_deno_100to500);

    // barrel numerator
    THStack * s_barrel_nume = new THStack("s_barrel_nume", "");
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_WW]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_WW]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_WW]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_WZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_WZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_WZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_ZZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_ZZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_ZZ]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_DY_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_DY_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_DY_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_QCDMuEnriched_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_QCDMuEnriched_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]);

    h_barrel_data_nume_50to70->Add(h_barrel_data_nume_70to100);
    h_barrel_data_nume_50to70->Add(h_barrel_data_nume_100to500);

    // endcap numerator
    THStack * s_endcap_nume = new THStack("s_endcap_nume", "");
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_WW]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_WW]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_WW]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_WZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_WZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_WZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_ZZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_ZZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_ZZ]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_QCDMuEnriched_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_QCDMuEnriched_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]);

    h_endcap_data_nume_50to70->Add(h_endcap_data_nume_70to100);
    h_endcap_data_nume_50to70->Add(h_endcap_data_nume_100to500);

// ------------------------------------- DRAWING -----------------------------------------------------
    myRatioPlot_t *RP_barrel_deno = new myRatioPlot_t("RP_barrel_deno", s_barrel_deno, h_barrel_data_deno_50to70);
    myRatioPlot_t *RP_endcap_deno = new myRatioPlot_t("RP_endcap_deno", s_endcap_deno, h_endcap_data_deno_50to70);
    myRatioPlot_t *RP_barrel_nume = new myRatioPlot_t("RP_barrel_nume", s_barrel_nume, h_barrel_data_nume_50to70);
    myRatioPlot_t *RP_endcap_nume = new myRatioPlot_t("RP_endcap_nume", s_endcap_nume, h_endcap_data_nume_50to70);
    RP_barrel_deno->SetPlots("p_{T} (#mu_{barrel}^{deno})", 0, 5);
    RP_endcap_deno->SetPlots("p_{T} (#mu_{endcap}^{deno})", 0, 5);
    RP_barrel_nume->SetPlots("p_{T} (#mu_{barrel}^{nume})", 0, 5);
    RP_endcap_nume->SetPlots("p_{T} (#mu_{endcap}^{nume})", 0, 5);

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
    RP_endcap_deno->ImportLegend(legend);
    RP_barrel_nume->ImportLegend(legend);
    RP_endcap_nume->ImportLegend(legend);
    RP_barrel_deno->Draw(1, 1e7, 0);
    RP_endcap_deno->Draw(1, 1e7, 0);
    RP_barrel_nume->Draw(1, 1e7, 0);
    RP_endcap_nume->Draw(1, 1e7, 0);
} // End of Fit_HistDrawer()



///----------------------- ESTIMATIONS ------------------------- ///
void E_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";

    f = new TFile(Dir+"QCDest_E.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special];
    TH1D *h_QCD_est;
    THStack * s_mass_wQCD = new THStack("s_mass_wQCD", "");
    THStack * s_mass_woQCD = new THStack("s_mass_woQCD", "");
    Color_t color = kBlack;

    // Loop over all processes (adding all histograms)
    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        h_mass[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
        }

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_QCDEMEnriched_Normal) color = kRed + 3;

        if (pr < _DoubleEG_B) // MC coloring
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

        // Adding hists to THStacks
        if (pr < _DoubleEG_B)
        {
            s_mass_wQCD->Add(h_mass[pr]);
        }
        if (pr < _QCDEMEnriched_20to30)
        {
            s_mass_woQCD->Add(h_mass[pr]);
        }

        // Adding up for convenience
        if (pr == _DY_10to50)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
        }
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
        }
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
        }
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
        }
        else if (pr == _QCDEMEnriched_20to30)
        {
            h_mass[_QCDEMEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDEMEnriched_Full")));
            h_mass[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_QCDEMEnriched_Normal)
        {
            h_mass[_QCDEMEnriched_Full]->Add(h_mass[pr]);
        }
        else if (pr == _DoubleEG_B)
        {
            h_mass[_DoubleEG_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DoubleEG_Full")));
            h_mass[_DoubleEG_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_SingleMuon_Normal)
        {
            h_mass[_DoubleEG_Full]->Add(h_mass[pr]);
        }


        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

    } // End of pr iteration

    // QCD estimation
    h_QCD_est = ((TH1D*)(h_mass[_DoubleEG_Full]->Clone("h_QCD_est")));
    Double_t err_data=0, int_data=0, err_qcd=0, int_qcd=0;
    int_data = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_data);
    h_QCD_est->SetTitle("");
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Add(h_mass[_DY_Full], -1);
    h_QCD_est->Add(h_mass[_ttbar_Full], -1);
    h_QCD_est->Add(h_mass[_VVnST], -1);
    h_QCD_est->Add(h_mass[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kBlack);
    h_QCD_est->SetMarkerStyle(0);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "QCD est events: " << int_qcd << "+-" << err_qcd << endl;

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass_wQCD = new myRatioPlot_t("c_mass_wQCD", s_mass_wQCD, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_DoubleEG_Full]);

    RP_mass_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend * legend_wQCD = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend_wQCD->AddEntry(h_mass[_DoubleEG_B], "Data", "pl");
    legend_wQCD->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_wQCD->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_wQCD->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend_wQCD->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend_wQCD->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend_wQCD->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend_wQCD->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend_wQCD->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    TLegend * legend_woQCD = ((TLegend*)(legend_wQCD->Clone()));
    legend_wQCD->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wQCD->ImportLegend(legend_wQCD);
    RP_mass_woQCD->ImportLegend(legend_woQCD);

    RP_mass_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD->Draw(1e-3, 1e3, 1);

    // Drawing estimated QCD
    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    TCanvas * c_QCD_est = new TCanvas("c_QCD_est", "QCD est", 750, 850);
    c_QCD_est->SetTopMargin(0.05);
    c_QCD_est->SetRightMargin(0.05);
    c_QCD_est->SetBottomMargin(0.15);
    c_QCD_est->SetLeftMargin(0.15);
    h_QCD_est->Draw("BAR");
    h_QCD_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_QCD_est->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est->GetXaxis()->SetMoreLogLabels();
    h_QCD_est->GetXaxis()->SetNoExponent();
    h_QCD_est->GetYaxis()->SetTitle("Number of events");
    h_QCD_est->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est->GetYaxis()->SetMoreLogLabels();
    h_QCD_est->GetYaxis()->SetNoExponent();
//    h_QCD_est->GetYaxis()->SetRangeUser(0, 300);

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_QCD_est->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_TEMPLATE_DOWN.root", "READ");
        TH1D *h_ratio, *h_up, *h_down, *h_fullsysterr, *h_fullerr;
        f_ratio->GetObject("h_QCD_est", h_ratio);
        f_up->GetObject("h_QCD_est", h_up);
        f_down->GetObject("h_QCD_est", h_down);
        h_fullsysterr = ((TH1D*)(h_ratio->Clone("h_QCD_fullsysterr")));
        h_fullerr = ((TH1D*)(h_ratio->Clone("h_QCD_fullerr")));

        for (Int_t i=1; i<h_ratio->GetSize()-1; i++)
        {
            // Systematic errors
            h_ratio->SetBinContent(i, fabs(h_ratio->GetBinContent(i)-h_QCD_est->GetBinContent(i)));
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_fullsysterr->SetBinContent(i, sqrt(h_ratio->GetBinContent(i)*h_ratio->GetBinContent(i)+h_up->GetBinContent(i)*h_up->GetBinContent(i)));
            // Statistical errors
            if (h_QCD_est->GetBinContent(i) > 0) h_QCD_est->SetBinError(i, sqrt(h_QCD_est->GetBinContent(i)));
            else h_QCD_est->SetBinError(i, 1);
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_QCD_est->GetBinError(i)*h_QCD_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
        }
        cout << "Systematic error: " << h_fullsysterr->Integral() << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_QCD_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_QCD_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);
        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);

        f_ratio->Close();
        f_up->Close();
        f_down->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullerr->Write();

        h_QCD_est->GetYaxis()->SetRangeUser(0, 200);
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
    }// End of if (systErr)

    l_QCD_est->Draw();
    c_QCD_est->SetLogx();
    c_QCD_est->SetGridx();
    c_QCD_est->SetGridy();
    c_QCD_est->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_E.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_E.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root COULD NOT BE CLOSED!\n" << endl;

} // End of E_QCDest_HistDrawer()


void E_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TString Dir = "/media/sf_DATA/FR/Electron/";
    TFile *f = new TFile(Dir+"WJETest_E.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special];
    TH1D *h_WJET_est;
    THStack * s_mass_wWJET = new THStack("s_mass_wWJET", "");
    THStack * s_mass_woWJET = new THStack("s_mass_woWJET", "");
    Color_t color = kBlack;

    // Getting data-driven QCD to subtract from data
    TFile *f_QCD = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root");
    TH1D *h_QCD_est;
    f_QCD->GetObject("h_QCD_est", h_QCD_est);
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Scale(2);
    s_mass_wWJET->Add(h_QCD_est);
    s_mass_woWJET->Add(h_QCD_est);

    for (Process_t pr=_DoubleEG_H; pr>=_DY_10to50; pr=previous(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        h_mass[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
        }

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_QCDEMEnriched_Normal) color = kRed + 3;

        if (pr < _DoubleEG_B) // MC coloring
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

        // Filling stack histograms
        if (pr < _EndOf_WJets_Normal)
        {
            s_mass_wWJET->Add(h_mass[pr]);
            if (pr < _WJets)
            {
                s_mass_woWJET->Add(h_mass[pr]);
            }
        }

        // Adding up for convenience
        if (pr == _DoubleEG_H)
        {
            h_mass[_DoubleEG_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_DoubleEG_Full]->SetDirectory(0);
        }
        else if (pr >= _DoubleEG_B)
        {
            h_mass[_DoubleEG_Full]->Add(h_mass[pr]);
        }
        else if (pr == _QCDEMEnriched_300toInf)
        {
            h_mass[_QCDEMEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDEMEnriched_Full")));
            h_mass[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else if (pr >= _QCDEMEnriched_20to30)
        {
            h_mass[_QCDEMEnriched_Full]->Add(h_mass[pr]);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
        }
        else if (pr >= _WJets)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
        }
        else if (pr == _WW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
        }
        else if (pr >= _tW)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
        }
        else if (pr == _ttbar_1000toInf)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr >= _ttbar)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
        }
        else if (pr == _DY_2000to3000)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
        }
        else if (pr >= _DY_10to50)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
        }

        if (pr == _QCDEMEnriched_20to30) pr = _EndOf_WJets_Normal; // next -- WJets_ext2v5
        if (pr == _ttbar) pr = _EndOf_DY_Normal; // next -- DY_2000to3000

    } // End of pr iteration

    // W+Jets estimation by subtraction
    Double_t int_data=0, err_data=0, int_wjet=0, err_wjet=0;
    h_WJET_est = ((TH1D*)(h_mass[_DoubleEG_Full]->Clone("h_WJET_est")));
    int_data = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_data);
    h_WJET_est->SetTitle("");
    h_WJET_est->SetDirectory(0);
    h_WJET_est->Add(h_mass[_DY_Full], -1);
    h_WJET_est->Add(h_mass[_ttbar_Full], -1);
    h_WJET_est->Add(h_mass[_VVnST], -1);
    h_WJET_est->Add(h_QCD_est, -1);
    removeNegativeBins(h_WJET_est);
    h_WJET_est->SetFillColor(kRed - 2);
    h_WJET_est->SetLineColor(kRed - 2);
    int_wjet = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_wjet);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "WJets est events (subraction): " << int_wjet << "+-" << err_wjet << endl;

    myRatioPlot_t *RP_mass_wWJET = new myRatioPlot_t("c_mass_wWJET", s_mass_wWJET, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woWJET = new myRatioPlot_t("c_mass_woWJET", s_mass_woWJET, h_mass[_DoubleEG_Full]);

    RP_mass_wWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend * legend_wWJET = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend_wWJET->AddEntry(h_mass[_DoubleEG_B], "Data", "pl");
    legend_wWJET->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_wWJET->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend_wWJET->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend_wWJET->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend_wWJET->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    TLegend * legend_woWJET = ((TLegend*)(legend_wWJET->Clone()));
    legend_wWJET->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend_wWJET->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend_woWJET->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);

    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");

    // Draw WJets from simple subtraction (opposite sign)
    TCanvas * c_WJET_est = new TCanvas("c_WJET_est", "W+Jets est", 750, 850);
    c_WJET_est->SetTopMargin(0.05);
    c_WJET_est->SetRightMargin(0.05);
    c_WJET_est->SetBottomMargin(0.15);
    c_WJET_est->SetLeftMargin(0.17);
    h_WJET_est->Draw("hist");
    h_WJET_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_WJET_est->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est->GetXaxis()->SetMoreLogLabels();
    h_WJET_est->GetXaxis()->SetNoExponent();
    h_WJET_est->GetYaxis()->SetTitle("Number of events");
    h_WJET_est->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est->GetYaxis()->SetMoreLogLabels();
    h_WJET_est->GetYaxis()->SetNoExponent();
//    h_WJET_est->GetYaxis()->SetRangeUser(0.01, 1e3);
    l_WJET_est->Draw();
    c_WJET_est->SetLogx();
    c_WJET_est->SetGridx();
    c_WJET_est->SetGridy();
    c_WJET_est->Update();


    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_WJET_est->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_TEMPLATE_DOWN.root", "READ");
        TH1D *h_ratio, *h_up, *h_down, *h_fullsysterr, *h_fullerr;
        f_ratio->GetObject("h_WJET_est", h_ratio);
        f_up->GetObject("h_WJET_est", h_up);
        f_down->GetObject("h_WJET_est", h_down);
        h_fullsysterr = ((TH1D*)(h_ratio->Clone("h_WJET_fullsysterr")));
        h_fullerr = ((TH1D*)(h_ratio->Clone("h_WJET_fullerr")));

        for (Int_t i=1; i<h_ratio->GetSize()-1; i++)
        {
            // Systematic errors
            h_ratio->SetBinContent(i, fabs(h_ratio->GetBinContent(i)-h_WJET_est->GetBinContent(i)));
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_fullsysterr->SetBinContent(i, sqrt(h_ratio->GetBinContent(i)*h_ratio->GetBinContent(i)+h_up->GetBinContent(i)*h_up->GetBinContent(i)));
            // Statistical errors
            if (h_WJET_est->GetBinContent(i) > 0) h_WJET_est->SetBinError(i, sqrt(h_WJET_est->GetBinContent(i)));
            else h_WJET_est->SetBinError(i, 1);
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_WJET_est->GetBinError(i)*h_WJET_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
        }
        cout << "Systematic error: " << h_fullsysterr->Integral() << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_WJET_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_WJET_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);
        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);
        f_ratio->Close();
        f_up->Close();
        f_down->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullerr->Write();

        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
    }// End of if (systErr)

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_E.root" << " has been closed successfully." << endl;
    else cout << "FILE " << Dir+"WJETest_E.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedEE/Histos/EstWJets_EE.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedEE/Histos/EstWJets_EE.root COULD NOT BE CLOSED!\n" << endl;

} // End of E_WJETest_HistDrawer()




void Mu_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";

    f = new TFile(Dir+"QCDest_Mu.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special];
    TH1D *h_QCD_est, *h_QCD_est_SS;
    THStack * s_mass_wQCD = new THStack("s_mass_wQCD", "");
    THStack * s_mass_woQCD = new THStack("s_mass_woQCD", "");
    THStack * s_mass_SS_wQCD = new THStack("s_mass_SS_wQCD", "");
    THStack * s_mass_SS_woQCD = new THStack("s_mass_SS_woQCD", "");
    TH1D * h_data_pT_lead;
    TH1D * h_data_pT_sublead;
    Color_t color = kBlack;

    // Loop over all processes (adding all histograms)
    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
        }

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
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
        }

        // Adding hists to THStacks
        if (pr < _SingleMuon_B)
        {
            s_mass_wQCD->Add(h_mass[pr]);
            s_mass_SS_wQCD->Add(h_mass_SS[pr]);
        }
        if (pr < _QCDMuEnriched_15to20)
        {
            s_mass_woQCD->Add(h_mass[pr]);
            s_mass_SS_woQCD->Add(h_mass_SS[pr]);
        }

        // Adding up for convenience
        if (pr == _DY_10to50)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
        }
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _QCDMuEnriched_15to20)
        {
            h_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDMuEnriched_Full")));
            h_mass[_QCDMuEnriched_Full]->SetDirectory(0);
            h_mass_SS[_QCDMuEnriched_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_QCDMuEnriched_Full")));
            h_mass_SS[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_QCDMuEnriched_Normal)
        {
            h_mass[_QCDMuEnriched_Full]->Add(h_mass[pr]);
            h_mass_SS[_QCDMuEnriched_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
            h_mass_SS[_SingleMuon_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_SingleMuon_Full")));
            h_mass_SS[_SingleMuon_Full]->SetDirectory(0);
            f->GetObject("h_pT_lead_"+Mgr.Procname[pr], h_data_pT_lead);
            f->GetObject("h_pT_sublead_"+Mgr.Procname[pr], h_data_pT_sublead);
            h_data_pT_lead->SetDirectory(0);
            h_data_pT_sublead->SetDirectory(0);
        }
        else if (pr < _EndOf_SingleMuon_Normal)
        {
            h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
            h_mass_SS[_SingleMuon_Full]->Add(h_mass_SS[pr]);
            TH1D * temp[2];
            f->GetObject("h_pT_lead_"+Mgr.Procname[pr], temp[0]);
            f->GetObject("h_pT_sublead_"+Mgr.Procname[pr], temp[1]);
            h_data_pT_lead->Add(temp[0]);
            h_data_pT_sublead->Add(temp[1]);
        }


        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

    } // End of pr iteration

    // QCD estimation
    h_QCD_est = ((TH1D*)(h_mass[_SingleMuon_Full]->Clone("h_QCD_est")));
    Double_t err_data=0, int_data=0, err_qcd=0, int_qcd=0;
    int_data = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_data);
    h_QCD_est->SetTitle("");
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Add(h_mass[_DY_Full], -1);
    h_QCD_est->Add(h_mass[_ttbar_Full], -1);
    h_QCD_est->Add(h_mass[_VVnST], -1);
    h_QCD_est->Add(h_mass[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kBlack);
    h_QCD_est->SetMarkerStyle(0);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "QCD est events: " << int_qcd << "+-" << err_qcd << endl;

    h_QCD_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_QCD_est_SS")));
    h_QCD_est_SS->SetTitle("");
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_QCD_est_SS->Add(h_mass_SS[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est_SS);
    h_QCD_est_SS->SetFillColor(kRed + 3);
    h_QCD_est_SS->SetLineColor(kRed + 3);

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass_wQCD = new myRatioPlot_t("c_mass_wQCD", s_mass_wQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_wQCD = new myRatioPlot_t("c_mass_SS_wQCD", s_mass_SS_wQCD, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_woQCD = new myRatioPlot_t("c_mass_SS_woQCD", s_mass_SS_woQCD, h_mass_SS[_SingleMuon_Full]);

    RP_mass_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);

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
    RP_mass_SS_wQCD->ImportLegend(legend_wQCD);
    RP_mass_SS_woQCD->ImportLegend(legend_woQCD);

    RP_mass_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD->Draw(1e-3, 1e3, 1);
    RP_mass_SS_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_SS_woQCD->Draw(1e-3, 1e3, 1);

    // Drawing estimated QCD
    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    TCanvas * c_QCD_est = new TCanvas("c_QCD_est", "QCD est", 750, 850);
    c_QCD_est->SetTopMargin(0.05);
    c_QCD_est->SetRightMargin(0.05);
    c_QCD_est->SetBottomMargin(0.15);
    c_QCD_est->SetLeftMargin(0.15);
    h_QCD_est->Draw("BAR");
    h_QCD_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_QCD_est->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est->GetXaxis()->SetMoreLogLabels();
    h_QCD_est->GetXaxis()->SetNoExponent();
    h_QCD_est->GetYaxis()->SetTitle("Number of events");
    h_QCD_est->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est->GetYaxis()->SetMoreLogLabels();
    h_QCD_est->GetYaxis()->SetNoExponent();
//    h_QCD_est->GetYaxis()->SetRangeUser(0, 300);

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_QCD_est->Write();
    h_QCD_est_SS->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu_TEMPLATE_DOWN.root", "READ");
        TH1D *h_ratio, *h_up, *h_down, *h_fullsysterr, *h_fullerr;
        f_ratio->GetObject("h_QCD_est", h_ratio);
        f_up->GetObject("h_QCD_est", h_up);
        f_down->GetObject("h_QCD_est", h_down);
        h_fullsysterr = ((TH1D*)(h_ratio->Clone("h_QCD_fullsysterr")));
        h_fullerr = ((TH1D*)(h_ratio->Clone("h_QCD_fullerr")));

        for (Int_t i=1; i<h_ratio->GetSize()-1; i++)
        {
            // Systematic errors
            h_ratio->SetBinContent(i, fabs(h_ratio->GetBinContent(i)-h_QCD_est->GetBinContent(i)));
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_fullsysterr->SetBinContent(i, sqrt(h_ratio->GetBinContent(i)*h_ratio->GetBinContent(i)+h_up->GetBinContent(i)*h_up->GetBinContent(i)));
            // Statistical errors
            if (h_QCD_est->GetBinContent(i) > 0) h_QCD_est->SetBinError(i, sqrt(h_QCD_est->GetBinContent(i)));
            else h_QCD_est->SetBinError(i, 1);
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_QCD_est->GetBinError(i)*h_QCD_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
        }
        cout << "Systematic error: " << h_fullsysterr->Integral() << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_QCD_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_QCD_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);
        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);

        f_ratio->Close();
        f_up->Close();
        f_down->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullerr->Write();

        h_QCD_est->GetYaxis()->SetRangeUser(0, 200);
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
    }// End of if (systErr)

    l_QCD_est->Draw();
    c_QCD_est->SetLogx();
    c_QCD_est->SetGridx();
    c_QCD_est->SetGridy();
    c_QCD_est->Update();

    TCanvas * c_QCD_est_SS = new TCanvas("c_QCD_SS_est", "QCD est (same sign)", 750, 850);
    c_QCD_est_SS->SetTopMargin(0.05);
    c_QCD_est_SS->SetRightMargin(0.05);
    c_QCD_est_SS->SetBottomMargin(0.15);
    c_QCD_est_SS->SetLeftMargin(0.15);
    h_QCD_est_SS->Draw("hist");
    h_QCD_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]");
    h_QCD_est_SS->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est_SS->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est_SS->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est_SS->GetXaxis()->SetMoreLogLabels();
    h_QCD_est_SS->GetXaxis()->SetNoExponent();
    h_QCD_est_SS->GetYaxis()->SetTitle("Number of events");
    h_QCD_est_SS->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est_SS->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est_SS->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est_SS->GetYaxis()->SetMoreLogLabels();
    h_QCD_est_SS->GetYaxis()->SetNoExponent();
//    h_QCD_est_SS->GetYaxis()->SetRangeUser(0, 300);
    l_QCD_est->Draw();
    c_QCD_est_SS->SetLogx();
    c_QCD_est_SS->SetGridx();
    c_QCD_est_SS->SetGridy();
    c_QCD_est_SS->Update();

/*
    TCanvas * c_pT_lead = new TCanvas("c_pT_lead", "pT lead", 750, 850);
//    TF1 * fit_lead = new TF1("fit_lead", "expo(0)", 1, 200);
    c_pT_lead->SetTopMargin(0.05);
    c_pT_lead->SetRightMargin(0.05);
    c_pT_lead->SetBottomMargin(0.15);
    c_pT_lead->SetLeftMargin(0.15);
//    h_data_pT_lead->Fit("fit_lead");
    h_data_pT_lead->Draw();
    h_data_pT_lead->GetXaxis()->SetTitle("p_{T} (#mu_{lead}) [GeV/c]");
    h_data_pT_lead->SetTitle("");
    h_data_pT_lead->GetXaxis()->SetTitleSize(0.062);
    h_data_pT_lead->GetXaxis()->SetTitleOffset(0.9);
    h_data_pT_lead->GetXaxis()->SetLabelSize(0.048);
    h_data_pT_lead->GetXaxis()->SetMoreLogLabels();
    h_data_pT_lead->GetXaxis()->SetNoExponent();
    h_data_pT_lead->GetYaxis()->SetTitle("Number of events");
    h_data_pT_lead->GetYaxis()->SetTitleSize(0.05);
    h_data_pT_lead->GetYaxis()->SetTitleOffset(1.2);
    h_data_pT_lead->GetYaxis()->SetLabelSize(0.043);
    c_pT_lead->SetLogy();
    c_pT_lead->SetGridx();
    c_pT_lead->SetGridy();
    c_pT_lead->Update();

    TCanvas * c_pT_sublead = new TCanvas("c_pT_sublead", "pT sublead", 750, 850);
//    TF1 * fit_sublead = new TF1("fit_sublead", "expo(0)", 1, 200);
    c_pT_sublead->SetTopMargin(0.05);
    c_pT_sublead->SetRightMargin(0.05);
    c_pT_sublead->SetBottomMargin(0.15);
    c_pT_sublead->SetLeftMargin(0.15);
//    h_data_pT_sublead->Fit("fit_sublead");
    h_data_pT_sublead->Draw("hist");
    h_data_pT_sublead->GetXaxis()->SetTitle("p_{T} (#mu_{sublead}) [GeV/c]");
    h_data_pT_sublead->SetTitle("");
    h_data_pT_sublead->GetXaxis()->SetTitleSize(0.062);
    h_data_pT_sublead->GetXaxis()->SetTitleOffset(0.9);
    h_data_pT_sublead->GetXaxis()->SetLabelSize(0.048);
    h_data_pT_sublead->GetXaxis()->SetMoreLogLabels();
    h_data_pT_sublead->GetXaxis()->SetNoExponent();
    h_data_pT_sublead->GetYaxis()->SetTitle("Number of events");
    h_data_pT_sublead->GetYaxis()->SetTitleSize(0.05);
    h_data_pT_sublead->GetYaxis()->SetTitleOffset(1.2);
    h_data_pT_sublead->GetYaxis()->SetLabelSize(0.043);
    c_pT_sublead->SetLogy();
    c_pT_sublead->SetGridx();
    c_pT_sublead->SetGridy();
    c_pT_sublead->Update();
*/
    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_Mu.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_Mu.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu.root COULD NOT BE CLOSED!\n" << endl;

} // End of Mu_QCDest_HistDrawer()


void Mu_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TString Dir = "/media/sf_DATA/FR/Muon/";
    TFile *f = new TFile(Dir+"WJETest_Mu.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special];
    TH1D *h_WJET_est, *h_WJET_est_SS;
    THStack * s_mass_wWJET = new THStack("s_mass_wWJET", "");
    THStack * s_mass_SS_wWJET = new THStack("s_mass_SS_wWJET", "");
    THStack * s_mass_woWJET = new THStack("s_mass_woWJET", "");
    THStack * s_mass_SS_woWJET = new THStack("s_mass_SS_woWJET", "");
    Color_t color = kBlack;

    // Getting data-driven QCD to subtract from data
    TFile *f_QCD = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu.root");
    TH1D *h_QCD_est, *h_QCD_est_SS;
    f_QCD->GetObject("h_QCD_est", h_QCD_est);
    f_QCD->GetObject("h_QCD_est_SS", h_QCD_est_SS);
    h_QCD_est->SetDirectory(0);
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est->Scale(2);
    cout << "Default Opposite sign QCD integral: " << h_QCD_est->Integral(1, 30) << endl;
    h_QCD_est->Scale(6.6329e+03 / h_QCD_est->Integral(1, 30));
    h_QCD_est_SS->Scale(2);
    s_mass_wWJET->Add(h_QCD_est);
    s_mass_woWJET->Add(h_QCD_est);
    s_mass_SS_wWJET->Add(h_QCD_est_SS);
    s_mass_SS_woWJET->Add(h_QCD_est_SS);

    for (Process_t pr=_SingleMuon_H; pr>=_DY_10to50; pr=previous(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
        }

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
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
        }

        // Filling stack histograms
        if (pr < _EndOf_WJets_Normal)
        {
            s_mass_wWJET->Add(h_mass[pr]);
            s_mass_SS_wWJET->Add(h_mass_SS[pr]);
            if (pr < _WJets)
            {
                s_mass_woWJET->Add(h_mass[pr]);
                s_mass_SS_woWJET->Add(h_mass_SS[pr]);
            }
        }

        // Adding up for convenience
        if (pr == _SingleMuon_H)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
            h_mass_SS[_SingleMuon_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_SingleMuon_Full")));
            h_mass_SS[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr >= _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
            h_mass_SS[_SingleMuon_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _QCDMuEnriched_1000toInf)
        {
            h_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDMuEnriched_Full")));
            h_mass[_QCDMuEnriched_Full]->SetDirectory(0);
            h_mass_SS[_QCDMuEnriched_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_QCDMuEnriched_Full")));
            h_mass_SS[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else if (pr >= _QCDMuEnriched_15to20)
        {
            h_mass[_QCDMuEnriched_Full]->Add(h_mass[pr]);
            h_mass_SS[_QCDMuEnriched_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
        }
        else if (pr >= _WJets)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _WW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
        }
        else if (pr >= _tW)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
        }
        else if (pr == _ttbar_1000toInf)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr >= _ttbar)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
        }
        else if (pr == _DY_2000to3000)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
        }
        else if (pr >= _DY_10to50)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
        }

        if (pr == _SingleMuon_B) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCD_Mu_1000toInf
        if (pr == _ttbar) pr = _EndOf_DY_Normal; // next -- DY_2000to3000

    } // End of pr iteration

    // W+Jets estimation by subtraction
    Double_t int_data=0, err_data=0, int_wjet=0, err_wjet=0;
    h_WJET_est = ((TH1D*)(h_mass[_SingleMuon_Full]->Clone("h_WJET_est")));
    int_data = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_data);
    h_WJET_est->SetTitle("");
    h_WJET_est->SetDirectory(0);
    h_WJET_est->Add(h_mass[_DY_Full], -1);
    h_WJET_est->Add(h_mass[_ttbar_Full], -1);
    h_WJET_est->Add(h_mass[_VVnST], -1);
    h_WJET_est->Add(h_QCD_est, -1);
    removeNegativeBins(h_WJET_est);
    h_WJET_est->SetFillColor(kRed - 2);
    h_WJET_est->SetLineColor(kRed - 2);  
    int_wjet = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_wjet);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "WJets est events (subraction): " << int_wjet << "+-" << err_wjet << endl;

    Double_t int_data_ss=0, err_data_ss=0, int_wjet_ss=0, err_wjet_ss=0;
    h_WJET_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_WJET_est_SS")));
    int_data_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_data_ss);
    h_WJET_est_SS->SetTitle("");
    h_WJET_est_SS->SetDirectory(0);
    h_WJET_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_WJET_est_SS->Add(h_QCD_est_SS, -1);
    removeNegativeBins(h_WJET_est_SS);
    h_WJET_est_SS->SetFillColor(kRed - 2);
    h_WJET_est_SS->SetLineColor(kRed - 2);
    int_wjet_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_wjet_ss);
    cout << "Data SS events: " << int_data_ss << "+-" << err_data_ss << endl;
    cout << "WJets SS est events (subraction): " << int_wjet_ss << "+-" << err_wjet_ss << endl;

    myRatioPlot_t *RP_mass_wWJET = new myRatioPlot_t("c_mass_wWJET", s_mass_wWJET, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woWJET = new myRatioPlot_t("c_mass_woWJET", s_mass_woWJET, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_wWJET = new myRatioPlot_t("c_mass_SS_wWJET", s_mass_SS_wWJET, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_woWJET = new myRatioPlot_t("c_mass_SS_woWJET", s_mass_SS_woWJET, h_mass_SS[_SingleMuon_Full]);

    RP_mass_wWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_wWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_woWJET->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);

    TLegend * legend_wWJET = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend_wWJET->AddEntry(h_mass[_SingleMuon_B], "Data", "pl");
    legend_wWJET->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_wWJET->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend_wWJET->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend_wWJET->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend_wWJET->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend_wWJET->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    TLegend * legend_woWJET = ((TLegend*)(legend_wWJET->Clone()));
    legend_wWJET->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend_wWJET->AddEntry(h_mass[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend_woWJET->AddEntry(h_mass[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);
    RP_mass_SS_wWJET->ImportLegend(legend_wWJET);
    RP_mass_SS_woWJET->ImportLegend(legend_woWJET);

    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);
    RP_mass_SS_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_SS_woWJET->Draw(1e-2, 1e5, 1);

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");

    // Draw WJets from simple subtraction (opposite sign)
    TCanvas * c_WJET_est = new TCanvas("c_WJET_est", "W+Jets est OS", 750, 850);
    c_WJET_est->SetTopMargin(0.05);
    c_WJET_est->SetRightMargin(0.05);
    c_WJET_est->SetBottomMargin(0.15);
    c_WJET_est->SetLeftMargin(0.17);
    h_WJET_est->Draw("hist");
    h_WJET_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_WJET_est->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est->GetXaxis()->SetMoreLogLabels();
    h_WJET_est->GetXaxis()->SetNoExponent();
    h_WJET_est->GetYaxis()->SetTitle("Number of events");
    h_WJET_est->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est->GetYaxis()->SetMoreLogLabels();
    h_WJET_est->GetYaxis()->SetNoExponent();
//    h_WJET_est->GetYaxis()->SetRangeUser(0.01, 1e3);
    l_WJET_est->Draw();
    c_WJET_est->SetLogx();
    c_WJET_est->SetGridx();
    c_WJET_est->SetGridy();
    c_WJET_est->Update();

    // Draw WJets from simple subtraction (same sign)
    TCanvas * c_WJET_est_SS = new TCanvas("c_WJET_est_SS", "W+Jets est SS", 750, 850);
    c_WJET_est_SS->SetTopMargin(0.05);
    c_WJET_est_SS->SetRightMargin(0.05);
    c_WJET_est_SS->SetBottomMargin(0.15);
    c_WJET_est_SS->SetLeftMargin(0.17);
    h_WJET_est_SS->Draw("hist");
    h_WJET_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]");
    h_WJET_est_SS->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est_SS->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est_SS->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est_SS->GetXaxis()->SetMoreLogLabels();
    h_WJET_est_SS->GetXaxis()->SetNoExponent();
    h_WJET_est_SS->GetYaxis()->SetTitle("Number of events");
    h_WJET_est_SS->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est_SS->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est_SS->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est_SS->GetYaxis()->SetMoreLogLabels();
    h_WJET_est_SS->GetYaxis()->SetNoExponent();
//    h_WJET_est_SS->GetYaxis()->SetRangeUser(0.01, 1e3);
    l_WJET_est->Draw();
    c_WJET_est_SS->SetLogx();
    c_WJET_est_SS->SetGridx();
    c_WJET_est_SS->SetGridy();
    c_WJET_est_SS->Update();

    // Draw WJets from template fitting (MC template)
    TH1D * h_WJET_fit_MC = ((TH1D*)(h_mass[_WJets_Full]->Clone("h_WJET_fit_MC")));
    h_WJET_fit_MC->SetTitle("");
//    h_WJET_fit_MC->Scale(6.6354e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 28, 17
    h_WJET_fit_MC->Scale(6.8312e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 28, 17; Mixed FR method with SF
//    h_WJET_fit_MC->Scale(7.0791e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 28, 17; Mixed FR method
//    h_WJET_fit_MC->Scale(1.2256e+02 / h_WJET_fit_MC->Integral()); // pT cuts: 52, 52
//    h_WJET_fit_MC->Scale(8.2256e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 52, 10
//    h_WJET_fit_MC->Scale(2.4489e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 52, 17
//    h_WJET_fit_MC->Scale(4.3978e+04 / h_WJET_fit_MC->Integral()); // pT cuts: 52, 0
    h_WJET_fit_MC->SetDirectory(0);
    h_WJET_fit_MC->SetFillColor(kRed - 2);
    h_WJET_fit_MC->SetLineColor(kRed - 2);
    TCanvas * c_WJET_fit_MC = new TCanvas("c_WJET_fit_MC", "W+Jets est fit MC", 750, 850);
    c_WJET_fit_MC->SetTopMargin(0.05);
    c_WJET_fit_MC->SetRightMargin(0.05);
    c_WJET_fit_MC->SetBottomMargin(0.15);
    c_WJET_fit_MC->SetLeftMargin(0.17);
    h_WJET_fit_MC->Draw("hist");
    h_WJET_fit_MC->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_WJET_fit_MC->GetXaxis()->SetTitleSize(0.062);
    h_WJET_fit_MC->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_fit_MC->GetXaxis()->SetLabelSize(0.048);
    h_WJET_fit_MC->GetXaxis()->SetMoreLogLabels();
    h_WJET_fit_MC->GetXaxis()->SetNoExponent();
    h_WJET_fit_MC->GetYaxis()->SetTitle("Number of events");
    h_WJET_fit_MC->GetYaxis()->SetTitleSize(0.05);
    h_WJET_fit_MC->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_fit_MC->GetYaxis()->SetLabelSize(0.043);
    h_WJET_fit_MC->GetYaxis()->SetMoreLogLabels();
    h_WJET_fit_MC->GetYaxis()->SetNoExponent();
//    h_WJET_fit_MC->GetYaxis()->SetRangeUser(0.01, 500);
    l_WJET_est->Draw();
    c_WJET_fit_MC->SetLogx();
    c_WJET_fit_MC->SetGridx();
    c_WJET_fit_MC->SetGridy();
    c_WJET_fit_MC->Update();

    // Draw WJets from template fitting (Same-sign data template BUT THE ESTIMATE IS FOR OPPOSITE SIGN!)
    TH1D * h_WJET_fit_SS = ((TH1D*)(h_WJET_est_SS->Clone("h_WJET_fit_SS")));
    h_WJET_fit_SS->SetTitle("");
//    h_WJET_fit_SS->Scale(5.3412e+03 / h_WJET_fit_SS->Integral()); // from full histogram fitting
    cout << "Same sign W+Jets integral: " << h_WJET_fit_SS->Integral(1,30) << endl;
//    h_WJET_fit_SS->Scale(2.2000e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (my FR 2)
//    h_WJET_fit_SS->Scale(2.7086e+04 / h_WJET_fit_SS->Integral(1,30)); // Test with (52, 2) GeV pT cuts
    h_WJET_fit_SS->Scale(2.5935e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (my FR 2)
//    h_WJET_fit_SS->Scale(3.3925e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (my FR 2)
//    h_WJET_fit_SS->Scale(3.5543e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS)
//    h_WJET_fit_SS->Scale(3.4181e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS (up))
//    h_WJET_fit_SS->Scale(3.7006e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS (down))
    Double_t int_wjet_est, err_wjet_est;
    int_wjet_est = h_WJET_fit_SS->IntegralAndError(1, h_WJET_fit_SS->GetSize()-2, err_wjet_est);
    cout << "WJets est events (FROM SS): " << int_wjet_est << "+-" << err_wjet_est << endl;
    h_WJET_fit_SS->SetDirectory(0);
    h_WJET_fit_SS->SetFillColor(kRed - 2);
    h_WJET_fit_SS->SetLineColor(kBlack);
    h_WJET_fit_SS->SetMarkerStyle(0);

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_WJET_est->Write();
    h_WJET_fit_MC->Write();
    h_WJET_fit_SS->Write();

    TCanvas * c_WJET_fit_SS = new TCanvas("c_WJET_fit_SS", "W+Jets est fit from SS", 750, 850);
    c_WJET_fit_SS->SetTopMargin(0.05);
    c_WJET_fit_SS->SetRightMargin(0.05);
    c_WJET_fit_SS->SetBottomMargin(0.15);
    c_WJET_fit_SS->SetLeftMargin(0.17);
    h_WJET_fit_SS->Draw("BAR");
    h_WJET_fit_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]");
    h_WJET_fit_SS->GetXaxis()->SetTitleSize(0.062);
    h_WJET_fit_SS->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_fit_SS->GetXaxis()->SetLabelSize(0.048);
    h_WJET_fit_SS->GetXaxis()->SetMoreLogLabels();
    h_WJET_fit_SS->GetXaxis()->SetNoExponent();
    h_WJET_fit_SS->GetYaxis()->SetTitle("Number of events");
    h_WJET_fit_SS->GetYaxis()->SetTitleSize(0.05);
    h_WJET_fit_SS->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_fit_SS->GetYaxis()->SetLabelSize(0.043);
    h_WJET_fit_SS->GetYaxis()->SetMoreLogLabels();
    h_WJET_fit_SS->GetYaxis()->SetNoExponent();
    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu_TEMPLATE_DOWN.root", "READ");
        TH1D *h_ratio, *h_up, *h_down, *h_fullsysterr, *h_fullerr;
        f_ratio->GetObject("h_WJET_fit_SS", h_ratio);
        f_up->GetObject("h_WJET_fit_SS", h_up);
        f_down->GetObject("h_WJET_fit_SS", h_down);
        h_fullsysterr = ((TH1D*)(h_ratio->Clone("h_WJET_fullsysterr")));
        h_fullerr = ((TH1D*)(h_ratio->Clone("h_WJET_fullerr")));

        for (Int_t i=1; i<h_ratio->GetSize()-1; i++)
        {
            // Systematic errors
            h_ratio->SetBinContent(i, fabs(h_ratio->GetBinContent(i)-h_WJET_fit_SS->GetBinContent(i)));
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_fullsysterr->SetBinContent(i, sqrt(h_ratio->GetBinContent(i)*h_ratio->GetBinContent(i)+h_up->GetBinContent(i)*h_up->GetBinContent(i)));
            // Statistical errors
            if (h_WJET_fit_SS->GetBinContent(i) > 0) h_WJET_fit_SS->SetBinError(i, sqrt(h_WJET_fit_SS->GetBinContent(i)));
            else h_WJET_fit_SS->SetBinError(i, 1);
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_WJET_fit_SS->GetBinError(i)*h_WJET_fit_SS->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
        }
        cout << "Systematic error: " << h_fullsysterr->Integral() << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_WJET_fit_SS->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_WJET_fit_SS->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);
        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);
        f_ratio->Close();
        f_up->Close();
        f_down->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullerr->Write();

        h_WJET_fit_SS->GetYaxis()->SetRangeUser(0, 550);
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
    }// End of if (systErr)
    l_WJET_est->Draw();
    c_WJET_fit_SS->SetLogx();
    c_WJET_fit_SS->SetGridx();
    c_WJET_fit_SS->SetGridy();
    c_WJET_fit_SS->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_Mu.root" << " has been closed successfully." << endl;
    else cout << "FILE " << Dir+"WJETest_Mu.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu.root COULD NOT BE CLOSED!\n" << endl;

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
