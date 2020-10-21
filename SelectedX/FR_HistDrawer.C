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
void E_HistDrawer_alt();
void E_HistDrawer_alt2();
void Mu_HistDrawer(Int_t type);
void Test_HistDrawer(Int_t type);
void Fit_HistDrawer();

void E_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr);
void E_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr);
void Mu_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr);
void Mu_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr);
void EMu_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr);
void EMu_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr);

void E_PR_HistDrawer();
void E_MatrixMethod_HistDrawer(Int_t remNegBins, Int_t systErr);

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

void FR_HistDrawer (TString WhichX = "", Int_t systErr = 0, Int_t type = 2)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    if (type < 1 || type > 2) // 1 -- draw histograms obtained from MakeSelectionForFR.C; 2 -- draw histograms obtained from FRgraphMaker.C
    {
        cout << "Wrong type!" << endl;
        return;
    }
    if (whichX.Contains("EMU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*******    EMu_QCDest_HistDrawer(" << type << ", " << systErr << ")   *******" << endl;
            EMu_QCDest_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*******    EMu_WJETest_HistDrawer(" << type << ", " << systErr << ")   *******" << endl;
            EMu_WJETest_HistDrawer(type, systErr);
        }
        else
        {
            cout << "Please specify by adding QCD or WJets" << endl;
            return;
        }
    }
    else if (whichX.Contains("MU"))
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
    else if (whichX.Contains("E"))
    {
        Xselected++;
        if (whichX.Contains("PR"))
        {
            cout << "\n*************     E_PR_HistDrawer()    *************" << endl;
            E_PR_HistDrawer();
        }
        else if (whichX.Contains("MATRIX") || whichX.Contains("MM"))
        {
            cout << "\n*******     E_MatrixMethod_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            E_MatrixMethod_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("QCD"))
        {
            cout << "\n*******     E_QCDest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            E_QCDest_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*******     E_WJETest_HistDrawer(" << type << ", " << systErr << ")    *******" << endl;
            E_WJETest_HistDrawer(type, systErr);
        }
        else if (whichX.Contains("ALT") && whichX.Contains("2"))
        {
            cout << "\n*******      E_HistDrawer_alt2()      *******" << endl;
            E_HistDrawer_alt2();
        }
        else if (whichX.Contains("ALT"))
        {
            cout << "\n*******      E_HistDrawer_alt()      *******" << endl;
            E_HistDrawer_alt();
        }
        else
        {
            cout << "\n*******      E_HistDrawer(" << type << ")      *******" << endl;
            E_HistDrawer(type);
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
    DYAnalyzer analyzer("Photon_OR");
    FileMgr fm;
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
    THStack *s_TrkIso_barrel_nume = new THStack("s_TrkIso_barrel_nume", "");
    THStack *s_TrkIso_endcap_nume = new THStack("s_TrkIso_endcap_nume", "");
    THStack *s_TrkIso_barrel_deno = new THStack("s_TrkIso_barrel_deno", "");
    THStack *s_TrkIso_endcap_deno = new THStack("s_TrkIso_endcap_deno", "");
    THStack *s_TrkIso_barrel_ctrl = new THStack("s_TrkIso_barrel_ctrl", "");
    THStack *s_TrkIso_endcap_ctrl = new THStack("s_TrkIso_endcap_ctrl", "");
    THStack *s_ECALiso_barrel_nume = new THStack("s_ECALiso_barrel_nume", "");
    THStack *s_ECALiso_endcap_nume = new THStack("s_ECALiso_endcap_nume", "");
    THStack *s_ECALiso_barrel_deno = new THStack("s_ECALiso_barrel_deno", "");
    THStack *s_ECALiso_endcap_deno = new THStack("s_ECALiso_endcap_deno", "");
    THStack *s_ECALiso_barrel_ctrl = new THStack("s_ECALiso_barrel_ctrl", "");
    THStack *s_ECALiso_endcap_ctrl = new THStack("s_ECALiso_endcap_ctrl", "");
    THStack *s_HCALiso_barrel_nume = new THStack("s_HCALiso_barrel_nume", "");
    THStack *s_HCALiso_endcap_nume = new THStack("s_HCALiso_endcap_nume", "");
    THStack *s_HCALiso_barrel_deno = new THStack("s_HCALiso_barrel_deno", "");
    THStack *s_HCALiso_endcap_deno = new THStack("s_HCALiso_endcap_deno", "");
    THStack *s_HCALiso_barrel_ctrl = new THStack("s_HCALiso_barrel_ctrl", "");
    THStack *s_HCALiso_endcap_ctrl = new THStack("s_HCALiso_endcap_ctrl", "");
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
    THStack *s_pT_barrel_deno_density = new THStack("s_pT_barrel_deno_density", "");
    THStack *s_pT_endcap_deno_density = new THStack("s_pT_endcap_deno_density", "");

    THStack *s_HoverE_barrel_template_int = new THStack("s_HoverE_barrel_template_int", "");
    THStack *s_HoverE_barrel_jetTemplate_int = new THStack("s_HoverE_barrel_jetTemplate_int", "");
    THStack *s_HoverE_endcap_template_int = new THStack("s_HoverE_endcap_template_int", "");
    THStack *s_HoverE_endcap_jetTemplate_int = new THStack("s_HoverE_endcap_jetTemplate_int", "");

    THStack *s_PFiso_Rho_barrel_separate = new THStack("s_PFiso_Rho_barrel_separate", "");
    THStack *s_PFiso_Rho_endcap_separate = new THStack("s_PFiso_Rho_endcap_separate", "");
    THStack *s_SigmaIEtaIEta_barrel_separate = new THStack("s_SigmaIEtaIEta_barrel_separate", "");
    THStack *s_SigmaIEtaIEta_endcap_separate = new THStack("s_SigmaIEtaIEta_endcap_separate", "");
    THStack *s_dEtaInSeed_barrel_separate = new THStack("s_dEtaInSeed_barrel_separate", "");
    THStack *s_dEtaInSeed_endcap_separate = new THStack("s_dEtaInSeed_endcap_separate", "");
    THStack *s_dPhiIn_barrel_separate = new THStack("s_dPhiIn_barrel_separate", "");
    THStack *s_dPhiIn_endcap_separate = new THStack("s_dPhiIn_endcap_separate", "");
    THStack *s_HoverE_barrel_separate = new THStack("s_HoverE_barrel_separate", "");
    THStack *s_HoverE_endcap_separate = new THStack("s_HoverE_endcap_separate", "");
    THStack *s_InvEminusInvP_barrel_separate = new THStack("s_InvEminusInvP_barrel_separate", "");
    THStack *s_InvEminusInvP_endcap_separate = new THStack("s_InvEminusInvP_endcap_separate", "");
    THStack *s_mHits_barrel_separate = new THStack("s_mHits_barrel_separate", "");
    THStack *s_mHits_endcap_separate = new THStack("s_mHits_endcap_separate", "");
    THStack *s_passConvVeto_barrel_separate = new THStack("s_passConvVeto_barrel_separate", "");
    THStack *s_passConvVeto_endcap_separate = new THStack("s_passConvVeto_endcap_separate", "");
    THStack *s_passMediumID_barrel_separate = new THStack("s_passMediumID_barrel_separate", "");
    THStack *s_passMediumID_endcap_separate = new THStack("s_passMediumID_endcap_separate", "");

    TH1D *h_PFiso_Rho_barrel_MC_nume[_EndOf_Data_Special], *h_PFiso_Rho_endcap_MC_nume[_EndOf_Data_Special],
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
         *h_TrkIso_barrel_MC_nume[_EndOf_Data_Special], *h_TrkIso_endcap_MC_nume[_EndOf_Data_Special],
         *h_TrkIso_barrel_MC_deno[_EndOf_Data_Special], *h_TrkIso_endcap_MC_deno[_EndOf_Data_Special],
         *h_TrkIso_barrel_MC_ctrl[_EndOf_Data_Special], *h_TrkIso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_ECALiso_barrel_MC_nume[_EndOf_Data_Special], *h_ECALiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_ECALiso_barrel_MC_deno[_EndOf_Data_Special], *h_ECALiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_ECALiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_ECALiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_HCALiso_barrel_MC_nume[_EndOf_Data_Special], *h_HCALiso_endcap_MC_nume[_EndOf_Data_Special],
         *h_HCALiso_barrel_MC_deno[_EndOf_Data_Special], *h_HCALiso_endcap_MC_deno[_EndOf_Data_Special],
         *h_HCALiso_barrel_MC_ctrl[_EndOf_Data_Special], *h_HCALiso_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_MET_MC[_EndOf_Data_Special], *h_eta_data, *h_nVTX_data,
         *h_MT_barrel_MC_nume[_EndOf_Data_Special], *h_MT_endcap_MC_nume[_EndOf_Data_Special],
         *h_MT_barrel_MC_deno[_EndOf_Data_Special], *h_MT_endcap_MC_deno[_EndOf_Data_Special],
         *h_MT_barrel_MC_ctrl[_EndOf_Data_Special], *h_MT_endcap_MC_ctrl[_EndOf_Data_Special],
         *h_eta_MC[_EndOf_Data_Special], *h_nVTX_MC[_EndOf_Data_Special], *h_mass_test_MC[_EndOf_Data_Special],
         *h_pT_barrel_MC_deno_density[_EndOf_Data_Special], *h_pT_endcap_MC_deno_density[_EndOf_Data_Special],            
         *h_PFiso_Rho_barrel_MC_separate[_EndOf_Data_Special], *h_PFiso_Rho_endcap_MC_separate[_EndOf_Data_Special],
         *h_SigmaIEtaIEta_barrel_MC_separate[_EndOf_Data_Special], *h_SigmaIEtaIEta_endcap_MC_separate[_EndOf_Data_Special],
         *h_dEtaInSeed_barrel_MC_separate[_EndOf_Data_Special], *h_dEtaInSeed_endcap_MC_separate[_EndOf_Data_Special],
         *h_dPhiIn_barrel_MC_separate[_EndOf_Data_Special], *h_dPhiIn_endcap_MC_separate[_EndOf_Data_Special],
         *h_HoverE_barrel_MC_separate[_EndOf_Data_Special], *h_HoverE_endcap_MC_separate[_EndOf_Data_Special],
         *h_InvEminusInvP_barrel_MC_separate[_EndOf_Data_Special], *h_InvEminusInvP_endcap_MC_separate[_EndOf_Data_Special],
         *h_mHits_barrel_MC_separate[_EndOf_Data_Special], *h_mHits_endcap_MC_separate[_EndOf_Data_Special],
         *h_passConvVeto_barrel_MC_separate[_EndOf_Data_Special], *h_passConvVeto_endcap_MC_separate[_EndOf_Data_Special],
         *h_passMediumID_barrel_MC_separate[_EndOf_Data_Special], *h_passMediumID_endcap_MC_separate[_EndOf_Data_Special],
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
         *h_TrkIso_barrel_data_nume, *h_TrkIso_endcap_data_nume,
         *h_TrkIso_barrel_data_deno, *h_TrkIso_endcap_data_deno,
         *h_TrkIso_barrel_data_ctrl, *h_TrkIso_endcap_data_ctrl,
         *h_ECALiso_barrel_data_nume, *h_ECALiso_endcap_data_nume,
         *h_ECALiso_barrel_data_deno, *h_ECALiso_endcap_data_deno,
         *h_ECALiso_barrel_data_ctrl, *h_ECALiso_endcap_data_ctrl,
         *h_HCALiso_barrel_data_nume, *h_HCALiso_endcap_data_nume,
         *h_HCALiso_barrel_data_deno, *h_HCALiso_endcap_data_deno,
         *h_HCALiso_barrel_data_ctrl, *h_HCALiso_endcap_data_ctrl,
         *h_MET_data,
         *h_MT_barrel_data_nume, *h_MT_endcap_data_nume,
         *h_MT_barrel_data_deno, *h_MT_endcap_data_deno,
         *h_MT_barrel_data_ctrl, *h_MT_endcap_data_ctrl,
         *h_mass_test_data,
         *h_pT_barrel_data_deno_density, *h_pT_endcap_data_deno_density,
         *h_HoverE_barrel_template_int_MC[_EndOf_Data_Special],    *h_HoverE_endcap_template_int_MC[_EndOf_Data_Special],
         *h_HoverE_barrel_jetTemplate_int_MC[_EndOf_Data_Special], *h_HoverE_endcap_jetTemplate_int_MC[_EndOf_Data_Special],
         *h_HoverE_barrel_template_int_data,    *h_HoverE_endcap_template_int_data,
         *h_HoverE_barrel_jetTemplate_int_data, *h_HoverE_endcap_jetTemplate_int_data,
         *h_PFiso_Rho_barrel_data_separate, *h_PFiso_Rho_endcap_data_separate,
         *h_SigmaIEtaIEta_barrel_data_separate, *h_SigmaIEtaIEta_endcap_data_separate,
         *h_dEtaInSeed_barrel_data_separate, *h_dEtaInSeed_endcap_data_separate,
         *h_dPhiIn_barrel_data_separate, *h_dPhiIn_endcap_data_separate,
         *h_HoverE_barrel_data_separate, *h_HoverE_endcap_data_separate,
         *h_InvEminusInvP_barrel_data_separate, *h_InvEminusInvP_endcap_data_separate,
         *h_mHits_barrel_data_separate, *h_mHits_endcap_data_separate,
         *h_passConvVeto_barrel_data_separate, *h_passConvVeto_endcap_data_separate,
         *h_passMediumID_barrel_data_separate, *h_passMediumID_endcap_data_separate;

    THStack *s_PFiso_Rho_barrel_template[nPtBin_ele], *s_PFiso_Rho_endcap_template[nPtBin_ele],
            *s_PFiso_Rho_barrel_jetTemplate[nPtBin_ele], *s_PFiso_Rho_endcap_jetTemplate[nPtBin_ele];

    cout << "Main hist containers created." << endl;

    TH1D *h_PFiso_Rho_barrel_template[_EndOf_Data_Special][nPtBin_ele], *h_PFiso_Rho_endcap_template[_EndOf_Data_Special][nPtBin_ele],
         *h_PFiso_Rho_barrel_jetTemplate[_EndOf_Data_Special][nPtBin_ele], *h_PFiso_Rho_endcap_jetTemplate[_EndOf_Data_Special][nPtBin_ele],
         *h_PFiso_Rho_barrel_data_template[nPtBin_ele], *h_PFiso_Rho_endcap_data_template[nPtBin_ele],
         *h_PFiso_Rho_barrel_data_jetTemplate[nPtBin_ele], *h_PFiso_Rho_endcap_data_jetTemplate[nPtBin_ele];

    cout << "All hist containers created." << endl;

    for (Int_t ih=0; ih<nPtBin_ele; ih++)
    {
        s_PFiso_Rho_barrel_template[ih] = new THStack("s_PFiso_Rho_barrel_template"+TString::Itoa(ih, 10), "");
        s_PFiso_Rho_barrel_jetTemplate[ih] = new THStack("s_PFiso_Rho_barrel_jetTemplate"+TString::Itoa(ih, 10), "");
        if (ih < nPtBin_ele)
        {
            s_PFiso_Rho_endcap_template[ih] = new THStack("s_PFiso_Rho_endcap_template"+TString::Itoa(ih, 10), "");
            s_PFiso_Rho_endcap_jetTemplate[ih] = new THStack("s_PFiso_Rho_endcap_jetTemplate"+TString::Itoa(ih, 10), "");
        }
    }

    cout << "THStacks created." << endl;

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr1]+".root", "READ");

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
        file->GetObject("h_TrkIso_barrel_nume", h_TrkIso_barrel_MC_nume[pr1]);
        file->GetObject("h_TrkIso_endcap_nume", h_TrkIso_endcap_MC_nume[pr1]);
        file->GetObject("h_TrkIso_barrel_deno", h_TrkIso_barrel_MC_deno[pr1]);
        file->GetObject("h_TrkIso_endcap_deno", h_TrkIso_endcap_MC_deno[pr1]);
        file->GetObject("h_TrkIso_barrel_ctrl", h_TrkIso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_TrkIso_endcap_ctrl", h_TrkIso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_ECALiso_barrel_nume", h_ECALiso_barrel_MC_nume[pr1]);
        file->GetObject("h_ECALiso_endcap_nume", h_ECALiso_endcap_MC_nume[pr1]);
        file->GetObject("h_ECALiso_barrel_deno", h_ECALiso_barrel_MC_deno[pr1]);
        file->GetObject("h_ECALiso_endcap_deno", h_ECALiso_endcap_MC_deno[pr1]);
        file->GetObject("h_ECALiso_barrel_ctrl", h_ECALiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_ECALiso_endcap_ctrl", h_ECALiso_endcap_MC_ctrl[pr1]);
        file->GetObject("h_HCALiso_barrel_nume", h_HCALiso_barrel_MC_nume[pr1]);
        file->GetObject("h_HCALiso_endcap_nume", h_HCALiso_endcap_MC_nume[pr1]);
        file->GetObject("h_HCALiso_barrel_deno", h_HCALiso_barrel_MC_deno[pr1]);
        file->GetObject("h_HCALiso_endcap_deno", h_HCALiso_endcap_MC_deno[pr1]);
        file->GetObject("h_HCALiso_barrel_ctrl", h_HCALiso_barrel_MC_ctrl[pr1]);
        file->GetObject("h_HCALiso_endcap_ctrl", h_HCALiso_endcap_MC_ctrl[pr1]);
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
        file->GetObject("h_HoverE_barrel_template_int", h_HoverE_barrel_template_int_MC[pr1]);
        file->GetObject("h_HoverE_barrel_jetTemplate_int", h_HoverE_barrel_jetTemplate_int_MC[pr1]);
        file->GetObject("h_HoverE_endcap_template_int", h_HoverE_endcap_template_int_MC[pr1]);
        file->GetObject("h_HoverE_endcap_jetTemplate_int", h_HoverE_endcap_jetTemplate_int_MC[pr1]);
        file->GetObject("h_PFiso_Rho_barrel_separate", h_PFiso_Rho_barrel_MC_separate[pr1]);
        file->GetObject("h_PFiso_Rho_endcap_separate", h_PFiso_Rho_endcap_MC_separate[pr1]);
        file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_SigmaIEtaIEta_barrel_MC_separate[pr1]);
        file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_SigmaIEtaIEta_endcap_MC_separate[pr1]);
        file->GetObject("h_dEtaInSeed_barrel_separate", h_dEtaInSeed_barrel_MC_separate[pr1]);
        file->GetObject("h_dEtaInSeed_endcap_separate", h_dEtaInSeed_endcap_MC_separate[pr1]);
        file->GetObject("h_dPhiIn_barrel_separate", h_dPhiIn_barrel_MC_separate[pr1]);
        file->GetObject("h_dPhiIn_endcap_separate", h_dPhiIn_endcap_MC_separate[pr1]);
        file->GetObject("h_HoverE_barrel_separate", h_HoverE_barrel_MC_separate[pr1]);
        file->GetObject("h_HoverE_endcap_separate", h_HoverE_endcap_MC_separate[pr1]);
        file->GetObject("h_InvEminusInvP_barrel_separate", h_InvEminusInvP_barrel_MC_separate[pr1]);
        file->GetObject("h_InvEminusInvP_endcap_separate", h_InvEminusInvP_endcap_MC_separate[pr1]);
        file->GetObject("h_mHits_barrel_separate", h_mHits_barrel_MC_separate[pr1]);
        file->GetObject("h_mHits_endcap_separate", h_mHits_endcap_MC_separate[pr1]);
        file->GetObject("h_passConvVeto_barrel_separate", h_passConvVeto_barrel_MC_separate[pr1]);
        file->GetObject("h_passConvVeto_endcap_separate", h_passConvVeto_endcap_MC_separate[pr1]);
        file->GetObject("h_passMediumID_barrel_separate", h_passMediumID_barrel_MC_separate[pr1]);
        file->GetObject("h_passMediumID_endcap_separate", h_passMediumID_endcap_MC_separate[pr1]);
        h_pT_barrel_MC_deno_density[pr1] = ((TH1D*)(h_pT_barrel_MC_deno[pr1]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr1])));
        h_pT_endcap_MC_deno_density[pr1] = ((TH1D*)(h_pT_endcap_MC_deno[pr1]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr1])));

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
        removeNegativeBins(h_TrkIso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_TrkIso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_TrkIso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_TrkIso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_TrkIso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_TrkIso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_ECALiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_ECALiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_ECALiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_ECALiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_ECALiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_ECALiso_endcap_MC_ctrl[pr1]);
        removeNegativeBins(h_HCALiso_barrel_MC_nume[pr1]);
        removeNegativeBins(h_HCALiso_endcap_MC_nume[pr1]);
        removeNegativeBins(h_HCALiso_barrel_MC_deno[pr1]);
        removeNegativeBins(h_HCALiso_endcap_MC_deno[pr1]);
        removeNegativeBins(h_HCALiso_barrel_MC_ctrl[pr1]);
        removeNegativeBins(h_HCALiso_endcap_MC_ctrl[pr1]);
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
        removeNegativeBins(h_HoverE_barrel_template_int_MC[pr1]);
        removeNegativeBins(h_HoverE_barrel_jetTemplate_int_MC[pr1]);
        removeNegativeBins(h_HoverE_endcap_template_int_MC[pr1]);
        removeNegativeBins(h_HoverE_endcap_jetTemplate_int_MC[pr1]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_separate[pr1]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_separate[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_separate[pr1]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_separate[pr1]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_separate[pr1]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_separate[pr1]);
        removeNegativeBins(h_dPhiIn_barrel_MC_separate[pr1]);
        removeNegativeBins(h_dPhiIn_endcap_MC_separate[pr1]);
        removeNegativeBins(h_HoverE_barrel_MC_separate[pr1]);
        removeNegativeBins(h_HoverE_endcap_MC_separate[pr1]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_separate[pr1]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_separate[pr1]);
        removeNegativeBins(h_mHits_barrel_MC_separate[pr1]);
        removeNegativeBins(h_mHits_endcap_MC_separate[pr1]);
        removeNegativeBins(h_passConvVeto_barrel_MC_separate[pr1]);
        removeNegativeBins(h_passConvVeto_endcap_MC_separate[pr1]);
        removeNegativeBins(h_passMediumID_barrel_MC_separate[pr1]);
        removeNegativeBins(h_passMediumID_endcap_MC_separate[pr1]);

        removeNegativeBins(h_pT_barrel_MC_deno_density[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr1]);

        for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr1]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr1]->GetBinContent(i_bin) /
                                                                  (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr1]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr1]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr1]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr1]->GetBinContent(i_bin) /
                                                                  (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr1]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr1]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
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
        h_TrkIso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_TrkIso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_TrkIso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_TrkIso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_TrkIso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_TrkIso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_ECALiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_ECALiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_ECALiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_ECALiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_ECALiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_ECALiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
        h_HCALiso_barrel_MC_nume[pr1]->SetFillColor(color);
        h_HCALiso_endcap_MC_nume[pr1]->SetFillColor(color);
        h_HCALiso_barrel_MC_deno[pr1]->SetFillColor(color);
        h_HCALiso_endcap_MC_deno[pr1]->SetFillColor(color);
        h_HCALiso_barrel_MC_ctrl[pr1]->SetFillColor(color);
        h_HCALiso_endcap_MC_ctrl[pr1]->SetFillColor(color);
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
        h_HoverE_barrel_template_int_MC[pr1]->SetFillColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr1]->SetFillColor(color);
        h_HoverE_endcap_template_int_MC[pr1]->SetFillColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr1]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr1]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr1]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr1]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr1]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr1]->SetFillColor(color);
        h_dPhiIn_barrel_MC_separate[pr1]->SetFillColor(color);
        h_dPhiIn_endcap_MC_separate[pr1]->SetFillColor(color);
        h_HoverE_barrel_MC_separate[pr1]->SetFillColor(color);
        h_HoverE_endcap_MC_separate[pr1]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr1]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr1]->SetFillColor(color);
        h_mHits_barrel_MC_separate[pr1]->SetFillColor(color);
        h_mHits_endcap_MC_separate[pr1]->SetFillColor(color);
        h_passConvVeto_barrel_MC_separate[pr1]->SetFillColor(color);
        h_passConvVeto_endcap_MC_separate[pr1]->SetFillColor(color);
        h_passMediumID_barrel_MC_separate[pr1]->SetFillColor(color);
        h_passMediumID_endcap_MC_separate[pr1]->SetFillColor(color);
        h_pT_barrel_MC_deno_density[pr1]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr1]->SetFillColor(color);

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
        h_TrkIso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_TrkIso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_TrkIso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_TrkIso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_TrkIso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_TrkIso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_ECALiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_ECALiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_ECALiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_ECALiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_ECALiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_ECALiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
        h_HCALiso_barrel_MC_nume[pr1]->SetLineColor(color);
        h_HCALiso_endcap_MC_nume[pr1]->SetLineColor(color);
        h_HCALiso_barrel_MC_deno[pr1]->SetLineColor(color);
        h_HCALiso_endcap_MC_deno[pr1]->SetLineColor(color);
        h_HCALiso_barrel_MC_ctrl[pr1]->SetLineColor(color);
        h_HCALiso_endcap_MC_ctrl[pr1]->SetLineColor(color);
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
        h_HoverE_barrel_template_int_MC[pr1]->SetLineColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr1]->SetLineColor(color);
        h_HoverE_endcap_template_int_MC[pr1]->SetLineColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr1]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr1]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr1]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr1]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr1]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr1]->SetLineColor(color);
        h_dPhiIn_barrel_MC_separate[pr1]->SetLineColor(color);
        h_dPhiIn_endcap_MC_separate[pr1]->SetLineColor(color);
        h_HoverE_barrel_MC_separate[pr1]->SetLineColor(color);
        h_HoverE_endcap_MC_separate[pr1]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr1]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr1]->SetLineColor(color);
        h_mHits_barrel_MC_separate[pr1]->SetLineColor(color);
        h_mHits_endcap_MC_separate[pr1]->SetLineColor(color);
        h_passConvVeto_barrel_MC_separate[pr1]->SetLineColor(color);
        h_passConvVeto_endcap_MC_separate[pr1]->SetLineColor(color);
        h_passMediumID_barrel_MC_separate[pr1]->SetLineColor(color);
        h_passMediumID_endcap_MC_separate[pr1]->SetLineColor(color);
        h_pT_barrel_MC_deno_density[pr1]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr1]->SetLineColor(color);

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
        h_TrkIso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_TrkIso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_TrkIso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_TrkIso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_TrkIso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_TrkIso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_ECALiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_ECALiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_ECALiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_ECALiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_ECALiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_ECALiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
        h_HCALiso_barrel_MC_nume[pr1]->SetDirectory(0);
        h_HCALiso_endcap_MC_nume[pr1]->SetDirectory(0);
        h_HCALiso_barrel_MC_deno[pr1]->SetDirectory(0);
        h_HCALiso_endcap_MC_deno[pr1]->SetDirectory(0);
        h_HCALiso_barrel_MC_ctrl[pr1]->SetDirectory(0);
        h_HCALiso_endcap_MC_ctrl[pr1]->SetDirectory(0);
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
        h_HoverE_barrel_template_int_MC[pr1]->SetDirectory(0);
        h_HoverE_barrel_jetTemplate_int_MC[pr1]->SetDirectory(0);
        h_HoverE_endcap_template_int_MC[pr1]->SetDirectory(0);
        h_HoverE_endcap_jetTemplate_int_MC[pr1]->SetDirectory(0);        
        h_PFiso_Rho_barrel_MC_separate[pr1]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_separate[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_separate[pr1]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_separate[pr1]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_separate[pr1]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_separate[pr1]->SetDirectory(0);
        h_dPhiIn_barrel_MC_separate[pr1]->SetDirectory(0);
        h_dPhiIn_endcap_MC_separate[pr1]->SetDirectory(0);
        h_HoverE_barrel_MC_separate[pr1]->SetDirectory(0);
        h_HoverE_endcap_MC_separate[pr1]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_separate[pr1]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_separate[pr1]->SetDirectory(0);
        h_mHits_barrel_MC_separate[pr1]->SetDirectory(0);
        h_mHits_endcap_MC_separate[pr1]->SetDirectory(0);
        h_passConvVeto_barrel_MC_separate[pr1]->SetDirectory(0);
        h_passConvVeto_endcap_MC_separate[pr1]->SetDirectory(0);
        h_passMediumID_barrel_MC_separate[pr1]->SetDirectory(0);
        h_passMediumID_endcap_MC_separate[pr1]->SetDirectory(0);
        h_pT_barrel_MC_deno_density[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr1]->SetDirectory(0);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_template[pr1][ih]);
            file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_jetTemplate[pr1][ih]);

            removeNegativeBins(h_PFiso_Rho_barrel_template[pr1][ih]);
            removeNegativeBins(h_PFiso_Rho_barrel_jetTemplate[pr1][ih]);

            h_PFiso_Rho_barrel_template[pr1][ih]   ->SetFillColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr1][ih]->SetFillColor(color);

            h_PFiso_Rho_barrel_template[pr1][ih]   ->SetLineColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr1][ih]->SetLineColor(color);

            h_PFiso_Rho_barrel_template[pr1][ih]   ->SetDirectory(0);
            h_PFiso_Rho_barrel_jetTemplate[pr1][ih]->SetDirectory(0);

            if (ih < nPtBin_ele)
            {
                file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_template[pr1][ih]);
                file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_jetTemplate[pr1][ih]);

                removeNegativeBins(h_PFiso_Rho_endcap_template[pr1][ih]);
                removeNegativeBins(h_PFiso_Rho_endcap_jetTemplate[pr1][ih]);

                h_PFiso_Rho_endcap_template[pr1][ih]   ->SetFillColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr1][ih]->SetFillColor(color);

                h_PFiso_Rho_endcap_template[pr1][ih]   ->SetLineColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr1][ih]->SetLineColor(color);

                h_PFiso_Rho_endcap_template[pr1][ih]   ->SetDirectory(0);
                h_PFiso_Rho_endcap_jetTemplate[pr1][ih]->SetDirectory(0);
            }
        }

        if (pr1 == _WJets)
        {
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
            h_TrkIso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_nume[pr1]->Clone("h_TrkIso_barrel_nume_WJets")));
            h_TrkIso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_nume[pr1]->Clone("h_TrkIso_endcap_nume_WJets")));
            h_TrkIso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_deno[pr1]->Clone("h_TrkIso_barrel_deno_WJets")));
            h_TrkIso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_deno[pr1]->Clone("h_TrkIso_endcap_deno_WJets")));
            h_TrkIso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_ctrl[pr1]->Clone("h_TrkIso_barrel_ctrl_WJets")));
            h_TrkIso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_ctrl[pr1]->Clone("h_TrkIso_endcap_ctrl_WJets")));
            h_ECALiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_nume[pr1]->Clone("h_ECALiso_barrel_nume_WJets")));
            h_ECALiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_nume[pr1]->Clone("h_ECALiso_endcap_nume_WJets")));
            h_ECALiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_deno[pr1]->Clone("h_ECALiso_barrel_deno_WJets")));
            h_ECALiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_deno[pr1]->Clone("h_ECALiso_endcap_deno_WJets")));
            h_ECALiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_ctrl[pr1]->Clone("h_ECALiso_barrel_ctrl_WJets")));
            h_ECALiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_ctrl[pr1]->Clone("h_ECALiso_endcap_ctrl_WJets")));
            h_HCALiso_barrel_MC_nume[_WJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_nume[pr1]->Clone("h_HCALiso_barrel_nume_WJets")));
            h_HCALiso_endcap_MC_nume[_WJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_nume[pr1]->Clone("h_HCALiso_endcap_nume_WJets")));
            h_HCALiso_barrel_MC_deno[_WJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_deno[pr1]->Clone("h_HCALiso_barrel_deno_WJets")));
            h_HCALiso_endcap_MC_deno[_WJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_deno[pr1]->Clone("h_HCALiso_endcap_deno_WJets")));
            h_HCALiso_barrel_MC_ctrl[_WJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_ctrl[pr1]->Clone("h_HCALiso_barrel_ctrl_WJets")));
            h_HCALiso_endcap_MC_ctrl[_WJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_ctrl[pr1]->Clone("h_HCALiso_endcap_ctrl_WJets")));
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
            h_HoverE_barrel_template_int_MC[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_template_int_MC[pr1]->Clone("h_HoverE_barrel_template_int_WJets")));
            h_HoverE_barrel_jetTemplate_int_MC[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_jetTemplate_int_MC[pr1]->Clone("h_HoverE_barrel_jetTemplate_int_WJets")));
            h_HoverE_endcap_template_int_MC[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_template_int_MC[pr1]->Clone("h_HoverE_endcap_template_int_WJets")));
            h_HoverE_endcap_jetTemplate_int_MC[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_jetTemplate_int_MC[pr1]->Clone("h_HoverE_endcap_jetTemplate_int_WJets")));           
            h_PFiso_Rho_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_separate[pr1]->Clone("h_PFiso_Rho_barrel_separate_WJets")));
            h_PFiso_Rho_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_separate[pr1]->Clone("h_PFiso_Rho_endcap_separate_WJets")));
            h_SigmaIEtaIEta_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_separate[pr1]->Clone("h_SigmaIEtaIEta_barrel_separate_WJets")));
            h_SigmaIEtaIEta_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_separate[pr1]->Clone("h_SigmaIEtaIEta_endcap_separate_WJets")));
            h_dEtaInSeed_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_separate[pr1]->Clone("h_dEtaInSeed_barrel_separate_WJets")));
            h_dEtaInSeed_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_separate[pr1]->Clone("h_dEtaInSeed_endcap_separate_WJets")));
            h_dPhiIn_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_separate[pr1]->Clone("h_dPhiIn_barrel_separate_WJets")));
            h_dPhiIn_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_separate[pr1]->Clone("h_dPhiIn_endcap_separate_WJets")));
            h_HoverE_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_separate[pr1]->Clone("h_HoverE_barrel_separate_WJets")));
            h_HoverE_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_separate[pr1]->Clone("h_HoverE_endcap_separate_WJets")));
            h_InvEminusInvP_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_separate[pr1]->Clone("h_InvEminusInvP_barrel_separate_WJets")));
            h_InvEminusInvP_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_separate[pr1]->Clone("h_InvEminusInvP_endcap_separate_WJets")));
            h_mHits_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_mHits_barrel_MC_separate[pr1]->Clone("h_mHits_barrel_MC_separate_WJets")));
            h_mHits_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_mHits_endcap_MC_separate[pr1]->Clone("h_mHits_endcap_MC_separate_WJets")));
            h_passConvVeto_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_passConvVeto_barrel_MC_separate[pr1]->Clone("h_passConvVeto_barrel_MC_separate_WJets")));
            h_passConvVeto_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_passConvVeto_endcap_MC_separate[pr1]->Clone("h_passConvVeto_endcap_MC_separate_WJets")));
            h_passMediumID_barrel_MC_separate[_WJets_Full] = ((TH1D*)(h_passMediumID_barrel_MC_separate[pr1]->Clone("h_passMediumID_barrel_MC_separate_WJets")));
            h_passMediumID_endcap_MC_separate[_WJets_Full] = ((TH1D*)(h_passMediumID_endcap_MC_separate[pr1]->Clone("h_passMediumID_endcap_MC_separate_WJets")));
            h_pT_barrel_MC_deno_density[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr1]->Clone("h_pT_barrel_deno_density_WJets_Full")));
            h_pT_endcap_MC_deno_density[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr1]->Clone("h_pT_endcap_deno_density_WJets_Full")));

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
            h_TrkIso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_nume[_WJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_nume[_WJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_deno[_WJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_deno[_WJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_ctrl[_WJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_ctrl[_WJets_Full]->SetDirectory(0);
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
            h_HoverE_barrel_template_int_MC[_WJets_Full]->SetDirectory(0);
            h_HoverE_barrel_jetTemplate_int_MC[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_template_int_MC[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_jetTemplate_int_MC[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_mHits_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_mHits_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_passConvVeto_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_passConvVeto_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_passMediumID_barrel_MC_separate[_WJets_Full]->SetDirectory(0);
            h_passMediumID_endcap_MC_separate[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno_density[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_WJets_Full]->SetDirectory(0);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_WJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_template[pr1][ih]->Clone("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10)+"_WJets")));
                h_PFiso_Rho_barrel_jetTemplate[_WJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_jetTemplate[pr1][ih]->Clone("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10)+"_WJets")));
                h_PFiso_Rho_barrel_template[_WJets_Full][ih]->SetDirectory(0);
                h_PFiso_Rho_barrel_jetTemplate[_WJets_Full][ih]->SetDirectory(0);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_WJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_template[pr1][ih]->Clone("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10)+"_WJets")));
                    h_PFiso_Rho_endcap_jetTemplate[_WJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_jetTemplate[pr1][ih]->Clone("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10)+"_WJets")));
                    h_PFiso_Rho_endcap_template[_WJets_Full][ih]->SetDirectory(0);
                    h_PFiso_Rho_endcap_jetTemplate[_WJets_Full][ih]->SetDirectory(0);
                }
            }
        }
        else if (pr1 == _WJets_ext2v5)
        {
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
            h_TrkIso_barrel_MC_nume[_WJets_Full]->Add(h_TrkIso_barrel_MC_nume[pr1]);
            h_TrkIso_endcap_MC_nume[_WJets_Full]->Add(h_TrkIso_endcap_MC_nume[pr1]);
            h_TrkIso_barrel_MC_deno[_WJets_Full]->Add(h_TrkIso_barrel_MC_deno[pr1]);
            h_TrkIso_endcap_MC_deno[_WJets_Full]->Add(h_TrkIso_endcap_MC_deno[pr1]);
            h_TrkIso_barrel_MC_ctrl[_WJets_Full]->Add(h_TrkIso_barrel_MC_ctrl[pr1]);
            h_TrkIso_endcap_MC_ctrl[_WJets_Full]->Add(h_TrkIso_endcap_MC_ctrl[pr1]);
            h_ECALiso_barrel_MC_nume[_WJets_Full]->Add(h_ECALiso_barrel_MC_nume[pr1]);
            h_ECALiso_endcap_MC_nume[_WJets_Full]->Add(h_ECALiso_endcap_MC_nume[pr1]);
            h_ECALiso_barrel_MC_deno[_WJets_Full]->Add(h_ECALiso_barrel_MC_deno[pr1]);
            h_ECALiso_endcap_MC_deno[_WJets_Full]->Add(h_ECALiso_endcap_MC_deno[pr1]);
            h_ECALiso_barrel_MC_ctrl[_WJets_Full]->Add(h_ECALiso_barrel_MC_ctrl[pr1]);
            h_ECALiso_endcap_MC_ctrl[_WJets_Full]->Add(h_ECALiso_endcap_MC_ctrl[pr1]);
            h_HCALiso_barrel_MC_nume[_WJets_Full]->Add(h_HCALiso_barrel_MC_nume[pr1]);
            h_HCALiso_endcap_MC_nume[_WJets_Full]->Add(h_HCALiso_endcap_MC_nume[pr1]);
            h_HCALiso_barrel_MC_deno[_WJets_Full]->Add(h_HCALiso_barrel_MC_deno[pr1]);
            h_HCALiso_endcap_MC_deno[_WJets_Full]->Add(h_HCALiso_endcap_MC_deno[pr1]);
            h_HCALiso_barrel_MC_ctrl[_WJets_Full]->Add(h_HCALiso_barrel_MC_ctrl[pr1]);
            h_HCALiso_endcap_MC_ctrl[_WJets_Full]->Add(h_HCALiso_endcap_MC_ctrl[pr1]);
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
            h_HoverE_barrel_template_int_MC[_WJets_Full]->Add(h_HoverE_barrel_template_int_MC[pr1]);
            h_HoverE_barrel_jetTemplate_int_MC[_WJets_Full]->Add(h_HoverE_barrel_jetTemplate_int_MC[pr1]);
            h_HoverE_endcap_template_int_MC[_WJets_Full]->Add(h_HoverE_endcap_template_int_MC[pr1]);
            h_HoverE_endcap_jetTemplate_int_MC[_WJets_Full]->Add(h_HoverE_endcap_jetTemplate_int_MC[pr1]);
            h_PFiso_Rho_barrel_MC_separate[_WJets_Full]->Add(h_PFiso_Rho_barrel_MC_separate[pr1]);
            h_PFiso_Rho_endcap_MC_separate[_WJets_Full]->Add(h_PFiso_Rho_endcap_MC_separate[pr1]);
            h_SigmaIEtaIEta_barrel_MC_separate[_WJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr1]);
            h_SigmaIEtaIEta_endcap_MC_separate[_WJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr1]);
            h_dEtaInSeed_barrel_MC_separate[_WJets_Full]->Add(h_dEtaInSeed_barrel_MC_separate[pr1]);
            h_dEtaInSeed_endcap_MC_separate[_WJets_Full]->Add(h_dEtaInSeed_endcap_MC_separate[pr1]);
            h_dPhiIn_barrel_MC_separate[_WJets_Full]->Add(h_dPhiIn_barrel_MC_separate[pr1]);
            h_dPhiIn_endcap_MC_separate[_WJets_Full]->Add(h_dPhiIn_endcap_MC_separate[pr1]);
            h_HoverE_barrel_MC_separate[_WJets_Full]->Add(h_HoverE_barrel_MC_separate[pr1]);
            h_HoverE_endcap_MC_separate[_WJets_Full]->Add(h_HoverE_endcap_MC_separate[pr1]);
            h_InvEminusInvP_barrel_MC_separate[_WJets_Full]->Add(h_InvEminusInvP_barrel_MC_separate[pr1]);
            h_InvEminusInvP_endcap_MC_separate[_WJets_Full]->Add(h_InvEminusInvP_endcap_MC_separate[pr1]);
            h_mHits_barrel_MC_separate[_WJets_Full]->Add(h_mHits_barrel_MC_separate[pr1]);
            h_mHits_endcap_MC_separate[_WJets_Full]->Add(h_mHits_endcap_MC_separate[pr1]);
            h_passConvVeto_barrel_MC_separate[_WJets_Full]->Add(h_passConvVeto_barrel_MC_separate[pr1]);
            h_passConvVeto_endcap_MC_separate[_WJets_Full]->Add(h_passConvVeto_endcap_MC_separate[pr1]);
            h_passMediumID_barrel_MC_separate[_WJets_Full]->Add(h_passMediumID_barrel_MC_separate[pr1]);
            h_passMediumID_endcap_MC_separate[_WJets_Full]->Add(h_passMediumID_endcap_MC_separate[pr1]);
            h_pT_barrel_MC_deno_density[_WJets_Full]->Add(h_pT_barrel_MC_deno_density[pr1]);
            h_pT_endcap_MC_deno_density[_WJets_Full]->Add(h_pT_endcap_MC_deno_density[pr1]);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_WJets_Full][ih]->Add(h_PFiso_Rho_barrel_template[pr1][ih]);
                h_PFiso_Rho_barrel_jetTemplate[_WJets_Full][ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr1][ih]);
                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_WJets_Full][ih]->Add(h_PFiso_Rho_endcap_template[pr1][ih]);
                    h_PFiso_Rho_endcap_jetTemplate[_WJets_Full][ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr1][ih]);
                }
            }
        }

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
        s_TrkIso_barrel_nume->Add(h_TrkIso_barrel_MC_nume[pr1]);
        s_TrkIso_endcap_nume->Add(h_TrkIso_endcap_MC_nume[pr1]);
        s_TrkIso_barrel_deno->Add(h_TrkIso_barrel_MC_deno[pr1]);
        s_TrkIso_endcap_deno->Add(h_TrkIso_endcap_MC_deno[pr1]);
        s_TrkIso_barrel_ctrl->Add(h_TrkIso_barrel_MC_ctrl[pr1]);
        s_TrkIso_endcap_ctrl->Add(h_TrkIso_endcap_MC_ctrl[pr1]);
        s_ECALiso_barrel_nume->Add(h_ECALiso_barrel_MC_nume[pr1]);
        s_ECALiso_endcap_nume->Add(h_ECALiso_endcap_MC_nume[pr1]);
        s_ECALiso_barrel_deno->Add(h_ECALiso_barrel_MC_deno[pr1]);
        s_ECALiso_endcap_deno->Add(h_ECALiso_endcap_MC_deno[pr1]);
        s_ECALiso_barrel_ctrl->Add(h_ECALiso_barrel_MC_ctrl[pr1]);
        s_ECALiso_endcap_ctrl->Add(h_ECALiso_endcap_MC_ctrl[pr1]);
        s_HCALiso_barrel_nume->Add(h_HCALiso_barrel_MC_nume[pr1]);
        s_HCALiso_endcap_nume->Add(h_HCALiso_endcap_MC_nume[pr1]);
        s_HCALiso_barrel_deno->Add(h_HCALiso_barrel_MC_deno[pr1]);
        s_HCALiso_endcap_deno->Add(h_HCALiso_endcap_MC_deno[pr1]);
        s_HCALiso_barrel_ctrl->Add(h_HCALiso_barrel_MC_ctrl[pr1]);
        s_HCALiso_endcap_ctrl->Add(h_HCALiso_endcap_MC_ctrl[pr1]);
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
        s_HoverE_barrel_template_int->Add(h_HoverE_barrel_template_int_MC[pr1]);
        s_HoverE_barrel_jetTemplate_int->Add(h_HoverE_barrel_jetTemplate_int_MC[pr1]);
        s_HoverE_endcap_template_int->Add(h_HoverE_endcap_template_int_MC[pr1]);
        s_HoverE_endcap_jetTemplate_int->Add(h_HoverE_endcap_jetTemplate_int_MC[pr1]);
        s_PFiso_Rho_barrel_separate->Add(h_PFiso_Rho_barrel_MC_separate[pr1]);
        s_PFiso_Rho_endcap_separate->Add(h_PFiso_Rho_endcap_MC_separate[pr1]);
        s_SigmaIEtaIEta_barrel_separate->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr1]);
        s_SigmaIEtaIEta_endcap_separate->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr1]);
        s_dEtaInSeed_barrel_separate->Add(h_dEtaInSeed_barrel_MC_separate[pr1]);
        s_dEtaInSeed_endcap_separate->Add(h_dEtaInSeed_endcap_MC_separate[pr1]);
        s_dPhiIn_barrel_separate->Add(h_dPhiIn_barrel_MC_separate[pr1]);
        s_dPhiIn_endcap_separate->Add(h_dPhiIn_endcap_MC_separate[pr1]);
        s_HoverE_barrel_separate->Add(h_HoverE_barrel_MC_separate[pr1]);
        s_HoverE_endcap_separate->Add(h_HoverE_endcap_MC_separate[pr1]);
        s_InvEminusInvP_barrel_separate->Add(h_InvEminusInvP_barrel_MC_separate[pr1]);
        s_InvEminusInvP_endcap_separate->Add(h_InvEminusInvP_endcap_MC_separate[pr1]);
        s_mHits_barrel_separate->Add(h_mHits_barrel_MC_separate[pr1]);
        s_mHits_endcap_separate->Add(h_mHits_endcap_MC_separate[pr1]);
        s_passConvVeto_barrel_separate->Add(h_passConvVeto_barrel_MC_separate[pr1]);
        s_passConvVeto_endcap_separate->Add(h_passConvVeto_endcap_MC_separate[pr1]);
        s_passMediumID_barrel_separate->Add(h_passMediumID_barrel_MC_separate[pr1]);
        s_passMediumID_endcap_separate->Add(h_passMediumID_endcap_MC_separate[pr1]);
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr1]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr1]);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            s_PFiso_Rho_barrel_template[ih]->Add(h_PFiso_Rho_barrel_template[pr1][ih]);
            s_PFiso_Rho_barrel_jetTemplate[ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr1][ih]);
            if (ih < nPtBin_ele)
            {
                s_PFiso_Rho_endcap_template[ih]->Add(h_PFiso_Rho_endcap_template[pr1][ih]);
                s_PFiso_Rho_endcap_jetTemplate[ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr1][ih]);
            }
        }

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

    cout << "Main bkg received" << endl;

    // Drell-Yan
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

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
        file->GetObject("h_TrkIso_barrel_nume", h_TrkIso_barrel_MC_nume[pr]);
        file->GetObject("h_TrkIso_endcap_nume", h_TrkIso_endcap_MC_nume[pr]);
        file->GetObject("h_TrkIso_barrel_deno", h_TrkIso_barrel_MC_deno[pr]);
        file->GetObject("h_TrkIso_endcap_deno", h_TrkIso_endcap_MC_deno[pr]);
        file->GetObject("h_TrkIso_barrel_ctrl", h_TrkIso_barrel_MC_ctrl[pr]);
        file->GetObject("h_TrkIso_endcap_ctrl", h_TrkIso_endcap_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_barrel_nume", h_ECALiso_barrel_MC_nume[pr]);
        file->GetObject("h_ECALiso_endcap_nume", h_ECALiso_endcap_MC_nume[pr]);
        file->GetObject("h_ECALiso_barrel_deno", h_ECALiso_barrel_MC_deno[pr]);
        file->GetObject("h_ECALiso_endcap_deno", h_ECALiso_endcap_MC_deno[pr]);
        file->GetObject("h_ECALiso_barrel_ctrl", h_ECALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_endcap_ctrl", h_ECALiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_barrel_nume", h_HCALiso_barrel_MC_nume[pr]);
        file->GetObject("h_HCALiso_endcap_nume", h_HCALiso_endcap_MC_nume[pr]);
        file->GetObject("h_HCALiso_barrel_deno", h_HCALiso_barrel_MC_deno[pr]);
        file->GetObject("h_HCALiso_endcap_deno", h_HCALiso_endcap_MC_deno[pr]);
        file->GetObject("h_HCALiso_barrel_ctrl", h_HCALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_endcap_ctrl", h_HCALiso_endcap_MC_ctrl[pr]);
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
        file->GetObject("h_HoverE_barrel_template_int", h_HoverE_barrel_template_int_MC[pr]);
        file->GetObject("h_HoverE_barrel_jetTemplate_int", h_HoverE_barrel_jetTemplate_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_template_int", h_HoverE_endcap_template_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_jetTemplate_int", h_HoverE_endcap_jetTemplate_int_MC[pr]);
        file->GetObject("h_PFiso_Rho_barrel_separate", h_PFiso_Rho_barrel_MC_separate[pr]);
        file->GetObject("h_PFiso_Rho_endcap_separate", h_PFiso_Rho_endcap_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_barrel_separate", h_dEtaInSeed_barrel_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_endcap_separate", h_dEtaInSeed_endcap_MC_separate[pr]);
        file->GetObject("h_dPhiIn_barrel_separate", h_dPhiIn_barrel_MC_separate[pr]);
        file->GetObject("h_dPhiIn_endcap_separate", h_dPhiIn_endcap_MC_separate[pr]);
        file->GetObject("h_HoverE_barrel_separate", h_HoverE_barrel_MC_separate[pr]);
        file->GetObject("h_HoverE_endcap_separate", h_HoverE_endcap_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_barrel_separate", h_InvEminusInvP_barrel_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_endcap_separate", h_InvEminusInvP_endcap_MC_separate[pr]);
        file->GetObject("h_mHits_barrel_separate", h_mHits_barrel_MC_separate[pr]);
        file->GetObject("h_mHits_endcap_separate", h_mHits_endcap_MC_separate[pr]);
        file->GetObject("h_passConvVeto_barrel_separate", h_passConvVeto_barrel_MC_separate[pr]);
        file->GetObject("h_passConvVeto_endcap_separate", h_passConvVeto_endcap_MC_separate[pr]);
        file->GetObject("h_passMediumID_barrel_separate", h_passMediumID_barrel_MC_separate[pr]);
        file->GetObject("h_passMediumID_endcap_separate", h_passMediumID_endcap_MC_separate[pr]);
        h_pT_barrel_MC_deno_density[pr] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr])));
        h_pT_endcap_MC_deno_density[pr] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr])));

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
        removeNegativeBins(h_TrkIso_barrel_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_ctrl[pr]);
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
        removeNegativeBins(h_HoverE_barrel_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_barrel_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_separate[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_separate[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_separate[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_separate[pr]);
        removeNegativeBins(h_mHits_barrel_MC_separate[pr]);
        removeNegativeBins(h_mHits_endcap_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_barrel_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_endcap_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_barrel_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_endcap_MC_separate[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr]);

        for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinError(i_bin) /
                                                               (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinError(i_bin) /
                                                               (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        }

        Color_t color = kOrange - 5;
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
        h_TrkIso_barrel_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetFillColor(color);
        h_mHits_barrel_MC_separate[pr]->SetFillColor(color);
        h_mHits_endcap_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetFillColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetLineColor(color);
        h_mHits_barrel_MC_separate[pr]->SetLineColor(color);
        h_mHits_endcap_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetLineColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
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
        h_HoverE_barrel_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetDirectory(0);      
        h_PFiso_Rho_barrel_MC_separate[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_separate[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_separate[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetDirectory(0);
        h_mHits_barrel_MC_separate[pr]->SetDirectory(0);
        h_mHits_endcap_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_barrel_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_endcap_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_barrel_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_endcap_MC_separate[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno_density[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr]->SetDirectory(0);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_template[pr][ih]);
            file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            removeNegativeBins(h_PFiso_Rho_barrel_template[pr][ih]);
            removeNegativeBins(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetFillColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetFillColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetLineColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetLineColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetDirectory(0);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetDirectory(0);

            if (ih < nPtBin_ele)
            {
                file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_template[pr][ih]);
                file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                removeNegativeBins(h_PFiso_Rho_endcap_template[pr][ih]);
                removeNegativeBins(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetFillColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetFillColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetLineColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetLineColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetDirectory(0);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetDirectory(0);
            }
        }

        if (pr == _DY_10to50)
        {
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
            h_TrkIso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_TrkIso_barrel_MC_nume[pr]->Clone("h_TrkIso_barrel_nume_DY")));
            h_TrkIso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_TrkIso_endcap_MC_nume[pr]->Clone("h_TrkIso_endcap_nume_DY")));
            h_TrkIso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_TrkIso_barrel_MC_deno[pr]->Clone("h_TrkIso_barrel_deno_DY")));
            h_TrkIso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_TrkIso_endcap_MC_deno[pr]->Clone("h_TrkIso_endcap_deno_DY")));
            h_TrkIso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_TrkIso_barrel_MC_ctrl[pr]->Clone("h_TrkIso_barrel_ctrl_DY")));
            h_TrkIso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_TrkIso_endcap_MC_ctrl[pr]->Clone("h_TrkIso_endcap_ctrl_DY")));
            h_ECALiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_ECALiso_barrel_MC_nume[pr]->Clone("h_ECALiso_barrel_nume_DY")));
            h_ECALiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_ECALiso_endcap_MC_nume[pr]->Clone("h_ECALiso_endcap_nume_DY")));
            h_ECALiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_ECALiso_barrel_MC_deno[pr]->Clone("h_ECALiso_barrel_deno_DY")));
            h_ECALiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_ECALiso_endcap_MC_deno[pr]->Clone("h_ECALiso_endcap_deno_DY")));
            h_ECALiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_ECALiso_barrel_MC_ctrl[pr]->Clone("h_ECALiso_barrel_ctrl_DY")));
            h_ECALiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_ECALiso_endcap_MC_ctrl[pr]->Clone("h_ECALiso_endcap_ctrl_DY")));
            h_HCALiso_barrel_MC_nume[_DY_Full] = ((TH1D*)(h_HCALiso_barrel_MC_nume[pr]->Clone("h_HCALiso_barrel_nume_DY")));
            h_HCALiso_endcap_MC_nume[_DY_Full] = ((TH1D*)(h_HCALiso_endcap_MC_nume[pr]->Clone("h_HCALiso_endcap_nume_DY")));
            h_HCALiso_barrel_MC_deno[_DY_Full] = ((TH1D*)(h_HCALiso_barrel_MC_deno[pr]->Clone("h_HCALiso_barrel_deno_DY")));
            h_HCALiso_endcap_MC_deno[_DY_Full] = ((TH1D*)(h_HCALiso_endcap_MC_deno[pr]->Clone("h_HCALiso_endcap_deno_DY")));
            h_HCALiso_barrel_MC_ctrl[_DY_Full] = ((TH1D*)(h_HCALiso_barrel_MC_ctrl[pr]->Clone("h_HCALiso_barrel_ctrl_DY")));
            h_HCALiso_endcap_MC_ctrl[_DY_Full] = ((TH1D*)(h_HCALiso_endcap_MC_ctrl[pr]->Clone("h_HCALiso_endcap_ctrl_DY")));
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
            h_HoverE_barrel_template_int_MC[_DY_Full] = ((TH1D*)(h_HoverE_barrel_template_int_MC[pr]->Clone("h_HoverE_barrel_template_int_DY")));
            h_HoverE_barrel_jetTemplate_int_MC[_DY_Full] = ((TH1D*)(h_HoverE_barrel_jetTemplate_int_MC[pr]->Clone("h_HoverE_barrel_jetTemplate_int_DY")));
            h_HoverE_endcap_template_int_MC[_DY_Full] = ((TH1D*)(h_HoverE_endcap_template_int_MC[pr]->Clone("h_HoverE_endcap_template_int_DY")));
            h_HoverE_endcap_jetTemplate_int_MC[_DY_Full] = ((TH1D*)(h_HoverE_endcap_jetTemplate_int_MC[pr]->Clone("h_HoverE_endcap_jetTemplate_int_DY")));
            h_PFiso_Rho_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_separate[pr]->Clone("h_PFiso_Rho_barrel_separate_DY")));
            h_PFiso_Rho_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_separate[pr]->Clone("h_PFiso_Rho_endcap_separate_DY")));
            h_SigmaIEtaIEta_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_separate[pr]->Clone("h_SigmaIEtaIEta_barrel_separate_DY")));
            h_SigmaIEtaIEta_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_separate[pr]->Clone("h_SigmaIEtaIEta_endcap_separate_DY")));
            h_dEtaInSeed_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_separate[pr]->Clone("h_dEtaInSeed_barrel_separate_DY")));
            h_dEtaInSeed_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_separate[pr]->Clone("h_dEtaInSeed_endcap_separate_DY")));
            h_dPhiIn_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_separate[pr]->Clone("h_dPhiIn_barrel_separate_DY")));
            h_dPhiIn_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_separate[pr]->Clone("h_dPhiIn_endcap_separate_DY")));
            h_HoverE_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_HoverE_barrel_MC_separate[pr]->Clone("h_HoverE_barrel_separate_DY")));
            h_HoverE_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_HoverE_endcap_MC_separate[pr]->Clone("h_HoverE_endcap_separate_DY")));
            h_InvEminusInvP_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_separate[pr]->Clone("h_InvEminusInvP_barrel_separate_DY")));
            h_InvEminusInvP_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_separate[pr]->Clone("h_InvEminusInvP_endcap_separate_DY")));
            h_mHits_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_mHits_barrel_MC_separate[pr]->Clone("h_mHits_barrel_separate_DY")));
            h_mHits_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_mHits_endcap_MC_separate[pr]->Clone("h_mHits_endcap_separate_DY")));
            h_passConvVeto_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_passConvVeto_barrel_MC_separate[pr]->Clone("h_passConvVeto_barrel_separate_DY")));
            h_passConvVeto_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_passConvVeto_endcap_MC_separate[pr]->Clone("h_passConvVeto_endcap_separate_DY")));
            h_passMediumID_barrel_MC_separate[_DY_Full] = ((TH1D*)(h_passMediumID_barrel_MC_separate[pr]->Clone("h_passMediumID_barrel_separate_DY")));
            h_passMediumID_endcap_MC_separate[_DY_Full] = ((TH1D*)(h_passMediumID_endcap_MC_separate[pr]->Clone("h_passMediumID_endcap_separate_DY")));
            h_pT_barrel_MC_deno_density[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr]->Clone("h_pT_barrel_deno_density_DY_Full")));
            h_pT_endcap_MC_deno_density[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr]->Clone("h_pT_endcap_deno_density_DY_Full")));

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
            h_TrkIso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_nume[_DY_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_nume[_DY_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_deno[_DY_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_deno[_DY_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_ctrl[_DY_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_ctrl[_DY_Full]->SetDirectory(0);
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
            h_HoverE_barrel_template_int_MC[_DY_Full]->SetDirectory(0);
            h_HoverE_barrel_jetTemplate_int_MC[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_template_int_MC[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_jetTemplate_int_MC[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_mHits_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_mHits_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_passConvVeto_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_passConvVeto_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_passMediumID_barrel_MC_separate[_DY_Full]->SetDirectory(0);
            h_passMediumID_endcap_MC_separate[_DY_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno_density[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_DY_Full]->SetDirectory(0);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_DY_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_template[pr][ih]->Clone("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10)+"_DY")));
                h_PFiso_Rho_barrel_jetTemplate[_DY_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10)+"_DY")));
                h_PFiso_Rho_barrel_template[_DY_Full][ih]->SetDirectory(0);
                h_PFiso_Rho_barrel_jetTemplate[_DY_Full][ih]->SetDirectory(0);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_DY_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_template[pr][ih]->Clone("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10)+"_DY")));
                    h_PFiso_Rho_endcap_jetTemplate[_DY_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10)+"_DY")));
                    h_PFiso_Rho_endcap_template[_DY_Full][ih]->SetDirectory(0);
                    h_PFiso_Rho_endcap_jetTemplate[_DY_Full][ih]->SetDirectory(0);
                }
            }
        }
        else
        {
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
            h_TrkIso_barrel_MC_nume[_DY_Full]->Add(h_TrkIso_barrel_MC_nume[pr]);
            h_TrkIso_endcap_MC_nume[_DY_Full]->Add(h_TrkIso_endcap_MC_nume[pr]);
            h_TrkIso_barrel_MC_deno[_DY_Full]->Add(h_TrkIso_barrel_MC_deno[pr]);
            h_TrkIso_endcap_MC_deno[_DY_Full]->Add(h_TrkIso_endcap_MC_deno[pr]);
            h_TrkIso_barrel_MC_ctrl[_DY_Full]->Add(h_TrkIso_barrel_MC_ctrl[pr]);
            h_TrkIso_endcap_MC_ctrl[_DY_Full]->Add(h_TrkIso_endcap_MC_ctrl[pr]);
            h_ECALiso_barrel_MC_nume[_DY_Full]->Add(h_ECALiso_barrel_MC_nume[pr]);
            h_ECALiso_endcap_MC_nume[_DY_Full]->Add(h_ECALiso_endcap_MC_nume[pr]);
            h_ECALiso_barrel_MC_deno[_DY_Full]->Add(h_ECALiso_barrel_MC_deno[pr]);
            h_ECALiso_endcap_MC_deno[_DY_Full]->Add(h_ECALiso_endcap_MC_deno[pr]);
            h_ECALiso_barrel_MC_ctrl[_DY_Full]->Add(h_ECALiso_barrel_MC_ctrl[pr]);
            h_ECALiso_endcap_MC_ctrl[_DY_Full]->Add(h_ECALiso_endcap_MC_ctrl[pr]);
            h_HCALiso_barrel_MC_nume[_DY_Full]->Add(h_HCALiso_barrel_MC_nume[pr]);
            h_HCALiso_endcap_MC_nume[_DY_Full]->Add(h_HCALiso_endcap_MC_nume[pr]);
            h_HCALiso_barrel_MC_deno[_DY_Full]->Add(h_HCALiso_barrel_MC_deno[pr]);
            h_HCALiso_endcap_MC_deno[_DY_Full]->Add(h_HCALiso_endcap_MC_deno[pr]);
            h_HCALiso_barrel_MC_ctrl[_DY_Full]->Add(h_HCALiso_barrel_MC_ctrl[pr]);
            h_HCALiso_endcap_MC_ctrl[_DY_Full]->Add(h_HCALiso_endcap_MC_ctrl[pr]);
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
            h_HoverE_barrel_template_int_MC[_DY_Full]->Add(h_HoverE_barrel_template_int_MC[pr]);
            h_HoverE_barrel_jetTemplate_int_MC[_DY_Full]->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
            h_HoverE_endcap_template_int_MC[_DY_Full]->Add(h_HoverE_endcap_template_int_MC[pr]);
            h_HoverE_endcap_jetTemplate_int_MC[_DY_Full]->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
            h_PFiso_Rho_barrel_MC_separate[_DY_Full]->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
            h_PFiso_Rho_endcap_MC_separate[_DY_Full]->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
            h_SigmaIEtaIEta_barrel_MC_separate[_DY_Full]->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
            h_SigmaIEtaIEta_endcap_MC_separate[_DY_Full]->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
            h_dEtaInSeed_barrel_MC_separate[_DY_Full]->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
            h_dEtaInSeed_endcap_MC_separate[_DY_Full]->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
            h_dPhiIn_barrel_MC_separate[_DY_Full]->Add(h_dPhiIn_barrel_MC_separate[pr]);
            h_dPhiIn_endcap_MC_separate[_DY_Full]->Add(h_dPhiIn_endcap_MC_separate[pr]);
            h_HoverE_barrel_MC_separate[_DY_Full]->Add(h_HoverE_barrel_MC_separate[pr]);
            h_HoverE_endcap_MC_separate[_DY_Full]->Add(h_HoverE_endcap_MC_separate[pr]);
            h_InvEminusInvP_barrel_MC_separate[_DY_Full]->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
            h_InvEminusInvP_endcap_MC_separate[_DY_Full]->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
            h_mHits_barrel_MC_separate[_DY_Full]->Add(h_mHits_barrel_MC_separate[pr]);
            h_mHits_endcap_MC_separate[_DY_Full]->Add(h_mHits_endcap_MC_separate[pr]);
            h_passConvVeto_barrel_MC_separate[_DY_Full]->Add(h_passConvVeto_barrel_MC_separate[pr]);
            h_passConvVeto_endcap_MC_separate[_DY_Full]->Add(h_passConvVeto_endcap_MC_separate[pr]);
            h_passMediumID_barrel_MC_separate[_DY_Full]->Add(h_passMediumID_barrel_MC_separate[pr]);
            h_passMediumID_endcap_MC_separate[_DY_Full]->Add(h_passMediumID_endcap_MC_separate[pr]);
            h_pT_barrel_MC_deno_density[_DY_Full]->Add(h_pT_barrel_MC_deno_density[pr]);
            h_pT_endcap_MC_deno_density[_DY_Full]->Add(h_pT_endcap_MC_deno_density[pr]);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_DY_Full][ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
                h_PFiso_Rho_barrel_jetTemplate[_DY_Full][ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_DY_Full][ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
                    h_PFiso_Rho_endcap_jetTemplate[_DY_Full][ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
                }
            }
        }

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
        s_TrkIso_barrel_nume->Add(h_TrkIso_barrel_MC_nume[pr]);
        s_TrkIso_endcap_nume->Add(h_TrkIso_endcap_MC_nume[pr]);
        s_TrkIso_barrel_deno->Add(h_TrkIso_barrel_MC_deno[pr]);
        s_TrkIso_endcap_deno->Add(h_TrkIso_endcap_MC_deno[pr]);
        s_TrkIso_barrel_ctrl->Add(h_TrkIso_barrel_MC_ctrl[pr]);
        s_TrkIso_endcap_ctrl->Add(h_TrkIso_endcap_MC_ctrl[pr]);
        s_ECALiso_barrel_nume->Add(h_ECALiso_barrel_MC_nume[pr]);
        s_ECALiso_endcap_nume->Add(h_ECALiso_endcap_MC_nume[pr]);
        s_ECALiso_barrel_deno->Add(h_ECALiso_barrel_MC_deno[pr]);
        s_ECALiso_endcap_deno->Add(h_ECALiso_endcap_MC_deno[pr]);
        s_ECALiso_barrel_ctrl->Add(h_ECALiso_barrel_MC_ctrl[pr]);
        s_ECALiso_endcap_ctrl->Add(h_ECALiso_endcap_MC_ctrl[pr]);
        s_HCALiso_barrel_nume->Add(h_HCALiso_barrel_MC_nume[pr]);
        s_HCALiso_endcap_nume->Add(h_HCALiso_endcap_MC_nume[pr]);
        s_HCALiso_barrel_deno->Add(h_HCALiso_barrel_MC_deno[pr]);
        s_HCALiso_endcap_deno->Add(h_HCALiso_endcap_MC_deno[pr]);
        s_HCALiso_barrel_ctrl->Add(h_HCALiso_barrel_MC_ctrl[pr]);
        s_HCALiso_endcap_ctrl->Add(h_HCALiso_endcap_MC_ctrl[pr]);
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
        s_HoverE_barrel_template_int->Add(h_HoverE_barrel_template_int_MC[pr]);
        s_HoverE_barrel_jetTemplate_int->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
        s_HoverE_endcap_template_int->Add(h_HoverE_endcap_template_int_MC[pr]);
        s_HoverE_endcap_jetTemplate_int->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
        s_PFiso_Rho_barrel_separate->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
        s_PFiso_Rho_endcap_separate->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
        s_SigmaIEtaIEta_barrel_separate->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        s_SigmaIEtaIEta_endcap_separate->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        s_dEtaInSeed_barrel_separate->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
        s_dEtaInSeed_endcap_separate->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
        s_dPhiIn_barrel_separate->Add(h_dPhiIn_barrel_MC_separate[pr]);
        s_dPhiIn_endcap_separate->Add(h_dPhiIn_endcap_MC_separate[pr]);
        s_HoverE_barrel_separate->Add(h_HoverE_barrel_MC_separate[pr]);
        s_HoverE_endcap_separate->Add(h_HoverE_endcap_MC_separate[pr]);
        s_InvEminusInvP_barrel_separate->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
        s_InvEminusInvP_endcap_separate->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
        s_mHits_barrel_separate->Add(h_mHits_barrel_MC_separate[pr]);
        s_mHits_endcap_separate->Add(h_mHits_endcap_MC_separate[pr]);
        s_passConvVeto_barrel_separate->Add(h_passConvVeto_barrel_MC_separate[pr]);
        s_passConvVeto_endcap_separate->Add(h_passConvVeto_endcap_MC_separate[pr]);
        s_passMediumID_barrel_separate->Add(h_passMediumID_barrel_MC_separate[pr]);
        s_passMediumID_endcap_separate->Add(h_passMediumID_endcap_MC_separate[pr]);
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr]);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            s_PFiso_Rho_barrel_template[ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
            s_PFiso_Rho_barrel_jetTemplate[ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);
            if (ih < nPtBin_ele)
            {
                s_PFiso_Rho_endcap_template[ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
                s_PFiso_Rho_endcap_jetTemplate[ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
            }
        }

        file->Close();
    }

    cout << "DY received" << endl;

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

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
        file->GetObject("h_TrkIso_barrel_nume", h_TrkIso_barrel_MC_nume[pr]);
        file->GetObject("h_TrkIso_endcap_nume", h_TrkIso_endcap_MC_nume[pr]);
        file->GetObject("h_TrkIso_barrel_deno", h_TrkIso_barrel_MC_deno[pr]);
        file->GetObject("h_TrkIso_endcap_deno", h_TrkIso_endcap_MC_deno[pr]);
        file->GetObject("h_TrkIso_barrel_ctrl", h_TrkIso_barrel_MC_ctrl[pr]);
        file->GetObject("h_TrkIso_endcap_ctrl", h_TrkIso_endcap_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_barrel_nume", h_ECALiso_barrel_MC_nume[pr]);
        file->GetObject("h_ECALiso_endcap_nume", h_ECALiso_endcap_MC_nume[pr]);
        file->GetObject("h_ECALiso_barrel_deno", h_ECALiso_barrel_MC_deno[pr]);
        file->GetObject("h_ECALiso_endcap_deno", h_ECALiso_endcap_MC_deno[pr]);
        file->GetObject("h_ECALiso_barrel_ctrl", h_ECALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_endcap_ctrl", h_ECALiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_barrel_nume", h_HCALiso_barrel_MC_nume[pr]);
        file->GetObject("h_HCALiso_endcap_nume", h_HCALiso_endcap_MC_nume[pr]);
        file->GetObject("h_HCALiso_barrel_deno", h_HCALiso_barrel_MC_deno[pr]);
        file->GetObject("h_HCALiso_endcap_deno", h_HCALiso_endcap_MC_deno[pr]);
        file->GetObject("h_HCALiso_barrel_ctrl", h_HCALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_endcap_ctrl", h_HCALiso_endcap_MC_ctrl[pr]);
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
        file->GetObject("h_HoverE_barrel_template_int", h_HoverE_barrel_template_int_MC[pr]);
        file->GetObject("h_HoverE_barrel_jetTemplate_int", h_HoverE_barrel_jetTemplate_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_template_int", h_HoverE_endcap_template_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_jetTemplate_int", h_HoverE_endcap_jetTemplate_int_MC[pr]);
        file->GetObject("h_PFiso_Rho_barrel_separate", h_PFiso_Rho_barrel_MC_separate[pr]);
        file->GetObject("h_PFiso_Rho_endcap_separate", h_PFiso_Rho_endcap_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_barrel_separate", h_dEtaInSeed_barrel_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_endcap_separate", h_dEtaInSeed_endcap_MC_separate[pr]);
        file->GetObject("h_dPhiIn_barrel_separate", h_dPhiIn_barrel_MC_separate[pr]);
        file->GetObject("h_dPhiIn_endcap_separate", h_dPhiIn_endcap_MC_separate[pr]);
        file->GetObject("h_HoverE_barrel_separate", h_HoverE_barrel_MC_separate[pr]);
        file->GetObject("h_HoverE_endcap_separate", h_HoverE_endcap_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_barrel_separate", h_InvEminusInvP_barrel_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_endcap_separate", h_InvEminusInvP_endcap_MC_separate[pr]);
        file->GetObject("h_mHits_barrel_separate", h_mHits_barrel_MC_separate[pr]);
        file->GetObject("h_mHits_endcap_separate", h_mHits_endcap_MC_separate[pr]);
        file->GetObject("h_passConvVeto_barrel_separate", h_passConvVeto_barrel_MC_separate[pr]);
        file->GetObject("h_passConvVeto_endcap_separate", h_passConvVeto_endcap_MC_separate[pr]);
        file->GetObject("h_passMediumID_barrel_separate", h_passMediumID_barrel_MC_separate[pr]);
        file->GetObject("h_passMediumID_endcap_separate", h_passMediumID_endcap_MC_separate[pr]);
        h_pT_barrel_MC_deno_density[pr] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr])));
        h_pT_endcap_MC_deno_density[pr] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr])));

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
        removeNegativeBins(h_TrkIso_barrel_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_ctrl[pr]);
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
        removeNegativeBins(h_HoverE_barrel_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_barrel_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_separate[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_separate[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_separate[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_separate[pr]);
        removeNegativeBins(h_mHits_barrel_MC_separate[pr]);
        removeNegativeBins(h_mHits_endcap_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_barrel_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_endcap_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_barrel_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_endcap_MC_separate[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr]);

        for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        }

        Color_t color = kRed + 3;
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
        h_TrkIso_barrel_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetFillColor(color);
        h_mHits_barrel_MC_separate[pr]->SetFillColor(color);
        h_mHits_endcap_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetFillColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetLineColor(color);
        h_mHits_barrel_MC_separate[pr]->SetLineColor(color);
        h_mHits_endcap_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetLineColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
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
        h_HoverE_barrel_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_separate[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_separate[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetDirectory(0);
        h_mHits_barrel_MC_separate[pr]->SetDirectory(0);
        h_mHits_endcap_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_barrel_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_endcap_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_barrel_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_endcap_MC_separate[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno_density[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr]->SetDirectory(0);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_template[pr][ih]);
            file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            removeNegativeBins(h_PFiso_Rho_barrel_template[pr][ih]);
            removeNegativeBins(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetFillColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetFillColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetLineColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetLineColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetDirectory(0);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetDirectory(0);

            if (ih < nPtBin_ele)
            {
                file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_template[pr][ih]);
                file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                removeNegativeBins(h_PFiso_Rho_endcap_template[pr][ih]);
                removeNegativeBins(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetFillColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetFillColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetLineColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetLineColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetDirectory(0);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetDirectory(0);

            }
        }

        if (pr == _QCDEMEnriched_20to30)
        {
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
            h_TrkIso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_barrel_MC_nume[pr]->Clone("h_TrkIso_barrel_nume_QCD")));
            h_TrkIso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_endcap_MC_nume[pr]->Clone("h_TrkIso_endcap_nume_QCD")));
            h_TrkIso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_barrel_MC_deno[pr]->Clone("h_TrkIso_barrel_deno_QCD")));
            h_TrkIso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_endcap_MC_deno[pr]->Clone("h_TrkIso_endcap_deno_QCD")));
            h_TrkIso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_barrel_MC_ctrl[pr]->Clone("h_TrkIso_barrel_ctrl_QCD")));
            h_TrkIso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_endcap_MC_ctrl[pr]->Clone("h_TrkIso_endcap_ctrl_QCD")));
            h_ECALiso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_barrel_MC_nume[pr]->Clone("h_ECALiso_barrel_nume_QCD")));
            h_ECALiso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_endcap_MC_nume[pr]->Clone("h_ECALiso_endcap_nume_QCD")));
            h_ECALiso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_barrel_MC_deno[pr]->Clone("h_ECALiso_barrel_deno_QCD")));
            h_ECALiso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_endcap_MC_deno[pr]->Clone("h_ECALiso_endcap_deno_QCD")));
            h_ECALiso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_barrel_MC_ctrl[pr]->Clone("h_ECALiso_barrel_ctrl_QCD")));
            h_ECALiso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_endcap_MC_ctrl[pr]->Clone("h_ECALiso_endcap_ctrl_QCD")));
            h_HCALiso_barrel_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_barrel_MC_nume[pr]->Clone("h_HCALiso_barrel_nume_QCD")));
            h_HCALiso_endcap_MC_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_endcap_MC_nume[pr]->Clone("h_HCALiso_endcap_nume_QCD")));
            h_HCALiso_barrel_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_barrel_MC_deno[pr]->Clone("h_HCALiso_barrel_deno_QCD")));
            h_HCALiso_endcap_MC_deno[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_endcap_MC_deno[pr]->Clone("h_HCALiso_endcap_deno_QCD")));
            h_HCALiso_barrel_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_barrel_MC_ctrl[pr]->Clone("h_HCALiso_barrel_ctrl_QCD")));
            h_HCALiso_endcap_MC_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_endcap_MC_ctrl[pr]->Clone("h_HCALiso_endcap_ctrl_QCD")));
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
            h_HoverE_barrel_template_int_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_template_int_MC[pr]->Clone("h_HoverE_barrel_template_int_QCD")));
            h_HoverE_barrel_jetTemplate_int_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_jetTemplate_int_MC[pr]->Clone("h_HoverE_barrel_jetTemplate_int_QCD")));
            h_HoverE_endcap_template_int_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_template_int_MC[pr]->Clone("h_HoverE_endcap_template_int_QCD")));
            h_HoverE_endcap_jetTemplate_int_MC[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_jetTemplate_int_MC[pr]->Clone("h_HoverE_endcap_jetTemplate_int_QCD")));            
            h_PFiso_Rho_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_separate[pr]->Clone("h_PFiso_Rho_barrel_separate_QCD")));
            h_PFiso_Rho_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_separate[pr]->Clone("h_PFiso_Rho_endcap_separate_QCD")));
            h_SigmaIEtaIEta_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_separate[pr]->Clone("h_SigmaIEtaIEta_barrel_separate_QCD")));
            h_SigmaIEtaIEta_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_separate[pr]->Clone("h_SigmaIEtaIEta_endcap_separate_QCD")));
            h_dEtaInSeed_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_separate[pr]->Clone("h_dEtaInSeed_barrel_separate_QCD")));
            h_dEtaInSeed_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_separate[pr]->Clone("h_dEtaInSeed_endcap_separate_QCD")));
            h_dPhiIn_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_separate[pr]->Clone("h_dPhiIn_barrel_separate_QCD")));
            h_dPhiIn_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_separate[pr]->Clone("h_dPhiIn_endcap_separate_QCD")));
            h_HoverE_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_barrel_MC_separate[pr]->Clone("h_HoverE_barrel_separate_QCD")));
            h_HoverE_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_HoverE_endcap_MC_separate[pr]->Clone("h_HoverE_endcap_separate_QCD")));
            h_InvEminusInvP_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_separate[pr]->Clone("h_InvEminusInvP_barrel_separate_QCD")));
            h_InvEminusInvP_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_separate[pr]->Clone("h_InvEminusInvP_endcap_separate_QCD")));
            h_mHits_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_mHits_barrel_MC_separate[pr]->Clone("h_mHits_barrel_separate_QCD")));
            h_mHits_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_mHits_endcap_MC_separate[pr]->Clone("h_mHits_endcap_separate_QCD")));
            h_passConvVeto_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_passConvVeto_barrel_MC_separate[pr]->Clone("h_passConvVeto_barrel_separate_QCD")));
            h_passConvVeto_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_passConvVeto_endcap_MC_separate[pr]->Clone("h_passConvVeto_endcap_separate_QCD")));
            h_passMediumID_barrel_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_passMediumID_barrel_MC_separate[pr]->Clone("h_passMediumID_barrel_separate_QCD")));
            h_passMediumID_endcap_MC_separate[_QCDEMEnriched_Full] = ((TH1D*)(h_passMediumID_endcap_MC_separate[pr]->Clone("h_passMediumID_endcap_separate_QCD")));
            h_pT_barrel_MC_deno_density[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr]->Clone("h_pT_barrel_deno_density_QCD")));
            h_pT_endcap_MC_deno_density[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr]->Clone("h_pT_endcap_deno_density_QCD")));

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
            h_TrkIso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_deno[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
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
            h_HoverE_barrel_template_int_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_barrel_jetTemplate_int_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_template_int_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_jetTemplate_int_MC[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mHits_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mHits_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_passConvVeto_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_passConvVeto_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_passMediumID_barrel_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_passMediumID_endcap_MC_separate[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno_density[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_QCDEMEnriched_Full]->SetDirectory(0);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_QCDEMEnriched_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_template[pr][ih]->Clone("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10)+"_QCD")));
                h_PFiso_Rho_barrel_jetTemplate[_QCDEMEnriched_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10)+"_QCD")));
                h_PFiso_Rho_barrel_template[_QCDEMEnriched_Full][ih]->SetDirectory(0);
                h_PFiso_Rho_barrel_jetTemplate[_QCDEMEnriched_Full][ih]->SetDirectory(0);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_QCDEMEnriched_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_template[pr][ih]->Clone("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10)+"_QCD")));
                    h_PFiso_Rho_endcap_jetTemplate[_QCDEMEnriched_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10)+"_QCD")));
                    h_PFiso_Rho_endcap_template[_QCDEMEnriched_Full][ih]->SetDirectory(0);
                    h_PFiso_Rho_endcap_jetTemplate[_QCDEMEnriched_Full][ih]->SetDirectory(0);
                }
            }
        }
        else
        {
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
            h_TrkIso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_TrkIso_barrel_MC_nume[pr]);
            h_TrkIso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_TrkIso_endcap_MC_nume[pr]);
            h_TrkIso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_TrkIso_barrel_MC_deno[pr]);
            h_TrkIso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_TrkIso_endcap_MC_deno[pr]);
            h_TrkIso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_TrkIso_barrel_MC_ctrl[pr]);
            h_TrkIso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_TrkIso_endcap_MC_ctrl[pr]);
            h_ECALiso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_ECALiso_barrel_MC_nume[pr]);
            h_ECALiso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_ECALiso_endcap_MC_nume[pr]);
            h_ECALiso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_ECALiso_barrel_MC_deno[pr]);
            h_ECALiso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_ECALiso_endcap_MC_deno[pr]);
            h_ECALiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_ECALiso_barrel_MC_ctrl[pr]);
            h_ECALiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_ECALiso_endcap_MC_ctrl[pr]);
            h_HCALiso_barrel_MC_nume[_QCDEMEnriched_Full]->Add(h_HCALiso_barrel_MC_nume[pr]);
            h_HCALiso_endcap_MC_nume[_QCDEMEnriched_Full]->Add(h_HCALiso_endcap_MC_nume[pr]);
            h_HCALiso_barrel_MC_deno[_QCDEMEnriched_Full]->Add(h_HCALiso_barrel_MC_deno[pr]);
            h_HCALiso_endcap_MC_deno[_QCDEMEnriched_Full]->Add(h_HCALiso_endcap_MC_deno[pr]);
            h_HCALiso_barrel_MC_ctrl[_QCDEMEnriched_Full]->Add(h_HCALiso_barrel_MC_ctrl[pr]);
            h_HCALiso_endcap_MC_ctrl[_QCDEMEnriched_Full]->Add(h_HCALiso_endcap_MC_ctrl[pr]);
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
            h_HoverE_barrel_template_int_MC[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_template_int_MC[pr]);
            h_HoverE_barrel_jetTemplate_int_MC[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
            h_HoverE_endcap_template_int_MC[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_template_int_MC[pr]);
            h_HoverE_endcap_jetTemplate_int_MC[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
            h_PFiso_Rho_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
            h_PFiso_Rho_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
            h_SigmaIEtaIEta_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
            h_SigmaIEtaIEta_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
            h_dEtaInSeed_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
            h_dEtaInSeed_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
            h_dPhiIn_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_dPhiIn_barrel_MC_separate[pr]);
            h_dPhiIn_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_dPhiIn_endcap_MC_separate[pr]);
            h_HoverE_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_HoverE_barrel_MC_separate[pr]);
            h_HoverE_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_HoverE_endcap_MC_separate[pr]);
            h_InvEminusInvP_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
            h_InvEminusInvP_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
            h_mHits_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_mHits_barrel_MC_separate[pr]);
            h_mHits_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_mHits_endcap_MC_separate[pr]);
            h_passConvVeto_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_passConvVeto_barrel_MC_separate[pr]);
            h_passConvVeto_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_passConvVeto_endcap_MC_separate[pr]);
            h_passMediumID_barrel_MC_separate[_QCDEMEnriched_Full]->Add(h_passMediumID_barrel_MC_separate[pr]);
            h_passMediumID_endcap_MC_separate[_QCDEMEnriched_Full]->Add(h_passMediumID_endcap_MC_separate[pr]);
            h_pT_barrel_MC_deno_density[_QCDEMEnriched_Full]->Add(h_pT_barrel_MC_deno_density[pr]);
            h_pT_endcap_MC_deno_density[_QCDEMEnriched_Full]->Add(h_pT_endcap_MC_deno_density[pr]);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
                h_PFiso_Rho_barrel_jetTemplate[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
                    h_PFiso_Rho_endcap_jetTemplate[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
                }
            }
        }

//        s_PFiso_Rho_barrel_nume->Add(h_PFiso_Rho_barrel_MC_nume[pr]);
//        s_PFiso_Rho_endcap_nume->Add(h_PFiso_Rho_endcap_MC_nume[pr]);
//        s_PFiso_Rho_barrel_deno->Add(h_PFiso_Rho_barrel_MC_deno[pr]);
//        s_PFiso_Rho_endcap_deno->Add(h_PFiso_Rho_endcap_MC_deno[pr]);
//        s_PFiso_Rho_barrel_ctrl->Add(h_PFiso_Rho_barrel_MC_ctrl[pr]);
//        s_PFiso_Rho_endcap_ctrl->Add(h_PFiso_Rho_endcap_MC_ctrl[pr]);
//        s_pT_barrel_nume->Add(h_pT_barrel_MC_nume[pr]);
//        s_pT_endcap_nume->Add(h_pT_endcap_MC_nume[pr]);
//        s_pT_barrel_deno->Add(h_pT_barrel_MC_deno[pr]);
//        s_pT_endcap_deno->Add(h_pT_endcap_MC_deno[pr]);
//        s_pT_barrel_ctrl->Add(h_pT_barrel_MC_ctrl[pr]);
//        s_pT_endcap_ctrl->Add(h_pT_endcap_MC_ctrl[pr]);
//        s_SigmaIEtaIEta_barrel_nume->Add(h_SigmaIEtaIEta_barrel_MC_nume[pr]);
//        s_SigmaIEtaIEta_endcap_nume->Add(h_SigmaIEtaIEta_endcap_MC_nume[pr]);
//        s_SigmaIEtaIEta_barrel_deno->Add(h_SigmaIEtaIEta_barrel_MC_deno[pr]);
//        s_SigmaIEtaIEta_endcap_deno->Add(h_SigmaIEtaIEta_endcap_MC_deno[pr]);
//        s_SigmaIEtaIEta_barrel_ctrl->Add(h_SigmaIEtaIEta_barrel_MC_ctrl[pr]);
//        s_SigmaIEtaIEta_endcap_ctrl->Add(h_SigmaIEtaIEta_endcap_MC_ctrl[pr]);
//        s_dEtaInSeed_barrel_nume->Add(h_dEtaInSeed_barrel_MC_nume[pr]);
//        s_dEtaInSeed_endcap_nume->Add(h_dEtaInSeed_endcap_MC_nume[pr]);
//        s_dEtaInSeed_barrel_deno->Add(h_dEtaInSeed_barrel_MC_deno[pr]);
//        s_dEtaInSeed_endcap_deno->Add(h_dEtaInSeed_endcap_MC_deno[pr]);
//        s_dEtaInSeed_barrel_ctrl->Add(h_dEtaInSeed_barrel_MC_ctrl[pr]);
//        s_dEtaInSeed_endcap_ctrl->Add(h_dEtaInSeed_endcap_MC_ctrl[pr]);
//        s_dPhiIn_barrel_nume->Add(h_dPhiIn_barrel_MC_nume[pr]);
//        s_dPhiIn_endcap_nume->Add(h_dPhiIn_endcap_MC_nume[pr]);
//        s_dPhiIn_barrel_deno->Add(h_dPhiIn_barrel_MC_deno[pr]);
//        s_dPhiIn_endcap_deno->Add(h_dPhiIn_endcap_MC_deno[pr]);
//        s_dPhiIn_barrel_ctrl->Add(h_dPhiIn_barrel_MC_ctrl[pr]);
//        s_dPhiIn_endcap_ctrl->Add(h_dPhiIn_endcap_MC_ctrl[pr]);
//        s_HoverE_barrel_nume->Add(h_HoverE_barrel_MC_nume[pr]);
//        s_HoverE_endcap_nume->Add(h_HoverE_endcap_MC_nume[pr]);
//        s_HoverE_barrel_deno->Add(h_HoverE_barrel_MC_deno[pr]);
//        s_HoverE_endcap_deno->Add(h_HoverE_endcap_MC_deno[pr]);
//        s_HoverE_barrel_ctrl->Add(h_HoverE_barrel_MC_ctrl[pr]);
//        s_HoverE_endcap_ctrl->Add(h_HoverE_endcap_MC_ctrl[pr]);
//        s_InvEminusInvP_barrel_nume->Add(h_InvEminusInvP_barrel_MC_nume[pr]);
//        s_InvEminusInvP_endcap_nume->Add(h_InvEminusInvP_endcap_MC_nume[pr]);
//        s_InvEminusInvP_barrel_deno->Add(h_InvEminusInvP_barrel_MC_deno[pr]);
//        s_InvEminusInvP_endcap_deno->Add(h_InvEminusInvP_endcap_MC_deno[pr]);
//        s_InvEminusInvP_barrel_ctrl->Add(h_InvEminusInvP_barrel_MC_ctrl[pr]);
//        s_InvEminusInvP_endcap_ctrl->Add(h_InvEminusInvP_endcap_MC_ctrl[pr]);
//        s_TrkIso_barrel_nume->Add(h_TrkIso_barrel_MC_nume[pr]);
//        s_TrkIso_endcap_nume->Add(h_TrkIso_endcap_MC_nume[pr]);
//        s_TrkIso_barrel_deno->Add(h_TrkIso_barrel_MC_deno[pr]);
//        s_TrkIso_endcap_deno->Add(h_TrkIso_endcap_MC_deno[pr]);
//        s_TrkIso_barrel_ctrl->Add(h_TrkIso_barrel_MC_ctrl[pr]);
//        s_TrkIso_endcap_ctrl->Add(h_TrkIso_endcap_MC_ctrl[pr]);
//        s_ECALiso_barrel_nume->Add(h_ECALiso_barrel_MC_nume[pr]);
//        s_ECALiso_endcap_nume->Add(h_ECALiso_endcap_MC_nume[pr]);
//        s_ECALiso_barrel_deno->Add(h_ECALiso_barrel_MC_deno[pr]);
//        s_ECALiso_endcap_deno->Add(h_ECALiso_endcap_MC_deno[pr]);
//        s_ECALiso_barrel_ctrl->Add(h_ECALiso_barrel_MC_ctrl[pr]);
//        s_ECALiso_endcap_ctrl->Add(h_ECALiso_endcap_MC_ctrl[pr]);
//        s_HCALiso_barrel_nume->Add(h_HCALiso_barrel_MC_nume[pr]);
//        s_HCALiso_endcap_nume->Add(h_HCALiso_endcap_MC_nume[pr]);
//        s_HCALiso_barrel_deno->Add(h_HCALiso_barrel_MC_deno[pr]);
//        s_HCALiso_endcap_deno->Add(h_HCALiso_endcap_MC_deno[pr]);
//        s_HCALiso_barrel_ctrl->Add(h_HCALiso_barrel_MC_ctrl[pr]);
//        s_HCALiso_endcap_ctrl->Add(h_HCALiso_endcap_MC_ctrl[pr]);
//        s_MET->Add(h_MET_MC[pr]);
//        s_MT_barrel_nume->Add(h_MT_barrel_MC_nume[pr]);
//        s_MT_endcap_nume->Add(h_MT_endcap_MC_nume[pr]);
//        s_MT_barrel_deno->Add(h_MT_barrel_MC_deno[pr]);
//        s_MT_endcap_deno->Add(h_MT_endcap_MC_deno[pr]);
//        s_MT_barrel_ctrl->Add(h_MT_barrel_MC_ctrl[pr]);
//        s_MT_endcap_ctrl->Add(h_MT_endcap_MC_ctrl[pr]);
//        s_eta->Add(h_eta_MC[pr]);
//        s_nVTX->Add(h_nVTX_MC[pr]);
//        s_mass_test->Add(h_mass_test_MC[pr]);
//        s_HoverE_barrel_template_int->Add(h_HoverE_barrel_template_int_MC[pr]);
//        s_HoverE_barrel_jetTemplate_int->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
//        s_HoverE_endcap_template_int->Add(h_HoverE_endcap_template_int_MC[pr]);
//        s_HoverE_endcap_jetTemplate_int->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
//        s_PFiso_Rho_barrel_separate->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
//        s_PFiso_Rho_endcap_separate->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
//        s_SigmaIEtaIEta_barrel_separate->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
//        s_SigmaIEtaIEta_endcap_separate->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
//        s_dEtaInSeed_barrel_separate->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
//        s_dEtaInSeed_endcap_separate->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
//        s_dPhiIn_barrel_separate->Add(h_dPhiIn_barrel_MC_separate[pr]);
//        s_dPhiIn_endcap_separate->Add(h_dPhiIn_endcap_MC_separate[pr]);
//        s_HoverE_barrel_separate->Add(h_HoverE_barrel_MC_separate[pr]);
//        s_HoverE_endcap_separate->Add(h_HoverE_endcap_MC_separate[pr]);
//        s_InvEminusInvP_barrel_separate->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
//        s_InvEminusInvP_endcap_separate->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
//        s_mHits_barrel_separate->Add(h_mHits_barrel_MC_separate[pr]);
//        s_mHits_endcap_separate->Add(h_mHits_endcap_MC_separate[pr]);
//        s_passConvVeto_barrel_separate->Add(h_passConvVeto_barrel_MC_separate[pr]);
//        s_passConvVeto_endcap_separate->Add(h_passConvVeto_endcap_MC_separate[pr]);
//        s_passMediumID_barrel_separate->Add(h_passMediumID_barrel_MC_separate[pr]);
//        s_passMediumID_endcap_separate->Add(h_passMediumID_endcap_MC_separate[pr]);

//        for (Int_t ih=0; ih<nPtBin_ele; ih++)
//        {
//            s_PFiso_Rho_barrel_template[ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
//            s_PFiso_Rho_barrel_jetTemplate[ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);
//            if (ih < nPtBin_ele)
//            {
//                s_PFiso_Rho_endcap_template[ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
//                s_PFiso_Rho_endcap_jetTemplate[ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
//            }
//        }

        file->Close();
    }

    cout << "QCD received" << endl;

    // GammaJets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

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
        file->GetObject("h_TrkIso_barrel_nume", h_TrkIso_barrel_MC_nume[pr]);
        file->GetObject("h_TrkIso_endcap_nume", h_TrkIso_endcap_MC_nume[pr]);
        file->GetObject("h_TrkIso_barrel_deno", h_TrkIso_barrel_MC_deno[pr]);
        file->GetObject("h_TrkIso_endcap_deno", h_TrkIso_endcap_MC_deno[pr]);
        file->GetObject("h_TrkIso_barrel_ctrl", h_TrkIso_barrel_MC_ctrl[pr]);
        file->GetObject("h_TrkIso_endcap_ctrl", h_TrkIso_endcap_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_barrel_nume", h_ECALiso_barrel_MC_nume[pr]);
        file->GetObject("h_ECALiso_endcap_nume", h_ECALiso_endcap_MC_nume[pr]);
        file->GetObject("h_ECALiso_barrel_deno", h_ECALiso_barrel_MC_deno[pr]);
        file->GetObject("h_ECALiso_endcap_deno", h_ECALiso_endcap_MC_deno[pr]);
        file->GetObject("h_ECALiso_barrel_ctrl", h_ECALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_ECALiso_endcap_ctrl", h_ECALiso_endcap_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_barrel_nume", h_HCALiso_barrel_MC_nume[pr]);
        file->GetObject("h_HCALiso_endcap_nume", h_HCALiso_endcap_MC_nume[pr]);
        file->GetObject("h_HCALiso_barrel_deno", h_HCALiso_barrel_MC_deno[pr]);
        file->GetObject("h_HCALiso_endcap_deno", h_HCALiso_endcap_MC_deno[pr]);
        file->GetObject("h_HCALiso_barrel_ctrl", h_HCALiso_barrel_MC_ctrl[pr]);
        file->GetObject("h_HCALiso_endcap_ctrl", h_HCALiso_endcap_MC_ctrl[pr]);
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
        file->GetObject("h_HoverE_barrel_template_int", h_HoverE_barrel_template_int_MC[pr]);
        file->GetObject("h_HoverE_barrel_jetTemplate_int", h_HoverE_barrel_jetTemplate_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_template_int", h_HoverE_endcap_template_int_MC[pr]);
        file->GetObject("h_HoverE_endcap_jetTemplate_int", h_HoverE_endcap_jetTemplate_int_MC[pr]);
        file->GetObject("h_PFiso_Rho_barrel_separate", h_PFiso_Rho_barrel_MC_separate[pr]);
        file->GetObject("h_PFiso_Rho_endcap_separate", h_PFiso_Rho_endcap_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_barrel_separate", h_dEtaInSeed_barrel_MC_separate[pr]);
        file->GetObject("h_dEtaInSeed_endcap_separate", h_dEtaInSeed_endcap_MC_separate[pr]);
        file->GetObject("h_dPhiIn_barrel_separate", h_dPhiIn_barrel_MC_separate[pr]);
        file->GetObject("h_dPhiIn_endcap_separate", h_dPhiIn_endcap_MC_separate[pr]);
        file->GetObject("h_HoverE_barrel_separate", h_HoverE_barrel_MC_separate[pr]);
        file->GetObject("h_HoverE_endcap_separate", h_HoverE_endcap_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_barrel_separate", h_InvEminusInvP_barrel_MC_separate[pr]);
        file->GetObject("h_InvEminusInvP_endcap_separate", h_InvEminusInvP_endcap_MC_separate[pr]);
        file->GetObject("h_mHits_barrel_separate", h_mHits_barrel_MC_separate[pr]);
        file->GetObject("h_mHits_endcap_separate", h_mHits_endcap_MC_separate[pr]);
        file->GetObject("h_passConvVeto_barrel_separate", h_passConvVeto_barrel_MC_separate[pr]);
        file->GetObject("h_passConvVeto_endcap_separate", h_passConvVeto_endcap_MC_separate[pr]);
        file->GetObject("h_passMediumID_barrel_separate", h_passMediumID_barrel_MC_separate[pr]);
        file->GetObject("h_passMediumID_endcap_separate", h_passMediumID_endcap_MC_separate[pr]);
        h_pT_barrel_MC_deno_density[pr] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_"+fm.Procname[pr])));
        h_pT_endcap_MC_deno_density[pr] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_"+fm.Procname[pr])));

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
        removeNegativeBins(h_TrkIso_barrel_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_nume[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_deno[pr]);
        removeNegativeBins(h_TrkIso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_TrkIso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_ECALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_ECALiso_endcap_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_nume[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_deno[pr]);
        removeNegativeBins(h_HCALiso_barrel_MC_ctrl[pr]);
        removeNegativeBins(h_HCALiso_endcap_MC_ctrl[pr]);
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
        removeNegativeBins(h_HoverE_barrel_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_barrel_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_template_int_MC[pr]);
        removeNegativeBins(h_HoverE_endcap_jetTemplate_int_MC[pr]);
        removeNegativeBins(h_PFiso_Rho_barrel_MC_separate[pr]);
        removeNegativeBins(h_PFiso_Rho_endcap_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        removeNegativeBins(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_barrel_MC_separate[pr]);
        removeNegativeBins(h_dEtaInSeed_endcap_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_barrel_MC_separate[pr]);
        removeNegativeBins(h_dPhiIn_endcap_MC_separate[pr]);
        removeNegativeBins(h_HoverE_barrel_MC_separate[pr]);
        removeNegativeBins(h_HoverE_endcap_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_barrel_MC_separate[pr]);
        removeNegativeBins(h_InvEminusInvP_endcap_MC_separate[pr]);
        removeNegativeBins(h_mHits_barrel_MC_separate[pr]);
        removeNegativeBins(h_mHits_endcap_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_barrel_MC_separate[pr]);
        removeNegativeBins(h_passConvVeto_endcap_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_barrel_MC_separate[pr]);
        removeNegativeBins(h_passMediumID_endcap_MC_separate[pr]);
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr]);

        for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                 (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
            h_pT_endcap_MC_deno_density[pr]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        }

        Color_t color = kYellow + 3;
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
        h_TrkIso_barrel_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetFillColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetFillColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetFillColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetFillColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetFillColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetFillColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetFillColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetFillColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetFillColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetFillColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetFillColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetFillColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetFillColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetFillColor(color);
        h_mHits_barrel_MC_separate[pr]->SetFillColor(color);
        h_mHits_endcap_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetFillColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetFillColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetFillColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetFillColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_nume[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_deno[pr]->SetLineColor(color);
        h_TrkIso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_TrkIso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_ECALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_ECALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_nume[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_deno[pr]->SetLineColor(color);
        h_HCALiso_barrel_MC_ctrl[pr]->SetLineColor(color);
        h_HCALiso_endcap_MC_ctrl[pr]->SetLineColor(color);
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
        h_HoverE_barrel_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_template_int_MC[pr]->SetLineColor(color);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetLineColor(color);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetLineColor(color);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetLineColor(color);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetLineColor(color);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_barrel_MC_separate[pr]->SetLineColor(color);
        h_dPhiIn_endcap_MC_separate[pr]->SetLineColor(color);
        h_HoverE_barrel_MC_separate[pr]->SetLineColor(color);
        h_HoverE_endcap_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetLineColor(color);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetLineColor(color);
        h_mHits_barrel_MC_separate[pr]->SetLineColor(color);
        h_mHits_endcap_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_barrel_MC_separate[pr]->SetLineColor(color);
        h_passConvVeto_endcap_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_barrel_MC_separate[pr]->SetLineColor(color);
        h_passMediumID_endcap_MC_separate[pr]->SetLineColor(color);
        h_pT_barrel_MC_deno_density[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetLineColor(color);

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
        h_TrkIso_barrel_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_nume[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_deno[pr]->SetDirectory(0);
        h_TrkIso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_TrkIso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_ECALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_ECALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_nume[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_deno[pr]->SetDirectory(0);
        h_HCALiso_barrel_MC_ctrl[pr]->SetDirectory(0);
        h_HCALiso_endcap_MC_ctrl[pr]->SetDirectory(0);
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
        h_HoverE_barrel_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_barrel_jetTemplate_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_template_int_MC[pr]->SetDirectory(0);
        h_HoverE_endcap_jetTemplate_int_MC[pr]->SetDirectory(0);
        h_PFiso_Rho_barrel_MC_separate[pr]->SetDirectory(0);
        h_PFiso_Rho_endcap_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_barrel_MC_separate[pr]->SetDirectory(0);
        h_SigmaIEtaIEta_endcap_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_barrel_MC_separate[pr]->SetDirectory(0);
        h_dEtaInSeed_endcap_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_barrel_MC_separate[pr]->SetDirectory(0);
        h_dPhiIn_endcap_MC_separate[pr]->SetDirectory(0);
        h_HoverE_barrel_MC_separate[pr]->SetDirectory(0);
        h_HoverE_endcap_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_barrel_MC_separate[pr]->SetDirectory(0);
        h_InvEminusInvP_endcap_MC_separate[pr]->SetDirectory(0);
        h_mHits_barrel_MC_separate[pr]->SetDirectory(0);
        h_mHits_endcap_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_barrel_MC_separate[pr]->SetDirectory(0);
        h_passConvVeto_endcap_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_barrel_MC_separate[pr]->SetDirectory(0);
        h_passMediumID_endcap_MC_separate[pr]->SetDirectory(0);
        h_pT_barrel_MC_deno_density[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr]->SetDirectory(0);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_template[pr][ih]);
            file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            removeNegativeBins(h_PFiso_Rho_barrel_template[pr][ih]);
            removeNegativeBins(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetFillColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetFillColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetLineColor(color);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetLineColor(color);

            h_PFiso_Rho_barrel_template[pr][ih]   ->SetDirectory(0);
            h_PFiso_Rho_barrel_jetTemplate[pr][ih]->SetDirectory(0);

            if (ih < nPtBin_ele)
            {
                file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_template[pr][ih]);
                file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                removeNegativeBins(h_PFiso_Rho_endcap_template[pr][ih]);
                removeNegativeBins(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetFillColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetFillColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetLineColor(color);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetLineColor(color);

                h_PFiso_Rho_endcap_template[pr][ih]   ->SetDirectory(0);
                h_PFiso_Rho_endcap_jetTemplate[pr][ih]->SetDirectory(0);
            }
        }

        if (pr == _GJets_20to100)
        {
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
            h_TrkIso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_nume[pr]->Clone("h_TrkIso_barrel_nume_GJets")));
            h_TrkIso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_nume[pr]->Clone("h_TrkIso_endcap_nume_GJets")));
            h_TrkIso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_deno[pr]->Clone("h_TrkIso_barrel_deno_GJets")));
            h_TrkIso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_deno[pr]->Clone("h_TrkIso_endcap_deno_GJets")));
            h_TrkIso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_TrkIso_barrel_MC_ctrl[pr]->Clone("h_TrkIso_barrel_ctrl_GJets")));
            h_TrkIso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_TrkIso_endcap_MC_ctrl[pr]->Clone("h_TrkIso_endcap_ctrl_GJets")));
            h_ECALiso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_nume[pr]->Clone("h_ECALiso_barrel_nume_GJets")));
            h_ECALiso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_nume[pr]->Clone("h_ECALiso_endcap_nume_GJets")));
            h_ECALiso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_deno[pr]->Clone("h_ECALiso_barrel_deno_GJets")));
            h_ECALiso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_deno[pr]->Clone("h_ECALiso_endcap_deno_GJets")));
            h_ECALiso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_ECALiso_barrel_MC_ctrl[pr]->Clone("h_ECALiso_barrel_ctrl_GJets")));
            h_ECALiso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_ECALiso_endcap_MC_ctrl[pr]->Clone("h_ECALiso_endcap_ctrl_GJets")));
            h_HCALiso_barrel_MC_nume[_GJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_nume[pr]->Clone("h_HCALiso_barrel_nume_GJets")));
            h_HCALiso_endcap_MC_nume[_GJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_nume[pr]->Clone("h_HCALiso_endcap_nume_GJets")));
            h_HCALiso_barrel_MC_deno[_GJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_deno[pr]->Clone("h_HCALiso_barrel_deno_GJets")));
            h_HCALiso_endcap_MC_deno[_GJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_deno[pr]->Clone("h_HCALiso_endcap_deno_GJets")));
            h_HCALiso_barrel_MC_ctrl[_GJets_Full] = ((TH1D*)(h_HCALiso_barrel_MC_ctrl[pr]->Clone("h_HCALiso_barrel_ctrl_GJets")));
            h_HCALiso_endcap_MC_ctrl[_GJets_Full] = ((TH1D*)(h_HCALiso_endcap_MC_ctrl[pr]->Clone("h_HCALiso_endcap_ctrl_GJets")));
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
            h_HoverE_barrel_template_int_MC[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_template_int_MC[pr]->Clone("h_HoverE_barrel_template_int_GJets")));
            h_HoverE_barrel_jetTemplate_int_MC[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_jetTemplate_int_MC[pr]->Clone("h_HoverE_barrel_jetTemplate_int_GJets")));
            h_HoverE_endcap_template_int_MC[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_template_int_MC[pr]->Clone("h_HoverE_endcap_template_int_GJets")));
            h_HoverE_endcap_jetTemplate_int_MC[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_jetTemplate_int_MC[pr]->Clone("h_HoverE_endcap_jetTemplate_int_GJets")));
            h_PFiso_Rho_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_barrel_MC_separate[pr]->Clone("h_PFiso_Rho_barrel_separate_GJets")));
            h_PFiso_Rho_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_PFiso_Rho_endcap_MC_separate[pr]->Clone("h_PFiso_Rho_endcap_separate_GJets")));
            h_SigmaIEtaIEta_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_barrel_MC_separate[pr]->Clone("h_SigmaIEtaIEta_barrel_separate_GJets")));
            h_SigmaIEtaIEta_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_SigmaIEtaIEta_endcap_MC_separate[pr]->Clone("h_SigmaIEtaIEta_endcap_separate_GJets")));
            h_dEtaInSeed_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_barrel_MC_separate[pr]->Clone("h_dEtaInSeed_barrel_separate_GJets")));
            h_dEtaInSeed_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_dEtaInSeed_endcap_MC_separate[pr]->Clone("h_dEtaInSeed_endcap_separate_GJets")));
            h_dPhiIn_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_dPhiIn_barrel_MC_separate[pr]->Clone("h_dPhiIn_barrel_separate_GJets")));
            h_dPhiIn_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_dPhiIn_endcap_MC_separate[pr]->Clone("h_dPhiIn_endcap_separate_GJets")));
            h_HoverE_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_HoverE_barrel_MC_separate[pr]->Clone("h_HoverE_barrel_separate_GJets")));
            h_HoverE_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_HoverE_endcap_MC_separate[pr]->Clone("h_HoverE_endcap_separate_GJets")));
            h_InvEminusInvP_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_barrel_MC_separate[pr]->Clone("h_InvEminusInvP_barrel_separate_GJets")));
            h_InvEminusInvP_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_InvEminusInvP_endcap_MC_separate[pr]->Clone("h_InvEminusInvP_endcap_separate_GJets")));
            h_mHits_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_mHits_barrel_MC_separate[pr]->Clone("h_mHits_barrel_separate_GJets")));
            h_mHits_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_mHits_endcap_MC_separate[pr]->Clone("h_mHits_endcap_separate_GJets")));
            h_passConvVeto_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_passConvVeto_barrel_MC_separate[pr]->Clone("h_passConvVeto_barrel_separate_GJets")));
            h_passConvVeto_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_passConvVeto_endcap_MC_separate[pr]->Clone("h_passConvVeto_endcap_separate_GJets")));
            h_passMediumID_barrel_MC_separate[_GJets_Full] = ((TH1D*)(h_passMediumID_barrel_MC_separate[pr]->Clone("h_passMediumID_barrel_separate_GJets")));
            h_passMediumID_endcap_MC_separate[_GJets_Full] = ((TH1D*)(h_passMediumID_endcap_MC_separate[pr]->Clone("h_passMediumID_endcap_separate_GJets")));
            h_pT_barrel_MC_deno_density[_GJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr]->Clone("h_pT_barrel_deno_density_GJets")));
            h_pT_endcap_MC_deno_density[_GJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr]->Clone("h_pT_endcap_deno_density_GJets")));

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
            h_TrkIso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_TrkIso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_TrkIso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_ECALiso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_ECALiso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_nume[_GJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_nume[_GJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_deno[_GJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_deno[_GJets_Full]->SetDirectory(0);
            h_HCALiso_barrel_MC_ctrl[_GJets_Full]->SetDirectory(0);
            h_HCALiso_endcap_MC_ctrl[_GJets_Full]->SetDirectory(0);
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
            h_HoverE_barrel_template_int_MC[_GJets_Full]->SetDirectory(0);
            h_HoverE_barrel_jetTemplate_int_MC[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_template_int_MC[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_jetTemplate_int_MC[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_PFiso_Rho_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_SigmaIEtaIEta_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_dEtaInSeed_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_dPhiIn_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_HoverE_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_HoverE_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_InvEminusInvP_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_mHits_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_mHits_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_passConvVeto_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_passConvVeto_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_passMediumID_barrel_MC_separate[_GJets_Full]->SetDirectory(0);
            h_passMediumID_endcap_MC_separate[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_MC_deno_density[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_GJets_Full]->SetDirectory(0);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_GJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_template[pr][ih]->Clone("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10)+"_GJets")));
                h_PFiso_Rho_barrel_jetTemplate[_GJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_barrel_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10)+"_GJets")));
                h_PFiso_Rho_barrel_template[_GJets_Full][ih]->SetDirectory(0);
                h_PFiso_Rho_barrel_jetTemplate[_GJets_Full][ih]->SetDirectory(0);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_GJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_template[pr][ih]->Clone("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10)+"_GJets")));
                    h_PFiso_Rho_endcap_jetTemplate[_GJets_Full][ih] = ((TH1D*)(h_PFiso_Rho_endcap_jetTemplate[pr][ih]->Clone("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10)+"_GJets")));
                    h_PFiso_Rho_endcap_template[_GJets_Full][ih]->SetDirectory(0);
                    h_PFiso_Rho_endcap_jetTemplate[_GJets_Full][ih]->SetDirectory(0);
                }
            }
        }
        else
        {
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
            h_TrkIso_barrel_MC_nume[_GJets_Full]->Add(h_TrkIso_barrel_MC_nume[pr]);
            h_TrkIso_endcap_MC_nume[_GJets_Full]->Add(h_TrkIso_endcap_MC_nume[pr]);
            h_TrkIso_barrel_MC_deno[_GJets_Full]->Add(h_TrkIso_barrel_MC_deno[pr]);
            h_TrkIso_endcap_MC_deno[_GJets_Full]->Add(h_TrkIso_endcap_MC_deno[pr]);
            h_TrkIso_barrel_MC_ctrl[_GJets_Full]->Add(h_TrkIso_barrel_MC_ctrl[pr]);
            h_TrkIso_endcap_MC_ctrl[_GJets_Full]->Add(h_TrkIso_endcap_MC_ctrl[pr]);
            h_ECALiso_barrel_MC_nume[_GJets_Full]->Add(h_ECALiso_barrel_MC_nume[pr]);
            h_ECALiso_endcap_MC_nume[_GJets_Full]->Add(h_ECALiso_endcap_MC_nume[pr]);
            h_ECALiso_barrel_MC_deno[_GJets_Full]->Add(h_ECALiso_barrel_MC_deno[pr]);
            h_ECALiso_endcap_MC_deno[_GJets_Full]->Add(h_ECALiso_endcap_MC_deno[pr]);
            h_ECALiso_barrel_MC_ctrl[_GJets_Full]->Add(h_ECALiso_barrel_MC_ctrl[pr]);
            h_ECALiso_endcap_MC_ctrl[_GJets_Full]->Add(h_ECALiso_endcap_MC_ctrl[pr]);
            h_HCALiso_barrel_MC_nume[_GJets_Full]->Add(h_HCALiso_barrel_MC_nume[pr]);
            h_HCALiso_endcap_MC_nume[_GJets_Full]->Add(h_HCALiso_endcap_MC_nume[pr]);
            h_HCALiso_barrel_MC_deno[_GJets_Full]->Add(h_HCALiso_barrel_MC_deno[pr]);
            h_HCALiso_endcap_MC_deno[_GJets_Full]->Add(h_HCALiso_endcap_MC_deno[pr]);
            h_HCALiso_barrel_MC_ctrl[_GJets_Full]->Add(h_HCALiso_barrel_MC_ctrl[pr]);
            h_HCALiso_endcap_MC_ctrl[_GJets_Full]->Add(h_HCALiso_endcap_MC_ctrl[pr]);
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
            h_HoverE_barrel_template_int_MC[_GJets_Full]->Add(h_HoverE_barrel_template_int_MC[pr]);
            h_HoverE_barrel_jetTemplate_int_MC[_GJets_Full]->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
            h_HoverE_endcap_template_int_MC[_GJets_Full]->Add(h_HoverE_endcap_template_int_MC[pr]);
            h_HoverE_endcap_jetTemplate_int_MC[_GJets_Full]->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
            h_PFiso_Rho_barrel_MC_separate[_GJets_Full]->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
            h_PFiso_Rho_endcap_MC_separate[_GJets_Full]->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
            h_SigmaIEtaIEta_barrel_MC_separate[_GJets_Full]->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
            h_SigmaIEtaIEta_endcap_MC_separate[_GJets_Full]->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
            h_dEtaInSeed_barrel_MC_separate[_GJets_Full]->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
            h_dEtaInSeed_endcap_MC_separate[_GJets_Full]->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
            h_dPhiIn_barrel_MC_separate[_GJets_Full]->Add(h_dPhiIn_barrel_MC_separate[pr]);
            h_dPhiIn_endcap_MC_separate[_GJets_Full]->Add(h_dPhiIn_endcap_MC_separate[pr]);
            h_HoverE_barrel_MC_separate[_GJets_Full]->Add(h_HoverE_barrel_MC_separate[pr]);
            h_HoverE_endcap_MC_separate[_GJets_Full]->Add(h_HoverE_endcap_MC_separate[pr]);
            h_InvEminusInvP_barrel_MC_separate[_GJets_Full]->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
            h_InvEminusInvP_endcap_MC_separate[_GJets_Full]->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
            h_mHits_barrel_MC_separate[_GJets_Full]->Add(h_mHits_barrel_MC_separate[pr]);
            h_mHits_endcap_MC_separate[_GJets_Full]->Add(h_mHits_endcap_MC_separate[pr]);
            h_passConvVeto_barrel_MC_separate[_GJets_Full]->Add(h_passConvVeto_barrel_MC_separate[pr]);
            h_passConvVeto_endcap_MC_separate[_GJets_Full]->Add(h_passConvVeto_endcap_MC_separate[pr]);
            h_passMediumID_barrel_MC_separate[_GJets_Full]->Add(h_passMediumID_barrel_MC_separate[pr]);
            h_passMediumID_endcap_MC_separate[_GJets_Full]->Add(h_passMediumID_endcap_MC_separate[pr]);
            h_pT_barrel_MC_deno_density[_GJets_Full]->Add(h_pT_barrel_MC_deno_density[pr]);
            h_pT_endcap_MC_deno_density[_GJets_Full]->Add(h_pT_endcap_MC_deno_density[pr]);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                h_PFiso_Rho_barrel_template[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
                h_PFiso_Rho_barrel_jetTemplate[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);

                if (ih < nPtBin_ele)
                {
                    h_PFiso_Rho_endcap_template[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
                    h_PFiso_Rho_endcap_jetTemplate[_QCDEMEnriched_Full][ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
                }
            }
        }

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
        s_TrkIso_barrel_nume->Add(h_TrkIso_barrel_MC_nume[pr]);
        s_TrkIso_endcap_nume->Add(h_TrkIso_endcap_MC_nume[pr]);
        s_TrkIso_barrel_deno->Add(h_TrkIso_barrel_MC_deno[pr]);
        s_TrkIso_endcap_deno->Add(h_TrkIso_endcap_MC_deno[pr]);
        s_TrkIso_barrel_ctrl->Add(h_TrkIso_barrel_MC_ctrl[pr]);
        s_TrkIso_endcap_ctrl->Add(h_TrkIso_endcap_MC_ctrl[pr]);
        s_ECALiso_barrel_nume->Add(h_ECALiso_barrel_MC_nume[pr]);
        s_ECALiso_endcap_nume->Add(h_ECALiso_endcap_MC_nume[pr]);
        s_ECALiso_barrel_deno->Add(h_ECALiso_barrel_MC_deno[pr]);
        s_ECALiso_endcap_deno->Add(h_ECALiso_endcap_MC_deno[pr]);
        s_ECALiso_barrel_ctrl->Add(h_ECALiso_barrel_MC_ctrl[pr]);
        s_ECALiso_endcap_ctrl->Add(h_ECALiso_endcap_MC_ctrl[pr]);
        s_HCALiso_barrel_nume->Add(h_HCALiso_barrel_MC_nume[pr]);
        s_HCALiso_endcap_nume->Add(h_HCALiso_endcap_MC_nume[pr]);
        s_HCALiso_barrel_deno->Add(h_HCALiso_barrel_MC_deno[pr]);
        s_HCALiso_endcap_deno->Add(h_HCALiso_endcap_MC_deno[pr]);
        s_HCALiso_barrel_ctrl->Add(h_HCALiso_barrel_MC_ctrl[pr]);
        s_HCALiso_endcap_ctrl->Add(h_HCALiso_endcap_MC_ctrl[pr]);
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
        s_HoverE_barrel_template_int->Add(h_HoverE_barrel_template_int_MC[pr]);
        s_HoverE_barrel_jetTemplate_int->Add(h_HoverE_barrel_jetTemplate_int_MC[pr]);
        s_HoverE_endcap_template_int->Add(h_HoverE_endcap_template_int_MC[pr]);
        s_HoverE_endcap_jetTemplate_int->Add(h_HoverE_endcap_jetTemplate_int_MC[pr]);
        s_PFiso_Rho_barrel_separate->Add(h_PFiso_Rho_barrel_MC_separate[pr]);
        s_PFiso_Rho_endcap_separate->Add(h_PFiso_Rho_endcap_MC_separate[pr]);
        s_SigmaIEtaIEta_barrel_separate->Add(h_SigmaIEtaIEta_barrel_MC_separate[pr]);
        s_SigmaIEtaIEta_endcap_separate->Add(h_SigmaIEtaIEta_endcap_MC_separate[pr]);
        s_dEtaInSeed_barrel_separate->Add(h_dEtaInSeed_barrel_MC_separate[pr]);
        s_dEtaInSeed_endcap_separate->Add(h_dEtaInSeed_endcap_MC_separate[pr]);
        s_dPhiIn_barrel_separate->Add(h_dPhiIn_barrel_MC_separate[pr]);
        s_dPhiIn_endcap_separate->Add(h_dPhiIn_endcap_MC_separate[pr]);
        s_HoverE_barrel_separate->Add(h_HoverE_barrel_MC_separate[pr]);
        s_HoverE_endcap_separate->Add(h_HoverE_endcap_MC_separate[pr]);
        s_InvEminusInvP_barrel_separate->Add(h_InvEminusInvP_barrel_MC_separate[pr]);
        s_InvEminusInvP_endcap_separate->Add(h_InvEminusInvP_endcap_MC_separate[pr]);
        s_mHits_barrel_separate->Add(h_mHits_barrel_MC_separate[pr]);
        s_mHits_endcap_separate->Add(h_mHits_endcap_MC_separate[pr]);
        s_passConvVeto_barrel_separate->Add(h_passConvVeto_barrel_MC_separate[pr]);
        s_passConvVeto_endcap_separate->Add(h_passConvVeto_endcap_MC_separate[pr]);
        s_passMediumID_barrel_separate->Add(h_passMediumID_barrel_MC_separate[pr]);
        s_passMediumID_endcap_separate->Add(h_passMediumID_endcap_MC_separate[pr]);
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr]);

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            s_PFiso_Rho_barrel_template[ih]->Add(h_PFiso_Rho_barrel_template[pr][ih]);
            s_PFiso_Rho_barrel_jetTemplate[ih]->Add(h_PFiso_Rho_barrel_jetTemplate[pr][ih]);
            if (ih < nPtBin_ele)
            {
                s_PFiso_Rho_endcap_template[ih]->Add(h_PFiso_Rho_endcap_template[pr][ih]);
                s_PFiso_Rho_endcap_jetTemplate[ih]->Add(h_PFiso_Rho_endcap_jetTemplate[pr][ih]);
            }
        }

        file->Close();
    }

    cout << "GJets received" << endl;

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");
        TH1D *h_temp[93];
        if (pr == _SinglePhoton_B)
        {
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
            file->GetObject("h_TrkIso_barrel_nume", h_TrkIso_barrel_data_nume);
            file->GetObject("h_TrkIso_endcap_nume", h_TrkIso_endcap_data_nume);
            file->GetObject("h_TrkIso_barrel_deno", h_TrkIso_barrel_data_deno);
            file->GetObject("h_TrkIso_endcap_deno", h_TrkIso_endcap_data_deno);
            file->GetObject("h_TrkIso_barrel_ctrl", h_TrkIso_barrel_data_ctrl);
            file->GetObject("h_TrkIso_endcap_ctrl", h_TrkIso_endcap_data_ctrl);
            file->GetObject("h_ECALiso_barrel_nume", h_ECALiso_barrel_data_nume);
            file->GetObject("h_ECALiso_endcap_nume", h_ECALiso_endcap_data_nume);
            file->GetObject("h_ECALiso_barrel_deno", h_ECALiso_barrel_data_deno);
            file->GetObject("h_ECALiso_endcap_deno", h_ECALiso_endcap_data_deno);
            file->GetObject("h_ECALiso_barrel_ctrl", h_ECALiso_barrel_data_ctrl);
            file->GetObject("h_ECALiso_endcap_ctrl", h_ECALiso_endcap_data_ctrl);
            file->GetObject("h_HCALiso_barrel_nume", h_HCALiso_barrel_data_nume);
            file->GetObject("h_HCALiso_endcap_nume", h_HCALiso_endcap_data_nume);
            file->GetObject("h_HCALiso_barrel_deno", h_HCALiso_barrel_data_deno);
            file->GetObject("h_HCALiso_endcap_deno", h_HCALiso_endcap_data_deno);
            file->GetObject("h_HCALiso_barrel_ctrl", h_HCALiso_barrel_data_ctrl);
            file->GetObject("h_HCALiso_endcap_ctrl", h_HCALiso_endcap_data_ctrl);
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
            file->GetObject("h_HoverE_barrel_template_int", h_HoverE_barrel_template_int_data);
            file->GetObject("h_HoverE_barrel_jetTemplate_int", h_HoverE_barrel_jetTemplate_int_data);
            file->GetObject("h_HoverE_endcap_template_int", h_HoverE_endcap_template_int_data);
            file->GetObject("h_HoverE_endcap_jetTemplate_int", h_HoverE_endcap_jetTemplate_int_data);
            file->GetObject("h_PFiso_Rho_barrel_separate", h_PFiso_Rho_barrel_data_separate);
            file->GetObject("h_PFiso_Rho_endcap_separate", h_PFiso_Rho_endcap_data_separate);
            file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_SigmaIEtaIEta_barrel_data_separate);
            file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_SigmaIEtaIEta_endcap_data_separate);
            file->GetObject("h_dEtaInSeed_barrel_separate", h_dEtaInSeed_barrel_data_separate);
            file->GetObject("h_dEtaInSeed_endcap_separate", h_dEtaInSeed_endcap_data_separate);
            file->GetObject("h_dPhiIn_barrel_separate", h_dPhiIn_barrel_data_separate);
            file->GetObject("h_dPhiIn_endcap_separate", h_dPhiIn_endcap_data_separate);
            file->GetObject("h_HoverE_barrel_separate", h_HoverE_barrel_data_separate);
            file->GetObject("h_HoverE_endcap_separate", h_HoverE_endcap_data_separate);
            file->GetObject("h_InvEminusInvP_barrel_separate", h_InvEminusInvP_barrel_data_separate);
            file->GetObject("h_InvEminusInvP_endcap_separate", h_InvEminusInvP_endcap_data_separate);
            file->GetObject("h_mHits_barrel_separate", h_mHits_barrel_data_separate);
            file->GetObject("h_mHits_endcap_separate", h_mHits_endcap_data_separate);
            file->GetObject("h_passConvVeto_barrel_separate", h_passConvVeto_barrel_data_separate);
            file->GetObject("h_passConvVeto_endcap_separate", h_passConvVeto_endcap_data_separate);
            file->GetObject("h_passMediumID_barrel_separate", h_passMediumID_barrel_data_separate);
            file->GetObject("h_passMediumID_endcap_separate", h_passMediumID_endcap_data_separate);
            h_pT_barrel_data_deno_density = ((TH1D*)(h_pT_barrel_data_deno->Clone("h_pT_barrel_data_deno_density")));
            h_pT_endcap_data_deno_density = ((TH1D*)(h_pT_endcap_data_deno->Clone("h_pT_endcap_data_deno_density")));

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
            removeNegativeBins(h_TrkIso_barrel_data_nume);
            removeNegativeBins(h_TrkIso_endcap_data_nume);
            removeNegativeBins(h_TrkIso_barrel_data_deno);
            removeNegativeBins(h_TrkIso_endcap_data_deno);
            removeNegativeBins(h_TrkIso_barrel_data_ctrl);
            removeNegativeBins(h_TrkIso_endcap_data_ctrl);

            removeNegativeBins(h_ECALiso_barrel_data_nume);
            removeNegativeBins(h_ECALiso_endcap_data_nume);
            removeNegativeBins(h_ECALiso_barrel_data_deno);
            removeNegativeBins(h_ECALiso_endcap_data_deno);
            removeNegativeBins(h_ECALiso_barrel_data_ctrl);
            removeNegativeBins(h_ECALiso_endcap_data_ctrl);

            removeNegativeBins(h_HCALiso_barrel_data_nume);
            removeNegativeBins(h_HCALiso_endcap_data_nume);
            removeNegativeBins(h_HCALiso_barrel_data_deno);
            removeNegativeBins(h_HCALiso_endcap_data_deno);
            removeNegativeBins(h_HCALiso_barrel_data_ctrl);
            removeNegativeBins(h_HCALiso_endcap_data_ctrl);
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
            removeNegativeBins(h_HoverE_barrel_template_int_data);
            removeNegativeBins(h_HoverE_barrel_jetTemplate_int_data);
            removeNegativeBins(h_HoverE_endcap_template_int_data);
            removeNegativeBins(h_HoverE_endcap_jetTemplate_int_data);
            removeNegativeBins(h_PFiso_Rho_barrel_data_separate);
            removeNegativeBins(h_PFiso_Rho_endcap_data_separate);
            removeNegativeBins(h_SigmaIEtaIEta_barrel_data_separate);
            removeNegativeBins(h_SigmaIEtaIEta_endcap_data_separate);
            removeNegativeBins(h_dEtaInSeed_barrel_data_separate);
            removeNegativeBins(h_dEtaInSeed_endcap_data_separate);
            removeNegativeBins(h_dPhiIn_barrel_data_separate);
            removeNegativeBins(h_dPhiIn_endcap_data_separate);
            removeNegativeBins(h_HoverE_barrel_data_separate);
            removeNegativeBins(h_HoverE_endcap_data_separate);
            removeNegativeBins(h_InvEminusInvP_barrel_data_separate);
            removeNegativeBins(h_InvEminusInvP_endcap_data_separate);
            removeNegativeBins(h_mHits_barrel_data_separate);
            removeNegativeBins(h_mHits_endcap_data_separate);
            removeNegativeBins(h_passConvVeto_barrel_data_separate);
            removeNegativeBins(h_passConvVeto_endcap_data_separate);
            removeNegativeBins(h_passMediumID_barrel_data_separate);
            removeNegativeBins(h_passMediumID_endcap_data_separate);
            removeNegativeBins(h_pT_barrel_data_deno_density);
            removeNegativeBins(h_pT_endcap_data_deno_density);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_data_template[ih]);
                file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_barrel_data_jetTemplate[ih]);
                removeNegativeBins(h_PFiso_Rho_barrel_data_template[ih]);
                removeNegativeBins(h_PFiso_Rho_barrel_data_jetTemplate[ih]);

                if (ih < nPtBin_ele)
                {
                    file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_data_template[ih]);
                    file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_PFiso_Rho_endcap_data_jetTemplate[ih]);
                    removeNegativeBins(h_PFiso_Rho_endcap_data_template[ih]);
                    removeNegativeBins(h_PFiso_Rho_endcap_data_jetTemplate[ih]);
                }
            }
        }
        else
        {
            file->GetObject("h_PFiso_Rho_barrel_nume", h_temp[1]);
            file->GetObject("h_PFiso_Rho_endcap_nume", h_temp[2]);
            file->GetObject("h_PFiso_Rho_barrel_deno", h_temp[3]);
            file->GetObject("h_PFiso_Rho_endcap_deno", h_temp[4]);
            file->GetObject("h_PFiso_Rho_barrel_ctrl", h_temp[5]);
            file->GetObject("h_PFiso_Rho_endcap_ctrl", h_temp[6]);
            file->GetObject("h_pT_barrel_nume", h_temp[7]);
            file->GetObject("h_pT_endcap_nume", h_temp[8]);
            file->GetObject("h_pT_barrel_deno", h_temp[9]);
            file->GetObject("h_pT_endcap_deno", h_temp[10]);
            file->GetObject("h_pT_barrel_ctrl", h_temp[11]);
            file->GetObject("h_pT_endcap_ctrl", h_temp[12]);
            file->GetObject("h_SigmaIEtaIEta_barrel_nume", h_temp[13]);
            file->GetObject("h_SigmaIEtaIEta_endcap_nume", h_temp[14]);
            file->GetObject("h_SigmaIEtaIEta_barrel_deno", h_temp[15]);
            file->GetObject("h_SigmaIEtaIEta_endcap_deno", h_temp[16]);
            file->GetObject("h_SigmaIEtaIEta_barrel_ctrl", h_temp[17]);
            file->GetObject("h_SigmaIEtaIEta_endcap_ctrl", h_temp[18]);
            file->GetObject("h_dEtaInSeed_barrel_nume", h_temp[19]);
            file->GetObject("h_dEtaInSeed_endcap_nume", h_temp[20]);
            file->GetObject("h_dEtaInSeed_barrel_deno", h_temp[21]);
            file->GetObject("h_dEtaInSeed_endcap_deno", h_temp[22]);
            file->GetObject("h_dEtaInSeed_barrel_ctrl", h_temp[23]);
            file->GetObject("h_dEtaInSeed_endcap_ctrl", h_temp[24]);
            file->GetObject("h_dPhiIn_barrel_nume", h_temp[25]);
            file->GetObject("h_dPhiIn_endcap_nume", h_temp[26]);
            file->GetObject("h_dPhiIn_barrel_deno", h_temp[27]);
            file->GetObject("h_dPhiIn_endcap_deno", h_temp[28]);
            file->GetObject("h_dPhiIn_barrel_ctrl", h_temp[29]);
            file->GetObject("h_dPhiIn_endcap_ctrl", h_temp[30]);
            file->GetObject("h_HoverE_barrel_nume", h_temp[31]);
            file->GetObject("h_HoverE_endcap_nume", h_temp[32]);
            file->GetObject("h_HoverE_barrel_deno", h_temp[33]);
            file->GetObject("h_HoverE_endcap_deno", h_temp[34]);
            file->GetObject("h_HoverE_barrel_ctrl", h_temp[35]);
            file->GetObject("h_HoverE_endcap_ctrl", h_temp[36]);
            file->GetObject("h_InvEminusInvP_barrel_nume", h_temp[37]);
            file->GetObject("h_InvEminusInvP_endcap_nume", h_temp[38]);
            file->GetObject("h_InvEminusInvP_barrel_deno", h_temp[39]);
            file->GetObject("h_InvEminusInvP_endcap_deno", h_temp[40]);
            file->GetObject("h_InvEminusInvP_barrel_ctrl", h_temp[41]);
            file->GetObject("h_InvEminusInvP_endcap_ctrl", h_temp[42]);
            file->GetObject("h_TrkIso_barrel_nume", h_temp[43]);
            file->GetObject("h_TrkIso_endcap_nume", h_temp[44]);
            file->GetObject("h_TrkIso_barrel_deno", h_temp[45]);
            file->GetObject("h_TrkIso_endcap_deno", h_temp[46]);
            file->GetObject("h_TrkIso_barrel_ctrl", h_temp[47]);
            file->GetObject("h_TrkIso_endcap_ctrl", h_temp[48]);
            file->GetObject("h_ECALiso_barrel_nume", h_temp[49]);
            file->GetObject("h_ECALiso_endcap_nume", h_temp[50]);
            file->GetObject("h_ECALiso_barrel_deno", h_temp[51]);
            file->GetObject("h_ECALiso_endcap_deno", h_temp[52]);
            file->GetObject("h_ECALiso_barrel_ctrl", h_temp[53]);
            file->GetObject("h_ECALiso_endcap_ctrl", h_temp[54]);
            file->GetObject("h_HCALiso_barrel_nume", h_temp[55]);
            file->GetObject("h_HCALiso_endcap_nume", h_temp[56]);
            file->GetObject("h_HCALiso_barrel_deno", h_temp[57]);
            file->GetObject("h_HCALiso_endcap_deno", h_temp[58]);
            file->GetObject("h_HCALiso_barrel_ctrl", h_temp[59]);
            file->GetObject("h_HCALiso_endcap_ctrl", h_temp[60]);
            file->GetObject("h_MET", h_temp[61]);
            file->GetObject("h_MT_barrel_nume", h_temp[62]);
            file->GetObject("h_MT_endcap_nume", h_temp[63]);
            file->GetObject("h_MT_barrel_deno", h_temp[64]);
            file->GetObject("h_MT_endcap_deno", h_temp[65]);
            file->GetObject("h_MT_barrel_ctrl", h_temp[66]);
            file->GetObject("h_MT_endcap_ctrl", h_temp[67]);
            file->GetObject("h_eta_deno", h_temp[68]);
            file->GetObject("h_nVTX", h_temp[69]);
            file->GetObject("h_mass_test", h_temp[70]);
            file->GetObject("h_HoverE_barrel_template_int", h_temp[71]);
            file->GetObject("h_HoverE_barrel_jetTemplate_int", h_temp[72]);
            file->GetObject("h_HoverE_endcap_template_int", h_temp[73]);
            file->GetObject("h_HoverE_endcap_jetTemplate_int", h_temp[74]);
            file->GetObject("h_PFiso_Rho_barrel_separate", h_temp[75]);
            file->GetObject("h_PFiso_Rho_endcap_separate", h_temp[76]);
            file->GetObject("h_SigmaIEtaIEta_barrel_separate", h_temp[77]);
            file->GetObject("h_SigmaIEtaIEta_endcap_separate", h_temp[78]);
            file->GetObject("h_dEtaInSeed_barrel_separate", h_temp[79]);
            file->GetObject("h_dEtaInSeed_endcap_separate", h_temp[80]);
            file->GetObject("h_dPhiIn_barrel_separate", h_temp[81]);
            file->GetObject("h_dPhiIn_endcap_separate", h_temp[82]);
            file->GetObject("h_HoverE_barrel_separate", h_temp[83]);
            file->GetObject("h_HoverE_endcap_separate", h_temp[84]);
            file->GetObject("h_InvEminusInvP_barrel_separate", h_temp[85]);
            file->GetObject("h_InvEminusInvP_endcap_separate", h_temp[86]);
            file->GetObject("h_mHits_barrel_separate", h_temp[87]);
            file->GetObject("h_mHits_endcap_separate", h_temp[88]);
            file->GetObject("h_passConvVeto_barrel_separate", h_temp[89]);
            file->GetObject("h_passConvVeto_endcap_separate", h_temp[90]);
            file->GetObject("h_passMediumID_barrel_separate", h_temp[91]);
            file->GetObject("h_passMediumID_endcap_separate", h_temp[92]);

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
            removeNegativeBins(h_temp[82]);
            removeNegativeBins(h_temp[83]);
            removeNegativeBins(h_temp[84]);
            removeNegativeBins(h_temp[85]);
            removeNegativeBins(h_temp[86]);
            removeNegativeBins(h_temp[87]);
            removeNegativeBins(h_temp[88]);
            removeNegativeBins(h_temp[89]);
            removeNegativeBins(h_temp[90]);
            removeNegativeBins(h_temp[91]);
            removeNegativeBins(h_temp[92]);

            h_PFiso_Rho_barrel_data_nume->Add(h_temp[1]);
            h_PFiso_Rho_endcap_data_nume->Add(h_temp[2]);
            h_PFiso_Rho_barrel_data_deno->Add(h_temp[3]);
            h_PFiso_Rho_endcap_data_deno->Add(h_temp[4]);
            h_PFiso_Rho_barrel_data_ctrl->Add(h_temp[5]);
            h_PFiso_Rho_endcap_data_ctrl->Add(h_temp[6]);
            h_pT_barrel_data_nume->Add(h_temp[7]);
            h_pT_endcap_data_nume->Add(h_temp[8]);
            h_pT_barrel_data_deno->Add(h_temp[9]);
            h_pT_endcap_data_deno->Add(h_temp[10]);
            h_pT_barrel_data_ctrl->Add(h_temp[11]);
            h_pT_endcap_data_ctrl->Add(h_temp[12]);
            h_SigmaIEtaIEta_barrel_data_nume->Add(h_temp[13]);
            h_SigmaIEtaIEta_endcap_data_nume->Add(h_temp[14]);
            h_SigmaIEtaIEta_barrel_data_deno->Add(h_temp[15]);
            h_SigmaIEtaIEta_endcap_data_deno->Add(h_temp[16]);
            h_SigmaIEtaIEta_barrel_data_ctrl->Add(h_temp[17]);
            h_SigmaIEtaIEta_endcap_data_ctrl->Add(h_temp[18]);
            h_dEtaInSeed_barrel_data_nume->Add(h_temp[19]);
            h_dEtaInSeed_endcap_data_nume->Add(h_temp[20]);
            h_dEtaInSeed_barrel_data_deno->Add(h_temp[21]);
            h_dEtaInSeed_endcap_data_deno->Add(h_temp[22]);
            h_dEtaInSeed_barrel_data_ctrl->Add(h_temp[23]);
            h_dEtaInSeed_endcap_data_ctrl->Add(h_temp[24]);
            h_dPhiIn_barrel_data_nume->Add(h_temp[25]);
            h_dPhiIn_endcap_data_nume->Add(h_temp[26]);
            h_dPhiIn_barrel_data_deno->Add(h_temp[27]);
            h_dPhiIn_endcap_data_deno->Add(h_temp[28]);
            h_dPhiIn_barrel_data_ctrl->Add(h_temp[29]);
            h_dPhiIn_endcap_data_ctrl->Add(h_temp[30]);
            h_HoverE_barrel_data_nume->Add(h_temp[31]);
            h_HoverE_endcap_data_nume->Add(h_temp[32]);
            h_HoverE_barrel_data_deno->Add(h_temp[33]);
            h_HoverE_endcap_data_deno->Add(h_temp[34]);
            h_HoverE_barrel_data_ctrl->Add(h_temp[35]);
            h_HoverE_endcap_data_ctrl->Add(h_temp[36]);
            h_InvEminusInvP_barrel_data_nume->Add(h_temp[37]);
            h_InvEminusInvP_endcap_data_nume->Add(h_temp[38]);
            h_InvEminusInvP_barrel_data_deno->Add(h_temp[39]);
            h_InvEminusInvP_endcap_data_deno->Add(h_temp[40]);
            h_InvEminusInvP_barrel_data_ctrl->Add(h_temp[41]);
            h_InvEminusInvP_endcap_data_ctrl->Add(h_temp[42]);
            h_TrkIso_barrel_data_nume->Add(h_temp[43]);
            h_TrkIso_endcap_data_nume->Add(h_temp[44]);
            h_TrkIso_barrel_data_deno->Add(h_temp[45]);
            h_TrkIso_endcap_data_deno->Add(h_temp[46]);
            h_TrkIso_barrel_data_ctrl->Add(h_temp[47]);
            h_TrkIso_endcap_data_ctrl->Add(h_temp[48]);
            h_ECALiso_barrel_data_nume->Add(h_temp[49]);
            h_ECALiso_endcap_data_nume->Add(h_temp[50]);
            h_ECALiso_barrel_data_deno->Add(h_temp[51]);
            h_ECALiso_endcap_data_deno->Add(h_temp[52]);
            h_ECALiso_barrel_data_ctrl->Add(h_temp[53]);
            h_ECALiso_endcap_data_ctrl->Add(h_temp[54]);
            h_HCALiso_barrel_data_nume->Add(h_temp[55]);
            h_HCALiso_endcap_data_nume->Add(h_temp[56]);
            h_HCALiso_barrel_data_deno->Add(h_temp[57]);
            h_HCALiso_endcap_data_deno->Add(h_temp[58]);
            h_HCALiso_barrel_data_ctrl->Add(h_temp[59]);
            h_HCALiso_endcap_data_ctrl->Add(h_temp[60]);
            h_MET_data->Add(h_temp[61]);
            h_MT_barrel_data_nume->Add(h_temp[62]);
            h_MT_endcap_data_nume->Add(h_temp[63]);
            h_MT_barrel_data_deno->Add(h_temp[64]);
            h_MT_endcap_data_deno->Add(h_temp[65]);
            h_MT_barrel_data_ctrl->Add(h_temp[66]);
            h_MT_endcap_data_ctrl->Add(h_temp[67]);
            h_eta_data->Add(h_temp[68]);
            h_nVTX_data->Add(h_temp[69]);
            h_mass_test_data->Add(h_temp[70]);
            h_HoverE_barrel_template_int_data->Add(h_temp[71]);
            h_HoverE_barrel_jetTemplate_int_data->Add(h_temp[72]);
            h_HoverE_endcap_template_int_data->Add(h_temp[73]);
            h_HoverE_endcap_jetTemplate_int_data->Add(h_temp[74]);
            h_PFiso_Rho_barrel_data_separate->Add(h_temp[75]);
            h_PFiso_Rho_endcap_data_separate->Add(h_temp[76]);
            h_SigmaIEtaIEta_barrel_data_separate->Add(h_temp[77]);
            h_SigmaIEtaIEta_endcap_data_separate->Add(h_temp[78]);
            h_dEtaInSeed_barrel_data_separate->Add(h_temp[79]);
            h_dEtaInSeed_endcap_data_separate->Add(h_temp[80]);
            h_dPhiIn_barrel_data_separate->Add(h_temp[81]);
            h_dPhiIn_endcap_data_separate->Add(h_temp[82]);
            h_HoverE_barrel_data_separate->Add(h_temp[83]);
            h_HoverE_endcap_data_separate->Add(h_temp[84]);
            h_InvEminusInvP_barrel_data_separate->Add(h_temp[85]);
            h_InvEminusInvP_endcap_data_separate->Add(h_temp[86]);
            h_mHits_barrel_data_separate->Add(h_temp[87]);
            h_mHits_endcap_data_separate->Add(h_temp[88]);
            h_passConvVeto_barrel_data_separate->Add(h_temp[89]);
            h_passConvVeto_endcap_data_separate->Add(h_temp[90]);
            h_passMediumID_barrel_data_separate->Add(h_temp[91]);
            h_passMediumID_endcap_data_separate->Add(h_temp[92]);
            h_pT_barrel_data_deno_density->Add(h_temp[9]);
            h_pT_endcap_data_deno_density->Add(h_temp[10]);

            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                TH1D *h_temp_template[4];

                file->GetObject("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10), h_temp_template[0]);
                file->GetObject("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_temp_template[1]);

                removeNegativeBins(h_temp_template[0]);
                removeNegativeBins(h_temp_template[1]);

                h_PFiso_Rho_barrel_data_template[ih]->Add(h_temp_template[0]);
                h_PFiso_Rho_barrel_data_jetTemplate[ih]->Add(h_temp_template[1]);

                if (ih < nPtBin_ele)
                {
                    file->GetObject("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10), h_temp_template[2]);
                    file->GetObject("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_temp_template[3]);

                    removeNegativeBins(h_temp_template[2]);
                    removeNegativeBins(h_temp_template[3]);

                    h_PFiso_Rho_endcap_data_template[ih]->Add(h_temp_template[2]);
                    h_PFiso_Rho_endcap_data_jetTemplate[ih]->Add(h_temp_template[3]);
                }
            }
        }
    }
    for (Int_t i_bin=1; i_bin<=nPtBin_ele; i_bin++)
    {
        h_pT_barrel_data_deno_density->SetBinContent(i_bin, h_pT_barrel_data_deno_density->GetBinContent(i_bin) /
                                                           (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        h_pT_barrel_data_deno_density->SetBinError(i_bin, h_pT_barrel_data_deno_density->GetBinError(i_bin) /
                                                         (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        h_pT_endcap_data_deno_density->SetBinContent(i_bin, h_pT_endcap_data_deno_density->GetBinContent(i_bin) /
                                                           (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
        h_pT_endcap_data_deno_density->SetBinError(i_bin, h_pT_endcap_data_deno_density->GetBinError(i_bin) /
                                                         (analyzer.ptbin_ele[i_bin]-analyzer.ptbin_ele[i_bin-1]));
    }

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
    h_TrkIso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_TrkIso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_TrkIso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_TrkIso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_TrkIso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_TrkIso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_ECALiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_barrel_data_nume->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_endcap_data_nume->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_barrel_data_deno->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_endcap_data_deno->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_barrel_data_ctrl->SetMarkerStyle(kFullDotLarge);
    h_HCALiso_endcap_data_ctrl->SetMarkerStyle(kFullDotLarge);
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
    h_HoverE_barrel_template_int_data->SetMarkerStyle(kFullDotLarge);
    h_HoverE_barrel_jetTemplate_int_data->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_template_int_data->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_jetTemplate_int_data->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_PFiso_Rho_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_SigmaIEtaIEta_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_dEtaInSeed_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_dPhiIn_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_HoverE_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_HoverE_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_InvEminusInvP_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_mHits_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_mHits_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_passConvVeto_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_passConvVeto_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_passMediumID_barrel_data_separate->SetMarkerStyle(kFullDotLarge);
    h_passMediumID_endcap_data_separate->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_data_deno_density->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_deno_density->SetMarkerStyle(kFullDotLarge);

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
    h_TrkIso_barrel_data_nume->SetMarkerColor(kBlack);
    h_TrkIso_endcap_data_nume->SetMarkerColor(kBlack);
    h_TrkIso_barrel_data_deno->SetMarkerColor(kBlack);
    h_TrkIso_endcap_data_deno->SetMarkerColor(kBlack);
    h_TrkIso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_TrkIso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_ECALiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_ECALiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_ECALiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_ECALiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_ECALiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_ECALiso_endcap_data_ctrl->SetMarkerColor(kBlack);
    h_HCALiso_barrel_data_nume->SetMarkerColor(kBlack);
    h_HCALiso_endcap_data_nume->SetMarkerColor(kBlack);
    h_HCALiso_barrel_data_deno->SetMarkerColor(kBlack);
    h_HCALiso_endcap_data_deno->SetMarkerColor(kBlack);
    h_HCALiso_barrel_data_ctrl->SetMarkerColor(kBlack);
    h_HCALiso_endcap_data_ctrl->SetMarkerColor(kBlack);
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
    h_HoverE_barrel_template_int_data->SetMarkerColor(kBlack);
    h_HoverE_barrel_jetTemplate_int_data->SetMarkerColor(kBlack);
    h_HoverE_endcap_template_int_data->SetMarkerColor(kBlack);
    h_HoverE_endcap_jetTemplate_int_data->SetMarkerColor(kBlack);
    h_PFiso_Rho_barrel_data_separate->SetMarkerColor(kBlack);
    h_PFiso_Rho_endcap_data_separate->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_separate->SetMarkerColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_separate->SetMarkerColor(kBlack);
    h_dEtaInSeed_barrel_data_separate->SetMarkerColor(kBlack);
    h_dEtaInSeed_endcap_data_separate->SetMarkerColor(kBlack);
    h_dPhiIn_barrel_data_separate->SetMarkerColor(kBlack);
    h_dPhiIn_endcap_data_separate->SetMarkerColor(kBlack);
    h_HoverE_barrel_data_separate->SetMarkerColor(kBlack);
    h_HoverE_endcap_data_separate->SetMarkerColor(kBlack);
    h_InvEminusInvP_barrel_data_separate->SetMarkerColor(kBlack);
    h_InvEminusInvP_endcap_data_separate->SetMarkerColor(kBlack);
    h_mHits_barrel_data_separate->SetMarkerColor(kBlack);
    h_mHits_endcap_data_separate->SetMarkerColor(kBlack);
    h_passConvVeto_barrel_data_separate->SetMarkerColor(kBlack);
    h_passConvVeto_endcap_data_separate->SetMarkerColor(kBlack);
    h_passMediumID_barrel_data_separate->SetMarkerColor(kBlack);
    h_passMediumID_endcap_data_separate->SetMarkerColor(kBlack);
    h_pT_barrel_data_deno_density->SetMarkerColor(kBlack);
    h_pT_endcap_data_deno_density->SetMarkerColor(kBlack);

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
    h_TrkIso_barrel_data_nume->SetLineColor(kBlack);
    h_TrkIso_endcap_data_nume->SetLineColor(kBlack);
    h_TrkIso_barrel_data_deno->SetLineColor(kBlack);
    h_TrkIso_endcap_data_deno->SetLineColor(kBlack);
    h_TrkIso_barrel_data_ctrl->SetLineColor(kBlack);
    h_TrkIso_endcap_data_ctrl->SetLineColor(kBlack);
    h_ECALiso_barrel_data_nume->SetLineColor(kBlack);
    h_ECALiso_endcap_data_nume->SetLineColor(kBlack);
    h_ECALiso_barrel_data_deno->SetLineColor(kBlack);
    h_ECALiso_endcap_data_deno->SetLineColor(kBlack);
    h_ECALiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_ECALiso_endcap_data_ctrl->SetLineColor(kBlack);
    h_HCALiso_barrel_data_nume->SetLineColor(kBlack);
    h_HCALiso_endcap_data_nume->SetLineColor(kBlack);
    h_HCALiso_barrel_data_deno->SetLineColor(kBlack);
    h_HCALiso_endcap_data_deno->SetLineColor(kBlack);
    h_HCALiso_barrel_data_ctrl->SetLineColor(kBlack);
    h_HCALiso_endcap_data_ctrl->SetLineColor(kBlack);
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
    h_HoverE_barrel_template_int_data->SetLineColor(kBlack);
    h_HoverE_barrel_jetTemplate_int_data->SetLineColor(kBlack);
    h_HoverE_endcap_template_int_data->SetLineColor(kBlack);
    h_HoverE_endcap_jetTemplate_int_data->SetLineColor(kBlack);
    h_PFiso_Rho_barrel_data_separate->SetLineColor(kBlack);
    h_PFiso_Rho_endcap_data_separate->SetLineColor(kBlack);
    h_SigmaIEtaIEta_barrel_data_separate->SetLineColor(kBlack);
    h_SigmaIEtaIEta_endcap_data_separate->SetLineColor(kBlack);
    h_dEtaInSeed_barrel_data_separate->SetLineColor(kBlack);
    h_dEtaInSeed_endcap_data_separate->SetLineColor(kBlack);
    h_dPhiIn_barrel_data_separate->SetLineColor(kBlack);
    h_dPhiIn_endcap_data_separate->SetLineColor(kBlack);
    h_HoverE_barrel_data_separate->SetLineColor(kBlack);
    h_HoverE_endcap_data_separate->SetLineColor(kBlack);
    h_InvEminusInvP_barrel_data_separate->SetLineColor(kBlack);
    h_InvEminusInvP_endcap_data_separate->SetLineColor(kBlack);
    h_mHits_barrel_data_separate->SetLineColor(kBlack);
    h_mHits_endcap_data_separate->SetLineColor(kBlack);
    h_passConvVeto_barrel_data_separate->SetLineColor(kBlack);
    h_passConvVeto_endcap_data_separate->SetLineColor(kBlack);
    h_passMediumID_barrel_data_separate->SetLineColor(kBlack);
    h_passMediumID_endcap_data_separate->SetLineColor(kBlack);
    h_pT_barrel_data_deno_density->SetLineColor(kBlack);
    h_pT_endcap_data_deno_density->SetLineColor(kBlack);

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
    h_TrkIso_barrel_data_nume->SetDirectory(0);
    h_TrkIso_endcap_data_nume->SetDirectory(0);
    h_TrkIso_barrel_data_deno->SetDirectory(0);
    h_TrkIso_endcap_data_deno->SetDirectory(0);
    h_TrkIso_barrel_data_ctrl->SetDirectory(0);
    h_TrkIso_endcap_data_ctrl->SetDirectory(0);
    h_ECALiso_barrel_data_nume->SetDirectory(0);
    h_ECALiso_endcap_data_nume->SetDirectory(0);
    h_ECALiso_barrel_data_deno->SetDirectory(0);
    h_ECALiso_endcap_data_deno->SetDirectory(0);
    h_ECALiso_barrel_data_ctrl->SetDirectory(0);
    h_ECALiso_endcap_data_ctrl->SetDirectory(0);
    h_HCALiso_barrel_data_nume->SetDirectory(0);
    h_HCALiso_endcap_data_nume->SetDirectory(0);
    h_HCALiso_barrel_data_deno->SetDirectory(0);
    h_HCALiso_endcap_data_deno->SetDirectory(0);
    h_HCALiso_barrel_data_ctrl->SetDirectory(0);
    h_HCALiso_endcap_data_ctrl->SetDirectory(0);
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
    h_HoverE_barrel_template_int_data->SetDirectory(0);
    h_HoverE_barrel_jetTemplate_int_data->SetDirectory(0);
    h_HoverE_endcap_template_int_data->SetDirectory(0);
    h_HoverE_endcap_jetTemplate_int_data->SetDirectory(0);
    h_PFiso_Rho_barrel_data_separate->SetDirectory(0);
    h_PFiso_Rho_endcap_data_separate->SetDirectory(0);
    h_SigmaIEtaIEta_barrel_data_separate->SetDirectory(0);
    h_SigmaIEtaIEta_endcap_data_separate->SetDirectory(0);
    h_dEtaInSeed_barrel_data_separate->SetDirectory(0);
    h_dEtaInSeed_endcap_data_separate->SetDirectory(0);
    h_dPhiIn_barrel_data_separate->SetDirectory(0);
    h_dPhiIn_endcap_data_separate->SetDirectory(0);
    h_HoverE_barrel_data_separate->SetDirectory(0);
    h_HoverE_endcap_data_separate->SetDirectory(0);
    h_InvEminusInvP_barrel_data_separate->SetDirectory(0);
    h_InvEminusInvP_endcap_data_separate->SetDirectory(0);
    h_mHits_barrel_data_separate->SetDirectory(0);
    h_mHits_endcap_data_separate->SetDirectory(0);
    h_passConvVeto_barrel_data_separate->SetDirectory(0);
    h_passConvVeto_endcap_data_separate->SetDirectory(0);
    h_passMediumID_barrel_data_separate->SetDirectory(0);
    h_passMediumID_endcap_data_separate->SetDirectory(0);
    h_pT_barrel_data_deno_density->SetDirectory(0);
    h_pT_endcap_data_deno_density->SetDirectory(0);

    for (Int_t ih=0; ih<nPtBin_ele; ih++)
    {
        h_PFiso_Rho_barrel_data_template[ih]->SetMarkerStyle(kFullDotLarge);
        h_PFiso_Rho_barrel_data_jetTemplate[ih]->SetMarkerStyle(kFullDotLarge);

        h_PFiso_Rho_barrel_data_template[ih]->SetMarkerColor(kBlack);
        h_PFiso_Rho_barrel_data_jetTemplate[ih]->SetMarkerColor(kBlack);

        h_PFiso_Rho_barrel_data_template[ih]->SetLineColor(kBlack);
        h_PFiso_Rho_barrel_data_jetTemplate[ih]->SetLineColor(kBlack);

        h_PFiso_Rho_barrel_data_template[ih]->SetDirectory(0);
        h_PFiso_Rho_barrel_data_jetTemplate[ih]->SetDirectory(0);

        if (ih < nPtBin_ele)
        {
            h_PFiso_Rho_endcap_data_template[ih]->SetMarkerStyle(kFullDotLarge);
            h_PFiso_Rho_endcap_data_jetTemplate[ih]->SetMarkerStyle(kFullDotLarge);

            h_PFiso_Rho_endcap_data_template[ih]->SetMarkerColor(kBlack);
            h_PFiso_Rho_endcap_data_jetTemplate[ih]->SetMarkerColor(kBlack);

            h_PFiso_Rho_endcap_data_template[ih]->SetLineColor(kBlack);
            h_PFiso_Rho_endcap_data_jetTemplate[ih]->SetLineColor(kBlack);

            h_PFiso_Rho_endcap_data_template[ih]->SetDirectory(0);
            h_PFiso_Rho_endcap_data_jetTemplate[ih]->SetDirectory(0);
        }
    }

    cout << "Data received" << endl;

//--------------------------------- Ratio Plots --------------------------------------

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
    myRatioPlot_t *RP_TrkIso_barrel_nume = new myRatioPlot_t("RP_TrkIso_barrel_nume", s_TrkIso_barrel_nume, h_TrkIso_barrel_data_nume);
    myRatioPlot_t *RP_TrkIso_endcap_nume = new myRatioPlot_t("RP_TrkIso_endcap_nume", s_TrkIso_endcap_nume, h_TrkIso_endcap_data_nume);
    myRatioPlot_t *RP_TrkIso_barrel_deno = new myRatioPlot_t("RP_TrkIso_barrel_deno", s_TrkIso_barrel_deno, h_TrkIso_barrel_data_deno);
    myRatioPlot_t *RP_TrkIso_endcap_deno = new myRatioPlot_t("RP_TrkIso_endcap_deno", s_TrkIso_endcap_deno, h_TrkIso_endcap_data_deno);
    myRatioPlot_t *RP_TrkIso_barrel_ctrl = new myRatioPlot_t("RP_TrkIso_barrel_ctrl", s_TrkIso_barrel_ctrl, h_TrkIso_barrel_data_ctrl);
    myRatioPlot_t *RP_TrkIso_endcap_ctrl = new myRatioPlot_t("RP_TrkIso_endcap_ctrl", s_TrkIso_endcap_ctrl, h_TrkIso_endcap_data_ctrl);
    myRatioPlot_t *RP_ECALiso_barrel_nume = new myRatioPlot_t("RP_ECALiso_barrel_nume", s_ECALiso_barrel_nume, h_ECALiso_barrel_data_nume);
    myRatioPlot_t *RP_ECALiso_endcap_nume = new myRatioPlot_t("RP_ECALiso_endcap_nume", s_ECALiso_endcap_nume, h_ECALiso_endcap_data_nume);
    myRatioPlot_t *RP_ECALiso_barrel_deno = new myRatioPlot_t("RP_ECALiso_barrel_deno", s_ECALiso_barrel_deno, h_ECALiso_barrel_data_deno);
    myRatioPlot_t *RP_ECALiso_endcap_deno = new myRatioPlot_t("RP_ECALiso_endcap_deno", s_ECALiso_endcap_deno, h_ECALiso_endcap_data_deno);
    myRatioPlot_t *RP_ECALiso_barrel_ctrl = new myRatioPlot_t("RP_ECALiso_barrel_ctrl", s_ECALiso_barrel_ctrl, h_ECALiso_barrel_data_ctrl);
    myRatioPlot_t *RP_ECALiso_endcap_ctrl = new myRatioPlot_t("RP_ECALiso_endcap_ctrl", s_ECALiso_endcap_ctrl, h_ECALiso_endcap_data_ctrl);
    myRatioPlot_t *RP_HCALiso_barrel_nume = new myRatioPlot_t("RP_HCALiso_barrel_nume", s_HCALiso_barrel_nume, h_HCALiso_barrel_data_nume);
    myRatioPlot_t *RP_HCALiso_endcap_nume = new myRatioPlot_t("RP_HCALiso_endcap_nume", s_HCALiso_endcap_nume, h_HCALiso_endcap_data_nume);
    myRatioPlot_t *RP_HCALiso_barrel_deno = new myRatioPlot_t("RP_HCALiso_barrel_deno", s_HCALiso_barrel_deno, h_HCALiso_barrel_data_deno);
    myRatioPlot_t *RP_HCALiso_endcap_deno = new myRatioPlot_t("RP_HCALiso_endcap_deno", s_HCALiso_endcap_deno, h_HCALiso_endcap_data_deno);
    myRatioPlot_t *RP_HCALiso_barrel_ctrl = new myRatioPlot_t("RP_HCALiso_barrel_ctrl", s_HCALiso_barrel_ctrl, h_HCALiso_barrel_data_ctrl);
    myRatioPlot_t *RP_HCALiso_endcap_ctrl = new myRatioPlot_t("RP_HCALiso_endcap_ctrl", s_HCALiso_endcap_ctrl, h_HCALiso_endcap_data_ctrl);
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
    myRatioPlot_t *RP_HoverE_barrel_template_int = new myRatioPlot_t("RP_HoverE_barrel_template_int", s_HoverE_barrel_template_int, h_HoverE_barrel_template_int_data);
    myRatioPlot_t *RP_HoverE_barrel_jetTemplate_int = new myRatioPlot_t("RP_HoverE_barrel_jetTemplate_int", s_HoverE_barrel_jetTemplate_int, h_HoverE_barrel_jetTemplate_int_data);
    myRatioPlot_t *RP_HoverE_endcap_template_int = new myRatioPlot_t("RP_HoverE_endcap_template_int", s_HoverE_endcap_template_int, h_HoverE_endcap_template_int_data);
    myRatioPlot_t *RP_HoverE_endcap_jetTemplate_int = new myRatioPlot_t("RP_HoverE_endcap_jetTemplate_int", s_HoverE_endcap_jetTemplate_int, h_HoverE_endcap_jetTemplate_int_data);   
    myRatioPlot_t *RP_PFiso_Rho_barrel_separate = new myRatioPlot_t("RP_PFiso_Rho_barrel_separate", s_PFiso_Rho_barrel_separate, h_PFiso_Rho_barrel_data_separate);
    myRatioPlot_t *RP_PFiso_Rho_endcap_separate = new myRatioPlot_t("RP_PFiso_Rho_endcap_separate", s_PFiso_Rho_endcap_separate, h_PFiso_Rho_endcap_data_separate);
    myRatioPlot_t *RP_SigmaIEtaIEta_barrel_separate = new myRatioPlot_t("RP_SigmaIEtaIEta_barrel_separate", s_SigmaIEtaIEta_barrel_separate, h_SigmaIEtaIEta_barrel_data_separate);
    myRatioPlot_t *RP_SigmaIEtaIEta_endcap_separate = new myRatioPlot_t("RP_SigmaIEtaIEta_endcap_separate", s_SigmaIEtaIEta_endcap_separate, h_SigmaIEtaIEta_endcap_data_separate);
    myRatioPlot_t *RP_dEtaInSeed_barrel_separate = new myRatioPlot_t("RP_dEtaInSeed_barrel_separate", s_dEtaInSeed_barrel_separate, h_dEtaInSeed_barrel_data_separate);
    myRatioPlot_t *RP_dEtaInSeed_endcap_separate = new myRatioPlot_t("RP_dEtaInSeed_endcap_separate", s_dEtaInSeed_endcap_separate, h_dEtaInSeed_endcap_data_separate);
    myRatioPlot_t *RP_dPhiIn_barrel_separate = new myRatioPlot_t("RP_dPhiIn_barrel_separate", s_dPhiIn_barrel_separate, h_dPhiIn_barrel_data_separate);
    myRatioPlot_t *RP_dPhiIn_endcap_separate = new myRatioPlot_t("RP_dPhiIn_endcap_separate", s_dPhiIn_endcap_separate, h_dPhiIn_endcap_data_separate);
    myRatioPlot_t *RP_HoverE_barrel_separate = new myRatioPlot_t("RP_HoverE_barrel_separate", s_HoverE_barrel_separate, h_HoverE_barrel_data_separate);
    myRatioPlot_t *RP_HoverE_endcap_separate = new myRatioPlot_t("RP_HoverE_endcap_separate", s_HoverE_endcap_separate, h_HoverE_endcap_data_separate);
    myRatioPlot_t *RP_InvEminusInvP_barrel_separate = new myRatioPlot_t("RP_InvEminusInvP_barrel_separate", s_InvEminusInvP_barrel_separate, h_InvEminusInvP_barrel_data_separate);
    myRatioPlot_t *RP_InvEminusInvP_endcap_separate = new myRatioPlot_t("RP_InvEminusInvP_endcap_separate", s_InvEminusInvP_endcap_separate, h_InvEminusInvP_endcap_data_separate);
    myRatioPlot_t *RP_mHits_barrel_separate = new myRatioPlot_t("RP_mHits_barrel_separate", s_mHits_barrel_separate, h_mHits_barrel_data_separate);
    myRatioPlot_t *RP_mHits_endcap_separate = new myRatioPlot_t("RP_mHits_endcap_separate", s_mHits_endcap_separate, h_mHits_endcap_data_separate);
    myRatioPlot_t *RP_passConvVeto_barrel_separate = new myRatioPlot_t("RP_passConvVeto_barrel_separate", s_passConvVeto_barrel_separate, h_passConvVeto_barrel_data_separate);
    myRatioPlot_t *RP_passConvVeto_endcap_separate = new myRatioPlot_t("RP_passConvVeto_endcap_separate", s_passConvVeto_endcap_separate, h_passConvVeto_endcap_data_separate);
    myRatioPlot_t *RP_passMediumID_barrel_separate = new myRatioPlot_t("RP_passMediumID_barrel_separate", s_passMediumID_barrel_separate, h_passMediumID_barrel_data_separate);
    myRatioPlot_t *RP_passMediumID_endcap_separate = new myRatioPlot_t("RP_passMediumID_endcap_separate", s_passMediumID_endcap_separate, h_passMediumID_endcap_data_separate);
    myRatioPlot_t *RP_pT_barrel_deno_density = new myRatioPlot_t("RP_pT_barrel_deno_density", s_pT_barrel_deno_density, h_pT_barrel_data_deno_density);
    myRatioPlot_t *RP_pT_endcap_deno_density = new myRatioPlot_t("RP_pT_endcap_deno_density", s_pT_endcap_deno_density, h_pT_endcap_data_deno_density);

    RP_PFiso_Rho_barrel_nume->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{nume})", 0, 0.2);
    RP_PFiso_Rho_endcap_nume->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{nume})", 0, 0.2);
    RP_PFiso_Rho_barrel_deno->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{deno})", 0, 5);
    RP_PFiso_Rho_endcap_deno->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{deno})", 0, 5);
    RP_PFiso_Rho_barrel_ctrl->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{barrel}}^{control})", 0, 5);
    RP_PFiso_Rho_endcap_ctrl->SetPlots("relPFiso_Rho (e_{#lower[-0.4]{endcap}}^{control})", 0, 5);
    RP_pT_barrel_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{nume}) [GeV/c]", 25, 800);
    RP_pT_endcap_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{nume}) [GeV/c]", 25, 800);
    RP_pT_barrel_deno->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{deno}) [GeV/c]", 25, 800);
    RP_pT_endcap_deno->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{deno}) [GeV/c]", 25, 800);
    RP_pT_barrel_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{control}) [GeV/c]", 25, 800);
    RP_pT_endcap_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{control}) [GeV/c]", 25, 800);
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
    RP_HoverE_barrel_deno->SetPlots("H/E (e_{#lower[-0.4]{barrel}}^{deno})", 0, 1);
    RP_HoverE_endcap_deno->SetPlots("H/E (e_{#lower[-0.4]{endcap}}^{deno})", 0, 1);
    RP_HoverE_barrel_ctrl->SetPlots("H/E (e_{#lower[-0.4]{barrel}}^{control})", 0, 1);
    RP_HoverE_endcap_ctrl->SetPlots("H/E (e_{#lower[-0.4]{endcap}}^{control})", 0, 1);
    RP_InvEminusInvP_barrel_nume->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{nume})", -0.5, 0.5);
    RP_InvEminusInvP_endcap_nume->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{nume})", -0.5, 0.5);
    RP_InvEminusInvP_barrel_deno->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{deno})", -6, 6);
    RP_InvEminusInvP_endcap_deno->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{deno})", -6, 6);
    RP_InvEminusInvP_barrel_ctrl->SetPlots("1/E-1/P (e_{#lower[-0.4]{barrel}}^{control})", -6, 6);
    RP_InvEminusInvP_endcap_ctrl->SetPlots("1/E-1/P (e_{#lower[-0.4]{endcap}}^{control})", -6, 6);
    RP_TrkIso_barrel_nume->SetPlots("relTrkIso (e_{#lower[-0.4]{barrel}}^{nume})", 0, 2);
    RP_TrkIso_endcap_nume->SetPlots("relTrkIso (e_{#lower[-0.4]{endcap}}^{nume})", 0, 2);
    RP_TrkIso_barrel_deno->SetPlots("relTrkIso (e_{#lower[-0.4]{barrel}}^{deno})", 0, 2);
    RP_TrkIso_endcap_deno->SetPlots("relTrkIso (e_{#lower[-0.4]{endcap}}^{deno})", 0, 2);
    RP_TrkIso_barrel_ctrl->SetPlots("relTrkIso (e_{#lower[-0.4]{barrel}}^{control})", 0, 2);
    RP_TrkIso_endcap_ctrl->SetPlots("relTrkIso (e_{#lower[-0.4]{endcap}}^{control})", 0, 2);
    RP_ECALiso_barrel_nume->SetPlots("relECALiso (e_{#lower[-0.4]{barrel}}^{nume})", 0, 2);
    RP_ECALiso_endcap_nume->SetPlots("relECALiso (e_{#lower[-0.4]{endcap}}^{nume})", 0, 2);
    RP_ECALiso_barrel_deno->SetPlots("relECALiso (e_{#lower[-0.4]{barrel}}^{deno})", 0, 2);
    RP_ECALiso_endcap_deno->SetPlots("relECALiso (e_{#lower[-0.4]{endcap}}^{deno})", 0, 2);
    RP_ECALiso_barrel_ctrl->SetPlots("relECALiso (e_{#lower[-0.4]{barrel}}^{control})", 0, 2);
    RP_ECALiso_endcap_ctrl->SetPlots("relECALiso (e_{#lower[-0.4]{endcap}}^{control})", 0, 2);
    RP_HCALiso_barrel_nume->SetPlots("relHCALiso (e_{#lower[-0.4]{barrel}}^{nume})", 0, 2);
    RP_HCALiso_endcap_nume->SetPlots("relHCALiso (e_{#lower[-0.4]{endcap}}^{nume})", 0, 2);
    RP_HCALiso_barrel_deno->SetPlots("relHCALiso (e_{#lower[-0.4]{barrel}}^{deno})", 0, 2);
    RP_HCALiso_endcap_deno->SetPlots("relHCALiso (e_{#lower[-0.4]{endcap}}^{deno})", 0, 2);
    RP_HCALiso_barrel_ctrl->SetPlots("relHCALiso (e_{#lower[-0.4]{barrel}}^{control})", 0, 2);
    RP_HCALiso_endcap_ctrl->SetPlots("relHCALiso (e_{#lower[-0.4]{endcap}}^{control})", 0, 2);
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
    RP_HoverE_barrel_template_int->SetPlots("H/E (barrel)", 0, 0.2);
    RP_HoverE_barrel_jetTemplate_int->SetPlots("H/E (barrel, jet selection)", 0, 0.2);
    RP_HoverE_endcap_template_int->SetPlots("H/E (endcap)", 0, 0.5);
    RP_HoverE_endcap_jetTemplate_int->SetPlots("H/E (endcap, jet selection)", 0, 0.5);
    RP_PFiso_Rho_barrel_separate->SetPlots("relPFiso_Rho", 0, 5);
    RP_PFiso_Rho_endcap_separate->SetPlots("relPFiso_Rho", 0, 5);
    RP_SigmaIEtaIEta_barrel_separate->SetPlots("#sigma_{i#etai#eta}", 0, 0.05);
    RP_SigmaIEtaIEta_endcap_separate->SetPlots("#sigma_{i#etai#eta}", 0, 0.1);
    RP_dEtaInSeed_barrel_separate->SetPlots("dEtaInSeed", -0.1, 0.1);
    RP_dEtaInSeed_endcap_separate->SetPlots("dEtaInSeed", -1, 1);
    RP_dPhiIn_barrel_separate->SetPlots("dPhiIn", -1, 1);
    RP_dPhiIn_endcap_separate->SetPlots("dPhiIn", -1, 1);
    RP_HoverE_barrel_separate->SetPlots("H/E", 0, 1);
    RP_HoverE_endcap_separate->SetPlots("H/E", 0, 1);
    RP_InvEminusInvP_barrel_separate->SetPlots("1/E-1/P", -6, 6);
    RP_InvEminusInvP_endcap_separate->SetPlots("1/E-1/P", -6, 6);
    RP_mHits_barrel_separate->SetPlots("# missing hits", 0-0.5, 3-0.5);
    RP_mHits_endcap_separate->SetPlots("# missing hits", 0-0.5, 3-0.5);
    RP_passConvVeto_barrel_separate->SetPlots("Pass Conversion Veto", 0-0.5, 2-0.5);
    RP_passConvVeto_endcap_separate->SetPlots("Pass Conversion Veto", 0-0.5, 2-0.5);
    RP_passMediumID_barrel_separate->SetPlots("Pass MediumID", 0-0.5, 2-0.5);
    RP_passMediumID_endcap_separate->SetPlots("Pass MediumID", 0-0.5, 2-0.5);
    RP_pT_barrel_deno_density->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{deno}) [GeV/c]", 25, 800);
    RP_pT_endcap_deno_density->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{deno}) [GeV/c]", 25, 800);


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
//    legend->AddEntry(h_eta_MC[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->AddEntry(h_eta_MC[_GJets_20to100], "#gamma+Jets", "f");
    legend->SetNColumns(2);

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
    RP_TrkIso_barrel_nume->ImportLegend(legend);
    RP_TrkIso_endcap_nume->ImportLegend(legend);
    RP_TrkIso_barrel_deno->ImportLegend(legend);
    RP_TrkIso_endcap_deno->ImportLegend(legend);
    RP_TrkIso_barrel_ctrl->ImportLegend(legend);
    RP_TrkIso_endcap_ctrl->ImportLegend(legend);
    RP_ECALiso_barrel_nume->ImportLegend(legend);
    RP_ECALiso_endcap_nume->ImportLegend(legend);
    RP_ECALiso_barrel_deno->ImportLegend(legend);
    RP_ECALiso_endcap_deno->ImportLegend(legend);
    RP_ECALiso_barrel_ctrl->ImportLegend(legend);
    RP_ECALiso_endcap_ctrl->ImportLegend(legend);
    RP_HCALiso_barrel_nume->ImportLegend(legend);
    RP_HCALiso_endcap_nume->ImportLegend(legend);
    RP_HCALiso_barrel_deno->ImportLegend(legend);
    RP_HCALiso_endcap_deno->ImportLegend(legend);
    RP_HCALiso_barrel_ctrl->ImportLegend(legend);
    RP_HCALiso_endcap_ctrl->ImportLegend(legend);
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
    RP_HoverE_barrel_template_int->ImportLegend(legend);
    RP_HoverE_barrel_jetTemplate_int->ImportLegend(legend);
    RP_HoverE_endcap_template_int->ImportLegend(legend);
    RP_HoverE_endcap_jetTemplate_int->ImportLegend(legend);
    RP_PFiso_Rho_barrel_separate->ImportLegend(legend);
    RP_PFiso_Rho_endcap_separate->ImportLegend(legend);
    RP_SigmaIEtaIEta_barrel_separate->ImportLegend(legend);
    RP_SigmaIEtaIEta_endcap_separate->ImportLegend(legend);
    RP_dEtaInSeed_barrel_separate->ImportLegend(legend);
    RP_dEtaInSeed_endcap_separate->ImportLegend(legend);
    RP_dPhiIn_barrel_separate->ImportLegend(legend);
    RP_dPhiIn_endcap_separate->ImportLegend(legend);
    RP_HoverE_barrel_separate->ImportLegend(legend);
    RP_HoverE_endcap_separate->ImportLegend(legend);
    RP_InvEminusInvP_barrel_separate->ImportLegend(legend);
    RP_InvEminusInvP_endcap_separate->ImportLegend(legend);
    RP_mHits_barrel_separate->ImportLegend(legend);
    RP_mHits_endcap_separate->ImportLegend(legend);
    RP_passConvVeto_barrel_separate->ImportLegend(legend);
    RP_passConvVeto_endcap_separate->ImportLegend(legend);
    RP_passMediumID_barrel_separate->ImportLegend(legend);
    RP_passMediumID_endcap_separate->ImportLegend(legend);
    RP_pT_barrel_deno_density->ImportLegend(legend);
    RP_pT_endcap_deno_density->ImportLegend(legend);

//    RP_PFiso_Rho_barrel_nume->Draw(1, 1e10, 0);
//    RP_PFiso_Rho_endcap_nume->Draw(1, 1e10, 0);
//    RP_PFiso_Rho_barrel_deno->Draw(1, 1e10, 0);
//    RP_PFiso_Rho_endcap_deno->Draw(1, 1e10, 0);
    RP_PFiso_Rho_barrel_ctrl->Draw(1, 1e10, 0);
    RP_PFiso_Rho_endcap_ctrl->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_barrel_nume->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_endcap_nume->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_barrel_deno->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_endcap_deno->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_SigmaIEtaIEta_endcap_ctrl->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_barrel_nume->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_endcap_nume->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_barrel_deno->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_endcap_deno->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_dEtaInSeed_endcap_ctrl->Draw(1, 1e10, 0);
//    RP_dPhiIn_barrel_nume->Draw(1, 1e10, 0);
//    RP_dPhiIn_endcap_nume->Draw(1, 1e10, 0);
//    RP_dPhiIn_barrel_deno->Draw(1, 1e10, 0);
//    RP_dPhiIn_endcap_deno->Draw(1, 1e10, 0);
//    RP_dPhiIn_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_dPhiIn_endcap_ctrl->Draw(1, 1e10, 0);
//    RP_HoverE_barrel_nume->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_nume->Draw(1, 1e10, 0);
//    RP_HoverE_barrel_deno->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_deno->Draw(1, 1e10, 0);
//    RP_HoverE_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_ctrl->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_barrel_nume->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_endcap_nume->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_barrel_deno->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_endcap_deno->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_InvEminusInvP_endcap_ctrl->Draw(1, 1e10, 0);
    RP_pT_barrel_nume->Draw(1, 1e10, 0);
    RP_pT_endcap_nume->Draw(1, 1e10, 0);
    RP_pT_barrel_deno->Draw(1, 1e10, 0);
    RP_pT_endcap_deno->Draw(1, 1e10, 0);
    RP_pT_barrel_ctrl->Draw(1, 1e10, 0);
    RP_pT_endcap_ctrl->Draw(1, 1e10, 0);
    RP_TrkIso_barrel_nume->Draw(1, 1e10, 0);
    RP_TrkIso_endcap_nume->Draw(1, 1e10, 0);
    RP_TrkIso_barrel_deno->Draw(1, 1e10, 0);
    RP_TrkIso_endcap_deno->Draw(1, 1e10, 0);
    RP_TrkIso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_TrkIso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_ECALiso_barrel_nume->Draw(1, 1e10, 0);
    RP_ECALiso_endcap_nume->Draw(1, 1e10, 0);
    RP_ECALiso_barrel_deno->Draw(1, 1e10, 0);
    RP_ECALiso_endcap_deno->Draw(1, 1e10, 0);
    RP_ECALiso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_ECALiso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_HCALiso_barrel_nume->Draw(1, 1e10, 0);
    RP_HCALiso_endcap_nume->Draw(1, 1e10, 0);
    RP_HCALiso_barrel_deno->Draw(1, 1e10, 0);
    RP_HCALiso_endcap_deno->Draw(1, 1e10, 0);
    RP_HCALiso_barrel_ctrl->Draw(1, 1e10, 0);
    RP_HCALiso_endcap_ctrl->Draw(1, 1e10, 0);
    RP_MET->Draw(1, 1e10, 0);
//    RP_MT_barrel_nume->Draw(1, 1e10, 0);
//    RP_MT_endcap_nume->Draw(1, 1e10, 0);
//    RP_MT_barrel_deno->Draw(1, 1e10, 0);
//    RP_MT_endcap_deno->Draw(1, 1e10, 0);
//    RP_MT_barrel_ctrl->Draw(1, 1e10, 0);
//    RP_MT_endcap_ctrl->Draw(1, 1e10, 0);
    RP_eta->Draw(1, 1e12, 0);
    RP_nVTX->Draw(1, 1e10, 0);
    RP_mass_test->Draw(1, 1e9, 1);
//    RP_HoverE_barrel_template_int->Draw(1, 1e10, 0);
//    RP_HoverE_barrel_jetTemplate_int->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_template_int->Draw(1, 1e10, 0);
//    RP_HoverE_endcap_jetTemplate_int->Draw(1, 1e10, 0);
    RP_PFiso_Rho_barrel_separate->Draw(1, 1e10, 0);
    RP_PFiso_Rho_endcap_separate->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_barrel_separate->Draw(1, 1e10, 0);
    RP_SigmaIEtaIEta_endcap_separate->Draw(1, 1e10, 0);
    RP_dEtaInSeed_barrel_separate->Draw(1, 1e10, 0);
    RP_dEtaInSeed_endcap_separate->Draw(1, 1e10, 0);
    RP_dPhiIn_barrel_separate->Draw(1, 1e10, 0);
    RP_dPhiIn_endcap_separate->Draw(1, 1e10, 0);
    RP_HoverE_barrel_separate->Draw(1, 1e10, 0);
    RP_HoverE_endcap_separate->Draw(1, 1e10, 0);
    RP_InvEminusInvP_barrel_separate->Draw(1, 1e10, 0);
    RP_InvEminusInvP_endcap_separate->Draw(1, 1e10, 0);
//    RP_mHits_barrel_separate->Draw(1, 1e10, 0);
//    RP_mHits_endcap_separate->Draw(1, 1e10, 0);
//    RP_passConvVeto_barrel_separate->Draw(1, 1e10, 0);
//    RP_passConvVeto_endcap_separate->Draw(1, 1e10, 0);
    RP_passMediumID_barrel_separate->Draw(1, 1e10, 0);
    RP_passMediumID_endcap_separate->Draw(1, 1e10, 0);
    RP_pT_barrel_deno_density->Draw(1, 1e10, 0, "HIST", "Number of events / GeV");
    RP_pT_endcap_deno_density->Draw(1, 1e10, 0, "HIST", "Number of events / GeV");

    cout << "Main histos drawn." << endl;

    // Check prompt/fake electron separation by different variables
    // Ammount of QCD
    Double_t n_pass_PFiso_barrel, n_pass_PFiso_endcap, n_fail_PFiso_barrel, n_fail_PFiso_endcap;
    Double_t n_pass_sigma_barrel, n_pass_sigma_endcap, n_fail_sigma_barrel, n_fail_sigma_endcap;
    Double_t n_pass_dEta_barrel,  n_pass_dEta_endcap,  n_fail_dEta_barrel,  n_fail_dEta_endcap;
    Double_t n_pass_dPhi_barrel,  n_pass_dPhi_endcap,  n_fail_dPhi_barrel,  n_fail_dPhi_endcap;
    Double_t n_pass_HovE_barrel,  n_pass_HovE_endcap,  n_fail_HovE_barrel,  n_fail_HovE_endcap;
    Double_t n_pass_iEmiP_barrel, n_pass_iEmiP_endcap, n_fail_iEmiP_barrel, n_fail_iEmiP_endcap;
    Double_t n_pass_pMID_barrel,  n_pass_pMID_endcap,  n_fail_pMID_barrel,  n_fail_pMID_endcap;
    // Ammount of data
    Double_t nd_pass_PFiso_barrel, nd_pass_PFiso_endcap, nd_fail_PFiso_barrel, nd_fail_PFiso_endcap;
    Double_t nd_pass_sigma_barrel, nd_pass_sigma_endcap, nd_fail_sigma_barrel, nd_fail_sigma_endcap;
    Double_t nd_pass_dEta_barrel,  nd_pass_dEta_endcap,  nd_fail_dEta_barrel,  nd_fail_dEta_endcap;
    Double_t nd_pass_dPhi_barrel,  nd_pass_dPhi_endcap,  nd_fail_dPhi_barrel,  nd_fail_dPhi_endcap;
    Double_t nd_pass_HovE_barrel,  nd_pass_HovE_endcap,  nd_fail_HovE_barrel,  nd_fail_HovE_endcap;
    Double_t nd_pass_iEmiP_barrel, nd_pass_iEmiP_endcap, nd_fail_iEmiP_barrel, nd_fail_iEmiP_endcap;
    Double_t nd_pass_pMID_barrel,  nd_pass_pMID_endcap,  nd_fail_pMID_barrel,  nd_fail_pMID_endcap;
    // Percentage of QCD
    Double_t p_pass_PFiso_barrel, p_pass_PFiso_endcap, p_fail_PFiso_barrel, p_fail_PFiso_endcap;
    Double_t p_pass_sigma_barrel, p_pass_sigma_endcap, p_fail_sigma_barrel, p_fail_sigma_endcap;
    Double_t p_pass_dEta_barrel,  p_pass_dEta_endcap,  p_fail_dEta_barrel,  p_fail_dEta_endcap;
    Double_t p_pass_dPhi_barrel,  p_pass_dPhi_endcap,  p_fail_dPhi_barrel,  p_fail_dPhi_endcap;
    Double_t p_pass_HovE_barrel,  p_pass_HovE_endcap,  p_fail_HovE_barrel,  p_fail_HovE_endcap;
    Double_t p_pass_iEmiP_barrel, p_pass_iEmiP_endcap, p_fail_iEmiP_barrel, p_fail_iEmiP_endcap;
    Double_t p_pass_pMID_barrel,  p_pass_pMID_endcap,  p_fail_pMID_barrel,  p_fail_pMID_endcap;

    nd_pass_PFiso_barrel = h_PFiso_Rho_barrel_data_separate->Integral(1,7) ;
    nd_pass_PFiso_endcap = h_PFiso_Rho_endcap_data_separate->Integral(1,8) ;
    nd_fail_PFiso_barrel = h_PFiso_Rho_barrel_data_separate->Integral(8,80);
    nd_fail_PFiso_endcap = h_PFiso_Rho_endcap_data_separate->Integral(9,80);

    nd_pass_sigma_barrel = h_SigmaIEtaIEta_barrel_data_separate->Integral(1,10) ;
    nd_pass_sigma_endcap = h_SigmaIEtaIEta_endcap_data_separate->Integral(1,30) ;
    nd_fail_sigma_barrel = h_SigmaIEtaIEta_barrel_data_separate->Integral(11,13);
    nd_fail_sigma_endcap = h_SigmaIEtaIEta_endcap_data_separate->Integral(31,35);

    nd_pass_dEta_barrel = h_dEtaInSeed_barrel_data_separate->Integral(1,3) ;
    nd_pass_dEta_endcap = h_dEtaInSeed_endcap_data_separate->Integral(1,6) ;
    nd_fail_dEta_barrel = h_dEtaInSeed_barrel_data_separate->Integral(4,10);
    nd_fail_dEta_endcap = h_dEtaInSeed_endcap_data_separate->Integral(7,30);

    nd_pass_dPhi_barrel = h_dPhiIn_barrel_data_separate->Integral(1,10);
    nd_pass_dPhi_endcap = h_dPhiIn_endcap_data_separate->Integral(1,5) ;
    nd_fail_dPhi_barrel = 0;
    nd_fail_dPhi_endcap = h_dPhiIn_endcap_data_separate->Integral(6,30);

    nd_pass_HovE_barrel = h_HoverE_barrel_data_separate->Integral(1,15);
    nd_pass_HovE_endcap = h_HoverE_endcap_data_separate->Integral(1,9) ;
    nd_fail_HovE_barrel = 0;
    nd_fail_HovE_endcap = h_HoverE_endcap_data_separate->Integral(10,15);

    nd_pass_iEmiP_barrel = h_InvEminusInvP_barrel_data_separate->Integral(1,1) ;
    nd_pass_iEmiP_endcap = h_InvEminusInvP_endcap_data_separate->Integral(1,1) ;
    nd_fail_iEmiP_barrel = h_InvEminusInvP_barrel_data_separate->Integral(2,10);
    nd_fail_iEmiP_endcap = h_InvEminusInvP_endcap_data_separate->Integral(2,10);

    nd_pass_pMID_barrel = h_passMediumID_barrel_data_separate->Integral(2,2);
    nd_pass_pMID_endcap = h_passMediumID_endcap_data_separate->Integral(2,2);
    nd_fail_pMID_barrel = h_passMediumID_barrel_data_separate->Integral(1,1);
    nd_fail_pMID_endcap = h_passMediumID_endcap_data_separate->Integral(1,1);


    n_pass_PFiso_barrel = h_PFiso_Rho_barrel_data_separate->Integral(1,7)  - ((TH1D*)s_PFiso_Rho_barrel_separate->GetStack()->Last())->Integral(1,7) ;
    n_pass_PFiso_endcap = h_PFiso_Rho_endcap_data_separate->Integral(1,8)  - ((TH1D*)s_PFiso_Rho_endcap_separate->GetStack()->Last())->Integral(1,8) ;
    n_fail_PFiso_barrel = h_PFiso_Rho_barrel_data_separate->Integral(8,80) - ((TH1D*)s_PFiso_Rho_barrel_separate->GetStack()->Last())->Integral(8,80);
    n_fail_PFiso_endcap = h_PFiso_Rho_endcap_data_separate->Integral(9,80) - ((TH1D*)s_PFiso_Rho_endcap_separate->GetStack()->Last())->Integral(9,80);

    n_pass_sigma_barrel = h_SigmaIEtaIEta_barrel_data_separate->Integral(1,10)  - ((TH1D*)s_SigmaIEtaIEta_barrel_separate->GetStack()->Last())->Integral(1,10) ;
    n_pass_sigma_endcap = h_SigmaIEtaIEta_endcap_data_separate->Integral(1,30)  - ((TH1D*)s_SigmaIEtaIEta_endcap_separate->GetStack()->Last())->Integral(1,30) ;
    n_fail_sigma_barrel = h_SigmaIEtaIEta_barrel_data_separate->Integral(11,13) - ((TH1D*)s_SigmaIEtaIEta_barrel_separate->GetStack()->Last())->Integral(11,13);
    n_fail_sigma_endcap = h_SigmaIEtaIEta_endcap_data_separate->Integral(31,35) - ((TH1D*)s_SigmaIEtaIEta_endcap_separate->GetStack()->Last())->Integral(31,35);

    n_pass_dEta_barrel = h_dEtaInSeed_barrel_data_separate->Integral(1,3)  - ((TH1D*)s_dEtaInSeed_barrel_separate->GetStack()->Last())->Integral(1,3) ;
    n_pass_dEta_endcap = h_dEtaInSeed_endcap_data_separate->Integral(1,6)  - ((TH1D*)s_dEtaInSeed_endcap_separate->GetStack()->Last())->Integral(1,6) ;
    n_fail_dEta_barrel = h_dEtaInSeed_barrel_data_separate->Integral(4,10) - ((TH1D*)s_dEtaInSeed_barrel_separate->GetStack()->Last())->Integral(4,10);
    n_fail_dEta_endcap = h_dEtaInSeed_endcap_data_separate->Integral(7,30) - ((TH1D*)s_dEtaInSeed_endcap_separate->GetStack()->Last())->Integral(7,30);

    n_pass_dPhi_barrel = h_dPhiIn_barrel_data_separate->Integral(1,10)  - ((TH1D*)s_dPhiIn_barrel_separate->GetStack()->Last())->Integral(1,10);
    n_pass_dPhi_endcap = h_dPhiIn_endcap_data_separate->Integral(1,5)   - ((TH1D*)s_dPhiIn_endcap_separate->GetStack()->Last())->Integral(1,5);
    n_fail_dPhi_barrel = 0;
    n_fail_dPhi_endcap = h_dPhiIn_endcap_data_separate->Integral(6,30)  - ((TH1D*)s_dPhiIn_endcap_separate->GetStack()->Last())->Integral(6,30);

    n_pass_HovE_barrel = h_HoverE_barrel_data_separate->Integral(1,15)  - ((TH1D*)s_HoverE_barrel_separate->GetStack()->Last())->Integral(1,15);
    n_pass_HovE_endcap = h_HoverE_endcap_data_separate->Integral(1,9)   - ((TH1D*)s_HoverE_endcap_separate->GetStack()->Last())->Integral(1,9);
    n_fail_HovE_barrel = 0;
    n_fail_HovE_endcap = h_HoverE_endcap_data_separate->Integral(10,15) - ((TH1D*)s_HoverE_endcap_separate->GetStack()->Last())->Integral(10,15);

    n_pass_iEmiP_barrel = h_InvEminusInvP_barrel_data_separate->Integral(1,1)  - ((TH1D*)s_InvEminusInvP_barrel_separate->GetStack()->Last())->Integral(1,1) ;
    n_pass_iEmiP_endcap = h_InvEminusInvP_endcap_data_separate->Integral(1,1)  - ((TH1D*)s_InvEminusInvP_endcap_separate->GetStack()->Last())->Integral(1,1) ;
    n_fail_iEmiP_barrel = h_InvEminusInvP_barrel_data_separate->Integral(2,10) - ((TH1D*)s_InvEminusInvP_barrel_separate->GetStack()->Last())->Integral(2,10);
    n_fail_iEmiP_endcap = h_InvEminusInvP_endcap_data_separate->Integral(2,10) - ((TH1D*)s_InvEminusInvP_endcap_separate->GetStack()->Last())->Integral(2,10);

    n_pass_pMID_barrel = h_passMediumID_barrel_data_separate->Integral(2,2) - ((TH1D*)s_passMediumID_barrel_separate->GetStack()->Last())->Integral(1,1);
    n_pass_pMID_endcap = h_passMediumID_endcap_data_separate->Integral(2,2) - ((TH1D*)s_passMediumID_endcap_separate->GetStack()->Last())->Integral(1,1);
    n_fail_pMID_barrel = h_passMediumID_barrel_data_separate->Integral(1,1) - ((TH1D*)s_passMediumID_barrel_separate->GetStack()->Last())->Integral(2,2);
    n_fail_pMID_endcap = h_passMediumID_endcap_data_separate->Integral(1,1) - ((TH1D*)s_passMediumID_endcap_separate->GetStack()->Last())->Integral(2,2);

    p_pass_PFiso_barrel = n_pass_PFiso_barrel / nd_pass_PFiso_barrel;
    p_pass_PFiso_endcap = n_pass_PFiso_endcap / nd_pass_PFiso_endcap;
    p_fail_PFiso_barrel = n_fail_PFiso_barrel / nd_fail_PFiso_barrel;
    p_fail_PFiso_endcap = n_fail_PFiso_endcap / nd_fail_PFiso_endcap;

    p_pass_sigma_barrel = n_pass_sigma_barrel / nd_pass_sigma_barrel;
    p_pass_sigma_endcap = n_pass_sigma_endcap / nd_pass_sigma_endcap;
    p_fail_sigma_barrel = n_fail_sigma_barrel / nd_fail_sigma_barrel;
    p_fail_sigma_endcap = n_fail_sigma_endcap / nd_fail_sigma_endcap;

    p_pass_dEta_barrel  = n_pass_dEta_barrel  / nd_pass_dEta_barrel ;
    p_pass_dEta_endcap  = n_pass_dEta_endcap  / nd_pass_dEta_endcap ;
    p_fail_dEta_barrel  = n_fail_dEta_barrel  / nd_fail_dEta_barrel ;
    p_fail_dEta_endcap  = n_fail_dEta_endcap  / nd_fail_dEta_endcap ;

    p_pass_dPhi_barrel  = n_pass_dPhi_barrel  / nd_pass_dPhi_barrel ;
    p_pass_dPhi_endcap  = n_pass_dPhi_endcap  / nd_pass_dPhi_endcap ;
    p_fail_dPhi_barrel  = n_fail_dPhi_barrel  / nd_fail_dPhi_barrel ;
    p_fail_dPhi_endcap  = n_fail_dPhi_endcap  / nd_fail_dPhi_endcap ;

    p_pass_HovE_barrel  = n_pass_HovE_barrel  / nd_pass_HovE_barrel ;
    p_pass_HovE_endcap  = n_pass_HovE_endcap  / nd_pass_HovE_endcap ;
    p_fail_HovE_barrel  = n_fail_HovE_barrel  / nd_fail_HovE_barrel ;
    p_fail_HovE_endcap  = n_fail_HovE_endcap  / nd_fail_HovE_endcap ;

    p_pass_iEmiP_barrel = n_pass_iEmiP_barrel / nd_pass_iEmiP_barrel;
    p_pass_iEmiP_endcap = n_pass_iEmiP_endcap / nd_pass_iEmiP_endcap;
    p_fail_iEmiP_barrel = n_fail_iEmiP_barrel / nd_fail_iEmiP_barrel;
    p_fail_iEmiP_endcap = n_fail_iEmiP_endcap / nd_fail_iEmiP_endcap;

    p_pass_pMID_barrel  = n_pass_pMID_barrel  / nd_pass_pMID_barrel ;
    p_pass_pMID_endcap  = n_pass_pMID_endcap  / nd_pass_pMID_endcap ;
    p_fail_pMID_barrel  = n_fail_pMID_barrel  / nd_fail_pMID_barrel ;
    p_fail_pMID_endcap  = n_fail_pMID_endcap  / nd_fail_pMID_endcap ;

    // Preliminary fake rates
    Double_t f_pass_PFiso_barrel = n_pass_PFiso_barrel / (n_pass_PFiso_barrel + n_fail_PFiso_barrel);
    Double_t f_pass_PFiso_endcap = n_pass_PFiso_endcap / (n_pass_PFiso_endcap + n_fail_PFiso_endcap);
    Double_t f_pass_sigma_barrel = n_pass_sigma_barrel / (n_pass_sigma_barrel + n_fail_sigma_barrel);
    Double_t f_pass_sigma_endcap = n_pass_sigma_endcap / (n_pass_sigma_endcap + n_fail_sigma_endcap);
    Double_t f_pass_dEta_barrel  = n_pass_dEta_barrel  / (n_pass_dEta_barrel  + n_fail_dEta_barrel);
    Double_t f_pass_dEta_endcap  = n_pass_dEta_endcap  / (n_pass_dEta_endcap  + n_fail_dEta_endcap);
    Double_t f_pass_dPhi_barrel  = n_pass_dPhi_barrel  / (n_pass_dPhi_barrel  + n_fail_dPhi_barrel);
    Double_t f_pass_dPhi_endcap  = n_pass_dPhi_endcap  / (n_pass_dPhi_endcap  + n_fail_dPhi_endcap);
    Double_t f_pass_HovE_barrel  = n_pass_HovE_barrel  / (n_pass_HovE_barrel  + n_fail_HovE_barrel);
    Double_t f_pass_HovE_endcap  = n_pass_HovE_endcap  / (n_pass_HovE_endcap  + n_fail_HovE_endcap);
    Double_t f_pass_iEmiP_barrel = n_pass_iEmiP_barrel / (n_pass_iEmiP_barrel + n_fail_iEmiP_barrel);
    Double_t f_pass_iEmiP_endcap = n_pass_iEmiP_endcap / (n_pass_iEmiP_endcap + n_fail_iEmiP_endcap);
    Double_t f_pass_pMID_barrel  = n_pass_pMID_barrel  / (n_pass_pMID_barrel  + n_fail_pMID_barrel);
    Double_t f_pass_pMID_endcap  = n_pass_pMID_endcap  / (n_pass_pMID_endcap  + n_fail_pMID_endcap);

    cout << "------- QCD percentage -------";
    cout << "\n  relPFiso_rho:\n   Barrel:\n    Pass: " << p_pass_PFiso_barrel << "  Fail: " << p_fail_PFiso_barrel << "  FR: " << f_pass_PFiso_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_PFiso_endcap << "  Fail: " << p_fail_PFiso_endcap << "  FR: " << f_pass_PFiso_endcap << endl;
    cout << "\n  SigmaIEtaIEta:\n   Barrel:\n    Pass: " << p_pass_sigma_barrel << "  Fail: " << p_fail_sigma_barrel << "  FR: " << f_pass_sigma_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_sigma_endcap << "  Fail: " << p_fail_sigma_endcap << "  FR: " << f_pass_sigma_endcap << endl;
    cout << "\n  dEtaInSeed:\n   Barrel:\n    Pass: " << p_pass_dEta_barrel << "  Fail: " << p_fail_dEta_barrel << "  FR: " << f_pass_dEta_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_dEta_endcap << "  Fail: " << p_fail_dEta_endcap << "  FR: " << f_pass_dEta_endcap << endl;
    cout << "\n  dPhiIn:\n   Barrel:\n    Pass: " << p_pass_dPhi_barrel << "  Fail: " << p_fail_dPhi_barrel << "  FR: " << f_pass_dPhi_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_dPhi_endcap << "  Fail: " << p_fail_dPhi_endcap << "  FR: " << f_pass_dPhi_endcap << endl;
    cout << "\n  HoverE:\n   Barrel:\n    Pass: " << p_pass_HovE_barrel << "  Fail: " << p_fail_HovE_barrel << "  FR: " << f_pass_HovE_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_HovE_endcap << "  Fail: " << p_fail_HovE_endcap << "  FR: " << f_pass_HovE_endcap << endl;
    cout << "\n  InvEminusInvP:\n   Barrel:\n    Pass: " << p_pass_iEmiP_barrel << "  Fail: " << p_fail_iEmiP_barrel << "  FR: " << f_pass_iEmiP_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_iEmiP_endcap << "  Fail: " << p_fail_iEmiP_endcap << "  FR: " << f_pass_iEmiP_endcap << endl;
    cout << "\n  passMediumID:\n   Barrel:\n    Pass: " << p_pass_pMID_barrel << "  Fail: " << p_fail_pMID_barrel << "  FR: " << f_pass_pMID_barrel << endl;
    cout << "   Endcap:\n    Pass: " << p_pass_pMID_endcap << "  Fail: " << p_fail_pMID_endcap << "  FR: " << f_pass_pMID_endcap << endl;
    cout << "------------------------------" << endl;

    myRatioPlot_t *RP_PFiso_Rho_barrel_template[nPtBin_ele], *RP_PFiso_Rho_endcap_template[nPtBin_ele],
                  *RP_PFiso_Rho_barrel_jetTemplate[nPtBin_ele], *RP_PFiso_Rho_endcap_jetTemplate[nPtBin_ele];

    for (Int_t ih=0; ih<nPtBin_ele; ih++)
    {
        RP_PFiso_Rho_barrel_template[ih] = new myRatioPlot_t("RP_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10),
                                                             s_PFiso_Rho_barrel_template[ih], h_PFiso_Rho_barrel_data_template[ih]);
        RP_PFiso_Rho_barrel_jetTemplate[ih] = new myRatioPlot_t("RP_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10),
                                                                s_PFiso_Rho_barrel_jetTemplate[ih], h_PFiso_Rho_barrel_data_jetTemplate[ih]);

        RP_PFiso_Rho_barrel_template[ih]->SetPlots("relPFiso_Rho", 0, 5);
        RP_PFiso_Rho_barrel_jetTemplate[ih]->SetPlots("relPFiso_Rho", 0, 5);

        RP_PFiso_Rho_barrel_template[ih]->ImportLegend(legend);
        RP_PFiso_Rho_barrel_jetTemplate[ih]->ImportLegend(legend);

//        if (ih % 5 == 0)
//        {
//            RP_PFiso_Rho_barrel_template[ih]->Draw(1, 1e9, 0);
//            RP_PFiso_Rho_barrel_jetTemplate[ih]->Draw(1, 1e9, 0);
//        }

        if (ih < nPtBin_ele)
        {
            RP_PFiso_Rho_endcap_template[ih] = new myRatioPlot_t("RP_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10),
                                                                 s_PFiso_Rho_endcap_template[ih], h_PFiso_Rho_endcap_data_template[ih]);
            RP_PFiso_Rho_endcap_jetTemplate[ih] = new myRatioPlot_t("RP_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10),
                                                                    s_PFiso_Rho_endcap_jetTemplate[ih], h_PFiso_Rho_endcap_data_jetTemplate[ih]);

            RP_PFiso_Rho_endcap_template[ih]->SetPlots("relPFiso_Rho", 0, 5);
            RP_PFiso_Rho_endcap_jetTemplate[ih]->SetPlots("relPFiso_Rho", 0, 5);

            RP_PFiso_Rho_endcap_template[ih]->ImportLegend(legend);
            RP_PFiso_Rho_endcap_jetTemplate[ih]->ImportLegend(legend);

//            if (ih % 5 == 0)
//            {
//                RP_PFiso_Rho_endcap_template[ih]->Draw(1, 1e9, 0);
//                RP_PFiso_Rho_endcap_jetTemplate[ih]->Draw(1, 1e9, 0);
//            }
        }
    }

    cout << "Template histos drawn." << endl;

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
    cout << "Gamma+Jets integral (deno): " << h_pT_barrel_MC_deno[_GJets_Full]->Integral()+h_pT_endcap_MC_deno[_GJets_Full]->Integral() << endl;
    cout << "DY integral (deno): " << h_pT_barrel_MC_deno[_DY_Full]->Integral()+h_pT_endcap_MC_deno[_DY_Full]->Integral() << endl;
    cout << "MC integral(nume): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() +
            ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume): " << h_pT_barrel_data_nume->Integral() + h_pT_endcap_data_nume->Integral() << endl;
    cout << "QCD integral (nume): " << h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume): " << h_pT_barrel_MC_nume[_GJets_Full]->Integral()+h_pT_endcap_MC_nume[_GJets_Full]->Integral() << endl;
    cout << "DY integral (nume): " << h_pT_barrel_MC_nume[_DY_Full]->Integral()+h_pT_endcap_MC_nume[_DY_Full]->Integral() << endl;

    cout << "\nMC integral (deno barrel): " << ((TH1D*)(s_pT_barrel_deno->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (deno barrel): " << h_pT_barrel_data_deno->Integral() << endl;
    cout << "QCD integral (deno barrel): " << h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (deno barrel): " << h_pT_barrel_MC_deno[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (nume barrel): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume barrel): " << h_pT_barrel_data_nume->Integral() << endl;
    cout << "QCD integral (nume barrel): " << h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume barrel): " << h_pT_barrel_MC_nume[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (deno endcap): " << ((TH1D*)(s_pT_endcap_deno->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (deno endcap): " << h_pT_endcap_data_deno->Integral() << endl;
    cout << "QCD integral (deno endcap): " << h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (deno endcap): " << h_pT_endcap_MC_deno[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (nume endcap): " << ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume endcap): " << h_pT_endcap_data_nume->Integral() << endl;
    cout << "QCD integral (nume endcap): " << h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume endcap): " << h_pT_endcap_MC_nume[_GJets_Full]->Integral() << endl;


    // ---- TEST OF MET CUTS ---- //
    Double_t QCD_full, WJets_full;
    Double_t QCD_red[100], WJets_red[100], cuts[100];
    Double_t SSB[100];
    QCD_full = h_MET_MC[_QCDEMEnriched_Full]->Integral();
    WJets_full = h_MET_MC[_WJets_Full]->Integral();
    for (Int_t i_bin=0; i_bin<100; i_bin++)
    {
        cuts[i_bin] = i_bin + 1;
        QCD_red[i_bin] = h_MET_MC[_QCDEMEnriched_Full]->Integral(1, i_bin+1) / QCD_full;
        WJets_red[i_bin] = h_MET_MC[_WJets_Full]->Integral(1, i_bin+1) / WJets_full;
        SSB[i_bin] = QCD_red[i_bin] / (QCD_red[i_bin] + WJets_red[i_bin]);
    }
    TGraph *g_QCD_cuts = new TGraph(100, cuts, QCD_red);
    TGraph *g_WJets_cuts = new TGraph(100, cuts, WJets_red);
    TGraph *g_SSB = new TGraph(100, cuts, SSB);
    g_QCD_cuts->SetLineWidth(3);
    g_WJets_cuts->SetLineWidth(3);
    g_SSB->SetLineWidth(3);
    g_WJets_cuts->SetLineColor(kRed);
    g_SSB->SetLineColor(kBlue);

    TLegend *l_cuts = new TLegend(0.7, 0.8, 0.95, 0.9);
    l_cuts->AddEntry(g_QCD_cuts, "QCD", "l");
    l_cuts->AddEntry(g_WJets_cuts, "W+Jets", "l");
    TLegend *l_SSB = new TLegend(0.65, 0.7, 0.95, 0.87);
    l_SSB->AddEntry(g_SSB, "QCD/(WJets+QCD)", "l");

    TCanvas *c_cuts = new TCanvas("c_cuts", "MET cuts", 800, 800);
    c_cuts->SetRightMargin(0.05);
    c_cuts->SetTopMargin(0.07);
    c_cuts->SetLeftMargin(0.13);
    c_cuts->SetBottomMargin(0.13);
    g_QCD_cuts->SetTitle("Numerator");
    g_QCD_cuts->GetXaxis()->SetTitle("E_{T}^{miss} cut [GeV]");
    g_QCD_cuts->GetYaxis()->SetTitle("Reduction percentage");
    g_QCD_cuts->GetXaxis()->SetTitleSize(0.05);
    g_QCD_cuts->GetYaxis()->SetTitleSize(0.05);
    g_QCD_cuts->Draw();
    g_QCD_cuts->GetXaxis()->SetRangeUser(0, 250);
    g_QCD_cuts->GetYaxis()->SetRangeUser(0, 1);
    g_WJets_cuts->Draw("same");
    l_cuts->Draw();
    c_cuts->SetGridx();
    c_cuts->SetGridy();
    c_cuts->Update();

    TCanvas *c_SSB = new TCanvas("c_SSB", "S/(S+B) MET cuts", 800, 800);
    c_SSB->SetRightMargin(0.05);
    c_SSB->SetTopMargin(0.13);
    c_SSB->SetLeftMargin(0.13);
    c_SSB->SetBottomMargin(0.13);
    g_SSB->SetTitle("#frac{QCD}{QCD+WJets}");
    g_SSB->GetXaxis()->SetTitle("E_{T}^{miss} cut [GeV]");
    g_SSB->GetYaxis()->SetTitle("QCD/(QCD+WJets)");
    g_SSB->GetXaxis()->SetTitleSize(0.05);
    g_SSB->GetYaxis()->SetTitleSize(0.05);
    g_SSB->GetXaxis()->SetRangeUser(0, 250);
    g_SSB->GetYaxis()->SetRangeUser(0, 1);
    g_SSB->Draw();
    l_SSB->Draw();
    c_SSB->SetGridx();
    c_SSB->SetGridy();
    c_SSB->Update();
} // End of EE_HistDrawer()


/// ############################################################################# ///
/// ------------------------------- Matrix method ------------------------------- ///
/// ############################################################################# ///
void E_MatrixMethod_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    DYAnalyzer analyzer("Photon_OR");
    FileMgr fm;
    THStack *s_mass_PP = new THStack("s_mass_PP", "");
    THStack *s_mass_PF = new THStack("s_mass_PF", "");
    THStack *s_mass_FF = new THStack("s_mass_FF", "");
    THStack *s_mass_PP_TT = new THStack("s_mass_PP_TT", "");
    THStack *s_mass_PF_TT = new THStack("s_mass_PF_TT", "");
    THStack *s_mass_FF_TT = new THStack("s_mass_FF_TT", "");
    THStack *s_massSS_PP = new THStack("s_mass_PP", "");
    THStack *s_massSS_PF = new THStack("s_mass_PF", "");
    THStack *s_massSS_FF = new THStack("s_mass_FF", "");
    THStack *s_massSS_PP_TT = new THStack("s_mass_PP_TT", "");
    THStack *s_massSS_PF_TT = new THStack("s_mass_PF_TT", "");
    THStack *s_massSS_FF_TT = new THStack("s_mass_FF_TT", "");

    TH1D *h_mass_PP[_EndOf_Data_Special],
         *h_mass_PF[_EndOf_Data_Special],
         *h_mass_FF[_EndOf_Data_Special],
         *h_mass_PP_TT[_EndOf_Data_Special],
         *h_mass_PF_TT[_EndOf_Data_Special],
         *h_mass_FF_TT[_EndOf_Data_Special],
         *h_massSS_PP[_EndOf_Data_Special],
         *h_massSS_PF[_EndOf_Data_Special],
         *h_massSS_FF[_EndOf_Data_Special],
         *h_massSS_PP_TT[_EndOf_Data_Special],
         *h_massSS_PF_TT[_EndOf_Data_Special],
         *h_massSS_FF_TT[_EndOf_Data_Special];

    TFile *file = new TFile("/media/sf_DATA/FR/Electron/MatrixMethod_E.root", "READ");


//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr1 = _WW;
    while (!stop)
    {
        file->GetObject("h_mass_PP_"     +fm.Procname[pr1], h_mass_PP     [pr1]);
        file->GetObject("h_mass_PF_"     +fm.Procname[pr1], h_mass_PF     [pr1]);
        file->GetObject("h_mass_FF_"     +fm.Procname[pr1], h_mass_FF     [pr1]);
        file->GetObject("h_mass_PP_TT_"  +fm.Procname[pr1], h_mass_PP_TT  [pr1]);
        file->GetObject("h_mass_PF_TT_"  +fm.Procname[pr1], h_mass_PF_TT  [pr1]);
        file->GetObject("h_mass_FF_TT_"  +fm.Procname[pr1], h_mass_FF_TT  [pr1]);
        file->GetObject("h_massSS_PP_"   +fm.Procname[pr1], h_massSS_PP   [pr1]);
        file->GetObject("h_massSS_PF_"   +fm.Procname[pr1], h_massSS_PF   [pr1]);
        file->GetObject("h_massSS_FF_"   +fm.Procname[pr1], h_massSS_FF   [pr1]);
        file->GetObject("h_massSS_PP_TT_"+fm.Procname[pr1], h_massSS_PP_TT[pr1]);
        file->GetObject("h_massSS_PF_TT_"+fm.Procname[pr1], h_massSS_PF_TT[pr1]);
        file->GetObject("h_massSS_FF_TT_"+fm.Procname[pr1], h_massSS_FF_TT[pr1]);

        removeNegativeBins(h_mass_PP[pr1]);
        removeNegativeBins(h_mass_PF[pr1]);
        removeNegativeBins(h_mass_FF[pr1]);
        removeNegativeBins(h_mass_PP_TT[pr1]);
        removeNegativeBins(h_mass_PF_TT[pr1]);
        removeNegativeBins(h_mass_FF_TT[pr1]);
        removeNegativeBins(h_massSS_PP[pr1]);
        removeNegativeBins(h_massSS_PF[pr1]);
        removeNegativeBins(h_massSS_FF[pr1]);
        removeNegativeBins(h_massSS_PP_TT[pr1]);
        removeNegativeBins(h_massSS_PF_TT[pr1]);
        removeNegativeBins(h_massSS_FF_TT[pr1]);

        TH1D *h_mass_FP      = (TH1D*)(file->Get("h_mass_FP_"     +fm.Procname[pr1]));
        TH1D *h_mass_FP_TT   = (TH1D*)(file->Get("h_mass_FP_TT_"  +fm.Procname[pr1]));
        TH1D *h_massSS_FP    = (TH1D*)(file->Get("h_massSS_FP_"   +fm.Procname[pr1]));
        TH1D *h_massSS_FP_TT = (TH1D*)(file->Get("h_massSS_FP_TT_"+fm.Procname[pr1]));
        removeNegativeBins(h_mass_FP);
        removeNegativeBins(h_mass_FP_TT);
        removeNegativeBins(h_massSS_FP);
        removeNegativeBins(h_massSS_FP_TT);
        h_mass_PF     [pr1]->Add(h_mass_FP);
        h_mass_PF_TT  [pr1]->Add(h_mass_FP_TT);
        h_massSS_PF   [pr1]->Add(h_massSS_FP);
        h_massSS_PF_TT[pr1]->Add(h_massSS_FP_TT);

        Color_t color = kBlack;
        if (pr1 == _WJets || pr1 == _WJets_ext2v5) color = kRed - 2;
        if (pr1 == _VVnST) color = kMagenta - 5;
        if (pr1 == _WW) color = kMagenta - 5;
        if (pr1 == _WZ) color = kMagenta - 2;
        if (pr1 == _ZZ) color = kMagenta - 6;
        if (pr1 == _tbarW) color = kGreen - 2;
        if (pr1 == _tW) color = kGreen + 2;
        if (pr1 == _ttbar || pr1 == _ttbar_700to1000 || pr1 == _ttbar_1000toInf) color = kCyan + 2;

        h_mass_PP     [pr1]->SetFillColor(color);
        h_mass_PF     [pr1]->SetFillColor(color);
        h_mass_FF     [pr1]->SetFillColor(color);
        h_mass_PP_TT  [pr1]->SetFillColor(color);
        h_mass_PF_TT  [pr1]->SetFillColor(color);
        h_mass_FF_TT  [pr1]->SetFillColor(color);
        h_massSS_PP   [pr1]->SetFillColor(color);
        h_massSS_PF   [pr1]->SetFillColor(color);
        h_massSS_FF   [pr1]->SetFillColor(color);
        h_massSS_PP_TT[pr1]->SetFillColor(color);
        h_massSS_PF_TT[pr1]->SetFillColor(color);
        h_massSS_FF_TT[pr1]->SetFillColor(color);

        h_mass_PP     [pr1]->SetLineColor(color);
        h_mass_PF     [pr1]->SetLineColor(color);
        h_mass_FF     [pr1]->SetLineColor(color);
        h_mass_PP_TT  [pr1]->SetLineColor(color);
        h_mass_PF_TT  [pr1]->SetLineColor(color);
        h_mass_FF_TT  [pr1]->SetLineColor(color);
        h_massSS_PP   [pr1]->SetLineColor(color);
        h_massSS_PF   [pr1]->SetLineColor(color);
        h_massSS_FF   [pr1]->SetLineColor(color);
        h_massSS_PP_TT[pr1]->SetLineColor(color);
        h_massSS_PF_TT[pr1]->SetLineColor(color);
        h_massSS_FF_TT[pr1]->SetLineColor(color);

        h_mass_PP     [pr1]->SetDirectory(0);
        h_mass_PF     [pr1]->SetDirectory(0);
        h_mass_FF     [pr1]->SetDirectory(0);
        h_mass_PP_TT  [pr1]->SetDirectory(0);
        h_mass_PF_TT  [pr1]->SetDirectory(0);
        h_mass_FF_TT  [pr1]->SetDirectory(0);
        h_massSS_PP   [pr1]->SetDirectory(0);
        h_massSS_PF   [pr1]->SetDirectory(0);
        h_massSS_FF   [pr1]->SetDirectory(0);
        h_massSS_PP_TT[pr1]->SetDirectory(0);
        h_massSS_PF_TT[pr1]->SetDirectory(0);
        h_massSS_FF_TT[pr1]->SetDirectory(0);

        if (pr1 == _WJets)
        {
            h_mass_PP     [_WJets_Full] = ((TH1D*)(h_mass_PP     [pr1]->Clone("h_mass_PP_WJets")));
            h_mass_PF     [_WJets_Full] = ((TH1D*)(h_mass_PF     [pr1]->Clone("h_mass_PF_WJets")));
            h_mass_FF     [_WJets_Full] = ((TH1D*)(h_mass_FF     [pr1]->Clone("h_mass_FF_WJets")));
            h_mass_PP_TT  [_WJets_Full] = ((TH1D*)(h_mass_PP_TT  [pr1]->Clone("h_mass_PP_TT_WJets")));
            h_mass_PF_TT  [_WJets_Full] = ((TH1D*)(h_mass_PF_TT  [pr1]->Clone("h_mass_PF_TT_WJets")));
            h_mass_FF_TT  [_WJets_Full] = ((TH1D*)(h_mass_FF_TT  [pr1]->Clone("h_mass_FF_TT_WJets")));
            h_massSS_PP   [_WJets_Full] = ((TH1D*)(h_massSS_PP   [pr1]->Clone("h_massSS_PP_WJets")));
            h_massSS_PF   [_WJets_Full] = ((TH1D*)(h_massSS_PF   [pr1]->Clone("h_massSS_PF_WJets")));
            h_massSS_FF   [_WJets_Full] = ((TH1D*)(h_massSS_FF   [pr1]->Clone("h_massSS_FF_WJets")));
            h_massSS_PP_TT[_WJets_Full] = ((TH1D*)(h_massSS_PP_TT[pr1]->Clone("h_massSS_PP_TT_WJets")));
            h_massSS_PF_TT[_WJets_Full] = ((TH1D*)(h_massSS_PF_TT[pr1]->Clone("h_massSS_PF_TT_WJets")));
            h_massSS_FF_TT[_WJets_Full] = ((TH1D*)(h_massSS_FF_TT[pr1]->Clone("h_massSS_FF_TT_WJets")));

            h_mass_PP     [_WJets_Full]->SetDirectory(0);
            h_mass_PF     [_WJets_Full]->SetDirectory(0);
            h_mass_FF     [_WJets_Full]->SetDirectory(0);
            h_mass_PP_TT  [_WJets_Full]->SetDirectory(0);
            h_mass_PF_TT  [_WJets_Full]->SetDirectory(0);
            h_mass_FF_TT  [_WJets_Full]->SetDirectory(0);
            h_massSS_PP   [_WJets_Full]->SetDirectory(0);
            h_massSS_PF   [_WJets_Full]->SetDirectory(0);
            h_massSS_FF   [_WJets_Full]->SetDirectory(0);
            h_massSS_PP_TT[_WJets_Full]->SetDirectory(0);
            h_massSS_PF_TT[_WJets_Full]->SetDirectory(0);
            h_massSS_FF_TT[_WJets_Full]->SetDirectory(0);
        }
        else if (pr1 == _WJets_ext2v5)
        {
            h_mass_PP     [_WJets_Full]->Add(h_mass_PP     [pr1]);
            h_mass_PF     [_WJets_Full]->Add(h_mass_PF     [pr1]);
            h_mass_FF     [_WJets_Full]->Add(h_mass_FF     [pr1]);
            h_mass_PP_TT  [_WJets_Full]->Add(h_mass_PP_TT  [pr1]);
            h_mass_PF_TT  [_WJets_Full]->Add(h_mass_PF_TT  [pr1]);
            h_mass_FF_TT  [_WJets_Full]->Add(h_mass_FF_TT  [pr1]);
            h_massSS_PP   [_WJets_Full]->Add(h_massSS_PP   [pr1]);
            h_massSS_PF   [_WJets_Full]->Add(h_massSS_PF   [pr1]);
            h_massSS_FF   [_WJets_Full]->Add(h_massSS_FF   [pr1]);
            h_massSS_PP_TT[_WJets_Full]->Add(h_massSS_PP_TT[pr1]);
            h_massSS_PF_TT[_WJets_Full]->Add(h_massSS_PF_TT[pr1]);
            h_massSS_FF_TT[_WJets_Full]->Add(h_massSS_FF_TT[pr1]);
        }

        s_mass_PP     ->Add(h_mass_PP     [pr1]);
        s_mass_PF     ->Add(h_mass_PF     [pr1]);
        s_mass_FF     ->Add(h_mass_FF     [pr1]);
        s_mass_PP_TT  ->Add(h_mass_PP_TT  [pr1]);
        s_mass_PF_TT  ->Add(h_mass_PF_TT  [pr1]);
        s_mass_FF_TT  ->Add(h_mass_FF_TT  [pr1]);
        s_massSS_PP   ->Add(h_massSS_PP   [pr1]);
        s_massSS_PF   ->Add(h_massSS_PF   [pr1]);
        s_massSS_FF   ->Add(h_massSS_FF   [pr1]);
        s_massSS_PP_TT->Add(h_massSS_PP_TT[pr1]);
        s_massSS_PF_TT->Add(h_massSS_PF_TT[pr1]);
        s_massSS_FF_TT->Add(h_massSS_FF_TT[pr1]);

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

    cout << "Main bkg received" << endl;
    cout << "WJets number of FP+PF events: " << h_mass_PF[_WJets_Full]->Integral() << endl;

    // Drell-Yan
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        file->GetObject("h_mass_PP_"     +fm.Procname[pr], h_mass_PP     [pr]);
        file->GetObject("h_mass_PF_"     +fm.Procname[pr], h_mass_PF     [pr]);
        file->GetObject("h_mass_FF_"     +fm.Procname[pr], h_mass_FF     [pr]);
        file->GetObject("h_mass_PP_TT_"  +fm.Procname[pr], h_mass_PP_TT  [pr]);
        file->GetObject("h_mass_PF_TT_"  +fm.Procname[pr], h_mass_PF_TT  [pr]);
        file->GetObject("h_mass_FF_TT_"  +fm.Procname[pr], h_mass_FF_TT  [pr]);
        file->GetObject("h_massSS_PP_"   +fm.Procname[pr], h_massSS_PP   [pr]);
        file->GetObject("h_massSS_PF_"   +fm.Procname[pr], h_massSS_PF   [pr]);
        file->GetObject("h_massSS_FF_"   +fm.Procname[pr], h_massSS_FF   [pr]);
        file->GetObject("h_massSS_PP_TT_"+fm.Procname[pr], h_massSS_PP_TT[pr]);
        file->GetObject("h_massSS_PF_TT_"+fm.Procname[pr], h_massSS_PF_TT[pr]);
        file->GetObject("h_massSS_FF_TT_"+fm.Procname[pr], h_massSS_FF_TT[pr]);

        removeNegativeBins(h_mass_PP     [pr]);
        removeNegativeBins(h_mass_PF     [pr]);
        removeNegativeBins(h_mass_FF     [pr]);
        removeNegativeBins(h_mass_PP_TT  [pr]);
        removeNegativeBins(h_mass_PF_TT  [pr]);
        removeNegativeBins(h_mass_FF_TT  [pr]);
        removeNegativeBins(h_massSS_PP   [pr]);
        removeNegativeBins(h_massSS_PF   [pr]);
        removeNegativeBins(h_massSS_FF   [pr]);
        removeNegativeBins(h_massSS_PP_TT[pr]);
        removeNegativeBins(h_massSS_PF_TT[pr]);
        removeNegativeBins(h_massSS_FF_TT[pr]);

        TH1D *h_mass_FP      = (TH1D*)(file->Get("h_mass_FP_"     +fm.Procname[pr]));
        TH1D *h_mass_FP_TT   = (TH1D*)(file->Get("h_mass_FP_TT_"  +fm.Procname[pr]));
        TH1D *h_massSS_FP    = (TH1D*)(file->Get("h_massSS_FP_"   +fm.Procname[pr]));
        TH1D *h_massSS_FP_TT = (TH1D*)(file->Get("h_massSS_FP_TT_"+fm.Procname[pr]));
        removeNegativeBins(h_mass_FP);
        removeNegativeBins(h_mass_FP_TT);
        removeNegativeBins(h_massSS_FP);
        removeNegativeBins(h_massSS_FP_TT);
        h_mass_PF     [pr]->Add(h_mass_FP);
        h_mass_PF_TT  [pr]->Add(h_mass_FP_TT);
        h_massSS_PF   [pr]->Add(h_massSS_FP);
        h_massSS_PF_TT[pr]->Add(h_massSS_FP_TT);

        Color_t color = kOrange - 5;
        h_mass_PP     [pr]->SetFillColor(color);
        h_mass_PF     [pr]->SetFillColor(color);
        h_mass_FF     [pr]->SetFillColor(color);
        h_mass_PP_TT  [pr]->SetFillColor(color);
        h_mass_PF_TT  [pr]->SetFillColor(color);
        h_mass_FF_TT  [pr]->SetFillColor(color);
        h_massSS_PP   [pr]->SetFillColor(color);
        h_massSS_PF   [pr]->SetFillColor(color);
        h_massSS_FF   [pr]->SetFillColor(color);
        h_massSS_PP_TT[pr]->SetFillColor(color);
        h_massSS_PF_TT[pr]->SetFillColor(color);
        h_massSS_FF_TT[pr]->SetFillColor(color);

        h_mass_PP     [pr]->SetLineColor(color);
        h_mass_PF     [pr]->SetLineColor(color);
        h_mass_FF     [pr]->SetLineColor(color);
        h_mass_PP_TT  [pr]->SetLineColor(color);
        h_mass_PF_TT  [pr]->SetLineColor(color);
        h_mass_FF_TT  [pr]->SetLineColor(color);
        h_massSS_PP   [pr]->SetLineColor(color);
        h_massSS_PF   [pr]->SetLineColor(color);
        h_massSS_FF   [pr]->SetLineColor(color);
        h_massSS_PP_TT[pr]->SetLineColor(color);
        h_massSS_PF_TT[pr]->SetLineColor(color);
        h_massSS_FF_TT[pr]->SetLineColor(color);

        h_mass_PP     [pr]->SetDirectory(0);
        h_mass_PF     [pr]->SetDirectory(0);
        h_mass_FF     [pr]->SetDirectory(0);
        h_mass_PP_TT  [pr]->SetDirectory(0);
        h_mass_PF_TT  [pr]->SetDirectory(0);
        h_mass_FF_TT  [pr]->SetDirectory(0);
        h_massSS_PP   [pr]->SetDirectory(0);
        h_massSS_PF   [pr]->SetDirectory(0);
        h_massSS_FF   [pr]->SetDirectory(0);
        h_massSS_PP_TT[pr]->SetDirectory(0);
        h_massSS_PF_TT[pr]->SetDirectory(0);
        h_massSS_FF_TT[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_mass_PP     [_DY_Full] = ((TH1D*)(h_mass_PP     [pr]->Clone("h_mass_PP_DY")));
            h_mass_PF     [_DY_Full] = ((TH1D*)(h_mass_PF     [pr]->Clone("h_mass_PF_DY")));
            h_mass_FF     [_DY_Full] = ((TH1D*)(h_mass_FF     [pr]->Clone("h_mass_FF_DY")));
            h_mass_PP_TT  [_DY_Full] = ((TH1D*)(h_mass_PP_TT  [pr]->Clone("h_mass_PP_TT_DY")));
            h_mass_PF_TT  [_DY_Full] = ((TH1D*)(h_mass_PF_TT  [pr]->Clone("h_mass_PF_TT_DY")));
            h_mass_FF_TT  [_DY_Full] = ((TH1D*)(h_mass_FF_TT  [pr]->Clone("h_mass_FF_TT_DY")));
            h_massSS_PP   [_DY_Full] = ((TH1D*)(h_massSS_PP   [pr]->Clone("h_massSS_PP_DY")));
            h_massSS_PF   [_DY_Full] = ((TH1D*)(h_massSS_PF   [pr]->Clone("h_massSS_PF_DY")));
            h_massSS_FF   [_DY_Full] = ((TH1D*)(h_massSS_FF   [pr]->Clone("h_massSS_FF_DY")));
            h_massSS_PP_TT[_DY_Full] = ((TH1D*)(h_massSS_PP_TT[pr]->Clone("h_massSS_PP_TT_DY")));
            h_massSS_PF_TT[_DY_Full] = ((TH1D*)(h_massSS_PF_TT[pr]->Clone("h_massSS_PF_TT_DY")));
            h_massSS_FF_TT[_DY_Full] = ((TH1D*)(h_massSS_FF_TT[pr]->Clone("h_massSS_FF_TT_DY")));

            h_mass_PP     [_DY_Full]->SetDirectory(0);
            h_mass_PF     [_DY_Full]->SetDirectory(0);
            h_mass_FF     [_DY_Full]->SetDirectory(0);
            h_mass_PP_TT  [_DY_Full]->SetDirectory(0);
            h_mass_PF_TT  [_DY_Full]->SetDirectory(0);
            h_mass_FF_TT  [_DY_Full]->SetDirectory(0);
            h_massSS_PP   [_DY_Full]->SetDirectory(0);
            h_massSS_PF   [_DY_Full]->SetDirectory(0);
            h_massSS_FF   [_DY_Full]->SetDirectory(0);
            h_massSS_PP_TT[_DY_Full]->SetDirectory(0);
            h_massSS_PF_TT[_DY_Full]->SetDirectory(0);
            h_massSS_FF_TT[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_PP     [_DY_Full]->Add(h_mass_PP     [pr]);
            h_mass_PF     [_DY_Full]->Add(h_mass_PF     [pr]);
            h_mass_FF     [_DY_Full]->Add(h_mass_FF     [pr]);
            h_mass_PP_TT  [_DY_Full]->Add(h_mass_PP_TT  [pr]);
            h_mass_PF_TT  [_DY_Full]->Add(h_mass_PF_TT  [pr]);
            h_mass_FF_TT  [_DY_Full]->Add(h_mass_FF_TT  [pr]);
            h_massSS_PP   [_DY_Full]->Add(h_massSS_PP   [pr]);
            h_massSS_PF   [_DY_Full]->Add(h_massSS_PF   [pr]);
            h_massSS_FF   [_DY_Full]->Add(h_massSS_FF   [pr]);
            h_massSS_PP_TT[_DY_Full]->Add(h_massSS_PP_TT[pr]);
            h_massSS_PF_TT[_DY_Full]->Add(h_massSS_PF_TT[pr]);
            h_massSS_FF_TT[_DY_Full]->Add(h_massSS_FF_TT[pr]);
        }

        s_mass_PP     ->Add(h_mass_PP     [pr]);
        s_mass_PF     ->Add(h_mass_PF     [pr]);
        s_mass_FF     ->Add(h_mass_FF     [pr]);
        s_mass_PP_TT  ->Add(h_mass_PP_TT  [pr]);
        s_mass_PF_TT  ->Add(h_mass_PF_TT  [pr]);
        s_mass_FF_TT  ->Add(h_mass_FF_TT  [pr]);
        s_massSS_PP   ->Add(h_massSS_PP   [pr]);
        s_massSS_PF   ->Add(h_massSS_PF   [pr]);
        s_massSS_FF   ->Add(h_massSS_FF   [pr]);
        s_massSS_PP_TT->Add(h_massSS_PP_TT[pr]);
        s_massSS_PF_TT->Add(h_massSS_PF_TT[pr]);
        s_massSS_FF_TT->Add(h_massSS_FF_TT[pr]);
    }

    cout << "DY received" << endl;

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        file->GetObject("h_mass_PP_"     +fm.Procname[pr], h_mass_PP     [pr]);
        file->GetObject("h_mass_PF_"     +fm.Procname[pr], h_mass_PF     [pr]);
        file->GetObject("h_mass_FF_"     +fm.Procname[pr], h_mass_FF     [pr]);
        file->GetObject("h_mass_PP_TT_"  +fm.Procname[pr], h_mass_PP_TT  [pr]);
        file->GetObject("h_mass_PF_TT_"  +fm.Procname[pr], h_mass_PF_TT  [pr]);
        file->GetObject("h_mass_FF_TT_"  +fm.Procname[pr], h_mass_FF_TT  [pr]);
        file->GetObject("h_massSS_PP_"   +fm.Procname[pr], h_massSS_PP   [pr]);
        file->GetObject("h_massSS_PF_"   +fm.Procname[pr], h_massSS_PF   [pr]);
        file->GetObject("h_massSS_FF_"   +fm.Procname[pr], h_massSS_FF   [pr]);
        file->GetObject("h_massSS_PP_TT_"+fm.Procname[pr], h_massSS_PP_TT[pr]);
        file->GetObject("h_massSS_PF_TT_"+fm.Procname[pr], h_massSS_PF_TT[pr]);
        file->GetObject("h_massSS_FF_TT_"+fm.Procname[pr], h_massSS_FF_TT[pr]);

        removeNegativeBins(h_mass_PP     [pr]);
        removeNegativeBins(h_mass_PF     [pr]);
        removeNegativeBins(h_mass_FF     [pr]);
        removeNegativeBins(h_mass_PP_TT  [pr]);
        removeNegativeBins(h_mass_PF_TT  [pr]);
        removeNegativeBins(h_mass_FF_TT  [pr]);
        removeNegativeBins(h_massSS_PP   [pr]);
        removeNegativeBins(h_massSS_PF   [pr]);
        removeNegativeBins(h_massSS_FF   [pr]);
        removeNegativeBins(h_massSS_PP_TT[pr]);
        removeNegativeBins(h_massSS_PF_TT[pr]);
        removeNegativeBins(h_massSS_FF_TT[pr]);

        TH1D *h_mass_FP      = (TH1D*)(file->Get("h_mass_FP_"     +fm.Procname[pr]));
        TH1D *h_mass_FP_TT   = (TH1D*)(file->Get("h_mass_FP_TT_"  +fm.Procname[pr]));
        TH1D *h_massSS_FP    = (TH1D*)(file->Get("h_massSS_FP_"   +fm.Procname[pr]));
        TH1D *h_massSS_FP_TT = (TH1D*)(file->Get("h_massSS_FP_TT_"+fm.Procname[pr]));
        removeNegativeBins(h_mass_FP);
        removeNegativeBins(h_mass_FP_TT);
        removeNegativeBins(h_massSS_FP);
        removeNegativeBins(h_massSS_FP_TT);
        h_mass_PF     [pr]->Add(h_mass_FP);
        h_mass_PF_TT  [pr]->Add(h_mass_FP_TT);
        h_massSS_PF   [pr]->Add(h_massSS_FP);
        h_massSS_PF_TT[pr]->Add(h_massSS_FP_TT);

        Color_t color = kRed + 3;
        h_mass_PP     [pr]->SetFillColor(color);
        h_mass_PF     [pr]->SetFillColor(color);
        h_mass_FF     [pr]->SetFillColor(color);
        h_mass_PP_TT  [pr]->SetFillColor(color);
        h_mass_PF_TT  [pr]->SetFillColor(color);
        h_mass_FF_TT  [pr]->SetFillColor(color);
        h_massSS_PP   [pr]->SetFillColor(color);
        h_massSS_PF   [pr]->SetFillColor(color);
        h_massSS_FF   [pr]->SetFillColor(color);
        h_massSS_PP_TT[pr]->SetFillColor(color);
        h_massSS_PF_TT[pr]->SetFillColor(color);
        h_massSS_FF_TT[pr]->SetFillColor(color);

        h_mass_PP     [pr]->SetLineColor(color);
        h_mass_PF     [pr]->SetLineColor(color);
        h_mass_FF     [pr]->SetLineColor(color);
        h_mass_PP_TT  [pr]->SetLineColor(color);
        h_mass_PF_TT  [pr]->SetLineColor(color);
        h_mass_FF_TT  [pr]->SetLineColor(color);
        h_massSS_PP   [pr]->SetLineColor(color);
        h_massSS_PF   [pr]->SetLineColor(color);
        h_massSS_FF   [pr]->SetLineColor(color);
        h_massSS_PP_TT[pr]->SetLineColor(color);
        h_massSS_PF_TT[pr]->SetLineColor(color);
        h_massSS_FF_TT[pr]->SetLineColor(color);

        h_mass_PP     [pr]->SetDirectory(0);
        h_mass_PF     [pr]->SetDirectory(0);
        h_mass_FF     [pr]->SetDirectory(0);
        h_mass_PP_TT  [pr]->SetDirectory(0);
        h_mass_PF_TT  [pr]->SetDirectory(0);
        h_mass_FF_TT  [pr]->SetDirectory(0);
        h_massSS_PP   [pr]->SetDirectory(0);
        h_massSS_PF   [pr]->SetDirectory(0);
        h_massSS_FF   [pr]->SetDirectory(0);
        h_massSS_PP_TT[pr]->SetDirectory(0);
        h_massSS_PF_TT[pr]->SetDirectory(0);
        h_massSS_FF_TT[pr]->SetDirectory(0);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_mass_PP     [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_PP     [pr]->Clone("h_mass_PP_QCD")));
            h_mass_PF     [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_PF     [pr]->Clone("h_mass_PF_QCD")));
            h_mass_FF     [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_FF     [pr]->Clone("h_mass_FF_QCD")));
            h_mass_PP_TT  [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_PP_TT  [pr]->Clone("h_mass_PP_TT_QCD")));
            h_mass_PF_TT  [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_PF_TT  [pr]->Clone("h_mass_PF_TT_QCD")));
            h_mass_FF_TT  [_QCDEMEnriched_Full] = ((TH1D*)(h_mass_FF_TT  [pr]->Clone("h_mass_FF_TT_QCD")));
            h_massSS_PP   [_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_PP   [pr]->Clone("h_massSS_PP_QCD")));
            h_massSS_PF   [_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_PF   [pr]->Clone("h_massSS_PF_QCD")));
            h_massSS_FF   [_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_FF   [pr]->Clone("h_massSS_FF_QCD")));
            h_massSS_PP_TT[_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_PP_TT[pr]->Clone("h_massSS_PP_TT_QCD")));
            h_massSS_PF_TT[_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_PF_TT[pr]->Clone("h_massSS_PF_TT_QCD")));
            h_massSS_FF_TT[_QCDEMEnriched_Full] = ((TH1D*)(h_massSS_FF_TT[pr]->Clone("h_massSS_FF_TT_QCD")));

            h_mass_PP     [_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_PF     [_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_FF     [_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_PP_TT  [_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_PF_TT  [_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_FF_TT  [_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_PP   [_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_PF   [_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_FF   [_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_PP_TT[_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_PF_TT[_QCDEMEnriched_Full]->SetDirectory(0);
            h_massSS_FF_TT[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_PP     [_QCDEMEnriched_Full]->Add(h_mass_PP     [pr]);
            h_mass_PF     [_QCDEMEnriched_Full]->Add(h_mass_PF     [pr]);
            h_mass_FF     [_QCDEMEnriched_Full]->Add(h_mass_FF     [pr]);
            h_mass_PP_TT  [_QCDEMEnriched_Full]->Add(h_mass_PP_TT  [pr]);
            h_mass_PF_TT  [_QCDEMEnriched_Full]->Add(h_mass_PF_TT  [pr]);
            h_mass_FF_TT  [_QCDEMEnriched_Full]->Add(h_mass_FF_TT  [pr]);
            h_massSS_PP   [_QCDEMEnriched_Full]->Add(h_massSS_PP   [pr]);
            h_massSS_PF   [_QCDEMEnriched_Full]->Add(h_massSS_PF   [pr]);
            h_massSS_FF   [_QCDEMEnriched_Full]->Add(h_massSS_FF   [pr]);
            h_massSS_PP_TT[_QCDEMEnriched_Full]->Add(h_massSS_PP_TT[pr]);
            h_massSS_PF_TT[_QCDEMEnriched_Full]->Add(h_massSS_PF_TT[pr]);
            h_massSS_FF_TT[_QCDEMEnriched_Full]->Add(h_massSS_FF_TT[pr]);
        }

        s_mass_PP     ->Add(h_mass_PP     [pr]);
        s_mass_PF     ->Add(h_mass_PF     [pr]);
        s_mass_FF     ->Add(h_mass_FF     [pr]);
        s_mass_PP_TT  ->Add(h_mass_PP_TT  [pr]);
        s_mass_PF_TT  ->Add(h_mass_PF_TT  [pr]);
        s_mass_FF_TT  ->Add(h_mass_FF_TT  [pr]);
        s_massSS_PP   ->Add(h_massSS_PP   [pr]);
        s_massSS_PF   ->Add(h_massSS_PF   [pr]);
        s_massSS_FF   ->Add(h_massSS_FF   [pr]);
        s_massSS_PP_TT->Add(h_massSS_PP_TT[pr]);
        s_massSS_PF_TT->Add(h_massSS_PF_TT[pr]);
        s_massSS_FF_TT->Add(h_massSS_FF_TT[pr]);
    }

    cout << "QCD received" << endl;

    // GammaJets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        file->GetObject("h_mass_PP_"     +fm.Procname[pr], h_mass_PP     [pr]);
        file->GetObject("h_mass_PF_"     +fm.Procname[pr], h_mass_PF     [pr]);
        file->GetObject("h_mass_FF_"     +fm.Procname[pr], h_mass_FF     [pr]);
        file->GetObject("h_mass_PP_TT_"  +fm.Procname[pr], h_mass_PP_TT  [pr]);
        file->GetObject("h_mass_PF_TT_"  +fm.Procname[pr], h_mass_PF_TT  [pr]);
        file->GetObject("h_mass_FF_TT_"  +fm.Procname[pr], h_mass_FF_TT  [pr]);
        file->GetObject("h_massSS_PP_"   +fm.Procname[pr], h_massSS_PP   [pr]);
        file->GetObject("h_massSS_PF_"   +fm.Procname[pr], h_massSS_PF   [pr]);
        file->GetObject("h_massSS_FF_"   +fm.Procname[pr], h_massSS_FF   [pr]);
        file->GetObject("h_massSS_PP_TT_"+fm.Procname[pr], h_massSS_PP_TT[pr]);
        file->GetObject("h_massSS_PF_TT_"+fm.Procname[pr], h_massSS_PF_TT[pr]);
        file->GetObject("h_massSS_FF_TT_"+fm.Procname[pr], h_massSS_FF_TT[pr]);

        removeNegativeBins(h_mass_PP     [pr]);
        removeNegativeBins(h_mass_PF     [pr]);
        removeNegativeBins(h_mass_FF     [pr]);
        removeNegativeBins(h_mass_PP_TT  [pr]);
        removeNegativeBins(h_mass_PF_TT  [pr]);
        removeNegativeBins(h_mass_FF_TT  [pr]);
        removeNegativeBins(h_massSS_PP   [pr]);
        removeNegativeBins(h_massSS_PF   [pr]);
        removeNegativeBins(h_massSS_FF   [pr]);
        removeNegativeBins(h_massSS_PP_TT[pr]);
        removeNegativeBins(h_massSS_PF_TT[pr]);
        removeNegativeBins(h_massSS_FF_TT[pr]);

        TH1D *h_mass_FP      = (TH1D*)(file->Get("h_mass_FP_"     +fm.Procname[pr]));
        TH1D *h_mass_FP_TT   = (TH1D*)(file->Get("h_mass_FP_TT_"  +fm.Procname[pr]));
        TH1D *h_massSS_FP    = (TH1D*)(file->Get("h_massSS_FP_"   +fm.Procname[pr]));
        TH1D *h_massSS_FP_TT = (TH1D*)(file->Get("h_massSS_FP_TT_"+fm.Procname[pr]));
        removeNegativeBins(h_mass_FP);
        removeNegativeBins(h_mass_FP_TT);
        removeNegativeBins(h_massSS_FP);
        removeNegativeBins(h_massSS_FP_TT);
        h_mass_PF     [pr]->Add(h_mass_FP);
        h_mass_PF_TT  [pr]->Add(h_mass_FP_TT);
        h_massSS_PF   [pr]->Add(h_massSS_FP);
        h_massSS_PF_TT[pr]->Add(h_massSS_FP_TT);

        Color_t color = kYellow + 3;
        h_mass_PP     [pr]->SetFillColor(color);
        h_mass_PF     [pr]->SetFillColor(color);
        h_mass_FF     [pr]->SetFillColor(color);
        h_mass_PP_TT  [pr]->SetFillColor(color);
        h_mass_PF_TT  [pr]->SetFillColor(color);
        h_mass_FF_TT  [pr]->SetFillColor(color);
        h_massSS_PP   [pr]->SetFillColor(color);
        h_massSS_PF   [pr]->SetFillColor(color);
        h_massSS_FF   [pr]->SetFillColor(color);
        h_massSS_PP_TT[pr]->SetFillColor(color);
        h_massSS_PF_TT[pr]->SetFillColor(color);
        h_massSS_FF_TT[pr]->SetFillColor(color);

        h_mass_PP     [pr]->SetLineColor(color);
        h_mass_PF     [pr]->SetLineColor(color);
        h_mass_FF     [pr]->SetLineColor(color);
        h_mass_PP_TT  [pr]->SetLineColor(color);
        h_mass_PF_TT  [pr]->SetLineColor(color);
        h_mass_FF_TT  [pr]->SetLineColor(color);
        h_massSS_PP   [pr]->SetLineColor(color);
        h_massSS_PF   [pr]->SetLineColor(color);
        h_massSS_FF   [pr]->SetLineColor(color);
        h_massSS_PP_TT[pr]->SetLineColor(color);
        h_massSS_PF_TT[pr]->SetLineColor(color);
        h_massSS_FF_TT[pr]->SetLineColor(color);

        h_mass_PP     [pr]->SetDirectory(0);
        h_mass_PF     [pr]->SetDirectory(0);
        h_mass_FF     [pr]->SetDirectory(0);
        h_mass_PP_TT  [pr]->SetDirectory(0);
        h_mass_PF_TT  [pr]->SetDirectory(0);
        h_mass_FF_TT  [pr]->SetDirectory(0);
        h_massSS_PP   [pr]->SetDirectory(0);
        h_massSS_PF   [pr]->SetDirectory(0);
        h_massSS_FF   [pr]->SetDirectory(0);
        h_massSS_PP_TT[pr]->SetDirectory(0);
        h_massSS_PF_TT[pr]->SetDirectory(0);
        h_massSS_FF_TT[pr]->SetDirectory(0);

        if (pr == _GJets_20to100)
        {
            h_mass_PP     [_GJets_Full] = ((TH1D*)(h_mass_PP     [pr]->Clone("h_mass_PP_GJets")));
            h_mass_PF     [_GJets_Full] = ((TH1D*)(h_mass_PF     [pr]->Clone("h_mass_PF_GJets")));
            h_mass_FF     [_GJets_Full] = ((TH1D*)(h_mass_FF     [pr]->Clone("h_mass_FF_GJets")));
            h_mass_PP_TT  [_GJets_Full] = ((TH1D*)(h_mass_PP_TT  [pr]->Clone("h_mass_PP_TT_GJets")));
            h_mass_PF_TT  [_GJets_Full] = ((TH1D*)(h_mass_PF_TT  [pr]->Clone("h_mass_PF_TT_GJets")));
            h_mass_FF_TT  [_GJets_Full] = ((TH1D*)(h_mass_FF_TT  [pr]->Clone("h_mass_FF_TT_GJets")));
            h_massSS_PP   [_GJets_Full] = ((TH1D*)(h_massSS_PP   [pr]->Clone("h_massSS_PP_GJets")));
            h_massSS_PF   [_GJets_Full] = ((TH1D*)(h_massSS_PF   [pr]->Clone("h_massSS_PF_GJets")));
            h_massSS_FF   [_GJets_Full] = ((TH1D*)(h_massSS_FF   [pr]->Clone("h_massSS_FF_GJets")));
            h_massSS_PP_TT[_GJets_Full] = ((TH1D*)(h_massSS_PP_TT[pr]->Clone("h_massSS_PP_TT_GJets")));
            h_massSS_PF_TT[_GJets_Full] = ((TH1D*)(h_massSS_PF_TT[pr]->Clone("h_massSS_PF_TT_GJets")));
            h_massSS_FF_TT[_GJets_Full] = ((TH1D*)(h_massSS_FF_TT[pr]->Clone("h_massSS_FF_TT_GJets")));

            h_mass_PP     [_GJets_Full]->SetDirectory(0);
            h_mass_PF     [_GJets_Full]->SetDirectory(0);
            h_mass_FF     [_GJets_Full]->SetDirectory(0);
            h_mass_PP_TT  [_GJets_Full]->SetDirectory(0);
            h_mass_PF_TT  [_GJets_Full]->SetDirectory(0);
            h_mass_FF_TT  [_GJets_Full]->SetDirectory(0);
            h_massSS_PP   [_GJets_Full]->SetDirectory(0);
            h_massSS_PF   [_GJets_Full]->SetDirectory(0);
            h_massSS_FF   [_GJets_Full]->SetDirectory(0);
            h_massSS_PP_TT[_GJets_Full]->SetDirectory(0);
            h_massSS_PF_TT[_GJets_Full]->SetDirectory(0);
            h_massSS_FF_TT[_GJets_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_PP     [_GJets_Full]->Add(h_mass_PP     [pr]);
            h_mass_PF     [_GJets_Full]->Add(h_mass_PF     [pr]);
            h_mass_FF     [_GJets_Full]->Add(h_mass_FF     [pr]);
            h_mass_PP_TT  [_GJets_Full]->Add(h_mass_PP_TT  [pr]);
            h_mass_PF_TT  [_GJets_Full]->Add(h_mass_PF_TT  [pr]);
            h_mass_FF_TT  [_GJets_Full]->Add(h_mass_FF_TT  [pr]);
            h_massSS_PP   [_GJets_Full]->Add(h_massSS_PP   [pr]);
            h_massSS_PF   [_GJets_Full]->Add(h_massSS_PF   [pr]);
            h_massSS_FF   [_GJets_Full]->Add(h_massSS_FF   [pr]);
            h_massSS_PP_TT[_GJets_Full]->Add(h_massSS_PP_TT[pr]);
            h_massSS_PF_TT[_GJets_Full]->Add(h_massSS_PF_TT[pr]);
            h_massSS_FF_TT[_GJets_Full]->Add(h_massSS_FF_TT[pr]);
        }

        s_mass_PP     ->Add(h_mass_PP     [pr]);
        s_mass_PF     ->Add(h_mass_PF     [pr]);
        s_mass_FF     ->Add(h_mass_FF     [pr]);
        s_mass_PP_TT  ->Add(h_mass_PP_TT  [pr]);
        s_mass_PF_TT  ->Add(h_mass_PF_TT  [pr]);
        s_mass_FF_TT  ->Add(h_mass_FF_TT  [pr]);
        s_massSS_PP   ->Add(h_massSS_PP   [pr]);
        s_massSS_PF   ->Add(h_massSS_PF   [pr]);
        s_massSS_FF   ->Add(h_massSS_FF   [pr]);
        s_massSS_PP_TT->Add(h_massSS_PP_TT[pr]);
        s_massSS_PF_TT->Add(h_massSS_PF_TT[pr]);
        s_massSS_FF_TT->Add(h_massSS_FF_TT[pr]);
    }

    cout << "GJets received" << endl;

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_DoubleEG_B; pr<=_DoubleEG_H; pr=next(pr))
    {
        file->GetObject("h_mass_PP_"     +fm.Procname[pr], h_mass_PP     [pr]);
        file->GetObject("h_mass_PF_"     +fm.Procname[pr], h_mass_PF     [pr]);
        file->GetObject("h_mass_FF_"     +fm.Procname[pr], h_mass_FF     [pr]);
        file->GetObject("h_mass_PP_TT_"  +fm.Procname[pr], h_mass_PP_TT  [pr]);
        file->GetObject("h_mass_PF_TT_"  +fm.Procname[pr], h_mass_PF_TT  [pr]);
        file->GetObject("h_mass_FF_TT_"  +fm.Procname[pr], h_mass_FF_TT  [pr]);
        file->GetObject("h_massSS_PP_"   +fm.Procname[pr], h_massSS_PP   [pr]);
        file->GetObject("h_massSS_PF_"   +fm.Procname[pr], h_massSS_PF   [pr]);
        file->GetObject("h_massSS_FF_"   +fm.Procname[pr], h_massSS_FF   [pr]);
        file->GetObject("h_massSS_PP_TT_"+fm.Procname[pr], h_massSS_PP_TT[pr]);
        file->GetObject("h_massSS_PF_TT_"+fm.Procname[pr], h_massSS_PF_TT[pr]);
        file->GetObject("h_massSS_FF_TT_"+fm.Procname[pr], h_massSS_FF_TT[pr]);

        removeNegativeBins(h_mass_PP     [pr]);
        removeNegativeBins(h_mass_PF     [pr]);
        removeNegativeBins(h_mass_FF     [pr]);
        removeNegativeBins(h_mass_PP_TT  [pr]);
        removeNegativeBins(h_mass_PF_TT  [pr]);
        removeNegativeBins(h_mass_FF_TT  [pr]);
        removeNegativeBins(h_massSS_PP   [pr]);
        removeNegativeBins(h_massSS_PF   [pr]);
        removeNegativeBins(h_massSS_FF   [pr]);
        removeNegativeBins(h_massSS_PP_TT[pr]);
        removeNegativeBins(h_massSS_PF_TT[pr]);
        removeNegativeBins(h_massSS_FF_TT[pr]);

        TH1D *h_mass_FP      = (TH1D*)(file->Get("h_mass_FP_"     +fm.Procname[pr]));
        TH1D *h_mass_FP_TT   = (TH1D*)(file->Get("h_mass_FP_TT_"  +fm.Procname[pr]));
        TH1D *h_massSS_FP    = (TH1D*)(file->Get("h_massSS_FP_"   +fm.Procname[pr]));
        TH1D *h_massSS_FP_TT = (TH1D*)(file->Get("h_massSS_FP_TT_"+fm.Procname[pr]));
        removeNegativeBins(h_mass_FP);
        removeNegativeBins(h_mass_FP_TT);
        removeNegativeBins(h_massSS_FP);
        removeNegativeBins(h_massSS_FP_TT);
        h_mass_PF     [pr]->Add(h_mass_FP);
        h_mass_PF_TT  [pr]->Add(h_mass_FP_TT);
        h_massSS_PF   [pr]->Add(h_massSS_FP);
        h_massSS_PF_TT[pr]->Add(h_massSS_FP_TT);

        h_mass_PP     [pr]->SetMarkerStyle(kFullDotLarge);
        h_mass_PF     [pr]->SetMarkerStyle(kFullDotLarge);
        h_mass_FF     [pr]->SetMarkerStyle(kFullDotLarge);
        h_mass_PP_TT  [pr]->SetMarkerStyle(kFullDotLarge);
        h_mass_PF_TT  [pr]->SetMarkerStyle(kFullDotLarge);
        h_mass_FF_TT  [pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_PP   [pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_PF   [pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_FF   [pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_PP_TT[pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_PF_TT[pr]->SetMarkerStyle(kFullDotLarge);
        h_massSS_FF_TT[pr]->SetMarkerStyle(kFullDotLarge);

        h_mass_PP     [pr]->SetMarkerColor(kBlack);
        h_mass_PF     [pr]->SetMarkerColor(kBlack);
        h_mass_FF     [pr]->SetMarkerColor(kBlack);
        h_mass_PP_TT  [pr]->SetMarkerColor(kBlack);
        h_mass_PF_TT  [pr]->SetMarkerColor(kBlack);
        h_mass_FF_TT  [pr]->SetMarkerColor(kBlack);
        h_massSS_PP   [pr]->SetMarkerColor(kBlack);
        h_massSS_PF   [pr]->SetMarkerColor(kBlack);
        h_massSS_FF   [pr]->SetMarkerColor(kBlack);
        h_massSS_PP_TT[pr]->SetMarkerColor(kBlack);
        h_massSS_PF_TT[pr]->SetMarkerColor(kBlack);
        h_massSS_FF_TT[pr]->SetMarkerColor(kBlack);

        h_mass_PP     [pr]->SetLineColor(kBlack);
        h_mass_PF     [pr]->SetLineColor(kBlack);
        h_mass_FF     [pr]->SetLineColor(kBlack);
        h_mass_PP_TT  [pr]->SetLineColor(kBlack);
        h_mass_PF_TT  [pr]->SetLineColor(kBlack);
        h_mass_FF_TT  [pr]->SetLineColor(kBlack);
        h_massSS_PP   [pr]->SetLineColor(kBlack);
        h_massSS_PF   [pr]->SetLineColor(kBlack);
        h_massSS_FF   [pr]->SetLineColor(kBlack);
        h_massSS_PP_TT[pr]->SetLineColor(kBlack);
        h_massSS_PF_TT[pr]->SetLineColor(kBlack);
        h_massSS_FF_TT[pr]->SetLineColor(kBlack);

        h_mass_PP     [pr]->SetDirectory(0);
        h_mass_PF     [pr]->SetDirectory(0);
        h_mass_FF     [pr]->SetDirectory(0);
        h_mass_PP_TT  [pr]->SetDirectory(0);
        h_mass_PF_TT  [pr]->SetDirectory(0);
        h_mass_FF_TT  [pr]->SetDirectory(0);
        h_massSS_PP   [pr]->SetDirectory(0);
        h_massSS_PF   [pr]->SetDirectory(0);
        h_massSS_FF   [pr]->SetDirectory(0);
        h_massSS_PP_TT[pr]->SetDirectory(0);
        h_massSS_PF_TT[pr]->SetDirectory(0);
        h_massSS_FF_TT[pr]->SetDirectory(0);

        if (pr == _DoubleEG_B)
        {
            h_mass_PP     [_DoubleEG_Full] = ((TH1D*)(h_mass_PP     [pr]->Clone("h_mass_PP_DoubleEG")));
            h_mass_PF     [_DoubleEG_Full] = ((TH1D*)(h_mass_PF     [pr]->Clone("h_mass_PF_DoubleEG")));
            h_mass_FF     [_DoubleEG_Full] = ((TH1D*)(h_mass_FF     [pr]->Clone("h_mass_FF_DoubleEG")));
            h_mass_PP_TT  [_DoubleEG_Full] = ((TH1D*)(h_mass_PP_TT  [pr]->Clone("h_mass_PP_TT_DoubleEG")));
            h_mass_PF_TT  [_DoubleEG_Full] = ((TH1D*)(h_mass_PF_TT  [pr]->Clone("h_mass_PF_TT_DoubleEG")));
            h_mass_FF_TT  [_DoubleEG_Full] = ((TH1D*)(h_mass_FF_TT  [pr]->Clone("h_mass_FF_TT_DoubleEG")));
            h_massSS_PP   [_DoubleEG_Full] = ((TH1D*)(h_massSS_PP   [pr]->Clone("h_massSS_PP_DoubleEG")));
            h_massSS_PF   [_DoubleEG_Full] = ((TH1D*)(h_massSS_PF   [pr]->Clone("h_massSS_PF_DoubleEG")));
            h_massSS_FF   [_DoubleEG_Full] = ((TH1D*)(h_massSS_FF   [pr]->Clone("h_massSS_FF_DoubleEG")));
            h_massSS_PP_TT[_DoubleEG_Full] = ((TH1D*)(h_massSS_PP_TT[pr]->Clone("h_massSS_PP_TT_DoubleEG")));
            h_massSS_PF_TT[_DoubleEG_Full] = ((TH1D*)(h_massSS_PF_TT[pr]->Clone("h_massSS_PF_TT_DoubleEG")));
            h_massSS_FF_TT[_DoubleEG_Full] = ((TH1D*)(h_massSS_FF_TT[pr]->Clone("h_massSS_FF_TT_DoubleEG")));

            h_mass_PP     [_DoubleEG_Full]->SetDirectory(0);
            h_mass_PF     [_DoubleEG_Full]->SetDirectory(0);
            h_mass_FF     [_DoubleEG_Full]->SetDirectory(0);
            h_mass_PP_TT  [_DoubleEG_Full]->SetDirectory(0);
            h_mass_PF_TT  [_DoubleEG_Full]->SetDirectory(0);
            h_mass_FF_TT  [_DoubleEG_Full]->SetDirectory(0);
            h_massSS_PP   [_DoubleEG_Full]->SetDirectory(0);
            h_massSS_PF   [_DoubleEG_Full]->SetDirectory(0);
            h_massSS_FF   [_DoubleEG_Full]->SetDirectory(0);
            h_massSS_PP_TT[_DoubleEG_Full]->SetDirectory(0);
            h_massSS_PF_TT[_DoubleEG_Full]->SetDirectory(0);
            h_massSS_FF_TT[_DoubleEG_Full]->SetDirectory(0);
        }
        else
        {
            h_mass_PP     [_DoubleEG_Full]->Add(h_mass_PP     [pr]);
            h_mass_PF     [_DoubleEG_Full]->Add(h_mass_PF     [pr]);
            h_mass_FF     [_DoubleEG_Full]->Add(h_mass_FF     [pr]);
            h_mass_PP_TT  [_DoubleEG_Full]->Add(h_mass_PP_TT  [pr]);
            h_mass_PF_TT  [_DoubleEG_Full]->Add(h_mass_PF_TT  [pr]);
            h_mass_FF_TT  [_DoubleEG_Full]->Add(h_mass_FF_TT  [pr]);
            h_massSS_PP   [_DoubleEG_Full]->Add(h_massSS_PP   [pr]);
            h_massSS_PF   [_DoubleEG_Full]->Add(h_massSS_PF   [pr]);
            h_massSS_FF   [_DoubleEG_Full]->Add(h_massSS_FF   [pr]);
            h_massSS_PP_TT[_DoubleEG_Full]->Add(h_massSS_PP_TT[pr]);
            h_massSS_PF_TT[_DoubleEG_Full]->Add(h_massSS_PF_TT[pr]);
            h_massSS_FF_TT[_DoubleEG_Full]->Add(h_massSS_FF_TT[pr]);
        }
    }

    file->Close();

    cout << "Data received" << endl;

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_mass_PP      = new myRatioPlot_t("RP_mass_PP     ", s_mass_PP     , h_mass_PP     [_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_PF      = new myRatioPlot_t("RP_mass_PF     ", s_mass_PF     , h_mass_PF     [_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_FF      = new myRatioPlot_t("RP_mass_FF     ", s_mass_FF     , h_mass_FF     [_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_PP_TT   = new myRatioPlot_t("RP_mass_PP_TT  ", s_mass_PP_TT  , h_mass_PP_TT  [_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_PF_TT   = new myRatioPlot_t("RP_mass_PF_TT  ", s_mass_PF_TT  , h_mass_PF_TT  [_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_FF_TT   = new myRatioPlot_t("RP_mass_FF_TT  ", s_mass_FF_TT  , h_mass_FF_TT  [_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_PP    = new myRatioPlot_t("RP_massSS_PP   ", s_massSS_PP   , h_massSS_PP   [_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_PF    = new myRatioPlot_t("RP_massSS_PF   ", s_massSS_PF   , h_massSS_PF   [_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_FF    = new myRatioPlot_t("RP_massSS_FF   ", s_massSS_FF   , h_massSS_FF   [_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_PP_TT = new myRatioPlot_t("RP_massSS_PP_TT", s_massSS_PP_TT, h_massSS_PP_TT[_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_PF_TT = new myRatioPlot_t("RP_massSS_PF_TT", s_massSS_PF_TT, h_massSS_PF_TT[_DoubleEG_Full]);
    myRatioPlot_t *RP_massSS_FF_TT = new myRatioPlot_t("RP_massSS_FF_TT", s_massSS_FF_TT, h_massSS_FF_TT[_DoubleEG_Full]);

    RP_mass_PP     ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_mass_PF     ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_mass_FF     ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_mass_PP_TT  ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_mass_PF_TT  ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_mass_FF_TT  ->SetPlots("m_{ee} [GeV/c^{2}]", 15, 3000);
    RP_massSS_PP   ->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);
    RP_massSS_PF   ->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);
    RP_massSS_FF   ->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);
    RP_massSS_PP_TT->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);
    RP_massSS_PF_TT->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);
    RP_massSS_FF_TT->SetPlots("m_{ee} [GeV/c^{2}] (same-sign)", 15, 3000);

    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_mass_PP[_DoubleEG_Full], "Data", "lp");
    legend->AddEntry(h_mass_PP[_DY_50to100], "DY", "f");
    legend->AddEntry(h_mass_PP[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_mass_PP[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_mass_PP[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_mass_PP[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_mass_PP[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_mass_PP[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend->AddEntry(h_mass_PP[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_mass_PP[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->AddEntry(h_mass_PP[_GJets_20to100], "#gamma+Jets", "f");
    legend->SetNColumns(2);

    RP_mass_PP     ->ImportLegend(legend);
    RP_mass_PF     ->ImportLegend(legend);
    RP_mass_FF     ->ImportLegend(legend);
    RP_mass_PP_TT  ->ImportLegend(legend);
    RP_mass_PF_TT  ->ImportLegend(legend);
    RP_mass_FF_TT  ->ImportLegend(legend);
    RP_massSS_PP   ->ImportLegend(legend);
    RP_massSS_PF   ->ImportLegend(legend);
    RP_massSS_FF   ->ImportLegend(legend);
    RP_massSS_PP_TT->ImportLegend(legend);
    RP_massSS_PF_TT->ImportLegend(legend);
    RP_massSS_FF_TT->ImportLegend(legend);

    RP_mass_PP     ->Draw(1, 1e9, 1, "HIST", "N_{PP}");
    RP_mass_PF     ->Draw(1, 1e9, 1, "HIST", "N_{PF+FP}");
    RP_mass_FF     ->Draw(1, 1e9, 1, "HIST", "N_{FF}");
    RP_mass_PP_TT  ->Draw(1, 1e9, 1, "HIST", "N_{PP}^{TT}");
    RP_mass_PF_TT  ->Draw(1, 1e9, 1, "HIST", "N_{PF+FP}^{TT}");
    RP_mass_FF_TT  ->Draw(1, 1e9, 1, "HIST", "N_{FF}^{TT}");
    RP_massSS_PP   ->Draw(1, 1e9, 1, "HIST", "N_{PP}");
    RP_massSS_PF   ->Draw(1, 1e9, 1, "HIST", "N_{PF+FP}");
    RP_massSS_FF   ->Draw(1, 1e9, 1, "HIST", "N_{FF}");
    RP_massSS_PP_TT->Draw(1, 1e9, 1, "HIST", "N_{PP}^{TT}");
    RP_massSS_PF_TT->Draw(1, 1e9, 1, "HIST", "N_{PF+FP}^{TT}");
    RP_massSS_FF_TT->Draw(1, 1e9, 1, "HIST", "N_{FF}^{TT}");

    cout << "Main histos drawn." << endl;


//    cout << "MC integral (deno): " << ((TH1D*)(s_pT_barrel_deno->GetStack()->Last()))->Integral() +
//                                     ((TH1D*)(s_pT_endcap_deno->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (deno): " << h_pT_barrel_data_deno->Integral() + h_pT_endcap_data_deno->Integral() << endl;
//    cout << "QCD integral (deno): " << h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (deno): " << h_pT_barrel_MC_deno[_GJets_Full]->Integral()+h_pT_endcap_MC_deno[_GJets_Full]->Integral() << endl;
//    cout << "DY integral (deno): " << h_pT_barrel_MC_deno[_DY_Full]->Integral()+h_pT_endcap_MC_deno[_DY_Full]->Integral() << endl;
//    cout << "MC integral(nume): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() +
//            ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (nume): " << h_pT_barrel_data_nume->Integral() + h_pT_endcap_data_nume->Integral() << endl;
//    cout << "QCD integral (nume): " << h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (nume): " << h_pT_barrel_MC_nume[_GJets_Full]->Integral()+h_pT_endcap_MC_nume[_GJets_Full]->Integral() << endl;
//    cout << "DY integral (nume): " << h_pT_barrel_MC_nume[_DY_Full]->Integral()+h_pT_endcap_MC_nume[_DY_Full]->Integral() << endl;

//    cout << "\nMC integral (deno barrel): " << ((TH1D*)(s_pT_barrel_deno->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (deno barrel): " << h_pT_barrel_data_deno->Integral() << endl;
//    cout << "QCD integral (deno barrel): " << h_pT_barrel_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (deno barrel): " << h_pT_barrel_MC_deno[_GJets_Full]->Integral() << endl;
//    cout << "\nMC integral (nume barrel): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (nume barrel): " << h_pT_barrel_data_nume->Integral() << endl;
//    cout << "QCD integral (nume barrel): " << h_pT_barrel_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (nume barrel): " << h_pT_barrel_MC_nume[_GJets_Full]->Integral() << endl;
//    cout << "\nMC integral (deno endcap): " << ((TH1D*)(s_pT_endcap_deno->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (deno endcap): " << h_pT_endcap_data_deno->Integral() << endl;
//    cout << "QCD integral (deno endcap): " << h_pT_endcap_MC_deno[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (deno endcap): " << h_pT_endcap_MC_deno[_GJets_Full]->Integral() << endl;
//    cout << "\nMC integral (nume endcap): " << ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
//    cout << "Data integral (nume endcap): " << h_pT_endcap_data_nume->Integral() << endl;
//    cout << "QCD integral (nume endcap): " << h_pT_endcap_MC_nume[_QCDEMEnriched_Full]->Integral() << endl;
//    cout << "Gamma+Jets integral (nume endcap): " << h_pT_endcap_MC_nume[_GJets_Full]->Integral() << endl;

} // End of E_MatrixMethod_HistDrawer()


void E_HistDrawer_alt()
{
    TFile *f = new TFile("/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root", "READ");

    TH1D *h_pT_barrel_lead_pass[4],
         *h_pT_endcap_lead_pass[4],
         *h_pT_barrel_lead_fail[4],
         *h_pT_endcap_lead_fail[4],
         *h_pT_barrel_sub_pass[4],
         *h_pT_endcap_sub_pass[4],
         *h_pT_barrel_sub_fail[4],
         *h_pT_endcap_sub_fail[4],
         *h_MET_fail[4],
         *h_MT_fail[4];
    THStack *s_pT_barrel_lead_pass = new THStack("s_pT_barrel_lead_pass", "");
    THStack *s_pT_endcap_lead_pass = new THStack("s_pT_endcap_lead_pass", "");
    THStack *s_pT_barrel_lead_fail = new THStack("s_pT_barrel_lead_fail", "");
    THStack *s_pT_endcap_lead_fail = new THStack("s_pT_endcap_lead_fail", "");
    THStack *s_pT_barrel_sub_pass = new THStack("s_pT_barrel_sub_pass", "");
    THStack *s_pT_endcap_sub_pass = new THStack("s_pT_endcap_sub_pass", "");
    THStack *s_pT_barrel_sub_fail = new THStack("s_pT_barrel_sub_fail", "");
    THStack *s_pT_endcap_sub_fail = new THStack("s_pT_endcap_sub_fail", "");
    THStack *s_MET_fail = new THStack("s_MET_fail", "");
    THStack *s_MT_fail = new THStack("s_MT_fail", "");
    Color_t color = kBlack;
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
        f->GetObject("h_MET_fail_"+type[i], h_MET_fail[i]);
        f->GetObject("h_MT_fail_"+type[i], h_MT_fail[i]);
        h_pT_barrel_lead_pass[i]->SetDirectory(0);
        h_pT_endcap_lead_pass[i]->SetDirectory(0);
        h_pT_barrel_lead_fail[i]->SetDirectory(0);
        h_pT_endcap_lead_fail[i]->SetDirectory(0);
        h_pT_barrel_sub_pass[i]->SetDirectory(0);
        h_pT_endcap_sub_pass[i]->SetDirectory(0);
        h_pT_barrel_sub_fail[i]->SetDirectory(0);
        h_pT_endcap_sub_fail[i]->SetDirectory(0);
        h_MET_fail[i]->SetDirectory(0);
        h_MT_fail[i]->SetDirectory(0);
        if (i == 0)
        {
            h_pT_barrel_lead_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_lead_pass[i]->SetMarkerColor(kBlack);
            h_pT_barrel_lead_pass[i]->SetLineColor(kBlack);
            h_pT_endcap_lead_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_lead_pass[i]->SetMarkerColor(kBlack);
            h_pT_endcap_lead_pass[i]->SetLineColor(kBlack);
            h_pT_barrel_lead_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_lead_fail[i]->SetMarkerColor(kBlack);
            h_pT_barrel_lead_fail[i]->SetLineColor(kBlack);
            h_pT_endcap_lead_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_lead_fail[i]->SetMarkerColor(kBlack);
            h_pT_endcap_lead_fail[i]->SetLineColor(kBlack);
            h_pT_barrel_sub_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_sub_pass[i]->SetMarkerColor(kBlack);
            h_pT_barrel_sub_pass[i]->SetLineColor(kBlack);
            h_pT_endcap_sub_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_sub_pass[i]->SetMarkerColor(kBlack);
            h_pT_endcap_sub_pass[i]->SetLineColor(kBlack);
            h_pT_barrel_sub_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_sub_fail[i]->SetMarkerColor(kBlack);
            h_pT_barrel_sub_fail[i]->SetLineColor(kBlack);
            h_pT_endcap_sub_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_sub_fail[i]->SetMarkerColor(kBlack);
            h_pT_endcap_sub_fail[i]->SetLineColor(kBlack);
            h_MET_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_MET_fail[i]->SetMarkerColor(kBlack);
            h_MET_fail[i]->SetLineColor(kBlack);
            h_MT_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_MT_fail[i]->SetMarkerColor(kBlack);
            h_MT_fail[i]->SetLineColor(kBlack);
        }
        else
        {
            if (i == 1) color = kOrange;
            else if (i == 2) color = kCyan + 2;
            else color = kRed + 3;
            h_pT_barrel_lead_pass[i]->SetFillColor(color);
            h_pT_barrel_lead_pass[i]->SetLineColor(color);
            h_pT_endcap_lead_pass[i]->SetFillColor(color);
            h_pT_endcap_lead_pass[i]->SetLineColor(color);
            h_pT_barrel_lead_fail[i]->SetFillColor(color);
            h_pT_barrel_lead_fail[i]->SetLineColor(color);
            h_pT_endcap_lead_fail[i]->SetFillColor(color);
            h_pT_endcap_lead_fail[i]->SetLineColor(color);
            h_pT_barrel_sub_pass[i]->SetFillColor(color);
            h_pT_barrel_sub_pass[i]->SetLineColor(color);
            h_pT_endcap_sub_pass[i]->SetFillColor(color);
            h_pT_endcap_sub_pass[i]->SetLineColor(color);
            h_pT_barrel_sub_fail[i]->SetFillColor(color);
            h_pT_barrel_sub_fail[i]->SetLineColor(color);
            h_pT_endcap_sub_fail[i]->SetFillColor(color);
            h_pT_endcap_sub_fail[i]->SetLineColor(color);
            h_MET_fail[i]->SetFillColor(color);
            h_MET_fail[i]->SetLineColor(color);
            h_MT_fail[i]->SetFillColor(color);
            h_MT_fail[i]->SetLineColor(color);

            s_pT_barrel_lead_pass->Add(h_pT_barrel_lead_pass[i]);
            s_pT_endcap_lead_pass->Add(h_pT_endcap_lead_pass[i]);
            s_pT_barrel_lead_fail->Add(h_pT_barrel_lead_fail[i]);
            s_pT_endcap_lead_fail->Add(h_pT_endcap_lead_fail[i]);
            s_pT_barrel_sub_pass->Add(h_pT_barrel_sub_pass[i]);
            s_pT_endcap_sub_pass->Add(h_pT_endcap_sub_pass[i]);
            s_pT_barrel_sub_fail->Add(h_pT_barrel_sub_fail[i]);
            s_pT_endcap_sub_fail->Add(h_pT_endcap_sub_fail[i]);
            s_MET_fail->Add(h_MET_fail[i]);
            s_MT_fail->Add(h_MT_fail[i]);
        }
    }


    // Creating and drawing ratio plots
    myRatioPlot_t *RP_pT_barrel_lead_pass = new myRatioPlot_t("c_pT_barrel_lead_pass", s_pT_barrel_lead_pass, h_pT_barrel_lead_pass[0]);
    myRatioPlot_t *RP_pT_endcap_lead_pass = new myRatioPlot_t("c_pT_endcap_lead_pass", s_pT_endcap_lead_pass, h_pT_endcap_lead_pass[0]);
    myRatioPlot_t *RP_pT_barrel_lead_fail = new myRatioPlot_t("c_pT_barrel_lead_fail", s_pT_barrel_lead_fail, h_pT_barrel_lead_fail[0]);
    myRatioPlot_t *RP_pT_endcap_lead_fail = new myRatioPlot_t("c_pT_endcap_lead_fail", s_pT_endcap_lead_fail, h_pT_endcap_lead_fail[0]);
    myRatioPlot_t *RP_pT_barrel_sub_pass = new myRatioPlot_t("c_pT_barrel_sub_pass", s_pT_barrel_sub_pass, h_pT_barrel_sub_pass[0]);
    myRatioPlot_t *RP_pT_endcap_sub_pass = new myRatioPlot_t("c_pT_endcap_sub_pass", s_pT_endcap_sub_pass, h_pT_endcap_sub_pass[0]);
    myRatioPlot_t *RP_pT_barrel_sub_fail = new myRatioPlot_t("c_pT_barrel_sub_fail", s_pT_barrel_sub_fail, h_pT_barrel_sub_fail[0]);
    myRatioPlot_t *RP_pT_endcap_sub_fail = new myRatioPlot_t("c_pT_endcap_sub_fail", s_pT_endcap_sub_fail, h_pT_endcap_sub_fail[0]);
    myRatioPlot_t *RP_MET_fail = new myRatioPlot_t("c_MET_fail", s_MET_fail, h_MET_fail[0]);
    myRatioPlot_t *RP_MT_fail = new myRatioPlot_t("c_MT_fail", s_MT_fail, h_MT_fail[0]);

    RP_pT_barrel_lead_pass->SetPlots("p_{#lower[-0.2]{T}}^{lead} [GeV/c]", 0, 1000);
    RP_pT_endcap_lead_pass->SetPlots("p_{#lower[-0.2]{T}}^{lead} [GeV/c]", 0, 1000);
    RP_pT_barrel_lead_fail->SetPlots("p_{#lower[-0.2]{T}}^{lead} [GeV/c]", 0, 1000);
    RP_pT_endcap_lead_fail->SetPlots("p_{#lower[-0.2]{T}}^{lead} [GeV/c]", 0, 1000);
    RP_pT_barrel_sub_pass->SetPlots("p_{#lower[-0.2]{T}}^{sub} [GeV/c]", 0, 1000);
    RP_pT_endcap_sub_pass->SetPlots("p_{#lower[-0.2]{T}}^{sub} [GeV/c]", 0, 1000);
    RP_pT_barrel_sub_fail->SetPlots("p_{#lower[-0.2]{T}}^{sub} [GeV/c]", 0, 1000);
    RP_pT_endcap_sub_fail->SetPlots("p_{#lower[-0.2]{T}}^{sub} [GeV/c]", 0, 1000);
    RP_MET_fail->SetPlots("E_{#lower[-0.2]{T}}^{miss} [GeV]", 0, 100);
    RP_MT_fail->SetPlots("m_{#lower[-0.2]{T}} [GeV/c^{2}]", 0, 200);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_pT_barrel_lead_pass[0], "Data", "pl");
    legend->AddEntry(h_pT_barrel_lead_pass[1], "DY", "f");
    legend->AddEntry(h_pT_barrel_lead_pass[2], "Bkg (real)", "f");
    legend->AddEntry(h_pT_barrel_lead_pass[3], "Bkg (fake)", "f");

    RP_pT_barrel_lead_pass->ImportLegend(legend);
    RP_pT_endcap_lead_pass->ImportLegend(legend);
    RP_pT_barrel_lead_fail->ImportLegend(legend);
    RP_pT_endcap_lead_fail->ImportLegend(legend);
    RP_pT_barrel_sub_pass->ImportLegend(legend);
    RP_pT_endcap_sub_pass->ImportLegend(legend);
    RP_pT_barrel_sub_fail->ImportLegend(legend);
    RP_pT_endcap_sub_fail->ImportLegend(legend);
    RP_MET_fail->ImportLegend(legend);
    RP_MT_fail->ImportLegend(legend);

    RP_pT_barrel_lead_pass->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_lead_pass->Draw(1e-1, 1e7, 1);
    RP_pT_barrel_lead_fail->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_lead_fail->Draw(1e-1, 1e7, 1);
    RP_pT_barrel_sub_pass->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_sub_pass->Draw(1e-1, 1e7, 1);
    RP_pT_barrel_sub_fail->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_sub_fail->Draw(1e-1, 1e7, 1);
    RP_MET_fail->Draw(1e-1, 1e7, 0);
    RP_MT_fail->Draw(1e-1, 1e7, 0);

    f->Close();
    if (!f->IsOpen()) cout << "File " << "/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " <<"/media/sf_DATA/FR/Electron/For_FakeRate_electron_alt.root" << " COULD NOT BE CLOSED!\n" << endl;

} // End of E_HistDrawer_alt()


void E_HistDrawer_alt2()
{
    DYAnalyzer analyzer("Ele23Ele12");
    FileMgr fm;
    THStack *s_pT_barrel_nume = new THStack("s_pT_barrel_nume", "");
    THStack *s_pT_endcap_nume = new THStack("s_pT_endcap_nume", "");
    THStack *s_pT_barrel_ctrl = new THStack("s_pT_barrel_ctrl", "");
    THStack *s_pT_endcap_ctrl = new THStack("s_pT_endcap_ctrl", "");
    THStack *s_eta_nume = new THStack("s_eta_nume", "");
    THStack *s_eta_ctrl = new THStack("s_eta_ctrl", "");
    THStack *s_MET = new THStack("s_MET", "");
    THStack *s_MT = new THStack("s_MT", "");
    THStack *s_nVTX = new THStack("s_nVTX", "");
    THStack *s_mass = new THStack("s_mass", "");

    TH1D *h_pT_barrel_nume[_EndOf_Data_Special], *h_pT_endcap_nume[_EndOf_Data_Special],
         *h_pT_barrel_ctrl[_EndOf_Data_Special], *h_pT_endcap_ctrl[_EndOf_Data_Special],
         *h_eta_nume[_EndOf_Data_Special], *h_eta_ctrl[_EndOf_Data_Special],
         *h_MET[_EndOf_Data_Special], *h_MT[_EndOf_Data_Special],
         *h_nVTX[_EndOf_Data_Special], *h_mass[_EndOf_Data_Special];

//----------------------------------- MC bkg -------------------------------------------------------

    // Other MC
    Int_t stop = 0;
    Process_t pr = _WW;
    while (!stop)
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_MET", h_MET[pr]);
        file->GetObject("h_MT", h_MT[pr]);
        file->GetObject("h_nVTX", h_nVTX[pr]);
        file->GetObject("h_mass", h_mass[pr]);
        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_MET[pr]);
        removeNegativeBins(h_MT[pr]);
        removeNegativeBins(h_nVTX[pr]);
        removeNegativeBins(h_mass[pr]);

        Color_t color = kBlack;
        if (pr == _WJets || pr == _WJets_ext2v5) color = kRed - 2;
        if (pr == _VVnST) color = kMagenta - 5;
        if (pr == _WW) color = kMagenta - 5;
        if (pr == _WZ) color = kMagenta - 2;
        if (pr == _ZZ) color = kMagenta - 6;
        if (pr == _tbarW) color = kGreen - 2;
        if (pr == _tW) color = kGreen + 2;
        if (pr == _ttbar || pr == _ttbar_700to1000 || pr == _ttbar_1000toInf) color = kCyan + 2;

        h_pT_barrel_nume[pr]->SetFillColor(color);
        h_pT_endcap_nume[pr]->SetFillColor(color);
        h_pT_barrel_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_ctrl[pr]->SetFillColor(color);
        h_eta_nume[pr]      ->SetFillColor(color);
        h_eta_ctrl[pr]      ->SetFillColor(color);
        h_MET[pr]           ->SetFillColor(color);
        h_MT[pr]            ->SetFillColor(color);
        h_nVTX[pr]          ->SetFillColor(color);
        h_mass[pr]          ->SetFillColor(color);

        h_pT_barrel_nume[pr]->SetLineColor(color);
        h_pT_endcap_nume[pr]->SetLineColor(color);
        h_pT_barrel_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_ctrl[pr]->SetLineColor(color);
        h_eta_nume[pr]      ->SetLineColor(color);
        h_eta_ctrl[pr]      ->SetLineColor(color);
        h_MET[pr]           ->SetLineColor(color);
        h_MT[pr]            ->SetLineColor(color);
        h_nVTX[pr]          ->SetLineColor(color);
        h_mass[pr]          ->SetLineColor(color);

        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]      ->SetDirectory(0);
        h_eta_ctrl[pr]      ->SetDirectory(0);
        h_MET[pr]           ->SetDirectory(0);
        h_MT[pr]            ->SetDirectory(0);
        h_nVTX[pr]          ->SetDirectory(0);
        h_mass[pr]          ->SetDirectory(0);

        if (pr == _WJets)
        {
            h_pT_barrel_nume[_WJets_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_WJets")));
            h_pT_endcap_nume[_WJets_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_WJets")));
            h_pT_barrel_ctrl[_WJets_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_WJets")));
            h_pT_endcap_ctrl[_WJets_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_WJets")));
            h_eta_nume[_WJets_Full]       = ((TH1D*)(h_eta_nume[pr]      ->Clone("h_eta_nume_WJets")));
            h_eta_ctrl[_WJets_Full]       = ((TH1D*)(h_eta_ctrl[pr]      ->Clone("h_eta_ctrl_WJets")));
            h_MET[_WJets_Full]            = ((TH1D*)(h_MET[pr]           ->Clone("h_MET_WJets")));
            h_MT[_WJets_Full]             = ((TH1D*)(h_MT[pr]            ->Clone("h_MT_WJets")));
            h_nVTX[_WJets_Full]           = ((TH1D*)(h_nVTX[pr]          ->Clone("h_nVTX_WJets")));
            h_mass[_WJets_Full]           = ((TH1D*)(h_mass[pr]          ->Clone("h_mass_WJets")));

            h_pT_barrel_nume[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_nume[_WJets_Full]->SetDirectory(0);
            h_pT_barrel_ctrl[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_WJets_Full]->SetDirectory(0);
            h_eta_nume[_WJets_Full]      ->SetDirectory(0);
            h_eta_ctrl[_WJets_Full]      ->SetDirectory(0);
            h_MET[_WJets_Full]           ->SetDirectory(0);
            h_MT[_WJets_Full]            ->SetDirectory(0);
            h_nVTX[_WJets_Full]          ->SetDirectory(0);
            h_mass[_WJets_Full]          ->SetDirectory(0);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_pT_barrel_nume[_WJets_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_WJets_Full]->Add(h_pT_endcap_nume[pr]);
            h_pT_barrel_ctrl[_WJets_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_WJets_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_eta_nume[_WJets_Full]      ->Add(h_eta_nume[pr]);
            h_eta_ctrl[_WJets_Full]      ->Add(h_eta_ctrl[pr]);
            h_MET[_WJets_Full]           ->Add(h_MET[pr]);
            h_MT[_WJets_Full]            ->Add(h_MT[pr]);
            h_nVTX[_WJets_Full]          ->Add(h_nVTX[pr]);
            h_mass[_WJets_Full]          ->Add(h_mass[pr]);
        }

        s_pT_barrel_nume->Add(h_pT_barrel_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_nume[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_ctrl[pr]);
        s_eta_nume      ->Add(h_eta_nume[pr]);
        s_eta_ctrl      ->Add(h_eta_ctrl[pr]);
        s_MET           ->Add(h_MET[pr]);
        s_MT            ->Add(h_MT[pr]);
        s_nVTX          ->Add(h_nVTX[pr]);
        s_mass          ->Add(h_mass[pr]);

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

    // Drell-Yan
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_MET", h_MET[pr]);
        file->GetObject("h_MT", h_MT[pr]);
        file->GetObject("h_nVTX", h_nVTX[pr]);
        file->GetObject("h_mass", h_mass[pr]);

        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_MET[pr]);
        removeNegativeBins(h_MT[pr]);
        removeNegativeBins(h_nVTX[pr]);
        removeNegativeBins(h_mass[pr]);

        Color_t color = kOrange - 5;
        h_pT_barrel_nume[pr]->SetFillColor(color);
        h_pT_endcap_nume[pr]->SetFillColor(color);
        h_pT_barrel_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_ctrl[pr]->SetFillColor(color);
        h_eta_nume[pr]      ->SetFillColor(color);
        h_eta_ctrl[pr]      ->SetFillColor(color);
        h_MET[pr]           ->SetFillColor(color);
        h_MT[pr]            ->SetFillColor(color);
        h_nVTX[pr]          ->SetFillColor(color);
        h_mass[pr]          ->SetFillColor(color);

        h_pT_barrel_nume[pr]->SetLineColor(color);
        h_pT_endcap_nume[pr]->SetLineColor(color);
        h_pT_barrel_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_ctrl[pr]->SetLineColor(color);
        h_eta_nume[pr]      ->SetLineColor(color);
        h_eta_ctrl[pr]      ->SetLineColor(color);
        h_MET[pr]           ->SetLineColor(color);
        h_MT[pr]            ->SetLineColor(color);
        h_nVTX[pr]          ->SetLineColor(color);
        h_mass[pr]          ->SetLineColor(color);

        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]      ->SetDirectory(0);
        h_eta_ctrl[pr]      ->SetDirectory(0);
        h_MET[pr]           ->SetDirectory(0);
        h_MT[pr]            ->SetDirectory(0);
        h_nVTX[pr]          ->SetDirectory(0);
        h_mass[pr]          ->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_pT_barrel_nume[_DY_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_DY")));
            h_pT_endcap_nume[_DY_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_DY")));
            h_pT_barrel_ctrl[_DY_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_DY")));
            h_pT_endcap_ctrl[_DY_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_DY")));
            h_eta_nume[_DY_Full]       = ((TH1D*)(h_eta_nume[pr]      ->Clone("h_eta_nume_DY")));
            h_eta_ctrl[_DY_Full]       = ((TH1D*)(h_eta_ctrl[pr]      ->Clone("h_eta_ctrl_DY")));
            h_MET[_DY_Full]            = ((TH1D*)(h_MET[pr]           ->Clone("h_MET_DY")));
            h_MT[_DY_Full]             = ((TH1D*)(h_MT[pr]            ->Clone("h_MT_DY")));
            h_nVTX[_DY_Full]           = ((TH1D*)(h_nVTX[pr]          ->Clone("h_nVTX_DY")));
            h_mass[_DY_Full]           = ((TH1D*)(h_mass[pr]          ->Clone("h_mass_DY")));

            h_pT_barrel_nume[_DY_Full]->SetDirectory(0);
            h_pT_endcap_nume[_DY_Full]->SetDirectory(0);
            h_pT_barrel_ctrl[_DY_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_DY_Full]->SetDirectory(0);
            h_eta_nume[_DY_Full]      ->SetDirectory(0);
            h_eta_ctrl[_DY_Full]      ->SetDirectory(0);
            h_MET[_DY_Full]           ->SetDirectory(0);
            h_MT[_DY_Full]            ->SetDirectory(0);
            h_nVTX[_DY_Full]          ->SetDirectory(0);
            h_mass[_DY_Full]          ->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_nume[_DY_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_DY_Full]->Add(h_pT_endcap_nume[pr]);
            h_pT_barrel_ctrl[_DY_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_DY_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_eta_nume[_DY_Full]      ->Add(h_eta_nume[pr]);
            h_eta_ctrl[_DY_Full]      ->Add(h_eta_ctrl[pr]);
            h_MET[_DY_Full]           ->Add(h_MET[pr]);
            h_MT[_DY_Full]            ->Add(h_MT[pr]);
            h_nVTX[_DY_Full]          ->Add(h_nVTX[pr]);
            h_mass[_DY_Full]          ->Add(h_mass[pr]);
        }

        s_pT_barrel_nume->Add(h_pT_barrel_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_nume[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_ctrl[pr]);
        s_eta_nume      ->Add(h_eta_nume[pr]);
        s_eta_ctrl      ->Add(h_eta_ctrl[pr]);
        s_MET           ->Add(h_MET[pr]);
        s_MT            ->Add(h_MT[pr]);
        s_nVTX          ->Add(h_nVTX[pr]);
        s_mass          ->Add(h_mass[pr]);

        file->Close();
    }

    // QCD
    for (Process_t pr = _QCDEMEnriched_20to30; pr <= _QCDEMEnriched_300toInf; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_MET", h_MET[pr]);
        file->GetObject("h_MT", h_MT[pr]);
        file->GetObject("h_nVTX", h_nVTX[pr]);
        file->GetObject("h_mass", h_mass[pr]);

        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_MET[pr]);
        removeNegativeBins(h_MT[pr]);
        removeNegativeBins(h_nVTX[pr]);
        removeNegativeBins(h_mass[pr]);

        Color_t color = kRed + 3;
        h_pT_barrel_nume[pr]->SetFillColor(color);
        h_pT_endcap_nume[pr]->SetFillColor(color);
        h_pT_barrel_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_ctrl[pr]->SetFillColor(color);
        h_eta_nume[pr]      ->SetFillColor(color);
        h_eta_ctrl[pr]      ->SetFillColor(color);
        h_MET[pr]           ->SetFillColor(color);
        h_MT[pr]            ->SetFillColor(color);
        h_nVTX[pr]          ->SetFillColor(color);
        h_mass[pr]          ->SetFillColor(color);

        h_pT_barrel_nume[pr]->SetLineColor(color);
        h_pT_endcap_nume[pr]->SetLineColor(color);
        h_pT_barrel_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_ctrl[pr]->SetLineColor(color);
        h_eta_nume[pr]      ->SetLineColor(color);
        h_eta_ctrl[pr]      ->SetLineColor(color);
        h_MET[pr]           ->SetLineColor(color);
        h_MT[pr]            ->SetLineColor(color);
        h_nVTX[pr]          ->SetLineColor(color);
        h_mass[pr]          ->SetLineColor(color);

        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]      ->SetDirectory(0);
        h_eta_ctrl[pr]      ->SetDirectory(0);
        h_MET[pr]           ->SetDirectory(0);
        h_MT[pr]            ->SetDirectory(0);
        h_nVTX[pr]          ->SetDirectory(0);
        h_mass[pr]          ->SetDirectory(0);

        if (pr == _QCDEMEnriched_20to30)
        {
            h_pT_barrel_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_QCD")));
            h_pT_endcap_nume[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_QCD")));
            h_pT_barrel_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_QCD")));
            h_pT_endcap_ctrl[_QCDEMEnriched_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_QCD")));
            h_eta_nume[_QCDEMEnriched_Full]       = ((TH1D*)(h_eta_nume[pr]      ->Clone("h_eta_nume")));
            h_eta_ctrl[_QCDEMEnriched_Full]       = ((TH1D*)(h_eta_ctrl[pr]      ->Clone("h_eta_ctrl")));
            h_MET[_QCDEMEnriched_Full]            = ((TH1D*)(h_MET[pr]           ->Clone("h_MET_QCD")));
            h_MT[_QCDEMEnriched_Full]             = ((TH1D*)(h_MT[pr]            ->Clone("h_MT_QCD")));
            h_nVTX[_QCDEMEnriched_Full]           = ((TH1D*)(h_nVTX[pr]          ->Clone("h_nVTX_QCD")));
            h_mass[_QCDEMEnriched_Full]           = ((TH1D*)(h_mass[pr]          ->Clone("h_mass_QCD")));

            h_pT_barrel_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_nume[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_barrel_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_QCDEMEnriched_Full]->SetDirectory(0);
            h_eta_nume[_QCDEMEnriched_Full]      ->SetDirectory(0);
            h_eta_ctrl[_QCDEMEnriched_Full]      ->SetDirectory(0);
            h_MET[_QCDEMEnriched_Full]           ->SetDirectory(0);
            h_MT[_QCDEMEnriched_Full]            ->SetDirectory(0);
            h_nVTX[_QCDEMEnriched_Full]          ->SetDirectory(0);
            h_mass[_QCDEMEnriched_Full]          ->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_nume[_QCDEMEnriched_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_QCDEMEnriched_Full]->Add(h_pT_endcap_nume[pr]);
            h_pT_barrel_ctrl[_QCDEMEnriched_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_QCDEMEnriched_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_eta_nume[_QCDEMEnriched_Full]      ->Add(h_eta_nume[pr]);
            h_eta_ctrl[_QCDEMEnriched_Full]      ->Add(h_eta_ctrl[pr]);
            h_MET[_QCDEMEnriched_Full]           ->Add(h_MET[pr]);
            h_MT[_QCDEMEnriched_Full]            ->Add(h_MT[pr]);
            h_nVTX[_QCDEMEnriched_Full]          ->Add(h_nVTX[pr]);
            h_mass[_QCDEMEnriched_Full]          ->Add(h_mass[pr]);
        }

        s_pT_barrel_nume->Add(h_pT_barrel_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_nume[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_ctrl[pr]);
        s_eta_nume      ->Add(h_eta_nume[pr]);
        s_eta_ctrl      ->Add(h_eta_ctrl[pr]);
        s_MET           ->Add(h_MET[pr]);
        s_MT            ->Add(h_MT[pr]);
        s_nVTX          ->Add(h_nVTX[pr]);
        s_mass          ->Add(h_mass[pr]);

        file->Close();
    }

    // GammaJets
    for (Process_t pr = _GJets_20to100; pr <= _GJets_2000to5000; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_MET", h_MET[pr]);
        file->GetObject("h_MT", h_MT[pr]);
        file->GetObject("h_nVTX", h_nVTX[pr]);
        file->GetObject("h_mass", h_mass[pr]);

        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_MET[pr]);
        removeNegativeBins(h_MT[pr]);
        removeNegativeBins(h_nVTX[pr]);
        removeNegativeBins(h_mass[pr]);

        Color_t color = kYellow + 3;
        h_pT_barrel_nume[pr]->SetFillColor(color);
        h_pT_endcap_nume[pr]->SetFillColor(color);
        h_pT_barrel_ctrl[pr]->SetFillColor(color);
        h_pT_endcap_ctrl[pr]->SetFillColor(color);
        h_eta_nume[pr]      ->SetFillColor(color);
        h_eta_ctrl[pr]      ->SetFillColor(color);
        h_MET[pr]           ->SetFillColor(color);
        h_MT[pr]            ->SetFillColor(color);
        h_nVTX[pr]          ->SetFillColor(color);
        h_mass[pr]          ->SetFillColor(color);

        h_pT_barrel_nume[pr]->SetLineColor(color);
        h_pT_endcap_nume[pr]->SetLineColor(color);
        h_pT_barrel_ctrl[pr]->SetLineColor(color);
        h_pT_endcap_ctrl[pr]->SetLineColor(color);
        h_eta_nume[pr]      ->SetLineColor(color);
        h_eta_ctrl[pr]      ->SetLineColor(color);
        h_MET[pr]           ->SetLineColor(color);
        h_MT[pr]            ->SetLineColor(color);
        h_nVTX[pr]          ->SetLineColor(color);
        h_mass[pr]          ->SetLineColor(color);

        h_pT_barrel_nume[pr]->SetDirectory(0);
        h_pT_endcap_nume[pr]->SetDirectory(0);
        h_pT_barrel_ctrl[pr]->SetDirectory(0);
        h_pT_endcap_ctrl[pr]->SetDirectory(0);
        h_eta_nume[pr]      ->SetDirectory(0);
        h_eta_ctrl[pr]      ->SetDirectory(0);
        h_MET[pr]           ->SetDirectory(0);
        h_MT[pr]            ->SetDirectory(0);
        h_nVTX[pr]          ->SetDirectory(0);
        h_mass[pr]          ->SetDirectory(0);

        if (pr == _GJets_20to100)
        {
            h_pT_barrel_nume[_GJets_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_GJets")));
            h_pT_endcap_nume[_GJets_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_GJets")));
            h_pT_barrel_ctrl[_GJets_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_GJets")));
            h_pT_endcap_ctrl[_GJets_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_GJets")));
            h_eta_nume[_GJets_Full]       = ((TH1D*)(h_eta_nume[pr]      ->Clone("h_eta_nume_GJets")));
            h_eta_ctrl[_GJets_Full]       = ((TH1D*)(h_eta_ctrl[pr]      ->Clone("h_eta_ctrl_GJets")));
            h_MET[_GJets_Full]            = ((TH1D*)(h_MET[pr]           ->Clone("h_MET_GJets")));
            h_MT[_GJets_Full]             = ((TH1D*)(h_MT[pr]            ->Clone("h_MT_GJets")));
            h_nVTX[_GJets_Full]           = ((TH1D*)(h_nVTX[pr]          ->Clone("h_nVTX_GJets")));
            h_mass[_GJets_Full]           = ((TH1D*)(h_mass[pr]          ->Clone("h_mass_GJets")));

            h_pT_barrel_nume[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_nume[_GJets_Full]->SetDirectory(0);
            h_pT_barrel_ctrl[_GJets_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_GJets_Full]->SetDirectory(0);
            h_eta_nume[_GJets_Full]      ->SetDirectory(0);
            h_eta_ctrl[_GJets_Full]      ->SetDirectory(0);
            h_MET[_GJets_Full]           ->SetDirectory(0);
            h_MT[_GJets_Full]            ->SetDirectory(0);
            h_nVTX[_GJets_Full]          ->SetDirectory(0);
            h_mass[_GJets_Full]          ->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_nume[_GJets_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_GJets_Full]->Add(h_pT_endcap_nume[pr]);
            h_pT_barrel_ctrl[_GJets_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_GJets_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_eta_nume[_GJets_Full]      ->Add(h_eta_nume[pr]);
            h_eta_ctrl[_GJets_Full]      ->Add(h_eta_ctrl[pr]);
            h_MET[_GJets_Full]           ->Add(h_MET[pr]);
            h_MT[_GJets_Full]            ->Add(h_MT[pr]);
            h_nVTX[_GJets_Full]          ->Add(h_nVTX[pr]);
            h_mass[_GJets_Full]          ->Add(h_mass[pr]);
        }

        s_pT_barrel_nume->Add(h_pT_barrel_nume[pr]);
        s_pT_endcap_nume->Add(h_pT_endcap_nume[pr]);
        s_pT_barrel_ctrl->Add(h_pT_barrel_ctrl[pr]);
        s_pT_endcap_ctrl->Add(h_pT_endcap_ctrl[pr]);
        s_eta_nume      ->Add(h_eta_nume[pr]);
        s_eta_ctrl      ->Add(h_eta_ctrl[pr]);
        s_MET           ->Add(h_MET[pr]);
        s_MT            ->Add(h_MT[pr]);
        s_nVTX          ->Add(h_nVTX[pr]);
        s_mass          ->Add(h_mass[pr]);

        file->Close();
    }

//--------------------------------------- DATA -----------------------------------------------------

    for (Process_t pr=_DoubleEG_B; pr<=_DoubleEG_H; pr=next(pr))
    {
        TFile *file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_alt_"+fm.Procname[pr]+".root", "READ");

        file->GetObject("h_pT_barrel_nume", h_pT_barrel_nume[pr]);
        file->GetObject("h_pT_endcap_nume", h_pT_endcap_nume[pr]);
        file->GetObject("h_pT_barrel_ctrl", h_pT_barrel_ctrl[pr]);
        file->GetObject("h_pT_endcap_ctrl", h_pT_endcap_ctrl[pr]);
        file->GetObject("h_eta_nume", h_eta_nume[pr]);
        file->GetObject("h_eta_ctrl", h_eta_ctrl[pr]);
        file->GetObject("h_MET", h_MET[pr]);
        file->GetObject("h_MT", h_MT[pr]);
        file->GetObject("h_nVTX", h_nVTX[pr]);
        file->GetObject("h_mass", h_mass[pr]);

        removeNegativeBins(h_pT_barrel_nume[pr]);
        removeNegativeBins(h_pT_endcap_nume[pr]);
        removeNegativeBins(h_pT_barrel_ctrl[pr]);
        removeNegativeBins(h_pT_endcap_ctrl[pr]);
        removeNegativeBins(h_eta_nume[pr]);
        removeNegativeBins(h_eta_ctrl[pr]);
        removeNegativeBins(h_MET[pr]);
        removeNegativeBins(h_MT[pr]);
        removeNegativeBins(h_nVTX[pr]);
        removeNegativeBins(h_mass[pr]);

        if (pr == _DoubleEG_B)
        {
            h_pT_barrel_nume[_DoubleEG_Full] = ((TH1D*)(h_pT_barrel_nume[pr]->Clone("h_pT_barrel_nume_data")));
            h_pT_endcap_nume[_DoubleEG_Full] = ((TH1D*)(h_pT_endcap_nume[pr]->Clone("h_pT_endcap_nume_data")));
            h_pT_barrel_ctrl[_DoubleEG_Full] = ((TH1D*)(h_pT_barrel_ctrl[pr]->Clone("h_pT_barrel_ctrl_data")));
            h_pT_endcap_ctrl[_DoubleEG_Full] = ((TH1D*)(h_pT_endcap_ctrl[pr]->Clone("h_pT_endcap_ctrl_data")));
            h_eta_nume[_DoubleEG_Full]       = ((TH1D*)(h_eta_nume[pr]      ->Clone("h_eta_nume_data")));
            h_eta_ctrl[_DoubleEG_Full]       = ((TH1D*)(h_eta_ctrl[pr]      ->Clone("h_eta_ctrl_data")));
            h_MET[_DoubleEG_Full]            = ((TH1D*)(h_MET[pr]           ->Clone("h_MET_data")));
            h_MT[_DoubleEG_Full]             = ((TH1D*)(h_MT[pr]            ->Clone("h_MT_data")));
            h_nVTX[_DoubleEG_Full]           = ((TH1D*)(h_nVTX[pr]          ->Clone("h_nVTX_data")));
            h_mass[_DoubleEG_Full]           = ((TH1D*)(h_mass[pr]          ->Clone("h_mass_data")));

            h_pT_barrel_nume[_DoubleEG_Full]->SetDirectory(0);
            h_pT_endcap_nume[_DoubleEG_Full]->SetDirectory(0);
            h_pT_barrel_ctrl[_DoubleEG_Full]->SetDirectory(0);
            h_pT_endcap_ctrl[_DoubleEG_Full]->SetDirectory(0);
            h_eta_nume[_DoubleEG_Full]      ->SetDirectory(0);
            h_eta_ctrl[_DoubleEG_Full]      ->SetDirectory(0);
            h_MET[_DoubleEG_Full]           ->SetDirectory(0);
            h_MT[_DoubleEG_Full]            ->SetDirectory(0);
            h_nVTX[_DoubleEG_Full]          ->SetDirectory(0);
            h_mass[_DoubleEG_Full]          ->SetDirectory(0);
        }
        else
        {
            h_pT_barrel_nume[_DoubleEG_Full]->Add(h_pT_barrel_nume[pr]);
            h_pT_endcap_nume[_DoubleEG_Full]->Add(h_pT_endcap_nume[pr]);
            h_pT_barrel_ctrl[_DoubleEG_Full]->Add(h_pT_barrel_ctrl[pr]);
            h_pT_endcap_ctrl[_DoubleEG_Full]->Add(h_pT_endcap_ctrl[pr]);
            h_eta_nume[_DoubleEG_Full]      ->Add(h_eta_nume[pr]);
            h_eta_ctrl[_DoubleEG_Full]      ->Add(h_eta_ctrl[pr]);
            h_MET[_DoubleEG_Full]           ->Add(h_MET[pr]);
            h_MT[_DoubleEG_Full]            ->Add(h_MT[pr]);
            h_nVTX[_DoubleEG_Full]          ->Add(h_nVTX[pr]);
            h_mass[_DoubleEG_Full]          ->Add(h_mass[pr]);
        }

    }

    h_pT_barrel_nume[_DoubleEG_Full]->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_nume[_DoubleEG_Full]->SetMarkerStyle(kFullDotLarge);
    h_pT_barrel_ctrl[_DoubleEG_Full]->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_ctrl[_DoubleEG_Full]->SetMarkerStyle(kFullDotLarge);
    h_eta_nume[_DoubleEG_Full]      ->SetMarkerStyle(kFullDotLarge);
    h_eta_ctrl[_DoubleEG_Full]      ->SetMarkerStyle(kFullDotLarge);
    h_MET[_DoubleEG_Full]           ->SetMarkerStyle(kFullDotLarge);
    h_MT[_DoubleEG_Full]            ->SetMarkerStyle(kFullDotLarge);
    h_nVTX[_DoubleEG_Full]          ->SetMarkerStyle(kFullDotLarge);
    h_mass[_DoubleEG_Full]          ->SetMarkerStyle(kFullDotLarge);

    h_pT_barrel_nume[_DoubleEG_Full]->SetMarkerColor(kBlack);
    h_pT_endcap_nume[_DoubleEG_Full]->SetMarkerColor(kBlack);
    h_pT_barrel_ctrl[_DoubleEG_Full]->SetMarkerColor(kBlack);
    h_pT_endcap_ctrl[_DoubleEG_Full]->SetMarkerColor(kBlack);
    h_eta_nume[_DoubleEG_Full]      ->SetMarkerColor(kBlack);
    h_eta_ctrl[_DoubleEG_Full]      ->SetMarkerColor(kBlack);
    h_MET[_DoubleEG_Full]           ->SetMarkerColor(kBlack);
    h_MT[_DoubleEG_Full]            ->SetMarkerColor(kBlack);
    h_nVTX[_DoubleEG_Full]          ->SetMarkerColor(kBlack);
    h_mass[_DoubleEG_Full]          ->SetMarkerColor(kBlack);

    h_pT_barrel_nume[_DoubleEG_Full]->SetLineColor(kBlack);
    h_pT_endcap_nume[_DoubleEG_Full]->SetLineColor(kBlack);
    h_pT_barrel_ctrl[_DoubleEG_Full]->SetLineColor(kBlack);
    h_pT_endcap_ctrl[_DoubleEG_Full]->SetLineColor(kBlack);
    h_eta_nume[_DoubleEG_Full]      ->SetLineColor(kBlack);
    h_eta_ctrl[_DoubleEG_Full]      ->SetLineColor(kBlack);
    h_MET[_DoubleEG_Full]           ->SetLineColor(kBlack);
    h_MT[_DoubleEG_Full]            ->SetLineColor(kBlack);
    h_nVTX[_DoubleEG_Full]          ->SetLineColor(kBlack);
    h_mass[_DoubleEG_Full]          ->SetLineColor(kBlack);

//--------------------------------- Ratio Plots --------------------------------------

    myRatioPlot_t *RP_pT_barrel_nume = new myRatioPlot_t("RP_pT_barrel_nume", s_pT_barrel_nume, h_pT_barrel_nume[_DoubleEG_Full]);
    myRatioPlot_t *RP_pT_endcap_nume = new myRatioPlot_t("RP_pT_endcap_nume", s_pT_endcap_nume, h_pT_endcap_nume[_DoubleEG_Full]);
    myRatioPlot_t *RP_pT_barrel_ctrl = new myRatioPlot_t("RP_pT_barrel_ctrl", s_pT_barrel_ctrl, h_pT_barrel_ctrl[_DoubleEG_Full]);
    myRatioPlot_t *RP_pT_endcap_ctrl = new myRatioPlot_t("RP_pT_endcap_ctrl", s_pT_endcap_ctrl, h_pT_endcap_ctrl[_DoubleEG_Full]);
    myRatioPlot_t *RP_eta_nume       = new myRatioPlot_t("RP_eta_nume",       s_eta_nume,       h_eta_nume[_DoubleEG_Full]);
    myRatioPlot_t *RP_eta_ctrl       = new myRatioPlot_t("RP_eta_ctrl",       s_eta_ctrl,       h_eta_ctrl[_DoubleEG_Full]);
    myRatioPlot_t *RP_MET            = new myRatioPlot_t("RP_MET",            s_MET,            h_MET[_DoubleEG_Full]);
    myRatioPlot_t *RP_MT             = new myRatioPlot_t("RP_MT",             s_MT,             h_MT[_DoubleEG_Full]);
    myRatioPlot_t *RP_nVTX           = new myRatioPlot_t("RP_nVTX",           s_nVTX,           h_nVTX[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass           = new myRatioPlot_t("RP_mass",           s_mass,           h_mass[_DoubleEG_Full]);

    RP_pT_barrel_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{signal}) [GeV/c]", 15, 10000);
    RP_pT_endcap_nume->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{signal}) [GeV/c]", 15, 10000);
    RP_pT_barrel_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{barrel}}^{non-signal}) [GeV/c]", 15, 10000);
    RP_pT_endcap_ctrl->SetPlots("p_{#lower[-0.25]{T}} (e_{#lower[-0.4]{endcap}}^{non-signal}) [GeV/c]", 15, 10000);
    RP_MET           ->SetPlots("E_{#lower[-0.25]{T}}^{miss} [GeV]", 0, 500);
    RP_MT            ->SetPlots("m_{#lower[-0.25]{T}} (e_{#lower[-0.4]{lead}}, E_{#lower[-0.25]{T}}^{miss}) [GeV/c^{2}]", 0, 500);
    RP_eta_nume      ->SetPlots("#eta (e_{#lower[-0.4]{signal}})", -3, 3);
    RP_eta_ctrl      ->SetPlots("#eta (e_{#lower[-0.4]{non-signal}})", -3, 3);
    RP_nVTX          ->SetPlots("N_{#lower[-0.25]{VTX}}", 0, 50);
    RP_mass          ->SetPlots("m_{#lower[-0.25]{ee}} [GeV/c^{2}]", 15, 3000);

    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_pT_barrel_nume[_DoubleEG_Full], "Data", "lp");
    legend->AddEntry(h_pT_barrel_nume[_DY_50to100], "DY", "f");
    legend->AddEntry(h_pT_barrel_nume[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend->AddEntry(h_pT_barrel_nume[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_pT_barrel_nume[_GJets_20to100], "#gamma+Jets", "f");
    legend->SetNColumns(2);

    RP_pT_barrel_nume->ImportLegend(legend);
    RP_pT_endcap_nume->ImportLegend(legend);
    RP_pT_barrel_ctrl->ImportLegend(legend);
    RP_pT_endcap_ctrl->ImportLegend(legend);
    RP_MET           ->ImportLegend(legend);
    RP_MT            ->ImportLegend(legend);
    RP_eta_nume      ->ImportLegend(legend);
    RP_eta_ctrl      ->ImportLegend(legend);
    RP_nVTX          ->ImportLegend(legend);
    RP_mass          ->ImportLegend(legend);

    RP_pT_barrel_nume->Draw(1e-1, 1e4, 1);
    RP_pT_endcap_nume->Draw(1e-1, 1e4, 1);
    RP_pT_barrel_ctrl->Draw(1e-1, 1e4, 1);
    RP_pT_endcap_ctrl->Draw(1e-1, 1e4, 1);
    RP_MET           ->Draw(1e-1, 1e4, 0);
    RP_MT            ->Draw(1e-1, 1e4, 0);
    RP_eta_nume      ->Draw(1e-1, 1e4, 0);
    RP_eta_ctrl      ->Draw(1e-1, 1e4, 0);
    RP_nVTX          ->Draw(1e-1, 1e4, 0);
    RP_mass          ->Draw(1e-1, 1e4, 1);

    cout << "MC integral (ctrl): " << ((TH1D*)(s_pT_barrel_ctrl->GetStack()->Last()))->Integral() +
                                     ((TH1D*)(s_pT_endcap_ctrl->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (ctrl): " << h_pT_barrel_ctrl[_DoubleEG_Full]->Integral() + h_pT_endcap_ctrl[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (ctrl): " << h_pT_barrel_ctrl[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_ctrl[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (ctrl): " << h_pT_barrel_ctrl[_GJets_Full]->Integral()+h_pT_endcap_ctrl[_GJets_Full]->Integral() << endl;
    cout << "DY integral (ctrl): " << h_pT_barrel_ctrl[_DY_Full]->Integral()+h_pT_endcap_ctrl[_DY_Full]->Integral() << endl;
    cout << "MC integral(nume): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() +
            ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume): " << h_pT_barrel_nume[_DoubleEG_Full]->Integral() + h_pT_endcap_nume[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (nume): " << h_pT_barrel_nume[_QCDEMEnriched_Full]->Integral()+h_pT_endcap_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume): " << h_pT_barrel_nume[_GJets_Full]->Integral()+h_pT_endcap_nume[_GJets_Full]->Integral() << endl;
    cout << "DY integral (nume): " << h_pT_barrel_nume[_DY_Full]->Integral()+h_pT_endcap_nume[_DY_Full]->Integral() << endl;

    cout << "\nMC integral (ctrl barrel): " << ((TH1D*)(s_pT_barrel_ctrl->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (ctrl barrel): " << h_pT_barrel_ctrl[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (ctrl barrel): " << h_pT_barrel_ctrl[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (ctrl barrel): " << h_pT_barrel_ctrl[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (nume barrel): " << ((TH1D*)(s_pT_barrel_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume barrel): " << h_pT_barrel_nume[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (nume barrel): " << h_pT_barrel_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume barrel): " << h_pT_barrel_nume[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (ctrl endcap): " << ((TH1D*)(s_pT_endcap_ctrl->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (ctrl endcap): " << h_pT_endcap_ctrl[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (ctrl endcap): " << h_pT_endcap_ctrl[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (ctrl endcap): " << h_pT_endcap_ctrl[_GJets_Full]->Integral() << endl;
    cout << "\nMC integral (nume endcap): " << ((TH1D*)(s_pT_endcap_nume->GetStack()->Last()))->Integral() << endl;
    cout << "Data integral (nume endcap): " << h_pT_endcap_nume[_DoubleEG_Full]->Integral() << endl;
    cout << "QCD integral (nume endcap): " << h_pT_endcap_nume[_QCDEMEnriched_Full]->Integral() << endl;
    cout << "Gamma+Jets integral (nume endcap): " << h_pT_endcap_nume[_GJets_Full]->Integral() << endl;

} // End of EE_HistDrawer_alt2()


/// ################################################################################## ///
/// -------------------------------- Muon Channel ------------------------------------ ///
/// ################################################################################## ///
void Mu_HistDrawer(Int_t type)
{
    FileMgr fm;
    DYAnalyzer analyzer("None");
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
    THStack *s_pT_barrel_deno_density = new THStack("s_pT_barrel_deno_density", "");
    THStack *s_pT_endcap_deno_density = new THStack("s_pT_endcap_deno_density", "");

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
         *h_pT_barrel_MC_deno_density[_EndOf_Data_Special], *h_pT_endcap_MC_deno_density[_EndOf_Data_Special],
         *h_PFiso_barrel_data_nume, *h_PFiso_endcap_data_nume, *h_PFiso_barrel_data_deno, *h_PFiso_endcap_data_deno,
         *h_PFiso_barrel_data_ctrl, *h_PFiso_endcap_data_ctrl, *h_pT_barrel_data_nume, *h_pT_endcap_data_nume,
         *h_pT_barrel_data_deno, *h_pT_endcap_data_deno, *h_pT_barrel_data_ctrl, *h_pT_endcap_data_ctrl,
         *h_MET_data, *h_MT_barrel_data_nume, *h_MT_endcap_data_nume,
         *h_MT_barrel_data_deno, *h_MT_endcap_data_deno, *h_MT_barrel_data_ctrl, *h_MT_endcap_data_ctrl,
         *h_pT_barrel_data_deno_density, *h_pT_endcap_data_deno_density;

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
        file->GetObject("h_nVTX", h_nVTX_MC[pr1]);
        h_pT_barrel_MC_deno_density[pr1] = ((TH1D*)(h_pT_barrel_MC_deno[pr1]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr1])));
        h_pT_endcap_MC_deno_density[pr1] = ((TH1D*)(h_pT_endcap_MC_deno[pr1]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr1])));
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
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr1]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr1]);
        for (Int_t i_bin=1; i_bin<=nPtBinBarrel; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr1]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr1]->GetBinContent(i_bin) /
                                                                  (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr1]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr1]->GetBinError(i_bin) /
                                                                (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            if (i_bin <= nPtBinEndcap)
            {
                h_pT_endcap_MC_deno_density[pr1]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr1]->GetBinContent(i_bin) /
                                                                      (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
                h_pT_endcap_MC_deno_density[pr1]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr1]->GetBinError(i_bin) /
                                                                    (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
            }
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
        h_pT_barrel_MC_deno_density[pr1]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr1]->SetFillColor(color);
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
        h_pT_barrel_MC_deno_density[pr1]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr1]->SetLineColor(color);
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
        h_pT_barrel_MC_deno_density[pr1]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr1]->SetDirectory(0);

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
            h_pT_barrel_MC_deno_density[_WJets_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr1]->Clone("h_pT_barrel_deno_density_WJets")));
            h_pT_endcap_MC_deno_density[_WJets_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr1]->Clone("h_pT_endcap_deno_density_WJets")));
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
            h_pT_barrel_MC_deno_density[_WJets_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_WJets_Full]->SetDirectory(0);
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
            h_pT_barrel_MC_deno_density[_WJets_Full]->Add(h_pT_barrel_MC_deno_density[pr1]);
            h_pT_endcap_MC_deno_density[_WJets_Full]->Add(h_pT_endcap_MC_deno_density[pr1]);
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
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr1]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr1]);

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
        h_pT_barrel_MC_deno_density[pr] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr])));
        h_pT_endcap_MC_deno_density[pr] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr])));
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
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr]);
        for (Int_t i_bin=1; i_bin<=nPtBinBarrel; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                  (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinError(i_bin) /
                                                               (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            if (i_bin <= nPtBinEndcap)
            {
                h_pT_endcap_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                     (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
                h_pT_endcap_MC_deno_density[pr]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                   (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
            }
        }

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
            h_pT_barrel_MC_deno_density[_DY_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr]->Clone("h_pT_barrel_deno_density_DY")));
            h_pT_endcap_MC_deno_density[_DY_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr]->Clone("h_pT_endcap_deno_density_DY")));
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
            h_pT_barrel_MC_deno_density[_DY_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_DY_Full]->SetDirectory(0);
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
            h_pT_barrel_MC_deno_density[_DY_Full]->Add(h_pT_barrel_MC_deno_density[pr]);
            h_pT_endcap_MC_deno_density[_DY_Full]->Add(h_pT_endcap_MC_deno_density[pr]);
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
        h_pT_barrel_MC_deno_density[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetFillColor(color);
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
        h_pT_barrel_MC_deno_density[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetLineColor(color);
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
        h_pT_barrel_MC_deno_density[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr]->SetDirectory(0);

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
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr]);

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
        h_pT_barrel_MC_deno_density[pr] = ((TH1D*)(h_pT_barrel_MC_deno[pr]->Clone("h_pT_barrel_MC_deno_density_"+fm.Procname[pr])));
        h_pT_endcap_MC_deno_density[pr] = ((TH1D*)(h_pT_endcap_MC_deno[pr]->Clone("h_pT_endcap_MC_deno_density_"+fm.Procname[pr])));
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
        removeNegativeBins(h_pT_barrel_MC_deno_density[pr]);
        removeNegativeBins(h_pT_endcap_MC_deno_density[pr]);
        for (Int_t i_bin=1; i_bin<=nPtBinBarrel; i_bin++)
        {
            h_pT_barrel_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                  (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            h_pT_barrel_MC_deno_density[pr]->SetBinError(i_bin, h_pT_barrel_MC_deno_density[pr]->GetBinError(i_bin) /
                                                               (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
            if (i_bin <= nPtBinEndcap)
            {
                h_pT_endcap_MC_deno_density[pr]->SetBinContent(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinContent(i_bin) /
                                                                     (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
                h_pT_endcap_MC_deno_density[pr]->SetBinError(i_bin, h_pT_endcap_MC_deno_density[pr]->GetBinError(i_bin) /
                                                                   (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
            }
        }

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
            h_pT_barrel_MC_deno_density[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_barrel_MC_deno_density[pr]->Clone("h_pT_barrel_deno_density_QCD")));
            h_pT_endcap_MC_deno_density[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_endcap_MC_deno_density[pr]->Clone("h_pT_endcap_deno_density_QCD")));
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
            h_pT_barrel_MC_deno_density[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_endcap_MC_deno_density[_QCDMuEnriched_Full]->SetDirectory(0);
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
            h_pT_barrel_MC_deno_density[_QCDMuEnriched_Full]->Add(h_pT_barrel_MC_deno_density[pr]);
            h_pT_endcap_MC_deno_density[_QCDMuEnriched_Full]->Add(h_pT_endcap_MC_deno_density[pr]);
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
        h_pT_barrel_MC_deno_density[pr]->SetFillColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetFillColor(color);
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
        h_pT_barrel_MC_deno_density[pr]->SetLineColor(color);
        h_pT_endcap_MC_deno_density[pr]->SetLineColor(color);
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
        h_pT_barrel_MC_deno_density[pr]->SetDirectory(0);
        h_pT_endcap_MC_deno_density[pr]->SetDirectory(0);

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
        s_pT_barrel_deno_density->Add(h_pT_barrel_MC_deno_density[pr]);
        s_pT_endcap_deno_density->Add(h_pT_endcap_MC_deno_density[pr]);

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
            h_pT_barrel_data_deno_density = ((TH1D*)(h_pT_barrel_data_deno->Clone("h_pT_barrel_data_deno_density")));
            h_pT_endcap_data_deno_density = ((TH1D*)(h_pT_endcap_data_deno->Clone("h_pT_endcap_data_deno_density")));
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
            removeNegativeBins(h_pT_barrel_data_deno_density);
            removeNegativeBins(h_pT_endcap_data_deno_density);
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
            h_pT_barrel_data_deno_density->Add(h_temp[8]);
            h_pT_endcap_data_deno_density->Add(h_temp[9]);
        }
    }
    for (Int_t i_bin=1; i_bin<=nPtBinBarrel; i_bin++)
    {
        h_pT_barrel_data_deno_density->SetBinContent(i_bin, h_pT_barrel_data_deno_density->GetBinContent(i_bin) /
                                                           (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
        h_pT_barrel_data_deno_density->SetBinError(i_bin, h_pT_barrel_data_deno_density->GetBinError(i_bin) /
                                                         (analyzer.ptbin_barrel[i_bin]-analyzer.ptbin_barrel[i_bin-1]));
        if (i_bin <= nPtBinEndcap)
        {
            h_pT_endcap_data_deno_density->SetBinContent(i_bin, h_pT_endcap_data_deno_density->GetBinContent(i_bin) /
                                                               (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
            h_pT_endcap_data_deno_density->SetBinError(i_bin, h_pT_endcap_data_deno_density->GetBinError(i_bin) /
                                                             (analyzer.ptbin_endcap[i_bin]-analyzer.ptbin_endcap[i_bin-1]));
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
    h_pT_barrel_data_deno_density->SetMarkerStyle(kFullDotLarge);
    h_pT_endcap_data_deno_density->SetMarkerStyle(kFullDotLarge);
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
    h_pT_barrel_data_deno_density->SetMarkerColor(kBlack);
    h_pT_endcap_data_deno_density->SetMarkerColor(kBlack);
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
    h_pT_barrel_data_deno_density->SetLineColor(kBlack);
    h_pT_endcap_data_deno_density->SetLineColor(kBlack);
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
    h_pT_barrel_data_deno_density->SetDirectory(0);
    h_pT_endcap_data_deno_density->SetDirectory(0);

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
    myRatioPlot_t *RP_pT_barrel_deno_density = new myRatioPlot_t("RP_pT_barrel_deno_density", s_pT_barrel_deno_density, h_pT_barrel_data_deno_density);
    myRatioPlot_t *RP_pT_endcap_deno_density = new myRatioPlot_t("RP_pT_endcap_deno_density", s_pT_endcap_deno_density, h_pT_endcap_data_deno_density);

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
    RP_pT_barrel_deno_density->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{barrel}}^{deno}) [GeV/c]", 52, 1000);
    RP_pT_endcap_deno_density->SetPlots("p_{#lower[-0.25]{T}} (#mu_{#lower[-0.4]{endcap}}^{deno}) [GeV/c]", 52, 1000);


    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_eta_data, "Data", "lp");
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
    RP_pT_barrel_deno_density->ImportLegend(legend);
    RP_pT_endcap_deno_density->ImportLegend(legend);

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

    RP_pT_barrel_deno_density->Draw(1, 1e8, 0, "HIST", "Number of events / GeV");
    RP_pT_endcap_deno_density->Draw(1, 1e8, 0, "HIST", "Number of events / GeV");

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

    TH1D *h_barrel_MC_ctrl_50to70  [_EndOf_Data_Special],
         *h_barrel_MC_nume_50to70  [_EndOf_Data_Special],
         *h_endcap_MC_ctrl_50to70  [_EndOf_Data_Special],
         *h_endcap_MC_nume_50to70  [_EndOf_Data_Special],
         *h_barrel_MC_ctrl_70to100 [_EndOf_Data_Special],
         *h_barrel_MC_nume_70to100 [_EndOf_Data_Special],
         *h_endcap_MC_ctrl_70to100 [_EndOf_Data_Special],
         *h_endcap_MC_nume_70to100 [_EndOf_Data_Special],
         *h_barrel_MC_ctrl_100to500[_EndOf_Data_Special],
         *h_barrel_MC_nume_100to500[_EndOf_Data_Special],
         *h_endcap_MC_ctrl_100to500[_EndOf_Data_Special],
         *h_endcap_MC_nume_100to500[_EndOf_Data_Special],
         *h_barrel_data_ctrl_50to70,
         *h_barrel_data_nume_50to70,
         *h_endcap_data_ctrl_50to70,
         *h_endcap_data_nume_50to70,
         *h_barrel_data_ctrl_70to100,
         *h_barrel_data_nume_70to100,
         *h_endcap_data_ctrl_70to100,
         *h_endcap_data_nume_70to100,
         *h_barrel_data_ctrl_100to500,
         *h_barrel_data_nume_100to500,
         *h_endcap_data_ctrl_100to500,
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
        file->GetObject("h_pT_barrel_ctrl_50to70",   h_barrel_MC_ctrl_50to70 [pr1]);
        file->GetObject("h_pT_endcap_ctrl_50to70",   h_endcap_MC_ctrl_50to70 [pr1]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr1]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr1]);
        file->GetObject("h_pT_barrel_ctrl_70to100",  h_barrel_MC_ctrl_70to100[pr1]);
        file->GetObject("h_pT_endcap_ctrl_70to100",  h_endcap_MC_ctrl_70to100[pr1]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr1]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr1]);
        file->GetObject("h_pT_barrel_ctrl_100to500", h_barrel_MC_ctrl_100to500[pr1]);
        file->GetObject("h_pT_endcap_ctrl_100to500", h_endcap_MC_ctrl_100to500[pr1]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr1]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr1]);

        removeNegativeBins(h_barrel_MC_ctrl_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_ctrl_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr1]);
        removeNegativeBins(h_barrel_MC_ctrl_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_ctrl_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr1]);
        removeNegativeBins(h_barrel_MC_ctrl_100to500[pr1]);
        removeNegativeBins(h_endcap_MC_ctrl_100to500[pr1]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr1]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr1]);

        h_barrel_MC_ctrl_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_ctrl_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr1]->SetDirectory(0);
        h_barrel_MC_ctrl_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_ctrl_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr1]->SetDirectory(0);
        h_barrel_MC_ctrl_100to500[pr1]->SetDirectory(0);
        h_endcap_MC_ctrl_100to500[pr1]->SetDirectory(0);
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

        h_barrel_MC_ctrl_50to70  [pr1]->SetFillColor(color);
        h_endcap_MC_ctrl_50to70  [pr1]->SetFillColor(color);
        h_barrel_MC_nume_50to70  [pr1]->SetFillColor(color);
        h_endcap_MC_nume_50to70  [pr1]->SetFillColor(color);
        h_barrel_MC_ctrl_70to100 [pr1]->SetFillColor(color);
        h_endcap_MC_ctrl_70to100 [pr1]->SetFillColor(color);
        h_barrel_MC_nume_70to100 [pr1]->SetFillColor(color);
        h_endcap_MC_nume_70to100 [pr1]->SetFillColor(color);
        h_barrel_MC_ctrl_100to500[pr1]->SetFillColor(color);
        h_endcap_MC_ctrl_100to500[pr1]->SetFillColor(color);
        h_barrel_MC_nume_100to500[pr1]->SetFillColor(color);
        h_endcap_MC_nume_100to500[pr1]->SetFillColor(color);

        h_barrel_MC_ctrl_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_ctrl_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_endcap_MC_nume_50to70  [pr1]->SetLineColor(color);
        h_barrel_MC_ctrl_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_ctrl_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_endcap_MC_nume_70to100 [pr1]->SetLineColor(color);
        h_barrel_MC_ctrl_100to500[pr1]->SetLineColor(color);
        h_endcap_MC_ctrl_100to500[pr1]->SetLineColor(color);
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
    h_barrel_MC_ctrl_50to70  [_ttbar]->Add(h_barrel_MC_ctrl_50to70  [_ttbar_700to1000]);
    h_endcap_MC_ctrl_50to70  [_ttbar]->Add(h_endcap_MC_ctrl_50to70  [_ttbar_700to1000]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_700to1000]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_700to1000]);
    h_barrel_MC_ctrl_70to100 [_ttbar]->Add(h_barrel_MC_ctrl_70to100 [_ttbar_700to1000]);
    h_endcap_MC_ctrl_70to100 [_ttbar]->Add(h_endcap_MC_ctrl_70to100 [_ttbar_700to1000]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_700to1000]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_700to1000]);
    h_barrel_MC_ctrl_100to500[_ttbar]->Add(h_barrel_MC_ctrl_100to500[_ttbar_700to1000]);
    h_endcap_MC_ctrl_100to500[_ttbar]->Add(h_endcap_MC_ctrl_100to500[_ttbar_700to1000]);
    h_barrel_MC_nume_100to500[_ttbar]->Add(h_barrel_MC_nume_100to500[_ttbar_700to1000]);
    h_endcap_MC_nume_100to500[_ttbar]->Add(h_endcap_MC_nume_100to500[_ttbar_700to1000]);

    h_barrel_MC_ctrl_50to70  [_ttbar]->Add(h_barrel_MC_ctrl_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_ctrl_50to70  [_ttbar]->Add(h_endcap_MC_ctrl_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_nume_50to70  [_ttbar]->Add(h_barrel_MC_nume_50to70  [_ttbar_1000toInf]);
    h_endcap_MC_nume_50to70  [_ttbar]->Add(h_endcap_MC_nume_50to70  [_ttbar_1000toInf]);
    h_barrel_MC_ctrl_70to100 [_ttbar]->Add(h_barrel_MC_ctrl_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_ctrl_70to100 [_ttbar]->Add(h_endcap_MC_ctrl_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_nume_70to100 [_ttbar]->Add(h_barrel_MC_nume_70to100 [_ttbar_1000toInf]);
    h_endcap_MC_nume_70to100 [_ttbar]->Add(h_endcap_MC_nume_70to100 [_ttbar_1000toInf]);
    h_barrel_MC_ctrl_100to500[_ttbar]->Add(h_barrel_MC_ctrl_100to500[_ttbar_1000toInf]);
    h_endcap_MC_ctrl_100to500[_ttbar]->Add(h_endcap_MC_ctrl_100to500[_ttbar_1000toInf]);
    h_barrel_MC_nume_100to500[_ttbar]->Add(h_barrel_MC_nume_100to500[_ttbar_1000toInf]);
    h_endcap_MC_nume_100to500[_ttbar]->Add(h_endcap_MC_nume_100to500[_ttbar_1000toInf]);

    h_barrel_MC_ctrl_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_ctrl_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_ctrl_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_ctrl_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_ctrl_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_ctrl_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_barrel_MC_nume_100to500[_ttbar]->SetFillColor(kCyan+2);
    h_endcap_MC_nume_100to500[_ttbar]->SetFillColor(kCyan+2);

    h_barrel_MC_ctrl_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_ctrl_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_50to70  [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_ctrl_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_ctrl_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_70to100 [_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_ctrl_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_ctrl_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_barrel_MC_nume_100to500[_ttbar]->SetLineColor(kCyan+2);
    h_endcap_MC_nume_100to500[_ttbar]->SetLineColor(kCyan+2);


    h_barrel_MC_ctrl_50to70  [_WJets]->Add(h_barrel_MC_ctrl_50to70  [_WJets_ext2v5]);
    h_endcap_MC_ctrl_50to70  [_WJets]->Add(h_endcap_MC_ctrl_50to70  [_WJets_ext2v5]);
    h_barrel_MC_nume_50to70  [_WJets]->Add(h_barrel_MC_nume_50to70  [_WJets_ext2v5]);
    h_endcap_MC_nume_50to70  [_WJets]->Add(h_endcap_MC_nume_50to70  [_WJets_ext2v5]);
    h_barrel_MC_ctrl_70to100 [_WJets]->Add(h_barrel_MC_ctrl_70to100 [_WJets_ext2v5]);
    h_endcap_MC_ctrl_70to100 [_WJets]->Add(h_endcap_MC_ctrl_70to100 [_WJets_ext2v5]);
    h_barrel_MC_nume_70to100 [_WJets]->Add(h_barrel_MC_nume_70to100 [_WJets_ext2v5]);
    h_endcap_MC_nume_70to100 [_WJets]->Add(h_endcap_MC_nume_70to100 [_WJets_ext2v5]);
    h_barrel_MC_ctrl_100to500[_WJets]->Add(h_barrel_MC_ctrl_100to500[_WJets_ext2v5]);
    h_endcap_MC_ctrl_100to500[_WJets]->Add(h_endcap_MC_ctrl_100to500[_WJets_ext2v5]);
    h_barrel_MC_nume_100to500[_WJets]->Add(h_barrel_MC_nume_100to500[_WJets_ext2v5]);
    h_endcap_MC_nume_100to500[_WJets]->Add(h_endcap_MC_nume_100to500[_WJets_ext2v5]);

    h_barrel_MC_ctrl_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_ctrl_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_ctrl_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_ctrl_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_ctrl_100to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_ctrl_100to500[_WJets]->SetFillColor(kRed-2);
    h_barrel_MC_nume_100to500[_WJets]->SetFillColor(kRed-2);
    h_endcap_MC_nume_100to500[_WJets]->SetFillColor(kRed-2);

    h_barrel_MC_ctrl_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_ctrl_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_50to70  [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_ctrl_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_ctrl_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_70to100 [_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_ctrl_100to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_ctrl_100to500[_WJets]->SetLineColor(kRed-2);
    h_barrel_MC_nume_100to500[_WJets]->SetLineColor(kRed-2);
    h_endcap_MC_nume_100to500[_WJets]->SetLineColor(kRed-2);

    // DY
    for (Process_t pr = _DY_10to50; pr <= _DY_2000to3000; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_pT_barrel_ctrl_50to70",   h_barrel_MC_ctrl_50to70  [pr]);
        file->GetObject("h_pT_endcap_ctrl_50to70",   h_endcap_MC_ctrl_50to70  [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70  [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70  [pr]);
        file->GetObject("h_pT_barrel_ctrl_70to100",  h_barrel_MC_ctrl_70to100 [pr]);
        file->GetObject("h_pT_endcap_ctrl_70to100",  h_endcap_MC_ctrl_70to100 [pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100 [pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100 [pr]);
        file->GetObject("h_pT_barrel_ctrl_100to500", h_barrel_MC_ctrl_100to500[pr]);
        file->GetObject("h_pT_endcap_ctrl_100to500", h_endcap_MC_ctrl_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

        removeNegativeBins(h_barrel_MC_ctrl_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_ctrl_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_ctrl_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_ctrl_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_ctrl_100to500[pr]);
        removeNegativeBins(h_endcap_MC_ctrl_100to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr]);

        h_barrel_MC_ctrl_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_ctrl_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_ctrl_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_ctrl_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);;
        h_barrel_MC_ctrl_100to500[pr]->SetDirectory(0);
        h_endcap_MC_ctrl_100to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to500[pr]->SetDirectory(0);

        if (pr == _DY_10to50)
        {
            h_barrel_MC_ctrl_50to70  [_DY_Full] = ((TH1D*)(h_barrel_MC_ctrl_50to70  [pr]->Clone("h_barrel_MC_ctrl_DY_50to70")));
            h_endcap_MC_ctrl_50to70  [_DY_Full] = ((TH1D*)(h_endcap_MC_ctrl_50to70  [pr]->Clone("h_endcap_MC_ctrl_DY_50to70")));
            h_barrel_MC_nume_50to70  [_DY_Full] = ((TH1D*)(h_barrel_MC_nume_50to70  [pr]->Clone("h_barrel_MC_nume_DY_50to70")));
            h_endcap_MC_nume_50to70  [_DY_Full] = ((TH1D*)(h_endcap_MC_nume_50to70  [pr]->Clone("h_endcap_MC_nume_DY_50to70")));
            h_barrel_MC_ctrl_70to100 [_DY_Full] = ((TH1D*)(h_barrel_MC_ctrl_70to100 [pr]->Clone("h_barrel_MC_ctrl_DY_70to100")));
            h_endcap_MC_ctrl_70to100 [_DY_Full] = ((TH1D*)(h_endcap_MC_ctrl_70to100 [pr]->Clone("h_endcap_MC_ctrl_DY_70to100")));
            h_barrel_MC_nume_70to100 [_DY_Full] = ((TH1D*)(h_barrel_MC_nume_70to100 [pr]->Clone("h_barrel_MC_nume_DY_70to100")));
            h_endcap_MC_nume_70to100 [_DY_Full] = ((TH1D*)(h_endcap_MC_nume_70to100 [pr]->Clone("h_endcap_MC_nume_DY_70to100")));
            h_barrel_MC_ctrl_100to500[_DY_Full] = ((TH1D*)(h_barrel_MC_ctrl_100to500[pr]->Clone("h_barrel_MC_ctrl_DY_100to500")));
            h_endcap_MC_ctrl_100to500[_DY_Full] = ((TH1D*)(h_endcap_MC_ctrl_100to500[pr]->Clone("h_endcap_MC_ctrl_DY_100to500")));
            h_barrel_MC_nume_100to500[_DY_Full] = ((TH1D*)(h_barrel_MC_nume_100to500[pr]->Clone("h_barrel_MC_nume_DY_100to500")));
            h_endcap_MC_nume_100to500[_DY_Full] = ((TH1D*)(h_endcap_MC_nume_100to500[pr]->Clone("h_endcap_MC_nume_DY_100to500")));

            h_barrel_MC_ctrl_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_DY_Full]->SetDirectory(0);
            h_barrel_MC_ctrl_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_DY_Full]->SetDirectory(0);
            h_barrel_MC_ctrl_100to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_100to500[_DY_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to500[_DY_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to500[_DY_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_ctrl_50to70  [_DY_Full]->Add(h_barrel_MC_ctrl_50to70  [pr]);
            h_endcap_MC_ctrl_50to70  [_DY_Full]->Add(h_endcap_MC_ctrl_50to70  [pr]);
            h_barrel_MC_nume_50to70  [_DY_Full]->Add(h_barrel_MC_nume_50to70  [pr]);
            h_endcap_MC_nume_50to70  [_DY_Full]->Add(h_endcap_MC_nume_50to70  [pr]);
            h_barrel_MC_ctrl_70to100 [_DY_Full]->Add(h_barrel_MC_ctrl_70to100 [pr]);
            h_endcap_MC_ctrl_70to100 [_DY_Full]->Add(h_endcap_MC_ctrl_70to100 [pr]);
            h_barrel_MC_nume_70to100 [_DY_Full]->Add(h_barrel_MC_nume_70to100 [pr]);
            h_endcap_MC_nume_70to100 [_DY_Full]->Add(h_endcap_MC_nume_70to100 [pr]);
            h_barrel_MC_ctrl_100to500[_DY_Full]->Add(h_barrel_MC_ctrl_100to500[pr]);
            h_endcap_MC_ctrl_100to500[_DY_Full]->Add(h_endcap_MC_ctrl_100to500[pr]);
            h_barrel_MC_nume_100to500[_DY_Full]->Add(h_barrel_MC_nume_100to500[pr]);
            h_endcap_MC_nume_100to500[_DY_Full]->Add(h_endcap_MC_nume_100to500[pr]);
        }
        file->Close();
    }

    h_barrel_MC_ctrl_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_ctrl_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_50to70  [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_ctrl_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_ctrl_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_70to100 [_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_ctrl_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_ctrl_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_barrel_MC_nume_100to500[_DY_Full]->SetFillColor(kOrange-5);
    h_endcap_MC_nume_100to500[_DY_Full]->SetFillColor(kOrange-5);

    h_barrel_MC_ctrl_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_ctrl_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_50to70  [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_ctrl_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_ctrl_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_70to100 [_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_ctrl_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_ctrl_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_barrel_MC_nume_100to500[_DY_Full]->SetLineColor(kOrange-5);
    h_endcap_MC_nume_100to500[_DY_Full]->SetLineColor(kOrange-5);

    // QCD
    for (Process_t pr = _QCDMuEnriched_15to20; pr <= _QCDMuEnriched_1000toInf; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Muon/FR_Hist_Mu_"+fm.Procname[pr]+".root", "READ");
        file->GetObject("h_pT_barrel_ctrl_50to70",   h_barrel_MC_ctrl_50to70 [pr]);
        file->GetObject("h_pT_endcap_ctrl_50to70",   h_endcap_MC_ctrl_50to70 [pr]);
        file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_MC_nume_50to70 [pr]);
        file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_MC_nume_50to70 [pr]);
        file->GetObject("h_pT_barrel_ctrl_70to100",  h_barrel_MC_ctrl_70to100[pr]);
        file->GetObject("h_pT_endcap_ctrl_70to100",  h_endcap_MC_ctrl_70to100[pr]);
        file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_MC_nume_70to100[pr]);
        file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_MC_nume_70to100[pr]);
        file->GetObject("h_pT_barrel_ctrl_100to500", h_barrel_MC_ctrl_100to500[pr]);
        file->GetObject("h_pT_endcap_ctrl_100to500", h_endcap_MC_ctrl_100to500[pr]);
        file->GetObject("h_pT_barrel_nume_100to500", h_barrel_MC_nume_100to500[pr]);
        file->GetObject("h_pT_endcap_nume_100to500", h_endcap_MC_nume_100to500[pr]);

        removeNegativeBins(h_barrel_MC_ctrl_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_ctrl_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_nume_50to70  [pr]);
        removeNegativeBins(h_endcap_MC_nume_50to70  [pr]);
        removeNegativeBins(h_barrel_MC_ctrl_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_ctrl_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_nume_70to100 [pr]);
        removeNegativeBins(h_endcap_MC_nume_70to100 [pr]);
        removeNegativeBins(h_barrel_MC_ctrl_100to500[pr]);
        removeNegativeBins(h_endcap_MC_ctrl_100to500[pr]);
        removeNegativeBins(h_barrel_MC_nume_100to500[pr]);
        removeNegativeBins(h_endcap_MC_nume_100to500[pr]);

        h_barrel_MC_ctrl_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_ctrl_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_nume_50to70  [pr]->SetDirectory(0);
        h_endcap_MC_nume_50to70  [pr]->SetDirectory(0);
        h_barrel_MC_ctrl_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_ctrl_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_nume_70to100 [pr]->SetDirectory(0);
        h_endcap_MC_nume_70to100 [pr]->SetDirectory(0);
        h_barrel_MC_ctrl_100to500[pr]->SetDirectory(0);
        h_endcap_MC_ctrl_100to500[pr]->SetDirectory(0);
        h_barrel_MC_nume_100to500[pr]->SetDirectory(0);
        h_endcap_MC_nume_100to500[pr]->SetDirectory(0);

        if (pr == _QCDMuEnriched_15to20)
        {
            h_barrel_MC_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_ctrl_50to70  [pr]->Clone("h_barrel_MC_ctrl_QCD_50to70")));
            h_endcap_MC_ctrl_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_ctrl_50to70  [pr]->Clone("h_endcap_MC_ctrl_QCD_50to70")));
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_50to70  [pr]->Clone("h_barrel_MC_nume_QCD_50to70")));
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_50to70  [pr]->Clone("h_endcap_MC_nume_QCD_50to70")));
            h_barrel_MC_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_ctrl_70to100 [pr]->Clone("h_barrel_MC_ctrl_QCD_70to100")));
            h_endcap_MC_ctrl_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_ctrl_70to100 [pr]->Clone("h_endcap_MC_ctrl_QCD_70to100")));
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_70to100 [pr]->Clone("h_barrel_MC_nume_QCD_70to100")));
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_70to100 [pr]->Clone("h_endcap_MC_nume_QCD_70to100")));
            h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_ctrl_100to500[pr]->Clone("h_barrel_MC_ctrl_QCD_100to500")));
            h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_ctrl_100to500[pr]->Clone("h_endcap_MC_ctrl_QCD_100to500")));
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_barrel_MC_nume_100to500[pr]->Clone("h_barrel_MC_nume_QCD_100to500")));
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full] = ((TH1D*)(h_endcap_MC_nume_100to500[pr]->Clone("h_endcap_MC_nume_QCD_100to500")));

            h_barrel_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else
        {
            h_barrel_MC_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_barrel_MC_ctrl_50to70 [pr]);
            h_endcap_MC_ctrl_50to70  [_QCDMuEnriched_Full]->Add(h_endcap_MC_ctrl_50to70 [pr]);
            h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_50to70 [pr]);
            h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_50to70 [pr]);
            h_barrel_MC_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_barrel_MC_ctrl_70to100[pr]);
            h_endcap_MC_ctrl_70to100 [_QCDMuEnriched_Full]->Add(h_endcap_MC_ctrl_70to100[pr]);
            h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_70to100[pr]);
            h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_70to100[pr]);
            h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_ctrl_100to500[pr]);
            h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_ctrl_100to500[pr]);
            h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->Add(h_barrel_MC_nume_100to500[pr]);
            h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->Add(h_endcap_MC_nume_100to500[pr]);
        }
        file->Close();
    }

    h_barrel_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);
    h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->SetFillColor(kRed+3);

    h_barrel_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_ctrl_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_50to70  [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_ctrl_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_nume_70to100 [_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
    h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]->SetLineColor(kRed+3);
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
            file->GetObject("h_pT_barrel_ctrl_50to70",   h_barrel_data_ctrl_50to70 );
            file->GetObject("h_pT_endcap_ctrl_50to70",   h_endcap_data_ctrl_50to70 );
            file->GetObject("h_pT_barrel_nume_50to70",   h_barrel_data_nume_50to70 );
            file->GetObject("h_pT_endcap_nume_50to70",   h_endcap_data_nume_50to70 );
            file->GetObject("h_pT_barrel_ctrl_70to100",  h_barrel_data_ctrl_70to100);
            file->GetObject("h_pT_endcap_ctrl_70to100",  h_endcap_data_ctrl_70to100);
            file->GetObject("h_pT_barrel_nume_70to100",  h_barrel_data_nume_70to100);
            file->GetObject("h_pT_endcap_nume_70to100",  h_endcap_data_nume_70to100);
            file->GetObject("h_pT_barrel_ctrl_100to500", h_barrel_data_ctrl_100to500);
            file->GetObject("h_pT_endcap_ctrl_100to500", h_endcap_data_ctrl_100to500);
            file->GetObject("h_pT_barrel_nume_100to500", h_barrel_data_nume_100to500);
            file->GetObject("h_pT_endcap_nume_100to500", h_endcap_data_nume_100to500);

            removeNegativeBins(h_barrel_data_ctrl_50to70);
            removeNegativeBins(h_endcap_data_ctrl_50to70);
            removeNegativeBins(h_barrel_data_nume_50to70);
            removeNegativeBins(h_endcap_data_nume_50to70);
            removeNegativeBins(h_barrel_data_ctrl_70to100);
            removeNegativeBins(h_endcap_data_ctrl_70to100);
            removeNegativeBins(h_barrel_data_nume_70to100);
            removeNegativeBins(h_endcap_data_nume_70to100);
            removeNegativeBins(h_barrel_data_ctrl_100to500);
            removeNegativeBins(h_endcap_data_ctrl_100to500);
            removeNegativeBins(h_barrel_data_nume_100to500);
            removeNegativeBins(h_endcap_data_nume_100to500);
        }
        else
        {
            file->GetObject("h_pT_barrel_ctrl_50to70",   h_temp[0]);
            file->GetObject("h_pT_endcap_ctrl_50to70",   h_temp[1]);
            file->GetObject("h_pT_barrel_nume_50to70",   h_temp[2]);
            file->GetObject("h_pT_endcap_nume_50to70",   h_temp[3]);
            file->GetObject("h_pT_barrel_ctrl_70to100",  h_temp[4]);
            file->GetObject("h_pT_endcap_ctrl_70to100",  h_temp[5]);
            file->GetObject("h_pT_barrel_nume_70to100",  h_temp[6]);
            file->GetObject("h_pT_endcap_nume_70to100",  h_temp[7]);
            file->GetObject("h_pT_barrel_ctrl_100to500", h_temp[8]);
            file->GetObject("h_pT_endcap_ctrl_100to500", h_temp[9]);
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

            h_barrel_data_ctrl_50to70  ->Add(h_temp[0]);
            h_endcap_data_ctrl_50to70  ->Add(h_temp[1]);
            h_barrel_data_nume_50to70  ->Add(h_temp[2]);
            h_endcap_data_nume_50to70  ->Add(h_temp[3]);
            h_barrel_data_ctrl_70to100 ->Add(h_temp[4]);
            h_endcap_data_ctrl_70to100 ->Add(h_temp[5]);
            h_barrel_data_nume_70to100 ->Add(h_temp[6]);
            h_endcap_data_nume_70to100 ->Add(h_temp[7]);
            h_barrel_data_ctrl_100to500->Add(h_temp[8]);
            h_endcap_data_ctrl_100to500->Add(h_temp[9]);
            h_barrel_data_nume_100to500->Add(h_temp[10]);
            h_endcap_data_nume_100to500->Add(h_temp[11]);
        }
    }

    h_barrel_data_ctrl_50to70  ->SetDirectory(0);
    h_endcap_data_ctrl_50to70  ->SetDirectory(0);
    h_barrel_data_nume_50to70  ->SetDirectory(0);
    h_endcap_data_nume_50to70  ->SetDirectory(0);
    h_barrel_data_ctrl_70to100 ->SetDirectory(0);
    h_endcap_data_ctrl_70to100 ->SetDirectory(0);
    h_barrel_data_nume_70to100 ->SetDirectory(0);
    h_endcap_data_nume_70to100 ->SetDirectory(0);
    h_barrel_data_ctrl_100to500->SetDirectory(0);
    h_endcap_data_ctrl_100to500->SetDirectory(0);
    h_barrel_data_nume_100to500->SetDirectory(0);
    h_endcap_data_nume_100to500->SetDirectory(0);

    h_barrel_data_ctrl_50to70  ->SetLineColor(kBlack);
    h_endcap_data_ctrl_50to70  ->SetLineColor(kBlack);
    h_barrel_data_nume_50to70  ->SetLineColor(kBlack);
    h_endcap_data_nume_50to70  ->SetLineColor(kBlack);
    h_barrel_data_ctrl_70to100 ->SetLineColor(kBlack);
    h_endcap_data_ctrl_70to100 ->SetLineColor(kBlack);
    h_barrel_data_nume_70to100 ->SetLineColor(kBlack);
    h_endcap_data_nume_70to100 ->SetLineColor(kBlack);
    h_barrel_data_ctrl_100to500->SetLineColor(kBlack);
    h_endcap_data_ctrl_100to500->SetLineColor(kBlack);
    h_barrel_data_nume_100to500->SetLineColor(kBlack);
    h_endcap_data_nume_100to500->SetLineColor(kBlack);

    h_barrel_data_ctrl_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_ctrl_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_endcap_data_nume_50to70  ->SetMarkerColor(kBlack);
    h_barrel_data_ctrl_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_ctrl_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_endcap_data_nume_70to100 ->SetMarkerColor(kBlack);
    h_barrel_data_ctrl_100to500->SetMarkerColor(kBlack);
    h_endcap_data_ctrl_100to500->SetMarkerColor(kBlack);
    h_barrel_data_nume_100to500->SetMarkerColor(kBlack);
    h_endcap_data_nume_100to500->SetMarkerColor(kBlack);

    h_barrel_data_ctrl_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_ctrl_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_50to70  ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_ctrl_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_ctrl_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_70to100 ->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_ctrl_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_ctrl_100to500->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_nume_100to500->SetMarkerStyle(kFullDotLarge);
    h_endcap_data_nume_100to500->SetMarkerStyle(kFullDotLarge);

//-------------------------------------- SCALING -----------------------------------------------------

    h_barrel_MC_ctrl_50to70[_DY_Full]             ->Scale(5.4747e+04 / h_barrel_MC_ctrl_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_ctrl_50to70[_WW]                  ->Scale(1.5728e+03 / h_barrel_MC_ctrl_50to70[_WW]                 ->Integral());
    h_barrel_MC_ctrl_50to70[_WZ]                  ->Scale(6.1305e+02 / h_barrel_MC_ctrl_50to70[_WZ]                 ->Integral());
    h_barrel_MC_ctrl_50to70[_ZZ]                  ->Scale(2.5470e+02 / h_barrel_MC_ctrl_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_ctrl_50to70[_tW]                  ->Scale(2.8822e+03 / h_barrel_MC_ctrl_50to70[_tW]                 ->Integral());
    h_barrel_MC_ctrl_50to70[_tbarW]               ->Scale(2.9110e+03 / h_barrel_MC_ctrl_50to70[_tbarW]              ->Integral());
    h_barrel_MC_ctrl_50to70[_ttbar]               ->Scale(9.0481e+04 / h_barrel_MC_ctrl_50to70[_ttbar]              ->Integral());
    h_barrel_MC_ctrl_50to70[_WJets]               ->Scale(1.0401e+05 / h_barrel_MC_ctrl_50to70[_WJets]              ->Integral());
    h_barrel_MC_ctrl_50to70[_QCDMuEnriched_Full]  ->Scale(1.4761e+07 / h_barrel_MC_ctrl_50to70[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_ctrl_70to100[_DY_Full]            ->Scale(1.5125e+04 / h_barrel_MC_ctrl_70to100[_DY_Full]            ->Integral());
    h_barrel_MC_ctrl_70to100[_WW]                 ->Scale(6.3972e+02 / h_barrel_MC_ctrl_70to100[_WW]                 ->Integral());
    h_barrel_MC_ctrl_70to100[_WZ]                 ->Scale(2.6216e+02 / h_barrel_MC_ctrl_70to100[_WZ]                 ->Integral());
    h_barrel_MC_ctrl_70to100[_ZZ]                 ->Scale(1.4548e+02 / h_barrel_MC_ctrl_70to100[_ZZ]                 ->Integral());
    h_barrel_MC_ctrl_70to100[_tW]                 ->Scale(1.5772e+03 / h_barrel_MC_ctrl_70to100[_tW]                 ->Integral());
    h_barrel_MC_ctrl_70to100[_tbarW]              ->Scale(1.5400e+03 / h_barrel_MC_ctrl_70to100[_tbarW]              ->Integral());
    h_barrel_MC_ctrl_70to100[_ttbar]              ->Scale(4.4142e+04 / h_barrel_MC_ctrl_70to100[_ttbar]              ->Integral());
    h_barrel_MC_ctrl_70to100[_WJets]              ->Scale(2.8937e+04 / h_barrel_MC_ctrl_70to100[_WJets]              ->Integral());
    h_barrel_MC_ctrl_70to100[_QCDMuEnriched_Full] ->Scale(4.1945e+06 / h_barrel_MC_ctrl_70to100[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_ctrl_100to500[_DY_Full]           ->Scale(4.8091e+03 / h_barrel_MC_ctrl_100to500[_DY_Full]            ->Integral());
    h_barrel_MC_ctrl_100to500[_WW]                ->Scale(3.0715e+02 / h_barrel_MC_ctrl_100to500[_WW]                 ->Integral());
    h_barrel_MC_ctrl_100to500[_WZ]                ->Scale(1.0903e+02 / h_barrel_MC_ctrl_100to500[_WZ]                 ->Integral());
    h_barrel_MC_ctrl_100to500[_ZZ]                ->Scale(5.2073e+01 / h_barrel_MC_ctrl_100to500[_ZZ]                 ->Integral());
    h_barrel_MC_ctrl_100to500[_tW]                ->Scale(1.1978e+03 / h_barrel_MC_ctrl_100to500[_tW]                 ->Integral());
    h_barrel_MC_ctrl_100to500[_tbarW]             ->Scale(1.2146e+03 / h_barrel_MC_ctrl_100to500[_tbarW]              ->Integral());
    h_barrel_MC_ctrl_100to500[_ttbar]             ->Scale(2.9037e+04 / h_barrel_MC_ctrl_100to500[_ttbar]              ->Integral());
    h_barrel_MC_ctrl_100to500[_WJets]             ->Scale(8.7085e+03 / h_barrel_MC_ctrl_100to500[_WJets]              ->Integral());
    h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]->Scale(8.8507e+05 / h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full] ->Integral());

    h_endcap_MC_ctrl_50to70[_DY_Full]             ->Scale(2.6401e+04 / h_endcap_MC_ctrl_50to70[_DY_Full]            ->Integral());
    h_endcap_MC_ctrl_50to70[_WW]                  ->Scale(8.2832e+02 / h_endcap_MC_ctrl_50to70[_WW]                 ->Integral());
    h_endcap_MC_ctrl_50to70[_WZ]                  ->Scale(2.8933e+02 / h_endcap_MC_ctrl_50to70[_WZ]                 ->Integral());
    h_endcap_MC_ctrl_50to70[_ZZ]                  ->Scale(1.4463e+02 / h_endcap_MC_ctrl_50to70[_ZZ]                 ->Integral());
    h_endcap_MC_ctrl_50to70[_tW]                  ->Scale(9.1301e+02 / h_endcap_MC_ctrl_50to70[_tW]                 ->Integral());
    h_endcap_MC_ctrl_50to70[_tbarW]               ->Scale(8.8836e+02 / h_endcap_MC_ctrl_50to70[_tbarW]              ->Integral());
    h_endcap_MC_ctrl_50to70[_ttbar]               ->Scale(2.8681e+04 / h_endcap_MC_ctrl_50to70[_ttbar]              ->Integral());
    h_endcap_MC_ctrl_50to70[_WJets]               ->Scale(4.9213e+04 / h_endcap_MC_ctrl_50to70[_WJets]              ->Integral());
    h_endcap_MC_ctrl_50to70[_QCDMuEnriched_Full]  ->Scale(6.8016e+06 / h_endcap_MC_ctrl_50to70[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_ctrl_70to100[_DY_Full]            ->Scale(7.7330e+03 / h_endcap_MC_ctrl_70to100[_DY_Full]            ->Integral());
    h_endcap_MC_ctrl_70to100[_WW]                 ->Scale(3.3516e+02 / h_endcap_MC_ctrl_70to100[_WW]                 ->Integral());
    h_endcap_MC_ctrl_70to100[_WZ]                 ->Scale(1.2805e+02 / h_endcap_MC_ctrl_70to100[_WZ]                 ->Integral());
    h_endcap_MC_ctrl_70to100[_ZZ]                 ->Scale(5.5725e+01 / h_endcap_MC_ctrl_70to100[_ZZ]                 ->Integral());
    h_endcap_MC_ctrl_70to100[_tW]                 ->Scale(4.4211e+02 / h_endcap_MC_ctrl_70to100[_tW]                 ->Integral());
    h_endcap_MC_ctrl_70to100[_tbarW]              ->Scale(4.4114e+02 / h_endcap_MC_ctrl_70to100[_tbarW]              ->Integral());
    h_endcap_MC_ctrl_70to100[_ttbar]              ->Scale(1.2446e+04 / h_endcap_MC_ctrl_70to100[_ttbar]              ->Integral());
    h_endcap_MC_ctrl_70to100[_WJets]              ->Scale(1.4291e+04 / h_endcap_MC_ctrl_70to100[_WJets]              ->Integral());
    h_endcap_MC_ctrl_70to100[_QCDMuEnriched_Full] ->Scale(1.7431e+06 / h_endcap_MC_ctrl_70to100[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_ctrl_100to500[_DY_Full]           ->Scale(3.4277e+03 / h_endcap_MC_ctrl_100to500[_DY_Full]            ->Integral());
    h_endcap_MC_ctrl_100to500[_WW]                ->Scale(1.7181e+02 / h_endcap_MC_ctrl_100to500[_WW]                 ->Integral());
    h_endcap_MC_ctrl_100to500[_WZ]                ->Scale(4.2215e+01 / h_endcap_MC_ctrl_100to500[_WZ]                 ->Integral());
    h_endcap_MC_ctrl_100to500[_ZZ]                ->Scale(3.3979e+01 / h_endcap_MC_ctrl_100to500[_ZZ]                 ->Integral());
    h_endcap_MC_ctrl_100to500[_tW]                ->Scale(2.7577e+02 / h_endcap_MC_ctrl_100to500[_tW]                 ->Integral());
    h_endcap_MC_ctrl_100to500[_tbarW]             ->Scale(2.5955e+02 / h_endcap_MC_ctrl_100to500[_tbarW]              ->Integral());
    h_endcap_MC_ctrl_100to500[_ttbar]             ->Scale(6.6796e+03 / h_endcap_MC_ctrl_100to500[_ttbar]              ->Integral());
    h_endcap_MC_ctrl_100to500[_WJets]             ->Scale(1.0192e+04 / h_endcap_MC_ctrl_100to500[_WJets]              ->Integral());
    h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]->Scale(2.9850e+05 / h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full] ->Integral());

    h_barrel_MC_nume_50to70[_DY_Full]             ->Scale(2.3263e+06 / h_barrel_MC_nume_50to70[_DY_Full]            ->Integral());
    h_barrel_MC_nume_50to70[_WW]                  ->Scale(5.1702e+04 / h_barrel_MC_nume_50to70[_WW]                 ->Integral());
    h_barrel_MC_nume_50to70[_WZ]                  ->Scale(1.7269e+04 / h_barrel_MC_nume_50to70[_WZ]                 ->Integral());
    h_barrel_MC_nume_50to70[_ZZ]                  ->Scale(4.7152e+03 / h_barrel_MC_nume_50to70[_ZZ]                 ->Integral());
    h_barrel_MC_nume_50to70[_tW]                  ->Scale(2.4993e+04 / h_barrel_MC_nume_50to70[_tW]                 ->Integral());
    h_barrel_MC_nume_50to70[_tbarW]               ->Scale(2.5035e+04 / h_barrel_MC_nume_50to70[_tbarW]              ->Integral());
    h_barrel_MC_nume_50to70[_ttbar]               ->Scale(5.8579e+05 / h_barrel_MC_nume_50to70[_ttbar]              ->Integral());
    h_barrel_MC_nume_50to70[_WJets]               ->Scale(9.5136e+06 / h_barrel_MC_nume_50to70[_WJets]              ->Integral());
    h_barrel_MC_nume_50to70[_QCDMuEnriched_Full]  ->Scale(1.0364e+06 / h_barrel_MC_nume_50to70[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_nume_70to100[_DY_Full]            ->Scale(1.1070e+06 / h_barrel_MC_nume_70to100[_DY_Full]            ->Integral());
    h_barrel_MC_nume_70to100[_WW]                 ->Scale(3.3367e+04 / h_barrel_MC_nume_70to100[_WW]                 ->Integral());
    h_barrel_MC_nume_70to100[_WZ]                 ->Scale(1.1559e+04 / h_barrel_MC_nume_70to100[_WZ]                 ->Integral());
    h_barrel_MC_nume_70to100[_ZZ]                 ->Scale(3.3983e+03 / h_barrel_MC_nume_70to100[_ZZ]                 ->Integral());
    h_barrel_MC_nume_70to100[_tW]                 ->Scale(2.3878e+04 / h_barrel_MC_nume_70to100[_tW]                 ->Integral());
    h_barrel_MC_nume_70to100[_tbarW]              ->Scale(2.3798e+04 / h_barrel_MC_nume_70to100[_tbarW]              ->Integral());
    h_barrel_MC_nume_70to100[_ttbar]              ->Scale(5.2500e+05 / h_barrel_MC_nume_70to100[_ttbar]              ->Integral());
    h_barrel_MC_nume_70to100[_WJets]              ->Scale(3.1498e+06 / h_barrel_MC_nume_70to100[_WJets]              ->Integral());
    h_barrel_MC_nume_70to100[_QCDMuEnriched_Full] ->Scale(2.6953e+05 / h_barrel_MC_nume_70to100[_QCDMuEnriched_Full] ->Integral());
    h_barrel_MC_nume_100to500[_DY_Full]           ->Scale(3.9714e+05 / h_barrel_MC_nume_100to500[_DY_Full]            ->Integral());
    h_barrel_MC_nume_100to500[_WW]                ->Scale(1.3940e+04 / h_barrel_MC_nume_100to500[_WW]                 ->Integral());
    h_barrel_MC_nume_100to500[_WZ]                ->Scale(4.6703e+03 / h_barrel_MC_nume_100to500[_WZ]                 ->Integral());
    h_barrel_MC_nume_100to500[_ZZ]                ->Scale(1.3270e+03 / h_barrel_MC_nume_100to500[_ZZ]                 ->Integral());
    h_barrel_MC_nume_100to500[_tW]                ->Scale(2.2911e+04 / h_barrel_MC_nume_100to500[_tW]                 ->Integral());
    h_barrel_MC_nume_100to500[_tbarW]             ->Scale(2.2654e+04 / h_barrel_MC_nume_100to500[_tbarW]              ->Integral());
    h_barrel_MC_nume_100to500[_ttbar]             ->Scale(3.6522e+05 / h_barrel_MC_nume_100to500[_ttbar]              ->Integral());
    h_barrel_MC_nume_100to500[_WJets]             ->Scale(1.2361e+06 / h_barrel_MC_nume_100to500[_WJets]              ->Integral());
    h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]->Scale(5.8852e+04 / h_barrel_MC_nume_100to500[_QCDMuEnriched_Full] ->Integral());

    h_endcap_MC_nume_50to70[_DY_Full]             ->Scale(1.9068e+06 / h_endcap_MC_nume_50to70[_DY_Full]            ->Integral());
    h_endcap_MC_nume_50to70[_WW]                  ->Scale(3.9296e+04 / h_endcap_MC_nume_50to70[_WW]                 ->Integral());
    h_endcap_MC_nume_50to70[_WZ]                  ->Scale(1.2833e+04 / h_endcap_MC_nume_50to70[_WZ]                 ->Integral());
    h_endcap_MC_nume_50to70[_ZZ]                  ->Scale(3.6118e+03 / h_endcap_MC_nume_50to70[_ZZ]                 ->Integral());
    h_endcap_MC_nume_50to70[_tW]                  ->Scale(1.1300e+04 / h_endcap_MC_nume_50to70[_tW]                 ->Integral());
    h_endcap_MC_nume_50to70[_tbarW]               ->Scale(1.1271e+04 / h_endcap_MC_nume_50to70[_tbarW]              ->Integral());
    h_endcap_MC_nume_50to70[_ttbar]               ->Scale(2.7070e+05 / h_endcap_MC_nume_50to70[_ttbar]              ->Integral());
    h_endcap_MC_nume_50to70[_WJets]               ->Scale(7.6807e+06 / h_endcap_MC_nume_50to70[_WJets]              ->Integral());
    h_endcap_MC_nume_50to70[_QCDMuEnriched_Full]  ->Scale(1.1949e+06 / h_endcap_MC_nume_50to70[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_nume_70to100[_DY_Full]            ->Scale(6.1464e+05 / h_endcap_MC_nume_70to100[_DY_Full]            ->Integral());
    h_endcap_MC_nume_70to100[_WW]                 ->Scale(2.4375e+04 / h_endcap_MC_nume_70to100[_WW]                 ->Integral());
    h_endcap_MC_nume_70to100[_WZ]                 ->Scale(8.2399e+03 / h_endcap_MC_nume_70to100[_WZ]                 ->Integral());
    h_endcap_MC_nume_70to100[_ZZ]                 ->Scale(2.5715e+03 / h_endcap_MC_nume_70to100[_ZZ]                 ->Integral());
    h_endcap_MC_nume_70to100[_tW]                 ->Scale(1.0180e+04 / h_endcap_MC_nume_70to100[_tW]                 ->Integral());
    h_endcap_MC_nume_70to100[_tbarW]              ->Scale(1.0047e+04 / h_endcap_MC_nume_70to100[_tbarW]              ->Integral());
    h_endcap_MC_nume_70to100[_ttbar]              ->Scale(2.3196e+05 / h_endcap_MC_nume_70to100[_ttbar]              ->Integral());
    h_endcap_MC_nume_70to100[_WJets]              ->Scale(2.7487e+06 / h_endcap_MC_nume_70to100[_WJets]              ->Integral());
    h_endcap_MC_nume_70to100[_QCDMuEnriched_Full] ->Scale(2.8630e+05 / h_endcap_MC_nume_70to100[_QCDMuEnriched_Full] ->Integral());
    h_endcap_MC_nume_100to500[_DY_Full]           ->Scale(2.9622e+05 / h_endcap_MC_nume_100to500[_DY_Full]            ->Integral());
    h_endcap_MC_nume_100to500[_WW]                ->Scale(8.9375e+03 / h_endcap_MC_nume_100to500[_WW]                 ->Integral());
    h_endcap_MC_nume_100to500[_WZ]                ->Scale(4.2826e+03 / h_endcap_MC_nume_100to500[_WZ]                 ->Integral());
    h_endcap_MC_nume_100to500[_ZZ]                ->Scale(1.1298e+03 / h_endcap_MC_nume_100to500[_ZZ]                 ->Integral());
    h_endcap_MC_nume_100to500[_tW]                ->Scale(7.9142e+03 / h_endcap_MC_nume_100to500[_tW]                 ->Integral());
    h_endcap_MC_nume_100to500[_tbarW]             ->Scale(7.9509e+03 / h_endcap_MC_nume_100to500[_tbarW]              ->Integral());
    h_endcap_MC_nume_100to500[_ttbar]             ->Scale(1.4191e+05 / h_endcap_MC_nume_100to500[_ttbar]              ->Integral());
    h_endcap_MC_nume_100to500[_WJets]             ->Scale(8.8404e+05 / h_endcap_MC_nume_100to500[_WJets]              ->Integral());
    h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]->Scale(6.0482e+04 / h_endcap_MC_nume_100to500[_QCDMuEnriched_Full] ->Integral());

    // barrel ctrlminator
    THStack * s_barrel_ctrl = new THStack("s_barrel_ctrl", "");
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_WW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_WW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_WW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_WZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_WZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_WZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_ZZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_ZZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_ZZ]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_tbarW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_tbarW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_tbarW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_tW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_tW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_tW]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_ttbar]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_ttbar]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_ttbar]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_WJets]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_WJets]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_WJets]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_DY_Full]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_DY_Full]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_DY_Full]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_50to70[_QCDMuEnriched_Full]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_70to100[_QCDMuEnriched_Full]);
    s_barrel_ctrl->Add(h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]);

    h_barrel_data_ctrl_50to70->Add(h_barrel_data_ctrl_70to100);
    h_barrel_data_ctrl_50to70->Add(h_barrel_data_ctrl_100to500);

    // endcap ctrlminator
    THStack * s_endcap_ctrl = new THStack("s_endcap_ctrl", "");
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_WW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_WW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_WW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_WZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_WZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_WZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_ZZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_ZZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_ZZ]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_tbarW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_tbarW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_tbarW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_tW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_tW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_tW]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_ttbar]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_ttbar]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_ttbar]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_WJets]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_WJets]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_WJets]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_DY_Full]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_DY_Full]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_DY_Full]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_50to70[_QCDMuEnriched_Full]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_70to100[_QCDMuEnriched_Full]);
    s_endcap_ctrl->Add(h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]);

    h_endcap_data_ctrl_50to70->Add(h_endcap_data_ctrl_70to100);
    h_endcap_data_ctrl_50to70->Add(h_endcap_data_ctrl_100to500);

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
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_tbarW]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_tW]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_ttbar]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_WJets]);
    s_barrel_nume->Add(h_barrel_MC_nume_50to70[_DY_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_70to100[_DY_Full]);
    s_barrel_nume->Add(h_barrel_MC_nume_100to500[_DY_Full]);
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
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_tbarW]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_tW]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_ttbar]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_WJets]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_DY_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_50to70[_QCDMuEnriched_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_70to100[_QCDMuEnriched_Full]);
    s_endcap_nume->Add(h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]);

    h_endcap_data_nume_50to70->Add(h_endcap_data_nume_70to100);
    h_endcap_data_nume_50to70->Add(h_endcap_data_nume_100to500);

    // barrel denominator
    THStack * s_barrel_deno = new THStack("s_barrel_deno", "");
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_WW]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_WW]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_WW]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_tW]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_tW]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_tW]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_nume_50to70[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_nume_70to100[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_nume_100to500[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_WW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_WW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_WW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_WZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_ZZ]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_tbarW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_tW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_tW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_tW]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_ttbar]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_WJets]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_DY_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_50to70[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_70to100[_QCDMuEnriched_Full]);
    s_barrel_deno->Add(h_barrel_MC_ctrl_100to500[_QCDMuEnriched_Full]);

    TH1D *h_barrel_data_deno = ((TH1D*)(h_barrel_data_nume_50to70->Clone("h_barrel_data_deno")));
    h_barrel_data_deno->Add(h_barrel_data_ctrl_50to70);

    // endcap numerator
    THStack * s_endcap_deno = new THStack("s_endcap_deno", "");
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_WW]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_WW]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_WW]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_tW]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_tW]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_tW]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_nume_50to70[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_nume_70to100[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_nume_100to500[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_WW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_WW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_WW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_WZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_ZZ]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_tbarW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_tW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_tW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_tW]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_ttbar]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_WJets]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_DY_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_50to70[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_70to100[_QCDMuEnriched_Full]);
    s_endcap_deno->Add(h_endcap_MC_ctrl_100to500[_QCDMuEnriched_Full]);

    TH1D *h_endcap_data_deno = ((TH1D*)(h_endcap_data_nume_50to70->Clone("h_endcap_data_deno")));
    h_endcap_data_deno->Add(h_endcap_data_ctrl_50to70);

// ------------------------------------- DRAWING -----------------------------------------------------
    myRatioPlot_t *RP_barrel_ctrl = new myRatioPlot_t("RP_barrel_ctrl", s_barrel_ctrl, h_barrel_data_ctrl_50to70);
    myRatioPlot_t *RP_endcap_ctrl = new myRatioPlot_t("RP_endcap_ctrl", s_endcap_ctrl, h_endcap_data_ctrl_50to70);
    myRatioPlot_t *RP_barrel_nume = new myRatioPlot_t("RP_barrel_nume", s_barrel_nume, h_barrel_data_nume_50to70);
    myRatioPlot_t *RP_endcap_nume = new myRatioPlot_t("RP_endcap_nume", s_endcap_nume, h_endcap_data_nume_50to70);
    myRatioPlot_t *RP_barrel_deno = new myRatioPlot_t("RP_barrel_deno", s_barrel_deno, h_barrel_data_deno);
    myRatioPlot_t *RP_endcap_deno = new myRatioPlot_t("RP_endcap_deno", s_endcap_deno, h_endcap_data_deno);
    RP_barrel_ctrl->SetPlots("p_{T} (#mu_{barrel}^{ctrl})", 0, 5);
    RP_endcap_ctrl->SetPlots("p_{T} (#mu_{endcap}^{ctrl})", 0, 5);
    RP_barrel_nume->SetPlots("p_{T} (#mu_{barrel}^{nume})", 0, 5);
    RP_endcap_nume->SetPlots("p_{T} (#mu_{endcap}^{nume})", 0, 5);
    RP_barrel_deno->SetPlots("p_{T} (#mu_{barrel}^{deno})", 0, 5);
    RP_endcap_deno->SetPlots("p_{T} (#mu_{endcap}^{deno})", 0, 5);

    TLegend *legend = new TLegend(0.5, 0.65, 0.95, 0.95);

    legend->AddEntry(h_barrel_data_ctrl_50to70, "Data", "lp");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_DY_Full], "DY", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_barrel_MC_ctrl_50to70[_QCDMuEnriched_Full], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend->SetNColumns(2);

    RP_barrel_ctrl->ImportLegend(legend);
    RP_endcap_ctrl->ImportLegend(legend);
    RP_barrel_nume->ImportLegend(legend);
    RP_endcap_nume->ImportLegend(legend);
    RP_barrel_deno->ImportLegend(legend);
    RP_endcap_deno->ImportLegend(legend);
    RP_barrel_ctrl->Draw(1, 1e8, 0);
    RP_endcap_ctrl->Draw(1, 1e8, 0);
    RP_barrel_nume->Draw(1, 1e8, 0);
    RP_endcap_nume->Draw(1, 1e8, 0);
    RP_barrel_deno->Draw(1, 1e8, 0);
    RP_endcap_deno->Draw(1, 1e8, 0);

    Double_t int_data_barrel, int_data_endcap, err_data_barrel, err_data_endcap, int_mc_barrel, int_mc_endcap, err_mc_barrel, err_mc_endcap;
    int_data_barrel = h_barrel_data_deno->IntegralAndError(1, h_barrel_data_deno->GetSize()-2, err_data_barrel);
    int_data_endcap = h_endcap_data_deno->IntegralAndError(1, h_endcap_data_deno->GetSize()-2, err_data_endcap);
    int_mc_barrel = ((TH1D*)(s_barrel_deno->GetStack()->Last()))->IntegralAndError(1, h_barrel_data_deno->GetSize()-2, err_mc_barrel);
    int_mc_endcap = ((TH1D*)(s_endcap_deno->GetStack()->Last()))->IntegralAndError(1, h_endcap_data_deno->GetSize()-2, err_mc_endcap);

    cout << "Data integral in barrel: " << int_data_barrel << "+-" << err_data_barrel << endl;
    cout << "Data integral in endcap: " << int_data_endcap << "+-" << err_data_endcap << endl;
    cout << "MC integral in barrel: " << int_mc_barrel << "+-" << err_mc_barrel << endl;
    cout << "MC integral in endcap: " << int_mc_endcap << "+-" << err_mc_endcap << endl;
} // End of Fit_HistDrawer()



///----------------------- ESTIMATIONS ------------------------- ///
void E_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";

    f = new TFile(Dir+"QCDest_E.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special],
         *h_mass_SS[_EndOf_Data_Special],
         *h_mass_temp[_EndOf_Data_Special],
         *h_mass_test[_EndOf_Data_Special],
         *h_pT[_EndOf_Data_Special], *h_rapi[_EndOf_Data_Special],
         *h_PFiso_barrel[_EndOf_Data_Special],
         *h_TrkIso_barrel[_EndOf_Data_Special],
         *h_ECALiso_barrel[_EndOf_Data_Special],
         *h_HCALiso_barrel[_EndOf_Data_Special],
         *h_PFiso_endcap[_EndOf_Data_Special],
         *h_TrkIso_endcap[_EndOf_Data_Special],
         *h_ECALiso_endcap[_EndOf_Data_Special],
         *h_HCALiso_endcap[_EndOf_Data_Special];
    TH1D *h_QCD_est, *h_QCD_est_SS, *h_QCD_est_temp;
    THStack *s_mass_wQCD = new THStack("s_mass_wQCD", "");
    THStack *s_mass_woQCD = new THStack("s_mass_woQCD_SS", "");
    THStack *s_mass_wQCD_SS = new THStack("s_mass_wQCD_SS", "");
    THStack *s_mass_woQCD_SS = new THStack("s_mass_woQCD", "");
    THStack *s_mass_wQCD_temp = new THStack("s_mass_wQCD_temp", "");
    THStack *s_mass_woQCD_temp = new THStack("s_mass_woQCD_temp", "");
    THStack *s_mass_wQCD_test = new THStack("s_mass_wQCD_test", "");
    THStack *s_mass_woQCD_test = new THStack("s_mass_woQCD_test", "");
    THStack *s_pT = new THStack("s_pT", "");
    THStack *s_rapi = new THStack("s_rapi", "");
    THStack *s_PFiso_barrel = new THStack("s_PFiso_barrel", "");
    THStack *s_TrkIso_barrel = new THStack("s_TrkIso_barrel", "");
    THStack *s_ECALiso_barrel = new THStack("s_ECALiso_barrel", "");
    THStack *s_HCALiso_barrel = new THStack("s_HCALiso_barrel", "");
    THStack *s_PFiso_endcap = new THStack("s_PFiso_endcap", "");
    THStack *s_TrkIso_endcap = new THStack("s_TrkIso_endcap", "");
    THStack *s_ECALiso_endcap = new THStack("s_ECALiso_endcap", "");
    THStack *s_HCALiso_endcap = new THStack("s_HCALiso_endcap", "");
    Color_t color = kBlack;

    // Loop over all processes (adding all histograms)
    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        f->GetObject("h_mass_template_"+Mgr.Procname[pr], h_mass_temp[pr]);
        f->GetObject("h_mass_test_"+Mgr.Procname[pr], h_mass_test[pr]);
        f->GetObject("h_pT_"+Mgr.Procname[pr], h_pT[pr]);
        f->GetObject("h_rapi_"+Mgr.Procname[pr], h_rapi[pr]);
        f->GetObject("h_PFiso_Rho_barrel_ctrl_"+Mgr.Procname[pr], h_PFiso_barrel[pr]);
        f->GetObject("h_TrkIso_barrel_ctrl_"+Mgr.Procname[pr], h_TrkIso_barrel[pr]);
        f->GetObject("h_ECALiso_barrel_ctrl_"+Mgr.Procname[pr], h_ECALiso_barrel[pr]);
        f->GetObject("h_HCALiso_barrel_ctrl_"+Mgr.Procname[pr], h_HCALiso_barrel[pr]);
        f->GetObject("h_PFiso_Rho_endcap_ctrl_"+Mgr.Procname[pr], h_PFiso_endcap[pr]);
        f->GetObject("h_TrkIso_endcap_ctrl_"+Mgr.Procname[pr], h_TrkIso_endcap[pr]);
        f->GetObject("h_ECALiso_endcap_ctrl_"+Mgr.Procname[pr], h_ECALiso_endcap[pr]);
        f->GetObject("h_HCALiso_endcap_ctrl_"+Mgr.Procname[pr], h_HCALiso_endcap[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);
        h_mass_temp[pr]->SetDirectory(0);
        h_mass_test[pr]->SetDirectory(0);
        h_pT[pr]->SetDirectory(0);
        h_rapi[pr]->SetDirectory(0);
        h_PFiso_barrel[pr]->SetDirectory(0);
        h_TrkIso_barrel[pr]->SetDirectory(0);
        h_ECALiso_barrel[pr]->SetDirectory(0);
        h_HCALiso_barrel[pr]->SetDirectory(0);
        h_PFiso_endcap[pr]->SetDirectory(0);
        h_TrkIso_endcap[pr]->SetDirectory(0);
        h_ECALiso_endcap[pr]->SetDirectory(0);
        h_HCALiso_endcap[pr]->SetDirectory(0);


        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
            removeNegativeBins(h_mass_temp[pr]);
            removeNegativeBins(h_mass_test[pr]);
            removeNegativeBins(h_pT[pr]);
            removeNegativeBins(h_rapi[pr]);
            removeNegativeBins(h_PFiso_barrel[pr]);
            removeNegativeBins(h_TrkIso_barrel[pr]);
            removeNegativeBins(h_ECALiso_barrel[pr]);
            removeNegativeBins(h_HCALiso_barrel[pr]);
            removeNegativeBins(h_PFiso_endcap[pr]);
            removeNegativeBins(h_TrkIso_endcap[pr]);
            removeNegativeBins(h_ECALiso_endcap[pr]);
            removeNegativeBins(h_HCALiso_endcap[pr]);
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
        else if (pr < _EndOf_GJets_Normal) color = kYellow + 3;

        if (pr < _DoubleEG_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
            h_mass_temp[pr]->SetFillColor(color);
            h_mass_temp[pr]->SetLineColor(color);
            h_mass_test[pr]->SetFillColor(color);
            h_mass_test[pr]->SetLineColor(color);
            h_pT[pr]->SetFillColor(color);
            h_pT[pr]->SetLineColor(color);
            h_rapi[pr]->SetFillColor(color);
            h_rapi[pr]->SetLineColor(color);
            h_PFiso_barrel[pr]->SetFillColor(color);
            h_TrkIso_barrel[pr]->SetFillColor(color);
            h_ECALiso_barrel[pr]->SetFillColor(color);
            h_HCALiso_barrel[pr]->SetFillColor(color);
            h_PFiso_endcap[pr]->SetFillColor(color);
            h_TrkIso_endcap[pr]->SetFillColor(color);
            h_ECALiso_endcap[pr]->SetFillColor(color);
            h_HCALiso_endcap[pr]->SetFillColor(color);
            h_PFiso_barrel[pr]->SetLineColor(color);
            h_TrkIso_barrel[pr]->SetLineColor(color);
            h_ECALiso_barrel[pr]->SetLineColor(color);
            h_HCALiso_barrel[pr]->SetLineColor(color);
            h_PFiso_endcap[pr]->SetLineColor(color);
            h_TrkIso_endcap[pr]->SetLineColor(color);
            h_ECALiso_endcap[pr]->SetLineColor(color);
            h_HCALiso_endcap[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
            h_mass_temp[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_temp[pr]->SetMarkerColor(kBlack);
            h_mass_temp[pr]->SetLineColor(kBlack);
            h_mass_test[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test[pr]->SetMarkerColor(kBlack);
            h_mass_test[pr]->SetLineColor(kBlack);
            h_pT[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT[pr]->SetMarkerColor(kBlack);
            h_pT[pr]->SetLineColor(kBlack);
            h_rapi[pr]->SetMarkerStyle(kFullDotLarge);
            h_rapi[pr]->SetMarkerColor(kBlack);
            h_rapi[pr]->SetLineColor(kBlack);
            h_PFiso_barrel[pr]->SetMarkerStyle(kFullDotLarge);
            h_TrkIso_barrel[pr]->SetMarkerStyle(kFullDotLarge);
            h_ECALiso_barrel[pr]->SetMarkerStyle(kFullDotLarge);
            h_HCALiso_barrel[pr]->SetMarkerStyle(kFullDotLarge);
            h_PFiso_endcap[pr]->SetMarkerStyle(kFullDotLarge);
            h_TrkIso_endcap[pr]->SetMarkerStyle(kFullDotLarge);
            h_ECALiso_endcap[pr]->SetMarkerStyle(kFullDotLarge);
            h_HCALiso_endcap[pr]->SetMarkerStyle(kFullDotLarge);
            h_PFiso_barrel[pr]->SetMarkerColor(kBlack);
            h_TrkIso_barrel[pr]->SetMarkerColor(kBlack);
            h_ECALiso_barrel[pr]->SetMarkerColor(kBlack);
            h_HCALiso_barrel[pr]->SetMarkerColor(kBlack);
            h_PFiso_endcap[pr]->SetMarkerColor(kBlack);
            h_TrkIso_endcap[pr]->SetMarkerColor(kBlack);
            h_ECALiso_endcap[pr]->SetMarkerColor(kBlack);
            h_HCALiso_endcap[pr]->SetMarkerColor(kBlack);
            h_PFiso_barrel[pr]->SetLineColor(kBlack);
            h_TrkIso_barrel[pr]->SetLineColor(kBlack);
            h_ECALiso_barrel[pr]->SetLineColor(kBlack);
            h_HCALiso_barrel[pr]->SetLineColor(kBlack);
            h_PFiso_endcap[pr]->SetLineColor(kBlack);
            h_TrkIso_endcap[pr]->SetLineColor(kBlack);
            h_ECALiso_endcap[pr]->SetLineColor(kBlack);
            h_HCALiso_endcap[pr]->SetLineColor(kBlack);
        }

        // Adding hists to THStacks
        if (pr < _DoubleEG_B)
        {
            s_mass_wQCD->Add(h_mass[pr]);
            s_mass_wQCD_SS->Add(h_mass_SS[pr]);
            s_mass_wQCD_temp->Add(h_mass_temp[pr]);
            s_mass_wQCD_test->Add(h_mass_test[pr]);
            s_pT->Add(h_pT[pr]);
            s_rapi->Add(h_rapi[pr]);            
        }
        if (pr < _QCDEMEnriched_20to30 /*|| (pr > _EndOf_QCDEMEnriched_Normal && pr < _DoubleEG_B)*/)
        {
            s_mass_woQCD->Add(h_mass[pr]);
            s_mass_woQCD_SS->Add(h_mass_SS[pr]);
            s_mass_woQCD_temp->Add(h_mass_temp[pr]);
            s_mass_woQCD_test->Add(h_mass_test[pr]);
        }

        // Adding up for convenience
        if (pr == _DY_10to50)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
            h_mass_temp[_DY_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DY_Full")));
            h_mass_temp[_DY_Full]->SetDirectory(0);
            h_mass_test[_DY_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DY_Full")));
            h_mass_test[_DY_Full]->SetDirectory(0);
            h_pT[_DY_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_DY_Full")));
            h_pT[_DY_Full]->SetDirectory(0);
            h_rapi[_DY_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_DY_Full")));
            h_rapi[_DY_Full]->SetDirectory(0);
            h_PFiso_barrel[_DY_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_DY_Full")));
            h_TrkIso_barrel[_DY_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_DY_Full")));
            h_ECALiso_barrel[_DY_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_DY_Full")));
            h_HCALiso_barrel[_DY_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_DY_Full")));
            h_PFiso_endcap[_DY_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_DY_Full")));
            h_TrkIso_endcap[_DY_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_DY_Full")));
            h_ECALiso_endcap[_DY_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_DY_Full")));
            h_HCALiso_endcap[_DY_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_DY_Full")));
            h_PFiso_barrel[_DY_Full]->SetDirectory(0);
            h_TrkIso_barrel[_DY_Full]->SetDirectory(0);
            h_ECALiso_barrel[_DY_Full]->SetDirectory(0);
            h_HCALiso_barrel[_DY_Full]->SetDirectory(0);
            h_PFiso_endcap[_DY_Full]->SetDirectory(0);
            h_TrkIso_endcap[_DY_Full]->SetDirectory(0);
            h_ECALiso_endcap[_DY_Full]->SetDirectory(0);
            h_HCALiso_endcap[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DY_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DY_Full]->Add(h_mass_test[pr]);
            h_pT[_DY_Full]->Add(h_pT[pr]);
            h_rapi[_DY_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_DY_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_DY_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_DY_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_DY_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_DY_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_DY_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_DY_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_DY_Full]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
            h_mass_temp[_ttbar_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_ttbar_Full")));
            h_mass_temp[_ttbar_Full]->SetDirectory(0);
            h_mass_test[_ttbar_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_ttbar_Full")));
            h_mass_test[_ttbar_Full]->SetDirectory(0);
            h_pT[_ttbar_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_ttbar_Full")));
            h_pT[_ttbar_Full]->SetDirectory(0);
            h_rapi[_ttbar_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_ttbar_Full")));
            h_rapi[_ttbar_Full]->SetDirectory(0);
            h_PFiso_barrel[_ttbar_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_ttbar_Full")));
            h_TrkIso_barrel[_ttbar_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_ttbar_Full")));
            h_ECALiso_barrel[_ttbar_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_ttbar_Full")));
            h_HCALiso_barrel[_ttbar_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_ttbar_Full")));
            h_PFiso_endcap[_ttbar_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_ttbar_Full")));
            h_TrkIso_endcap[_ttbar_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_ttbar_Full")));
            h_ECALiso_endcap[_ttbar_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_ttbar_Full")));
            h_HCALiso_endcap[_ttbar_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_ttbar_Full")));
            h_PFiso_barrel[_ttbar_Full]->SetDirectory(0);
            h_TrkIso_barrel[_ttbar_Full]->SetDirectory(0);
            h_ECALiso_barrel[_ttbar_Full]->SetDirectory(0);
            h_HCALiso_barrel[_ttbar_Full]->SetDirectory(0);
            h_PFiso_endcap[_ttbar_Full]->SetDirectory(0);
            h_TrkIso_endcap[_ttbar_Full]->SetDirectory(0);
            h_ECALiso_endcap[_ttbar_Full]->SetDirectory(0);
            h_HCALiso_endcap[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_ttbar_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_ttbar_Full]->Add(h_mass_test[pr]);
            h_pT[_ttbar_Full]->Add(h_pT[pr]);
            h_rapi[_ttbar_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_ttbar_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_ttbar_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_ttbar_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_ttbar_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_ttbar_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_ttbar_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_ttbar_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_ttbar_Full]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
            h_mass_temp[_VVnST] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_VVnST")));
            h_mass_temp[_VVnST]->SetDirectory(0);
            h_mass_test[_VVnST] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_VVnST")));
            h_mass_test[_VVnST]->SetDirectory(0);
            h_pT[_VVnST] = ((TH1D*)(h_pT[pr]->Clone("h_pT_VVnST")));
            h_pT[_VVnST]->SetDirectory(0);
            h_rapi[_VVnST] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_VVnST")));
            h_rapi[_VVnST]->SetDirectory(0);
            h_PFiso_barrel[_VVnST] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_VVnST")));
            h_TrkIso_barrel[_VVnST] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_VVnST")));
            h_ECALiso_barrel[_VVnST] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_VVnST")));
            h_HCALiso_barrel[_VVnST] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_VVnST")));
            h_PFiso_endcap[_VVnST] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_VVnST")));
            h_TrkIso_endcap[_VVnST] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_VVnST")));
            h_ECALiso_endcap[_VVnST] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_VVnST")));
            h_HCALiso_endcap[_VVnST] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_VVnST")));
            h_PFiso_barrel[_VVnST]->SetDirectory(0);
            h_TrkIso_barrel[_VVnST]->SetDirectory(0);
            h_ECALiso_barrel[_VVnST]->SetDirectory(0);
            h_HCALiso_barrel[_VVnST]->SetDirectory(0);
            h_PFiso_endcap[_VVnST]->SetDirectory(0);
            h_TrkIso_endcap[_VVnST]->SetDirectory(0);
            h_ECALiso_endcap[_VVnST]->SetDirectory(0);
            h_HCALiso_endcap[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
            h_mass_temp[_VVnST]->Add(h_mass_temp[pr]);
            h_mass_test[_VVnST]->Add(h_mass_test[pr]);
            h_pT[_VVnST]->Add(h_pT[pr]);
            h_rapi[_VVnST]->Add(h_rapi[pr]);
            h_PFiso_barrel[_VVnST]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_VVnST]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_VVnST]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_VVnST]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_VVnST]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_VVnST]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_VVnST]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_VVnST]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
            h_mass_temp[_WJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_WJets_Full")));
            h_mass_temp[_WJets_Full]->SetDirectory(0);
            h_mass_test[_WJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_WJets_Full")));
            h_mass_test[_WJets_Full]->SetDirectory(0);
            h_pT[_WJets_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_WJets_Full")));
            h_pT[_WJets_Full]->SetDirectory(0);
            h_rapi[_WJets_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_WJets_Full")));
            h_rapi[_WJets_Full]->SetDirectory(0);
            h_PFiso_barrel[_WJets_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_WJets_Full")));
            h_TrkIso_barrel[_WJets_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_WJets_Full")));
            h_ECALiso_barrel[_WJets_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_WJets_Full")));
            h_HCALiso_barrel[_WJets_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_WJets_Full")));
            h_PFiso_endcap[_WJets_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_WJets_Full")));
            h_TrkIso_endcap[_WJets_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_WJets_Full")));
            h_ECALiso_endcap[_WJets_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_WJets_Full")));
            h_HCALiso_endcap[_WJets_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_WJets_Full")));
            h_PFiso_barrel[_WJets_Full]->SetDirectory(0);
            h_TrkIso_barrel[_WJets_Full]->SetDirectory(0);
            h_ECALiso_barrel[_WJets_Full]->SetDirectory(0);
            h_HCALiso_barrel[_WJets_Full]->SetDirectory(0);
            h_PFiso_endcap[_WJets_Full]->SetDirectory(0);
            h_TrkIso_endcap[_WJets_Full]->SetDirectory(0);
            h_ECALiso_endcap[_WJets_Full]->SetDirectory(0);
            h_HCALiso_endcap[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_WJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_WJets_Full]->Add(h_mass_test[pr]);
            h_pT[_WJets_Full]->Add(h_pT[pr]);
            h_rapi[_WJets_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_WJets_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_WJets_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_WJets_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_WJets_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_WJets_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_WJets_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_WJets_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_WJets_Full]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _QCDEMEnriched_20to30)
        {
            h_mass[_QCDEMEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDEMEnriched_Full")));
            h_mass[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_SS[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_QCDEMEnriched_Full")));
            h_mass_SS[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_temp[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_QCDEMEnriched_Full")));
            h_mass_temp[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_test[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_QCDEMEnriched_Full")));
            h_mass_test[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT[_QCDEMEnriched_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_QCDEMEnriched_Full")));
            h_pT[_QCDEMEnriched_Full]->SetDirectory(0);
            h_rapi[_QCDEMEnriched_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_QCDEMEnriched_Full")));
            h_rapi[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_barrel[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_QCDEMEnriched_Full")));
            h_TrkIso_barrel[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_QCDEMEnriched_Full")));
            h_ECALiso_barrel[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_QCDEMEnriched_Full")));
            h_HCALiso_barrel[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_QCDEMEnriched_Full")));
            h_PFiso_endcap[_QCDEMEnriched_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_QCDEMEnriched_Full")));
            h_TrkIso_endcap[_QCDEMEnriched_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_QCDEMEnriched_Full")));
            h_ECALiso_endcap[_QCDEMEnriched_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_QCDEMEnriched_Full")));
            h_HCALiso_endcap[_QCDEMEnriched_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_QCDEMEnriched_Full")));
            h_PFiso_barrel[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_barrel[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_barrel[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_barrel[_QCDEMEnriched_Full]->SetDirectory(0);
            h_PFiso_endcap[_QCDEMEnriched_Full]->SetDirectory(0);
            h_TrkIso_endcap[_QCDEMEnriched_Full]->SetDirectory(0);
            h_ECALiso_endcap[_QCDEMEnriched_Full]->SetDirectory(0);
            h_HCALiso_endcap[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_QCDEMEnriched_Normal)
        {
            h_mass[_QCDEMEnriched_Full]->Add(h_mass[pr]);
            h_mass_SS[_QCDEMEnriched_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_QCDEMEnriched_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_QCDEMEnriched_Full]->Add(h_mass_test[pr]);
            h_pT[_QCDEMEnriched_Full]->Add(h_pT[pr]);
            h_rapi[_QCDEMEnriched_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_QCDEMEnriched_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_QCDEMEnriched_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_QCDEMEnriched_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_QCDEMEnriched_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_QCDEMEnriched_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_QCDEMEnriched_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_QCDEMEnriched_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_QCDEMEnriched_Full]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _GJets_20to100)
        {
            h_mass[_GJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_GJets_Full")));
            h_mass[_GJets_Full]->SetDirectory(0);
            h_mass_SS[_GJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_GJets_Full")));
            h_mass_SS[_GJets_Full]->SetDirectory(0);
            h_mass_temp[_GJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_GJets_Full")));
            h_mass_temp[_GJets_Full]->SetDirectory(0);
            h_mass_test[_GJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_GJets_Full")));
            h_mass_test[_GJets_Full]->SetDirectory(0);
            h_pT[_GJets_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_GJets_Full")));
            h_pT[_GJets_Full]->SetDirectory(0);
            h_rapi[_GJets_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_GJets_Full")));
            h_rapi[_GJets_Full]->SetDirectory(0);
            h_PFiso_barrel[_GJets_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_GJets_Full")));
            h_TrkIso_barrel[_GJets_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_GJets_Full")));
            h_ECALiso_barrel[_GJets_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_GJets_Full")));
            h_HCALiso_barrel[_GJets_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_GJets_Full")));
            h_PFiso_endcap[_GJets_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_GJets_Full")));
            h_TrkIso_endcap[_GJets_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_GJets_Full")));
            h_ECALiso_endcap[_GJets_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_GJets_Full")));
            h_HCALiso_endcap[_GJets_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_GJets_Full")));
            h_PFiso_barrel[_GJets_Full]->SetDirectory(0);
            h_TrkIso_barrel[_GJets_Full]->SetDirectory(0);
            h_ECALiso_barrel[_GJets_Full]->SetDirectory(0);
            h_HCALiso_barrel[_GJets_Full]->SetDirectory(0);
            h_PFiso_endcap[_GJets_Full]->SetDirectory(0);
            h_TrkIso_endcap[_GJets_Full]->SetDirectory(0);
            h_ECALiso_endcap[_GJets_Full]->SetDirectory(0);
            h_HCALiso_endcap[_GJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_GJets_Normal)
        {
            h_mass[_GJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_GJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_GJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_GJets_Full]->Add(h_mass_test[pr]);
            h_pT[_GJets_Full]->Add(h_pT[pr]);
            h_rapi[_GJets_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_GJets_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_GJets_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_GJets_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_GJets_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_GJets_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_GJets_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_GJets_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_GJets_Full]->Add(h_HCALiso_endcap[pr]);
        }
        else if (pr == _DoubleEG_B)
        {
            h_mass[_DoubleEG_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DoubleEG_Full")));
            h_mass[_DoubleEG_Full]->SetDirectory(0);
            h_mass_SS[_DoubleEG_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DoubleEG_Full")));
            h_mass_SS[_DoubleEG_Full]->SetDirectory(0);
            h_mass_temp[_DoubleEG_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DoubleEG_Full")));
            h_mass_temp[_DoubleEG_Full]->SetDirectory(0);
            h_mass_test[_DoubleEG_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DoubleEG_Full")));
            h_mass_test[_DoubleEG_Full]->SetDirectory(0);
            h_pT[_DoubleEG_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_DoubleEG_Full")));
            h_pT[_DoubleEG_Full]->SetDirectory(0);
            h_rapi[_DoubleEG_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_DoubleEG_Full")));
            h_rapi[_DoubleEG_Full]->SetDirectory(0);
            h_PFiso_barrel[_DoubleEG_Full] = ((TH1D*)(h_PFiso_barrel[pr]->Clone("h_PFiso_barrel_DoubleEG_Full")));
            h_TrkIso_barrel[_DoubleEG_Full] = ((TH1D*)(h_TrkIso_barrel[pr]->Clone("h_TrkIso_barrel_DoubleEG_Full")));
            h_ECALiso_barrel[_DoubleEG_Full] = ((TH1D*)(h_ECALiso_barrel[pr]->Clone("h_ECALiso_barrel_DoubleEG_Full")));
            h_HCALiso_barrel[_DoubleEG_Full] = ((TH1D*)(h_HCALiso_barrel[pr]->Clone("h_HCALiso_barrel_DoubleEG_Full")));
            h_PFiso_endcap[_DoubleEG_Full] = ((TH1D*)(h_PFiso_endcap[pr]->Clone("h_PFiso_endcap_DoubleEG_Full")));
            h_TrkIso_endcap[_DoubleEG_Full] = ((TH1D*)(h_TrkIso_endcap[pr]->Clone("h_TrkIso_endcap_DoubleEG_Full")));
            h_ECALiso_endcap[_DoubleEG_Full] = ((TH1D*)(h_ECALiso_endcap[pr]->Clone("h_ECALiso_endcap_DoubleEG_Full")));
            h_HCALiso_endcap[_DoubleEG_Full] = ((TH1D*)(h_HCALiso_endcap[pr]->Clone("h_HCALiso_endcap_DoubleEG_Full")));
            h_PFiso_barrel[_DoubleEG_Full]->SetDirectory(0);
            h_TrkIso_barrel[_DoubleEG_Full]->SetDirectory(0);
            h_ECALiso_barrel[_DoubleEG_Full]->SetDirectory(0);
            h_HCALiso_barrel[_DoubleEG_Full]->SetDirectory(0);
            h_PFiso_endcap[_DoubleEG_Full]->SetDirectory(0);
            h_TrkIso_endcap[_DoubleEG_Full]->SetDirectory(0);
            h_ECALiso_endcap[_DoubleEG_Full]->SetDirectory(0);
            h_HCALiso_endcap[_DoubleEG_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_SingleMuon_Normal)
        {
            h_mass[_DoubleEG_Full]->Add(h_mass[pr]);
            h_mass_SS[_DoubleEG_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DoubleEG_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DoubleEG_Full]->Add(h_mass_test[pr]);
            h_pT[_DoubleEG_Full]->Add(h_pT[pr]);
            h_rapi[_DoubleEG_Full]->Add(h_rapi[pr]);
            h_PFiso_barrel[_DoubleEG_Full]->Add(h_PFiso_barrel[pr]);
            h_TrkIso_barrel[_DoubleEG_Full]->Add(h_TrkIso_barrel[pr]);
            h_ECALiso_barrel[_DoubleEG_Full]->Add(h_ECALiso_barrel[pr]);
            h_HCALiso_barrel[_DoubleEG_Full]->Add(h_HCALiso_barrel[pr]);
            h_PFiso_endcap[_DoubleEG_Full]->Add(h_PFiso_endcap[pr]);
            h_TrkIso_endcap[_DoubleEG_Full]->Add(h_TrkIso_endcap[pr]);
            h_ECALiso_endcap[_DoubleEG_Full]->Add(h_ECALiso_endcap[pr]);
            h_HCALiso_endcap[_DoubleEG_Full]->Add(h_HCALiso_endcap[pr]);
        }

        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

    } // End of pr iteration

    s_PFiso_barrel->Add(h_PFiso_barrel[_WW]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_WW]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_WW]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_WW]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_WW]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_WW]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_WW]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_WW]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_WZ]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_WZ]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_WZ]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_WZ]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_WZ]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_WZ]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_WZ]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_WZ]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_ZZ]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_ZZ]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_ZZ]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_ZZ]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_ZZ]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_ZZ]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_ZZ]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_ZZ]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_tbarW]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_tbarW]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_tbarW]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_tbarW]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_tbarW]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_tbarW]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_tbarW]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_tbarW]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_tW]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_tW]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_tW]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_tW]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_tW]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_tW]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_tW]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_tW]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_ttbar_Full]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_ttbar_Full]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_ttbar_Full]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_ttbar_Full]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_ttbar_Full]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_ttbar_Full]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_ttbar_Full]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_ttbar_Full]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_WJets_Full]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_WJets_Full]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_WJets_Full]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_WJets_Full]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_WJets_Full]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_WJets_Full]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_WJets_Full]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_WJets_Full]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_DY_Full]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_DY_Full]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_DY_Full]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_DY_Full]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_DY_Full]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_DY_Full]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_DY_Full]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_DY_Full]);
    s_PFiso_barrel->Add(h_PFiso_barrel[_GJets_Full]);
    s_TrkIso_barrel->Add(h_TrkIso_barrel[_GJets_Full]);
    s_ECALiso_barrel->Add(h_ECALiso_barrel[_GJets_Full]);
    s_HCALiso_barrel->Add(h_HCALiso_barrel[_GJets_Full]);
    s_PFiso_endcap->Add(h_PFiso_endcap[_GJets_Full]);
    s_TrkIso_endcap->Add(h_TrkIso_endcap[_GJets_Full]);
    s_ECALiso_endcap->Add(h_ECALiso_endcap[_GJets_Full]);
    s_HCALiso_endcap->Add(h_HCALiso_endcap[_GJets_Full]);

    // QCD estimation
    h_QCD_est = ((TH1D*)(h_mass[_DoubleEG_Full]->Clone("h_QCD_est")));
    Double_t err_data=0, int_data=0, err_qcd=0, int_qcd=0, err_dy=0, int_dy=0, err_tt=0, int_tt=0, err_vvnst=0, int_vvnst=0, err_wj=0, int_wj=0;
    int_data = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_data);
    h_QCD_est->SetTitle("");
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Add(h_mass[_DY_Full], -1);
    h_QCD_est->Add(h_mass[_ttbar_Full], -1);
    h_QCD_est->Add(h_mass[_VVnST], -1);
    h_QCD_est->Add(h_mass[_WJets_Full], -1);
//    h_QCD_est->Add(h_mass[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est->GetSize()-1; i_bin++)
    {
        h_QCD_est->SetBinError(i_bin, sqrt(h_QCD_est->GetBinContent(i_bin)));
        if (h_QCD_est->GetBinContent(i_bin) == 0)
        {
            h_QCD_est->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kBlack);
    h_QCD_est->SetMarkerStyle(0);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    int_dy = h_mass[_DY_Full]->IntegralAndError(1, h_mass[_DY_Full]->GetSize()-2, err_dy);
    int_tt = h_mass[_ttbar_Full]->IntegralAndError(1, h_mass[_ttbar_Full]->GetSize()-2, err_tt);
    int_vvnst = h_mass[_VVnST]->IntegralAndError(1, h_mass[_VVnST]->GetSize()-2, err_vvnst);
    int_wj = h_mass[_WJets_Full]->IntegralAndError(1, h_mass[_WJets_Full]->GetSize()-2, err_wj);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "QCD est events: " << int_qcd << "+-" << err_qcd << "(" << int_qcd/int_data << ")" << endl;
    cout << "DY events: " << int_dy << "+-" << err_dy << "(" << int_dy/int_data << ")" << endl;
    cout << "ttbar events: " << int_tt << "+-" << err_tt << "(" << int_tt/int_data << ")" << endl;
    cout << "VVnST events: " << int_vvnst << "+-" << err_vvnst << "(" << int_vvnst/int_data << ")" << endl;
    cout << "WJets events: " << int_wj << "+-" << err_wj << "(" << int_wj/int_data << ")" << endl;

    // Same-sign
    h_QCD_est_SS = ((TH1D*)(h_mass_SS[_DoubleEG_Full]->Clone("h_QCD_est_SS")));
    Double_t err_data_ss=0, int_data_ss=0, err_qcd_ss=0, int_qcd_ss=0, err_dy_ss=0, int_dy_ss=0, err_tt_ss=0, int_tt_ss=0,
            err_vvnst_ss=0, int_vvnst_ss=0, err_wj_ss=0, int_wj_ss=0;
    int_data_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_data_ss);
    h_QCD_est_SS->SetTitle("");
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_QCD_est_SS->Add(h_mass_SS[_WJets_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est_SS);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est_SS->GetSize()-1; i_bin++)
    {
        h_QCD_est_SS->SetBinError(i_bin, sqrt(h_QCD_est_SS->GetBinContent(i_bin)));
        if (h_QCD_est_SS->GetBinContent(i_bin) == 0)
        {
            h_QCD_est_SS->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est_SS->SetFillColor(kRed + 3);
    h_QCD_est_SS->SetLineColor(kBlack);
    h_QCD_est_SS->SetMarkerStyle(0);
    int_qcd_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_qcd_ss);
    int_dy_ss = h_mass_SS[_DY_Full]->IntegralAndError(1, h_mass[_DY_Full]->GetSize()-2, err_dy_ss);
    int_tt_ss = h_mass_SS[_ttbar_Full]->IntegralAndError(1, h_mass[_ttbar_Full]->GetSize()-2, err_tt_ss);
    int_vvnst_ss = h_mass_SS[_VVnST]->IntegralAndError(1, h_mass[_VVnST]->GetSize()-2, err_vvnst_ss);
    int_wj_ss = h_mass_SS[_WJets_Full]->IntegralAndError(1, h_mass[_WJets_Full]->GetSize()-2, err_wj_ss);
    cout << "QCD est events (same-sign): " << int_qcd_ss << "+-" << err_qcd_ss << endl;
    cout << "DY events (same-sign): " << int_dy_ss << "+-" << err_dy_ss << "(" << int_dy_ss/int_data_ss << ")" << endl;
    cout << "ttbar events (same-sign): " << int_tt_ss << "+-" << err_tt_ss << "(" << int_tt_ss/int_data_ss << ")" << endl;
    cout << "VVnST events (same-sign): " << int_vvnst_ss << "+-" << err_vvnst_ss << "(" << int_vvnst_ss/int_data_ss << ")" << endl;
    cout << "WJets events (same-sign): " << int_wj_ss << "+-" << err_wj_ss << "(" << int_wj_ss/int_data_ss << ")" << endl;

    // Template
    h_QCD_est_temp = ((TH1D*)(h_mass_temp[_DoubleEG_Full]->Clone("h_QCD_est_template")));
    h_QCD_est_temp->SetTitle("");
    h_QCD_est_temp->SetDirectory(0);
    h_QCD_est_temp->Add(h_mass_temp[_DY_Full], -1);
    h_QCD_est_temp->Add(h_mass_temp[_ttbar_Full], -1);
    h_QCD_est_temp->Add(h_mass_temp[_VVnST], -1);
    h_QCD_est_temp->Add(h_mass_temp[_WJets_Full], -1);
    h_QCD_est_temp->Add(h_mass_temp[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est_temp);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est_temp->GetSize()-1; i_bin++)
    {
        h_QCD_est_temp->SetBinError(i_bin, sqrt(h_QCD_est_temp->GetBinContent(i_bin)));
        if (h_QCD_est_temp->GetBinContent(i_bin) == 0)
        {
            h_QCD_est_temp->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est_temp->SetFillColor(kRed + 3);
    h_QCD_est_temp->SetLineColor(kBlack);
    h_QCD_est_temp->SetMarkerStyle(0);

    // QCD estimation from same-sign template fit
    Double_t int_QCD_fit=0, err_QCD_fit=0;
    TH1D *h_QCD_est_fit = ((TH1D*)(h_QCD_est_SS->Clone("h_QCD_est_fit")));
    h_QCD_est_fit->SetDirectory(0);
    TH1D *h_DY_fit = ((TH1D*)(h_mass[_DY_Full]->Clone("h_DY_fit")));
    TH1D *h_ttbar_fit = ((TH1D*)(h_mass[_ttbar_Full]->Clone("h_ttbar_fit")));
    TH1D *h_VVnST_fit = ((TH1D*)(h_mass[_VVnST]->Clone("h_VVnST_fit")));
    TH1D *h_WJets_fit = ((TH1D*)(h_mass[_WJets_Full]->Clone("h_WJets_fit")));    

    // Old FR (numerator = full mediumID)
//    cout << "QCD fit factor: " << 7.8129e+03/*8.4683e+03*/ / h_QCD_est_fit->Integral() << endl;
//    cout << "W+Jets fit factor: " << 1.2812e+02 / h_WJets_fit->Integral() << endl;
//    cout << "DY fit factor: " << 1.5534e+03 / h_DY_fit->Integral() << endl;
//    cout << "ttbar fit factor: " << 4.8172e+01 / h_ttbar_fit->Integral() << endl;
//    cout << "VVnST fit factor: " << 8.3653e+00 / h_VVnST_fit->Integral() << endl;
//    h_QCD_est_fit->Scale(7.8129e+03/*8.4683e+03*/ / h_QCD_est_fit->Integral());
//    h_WJets_fit->Scale(1.2812e+02 / h_WJets_fit->Integral());
//    h_DY_fit->Scale(1.5534e+03 / h_DY_fit->Integral());
//    h_ttbar_fit->Scale(4.8172e+01 / h_ttbar_fit->Integral());
//    h_VVnST_fit->Scale(8.3653e+00 / h_VVnST_fit->Integral());
    // New FR (denominator = almost mediumID (minus relPFiso), numerator = full mediumID
    cout << "QCD fit factor: " << 9.9033e+03 / h_QCD_est_fit->Integral() << endl;
    cout << "W+Jets fit factor: " << 2.4294e+02 / h_WJets_fit->Integral() << endl;
    cout << "DY fit factor: " << 7.9040e+03 / h_DY_fit->Integral() << endl;
    cout << "ttbar fit factor: " << 2.4973e+02 / h_ttbar_fit->Integral() << endl;
    cout << "VVnST fit factor: " << 3.6332e+01 / h_VVnST_fit->Integral() << endl;
    h_QCD_est_fit->Scale(9.9033e+03 / h_QCD_est_fit->Integral());
    h_WJets_fit->Scale(2.4294e+02 / h_WJets_fit->Integral());
    h_DY_fit->Scale(7.9040e+03 / h_DY_fit->Integral());
    h_ttbar_fit->Scale(2.4973e+02 / h_ttbar_fit->Integral());
    h_VVnST_fit->Scale(3.6332e+01 / h_VVnST_fit->Integral());

    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est_fit->GetSize()-1; i_bin++)
    {
        h_QCD_est_fit->SetBinError(i_bin, sqrt(h_QCD_est_fit->GetBinContent(i_bin)));
        if (h_QCD_est_fit->GetBinContent(i_bin) == 0)
        {
            h_QCD_est_fit->SetBinError(i_bin, 1);
        }
    }
    THStack *s_mass_fit = new THStack("s_mass_fit", "");
    s_mass_fit->Add(h_VVnST_fit);
    s_mass_fit->Add(h_WJets_fit);
    s_mass_fit->Add(h_ttbar_fit);
    s_mass_fit->Add(h_DY_fit);
    s_mass_fit->Add(h_QCD_est_fit);
    int_QCD_fit = h_QCD_est_fit->IntegralAndError(1, h_QCD_est_fit->GetSize()-2, err_QCD_fit);
    cout << "QCD est events (template fit): " << int_QCD_fit << "+-" << err_QCD_fit << endl;


    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass_wQCD = new myRatioPlot_t("c_mass_wQCD", s_mass_wQCD, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wQCD_SS = new myRatioPlot_t("c_mass_wQCD_SS", s_mass_wQCD_SS, h_mass_SS[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woQCD_SS = new myRatioPlot_t("c_mass_woQCD_SS", s_mass_woQCD_SS, h_mass_SS[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wQCD_temp = new myRatioPlot_t("c_mass_wQCD_template", s_mass_wQCD_temp, h_mass_temp[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woQCD_temp = new myRatioPlot_t("c_mass_woQCD_template", s_mass_woQCD_temp, h_mass_temp[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wQCD_test = new myRatioPlot_t("c_mass_wQCD_test", s_mass_wQCD_test, h_mass_test[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woQCD_test = new myRatioPlot_t("c_mass_woQCD_test", s_mass_woQCD_test, h_mass_test[_DoubleEG_Full]);
    myRatioPlot_t *RP_pT = new myRatioPlot_t("c_pT", s_pT, h_pT[_DoubleEG_Full]);
    myRatioPlot_t *RP_rapi = new myRatioPlot_t("c_rapi", s_rapi, h_rapi[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_fit = new myRatioPlot_t("c_mass_fit", s_mass_fit, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_PFiso_barrel = new myRatioPlot_t("c_PFiso_barrel", s_PFiso_barrel, h_PFiso_barrel[_DoubleEG_Full]);
    myRatioPlot_t *RP_TrkIso_barrel = new myRatioPlot_t("c_TrkIso_barrel", s_TrkIso_barrel, h_TrkIso_barrel[_DoubleEG_Full]);
    myRatioPlot_t *RP_ECALiso_barrel = new myRatioPlot_t("c_ECALiso_barrel", s_ECALiso_barrel, h_ECALiso_barrel[_DoubleEG_Full]);
    myRatioPlot_t *RP_HCALiso_barrel = new myRatioPlot_t("c_HCALiso_barrel", s_HCALiso_barrel, h_HCALiso_barrel[_DoubleEG_Full]);
    myRatioPlot_t *RP_PFiso_endcap = new myRatioPlot_t("c_PFiso_endcap", s_PFiso_endcap, h_PFiso_endcap[_DoubleEG_Full]);
    myRatioPlot_t *RP_TrkIso_endcap = new myRatioPlot_t("c_TrkIso_endcap", s_TrkIso_endcap, h_TrkIso_endcap[_DoubleEG_Full]);
    myRatioPlot_t *RP_ECALiso_endcap = new myRatioPlot_t("c_ECALiso_endcap", s_ECALiso_endcap, h_ECALiso_endcap[_DoubleEG_Full]);
    myRatioPlot_t *RP_HCALiso_endcap = new myRatioPlot_t("c_HCALiso_endcap", s_HCALiso_endcap, h_HCALiso_endcap[_DoubleEG_Full]);

    RP_mass_wQCD->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_wQCD_SS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_SS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_wQCD_temp->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_temp->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_wQCD_test->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_test->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_pT->SetPlots("p_{#lower[-0.2]{#font[12]{#scale[1.2]{ee} T}}} [GeV/c]", 0, 1000);
    RP_rapi->SetPlots("y_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}}", -4, 4);
    RP_mass_fit->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_PFiso_barrel->SetPlots(/*"I_{PF}^{rel.}"*/"relPFiso_Rho", 0, 5);
    RP_TrkIso_barrel->SetPlots("I_{Trk}^{rel.}", 0, 2);
    RP_ECALiso_barrel->SetPlots("I_{ECAL}^{rel.}", 0, 2);
    RP_HCALiso_barrel->SetPlots("I_{HCAL}^{rel.}", 0, 2);
    RP_PFiso_endcap->SetPlots(/*"I_{PF}^{rel.}"*/"relPFiso_Rho", 0, 5);
    RP_TrkIso_endcap->SetPlots("I_{Trk}^{rel.}", 0, 2);
    RP_ECALiso_endcap->SetPlots("I_{ECAL}^{rel.}", 0, 2);
    RP_HCALiso_endcap->SetPlots("I_{HCAL}^{rel.}", 0, 2);

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
    legend_wQCD->AddEntry(h_mass[_GJets_Full], "#gamma+Jets", "f");
    TLegend * legend_woQCD = ((TLegend*)(legend_wQCD->Clone()));
    legend_wQCD->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    TLegend *legend_fit = new TLegend(0.7, 0.64, 0.95, 0.95);
    legend_fit->AddEntry(h_mass[_DoubleEG_B], "Data", "pl");
    legend_fit->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    legend_fit->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_fit->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_fit->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend_fit->AddEntry(h_mass[_tbarW], "#font[12]{#scale[1.1]{tW+#bar{t}W+VV}}", "f");
    legend_woQCD->SetNColumns(2);

    RP_mass_wQCD->ImportLegend(legend_wQCD);
    RP_mass_woQCD->ImportLegend(legend_woQCD);
    RP_mass_wQCD_SS->ImportLegend(legend_wQCD);
    RP_mass_woQCD_SS->ImportLegend(legend_woQCD);
    RP_mass_wQCD_temp->ImportLegend(legend_wQCD);
    RP_mass_woQCD_temp->ImportLegend(legend_woQCD);
    RP_mass_wQCD_test->ImportLegend(legend_wQCD);
    RP_mass_woQCD_test->ImportLegend(legend_woQCD);
    RP_pT->ImportLegend(legend_wQCD);
    RP_rapi->ImportLegend(legend_wQCD);
    RP_mass_fit->ImportLegend(legend_fit);
    RP_PFiso_barrel->ImportLegend(legend_woQCD);
    RP_TrkIso_barrel->ImportLegend(legend_woQCD);
    RP_ECALiso_barrel->ImportLegend(legend_woQCD);
    RP_HCALiso_barrel->ImportLegend(legend_woQCD);
    RP_PFiso_endcap->ImportLegend(legend_woQCD);
    RP_TrkIso_endcap->ImportLegend(legend_woQCD);
    RP_ECALiso_endcap->ImportLegend(legend_woQCD);
    RP_HCALiso_endcap->ImportLegend(legend_woQCD);

    RP_mass_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD->Draw(1e-3, 1e3, 1);
    RP_mass_wQCD_SS->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD_SS->Draw(1e-3, 1e3, 1);
    RP_mass_wQCD_temp->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD_temp->Draw(1e-3, 1e3, 1);
    RP_mass_wQCD_test->Draw(1e-2, 1e6, 1);
    RP_mass_woQCD_test->Draw(1e-2, 1e6, 1);
    RP_pT->Draw(1e-2, 1e6, 0);
    RP_rapi->Draw(1e-2, 1e6, 0);
    RP_mass_fit->Draw(1e-2, 1e4, 1);
    RP_PFiso_barrel->Draw(1e-1, 1e8, 0);
    RP_TrkIso_barrel->Draw(1e-1, 1e8, 0);
    RP_ECALiso_barrel->Draw(1e-1, 1e8, 0);
    RP_HCALiso_barrel->Draw(1e-1, 1e8, 0);
    RP_PFiso_endcap->Draw(1e-1, 1e8, 0);
    RP_TrkIso_endcap->Draw(1e-1, 1e8, 0);
    RP_ECALiso_endcap->Draw(1e-1, 1e8, 0);
    RP_HCALiso_endcap->Draw(1e-1, 1e8, 0);

    // Same-sign vs opposite-sign plots
    TH1D *h_mass_QCD_MC = ((TH1D*)(h_mass[_QCDEMEnriched_Full]->Clone("h_mass_QCD_MC")));
    TH1D *h_mass_SS_QCD_MC = ((TH1D*)(h_mass_SS[_QCDEMEnriched_Full]->Clone("h_mass_SS_QCD_MC")));
    TH1D *h_mass_QCD = ((TH1D*)(h_QCD_est->Clone("h_mass_QCD")));
    TH1D *h_mass_SS_QCD = ((TH1D*)(h_QCD_est_fit->Clone("h_mass_SS_QCD")));
    h_mass_QCD_MC->SetMarkerStyle(kFullDotLarge);
    h_mass_QCD->SetMarkerStyle(kFullDotLarge);
    h_mass_QCD_MC->SetMarkerColor(kBlack);
    h_mass_QCD->SetMarkerColor(kBlack);
    h_mass_QCD_MC->SetLineColor(kBlack);
    h_mass_QCD->SetLineColor(kBlack);
    h_mass_SS_QCD_MC->SetMarkerStyle(kFullSquare);
    h_mass_SS_QCD->SetMarkerStyle(kFullSquare);
    h_mass_SS_QCD_MC->SetMarkerColor(kRed);
    h_mass_SS_QCD->SetMarkerColor(kRed);
    h_mass_SS_QCD_MC->SetLineColor(kRed);
    h_mass_SS_QCD->SetLineColor(kRed);

    myRatioPlot_t *RP_mass_QCD_SSvsOS = new myRatioPlot_t("c_mass_QCD_SSvsOS", h_mass_QCD_MC, h_mass_SS_QCD_MC);
    myRatioPlot_t *RP_mass_QCDest_SSvsOS = new myRatioPlot_t("c_mass_QCDest_SSvsOS", h_mass_QCD, h_mass_SS_QCD);
    RP_mass_QCD_SSvsOS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "SS/(OS+SS)");
    RP_mass_QCDest_SSvsOS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "Temp./Orig.");

    TLegend *l_SSvsOS_MC = new TLegend(0.5, 0.8, 0.95, 0.95);
    TLegend *l_SSvsOS = new TLegend(0.5, 0.8, 0.95, 0.95);
    l_SSvsOS_MC->AddEntry(h_mass_QCD_MC, "Any-sign QCD MC", "lp");
    l_SSvsOS_MC->AddEntry(h_mass_SS_QCD_MC, "Same-sign QCD MC", "lp");
    l_SSvsOS->AddEntry(h_mass_QCD, "Any-sign QCD estimation", "lp");
    l_SSvsOS->AddEntry(h_mass_SS_QCD, "Same-sign QCD estimation (scaled)", "lp");
    RP_mass_QCD_SSvsOS->ImportLegend(l_SSvsOS_MC);
    RP_mass_QCDest_SSvsOS->ImportLegend(l_SSvsOS);

    RP_mass_QCD_SSvsOS->Draw(5, 1e3, 1, "E");
    RP_mass_QCDest_SSvsOS->Draw(5, 1e3, 1, "E");

    // Drawing estimated QCD
    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (FR)", "f");
    TCanvas * c_QCD_est = new TCanvas("c_QCD_est", "QCD est", 750, 850);
    c_QCD_est->SetTopMargin(0.05);
    c_QCD_est->SetRightMargin(0.05);
    c_QCD_est->SetBottomMargin(0.15);
    c_QCD_est->SetLeftMargin(0.15);
    h_QCD_est->Draw("BAR");
    h_QCD_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]");
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
    l_QCD_est->Draw();
    c_QCD_est->SetLogx();
    c_QCD_est->SetGridx();
    c_QCD_est->SetGridy();
    c_QCD_est->Update();

    // Drawing estimated QCD same-sign template
    TCanvas * c_QCD_est_SS = new TCanvas("c_QCD_est_SS", "QCD est (same-sign template)", 750, 850);
    c_QCD_est_SS->SetTopMargin(0.05);
    c_QCD_est_SS->SetRightMargin(0.05);
    c_QCD_est_SS->SetBottomMargin(0.15);
    c_QCD_est_SS->SetLeftMargin(0.15);
    h_QCD_est_SS->Draw("BAR");
    h_QCD_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]");
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
    l_QCD_est->Draw();
    c_QCD_est_SS->SetLogx();
    c_QCD_est_SS->SetGridx();
    c_QCD_est_SS->SetGridy();
    c_QCD_est_SS->Update();

    // Drawing estimated QCD same sign (with tight charge requirements) template
    TCanvas * c_QCD_est_temp = new TCanvas("c_QCD_est_temp", "QCD est (tight same-sign template)", 750, 850);
    c_QCD_est_temp->SetTopMargin(0.05);
    c_QCD_est_temp->SetRightMargin(0.05);
    c_QCD_est_temp->SetBottomMargin(0.15);
    c_QCD_est_temp->SetLeftMargin(0.15);
    h_QCD_est_temp->Draw("BAR");
    h_QCD_est_temp->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]");
    h_QCD_est_temp->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est_temp->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est_temp->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est_temp->GetXaxis()->SetMoreLogLabels();
    h_QCD_est_temp->GetXaxis()->SetNoExponent();
    h_QCD_est_temp->GetYaxis()->SetTitle("Number of events");
    h_QCD_est_temp->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est_temp->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est_temp->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est_temp->GetYaxis()->SetMoreLogLabels();
    h_QCD_est_temp->GetYaxis()->SetNoExponent();
    l_QCD_est->Draw();
    c_QCD_est_temp->SetLogx();
    c_QCD_est_temp->SetGridx();
    c_QCD_est_temp->SetGridy();
    c_QCD_est_temp->Update();

    // Drawing estimated QCD from template fit
    TCanvas * c_QCD_est_fit = new TCanvas("c_QCD_est_fit", "QCD est (template fit)", 750, 850);
    c_QCD_est_fit->SetTopMargin(0.05);
    c_QCD_est_fit->SetRightMargin(0.05);
    c_QCD_est_fit->SetBottomMargin(0.15);
    c_QCD_est_fit->SetLeftMargin(0.15);
    h_QCD_est_fit->Draw("BAR");
    h_QCD_est_fit->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]");
    h_QCD_est_fit->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est_fit->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est_fit->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est_fit->GetXaxis()->SetMoreLogLabels();
    h_QCD_est_fit->GetXaxis()->SetNoExponent();
    h_QCD_est_fit->GetYaxis()->SetTitle("Number of events");
    h_QCD_est_fit->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est_fit->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est_fit->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est_fit->GetYaxis()->SetMoreLogLabels();
    h_QCD_est_fit->GetYaxis()->SetNoExponent();

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_QCD_est->Write();
    h_QCD_est_SS->Write();
    h_QCD_est_temp->Write();
    h_QCD_est_fit->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_DOWN.root", "READ");
        TFile *f_alt = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE_PFiso.root", "READ");
        TH1D *h_up, *h_down, *h_up_temp, *h_down_temp, *h_alt, *h_alt_temp, *h_tempVSsub, *h_fitUncer,
             *h_fullsysterr, *h_fullsysterr_temp, *h_fullerr, *h_fullerr_temp;
        f_up->GetObject("h_QCD_est", h_up);
        f_up->GetObject("h_QCD_est_fit", h_up_temp);
        f_down->GetObject("h_QCD_est", h_down);
        f_down->GetObject("h_QCD_est_fit", h_down_temp);
        f_alt->GetObject("h_QCD_est", h_alt);
        f_alt->GetObject("h_QCD_est_fit", h_alt_temp);
        h_tempVSsub = ((TH1D*)(h_QCD_est->Clone("h_QCD_tempVSsub")));
        h_fitUncer = ((TH1D*)(h_QCD_est_fit->Clone("h_fitUncer")));
        h_fullsysterr = ((TH1D*)(h_up->Clone("h_QCD_fullsysterr")));
        h_fullerr = ((TH1D*)(h_up->Clone("h_QCD_fullerr")));
        h_fullsysterr_temp = ((TH1D*)(h_up_temp->Clone("h_QCD_fullsysterr_temp")));
        h_fullerr_temp = ((TH1D*)(h_up_temp->Clone("h_QCD_fullerr_temp")));

        for (Int_t i=1; i<h_up->GetSize()-1; i++)
        {
            // Systematic errors
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_up_temp->SetBinContent(i, fabs(h_up_temp->GetBinContent(i)-h_down_temp->GetBinContent(i))/2);
            h_alt->SetBinContent(i, fabs(h_alt->GetBinContent(i)-h_QCD_est->GetBinContent(i)));
            h_alt_temp->SetBinContent(i, fabs(h_alt_temp->GetBinContent(i)-h_QCD_est_fit->GetBinContent(i)));
            h_tempVSsub->SetBinContent(i, fabs(h_tempVSsub->GetBinContent(i)-h_QCD_est_fit->GetBinContent(i)));
            h_fitUncer->SetBinContent(i, 0.00121*h_fitUncer->GetBinContent(i));
            Double_t err = h_up->GetBinContent(i) * h_up->GetBinContent(i) + /*h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +*/
                           h_alt->GetBinContent(i) * h_alt->GetBinContent(i);
            Double_t err_temp = h_up_temp->GetBinContent(i) * h_up_temp->GetBinContent(i) + h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +
                                h_alt_temp->GetBinContent(i) * h_alt_temp->GetBinContent(i) + h_fitUncer->GetBinContent(i) * h_fitUncer->GetBinContent(i);
            if (sqrt(err) < h_QCD_est->GetBinContent(i))
                h_fullsysterr->SetBinContent(i, sqrt(err));
            else
                h_fullsysterr->SetBinContent(i, h_QCD_est->GetBinContent(i));

            if (sqrt(err_temp) < h_QCD_est_fit->GetBinContent(i))
                h_fullsysterr_temp->SetBinContent(i, sqrt(err_temp));
            else
                h_fullsysterr_temp->SetBinContent(i, h_QCD_est_fit->GetBinContent(i));
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_QCD_est->GetBinError(i)*h_QCD_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
            h_fullerr_temp->SetBinContent(i, sqrt(h_QCD_est_fit->GetBinError(i)*h_QCD_est_fit->GetBinError(i)+
                                                  h_fullsysterr_temp->GetBinContent(i)*h_fullsysterr_temp->GetBinContent(i)));
        }
        cout << "Estimated QCD events (from template fit): " << int_QCD_fit << "+-" << err_QCD_fit << "+-" <<
                h_fullsysterr_temp->Integral() << "   (+-" << h_fullerr_temp << ")" << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_QCD_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_QCD_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetLineColor(39);
        h_draw_2->SetLineColor(39);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);

        TH1D *h_draw_1_temp = ((TH1D*)(h_QCD_est_fit->Clone("h_draw_1_temp")));
        TH1D *h_draw_2_temp = ((TH1D*)(h_QCD_est_fit->Clone("h_draw_2_temp")));
        h_draw_1_temp->Add(h_fullerr_temp, 1);
        h_draw_2_temp->Add(h_fullerr_temp, -1);
        h_draw_1_temp->SetMarkerStyle(0);
        h_draw_1_temp->SetFillColor(38);
        h_draw_1_temp->SetFillStyle(3244);
        h_draw_1_temp->SetLineColor(39);
        h_draw_2_temp->SetLineColor(39);
        h_draw_1_temp->SetDirectory(0);
        h_draw_2_temp->SetDirectory(0);

        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);
        h_fullsysterr_temp->SetDirectory(0);
        h_fullerr_temp->SetDirectory(0);
        h_tempVSsub->SetDirectory(0);

        f_up->Close();
        f_down->Close();
        f_alt->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullsysterr_temp->Write();
        h_fullerr->Write();
        h_fullerr_temp->Write();

        c_QCD_est->cd();
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
        c_QCD_est->Update();

        c_QCD_est_fit->cd();
        h_draw_1_temp->Draw("samehist");
        h_draw_2_temp->Draw("samehist");
        c_QCD_est_fit->Update();
    }// End of if (systErr)

    l_QCD_est->Draw();
    c_QCD_est_fit->SetLogx();
    c_QCD_est_fit->SetGridx();
    c_QCD_est_fit->SetGridy();
    c_QCD_est_fit->Update();

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

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special], *h_mass_temp[_EndOf_Data_Special], *h_mass_test[_EndOf_Data_Special],
         *h_MET[_EndOf_Data_Special], *h_pT[_EndOf_Data_Special], *h_rapi[_EndOf_Data_Special];
    TH1D *h_WJET_est, *h_WJET_est_SS, *h_WJET_est_temp, *h_WJET_est_fit, *h_WJET_est_fit_2;
    THStack *s_mass_wWJET = new THStack("s_mass_wWJET", "");
    THStack *s_mass_woWJET = new THStack("s_mass_woWJET", "");
    THStack *s_mass_wWJET_SS = new THStack("s_mass_wWJET_SS", "");
    THStack *s_mass_woWJET_SS = new THStack("s_mass_woWJET_SS", "");
    THStack *s_mass_wWJET_temp = new THStack("s_mass_wWJET_template", "");
    THStack *s_mass_woWJET_temp = new THStack("s_mass_woWJET_template", "");
    THStack *s_mass_test = new THStack("s_mass_wWJET_test", "");
    THStack *s_MET = new THStack("s_MET", "");
    THStack *s_pT = new THStack("s_pT", "");
    THStack *s_rapi = new THStack("s_rapi", "");
    Color_t color = kBlack;

    // Getting data-driven QCD to subtract from data
    TFile *f_QCD = new TFile("/media/sf_DATA/SelectedEE/Histos/EstQCD_EE.root");
    TH1D *h_QCD_est, *h_QCD_est_SS, *h_QCD_est_temp, *h_QCD_est_fit;

    f_QCD->GetObject("h_QCD_est", h_QCD_est);
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Scale(2);
    s_mass_wWJET->Add(h_QCD_est);
    s_mass_woWJET->Add(h_QCD_est);

    f_QCD->GetObject("h_QCD_est_SS", h_QCD_est_SS);
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->Scale(2);
    s_mass_wWJET_SS->Add(h_QCD_est_SS);
    s_mass_woWJET_SS->Add(h_QCD_est_SS);

    f_QCD->GetObject("h_QCD_est_template", h_QCD_est_temp);
    h_QCD_est_temp->SetDirectory(0);
    h_QCD_est_temp->Scale(2);
    s_mass_wWJET_temp->Add(h_QCD_est_temp);
    s_mass_woWJET_temp->Add(h_QCD_est_temp);

    f_QCD->GetObject("h_QCD_est_fit", h_QCD_est_fit);
    h_QCD_est_fit->SetDirectory(0);
    h_QCD_est_fit->Scale(2);

    for (Process_t pr=_DoubleEG_H; pr>=_DY_10to50; pr=previous(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        f->GetObject("h_mass_template_"+Mgr.Procname[pr], h_mass_temp[pr]);
        f->GetObject("h_mass_test_"+Mgr.Procname[pr], h_mass_test[pr]);
        f->GetObject("h_MET_"+Mgr.Procname[pr], h_MET[pr]);
        f->GetObject("h_pT_"+Mgr.Procname[pr], h_pT[pr]);
        f->GetObject("h_rapi_"+Mgr.Procname[pr], h_rapi[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);
        h_mass_temp[pr]->SetDirectory(0);
        h_mass_test[pr]->SetDirectory(0);
        h_MET[pr]->SetDirectory(0);
        h_pT[pr]->SetDirectory(0);
        h_rapi[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
            removeNegativeBins(h_mass_temp[pr]);
            removeNegativeBins(h_mass_test[pr]);
            removeNegativeBins(h_MET[pr]);
            removeNegativeBins(h_pT[pr]);
            removeNegativeBins(h_rapi[pr]);
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
        else if (pr < _EndOf_GJets_Normal) color = kYellow + 3;

        if (pr < _DoubleEG_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
            h_mass_temp[pr]->SetFillColor(color);
            h_mass_temp[pr]->SetLineColor(color);
            h_mass_test[pr]->SetFillColor(color);
            h_mass_test[pr]->SetLineColor(color);
            h_MET[pr]->SetFillColor(color);
            h_MET[pr]->SetLineColor(color);
            h_pT[pr]->SetFillColor(color);
            h_pT[pr]->SetLineColor(color);
            h_rapi[pr]->SetFillColor(color);
            h_rapi[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
            h_mass_temp[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_temp[pr]->SetMarkerColor(kBlack);
            h_mass_temp[pr]->SetLineColor(kBlack);
            h_mass_test[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test[pr]->SetMarkerColor(kBlack);
            h_mass_test[pr]->SetLineColor(kBlack);
            h_MET[pr]->SetMarkerStyle(kFullDotLarge);
            h_MET[pr]->SetMarkerColor(kBlack);
            h_MET[pr]->SetLineColor(kBlack);
            h_pT[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT[pr]->SetMarkerColor(kBlack);
            h_pT[pr]->SetLineColor(kBlack);
            h_rapi[pr]->SetMarkerStyle(kFullDotLarge);
            h_rapi[pr]->SetMarkerColor(kBlack);
            h_rapi[pr]->SetLineColor(kBlack);
        }

        // Filling stack histograms
        if (pr < _EndOf_WJets_Normal /*|| (pr >= _GJets_20to100 && pr <= _GJets_2000to5000)*/)
        {
            s_mass_wWJET->Add(h_mass[pr]);
            s_mass_wWJET_SS->Add(h_mass_SS[pr]);
            s_mass_wWJET_temp->Add(h_mass_temp[pr]);
            if (pr < _WJets /*|| (pr >= _GJets_20to100 && pr <= _GJets_2000to5000)*/)
            {
                s_mass_woWJET->Add(h_mass[pr]);
                s_mass_woWJET_SS->Add(h_mass_SS[pr]);
                s_mass_woWJET_temp->Add(h_mass_temp[pr]);
            }
        }
        if (pr < _DoubleEG_B)
        {
            s_MET->Add(h_MET[pr]);
            s_mass_test->Add(h_mass_test[pr]);
            s_pT->Add(h_pT[pr]);
            s_rapi->Add(h_rapi[pr]);
        }

        // Adding up for convenience
        if (pr == _DoubleEG_H)
        {
            h_mass[_DoubleEG_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DoubleEG_Full")));
            h_mass[_DoubleEG_Full]->SetDirectory(0);
            h_mass_SS[_DoubleEG_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DoubleEG_Full")));
            h_mass_SS[_DoubleEG_Full]->SetDirectory(0);
            h_mass_temp[_DoubleEG_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DoubleEG_Full")));
            h_mass_temp[_DoubleEG_Full]->SetDirectory(0);
            h_mass_test[_DoubleEG_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DoubleEG_Full")));
            h_mass_test[_DoubleEG_Full]->SetDirectory(0);
            h_MET[_DoubleEG_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_DoubleEG_Full")));
            h_MET[_DoubleEG_Full]->SetDirectory(0);
            h_pT[_DoubleEG_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_DoubleEG_Full")));
            h_pT[_DoubleEG_Full]->SetDirectory(0);
            h_rapi[_DoubleEG_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_DoubleEG_Full")));
            h_rapi[_DoubleEG_Full]->SetDirectory(0);
        }
        else if (pr >= _DoubleEG_B)
        {
            h_mass[_DoubleEG_Full]->Add(h_mass[pr]);
            h_mass_SS[_DoubleEG_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DoubleEG_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DoubleEG_Full]->Add(h_mass_test[pr]);
            h_MET[_DoubleEG_Full]->Add(h_MET[pr]);
            h_pT[_DoubleEG_Full]->Add(h_pT[pr]);
            h_rapi[_DoubleEG_Full]->Add(h_rapi[pr]);
        }
        else if (pr == _GJets_2000to5000)
        {
            h_mass[_GJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_GJets_Full")));
            h_mass[_GJets_Full]->SetDirectory(0);
            h_mass_SS[_GJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_GJets_Full")));
            h_mass_SS[_GJets_Full]->SetDirectory(0);
            h_mass_temp[_GJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_GJets_Full")));
            h_mass_temp[_GJets_Full]->SetDirectory(0);
            h_mass_test[_GJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_GJets_Full")));
            h_mass_test[_GJets_Full]->SetDirectory(0);
            h_MET[_GJets_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_GJets_Full")));
            h_MET[_GJets_Full]->SetDirectory(0);
            h_pT[_GJets_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_GJets_Full")));
            h_pT[_GJets_Full]->SetDirectory(0);
            h_rapi[_GJets_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_GJets_Full")));
            h_rapi[_GJets_Full]->SetDirectory(0);
        }
        else if (pr >= _GJets_20to100)
        {
            h_mass[_GJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_GJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_GJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_GJets_Full]->Add(h_mass_test[pr]);
            h_MET[_GJets_Full]->Add(h_MET[pr]);
            h_pT[_GJets_Full]->Add(h_pT[pr]);
            h_rapi[_GJets_Full]->Add(h_rapi[pr]);
        }
        else if (pr == _QCDEMEnriched_300toInf)
        {
            h_mass[_QCDEMEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDEMEnriched_Full")));
            h_mass[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_SS[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_QCDEMEnriched_Full")));
            h_mass_SS[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_temp[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_QCDEMEnriched_Full")));
            h_mass_temp[_QCDEMEnriched_Full]->SetDirectory(0);
            h_mass_test[_QCDEMEnriched_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_QCDEMEnriched_Full")));
            h_mass_test[_QCDEMEnriched_Full]->SetDirectory(0);
            h_MET[_QCDEMEnriched_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_QCDEMEnriched_Full")));
            h_MET[_QCDEMEnriched_Full]->SetDirectory(0);
            h_pT[_QCDEMEnriched_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_QCDEMEnriched_Full")));
            h_pT[_QCDEMEnriched_Full]->SetDirectory(0);
            h_rapi[_QCDEMEnriched_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_QCDEMEnriched_Full")));
            h_rapi[_QCDEMEnriched_Full]->SetDirectory(0);
        }
        else if (pr >= _QCDEMEnriched_20to30)
        {
            h_mass[_QCDEMEnriched_Full]->Add(h_mass[pr]);
            h_mass_SS[_QCDEMEnriched_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_QCDEMEnriched_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_QCDEMEnriched_Full]->Add(h_mass_test[pr]);
            h_MET[_QCDEMEnriched_Full]->Add(h_MET[pr]);
            h_pT[_QCDEMEnriched_Full]->Add(h_pT[pr]);
            h_rapi[_QCDEMEnriched_Full]->Add(h_rapi[pr]);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
            h_mass_temp[_WJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_WJets_Full")));
            h_mass_temp[_WJets_Full]->SetDirectory(0);
            h_mass_test[_WJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_WJets_Full")));
            h_mass_test[_WJets_Full]->SetDirectory(0);
            h_MET[_WJets_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_WJets_Full")));
            h_MET[_WJets_Full]->SetDirectory(0);
            h_pT[_WJets_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_WJets_Full")));
            h_pT[_WJets_Full]->SetDirectory(0);
            h_rapi[_WJets_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_WJets_Full")));
            h_rapi[_WJets_Full]->SetDirectory(0);
        }
        else if (pr >= _WJets)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_WJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_WJets_Full]->Add(h_mass_test[pr]);
            h_MET[_WJets_Full]->Add(h_MET[pr]);
            h_pT[_WJets_Full]->Add(h_pT[pr]);
            h_rapi[_WJets_Full]->Add(h_rapi[pr]);
        }
        else if (pr == _WW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
            h_mass_temp[_VVnST] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_VVnST")));
            h_mass_temp[_VVnST]->SetDirectory(0);
            h_mass_test[_VVnST] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_VVnST")));
            h_mass_test[_VVnST]->SetDirectory(0);
            h_MET[_VVnST] = ((TH1D*)(h_MET[pr]->Clone("h_MET_VVnST")));
            h_MET[_VVnST]->SetDirectory(0);
            h_pT[_VVnST] = ((TH1D*)(h_pT[pr]->Clone("h_pT_VVnST")));
            h_pT[_VVnST]->SetDirectory(0);
            h_rapi[_VVnST] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_VVnST")));
            h_rapi[_VVnST]->SetDirectory(0);
        }
        else if (pr >= _tW)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
            h_mass_temp[_VVnST]->Add(h_mass_temp[pr]);
            h_mass_test[_VVnST]->Add(h_mass_test[pr]);
            h_MET[_VVnST]->Add(h_MET[pr]);
            h_pT[_VVnST]->Add(h_pT[pr]);
            h_rapi[_VVnST]->Add(h_rapi[pr]);
        }
        else if (pr == _ttbar_1000toInf)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
            h_mass_temp[_ttbar_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_ttbar_Full")));
            h_mass_temp[_ttbar_Full]->SetDirectory(0);
            h_mass_test[_ttbar_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_ttbar_Full")));
            h_mass_test[_ttbar_Full]->SetDirectory(0);
            h_MET[_ttbar_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_ttbar_Full")));
            h_MET[_ttbar_Full]->SetDirectory(0);
            h_pT[_ttbar_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_ttbar_Full")));
            h_pT[_ttbar_Full]->SetDirectory(0);
            h_rapi[_ttbar_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_ttbar_Full")));
            h_rapi[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr >= _ttbar)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_ttbar_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_ttbar_Full]->Add(h_mass_test[pr]);
            h_MET[_ttbar_Full]->Add(h_MET[pr]);
            h_pT[_ttbar_Full]->Add(h_pT[pr]);
            h_rapi[_ttbar_Full]->Add(h_rapi[pr]);
        }
        else if (pr == _DY_2000to3000)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
            h_mass_temp[_DY_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DY_Full")));
            h_mass_temp[_DY_Full]->SetDirectory(0);
            h_mass_test[_DY_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DY_Full")));
            h_mass_test[_DY_Full]->SetDirectory(0);
            h_MET[_DY_Full] = ((TH1D*)(h_MET[pr]->Clone("h_MET_DY_Full")));
            h_MET[_DY_Full]->SetDirectory(0);
            h_pT[_DY_Full] = ((TH1D*)(h_pT[pr]->Clone("h_pT_DY_Full")));
            h_pT[_DY_Full]->SetDirectory(0);
            h_rapi[_DY_Full] = ((TH1D*)(h_rapi[pr]->Clone("h_rapi_DY_Full")));
            h_rapi[_DY_Full]->SetDirectory(0);
        }
        else if (pr >= _DY_10to50)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DY_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DY_Full]->Add(h_mass_test[pr]);
            h_MET[_DY_Full]->Add(h_MET[pr]);
            h_pT[_DY_Full]->Add(h_pT[pr]);
            h_rapi[_DY_Full]->Add(h_rapi[pr]);
        }

        if (pr == _QCDEMEnriched_20to30) pr = _EndOf_WJets_Normal; // next -- WJets_ext2v5
        if (pr == _ttbar) pr = _EndOf_DY_Normal; // next -- DY_2000to3000

    } // End of pr iteration

    // W+Jets estimation by subtraction
    Double_t int_data=0, err_data=0, int_wjet=0, err_wjet=0, int_DY=0, err_DY=0, int_tt=0, err_tt=0, int_vvnst=0, err_vvnst=0, int_qcd=0, err_qcd=0;
    h_WJET_est = ((TH1D*)(h_mass[_DoubleEG_Full]->Clone("h_WJET_est")));
    int_data = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_data);
    h_WJET_est->SetTitle("");
    h_WJET_est->SetDirectory(0);
    h_WJET_est->Add(h_mass[_DY_Full], -1);
    h_WJET_est->Add(h_mass[_ttbar_Full], -1);
    h_WJET_est->Add(h_mass[_VVnST], -1);
//    h_WJET_est->Add(h_mass[_GJets_Full], -1);
    h_WJET_est->Add(h_QCD_est, -1);
    removeNegativeBins(h_WJET_est);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est->GetSize()-1; i_bin++)
    {
        h_WJET_est->SetBinError(i_bin, sqrt(h_WJET_est->GetBinContent(i_bin)));
        if (h_WJET_est->GetBinContent(i_bin) == 0)
        {
            h_WJET_est->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est->SetFillColor(kRed - 2);
    h_WJET_est->SetLineColor(kBlack);
    h_WJET_est->SetMarkerStyle(0);
    int_wjet = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_wjet);
    int_DY = h_mass[_DY_Full]->IntegralAndError(1, h_mass[_DY_Full]->GetSize()-2, err_DY);
    int_tt = h_mass[_ttbar_Full]->IntegralAndError(1, h_mass[_ttbar_Full]->GetSize()-2, err_tt);
    int_vvnst = h_mass[_VVnST]->IntegralAndError(1, h_mass[_VVnST]->GetSize()-2, err_vvnst);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "WJets est events (subraction): " << int_wjet << "+-" << err_wjet << "(" << int_wjet/int_data << ")" << endl;
    cout << "DY events (subraction): " << int_DY << "+-" << err_DY << "(" << int_DY/int_data << ")" << endl;
    cout << "ttbar events (subraction): " << int_tt << "+-" << err_tt << "(" << int_tt/int_data << ")" << endl;
    cout << "VVnST events (subraction): " << int_vvnst << "+-" << err_vvnst << "(" << int_vvnst/int_data << ")" << endl;
    cout << "QCD events (subraction): " << int_qcd << "+-" << err_qcd << "(" << int_qcd/int_data << ")" << endl;

    // Same-sign template
    Double_t int_data_ss=0, err_data_ss=0, int_wjet_ss=0, err_wjet_ss=0, int_DY_ss=0, err_DY_ss=0, int_tt_ss=0, err_tt_ss=0,
            int_vvnst_ss=0, err_vvnst_ss=0, int_qcd_ss=0, err_qcd_ss=0;
    h_WJET_est_SS = ((TH1D*)(h_mass_SS[_DoubleEG_Full]->Clone("h_WJET_est_SS")));
    int_data_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_data_ss);
    h_WJET_est_SS->SetTitle("");
    h_WJET_est_SS->SetDirectory(0);
    h_WJET_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_VVnST], -1);
//    h_WJET_est_SS->Add(h_mass_SS[_GJets_Full], -1);
    h_WJET_est_SS->Add(h_QCD_est_SS, -1);
    removeNegativeBins(h_WJET_est_SS);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_SS->GetSize()-1; i_bin++)
    {
        h_WJET_est_SS->SetBinError(i_bin, sqrt(h_WJET_est_SS->GetBinContent(i_bin)));
        if (h_WJET_est_SS->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_SS->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est_SS->SetFillColor(kRed - 2);
    h_WJET_est_SS->SetLineColor(kBlack);
    h_WJET_est_SS->SetMarkerStyle(0);
    int_wjet_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_wjet_ss);
    int_DY_ss = h_mass_SS[_DY_Full]->IntegralAndError(1, h_mass_SS[_DY_Full]->GetSize()-2, err_DY_ss);
    int_tt_ss = h_mass_SS[_ttbar_Full]->IntegralAndError(1, h_mass_SS[_ttbar_Full]->GetSize()-2, err_tt_ss);
    int_vvnst_ss = h_mass_SS[_VVnST]->IntegralAndError(1, h_mass_SS[_VVnST]->GetSize()-2, err_vvnst_ss);
    int_qcd_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_qcd_ss);
    cout << "Data events: " << int_data << "+-" << err_data_ss << endl;
    cout << "WJets est events (same-sign): " << int_wjet_ss << "+-" << err_wjet_ss << "(" << int_wjet_ss/int_data_ss << ")" << endl;
    cout << "DY events (same-sign): " << int_DY_ss << "+-" << err_DY_ss << "(" << int_DY_ss/int_data_ss << ")" << endl;
    cout << "ttbar events (same-sign): " << int_tt_ss << "+-" << err_tt_ss << "(" << int_tt_ss/int_data_ss << ")" << endl;
    cout << "VVnST events (same-sign): " << int_vvnst_ss << "+-" << err_vvnst_ss << "(" << int_vvnst_ss/int_data_ss << ")" << endl;
    cout << "QCD events (same-sign): " << int_qcd_ss << "+-" << err_qcd_ss << "(" << int_qcd_ss/int_data_ss << ")" << endl;


    // Same-sign template with tight charge requirement
    Double_t int_data_temp=0, err_data_temp=0, int_wjet_temp=0, err_wjet_temp=0, int_DY_temp=0, err_DY_temp=0, int_tt_temp=0, err_tt_temp=0,
             int_vvnst_temp=0, err_vvnst_temp=0, int_qcd_temp=0, err_qcd_temp=0;
    h_WJET_est_temp = ((TH1D*)(h_mass_temp[_DoubleEG_Full]->Clone("h_WJET_est_template")));
    int_data_temp = h_WJET_est_temp->IntegralAndError(1, h_WJET_est_temp->GetSize()-2, err_data_temp);
    h_WJET_est_temp->SetTitle("");
    h_WJET_est_temp->SetDirectory(0);
    h_WJET_est_temp->Add(h_mass_temp[_DY_Full], -1);
    h_WJET_est_temp->Add(h_mass_temp[_ttbar_Full], -1);
    h_WJET_est_temp->Add(h_mass_temp[_VVnST], -1);
//    h_WJET_est_temp->Add(h_mass_temp[_GJets_Full], -1);
    h_WJET_est_temp->Add(h_QCD_est_temp, -1);
    removeNegativeBins(h_WJET_est_temp);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_temp->GetSize()-1; i_bin++)
    {
        h_WJET_est_temp->SetBinError(i_bin, sqrt(h_WJET_est_temp->GetBinContent(i_bin)));
        if (h_WJET_est_temp->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_temp->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est_temp->SetFillColor(kRed - 2);
    h_WJET_est_temp->SetLineColor(kBlack);
    h_WJET_est_temp->SetMarkerStyle(0);
    int_wjet_temp = h_WJET_est_temp->IntegralAndError(1, h_WJET_est_temp->GetSize()-2, err_wjet_temp);
    int_DY_temp = h_mass_temp[_DY_Full]->IntegralAndError(1, h_mass_temp[_DY_Full]->GetSize()-2, err_DY_temp);
    int_tt_temp = h_mass_temp[_ttbar_Full]->IntegralAndError(1, h_mass_temp[_ttbar_Full]->GetSize()-2, err_tt_temp);
    int_vvnst_temp = h_mass_temp[_VVnST]->IntegralAndError(1, h_mass_temp[_VVnST]->GetSize()-2, err_vvnst_temp);
    int_qcd_temp = h_QCD_est_temp->IntegralAndError(1, h_QCD_est_temp->GetSize()-2, err_qcd_temp);
    cout << "Data events: " << int_data << "+-" << err_data_ss << endl;
    cout << "WJets est events (tight same-sign): " << int_wjet_temp << "+-" << err_wjet_temp << "(" << int_wjet_temp/int_data_temp << ")" << endl;
    cout << "DY events (tight same-sign): " << int_DY_temp << "+-" << err_DY_temp << "(" << int_DY_temp/int_data_temp << ")" << endl;
    cout << "ttbar events (tight same-sign): " << int_tt_temp << "+-" << err_tt_temp << "(" << int_tt_temp/int_data_temp << ")" << endl;
    cout << "VVnST events (tight same-sign): " << int_vvnst_temp << "+-" << err_vvnst_temp << "(" << int_vvnst_temp/int_data_temp << ")" << endl;
    cout << "QCD events (tight same-sign): " << int_qcd_temp << "+-" << err_qcd_temp << "(" << int_qcd_temp/int_data_temp << ")" << endl;


    // W+Jets estimation from same-sign template fit
    cout << "---- estimation from same-sign template fit ----" << endl;
    Double_t int_wjet_fit=0, err_wjet_fit=0;
//    h_WJET_est_fit = ((TH1D*)(h_WJET_est_temp->Clone("h_WJET_est_fit")));
    h_WJET_est_fit = ((TH1D*)(h_WJET_est_SS->Clone("h_WJET_est_fit")));
    TH1D *h_DY_fit = ((TH1D*)(h_mass[_DY_Full]->Clone("h_DY_fit")));
    TH1D *h_ttbar_fit = ((TH1D*)(h_mass[_ttbar_Full]->Clone("h_ttbar_fit")));
    TH1D *h_VVnST_fit = ((TH1D*)(h_mass[_VVnST]->Clone("h_VVnST_fit")));
    TH1D *h_QCD_fit = ((TH1D*)(h_QCD_est_fit->Clone("h_QCD_fit")));

    // Old FR (numerator = full mediumID)
//    cout << "W+Jets fit factor: " << 1.1090e+04/*1.2992e+04*/ / h_WJET_est_fit->Integral() << endl;
//    cout << "DY fit factor: " << 2.9014e+05 / h_DY_fit->Integral() << endl;
//    cout << "ttbar fit factor: " << 3.9280e+03 / h_ttbar_fit->Integral() << endl;
//    cout << "VVnST fit factor: " << 1.3173e+03 / h_VVnST_fit->Integral() << endl;
//    cout << "QCD fit factor: " << 1.5424e+04 / h_QCD_fit->Integral() << endl;
//    h_WJET_est_fit->Scale(1.1090e+04/*1.2992e+04*/ / h_WJET_est_fit->Integral());
//    h_DY_fit->Scale(2.9014e+05 / h_DY_fit->Integral());
//    h_ttbar_fit->Scale(3.9280e+03 / h_ttbar_fit->Integral());
//    h_VVnST_fit->Scale(1.3173e+03 / h_VVnST_fit->Integral());
//    h_QCD_fit->Scale(1.5424e+04 / h_QCD_fit->Integral());
    // New FR (denominator = almost mediumID (minus relPFiso), numerator = mediumID)
    cout << "W+Jets fit factor: " << 7.4789e+03 / h_WJET_est_fit->Integral() << endl;
    cout << "DY fit factor: " << 2.6971e+05 / h_DY_fit->Integral() << endl;
    cout << "ttbar fit factor: " << 3.7000e+03 / h_ttbar_fit->Integral() << endl;
    cout << "VVnST fit factor: " << 1.0802e+03 / h_VVnST_fit->Integral() << endl;
    cout << "QCD fit factor: " << 1.9583e+03 / h_QCD_fit->Integral() << endl;
    h_WJET_est_fit->Scale(7.4789e+03 / h_WJET_est_fit->Integral());
    h_DY_fit->Scale(2.6971e+05 / h_DY_fit->Integral());
    h_ttbar_fit->Scale(3.7000e+03 / h_ttbar_fit->Integral());
    h_VVnST_fit->Scale(1.0802e+03 / h_VVnST_fit->Integral());
    h_QCD_fit->Scale(1.9583e+03 / h_QCD_fit->Integral());

    h_VVnST_fit->SetFillColor(kGreen+2);
    h_VVnST_fit->SetLineColor(kGreen+2);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_fit->GetSize()-1; i_bin++)
    {
        h_WJET_est_fit->SetBinError(i_bin, sqrt(h_WJET_est_fit->GetBinContent(i_bin)));
        if (h_WJET_est_fit->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_fit->SetBinError(i_bin, 1);
        }
    }
    THStack *s_mass_fit = new THStack("s_mass_fit", "");
    s_mass_fit->Add(h_QCD_fit);
    s_mass_fit->Add(h_VVnST_fit);
    s_mass_fit->Add(h_ttbar_fit);
    s_mass_fit->Add(h_WJET_est_fit);
    s_mass_fit->Add(h_DY_fit);
    int_wjet_fit = h_WJET_est_fit->IntegralAndError(1, h_WJET_est_fit->GetSize()-2, err_wjet_fit);
    cout << "WJets est events (same-sign template fit): " << int_wjet_fit << "+-" << err_wjet_fit << endl;


    // W+Jets estimation from TIGHT same-sign template fit
    cout << "---- estimation from tight same-sign template fit ----" << endl;
    Double_t int_wjet_fit_2=0, err_wjet_fit_2=0;
    h_WJET_est_fit_2 = ((TH1D*)(h_WJET_est_temp->Clone("h_WJET_est_fit_2")));
    TH1D *h_DY_fit_2 = ((TH1D*)(h_mass[_DY_Full]->Clone("h_DY_fit_2")));
    TH1D *h_ttbar_fit_2 = ((TH1D*)(h_mass[_ttbar_Full]->Clone("h_ttbar_fit_2")));
    TH1D *h_VVnST_fit_2 = ((TH1D*)(h_mass[_VVnST]->Clone("h_VVnST_fit_2")));
    TH1D *h_QCD_fit_2 = ((TH1D*)(h_QCD_est_fit->Clone("h_QCD_fit_2")));

    // Old FR (numerator = full mediumID)
    // empty

    // New FR (denominator = almost mediumID (minus relPFiso), numerator = mediumID)
    cout << "W+Jets fit factor: " << 6.1734e+03 / h_WJET_est_fit_2->Integral() << endl;
    cout << "DY fit factor: " << 2.7116e+05 / h_DY_fit_2->Integral() << endl;
    cout << "ttbar fit factor: " << 3.7547e+03 / h_ttbar_fit_2->Integral() << endl;
    cout << "VVnST fit factor: " << 9.1039e+02 / h_VVnST_fit_2->Integral() << endl;
    cout << "QCD fit factor: " << 1.9875e+03 / h_QCD_fit_2->Integral() << endl;
    h_WJET_est_fit_2->Scale(6.1734e+03 / h_WJET_est_fit_2->Integral());
    h_DY_fit_2->Scale(2.7116e+05 / h_DY_fit_2->Integral());
    h_ttbar_fit_2->Scale(3.7547e+03 / h_ttbar_fit_2->Integral());
    h_VVnST_fit_2->Scale(9.1039e+02 / h_VVnST_fit_2->Integral());
    h_QCD_fit_2->Scale(1.9875e+03 / h_QCD_fit_2->Integral());

    h_VVnST_fit_2->SetFillColor(kGreen+2);
    h_VVnST_fit_2->SetLineColor(kGreen+2);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_fit_2->GetSize()-1; i_bin++)
    {
        h_WJET_est_fit_2->SetBinError(i_bin, sqrt(h_WJET_est_fit_2->GetBinContent(i_bin)));
        if (h_WJET_est_fit_2->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_fit_2->SetBinError(i_bin, 1);
        }
    }
    THStack *s_mass_fit_2 = new THStack("s_mass_fit_2", "");
    s_mass_fit_2->Add(h_QCD_fit_2);
    s_mass_fit_2->Add(h_VVnST_fit_2);
    s_mass_fit_2->Add(h_ttbar_fit_2);
    s_mass_fit_2->Add(h_WJET_est_fit_2);
    s_mass_fit_2->Add(h_DY_fit_2);
    int_wjet_fit_2 = h_WJET_est_fit_2->IntegralAndError(1, h_WJET_est_fit_2->GetSize()-2, err_wjet_fit_2);
    cout << "WJets est events (tight same-sign template fit): " << int_wjet_fit_2 << "+-" << err_wjet_fit_2 << endl;


    myRatioPlot_t *RP_MET = new myRatioPlot_t("c_MET", s_MET, h_MET[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wWJET = new myRatioPlot_t("c_mass_wWJET", s_mass_wWJET, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woWJET = new myRatioPlot_t("c_mass_woWJET", s_mass_woWJET, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wWJET_SS = new myRatioPlot_t("c_mass_wWJET_SS", s_mass_wWJET_SS, h_mass_SS[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woWJET_SS = new myRatioPlot_t("c_mass_woWJET_SS", s_mass_woWJET_SS, h_mass_SS[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_wWJET_temp = new myRatioPlot_t("c_mass_wWJET_template", s_mass_wWJET_temp, h_mass_temp[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_woWJET_temp = new myRatioPlot_t("c_mass_woWJET_template", s_mass_woWJET_temp, h_mass_temp[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_test = new myRatioPlot_t("c_mass_test", s_mass_test, h_mass_test[_DoubleEG_Full]);
    myRatioPlot_t *RP_pT = new myRatioPlot_t("c_pT", s_pT, h_pT[_DoubleEG_Full]);
    myRatioPlot_t *RP_rapi = new myRatioPlot_t("c_rapi", s_rapi, h_rapi[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_fit = new myRatioPlot_t("c_mass_fit", s_mass_fit, h_mass[_DoubleEG_Full]);
    myRatioPlot_t *RP_mass_fit_2 = new myRatioPlot_t("c_mass_fit_2", s_mass_fit_2, h_mass[_DoubleEG_Full]);

    RP_MET->SetPlots("p_{T}^{miss} [GeV/c]", 0, 1000);
    RP_mass_wWJET->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_wWJET_SS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET_SS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_wWJET_temp->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET_temp->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_test->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_pT->SetPlots("p_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}} T}} [GeV/c]", 15, 3000);
    RP_rapi->SetPlots("y_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}}", 15, 3000);
    RP_mass_fit->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_fit_2->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000);

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
//    legend_wWJET->AddEntry(h_mass[_GJets_20to100], "#gamma+Jets", "f");
    legend_wWJET->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
//    legend_woWJET->AddEntry(h_mass[_GJets_20to100], "#gamma+Jets", "f");
    legend_woWJET->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
    TLegend * legend_fit = new TLegend(0.7, 0.64, 0.95, 0.95);
    legend_fit->AddEntry(h_mass[_DoubleEG_B], "Data", "pl");
    legend_fit->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_fit->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend_fit->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_fit->AddEntry(h_mass[_tbarW], "#font[12]{#scale[1.1]{tW+#bar{t}W+VV}}", "f");
    legend_fit->AddEntry(h_mass[_QCDEMEnriched_20to30], "#font[12]{#scale[1.1]{QCD}}", "f");
//    legend_fit->AddEntry(h_mass[_GJets_20to100], "#gamma+Jets", "f");

    RP_MET->ImportLegend(legend_wWJET);
    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);
    RP_mass_wWJET_SS->ImportLegend(legend_wWJET);
    RP_mass_woWJET_SS->ImportLegend(legend_woWJET);
    RP_mass_wWJET_temp->ImportLegend(legend_wWJET);
    RP_mass_woWJET_temp->ImportLegend(legend_woWJET);
    RP_mass_test->ImportLegend(legend_wWJET);
    RP_pT->ImportLegend(legend_wWJET);
    RP_rapi->ImportLegend(legend_wWJET);
    RP_mass_fit->ImportLegend(legend_fit);
    RP_mass_fit_2->ImportLegend(legend_fit);

    RP_MET->Draw(1e-2, 1e5, 0);
    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);
    RP_mass_wWJET_SS->Draw(1e-2, 1e4, 1);
    RP_mass_woWJET_SS->Draw(1e-2, 1e4, 1);
    RP_mass_wWJET_temp->Draw(1e-2, 1e4, 1);
    RP_mass_woWJET_temp->Draw(1e-2, 1e4, 1);
    RP_mass_test->Draw(1e-2, 1e6, 1);
    RP_pT->Draw(1e-2, 1e6, 0);
    RP_rapi->Draw(1e-2, 1e6, 0);
    RP_mass_fit->Draw(10, 1e5, 1);
    RP_mass_fit_2->Draw(10, 1e5, 1);

    // Same-sign vs opposite-sign plots
    TH1D *h_mass_WJET_MC = ((TH1D*)(h_mass[_WJets_Full]->Clone("h_mass_WJET_MC")));
    TH1D *h_mass_SS_WJET_MC = ((TH1D*)(h_mass_SS[_WJets_Full]->Clone("h_mass_SS_WJET_MC")));
    TH1D *h_mass_temp_WJET_MC = ((TH1D*)(h_mass_temp[_WJets_Full]->Clone("h_mass_temp_WJET_MC")));
    TH1D *h_mass_SS_WJET = ((TH1D*)(h_WJET_est_SS->Clone("h_mass_SS_WJET")));
    TH1D *h_mass_temp_WJET = ((TH1D*)(h_WJET_est_temp->Clone("h_mass_temp_WJET")));
    h_mass_WJET_MC->SetMarkerStyle(kFullDotLarge);
    h_mass_SS_WJET->SetMarkerStyle(kFullDotLarge);
    h_mass_temp_WJET->SetMarkerStyle(kFullDotLarge);
    h_mass_WJET_MC->SetMarkerColor(kBlack);
    h_mass_SS_WJET->SetMarkerColor(kBlack);
    h_mass_temp_WJET->SetMarkerColor(kBlack);
    h_mass_WJET_MC->SetLineColor(kBlack);
    h_mass_SS_WJET->SetLineColor(kBlack);
    h_mass_temp_WJET->SetLineColor(kBlack);
    h_mass_SS_WJET_MC->SetMarkerStyle(kFullSquare);
    h_mass_temp_WJET_MC->SetMarkerStyle(kFullSquare);
    h_mass_SS_WJET_MC->SetMarkerColor(kRed);
    h_mass_temp_WJET_MC->SetMarkerColor(kRed);
    h_mass_SS_WJET_MC->SetLineColor(kRed);
    h_mass_temp_WJET_MC->SetLineColor(kRed);

    myRatioPlot_t *RP_mass_WJET_tSSvsOS = new myRatioPlot_t("RP_mass_WJET_tSSvsOS", h_mass_WJET_MC, h_mass_temp_WJET_MC);
    myRatioPlot_t *RP_mass_WJET_SSvsOS = new myRatioPlot_t("RP_mass_WJET_SSvsOS", h_mass_WJET_MC, h_mass_SS_WJET_MC);
    myRatioPlot_t *RP_mass_WJET_SSMCvsEst = new myRatioPlot_t("RP_mass_WJET_SSMCvsEst", h_mass_SS_WJET_MC, h_mass_SS_WJET);
    myRatioPlot_t *RP_mass_WJET_tSSMCvsEst = new myRatioPlot_t("RP_mass_WJET_SSMCvsEst", h_mass_temp_WJET_MC, h_mass_temp_WJET);
    RP_mass_WJET_tSSvsOS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "SS/(OS+SS)");
    RP_mass_WJET_SSvsOS->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "SS/(OS+SS)");
    RP_mass_WJET_SSMCvsEst->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "MC/Est");
    RP_mass_WJET_tSSMCvsEst->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 15, 3000, "MC/Est");

    TLegend *l_tSSvsOS = new TLegend(0.5, 0.8, 0.95, 0.95);
    TLegend *l_SSvsOS = new TLegend(0.5, 0.8, 0.95, 0.95);
    TLegend *l_MCvsEst = new TLegend(0.5, 0.8, 0.95, 0.95);
    TLegend *l_tSSMCvsEst = new TLegend(0.5, 0.8, 0.95, 0.95);
    l_tSSvsOS->AddEntry(h_mass_WJET_MC, "WJets MC", "lp");
    l_tSSvsOS->AddEntry(h_mass_temp_WJET_MC, "WJets MC (Tight same-sign template)", "lp");
    l_SSvsOS->AddEntry(h_mass_WJET_MC, "Any-sign W+Jets MC", "lp");
    l_SSvsOS->AddEntry(h_mass_SS_WJET_MC, "Same-sign W+Jets MC", "lp");
    l_MCvsEst->AddEntry(h_mass_SS_WJET, "Same-sign WJets (data)", "lp");
    l_MCvsEst->AddEntry(h_mass_SS_WJET_MC, "Same-sign WJets (MC)", "lp");
    l_tSSMCvsEst->AddEntry(h_mass_temp_WJET, "Tight same-sign WJets (data)", "lp");
    l_tSSMCvsEst->AddEntry(h_mass_temp_WJET_MC, "Tight same-sign WJets (MC)", "lp");
    RP_mass_WJET_tSSvsOS->ImportLegend(l_tSSvsOS);
    RP_mass_WJET_SSvsOS->ImportLegend(l_SSvsOS);
    RP_mass_WJET_SSMCvsEst->ImportLegend(l_MCvsEst);
    RP_mass_WJET_tSSMCvsEst->ImportLegend(l_tSSMCvsEst);

    RP_mass_WJET_tSSvsOS->Draw(5, 1e3, 1, "E");
    RP_mass_WJET_SSvsOS->Draw(5, 1e3, 1, "E");
    RP_mass_WJET_SSMCvsEst->Draw(5, 1e3, 1, "E");
    RP_mass_WJET_tSSMCvsEst->Draw(5, 1e3, 1, "E");

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (FR)", "f");

    // Draw WJets from simple subtraction
    TCanvas * c_WJET_est = new TCanvas("c_WJET_est", "W+Jets est", 750, 850);
    c_WJET_est->SetTopMargin(0.05);
    c_WJET_est->SetRightMargin(0.05);
    c_WJET_est->SetBottomMargin(0.15);
    c_WJET_est->SetLeftMargin(0.17);
    h_WJET_est->Draw("BAR");
    h_WJET_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]");
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

    // Draw WJets from simple subtraction (same-sign template)
    TCanvas * c_WJET_est_SS = new TCanvas("c_WJET_est_SS", "W+Jets est (same-sign template)", 750, 850);
    c_WJET_est_SS->SetTopMargin(0.05);
    c_WJET_est_SS->SetRightMargin(0.05);
    c_WJET_est_SS->SetBottomMargin(0.15);
    c_WJET_est_SS->SetLeftMargin(0.17);
    h_WJET_est_SS->Draw("BAR");
    h_WJET_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (same-sign) [GeV/c^{2}]");
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
    l_WJET_est->Draw();
    c_WJET_est_SS->SetLogx();
    c_WJET_est_SS->SetGridx();
    c_WJET_est_SS->SetGridy();
    c_WJET_est_SS->Update();

    // Draw WJets from simple subtraction (tight same-sign template)
    TCanvas * c_WJET_est_temp = new TCanvas("c_WJET_est_temp", "W+Jets est (tight same-sign template)", 750, 850);
    c_WJET_est_temp->SetTopMargin(0.05);
    c_WJET_est_temp->SetRightMargin(0.05);
    c_WJET_est_temp->SetBottomMargin(0.15);
    c_WJET_est_temp->SetLeftMargin(0.17);
    h_WJET_est_temp->Draw("BAR");
    h_WJET_est_temp->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} (tight same-sign) [GeV/c^{2}]");
    h_WJET_est_temp->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est_temp->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est_temp->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est_temp->GetXaxis()->SetMoreLogLabels();
    h_WJET_est_temp->GetXaxis()->SetNoExponent();
    h_WJET_est_temp->GetYaxis()->SetTitle("Number of events");
    h_WJET_est_temp->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est_temp->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est_temp->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est_temp->GetYaxis()->SetMoreLogLabels();
    h_WJET_est_temp->GetYaxis()->SetNoExponent();
    l_WJET_est->Draw();
    c_WJET_est_temp->SetLogx();
    c_WJET_est_temp->SetGridx();
    c_WJET_est_temp->SetGridy();
    c_WJET_est_temp->Update();

    // Draw WJets from same-sign template fit
    TCanvas * c_WJET_fit = new TCanvas("c_WJET_fit", "W+Jets fit (same-sign template)", 750, 850);
    c_WJET_fit->SetTopMargin(0.05);
    c_WJET_fit->SetRightMargin(0.05);
    c_WJET_fit->SetBottomMargin(0.15);
    c_WJET_fit->SetLeftMargin(0.17);
    h_WJET_est_fit->Draw("BAR");
    h_WJET_est_fit->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]");
    h_WJET_est_fit->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est_fit->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est_fit->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est_fit->GetXaxis()->SetMoreLogLabels();
    h_WJET_est_fit->GetXaxis()->SetNoExponent();
    h_WJET_est_fit->GetYaxis()->SetTitle("Number of events");
    h_WJET_est_fit->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est_fit->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est_fit->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est_fit->GetYaxis()->SetMoreLogLabels();
    h_WJET_est_fit->GetYaxis()->SetNoExponent();    

    // Draw WJets from tight same-sign template fit
    TCanvas * c_WJET_fit_2 = new TCanvas("c_WJET_fit_2", "W+Jets fit (tight same-sign template)", 750, 850);
    c_WJET_fit_2->SetTopMargin(0.05);
    c_WJET_fit_2->SetRightMargin(0.05);
    c_WJET_fit_2->SetBottomMargin(0.15);
    c_WJET_fit_2->SetLeftMargin(0.17);
    h_WJET_est_fit_2->Draw("BAR");
    h_WJET_est_fit_2->GetXaxis()->SetTitle("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]");
    h_WJET_est_fit_2->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est_fit_2->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est_fit_2->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est_fit_2->GetXaxis()->SetMoreLogLabels();
    h_WJET_est_fit_2->GetXaxis()->SetNoExponent();
    h_WJET_est_fit_2->GetYaxis()->SetTitle("Number of events");
    h_WJET_est_fit_2->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est_fit_2->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est_fit_2->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est_fit_2->GetYaxis()->SetMoreLogLabels();
    h_WJET_est_fit_2->GetYaxis()->SetNoExponent();
    l_WJET_est->Draw();
    c_WJET_fit_2->SetLogx();
    c_WJET_fit_2->SetGridx();
    c_WJET_fit_2->SetGridy();
    c_WJET_fit_2->Update();

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_WJET_est->Write();
    h_WJET_est_SS->Write();
    h_WJET_est_temp->Write();
    h_WJET_est_fit->Write();
    h_WJET_est_fit_2->Write();

    if (systErr > 0) // UPDATE FOR fit_2
    {   // Errors
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_DOWN.root", "READ");
        TFile *f_alt = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE_PFiso.root", "READ");
        TH1D *h_up, *h_down, *h_up_temp, *h_down_temp, *h_alt, *h_alt_temp, *h_tempVSsub, *h_fitUncer,
             *h_fullsysterr, *h_fullsysterr_temp, *h_fullerr, *h_fullerr_temp;
        f_up->GetObject("h_WJET_est", h_up);
        f_up->GetObject("h_WJET_est_fit", h_up_temp);
        f_down->GetObject("h_WJET_est", h_down);
        f_down->GetObject("h_WJET_est_fit", h_down_temp);
        f_alt->GetObject("h_WJET_est", h_alt);
        f_alt->GetObject("h_WJET_est_fit", h_alt_temp);
        h_tempVSsub = ((TH1D*)(h_WJET_est->Clone("h_WJET_tempVSsub")));
        h_fitUncer = ((TH1D*)(h_WJET_est_fit->Clone("h_fitUncer")));
        h_fullsysterr = ((TH1D*)(h_up->Clone("h_WJET_fullsysterr")));
        h_fullerr = ((TH1D*)(h_up->Clone("h_WJET_fullerr")));
        h_fullsysterr_temp = ((TH1D*)(h_up_temp->Clone("h_WJET_fullsysterr_temp")));
        h_fullerr_temp = ((TH1D*)(h_up_temp->Clone("h_WJET_fullerr_temp")));

        for (Int_t i=1; i<h_up->GetSize()-1; i++)
        {
            // Systematic errors
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_up_temp->SetBinContent(i, fabs(h_up_temp->GetBinContent(i)-h_down_temp->GetBinContent(i))/2);
            h_alt->SetBinContent(i, fabs(h_alt->GetBinContent(i)-h_WJET_est->GetBinContent(i)));
            h_alt_temp->SetBinContent(i, fabs(h_alt_temp->GetBinContent(i)-h_WJET_est_fit->GetBinContent(i)));
            h_tempVSsub->SetBinContent(i, fabs(h_tempVSsub->GetBinContent(i)-h_WJET_est_fit->GetBinContent(i)));
            h_fitUncer->SetBinContent(i, 0.01443*h_fitUncer->GetBinContent(i));
            Double_t err = h_up->GetBinContent(i) * h_up->GetBinContent(i) + /*h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +*/
                           h_alt->GetBinContent(i) * h_alt->GetBinContent(i);
            Double_t err_temp = h_up_temp->GetBinContent(i) * h_up_temp->GetBinContent(i) + h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +
                                h_alt_temp->GetBinContent(i) * h_alt_temp->GetBinContent(i) + h_fitUncer->GetBinContent(i) * h_fitUncer->GetBinContent(i);
            if (sqrt(err) < h_WJET_est->GetBinContent(i))
                h_fullsysterr->SetBinContent(i, sqrt(err));
            else
                h_fullsysterr->SetBinContent(i, h_WJET_est->GetBinContent(i));

            if (sqrt(err_temp) < h_WJET_est_fit->GetBinContent(i))
                h_fullsysterr_temp->SetBinContent(i, sqrt(err_temp));
            else
                h_fullsysterr_temp->SetBinContent(i, h_WJET_est_fit->GetBinContent(i));
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_WJET_est->GetBinError(i)*h_WJET_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
            h_fullerr_temp->SetBinContent(i, sqrt(h_WJET_est_fit->GetBinError(i)*h_WJET_est_fit->GetBinError(i)+
                                                  h_fullsysterr_temp->GetBinContent(i)*h_fullsysterr_temp->GetBinContent(i)));
        }
        cout << "Estimated W+Jets events (from template fit): " << int_wjet_fit << "+-" << err_wjet_fit << "+-" <<
                h_fullsysterr_temp->Integral() << "   (+-" << h_fullerr_temp->Integral() << ")" << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_WJET_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_WJET_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetLineColor(39);
        h_draw_2->SetLineColor(39);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);

        TH1D *h_draw_1_temp = ((TH1D*)(h_WJET_est_fit->Clone("h_draw_1_temp")));
        TH1D *h_draw_2_temp = ((TH1D*)(h_WJET_est_fit->Clone("h_draw_2_temp")));
        h_draw_1_temp->Add(h_fullerr_temp, 1);
        h_draw_2_temp->Add(h_fullerr_temp, -1);
        h_draw_1_temp->SetMarkerStyle(0);
        h_draw_1_temp->SetFillColor(38);
        h_draw_1_temp->SetFillStyle(3244);
        h_draw_1_temp->SetLineColor(39);
        h_draw_2_temp->SetLineColor(39);
        h_draw_1_temp->SetDirectory(0);
        h_draw_2_temp->SetDirectory(0);

        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);
        h_fullsysterr_temp->SetDirectory(0);
        h_fullerr_temp->SetDirectory(0);
        h_tempVSsub->SetDirectory(0);

        f_up->Close();
        f_down->Close();
        f_alt->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullsysterr_temp->Write();
        h_fullerr->Write();
        h_fullerr_temp->Write();

        c_WJET_est->cd();
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
        c_WJET_est->Update();

        c_WJET_fit->cd();
        h_draw_1_temp->Draw("samehist");
        h_draw_2_temp->Draw("samehist");
        c_WJET_fit->Update();
    }// End of if (systErr)
    l_WJET_est->Draw();
    c_WJET_fit->SetLogx();
    c_WJET_fit->SetGridx();
    c_WJET_fit->SetGridy();
    c_WJET_fit->Update();

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

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special], *h_pT_lead[_EndOf_Data_Special],
         *h_pT_sublead[_EndOf_Data_Special], *h_HLT_pT[_EndOf_Data_Special];
    TH1D *h_QCD_est, *h_QCD_est_SS;
    THStack * s_mass_wQCD = new THStack("s_mass_wQCD", "");
    THStack * s_mass_woQCD = new THStack("s_mass_woQCD", "");
    THStack * s_mass_SS_wQCD = new THStack("s_mass_SS_wQCD", "");
    THStack * s_mass_SS_woQCD = new THStack("s_mass_SS_woQCD", "");
    THStack * s_pT_lead = new THStack("s_pT_lead", "");
    THStack * s_pT_sublead = new THStack("s_pT_sublead", "");
    THStack * s_HLT_pT = new THStack("s_HLT_pT", "");
    Color_t color = kBlack;

    // Loop over all processes (adding all histograms)
    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        f->GetObject("h_pT_lead_"+Mgr.Procname[pr], h_pT_lead[pr]);
        f->GetObject("h_pT_sublead_"+Mgr.Procname[pr], h_pT_sublead[pr]);
        f->GetObject("h_HLT_pT_"+Mgr.Procname[pr], h_HLT_pT[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
            removeNegativeBins(h_pT_lead[pr]);
            removeNegativeBins(h_pT_sublead[pr]);
            removeNegativeBins(h_HLT_pT[pr]);
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
            h_pT_lead[pr]->SetFillColor(color);
            h_pT_lead[pr]->SetLineColor(color);
            h_pT_sublead[pr]->SetFillColor(color);
            h_pT_sublead[pr]->SetLineColor(color);
            h_HLT_pT[pr]->SetFillColor(color);
            h_HLT_pT[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
            h_pT_lead[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_lead[pr]->SetMarkerColor(kBlack);
            h_pT_lead[pr]->SetLineColor(kBlack);
            h_pT_sublead[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_sublead[pr]->SetMarkerColor(kBlack);
            h_pT_sublead[pr]->SetLineColor(kBlack);
            h_HLT_pT[pr]->SetMarkerStyle(kFullDotLarge);
            h_HLT_pT[pr]->SetMarkerColor(kBlack);
            h_HLT_pT[pr]->SetLineColor(kBlack);
        }

        // Adding hists to THStacks
        if (pr < _SingleMuon_B)
        {
            s_mass_wQCD->Add(h_mass[pr]);
            s_mass_SS_wQCD->Add(h_mass_SS[pr]);
            s_pT_lead->Add(h_pT_lead[pr]);
            s_pT_sublead->Add(h_pT_sublead[pr]);
            s_HLT_pT->Add(h_HLT_pT[pr]);
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
            h_pT_lead[_DY_Full] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_DY_Full")));
            h_pT_lead[_DY_Full]->SetDirectory(0);
            h_pT_sublead[_DY_Full] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_DY_Full")));
            h_pT_sublead[_DY_Full]->SetDirectory(0);
            h_HLT_pT[_DY_Full] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_DY_Full")));
            h_HLT_pT[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
            h_pT_lead[_DY_Full]->Add(h_pT_lead[pr]);
            h_pT_sublead[_DY_Full]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_DY_Full]->Add(h_HLT_pT[pr]);
        }
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
            h_pT_lead[_ttbar_Full] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_ttbar_Full")));
            h_pT_lead[_ttbar_Full]->SetDirectory(0);
            h_pT_sublead[_ttbar_Full] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_ttbar_Full")));
            h_pT_sublead[_ttbar_Full]->SetDirectory(0);
            h_HLT_pT[_ttbar_Full] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_ttbar_Full")));
            h_HLT_pT[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
            h_pT_lead[_ttbar_Full]->Add(h_pT_lead[pr]);
            h_pT_sublead[_ttbar_Full]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_ttbar_Full]->Add(h_HLT_pT[pr]);
        }
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
            h_pT_lead[_VVnST] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_VVnST")));
            h_pT_lead[_VVnST]->SetDirectory(0);
            h_pT_sublead[_VVnST] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_VVnST")));
            h_pT_sublead[_VVnST]->SetDirectory(0);
            h_HLT_pT[_VVnST] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_VVnST")));
            h_HLT_pT[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
            h_pT_lead[_VVnST]->Add(h_pT_lead[pr]);
            h_pT_sublead[_VVnST]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_VVnST]->Add(h_HLT_pT[pr]);
        }
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
            h_pT_lead[_WJets_Full] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_WJets_Full")));
            h_pT_lead[_WJets_Full]->SetDirectory(0);
            h_pT_sublead[_WJets_Full] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_WJets_Full")));
            h_pT_sublead[_WJets_Full]->SetDirectory(0);
            h_HLT_pT[_WJets_Full] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_WJets_Full")));
            h_HLT_pT[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
            h_pT_lead[_WJets_Full]->Add(h_pT_lead[pr]);
            h_pT_sublead[_WJets_Full]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_WJets_Full]->Add(h_HLT_pT[pr]);
        }
        else if (pr == _QCDMuEnriched_15to20)
        {
            h_mass[_QCDMuEnriched_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_QCDMuEnriched_Full")));
            h_mass[_QCDMuEnriched_Full]->SetDirectory(0);
            h_mass_SS[_QCDMuEnriched_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_QCDMuEnriched_Full")));
            h_mass_SS[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_lead[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_QCDMuEnriched_Full")));
            h_pT_lead[_QCDMuEnriched_Full]->SetDirectory(0);
            h_pT_sublead[_QCDMuEnriched_Full] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_QCDMuEnriched_Full")));
            h_pT_sublead[_QCDMuEnriched_Full]->SetDirectory(0);
            h_HLT_pT[_QCDMuEnriched_Full] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_QCDMuEnriched_Full")));
            h_HLT_pT[_QCDMuEnriched_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_QCDMuEnriched_Normal)
        {
            h_mass[_QCDMuEnriched_Full]->Add(h_mass[pr]);
            h_mass_SS[_QCDMuEnriched_Full]->Add(h_mass_SS[pr]);
            h_pT_lead[_QCDMuEnriched_Full]->Add(h_pT_lead[pr]);
            h_pT_sublead[_QCDMuEnriched_Full]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_QCDMuEnriched_Full]->Add(h_HLT_pT[pr]);
        }
        else if (pr == _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_SingleMuon_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
            h_mass_SS[_SingleMuon_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_SingleMuon_Full")));
            h_mass_SS[_SingleMuon_Full]->SetDirectory(0);
            h_pT_lead[_SingleMuon_Full] = ((TH1D*)(h_pT_lead[pr]->Clone("h_pT_lead_SingleMuon_Full")));
            h_pT_lead[_SingleMuon_Full]->SetDirectory(0);
            h_pT_sublead[_SingleMuon_Full] = ((TH1D*)(h_pT_sublead[pr]->Clone("h_pT_sublead_SingleMuon_Full")));
            h_pT_sublead[_SingleMuon_Full]->SetDirectory(0);
            h_HLT_pT[_SingleMuon_Full] = ((TH1D*)(h_HLT_pT[pr]->Clone("h_HLT_pT_SingleMuon_Full")));
            h_HLT_pT[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_SingleMuon_Normal)
        {
            h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
            h_mass_SS[_SingleMuon_Full]->Add(h_mass_SS[pr]);
            h_pT_lead[_SingleMuon_Full]->Add(h_pT_lead[pr]);
            h_pT_sublead[_SingleMuon_Full]->Add(h_pT_sublead[pr]);
            h_HLT_pT[_SingleMuon_Full]->Add(h_HLT_pT[pr]);
        }


        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

    } // End of pr iteration

    // QCD estimation
    h_QCD_est = ((TH1D*)(h_mass[_SingleMuon_Full]->Clone("h_QCD_est")));
    Double_t err_data=0, int_data=0, err_qcd=0, int_qcd=0, err_dy=0, int_dy=0, err_tt=0, int_tt=0, err_vvnst=0, int_vvnst=0, err_wj=0, int_wj=0;
    int_data = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_data);
    h_QCD_est->SetTitle("");
    h_QCD_est->SetDirectory(0);
    h_QCD_est->Add(h_mass[_DY_Full], -1);
    h_QCD_est->Add(h_mass[_ttbar_Full], -1);
    h_QCD_est->Add(h_mass[_VVnST], -1);
    h_QCD_est->Add(h_mass[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    // Setting proper errors
    for (Int_t i_bin=1; i_bin<h_QCD_est->GetSize()-1; i_bin++)
    {
        h_QCD_est->SetBinError(i_bin, sqrt(h_QCD_est->GetBinContent(i_bin)));
        if (h_QCD_est->GetBinContent(i_bin) == 0)
            h_QCD_est->SetBinError(i_bin, 1);
    }
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kBlack);
    h_QCD_est->SetMarkerStyle(0);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    int_dy = h_mass[_DY_Full]->IntegralAndError(1, h_mass[_DY_Full]->GetSize()-2, err_dy);
    int_tt = h_mass[_ttbar_Full]->IntegralAndError(1, h_mass[_ttbar_Full]->GetSize()-2, err_tt);
    int_vvnst = h_mass[_VVnST]->IntegralAndError(1, h_mass[_VVnST]->GetSize()-2, err_vvnst);
    int_wj = h_mass[_WJets_Full]->IntegralAndError(1, h_mass[_WJets_Full]->GetSize()-2, err_wj);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "QCD est events: " << int_qcd << "+-" << err_qcd << endl;
    cout << "DY events: " << int_dy << "+-" << err_dy << "(" << int_dy/int_data << ")" << endl;
    cout << "ttbar events: " << int_tt << "+-" << err_tt << "(" << int_tt/int_data << ")" << endl;
    cout << "VVnST events: " << int_vvnst << "+-" << err_vvnst << "(" << int_vvnst/int_data << ")" << endl;
    cout << "WJets events: " << int_wj << "+-" << err_wj << "(" << int_wj/int_data << ")" << endl;

    Double_t err_data_ss=0, int_data_ss=0, err_qcd_ss=0, int_qcd_ss=0, err_dy_ss=0, int_dy_ss=0, err_tt_ss=0, int_tt_ss=0,
            err_vvnst_ss=0, int_vvnst_ss=0, err_wj_ss=0, int_wj_ss=0;
    h_QCD_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_QCD_est_SS")));
    int_data_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_data_ss);
    h_QCD_est_SS->SetTitle("");
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_QCD_est_SS->Add(h_mass_SS[_WJets_Full], -1);
    removeNegativeBins(h_QCD_est_SS);
    h_QCD_est_SS->SetFillColor(kRed + 3);
    h_QCD_est_SS->SetLineColor(kRed + 3);
    int_qcd_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_qcd_ss);
    int_dy_ss = h_mass_SS[_DY_Full]->IntegralAndError(1, h_mass_SS[_DY_Full]->GetSize()-2, err_dy_ss);
    int_tt_ss = h_mass_SS[_ttbar_Full]->IntegralAndError(1, h_mass_SS[_ttbar_Full]->GetSize()-2, err_tt_ss);
    int_vvnst_ss = h_mass_SS[_VVnST]->IntegralAndError(1, h_mass_SS[_VVnST]->GetSize()-2, err_vvnst_ss);
    int_wj_ss = h_mass_SS[_WJets_Full]->IntegralAndError(1, h_mass_SS[_WJets_Full]->GetSize()-2, err_wj_ss);
    cout << "Data events (same-sign): " << int_data_ss << "+-" << err_data_ss << endl;
    cout << "QCD est events (same-sign): " << int_qcd_ss << "+-" << err_qcd_ss << endl;
    cout << "DY events (same-sign): " << int_dy_ss << "+-" << err_dy_ss << "(" << int_dy_ss/int_data_ss << ")" << endl;
    cout << "ttbar events (same-sign): " << int_tt_ss << "+-" << err_tt_ss << "(" << int_tt_ss/int_data_ss << ")" << endl;
    cout << "VVnST events (same-sign): " << int_vvnst_ss << "+-" << err_vvnst_ss << "(" << int_vvnst_ss/int_data_ss << ")" << endl;
    cout << "WJets events (same-sign): " << int_wj_ss << "+-" << err_wj_ss << "(" << int_wj_ss/int_data_ss << ")" << endl;

    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass_wQCD = new myRatioPlot_t("c_mass_wQCD", s_mass_wQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_wQCD = new myRatioPlot_t("c_mass_SS_wQCD", s_mass_SS_wQCD, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_SS_woQCD = new myRatioPlot_t("c_mass_SS_woQCD", s_mass_SS_woQCD, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_lead = new myRatioPlot_t("c_pT_lead", s_pT_lead, h_pT_lead[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_sublead = new myRatioPlot_t("c_pT_sublead", s_pT_sublead, h_pT_sublead[_SingleMuon_Full]);
    myRatioPlot_t *RP_HLT_pT = new myRatioPlot_t("c_HLT_pT", s_HLT_pT, h_HLT_pT[_SingleMuon_Full]);

    RP_mass_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_wQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_SS_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} (same sign) [GeV/c^{2}]", 15, 3000);
    RP_pT_lead->SetPlots("p_{#lower[-0.2]{T}} (#mu_{lead}) [GeV/c]", 0, 1000);
    RP_pT_sublead->SetPlots("p_{#lower[-0.2]{T}} (#mu_{sublead}) [GeV/c]", 0, 1000);
    RP_HLT_pT->SetPlots("HLT object p_{#lower[-0.2]{T}} [GeV/c]", 20, 200);

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
    RP_pT_lead->ImportLegend(legend_wQCD);
    RP_pT_sublead->ImportLegend(legend_wQCD);
    RP_HLT_pT->ImportLegend(legend_wQCD);

    RP_mass_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_woQCD->Draw(1e-3, 1e3, 1);
    RP_mass_SS_wQCD->Draw(1e-3, 1e3, 1);
    RP_mass_SS_woQCD->Draw(1e-3, 1e3, 1);
    RP_pT_lead->Draw(1e-1, 1e5, 0);
    RP_pT_sublead->Draw(1e-1, 1e5, 0);
    RP_HLT_pT->Draw(10, 1e7, 0);

    // Drawing estimated QCD
    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (FR)", "f");
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
        TFile *f_up = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstQCD_MuMu_DOWN.root", "READ");
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
        cout << "Estimated QCD events: " << int_qcd << "+-" << err_qcd << "+-" << h_fullsysterr->Integral() << "   (+-" <<
                h_fullerr->Integral() << ")" << endl;
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
    Double_t int_data=0, err_data=0, int_wjet=0, err_wjet=0, int_dy=0, err_dy=0, int_tt=0, err_tt=0, int_vvnst=0, err_vvnst=0, int_qcd=0, err_qcd=0;
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
    int_dy = h_mass[_DY_Full]->IntegralAndError(1, h_mass[_DY_Full]->GetSize()-2, err_dy);
    int_tt = h_mass[_ttbar_Full]->IntegralAndError(1, h_mass[_ttbar_Full]->GetSize()-2, err_tt);
    int_vvnst = h_mass[_VVnST]->IntegralAndError(1, h_mass[_VVnST]->GetSize()-2, err_vvnst);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "WJets est events (subraction): " << int_wjet << "+-" << err_wjet << endl;
    cout << "DY events: " << int_dy << "+-" << err_dy << "(" << int_dy/int_data << ")" << endl;
    cout << "ttbar events: " << int_tt << "+-" << err_tt << "(" << int_tt/int_data << ")" << endl;
    cout << "VVnST events: " << int_vvnst << "+-" << err_vvnst << "(" << int_vvnst/int_data << ")" << endl;
    cout << "QCD events: " << int_qcd << "+-" << err_qcd << "(" << int_qcd/int_data << ")" << endl;

    Double_t int_data_ss=0, err_data_ss=0, int_wjet_ss=0, err_wjet_ss=0, int_dy_ss=0, err_dy_ss=0, int_tt_ss=0, err_tt_ss=0,
            int_vvnst_ss=0, err_vvnst_ss=0, int_qcd_ss=0, err_qcd_ss=0;
    h_WJET_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_WJET_est_SS")));
    int_data_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_data_ss);
    h_WJET_est_SS->SetTitle("");
    h_WJET_est_SS->SetDirectory(0);
    h_WJET_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
//    h_WJET_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_WJET_est_SS->Add(h_QCD_est_SS, -1);
    removeNegativeBins(h_WJET_est_SS);
    h_WJET_est_SS->SetFillColor(kRed - 2);
    h_WJET_est_SS->SetLineColor(kRed - 2);
    int_wjet_ss = h_WJET_est_SS->IntegralAndError(1, h_WJET_est_SS->GetSize()-2, err_wjet_ss);
    int_dy_ss = h_mass_SS[_DY_Full]->IntegralAndError(1, h_mass_SS[_DY_Full]->GetSize()-2, err_dy_ss);
    int_tt_ss = h_mass_SS[_ttbar_Full]->IntegralAndError(1, h_mass_SS[_ttbar_Full]->GetSize()-2, err_tt_ss);
    int_vvnst_ss = h_mass_SS[_VVnST]->IntegralAndError(1, h_mass_SS[_VVnST]->GetSize()-2, err_vvnst_ss);
    int_qcd_ss = h_QCD_est_SS->IntegralAndError(1, h_QCD_est_SS->GetSize()-2, err_qcd_ss);
    cout << "Data SS events: " << int_data_ss << "+-" << err_data_ss << endl;
    cout << "WJets SS est events (subraction): " << int_wjet_ss << "+-" << err_wjet_ss << endl;
    cout << "DY SS events: " << int_dy_ss << "+-" << err_dy_ss << "(" << int_dy_ss/int_data_ss << ")" << endl;
    cout << "ttbar SS events: " << int_tt_ss << "+-" << err_tt_ss << "(" << int_tt_ss/int_data_ss << ")" << endl;
    cout << "VVnST SS events: " << int_vvnst_ss << "+-" << err_vvnst_ss << "(" << int_vvnst_ss/int_data_ss << ")" << endl;
    cout << "QCD SS events: " << int_qcd_ss << "+-" << err_qcd_ss << "(" << int_qcd_ss/int_data_ss << ")" << endl;

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
    TLegend * legend_fit = new TLegend(0.7, 0.64, 0.95, 0.95);
    legend_fit->AddEntry(h_mass[_SingleMuon_B], "Data", "pl");
    legend_fit->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend_fit->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend_fit->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend_fit->AddEntry(h_mass[_tbarW], "#font[12]{#scale[1.1]{tW+#bar{t}W+VV}}", "f");
    legend_fit->AddEntry(h_mass[_QCDMuEnriched_15to20], "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);
    RP_mass_SS_wWJET->ImportLegend(legend_wWJET);
    RP_mass_SS_woWJET->ImportLegend(legend_woWJET);

    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);
    RP_mass_SS_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_SS_woWJET->Draw(1e-2, 1e5, 1);

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (FR)", "f");

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
//    h_WJET_fit_MC->Scale(3.8159e+03 / h_WJET_fit_MC->Integral()); // pT cuts: 52, 0    
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
    cout << "W+Jets fit factor: " << 3.1476e+03 / h_WJET_fit_SS->Integral() << endl;
    cout << "Same sign W+Jets integral: " << h_WJET_fit_SS->Integral(1,30) << endl;
//    h_WJET_fit_SS->Scale(3.5028e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (my FR 2)
    h_WJET_fit_SS->Scale(3.1476e+03 / h_WJET_fit_SS->Integral()); // from fitting with 5 GeV bins (my FR 2)
//    h_WJET_fit_SS->Scale(3.2959e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS)
//    h_WJET_fit_SS->Scale(3.1310e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS (up))
//    h_WJET_fit_SS->Scale(3.1729e+03 / h_WJET_fit_SS->Integral(1,30)); // from fitting with 5 GeV bins (FOR SYSTEMATIC ERRORS (down))
    // Setting proper errors
    for (Int_t i_bin=1; i_bin < h_WJET_fit_SS->GetSize()-1; i_bin++)
    {
        h_WJET_fit_SS->SetBinError(i_bin, sqrt(h_WJET_fit_SS->GetBinContent(i_bin)));
        if (h_WJET_fit_SS->GetBinContent(i_bin) == 0)
            h_WJET_fit_SS->SetBinError(i_bin, 1);
    }
    Double_t int_wjet_est, err_wjet_est;
    int_wjet_est = h_WJET_fit_SS->IntegralAndError(1, h_WJET_fit_SS->GetSize()-2, err_wjet_est);
    cout << "WJets est events (FROM SS): " << int_wjet_est << "+-" << err_wjet_est << endl;
    h_WJET_fit_SS->SetDirectory(0);
    h_WJET_fit_SS->SetFillColor(kRed - 2);
    h_WJET_fit_SS->SetLineColor(kBlack);
    h_WJET_fit_SS->SetMarkerStyle(0);

    // Full stack histogram from template fit
    Double_t int_wjet_fit=0, err_wjet_fit=0;
    TH1D *h_DY_fit = ((TH1D*)(h_mass[_DY_Full]->Clone("h_DY_fit")));
    TH1D *h_ttbar_fit = ((TH1D*)(h_mass[_ttbar_Full]->Clone("h_ttbar_fit")));
    TH1D *h_VVnST_fit = ((TH1D*)(h_mass[_VVnST]->Clone("h_VVnST_fit")));
    TH1D *h_QCD_fit = ((TH1D*)(h_QCD_est->Clone("h_QCD_fit")));
    cout << "DY fit factor: " << 3.0397e+05 / h_DY_fit->Integral() << endl;
    cout << "ttbar fit factor: " << 1.0230e+04 / h_ttbar_fit->Integral() << endl;
    cout << "VVnST fit factor: " << 2.2971e+03 / h_VVnST_fit->Integral() << endl;
    cout << "QCD fit factor: " << 6.7421e+03 / h_QCD_fit->Integral() << endl;
//    h_DY_fit->Scale(3.0356e+05 / h_DY_fit->Integral(1,30));
//    h_ttbar_fit->Scale(3.9280e+03 / h_ttbar_fit->Integral(1,30));
//    h_VVnST_fit->Scale(2.2021e+03 / h_VVnST_fit->Integral(1,30));
//    h_QCD_fit->Scale(5.9731e+03 / h_QCD_fit->Integral(1,30));
    h_DY_fit->Scale(3.0397e+05 / h_DY_fit->Integral());
    h_ttbar_fit->Scale(1.0230e+04 / h_ttbar_fit->Integral());
    h_VVnST_fit->Scale(2.2971e+03 / h_VVnST_fit->Integral());
    h_QCD_fit->Scale(6.7421e+03 / h_QCD_fit->Integral());
    h_VVnST_fit->SetFillColor(kGreen+2);
    h_VVnST_fit->SetLineColor(kGreen+2);
    THStack *s_mass_fit = new THStack("s_mass_fit", "");
    s_mass_fit->Add(h_QCD_fit);
    s_mass_fit->Add(h_VVnST_fit);
    s_mass_fit->Add(h_ttbar_fit);
    s_mass_fit->Add(h_WJET_fit_SS);
    s_mass_fit->Add(h_DY_fit);
    int_wjet_fit = h_WJET_fit_SS->IntegralAndError(1, h_WJET_fit_SS->GetSize()-2, err_wjet_fit);
    cout << "WJets est events (template fit): " << int_wjet_fit << "+-" << err_wjet_fit << endl;

    myRatioPlot_t *RP_mass_fit = new myRatioPlot_t("c_mass_fit", s_mass_fit, h_mass[_SingleMuon_Full]);
    RP_mass_fit->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#mu#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_fit->ImportLegend(legend_fit);
    RP_mass_fit->Draw(1e-2, 1e5, 1);

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
        TFile *f_up = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu_DOWN.root", "READ");
        TH1D *h_up, *h_down, *h_up_temp, *h_down_temp, *h_ratio, *h_ratio_temp, *h_tempVSsub, *h_fitUncer,
             *h_fullsysterr, *h_fullsysterr_temp, *h_fullerr, *h_fullerr_temp;
        f_up->GetObject("h_WJET_est", h_up);
        f_up->GetObject("h_WJET_fit_SS", h_up_temp);
        f_down->GetObject("h_WJET_est", h_down);
        f_down->GetObject("h_WJET_fit_SS", h_down_temp);
        f_ratio->GetObject("h_WJET_est", h_ratio);
        f_ratio->GetObject("h_WJET_fit_SS", h_ratio_temp);
        h_tempVSsub = ((TH1D*)(h_WJET_est->Clone("h_WJET_tempVSsub")));
        h_fitUncer = ((TH1D*)(h_WJET_fit_SS->Clone("h_fitUncer")));
        h_fullsysterr = ((TH1D*)(h_up->Clone("h_WJET_fullsysterr")));
        h_fullerr = ((TH1D*)(h_up->Clone("h_WJET_fullerr")));
        h_fullsysterr_temp = ((TH1D*)(h_up_temp->Clone("h_WJET_fullsysterr_temp")));
        h_fullerr_temp = ((TH1D*)(h_up_temp->Clone("h_WJET_fullerr_temp")));

        for (Int_t i=1; i<h_ratio->GetSize()-1; i++)
        {
            // Systematic errors
            h_up->SetBinContent(i, fabs(h_up->GetBinContent(i)-h_down->GetBinContent(i))/2);
            h_up_temp->SetBinContent(i, fabs(h_up_temp->GetBinContent(i)-h_down_temp->GetBinContent(i))/2);
            h_ratio->SetBinContent(i, fabs(h_ratio->GetBinContent(i)-h_WJET_est->GetBinContent(i)));
            h_ratio_temp->SetBinContent(i, fabs(h_ratio_temp->GetBinContent(i)-h_WJET_fit_SS->GetBinContent(i)));
            h_tempVSsub->SetBinContent(i, fabs(h_tempVSsub->GetBinContent(i)-h_WJET_fit_SS->GetBinContent(i)));
            h_fitUncer->SetBinContent(i, 0.10929*h_fitUncer->GetBinContent(i));
            Double_t err = h_up->GetBinContent(i) * h_ratio->GetBinContent(i) + /*h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +*/
                           h_ratio->GetBinContent(i) * h_ratio->GetBinContent(i);
            Double_t err_temp = h_up_temp->GetBinContent(i) * h_up_temp->GetBinContent(i) + h_tempVSsub->GetBinContent(i) * h_tempVSsub->GetBinContent(i) +
                                h_ratio_temp->GetBinContent(i) * h_ratio_temp->GetBinContent(i) + h_fitUncer->GetBinContent(i) * h_fitUncer->GetBinContent(i);
            if (sqrt(err) < h_WJET_est->GetBinContent(i))
                h_fullsysterr->SetBinContent(i, sqrt(err));
            else
                h_fullsysterr->SetBinContent(i, h_WJET_est->GetBinContent(i));

            if (sqrt(err_temp) < h_WJET_fit_SS->GetBinContent(i))
                h_fullsysterr_temp->SetBinContent(i, sqrt(err_temp));
            else
                h_fullsysterr_temp->SetBinContent(i, h_WJET_fit_SS->GetBinContent(i));
            // Combined error
            h_fullerr->SetBinContent(i, sqrt(h_WJET_est->GetBinError(i)*h_WJET_est->GetBinError(i)+
                                             h_fullsysterr->GetBinContent(i)*h_fullsysterr->GetBinContent(i)));
            h_fullerr_temp->SetBinContent(i, sqrt(h_WJET_fit_SS->GetBinError(i)*h_WJET_fit_SS->GetBinError(i)+
                                                  h_fullsysterr_temp->GetBinContent(i)*h_fullsysterr_temp->GetBinContent(i)));
        }
        cout << "Estimated W+Jets events (from template fit): " << int_wjet_est << "+-" << err_wjet_est << "+-" <<
                h_fullsysterr_temp->Integral() << "   (+-" << h_fullerr_temp->Integral() << ")" << endl;
        TH1D *h_draw_1 = ((TH1D*)(h_WJET_est->Clone("h_draw_1")));
        TH1D *h_draw_2 = ((TH1D*)(h_WJET_est->Clone("h_draw_2")));
        h_draw_1->Add(h_fullerr, 1);
        h_draw_2->Add(h_fullerr, -1);
        h_draw_1->SetMarkerStyle(0);
        h_draw_1->SetFillColor(38);
        h_draw_1->SetFillStyle(3244);
        h_draw_1->SetLineColor(39);
        h_draw_2->SetLineColor(39);
        h_draw_1->SetDirectory(0);
        h_draw_2->SetDirectory(0);

        TH1D *h_draw_1_temp = ((TH1D*)(h_WJET_fit_SS->Clone("h_draw_1_temp")));
        TH1D *h_draw_2_temp = ((TH1D*)(h_WJET_fit_SS->Clone("h_draw_2_temp")));
        h_draw_1_temp->Add(h_fullerr_temp, 1);
        h_draw_2_temp->Add(h_fullerr_temp, -1);
        h_draw_1_temp->SetMarkerStyle(0);
        h_draw_1_temp->SetFillColor(38);
        h_draw_1_temp->SetFillStyle(3244);
        h_draw_1_temp->SetLineColor(39);
        h_draw_2_temp->SetLineColor(39);
        h_draw_1_temp->SetDirectory(0);
        h_draw_2_temp->SetDirectory(0);

        h_fullsysterr->SetDirectory(0);
        h_fullerr->SetDirectory(0);
        h_fullsysterr_temp->SetDirectory(0);
        h_fullerr_temp->SetDirectory(0);
        h_tempVSsub->SetDirectory(0);

        f_up->Close();
        f_down->Close();
        f_ratio->Close();
        f_out->cd();
        h_fullsysterr->Write();
        h_fullsysterr_temp->Write();
        h_fullerr->Write();
        h_fullerr_temp->Write();

        c_WJET_est->cd();
        h_draw_1->Draw("samehist");
        h_draw_2->Draw("samehist");
        c_WJET_est->Update();

        c_WJET_fit_SS->cd();
        h_draw_1_temp->Draw("samehist");
        h_draw_2_temp->Draw("samehist");
        c_WJET_fit_SS->Update();

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



void EMu_QCDest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/EMu/";

    f = new TFile(Dir+"QCDest_EMu.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special], *h_mass_temp[_EndOf_Data_Special], *h_mass_test[_EndOf_Data_Special],
         *h_pT_e[_EndOf_Data_Special], *h_pT_mu[_EndOf_Data_Special];
    TH1D *h_QCD_est, *h_QCD_est_SS, *h_QCD_est_temp;
    THStack * s_mass_woQCD = new THStack("s_mass_woQCD", "");
    THStack * s_mass_woQCD_SS = new THStack("s_mass_woQCD_SS", "");
    THStack * s_mass_woQCD_temp = new THStack("s_mass_woQCD_template", "");
    THStack * s_mass_woQCD_test = new THStack("s_mass_woQCD_test", "");
    THStack * s_pT_e = new THStack("s_pT_e", "");
    THStack * s_pT_mu = new THStack("s_pT_mu", "");
    Color_t color = kBlack;

    // Loop over all processes (adding all histograms)
    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        f->GetObject("h_mass_template_"+Mgr.Procname[pr], h_mass_temp[pr]);
        f->GetObject("h_mass_test_"+Mgr.Procname[pr], h_mass_test[pr]);
        f->GetObject("h_pT_e_"+Mgr.Procname[pr], h_pT_e[pr]);
        f->GetObject("h_pT_mu_"+Mgr.Procname[pr], h_pT_mu[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);
        h_mass_temp[pr]->SetDirectory(0);
        h_mass_test[pr]->SetDirectory(0);
        h_pT_e[pr]->SetDirectory(0);
        h_pT_mu[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
            removeNegativeBins(h_mass_temp[pr]);
            removeNegativeBins(h_mass_test[pr]);
            removeNegativeBins(h_pT_e[pr]);
            removeNegativeBins(h_pT_mu[pr]);
        }

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_GJets_Normal) color = kYellow + 3;

        if (pr < _DoubleEG_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
            h_mass_temp[pr]->SetFillColor(color);
            h_mass_temp[pr]->SetLineColor(color);
            h_mass_test[pr]->SetFillColor(color);
            h_mass_test[pr]->SetLineColor(color);
            h_pT_e[pr]->SetFillColor(color);
            h_pT_e[pr]->SetLineColor(color);
            h_pT_mu[pr]->SetFillColor(color);
            h_pT_mu[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
            h_mass_temp[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_temp[pr]->SetMarkerColor(kBlack);
            h_mass_temp[pr]->SetLineColor(kBlack);
            h_mass_test[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test[pr]->SetMarkerColor(kBlack);
            h_mass_test[pr]->SetLineColor(kBlack);
            h_pT_e[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_e[pr]->SetMarkerColor(kBlack);
            h_pT_e[pr]->SetLineColor(kBlack);
            h_pT_mu[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_mu[pr]->SetMarkerColor(kBlack);
            h_pT_mu[pr]->SetLineColor(kBlack);
        }

        // Adding hists to THStacks
        if (pr < _DoubleEG_B)
        {
            s_mass_woQCD->Add(h_mass[pr]);
            s_mass_woQCD_SS->Add(h_mass_SS[pr]);
            s_mass_woQCD_temp->Add(h_mass_temp[pr]);
            s_mass_woQCD_test->Add(h_mass_test[pr]);
            s_pT_e->Add(h_pT_e[pr]);
            s_pT_mu->Add(h_pT_mu[pr]);
        }

        // Adding up for convenience
        if (pr == _DY_10to50)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
            h_mass_temp[_DY_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DY_Full")));
            h_mass_temp[_DY_Full]->SetDirectory(0);
            h_mass_test[_DY_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DY_Full")));
            h_mass_test[_DY_Full]->SetDirectory(0);
            h_pT_e[_DY_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_DY_Full")));
            h_pT_e[_DY_Full]->SetDirectory(0);
            h_pT_mu[_DY_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_DY_Full")));
            h_pT_mu[_DY_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_DY_Normal)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DY_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DY_Full]->Add(h_mass_test[pr]);
            h_pT_e[_DY_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_DY_Full]->Add(h_pT_e[pr]);
        }
        else if (pr == _ttbar)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
            h_mass_temp[_ttbar_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_ttbar_Full")));
            h_mass_temp[_ttbar_Full]->SetDirectory(0);
            h_mass_test[_ttbar_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_ttbar_Full")));
            h_mass_test[_ttbar_Full]->SetDirectory(0);
            h_pT_e[_ttbar_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_ttbar_Full")));
            h_pT_e[_ttbar_Full]->SetDirectory(0);
            h_pT_mu[_ttbar_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_ttbar_Full")));
            h_pT_mu[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_ttbar_Normal)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_ttbar_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_ttbar_Full]->Add(h_mass_test[pr]);
            h_pT_e[_ttbar_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_ttbar_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _tW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
            h_mass_temp[_VVnST] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_VVnST")));
            h_mass_temp[_VVnST]->SetDirectory(0);
            h_mass_test[_VVnST] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_VVnST")));
            h_mass_test[_VVnST]->SetDirectory(0);
            h_pT_e[_VVnST] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_VVnST")));
            h_pT_e[_VVnST]->SetDirectory(0);
            h_pT_mu[_VVnST] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_VVnST")));
            h_pT_mu[_VVnST]->SetDirectory(0);
        }
        else if (pr < _EndOf_VVnST_Normal)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
            h_mass_temp[_VVnST]->Add(h_mass_temp[pr]);
            h_mass_test[_VVnST]->Add(h_mass_test[pr]);
            h_pT_e[_VVnST]->Add(h_pT_e[pr]);
            h_pT_mu[_VVnST]->Add(h_pT_mu[pr]);
        }
        else if (pr == _WJets)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
            h_mass_temp[_WJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_WJets_Full")));
            h_mass_temp[_WJets_Full]->SetDirectory(0);
            h_mass_test[_WJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_WJets_Full")));
            h_mass_test[_WJets_Full]->SetDirectory(0);
            h_pT_e[_WJets_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_WJets_Full")));
            h_pT_e[_WJets_Full]->SetDirectory(0);
            h_pT_mu[_WJets_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_WJets_Full")));
            h_pT_mu[_WJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_WJets_Normal)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_WJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_WJets_Full]->Add(h_mass_test[pr]);
            h_pT_e[_WJets_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_WJets_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _GJets_20to100)
        {
            h_mass[_GJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_GJets_Full")));
            h_mass[_GJets_Full]->SetDirectory(0);
            h_mass_SS[_GJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_GJets_Full")));
            h_mass_SS[_GJets_Full]->SetDirectory(0);
            h_mass_temp[_GJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_GJets_Full")));
            h_mass_temp[_GJets_Full]->SetDirectory(0);
            h_mass_test[_GJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_GJets_Full")));
            h_mass_test[_GJets_Full]->SetDirectory(0);
            h_pT_e[_GJets_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_GJets_Full")));
            h_pT_e[_GJets_Full]->SetDirectory(0);
            h_pT_mu[_GJets_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_GJets_Full")));
            h_pT_mu[_GJets_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_GJets_Normal)
        {
            h_mass[_GJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_GJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_GJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_GJets_Full]->Add(h_mass_test[pr]);
            h_pT_e[_GJets_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_GJets_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DoubleEG_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
            h_mass_SS[_SingleMuon_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DoubleEG_Full")));
            h_mass_SS[_SingleMuon_Full]->SetDirectory(0);
            h_mass_temp[_SingleMuon_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DoubleEG_Full")));
            h_mass_temp[_SingleMuon_Full]->SetDirectory(0);
            h_mass_test[_SingleMuon_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DoubleEG_Full")));
            h_mass_test[_SingleMuon_Full]->SetDirectory(0);
            h_pT_e[_SingleMuon_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_DoubleEG_Full")));
            h_pT_e[_SingleMuon_Full]->SetDirectory(0);
            h_pT_mu[_SingleMuon_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_DoubleEG_Full")));
            h_pT_mu[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr < _EndOf_SingleMuon_Normal)
        {
            h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
            h_mass_SS[_SingleMuon_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_SingleMuon_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_SingleMuon_Full]->Add(h_mass_test[pr]);
            h_pT_e[_SingleMuon_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_SingleMuon_Full]->Add(h_pT_mu[pr]);
        }

        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDEMEnriched_Normal; // next -- GJets_20to100
        if (pr == _GJets_2000to5000) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuonB

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
//    h_QCD_est->Add(h_mass[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est->GetSize()-1; i_bin++)
    {
        h_QCD_est->SetBinError(i_bin, sqrt(h_QCD_est->GetBinContent(i_bin)));
        if (h_QCD_est->GetBinContent(i_bin) == 0)
        {
            h_QCD_est->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est->SetFillColor(kRed + 3);
    h_QCD_est->SetLineColor(kBlack);
    h_QCD_est->SetMarkerStyle(0);
    int_qcd = h_QCD_est->IntegralAndError(1, h_QCD_est->GetSize()-2, err_qcd);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "QCD est events: " << int_qcd << "+-" << err_qcd << endl;

    // Same-sign template
    h_QCD_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_QCD_est_SS")));
    h_QCD_est_SS->SetTitle("");
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
    h_QCD_est_SS->Add(h_mass_SS[_VVnST], -1);
    h_QCD_est_SS->Add(h_mass_SS[_WJets_Full], -1);
//    h_QCD_est_SS->Add(h_mass_SS[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est_SS);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est_SS->GetSize()-1; i_bin++)
    {
        h_QCD_est_SS->SetBinError(i_bin, sqrt(h_QCD_est_SS->GetBinContent(i_bin)));
        if (h_QCD_est_SS->GetBinContent(i_bin) == 0)
        {
            h_QCD_est_SS->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est_SS->SetFillColor(kRed + 3);
    h_QCD_est_SS->SetLineColor(kBlack);
    h_QCD_est_SS->SetMarkerStyle(0);

    // Possible template for ee
    h_QCD_est_temp = ((TH1D*)(h_mass_temp[_SingleMuon_Full]->Clone("h_QCD_est_template")));
    h_QCD_est_temp->SetTitle("");
    h_QCD_est_temp->SetDirectory(0);
    h_QCD_est_temp->Add(h_mass_temp[_DY_Full], -1);
    h_QCD_est_temp->Add(h_mass_temp[_ttbar_Full], -1);
    h_QCD_est_temp->Add(h_mass_temp[_VVnST], -1);
    h_QCD_est_temp->Add(h_mass_temp[_WJets_Full], -1);
//    h_QCD_est_temp->Add(h_mass_temp[_GJets_Full], -1);
    removeNegativeBins(h_QCD_est_temp);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_QCD_est_temp->GetSize()-1; i_bin++)
    {
        h_QCD_est_temp->SetBinError(i_bin, sqrt(h_QCD_est_temp->GetBinContent(i_bin)));
        if (h_QCD_est_temp->GetBinContent(i_bin) == 0)
        {
            h_QCD_est_temp->SetBinError(i_bin, 1);
        }
    }
    h_QCD_est_temp->SetFillColor(kRed + 3);
    h_QCD_est_temp->SetLineColor(kBlack);
    h_QCD_est_temp->SetMarkerStyle(0);


    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass_woQCD = new myRatioPlot_t("c_mass_woQCD", s_mass_woQCD, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD_SS = new myRatioPlot_t("c_mass_woQCD_SS", s_mass_woQCD_SS, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD_temp = new myRatioPlot_t("c_mass_woQCD_temp", s_mass_woQCD_temp, h_mass_temp[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woQCD_test = new myRatioPlot_t("c_mass_woQCD_test", s_mass_woQCD_test, h_mass_test[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_e = new myRatioPlot_t("c_pT_e", s_pT_e, h_pT_e[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_mu = new myRatioPlot_t("c_pT_mu", s_pT_mu, h_pT_mu[_SingleMuon_Full]);

    RP_mass_woQCD->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_SS->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} (same-sign) [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_temp->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woQCD_test->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_pT_e->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}) [GeV/c]", 0, 1000);
    RP_pT_mu->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{#mu}) [GeV/c]", 0, 1000);

    TLegend *legend = new TLegend(0.8, 0.52, 0.95, 0.95);
    legend->AddEntry(h_mass[_SingleMuon_B], "Data", "pl");
    legend->AddEntry(h_mass[_DY_10to50], "DY", "f");
    legend->AddEntry(h_mass[_ttbar], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_mass[_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_mass[_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_mass[_WW], "#font[12]{#scale[1.1]{WW}}", "f");
    legend->AddEntry(h_mass[_WZ], "#font[12]{#scale[1.1]{WZ}}", "f");
    legend->AddEntry(h_mass[_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}}", "f");
    legend->AddEntry(h_mass[_WJets], "#font[12]{#scale[1.1]{W}}+Jets", "f");
    legend->AddEntry(h_mass[_GJets_Full], "#gamma+Jets", "f");

    RP_mass_woQCD->ImportLegend(legend);
    RP_mass_woQCD_SS->ImportLegend(legend);
    RP_mass_woQCD_temp->ImportLegend(legend);
    RP_mass_woQCD_test->ImportLegend(legend);
    RP_pT_e->ImportLegend(legend);
    RP_pT_mu->ImportLegend(legend);

    RP_mass_woQCD->Draw(1e-2, 1e4, 1);
    RP_mass_woQCD_SS->Draw(1e-2, 1e4, 1);
    RP_mass_woQCD_temp->Draw(1e-2, 1e4, 1);
    RP_mass_woQCD_test->Draw(1e-1, 1e6, 1);
    RP_pT_e->Draw(1e-1, 1e6, 0);
    RP_pT_mu->Draw(1e-1, 1e6, 0);

    // Drawing estimated QCD
    TLegend * l_QCD_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_QCD_est->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}} (FR)", "f");
    TCanvas * c_QCD_est = new TCanvas("c_QCD_est", "QCD est", 750, 850);
    c_QCD_est->SetTopMargin(0.05);
    c_QCD_est->SetRightMargin(0.05);
    c_QCD_est->SetBottomMargin(0.15);
    c_QCD_est->SetLeftMargin(0.15);
    h_QCD_est->Draw("BAR");
    h_QCD_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]");
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
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_QCD_est->Write();
    h_QCD_est_SS->Write();
    h_QCD_est_temp->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu_TEMPLATE_DOWN.root", "READ");
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

    // Drawing estimated QCD same-sign template
    TCanvas * c_QCD_est_SS = new TCanvas("c_QCD_est_SS", "QCD est (same-sign template)", 750, 850);
    c_QCD_est_SS->SetTopMargin(0.05);
    c_QCD_est_SS->SetRightMargin(0.05);
    c_QCD_est_SS->SetBottomMargin(0.15);
    c_QCD_est_SS->SetLeftMargin(0.15);
    h_QCD_est_SS->Draw("BAR");
    h_QCD_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} (same-sign) [GeV/c^{2}]");
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
    l_QCD_est->Draw();
    c_QCD_est_SS->SetLogx();
    c_QCD_est_SS->SetGridx();
    c_QCD_est_SS->SetGridy();
    c_QCD_est_SS->Update();    

    // Drawing estimated QCD template for ee
    TCanvas * c_QCD_est_temp = new TCanvas("c_QCD_est_temp", "QCD est (template for ee)", 750, 850);
    c_QCD_est_temp->SetTopMargin(0.05);
    c_QCD_est_temp->SetRightMargin(0.05);
    c_QCD_est_temp->SetBottomMargin(0.15);
    c_QCD_est_temp->SetLeftMargin(0.15);
    h_QCD_est_temp->Draw("BAR");
    h_QCD_est_temp->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]");
    h_QCD_est_temp->GetXaxis()->SetTitleSize(0.062);
    h_QCD_est_temp->GetXaxis()->SetTitleOffset(0.9);
    h_QCD_est_temp->GetXaxis()->SetLabelSize(0.048);
    h_QCD_est_temp->GetXaxis()->SetMoreLogLabels();
    h_QCD_est_temp->GetXaxis()->SetNoExponent();
    h_QCD_est_temp->GetYaxis()->SetTitle("Number of events");
    h_QCD_est_temp->GetYaxis()->SetTitleSize(0.05);
    h_QCD_est_temp->GetYaxis()->SetTitleOffset(1.4);
    h_QCD_est_temp->GetYaxis()->SetLabelSize(0.043);
    h_QCD_est_temp->GetYaxis()->SetMoreLogLabels();
    h_QCD_est_temp->GetYaxis()->SetNoExponent();
    l_QCD_est->Draw();
    c_QCD_est_temp->SetLogx();
    c_QCD_est_temp->SetGridx();
    c_QCD_est_temp->SetGridy();
    c_QCD_est_temp->Update();

    // Same-sign vs opposite-sign plots
    TH1D *h_mass_QCD = ((TH1D*)(h_QCD_est->Clone("h_mass_QCD")));
    TH1D *h_mass_SS_QCD = ((TH1D*)(h_QCD_est_SS->Clone("h_mass_SS_QCD")));
    h_mass_QCD->SetMarkerStyle(kFullDotLarge);
    h_mass_QCD->SetMarkerColor(kBlack);
    h_mass_QCD->SetLineColor(kBlack);
    h_mass_SS_QCD->SetMarkerStyle(kFullSquare);
    h_mass_SS_QCD->SetMarkerColor(kRed);
    h_mass_SS_QCD->SetLineColor(kRed);

    myRatioPlot_t *RP_mass_QCD_SSvsOS = new myRatioPlot_t("c_mass_QCD_SSvsOS", h_mass_QCD, h_mass_SS_QCD);
    RP_mass_QCD_SSvsOS->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000, "OS/SS");
    TLegend *l_SSvsOS = new TLegend(0.5, 0.8, 0.95, 0.95);
    l_SSvsOS->AddEntry(h_mass_QCD, "Opposite-sign QCD estimation", "lp");
    l_SSvsOS->AddEntry(h_mass_SS_QCD, "Same-sign QCD estimation", "lp");
    RP_mass_QCD_SSvsOS->ImportLegend(l_SSvsOS);
    RP_mass_QCD_SSvsOS->Draw(5, 1e3, 1, "E");

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_EMu.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_EMu.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu.root COULD NOT BE CLOSED!\n" << endl;

} // End of EMu_QCDest_HistDrawer()


void EMu_WJETest_HistDrawer(Int_t remNegBins, Int_t systErr)
{
    FileMgr Mgr;

    TString Dir = "/media/sf_DATA/FR/EMu/";
    TFile *f = new TFile(Dir+"WJETest_EMu.root", "READ");

    TH1D *h_mass[_EndOf_Data_Special], *h_mass_SS[_EndOf_Data_Special], *h_mass_temp[_EndOf_Data_Special],
         *h_mass_test[_EndOf_Data_Special], *h_mass_test_ef[_EndOf_Data_Special], *h_mass_test_mf[_EndOf_Data_Special],
         *h_pT_e[_EndOf_Data_Special], *h_pT_mu[_EndOf_Data_Special];
    TH1D *h_WJET_est, *h_WJET_est_SS, *h_WJET_est_temp, *h_WJET_est_fit;
    THStack * s_mass_wWJET = new THStack("s_mass_wWJET", "");
    THStack * s_mass_woWJET = new THStack("s_mass_woWJET", "");
    THStack * s_mass_wWJET_SS = new THStack("s_mass_wWJET_SS", "");
    THStack * s_mass_woWJET_SS = new THStack("s_mass_woWJET_SS", "");
    THStack * s_mass_wWJET_temp = new THStack("s_mass_wWJET_template", "");
    THStack * s_mass_woWJET_temp = new THStack("s_mass_woWJET_template", "");
    THStack * s_mass_test = new THStack("s_mass_test", "");
    THStack * s_mass_test_ef = new THStack("s_mass_test_ef", "");
    THStack * s_mass_test_mf = new THStack("s_mass_test_mf", "");
    THStack * s_pT_e = new THStack("s_pT_e", "");
    THStack * s_pT_mu = new THStack("s_pT_mu", "");
    Color_t color = kBlack;

    // Getting data-driven QCD to subtract from data
    TFile *f_QCD = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstQCD_EMu.root");
    TH1D *h_QCD_est, *h_QCD_est_SS, *h_QCD_est_temp;

    f_QCD->GetObject("h_QCD_est", h_QCD_est);
    h_QCD_est->SetDirectory(0);
    h_QCD_est->SetFillColor(kRed+3);
    h_QCD_est->SetLineColor(kRed+3);
//    h_QCD_est->Scale(2);
    s_mass_wWJET->Add(h_QCD_est);
    s_mass_woWJET->Add(h_QCD_est);

    f_QCD->GetObject("h_QCD_est_SS", h_QCD_est_SS);
    h_QCD_est_SS->SetDirectory(0);
    h_QCD_est_SS->SetFillColor(kRed+3);
    h_QCD_est_SS->SetLineColor(kRed+3);
//    h_QCD_est_SS->Scale(2);
    s_mass_wWJET_SS->Add(h_QCD_est_SS);
    s_mass_woWJET_SS->Add(h_QCD_est_SS);

    f_QCD->GetObject("h_QCD_est_template", h_QCD_est_temp);
    h_QCD_est_temp->SetDirectory(0);
    h_QCD_est_temp->SetFillColor(kRed+3);
    h_QCD_est_temp->SetLineColor(kRed+3);
//    h_QCD_est_temp->Scale(2);
    s_mass_wWJET_temp->Add(h_QCD_est_temp);
    s_mass_woWJET_temp->Add(h_QCD_est_temp);

    for (Process_t pr=_SingleMuon_H; pr>=_DY_10to50; pr=previous(pr))
    {
        Mgr.SetProc(pr);

        f->GetObject("h_mass_"+Mgr.Procname[pr], h_mass[pr]);
        f->GetObject("h_mass_SS_"+Mgr.Procname[pr], h_mass_SS[pr]);
        f->GetObject("h_mass_template_"+Mgr.Procname[pr], h_mass_temp[pr]);
        f->GetObject("h_mass_test_"+Mgr.Procname[pr], h_mass_test[pr]);
        f->GetObject("h_mass_test_elefail_"+Mgr.Procname[pr], h_mass_test_ef[pr]);
        f->GetObject("h_mass_test_mufail_"+Mgr.Procname[pr], h_mass_test_mf[pr]);
        f->GetObject("h_pT_e_"+Mgr.Procname[pr], h_pT_e[pr]);
        f->GetObject("h_pT_mu_"+Mgr.Procname[pr], h_pT_mu[pr]);
        h_mass[pr]->SetDirectory(0);
        h_mass_SS[pr]->SetDirectory(0);
        h_mass_temp[pr]->SetDirectory(0);
        h_mass_test[pr]->SetDirectory(0);
        h_mass_test_ef[pr]->SetDirectory(0);
        h_mass_test_mf[pr]->SetDirectory(0);
        h_pT_e[pr]->SetDirectory(0);
        h_pT_mu[pr]->SetDirectory(0);

        // 2 -- use negative bin removal, 1 -- do not use
        if (remNegBins > 1)
        {
            removeNegativeBins(h_mass[pr]);
            removeNegativeBins(h_mass_SS[pr]);
            removeNegativeBins(h_mass_temp[pr]);
            removeNegativeBins(h_mass_test[pr]);
            removeNegativeBins(h_mass_test_ef[pr]);
            removeNegativeBins(h_mass_test_mf[pr]);
            removeNegativeBins(h_pT_e[pr]);
            removeNegativeBins(h_pT_mu[pr]);
        }

        if (pr < _EndOf_DY_Normal) color = kOrange - 5;
        else if (pr < _EndOf_ttbar_Normal) color = kCyan + 2;
        else if (pr == _tW) color = kGreen + 2;
        else if (pr == _tbarW) color = kGreen - 2;
        else if (pr == _WW) color = kMagenta - 5;
        else if (pr == _WZ) color = kMagenta - 2;
        else if (pr == _ZZ) color = kMagenta - 6;
        else if (pr < _EndOf_WJets_Normal) color = kRed - 2;
        else if (pr < _EndOf_GJets_Normal) color = kYellow + 3;

        if (pr < _DoubleEG_B) // MC coloring
        {
            h_mass[pr]->SetFillColor(color);
            h_mass[pr]->SetLineColor(color);
            h_mass_SS[pr]->SetFillColor(color);
            h_mass_SS[pr]->SetLineColor(color);
            h_mass_temp[pr]->SetFillColor(color);
            h_mass_temp[pr]->SetLineColor(color);
            h_mass_test[pr]->SetFillColor(color);
            h_mass_test[pr]->SetLineColor(color);
            h_mass_test_ef[pr]->SetFillColor(color);
            h_mass_test_ef[pr]->SetLineColor(color);
            h_mass_test_mf[pr]->SetFillColor(color);
            h_mass_test_mf[pr]->SetLineColor(color);
            h_pT_e[pr]->SetFillColor(color);
            h_pT_e[pr]->SetLineColor(color);
            h_pT_mu[pr]->SetFillColor(color);
            h_pT_mu[pr]->SetLineColor(color);
        }
        else // Data coloring
        {
            h_mass[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass[pr]->SetMarkerColor(kBlack);
            h_mass[pr]->SetLineColor(kBlack);
            h_mass_SS[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_SS[pr]->SetMarkerColor(kBlack);
            h_mass_SS[pr]->SetLineColor(kBlack);
            h_mass_temp[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_temp[pr]->SetMarkerColor(kBlack);
            h_mass_temp[pr]->SetLineColor(kBlack);
            h_mass_test[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test[pr]->SetMarkerColor(kBlack);
            h_mass_test[pr]->SetLineColor(kBlack);
            h_mass_test_ef[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test_ef[pr]->SetMarkerColor(kBlack);
            h_mass_test_ef[pr]->SetLineColor(kBlack);
            h_mass_test_mf[pr]->SetMarkerStyle(kFullDotLarge);
            h_mass_test_mf[pr]->SetMarkerColor(kBlack);
            h_mass_test_mf[pr]->SetLineColor(kBlack);
            h_pT_e[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_e[pr]->SetMarkerColor(kBlack);
            h_pT_e[pr]->SetLineColor(kBlack);
            h_pT_mu[pr]->SetMarkerStyle(kFullDotLarge);
            h_pT_mu[pr]->SetMarkerColor(kBlack);
            h_pT_mu[pr]->SetLineColor(kBlack);
        }

        // Filling stack histograms
        if (pr < _EndOf_WJets_Normal || (pr >= _GJets_20to100 && pr <= _GJets_2000to5000))
        {
            s_mass_wWJET->Add(h_mass[pr]);
            s_mass_wWJET_SS->Add(h_mass_SS[pr]);
            s_mass_wWJET_temp->Add(h_mass_temp[pr]);
            if (pr < _WJets || (pr >= _GJets_20to100 && pr <= _GJets_2000to5000))
            {
                s_mass_woWJET->Add(h_mass[pr]);
                s_mass_woWJET_SS->Add(h_mass_SS[pr]);
                s_mass_woWJET_temp->Add(h_mass_temp[pr]);
            }
            if (pr < _SingleMuon_B)
            {
                s_mass_test->Add(h_mass_test[pr]);
                s_mass_test_ef->Add(h_mass_test_ef[pr]);
                s_mass_test_mf->Add(h_mass_test_mf[pr]);
                s_pT_e->Add(h_pT_e[pr]);
                s_pT_mu->Add(h_pT_mu[pr]);
            }
        }

        // Adding up for convenience
        if (pr == _SingleMuon_H)
        {
            h_mass[_SingleMuon_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_data_Full")));
            h_mass[_SingleMuon_Full]->SetDirectory(0);
            h_mass_SS[_SingleMuon_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_data_Full")));
            h_mass_SS[_SingleMuon_Full]->SetDirectory(0);
            h_mass_temp[_SingleMuon_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_data_Full")));
            h_mass_temp[_SingleMuon_Full]->SetDirectory(0);
            h_mass_test[_SingleMuon_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_data_Full")));
            h_mass_test[_SingleMuon_Full]->SetDirectory(0);
            h_mass_test_ef[_SingleMuon_Full] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_data_Full")));
            h_mass_test_ef[_SingleMuon_Full]->SetDirectory(0);
            h_mass_test_mf[_SingleMuon_Full] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_data_Full")));
            h_mass_test_mf[_SingleMuon_Full]->SetDirectory(0);
            h_pT_e[_SingleMuon_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_data_Full")));
            h_pT_e[_SingleMuon_Full]->SetDirectory(0);
            h_pT_mu[_SingleMuon_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_data_Full")));
            h_pT_mu[_SingleMuon_Full]->SetDirectory(0);
        }
        else if (pr >= _SingleMuon_B)
        {
            h_mass[_SingleMuon_Full]->Add(h_mass[pr]);
            h_mass_SS[_SingleMuon_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_SingleMuon_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_SingleMuon_Full]->Add(h_mass_test[pr]);
            h_mass_test_ef[_SingleMuon_Full]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_SingleMuon_Full]->Add(h_mass_test_mf[pr]);
            h_pT_e[_SingleMuon_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_SingleMuon_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _GJets_2000to5000)
        {
            h_mass[_GJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_GJets_Full")));
            h_mass[_GJets_Full]->SetDirectory(0);
            h_mass_SS[_GJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_GJets_Full")));
            h_mass_SS[_GJets_Full]->SetDirectory(0);
            h_mass_temp[_GJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_GJets_Full")));
            h_mass_temp[_GJets_Full]->SetDirectory(0);
            h_mass_test[_GJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_GJets_Full")));
            h_mass_test[_GJets_Full]->SetDirectory(0);
            h_mass_test_ef[_GJets_Full] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_GJets_Full")));
            h_mass_test_ef[_GJets_Full]->SetDirectory(0);
            h_mass_test_mf[_GJets_Full] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_GJets_Full")));
            h_mass_test_mf[_GJets_Full]->SetDirectory(0);
            h_pT_e[_GJets_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_GJets_Full")));
            h_pT_e[_GJets_Full]->SetDirectory(0);
            h_pT_mu[_GJets_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_GJets_Full")));
            h_pT_mu[_GJets_Full]->SetDirectory(0);
        }
        else if (pr >= _GJets_20to100)
        {
            h_mass[_GJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_GJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_GJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_GJets_Full]->Add(h_mass_test[pr]);
            h_mass_test_ef[_GJets_Full]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_GJets_Full]->Add(h_mass_test_mf[pr]);
            h_pT_e[_GJets_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_GJets_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _WJets_ext2v5)
        {
            h_mass[_WJets_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_WJets_Full")));
            h_mass[_WJets_Full]->SetDirectory(0);
            h_mass_SS[_WJets_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_WJets_Full")));
            h_mass_SS[_WJets_Full]->SetDirectory(0);
            h_mass_temp[_WJets_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_WJets_Full")));
            h_mass_temp[_WJets_Full]->SetDirectory(0);
            h_mass_test[_WJets_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_WJets_Full")));
            h_mass_test[_WJets_Full]->SetDirectory(0);
            h_mass_test_ef[_WJets_Full] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_WJets_Full")));
            h_mass_test_ef[_WJets_Full]->SetDirectory(0);
            h_mass_test_mf[_WJets_Full] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_WJets_Full")));
            h_mass_test_mf[_WJets_Full]->SetDirectory(0);
            h_pT_e[_WJets_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_WJets_Full")));
            h_pT_e[_WJets_Full]->SetDirectory(0);
            h_pT_mu[_WJets_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_WJets_Full")));
            h_pT_mu[_WJets_Full]->SetDirectory(0);
        }
        else if (pr >= _WJets)
        {
            h_mass[_WJets_Full]->Add(h_mass[pr]);
            h_mass_SS[_WJets_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_WJets_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_WJets_Full]->Add(h_mass_test[pr]);
            h_mass_test_ef[_WJets_Full]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_WJets_Full]->Add(h_mass_test_mf[pr]);
            h_pT_e[_WJets_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_WJets_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _WW)
        {
            h_mass[_VVnST] = ((TH1D*)(h_mass[pr]->Clone("h_mass_VVnST")));
            h_mass[_VVnST]->SetDirectory(0);
            h_mass_SS[_VVnST] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_VVnST")));
            h_mass_SS[_VVnST]->SetDirectory(0);
            h_mass_temp[_VVnST] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_VVnST")));
            h_mass_temp[_VVnST]->SetDirectory(0);
            h_mass_test[_VVnST] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_VVnST")));
            h_mass_test[_VVnST]->SetDirectory(0);
            h_mass_test_ef[_VVnST] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_VVnST")));
            h_mass_test_ef[_VVnST]->SetDirectory(0);
            h_mass_test_mf[_VVnST] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_VVnST")));
            h_mass_test_mf[_VVnST]->SetDirectory(0);
            h_pT_e[_VVnST] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_VVnST")));
            h_pT_e[_VVnST]->SetDirectory(0);
            h_pT_mu[_VVnST] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_VVnST")));
            h_pT_mu[_VVnST]->SetDirectory(0);
        }
        else if (pr >= _tW)
        {
            h_mass[_VVnST]->Add(h_mass[pr]);
            h_mass_SS[_VVnST]->Add(h_mass_SS[pr]);
            h_mass_temp[_VVnST]->Add(h_mass_temp[pr]);
            h_mass_test[_VVnST]->Add(h_mass_test[pr]);
            h_mass_test_ef[_VVnST]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_VVnST]->Add(h_mass_test_mf[pr]);
            h_pT_e[_VVnST]->Add(h_pT_e[pr]);
            h_pT_mu[_VVnST]->Add(h_pT_mu[pr]);
        }
        else if (pr == _ttbar_1000toInf)
        {
            h_mass[_ttbar_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_ttbar_Full")));
            h_mass[_ttbar_Full]->SetDirectory(0);
            h_mass_SS[_ttbar_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_ttbar_Full")));
            h_mass_SS[_ttbar_Full]->SetDirectory(0);
            h_mass_temp[_ttbar_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_ttbar_Full")));
            h_mass_temp[_ttbar_Full]->SetDirectory(0);
            h_mass_test[_ttbar_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_ttbar_Full")));
            h_mass_test[_ttbar_Full]->SetDirectory(0);
            h_mass_test_ef[_ttbar_Full] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_ttbar_Full")));
            h_mass_test_ef[_ttbar_Full]->SetDirectory(0);
            h_mass_test_mf[_ttbar_Full] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_ttbar_Full")));
            h_mass_test_mf[_ttbar_Full]->SetDirectory(0);
            h_pT_e[_ttbar_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_ttbar_Full")));
            h_pT_e[_ttbar_Full]->SetDirectory(0);
            h_pT_mu[_ttbar_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_ttbar_Full")));
            h_pT_mu[_ttbar_Full]->SetDirectory(0);
        }
        else if (pr >= _ttbar)
        {
            h_mass[_ttbar_Full]->Add(h_mass[pr]);
            h_mass_SS[_ttbar_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_ttbar_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_ttbar_Full]->Add(h_mass_test[pr]);
            h_mass_test_ef[_ttbar_Full]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_ttbar_Full]->Add(h_mass_test_mf[pr]);
            h_pT_e[_ttbar_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_ttbar_Full]->Add(h_pT_mu[pr]);
        }
        else if (pr == _DY_2000to3000)
        {
            h_mass[_DY_Full] = ((TH1D*)(h_mass[pr]->Clone("h_mass_DY_Full")));
            h_mass[_DY_Full]->SetDirectory(0);
            h_mass_SS[_DY_Full] = ((TH1D*)(h_mass_SS[pr]->Clone("h_mass_SS_DY_Full")));
            h_mass_SS[_DY_Full]->SetDirectory(0);
            h_mass_temp[_DY_Full] = ((TH1D*)(h_mass_temp[pr]->Clone("h_mass_template_DY_Full")));
            h_mass_temp[_DY_Full]->SetDirectory(0);
            h_mass_test[_DY_Full] = ((TH1D*)(h_mass_test[pr]->Clone("h_mass_test_DY_Full")));
            h_mass_test[_DY_Full]->SetDirectory(0);
            h_mass_test_ef[_DY_Full] = ((TH1D*)(h_mass_test_ef[pr]->Clone("h_mass_test_elefail_DY_Full")));
            h_mass_test_ef[_DY_Full]->SetDirectory(0);
            h_mass_test_mf[_DY_Full] = ((TH1D*)(h_mass_test_mf[pr]->Clone("h_mass_test_mufail_DY_Full")));
            h_mass_test_mf[_DY_Full]->SetDirectory(0);
            h_pT_e[_DY_Full] = ((TH1D*)(h_pT_e[pr]->Clone("h_pT_e_DY_Full")));
            h_pT_e[_DY_Full]->SetDirectory(0);
            h_pT_mu[_DY_Full] = ((TH1D*)(h_pT_mu[pr]->Clone("h_pT_mu_DY_Full")));
            h_pT_mu[_DY_Full]->SetDirectory(0);
        }
        else if (pr >= _DY_10to50)
        {
            h_mass[_DY_Full]->Add(h_mass[pr]);
            h_mass_SS[_DY_Full]->Add(h_mass_SS[pr]);
            h_mass_temp[_DY_Full]->Add(h_mass_temp[pr]);
            h_mass_test[_DY_Full]->Add(h_mass_test[pr]);
            h_mass_test_ef[_DY_Full]->Add(h_mass_test_ef[pr]);
            h_mass_test_mf[_DY_Full]->Add(h_mass_test_mf[pr]);
            h_pT_e[_DY_Full]->Add(h_pT_e[pr]);
            h_pT_mu[_DY_Full]->Add(h_pT_mu[pr]);
        }

        if (pr == _SingleMuon_B) pr = _EndOf_GJets_Normal; // next -- GJets_2000to5000
        if (pr == _GJets_20to100) pr = _EndOf_WJets_Normal; // next -- WJets_ext2v5
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
//    h_WJET_est->Add(h_mass[_GJets_Full], -1);
    h_WJET_est->Add(h_QCD_est, -1);
    removeNegativeBins(h_WJET_est);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est->GetSize()-1; i_bin++)
    {
        h_WJET_est->SetBinError(i_bin, sqrt(h_WJET_est->GetBinContent(i_bin)));
        if (h_WJET_est->GetBinContent(i_bin) == 0)
        {
            h_WJET_est->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est->SetFillColor(kRed - 2);
    h_WJET_est->SetLineColor(kBlack);
    h_WJET_est->SetMarkerStyle(0);
    int_wjet = h_WJET_est->IntegralAndError(1, h_WJET_est->GetSize()-2, err_wjet);
    cout << "Data events: " << int_data << "+-" << err_data << endl;
    cout << "WJets est events (subraction): " << int_wjet << "+-" << err_wjet << endl;

    // Same-sign template
    h_WJET_est_SS = ((TH1D*)(h_mass_SS[_SingleMuon_Full]->Clone("h_WJET_est_SS")));
    h_WJET_est_SS->SetTitle("");
    h_WJET_est_SS->SetDirectory(0);
    h_WJET_est_SS->Add(h_mass_SS[_DY_Full], -1);
    h_WJET_est_SS->Add(h_mass_SS[_ttbar_Full], -1);
//    h_WJET_est_SS->Add(h_mass_SS[_WJets_Full], -1); // REMOVE AFTER TEST
    h_WJET_est_SS->Add(h_mass_SS[_VVnST], -1);
//    h_WJET_est_SS->Add(h_mass_SS[_GJets_Full], -1);
    h_WJET_est_SS->Add(h_QCD_est_SS, -1);
    removeNegativeBins(h_WJET_est_SS);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_SS->GetSize()-1; i_bin++)
    {
        h_WJET_est_SS->SetBinError(i_bin, sqrt(h_WJET_est_SS->GetBinContent(i_bin)));
        if (h_WJET_est_SS->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_SS->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est_SS->SetFillColor(kRed - 2);
    h_WJET_est_SS->SetLineColor(kBlack);
    h_WJET_est_SS->SetMarkerStyle(0);

    // Template for fit
    h_WJET_est_temp = ((TH1D*)(h_mass_temp[_SingleMuon_Full]->Clone("h_WJET_est_template")));
    h_WJET_est_temp->SetTitle("");
    h_WJET_est_temp->SetDirectory(0);
    h_WJET_est_temp->Add(h_mass_temp[_DY_Full], -1);
    h_WJET_est_temp->Add(h_mass_temp[_ttbar_Full], -1);
//    h_WJET_est_temp->Add(h_mass_temp[_WJets_Full], -1); // REMOVE AFTER TEST
    h_WJET_est_temp->Add(h_mass_temp[_VVnST], -1);
//    h_WJET_est_temp->Add(h_mass_temp[_GJets_Full], -1);
    h_WJET_est_temp->Add(h_QCD_est_temp, -1);
    removeNegativeBins(h_WJET_est_temp);
    // Setting errors
    for (Int_t i_bin=1; i_bin<h_WJET_est_temp->GetSize()-1; i_bin++)
    {
        h_WJET_est_temp->SetBinError(i_bin, sqrt(h_WJET_est_temp->GetBinContent(i_bin)));
        if (h_WJET_est_temp->GetBinContent(i_bin) == 0)
        {
            h_WJET_est_temp->SetBinError(i_bin, 1);
        }
    }
    h_WJET_est_temp->SetFillColor(kRed - 2);
    h_WJET_est_temp->SetLineColor(kBlack);
    h_WJET_est_temp->SetMarkerStyle(0);

    myRatioPlot_t *RP_mass_wWJET = new myRatioPlot_t("c_mass_wWJET", s_mass_wWJET, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woWJET = new myRatioPlot_t("c_mass_woWJET", s_mass_woWJET, h_mass[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_wWJET_SS = new myRatioPlot_t("c_mass_wWJET_SS", s_mass_wWJET_SS, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woWJET_SS = new myRatioPlot_t("c_mass_woWJET_SS", s_mass_woWJET_SS, h_mass_SS[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_wWJET_temp = new myRatioPlot_t("c_mass_wWJET_template", s_mass_wWJET_temp, h_mass_temp[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_woWJET_temp = new myRatioPlot_t("c_mass_woWJET_template", s_mass_woWJET_temp, h_mass_temp[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_test = new myRatioPlot_t("c_mass_test", s_mass_test, h_mass_test[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_test_ef = new myRatioPlot_t("c_mass_test_elefail", s_mass_test_ef, h_mass_test_ef[_SingleMuon_Full]);
    myRatioPlot_t *RP_mass_test_mf = new myRatioPlot_t("c_mass_test_mufail", s_mass_test_mf, h_mass_test_mf[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_e = new myRatioPlot_t("c_pT_e", s_pT_e, h_pT_e[_SingleMuon_Full]);
    myRatioPlot_t *RP_pT_mu = new myRatioPlot_t("c_pT_mu", s_pT_mu, h_pT_mu[_SingleMuon_Full]);

    RP_mass_wWJET->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_wWJET_SS->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET_SS->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_wWJET_temp->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_woWJET_temp->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_test->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_test_ef->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_mass_test_mf->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000);
    RP_pT_e->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{e}) [GeV/c]", 0, 1000);
    RP_pT_mu->SetPlots("p_{#lower[-0.2]{T}} (#font[12]{#mu}) [GeV/c]", 0, 1000);

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
    legend_wWJET->AddEntry(h_mass[_GJets_20to100], "#gamma+Jets", "f");
    legend_wWJET->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}}", "f");
    legend_woWJET->AddEntry(h_mass[_GJets_20to100], "#gamma+Jets", "f");
    legend_woWJET->AddEntry(h_QCD_est, "#font[12]{#scale[1.1]{QCD}}", "f");

    RP_mass_wWJET->ImportLegend(legend_wWJET);
    RP_mass_woWJET->ImportLegend(legend_woWJET);
    RP_mass_wWJET_SS->ImportLegend(legend_wWJET);
    RP_mass_woWJET_SS->ImportLegend(legend_woWJET);
    RP_mass_wWJET_temp->ImportLegend(legend_wWJET);
    RP_mass_woWJET_temp->ImportLegend(legend_woWJET);
    RP_mass_test->ImportLegend(legend_wWJET);
    RP_mass_test_ef->ImportLegend(legend_wWJET);
    RP_mass_test_mf->ImportLegend(legend_wWJET);
    RP_pT_e->ImportLegend(legend_wWJET);
    RP_pT_mu->ImportLegend(legend_wWJET);

    RP_mass_wWJET->Draw(1e-2, 1e5, 1);
    RP_mass_woWJET->Draw(1e-2, 1e5, 1);
    RP_mass_wWJET_SS->Draw(1e-2, 1e4, 1);
    RP_mass_woWJET_SS->Draw(1e-2, 1e4, 1);
    RP_mass_wWJET_temp->Draw(1e-2, 1e4, 1);
    RP_mass_woWJET_temp->Draw(1e-2, 1e4, 1);
    RP_mass_test->Draw(1e-2, 1e5, 1);
    RP_mass_test_ef->Draw(10, 1e6, 1);
    RP_mass_test_mf->Draw(1e-1, 1e5, 1);
    RP_pT_e->Draw(1e-1, 1e7, 0);
    RP_pT_mu->Draw(1e-1, 1e7, 0);

    TLegend * l_WJET_est = new TLegend(0.7, 0.88, 0.95, 0.95);
    l_WJET_est->AddEntry(h_WJET_est, "#font[12]{#scale[1.1]{W}}+Jets (FR)", "f");

    // Draw WJets from simple subtraction
    TCanvas * c_WJET_est = new TCanvas("c_WJET_est", "W+Jets est", 750, 850);
    c_WJET_est->SetTopMargin(0.05);
    c_WJET_est->SetRightMargin(0.05);
    c_WJET_est->SetBottomMargin(0.15);
    c_WJET_est->SetLeftMargin(0.17);
    h_WJET_est->Draw("BAR");
    h_WJET_est->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]");
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

    // Draw WJets from simple subtraction (same-sign template)
    TCanvas * c_WJET_est_SS = new TCanvas("c_WJET_est_SS", "W+Jets est (same-sign template)", 750, 850);
    c_WJET_est_SS->SetTopMargin(0.05);
    c_WJET_est_SS->SetRightMargin(0.05);
    c_WJET_est_SS->SetBottomMargin(0.15);
    c_WJET_est_SS->SetLeftMargin(0.17);
    h_WJET_est_SS->Draw("BAR");
    h_WJET_est_SS->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} (same-sign) [GeV/c^{2}]");
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
    l_WJET_est->Draw();
    c_WJET_est_SS->SetLogx();
    c_WJET_est_SS->SetGridx();
    c_WJET_est_SS->SetGridy();
    c_WJET_est_SS->Update();

    // Draw WJets from simple subtraction (template for ee)
    TCanvas * c_WJET_est_temp = new TCanvas("c_WJET_est_temp", "W+Jets est (template for ee)", 750, 850);
    c_WJET_est_temp->SetTopMargin(0.05);
    c_WJET_est_temp->SetRightMargin(0.05);
    c_WJET_est_temp->SetBottomMargin(0.15);
    c_WJET_est_temp->SetLeftMargin(0.17);
    h_WJET_est_temp->Draw("BAR");
    h_WJET_est_temp->GetXaxis()->SetTitle("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]");
    h_WJET_est_temp->GetXaxis()->SetTitleSize(0.062);
    h_WJET_est_temp->GetXaxis()->SetTitleOffset(0.9);
    h_WJET_est_temp->GetXaxis()->SetLabelSize(0.048);
    h_WJET_est_temp->GetXaxis()->SetMoreLogLabels();
    h_WJET_est_temp->GetXaxis()->SetNoExponent();
    h_WJET_est_temp->GetYaxis()->SetTitle("Number of events");
    h_WJET_est_temp->GetYaxis()->SetTitleSize(0.05);
    h_WJET_est_temp->GetYaxis()->SetTitleOffset(1.6);
    h_WJET_est_temp->GetYaxis()->SetLabelSize(0.043);
    h_WJET_est_temp->GetYaxis()->SetMoreLogLabels();
    h_WJET_est_temp->GetYaxis()->SetNoExponent();
    l_WJET_est->Draw();
    c_WJET_est_temp->SetLogx();
    c_WJET_est_temp->SetGridx();
    c_WJET_est_temp->SetGridy();
    c_WJET_est_temp->Update();

    // Starting writing
    TFile * f_out = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu.root", "RECREATE");
    if (f_out->IsOpen()) cout << "Writing output file..." << endl;
    else cout << "Error while writing output!" << endl;
    h_WJET_est->Write();
    h_WJET_est_SS->Write();
    h_WJET_est_temp->Write();

    if (systErr > 0)
    {   // Errors
        TFile *f_ratio = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu_RATIO.root", "READ");
        TFile *f_up = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu_TEMPLATE_UP.root", "READ");
        TFile *f_down = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu_TEMPLATE_DOWN.root", "READ");
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

    // Same-sign vs opposite-sign plots
    TH1D *h_mass_WJET = ((TH1D*)(h_WJET_est->Clone("h_mass_WJets")));
    TH1D *h_mass_SS_WJET = ((TH1D*)(h_WJET_est_SS->Clone("h_mass_SS_WJets")));
    h_mass_WJET->SetMarkerStyle(kFullDotLarge);
    h_mass_WJET->SetMarkerColor(kBlack);
    h_mass_WJET->SetLineColor(kBlack);
    h_mass_SS_WJET->SetMarkerStyle(kFullSquare);
    h_mass_SS_WJET->SetMarkerColor(kRed);
    h_mass_SS_WJET->SetLineColor(kRed);

    myRatioPlot_t *RP_mass_WJET_SSvsOS = new myRatioPlot_t("RP_mass_WJET_SSvsOS", h_mass_WJET, h_mass_SS_WJET);
    RP_mass_WJET_SSvsOS->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#font[12]{e#mu}}}} [GeV/c^{2}]", 15, 3000, "OS/SS");
    TLegend *l_SSvsOS = new TLegend(0.5, 0.8, 0.95, 0.95);
    l_SSvsOS->AddEntry(h_mass_WJET, "Opposite-sign W+Jets estimation", "lp");
    l_SSvsOS->AddEntry(h_mass_SS_WJET, "Same-sign W+Jets estimation", "lp");
    RP_mass_WJET_SSvsOS->ImportLegend(l_SSvsOS);
    RP_mass_WJET_SSvsOS->Draw(5, 1e3, 1, "");

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_EMu.root" << " has been closed successfully." << endl;
    else cout << "FILE " << Dir+"WJETest_EMu.root" << " COULD NOT BE CLOSED!\n" << endl;

    f_out->Close();
    if (!f_out->IsOpen()) cout << "File /media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu.root has been closed successfully.\n" << endl;
    else cout << "FILE /media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu.root COULD NOT BE CLOSED!\n" << endl;

} // End of EMu_WJETest_HistDrawer()


void E_PR_HistDrawer()
{
    TFile *f = new TFile("/media/sf_DATA/SelectedEE/Histos/DYefficiency.root", "READ");

    TH1D *h_mass[4],
         *h_pT_barrel_pass[4],
         *h_pT_endcap_pass[4],
         *h_pT_barrel_fail[4],
         *h_pT_endcap_fail[4],
         *h_MET_fail[4],
         *h_MT_fail[4];
    THStack *s_mass = new THStack("s_mass", "");
    THStack *s_pT_barrel_pass = new THStack("s_pT_barrel_pass", "");
    THStack *s_pT_endcap_pass = new THStack("s_pT_endcap_pass", "");
    THStack *s_pT_barrel_fail = new THStack("s_pT_barrel_fail", "");
    THStack *s_pT_endcap_fail = new THStack("s_pT_endcap_fail", "");
    THStack *s_MET_fail = new THStack("s_MET_fail", "");
    THStack *s_MT_fail = new THStack("s_MT_fail", "");
    Color_t color = kBlack;
    TString type[4] = {"data", "DY", "bkgr", "bkgf"};

    for (Int_t i=3; i>=0; i--)
    {
        f->GetObject("h_mass_"+type[i], h_mass[i]);
        f->GetObject("h_pT_barrel_pass_"+type[i], h_pT_barrel_pass[i]);
        f->GetObject("h_pT_endcap_pass_"+type[i], h_pT_endcap_pass[i]);
        f->GetObject("h_pT_barrel_fail_"+type[i], h_pT_barrel_fail[i]);
        f->GetObject("h_pT_endcap_fail_"+type[i], h_pT_endcap_fail[i]);
        f->GetObject("h_MET_fail_"+type[i], h_MET_fail[i]);
        f->GetObject("h_MT_fail_"+type[i], h_MT_fail[i]);
        h_mass[i]->SetDirectory(0);
        h_pT_barrel_pass[i]->SetDirectory(0);
        h_pT_endcap_pass[i]->SetDirectory(0);
        h_pT_barrel_fail[i]->SetDirectory(0);
        h_pT_endcap_fail[i]->SetDirectory(0);
        h_MET_fail[i]->SetDirectory(0);
        h_MT_fail[i]->SetDirectory(0);
        if (i == 0)
        {
            h_mass[i]->SetMarkerStyle(kFullDotLarge);
            h_mass[i]->SetMarkerColor(kBlack);
            h_mass[i]->SetLineColor(kBlack);
            h_pT_barrel_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_pass[i]->SetMarkerColor(kBlack);
            h_pT_barrel_pass[i]->SetLineColor(kBlack);
            h_pT_endcap_pass[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_pass[i]->SetMarkerColor(kBlack);
            h_pT_endcap_pass[i]->SetLineColor(kBlack);
            h_pT_barrel_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_barrel_fail[i]->SetMarkerColor(kBlack);
            h_pT_barrel_fail[i]->SetLineColor(kBlack);
            h_pT_endcap_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_pT_endcap_fail[i]->SetMarkerColor(kBlack);
            h_pT_endcap_fail[i]->SetLineColor(kBlack);
            h_MET_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_MET_fail[i]->SetMarkerColor(kBlack);
            h_MET_fail[i]->SetLineColor(kBlack);
            h_MT_fail[i]->SetMarkerStyle(kFullDotLarge);
            h_MT_fail[i]->SetMarkerColor(kBlack);
            h_MT_fail[i]->SetLineColor(kBlack);
        }
        else
        {
            if (i == 1) color = kOrange;
            else if (i == 2) color = kCyan + 2;
            else color = kRed + 3;
            h_mass[i]->SetFillColor(color);
            h_mass[i]->SetLineColor(color);
            h_pT_barrel_pass[i]->SetFillColor(color);
            h_pT_barrel_pass[i]->SetLineColor(color);
            h_pT_endcap_pass[i]->SetFillColor(color);
            h_pT_endcap_pass[i]->SetLineColor(color);
            h_pT_barrel_fail[i]->SetFillColor(color);
            h_pT_barrel_fail[i]->SetLineColor(color);
            h_pT_endcap_fail[i]->SetFillColor(color);
            h_pT_endcap_fail[i]->SetLineColor(color);
            h_MET_fail[i]->SetFillColor(color);
            h_MET_fail[i]->SetLineColor(color);
            h_MT_fail[i]->SetFillColor(color);
            h_MT_fail[i]->SetLineColor(color);


            s_mass->Add(h_mass[i]);
            s_pT_barrel_pass->Add(h_pT_barrel_pass[i]);
            s_pT_endcap_pass->Add(h_pT_endcap_pass[i]);
            s_pT_barrel_fail->Add(h_pT_barrel_fail[i]);
            s_pT_endcap_fail->Add(h_pT_endcap_fail[i]);
            s_MET_fail->Add(h_MET_fail[i]);
            s_MT_fail->Add(h_MT_fail[i]);
        }
    }


    // Creating and drawing ratio plots
    myRatioPlot_t *RP_mass = new myRatioPlot_t("c_mass", s_mass, h_mass[0]);
    myRatioPlot_t *RP_pT_barrel_pass = new myRatioPlot_t("c_pT_barrel_pass", s_pT_barrel_pass, h_pT_barrel_pass[0]);
    myRatioPlot_t *RP_pT_endcap_pass = new myRatioPlot_t("c_pT_endcap_pass", s_pT_endcap_pass, h_pT_endcap_pass[0]);
    myRatioPlot_t *RP_pT_barrel_fail = new myRatioPlot_t("c_pT_barrel_fail", s_pT_barrel_fail, h_pT_barrel_fail[0]);
    myRatioPlot_t *RP_pT_endcap_fail = new myRatioPlot_t("c_pT_endcap_fail", s_pT_endcap_fail, h_pT_endcap_fail[0]);
    myRatioPlot_t *RP_MET_fail = new myRatioPlot_t("c_MET_fail", s_MET_fail, h_MET_fail[0]);
    myRatioPlot_t *RP_MT_fail = new myRatioPlot_t("c_MT_fail", s_MT_fail, h_MT_fail[0]);

    RP_mass->SetPlots("m_{#lower[-0.2]{#font[12]{#scale[1.2]{ee}}}} [GeV/c^{2}]", 81, 101);
    RP_pT_barrel_pass->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 1000);
    RP_pT_endcap_pass->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 1000);
    RP_pT_barrel_fail->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 1000);
    RP_pT_endcap_fail->SetPlots("p_{#lower[-0.2]{T}} [GeV/c]", 0, 1000);
    RP_MET_fail->SetPlots("E_{#lower[-0.2]{T}}^{miss} [GeV]", 0, 100);
    RP_MT_fail->SetPlots("m_{#lower[-0.2]{T}} [GeV/c^{2}]", 0, 200);

    TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->AddEntry(h_mass[0], "Data", "pl");
    legend->AddEntry(h_mass[1], "DY", "f");
    legend->AddEntry(h_mass[2], "Bkg (real)", "f");
    legend->AddEntry(h_mass[3], "Bkg (fake)", "f");

    RP_mass->ImportLegend(legend);
    RP_pT_barrel_pass->ImportLegend(legend);
    RP_pT_endcap_pass->ImportLegend(legend);
    RP_pT_barrel_fail->ImportLegend(legend);
    RP_pT_endcap_fail->ImportLegend(legend);
    RP_MET_fail->ImportLegend(legend);
    RP_MT_fail->ImportLegend(legend);

    RP_mass->Draw(1e-1, 1e8, 1);
    RP_pT_barrel_pass->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_pass->Draw(1e-1, 1e7, 1);
    RP_pT_barrel_fail->Draw(1e-1, 1e7, 1);
    RP_pT_endcap_fail->Draw(1e-1, 1e7, 1);
    RP_MET_fail->Draw(1e-1, 1e7, 0);
    RP_MT_fail->Draw(1e-1, 1e7, 0);

    f->Close();
    if (!f->IsOpen()) cout << "File " << "/media/sf_DATA/SelectedEE/Histos/DYefficiency.root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " <<"/media/sf_DATA/SelectedEE/Histos/DYefficiency.root" << " COULD NOT BE CLOSED!\n" << endl;

} // End of E_DYefficiency()


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
