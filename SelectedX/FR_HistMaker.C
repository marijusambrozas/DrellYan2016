#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"

void E_FR_HistMaker (Bool_t DEBUG);
void Mu_FR_HistMaker (Bool_t DEBUG);

void E_QCD_HistMaker (Bool_t DEBUG);
void E_WJET_HistMaker (Bool_t DEBUG);
void Mu_QCD_HistMaker (Bool_t DEBUG, Int_t type);
void Mu_WJET_HistMaker (Bool_t DEBUG, Int_t type);
void EMu_QCD_HistMaker (Bool_t DEBUG);
void EMu_WJET_HistMaker (Bool_t DEBUG);

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

const Double_t etabins[51] = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
                              0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};

// -- RelPFiso bins for electron FR templates -- //
const Int_t binnum_ele_template = 5;
const Double_t bins_ele_template_barrel[6] = {0, 0.0695, 0.2, 0.4, 0.8, 3};
const Double_t bins_ele_template_endcap[6] = {0, 0.0821, 0.2, 0.4, 0.8, 3};
const Double_t bins_ele_HE_barrel[6] = {0, 0.0253, 0.026, 0.03, 0.04, 0.2};
const Double_t bins_ele_HE_endcap[6] = {0, 0.0878, 0.09, 0.092, 0.095, 0.5};
const Double_t bins_ele_sigma_barrel[6] = {0, 0.00998, 0.02, 0.03, 0.05, 0.1};
const Double_t bins_ele_sigma_endcap[6] = {0, 0.0298, 0.05, 0.07, 0.09, 0.2};
const Double_t bins_ele_dEta_barrel[6] = {0, 0.00311, 0.01, 0.03, 0.05, 0.1};
const Double_t bins_ele_dEta_endcap[6] = {0, 0.00609, 0.02, 0.05, 0.1, 0.5};
const Double_t bins_ele_dPhi_barrel[6] = {0, 0.103, 0.11, 0.12, 0.15, 0.2};
const Double_t bins_ele_dPhi_endcap[6] = {0, 0.045, 0.05, 0.07, 0.09, 0.2};
const Double_t bins_ele_EP_barrel[6] = {0, 0.134, 0.3, 0.5, 0.8, 1.2};
const Double_t bins_ele_EP_endcap[6] = {0, 0.13, 0.3, 0.5, 0.8, 1.2};

void FR_HistMaker (TString WhichX = "", Int_t type=1)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    Bool_t DEBUG = kFALSE;
    if (whichX.Contains("DEBUG"))
    {
        DEBUG = kTRUE;
        cout << "**** DEBUG MODE: Running with 100 events only. ****" << endl;
    }
    if (whichX.Contains("EMU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*****  Mu_QCD_HistMaker()  *****" << endl;
            EMu_QCD_HistMaker(DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  Mu_WJET_HistMaker()  *****" << endl;
            EMu_WJET_HistMaker(DEBUG);
        }
        else
        {
            cout << "Please specify by adding QCD or WJET" << endl;
            return;
        }
    }
    else if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*****  Mu_QCD_HistMaker(" << type << ")  *****" << endl;
            Mu_QCD_HistMaker(DEBUG, type);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  Mu_WJET_HistMaker(" << type << ")  *****" << endl;
            Mu_WJET_HistMaker(DEBUG, type);
        }
        else
        {
            cout << "\n*****  Mu_FR_HistMaker  *****" << endl;
            Mu_FR_HistMaker(DEBUG);
        }
    }
    else if (whichX.Contains("E"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*****  E_QCD_HistMaker  *****" << endl;
            E_QCD_HistMaker(DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  E_WJET_HistMaker  *****" << endl;
            E_WJET_HistMaker(DEBUG);
        }
        else
        {
            cout << "\n*****    E_FR_HistMaker    *****" << endl;
            E_FR_HistMaker(DEBUG);
        }
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FR_HistMaker.C+(\"whichX\", type)" << endl;

} // End of HistMaker()


/// ----------------------------- Electron Channel ------------------------------ ///
void E_FR_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    UInt_t n2MC=0, n2Data=0;

    for (Process_t pr=_DY_10to50; pr<_EndOf_SinglePhoton_Normal; pr=next(pr))
//    for (Process_t pr=_QCDEMEnriched_20to30; pr<_EndOf_QCDEMEnriched_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_E_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        /*///////////////////////////////////////
        if (pr == _WJets_ext2v5) // REMOVE LATER
        {
            f->Close();
            continue;
        }//////////////////////////////////////*/

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer("Photon_OR");

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_electron();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_endcap_deno->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBin_ele, analyzer->ptbin_ele /*250, 0, 500*/); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 50, etabins); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 50, etabins); h_eta_deno->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 50, etabins); h_eta_ctrl->Sumw2();
        TH1D* h_PFiso_dBeta_barrel_nume = new TH1D("h_PFiso_dBeta_barrel_nume", "h_PFiso_dBeta_barrel_nume", 30, 0, 0.3); h_PFiso_dBeta_barrel_nume->Sumw2();
        TH1D* h_PFiso_dBeta_endcap_nume = new TH1D("h_PFiso_dBeta_endcap_nume", "h_PFiso_dBeta_endcap_nume", 30, 0, 0.3); h_PFiso_dBeta_endcap_nume->Sumw2();
        TH1D* h_PFiso_dBeta_barrel_deno = new TH1D("h_PFiso_dBeta_barrel_deno", "h_PFiso_dBeta_barrel_deno", 50, 0, 5); h_PFiso_dBeta_barrel_deno->Sumw2();
        TH1D* h_PFiso_dBeta_endcap_deno = new TH1D("h_PFiso_dBeta_endcap_deno", "h_PFiso_dBeta_endcap_deno", 50, 0, 5); h_PFiso_dBeta_endcap_deno->Sumw2();
        TH1D* h_PFiso_dBeta_barrel_ctrl = new TH1D("h_PFiso_dBeta_barrel_ctrl", "h_PFiso_dBeta_barrel_ctrl", 50, 0, 5); h_PFiso_dBeta_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_dBeta_endcap_ctrl = new TH1D("h_PFiso_dBeta_endcap_ctrl", "h_PFiso_dBeta_endcap_ctrl", 50, 0, 5); h_PFiso_dBeta_endcap_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_barrel_nume = new TH1D("h_PFiso_Rho_barrel_nume", "h_PFiso_Rho_barrel_nume", 20, 0, 0.2); h_PFiso_Rho_barrel_nume->Sumw2();
        TH1D* h_PFiso_Rho_endcap_nume = new TH1D("h_PFiso_Rho_endcap_nume", "h_PFiso_Rho_endcap_nume", 20, 0, 0.2); h_PFiso_Rho_endcap_nume->Sumw2();
        TH1D* h_PFiso_Rho_barrel_deno = new TH1D("h_PFiso_Rho_barrel_deno", "h_PFiso_Rho_barrel_deno", 50, 0, 5); h_PFiso_Rho_barrel_deno->Sumw2();
        TH1D* h_PFiso_Rho_endcap_deno = new TH1D("h_PFiso_Rho_endcap_deno", "h_PFiso_Rho_endcap_deno", 50, 0, 5); h_PFiso_Rho_endcap_deno->Sumw2();
        TH1D* h_PFiso_Rho_barrel_ctrl = new TH1D("h_PFiso_Rho_barrel_ctrl", "h_PFiso_Rho_barrel_ctrl", 50, 0, 5); h_PFiso_Rho_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_endcap_ctrl = new TH1D("h_PFiso_Rho_endcap_ctrl", "h_PFiso_Rho_endcap_ctrl", 50, 0, 5); h_PFiso_Rho_endcap_ctrl->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_nume = new TH1D("h_SigmaIEtaIEta_barrel_nume", "h_SigmaIEtaIEta_barrel_nume", 10, 0, 0.01); h_SigmaIEtaIEta_barrel_nume->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_nume = new TH1D("h_SigmaIEtaIEta_endcap_nume", "h_SigmaIEtaIEta_endcap_nume", 30, 0, 0.03); h_SigmaIEtaIEta_endcap_nume->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_deno = new TH1D("h_SigmaIEtaIEta_barrel_deno", "h_SigmaIEtaIEta_barrel_deno", 50, 0, 0.05); h_SigmaIEtaIEta_barrel_deno->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_deno = new TH1D("h_SigmaIEtaIEta_endcap_deno", "h_SigmaIEtaIEta_endcap_deno", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap_deno->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_ctrl = new TH1D("h_SigmaIEtaIEta_barrel_ctrl", "h_SigmaIEtaIEta_barrel_ctrl", 50, 0, 0.05); h_SigmaIEtaIEta_barrel_ctrl->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_ctrl = new TH1D("h_SigmaIEtaIEta_endcap_ctrl", "h_SigmaIEtaIEta_endcap_ctrl", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap_ctrl->Sumw2();
        TH1D* h_dEtaInSeed_barrel_nume = new TH1D("h_dEtaInSeed_barrel_nume", "h_dEtaInSeed_barrel_nume", 20, -0.1, 0.1); h_dEtaInSeed_barrel_nume->Sumw2();
        TH1D* h_dEtaInSeed_endcap_nume = new TH1D("h_dEtaInSeed_endcap_nume", "h_dEtaInSeed_endcap_nume", 20, -0.1, 0.1); h_dEtaInSeed_endcap_nume->Sumw2();
        TH1D* h_dEtaInSeed_barrel_deno = new TH1D("h_dEtaInSeed_barrel_deno", "h_dEtaInSeed_barrel_deno", 100, -1, 1); h_dEtaInSeed_barrel_deno->Sumw2();
        TH1D* h_dEtaInSeed_endcap_deno = new TH1D("h_dEtaInSeed_endcap_deno", "h_dEtaInSeed_endcap_deno", 100, -1, 1); h_dEtaInSeed_endcap_deno->Sumw2();
        TH1D* h_dEtaInSeed_barrel_ctrl = new TH1D("h_dEtaInSeed_barrel_ctrl", "h_dEtaInSeed_barrel_ctrl", 100, -1, 1); h_dEtaInSeed_barrel_ctrl->Sumw2();
        TH1D* h_dEtaInSeed_endcap_ctrl = new TH1D("h_dEtaInSeed_endcap_ctrl", "h_dEtaInSeed_endcap_ctrl", 100, -1, 1); h_dEtaInSeed_endcap_ctrl->Sumw2();
        TH1D* h_dPhiIn_barrel_nume = new TH1D("h_dPhiIn_barrel_nume", "h_dPhiIn_barrel_nume", 20, -0.1, 0.1); h_dPhiIn_barrel_nume->Sumw2();
        TH1D* h_dPhiIn_endcap_nume = new TH1D("h_dPhiIn_endcap_nume", "h_dPhiIn_endcap_nume", 20, -0.1, 0.1); h_dPhiIn_endcap_nume->Sumw2();
        TH1D* h_dPhiIn_barrel_deno = new TH1D("h_dPhiIn_barrel_deno", "h_dPhiIn_barrel_deno", 20, -0.1, 0.1); h_dPhiIn_barrel_deno->Sumw2();
        TH1D* h_dPhiIn_endcap_deno = new TH1D("h_dPhiIn_endcap_deno", "h_dPhiIn_endcap_deno", 100, -1, 1);    h_dPhiIn_endcap_deno->Sumw2();
        TH1D* h_dPhiIn_barrel_ctrl = new TH1D("h_dPhiIn_barrel_ctrl", "h_dPhiIn_barrel_ctrl", 20, -0.1, 0.1); h_dPhiIn_barrel_ctrl->Sumw2();
        TH1D* h_dPhiIn_endcap_ctrl = new TH1D("h_dPhiIn_endcap_ctrl", "h_dPhiIn_endcap_ctrl", 100, -1, 1);    h_dPhiIn_endcap_ctrl->Sumw2();
        TH1D* h_HoverE_barrel_nume = new TH1D("h_HoverE_barrel_nume", "h_HoverE_barrel_nume", 20, 0, 0.1); h_HoverE_barrel_nume->Sumw2();
        TH1D* h_HoverE_endcap_nume = new TH1D("h_HoverE_endcap_nume", "h_HoverE_endcap_nume", 30, 0, 0.15); h_HoverE_endcap_nume->Sumw2();
        TH1D* h_HoverE_barrel_deno = new TH1D("h_HoverE_barrel_deno", "h_HoverE_barrel_deno", 100, 0, 1);  h_HoverE_barrel_deno->Sumw2();
        TH1D* h_HoverE_endcap_deno = new TH1D("h_HoverE_endcap_deno", "h_HoverE_endcap_deno", 100, 0, 1); h_HoverE_endcap_deno->Sumw2();
        TH1D* h_HoverE_barrel_ctrl = new TH1D("h_HoverE_barrel_ctrl", "h_HoverE_barrel_ctrl", 100, 0, 1);  h_HoverE_barrel_ctrl->Sumw2();
        TH1D* h_HoverE_endcap_ctrl = new TH1D("h_HoverE_endcap_ctrl", "h_HoverE_endcap_ctrl", 100, 0, 1); h_HoverE_endcap_ctrl->Sumw2();
        TH1D* h_InvEminusInvP_barrel_nume = new TH1D("h_InvEminusInvP_barrel_nume", "h_InvEminusInvP_barrel_nume", 100, -0.5, 0.5); h_InvEminusInvP_barrel_nume->Sumw2();
        TH1D* h_InvEminusInvP_endcap_nume = new TH1D("h_InvEminusInvP_endcap_nume", "h_InvEminusInvP_endcap_nume", 100, -0.5, 0.5); h_InvEminusInvP_endcap_nume->Sumw2();
        TH1D* h_InvEminusInvP_barrel_deno = new TH1D("h_InvEminusInvP_barrel_deno", "h_InvEminusInvP_barrel_deno", 120, -6, 6); h_InvEminusInvP_barrel_deno->Sumw2();
        TH1D* h_InvEminusInvP_endcap_deno = new TH1D("h_InvEminusInvP_endcap_deno", "h_InvEminusInvP_endcap_deno", 120, -6, 6); h_InvEminusInvP_endcap_deno->Sumw2();
        TH1D* h_InvEminusInvP_barrel_ctrl = new TH1D("h_InvEminusInvP_barrel_ctrl", "h_InvEminusInvP_barrel_ctrl", 120, -6, 6); h_InvEminusInvP_barrel_ctrl->Sumw2();
        TH1D* h_InvEminusInvP_endcap_ctrl = new TH1D("h_InvEminusInvP_endcap_ctrl", "h_InvEminusInvP_endcap_ctrl", 120, -6, 6); h_InvEminusInvP_endcap_ctrl->Sumw2();
        TH1D* h_chiso_barrel_nume = new TH1D("h_chiso_barrel_nume", "h_chiso_barrel_nume", 20, 0, 0.2); h_chiso_barrel_nume->Sumw2();
        TH1D* h_chiso_endcap_nume = new TH1D("h_chiso_endcap_nume", "h_chiso_endcap_nume", 20, 0, 0.2); h_chiso_endcap_nume->Sumw2();
        TH1D* h_chiso_barrel_deno = new TH1D("h_chiso_barrel_deno", "h_chiso_barrel_deno", 100, 0, 10); h_chiso_barrel_deno->Sumw2();
        TH1D* h_chiso_endcap_deno = new TH1D("h_chiso_endcap_deno", "h_chiso_endcap_deno", 100, 0, 10); h_chiso_endcap_deno->Sumw2();
        TH1D* h_chiso_barrel_ctrl = new TH1D("h_chiso_barrel_ctrl", "h_chiso_barrel_ctrl", 100, 0, 10); h_chiso_barrel_ctrl->Sumw2();
        TH1D* h_chiso_endcap_ctrl = new TH1D("h_chiso_endcap_ctrl", "h_chiso_endcap_ctrl", 100, 0, 10); h_chiso_endcap_ctrl->Sumw2();
        TH1D* h_nhiso_barrel_nume = new TH1D("h_nhiso_barrel_nume", "h_nhiso_barrel_nume", 20, 0, 0.2); h_nhiso_barrel_nume->Sumw2();
        TH1D* h_nhiso_endcap_nume = new TH1D("h_nhiso_endcap_nume", "h_nhiso_endcap_nume", 20, 0, 0.2); h_nhiso_endcap_nume->Sumw2();
        TH1D* h_nhiso_barrel_deno = new TH1D("h_nhiso_barrel_deno", "h_nhiso_barrel_deno", 100, 0, 10); h_nhiso_barrel_deno->Sumw2();
        TH1D* h_nhiso_endcap_deno = new TH1D("h_nhiso_endcap_deno", "h_nhiso_endcap_deno", 100, 0, 10); h_nhiso_endcap_deno->Sumw2();
        TH1D* h_nhiso_barrel_ctrl = new TH1D("h_nhiso_barrel_ctrl", "h_nhiso_barrel_ctrl", 100, 0, 10); h_nhiso_barrel_ctrl->Sumw2();
        TH1D* h_nhiso_endcap_ctrl = new TH1D("h_nhiso_endcap_ctrl", "h_nhiso_endcap_ctrl", 100, 0, 10); h_nhiso_endcap_ctrl->Sumw2();
        TH1D* h_phiso_barrel_nume = new TH1D("h_phiso_barrel_nume", "h_phiso_barrel_nume", 20, 0, 0.2); h_phiso_barrel_nume->Sumw2();
        TH1D* h_phiso_endcap_nume = new TH1D("h_phiso_endcap_nume", "h_phiso_endcap_nume", 20, 0, 0.2); h_phiso_endcap_nume->Sumw2();
        TH1D* h_phiso_barrel_deno = new TH1D("h_phiso_barrel_deno", "h_phiso_barrel_deno", 100, 0, 10); h_phiso_barrel_deno->Sumw2();
        TH1D* h_phiso_endcap_deno = new TH1D("h_phiso_endcap_deno", "h_phiso_endcap_deno", 100, 0, 10); h_phiso_endcap_deno->Sumw2();
        TH1D* h_phiso_barrel_ctrl = new TH1D("h_phiso_barrel_ctrl", "h_phiso_barrel_ctrl", 100, 0, 10); h_phiso_barrel_ctrl->Sumw2();
        TH1D* h_phiso_endcap_ctrl = new TH1D("h_phiso_endcap_ctrl", "h_phiso_endcap_ctrl", 100, 0, 10); h_phiso_endcap_ctrl->Sumw2();
        TH1D* h_chisoPU_barrel_nume = new TH1D("h_chisoPU_barrel_nume", "h_chisoPU_barrel_nume", 10, 0, 1); h_chisoPU_barrel_nume->Sumw2();
        TH1D* h_chisoPU_endcap_nume = new TH1D("h_chisoPU_endcap_nume", "h_chisoPU_endcap_nume", 10, 0, 1); h_chisoPU_endcap_nume->Sumw2();
        TH1D* h_chisoPU_barrel_deno = new TH1D("h_chisoPU_barrel_deno", "h_chisoPU_barrel_deno", 50, 0, 5); h_chisoPU_barrel_deno->Sumw2();
        TH1D* h_chisoPU_endcap_deno = new TH1D("h_chisoPU_endcap_deno", "h_chisoPU_endcap_deno", 50, 0, 5); h_chisoPU_endcap_deno->Sumw2();
        TH1D* h_chisoPU_barrel_ctrl = new TH1D("h_chisoPU_barrel_ctrl", "h_chisoPU_barrel_ctrl", 50, 0, 5); h_chisoPU_barrel_ctrl->Sumw2();
        TH1D* h_chisoPU_endcap_ctrl = new TH1D("h_chisoPU_endcap_ctrl", "h_chisoPU_endcap_ctrl", 50, 0, 5); h_chisoPU_endcap_ctrl->Sumw2();
        TH1D* h_MET = new TH1D("h_MET", "h_MET", 100, 0, 100); h_MET->Sumw2();
        TH1D* h_MT_barrel_nume = new TH1D("h_MT_barrel_nume", "h_MT_barrel_nume", 500, 0, 1000); h_MT_barrel_nume->Sumw2();
        TH1D* h_MT_endcap_nume = new TH1D("h_MT_endcap_nume", "h_MT_endcap_nume", 500, 0, 1000); h_MT_endcap_nume->Sumw2();
        TH1D* h_MT_barrel_deno = new TH1D("h_MT_barrel_deno", "h_MT_barrel_deno", 500, 0, 1000); h_MT_barrel_deno->Sumw2();
        TH1D* h_MT_endcap_deno = new TH1D("h_MT_endcap_deno", "h_MT_endcap_deno", 500, 0, 1000); h_MT_endcap_deno->Sumw2();
        TH1D* h_MT_barrel_ctrl = new TH1D("h_MT_barrel_ctrl", "h_MT_barrel_ctrl", 500, 0, 1000); h_MT_barrel_ctrl->Sumw2();
        TH1D* h_MT_endcap_ctrl = new TH1D("h_MT_endcap_ctrl", "h_MT_endcap_ctrl", 500, 0, 1000); h_MT_endcap_ctrl->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();

        TH1D* h_HoverE_barrel_template_int = new TH1D("h_HoverE_barrel_template_int", "h_HoverE_barrel_template_int", binnum_ele_template, bins_ele_HE_barrel); h_HoverE_barrel_template_int->Sumw2();
        TH1D* h_HoverE_endcap_template_int = new TH1D("h_HoverE_endcap_template_int", "h_HoverE_endcap_template_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap_template_int->Sumw2();
        TH1D* h_HoverE_barrel_jetTemplate_int = new TH1D("h_HoverE_barrel_jetTemplate_int", "h_HoverE_barrel_jetTemplate_int", binnum_ele_template, bins_ele_HE_barrel); h_HoverE_barrel_jetTemplate_int->Sumw2();
        TH1D* h_HoverE_endcap_jetTemplate_int = new TH1D("h_HoverE_endcap_jetTemplate_int", "h_HoverE_endcap_jetTemplate_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap_jetTemplate_int->Sumw2();

        TH2D *h2_HEvsIso_barrel[nPtBin_ele],     *h2_HEvsIso_endcap[nPtBin_ele],
             *h2_ISOvsSigma_barrel[nPtBin_ele],  *h2_ISOvsSigma_endcap[nPtBin_ele],
             *h2_HEvsSigma_barrel[nPtBin_ele],   *h2_HEvsSigma_endcap[nPtBin_ele],
             *h2_ISOvsdEta_barrel[nPtBin_ele],   *h2_ISOvsdEta_endcap[nPtBin_ele],
             *h2_HEvsdEta_barrel[nPtBin_ele],    *h2_HEvsdEta_endcap[nPtBin_ele],
             *h2_ISOvsEP_barrel[nPtBin_ele],     *h2_ISOvsEP_endcap[nPtBin_ele],
             *h2_HEvsEP_barrel[nPtBin_ele],      *h2_HEvsEP_endcap[nPtBin_ele],
             *h2_ISOvsdPhi_barrel[nPtBin_ele],   *h2_ISOvsdPhi_endcap[nPtBin_ele],
             *h2_HEvsdPhi_barrel[nPtBin_ele],    *h2_HEvsdPhi_endcap[nPtBin_ele],
             *h2_SigmavsdEta_barrel[nPtBin_ele], *h2_SigmavsdEta_endcap[nPtBin_ele],
             *h2_SigmavsdPhi_barrel[nPtBin_ele], *h2_SigmavsdPhi_endcap[nPtBin_ele],
             *h2_SigmavsEP_barrel[nPtBin_ele],   *h2_SigmavsEP_endcap[nPtBin_ele],
             *h2_dEtavsdPhi_barrel[nPtBin_ele],  *h2_dEtavsdPhi_endcap[nPtBin_ele],
             *h2_dEtavsEP_barrel[nPtBin_ele],    *h2_dEtavsEP_endcap[nPtBin_ele],
             *h2_dPhivsEP_barrel[nPtBin_ele],    *h2_dPhivsEP_endcap[nPtBin_ele];

        TH2D* h2_HEvsPt_barrel = new TH2D("h2_HEvsPt_barrel", "h2_HEvsPt_barrel", binnum_ele_template, bins_ele_HE_barrel, nPtBin_ele, analyzer->ptbin_ele);
        TH2D* h2_dEtavsPt_barrel = new TH2D("h2_dEtavsPt_barrel", "h2_dEtavsPt_barrel", binnum_ele_template, bins_ele_dEta_barrel, nPtBin_ele, analyzer->ptbin_ele);
        TH2D* h2_SigmavsPt_barrel = new TH2D("h2_SigmavsPt_barrel", "h2_SigmavsPt_barrel", binnum_ele_template, bins_ele_sigma_barrel, nPtBin_ele, analyzer->ptbin_ele);

        TH2D* h2_HEvsPt_endcap = new TH2D("h2_HEvsPt_endcap", "h2_HEvsPt_endcap", binnum_ele_template, bins_ele_HE_endcap, nPtBin_ele, analyzer->ptbin_ele);
        TH2D* h2_dEtavsPt_endcap = new TH2D("h2_dEtavsPt_endcap", "h2_dEtavsPt_endcap", binnum_ele_template, bins_ele_dEta_endcap, nPtBin_ele, analyzer->ptbin_ele);
        TH2D* h2_SigmavsPt_endcap = new TH2D("h2_SigmavsPt_endcap", "h2_SigmavsPt_endcap", binnum_ele_template, bins_ele_sigma_endcap, nPtBin_ele, analyzer->ptbin_ele);

        TH1D *h_PFiso_Rho_barrel_template[nPtBin_ele],       *h_PFiso_Rho_endcap_template[nPtBin_ele],
             *h_PFiso_Rho_barrel_jetTemplate[nPtBin_ele],    *h_PFiso_Rho_endcap_jetTemplate[nPtBin_ele],
             *h_PFiso_Rho_barrel_badJetTemplate[nPtBin_ele], *h_PFiso_Rho_endcap_badJetTemplate[nPtBin_ele],
             *h_HoverE_barrel_template[nPtBin_ele],       *h_HoverE_endcap_template[nPtBin_ele],
             *h_HoverE_barrel_jetTemplate[nPtBin_ele],    *h_HoverE_endcap_jetTemplate[nPtBin_ele];
        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            h_PFiso_Rho_barrel_template[ih] = new TH1D("h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10),
                                                       "h_PFiso_Rho_barrel_template_"+TString::Itoa(ih, 10),
                                                       binnum_ele_template, bins_ele_template_barrel);
            h_PFiso_Rho_barrel_jetTemplate[ih] = new TH1D("h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10),
                                                          "h_PFiso_Rho_barrel_jetTemplate_"+TString::Itoa(ih, 10),
                                                          binnum_ele_template, bins_ele_template_barrel);
            h_PFiso_Rho_barrel_badJetTemplate[ih] = new TH1D("h_PFiso_Rho_barrel_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             "h_PFiso_Rho_barrel_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             binnum_ele_template, bins_ele_template_barrel);
            h_HoverE_barrel_template[ih] = new TH1D("h_HoverE_barrel_template_"+TString::Itoa(ih, 10),
                                                    "h_HoverE_barrel_template_"+TString::Itoa(ih, 10),
                                                    binnum_ele_template, bins_ele_HE_barrel);
            h_HoverE_barrel_jetTemplate[ih] = new TH1D("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10),
                                                       "h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10),
                                                       binnum_ele_template, bins_ele_HE_barrel);
            h_PFiso_Rho_endcap_template[ih] = new TH1D("h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10),
                                                       "h_PFiso_Rho_endcap_template_"+TString::Itoa(ih, 10),
                                                       binnum_ele_template, bins_ele_template_endcap);
            h_PFiso_Rho_endcap_jetTemplate[ih] = new TH1D("h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10),
                                                          "h_PFiso_Rho_endcap_jetTemplate_"+TString::Itoa(ih, 10),
                                                          binnum_ele_template, bins_ele_template_endcap);
            h_PFiso_Rho_endcap_badJetTemplate[ih] = new TH1D("h_PFiso_Rho_endcap_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             "h_PFiso_Rho_endcap_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             binnum_ele_template, bins_ele_template_endcap);
            h_HoverE_endcap_template[ih] = new TH1D("h_HoverE_endcap_template_"+TString::Itoa(ih, 10),
                                                    "h_HoverE_endcap_template_"+TString::Itoa(ih, 10),
                                                    binnum_ele_template, bins_ele_HE_endcap);
            h_HoverE_endcap_jetTemplate[ih] = new TH1D("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10),
                                                       "h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10),
                                                       binnum_ele_template, bins_ele_HE_endcap);
            h_PFiso_Rho_barrel_template[ih]->Sumw2();
            h_PFiso_Rho_barrel_jetTemplate[ih]->Sumw2();
            h_PFiso_Rho_barrel_badJetTemplate[ih]->Sumw2();
            h_HoverE_barrel_template[ih]->Sumw2();
            h_HoverE_barrel_jetTemplate[ih]->Sumw2();
            h_PFiso_Rho_endcap_template[ih]->Sumw2();
            h_PFiso_Rho_endcap_jetTemplate[ih]->Sumw2();
            h_PFiso_Rho_endcap_badJetTemplate[ih]->Sumw2();
            h_HoverE_endcap_template[ih]->Sumw2();
            h_HoverE_endcap_jetTemplate[ih]->Sumw2();

            h2_HEvsIso_barrel[ih] = new TH2D("h2_HEvsIso_barrel_"+TString::Itoa(ih, 10), "h2_HEvsIso_barrel_"+TString::Itoa(ih, 10),
                                             binnum_ele_template, bins_ele_HE_barrel, binnum_ele_template, bins_ele_template_barrel);
            h2_ISOvsSigma_barrel[ih] = new TH2D("h2_ISOvsSigma_barrel_"+TString::Itoa(ih, 10), "h2_ISOvsSigma_barrel_"+TString::Itoa(ih, 10),
                                                binnum_ele_template, bins_ele_template_barrel, binnum_ele_template, bins_ele_sigma_barrel);
            h2_HEvsSigma_barrel[ih] = new TH2D("h2_HEvsSigma_barrel_"+TString::Itoa(ih, 10), "h2_HEvsIso_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_HE_barrel, binnum_ele_template, bins_ele_sigma_barrel);
            h2_ISOvsdEta_barrel[ih] = new TH2D("h2_ISOvsdEta_barrel_"+TString::Itoa(ih, 10), "h2_ISOvsdEta_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_template_barrel, binnum_ele_template, bins_ele_dEta_barrel);
            h2_HEvsdEta_barrel[ih] = new TH2D("h2_HEvsdEta_barrel_"+TString::Itoa(ih, 10), "h2_HEvsdEta_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_HE_barrel, binnum_ele_template, bins_ele_dEta_barrel);
            h2_ISOvsEP_barrel[ih] = new TH2D("h2_ISOvsEP_barrel_"+TString::Itoa(ih, 10), "h2_ISOvsEP_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_template_barrel, binnum_ele_template, bins_ele_EP_barrel);
            h2_HEvsEP_barrel[ih] = new TH2D("h2_HEvsEP_barrel_"+TString::Itoa(ih, 10), "h2_HEvsEP_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_HE_barrel, binnum_ele_template, bins_ele_EP_barrel);          
            h2_ISOvsdPhi_barrel[ih] = new TH2D("h2_ISOvsdPhi_barrel_"+TString::Itoa(ih, 10), "h2_ISOvsdPhi_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_template_barrel, binnum_ele_template, bins_ele_dPhi_barrel);
            h2_HEvsdPhi_barrel[ih] = new TH2D("h2_HEvsdPhi_barrel_"+TString::Itoa(ih, 10), "h2_HEvsdPhi_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_HE_barrel, binnum_ele_template, bins_ele_dPhi_barrel);
            h2_SigmavsdEta_barrel[ih] = new TH2D("h2_SigmavsdEta_barrel_"+TString::Itoa(ih, 10), "h2_SigmavsdEta_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_sigma_barrel, binnum_ele_template, bins_ele_dEta_barrel);
            h2_SigmavsdPhi_barrel[ih] = new TH2D("h2_SigmavsdPhi_barrel_"+TString::Itoa(ih, 10), "h2_SigmavsdPhi_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_sigma_barrel, binnum_ele_template, bins_ele_dPhi_barrel);
            h2_SigmavsEP_barrel[ih] = new TH2D("h2_SigmavsEP_barrel_"+TString::Itoa(ih, 10), "h2_SigmavsEP_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_sigma_barrel, binnum_ele_template, bins_ele_EP_barrel);
            h2_dEtavsdPhi_barrel[ih] = new TH2D("h2_dEtavsdPhi_barrel_"+TString::Itoa(ih, 10), "h2_dEtavsdPhi_barrel_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_dEta_barrel, binnum_ele_template, bins_ele_dPhi_barrel);
            h2_dEtavsEP_barrel[ih] = new TH2D("h2_dEtavsEP_barrel_"+TString::Itoa(ih, 10), "h2_dEtavsEP_barrel_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_dEta_barrel, binnum_ele_template, bins_ele_EP_barrel);
            h2_dPhivsEP_barrel[ih] = new TH2D("h2_dPhivsEP_barrel_"+TString::Itoa(ih, 10), "h2_dPhivsEP_barrel_"+TString::Itoa(ih, 10),
                                             binnum_ele_template, bins_ele_dPhi_barrel, binnum_ele_template, bins_ele_EP_barrel);

            h2_HEvsIso_endcap[ih] = new TH2D("h2_HEvsIso_endcap_"+TString::Itoa(ih, 10), "h2_HEvsIso_endcap_"+TString::Itoa(ih, 10),
                                             binnum_ele_template, bins_ele_HE_endcap, binnum_ele_template, bins_ele_template_endcap);
            h2_ISOvsSigma_endcap[ih] = new TH2D("h2_ISOvsSigma_endcap_"+TString::Itoa(ih, 10), "h2_ISOvsSigma_endcap_"+TString::Itoa(ih, 10),
                                                binnum_ele_template, bins_ele_template_endcap, binnum_ele_template, bins_ele_sigma_endcap);
            h2_HEvsSigma_endcap[ih] = new TH2D("h2_HEvsSigma_endcap_"+TString::Itoa(ih, 10), "h2_HEvsSigma_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_HE_endcap, binnum_ele_template, bins_ele_sigma_endcap);
            h2_ISOvsdEta_endcap[ih] = new TH2D("h2_ISOvsdEta_endcap_"+TString::Itoa(ih, 10), "h2_ISOvsdEta_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_template_endcap, binnum_ele_template, bins_ele_dEta_endcap);
            h2_HEvsdEta_endcap[ih] = new TH2D("h2_HEvsdEta_endcap_"+TString::Itoa(ih, 10), "h2_HEvsdEta_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_HE_endcap, binnum_ele_template, bins_ele_dEta_endcap);
            h2_ISOvsEP_endcap[ih] = new TH2D("h2_ISOvsEP_endcap_"+TString::Itoa(ih, 10), "h2_ISOvsEP_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_template_endcap, binnum_ele_template, bins_ele_EP_endcap);
            h2_HEvsEP_endcap[ih] = new TH2D("h2_HEvsEP_endcap_"+TString::Itoa(ih, 10), "h2_HEvsEP_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_HE_endcap, binnum_ele_template, bins_ele_EP_endcap);
            h2_ISOvsdPhi_endcap[ih] = new TH2D("h2_ISOvsdPhi_endcap_"+TString::Itoa(ih, 10), "h2_ISOvsdPhi_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_template_endcap, binnum_ele_template, bins_ele_dPhi_endcap);
            h2_HEvsdPhi_endcap[ih] = new TH2D("h2_HEvsdPhi_endcap_"+TString::Itoa(ih, 10), "h2_HEvsdPhi_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_HE_endcap, binnum_ele_template, bins_ele_dPhi_endcap);
            h2_SigmavsdEta_endcap[ih] = new TH2D("h2_SigmavsdEta_endcap_"+TString::Itoa(ih, 10), "h2_SigmavsdEta_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_sigma_endcap, binnum_ele_template, bins_ele_dEta_endcap);
            h2_SigmavsdPhi_endcap[ih] = new TH2D("h2_SigmavsdPhi_endcap_"+TString::Itoa(ih, 10), "h2_SigmavsdPhi_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_sigma_endcap, binnum_ele_template, bins_ele_dPhi_endcap);
            h2_SigmavsEP_endcap[ih] = new TH2D("h2_SigmavsEP_endcap_"+TString::Itoa(ih, 10), "h2_SigmavsEP_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_sigma_endcap, binnum_ele_template, bins_ele_EP_endcap);
            h2_dEtavsdPhi_endcap[ih] = new TH2D("h2_dEtavsdPhi_endcap_"+TString::Itoa(ih, 10), "h2_dEtavsdPhi_endcap_"+TString::Itoa(ih, 10),
                                              binnum_ele_template, bins_ele_dEta_endcap, binnum_ele_template, bins_ele_dPhi_endcap);
            h2_dEtavsEP_endcap[ih] = new TH2D("h2_dEtavsEP_endcap_"+TString::Itoa(ih, 10), "h2_dEtavsEP_endcap_"+TString::Itoa(ih, 10),
                                               binnum_ele_template, bins_ele_dEta_endcap, binnum_ele_template, bins_ele_EP_endcap);
            h2_dPhivsEP_endcap[ih] = new TH2D("h2_dPhivsEP_endcap_"+TString::Itoa(ih, 10), "h2_dPhivsEP_endcap_"+TString::Itoa(ih, 10),
                                             binnum_ele_template, bins_ele_dPhi_endcap, binnum_ele_template, bins_ele_EP_endcap);

        }

        TH1D* h_PFiso_Rho_barrel_template_PU10to20 = new TH1D("h_PFiso_Rho_barrel_template_PU10to20", "h_PFiso_Rho_barrel_template_PU10to20", 50, 0, 5); h_PFiso_Rho_barrel_template_PU10to20->Sumw2();
        TH1D* h_PFiso_Rho_endcap_template_PU10to20 = new TH1D("h_PFiso_Rho_endcap_template_PU10to20", "h_PFiso_Rho_endcap_template_PU10to20", 50, 0, 5); h_PFiso_Rho_endcap_template_PU10to20->Sumw2();
        TH1D* h_PFiso_Rho_barrel_template_PU20to30 = new TH1D("h_PFiso_Rho_barrel_template_PU20to30", "h_PFiso_Rho_barrel_template_PU20to30", 50, 0, 5); h_PFiso_Rho_barrel_template_PU20to30->Sumw2();
        TH1D* h_PFiso_Rho_endcap_template_PU20to30 = new TH1D("h_PFiso_Rho_endcap_template_PU20to30", "h_PFiso_Rho_endcap_template_PU20to30", 50, 0, 5); h_PFiso_Rho_endcap_template_PU20to30->Sumw2();
        TH1D* h_PFiso_Rho_barrel_template_PU30to50 = new TH1D("h_PFiso_Rho_barrel_template_PU30to50", "h_PFiso_Rho_barrel_template_PU30to50", 50, 0, 5); h_PFiso_Rho_barrel_template_PU30to50->Sumw2();
        TH1D* h_PFiso_Rho_endcap_template_PU30to50 = new TH1D("h_PFiso_Rho_endcap_template_PU30to50", "h_PFiso_Rho_endcap_template_PU30to50", 50, 0, 5); h_PFiso_Rho_endcap_template_PU30to50->Sumw2();

        TH1D* h_mass_test = new TH1D("h_mass_test", "h_mass_test", binnum, massbins); h_mass_test->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<double> *chIso03 = new std::vector<double>;
        std::vector<double> *nhIso03 = new std::vector<double>;
        std::vector<double> *phIso03 = new std::vector<double>;
        std::vector<double> *ChIso03FromPU = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<double> *relPFiso_dBeta = new std::vector<double>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t runNum;
        Int_t lumiBlock;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("InvEminusInvP", 1);
        chain->SetBranchStatus("chIso03", 1);
        chain->SetBranchStatus("nhIso03", 1);
        chain->SetBranchStatus("phIso03", 1);
        chain->SetBranchStatus("ChIso03FromPU", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("relPFiso_dBeta", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("prescale_factor", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("chIso03", &chIso03);
        chain->SetBranchAddress("nhIso03", &nhIso03);
        chain->SetBranchAddress("phIso03", &phIso03);
        chain->SetBranchAddress("ChIso03FromPU", &ChIso03FromPU);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("relPFiso_dBeta", &relPFiso_dBeta);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE){
                cout << "\nEvt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
                cout << "eta[1] = " << eta->at(0) << endl;
                cout << "phi[1] = " << phi->at(0) << endl;
            }

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE) TopPtWeight = top_weight;

            // -- Normalization -- //
            Double_t TotWeight = gen_weight;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            /*///////////////////////////
            if (pr == _WJets) TotWeight = (Lumi * Mgr.Xsec[0] / (16433848+161144203)) * gen_weight; // REMOVE LATER
            ///////////////////////////*/
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl;

            if (Mgr.isMC == kTRUE && p_T->size() > 1) n2MC += TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight;
            if (Mgr.isMC == kFALSE && p_T->size() > 1) n2Data++;

            // For finding the leading electron
            TLorentzVector ele_lead;
            ele_lead.SetPtEtaPhiM(0, 0, 0, M_Elec);
            Double_t iso_lead = -9999;

            if (p_T->size() != passMediumID->size())
            {
                cout << "ERROR: vector sizes do not match!" << endl;
                break;
            }

            Double_t med_count = 0;
            Int_t i_lead = -1;
            for (UInt_t i_ele=0; i_ele<p_T->size(); i_ele++)
            {
                if (passMediumID->at(i_ele)) med_count++;
                if (p_T->at(i_ele) <= 25) continue;
//                if (fabs(eta->at(i_ele)) < 1.4442)
//                {
//                    if (Full5x5_SigmaIEtaIEta->at(i_ele) > 0.013) continue;
//                    if (HoverE->at(i_ele) > 0.13) continue;
//                    if (fabs(dEtaInSeed->at(i_ele)) > 0.01) continue;
//                    if (fabs(dPhiIn->at(i_ele)) > 0.07) continue;
//                }
//                else if (fabs(eta->at(i_ele)) > 1.566)
//                {
//                    if (Full5x5_SigmaIEtaIEta->at(i_ele) > 0.035) continue;
//                    if (HoverE->at(i_ele) > 0.13) continue;
//                }

                // Selecting leading electron (could also try finding a muon with the best isolation)
                if (p_T->at(i_ele) > ele_lead.Pt())
                {
                    Int_t matched = 0;
                    for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                    {
                        if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
                        {
                            matched = 1;
                            i_lead = i_tr;
                        }
                    }
                    if (matched == 0) continue;
                    ele_lead.SetPtEtaPhiM(p_T->at(i_ele), eta->at(i_ele), phi->at(i_ele), M_Elec);
                    iso_lead = relPFiso_dBeta->at(i_ele);
                }
            }
            if (i_lead<0) continue;

            Double_t dTheta = ele_lead.Phi() - MET_phi;
            Double_t MT = sqrt(2 * ele_lead.Pt() * MET_pT * (1 - cos(dTheta)));
//            if (MT >= 60) continue;
//            if (MET_pT >= 20) continue;

            Double_t prescale_lead = analyzer->getPrescale_alt(trig_pT->at(i_lead));
            if (Mgr.isMC == kTRUE) prescale_lead = 1;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);

            if (p_T->size() == 2 && passMediumID->at(0) && passMediumID->at(1))
            {
                TLorentzVector ele1, ele2;
                ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
                ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
                Electron e1, e2;
                e1.Pt  = p_T->at(0);
                e1.etaSC = eta->at(0);
                e2.Pt  = p_T->at(1);
                e2.etaSC = eta->at(1);
//                effweight = analyzer->EfficiencySF_EventWeight_electron(e1, e2);
                h_mass_test->Fill((ele1+ele2).M(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
            }

            if (med_count > 1) continue;

            if (DEBUG == kTRUE)
            {
                cout << "Triggers:" << endl;
                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
                    cout << "Photon" << trig_fired->at(i_tr) << "  p_T: " << trig_pT->at(i_tr) << "  matched to " << trig_matched->at(i_tr) << endl;
                }
            }

            for (UInt_t i_ele=0; i_ele<p_T->size(); i_ele++)
            {
                if (p_T->at(i_ele) != p_T->at(i_ele))
                {
                    cout << p_T->at(i_ele) << " " << eta->at(i_ele) << " " << phi->at(i_ele) << " " << relPFiso_dBeta->at(i_ele) << endl;
                    continue;
                }
                if (p_T->at(i_ele) <= 25) continue;
//                if (fabs(eta->at(i_ele)) < 1.4442)
//                {
//                    if (Full5x5_SigmaIEtaIEta->at(i_ele) > 0.013) continue;
//                    if (HoverE->at(i_ele) > 0.13) continue;
//                    if (fabs(dEtaInSeed->at(i_ele)) > 0.01) continue;
//                    if (fabs(dPhiIn->at(i_ele)) > 0.07) continue;
//                }
//                else if (fabs(eta->at(i_ele)) > 1.566)
//                {
//                    if (Full5x5_SigmaIEtaIEta->at(i_ele) > 0.035) continue;
//                    if (HoverE->at(i_ele) > 0.13) continue;
//                }

                if (DEBUG == kTRUE) cout << "i_ele = " << i_ele << endl;

                Int_t matched = 0;
                Int_t matched22=0, matched30=0, matched36=0, matched50=0, matched75=0, matched90=0, matched120=0, matched175=0;
                Int_t i_22=-1, i_30=-1, i_36=-1, i_50=-1, i_75=-1, i_90=-1, i_120=-1, i_175=-1;
                Double_t prescale_alt = 1.;

                if (Mgr.isMC == kFALSE)
                {
                    for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                    {
//                        if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
//                        {
    //                        matched = 1;
    //                        if (Mgr.isMC == kFALSE) prescale_alt = analyzer->getPrescale_alt(trig_pT->at(i_tr));
//                        }
                        if (trig_fired->at(i_tr) < 22) continue;
                        if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
                        {
                            if (trig_fired->at(i_tr) == 22 && trig_pT->at(i_tr)<30)
                            {
                                i_22 = i_tr;
                                matched22 = 1;
                            }
                            if (trig_fired->at(i_tr) == 30 && trig_pT->at(i_tr)<36)
                            {
                                i_30 = i_tr;
                                matched30 = 1;
                            }
                            if (trig_fired->at(i_tr) == 36 && trig_pT->at(i_tr)<50)
                            {
                                i_36 = i_tr;
                                matched36 = 1;
                            }
                            if (trig_fired->at(i_tr) == 50 && trig_pT->at(i_tr)<75)
                            {
                                i_50 = i_tr;
                                matched50 = 1;
                            }
                            if (trig_fired->at(i_tr) == 75 && trig_pT->at(i_tr)<90)
                            {
                                i_75 = i_tr;
                                matched75 = 1;
                            }
                            if (trig_fired->at(i_tr) == 90 && trig_pT->at(i_tr)<120)
                            {
                                i_90 = i_tr;
                                matched90 = 1;
                            }
                            if (trig_fired->at(i_tr) == 120 && trig_pT->at(i_tr)<175)
                            {
                                i_120 = i_tr;
                                matched120 = 1;
                            }
                            else if (trig_fired->at(i_tr) == 175)
                            {
                                i_175 = i_tr;
                                matched175 = 1;
                            }
                        }

    //                if (matched == 0) continue;
                    }
                    if (matched22==0 && matched30==0 && matched36==0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0) continue;
                    else if (matched22 == 1 && matched30 == 0 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 813./18621470.;
                    }
                    else if (matched30 == 1 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 3211./18621470.;
                    }
                    else if (matched36 == 1 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 6372./18621470.;
                    }
                    else if (matched50 == 1 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 12648./18621470.;
                    }
                    else if (matched75 == 1 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 63170./18621470.;
                    }
                    else if (matched90 == 1 && matched120 == 0 && matched175 == 0)
                    {
                        prescale_alt = 126981./18621470.;
                    }
                    else if (matched120 == 1 && matched175 == 0)
                    {
                        prescale_alt = 260278./18621470.;
                    }
                    else if (matched175 == 1)
                    {
                        prescale_alt = 1.;
                    }
                }


                if (passMediumID->at(i_ele)) // Signal/Numerator
                {
                    h_eta_nume->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    if (fabs(eta->at(i_ele)) < 1.4442) // Barrel
                    {
                        h_pT_barrel_nume->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_barrel_nume->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_barrel_nume->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_barrel_nume->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_barrel_nume->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_barrel_nume->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_barrel_nume->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_barrel_nume->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_barrel_nume->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_barrel_nume->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_barrel_nume->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_barrel_nume->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) > 1.566) // Endcap
                    {
                        h_pT_endcap_nume->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_endcap_nume->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_endcap_nume->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_endcap_nume->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_endcap_nume->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_endcap_nume->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_endcap_nume->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_endcap_nume->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_endcap_nume->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_endcap_nume->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_endcap_nume->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_endcap_nume->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                } // End of if(Signal/Numerator)
                else // Control
                {
                    h_eta_ctrl->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    if (fabs(eta->at(i_ele)) < 1.4442 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.013 && HoverE->at(i_ele) < 0.13 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.01 && fabs(dPhiIn->at(i_ele)) < 0.07 && mHits->at(i_ele) <= 1) // Barrel
                    {
                        h_pT_barrel_ctrl->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_barrel_ctrl->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_barrel_ctrl->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_barrel_ctrl->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_barrel_ctrl->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_barrel_ctrl->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_barrel_ctrl->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_barrel_ctrl->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_barrel_ctrl->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_barrel_ctrl->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_barrel_ctrl->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_barrel_ctrl->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) > 1.566 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 && HoverE->at(i_ele) < 0.13 && mHits->at(i_ele) <= 1) // Endcap
                    {
                        h_pT_endcap_ctrl->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_endcap_ctrl->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_endcap_ctrl->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_endcap_ctrl->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_endcap_ctrl->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_endcap_ctrl->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_endcap_ctrl->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_endcap_ctrl->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_endcap_ctrl->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_endcap_ctrl->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_endcap_ctrl->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_endcap_ctrl->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                }// End of if(Control)
                // Denominator
                if (fabs(eta->at(i_ele)) < 1.4442) // Barrel
                {
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.013 && HoverE->at(i_ele) < 0.13 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.01 && fabs(dPhiIn->at(i_ele)) < 0.07 && mHits->at(i_ele) <= 1)
                    {
                        h_eta_deno->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_pT_barrel_deno->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_barrel_deno->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_barrel_deno->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_barrel_deno->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_barrel_deno->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_barrel_deno->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_barrel_deno->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_barrel_deno->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_barrel_deno->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_barrel_deno->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_barrel_deno->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_barrel_deno->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        HoverE->at(i_ele) < 0.253 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_template[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (HoverE->at(i_ele) < 0.253 && (fabs(InvEminusInvP->at(i_ele)) >= 0.134 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.00998 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00311 || fabs(dPhiIn->at(i_ele)) >= 0.103 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_jetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                             HoverE->at(i_ele) >= 0.253 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_badJetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (relPFiso_Rho->at(i_ele) < 0.0695 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 &&
                        fabs(dPhiIn->at(i_ele)) < 0.103 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h_HoverE_barrel_template_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h2_HEvsPt_barrel->Fill(HoverE->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_barrel_template[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (relPFiso_Rho->at(i_ele) < 0.0695 && (fabs(InvEminusInvP->at(i_ele)) >= 0.134 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.00998 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00311 || fabs(dPhiIn->at(i_ele)) >= 0.103 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        h_HoverE_barrel_jetTemplate_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_barrel_jetTemplate[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (HoverE->at(i_ele) < 0.253 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        relPFiso_Rho->at(i_ele) < 0.0695 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h2_SigmavsPt_barrel->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    if (HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        relPFiso_Rho->at(i_ele) < 0.0695 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h2_dEtavsPt_barrel->Fill(dEtaInSeed->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsIso_barrel[ih]->Fill(HoverE->at(i_ele), relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.253 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsSigma_barrel[ih]->Fill(relPFiso_Rho->at(i_ele), Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsSigma_barrel[ih]->Fill(HoverE->at(i_ele), Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsdEta_barrel[ih]->Fill(relPFiso_Rho->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsdEta_barrel[ih]->Fill(HoverE->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00311 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsEP_barrel[ih]->Fill(relPFiso_Rho->at(i_ele), InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00311 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsEP_barrel[ih]->Fill(HoverE->at(i_ele), InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }           
                    if (HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsdPhi_barrel[ih]->Fill(relPFiso_Rho->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsdPhi_barrel[ih]->Fill(HoverE->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.253 && relPFiso_Rho->at(i_ele) < 0.0695 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsdEta_barrel[ih]->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && HoverE->at(i_ele) < 0.253 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsdPhi_barrel[ih]->Fill(HoverE->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.253 && relPFiso_Rho->at(i_ele) < 0.0695 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00311 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsEP_barrel[ih]->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dEtavsdPhi_barrel[ih]->Fill(fabs(dEtaInSeed->at(i_ele)), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 &&
                        fabs(dPhiIn->at(i_ele)) < 0.103 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dEtavsEP_barrel[ih]->Fill(fabs(dEtaInSeed->at(i_ele)), fabs(InvEminusInvP->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0695 && HoverE->at(i_ele) < 0.253 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00311 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dPhivsEP_barrel[ih]->Fill(fabs(dPhiIn->at(i_ele)), fabs(InvEminusInvP->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }


                    if (nVTX > 10 && nVTX < 20)
                    {
                        h_PFiso_Rho_barrel_template_PU10to20->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (nVTX > 20 && nVTX < 30)
                    {
                       h_PFiso_Rho_barrel_template_PU20to30->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (nVTX > 30 && nVTX < 50)
                    {
                        h_PFiso_Rho_barrel_template_PU30to50->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                }
                else if (fabs(eta->at(i_ele)) > 1.566) // Endcap
                {
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 && HoverE->at(i_ele) < 0.13 && mHits->at(i_ele) <= 1)
                    {
                        h_eta_deno->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_pT_endcap_deno->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_dBeta_endcap_deno->Fill(relPFiso_dBeta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_PFiso_Rho_endcap_deno->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_SigmaIEtaIEta_endcap_deno->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dEtaInSeed_endcap_deno->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_dPhiIn_endcap_deno->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_HoverE_endcap_deno->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_InvEminusInvP_endcap_deno->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chiso_endcap_deno->Fill(chIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_nhiso_endcap_deno->Fill(nhIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_phiso_endcap_deno->Fill(phIso03->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h_chisoPU_endcap_deno->Fill(ChIso03FromPU->at(i_ele)/p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_template[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (HoverE->at(i_ele) < 0.0878 && (fabs(InvEminusInvP->at(i_ele)) >= 0.13 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00609 || fabs(dPhiIn->at(i_ele)) >= 0.045 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_jetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                             HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) >= 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_badJetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        relPFiso_Rho->at(i_ele) < 0.0821 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h_HoverE_endcap_template_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        h2_HEvsPt_endcap->Fill(HoverE->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap_template[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    else if (relPFiso_Rho->at(i_ele) < 0.0821 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 || fabs(dEtaInSeed->at(i_ele)) >= 0.00609 ||
                             fabs(dPhiIn->at(i_ele)) >= 0.045 || fabs(InvEminusInvP->at(i_ele)) >= 0.13 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        h_HoverE_endcap_jetTemplate_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap_jetTemplate[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (HoverE->at(i_ele) < 0.0878 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        relPFiso_Rho->at(i_ele) < 0.0821 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h2_SigmavsPt_barrel->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    if (HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        relPFiso_Rho->at(i_ele) < 0.0821 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h2_dEtavsPt_barrel->Fill(dEtaInSeed->at(i_ele), p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0045103 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsIso_endcap[ih]->Fill(HoverE->at(i_ele), relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.0878 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsSigma_endcap[ih]->Fill(relPFiso_Rho->at(i_ele), Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsSigma_endcap[ih]->Fill(HoverE->at(i_ele), Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsdEta_endcap[ih]->Fill(relPFiso_Rho->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsdEta_endcap[ih]->Fill(HoverE->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00609 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsEP_endcap[ih]->Fill(relPFiso_Rho->at(i_ele), InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00609 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsEP_endcap[ih]->Fill(HoverE->at(i_ele), InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_ISOvsdPhi_endcap[ih]->Fill(relPFiso_Rho->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_HEvsdPhi_endcap[ih]->Fill(HoverE->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (HoverE->at(i_ele) < 0.0878 && relPFiso_Rho->at(i_ele) < 0.0821 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsdEta_endcap[ih]->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && HoverE->at(i_ele) < 0.0878 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsdPhi_endcap[ih]->Fill(HoverE->at(i_ele), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (HoverE->at(i_ele) < 0.0878 && relPFiso_Rho->at(i_ele) < 0.0821 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00609 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_SigmavsEP_endcap[ih]->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 &&
                        fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dEtavsdPhi_endcap[ih]->Fill(fabs(dEtaInSeed->at(i_ele)), fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 &&
                        fabs(dPhiIn->at(i_ele)) < 0.045 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dEtavsEP_endcap[ih]->Fill(fabs(dEtaInSeed->at(i_ele)), fabs(InvEminusInvP->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }
                    if (relPFiso_Rho->at(i_ele) < 0.0821 && HoverE->at(i_ele) < 0.0878 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.00609 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h2_dPhivsEP_endcap[ih]->Fill(fabs(dPhiIn->at(i_ele)), fabs(InvEminusInvP->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                        }
                    }

                    if (nVTX > 10 && nVTX < 20)
                    {
                        h_PFiso_Rho_endcap_template_PU10to20->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (nVTX > 20 && nVTX < 30)
                    {
                       h_PFiso_Rho_endcap_template_PU20to30->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                    else if (nVTX > 30 && nVTX < 50)
                    {
                        h_PFiso_Rho_endcap_template_PU30to50->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    }
                }
            }// End of i_ele iteration

            if (fabs(ele_lead.Eta()) < 1.4442) // Barrel
            {
                h_MT_barrel_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
                if (iso_lead < 0.15) h_MT_barrel_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
                else h_MT_barrel_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
            }
            else if (fabs(ele_lead.Eta()) > 1.566) // Endcap
            {
                h_MT_endcap_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
                if (iso_lead < 0.15) h_MT_endcap_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
                else h_MT_endcap_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_lead);
            }

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

        if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);

        f->cd();
        cout << "\tWriting into file...";

        h_pT_barrel_nume->Write();
        h_pT_endcap_nume->Write();
        h_pT_barrel_deno->Write();
        h_pT_endcap_deno->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_eta_nume->Write();
        h_eta_deno->Write();
        h_eta_ctrl->Write();
        h_PFiso_dBeta_barrel_nume->Write();
        h_PFiso_dBeta_endcap_nume->Write();
        h_PFiso_dBeta_barrel_deno->Write();
        h_PFiso_dBeta_endcap_deno->Write();
        h_PFiso_dBeta_barrel_ctrl->Write();
        h_PFiso_dBeta_endcap_ctrl->Write();
        h_PFiso_Rho_barrel_nume->Write();
        h_PFiso_Rho_endcap_nume->Write();
        h_PFiso_Rho_barrel_deno->Write();
        h_PFiso_Rho_endcap_deno->Write();
        h_PFiso_Rho_barrel_ctrl->Write();
        h_PFiso_Rho_endcap_ctrl->Write();
        h_SigmaIEtaIEta_barrel_nume->Write();
        h_SigmaIEtaIEta_endcap_nume->Write();
        h_SigmaIEtaIEta_barrel_deno->Write();
        h_SigmaIEtaIEta_endcap_deno->Write();
        h_SigmaIEtaIEta_barrel_ctrl->Write();
        h_SigmaIEtaIEta_endcap_ctrl->Write();
        h_dEtaInSeed_barrel_nume->Write();
        h_dEtaInSeed_endcap_nume->Write();
        h_dEtaInSeed_barrel_deno->Write();
        h_dEtaInSeed_endcap_deno->Write();
        h_dEtaInSeed_barrel_ctrl->Write();
        h_dEtaInSeed_endcap_ctrl->Write();
        h_dPhiIn_barrel_nume->Write();
        h_dPhiIn_endcap_nume->Write();
        h_dPhiIn_barrel_deno->Write();
        h_dPhiIn_endcap_deno->Write();
        h_dPhiIn_barrel_ctrl->Write();
        h_dPhiIn_endcap_ctrl->Write();
        h_HoverE_barrel_nume->Write();
        h_HoverE_endcap_nume->Write();
        h_HoverE_barrel_deno->Write();
        h_HoverE_endcap_deno->Write();
        h_HoverE_barrel_ctrl->Write();
        h_HoverE_endcap_ctrl->Write();
        h_InvEminusInvP_barrel_nume->Write();
        h_InvEminusInvP_endcap_nume->Write();
        h_InvEminusInvP_barrel_deno->Write();
        h_InvEminusInvP_endcap_deno->Write();
        h_InvEminusInvP_barrel_ctrl->Write();
        h_InvEminusInvP_endcap_ctrl->Write();
        h_chiso_barrel_nume->Write();
        h_chiso_endcap_nume->Write();
        h_chiso_barrel_deno->Write();
        h_chiso_endcap_deno->Write();
        h_chiso_barrel_ctrl->Write();
        h_chiso_endcap_ctrl->Write();
        h_nhiso_barrel_nume->Write();
        h_nhiso_endcap_nume->Write();
        h_nhiso_barrel_deno->Write();
        h_nhiso_endcap_deno->Write();
        h_nhiso_barrel_ctrl->Write();
        h_nhiso_endcap_ctrl->Write();
        h_phiso_barrel_nume->Write();
        h_phiso_endcap_nume->Write();
        h_phiso_barrel_deno->Write();
        h_phiso_endcap_deno->Write();
        h_phiso_barrel_ctrl->Write();
        h_phiso_endcap_ctrl->Write();
        h_chisoPU_barrel_nume->Write();
        h_chisoPU_endcap_nume->Write();
        h_chisoPU_barrel_deno->Write();
        h_chisoPU_endcap_deno->Write();
        h_chisoPU_barrel_ctrl->Write();
        h_chisoPU_endcap_ctrl->Write();
        h_MET->Write();
        h_MT_barrel_nume->Write();
        h_MT_endcap_nume->Write();
        h_MT_barrel_deno->Write();
        h_MT_endcap_deno->Write();
        h_MT_barrel_ctrl->Write();
        h_MT_endcap_ctrl->Write();
        h_nVTX->Write();

        h_HoverE_barrel_template_int->Write();
        h_HoverE_barrel_jetTemplate_int->Write();
        h_HoverE_endcap_template_int->Write();
        h_HoverE_endcap_jetTemplate_int->Write();
        h2_HEvsPt_barrel->Write();
        h2_SigmavsPt_barrel->Write();
        h2_dEtavsPt_barrel->Write();
        h2_HEvsPt_endcap->Write();
        h2_SigmavsPt_endcap->Write();
        h2_dEtavsPt_endcap->Write();

        for (Int_t ih=0; ih<nPtBin_ele; ih++)
        {
            h_PFiso_Rho_barrel_template[ih]->Write();
            h_PFiso_Rho_barrel_jetTemplate[ih]->Write();
            h_PFiso_Rho_barrel_badJetTemplate[ih]->Write();
            h_HoverE_barrel_template[ih]->Write();
            h_HoverE_barrel_jetTemplate[ih]->Write();

            h_PFiso_Rho_endcap_template[ih]->Write();
            h_PFiso_Rho_endcap_jetTemplate[ih]->Write();
            h_PFiso_Rho_endcap_badJetTemplate[ih]->Write();
            h_HoverE_endcap_template[ih]->Write();
            h_HoverE_endcap_jetTemplate[ih]->Write();

            h2_HEvsIso_barrel[ih]->Write();
            h2_ISOvsSigma_barrel[ih]->Write();
            h2_HEvsSigma_barrel[ih]->Write();
            h2_ISOvsdEta_barrel[ih]->Write();
            h2_HEvsdEta_barrel[ih]->Write();
            h2_HEvsEP_barrel[ih]->Write();
            h2_ISOvsEP_barrel[ih]->Write();
            h2_HEvsdPhi_barrel[ih]->Write();
            h2_ISOvsdPhi_barrel[ih]->Write();
            h2_SigmavsdEta_barrel[ih]->Write();
            h2_SigmavsdPhi_barrel[ih]->Write();
            h2_SigmavsEP_barrel[ih]->Write();
            h2_dEtavsdPhi_barrel[ih]->Write();
            h2_dEtavsEP_barrel[ih]->Write();
            h2_dPhivsEP_barrel[ih]->Write();

            h2_HEvsIso_endcap[ih]->Write();
            h2_ISOvsSigma_endcap[ih]->Write();
            h2_HEvsSigma_endcap[ih]->Write();
            h2_ISOvsdEta_endcap[ih]->Write();
            h2_HEvsdEta_endcap[ih]->Write();
            h2_HEvsEP_endcap[ih]->Write();
            h2_ISOvsEP_endcap[ih]->Write();
            h2_HEvsdPhi_endcap[ih]->Write();
            h2_ISOvsdPhi_endcap[ih]->Write();
            h2_SigmavsdEta_endcap[ih]->Write();
            h2_SigmavsdPhi_endcap[ih]->Write();
            h2_SigmavsEP_endcap[ih]->Write();
            h2_dEtavsdPhi_endcap[ih]->Write();
            h2_dEtavsEP_endcap[ih]->Write();
            h2_dPhivsEP_endcap[ih]->Write();
        }

        h_PFiso_Rho_barrel_template_PU10to20->Write();
        h_PFiso_Rho_endcap_template_PU10to20->Write();
        h_PFiso_Rho_barrel_template_PU20to30->Write();
        h_PFiso_Rho_endcap_template_PU20to30->Write();
        h_PFiso_Rho_barrel_template_PU30to50->Write();
        h_PFiso_Rho_endcap_template_PU30to50->Write();

        h_mass_test->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30
        if (pr == _GJets_2000to5000) pr = _EndOf_SingleElectron_Normal; // next -- SinglePhoton_B

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Total weighted events with 2 or more denominator electrons in MC: " << n2MC << endl;
    cout << "Total events with 2 or more denominator electrons in Data: " << n2Data << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FR_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_FR_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    UInt_t n2MC=0, n2Data=0;

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 48, -2.4, 2.4); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 48, -2.4, 2.4); h_eta_deno->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 48, -2.4, 2.4); h_eta_ctrl->Sumw2();
        TH1D* h_PFiso_barrel_nume = new TH1D("h_PFiso_barrel_nume", "h_PFiso_barrel_nume", 30, 0, 0.15); h_PFiso_barrel_nume->Sumw2();
        TH1D* h_PFiso_endcap_nume = new TH1D("h_PFiso_endcap_nume", "h_PFiso_endcap_nume", 30, 0, 0.15); h_PFiso_endcap_nume->Sumw2();
        TH1D* h_PFiso_barrel_deno = new TH1D("h_PFiso_barrel_deno", "h_PFiso_barrel_deno", 100, 0, 5); h_PFiso_barrel_deno->Sumw2();
        TH1D* h_PFiso_endcap_deno = new TH1D("h_PFiso_endcap_deno", "h_PFiso_endcap_deno", 100, 0, 5); h_PFiso_endcap_deno->Sumw2();
        TH1D* h_PFiso_barrel_ctrl = new TH1D("h_PFiso_barrel_ctrl", "h_PFiso_barrel_ctrl", 100, 0, 5); h_PFiso_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_endcap_ctrl = new TH1D("h_PFiso_endcap_ctrl", "h_PFiso_endcap_ctrl", 100, 0, 5); h_PFiso_endcap_ctrl->Sumw2();
        TH1D* h_TRKiso_barrel_nume = new TH1D("h_TRKiso_barrel_nume", "h_TRKiso_barrel_nume", 50, 0, 1); h_TRKiso_barrel_nume->Sumw2();
        TH1D* h_TRKiso_endcap_nume = new TH1D("h_TRKiso_endcap_nume", "h_TRKiso_endcap_nume", 50, 0, 1); h_TRKiso_endcap_nume->Sumw2();
        TH1D* h_TRKiso_barrel_deno = new TH1D("h_TRKiso_barrel_deno", "h_TRKiso_barrel_deno", 100, 0, 5); h_TRKiso_barrel_deno->Sumw2();
        TH1D* h_TRKiso_endcap_deno = new TH1D("h_TRKiso_endcap_deno", "h_TRKiso_endcap_deno", 100, 0, 5); h_TRKiso_endcap_deno->Sumw2();
        TH1D* h_TRKiso_barrel_ctrl = new TH1D("h_TRKiso_barrel_ctrl", "h_TRKiso_barrel_ctrl", 50, 0.15, 5); h_TRKiso_barrel_deno->Sumw2();
        TH1D* h_TRKiso_endcap_ctrl = new TH1D("h_TRKiso_endcap_ctrl", "h_TRKiso_endcap_ctrl", 50, 0.15, 5); h_TRKiso_endcap_deno->Sumw2();
        TH1D* h_MET = new TH1D("h_MET", "h_MET", 100, 0, 1000); h_MET->Sumw2();
        TH1D* h_MT_barrel_nume = new TH1D("h_MT_barrel_nume", "h_MT_barrel_nume", 500, 0, 1000); h_MT_barrel_nume->Sumw2();
        TH1D* h_MT_endcap_nume = new TH1D("h_MT_endcap_nume", "h_MT_endcap_nume", 500, 0, 1000); h_MT_endcap_nume->Sumw2();
        TH1D* h_MT_barrel_deno = new TH1D("h_MT_barrel_deno", "h_MT_barrel_deno", 500, 0, 1000); h_MT_barrel_deno->Sumw2();
        TH1D* h_MT_endcap_deno = new TH1D("h_MT_endcap_deno", "h_MT_endcap_deno", 500, 0, 1000); h_MT_endcap_deno->Sumw2();
        TH1D* h_MT_barrel_ctrl = new TH1D("h_MT_barrel_ctrl", "h_MT_barrel_ctrl", 500, 0, 1000); h_MT_barrel_ctrl->Sumw2();
        TH1D* h_MT_endcap_ctrl = new TH1D("h_MT_endcap_ctrl", "h_MT_endcap_ctrl", 500, 0, 1000); h_MT_endcap_ctrl->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();

        TH1D* h_PFiso_barrel_nume_50to70   = new TH1D("h_PFiso_barrel_nume_50to70",   "h_PFiso_barrel_nume_50to70",   15, 0, 0.15); h_PFiso_barrel_nume_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_nume_50to70   = new TH1D("h_PFiso_endcap_nume_50to70",   "h_PFiso_endcap_nume_50to70",   15, 0, 0.15); h_PFiso_endcap_nume_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_deno_50to70   = new TH1D("h_PFiso_barrel_deno_50to70",   "h_PFiso_barrel_deno_50to70",   50, 0, 5);    h_PFiso_barrel_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_deno_50to70   = new TH1D("h_PFiso_endcap_deno_50to70",   "h_PFiso_endcap_deno_50to70",   50, 0, 5);    h_PFiso_endcap_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_50to70   = new TH1D("h_PFiso_barrel_ctrl_50to70",   "h_PFiso_barrel_ctrl_50to70",   50, 0.15, 5);    h_PFiso_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_50to70   = new TH1D("h_PFiso_endcap_ctrl_50to70",   "h_PFiso_endcap_ctrl_50to70",   50, 0.15, 5);    h_PFiso_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_nume_70to100  = new TH1D("h_PFiso_barrel_nume_70to100",  "h_PFiso_barrel_nume_70to100",  15, 0, 0.15); h_PFiso_barrel_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_nume_70to100  = new TH1D("h_PFiso_endcap_nume_70to100",  "h_PFiso_endcap_nume_70to100",  15, 0, 0.15); h_PFiso_endcap_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_deno_70to100  = new TH1D("h_PFiso_barrel_deno_70to100",  "h_PFiso_barrel_deno_70to100",  50, 0, 5);    h_PFiso_barrel_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_deno_70to100  = new TH1D("h_PFiso_endcap_deno_70to100",  "h_PFiso_endcap_deno_70to100",  50, 0, 5);    h_PFiso_endcap_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_70to100  = new TH1D("h_PFiso_barrel_ctrl_70to100",  "h_PFiso_barrel_ctrl_70to100",  50, 0.15, 5);    h_PFiso_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_70to100  = new TH1D("h_PFiso_endcap_ctrl_70to100",  "h_PFiso_endcap_ctrl_70to100",  50, 0.15, 5);    h_PFiso_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_nume_100to500 = new TH1D("h_PFiso_barrel_nume_100to500", "h_PFiso_barrel_nume_100to500", 15, 0, 0.15); h_PFiso_barrel_nume_100to500->Sumw2();
        TH1D* h_PFiso_endcap_nume_100to500 = new TH1D("h_PFiso_endcap_nume_100to500", "h_PFiso_endcap_nume_100to500", 15, 0, 0.15); h_PFiso_endcap_nume_100to500->Sumw2();
        TH1D* h_PFiso_barrel_deno_100to500 = new TH1D("h_PFiso_barrel_deno_100to500", "h_PFiso_barrel_deno_100to500", 50, 0, 5);    h_PFiso_barrel_deno_100to500->Sumw2();
        TH1D* h_PFiso_endcap_deno_100to500 = new TH1D("h_PFiso_endcap_deno_100to500", "h_PFiso_endcap_deno_100to500", 50, 0, 5);    h_PFiso_endcap_deno_100to500->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_100to500 = new TH1D("h_PFiso_barrel_ctrl_100to500", "h_PFiso_barrel_ctrl_100to500", 50, 0.15, 5);    h_PFiso_barrel_ctrl_100to500->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_100to500 = new TH1D("h_PFiso_endcap_ctrl_100to500", "h_PFiso_endcap_ctrl_100to500", 50, 0.15, 5);    h_PFiso_endcap_ctrl_100to500->Sumw2();

        TH1D* h_pT_barrel_nume_50to70   = new TH1D("h_pT_barrel_nume_50to70 ",  "h_pT_barrel_nume_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_50to70  ->Sumw2();
        TH1D* h_pT_endcap_nume_50to70   = new TH1D("h_pT_endcap_nume_50to70 ",  "h_pT_endcap_nume_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_50to70  ->Sumw2();
        TH1D* h_pT_barrel_deno_50to70   = new TH1D("h_pT_barrel_deno_50to70 ",  "h_pT_barrel_deno_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_50to70  ->Sumw2();
        TH1D* h_pT_endcap_deno_50to70   = new TH1D("h_pT_endcap_deno_50to70 ",  "h_pT_endcap_deno_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_50to70  ->Sumw2();
        TH1D* h_pT_barrel_ctrl_50to70   = new TH1D("h_pT_barrel_ctrl_50to70 ",  "h_pT_barrel_ctrl_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_endcap_ctrl_50to70   = new TH1D("h_pT_endcap_ctrl_50to70 ",  "h_pT_endcap_ctrl_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_barrel_nume_70to100  = new TH1D("h_pT_barrel_nume_70to100",  "h_pT_barrel_nume_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_70to100 ->Sumw2();
        TH1D* h_pT_endcap_nume_70to100  = new TH1D("h_pT_endcap_nume_70to100",  "h_pT_endcap_nume_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_70to100 ->Sumw2();
        TH1D* h_pT_barrel_deno_70to100  = new TH1D("h_pT_barrel_deno_70to100",  "h_pT_barrel_deno_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_70to100 ->Sumw2();
        TH1D* h_pT_endcap_deno_70to100  = new TH1D("h_pT_endcap_deno_70to100",  "h_pT_endcap_deno_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_70to100 ->Sumw2();
        TH1D* h_pT_barrel_ctrl_70to100  = new TH1D("h_pT_barrel_ctrl_70to100",  "h_pT_barrel_ctrl_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_endcap_ctrl_70to100  = new TH1D("h_pT_endcap_ctrl_70to100",  "h_pT_endcap_ctrl_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_barrel_nume_100to500 = new TH1D("h_pT_barrel_nume_100to500", "h_pT_barrel_nume_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_100to500->Sumw2();
        TH1D* h_pT_endcap_nume_100to500 = new TH1D("h_pT_endcap_nume_100to500", "h_pT_endcap_nume_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_100to500->Sumw2();
        TH1D* h_pT_barrel_deno_100to500 = new TH1D("h_pT_barrel_deno_100to500", "h_pT_barrel_deno_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_100to500->Sumw2();
        TH1D* h_pT_endcap_deno_100to500 = new TH1D("h_pT_endcap_deno_100to500", "h_pT_endcap_deno_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_100to500->Sumw2();
        TH1D* h_pT_barrel_ctrl_100to500 = new TH1D("h_pT_barrel_ctrl_100to500", "h_pT_barrel_ctrl_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_100to500->Sumw2();
        TH1D* h_pT_endcap_ctrl_100to500 = new TH1D("h_pT_endcap_ctrl_100to500", "h_pT_endcap_ctrl_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_100to500->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Double_t MET_pT, MET_phi, MET_sumEt;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
//        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("MET_sumEt", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("MET_sumEt", &MET_sumEt);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE){
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
                cout << "eta[1] = " << eta->at(0) << endl;
                cout << "phi[1] = " << phi->at(0) << endl;
            }

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            // -- Normalization -- //
            Double_t TotWeight = gen_weight;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl;

            if (Mgr.isMC == kTRUE && p_T->size() > 1) n2MC += TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight;
            if (Mgr.isMC == kFALSE && p_T->size() > 1) n2Data++;

            // For finding the leading muon
            TLorentzVector mu_lead;
            mu_lead.SetPtEtaPhiM(0, 0, 0, M_Mu);
            Double_t iso_lead = -9999;

            Double_t iso_count = 0;
            for (UInt_t i_mu=0; i_mu<p_T->size(); i_mu++)
            {
                if (p_T->at(i_mu) <= 52) continue;
                if (relPFiso->at(i_mu) < 0.15) iso_count++;

                // Selecting leading muon (could also try finding a muon with the best isolation)
                if (p_T->at(i_mu) > mu_lead.Pt())
                {
                    mu_lead.SetPtEtaPhiM(p_T->at(i_mu), eta->at(i_mu), phi->at(i_mu), M_Mu);
                    iso_lead = relPFiso->at(i_mu);
                }
            }
            if (iso_count > 1) continue;
            Double_t dTheta = mu_lead.Phi() - MET_phi;
            Double_t MT = sqrt(2 * mu_lead.Pt() * MET_pT * (1 - cos(dTheta)));
//            if (MT >= 60) continue;
//            if (MET_pT >= 50) continue;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            for (UInt_t i_mu=0; i_mu<p_T->size(); i_mu++)
            {
                if (p_T->at(i_mu) != p_T->at(i_mu))
                {
                    cout << p_T->at(i_mu) << " " << eta->at(i_mu) << " " << phi->at(i_mu) << " " << charge->at(i_mu) << " " << relPFiso->at(i_mu) << endl;
                    continue;
                }
                if (p_T->at(i_mu) <= 52) continue;
                if (DEBUG == kTRUE) cout << "i_mu = " << i_mu << endl;

                // -- Efficiency scale factor -- //
//                if(Mgr.isMC == kTRUE)
//                {
//                    weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF_new(MuMu);
//                    weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH_new(MuMu);
//                    effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
//                }

                if (relPFiso->at(i_mu) < 0.15) // Signal/Numerator
                {
                    h_eta_nume->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                    {
                        h_pT_barrel_nume->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_nume->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_barrel_nume->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_barrel_nume_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_barrel_nume_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_barrel_nume_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    else // Endcap
                    {
                        h_pT_endcap_nume->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_nume->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_endcap_nume->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap_nume_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap_nume_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap_nume_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                } // End of if(Signal/Numerator)
                else // Control
                {
                    h_eta_ctrl->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                    {
                        h_pT_barrel_ctrl->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_ctrl->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_barrel_ctrl->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_barrel_ctrl_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_barrel_ctrl_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_barrel_ctrl_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    else // Endcap
                    {
                        h_pT_endcap_ctrl->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_ctrl->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_endcap_ctrl->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap_ctrl_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap_ctrl_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap_ctrl_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                }// End of if(Control)
                // Denominator
                h_eta_deno->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                {
                    h_pT_barrel_deno->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_PFiso_barrel_deno->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_TRKiso_barrel_deno->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (p_T->at(i_mu) < 70)
                    {
                        h_pT_barrel_deno_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (p_T->at(i_mu) < 100)
                    {
                        h_pT_barrel_deno_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else //if (p_T->at(i_mu) < 500)
                    {
                        h_pT_barrel_deno_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                else // Endcap
                {
                    h_pT_endcap_deno->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_PFiso_endcap_deno->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_TRKiso_endcap_deno->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (p_T->at(i_mu) < 70)
                    {
                        h_pT_endcap_deno_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (p_T->at(i_mu) < 100)
                    {
                        h_pT_endcap_deno_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else //if (p_T->at(i_mu) < 500)
                    {
                        h_pT_endcap_deno_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
            }// End of i_mu iteration

            if (fabs(mu_lead.Eta()) < 1.2)
            {
                h_MT_barrel_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (iso_lead < 0.15) h_MT_barrel_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                else h_MT_barrel_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_MT_endcap_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (iso_lead < 0.15) h_MT_endcap_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                else h_MT_endcap_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }

            bar.Draw(i);

        }// End of event iteration

        if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);

        f->cd();
        cout << "\tWriting into file...";

        h_pT_barrel_nume->Write();
        h_pT_endcap_nume->Write();
        h_pT_barrel_deno->Write();
        h_pT_endcap_deno->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_eta_nume->Write();
        h_eta_deno->Write();
        h_eta_ctrl->Write();
        h_PFiso_barrel_nume->Write();
        h_PFiso_endcap_nume->Write();
        h_PFiso_barrel_deno->Write();
        h_PFiso_endcap_deno->Write();
        h_PFiso_barrel_ctrl->Write();
        h_PFiso_endcap_ctrl->Write();
        h_TRKiso_barrel_nume->Write();
        h_TRKiso_endcap_nume->Write();
        h_TRKiso_barrel_deno->Write();
        h_TRKiso_endcap_deno->Write();
        h_TRKiso_barrel_ctrl->Write();
        h_TRKiso_endcap_ctrl->Write();
        h_MET->Write();
        h_MT_barrel_nume->Write();
        h_MT_endcap_nume->Write();
        h_MT_barrel_deno->Write();
        h_MT_endcap_deno->Write();
        h_MT_barrel_ctrl->Write();
        h_MT_endcap_ctrl->Write();
        h_nVTX->Write();

        h_pT_barrel_nume_50to70->Write();
        h_pT_endcap_nume_50to70->Write();
        h_pT_barrel_deno_50to70->Write();
        h_pT_endcap_deno_50to70->Write();
        h_pT_barrel_ctrl_50to70->Write();
        h_pT_endcap_ctrl_50to70->Write();
        h_pT_barrel_nume_70to100->Write();
        h_pT_endcap_nume_70to100->Write();
        h_pT_barrel_deno_70to100->Write();
        h_pT_endcap_deno_70to100->Write();
        h_pT_barrel_ctrl_70to100->Write();
        h_pT_endcap_ctrl_70to100->Write();
        h_pT_barrel_nume_100to500->Write();
        h_pT_endcap_nume_100to500->Write();
        h_pT_barrel_deno_100to500->Write();
        h_pT_endcap_deno_100to500->Write();
        h_pT_barrel_ctrl_100to500->Write();
        h_pT_endcap_ctrl_100to500->Write();

        h_PFiso_barrel_nume_50to70->Write();
        h_PFiso_endcap_nume_50to70->Write();
        h_PFiso_barrel_deno_50to70->Write();
        h_PFiso_endcap_deno_50to70->Write();
        h_PFiso_barrel_ctrl_50to70->Write();
        h_PFiso_endcap_ctrl_50to70->Write();
        h_PFiso_barrel_nume_70to100->Write();
        h_PFiso_endcap_nume_70to100->Write();
        h_PFiso_barrel_deno_70to100->Write();
        h_PFiso_endcap_deno_70to100->Write();
        h_PFiso_barrel_ctrl_70to100->Write();
        h_PFiso_endcap_ctrl_70to100->Write();
        h_PFiso_barrel_nume_100to500->Write();
        h_PFiso_endcap_nume_100to500->Write();
        h_PFiso_barrel_deno_100to500->Write();
        h_PFiso_endcap_deno_100to500->Write();
        h_PFiso_barrel_ctrl_100to500->Write();
        h_PFiso_endcap_ctrl_100to500->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Total weighted events with 2 or more denominator muons in MC: " << n2MC << endl;
    cout << "Total events with 2 or more denominator muons in Data: " << n2Data << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_FR_HistMaker()


/// -------------------------------- Electron Channel ------------------------------------ ///
void E_QCD_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"QCDest_E"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.01);

    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX_"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_temp = new TH1D("h_mass_template_"+Mgr.Procname[pr], "h_mass_template_"+Mgr.Procname[pr], binnum, massbins); h_mass_temp->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<double> *chIso03 = new std::vector<double>;
        std::vector<double> *nhIso03 = new std::vector<double>;
        std::vector<double> *phIso03 = new std::vector<double>;
        std::vector<double> *ChIso03FromPU = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<double> *relPFiso_dBeta = new std::vector<double>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("etaSC", 1);
        chain->SetBranchStatus("charge", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("InvEminusInvP", 1);
        chain->SetBranchStatus("chIso03", 1);
        chain->SetBranchStatus("nhIso03", 1);
        chain->SetBranchStatus("phIso03", 1);
        chain->SetBranchStatus("ChIso03FromPU", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("relPFiso_dBeta", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("etaSC", &etaSC);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("chIso03", &chIso03);
        chain->SetBranchAddress("nhIso03", &nhIso03);
        chain->SetBranchAddress("phIso03", &phIso03);
        chain->SetBranchAddress("ChIso03FromPU", &ChIso03FromPU);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_dBeta", &relPFiso_dBeta);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassEle=0, nFailEle=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            for (UInt_t e=0; e<passMediumID->size(); e++)
            {
                if (passMediumID->at(e) == 1) nPassEle++;
                else nFailEle++;
            }
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
            if ((etaSC->at(0) > 1.4442 && etaSC->at(0) < 1.566) || (etaSC->at(1) > 1.4442 && etaSC->at(1) < 1.566)) continue;
            if (etaSC->at(0) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                          fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (etaSC->at(1) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                          fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (etaSC->at(0) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (etaSC->at(1) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;

            if (passMediumID->at(0) == 1 || passMediumID->at(1)  == 1) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nElectrons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele1, ele2, ele1_SF, ele2_SF;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            ele1_SF.SetPtEtaPhiM(p_T->at(0), etaSC->at(0), phi->at(0), M_Elec);
            ele2_SF.SetPtEtaPhiM(p_T->at(1), etaSC->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;
            if (Mgr.isMC == kTRUE)
            {
                effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, -1);
            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR1, FR2;
            FR1 = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
            FR2 = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (MET_pT > 75)
                h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
//            if (charge->at(0) == charge->at(1))
//                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }
        cout << "\t # passed electrons: " << nPassEle << endl;
        cout << "\t # failed electrons: " << nFailEle << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_SS->Write();
        h_mass_temp->Write();
        h_nVTX->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_E"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_E"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_QCD_HistMaker()


void E_WJET_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"WJETest_E"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For W+Jets estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.5);

    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_MET = new TH1D("h_MET_"+Mgr.Procname[pr], "h_MET"+Mgr.Procname[pr], 100, 0, 1000); h_MET->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_temp = new TH1D("h_mass_template_"+Mgr.Procname[pr], "h_mass_template_"+Mgr.Procname[pr], binnum, massbins); h_mass_temp->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<double> *chIso03 = new std::vector<double>;
        std::vector<double> *nhIso03 = new std::vector<double>;
        std::vector<double> *phIso03 = new std::vector<double>;
        std::vector<double> *ChIso03FromPU = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<double> *relPFiso_dBeta = new std::vector<double>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("etaSC", 1);
        chain->SetBranchStatus("charge", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("InvEminusInvP", 1);
        chain->SetBranchStatus("chIso03", 1);
        chain->SetBranchStatus("nhIso03", 1);
        chain->SetBranchStatus("phIso03", 1);
        chain->SetBranchStatus("ChIso03FromPU", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("relPFiso_dBeta", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("etaSC", &etaSC);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("chIso03", &chIso03);
        chain->SetBranchAddress("nhIso03", &nhIso03);
        chain->SetBranchAddress("phIso03", &phIso03);
        chain->SetBranchAddress("ChIso03FromPU", &ChIso03FromPU);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_dBeta", &relPFiso_dBeta);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // W+Jets selection
            if (p_T->size() != 2) continue;
            if (passMediumID->at(0) == 1 && passMediumID->at(1) == 1) continue;
            if (passMediumID->at(0) == 0 && passMediumID->at(1) == 0) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if ((etaSC->at(0) > 1.4442 && etaSC->at(0) < 1.566) || (etaSC->at(1) > 1.4442 && etaSC->at(1) < 1.566)) continue;
            if (etaSC->at(0) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                          fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (etaSC->at(1) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                          fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (etaSC->at(0) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (etaSC->at(1) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele1, ele2, ele1_SF, ele2_SF;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            ele1_SF.SetPtEtaPhiM(p_T->at(0), etaSC->at(0), phi->at(0), M_Elec);
            ele2_SF.SetPtEtaPhiM(p_T->at(1), etaSC->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR (and SF) WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR;
            if (passMediumID->at(0) == 0 && passMediumID->at(1) == 1) // First fails, second passes
            {
                FR = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
                if (Mgr.isMC == kTRUE)
                {
                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 2);
                }

            }
            else // Second fails, first passes
            {
                FR = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
                if (Mgr.isMC == kTRUE)
                {
                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 1);
                }
            }
            FRweight = FR / (1 - FR);
            if (DEBUG == kTRUE) cout << "FR = " << FR << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- Histogram filling -- //
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (MET_pT > 75)
                h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
//            if (charge->at(0) == charge->at(1))
//                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();        
        h_mass_SS->Write();
        h_mass_temp->Write();
        h_MET->Write();
        h_nVTX->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_E"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_E"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_WJET_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_QCD_HistMaker (Bool_t DEBUG, Int_t type=1)
// type=0 uses files with Mu50 trigger
// type=1 -- no trigger
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"QCDest_Mu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", "sigCtrl_template");

    TH1D *h_mass_test[_EndOf_SingleMuon_Normal];
    TH1D *h_mass_test_SS[_EndOf_SingleMuon_Normal];

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.4);

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_forFit = new TH1D("h_mass_forFit_"+Mgr.Procname[pr], "h_mass_forFit_"+Mgr.Procname[pr], 37, 15, 200); h_mass_forFit->Sumw2();
        TH1D* h_mass_SS_forFit = new TH1D("h_mass_SS_forFit_"+Mgr.Procname[pr], "h_mass_SS_forFit_"+Mgr.Procname[pr], 37, 15, 200); h_mass_SS_forFit->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX_"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_pT_lead = new TH1D("h_pT_lead_"+Mgr.Procname[pr], "h_pT_lead_"+Mgr.Procname[pr], 100, 0, 1000); h_pT_lead->Sumw2();
        TH1D* h_pT_sublead = new TH1D("h_pT_sublead_"+Mgr.Procname[pr], "h_pT_sublead_"+Mgr.Procname[pr], 100, 0, 1000); h_pT_sublead->Sumw2();
        TH2D* h2_pT = new TH2D("h2_pT_"+Mgr.Procname[pr], "h2_pT_"+Mgr.Procname[pr], 49, 10, 500, 49, 10, 500);
        h_mass_test[pr] = new TH1D("h_test_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test[pr]->Sumw2();
        h_mass_test[pr]->SetDirectory(0);
        h_mass_test_SS[pr] = new TH1D("h_test_SS_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test_SS[pr]->Sumw2();
        h_mass_test_SS[pr]->SetDirectory(0);

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        if (type == 0) // Mu50 files
        {
            chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 1) // Triggerless files
        {
            chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else
        {
            cout << "Wrong type! Select 0 (Mu50) or 1 (Triggerless)" << endl;
            return;
        }
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassMu=0, nFailMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            for (UInt_t m=0; m<relPFiso->size(); m++)
            {
                if (relPFiso->at(m) < 0.15) nPassMu++;
                else nFailMu++;
            }
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
//            if (charge->at(0) == charge->at(1)) continue;
            if (relPFiso->at(0) < 0.15 || relPFiso->at(1) < 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
//            if (p_T->at(0) < 2 || p_T->at(1) < 2) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;
            if (Mgr.isMC == kTRUE)
            {
                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(mu1, mu2);
                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(mu1, mu2);
                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR1, FR2;
            FR1 = analyzer->FakeRate(p_T->at(0), eta->at(0));
            FR2 = analyzer->FakeRate(p_T->at(1), eta->at(1));
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (charge->at(0) != charge->at(1))
            {
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_SS_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test_SS[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            if (mu1.Pt() > mu2.Pt())
            {
                h2_pT->Fill(mu1.Pt(), mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_lead->Fill(mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_sublead->Fill(mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }
            else
            {
                h2_pT->Fill(mu2.Pt(), mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_lead->Fill(mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_sublead->Fill(mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }
        cout << "\t # passed muons: " << nPassMu << endl;
        cout << "\t # failed muons: " << nFailMu << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_forFit->Write();
        h_mass_SS->Write();
        h_mass_SS_forFit->Write();
        h_nVTX->Write();
        h_pT_lead->Write();
        h_pT_sublead->Write();
        h2_pT->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_data = new TCanvas("data","data", 800, 800);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_C]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_D]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_E]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_F]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_G]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_H]);
    h_mass_test[_SingleMuon_B]->Draw("hist");
    c_data->SetLogx();
    c_data->Update();

    TCanvas *c_data_SS = new TCanvas("data_SS","data SS", 800, 800);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_C]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_D]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_E]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_F]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_G]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_H]);
    h_mass_test_SS[_SingleMuon_B]->Draw("hist");
    c_data_SS->SetLogx();
    c_data_SS->Update();

    TCanvas *c_DY = new TCanvas("dy","dy", 800, 800);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_50to100]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_100to200]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_200to400]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_400to500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_500to700]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_700to800]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_800to1000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1000to1500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1500to2000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_2000to3000]);
    h_mass_test[_DY_10to50]->Draw("hist");
    c_DY->SetLogx();
    c_DY->Update();

    TCanvas *c_DY_SS = new TCanvas("dy_SS","dy SS", 800, 800);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_50to100]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_100to200]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_200to400]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_400to500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_500to700]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_700to800]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_800to1000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1000to1500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1500to2000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_2000to3000]);
    h_mass_test_SS[_DY_10to50]->Draw("hist");
    c_DY_SS->SetLogx();
    c_DY_SS->Update();

    TCanvas *c_tt = new TCanvas("ttbar","ttbar", 800, 800);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_700to1000]);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_1000toInf]);
    h_mass_test[_ttbar]->Draw("hist");
    c_tt->SetLogx();
    c_tt->Update();

    TCanvas *c_tt_SS = new TCanvas("ttbar_SS","ttbar SS", 800, 800);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_700to1000]);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_1000toInf]);
    h_mass_test_SS[_ttbar]->Draw("hist");
    c_tt_SS->SetLogx();
    c_tt_SS->Update();

    TCanvas *c_wjets = new TCanvas("wjets","wjets", 800, 800);
    h_mass_test[_WJets]->Add(h_mass_test[_WJets_ext2v5]);
    h_mass_test[_WJets]->Draw("hist");
    c_wjets->SetLogx();
    c_wjets->Update();

    TCanvas *c_wjets_SS = new TCanvas("wjets_SS","wjets SS", 800, 800);
    h_mass_test_SS[_WJets]->Add(h_mass_test_SS[_WJets_ext2v5]);
    h_mass_test_SS[_WJets]->Draw("hist");
    c_wjets_SS->SetLogx();
    c_wjets_SS->Update();

    TH1D *h_wjets_SS_scaled = (TH1D*)h_mass_test_SS[_WJets]->Clone("h_wjets_SS_scaled");
    h_wjets_SS_scaled->Rebin(8);
    h_wjets_SS_scaled->Scale(3);
    h_wjets_SS_scaled->SetLineColor(kRed);
    h_wjets_SS_scaled->SetDirectory(0);
    TH1D *h_wjets_OS_compare = (TH1D*)h_mass_test[_WJets]->Clone("h_wjets_OS_compare");
    h_wjets_OS_compare->SetDirectory(0);
    h_wjets_OS_compare->Rebin(8);
    TCanvas *c_wjets_compare = new TCanvas("wjets_compare", "WJets OS vs SSx3", 800, 800);
    h_wjets_OS_compare->Draw("hist");
    h_wjets_SS_scaled->Draw("histsame");
    c_wjets_compare->SetLogx();
    c_wjets_compare->Update();
    cout << "Ratio of W+Jets opposite-sign and same-sign integrals: " << h_mass_test[_WJets]->Integral()/h_mass_test_SS[_WJets]->Integral() << endl;

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_Mu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_Mu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_QCD_HistMaker()


void Mu_WJET_HistMaker (Bool_t DEBUG, Int_t type=2)
// type=0 uses files with Mu50 trigger
// type=1 -- IsoMu24_OR_IsoTkMu24
// type=2 -- no trigger
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"WJETest_Mu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For W+Jets estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", "sigCtrl_template");

    TH1D *h_mass_test[_EndOf_SingleMuon_Normal];
    TH1D *h_mass_test_SS[_EndOf_SingleMuon_Normal];

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_forFit = new TH1D("h_mass_forFit_"+Mgr.Procname[pr], "h_mass_forFit_"+Mgr.Procname[pr], 37,15,200); h_mass_forFit->Sumw2();
        TH1D* h_mass_SS_forFit = new TH1D("h_mass_SS_forFit_"+Mgr.Procname[pr], "h_mass_SS_forFit_"+Mgr.Procname[pr], 37,15,200); h_mass_SS_forFit->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_eta_diff = new TH1D("h_eta_diff_"+Mgr.Procname[pr], "h_eta_diff_"+Mgr.Procname[pr], 50, 0, 5); h_eta_diff->Sumw2();
        TH1D* h_phi_diff = new TH1D("h_phi_diff_"+Mgr.Procname[pr], "h_phi_diff_"+Mgr.Procname[pr], 63, -3.15, 3.15); h_phi_diff->Sumw2();
        h_mass_test[pr] = new TH1D("h_test_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test[pr]->Sumw2();
        h_mass_test[pr]->SetDirectory(0);
        h_mass_test_SS[pr] = new TH1D("h_test_SS_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test_SS[pr]->Sumw2();
        h_mass_test_SS[pr]->SetDirectory(0);

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        if (type == 0) // Mu50 files
        {
            chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 1) // IsoMu24 files
        {
            chain->Add(Dir+"SelectedForWJETest_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForWJETest_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 2) // Triggerless files
        {
            chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else
        {
            cout << "Wrong type! Select 0 (Mu50), 1 (IsoMu24) or 2 (Triggerless)" << endl;
            return;
        }
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // W+Jets selection
            if (p_T->size() != 2) continue;
//            if (charge->at(0) == charge->at(1)) continue;
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1) < 0.15) continue;
            if (relPFiso->at(0) > 0.15 && relPFiso->at(1) > 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
//            if (p_T->at(0) < 2 || p_T->at(1) < 2) continue;

            nPass++;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);

            Double_t mass = (mu1+mu2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;
            if (Mgr.isMC == kTRUE)
            {
                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(mu1, mu2);
                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(mu1, mu2);
                if (weight1<0 || weight2<0)
                    cout << "weights: " << weight1 << " " << weight2 << "   Pt: " << p_T->at(0) << " " << mu1.Pt() << " " << p_T->at(1)
                         << " " << mu2.Pt() << "   Eta: " << eta->at(0) << " " << mu1.Eta() << " " << eta->at(1) << " " << mu2.Eta() << endl;
                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            if (relPFiso->at(0) >= 0.15 && relPFiso->at(1) >= 0.15) // Both failing -- skip
            {
                cout << "Both failed" << endl; // Comment this when not needed
                continue;
            }
            else // Only one failing
            {
                Double_t FR;
                if (relPFiso->at(0) > relPFiso->at(1)) // First fails, second passes
                    FR = analyzer->FakeRate(p_T->at(0), eta->at(0));
                else // Second fails, first passes
                    FR = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FRweight = FR / (1 - FR);
                if (DEBUG == kTRUE) cout << "FR = " << FR << "   FRweight = " << FRweight << endl;
                avgFRweight += FRweight;
            }

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- Histogram filling -- //
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_eta_diff->Fill(fabs(eta->at(0)-eta->at(1)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_phi_diff->Fill(phi->at(0)-phi->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (charge->at(0) != charge->at(1))
            {
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_SS_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test_SS[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }


        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_forFit->Write();
        h_mass_SS->Write();
        h_mass_SS_forFit->Write();
        h_nVTX->Write();
        h_eta_diff->Write();
        h_phi_diff->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_data = new TCanvas("data","data", 800, 800);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_C]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_D]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_E]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_F]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_G]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_H]);
    h_mass_test[_SingleMuon_B]->Draw("hist");
    c_data->SetLogx();
    c_data->Update();

    TCanvas *c_data_SS = new TCanvas("data_SS","data SS", 800, 800);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_C]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_D]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_E]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_F]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_G]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_H]);
    h_mass_test_SS[_SingleMuon_B]->Draw("hist");
    c_data_SS->SetLogx();
    c_data_SS->Update();

    TCanvas *c_DY = new TCanvas("dy","dy", 800, 800);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_50to100]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_100to200]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_200to400]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_400to500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_500to700]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_700to800]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_800to1000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1000to1500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1500to2000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_2000to3000]);
    h_mass_test[_DY_10to50]->Draw("hist");
    c_DY->SetLogx();
    c_DY->Update();

    TCanvas *c_DY_SS = new TCanvas("dy_SS","dy SS", 800, 800);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_50to100]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_100to200]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_200to400]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_400to500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_500to700]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_700to800]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_800to1000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1000to1500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1500to2000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_2000to3000]);
    h_mass_test_SS[_DY_10to50]->Draw("hist");
    c_DY_SS->SetLogx();
    c_DY_SS->Update();

    TCanvas *c_tt = new TCanvas("ttbar","ttbar", 800, 800);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_700to1000]);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_1000toInf]);
    h_mass_test[_ttbar]->Draw("hist");
    c_tt->SetLogx();
    c_tt->Update();

    TCanvas *c_tt_SS = new TCanvas("ttbar_SS","ttbar SS", 800, 800);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_700to1000]);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_1000toInf]);
    h_mass_test_SS[_ttbar]->Draw("hist");
    c_tt_SS->SetLogx();
    c_tt_SS->Update();

    TCanvas *c_wjets = new TCanvas("wjets","wjets", 800, 800);
    h_mass_test[_WJets]->Add(h_mass_test[_WJets_ext2v5]);
    h_mass_test[_WJets]->Draw("hist");
    c_wjets->SetLogx();
    c_wjets->Update();

    TCanvas *c_wjets_SS = new TCanvas("wjets_SS","wjets SS", 800, 800);
    h_mass_test_SS[_WJets]->Add(h_mass_test_SS[_WJets_ext2v5]);
    h_mass_test_SS[_WJets]->Draw("hist");
    c_wjets_SS->SetLogx();
    c_wjets_SS->Update();

    TH1D *h_wjets_SS_scaled = (TH1D*)h_mass_test_SS[_WJets]->Clone("h_wjets_SS_scaled");
    h_wjets_SS_scaled->Rebin(8);
    h_wjets_SS_scaled->Scale(3);
    h_wjets_SS_scaled->SetLineColor(kRed);
    h_wjets_SS_scaled->SetDirectory(0);
    TH1D *h_wjets_OS_compare = (TH1D*)h_mass_test[_WJets]->Clone("h_wjets_OS_compare");
    h_wjets_OS_compare->SetDirectory(0);
    h_wjets_OS_compare->Rebin(8);
    TCanvas *c_wjets_compare = new TCanvas("wjets_compare", "WJets OS vs SSx3", 800, 800);
    h_wjets_OS_compare->Draw("hist");
    h_wjets_SS_scaled->Draw("histsame");
    c_wjets_compare->SetLogx();
    c_wjets_compare->Update();
    cout << "Ratio of W+Jets opposite-sign and same-sign integrals: " << h_mass_test[_WJets]->Integral()/h_mass_test_SS[_WJets]->Integral() << endl;

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_Mu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_Mu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_WJET_HistMaker()


/// -------------------------------- Electron Channel ------------------------------------ ///
void EMu_QCD_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/EMu/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"QCDest_EMu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("IsoMu24_OR_IsoTkMu24");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele("/media/sf_DATA/FR/Electron/FakeRate_electron.root", "subtract");
    analyzer->SetupFRvalues("/media/sf_DATA/FR/Muon/FakeRate_muon.root", "sigCtrl_template");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.01);

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX_"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();

        std::vector<double> *e_p_T = new std::vector<double>;
        std::vector<double> *e_eta = new std::vector<double>;
        std::vector<double> *e_etaSC = new std::vector<double>;
        std::vector<double> *e_phi = new std::vector<double>;
        std::vector<int> *e_charge = new std::vector<int>;
        std::vector<double> *e_Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *e_dEtaInSeed = new std::vector<double>;
        std::vector<double> *e_dPhiIn = new std::vector<double>;
        std::vector<double> *e_HoverE = new std::vector<double>;
        std::vector<double> *e_InvEminusInvP = new std::vector<double>;
        std::vector<double> *e_chIso03 = new std::vector<double>;
        std::vector<double> *e_nhIso03 = new std::vector<double>;
        std::vector<double> *e_phIso03 = new std::vector<double>;
        std::vector<double> *e_ChIso03FromPU = new std::vector<double>;
        std::vector<int> *e_mHits = new std::vector<int>;
        std::vector<int> *e_passConvVeto = new std::vector<int>;
        std::vector<double> *e_relPFiso_dBeta = new std::vector<double>;
        std::vector<double> *e_relPFiso_Rho = new std::vector<double>;
        std::vector<int> *e_passMediumID = new std::vector<int>;
        std::vector<double> *mu_p_T = new std::vector<double>;
        std::vector<double> *mu_eta = new std::vector<double>;
        std::vector<double> *mu_phi = new std::vector<double>;
        std::vector<int> *mu_charge = new std::vector<int>;
        std::vector<double> *mu_relPFiso_dBeta = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("e_p_T", 1);
        chain->SetBranchStatus("e_eta", 1);
        chain->SetBranchStatus("e_etaSC", 1);
        chain->SetBranchStatus("e_phi", 1);
        chain->SetBranchStatus("e_charge", 1);
        chain->SetBranchStatus("e_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("e_dEtaInSeed", 1);
        chain->SetBranchStatus("e_dPhiIn", 1);
        chain->SetBranchStatus("e_HoverE", 1);
        chain->SetBranchStatus("e_InvEminusInvP", 1);
        chain->SetBranchStatus("e_chIso03", 1);
        chain->SetBranchStatus("e_nhIso03", 1);
        chain->SetBranchStatus("e_phIso03", 1);
        chain->SetBranchStatus("e_ChIso03FromPU", 1);
        chain->SetBranchStatus("e_mHits", 1);
        chain->SetBranchStatus("e_passConvVeto", 1);
        chain->SetBranchStatus("e_relPFiso_dBeta", 1);
        chain->SetBranchStatus("e_relPFiso_Rho", 1);
        chain->SetBranchStatus("e_passMediumID", 1);
        chain->SetBranchStatus("mu_p_T", 1);
        chain->SetBranchStatus("mu_eta", 1);
        chain->SetBranchStatus("mu_phi", 1);
        chain->SetBranchStatus("mu_charge", 1);
        chain->SetBranchStatus("mu_relPFiso_dBeta", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);

        chain->SetBranchAddress("e_p_T", &e_p_T);
        chain->SetBranchAddress("e_eta", &e_eta);
        chain->SetBranchAddress("e_etaSC", &e_etaSC);
        chain->SetBranchAddress("e_phi", &e_phi);
        chain->SetBranchAddress("e_charge", &e_charge);
        chain->SetBranchAddress("e_Full5x5_SigmaIEtaIEta", &e_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("e_dEtaInSeed", &e_dEtaInSeed);
        chain->SetBranchAddress("e_dPhiIn", &e_dPhiIn);
        chain->SetBranchAddress("e_HoverE", &e_HoverE);
        chain->SetBranchAddress("e_InvEminusInvP", &e_InvEminusInvP);
        chain->SetBranchAddress("e_chIso03", &e_chIso03);
        chain->SetBranchAddress("e_nhIso03", &e_nhIso03);
        chain->SetBranchAddress("e_phIso03", &e_phIso03);
        chain->SetBranchAddress("e_ChIso03FromPU", &e_ChIso03FromPU);
        chain->SetBranchAddress("e_mHits", &e_mHits);
        chain->SetBranchAddress("e_passConvVeto", &e_passConvVeto);
        chain->SetBranchAddress("e_relPFiso_dBeta", &e_relPFiso_dBeta);
        chain->SetBranchAddress("e_relPFiso_Rho", &e_relPFiso_Rho);
        chain->SetBranchAddress("e_passMediumID", &e_passMediumID);
        chain->SetBranchAddress("mu_p_T", &mu_p_T);
        chain->SetBranchAddress("mu_eta", &mu_eta);
        chain->SetBranchAddress("mu_phi", &mu_phi);
        chain->SetBranchAddress("mu_charge", &mu_charge);
        chain->SetBranchAddress("mu_relPFiso_dBeta", &mu_relPFiso_dBeta);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassEle=0, nFailEle=0, nPassMu=0, nFailMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            if (e_passMediumID->at(0)) nPassEle++;
            else nFailEle++;
            if (mu_relPFiso_dBeta->at(0) < 0.15) nPassMu++;
            else nFailMu++;

            // QCD selection
            if (e_passMediumID->at(0) || mu_relPFiso_dBeta->at(0)  < 0.15) continue;
            if (e_p_T->at(0) < 17 || mu_p_T->at(0) < 17) continue;
            if (e_p_T->at(0) < 28 && mu_p_T->at(0) < 28) continue;

            if (e_p_T->at(0) != e_p_T->at(0)) cout << e_p_T->at(0) << " " << e_eta->at(0) << " " << e_phi->at(0) << " " << e_charge->at(0) << " " << e_relPFiso_Rho->at(0) << endl;
            if (mu_p_T->at(0) != mu_p_T->at(0)) cout << mu_p_T->at(0) << " " << mu_eta->at(0) << " " << mu_phi->at(0) << " " << mu_charge->at(0) << " " << mu_relPFiso_dBeta->at(0) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T->at(0);
                cout << "\teta = " << e_eta->at(0);
                cout << "\tphi = " << e_phi->at(0) << endl;
                cout << "\tpassMediumID = " << e_passMediumID->at(0) << endl;
                cout << "Muon p_T = " << mu_p_T->at(0);
                cout << "\teta = " << mu_eta->at(0);
                cout << "\tphi = " << mu_phi->at(0) << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta->at(0) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu, ele_SF;
            ele.SetPtEtaPhiM(e_p_T->at(0), e_eta->at(0), e_phi->at(0), M_Elec);
            mu.SetPtEtaPhiM(mu_p_T->at(0), mu_eta->at(0), mu_phi->at(0), M_Mu);
            ele_SF.SetPtEtaPhiM(e_p_T->at(0), e_etaSC->at(0), e_phi->at(0), M_Elec);
            Double_t mass = (ele+mu).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;
//            if (Mgr.isMC == kTRUE)
//            {
//                effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, -1);
//            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR1, FR2;
            FR1 = analyzer->FakeRate_ele(e_p_T->at(0), e_etaSC->at(0));
            FR2 = analyzer->FakeRate(mu_p_T->at(0), mu_eta->at(0));
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (e_charge->at(0) != mu_charge->at(0))
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            else
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }
        cout << "\t # passed electrons: " << nPassEle << endl;
        cout << "\t # failed electrons: " << nFailEle << endl;
        cout << "\t # passed muons: " << nPassMu << endl;
        cout << "\t # failed muons: " << nFailMu << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_nVTX->Write();
        h_mass_SS->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDEMEnriched_Normal; // next -- GJets_20to100
        if (pr == _GJets_2000to5000) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_EMu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_EMu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_QCD_HistMaker()


void EMu_WJET_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/EMu/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"WJETest_EMu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("IsoMu24_OR_IsoTkMu24");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For W+Jets estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele("/media/sf_DATA/FR/Electron/FakeRate_electron.root", "subtract");
    analyzer->SetupFRvalues("/media/sf_DATA/FR/Muon/FakeRate_muon.root", "sigCtrl_template");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.5);

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_temp = new TH1D("h_mass_template_"+Mgr.Procname[pr], "h_mass_template_"+Mgr.Procname[pr], binnum, massbins); h_mass_temp->Sumw2();

        std::vector<double> *e_p_T = new std::vector<double>;
        std::vector<double> *e_eta = new std::vector<double>;
        std::vector<double> *e_etaSC = new std::vector<double>;
        std::vector<double> *e_phi = new std::vector<double>;
        std::vector<int> *e_charge = new std::vector<int>;
        std::vector<double> *e_Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *e_dEtaInSeed = new std::vector<double>;
        std::vector<double> *e_dPhiIn = new std::vector<double>;
        std::vector<double> *e_HoverE = new std::vector<double>;
        std::vector<double> *e_InvEminusInvP = new std::vector<double>;
        std::vector<double> *e_chIso03 = new std::vector<double>;
        std::vector<double> *e_nhIso03 = new std::vector<double>;
        std::vector<double> *e_phIso03 = new std::vector<double>;
        std::vector<double> *e_ChIso03FromPU = new std::vector<double>;
        std::vector<int> *e_mHits = new std::vector<int>;
        std::vector<int> *e_passConvVeto = new std::vector<int>;
        std::vector<double> *e_relPFiso_dBeta = new std::vector<double>;
        std::vector<double> *e_relPFiso_Rho = new std::vector<double>;
        std::vector<int> *e_passMediumID = new std::vector<int>;
        std::vector<double> *mu_p_T = new std::vector<double>;
        std::vector<double> *mu_eta = new std::vector<double>;
        std::vector<double> *mu_phi = new std::vector<double>;
        std::vector<int> *mu_charge = new std::vector<int>;
        std::vector<double> *mu_relPFiso_dBeta = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("e_p_T", 1);
        chain->SetBranchStatus("e_eta", 1);
        chain->SetBranchStatus("e_etaSC", 1);
        chain->SetBranchStatus("e_phi", 1);
        chain->SetBranchStatus("e_charge", 1);
        chain->SetBranchStatus("e_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("e_dEtaInSeed", 1);
        chain->SetBranchStatus("e_dPhiIn", 1);
        chain->SetBranchStatus("e_HoverE", 1);
        chain->SetBranchStatus("e_InvEminusInvP", 1);
        chain->SetBranchStatus("e_chIso03", 1);
        chain->SetBranchStatus("e_nhIso03", 1);
        chain->SetBranchStatus("e_phIso03", 1);
        chain->SetBranchStatus("e_ChIso03FromPU", 1);
        chain->SetBranchStatus("e_mHits", 1);
        chain->SetBranchStatus("e_passConvVeto", 1);
        chain->SetBranchStatus("e_relPFiso_dBeta", 1);
        chain->SetBranchStatus("e_relPFiso_Rho", 1);
        chain->SetBranchStatus("e_passMediumID", 1);
        chain->SetBranchStatus("mu_p_T", 1);
        chain->SetBranchStatus("mu_eta", 1);
        chain->SetBranchStatus("mu_phi", 1);
        chain->SetBranchStatus("mu_charge", 1);
        chain->SetBranchStatus("mu_relPFiso_dBeta", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);

        chain->SetBranchAddress("e_p_T", &e_p_T);
        chain->SetBranchAddress("e_eta", &e_eta);
        chain->SetBranchAddress("e_etaSC", &e_etaSC);
        chain->SetBranchAddress("e_phi", &e_phi);
        chain->SetBranchAddress("e_charge", &e_charge);
        chain->SetBranchAddress("e_Full5x5_SigmaIEtaIEta", &e_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("e_dEtaInSeed", &e_dEtaInSeed);
        chain->SetBranchAddress("e_dPhiIn", &e_dPhiIn);
        chain->SetBranchAddress("e_HoverE", &e_HoverE);
        chain->SetBranchAddress("e_InvEminusInvP", &e_InvEminusInvP);
        chain->SetBranchAddress("e_chIso03", &e_chIso03);
        chain->SetBranchAddress("e_nhIso03", &e_nhIso03);
        chain->SetBranchAddress("e_phIso03", &e_phIso03);
        chain->SetBranchAddress("e_ChIso03FromPU", &e_ChIso03FromPU);
        chain->SetBranchAddress("e_mHits", &e_mHits);
        chain->SetBranchAddress("e_passConvVeto", &e_passConvVeto);
        chain->SetBranchAddress("e_relPFiso_dBeta", &e_relPFiso_dBeta);
        chain->SetBranchAddress("e_relPFiso_Rho", &e_relPFiso_Rho);
        chain->SetBranchAddress("e_passMediumID", &e_passMediumID);
        chain->SetBranchAddress("mu_p_T", &mu_p_T);
        chain->SetBranchAddress("mu_eta", &mu_eta);
        chain->SetBranchAddress("mu_phi", &mu_phi);
        chain->SetBranchAddress("mu_charge", &mu_charge);
        chain->SetBranchAddress("mu_relPFiso_dBeta", &mu_relPFiso_dBeta);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);


        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // W+Jets selection
            if (e_passMediumID->at(0) == 1 && mu_relPFiso_dBeta->at(0) < 0.15) continue;
            if (e_passMediumID->at(0) == 0 && mu_relPFiso_dBeta->at(0) >= 0.15) continue;
            if (e_p_T->at(0) < 17 || mu_p_T->at(0) < 17) continue;
            if (e_p_T->at(0) < 28 && mu_p_T->at(0) < 28) continue;

            if (e_p_T->at(0) != e_p_T->at(0)) cout << e_p_T->at(0) << " " << e_eta->at(0) << " " << e_phi->at(0) << " " << e_charge->at(0) << " " << e_relPFiso_Rho->at(0) << endl;
            if (mu_p_T->at(0) != mu_p_T->at(0)) cout << mu_p_T->at(0) << " " << mu_eta->at(0) << " " << mu_phi->at(0) << " " << mu_charge->at(0) << " " << mu_relPFiso_dBeta->at(0) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T->at(0);
                cout << "\teta = " << e_eta->at(0);
                cout << "\tphi = " << e_phi->at(0) << endl;
                cout << "\tpassMediumID = " << e_passMediumID->at(0) << endl;
                cout << "Muon p_T = " << mu_p_T->at(0);
                cout << "\teta = " << mu_eta->at(0);
                cout << "\tphi = " << mu_phi->at(0) << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta->at(0) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu, ele_SF;
            ele.SetPtEtaPhiM(e_p_T->at(0), e_eta->at(0), e_phi->at(0), M_Elec);
            mu.SetPtEtaPhiM(mu_p_T->at(0), mu_eta->at(0), mu_phi->at(0), M_Mu);
            ele_SF.SetPtEtaPhiM(e_p_T->at(0), e_etaSC->at(0), e_phi->at(0), M_Elec);
            Double_t mass = (ele+mu).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- FR (and SF) WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR = -1;
            if (e_passMediumID->at(0) == 0 && mu_relPFiso_dBeta->at(0) < 0.15) // Electron fails, muon passes
            {
                FR = analyzer->FakeRate_ele(e_p_T->at(0), e_etaSC->at(0));
                FRweight = FR / (1 - FR);
//                if (Mgr.isMC == kTRUE)
//                {
//                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 2);
//                }
                h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }
            else if (e_passMediumID->at(0) == 1 && mu_relPFiso_dBeta->at(0) >= 0.15) // Muon fails, electron passes
            {
                FR = analyzer->FakeRate(mu_p_T->at(0), mu_eta->at(0));
                FRweight = FR / (1 - FR);
//                if (Mgr.isMC == kTRUE)
//                {
//                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 1);
//                }
            }
            if (DEBUG == kTRUE) cout << "FR = " << FR << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Histogram filling -- //
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (e_charge->at(0) != mu_charge->at(0))
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            else
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_nVTX->Write();
        h_mass_SS->Write();
        h_mass_temp->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDEMEnriched_Normal; // next -- GJets_20to100
        if (pr == _GJets_2000to5000) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_EMu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_EMu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_WJET_HistMaker()
