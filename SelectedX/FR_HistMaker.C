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
#include "./header/LocalFileMgr.h"
#include "./header/PrescaleProvider.h"

void E_FR_HistMaker (Bool_t DEBUG);
void E_FR_HistMaker_alt (Bool_t DEBUG);
void E_FR_HistMaker_alt2 (Bool_t DEBUG);
void Mu_FR_HistMaker (Bool_t DEBUG);

void E_QCD_HistMaker (Bool_t DEBUG);
void E_WJET_HistMaker (Bool_t DEBUG);
void Mu_QCD_HistMaker (Bool_t DEBUG);
void Mu_WJET_HistMaker (Bool_t DEBUG);
void EMu_QCD_HistMaker (Bool_t DEBUG);
void EMu_WJET_HistMaker (Bool_t DEBUG);
void E_PR_HistMaker (Bool_t DEBUG); // PROMPT RATE
void Mu_PR_HistMaker (Bool_t DEBUG);

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

const Double_t etabins[67] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2.0,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.566,
                              -1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
                              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,
                              1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
const Double_t etabins2[6] = {0, 0.7, 1.4442, 1.566, 2.0, 2.4};

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
            cout << "\n*****  EMu_QCD_HistMaker()  *****" << endl;
            EMu_QCD_HistMaker(DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  EMu_WJET_HistMaker()  *****" << endl;
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
        if (whichX.Contains("PR"))
        {
            cout << "\n*****  Mu_PR_HistMaker()  *****" << endl;
            Mu_PR_HistMaker(DEBUG);
        }
        else if (whichX.Contains("QCD"))
        {
            cout << "\n*****  Mu_QCD_HistMaker()  *****" << endl;
            Mu_QCD_HistMaker(DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  Mu_WJET_HistMaker()  *****" << endl;
            Mu_WJET_HistMaker(DEBUG);
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
        if (whichX.Contains("PR"))
        {
            cout << "\n*****  E_PR_HistMaker  *****" << endl;
            E_PR_HistMaker(DEBUG);
        }
        else if (whichX.Contains("QCD"))
        {
            cout << "\n*****  E_QCD_HistMaker  *****" << endl;
            E_QCD_HistMaker(DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  E_WJET_HistMaker  *****" << endl;
            E_WJET_HistMaker(DEBUG);
        }
        else if (whichX.Contains("ALT") && whichX.Contains("2"))
        {
            cout << "\n*****    E_FR_HistMaker_alt2    *****" << endl;
            E_FR_HistMaker_alt2(DEBUG);
        }
        else if (whichX.Contains("ALT"))
        {
            cout << "\n*****    E_FR_HistMaker_alt    *****" << endl;
            E_FR_HistMaker_alt(DEBUG);
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
    PrescaleProvider pp("etc/prescale/triggerData2016");

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
        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", nPtBin_ele, analyzer->ptbin_ele); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_endcap2_nume = new TH1D("h_pT_endcap2_nume", "h_pT_endcap2_nume", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap2_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", nPtBin_ele, analyzer->ptbin_ele); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap_deno->Sumw2();
        TH1D* h_pT_endcap2_deno = new TH1D("h_pT_endcap2_deno", "h_pT_endcap2_deno", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap2_deno->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBin_ele, analyzer->ptbin_ele); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_pT_endcap2_ctrl = new TH1D("h_pT_endcap2_ctrl", "h_pT_endcap2_ctrl", nPtBin_ele, analyzer->ptbin_ele); h_pT_endcap2_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 66, etabins); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 66, etabins); h_eta_deno->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 66, etabins); h_eta_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_barrel_nume = new TH1D("h_PFiso_Rho_barrel_nume", "h_PFiso_Rho_barrel_nume", 20, 0, 0.2); h_PFiso_Rho_barrel_nume->Sumw2();
        TH1D* h_PFiso_Rho_endcap_nume = new TH1D("h_PFiso_Rho_endcap_nume", "h_PFiso_Rho_endcap_nume", 20, 0, 0.2); h_PFiso_Rho_endcap_nume->Sumw2();
        TH1D* h_PFiso_Rho_endcap2_nume = new TH1D("h_PFiso_Rho_endcap2_nume", "h_PFiso_Rho_endcap2_nume", 20, 0, 0.2); h_PFiso_Rho_endcap2_nume->Sumw2();
        TH1D* h_PFiso_Rho_barrel_deno = new TH1D("h_PFiso_Rho_barrel_deno", "h_PFiso_Rho_barrel_deno", 50, 0, 5); h_PFiso_Rho_barrel_deno->Sumw2();
        TH1D* h_PFiso_Rho_endcap_deno = new TH1D("h_PFiso_Rho_endcap_deno", "h_PFiso_Rho_endcap_deno", 50, 0, 5); h_PFiso_Rho_endcap_deno->Sumw2();
        TH1D* h_PFiso_Rho_endcap2_deno = new TH1D("h_PFiso_Rho_endcap2_deno", "h_PFiso_Rho_endcap2_deno", 50, 0, 5); h_PFiso_Rho_endcap2_deno->Sumw2();
        TH1D* h_PFiso_Rho_barrel_ctrl = new TH1D("h_PFiso_Rho_barrel_ctrl", "h_PFiso_Rho_barrel_ctrl", 50, 0, 5); h_PFiso_Rho_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_endcap_ctrl = new TH1D("h_PFiso_Rho_endcap_ctrl", "h_PFiso_Rho_endcap_ctrl", 50, 0, 5); h_PFiso_Rho_endcap_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_endcap2_ctrl = new TH1D("h_PFiso_Rho_endcap2_ctrl", "h_PFiso_Rho_endcap2_ctrl", 50, 0, 5); h_PFiso_Rho_endcap2_ctrl->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_nume = new TH1D("h_SigmaIEtaIEta_barrel_nume", "h_SigmaIEtaIEta_barrel_nume", 10, 0, 0.01); h_SigmaIEtaIEta_barrel_nume->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_nume = new TH1D("h_SigmaIEtaIEta_endcap_nume", "h_SigmaIEtaIEta_endcap_nume", 30, 0, 0.03); h_SigmaIEtaIEta_endcap_nume->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap2_nume = new TH1D("h_SigmaIEtaIEta_endcap2_nume", "h_SigmaIEtaIEta_endcap2_nume", 30, 0, 0.03); h_SigmaIEtaIEta_endcap2_nume->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_deno = new TH1D("h_SigmaIEtaIEta_barrel_deno", "h_SigmaIEtaIEta_barrel_deno", 50, 0, 0.05); h_SigmaIEtaIEta_barrel_deno->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_deno = new TH1D("h_SigmaIEtaIEta_endcap_deno", "h_SigmaIEtaIEta_endcap_deno", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap_deno->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap2_deno = new TH1D("h_SigmaIEtaIEta_endcap2_deno", "h_SigmaIEtaIEta_endcap2_deno", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap2_deno->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_ctrl = new TH1D("h_SigmaIEtaIEta_barrel_ctrl", "h_SigmaIEtaIEta_barrel_ctrl", 50, 0, 0.05); h_SigmaIEtaIEta_barrel_ctrl->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_ctrl = new TH1D("h_SigmaIEtaIEta_endcap_ctrl", "h_SigmaIEtaIEta_endcap_ctrl", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap_ctrl->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap2_ctrl = new TH1D("h_SigmaIEtaIEta_endcap2_ctrl", "h_SigmaIEtaIEta_endcap2_ctrl", 50, 0, 0.1);  h_SigmaIEtaIEta_endcap2_ctrl->Sumw2();
        TH1D* h_dEtaInSeed_barrel_nume = new TH1D("h_dEtaInSeed_barrel_nume", "h_dEtaInSeed_barrel_nume", 20, -0.1, 0.1); h_dEtaInSeed_barrel_nume->Sumw2();
        TH1D* h_dEtaInSeed_endcap_nume = new TH1D("h_dEtaInSeed_endcap_nume", "h_dEtaInSeed_endcap_nume", 20, -0.1, 0.1); h_dEtaInSeed_endcap_nume->Sumw2();
        TH1D* h_dEtaInSeed_endcap2_nume = new TH1D("h_dEtaInSeed_endcap2_nume", "h_dEtaInSeed_endcap2_nume", 20, -0.1, 0.1); h_dEtaInSeed_endcap2_nume->Sumw2();
        TH1D* h_dEtaInSeed_barrel_deno = new TH1D("h_dEtaInSeed_barrel_deno", "h_dEtaInSeed_barrel_deno", 100, -1, 1); h_dEtaInSeed_barrel_deno->Sumw2();
        TH1D* h_dEtaInSeed_endcap_deno = new TH1D("h_dEtaInSeed_endcap_deno", "h_dEtaInSeed_endcap_deno", 100, -1, 1); h_dEtaInSeed_endcap_deno->Sumw2();
        TH1D* h_dEtaInSeed_endcap2_deno = new TH1D("h_dEtaInSeed_endcap2_deno", "h_dEtaInSeed_endcap2_deno", 100, -1, 1); h_dEtaInSeed_endcap2_deno->Sumw2();
        TH1D* h_dEtaInSeed_barrel_ctrl = new TH1D("h_dEtaInSeed_barrel_ctrl", "h_dEtaInSeed_barrel_ctrl", 100, -1, 1); h_dEtaInSeed_barrel_ctrl->Sumw2();
        TH1D* h_dEtaInSeed_endcap_ctrl = new TH1D("h_dEtaInSeed_endcap_ctrl", "h_dEtaInSeed_endcap_ctrl", 100, -1, 1); h_dEtaInSeed_endcap_ctrl->Sumw2();
        TH1D* h_dEtaInSeed_endcap2_ctrl = new TH1D("h_dEtaInSeed_endcap2_ctrl", "h_dEtaInSeed_endcap2_ctrl", 100, -1, 1); h_dEtaInSeed_endcap2_ctrl->Sumw2();
        TH1D* h_dPhiIn_barrel_nume = new TH1D("h_dPhiIn_barrel_nume", "h_dPhiIn_barrel_nume", 20, -0.1, 0.1); h_dPhiIn_barrel_nume->Sumw2();
        TH1D* h_dPhiIn_endcap_nume = new TH1D("h_dPhiIn_endcap_nume", "h_dPhiIn_endcap_nume", 20, -0.1, 0.1); h_dPhiIn_endcap_nume->Sumw2();
        TH1D* h_dPhiIn_endcap2_nume = new TH1D("h_dPhiIn_endcap2_nume", "h_dPhiIn_endcap2_nume", 20, -0.1, 0.1); h_dPhiIn_endcap2_nume->Sumw2();
        TH1D* h_dPhiIn_barrel_deno = new TH1D("h_dPhiIn_barrel_deno", "h_dPhiIn_barrel_deno", 20, -0.1, 0.1); h_dPhiIn_barrel_deno->Sumw2();
        TH1D* h_dPhiIn_endcap_deno = new TH1D("h_dPhiIn_endcap_deno", "h_dPhiIn_endcap_deno", 100, -1, 1);    h_dPhiIn_endcap_deno->Sumw2();
        TH1D* h_dPhiIn_endcap2_deno = new TH1D("h_dPhiIn_endcap2_deno", "h_dPhiIn_endcap2_deno", 100, -1, 1);    h_dPhiIn_endcap2_deno->Sumw2();
        TH1D* h_dPhiIn_barrel_ctrl = new TH1D("h_dPhiIn_barrel_ctrl", "h_dPhiIn_barrel_ctrl", 20, -0.1, 0.1); h_dPhiIn_barrel_ctrl->Sumw2();
        TH1D* h_dPhiIn_endcap_ctrl = new TH1D("h_dPhiIn_endcap_ctrl", "h_dPhiIn_endcap_ctrl", 100, -1, 1);    h_dPhiIn_endcap_ctrl->Sumw2();
        TH1D* h_dPhiIn_endcap2_ctrl = new TH1D("h_dPhiIn_endcap2_ctrl", "h_dPhiIn_endcap2_ctrl", 100, -1, 1);    h_dPhiIn_endcap2_ctrl->Sumw2();
        TH1D* h_HoverE_barrel_nume = new TH1D("h_HoverE_barrel_nume", "h_HoverE_barrel_nume", 20, 0, 0.1); h_HoverE_barrel_nume->Sumw2();
        TH1D* h_HoverE_endcap_nume = new TH1D("h_HoverE_endcap_nume", "h_HoverE_endcap_nume", 30, 0, 0.15); h_HoverE_endcap_nume->Sumw2();
        TH1D* h_HoverE_endcap2_nume = new TH1D("h_HoverE_endcap2_nume", "h_HoverE_endcap2_nume", 30, 0, 0.15); h_HoverE_endcap2_nume->Sumw2();
        TH1D* h_HoverE_barrel_deno = new TH1D("h_HoverE_barrel_deno", "h_HoverE_barrel_deno", 100, 0, 1);  h_HoverE_barrel_deno->Sumw2();
        TH1D* h_HoverE_endcap_deno = new TH1D("h_HoverE_endcap_deno", "h_HoverE_endcap_deno", 100, 0, 1); h_HoverE_endcap_deno->Sumw2();
        TH1D* h_HoverE_endcap2_deno = new TH1D("h_HoverE_endcap2_deno", "h_HoverE_endcap2_deno", 100, 0, 1); h_HoverE_endcap2_deno->Sumw2();
        TH1D* h_HoverE_barrel_ctrl = new TH1D("h_HoverE_barrel_ctrl", "h_HoverE_barrel_ctrl", 100, 0, 1);  h_HoverE_barrel_ctrl->Sumw2();
        TH1D* h_HoverE_endcap_ctrl = new TH1D("h_HoverE_endcap_ctrl", "h_HoverE_endcap_ctrl", 100, 0, 1); h_HoverE_endcap_ctrl->Sumw2();
        TH1D* h_HoverE_endcap2_ctrl = new TH1D("h_HoverE_endcap2_ctrl", "h_HoverE_endcap2_ctrl", 100, 0, 1); h_HoverE_endcap2_ctrl->Sumw2();
        TH1D* h_InvEminusInvP_barrel_nume = new TH1D("h_InvEminusInvP_barrel_nume", "h_InvEminusInvP_barrel_nume", 100, -0.5, 0.5); h_InvEminusInvP_barrel_nume->Sumw2();
        TH1D* h_InvEminusInvP_endcap_nume = new TH1D("h_InvEminusInvP_endcap_nume", "h_InvEminusInvP_endcap_nume", 100, -0.5, 0.5); h_InvEminusInvP_endcap_nume->Sumw2();
        TH1D* h_InvEminusInvP_endcap2_nume = new TH1D("h_InvEminusInvP_endcap2_nume", "h_InvEminusInvP_endcap2_nume", 100, -0.5, 0.5); h_InvEminusInvP_endcap2_nume->Sumw2();
        TH1D* h_InvEminusInvP_barrel_deno = new TH1D("h_InvEminusInvP_barrel_deno", "h_InvEminusInvP_barrel_deno", 120, -6, 6); h_InvEminusInvP_barrel_deno->Sumw2();
        TH1D* h_InvEminusInvP_endcap_deno = new TH1D("h_InvEminusInvP_endcap_deno", "h_InvEminusInvP_endcap_deno", 120, -6, 6); h_InvEminusInvP_endcap_deno->Sumw2();
        TH1D* h_InvEminusInvP_endcap2_deno = new TH1D("h_InvEminusInvP_endcap2_deno", "h_InvEminusInvP_endcap2_deno", 120, -6, 6); h_InvEminusInvP_endcap2_deno->Sumw2();
        TH1D* h_InvEminusInvP_barrel_ctrl = new TH1D("h_InvEminusInvP_barrel_ctrl", "h_InvEminusInvP_barrel_ctrl", 120, -6, 6); h_InvEminusInvP_barrel_ctrl->Sumw2();
        TH1D* h_InvEminusInvP_endcap_ctrl = new TH1D("h_InvEminusInvP_endcap_ctrl", "h_InvEminusInvP_endcap_ctrl", 120, -6, 6); h_InvEminusInvP_endcap_ctrl->Sumw2();
        TH1D* h_InvEminusInvP_endcap2_ctrl = new TH1D("h_InvEminusInvP_endcap2_ctrl", "h_InvEminusInvP_endcap2_ctrl", 120, -6, 6); h_InvEminusInvP_endcap2_ctrl->Sumw2();
        TH1D* h_TrkIso_barrel_nume = new TH1D("h_TrkIso_barrel_nume", "h_TrkIso_barrel_nume", 50, 0, 5); h_TrkIso_barrel_nume->Sumw2();
        TH1D* h_TrkIso_endcap_nume = new TH1D("h_TrkIso_endcap_nume", "h_TrkIso_endcap_nume", 50, 0, 5); h_TrkIso_endcap_nume->Sumw2();
        TH1D* h_TrkIso_endcap2_nume = new TH1D("h_TrkIso_endcap2_nume", "h_TrkIso_endcap2_nume", 50, 0, 5); h_TrkIso_endcap2_nume->Sumw2();
        TH1D* h_TrkIso_barrel_deno = new TH1D("h_TrkIso_barrel_deno", "h_TrkIso_barrel_deno", 50, 0, 5); h_TrkIso_barrel_deno->Sumw2();
        TH1D* h_TrkIso_endcap_deno = new TH1D("h_TrkIso_endcap_deno", "h_TrkIso_endcap_deno", 50, 0, 5); h_TrkIso_endcap_deno->Sumw2();
        TH1D* h_TrkIso_endcap2_deno = new TH1D("h_TrkIso_endcap2_deno", "h_TrkIso_endcap2_deno", 50, 0, 5); h_TrkIso_endcap2_deno->Sumw2();
        TH1D* h_TrkIso_barrel_ctrl = new TH1D("h_TrkIso_barrel_ctrl", "h_TrkIso_barrel_ctrl", 50, 0, 5); h_TrkIso_barrel_ctrl->Sumw2();
        TH1D* h_TrkIso_endcap_ctrl = new TH1D("h_TrkIso_endcap_ctrl", "h_TrkIso_endcap_ctrl", 50, 0, 5); h_TrkIso_endcap_ctrl->Sumw2();
        TH1D* h_TrkIso_endcap2_ctrl = new TH1D("h_TrkIso_endcap2_ctrl", "h_TrkIso_endcap2_ctrl", 50, 0, 5); h_TrkIso_endcap2_ctrl->Sumw2();
        TH1D* h_ECALiso_barrel_nume = new TH1D("h_ECALiso_barrel_nume", "h_ECALiso_barrel_nume", 50, 0, 5); h_ECALiso_barrel_nume->Sumw2();
        TH1D* h_ECALiso_endcap_nume = new TH1D("h_ECALiso_endcap_nume", "h_ECALiso_endcap_nume", 50, 0, 5); h_ECALiso_endcap_nume->Sumw2();
        TH1D* h_ECALiso_endcap2_nume = new TH1D("h_ECALiso_endcap2_nume", "h_ECALiso_endcap2_nume", 50, 0, 5); h_ECALiso_endcap2_nume->Sumw2();
        TH1D* h_ECALiso_barrel_deno = new TH1D("h_ECALiso_barrel_deno", "h_ECALiso_barrel_deno", 50, 0, 5); h_ECALiso_barrel_deno->Sumw2();
        TH1D* h_ECALiso_endcap_deno = new TH1D("h_ECALiso_endcap_deno", "h_ECALiso_endcap_deno", 50, 0, 5); h_ECALiso_endcap_deno->Sumw2();
        TH1D* h_ECALiso_endcap2_deno = new TH1D("h_ECALiso_endcap2_deno", "h_ECALiso_endcap2_deno", 50, 0, 5); h_ECALiso_endcap2_deno->Sumw2();
        TH1D* h_ECALiso_barrel_ctrl = new TH1D("h_ECALiso_barrel_ctrl", "h_ECALiso_barrel_ctrl", 50, 0, 5); h_ECALiso_barrel_ctrl->Sumw2();
        TH1D* h_ECALiso_endcap_ctrl = new TH1D("h_ECALiso_endcap_ctrl", "h_ECALiso_endcap_ctrl", 50, 0, 5); h_ECALiso_endcap_ctrl->Sumw2();
        TH1D* h_ECALiso_endcap2_ctrl = new TH1D("h_ECALiso_endcap2_ctrl", "h_ECALiso_endcap2_ctrl", 50, 0, 5); h_ECALiso_endcap2_ctrl->Sumw2();
        TH1D* h_HCALiso_barrel_nume = new TH1D("h_HCALiso_barrel_nume", "h_HCALiso_barrel_nume", 50, 0, 5); h_HCALiso_barrel_nume->Sumw2();
        TH1D* h_HCALiso_endcap_nume = new TH1D("h_HCALiso_endcap_nume", "h_HCALiso_endcap_nume", 50, 0, 5); h_HCALiso_endcap_nume->Sumw2();
        TH1D* h_HCALiso_endcap2_nume = new TH1D("h_HCALiso_endcap2_nume", "h_HCALiso_endcap2_nume", 50, 0, 5); h_HCALiso_endcap2_nume->Sumw2();
        TH1D* h_HCALiso_barrel_deno = new TH1D("h_HCALiso_barrel_deno", "h_HCALiso_barrel_deno", 50, 0, 5); h_HCALiso_barrel_deno->Sumw2();
        TH1D* h_HCALiso_endcap_deno = new TH1D("h_HCALiso_endcap_deno", "h_HCALiso_endcap_deno", 50, 0, 5); h_HCALiso_endcap_deno->Sumw2();
        TH1D* h_HCALiso_endcap2_deno = new TH1D("h_HCALiso_endcap2_deno", "h_HCALiso_endcap2_deno", 50, 0, 5); h_HCALiso_endcap2_deno->Sumw2();
        TH1D* h_HCALiso_barrel_ctrl = new TH1D("h_HCALiso_barrel_ctrl", "h_HCALiso_barrel_ctrl", 50, 0, 5); h_HCALiso_barrel_ctrl->Sumw2();
        TH1D* h_HCALiso_endcap_ctrl = new TH1D("h_HCALiso_endcap_ctrl", "h_HCALiso_endcap_ctrl", 50, 0, 5); h_HCALiso_endcap_ctrl->Sumw2();
        TH1D* h_HCALiso_endcap2_ctrl = new TH1D("h_HCALiso_endcap2_ctrl", "h_HCALiso_endcap2_ctrl", 50, 0, 5); h_HCALiso_endcap2_ctrl->Sumw2();
        TH1D* h_MET = new TH1D("h_MET", "h_MET", 100, 0, 100); h_MET->Sumw2();
        TH1D* h_MT_barrel_nume = new TH1D("h_MT_barrel_nume", "h_MT_barrel_nume", 500, 0, 1000); h_MT_barrel_nume->Sumw2();
        TH1D* h_MT_endcap_nume = new TH1D("h_MT_endcap_nume", "h_MT_endcap_nume", 500, 0, 1000); h_MT_endcap_nume->Sumw2();
        TH1D* h_MT_endcap2_nume = new TH1D("h_MT_endcap2_nume", "h_MT_endcap2_nume", 500, 0, 1000); h_MT_endcap2_nume->Sumw2();
        TH1D* h_MT_barrel_deno = new TH1D("h_MT_barrel_deno", "h_MT_barrel_deno", 500, 0, 1000); h_MT_barrel_deno->Sumw2();
        TH1D* h_MT_endcap_deno = new TH1D("h_MT_endcap_deno", "h_MT_endcap_deno", 500, 0, 1000); h_MT_endcap_deno->Sumw2();
        TH1D* h_MT_endcap2_deno = new TH1D("h_MT_endcap2_deno", "h_MT_endcap2_deno", 500, 0, 1000); h_MT_endcap2_deno->Sumw2();
        TH1D* h_MT_barrel_ctrl = new TH1D("h_MT_barrel_ctrl", "h_MT_barrel_ctrl", 500, 0, 1000); h_MT_barrel_ctrl->Sumw2();
        TH1D* h_MT_endcap_ctrl = new TH1D("h_MT_endcap_ctrl", "h_MT_endcap_ctrl", 500, 0, 1000); h_MT_endcap_ctrl->Sumw2();
        TH1D* h_MT_endcap2_ctrl = new TH1D("h_MT_endcap2_ctrl", "h_MT_endcap2_ctrl", 500, 0, 1000); h_MT_endcap2_ctrl->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();

        TH1D* h_HoverE_barrel_template_int = new TH1D("h_HoverE_barrel_template_int", "h_HoverE_barrel_template_int", binnum_ele_template, bins_ele_HE_barrel); h_HoverE_barrel_template_int->Sumw2();
        TH1D* h_HoverE_endcap_template_int = new TH1D("h_HoverE_endcap_template_int", "h_HoverE_endcap_template_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap_template_int->Sumw2();
        TH1D* h_HoverE_endcap2_template_int = new TH1D("h_HoverE_endcap2_template_int", "h_HoverE_endcap2_template_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap2_template_int->Sumw2();
        TH1D* h_HoverE_barrel_jetTemplate_int = new TH1D("h_HoverE_barrel_jetTemplate_int", "h_HoverE_barrel_jetTemplate_int", binnum_ele_template, bins_ele_HE_barrel); h_HoverE_barrel_jetTemplate_int->Sumw2();
        TH1D* h_HoverE_endcap_jetTemplate_int = new TH1D("h_HoverE_endcap_jetTemplate_int", "h_HoverE_endcap_jetTemplate_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap_jetTemplate_int->Sumw2();
        TH1D* h_HoverE_endcap2_jetTemplate_int = new TH1D("h_HoverE_endcap2_jetTemplate_int", "h_HoverE_endcap2_jetTemplate_int", binnum_ele_template, bins_ele_HE_endcap); h_HoverE_endcap2_jetTemplate_int->Sumw2();

        TH1D *h_PFiso_Rho_barrel_template[nPtBin_ele],       *h_PFiso_Rho_endcap_template[nPtBin_ele],       *h_PFiso_Rho_endcap2_template[nPtBin_ele],
             *h_PFiso_Rho_barrel_jetTemplate[nPtBin_ele],    *h_PFiso_Rho_endcap_jetTemplate[nPtBin_ele],    *h_PFiso_Rho_endcap2_jetTemplate[nPtBin_ele],
             *h_PFiso_Rho_barrel_badJetTemplate[nPtBin_ele], *h_PFiso_Rho_endcap_badJetTemplate[nPtBin_ele], *h_PFiso_Rho_endcap2_badJetTemplate[nPtBin_ele],
             *h_HoverE_barrel_template[nPtBin_ele],          *h_HoverE_endcap_template[nPtBin_ele],          *h_HoverE_endcap2_template[nPtBin_ele],
             *h_HoverE_barrel_jetTemplate[nPtBin_ele],       *h_HoverE_endcap_jetTemplate[nPtBin_ele],       *h_HoverE_endcap2_jetTemplate[nPtBin_ele];
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
            h_PFiso_Rho_endcap2_template[ih] = new TH1D("h_PFiso_Rho_endcap2_template_"+TString::Itoa(ih, 10),
                                                       "h_PFiso_Rho_endcap2_template_"+TString::Itoa(ih, 10),
                                                       binnum_ele_template, bins_ele_template_endcap);
            h_PFiso_Rho_endcap2_jetTemplate[ih] = new TH1D("h_PFiso_Rho_endcap2_jetTemplate_"+TString::Itoa(ih, 10),
                                                          "h_PFiso_Rho_endcap2_jetTemplate_"+TString::Itoa(ih, 10),
                                                          binnum_ele_template, bins_ele_template_endcap);
            h_PFiso_Rho_endcap2_badJetTemplate[ih] = new TH1D("h_PFiso_Rho_endcap2_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             "h_PFiso_Rho_endcap2_badJetTemplate_"+TString::Itoa(ih, 10),
                                                             binnum_ele_template, bins_ele_template_endcap);
            h_HoverE_endcap2_template[ih] = new TH1D("h_HoverE_endcap2_template_"+TString::Itoa(ih, 10),
                                                    "h_HoverE_endcap2_template_"+TString::Itoa(ih, 10),
                                                    binnum_ele_template, bins_ele_HE_endcap);
            h_HoverE_endcap2_jetTemplate[ih] = new TH1D("h_HoverE_endcap2_jetTemplate_"+TString::Itoa(ih, 10),
                                                       "h_HoverE_endcap2_jetTemplate_"+TString::Itoa(ih, 10),
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
            h_PFiso_Rho_endcap2_template[ih]->Sumw2();
            h_PFiso_Rho_endcap2_jetTemplate[ih]->Sumw2();
            h_PFiso_Rho_endcap2_badJetTemplate[ih]->Sumw2();
            h_HoverE_endcap2_template[ih]->Sumw2();
            h_HoverE_endcap2_jetTemplate[ih]->Sumw2();
        }

        TH1D* h_PFiso_Rho_barrel_separate = new TH1D("h_PFiso_Rho_barrel_separate", "h_PFiso_Rho_barrel_separate", 80, 0, 0.8); h_PFiso_Rho_barrel_separate->Sumw2();
        TH1D* h_PFiso_Rho_endcap_separate = new TH1D("h_PFiso_Rho_endcap_separate", "h_PFiso_Rho_endcap_separate", 80, 0, 0.8); h_PFiso_Rho_endcap_separate->Sumw2();
        TH1D* h_PFiso_Rho_endcap2_separate = new TH1D("h_PFiso_Rho_endcap2_separate", "h_PFiso_Rho_endcap2_separate", 80, 0, 0.8); h_PFiso_Rho_endcap2_separate->Sumw2();
        TH1D* h_SigmaIEtaIEta_barrel_separate = new TH1D("h_SigmaIEtaIEta_barrel_separate", "h_SigmaIEtaIEta_barrel_separate", 13, 0, 0.013); h_SigmaIEtaIEta_barrel_separate->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap_separate = new TH1D("h_SigmaIEtaIEta_endcap_separate", "h_SigmaIEtaIEta_endcap_separate", 35, 0, 0.035);  h_SigmaIEtaIEta_endcap_separate->Sumw2();
        TH1D* h_SigmaIEtaIEta_endcap2_separate = new TH1D("h_SigmaIEtaIEta_endcap2_separate", "h_SigmaIEtaIEta_endcap2_separate", 35, 0, 0.035);  h_SigmaIEtaIEta_endcap2_separate->Sumw2();
        TH1D* h_dEtaInSeed_barrel_separate = new TH1D("h_dEtaInSeed_barrel_separate", "h_dEtaInSeed_barrel_separate", 10, 0, 0.01); h_dEtaInSeed_barrel_separate->Sumw2();
        TH1D* h_dEtaInSeed_endcap_separate = new TH1D("h_dEtaInSeed_endcap_separate", "h_dEtaInSeed_endcap_separate", 30, 0, 0.03); h_dEtaInSeed_endcap_separate->Sumw2();
        TH1D* h_dEtaInSeed_endcap2_separate = new TH1D("h_dEtaInSeed_endcap2_separate", "h_dEtaInSeed_endcap2_separate", 30, 0, 0.03); h_dEtaInSeed_endcap2_separate->Sumw2();
        TH1D* h_dPhiIn_barrel_separate = new TH1D("h_dPhiIn_barrel_separate", "h_dPhiIn_barrel_separate", 10, 0, 0.1); h_dPhiIn_barrel_separate->Sumw2();
        TH1D* h_dPhiIn_endcap_separate = new TH1D("h_dPhiIn_endcap_separate", "h_dPhiIn_endcap_separate", 30, 0, 0.3);    h_dPhiIn_endcap_separate->Sumw2();
        TH1D* h_dPhiIn_endcap2_separate = new TH1D("h_dPhiIn_endcap2_separate", "h_dPhiIn_endcap2_separate", 30, 0, 0.3);    h_dPhiIn_endcap2_separate->Sumw2();
        TH1D* h_HoverE_barrel_separate = new TH1D("h_HoverE_barrel_separate", "h_HoverE_barrel_separate", 15, 0, 0.15);  h_HoverE_barrel_separate->Sumw2();
        TH1D* h_HoverE_endcap_separate = new TH1D("h_HoverE_endcap_separate", "h_HoverE_endcap_separate", 15, 0, 0.15); h_HoverE_endcap_separate->Sumw2();
        TH1D* h_HoverE_endcap2_separate = new TH1D("h_HoverE_endcap2_separate", "h_HoverE_endcap2_separate", 15, 0, 0.15); h_HoverE_endcap2_separate->Sumw2();
        TH1D* h_InvEminusInvP_barrel_separate = new TH1D("h_InvEminusInvP_barrel_separate", "h_InvEminusInvP_barrel_separate", 10, 0, 1); h_InvEminusInvP_barrel_separate->Sumw2();
        TH1D* h_InvEminusInvP_endcap_separate = new TH1D("h_InvEminusInvP_endcap_separate", "h_InvEminusInvP_endcap_separate", 10, 0, 1); h_InvEminusInvP_endcap_separate->Sumw2();
        TH1D* h_InvEminusInvP_endcap2_separate = new TH1D("h_InvEminusInvP_endcap2_separate", "h_InvEminusInvP_endcap2_separate", 10, 0, 1); h_InvEminusInvP_endcap2_separate->Sumw2();
        TH1D* h_mHits_barrel_separate = new TH1D("h_mHits_barrel_separate", "h_mHits_barrel_separate", 3, 0-0.5, 3-0.5); h_mHits_barrel_separate->Sumw2();
        TH1D* h_mHits_endcap_separate = new TH1D("h_mHits_endcap_separate", "h_mHits_endcap_separate", 3, 0-0.5, 3-0.5); h_mHits_endcap_separate->Sumw2();
        TH1D* h_mHits_endcap2_separate = new TH1D("h_mHits_endcap2_separate", "h_mHits_endcap2_separate", 3, 0-0.5, 3-0.5); h_mHits_endcap2_separate->Sumw2();
        TH1D* h_passConvVeto_barrel_separate = new TH1D("h_passConvVeto_barrel_separate", "h_passConvVeto_barrel_separate", 2, 0-0.5, 2-0.5); h_passConvVeto_barrel_separate->Sumw2();
        TH1D* h_passConvVeto_endcap_separate = new TH1D("h_passConvVeto_endcap_separate", "h_passConvVeto_endcap_separate", 2, 0-0.5, 2-0.5); h_passConvVeto_endcap_separate->Sumw2();
        TH1D* h_passConvVeto_endcap2_separate = new TH1D("h_passConvVeto_endcap2_separate", "h_passConvVeto_endcap2_separate", 2, 0-0.5, 2-0.5); h_passConvVeto_endcap2_separate->Sumw2();
        TH1D* h_passMediumID_barrel_separate = new TH1D("h_passMediumID_barrel_separate", "h_passMediumID_barrel_separate", 2, 0-0.5, 2-0.5); h_passMediumID_barrel_separate->Sumw2();
        TH1D* h_passMediumID_endcap_separate = new TH1D("h_passMediumID_endcap_separate", "h_passMediumID_endcap_separate", 2, 0-0.5, 2-0.5); h_passMediumID_endcap_separate->Sumw2();
        TH1D* h_passMediumID_endcap2_separate = new TH1D("h_passMediumID_endcap2_separate", "h_passMediumID_endcap2_separate", 2, 0-0.5, 2-0.5); h_passMediumID_endcap2_separate->Sumw2();

        TH1D* h_mass_test = new TH1D("h_mass_test", "h_mass_test", binnum, massbins); h_mass_test->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
        std::vector<double> *relTrkIso = new std::vector<double>;
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
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("relTrkIso", 1);
        chain->SetBranchStatus("relECALiso", 1);
        chain->SetBranchStatus("relHCALiso", 1);
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
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("relTrkIso", &relTrkIso);
        chain->SetBranchAddress("relECALiso", &relECALiso);
        chain->SetBranchAddress("relHCALiso", &relHCALiso);
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
            if (DEBUG == kFALSE) bar.Draw(i);
            if (MET_pT >= 25) continue;
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
                if (fabs(eta->at(i_ele)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.013 || HoverE->at(i_ele) >= 0.13 ||
                                                      fabs(dEtaInSeed->at(i_ele)) >= 0.01 || fabs(dPhiIn->at(i_ele)) >= 0.07)) continue;
                else if (fabs(eta->at(i_ele)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.035 || HoverE->at(i_ele) >= 0.13)) continue;
                if (mHits->at(i_ele) > 1) continue;


                // Selecting leading electron (could also try finding a electron with the best isolation)
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
                    iso_lead = relPFiso_Rho->at(i_ele);
                }
            }
            if (i_lead<0) continue;

            Double_t dTheta = ele_lead.Phi() - MET_phi;
            Double_t MT = sqrt(2 * ele_lead.Pt() * MET_pT * (1 - cos(dTheta)));

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
                    cout << p_T->at(i_ele) << " " << eta->at(i_ele) << " " << phi->at(i_ele) << " " << relPFiso_Rho->at(i_ele) << endl;
                    continue;
                }
                if (p_T->at(i_ele) <= 25) continue;
                if (fabs(eta->at(i_ele)) >= 1.4442 && fabs(eta->at(i_ele)) <= 1.566) continue;
                if (fabs(eta->at(i_ele)) >= 2.4) continue;

                if (fabs(eta->at(i_ele)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.013 || HoverE->at(i_ele) >= 0.13 ||
                                                      fabs(dEtaInSeed->at(i_ele)) >= 0.01 || fabs(dPhiIn->at(i_ele)) >= 0.07)) continue;
                else if (fabs(eta->at(i_ele)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.035 || HoverE->at(i_ele) >= 0.13)) continue;

//                if (fabs(eta->at(i_ele)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.00998 || fabs(dEtaInSeed->at(i_ele)) >= 0.00311 ||
//                    fabs(dPhiIn->at(i_ele)) >= 0.07/*0.103*/ || HoverE->at(i_ele) >= 0.13/*0.253*/ || InvEminusInvP->at(i_ele) >= 0.134))
//                    continue;
//                if (fabs(eta->at(i_ele)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 || fabs(dEtaInSeed->at(i_ele)) >= 0.00609 ||
//                    fabs(dPhiIn->at(i_ele)) >= 0.045 || HoverE->at(i_ele) >= 0.0878 || InvEminusInvP->at(i_ele) >= 0.13 ))
//                    continue;
                if (mHits->at(i_ele) > 1 /*|| !passConvVeto->at(i_ele)*/ || relECALiso->at(i_ele) >= 0.5 || relHCALiso->at(i_ele) >= 0.5 || relTrkIso->at(i_ele) >= 0.2) continue;

                if (DEBUG == kTRUE) cout << "i_ele = " << i_ele << endl;

                Int_t matched22=0, matched30=0, matched36=0, matched50=0, matched75=0, matched90=0, matched120=0, matched175=0;
                Int_t i_22=-1, i_30=-1, i_36=-1, i_50=-1, i_75=-1, i_90=-1, i_120=-1, i_175=-1;
                Double_t prescale_alt = 1.;

                if (Mgr.isMC == kFALSE)
                {                        
                    for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                    {
                        if (trig_pT->at(i_tr) < 22) continue;
                        if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
                        {
                            if (trig_fired->at(i_tr) == 22 && trig_pT->at(i_tr) > 22 && trig_pT->at(i_tr) <= 30)
                            {
                                i_22 = i_tr;
                                matched22 = 1;
                            }
                            if (trig_fired->at(i_tr) == 30 && trig_pT->at(i_tr) > 30 && trig_pT->at(i_tr) <= 36)
                            {
                                i_30 = i_tr;
                                matched30 = 1;
                            }
                            if (trig_fired->at(i_tr) == 36 && trig_pT->at(i_tr) > 36 && trig_pT->at(i_tr) <= 50)
                            {
                                i_36 = i_tr;
                                matched36 = 1;
                            }
                            if (trig_fired->at(i_tr) == 50 && trig_pT->at(i_tr) > 50 && trig_pT->at(i_tr) <= 75)
                            {
                                i_50 = i_tr;
                                matched50 = 1;
                            }
                            if (trig_fired->at(i_tr) == 75 && trig_pT->at(i_tr) > 75 && trig_pT->at(i_tr) <= 90)
                            {
                                i_75 = i_tr;
                                matched75 = 1;
                            }
                            if (trig_fired->at(i_tr) == 90 && trig_pT->at(i_tr) > 90 && trig_pT->at(i_tr) <= 120)
                            {
                                i_90 = i_tr;
                                matched90 = 1;
                            }
                            if (trig_fired->at(i_tr) == 120 && trig_pT->at(i_tr) > 120 && trig_pT->at(i_tr) <= 175)
                            {
                                i_120 = i_tr;
                                matched120 = 1;
                            }
                            if (trig_fired->at(i_tr) == 175 && trig_pT->at(i_tr) > 175)
                            {
                                i_175 = i_tr;
                                matched175 = 1;
                            }
                        } // End of if (trig matched with i_ele)
                    } // End of for(i_tr)

                    if (matched22==0 && matched30==0 && matched36==0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0) continue;

                    if (matched22) // Matched an electron to HLT_Photon22 with 22<HLT_pT<30
                    { // get the prescale
                        prescale_alt = (double)(pp.hltPrescale("HLT_Photon22_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG18", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched30) // Matched an electron to HLT_Photon22 with 30<HLT_pT<36
                    {
                        prescale_alt = (double)(pp.hltPrescale("HLT_Photon30_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched36)
                    {
                        prescale_alt = (double)(pp.hltPrescale("HLT_Photon36_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched50)
                    {
                        prescale_alt = 0;
                        if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                            prescale_alt = pp.hltPrescale("HLT_Photon50_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                        if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                            prescale_alt = pp.hltPrescale("HLT_Photon50_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched75)
                    {
                        prescale_alt = 0.0;
                        if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon75_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock));
                        if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon75_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched90)
                    {
                        prescale_alt = 0;
                        if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon90_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock));
                        if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon90_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched120)
                    {
                        prescale_alt = 0;
                        if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon120_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock));
                        if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                            prescale_alt = (double)(pp.hltPrescale("HLT_Photon120_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock));
                        if (Mgr.CurrentProc == _SinglePhoton_B) prescale_alt *= 1.37;
                    }
                    else if (matched175)
                    {
                         prescale_alt = 0;
                         if (pp.l1Prescale("L1_SingleEG30", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG30", runNum, lumiBlock));
                         if (pp.l1Prescale("L1_SingleEG32", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG32", runNum, lumiBlock));
                         if (pp.l1Prescale("L1_SingleEG34", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG34", runNum, lumiBlock));
                         if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock));
                         if (pp.l1Prescale("L1_SingleEG38", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG38", runNum, lumiBlock));
                         if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                             prescale_alt = (double)(pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock));
                    }
                    else continue;
                }// End of if(MC=False)

                if (passMediumID->at(i_ele)) // Signal/Numerator
                {
                    h_eta_nume->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    if (fabs(eta->at(i_ele)) < 1.4442) // Barrel
                    {
                        h_pT_barrel_nume->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_barrel_nume->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_barrel_nume->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_barrel_nume->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_barrel_nume->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_barrel_nume->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_barrel_nume->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_barrel_nume->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_barrel_nume->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_barrel_nume->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) > 1.566 && fabs(eta->at(i_ele)) < 2.2) // Endcap
                    {
                        h_pT_endcap_nume->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap_nume->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap_nume->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap_nume->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap_nume->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap_nume->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap_nume->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap_nume->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap_nume->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap_nume->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) >= 2.2 && fabs(eta->at(i_ele)) < 2.4) // Far endcap
                    {
                        h_pT_endcap2_nume->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap2_nume->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap2_nume->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap2_nume->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap2_nume->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap2_nume->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap2_nume->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap2_nume->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap2_nume->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap2_nume->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                } // End of if(Signal/Numerator)
                else // Control
                {
                    h_eta_ctrl->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    if (fabs(eta->at(i_ele)) < 1.4442) // Barrel
                    {
                        h_pT_barrel_ctrl->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_barrel_ctrl->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_barrel_ctrl->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_barrel_ctrl->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_barrel_ctrl->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight  * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_barrel_ctrl->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_barrel_ctrl->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_barrel_ctrl->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_barrel_ctrl->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_barrel_ctrl->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) > 1.566 && fabs(eta->at(i_ele)) < 2.2) // Endcap
                    {
                        h_pT_endcap_ctrl->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap_ctrl->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap_ctrl->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap_ctrl->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap_ctrl->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap_ctrl->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap_ctrl->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap_ctrl->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap_ctrl->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap_ctrl->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                    else if (fabs(eta->at(i_ele)) >= 2.2 && fabs(eta->at(i_ele)) < 2.4) // Far endcap
                    {
                        h_pT_endcap2_ctrl->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap2_ctrl->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap2_ctrl->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap2_ctrl->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap2_ctrl->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap2_ctrl->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap2_ctrl->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap2_ctrl->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap2_ctrl->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap2_ctrl->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }
                }// End of if(Control)
                // Denominator
                if (fabs(eta->at(i_ele)) < 1.4442) // Barrel
                {
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.013 && HoverE->at(i_ele) < 0.13 &&
                        fabs(dEtaInSeed->at(i_ele)) < 0.01 && fabs(dPhiIn->at(i_ele)) < 0.07 && mHits->at(i_ele) <= 1)
                    {
                        h_eta_deno->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_pT_barrel_deno->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_barrel_deno->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_barrel_deno->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_barrel_deno->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_barrel_deno->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_barrel_deno->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_barrel_deno->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_barrel_deno->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_barrel_deno->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_barrel_deno->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }

                    /// TEST PART ///
                    // MediumID minus SigmaIEtaIEta
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.013 // relaxed
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_SigmaIEtaIEta_barrel_separate->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dEtaInSeed
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.01 // relaxed
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dEtaInSeed_barrel_separate->Fill(fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dPhiIn
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dPhiIn_barrel_separate->Fill(fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus RelPFIso
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_PFiso_Rho_barrel_separate->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus HoverE
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_HoverE_barrel_separate->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus InvEminusInvP
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_InvEminusInvP_barrel_separate->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus mHits
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && passConvVeto->at(i_ele))
                        h_mHits_barrel_separate->Fill(mHits->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus passConvVeto
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00311
                        && fabs(dPhiIn->at(i_ele)) < 0.07 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.13 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0695
                        && InvEminusInvP->at(i_ele) < 0.134
                        && mHits->at(i_ele) <= 1)
                        h_passConvVeto_barrel_separate->Fill(passConvVeto->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID yes or no
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.013
                        && HoverE->at(i_ele) < 0.13
                        && fabs(dEtaInSeed->at(i_ele)) < 0.01
                        && fabs(dPhiIn->at(i_ele)) < 0.07
                        && mHits->at(i_ele) <= 1)
                        h_passMediumID_barrel_separate->Fill(passMediumID->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    /// End test part ///

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                        HoverE->at(i_ele) < 0.253 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_template[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (HoverE->at(i_ele) < 0.253 && (fabs(InvEminusInvP->at(i_ele)) >= 0.134 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.00998 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00311 || fabs(dPhiIn->at(i_ele)) >= 0.103 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_jetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 && fabs(dPhiIn->at(i_ele)) < 0.103 &&
                             HoverE->at(i_ele) >= 0.253 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_barrel_badJetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }

                    if (relPFiso_Rho->at(i_ele) < 0.0695 && Full5x5_SigmaIEtaIEta->at(i_ele) < 0.00998 && fabs(dEtaInSeed->at(i_ele)) < 0.00311 &&
                        fabs(dPhiIn->at(i_ele)) < 0.103 && fabs(InvEminusInvP->at(i_ele)) < 0.134 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h_HoverE_barrel_template_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_barrel_template[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (relPFiso_Rho->at(i_ele) < 0.0695 && (fabs(InvEminusInvP->at(i_ele)) >= 0.134 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.00998 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00311 || fabs(dPhiIn->at(i_ele)) >= 0.103 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        h_HoverE_barrel_jetTemplate_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_barrel_jetTemplate[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                }// End of if(barrel)
                else if (fabs(eta->at(i_ele)) > 1.566 && fabs(eta->at(i_ele)) < 2.2) // Endcap
                {
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 && HoverE->at(i_ele) < 0.13 && mHits->at(i_ele) <= 1)
                    {
                        h_eta_deno->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_pT_endcap_deno->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap_deno->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap_deno->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap_deno->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap_deno->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap_deno->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap_deno->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap_deno->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap_deno->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap_deno->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }

                    /// TEST PART ///
                    // MediumID minus SigmaIEtaIEta
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 // relaxed
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_SigmaIEtaIEta_endcap_separate->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dEtaInSeed
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dEtaInSeed_endcap_separate->Fill(fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dPhiIn
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dPhiIn_endcap_separate->Fill(fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus HoverE
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && HoverE->at(i_ele) < 0.13 // relaxed (Ele23Ele12 threshold)
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_HoverE_endcap_separate->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus RelPFIso
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_PFiso_Rho_endcap_separate->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus InvEminusInvP
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.0878 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_InvEminusInvP_endcap_separate->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus mHits
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && passConvVeto->at(i_ele))
                        h_mHits_endcap_separate->Fill(mHits->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus passConvVeto
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1)
                        h_passConvVeto_endcap_separate->Fill(passConvVeto->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID yes or no
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035
                        && HoverE->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1)
                        h_passMediumID_endcap_separate->Fill(passMediumID->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    /// End test part ///

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_template[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (HoverE->at(i_ele) < 0.0878 && (fabs(InvEminusInvP->at(i_ele)) >= 0.13 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00609 || fabs(dPhiIn->at(i_ele)) >= 0.045 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_jetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                             HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) >= 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap_badJetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        relPFiso_Rho->at(i_ele) < 0.0821 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h_HoverE_endcap_template_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap_template[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (relPFiso_Rho->at(i_ele) < 0.0821 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 || fabs(dEtaInSeed->at(i_ele)) >= 0.00609 ||
                             fabs(dPhiIn->at(i_ele)) >= 0.045 || fabs(InvEminusInvP->at(i_ele)) >= 0.13 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        h_HoverE_endcap_jetTemplate_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap_jetTemplate[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                }// End of if(endcap)
                else if (fabs(eta->at(i_ele)) >= 2.2 && fabs(eta->at(i_ele)) < 2.4) // Far endcap
                {
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 && HoverE->at(i_ele) < 0.13 && mHits->at(i_ele) <= 1)
                    {
                        h_eta_deno->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_pT_endcap2_deno->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_PFiso_Rho_endcap2_deno->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_SigmaIEtaIEta_endcap2_deno->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dEtaInSeed_endcap2_deno->Fill(dEtaInSeed->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_dPhiIn_endcap2_deno->Fill(dPhiIn->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HoverE_endcap2_deno->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_InvEminusInvP_endcap2_deno->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_TrkIso_endcap2_deno->Fill(relTrkIso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_ECALiso_endcap2_deno->Fill(relECALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        h_HCALiso_endcap2_deno->Fill(relHCALiso->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                    }

                    /// TEST PART ///
                    // MediumID minus SigmaIEtaIEta
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035 // relaxed
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_SigmaIEtaIEta_endcap2_separate->Fill(Full5x5_SigmaIEtaIEta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dEtaInSeed
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dEtaInSeed_endcap2_separate->Fill(fabs(dEtaInSeed->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus dPhiIn
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_dPhiIn_endcap2_separate->Fill(fabs(dPhiIn->at(i_ele)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus HoverE
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && HoverE->at(i_ele) < 0.13 // relaxed (Ele23Ele12 threshold)
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_HoverE_endcap2_separate->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus RelPFIso
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_PFiso_Rho_endcap2_separate->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus InvEminusInvP
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045 // Tighter by Ele23Ele12
                        && HoverE->at(i_ele) < 0.0878 // Tighter by Ele23Ele12
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && mHits->at(i_ele) <= 1
                        && passConvVeto->at(i_ele))
                        h_InvEminusInvP_endcap2_separate->Fill(InvEminusInvP->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus mHits
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && passConvVeto->at(i_ele))
                        h_mHits_endcap2_separate->Fill(mHits->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID minus passConvVeto
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298
                        && fabs(dEtaInSeed->at(i_ele)) < 0.00609
                        && fabs(dPhiIn->at(i_ele)) < 0.045
                        && HoverE->at(i_ele) < 0.0878
                        && relPFiso_Rho->at(i_ele) < 0.0821
                        && InvEminusInvP->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1)
                        h_passConvVeto_endcap2_separate->Fill(passConvVeto->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    // MediumID yes or no
                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.035
                        && HoverE->at(i_ele) < 0.13
                        && mHits->at(i_ele) <= 1)
                        h_passMediumID_endcap2_separate->Fill(passMediumID->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);

                    /// End test part ///

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap2_template[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (HoverE->at(i_ele) < 0.0878 && (fabs(InvEminusInvP->at(i_ele)) >= 0.13 || Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 ||
                             fabs(dEtaInSeed->at(i_ele)) >= 0.00609 || fabs(dPhiIn->at(i_ele)) >= 0.045 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap2_jetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                             HoverE->at(i_ele) < 0.0878 && fabs(InvEminusInvP->at(i_ele)) >= 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_PFiso_Rho_endcap2_badJetTemplate[ih]->Fill(relPFiso_Rho->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }

                    if (Full5x5_SigmaIEtaIEta->at(i_ele) < 0.0298 && fabs(dEtaInSeed->at(i_ele)) < 0.00609 && fabs(dPhiIn->at(i_ele)) < 0.045 &&
                        relPFiso_Rho->at(i_ele) < 0.0821 && fabs(InvEminusInvP->at(i_ele)) < 0.13 && passConvVeto->at(i_ele) && mHits->at(i_ele) <= 1)
                    {
                        h_HoverE_endcap_template_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap2_template[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                    else if (relPFiso_Rho->at(i_ele) < 0.0821 && (Full5x5_SigmaIEtaIEta->at(i_ele) >= 0.0298 || fabs(dEtaInSeed->at(i_ele)) >= 0.00609 ||
                             fabs(dPhiIn->at(i_ele)) >= 0.045 || fabs(InvEminusInvP->at(i_ele)) >= 0.13 || !passConvVeto->at(i_ele) || mHits->at(i_ele) > 1))
                    {
                        h_HoverE_endcap_jetTemplate_int->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        for (Int_t ih=0; ih<nPtBin_ele; ih++)
                        {
                            if (p_T->at(i_ele) > analyzer->ptbin_ele[ih] && p_T->at(i_ele) < analyzer->ptbin_ele[ih+1])
                                h_HoverE_endcap2_jetTemplate[ih]->Fill(HoverE->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_alt);
                        }
                    }
                }// End of if(far endcap)
            }// End of i_ele iteration

            if (fabs(ele_lead.Eta()) < 1.4442) // Barrel
            {
                h_MT_barrel_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                if (iso_lead < 0.15) h_MT_barrel_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                else h_MT_barrel_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
            }
            else if (fabs(ele_lead.Eta()) > 1.566 && fabs(ele_lead.Eta()) < 2.2) // Endcap
            {
                h_MT_endcap_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                if (iso_lead < 0.15) h_MT_endcap_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                else h_MT_endcap_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
            }
            else if (fabs(ele_lead.Eta()) >= 2.2 && fabs(ele_lead.Eta()) < 2.4) // Far endcap
            {
                h_MT_endcap2_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                if (iso_lead < 0.15) h_MT_endcap2_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
                else h_MT_endcap2_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * prescale_lead);
            }

        }// End of event iteration

        if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);

        f->cd();
        cout << "\tWriting into file...";

        h_pT_barrel_nume->Write();
        h_pT_endcap_nume->Write();
        h_pT_endcap2_nume->Write();
        h_pT_barrel_deno->Write();
        h_pT_endcap_deno->Write();
        h_pT_endcap2_deno->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_pT_endcap2_ctrl->Write();
        h_eta_nume->Write();
        h_eta_deno->Write();
        h_eta_ctrl->Write();
        h_PFiso_Rho_barrel_nume->Write();
        h_PFiso_Rho_endcap_nume->Write();
        h_PFiso_Rho_endcap2_nume->Write();
        h_PFiso_Rho_barrel_deno->Write();
        h_PFiso_Rho_endcap_deno->Write();
        h_PFiso_Rho_endcap2_deno->Write();
        h_PFiso_Rho_barrel_ctrl->Write();
        h_PFiso_Rho_endcap_ctrl->Write();
        h_PFiso_Rho_endcap2_ctrl->Write();
        h_SigmaIEtaIEta_barrel_nume->Write();
        h_SigmaIEtaIEta_endcap_nume->Write();
        h_SigmaIEtaIEta_endcap2_nume->Write();
        h_SigmaIEtaIEta_barrel_deno->Write();
        h_SigmaIEtaIEta_endcap_deno->Write();
        h_SigmaIEtaIEta_endcap2_deno->Write();
        h_SigmaIEtaIEta_barrel_ctrl->Write();
        h_SigmaIEtaIEta_endcap_ctrl->Write();
        h_SigmaIEtaIEta_endcap2_ctrl->Write();
        h_dEtaInSeed_barrel_nume->Write();
        h_dEtaInSeed_endcap_nume->Write();
        h_dEtaInSeed_endcap2_nume->Write();
        h_dEtaInSeed_barrel_deno->Write();
        h_dEtaInSeed_endcap_deno->Write();
        h_dEtaInSeed_endcap2_deno->Write();
        h_dEtaInSeed_barrel_ctrl->Write();
        h_dEtaInSeed_endcap_ctrl->Write();
        h_dEtaInSeed_endcap2_ctrl->Write();
        h_dPhiIn_barrel_nume->Write();
        h_dPhiIn_endcap_nume->Write();
        h_dPhiIn_endcap2_nume->Write();
        h_dPhiIn_barrel_deno->Write();
        h_dPhiIn_endcap_deno->Write();
        h_dPhiIn_endcap2_deno->Write();
        h_dPhiIn_barrel_ctrl->Write();
        h_dPhiIn_endcap_ctrl->Write();
        h_dPhiIn_endcap2_ctrl->Write();
        h_HoverE_barrel_nume->Write();
        h_HoverE_endcap_nume->Write();
        h_HoverE_endcap2_nume->Write();
        h_HoverE_barrel_deno->Write();
        h_HoverE_endcap_deno->Write();
        h_HoverE_endcap2_deno->Write();
        h_HoverE_barrel_ctrl->Write();
        h_HoverE_endcap_ctrl->Write();
        h_HoverE_endcap2_ctrl->Write();
        h_InvEminusInvP_barrel_nume->Write();
        h_InvEminusInvP_endcap_nume->Write();
        h_InvEminusInvP_endcap2_nume->Write();
        h_InvEminusInvP_barrel_deno->Write();
        h_InvEminusInvP_endcap_deno->Write();
        h_InvEminusInvP_endcap2_deno->Write();
        h_InvEminusInvP_barrel_ctrl->Write();
        h_InvEminusInvP_endcap_ctrl->Write();
        h_InvEminusInvP_endcap2_ctrl->Write();
        h_TrkIso_barrel_nume->Write();
        h_TrkIso_endcap_nume->Write();
        h_TrkIso_endcap2_nume->Write();
        h_TrkIso_barrel_deno->Write();
        h_TrkIso_endcap_deno->Write();
        h_TrkIso_endcap2_deno->Write();
        h_TrkIso_barrel_ctrl->Write();
        h_TrkIso_endcap_ctrl->Write();
        h_TrkIso_endcap2_ctrl->Write();
        h_ECALiso_barrel_nume->Write();
        h_ECALiso_endcap_nume->Write();
        h_ECALiso_endcap2_nume->Write();
        h_ECALiso_barrel_deno->Write();
        h_ECALiso_endcap_deno->Write();
        h_ECALiso_endcap2_deno->Write();
        h_ECALiso_barrel_ctrl->Write();
        h_ECALiso_endcap_ctrl->Write();
        h_ECALiso_endcap2_ctrl->Write();
        h_HCALiso_barrel_nume->Write();
        h_HCALiso_endcap_nume->Write();
        h_HCALiso_endcap2_nume->Write();
        h_HCALiso_barrel_deno->Write();
        h_HCALiso_endcap_deno->Write();
        h_HCALiso_endcap2_deno->Write();
        h_HCALiso_barrel_ctrl->Write();
        h_HCALiso_endcap_ctrl->Write();
        h_HCALiso_endcap2_ctrl->Write();
        h_MET->Write();
        h_MT_barrel_nume->Write();
        h_MT_endcap_nume->Write();
        h_MT_endcap2_nume->Write();
        h_MT_barrel_deno->Write();
        h_MT_endcap_deno->Write();
        h_MT_endcap2_deno->Write();
        h_MT_barrel_ctrl->Write();
        h_MT_endcap_ctrl->Write();
        h_MT_endcap2_ctrl->Write();
        h_nVTX->Write();

        h_HoverE_barrel_template_int->Write();
        h_HoverE_barrel_jetTemplate_int->Write();
        h_HoverE_endcap_template_int->Write();
        h_HoverE_endcap_jetTemplate_int->Write();
        h_HoverE_endcap2_template_int->Write();
        h_HoverE_endcap2_jetTemplate_int->Write();

        h_PFiso_Rho_barrel_separate->Write();
        h_PFiso_Rho_endcap_separate->Write();
        h_PFiso_Rho_endcap2_separate->Write();
        h_SigmaIEtaIEta_barrel_separate->Write();
        h_SigmaIEtaIEta_endcap_separate->Write();
        h_SigmaIEtaIEta_endcap2_separate->Write();
        h_dEtaInSeed_barrel_separate->Write();
        h_dEtaInSeed_endcap_separate->Write();
        h_dEtaInSeed_endcap2_separate->Write();
        h_dPhiIn_barrel_separate->Write();
        h_dPhiIn_endcap_separate->Write();
        h_dPhiIn_endcap2_separate->Write();
        h_HoverE_barrel_separate->Write();
        h_HoverE_endcap_separate->Write();
        h_HoverE_endcap2_separate->Write();
        h_InvEminusInvP_barrel_separate->Write();
        h_InvEminusInvP_endcap_separate->Write();
        h_InvEminusInvP_endcap2_separate->Write();
        h_mHits_barrel_separate->Write();
        h_mHits_endcap_separate->Write();
        h_mHits_endcap2_separate->Write();
        h_passConvVeto_barrel_separate->Write();
        h_passConvVeto_endcap_separate->Write();
        h_passConvVeto_endcap2_separate->Write();
        h_passMediumID_barrel_separate->Write();
        h_passMediumID_endcap_separate->Write();
        h_passMediumID_endcap2_separate->Write();

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

            h_PFiso_Rho_endcap2_template[ih]->Write();
            h_PFiso_Rho_endcap2_jetTemplate[ih]->Write();
            h_PFiso_Rho_endcap2_badJetTemplate[ih]->Write();
            h_HoverE_endcap2_template[ih]->Write();
            h_HoverE_endcap2_jetTemplate[ih]->Write();
        }

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
        TH1D* h_pT_endcap2_nume = new TH1D("h_pT_endcap2_nume", "h_pT_endcap2_nume", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno->Sumw2();
        TH1D* h_pT_endcap2_deno = new TH1D("h_pT_endcap2_deno", "h_pT_endcap2_deno", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_deno->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_pT_endcap2_ctrl = new TH1D("h_pT_endcap2_ctrl", "h_pT_endcap2_ctrl", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 96, -2.4, 2.4); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 96, -2.4, 2.4); h_eta_deno->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 96, -2.4, 2.4); h_eta_ctrl->Sumw2();
        TH1D* h_PFiso_barrel_nume = new TH1D("h_PFiso_barrel_nume", "h_PFiso_barrel_nume", 30, 0, 0.15); h_PFiso_barrel_nume->Sumw2();
        TH1D* h_PFiso_endcap_nume = new TH1D("h_PFiso_endcap_nume", "h_PFiso_endcap_nume", 30, 0, 0.15); h_PFiso_endcap_nume->Sumw2();
        TH1D* h_PFiso_endcap2_nume = new TH1D("h_PFiso_endcap2_nume", "h_PFiso_endcap2_nume", 30, 0, 0.15); h_PFiso_endcap2_nume->Sumw2();
        TH1D* h_PFiso_barrel_deno = new TH1D("h_PFiso_barrel_deno", "h_PFiso_barrel_deno", 100, 0, 5); h_PFiso_barrel_deno->Sumw2();
        TH1D* h_PFiso_endcap_deno = new TH1D("h_PFiso_endcap_deno", "h_PFiso_endcap_deno", 100, 0, 5); h_PFiso_endcap_deno->Sumw2();
        TH1D* h_PFiso_endcap2_deno = new TH1D("h_PFiso_endcap2_deno", "h_PFiso_endcap2_deno", 100, 0, 5); h_PFiso_endcap2_deno->Sumw2();
        TH1D* h_PFiso_barrel_ctrl = new TH1D("h_PFiso_barrel_ctrl", "h_PFiso_barrel_ctrl", 100, 0, 5); h_PFiso_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_endcap_ctrl = new TH1D("h_PFiso_endcap_ctrl", "h_PFiso_endcap_ctrl", 100, 0, 5); h_PFiso_endcap_ctrl->Sumw2();
        TH1D* h_PFiso_endcap2_ctrl = new TH1D("h_PFiso_endcap2_ctrl", "h_PFiso_endcap2_ctrl", 100, 0, 5); h_PFiso_endcap2_ctrl->Sumw2();
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
        TH1D* h_PFiso_endcap2_nume_50to70   = new TH1D("h_PFiso_endcap2_nume_50to70",   "h_PFiso_endcap2_nume_50to70",   15, 0, 0.15); h_PFiso_endcap2_nume_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_deno_50to70   = new TH1D("h_PFiso_barrel_deno_50to70",   "h_PFiso_barrel_deno_50to70",   50, 0, 5);    h_PFiso_barrel_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_deno_50to70   = new TH1D("h_PFiso_endcap_deno_50to70",   "h_PFiso_endcap_deno_50to70",   50, 0, 5);    h_PFiso_endcap_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap2_deno_50to70   = new TH1D("h_PFiso_endcap2_deno_50to70",   "h_PFiso_endcap2_deno_50to70",   50, 0, 5);    h_PFiso_endcap2_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_50to70   = new TH1D("h_PFiso_barrel_ctrl_50to70",   "h_PFiso_barrel_ctrl_50to70",   50, 0.15, 5);    h_PFiso_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_50to70   = new TH1D("h_PFiso_endcap_ctrl_50to70",   "h_PFiso_endcap_ctrl_50to70",   50, 0.15, 5);    h_PFiso_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap2_ctrl_50to70   = new TH1D("h_PFiso_endcap2_ctrl_50to70",   "h_PFiso_endcap2_ctrl_50to70",   50, 0.15, 5);    h_PFiso_endcap2_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_nume_70to100  = new TH1D("h_PFiso_barrel_nume_70to100",  "h_PFiso_barrel_nume_70to100",  15, 0, 0.15); h_PFiso_barrel_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_nume_70to100  = new TH1D("h_PFiso_endcap_nume_70to100",  "h_PFiso_endcap_nume_70to100",  15, 0, 0.15); h_PFiso_endcap_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap2_nume_70to100  = new TH1D("h_PFiso_endcap2_nume_70to100",  "h_PFiso_endcap2_nume_70to100",  15, 0, 0.15); h_PFiso_endcap2_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_deno_70to100  = new TH1D("h_PFiso_barrel_deno_70to100",  "h_PFiso_barrel_deno_70to100",  50, 0, 5);    h_PFiso_barrel_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_deno_70to100  = new TH1D("h_PFiso_endcap_deno_70to100",  "h_PFiso_endcap_deno_70to100",  50, 0, 5);    h_PFiso_endcap_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap2_deno_70to100  = new TH1D("h_PFiso_endcap2_deno_70to100",  "h_PFiso_endcap2_deno_70to100",  50, 0, 5);    h_PFiso_endcap2_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_70to100  = new TH1D("h_PFiso_barrel_ctrl_70to100",  "h_PFiso_barrel_ctrl_70to100",  50, 0.15, 5);    h_PFiso_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_70to100  = new TH1D("h_PFiso_endcap_ctrl_70to100",  "h_PFiso_endcap_ctrl_70to100",  50, 0.15, 5);    h_PFiso_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap2_ctrl_70to100  = new TH1D("h_PFiso_endcap2_ctrl_70to100",  "h_PFiso_endcap2_ctrl_70to100",  50, 0.15, 5);    h_PFiso_endcap2_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_nume_100to500 = new TH1D("h_PFiso_barrel_nume_100to500", "h_PFiso_barrel_nume_100to500", 15, 0, 0.15); h_PFiso_barrel_nume_100to500->Sumw2();
        TH1D* h_PFiso_endcap_nume_100to500 = new TH1D("h_PFiso_endcap_nume_100to500", "h_PFiso_endcap_nume_100to500", 15, 0, 0.15); h_PFiso_endcap_nume_100to500->Sumw2();
        TH1D* h_PFiso_endcap2_nume_100to500 = new TH1D("h_PFiso_endcap2_nume_100to500", "h_PFiso_endcap2_nume_100to500", 15, 0, 0.15); h_PFiso_endcap2_nume_100to500->Sumw2();
        TH1D* h_PFiso_barrel_deno_100to500 = new TH1D("h_PFiso_barrel_deno_100to500", "h_PFiso_barrel_deno_100to500", 50, 0, 5);    h_PFiso_barrel_deno_100to500->Sumw2();
        TH1D* h_PFiso_endcap_deno_100to500 = new TH1D("h_PFiso_endcap_deno_100to500", "h_PFiso_endcap_deno_100to500", 50, 0, 5);    h_PFiso_endcap_deno_100to500->Sumw2();
        TH1D* h_PFiso_endcap2_deno_100to500 = new TH1D("h_PFiso_endcap2_deno_100to500", "h_PFiso_endcap2_deno_100to500", 50, 0, 5);    h_PFiso_endcap2_deno_100to500->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_100to500 = new TH1D("h_PFiso_barrel_ctrl_100to500", "h_PFiso_barrel_ctrl_100to500", 50, 0.15, 5);    h_PFiso_barrel_ctrl_100to500->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_100to500 = new TH1D("h_PFiso_endcap_ctrl_100to500", "h_PFiso_endcap_ctrl_100to500", 50, 0.15, 5);    h_PFiso_endcap_ctrl_100to500->Sumw2();
        TH1D* h_PFiso_endcap2_ctrl_100to500 = new TH1D("h_PFiso_endcap2_ctrl_100to500", "h_PFiso_endcap2_ctrl_100to500", 50, 0.15, 5);    h_PFiso_endcap2_ctrl_100to500->Sumw2();

        TH1D* h_pT_barrel_nume_50to70   = new TH1D("h_pT_barrel_nume_50to70 ",  "h_pT_barrel_nume_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_50to70  ->Sumw2();
        TH1D* h_pT_endcap_nume_50to70   = new TH1D("h_pT_endcap_nume_50to70 ",  "h_pT_endcap_nume_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_50to70  ->Sumw2();
        TH1D* h_pT_endcap2_nume_50to70   = new TH1D("h_pT_endcap2_nume_50to70 ",  "h_pT_endcap2_nume_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_nume_50to70  ->Sumw2();
        TH1D* h_pT_barrel_deno_50to70   = new TH1D("h_pT_barrel_deno_50to70 ",  "h_pT_barrel_deno_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_50to70  ->Sumw2();
        TH1D* h_pT_endcap_deno_50to70   = new TH1D("h_pT_endcap_deno_50to70 ",  "h_pT_endcap_deno_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_50to70  ->Sumw2();
        TH1D* h_pT_endcap2_deno_50to70   = new TH1D("h_pT_endcap2_deno_50to70 ",  "h_pT_endcap2_deno_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_deno_50to70  ->Sumw2();
        TH1D* h_pT_barrel_ctrl_50to70   = new TH1D("h_pT_barrel_ctrl_50to70 ",  "h_pT_barrel_ctrl_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_endcap_ctrl_50to70   = new TH1D("h_pT_endcap_ctrl_50to70 ",  "h_pT_endcap_ctrl_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_endcap2_ctrl_50to70   = new TH1D("h_pT_endcap2_ctrl_50to70 ",  "h_pT_endcap2_ctrl_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_barrel_nume_70to100  = new TH1D("h_pT_barrel_nume_70to100",  "h_pT_barrel_nume_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_70to100 ->Sumw2();
        TH1D* h_pT_endcap_nume_70to100  = new TH1D("h_pT_endcap_nume_70to100",  "h_pT_endcap_nume_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_70to100 ->Sumw2();
        TH1D* h_pT_endcap2_nume_70to100  = new TH1D("h_pT_endcap2_nume_70to100",  "h_pT_endcap2_nume_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_nume_70to100 ->Sumw2();
        TH1D* h_pT_barrel_deno_70to100  = new TH1D("h_pT_barrel_deno_70to100",  "h_pT_barrel_deno_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_70to100 ->Sumw2();
        TH1D* h_pT_endcap_deno_70to100  = new TH1D("h_pT_endcap_deno_70to100",  "h_pT_endcap_deno_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_70to100 ->Sumw2();
        TH1D* h_pT_endcap2_deno_70to100  = new TH1D("h_pT_endcap2_deno_70to100",  "h_pT_endcap2_deno_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_deno_70to100 ->Sumw2();
        TH1D* h_pT_barrel_ctrl_70to100  = new TH1D("h_pT_barrel_ctrl_70to100",  "h_pT_barrel_ctrl_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_endcap_ctrl_70to100  = new TH1D("h_pT_endcap_ctrl_70to100",  "h_pT_endcap_ctrl_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_endcap2_ctrl_70to100  = new TH1D("h_pT_endcap2_ctrl_70to100",  "h_pT_endcap2_ctrl_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_barrel_nume_100to500 = new TH1D("h_pT_barrel_nume_100to500", "h_pT_barrel_nume_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_100to500->Sumw2();
        TH1D* h_pT_endcap_nume_100to500 = new TH1D("h_pT_endcap_nume_100to500", "h_pT_endcap_nume_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_100to500->Sumw2();
        TH1D* h_pT_endcap2_nume_100to500 = new TH1D("h_pT_endcap2_nume_100to500", "h_pT_endcap2_nume_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_nume_100to500->Sumw2();
        TH1D* h_pT_barrel_deno_100to500 = new TH1D("h_pT_barrel_deno_100to500", "h_pT_barrel_deno_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_100to500->Sumw2();
        TH1D* h_pT_endcap_deno_100to500 = new TH1D("h_pT_endcap_deno_100to500", "h_pT_endcap_deno_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_100to500->Sumw2();
        TH1D* h_pT_endcap2_deno_100to500 = new TH1D("h_pT_endcap2_deno_100to500", "h_pT_endcap2_deno_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_deno_100to500->Sumw2();
        TH1D* h_pT_barrel_ctrl_100to500 = new TH1D("h_pT_barrel_ctrl_100to500", "h_pT_barrel_ctrl_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_100to500->Sumw2();
        TH1D* h_pT_endcap_ctrl_100to500 = new TH1D("h_pT_endcap_ctrl_100to500", "h_pT_endcap_ctrl_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_100to500->Sumw2();
        TH1D* h_pT_endcap2_ctrl_100to500 = new TH1D("h_pT_endcap2_ctrl_100to500", "h_pT_endcap2_ctrl_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap2_ctrl_100to500->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Double_t MET_pT, MET_phi;
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
                    else if (fabs(eta->at(i_mu)) >= 1.2 && fabs(eta->at(i_mu)) < 1.8) // Endcap
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
                    else if (fabs(eta->at(i_mu)) >= 1.8 && fabs(eta->at(i_mu)) < 2.4) // Far endcap
                    {
                        h_pT_endcap2_nume->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap2_nume->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap2_nume_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_nume_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap2_nume_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_nume_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap2_nume_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_nume_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
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
                    else if (fabs(eta->at(i_mu)) >= 1.2 && fabs(eta->at(i_mu)) < 1.8) // Endcap
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
                    else if (fabs(eta->at(i_mu)) >= 1.8 && fabs(eta->at(i_mu)) < 2.4) // Far endcap
                    {
                        h_pT_endcap2_ctrl->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap2_ctrl->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap2_ctrl_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_ctrl_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap2_ctrl_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_ctrl_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap2_ctrl_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap2_ctrl_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
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
                else if (fabs(eta->at(i_mu)) >= 1.2 && fabs(eta->at(i_mu)) < 1.8) // Endcap
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
                else if (fabs(eta->at(i_mu)) >= 1.8 && fabs(eta->at(i_mu)) < 2.4) // Endcap
                {
                    h_pT_endcap2_deno->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_PFiso_endcap2_deno->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (p_T->at(i_mu) < 70)
                    {
                        h_pT_endcap2_deno_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap2_deno_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (p_T->at(i_mu) < 100)
                    {
                        h_pT_endcap2_deno_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap2_deno_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else //if (p_T->at(i_mu) < 500)
                    {
                        h_pT_endcap2_deno_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap2_deno_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
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
        h_pT_endcap2_nume->Write();
        h_pT_barrel_deno->Write();
        h_pT_endcap_deno->Write();
        h_pT_endcap2_deno->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_pT_endcap2_ctrl->Write();
        h_eta_nume->Write();
        h_eta_deno->Write();
        h_eta_ctrl->Write();
        h_PFiso_barrel_nume->Write();
        h_PFiso_endcap_nume->Write();
        h_PFiso_endcap2_nume->Write();
        h_PFiso_barrel_deno->Write();
        h_PFiso_endcap_deno->Write();
        h_PFiso_endcap2_deno->Write();
        h_PFiso_barrel_ctrl->Write();
        h_PFiso_endcap_ctrl->Write();
        h_PFiso_endcap2_ctrl->Write();
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
        h_pT_endcap2_nume_50to70->Write();
        h_pT_barrel_deno_50to70->Write();
        h_pT_endcap_deno_50to70->Write();
        h_pT_endcap2_deno_50to70->Write();
        h_pT_barrel_ctrl_50to70->Write();
        h_pT_endcap_ctrl_50to70->Write();
        h_pT_endcap2_ctrl_50to70->Write();
        h_pT_barrel_nume_70to100->Write();
        h_pT_endcap_nume_70to100->Write();
        h_pT_endcap2_nume_70to100->Write();
        h_pT_barrel_deno_70to100->Write();
        h_pT_endcap_deno_70to100->Write();
        h_pT_endcap2_deno_70to100->Write();
        h_pT_barrel_ctrl_70to100->Write();
        h_pT_endcap_ctrl_70to100->Write();
        h_pT_endcap2_ctrl_70to100->Write();
        h_pT_barrel_nume_100to500->Write();
        h_pT_endcap_nume_100to500->Write();
        h_pT_endcap2_nume_100to500->Write();
        h_pT_barrel_deno_100to500->Write();
        h_pT_endcap_deno_100to500->Write();
        h_pT_endcap2_deno_100to500->Write();
        h_pT_barrel_ctrl_100to500->Write();
        h_pT_endcap_ctrl_100to500->Write();
        h_pT_endcap2_ctrl_100to500->Write();

        h_PFiso_barrel_nume_50to70->Write();
        h_PFiso_endcap_nume_50to70->Write();
        h_PFiso_endcap2_nume_50to70->Write();
        h_PFiso_barrel_deno_50to70->Write();
        h_PFiso_endcap_deno_50to70->Write();
        h_PFiso_endcap2_deno_50to70->Write();
        h_PFiso_barrel_ctrl_50to70->Write();
        h_PFiso_endcap_ctrl_50to70->Write();
        h_PFiso_endcap2_ctrl_50to70->Write();
        h_PFiso_barrel_nume_70to100->Write();
        h_PFiso_endcap_nume_70to100->Write();
        h_PFiso_endcap2_nume_70to100->Write();
        h_PFiso_barrel_deno_70to100->Write();
        h_PFiso_endcap_deno_70to100->Write();
        h_PFiso_endcap2_deno_70to100->Write();
        h_PFiso_barrel_ctrl_70to100->Write();
        h_PFiso_endcap_ctrl_70to100->Write();
        h_PFiso_endcap2_ctrl_70to100->Write();
        h_PFiso_barrel_nume_100to500->Write();
        h_PFiso_endcap_nume_100to500->Write();
        h_PFiso_endcap2_nume_100to500->Write();
        h_PFiso_barrel_deno_100to500->Write();
        h_PFiso_endcap_deno_100to500->Write();
        h_PFiso_endcap2_deno_100to500->Write();
        h_PFiso_barrel_ctrl_100to500->Write();
        h_PFiso_endcap_ctrl_100to500->Write();
        h_PFiso_endcap2_ctrl_100to500->Write();

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
//    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract", 1);
    analyzer->SetupFRvalues_ele_fit(Dir+"FakeRate_electron.root");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.1);

    // Sign-missasignment correction
    TFile *f_ch = new TFile("~/DrellYan2016/SelectedX/etc/misid_majority.root");
    TH1D *h_misid = ((TH1D*)(f_ch->Get("data")));
    h_misid->SetDirectory(0);
    f_ch->Close();
    Double_t xAxis[13] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.44, 1.57, 1.8, 2, 2.2, 2.5};

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
        TH1D* h_mass_test = new TH1D("h_mass_test_"+Mgr.Procname[pr], "h_mass_test_"+Mgr.Procname[pr], binnum, massbins); h_mass_test->Sumw2();
        TH2D* h_zpeak = new TH2D("h_zpeak_"+Mgr.Procname[pr], "h_zpeak_"+Mgr.Procname[pr], nPtBin_ele, analyzer->ptbin_ele, 5, etabins2);
        TH1D* h_pT = new TH1D("h_pT_"+Mgr.Procname[pr], "h_pT_"+Mgr.Procname[pr], 100, 0, 1000); h_pT->Sumw2();
        TH1D* h_rapi = new TH1D("h_rapi_"+Mgr.Procname[pr], "h_rapi_"+Mgr.Procname[pr], 80, -4, 4); h_pT->Sumw2();

        TH1D* h_PFiso_Rho_barrel_ctrl = new TH1D("h_PFiso_Rho_barrel_ctrl_"+Mgr.Procname[pr], "h_PFiso_Rho_barrel_ctrl_"+Mgr.Procname[pr], 50, 0, 5); h_PFiso_Rho_barrel_ctrl->Sumw2();
        TH1D* h_TrkIso_barrel_ctrl = new TH1D("h_TrkIso_barrel_ctrl_"+Mgr.Procname[pr], "h_TrkIso_barrel_ctrl_"+Mgr.Procname[pr], 100, 0, 10); h_TrkIso_barrel_ctrl->Sumw2();
        TH1D* h_ECALiso_barrel_ctrl = new TH1D("h_ECALiso_barrel_ctrl_"+Mgr.Procname[pr], "h_ECALiso_barrel_ctrl_"+Mgr.Procname[pr], 100, 0, 10); h_ECALiso_barrel_ctrl->Sumw2();
        TH1D* h_HCALiso_barrel_ctrl = new TH1D("h_HCALiso_barrel_ctrl_"+Mgr.Procname[pr], "h_HCALiso_barrel_ctrl_"+Mgr.Procname[pr], 100, 0, 10); h_HCALiso_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_Rho_endcap_ctrl = new TH1D("h_PFiso_Rho_endcap_ctrl_"+Mgr.Procname[pr], "h_PFiso_Rho_endcap_ctrl_"+Mgr.Procname[pr], 50, 0, 5); h_PFiso_Rho_endcap_ctrl->Sumw2();
        TH1D* h_TrkIso_endcap_ctrl = new TH1D("h_TrkIso_endcap_ctrl_"+Mgr.Procname[pr], "h_TrkIso_endcap_ctrl_"+Mgr.Procname[pr], 20, 0, 2); h_TrkIso_endcap_ctrl->Sumw2();
        TH1D* h_ECALiso_endcap_ctrl = new TH1D("h_ECALiso_endcap_ctrl_"+Mgr.Procname[pr], "h_ECALiso_endcap_ctrl_"+Mgr.Procname[pr], 20, 0, 2); h_ECALiso_endcap_ctrl->Sumw2();
        TH1D* h_HCALiso_endcap_ctrl = new TH1D("h_HCALiso_endcap_ctrl_"+Mgr.Procname[pr], "h_HCALiso_endcap_ctrl_"+Mgr.Procname[pr], 20, 0, 2); h_HCALiso_endcap_ctrl->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<int> *scPixCharge = new std::vector<int>;
        std::vector<int> *isGsfCtfScPixConsistent = new std::vector<int>;
        std::vector<int> *isGsfScPixConsistent = new std::vector<int>;
        std::vector<int> *isGsfCtfConsistent = new std::vector<int>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<double> *relTrkIso = new std::vector<double>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
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
        chain->SetBranchStatus("scPixCharge", 1);
        chain->SetBranchStatus("isGsfCtfScPixConsistent", 1);
        chain->SetBranchStatus("isGsfScPixConsistent", 1);
        chain->SetBranchStatus("isGsfCtfConsistent", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("InvEminusInvP", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("relTrkIso", 1);
        chain->SetBranchStatus("relECALiso", 1);
        chain->SetBranchStatus("relHCALiso", 1);
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
        chain->SetBranchAddress("scPixCharge", &scPixCharge);
        chain->SetBranchAddress("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        chain->SetBranchAddress("isGsfScPixConsistent", &isGsfScPixConsistent);
        chain->SetBranchAddress("isGsfCtfConsistent", &isGsfCtfConsistent);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("relTrkIso", &relTrkIso);
        chain->SetBranchAddress("relECALiso", &relECALiso);
        chain->SetBranchAddress("relHCALiso", &relHCALiso);
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
        UInt_t nSameSign=0, nSameSignTight=0;
        UInt_t n_charge_tightCharge_match=0;
        UInt_t nEle=0;

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

            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (fabs(etaSC->at(0)) >= 2.4 || fabs(etaSC->at(1)) >= 2.4) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;
            if (passMediumID->at(0) == 1 || passMediumID->at(1)  == 1) continue;
            if (mHits->at(0) > 1 /*|| !passConvVeto->at(0)*/ || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 /*|| !passConvVeto->at(1)*/ || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;

            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;

//            // ALTERNATIVE FR (allowing only relPFiso to fail)        
//            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.00998 || fabs(dEtaInSeed->at(0)) >= 0.00311 ||
//                fabs(dPhiIn->at(0)) >= 0.07/*0.103*/ || HoverE->at(0) >= 0.13/*0.253*/ || InvEminusInvP->at(0) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.0298 || fabs(dEtaInSeed->at(0)) >= 0.00609 ||
//                fabs(dPhiIn->at(0)) >= 0.045 || HoverE->at(0) >= 0.0878 || InvEminusInvP->at(0) >= 0.13 ))
//                continue;
//            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.00998 || fabs(dEtaInSeed->at(1)) >= 0.00311 ||
//                fabs(dPhiIn->at(1)) >= 0.07/*0.103*/ || HoverE->at(1) >= 0.13/*0.253*/ || InvEminusInvP->at(1) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.0298 || fabs(dEtaInSeed->at(1)) >= 0.00609 ||
//                fabs(dPhiIn->at(1)) >= 0.045 || HoverE->at(1) >= 0.0878 || InvEminusInvP->at(1) >= 0.13))
//                continue;

            if (charge->at(0) == charge->at(1)) nSameSign++;
//            if (scPixCharge->at(0) == scPixCharge->at(1) && isGsfCtfScPixConsistent->at(0) && isGsfCtfScPixConsistent->at(1)) nSameSignTight++;
//            if (charge->at(0) == scPixCharge->at(0) && isGsfCtfScPixConsistent->at(0)) n_charge_tightCharge_match++;
//            if (charge->at(1) == scPixCharge->at(1) && isGsfCtfScPixConsistent->at(1)) n_charge_tightCharge_match++;
            nEle += 2;

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
            Double_t PT = (ele1+ele2).Pt();
            Double_t rapi = (ele1+ele2).Rapidity();

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
            FR1 = analyzer->FakeRate_ele_fit(p_T->at(0), etaSC->at(0));
            FR2 = analyzer->FakeRate_ele_fit(p_T->at(1), etaSC->at(1));
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            Double_t nFails1=0, nFails2=0;
            if (fabs(etaSC->at(0)) < 1.4442)
            {
                if (Full5x5_SigmaIEtaIEta->at(0) >= 0.00998) nFails1++;
                if (fabs(dEtaInSeed->at(0)) >= 0.00311) nFails1++;
                if (fabs(dPhiIn->at(0)) >= 0.103) nFails1++;
                if (HoverE->at(0) >= 0.253) nFails1++;
                if (relPFiso_Rho->at(0) >= 0.0695) nFails1++;
                if (InvEminusInvP->at(0) >= 0.134) nFails1++;
                if (mHits->at(0) > 1) nFails1++;
                if (passConvVeto->at(0) == 0) nFails1++;
            }
            else if (fabs(etaSC->at(0)) > 1.566)
            {
                if (Full5x5_SigmaIEtaIEta->at(0) >= 0.0298) nFails1++;
                if (fabs(dEtaInSeed->at(0)) >= 0.00609) nFails1++;
                if (fabs(dPhiIn->at(0)) >= 0.045) nFails1++;
                if (HoverE->at(0) >= 0.0878) nFails1++;
                if (relPFiso_Rho->at(0) >= 0.0821) nFails1++;
                if (InvEminusInvP->at(0) >= 0.13) nFails1++;
                if (mHits->at(0) > 1) nFails1++;
                if (passConvVeto->at(0) == 0) nFails1++;
            }
//            if (fabs(etaSC->at(1)) < 1.4442)
//            {
//                if (Full5x5_SigmaIEtaIEta->at(1) >= 0.00998) nFails2++;
//                if (fabs(dEtaInSeed->at(1)) >= 0.00311) nFails2++;
//                if (fabs(dPhiIn->at(1)) >= 0.103) nFails2++;
//                if (HoverE->at(1) >= 0.253) nFails2++;
//                if (relPFiso_Rho->at(1) >= 0.0695) nFails2++;
//                if (InvEminusInvP->at(1) >= 0.134) nFails2++;
//                if (mHits->at(1) > 1) nFails2++;
//                if (passConvVeto->at(1) == 0) nFails2++;
//            }
//            else if (fabs(etaSC->at(1)) > 1.566)
//            {
//                if (Full5x5_SigmaIEtaIEta->at(1) >= 0.0298) nFails2++;
//                if (fabs(dEtaInSeed->at(1)) >= 0.00609) nFails2++;
//                if (fabs(dPhiIn->at(1)) >= 0.045) nFails2++;
//                if (HoverE->at(1) >= 0.0878) nFails2++;
//                if (relPFiso_Rho->at(1) >= 0.0821) nFails2++;
//                if (InvEminusInvP->at(1) >= 0.13) nFails2++;
//                if (mHits->at(1) > 1) nFails2++;
//                if (passConvVeto->at(1) == 0) nFails2++;
//            }

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass_test->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_pT->Fill(PT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_rapi->Fill(rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (charge->at(0) == charge->at(1))
            {
//                Double_t f1=0, f2=0, weight=0;
//                for (Int_t i_eta=0; i_eta<12; i_eta++)
//                {
//                    if (fabs(eta->at(0)) >= xAxis[i_eta] && fabs(eta->at(0)) < xAxis[i_eta+1])
//                        f1 = 1 - h_misid->GetBinContent(i_eta+1);
//                    if (fabs(eta->at(1)) >= xAxis[i_eta] && fabs(eta->at(1)) < xAxis[i_eta+1])
//                        f2 = 1 - h_misid->GetBinContent(i_eta+1);
//                }
//                f1 = 0.01; f2 = 0.01;
//                Double_t part1 = (f1 * (1-f2) + f2 * (1-f1)) * (f1 * (1-f2) + f2 * (1-f1)) / ((1-f1) * (1-f2) * (1-f1) * (1-f2));
//                weight = 1 / ((1-f1) * (1-f2)) + part1;
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight/* * weight*/);
                if (scPixCharge->at(0) == charge->at(0) && scPixCharge->at(1) == charge->at(1) && isGsfCtfScPixConsistent->at(0) && isGsfCtfScPixConsistent->at(1))
                    h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight/* * weight*/);
            }
//            else
//            {
//                Double_t f1=0, f2=0, weight=0;
//                for (Int_t i_eta=0; i_eta<12; i_eta++)
//                {
//                    if (fabs(eta->at(0)) >= xAxis[i_eta] && fabs(eta->at(0)) < xAxis[i_eta+1])
//                        f1 = 1 - h_misid->GetBinContent(i_eta+1);
//                    if (fabs(eta->at(1)) >= xAxis[i_eta] && fabs(eta->at(1)) < xAxis[i_eta+1])
//                        f2 = 1 - h_misid->GetBinContent(i_eta+1);
//                }
//                f1 = 0.01; f2 = 0.01;
//                weight = 0 - (f1 * (1-f2) + f2 * (1-f1)) / ((1-f1) * (1-f2) * (1-f1) * (1-f2));
//                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight * weight);
//            }
            if (mass < 60)
            {
                h_zpeak->Fill(p_T->at(0), fabs(eta->at(0)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_zpeak->Fill(p_T->at(1), fabs(eta->at(1)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }
            if (fabs(eta->at(0)) < 1.4442)
            {
                h_PFiso_Rho_barrel_ctrl->Fill(relPFiso_Rho->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_TrkIso_barrel_ctrl->Fill(relTrkIso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_ECALiso_barrel_ctrl->Fill(relECALiso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_HCALiso_barrel_ctrl->Fill(relHCALiso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else if (fabs(eta->at(0)) > 1.566)
            {
                h_PFiso_Rho_endcap_ctrl->Fill(relPFiso_Rho->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_TrkIso_endcap_ctrl->Fill(relTrkIso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_ECALiso_endcap_ctrl->Fill(relECALiso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_HCALiso_endcap_ctrl->Fill(relHCALiso->at(0), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            if (fabs(eta->at(1)) < 1.4442)
            {
                h_PFiso_Rho_barrel_ctrl->Fill(relPFiso_Rho->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_TrkIso_barrel_ctrl->Fill(relTrkIso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_ECALiso_barrel_ctrl->Fill(relECALiso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_HCALiso_barrel_ctrl->Fill(relHCALiso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else if (fabs(eta->at(1)) > 1.566)
            {
                h_PFiso_Rho_endcap_ctrl->Fill(relPFiso_Rho->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_TrkIso_endcap_ctrl->Fill(relTrkIso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_ECALiso_endcap_ctrl->Fill(relECALiso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_HCALiso_endcap_ctrl->Fill(relHCALiso->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
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
        cout << "\t # passed electrons: " << nPassEle << endl;
        cout << "\t # failed electrons: " << nFailEle << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_SS->Write();
        h_mass_temp->Write();
        h_mass_test->Write();
        h_nVTX->Write();
        h_zpeak->Write();
        h_pT->Write();
        h_rapi->Write();
        h_PFiso_Rho_barrel_ctrl->Write();
        h_TrkIso_barrel_ctrl->Write();
        h_ECALiso_barrel_ctrl->Write();
        h_HCALiso_barrel_ctrl->Write();
        h_PFiso_Rho_endcap_ctrl->Write();
        h_TrkIso_endcap_ctrl->Write();
        h_ECALiso_endcap_ctrl->Write();
        h_HCALiso_endcap_ctrl->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

        cout << "nSameSign = " << nSameSign << endl;
        cout << "nSameSignTight = " << nSameSignTight << endl;
        Double_t Ratio = (Double_t)nSameSignTight/(Double_t)nSameSign;
        cout << "Ratio = " << Ratio << endl;
        cout << "Matches between charge and tight charge: " << n_charge_tightCharge_match << endl;
        Double_t rMatches = (Double_t)n_charge_tightCharge_match/(Double_t)nEle;
        cout << "Ratio = " << rMatches << endl;

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

    // -- TEST!! DY efficiency -- //
//    analyzer->SetupDYeff();

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For W+Jets estimation from Fake Rate -- //
//    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract", 1);
    analyzer->SetupFRvalues_ele_fit(Dir+"FakeRate_electron.root");

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.5);
    TH1D* h_DYweight = new TH1D("h_DYweight", "DY eff weights", 100, 0.5, 1.5);

    // Sign-missasignment correction
    TFile *f_ch = new TFile("~/DrellYan2016/SelectedX/etc/misid_majority.root");
    TH1D *h_misid = ((TH1D*)(f_ch->Get("data")));
    h_misid->SetDirectory(0);
    f_ch->Close();
    Double_t xAxis[13] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.44, 1.57, 1.8, 2, 2.2, 2.5};

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
        TH1D* h_mass_test = new TH1D("h_mass_test_"+Mgr.Procname[pr], "h_mass_test_"+Mgr.Procname[pr], binnum, massbins); h_mass_test->Sumw2();
        TH2D* h_zpeak = new TH2D("h_zpeak_"+Mgr.Procname[pr], "h_zpeak_"+Mgr.Procname[pr], nPtBin_ele, analyzer->ptbin_ele, 5, etabins2);
        TH1D* h_pT = new TH1D("h_pT_"+Mgr.Procname[pr], "h_pT_"+Mgr.Procname[pr], 100, 0, 1000); h_pT->Sumw2();
        TH1D* h_rapi = new TH1D("h_rapi_"+Mgr.Procname[pr], "h_rapi_"+Mgr.Procname[pr], 80, -4, 4); h_rapi->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<int> *scPixCharge = new std::vector<int>;
        std::vector<int> *isGsfCtfScPixConsistent = new std::vector<int>;
        std::vector<int> *isGsfScPixConsistent = new std::vector<int>;
        std::vector<int> *isGsfCtfConsistent = new std::vector<int>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *InvEminusInvP = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<double> *relTrkIso = new std::vector<double>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG) cout << Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("etaSC", 1);
        chain->SetBranchStatus("charge", 1);
        chain->SetBranchStatus("scPixCharge", 1);
        chain->SetBranchStatus("isGsfCtfScPixConsistent", 1);
        chain->SetBranchStatus("isGsfScPixConsistent", 1);
        chain->SetBranchStatus("isGsfCtfConsistent", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("InvEminusInvP", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("relTrkIso", 1);
        chain->SetBranchStatus("relECALiso", 1);
        chain->SetBranchStatus("relHCALiso", 1);
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
        chain->SetBranchAddress("scPixCharge", &scPixCharge);
        chain->SetBranchAddress("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        chain->SetBranchAddress("isGsfScPixConsistent", &isGsfScPixConsistent);
        chain->SetBranchAddress("isGsfCtfConsistent", &isGsfCtfConsistent);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("relTrkIso", &relTrkIso);
        chain->SetBranchAddress("relECALiso", &relECALiso);
        chain->SetBranchAddress("relHCALiso", &relHCALiso);
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
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || (fabs(etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;
            if (fabs(etaSC->at(0)) >= 2.4 || fabs(etaSC->at(0)) >= 2.4) continue;
            if (passMediumID->at(0) == 1 && passMediumID->at(1) == 1) continue;
            if (passMediumID->at(0) == 0 && passMediumID->at(1) == 0) continue;

            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;

            if (mHits->at(0) > 1 /*|| !passConvVeto->at(0)*/ || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 /*|| !passConvVeto->at(1)*/ || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;

            // ALTERNATIVE FR (allowing only relPFiso to fail)
//            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.00998 || fabs(dEtaInSeed->at(0)) >= 0.00311 ||
//                fabs(dPhiIn->at(0)) >= 0.07/*0.103*/ || HoverE->at(0) >= 0.13/*0.253*/ || InvEminusInvP->at(0) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.0298 || fabs(dEtaInSeed->at(0)) >= 0.00609 ||
//                fabs(dPhiIn->at(0)) >= 0.045 || HoverE->at(0) >= 0.0878 || InvEminusInvP->at(0) >= 0.13))
//                continue;
//            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.00998 || fabs(dEtaInSeed->at(1)) >= 0.00311 ||
//                fabs(dPhiIn->at(1)) >= 0.07/*0.103*/ || HoverE->at(1) >= 0.13/*0.253*/ || InvEminusInvP->at(1) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.0298 || fabs(dEtaInSeed->at(1)) >= 0.00609 ||
//                fabs(dPhiIn->at(1)) >= 0.045 || HoverE->at(1) >= 0.0878 || InvEminusInvP->at(1) >= 0.13))
//                continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;

            nPass++;

            if (DEBUG)
            {
                if (nPass >= 1000) break;
                cout << "Evt " << i << endl;
                cout << "nElectrons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tcharge[0] = " << charge->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tcharge[1] = " << charge->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "\nMET = " << MET_pT << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele1, ele2, ele1_SF, ele2_SF;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            ele1_SF.SetPtEtaPhiM(p_T->at(0), etaSC->at(0), phi->at(0), M_Elec);
            ele2_SF.SetPtEtaPhiM(p_T->at(1), etaSC->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();
            Double_t PT = (ele1+ele2).Pt();
            Double_t rapi = (ele1+ele2).Rapidity();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            // -- TEST!!! DY eff weights -- //
            Double_t DYeffWeight = 1;
//            if (Mgr.isMC && !Mgr.Tag[0].Contains("QCD") && !Mgr.Tag[0].Contains("Jet"))
//            {
//                DYeffWeight = analyzer->DYeff_evtWeight(p_T->at(0), eta->at(0), passMediumID->at(0), p_T->at(1), eta->at(1), passMediumID->at(1));
//                h_DYweight->Fill(DYeffWeight);
//            }

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << "\tDY eff weight (TEST): " << DYeffWeight << endl;

            // -- FR (and SF) WEIGHTS -- //
            Int_t nFails = 0;
            Double_t FRweight = 1;
            Double_t FR;
            if (passMediumID->at(0) == 0 && passMediumID->at(1) == 1) // First fails, second passes
            {
                FR = analyzer->FakeRate_ele_fit(p_T->at(0), etaSC->at(0));
                if (Mgr.isMC == kTRUE)
                {
                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 2);
                }
                if (fabs(etaSC->at(0)) < 1.4442)
                {
                    if (Full5x5_SigmaIEtaIEta->at(0) >= 0.00998) nFails++;
                    if (fabs(dEtaInSeed->at(0)) >= 0.00311) nFails++;
                    if (fabs(dPhiIn->at(0)) >= 0.103) nFails++;
                    if (HoverE->at(0) >= 0.253) nFails++;
                    if (relPFiso_Rho->at(0) >= 0.0695) nFails++;
                    if (InvEminusInvP->at(0) >= 0.134) nFails++;
                    if (mHits->at(0) > 1) nFails++;
                    if (passConvVeto->at(0) == 0) nFails++;
                }
                else if (fabs(etaSC->at(0)) > 1.566)
                {
                    if (Full5x5_SigmaIEtaIEta->at(0) >= 0.0298) nFails++;
                    if (fabs(dEtaInSeed->at(0)) >= 0.00609) nFails++;
                    if (fabs(dPhiIn->at(0)) >= 0.045) nFails++;
                    if (HoverE->at(0) >= 0.0878) nFails++;
                    if (relPFiso_Rho->at(0) >= 0.0821) nFails++;
                    if (InvEminusInvP->at(0) >= 0.13) nFails++;
                    if (mHits->at(0) > 1) nFails++;
                    if (passConvVeto->at(0) == 0) nFails++;
                }
            }
            else // Second fails, first passes
            {
                FR = analyzer->FakeRate_ele_fit(p_T->at(1), etaSC->at(1));
                if (Mgr.isMC == kTRUE)
                {
                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 1);
                }

                if (fabs(etaSC->at(1)) < 1.4442)
                {
                    if (Full5x5_SigmaIEtaIEta->at(1) >= 0.00998) nFails++;
                    if (fabs(dEtaInSeed->at(1)) >= 0.00311) nFails++;
                    if (fabs(dPhiIn->at(1)) >= 0.103) nFails++;
                    if (HoverE->at(1) >= 0.253) nFails++;
                    if (relPFiso_Rho->at(1) >= 0.0695) nFails++;
                    if (InvEminusInvP->at(1) >= 0.134) nFails++;
                    if (mHits->at(1) > 1) nFails++;
                    if (passConvVeto->at(1) == 0) nFails++;
                }
                else if (fabs(etaSC->at(1)) > 1.566)
                {
                    if (Full5x5_SigmaIEtaIEta->at(1) >= 0.0298) nFails++;
                    if (fabs(dEtaInSeed->at(1)) >= 0.00609) nFails++;
                    if (fabs(dPhiIn->at(1)) >= 0.045) nFails++;
                    if (HoverE->at(1) >= 0.0878) nFails++;
                    if (relPFiso_Rho->at(1) >= 0.0821) nFails++;
                    if (InvEminusInvP->at(1) >= 0.13) nFails++;
                    if (mHits->at(1) > 1) nFails++;
                    if (passConvVeto->at(1) == 0) nFails++;
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
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
            h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
            h_mass_test->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * DYeffWeight * TopPtWeight);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
            h_pT->Fill(PT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
            h_rapi->Fill(rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);

            if (charge->at(0) == charge->at(1))
            {
//                Double_t f1=0, f2=0, weight=0;
//                for (Int_t i_eta=0; i_eta<12; i_eta++)
//                {
//                    if (fabs(eta->at(0)) >= xAxis[i_eta] && fabs(eta->at(0)) < xAxis[i_eta+1])
//                        f1 = 1 - h_misid->GetBinContent(i_eta+1);
//                    if (fabs(eta->at(1)) >= xAxis[i_eta] && fabs(eta->at(1)) < xAxis[i_eta+1])
//                        f2 = 1 - h_misid->GetBinContent(i_eta+1);
//                }
//                Double_t part1 = (f1 * (1-f2) + f2 * (1-f1)) * (f1 * (1-f2) + f2 * (1-f1)) / ((1-f1) * (1-f2) * (1-f1) * (1-f2));
//                weight = 1 / ((1-f1) * (1-f2)) + part1;
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight/* * weight*/);
                if (scPixCharge->at(0) == charge->at(0) && scPixCharge->at(1) == charge->at(1) && isGsfCtfScPixConsistent->at(0) && isGsfCtfScPixConsistent->at(1))
                    h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight/* * weight*/);
            }
//            else
//            {
//                Double_t f1=0, f2=0, weight=0;
//                for (Int_t i_eta=0; i_eta<12; i_eta++)
//                {
//                    if (fabs(eta->at(0)) >= xAxis[i_eta] && fabs(eta->at(0)) < xAxis[i_eta+1])
//                        f1 = 1 - h_misid->GetBinContent(i_eta+1);
//                    if (fabs(eta->at(1)) >= xAxis[i_eta] && fabs(eta->at(1)) < xAxis[i_eta+1])
//                        f2 = 1 - h_misid->GetBinContent(i_eta+1);
//                }
//                weight = 0 - (f1 * (1-f2) + f2 * (1-f1)) / ((1-f1) * (1-f2) * (1-f1) * (1-f2));
//                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight * weight);
//            }
            if (mass > 75 && mass < 105)
            {
                if (passMediumID->at(0) == 0 && passMediumID->at(1) == 1)
                    h_zpeak->Fill(p_T->at(0), fabs(eta->at(0)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
                if (passMediumID->at(0) == 1 && passMediumID->at(1) == 0)
                    h_zpeak->Fill(p_T->at(1), fabs(eta->at(1)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * DYeffWeight * FRweight);
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
        h_mass_SS->Write();
        h_mass_temp->Write();
        h_mass_test->Write();
        h_MET->Write();
        h_nVTX->Write();
        h_zpeak->Write();
        h_pT->Write();
        h_rapi->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    TCanvas *c_DYweight = new TCanvas("DYweight","DY weights", 800, 800);
    h_DYweight->Draw();
    c_DYweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_E"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_E"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_WJET_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_QCD_HistMaker (Bool_t DEBUG)
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
        TH1D* h_pT_lead = new TH1D("h_pT_lead_"+Mgr.Procname[pr], "h_pT_lead_"+Mgr.Procname[pr], 175, 25, 200); h_pT_lead->Sumw2();
        TH1D* h_pT_sublead = new TH1D("h_pT_sublead_"+Mgr.Procname[pr], "h_pT_sublead_"+Mgr.Procname[pr], 85, 15, 100); h_pT_sublead->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT_"+Mgr.Procname[pr], "h_HLT_pT_"+Mgr.Procname[pr], 100, 0, 200); h_HLT_pT->Sumw2();
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
        std::vector<string> *trig_name = new std::vector<string>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Double_t MET_pT, MET_Phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);        
        chain->SetBranchStatus("trig_name", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("prescale_factor", 1);
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
        chain->SetBranchAddress("trig_name", &trig_name);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
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
            if (p_T->at(0) < 52 && p_T->at(1) < 52) continue;
//            if (p_T->at(0) < 2 || p_T->at(1) < 2) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            nPass++;
            Int_t nTrig = (Int_t)trig_fired->size();

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
                if (nTrig)
                {
                    for (Int_t i_tr=0; i_tr<nTrig; i_tr++)
                    {
                        if (trig_matched->at(i_tr) == 0)
                        {
                            h_HLT_pT->Fill(trig_pT->at(i_tr), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            break;
                        }
                    }
                }
            }
            else
            {
                h2_pT->Fill(mu2.Pt(), mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_lead->Fill(mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_sublead->Fill(mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                if (nTrig)
                {
                    for (Int_t i_tr=0; i_tr<nTrig; i_tr++)
                    {
                        if (trig_matched->at(i_tr) == 1)
                        {
                            h_HLT_pT->Fill(trig_pT->at(i_tr), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            break;
                        }
                    }
                }
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
        h_HLT_pT->Write();
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


void Mu_WJET_HistMaker (Bool_t DEBUG)
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
        chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

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
        TH1D* h_mass_temp = new TH1D("h_mass_template_"+Mgr.Procname[pr], "h_mass_template_"+Mgr.Procname[pr], binnum, massbins); h_mass_temp->Sumw2();
        TH1D* h_mass_test = new TH1D("h_mass_test_"+Mgr.Procname[pr], "h_mass_test_"+Mgr.Procname[pr], binnum, massbins); h_mass_test->Sumw2();
        TH1D* h_pT_e = new TH1D("h_pT_e_"+Mgr.Procname[pr], "h_pT_e_"+Mgr.Procname[pr], 50, 0, 1000); h_pT_e->Sumw2();
        TH1D* h_pT_mu = new TH1D("h_pT_mu_"+Mgr.Procname[pr], "h_pT_mu_"+Mgr.Procname[pr], 50, 0, 1000); h_pT_mu->Sumw2();

        Double_t e_p_T, e_eta, e_etaSC, e_phi;
        Int_t e_charge;
        Double_t e_Full5x5_SigmaIEtaIEta, e_dEtaInSeed, e_dPhiIn, e_HoverE, e_InvEminusInvP;
        Double_t e_chIso03, e_nhIso03, e_phIso03, e_ChIso03FromPU;
        Double_t e_relPFiso_dBeta, e_relPFiso_Rho;
        Int_t e_mHits;
        Int_t e_passConvVeto, e_passMediumID;
        Double_t mu_p_T, mu_eta, mu_phi;
        Int_t mu_charge;
        Double_t mu_relPFiso_dBeta;
        Int_t nPU, nVTX;
        Double_t PVz, gen_weight, top_weight, prefiring_weight, prefiring_weight_up, prefiring_weight_down;

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

            if (e_passMediumID) nPassEle++;
            else nFailEle++;
            if (mu_relPFiso_dBeta < 0.15) nPassMu++;
            else nFailMu++;

            // QCD selection
            if (e_passMediumID || mu_relPFiso_dBeta < 0.15) continue;
            if (e_p_T < 17 || mu_p_T < 17) continue;
            if (e_p_T < 28 && mu_p_T < 28) continue;
            if (fabs(e_etaSC) > 1.4442 && fabs(e_etaSC) < 1.566) continue;
            if (fabs(e_etaSC) < 1.4442 && (e_Full5x5_SigmaIEtaIEta >= 0.013 || e_HoverE >= 0.13 ||
                                           fabs(e_dEtaInSeed) >= 0.01 || fabs(e_dPhiIn) >= 0.07 || e_mHits > 1)) continue;
            if (fabs(e_etaSC) > 1.566 && (e_Full5x5_SigmaIEtaIEta >= 0.035 || e_HoverE >= 0.13 || e_mHits > 1)) continue;

            if (e_p_T != e_p_T) cout << e_p_T << " " << e_eta << " " << e_phi << " " << e_charge << " " << e_relPFiso_Rho << endl;
            if (mu_p_T != mu_p_T) cout << mu_p_T << " " << mu_eta << " " << mu_phi << " " << mu_charge << " " << mu_relPFiso_dBeta << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T;
                cout << "\teta = " << e_eta;
                cout << "\tphi = " << e_phi << endl;
                cout << "\tpassMediumID = " << e_passMediumID << endl;
                cout << "Muon p_T = " << mu_p_T;
                cout << "\teta = " << mu_eta;
                cout << "\tphi = " << mu_phi << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu, ele_SF;
            ele.SetPtEtaPhiM(e_p_T, e_eta, e_phi, M_Elec);
            mu.SetPtEtaPhiM(mu_p_T, mu_eta, mu_phi, M_Mu);
            ele_SF.SetPtEtaPhiM(e_p_T, e_etaSC, e_phi, M_Elec);
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
            FR1 = analyzer->FakeRate_ele(e_p_T, e_etaSC);
            FR2 = analyzer->FakeRate(mu_p_T, mu_eta);
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (e_charge != mu_charge)
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            else
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass_test->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_pT_e->Fill(e_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_pT_mu->Fill(mu_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

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
        h_mass_temp->Write();
        h_mass_test->Write();
        h_pT_e->Write();
        h_pT_mu->Write();

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

    Int_t n_efail=0, n_mufail=0;
    Double_t nw_efail=0, nw_mufail=0;

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        if (DEBUG == kTRUE && pr != _WJets_ext2v5) continue;
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
        TH1D* h_mass_test = new TH1D("h_mass_test_"+Mgr.Procname[pr], "h_mass_test_"+Mgr.Procname[pr], binnum, massbins); h_mass_test->Sumw2();
        TH1D* h_mass_test_ef = new TH1D("h_mass_test_elefail_"+Mgr.Procname[pr], "h_mass_test_elefail_"+Mgr.Procname[pr], binnum, massbins); h_mass_test_ef->Sumw2();
        TH1D* h_mass_test_mf = new TH1D("h_mass_test_mufail_"+Mgr.Procname[pr], "h_mass_test_mufail_"+Mgr.Procname[pr], binnum, massbins); h_mass_test_mf->Sumw2();
        TH1D* h_pT_e = new TH1D("h_pT_e_"+Mgr.Procname[pr], "h_pT_e_"+Mgr.Procname[pr], 50, 0, 1000); h_pT_e->Sumw2();
        TH1D* h_pT_mu = new TH1D("h_pT_mu_"+Mgr.Procname[pr], "h_pT_mu_"+Mgr.Procname[pr], 50, 0, 1000); h_pT_mu->Sumw2();

        Double_t e_p_T, e_eta, e_etaSC, e_phi;
        Int_t e_charge;
        Double_t e_Full5x5_SigmaIEtaIEta, e_dEtaInSeed, e_dPhiIn, e_HoverE, e_InvEminusInvP;
        Double_t e_chIso03, e_nhIso03, e_phIso03, e_ChIso03FromPU;
        Double_t e_relPFiso_dBeta, e_relPFiso_Rho;
        Int_t e_mHits;
        Int_t e_passConvVeto, e_passMediumID;
        Double_t mu_p_T, mu_eta, mu_phi;
        Int_t mu_charge;
        Double_t mu_relPFiso_dBeta;
        Int_t nPU, nVTX;
        Double_t PVz, gen_weight, top_weight, prefiring_weight, prefiring_weight_up, prefiring_weight_down;

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
            if (e_passMediumID == 1 && mu_relPFiso_dBeta < 0.15) continue;
            if (e_passMediumID == 0 && mu_relPFiso_dBeta >= 0.15) continue;
            if (e_p_T < 17 || mu_p_T < 17) continue;
            if (e_p_T < 28 && mu_p_T < 28) continue;
            if (fabs(e_etaSC) > 1.4442 && fabs(e_etaSC) < 1.566) continue;
            if (fabs(e_etaSC) < 1.4442 && (e_Full5x5_SigmaIEtaIEta >= 0.013 || e_HoverE >= 0.13 ||
                                           fabs(e_dEtaInSeed) >= 0.01 || fabs(e_dPhiIn) >= 0.07 || e_mHits > 1)) continue;
            if (fabs(e_etaSC) > 1.566 && (e_Full5x5_SigmaIEtaIEta >= 0.035 || e_HoverE >= 0.13 || e_mHits > 1)) continue;

            if (e_p_T != e_p_T) cout << e_p_T << " " << e_eta << " " << e_phi << " " << e_charge << " " << e_relPFiso_Rho << endl;
            if (mu_p_T != mu_p_T) cout << mu_p_T << " " << mu_eta << " " << mu_phi << " " << mu_charge << " " << mu_relPFiso_dBeta << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 50) break;
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T;
                cout << "\teta = " << e_eta;
                cout << "\tphi = " << e_phi << endl;
                cout << "\tpassMediumID = " << e_passMediumID << endl;
                cout << "Muon p_T = " << mu_p_T;
                cout << "\teta = " << mu_eta;
                cout << "\tphi = " << mu_phi << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu, ele_SF;
            ele.SetPtEtaPhiM(e_p_T, e_eta, e_phi, M_Elec);
            mu.SetPtEtaPhiM(mu_p_T, mu_eta, mu_phi, M_Mu);
            ele_SF.SetPtEtaPhiM(e_p_T, e_etaSC, e_phi, M_Elec);
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
            if (e_passMediumID == 0 && mu_relPFiso_dBeta < 0.15) // Electron fails, muon passes
            {                
                FR = analyzer->FakeRate_ele(e_p_T, e_etaSC);
                FRweight = FR / (1 - FR);
//                if (Mgr.isMC == kTRUE)
//                {
//                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 2);
//                }
                if (e_charge != mu_charge)
                    h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                else
                    h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                if (Mgr.isMC == kFALSE)
                {
                    n_efail++;
                    nw_efail += FRweight;
                }
                h_mass_test_ef->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_pT_mu->Fill(mu_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                h_pT_e->Fill(e_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else if (e_passMediumID == 1 && mu_relPFiso_dBeta >= 0.15) // Muon fails, electron passes
            {                
                FR = analyzer->FakeRate(mu_p_T, mu_eta);
                FRweight = FR / (1 - FR);
//                if (Mgr.isMC == kTRUE)
//                {
//                    effweight = analyzer->EfficiencySF_EventWeight_electronFR(ele1_SF, ele2_SF, 1);
//                }
//                h_mass_temp->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
//                if (e_charge != mu_charge)
//                    h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
//                else
//                    h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                if (Mgr.isMC == kFALSE)
                {
                    n_mufail++;
                    nw_mufail += FRweight;
                }
                h_mass_test_mf->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                if (mu_p_T > 50)
//                {
//                    h_pT_mu->Fill(mu_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_pT_e->Fill(e_p_T, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                }
            }
            if (DEBUG == kTRUE) cout << "FR = " << FR << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Histogram filling -- //
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_mass_test->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

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
        h_mass_test->Write();
        h_mass_test_ef->Write();
        h_mass_test_mf->Write();
        h_pT_e->Write();
        h_pT_mu->Write();

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

    cout << "Events with failed electron: " << n_efail << endl;
    cout << "Events with failed muon: " << n_mufail << endl;
    cout << "FR weighted events with failed electron: " << nw_efail << endl;
    cout << "FR weighted events with failed muon: " << nw_mufail << endl;


    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_WJET_HistMaker()


/// --------------------------- ELECTRON PROMPT RATE ---------------------------- ///
void E_PR_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir1 = "/media/sf_DATA/FR/Electron/";
    TString Dir2 = "/media/sf_DATA/SelectedEE/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir1+"PR_Hist_E"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    TH1D *h_mass_data = new TH1D ("h_mass_data", "h_mass_data", 20, 81, 101);
    TH1D *h_mass_DY = new TH1D ("h_mass_DY", "h_mass_DY", 20, 81, 101);
    TH1D *h_mass_bkgr = new TH1D ("h_mass_bkgr", "h_mass_bkgr", 20, 81, 101);
    TH1D *h_mass_bkgf = new TH1D ("h_mass_bkgf", "h_mass_bkgf", 20, 81, 101);
    TH1D *h_pT_barrel_pass_data = new TH1D ("h_pT_barrel_pass_data", "h_pT_barrel_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_DY = new TH1D ("h_pT_barrel_pass_DY", "h_pT_barrel_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_bkgr = new TH1D ("h_pT_barrel_pass_bkgr", "h_pT_barrel_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_bkgf = new TH1D ("h_pT_barrel_pass_bkgf", "h_pT_barrel_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_data = new TH1D ("h_pT_endcap_pass_data", "h_pT_endcap_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_DY = new TH1D ("h_pT_endcap_pass_DY", "h_pT_endcap_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_bkgr = new TH1D ("h_pT_endcap_pass_bkgr", "h_pT_endcap_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_bkgf = new TH1D ("h_pT_endcap_pass_bkgf", "h_pT_endcap_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_data = new TH1D ("h_pT_endcap2_pass_data", "h_pT_endcap2_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_DY = new TH1D ("h_pT_endcap2_pass_DY", "h_pT_endcap2_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_bkgr = new TH1D ("h_pT_endcap2_pass_bkgr", "h_pT_endcap2_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_bkgf = new TH1D ("h_pT_endcap2_pass_bkgf", "h_pT_endcap2_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D* h_eta_pass_data = new TH1D("h_eta_pass_data", "h_eta_pass_data", 66, etabins);
    TH1D* h_eta_pass_DY = new TH1D("h_eta_pass_DY", "h_eta_pass_DY", 66, etabins);
    TH1D* h_eta_pass_bkgr = new TH1D("h_eta_pass_bkgr", "h_eta_pass_bkgr", 66, etabins);
    TH1D* h_eta_pass_bkgf = new TH1D("h_eta_pass_bkgf", "h_eta_pass_bkgf", 66, etabins);
    TH1D *h_pT_barrel_fail_data = new TH1D ("h_pT_barrel_fail_data", "h_pT_barrel_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_DY = new TH1D ("h_pT_barrel_fail_DY", "h_pT_barrel_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_bkgr = new TH1D ("h_pT_barrel_fail_bkgr", "h_pT_barrel_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_bkgf = new TH1D ("h_pT_barrel_fail_bkgf", "h_pT_barrel_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_data = new TH1D ("h_pT_endcap_fail_data", "h_pT_endcap_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_DY = new TH1D ("h_pT_endcap_fail_DY", "h_pT_endcap_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_bkgr = new TH1D ("h_pT_endcap_fail_bkgr", "h_pT_endcap_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_bkgf = new TH1D ("h_pT_endcap_fail_bkgf", "h_pT_endcap_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_data = new TH1D ("h_pT_endcap2_fail_data", "h_pT_endcap2_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_DY = new TH1D ("h_pT_endcap2_fail_DY", "h_pT_endcap2_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_bkgr = new TH1D ("h_pT_endcap2_fail_bkgr", "h_pT_endcap2_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_bkgf = new TH1D ("h_pT_endcap2_fail_bkgf", "h_pT_endcap2_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D* h_eta_fail_data = new TH1D("h_eta_fail_data", "h_eta_fail_data", 66, etabins);
    TH1D* h_eta_fail_DY = new TH1D("h_eta_fail_DY", "h_eta_fail_DY", 66, etabins);
    TH1D* h_eta_fail_bkgr = new TH1D("h_eta_fail_bkgr", "h_eta_fail_bkgr", 66, etabins);
    TH1D* h_eta_fail_bkgf = new TH1D("h_eta_fail_bkgf", "h_eta_fail_bkgf", 66, etabins);
    TH1D *h_MET_fail_data = new TH1D ("h_MET_fail_data", "h_MET_fail_data", 50, 0, 100);
    TH1D *h_MET_fail_DY = new TH1D ("h_MET_fail_DY", "h_MET_fail_DY", 50, 0, 100);
    TH1D *h_MET_fail_bkgr = new TH1D ("h_MET_fail_bkgr", "h_MET_fail_bkgr", 50, 0, 100);
    TH1D *h_MET_fail_bkgf = new TH1D ("h_MET_fail_bkgf", "h_MET_fail_bkgf", 50, 0, 100);
    TH1D *h_MT_fail_data = new TH1D ("h_MT_fail_data", "h_MT_fail_data", 100, 0, 200);
    TH1D *h_MT_fail_DY = new TH1D ("h_MT_fail_DY", "h_MT_fail_DY", 100, 0, 200);
    TH1D *h_MT_fail_bkgr = new TH1D ("h_MT_fail_bkgr", "h_MT_fail_bkgr", 100, 0, 200);
    TH1D *h_MT_fail_bkgf = new TH1D ("h_MT_fail_bkgf", "h_MT_fail_bkgf", 100, 0, 200);

    // -- For PU re-weighting -- //
    analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

    // -- For PVz reweighting -- //
    analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

    // FR part (at least one electron fails the selection)
    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir1 << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;

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
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passConvVeto = new std::vector<int>;
        std::vector<double> *relPFiso_Rho = new std::vector<double>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
        std::vector<double> *relTrkIso = new std::vector<double>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir1+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir1+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

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
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passConvVeto", 1);
        chain->SetBranchStatus("relPFiso_Rho", 1);
        chain->SetBranchStatus("relECALiso", 1);
        chain->SetBranchStatus("relHCALiso", 1);
        chain->SetBranchStatus("relTrkIso", 1);
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
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("relECALiso", &relECALiso);
        chain->SetBranchAddress("relHCALiso", &relHCALiso);
        chain->SetBranchAddress("relTrkIso", &relTrkIso);
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

            // Preselection
            if (p_T->size() != 2) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (fabs(etaSC->at(0)) >= 2.4 || fabs(etaSC->at(1)) >= 2.4) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;

            if (fabs(eta->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                              fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07)) continue;
            else if (fabs(eta->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13)) continue;
            if (fabs(eta->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                              fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07)) continue;
            else if (fabs(eta->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13)) continue;

//            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.00998 || fabs(dEtaInSeed->at(0)) >= 0.00311 ||
//                fabs(dPhiIn->at(0)) >= 0.07/*0.103*/ || HoverE->at(0) >= 0.13/*0.253*/ || InvEminusInvP->at(0) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.0298 || fabs(dEtaInSeed->at(0)) >= 0.00609 ||
//                fabs(dPhiIn->at(0)) >= 0.045 || HoverE->at(0) >= 0.0878 || InvEminusInvP->at(0) >= 0.13))
//                continue;
//            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.00998 || fabs(dEtaInSeed->at(1)) >= 0.00311 ||
//                fabs(dPhiIn->at(1)) >= 0.07/*0.103*/ || HoverE->at(1) >= 0.13/*0.253*/ || InvEminusInvP->at(1) >= 0.134))
//                continue;
//            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.0298 || fabs(dEtaInSeed->at(1)) >= 0.00609 ||
//                fabs(dPhiIn->at(1)) >= 0.045 || HoverE->at(1) >= 0.0878 || InvEminusInvP->at(1) >= 0.13))
//                continue;

            if (mHits->at(0) > 1 /*|| !passConvVeto->at(0)*/ || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 /*|| !passConvVeto->at(1)*/ || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;

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

            // Look for tight electron:
            Int_t tight0=0, tight1=0, i_fill=-1;
            if (fabs(etaSC->at(0)) < 1.4442 && Full5x5_SigmaIEtaIEta->at(0) < 0.00998 && fabs(dEtaInSeed->at(0)) < 0.00308 &&
                fabs(dPhiIn->at(0)) < 0.07/*0.0816*/ && HoverE->at(0) < 0.0414 && relPFiso_Rho->at(0) < 0.0588 &&
                InvEminusInvP->at(0) < 0.0129 && mHits->at(0) <= 1 && passConvVeto->at(0))
                tight0 = 1;
            if (fabs(etaSC->at(0)) > 1.566 && Full5x5_SigmaIEtaIEta->at(0) < 0.0292 && fabs(dEtaInSeed->at(0)) < 0.00605 &&
                fabs(dPhiIn->at(0)) < 0.0394 && HoverE->at(0) < 0.0641 && relPFiso_Rho->at(0) < 0.0571 &&
                InvEminusInvP->at(0) < 0.0129 && mHits->at(0) <= 1 && passConvVeto->at(0))
                tight0 = 1;
            if (fabs(etaSC->at(1)) < 1.4442 && Full5x5_SigmaIEtaIEta->at(1) < 0.00998 && fabs(dEtaInSeed->at(1)) < 0.00308 &&
                fabs(dPhiIn->at(1)) < 0.07/*0.0816*/ && HoverE->at(1) < 0.0414 && relPFiso_Rho->at(1) < 0.0588 &&
                InvEminusInvP->at(1) < 0.0129 && mHits->at(1) <= 1 && passConvVeto->at(1))
                tight1 = 1;
            if (fabs(etaSC->at(1)) > 1.566 && Full5x5_SigmaIEtaIEta->at(1) < 0.0292 && fabs(dEtaInSeed->at(1)) < 0.00605 &&
                fabs(dPhiIn->at(1)) < 0.0394 && HoverE->at(1) < 0.0641 && relPFiso_Rho->at(1) < 0.0571 &&
                InvEminusInvP->at(1) < 0.0129 && mHits->at(1) <= 1 && passConvVeto->at(1))
                tight1 = 1;

            if (!tight0 && !tight1) continue;
            else if (tight0 && tight1)
            {
                if (p_T->at(0) > p_T->at(1)) i_fill = 1;
                else i_fill = 0;
            }
            else if (tight0) i_fill = 1;
            else if (tight1) i_fill = 0;

            nPass++;            

            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();
            if (mass < 76 || mass >106) continue;
            Double_t dTheta = phi->at(1-i_fill) - MET_phi;
            Double_t MT = sqrt(2 * p_T->at(1-i_fill) * MET_pT * (1 - cos(dTheta)));
//            if (MT >= 60) continue;

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << "\tTop pT weight: " << TopPtWeight << endl;

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            if (Mgr.Type == "SIGNAL")
            {
                h_mass_DY->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (passMediumID->at(i_fill))
                {
                    h_eta_pass_DY->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(etaSC->at(i_fill)) < 1.4442) // Barrel
                        h_pT_barrel_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1) // Endcap
                        h_pT_endcap_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4) // Far endcap
                        h_pT_endcap2_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
                else
                {
                    h_eta_fail_DY->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MET_fail_DY->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MT_fail_DY->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(etaSC->at(i_fill)) < 1.4442) // Barrel
                        h_pT_barrel_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1) // Endcap
                        h_pT_endcap_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4) // Far endcap
                        h_pT_endcap2_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }

            }
            else if (Mgr.Type == "BKG")
            {
                if (Mgr.Tag[0].Contains("QCD") || Mgr.Tag[0].Contains("Jet"))
                {
                    h_mass_bkgf->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                    if (passMediumID->at(i_fill))
                    {
                        h_eta_pass_bkgf->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(etaSC->at(i_fill)) < 1.4442)
                            h_pT_barrel_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                            h_pT_endcap_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                            h_pT_endcap2_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else
                    {
                        h_eta_fail_bkgf->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_bkgf->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_bkgf->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                        if (fabs(etaSC->at(i_fill)) < 1.4442)
                            h_pT_barrel_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                            h_pT_endcap_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                            h_pT_endcap2_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }

                }
                else
                {
                    h_mass_bkgr->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (passMediumID->at(i_fill))
                    {
                        h_eta_pass_bkgr->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(etaSC->at(i_fill)) < 1.4442)
                            h_pT_barrel_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                            h_pT_endcap_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                            h_pT_endcap2_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else
                    {
                        h_eta_fail_bkgr->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_bkgr->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_bkgr->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                        if (fabs(etaSC->at(i_fill)) < 1.4442)
                            h_pT_barrel_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                            h_pT_endcap_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                            h_pT_endcap2_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
            }
            else if (Mgr.Type == "DATA")
            {
                h_mass_data->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (passMediumID->at(i_fill))
                {
                    h_eta_pass_data->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(etaSC->at(i_fill)) < 1.4442)
                        h_pT_barrel_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                        h_pT_endcap_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                        h_pT_endcap2_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
                else
                {
                    h_eta_fail_data->Fill(etaSC->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MET_fail_data->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MT_fail_data->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                    if (fabs(etaSC->at(i_fill)) < 1.4442 && !passMediumID->at(i_fill))
                        h_pT_barrel_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) > 1.566 && fabs(etaSC->at(i_fill)) < 2.1)
                        h_pT_endcap_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(etaSC->at(i_fill)) >= 2.1 && fabs(etaSC->at(i_fill)) < 2.4)
                        h_pT_endcap2_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
            }

        }// End of event iteration

        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    f->cd();
    cout << "\tWriting into file...";

    h_mass_DY->Write();
    h_mass_bkgr->Write();
    h_mass_bkgf->Write();
    h_mass_data->Write();
    h_pT_barrel_pass_DY->Write();
    h_pT_barrel_pass_bkgr->Write();
    h_pT_barrel_pass_bkgf->Write();
    h_pT_barrel_pass_data->Write();
    h_pT_endcap_pass_DY->Write();
    h_pT_endcap_pass_bkgr->Write();
    h_pT_endcap_pass_bkgf->Write();
    h_pT_endcap_pass_data->Write();
    h_pT_endcap2_pass_DY->Write();
    h_pT_endcap2_pass_bkgr->Write();
    h_pT_endcap2_pass_bkgf->Write();
    h_pT_endcap2_pass_data->Write();
    h_eta_pass_DY->Write();
    h_eta_pass_bkgr->Write();
    h_eta_pass_bkgf->Write();
    h_eta_pass_data->Write();
    h_pT_barrel_fail_DY->Write();
    h_pT_barrel_fail_bkgr->Write();
    h_pT_barrel_fail_bkgf->Write();
    h_pT_barrel_fail_data->Write();
    h_pT_endcap_fail_DY->Write();
    h_pT_endcap_fail_bkgr->Write();
    h_pT_endcap_fail_bkgf->Write();
    h_pT_endcap_fail_data->Write();
    h_pT_endcap2_fail_DY->Write();
    h_pT_endcap2_fail_bkgr->Write();
    h_pT_endcap2_fail_bkgf->Write();
    h_pT_endcap2_fail_data->Write();
    h_eta_fail_DY->Write();
    h_eta_fail_bkgr->Write();
    h_eta_fail_bkgf->Write();
    h_eta_fail_data->Write();
    h_MET_fail_DY->Write();
    h_MET_fail_bkgr->Write();
    h_MET_fail_bkgf->Write();
    h_MET_fail_data->Write();
    h_MT_fail_DY->Write();
    h_MT_fail_bkgr->Write();
    h_MT_fail_bkgf->Write();
    h_MT_fail_data->Write();

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_PR_HistDrawer()


/// --------------------------- MUON PROMPT RATE ---------------------------- ///
void Mu_PR_HistMaker (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir1 = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir1+"PR_Hist_Mu"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    TH1D *h_mass_data = new TH1D ("h_mass_data", "h_mass_data", 20, 81, 101);
    TH1D *h_mass_DY = new TH1D ("h_mass_DY", "h_mass_DY", 20, 81, 101);
    TH1D *h_mass_bkgr = new TH1D ("h_mass_bkgr", "h_mass_bkgr", 20, 81, 101);
    TH1D *h_mass_bkgf = new TH1D ("h_mass_bkgf", "h_mass_bkgf", 20, 81, 101);
    TH1D *h_pT_barrel_pass_data = new TH1D ("h_pT_barrel_pass_data", "h_pT_barrel_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_DY = new TH1D ("h_pT_barrel_pass_DY", "h_pT_barrel_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_bkgr = new TH1D ("h_pT_barrel_pass_bkgr", "h_pT_barrel_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_pass_bkgf = new TH1D ("h_pT_barrel_pass_bkgf", "h_pT_barrel_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_pass_data = new TH1D ("h_pT_barrel2_pass_data", "h_pT_barrel2_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_pass_DY = new TH1D ("h_pT_barrel2_pass_DY", "h_pT_barrel2_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_pass_bkgr = new TH1D ("h_pT_barrel2_pass_bkgr", "h_pT_barrel2_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_pass_bkgf = new TH1D ("h_pT_barrel2_pass_bkgf", "h_pT_barrel2_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_data = new TH1D ("h_pT_endcap_pass_data", "h_pT_endcap_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_DY = new TH1D ("h_pT_endcap_pass_DY", "h_pT_endcap_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_bkgr = new TH1D ("h_pT_endcap_pass_bkgr", "h_pT_endcap_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_pass_bkgf = new TH1D ("h_pT_endcap_pass_bkgf", "h_pT_endcap_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_data = new TH1D ("h_pT_endcap2_pass_data", "h_pT_endcap2_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_DY = new TH1D ("h_pT_endcap2_pass_DY", "h_pT_endcap2_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_bkgr = new TH1D ("h_pT_endcap2_pass_bkgr", "h_pT_endcap2_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_pass_bkgf = new TH1D ("h_pT_endcap2_pass_bkgf", "h_pT_endcap2_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_eta_pass_data = new TH1D("h_eta_pass_data", "h_eta_pass_data", 96, -2.4, 2.4);
    TH1D *h_eta_pass_DY = new TH1D("h_eta_pass_DY", "h_eta_pass_DY", 96, -2.4, 2.4);
    TH1D *h_eta_pass_bkgr = new TH1D("h_eta_pass_bkgr", "h_eta_pass_bkgr", 96, -2.4, 2.4);
    TH1D *h_eta_pass_bkgf = new TH1D("h_eta_pass_bkgf", "h_eta_pass_bkgf", 96, -2.4, 2.4);
    TH1D *h_pT_barrel_fail_data = new TH1D ("h_pT_barrel_fail_data", "h_pT_barrel_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_DY = new TH1D ("h_pT_barrel_fail_DY", "h_pT_barrel_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_bkgr = new TH1D ("h_pT_barrel_fail_bkgr", "h_pT_barrel_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_fail_bkgf = new TH1D ("h_pT_barrel_fail_bkgf", "h_pT_barrel_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_fail_data = new TH1D ("h_pT_barrel2_fail_data", "h_pT_barrel2_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_fail_DY = new TH1D ("h_pT_barrel2_fail_DY", "h_pT_barrel2_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_fail_bkgr = new TH1D ("h_pT_barrel2_fail_bkgr", "h_pT_barrel2_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel2_fail_bkgf = new TH1D ("h_pT_barrel2_fail_bkgf", "h_pT_barrel2_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_data = new TH1D ("h_pT_endcap_fail_data", "h_pT_endcap_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_DY = new TH1D ("h_pT_endcap_fail_DY", "h_pT_endcap_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_bkgr = new TH1D ("h_pT_endcap_fail_bkgr", "h_pT_endcap_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_fail_bkgf = new TH1D ("h_pT_endcap_fail_bkgf", "h_pT_endcap_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_data = new TH1D ("h_pT_endcap2_fail_data", "h_pT_endcap2_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_DY = new TH1D ("h_pT_endcap2_fail_DY", "h_pT_endcap2_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_bkgr = new TH1D ("h_pT_endcap2_fail_bkgr", "h_pT_endcap2_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap2_fail_bkgf = new TH1D ("h_pT_endcap2_fail_bkgf", "h_pT_endcap2_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_eta_fail_data = new TH1D("h_eta_fail_data", "h_eta_fail_data", 96, -2.4, 2.4);
    TH1D *h_eta_fail_DY = new TH1D("h_eta_fail_DY", "h_eta_fail_DY", 96, -2.4, 2.4);
    TH1D *h_eta_fail_bkgr = new TH1D("h_eta_fail_bkgr", "h_eta_fail_bkgr", 96, -2.4, 2.4);
    TH1D *h_eta_fail_bkgf = new TH1D("h_eta_fail_bkgf", "h_eta_fail_bkgf", 96, -2.4, 2.4);
    TH1D *h_MET_fail_data = new TH1D ("h_MET_fail_data", "h_MET_fail_data", 50, 0, 100);
    TH1D *h_MET_fail_DY = new TH1D ("h_MET_fail_DY", "h_MET_fail_DY", 50, 0, 100);
    TH1D *h_MET_fail_bkgr = new TH1D ("h_MET_fail_bkgr", "h_MET_fail_bkgr", 50, 0, 100);
    TH1D *h_MET_fail_bkgf = new TH1D ("h_MET_fail_bkgf", "h_MET_fail_bkgf", 50, 0, 100);
    TH1D *h_MT_fail_data = new TH1D ("h_MT_fail_data", "h_MT_fail_data", 100, 0, 200);
    TH1D *h_MT_fail_DY = new TH1D ("h_MT_fail_DY", "h_MT_fail_DY", 100, 0, 200);
    TH1D *h_MT_fail_bkgr = new TH1D ("h_MT_fail_bkgr", "h_MT_fail_bkgr", 100, 0, 200);
    TH1D *h_MT_fail_bkgf = new TH1D ("h_MT_fail_bkgf", "h_MT_fail_bkgf", 100, 0, 200);

    // -- For PU re-weighting -- //
    analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

    // -- For PVz reweighting -- //
    analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir1 << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        std::vector<string> *trig_name = new std::vector<string>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t MET_pT, MET_phi;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir1+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir1+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("trig_name", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("prescale_factor", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
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
        chain->SetBranchAddress("trig_name", &trig_name);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
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

            if (p_T->size() != 2) continue;
            if (charge->at(0) == charge->at(1)) continue;
            if (p_T->at(0) <= 17 || p_T->at(1) <= 17) continue;
            if (p_T->at(0) <= 52 && p_T->at(1) <= 52) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            // Look for tight electron:
            Int_t tight0=0, tight1=0, i_fill=-1;
            if (relPFiso->at(0) < 0.10)
                tight0 = 1;
            if (relPFiso->at(1) < 0.10)
                tight1 = 1;

            if (!tight0 && !tight1) continue;
            else if (tight0 && tight1)
            {
                if (p_T->at(0) > p_T->at(1)) i_fill = 1;
                else i_fill = 0;
            }
            else if (tight0) i_fill = 1;
            else if (tight1) i_fill = 0;

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (mu1+mu2).M();
            if (mass < 76 || mass >106) continue;
            Double_t dTheta = phi->at(1-i_fill) - MET_phi;
            Double_t MT = sqrt(2 * p_T->at(1-i_fill) * MET_pT * (1 - cos(dTheta)));
//            if (MT >= 60) continue;
            nPass++;

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << "\tTop pT weight: " << TopPtWeight << endl;

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            if (Mgr.Type == "SIGNAL")
            {
                h_mass_DY->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (relPFiso->at(i_fill) < 0.15)
                {
                    h_eta_pass_DY->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                        h_pT_barrel_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                        h_pT_barrel2_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8) // Endcap
                        h_pT_endcap_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4) // Far endcap
                        h_pT_endcap2_pass_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
                else
                {
                    h_eta_fail_DY->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MET_fail_DY->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MT_fail_DY->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                        h_pT_barrel_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                        h_pT_barrel2_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8) // Endcap
                        h_pT_endcap_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4) // Far endcap
                        h_pT_endcap2_fail_DY->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }

            }
            else if (Mgr.Type == "BKG")
            {
                if (Mgr.Tag[0].Contains("QCD") || Mgr.Tag[0].Contains("Jet"))
                {
                    h_mass_bkgf->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                    if (relPFiso->at(i_fill) < 0.15)
                    {
                        h_eta_pass_bkgf->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                            h_pT_barrel_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                            h_pT_barrel2_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                            h_pT_endcap_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                            h_pT_endcap2_pass_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else
                    {
                        h_eta_fail_bkgf->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_bkgf->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_bkgf->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                        if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                            h_pT_barrel_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                            h_pT_barrel2_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                            h_pT_endcap_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                            h_pT_endcap2_fail_bkgf->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }

                }
                else
                {
                    h_mass_bkgr->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (relPFiso->at(i_fill) < 0.15)
                    {
                        h_eta_pass_bkgr->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                            h_pT_barrel_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                            h_pT_barrel2_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                            h_pT_endcap_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                            h_pT_endcap2_pass_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else
                    {
                        h_eta_fail_bkgr->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_bkgr->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_bkgr->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                        if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                            h_pT_barrel_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                            h_pT_barrel2_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                            h_pT_endcap_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                            h_pT_endcap2_fail_bkgr->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
            }
            else if (Mgr.Type == "DATA")
            {
                h_mass_data->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (relPFiso->at(i_fill) < 0.15)
                {
                    h_eta_pass_data->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                        h_pT_barrel_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                        h_pT_barrel2_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                        h_pT_endcap_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                        h_pT_endcap2_pass_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
                else
                {
                    h_eta_fail_data->Fill(eta->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MET_fail_data->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    h_MT_fail_data->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

                    if (fabs(eta->at(i_fill)) < 0.7) // Barrel
                        h_pT_barrel_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 0.7 && fabs(eta->at(i_fill)) < 1.2) // Far barrel
                        h_pT_barrel2_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.2 && fabs(eta->at(i_fill)) < 1.8)
                        h_pT_endcap_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_fill)) >= 1.8 && fabs(eta->at(i_fill)) < 2.4)
                        h_pT_endcap2_fail_data->Fill(p_T->at(i_fill), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
            }

        }// End of event iteration

        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    f->cd();
    cout << "\tWriting into file...";

    h_mass_DY->Write();
    h_mass_bkgr->Write();
    h_mass_bkgf->Write();
    h_mass_data->Write();
    h_pT_barrel_pass_DY->Write();
    h_pT_barrel_pass_bkgr->Write();
    h_pT_barrel_pass_bkgf->Write();
    h_pT_barrel_pass_data->Write();
    h_pT_barrel2_pass_DY->Write();
    h_pT_barrel2_pass_bkgr->Write();
    h_pT_barrel2_pass_bkgf->Write();
    h_pT_barrel2_pass_data->Write();
    h_pT_endcap_pass_DY->Write();
    h_pT_endcap_pass_bkgr->Write();
    h_pT_endcap_pass_bkgf->Write();
    h_pT_endcap_pass_data->Write();
    h_pT_endcap2_pass_DY->Write();
    h_pT_endcap2_pass_bkgr->Write();
    h_pT_endcap2_pass_bkgf->Write();
    h_pT_endcap2_pass_data->Write();
    h_eta_pass_DY->Write();
    h_eta_pass_bkgr->Write();
    h_eta_pass_bkgf->Write();
    h_eta_pass_data->Write();
    h_pT_barrel_fail_DY->Write();
    h_pT_barrel_fail_bkgr->Write();
    h_pT_barrel_fail_bkgf->Write();
    h_pT_barrel_fail_data->Write();
    h_pT_barrel2_fail_DY->Write();
    h_pT_barrel2_fail_bkgr->Write();
    h_pT_barrel2_fail_bkgf->Write();
    h_pT_barrel2_fail_data->Write();
    h_pT_endcap_fail_DY->Write();
    h_pT_endcap_fail_bkgr->Write();
    h_pT_endcap_fail_bkgf->Write();
    h_pT_endcap_fail_data->Write();
    h_pT_endcap2_fail_DY->Write();
    h_pT_endcap2_fail_bkgr->Write();
    h_pT_endcap2_fail_bkgf->Write();
    h_pT_endcap2_fail_data->Write();
    h_eta_fail_DY->Write();
    h_eta_fail_bkgr->Write();
    h_eta_fail_bkgf->Write();
    h_eta_fail_data->Write();
    h_MET_fail_DY->Write();
    h_MET_fail_bkgr->Write();
    h_MET_fail_bkgf->Write();
    h_MET_fail_data->Write();
    h_MT_fail_DY->Write();
    h_MT_fail_bkgr->Write();
    h_MT_fail_bkgf->Write();
    h_MT_fail_data->Write();

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_PR_HistDrawer()


/// --------------------------- DY MediumID efficiency check ---------------------------- ///
void E_FR_HistMaker_alt (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir1 = "/media/sf_DATA/FR/Electron/";
    TString Dir2 = "/media/sf_DATA/SelectedEE/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir1+"For_FakeRate_electron_alt"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    TH1D *h_pT_barrel_lead_pass_data = new TH1D ("h_pT_barrel_lead_pass_data", "h_pT_barrel_lead_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_pass_DY = new TH1D ("h_pT_barrel_lead_pass_DY", "h_pT_barrel_lead_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_pass_bkgr = new TH1D ("h_pT_barrel_lead_pass_bkgr", "h_pT_barrel_lead_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_pass_bkgf = new TH1D ("h_pT_barrel_lead_pass_bkgf", "h_pT_barrel_lead_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_pass_data = new TH1D ("h_pT_endcap_lead_pass_data", "h_pT_endcap_lead_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_pass_DY = new TH1D ("h_pT_endcap_lead_pass_DY", "h_pT_endcap_lead_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_pass_bkgr = new TH1D ("h_pT_endcap_lead_pass_bkgr", "h_pT_endcap_lead_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_pass_bkgf = new TH1D ("h_pT_endcap_lead_pass_bkgf", "h_pT_endcap_lead_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_fail_data = new TH1D ("h_pT_barrel_lead_fail_data", "h_pT_barrel_lead_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_fail_DY = new TH1D ("h_pT_barrel_lead_fail_DY", "h_pT_barrel_lead_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_fail_bkgr = new TH1D ("h_pT_barrel_lead_fail_bkgr", "h_pT_barrel_lead_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_lead_fail_bkgf = new TH1D ("h_pT_barrel_lead_fail_bkgf", "h_pT_barrel_lead_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_fail_data = new TH1D ("h_pT_endcap_lead_fail_data", "h_pT_endcap_lead_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_fail_DY = new TH1D ("h_pT_endcap_lead_fail_DY", "h_pT_endcap_lead_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_fail_bkgr = new TH1D ("h_pT_endcap_lead_fail_bkgr", "h_pT_endcap_lead_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_lead_fail_bkgf = new TH1D ("h_pT_endcap_lead_fail_bkgf", "h_pT_endcap_lead_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);

    TH1D *h_pT_barrel_sub_pass_data = new TH1D ("h_pT_barrel_sub_pass_data", "h_pT_barrel_sub_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_pass_DY = new TH1D ("h_pT_barrel_sub_pass_DY", "h_pT_barrel_sub_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_pass_bkgr = new TH1D ("h_pT_barrel_sub_pass_bkgr", "h_pT_barrel_sub_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_pass_bkgf = new TH1D ("h_pT_barrel_sub_pass_bkgf", "h_pT_barrel_sub_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_pass_data = new TH1D ("h_pT_endcap_sub_pass_data", "h_pT_endcap_sub_pass_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_pass_DY = new TH1D ("h_pT_endcap_sub_pass_DY", "h_pT_endcap_sub_pass_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_pass_bkgr = new TH1D ("h_pT_endcap_sub_pass_bkgr", "h_pT_endcap_sub_pass_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_pass_bkgf = new TH1D ("h_pT_endcap_sub_pass_bkgf", "h_pT_endcap_sub_pass_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_fail_data = new TH1D ("h_pT_barrel_sub_fail_data", "h_pT_barrel_sub_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_fail_DY = new TH1D ("h_pT_barrel_sub_fail_DY", "h_pT_barrel_sub_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_fail_bkgr = new TH1D ("h_pT_barrel_sub_fail_bkgr", "h_pT_barrel_sub_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_barrel_sub_fail_bkgf = new TH1D ("h_pT_barrel_sub_fail_bkgf", "h_pT_barrel_sub_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_fail_data = new TH1D ("h_pT_endcap_sub_fail_data", "h_pT_endcap_sub_fail_data", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_fail_DY = new TH1D ("h_pT_endcap_sub_fail_DY", "h_pT_endcap_sub_fail_DY", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_fail_bkgr = new TH1D ("h_pT_endcap_sub_fail_bkgr", "h_pT_endcap_sub_fail_bkgr", nPtBin_DY, analyzer->ptbin_DY);
    TH1D *h_pT_endcap_sub_fail_bkgf = new TH1D ("h_pT_endcap_sub_fail_bkgf", "h_pT_endcap_sub_fail_bkgf", nPtBin_DY, analyzer->ptbin_DY);

    TH1D *h_MET_fail_data = new TH1D ("h_MET_fail_data", "h_MET_fail_data", 50, 0, 100);
    TH1D *h_MET_fail_DY = new TH1D ("h_MET_fail_DY", "h_MET_fail_DY", 50, 0, 100);
    TH1D *h_MET_fail_bkgr = new TH1D ("h_MET_fail_bkgr", "h_MET_fail_bkgr", 50, 0, 100);
    TH1D *h_MET_fail_bkgf = new TH1D ("h_MET_fail_bkgf", "h_MET_fail_bkgf", 50, 0, 100);
    TH1D *h_MT_fail_data = new TH1D ("h_MT_fail_data", "h_MT_fail_data", 100, 0, 200);
    TH1D *h_MT_fail_DY = new TH1D ("h_MT_fail_DY", "h_MT_fail_DY", 100, 0, 200);
    TH1D *h_MT_fail_bkgr = new TH1D ("h_MT_fail_bkgr", "h_MT_fail_bkgr", 100, 0, 200);
    TH1D *h_MT_fail_bkgf = new TH1D ("h_MT_fail_bkgf", "h_MT_fail_bkgf", 100, 0, 200);

    // -- For PU re-weighting -- //
    analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

    // -- For PVz reweighting -- //
    analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

    // FR part (at least one electron fails the selection)
    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir1 << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;

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
        chain->Add(Dir1+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir1+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

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

            // QCD selection
            if (p_T->size() != 2) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;
//            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
//                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
//            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
//                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
//            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
//            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (!passMediumID->at(0) && !passMediumID->at(1)) continue;
            if (charge->at(0) != charge->at(1)) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;


            Int_t i_lead = 0;
            if (p_T->at(1) > p_T->at(0)) i_lead = 1;
//            if (passMediumID->at(1)) i_lead = 1;
//            if (MET_pT > 20) continue;

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

            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();
            if (mass > 76 && mass <106) continue;
            Double_t dTheta = phi->at(i_lead) - MET_phi;
            Double_t MT = sqrt(2 * p_T->at(i_lead) * MET_pT * (1 - cos(dTheta)));
            if (MT < 30) continue;

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << "\tTop pT weight: " << TopPtWeight << endl;

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            if (Mgr.Type == "SIGNAL")
            {
                if (passMediumID->at(i_lead))
                {
                    if (fabs(eta->at(1-i_lead)) < 1.4442 && passMediumID->at(1-i_lead))
                        h_pT_barrel_sub_pass_DY->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(1-i_lead)) > 1.566 && passMediumID->at(1-i_lead))
                        h_pT_endcap_sub_pass_DY->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(1-i_lead)) < 1.4442 && !passMediumID->at(1-i_lead))
                    {
                        h_pT_barrel_sub_fail_DY->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_DY->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_DY->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (fabs(eta->at(1-i_lead)) > 1.566 && !passMediumID->at(1-i_lead))
                    {
                        h_pT_endcap_sub_fail_DY->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_DY->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_DY->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                if (passMediumID->at(1-i_lead))
                {
                    if (fabs(eta->at(i_lead)) < 1.4442 && passMediumID->at(i_lead))
                        h_pT_barrel_lead_pass_DY->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) > 1.566 && passMediumID->at(i_lead))
                        h_pT_endcap_lead_pass_DY->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) < 1.4442 && !passMediumID->at(i_lead))
                        h_pT_barrel_lead_fail_DY->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) > 1.566 && !passMediumID->at(i_lead))
                        h_pT_endcap_lead_fail_DY->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
            }
            else if (Mgr.Type == "BKG")
            {
                if (Mgr.Tag[0].Contains("QCD") || Mgr.Tag[0].Contains("Jet"))
                {
                    if (passMediumID->at(i_lead))
                    {
                        if (fabs(eta->at(1-i_lead)) < 1.4442 && passMediumID->at(1-i_lead))
                            h_pT_barrel_sub_pass_bkgf->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(1-i_lead)) > 1.566 && passMediumID->at(1-i_lead))
                            h_pT_endcap_sub_pass_bkgf->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(1-i_lead)) < 1.4442 && !passMediumID->at(1-i_lead))
                        {
                            h_pT_barrel_sub_fail_bkgf->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MET_fail_bkgf->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MT_fail_bkgf->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (fabs(eta->at(1-i_lead)) > 1.566 && !passMediumID->at(1-i_lead))
                        {
                            h_pT_endcap_sub_fail_bkgf->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MET_fail_bkgf->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MT_fail_bkgf->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    if (passMediumID->at(1-i_lead))
                    {
                        if (fabs(eta->at(i_lead)) < 1.4442 && passMediumID->at(i_lead))
                            h_pT_barrel_lead_pass_bkgf->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) > 1.566 && passMediumID->at(i_lead))
                            h_pT_endcap_lead_pass_bkgf->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) < 1.4442 && !passMediumID->at(i_lead))
                            h_pT_barrel_lead_fail_bkgf->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) > 1.566 && !passMediumID->at(i_lead))
                            h_pT_endcap_lead_fail_bkgf->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                else
                {
                    if (passMediumID->at(i_lead))
                    {
                        if (fabs(eta->at(1-i_lead)) < 1.4442 && passMediumID->at(1-i_lead))
                            h_pT_barrel_sub_pass_bkgr->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(1-i_lead)) > 1.566 && passMediumID->at(1-i_lead))
                            h_pT_endcap_sub_pass_bkgr->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(1-i_lead)) < 1.4442 && !passMediumID->at(1-i_lead))
                        {
                            h_pT_barrel_sub_fail_bkgr->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MET_fail_bkgr->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MT_fail_bkgr->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (fabs(eta->at(1-i_lead)) > 1.566 && !passMediumID->at(1-i_lead))
                        {
                            h_pT_endcap_sub_fail_bkgr->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MET_fail_bkgr->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                            h_MT_fail_bkgr->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    if (passMediumID->at(1-i_lead))
                    {
                        if (fabs(eta->at(i_lead)) < 1.4442 && passMediumID->at(i_lead))
                            h_pT_barrel_lead_pass_bkgr->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) > 1.566 && passMediumID->at(i_lead))
                            h_pT_endcap_lead_pass_bkgr->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) < 1.4442 && !passMediumID->at(i_lead))
                            h_pT_barrel_lead_fail_bkgr->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(eta->at(i_lead)) > 1.566 && !passMediumID->at(i_lead))
                            h_pT_endcap_lead_fail_bkgr->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
            }
            else if (Mgr.Type == "DATA")
            {
                if (passMediumID->at(i_lead))
                {
                    if (fabs(eta->at(1-i_lead)) < 1.4442 && passMediumID->at(1-i_lead))
                        h_pT_barrel_sub_pass_data->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(1-i_lead)) > 1.566 && passMediumID->at(1-i_lead))
                        h_pT_endcap_sub_pass_data->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(1-i_lead)) < 1.4442 && !passMediumID->at(1-i_lead))
                    {
                        h_pT_barrel_sub_fail_data->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_data->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_data->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (fabs(eta->at(1-i_lead)) > 1.566 && !passMediumID->at(1-i_lead))
                    {
                        h_pT_endcap_sub_fail_data->Fill(p_T->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MET_fail_data->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        h_MT_fail_data->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                if (passMediumID->at(1-i_lead))
                {
                    if (fabs(eta->at(i_lead)) < 1.4442 && passMediumID->at(i_lead))
                        h_pT_barrel_lead_pass_data->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) > 1.566 && passMediumID->at(i_lead))
                        h_pT_endcap_lead_pass_data->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) < 1.4442 && !passMediumID->at(i_lead))
                        h_pT_barrel_lead_fail_data->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(eta->at(i_lead)) > 1.566 && !passMediumID->at(i_lead))
                        h_pT_endcap_lead_fail_data->Fill(p_T->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
            }

        }// End of event iteration

        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    // Signal part (both electrons pas the selection)
    LocalFileMgr LMgr;
    for (SelProc_t pr=_EE_DY_10to50; pr<_EndOf_EE_DoubleEG_Normal; pr=next(pr))
    {
        LMgr.SetProc(pr);

        cout << "Process: " << LMgr.Procname[pr] << endl;
        cout << "Type: " << LMgr.Type << endl;
        cout << "DATA location: " << LMgr.BaseLocation << endl;

        //Loop for all samples
        const Int_t Ntup = LMgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            cout << "\t<" << LMgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(LMgr.TreeName[i_tup]);
            chain->Add(LMgr.FullLocation[i_tup]);

            SelectedEE_t *EE = new SelectedEE_t();
            EE->CreateFromChain(chain);

            Int_t NEvents = chain->GetEntries();
            if (NEvents != LMgr.nEvents[i_tup]) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights:: " << LMgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            if (DEBUG == kTRUE) NEvents = 10;
            myProgressBar_t bar(NEvents);

            for (Int_t i=0; i<NEvents; i++)
            {
                EE->GetEvent(i);
                if (!DEBUG) bar.Draw(i);

                if (!EE->isSelPassed) continue;
                if (EE->Electron_charge->at(0) != EE->Electron_charge->at(1)) continue;
//                if (EE->MET_pT > 20) continue;
                Int_t i_lead = 0;
                if (EE->Electron_pT->at(1) > EE->Electron_pT->at(0)) i_lead = 1;

                TLorentzVector ele1, ele2;
                ele1.SetPtEtaPhiE(EE->Electron_pT->at(0), EE->Electron_eta->at(0), EE->Electron_phi->at(0), EE->Electron_Energy->at(0));
                ele2.SetPtEtaPhiE(EE->Electron_pT->at(1), EE->Electron_eta->at(1), EE->Electron_phi->at(1), EE->Electron_Energy->at(1));
                Double_t mass = (ele1 + ele2).M();
                if (mass > 76 && mass < 106) continue;

                Double_t dTheta = EE->Electron_phi->at(0) - EE->MET_phi;
                Double_t MT = sqrt(2 * EE->Electron_pT->at(0) * EE->MET_pT * (1 - cos(dTheta)));
                if (MT < 30) continue;

                // -- Bring the weights -- //
                Double_t GenWeight = EE->GENEvt_weight;
                if (GenWeight > 1) GenWeight = 1;
                else if (GenWeight < -1) GenWeight = -1;

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if(LMgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(EE->nPileUp);

                // -- PVz weights -- //
                Double_t PVzWeight = 1;
                if(LMgr.isMC == kTRUE /*&& pr < _EE_QCDEMEnriched_20to30*/) PVzWeight = analyzer->PVzWeightValue(EE->PVz);

                // -- L1 prefiring weights -- //
                Double_t L1weight = 1;
                if (LMgr.isMC == kTRUE /*&& pr < _EE_QCDEMEnriched_20to30*/) L1weight = EE->_prefiringweight;

                // -- Top Pt weights -- //
                Double_t TopPtWeight = 1;
                if (LMgr.isMC == kTRUE /*&& pr < _EE_QCDEMEnriched_20to30*/) TopPtWeight = EE->_topPtWeight;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if (LMgr.isMC == kTRUE) TotWeight = (Lumi * LMgr.Xsec[i_tup] / LMgr.Wsum[i_tup]) * GenWeight;
                if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;


                if (LMgr.Type == "SIGNAL")
                {                   
                    if (fabs(EE->Electron_eta->at(i_lead)) < 1.4442)
                        h_pT_barrel_lead_pass_DY->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(EE->Electron_eta->at(i_lead)) > 1.566)
                        h_pT_endcap_lead_pass_DY->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(EE->Electron_eta->at(1-i_lead)) < 1.4442)
                        h_pT_barrel_sub_pass_DY->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(EE->Electron_eta->at(1-i_lead)) > 1.566)
                        h_pT_endcap_sub_pass_DY->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (LMgr.Type == "BKG")
                {
                    if (LMgr.Tag[0].Contains("QCD") || LMgr.Tag[0].Contains("Jet"))
                    {                        
                        if (fabs(EE->Electron_eta->at(i_lead)) < 1.4442)
                            h_pT_barrel_lead_pass_bkgf->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(EE->Electron_eta->at(i_lead)) > 1.566)
                            h_pT_endcap_lead_pass_bkgf->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(EE->Electron_eta->at(1-i_lead)) < 1.4442)
                            h_pT_barrel_sub_pass_bkgf->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(EE->Electron_eta->at(1-i_lead)) > 1.566)
                            h_pT_endcap_sub_pass_bkgf->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else
                    {                        
                        if (fabs(EE->Electron_eta->at(i_lead)) < 1.4442)
                            h_pT_barrel_lead_pass_bkgr->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(EE->Electron_eta->at(i_lead)) > 1.566)
                            h_pT_endcap_lead_pass_bkgr->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        if (fabs(EE->Electron_eta->at(1-i_lead)) < 1.4442)
                            h_pT_barrel_sub_pass_bkgr->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                        else if (fabs(EE->Electron_eta->at(1-i_lead)) > 1.566)
                            h_pT_endcap_sub_pass_bkgr->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                else if (LMgr.Type == "DATA")
                {                    
                    if (fabs(EE->Electron_eta->at(i_lead)) < 1.4442)
                        h_pT_barrel_lead_pass_data->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(EE->Electron_eta->at(i_lead)) > 1.566)
                        h_pT_endcap_lead_pass_data->Fill(EE->Electron_pT->at(i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(EE->Electron_eta->at(1-i_lead)) < 1.4442)
                        h_pT_barrel_sub_pass_data->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                    else if (fabs(EE->Electron_eta->at(1-i_lead)) > 1.566)
                        h_pT_endcap_sub_pass_data->Fill(EE->Electron_pT->at(1-i_lead), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                }

            } // End of event iteration

            if(LMgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*LMgr.Xsec[i_tup]/LMgr.Wsum[i_tup]);

        }// End of i_tup iteration

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _EE_DY_10to50) break;

        cout << "===========================================================\n" << endl;

    }// End of pr iteration

    f->cd();
    cout << "\tWriting into file...";

    h_pT_barrel_lead_pass_DY->Write();
    h_pT_barrel_lead_pass_bkgr->Write();
    h_pT_barrel_lead_pass_bkgf->Write();
    h_pT_barrel_lead_pass_data->Write();
    h_pT_endcap_lead_pass_DY->Write();
    h_pT_endcap_lead_pass_bkgr->Write();
    h_pT_endcap_lead_pass_bkgf->Write();
    h_pT_endcap_lead_pass_data->Write();
    h_pT_barrel_lead_fail_DY->Write();
    h_pT_barrel_lead_fail_bkgr->Write();
    h_pT_barrel_lead_fail_bkgf->Write();
    h_pT_barrel_lead_fail_data->Write();
    h_pT_endcap_lead_fail_DY->Write();
    h_pT_endcap_lead_fail_bkgr->Write();
    h_pT_endcap_lead_fail_bkgf->Write();
    h_pT_endcap_lead_fail_data->Write();

    h_pT_barrel_sub_pass_DY->Write();
    h_pT_barrel_sub_pass_bkgr->Write();
    h_pT_barrel_sub_pass_bkgf->Write();
    h_pT_barrel_sub_pass_data->Write();
    h_pT_endcap_sub_pass_DY->Write();
    h_pT_endcap_sub_pass_bkgr->Write();
    h_pT_endcap_sub_pass_bkgf->Write();
    h_pT_endcap_sub_pass_data->Write();
    h_pT_barrel_sub_fail_DY->Write();
    h_pT_barrel_sub_fail_bkgr->Write();
    h_pT_barrel_sub_fail_bkgf->Write();
    h_pT_barrel_sub_fail_data->Write();
    h_pT_endcap_sub_fail_DY->Write();
    h_pT_endcap_sub_fail_bkgr->Write();
    h_pT_endcap_sub_fail_bkgf->Write();
    h_pT_endcap_sub_fail_data->Write();

    h_MET_fail_DY->Write();
    h_MET_fail_bkgr->Write();
    h_MET_fail_bkgf->Write();
    h_MET_fail_data->Write();
    h_MT_fail_DY->Write();
    h_MT_fail_bkgr->Write();
    h_MT_fail_bkgf->Write();
    h_MT_fail_data->Write();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir1+"For_FakeRate_electron_alt"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir1+"For_FakeRate_electron_alt"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FR_HistMaker_alt()


void E_FR_HistMaker_alt2 (Bool_t DEBUG) // UPDATE LATER: FR from 3rd "electron" in DY event
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


    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For PU re-weighting -- //
    analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

    // -- For PVz reweighting -- //
    analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

    for (Process_t pr=_DY_10to50; pr<_EndOf_DoubleEG_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_E_alt_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", nPtBin_ele_alt, analyzer->ptbin_ele_alt); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", nPtBin_ele_alt, analyzer->ptbin_ele_alt); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBin_ele_alt, analyzer->ptbin_ele_alt); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBin_ele_alt, analyzer->ptbin_ele_alt); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 50, etabins); h_eta_nume->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 50, etabins); h_eta_ctrl->Sumw2();
        TH1D* h_MET = new TH1D("h_MET", "h_MET", 100, 0, 100); h_MET->Sumw2();
        TH1D* h_MT = new TH1D("h_MT", "h_MT", 500, 0, 1000); h_MT->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_mass = new TH1D("h_mass", "h_mass", nMassBin, massbins); h_mass->Sumw2();

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
        chain->Add(Dir+"SelectedForFR_E_alt_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_E_alt_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

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
        Int_t nPass = 0;
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // FR selection
            if (p_T->size() != 3) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566) ||
                (fabs(etaSC->at(2)) > 1.4442 && fabs(etaSC->at(2)) < 1.566)) continue;
            if (fabs(etaSC->at(2)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(2) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(2)) >= 0.01 || fabs(dPhiIn->at(2)) >= 0.07 || mHits->at(2) > 1)) continue;
            if (fabs(etaSC->at(2)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(2) >= 0.035 || HoverE->at(2) >= 0.13 || mHits->at(2) > 1)) continue;
            if (p_T->at(0) < 28 || p_T->at(1) < 17 || p_T->at(2) < 15) continue;
            if (!passMediumID->at(0) || !passMediumID->at(1)) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso_Rho->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso_Rho->at(1) << endl;
            if (p_T->at(2) != p_T->at(2)) cout << p_T->at(2) << " " << eta->at(2) << " " << phi->at(2) << " " << charge->at(2) << " " << relPFiso_Rho->at(2) << endl;

//            if (MET_pT > 25) continue;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "pT size = " << p_T->size() << endl;
                cout << "eta size = " << eta->size() << endl;
                cout << "etaSC size = " << etaSC->size() << endl;
                cout << "phi size = " << phi->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "p_T[2] = " << p_T->at(2);
                cout << "\teta[2] = " << eta->at(2);
                cout << "\tphi[2] = " << phi->at(2) << endl;
                cout << "\tpassMediumID[2] = " << passMediumID->at(2) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            // -- Invariant mass and transverse mass -- //
            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();
            if (mass < 76.0 || mass > 106.0) {cout << "MASS = " << mass << endl; continue;}

            Double_t dTheta = phi->at(0) - MET_phi;
            Double_t MT = sqrt(2 * p_T->at(0) * MET_pT * (1 - cos(dTheta)));
//            if (MT < 30) continue;

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << "\tTop pT weight: " << TopPtWeight << endl;

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- Filling the histograms -- //
            h_mass->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
            h_MT->Fill(MT, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);

            if (passMediumID->at(2)) // Numerator
            {
                h_eta_nume->Fill(etaSC->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (etaSC->at(2) < 1.4442) // Barrel
                    h_pT_barrel_nume->Fill(p_T->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                else if (etaSC->at(2) > 1.566) // Endcap
                    h_pT_endcap_nume->Fill(p_T->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
            }
            else // Non-signal (failed numerator or Deno-nume)
            {
                h_eta_ctrl->Fill(etaSC->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                if (etaSC->at(2) < 1.4442) // Barrel
                    h_pT_barrel_ctrl->Fill(p_T->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
                else if (etaSC->at(2) > 1.566) // Endcap
                    h_pT_endcap_ctrl->Fill(p_T->at(2), TotWeight * PUWeight * PVzWeight * L1weight * TopPtWeight);
            }

            nPass++;

        }// End of event iteration

        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        cout << " Finished.\n" << endl;

        // -- Writing into files -- //
        f->cd();
        cout << "\tWriting into file...";

        h_pT_barrel_nume->Write();
        h_pT_endcap_nume->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_eta_nume->Write();
        h_eta_ctrl->Write();
        h_MET->Write();
        h_MT->Write();
        h_nVTX->Write();
        h_mass->Write();

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_E_alt_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_E_alt_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_50to100) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched_20to30

        cout << "===========================================================\n" << endl;
    } // End of pr iteration  

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FR_HistMaker_alt2()
