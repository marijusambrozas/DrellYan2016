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

void E_FRgraphMaker (Bool_t DEBUG);
void Mu_FRgraphMaker (Bool_t DEBUG);

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


void FRgraphMaker_TEST (TString WhichX = "")
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
//    if (whichX.Contains("E"))
//    {
//        Xselected++;
//        if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
//        else HLT = HLTname;
//        cout << "\n*******      E_FRgraphMaker      *******" << endl;
//        E_FRgraphMaker(DEBUG);
//    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*****  Mu_FRgraphMaker  *****" << endl;
        Mu_FRgraphMaker(DEBUG);
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FRgraphMaker_TEST.C+(\"whichX\")" << endl;

} // End of HistMaker()


/// ----------------------------- Electron Channel ------------------------------ ///
void E_FRgraphMaker (Bool_t DEBUG)
{
    return;
//    if (!type.Length())
//    {
//        cout << "Error: no type specified!" << endl;
//        return;
//    }

//    TTimeStamp ts_start;
//    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
//    TStopwatch totaltime;
//    totaltime.Start();

//    LocalFileMgr Mgr;
//    vector<SelProc_t> Processes = Mgr.FindProc(type);
//    TFile *f;
//    TString OutputDir;
//    TString debug = "";
//    if (DEBUG == kTRUE) debug = "_DEBUG";
//    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)

//    if (!Processes.size())
//    {
//        cout << "Error: no processes!" << endl;
//        return;
//    }

//    if (Processes[0] == _EE_Bkg_Full)
//    {
//        isBkgFull = kTRUE;
//        Mgr.SetProc(_EE_Bkg_Full);
//        // -- Output ROOTFile -- //
//        OutputDir = Mgr.HistLocation;
//        f = new TFile(OutputDir+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+debug+".root", "RECREATE");
//        Processes.clear();
//        Processes.push_back(_EE_DYTauTau_Full);
//        Processes.push_back(_EE_ttbar_Full);
//        Processes.push_back(_EE_tW);
//        Processes.push_back(_EE_tbarW);
//        Processes.push_back(_EE_WW);
//        Processes.push_back(_EE_WZ);
//        Processes.push_back(_EE_ZZ);
//        Processes.push_back(_EE_WJets_Full);
//        Processes.push_back(_EE_QCDEMEnriched_Full);
//    }

//    for (Int_t i_proc=0; i_proc<((int)(Processes.size())); i_proc++)
//    {
//        if (Processes[i_proc] <= _EndOf_MuMu || Processes[i_proc] > _EndOf_EE)
//        {
//            cout << "Error: process " << Mgr.Procname[Processes[i_proc]] << " is not EE!" << endl;
//            continue;
//        }

//        Mgr.SetProc(Processes[i_proc]);

//        // -- Output ROOTFile -- //
//        if (isBkgFull == kFALSE)
//        {
//            OutputDir = Mgr.HistLocation;
//            f = new TFile(OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");
//        }

//        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
//        cout << "Type: " << Mgr.Type << endl;
//        cout << "DATA location: " << Mgr.BaseLocation << endl;
//        cout << "Output directory: " << OutputDir << endl;

//        TStopwatch totaltime;
//        totaltime.Start();

//        DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

//        // -- For PU re-weighting -- //
//        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

//        // -- For efficiency SF -- //
//        analyzer->SetupEfficiencyScaleFactor_electron();

//        // -- For PVz reweighting -- //
////        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");
//        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");


//        // -- Creating Histograms -- //
//        TH1D *h_mass_before_PUCorr = new TH1D("h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass_before_EffCorr = new TH1D("h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass_before_PVzCorr = new TH1D("h_mass_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass_before_L1Corr = new TH1D("h_mass_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass_before_TopPtCorr = new TH1D("h_mass_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000);
//        TH1D *h_mass = new TH1D("h_mass_"+Mgr.Procname[Mgr.CurrentProc], "", binnum, massbins);
//        TH1D *h_mass2 = new TH1D("h_mass2_"+Mgr.Procname[Mgr.CurrentProc], "", binnum2, massbins2);
//        TH1D *h_pT_before_PUCorr = new TH1D("h_pT_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_before_EffCorr = new TH1D("h_pT_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_before_PVzCorr = new TH1D("h_pT_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_before_L1Corr = new TH1D("h_pT_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_before_TopPtCorr = new TH1D("h_pT_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT = new TH1D("h_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_rapi_before_PUCorr = new TH1D("h_rapi_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
//        TH1D *h_rapi_before_EffCorr = new TH1D("h_rapi_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
//        TH1D *h_rapi_before_PVzCorr = new TH1D("h_rapi_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
//        TH1D *h_rapi_before_L1Corr = new TH1D("h_rapi_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
//        TH1D *h_rapi_before_TopPtCorr = new TH1D("h_rapi_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
//        TH1D *h_rapi = new TH1D("h_rapi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);

//        TH1D *h_nVTX_before_PUCorr = new TH1D("h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75);
//        TH1D *h_nVTX_before_EffCorr = new TH1D("h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75);
//        TH1D *h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75);

//        TH1D *h_pT_lead_before_PUCorr = new TH1D("h_pT_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead_before_PUCorr = new TH1D("h_pT_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead_before_PUCorr = new TH1D("h_eta_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead_before_PUCorr = new TH1D("h_eta_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead_before_PUCorr = new TH1D("h_phi_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead_before_PUCorr = new TH1D("h_phi_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        TH1D *h_pT_lead_before_EffCorr = new TH1D("h_pT_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead_before_EffCorr = new TH1D("h_pT_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead_before_EffCorr = new TH1D("h_eta_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead_before_EffCorr = new TH1D("h_eta_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead_before_EffCorr = new TH1D("h_phi_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead_before_EffCorr = new TH1D("h_phi_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        TH1D *h_pT_lead_before_PVzCorr = new TH1D("h_pT_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead_before_PVzCorr = new TH1D("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead_before_PVzCorr = new TH1D("h_eta_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead_before_PVzCorr = new TH1D("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead_before_PVzCorr = new TH1D("h_phi_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead_before_PVzCorr = new TH1D("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        TH1D *h_pT_lead_before_L1Corr = new TH1D("h_pT_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead_before_L1Corr = new TH1D("h_pT_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead_before_L1Corr = new TH1D("h_eta_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead_before_L1Corr = new TH1D("h_eta_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead_before_L1Corr = new TH1D("h_phi_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead_before_L1Corr = new TH1D("h_phi_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        TH1D *h_pT_lead_before_TopPtCorr = new TH1D("h_pT_lead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead_before_TopPtCorr = new TH1D("h_pT_sublead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead_before_TopPtCorr = new TH1D("h_eta_lead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead_before_TopPtCorr = new TH1D("h_eta_sublead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead_before_TopPtCorr = new TH1D("h_phi_lead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead_before_TopPtCorr = new TH1D("h_phi_sublead_before_TopPtCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        TH1D *h_pT_lead = new TH1D("h_pT_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_pT_sublead = new TH1D("h_pT_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
//        TH1D *h_eta_lead = new TH1D("h_eta_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_eta_sublead = new TH1D("h_eta_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_lead = new TH1D("h_phi_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
//        TH1D *h_phi_sublead = new TH1D("h_phi_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

//        //Loop for all samples
//        const Int_t Ntup = Mgr.FileLocation.size();
//        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
//        {
//            TStopwatch looptime;
//            looptime.Start();

//            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

//            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
//            chain->Add(Mgr.FullLocation[i_tup]);

//            SelectedEE_t *EE = new SelectedEE_t();
//            EE->CreateFromChain(chain);

//            Int_t NEvents = chain->GetEntries();
//            if (NEvents != Mgr.nEvents[i_tup]) cout << "\tEvent numbers do not match!!!" << endl;
//            cout << "\t[Sum of weights:: " << Mgr.Wsum[i_tup] << "]" << endl;
//            cout << "\t[Selected Events: " << NEvents << "]" << endl;

//            if (DEBUG == kTRUE) NEvents = 10;
//            myProgressBar_t bar(NEvents);

//            for (Int_t i=0; i<NEvents; i++)
//            {
//                EE->GetEvent(i);

//                // -- Bring the weights -- //
//                Double_t GenWeight = EE->GENEvt_weight;
//                if (GenWeight > 1) GenWeight = 1;
//                else if (GenWeight < -1) GenWeight = -1;


//                // -- Pileup-Reweighting -- //
//                Double_t PUWeight = 1;
//                if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(EE->nPileUp);

//                // -- efficiency weights -- //
//                Double_t effweight = 1;

//                // -- PVz weights -- //
//                Double_t PVzWeight = 1;
//                if(Mgr.isMC == kTRUE && !Mgr.Tag[i_tup].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(EE->PVz);

//                // -- L1 prefiring weights -- //
//                Double_t L1weight = 1;
//                if (Mgr.isMC == kTRUE && !Mgr.Tag[i_tup].Contains("QCD")) L1weight = EE->_prefiringweight;
////                if (Mgr.isMC == kTRUE) L1weight = EE->_prefiringweightup;
////                if (Mgr.isMC == kTRUE) L1weight = EE->_prefiringweightdown;

//                // -- Top Pt weights -- //
//                Double_t TopPtWeight = 1;
//                if (Mgr.isMC == kTRUE && !Mgr.Tag[i_tup].Contains("QCD")) TopPtWeight = EE->_topPtWeight;

//                // -- Normalization -- //
//                Double_t TotWeight = GenWeight;
//                if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup]) * GenWeight;

//                if(EE->isSelPassed == kTRUE)
//                {
//                    TLorentzVector ele1, ele2;
//                    ele1.SetPtEtaPhiE(EE->Electron_pT->at(0), EE->Electron_eta->at(0), EE->Electron_phi->at(0), EE->Electron_Energy->at(0));
//                    ele2.SetPtEtaPhiE(EE->Electron_pT->at(1), EE->Electron_eta->at(1), EE->Electron_phi->at(1), EE->Electron_Energy->at(1));
//                    Double_t reco_Pt = (ele1 + ele2).Pt();
//                    Double_t reco_rapi = (ele1 + ele2).Rapidity();

//                    // -- Apply efficiency correcion -- //
//                    if(Mgr.isMC == kTRUE)
//                            effweight = analyzer->EfficiencySF_EventWeight_electron(EE);
////                    if(i%10000 == 0) cout << effweight << endl;

//                    h_mass_before_PUCorr->Fill(EE->Electron_InvM, TotWeight);
//                    h_mass_before_EffCorr->Fill(EE->Electron_InvM, TotWeight * PUWeight);
//                    h_mass_before_PVzCorr->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight);
//                    h_mass_before_L1Corr->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight);
//                    h_mass_before_TopPtCorr->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_mass_fine->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_mass->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_mass2->Fill(EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_pT_before_PUCorr->Fill(reco_Pt, TotWeight);
//                    h_pT_before_EffCorr->Fill(reco_Pt, TotWeight * PUWeight);
//                    h_pT_before_PVzCorr->Fill(reco_Pt, TotWeight * PUWeight * effweight);
//                    h_pT_before_L1Corr->Fill(reco_Pt, TotWeight * PUWeight * effweight * PVzWeight);
//                    h_pT_before_TopPtCorr->Fill(reco_Pt, TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_pT->Fill(reco_Pt, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_rapi_before_PUCorr->Fill(reco_rapi, TotWeight);
//                    h_rapi_before_EffCorr->Fill(reco_rapi, TotWeight * PUWeight);
//                    h_rapi_before_PVzCorr->Fill(reco_rapi, TotWeight * PUWeight * effweight);
//                    h_rapi_before_L1Corr->Fill(reco_rapi, TotWeight * PUWeight * effweight * PVzWeight);
//                    h_rapi_before_TopPtCorr->Fill(reco_rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_rapi->Fill(reco_rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

//                    h_nVTX_before_PUCorr->Fill(EE->nVertices, TotWeight);
//                    h_nVTX_before_EffCorr->Fill(EE->nVertices, TotWeight * PUWeight);
//                    h_nVTX->Fill(EE->nVertices, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

//                    int lead=0, sublead=1;
//                    if (EE->Electron_pT->at(0) < EE->Electron_pT->at(1))
//                    {
//                        lead = 1;
//                        sublead = 0;
//                    }

//                    h_pT_lead_before_PUCorr->Fill(EE->Electron_pT->at(lead), TotWeight);
//                    h_pT_sublead_before_PUCorr->Fill(EE->Electron_pT->at(sublead), TotWeight);
////                    h_eta_lead_before_PUCorr->Fill(EE->Electron_eta->at(lead), TotWeight);
////                    h_eta_sublead_before_PUCorr->Fill(EE->Electron_eta->at(sublead), TotWeight);
//                    h_eta_lead_before_PUCorr->Fill(EE->Electron_etaSC->at(lead), TotWeight);
//                    h_eta_sublead_before_PUCorr->Fill(EE->Electron_etaSC->at(sublead), TotWeight);
//                    h_phi_lead_before_PUCorr->Fill(EE->Electron_phi->at(lead), TotWeight);
//                    h_phi_sublead_before_PUCorr->Fill(EE->Electron_phi->at(sublead), TotWeight);

//                    h_pT_lead_before_EffCorr->Fill(EE->Electron_pT->at(lead), TotWeight * PUWeight);
//                    h_pT_sublead_before_EffCorr->Fill(EE->Electron_pT->at(sublead), TotWeight * PUWeight);
////                    h_eta_lead_before_EffCorr->Fill(EE->Electron_eta->at(lead), TotWeight * PUWeight);
////                    h_eta_sublead_before_EffCorr->Fill(EE->Electron_eta->at(sublead), TotWeight * PUWeight);
//                    h_eta_lead_before_EffCorr->Fill(EE->Electron_etaSC->at(lead), TotWeight * PUWeight);
//                    h_eta_sublead_before_EffCorr->Fill(EE->Electron_etaSC->at(sublead), TotWeight * PUWeight);
//                    h_phi_lead_before_EffCorr->Fill(EE->Electron_phi->at(lead), TotWeight * PUWeight);
//                    h_phi_sublead_before_EffCorr->Fill(EE->Electron_phi->at(sublead), TotWeight * PUWeight);

//                    h_pT_lead_before_PVzCorr->Fill(EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight);
//                    h_pT_sublead_before_PVzCorr->Fill(EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight);
////                    h_eta_lead_before_PVzCorr->Fill(EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight);
////                    h_eta_sublead_before_PVzCorr->Fill(EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight);
//                    h_eta_lead_before_PVzCorr->Fill(EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight);
//                    h_eta_sublead_before_PVzCorr->Fill(EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight);
//                    h_phi_lead_before_PVzCorr->Fill(EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight);
//                    h_phi_sublead_before_PVzCorr->Fill(EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight);

//                    h_pT_lead_before_L1Corr->Fill(EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight);
//                    h_pT_sublead_before_L1Corr->Fill(EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight);
////                    h_eta_lead_before_L1Corr->Fill(EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight);
////                    h_eta_sublead_before_L1Corr->Fill(EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight);
//                    h_eta_lead_before_L1Corr->Fill(EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight * PVzWeight);
//                    h_eta_sublead_before_L1Corr->Fill(EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight * PVzWeight);
//                    h_phi_lead_before_L1Corr->Fill(EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight);
//                    h_phi_sublead_before_L1Corr->Fill(EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight);

//                    h_pT_lead_before_TopPtCorr->Fill(EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_pT_sublead_before_TopPtCorr->Fill(EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
////                    h_eta_lead_before_TopPtCorr->Fill(EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
////                    h_eta_sublead_before_TopPtCorr->Fill(EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_eta_lead_before_TopPtCorr->Fill(EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_eta_sublead_before_TopPtCorr->Fill(EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_phi_lead_before_TopPtCorr->Fill(EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);
//                    h_phi_sublead_before_TopPtCorr->Fill(EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight);

//                    h_pT_lead->Fill(EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_pT_sublead->Fill(EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
////                    h_eta_lead->Fill(EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
////                    h_eta_sublead->Fill(EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_eta_lead->Fill(EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_eta_sublead->Fill(EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_phi_lead->Fill(EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//                    h_phi_sublead->Fill(EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

//                } // End of event selection
//                bar.Draw(i);

//            } // End of event iteration

//            if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);

//            Double_t LoopRunTime = looptime.CpuTime();
//            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

//        }// End of i_tup iteration

//        f->cd();
//        cout << "\tWriting into file...";

//        h_mass_before_PUCorr->Write();
//        h_mass_before_EffCorr->Write();
//        h_mass_before_PVzCorr->Write();
//        h_mass_before_L1Corr->Write();
//        h_mass_before_TopPtCorr->Write();
//        h_mass_fine->Write();
//        h_mass->Write();
//        h_mass2->Write();
//        Double_t err = 0;
//        Double_t integ = h_mass->IntegralAndError(1, h_mass->GetSize()-2, err);
//        cout << "Number of events in mass histogram: " << integ << " +- " << err << endl;

//        h_pT_before_PUCorr->Write();
//        h_pT_before_EffCorr->Write();
//        h_pT_before_PVzCorr->Write();
//        h_pT_before_L1Corr->Write();
//        h_pT_before_TopPtCorr->Write();
//        h_pT->Write();

//        h_rapi_before_PUCorr->Write();
//        h_rapi_before_EffCorr->Write();
//        h_rapi_before_PVzCorr->Write();
//        h_rapi_before_L1Corr->Write();
//        h_rapi_before_TopPtCorr->Write();
//        h_rapi->Write();

//        h_nVTX_before_PUCorr->Write();
//        h_nVTX_before_EffCorr->Write();
//        h_nVTX->Write();

//        h_pT_lead_before_PUCorr->Write();
//        h_pT_sublead_before_PUCorr->Write();
//        h_eta_lead_before_PUCorr->Write();
//        h_eta_sublead_before_PUCorr->Write();
//        h_phi_lead_before_PUCorr->Write();
//        h_phi_sublead_before_PUCorr->Write();

//        h_pT_lead_before_EffCorr->Write();
//        h_pT_sublead_before_EffCorr->Write();
//        h_eta_lead_before_EffCorr->Write();
//        h_eta_sublead_before_EffCorr->Write();
//        h_phi_lead_before_EffCorr->Write();
//        h_phi_sublead_before_EffCorr->Write();

//        h_pT_lead_before_PVzCorr->Write();
//        h_pT_sublead_before_PVzCorr->Write();
//        h_eta_lead_before_PVzCorr->Write();
//        h_eta_sublead_before_PVzCorr->Write();
//        h_phi_lead_before_PVzCorr->Write();
//        h_phi_sublead_before_PVzCorr->Write();

//        h_pT_lead_before_L1Corr->Write();
//        h_pT_sublead_before_L1Corr->Write();
//        h_eta_lead_before_L1Corr->Write();
//        h_eta_sublead_before_L1Corr->Write();
//        h_phi_lead_before_L1Corr->Write();
//        h_phi_sublead_before_L1Corr->Write();

//        h_pT_lead_before_TopPtCorr->Write();
//        h_pT_sublead_before_TopPtCorr->Write();
//        h_eta_lead_before_TopPtCorr->Write();
//        h_eta_sublead_before_TopPtCorr->Write();
//        h_phi_lead_before_TopPtCorr->Write();
//        h_phi_sublead_before_TopPtCorr->Write();

//        h_pT_lead->Write();
//        h_pT_sublead->Write();
//        h_eta_lead->Write();
//        h_eta_sublead->Write();
//        h_phi_lead->Write();
//        h_phi_sublead->Write();

//        cout << " Finished." << endl;

//    }// End of i_proc iteration

//    f->Close();
//    if (isBkgFull == kTRUE)
//    {
//        if (!f->IsOpen()) cout << "File Hist_" << Mgr.Procname[_EE_Bkg_Full]+debug << ".root has been closed successfully.\n" << endl;
//        else cout << "FILE Hist_" << Mgr.Procname[_EE_Bkg_Full]+debug << ".root COULD NOT BE CLOSED!\n" << endl;
//    }
//    else
//    {
//        if (!f->IsOpen()) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc]+debug << ".root has been closed successfully.\n" << endl;
//        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc]+debug << ".root COULD NOT BE CLOSED!\n" << endl;
//    }

//    Double_t TotalRunTime = totaltime.CpuTime();
//    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

//    TTimeStamp ts_end;
//    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FRgraphMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_FRgraphMaker (Bool_t DEBUG)
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

    for (Process_t pr=_DY_10to50; pr<_EndOf_SinglMuon_Normal; pr=next(pr))
    {
        // For counting how many events make it into the histograms (DEBUGGING)
        Int_t Npass_raw = 0;
        Double_t Npass_weighed = 0;

        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+"_TEST.root", "RECREATE");

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
//        analyzer->SetupEfficiencyScaleFactor_BtoF();
//        analyzer->SetupEfficiencyScaleFactor_GtoH();
//        analyzer->SetupEfficiencyScaleFactor_BtoF_new();
//        analyzer->SetupEfficiencyScaleFactor_GtoH_new();

        // -- For PVz reweighting -- //
//        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass", "h_mass", binnum, massbins); h_mass->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t evt_weight;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("evt_weight", &evt_weight);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            bar.Draw(i);
            if (DEBUG == kTRUE){
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
                cout << "eta[1] = " << eta->at(0) << endl;
                cout << "phi[1] = " << phi->at(0) << endl;
            }

            // Added test part that mimics DY->mumu selection
            if (p_T->size() != 2) continue;
            if (relPFiso->at(0) > 0.15 || relPFiso->at(1) > 0.15) continue;

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();


            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
//            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(MuMu->PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
//            if (Mgr.isMC == kTRUE) L1weight = MuMu->_prefiringweight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
//            if (Mgr.isMC == kTRUE) TopPtWeight = MuMu->_topPtWeight;

            // -- Normalization -- //
            Double_t TotWeight = evt_weight;
//                if(Mgr.isMC == kTRUE) TotWeight = (Lumi_GtoH * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup]) * GenWeight;
            if(Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * evt_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            Npass_raw++;
            Npass_weighed += TotWeight * PUWeight;
        }// End of event iteration

        cout << "\t *** " << Npass_raw << " raw events have passed the selection" << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** " << Npass_weighed << " normalised events have passed the selection" << endl;
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_nVTX->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+"_TEST.root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+"_TEST.root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MuMu_HistMaker()
