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

void MatrixMethod_ele (Bool_t DEBUG);
void MatrixMethod_ele_MC (Bool_t DEBUG);
void MatrixMethod_mu (Bool_t DEBUG);
void MatrixMethod_mu_MC (Bool_t DEBUG);
void MatrixMethod_singleMu_MC (Bool_t DEBUG);
void MatrixMethod_EMu (Bool_t DEBUG);

void MatrixMethod (TString WhichX = "", Bool_t DEBUG=kFALSE)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    if (whichX.Contains("DEBUG")) DEBUG = kTRUE;
    Int_t Xselected = 0;

    if (whichX.Contains("EMU"))
    {
        Xselected++;
        cout << "\n*****  MatrixMethod_EMu  *****" << endl;
        MatrixMethod_EMu(DEBUG);
    }
    else if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("MC"))
        {
            if (whichX.Contains("SINGLE"))
            {
                cout << "\n*****  MatrixMethod_singleMu_MC  *****" << endl;
                MatrixMethod_singleMu_MC(DEBUG);
            }
            else
            {
                cout << "\n*****  MatrixMethod_mu_MC  *****" << endl;
                MatrixMethod_mu_MC(DEBUG);
            }
        }
        else
        {
            cout << "\n*****  MatrixMethod_mu  *****" << endl;
            MatrixMethod_mu(DEBUG);
        }
    }
    else if (whichX.Contains("E"))
    {
        Xselected++;
        if (whichX.Contains("MC"))
        {
            cout << "\n*****  MatrixMethod_ele_MC  *****" << endl;
            MatrixMethod_ele_MC(DEBUG);
        }
        else
        {
            cout << "\n*****    MatrixMethod_ele    *****" << endl;
            MatrixMethod_ele(DEBUG);
        }
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FR_HistMaker.C+(\"whichX\", type)" << endl;

} // End of MatrixMethod()


/// -------------------------------- Electron Channel ------------------------------------ ///
void MatrixMethod_ele (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_E"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract");
    analyzer->SetupPRvalues_ele(Dir+"PromptRate_electron_alt.root", "subtract");

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
        TH1D* h_mass_PP = new TH1D("h_mass_PP_"+Mgr.Procname[pr], "h_mass_PP_"+Mgr.Procname[pr], binnum, massbins); h_mass_PP->Sumw2();
        TH1D* h_mass_PF = new TH1D("h_mass_PF_"+Mgr.Procname[pr], "h_mass_PF_"+Mgr.Procname[pr], binnum, massbins); h_mass_PF->Sumw2();
        TH1D* h_mass_FP = new TH1D("h_mass_FP_"+Mgr.Procname[pr], "h_mass_FP_"+Mgr.Procname[pr], binnum, massbins); h_mass_FP->Sumw2();
        TH1D* h_mass_FF = new TH1D("h_mass_FF_"+Mgr.Procname[pr], "h_mass_FF_"+Mgr.Procname[pr], binnum, massbins); h_mass_FF->Sumw2();
        TH1D* h_mass_PP_TT = new TH1D("h_mass_PP_TT_"+Mgr.Procname[pr], "h_mass_PP_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_PP_TT->Sumw2();
        TH1D* h_mass_PF_TT = new TH1D("h_mass_PF_TT_"+Mgr.Procname[pr], "h_mass_PF_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_PF_TT->Sumw2();
        TH1D* h_mass_FP_TT = new TH1D("h_mass_FP_TT_"+Mgr.Procname[pr], "h_mass_FP_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_FP_TT->Sumw2();
        TH1D* h_mass_FF_TT = new TH1D("h_mass_FF_TT_"+Mgr.Procname[pr], "h_mass_FF_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_FF_TT->Sumw2();
        TH1D* h_mass_TT = new TH1D("h_mass_TT_"+Mgr.Procname[pr], "h_mass_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_TT->Sumw2();
        TH1D* h_mass_TN = new TH1D("h_mass_TN_"+Mgr.Procname[pr], "h_mass_TN_"+Mgr.Procname[pr], binnum, massbins); h_mass_TN->Sumw2();
        TH1D* h_mass_NT = new TH1D("h_mass_NT_"+Mgr.Procname[pr], "h_mass_NT_"+Mgr.Procname[pr], binnum, massbins); h_mass_NT->Sumw2();
        TH1D* h_mass_NN = new TH1D("h_mass_NN_"+Mgr.Procname[pr], "h_mass_NN_"+Mgr.Procname[pr], binnum, massbins); h_mass_NN->Sumw2();

        TH1D* h_massSS_PP = new TH1D("h_massSS_PP_"+Mgr.Procname[pr], "h_massSS_PP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP->Sumw2();
        TH1D* h_massSS_PF = new TH1D("h_massSS_PF_"+Mgr.Procname[pr], "h_massSS_PF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF->Sumw2();
        TH1D* h_massSS_FP = new TH1D("h_massSS_FP_"+Mgr.Procname[pr], "h_massSS_FP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP->Sumw2();
        TH1D* h_massSS_FF = new TH1D("h_massSS_FF_"+Mgr.Procname[pr], "h_massSS_FF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF->Sumw2();
        TH1D* h_massSS_PP_TT = new TH1D("h_massSS_PP_TT_"+Mgr.Procname[pr], "h_massSS_PP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP_TT->Sumw2();
        TH1D* h_massSS_PF_TT = new TH1D("h_massSS_PF_TT_"+Mgr.Procname[pr], "h_massSS_PF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF_TT->Sumw2();
        TH1D* h_massSS_FP_TT = new TH1D("h_massSS_FP_TT_"+Mgr.Procname[pr], "h_massSS_FP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP_TT->Sumw2();
        TH1D* h_massSS_FF_TT = new TH1D("h_massSS_FF_TT_"+Mgr.Procname[pr], "h_massSS_FF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF_TT->Sumw2();

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
//            if (pr == _ttbar) continue;
            chain->GetEntry(i);
            for (UInt_t e=0; e<passMediumID->size(); e++)
            {
                if (passMediumID->at(e) == 1) nPassEle++;
                else nFailEle++;
            }
            if (!DEBUG) bar.Draw(i);

            // pre-selection
            if (p_T->size() != 2) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (fabs(etaSC->at(0)) >= 2.4 || fabs(etaSC->at(1)) >= 2.4) continue;
            if (fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) continue;
            if (fabs(etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566) continue;

            if (mHits->at(0) > 1 || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;

            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07)) continue;
            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07)) continue;
            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13)) continue;
            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13)) continue;

            // ALTERNATIVE FR (allowing only relPFiso to fail)
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
            if (scPixCharge->at(0) == scPixCharge->at(1) && isGsfCtfScPixConsistent->at(0) && isGsfCtfScPixConsistent->at(1)) nSameSignTight++;
            if (charge->at(0) == scPixCharge->at(0) && isGsfCtfScPixConsistent->at(0)) n_charge_tightCharge_match++;
            if (charge->at(1) == scPixCharge->at(1) && isGsfCtfScPixConsistent->at(1)) n_charge_tightCharge_match++;
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
            std::vector<double> MatrixMethod_weights(4, 1.0);
            Double_t FR1=0, FR2=0, PR1=1, PR2=1;
            if (p_T->at(0) >= p_T->at(1))
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight_ele(p_T->at(0), etaSC->at(0), passMediumID->at(0), p_T->at(1), etaSC->at(1), passMediumID->at(1));
                FR1 = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
                FR2 = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
                PR1 = analyzer->PromptRate_ele(p_T->at(0), etaSC->at(0));
                PR2 = analyzer->PromptRate_ele(p_T->at(1), etaSC->at(1));
            }
            else
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight_ele(p_T->at(1), etaSC->at(1), passMediumID->at(1), p_T->at(0), etaSC->at(0), passMediumID->at(0));
                FR1 = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
                FR2 = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
                PR1 = analyzer->PromptRate_ele(p_T->at(1), etaSC->at(1));
                PR2 = analyzer->PromptRate_ele(p_T->at(0), etaSC->at(0));
            }

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_mass_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
            h_mass_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
            h_mass_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
            h_mass_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
            h_mass_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PR1 * PR2);
            h_mass_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PR1 * FR2);
            h_mass_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FR1 * PR2);
            h_mass_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FR1 * FR2);
            if (charge->at(0) == charge->at(1))
            {
                h_massSS_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
                h_massSS_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
                h_massSS_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
                h_massSS_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
                h_massSS_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PR1 * PR2);
                h_massSS_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PR1 * FR2);
                h_massSS_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FR1 * PR2);
                h_massSS_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FR1 * FR2);
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

        h_mass_PP->Write();
        h_mass_PF->Write();
        h_mass_FP->Write();
        h_mass_FF->Write();
        h_mass_PP_TT->Write();
        h_mass_PF_TT->Write();
        h_mass_FP_TT->Write();
        h_mass_FF_TT->Write();
        h_massSS_PP->Write();
        h_massSS_PF->Write();
        h_massSS_FP->Write();
        h_massSS_FF->Write();
        h_massSS_PP_TT->Write();
        h_massSS_PF_TT->Write();
        h_massSS_FP_TT->Write();
        h_massSS_FF_TT->Write();

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

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_ele()


void MatrixMethod_ele_MC (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString Dir1 = "/media/sf_DownloadData/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_E_MC"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Photon_OR");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRandPRvalues_ele_MC();

    // Histogram containers
    TH1D* h_mass_PP[4];
    TH1D* h_mass_PF[4];
    TH1D* h_mass_FP[4];
    TH1D* h_mass_FF[4];
    TH1D* h_mass_TT[4];
    TH1D* h_mass_TN[4];
    TH1D* h_mass_NT[4];
    TH1D* h_mass_NN[4];
    TH1D* h_mass_PP_true[4];
    TH1D* h_mass_PF_true[4];
    TH1D* h_mass_FP_true[4];
    TH1D* h_mass_FF_true[4];

    TString names[4] = {"DY", "ttbar", "WJets", "QCD"};
    // -- Creating Histograms -- //
    for (Int_t i_type=0; i_type<4; i_type++)
    {
        h_mass_PP[i_type] = new TH1D("h_mass_PP_"+names[i_type], "", binnum, massbins);
        h_mass_PF[i_type] = new TH1D("h_mass_PF_"+names[i_type], "", binnum, massbins);
        h_mass_FP[i_type] = new TH1D("h_mass_FP_"+names[i_type], "", binnum, massbins);
        h_mass_FF[i_type] = new TH1D("h_mass_FF_"+names[i_type], "", binnum, massbins);
        h_mass_TT[i_type] = new TH1D("h_mass_TT_"+names[i_type], "", binnum, massbins);
        h_mass_TN[i_type] = new TH1D("h_mass_TN_"+names[i_type], "", binnum, massbins);
        h_mass_NT[i_type] = new TH1D("h_mass_NT_"+names[i_type], "", binnum, massbins);
        h_mass_NN[i_type] = new TH1D("h_mass_NN_"+names[i_type], "", binnum, massbins);
        h_mass_PP_true[i_type] = new TH1D("h_mass_PP_true_"+names[i_type], "", binnum, massbins);
        h_mass_PF_true[i_type] = new TH1D("h_mass_PF_true_"+names[i_type], "", binnum, massbins);
        h_mass_FP_true[i_type] = new TH1D("h_mass_FP_true_"+names[i_type], "", binnum, massbins);
        h_mass_FF_true[i_type] = new TH1D("h_mass_FF_true_"+names[i_type], "", binnum, massbins);
    }
    Int_t i_type = 0;

    Process_t pr_first = _DY_10to50, pr_last = _EndOf_QCDEMEnriched_Normal;
    if (DEBUG) pr_first = _ttbar, pr_last = _ttbar_700to1000;
    for (Process_t pr=pr_first; pr<pr_last; pr=next(pr))
    {
        Mgr.SetProc(pr);

        if (pr < _EndOf_DY_Normal) i_type = 0;
        else if (pr < _EndOf_ttbar_Normal) i_type = 1;
        else if (pr < _EndOf_WJets_Normal) i_type = 2;
        else if (pr < _EndOf_QCDEMEnriched_Normal) i_type = 3;

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

        Int_t nElectrons;
        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
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
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
        std::vector<double> *relTrkIso = new std::vector<double>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
        Int_t nGenElectrons;
        std::vector<double> *genEle_pT = new std::vector<double>;
        std::vector<double> *genEle_eta = new std::vector<double>;
        std::vector<double> *genEle_phi = new std::vector<double>;
        std::vector<int> *genEle_charge = new std::vector<int>;
        std::vector<int> *genEle_isPrompt = new std::vector<int>;
        std::vector<int> *genEle_isPromptFinalState = new std::vector<int>;
        std::vector<int> *genEle_isTauDecayProduct = new std::vector<int>;
        std::vector<int> *genEle_isPromptTauDecayProduct = new std::vector<int>;
        std::vector<int> *genEle_isDirectPromptTauDecayProductFinalState = new std::vector<int>;
        std::vector<int> *genEle_isHardProcess = new std::vector<int>;
        std::vector<int> *genEle_isLastCopy = new std::vector<int>;
        std::vector<int> *genEle_isLastCopyBeforeFSR = new std::vector<int>;
        std::vector<int> *genEle_fromHardProcessBeforeFSR = new std::vector<int>;
        std::vector<int> *genEle_fromHardProcessFinalState = new std::vector<int>;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir1+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir1+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("nElectrons", 1);
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
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
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("relECALiso", 1);
        chain->SetBranchStatus("relHCALiso", 1);
        chain->SetBranchStatus("relTrkIso", 1);
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
        chain->SetBranchStatus("nGenElectrons", 1);
        chain->SetBranchStatus("genEle_pT", 1);
        chain->SetBranchStatus("genEle_eta", 1);
        chain->SetBranchStatus("genEle_phi", 1);
        chain->SetBranchStatus("genEle_charge", 1);
        chain->SetBranchStatus("genEle_isPrompt", 1);
        chain->SetBranchStatus("genEle_isPromptFinalState", 1);
        chain->SetBranchStatus("genEle_isTauDecayProduct", 1);
        chain->SetBranchStatus("genEle_isPromptTauDecayProduct", 1);
        chain->SetBranchStatus("genEle_isDirectPromptTauDecayProductFinalState", 1);
        chain->SetBranchStatus("genEle_isHardProcess", 1);
        chain->SetBranchStatus("genEle_isLastCopy", 1);
        chain->SetBranchStatus("genEle_isLastCopyBeforeFSR", 1);
        chain->SetBranchStatus("genEle_fromHardProcessBeforeFSR", 1);
        chain->SetBranchStatus("genEle_fromHardProcessFinalState", 1);

        chain->SetBranchAddress("nElectrons", &nElectrons);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("scPixCharge", &scPixCharge);
        chain->SetBranchAddress("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        chain->SetBranchAddress("isGsfScPixConsistent", &isGsfScPixConsistent);
        chain->SetBranchAddress("isGsfCtfConsistent", &isGsfCtfConsistent);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("dEtaInSeed", dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("InvEminusInvP", &InvEminusInvP);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passConvVeto", &passConvVeto);
        chain->SetBranchAddress("relPFiso_Rho", &relPFiso_Rho);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("relECALiso", &relECALiso);
        chain->SetBranchAddress("relHCALiso", &relHCALiso);
        chain->SetBranchAddress("relTrkIso", &relTrkIso);
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
        chain->SetBranchAddress("nGenElectrons", &nGenElectrons);
        chain->SetBranchAddress("genEle_pT", &genEle_pT);
        chain->SetBranchAddress("genEle_eta", &genEle_eta);
        chain->SetBranchAddress("genEle_phi", &genEle_phi);
        chain->SetBranchAddress("genEle_charge", &genEle_charge);
        chain->SetBranchAddress("genEle_isPrompt", &genEle_isPrompt);
        chain->SetBranchAddress("genEle_isPromptFinalState", &genEle_isPromptFinalState);
        chain->SetBranchAddress("genEle_isTauDecayProduct", &genEle_isTauDecayProduct);
        chain->SetBranchAddress("genEle_isPromptTauDecayProduct", &genEle_isPromptTauDecayProduct);
        chain->SetBranchAddress("genEle_isDirectPromptTauDecayProductFinalState", &genEle_isDirectPromptTauDecayProductFinalState);
        chain->SetBranchAddress("genEle_isHardProcess", &genEle_isHardProcess);
        chain->SetBranchAddress("genEle_isLastCopy", &genEle_isLastCopy);
        chain->SetBranchAddress("genEle_isLastCopyBeforeFSR", &genEle_isLastCopyBeforeFSR);
        chain->SetBranchAddress("genEle_fromHardProcessBeforeFSR", &genEle_fromHardProcessBeforeFSR);
        chain->SetBranchAddress("genEle_fromHardProcessFinalState", &genEle_fromHardProcessFinalState);

        Int_t NEvents = chain->GetEntries();
        if (DEBUG) NEvents = 1000;
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassEle=0, nFailEle=0;
        UInt_t nSameSign=0;
        UInt_t nEle=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            for (UInt_t ele=0; ele<passMediumID->size(); ele++)
            {
                if (passMediumID->at(ele)) nPassEle++;
                else nFailEle++;
            }
            if (!DEBUG) bar.Draw(i);

            // pre-selection
            if (p_T->size() != 2) continue;
//            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 25 || p_T->at(1) < 25) continue;
            if (fabs(eta->at(0)) >= 2.4 || fabs(eta->at(1)) >= 2.4) continue;
            if (fabs(eta->at(0)) >= 1.4442 && fabs(eta->at(0)) <= 1.566) continue;
            if (fabs(eta->at(1)) >= 1.4442 && fabs(eta->at(1)) <= 1.566) continue;
            if (mHits->at(0) > 1 || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;
            if (fabs(eta->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                              fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07)) continue;
            else if (fabs(eta->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13)) continue;
            if (fabs(eta->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                              fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07)) continue;
            else if (fabs(eta->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13)) continue;

            if (charge->at(0) == charge->at(1)) nSameSign++;
            nEle += 2;

            if (DEBUG)
            {
                if (nPass >= 10) break;
                cout << "Evt " << i << endl;
            }

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << passMediumID->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << passMediumID->at(1) << endl;

            Int_t i_lead = (p_T->at(0) >= p_T->at(1)) ? 0 : 1;
            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
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

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            std::vector<int> matches_lead;
            std::vector<int> matches_sub;
            Double_t dRmin_lead = 9999;
            Double_t dRmin_sub = 9999;
            Double_t rel_dpTmin_lead = 9999;
            Double_t rel_dpTmin_sub = 9999;
            Int_t i_bestmatch_lead = -1;
            Int_t i_bestmatch_sub = -1;
            Int_t veto_lead = 0; // used to find fake leptons
            Int_t veto_sub = 0; // used to find fake leptons
            if (DEBUG) cout << "  nGenElectrons = " << genEle_pT->size() << endl;
            for (UInt_t i_gen=0; i_gen<genEle_pT->size(); i_gen++)
            {
                if (DEBUG)
                {
                    cout << "   GenElectron" << i_gen << ":  p_T=" << genEle_pT->at(i_gen);
                    cout << "  eta=" << genEle_eta->at(i_gen) << "  phi=" << genEle_phi->at(i_gen) << endl;
                    cout << "    isPrompt=" << genEle_isPrompt->at(i_gen) << "  isLastCopy=" << genEle_isLastCopy->at(i_gen);
                    cout << "  iLCbeforeFSR=" << genEle_isLastCopyBeforeFSR->at(i_gen) << endl;
                }

                Double_t dEta_lead = eta->at(i_lead) - genEle_eta->at(i_gen);
                Double_t dPhi_lead = phi->at(i_lead) - genEle_phi->at(i_gen);
                Double_t dR_lead = sqrt(dEta_lead*dEta_lead + dPhi_lead*dPhi_lead);
                Double_t rel_dpT_lead = fabs(p_T->at(i_lead) - genEle_pT->at(i_gen)) / p_T->at(i_lead);

                Double_t dEta_sub = eta->at(1-i_lead) - genEle_eta->at(i_gen);
                Double_t dPhi_sub = phi->at(1-i_lead) - genEle_phi->at(i_gen);
                Double_t dR_sub = sqrt(dEta_sub*dEta_sub + dPhi_sub*dPhi_sub);
                Double_t rel_dpT_sub = fabs(p_T->at(1-i_lead) - genEle_pT->at(i_gen)) / p_T->at(1-i_lead);
                if (dR_lead < 0.3)
                {
                    if (genEle_isPrompt->at(i_gen) || genEle_isPromptFinalState->at(i_gen) || genEle_isTauDecayProduct->at(i_gen) ||
                        genEle_isPromptTauDecayProduct->at(i_gen) || genEle_isDirectPromptTauDecayProductFinalState->at(i_gen) ||
                        genEle_isHardProcess->at(i_gen) || genEle_fromHardProcessBeforeFSR->at(i_gen) || genEle_fromHardProcessFinalState->at(i_gen))
                        veto_lead = 1;
                    if (rel_dpT_lead < 0.5)
                    {
                        matches_lead.push_back(i_gen);
                        if (genEle_isLastCopy->at(i_gen) && (genEle_isPrompt->at(i_gen) || genEle_isPromptFinalState->at(i_gen)))
                        {
                            if (dR_lead < dRmin_lead || (dR_lead == dRmin_lead && rel_dpT_lead <= rel_dpTmin_lead))
                            {
                                dRmin_lead = dR_lead;
                                rel_dpTmin_lead = rel_dpT_lead;
                                i_bestmatch_lead = i_gen;
                            }
                        }
                    }
                }
                if (dR_sub < 0.3)
                {
                    if (genEle_isPrompt->at(i_gen) || genEle_isPromptFinalState->at(i_gen) || genEle_isTauDecayProduct->at(i_gen) ||
                        genEle_isPromptTauDecayProduct->at(i_gen) || genEle_isDirectPromptTauDecayProductFinalState->at(i_gen) ||
                        genEle_isHardProcess->at(i_gen) || genEle_fromHardProcessBeforeFSR->at(i_gen) || genEle_fromHardProcessFinalState->at(i_gen))
                        veto_sub = 1;
                    if (rel_dpT_sub < 0.5)
                    {
                        matches_sub.push_back(i_gen);
                        if (genEle_isLastCopy->at(i_gen) && (genEle_isPrompt->at(i_gen) || genEle_isPromptFinalState->at(i_gen)))
                        {
                            if (dR_sub < dRmin_sub || (dR_sub == dRmin_sub && rel_dpT_sub <= rel_dpTmin_sub))
                            {
                                dRmin_sub = dR_sub;
                                rel_dpTmin_sub = rel_dpT_sub;
                                i_bestmatch_sub = i_gen;
                            }
                        }
                    }
                }
            }

            if (DEBUG)
            {
                cout << "\n  nElectrons = " << p_T->size() << endl;
                cout << "   p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0);
                Int_t printval = (i_lead == 0)?i_bestmatch_lead:i_bestmatch_sub;
                cout << "\tMatched to " << printval;
                if ((i_lead == 0 && veto_lead) || (i_lead != 0 && veto_sub))
                    cout << "\tVETOED";
                cout << endl;
                cout << "   p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1);
                printval = (i_lead == 1)?i_bestmatch_lead:i_bestmatch_sub;
                cout << "\tMatched to " << printval;
                if ((i_lead == 1 && veto_lead) || (i_lead != 1 && veto_sub))
                    cout << "\tVETOED";
                cout << endl;
            }

            if ((i_bestmatch_lead < 0 && veto_lead) || (i_bestmatch_sub < 0 && veto_sub)) continue;
            nPass++;

            if (i_bestmatch_lead >= 0 && i_bestmatch_sub >= 0)
                h_mass_PP_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if (i_bestmatch_lead >= 0 && (i_bestmatch_sub < 0 && !veto_sub))
                h_mass_PF_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if ((i_bestmatch_lead < 0 && !veto_lead) && i_bestmatch_sub >= 0)
                h_mass_FP_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if ((i_bestmatch_lead < 0 && !veto_lead) && (i_bestmatch_sub < 0 && !veto_sub))
                h_mass_FF_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            // -- FR WEIGHTS -- //
            std::vector<double> MMweights(4, 1.0);
            Double_t FR1=0, FR2=0, PR1=1, PR2=1;

            MMweights = analyzer->MatrixMethod_ele_MC(p_T->at(i_lead), eta->at(i_lead), passMediumID->at(i_lead),
                                                      p_T->at(1-i_lead), eta->at(1-i_lead), passMediumID->at(1-i_lead), "DY", "QCD");
            FR1 = analyzer->FakeRate_ele_MC(p_T->at(i_lead), eta->at(i_lead), "QCD");
            FR2 = analyzer->FakeRate_ele_MC(p_T->at(1-i_lead), eta->at(1-i_lead), "QCD");
            PR1 = analyzer->PromptRate_ele_MC(p_T->at(i_lead), eta->at(i_lead), "DY");
            PR2 = analyzer->PromptRate_ele_MC(p_T->at(1-i_lead), eta->at(1-i_lead), "DY");
//            MMweights = analyzer->MatrixMethod_ele_eta_MC(eta->at(i_lead), passMediumID->at(i_lead),
//                                                          eta->at(1-i_lead), passMediumID->at(1-i_lead), "ttbar", "ttbar");
//            FR1 = analyzer->FakeRate_ele_eta_MC(eta->at(i_lead), "QCD");
//            FR2 = analyzer->FakeRate_ele_eta_MC(eta->at(1-i_lead), "QCD");
//            PR1 = analyzer->PromptRate_ele_eta_MC(eta->at(i_lead), "DY");
//            PR2 = analyzer->PromptRate_ele_eta_MC(eta->at(1-i_lead), "DY");

            // Filling histos
            h_mass_PP[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[0]);
            h_mass_PF[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[1]);
            h_mass_FP[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[2]);
            h_mass_FF[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[3]);

            if(passMediumID->at(i_lead) && passMediumID->at(1-i_lead))
                h_mass_TT[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(passMediumID->at(i_lead) && !passMediumID->at(1-i_lead))
                h_mass_TN[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(!passMediumID->at(i_lead) && passMediumID->at(1-i_lead))
                h_mass_NT[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(!passMediumID->at(i_lead) && !passMediumID->at(1-i_lead))
                h_mass_NN[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

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

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _ttbar_1000toInf) pr = _EndOf_VVnST_Normal; // next -- WJets
        if (pr == _EndOf_WJets_Normal) pr = _QCDEMEnriched_20to30; // next -- QCDEMEnriched

        cout << "nSameSign = " << nSameSign << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    for (Int_t i=0; i<4; i++)
    {
        h_mass_PP[i]->Write();
        h_mass_PF[i]->Write();
        h_mass_FP[i]->Write();
        h_mass_FF[i]->Write();
        h_mass_TT[i]->Write();
        h_mass_TN[i]->Write();
        h_mass_NT[i]->Write();
        h_mass_NN[i]->Write();
        h_mass_PP_true[i]->Write();
        h_mass_PF_true[i]->Write();
        h_mass_FP_true[i]->Write();
        h_mass_FF_true[i]->Write();
    }

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_ele_MC()


/// -------------------------------- Muon Channel ------------------------------------ ///
void MatrixMethod_mu (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_Mu"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", "sigCtrl_template");
    analyzer->SetupPRvalues(Dir+"PromptRate_muon_alt.root", "subtract", 1);

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
        TH1D* h_mass_PP = new TH1D("h_mass_PP_"+Mgr.Procname[pr], "h_mass_PP_"+Mgr.Procname[pr], binnum, massbins); h_mass_PP->Sumw2();
        TH1D* h_mass_PF = new TH1D("h_mass_PF_"+Mgr.Procname[pr], "h_mass_PF_"+Mgr.Procname[pr], binnum, massbins); h_mass_PF->Sumw2();
        TH1D* h_mass_FP = new TH1D("h_mass_FP_"+Mgr.Procname[pr], "h_mass_FP_"+Mgr.Procname[pr], binnum, massbins); h_mass_FP->Sumw2();
        TH1D* h_mass_FF = new TH1D("h_mass_FF_"+Mgr.Procname[pr], "h_mass_FF_"+Mgr.Procname[pr], binnum, massbins); h_mass_FF->Sumw2();
        TH1D* h_mass_PP_TT = new TH1D("h_mass_PP_TT_"+Mgr.Procname[pr], "h_mass_PP_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_PP_TT->Sumw2();
        TH1D* h_mass_PF_TT = new TH1D("h_mass_PF_TT_"+Mgr.Procname[pr], "h_mass_PF_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_PF_TT->Sumw2();
        TH1D* h_mass_FP_TT = new TH1D("h_mass_FP_TT_"+Mgr.Procname[pr], "h_mass_FP_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_FP_TT->Sumw2();
        TH1D* h_mass_FF_TT = new TH1D("h_mass_FF_TT_"+Mgr.Procname[pr], "h_mass_FF_TT_"+Mgr.Procname[pr], binnum, massbins); h_mass_FF_TT->Sumw2();

        TH1D* h_massSS_PP = new TH1D("h_massSS_PP_"+Mgr.Procname[pr], "h_massSS_PP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP->Sumw2();
        TH1D* h_massSS_PF = new TH1D("h_massSS_PF_"+Mgr.Procname[pr], "h_massSS_PF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF->Sumw2();
        TH1D* h_massSS_FP = new TH1D("h_massSS_FP_"+Mgr.Procname[pr], "h_massSS_FP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP->Sumw2();
        TH1D* h_massSS_FF = new TH1D("h_massSS_FF_"+Mgr.Procname[pr], "h_massSS_FF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF->Sumw2();
        TH1D* h_massSS_PP_TT = new TH1D("h_massSS_PP_TT_"+Mgr.Procname[pr], "h_massSS_PP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP_TT->Sumw2();
        TH1D* h_massSS_PF_TT = new TH1D("h_massSS_PF_TT_"+Mgr.Procname[pr], "h_massSS_PF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF_TT->Sumw2();
        TH1D* h_massSS_FP_TT = new TH1D("h_massSS_FP_TT_"+Mgr.Procname[pr], "h_massSS_FP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP_TT->Sumw2();
        TH1D* h_massSS_FF_TT = new TH1D("h_massSS_FF_TT_"+Mgr.Procname[pr], "h_massSS_FF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF_TT->Sumw2();

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
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

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

        UInt_t nPassMu=0, nFailMu=0;
        UInt_t nSameSign=0;
        UInt_t nMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
//            if (pr == _ttbar) continue;
            chain->GetEntry(i);
            for (UInt_t mu=0; mu<relPFiso->size(); mu++)
            {
                if (relPFiso->at(mu) < 0.15) nPassMu++;
                else nFailMu++;
            }
            if (!DEBUG) bar.Draw(i);

            // pre-selection
            if (p_T->size() != 2) continue;
//            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 52 || p_T->at(1) < 52) continue;
            if (fabs(eta->at(0)) >= 2.4 || fabs(eta->at(1)) >= 2.4) continue;

            if (charge->at(0) == charge->at(1)) nSameSign++;
            nMu += 2;

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
                cout << "\trelPFiso[0] = " << relPFiso->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\trelPFiso[1] = " << relPFiso->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();
            Double_t pTll = (mu1+mu2).Pt();
            Double_t etall = (mu1+mu2).Eta();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;
            // -- Efficiency scale factor -- //
//            if(Mgr.isMC == kTRUE)
//            {
//                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(MuMu);
//                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(MuMu);
//                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
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
            std::vector<double> MatrixMethod_weights(4, 1.0);
            Double_t FR1=0, FR2=0, PR1=1, PR2=1;
            if (p_T->at(0) >= p_T->at(1))
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight(p_T->at(0), eta->at(0), relPFiso->at(0), p_T->at(1), eta->at(1), relPFiso->at(1), 1);
                FR1 = analyzer->FakeRate(p_T->at(0), eta->at(0));
                FR2 = analyzer->FakeRate(p_T->at(1), eta->at(1));
                PR1 = analyzer->PromptRate(p_T->at(0), eta->at(0), 1);
                PR2 = analyzer->PromptRate(p_T->at(1), eta->at(1), 1);
            }
            else
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight(p_T->at(1), eta->at(1), relPFiso->at(1), p_T->at(0), eta->at(0), relPFiso->at(0), 1);
                FR1 = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FR2 = analyzer->FakeRate(p_T->at(0), eta->at(0));
                PR1 = analyzer->PromptRate(p_T->at(1), eta->at(1), 1);
                PR2 = analyzer->PromptRate(p_T->at(0), eta->at(0), 1);
            }

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            if (charge->at(0) != charge->at(1))
            {
                h_mass_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
                h_mass_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
                h_mass_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
                h_mass_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
                h_mass_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PR1 * PR2);
                h_mass_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PR1 * FR2);
                h_mass_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FR1 * PR2);
                h_mass_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FR1 * FR2);
            }
            else
            {
                h_massSS_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
                h_massSS_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
                h_massSS_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
                h_massSS_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
                h_massSS_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PR1 * PR2);
                h_massSS_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PR1 * FR2);
                h_massSS_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FR1 * PR2);
                h_massSS_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FR1 * FR2);
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

        h_mass_PP->Write();
        h_mass_PF->Write();
        h_mass_FP->Write();
        h_mass_FF->Write();
        h_mass_PP_TT->Write();
        h_mass_PF_TT->Write();
        h_mass_FP_TT->Write();
        h_mass_FF_TT->Write();
        h_massSS_PP->Write();
        h_massSS_PF->Write();
        h_massSS_FP->Write();
        h_massSS_FF->Write();
        h_massSS_PP_TT->Write();
        h_massSS_PF_TT->Write();
        h_massSS_FP_TT->Write();
        h_massSS_FF_TT->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "nSameSign = " << nSameSign << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_mu()


void MatrixMethod_mu_MC (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_Mu_MC"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRandPRvalues_MC("NONE");

    // Histogram containers
    TH1D* h_mass_PP[4];
    TH1D* h_mass_PF[4];
    TH1D* h_mass_FP[4];
    TH1D* h_mass_FF[4];
    TH1D* h_mass_TT[4];
    TH1D* h_mass_TN[4];
    TH1D* h_mass_NT[4];
    TH1D* h_mass_NN[4];
    TH1D* h_mass_PP_true[4];
    TH1D* h_mass_PF_true[4];
    TH1D* h_mass_FP_true[4];
    TH1D* h_mass_FF_true[4];

    TString names[4] = {"DY", "ttbar", "WJets", "QCD"};
    // -- Creating Histograms -- //
    for (Int_t i_type=0; i_type<4; i_type++)
    {
        h_mass_PP[i_type] = new TH1D("h_mass_PP_"+names[i_type], "", binnum, massbins);
        h_mass_PF[i_type] = new TH1D("h_mass_PF_"+names[i_type], "", binnum, massbins);
        h_mass_FP[i_type] = new TH1D("h_mass_FP_"+names[i_type], "", binnum, massbins);
        h_mass_FF[i_type] = new TH1D("h_mass_FF_"+names[i_type], "", binnum, massbins);
        h_mass_TT[i_type] = new TH1D("h_mass_TT_"+names[i_type], "", binnum, massbins);
        h_mass_TN[i_type] = new TH1D("h_mass_TN_"+names[i_type], "", binnum, massbins);
        h_mass_NT[i_type] = new TH1D("h_mass_NT_"+names[i_type], "", binnum, massbins);
        h_mass_NN[i_type] = new TH1D("h_mass_NN_"+names[i_type], "", binnum, massbins);
        h_mass_PP_true[i_type] = new TH1D("h_mass_PP_true_"+names[i_type], "", binnum, massbins);
        h_mass_PF_true[i_type] = new TH1D("h_mass_PF_true_"+names[i_type], "", binnum, massbins);
        h_mass_FP_true[i_type] = new TH1D("h_mass_FP_true_"+names[i_type], "", binnum, massbins);
        h_mass_FF_true[i_type] = new TH1D("h_mass_FF_true_"+names[i_type], "", binnum, massbins);
    }
    Int_t i_type = 0;

    Process_t pr_first = _DY_10to50, pr_last = _EndOf_QCDMuEnriched_Normal;
    if (DEBUG) pr_first = _ttbar, pr_last = _ttbar_700to1000;
    for (Process_t pr=pr_first; pr<pr_last; pr=next(pr))
    {
        Mgr.SetProc(pr);

        if (pr < _EndOf_DY_Normal) i_type = 0;
        else if (pr < _EndOf_ttbar_Normal) i_type = 1;
        else if (pr < _EndOf_WJets_Normal) i_type = 2;
        else if (pr < _EndOf_QCDMuEnriched_Normal) i_type = 3;

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

        Int_t nMuons;
        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nGenMuons;
        std::vector<double> *genMu_pT = new std::vector<double>;
        std::vector<double> *genMu_eta = new std::vector<double>;
        std::vector<double> *genMu_phi = new std::vector<double>;
        std::vector<int> *genMu_charge = new std::vector<int>;
        std::vector<int> *genMu_isPrompt = new std::vector<int>;
        std::vector<int> *genMu_isPromptFinalState = new std::vector<int>;
        std::vector<int> *genMu_isTauDecayProduct = new std::vector<int>;
        std::vector<int> *genMu_isPromptTauDecayProduct = new std::vector<int>;
        std::vector<int> *genMu_isDirectPromptTauDecayProductFinalState = new std::vector<int>;
        std::vector<int> *genMu_isHardProcess = new std::vector<int>;
        std::vector<int> *genMu_isLastCopy = new std::vector<int>;
        std::vector<int> *genMu_isLastCopyBeforeFSR = new std::vector<int>;
        std::vector<int> *genMu_isPromptDecayed = new std::vector<int>;
        std::vector<int> *genMu_isDecayedLeptonHadron = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessBeforeFSR = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessDecayed = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessFinalState = new std::vector<int>;
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

        chain->SetBranchStatus("nMuons", 1);
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
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchStatus("nGenMuons", 1);
        chain->SetBranchStatus("genMu_pT", 1);
        chain->SetBranchStatus("genMu_eta", 1);
        chain->SetBranchStatus("genMu_phi", 1);
        chain->SetBranchStatus("genMu_charge", 1);
        chain->SetBranchStatus("genMu_isPrompt", 1);
        chain->SetBranchStatus("genMu_isPromptFinalState", 1);
        chain->SetBranchStatus("genMu_isTauDecayProduct", 1);
        chain->SetBranchStatus("genMu_isPromptTauDecayProduct", 1);
        chain->SetBranchStatus("genMu_isDirectPromptTauDecayProductFinalState", 1);
        chain->SetBranchStatus("genMu_isHardProcess", 1);
        chain->SetBranchStatus("genMu_isLastCopy", 1);
        chain->SetBranchStatus("genMu_isLastCopyBeforeFSR", 1);
        chain->SetBranchStatus("genMu_isPromptDecayed", 1);
        chain->SetBranchStatus("genMu_isDecayedLeptonHadron", 1);
        chain->SetBranchStatus("genMu_fromHardProcessBeforeFSR", 1);
        chain->SetBranchStatus("genMu_fromHardProcessDecayed", 1);
        chain->SetBranchStatus("genMu_fromHardProcessFinalState", 1);

        chain->SetBranchAddress("nMuons", &nMuons);
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
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);
        chain->SetBranchAddress("genMu_pT", &genMu_pT);
        chain->SetBranchAddress("genMu_eta", &genMu_eta);
        chain->SetBranchAddress("genMu_phi", &genMu_phi);
        chain->SetBranchAddress("genMu_charge", &genMu_charge);
        chain->SetBranchAddress("genMu_isPrompt", &genMu_isPrompt);
        chain->SetBranchAddress("genMu_isPromptFinalState", &genMu_isPromptFinalState);
        chain->SetBranchAddress("genMu_isTauDecayProduct", &genMu_isTauDecayProduct);
        chain->SetBranchAddress("genMu_isPromptTauDecayProduct", &genMu_isPromptTauDecayProduct);
        chain->SetBranchAddress("genMu_isDirectPromptTauDecayProductFinalState", &genMu_isDirectPromptTauDecayProductFinalState);
        chain->SetBranchAddress("genMu_isHardProcess", &genMu_isHardProcess);
        chain->SetBranchAddress("genMu_isLastCopy", &genMu_isLastCopy);
        chain->SetBranchAddress("genMu_isLastCopyBeforeFSR", &genMu_isLastCopyBeforeFSR);
        chain->SetBranchAddress("genMu_isPromptDecayed", &genMu_isPromptDecayed);
        chain->SetBranchAddress("genMu_isDecayedLeptonHadron", &genMu_isDecayedLeptonHadron);
        chain->SetBranchAddress("genMu_fromHardProcessBeforeFSR", &genMu_fromHardProcessBeforeFSR);
        chain->SetBranchAddress("genMu_fromHardProcessDecayed", &genMu_fromHardProcessDecayed);
        chain->SetBranchAddress("genMu_fromHardProcessFinalState", &genMu_fromHardProcessFinalState);

        Int_t NEvents = chain->GetEntries();
        if (DEBUG) NEvents = 1000;
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassMu=0, nFailMu=0;
        UInt_t nSameSign=0;
        UInt_t nMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            for (UInt_t mu=0; mu<relPFiso->size(); mu++)
            {
                if (relPFiso->at(mu) < 0.15) nPassMu++;
                else nFailMu++;
            }
            if (!DEBUG) bar.Draw(i);

            // pre-selection
            if (p_T->size() != 2) continue;
//            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 52 || p_T->at(1) < 52) continue;
            if (fabs(eta->at(0)) >= 2.4 || fabs(eta->at(1)) >= 2.4) continue;

            if (charge->at(0) == charge->at(1)) {nSameSign++; continue;}
            nMu += 2;

            if (DEBUG)
            {
                if (nPass >= 10) break;
                cout << "Evt " << i << endl;
            }

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            Int_t i_lead = (p_T->at(0) >= p_T->at(1)) ? 0 : 1;
            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();

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

            std::vector<int> matches_lead;
            std::vector<int> matches_sub;
            Double_t dRmin_lead = 9999;
            Double_t dRmin_sub = 9999;
            Double_t rel_dpTmin_lead = 9999;
            Double_t rel_dpTmin_sub = 9999;
            Int_t i_bestmatch_lead = -1;
            Int_t i_bestmatch_sub = -1;
            Int_t veto_lead = 0; // used to find fake leptons
            Int_t veto_sub = 0; // used to find fake leptons
            if (DEBUG) cout << "  nGenMuons = " << genMu_pT->size() << endl;
            for (UInt_t i_gen=0; i_gen<genMu_pT->size(); i_gen++)
            {
                if (DEBUG)
                {
                    cout << "   GenMuon" << i_gen << ":  p_T=" << genMu_pT->at(i_gen);
                    cout << "  eta=" << genMu_eta->at(i_gen) << "  phi=" << genMu_phi->at(i_gen) << endl;
                    cout << "    isPrompt=" << genMu_isPrompt->at(i_gen) << "  isLastCopy=" << genMu_isLastCopy->at(i_gen);
                    cout << "  iLCbeforeFSR=" << genMu_isLastCopyBeforeFSR->at(i_gen) << endl;
                }

                Double_t dEta_lead = eta->at(i_lead) - genMu_eta->at(i_gen);
                Double_t dPhi_lead = phi->at(i_lead) - genMu_phi->at(i_gen);
                Double_t dR_lead = sqrt(dEta_lead*dEta_lead + dPhi_lead*dPhi_lead);
                Double_t rel_dpT_lead = fabs(p_T->at(i_lead) - genMu_pT->at(i_gen)) / p_T->at(i_lead);

                Double_t dEta_sub = eta->at(1-i_lead) - genMu_eta->at(i_gen);
                Double_t dPhi_sub = phi->at(1-i_lead) - genMu_phi->at(i_gen);
                Double_t dR_sub = sqrt(dEta_sub*dEta_sub + dPhi_sub*dPhi_sub);
                Double_t rel_dpT_sub = fabs(p_T->at(1-i_lead) - genMu_pT->at(i_gen)) / p_T->at(1-i_lead);
                if (dR_lead < 0.3)
                {
                    if (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen) || genMu_isTauDecayProduct->at(i_gen) ||
                        genMu_isPromptTauDecayProduct->at(i_gen) || genMu_isDirectPromptTauDecayProductFinalState->at(i_gen) ||
                        genMu_isHardProcess->at(i_gen) || genMu_fromHardProcessBeforeFSR->at(i_gen) || genMu_fromHardProcessFinalState->at(i_gen))
                        veto_lead = 1;
                    if (rel_dpT_lead < 0.5)
                    {
                        matches_lead.push_back(i_gen);
                        if (genMu_isLastCopy->at(i_gen) && (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen)))
                        {
                            if (dR_lead < dRmin_lead || (dR_lead == dRmin_lead && rel_dpT_lead <= rel_dpTmin_lead))
                            {
                                dRmin_lead = dR_lead;
                                rel_dpTmin_lead = rel_dpT_lead;
                                i_bestmatch_lead = i_gen;
                            }
                        }
                    }
                }
                if (dR_sub < 0.3)
                {
                    if (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen) || genMu_isTauDecayProduct->at(i_gen) ||
                        genMu_isPromptTauDecayProduct->at(i_gen) || genMu_isDirectPromptTauDecayProductFinalState->at(i_gen) ||
                        genMu_isHardProcess->at(i_gen) || genMu_fromHardProcessBeforeFSR->at(i_gen) || genMu_fromHardProcessFinalState->at(i_gen))
                        veto_sub = 1;
                    if (rel_dpT_sub < 0.5)
                    {
                        matches_sub.push_back(i_gen);
                        if (genMu_isLastCopy->at(i_gen) && (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen)))
                        {
                            if (dR_sub < dRmin_sub || (dR_sub == dRmin_sub && rel_dpT_sub <= rel_dpTmin_sub))
                            {
                                dRmin_sub = dR_sub;
                                rel_dpTmin_sub = rel_dpT_sub;
                                i_bestmatch_sub = i_gen;
                            }
                        }
                    }
                }
            }

            if (DEBUG)
            {
                cout << "\n  nMuons = " << p_T->size() << endl;
                cout << "   p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\trelPFiso[0] = " << relPFiso->at(0);
                Int_t printval = (i_lead == 0)?i_bestmatch_lead:i_bestmatch_sub;
                cout << "\tMatched to " << printval;
                if ((i_lead == 0 && veto_lead) || (i_lead != 0 && veto_sub))
                    cout << "\tVETOED";
                cout << endl;
                cout << "   p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\trelPFiso[1] = " << relPFiso->at(1);
                printval = (i_lead == 1)?i_bestmatch_lead:i_bestmatch_sub;
                cout << "\tMatched to " << printval;
                if ((i_lead == 1 && veto_lead) || (i_lead != 1 && veto_sub))
                    cout << "\tVETOED";
                cout << endl;
            }

            if ((i_bestmatch_lead < 0 && veto_lead) || (i_bestmatch_sub < 0 && veto_sub)) continue;
            nPass++;

            if (i_bestmatch_lead >= 0 && i_bestmatch_sub >= 0)
                h_mass_PP_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if (i_bestmatch_lead >= 0 && (i_bestmatch_sub < 0 && !veto_sub))
                h_mass_PF_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if ((i_bestmatch_lead < 0 && !veto_lead) && i_bestmatch_sub >= 0)
                h_mass_FP_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            if ((i_bestmatch_lead < 0 && !veto_lead) && (i_bestmatch_sub < 0 && !veto_sub))
                h_mass_FF_true[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            // -- FR WEIGHTS -- //
            std::vector<double> MMweights(4, 1.0);
            Double_t FR1=0, FR2=0, PR1=1, PR2=1;

            MMweights = analyzer->MatrixMethod_MC(p_T->at(i_lead), eta->at(i_lead), relPFiso->at(i_lead),
                                                  p_T->at(1-i_lead), eta->at(1-i_lead), relPFiso->at(1-i_lead), "ttbar", "ttbar");
            FR1 = analyzer->FakeRate_MC(p_T->at(i_lead), eta->at(i_lead), "QCD");
            FR2 = analyzer->FakeRate_MC(p_T->at(1-i_lead), eta->at(1-i_lead), "QCD");
            PR1 = analyzer->PromptRate_MC(p_T->at(i_lead), eta->at(i_lead), "DY");
            PR2 = analyzer->PromptRate_MC(p_T->at(1-i_lead), eta->at(1-i_lead), "DY");
//            MMweights = analyzer->MatrixMethod_eta_MC(eta->at(i_lead), relPFiso->at(i_lead),
//                                                      eta->at(1-i_lead), relPFiso->at(1-i_lead), "ttbar", "ttbar");
//            FR1 = analyzer->FakeRate_eta_MC(eta->at(i_lead), "QCD");
//            FR2 = analyzer->FakeRate_eta_MC(eta->at(1-i_lead), "QCD");
//            PR1 = analyzer->PromptRate_eta_MC(eta->at(i_lead), "DY");
//            PR2 = analyzer->PromptRate_eta_MC(eta->at(1-i_lead), "DY");

            // Filling histos
            h_mass_PP[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[0]);
            h_mass_PF[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[1]);
            h_mass_FP[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[2]);
            h_mass_FF[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MMweights[3]);

            if(relPFiso->at(i_lead) < 0.15 && relPFiso->at(1-i_lead) < 0.15)
                h_mass_TT[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(relPFiso->at(i_lead) < 0.15 && relPFiso->at(1-i_lead) >= 0.15)
                h_mass_TN[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(relPFiso->at(i_lead) >= 0.15 && relPFiso->at(1-i_lead) < 0.15)
                h_mass_NT[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            else if(relPFiso->at(i_lead) >= 0.15 && relPFiso->at(1-i_lead) >= 0.15)
                h_mass_NN[i_type]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

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

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _ttbar_1000toInf) pr = _EndOf_VVnST_Normal; // next -- WJets

        cout << "nSameSign = " << nSameSign << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    for (Int_t i=0; i<4; i++)
    {
        h_mass_PP[i]->Write();
        h_mass_PF[i]->Write();
        h_mass_FP[i]->Write();
        h_mass_FF[i]->Write();
        h_mass_TT[i]->Write();
        h_mass_TN[i]->Write();
        h_mass_NT[i]->Write();
        h_mass_NN[i]->Write();
        h_mass_PP_true[i]->Write();
        h_mass_PF_true[i]->Write();
        h_mass_FP_true[i]->Write();
        h_mass_FF_true[i]->Write();
    }

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_mu_MC()


void MatrixMethod_singleMu_MC (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_singleMu_MC"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRandPRvalues_MC();

    // Histogram containers
    TH1D* h_pT_P[4];
    TH1D* h_pT_F[4];
    TH1D* h_pT_T[4];
    TH1D* h_pT_N[4];
    TH1D* h_pT_P_true[4];
    TH1D* h_pT_F_true[4];

    TString names[4] = {"DY", "ttbar", "WJets", "QCD"};
    // -- Creating Histograms -- //
    for (Int_t i_type=0; i_type<4; i_type++)
    {
        h_pT_P[i_type] = new TH1D("h_pT_P_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
        h_pT_F[i_type] = new TH1D("h_pT_F_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
        h_pT_T[i_type] = new TH1D("h_pT_T_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
        h_pT_N[i_type] = new TH1D("h_pT_N_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
        h_pT_P_true[i_type] = new TH1D("h_pT_P_true_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
        h_pT_F_true[i_type] = new TH1D("h_pT_F_true_"+names[i_type], "", nPtBinEndcap, analyzer->ptbin_endcap);
    }
    Int_t i_type = 0;

    Process_t pr_first = _DY_10to50, pr_last = _EndOf_QCDMuEnriched_Normal;
    if (DEBUG) pr_first = _ttbar, pr_last = _ttbar_700to1000;
    for (Process_t pr=pr_first; pr<pr_last; pr=next(pr))
    {
        Mgr.SetProc(pr);

        if (pr < _EndOf_DY_Normal) i_type = 0;
        else if (pr < _EndOf_ttbar_Normal) i_type = 1;
        else if (pr < _EndOf_WJets_Normal) i_type = 2;
        else if (pr < _EndOf_QCDMuEnriched_Normal) i_type = 3;

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

        Int_t nMuons;
        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nGenMuons;
        std::vector<double> *genMu_pT = new std::vector<double>;
        std::vector<double> *genMu_eta = new std::vector<double>;
        std::vector<double> *genMu_phi = new std::vector<double>;
        std::vector<int> *genMu_charge = new std::vector<int>;
        std::vector<int> *genMu_isPrompt = new std::vector<int>;
        std::vector<int> *genMu_isPromptFinalState = new std::vector<int>;
        std::vector<int> *genMu_isTauDecayProduct = new std::vector<int>;
        std::vector<int> *genMu_isPromptTauDecayProduct = new std::vector<int>;
        std::vector<int> *genMu_isDirectPromptTauDecayProductFinalState = new std::vector<int>;
        std::vector<int> *genMu_isHardProcess = new std::vector<int>;
        std::vector<int> *genMu_isLastCopy = new std::vector<int>;
        std::vector<int> *genMu_isLastCopyBeforeFSR = new std::vector<int>;
        std::vector<int> *genMu_isPromptDecayed = new std::vector<int>;
        std::vector<int> *genMu_isDecayedLeptonHadron = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessBeforeFSR = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessDecayed = new std::vector<int>;
        std::vector<int> *genMu_fromHardProcessFinalState = new std::vector<int>;
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

        chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("nMuons", 1);
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
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchStatus("nGenMuons", 1);
        chain->SetBranchStatus("genMu_pT", 1);
        chain->SetBranchStatus("genMu_eta", 1);
        chain->SetBranchStatus("genMu_phi", 1);
        chain->SetBranchStatus("genMu_charge", 1);
        chain->SetBranchStatus("genMu_isPrompt", 1);
        chain->SetBranchStatus("genMu_isPromptFinalState", 1);
        chain->SetBranchStatus("genMu_isTauDecayProduct", 1);
        chain->SetBranchStatus("genMu_isPromptTauDecayProduct", 1);
        chain->SetBranchStatus("genMu_isDirectPromptTauDecayProductFinalState", 1);
        chain->SetBranchStatus("genMu_isHardProcess", 1);
        chain->SetBranchStatus("genMu_isLastCopy", 1);
        chain->SetBranchStatus("genMu_isLastCopyBeforeFSR", 1);
        chain->SetBranchStatus("genMu_isPromptDecayed", 1);
        chain->SetBranchStatus("genMu_isDecayedLeptonHadron", 1);
        chain->SetBranchStatus("genMu_fromHardProcessBeforeFSR", 1);
        chain->SetBranchStatus("genMu_fromHardProcessDecayed", 1);
        chain->SetBranchStatus("genMu_fromHardProcessFinalState", 1);

        chain->SetBranchAddress("nMuons", &nMuons);
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
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);
        chain->SetBranchAddress("nGenMuons", &nGenMuons);
        chain->SetBranchAddress("genMu_pT", &genMu_pT);
        chain->SetBranchAddress("genMu_eta", &genMu_eta);
        chain->SetBranchAddress("genMu_phi", &genMu_phi);
        chain->SetBranchAddress("genMu_charge", &genMu_charge);
        chain->SetBranchAddress("genMu_isPrompt", &genMu_isPrompt);
        chain->SetBranchAddress("genMu_isPromptFinalState", &genMu_isPromptFinalState);
        chain->SetBranchAddress("genMu_isTauDecayProduct", &genMu_isTauDecayProduct);
        chain->SetBranchAddress("genMu_isPromptTauDecayProduct", &genMu_isPromptTauDecayProduct);
        chain->SetBranchAddress("genMu_isDirectPromptTauDecayProductFinalState", &genMu_isDirectPromptTauDecayProductFinalState);
        chain->SetBranchAddress("genMu_isHardProcess", &genMu_isHardProcess);
        chain->SetBranchAddress("genMu_isLastCopy", &genMu_isLastCopy);
        chain->SetBranchAddress("genMu_isLastCopyBeforeFSR", &genMu_isLastCopyBeforeFSR);
        chain->SetBranchAddress("genMu_isPromptDecayed", &genMu_isPromptDecayed);
        chain->SetBranchAddress("genMu_isDecayedLeptonHadron", &genMu_isDecayedLeptonHadron);
        chain->SetBranchAddress("genMu_fromHardProcessBeforeFSR", &genMu_fromHardProcessBeforeFSR);
        chain->SetBranchAddress("genMu_fromHardProcessDecayed", &genMu_fromHardProcessDecayed);
        chain->SetBranchAddress("genMu_fromHardProcessFinalState", &genMu_fromHardProcessFinalState);

        Int_t NEvents = chain->GetEntries();
        if (DEBUG) NEvents = 1000;
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassMu=0, nFailMu=0;
        UInt_t nSameSign=0;
        UInt_t nMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);
            else
            {
                cout << "Evt " << i << endl;
                cout << "  mMuons=" << p_T->size() << endl;
            }

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

            for (UInt_t i_mu=0; i_mu<relPFiso->size(); i_mu++)
            {
                if (relPFiso->at(i_mu) < 0.15) nPassMu++;
                else nFailMu++;

                //Preselection
                if (p_T->at(i_mu) < 52 || fabs(eta->at(i_mu)) >= 2.4) continue;
                if (p_T->at(i_mu) != p_T->at(i_mu)) cout << p_T->at(i_mu) << " " << eta->at(i_mu) << " " << phi->at(i_mu) <<
                                                            " " << charge->at(i_mu) << " " << relPFiso->at(i_mu) << endl;

                std::vector<int> matches;
                Double_t dRmin = 9999;
                Double_t rel_dpTmin = 9999;
                Int_t i_bestmatch = -1;
                Int_t veto = 0; // used to find fake leptons
                for (UInt_t i_gen=0; i_gen<genMu_pT->size(); i_gen++)
                {
                    Double_t dEta = eta->at(i_mu) - genMu_eta->at(i_gen);
                    Double_t dPhi = phi->at(i_mu) - genMu_phi->at(i_gen);
                    Double_t dR = sqrt(dEta*dEta + dPhi*dPhi);
                    Double_t rel_dpT = fabs(p_T->at(i_mu) - genMu_pT->at(i_gen)) / p_T->at(i_mu);
                    if (dR < 0.3)
                    {
                        if (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen) || genMu_isTauDecayProduct->at(i_gen) ||
                            genMu_isPromptTauDecayProduct->at(i_gen) || genMu_isDirectPromptTauDecayProductFinalState->at(i_gen) ||
                            genMu_isHardProcess->at(i_gen) || genMu_fromHardProcessBeforeFSR->at(i_gen) || genMu_fromHardProcessFinalState->at(i_gen))
                            veto = 1;
                        if (rel_dpT < 0.5)
                        {
                            matches.push_back(i_gen);
                            if (genMu_isLastCopy->at(i_gen) && (genMu_isPrompt->at(i_gen) || genMu_isPromptFinalState->at(i_gen)))
                            {
                                if (dR < dRmin || (dR == dRmin && rel_dpT <= rel_dpTmin))
                                {
                                    dRmin = dR;
                                    rel_dpTmin = rel_dpT;
                                    i_bestmatch = i_gen;
                                }
                            }
                        }
                    }
                }

                Double_t FR = analyzer->FakeRate_MC(p_T->at(i_mu), eta->at(i_mu), "ttbar");
                Double_t PR = analyzer->PromptRate_MC(p_T->at(i_mu), eta->at(i_mu), "ttbar");

                if (i_bestmatch >= 0)
                {
                    if (DEBUG)
                    {
                        cout << "  isGenMatched = 1\n  Matched with " << matches.size() <<" gen leptons. Best match = " << i_bestmatch << endl;
                        for (UInt_t i_match=0; i_match<matches.size(); i_match++)
                        {
                            cout << "    match" << i_match << ":  pT=" << genMu_pT->at(matches[i_match]);
                            cout << "  eta=" << genMu_eta->at(matches[i_match]) << "  phi=" << genMu_phi->at(matches[i_match]);
                            cout << "  charge=" << genMu_charge->at(matches[i_match]) << endl;
                            cout << "     isPrompt = " << genMu_isPrompt->at(matches[i_match]) << endl;
                            cout << "     isPromptFinalState = " << genMu_isPromptFinalState->at(matches[i_match]) << endl;
                            cout << "     isTauDecayProduct = " << genMu_isTauDecayProduct->at(matches[i_match]) << endl;
                            cout << "     isPromptTauDecayProduct = " << genMu_isPromptTauDecayProduct->at(matches[i_match]) << endl;
                            cout << "     isDirectPromptTauDecayProductFinalState = " << genMu_isDirectPromptTauDecayProductFinalState->at(matches[i_match]) << endl;
                            cout << "     isHardProcess = " << genMu_isHardProcess->at(matches[i_match]) << endl;
                            cout << "     isLastCopy = " << genMu_isLastCopy->at(matches[i_match]) << endl;
                            cout << "     isLastCopyBeforeFSR = " << genMu_isLastCopyBeforeFSR->at(matches[i_match]) << endl;
                            cout << "     isPromptDecayed = " << genMu_isPromptDecayed->at(matches[i_match]) << endl;
                            cout << "     isDecayedLeptonHadron = " << genMu_isDecayedLeptonHadron->at(matches[i_match]) << endl;
                            cout << "     fromHardProcessBeforeFSR = " << genMu_fromHardProcessBeforeFSR->at(matches[i_match]) << endl;
                            cout << "     fromHardProcessDecayed = " << genMu_fromHardProcessDecayed->at(matches[i_match]) << endl;
                            cout << "     fromHardProcessFinalState = " << genMu_fromHardProcessFinalState->at(matches[i_match]) << endl;
                        }
                    }
                    h_pT_P_true[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                }
                else if (!veto)
                {
                    h_pT_F_true[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (DEBUG) cout << "  isGenMatched = 0" << endl;
                }
                else
                {
                    if (DEBUG) cout << " isGenMatched = -1 (not matched but also vetoed as fake)" << endl;
                    continue;
                }

                if (relPFiso->at(i_mu) < 0.15)
                {
                    Double_t MM_weight1 = (1-FR) / (PR - FR); // contribution to prompt events
                    Double_t MM_weight2 = (PR-1) / (PR - FR); // contribution to fake events
                    h_pT_P[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MM_weight1);
                    h_pT_F[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MM_weight2);
                    h_pT_T[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else
                {
                    Double_t MM_weight1 = -FR / (PR - FR); // contribution to prompt events
                    Double_t MM_weight2 =  PR / (PR - FR); // contribution to fake events
                    h_pT_P[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MM_weight1);
                    h_pT_F[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MM_weight2);
                    h_pT_N[i_type]->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }

            }

        }// End of event iteration

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

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _ttbar_1000toInf) pr = _EndOf_VVnST_Normal; // next -- WJets

        cout << "nSameSign = " << nSameSign << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    for (Int_t i=0; i<4; i++)
    {
        h_pT_P[i]->Write();
        h_pT_F[i]->Write();
        h_pT_T[i]->Write();
        h_pT_N[i]->Write();
        h_pT_P_true[i]->Write();
        h_pT_F_true[i]->Write();
    }

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_singleMu_MC()


/// -------------------------------- EMu Channel ------------------------------------ ///
void MatrixMethod_EMu (Bool_t DEBUG=kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/EMu/";
    TString Dir_e = "/media/sf_DATA/FR/Electron/";
    TString Dir_mu = "/media/sf_DATA/FR/Muon/";
    TString Dir1 = "/media/sf_DownloadData/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    TString outName = Dir+"MatrixMethod_EMu"+debug+".root";
    f = new TFile(outName, "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();
    analyzer->SetupEfficiencyScaleFactor_electron();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir_mu+"FakeRate_muon.root");
    analyzer->SetupPRvalues(Dir_mu+"PromptRate_muon_alt.root", "ratio", 1);
    analyzer->SetupFRvalues_ele(Dir_e+"FakeRate_electron.root");
    analyzer->SetupPRvalues_ele(Dir_e+"PromptRate_electron_alt.root");

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
        TH1D* h_mass_PP = new TH1D("h_mass_PP_"+Mgr.Procname[pr], "h_mass_PP_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_PP->Sumw2();
        TH1D* h_mass_PF = new TH1D("h_mass_PF_"+Mgr.Procname[pr], "h_mass_PF_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_PF->Sumw2();
        TH1D* h_mass_FP = new TH1D("h_mass_FP_"+Mgr.Procname[pr], "h_mass_FP_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_FP->Sumw2();
        TH1D* h_mass_FF = new TH1D("h_mass_FF_"+Mgr.Procname[pr], "h_mass_FF_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_FF->Sumw2();
        TH1D* h_mass_PP_TT = new TH1D("h_mass_PP_TT_"+Mgr.Procname[pr], "h_mass_PP_TT_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_PP_TT->Sumw2();
        TH1D* h_mass_PF_TT = new TH1D("h_mass_PF_TT_"+Mgr.Procname[pr], "h_mass_PF_TT_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_PF_TT->Sumw2();
        TH1D* h_mass_FP_TT = new TH1D("h_mass_FP_TT_"+Mgr.Procname[pr], "h_mass_FP_TT_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_FP_TT->Sumw2();
        TH1D* h_mass_FF_TT = new TH1D("h_mass_FF_TT_"+Mgr.Procname[pr], "h_mass_FF_TT_"+Mgr.Procname[pr], 48,-2.4,2.4/*binnum, massbins*/); h_mass_FF_TT->Sumw2();

        TH1D* h_massSS_PP = new TH1D("h_massSS_PP_"+Mgr.Procname[pr], "h_massSS_PP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP->Sumw2();
        TH1D* h_massSS_PF = new TH1D("h_massSS_PF_"+Mgr.Procname[pr], "h_massSS_PF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF->Sumw2();
        TH1D* h_massSS_FP = new TH1D("h_massSS_FP_"+Mgr.Procname[pr], "h_massSS_FP_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP->Sumw2();
        TH1D* h_massSS_FF = new TH1D("h_massSS_FF_"+Mgr.Procname[pr], "h_massSS_FF_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF->Sumw2();
        TH1D* h_massSS_PP_TT = new TH1D("h_massSS_PP_TT_"+Mgr.Procname[pr], "h_massSS_PP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PP_TT->Sumw2();
        TH1D* h_massSS_PF_TT = new TH1D("h_massSS_PF_TT_"+Mgr.Procname[pr], "h_massSS_PF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_PF_TT->Sumw2();
        TH1D* h_massSS_FP_TT = new TH1D("h_massSS_FP_TT_"+Mgr.Procname[pr], "h_massSS_FP_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FP_TT->Sumw2();
        TH1D* h_massSS_FF_TT = new TH1D("h_massSS_FF_TT_"+Mgr.Procname[pr], "h_massSS_FF_TT_"+Mgr.Procname[pr], binnum, massbins); h_massSS_FF_TT->Sumw2();

        Double_t e_p_T, e_eta, e_etaSC, e_phi;
        Int_t e_charge, e_scPixCharge, e_isGsfCtfScPixConsistent, e_isGsfScPixConsistent, e_isGsfCtfConsistent;
        Double_t e_Full5x5_SigmaIEtaIEta, e_dEtaInSeed, e_dPhiIn, e_HoverE, e_InvEminusInvP;
        Double_t e_relPFiso;
        Double_t e_relECALiso, e_relHCALiso, e_relTrkIso;
        Int_t e_mHits, e_passConvVeto, e_passMediumID;
        Double_t mu_p_T, mu_eta, mu_phi;
        Int_t mu_charge;
        Double_t mu_relPFiso;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU, nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
        Int_t runNum;
        Int_t lumiBlock;


        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("e_p_T", 1);
        chain->SetBranchStatus("e_eta", 1);
        chain->SetBranchStatus("e_etaSC", 1);
        chain->SetBranchStatus("e_phi", 1);
        chain->SetBranchStatus("e_charge", 1);
        chain->SetBranchStatus("e_scPixCharge", 1);
        chain->SetBranchStatus("e_isGsfCtfScPixConsistent", 1);
        chain->SetBranchStatus("e_isGsfScPixConsistent", 1);
        chain->SetBranchStatus("e_isGsfCtfConsistent", 1);
        chain->SetBranchStatus("e_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("e_dEtaInSeed", 1);
        chain->SetBranchStatus("e_dPhiIn", 1);
        chain->SetBranchStatus("e_HoverE", 1);
        chain->SetBranchStatus("e_InvEminusInvP", 1);
        chain->SetBranchStatus("e_relECALiso", 1);
        chain->SetBranchStatus("e_relHCALiso", 1);
        chain->SetBranchStatus("e_relTrkIso", 1);
        chain->SetBranchStatus("e_mHits", 1);
        chain->SetBranchStatus("e_passConvVeto", 1);
        chain->SetBranchStatus("e_relPFiso_Rho", 1);
        chain->SetBranchStatus("e_passMediumID", 1);
        chain->SetBranchStatus("mu_p_T", 1);
        chain->SetBranchStatus("mu_eta", 1);
        chain->SetBranchStatus("mu_phi", 1);
        chain->SetBranchStatus("mu_charge", 1);
        chain->SetBranchStatus("mu_relPFiso_dBeta", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("prescale_factor", 1);
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
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);

        chain->SetBranchAddress("e_p_T", &e_p_T);
        chain->SetBranchAddress("e_eta", &e_eta);
        chain->SetBranchAddress("e_etaSC", &e_etaSC);
        chain->SetBranchAddress("e_phi", &e_phi);
        chain->SetBranchAddress("e_charge", &e_charge);
        chain->SetBranchAddress("e_scPixCharge", &e_scPixCharge);
        chain->SetBranchAddress("e_isGsfCtfScPixConsistent", &e_isGsfCtfScPixConsistent);
        chain->SetBranchAddress("e_isGsfScPixConsistent", &e_isGsfScPixConsistent);
        chain->SetBranchAddress("e_isGsfCtfConsistent", &e_isGsfCtfConsistent);
        chain->SetBranchAddress("e_Full5x5_SigmaIEtaIEta", &e_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("e_dEtaInSeed", &e_dEtaInSeed);
        chain->SetBranchAddress("e_dPhiIn", &e_dPhiIn);
        chain->SetBranchAddress("e_HoverE", &e_HoverE);
        chain->SetBranchAddress("e_InvEminusInvP", &e_InvEminusInvP);
        chain->SetBranchAddress("e_relECALiso", &e_relECALiso);
        chain->SetBranchAddress("e_relHCALiso", &e_relHCALiso);
        chain->SetBranchAddress("e_relTrkIso", &e_relTrkIso);
        chain->SetBranchAddress("e_mHits", &e_mHits);
        chain->SetBranchAddress("e_passConvVeto", &e_passConvVeto);
        chain->SetBranchAddress("e_relPFiso_Rho", &e_relPFiso);
        chain->SetBranchAddress("e_passMediumID", &e_passMediumID);
        chain->SetBranchAddress("mu_p_T", &mu_p_T);
        chain->SetBranchAddress("mu_eta", &mu_eta);
        chain->SetBranchAddress("mu_phi", &mu_phi);
        chain->SetBranchAddress("mu_charge", &mu_charge);
        chain->SetBranchAddress("mu_relPFiso_dBeta", &mu_relPFiso);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
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
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        myProgressBar_t bar(NEvents);
        Int_t nSameSign = 0;

        for(Int_t i=0; i<NEvents; i++)
        {
//            if (pr == _ttbar) continue;
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // pre-selection
            if (mu_p_T < 52 || e_p_T < 25) continue;
            if (fabs(mu_eta) >= 2.4 || fabs(e_eta) >= 2.4 || (fabs(e_eta) >= 1.4442 && fabs(e_eta) <= 1.566)) continue;
            if (fabs(e_eta) < 1.4442 && (e_Full5x5_SigmaIEtaIEta >= 0.013 || e_HoverE >= 0.13 ||
                                         fabs(e_dEtaInSeed) >= 0.01 || fabs(e_dPhiIn) >= 0.07)) continue;
            else if (fabs(e_eta) > 1.566 && (e_Full5x5_SigmaIEtaIEta >= 0.035 || e_HoverE >= 0.13)) continue;
            if (e_mHits > 1 /*|| e_relECALiso >= 0.5 || e_relHCALiso >= 0.3 || e_relTrkIso >= 0.2*/) continue;

            if (mu_charge == e_charge) nSameSign++;

            if (mu_p_T != mu_p_T) cout << mu_p_T << " " << mu_eta << " " << mu_phi << " " << mu_charge << " " << mu_relPFiso << endl;
            if (e_p_T != e_p_T) cout << e_p_T << " " << e_eta << " " << e_phi << " " << e_charge << " " << e_passMediumID << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "p_T[mu] = " << mu_p_T;
                cout << "\teta[mu] = " << mu_eta;
                cout << "\tphi[mu] = " << mu_phi << endl;
                cout << "\trelPFiso[mu] = " << mu_relPFiso << endl;
                cout << "p_T[e] = " << e_p_T;
                cout << "\teta[e] = " << e_eta;
                cout << "\tphi[e] = " << e_phi << endl;
                cout << "\tpassMediumID[e] = " << e_passMediumID << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu, ele;
            mu.SetPtEtaPhiM(mu_p_T, mu_eta, mu_phi, M_Mu);
            ele.SetPtEtaPhiM(e_p_T, e_eta, e_phi, M_Elec);
            Double_t mass = (mu+ele).M();
            Double_t pTll = (mu+ele).Pt();
            Double_t etall = (mu+ele).Eta();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;
            // -- Efficiency scale factor -- //
//            if(Mgr.isMC == kTRUE)
//            {
//                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(MuMu);
//                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(MuMu);
//                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
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
            std::vector<double> MatrixMethod_weights(4, 1.0);
            Double_t FRmu=0, FRe=0, PRmu=1, PRe=1;
            MatrixMethod_weights = analyzer->MatrixMethod_evtWeight_EMu(mu_p_T, mu_eta, mu_relPFiso, e_p_T, e_eta, e_passMediumID);
            FRmu = analyzer->FakeRate(mu_p_T, mu_eta);
            FRe  = analyzer->FakeRate_ele(e_p_T, e_eta);
            PRmu = analyzer->PromptRate(mu_p_T, mu_eta, 1);
            PRe  = analyzer->PromptRate_ele(e_p_T, e_eta);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            if (mu_charge != e_charge)
            {
                h_mass_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
                h_mass_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
                h_mass_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
                h_mass_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
                h_mass_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PRmu * PRe);
                h_mass_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PRmu * FRe);
                h_mass_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FRmu * PRe);
                h_mass_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FRmu * FRe);
            }
            else
            {
                h_massSS_PP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0]);
                h_massSS_PF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1]);
                h_massSS_FP->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2]);
                h_massSS_FF->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3]);
                h_massSS_PP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[0] * PRmu * PRe);
                h_massSS_PF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[1] * PRmu * FRe);
                h_massSS_FP_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[2] * FRmu * PRe);
                h_massSS_FF_TT->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * MatrixMethod_weights[3] * FRmu * FRe);
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

        h_mass_PP->Write();
        h_mass_PF->Write();
        h_mass_FP->Write();
        h_mass_FF->Write();
        h_mass_PP_TT->Write();
        h_mass_PF_TT->Write();
        h_mass_FP_TT->Write();
        h_mass_FF_TT->Write();
        h_massSS_PP->Write();
        h_massSS_PF->Write();
        h_massSS_FP->Write();
        h_massSS_FF->Write();
        h_massSS_PP_TT->Write();
        h_massSS_PF_TT->Write();
        h_massSS_FP_TT->Write();
        h_massSS_FF_TT->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "nSameSign = " << nSameSign << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    f->Close();
    if (!f->IsOpen()) cout << "File " << outName << " has been closed successfully.\n" << endl;
    else cout << "FILE " << outName << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MatrixMethod_EMu()
