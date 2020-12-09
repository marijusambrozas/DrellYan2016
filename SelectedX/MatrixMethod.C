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
void MatrixMethod_mu (Bool_t DEBUG);

void MatrixMethod (TString WhichX = "", Bool_t DEBUG=kFALSE)
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;

    if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*****  MatrixMethod_mu  *****" << endl;
        MatrixMethod_mu(DEBUG);
    }
    else if (whichX.Contains("E"))
    {
        Xselected++;
        cout << "\n*****    MatrixMethod_ele    *****" << endl;
        MatrixMethod_ele(DEBUG);
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
//    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", "subtract");
    analyzer->SetupFRvalues_ele_fit(Dir+"FakeRate_electron.root");
//    analyzer->SetupPRvalues_ele(Dir+"PromptRate_electron.root", "subtract");
    analyzer->SetupPRvalues_ele_fit(Dir+"PromptRate_electron.root");

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

            if (mHits->at(0) > 1 /*|| !passConvVeto->at(0)*/ || relECALiso->at(0) >= 0.5 || relHCALiso->at(0) >= 0.5 || relTrkIso->at(0) >= 0.2) continue;
            if (mHits->at(1) > 1 /*|| !passConvVeto->at(1)*/ || relECALiso->at(1) >= 0.5 || relHCALiso->at(1) >= 0.5 || relTrkIso->at(1) >= 0.2) continue;

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
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight_ele_fit(p_T->at(0), etaSC->at(0), passMediumID->at(0), p_T->at(1), etaSC->at(1), passMediumID->at(1));
                FR1 = analyzer->FakeRate_ele_fit(p_T->at(0), etaSC->at(0));
                FR2 = analyzer->FakeRate_ele_fit(p_T->at(1), etaSC->at(1));
                PR1 = analyzer->PromptRate_ele_fit(p_T->at(0), etaSC->at(0));
                PR2 = analyzer->PromptRate_ele_fit(p_T->at(1), etaSC->at(1));
            }
            else
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight_ele_fit(p_T->at(1), etaSC->at(1), passMediumID->at(1), p_T->at(0), etaSC->at(0), passMediumID->at(0));
                FR1 = analyzer->FakeRate_ele_fit(p_T->at(1), etaSC->at(1));
                FR2 = analyzer->FakeRate_ele_fit(p_T->at(0), etaSC->at(0));
                PR1 = analyzer->PromptRate_ele_fit(p_T->at(1), etaSC->at(1));
                PR2 = analyzer->PromptRate_ele_fit(p_T->at(0), etaSC->at(0));
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


/// -------------------------------- Electron Channel ------------------------------------ ///
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
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root");
    analyzer->SetupPRvalues(Dir+"PromptRate_muon.root");

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
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 52 && p_T->at(1) < 52) continue;
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
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (mu1+mu2).M();

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
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight(p_T->at(0), eta->at(0), relPFiso->at(0), p_T->at(1), eta->at(1), relPFiso->at(1));
                FR1 = analyzer->FakeRate(p_T->at(0), eta->at(0));
                FR2 = analyzer->FakeRate(p_T->at(1), eta->at(1));
                PR1 = analyzer->PromptRate(p_T->at(0), eta->at(0));
                PR2 = analyzer->PromptRate(p_T->at(1), eta->at(1));
            }
            else
            {
                MatrixMethod_weights = analyzer->MatrixMethod_evtWeight(p_T->at(1), eta->at(1), relPFiso->at(1), p_T->at(0), eta->at(0), relPFiso->at(0));
                FR1 = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FR2 = analyzer->FakeRate(p_T->at(0), eta->at(0));
                PR1 = analyzer->PromptRate(p_T->at(1), eta->at(1));
                PR2 = analyzer->PromptRate(p_T->at(0), eta->at(0));
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
