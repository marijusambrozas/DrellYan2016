#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <tuple>


// -- Macro for making new data files with only selection-passing events  -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.cc"
#include "./header/FileMgr.h"
#include "./etc/RoccoR/RoccoR.cc"

void MakeSelectionForFR_E (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForFR_E_alt (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForFR_Mu (TString type, TString HLTname, Bool_t Debug);
// PROMPT RATE
void MakeSelectionForPR_Mu (TString type, TString HLTname, Bool_t Debug);

void MakeSelectionForBKGest_EE (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForBKGest_MuMu (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForBKGest_EMu (TString type, TString HLTname, Bool_t Debug);
void CountObjectsInAcceptance (TString type, TString HLTname, Bool_t Debug);

Int_t skipLumiSection (Int_t runNo, Int_t lumiSection);
Int_t skipRun (Int_t runNo);

void MakeSelectionForFR (TString WhichX, TString type = "", TString HLTname = "DEFAULT")
{
    TString whichX = WhichX;
    whichX.ToUpper();
    TString HLT;
    Int_t Xselected = 0;
    Bool_t Debug = kFALSE;

    // for non-interactive status output
    ofstream fs;
    fs.open("./Status/StatusFR_"+whichX+"_"+type+".txt");
    fs << "Process MakeSelectionForFR (" << whichX << ", " << type << ", " << HLTname << ") initiated.\n";

    if (whichX.Contains("DEBUG"))
    {
        Debug = kTRUE;
        cout << "DEBUG MODE ON. Running on 100 events only" << endl;
    }
    if (whichX.Contains("COUNT"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "None";
        else HLT = HLTname;
        cout << "\n*******      CountObjectsInAcceptance (" << type << ", " << HLT << ")      *******" << endl;
        CountObjectsInAcceptance(type, HLTname, Debug);
    }
    else if (whichX.Contains("EMU") && whichX.Contains("EST"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "Mu_OR";
        else HLT = HLTname;
        cout << "\n*****  MakeSelectionForBKGest_EMu (" << type << ", " << HLT << ")  *****" << endl;
        MakeSelectionForBKGest_EMu(type, HLT, Debug);
    }
    else if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("EST"))
        {            
            if (HLTname == "DEFAULT") HLT = "Mu20_OR_Mu27_OR_Mu50";
            else HLT = HLTname;
            cout << "\n*****  MakeSelectionForBKGest_MuMu (" << type << ", " << HLT << ")  *****" << endl;
            MakeSelectionForBKGest_MuMu(type, HLT, Debug);
        }
        else if (whichX.Contains("PR"))
        {
            if (HLTname == "DEFAULT") HLT = "Mu50";
            else HLT = HLTname;
            cout << "\n*****  MakeSelectionForPR_Mu (" << type << ", " << HLT << ")  *****" << endl;
            MakeSelectionForPR_Mu(type, HLT, Debug);
        }
        else
        {
            if (HLTname == "DEFAULT") HLT = "Mu50";
            else HLT = HLTname;
            cout << "\n*****  MakeSelectionForFR_Mu (" << type << ", " << HLT << ")  *****" << endl;
            MakeSelectionForFR_Mu(type, HLT, Debug);
        }
    }
    else if (whichX.Contains("E"))
    {
        if (whichX.Contains("EST"))
        {
            Xselected++;
            if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
            else HLT = HLTname;
            cout << "\n*******      MakeSelectionForBKGest_EE (" << type << ", " << HLT << ")      *******" << endl;
            MakeSelectionForBKGest_EE(type, HLT, Debug);
        }
        else if (whichX.Contains("ALT"))
        {
            Xselected++;
            if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
            else HLT = HLTname;
            cout << "\n*******      MakeSelectionForFR_E_alt (" << type << ", " << HLT << ")      *******" << endl;
            MakeSelectionForFR_E_alt(type, HLT, Debug);
        }
        else
        {
            Xselected++;
            if (HLTname == "DEFAULT") HLT = "Photon_OR";
            else HLT = HLTname;
            cout << "\n*******      MakeSelectionForFR_E (" << type << ", " << HLT << ")      *******" << endl;
            MakeSelectionForFR_E(type, HLT, Debug);
        }
    }

    if (Xselected == 0) {
        cout << "Wrong arument!" << endl;
        fs << "Process MakeSelectionForFR (" << whichX << ", " << type << ", " << HLTname << ") FAILED. Wrong arguments!\n";
    }
    else fs << "Process MakeSelectionForFR (" << whichX << ", " << type << ", " << HLTname << ") FINISHED.\n";

    fs.close();

} // End of MakeSelectionForFR()

/////////////////////////////////////////////////////////////////////////////
/// --------------------------- FR estimation --------------------------- ///
/// /////////////////////////////////////////////////////////////////////////

/// ----------------------------- Electron ------------------------------ ///
void MakeSelectionForFR_E (TString type, TString HLTname , Bool_t Debug)
{
    gInterpreter->GenerateDictionary("std::vector<std::vector<std::pair<std::string,int>>>", "vector");

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;   

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* ElectronFile;
        out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
//        out_base = "~/Desktop/";
        out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];

        if (Mgr.Type == "TEST")
            out_base = "/media/sf_DATA/test/";

        if (Debug == kTRUE)
            ElectronFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            ElectronFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        ElectronFile->cd();

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
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *HLT_trigPS = new std::vector<int>;
        std::vector<int> *L1_trigPS = new std::vector<int>;
//        std::vector<std::vector<std::pair<std::string, int>>> L1seed_trigPSinDetail;
        Double_t MET_pT, MET_phi;
        Int_t runNum;
        Int_t lumiBlock;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        // -- Creating electron variables to assign branches -- //
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
        ElectronTree->Branch("scPixCharge", &scPixCharge);
        ElectronTree->Branch("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        ElectronTree->Branch("isGsfScPixConsistent", &isGsfScPixConsistent);
        ElectronTree->Branch("isGsfCtfConsistent", &isGsfCtfConsistent);
        ElectronTree->Branch("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        ElectronTree->Branch("dEtaInSeed", dEtaInSeed);
        ElectronTree->Branch("dPhiIn", &dPhiIn);
        ElectronTree->Branch("HoverE", &HoverE);
        ElectronTree->Branch("InvEminusInvP", &InvEminusInvP);
        ElectronTree->Branch("mHits", &mHits);
        ElectronTree->Branch("passConvVeto", &passConvVeto);
        ElectronTree->Branch("relPFiso_Rho", &relPFiso_Rho);
        ElectronTree->Branch("passMediumID", &passMediumID);
        ElectronTree->Branch("relECALiso", &relECALiso);
        ElectronTree->Branch("relHCALiso", &relHCALiso);
        ElectronTree->Branch("relTrkIso", &relTrkIso);
        ElectronTree->Branch("trig_fired", &trig_fired);
        ElectronTree->Branch("trig_matched", &trig_matched);
        ElectronTree->Branch("trig_pT", &trig_pT);
        ElectronTree->Branch("MET_pT", &MET_pT);
        ElectronTree->Branch("MET_phi", &MET_phi);
        ElectronTree->Branch("runNum", &runNum);
        ElectronTree->Branch("lumiBlock", &lumiBlock);
        ElectronTree->Branch("nPU", &nPU);
        ElectronTree->Branch("nVTX", &nVTX);
        ElectronTree->Branch("PVz", &PVz);
        ElectronTree->Branch("gen_weight", &gen_weight);
        ElectronTree->Branch("top_weight", &top_weight);
        ElectronTree->Branch("prefiring_weight", &prefiring_weight);
        ElectronTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        ElectronTree->Branch("prefiring_weight_down", &prefiring_weight_down);
        ElectronTree->Branch("HLT_trigPS", &HLT_trigPS);
        ElectronTree->Branch("L1_trigPS", &L1_trigPS);
//        ElectronTree->Branch("L1seed_trigPSinDetail", &L1seed_trigPSinDetail);

        Int_t currentRunNo=-999, currentLumiSection=-999;
//        std::map<std::pair<int,int>,int> lumis; // to search present run numbers and lumi sections

        // to count #events in different runs
        std::map<int, int> nEvtInRun;
        Int_t nEvt=0, nEvtG=0, nLumis=0, nEvtInLumi=0, minEvtInLumi=99999999, maxEvtInLumi=-9999999;

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            // TEST
            chain->Add(Mgr.FullLocation[i_tup]);
//            Mgr.SetupChain(i_tup, chain);

//            std::map<std::tuple<int,int,unsigned long long>, int> repeats; // to search for repeating events

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            // Loop for all events in the chain
//            for (Int_t i=0; i<NEvents; i++)
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);
/*
                //TEST/////////////////
                Int_t pass = 0;
                trig_fired->clear();
                trig_pT->clear();
                for( Int_t k = 0; k < ntuple->HLT_ntrig; k++ )
                {
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon22_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(22);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon30_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(30);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon36_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(36);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon50_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(50);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon75_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(75);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon90_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(90);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon120_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(120);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon175_v*" && ntuple->HLT_trigFired[k] )
                    {
                        pass = 1;
                        trig_fired->push_back(175);
                        trig_pT->push_back(ntuple->HLT_trigPt[k]);
                    }
                    if (pass)
                    {
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        ElectronTree->Fill();
                    }
                }
                // END TEST /////////////////////////////////////
*/

//                nEvt++;
//                Int_t skipLumi = skipLumiSection(ntuple->runNum, ntuple->lumiBlock);
//                Int_t skiprun = skipRun(ntuple->runNum);
//                if (skipLumi) continue;
//                if (skiprun) continue;
//                if (ntuple->runNum < 275900) continue;
//                nEvtG++;
                nEvtInRun[ntuple->runNum]++;

                if (ntuple->runNum != currentRunNo)
                {
                    currentRunNo = ntuple->runNum;
                    currentLumiSection = ntuple->lumiBlock;
                    nLumis++;
                    if (nEvtInLumi > maxEvtInLumi && nEvtInLumi > 0) maxEvtInLumi = nEvtInLumi;
                    if (nEvtInLumi < minEvtInLumi && nEvtInLumi > 0) minEvtInLumi = nEvtInLumi;
                    nEvtInLumi = 0;
//                    lumis[std::make_pair(currentRunNo, currentLumiSection)]++;
                }
                else if (ntuple->lumiBlock != currentLumiSection)
                {
                    currentLumiSection = ntuple->lumiBlock;
                    nLumis++;
                    if (nEvtInLumi > maxEvtInLumi && nEvtInLumi > 0) maxEvtInLumi = nEvtInLumi;
                    if (nEvtInLumi < minEvtInLumi && nEvtInLumi > 0) minEvtInLumi = nEvtInLumi;
                    nEvtInLumi = 0;
//                    lumis[std::make_pair(currentRunNo, currentLumiSection)]++;
                }
                nEvtInLumi++;

//                if (ntuple->pfMET_pT >= 20) continue;

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;
                if (Processes[i_proc] >= _GJets_20to100 && Processes[i_proc] <= _GJets_2000to5000)
                    gen_weight = ntuple->GENEvt_weight; // Resetting to +-1 affects normalization

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;


                Bool_t TriggerFlag = kFALSE;
                TString triggername;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT, &triggername);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Debug == kTRUE && triggername == "HLT_Photon175_v*") {NEvents++; continue;} // FOR TEST
                    if (Debug == kTRUE)
                    {
                        for (Int_t i_tr=0; i_tr<ntuple->HLT_ntrig; i_tr++)
                        {
                            cout << ntuple->HLT_trigName->at(i_tr) << "   prescale=" << ntuple->HLT_trigPS->at(i_tr) << endl;
                        }
                        cout << endl;
                    }
                    // -- Reco level selection -- //
                    vector< Electron > ElectronCollection;
                    Int_t NLeptons = ntuple->Nelectrons;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_reco);
                        ElectronCollection.push_back(ele);
                    }

                    // -- Event Selection -- //
                    vector< Electron > SelectedElectronCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FR(ElectronCollection, ntuple, &SelectedElectronCollection);

                    if (isPassEventSelection == kTRUE)
                    {
//                        repeats[std::make_tuple(ntuple->runNum, ntuple->lumiBlock, ntuple->evtNum)]++;
//                        if (repeats[std::make_tuple(ntuple->runNum, ntuple->lumiBlock, ntuple->evtNum)] > 1)
//                            cout << "Evt " << ntuple->runNum << "; " << ntuple->lumiBlock << "; " << ntuple->evtNum << " repeated " <<
//                                    repeats[std::make_tuple(ntuple->runNum, ntuple->lumiBlock, ntuple->evtNum)] << " times." << endl;

                        if (Debug == kTRUE) cout << "\nEvent " << i << endl << triggername << endl;
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
                        scPixCharge->clear();
                        isGsfCtfScPixConsistent->clear();
                        isGsfScPixConsistent->clear();
                        isGsfCtfConsistent->clear();
                        Full5x5_SigmaIEtaIEta->clear();
                        dEtaInSeed->clear();
                        dPhiIn->clear();
                        HoverE->clear();
                        InvEminusInvP->clear();
                        mHits->clear();
                        passConvVeto->clear();
                        relPFiso_Rho->clear();
                        passMediumID->clear();
                        relECALiso->clear();
                        relHCALiso->clear();
                        relTrkIso->clear();
                        trig_fired->clear();
                        trig_matched->clear();
                        trig_pT->clear();
                        HLT_trigPS->clear();
                        L1_trigPS->clear();
//                        L1seed_trigPSinDetail.clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- MET information -- //
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;

                        // -- Information for various other reweightings -- //
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        Int_t triggered = 0;
                        std::vector<int> trig_index;
                        triggered = analyzer->FindTrigger(SelectedElectronCollection, ntuple, &trig_index, trig_fired, trig_matched);

                        // -- Trigger vector filling -- //
                        if (!triggered) continue;
                        for (UInt_t i_tr=0; i_tr<trig_index.size(); i_tr++)
                        {
                            trig_pT->push_back(ntuple->HLT_trigPt[trig_index[i_tr]]);
                            HLT_trigPS->push_back(ntuple->HLT_trigPS->at(trig_index[i_tr]));
//                            L1seed_trigPS->push_back(ntuple->L1seed_trigPS->at(trig_index[i_tr]));
//                            L1seed_trigPSinDetail.push_back(ntuple->L1seed_trigPSinDetail->at(trig_index[i_tr]));
                        }
                        for (Int_t i_L1=0; i_L1<4; i_L1++)
                            L1_trigPS->push_back(ntuple->L1_trigPS->at(i_L1));

                        // -- Electron vector filling -- //
                        for (UInt_t i_ele=0; i_ele<SelectedElectronCollection.size(); i_ele++)
                        {
                            p_T->push_back(SelectedElectronCollection[i_ele].Pt);
                            eta->push_back(SelectedElectronCollection[i_ele].etaSC);
                            phi->push_back(SelectedElectronCollection[i_ele].phi);
                            charge->push_back(SelectedElectronCollection[i_ele].charge);
                            scPixCharge->push_back(SelectedElectronCollection[i_ele].scPixCharge);
                            isGsfCtfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfScPixConsistent);
                            isGsfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfScPixConsistent);
                            isGsfCtfConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfConsistent);
                            Full5x5_SigmaIEtaIEta->push_back(SelectedElectronCollection[i_ele].Full5x5_SigmaIEtaIEta);
                            dEtaInSeed->push_back(SelectedElectronCollection[i_ele].dEtaInSeed);
                            dPhiIn->push_back(SelectedElectronCollection[i_ele].dPhiIn);
                            HoverE->push_back(SelectedElectronCollection[i_ele].HoverE);
                            InvEminusInvP->push_back(SelectedElectronCollection[i_ele].InvEminusInvP);
                            mHits->push_back(SelectedElectronCollection[i_ele].mHits);
                            passConvVeto->push_back(SelectedElectronCollection[i_ele].passConvVeto);
                            relPFiso_Rho->push_back(SelectedElectronCollection[i_ele].RelPFIso_Rho);
                            passMediumID->push_back(SelectedElectronCollection[i_ele].passMediumID);
                            relECALiso->push_back(SelectedElectronCollection[i_ele].ecalIso03 / SelectedElectronCollection[i_ele].Pt);
                            relHCALiso->push_back(SelectedElectronCollection[i_ele].hcalIso03 / SelectedElectronCollection[i_ele].Pt);
                            relTrkIso->push_back(SelectedElectronCollection[i_ele].trkIso03 / SelectedElectronCollection[i_ele].Pt);

                            if (Debug == kTRUE)
                            {
                                cout << "Passing electron " << i_ele << ": p_T = " << SelectedElectronCollection[i_ele].Pt;
                                cout << "   eta = " << SelectedElectronCollection[i_ele].etaSC;
                                cout << "   charge = " << SelectedElectronCollection[i_ele].charge;
                                if (SelectedElectronCollection[i_ele].passMediumID == 1) cout << "   MediumID";
                                cout << endl;
                                if (i_ele < trig_fired->size())
                                    cout << "Trigger Photon" << trig_fired->at(i_ele) << " HLT pT: " << trig_pT->at(i_ele) << endl;
                            }
                        }
                        ElectronTree->Fill();

                    } // End of event selection

                } // End of if(isTriggered)

                if (!Debug) bar.Draw(i);
            } // End of event iteration

            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

//            Int_t noreps=0, reps=0;
//            for (std::map<std::tuple<int,int,unsigned long long>,int>::const_iterator it=repeats.begin(); it!=repeats.end(); it++)
//            {
//                if (it->second == 1)noreps++;
//                else if (it->second == 0) cout << "0 repetitions.. What??" << endl;
//                else reps++;
//                if (Debug)
//                {
//                    cout << "RunNo=" << std::get<0>(it->first) << "  LumiSec=" << std::get<1>(it->first) << "  EvtNo=" << std::get<2>(it->first);
//                    cout << "   repeats " << it->second << " times." << endl;
//                }
//            }
//            cout << noreps << " events without repetitions\n" << reps << " events with repetitions" << endl;

        } // End of i_tup iteration

//        cout << "Full #events: " << nEvt << endl;
//        cout << "Good #events: " << nEvtG << endl;
//        cout << "#events by run:" << endl;
//        for (std::map<int,int>::const_iterator it=nEvtInRun.begin(); it!=nEvtInRun.end(); it++)
//        {
//            cout << "   Run " << it->first << ":  " << it->second << " events" << endl;
//        }
//        cout << "Avg #events in lumisection: " << ((float)nEvtG)/((float)nLumis) << endl;
//        cout << "Min #events in lumisection: " << minEvtInLumi << endl;
//        cout << "Max #events in lumisection: " << maxEvtInLumi << "\n" << endl;

//        ofstream runOutput;
//        runOutput.open("Lumi_"+Mgr.Procname[Mgr.CurrentProc]+".txt");
//        runOutput << "RunNo\tLumiSection\n";

//        for (Int_t rn=270000; rn<290000; rn++)
//        {
//            Int_t printed = 0;
//            for (Int_t ls=1; ls<4000; ls++)
//            {
//                for (std::map<std::pair<int,int>,int>::const_iterator it=lumis.begin(); it!=lumis.end(); it++)
//                {
//                    if (std::get<0>(it->first) != rn) continue;
//                    if (std::get<1>(it->first) != ls) continue;
//                    if (!printed) {printed++; runOutput << rn << "\n";}
//                    runOutput << "\t" << ls << "\n";
//                }
//            }
//        }

//        runOutput.close();

        // Writing
        cout << "Writing into file...";
        ElectronFile->cd();
        Int_t write;
        write = ElectronTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            ElectronFile->Close();
            if (!ElectronFile->IsOpen()) cout << "File SelectedForFR_E_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForFR_E_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!\n" << endl;
            ElectronFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForFR_E


void MakeSelectionForFR_E_alt (TString type, TString HLTname , Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* ElectronFile;
        out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
        out_dir = "SelectedForFR_E_alt_"+Mgr.Procname[Mgr.CurrentProc];

        if (Debug == kTRUE)
            ElectronFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            ElectronFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        ElectronFile->cd();

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
        Int_t runNum;
        Int_t lumiBlock;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        // -- Creating electron variables to assign branches -- //
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("etaSC", &etaSC);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
        ElectronTree->Branch("scPixCharge", &scPixCharge);
        ElectronTree->Branch("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        ElectronTree->Branch("isGsfScPixConsistent", &isGsfScPixConsistent);
        ElectronTree->Branch("isGsfCtfConsistent", &isGsfCtfConsistent);
        ElectronTree->Branch("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        ElectronTree->Branch("dEtaInSeed", dEtaInSeed);
        ElectronTree->Branch("dPhiIn", &dPhiIn);
        ElectronTree->Branch("HoverE", &HoverE);
        ElectronTree->Branch("InvEminusInvP", &InvEminusInvP);
        ElectronTree->Branch("chIso03", &chIso03);
        ElectronTree->Branch("nhIso03", &nhIso03);
        ElectronTree->Branch("phIso03", &phIso03);
        ElectronTree->Branch("ChIso03FromPU", &ChIso03FromPU);
        ElectronTree->Branch("mHits", &mHits);
        ElectronTree->Branch("passConvVeto", &passConvVeto);
        ElectronTree->Branch("relPFiso_dBeta", &relPFiso_dBeta);
        ElectronTree->Branch("relPFiso_Rho", &relPFiso_Rho);
        ElectronTree->Branch("passMediumID", &passMediumID);
        ElectronTree->Branch("MET_pT", &MET_pT);
        ElectronTree->Branch("MET_phi", &MET_phi);
        ElectronTree->Branch("runNum", &runNum);
        ElectronTree->Branch("lumiBlock", &lumiBlock);
        ElectronTree->Branch("nPU", &nPU);
        ElectronTree->Branch("nVTX", &nVTX);
        ElectronTree->Branch("PVz", &PVz);
        ElectronTree->Branch("gen_weight", &gen_weight);
        ElectronTree->Branch("top_weight", &top_weight);
        ElectronTree->Branch("prefiring_weight", &prefiring_weight);
        ElectronTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        ElectronTree->Branch("prefiring_weight_down", &prefiring_weight_down);

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;
                if (Processes[i_proc] >= _GJets_20to100 && Processes[i_proc] <= _GJets_2000to5000)
                    gen_weight = ntuple->GENEvt_weight; // Resetting to +-1 affects normalization

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kFALSE;
                TString triggername;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT, &triggername);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    // -- Reco level selection -- //
                    vector< Electron > ElectronCollection;
                    Int_t NLeptons = ntuple->Nelectrons;
                    if (NLeptons < 3) continue;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_reco);
                        ElectronCollection.push_back(ele);
                    }

                    // -- Event Selection -- //
                    vector< Electron > SelectedElectronCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FR_alt(ElectronCollection, ntuple, &SelectedElectronCollection);

                    if (isPassEventSelection == kTRUE)
                    {
                        if (Debug == kTRUE) cout << "\nEvent " << i << endl << triggername << endl;
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        etaSC->clear();
                        phi->clear();
                        charge->clear();
                        scPixCharge->clear();
                        isGsfCtfScPixConsistent->clear();
                        isGsfScPixConsistent->clear();
                        isGsfCtfConsistent->clear();
                        Full5x5_SigmaIEtaIEta->clear();
                        dEtaInSeed->clear();
                        dPhiIn->clear();
                        HoverE->clear();
                        InvEminusInvP->clear();
                        chIso03->clear();
                        nhIso03->clear();
                        phIso03->clear();
                        ChIso03FromPU->clear();
                        mHits->clear();
                        passConvVeto->clear();
                        relPFiso_dBeta->clear();
                        relPFiso_Rho->clear();
                        passMediumID->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- MET information -- //
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;

                        // -- Information for various other reweightings -- //
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        // -- Vector filling -- //
                        for (UInt_t i_ele=0; i_ele<SelectedElectronCollection.size(); i_ele++)
                        {
                            p_T->push_back(SelectedElectronCollection[i_ele].Pt);
                            eta->push_back(SelectedElectronCollection[i_ele].eta);
                            etaSC->push_back(SelectedElectronCollection[i_ele].etaSC);
                            phi->push_back(SelectedElectronCollection[i_ele].phi);
                            charge->push_back(SelectedElectronCollection[i_ele].charge);
                            scPixCharge->push_back(SelectedElectronCollection[i_ele].scPixCharge);
                            isGsfCtfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfScPixConsistent);
                            isGsfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfScPixConsistent);
                            isGsfCtfConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfConsistent);
                            Full5x5_SigmaIEtaIEta->push_back(SelectedElectronCollection[i_ele].Full5x5_SigmaIEtaIEta);
                            dEtaInSeed->push_back(SelectedElectronCollection[i_ele].dEtaInSeed);
                            dPhiIn->push_back(SelectedElectronCollection[i_ele].dPhiIn);
                            HoverE->push_back(SelectedElectronCollection[i_ele].HoverE);
                            InvEminusInvP->push_back(SelectedElectronCollection[i_ele].InvEminusInvP);
                            chIso03->push_back(SelectedElectronCollection[i_ele].chIso03);
                            nhIso03->push_back(SelectedElectronCollection[i_ele].nhIso03);
                            phIso03->push_back(SelectedElectronCollection[i_ele].phIso03);
                            ChIso03FromPU->push_back(SelectedElectronCollection[i_ele].ChIso03FromPU);
                            mHits->push_back(SelectedElectronCollection[i_ele].mHits);
                            passConvVeto->push_back(SelectedElectronCollection[i_ele].passConvVeto);
                            relPFiso_dBeta->push_back(SelectedElectronCollection[i_ele].RelPFIso_dBeta);
                            relPFiso_Rho->push_back(SelectedElectronCollection[i_ele].RelPFIso_Rho);
                            passMediumID->push_back(SelectedElectronCollection[i_ele].passMediumID);

                            if (Debug == kTRUE)
                            {
                                cout << "Passing electron " << i_ele << ": p_T = " << SelectedElectronCollection[i_ele].Pt;
                                cout << "   eta = " << SelectedElectronCollection[i_ele].eta;
                                cout << "   etaSC = " << SelectedElectronCollection[i_ele].etaSC;
                                cout << "   charge = " << SelectedElectronCollection[i_ele].charge;
                                if (SelectedElectronCollection[i_ele].passMediumID == 1) cout << "   MediumID";
                                cout << endl;
                            }
                        }
                        if (Debug == kTRUE)
                        {
                            Double_t mass = (SelectedElectronCollection[0].Momentum + SelectedElectronCollection[1].Momentum).M();
                            cout << "mass = " << mass << endl;
                        }
                        ElectronTree->Fill();

                    } // End of event selection

                } // End of if(isTriggered)

                if (!Debug) bar.Draw(i);
            } // End of event iteration

            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into file...";
        ElectronFile->cd();
        Int_t write;
        write = ElectronTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            ElectronFile->Close();
            if (!ElectronFile->IsOpen()) cout << "File SelectedForFR_E_alt_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForFR_E_alt_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!\n" << endl;
            ElectronFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForFR_E_alt


/// -------------------------------- Muon ------------------------------------ ///
void MakeSelectionForFR_Mu (TString type, TString HLTname, Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);
    // -- For Rochester correction -- //
    TRandom3 *r1 = new TRandom3(0);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* MuonFile;
        if (Mgr.Type == "DATA")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "SIGNAL")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "BKG")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "TEST")
        {
            out_base = "/media/sf_DATA/test/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else
        {
            cout << "Problems with TYPE." << endl;
            return;
        }

        if (Debug == kTRUE)
            MuonFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            MuonFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        MuonFile->cd();

        Int_t nMuons;
        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        std::vector<int> *isGenMatched = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
        Int_t nElectrons;
        std::vector<double> *ele_pT = new std::vector<double>;
        std::vector<double> *ele_eta = new std::vector<double>;
        std::vector<double> *ele_etaSC = new std::vector<double>;
        std::vector<double> *ele_phi = new std::vector<double>;
        std::vector<int> *ele_charge = new std::vector<int>;
        Int_t nJets;
        std::vector<double> *jet_pT = new std::vector<double>;
        std::vector<double> *jet_eta = new std::vector<double>;
        std::vector<double> *jet_phi = new std::vector<double>;
        std::vector<int> *jet_charge = new std::vector<int>;

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating muon variables to assign branches -- //
        MuonTree->Branch("nMuons", &nMuons);
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("phi", &phi);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("TRKiso", &TRKiso);
        MuonTree->Branch("isGenMatched", &isGenMatched);
        MuonTree->Branch("MET_pT", &MET_pT);
        MuonTree->Branch("MET_phi", &MET_phi);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);
        MuonTree->Branch("nElectrons", &nElectrons);
        MuonTree->Branch("ele_pT", &ele_pT);
        MuonTree->Branch("ele_eta", &ele_eta);
        MuonTree->Branch("ele_etaSC", &ele_etaSC);
        MuonTree->Branch("ele_phi", &ele_phi);
        MuonTree->Branch("ele_charge", &ele_charge);
        MuonTree->Branch("nJets", &nJets);
        MuonTree->Branch("jet_pT", &jet_pT);
        MuonTree->Branch("jet_eta", &jet_eta);
        MuonTree->Branch("jet_phi", &jet_phi);
        MuonTree->Branch("jet_charge", &jet_charge);

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
//            chain->Add(Mgr.FullLocation[i_tup]+"*.root");
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Muon();
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_Jet();

            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 100; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            std::string RCaddress;
            if (Mgr.Type == "TEST")
                RCaddress = "./etc/RoccoR/rcdata.2016.v3";
            else RCaddress = "/cms/ldap_home/mambroza/DrellYan2016/SelectedX/etc/RoccoR/rcdata.2016.v3";
            RoccoR rc(RCaddress);

            std::vector<string> HLTs;

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                if (Debug)
                {
                    for(Int_t i_trig=0; i_trig<ntuple->HLT_ntrig; i_trig++)
                    {
                        Int_t written_already = 0;
                        if (HLTs.size())
                        {
                            for(UInt_t i_tr=0; i_tr<HLTs.size(); i_tr++)
                            {
                                if (HLTs[i_tr] == ntuple->HLT_trigName->at(i_trig))
                                    written_already = 1;
                            }
                        }
                        if (!written_already)
                            HLTs.push_back(ntuple->HLT_trigName->at(i_trig));
                    }
                }

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    // -- Reco level selection -- //
                    vector< Muon > MuonCollection;
                    Int_t NLeptons = ntuple->nMuon;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("finalState_OR_hadronDecay", ntuple, mu);
                            if (genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                            if (mu.Pt != mu.Pt || SF != SF || SF*mu.Pt != SF*mu.Pt)
                                cout << "\nGenPt: " << genPt << "  Pt: " << mu.Pt << "  Corr Pt: " << SF*mu.Pt << "  multiplier: " << SF
                                     << "\nEta: " << mu.eta << "   Phi: " << mu.phi << "  Charge: " << mu.charge << "\ntrLayers: " << mu.trackerLayers
                                     << "  PFiso: " << mu.relPFiso << endl;
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);

                    } // End of i_reco iteration

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection_nume, SelectedMuonCollection_deno;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FR(MuonCollection, ntuple, &SelectedMuonCollection_nume, &SelectedMuonCollection_deno);

                    if (isPassEventSelection == kTRUE)
                    {
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
                        relPFiso->clear();
                        TRKiso->clear();
                        isGenMatched->clear();
                        ele_pT->clear();
                        ele_eta->clear();
                        ele_etaSC->clear();
                        ele_phi->clear();
                        ele_charge->clear();
                        jet_pT->clear();
                        jet_eta->clear();
                        jet_phi->clear();
                        jet_charge->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- MET information -- //
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;

                        // -- Information for various other reweightings -- //
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        nMuons = SelectedMuonCollection_deno.size();
                        // -- Vector filling -- //
                        for (UInt_t i=0; i<SelectedMuonCollection_deno.size(); i++)
                        {
                            p_T->push_back(SelectedMuonCollection_deno[i].Pt);
                            eta->push_back(SelectedMuonCollection_deno[i].eta);
                            phi->push_back(SelectedMuonCollection_deno[i].phi);
                            charge->push_back(SelectedMuonCollection_deno[i].charge);
                            relPFiso->push_back(SelectedMuonCollection_deno[i].RelPFIso_dBeta);
                            TRKiso->push_back(SelectedMuonCollection_deno[i].trkiso);
                            isGenMatched->push_back(analyzer->isGenMatched(SelectedMuonCollection_deno[i], "fromHardProcess", ntuple));
                        }
                        nElectrons = ntuple->Nelectrons;
                        for (Int_t i_ele=0; i_ele<nElectrons; i_ele++)
                        {
                            ele_pT->push_back(ntuple->Electron_pT[i_ele]);
                            ele_eta->push_back(ntuple->Electron_eta[i_ele]);
                            ele_etaSC->push_back(ntuple->Electron_etaSC[i_ele]);
                            ele_phi->push_back(ntuple->Electron_phi[i_ele]);
                            ele_charge->push_back(ntuple->Electron_charge[i_ele]);
                        }
                        nJets = ntuple->Njets;
                        for (Int_t i_jet=0; i_jet<nJets; i_jet++)
                        {
                            jet_pT->push_back(ntuple->Jet_pT[i_jet]);
                            jet_eta->push_back(ntuple->Jet_eta[i_jet]);
                            jet_phi->push_back(ntuple->Jet_phi[i_jet]);
                            jet_charge->push_back(ntuple->Jet_Charge[i_jet]);
                        }

                        MuonTree->Fill();
                    } // End of isPassEvtSelection

                } // End of if(isTriggered)

                bar.Draw(i);

            } // End of event iteration
            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            if (Debug)
            {
                cout << "Available trigger names:" << endl;
                if (HLTs.size())
                {
                    for(UInt_t i_tr=0; i_tr<HLTs.size(); i_tr++)
                        cout << HLTs[i_tr] << endl;
                    cout << endl;
                }
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into files...";
        MuonFile->cd();
        Int_t write;
        write = MuonTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            MuonFile->Close();
            if (!MuonFile->IsOpen()) cout << "File SelectedForFR_Mu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForFR_Mu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            MuonFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForFR_Mu


/////////////////////////////////////////////////////////////////////////////
/// -------------------------- BKG estimation --------------------------- ///
/// /////////////////////////////////////////////////////////////////////////

/// Electron QCD+WJET estimation
void MakeSelectionForBKGest_EE (TString type, TString HLTname, Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* ElectronFile;
        if (Mgr.Type == "DATA")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "SIGNAL")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "BKG")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "TEST")
        {
            out_base = "/media/sf_DATA/test/";
            out_dir = "SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else
        {
            cout << "Problems with TYPE." << endl;
            return;
        }

        if (Debug == kTRUE)
            ElectronFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            ElectronFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        ElectronFile->cd();

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
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<double> *relTrkIso = new std::vector<double>;
        std::vector<double> *relECALiso = new std::vector<double>;
        std::vector<double> *relHCALiso = new std::vector<double>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("etaSC", &etaSC);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
        ElectronTree->Branch("scPixCharge", &scPixCharge);
        ElectronTree->Branch("isGsfCtfScPixConsistent", &isGsfCtfScPixConsistent);
        ElectronTree->Branch("isGsfScPixConsistent", &isGsfScPixConsistent);
        ElectronTree->Branch("isGsfCtfConsistent", &isGsfCtfConsistent);
        ElectronTree->Branch("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        ElectronTree->Branch("dEtaInSeed", &dEtaInSeed);
        ElectronTree->Branch("dPhiIn", &dPhiIn);
        ElectronTree->Branch("HoverE", &HoverE);
        ElectronTree->Branch("InvEminusInvP", &InvEminusInvP);
        ElectronTree->Branch("mHits", &mHits);
        ElectronTree->Branch("passConvVeto", &passConvVeto);
        ElectronTree->Branch("relPFiso_Rho", &relPFiso_Rho);
        ElectronTree->Branch("passMediumID", &passMediumID);
        ElectronTree->Branch("relTrkIso", &relTrkIso);
        ElectronTree->Branch("relECALiso", &relECALiso);
        ElectronTree->Branch("relHCALiso", &relHCALiso);
        ElectronTree->Branch("MET_pT", &MET_pT);
        ElectronTree->Branch("MET_phi", &MET_phi);
        ElectronTree->Branch("nPU", &nPU);
        ElectronTree->Branch("nVTX", &nVTX);
        ElectronTree->Branch("PVz", &PVz);
        ElectronTree->Branch("gen_weight", &gen_weight);
        ElectronTree->Branch("top_weight", &top_weight);
        ElectronTree->Branch("prefiring_weight", &prefiring_weight);
        ElectronTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        ElectronTree->Branch("prefiring_weight_down", &prefiring_weight_down);

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 100; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;
                if (Processes[i_proc] >= _GJets_20to100 && Processes[i_proc] <= _GJets_2000to5000)
                    gen_weight = ntuple->GENEvt_weight; // Resetting to +-1 affects normalization

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    // -- Reco level selection -- //
                    vector< Electron > ElectronCollection;
                    Int_t NLeptons = ntuple->Nelectrons;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_reco);
                        ElectronCollection.push_back(ele);

                    } // End of i_reco iteration

                    // -- Event Selection -- //
                    vector< Electron > SelectedElectronCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FakeElectrons(ElectronCollection, ntuple, &SelectedElectronCollection);

                    if (isPassEventSelection == kTRUE)
                    {
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        etaSC->clear();
                        phi->clear();
                        charge->clear();
                        scPixCharge->clear();
                        isGsfCtfScPixConsistent->clear();
                        isGsfScPixConsistent->clear();
                        isGsfCtfConsistent->clear();
                        Full5x5_SigmaIEtaIEta->clear();
                        dEtaInSeed->clear();
                        dPhiIn->clear();
                        HoverE->clear();
                        InvEminusInvP->clear();
                        mHits->clear();
                        passConvVeto->clear();
                        relPFiso_Rho->clear();
                        passMediumID->clear();
                        relTrkIso->clear();
                        relECALiso->clear();
                        relHCALiso->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // MET information
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;

                        // Information for various other reweightings
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        // -- Vector filling -- //
                        for (UInt_t i_ele=0; i_ele<SelectedElectronCollection.size(); i_ele++)
                        {
                            p_T->push_back(SelectedElectronCollection[i_ele].Pt);
                            eta->push_back(SelectedElectronCollection[i_ele].eta);
                            etaSC->push_back(SelectedElectronCollection[i_ele].etaSC);
                            phi->push_back(SelectedElectronCollection[i_ele].phi);
                            charge->push_back(SelectedElectronCollection[i_ele].charge);
                            scPixCharge->push_back(SelectedElectronCollection[i_ele].scPixCharge);
                            isGsfCtfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfScPixConsistent);
                            isGsfScPixConsistent->push_back(SelectedElectronCollection[i_ele].isGsfScPixConsistent);
                            isGsfCtfConsistent->push_back(SelectedElectronCollection[i_ele].isGsfCtfConsistent);
                            Full5x5_SigmaIEtaIEta->push_back(SelectedElectronCollection[i_ele].Full5x5_SigmaIEtaIEta);
                            dEtaInSeed->push_back(SelectedElectronCollection[i_ele].dEtaInSeed);
                            dPhiIn->push_back(SelectedElectronCollection[i_ele].dPhiIn);
                            HoverE->push_back(SelectedElectronCollection[i_ele].HoverE);
                            InvEminusInvP->push_back(SelectedElectronCollection[i_ele].InvEminusInvP);
                            mHits->push_back(SelectedElectronCollection[i_ele].mHits);
                            passConvVeto->push_back(SelectedElectronCollection[i_ele].passConvVeto);
                            relPFiso_Rho->push_back(SelectedElectronCollection[i_ele].RelPFIso_Rho);
                            passMediumID->push_back(SelectedElectronCollection[i_ele].passMediumID);
                            relTrkIso->push_back(SelectedElectronCollection[i_ele].trkIso03 / SelectedElectronCollection[i_ele].Pt);
                            relECALiso->push_back(SelectedElectronCollection[i_ele].ecalIso03 / SelectedElectronCollection[i_ele].Pt);
                            relHCALiso->push_back(SelectedElectronCollection[i_ele].hcalIso03 / SelectedElectronCollection[i_ele].Pt);
                        }
                        ElectronTree->Fill();
                    } // End of isPassEvtSelection

                } // End of if(isTriggered)

                bar.Draw(i);

            } // End of event iteration
            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into files...";
        ElectronFile->cd();
        Int_t write;
        write = ElectronTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            ElectronFile->Close();
            if (!ElectronFile->IsOpen()) cout << "File SelectedForBKGest_E_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForBKGest_E_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            ElectronFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForBKGest_EE

/// Muon QCD+WJET estimation
void MakeSelectionForBKGest_MuMu (TString type, TString HLTname, Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);
    // -- For Rochester correction -- //
    TRandom3 *r1 = new TRandom3(0);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* MuonFile;
        if (Mgr.Type == "DATA")
        {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "Data/SelectedForBKGest_Mu_Triggerless_"+Mgr.Tag[i_tup];

            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "SIGNAL")
        {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_signal/SelectedForBKGest_Mu_Triggerless_"+Mgr.Tag[i_tup];

            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "BKG")
        {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_bkg/SelectedForBKGest_Mu_Triggerless_"+Mgr.Tag[i_tup];

            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "TEST")
        {
            out_base = "/media/sf_DATA/test/";
            out_dir = "SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else
        {
            cout << "Problems with TYPE." << endl;
            return;
        }

        if (Debug == kTRUE)
            MuonFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            MuonFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        MuonFile->cd();

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
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
        Int_t runNum;
        Int_t lumiBlock;

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("phi", &phi);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("TRKiso", &TRKiso);
        MuonTree->Branch("trig_name", &trig_name);
        MuonTree->Branch("trig_fired", &trig_fired);
        MuonTree->Branch("trig_matched", &trig_matched);
        MuonTree->Branch("trig_pT", &trig_pT);
        MuonTree->Branch("prescale_factor", &prescale_factor);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);
        MuonTree->Branch("runNum", &runNum);
        MuonTree->Branch("lumiBlock", &lumiBlock);

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Muon();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 328; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            std::string RCaddress;
            if (Mgr.Type == "TEST")
                RCaddress = "./etc/RoccoR/rcdata.2016.v3";
            else RCaddress = "/cms/ldap_home/mambroza/DrellYan2016/SelectedX/etc/RoccoR/rcdata.2016.v3";
            RoccoR rc(RCaddress);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                if (!Debug) bar.Draw(i);
                ntuple->GetEvent(i);
                if (ntuple->nMuon < 2) continue;
                Int_t nTight = 0;
                for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                {
                    if (ntuple->Muon_passTightID[i_mu] && ntuple->Muon_pT[i_mu] > 15) nTight++;
                }
                if (nTight < 2) continue;

                if (Debug == kTRUE)
                {
                    cout << "\nEvt " << i << endl;
                    cout << "Nmuons: " << ntuple->nMuon << ",  Ntrig: " << ntuple->HLT_ntrig << endl;
                    for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                    {
                        cout << "\tmu" << i_mu << ": pT=" << ntuple->Muon_pT[i_mu] << "  eta=" << ntuple->Muon_eta[i_mu];
                        Double_t PFiso = (ntuple->Muon_PfChargedHadronIsoR04[i_mu] +
                                 max(0.0, ntuple->Muon_PfNeutralHadronIsoR04[i_mu] + ntuple->Muon_PfGammaIsoR04[i_mu] - 0.5*ntuple->Muon_PFSumPUIsoR04[i_mu]))
                                 / ntuple->Muon_pT[i_mu];
                        cout << "  isTight=" << ntuple->Muon_passTightID[i_mu] << "  PFiso=" << PFiso << endl;
                    }
                }

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kTRUE;
//                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Debug == kTRUE)
                        cout << "PASS" << endl;
                    // -- Reco level selection -- //
                    vector< Muon > MuonCollection;
                    Int_t NLeptons = ntuple->nMuon;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("finalState_OR_hadronDecay", ntuple, mu);
                            if (genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                            if (mu.Pt != mu.Pt || SF != SF || SF*mu.Pt != SF*mu.Pt)
                                cout << "\nGenPt: " << genPt << "  Pt: " << mu.Pt << "  Corr Pt: " << SF*mu.Pt << "  multiplier: " << SF
                                     << "\nEta: " << mu.eta << "   Phi: " << mu.phi << "  Charge: " << mu.charge << "\ntrLayers: " << mu.trackerLayers
                                     << "  PFiso: " << mu.relPFiso << endl;
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);

                    } // End of i_reco iteration
                    if (Debug == kTRUE) cout << "Muons in collection: " << MuonCollection.size() << endl;

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FakeMuons(MuonCollection, ntuple, &SelectedMuonCollection);

                    if (isPassEventSelection == kTRUE)
                    {
                        if (Debug == kTRUE) cout << "Selection passed" << endl;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
                        relPFiso->clear();
                        TRKiso->clear();
                        trig_name->clear();
                        trig_fired->clear();
                        trig_matched->clear();
                        trig_pT->clear();
                        prescale_factor->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- Information for various other reweightings -- //
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        Int_t triggered = 0;
                        triggered = analyzer->FindTriggerAndPrescale2(SelectedMuonCollection, ntuple, trig_name, trig_fired, prescale_factor,
                                                                      trig_matched, trig_pT, Debug);
                        if (Debug == kTRUE) cout << "triggered=" << triggered << endl;
                        if (!triggered) continue;
                        timesPassed++;

                        // -- Vector filling -- //
                        for (UInt_t i=0; i<SelectedMuonCollection.size(); i++)
                        {
                            p_T->push_back(SelectedMuonCollection[i].Pt);
                            eta->push_back(SelectedMuonCollection[i].eta);
                            phi->push_back(SelectedMuonCollection[i].phi);
                            charge->push_back(SelectedMuonCollection[i].charge);
                            relPFiso->push_back(SelectedMuonCollection[i].RelPFIso_dBeta);
                            TRKiso->push_back(SelectedMuonCollection[i].trkiso);
                        }
                        MuonTree->Fill();
                    } // End of isPassEvtSelection

                } // End of if(isTriggered)                

            } // End of event iteration
            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into files...";
        MuonFile->cd();
        Int_t write;
        write = MuonTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            MuonFile->Close();
            if (!MuonFile->IsOpen()) cout << "File SelectedForBKGest_Mu_Triggerless_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForBKGest_Mu_Triggerless_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            MuonFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForBKGest_MuMu


void MakeSelectionForBKGest_EMu (TString type, TString HLTname, Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);
    // -- For Rochester correction -- //
    TRandom3 *r1 = new TRandom3(0);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
        TString out_dir = "SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc];
        TFile* MuonFile;

        if (Debug == kTRUE)
            MuonFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            MuonFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        MuonFile->cd();

        Double_t e_p_T, e_eta, e_etaSC, e_phi;
        Int_t e_charge;
        Double_t e_Full5x5_SigmaIEtaIEta, e_dEtaInSeed, e_dPhiIn, e_HoverE, e_InvEminusInvP;
        Double_t e_relPFiso_dBeta, e_relPFiso_Rho;
        Double_t e_chIso03, e_nhIso03, e_phIso03, e_ChIso03FromPU;
        Int_t e_mHits, e_passConvVeto, e_passMediumID;
        Double_t mu_p_T, mu_eta, mu_phi;
        Int_t mu_charge;
        Double_t mu_relPFiso_dBeta;
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

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        MuonTree->Branch("e_p_T", &e_p_T);
        MuonTree->Branch("e_eta", &e_eta);
        MuonTree->Branch("e_etaSC", &e_etaSC);
        MuonTree->Branch("e_phi", &e_phi);
        MuonTree->Branch("e_charge", &e_charge);
        MuonTree->Branch("e_Full5x5_SigmaIEtaIEta", &e_Full5x5_SigmaIEtaIEta);
        MuonTree->Branch("e_dEtaInSeed", &e_dEtaInSeed);
        MuonTree->Branch("e_dPhiIn", &e_dPhiIn);
        MuonTree->Branch("e_HoverE", &e_HoverE);
        MuonTree->Branch("e_InvEminusInvP", &e_InvEminusInvP);
        MuonTree->Branch("e_chIso03", &e_chIso03);
        MuonTree->Branch("e_nhIso03", &e_nhIso03);
        MuonTree->Branch("e_phIso03", &e_phIso03);
        MuonTree->Branch("e_ChIso03FromPU", &e_ChIso03FromPU);
        MuonTree->Branch("e_mHits", &e_mHits);
        MuonTree->Branch("e_passConvVeto", &e_passConvVeto);
        MuonTree->Branch("e_relPFiso_dBeta", &e_relPFiso_dBeta);
        MuonTree->Branch("e_relPFiso_Rho", &e_relPFiso_Rho);
        MuonTree->Branch("e_passMediumID", &e_passMediumID);
        MuonTree->Branch("mu_p_T", &mu_p_T);
        MuonTree->Branch("mu_eta", &mu_eta);
        MuonTree->Branch("mu_phi", &mu_phi);
        MuonTree->Branch("mu_charge", &mu_charge);
        MuonTree->Branch("mu_relPFiso_dBeta", &mu_relPFiso_dBeta);
        MuonTree->Branch("trig_fired", &trig_fired);
        MuonTree->Branch("trig_pT", &trig_pT);
        MuonTree->Branch("prescale_factor", &prescale_factor);
        MuonTree->Branch("MET_pT", &MET_pT);
        MuonTree->Branch("MET_phi", &MET_phi);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);
        MuonTree->Branch("runNum", &runNum);
        MuonTree->Branch("lumiBlock", &lumiBlock);

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_Muon();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 100; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            std::string RCaddress;
            if (Mgr.Type == "TEST")
                RCaddress = "./etc/RoccoR/rcdata.2016.v3";
            else RCaddress = "/cms/ldap_home/mambroza/DrellYan2016/SelectedX/etc/RoccoR/rcdata.2016.v3";
            RoccoR rc(RCaddress);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;
                if (Processes[i_proc] >= _GJets_20to100 && Processes[i_proc] <= _GJets_2000to5000)
                    gen_weight = ntuple->GENEvt_weight; // Resetting to +-1 affects normalization

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kTRUE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    // -- Reco level selection -- //
                    vector<Muon> MuonCollection;
                    vector<Electron> ElectronCollection;
                    for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_mu);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("finalState_OR_hadronDecay", ntuple, mu);
                            if (genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                            if (mu.Pt != mu.Pt || SF != SF || SF*mu.Pt != SF*mu.Pt)
                                cout << "\nGenPt: " << genPt << "  Pt: " << mu.Pt << "  Corr Pt: " << SF*mu.Pt << "  multiplier: " << SF
                                     << "\nEta: " << mu.eta << "   Phi: " << mu.phi << "  Charge: " << mu.charge << "\ntrLayers: " << mu.trackerLayers
                                     << "  PFiso: " << mu.relPFiso << endl;
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);
                    } // End of i_mu iteration

                    for (Int_t i_ele=0; i_ele<ntuple->Nelectrons; i_ele++)
                    {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_ele);
                        ElectronCollection.push_back(ele);
                    } // End of i_ele iteration

                    // -- Event Selection -- //
                    Muon SelectedMuon;
                    Electron SelectedElectron;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FakeEMu(ElectronCollection, MuonCollection, ntuple, &SelectedElectron, &SelectedMuon);

                    if (isPassEventSelection == kTRUE)
                    {
                        if (SelectedElectron.passMediumID && SelectedMuon.RelPFIso_dBeta < 0.15) continue;

                        timesPassed++;
                        trig_fired->clear();
                        prescale_factor->clear();
                        trig_pT->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // Trigger info
                        Int_t triggered = 0;
                        vector<Muon> SelectedMuonCollection;
                        SelectedMuonCollection.push_back(SelectedMuon);
                        std::vector<int> *trig_m = new std::vector<int>;
                        triggered = analyzer->FindTriggerAndPrescale(SelectedMuonCollection, ntuple, trig_fired, prescale_factor, trig_m, trig_pT);
                        if (!triggered) continue;

                        // Information for reweightings
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;
                        // MET information
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;
                        // EMu information
                        e_p_T = SelectedElectron.Pt;
                        e_eta = SelectedElectron.eta;
                        e_etaSC = SelectedElectron.etaSC;
                        e_phi = SelectedElectron.phi;
                        e_charge = SelectedElectron.charge;
                        e_Full5x5_SigmaIEtaIEta = SelectedElectron.Full5x5_SigmaIEtaIEta;
                        e_dEtaInSeed = SelectedElectron.dEtaInSeed;
                        e_dPhiIn = SelectedElectron.dPhiIn;
                        e_HoverE = SelectedElectron.HoverE;
                        e_InvEminusInvP = SelectedElectron.InvEminusInvP;
                        e_chIso03 = SelectedElectron.chIso03;
                        e_nhIso03 = SelectedElectron.nhIso03;
                        e_phIso03 = SelectedElectron.phIso03;
                        e_ChIso03FromPU = SelectedElectron.ChIso03FromPU;
                        e_mHits = SelectedElectron.mHits;
                        e_passConvVeto = SelectedElectron.passConvVeto;
                        e_relPFiso_dBeta = SelectedElectron.RelPFIso_dBeta;
                        e_relPFiso_Rho = SelectedElectron.RelPFIso_Rho;
                        e_passMediumID = SelectedElectron.passMediumID;
                        mu_p_T = SelectedMuon.Pt;
                        mu_eta = SelectedMuon.eta;
                        mu_phi = SelectedMuon.phi;
                        mu_charge = SelectedMuon.charge;
                        mu_relPFiso_dBeta = SelectedMuon.RelPFIso_dBeta;

                        MuonTree->Fill();
                    } // End of isPassEvtSelection

                } // End of if(isTriggered)

                bar.Draw(i);

            } // End of event iteration
            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into files...";
        MuonFile->cd();
        Int_t write;
        write = MuonTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            MuonFile->Close();
            if (!MuonFile->IsOpen()) cout << "File SelectedForBKGest_EMu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForBKGest_EMu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            MuonFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForBKGest_EMu


/// Muon PROMPT RATE estimation
void MakeSelectionForPR_Mu (TString type, TString HLTname, Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);
    // -- For Rochester correction -- //
    TRandom3 *r1 = new TRandom3(0);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
        TString out_dir = "SelectedForPR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        TFile* MuonFile;

        if (Debug == kTRUE)
            MuonFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
        else
            MuonFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
        MuonFile->cd();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        std::vector<int> *isGenMatched = new std::vector<int>;
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
        Int_t runNum;
        Int_t lumiBlock;
        Int_t nElectrons;
        std::vector<double> *ele_pT = new std::vector<double>;
        std::vector<double> *ele_eta = new std::vector<double>;
        std::vector<double> *ele_etaSC = new std::vector<double>;
        std::vector<double> *ele_phi = new std::vector<double>;
        std::vector<int> *ele_charge = new std::vector<int>;
        Int_t nJets;
        std::vector<double> *jet_pT = new std::vector<double>;
        std::vector<double> *jet_eta = new std::vector<double>;
        std::vector<double> *jet_phi = new std::vector<double>;
        std::vector<int> *jet_charge = new std::vector<int>;

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("phi", &phi);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("TRKiso", &TRKiso);
        MuonTree->Branch("isGenMatched", &isGenMatched);
        MuonTree->Branch("trig_name", &trig_name);
        MuonTree->Branch("trig_fired", &trig_fired);
        MuonTree->Branch("trig_matched", &trig_matched);
        MuonTree->Branch("trig_pT", &trig_pT);
        MuonTree->Branch("prescale_factor", &prescale_factor);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("MET_pT", &MET_pT);
        MuonTree->Branch("MET_phi", &MET_phi);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);
        MuonTree->Branch("runNum", &runNum);
        MuonTree->Branch("lumiBlock", &lumiBlock);
        MuonTree->Branch("nElectrons", &nElectrons);
        MuonTree->Branch("ele_pT", &ele_pT);
        MuonTree->Branch("ele_eta", &ele_eta);
        MuonTree->Branch("ele_etaSC", &ele_etaSC);
        MuonTree->Branch("ele_phi", &ele_phi);
        MuonTree->Branch("ele_charge", &ele_charge);
        MuonTree->Branch("nJets", &nJets);
        MuonTree->Branch("jet_pT", &jet_pT);
        MuonTree->Branch("jet_eta", &jet_eta);
        MuonTree->Branch("jet_phi", &jet_phi);
        MuonTree->Branch("jet_charge", &jet_charge);

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Muon();
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_Jet();
            ntuple->TurnOnBranches_MET();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

            std::string RCaddress;
            if (Mgr.Type == "TEST")
                RCaddress = "./etc/RoccoR/rcdata.2016.v3";
            else RCaddress = "/cms/ldap_home/mambroza/DrellYan2016/SelectedX/etc/RoccoR/rcdata.2016.v3";
            RoccoR rc(RCaddress);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                if (!Debug) bar.Draw(i);
                ntuple->GetEvent(i);
                if (ntuple->nMuon < 2) continue;

                if (Debug == kTRUE)
                {
                    cout << "\nEvt " << i << endl;
                    cout << "Nmuons: " << ntuple->nMuon << ",  Ntrig: " << ntuple->HLT_ntrig << endl;
                    for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                    {
                        cout << "\tmu" << i_mu << ": pT=" << ntuple->Muon_pT[i_mu] << "  eta=" << ntuple->Muon_eta[i_mu];
                        Double_t PFiso = (ntuple->Muon_PfChargedHadronIsoR04[i_mu] +
                                 max(0.0, ntuple->Muon_PfNeutralHadronIsoR04[i_mu] + ntuple->Muon_PfGammaIsoR04[i_mu] - 0.5*ntuple->Muon_PFSumPUIsoR04[i_mu]))
                                 / ntuple->Muon_pT[i_mu];
                        cout << "  isTight=" << ntuple->Muon_passTightID[i_mu] << "  PFiso=" << PFiso << endl;
                    }
                }

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                SumWeight += gen_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kTRUE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Debug == kTRUE)
                        cout << "PASS" << endl;
                    // -- Reco level selection -- //
                    vector< Muon > MuonCollection;
                    Int_t NLeptons = ntuple->nMuon;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("finalState_OR_hadronDecay", ntuple, mu);
                            if (genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                            if (mu.Pt != mu.Pt || SF != SF || SF*mu.Pt != SF*mu.Pt)
                                cout << "\nGenPt: " << genPt << "  Pt: " << mu.Pt << "  Corr Pt: " << SF*mu.Pt << "  multiplier: " << SF
                                     << "\nEta: " << mu.eta << "   Phi: " << mu.phi << "  Charge: " << mu.charge << "\ntrLayers: " << mu.trackerLayers
                                     << "  PFiso: " << mu.relPFiso << endl;
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);

                    } // End of i_reco iteration
                    if (Debug == kTRUE) cout << "Muons in collection: " << MuonCollection.size() << endl;

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_PR(MuonCollection, ntuple, &SelectedMuonCollection);

                    if (isPassEventSelection == kTRUE)
                    {
                        if (Debug == kTRUE) cout << "Selection passed" << endl;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
                        relPFiso->clear();
                        TRKiso->clear();
                        isGenMatched->clear();
                        trig_name->clear();
                        trig_fired->clear();
                        trig_matched->clear();
                        trig_pT->clear();
                        prescale_factor->clear();
                        ele_pT->clear();
                        ele_eta->clear();
                        ele_etaSC->clear();
                        ele_phi->clear();
                        ele_charge->clear();
                        jet_pT->clear();
                        jet_eta->clear();
                        jet_phi->clear();
                        jet_charge->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- MET information -- //
                        MET_pT = ntuple->pfMET_pT;
                        MET_phi = ntuple->pfMET_phi;

                        // -- Information for various other reweightings -- //
                        runNum = ntuple->runNum;
                        lumiBlock = ntuple->lumiBlock;
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        Int_t triggered = 0;
                        triggered = analyzer->FindTriggerAndPrescale2(SelectedMuonCollection, ntuple, trig_name, trig_fired, prescale_factor,
                                                                      trig_matched, trig_pT, Debug);
                        if (Debug == kTRUE) cout << "triggered=" << triggered << endl;
                        if (!triggered) continue;
                        timesPassed++;

                        // -- Vector filling -- //
                        for (UInt_t i=0; i<SelectedMuonCollection.size(); i++)
                        {
                            p_T->push_back(SelectedMuonCollection[i].Pt);
                            eta->push_back(SelectedMuonCollection[i].eta);
                            phi->push_back(SelectedMuonCollection[i].phi);
                            charge->push_back(SelectedMuonCollection[i].charge);
                            relPFiso->push_back(SelectedMuonCollection[i].RelPFIso_dBeta);
                            TRKiso->push_back(SelectedMuonCollection[i].trkiso);
                            isGenMatched->push_back(analyzer->isGenMatched(SelectedMuonCollection[i], "fromHardProcess", ntuple));
                        }                        
                        nElectrons = ntuple->Nelectrons;
                        for (Int_t i_ele=0; i_ele<nElectrons; i_ele++)
                        {
                            ele_pT->push_back(ntuple->Electron_pT[i_ele]);
                            ele_eta->push_back(ntuple->Electron_eta[i_ele]);
                            ele_etaSC->push_back(ntuple->Electron_etaSC[i_ele]);
                            ele_phi->push_back(ntuple->Electron_phi[i_ele]);
                            ele_charge->push_back(ntuple->Electron_charge[i_ele]);
                        }
                        nJets = ntuple->Njets;
                        for (Int_t i_jet=0; i_jet<nJets; i_jet++)
                        {
                            jet_pT->push_back(ntuple->Jet_pT[i_jet]);
                            jet_eta->push_back(ntuple->Jet_eta[i_jet]);
                            jet_phi->push_back(ntuple->Jet_phi[i_jet]);
                            jet_charge->push_back(ntuple->Jet_Charge[i_jet]);
                        }
                        MuonTree->Fill();
                    } // End of isPassEvtSelection

                } // End of if(isTriggered)

            } // End of event iteration
            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

            if (Mgr.isMC == kTRUE)
            {
                printf("\tTotal sum of weights: %.1lf\n", SumWeight);
                printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
                printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into files...";
        MuonFile->cd();
        Int_t write;
        write = MuonTree->Write();

        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        if (write)
        {
            cout << " Tree writing finished." << endl << "Closing a file..." << endl;
            MuonFile->Close();
            if (!MuonFile->IsOpen()) cout << "File SelectedForPR_Mu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
            else cout << "FILE SelectedForPR_Mu_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            MuonFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectionForPR_Mu


void CountObjectsInAcceptance (TString type, TString HLTname , Bool_t Debug)
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();
        Double_t n_pass = 0;

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();
            ntuple->TurnOnBranches_Photon();
            ntuple->TurnOnBranches_Muon();
            ntuple->TurnOnBranches_Jet();
            ntuple->TurnOnBranches_MET();

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);
                Int_t n_28=0, n_17=0;

                Double_t gen_weight;
                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;
                if (Processes[i_proc] >= _GJets_20to100 && Processes[i_proc] <= _GJets_2000to5000)
                    gen_weight = ntuple->GENEvt_weight; // Resetting to +-1 affects normalization

                for (Int_t i_ele=0; i_ele<ntuple->Nelectrons; i_ele++)
                {
                    if (ntuple->Electron_pT[i_ele] > 28 && fabs(ntuple->Electron_etaSC[i_ele]) < 2.4 &&
                        (fabs(ntuple->Electron_etaSC[i_ele]) < 1.4442 || fabs(ntuple->Electron_etaSC[i_ele]) > 1.566))
                        n_28++;
                    else if (ntuple->Electron_pT[i_ele] > 17 && fabs(ntuple->Electron_etaSC[i_ele]) < 2.4 &&
                        (fabs(ntuple->Electron_etaSC[i_ele]) < 1.4442 || fabs(ntuple->Electron_etaSC[i_ele]) > 1.566))
                        n_17++;
                }

                for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                {
                    if (ntuple->Muon_pT[i_mu] > 28 && ntuple->Muon_eta[i_mu] < 2.4)
                        n_28++;
                    else if (ntuple->Muon_pT[i_mu] > 17 && ntuple->Muon_eta[i_mu] < 2.4)
                        n_17++;
                }

                for (Int_t i_jet=0; i_jet<ntuple->Njets; i_jet++)
                {
                    if (ntuple->Jet_pT[i_jet] > 28 && ntuple->Jet_eta[i_jet] < 2.4)
                        n_28++;
                    else if (ntuple->Jet_pT[i_jet] > 17 && ntuple->Jet_eta[i_jet] < 2.4)
                        n_17++;
                }

                for (Int_t i_pho=0; i_pho<ntuple->nPhotons; i_pho++)
                {
                    if (ntuple->Photon_pT[i_pho] > 28 && fabs(ntuple->Photon_etaSC[i_pho]) < 2.4 &&
                        (fabs(ntuple->Photon_etaSC[i_pho]) < 1.4442 || fabs(ntuple->Photon_etaSC[i_pho]) > 1.566))
                        n_28++;
                    else if (ntuple->Photon_pT[i_pho] > 17 && fabs(ntuple->Photon_etaSC[i_pho]) < 2.4 &&
                        (fabs(ntuple->Photon_etaSC[i_pho]) < 1.4442 || fabs(ntuple->Photon_etaSC[i_pho]) > 1.566))
                        n_17++;
                }

                if (n_28 > 0 && n_17 > 0)
                    n_pass += gen_weight * Lumi * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup];

                if (!Debug) bar.Draw(i);
            } // End of event iteration


            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        cout << "\t" << n_pass << " events have passed the event selection." << endl;
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of CountObjectsInAcceptance


Int_t skipLumiSection(Int_t runNo, Int_t lumiSection)
{
    Int_t skip = 0;

         if      (runNo == 273158 && (lumiSection == 1   || lumiSection == 1279)) skip = 1;
         else if (runNo == 273302 && (lumiSection == 1   || lumiSection == 459 )) skip = 1;
         else if (runNo == 273402 && (lumiSection == 100 || lumiSection == 292 )) skip = 1;
         else if (runNo == 273403 && (lumiSection == 1   || lumiSection == 53  )) skip = 1;
         else if (runNo == 273404 && (lumiSection == 1   || lumiSection == 18  )) skip = 1;
         else if (runNo == 273405 && (lumiSection == 2   || lumiSection == 25  )) skip = 1;
         else if (runNo == 273406 && (lumiSection == 1   || lumiSection == 112 )) skip = 1;
         else if (runNo == 273408 && (lumiSection == 1   || lumiSection == 6   )) skip = 1;
         else if (runNo == 273409 && (lumiSection == 1   || lumiSection == 309 )) skip = 1;
         else if (runNo == 273410 && (lumiSection == 1   || lumiSection == 90  )) skip = 1;
         else if (runNo == 273411 && (lumiSection == 1   || lumiSection == 29  )) skip = 1;
         else if (runNo == 273425 && (lumiSection == 62  || lumiSection == 733 )) skip = 1;
         else if (runNo == 273446 && (lumiSection == 1   || lumiSection == 33  )) skip = 1;
         else if (runNo == 273447 && (lumiSection == 1   || lumiSection == 412 )) skip = 1;
         else if (runNo == 273448 && (lumiSection == 1   || lumiSection == 391 )) skip = 1;
         else if (runNo == 273449 && (lumiSection == 1   || lumiSection == 214 )) skip = 1;
         else if (runNo == 273450 && (lumiSection == 1   || lumiSection == 647 )) skip = 1;
         else if (runNo == 273492 && (lumiSection == 71  || lumiSection == 338 )) skip = 1;
         else if (runNo == 273493 && (lumiSection == 1   || lumiSection == 233 )) skip = 1;
         else if (runNo == 273494 && (lumiSection == 1   || lumiSection == 192 )) skip = 1;
         else if (runNo == 273502 && (lumiSection == 73  || lumiSection == 1064)) skip = 1;
         else if (runNo == 273503 && (lumiSection == 1   || lumiSection == 598 )) skip = 1;
         else if (runNo == 273554 && (lumiSection == 77  || lumiSection == 437 )) skip = 1;
         else if (runNo == 273555 && (lumiSection == 1   || lumiSection == 173 )) skip = 1;
         else if (runNo == 273725 && (lumiSection == 83  || lumiSection == 2545)) skip = 1;
         else if (runNo == 273728 && (lumiSection == 1   || lumiSection == 100 )) skip = 1;
         else if (runNo == 273730 && (lumiSection == 1   || lumiSection == 2126)) skip = 1;
         else if (runNo == 274094 && (lumiSection == 108 || lumiSection == 332 )) skip = 1;
         else if (runNo == 274146 && (lumiSection == 1   || lumiSection == 67  )) skip = 1;
         else if (runNo == 274157 && (lumiSection == 105 || lumiSection == 534 )) skip = 1;
         else if (runNo == 274159 && (lumiSection == 1   || lumiSection == 43  )) skip = 1;
         else if (runNo == 274160 && (lumiSection == 1   || lumiSection == 207 )) skip = 1;
         else if (runNo == 274161 && (lumiSection == 1   || lumiSection == 516 )) skip = 1;
         else if (runNo == 274172 && (lumiSection == 31  || lumiSection == 95  )) skip = 1;
         else if (runNo == 274198 && (lumiSection == 81  || lumiSection == 191 )) skip = 1;
         else if (runNo == 274199 && (lumiSection == 1   || lumiSection == 623 )) skip = 1;
         else if (runNo == 274200 && (lumiSection == 1   || lumiSection == 678 )) skip = 1;
         else if (runNo == 274240 && (lumiSection == 1   || lumiSection == 82  )) skip = 1;
         else if (runNo == 274241 && (lumiSection == 1   || lumiSection == 1176)) skip = 1;
         else if (runNo == 274244 && (lumiSection == 1   || lumiSection == 607 )) skip = 1;
         else if (runNo == 274250 && (lumiSection == 1   || lumiSection == 701 )) skip = 1;
         else if (runNo == 274251 && (lumiSection == 1   || lumiSection == 546 )) skip = 1;
         else if (runNo == 274283 && (lumiSection == 2   || lumiSection == 19  )) skip = 1;
         else if (runNo == 274284 && (lumiSection == 1   || lumiSection == 210 )) skip = 1;
         else if (runNo == 274286 && (lumiSection == 1   || lumiSection == 154 )) skip = 1;
         else if (runNo == 274314 && (lumiSection == 97  || lumiSection == 158 )) skip = 1;
         else if (runNo == 274315 && (lumiSection == 1   || lumiSection == 424 )) skip = 1;
         else if (runNo == 274316 && (lumiSection == 1   || lumiSection == 959 )) skip = 1;
         else if (runNo == 274317 && (lumiSection == 1   || lumiSection == 3   )) skip = 1;
         else if (runNo == 274319 && (lumiSection == 1   || lumiSection == 225 )) skip = 1;
         else if (runNo == 274335 && (lumiSection == 60  || lumiSection == 1003)) skip = 1;
         else if (runNo == 274336 && (lumiSection == 1   || lumiSection == 14  )) skip = 1;
         else if (runNo == 274337 && (lumiSection == 3   || lumiSection == 17  )) skip = 1;
         else if (runNo == 274338 && (lumiSection == 1   || lumiSection == 698 )) skip = 1;
         else if (runNo == 274339 && (lumiSection == 1   || lumiSection == 93  )) skip = 1;
         else if (runNo == 274344 && (lumiSection == 1   || lumiSection == 632 )) skip = 1;
         else if (runNo == 274345 && (lumiSection == 1   || lumiSection == 170 )) skip = 1;
         else if (runNo == 274382 && (lumiSection == 94  || lumiSection == 144 )) skip = 1;
         else if (runNo == 274387 && (lumiSection == 88  || lumiSection == 439 )) skip = 1;
         else if (runNo == 274388 && (lumiSection == 1   || lumiSection == 1820)) skip = 1;
         else if (runNo == 274420 && (lumiSection == 94  || lumiSection == 268 )) skip = 1;
         else if (runNo == 274421 && (lumiSection == 1   || lumiSection == 342 )) skip = 1;
         else if (runNo == 274422 && (lumiSection == 1   || lumiSection == 2207)) skip = 1;
         else if (runNo == 274440 && (lumiSection == 92  || lumiSection == 493 )) skip = 1;
         else if (runNo == 274441 && (lumiSection == 1   || lumiSection == 431 )) skip = 1;
         else if (runNo == 274442 && (lumiSection == 1   || lumiSection == 752 )) skip = 1;
         else if (runNo == 274954 && (lumiSection == 37  || lumiSection == 57  )) skip = 1;
         else if (runNo == 274955 && (lumiSection == 1   || lumiSection == 91  )) skip = 1;
         else if (runNo == 274968 && (lumiSection == 1   || lumiSection == 1192)) skip = 1;
         else if (runNo == 274969 && (lumiSection == 1   || lumiSection == 1003)) skip = 1;
         else if (runNo == 274970 && (lumiSection == 1   || lumiSection == 47  )) skip = 1;
         else if (runNo == 274971 && (lumiSection == 1   || lumiSection == 905 )) skip = 1;
         else if (runNo == 274998 && (lumiSection == 64  || lumiSection == 782 )) skip = 1;
         else if (runNo == 274999 && (lumiSection == 1   || lumiSection == 1241)) skip = 1;
         else if (runNo == 275000 && (lumiSection == 1   || lumiSection == 136 )) skip = 1;
         else if (runNo == 275001 && (lumiSection == 1   || lumiSection == 2061)) skip = 1;
         else if (runNo == 275059 && (lumiSection == 78  || lumiSection == 137 )) skip = 1;
         else if (runNo == 275066 && (lumiSection == 1   || lumiSection == 96  )) skip = 1;
         else if (runNo == 275067 && (lumiSection == 1   || lumiSection == 392 )) skip = 1;
         else if (runNo == 275068 && (lumiSection == 1   || lumiSection == 915 )) skip = 1;
         else if (runNo == 275073 && (lumiSection == 1   || lumiSection == 517 )) skip = 1;
         else if (runNo == 275074 && (lumiSection == 1   || lumiSection == 647 )) skip = 1;
         else if (runNo == 275124 && (lumiSection == 106 || lumiSection == 431 )) skip = 1;
         else if (runNo == 275125 && (lumiSection == 1   || lumiSection == 989 )) skip = 1;
         else if (runNo == 275282 && (lumiSection == 91  || lumiSection == 180 )) skip = 1;
         else if (runNo == 275283 && (lumiSection == 1   || lumiSection == 132 )) skip = 1;
         else if (runNo == 275284 && (lumiSection == 1   || lumiSection == 74  )) skip = 1;
         else if (runNo == 275290 && (lumiSection == 96  || lumiSection == 143 )) skip = 1;
         else if (runNo == 275291 && (lumiSection == 1   || lumiSection == 347 )) skip = 1;
         else if (runNo == 275292 && (lumiSection == 1   || lumiSection == 121 )) skip = 1;
         else if (runNo == 275293 && (lumiSection == 1   || lumiSection == 201 )) skip = 1;
         else if (runNo == 275309 && (lumiSection == 55  || lumiSection == 617 )) skip = 1;
         else if (runNo == 275310 && (lumiSection == 1   || lumiSection == 1929)) skip = 1;
         else if (runNo == 275311 && (lumiSection == 1   || lumiSection == 1253)) skip = 1;
         else if (runNo == 275319 && (lumiSection == 141 || lumiSection == 282 )) skip = 1;
         else if (runNo == 275337 && (lumiSection == 1   || lumiSection == 427 )) skip = 1;
         else if (runNo == 275338 && (lumiSection == 1   || lumiSection == 520 )) skip = 1;
         else if (runNo == 275344 && (lumiSection == 76  || lumiSection == 356 )) skip = 1;
         else if (runNo == 275345 && (lumiSection == 1   || lumiSection == 353 )) skip = 1;
         else if (runNo == 275370 && (lumiSection == 81  || lumiSection == 365 )) skip = 1;
         else if (runNo == 275371 && (lumiSection == 1   || lumiSection == 569 )) skip = 1;
         else if (runNo == 275375 && (lumiSection == 127 || lumiSection == 1449)) skip = 1;
         else if (runNo == 275376 && (lumiSection == 1   || lumiSection == 3096)) skip = 1;
         else if (runNo == 275657 && (lumiSection == 1   || lumiSection == 105 )) skip = 1;
         else if (runNo == 275658 && (lumiSection == 1   || lumiSection == 337 )) skip = 1;
         else if (runNo == 275659 && (lumiSection == 1   || lumiSection == 17  )) skip = 1;
         else if (runNo == 275761 && (lumiSection == 1   || lumiSection == 9   )) skip = 1;
         else if (runNo == 275767 && (lumiSection == 1   || lumiSection == 4   )) skip = 1;
         else if (runNo == 275772 && (lumiSection == 1   || lumiSection == 56  )) skip = 1;
         else if (runNo == 275773 && (lumiSection == 1   || lumiSection == 7   )) skip = 1;
         else if (runNo == 275774 && (lumiSection == 1   || lumiSection == 315 )) skip = 1;
         else if (runNo == 275776 && (lumiSection == 1   || lumiSection == 140 )) skip = 1;
         else if (runNo == 275777 && (lumiSection == 1   || lumiSection == 300 )) skip = 1;
         else if (runNo == 275778 && (lumiSection == 1   || lumiSection == 305 )) skip = 1;
         else if (runNo == 275782 && (lumiSection == 1   || lumiSection == 762 )) skip = 1;
         else if (runNo == 275832 && (lumiSection == 1   || lumiSection == 367 )) skip = 1;
         else if (runNo == 275833 && (lumiSection == 1   || lumiSection == 251 )) skip = 1;
         else if (runNo == 275834 && (lumiSection == 1   || lumiSection == 297 )) skip = 1;
         else if (runNo == 275835 && (lumiSection == 1   || lumiSection == 13  )) skip = 1;
         else if (runNo == 275836 && (lumiSection == 1   || lumiSection == 1293)) skip = 1;
         else if (runNo == 275837 && (lumiSection == 1   || lumiSection == 726 )) skip = 1;
         else if (runNo == 275847 && (lumiSection == 1   || lumiSection == 2263)) skip = 1;
         else if (runNo == 275886 && (lumiSection == 73  || lumiSection == 109 )) skip = 1;
         else if (runNo == 275890 && (lumiSection == 1   || lumiSection == 1393)) skip = 1;
         else if (runNo == 275911 && (lumiSection == 62  || lumiSection == 440 )) skip = 1;
         else if (runNo == 275912 && (lumiSection == 1   || lumiSection == 289 )) skip = 1;
         else if (runNo == 275913 && (lumiSection == 1   || lumiSection == 475 )) skip = 1;
         else if (runNo == 275918 && (lumiSection == 1   || lumiSection == 361 )) skip = 1;
         else if (runNo == 275920 && (lumiSection == 5   || lumiSection == 463 )) skip = 1;
         else if (runNo == 275921 && (lumiSection == 1   || lumiSection == 20  )) skip = 1;
         else if (runNo == 275923 && (lumiSection == 3   || lumiSection == 126 )) skip = 1;
         else if (runNo == 275931 && (lumiSection == 1   || lumiSection == 89  )) skip = 1;
         else if (runNo == 275963 && (lumiSection == 82  || lumiSection == 172 )) skip = 1;
         else if (runNo == 276092 && (lumiSection == 74  || lumiSection == 149 )) skip = 1;
         else if (runNo == 276097 && (lumiSection == 1   || lumiSection == 507 )) skip = 1;
         else if (runNo == 276242 && (lumiSection == 1   || lumiSection == 1664)) skip = 1;
         else if (runNo == 276243 && (lumiSection == 1   || lumiSection == 611 )) skip = 1;
         else if (runNo == 276244 && (lumiSection == 3   || lumiSection == 1202)) skip = 1;
         else if (runNo == 276282 && (lumiSection == 75  || lumiSection == 1142)) skip = 1;
         else if (runNo == 276283 && (lumiSection == 3   || lumiSection == 1087)) skip = 1;
         else if (runNo == 276315 && (lumiSection == 40  || lumiSection == 217 )) skip = 1;
         else if (runNo == 276317 && (lumiSection == 3   || lumiSection == 138 )) skip = 1;
         else if (runNo == 276318 && (lumiSection == 3   || lumiSection == 570 )) skip = 1;
         else if (runNo == 276355 && (lumiSection == 1   || lumiSection == 33  )) skip = 1;
         else if (runNo == 276361 && (lumiSection == 1   || lumiSection == 833 )) skip = 1;
         else if (runNo == 276363 && (lumiSection == 1   || lumiSection == 1482)) skip = 1;
         else if (runNo == 276384 && (lumiSection == 2   || lumiSection == 1117)) skip = 1;
         else if (runNo == 276437 && (lumiSection == 63  || lumiSection == 2190)) skip = 1;
         else if (runNo == 276454 && (lumiSection == 1   || lumiSection == 527 )) skip = 1;
         else if (runNo == 276458 && (lumiSection == 1   || lumiSection == 341 )) skip = 1;
         else if (runNo == 276495 && (lumiSection == 87  || lumiSection == 268 )) skip = 1;
         else if (runNo == 276501 && (lumiSection == 4   || lumiSection == 2547)) skip = 1;
         else if (runNo == 276502 && (lumiSection == 2   || lumiSection == 741 )) skip = 1;
         else if (runNo == 276525 && (lumiSection == 88  || lumiSection == 2893)) skip = 1;
         else if (runNo == 276527 && (lumiSection == 1   || lumiSection == 214 )) skip = 1;
         else if (runNo == 276528 && (lumiSection == 4   || lumiSection == 394 )) skip = 1;
         else if (runNo == 276542 && (lumiSection == 74  || lumiSection == 857 )) skip = 1;
         else if (runNo == 276543 && (lumiSection == 1   || lumiSection == 952 )) skip = 1;
         else if (runNo == 276544 && (lumiSection == 2   || lumiSection == 161 )) skip = 1;
         else if (runNo == 276545 && (lumiSection == 2   || lumiSection == 213 )) skip = 1;
         else if (runNo == 276581 && (lumiSection == 79  || lumiSection == 444 )) skip = 1;
         else if (runNo == 276582 && (lumiSection == 1   || lumiSection == 871 )) skip = 1;
         else if (runNo == 276583 && (lumiSection == 1   || lumiSection == 52  )) skip = 1;
         else if (runNo == 276584 && (lumiSection == 1   || lumiSection == 2   )) skip = 1;
         else if (runNo == 276585 && (lumiSection == 1   || lumiSection == 246 )) skip = 1;
         else if (runNo == 276586 && (lumiSection == 2   || lumiSection == 773 )) skip = 1;
         else if (runNo == 276587 && (lumiSection == 1   || lumiSection == 1006)) skip = 1;
         else if (runNo == 276653 && (lumiSection == 72  || lumiSection == 550 )) skip = 1;
         else if (runNo == 276655 && (lumiSection == 1   || lumiSection == 1106)) skip = 1;
         else if (runNo == 276659 && (lumiSection == 1   || lumiSection == 252 )) skip = 1;
         else if (runNo == 276775 && (lumiSection == 96  || lumiSection == 1260)) skip = 1;
         else if (runNo == 276776 && (lumiSection == 1   || lumiSection == 1823)) skip = 1;
         else if (runNo == 276794 && (lumiSection == 1   || lumiSection == 885 )) skip = 1;
         else if (runNo == 276807 && (lumiSection == 66  || lumiSection == 220 )) skip = 1;
         else if (runNo == 276808 && (lumiSection == 1   || lumiSection == 875 )) skip = 1;
         else if (runNo == 276810 && (lumiSection == 1   || lumiSection == 287 )) skip = 1;
         else if (runNo == 276811 && (lumiSection == 1   || lumiSection == 2563)) skip = 1;
         else if (runNo == 276831 && (lumiSection == 64  || lumiSection == 2702)) skip = 1;
         else if (runNo == 276834 && (lumiSection == 1   || lumiSection == 720 )) skip = 1;
         else if (runNo == 276870 && (lumiSection == 78  || lumiSection == 3484)) skip = 1;
         else if (runNo == 276935 && (lumiSection == 79  || lumiSection == 906 )) skip = 1;
         else if (runNo == 276940 && (lumiSection == 70  || lumiSection == 213 )) skip = 1;
         else if (runNo == 276946 && (lumiSection == 1   || lumiSection == 27  )) skip = 1;
         else if (runNo == 276947 && (lumiSection == 1   || lumiSection == 141 )) skip = 1;
         else if (runNo == 276948 && (lumiSection == 1   || lumiSection == 474 )) skip = 1;
         else if (runNo == 276950 && (lumiSection == 1   || lumiSection == 2353)) skip = 1;
         else if (runNo == 277069 && (lumiSection == 81  || lumiSection == 390 )) skip = 1;
         else if (runNo == 277070 && (lumiSection == 1   || lumiSection == 1059)) skip = 1;
         else if (runNo == 277071 && (lumiSection == 1   || lumiSection == 178 )) skip = 1;
         else if (runNo == 277072 && (lumiSection == 1   || lumiSection == 466 )) skip = 1;
         else if (runNo == 277073 && (lumiSection == 1   || lumiSection == 90  )) skip = 1;
         else if (runNo == 277076 && (lumiSection == 1   || lumiSection == 1037)) skip = 1;
         else if (runNo == 277087 && (lumiSection == 204 || lumiSection == 1191)) skip = 1;
         else if (runNo == 277094 && (lumiSection == 1   || lumiSection == 584 )) skip = 1;
         else if (runNo == 277096 && (lumiSection == 1   || lumiSection == 2086)) skip = 1;
         else if (runNo == 277112 && (lumiSection == 1   || lumiSection == 155 )) skip = 1;
         else if (runNo == 277126 && (lumiSection == 42  || lumiSection == 59  )) skip = 1;
         else if (runNo == 277127 && (lumiSection == 1   || lumiSection == 902 )) skip = 1;
         else if (runNo == 277148 && (lumiSection == 83  || lumiSection == 700 )) skip = 1;
         else if (runNo == 277166 && (lumiSection == 77  || lumiSection == 431 )) skip = 1;
         else if (runNo == 277168 && (lumiSection == 1   || lumiSection == 2223)) skip = 1;
         else if (runNo == 277180 && (lumiSection == 88  || lumiSection == 228 )) skip = 1;
         else if (runNo == 277194 && (lumiSection == 113 || lumiSection == 2070)) skip = 1;
         else if (runNo == 277305 && (lumiSection == 62  || lumiSection == 744 )) skip = 1;
         else if (runNo == 277420 && (lumiSection == 84  || lumiSection == 346 )) skip = 1;
         else if (runNo == 277981 && (lumiSection == 82  || lumiSection == 163 )) skip = 1;
         else if (runNo == 277991 && (lumiSection == 1   || lumiSection == 98  )) skip = 1;
         else if (runNo == 277992 && (lumiSection == 1   || lumiSection == 312 )) skip = 1;
         else if (runNo == 278017 && (lumiSection == 77  || lumiSection == 589 )) skip = 1;
         else if (runNo == 278018 && (lumiSection == 1   || lumiSection == 1181)) skip = 1;
         else if (runNo == 278167 && (lumiSection == 87  || lumiSection == 2258)) skip = 1;
         else if (runNo == 278175 && (lumiSection == 1   || lumiSection == 88  )) skip = 1;
         else if (runNo == 278193 && (lumiSection == 77  || lumiSection == 231 )) skip = 1;
         else if (runNo == 278239 && (lumiSection == 76  || lumiSection == 740 )) skip = 1;
         else if (runNo == 278240 && (lumiSection == 1   || lumiSection == 1309)) skip = 1;
         else if (runNo == 278273 && (lumiSection == 75  || lumiSection == 110 )) skip = 1;
         else if (runNo == 278274 && (lumiSection == 1   || lumiSection == 85  )) skip = 1;
         else if (runNo == 278288 && (lumiSection == 67  || lumiSection == 81  )) skip = 1;
         else if (runNo == 278289 && (lumiSection == 1   || lumiSection == 52  )) skip = 1;
         else if (runNo == 278290 && (lumiSection == 1   || lumiSection == 11  )) skip = 1;
         else if (runNo == 278308 && (lumiSection == 87  || lumiSection == 1880)) skip = 1;
         else if (runNo == 278310 && (lumiSection == 1   || lumiSection == 709 )) skip = 1;
         else if (runNo == 278315 && (lumiSection == 73  || lumiSection == 767 )) skip = 1;
         else if (runNo == 278345 && (lumiSection == 84  || lumiSection == 831 )) skip = 1;
         else if (runNo == 278346 && (lumiSection == 1   || lumiSection == 117 )) skip = 1;
         else if (runNo == 278349 && (lumiSection == 1   || lumiSection == 633 )) skip = 1;
         else if (runNo == 278366 && (lumiSection == 1   || lumiSection == 453 )) skip = 1;
         else if (runNo == 278406 && (lumiSection == 85  || lumiSection == 1682)) skip = 1;
         else if (runNo == 278509 && (lumiSection == 91  || lumiSection == 1557)) skip = 1;
         else if (runNo == 278769 && (lumiSection == 75  || lumiSection == 104 )) skip = 1;
         else if (runNo == 278770 && (lumiSection == 1   || lumiSection == 767 )) skip = 1;
         else if (runNo == 278801 && (lumiSection == 48  || lumiSection == 85  )) skip = 1;
         else if (runNo == 278802 && (lumiSection == 1   || lumiSection == 17  )) skip = 1;
         else if (runNo == 278803 && (lumiSection == 1   || lumiSection == 323 )) skip = 1;
         else if (runNo == 278804 && (lumiSection == 1   || lumiSection == 4   )) skip = 1;
         else if (runNo == 278805 && (lumiSection == 3   || lumiSection == 288 )) skip = 1;
         else if (runNo == 278808 && (lumiSection == 1   || lumiSection == 1793)) skip = 1;
         else if (runNo == 278820 && (lumiSection == 17  || lumiSection == 1533)) skip = 1;
         else if (runNo == 278822 && (lumiSection == 1   || lumiSection == 1627)) skip = 1;
         else if (runNo == 278873 && (lumiSection == 70  || lumiSection == 129 )) skip = 1;
         else if (runNo == 278874 && (lumiSection == 1   || lumiSection == 478 )) skip = 1;
         else if (runNo == 278875 && (lumiSection == 1   || lumiSection == 834 )) skip = 1;
         else if (runNo == 278923 && (lumiSection == 55  || lumiSection == 467 )) skip = 1;
         else if (runNo == 278957 && (lumiSection == 79  || lumiSection == 227 )) skip = 1;
         else if (runNo == 278962 && (lumiSection == 68  || lumiSection == 408 )) skip = 1;
         else if (runNo == 278963 && (lumiSection == 1   || lumiSection == 175 )) skip = 1;
         else if (runNo == 278969 && (lumiSection == 70  || lumiSection == 1460)) skip = 1;
         else if (runNo == 278975 && (lumiSection == 1   || lumiSection == 850 )) skip = 1;
         else if (runNo == 278976 && (lumiSection == 1   || lumiSection == 20  )) skip = 1;
         else if (runNo == 278986 && (lumiSection == 71  || lumiSection == 199 )) skip = 1;
         else if (runNo == 279024 && (lumiSection == 82  || lumiSection == 382 )) skip = 1;
         else if (runNo == 279029 && (lumiSection == 1   || lumiSection == 260 )) skip = 1;
         else if (runNo == 279071 && (lumiSection == 71  || lumiSection == 244 )) skip = 1;
         else if (runNo == 279080 && (lumiSection == 68  || lumiSection == 224 )) skip = 1;
         else if (runNo == 279115 && (lumiSection == 118 || lumiSection == 524 )) skip = 1;
         else if (runNo == 279116 && (lumiSection == 38  || lumiSection == 485 )) skip = 1;
         else if (runNo == 279479 && (lumiSection == 86  || lumiSection == 190 )) skip = 1;
         else if (runNo == 279588 && (lumiSection == 100 || lumiSection == 1259)) skip = 1;
         else if (runNo == 279653 && (lumiSection == 77  || lumiSection == 261 )) skip = 1;
         else if (runNo == 279654 && (lumiSection == 1   || lumiSection == 1299)) skip = 1;
         else if (runNo == 279656 && (lumiSection == 1   || lumiSection == 43  )) skip = 1;
         else if (runNo == 279658 && (lumiSection == 1   || lumiSection == 713 )) skip = 1;
         else if (runNo == 279667 && (lumiSection == 68  || lumiSection == 1033)) skip = 1;
         else if (runNo == 279681 && (lumiSection == 77  || lumiSection == 104 )) skip = 1;
         else if (runNo == 279682 && (lumiSection == 1   || lumiSection == 38  )) skip = 1;
         else if (runNo == 279683 && (lumiSection == 1   || lumiSection == 26  )) skip = 1;
         else if (runNo == 279684 && (lumiSection == 1   || lumiSection == 22  )) skip = 1;
         else if (runNo == 279685 && (lumiSection == 1   || lumiSection == 209 )) skip = 1;
         else if (runNo == 279691 && (lumiSection == 71  || lumiSection == 113 )) skip = 1;
         else if (runNo == 279694 && (lumiSection == 1   || lumiSection == 2235)) skip = 1;
         else if (runNo == 279715 && (lumiSection == 71  || lumiSection == 691 )) skip = 1;
         else if (runNo == 279716 && (lumiSection == 1   || lumiSection == 1653)) skip = 1;
         else if (runNo == 279760 && (lumiSection == 68  || lumiSection == 728 )) skip = 1;
         else if (runNo == 279766 && (lumiSection == 1   || lumiSection == 1689)) skip = 1;
         else if (runNo == 279767 && (lumiSection == 1   || lumiSection == 776 )) skip = 1;
         else if (runNo == 279794 && (lumiSection == 77  || lumiSection == 1100)) skip = 1;
         else if (runNo == 279823 && (lumiSection == 61  || lumiSection == 395 )) skip = 1;
         else if (runNo == 279841 && (lumiSection == 75  || lumiSection == 2122)) skip = 1;
         else if (runNo == 279844 && (lumiSection == 72  || lumiSection == 295 )) skip = 1;
         else if (runNo == 279887 && (lumiSection == 79  || lumiSection == 397 )) skip = 1;
         else if (runNo == 279931 && (lumiSection == 84  || lumiSection == 3022)) skip = 1;
         else if (runNo == 279966 && (lumiSection == 79  || lumiSection == 441 )) skip = 1;
         else if (runNo == 279975 && (lumiSection == 70  || lumiSection == 1121)) skip = 1;
         else if (runNo == 279993 && (lumiSection == 85  || lumiSection == 156 )) skip = 1;
         else if (runNo == 279994 && (lumiSection == 1   || lumiSection == 47  )) skip = 1;
         else if (runNo == 280013 && (lumiSection == 1   || lumiSection == 25  )) skip = 1;
         else if (runNo == 280015 && (lumiSection == 1   || lumiSection == 580 )) skip = 1;
         else if (runNo == 280016 && (lumiSection == 1   || lumiSection == 149 )) skip = 1;
         else if (runNo == 280017 && (lumiSection == 1   || lumiSection == 608 )) skip = 1;
         else if (runNo == 280018 && (lumiSection == 1   || lumiSection == 1281)) skip = 1;
         else if (runNo == 280020 && (lumiSection == 1   || lumiSection == 45  )) skip = 1;
         else if (runNo == 280024 && (lumiSection == 1   || lumiSection == 427 )) skip = 1;
         else if (runNo == 280187 && (lumiSection == 4   || lumiSection == 60  )) skip = 1;
         else if (runNo == 280188 && (lumiSection == 1   || lumiSection == 245 )) skip = 1;
         else if (runNo == 280191 && (lumiSection == 1   || lumiSection == 900 )) skip = 1;
         else if (runNo == 280194 && (lumiSection == 1   || lumiSection == 238 )) skip = 1;
         else if (runNo == 280242 && (lumiSection == 1   || lumiSection == 627 )) skip = 1;
         else if (runNo == 280249 && (lumiSection == 1   || lumiSection == 1433)) skip = 1;
         else if (runNo == 280251 && (lumiSection == 1   || lumiSection == 372 )) skip = 1;
         else if (runNo == 280327 && (lumiSection == 49  || lumiSection == 85  )) skip = 1;
         else if (runNo == 280330 && (lumiSection == 1   || lumiSection == 857 )) skip = 1;
         else if (runNo == 280349 && (lumiSection == 1   || lumiSection == 626 )) skip = 1;
         else if (runNo == 280363 && (lumiSection == 1   || lumiSection == 359 )) skip = 1;
         else if (runNo == 280364 && (lumiSection == 1   || lumiSection == 1363)) skip = 1;
         else if (runNo == 280383 && (lumiSection == 64  || lumiSection == 65  )) skip = 1;
         else if (runNo == 280384 && (lumiSection == 2   || lumiSection == 34  )) skip = 1;
         else if (runNo == 280385 && (lumiSection == 1   || lumiSection == 2022)) skip = 1;
         else if (runNo == 281613 && (lumiSection == 101 || lumiSection == 903 )) skip = 1;
         else if (runNo == 281639 && (lumiSection == 1   || lumiSection == 132 )) skip = 1;
         else if (runNo == 281641 && (lumiSection == 1   || lumiSection == 319 )) skip = 1;
         else if (runNo == 281693 && (lumiSection == 1   || lumiSection == 2191)) skip = 1;
         else if (runNo == 281707 && (lumiSection == 99  || lumiSection == 1065)) skip = 1;
         else if (runNo == 281726 && (lumiSection == 1   || lumiSection == 288 )) skip = 1;
         else if (runNo == 281727 && (lumiSection == 1   || lumiSection == 1605)) skip = 1;
         else if (runNo == 281797 && (lumiSection == 125 || lumiSection == 2176)) skip = 1;
         else if (runNo == 281975 && (lumiSection == 1   || lumiSection == 215 )) skip = 1;
         else if (runNo == 281976 && (lumiSection == 1   || lumiSection == 2166)) skip = 1;
         else if (runNo == 282033 && (lumiSection == 82  || lumiSection == 117 )) skip = 1;
         else if (runNo == 282034 && (lumiSection == 1   || lumiSection == 33  )) skip = 1;
         else if (runNo == 282035 && (lumiSection == 1   || lumiSection == 40  )) skip = 1;
         else if (runNo == 282037 && (lumiSection == 1   || lumiSection == 1862)) skip = 1;
         else if (runNo == 282092 && (lumiSection == 92  || lumiSection == 2276)) skip = 1;
         else if (runNo == 282708 && (lumiSection == 1   || lumiSection == 8   )) skip = 1;
         else if (runNo == 282710 && (lumiSection == 1   || lumiSection == 8   )) skip = 1;
         else if (runNo == 282712 && (lumiSection == 1   || lumiSection == 68  )) skip = 1;
         else if (runNo == 282730 && (lumiSection == 89  || lumiSection == 164))  skip = 1;
         else if (runNo == 282731 && (lumiSection == 1   || lumiSection == 172 )) skip = 1;
         else if (runNo == 282732 && (lumiSection == 1   || lumiSection == 69  )) skip = 1;
         else if (runNo == 282733 && (lumiSection == 1   || lumiSection == 177 )) skip = 1;
         else if (runNo == 282734 && (lumiSection == 1   || lumiSection == 327 )) skip = 1;
         else if (runNo == 282735 && (lumiSection == 1   || lumiSection == 1823)) skip = 1;
         else if (runNo == 282800 && (lumiSection == 1   || lumiSection == 377 )) skip = 1;
         else if (runNo == 282807 && (lumiSection == 1   || lumiSection == 326 )) skip = 1;
         else if (runNo == 282814 && (lumiSection == 1   || lumiSection == 1843)) skip = 1;
         else if (runNo == 282842 && (lumiSection == 1   || lumiSection == 80  )) skip = 1;
         else if (runNo == 282917 && (lumiSection == 117 || lumiSection == 191))  skip = 1;
         else if (runNo == 282918 && (lumiSection == 1   || lumiSection == 51  )) skip = 1;
         else if (runNo == 282919 && (lumiSection == 1   || lumiSection == 243 )) skip = 1;
         else if (runNo == 282922 && (lumiSection == 1   || lumiSection == 131 )) skip = 1;
         else if (runNo == 282923 && (lumiSection == 1   || lumiSection == 224 )) skip = 1;
         else if (runNo == 283042 && (lumiSection == 1   || lumiSection == 6   )) skip = 1;
         else if (runNo == 283043 && (lumiSection == 1   || lumiSection == 519 )) skip = 1;
         else if (runNo == 283049 && (lumiSection == 82  || lumiSection == 93  )) skip = 1;
         else if (runNo == 283050 && (lumiSection == 1   || lumiSection == 212 )) skip = 1;
         else if (runNo == 283052 && (lumiSection == 1   || lumiSection == 111 )) skip = 1;
         else if (runNo == 283059 && (lumiSection == 1   || lumiSection == 451 )) skip = 1;
         else if (runNo == 283270 && (lumiSection == 76  || lumiSection == 1912)) skip = 1;
         else if (runNo == 283283 && (lumiSection == 4   || lumiSection == 1748)) skip = 1;
         else if (runNo == 283305 && (lumiSection == 79  || lumiSection == 85  )) skip = 1;
         else if (runNo == 283306 && (lumiSection == 1   || lumiSection == 289 )) skip = 1;
         else if (runNo == 283307 && (lumiSection == 1   || lumiSection == 456 )) skip = 1;
         else if (runNo == 283308 && (lumiSection == 1   || lumiSection == 948 )) skip = 1;
         else if (runNo == 283353 && (lumiSection == 80  || lumiSection == 822 )) skip = 1;
         else if (runNo == 283358 && (lumiSection == 1   || lumiSection == 981 )) skip = 1;
         else if (runNo == 283359 && (lumiSection == 1   || lumiSection == 428 )) skip = 1;
         else if (runNo == 283407 && (lumiSection == 82  || lumiSection == 114 )) skip = 1;
         else if (runNo == 283408 && (lumiSection == 1   || lumiSection == 2542)) skip = 1;
         else if (runNo == 283416 && (lumiSection == 49  || lumiSection == 245 )) skip = 1;
         else if (runNo == 283453 && (lumiSection == 83  || lumiSection == 537 )) skip = 1;
         else if (runNo == 283469 && (lumiSection == 74  || lumiSection == 74  )) skip = 1;
         else if (runNo == 283478 && (lumiSection == 76  || lumiSection == 969 )) skip = 1;
         else if (runNo == 283548 && (lumiSection == 145 || lumiSection == 288 )) skip = 1;
         else if (runNo == 283680 && (lumiSection == 1   || lumiSection == 81  )) skip = 1;
         else if (runNo == 283681 && (lumiSection == 1   || lumiSection == 17  )) skip = 1;
         else if (runNo == 283682 && (lumiSection == 1   || lumiSection == 384 )) skip = 1;
         else if (runNo == 283685 && (lumiSection == 1   || lumiSection == 314 )) skip = 1;
         else if (runNo == 283820 && (lumiSection == 67  || lumiSection == 1548)) skip = 1;
         else if (runNo == 283830 && (lumiSection == 1   || lumiSection == 722 )) skip = 1;
         else if (runNo == 283834 && (lumiSection == 1   || lumiSection == 82  )) skip = 1;
         else if (runNo == 283835 && (lumiSection == 1   || lumiSection == 112 )) skip = 1;
         else if (runNo == 283865 && (lumiSection == 1   || lumiSection == 1177)) skip = 1;
         else if (runNo == 283876 && (lumiSection == 65  || lumiSection == 724 )) skip = 1;
         else if (runNo == 283877 && (lumiSection == 1   || lumiSection == 1496)) skip = 1;
         else if (runNo == 283884 && (lumiSection == 349 || lumiSection == 756 )) skip = 1;
         else if (runNo == 283885 && (lumiSection == 1   || lumiSection == 1723)) skip = 1;
         else if (runNo == 283933 && (lumiSection == 88  || lumiSection == 232 )) skip = 1;
         else if (runNo == 283934 && (lumiSection == 1   || lumiSection == 1291)) skip = 1;
         else if (runNo == 283946 && (lumiSection == 85  || lumiSection == 1462)) skip = 1;
         else if (runNo == 283964 && (lumiSection == 1   || lumiSection == 388 )) skip = 1;
         else if (runNo == 284006 && (lumiSection == 73  || lumiSection == 390 )) skip = 1;
         else if (runNo == 284014 && (lumiSection == 1   || lumiSection == 266 )) skip = 1;
         else if (runNo == 284025 && (lumiSection == 110 || lumiSection == 157 )) skip = 1;
         else if (runNo == 284029 && (lumiSection == 1   || lumiSection == 112 )) skip = 1;
         else if (runNo == 284035 && (lumiSection == 1   || lumiSection == 360 )) skip = 1;
         else if (runNo == 284036 && (lumiSection == 1   || lumiSection == 348 )) skip = 1;
         else if (runNo == 284037 && (lumiSection == 1   || lumiSection == 340 )) skip = 1;
         else if (runNo == 284038 && (lumiSection == 1   || lumiSection == 55  )) skip = 1;
         else if (runNo == 284039 && (lumiSection == 1   || lumiSection == 30  )) skip = 1;
         else if (runNo == 284040 && (lumiSection == 1   || lumiSection == 33  )) skip = 1;
         else if (runNo == 284041 && (lumiSection == 1   || lumiSection == 44  )) skip = 1;
         else if (runNo == 284042 && (lumiSection == 1   || lumiSection == 129 )) skip = 1;
         else if (runNo == 284043 && (lumiSection == 1   || lumiSection == 224 )) skip = 1;
         else if (runNo == 284044 && (lumiSection == 1   || lumiSection == 30  )) skip = 1;

    return skip;
}

Int_t skipRun(Int_t runNo) // Option to skip runs which have at least one of the lumisections missing (compared to json file)
{
    Int_t skip = 0;
    // RunB
    if      (runNo == 273158) skip = 1;
    else if (runNo == 273402) skip = 1;
    else if (runNo == 273425) skip = 1;
    else if (runNo == 273449) skip = 1;
    else if (runNo == 273492) skip = 1;
    else if (runNo == 273493) skip = 1;
    else if (runNo == 273554) skip = 1;
    else if (runNo == 273728) skip = 1;
    else if (runNo == 273730) skip = 1;
    else if (runNo == 274146) skip = 1;
    else if (runNo == 274157) skip = 1;
    else if (runNo == 274198) skip = 1;
    else if (runNo == 274200) skip = 1;
    else if (runNo == 274244) skip = 1;
    else if (runNo == 274250) skip = 1;
    else if (runNo == 274314) skip = 1;
    else if (runNo == 274339) skip = 1;
    else if (runNo == 274344) skip = 1;
    else if (runNo == 274387) skip = 1;
    else if (runNo == 274440) skip = 1;
    else if (runNo == 274954) skip = 1;
    else if (runNo == 274969) skip = 1;
    else if (runNo == 275124) skip = 1;
    else if (runNo == 275309) skip = 1;
    else if (runNo == 275319) skip = 1;
    else if (runNo == 275370) skip = 1;
    // RunC
    else if (runNo == 275837) skip = 1;
    else if (runNo == 275911) skip = 1;
    else if (runNo == 276092) skip = 1;
    else if (runNo == 276282) skip = 1;
    // RunD
    else if (runNo == 276437) skip = 1;
    else if (runNo == 276495) skip = 1;
    else if (runNo == 276581) skip = 1;
    else if (runNo == 276775) skip = 1;
    else if (runNo == 275807) skip = 1;
    // RunE
    else if (runNo == 276831) skip = 1;
    else if (runNo == 276834) skip = 1;
    else if (runNo == 276870) skip = 1;
    else if (runNo == 276935) skip = 1;
    else if (runNo == 276940) skip = 1;
    else if (runNo == 276947) skip = 1;
    else if (runNo == 277071) skip = 1;
    else if (runNo == 277076) skip = 1;
    else if (runNo == 277148) skip = 1;
    else if (runNo == 277166) skip = 1;
    else if (runNo == 277168) skip = 1;
    else if (runNo == 277180) skip = 1;
    else if (runNo == 277305) skip = 1;
    else if (runNo == 277420) skip = 1;
    // RunF
    else if (runNo == 277981) skip = 1;
    else if (runNo == 277991) skip = 1;
    else if (runNo == 278017) skip = 1;
    else if (runNo == 278167) skip = 1;
    else if (runNo == 278239) skip = 1;
    else if (runNo == 278273) skip = 1;
    else if (runNo == 278288) skip = 1;
    else if (runNo == 278349) skip = 1;
    else if (runNo == 278406) skip = 1;
    else if (runNo == 278509) skip = 1;
    else if (runNo == 278769) skip = 1;
    else if (runNo == 278801) skip = 1;
    // RunG
    else if (runNo == 278873) skip = 1;
    else if (runNo == 278986) skip = 1;
    else if (runNo == 279080) skip = 1;
    else if (runNo == 279667) skip = 1;
    else if (runNo == 279760) skip = 1;
    else if (runNo == 279766) skip = 1;
    else if (runNo == 279794) skip = 1;
    else if (runNo == 279841) skip = 1;
    else if (runNo == 279887) skip = 1;
    else if (runNo == 279966) skip = 1;
    else if (runNo == 279975) skip = 1;
    else if (runNo == 279993) skip = 1;
    else if (runNo == 280187) skip = 1;
    else if (runNo == 280327) skip = 1;
    else if (runNo == 280349) skip = 1;
    else if (runNo == 280364) skip = 1;
    // RunH
    else if (runNo == 281613) skip = 1;
    else if (runNo == 281641) skip = 1;
    else if (runNo == 281693) skip = 1;
    else if (runNo == 281707) skip = 1;
    else if (runNo == 281727) skip = 1;
    else if (runNo == 281797) skip = 1;
    else if (runNo == 282092) skip = 1;
    else if (runNo == 282710) skip = 1;
    else if (runNo == 282917) skip = 1;
    else if (runNo == 282923) skip = 1;
    else if (runNo == 283049) skip = 1;
    else if (runNo == 283059) skip = 1;
    else if (runNo == 283270) skip = 1;
    else if (runNo == 283283) skip = 1;
    else if (runNo == 283305) skip = 1;
    else if (runNo == 283353) skip = 1;
    else if (runNo == 283359) skip = 1;
    else if (runNo == 283407) skip = 1;
    else if (runNo == 283408) skip = 1;
    else if (runNo == 283416) skip = 1;
    else if (runNo == 283453) skip = 1;
    else if (runNo == 283469) skip = 1; // This has only 1 valid LS and does not exist in SinglePhoton_H file
    else if (runNo == 283478) skip = 1;
    else if (runNo == 283685) skip = 1;
    else if (runNo == 283877) skip = 1;
    else if (runNo == 283933) skip = 1;
    else if (runNo == 283946) skip = 1;
    else if (runNo == 284006) skip = 1;
    else if (runNo == 284044) skip = 1;

    return skip;
}
