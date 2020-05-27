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


// -- Macro for making new data files with only selection-passing events  -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.cc"
#include "./header/FileMgr.h"
#include "./etc/RoccoR/RoccoR.cc"

void MakeSelectionForFR_E (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForFR_Mu (TString type, TString HLTname, Bool_t Debug);

void MakeSelectionForBKGest_EE (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForBKGest_MuMu (TString type, TString HLTname, Bool_t Debug);
void MakeSelectionForBKGest_EMu (TString type, TString HLTname, Bool_t Debug);
void CountObjectsInAcceptance (TString type, TString HLTname, Bool_t Debug);

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
    if (whichX.Contains("EMU") && whichX.Contains("EST"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MakeSelectionForBKGest_EMu (" << type << ", " << HLT << ")  *****" << endl;
        MakeSelectionForBKGest_EMu(type, HLT, Debug);
    }
    else if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("EST"))
        {            
            if (HLTname == "DEFAULT") HLT = "IsoMu24_OR_IsoTkMu24";
            else HLT = HLTname;
            cout << "\n*****  MakeSelectionForBKGest_MuMu (" << type << ", " << HLT << ")  *****" << endl;
            MakeSelectionForBKGest_MuMu(type, HLT, Debug);
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
        else
        {
            Xselected++;
            if (HLTname == "DEFAULT") HLT = "Photon_OR";
            else HLT = HLTname;
            cout << "\n*******      MakeSelectionForFR_E (" << type << ", " << HLT << ")      *******" << endl;
            MakeSelectionForFR_E(type, HLT, Debug);
        }
    }
    else if (whichX.Contains("COUNT"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "None";
        else HLT = HLTname;
        cout << "\n*******      CountObjectsInAcceptance (" << type << ", " << HLT << ")      *******" << endl;
        CountObjectsInAcceptance(type, HLTname, Debug);
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
            out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "SIGNAL")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "BKG")
        {
            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "TEST")
        {
            out_base = "/media/sf_DATA/test/";
            out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];
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

        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        // -- Creating electron variables to assign branches -- //
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
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
        ElectronTree->Branch("prescale_factor", &prescale_factor);

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            std::map<std::pair<int,int>, int> repeats; // to search for repeating events

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
                        for (Int_t i_tr=0; i_tr<ntuple->/*HLT_trigPS->size()*/HLT_ntrig; i_tr++)
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
                        repeats[std::make_pair(ntuple->runNum, ntuple->evtNum)]++;
                        if (repeats[std::make_pair(ntuple->runNum, ntuple->evtNum)] > 1)
                            cout << "Evt " << ntuple->runNum << "; " << ntuple->evtNum << " repeated " <<
                                    repeats[std::make_pair(ntuple->runNum, ntuple->evtNum)] << " times." << endl;

                        if (Debug == kTRUE) cout << "\nEvent " << i << endl << triggername << endl;
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
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

//                        if (analyzer->HLT == "HLT_Photon*")
//                            prescale_factor = analyzer->PrescaleFactor(SelectedElectronCollection, ntuple, trig_fired, trig_matched, trig_pT);
//                        if (analyzer->HLT == "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")
//                            prescale_factor = analyzer->PrescaleFactor2(SelectedElectronCollection, ntuple, trig_fired, trig_matched, trig_pT);
//                        if (analyzer->HLT == "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v* && HLT_Photon*")
//                            prescale_factor = analyzer->PrescaleFactor3(SelectedElectronCollection, ntuple, trig_fired, trig_matched, trig_pT);

//                        if (prescale_factor < 0) continue; // If no trigger match between selected electrons

                        Int_t triggered = 0;
                        triggered = analyzer->FindTriggerAndPrescale(SelectedElectronCollection, ntuple, trig_fired, prescale_factor, trig_matched, trig_pT);

                        // -- Vector filling -- //
                        for (UInt_t i_ele=0; i_ele<SelectedElectronCollection.size(); i_ele++)
                        {
                            p_T->push_back(SelectedElectronCollection[i_ele].Pt);
                            eta->push_back(SelectedElectronCollection[i_ele].etaSC);
                            phi->push_back(SelectedElectronCollection[i_ele].phi);
                            charge->push_back(SelectedElectronCollection[i_ele].charge);
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

            Int_t noreps, reps;
            for (std::map<std::pair<int,int>,int>::iterator it=repeats.begin(); it!=repeats.end(); it++)
            {
                if (it->second == 1)noreps++;
                else if (it->second == 0) cout << "0 repetitions.. What??" << endl;
                else reps++;
            }
            cout << noreps << " events without repetitions\n" << reps << " events with repetitions" << endl;

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
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating muon variables to assign branches -- //
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("phi", &phi);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("TRKiso", &TRKiso);
        MuonTree->Branch("MET_pT", &MET_pT);
        MuonTree->Branch("MET_phi", &MET_phi);
        MuonTree->Branch("MET_sumEt", &MET_sumEt);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);

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
                        MET_sumEt = ntuple->pfMET_SumEt;

                        // -- Information for various other reweightings -- //
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

                        // -- Vector filling -- //
                        for (UInt_t i=0; i<SelectedMuonCollection_deno.size(); i++)
                        {
                            p_T->push_back(SelectedMuonCollection_deno[i].Pt);
                            eta->push_back(SelectedMuonCollection_deno[i].eta);
                            phi->push_back(SelectedMuonCollection_deno[i].phi);
                            charge->push_back(SelectedMuonCollection_deno[i].charge);
                            relPFiso->push_back(SelectedMuonCollection_deno[i].RelPFIso_dBeta);
                            TRKiso->push_back(SelectedMuonCollection_deno[i].trkiso);
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

        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("etaSC", &etaSC);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
        ElectronTree->Branch("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        ElectronTree->Branch("dEtaInSeed", &dEtaInSeed);
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

/// Muon QCD+WJET estimation WITHOUT TRIGGER
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
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("phi", &phi);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("TRKiso", &TRKiso);
        MuonTree->Branch("nPU", &nPU);
        MuonTree->Branch("nVTX", &nVTX);
        MuonTree->Branch("PVz", &PVz);
        MuonTree->Branch("gen_weight", &gen_weight);
        MuonTree->Branch("top_weight", &top_weight);
        MuonTree->Branch("prefiring_weight", &prefiring_weight);
        MuonTree->Branch("prefiring_weight_up", &prefiring_weight_up);
        MuonTree->Branch("prefiring_weight_down", &prefiring_weight_down);

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

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Get GenTopCollection -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += gen_weight;

                Bool_t TriggerFlag = kTRUE; // We don't use trigger here
//                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

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
                    vector< Muon > SelectedMuonCollection;
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FakeMuons_Triggerless(MuonCollection, ntuple, &SelectedMuonCollection);

                    if (isPassEventSelection == kTRUE)
                    {
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        phi->clear();
                        charge->clear();
                        relPFiso->clear();
                        TRKiso->clear();

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // -- Information for various other reweightings -- //
                        nPU = ntuple->nPileUp;
                        nVTX = ntuple->nVertices;
                        PVz = ntuple->PVz;
                        prefiring_weight = ntuple->_prefiringweight;
                        prefiring_weight_up = ntuple->_prefiringweightup;
                        prefiring_weight_down = ntuple->_prefiringweightdown;

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
        Double_t MET_pT, MET_phi;
        Int_t nPU, nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

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

                Bool_t TriggerFlag = kTRUE; // We don't use trigger here

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

                        // -- Top pT reweighting -- //
                        top_weight = 1;
                        if (Mgr.Tag[i_tup].Contains("ttbar"))
                        {
                            Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                            Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                            top_weight = sqrt(SF0 * SF1);
                        }

                        // Information for reweightings
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

            Double_t n_pass = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);
            Int_t timesPassed = 0;

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
                    if (ntuple->Electron_pT[i_ele] > 17 && fabs(ntuple->Electron_etaSC[i_ele]) < 2.4 &&
                        (fabs(ntuple->Electron_etaSC[i_ele]) < 1.4442 || fabs(ntuple->Electron_etaSC[i_ele]) > 1.566))
                        n_17++;
                }

                for (Int_t i_mu=0; i_mu<ntuple->nMuon; i_mu++)
                {
                    if (ntuple->Muon_pT[i_mu] > 28 && ntuple->Muon_eta[i_mu] < 2.4)
                        n_28++;
                    if (ntuple->Muon_pT[i_mu] > 17 && ntuple->Muon_eta[i_mu] < 2.4)
                        n_17++;
                }

                for (Int_t i_jet=0; i_jet<ntuple->Njets; i_jet++)
                {
                    if (ntuple->Jet_pT[i_jet] > 28 && ntuple->Jet_eta[i_jet] < 2.4)
                        n_28++;
                    if (ntuple->Jet_pT[i_jet] > 17 && ntuple->Jet_eta[i_jet] < 2.4)
                        n_17++;
                }

                for (Int_t i_pho=0; i_pho<ntuple->nPhotons; i_pho++)
                {
                    if (ntuple->Photon_pT[i_ele] > 28 && fabs(ntuple->Photon_etaSC[i_ele]) < 2.4 &&
                        (fabs(ntuple->Photon_etaSC[i_ele]) < 1.4442 || fabs(ntuple->Photon_etaSC[i_ele]) > 1.566))
                        n_28++;
                    if (ntuple->Photon_pT[i_ele] > 17 && fabs(ntuple->Photon_etaSC[i_ele]) < 2.4 &&
                        (fabs(ntuple->Photon_etaSC[i_ele]) < 1.4442 || fabs(ntuple->Photon_etaSC[i_ele]) > 1.566))
                        n_17++;
                }

                if (n_28 > 0 && n_17 > 0)
                    n_pass += gen_weight * Lumi * Mgr.Xsec / Mgr.Wsum;

                if (!Debug) bar.Draw(i);
            } // End of event iteration

            cout << "\t" << n_pass << " events have passed the event selection." << endl;

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of CountObjectsInAcceptance
