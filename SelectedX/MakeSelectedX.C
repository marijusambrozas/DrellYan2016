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

void MakeSelectedEE (TString type, TString HLTname, Bool_t Debug);
void MakeSelectedMuMu (TString type, TString HLTname, Bool_t Debug);
void MakeSelectedEMu (TString type, TString HLTname, Bool_t Debug);

void MakeSelectedQCDEM_120to170 (TString HLTname, Int_t name, Bool_t Debug);
void MakeSelectedQCDEM_120to170_merged();


void MakeSelectedX (TString WhichX, TString type = "", TString HLTname = "DEFAULT")
{
    TString whichX = WhichX;
    whichX.ToUpper();
    TString HLT;
    Int_t Xselected = 0;
    Bool_t Debug = kFALSE;

    // for non-interactive status output
    ofstream fs;
    fs.open("./Status/Status_"+whichX+"_"+type+".txt");
    fs << "Process MakeSelectedX (" << whichX << ", " << type << ", " << HLTname << ") initiated.\n";

    if (whichX.Contains("DEBUG"))
    {
        Debug = kTRUE;
        cout << "DEBUG MODE ON. Running on 100 events only" << endl;
    }
    if (whichX.Contains("EE"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******      MakeSelectedEE (" << type << ", " << HLT << ")      *******" << endl;
        MakeSelectedEE(type, HLT, Debug);
    }
    if (whichX.Contains("MUMU"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MakeSelectedMuMu (" << type << ", " << HLT << ")  *****" << endl;
        MakeSelectedMuMu(type, HLT, Debug);
    }
    if (whichX.Contains("EMU"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   MakeSelectedEMu (" << type << ", " << HLT << ")  *****" << endl;
        MakeSelectedEMu(type, HLT, Debug);
    }
    if (whichX.Contains("QCDFAIL"))   // To run through QCDEMEnriched_pT120to170 file that crashes
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
        else HLT = HLTname;      
        for (Int_t name = 1; name <= 316; name++)
        {
            if (name == 116) continue;
            cout << "\n** MakeSelectedQCDEM_120to170 (skim_" << name << ", " << HLT << ") **" << endl;
            MakeSelectedQCDEM_120to170(HLT, name, Debug);
        }
    }
    if (whichX.Contains("QCDMERGE")) // to merge 316-1(that fails) selected QCD files
    {
        Xselected++;
        cout << "\n****   MakeSelectedQCDEM_120to170()   ****" << endl;
        MakeSelectedQCDEM_120to170_merged();
    }

    if (Xselected == 0) {
        cout << "Wrong arument!" << endl;
        fs << "Process MakeSelectedX (" << whichX << ", " << type << ", " << HLTname << ") FAILED. Wrong arguments!\n";
    }
    else fs << "Process MakeSelectedX (" << whichX << ", " << type << ", " << HLTname << ") FINISHED.\n";

    fs.close();

} // End of MakeSelectedX()


/// ----------------------------- Electron Channel ------------------------------ ///
void MakeSelectedEE (TString type, TString HLTname , Bool_t Debug)
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
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {           
//            if (Mgr.CurrentProc == _WJets) i_tup = 2;

            TStopwatch looptime;
            looptime.Start();

            if (Mgr.Tag[i_tup] == "QCDEMEnriched_Pt120to170") continue; // One of the skims crashes here

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            //Creating a file
            TString out_base;
            TString out_dir;
            TFile* ElectronFile;
            if (Mgr.Type == "DATA")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/";
//                out_dir = "Data/SelectedEE_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedEE_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "SIGNAL")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/";
//                out_dir = "MC_signal/SelectedEE_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedEE_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "BKG")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/";
//                out_dir = "MC_bkg/SelectedEE_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedEE_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "TEST")
            {
                out_base = "/media/sf_DATA/test/";
                out_dir = "SelectedEE_"+Mgr.Tag[i_tup];
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


            TTree* ElectronTree = new TTree("DYTree", "DYTree");
            // -- Creating LongSelectedEE variables to assign branches -- //
            SelectedEE_t EE; EE.CreateNew();
            EE.MakeBranches(ElectronTree);

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            if (Mgr.isMC == kTRUE)
            {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
            }
            ntuple->TurnOnBranches_Electron();

            Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            Int_t timesPassed = 0;           
            Int_t isClear = 0; // vectors that are writen into files should be cleared afer each event

            myProgressBar_t bar(NEvents);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? EE.GENEvt_weight = -1 : EE.GENEvt_weight = 1;
                SumWeight += EE.GENEvt_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Separate ttbar samples -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += EE.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Mgr.Tag[i_tup].Contains("ttbar"))
                    {
                        // -- Top pT reweighting -- //
                        Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                        Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                        EE._topPtWeight = sqrt(SF0 * SF1);
                    }
                    else EE._topPtWeight = 1;

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
                    EE.isSelPassed = 0;
                    isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection);

                    if (isPassEventSelection == kTRUE)
                    {                       
                        timesPassed++;
                        Electron ele1 = SelectedElectronCollection[0];
                        Electron ele2 = SelectedElectronCollection[1];

                        EE.isSelPassed = 1;
                        EE.nVertices = ntuple->nVertices;
                        EE.nPileUp = ntuple->nPileUp;
                        EE._prefiringweight = ntuple->_prefiringweight;
                        EE._prefiringweightup = ntuple->_prefiringweightup;
                        EE._prefiringweightdown = ntuple->_prefiringweightdown;
                        EE.PVz = ntuple->PVz;
                        EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

                        EE.Electron_pT->push_back(ele1.Pt);
                        EE.Electron_pT->push_back(ele2.Pt);
                        EE.Electron_eta->push_back(ele1.eta);
                        EE.Electron_eta->push_back(ele2.eta);
                        EE.Electron_phi->push_back(ele1.phi);
                        EE.Electron_phi->push_back(ele2.phi);
                        EE.Electron_Energy->push_back(ele1.Energy);
                        EE.Electron_Energy->push_back(ele2.Energy);
                        EE.Electron_charge->push_back(ele1.charge);
                        EE.Electron_charge->push_back(ele2.charge);
                        EE.Electron_etaSC->push_back(ele1.etaSC);
                        EE.Electron_etaSC->push_back(ele2.etaSC);
                        EE.Electron_phiSC->push_back(ele1.phiSC);
                        EE.Electron_phiSC->push_back(ele2.phiSC);
                        EE.Electron_Energy_uncorr->push_back(ele1.Energy_uncorr);
                        EE.Electron_Energy_uncorr->push_back(ele2.Energy_uncorr);

                        ElectronTree->Fill();
                        isClear = EE.ClearVectors();
                        if (!isClear)
                        {
                            cout << "======== ERROR: The vectors were not cleared ========" << endl;
                            break;
                        }

                    } // End of event selection

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

            // Writing
            cout << "Writing into file...";
            Int_t write;
            write = ElectronTree->Write();
            if (write)
            {
                cout << " Finished." << endl << "Closing a file..." << endl;
                TString addition = "";
                if (Debug == kTRUE) addition = "_DEBUG";
                ElectronFile->Close();
                if (!ElectronFile->IsOpen()) cout << "File SelectedEE_" << Mgr.Tag[i_tup]+addition << ".root has been closed successfully.\n" << endl;
                else cout << "FILE SelectedEE_" << Mgr.Tag[i_tup]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
            }
            else
            {
                cout << " Writing was NOT successful!\n" << endl;
                ElectronFile->Close();
            }

        } // End of i_tup iteration

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectedEE


/// -------------------------------- Muon Channel ------------------------------------ ///
void MakeSelectedMuMu (TString type, TString HLTname, Bool_t Debug)
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
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
//            if (Mgr.CurrentProc == _WJets) i_tup = 2;

            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            //Creating a file
            TString out_base;
            TString out_dir;
            TFile* MuonFile;
            if (Mgr.Type == "DATA")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "Data/SelectedMuMu_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedMuMu_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "SIGNAL")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_signal/SelectedMuMu_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedMuMu_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "BKG")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_bkg/SelectedMuMu_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedMuMu_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "TEST")
            {
                out_base = "/media/sf_DATA/test/";
                out_dir = "SelectedMuMu_"+Mgr.Tag[i_tup];
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


            TTree* MuonTree = new TTree("DYTree", "DYTree");
            // -- Creating SelectedMuMu variables to assign branches -- //
            SelectedMuMu_t MuMu; MuMu.CreateNew();
            MuMu.MakeBranches(MuonTree);

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
            Int_t isClear = 0; // vectors that are writen into files should be cleared afer each event

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
                ntuple->GENEvt_weight < 0 ? MuMu.GENEvt_weight = -1 : MuMu.GENEvt_weight = 1;
                SumWeight += MuMu.GENEvt_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Separate ttbar samples -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += MuMu.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Mgr.Tag[i_tup].Contains("ttbar"))
                    {
                        // -- Top pT reweighting -- //
                        Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                        Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                        MuMu._topPtWeight = sqrt(SF0 * SF1);
                    }
                    else MuMu._topPtWeight = 1;

                    // -- Reco level selection -- //
                    vector< Muon > MuonCollection;
                    vector< Muon > MuonCollection_noRocCorr;
                    Int_t NLeptons = ntuple->nMuon;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

                        MuonCollection_noRocCorr.push_back(mu);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("fromHardProcessFinalState", ntuple, mu);
                            if (genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);

                    } // End of i_reco iteration

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection, SelectedMuonCollection_noRocCorr;
                    Bool_t isPassEventSelection = kFALSE;
                    Bool_t isPassEventSelection_noRocCorr = kFALSE;
                    MuMu.isSelPassed = 0;
                    isPassEventSelection = analyzer->EventSelection(MuonCollection, ntuple, &SelectedMuonCollection);
                    isPassEventSelection_noRocCorr = analyzer->EventSelection(MuonCollection_noRocCorr, ntuple, &SelectedMuonCollection_noRocCorr);

                    if (isPassEventSelection == kTRUE || isPassEventSelection_noRocCorr == kTRUE)
                    {                     
                        if (isPassEventSelection == kTRUE && isPassEventSelection_noRocCorr == kTRUE)
                        {
                            timesPassed++;
                            MuMu.isSelPassed = 1;

                            Muon mu1 = SelectedMuonCollection[0];
                            Muon mu2 = SelectedMuonCollection[1];
                            Muon mu1_noRC = SelectedMuonCollection_noRocCorr[0];
                            Muon mu2_noRC = SelectedMuonCollection_noRocCorr[1];

                            MuMu.Muon_pT->push_back(mu1.Pt);
                            MuMu.Muon_pT->push_back(mu2.Pt);
                            MuMu.Muon_Energy->push_back(mu1.Momentum.E());
                            MuMu.Muon_Energy->push_back(mu2.Momentum.E());
                            MuMu.Muon_eta->push_back(mu1.eta);
                            MuMu.Muon_eta->push_back(mu2.eta);
                            MuMu.Muon_phi->push_back(mu1.phi);
                            MuMu.Muon_phi->push_back(mu2.phi);
                            MuMu.Muon_charge->push_back(mu1.charge);
                            MuMu.Muon_charge->push_back(mu2.charge);
                            MuMu.Muon_TuneP_pT->push_back(mu1.TuneP_pT);
                            MuMu.Muon_TuneP_pT->push_back(mu2.TuneP_pT);
                            MuMu.Muon_pT_uncorr->push_back(mu1_noRC.Pt);
                            MuMu.Muon_pT_uncorr->push_back(mu2_noRC.Pt);
                            MuMu.Muon_eta_uncorr->push_back(mu1_noRC.eta);
                            MuMu.Muon_eta_uncorr->push_back(mu2_noRC.eta);
                            MuMu.Muon_phi_uncorr->push_back(mu1_noRC.phi);
                            MuMu.Muon_phi_uncorr->push_back(mu2_noRC.phi);
                            MuMu.Muon_charge_uncorr->push_back(mu1_noRC.charge);
                            MuMu.Muon_charge_uncorr->push_back(mu2_noRC.charge);

                            MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();
                        }
                        else if (isPassEventSelection == kTRUE && isPassEventSelection_noRocCorr == kFALSE)
                        {
                            timesPassed++;
                            MuMu.isSelPassed = 2;

                            Muon mu1 = SelectedMuonCollection[0];
                            Muon mu2 = SelectedMuonCollection[1];

                            MuMu.Muon_pT->push_back(mu1.Pt);
                            MuMu.Muon_pT->push_back(mu2.Pt);
                            MuMu.Muon_Energy->push_back(mu1.Momentum.E());
                            MuMu.Muon_Energy->push_back(mu2.Momentum.E());
                            MuMu.Muon_eta->push_back(mu1.eta);
                            MuMu.Muon_eta->push_back(mu2.eta);
                            MuMu.Muon_phi->push_back(mu1.phi);
                            MuMu.Muon_phi->push_back(mu2.phi);
                            MuMu.Muon_charge->push_back(mu1.charge);
                            MuMu.Muon_charge->push_back(mu2.charge);
                            MuMu.Muon_TuneP_pT->push_back(mu1.TuneP_pT);
                            MuMu.Muon_TuneP_pT->push_back(mu2.TuneP_pT);

                            MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();
                        }
                        else if (isPassEventSelection == kFALSE && isPassEventSelection_noRocCorr == kTRUE)
                        {
                            MuMu.isSelPassed = 3;
                            Muon mu1_noRC = SelectedMuonCollection_noRocCorr[0];
                            Muon mu2_noRC = SelectedMuonCollection_noRocCorr[1];

                            MuMu.Muon_pT_uncorr->push_back(mu1_noRC.Pt);
                            MuMu.Muon_pT_uncorr->push_back(mu2_noRC.Pt);
                            MuMu.Muon_eta_uncorr->push_back(mu1_noRC.eta);
                            MuMu.Muon_eta_uncorr->push_back(mu2_noRC.eta);
                            MuMu.Muon_phi_uncorr->push_back(mu1_noRC.phi);
                            MuMu.Muon_phi_uncorr->push_back(mu2_noRC.phi);
                            MuMu.Muon_charge_uncorr->push_back(mu1_noRC.charge);
                            MuMu.Muon_charge_uncorr->push_back(mu2_noRC.charge);
                        }

                        MuMu.nVertices = ntuple->nVertices;
                        MuMu.nPileUp = ntuple->nPileUp;
                        MuMu._prefiringweight = ntuple->_prefiringweight;
                        MuMu._prefiringweightup = ntuple->_prefiringweightup;
                        MuMu._prefiringweightdown = ntuple->_prefiringweightdown;
                        MuMu.PVz = ntuple->PVz;

                        MuonTree->Fill();
                        isClear = MuMu.ClearVectors();
                        if (!isClear)
                        {
                            cout << "======== ERROR: The vectors were not cleared ========" << endl;
                            break;
                        }
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
                printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.nEvents[i_tup]);
            }

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

            // Writing
            cout << "Writing into file...";
            Int_t write;
            write = MuonTree->Write();
            if (write)
            {
                cout << " Finished." << endl << "Closing a file..." << endl;
                TString addition = "";
                if (Debug == kTRUE) addition = "_DEBUG";
                MuonFile->Close();
                if (!MuonFile->IsOpen()) cout << "File SelectedMuMu_" << Mgr.Tag[i_tup]+addition << ".root has been closed successfully.\n" << endl;
                else cout << "FILE SelectedMuMu_" << Mgr.Tag[i_tup]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
            }
            else
            {
                cout << " Writing was NOT successful!\n" << endl;
                MuonFile->Close();
            }

        } // End of i_tup iteration

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectedMuMu


/// --------------------------------- EMu events --------------------------------- ///
void MakeSelectedEMu (TString type, TString HLTname, Bool_t Debug)
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
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
//            if (Mgr.CurrentProc == _WJets) i_tup = 2;

            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            //Creating a file
            TString out_base;
            TString out_dir;
            TFile* EMuFile;
            if (Mgr.Type == "DATA")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEMu/";
//                out_dir = "Data/SelectedEMu_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedEMu_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "BKG")
            {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEMu/";
//                out_dir = "MC_bkg/SelectedEMu_"+Mgr.Tag[i_tup];

                out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
                out_dir = "SelectedEMu_"+Mgr.Tag[i_tup];
            }
            else if (Mgr.Type == "TEST")
            {
                out_base = "/media/sf_DATA/test/";
                out_dir = "SelectedEMu_"+Mgr.Tag[i_tup];
            }
            else
            {
                cout << "Problems with TYPE." << endl;
                return;
            }

            if (Debug == kTRUE)
                EMuFile = TFile::Open(out_base+out_dir+"_DEBUG.root", "RECREATE");
            else
                EMuFile = TFile::Open(out_base+out_dir+".root", "RECREATE");
            EMuFile->cd();


            TTree* EMuTree = new TTree("DYTree", "DYTree");
            // -- Creating SelectedMuMu variables to assign branches -- //
            SelectedEMu_t EMu; EMu.CreateNew();
            EMu.MakeBranches(EMuTree);

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

            RoccoR rc("./etc/RoccoR/rcdata.2016.v3");

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                ntuple->GENEvt_weight < 0 ? EMu.GENEvt_weight = -1 : EMu.GENEvt_weight = 1;
                SumWeight += EMu.GENEvt_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kFALSE;
                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Separate ttbar samples -- //
                Bool_t GenFlag_top = kFALSE;
                vector<GenOthers> GenTopCollection;
                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += EMu.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
                    if (Mgr.Tag[i_tup].Contains("ttbar"))
                    {
                        // -- Top pT reweighting -- //
                        Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
                        Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
                        EMu._topPtWeight = sqrt(SF0 * SF1);
                    }
                    else EMu._topPtWeight = 1;

                    // -- Reco level selection -- //
                    // Muons
                    vector<Muon> MuonCollection;
                    vector<Muon> MuonCollection_noRocCorr;
                    Int_t N_Muons = ntuple->nMuon;
                    for(Int_t i_reco=0; i_reco<N_Muons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

                        MuonCollection_noRocCorr.push_back(mu);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.isMC == kFALSE)
                            SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
                        else
                        {
                            Double_t genPt = analyzer->GenMuonPt("fromHardProcessFinalState", ntuple, mu);
                            if(genPt > 0)
                                SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
                            else
                                SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
                        }
                        mu.Pt = SF*mu.Pt;
                        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);

                        MuonCollection.push_back(mu);
                    }// End of muon i_reco iteration

                    // Electrons
                    vector< Electron > ElectronCollection;
                    Int_t N_Electrons = ntuple->Nelectrons;
                    for (Int_t i_reco=0; i_reco<N_Electrons; i_reco++)
                    {
                            Electron ele;
                            ele.FillFromNtuple(ntuple, i_reco);

                            ElectronCollection.push_back(ele);
                    }

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection, SelectedMuonCollection_noRocCorr;
                    vector< Electron > SelectedElectronCollection, SelectedElectronCollection_noRocCorr;
                    Bool_t isPassEventSelection = kFALSE;
                    Bool_t isPassEventSelection_noRocCorr = kFALSE;
                    EMu.isSelPassed = 0;
                    isPassEventSelection = analyzer->EventSelection_emu_method(MuonCollection, ElectronCollection, ntuple, &SelectedMuonCollection,
                                                                               &SelectedElectronCollection);
                    isPassEventSelection_noRocCorr = analyzer->EventSelection_emu_method(MuonCollection_noRocCorr, ElectronCollection, ntuple,
                                                                                         &SelectedMuonCollection_noRocCorr, &SelectedElectronCollection_noRocCorr);

                    if (isPassEventSelection == kTRUE || isPassEventSelection_noRocCorr == kTRUE)
                    {
                        if (isPassEventSelection == kTRUE && isPassEventSelection_noRocCorr == kTRUE)
                        {
                            timesPassed++;
                            EMu.isSelPassed = 1;

                            Muon mu = SelectedMuonCollection[0];
                            Muon mu_noRC = SelectedMuonCollection_noRocCorr[0];
                            Electron ele = SelectedElectronCollection[0];
                            Electron ele_noRC = SelectedElectronCollection_noRocCorr[0];

                            EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();

                            EMu.Muon_pT = mu.Pt;
                            EMu.Muon_eta = mu.eta;
                            EMu.Muon_phi = mu.phi;
                            EMu.Muon_Energy = mu.Momentum.E();
                            EMu.Muon_charge = mu.charge;
                            EMu.Muon_TuneP_pT = mu.TuneP_pT;

                            EMu.Muon_pT_uncorr = mu_noRC.Pt;
                            EMu.Muon_eta_uncorr = mu_noRC.eta;
                            EMu.Muon_phi_uncorr = mu_noRC.phi;
                            EMu.Muon_charge_uncorr = mu_noRC.charge;

                            EMu.Electron_pT = ele.Pt;
                            EMu.Electron_eta = ele.eta;
                            EMu.Electron_phi = ele.phi;
                            EMu.Electron_Energy = ele.Energy;
                            EMu.Electron_charge = ele.charge;
                            EMu.Electron_etaSC = ele.etaSC;
                            EMu.Electron_phiSC = ele.phiSC;

                            EMu.Electron_Energy_uncorr = ele.Energy_uncorr;
                        }
                        else if (isPassEventSelection == kTRUE && isPassEventSelection_noRocCorr == kFALSE)
                        {
                            timesPassed++;
                            EMu.isSelPassed = 2;

                            Muon mu = SelectedMuonCollection[0];
                            Electron ele = SelectedElectronCollection[0];

                            EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();

                            EMu.Muon_pT = mu.Pt;
                            EMu.Muon_eta = mu.eta;
                            EMu.Muon_phi = mu.phi;
                            EMu.Muon_Energy = mu.Momentum.E();
                            EMu.Muon_charge = mu.charge;
                            EMu.Muon_TuneP_pT = mu.TuneP_pT;

                            EMu.Muon_pT_uncorr = -9999;
                            EMu.Muon_eta_uncorr = -9999;
                            EMu.Muon_phi_uncorr = -9999;
                            EMu.Muon_charge_uncorr = -9999;

                            EMu.Electron_pT = ele.Pt;
                            EMu.Electron_eta = ele.eta;
                            EMu.Electron_phi = ele.phi;
                            EMu.Electron_Energy = ele.Energy;
                            EMu.Electron_charge = ele.charge;
                            EMu.Electron_etaSC = ele.etaSC;
                            EMu.Electron_phiSC = ele.phiSC;

                            EMu.Electron_Energy_uncorr = ele.Energy_uncorr;
                        }
                        else if (isPassEventSelection == kFALSE && isPassEventSelection_noRocCorr == kTRUE)
                        {
                            EMu.isSelPassed = 3;

                            Muon mu_noRC = SelectedMuonCollection_noRocCorr[0];
                            Electron ele_noRC = SelectedElectronCollection_noRocCorr[0];

                            EMu.EMu_InvM = (mu_noRC.Momentum + ele_noRC.Momentum).M();

                            EMu.Muon_pT = -9999;
                            EMu.Muon_eta = -9999;
                            EMu.Muon_phi = -9999;
                            EMu.Muon_Energy = -9999;
                            EMu.Muon_charge = -9999;
                            EMu.Muon_TuneP_pT = -9999;

                            EMu.Muon_pT_uncorr = mu_noRC.Pt;
                            EMu.Muon_eta_uncorr = mu_noRC.eta;
                            EMu.Muon_phi_uncorr = mu_noRC.phi;
                            EMu.Muon_charge_uncorr = mu_noRC.charge;

                            EMu.Electron_pT = ele_noRC.Pt;
                            EMu.Electron_eta = ele_noRC.eta;
                            EMu.Electron_phi = ele_noRC.phi;
                            EMu.Electron_Energy = ele_noRC.Energy;
                            EMu.Electron_charge = ele_noRC.charge;
                            EMu.Electron_etaSC = ele_noRC.etaSC;
                            EMu.Electron_phiSC = ele_noRC.phiSC;

                            EMu.Electron_Energy_uncorr = ele_noRC.Energy_uncorr;
                        }

                        EMu._prefiringweight = ntuple->_prefiringweight;
                        EMu._prefiringweightup = ntuple->_prefiringweightup;
                        EMu._prefiringweightdown = ntuple->_prefiringweightdown;
                        EMu.PVz = ntuple->PVz;

                        EMu.nVertices = ntuple->nVertices;
                        EMu.nPileUp = ntuple->nPileUp;

                        EMuTree->Fill();

                    } // End of event selection

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

            // Writing
            cout << "Writing into file...";
            Int_t write;
            write = EMuTree->Write();
            if (write)
            {
                cout << " Finished." << endl << "Closing a file..." << endl;
                TString addition = "";
                if (Debug == kTRUE) addition = "_DEBUG";
                EMuFile->Close();
                if (!EMuFile->IsOpen()) cout << "File SelectedEMu_" << Mgr.Tag[i_tup]+addition << ".root has been closed successfully.\n" << endl;
                else cout << "FILE SelectedEMu_" << Mgr.Tag[i_tup]+addition << ".root COULD NOT BE CLOSED!\n" << endl;
            }
            else
            {
                cout << " Writing was NOT successful!\n" << endl;
                EMuFile->Close();
            }

        }// End of i_tup iteration

    }// End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// End of MakeSelectedEMu


/// ----------------------------- For QCD dataset that fails ------------------------------ ///
void MakeSelectedQCDEM_120to170 (TString HLTname, Int_t name, Bool_t Debug)
{
    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer(HLTname);

    FileMgr Mgr(_QCDEMEnriched_120to170);

    cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
    cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

    cout << "\t<" << Mgr.Tag[0] << ">" << endl;
    cout << "\tntuple_skim_" << name << ".root" << endl;

    stringstream ss;
    ss << name;
    TString Name = ss.str();

    //Creating a file
    TFile* ElectronFile;
    if (Debug == kTRUE)
        ElectronFile = TFile::Open("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/QCDfail/SelectedEE_"
                                    +Mgr.Tag[0]+"_"+Name+"_DEBUG.root", "RECREATE");
    else
        ElectronFile = TFile::Open("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/QCDfail/SelectedEE_"
                                    +Mgr.Tag[0]+"_"+Name+".root", "RECREATE");
    ElectronFile->cd();

    TTree* ElectronTree = new TTree("DYTree", "DYTree");
    // -- Creating LongSelectedEE variables to assign branches -- //
    SelectedEE_t EE; EE.CreateNew();

    ElectronTree->Branch("isSelPassed", &EE.isSelPassed);
    ElectronTree->Branch("nVertices", &EE.nVertices);
    ElectronTree->Branch("nPileUp", &EE.nPileUp);
    ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
    ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
    ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
    ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
    ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
    ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
    ElectronTree->Branch("Electron_charge", &EE.Electron_charge);
    ElectronTree->Branch("Electron_etaSC", &EE.Electron_etaSC);
    ElectronTree->Branch("Electron_phiSC", &EE.Electron_phiSC);

    TChain *chain = new TChain(Mgr.TreeName[0]);
    chain->Add(Mgr.BaseLocation+"QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000/ntuple_skim_"+Name+".root");

    NtupleHandle *ntuple = new NtupleHandle(chain);
    if (Mgr.isMC == kTRUE)
    {
        ntuple->TurnOnBranches_GenLepton(); // for all leptons
        ntuple->TurnOnBranches_GenOthers(); // for quarks
    }
    ntuple->TurnOnBranches_Electron();

    Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

    Int_t timesPassed = 0;    

    Int_t NEvents = chain->GetEntries();
    if (Debug == kTRUE) NEvents = 100; // using few events for debugging

    myProgressBar_t bar(NEvents);
    cout << "\tNumber of events: " << NEvents << endl;

    // Loop for all events in the chain
    for (Int_t i=0; i<NEvents; i++)
    {
        ntuple->GetEvent(i);

        // -- Positive/Negative Gen-weights -- //
        ntuple->GENEvt_weight < 0 ? EE.GENEvt_weight = -1 : EE.GENEvt_weight = 1;
        SumWeight += EE.GENEvt_weight;
        SumWeightRaw += ntuple->GENEvt_weight;

        // -- Separate DYLL samples -- //
        Bool_t GenFlag = kFALSE;
        GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[0], ntuple);

        // -- Separate ttbar samples -- //
        Bool_t GenFlag_top = kFALSE;
        vector<GenOthers> GenTopCollection;
        GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[0], ntuple, &GenTopCollection);

        if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += EE.GENEvt_weight;

        Bool_t TriggerFlag = kFALSE;
        TriggerFlag = ntuple->isTriggered(analyzer->HLT);

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
            }

            // -- Event Selection -- //
            vector< Electron > SelectedElectronCollection;
            vector< Int_t > Sel_Index; // Ntuple indexes of electrons that passed the selection
            Bool_t isPassEventSelection = kFALSE;
            isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection, &Sel_Index);

            if (isPassEventSelection == kTRUE)
            {
                timesPassed++;
                Electron ele1 = SelectedElectronCollection[0];
                Electron ele2 = SelectedElectronCollection[1];

                EE.isSelPassed = kTRUE;
                EE.nVertices = ntuple->nVertices;
                EE.nPileUp = ntuple->nPileUp;
                EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

                if (Sel_Index.size() != 2) cout << "======== ERROR: The number of electrons saved is not 2 ========" << endl;
                else
                {
                    for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                    {
                        Int_t index = Sel_Index[iter];

                        EE.Electron_pT->push_back(ntuple->Electron_pT[index]);
                        EE.Electron_eta->push_back(ntuple->Electron_eta[index]);
                        EE.Electron_phi->push_back(ntuple->Electron_phi[index]);
                        EE.Electron_Energy->push_back(ntuple->Electron_Energy[index]);
                        EE.Electron_charge->push_back(ntuple->Electron_charge[index]);
                        EE.Electron_etaSC->push_back(ntuple->Electron_etaSC[index]);
                        EE.Electron_phiSC->push_back(ntuple->Electron_phiSC[index]);

                    } // End of vector filling

                } // End of else()

                ElectronTree->Fill();

                EE.Electron_pT->clear();
                EE.Electron_eta->clear();
                EE.Electron_phi->clear();
                EE.Electron_Energy->clear();
                EE.Electron_charge->clear();
                EE.Electron_etaSC->clear();
                EE.Electron_phiSC->clear();

            } // End of event selection

        } // End of if(isTriggered)

        bar.Draw(i);
    } // End of event iteration

    cout << "\t" << timesPassed << " events have passed the event selection." << endl;

    // Writing
    cout << "\tWriting into file...";
    Int_t write;
    write = ElectronTree->Write();
    if (write)
    {
        cout << " Finished." << endl << "\tClosing a file..." << endl;
        TString addition = "";
        if (Debug == kTRUE) addition = "_DEBUG";
        ElectronFile->Close();
        if (!ElectronFile->IsOpen()) cout << "\tFile SelectedEE_" << Mgr.Tag[0]+addition << "_" << Name << ".root has been closed successfully.\n" << endl;
        else cout << "\tFILE SelectedEE_" << Mgr.Tag[0]+addition << "_" << Name << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!\n" << endl;
        ElectronFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "RunTime: " << TotalRunTime << " seconds" << endl;

} // End of MakeSelectedQCDEM_120to170

/// ----------------------------- To merge small selected QCD files ------------------------------ ///
void MakeSelectedQCDEM_120to170_merged()
{
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr(_QCDEMEnriched_120to170);

    cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
    cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

    cout << "\t<" << Mgr.Tag[0] << ">" << endl;
    cout << "\tMerging all selected events from ntuples that didn't fail into a single file:" << endl;

    //Creating a file
    TFile* ElectronFile = TFile::Open("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_bkg/SelectedEE_"+Mgr.Tag[0]+".root", "RECREATE");
    ElectronFile->cd();

    TTree* ElectronTree = new TTree("DYTree", "DYTree");
    // -- Creating LongSelectedEE variables to assign branches -- //
    SelectedEE_t EE; EE.CreateNew();

    ElectronTree->Branch("isSelPassed", &EE.isSelPassed);
    ElectronTree->Branch("nVertices", &EE.nVertices);
    ElectronTree->Branch("nPileUp", &EE.nPileUp);
    ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
    ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
    ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
    ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
    ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
    ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
    ElectronTree->Branch("Electron_charge", &EE.Electron_charge);
    ElectronTree->Branch("Electron_etaSC", &EE.Electron_etaSC);
    ElectronTree->Branch("Electron_phiSC", &EE.Electron_phiSC);

    TChain *chain = new TChain("DYTree");
    chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedEE/QCDfail/*.root");

    SelectedEE_t* QCD_EE = new SelectedEE_t();
    QCD_EE->CreateFromChain(chain);

    Int_t NEvents = chain->GetEntries();
    myProgressBar_t bar(NEvents);
    cout << "\tNumber of events: " << NEvents << endl;

    // Loop for all events in the chain
    for (Int_t i=0; i<NEvents; i++)
    {
        QCD_EE->GetEvent(i);

        if (QCD_EE->isSelPassed == kTRUE)
        {
            EE.isSelPassed = kTRUE;
            EE.GENEvt_weight = QCD_EE->GENEvt_weight;
            EE.nVertices = QCD_EE->nVertices;
            EE.nPileUp = QCD_EE->nPileUp;
            EE.Electron_InvM = QCD_EE->Electron_InvM;

            if (QCD_EE->Electron_charge->size() != 2) cout << "======== ERROR: Vector sizes are not 2 ========" << endl;
            else
            {
                for (UInt_t iter=0; iter<2; iter++)
                {
                    EE.Electron_pT->push_back(QCD_EE->Electron_pT->at(iter));
                    EE.Electron_eta->push_back(QCD_EE->Electron_eta->at(iter));
                    EE.Electron_phi->push_back(QCD_EE->Electron_phi->at(iter));
                    EE.Electron_Energy->push_back(QCD_EE->Electron_Energy->at(iter));
                    EE.Electron_charge->push_back(QCD_EE->Electron_charge->at(iter));
                    EE.Electron_etaSC->push_back(QCD_EE->Electron_etaSC->at(iter));
                    EE.Electron_phiSC->push_back(QCD_EE->Electron_phiSC->at(iter));

                } // End of vector filling

            } // End of else()

            ElectronTree->Fill();

            EE.Electron_pT->clear();
            EE.Electron_eta->clear();
            EE.Electron_phi->clear();
            EE.Electron_Energy->clear();
            EE.Electron_charge->clear();
            EE.Electron_etaSC->clear();
            EE.Electron_phiSC->clear();

        } // End of event selection

        bar.Draw(i);
    } // End of event iteration

    // Writing
    cout << "\tWriting into file...";
    Int_t write;
    write = ElectronTree->Write();
    if (write)
    {
        cout << " Finished." << endl << "\tClosing a file..." << endl;        
        ElectronFile->Close();
        if (!ElectronFile->IsOpen()) cout << "\tFile SelectedEE_" << Mgr.Tag[0] << ".root has been closed successfully.\n" << endl;
        else cout << "\tFILE SelectedEE_" << Mgr.Tag[0] << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!\n" << endl;
        ElectronFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "RunTime: " << TotalRunTime << " seconds" << endl;

} // End of MakeSelectedQCDEM_120to170_merged
