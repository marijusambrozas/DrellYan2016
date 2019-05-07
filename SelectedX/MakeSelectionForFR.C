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
//    if (whichX.Contains("E"))
//    {
//        Xselected++;
//        if (HLTname == "DEFAULT") HLT = "Ele23Ele12";
//        else HLT = HLTname;
//        cout << "\n*******      MakeSelectionForFR_E (" << type << ", " << HLT << ")      *******" << endl;
//        MakeSelectionForFR_E(type, HLT, Debug);
//    }
    if (whichX.Contains("MU"))
    {
        Xselected++;
        if (HLTname == "DEFAULT") HLT = "Mu50";
        else HLT = HLTname;
        cout << "\n*****  MakeSelectionForFR_Mu (" << type << ", " << HLT << ")  *****" << endl;
        MakeSelectionForFR_Mu(type, HLT, Debug);
    }

    if (Xselected == 0) {
        cout << "Wrong arument!" << endl;
        fs << "Process MakeSelectionForFR (" << whichX << ", " << type << ", " << HLTname << ") FAILED. Wrong arguments!\n";
    }
    else fs << "Process MakeSelectionForFR (" << whichX << ", " << type << ", " << HLTname << ") FINISHED.\n";

    fs.close();

} // End of MakeSelectedX()


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
            chain->Add(Mgr.FullLocation[i_tup]);

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
//                out_dir = "Data/SelectedMuMu_"+Mgr.Tag[i_tup];

            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "SIGNAL")
        {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_signal/SelectedMuMu_"+Mgr.Tag[i_tup];

            out_base = "/cms/ldap_home/mambroza/DrellYan2016/";
            out_dir = "SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc];
        }
        else if (Mgr.Type == "BKG")
        {
//                out_base = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/mambroza/SelectedX_v1/SelectedMuMu/";
//                out_dir = "MC_bkg/SelectedMuMu_"+Mgr.Tag[i_tup];

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
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        double evt_weight;

        const int ptbinnum_endcap = 9;
        double ptbin_endcap[ptbinnum_endcap+1] = {47,52,60,70,80,90,100,150,200,500};
        const int ptbinnum = 17;
        double ptbin[ptbinnum+1] = {47,52,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500};

        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", ptbinnum, ptbin); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", ptbinnum_endcap, ptbin_endcap); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", ptbinnum, ptbin); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", ptbinnum_endcap, ptbin_endcap); h_pT_endcap_deno->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 48, -2.4, 2.4); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 48, -2.4, 2.4); h_eta_deno->Sumw2();
        TH1D* h_iso_barrel_nume = new TH1D("h_iso_barrel_nume", "h_iso_barrel_nume", 100, 0, 5); h_iso_barrel_nume->Sumw2();
        TH1D* h_iso_endcap_nume = new TH1D("h_iso_endcap_nume", "h_iso_endcap_nume", 100, 0, 5); h_iso_endcap_nume->Sumw2();
        TH1D* h_iso_barrel_deno = new TH1D("h_iso_barrel_deno", "h_iso_barrel_deno", 100, 0, 5); h_iso_barrel_deno->Sumw2();
        TH1D* h_iso_endcap_deno = new TH1D("h_iso_endcap_deno", "h_iso_endcap_deno", 100, 0, 5); h_iso_endcap_deno->Sumw2();

        TTree* MuonTree = new TTree("FRTree", "FRTree");
        // -- Creating SelectedMuMu variables to assign branches -- //
        MuonTree->Branch("p_T", &p_T);
        MuonTree->Branch("eta", &eta);
        MuonTree->Branch("charge", &charge);
        MuonTree->Branch("relPFiso", &relPFiso);
        MuonTree->Branch("evt_weight", &evt_weight);

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
//            if (Mgr.CurrentProc == _WJets) i_tup = 2;

            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            chain->Add(Mgr.FullLocation[i_tup]);

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
                ntuple->GENEvt_weight < 0 ? evt_weight = -1 : evt_weight = 1;
                SumWeight += evt_weight;
                SumWeightRaw += ntuple->GENEvt_weight;

                // -- Separate DYLL samples -- //
                Bool_t GenFlag = kTRUE;//kFALSE;
//                GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Mgr.Tag[i_tup], ntuple);

                // -- Separate ttbar samples -- //
                Bool_t GenFlag_top = kTRUE;//kFALSE;
//                vector<GenOthers> GenTopCollection;
//                GenFlag_top = analyzer->Separate_ttbarSample(Mgr.Tag[i_tup], ntuple, &GenTopCollection);

//                if (GenFlag == kTRUE && GenFlag_top == kTRUE) SumWeight_Separated += MuMu.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = ntuple->isTriggered(analyzer->HLT);

                if (TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE)
                {
//                    if (Mgr.Tag[i_tup].Contains("ttbar"))
//                    {
//                        // -- Top pT reweighting -- //
//                        Double_t SF0 = exp(0.0615 - (0.0005 * GenTopCollection[0].Pt));
//                        Double_t SF1 = exp(0.0615 - (0.0005 * GenTopCollection[1].Pt));
//                        evt_weight *= sqrt(SF0 * SF1);
//                    }

                    // -- Reco level selection -- //
                    vector< Muon > MuonCollection;
//                    vector< Muon > MuonCollection_noRocCorr;
                    Int_t NLeptons = ntuple->nMuon;
                    for (Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Muon mu;
                        mu.FillFromNtuple(ntuple, i_reco);

//                        MuonCollection_noRocCorr.push_back(mu);

                        // -- Rochester correction -- //
                        Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                        Int_t s, m;
                        if(Mgr.Tag[i_tup] == "DATA")
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
                    vector< Muon > SelectedMuonCollection_nume, SelectedMuonCollection_deno, SelectedMuonCollection_noRocCorr_nume, SelectedMuonCollection_noRocCorr_deno;
                    Bool_t isPassEventSelection = kFALSE;
//                    Bool_t isPassEventSelection_noRocCorr = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_FR(MuonCollection, ntuple, &SelectedMuonCollection_nume, &SelectedMuonCollection_deno);
//                    isPassEventSelection_noRocCorr = analyzer->EventSelection(MuonCollection_noRocCorr, ntuple, &SelectedMuonCollection_noRocCorr_nume, &SelectedMuonCollection_noRocCorr_deno);

                    if (isPassEventSelection == kTRUE)
                    {
                        timesPassed++;
                        p_T->clear();
                        eta->clear();
                        charge->clear();
                        relPFiso->clear();
                        for (UInt_t i=0; i<SelectedMuonCollection_deno.size(); i++)
                        {
                            p_T->push_back(SelectedMuonCollection_deno[i].Pt);
                            eta->push_back(SelectedMuonCollection_deno[i].eta);
                            charge->push_back(SelectedMuonCollection_deno[i].charge);
                            relPFiso->push_back(SelectedMuonCollection_deno[i].RelPFIso_dBeta);

                            int weight = 1;
                            if (Mgr.isMC)
                                weight = evt_weight*L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup];
                            cout << evt_weight << endl;

                            h_eta_deno->Fill(SelectedMuonCollection_deno[i].eta, weight);
                            if (fabs(SelectedMuonCollection_deno[i].eta) < 1.2)
                            {
                                h_pT_barrel_deno->Fill(SelectedMuonCollection_deno[i].Pt, weight);
                                h_iso_barrel_deno->Fill(SelectedMuonCollection_deno[i].RelPFIso_dBeta, weight);
                            }
                            else
                            {
                                h_pT_endcap_deno->Fill(SelectedMuonCollection_deno[i].Pt, weight);
                                h_iso_endcap_deno->Fill(SelectedMuonCollection_deno[i].RelPFIso_dBeta, weight);
                            }
                            if (i < SelectedMuonCollection_nume.size())
                            {
                                h_eta_nume->Fill(SelectedMuonCollection_nume[i].eta, weight);
                                if (fabs(SelectedMuonCollection_nume[i].eta) < 1.2)
                                {
                                    h_pT_barrel_nume->Fill(SelectedMuonCollection_nume[i].Pt, weight);
                                    h_iso_barrel_nume->Fill(SelectedMuonCollection_nume[i].RelPFIso_dBeta, weight);
                                }
                                else
                                {
                                    h_pT_endcap_nume->Fill(SelectedMuonCollection_nume[i].Pt, weight);
                                    h_iso_endcap_nume->Fill(SelectedMuonCollection_nume[i].RelPFIso_dBeta, weight);
                                }
                            }
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
        Int_t write;
        write = MuonTree->Write();

        h_pT_barrel_deno->SetDirectory(0); h_pT_barrel_deno->Write();
        h_pT_endcap_deno->SetDirectory(0); h_pT_endcap_deno->Write();
        h_pT_barrel_nume->SetDirectory(0); h_pT_barrel_nume->Write();
        h_pT_endcap_nume->SetDirectory(0); h_pT_endcap_nume->Write();
        h_iso_barrel_deno->SetDirectory(0); h_iso_barrel_deno->Write();
        h_iso_endcap_deno->SetDirectory(0); h_iso_endcap_deno->Write();
        h_iso_barrel_nume->SetDirectory(0); h_iso_barrel_nume->Write();
        h_iso_endcap_nume->SetDirectory(0); h_iso_endcap_nume->Write();
        h_eta_deno->SetDirectory(0); h_eta_deno->Write();
        h_eta_nume->SetDirectory(0); h_eta_nume->Write();

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
            cout << " Writing was NOT successful!\n" << endl;
            MuonFile->Close();
        }

//        TFile * file_h = new TFile(out_base+"FRhists_"+Mgr.Procname[Mgr.CurrentProc]+addition+".root", "RECREATE");
//        file_h->cd();
//        h_pT_barrel_deno->SetDirectory(0); h_pT_barrel_deno->Write();
//        h_pT_endcap_deno->SetDirectory(0); h_pT_endcap_deno->Write();
//        h_pT_barrel_nume->SetDirectory(0); h_pT_barrel_nume->Write();
//        h_pT_endcap_nume->SetDirectory(0); h_pT_endcap_nume->Write();
//        h_iso_barrel_deno->SetDirectory(0); h_iso_barrel_deno->Write();
//        h_iso_endcap_deno->SetDirectory(0); h_iso_endcap_deno->Write();
//        h_iso_barrel_nume->SetDirectory(0); h_iso_barrel_nume->Write();
//        h_iso_endcap_nume->SetDirectory(0); h_iso_endcap_nume->Write();
//        h_eta_deno->SetDirectory(0); h_eta_deno->Write();
//        h_eta_nume->SetDirectory(0); h_eta_nume->Write();

//        cout << " Histogram writing finished." << endl << "Closing a file..." << endl;
//        file_h->Close();
//        if (!file_h->IsOpen()) cout << "File FRhists_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root has been closed successfully.\n" << endl;
//        else cout << "FILE FRhists_" << Mgr.Procname[Mgr.CurrentProc]+addition << ".root COULD NOT BE CLOSED!\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MakeSelectedMuMu
