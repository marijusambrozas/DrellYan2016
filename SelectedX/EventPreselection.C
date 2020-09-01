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

#include "./header/myProgressBar_t.cc"
#include "./header/FileMgr.h"
#include "./etc/RoccoR/RoccoR.cc"

// "type" is a string with approximate file name, for example type="SinglePhoton_runB"
// will provide ntuples corresponding to SinglePhoton run 2016B dataset
//
// Set selection to kFALSE if you want to save all the events with at least one
// single photon trigger fired and no further event selection (but keeps only HLT info)
void EventPreselection (TString type, Bool_t selection = kTRUE, Bool_t Debug = kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;   

    TStopwatch totaltime;
    totaltime.Start();

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
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        //Creating a file
        TString out_base;
        TString out_dir;
        TFile* ElectronFile;
        out_base = "/media/sf_DATA/test/";
        out_dir = "SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc];
        TString addition = "";
        if (Debug) addition = "_DEBUG";
        ElectronFile = TFile::Open(out_base+out_dir+addition+".root", "RECREATE");
        ElectronFile->cd();

        // Variables to fill with event data
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
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<int> *trig_fired = new std::vector<int>; // This keeps the "trigger name": 22 if HLT_Photon22_v* has fired and so on
        std::vector<int> *trig_matched = new std::vector<int>; // This keeps the number of electron that was matched to the HLT object
        std::vector<double> *trig_pT = new std::vector<double>;
        Int_t runNum;
        Int_t lumiBlock;

        // Output tree
        TTree* ElectronTree = new TTree("FRTree", "FRTree");
        ElectronTree->Branch("p_T", &p_T);
        ElectronTree->Branch("eta", &eta);
        ElectronTree->Branch("etaSC", &etaSC);
        ElectronTree->Branch("phi", &phi);
        ElectronTree->Branch("charge", &charge);
        ElectronTree->Branch("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        ElectronTree->Branch("dEtaInSeed", dEtaInSeed);
        ElectronTree->Branch("dPhiIn", &dPhiIn);
        ElectronTree->Branch("HoverE", &HoverE);
        ElectronTree->Branch("InvEminusInvP", &InvEminusInvP);
        ElectronTree->Branch("mHits", &mHits);
        ElectronTree->Branch("passConvVeto", &passConvVeto);
        ElectronTree->Branch("relPFiso_Rho", &relPFiso_Rho);
        ElectronTree->Branch("passMediumID", &passMediumID);
        ElectronTree->Branch("trig_fired", &trig_fired);
        ElectronTree->Branch("trig_matched", &trig_matched);
        ElectronTree->Branch("trig_pT", &trig_pT);
        ElectronTree->Branch("runNum", &runNum);
        ElectronTree->Branch("lumiBlock", &lumiBlock);

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]); // Setting up the chain
            Mgr.SetupChain(i_tup, chain);
            NtupleHandle *ntuple = new NtupleHandle(chain); // Ntuple handle helps getting the info from ntuples
            ntuple->TurnOnBranches_Electron();

            Int_t NEvents = chain->GetEntries();
            Int_t timesPassed = 0;
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging
            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);
                if (!Debug) bar.Draw(i);

                if (!selection) // No event preselection, saving only HLT info
                {
                    Bool_t pass = kFALSE;
                    trig_fired->clear();
                    trig_pT->clear();
                    for( Int_t k = 0; k < ntuple->HLT_ntrig; k++ )
                    {
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon22_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(22);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon30_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(30);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon36_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(36);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon50_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(50);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon75_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(75);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon90_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(90);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon120_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(120);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }
                        if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon175_v*" && ntuple->HLT_trigFired[k] )
                        {
                            pass = kTRUE;
                            trig_fired->push_back(175);
                            trig_pT->push_back(ntuple->HLT_trigPt[k]);
                        }

                        if (pass)
                        {
                            runNum = ntuple->runNum;
                            lumiBlock = ntuple->lumiBlock;
                            ElectronTree->Fill();
                        }

                    } // End of for(triggers)

                } // End of if(!selection)

                else // Performing some event preselection
                {
                    Bool_t TriggerFlag = kFALSE;
                    for (Int_t k = 0; k < ntuple->HLT_ntrig; k++) // Checking if any of the triggers has fired
                    {
                        if ( ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon22_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon30_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon36_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon50_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon75_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon90_v*"  && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon120_v*" && ntuple->HLT_trigFired[k] ) ||
                             ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon175_v*" && ntuple->HLT_trigFired[k] ) )
                        {
                            TriggerFlag = kTRUE;
                            break;
                        }
                    } // End of for(triggers)

                    if (TriggerFlag == kTRUE)
                    {
                        // Making the vectors empty for new event
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
                        mHits->clear();
                        passConvVeto->clear();
                        relPFiso_Rho->clear();
                        passMediumID->clear();
                        trig_fired->clear();
                        trig_matched->clear();
                        trig_pT->clear();

                        // -- Event Selection -- //
                        // Go through all electrons, select only trigger matched ones that surpass minimum thresholds
                        Bool_t isPassEventSelection = kFALSE;
                        for (Int_t j=0; j<ntuple->Nelectrons; j++)
                        {
                            if(ntuple->Electron_pT[j] > 25 && fabs(ntuple->Electron_etaSC[j]) < 2.4 && ntuple->Electron_mHits[j] <= 3 &&
                               !(fabs(ntuple->Electron_etaSC[j]) > 1.4442 && fabs(ntuple->Electron_etaSC[j]) < 1.566)) // very relaxed preselection
                            {
                                Bool_t pass = kFALSE;
                                for (Int_t k = 0; k < ntuple->HLT_ntrig; k++) // HLT matching
                                {
                                    if ( ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon22_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon30_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon36_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon50_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon75_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon90_v*"  && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon120_v*" && ntuple->HLT_trigFired[k] ) ||
                                         ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon175_v*" && ntuple->HLT_trigFired[k] ) )
                                    {
                                        Double_t dR = sqrt((ntuple->Electron_eta[j] - ntuple->HLT_trigEta[k])*(ntuple->Electron_eta[j] - ntuple->HLT_trigEta[k]) +
                                                           (ntuple->Electron_phi[j] - ntuple->HLT_trigPhi[k])*(ntuple->Electron_phi[j] - ntuple->HLT_trigPhi[k]));
                                        if (dR < 0.3)
                                        { // Saving all HLTs that have been matched with selected electron
                                            trig_matched->push_back(k);
                                            trig_pT->push_back(ntuple->HLT_trigPt);
                                            // Reducing output size: storing only trigger threshold value instead of full name
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon22_v*" ) trig_fired->push_back(22);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon30_v*" ) trig_fired->push_back(30);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon36_v*" ) trig_fired->push_back(36);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon50_v*" ) trig_fired->push_back(50);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon75_v*" ) trig_fired->push_back(75);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon90_v*" ) trig_fired->push_back(90);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon120_v*") trig_fired->push_back(120);
                                            if ( (ntuple->HLT_trigName->at((unsigned int)k)) == "HLT_Photon175_v*") trig_fired->push_back(175);
                                            pass = kTRUE;
                                        }
                                    } // End of if(HLT fired)
                                } // End of for(triggers)

                                if (pass) // Fill electron vectors
                                {
                                    p_T->push_back(ntuple->Electron_pT[j]);
                                    eta->push_back(ntuple->Electron_eta[j]);
                                    etaSC->push_back(ntuple->Electron_etaSC[j]);
                                    phi->push_back(ntuple->Electron_phi[j]);
                                    charge->push_back(ntuple->Electron_charge[j]);
                                    Full5x5_SigmaIEtaIEta->push_back(ntuple->Electron_Full5x5_SigmaIEtaIEta[j]);
                                    dEtaInSeed->push_back(ntuple->Electron_dEtaInSeed[j]);
                                    dPhiIn->push_back(ntuple->Electron_dPhiIn[j]);
                                    HoverE->push_back(ntuple->Electron_HoverE[j]);
                                    InvEminusInvP->push_back(ntuple->Electron_InvEminusInvP[j]);
                                    mHits->push_back(ntuple->Electron_mHits[j]);
                                    passConvVeto->push_back(ntuple->Electron_passConvVeto[j]);
                                    relPFiso_Rho->push_back(ntuple->Electron_RelPFIso_Rho[j]);
                                    passMediumID->push_back(ntuple->Electron_passMediumID[j]);

                                    isPassEventSelection = kTRUE;
                                }
                            } // End of if(pass preselection)
                        } // End of for(electrons)

                        if (isPassEventSelection == kTRUE)
                        {
                            timesPassed++;
                            runNum = ntuple->runNum;
                            lumiBlock = ntuple->lumiBlock;
                            ElectronTree->Fill();
                        }
                    } // End of if(isTriggered)
                } // End of if(selection)
            } // End of event loop

            cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into file...";
        ElectronFile->cd();
        Int_t write;
        write = ElectronTree->Write();
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

} // End of EventPreselection
