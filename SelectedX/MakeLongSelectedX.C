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

// -- Customized Analyzer for EE selection -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"

void MakeLongSelectedEE ( TString HLTname );
void MakeLongSelectedMuMu ( TString HLTname );
void MakeLongSelectedEMu ( TString HLTname );

void MakeLongSelectedX ( TString whichX, TString HLTname = "DEFAULT" )
{
    TString HLT;
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******      MakeLongSelectedEE ( " << HLT << " )      *******" << endl;
        MakeLongSelectedEE(HLT);
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MakeLongSelectedMuMu ( " << HLT << " )  *****" << endl;
        MakeLongSelectedMuMu(HLT);
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   MakeLongSelectedEMu ( " << HLT << " )  *****" << endl;
        MakeLongSelectedEMu(HLT);
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;
}

/// ----------------------------- Electron Channel ------------------------------ ///
//void MakeLongSelectedEE ( Int_t type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname )
void MakeLongSelectedEE ( TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "LongSelectedEE_ZToEE_M4500to6000_2.root";

//    TString DataType, DataLocation, DataLocation2, Type;

//    if( type == 1 ) {
//            DataType = "B";
//            DataLocation = "DoubleEG_Run2016B";
//    }
//    if( type == 2 ) {
//            DataType = "C";
//            DataLocation = "DoubleEG_Run2016C";
//    }
//    if( type == 3 ) {
//            DataType = "D";
//            DataLocation = "DoubleEG_Run2016D";
//    }
//    if( type == 4 ) {
//            DataType = "E";
//            DataLocation = "DoubleEG_Run2016E";
//    }
//    if( type == 5 ) {
//            DataType = "F";
//            DataLocation = "DoubleEG_Run2016F";
//    }
//    if( type == 6 ) {
//            DataType = "G";
//            DataLocation = "DoubleEG_Run2016G";
//    }
//    if( type == 7 ) {
//            DataType = "H";
//            DataLocation = "DoubleEG_Run2016Hver2";
//            DataLocation2 = "DoubleEG_Run2016Hver3";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
//    // -- Signal MC samples -- //
//    if( type == 11 ) Type = "DYEE_M10to50";
//    if( type == 12 ) Type = "DYEE_M50to100";
//    if( type == 13 ) Type = "DYEE_M100toInf";
//    // -- Background MC samples -- //
//    if( type == 21 ) Type = "ttbar";
//    if( type == 22 ) Type = "ttbarBackup";
//    if( type == 23 ) Type = "ttbar_M700toInf";
//    if( type == 31 ) Type = "DYTauTau_M10to50";
//    if( type == 32 ) Type = "DYTauTau_M50toInf";
//    if( type == 41 ) Type = "VVnST";
//    if( type == 51 ) Type = "WJetsToLNu";
//    //if( type == 61 ) Type = "QCDEMEnriched";

    //Creating a file
    TFile* ElectronFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
    TTree* ElectronTree = new TTree("DYTree", "DYTree");


    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
//    cout << "Type: " << Type << endl;
//    if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

//    TString BaseLocation;
//    if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
//    else BaseLocation = "/data9/DATA/DYntuple/v2.1";
//    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("ZToEE_M4500to6000"); // There is no such type mentioned in DYAnalyzer
    nEvents.push_back(74606);       // The actual number of events in file is taken
    Xsec.push_back(4.56E-07);       // Copied from ZToMuMu_M4500to6000
//    if( Type == "Data" ) {
//            ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating LongSelectedEE variables to assign branches -- //
    LongSelectedEE_t EE; EE.CreateNew();

    ElectronTree->Branch("nVertices", &EE.nVertices);
    ElectronTree->Branch("runNum", &EE.runNum);
    ElectronTree->Branch("lumiBlock", &EE.lumiBlock);
    ElectronTree->Branch("evtNum", &EE.evtNum);
    ElectronTree->Branch("nPileUp", &EE.nPileUp);
    ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
    ElectronTree->Branch("HLT_ntrig", &EE.HLT_ntrig);
    ElectronTree->Branch("HLT_trigFired", &EE.HLT_trigFired);
    ElectronTree->Branch("HLT_trigName", &EE.HLT_trigName);
//    ElectronTree->Branch("HLT_trigPt", &EE.HLT_trigPt);
    ElectronTree->Branch("HLT_trigEta", &EE.HLT_trigEta);
    ElectronTree->Branch("HLT_trigPhi", &EE.HLT_trigPhi);
    ElectronTree->Branch("isHardProcess", &EE.isHardProcess);
    ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
    ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
    ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
    ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
    ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
    ElectronTree->Branch("Electron_charge", &EE.Electron_charge);
    ElectronTree->Branch("Electron_gsfpT", &EE.Electron_gsfpT);
    ElectronTree->Branch("Electron_gsfPx", &EE.Electron_gsfPx);
    ElectronTree->Branch("Electron_gsfPy", &EE.Electron_gsfPy);
    ElectronTree->Branch("Electron_gsfPz", &EE.Electron_gsfPz);
    ElectronTree->Branch("Electron_gsfEta", &EE.Electron_gsfEta);
    ElectronTree->Branch("Electron_gsfPhi", &EE.Electron_gsfEta);
    ElectronTree->Branch("Electron_gsfCharge", &EE.Electron_gsfCharge);
    ElectronTree->Branch("Electron_etaSC", &EE.Electron_etaSC);
    ElectronTree->Branch("Electron_phiSC", &EE.Electron_phiSC);
    ElectronTree->Branch("Electron_etaWidth", &EE.Electron_etaWidth);
    ElectronTree->Branch("Electron_phiWidth", &EE.Electron_phiWidth);
    ElectronTree->Branch("Electron_dEtaIn", &EE.Electron_dEtaIn);
    ElectronTree->Branch("Electron_dEtaInSeed", &EE.Electron_dEtaInSeed);
    ElectronTree->Branch("Electron_dPhiIn", &EE.Electron_dPhiIn);
    ElectronTree->Branch("Electron_sigmaIEtaIEta", &EE.Electron_sigmaIEtaIEta);
    ElectronTree->Branch("Electron_Full5x5_SigmaIEtaIEta", &EE.Electron_Full5x5_SigmaIEtaIEta);
    ElectronTree->Branch("Electron_HoverE", &EE.Electron_HoverE);
    ElectronTree->Branch("Electron_fbrem", &EE.Electron_fbrem);
    ElectronTree->Branch("Electron_eOverP", &EE.Electron_eOverP);
    ElectronTree->Branch("Electron_InvEminusInvP", &EE.Electron_InvEminusInvP);
    ElectronTree->Branch("Electron_dxyVTX", &EE.Electron_dxyVTX);
    ElectronTree->Branch("Electron_dzVTX", &EE.Electron_dzVTX);
    ElectronTree->Branch("Electron_dxy", &EE.Electron_dxy);
    ElectronTree->Branch("Electron_dz", &EE.Electron_dz);
    ElectronTree->Branch("Electron_dxyBS", &EE.Electron_dxyBS);
    ElectronTree->Branch("Electron_dzBS", &EE.Electron_dxyBS);
    ElectronTree->Branch("Electron_chIso03", &EE.Electron_chIso03);
    ElectronTree->Branch("Electron_nhIso03", &EE.Electron_nhIso03);
    ElectronTree->Branch("Electron_phIso03", &EE.Electron_phIso03);
    ElectronTree->Branch("Electron_ChIso03FromPU", &EE.Electron_ChIso03FromPU);
    ElectronTree->Branch("Electron_mHits", &EE.Electron_mHits);
    ElectronTree->Branch("Electron_EnergySC", &EE.Electron_EnergySC);
    ElectronTree->Branch("Electron_preEnergySC", &EE.Electron_preEnergySC);
    ElectronTree->Branch("Electron_rawEnergySC", &EE.Electron_rawEnergySC);
    ElectronTree->Branch("Electron_etSC", &EE.Electron_etSC);
    ElectronTree->Branch("Electron_E15", &EE.Electron_E15);
    ElectronTree->Branch("Electron_E25", &EE.Electron_E25);
    ElectronTree->Branch("Electron_E55", &EE.Electron_E55);
    ElectronTree->Branch("Electron_RelPFIso_dBeta", &EE.Electron_RelPFIso_dBeta);
    ElectronTree->Branch("Electron_RelPFIso_Rho", &EE.Electron_RelPFIso_Rho);
    ElectronTree->Branch("Electron_r9", &EE.Electron_r9);
    ElectronTree->Branch("Electron_ecalDriven", &EE.Electron_ecalDriven);
    ElectronTree->Branch("Electron_passConvVeto", &EE.Electron_passConvVeto);
//    ElectronTree->Branch("Electron_passLooseID", &EE.Electron_passLooseID);
    ElectronTree->Branch("Electron_passMediumID", &EE.Electron_passMediumID);
//    ElectronTree->Branch("Electron_passTightID", &EE.Electron_passTightID);
//    ElectronTree->Branch("Electron_passMVAID_WP80", &EE.Electron_passMVAID_WP80);
//    ElectronTree->Branch("Electron_passMVAID_WP90", &EE.Electron_passMVAID_WP90);
//    ElectronTree->Branch("Electron_passHEEPID", &EE.Electron_passHEEPID);

    //Loop for all samples
//	const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;       //Using just 1 ntuple for test
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << Tag[i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
        chain->Add("/media/sf_DATA/test/ZToEE_M4500to6000_2.root");

//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//                chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//                if( type == 7 ) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;		// test using few events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);
        Int_t timesPassed = 0;

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? EE.GENEvt_weight = -1 : EE.GENEvt_weight = 1;
            SumWeight += EE.GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kTRUE;
            //Bool_t GenFlag_top = kFALSE;
            //vector<GenOthers> GenTopCollection;
            //GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

            // -- Normalization -- //
            Double_t TotWeight = EE.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EE.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                EE.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                vector< Electron > ElectronCollection;
                Int_t NLeptons = ntuple->Nelectrons;
                for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                {
                    Electron ele;
                    ele.FillFromNtuple(ntuple, i_reco);
                    ElectronCollection.push_back( ele );
                }

                // -- Event Selection -- //
                vector< Electron > SelectedElectronCollection;
                vector< Int_t > Sel_Index; // Ntuple indexes of electrons that passed the selection
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection, &Sel_Index);

                if( isPassEventSelection == kTRUE )
                {
                    timesPassed++;
                    Electron ele1 = SelectedElectronCollection[0];
                    Electron ele2 = SelectedElectronCollection[1];

                    EE.nVertices = ntuple->nVertices;
                    EE.runNum = ntuple->runNum;
                    EE.lumiBlock = ntuple->lumiBlock;
                    EE.evtNum = ntuple->evtNum;
                    EE.nPileUp = ntuple->nPileUp;
                    EE.HLT_ntrig = ntuple->HLT_ntrig;

                    EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

                    Int_t zero_count = 0; // resolving if there is no more information in arrays

//                    for (Int_t iter=0; iter<1000; iter++)
//                    {
//                        if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
//                        {
//                            EE.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
//                            EE.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
//                            EE.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
//                            if(iter<ntuple->HLT_ntrig) EE.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
//                        }
//                        else
//                        {
//                            zero_count++;
//                            if(zero_count>1) break;
//                        }
//                    }

                    for (Int_t iter=0; iter<ntuple->HLT_ntrig; iter++)
                    {
                        EE.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                        EE.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                        EE.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                        EE.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                    }

                    if(Sel_Index.size()!=2) cout << "======== ERROR: The number of electrons saved is not 2 ========" << endl;
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
                            EE.Electron_gsfpT->push_back(ntuple->Electron_gsfpT[index]);
                            EE.Electron_gsfPx->push_back(ntuple->Electron_gsfPx[index]);
                            EE.Electron_gsfPy->push_back(ntuple->Electron_gsfPy[index]);
                            EE.Electron_gsfPz->push_back(ntuple->Electron_gsfPz[index]);
                            EE.Electron_gsfEta->push_back(ntuple->Electron_gsfEta[index]);
                            EE.Electron_gsfPhi->push_back(ntuple->Electron_gsfPhi[index]);
                            EE.Electron_gsfCharge->push_back(ntuple->Electron_gsfCharge[index]);
                            EE.Electron_etaSC->push_back(ntuple->Electron_etaSC[index]);
                            EE.Electron_phiSC->push_back(ntuple->Electron_phiSC[index]);
                            EE.Electron_etaWidth->push_back(ntuple->Electron_etaWidth[index]);
                            EE.Electron_phiWidth->push_back(ntuple->Electron_phiWidth[index]);
                            EE.Electron_dEtaIn->push_back(ntuple->Electron_dEtaIn[index]);
                            EE.Electron_dEtaInSeed->push_back(ntuple->Electron_dEtaInSeed[index]);
                            EE.Electron_dPhiIn->push_back(ntuple->Electron_dPhiIn[index]);
                            EE.Electron_sigmaIEtaIEta->push_back(ntuple->Electron_sigmaIEtaIEta[index]);
                            EE.Electron_Full5x5_SigmaIEtaIEta->push_back(ntuple->Electron_Full5x5_SigmaIEtaIEta[index]);
                            EE.Electron_HoverE->push_back(ntuple->Electron_HoverE[index]);
                            EE.Electron_fbrem->push_back(ntuple->Electron_fbrem[index]);
                            EE.Electron_eOverP->push_back(ntuple->Electron_eOverP[index]);
                            EE.Electron_InvEminusInvP->push_back(ntuple->Electron_InvEminusInvP[index]);
                            EE.Electron_dxyVTX->push_back(ntuple->Electron_dxyVTX[index]);
                            EE.Electron_dzVTX->push_back(ntuple->Electron_dzVTX[index]);
                            EE.Electron_dxy->push_back(ntuple->Electron_dxy[index]);
                            EE.Electron_dz->push_back(ntuple->Electron_dz[index]);
                            EE.Electron_dxyBS->push_back(ntuple->Electron_dxyBS[index]);
                            EE.Electron_dzBS->push_back(ntuple->Electron_dzBS[index]);
                            EE.Electron_chIso03->push_back(ntuple->Electron_chIso03[index]);
                            EE.Electron_nhIso03->push_back(ntuple->Electron_nhIso03[index]);
                            EE.Electron_phIso03->push_back(ntuple->Electron_phIso03[index]);
                            EE.Electron_ChIso03FromPU->push_back(ntuple->Electron_ChIso03FromPU[index]);
                            EE.Electron_mHits->push_back(ntuple->Electron_mHits[index]);
                            EE.Electron_EnergySC->push_back(ntuple->Electron_EnergySC[index]);
                            EE.Electron_preEnergySC->push_back(ntuple->Electron_preEnergySC[index]);
                            EE.Electron_rawEnergySC->push_back(ntuple->Electron_rawEnergySC[index]);
                            EE.Electron_etSC->push_back(ntuple->Electron_etSC[index]);
                            EE.Electron_E15->push_back(ntuple->Electron_E15[index]);
                            EE.Electron_E25->push_back(ntuple->Electron_E25[index]);
                            EE.Electron_E55->push_back(ntuple->Electron_E55[index]);
                            EE.Electron_RelPFIso_dBeta->push_back(ntuple->Electron_RelPFIso_dBeta[index]);
                            EE.Electron_RelPFIso_Rho->push_back(ntuple->Electron_RelPFIso_Rho[index]);
                            EE.Electron_r9->push_back(ntuple->Electron_r9[index]);
                            EE.Electron_ecalDriven->push_back(ntuple->Electron_ecalDriven[index]);
                            EE.Electron_passConvVeto->push_back(ntuple->Electron_passConvVeto[index]);
//                            EE.Electron_passLooseID->push_back(ntuple->Electron_passLooseID[index]);
                            EE.Electron_passMediumID->push_back(ntuple->Electron_passMediumID[index]);
//                            EE.Electron_passTightID->push_back(ntuple->Electron_passTightID[index]);
//                            EE.Electron_passMVAID_WP80->push_back(ntuple->Electron_passMVAID_WP80[index]);
//                            EE.Electron_passMVAID_WP90->push_back(ntuple->Electron_passMVAID_WP90[index]);
//                            EE.Electron_passHEEPID->push_back(ntuple->Electron_passHEEPID[index]);

                        } // End of vector filling

                    } // End of else()

                    ElectronTree->Fill();

                    EE.HLT_trigFired->clear();
                    EE.HLT_trigEta->clear();
                    EE.HLT_trigPhi->clear();
                    EE.HLT_trigName->clear();
                    EE.Electron_pT->clear();
                    EE.Electron_eta->clear();
                    EE.Electron_phi->clear();
                    EE.Electron_Energy->clear();
                    EE.Electron_charge->clear();
                    EE.Electron_gsfpT->clear();
                    EE.Electron_gsfPx->clear();
                    EE.Electron_gsfPy->clear();
                    EE.Electron_gsfPz->clear();
                    EE.Electron_gsfEta->clear();
                    EE.Electron_gsfPhi->clear();
                    EE.Electron_gsfCharge->clear();
                    EE.Electron_etaSC->clear();
                    EE.Electron_phiSC->clear();
                    EE.Electron_etaWidth->clear();
                    EE.Electron_phiWidth->clear();
                    EE.Electron_dEtaIn->clear();
                    EE.Electron_dEtaInSeed->clear();
                    EE.Electron_dPhiIn->clear();
                    EE.Electron_sigmaIEtaIEta->clear();
                    EE.Electron_Full5x5_SigmaIEtaIEta->clear();
                    EE.Electron_HoverE->clear();
                    EE.Electron_fbrem->clear();
                    EE.Electron_eOverP->clear();
                    EE.Electron_InvEminusInvP->clear();
                    EE.Electron_dxyVTX->clear();
                    EE.Electron_dzVTX->clear();
                    EE.Electron_dxy->clear();
                    EE.Electron_dz->clear();
                    EE.Electron_dxyBS->clear();
                    EE.Electron_dzBS->clear();
                    EE.Electron_chIso03->clear();
                    EE.Electron_nhIso03->clear();
                    EE.Electron_phIso03->clear();
                    EE.Electron_ChIso03FromPU->clear();
                    EE.Electron_mHits->clear();
                    EE.Electron_EnergySC->clear();
                    EE.Electron_preEnergySC->clear();
                    EE.Electron_rawEnergySC->clear();
                    EE.Electron_etSC->clear();
                    EE.Electron_E15->clear();
                    EE.Electron_E25->clear();
                    EE.Electron_E55->clear();
                    EE.Electron_RelPFIso_dBeta->clear();
                    EE.Electron_RelPFIso_Rho->clear();
                    EE.Electron_r9->clear();
                    EE.Electron_ecalDriven->clear();
                    EE.Electron_passConvVeto->clear();
//                    EE.Electron_passLooseID->clear();
                    EE.Electron_passMediumID->clear();
//                    EE.Electron_passTightID->clear();
//                    EE.Electron_passMVAID_WP80->clear();
//                    EE.Electron_passMVAID_WP90->clear();
//                    EE.Electron_passHEEPID->clear();

//                    Double_t reco_Pt = (ele1.Momentum + ele2.Momentum).Pt();
//                    Double_t reco_rapi = (ele1.Momentum + ele2.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration
        cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if ( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    ElectronFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = ElectronTree->Write();
    if ( write )
    {
        cout  << " Finished." << endl << "Closing a file..." << endl;
        ElectronFile->Close();
        if ( !ElectronFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        ElectronFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

/// -------------------------------- Muon Channel ------------------------------------ ///
//void MakeLongSelectedMuMu ( Int_t type, TString HLTname )
void MakeLongSelectedMuMu ( TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "LongSelectedMuMu_ZToMuMu_M4500to6000_4.root";

//    TString DataType, DataLocation, DataLocation2, Type, OutputName;

//    if( type == 1 ) {
//        DataType = "B";
//        DataLocation = "SingleMuon_Run2016B";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016B.root";
//    }
//    if( type == 2 ) {
//        DataType = "C";
//        DataLocation = "SingleMuon_Run2016C";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016C.root";
//    }
//    if( type == 3 ) {
//        DataType = "D";
//        DataLocation = "SingleMuon_Run2016D";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016D.root";
//    }
//    if( type == 4 ) {
//        DataType = "E";
//        DataLocation = "SingleMuon_Run2016E";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016E.root";
//    }
//    if( type == 5 ) {
//        DataType = "F";
//        DataLocation = "SingleMuon_Run2016F";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016F.root";
//    }
//    if( type == 6 ) {
//        DataType = "G";
//        DataLocation = "SingleMuon_Run2016G";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016G.root";
//    }
//    if( type == 7 ) {
//        DataType = "H";
//        DataLocation = "SingleMuon_Run2016Hver2";
//        DataLocation2 = "SingleMuon_Run2016Hver3";
//        OutputName = "LongSelectedMuMu_SingleMuon_Run2016H.root";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
//    // -- Signal MC samples -- //
//    if( type == 11 ) {
//        Type = "DYMuMu_M10to50";
//        OutputName = "LongSelectedMuMu_DYMuMu_M10to50.root";
//    }
//    if( type == 12 ) {
//        Type = "DYMuMu_M50to100";
//        OutputName = "LongSelectedMuMu_DYMuMu_M50to100.root";
//    }
//    if( type == 13 ) {
//        Type = "DYMuMu_M100toInf";
//        OutputName = "LongSelectedMuMu_DYMuMu_M100toInf.root";
//    }
//    // -- Background MC samples -- //
//    if( type == 21 ) {
//        Type = "ttbar";
//        OutputName = "LongSelectedMuMu_ttbar.root";
//    }
//    if( type == 22 ) {
//        Type = "ttbarBackup";
//        OutputName = "LongSelectedMuMu_ttbarBackup.root";
//    }
//    if( type == 23 ) {
//        Type = "ttbar_M700toInf";
//        OutputName = "LongSelectedMuMu_ttbar_M700toInf.root";
//    }
//    if( type == 31 ) {
//        Type = "DYTauTau_M10to50";
//        OutputName = "LongSelectedMuMu_DYTauTau_M10to50.root";
//    }
//    if( type == 32 ) {
//        Type = "DYTauTau_M50toInf";
//        OutputName = "LongSelectedMuMu_DYTauTau_M50toInf.root";
//    }
//    if( type == 41 ) {
//        Type = "VVnST";
//        OutputName = "LongSelectedMuMu_VVnST.root";
//    }
//    if( type == 51 ) {
//        Type = "WJetsToLNu";
//        OutputName = "LongSelectedMuMu_WJetsToLNu.root";
//    }
//    if( type == 61 ) {
//        Type = "QCDMuEnriched";
//        OutputName = "LongSelectedMuMu_QCDMuEnriched.root";
//    }

    //Creating a file
    TFile* MuonFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
    TTree* MuonTree = new TTree("DYTree", "DYTree");

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
//    cout << "Type: " << Type << endl;
//    if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

//    TString BaseLocation;
//    if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
//    else BaseLocation = "/data9/DATA/DYntuple/v2.1";
//    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("ZMuMu_M4500to6000");
    nEvents.push_back(100000);
    Xsec.push_back(4.56E-07);
//    if( Type == "Data" ) {
//        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating LongSelectedMuMu variables to assign branches -- //
    LongSelectedMuMu_t MuMu; MuMu.CreateNew();

    MuonTree->Branch("nVertices", &MuMu.nVertices);
    MuonTree->Branch("runNum", &MuMu.runNum);
    MuonTree->Branch("lumiBlock", &MuMu.lumiBlock);
    MuonTree->Branch("evtNum", &MuMu.evtNum);
    MuonTree->Branch("nPileUp", &MuMu.nPileUp);
    MuonTree->Branch("GENEvt_weight", &MuMu.GENEvt_weight);
    MuonTree->Branch("HLT_ntrig", &MuMu.HLT_ntrig);
    MuonTree->Branch("HLT_trigFired", &MuMu.HLT_trigFired);
    MuonTree->Branch("HLT_trigName", &MuMu.HLT_trigName);
//    MuonTree->Branch("HLT_trigPt", &MuMu.HLT_trigPt);
    MuonTree->Branch("HLT_trigEta", &MuMu.HLT_trigEta);
    MuonTree->Branch("HLT_trigPhi", &MuMu.HLT_trigPhi);
    MuonTree->Branch("isHardProcess", &MuMu.isHardProcess);
    MuonTree->Branch("Muon_pT", &MuMu.Muon_pT);
    MuonTree->Branch("Muon_eta", &MuMu.Muon_eta);
    MuonTree->Branch("Muon_phi", &MuMu.Muon_phi);
    MuonTree->Branch("isGLBmuon", &MuMu.isGLBmuon);
    MuonTree->Branch("isPFmuon", &MuMu.isPFmuon);
    MuonTree->Branch("isTRKmuon", &MuMu.isTRKmuon);
    MuonTree->Branch("Muon_charge", &MuMu.Muon_charge);
    MuonTree->Branch("Muon_chi2dof", &MuMu.Muon_chi2dof);
    MuonTree->Branch("Muon_muonHits", &MuMu.Muon_muonHits);
    MuonTree->Branch("Muon_nSegments", &MuMu.Muon_nSegments);
    MuonTree->Branch("Muon_nMatches", &MuMu.Muon_nMatches);
    MuonTree->Branch("Muon_trackerLayers", &MuMu.Muon_trackerLayers);
    MuonTree->Branch("Muon_pixelHits", &MuMu.Muon_pixelHits);
    MuonTree->Branch("Muon_dxyVTX", &MuMu.Muon_dxyVTX);
    MuonTree->Branch("Muon_dzVTX", &MuMu.Muon_dzVTX);
    MuonTree->Branch("Muon_trkiso", &MuMu.Muon_trkiso);
    MuonTree->Branch("Muon_PfChargedHadronIsoR04", &MuMu.Muon_PfChargedHadronIsoR04);
    MuonTree->Branch("Muon_PfNeutralHadronIsoR04", &MuMu.Muon_PfNeutralHadronIsoR04);
    MuonTree->Branch("Muon_PfGammaIsoR04", &MuMu.Muon_PfGammaIsoR04);
    MuonTree->Branch("Muon_PFSumPUIsoR04", &MuMu.Muon_PFSumPUIsoR04);
    MuonTree->Branch("Muon_Px", &MuMu.Muon_Px);
    MuonTree->Branch("Muon_Py", &MuMu.Muon_Py);
    MuonTree->Branch("Muon_Pz", &MuMu.Muon_Pz);
    MuonTree->Branch("Muon_Energy", &MuMu.Muon_Energy);
    MuonTree->Branch("Muon_InvM", &MuMu.Muon_InvM);
    MuonTree->Branch("Muon_Best_pT", &MuMu.Muon_Best_pT);
    MuonTree->Branch("Muon_Best_pTError", &MuMu.Muon_Best_pTError);
    MuonTree->Branch("Muon_Best_Px", &MuMu.Muon_Best_Px);
    MuonTree->Branch("Muon_Best_Py", &MuMu.Muon_Best_Py);
    MuonTree->Branch("Muon_Best_Pz", &MuMu.Muon_Best_Pz);
    MuonTree->Branch("Muon_Best_eta", &MuMu.Muon_Best_eta);
    MuonTree->Branch("Muon_Best_phi", &MuMu.Muon_Best_phi);
    MuonTree->Branch("Muon_Inner_pT", &MuMu.Muon_Inner_pT);
    MuonTree->Branch("Muon_Inner_pTError", &MuMu.Muon_Inner_pTError);
    MuonTree->Branch("Muon_Inner_Px", &MuMu.Muon_Inner_Px);
    MuonTree->Branch("Muon_Inner_Py", &MuMu.Muon_Inner_Py);
    MuonTree->Branch("Muon_Inner_Pz", &MuMu.Muon_Inner_Pz);
    MuonTree->Branch("Muon_Inner_eta", &MuMu.Muon_Inner_eta);
    MuonTree->Branch("Muon_Inner_phi", &MuMu.Muon_Inner_phi);
    MuonTree->Branch("Muon_Outer_pT", &MuMu.Muon_Outer_pT);
    MuonTree->Branch("Muon_Outer_pTError", &MuMu.Muon_Outer_pTError);
    MuonTree->Branch("Muon_Outer_Px", &MuMu.Muon_Outer_Px);
    MuonTree->Branch("Muon_Outer_Py", &MuMu.Muon_Outer_Py);
    MuonTree->Branch("Muon_Outer_Pz", &MuMu.Muon_Outer_Pz);
    MuonTree->Branch("Muon_Outer_eta", &MuMu.Muon_Outer_eta);
    MuonTree->Branch("Muon_Outer_phi", &MuMu.Muon_Outer_phi);
    MuonTree->Branch("Muon_GLB_pT", &MuMu.Muon_GLB_pT);
    MuonTree->Branch("Muon_GLB_pTError", &MuMu.Muon_GLB_pTError);
    MuonTree->Branch("Muon_GLB_Px", &MuMu.Muon_GLB_Px);
    MuonTree->Branch("Muon_GLB_Py", &MuMu.Muon_GLB_Py);
    MuonTree->Branch("Muon_GLB_Pz", &MuMu.Muon_GLB_Pz);
    MuonTree->Branch("Muon_GLB_eta", &MuMu.Muon_GLB_eta);
    MuonTree->Branch("Muon_GLB_phi", &MuMu.Muon_GLB_phi);
    MuonTree->Branch("Muon_TuneP_pT", &MuMu.Muon_TuneP_pT);
    MuonTree->Branch("Muon_TuneP_pTError", &MuMu.Muon_TuneP_pTError);
    MuonTree->Branch("Muon_TuneP_Px", &MuMu.Muon_TuneP_Px);
    MuonTree->Branch("Muon_TuneP_Py", &MuMu.Muon_TuneP_Py);
    MuonTree->Branch("Muon_TuneP_Pz", &MuMu.Muon_TuneP_Pz);
    MuonTree->Branch("Muon_TuneP_eta", &MuMu.Muon_TuneP_eta);
    MuonTree->Branch("Muon_TuneP_phi", &MuMu.Muon_TuneP_phi);
    MuonTree->Branch("CosAngle", &MuMu.CosAngle);
    MuonTree->Branch("vtxTrkChi2", &MuMu.vtxTrkChi2);
    MuonTree->Branch("vtxTrkProb", &MuMu.vtxTrkProb);
    MuonTree->Branch("vtxTrkNdof", &MuMu.vtxTrkNdof);
    MuonTree->Branch("vtxTrkCkt1Pt", &MuMu.vtxTrkCkt1Pt);
    MuonTree->Branch("vtxTrkCkt2Pt", &MuMu.vtxTrkCkt2Pt);

    //Loop for all samples
//    const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;       //Using just 1 ntuple for test
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//            chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//            if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }
        chain->Add("/media/sf_DATA/test/ZToMuMu_M4500to6000_4.root");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Muon();

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? MuMu.GENEvt_weight = -1 : MuMu.GENEvt_weight = 1;
            SumWeight += MuMu.GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kTRUE;

            // -- Normalization -- //
            Double_t TotWeight = MuMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                MuMu.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                vector< Muon > MuonCollection;
                Int_t NLeptons = ntuple->nMuon;
                for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                {
                    Muon mu;
                    mu.FillFromNtuple(ntuple, i_reco);

                    // -- Convert to TuneP variables -- //
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );
                }

                // -- Event Selection -- //
                vector< Muon > SelectedMuonCollection;
                vector< Int_t > Sel_Index;  // Ntuple indexes of muons that passed the selection
                Int_t IndexDi;      // Index of the member of the dimuon vectors that describes the 2 selected muons
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_Zdiff_13TeV_HighPt(MuonCollection, ntuple, &SelectedMuonCollection, &Sel_Index, IndexDi);

                if( isPassEventSelection == kTRUE )
                {
                    Muon mu1 = SelectedMuonCollection[0];
                    Muon mu2 = SelectedMuonCollection[1];

                    MuMu.nVertices = ntuple->nVertices;
                    MuMu.runNum = ntuple->runNum;
                    MuMu.lumiBlock = ntuple->lumiBlock;
                    MuMu.evtNum = ntuple->evtNum;
                    MuMu.nPileUp = ntuple->nPileUp;
                    MuMu.HLT_ntrig = ntuple->HLT_ntrig;

                    MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();

                    if( IndexDi!=-1 )
                    {
                        MuMu.CosAngle->push_back(ntuple->CosAngle->at(IndexDi));
                        MuMu.vtxTrkChi2->push_back(ntuple->vtxTrkChi2->at(IndexDi));
                        MuMu.vtxTrkProb->push_back(ntuple->vtxTrkProb->at(IndexDi));
                        MuMu.vtxTrkNdof->push_back(ntuple->vtxTrkNdof->at(IndexDi));
                        MuMu.vtxTrkCkt1Pt->push_back(ntuple->vtxTrkCkt1Pt->at(IndexDi));
                        MuMu.vtxTrkCkt2Pt->push_back(ntuple->vtxTrkCkt2Pt->at(IndexDi));
                    }
                    else cout << "== ERROR: Event selection was passed but no dimuon index registered. ==" << endl;

                    Int_t zero_count = 0; // resolving if there is no more information in arrays

//                    for (Int_t iter=0; iter<1000; iter++)
//                    {
//                        if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
//                        {
//                            MuMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
//                            MuMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
//                            MuMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
//                            if(iter<ntuple->HLT_ntrig) MuMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
//                        }
//                        else
//                        {
//                            zero_count++;
//                            if(zero_count>1) break;
//                        }
//                    }

                    for (Int_t iter=0; iter<ntuple->HLT_ntrig; iter++)      // There are more trigEta and trigPhi values than nTrig suggests
                    {
                        MuMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                        MuMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                        MuMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                        MuMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                    }

                    if(Sel_Index.size()!=2) cout << "======== ERROR: The number of muons saved is not 2 ========" << endl;
                    else
                    {
                        for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                        {
                            Int_t index = Sel_Index[iter];

                            MuMu.Muon_pT->push_back(ntuple->Muon_pT[index]);
                            MuMu.Muon_eta->push_back(ntuple->Muon_eta[index]);
                            MuMu.Muon_phi->push_back(ntuple->Muon_phi[index]);
                            MuMu.isGLBmuon->push_back(ntuple->isGLBmuon[index]);
                            MuMu.isPFmuon->push_back(ntuple->isPFmuon[index]);
                            MuMu.isTRKmuon->push_back(ntuple->isTRKmuon[index]);
                            MuMu.Muon_charge->push_back(ntuple->Muon_charge[index]);
                            MuMu.Muon_chi2dof->push_back(ntuple->Muon_chi2dof[index]);
                            MuMu.Muon_muonHits->push_back(ntuple->Muon_muonHits[index]);
                            MuMu.Muon_nSegments->push_back(ntuple->Muon_nSegments[index]);
                            MuMu.Muon_nMatches->push_back(ntuple->Muon_nMatches[index]);
                            MuMu.Muon_trackerLayers->push_back(ntuple->Muon_trackerLayers[index]);
                            MuMu.Muon_pixelHits->push_back(ntuple->Muon_pixelHits[index]);
                            MuMu.Muon_dxyVTX->push_back(ntuple->Muon_dxyVTX[index]);
                            MuMu.Muon_dzVTX->push_back(ntuple->Muon_dzVTX[index]);
                            MuMu.Muon_trkiso->push_back(ntuple->Muon_trkiso[index]);
                            MuMu.Muon_PfChargedHadronIsoR04->push_back(ntuple->Muon_PfChargedHadronIsoR04[index]);
                            MuMu.Muon_PfNeutralHadronIsoR04->push_back(ntuple->Muon_PfNeutralHadronIsoR04[index]);
                            MuMu.Muon_PfGammaIsoR04->push_back(ntuple->Muon_PfGammaIsoR04[index]);
                            MuMu.Muon_PFSumPUIsoR04->push_back(ntuple->Muon_PFSumPUIsoR04[index]);
                            MuMu.Muon_Px->push_back(ntuple->Muon_Px[index]);
                            MuMu.Muon_Py->push_back(ntuple->Muon_Py[index]);
                            MuMu.Muon_Pz->push_back(ntuple->Muon_Pz[index]);

                            Double_t Mu_E = sqrt( ntuple->Muon_Px[index]*ntuple->Muon_Px[index] + ntuple->Muon_Py[index]*ntuple->Muon_Py[index]
                                                  + ntuple->Muon_Pz[index]*ntuple->Muon_Pz[index] + M_Mu*M_Mu );

                            MuMu.Muon_Energy->push_back(Mu_E);

                            MuMu.Muon_Best_pT->push_back(ntuple->Muon_Best_pT[index]);
                            MuMu.Muon_Best_pTError->push_back(ntuple->Muon_Best_pTError[index]);
                            MuMu.Muon_Best_Px->push_back(ntuple->Muon_Best_Px[index]);
                            MuMu.Muon_Best_Py->push_back(ntuple->Muon_Best_Py[index]);
                            MuMu.Muon_Best_Pz->push_back(ntuple->Muon_Best_Pz[index]);
                            MuMu.Muon_Best_eta->push_back(ntuple->Muon_Best_eta[index]);
                            MuMu.Muon_Best_phi->push_back(ntuple->Muon_Best_phi[index]);
                            MuMu.Muon_Inner_pT->push_back(ntuple->Muon_Inner_pT[index]);
                            MuMu.Muon_Inner_pTError->push_back(ntuple->Muon_Inner_pTError[index]);
                            MuMu.Muon_Inner_Px->push_back(ntuple->Muon_Inner_Px[index]);
                            MuMu.Muon_Inner_Py->push_back(ntuple->Muon_Inner_Py[index]);
                            MuMu.Muon_Inner_Pz->push_back(ntuple->Muon_Inner_Pz[index]);
                            MuMu.Muon_Inner_eta->push_back(ntuple->Muon_Inner_eta[index]);
                            MuMu.Muon_Inner_phi->push_back(ntuple->Muon_Inner_phi[index]);
                            MuMu.Muon_Outer_pT->push_back(ntuple->Muon_Outer_pT[index]);
                            MuMu.Muon_Outer_pTError->push_back(ntuple->Muon_Outer_pTError[index]);
                            MuMu.Muon_Outer_Px->push_back(ntuple->Muon_Outer_Px[index]);
                            MuMu.Muon_Outer_Py->push_back(ntuple->Muon_Outer_Py[index]);
                            MuMu.Muon_Outer_Pz->push_back(ntuple->Muon_Outer_Pz[index]);
                            MuMu.Muon_Outer_eta->push_back(ntuple->Muon_Outer_eta[index]);
                            MuMu.Muon_Outer_phi->push_back(ntuple->Muon_Outer_phi[index]);
                            MuMu.Muon_GLB_pT->push_back(ntuple->Muon_GLB_pT[index]);
                            MuMu.Muon_GLB_pTError->push_back(ntuple->Muon_GLB_pTError[index]);
                            MuMu.Muon_GLB_Px->push_back(ntuple->Muon_GLB_Px[index]);
                            MuMu.Muon_GLB_Py->push_back(ntuple->Muon_GLB_Py[index]);
                            MuMu.Muon_GLB_Pz->push_back(ntuple->Muon_GLB_Pz[index]);
                            MuMu.Muon_GLB_eta->push_back(ntuple->Muon_GLB_eta[index]);
                            MuMu.Muon_GLB_phi->push_back(ntuple->Muon_GLB_phi[index]);
                            MuMu.Muon_TuneP_pT->push_back(ntuple->Muon_TuneP_pT[index]);
                            MuMu.Muon_TuneP_pTError->push_back(ntuple->Muon_TuneP_pTError[index]);
                            MuMu.Muon_TuneP_Px->push_back(ntuple->Muon_TuneP_Px[index]);
                            MuMu.Muon_TuneP_Py->push_back(ntuple->Muon_TuneP_Py[index]);
                            MuMu.Muon_TuneP_Pz->push_back(ntuple->Muon_TuneP_Pz[index]);
                            MuMu.Muon_TuneP_eta->push_back(ntuple->Muon_TuneP_eta[index]);
                            MuMu.Muon_TuneP_phi->push_back(ntuple->Muon_TuneP_phi[index]);
                        } // End of vector filling

                    } // End of else()

                    MuonTree->Fill();

                    MuMu.HLT_trigFired->clear();
                    MuMu.HLT_trigName->clear();
                    MuMu.HLT_trigEta->clear();
                    MuMu.HLT_trigPhi->clear();
                    MuMu.Muon_pT->clear();
                    MuMu.Muon_eta->clear();
                    MuMu.Muon_phi->clear();
                    MuMu.isGLBmuon->clear();
                    MuMu.isPFmuon->clear();
                    MuMu.isTRKmuon->clear();
                    MuMu.Muon_charge->clear();
                    MuMu.Muon_chi2dof->clear();
                    MuMu.Muon_muonHits->clear();
                    MuMu.Muon_nSegments->clear();
                    MuMu.Muon_nMatches->clear();
                    MuMu.Muon_trackerLayers->clear();
                    MuMu.Muon_pixelHits->clear();
                    MuMu.Muon_dxyVTX->clear();
                    MuMu.Muon_dzVTX->clear();
                    MuMu.Muon_trkiso->clear();
                    MuMu.Muon_PfChargedHadronIsoR04->clear();
                    MuMu.Muon_PfNeutralHadronIsoR04->clear();
                    MuMu.Muon_PfGammaIsoR04->clear();
                    MuMu.Muon_PFSumPUIsoR04->clear();
                    MuMu.Muon_Px->clear();
                    MuMu.Muon_Py->clear();
                    MuMu.Muon_Pz->clear();
                    MuMu.Muon_Energy->clear();
                    MuMu.Muon_Best_pT->clear();
                    MuMu.Muon_Best_pTError->clear();
                    MuMu.Muon_Best_Px->clear();
                    MuMu.Muon_Best_Py->clear();
                    MuMu.Muon_Best_Pz->clear();
                    MuMu.Muon_Best_eta->clear();
                    MuMu.Muon_Best_phi->clear();
                    MuMu.Muon_Inner_pT->clear();
                    MuMu.Muon_Inner_pTError->clear();
                    MuMu.Muon_Inner_Px->clear();
                    MuMu.Muon_Inner_Py->clear();
                    MuMu.Muon_Inner_Pz->clear();
                    MuMu.Muon_Inner_eta->clear();
                    MuMu.Muon_Inner_phi->clear();
                    MuMu.Muon_Outer_pT->clear();
                    MuMu.Muon_Outer_pTError->clear();
                    MuMu.Muon_Outer_Px->clear();
                    MuMu.Muon_Outer_Py->clear();
                    MuMu.Muon_Outer_Pz->clear();
                    MuMu.Muon_Outer_eta->clear();
                    MuMu.Muon_Outer_phi->clear();
                    MuMu.Muon_GLB_pT->clear();
                    MuMu.Muon_GLB_pTError->clear();
                    MuMu.Muon_GLB_Px->clear();
                    MuMu.Muon_GLB_Py->clear();
                    MuMu.Muon_GLB_Pz->clear();
                    MuMu.Muon_GLB_eta->clear();
                    MuMu.Muon_GLB_phi->clear();
                    MuMu.Muon_TuneP_pT->clear();;
                    MuMu.Muon_TuneP_pTError->clear();
                    MuMu.Muon_TuneP_Px->clear();
                    MuMu.Muon_TuneP_Py->clear();
                    MuMu.Muon_TuneP_Pz->clear();
                    MuMu.Muon_TuneP_eta->clear();
                    MuMu.Muon_TuneP_phi->clear();
                    MuMu.CosAngle->clear();
                    MuMu.vtxTrkChi2->clear();
                    MuMu.vtxTrkProb->clear();
                    MuMu.vtxTrkNdof->clear();
                    MuMu.vtxTrkCkt1Pt->clear();
                    MuMu.vtxTrkCkt2Pt->clear();

//                    Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
//                    Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    MuonFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = MuonTree->Write();
    if ( write )
    {
        cout << " Finished." << endl << "Closing a file..." << endl;
        MuonFile->Close();
        if ( !MuonFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        MuonFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

/// --------------------------------- EMu events --------------------------------- ///
//void MakeLongSelectedEMu ( Int_t type, TString HLTname )
void MakeLongSelectedEMu ( TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "LongSelectedEMu_WW_34.root";

//    TString DataType, DataLocation, DataLocation2, Type, OutputName;

//    if( type == 1 ) {
//            DataType = "B";
//            DataLocation = "SingleMuon_Run2016B";
//    }
//    else if( type == 2 ) {
//            DataType = "C";
//            DataLocation = "SingleMuon_Run2016C";
//    }
//    else if( type == 3 ) {
//            DataType = "D";
//            DataLocation = "SingleMuon_Run2016D";
//    }
//    else if( type == 4 ) {
//            DataType = "E";
//            DataLocation = "SingleMuon_Run2016E";
//    }
//    else if( type == 5 ) {
//            DataType = "F";
//            DataLocation = "SingleMuon_Run2016F";
//    }
//    else if( type == 6 ) {
//            DataType = "G";
//            DataLocation = "SingleMuon_Run2016G";
//    }
//    else if( type == 7 ) {
//            DataType = "H";
//            DataLocation = "SingleMuon_Run2016Hver2";
//            DataLocation2 = "SingleMuon_Run2016Hver3";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) Type = "Data";
//    // -- Background MC samples -- //
//    else if( type == 21 ) Type = "ttbar";
//    else if( type == 22 ) Type = "ttbarBackup";
//    else if( type == 23 ) Type = "ttbar_M700toInf";
//    else if( type == 31 ) Type = "DYTauTau_M10to50";
//    else if( type == 32 ) Type = "DYTauTau_M50toInf";
//    else if( type == 41 ) Type = "VVnST";

    //Creating a file
    TFile* EMuFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
    TTree* EMuTree = new TTree("DYTree", "DYTree");

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
//    cout << "Type: " << Type << endl;
//    if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

//    TString BaseLocation;
//    if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
//    else BaseLocation = "/data9/DATA/DYntuple/v2.1";
//    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("WW");
    nEvents.push_back(6987123);
    Xsec.push_back(118.7);
//    if( Type == "Data" ) {
//        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating LongSelectedMuMu variables to assign branches -- //
    LongSelectedEMu_t EMu; EMu.CreateNew();

    EMuTree->Branch("nVertices", &EMu.nVertices);
    EMuTree->Branch("runNum", &EMu.runNum);
    EMuTree->Branch("lumiBlock", &EMu.lumiBlock);
    EMuTree->Branch("evtNum", &EMu.evtNum);
    EMuTree->Branch("nPileUp", &EMu.nPileUp);
    EMuTree->Branch("GENEvt_weight", &EMu.GENEvt_weight);
    EMuTree->Branch("HLT_ntrig", &EMu.HLT_ntrig);
    EMuTree->Branch("HLT_trigFired", &EMu.HLT_trigFired);
    EMuTree->Branch("HLT_trigName", &EMu.HLT_trigName);
//    EMuTree->Branch("HLT_trigPt", &EMu.HLT_trigPt);
    EMuTree->Branch("HLT_trigEta", &EMu.HLT_trigEta);
    EMuTree->Branch("HLT_trigPhi", &EMu.HLT_trigPhi);
    EMuTree->Branch("isHardProcess", &EMu.isHardProcess);
    EMuTree->Branch("EMu_InvM", &EMu.EMu_InvM);
    EMuTree->Branch("Muon_pT", &EMu.Muon_pT);
    EMuTree->Branch("Muon_eta", &EMu.Muon_eta);
    EMuTree->Branch("Muon_phi", &EMu.Muon_phi);
    EMuTree->Branch("isGLBmuon", &EMu.isGLBmuon);
    EMuTree->Branch("isPFmuon", &EMu.isPFmuon);
    EMuTree->Branch("isTRKmuon", &EMu.isTRKmuon);
    EMuTree->Branch("Muon_charge", &EMu.Muon_charge);
    EMuTree->Branch("Muon_chi2dof", &EMu.Muon_chi2dof);
    EMuTree->Branch("Muon_muonHits", &EMu.Muon_muonHits);
    EMuTree->Branch("Muon_nSegments", &EMu.Muon_nSegments);
    EMuTree->Branch("Muon_nMatches", &EMu.Muon_nMatches);
    EMuTree->Branch("Muon_trackerLayers", &EMu.Muon_trackerLayers);
    EMuTree->Branch("Muon_pixelHits", &EMu.Muon_pixelHits);
    EMuTree->Branch("Muon_dxyVTX", &EMu.Muon_dxyVTX);
    EMuTree->Branch("Muon_dzVTX", &EMu.Muon_dzVTX);
    EMuTree->Branch("Muon_trkiso", &EMu.Muon_trkiso);
    EMuTree->Branch("Muon_PfChargedHadronIsoR04", &EMu.Muon_PfChargedHadronIsoR04);
    EMuTree->Branch("Muon_PfNeutralHadronIsoR04", &EMu.Muon_PfNeutralHadronIsoR04);
    EMuTree->Branch("Muon_PfGammaIsoR04", &EMu.Muon_PfGammaIsoR04);
    EMuTree->Branch("Muon_PFSumPUIsoR04", &EMu.Muon_PFSumPUIsoR04);
    EMuTree->Branch("Muon_Px", &EMu.Muon_Px);
    EMuTree->Branch("Muon_Py", &EMu.Muon_Py);
    EMuTree->Branch("Muon_Pz", &EMu.Muon_Pz);
    EMuTree->Branch("Muon_Energy", &EMu.Muon_Energy);
    EMuTree->Branch("Muon_Best_pT", &EMu.Muon_Best_pT);
    EMuTree->Branch("Muon_Best_pTError", &EMu.Muon_Best_pTError);
    EMuTree->Branch("Muon_Best_Px", &EMu.Muon_Best_Px);
    EMuTree->Branch("Muon_Best_Py", &EMu.Muon_Best_Py);
    EMuTree->Branch("Muon_Best_Pz", &EMu.Muon_Best_Pz);
    EMuTree->Branch("Muon_Best_eta", &EMu.Muon_Best_eta);
    EMuTree->Branch("Muon_Best_phi", &EMu.Muon_Best_phi);
    EMuTree->Branch("Muon_Inner_pT", &EMu.Muon_Inner_pT);
    EMuTree->Branch("Muon_Inner_pTError", &EMu.Muon_Inner_pTError);
    EMuTree->Branch("Muon_Inner_Px", &EMu.Muon_Inner_Px);
    EMuTree->Branch("Muon_Inner_Py", &EMu.Muon_Inner_Py);
    EMuTree->Branch("Muon_Inner_Pz", &EMu.Muon_Inner_Pz);
    EMuTree->Branch("Muon_Inner_eta", &EMu.Muon_Inner_eta);
    EMuTree->Branch("Muon_Inner_phi", &EMu.Muon_Inner_phi);
    EMuTree->Branch("Muon_Outer_pT", &EMu.Muon_Outer_pT);
    EMuTree->Branch("Muon_Outer_pTError", &EMu.Muon_Outer_pTError);
    EMuTree->Branch("Muon_Outer_Px", &EMu.Muon_Outer_Px);
    EMuTree->Branch("Muon_Outer_Py", &EMu.Muon_Outer_Py);
    EMuTree->Branch("Muon_Outer_Pz", &EMu.Muon_Outer_Pz);
    EMuTree->Branch("Muon_Outer_eta", &EMu.Muon_Outer_eta);
    EMuTree->Branch("Muon_Outer_phi", &EMu.Muon_Outer_phi);
    EMuTree->Branch("Muon_GLB_pT", &EMu.Muon_GLB_pT);
    EMuTree->Branch("Muon_GLB_pTError", &EMu.Muon_GLB_pTError);
    EMuTree->Branch("Muon_GLB_Px", &EMu.Muon_GLB_Px);
    EMuTree->Branch("Muon_GLB_Py", &EMu.Muon_GLB_Py);
    EMuTree->Branch("Muon_GLB_Pz", &EMu.Muon_GLB_Pz);
    EMuTree->Branch("Muon_GLB_eta", &EMu.Muon_GLB_eta);
    EMuTree->Branch("Muon_GLB_phi", &EMu.Muon_GLB_phi);
    EMuTree->Branch("Muon_TuneP_pT", &EMu.Muon_TuneP_pT);
    EMuTree->Branch("Muon_TuneP_pTError", &EMu.Muon_TuneP_pTError);
    EMuTree->Branch("Muon_TuneP_Px", &EMu.Muon_TuneP_Px);
    EMuTree->Branch("Muon_TuneP_Py", &EMu.Muon_TuneP_Py);
    EMuTree->Branch("Muon_TuneP_Pz", &EMu.Muon_TuneP_Pz);
    EMuTree->Branch("Muon_TuneP_eta", &EMu.Muon_TuneP_eta);
    EMuTree->Branch("Muon_TuneP_phi", &EMu.Muon_TuneP_phi);
    EMuTree->Branch("Electron_pT", &EMu.Electron_pT);
    EMuTree->Branch("Electron_eta", &EMu.Electron_eta);
    EMuTree->Branch("Electron_phi", &EMu.Electron_phi);
    EMuTree->Branch("Electron_Energy", &EMu.Electron_Energy);
    EMuTree->Branch("Electron_charge", &EMu.Electron_charge);
    EMuTree->Branch("Electron_gsfpT", &EMu.Electron_gsfpT);
    EMuTree->Branch("Electron_gsfPx", &EMu.Electron_gsfPx);
    EMuTree->Branch("Electron_gsfPy", &EMu.Electron_gsfPy);
    EMuTree->Branch("Electron_gsfPz", &EMu.Electron_gsfPz);
    EMuTree->Branch("Electron_gsfEta", &EMu.Electron_gsfEta);
    EMuTree->Branch("Electron_gsfPhi", &EMu.Electron_gsfEta);
    EMuTree->Branch("Electron_gsfCharge", &EMu.Electron_gsfCharge);
    EMuTree->Branch("Electron_etaSC", &EMu.Electron_etaSC);
    EMuTree->Branch("Electron_phiSC", &EMu.Electron_phiSC);
    EMuTree->Branch("Electron_etaWidth", &EMu.Electron_etaWidth);
    EMuTree->Branch("Electron_phiWidth", &EMu.Electron_phiWidth);
    EMuTree->Branch("Electron_dEtaIn", &EMu.Electron_dEtaIn);
    EMuTree->Branch("Electron_dEtaInSeed", &EMu.Electron_dEtaInSeed);
    EMuTree->Branch("Electron_dPhiIn", &EMu.Electron_dPhiIn);
    EMuTree->Branch("Electron_sigmaIEtaIEta", &EMu.Electron_sigmaIEtaIEta);
    EMuTree->Branch("Electron_Full5x5_SigmaIEtaIEta", &EMu.Electron_Full5x5_SigmaIEtaIEta);
    EMuTree->Branch("Electron_HoverE", &EMu.Electron_HoverE);
    EMuTree->Branch("Electron_fbrem", &EMu.Electron_fbrem);
    EMuTree->Branch("Electron_eOverP", &EMu.Electron_eOverP);
    EMuTree->Branch("Electron_InvEminusInvP", &EMu.Electron_InvEminusInvP);
    EMuTree->Branch("Electron_dxyVTX", &EMu.Electron_dxyVTX);
    EMuTree->Branch("Electron_dzVTX", &EMu.Electron_dzVTX);
    EMuTree->Branch("Electron_dxy", &EMu.Electron_dxy);
    EMuTree->Branch("Electron_dz", &EMu.Electron_dz);
    EMuTree->Branch("Electron_dxyBS", &EMu.Electron_dxyBS);
    EMuTree->Branch("Electron_dzBS", &EMu.Electron_dxyBS);
    EMuTree->Branch("Electron_chIso03", &EMu.Electron_chIso03);
    EMuTree->Branch("Electron_nhIso03", &EMu.Electron_nhIso03);
    EMuTree->Branch("Electron_phIso03", &EMu.Electron_phIso03);
    EMuTree->Branch("Electron_ChIso03FromPU", &EMu.Electron_ChIso03FromPU);
    EMuTree->Branch("Electron_mHits", &EMu.Electron_mHits);
    EMuTree->Branch("Electron_EnergySC", &EMu.Electron_EnergySC);
    EMuTree->Branch("Electron_preEnergySC", &EMu.Electron_preEnergySC);
    EMuTree->Branch("Electron_rawEnergySC", &EMu.Electron_rawEnergySC);
    EMuTree->Branch("Electron_etSC", &EMu.Electron_etSC);
    EMuTree->Branch("Electron_E15", &EMu.Electron_E15);
    EMuTree->Branch("Electron_E25", &EMu.Electron_E25);
    EMuTree->Branch("Electron_E55", &EMu.Electron_E55);
    EMuTree->Branch("Electron_RelPFIso_dBeta", &EMu.Electron_RelPFIso_dBeta);
    EMuTree->Branch("Electron_RelPFIso_Rho", &EMu.Electron_RelPFIso_Rho);
    EMuTree->Branch("Electron_r9", &EMu.Electron_r9);
    EMuTree->Branch("Electron_ecalDriven", &EMu.Electron_ecalDriven);
    EMuTree->Branch("Electron_passConvVeto", &EMu.Electron_passConvVeto);
//    EMuTree->Branch("Electron_passLooseID", &EMu.Electron_passLooseID);
    EMuTree->Branch("Electron_passMediumID", &EMu.Electron_passMediumID);
//    EMuTree->Branch("Electron_passTightID", &EMu.Electron_passTightID);
//    EMuTree->Branch("Electron_passMVAID_WP80", &EMu.Electron_passMVAID_WP80);
//    EMuTree->Branch("Electron_passMVAID_WP90", &EMu.Electron_passMVAID_WP90);
//    EMuTree->Branch("Electron_passHEEPID", &EMu.Electron_passHEEPID);

    //Loop for all samples
//    const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;       //Using just 1 ntuple for test
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//            chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//            if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }
        chain->Add("/media/sf_DATA/test/WW_34.root");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();
        ntuple->TurnOnBranches_Muon();

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? EMu.GENEvt_weight = -1 : EMu.GENEvt_weight = 1;
            SumWeight += EMu.GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kTRUE;

            // -- Normalization -- //
            Double_t TotWeight = EMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                EMu.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                // Muons
                vector< Muon > MuonCollection;
                Int_t N_Muons = ntuple->nMuon;
                for(Int_t i_reco=0; i_reco<N_Muons; i_reco++)
                {
                    Muon mu;
                    mu.FillFromNtuple(ntuple, i_reco);

                    // -- Convert to TuneP variables -- //
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );
                }

                // Electrons
                vector< Electron > ElectronCollection;
                Int_t N_Electrons = ntuple->Nelectrons;
                for(Int_t i_reco=0; i_reco<N_Electrons; i_reco++)
                {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_reco);

                        ElectronCollection.push_back( ele );
                }

                // -- Event Selection -- //
                vector< Muon > SelectedMuonCollection;
                vector< Electron > SelectedElectronCollection;
                Int_t Sel_Index_Mu, Sel_Index_Ele;  // Ntuple indexes of electron and muon that passed the selection
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_emu_method_test(MuonCollection, ElectronCollection, ntuple, &SelectedMuonCollection,
                                                                                &SelectedElectronCollection, Sel_Index_Mu, Sel_Index_Ele);

                if( isPassEventSelection == kTRUE && Sel_Index_Mu != -1 && Sel_Index_Ele != -1 )
                {
                    Muon mu = SelectedMuonCollection[0];
                    Electron ele = SelectedElectronCollection[0];

                    EMu.nVertices = ntuple->nVertices;
                    EMu.runNum = ntuple->runNum;
                    EMu.lumiBlock = ntuple->lumiBlock;
                    EMu.evtNum = ntuple->evtNum;
                    EMu.nPileUp = ntuple->nPileUp;
                    EMu.HLT_ntrig = ntuple->HLT_ntrig;

                    EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();

                    Int_t zero_count = 0; // resolving if there is no more information in arrays

//                    for (Int_t iter=0; iter<1000; iter++)
//                    {
//                        if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
//                        {
//                            EMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
//                            EMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
//                            EMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
//                            if(iter<ntuple->HLT_ntrig) EMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
//                        }
//                        else
//                        {
//                            zero_count++;
//                            if(zero_count>1) break;
//                        }
//                    }

                    for (Int_t iter=0; iter<ntuple->HLT_ntrig; iter++)      // There are more trigEta and trigPhi values than nTrig suggests
                    {
                        EMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                        EMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                        EMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                        EMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                    }

                    EMu.Muon_pT = ntuple->Muon_pT[Sel_Index_Mu];
                    EMu.Muon_eta = ntuple->Muon_eta[Sel_Index_Mu];
                    EMu.Muon_phi = ntuple->Muon_phi[Sel_Index_Mu];
                    EMu.isGLBmuon = ntuple->isGLBmuon[Sel_Index_Mu];
                    EMu.isPFmuon = ntuple->isPFmuon[Sel_Index_Mu];
                    EMu.isTRKmuon = ntuple->isTRKmuon[Sel_Index_Mu];
                    EMu.Muon_charge = ntuple->Muon_charge[Sel_Index_Mu];
                    EMu.Muon_chi2dof = ntuple->Muon_chi2dof[Sel_Index_Mu];
                    EMu.Muon_muonHits = ntuple->Muon_muonHits[Sel_Index_Mu];
                    EMu.Muon_nSegments = ntuple->Muon_nSegments[Sel_Index_Mu];
                    EMu.Muon_nMatches = ntuple->Muon_nMatches[Sel_Index_Mu];
                    EMu.Muon_trackerLayers = ntuple->Muon_trackerLayers[Sel_Index_Mu];
                    EMu.Muon_pixelHits = ntuple->Muon_pixelHits[Sel_Index_Mu];
                    EMu.Muon_dxyVTX = ntuple->Muon_dxyVTX[Sel_Index_Mu];
                    EMu.Muon_dzVTX = ntuple->Muon_dzVTX[Sel_Index_Mu];
                    EMu.Muon_trkiso = ntuple->Muon_trkiso[Sel_Index_Mu];
                    EMu.Muon_PfChargedHadronIsoR04 = ntuple->Muon_PfChargedHadronIsoR04[Sel_Index_Mu];
                    EMu.Muon_PfNeutralHadronIsoR04 = ntuple->Muon_PfNeutralHadronIsoR04[Sel_Index_Mu];
                    EMu.Muon_PfGammaIsoR04 = ntuple->Muon_PfGammaIsoR04[Sel_Index_Mu];
                    EMu.Muon_PFSumPUIsoR04 = ntuple->Muon_PFSumPUIsoR04[Sel_Index_Mu];
                    EMu.Muon_Px = ntuple->Muon_Px[Sel_Index_Mu];
                    EMu.Muon_Py = ntuple->Muon_Py[Sel_Index_Mu];
                    EMu.Muon_Pz = ntuple->Muon_Pz[Sel_Index_Mu];
                    EMu.Muon_Energy = sqrt( ntuple->Muon_Px[Sel_Index_Mu]*ntuple->Muon_Px[Sel_Index_Mu] + ntuple->Muon_Py[Sel_Index_Mu]*ntuple->Muon_Py[Sel_Index_Mu]
                                          + ntuple->Muon_Pz[Sel_Index_Mu]*ntuple->Muon_Pz[Sel_Index_Mu] + M_Mu*M_Mu );
                    EMu.Muon_Best_pT = ntuple->Muon_Best_pT[Sel_Index_Mu];
                    EMu.Muon_Best_pTError = ntuple->Muon_Best_pTError[Sel_Index_Mu];
                    EMu.Muon_Best_Px = ntuple->Muon_Best_Px[Sel_Index_Mu];
                    EMu.Muon_Best_Py = ntuple->Muon_Best_Py[Sel_Index_Mu];
                    EMu.Muon_Best_Pz = ntuple->Muon_Best_Pz[Sel_Index_Mu];
                    EMu.Muon_Best_eta = ntuple->Muon_Best_eta[Sel_Index_Mu];
                    EMu.Muon_Best_phi = ntuple->Muon_Best_phi[Sel_Index_Mu];
                    EMu.Muon_Inner_pT = ntuple->Muon_Inner_pT[Sel_Index_Mu];
                    EMu.Muon_Inner_pTError = ntuple->Muon_Inner_pTError[Sel_Index_Mu];
                    EMu.Muon_Inner_Px = ntuple->Muon_Inner_Px[Sel_Index_Mu];
                    EMu.Muon_Inner_Py = ntuple->Muon_Inner_Py[Sel_Index_Mu];
                    EMu.Muon_Inner_Pz = ntuple->Muon_Inner_Pz[Sel_Index_Mu];
                    EMu.Muon_Inner_eta = ntuple->Muon_Inner_eta[Sel_Index_Mu];
                    EMu.Muon_Inner_phi = ntuple->Muon_Inner_phi[Sel_Index_Mu];
                    EMu.Muon_Outer_pT = ntuple->Muon_Outer_pT[Sel_Index_Mu];
                    EMu.Muon_Outer_pTError = ntuple->Muon_Outer_pTError[Sel_Index_Mu];
                    EMu.Muon_Outer_Px = ntuple->Muon_Outer_Px[Sel_Index_Mu];
                    EMu.Muon_Outer_Py = ntuple->Muon_Outer_Py[Sel_Index_Mu];
                    EMu.Muon_Outer_Pz = ntuple->Muon_Outer_Pz[Sel_Index_Mu];
                    EMu.Muon_Outer_eta = ntuple->Muon_Outer_eta[Sel_Index_Mu];
                    EMu.Muon_Outer_phi = ntuple->Muon_Outer_phi[Sel_Index_Mu];
                    EMu.Muon_GLB_pT = ntuple->Muon_GLB_pT[Sel_Index_Mu];
                    EMu.Muon_GLB_pTError = ntuple->Muon_GLB_pTError[Sel_Index_Mu];
                    EMu.Muon_GLB_Px = ntuple->Muon_GLB_Px[Sel_Index_Mu];
                    EMu.Muon_GLB_Py = ntuple->Muon_GLB_Py[Sel_Index_Mu];
                    EMu.Muon_GLB_Pz = ntuple->Muon_GLB_Pz[Sel_Index_Mu];
                    EMu.Muon_GLB_eta = ntuple->Muon_GLB_eta[Sel_Index_Mu];
                    EMu.Muon_GLB_phi = ntuple->Muon_GLB_phi[Sel_Index_Mu];
                    EMu.Muon_TuneP_pT = ntuple->Muon_TuneP_pT[Sel_Index_Mu];
                    EMu.Muon_TuneP_pTError = ntuple->Muon_TuneP_pTError[Sel_Index_Mu];
                    EMu.Muon_TuneP_Px = ntuple->Muon_TuneP_Px[Sel_Index_Mu];
                    EMu.Muon_TuneP_Py = ntuple->Muon_TuneP_Py[Sel_Index_Mu];
                    EMu.Muon_TuneP_Pz = ntuple->Muon_TuneP_Pz[Sel_Index_Mu];
                    EMu.Muon_TuneP_eta = ntuple->Muon_TuneP_eta[Sel_Index_Mu];
                    EMu.Muon_TuneP_phi = ntuple->Muon_TuneP_phi[Sel_Index_Mu];
                    EMu.Electron_pT = ntuple->Electron_pT[Sel_Index_Ele];
                    EMu.Electron_eta = ntuple->Electron_eta[Sel_Index_Ele];
                    EMu.Electron_phi = ntuple->Electron_phi[Sel_Index_Ele];
                    EMu.Electron_Energy = ntuple->Electron_Energy[Sel_Index_Ele];
                    EMu.Electron_charge = ntuple->Electron_charge[Sel_Index_Ele];
                    EMu.Electron_gsfpT = ntuple->Electron_gsfpT[Sel_Index_Ele];
                    EMu.Electron_gsfPx = ntuple->Electron_gsfPx[Sel_Index_Ele];
                    EMu.Electron_gsfPy = ntuple->Electron_gsfPy[Sel_Index_Ele];
                    EMu.Electron_gsfPz = ntuple->Electron_gsfPz[Sel_Index_Ele];
                    EMu.Electron_gsfEta = ntuple->Electron_gsfEta[Sel_Index_Ele];
                    EMu.Electron_gsfPhi = ntuple->Electron_gsfPhi[Sel_Index_Ele];
                    EMu.Electron_gsfCharge = ntuple->Electron_gsfCharge[Sel_Index_Ele];
                    EMu.Electron_etaSC = ntuple->Electron_etaSC[Sel_Index_Ele];
                    EMu.Electron_phiSC = ntuple->Electron_phiSC[Sel_Index_Ele];
                    EMu.Electron_etaWidth = ntuple->Electron_etaWidth[Sel_Index_Ele];
                    EMu.Electron_phiWidth = ntuple->Electron_phiWidth[Sel_Index_Ele];
                    EMu.Electron_dEtaIn = ntuple->Electron_dEtaIn[Sel_Index_Ele];
                    EMu.Electron_dEtaInSeed = ntuple->Electron_dEtaInSeed[Sel_Index_Ele];
                    EMu.Electron_dPhiIn = ntuple->Electron_dPhiIn[Sel_Index_Ele];
                    EMu.Electron_sigmaIEtaIEta = ntuple->Electron_sigmaIEtaIEta[Sel_Index_Ele];
                    EMu.Electron_Full5x5_SigmaIEtaIEta = ntuple->Electron_Full5x5_SigmaIEtaIEta[Sel_Index_Ele];
                    EMu.Electron_HoverE = ntuple->Electron_HoverE[Sel_Index_Ele];
                    EMu.Electron_fbrem = ntuple->Electron_fbrem[Sel_Index_Ele];
                    EMu.Electron_eOverP = ntuple->Electron_eOverP[Sel_Index_Ele];
                    EMu.Electron_InvEminusInvP = ntuple->Electron_InvEminusInvP[Sel_Index_Ele];
                    EMu.Electron_dxyVTX = ntuple->Electron_dxyVTX[Sel_Index_Ele];
                    EMu.Electron_dzVTX = ntuple->Electron_dzVTX[Sel_Index_Ele];
                    EMu.Electron_dxy = ntuple->Electron_dxy[Sel_Index_Ele];
                    EMu.Electron_dz = ntuple->Electron_dz[Sel_Index_Ele];
                    EMu.Electron_dxyBS = ntuple->Electron_dxyBS[Sel_Index_Ele];
                    EMu.Electron_dzBS = ntuple->Electron_dzBS[Sel_Index_Ele];
                    EMu.Electron_chIso03 = ntuple->Electron_chIso03[Sel_Index_Ele];
                    EMu.Electron_nhIso03 = ntuple->Electron_nhIso03[Sel_Index_Ele];
                    EMu.Electron_phIso03 = ntuple->Electron_phIso03[Sel_Index_Ele];
                    EMu.Electron_ChIso03FromPU = ntuple->Electron_ChIso03FromPU[Sel_Index_Ele];
                    EMu.Electron_mHits = ntuple->Electron_mHits[Sel_Index_Ele];
                    EMu.Electron_EnergySC = ntuple->Electron_EnergySC[Sel_Index_Ele];
                    EMu.Electron_preEnergySC = ntuple->Electron_preEnergySC[Sel_Index_Ele];
                    EMu.Electron_rawEnergySC = ntuple->Electron_rawEnergySC[Sel_Index_Ele];
                    EMu.Electron_etSC = ntuple->Electron_etSC[Sel_Index_Ele];
                    EMu.Electron_E15 = ntuple->Electron_E15[Sel_Index_Ele];
                    EMu.Electron_E25 = ntuple->Electron_E25[Sel_Index_Ele];
                    EMu.Electron_E55 = ntuple->Electron_E55[Sel_Index_Ele];
                    EMu.Electron_RelPFIso_dBeta = ntuple->Electron_RelPFIso_dBeta[Sel_Index_Ele];
                    EMu.Electron_RelPFIso_Rho = ntuple->Electron_RelPFIso_Rho[Sel_Index_Ele];
                    EMu.Electron_r9 = ntuple->Electron_r9[Sel_Index_Ele];
                    EMu.Electron_ecalDriven = ntuple->Electron_ecalDriven[Sel_Index_Ele];
                    EMu.Electron_passConvVeto = ntuple->Electron_passConvVeto[Sel_Index_Ele];
//                            EMu.Electron_passLooseID = ntuple->Electron_passLooseID[Sel_Index_Ele];
                    EMu.Electron_passMediumID = ntuple->Electron_passMediumID[Sel_Index_Ele];
//                            EMu.Electron_passTightID = ntuple->Electron_passTightID[Sel_Index_Ele];
//                            EMu.Electron_passMVAID_WP80 = ntuple->Electron_passMVAID_WP80[Sel_Index_Ele];
//                            EMu.Electron_passMVAID_WP90 = ntuple->Electron_passMVAID_WP90[Sel_Index_Ele];
//                            EMu.Electron_passHEEPID = ntuple->Electron_passHEEPID[Sel_Index_Ele];

                    EMuTree->Fill();

                    EMu.HLT_trigFired->clear();
                    EMu.HLT_trigName->clear();
                    EMu.HLT_trigEta->clear();
                    EMu.HLT_trigPhi->clear();


//                    Double_t reco_Pt = (mu.Momentum + ele.Momentum).Pt();
//                    Double_t reco_rapi = (mu.Momentum + ele.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if( isMC == kTRUE && nEvents.size()>0 ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    EMuFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = EMuTree->Write();
    if ( write )
    {
        cout << " Finished." << endl << "Closing a file..." << endl;
        EMuFile->Close();
        if ( !EMuFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        EMuFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
