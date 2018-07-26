#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TAttMarker.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for MuMu selection -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "header/myProgressBar_t.h"

#define M_Mu 0.1056583715 // -- GeV -- //

static inline void loadBar(int x, int n, int r, int w);

// -- Muon Channel -- //
//void MakeSelectedMuMu(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeReSelectedMuMu(TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
	// -- Run2016 luminosity [/pb] -- //
	Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
	L = L_B2H;
        TString OutputName = "ReSelectedMuMu_ZToMuMu_M4500to6000_4.root";

//        TString DataType, DataLocation, DataLocation2, Type, OutputName;

//	if( type == 1 ) {
//		DataType = "B";
//		DataLocation = "SingleMuon_Run2016B";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016B.root";
//	}
//	if( type == 2 ) {
//		DataType = "C";
//		DataLocation = "SingleMuon_Run2016C";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016C.root";
//	}
//	if( type == 3 ) {
//		DataType = "D";
//		DataLocation = "SingleMuon_Run2016D";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016D.root";
//	}
//	if( type == 4 ) {
//		DataType = "E";
//		DataLocation = "SingleMuon_Run2016E";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016E.root";
//	}
//	if( type == 5 ) {
//		DataType = "F";
//		DataLocation = "SingleMuon_Run2016F";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016F.root";
//	}
//	if( type == 6 ) {
//		DataType = "G";
//		DataLocation = "SingleMuon_Run2016G";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016G.root";
//	}
//	if( type == 7 ) {
//		DataType = "H";
//		DataLocation = "SingleMuon_Run2016Hver2";
//		DataLocation2 = "SingleMuon_Run2016Hver3";
//                OutputName = "SelectedMuMu_SingleMuon_Run2016H.root";
//	}

        Bool_t isMC = kTRUE;
//	if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
//	// -- Signal MC samples -- //
//        if( type == 11 ) {
//            Type = "DYMuMu_M10to50";
//            OutputName = "SelectedMuMu_DYMuMu_M10to50.root";
//        }
//        if( type == 12 ) {
//            Type = "DYMuMu_M50to100";
//            OutputName = "SelectedMuMu_DYMuMu_M50to100.root";
//        }
//        if( type == 13 ) {
//            Type = "DYMuMu_M100toInf";
//            OutputName = "SelectedMuMu_DYMuMu_M100toInf.root";
//        }
//	// -- Background MC samples -- //
//        if( type == 21 ) {
//            Type = "ttbar";
//            OutputName = "SelectedMuMu_ttbar.root";
//        }
//        if( type == 22 ) {
//            Type = "ttbarBackup";
//            OutputName = "SelectedMuMu_ttbarBackup.root";
//        }
//        if( type == 23 ) {
//            Type = "ttbar_M700toInf";
//            OutputName = "SelectedMuMu_ttbar_M700toInf.root";
//        }
//        if( type == 31 ) {
//            Type = "DYTauTau_M10to50";
//            OutputName = "SelectedMuMu_DYTauTau_M10to50.root";
//        }
//        if( type == 32 ) {
//            Type = "DYTauTau_M50toInf";
//            OutputName = "SelectedMuMu_DYTauTau_M50toInf.root";
//        }
//        if( type == 41 ) {
//            Type = "VVnST";
//            OutputName = "SelectedMuMu_VVnST.root";
//        }
//        if( type == 51 ) {
//            Type = "WJetsToLNu";
//            OutputName = "SelectedMuMu_WJetsToLNu.root";
//        }
        //if( type == 61 ) {
//            Type = "QCDMuEnriched";
//            OutputName = "SelectedMuMu_QCDMuEnriched.root";
//        }

        //Creating a file
        TFile* MuonFile = new TFile("/media/sf_DATA/"+OutputName, "RECREATE");
        TTree* MuonTree = new TTree("DYTree", "DYTree");

        TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
//	cout << "Type: " << Type << endl;
//	if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

//	TString BaseLocation;
//	if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
//	else BaseLocation = "/data9/DATA/DYntuple/v2.1";
//	cout << "DATA location: " << BaseLocation << endl;

	TStopwatch totaltime;
	totaltime.Start();

	DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

	// -- Each ntuple directory & corresponding Tags -- //
//	vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
        vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
        Tag.push_back("ZMuMu_M4500to6000");
        nEvents.push_back(100000);
        Xsec.push_back(4.56E-07);
//	if( Type == "Data" ) {
//		ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//	}
//	else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

        // -- Creating SelectedMuMu variables to assign branches -- //
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
//        MuonTree->Branch("HLT_trigPt", &MuMu.HLT_trigPt);
        MuonTree->Branch("HLT_trigEta", &MuMu.HLT_trigEta);
        MuonTree->Branch("HLT_trigPhi", &MuMu.HLT_trigPhi);
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
        MuonTree->Branch("Muon_E", &MuMu.Muon_E);
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
//	const Int_t Ntup = ntupleDirectory.size();
        const Int_t Ntup = 1;
	for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
	{
		TStopwatch looptime;
		looptime.Start();

//		cout << "\t<" << [i_tup] << ">" << endl;

                TChain *chain = new TChain("DYTree");
//		//Set MC chain
//		if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//		//Set Data chain
//		else {
//			chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//			if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//		}
                chain->Add("/media/sf_DATA/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");

                LongSelectedMuMu_t * MuMuInput = new LongSelectedMuMu_t();
                MuMuInput->CreateFromChain(chain);

		Double_t SumWeight = 0, SumWeight_Separated = 0;

                Int_t NEvents = chain->GetEntries();
//		Int_t NEvents = 10000;						// test using small events
		cout << "\t[Total Events: " << NEvents << "]" << endl;
                myProgressBar_t bar(NEvents);

		for(Int_t i=0; i<NEvents; i++)		
		{	                       
                        MuMuInput->GetEvent(i);
			
			// -- Positive/Negative Gen-weights -- //
                        MuMuInput->GENEvt_weight < 0 ? MuMu.GENEvt_weight = -1 : MuMu.GENEvt_weight = 1;
                        SumWeight += MuMu.GENEvt_weight;

			// -- Normalization -- //
                        Double_t TotWeight = MuMu.GENEvt_weight;
                        if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu.GENEvt_weight;

			Bool_t TriggerFlag = kFALSE;
                        TriggerFlag = MuMuInput->isTriggered( analyzer->HLT );

                        if( TriggerFlag == kTRUE )
			{
				////////////////////////////////
				// -- Reco level selection -- //
				////////////////////////////////
                                vector< Muon > MuonCollection;
                                Int_t NLeptons = MuMuInput->Muon_pT->size();
				for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
				{
                                        Muon mu;
                                        mu.FillFromSelectedX(MuMuInput, i_reco);

					// -- Convert to TuneP variables -- //
					analyzer->ConvertToTunePInfo( mu );
					MuonCollection.push_back( mu );
				}	

				// -- Event Selection -- //
				vector< Muon > SelectedMuonCollection;
                                vector< Int_t > Sel_Index;
                                Int_t IndexDi;
				Bool_t isPassEventSelection = kFALSE;
                                isPassEventSelection = analyzer->EventSelection_Zdiff_13TeV_HighPt(MuonCollection, MuMuInput, &SelectedMuonCollection, &Sel_Index, IndexDi);

				if( isPassEventSelection == kTRUE )
				{
                                        Muon mu1 = SelectedMuonCollection[0];
                                        Muon mu2 = SelectedMuonCollection[1];

                                        MuMu.nVertices = MuMuInput->nVertices;
                                        MuMu.runNum = MuMuInput->runNum;
                                        MuMu.lumiBlock = MuMuInput->lumiBlock;
                                        MuMu.evtNum = MuMuInput->evtNum;
                                        MuMu.nPileUp = MuMuInput->nPileUp;
                                        MuMu.HLT_ntrig = MuMuInput->HLT_ntrig;

                                        MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();

                                        if( IndexDi!=-1 )
                                        {
                                            MuMu.CosAngle->push_back(MuMuInput->CosAngle->at(IndexDi));
                                            MuMu.vtxTrkChi2->push_back(MuMuInput->vtxTrkChi2->at(IndexDi));
                                            MuMu.vtxTrkProb->push_back(MuMuInput->vtxTrkProb->at(IndexDi));
                                            MuMu.vtxTrkNdof->push_back(MuMuInput->vtxTrkNdof->at(IndexDi));
                                            MuMu.vtxTrkCkt1Pt->push_back(MuMuInput->vtxTrkCkt1Pt->at(IndexDi));
                                            MuMu.vtxTrkCkt2Pt->push_back(MuMuInput->vtxTrkCkt2Pt->at(IndexDi));
                                        }
                                        else cout << "== ERROR: Event selection was passed but no dimuon index registered. ==" << endl;

                                        for (UInt_t iter=0; iter<MuMuInput->HLT_trigEta->size(); iter++)
                                        {
                                            if(MuMuInput->HLT_trigEta->at(iter) && MuMuInput->HLT_trigPhi->at(iter))
                                            {
                                                MuMu.HLT_trigFired->push_back(MuMuInput->HLT_trigFired->at(iter));
                                                MuMu.HLT_trigEta->push_back(MuMuInput->HLT_trigEta->at(iter));
                                                MuMu.HLT_trigPhi->push_back(MuMuInput->HLT_trigPhi->at(iter));
                                                if(((Int_t)(iter))<MuMuInput->HLT_ntrig) MuMu.HLT_trigName->push_back(MuMuInput->HLT_trigName->at(iter));
                                            }
                                            else
                                            {
                                                cout << "HLT_trigEta and HLT_trigPhi do not exist. Breaking at " << iter << " iteration." << endl;
                                                break;
                                            }
                                        }

                                        if(Sel_Index.size()!=2) cout << "=========== ERROR: The number of muons saved is not 2 ===========" << endl;
                                        else
                                        {
                                            for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                                            {
                                                Int_t index = Sel_Index[iter];

                                                MuMu.Muon_pT->push_back(MuMuInput->Muon_pT->at(index));
                                                MuMu.Muon_eta->push_back(MuMuInput->Muon_eta->at(index));
                                                MuMu.Muon_phi->push_back(MuMuInput->Muon_phi->at(index));
                                                MuMu.isGLBmuon->push_back(MuMuInput->isGLBmuon->at(index));
                                                MuMu.isPFmuon->push_back(MuMuInput->isPFmuon->at(index));
                                                MuMu.isTRKmuon->push_back(MuMuInput->isTRKmuon->at(index));
                                                MuMu.Muon_charge->push_back(MuMuInput->Muon_charge->at(index));
                                                MuMu.Muon_chi2dof->push_back(MuMuInput->Muon_chi2dof->at(index));
                                                MuMu.Muon_muonHits->push_back(MuMuInput->Muon_muonHits->at(index));
                                                MuMu.Muon_nSegments->push_back(MuMuInput->Muon_nSegments->at(index));
                                                MuMu.Muon_nMatches->push_back(MuMuInput->Muon_nMatches->at(index));
                                                MuMu.Muon_trackerLayers->push_back(MuMuInput->Muon_trackerLayers->at(index));
                                                MuMu.Muon_pixelHits->push_back(MuMuInput->Muon_pixelHits->at(index));
                                                MuMu.Muon_dxyVTX->push_back(MuMuInput->Muon_dxyVTX->at(index));
                                                MuMu.Muon_dzVTX->push_back(MuMuInput->Muon_dzVTX->at(index));
                                                MuMu.Muon_trkiso->push_back(MuMuInput->Muon_trkiso->at(index));
                                                MuMu.Muon_PfChargedHadronIsoR04->push_back(MuMuInput->Muon_PfChargedHadronIsoR04->at(index));
                                                MuMu.Muon_PfNeutralHadronIsoR04->push_back(MuMuInput->Muon_PfNeutralHadronIsoR04->at(index));
                                                MuMu.Muon_PfGammaIsoR04->push_back(MuMuInput->Muon_PfGammaIsoR04->at(index));
                                                MuMu.Muon_PFSumPUIsoR04->push_back(MuMuInput->Muon_PFSumPUIsoR04->at(index));
                                                MuMu.Muon_Px->push_back(MuMuInput->Muon_Px->at(index));
                                                MuMu.Muon_Py->push_back(MuMuInput->Muon_Py->at(index));
                                                MuMu.Muon_Pz->push_back(MuMuInput->Muon_Pz->at(index));

                                                Double_t Mu_E = sqrt( MuMuInput->Muon_Px->at(index)*MuMuInput->Muon_Px->at(index) + MuMuInput->Muon_Py->at(index)
                                                                      *MuMuInput->Muon_Py->at(index)+ MuMuInput->Muon_Pz->at(index)*MuMuInput->Muon_Pz->at(index)
                                                                      + M_Mu*M_Mu );

                                                MuMu.Muon_E->push_back(Mu_E);                                                

                                                MuMu.Muon_Best_pT->push_back(MuMuInput->Muon_Best_pT->at(index));
                                                MuMu.Muon_Best_pTError->push_back(MuMuInput->Muon_Best_pTError->at(index));
                                                MuMu.Muon_Best_Px->push_back(MuMuInput->Muon_Best_Px->at(index));
                                                MuMu.Muon_Best_Py->push_back(MuMuInput->Muon_Best_Py->at(index));
                                                MuMu.Muon_Best_Pz->push_back(MuMuInput->Muon_Best_Pz->at(index));
                                                MuMu.Muon_Best_eta->push_back(MuMuInput->Muon_Best_eta->at(index));
                                                MuMu.Muon_Best_phi->push_back(MuMuInput->Muon_Best_phi->at(index));
                                                MuMu.Muon_Inner_pT->push_back(MuMuInput->Muon_Inner_pT->at(index));
                                                MuMu.Muon_Inner_pTError->push_back(MuMuInput->Muon_Inner_pTError->at(index));
                                                MuMu.Muon_Inner_Px->push_back(MuMuInput->Muon_Inner_Px->at(index));
                                                MuMu.Muon_Inner_Py->push_back(MuMuInput->Muon_Inner_Py->at(index));
                                                MuMu.Muon_Inner_Pz->push_back(MuMuInput->Muon_Inner_Pz->at(index));
                                                MuMu.Muon_Inner_eta->push_back(MuMuInput->Muon_Inner_eta->at(index));
                                                MuMu.Muon_Inner_phi->push_back(MuMuInput->Muon_Inner_phi->at(index));
                                                MuMu.Muon_Outer_pT->push_back(MuMuInput->Muon_Outer_pT->at(index));
                                                MuMu.Muon_Outer_pTError->push_back(MuMuInput->Muon_Outer_pTError->at(index));
                                                MuMu.Muon_Outer_Px->push_back(MuMuInput->Muon_Outer_Px->at(index));
                                                MuMu.Muon_Outer_Py->push_back(MuMuInput->Muon_Outer_Py->at(index));
                                                MuMu.Muon_Outer_Pz->push_back(MuMuInput->Muon_Outer_Pz->at(index));
                                                MuMu.Muon_Outer_eta->push_back(MuMuInput->Muon_Outer_eta->at(index));
                                                MuMu.Muon_Outer_phi->push_back(MuMuInput->Muon_Outer_phi->at(index));
                                                MuMu.Muon_GLB_pT->push_back(MuMuInput->Muon_GLB_pT->at(index));
                                                MuMu.Muon_GLB_pTError->push_back(MuMuInput->Muon_GLB_pTError->at(index));
                                                MuMu.Muon_GLB_Px->push_back(MuMuInput->Muon_GLB_Px->at(index));
                                                MuMu.Muon_GLB_Py->push_back(MuMuInput->Muon_GLB_Py->at(index));
                                                MuMu.Muon_GLB_Pz->push_back(MuMuInput->Muon_GLB_Pz->at(index));
                                                MuMu.Muon_GLB_eta->push_back(MuMuInput->Muon_GLB_eta->at(index));
                                                MuMu.Muon_GLB_phi->push_back(MuMuInput->Muon_GLB_phi->at(index));
                                                MuMu.Muon_TuneP_pT->push_back(MuMuInput->Muon_TuneP_pT->at(index));
                                                MuMu.Muon_TuneP_pTError->push_back(MuMuInput->Muon_TuneP_pTError->at(index));
                                                MuMu.Muon_TuneP_Px->push_back(MuMuInput->Muon_TuneP_Px->at(index));
                                                MuMu.Muon_TuneP_Py->push_back(MuMuInput->Muon_TuneP_Py->at(index));
                                                MuMu.Muon_TuneP_Pz->push_back(MuMuInput->Muon_TuneP_Pz->at(index));
                                                MuMu.Muon_TuneP_eta->push_back(MuMuInput->Muon_TuneP_eta->at(index));
                                                MuMu.Muon_TuneP_phi->push_back(MuMuInput->Muon_TuneP_phi->at(index));
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
                                        MuMu.Muon_E->clear();
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

//					Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
//					Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

				} // End of event selection

			} //End of if( isTriggered )

                        bar.Draw(i);
		} //End of event iteration

		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);                
		Double_t LoopRunTime = looptime.CpuTime();
                cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
                cout << "============================================================\n" << endl;

                cout << "Checking number of entries in output and input files:   ";
                if(chain->GetEntries()==MuonTree->GetEntries()) cout << "OK! (" << MuonTree->GetEntries() << ")" << endl;
                else cout << "Input: " << chain->GetEntries() << "     Output: " << MuonTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

	} //end of i_tup iteration

        MuonFile->cd();
        cout << "Writing into file...";
        MuonTree->Write();
        cout << "Finished." << endl << "Closing a file..." << endl;
        MuonFile->Close();
        if (!MuonFile->IsOpen()) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;

	Double_t TotalRunTime = totaltime.CpuTime();
	cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

	TTimeStamp ts_end;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
