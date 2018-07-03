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

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"

#define M_Mu 0.1056583715 // -- GeV -- //

static inline void loadBar(int x, int n, int r, int w);

// -- Muon Channel -- //
//void MakeSelectedMuMu(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeSelectedMuMu(TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
	// -- Run2016 luminosity [/pb] -- //
	Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
	L = L_B2H;
        TString OutputName = "SelectedMuMu_ZToMuMu_M4500to6000_4.root";

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
        TTree* MuonTree = new TTree("DYTree", "/recoTree");

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
        Int_t nVertices, runNum, lumiBlock, evtNum, nPileUp;
        Double_t GenWeight;
        Int_t HLT_ntrig;
        vector<int> *HLT_trigFired = new vector<int>;
        vector<string> *HLT_trigName = new vector<string>;
//        vector<double> *HLT_trigPt = new vecotr<double>;
        vector<double> *HLT_trigEta = new vector<double>;
        vector<double> *HLT_trigPhi = new vector<double>;
        vector<double> *Muon_pT = new vector<double>;
        vector<double> *Muon_eta = new vector<double>;
        vector<double> *Muon_phi = new vector<double>;
        vector<int> *isGLBmuon = new vector<int>;
        vector<int> *isPFmuon = new vector<int>;
        vector<int> *isTRKmuon = new vector<int>;
        vector<int> *Muon_charge = new vector<int>;
        vector<double> *Muon_chi2dof = new vector<double>;
        vector<int> *Muon_muonHits = new vector<int>;
        vector<int> *Muon_nSegments = new vector<int>;
        vector<int> *Muon_nMatches = new vector<int>;
        vector<int> *Muon_trackerLayers = new vector<int>;
        vector<int> *Muon_pixelHits = new vector<int>;
        vector<double> *Muon_dxyVTX = new vector<double>;
        vector<double> *Muon_dzVTX = new vector<double>;
        vector<double> *Muon_trkiso = new vector<double>;
        vector<double> *Muon_PfChargedHadronIsoR04 = new vector<double>;
        vector<double> *Muon_PfNeutralHadronIsoR04 = new vector<double>;
        vector<double> *Muon_PfGammaIsoR04 = new vector<double>;
        vector<double> *Muon_PFSumPUIsoR04 = new vector<double>;
        vector<double> *Muon_Px = new vector<double>;
        vector<double> *Muon_Py = new vector<double>;
        vector<double> *Muon_Pz = new vector<double>;
        vector<double> *Muon_E = new vector<double>;
        Double_t Muon_InvM;
        vector<double> *Muon_Best_pT = new vector<double>;
        vector<double> *Muon_Best_pTError = new vector<double>;
        vector<double> *Muon_Best_Px = new vector<double>;
        vector<double> *Muon_Best_Py = new vector<double>;
        vector<double> *Muon_Best_Pz = new vector<double>;
        vector<double> *Muon_Best_eta = new vector<double>;
        vector<double> *Muon_Best_phi = new vector<double>;
        vector<double> *Muon_Inner_pT = new vector<double>;
        vector<double> *Muon_Inner_pTError = new vector<double>;
        vector<double> *Muon_Inner_Px = new vector<double>;
        vector<double> *Muon_Inner_Py = new vector<double>;
        vector<double> *Muon_Inner_Pz = new vector<double>;
        vector<double> *Muon_Inner_eta = new vector<double>;
        vector<double> *Muon_Inner_phi = new vector<double>;
        vector<double> *Muon_Outer_pT = new vector<double>;
        vector<double> *Muon_Outer_pTError = new vector<double>;
        vector<double> *Muon_Outer_Px = new vector<double>;
        vector<double> *Muon_Outer_Py = new vector<double>;
        vector<double> *Muon_Outer_Pz = new vector<double>;
        vector<double> *Muon_Outer_eta = new vector<double>;
        vector<double> *Muon_Outer_phi = new vector<double>;
        vector<double> *Muon_GLB_pT = new vector<double>;
        vector<double> *Muon_GLB_pTError = new vector<double>;
        vector<double> *Muon_GLB_Px = new vector<double>;
        vector<double> *Muon_GLB_Py = new vector<double>;
        vector<double> *Muon_GLB_Pz = new vector<double>;
        vector<double> *Muon_GLB_eta = new vector<double>;
        vector<double> *Muon_GLB_phi = new vector<double>;
        vector<double> *Muon_TuneP_pT = new vector<double>;
        vector<double> *Muon_TuneP_pTError = new vector<double>;
        vector<double> *Muon_TuneP_Px = new vector<double>;
        vector<double> *Muon_TuneP_Py = new vector<double>;
        vector<double> *Muon_TuneP_Pz = new vector<double>;
        vector<double> *Muon_TuneP_eta = new vector<double>;
        vector<double> *Muon_TuneP_phi = new vector<double>;

        MuonTree->Branch("nVertices", &nVertices);
        MuonTree->Branch("runNum", &runNum);
        MuonTree->Branch("lumiBlock", &lumiBlock);
        MuonTree->Branch("evtNum", &evtNum);
        MuonTree->Branch("nPileUp", &nPileUp);
        MuonTree->Branch("GENEvt_weight", &GenWeight);
        MuonTree->Branch("HLT_ntrig", &HLT_ntrig);
        MuonTree->Branch("HLT_trigFired", &HLT_trigFired);
        MuonTree->Branch("HLT_trigName", &HLT_trigName);
//        MuonTree->Branch("HLT_trigPt", &HLT_trigPt);
        MuonTree->Branch("HLT_trigEta", &HLT_trigEta);
        MuonTree->Branch("HLT_trigPhi", &HLT_trigPhi);
        MuonTree->Branch("Muon_pT", &Muon_pT);
        MuonTree->Branch("Muon_eta", &Muon_eta);
        MuonTree->Branch("Muon_phi", &Muon_phi);
        MuonTree->Branch("isGLBmuon", &isGLBmuon);
        MuonTree->Branch("isPFmuon", &isPFmuon);
        MuonTree->Branch("isTRKmuon", &isTRKmuon);
        MuonTree->Branch("Muon_charge", &Muon_charge);
        MuonTree->Branch("Muon_chi2dof", &Muon_chi2dof);
        MuonTree->Branch("Muon_muonHits", &Muon_muonHits);
        MuonTree->Branch("Muon_nSegments", &Muon_nSegments);
        MuonTree->Branch("Muon_nMatches", &Muon_nMatches);
        MuonTree->Branch("Muon_trackerLayers", &Muon_trackerLayers);
        MuonTree->Branch("Muon_pixelHits", &Muon_pixelHits);
        MuonTree->Branch("Muon_dxyVTX", &Muon_dxyVTX);
        MuonTree->Branch("Muon_dzVTX", &Muon_dzVTX);
        MuonTree->Branch("Muon_trkiso", &Muon_trkiso);
        MuonTree->Branch("Muon_PfChargedHadronIsoR04", &Muon_PfChargedHadronIsoR04);
        MuonTree->Branch("Muon_PfNeutralHadronIsoR04", &Muon_PfNeutralHadronIsoR04);
        MuonTree->Branch("Muon_PfGammaIsoR04", &Muon_PfGammaIsoR04);
        MuonTree->Branch("Muon_PFSumPUIsoR04", &Muon_PFSumPUIsoR04);
        MuonTree->Branch("Muon_Px", &Muon_Px);
        MuonTree->Branch("Muon_Py", &Muon_Py);
        MuonTree->Branch("Muon_Pz", &Muon_Pz);
        MuonTree->Branch("Muon_E", &Muon_E);
        MuonTree->Branch("Muon_InvM", &Muon_InvM);
        MuonTree->Branch("Muon_Best_pT", &Muon_Best_pT);
        MuonTree->Branch("Muon_Best_pTError", &Muon_Best_pTError);
        MuonTree->Branch("Muon_Best_Px", &Muon_Best_Px);
        MuonTree->Branch("Muon_Best_Py", &Muon_Best_Py);
        MuonTree->Branch("Muon_Best_Pz", &Muon_Best_Pz);
        MuonTree->Branch("Muon_Best_eta", &Muon_Best_eta);
        MuonTree->Branch("Muon_Best_phi", &Muon_Best_phi);
        MuonTree->Branch("Muon_Inner_pT", &Muon_Inner_pT);
        MuonTree->Branch("Muon_Inner_pTError", &Muon_Inner_pTError);
        MuonTree->Branch("Muon_Inner_Px", &Muon_Inner_Px);
        MuonTree->Branch("Muon_Inner_Py", &Muon_Inner_Py);
        MuonTree->Branch("Muon_Inner_Pz", &Muon_Inner_Pz);
        MuonTree->Branch("Muon_Inner_eta", &Muon_Inner_eta);
        MuonTree->Branch("Muon_Inner_phi", &Muon_Inner_phi);
        MuonTree->Branch("Muon_Outer_pT", &Muon_Outer_pT);
        MuonTree->Branch("Muon_Outer_pTError", &Muon_Outer_pTError);
        MuonTree->Branch("Muon_Outer_Px", &Muon_Outer_Px);
        MuonTree->Branch("Muon_Outer_Py", &Muon_Outer_Py);
        MuonTree->Branch("Muon_Outer_Pz", &Muon_Outer_Pz);
        MuonTree->Branch("Muon_Outer_eta", &Muon_Outer_eta);
        MuonTree->Branch("Muon_Outer_phi", &Muon_Outer_phi);
        MuonTree->Branch("Muon_GLB_pT", &Muon_GLB_pT);
        MuonTree->Branch("Muon_GLB_pTError", &Muon_GLB_pTError);
        MuonTree->Branch("Muon_GLB_Px", &Muon_GLB_Px);
        MuonTree->Branch("Muon_GLB_Py", &Muon_GLB_Py);
        MuonTree->Branch("Muon_GLB_Pz", &Muon_GLB_Pz);
        MuonTree->Branch("Muon_GLB_eta", &Muon_GLB_eta);
        MuonTree->Branch("Muon_GLB_phi", &Muon_GLB_phi);
        MuonTree->Branch("Muon_TuneP_pT", &Muon_TuneP_pT);
        MuonTree->Branch("Muon_TuneP_pTError", &Muon_TuneP_pTError);
        MuonTree->Branch("Muon_TuneP_Px", &Muon_TuneP_Px);
        MuonTree->Branch("Muon_TuneP_Py", &Muon_TuneP_Py);
        MuonTree->Branch("Muon_TuneP_Pz", &Muon_TuneP_Pz);
        MuonTree->Branch("Muon_TuneP_eta", &Muon_TuneP_eta);
        MuonTree->Branch("Muon_TuneP_phi", &Muon_TuneP_phi);

	//Loop for all samples
//	const Int_t Ntup = ntupleDirectory.size();
        const Int_t Ntup = 1;
	for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
	{
		TStopwatch looptime;
		looptime.Start();

//		cout << "\t<" << [i_tup] << ">" << endl;

                TChain *chain = new TChain("recoTree/DYTree");
//		//Set MC chain
//		if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//		//Set Data chain
//		else {
//			chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//			if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//		}
                chain->Add("/media/sf_DATA/ZToMuMu_M4500to6000_4.root/recoTree/DYTree;2"); // NEED A WAY TO TELL THE NUMBER OF CYCLES AND THEIR EXTENTION NAMES

                chain->Add("/media/sf_DATA/ZToMuMu_M4500to6000_4.root/recoTree/DYTree;3");

		NtupleHandle *ntuple = new NtupleHandle( chain );
		if( isMC == kTRUE ) {
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			ntuple->TurnOnBranches_GenOthers(); // for quarks
		}
		ntuple->TurnOnBranches_Muon();

		Double_t SumWeight = 0, SumWeight_Separated = 0;

                Int_t NEvents = chain->GetEntries();
//		Int_t NEvents = 10000;						// test using small events
		cout << "\t[Total Events: " << NEvents << "]" << endl;
		for(Int_t i=0; i<NEvents; i++)		
		{	                       
                        loadBar(i+1, NEvents, 50, 50);
		
			ntuple->GetEvent(i);
			
			// -- Positive/Negative Gen-weights -- //
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;
			SumWeight += GenWeight;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			// -- Separate ttbar samples -- //
			Bool_t GenFlag_top = kTRUE;

			// -- Normalization -- //
			Double_t TotWeight = GenWeight;
			if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*GenWeight;

			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				////////////////////////////////
				// -- Reco level selection -- //
				////////////////////////////////
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
                                vector< Int_t > Sel_Index;
				Bool_t isPassEventSelection = kFALSE;
                                isPassEventSelection = analyzer->EventSelection_Zdiff_13TeV_HighPt(MuonCollection, ntuple, &SelectedMuonCollection, &Sel_Index);

				if( isPassEventSelection == kTRUE )
				{
                                        Muon mu1 = SelectedMuonCollection[0];
                                        Muon mu2 = SelectedMuonCollection[1];

                                        nVertices = ntuple->nVertices;
                                        runNum = ntuple->runNum;
                                        lumiBlock = ntuple->lumiBlock;
                                        evtNum = ntuple->evtNum;
                                        nPileUp = ntuple->nPileUp;
                                        HLT_ntrig = ntuple->HLT_ntrig;

                                        Int_t zero_count = 0; // resolving if there is no more information in arrays

                                        for (Int_t iter=0; iter<1000; iter++)
                                        {
                                            if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
                                            {
                                                HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                                                HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                                                HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                                                if(iter<ntuple->HLT_ntrig) HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                                            }
                                            else
                                            {
                                                zero_count++;
                                                if(zero_count>1) break;
                                            }
                                        }

                                        if(Sel_Index.size()!=2) cout << "========== ERROR: More than 2 muons saved ==========" << endl;
                                        else
                                        {
                                            for (Int_t iter=0; iter<Sel_Index.size(); iter++)
                                            {
                                                Int_t i = Sel_Index[iter];

                                                Muon_pT->push_back(ntuple->Muon_pT[i]);
                                                Muon_eta->push_back(ntuple->Muon_eta[i]);
                                                Muon_phi->push_back(ntuple->Muon_phi[i]);
                                                isGLBmuon->push_back(ntuple->isGLBmuon[i]);
                                                isPFmuon->push_back(ntuple->isPFmuon[i]);
                                                isTRKmuon->push_back(ntuple->isTRKmuon[i]);
                                                Muon_charge->push_back(ntuple->Muon_charge[i]);
                                                Muon_chi2dof->push_back(ntuple->Muon_chi2dof[i]);
                                                Muon_muonHits->push_back(ntuple->Muon_muonHits[i]);
                                                Muon_nSegments->push_back(ntuple->Muon_nSegments[i]);
                                                Muon_nMatches->push_back(ntuple->Muon_nMatches[i]);
                                                Muon_trackerLayers->push_back(ntuple->Muon_trackerLayers[i]);
                                                Muon_pixelHits->push_back(ntuple->Muon_pixelHits[i]);
                                                Muon_dxyVTX->push_back(ntuple->Muon_dxyVTX[i]);
                                                Muon_dzVTX->push_back(ntuple->Muon_dzVTX[i]);
                                                Muon_trkiso->push_back(ntuple->Muon_trkiso[i]);
                                                Muon_PfChargedHadronIsoR04->push_back(ntuple->Muon_PfChargedHadronIsoR04[i]);
                                                Muon_PfNeutralHadronIsoR04->push_back(ntuple->Muon_PfNeutralHadronIsoR04[i]);
                                                Muon_PfGammaIsoR04->push_back(ntuple->Muon_PfGammaIsoR04[i]);
                                                Muon_PFSumPUIsoR04->push_back(ntuple->Muon_PFSumPUIsoR04[i]);
                                                Muon_Px->push_back(ntuple->Muon_Px[i]);
                                                Muon_Py->push_back(ntuple->Muon_Py[i]);
                                                Muon_Pz->push_back(ntuple->Muon_Pz[i]);

                                                Double_t Mu_E = sqrt( ntuple->Muon_Px[i]*ntuple->Muon_Px[i] + ntuple->Muon_Py[i]*ntuple->Muon_Py[i]
                                                                      + ntuple->Muon_Pz[i]*ntuple->Muon_Pz[i] + M_Mu*M_Mu );

                                                Muon_E->push_back(Mu_E);
                                                Muon_InvM = (mu1.Momentum + mu2.Momentum).M();

                                                Muon_Best_pT->push_back(ntuple->Muon_Best_pT[i]);
                                                Muon_Best_pTError->push_back(ntuple->Muon_Best_pTError[i]);
                                                Muon_Best_Px->push_back(ntuple->Muon_Best_Px[i]);
                                                Muon_Best_Py->push_back(ntuple->Muon_Best_Py[i]);
                                                Muon_Best_Pz->push_back(ntuple->Muon_Best_Pz[i]);
                                                Muon_Best_eta->push_back(ntuple->Muon_Best_eta[i]);
                                                Muon_Best_phi->push_back(ntuple->Muon_Best_phi[i]);
                                                Muon_Inner_pT->push_back(ntuple->Muon_Inner_pT[i]);
                                                Muon_Inner_pTError->push_back(ntuple->Muon_Inner_pTError[i]);
                                                Muon_Inner_Px->push_back(ntuple->Muon_Inner_Px[i]);
                                                Muon_Inner_Py->push_back(ntuple->Muon_Inner_Py[i]);
                                                Muon_Inner_Pz->push_back(ntuple->Muon_Inner_Pz[i]);
                                                Muon_Inner_eta->push_back(ntuple->Muon_Inner_eta[i]);
                                                Muon_Inner_phi->push_back(ntuple->Muon_Inner_phi[i]);
                                                Muon_Outer_pT->push_back(ntuple->Muon_Outer_pT[i]);
                                                Muon_Outer_pTError->push_back(ntuple->Muon_Outer_pTError[i]);
                                                Muon_Outer_Px->push_back(ntuple->Muon_Outer_Px[i]);
                                                Muon_Outer_Py->push_back(ntuple->Muon_Outer_Py[i]);
                                                Muon_Outer_Pz->push_back(ntuple->Muon_Outer_Pz[i]);
                                                Muon_Outer_eta->push_back(ntuple->Muon_Outer_eta[i]);
                                                Muon_Outer_phi->push_back(ntuple->Muon_Outer_phi[i]);
                                                Muon_GLB_pT->push_back(ntuple->Muon_GLB_pT[i]);
                                                Muon_GLB_pTError->push_back(ntuple->Muon_GLB_pTError[i]);
                                                Muon_GLB_Px->push_back(ntuple->Muon_GLB_Px[i]);
                                                Muon_GLB_Py->push_back(ntuple->Muon_GLB_Py[i]);
                                                Muon_GLB_Pz->push_back(ntuple->Muon_GLB_Pz[i]);
                                                Muon_GLB_eta->push_back(ntuple->Muon_GLB_eta[i]);
                                                Muon_GLB_phi->push_back(ntuple->Muon_GLB_phi[i]);
                                                Muon_TuneP_pT->push_back(ntuple->Muon_TuneP_pT[i]);
                                                Muon_TuneP_pTError->push_back(ntuple->Muon_TuneP_pTError[i]);
                                                Muon_TuneP_Px->push_back(ntuple->Muon_TuneP_Px[i]);
                                                Muon_TuneP_Py->push_back(ntuple->Muon_TuneP_Py[i]);
                                                Muon_TuneP_Pz->push_back(ntuple->Muon_TuneP_Pz[i]);
                                                Muon_TuneP_eta->push_back(ntuple->Muon_TuneP_eta[i]);
                                                Muon_TuneP_phi->push_back(ntuple->Muon_TuneP_phi[i]);
                                            } // End of vector filling
                                        } // End of else()

                                        MuonTree->Fill();

                                        HLT_trigFired->clear();
                                        HLT_trigName->clear();
                                        HLT_trigEta->clear();
                                        HLT_trigPhi->clear();
                                        Muon_pT->clear();
                                        Muon_eta->clear();
                                        Muon_phi->clear();
                                        isGLBmuon->clear();
                                        isPFmuon->clear();
                                        isTRKmuon->clear();
                                        Muon_charge->clear();
                                        Muon_chi2dof->clear();
                                        Muon_muonHits->clear();
                                        Muon_nSegments->clear();
                                        Muon_nMatches->clear();
                                        Muon_trackerLayers->clear();
                                        Muon_pixelHits->clear();
                                        Muon_dxyVTX->clear();
                                        Muon_dzVTX->clear();
                                        Muon_trkiso->clear();
                                        Muon_PfChargedHadronIsoR04->clear();
                                        Muon_PfNeutralHadronIsoR04->clear();
                                        Muon_PfGammaIsoR04->clear();
                                        Muon_PFSumPUIsoR04->clear();
                                        Muon_Px->clear();
                                        Muon_Py->clear();
                                        Muon_Pz->clear();
                                        Muon_E->clear();
                                        Muon_Best_pT->clear();
                                        Muon_Best_pTError->clear();
                                        Muon_Best_Px->clear();
                                        Muon_Best_Py->clear();
                                        Muon_Best_Pz->clear();
                                        Muon_Best_eta->clear();
                                        Muon_Best_phi->clear();
                                        Muon_Inner_pT->clear();
                                        Muon_Inner_pTError->clear();
                                        Muon_Inner_Px->clear();
                                        Muon_Inner_Py->clear();
                                        Muon_Inner_Pz->clear();
                                        Muon_Inner_eta->clear();
                                        Muon_Inner_phi->clear();
                                        Muon_Outer_pT->clear();
                                        Muon_Outer_pTError->clear();
                                        Muon_Outer_Px->clear();
                                        Muon_Outer_Py->clear();
                                        Muon_Outer_Pz->clear();
                                        Muon_Outer_eta->clear();
                                        Muon_Outer_phi->clear();
                                        Muon_GLB_pT->clear();
                                        Muon_GLB_pTError->clear();
                                        Muon_GLB_Px->clear();
                                        Muon_GLB_Py->clear();
                                        Muon_GLB_Pz->clear();
                                        Muon_GLB_eta->clear();
                                        Muon_GLB_phi->clear();
                                        Muon_TuneP_pT->clear();;
                                        Muon_TuneP_pTError->clear();
                                        Muon_TuneP_Px->clear();
                                        Muon_TuneP_Py->clear();
                                        Muon_TuneP_Pz->clear();
                                        Muon_TuneP_eta->clear();
                                        Muon_TuneP_phi->clear();

//					Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
//					Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

				} // End of event selection

			} //End of if( isTriggered )

		} //End of event iteration

		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

		Double_t LoopRunTime = looptime.CpuTime();
		cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

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

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
    	cout << endl;

    if ( x % (n/r +1) != 0 ) return;

 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";
 
    for (int x=c; x<w; x++) cout << " ";
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
	cout << "]\r" << flush;

}

