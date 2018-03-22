#include <TChain.h>
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

// -- for Rochester Muon momentum correction -- //
#include <./etc/RoccoR/RoccoR.cc>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


static inline void loadBar(int x, int n, int r, int w);

// -- EMu Method for Background Estimation -- //
void EMuSelection(Int_t type, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
	// -- Run2016 luminosity [/pb] -- //
	Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
	L = L_B2H;

	TString DataType, DataLocation, DataLocation2, Type;

	if( type == 1 ) {
		DataType = "B";
		DataLocation = "SingleMuon_Run2016B";
	}
	else if( type == 2 ) {
		DataType = "C";
		DataLocation = "SingleMuon_Run2016C";
	}
	else if( type == 3 ) {
		DataType = "D";
		DataLocation = "SingleMuon_Run2016D";
	}
	else if( type == 4 ) {
		DataType = "E";
		DataLocation = "SingleMuon_Run2016E";
	}
	else if( type == 5 ) {
		DataType = "F";
		DataLocation = "SingleMuon_Run2016F";
	}
	else if( type == 6 ) {
		DataType = "G";
		DataLocation = "SingleMuon_Run2016G";
	}
	else if( type == 7 ) {
		DataType = "H";
		DataLocation = "SingleMuon_Run2016Hver2";
		DataLocation2 = "SingleMuon_Run2016Hver3";
	}

	if( type < 10  ) Type = "Data";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 23 ) Type = "ttbar_M700toInf";
	else if( type == 31 ) Type = "DYTauTau_M10to50";
	else if( type == 32 ) Type = "DYTauTau_M50toInf";
	else if( type == 41 ) Type = "VVnST";

	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	cout << "Type: " << Type << endl;
	if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

	TString BaseLocation;
	if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
	else BaseLocation = "/data9/DATA/DYntuple/v2.1";
	cout << "DATA location: " << BaseLocation << endl;

	TStopwatch totaltime;
	totaltime.Start();

	DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

	// -- For PU re-weighting -- //
	analyzer->SetupPileUpReWeighting_80X( isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

	// -- For Rochester correction -- //
	TRandom3 *r1 = new TRandom3(0);
	
	// -- For efficiency SF -- //
	analyzer->SetupEfficiencyScaleFactor_BtoF();
	analyzer->SetupEfficiencyScaleFactor_GtoH();
	analyzer->SetupEfficiencyScaleFactor_electron();

	TString OutputDir = "./result";
	TFile *f = new TFile(OutputDir+"/ROOTFile_EMuSelection_"+TString::Itoa(type,10)+"_"+TString::Itoa(isTopPtReweighting,10)+".root", "RECREATE");

	// -- Each ntuple directory & corresponding Tags -- //
	vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
	if( Type == "Data" ) {
		ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
	}
	else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

	//Loop for all samples
	const Int_t Ntup = ntupleDirectory.size();
	for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
	{
		TStopwatch looptime;
		looptime.Start();

		cout << "\t<" << Tag[i_tup] << ">" << endl;

		TChain *chain = new TChain("recoTree/DYTree");
		//Set MC chain
		if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
		//Set Data chain
		else {
			chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
			if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		if( isMC == kTRUE ) {
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			ntuple->TurnOnBranches_GenOthers(); // for quarks
		}
		ntuple->TurnOnBranches_Muon();
		ntuple->TurnOnBranches_Electron();
	
		RoccoR rc("./etc/RoccoR/rcdata.2016.v3");

		// -- Making Histogram -- //
                TH1D *h_emu_mass = new TH1D("h_emu_mass_"+Tag[i_tup], "", 43, massbins);
                TH1D *h_emuSS_mass = new TH1D("h_emuSS_mass_"+Tag[i_tup], "", 43, massbins);
                TH1D *h_emu_mass_fine = new TH1D("h_emu_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_emuSS_mass_fine = new TH1D("h_emuSS_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);

                TH1D *h_mu_pT = new TH1D("h_mu_pT_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_mu_eta = new TH1D("h_mu_eta_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_mu_phi = new TH1D("h_mu_phi_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_muSS_pT = new TH1D("h_muSS_pT_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_muSS_eta = new TH1D("h_muSS_eta_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_muSS_phi = new TH1D("h_muSS_phi_"+Tag[i_tup], "", 100, -5, 5);

                TH1D *h_ele_pT = new TH1D("h_ele_pT_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_ele_eta = new TH1D("h_ele_eta_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_ele_phi = new TH1D("h_ele_phi_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_eleSS_pT = new TH1D("h_eleSS_pT_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_eleSS_eta = new TH1D("h_eleSS_eta_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_eleSS_phi = new TH1D("h_eleSS_phi_"+Tag[i_tup], "", 100, -5, 5);

		Double_t SumWeight = 0, SumWeight_Separated = 0;

		//Int_t NEvents = chain->GetEntries();
		int NEvents = 10000;						// test using few samples
		cout << "\t[Total Events: " << NEvents << "]" << endl;
		for(Int_t i=0; i<NEvents; i++)		
		{	
			loadBar(i+1, NEvents, 100, 100);
		
			ntuple->GetEvent(i);

			/////////////////////////////
			// -- Bring the weights -- //
			/////////////////////////////
			
			// -- Positive/Negative Gen-weights -- //
			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;
			SumWeight += GenWeight;

			// -- Pileup-Reweighting -- //
			Double_t PUWeight = 1;
			if( isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( ntuple->nPileUp );

			// -- efficiency weights -- //
			Double_t weight1 = 0, weight2 = 0, effweight = 0;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			// -- Separate ttbar samples -- //
			Bool_t GenFlag_top = kTRUE;
			//Bool_t GenFlag_top = kFALSE;
			//vector<GenOthers> GenTopCollection;
			//GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

			if( GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				SumWeight_Separated += GenWeight;

				// -- Top Pt Reweighting -- //
				/*if( isTopPtReweighting == 1 && Tag[i_tup].Contains("ttbar") )
				{
					GenOthers t1 = GenTopCollection[0];
					GenOthers t2 = GenTopCollection[1];

					Double_t SF1 = exp(0.0615 - 0.0005*(t1.Pt));
					Double_t SF2 = exp(0.0615 - 0.0005*(t2.Pt));
					GenWeight = GenWeight*sqrt(SF1*SF2);
				}*/
			}

			// -- Normalization -- //
			Double_t TotWeight = GenWeight;
			if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*GenWeight;

			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				///////////////////////////////////////
				// -- Reco level selection : Muon -- //
				///////////////////////////////////////
				vector< Muon > MuonCollection;
				Int_t N_Muons = ntuple->nMuon;
				for(Int_t i_reco=0; i_reco<N_Muons; i_reco++)
				{
					Muon mu;
					mu.FillFromNtuple(ntuple, i_reco);

					////////////////////////////////
					// -- Rochester correction -- //
					////////////////////////////////
					Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
					Int_t s, m;
						
					if( Tag[i_tup] == "Data" )
						SF = rc.kScaleDT(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, s=0, m=0);
					else
						SF = rc.kScaleAndSmearMC(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);

					mu.TuneP_pT = SF*mu.TuneP_pT;
				
					// -- Convert to TuneP variables -- //
					analyzer->ConvertToTunePInfo( mu );

					MuonCollection.push_back( mu );
				}	

				///////////////////////////////////////////
				// -- Reco level selection : Electron -- //
				///////////////////////////////////////////
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
				Bool_t isPassEventSelection = kFALSE;

				isPassEventSelection = analyzer->EventSelection_emu_method_test(MuonCollection, ElectronCollection, ntuple,
												&SelectedMuonCollection, &SelectedElectronCollection);
				if( isPassEventSelection == kTRUE )
				{
					Muon mu = SelectedMuonCollection[0];
					Electron ele = SelectedElectronCollection[0];

					Double_t reco_M = (mu.Momentum + ele.Momentum).M();

					// -- Apply efficiency correcion -- //
					if( isMC == kTRUE )
					{
						weight1 = analyzer->EfficiencySF_EventWeight_emu_BtoF( mu, ele );
						weight2 = analyzer->EfficiencySF_EventWeight_emu_GtoH( mu, ele );
						effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;
					}

					if( mu.charge != ele.charge )
					{
						h_emu_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_emu_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );

						h_mu_pT->Fill( mu.Pt, TotWeight * PUWeight * effweight );
						h_mu_eta->Fill( mu.eta, TotWeight * PUWeight * effweight );
						h_mu_phi->Fill( mu.phi, TotWeight * PUWeight * effweight );

						h_ele_pT->Fill( ele.Pt, TotWeight * PUWeight * effweight );
						h_ele_eta->Fill( ele.eta, TotWeight * PUWeight * effweight );
						h_ele_phi->Fill( ele.phi, TotWeight * PUWeight * effweight );
					}
					else
					{
						h_emuSS_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_emuSS_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );

						h_muSS_pT->Fill( mu.Pt, TotWeight * PUWeight * effweight );
						h_muSS_eta->Fill( mu.eta, TotWeight * PUWeight * effweight );
						h_muSS_phi->Fill( mu.phi, TotWeight * PUWeight * effweight );

						h_eleSS_pT->Fill( ele.Pt, TotWeight * PUWeight * effweight );
						h_eleSS_eta->Fill( ele.eta, TotWeight * PUWeight * effweight );
						h_eleSS_phi->Fill( ele.phi, TotWeight * PUWeight * effweight );
					}
				} // End of event selection

			} //End of if( isTriggered )

		} //End of event iteration

		h_emu_mass->Write();
		h_emu_mass_fine->Write();
		h_emuSS_mass->Write();
		h_emuSS_mass_fine->Write();

		h_mu_pT->Write();
		h_mu_eta->Write();
		h_mu_phi->Write();
		h_muSS_pT->Write();
		h_muSS_eta->Write();
		h_muSS_phi->Write();

		h_ele_pT->Write();
		h_ele_eta->Write();
		h_ele_phi->Write();
		h_eleSS_pT->Write();
		h_eleSS_eta->Write();
		h_eleSS_phi->Write();

		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

		Double_t LoopRunTime = looptime.CpuTime();
		cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

	} //end of i_tup iteration

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

