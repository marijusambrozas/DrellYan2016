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

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


static inline void loadBar(int x, int n, int r, int w);

// -- Electron Channel -- //
void EESelection(Int_t type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
{
	// -- Run2016 luminosity [/pb] -- //
	Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
	L = L_B2H;

	TString DataType, DataLocation, DataLocation2, Type;

	if( type == 1 ) {
		DataType = "B";
		DataLocation = "DoubleEG_Run2016B";
	}
	if( type == 2 ) {
		DataType = "C";
		DataLocation = "DoubleEG_Run2016C";
	}
	if( type == 3 ) {
		DataType = "D";
		DataLocation = "DoubleEG_Run2016D";
	}
	if( type == 4 ) {
		DataType = "E";
		DataLocation = "DoubleEG_Run2016E";
	}
	if( type == 5 ) {
		DataType = "F";
		DataLocation = "DoubleEG_Run2016F";
	}
	if( type == 6 ) {
		DataType = "G";
		DataLocation = "DoubleEG_Run2016G";
	}
	if( type == 7 ) {
		DataType = "H";
		DataLocation = "DoubleEG_Run2016Hver2";
		DataLocation2 = "DoubleEG_Run2016Hver3";
	}

	Bool_t isMC = kTRUE;
	if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
	// -- Signal MC samples -- //
	if( type == 11 ) Type = "DYEE_M10to50";
	if( type == 12 ) Type = "DYEE_M50to100";
	if( type == 13 ) Type = "DYEE_M100toInf";
	// -- Background MC samples -- //
	if( type == 21 ) Type = "ttbar";
	if( type == 22 ) Type = "ttbarBackup";
	if( type == 23 ) Type = "ttbar_M700toInf";
	if( type == 31 ) Type = "DYTauTau_M10to50";
	if( type == 32 ) Type = "DYTauTau_M50toInf";
	if( type == 41 ) Type = "VVnST";
	if( type == 51 ) Type = "WJetsToLNu";
	//if( type == 61 ) Type = "QCDEMEnriched";

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

	// -- For efficiency SF -- //
	analyzer->SetupEfficiencyScaleFactor_electron();

	// -- Output ROOTFile -- //	
	TString OutputDir = "./result";
	TFile *f = new TFile(OutputDir+"/ROOTFile_EESelection_"+TString::Itoa(type,10)+"_"+TString::Itoa(isTopPtReweighting,10)+".root", "RECREATE");

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
			if( type == 7 ) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		if( isMC == kTRUE ) {
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			ntuple->TurnOnBranches_GenOthers(); // for quarks
		}
		ntuple->TurnOnBranches_Electron();
	
		// -- Making Histogram -- //
                TH1D *h_mass_fine_before_PUCorr = new TH1D("h_mass_fine_before_PUCorr_"+Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass_fine_before_EffCorr = new TH1D("h_mass_fine_before_EffCorr_"+Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass = new TH1D("h_mass_"+Tag[i_tup], "", 43, massbins);
                TH1D *h_Pt = new TH1D("h_Pt_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_rapi = new TH1D("h_rapi_"+Tag[i_tup], "", 100, -5, 5);

                TH1D *h_pT = new TH1D("h_pT_"+Tag[i_tup], "", 300, 0, 600);
                TH1D *h_eta = new TH1D("h_eta_"+Tag[i_tup], "", 100, -5, 5);
                TH1D *h_phi = new TH1D("h_phi_"+Tag[i_tup], "", 100, -5, 5);

		Double_t SumWeight = 0, SumWeight_Separated = 0;

		//Int_t NEvents = chain->GetEntries();
		Int_t NEvents = 10000;					// test using small events
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
			Double_t effweight = 1;

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
				////////////////////////////////
				// -- Reco level selection -- //
				////////////////////////////////
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
				Bool_t isPassEventSelection = kFALSE;
				isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection);

				if( isPassEventSelection == kTRUE )
				{
					Electron ele1 = SelectedElectronCollection[0];
					Electron ele2 = SelectedElectronCollection[1];

					Double_t reco_M = (ele1.Momentum + ele2.Momentum).M();
					Double_t reco_Pt = (ele1.Momentum + ele2.Momentum).Pt();
					Double_t reco_rapi = (ele1.Momentum + ele2.Momentum).Rapidity();

					// -- Apply efficiency correcion -- //
					if( isMC == kTRUE )
						effweight = analyzer->EfficiencySF_EventWeight_electron( ele1, ele2 );

					h_mass_fine_before_PUCorr->Fill( reco_M, TotWeight );
					h_mass_fine_before_EffCorr->Fill( reco_M, TotWeight * PUWeight );
					h_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
					h_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
					h_Pt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
					h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

					h_pT->Fill( ele1.Pt, TotWeight * PUWeight * effweight ); h_pT->Fill( ele2.Pt, TotWeight * PUWeight * effweight );
					h_eta->Fill( ele1.eta, TotWeight * PUWeight * effweight ); h_eta->Fill( ele2.eta, TotWeight * PUWeight * effweight );
					h_phi->Fill( ele1.phi, TotWeight * PUWeight * effweight ); h_phi->Fill( ele2.phi, TotWeight * PUWeight * effweight );

				} // End of event selection

			} //End of if( isTriggered )

		} //End of event iteration

		h_mass_fine_before_PUCorr->Write();
		h_mass_fine_before_EffCorr->Write();
		h_mass_fine->Write();
		h_mass->Write();
		h_Pt->Write();
		h_rapi->Write();

		h_pT->Write();
		h_eta->Write();
		h_phi->Write();


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

