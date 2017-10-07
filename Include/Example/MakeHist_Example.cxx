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
#include <TEfficiency.h>

#include <vector>

#include <Include/DYAnalyzer.h>

class HistogramProducer
{
public:
	TString FileName_ROOTFileList;
	TString Tag;
	Int_t IsMC;
	Double_t NormFactor;

	TH1D* h_mass;

	HistogramProducer()
	{
		this->Init();
	}

	HistogramProducer( TString _FileName_ROOTFileList, TString _Tag, Int_t _IsMC ): HistogramProducer()
	{
		this->FileName_ROOTFileList = _FileName_ROOTFileList;
		this->Tag = _Tag;
		this->IsMC = _IsMC;
	}

	void Set_NormFactor( Double_t _value )
	{
		this->NormFactor = _value;

		printf("===============[HistogramProducer]===============\n");
		printf("[IsMC = %d] Tag = %s, NormFactor = %e\n", this->IsMC, this->Tag.Data(), this->NormFactor);
		printf("=================================================\n");
	}

	void Save( TFile *f_output )
	{
		f_output->cd();
		this->h_mass->Write();

		TH1D* h_mass_10GeV = (TH1D*)this->h_mass->Clone("h_mass_10GeV");
		h_mass_10GeV->Rebin(10);
		h_mass_10GeV->Write();

		TH1D* h_mass_100GeV = (TH1D*)this->h_mass->Clone("h_mass_100GeV");
		h_mass_100GeV->Rebin(100);
		h_mass_100GeV->Write();
	}

	void Produce()
	{
		TStopwatch totaltime;
		totaltime.Start();

		DYAnalyzer *analyzer = new DYAnalyzer( "IsoMu24_OR_IsoTkMu24" );

		// -- make chain -- //
		TChain *chain = new TChain("recoTree/DYTree");
		analyzer->MakeTChain_fromTextFile( chain, FileName_ROOTFileList );

		// -- turn-on ntuple -- //		
		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		if( this->IsMC )
			ntuple->TurnOnBranches_GenLepton();
		ntuple->Ready();

		Int_t nEvent = chain->GetEntries();
		cout << "\t[Total Events: " << nEvent << "]" << endl;

		for(Int_t i=0; i<nEvent; i++)
		{
			loadBar(i+1, nEvent, 100, 100);
			
			ntuple->GetEvent(i);

			//Bring weights for NLO MC events
			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;

			Double_t TotWeight = this->NormFactor*GenWeight;

			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(this->Tag, ntuple); // -- for DY->tautau process -- //

			if( GenFlag == kTRUE )
			{
				////////////////////////////////
				// -- reco-level selection -- //
				////////////////////////////////
				if( ntuple->isTriggered( analyzer->HLT ) )
				{
					vector< Muon > MuonCollection;
					Int_t nLepton = ntuple->nMuon;
					for(Int_t i_reco=0; i_reco<nLepton; i_reco++)
					{
						Muon mu(ntuple, i_reco);
						MuonCollection.push_back( mu );
					}

					MuPair SelectedPair = analyzer->EventSelection_MuonChannel( MuonCollection, ntuple );

					if( SelectedPair.Flag_IsNonNull )
					{
						this->h_mass->Fill( SelectedPair.M, TotWeight );
					}
				} // -- end of if( ntuple->isTriggered( analyzer->HLT ) ) -- //
				
			} // -- end of if( GenFlag == kTRUE ) -- //

		} // -- end of event iteration -- //

		Double_t TotalRunTime = totaltime.CpuTime();
		cout << "\tTotal RunTime(" << this->Tag << "): " << TotalRunTime << " seconds\n" << endl;

		printf("============================\nProducer() is finished\n============================\n\n");

	}

private:
	void Init()
	{
		this->FileName_ROOTFileList = "";
		this->Tag = "";
		Int_t IsMC = 0;
		this->NormFactor = 0;

		this->h_mass = new TH1D("h_mass_"+this->Tag, "", 10000, 0, 10000);
	}
};

void MakeHist_Example( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	HistogramProducer *producer = new HistogramProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_Example.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}