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

class HistProducer
{
public:
	TString FileName_txt;
	TString Tag;
	Int_t IsMC;

	Double_t NormFactor;

	TH1D* h_SumWeight;

	HistProducer(TString _FileName_ROOTFileList, TString _Tag, Int_t _IsMC )
	{
		this->FileName_txt = _FileName_ROOTFileList;
		this->Tag = _Tag;
		this->IsMC = _IsMC;

		this->h_SumWeight = new TH1D("h_SumWeight", "", 1, 0, 1);
	}

	void Set_NormFactor( Double_t _NormFactor )
	{
		this->NormFactor = _NormFactor;
	}

	void Save( TFile *f_output )
	{
		f_output->cd();
		this->h_SumWeight->Write();
	}

	void Produce()
	{
		TStopwatch totaltime;
		totaltime.Start();

		DYAnalyzer *analyzer = new DYAnalyzer( "None" );

		// -- make chain -- //
		TChain *chain = new TChain("recoTree/DYTree");
		analyzer->MakeTChain_fromTextFile( chain, this->FileName_txt );

		// -- turn-on ntuple -- //		
		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		if( this->IsMC )
			ntuple->TurnOnBranches_GenLepton();
		ntuple->Ready();

		Int_t nEvent = chain->GetEntries();
		cout << "\t[Total Events: " << nEvent << "]" << endl;

		Double_t SumWeight = 0;

		for(Int_t i=0; i<nEvent; i++)
		{
			loadBar(i+1, nEvent, 100, 100);
			
			ntuple->GetEvent(i);

			//Bring weights for NLO MC events
			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;

			Double_t TotWeight = this->NormFactor*GenWeight;

			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(this->Tag, ntuple);

			if( GenFlag	)
			{
				SumWeight += GenWeight;
			}

		} // -- end of event iteration -- //

		this->h_SumWeight->SetBinContent(1, SumWeight);
		this->h_SumWeight->SetBinError(1, 0);
	}

};



void Calc_SumWeight( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_Calc_SumWeight.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}