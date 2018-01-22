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

const Int_t nBin = 9;
const Double_t arr_MassBinEdge[nBin+1] = 
{
	100, 200, 400, 500, 700, 800, 1000, 1500, 2000, 3000
};

class HistContainerMassBin
{
public:
	Double_t lowerEdge;
	Double_t upperEdge;

	vector<TH1D*> vec_Hist;
	// -- single lepton variables -- //
	TH1D* h_pt_preFSR;
	TH1D* h_eta_preFSR;
	TH1D* h_phi_preFSR;
	TH1D* h_ptLead_preFSR;
	TH1D* h_etaLead_preFSR;
	TH1D* h_phiLead_preFSR;
	TH1D* h_ptSub_preFSR;
	TH1D* h_etaSub_preFSR;
	TH1D* h_phiSub_preFSR;

	// -- di-lepton variables -- //
	TH1D* h_mass_preFSR;
	TH1D* h_diPt_preFSR;
	TH1D* h_diRap_preFSR;

	HistContainerMassBin( Double_t _lowerEdge, Double_t _upperEdge )
	{
		this->lowerEdge = _lowerEdge;
		this->upperEdge = _upperEdge;
		this->Init();
	}

	void Fill(const GenPair& genpair_preFSR, Double_t weight)
	{
		this->h_pt_preFSR->Fill( genpair_preFSR.First.Pt, weight );
		this->h_eta_preFSR->Fill( genpair_preFSR.First.eta, weight );
		this->h_phi_preFSR->Fill( genpair_preFSR.First.phi, weight );
		this->h_pt_preFSR->Fill( genpair_preFSR.Second.Pt, weight );
		this->h_eta_preFSR->Fill( genpair_preFSR.Second.eta, weight );
		this->h_phi_preFSR->Fill( genpair_preFSR.Second.phi, weight );

		this->h_ptLead_preFSR->Fill( genpair_preFSR.First.Pt, weight );
		this->h_etaLead_preFSR->Fill( genpair_preFSR.First.eta, weight );
		this->h_phiLead_preFSR->Fill( genpair_preFSR.First.phi, weight );

		this->h_ptSub_preFSR->Fill( genpair_preFSR.Second.Pt, weight );
		this->h_etaSub_preFSR->Fill( genpair_preFSR.Second.eta, weight );
		this->h_phiSub_preFSR->Fill( genpair_preFSR.Second.phi, weight );

		this->h_mass_preFSR->Fill( genpair_preFSR.M, weight );
		this->h_diPt_preFSR->Fill( genpair_preFSR.Pt, weight );
		this->h_diRap_preFSR->Fill( genpair_preFSR.Rap, weight );
	}

	void Save( TFile *f_output )
	{
		f_output->cd();

		for(const auto& h : vec_Hist )
			h->Write();
	}

private:
	void Init()
	{
		TString TStr_massRange = TString::Format("M%.0lfto%.0lf", this->lowerEdge, this->upperEdge);

		this->h_pt_preFSR = new TH1D("h_pt_preFSR_"+TStr_massRange, "", 10000, 0, 10000);
		vec_Hist.push_back( this->h_pt_preFSR );

		this->h_eta_preFSR = new TH1D("h_eta_preFSR_"+TStr_massRange, "", 4000, -20, 20);
		vec_Hist.push_back( this->h_eta_preFSR );

		this->h_phi_preFSR = new TH1D("h_phi_preFSR_"+TStr_massRange, "", 80, -4, 4);
		vec_Hist.push_back( this->h_phi_preFSR );

		this->h_ptLead_preFSR = new TH1D("h_ptLead_preFSR_"+TStr_massRange, "", 10000, 0, 10000);
		vec_Hist.push_back( this->h_ptLead_preFSR );

		this->h_etaLead_preFSR = new TH1D("h_etaLead_preFSR_"+TStr_massRange, "", 4000, -20, 20);
		vec_Hist.push_back( this->h_etaLead_preFSR );

		this->h_phiLead_preFSR = new TH1D("h_phiLead_preFSR_"+TStr_massRange, "", 80, -4, 4);
		vec_Hist.push_back( this->h_phiLead_preFSR );

		this->h_ptSub_preFSR = new TH1D("h_ptSub_preFSR_"+TStr_massRange, "", 10000, 0, 10000);
		vec_Hist.push_back( this->h_ptSub_preFSR );

		this->h_etaSub_preFSR = new TH1D("h_etaSub_preFSR_"+TStr_massRange, "", 4000, -20, 20);
		vec_Hist.push_back( this->h_etaSub_preFSR );

		this->h_phiSub_preFSR = new TH1D("h_phiSub_preFSR_"+TStr_massRange, "", 80, -4, 4);
		vec_Hist.push_back( this->h_phiSub_preFSR );

		this->h_mass_preFSR = new TH1D("h_mass_preFSR_"+TStr_massRange, "", 10000, 0, 10000);
		vec_Hist.push_back( this->h_mass_preFSR );

		this->h_diPt_preFSR = new TH1D("h_diPt_preFSR_"+TStr_massRange, "", 10000, 0, 10000);
		vec_Hist.push_back( this->h_diPt_preFSR );

		this->h_diRap_preFSR = new TH1D("h_diRap_preFSR_"+TStr_massRange, "", 4000, -20, 20);
		vec_Hist.push_back( this->h_diRap_preFSR );

		for( const auto &hist : vec_Hist ) 	hist->Sumw2();
	}
};

class HistContainer
{
public:
	vector< HistContainerMassBin* > vec_Hist;

	HistContainer()
	{
		this->Init();
	}

	void Fill( const GenPair& genpair_preFSR, Double_t weight )
	{
		Double_t mass = genpair_preFSR.M();

		for(Int_t i=0; i<nBin; i++)
		{
			if( mass > arr_MassBinEdge[i] && mass < arr_MassBinEdge[i+1] )
			{
				vec_Hist[i]->Fill( genpair_preFSR, weight );
				break;
			}
		}
	}

	void Save( TFile *f_output )
	{
		f_output->cd();

		for(const auto& h : vec_Hist )
			h->Save( f_output );
	}

private:
	void Init()
	{
		for(Int_t i=0; i<nBin; i++)
		{
			Double_t lowerEdge = arr_MassBinEdge[i];
			Double_t upperEdge = arr_MassBinEdge[i+1];

			HistContainerMassBin *temp = new HistContainerMassBin( lowerEdge, upperEdge );
			this->vec_Hist.push_back( HistContainerMassBin );
		}
	}
}

class HistProducer
{
public:
	TString FileName_ROOTFileList;
	TString Tag;
	Int_t IsMC;
	Double_t NormFactor;

	TString ChannelType;
	Int_t LeptonID;

	HistContainer* Hists;


	HistProducer()
	{

	}

	HistProducer( TString _FileName_ROOTFileList, TString _Tag, Int_t _IsMC ): HistProducer()
	{
		this->FileName_ROOTFileList = _FileName_ROOTFileList;
		this->Tag = _Tag;
		this->IsMC = _IsMC;

		// -- histograms -- //
		this->Hists = new HistContainer();
	}

	void Set_NormFactor( Double_t _value )
	{
		this->NormFactor = _value;

		printf("===============[HistProducer]===============\n");
		printf("[IsMC = %d] Tag = %s, NormFactor = %e\n", this->IsMC, this->Tag.Data(), this->NormFactor);
		printf("=================================================\n");
	}

	void Set_ChannelType( TString _Channel )
	{
		this->ChannelType = _Channel;
		this->Setup_Channel();
	}

	void Save( TFile *f_output )
	{
		this->Hists->Save( f_output );
	}

	void Produce()
	{
		TStopwatch totaltime;
		totaltime.Start();

		DYAnalyzer *analyzer = new DYAnalyzer( "None" );

		// -- make chain -- //
		TChain *chain = new TChain("recoTree/DYTree");
		analyzer->MakeTChain_fromTextFile( chain, FileName_ROOTFileList );

		// -- turn-on ntuple -- //		
		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		ntuple->TurnOnBranches_Electron();
		if( this->IsMC )
			ntuple->TurnOnBranches_GenLepton();
		ntuple->Ready();

		Int_t nEvent = chain->GetEntries();
		cout << "\t[Total Events: " << nEvent << "]" << endl;

		for(Int_t i=0; i<nEvent; i++)
		{
			loadBar(i+1, nEvent, 100, 100);
			
			ntuple->GetEvent(i);

			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;

			Double_t TotWeight = this->NormFactor*GenWeight;

			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(this->Tag, ntuple);

			// -- only DY->mumu or DY->ee events according to its tag name -- //
			if( GenFlag == kTRUE )
			{
				///////////////////////////////
				// -- gen-level selection -- //
				///////////////////////////////
				vector< GenLepton > vec_GenLepton_HP;
				for(Int_t i_gen=0; i_gen<ntuple->gnpair; i_gen++)
				{
					GenLepton genlep( ntuple, i_gen );

					if( fabs(genlep.ID) == this->LeptonID && genlep.isHardProcess )
					{
						vec_GenLepton_HP.push_back( genlep );
					}
				}

				if( vec_GenLepton_HP.size() != 2 )
				{
					cout << "# gen leptons with isHardProcess is not 2! ... # gen lepton = " << vec_GenLepton_HP.size() << endl;
					break;
				}

				GenPair genpair_HP( vec_GenLepton_HP[0], vec_GenLepton_HP[1] );
				this->Hists->Fill( genpair_HP, TotWeight );
				
			} // -- end of if( GenFlag == kTRUE ) -- //

		} // -- end of event iteration -- //

		Double_t TotalRunTime = totaltime.CpuTime();
		cout << "\tTotal RunTime(" << this->Tag << "): " << TotalRunTime << " seconds\n" << endl;

		printf("============================\nProducer() is finished\n============================\n\n");
	}

private:
	void Setup_Channel()
	{
		if( this->ChannelType == "MuMu" )
		{
			this->LeptonID = 13;
		}
		else if( this->ChannelType == "EE" )
		{
			this->LeptonID = 11;
		}
		else
		{
			cout << this->ChannelType << " is wrong type! ... it should be MuMu or EE!" << endl;
		}
	}
};