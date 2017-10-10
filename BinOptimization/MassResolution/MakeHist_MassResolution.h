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

// -- binning for calculating mass resolution -- //
const Int_t nBin = 34;
const Int_t MassBinEdges[nBin+1] =
{
	10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 
	60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 
	110, 115, 120, 200, 300, 400, 500, 600, 700, 800, 
	900, 1000, 1500, 2000, 3000
};

class HistContainer
{
public:
	TH1D* h_RelDiff_All;
	vector<TH1D*> vec_Hist_MassBinned;

	HistContainer()
	{
		this->Init();
	}

	void Fill()
	{

	}

private:
	void Init()
	{
		this->h_RelDiff_All = new TH1D("h_RelDiff_All", "", 10000, -5, 5);

		for(Int_t i=0; i<nBin; i++)
		{
			TString TStr_MassRange = TString::Format("M%dto%d", MassBinEdges[i], MassBinEdges[i+1]);
			TH1D* h_temp = new TH1D("h_RelDiff_"+TStr_MassRange, "", 10000, -5, 5);
			this->vec_Hist_MassBinned.push_bacK( h_temp );
		}
	}

	void Save( TFile *f_output )
	{
		f_output->cd();

		h_RelDiff_All->Write();
		for(const auto& h : vec_Hist_MassBinned )
			h->Write();
	}
};

class HistProducer
{
public:
	TString FileName_ROOTFileList;
	TString Tag;
	Int_t IsMC;
	Double_t NormFactor;

	TString ChannelType;
	TString TriggerName;
	Bool_t Exclude_ECALGAP;
	Int_t LeptonID;


	HistProducer()
	{

	}

	HistProducer( TString _FileName_ROOTFileList, TString _Tag, Int_t _IsMC ): HistProducer()
	{
		this->FileName_ROOTFileList = _FileName_ROOTFileList;
		this->Tag = _Tag;
		this->IsMC = _IsMC;
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
		f_output->cd();
	}

	void Produce()
	{
		TStopwatch totaltime;
		totaltime.Start();

		DYAnalyzer *analyzer = new DYAnalyzer( this->TriggerName );

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
				vector< GenLepton > vec_GenLepton_HPFS;
				for(Int_t i_gen=0; i_gen<ntuple->gnpair; i_gen++)
				{
					GenLepton genlep( ntuple, i_gen );

					if( fabs(genlep.ID) == this->LeptonID )
					{
						if( GenLepton.isHardProcess )
							vec_GenLepton_HP.push_back( genlep );

						if( GenLepton.fromHardProcessFinalState )
							vec_GenLepton_HPFS.push_back( genlep );
					}
				}

				if( vec_GenLepton_HP.size() != 2 )
				{
					cout << "# gen leptons with isHardProcess is not 2! ... # gen lepton = " << vec_GenLepton_HP.size() << endl;
					break;
				}

				if( vec_GenLepton_HPFS.size() != 2 )
				{
					cout << "# gen leptons with fromHardProcessFinalState is not 2! ... # gen lepton = " << vec_GenLepton_HPFS.size() << endl;
					break;
				}

				GenPair genpair_HP( vec_GenLepton_HP[0], vec_GenLepton_HP[1] );
				GenPair genpair_HPFS( vec_GenLepton_HPFS[0], vec_GenLepton_HPFS[1] );

				Bool_t Flag_GenPassAcc = analyzer->Flag_PassAcc_Dilepton( 
					genpair_HPFS.First.Pt, genpair_HPFS.First.eta, 
					genpair_HPFS.Second.Pt, genpair_HPFS.Second.eta, this->Exclude_ECALGAP );

				if( Flag_GenPassAcc )
				{
					////////////////////////////////
					// -- reco-level selection -- //
					////////////////////////////////
					if( ntuple->isTriggered( analyzer->HLT ) )
					{
						if( this->ChannelType == "MuMu" )
						{
							MuPair mupair;
							Bool_t Flag_Matched = this->GenRecoMatching_MuPair( ntuple, analyzer, genpair_HPFS, mupair );

							if( Flag_Matched )
							{
								PrintOutDoubleMuInfo( mupair.First, mupair.Second );
								// Hists->Fill( genpair_HP, genpair_HPFS, mupair );
							}
						}
						else if( this->ChannelType == "EE" )
						{
							// ElecPair elecpair;
							// Bool_t Flag_Matched = this->GenRecoMatching_ElecPair( ntuple, analyzer, genpair_HPFS, elecpair );
							// Hists->Fill( genpair_HP, genpair_HPFS, elecpair );
						}

					} // -- end of if( ntuple->isTriggered( analyzer->HLT ) ) -- //

				} // -- end of if( Flag_GenPassAcc ) -- //
				
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
			this->TriggerName = "IsoMu24_OR_IsoTkMu24";
			this->Exclude_ECALGAP = kFALSE;
			this->LeptonID = 13;
		}
		else if( this->ChannelType == "EE" )
		{
			this->TriggerName = "Ele23_Ele12";
			this->Exclude_ECALGAP = kTRUE;
			this->LeptonID = 11;
		}
		else
		{
			cout << this->ChannelType << " is wrong type! ... it should be MuMu or EE!" << endl;
		}
	}

	Bool_t GenRecoMatching_MuPair( NtupleHandle* ntuple, DYAnalyzer *analyzer, GenPair genpair, MuPair &mupair )
	{
		Bool_t Flag_Pass = kFALSE;

		vector< GenLepton > vec_GenLepton;
		vec_GenLepton.push_back( genpair.First );
		vec_GenLepton.push_back( genpair.Second );

		vector< Muon > vec_MatchedRecoLepton;
		// -- for each gen-lepton, find reco lepton that has smallest dR -- //
		for( const auto& genlep : vec_GenLepton )
		{
			Double_t dR_Smallest = 999;
			Muon Lepton_SmallestdR;

			for(Int_t i_reco=0; i_reco<ntuple->nMuon; i_reco++)
			{
				Muon lepton(ntuple, i_reco);
				Double_t dR = genlep.Momentum.DeltaR( lepton.Momentum );

				if( dR < dR_Smallest )
				{
					dR_Smallest = dR;
					Lepton_SmallestdR = lepton;
				}

			} // -- iteration over all leptons -- //

			// -- at least this dR should be less than 0.3 -- //
			if( dR_Smallest < 0.3 )
				vec_MatchedRecoLepton.push_back( Lepton_SmallestdR );
		}

		// -- all matched leptons are found -- //
		if( vec_MatchedRecoLepton.size() == 2 )
		{
			Double_t LeadPt, LeadEta, SubPt, SubEta;
			if( vec_MatchedRecoLepton[0].Pt > vec_MatchedRecoLepton[1].Pt )
			{
				LeadPt = vec_MatchedRecoLepton[0].Pt;
				LeadEta = vec_MatchedRecoLepton[0].eta;
				SubPt = vec_MatchedRecoLepton[1].Pt;
				SubEta = vec_MatchedRecoLepton[1].eta;
			}
			else
			{
				LeadPt = vec_MatchedRecoLepton[1].Pt;
				LeadEta = vec_MatchedRecoLepton[1].eta;
				SubPt = vec_MatchedRecoLepton[0].Pt;
				SubEta = vec_MatchedRecoLepton[0].eta;
			}
			Bool_t Flag_PassAcc = analyzer->Flag_PassAcc_Dilepton( LeadPt, LeadEta, SubPt, SubEta, this->Exclude_ECALGAP );

			if( Flag_PassAcc )
			{
				Flag_Pass = kTRUE;
				mupair.Set( vec_MatchedRecoLepton[0], vec_MatchedRecoLepton[1] );
			}
		}

		return Flag_Pass;
	}

	Bool_t GenRecoMatching_ElecPair( NtupleHandle* ntuple, DYAnalyzer *analyzer, GenPair genpair, ElecPair &elecpair )
	{

	}

};