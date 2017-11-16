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
#include <Include/RoccoR/RoccoR.cc>

#define nMassBin 28

Double_t arr_MassBinEdge[nMassBin+1] = 
{
	15, 17, 19, 21, 24, 27, 31, 35, 40, 45, 
	51, 58, 65, 73, 82, 91, 100, 112, 124, 137, 
	150, 163, 178, 193, 210, 228, 249, 272, 300
};

class HistContainer
{
public:
	Int_t binSize;

	vector< TH1D* > vec_Hist;
	TH1D* h_mass_preFSR;
	TH1D* h_mass_dressed;
	TH1D* h_mass_postFSR;
	TH1D* h_mass_reco;

	vector< TH2D* > vec_2DHist;
	TH2D* h2D_preFSR_reco;
	TH2D* h2D_dressed_reco;
	TH2D* h2D_postFSR_reco;

	HistContainer()
	{
		this->Init();
	}

	void Save( TFile *f_output )
	{
		f_output->cd();
		for( const auto& h : this->vec_Hist )
			h->Write();

		// -- too large to save -- //
		for( const auto& h : this->vec_2DHist )
			h->Write();
	}

	void Fill( GenPair pair_preFSR, GenPair pair_dressed, GenPair pair_postFSR, Double_t M_reco, Double_t weight )
	{
		Double_t M_preFSR = pair_preFSR.M;
		Double_t M_dressed = pair_dressed.M;
		Double_t M_postFSR = pair_postFSR.M;

		this->h_mass_preFSR->Fill( M_preFSR, weight );
		this->h_mass_dressed->Fill( M_dressed, weight );
		this->h_mass_postFSR->Fill( M_postFSR, weight );
		this->h_mass_reco->Fill( M_reco, weight );

		this->h2D_preFSR_reco->Fill( M_preFSR, M_reco, weight );
		this->h2D_dressed_reco->Fill( M_dressed, M_reco, weight );
		this->h2D_postFSR_reco->Fill( M_postFSR, M_reco, weight );
	}

private:
	void Init()
	{
		this->h_mass_preFSR = new TH1D("h_mass_preFSR", "", nMassBin, arr_MassBinEdge );
		this->vec_Hist.push_back( this->h_mass_preFSR );

		this->h_mass_dressed = new TH1D("h_mass_dressed", "", nMassBin, arr_MassBinEdge );
		this->vec_Hist.push_back( this->h_mass_dressed );

		this->h_mass_postFSR = new TH1D("h_mass_postFSR", "", nMassBin, arr_MassBinEdge );
		this->vec_Hist.push_back( this->h_mass_postFSR );

		this->h_mass_reco = new TH1D("h_mass_reco", "", nMassBin, arr_MassBinEdge );
		this->vec_Hist.push_back( this->h_mass_reco );

		// -- 2D hists -- //
		this->h2D_preFSR_reco = new TH2D("h2D_preFSR_reco", "", nMassBin, arr_MassBinEdge, nMassBin, arr_MassBinEdge );
		this->vec_2DHist.push_back( this->h2D_preFSR_reco );

		this->h2D_dressed_reco = new TH2D("h2D_dressed_reco", "", nMassBin, arr_MassBinEdge, nMassBin, arr_MassBinEdge );
		this->vec_2DHist.push_back( this->h2D_dressed_reco );

		this->h2D_postFSR_reco = new TH2D("h2D_postFSR_reco", "", nMassBin, arr_MassBinEdge, nMassBin, arr_MassBinEdge );
		this->vec_2DHist.push_back( this->h2D_postFSR_reco );
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

		DYAnalyzer *analyzer = new DYAnalyzer( this->TriggerName );
		analyzer->SetupPileUpReWeighting( this->IsMC );

		// -- setup for Rochester muon momentum correction -- //
		TRandom *r1 = new TRandom();
		TRandom *r2 = new TRandom();
		TString IncludePath = gSystem->Getenv("KP_INCLUDE_PATH");
		TString MomCorrPath = TString::Format("%s/RoccoR/rcdata.2016.v3", IncludePath.Data());
		RoccoR rc( MomCorrPath.Data() );

		// -- make chain -- //
		TChain *chain = new TChain("recoTree/DYTree");
		analyzer->MakeTChain_fromTextFile( chain, FileName_ROOTFileList );

		// -- turn-on ntuple -- //		
		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		ntuple->TurnOnBranches_Electron();
		if( this->IsMC )
		{			
			ntuple->TurnOnBranches_GenLepton();
			ntuple->TurnOnBranches_GenOthers();
		}
		ntuple->Ready();

		Int_t nEvent = chain->GetEntries();
		cout << "\t[Total Events: " << nEvent << "]" << endl;

		for(Int_t i=0; i<nEvent; i++)
		{
			loadBar(i+1, nEvent, 100, 100);
			
			ntuple->GetEvent(i);

			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;

			Double_t PUWeight = this->IsMC ? analyzer->PileUpWeightValue(ntuple->nPileUp) : 1.0;

			Double_t TotWeight = this->NormFactor*GenWeight*PUWeight;

			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(this->Tag, ntuple);

			// -- only DY->mumu or DY->ee events according to its tag name -- //
			if( GenFlag == kTRUE )
			{
				////////////////////////////////
				// -- reco-level selection -- //
				////////////////////////////////
				Bool_t flag_PassSel = kFALSE;
				Double_t M_reco = -1;

				if( ntuple->isTriggered( analyzer->HLT ) )
				{
					if( this->ChannelType == "MuMu" )
					{
						MuPair muPair;

						vector< Muon > vec_Muon;
						Int_t nLepton = ntuple->nMuon;
						for(Int_t i_reco=0; i_reco<nLepton; i_reco++)
						{
							Muon mu(ntuple, i_reco);
							// -- Obtain Rochester momentum scale correction -- //
							double u1 = r1->Rndm();
							double u2 = r2->Rndm();
							double SF = 0; int s; int m;
								
							if( !this->IsMC )
								SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
							else
								SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, u1, u2, s=0, m=0);

							// -- Change Muon momentum with updated(corrected) one -- //
							mu.UpdateKinematicVariable_UsingNewPt( SF*mu.Pt );
							vec_Muon.push_back( mu );
						}
						flag_PassSel = analyzer->EventSelection_MuonChannel( vec_Muon, ntuple, muPair);
						if( flag_PassSel )
							M_reco = muPair.M;
					}
					else if( this->ChannelType == "EE" )
					{
						ElecPair elecPair;

						vector< Electron > vec_Electron;
						Int_t nLepton = ntuple->Nelectrons;
						for(Int_t i_reco=0; i_reco<nLepton; i_reco++)
						{
							Electron electron(ntuple, i_reco);
							vec_Electron.push_back( electron );
						}
						flag_PassSel = analyzer->EventSelection_ElectronChannel( vec_Electron, ntuple, elecPair );
						if( flag_PassSel )
							M_reco = elecPair.M;
					}
				}

				// -- only when reco-level event passes full selection -- //
				if( flag_PassSel )
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
							if( genlep.isHardProcess )
								vec_GenLepton_HP.push_back( genlep );

							if( genlep.fromHardProcessFinalState )
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

					// -- dressed leptons -- //
					vector< GenOthers > vec_GenPhoton1;
					vector< GenOthers > vec_GenPhoton2;
					GenLepton genlep_dressed1 = analyzer->DressingLepton( ntuple, vec_GenLepton_HPFS[0], 0.1, &vec_GenPhoton1 );
					GenLepton genlep_dressed2 = analyzer->DressingLepton( ntuple, vec_GenLepton_HPFS[1], 0.1, &vec_GenPhoton2 );
					GenPair genpair_dressed( genlep_dressed1, genlep_dressed2 );

					// Bool_t Flag_GenPassAcc = analyzer->Flag_PassAcc_Dilepton( 
					// 	genpair_HPFS.First.Pt, genpair_HPFS.First.eta, 
					// 	genpair_HPFS.Second.Pt, genpair_HPFS.Second.eta, this->Exclude_ECALGAP );

					// -- finally fill the histograms -- //
					Hists->Fill( genpair_HP, genpair_dressed, genpair_HPFS, M_reco, TotWeight );

				} // -- end of if( flag_PassSel ) -- //
			}

		}

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
};