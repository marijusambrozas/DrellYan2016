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


class HistContainer
{
public:
	TString tag;

	vector<TH2D*> vec_Hist2D;
	TH2D* h2D_x1Q2;
	TH2D* h2D_x2Q2;

	TH2D* h2D_x1Q2_M60to120;
	TH2D* h2D_x2Q2_M60to120;

	HistContainer( TString _tag )
	{
		this->tag = _tag;
		this->Init();
	}

	void Fill( NtupleHandle* ntuple, Double_t M, Double_t weight )
	{
		Double_t Q = ntuple->GENEvt_QScale;
		Double_t Q2 = Q*Q;

		Double_t x1 = ntuple->GENEvt_x1;
		Double_t x2 = ntuple->GENEvt_x2;

		h2D_x1Q2->Fill( x1, Q2, weight );
		h2D_x2Q2->Fill( x2, Q2, weight );

		if( M > 60 && M < 120 )
		{
			h2D_x1Q2_M60to120->Fill( x1, Q2, weight );
			h2D_x2Q2_M60to120->Fill( x2, Q2, weight );
		}
	}

	void Save( TFile *f_output )
	{
		f_output->cd();
		for(const &auto h2D : this->vec_Hist2D )
			h2D->Write();
	}

public:
	void Init()
	{
		// -- make array containing bin edges for log-scale axis -- //
		Double_t minX = 1e-10;
		Double_t maxX = 1;
		Int_t nBinX = 0;
		Double_t *arr_BinEdgeX;
		this->SetBinEdges_LogScale( minX, maxX, nBinX, arr_BinEdgeX );
		for(Int_t i=0; i<nBinX+1; i++)
			printf("[%02d-th X bin edge] %.1e\n", i, arr_BinEdgeX[i]);

		Double_t minY = 1;
		Double_t maxY = 1e10;
		Int_t nBinY = 0;
		Double_t *arr_BinEdgeY;
		this->SetBinEdges_LogScale( minY, maxY, nBinY, arr_BinEdgeY );
		for(Int_t i=0; i<nBinY+1; i++)
			printf("[%02d-th Y bin edge] %.1e\n", i, arr_BinEdgeY[i]);

		// -- inititialization -- //
		this->h2D_x1Q2 = new TH2D("h2D_x1Q2_"+this->tag, "", nBinX, arr_BinEdgeX, nBinY, arr_BinEdgeY );
		this->vec_Hist2D.push_back( h2D_x1Q2 );

		this->h2D_x2Q2 = new TH2D("h2D_x2Q2_"+this->tag, "", nBinX, arr_BinEdgeX, nBinY, arr_BinEdgeY );
		this->vec_Hist2D.push_back( h2D_x2Q2 );

		this->h2D_x1Q2_M60to120 = new TH2D("h2D_x1Q2_M60to120_"+this->tag, "", nBinX, arr_BinEdgeX, nBinY, arr_BinEdgeY );
		this->vec_Hist2D.push_back( h2D_x1Q2_M60to120 );

		this->h2D_x2Q2_M60to120 = new TH2D("h2D_x2Q2_M60to120_"+this->tag, "", nBinX, arr_BinEdgeX, nBinY, arr_BinEdgeY );
		this->vec_Hist2D.push_back( h2D_x2Q2_M60to120 );
	}

	void SetBinEdges_LogScale(Double_t minX, Double_t maxX, Int_t& nBin, Double_t*& arr_BinEdge)
	{
		vector<Double_t> vec_BinEdge;
		Int_t i = 1;
		Double_t order = 1e0;
		while( 1 )
		{
			Double_t binEdge = minX * i * order;
			vec_BinEdge.push_back( binEdge );
			i++;
			if( i % 10 == 0 )
			{
				order = order * i;
				i = 1;
			}

			if( binEdge == maxX ) break;
		}

		nBin = (Int_t)vec_BinEdge.size()-1;

		arr_BinEdge = new Double_t[nBin+1];
		for(Int_t i=0; i<nBin+1; i++)
		{
			arr_BinEdge[i] = vec_BinEdge[i];
		}
	}
};

class HistProducer
{
public:
	TString txtFileName;
	TString tag;
	Int_t isMC;
	Double_t normFactor;

	TString channelType;
	TString triggerName;
	Bool_t exclude_ECALGAP;
	Int_t leptonID;

	HistContainer* Hists_GEN;
	HistContainer* Hists_RECO;

	HistProducer()
	{

	}

	HistProducer( TString _txtFileName, TString _tag, Int_t _isMC ): HistProducer()
	{
		this->txtFileName = _txtFileName;
		this->tag = _tag;
		this->isMC = _isMC;

		// -- histograms -- //
		this->Hists_GEN = new HistContainer( "GEN" );
		this->Hists_RECO = new HistContainer( "RECO" );
	}

	void SetNormFactor( Double_t _value )
	{
		this->normFactor = _value;

		printf("===============[HistProducer]===============\n");
		printf("[isMC = %d] tag = %s, normFactor = %e\n", this->isMC, this->tag.Data(), this->normFactor);
		printf("=================================================\n");
	}

	void SetChannelType( TString _Channel )
	{
		this->channelType = _Channel;
		this->SetupChannel();
	}

	void Save( TFile *f_output )
	{
		this->Hists_GEN->Save( f_output );
		this->Hists_RECO->Save( f_output );
	}

	void Produce()
	{
		TStopwatch totaltime;
		totaltime.Start();

		DYAnalyzer *analyzer = new DYAnalyzer( this->triggerName );
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
				///////////////////////////////
				// -- gen-level selection -- //
				///////////////////////////////
				vector< GenLepton > vec_GenLepton_HP;
				vector< GenLepton > vec_GenLepton_HPFS;
				for(Int_t i_gen=0; i_gen<ntuple->gnpair; i_gen++)
				{
					GenLepton genlep( ntuple, i_gen );

					if( fabs(genlep.ID) == this->leptonID )
					{
						if( genlep.isHardProcess )
							vec_GenLepton_HP.push_back( genlep );

						if( genlep.fromHardProcessFinalState )
							vec_GenLepton_HPFS.push_back( genlep );
					}
				}

				GenPair genpair_HP( vec_GenLepton_HP[0], vec_GenLepton_HP[1] );
				GenPair genpair_HPFS( vec_GenLepton_HPFS[0], vec_GenLepton_HPFS[1] );

				Bool_t Flag_GenPassAcc = analyzer->Flag_PassAcc_Dilepton( 
					genpair_HPFS.First.Pt, genpair_HPFS.First.eta, 
					genpair_HPFS.Second.Pt, genpair_HPFS.Second.eta, this->exclude_ECALGAP );

				if( Flag_GenPassAcc )
				{
					Double_t mass = genpair_HPFS.M;
					if( mass > 15 && mass < 3000 )
						this->Hists_GEN->Fill( ntuple, mass, this->NormFactor*GenWeight ); // -- no PU weights -- //
				}

				////////////////////////////////
				// -- reco-level selection -- //
				////////////////////////////////
				if( ntuple->isTriggered( analyzer->HLT ) )
				{
					Bool_t flag_PassSel = kFALSE;
					Double_t M_reco = -999;

					// -- selections corresponding each channel -- //
					if( this->channelType == "MuMu" )
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

						flag_PassSel = analyzer->EventSelection_MuonChannel( vec_Muon, ntuple, muPair );
						if( flag_PassSel ) M_reco = muPair.M;
					}
					else if( this->channelType == "EE" )
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
						if( flag_PassSel ) M_reco = elecPair.M;

					}
					// -- end of selections corresponding each channel -- //

					if( flag_PassSel )
					{
						Double_t mass = M_reco;
						if( mass > 15 && mass < 3000 )
							this->Hists_RECO->Fill( ntuple, mass, TotWeight );
					}

				} // -- end of if( ntuple->isTriggered( analyzer->HLT ) ) -- //

				// -- dressed leptons -- //
				// vector< GenOthers > vec_GenPhoton1;
				// vector< GenOthers > vec_GenPhoton2;
				// GenLepton genlep_dressed1 = analyzer->DressingLepton( ntuple, vec_GenLepton_HPFS[0], 0.1, &vec_GenPhoton1 );
				// GenLepton genlep_dressed2 = analyzer->DressingLepton( ntuple, vec_GenLepton_HPFS[1], 0.1, &vec_GenPhoton2 );
				// GenPair genpair_dressed( genlep_dressed1, genlep_dressed2 );
			} // -- end of GenFlag -- //

		} // -- end of event iteration -- //

		Double_t TotalRunTime = totaltime.CpuTime();
		cout << "\tTotal RunTime(" << this->Tag << "): " << TotalRunTime << " seconds\n" << endl;

		printf("============================\nProducer() is finished\n============================\n\n");
	}

private:
	void SetupChannel()
	{
		if( this->channelType == "MuMu" )
		{
			this->triggerName = "IsoMu24_OR_IsoTkMu24";
			this->exclude_ECALGAP = kFALSE;
			this->leptonID = 13;
		}
		else if( this->channelType == "EE" )
		{
			this->triggerName = "Ele23_Ele12";
			this->exclude_ECALGAP = kTRUE;
			this->leptonID = 11;
		}
		else
		{
			cout << this->channelType << " is wrong type! ... it should be MuMu or EE!" << endl;
		}
	}
};