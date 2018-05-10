#include <Include/PlotTools.h>

// const Int_t nBin = 9;
// const Double_t arr_MassBinEdge[nBin+1] = 
// {
// 	100, 200, 400, 500, 700, 800, 1000, 1500, 2000, 3000
// };

const Int_t nBin = 5;
const Double_t arr_MassBinEdge[nBin+1] = 
{
	100, 200, 400, 500, 700, 800
};

class HistDrawer
{
public:
	TString channelType;
	SampleInfo *siInc;
	SampleInfo *siBinned;

	vector< HistInfo* > vec_Hist;

	HistDrawer( TString _channelType )
	{
		this->channelType = _channelType;
		this->Init();
	}

	void DrawAll()
	{
		for(Int_t i=0; i<nBin; i++)
		{
			Double_t lowerEdge = arr_MassBinEdge[i];
			Double_t upperEdge = arr_MassBinEdge[i+1];
			for( const auto& hist : this->vec_Hist )
			{
				this->Draw( hist, lowerEdge, upperEdge );
			}
		}


	}
private:
	void Init()
	{
		// -- histogram information to be drawn -- //
		HistInfo *histMass10GeV = new HistInfo( "h_mass_preFSR", "m [GeV]", "Entries / 10 GeV");
		histMass10GeV->SetFullName( "DY"+this->channelType+"_Mass10GeV");
		histMass10GeV->SetRebinX( 10 );
		// histMass10GeV->SetXRange( 50, 3000 );
		this->vec_Hist.push_back( histMass10GeV );

		HistInfo *histMass20GeV = new HistInfo( "h_mass_preFSR", "m [GeV]", "Entries / 20 GeV");
		histMass20GeV->SetFullName( "DY"+this->channelType+"_Mass20GeV");
		histMass20GeV->SetRebinX( 20 );
		// histMass20GeV->SetXRange( 50, 3000 );
		this->vec_Hist.push_back( histMass20GeV );

		HistInfo *histDiRap = new HistInfo( "h_diRap_preFSR", "y(ll)", "Entries / 0.1");
		histDiRap->SetFullName( "DY"+this->channelType+"_DiRap");
		histDiRap->SetRebinX( 10 );
		histDiRap->SetXRange( -5, 5 );
		this->vec_Hist.push_back( histDiRap );

		HistInfo *histDiPt = new HistInfo( "h_diPt_preFSR", "P_{T}(ll) [GeV]", "Entries / 10 GeV");
		histDiPt->SetFullName( "DY"+this->channelType+"_DiPt");
		histDiPt->SetRebinX( 10 );
		histDiPt->SetXRange( 0, 250 );
		this->vec_Hist.push_back( histDiPt );

		HistInfo *histPt = new HistInfo( "h_pt_preFSR", "P_{T}(l) [GeV]", "Entries / 10 GeV");
		histPt->SetFullName( "DY"+this->channelType+"_Pt");
		histPt->SetRebinX( 10 );
		histPt->SetXRange( 0, 500 );
		this->vec_Hist.push_back( histPt );

		HistInfo *histPtLead = new HistInfo( "h_ptLead_preFSR", "P_{T}(l, leading) [GeV]", "Entries / 10 GeV");
		histPtLead->SetFullName( "DY"+this->channelType+"_PtLead");
		histPtLead->SetRebinX( 10 );
		histPtLead->SetXRange( 0, 500 );
		this->vec_Hist.push_back( histPtLead );

		HistInfo *histPtSub = new HistInfo( "h_ptSub_preFSR", "P_{T}(l, sub-leading) [GeV]", "Entries / 10 GeV");
		histPtSub->SetFullName( "DY"+this->channelType+"_PtSub");
		histPtSub->SetRebinX( 10 );
		histPtSub->SetXRange( 0, 500 );
		this->vec_Hist.push_back( histPtSub );

		HistInfo *histEta = new HistInfo( "h_eta_preFSR", "#eta(l)", "Entries / 0.1");
		histEta->SetFullName( "DY"+this->channelType+"_Eta");
		histEta->SetRebinX( 10 );
		histEta->SetXRange( -10, 10 );
		this->vec_Hist.push_back( histEta );

		HistInfo *histEtaLead = new HistInfo( "h_etaLead_preFSR", "#eta(l, leading)", "Entries / 0.1");
		histEtaLead->SetFullName( "DY"+this->channelType+"_EtaLead");
		histEtaLead->SetRebinX( 10 );
		histEtaLead->SetXRange( -10, 10 );
		this->vec_Hist.push_back( histEtaLead );

		HistInfo *histEtaSub = new HistInfo( "h_etaSub_preFSR", "#eta(l, sub-leading)", "Entries / 0.1");
		histEtaSub->SetFullName( "DY"+this->channelType+"_EtaSub");
		histEtaSub->SetRebinX( 10 );
		histEtaSub->SetXRange( -10, 10 );
		this->vec_Hist.push_back( histEtaSub );

		HistInfo *histPhi = new HistInfo( "h_phi_preFSR", "#phi(l)", "Entries / 0.1");
		histPhi->SetFullName( "DY"+this->channelType+"_Phi");
		this->vec_Hist.push_back( histPhi );

		HistInfo *histPhiLead = new HistInfo( "h_phiLead_preFSR", "#phi(l, leading)", "Entries / 0.1");
		histPhiLead->SetFullName( "DY"+this->channelType+"_PhiLead");
		this->vec_Hist.push_back( histPhiLead );

		HistInfo *histPhiSub = new HistInfo( "h_phiSub_preFSR", "#phi(l, sub-leading)", "Entries / 0.1");
		histPhiSub->SetFullName( "DY"+this->channelType+"_PhiSub");
		this->vec_Hist.push_back( histPhiSub );

		HistInfo *histDiPt_0jet = new HistInfo( "h_diPt_0jet_preFSR", "P_{T}(ll) [GeV] (0-jet)", "Entries / 10 GeV");
		histDiPt_0jet->SetFullName( "DY"+this->channelType+"_DiPt_0jet");
		histDiPt_0jet->SetRebinX( 10 );
		histDiPt_0jet->SetXRange( 0, 250 );
		this->vec_Hist.push_back( histDiPt_0jet );

		HistInfo *histDiPt_1jet = new HistInfo( "h_diPt_1jet_preFSR", "P_{T}(ll) [GeV] (1-jet)", "Entries / 10 GeV");
		histDiPt_1jet->SetFullName( "DY"+this->channelType+"_DiPt_1jet");
		histDiPt_1jet->SetRebinX( 10 );
		histDiPt_1jet->SetXRange( 0, 250 );
		this->vec_Hist.push_back( histDiPt_1jet );

		HistInfo *histDiPt_2jet = new HistInfo( "h_diPt_2jet_preFSR", "P_{T}(ll) [GeV] (2-jet)", "Entries / 10 GeV");
		histDiPt_2jet->SetFullName( "DY"+this->channelType+"_DiPt_2jet");
		histDiPt_2jet->SetRebinX( 10 );
		histDiPt_2jet->SetXRange( 0, 250 );
		this->vec_Hist.push_back( histDiPt_2jet );

	}

	void Draw( HistInfo* hist, Double_t lowerEdge, Double_t upperEdge )
	{
		// -- setup SampleInfo -- //
		TString fileNameInc = TString::Format("Local/ROOTFile_MakeHist_Inclusive_%s.root", this->channelType.Data());
		SampleInfo *siInc = new SampleInfo( 0, "Inclusive", "Inclusive (M50-inf)" );
		siInc->SetFileName( fileNameInc );
		siInc->SetColor( kRed );

		TString fileNameBinned = TString::Format("Local/ROOTFile_MakeHist_GenLevel_%s_MassBinned.root", this->channelType.Data());
		SampleInfo *siBinned = new SampleInfo( 0, "Mass-binned", "Mass-binned" );
		siBinned->SetFileName( fileNameBinned );
		siBinned->SetColor( kBlack );

		// -- setup histInfo -- //
		if( hist->name.Contains("mass") ) hist->SetXRange( lowerEdge, upperEdge );

		Double_t ratioMin = 0.9;
		Double_t ratioMax = 1.1;
		this->SetRange_EachMassRange( hist, lowerEdge, upperEdge, ratioMin, ratioMax );

		// -- get histogram -- //
		TString histNameInc = TString::Format("DY%s/%s_M%.0lfto%.0lf", this->channelType.Data(), hist->name.Data(), lowerEdge, upperEdge);
		TH1D* hInc = Get_Hist( siInc->fileName, histNameInc );

		TString histNameBinned = TString::Format("DY%s_M%.0lfto%.0lf/%s", this->channelType.Data(), lowerEdge, upperEdge, hist->name.Data() );
		TH1D* hBinned = Get_Hist( siBinned->fileName, histNameBinned );

		// -- setup TH1Ext -- //
		TH1Ext *hextInc = new TH1Ext( siInc, hist, hInc );
		hextInc->CalcRatio_DEN( hextInc->h ); // -- dummy -- //
		TH1Ext *hextBinned = new TH1Ext( siBinned, hist, hBinned );
		hextBinned->CalcRatio_DEN( hextInc->h ); // -- dummy -- //

		TString CanvasName = "Local/c" + TString::Format("_M%.0lfto%.0lf_", lowerEdge, upperEdge) + hist->fullName;
		TCanvas *c; TPad *TopPad; TPad *BottomPad;
		SetCanvas_Ratio( c, CanvasName, TopPad, BottomPad, 0, 1 );
		c->cd();
		TopPad->cd();

		hextInc->DrawAndSet( "EPSAME" );
		hextBinned->DrawAndSet( "EPSAME" );

		TLegend *legend;
		SetLegend( legend, 0.50, 0.82, 0.97, 0.95 );
		hextInc->AddToLegend( legend );
		hextBinned->AddToLegend( legend );
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );

		TString channelInfo = "";
		if( this->channelType == "MuMu" ) channelInfo = "#mu channel";
		if( this->channelType == "EE" ) channelInfo = "e channel";
		latex.DrawLatexNDC(0.16, 0.91, "#scale[0.8]{#font[62]{"+channelInfo+"}}");

		TString massRangeInfo = TString::Format("%.0lf < M < %.0lf GeV", lowerEdge, upperEdge);
		latex.DrawLatexNDC(0.16, 0.875, "#scale[0.6]{#font[42]{"+massRangeInfo+"}}");


		c->cd();
		BottomPad->cd();
		hextBinned->DrawRatioAndSet( "EPSAME", "Binned/Inc.", ratioMin, ratioMax );
		hextBinned->h_ratio->GetYaxis()->SetNdivisions(006);

		TF1 *f_line;
		DrawLine( f_line );

		c->SaveAs(".pdf");
	}

	void SetRange_EachMassRange( HistInfo* hist, Double_t lowerEdge, Double_t upperEdge, Double_t &ratioMin, Double_t &ratioMax )
	{
		if( lowerEdge == 100 && upperEdge == 200 )
		{
			if( hist->name.Contains("pt") )
			{
				hist->SetYRange(1, 1e7);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("eta") ) hist->SetYRange(1, 1e6);
			if( hist->name.Contains("mass") ) hist->SetYRange(1e4, 1e7);
			if( hist->name.Contains("diPt") )
			{
				hist->SetYRange(1e2, 1e7);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("diRap") ) hist->SetYRange(1, 1e6);
		}
		else if( lowerEdge == 200 && upperEdge == 400 )
		{
			if( hist->name.Contains("pt") )
			{
				hist->SetYRange(1, 1e5);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("eta") ) hist->SetYRange(1, 1e4);
			if( hist->name.Contains("mass") ) hist->SetYRange(1e2, 1e6);
			if( hist->name.Contains("diPt") )
			{
				hist->SetYRange(1e2, 1e6);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("diRap") ) hist->SetYRange(1, 1e4);
		}
		else if( lowerEdge == 400 && upperEdge == 500 )
		{
			if( hist->name.Contains("pt") )
			{
				hist->SetYRange(1, 1e4);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("eta") ) hist->SetYRange(1, 1e3);
			if( hist->name.Contains("mass") ) hist->SetYRange(1e2, 1e4);
			if( hist->name.Contains("diPt") )
			{
				hist->SetYRange(1, 1e5);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("diRap") ) hist->SetYRange(1, 1e3);
		}
		else if( lowerEdge == 500 && upperEdge == 700 )
		{
			if( hist->name.Contains("pt") )
			{
				hist->SetYRange(1, 1e4);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("eta") ) hist->SetYRange(1, 1e3);
			if( hist->name.Contains("mass") ) hist->SetYRange(10, 1e4);
			if( hist->name.Contains("diPt") )
			{
				hist->SetYRange(1, 1e5);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("diRap") ) hist->SetYRange(1, 1e3);
		}
		else if( lowerEdge == 700 && upperEdge == 800 )
		{
			if( hist->name.Contains("pt") )
			{
				hist->SetYRange(1e-1, 1e3);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("eta") ) hist->SetYRange(1, 1e2);
			if( hist->name.Contains("mass") ) hist->SetYRange(1, 1e3);
			if( hist->name.Contains("diPt") )
			{
				hist->SetYRange(1e-1, 1e3);
				ratioMin = 0.4;
				ratioMax = 1.6;
			}
			if( hist->name.Contains("diRap") ) hist->SetYRange(1, 1e2);
		}

	}
};

void DrawHist_MassRange()
{
	HistDrawer* drawerMuMu = new HistDrawer( "MuMu" );
	drawerMuMu->DrawAll();

	HistDrawer* drawerEE = new HistDrawer( "EE" );
	drawerEE->DrawAll();
}