#include <Include/PlotTools.h>

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
		for( const auto& hist : this->vec_Hist )
		{
			this->Draw( hist );
		}
	}

private:
	void Init()
	{
		TString fileNameInc = TString::Format("Local/ROOTFile_MakeHist_GenLevel_%s_Inclusive.root", this->channelType.Data());
		// TString fileNameInc = TString::Format("Local/ROOTFile_MakeHist_GenLevel_%s_MassBinned.root", this->channelType.Data());
		this->siInc = new SampleInfo( 0, "Inclusive", "Inclusive (M50-inf)" );
		this->siInc->SetFileName( fileNameInc );
		this->siInc->SetColor( kRed );

		TString fileNameBinned = TString::Format("Local/ROOTFile_MakeHist_GenLevel_%s_MassBinned.root", this->channelType.Data());
		this->siBinned = new SampleInfo( 0, "Mass-binned", "Mass-binned" );
		this->siBinned->SetFileName( fileNameBinned );
		this->siBinned->SetColor( kBlack );

		// -- histogram information to be drawn -- //
		HistInfo *histMass10GeV = new HistInfo( "DY"+this->channelType+"/h_mass_preFSR", "m [GeV]", "Entries / 10 GeV");
		histMass10GeV->SetFullName( "DY"+this->channelType+"_Mass10GeV");
		histMass10GeV->SetRebinX( 10 );
		histMass10GeV->SetXRange( 50, 3000 );
		histMass10GeV->SetYRange( 1e-3, 1e9 );
		this->vec_Hist.push_back( histMass10GeV );

		HistInfo *histMass50GeV = new HistInfo( "DY"+this->channelType+"/h_mass_preFSR", "m [GeV]", "Entries / 50 GeV");
		histMass50GeV->SetFullName( "DY"+this->channelType+"_Mass50GeV");
		histMass50GeV->SetRebinX( 50 );
		histMass50GeV->SetXRange( 50, 3000 );
		histMass50GeV->SetYRange( 1e-3, 1e9 );
		this->vec_Hist.push_back( histMass50GeV );

		HistInfo *histMass10GeVZoomIn = new HistInfo( "DY"+this->channelType+"/h_mass_preFSR", "m [GeV]", "Entries / 10 GeV");
		histMass10GeVZoomIn->SetFullName( "DY"+this->channelType+"_Mass10GeVZoomIn");
		histMass10GeVZoomIn->SetRebinX( 10 );
		histMass10GeVZoomIn->SetXRange( 50, 500 );
		histMass10GeVZoomIn->SetYRange( 1e2, 1e9 );
		this->vec_Hist.push_back( histMass10GeVZoomIn );

		HistInfo *histMass50GeVZoomIn = new HistInfo( "DY"+this->channelType+"/h_mass_preFSR", "m [GeV]", "Entries / 50 GeV");
		histMass50GeVZoomIn->SetFullName( "DY"+this->channelType+"_Mass50GeVZoomIn");
		histMass50GeVZoomIn->SetRebinX( 50 );
		histMass50GeVZoomIn->SetXRange( 50, 500 );
		histMass50GeVZoomIn->SetYRange( 1e2, 1e9 );
		this->vec_Hist.push_back( histMass50GeVZoomIn );

		HistInfo *histMass50GeVZoomIn2 = new HistInfo( "DY"+this->channelType+"/h_mass_preFSR", "m [GeV]", "Entries / 50 GeV");
		histMass50GeVZoomIn2->SetFullName( "DY"+this->channelType+"_Mass50GeVZoomIn2");
		histMass50GeVZoomIn2->SetRebinX( 50 );
		histMass50GeVZoomIn2->SetXRange( 50, 1000 );
		histMass50GeVZoomIn2->SetYRange( 1, 1e9 );
		this->vec_Hist.push_back( histMass50GeVZoomIn2 );
	}

	void Draw( HistInfo* hist )
	{
		TH1Ext *hextInc = new TH1Ext( this->siInc, hist );
		hextInc->CalcRatio_DEN( hextInc->h ); // -- dummy -- //
		TH1Ext *hextBinned = new TH1Ext( this->siBinned, hist );
		hextBinned->CalcRatio_DEN( hextInc->h );


		TString CanvasName = "Local/c_"+hist->fullName;
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

		c->cd();
		BottomPad->cd();
		hextBinned->DrawRatioAndSet( "EPSAME", "Binned/Inc.", 0.9, 1.1 );
		hextBinned->h_ratio->GetYaxis()->SetNdivisions(006);

		TF1 *f_line;
		DrawLine( f_line );

		c->SaveAs(".pdf");
	}
};

void DrawHist()
{
	HistDrawer* drawerMuMu = new HistDrawer( "MuMu" );
	drawerMuMu->DrawAll();

	HistDrawer* drawerEE = new HistDrawer( "EE" );
	drawerEE->DrawAll();
}