#include <Include/PlotTools.h>

class DrawingTool
{
public:
	TString channelType;
	DrawingTool( TString _channelType )
	{
		this->channelType = _channelType;
	}

	void Draw()
	{
		// this->Draw_Full_Zpeak();
		// this->Draw_Full_COLZ();
		this->Draw_ExpectZpeak_COLZ();
	}

	void Draw_x1Q2_x2Q2()
	{

	}

	void Draw_Full_Zpeak()
	{
		SampleInfo* sampleInfo_FullRange = new SampleInfo( 0, "FullRange", "15 < M < 3000 GeV");
		sampleInfo_FullRange->SetFileName( "./Local/ROOTFile_MakeHist_xQ2_"+this->channelType+".root");
		sampleInfo_FullRange->SetColor( kBlack );

		SampleInfo* sampleInfo_Zpeak = new SampleInfo( 0, "Zpeak", "60 < M < 120 GeV");
		sampleInfo_Zpeak->SetFileName( "./Local/ROOTFile_MakeHist_xQ2_"+this->channelType+".root");
		sampleInfo_Zpeak->SetColor( kRed );


		HistInfo* histInfo_x1Q2_GEN = new HistInfo( "h2D_x1Q2_GEN", "x1", "Q^{2} [GeV^{2}]" );
		histInfo_x1Q2_GEN->SetXRange( 1e-8, 1 );

		HistInfo* histInfo_x2Q2_GEN = new HistInfo( "h2D_x2Q2_GEN", "x2", "Q^{2} [GeV^{2}]" );
		histInfo_x2Q2_GEN->SetXRange( 1e-8, 1 );

		vector<HistInfo*> vec_HistInfo;
		vec_HistInfo.push_back( histInfo_x1Q2_GEN );
		vec_HistInfo.push_back( histInfo_x2Q2_GEN );

		for(const auto& histInfo : vec_HistInfo )
		{
			TString histName_FullRange = "DY"+this->channelType+"/"+histInfo->name;
			TH2Ext* h_FullRange = new TH2Ext( sampleInfo_FullRange, histInfo, histName_FullRange );

			TString histName_Zpeak = histName_FullRange;
			histName_Zpeak.ReplaceAll("GEN", "M60to120_GEN");
			TH2Ext* h_Zpeak = new TH2Ext( sampleInfo_Zpeak, histInfo, histName_Zpeak );

			TCanvas *c;
			TString CanvasName = histInfo->name;
			CanvasName.ReplaceAll("h2D_", "Local/c2D_"+this->channelType+"_");
			SetCanvas_Square( c, CanvasName, 1, 1 );

			h_FullRange->DrawAndSet( "SCAT" );
			h_Zpeak->DrawAndSet( "SCATSAME" );
			h_Zpeak->h->SetMarkerStyle( 24 );

			TLegend *legend;
			SetLegend( legend );
			h_FullRange->AddToLegend( legend );
			h_Zpeak->AddToLegend( legend );
			legend->Draw();

			c->SaveAs(".png");
			// c->SaveAs(".pdf");
		}
	}

	void Draw_Full_COLZ()
	{
		SampleInfo* sampleInfo_FullRange = new SampleInfo( 0, "FullRange", "15 < M^{post-FSR}_{ll} < 3000 GeV");
		sampleInfo_FullRange->SetFileName( "./Local/ROOTFile_MakeHist_xQ2_"+this->channelType+".root");
		sampleInfo_FullRange->SetColor( kBlack );

		vector<HistInfo*> vec_HistInfo;

		HistInfo* histInfo_x1Q2_GEN = new HistInfo( "h2D_x1Q2_GEN", "x_{1}", "Q^{2} [GeV^{2}]" );
		histInfo_x1Q2_GEN->SetXRange( 1e-5, 1 );
		histInfo_x1Q2_GEN->SetYRange( 1, 1e8 );
		histInfo_x1Q2_GEN->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x1Q2_GEN );

		HistInfo* histInfo_x2Q2_GEN = new HistInfo( "h2D_x2Q2_GEN", "x_{2}", "Q^{2} [GeV^{2}]" );
		histInfo_x2Q2_GEN->SetXRange( 1e-5, 1 );
		histInfo_x2Q2_GEN->SetYRange( 1, 1e8 );
		histInfo_x2Q2_GEN->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x2Q2_GEN );

		HistInfo* histInfo_x1Q2_M60to120_GEN = new HistInfo( "h2D_x1Q2_M60to120_GEN", "x_{1}", "Q^{2} [GeV^{2}]" );
		histInfo_x1Q2_M60to120_GEN->SetXRange( 1e-5, 1 );
		histInfo_x1Q2_M60to120_GEN->SetYRange( 1, 1e8 );
		histInfo_x1Q2_M60to120_GEN->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x1Q2_M60to120_GEN );

		HistInfo* histInfo_x2Q2_M60to120_GEN = new HistInfo( "h2D_x2Q2_M60to120_GEN", "x_{2}", "Q^{2} [GeV^{2}]" );
		histInfo_x2Q2_M60to120_GEN->SetXRange( 1e-5, 1 );
		histInfo_x2Q2_M60to120_GEN->SetYRange( 1, 1e8 );
		histInfo_x2Q2_M60to120_GEN->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x2Q2_M60to120_GEN );

		HistInfo* histInfo_x1Q2_RECO = new HistInfo( "h2D_x1Q2_RECO", "x_{1}", "Q^{2} [GeV^{2}]" );
		histInfo_x1Q2_RECO->SetXRange( 1e-5, 1 );
		histInfo_x1Q2_RECO->SetYRange( 1, 1e8 );
		histInfo_x1Q2_RECO->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x1Q2_RECO );

		HistInfo* histInfo_x2Q2_RECO = new HistInfo( "h2D_x2Q2_RECO", "x_{2}", "Q^{2} [GeV^{2}]" );
		histInfo_x2Q2_RECO->SetXRange( 1e-5, 1 );
		histInfo_x2Q2_RECO->SetYRange( 1, 1e8 );
		histInfo_x2Q2_RECO->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x2Q2_RECO );

		for(const auto& histInfo : vec_HistInfo )
		{
			TString histName_FullRange = "DY"+this->channelType+"/"+histInfo->name;
			TH2Ext* h_FullRange = new TH2Ext( sampleInfo_FullRange, histInfo, histName_FullRange );

			TCanvas *c;
			TString CanvasName = histInfo->name;
			CanvasName.ReplaceAll("h2D_", "Local/c2D_"+this->channelType+"_COLZ_");
			SetCanvas_Square2D( c, CanvasName, 1, 1 );
			c->SetLogz();

			h_FullRange->DrawAndSet( "COLZ" );
			h_FullRange->h->SetMinimum(1);
			gStyle->SetPalette(1);

			TLatex latex;
			Latex_Simulation( latex );

			TString channelInfo = "#mu channel";
			if( this->channelType == "EE" ) channelInfo = "e channel";
			latex.DrawLatexNDC(0.16, 0.91, "#font[62]{#scale[0.6]{"+channelInfo+"}}");

			TString massRangeInfo = sampleInfo_FullRange->fullName;
			if( histInfo->name.Contains("M60to120") ) massRangeInfo = "60 < M^{post-FSR}_{ll} < 120 GeV";
			if( histInfo->name.Contains("RECO") ) massRangeInfo = "15 < M^{RECO}_{ll} < 3000 GeV";
			latex.DrawLatexNDC(0.16, 0.87, "#font[42]{#scale[0.5]{"+massRangeInfo+"}}");


			// c->SaveAs(".png");
			c->SaveAs(".pdf");
		}
	}

	void Draw_ExpectZpeak_COLZ()
	{
		SampleInfo* sampleInfo = new SampleInfo( 0, "FullRange", "15 < M^{post-FSR}_{ll} < 3000 GeV except Z-peak");
		sampleInfo->SetFileName( "./Local/ROOTFile_MakeHist_xQ2_"+this->channelType+".root");
		sampleInfo->SetColor( kBlack );

		vector<HistInfo*> vec_HistInfo;

		HistInfo* histInfo_x1Q2_GEN = new HistInfo( "h2D_x1Q2_GEN", "x_{1}", "Q^{2} [GeV^{2}]" );
		histInfo_x1Q2_GEN->SetXRange( 1e-5, 1 );
		histInfo_x1Q2_GEN->SetYRange( 1, 1e8 );
		histInfo_x1Q2_GEN->SetZRange( 1, 2e6 );
		vec_HistInfo.push_back( histInfo_x1Q2_GEN );

		for(const auto& histInfo : vec_HistInfo )
		{
			TString histName_FullRange = "DY"+this->channelType+"/"+histInfo->name;
			TString histName_Zpeak = "DY"+this->channelType+"/"+histInfo->name;
			histName_Zpeak.ReplaceAll("_GEN", "_M60to120_GEN");

			TH2D* h_FullRange = Get_Hist_2D( sampleInfo->fileName, histName_FullRange );
			TH2D* h_Zpeak = Get_Hist_2D( sampleInfo->fileName, histName_Zpeak );
			TH2D* h2D_subtract = (TH2D*)h_FullRange->Clone(histInfo->name+"_SubtractZpeak");
			h2D_subtract->Add( h_Zpeak, -1 );

			TH2Ext* h2D = new TH2Ext( sampleInfo, histInfo, h2D_subtract );

			TCanvas *c;
			TString CanvasName = h2D->h->GetName();
			CanvasName.ReplaceAll("h2D_", "Local/c2D_"+this->channelType+"_COLZ_");
			SetCanvas_Square2D( c, CanvasName, 1, 1 );
			c->SetLogz();

			h2D->DrawAndSet( "COLZ" );
			h2D->h->SetMinimum(1);
			gStyle->SetPalette(1);

			TLatex latex;
			Latex_Simulation( latex );

			TString channelInfo = "#mu channel";
			if( this->channelType == "EE" ) channelInfo = "e channel";
			latex.DrawLatexNDC(0.16, 0.91, "#font[62]{#scale[0.6]{"+channelInfo+"}}");

			TString massRangeInfo = sampleInfo->fullName;
			latex.DrawLatexNDC(0.16, 0.87, "#font[42]{#scale[0.5]{"+massRangeInfo+"}}");


			// c->SaveAs(".png");
			c->SaveAs(".pdf");
		}


	}

private:
	void Init()
	{

	}
};

void DrawHist()
{
	DrawingTool *tool_MuMu = new DrawingTool( "MuMu" );
	tool_MuMu->Draw();

	DrawingTool *tool_EE = new DrawingTool( "EE" );
	tool_EE->Draw();
}