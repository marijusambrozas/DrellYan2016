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
		this->Draw_Full_Zpeak();
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
			SetCanvas_Square( c, CanvasName );

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

private:
	void Init()
	{

	}
};

void DrawHist()
{
	DrawingTool *tool = new DrawingTool( "MuMu" );
	tool->Draw();

}