#include <Include/PlotTools.h>

class DrawingTool
{
public:
	TString fileName;
	vector<TString> vec_massRange;
	vector<Int_t> vec_color;
	DrawingTool()
	{
		this->Init();
	}

	void DrawAll()
	{
		for( const auto& massRange : this->vec_massRange )
		{
			this->Draw( massRange );

			vector<Double_t> vec_diRapBinEdge = this->DiRapBinEdge( massRange );
			Int_t nDiRapBin = (Int_t)vec_diRapBinEdge.size() - 1;

			for(Int_t i_diRap=0; i_diRap<nDiRapBin; i_diRap++)
				this->DrawRelUncPlot( massRange, i_diRap, vec_diRapBinEdge[i_diRap], vec_diRapBinEdge[i_diRap+1] );
		}
	}

private:
	void Draw( TString massRange )
	{
		// HistInfo* histInfo = new HistInfo("", "P_{T}(ll) [GeV]", "d^{2}#sigma/dP_{T}d|y| [pb/GeV]");
		HistInfo* histInfo = new HistInfo("", "P_{T}(ll) [GeV]", "d^{2}#sigma/dP_{T} [pb/GeV]");
		histInfo->SetFullName( massRange+"_diPt" );
		// histInfo->SetXRange( 0, 1000 );
		histInfo->SetYRange( 1e-6, 1e+3 );

		vector<Double_t> vec_diRapBinEdge = this->DiRapBinEdge( massRange );
		Int_t nDiRapBin = (Int_t)vec_diRapBinEdge.size() - 1;

		vector< TH1Ext* > vec_hext;
		for(Int_t i_diRap=0; i_diRap<nDiRapBin; i_diRap++)
		{
			TString histName = TString::Format("%s_DiRapBin%02d", massRange.Data(), i_diRap);
			cout << "histName: " << histName << endl;
			TH1D* h = Get_Hist( this->fileName, histName );
			// Double_t diRapBinSize = vec_diRapBinEdge[i_diRap+1] - vec_diRapBinEdge[i_diRap];
			// h->Scale( 1 / diRapBinSize );
			h = DivideEachBin_ByBinWidth( h );

			SampleInfo *si = new SampleInfo( 0, TString::Format("DiRapBin%02d", i_diRap), TString::Format("%.1lf < |y| < %.1lf", vec_diRapBinEdge[i_diRap], vec_diRapBinEdge[i_diRap+1]) );
			si->SetColor( this->vec_color[i_diRap] );

			TH1Ext *hext = new TH1Ext( si, histInfo, h );
			vec_hext.push_back( hext );
		}

		TString CanvasName = histInfo->fullName;
		TCanvas *c;
		SetCanvas_Square( c, CanvasName, 0, 1 );
		c->cd();

		TLegend *legend;
		SetLegend( legend, 0.45, 0.80, 0.95, 0.95 );
		legend->SetNColumns(2);
		for(const auto& h : vec_hext )
		{
			h->DrawAndSet("EPSAME");
			h->AddToLegend( legend );
			// h->h->GetXaxis()->SetRangeUser(5, 1000);
		}
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );
		latex.DrawLatexNDC(0.16, 0.91, "#font[42]{#scale[0.7]{FEWZ (NNPDF 3.1)}}");
		TString massRangeInfo = "";
		if( massRange == "M15to64" ) massRangeInfo = "15 < M < 64 GeV";
		if( massRange == "M64to106" ) massRangeInfo = "64 < M < 106 GeV";
		if( massRange == "M106to3000" ) massRangeInfo = "106 < M < 3000 GeV";
		latex.DrawLatexNDC(0.16, 0.87, "#font[42]{#scale[0.6]{"+massRangeInfo+"}}");

		c->SaveAs(".pdf");
	}

	void DrawRelUncPlot( TString massRange, Int_t i_diRap, Double_t lowerDiRap, Double_t upperDiRap )
	{
		SampleInfo *si_integ = new SampleInfo( 0, "relUnc_integ", "Integration error" );
		si_integ->SetColor( kBlue );
		SampleInfo *si_pdf = new SampleInfo( 0, "relUnc_pdf", "PDF error" );
		si_pdf->SetColor( kGreen+2 );
		SampleInfo *si_tot = new SampleInfo( 0, "relUnc_tot", "Total error" );
		si_tot->SetColor( kRed );

		TString histName_base = TString::Format("%s_DiRapBin%02d", massRange.Data(), i_diRap);

		HistInfo* histInfo = new HistInfo("", "P_{T}(ll) [GeV]", "Rel. unc. (%)");
		histInfo->SetFullName( histName_base+"_relUnc" );
		histInfo->SetYRange( 1e-2, 1e2 );

		TH1D* h_relUnc_integ = Get_Hist(  this->fileName, histName_base+"_relUnc_integ" );
		h_relUnc_integ->Scale( 100 );
		TH1D* h_relUnc_pdf = Get_Hist(  this->fileName, histName_base+"_relUnc_pdf" );
		h_relUnc_pdf->Scale( 100 );
		TH1D* h_relUnc_tot = Get_Hist(  this->fileName, histName_base+"_relUnc_tot" );
		h_relUnc_tot->Scale( 100 );

		TH1Ext* hext_relUnc_integ = new TH1Ext( si_integ, histInfo, h_relUnc_integ );
		TH1Ext* hext_relUnc_pdf = new TH1Ext( si_pdf, histInfo, h_relUnc_pdf );
		TH1Ext* hext_relUnc_tot = new TH1Ext( si_tot, histInfo, h_relUnc_tot );

		vector< TH1Ext* > vec_hext = {hext_relUnc_integ, hext_relUnc_pdf, hext_relUnc_tot};

		TString CanvasName = histInfo->fullName;
		TCanvas *c;
		SetCanvas_Square( c, CanvasName, 0, 1 );
		c->cd();

		TLegend *legend;
		SetLegend( legend, 0.60, 0.80, 0.95, 0.95 );
		// legend->SetNColumns(2);
		for(const auto& h : vec_hext )
		{
			h->DrawAndSet("LPSAME");
			h->AddToLegend( legend );
			// h->h->GetXaxis()->SetRangeUser(5, 1000);
		}
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );
		latex.DrawLatexNDC(0.16, 0.91, "#font[42]{#scale[0.7]{FEWZ (NNPDF 3.1)}}");
		TString massRangeInfo = "";
		if( massRange == "M15to64" ) massRangeInfo = "15 < M < 64 GeV";
		if( massRange == "M64to106" ) massRangeInfo = "64 < M < 106 GeV";
		if( massRange == "M106to3000" ) massRangeInfo = "106 < M < 3000 GeV";
		TString diRapInfo = TString::Format("%.1lf < |y| < %.1lf", lowerDiRap, upperDiRap);
		TString fullInfo = massRangeInfo + ", " + diRapInfo;
		latex.DrawLatexNDC(0.16, 0.87, "#font[42]{#scale[0.6]{"+fullInfo+"}}");

		c->SaveAs(".pdf");
	}

	void Init()
	{
		this->fileName = "ROOTFile_FEWZ_diPt.root";

		this->vec_massRange.push_back( "M15to64" );
		this->vec_massRange.push_back( "M64to106" );
		this->vec_massRange.push_back( "M106to3000" );

		this->vec_color.push_back( kBlack );
		this->vec_color.push_back( kRed );
		this->vec_color.push_back( kBlue );
		this->vec_color.push_back( kGreen+2 );
		this->vec_color.push_back( kViolet );
	}

	vector<Double_t> DiRapBinEdge( TString massRange )
	{
		vector<Double_t> vec_BinEdge;
		if( massRange == "M15to64" )
		{
			vec_BinEdge.push_back( 0 );
			vec_BinEdge.push_back( 0.7 );
			vec_BinEdge.push_back( 1.1 );
			vec_BinEdge.push_back( 1.9 );
			vec_BinEdge.push_back( 2.4 );
			vec_BinEdge.push_back( 1000 );
		}
		if( massRange == "M64to106" )
		{
			vec_BinEdge.push_back( 0 );
			vec_BinEdge.push_back( 0.7 );
			vec_BinEdge.push_back( 1.9 );
			vec_BinEdge.push_back( 1000 );
		}
		if( massRange == "M106to3000" )
		{
			vec_BinEdge.push_back( 0 );
			vec_BinEdge.push_back( 0.7 );
			vec_BinEdge.push_back( 1000 );
		}

		return vec_BinEdge;
	}
};

void DrawHist()
{
	DrawingTool *tool = new DrawingTool();
	tool->DrawAll();
}