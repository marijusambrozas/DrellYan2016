#include <Include/PlotTools.h>

class PrintingTool
{
public:
	vector< TH1D* > vec_Hist;
	TString channelType;
	TString fileName;

	PrintingTool(TString _channelType, TString _fileName = "ROOTFile_GetResolution.root")
	{
		this->channelType = _channelType;
		this->fileName = _fileName;
		this->Init();
	}

	void PrintAll()
	{
		for( const auto& h : this->vec_Hist )
		{
			this->PrintHistogram( h );
		}
	}

private:
	void Init()
	{
		TH1D* h_Mean = Get_Hist(this->fileName, "h_Mean_"+this->channelType);
		this->vec_Hist.push_back( h_Mean );

		TH1D* h_Resolution = Get_Hist(this->fileName, "h_Resolution_"+this->channelType);
		this->vec_Hist.push_back( h_Resolution );

		TH1D* h_RMS = Get_Hist(this->fileName, "h_RMS_"+this->channelType);
		this->vec_Hist.push_back( h_RMS );

		TH1D* h_NormChi2 = Get_Hist(this->fileName, "h_NormChi2_"+this->channelType);
		this->vec_Hist.push_back( h_NormChi2 );
	}

	void PrintHistogram( TH1D* h )
	{
		TCanvas *c;
		TString CanvasName = h->GetName();
		CanvasName.ReplaceAll("h_", "Local/c_"+this->channelType+"_");

		HistInfo *Hist = new HistInfo( kBlack, "temp", h );

		SetCanvas_Square( c, CanvasName, 1, 0 );
		c->cd();
		Hist->Draw("HISTLPSAME");

		TString HistName = h->GetName();
		TString titleY = this->GetTitle_Y( HistName );
		SetHistFormat_SinglePad(Hist->h, "m [GeV]", titleY);
		if( HistName.Contains("Resolution") )
			Hist->h->GetYaxis()->SetRangeUser(0, 0.05);

		TLatex latex;
		Latex_Simulation( latex );

		c->SaveAs(".pdf");
	}

	TString GetTitle_Y( TString HistName )
	{
		TString titleY = "";
		if( HistName.Contains("NormChi2") )
			titleY = "Normalized #chi^{2}";
		else if( HistName.Contains("RMS") )
			titleY = "RMS";
		else if( HistName.Contains("Mean") )
			titleY = "Mean";
		else if( HistName.Contains("Resolution") )
			titleY = "Resolution";

		return titleY;
	}

};
void PrintHistogram()
{
	PrintingTool* tool = new PrintingTool("EE");
	tool->PrintAll();

}