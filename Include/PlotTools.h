#pragma once

#include <TH1D.h>
#include <TH2D.h>
#include <TColor.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TFile.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>

#include <vector>


TH1D* Get_Hist(TString FileName, TString HistName, TString HistName_New = "" )
{
	TH1::AddDirectory(kFALSE);

	TFile *f_input = TFile::Open( FileName );
	TH1D* h_temp = (TH1D*)f_input->Get(HistName)->Clone();
	if( HistName_New != "" )
		h_temp->SetName( HistName_New );

	f_input->Close();
	// delete f_input;

	return h_temp;
}

TH2D* Get_Hist_2D(TString FileName, TString HistName, TString HistName_New = "" )
{
	TH1::AddDirectory(kFALSE);

	TFile *f_input = TFile::Open( FileName );
	TH2D* h_temp = (TH2D*)f_input->Get(HistName)->Clone();
	if( HistName_New != "" )
		h_temp->SetName( HistName_New );

	f_input->Close();
	// delete f_input;

	return h_temp;
}

TGraphAsymmErrors* Get_Graph(TString FileName, TString GraphName, TString GraphName_New = "" )
{
	TFile *f_input = TFile::Open( FileName );
	TGraphAsymmErrors* g_temp = (TGraphAsymmErrors*)f_input->Get(GraphName)->Clone();
	if( GraphName_New != "" )
		g_temp->SetName( GraphName_New );

	f_input->Close();

	return g_temp;
}

void SetAxis_SinglePad( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
{
	X_axis->SetTitle( XTitle );
	X_axis->SetLabelSize(0.04);
	X_axis->SetTitleOffset(1.1);
	X_axis->SetTitleSize(0.05);
	X_axis->SetNoExponent();
	X_axis->SetMoreLogLabels();

	Y_axis->SetTitle( YTitle );
	Y_axis->SetTitleSize(0.05);
	Y_axis->SetTitleOffset(1.2);
	Y_axis->SetLabelSize(0.045);
}

void SetAxis_TopPad( TAxis *X_axis, TAxis *Y_axis, TString YTitle )
{
	X_axis->SetLabelFont(42);
	X_axis->SetLabelSize(0.000);
	X_axis->SetTitleSize(0.000);

	Y_axis->SetTitleFont(42);
	Y_axis->SetTitle( YTitle );
	Y_axis->SetTitleSize(0.05);
	Y_axis->SetTitleFont(42);
	Y_axis->SetTitleOffset(1.25);
	Y_axis->SetLabelFont(42);
	Y_axis->SetLabelSize(0.04);
}

void SetAxis_BottomPad( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle, Double_t yMin = 0.51, Double_t yMax = 1.49 )
{
	X_axis->SetMoreLogLabels();
	X_axis->SetNoExponent();
	X_axis->SetTitle( XTitle );
	X_axis->SetTitleOffset( 0.85 );
	X_axis->SetTitleSize( 0.2 );
	X_axis->SetLabelColor(1);
	X_axis->SetLabelFont(42);
	X_axis->SetLabelOffset(0.01);
	X_axis->SetLabelSize(0.13);

	Y_axis->SetTitle( YTitle );
	Y_axis->SetTitleOffset( 0.55 );
	Y_axis->SetTitleSize( 0.12);
	Y_axis->SetRangeUser( yMin, yMax );
	Y_axis->SetLabelSize( 0.10 );
}

void SetAxis_2D( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
{
	X_axis->SetTitle( XTitle );
	X_axis->SetLabelSize(0.03);
	X_axis->SetTitleOffset(1.35);
	X_axis->SetTitleSize(0.04);
	// X_axis->SetNoExponent();
	// X_axis->SetMoreLogLabels();

	Y_axis->SetTitle( YTitle );
	Y_axis->SetTitleSize(0.04);
	Y_axis->SetTitleOffset(1.5);
	Y_axis->SetLabelSize(0.03);
}

void SetCanvas_Square( TCanvas*& c, TString CanvasName, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t SizeX = 800, Double_t SizeY = 800 )
{
	c = new TCanvas(CanvasName, "", SizeX, SizeY);
	c->cd();
	
	c->SetTopMargin(0.05);
	c->SetLeftMargin(0.13);
	c->SetRightMargin(0.045);
	c->SetBottomMargin(0.13);

	if( isLogx )
		c->SetLogx();
	if( isLogy )
		c->SetLogy();
}

void SetCanvas_Square2D( TCanvas*& c, TString CanvasName, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t SizeX = 800, Double_t SizeY = 800 )
{
	SetCanvas_Square( c, CanvasName, isLogx, isLogy, SizeX, SizeY );
	c->SetRightMargin(0.12);
}

void SetCanvas_Ratio( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
{
	c = new TCanvas(CanvasName, "", 800, 800);
	c->cd();

	TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
	TopPad->Draw();
	TopPad->cd();

	TopPad->SetTopMargin(0.05);
	TopPad->SetLeftMargin(0.13);
	TopPad->SetRightMargin(0.045);
	TopPad->SetBottomMargin(0.3);

	if( isLogx )
		TopPad->SetLogx();
	if( isLogy )
		TopPad->SetLogy();

	c->cd();
	BottomPad = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.29 );
	BottomPad->Draw();
	BottomPad->cd();
	BottomPad->SetGridx();
	BottomPad->SetGridy();
	BottomPad->SetTopMargin(0.05);
	BottomPad->SetBottomMargin(0.4);
	BottomPad->SetRightMargin(0.045);
	BottomPad->SetLeftMargin(0.13);

	if( isLogx )
		BottomPad->SetLogx();
}

void SetLegend( TLegend *& legend, Double_t xMin = 0.75, Double_t yMin = 0.75, Double_t xMax = 0.95, Double_t yMax = 0.95 )
{
	legend = new TLegend( xMin, yMin, xMax, yMax );
	legend->SetFillStyle(0);
	legend->SetBorderSize(0);
	legend->SetTextFont( 62 );
}

TH1D* QuadSum_NoError( TH1D* h1, TH1D* h2 )
{
	TH1D* h_QuadSum = (TH1D*)h1->Clone( "h_QuadSum" );
	Int_t nBin = h1->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;

		Double_t value1 = h1->GetBinContent(i_bin);
		Double_t value2 = h2->GetBinContent(i_bin);

		Double_t QuadSum = sqrt( value1*value1 + value2*value2 );

		h_QuadSum->SetBinContent( i_bin, QuadSum );
		h_QuadSum->SetBinError( i_bin, 0 );
	}
	return h_QuadSum;
}

void AssignErrors( TH1D* h_cv, TH1D* h_RelUnc, Bool_t isPercent = kFALSE )
{
	if( h_cv->GetNbinsX() != h_RelUnc->GetNbinsX() )
	{
		printf("# bins for central value and relative uncertainty histograms are not same! .. %d != %d\n", 
			h_cv->GetNbinsX(), h_RelUnc->GetNbinsX() );

		return;
	}

	Int_t nBin = h_cv->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;

		Double_t value = h_cv->GetBinContent(i_bin);
		Double_t RelUnc = h_RelUnc->GetBinContent(i_bin);
		if( isPercent ) RelUnc = RelUnc / 100.0;

		Double_t AbsUnc = value * RelUnc;

		h_cv->SetBinError(i_bin, AbsUnc);
	}
}

void AssignErrors_Graph( TGraphAsymmErrors* g, TH1D* h_RelUnc, Bool_t isPercent = kFALSE )
{
	Int_t nPoint = g->GetN();
	Int_t nBin = h_RelUnc->GetNbinsX();
	if( nPoint != nBin )
	{
		printf("[nPoint != nBin: %d, %d]\n", nPoint, nBin);
		return;
	}

	for(Int_t i=0; i<nPoint; i++)
	{
		Double_t x = 0;
		Double_t y = 0;
		g->GetPoint( i, x, y );

		Int_t i_bin = i+1;
		Double_t RelUnc = h_RelUnc->GetBinContent(i_bin);
		if( isPercent ) RelUnc = RelUnc / 100.0;
		
		Double_t AbsUnc = y * RelUnc;

		g->SetPointEYlow(i, AbsUnc );
		g->SetPointEYhigh(i, AbsUnc );
	}
}

TH1D* Convert_GraphToHist( TGraphAsymmErrors *g )
{
	const Int_t nBin = g->GetN();
	Double_t *BinEdges = new Double_t[nBin+1];
	Double_t *value = new Double_t[nBin];
	Double_t *error = new Double_t[nBin];

	for(Int_t i=0; i<nBin; i++)
	{
		Double_t x, y;
		g->GetPoint(i, x, y);

		// -- make BinEdges array -- //
		Double_t ErrX_Low = g->GetErrorXlow(i);
		Double_t ErrX_High = g->GetErrorXhigh(i);

		if( i == nBin-1 )
		{
			BinEdges[i] = x - ErrX_Low;
			BinEdges[i+1] = x + ErrX_High;
		}
		else
			BinEdges[i] = x - ErrX_Low;


		// -- store graph information -- //
		value[i] = y;

		Double_t ErrY_Low = g->GetErrorYlow(i);
		Double_t ErrY_High = g->GetErrorYhigh(i);

		// -- take the larger one -- //
		error[i] = ErrY_Low > ErrY_High ? ErrY_Low : ErrY_High;
	}

	TString GraphName = g->GetName();
	TH1D* h_temp = new TH1D( "h_"+GraphName, "", nBin, BinEdges );

	// -- fill this histogram using graph information -- //
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;
		h_temp->SetBinContent( i_bin, value[i] );
		h_temp->SetBinError( i_bin, error[i] );
	}

	return h_temp;
}

TH1D* Extract_RelUnc( TH1D* h, TString HistName = "", Bool_t ConvertToPercent = kFALSE )
{
	TH1D* h_RelUnc = (TH1D*)h->Clone();
	if( HistName != "" )
		h_RelUnc->SetName(HistName);

	Int_t nBin = h->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;

		Double_t value = h->GetBinContent(i_bin);
		Double_t error = h->GetBinError(i_bin);

		Double_t RelUnc = error / value;
		if( ConvertToPercent )
			RelUnc = RelUnc * 100;

		h_RelUnc->SetBinContent(i_bin, RelUnc );
		h_RelUnc->SetBinError(i_bin, 0);
	}

	return h_RelUnc;
}

TH1D* ConvertHist_AbsToRel( TH1D* h_CenV, TH1D* h_AbsUnc, Bool_t ConvertToPercent = kFALSE )
{
	TH1D* h_RelUnc = (TH1D*)h_AbsUnc->Clone();

	Int_t nBin = h_CenV->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;

		Double_t CenV = h_CenV->GetBinContent(i_bin);
		Double_t AbsUnc = h_AbsUnc->GetBinContent(i_bin);
		Double_t RelUnc = AbsUnc / CenV;
		if( ConvertToPercent )
			RelUnc = RelUnc * 100;

		h_RelUnc->SetBinContent(i_bin, RelUnc );
		h_RelUnc->SetBinError(i_bin, 0);
	}

	return h_RelUnc;
}

TH1D* DivideEachBin_ByBinWidth( TH1D* h, TString HistName = "" )
{
	TH1D* h_return = (TH1D*)h->Clone();
	if( HistName != "" )
		h_return->SetName(HistName);

	Int_t nBin = h->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;
		Double_t Entry_before = h->GetBinContent(i_bin);
		Double_t Error_before = h->GetBinError(i_bin);
		Double_t BinWidth = h->GetBinWidth(i_bin);

		Double_t Entry_after = Entry_before / BinWidth;
		Double_t Error_after = Error_before / BinWidth;

		h_return->SetBinContent(i_bin, Entry_after);
		h_return->SetBinError(i_bin, Error_after);
	}

	return h_return;
}

TH1D* MultiplyEachBin_byBinWidth( TH1D* h, TString HistName = "" )
{
	TH1D* h_return = (TH1D*)h->Clone();
	if( HistName != "" )
		h_return->SetName(HistName);

	Int_t nBin = h->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;
		Double_t Entry_before = h->GetBinContent(i_bin);
		Double_t Error_before = h->GetBinError(i_bin);
		Double_t BinWidth = h->GetBinWidth(i_bin);

		Double_t Entry_after = Entry_before * BinWidth;
		Double_t Error_after = Error_before * BinWidth;

		h_return->SetBinContent(i_bin, Entry_after);
		h_return->SetBinError(i_bin, Error_after);
	}

	return h_return;
}

void SaveAsHist_OneContent( Double_t content, TString HistName, TFile *f_output )
{
	TH1D *h = new TH1D(HistName, "", 1, 0, 1);
	h->SetBinContent(1, content);
	h->SetBinError(1, 0);
	f_output->cd();
	h->Write();
}

void SaveAsHist_OneContent_WithError( Double_t content, Double_t error, TString HistName, TFile *f_output )
{
	TH1D *h = new TH1D(HistName, "", 1, 0, 1);
	h->SetBinContent(1, content);
	h->SetBinError(1, error);
	f_output->cd();
	h->Write();
}

Double_t GetContent_OneBinHist( TString FileName, TString HistName )
{
	TH1D* h_temp = Get_Hist( FileName, HistName );
	if( h_temp->GetNbinsX() != 1 )
	{
		cout << "This histogram has more than 1 bin! ... please check. Return -999" << endl;
		return -999;
	}

	return h_temp->GetBinContent(1);
}

void Print_Histogram( TH1D* h, Bool_t NegativeCheck = kFALSE )
{
	h->Print();

	// -- underflow -- //
	Double_t value_uf = h->GetBinContent(0);
	Double_t errorAbs_uf = h->GetBinError(0);
	Double_t errorRel_uf = value_uf == 0 ? 0 : errorAbs_uf / value_uf;

	printf( "Underflow: (value, error) = (%lf, %lf (%7.3lf %%))\n", 
		     value_uf, errorAbs_uf, errorRel_uf*100 );

	if( NegativeCheck && value_uf < 0 )
		printf("################## NEGATIVE BIN ##################");

	Int_t nBin = h->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;
		Double_t LowerEdge = h->GetBinLowEdge(i_bin);
		Double_t UpperEdge = h->GetBinLowEdge(i_bin+1);

		Double_t value = h->GetBinContent(i_bin);
		Double_t errorAbs = h->GetBinError(i_bin);
		Double_t errorRel;
		if( value != 0 )
			errorRel = errorAbs / value;
		else
			errorRel = 0;

		printf( "%02d bin: [%6.1lf, %6.1lf] (value, error) = (%lf, %lf (%7.3lf %%))\n", 
			     i_bin, LowerEdge, UpperEdge, value, errorAbs, errorRel*100 );
		
		if( NegativeCheck && value < 0 )
			printf("################## NEGATIVE BIN ##################");
	}

	// -- overflow -- //
	Double_t value_of = h->GetBinContent(nBin+1);
	Double_t errorAbs_of = h->GetBinError(nBin+1);
	Double_t errorRel_of = value_of == 0 ? 0 : errorAbs_of / value_of;

	printf( "Overflow: (value, error) = (%lf, %lf (%7.3lf %%))\n", 
		     value_of, errorAbs_of, errorRel_of*100 );

	if( NegativeCheck && value_of < 0 )
		printf("################## NEGATIVE BIN ##################");

	printf("\n\n");
}

void Print_Graph( TGraphAsymmErrors* g )
{
	TString GraphName = g->GetName();
	printf("[GraphName: %s]\n", GraphName.Data());
	Int_t nPoint = g->GetN();
	for(Int_t i=0; i<nPoint; i++)
	{
		Double_t x, y;
		g->GetPoint(i, x, y);

		Double_t xErrLow = g->GetErrorXlow(i);
		Double_t xErrHigh = g->GetErrorXhigh(i);
		Double_t LowerEdge = x - xErrLow;
		Double_t UpperEdge = x + xErrHigh;

		Double_t yErrLow = g->GetErrorYlow(i);
		Double_t yRelErrLow = yErrLow / y;
		Double_t yErrHigh = g->GetErrorYhigh(i);
		Double_t yRelErrHigh = yErrHigh / y;

		printf( "%02d point: [%6.1lf, %6.1lf] (value, errorLow, errorHigh) = (%lf, %lf (%7.3lf %%), %lf (%7.3lf %%))\n", 
			     i, LowerEdge, UpperEdge, y, yErrLow, yRelErrLow*100, yErrHigh, yRelErrHigh*100 );
	}
	printf("\n\n");
}

// -- http://igotit.tistory.com/entry/C-함수-인자로-포인터-전달하고-함수내에서-동적-메모리-할당-받기-2가지-방식 -- //
void DrawLine( TF1*& f_line, Int_t color = kRed )
{
	f_line = new TF1("f_line", "1", -10000, 10000);
	f_line->SetLineColor(color);
	f_line->SetLineWidth(1);
	f_line->Draw("PSAME");
}

void Latex_Preliminary_NoDataInfo( TLatex &latex )
{
	latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Preliminary}}}");
}

void Latex_Preliminary( TLatex &latex, Double_t lumi  )
{
	Latex_Preliminary_NoDataInfo( latex );
	latex.DrawLatexNDC(0.70, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%.2lf", lumi)+" fb^{-1} (13 TeV)}}");
}

void Latex_Preliminary( TLatex &latex, Double_t lumi, Int_t E_CM  )
{
	Latex_Preliminary_NoDataInfo( latex );
	latex.DrawLatexNDC(0.69, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%.1lf fb^{-1} (%d TeV)", lumi, E_CM)+"}}");
}

void Latex_Simulation( TLatex &latex )
{
	latex.DrawLatexNDC(0.82, 0.96, "#font[42]{#scale[0.8]{13 TeV}}");
	latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Simulation}}}");
}

void Latex_Preliminary_EffPlotStyle( TLatex &latex, Int_t year, Int_t E_CM = 13 )
{
	latex.DrawLatexNDC(0.70, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%d, %d TeV", year, E_CM)+"}}");
	latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}");
	latex.DrawLatexNDC(0.24, 0.96, "#font[42]{#it{#scale[0.8]{Preliminary}}}");
}


class SampleInfo
{
public:
	TString name; // -- short name for the convenience in the code: ex> DYMuMu_M50 -- //
	TString fullName; // -- for legend: ex> DY #rightarrow Z/#gamma* -- //

	Bool_t isRealData;
	TString fileName;
	Int_t color;

	Double_t xSec;
	Double_t sumWeight; // -- same with # events for the samples without negative weights -- //
	Double_t normFactor;

	SampleInfo()
	{
		this->Init();
	}

	SampleInfo( Bool_t _isRealData, TString _name, TString _fullName )
	{
		this->isRealData = _isRealData;
		this->SetName( _name, _fullName );
	}

	void SetName( TString _name, TString _fullName )
	{
		this->name = _name;
		this->fullName = _fullName;
	}

	void SetFileName( TString _name )
	{
		this->fileName = _name;
	}
	
	void SetColor( Int_t _color )
	{
		this->color = _color;
	}

	void SetNormFactor( Double_t lumi, Double_t _xSec, Double_t _sumWeight )
	{
		if( this->isRealData )
		{
			cout << "[SetNormFactor] This is real data! ... something goes wrong" << endl;
			return;
		}

		this->xSec = _xSec;
		this->sumWeight = _sumWeight;
		this->normFactor = (lumi * this->xSec) / this->sumWeight;
		printf("[SetNormFactor] Sample: %s\n", name.Data() );
		printf("\tNormalization factor = (%.3lf * %.3e) / (%.1lf) = %.3e\n", lumi, this->xSec, this->sumWeight, this->normFactor);
	}

private:
	void Init()
	{
		this->name = "";
		this->fullName = "";

		this->isRealData = kFALSE;
		this->fileName = "";
		this->color = 0;

		this->xSec = -999;
		this->sumWeight = -999;
		this->normFactor = -999;
	}
};


class HistInfo
{
public:
	TString name;
	TString titleX;
	TString titleY;

	Bool_t hasXRange;
	Double_t minX;
	Double_t maxX;

	Bool_t hasYRange;
	Double_t minY;
	Double_t maxY;

	Bool_t hasZRange;
	Double_t minZ;
	Double_t maxZ;

	Double_t hasRebinX;
	Double_t nRebinX;

	Double_t hasRebinY;
	Double_t nRebinY;

	Bool_t isFilled;

	// -- log axis: property of canvases, not histograms -- //
	// Bool_t isLogX;
	// Bool_t isLogY;

	HistInfo()
	{
		this->Init();
	}

	HistInfo( TString _name ): HistInfo()
	{
		this->name = _name;
	}

	HistInfo( TString _name, TString _titleX, TString _titleY ): HistInfo()
	{
		this->name = _name;
		this->SetTitle( _titleX, _titleY );
	}

	HistInfo( TString _titleX, TString _titleY ): HistInfo()
	{
		this->SetTitle( _titleX, _titleY );
	}

	void SetTitle( TString X, TString Y )
	{
		this->titleX = X;
		this->titleY = Y;
	}

	// void SetLogAxis( Bool_t X, Bool_t Y )
	// {
	// 	this->isLogX = X;
	// 	this->isLogY = Y;
	// }

	void SetXRange( Double_t min, Double_t max )
	{
		this->hasXRange = kTRUE;
		this->minX = min;
		this->maxX = max;
	}

	void SetYRange( Double_t min, Double_t max )
	{
		this->hasYRange = kTRUE;
		this->minY = min;
		this->maxY = max;
	}

	void SetZRange( Double_t min, Double_t max )
	{
		this->hasZRange = kTRUE;
		this->minZ = min;
		this->maxZ = max;
	}

	void SetRebinX( Int_t _nRebin )
	{
		this->hasRebinX = kTRUE;
		this->nRebinX = _nRebin;
	}

	void SetRebinY( Int_t _nRebin )
	{
		this->hasRebinY = kTRUE;
		this->nRebinY = _nRebin;
	}

	void IsFilled( Bool_t _isFilled = kTRUE )
	{
		this->isFilled = _isFilled;
	}

private:
	void Init()
	{
		this->name = "";
		this->titleX = "";
		this->titleY = "";

		this->hasXRange = kFALSE;
		this->minX = -999;
		this->maxX = -999;

		this->hasYRange = kFALSE;
		this->minY = -999;
		this->maxY = -999;

		this->hasZRange = kFALSE;
		this->minZ = -999;
		this->maxZ = -999;

		this->hasRebinX = kFALSE;
		this->nRebinX = 1;

		this->hasRebinY = kFALSE;
		this->nRebinY = 1;

		this->isFilled = kFALSE;

		// this->isLogX = kFALSE;
		// this->isLogY = kFALSE;
	}

};

// -- TH1 extension -- //
class TH1Ext
{
public:
	TH1D* h;

	Bool_t hasRatio;
	TH1D* h_ratio;

	SampleInfo* sampleInfo;
	HistInfo* histInfo;

	TH1Ext()
	{
		TH1::AddDirectory( kFALSE );

		this->h = NULL;

		this->hasRatio = kFALSE;
		this->h_ratio = NULL;
	}

	TH1Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString histName = "" ): TH1Ext()
	{
		this->sampleInfo = _sampleInfo;
		this->histInfo = _histInfo;

		if( histName == "" )
			this->h = Get_Hist( this->sampleInfo->fileName, this->histInfo->name );
		else
			this->h = Get_Hist( this->sampleInfo->fileName, histName );

		// -- it should be done earlier: to be consistent with the ratio calculation -- //
		if( this->histInfo->hasRebinX )
			h->Rebin( this->histInfo->nRebinX );
	}

	TH1Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TH1D* _h ): TH1Ext()
	{
		this->sampleInfo = _sampleInfo;
		this->histInfo = _histInfo;

		this->h = (TH1D*)_h->Clone();

		// -- it should be done earlier: to be consistent with the ratio calculation -- //
		if( this->histInfo->hasRebinX )
			h->Rebin( this->histInfo->nRebinX );
	}

	void DrawAndSet( TString drawOp )
	{
		this->h->Draw( drawOp );
		this->SetAttributes(); // -- setting after drawing: to be consistent with TGraphExt case -- //
	}

	void DrawRatioAndSet( TString DrawOp, TString ratioTitle, Double_t minRatio = 0.5, Double_t maxRatio = 1.5 )
	{
		this->h_ratio->Draw( DrawOp );
		this->SetAttributesRatio(ratioTitle, minRatio, maxRatio);
	}

	void AddToLegend( TLegend *legend )
	{
		legend->AddEntry( this->h, this->sampleInfo->fullName );
	}

	void CalcRatio_DEN( TH1D* h_DEN )
	{
		this->hasRatio = kTRUE;

		if( h == NULL )
		{
			cout << "Histogram is not assigned yet!" << endl;
			return;
		}

		h->Sumw2();
		h_DEN->Sumw2();

		this->h_ratio = (TH1D*)this->h->Clone();
		h_ratio->Divide( this->h, h_DEN );
	}

	void CalcRatio_NUM( TH1D* h_NUM )
	{
		this->hasRatio = kTRUE;

		if( h == NULL )
		{
			cout << "Histogram is not assigned yet!" << endl;
			return;
		}

		h->Sumw2();
		h_NUM->Sumw2();

		this->h_ratio = (TH1D*)this->h->Clone();
		h_ratio->Divide( h_NUM, this->h );
	}

protected:
	void SetAttributes()
	{
		this->h->SetTitle("");
		this->h->SetStats(kFALSE);

		this->h->SetLineColor( this->sampleInfo->color );
		this->h->SetFillColorAlpha( kWhite, 0 );
		this->h->SetMarkerStyle( 20 );
		this->h->SetMarkerColor( this->sampleInfo->color );
		if( this->histInfo->isFilled )
		{
			this->h->SetLineColorAlpha( kWhite, 0 );
			this->h->SetMarkerColorAlpha( kWhite, 0 );
			this->h->SetFillColorAlpha( this->sampleInfo->color, 1 );
		}

		if( this->histInfo->hasXRange )
			h->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );

		if( this->histInfo->hasYRange )
			h->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );

		if( this->histInfo->hasZRange )
			h->GetZaxis()->SetRangeUser( this->histInfo->minZ, this->histInfo->maxZ );

		if( this->hasRatio )
			SetAxis_TopPad( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleY );
		else
			SetAxis_SinglePad( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
	}

	void SetAttributesRatio(TString ratioTitle, Double_t minRatio, Double_t maxRatio)
	{
		if( this->h_ratio == NULL ) return;

		this->h_ratio->SetTitle("");
		this->h_ratio->SetStats(kFALSE);

		this->h_ratio->SetLineColor( this->sampleInfo->color );
		this->h_ratio->SetFillColorAlpha( kWhite, 0 );
		this->h_ratio->SetMarkerStyle( 20 );
		this->h_ratio->SetMarkerColor( this->sampleInfo->color );

		SetAxis_BottomPad( this->h_ratio->GetXaxis(), this->h_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
	}
};

// -- TH2 extension -- //
class TH2Ext
{
public:
	TH2D* h;

	SampleInfo* sampleInfo;
	HistInfo* histInfo;

	TH2Ext()
	{
		TH1::AddDirectory( kFALSE );

		this->h = NULL;
	}

	TH2Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString histName = "" ): TH2Ext()
	{
		this->sampleInfo = _sampleInfo;
		this->histInfo = _histInfo;

		if( histName == "" )
			this->h = Get_Hist_2D( this->sampleInfo->fileName, this->histInfo->name );
		else
			this->h = Get_Hist_2D( this->sampleInfo->fileName, histName );
	}

	// -- when the histogram already exists -- //
	TH2Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TH2D* _h ): TH2Ext()
	{
		this->sampleInfo = _sampleInfo;
		this->histInfo = _histInfo;
		this->h = (TH2D*)_h->Clone();
	}

	void DrawAndSet( TString drawOp )
	{
		this->h->Draw( drawOp );
		this->SetAttributes( drawOp ); // -- setting after drawing: to be consistent with TGraphExt case -- //
	}

	void AddToLegend( TLegend *legend )
	{
		legend->AddEntry( this->h, this->sampleInfo->fullName );
	}

protected:
	void SetAttributes( TString drawOp )
	{
		this->h->SetTitle("");
		this->h->SetStats(kFALSE);

		if( drawOp.Contains("SCAT") )
		{
			this->h->SetMarkerStyle( 20 );
			this->h->SetLineColorAlpha( kWhite, 0 );
			this->h->SetMarkerColor( this->sampleInfo->color );
		}

		if( this->histInfo->hasXRange )
			h->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );

		if( this->histInfo->hasYRange )
			h->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );

		if( this->histInfo->hasZRange )
			h->GetZaxis()->SetRangeUser( this->histInfo->minZ, this->histInfo->maxZ );

		if( this->histInfo->hasRebinX )
			h->RebinX( this->histInfo->nRebinX );

		if( this->histInfo->hasRebinY )
			h->RebinY( this->histInfo->nRebinY );

		SetAxis_2D( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
	}
};

class TGraphExt
{
public:
	TGraphAsymmErrors* g;

	Bool_t hasRatio;
	TGraphAsymmErrors* g_ratio;

	SampleInfo* sampleInfo;
	HistInfo* histInfo;

	TGraphExt()
	{
		this->g = NULL;

		this->hasRatio = kFALSE;
		this->g_ratio = NULL;
	}

	TGraphExt( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString graphName = "" ): TGraphExt()
	{
		this->sampleInfo = _sampleInfo;
		this->histInfo = _histInfo;

		if( graphName == "" )
			this->g = Get_Graph( this->sampleInfo->fileName, this->histInfo->name );
		else
			this->g = Get_Graph( this->sampleInfo->fileName, graphName );
	}

	void DrawAndSet( TString drawOp )
	{
		this->g->Draw( drawOp );
		this->SetAttributes();
	}

	void CalcRatio_DEN( TGraphAsymmErrors* g_DEN )
	{
		this->g_ratio = this->MakeRatioGraph( g, g_DEN );
	}

	void CalcRatio_NUM( TGraphAsymmErrors* g_NUM )
	{
		this->g_ratio = this->MakeRatioGraph( g_NUM, g );
	}

	void DrawRatioAndSet( TString DrawOp, TString ratioTitle, Double_t minRatio = 0.5, Double_t maxRatio = 1.5 )
	{
		this->g_ratio->Draw( DrawOp );
		this->SetAttributesRatio(ratioTitle, minRatio, maxRatio);
	}

	void AddToLegend( TLegend *legend )
	{
		legend->AddEntry( this->g, this->sampleInfo->fullName );
	}

protected:
	void SetAttributes()
	{
		this->g->SetTitle("");
		this->g->SetLineColor( this->sampleInfo->color );
		this->g->SetLineWidth( 1 );
		this->g->SetMarkerStyle( 20 );
		this->g->SetMarkerSize( 1.3 );
		this->g->SetMarkerColor( this->sampleInfo->color );
		this->g->SetFillColorAlpha( kWhite, 0 );
		this->g->GetXaxis()->SetTitleFont(42);
		this->g->GetYaxis()->SetTitleFont(42);
		this->g->GetXaxis()->SetLabelFont(42);
		this->g->GetYaxis()->SetLabelFont(42);
		this->g->GetXaxis()->SetLabelSize(0.04);
		this->g->GetYaxis()->SetLabelSize(0.04);

		if( this->histInfo->hasXRange )
			g->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );

		if( this->histInfo->hasYRange )
			g->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );

		if( this->hasRatio )
			SetAxis_TopPad( this->g->GetXaxis(), this->g->GetYaxis(), this->histInfo->titleY );
		else
			SetAxis_SinglePad( this->g->GetXaxis(), this->g->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
	}

	void SetAttributesRatio(TString ratioTitle, Double_t minRatio, Double_t maxRatio)
	{
		if( this->g_ratio == NULL ) return;

		this->g_ratio->SetTitle( "" );

		this->g_ratio->SetLineColor( this->sampleInfo->color );
		this->g_ratio->SetMarkerStyle( 20 );
		this->g_ratio->SetMarkerColor( this->sampleInfo->color );
		this->g_ratio->SetFillColorAlpha( kWhite, 0 );

		SetAxis_BottomPad( this->g_ratio->GetXaxis(), this->g_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
	}

	TGraphAsymmErrors* MakeRatioGraph(TGraphAsymmErrors *g_Type1, TGraphAsymmErrors *g_Type2)
	{
		g_ratio = (TGraphAsymmErrors*)g_Type2->Clone();
		g_ratio->Set(0); // --Remove all points (reset) -- //

		Int_t nPoint = g_Type1->GetN();
		Int_t nPoint_2 = g_Type2->GetN();
		if( nPoint != nPoint_2 )
			printf("# points is different bewteen two graph...be careful for the ratio plot\n");

		for(Int_t i_p=0; i_p<nPoint; i_p++)
		{
			// cout << i_p << "th Point" << endl;
			//Get Type1 point
			Double_t x_Type1, y_Type1;
			g_Type1->GetPoint(i_p, x_Type1, y_Type1);
			Double_t error_Type1 = this->ReturnLargerValue( g_Type1->GetErrorYhigh(i_p), g_Type1->GetErrorYlow(i_p) );
			// cout << "x_Type1: " << x_Type1 << " y_Type1: " << y_Type1 << " error_Type1: " << error_Type1 << " g_Type1->GetErrorYhigh: " << g_Type1->GetErrorYhigh(i_p) << " g_Type1->GetErrorYlow: " << g_Type1->GetErrorYlow(i_p) << endl;

			//Get Type2 point
			Double_t x_Type2, y_Type2;
			g_Type2->GetPoint(i_p, x_Type2, y_Type2);
			Double_t error_Type2 = this->ReturnLargerValue( g_Type2->GetErrorYhigh(i_p), g_Type2->GetErrorYlow(i_p) );
			// cout << "x_Type2: " << x_Type2 << " y_Type2: " << y_Type2 << " error_Type2: " << error_Type2 << " g_Type2->GetErrorYhigh: " << g_Type2->GetErrorYhigh(i_p) << " g_Type2->GetErrorYlow: " << g_Type2->GetErrorYlow(i_p) << endl;

			Double_t ratio;
			Double_t ratio_error;
			if( (nPoint != nPoint_2) && i_p >= nPoint_2)
			{
				ratio = 0;
				ratio_error = 0;
			}
			// else if(y_Type1 != 0 && error_Type1 != 0 && y_Type2 != 0 && error_Type2 != 0)
			else if(y_Type2 != 0)
			{
				ratio = y_Type1 / y_Type2;
				ratio_error = this->Error_PropagatedAoverB(y_Type1, error_Type1, y_Type2, error_Type2);
				//calculate Scale Factor(Type1/Type2) & error

				// cout << "ratio: " << ratio << " ratio_error: " << ratio_error << endl;
			}
			else
			{
				cout << "Denominator is 0! ... ratio and its error are set as 0" << endl;
				ratio = 0;
				ratio_error = 0;
			}

			//Set Central value
			g_ratio->SetPoint(i_p, x_Type1, ratio);

			//Set the error
			Double_t error_XLow = g_Type1->GetErrorXlow(i_p);
			Double_t error_Xhigh = g_Type1->GetErrorXhigh(i_p);
			g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);

			// cout << endl;
		}

		return g_ratio;
	}

	Double_t Error_PropagatedAoverB(Double_t A, Double_t sigma_A, Double_t B, Double_t sigma_B)
	{
		Double_t ratio_A = (sigma_A) / A;
		Double_t ratio_B = (sigma_B) / B;

		Double_t errorSquare = ratio_A * ratio_A + ratio_B * ratio_B;

		return (A/B) * sqrt(errorSquare);
	}

	Double_t ReturnLargerValue(Double_t a, Double_t b)
	{
		if( a > b )
			return a;
		else
			return b;
	}
};