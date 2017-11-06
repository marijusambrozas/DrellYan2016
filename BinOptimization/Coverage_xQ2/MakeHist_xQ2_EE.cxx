#include <BinOptimization/Coverage_xQ2/MakeHist_xQ2.h>

void MakeHist_xQ2_EE( TString txtFileName, TString tag, Int_t isMC, Double_t normFactor)
{
	HistProducer *producer = new HistProducer(txtFileName, tag, isMC);

	producer->SetNormFactor( normFactor );
	producer->SetChannelType( "EE" );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_xQ2_EE.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}