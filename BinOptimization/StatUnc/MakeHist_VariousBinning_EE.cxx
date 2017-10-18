#include <BinOptimization/StatUnc/MakeHist_VariousBinning.h>

void MakeHist_VariousBinning_EE( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Set_ChannelType( "EE" );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_VariousBinning_EE.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}