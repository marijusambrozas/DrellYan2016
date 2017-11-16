#include <BinOptimization/BinTest/MakeHist_FixBinning.h>

void MakeHist_FixBinnning_MuMu( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Set_ChannelType( "MuMu" );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_FixBinning_MuMu.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}