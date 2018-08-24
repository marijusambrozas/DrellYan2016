#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

// -- Electron Channel -- //
void EE_HistMaker ( TString type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12" )
{
	// -- Run2016 luminosity [/pb] -- //
	Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
	L = L_B2H;

	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

        LocalFileMgr Mgr;
        vector<SelProc_t> Processes = Mgr.FindProc( type );

        if ( !Processes.size() )
        {
            cout << "Error: no processes!" << endl;
            return;
        }

        for ( Int_t i_proc; i_proc<Processes.size(); i_proc++ )
        {
            Mgr.GetProc( Processes[i_proc] );

            cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
            cout << "Type: " << Mgr.Type << endl;
            cout << "DATA location: " << Mgr.BaseLocation << endl;

            TStopwatch totaltime;
            totaltime.Start();

            DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

            // -- For PU re-weighting -- //
            analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

            // -- For efficiency SF -- //
            analyzer->SetupEfficiencyScaleFactor_electron();

            // -- Output ROOTFile -- //
            TString OutputDir = Mgr.BaseLocation+"SelectedEE/Histos/";
            TFile *f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );

            //Loop for all samples
            const Int_t Ntup = Mgr.FileLocation.size();
            for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
            {
                TStopwatch looptime;
                looptime.Start();

                cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

                TChain *chain = new TChain( "DYTree" );
                chain->Add( Mgr.FullLocation[i_tup] );

                SelectedEE_t *EE = new SelectedEE_t();
                EE->CreateFromChain( chain );

                // -- Making Histogram -- //
                TH1D *h_mass_fine_before_PUCorr = new TH1D("h_mass_fine_before_PUCorr_"+Mgr.Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass_fine_before_EffCorr = new TH1D("h_mass_fine_before_EffCorr_"+Mgr.Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Mgr.Tag[i_tup], "", 10000, 0, 10000);
                TH1D *h_mass = new TH1D("h_mass_"+Mgr.Tag[i_tup], "", 43, massbins);
                TH1D *h_Pt = new TH1D("h_Pt_"+Mgr.Tag[i_tup], "", 300, 0, 600);
                TH1D *h_rapi = new TH1D("h_rapi_"+Mgr.Tag[i_tup], "", 100, -5, 5);

                TH1D *h_pT = new TH1D("h_pT_"+Mgr.Tag[i_tup], "", 300, 0, 600);
                TH1D *h_eta = new TH1D("h_eta_"+Mgr.Tag[i_tup], "", 100, -5, 5);
                TH1D *h_phi = new TH1D("h_phi_"+Mgr.Tag[i_tup], "", 100, -5, 5);

                Int_t NEvents = chain->GetEntries();
                if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
                cout << "\t[Total events: " << Mgr.Wsum << "]" << endl;
                cout << "\t[Selected Events: " << NEvents << "]" << endl;

                myProgressBar_t bar( NEvents );

                for ( Int_t i=0; i<NEvents; i++ )
                {
                    EE->GetEvent(i);

                    // -- Bring the weights -- //
                    Double_t GenWeight = EE->GENEvt_weight;
                    SumWeight += GenWeight;

                    // -- Pileup-Reweighting -- //
                    Double_t PUWeight = 1;
                    if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( ntuple->nPileUp );

                    // -- efficiency weights -- //
                    Double_t effweight = 1;

                    // -- Normalization -- //
                    Double_t TotWeight = GenWeight;
                    if( Mgr.isMC == kTRUE ) TotWeight = ( L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] ) * GenWeight;

                    if( EE->isSelPassed == kTRUE )
                    {
                        TLorentzVector ele1, ele2;
                        ele1.SetPtEtaPhiE( EE->Electron_pT->at(0), EE->Electron_eta->at(0), EE->Electron_phi->at(0), EE->Electron_Energy->at(0) );
                        ele1.SetPtEtaPhiE( EE->Electron_pT->at(1), EE->Electron_eta->at(1), EE->Electron_phi->at(1), EE->Electron_Energy->at(1) );
                        Double_t reco_Pt = ( ele1 + ele2 ).Pt();
                        Double_t reco_rapi = ( ele1 + ele2 ).Rapidity();

                        // -- Apply efficiency correcion -- //
                        if( Mgr.isMC == kTRUE )
                                effweight = analyzer->EfficiencySF_EventWeight_electron( EE );

                        h_mass_fine_before_PUCorr->Fill( EE->Electron_InvM, TotWeight );
                        h_mass_fine_before_EffCorr->Fill( EE->Electron_InvM, TotWeight * PUWeight );
                        h_mass_fine->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight );
                        h_mass->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight );
                        h_Pt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
                        h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

                        h_pT->Fill( EE->Electron_pT->at(0), TotWeight * PUWeight * effweight );
                        h_pT->Fill( EE->Electron_pT->at(1), TotWeight * PUWeight * effweight );
                        h_eta->Fill( EE->Electron_eta->at(0), TotWeight * PUWeight * effweight );
                        h_eta->Fill( EE->Electron_eta->at(1), TotWeight * PUWeight * effweight );
                        h_phi->Fill( EE->Electron_phi->at(0), TotWeight * PUWeight * effweight );
                        h_phi->Fill( EE->Electron_phi->at(1), TotWeight * PUWeight * effweight );

                    } // End of event selection
                    bar.Draw(i);

                } //End of event iteration

                f->cd();
                cout << "\tWriting into file...";

                h_mass_fine_before_PUCorr->Write();
                h_mass_fine_before_EffCorr->Write();
                h_mass_fine->Write();
                h_mass->Write();
                h_Pt->Write();
                h_rapi->Write();

                h_pT->Write();
                h_eta->Write();
                h_phi->Write();

                cout << " Finished." << endl;
                if( Mgr.isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup]);

                Double_t LoopRunTime = looptime.CpuTime();
                cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

            }// End of i_tup iteration

            f->Close();
            if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root has been closed successfully.\n" << endl;
            else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root COULD NOT BE CLOSED!\n" << endl;

        }// End of i_proc iteration

	Double_t TotalRunTime = totaltime.CpuTime();
	cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

	TTimeStamp ts_end;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
