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

// -- for Rochester Muon momentum correction -- //
#include "./etc/RoccoR/RoccoR.cc"

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

// -- Muon Channel -- //
void MuMu_HistMaker ( TString type = "", TString HLTname = "IsoMu24_OR_IsoTkMu24" )
{
    if ( !type.Length() )
    {
        cout << "Error: no type specified!" << endl;
        return;
    }

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;
    vector<SelProc_t> Processes = Mgr.FindProc( type );

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    for ( Int_t i_proc; i_proc<((int)(Processes.size())); i_proc++ )
    {
        if ( Processes[i_proc] > _EndOf_MuMu )
        {
            cout << "Error: process " << Mgr.Procname[Processes[i_proc]] << " is not MuMu!" << endl;
            continue;
        }

        Mgr.GetProc( Processes[i_proc] );

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "DATA location: " << Mgr.BaseLocation << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For Rochester correction -- //
        TRandom3 *r1 = new TRandom3(0);

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- Output ROOTFile -- //
        TString OutputDir = Mgr.HistLocation;
        TFile *f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedMuMu_t *MuMu = new SelectedMuMu_t();
            MuMu->CreateFromChain( chain );

            RoccoR rc("./etc/RoccoR/rcdata.2016.v3");

            // -- Making Histogram -- //
            TH1D *h_mass_fine_before_PUCorr = new TH1D( "h_mass_fine_before_PUCorr_"+Mgr.Tag[i_tup], "", 10000, 0, 10000 );
            TH1D *h_mass_fine_before_RoccoR = new TH1D( "h_mass_fine_before_RoccoR_"+Mgr.Tag[i_tup], "", 10000, 0, 10000 );
            TH1D *h_mass_fine_before_EffCorr = new TH1D( "h_mass_fine_before_EffCorr_"+Mgr.Tag[i_tup], "", 10000, 0, 10000 );
            TH1D *h_mass_fine = new TH1D( "h_mass_fine_"+Mgr.Tag[i_tup], "", 10000, 0, 10000 );
            TH1D *h_mass = new TH1D( "h_mass_"+Mgr.Tag[i_tup], "", 43, massbins );
            TH1D *h_Pt = new TH1D( "h_Pt_"+Mgr.Tag[i_tup], "", 300, 0, 600 );
            TH1D *h_rapi = new TH1D( "h_rapi_"+Mgr.Tag[i_tup], "", 100, -5, 5 );

            TH1D *h_nPU_beforePUCorr = new TH1D( "h_nPU_before_PUCorr"+Mgr.Tag[i_tup], "", 50, 0, 50 );
            TH1D *h_nPU_beforeEffCorr = new TH1D( "h_nPU_before_EffCorr"+Mgr.Tag[i_tup], "", 50, 0, 50 );
            TH1D *h_nPU = new TH1D( "h_nPU"+Mgr.Tag[i_tup], "", 50, 0, 50 );

            TH1D *h_pT = new TH1D( "h_pT_"+Mgr.Tag[i_tup], "", 300, 0, 600 );
            TH1D *h_eta = new TH1D( "h_eta_"+Mgr.Tag[i_tup], "", 100, -5, 5 );
            TH1D *h_phi = new TH1D( "h_phi_"+Mgr.Tag[i_tup], "", 100, -5, 5 );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Total events: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {
                MuMu->GetEvent(i);

                Double_t GenWeight = MuMu->GENEvt_weight;

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( MuMu->nPileUp );

                // -- efficiency weights -- //
                Double_t weight1 = 0, weight2 = 0, effweight = 1;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                // -- Rochester correction -- //
                for ( Int_t iter = 0; i<=1; i++ )
                {
                    Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                    Int_t s, m;

                    if( Mgr.Type == "DATA" )
                            SF = rc.kScaleDT(MuMu->Muon_charge->at(iter), MuMu->Muon_TuneP_pT->at(iter), MuMu->Muon_TuneP_eta->at(iter),
                                             MuMu->Muon_TuneP_phi->at(iter), s=0, m=0);
                    else
                            SF = rc.kScaleAndSmearMC(MuMu->Muon_charge->at(iter), MuMu->Muon_TuneP_pT->at(iter),
                                                     MuMu->Muon_TuneP_eta->at(iter), MuMu->Muon_TuneP_phi->at(iter),
                                                     MuMu->Muon_trackerLayers->at(iter), rndm[0], rndm[1], s=0, m=0);

                    MuMu->Muon_TuneP_pT->at(iter) = SF*MuMu->Muon_TuneP_pT->at(iter);
                }

                if( MuMu->isSelPassed == kTRUE )
                {
                    TLorentzVector mu1, mu2;
                    mu1.SetPtEtaPhiE( MuMu->Muon_TuneP_pT->at(0), MuMu->Muon_TuneP_eta->at(0), MuMu->Muon_TuneP_phi->at(0), MuMu->Muon_Energy->at(0) );
                    mu2.SetPtEtaPhiE( MuMu->Muon_TuneP_pT->at(1), MuMu->Muon_TuneP_eta->at(1), MuMu->Muon_TuneP_phi->at(1), MuMu->Muon_Energy->at(1) );
                    // -- Apply efficiency scale factor -- //
                    if( Mgr.isMC == kTRUE )
                    {
                        weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF( MuMu );
                        weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH( MuMu );
                        effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;
                    }

                    Double_t reco_Pt = ( mu1 + mu2 ).Pt();
                    Double_t reco_rapi = ( mu1 + mu2 ).Rapidity();

                    h_mass_fine_before_PUCorr->Fill( MuMu->Muon_InvM, TotWeight );
                    h_mass_fine_before_EffCorr->Fill( MuMu->Muon_InvM, TotWeight * PUWeight );
                    h_mass_fine->Fill( MuMu->Muon_InvM, TotWeight * PUWeight * effweight );
                    h_mass->Fill( MuMu->Muon_InvM, TotWeight * PUWeight * effweight );
                    h_Pt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
                    h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

                    h_nPU_beforePUCorr->Fill( MuMu->nPileUp, TotWeight );
                    h_nPU_beforeEffCorr->Fill( MuMu->nPileUp, TotWeight * PUWeight );
                    h_nPU->Fill( MuMu->nPileUp, TotWeight * PUWeight * effweight );

                    h_pT->Fill( MuMu->Muon_TuneP_pT->at(0), TotWeight * PUWeight * effweight );
                    h_pT->Fill( MuMu->Muon_TuneP_pT->at(1), TotWeight * PUWeight * effweight );
                    h_eta->Fill( MuMu->Muon_TuneP_eta->at(0), TotWeight * PUWeight * effweight );
                    h_eta->Fill( MuMu->Muon_TuneP_eta->at(1), TotWeight * PUWeight * effweight );
                    h_phi->Fill( MuMu->Muon_TuneP_phi->at(0), TotWeight * PUWeight * effweight );
                    h_phi->Fill( MuMu->Muon_TuneP_phi->at(1), TotWeight * PUWeight * effweight );

                }// End of event selection
                bar.Draw(i);

            }// End of event iteration

            f->cd();
            cout << "\tWriting into file...";

            h_mass_fine_before_PUCorr->Write();
            h_mass_fine_before_RoccoR->Write();
            h_mass_fine_before_EffCorr->Write();
            h_mass_fine->Write();
            h_mass->Write();
            h_Pt->Write();
            h_rapi->Write();

            h_nPU_beforePUCorr->Write();
            h_nPU_beforeEffCorr->Write();
            h_nPU->Write();

            h_pT->Write();
            h_eta->Write();
            h_phi->Write();

            cout << " Finished." << endl;
            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

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
