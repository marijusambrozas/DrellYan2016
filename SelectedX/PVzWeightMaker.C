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

#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"

void PVzWeightMaker ()
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;

    TFile *f = new TFile( "./etc/PVzWeights.root", "RECREATE" );

    TH1D* h_PVz_MC_ee = new TH1D( "h_PVz_MC_ee", "", 80, -20, 20 );
    TH1D* h_PVz_Data_ee = new TH1D( "h_PVz_Data_ee", "", 80, -20, 20 );
    TH1D* h_PVz_MC_mumu = new TH1D( "h_PVz_MC_mumu", "", 80, -20, 20 );
    TH1D* h_PVz_Data_mumu = new TH1D( "h_PVz_Data_mumu", "", 80, -20, 20 );
    TH1D* h_PVz_MC_emu = new TH1D( "h_PVz_MC_emu", "", 80, -20, 20 );
    TH1D* h_PVz_Data_emu = new TH1D( "h_PVz_Data_emu", "", 80, -20, 20 );


///////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------ EE PART -----------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////////

// -------------------------------------- MC --------------------------------------------//
/*    for ( SelProc_t pr=_EE_DY_10to50; pr<_EndOf_EE_WJets; pr=Next(pr) )
    {
        Mgr.SetProc( pr );
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        DYAnalyzer *analyzer = new DYAnalyzer( "Ele23Ele12" );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_electron();        

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEE_t *EE = new SelectedEE_t();
            EE->CreateFromChain( chain );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights:: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for ( Int_t i=0; i<NEvents; i++ )
            {
                EE->GetEvent(i);

                if( EE->isSelPassed == kTRUE )
                {
                    // -- Normalization -- //
                    Double_t TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * EE->GENEvt_weight;

                    // -- Pileup-Reweighting -- //
                    Double_t PUWeight = analyzer->PileUpWeightValue_80X( EE->nPileUp );

                    // -- efficiency weights -- //
                    Double_t effweight = analyzer->EfficiencySF_EventWeight_electron( EE );

                    // -- Filling -- //
                    h_PVz_MC_ee->Fill( EE->PVz, TotWeight * PUWeight * effweight );
                } // End of event selection
                bar.Draw(i);

            } // End of event iteration
        }// End of i_tup iteration
    }// End of process iteration


// ------------------------------------- Data -------------------------------------------//
    for ( SelProc_t pr=_EE_DoubleEG_B; pr<_EndOf_EE_DoubleEG_Normal; pr=Next(pr) )
    {
        Mgr.SetProc( pr );
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for( Int_t i_tup = 0; i_tup<Ntup; i_tup++ )
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEE_t *EE = new SelectedEE_t();
            EE->CreateFromChain( chain );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for ( Int_t i=0; i<NEvents; i++ )
            {
                EE->GetEvent(i);
                if( EE->isSelPassed == kTRUE )
                    h_PVz_Data_ee->Fill( EE->PVz );
                bar.Draw(i);
            } // End of event iteration
        }// End of i_tup iteration
    }// End of process iteration
*/

///////////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------- MuMu PART ----------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////////

// -------------------------------------- MC --------------------------------------------//

    for ( SelProc_t pr=_MuMu_DY_10to50; pr<_EndOf_MuMu_WJets; pr=next(pr) )
    {
        Mgr.SetProc(pr);

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        DYAnalyzer *analyzer = new DYAnalyzer( "IsoMu24_OR_IsoTkMu24" );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for( Int_t i_tup = 0; i_tup<Ntup; i_tup++ )
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedMuMu_t *MuMu = new SelectedMuMu_t();
            MuMu->CreateFromChain( chain );         

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for( Int_t i=0; i<NEvents; i++ )
            {              
                MuMu->GetEvent(i);

                if( MuMu->isSelPassed == kTRUE )
                {
                    // -- Normalization -- //
                    Double_t TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * MuMu->GENEvt_weight;

                    // -- Pileup-Reweighting -- //
                    Double_t PUWeight = analyzer->PileUpWeightValue_80X( MuMu->nPileUp );

                    // -- efficiency weights -- //
                    Double_t weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF( MuMu );
                    Double_t weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH( MuMu );
                    Double_t effweight = ( L_B2F * weight1 + L_G2H * weight2 ) / L_B2H;

                    h_PVz_MC_mumu->Fill( MuMu->PVz, TotWeight * PUWeight * effweight );
                }
                bar.Draw(i);
            }// End of event iteration
        }// End of i_tup iteration
    } // End of process iteration


// ------------------------------------- Data -------------------------------------------//
    for ( SelProc_t pr=_MuMu_SingleMuon_B; pr<_EndOf_MuMu_Data_Normal; pr=next(pr) )
    {
        Mgr.SetProc(pr);

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for( Int_t i_tup = 0; i_tup<Ntup; i_tup++ )
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedMuMu_t *MuMu = new SelectedMuMu_t();
            MuMu->CreateFromChain( chain );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for( Int_t i=0; i<NEvents; i++ )
            {
                MuMu->GetEvent(i);
                if( MuMu->isSelPassed == kTRUE )
                    h_PVz_Data_mumu->Fill( MuMu->PVz );
                bar.Draw(i);
            }// End of event iteration
        }// End of i_tup iteration
    } // End of process iteration


///////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------ EMu PART ----------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////////

// -------------------------------------- MC --------------------------------------------//

    for ( SelProc_t pr=_EMu_DYTauTau_10to50; pr<_EndOf_EMu_WJets; pr=next(pr) )
    {
        Mgr.SetProc( pr );

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        DYAnalyzer *analyzer = new DYAnalyzer( "IsoMu24_OR_IsoTkMu24" );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();
        analyzer->SetupEfficiencyScaleFactor_electron();

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEMu_t *EMu = new SelectedEMu_t();
            EMu->CreateFromChain( chain );        

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {
                EMu->GetEvent(i);                           

                if( EMu->isSelPassed == kTRUE )
                {
                    // -- Normalization -- //
                    Double_t TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * EMu->GENEvt_weight;

                    // -- Pileup-Reweighting -- //
                    Double_t PUWeight = analyzer->PileUpWeightValue_80X( EMu->nPileUp );

                    // -- efficiency weights -- //
                    Double_t weight1 = analyzer->EfficiencySF_EventWeight_emu_BtoF( EMu );
                    Double_t weight2 = analyzer->EfficiencySF_EventWeight_emu_GtoH( EMu );
                    Double_t effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;

                    h_PVz_MC_emu->Fill( EMu->PVz, TotWeight * PUWeight * effweight );
                }// End of event selection
                bar.Draw(i);

            }// End of event iteration
        }// End of i_tup iteration      
    }// End of process iteration

// ------------------------------------- Data -------------------------------------------//

    for ( SelProc_t pr=_EMu_SingleMuon_B; pr<_EndOf_EMu_Data_Normal; pr=next(pr) )
    {
        Mgr.SetProc( pr );

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEMu_t *EMu = new SelectedEMu_t();
            EMu->CreateFromChain( chain );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {
                EMu->GetEvent(i);
                if( EMu->isSelPassed == kTRUE )
                    h_PVz_Data_emu->Fill( EMu->PVz );
                bar.Draw(i);
            }// End of event iteration
        }// End of i_tup iteration
    }// End of process iteration

// ------------------------------------ FINAL STEPS -------------------------------------//

    // -- Scaling histograms to have integral equal to 1 -- //
//    h_PVz_MC_ee->Scale( 1/h_PVz_MC_ee->Integral(1, h_PVz_MC_ee->GetSize()-2) );
//    h_PVz_Data_ee->Scale( 1/h_PVz_Data_ee->Integral(1, h_PVz_Data_ee->GetSize()-2) );
    h_PVz_MC_mumu->Scale( 1/h_PVz_MC_mumu->Integral(1, h_PVz_MC_mumu->GetSize()-2) );
    h_PVz_Data_mumu->Scale( 1/h_PVz_Data_mumu->Integral(1, h_PVz_Data_mumu->GetSize()-2) );
    h_PVz_MC_emu->Scale( 1/h_PVz_MC_emu->Integral(1, h_PVz_MC_emu->GetSize()-2) );
    h_PVz_Data_emu->Scale( 1/h_PVz_Data_emu->Integral(1, h_PVz_Data_emu->GetSize()-2) );

    // -- Creating weight histograms -- //
//    TH1D* h_PVzWeights_ee = ( (TH1D*)(h_PVz_Data_ee->Clone("h_PVzWeights_ee")) );
    TH1D* h_PVzWeights_mumu = ( (TH1D*)(h_PVz_Data_mumu->Clone("h_PVzWeights_mumu")) );
    TH1D* h_PVzWeights_emu = ( (TH1D*)(h_PVz_Data_emu->Clone("h_PVzWeights_emu")) );

//    h_PVzWeights_ee->Divide( h_PVz_MC_ee );
    h_PVzWeights_mumu->Divide( h_PVz_MC_mumu );
    h_PVzWeights_emu->Divide( h_PVz_MC_emu );

    // -- Drawing -- //
//    TCanvas *c_ee = new TCanvas("ee", "ee", 800, 800);
//    h_PVz_MC_ee->SetMarkerStyle( kFullDotLarge );
//    h_PVz_Data_ee->SetMarkerStyle( kFullDotLarge );
//    h_PVz_MC_ee->SetMarkerColor( kRed );
//    h_PVz_Data_ee->SetMarkerStyle( kBlack );
//    h_PVz_MC_ee->Draw();
//    h_PVz_Data_ee->Draw("SAME");
//    c_ee->Update();

    TCanvas *c_mumu = new TCanvas("mumu", "mumu", 800, 800);
    h_PVz_MC_mumu->SetMarkerStyle( kFullDotLarge );
    h_PVz_Data_mumu->SetMarkerStyle( kFullDotLarge );
    h_PVz_MC_mumu->SetMarkerColor( kRed );
    h_PVz_Data_mumu->SetMarkerStyle( kBlack );
    h_PVz_MC_mumu->Draw();
    h_PVz_Data_mumu->Draw("SAME");
    c_mumu->Update();

    TCanvas *c_emu = new TCanvas("emu", "emu", 800, 800);
    h_PVz_MC_emu->SetMarkerStyle( kFullDotLarge );
    h_PVz_Data_emu->SetMarkerStyle( kFullDotLarge );
    h_PVz_MC_emu->SetMarkerColor( kRed );
    h_PVz_Data_emu->SetMarkerStyle( kBlack );
    h_PVz_MC_emu->Draw();
    h_PVz_Data_emu->Draw("SAME");
    c_emu->Update();

    TCanvas *c_weights = new TCanvas("weights", "weights", 800, 800);
//    h_PVzWeights_ee->SetMarkerStyle( kFullDotLarge );
    h_PVzWeights_mumu->SetMarkerStyle( kFullDotLarge );
    h_PVzWeights_emu->SetMarkerStyle( kFullDotLarge );
//    h_PVzWeights_ee->SetMarkerColor( kBlack );
    h_PVzWeights_mumu->SetMarkerColor( kRed );
    h_PVzWeights_emu->SetMarkerColor( kBlue );
//    h_PVzWeights_ee->Draw();
    h_PVzWeights_mumu->Draw();
    h_PVzWeights_emu->Draw("SAME");
    c_weights->Update();

//    h_PVzWeights_ee->SetDirectory(0);
    h_PVzWeights_mumu->SetDirectory(0);
    h_PVzWeights_emu->SetDirectory(0);
//    h_PVz_MC_ee->SetDirectory(0);
//    h_PVz_Data_ee->SetDirectory(0);
    h_PVz_MC_mumu->SetDirectory(0);
    h_PVz_Data_mumu->SetDirectory(0);
    h_PVz_MC_emu->SetDirectory(0);
    h_PVz_Data_emu->SetDirectory(0);

    // -- Writing -- //
    f->cd();
//    h_PVz_MC_ee->Write();
//    h_PVz_Data_ee->Write();
    h_PVz_MC_mumu->Write();
    h_PVz_Data_mumu->Write();
    h_PVz_MC_emu->Write();
    h_PVz_Data_emu->Write();
//    h_PVzWeights_ee->Write();
    h_PVzWeights_mumu->Write();
    h_PVzWeights_emu->Write();
    f->Close();

    if ( !f->IsOpen() ) cout << "File PVzWeights.root has been closed successfully.\n" << endl;
    else cout << "FILE PVzWeights.root COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of PVzWeightMaker.C
