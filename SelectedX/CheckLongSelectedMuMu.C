#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TAttMarker.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>

// -- Macro to check if selected and reselected events are the same -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "header/myProgressBar_t.h"

void CheckLongSelectedMuMu(Bool_t DrawHistos = kTRUE)
{
	// -- Run2016 luminosity [/pb] -- //
        Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0; // Need checking
	L = L_B2H;

        TTimeStamp ts_start;
        cout << "\n[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

	TStopwatch totaltime;
	totaltime.Start();

	// -- Each ntuple directory & corresponding Tags -- //
//	vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
        vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
        Tag.push_back("ZMuMu_M4500to6000");
        nEvents.push_back(100000);  // Primary dataset actually contained 19672 events.
        Xsec.push_back(4.56E-07);   // Needs checking
        const Int_t Ntup = 1;
	for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
	{
            TStopwatch looptime;
            looptime.Start();

//		cout << "\t<" << [i_tup] << ">" << endl;

            TChain *ch_s = new TChain("DYTree");
            TChain *ch_re = new TChain("DYTree");
            ch_s->Add("/media/sf_DATA/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");
            ch_re->Add("/media/sf_DATA/ReSelectedMuMu_ZToMuMu_M4500to6000_4.root");

            LongSelectedMuMu_t * MuMu_s = new LongSelectedMuMu_t();
            LongSelectedMuMu_t *MuMu_re = new LongSelectedMuMu_t();
            MuMu_s->CreateFromChain(ch_s);
            MuMu_re->CreateFromChain(ch_re);

            TH1D* h1_pT_s = new TH1D("h1pts", "h1_pT_s", 400, 0, 4000);
            TH1D* h1_pT_re = new TH1D("h1ptre", "h1_pT_re", 400, 0, 4000);
            TH1D* h1_TuneP_pT_s = new TH1D("h1tuneppts", "h1_TuneP_pT_s", 400, 0, 4000);
            TH1D* h1_TuneP_pT_re = new TH1D("h1tunepptre", "h1_TuneP_pT_re", 400, 0, 4000);
            TH1D* h1_eta_s = new TH1D("h1etas", "h1_eta_s", 120, -3, 3);
            TH1D* h1_eta_re = new TH1D("h1etare", "h1_eta_re", 120, -3, 3);
            TH1D* h1_TuneP_eta_s = new TH1D("h1tunepetas", "h1_TuneP_eta_s", 120, -3, 3);
            TH1D* h1_TuneP_eta_re = new TH1D("h1tunepetare", "h1_TuneP_eta_re", 120, -3, 3);
            TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
            TH1D* h1_invm_re = new TH1D("h1invmre", "h1_invm_re", 350, 0, 7000);
            TH1D* h1_trigPhi_s = new TH1D("h1trigphis", "h1_trigPhi_s", 400, -4, 4);
            TH1D* h1_trigPhi_re = new TH1D("h1trigphire", "h1_trigPhi_re", 400, -4, 4);
            TH1D* h1_dzVTX_s = new TH1D("h1dzvtxs", "h1_dzVTX_s", 100, -0.05, 0.05);
            TH1D* h1_dzVTX_re = new TH1D("h1dzvtxre", "h1_dzVTX_re", 100, -0.05, 0.05);
            TH1D* h1_vtxTrkProb_s = new TH1D("h1vtxtrkprobs", "h1_vtxTrkProb_s", 55, 0, 1.1);
            TH1D* h1_vtxTrkProb_re = new TH1D("h1vtxtrkprobre", "h1_vtxTrkProb_re", 55, 0, 1.1);

            Bool_t AllOk = kTRUE;

            Int_t NEvents_s = ch_s->GetEntries();
            Int_t NEvents_re = ch_re->GetEntries();
            cout << "\tTotal selected Events: " << NEvents_s << endl;
            cout << "\tTotal reselected Events: " << NEvents_re << endl;
            if(NEvents_s!=NEvents_re)
            {
                cout << "\tNUMBER IS NOT THE SAME!!" << endl;
                AllOk = kFALSE;
            }
            else
            {
                myProgressBar_t bar(NEvents_s);
                cout << "\tChecking separate events:" << endl;
                for(Int_t i=0; i<NEvents_s; i++)
                {
                    MuMu_s->GetEvent(i);
                    MuMu_re->GetEvent(i);

                    if ( MuMu_s->Muon_pT->size()!=2 || MuMu_re->Muon_pT->size()!=2 )
                    {
                        cout << "Event " << i << ": pT vector size is is not 2!! Breaking.." << endl;
                        AllOk = kFALSE;
                        break;
                    }
                    if( MuMu_s->HLT_trigPhi->size()!=MuMu_re->HLT_trigPhi->size() )
                    {
                        cout << "Event " << i << ": HLT_trigPhi vector sizes do not match!! Breaking.." << endl;
                        AllOk = kFALSE;
                        break;
                    }
                    if ( MuMu_s->evtNum!=MuMu_re->evtNum ) {
                        cout << "Event " << i << ": Event numbers do not match." << endl;
                        AllOk = kFALSE;
                    }
                    else if (AllOk==kTRUE)
                    {
                        if ( MuMu_s->Muon_charge->at(0)!=MuMu_re->Muon_charge->at(0) || MuMu_s->Muon_charge->at(1)!=MuMu_re->Muon_charge->at(1) )
                        {
                            cout << "Event " << i << ": Muon charges do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->Muon_pT->at(0)!=MuMu_re->Muon_pT->at(0) || MuMu_s->Muon_pT->at(1)!=MuMu_re->Muon_pT->at(1) )
                        {
                            cout << "Event " << i << ": Muon_pT do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->Muon_Px->at(0)!=MuMu_re->Muon_Px->at(0) || MuMu_s->Muon_Px->at(1)!=MuMu_re->Muon_Px->at(1) )
                        {
                            cout << "Event " << i << ": Muon_Px do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->Muon_TuneP_phi->at(0)!=MuMu_re->Muon_TuneP_phi->at(0) || MuMu_s->Muon_TuneP_phi->at(1)!=MuMu_re->Muon_TuneP_phi->at(1) )
                        {
                            cout << "Event " << i << ": Muon_TuneP_phi do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->HLT_trigEta->at(0)!=MuMu_re->HLT_trigEta->at(0) ||
                             MuMu_s->HLT_trigEta->at(MuMu_s->HLT_trigEta->size()-1)!=MuMu_re->HLT_trigEta->at(MuMu_re->HLT_trigEta->size()-1) )
                        {
                            cout << "Event " << i << ": HLT_trigEta do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->HLT_ntrig!=MuMu_re->HLT_ntrig )
                        {
                            cout << "Event " << i << ": Trigger numbers do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->CosAngle->at(0)!=MuMu_re->CosAngle->at(0) || MuMu_s->vtxTrkChi2->at(0)!=MuMu_re->vtxTrkChi2->at(0) )
                        {
                            cout << "Event " << i << ": Muon_TuneP_phi do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                        if ( MuMu_s->GENEvt_weight!=MuMu_re->GENEvt_weight )
                        {
                            cout << "Event " << i << ": Event weights do not match!!" << endl;
                            AllOk = kFALSE;
                        }
                    } // End of else if(AllOk)

                    // -- Normalization -- //
//                            Double_t TotWeight_s = MuMu_s->GENEvt_weight;
//                            Double_t TotWeight_re = MuMu_re->GENEvt_weight;
//                            TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                            TotWeight_re = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_re->GENEvt_weight;

                    // -- Histogram filling -- //
                    for (Int_t j=0; j<2; j++)
                    {
                        h1_pT_s->Fill(MuMu_s->Muon_pT->at(j), MuMu_s->GENEvt_weight);
                        h1_pT_re->Fill(MuMu_re->Muon_pT->at(j), MuMu_re->GENEvt_weight);
                        h1_TuneP_pT_s->Fill(MuMu_s->Muon_TuneP_pT->at(j), MuMu_s->GENEvt_weight);
                        h1_TuneP_pT_re->Fill(MuMu_re->Muon_TuneP_pT->at(j), MuMu_re->GENEvt_weight);
                        h1_eta_s->Fill(MuMu_s->Muon_eta->at(j), MuMu_s->GENEvt_weight);
                        h1_eta_re->Fill(MuMu_re->Muon_eta->at(j), MuMu_re->GENEvt_weight);
                        h1_TuneP_eta_s->Fill(MuMu_s->Muon_TuneP_eta->at(j), MuMu_s->GENEvt_weight);
                        h1_TuneP_eta_re->Fill(MuMu_re->Muon_TuneP_eta->at(j), MuMu_re->GENEvt_weight);
                        h1_dzVTX_s->Fill(MuMu_s->Muon_dzVTX->at(j), MuMu_s->GENEvt_weight);
                        h1_dzVTX_re->Fill(MuMu_re->Muon_dzVTX->at(j), MuMu_re->GENEvt_weight);
                    }
                    h1_invm_s->Fill(MuMu_s->Muon_InvM, MuMu_s->GENEvt_weight);
                    h1_invm_re->Fill(MuMu_re->Muon_InvM, MuMu_re->GENEvt_weight);
                    h1_vtxTrkProb_s->Fill(MuMu_s->vtxTrkProb->at(0), MuMu_s->GENEvt_weight);
                    h1_vtxTrkProb_re->Fill(MuMu_re->vtxTrkProb->at(0), MuMu_re->GENEvt_weight);
                    for (UInt_t j=0; j<MuMu_s->HLT_trigPhi->size(); j++)
                    {
                        h1_trigPhi_s->Fill(MuMu_s->HLT_trigPhi->at(j), MuMu_s->GENEvt_weight);
                        h1_trigPhi_re->Fill(MuMu_re->HLT_trigPhi->at(j), MuMu_re->GENEvt_weight);
                    }

                    bar.Draw(i);
                } //End of event iteration
            } //End of if( same event numbers )

            if ( AllOk==kTRUE )
            {
                cout << "\tNo problems found so far." << endl;
            }
            cout << "\tChecking histograms:  ";

            // -- Histogram checking -- //
            for (Int_t i=1; i<401; i++)
            {
                if(h1_pT_s->GetBinContent(i)!=h1_pT_re->GetBinContent(i))
                {
                    cout << "h1_pT bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(h1_TuneP_pT_s->GetBinContent(i)!=h1_TuneP_pT_s->GetBinContent(i))
                {
                    cout << "h1_TuneP_pT bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(h1_trigPhi_s->GetBinContent(i)!=h1_trigPhi_re->GetBinContent(i))
                {
                    cout << "h1_trigPhi bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(i<351 && h1_invm_s->GetBinContent(i)!=h1_invm_re->GetBinContent(i))
                {
                    cout << "h1_invm bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(i<121 && h1_eta_s->GetBinContent(i)!=h1_eta_re->GetBinContent(i))
                {
                    cout << "h1_eta bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(i<121 && h1_TuneP_eta_s->GetBinContent(i)!=h1_TuneP_eta_re->GetBinContent(i))
                {
                    cout << "h1_TuneP_eta bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(i<101 && h1_dzVTX_s->GetBinContent(i)!=h1_dzVTX_re->GetBinContent(i))
                {
                    cout << "h1_dzVTX bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
                if(i<56 && h1_vtxTrkProb_s->GetBinContent(i)!=h1_vtxTrkProb_re->GetBinContent(i))
                {
                    cout << "h1_vtxTrkProb bin " << i << ": Bin contents do not match!! Breaking.." << endl;
                    AllOk = kFALSE;
                    break;
                }
            }// End of for(histogram bins)

            // -- Drawing histos -- //
            if ( DrawHistos==kTRUE )
            {
                TCanvas* c_pT = new TCanvas ("pT", "pT", 1000, 1000);
                TCanvas* c_TuneP_pT = new TCanvas ("TuneP_pT", "TuneP_pT", 1000, 1000);
                TCanvas* c_eta = new TCanvas ("eta", "eta", 1000, 1000);
                TCanvas* c_TuneP_eta = new TCanvas ("TuneP_eta", "TuneP_eta", 1000, 1000);
                TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
                TCanvas* c_trigPhi = new TCanvas ("trigPhi", "trigPhi", 1000, 1000);
                TCanvas* c_dzVTX = new TCanvas ("dzVTX", "dzVTX", 1000, 1000);
                TCanvas* c_vtxTrkProb = new TCanvas ("vtxTrkProb", "vtxTrkProb", 1000, 1000);
                TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

                h1_pT_s->SetLineColor(kBlue);
                h1_pT_s->SetLineWidth(2);
                h1_pT_re->SetLineColor(kRed);
                h1_TuneP_pT_s->SetLineColor(kBlue);
                h1_TuneP_pT_s->SetLineWidth(2);
                h1_TuneP_pT_re->SetLineColor(kRed);
                h1_eta_s->SetLineColor(kBlue);
                h1_eta_s->SetLineWidth(2);
                h1_eta_re->SetLineColor(kRed);
                h1_TuneP_eta_s->SetLineColor(kBlue);
                h1_TuneP_eta_s->SetLineWidth(2);
                h1_TuneP_eta_re->SetLineColor(kRed);
                h1_invm_s->SetLineColor(kBlue);
                h1_invm_s->SetLineWidth(2);
                h1_invm_re->SetLineColor(kRed);
                h1_trigPhi_s->SetLineColor(kBlue);
                h1_trigPhi_s->SetLineWidth(2);
                h1_trigPhi_re->SetLineColor(kRed);
                h1_dzVTX_s->SetLineColor(kBlue);
                h1_dzVTX_s->SetLineWidth(2);
                h1_dzVTX_re->SetLineColor(kRed);
                h1_vtxTrkProb_s->SetLineColor(kBlue);
                h1_vtxTrkProb_s->SetLineWidth(2);
                h1_vtxTrkProb_re->SetLineColor(kRed);
                leg->AddEntry(h1_pT_s, "Selected MuMu", "l");
                leg->AddEntry(h1_pT_re, "Reselected MuMu", "l");

                c_pT->cd(); h1_pT_s->Draw(); h1_pT_re->Draw("SAME"); leg->Draw(); c_pT->Update();
                c_TuneP_pT->cd(); h1_TuneP_pT_s->Draw(); h1_TuneP_pT_re->Draw("SAME"); leg->Draw(); c_TuneP_pT->Update();
                c_eta->cd(); h1_eta_s->Draw(); h1_eta_re->Draw("SAME"); leg->Draw(); c_eta->Update();
                c_TuneP_eta->cd(); h1_TuneP_eta_s->Draw(); h1_TuneP_eta_re->Draw("SAME"); leg->Draw(); c_TuneP_eta->Update();
                c_invm->cd(); h1_invm_s->Draw(); h1_invm_re->Draw("SAME"); leg->Draw(); c_invm->Update();
                c_trigPhi->cd(); h1_trigPhi_s->Draw(); h1_trigPhi_re->Draw("SAME"); leg->Draw(); c_trigPhi->Update();
                c_dzVTX->cd(); h1_dzVTX_s->Draw(); h1_dzVTX_re->Draw("SAME"); leg->Draw(); c_dzVTX->Update();
                c_vtxTrkProb->cd(); h1_vtxTrkProb_s->Draw(); h1_vtxTrkProb_re->Draw("SAME"); leg->Draw(); c_vtxTrkProb->Update();
            }

            if ( AllOk==kTRUE )
            {
                cout << "All bin values match.\nSelected and reselected events match perfectly! Hooray!" << endl;
            }
	} //end of i_tup iteration

	Double_t TotalRunTime = totaltime.CpuTime();
        cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

	TTimeStamp ts_end;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
