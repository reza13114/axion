#include <stdio.h>
#include "iostream"
#include "TH1F.h"
#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "stdlib.h"
#include "math.h"
#include "TPaveText.h"
#include <string>
#include <vector>
#include "TVirtualFitter.h"
#include "TBranch.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include"THStack.h"
#include"TPad.h"
#include"TLegend.h"
#include <vector>
#include "myclass.C"
//#if not defined(__CINT__) || defined(__MAKECINT__)

// define signal histograms

//    Double_t   sigma_MG [] = { tbr13TeV,  tbi13TeV,    ep2bbjvbkg1,   ep2ccjvbkg2,  ep2vjz2bbbkg3,    ep2vjz2ccbkg4,  ep2vjz2jjbkg5,    ep2jt2w2jjbkg6,   ep2vjh2bbSMbkg7}; 
//OLD: 1.3 TeV sigma_MG [9] = {0.42564,     0.91731,     0.16329,       0.25351,      0.07348,           0.0580833,     0.207184,         0.5007655,          0.07667}; //pb, old
const Double_t sigma_MG [9] = {0.1359855, 0.1188306, 0.1641728379,   2.65558E-01,     7.36076E-02,       5.80564E-02,     2.06976E-01,      5.00882E-01,      1.15423E-01}; //1.3 TeV (pb)
//const Double_t sigma_MG [9] = {6.05582E-01, 5.95730E-01, 6.73908E-01, ???????,     2.91544E-01,       2.30013E-01,     8.19878E-01,      5.17822E+00,      3.92817E-01}; //3.4 TeV (pb)
//const Double_t sigma_MG [9] = {1.80681E-01}, //1.3 TeV (pb) 2-dim signals (TR and TI)
//const Double_t sigma_MG [9] = {6.14823E-01}, //3.4 TeV (pb)  2-dim signals (TR and TI)
Double_t            sigma;
//Double_t            sigma_pre;
Double_t            N_total;
Double_t            N_accepted;
//Double_t            N_accepted_pre;
//Double_t            efficiency_pre;
Double_t            efficiency;

// Arrays for groups
int N_accepted_group[6] = {0};   // 0=hz,1=tt,2=ww,3=zz,4=jets, 5=others
int N_total_group[6] = {0};
Double_t sigma_MG_group [6] = {1, // bkg_hz
                                 1,         // bkg_tt
                                 1,                    // bkg_ww
                                 1,                        // bkg_zz
                                 1,                    // bkg_jets
                                 1.0};                                       // others (signal)

TString group_name[8] = {"bkg_hz", "bkg_ll", "bkg_tautau ", "bkg_ww", "bkg_zz", "bkg_tt", "bkg_za", "bkg_jets"};
int group = 0;


//char *samplename[] = {"xtar", "tbr13TeV"};
//vector  <TString> name = {xtar, tbr13TeV, xtr, xtai, tbi13TeV, xti, ep2bbjvbkg1, ep2ccjvbkg2, ep2vjz2bbbkg3, ep2vjz2ccbkg4, ep2vjz2jjbkg5, ep2jt2w2jjbkg6, ep2vjh2bbSMbkg7};

TH1F *hjet1_pt[9], *hjet2_pt[9], *hjet1_eta[9], *hjet2_eta[9], *hjet1_P[9], *hjet2_P[9];
TH1F *hjetsInvariantMass[9], *hMET[9], *hHT[9], *hdeltaeta_jets[9];
TH1F *hdeltaphi_jets[9], *hdeltaR_jets[9];
TH1F *hcos_jets[9], *hjetmultiplicity[9];

TFile *f;

//----------------------------------------------------- END 

void SG1plots(){ // jets == 2 (no tau-tagged jets) 


//======================================================= sample loop

    //for (int r = 0; r < 18; r++)
    for (int r = 11; r < 16; r++) { // r: sample number

    
        cout << endl << "sample " << r << ": " << endl;

        hjet1_pt[r] = new TH1F("hjet1_pt","", 17, 0 ,340);
        hjet2_pt[r]  = new TH1F("hjet2_pt","", 10, 0 ,200);

        hjet1_eta[r]  = new TH1F("hjet1_eta","", 20, -2 ,6);
        hjet2_eta[r]  = new TH1F("hjet2_eta","", 15, -2 ,4);

        hjet1_P[r]  = new TH1F("hjet1_P","", 40, 0, 3200);
        hjet2_P[r]  = new TH1F("hjet2_P","", 20, 0, 800);

        hjetsInvariantMass[r]  = new TH1F("hjetsInvariantMass","", 70, 0 ,700);
        hMET[r]  = new TH1F("hMET","", 25, 0 ,500);
        hHT[r]  = new TH1F("hHT","", 50, 0 ,1000);

        hdeltaeta_jets[r]  = new TH1F("hdeltaeta_jets","", 25, -5 ,5);
        hdeltaphi_jets[r]  = new TH1F("hdeltaphi_jets","", 20, -8 ,8);
        hdeltaR_jets[r]  = new TH1F("hdeltaR_jets","", 15, 0 ,6);

        hcos_jets[r]  = new TH1F("hcos_jets","", 15, -1.5, 1.5);
        hjetmultiplicity[r]  = new TH1F("hjetmultiplicity","", 30, 0, 30);

        hjet1_pt[r]->SetDirectory(0);
        hjet2_pt[r]->SetDirectory(0);

        hjet1_eta[r]->SetDirectory(0);
        hjet2_eta[r]->SetDirectory(0);

        hjet1_P[r]->SetDirectory(0);
        hjet2_P[r]->SetDirectory(0);

        hjetsInvariantMass[r]->SetDirectory(0);
        hMET[r]->SetDirectory(0);
        hHT[r]->SetDirectory(0);

        hdeltaeta_jets[r]->SetDirectory(0);
        hdeltaphi_jets[r]->SetDirectory(0);
        hdeltaR_jets[r]->SetDirectory(0);
        hcos_jets[r]->SetDirectory(0);
        hjetmultiplicity[r]->SetDirectory(0);

        if (r==0)  f = new TFile("sig1.root","READ");

        if (r==1) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/hz_background/Events/bkg_hz1.root","READ");
        if (r==2) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/hz_background/Events/bkg_hz2.root","READ");
        if (r==3) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/hz_background/Events/bkg_hz3.root","READ");
        if (r==4) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/hz_background/Events/bkg_hz4.root","READ");

        if (r==5) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/l_background/Events/bkg_ll.root","READ");

        if (r==6) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/ta_background/Events/bkg_tautau.root","READ");

        if (r==7) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/w_background/Events/bkg_ww1.root","READ");
        if (r==8) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/w_background/Events/bkg_ww2.root","READ");

        if (r==9) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/z_background/Events/bkg_zz1.root","READ");
        if (r==10) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/z_background/Events/bkg_zz2.root","READ");

        if (r==11) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/t_background/Events/bkg_tt1.root","READ");
        if (r==12) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/t_background/Events/bkg_tt2.root","READ");
        if (r==13) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/t_background/Events/bkg_tt3.root","READ");
        if (r==14) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/t_background/Events/bkg_tt4.root","READ");

        if (r==15) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/za_background/Events/bkg_za.root","READ");
       
        if (r==16) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/bkg_jets/bkg_jets1.root","READ");
        if (r==17) f = new TFile("/cmsdata4/haghighat/axion/MG5_aMC_v2_6_7/bkg_jets/bkg_jets2.root","READ");


        TTree *tree_sample;
        f->GetObject("Delphes",tree_sample);   // Get objects from root file
        myclass     my(tree_sample);   // Define an object for Getting Entries
 /*
        TTree *tree_sample = (TTree*) f->Get( "Delphes" );
        myclass my(tree_sample);
  */ 
        TLorentzVector                   jet_4v;
        vector <TLorentzVector>          jet;

        TLorentzVector                   jet1_4v;
        vector <TLorentzVector>          jet1;


        Double_t                         jet1_pt;
        Double_t                         jet2_pt;

        Double_t                         jet1_eta;
        Double_t                         jet2_eta;

        TVector3                         jet1_Pvec;
        TVector3                         jet2_Pvec;


        Double_t                         jet1_P;
        Double_t                         jet2_P;

        Double_t                         jetsInvariantMass;
        Double_t                         MET;
        Double_t                         HT;

        Double_t                         deltaeta_jets;
        Double_t                         deltaphi_jets;
        Double_t                         deltaR_jets;

        Double_t                         cos_jets;
        Int_t                            jetmultiplicity;

        Double_t                         alpha;
        Double_t                         beta1;
        Double_t                         beta2;

//======================================================= event loop
        N_accepted = 0;
        //N_accepted_pre = 0;

        //for (int i = 0; i < my.fChain->GetEntriesFast(); i++) {
        for (int i = 0; i<1000; i++) { // i: event number

            if (i%200 == 0) std::cout << "sample number " << r << ": ---------- ... Processing event: " << i << std::endl;
                   
            my.GetEntry(i);

            //if (my.Jet_size < 3) continue;

            for (int m = 0; m < my.Jet_size; m++)
            {
                if (/*my.Jet_BTag[m] == 0 &&*/ // Non-BTagged jets
                    my.Jet_PT[m] >= 30 &&
                    fabs(my.Jet_Eta[m]) <= 2.5)
                {
                    jet_4v.SetPtEtaPhiM(my.Jet_PT[m],my.Jet_Eta[m],my.Jet_Phi[m],0);
                    jet.push_back(jet_4v);                    
                }

                // if (my.Jet_TauTag[m] == 1 && // BTag == 1
                //     my.Jet_PT[m] >= 30 && 
                //     fabs(my.Jet_Eta[m]) <= 2.5)
                // { 
                //     tau_4v.SetPtEtaPhiM(my.Jet_PT[m],my.Jet_Eta[m],my.Jet_Phi[m],0);
                //     tau.push_back(tau_4v);
                // }
          }
//======================================= event selection =============================================

            if (jet.size() == 2 &&
                jet.at(0).DeltaR(jet.at(1)) > 0.5)
            {
                //N_accepted_pre++;

                //if (jet[0].Pt() < jet[1].Pt())     swap(jet[0], jet[1]);

                jet1_pt = jet[0].Pt();
                jet2_pt = jet[1].Pt(); 

                jet1_eta = jet[0].Eta();
                jet2_eta = jet[1].Eta(); 

                jetsInvariantMass = (jet.at(0) + jet.at(1)).M();
                MET = my.MissingET_MET[0];
                HT = jet[0].Pt() + jet[0].Pt();

                deltaeta_jets = jet[0].Eta() - jet[1].Eta();
                deltaphi_jets = jet[0].Phi() - jet[1].Phi();
                deltaR_jets = jet.at(0).DeltaR(jet.at(1));

                //jetmultiplicity = my.Jet_size;
                jetmultiplicity = jet.size();

                jet1_Pvec = jet[0].Vect();
                jet2_Pvec = jet[1].Vect(); 
                jet1_P = jet1_Pvec.Mag();
                jet2_P = jet2_Pvec.Mag(); 
/*
		beta1 = acos(( jet[0].Px()*jet[1].Px() +
			       jet[0].Py()*jet[1].Py() + 
			       jet[0].Pz()*jet[1].Pz()  )/(
			  sqrt(jet[0].Px()*jet[0].Px() +
			       jet[0].Py()*jet[0].Py() + 
			       jet[0].Pz()*jet[0].Pz())*
			  sqrt(jet[1].Px()*jet[1].Px() + 
			       jet[1].Py()*jet[1].Py() + 
			       jet[1].Pz()*jet[1].Pz()) ));

		beta2 = acos( (jet_Pvec.Dot(jet2_Pvec)) / (jet_Pvec.Mag() * jet2_Pvec.Mag()) );
*/
		//alpha = jet_Pvec.Angle(jet2_Pvec); // in lab frame 
		
//cout << "alpha = " << alpha << " vs " << beta1 << " vs " << beta2 << endl;

		// alpha_jet1b1 = jet_Pvec.Angle(jet1_Pvec); // in lab frame 
		// alpha_jet1b2 = jet1_Pvec.Angle(jet2_Pvec); // in lab frame 
		//cout << "alpha = " << alpha << endl;
		// OLD: alphaZMF = acos ((jet_PvecZMF.Dot(jet2_PvecZMF)) / (jet_PvecZMF.Mag() * jet2_PvecZMF.Mag())); // in zero mumentum frame (ZMF)
		//cout << "alphaZMF = " << alphaZMF << endl;

                cos_jets = (jet1_Pvec.Dot(jet2_Pvec)) / (jet1_Pvec.Mag() * jet2_Pvec.Mag());

	//--------------inputs for alphaZMF:
		// double mass_1 = 125.0;
		// double Px_1 = jet[0].Px()+jet[1].Px();
		// double Py_1 = jet[0].Py()+jet[1].Py();
		// double Pz_1 = jet[0].Pz()+jet[1].Pz();

		// double mass_2 = 4.70;
		// double Px_2 = jet[0].Px();
		// double Py_2 = jet[0].Py();
		// double Pz_2 = jet[0].Pz();

		// double mass_3 = 4.70;
		// //double mass_3 = jet[1].M(); or Mass();
		// double Px_3 = jet[1].Px();
		// double Py_3 = jet[1].Py();
		// double Pz_3 = jet[1].Pz();

		// alphaZMF = angleinPRF(mass_1, Px_1, Py_1, Pz_1, mass_2, Px_2, Py_2, Pz_2, mass_3, Px_3, Py_3, Pz_3);

                 //======================================= Fill histograms

//if (MET > 20.0) {

                hjet1_pt[r]->Fill(jet1_pt); 
                hjet2_pt[r]->Fill(jet2_pt);

                hjet1_eta[r]->Fill(jet1_eta); 
                hjet2_eta[r]->Fill(jet2_eta); 

                hjet1_P[r]->Fill(jet1_P); 
                hjet2_P[r]->Fill(jet2_P); 

                hjetsInvariantMass[r]->Fill(jetsInvariantMass); 
                hMET[r]->Fill(MET); 
                hHT[r]->Fill(HT); 

                hdeltaeta_jets[r]->Fill(deltaeta_jets); 
                hdeltaphi_jets[r]->Fill(deltaphi_jets); 
                hdeltaR_jets[r]->Fill(deltaR_jets); 
                hcos_jets[r]->Fill(cos_jets); 

                hjetmultiplicity[r]->Fill(jetmultiplicity); 
//}

                  //======================================= final cuts

                //if (MET > 0.0 && HT > 0.0 && jetsInvariantMass > 0.0)
                 {
                    N_accepted++;
                 } 

            } //======================================= event selection 
               
                jet.clear();
  
        } //======================================================= event loop

//============================================================ calculation of sigma and efficiency :


// Decide which group the current sample belongs to
if(r>=1 && r<=4) group=0;   // bkg_hz1-4
if(r==5) group=1;           // bkg_ll
if(r==6) group=2;           // bkg_tautau
if(r==7 || r==8) group=3;   // bkg_ww1-2
if(r==9 || r==10) group=4;  // bkg_zz1-2
if(r>=11 && r<=14) group=5; // bkg_tt1-4
if(r==15) group=6;           // bkg_za
if(r==16 || r==17) group=7; // bkg_jets1-2

N_accepted_group[group] += N_accepted;
N_total_group[group]    += my.fChain->GetEntriesFast();

//cout << "test1" << endl;
} //======================================================= sample loop

////cout << "test2" << endl;


for(int g=0; g<8; g++) {

//cout << "test2 = " << g << endl;
 //if (g != group) continue; 
    if (N_total_group[g] == 0) continue;


    double eff_combined = (double)N_accepted_group[g] / N_total_group[g];
    double sigma_combined = eff_combined * sigma_MG_group[g] * 1000; // fb

    cout << "\n N_total for "<< group_name[g] << " = " << N_total_group[g] <<  endl;
    cout << "N_accepted for "<< group_name[g] << " = " << N_accepted_group[g] <<  endl;
    cout << "efficiency for " << group_name[g] << " = " << eff_combined << endl;
    cout << "sigma = " << sigma_combined << " fb ===" << endl;
}

    // sigma_pre = efficiency_pre * sigma_MG [r] * 1000; //fb 
    // sigma = efficiency * sigma_MG [r] * 1000; //fb 

    // cout << "N_total = " << my.fChain->GetEntriesFast() << endl;

    // cout << "N_accepted_pre = " << N_accepted_pre << endl;
    // cout << "N_accepted = " << N_accepted << endl;

    // cout << "efficiency_pre = " << efficiency_pre << endl;
    // cout << "efficiency = " << efficiency << endl;

    // cout << "sigma_MG = " << sigma_MG [r] << " fb" << endl; //fb 
    // cout << "sigma_pre = " << sigma_pre << " fb" << endl; //fb 
    // cout << "sigma = " << sigma << " fb" << endl; //fb 


//     //======================================= Draw histograms

//     TCanvas *c1 = new TCanvas("c1","p_{T}(jet1)"); c1->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c3 = new TCanvas("c3","p_{T}(jet2)");    c3->Divide(1,1); gStyle->SetOptStat(0);

//     TCanvas *c4 = new TCanvas("c4","#eta(jet1)");    c4->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c6 = new TCanvas("c6","#eta(jet2)");    c6->Divide(1,1); gStyle->SetOptStat(0);

//     TCanvas *c7 = new TCanvas("c7","E(jet1)");    c7->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c9 = new TCanvas("c9","E(jet2)");    c9->Divide(1,1); gStyle->SetOptStat(0);

//     TCanvas *c10 = new TCanvas("c10","IM_{jets}");    c10->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c11 = new TCanvas("c11","MET");    c11->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c12 = new TCanvas("c12","HT");    c12->Divide(1,1); gStyle->SetOptStat(0);

//     TCanvas *c13 = new TCanvas("c13","deltaeta_jets");    c13->Divide(1,1); gStyle->SetOptStat(0);
//     TCanvas *c16 = new TCanvas("c16","deltaphi_jets");    c16->Divide(1,1); gStyle->SetOptStat(0);


//     TCanvas *c19 = new TCanvas("c19","deltaR_jets");    c19->Divide(1,1); gStyle->SetOptStat(0);

//     TCanvas *c22 = new TCanvas("c22","cos_jets");    c22->Divide(1,1); gStyle->SetOptStat(0);


//     TCanvas *c25 = new TCanvas("c25","jetmultiplicity");    c25->Divide(1,1); gStyle->SetOptStat(0);

//     // TCanvas *c26 = new TCanvas("c26","alpha");    c26->Divide(1,1); gStyle->SetOptStat(0);
//     // TCanvas *c27 = new TCanvas("c27","alphaZMF");    c27->Divide(1,1); gStyle->SetOptStat(0);
//     // TCanvas *c28 = new TCanvas("c28","alpha_jet1b1");    c28->Divide(1,1); gStyle->SetOptStat(0);
//     // TCanvas *c29 = new TCanvas("c29","alpha_jet1b2");    c29->Divide(1,1); gStyle->SetOptStat(0);

//     gStyle->SetLineStyleString(11,"100 15");
//     gStyle->SetLineStyleString(12,"150 15");
//     gStyle->SetLineStyleString(13,"200 15");

//     for (int r = 0; r < 18; r++) // other loop over samples
//     {

//         hjet1_pt[r]->Scale(1 /(hjet1_pt[r]->Integral())); 
//         hjet1_pt[r]->SetLineWidth(2);
//         if (r < 2){ hjet1_pt[r]->SetLineColor(600-1-r);  hjet1_pt[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet1_pt[r]->SetLineStyle(2+r);
//             //hjet1_pt[r]->SetLineColor(632-1-r); 
//             hjet1_pt[2]->SetLineColor(2);  
//             hjet1_pt[3]->SetLineColor(3);  
//             hjet1_pt[4]->SetLineColor(7);  
//             hjet1_pt[5]->SetLineColor(6);  
//             hjet1_pt[6]->SetLineColor(kOrange-3);  
//             hjet1_pt[7]->SetLineColor(kPink+2);
//             hjet1_pt[8]->SetLineColor(8); 
//             hjet1_pt[8]->SetLineStyle(1); 
//         }
//         hjet1_pt[r]->GetXaxis()->SetTitle("p_{T}^{light jet} [GeV]");
//         hjet1_pt[r]->GetXaxis()->SetLabelFont(132);
//         hjet1_pt[r]->GetXaxis()->SetTitleFont(22);
//         hjet1_pt[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet1_pt[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet1_pt[r]->GetYaxis()->SetLabelFont(132);
//         hjet1_pt[r]->GetYaxis()->SetTitleFont(22);
//         hjet1_pt[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet1_pt[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hjet_pt[r]->Scale(1 /(hjet_pt[r]->Integral()));
//         hjet_pt[r]->SetLineStyle(2+r);
//         hjet_pt[r]->SetLineWidth(2);
//         if (r < 2){ hjet_pt[r]->SetLineColor(600-1-r);  hjet_pt[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet_pt[r]->SetLineStyle(2+r);
//             //hjet_pt[r]->SetLineColor(632-1-r); 
//             hjet_pt[2]->SetLineColor(2);  
//             hjet_pt[3]->SetLineColor(3);  
//             hjet_pt[4]->SetLineColor(7);  
//             hjet_pt[5]->SetLineColor(6);  
//             hjet_pt[6]->SetLineColor(kOrange-3);  
//             hjet_pt[7]->SetLineColor(kPink+2);  
//             hjet_pt[8]->SetLineColor(8); 
//             hjet_pt[8]->SetLineStyle(1);  
//         }
//         hjet_pt[r]->GetXaxis()->SetTitle("p_{T}^{b_{1}} [GeV]");
//         hjet_pt[r]->GetXaxis()->SetLabelFont(132);
//         hjet_pt[r]->GetXaxis()->SetTitleFont(22);
//         hjet_pt[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet_pt[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet_pt[r]->GetYaxis()->SetLabelFont(132);
//         hjet_pt[r]->GetYaxis()->SetTitleFont(22);
//         hjet_pt[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet_pt[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hjet2_pt[r]->Scale(1 /(hjet2_pt[r]->Integral()));
//         hjet2_pt[r]->SetLineStyle(2+r);
//         hjet2_pt[r]->SetLineWidth(2);
//         if (r < 2){ hjet2_pt[r]->SetLineColor(600-1-r);  hjet2_pt[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet2_pt[r]->SetLineStyle(2+r);
//             //hjet2_pt[r]->SetLineColor(632-1-r); 
//             hjet2_pt[2]->SetLineColor(2);  
//             hjet2_pt[3]->SetLineColor(3);  
//             hjet2_pt[4]->SetLineColor(7);  
//             hjet2_pt[5]->SetLineColor(6);  
//             hjet2_pt[6]->SetLineColor(kOrange-3);  
//             hjet2_pt[7]->SetLineColor(kPink+2);  
//             hjet2_pt[8]->SetLineColor(8);  
//             hjet2_pt[8]->SetLineStyle(1); 
//         }
//         hjet2_pt[r]->GetXaxis()->SetTitle("p_{T}^{b_{2}} [GeV]");
//         hjet2_pt[r]->GetXaxis()->SetLabelFont(132);
//         hjet2_pt[r]->GetXaxis()->SetTitleFont(22);
//         hjet2_pt[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet2_pt[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet2_pt[r]->GetYaxis()->SetLabelFont(132);
//         hjet2_pt[r]->GetYaxis()->SetTitleFont(22);
//         hjet2_pt[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet2_pt[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hjet1_eta[r]->Scale(1 /(hjet1_eta[r]->Integral()));
//         hjet1_eta[r]->SetLineStyle(2+r);
//         hjet1_eta[r]->SetLineWidth(2);
//         if (r < 2){ hjet1_eta[r]->SetLineColor(600-1-r);  hjet1_eta[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet1_eta[r]->SetLineStyle(2+r);
//             //hjet1_eta[r]->SetLineColor(632-1-r); 
//             hjet1_eta[2]->SetLineColor(2);  
//             hjet1_eta[3]->SetLineColor(3);  
//             hjet1_eta[4]->SetLineColor(7);  
//             hjet1_eta[5]->SetLineColor(6);  
//             hjet1_eta[6]->SetLineColor(kOrange-3);  
//             hjet1_eta[7]->SetLineColor(kPink+2);  
//             hjet1_eta[8]->SetLineColor(8);  
//             hjet1_eta[8]->SetLineStyle(1); 
//         }
//         hjet1_eta[r]->GetXaxis()->SetTitle("#eta_{light jet}");
//         hjet1_eta[r]->GetXaxis()->SetLabelFont(132);
//         hjet1_eta[r]->GetXaxis()->SetTitleFont(22);
//         hjet1_eta[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet1_eta[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet1_eta[r]->GetYaxis()->SetLabelFont(132);
//         hjet1_eta[r]->GetYaxis()->SetTitleFont(22);
//         hjet1_eta[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet1_eta[r]->GetYaxis()->SetRangeUser(0.0, .3);

//         hjet_eta[r]->Scale(1 /(hjet_eta[r]->Integral()));
//         hjet_eta[r]->SetLineStyle(2+r);
//         hjet_eta[r]->SetLineWidth(2);
//         if (r < 2){ hjet_eta[r]->SetLineColor(600-1-r);  hjet_eta[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet_eta[r]->SetLineStyle(2+r);
//             //hjet_eta[r]->SetLineColor(632-1-r); 
//             hjet_eta[2]->SetLineColor(2);  
//             hjet_eta[3]->SetLineColor(3);  
//             hjet_eta[4]->SetLineColor(7);  
//             hjet_eta[5]->SetLineColor(6);  
//             hjet_eta[6]->SetLineColor(kOrange-3);  
//             hjet_eta[7]->SetLineColor(kPink+2);  
//             hjet_eta[8]->SetLineColor(8);  
//             hjet_eta[8]->SetLineStyle(1); 
//         }
//         hjet_eta[r]->GetXaxis()->SetTitle("#eta_{b_{1}}");
//         hjet_eta[r]->GetXaxis()->SetLabelFont(132);
//         hjet_eta[r]->GetXaxis()->SetTitleFont(22);
//         hjet_eta[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet_eta[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet_eta[r]->GetYaxis()->SetLabelFont(132);
//         hjet_eta[r]->GetYaxis()->SetTitleFont(22);
//         hjet_eta[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet_eta[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hjet2_eta[r]->Scale(1 /(hjet2_eta[r]->Integral()));
//         hjet2_eta[r]->SetLineStyle(2+r);
//         hjet2_eta[r]->SetLineWidth(2);
//         if (r < 2){ hjet2_eta[r]->SetLineColor(600-1-r);  hjet2_eta[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet2_eta[r]->SetLineStyle(2+r);
//             //hjet2_eta[r]->SetLineColor(632-1-r); 
//             hjet2_eta[2]->SetLineColor(2);  
//             hjet2_eta[3]->SetLineColor(3);  
//             hjet2_eta[4]->SetLineColor(7);  
//             hjet2_eta[5]->SetLineColor(6);  
//             hjet2_eta[6]->SetLineColor(kOrange-3);  
//             hjet2_eta[7]->SetLineColor(kPink+2);  
//             hjet2_eta[8]->SetLineColor(8);  
//             hjet2_eta[8]->SetLineStyle(1); 
//         }
//         hjet2_eta[r]->GetXaxis()->SetTitle("#eta_{b_{2}}");
//         hjet2_eta[r]->GetXaxis()->SetLabelFont(132);
//         hjet2_eta[r]->GetXaxis()->SetTitleFont(22);
//         hjet2_eta[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet2_eta[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet2_eta[r]->GetYaxis()->SetLabelFont(132);
//         hjet2_eta[r]->GetYaxis()->SetTitleFont(22);
//         hjet2_eta[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet2_eta[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hjet1_P[r]->Scale(1 /(hjet1_P[r]->Integral()));
//         hjet1_P[r]->SetLineStyle(2+r);
//         hjet1_P[r]->SetLineWidth(2);
//         if (r < 2){ hjet1_P[r]->SetLineColor(600-1-r);  hjet1_P[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet1_P[r]->SetLineStyle(2+r);
//             //hjet1_P[r]->SetLineColor(632-1-r); 
//             hjet1_P[2]->SetLineColor(2);  
//             hjet1_P[3]->SetLineColor(3);  
//             hjet1_P[4]->SetLineColor(7);  
//             hjet1_P[5]->SetLineColor(6);  
//             hjet1_P[6]->SetLineColor(kOrange-3);  
//             hjet1_P[7]->SetLineColor(kPink+2);  
//             hjet1_P[8]->SetLineColor(8);  
//             hjet1_P[8]->SetLineStyle(1); 
//         }
//         hjet1_P[r]->GetXaxis()->SetTitle("P_{light jet} [GeV]");
//         hjet1_P[r]->GetXaxis()->SetLabelFont(132);
//         hjet1_P[r]->GetXaxis()->SetTitleFont(22);
//         hjet1_P[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet1_P[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet1_P[r]->GetYaxis()->SetLabelFont(132);
//         hjet1_P[r]->GetYaxis()->SetTitleFont(22);
//         hjet1_P[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet1_P[r]->GetYaxis()->SetRangeUser(0.0, .05);

//         hjet_P[r]->Scale(1 /(hjet_P[r]->Integral()));
//         hjet_P[r]->SetLineStyle(2+r);
//         hjet_P[r]->SetLineWidth(2);
//         if (r < 2){ hjet_P[r]->SetLineColor(600-1-r);  hjet_P[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet_P[r]->SetLineStyle(2+r);
//             //hjet_P[r]->SetLineColor(632-1-r); 
//             hjet_P[2]->SetLineColor(2);  
//             hjet_P[3]->SetLineColor(3);  
//             hjet_P[4]->SetLineColor(7);  
//             hjet_P[5]->SetLineColor(6);  
//             hjet_P[6]->SetLineColor(kOrange-3);  
//             hjet_P[7]->SetLineColor(kPink+2);  
//             hjet_P[8]->SetLineColor(8);  
//             hjet_P[8]->SetLineStyle(1); 
//         }
//         hjet_P[r]->GetXaxis()->SetTitle("P_{b_{1}} [GeV]");
//         hjet_P[r]->GetXaxis()->SetLabelFont(132);
//         hjet_P[r]->GetXaxis()->SetTitleFont(22);
//         hjet_P[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet_P[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet_P[r]->GetYaxis()->SetLabelFont(132);
//         hjet_P[r]->GetYaxis()->SetTitleFont(22);
//         hjet_P[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet_P[r]->GetYaxis()->SetRangeUser(0.0, .07);

//         hjet2_P[r]->Scale(1 /(hjet2_P[r]->Integral()));
//         hjet2_P[r]->SetLineStyle(2+r);
//         hjet2_P[r]->SetLineWidth(2);
//         if (r < 2){ hjet2_P[r]->SetLineColor(600-1-r);  hjet2_P[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjet2_P[r]->SetLineStyle(2+r);
//             //hjet2_P[r]->SetLineColor(632-1-r); 
//             hjet2_P[2]->SetLineColor(2);  
//             hjet2_P[3]->SetLineColor(3);  
//             hjet2_P[4]->SetLineColor(7);  
//             hjet2_P[5]->SetLineColor(6);  
//             hjet2_P[6]->SetLineColor(kOrange-3);  
//             hjet2_P[7]->SetLineColor(kPink+2);  
//             hjet2_P[8]->SetLineColor(8);  
//             hjet2_P[8]->SetLineStyle(1); 
//         }
//         hjet2_P[r]->GetXaxis()->SetTitle("P_{b_{2}} [GeV]");
//         hjet2_P[r]->GetXaxis()->SetLabelFont(132);
//         hjet2_P[r]->GetXaxis()->SetTitleFont(22);
//         hjet2_P[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjet2_P[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjet2_P[r]->GetYaxis()->SetLabelFont(132);
//         hjet2_P[r]->GetYaxis()->SetTitleFont(22);
//         hjet2_P[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjet2_P[r]->GetYaxis()->SetRangeUser(0.0, .12);

//         hjetsInvariantMass[r]->Scale(1 /(hjetsInvariantMass[r]->Integral()));
//         hjetsInvariantMass[r]->SetLineStyle(2+r);
//         hjetsInvariantMass[r]->SetLineWidth(2);
//         if (r < 2){ hjetsInvariantMass[r]->SetLineColor(600-1-r);  hjetsInvariantMass[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjetsInvariantMass[r]->SetLineStyle(2+r);
//             //hjetsInvariantMass[r]->SetLineColor(632-1-r); 
//             hjetsInvariantMass[2]->SetLineColor(2);  
//             hjetsInvariantMass[3]->SetLineColor(3);  
//             hjetsInvariantMass[4]->SetLineColor(7);  
//             hjetsInvariantMass[5]->SetLineColor(6);  
//             hjetsInvariantMass[6]->SetLineColor(kOrange-3);  
//             hjetsInvariantMass[7]->SetLineColor(kPink+2);  
//             hjetsInvariantMass[8]->SetLineColor(8);  
//             hjetsInvariantMass[8]->SetLineStyle(1); 
//         }
//         hjetsInvariantMass[r]->GetXaxis()->SetTitle("M_{b_{1},b_{2}} [GeV]");
//         hjetsInvariantMass[r]->GetXaxis()->SetLabelFont(132);
//         hjetsInvariantMass[r]->GetXaxis()->SetTitleFont(22);
//         hjetsInvariantMass[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjetsInvariantMass[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjetsInvariantMass[r]->GetYaxis()->SetLabelFont(132);
//         hjetsInvariantMass[r]->GetYaxis()->SetTitleFont(22);
//         hjetsInvariantMass[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjetsInvariantMass[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hMET[r]->Scale(1 /(hMET[r]->Integral()));
//         hMET[r]->SetLineStyle(2+r);
//         hMET[r]->SetLineWidth(2);
//         if (r < 2){ hMET[r]->SetLineColor(600-1-r);  hMET[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hMET[r]->SetLineStyle(2+r);
//             //hMET[r]->SetLineColor(632-1-r); 
//             hMET[2]->SetLineColor(2);  
//             hMET[3]->SetLineColor(3);  
//             hMET[4]->SetLineColor(7);  
//             hMET[5]->SetLineColor(6);  
//             hMET[6]->SetLineColor(kOrange-3);  
//             hMET[7]->SetLineColor(kPink+2);  
//             hMET[8]->SetLineColor(8);  
//             hMET[8]->SetLineStyle(1); 
//         }
//         hMET[r]->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
//         hMET[r]->GetXaxis()->SetLabelFont(132);
//         hMET[r]->GetXaxis()->SetTitleFont(22);
//         hMET[r]->GetXaxis()->SetTitleOffset(1.10);
//         hMET[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hMET[r]->GetYaxis()->SetLabelFont(132);
//         hMET[r]->GetYaxis()->SetTitleFont(22);
//         hMET[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hMET[r]->GetYaxis()->SetRangeUser(0.0, .1);

//         hHT[r]->Scale(1 /(hHT[r]->Integral()));
//         hHT[r]->SetLineStyle(2+r);
//         hHT[r]->SetLineWidth(2);
//         if (r < 2){ hHT[r]->SetLineColor(600-1-r);  hHT[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hHT[r]->SetLineStyle(2+r);
//             //hHT[r]->SetLineColor(632-1-r); 
//             hHT[2]->SetLineColor(2);  
//             hHT[3]->SetLineColor(3);  
//             hHT[4]->SetLineColor(7);  
//             hHT[5]->SetLineColor(6);  
//             hHT[6]->SetLineColor(kOrange-3);  
//             hHT[7]->SetLineColor(kPink+2);  
//             hHT[8]->SetLineColor(8);  
//             hHT[8]->SetLineStyle(1); 
//         }
//         hHT[r]->GetXaxis()->SetTitle("H_{T} [GeV]");
//         hHT[r]->GetXaxis()->SetLabelFont(132);
//         hHT[r]->GetXaxis()->SetTitleFont(22);
//         hHT[r]->GetXaxis()->SetTitleOffset(1.10);
//         hHT[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hHT[r]->GetYaxis()->SetLabelFont(132);
//         hHT[r]->GetYaxis()->SetTitleFont(22);
//         hHT[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hHT[r]->GetYaxis()->SetRangeUser(0.0, .07);

//         hdeltaeta_jets[r]->Scale(1 /(hdeltaeta_jets[r]->Integral()));
//         hdeltaeta_jets[r]->SetLineStyle(2+r);
//         hdeltaeta_jets[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaeta_jets[r]->SetLineColor(600-1-r);  hdeltaeta_jets[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaeta_jets[r]->SetLineStyle(2+r);
//             //hdeltaeta_jets[r]->SetLineColor(632-1-r); 
//             hdeltaeta_jets[2]->SetLineColor(2);  
//             hdeltaeta_jets[3]->SetLineColor(3);  
//             hdeltaeta_jets[4]->SetLineColor(7);  
//             hdeltaeta_jets[5]->SetLineColor(6);  
//             hdeltaeta_jets[6]->SetLineColor(kOrange-3);  
//             hdeltaeta_jets[7]->SetLineColor(kPink+2);  
//             hdeltaeta_jets[8]->SetLineColor(8); 
//             hdeltaeta_jets[8]->SetLineStyle(1);  
//         }
//         hdeltaeta_jets[r]->GetXaxis()->SetTitle("#Delta#eta(b_{1},b_{2})");
//         hdeltaeta_jets[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaeta_jets[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaeta_jets[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaeta_jets[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaeta_jets[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaeta_jets[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaeta_jets[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaeta_jets[r]->GetYaxis()->SetRangeUser(0.0, .18);

//         hdeltaeta_jet1jet[r]->Scale(1 /(hdeltaeta_jet1jet[r]->Integral()));
//         hdeltaeta_jet1jet[r]->SetLineStyle(2+r);
//         hdeltaeta_jet1jet[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaeta_jet1jet[r]->SetLineColor(600-1-r);  hdeltaeta_jet1jet[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaeta_jet1jet[r]->SetLineStyle(2+r);
//             //hdeltaeta_jet1jet[r]->SetLineColor(632-1-r); 
//             hdeltaeta_jet1jet[2]->SetLineColor(2);  
//             hdeltaeta_jet1jet[3]->SetLineColor(3);  
//             hdeltaeta_jet1jet[4]->SetLineColor(7);  
//             hdeltaeta_jet1jet[5]->SetLineColor(6);  
//             hdeltaeta_jet1jet[6]->SetLineColor(kOrange-3);  
//             hdeltaeta_jet1jet[7]->SetLineColor(kPink+2);  
//             hdeltaeta_jet1jet[8]->SetLineColor(8);  
//             hdeltaeta_jet1jet[8]->SetLineStyle(1); 
//         }
//         hdeltaeta_jet1jet[r]->GetXaxis()->SetTitle("#Delta#eta(b_{1},l-jet)");
//         hdeltaeta_jet1jet[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaeta_jet1jet[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaeta_jet1jet[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaeta_jet1jet[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaeta_jet1jet[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaeta_jet1jet[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaeta_jet1jet[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaeta_jet1jet[r]->GetYaxis()->SetRangeUser(0.0, .1);

//         hdeltaeta_jet1jet2[r]->Scale(1 /(hdeltaeta_jet1jet2[r]->Integral()));
//         hdeltaeta_jet1jet2[r]->SetLineStyle(2+r);
//         hdeltaeta_jet1jet2[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaeta_jet1jet2[r]->SetLineColor(600-1-r);  hdeltaeta_jet1jet2[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaeta_jet1jet2[r]->SetLineStyle(2+r);
//             //hdeltaeta_jet1jet2[r]->SetLineColor(632-1-r); 
//             hdeltaeta_jet1jet2[2]->SetLineColor(2);  
//             hdeltaeta_jet1jet2[3]->SetLineColor(3);  
//             hdeltaeta_jet1jet2[4]->SetLineColor(7);  
//             hdeltaeta_jet1jet2[5]->SetLineColor(6);  
//             hdeltaeta_jet1jet2[6]->SetLineColor(kOrange-3);  
//             hdeltaeta_jet1jet2[7]->SetLineColor(kPink+2);  
//             hdeltaeta_jet1jet2[8]->SetLineColor(8);  
//             hdeltaeta_jet1jet2[8]->SetLineStyle(1); 
//         }
//         hdeltaeta_jet1jet2[r]->GetXaxis()->SetTitle("#Delta#eta(l-jet,b_{2})");
//         hdeltaeta_jet1jet2[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaeta_jet1jet2[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaeta_jet1jet2[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaeta_jet1jet2[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaeta_jet1jet2[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaeta_jet1jet2[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaeta_jet1jet2[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaeta_jet1jet2[r]->GetYaxis()->SetRangeUser(0.0, .1);

//         hdeltaphi_jets[r]->Scale(1 /(hdeltaphi_jets[r]->Integral()));
//         hdeltaphi_jets[r]->SetLineStyle(2+r);
//         hdeltaphi_jets[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaphi_jets[r]->SetLineColor(600-1-r);  hdeltaphi_jets[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaphi_jets[r]->SetLineStyle(2+r);
//             //hdeltaphi_jets[r]->SetLineColor(632-1-r); 
//             hdeltaphi_jets[2]->SetLineColor(2);  
//             hdeltaphi_jets[3]->SetLineColor(3);  
//             hdeltaphi_jets[4]->SetLineColor(7);  
//             hdeltaphi_jets[5]->SetLineColor(6);  
//             hdeltaphi_jets[6]->SetLineColor(kOrange-3);  
//             hdeltaphi_jets[7]->SetLineColor(kPink+2);  
//             hdeltaphi_jets[8]->SetLineColor(8);  
//             hdeltaphi_jets[8]->SetLineStyle(1); 
//         }
//         hdeltaphi_jets[r]->GetXaxis()->SetTitle("#Delta#phi(b_{1},b_{2})");
//         hdeltaphi_jets[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaphi_jets[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaphi_jets[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaphi_jets[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaphi_jets[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaphi_jets[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaphi_jets[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaphi_jets[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hdeltaphi_jet1jet[r]->Scale(1 /(hdeltaphi_jet1jet[r]->Integral()));
//         hdeltaphi_jet1jet[r]->SetLineStyle(2+r);
//         hdeltaphi_jet1jet[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaphi_jet1jet[r]->SetLineColor(600-1-r);  hdeltaphi_jet1jet[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaphi_jet1jet[r]->SetLineStyle(2+r);
//             //hdeltaphi_jet1jet[r]->SetLineColor(632-1-r); 
//             hdeltaphi_jet1jet[2]->SetLineColor(2);  
//             hdeltaphi_jet1jet[3]->SetLineColor(3);  
//             hdeltaphi_jet1jet[4]->SetLineColor(7);  
//             hdeltaphi_jet1jet[5]->SetLineColor(6);  
//             hdeltaphi_jet1jet[6]->SetLineColor(kOrange-3);  
//             hdeltaphi_jet1jet[7]->SetLineColor(kPink+2);  
//             hdeltaphi_jet1jet[8]->SetLineColor(8);  
//             hdeltaphi_jet1jet[8]->SetLineStyle(1); 
//         }
//         hdeltaphi_jet1jet[r]->GetXaxis()->SetTitle("#Delta#phi(b_{1},l-jet)");
//         hdeltaphi_jet1jet[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaphi_jet1jet[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaphi_jet1jet[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaphi_jet1jet[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaphi_jet1jet[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaphi_jet1jet[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaphi_jet1jet[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaphi_jet1jet[r]->GetYaxis()->SetRangeUser(0.0, .1);

//         hdeltaphi_jet1jet2[r]->Scale(1 /(hdeltaphi_jet1jet2[r]->Integral()));
//         hdeltaphi_jet1jet2[r]->SetLineStyle(2+r);
//         hdeltaphi_jet1jet2[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaphi_jet1jet2[r]->SetLineColor(600-1-r);  hdeltaphi_jet1jet2[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaphi_jet1jet2[r]->SetLineStyle(2+r);
//             //hdeltaphi_jet1jet2[r]->SetLineColor(632-1-r); 
//             hdeltaphi_jet1jet2[2]->SetLineColor(2);  
//             hdeltaphi_jet1jet2[3]->SetLineColor(3);  
//             hdeltaphi_jet1jet2[4]->SetLineColor(7);  
//             hdeltaphi_jet1jet2[5]->SetLineColor(6);  
//             hdeltaphi_jet1jet2[6]->SetLineColor(kOrange-3);  
//             hdeltaphi_jet1jet2[7]->SetLineColor(kPink+2);  
//             hdeltaphi_jet1jet2[8]->SetLineColor(8);  
//             hdeltaphi_jet1jet2[8]->SetLineStyle(1); 
//         }
//         hdeltaphi_jet1jet2[r]->GetXaxis()->SetTitle("#Delta#phi(l-jet,b_{2})");
//         hdeltaphi_jet1jet2[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaphi_jet1jet2[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaphi_jet1jet2[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaphi_jet1jet2[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaphi_jet1jet2[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaphi_jet1jet2[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaphi_jet1jet2[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaphi_jet1jet2[r]->GetYaxis()->SetRangeUser(0.0, .1);

//         hdeltaR_jets[r]->Scale(1 /(hdeltaR_jets[r]->Integral()));
//         hdeltaR_jets[r]->SetLineStyle(2+r);
//         hdeltaR_jets[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaR_jets[r]->SetLineColor(600-1-r);  hdeltaR_jets[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaR_jets[r]->SetLineStyle(2+r);
//             //hdeltaR_jets[r]->SetLineColor(632-1-r); 
//             hdeltaR_jets[2]->SetLineColor(2);  
//             hdeltaR_jets[3]->SetLineColor(3);  
//             hdeltaR_jets[4]->SetLineColor(7);  
//             hdeltaR_jets[5]->SetLineColor(6);  
//             hdeltaR_jets[6]->SetLineColor(kOrange-3);  
//             hdeltaR_jets[7]->SetLineColor(kPink+2);  
//             hdeltaR_jets[8]->SetLineColor(8);  
//             hdeltaR_jets[8]->SetLineStyle(1); 
//         }
//         hdeltaR_jets[r]->GetXaxis()->SetTitle("#DeltaR(b_{1},b_{2})");
//         hdeltaR_jets[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaR_jets[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaR_jets[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaR_jets[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaR_jets[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaR_jets[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaR_jets[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaR_jets[r]->GetYaxis()->SetRangeUser(0.0, .3);

//         hdeltaR_jet1jet[r]->Scale(1 /(hdeltaR_jet1jet[r]->Integral()));
//         hdeltaR_jet1jet[r]->SetLineStyle(2+r);
//         hdeltaR_jet1jet[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaR_jet1jet[r]->SetLineColor(600-1-r);  hdeltaR_jet1jet[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaR_jet1jet[r]->SetLineStyle(2+r);
//             //hdeltaR_jet1jet[r]->SetLineColor(632-1-r); 
//             hdeltaR_jet1jet[2]->SetLineColor(2);  
//             hdeltaR_jet1jet[3]->SetLineColor(3);  
//             hdeltaR_jet1jet[4]->SetLineColor(7);  
//             hdeltaR_jet1jet[5]->SetLineColor(6);  
//             hdeltaR_jet1jet[6]->SetLineColor(kOrange-3);  
//             hdeltaR_jet1jet[7]->SetLineColor(kPink+2);  
//             hdeltaR_jet1jet[8]->SetLineColor(8);  
//             hdeltaR_jet1jet[8]->SetLineStyle(1); 
//         }
//         hdeltaR_jet1jet[r]->GetXaxis()->SetTitle("#DeltaR(l-jet,b_{1})");
//         hdeltaR_jet1jet[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaR_jet1jet[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaR_jet1jet[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaR_jet1jet[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaR_jet1jet[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaR_jet1jet[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaR_jet1jet[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaR_jet1jet[r]->GetYaxis()->SetRangeUser(0.0, .15);

//         hdeltaR_jet1jet2[r]->Scale(1 /(hdeltaR_jet1jet2[r]->Integral()));
//         hdeltaR_jet1jet2[r]->SetLineStyle(2+r);
//         hdeltaR_jet1jet2[r]->SetLineWidth(2);
//         if (r < 2){ hdeltaR_jet1jet2[r]->SetLineColor(600-1-r);  hdeltaR_jet1jet2[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hdeltaR_jet1jet2[r]->SetLineStyle(2+r);
//             //hdeltaR_jet1jet2[r]->SetLineColor(632-1-r); 
//             hdeltaR_jet1jet2[2]->SetLineColor(2);  
//             hdeltaR_jet1jet2[3]->SetLineColor(3);  
//             hdeltaR_jet1jet2[4]->SetLineColor(7);  
//             hdeltaR_jet1jet2[5]->SetLineColor(6);  
//             hdeltaR_jet1jet2[6]->SetLineColor(kOrange-3);  
//             hdeltaR_jet1jet2[7]->SetLineColor(kPink+2);  
//             hdeltaR_jet1jet2[8]->SetLineColor(8);  
//             hdeltaR_jet1jet2[8]->SetLineStyle(1); 
//         }
//         hdeltaR_jet1jet2[r]->GetXaxis()->SetTitle("#DeltaR(l-jet,b_{2})");
//         hdeltaR_jet1jet2[r]->GetXaxis()->SetLabelFont(132);
//         hdeltaR_jet1jet2[r]->GetXaxis()->SetTitleFont(22);
//         hdeltaR_jet1jet2[r]->GetXaxis()->SetTitleOffset(1.10);
//         hdeltaR_jet1jet2[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hdeltaR_jet1jet2[r]->GetYaxis()->SetLabelFont(132);
//         hdeltaR_jet1jet2[r]->GetYaxis()->SetTitleFont(22);
//         hdeltaR_jet1jet2[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hdeltaR_jet1jet2[r]->GetYaxis()->SetRangeUser(0.0, .15);

//         hcos_jets[r]->Scale(1 /(hcos_jets[r]->Integral()));
//         hcos_jets[r]->SetLineStyle(2+r);
//         hcos_jets[r]->SetLineWidth(2);
//         if (r < 2){ hcos_jets[r]->SetLineColor(600-1-r);  hcos_jets[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hcos_jets[r]->SetLineStyle(2+r);
//             //hcos_jets[r]->SetLineColor(632-1-r); 
//             hcos_jets[2]->SetLineColor(2);  
//             hcos_jets[3]->SetLineColor(3);  
//             hcos_jets[4]->SetLineColor(7);  
//             hcos_jets[5]->SetLineColor(6);  
//             hcos_jets[6]->SetLineColor(kOrange-3);  
//             hcos_jets[7]->SetLineColor(kPink+2);  
//             hcos_jets[8]->SetLineColor(8);  
//             hcos_jets[8]->SetLineStyle(1); 
//         }
//         hcos_jets[r]->GetXaxis()->SetTitle("cos (b_{1},b_{2})");
//         hcos_jets[r]->GetXaxis()->SetLabelFont(132);
//         hcos_jets[r]->GetXaxis()->SetTitleFont(22);
//         hcos_jets[r]->GetXaxis()->SetTitleOffset(1.10);
//         hcos_jets[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hcos_jets[r]->GetYaxis()->SetLabelFont(132);
//         hcos_jets[r]->GetYaxis()->SetTitleFont(22);
//         hcos_jets[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hcos_jets[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         hcos_jet1jet[r]->Scale(1 /(hcos_jet1jet[r]->Integral()));
//         hcos_jet1jet[r]->SetLineStyle(2+r);
//         hcos_jet1jet[r]->SetLineWidth(2);
//         if (r < 2){ hcos_jet1jet[r]->SetLineColor(600-1-r);  hcos_jet1jet[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hcos_jet1jet[r]->SetLineStyle(2+r);
//             //hcos_jet1jet[r]->SetLineColor(632-1-r); 
//             hcos_jet1jet[2]->SetLineColor(2);  
//             hcos_jet1jet[3]->SetLineColor(3);  
//             hcos_jet1jet[4]->SetLineColor(7);  
//             hcos_jet1jet[5]->SetLineColor(6);  
//             hcos_jet1jet[6]->SetLineColor(kOrange-3);  
//             hcos_jet1jet[7]->SetLineColor(kPink+2);  
//             hcos_jet1jet[8]->SetLineColor(8);  
//             hcos_jet1jet[8]->SetLineStyle(1); 
//         }
//         hcos_jet1jet[r]->GetXaxis()->SetTitle("cos (l-jet,b_{1})");
//         hcos_jet1jet[r]->GetXaxis()->SetLabelFont(132);
//         hcos_jet1jet[r]->GetXaxis()->SetTitleFont(22);
//         hcos_jet1jet[r]->GetXaxis()->SetTitleOffset(1.10);
//         hcos_jet1jet[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hcos_jet1jet[r]->GetYaxis()->SetLabelFont(132);
//         hcos_jet1jet[r]->GetYaxis()->SetTitleFont(22);
//         hcos_jet1jet[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hcos_jet1jet[r]->GetYaxis()->SetRangeUser(0.0, .45);

//         hcos_jet1jet2[r]->Scale(1 /(hcos_jet1jet2[r]->Integral()));
//         hcos_jet1jet2[r]->SetLineStyle(2+r);
//         hcos_jet1jet2[r]->SetLineWidth(2);
//         if (r < 2){ hcos_jet1jet2[r]->SetLineColor(600-1-r);  hcos_jet1jet2[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hcos_jet1jet2[r]->SetLineStyle(2+r);
//             //hcos_jet1jet2[r]->SetLineColor(632-1-r); 
//             hcos_jet1jet2[2]->SetLineColor(2);  
//             hcos_jet1jet2[3]->SetLineColor(3);  
//             hcos_jet1jet2[4]->SetLineColor(7);  
//             hcos_jet1jet2[5]->SetLineColor(6);  
//             hcos_jet1jet2[6]->SetLineColor(kOrange-3);  
//             hcos_jet1jet2[7]->SetLineColor(kPink+2);  
//             hcos_jet1jet2[8]->SetLineColor(8);  
//             hcos_jet1jet2[8]->SetLineStyle(1); 
//         }
//         hcos_jet1jet2[r]->GetXaxis()->SetTitle("cos (l-jet,b_{2})");
//         hcos_jet1jet2[r]->GetXaxis()->SetLabelFont(132);
//         hcos_jet1jet2[r]->GetXaxis()->SetTitleFont(22);
//         hcos_jet1jet2[r]->GetXaxis()->SetTitleOffset(1.10);
//         hcos_jet1jet2[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hcos_jet1jet2[r]->GetYaxis()->SetLabelFont(132);
//         hcos_jet1jet2[r]->GetYaxis()->SetTitleFont(22);
//         hcos_jet1jet2[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hcos_jet1jet2[r]->GetYaxis()->SetRangeUser(0.0, .4);

//         hjetmultiplicity[r]->Scale(1 /(hjetmultiplicity[r]->Integral()));
//         hjetmultiplicity[r]->SetLineStyle(2+r);
//         hjetmultiplicity[r]->SetLineWidth(2);
//         if (r < 2){ hjetmultiplicity[r]->SetLineColor(600-1-r);  hjetmultiplicity[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             hjetmultiplicity[r]->SetLineStyle(2+r);
//             //hjetmultiplicity[r]->SetLineColor(632-1-r); 
//             hjetmultiplicity[2]->SetLineColor(2);  
//             hjetmultiplicity[3]->SetLineColor(3);  
//             hjetmultiplicity[4]->SetLineColor(7);  
//             hjetmultiplicity[5]->SetLineColor(6);  
//             hjetmultiplicity[6]->SetLineColor(kOrange-3);  
//             hjetmultiplicity[7]->SetLineColor(kPink+2);  
//             hjetmultiplicity[8]->SetLineColor(8);  
//             hjetmultiplicity[8]->SetLineStyle(1); 
//         }
//         hjetmultiplicity[r]->GetXaxis()->SetTitle("Jet multiplicity");
//         hjetmultiplicity[r]->GetXaxis()->SetLabelFont(132);
//         hjetmultiplicity[r]->GetXaxis()->SetTitleFont(22);
//         hjetmultiplicity[r]->GetXaxis()->SetTitleOffset(1.10);
//         hjetmultiplicity[r]->GetYaxis()->SetTitle("Normalized distribution");
//         hjetmultiplicity[r]->GetYaxis()->SetLabelFont(132);
//         hjetmultiplicity[r]->GetYaxis()->SetTitleFont(22);
//         hjetmultiplicity[r]->GetYaxis()->SetTitleOffset(1.10);
//         //hjetmultiplicity[r]->GetYaxis()->SetRangeUser(0.0, .2);

//         halpha[r]->Scale(1 /(halpha[r]->Integral()));
//         halpha[r]->SetLineStyle(2+r);
//         halpha[r]->SetLineWidth(2);
//         if (r < 2){ halpha[r]->SetLineColor(600-1-r);  halpha[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             halpha[r]->SetLineStyle(2+r);
//             //halpha[r]->SetLineColor(632-1-r); 
//             halpha[2]->SetLineColor(2);  
//             halpha[3]->SetLineColor(3);  
//             halpha[4]->SetLineColor(7);  
//             halpha[5]->SetLineColor(6);  
//             halpha[6]->SetLineColor(kOrange-3);  
//             halpha[7]->SetLineColor(kPink+2);  
//             halpha[8]->SetLineColor(8);  
//             halpha[8]->SetLineStyle(1); 
//         }
//         halpha[r]->GetXaxis()->SetTitle("#alpha(b_{1},b_{2}) [Rad]");
//         halpha[r]->GetXaxis()->SetLabelFont(132);
//         halpha[r]->GetXaxis()->SetTitleFont(22);
//         halpha[r]->GetXaxis()->SetTitleOffset(1.10);
//         halpha[r]->GetYaxis()->SetTitle("Normalized distribution");
//         halpha[r]->GetYaxis()->SetLabelFont(132);
//         halpha[r]->GetYaxis()->SetTitleFont(22);
//         halpha[r]->GetYaxis()->SetTitleOffset(1.10);
//         //halpha[r]->GetYaxis()->SetRangeUser(0.0, .2);
// //-----------------
//         halphaZMF[r]->Scale(1 /(halphaZMF[r]->Integral()));
//         halphaZMF[r]->SetLineStyle(2+r);
//         halphaZMF[r]->SetLineWidth(2);
//         if (r < 2){ halphaZMF[r]->SetLineColor(600-1-r);  halphaZMF[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             halphaZMF[r]->SetLineStyle(2+r);
//             //halphaZMF[r]->SetLineColor(632-1-r); 
//             halphaZMF[2]->SetLineColor(2);  
//             halphaZMF[3]->SetLineColor(3);  
//             halphaZMF[4]->SetLineColor(7);  
//             halphaZMF[5]->SetLineColor(6);  
//             halphaZMF[6]->SetLineColor(kOrange-3);  
//             halphaZMF[7]->SetLineColor(kPink+2);  
//             halphaZMF[8]->SetLineColor(8);  
//             halphaZMF[8]->SetLineStyle(1); 
//         }
//         halphaZMF[r]->GetXaxis()->SetTitle("#alpha_{ZMF}(b_{1},b_{2}) [Rad]");
//         halphaZMF[r]->GetXaxis()->SetLabelFont(132);
//         halphaZMF[r]->GetXaxis()->SetTitleFont(22);
//         halphaZMF[r]->GetXaxis()->SetTitleOffset(1.10);
//         halphaZMF[r]->GetYaxis()->SetTitle("Normalized distribution");
//         halphaZMF[r]->GetYaxis()->SetLabelFont(132);
//         halphaZMF[r]->GetYaxis()->SetTitleFont(22);
//         halphaZMF[r]->GetYaxis()->SetTitleOffset(1.10);
//         //halphaZMF[r]->GetYaxis()->SetRangeUser(0.0, .2);
// //-----------------
//         halpha_jet1b1[r]->Scale(1 /(halpha_jet1b1[r]->Integral()));
//         halpha_jet1b1[r]->SetLineStyle(2+r);
//         halpha_jet1b1[r]->SetLineWidth(2);
//         if (r < 2){ halpha_jet1b1[r]->SetLineColor(600-1-r);  halpha_jet1b1[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             halpha_jet1b1[r]->SetLineStyle(2+r);
//             //halpha_jet1b1[r]->SetLineColor(632-1-r); 
//             halpha_jet1b1[2]->SetLineColor(2);  
//             halpha_jet1b1[3]->SetLineColor(3);  
//             halpha_jet1b1[4]->SetLineColor(7);  
//             halpha_jet1b1[5]->SetLineColor(6);  
//             halpha_jet1b1[6]->SetLineColor(kOrange-3);  
//             halpha_jet1b1[7]->SetLineColor(kPink+2);  
//             halpha_jet1b1[8]->SetLineColor(8);  
//             halpha_jet1b1[8]->SetLineStyle(1); 
//         }
//         halpha_jet1b1[r]->GetXaxis()->SetTitle("#alpha(b_{1},l-jet) [Rad]");
//         halpha_jet1b1[r]->GetXaxis()->SetLabelFont(132);
//         halpha_jet1b1[r]->GetXaxis()->SetTitleFont(22);
//         halpha_jet1b1[r]->GetXaxis()->SetTitleOffset(1.10);
//         halpha_jet1b1[r]->GetYaxis()->SetTitle("Normalized distribution");
//         halpha_jet1b1[r]->GetYaxis()->SetLabelFont(132);
//         halpha_jet1b1[r]->GetYaxis()->SetTitleFont(22);
//         halpha_jet1b1[r]->GetYaxis()->SetTitleOffset(1.10);
//         //halpha_jet1b1[r]->GetYaxis()->SetRangeUser(0.0, .2);
// //-----------------
//         halpha_jet1b2[r]->Scale(1 /(halpha_jet1b2[r]->Integral()));
//         halpha_jet1b2[r]->SetLineStyle(2+r);
//         halpha_jet1b2[r]->SetLineWidth(2);
//         if (r < 2){ halpha_jet1b2[r]->SetLineColor(600-1-r);  halpha_jet1b2[r]->SetLineStyle(2+r);}
//         if (r > 1 && r < 9){
//             halpha_jet1b2[r]->SetLineStyle(2+r);
//             //halpha_jet1b2[r]->SetLineColor(632-1-r); 
//             halpha_jet1b2[2]->SetLineColor(2);  
//             halpha_jet1b2[3]->SetLineColor(3);  
//             halpha_jet1b2[4]->SetLineColor(7);  
//             halpha_jet1b2[5]->SetLineColor(6);  
//             halpha_jet1b2[6]->SetLineColor(kOrange-3);  
//             halpha_jet1b2[7]->SetLineColor(kPink+2);  
//             halpha_jet1b2[8]->SetLineColor(8);  
//             halpha_jet1b2[8]->SetLineStyle(1); 
//         }
//         halpha_jet1b2[r]->GetXaxis()->SetTitle("#alpha(b_{2},l-jet) [Rad]");
//         halpha_jet1b2[r]->GetXaxis()->SetLabelFont(132);
//         halpha_jet1b2[r]->GetXaxis()->SetTitleFont(22);
//         halpha_jet1b2[r]->GetXaxis()->SetTitleOffset(1.10);
//         halpha_jet1b2[r]->GetYaxis()->SetTitle("Normalized distribution");
//         halpha_jet1b2[r]->GetYaxis()->SetLabelFont(132);
//         halpha_jet1b2[r]->GetYaxis()->SetTitleFont(22);
//         halpha_jet1b2[r]->GetYaxis()->SetTitleOffset(1.10);
//         //halpha_jet1b2[r]->GetYaxis()->SetRangeUser(0.0, .2);


//         if (r == 0)
//         {
//             c1->cd(1);    hjet1_pt[r]->Draw("hist");    
//             c2->cd(1);    hjet_pt[r]->Draw("hist");
//             c3->cd(1);    hjet2_pt[r]->Draw("hist");

//             c4->cd(1);    hjet1_eta[r]->Draw("hist");
//             c5->cd(1);    hjet_eta[r]->Draw("hist");
//             c6->cd(1);    hjet2_eta[r]->Draw("hist");

//             c7->cd(1);    hjet1_P[r]->Draw("hist");
//             c8->cd(1);    hjet_P[r]->Draw("hist");
//             c9->cd(1);    hjet2_P[r]->Draw("hist");

//             c10->cd(1);    hjetsInvariantMass[r]->Draw("hist");
//             c11->cd(1);    hMET[r]->Draw("hist");
//             c12->cd(1);    hHT[r]->Draw("hist");

//             c13->cd(1);    hdeltaeta_jets[r]->Draw("hist");
//             c14->cd(1);    hdeltaeta_jet1jet[r]->Draw("hist");
//             c15->cd(1);    hdeltaeta_jet1jet2[r]->Draw("hist");

//             c16->cd(1);    hdeltaphi_jets[r]->Draw("hist");
//             c17->cd(1);    hdeltaphi_jet1jet[r]->Draw("hist");
//             c18->cd(1);    hdeltaphi_jet1jet2[r]->Draw("hist");

//             c19->cd(1);    hdeltaR_jets[r]->Draw("hist");
//             c20->cd(1);    hdeltaR_jet1jet[r]->Draw("hist");
//             c21->cd(1);    hdeltaR_jet1jet2[r]->Draw("hist");

//             c22->cd(1);    hcos_jets[r]->Draw("hist");
//             c23->cd(1);    hcos_jet1jet[r]->Draw("hist");
//             c24->cd(1);    hcos_jet1jet2[r]->Draw("hist");

//             c25->cd(1);    hjetmultiplicity[r]->Draw("hist");

//             c26->cd(1);    halpha[r]->Draw("hist");
//             c27->cd(1);    halphaZMF[r]->Draw("hist");
//             c28->cd(1);    halpha_jet1b1[r]->Draw("hist");
//             c29->cd(1);    halpha_jet1b2[r]->Draw("hist");

//         }
//         else
//         {
//             c1->cd(1);    hjet1_pt[r]->Draw("hist same");
//             c2->cd(1);    hjet_pt[r]->Draw("hist same");
//             c3->cd(1);    hjet2_pt[r]->Draw("hist same");

//             c4->cd(1);    hjet1_eta[r]->Draw("hist same");
//             c5->cd(1);    hjet_eta[r]->Draw("hist same");
//             c6->cd(1);    hjet2_eta[r]->Draw("hist same");

//             c7->cd(1);    hjet1_P[r]->Draw("hist same");
//             c8->cd(1);    hjet_P[r]->Draw("hist same");
//             c9->cd(1);    hjet2_P[r]->Draw("hist same");

//             c10->cd(1);    hjetsInvariantMass[r]->Draw("hist same");
//             c11->cd(1);    hMET[r]->Draw("hist same");
//             c12->cd(1);    hHT[r]->Draw("hist same");

//             c13->cd(1);    hdeltaeta_jets[r]->Draw("hist same");
//             c14->cd(1);    hdeltaeta_jet1jet[r]->Draw("hist same");
//             c15->cd(1);    hdeltaeta_jet1jet2[r]->Draw("hist same");

//             c16->cd(1);    hdeltaphi_jets[r]->Draw("hist same");
//             c17->cd(1);    hdeltaphi_jet1jet[r]->Draw("hist same");
//             c18->cd(1);    hdeltaphi_jet1jet2[r]->Draw("hist same");

//             c19->cd(1);    hdeltaR_jets[r]->Draw("hist same");
//             c20->cd(1);    hdeltaR_jet1jet[r]->Draw("hist same");
//             c21->cd(1);    hdeltaR_jet1jet2[r]->Draw("hist same");

//             c22->cd(1);    hcos_jets[r]->Draw("hist same");
//             c23->cd(1);    hcos_jet1jet[r]->Draw("hist same");
//             c24->cd(1);    hcos_jet1jet2[r]->Draw("hist same");

//             c25->cd(1);    hjetmultiplicity[r]->Draw("hist same");

//             c26->cd(1);    halpha[r]->Draw("hist same");
//             c27->cd(1);    halphaZMF[r]->Draw("hist same");
//             c28->cd(1);    halpha_jet1b1[r]->Draw("hist same");
//             c29->cd(1);    halpha_jet1b2[r]->Draw("hist same");

//         }

//     } // other loop over samples

//         c1->cd(1);
//         TLegend *legend1 = new TLegend(0.65,0.4,0.83,0.87);// or: (0.73,0.6,0.87,0.87)
//         legend1->SetMargin(0.5); 
//         legend1->SetNColumns(1);
//         //legend1->AddEntry(hjet1_pt[0], "#font[12]{xtar}", "l"); //#font[132]
//         legend1->AddEntry(hjet1_pt[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend1->AddEntry(hjet1_pt[2], "#font[12]{xtr}", "l");
//         //legend1->AddEntry(hjet1_pt[3], "#font[12]{xtai}", "l");
//         legend1->AddEntry(hjet1_pt[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend1->AddEntry(hjet1_pt[5], "#font[12]{xti}", "l");
//         legend1->AddEntry(hjet1_pt[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[7], "#font[12]{tj#nu_{e}}", "l");
//         legend1->AddEntry(hjet1_pt[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend1->SetLineStyle(1);
//         legend1->SetLineWidth(1);
//         legend1->SetTextSize(0.035);
//         //legend1->SetTextSize(0.045);
//         legend1->SetLineColor(10);
//         legend1->SetFillColor(10);
//         legend1->Draw("hist");

//         c2->cd(1);
//         TLegend *legend2 = new TLegend(0.65,0.4,0.83,0.87);
//         legend2->SetMargin(0.5);
//         legend2->SetNColumns(1);
//         //legend2->AddEntry(hjet_pt[0], "#font[12]{xtar}", "l");
//         legend2->AddEntry(hjet_pt[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend2->AddEntry(hjet_pt[2], "#font[12]{xtr}", "l");
//         //legend2->AddEntry(hjet_pt[3], "#font[12]{xtai}", "l");
//         legend2->AddEntry(hjet_pt[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend2->AddEntry(hjet_pt[5], "#font[12]{xti}", "l");
//         legend2->AddEntry(hjet_pt[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[7], "#font[12]{tj#nu_{e}}", "l");
//         legend2->AddEntry(hjet_pt[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend2->SetLineStyle(1);
//         legend2->SetLineWidth(1);
//         legend2->SetTextSize(0.035);
//         legend2->SetLineColor(10);
//         legend2->SetFillColor(10);
//         legend2->Draw("hist");

//         c3->cd(1);
//         TLegend *legend3 = new TLegend(0.65,0.4,0.83,0.87);
//         legend3->SetMargin(0.5);
//         legend3->SetNColumns(1);
//         //legend3->AddEntry(hjet2_pt[0], "#font[12]{xtar}", "l");
//         legend3->AddEntry(hjet2_pt[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend3->AddEntry(hjet2_pt[2], "#font[12]{xtr}", "l");
//         //legend3->AddEntry(hjet2_pt[3], "#font[12]{xtai}", "l");
//         legend3->AddEntry(hjet2_pt[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend3->AddEntry(hjet2_pt[5], "#font[12]{xti}", "l");
//         legend3->AddEntry(hjet2_pt[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[7], "#font[12]{tj#nu_{e}}", "l");
//         legend3->AddEntry(hjet2_pt[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend3->SetLineStyle(1);
//         legend3->SetLineWidth(1);
//         legend3->SetTextSize(0.035);
//         legend3->SetLineColor(10);
//         legend3->SetFillColor(10);
//         legend3->Draw("hist");

//         c4->cd(1);
//         TLegend *legend4 = new TLegend(0.65,0.4,0.83,0.87);
//         legend4->SetMargin(0.5);
//         legend4->SetNColumns(1);
//         //legend4->AddEntry(hjet1_eta[0], "#font[12]{xtar}", "l");
//         legend4->AddEntry(hjet1_eta[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend4->AddEntry(hjet1_eta[2], "#font[12]{xtr}", "l");
//         //legend4->AddEntry(hjet1_eta[3], "#font[12]{xtai}", "l");
//         legend4->AddEntry(hjet1_eta[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend4->AddEntry(hjet1_eta[5], "#font[12]{xti}", "l");
//         legend4->AddEntry(hjet1_eta[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[7], "#font[12]{tj#nu_{e}}", "l");
//         legend4->AddEntry(hjet1_eta[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend4->SetLineStyle(1);
//         legend4->SetLineWidth(1);
//         legend4->SetTextSize(0.035);
//         legend4->SetLineColor(10);
//         legend4->SetFillColor(10);
//         legend4->Draw("hist");

//         c5->cd(1);
//         TLegend *legend5 = new TLegend(0.65,0.4,0.83,0.87);
//         legend5->SetMargin(0.5);
//         legend5->SetNColumns(1);
//         //legend5->AddEntry(hjet_eta[0], "#font[12]{xtar}", "l");
//         legend5->AddEntry(hjet_eta[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend5->AddEntry(hjet_eta[2], "#font[12]{xtr}", "l");
//         //legend5->AddEntry(hjet_eta[3], "#font[12]{xtai}", "l");
//         legend5->AddEntry(hjet_eta[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend5->AddEntry(hjet_eta[5], "#font[12]{xti}", "l");
//         legend5->AddEntry(hjet_eta[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[7], "#font[12]{tj#nu_{e}}", "l");
//         legend5->AddEntry(hjet_eta[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend5->SetLineStyle(1);
//         legend5->SetLineWidth(1);
//         legend5->SetTextSize(0.035);
//         legend5->SetLineColor(10);
//         legend5->SetFillColor(10);
//         legend5->Draw("hist");

//         c6->cd(1);
//         TLegend *legend6 = new TLegend(0.65,0.4,0.83,0.87);
//         legend6->SetMargin(0.5);
//         legend6->SetNColumns(1);
//         //legend6->AddEntry(hjet2_eta[0], "#font[12]{xtar}", "l");
//         legend6->AddEntry(hjet2_eta[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend6->AddEntry(hjet2_eta[2], "#font[12]{xtr}", "l");
//         //legend6->AddEntry(hjet2_eta[3], "#font[12]{xtai}", "l");
//         legend6->AddEntry(hjet2_eta[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend6->AddEntry(hjet2_eta[5], "#font[12]{xti}", "l");
//         legend6->AddEntry(hjet2_eta[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[7], "#font[12]{tj#nu_{e}}", "l");
//         legend6->AddEntry(hjet2_eta[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend6->SetLineStyle(1);
//         legend6->SetLineWidth(1);
//         legend6->SetTextSize(0.035);
//         legend6->SetLineColor(10);
//         legend6->SetFillColor(10);
//         legend6->Draw("hist");


//         c7->cd(1);
//         TLegend *legend7 = new TLegend(0.65,0.4,0.83,0.87);
//         legend7->SetMargin(0.5);
//         legend7->SetNColumns(1);
//         //legend7->AddEntry(hjet1_P[0], "#font[12]{xtar}", "l");
//         legend7->AddEntry(hjet1_P[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend7->AddEntry(hjet1_P[2], "#font[12]{xtr}", "l");
//         //legend7->AddEntry(hjet1_P[3], "#font[12]{xtai}", "l");
//         legend7->AddEntry(hjet1_P[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend7->AddEntry(hjet1_P[5], "#font[12]{xti}", "l");
//         legend7->AddEntry(hjet1_P[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[7], "#font[12]{tj#nu_{e}}", "l");
//         legend7->AddEntry(hjet1_P[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend7->SetLineStyle(1);
//         legend7->SetLineWidth(1);
//         legend7->SetTextSize(0.035);
//         legend7->SetLineColor(10);
//         legend7->SetFillColor(10);
//         legend7->Draw("hist");

//         c8->cd(1);
//         TLegend *legend8 = new TLegend(0.65,0.4,0.83,0.87);
//         legend8->SetMargin(0.5);
//         legend8->SetNColumns(1);
//         //legend8->AddEntry(hjet_P[0], "#font[12]{xtar}", "l");
//         legend8->AddEntry(hjet_P[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend8->AddEntry(hjet_P[2], "#font[12]{xtr}", "l");
//         //legend8->AddEntry(hjet_P[3], "#font[12]{xtai}", "l");
//         legend8->AddEntry(hjet_P[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend8->AddEntry(hjet_P[5], "#font[12]{xti}", "l");
//         legend8->AddEntry(hjet_P[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[7], "#font[12]{tj#nu_{e}}", "l");
//         legend8->AddEntry(hjet_P[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend8->SetLineStyle(1);
//         legend8->SetLineWidth(1);
//         legend8->SetTextSize(0.035);
//         legend8->SetLineColor(10);
//         legend8->SetFillColor(10);
//         legend8->Draw("hist");

//         c9->cd(1);
//         TLegend *legend9 = new TLegend(0.65,0.4,0.83,0.87);
//         legend9->SetMargin(0.5);
//         legend9->SetNColumns(1);
//         //legend9->AddEntry(hjet2_P[0], "#font[12]{xtar}", "l");
//         legend9->AddEntry(hjet2_P[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend9->AddEntry(hjet2_P[2], "#font[12]{xtr}", "l");
//         //legend9->AddEntry(hjet2_P[3], "#font[12]{xtai}", "l");
//         legend9->AddEntry(hjet2_P[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend9->AddEntry(hjet2_P[5], "#font[12]{xti}", "l");
//         legend9->AddEntry(hjet2_P[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[7], "#font[12]{tj#nu_{e}}", "l");
//         legend9->AddEntry(hjet2_P[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend9->SetLineStyle(1);
//         legend9->SetLineWidth(1);
//         legend9->SetTextSize(0.035);
//         legend9->SetLineColor(10);
//         legend9->SetFillColor(10);
//         legend9->Draw("hist");

//         c10->cd(1);
//         TLegend *legend10 = new TLegend(0.65,0.4,0.83,0.87);
//         legend10->SetMargin(0.5);
//         legend10->SetNColumns(1);
//         //legend10->AddEntry(hjetsInvariantMass[0], "#font[12]{xtar}", "l");
//         legend10->AddEntry(hjetsInvariantMass[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend10->AddEntry(hjetsInvariantMass[2], "#font[12]{xtr}", "l");
//         //legend10->AddEntry(hjetsInvariantMass[3], "#font[12]{xtai}", "l");
//         legend10->AddEntry(hjetsInvariantMass[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend10->AddEntry(hjetsInvariantMass[5], "#font[12]{xti}", "l");
//         legend10->AddEntry(hjetsInvariantMass[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[7], "#font[12]{tj#nu_{e}}", "l");
//         legend10->AddEntry(hjetsInvariantMass[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend10->SetLineStyle(1);
//         legend10->SetLineWidth(1);
//         legend10->SetTextSize(0.035);
//         legend10->SetLineColor(10);
//         legend10->SetFillColor(10);
//         legend10->Draw("hist");

//         c11->cd(1);
//         TLegend *legend11 = new TLegend(0.65,0.4,0.83,0.87);
//         legend11->SetMargin(0.5);
//         legend11->SetNColumns(1);
//         //legend11->AddEntry(hMET[0], "#font[12]{xtar}", "l");
//         legend11->AddEntry(hMET[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend11->AddEntry(hMET[2], "#font[12]{xtr}", "l");
//         //legend11->AddEntry(hMET[3], "#font[12]{xtai}", "l");
//         legend11->AddEntry(hMET[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend11->AddEntry(hMET[5], "#font[12]{xti}", "l");
//         legend11->AddEntry(hMET[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend11->AddEntry(hMET[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend11->AddEntry(hMET[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend11->AddEntry(hMET[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend11->AddEntry(hMET[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend11->AddEntry(hMET[7], "#font[12]{tj#nu_{e}}", "l");
//         legend11->AddEntry(hMET[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend11->SetLineStyle(1);
//         legend11->SetLineWidth(1);
//         legend11->SetTextSize(0.035);
//         legend11->SetLineColor(10);
//         legend11->SetFillColor(10);
//         legend11->Draw("hist");


//         c12->cd(1);
//         TLegend *legend12 = new TLegend(0.65,0.4,0.83,0.87);
//         legend12->SetMargin(0.5);
//         legend12->SetNColumns(1);
//         //legend12->AddEntry(hHT[0], "#font[12]{xtar}", "l");
//         legend12->AddEntry(hHT[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend12->AddEntry(hHT[2], "#font[12]{xtr}", "l");
//         //legend12->AddEntry(hHT[3], "#font[12]{xtai}", "l");
//         legend12->AddEntry(hHT[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend12->AddEntry(hHT[5], "#font[12]{xti}", "l");
//         legend12->AddEntry(hHT[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend12->AddEntry(hHT[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend12->AddEntry(hHT[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend12->AddEntry(hHT[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend12->AddEntry(hHT[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend12->AddEntry(hHT[7], "#font[12]{tj#nu_{e}}", "l");
//         legend12->AddEntry(hHT[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend12->SetLineStyle(1);
//         legend12->SetLineWidth(1);
//         legend12->SetTextSize(0.035);
//         legend12->SetLineColor(10);
//         legend12->SetFillColor(10);
//         legend12->Draw("hist");

//         c13->cd(1);
//         TLegend *legend13 = new TLegend(0.65,0.4,0.83,0.87);
//         legend13->SetMargin(0.5);
//         legend13->SetNColumns(1);
//         //legend13->AddEntry(hdeltaeta_jets[0], "#font[12]{xtar}", "l");
//         legend13->AddEntry(hdeltaeta_jets[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend13->AddEntry(hdeltaeta_jets[2], "#font[12]{xtr}", "l");
//         //legend13->AddEntry(hdeltaeta_jets[3], "#font[12]{xtai}", "l");
//         legend13->AddEntry(hdeltaeta_jets[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend13->AddEntry(hdeltaeta_jets[5], "#font[12]{xti}", "l");
//         legend13->AddEntry(hdeltaeta_jets[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[7], "#font[12]{tj#nu_{e}}", "l");
//         legend13->AddEntry(hdeltaeta_jets[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend13->SetLineStyle(1);
//         legend13->SetLineWidth(1);
//         legend13->SetTextSize(0.035);
//         legend13->SetLineColor(10);
//         legend13->SetFillColor(10);
//         legend13->Draw("hist");

//         c14->cd(1);
//         TLegend *legend14 = new TLegend(0.65,0.4,0.83,0.87);
//         legend14->SetMargin(0.5);
//         legend14->SetNColumns(1);
//         //legend14->AddEntry(hdeltaeta_jet1jet[0], "#font[12]{xtar}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend14->AddEntry(hdeltaeta_jet1jet[2], "#font[12]{xtr}", "l");
//         //legend14->AddEntry(hdeltaeta_jet1jet[3], "#font[12]{xtai}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend14->AddEntry(hdeltaeta_jet1jet[5], "#font[12]{xti}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[7], "#font[12]{tj#nu_{e}}", "l");
//         legend14->AddEntry(hdeltaeta_jet1jet[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend14->SetLineStyle(1);
//         legend14->SetLineWidth(1);
//         legend14->SetTextSize(0.035);
//         legend14->SetLineColor(10);
//         legend14->SetFillColor(10);
//         legend14->Draw("hist");

//         c15->cd(1);
//         TLegend *legend15 = new TLegend(0.65,0.4,0.83,0.87);
//         legend15->SetMargin(0.5);
//         legend15->SetNColumns(1);
//         //legend15->AddEntry(hdeltaeta_jet1jet2[0], "#font[12]{xtar}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend15->AddEntry(hdeltaeta_jet1jet2[2], "#font[12]{xtr}", "l");
//         //legend15->AddEntry(hdeltaeta_jet1jet2[3], "#font[12]{xtai}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend15->AddEntry(hdeltaeta_jet1jet2[5], "#font[12]{xti}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[7], "#font[12]{tj#nu_{e}}", "l");
//         legend15->AddEntry(hdeltaeta_jet1jet2[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend15->SetLineStyle(1);
//         legend15->SetLineWidth(1);
//         legend15->SetTextSize(0.035);
//         legend15->SetLineColor(10);
//         legend15->SetFillColor(10);
//         legend15->Draw("hist");

//         c16->cd(1);
//         TLegend *legend16 = new TLegend(0.65,0.4,0.83,0.87);
//         legend16->SetMargin(0.5);
//         legend16->SetNColumns(1);
//         //legend16->AddEntry(hdeltaphi_jets[0], "#font[12]{xtar}", "l");
//         legend16->AddEntry(hdeltaphi_jets[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend16->AddEntry(hdeltaphi_jets[2], "#font[12]{xtr}", "l");
//         //legend16->AddEntry(hdeltaphi_jets[3], "#font[12]{xtai}", "l");
//         legend16->AddEntry(hdeltaphi_jets[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend16->AddEntry(hdeltaphi_jets[5], "#font[12]{xti}", "l");
//         legend16->AddEntry(hdeltaphi_jets[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[7], "#font[12]{tj#nu_{e}}", "l");
//         legend16->AddEntry(hdeltaphi_jets[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend16->SetLineStyle(1);
//         legend16->SetLineWidth(1);
//         legend16->SetTextSize(0.035);
//         legend16->SetLineColor(10);
//         legend16->SetFillColor(10);
//         legend16->Draw("hist");


//         c17->cd(1);
//         TLegend *legend17 = new TLegend(0.65,0.4,0.83,0.87);
//         legend17->SetMargin(0.5);
//         legend17->SetNColumns(1);
//         //legend17->AddEntry(hdeltaphi_jet1jet[0], "#font[12]{xtar}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend17->AddEntry(hdeltaphi_jet1jet[2], "#font[12]{xtr}", "l");
//         //legend17->AddEntry(hdeltaphi_jet1jet[3], "#font[12]{xtai}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend17->AddEntry(hdeltaphi_jet1jet[5], "#font[12]{xti}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[7], "#font[12]{tj#nu_{e}}", "l");
//         legend17->AddEntry(hdeltaphi_jet1jet[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend17->SetLineStyle(1);
//         legend17->SetLineWidth(1);
//         legend17->SetTextSize(0.035);
//         legend17->SetLineColor(10);
//         legend17->SetFillColor(10);
//         legend17->Draw("hist");

//         c18->cd(1);
//         TLegend *legend18 = new TLegend(0.65,0.4,0.83,0.87);
//         legend18->SetMargin(0.5);
//         legend18->SetNColumns(1);
//         //legend18->AddEntry(hdeltaphi_jet1jet2[0], "#font[12]{xtar}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend18->AddEntry(hdeltaphi_jet1jet2[2], "#font[12]{xtr}", "l");
//         //legend18->AddEntry(hdeltaphi_jet1jet2[3], "#font[12]{xtai}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend18->AddEntry(hdeltaphi_jet1jet2[5], "#font[12]{xti}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[7], "#font[12]{tj#nu_{e}}", "l");
//         legend18->AddEntry(hdeltaphi_jet1jet2[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend18->SetLineStyle(1);
//         legend18->SetLineWidth(1);
//         legend18->SetTextSize(0.035);
//         legend18->SetLineColor(10);
//         legend18->SetFillColor(10);
//         legend18->Draw("hist");

//         c19->cd(1);
//         TLegend *legend19 = new TLegend(0.65,0.4,0.83,0.87);
//         legend19->SetMargin(0.5);
//         legend19->SetNColumns(1);
//         //legend19->AddEntry(hdeltaR_jets[0], "#font[12]{xtar}", "l");
//         legend19->AddEntry(hdeltaR_jets[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend19->AddEntry(hdeltaR_jets[2], "#font[12]{xtr}", "l");
//         //legend19->AddEntry(hdeltaR_jets[3], "#font[12]{xtai}", "l");
//         legend19->AddEntry(hdeltaR_jets[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend19->AddEntry(hdeltaR_jets[5], "#font[12]{xti}", "l");
//         legend19->AddEntry(hdeltaR_jets[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[7], "#font[12]{tj#nu_{e}}", "l");
//         legend19->AddEntry(hdeltaR_jets[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend19->SetLineStyle(1);
//         legend19->SetLineWidth(1);
//         legend19->SetTextSize(0.035);
//         legend19->SetLineColor(10);
//         legend19->SetFillColor(10);
//         legend19->Draw("hist");

//         c20->cd(1);
//         TLegend *legend20 = new TLegend(0.65,0.4,0.83,0.87);
//         legend20->SetMargin(0.5);
//         legend20->SetNColumns(1);
//         //legend20->AddEntry(hdeltaR_jet1jet[0], "#font[12]{xtar}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend20->AddEntry(hdeltaR_jet1jet[2], "#font[12]{xtr}", "l");
//         //legend20->AddEntry(hdeltaR_jet1jet[3], "#font[12]{xtai}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend20->AddEntry(hdeltaR_jet1jet[5], "#font[12]{xti}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[7], "#font[12]{tj#nu_{e}}", "l");
//         legend20->AddEntry(hdeltaR_jet1jet[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend20->SetLineStyle(1);
//         legend20->SetLineWidth(1);
//         legend20->SetTextSize(0.035);
//         legend20->SetLineColor(10);
//         legend20->SetFillColor(10);
//         legend20->Draw("hist");

//         c21->cd(1);
//         TLegend *legend21 = new TLegend(0.65,0.4,0.83,0.87);
//         legend21->SetMargin(0.5);
//         legend21->SetNColumns(1);
//         //legend21->AddEntry(hdeltaR_jet1jet2[0], "#font[12]{xtar}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend21->AddEntry(hdeltaR_jet1jet2[2], "#font[12]{xtr}", "l");
//         //legend21->AddEntry(hdeltaR_jet1jet2[3], "#font[12]{xtai}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend21->AddEntry(hdeltaR_jet1jet2[5], "#font[12]{xti}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[7], "#font[12]{tj#nu_{e}}", "l");
//         legend21->AddEntry(hdeltaR_jet1jet2[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend21->SetLineStyle(1);
//         legend21->SetLineWidth(1);
//         legend21->SetTextSize(0.035);
//         legend21->SetLineColor(10);
//         legend21->SetFillColor(10);
//         legend21->Draw("hist");


//         c22->cd(1);
//         TLegend *legend22 = new TLegend(0.65,0.4,0.83,0.87);
//         legend22->SetMargin(0.5);
//         legend22->SetNColumns(1);
//         //legend22->AddEntry(hcos_jets[0], "#font[12]{xtar}", "l");
//         legend22->AddEntry(hcos_jets[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend22->AddEntry(hcos_jets[2], "#font[12]{xtr}", "l");
//         //legend22->AddEntry(hcos_jets[3], "#font[12]{xtai}", "l");
//         legend22->AddEntry(hcos_jets[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend22->AddEntry(hcos_jets[5], "#font[12]{xti}", "l");
//         legend22->AddEntry(hcos_jets[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[7], "#font[12]{tj#nu_{e}}", "l");
//         legend22->AddEntry(hcos_jets[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend22->SetLineStyle(1);
//         legend22->SetLineWidth(1);
//         legend22->SetTextSize(0.035);
//         legend22->SetLineColor(10);
//         legend22->SetFillColor(10);
//         legend22->Draw("hist");

//         c23->cd(1);
//         TLegend *legend23 = new TLegend(0.65,0.4,0.83,0.87);
//         legend23->SetMargin(0.5);
//         legend23->SetNColumns(1);
//         //legend23->AddEntry(hcos_jet1jet[0], "#font[12]{xtar}", "l");
//         legend23->AddEntry(hcos_jet1jet[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend23->AddEntry(hcos_jet1jet[2], "#font[12]{xtr}", "l");
//         //legend23->AddEntry(hcos_jet1jet[3], "#font[12]{xtai}", "l");
//         legend23->AddEntry(hcos_jet1jet[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend23->AddEntry(hcos_jet1jet[5], "#font[12]{xti}", "l");
//         legend23->AddEntry(hcos_jet1jet[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[7], "#font[12]{tj#nu_{e}}", "l");
//         legend23->AddEntry(hcos_jet1jet[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend23->SetLineStyle(1);
//         legend23->SetLineWidth(1);
//         legend23->SetTextSize(0.035);
//         legend23->SetLineColor(10);
//         legend23->SetFillColor(10);
//         legend23->Draw("hist");

//         c24->cd(1);
//         TLegend *legend24 = new TLegend(0.65,0.4,0.83,0.87);
//         legend24->SetMargin(0.5);
//         legend24->SetNColumns(1);
//         //legend24->AddEntry(hcos_jet1jet2[0], "#font[12]{xtar}", "l");
//         legend24->AddEntry(hcos_jet1jet2[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend24->AddEntry(hcos_jet1jet2[2], "#font[12]{xtr}", "l");
//         //legend24->AddEntry(hcos_jet1jet2[3], "#font[12]{xtai}", "l");
//         legend24->AddEntry(hcos_jet1jet2[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend24->AddEntry(hcos_jet1jet2[5], "#font[12]{xti}", "l");
//         legend24->AddEntry(hcos_jet1jet2[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[7], "#font[12]{tj#nu_{e}}", "l");
//         legend24->AddEntry(hcos_jet1jet2[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend24->SetLineStyle(1);
//         legend24->SetLineWidth(1);
//         legend24->SetTextSize(0.035);
//         legend24->SetLineColor(10);
//         legend24->SetFillColor(10);
//         legend24->Draw("hist");

//         c25->cd(1);
//         TLegend *legend25 = new TLegend(0.65,0.4,0.83,0.87);
//         legend25->SetMargin(0.5);
//         legend25->SetNColumns(1);
//         //legend25->AddEntry(hjetmultiplicity[0], "#font[12]{xtar}", "l");
//         legend25->AddEntry(hjetmultiplicity[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend25->AddEntry(hjetmultiplicity[2], "#font[12]{xtr}", "l");
//         //legend25->AddEntry(hjetmultiplicity[3], "#font[12]{xtai}", "l");
//         legend25->AddEntry(hjetmultiplicity[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend25->AddEntry(hjetmultiplicity[5], "#font[12]{xti}", "l");
//         legend25->AddEntry(hjetmultiplicity[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[7], "#font[12]{tj#nu_{e}}", "l");
//         legend25->AddEntry(hjetmultiplicity[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend25->SetLineStyle(1);
//         legend25->SetLineWidth(1);
//         legend25->SetTextSize(0.035);
//         legend25->SetLineColor(10);
//         legend25->SetFillColor(10);
//         legend25->Draw("hist");

//         c26->cd(1);
//         TLegend *legend26 = new TLegend(0.65,0.4,0.83,0.87);
//         legend26->SetMargin(0.5);
//         legend26->SetNColumns(1);
//         //legend26->AddEntry(halpha[0], "#font[12]{xtar}", "l");
//         legend26->AddEntry(halpha[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend26->AddEntry(halpha[2], "#font[12]{xtr}", "l");
//         //legend26->AddEntry(halpha[3], "#font[12]{xtai}", "l");
//         legend26->AddEntry(halpha[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend26->AddEntry(halpha[5], "#font[12]{xti}", "l");
//         legend26->AddEntry(halpha[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend26->AddEntry(halpha[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend26->AddEntry(halpha[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend26->AddEntry(halpha[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend26->AddEntry(halpha[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend26->AddEntry(halpha[7], "#font[12]{tj#nu_{e}}", "l");
//         legend26->AddEntry(halpha[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend26->SetLineStyle(1);
//         legend26->SetLineWidth(1);
//         legend26->SetTextSize(0.035);
//         legend26->SetLineColor(10);
//         legend26->SetFillColor(10);
//         legend26->Draw("hist");

//         c27->cd(1);
//         TLegend *legend27 = new TLegend(0.65,0.4,0.83,0.87);
//         legend27->SetMargin(0.5);
//         legend27->SetNColumns(1);
//         //legend27->AddEntry(halphaZMF[0], "#font[12]{xtar}", "l");
//         legend27->AddEntry(halphaZMF[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend27->AddEntry(halphaZMF[2], "#font[12]{xtr}", "l");
//         //legend27->AddEntry(halphaZMF[3], "#font[12]{xtai}", "l");
//         legend27->AddEntry(halphaZMF[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend27->AddEntry(halphaZMF[5], "#font[12]{xti}", "l");
//         legend27->AddEntry(halphaZMF[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[7], "#font[12]{tj#nu_{e}}", "l");
//         legend27->AddEntry(halphaZMF[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend27->SetLineStyle(1);
//         legend27->SetLineWidth(1);
//         legend27->SetTextSize(0.035);
//         legend27->SetLineColor(10);
//         legend27->SetFillColor(10);
//         legend27->Draw("hist");

//         c28->cd(1);
//         TLegend *legend28 = new TLegend(0.65,0.4,0.83,0.87);
//         legend28->SetMargin(0.5);
//         legend28->SetNColumns(1);
//         //legend28->AddEntry(halpha_jet1b1[0], "#font[12]{xtar}", "l");
//         legend28->AddEntry(halpha_jet1b1[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend28->AddEntry(halpha_jet1b1[2], "#font[12]{xtr}", "l");
//         //legend28->AddEntry(halpha_jet1b1[3], "#font[12]{xtai}", "l");
//         legend28->AddEntry(halpha_jet1b1[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend28->AddEntry(halpha_jet1b1[5], "#font[12]{xti}", "l");
//         legend28->AddEntry(halpha_jet1b1[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[7], "#font[12]{tj#nu_{e}}", "l");
//         legend28->AddEntry(halpha_jet1b1[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend28->SetLineStyle(1);
//         legend28->SetLineWidth(1);
//         legend28->SetTextSize(0.035);
//         legend28->SetLineColor(10);
//         legend28->SetFillColor(10);
//         legend28->Draw("hist");

//         c29->cd(1);
//         TLegend *legend29 = new TLegend(0.65,0.4,0.83,0.87);
//         legend29->SetMargin(0.5);
//         legend29->SetNColumns(1);
//         //legend29->AddEntry(halpha_jet1b2[0], "#font[12]{xtar}", "l");
//         legend29->AddEntry(halpha_jet1b2[0], "#font[12]{H+j#nu_{e}(T_{R})}", "l");
//         //legend29->AddEntry(halpha_jet1b2[2], "#font[12]{xtr}", "l");
//         //legend29->AddEntry(halpha_jet1b2[3], "#font[12]{xtai}", "l");
//         legend29->AddEntry(halpha_jet1b2[1], "#font[12]{H+j#nu_{e}(T_{I})}", "l");
//         //legend29->AddEntry(halpha_jet1b2[5], "#font[12]{xti}", "l");
//         legend29->AddEntry(halpha_jet1b2[2], "#font[12]{b#bar{b}+j#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[3], "#font[12]{c#bar{c}+j#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[4], "#font[12]{Z(b#bar{b})+j#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[5], "#font[12]{Z(c#bar{c})+j#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[6], "#font[12]{Z(j_{l}#bar{j_{l}})+j#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[7], "#font[12]{tj#nu_{e}}", "l");
//         legend29->AddEntry(halpha_jet1b2[8], "#font[12]{H+j#nu_{e}}", "l");
//         legend29->SetLineStyle(1);
//         legend29->SetLineWidth(1);
//         legend29->SetTextSize(0.035);
//         legend29->SetLineColor(10);
//         legend29->SetFillColor(10);
//         legend29->Draw("hist");


// c1->SaveAs("plots1/jet1_pt.pdf");c1->SaveAs("plots1/jet1_pt.C");
// c2->SaveAs("plots1/jet_pt.pdf");c2->SaveAs("plots1/jet_pt.C");
// c3->SaveAs("plots1/jet2_pt.pdf");c3->SaveAs("plots1/jet2_pt.C");
// c4->SaveAs("plots1/jet1_eta.pdf");c4->SaveAs("plots1/jet1_eta.C");
// c5->SaveAs("plots1/jet_eta.pdf");c5->SaveAs("plots1/jet_eta.C");
// c6->SaveAs("plots1/jet2_eta.pdf");c6->SaveAs("plots1/jet2_eta.C");
// c7->SaveAs("plots1/jet1_P.pdf");c7->SaveAs("plots1/jet1_P.C");
// c8->SaveAs("plots1/jet_P.pdf");c8->SaveAs("plots1/jet_P.C");
// c9->SaveAs("plots1/jet2_P.pdf");c9->SaveAs("plots1/jet2_P.C");
// c10->SaveAs("plots1/jetsInvariantMass.pdf");c10->SaveAs("plots1/jetsInvariantMass.C");
// c11->SaveAs("plots1/MET.pdf");c11->SaveAs("plots1/MET.C");
// c12->SaveAs("plots1/HT.pdf");c12->SaveAs("plots1/HT.C");
// c13->SaveAs("plots1/deltaeta_jets.pdf");c13->SaveAs("plots1/deltaeta_jets.C");
// c14->SaveAs("plots1/deltaeta_jet1jet.pdf");c14->SaveAs("plots1/deltaeta_jet1jet.C");
// c15->SaveAs("plots1/deltaeta_jet1jet2.pdf");c15->SaveAs("plots1/deltaeta_jet1jet2.C");
// c16->SaveAs("plots1/deltaphi_jets.pdf");c16->SaveAs("plots1/deltaphi_jets.C");
// c17->SaveAs("plots1/deltaphi_jet1jet.pdf");c17->SaveAs("plots1/deltaphi_jet1jet.C");
// c18->SaveAs("plots1/deltaphi_jet1jet2.pdf");c18->SaveAs("plots1/deltaphi_jet1jet2.C");
// c19->SaveAs("plots1/deltaR_jets.pdf");c19->SaveAs("plots1/deltaR_jets.C");
// c20->SaveAs("plots1/deltaR_jet1jet.pdf");c20->SaveAs("plots1/deltaR_jet1jet.C");
// c21->SaveAs("plots1/deltaR_jet1jet2.pdf");c21->SaveAs("plots1/deltaR_jet1jet2.C");
// c22->SaveAs("plots1/cos_jets.pdf");c22->SaveAs("plots1/cos_jets.C");
// c23->SaveAs("plots1/cos_jet1jet.pdf");c23->SaveAs("plots1/cos_jet1jet.C");
// c24->SaveAs("plots1/cos_jet1jet2.pdf");c24->SaveAs("plots1/cos_jet1jet2.C");
// c25->SaveAs("plots1/jetmultiplicity.pdf");c25->SaveAs("plots1/jetmultiplicity.C");
// c26->SaveAs("plots1/alpha.pdf");c26->SaveAs("plots1/alpha.C");
// c27->SaveAs("plots1/alphaCOM.pdf");c27->SaveAs("plots1/alphaCOM.C");
// c28->SaveAs("plots1/alpha_jet1b1.pdf");c28->SaveAs("plots1/alpha_jet1b1.C");
// c29->SaveAs("plots1/alpha_jet1b2.pdf");c29->SaveAs("plots1/alpha_jet1b2.C");


cout << endl;
cout << "=================================" << endl;
cout << "============= Finish ============" << endl;
cout << "=================================" << endl;

}//======================================================= SG scope



