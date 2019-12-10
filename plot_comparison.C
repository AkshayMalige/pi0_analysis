
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"


#include <iostream>
#include <map>
#include <stdio.h>

//
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLine.h"
#include "/u/sudol/unigen/v2.1/drawStyle.C"

void plot_comparison()
{
   drawStyle();

    gStyle->SetOptTitle(0);

    TFile *in = new TFile("diphoton_massEXP_gen2_opAng10_50.root","read");


    TFile *out    = new TFile("pi0sec_2_4comparision.root","RECREATE");
    out->cd();


    // ------>  some 2D histograms which are showing phps of gamma e+e-


    //  ------>   histograms for mass spectra in pt-y bins  <--------------- //

    Int_t ptbins = 12;
    Int_t ptmin  = 200;
    Int_t ptmax  = 800;

    Int_t ybins    = 8;
    Double_t ymin  = 0.4;
    Double_t ymax  = 2.;
    
    Double_t PtBinSize = (ptmax-ptmin)/ptbins; // MeV/bin
    Double_t YBinSize  = (ymax-ymin)/ybins;  // rapidity/bin

    TH1F *MassPtY[2][12][8];
    TH1F *MassPtY_Less2[2][12][8];
    TH1F *MassPtY_More2[2][12][8];



   // MassPtY_More2

    Int_t binMin1;
    Int_t binMax1;

    for (Int_t c = 0; c < 2; c++)
    {
        for (Int_t pt = 0; pt < 12; pt++)
        {
            for (Int_t y = 0; y < 8; y++)
            {
                MassPtY[c][pt][y]   = (TH1F*) in->Get(Form("Mass_pt%i_y%i_c%i",pt,y,c));
            	setTH1(MassPtY[c][pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","#counts");
                MassPtY[c][pt][y]->SetLineColor(1);
                MassPtY[c][pt][y]->SetLineWidth(2);

                MassPtY_More2[c][pt][y] = (TH1F*) in->Get(Form("Mass_More_pt%i_y%i_c%i",pt,y,c));
                MassPtY_More2[c][pt][y]->SetLineColor(2);
                MassPtY_More2[c][pt][y]->SetLineWidth(2);
                
                MassPtY_Less2[c][pt][y] = (TH1F*) in->Get(Form("Mass_Less_pt%i_y%i_c%i",pt,y,c));
                MassPtY_Less2[c][pt][y]->SetLineColor(2);
		MassPtY_Less2[c][pt][y]->SetLineWidth(2);



            }
        }
    }

    // ----> rysowanie

    // add c
  // drawStyle();

  //  gStyle->SetOptTitle(0);

    TCanvas *canPtY[2][8];
  //  TCanvas *canPtYmore[2][8];
  //  TCanvas *canPtYless[2][8];

   //  TCanvas *canPtYmore[8];

    for(Int_t c=0; c<2; c++)
    {

        for(Int_t i=0; i<8; i++)
        {

            canPtY[c][i] = new TCanvas(Form("canc%iPtY%i",c,i),Form("can%iPtY%i",c,i),10,10,1000,1000);
            canPtY[c][i]->Divide(4,3);
            
        //    canPtYmore[c][i] = new TCanvas(Form("canc_more%iPtY%i",c,i),Form("can_more%iPtY%i",c,i),10,10,1000,1000);
        //    canPtYmore[c][i]->Divide(4,3);
            
        //    canPtYless[c][i] = new TCanvas(Form("canc_less%iPtY%i",c,i),Form("can_less%iPtY%i",c,i),10,10,1000,1000);
        //    canPtYless[c][i]->Divide(4,3);

            for(Int_t j=0; j<12; j++)
            {
                canPtY[c][i]->cd(j+1);
                MassPtY[c][j][i]->GetXaxis()->SetRangeUser(0,600);
                MassPtY[c][j][i]->SetLineColor(kBlack);
                MassPtY[c][j][i]->Draw("hist");
                
              //  canPtY[c][i]->cd(j+1);
             //   MassPtY_More2[c][j][i]->GetXaxis()->SetRangeUser(0,600);
                MassPtY_More2[c][j][i]->SetLineColor(kBlue);
                MassPtY_More2[c][j][i]->Draw("same");
                
              //  canPtY[c][i]->cd(j+1);
              //  MassPtY_Less2[c][j][i]->GetXaxis()->SetRangeUser(0,600);
		MassPtY_Less2[c][j][i]->SetLineColor(kRed);
		MassPtY_Less2[c][j][i]->Draw("same");

                setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(200+j*PtBinSize),(200+((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);



            }
        }
    }
   drawStyle();

    gStyle->SetOptTitle(0);

    out->cd();
    out->GetList()->Write();

    for (Int_t a=0; a<2; a++){
        for(Int_t b=0; b<8; b++){
            canPtY[a][b]->Write();
       //     canPtYmore[a][b]->Write();
       //     canPtYless[a][b]->Write();
        }
    }
    out->Save();

}
