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

    TFile *inMIX = new TFile("diphoton_massEXP_gen2_opAng5_MIX.root","read");
    TFile *in    = new TFile("diphoton_massEXP_gen2_opAng5.root","read");

    TFile *out    = new TFile("pi0sec_2_4comparision_MIX.root","RECREATE");
    out->cd();

    // ------>  some 2D histograms which are showing phps of gamma e+e-

    //  ------>   histograms for mass spectra in pt-y bins  <--------------- //

    Int_t ptbins = 20;
    Int_t ptmin  = 0;
    Int_t ptmax  = 1000;

    Int_t ybins    = 8;
    Double_t ymin  = 0.4;
    Double_t ymax  = 2.;
    
    Double_t PtBinSize = (ptmax-ptmin)/ptbins; // MeV/bin
    Double_t YBinSize  = (ymax-ymin)/ybins;  // rapidity/bin

    TH1F *MassPtY[2][20][8];
    TH1F *MassPtY_Less2[2][20][8];
    TH1F *MassPtY_More2[2][20][8];

    TH2F *GEPMassPt[2][20][8];
    TH2F *GEPMassPt_Less2[2][20][8];
    TH2F *GEPMassPt_More2[2][20][8];

    TH2F *GEPMassY[2][20][8];
    TH2F *GEPMassY_Less2[2][20][8];
    TH2F *GEPMassY_More2[2][20][8];

    TH1F *MassPtY_all[20][8];
    TH1F *MassPtY_mix[20][8];
    TH1F *MassPtY_sig[20][8];

    TH2F *PtY_sig =  new TH2F("PtY_sig",";p_{t} [MeV/c];y",ptbins,ptmin,ptmax,ybins,ymin,ymax);

   // MassPtY_More2

    Int_t binMin;
    Int_t binMax;

     for (Int_t pt = 0; pt < 20; pt++)
        {
            for (Int_t y = 0; y < 8; y++)
	    {
		MassPtY_mix[pt][y]   = (TH1F*) inMIX->Get(Form("Mass_pt_%i_y%i_all",pt,y));
		setTH1(MassPtY_mix[pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","#counts");
		MassPtY_mix[pt][y]->SetLineColor(4);
		MassPtY_mix[pt][y]->SetLineWidth(2);
		MassPtY_mix[pt][y]->SetNameTitle(Form("MassMIX_pt_%i_y%i_all",pt,y),Form("MassMIX_pt_%i_y%i_all",pt,y));

		MassPtY_all[pt][y]   = (TH1F*) in->Get(Form("Mass_pt_%i_y%i_all",pt,y));
		setTH1(MassPtY_all[pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","#counts");
		MassPtY_all[pt][y]->SetLineColor(1);
		MassPtY_all[pt][y]->SetLineWidth(2);

                binMin = MassPtY_all[pt][y]->FindBin(250);
                binMax = MassPtY_all[pt][y]->FindBin(600);

                Float_t int_all = MassPtY_all[pt][y]->Integral(binMin,binMax);
		Float_t int_mix = MassPtY_mix[pt][y]->Integral(binMin,binMax);

	       // cout<<"int_all : "<<int_all<<"  int_mix  "<<int_mix<<endl;

		MassPtY_mix[pt][y]->Scale(int_all/int_mix);

		MassPtY_sig[pt][y] = (TH1F*)MassPtY_all[pt][y]->Clone();
		MassPtY_sig[pt][y]->SetNameTitle(Form("MassSig_pt_%i_y%i_all",pt,y),Form("MassSig_pt_%i_y%i_all",pt,y));
		MassPtY_sig[pt][y]->Add(MassPtY_mix[pt][y],-1);
		MassPtY_sig[pt][y]->SetLineColor(2);
		MassPtY_sig[pt][y]->SetLineWidth(2);


	    }
	}

    for (Int_t c = 0; c < 2; c++)
    {
        for (Int_t pt = 0; pt < 20; pt++)
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
                MassPtY_Less2[c][pt][y]->SetLineColor(3);
		MassPtY_Less2[c][pt][y]->SetLineWidth(2);


                GEPMassPt[c][pt][y]   = (TH2F*) in->Get(Form("MassGEP_Pt_pt%i_y%i_c%i",pt,y,c));
            	setTH2(GEPMassPt[c][pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","Pt");
                GEPMassPt_More2[c][pt][y] = (TH2F*) in->Get(Form("MoreMassGEP_Pt_pt%i_y%i_c%i",pt,y,c));
		GEPMassPt_Less2[c][pt][y] = (TH2F*) in->Get(Form("LessMassGEP_Pt_pt%i_y%i_c%i",pt,y,c));

                GEPMassY[c][pt][y]   = (TH2F*) in->Get(Form("MassGEP_Y_pt%i_y%i_c%i",pt,y,c));
            	setTH2(GEPMassY[c][pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","Y");
                GEPMassY_More2[c][pt][y] = (TH2F*) in->Get(Form("MoreMassGEP_Y_pt%i_y%i_c%i",pt,y,c));
                GEPMassY_Less2[c][pt][y] = (TH2F*) in->Get(Form("LessMassGEP_Y_pt%i_y%i_c%i",pt,y,c));




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
            canPtY[c][i]->Divide(5,4);
            
        //    canPtYmore[c][i] = new TCanvas(Form("canc_more%iPtY%i",c,i),Form("can_more%iPtY%i",c,i),10,10,1000,1000);
        //    canPtYmore[c][i]->Divide(4,3);
            
        //    canPtYless[c][i] = new TCanvas(Form("canc_less%iPtY%i",c,i),Form("can_less%iPtY%i",c,i),10,10,1000,1000);
        //    canPtYless[c][i]->Divide(4,3);

            for(Int_t j=0; j<20; j++)
            {
                canPtY[c][i]->cd(j+1);
                MassPtY[c][j][i]->GetXaxis()->SetRangeUser(0,600);
                MassPtY[c][j][i]->SetLineColor(kBlack);
                MassPtY[c][j][i]->Draw("hist");
                
              //  canPtY[c][i]->cd(j+1);
             //   MassPtY_More2[c][j][i]->GetXaxis()->SetRangeUser(0,600);
                MassPtY_More2[c][j][i]->Draw("same");
                
              //  canPtY[c][i]->cd(j+1);
              //  MassPtY_Less2[c][j][i]->GetXaxis()->SetRangeUser(0,600);
		MassPtY_Less2[c][j][i]->Draw("same");

                setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.7,0.9,0.999,1,0.2);



            }
        }
    }
    Double_t sig_binMin=0.;
    Double_t sig_binMax=0.;
    Double_t sig_count=0.;
    Double_t sig_peak=0;
//   MassPtY_all[j][i]->Scale(1./2);
//   MassPtY_mix[j][i]->Rebin(2);
    TCanvas *canPtY_all[8];
    for(Int_t i=0; i<8; i++) // rapidity
    {

	canPtY_all[i] = new TCanvas(Form("cancPtY%i",i),Form("canPtY%i",i),10,10,1000,1000);
	canPtY_all[i]->Divide(5,4);

	//    canPtYmore[c][i] = new TCanvas(Form("canc_more%iPtY%i",c,i),Form("can_more%iPtY%i",c,i),10,10,1000,1000);
	//    canPtYmore[c][i]->Divide(4,3);

	//    canPtYless[c][i] = new TCanvas(Form("canc_less%iPtY%i",c,i),Form("can_less%iPtY%i",c,i),10,10,1000,1000);
	//    canPtYless[c][i]->Divide(4,3);

	for(Int_t j=0; j<20; j++) // pt
	{

	    canPtY_all[i]->cd(j+1);
 
	   // canPtY_all[i]->cd(j+1);
	    MassPtY_all[j][i]->GetXaxis()->SetRangeUser(0,600);
	    MassPtY_all[j][i]->Draw("");
	    setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);

	    MassPtY_mix[j][i]->GetXaxis()->SetRangeUser(0,600);
	    MassPtY_mix[j][i]->Draw("same");
	    setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);

	//    MassPtY_all[j][i]->Rebin(2);
	//    MassPtY_mix[j][i]->Rebin(2);


	    MassPtY_sig[j][i]->GetXaxis()->SetRangeUser(0,600);
	    MassPtY_sig[j][i]->Draw("same");
	    setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);

            sig_peak = MassPtY_sig[j][i]->GetMaximumBin();

	    sig_binMin = MassPtY_sig[j][i]->FindBin(100);
	    sig_binMax = MassPtY_sig[j][i]->FindBin(200);

	    // cout<<j<<"\t"<<i<<endl ;
         //   gStyle->SetPaintTextFormat("4.0f");
	    sig_count =   MassPtY_sig[j][i]->Integral(sig_peak-3,sig_peak+3);

	    if (sig_count>0)
	    {
	       // cout<<"j :"<<j*PtBinSize<<"\ti :"<<0.4+(i*YBinSize)<<"\tintg :"<< MassPtY_sig[j][i]->Integral(sig_binMin,sig_binMax)<<endl;

		PtY_sig->SetBinContent(j+1,i+1,sig_count);

                cout<<PtY_sig->GetBinContent(j+1,i+1)<<endl;
	    }
	    else             PtY_sig->SetBinContent(j+1,i+1,0);
	   // PtY_sig->SetMarkerColor(kRed);
	  //  PtY_sig->Draw("COLZTEXT");
//	    PtY_sig->Draw("TEXT SAME");
//	    if (sig_count>0) PtY_sig->Fill(j*PtBinSize,0.4+(i*YBinSize),sig_count);
//	    else             PtY_sig->Fill(j*PtBinSize,0.4+(i*YBinSize),0);
	}
    }

  // drawStyle();

   // gStyle->SetOptTitle(0);

    out->cd();
    out->GetList()->Write();

    for (Int_t a=0; a<2; a++){
        for(Int_t b=0; b<8; b++){
            canPtY[a][b]->Write();
           // canPtY_all[a][b]->Write();
       //     canPtYless[a][b]->Write();
        }
    }
    for (Int_t a =0; a<8; a++)
    {
	canPtY_all[a]->Write();
    }

    out->Save();

}
