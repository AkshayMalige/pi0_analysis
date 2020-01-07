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

    TFile *inMIX = new TFile("diphoton_massEXP_gen2_opAng5_C10_50_MIX.root","read");
    TFile *in    = new TFile("diphoton_massEXP_gen2_opAng5_C10_50.root","read");

    TFile *out    = new TFile("pi0sec_2_4_opAng5_C10_50_sig.root","RECREATE");
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

//////////////////////////////////////////////
    TH1F *MassPtY[2][20][8];  // MassPtY[2][20][8]->Sum2w();

    TH1F *MassPtY_Less2_all[2][20][8];// MassPtY_Less2_all[2][20][8]->Sum2w();
    TH1F *MassPtY_More2_all[2][20][8]; //MassPtY_More2_all[2][20][8]->Sum2w();

    TH1F *MassPtY_Less2_mix[2][20][8]; //MassPtY_Less2_mix[2][20][8]->Sum2w();
    TH1F *MassPtY_More2_mix[2][20][8]; //MassPtY_More2_mix[2][20][8]->Sum2w();

    TH1F *MassPtY_Less2_sig[2][20][8];// MassPtY_Less2_sig[2][20][8]->Sum2w();
    TH1F *MassPtY_More2_sig[2][20][8]; //MassPtY_More2_sig[2][20][8]->Sum2w();

    TH2F *PtY_more_sig[2];//  PtY_more_sig[2]->Sumw2();
    TH2F *PtY_less_sig[2]; // PtY_less_sig[2]->Sumw2();

    for (Int_t t=0; t<2; t++)
    {

	PtY_more_sig[t] =  new TH2F ( Form("PtY_more_sig_cent%i",t) , Form("PtY_more_sig_cent%i;p_{t} [MeV/c];y",t),ptbins,ptmin,ptmax,ybins,ymin,ymax);
	PtY_less_sig[t] =  new TH2F ( Form("PtY_less_sig_cent%i",t) , Form("PtY_less_sig_cent%i;p_{t} [MeV/c];y",t),ptbins,ptmin,ptmax,ybins,ymin,ymax);
     //   PtY_more_sig[t]->Sumw2();
      //  PtY_less_sig[t]->Sumw2();

    }

///////////////////////////////////////////////


    TH2F *GEPMassPt[2][20][8];
    TH2F *GEPMassPt_Less2[2][20][8];
    TH2F *GEPMassPt_More2[2][20][8];
    TH2F *GEPMassY[2][20][8];
    TH2F *GEPMassY_Less2[2][20][8];
    TH2F *GEPMassY_More2[2][20][8];

////////////////////////////////////////////////
    TH1F *MassPtY_all[20][8];
    TH1F *MassPtY_mix[20][8];
    TH1F *MassPtY_sig[20][8];
    TH2F *PtY_sig =  new TH2F("PtY_sig",";p_{t} [MeV/c];y",ptbins,ptmin,ptmax,ybins,ymin,ymax);
//////////////////////////////////////////////
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
              //  MassPtY_mix[pt][y]->Sumw2();

		MassPtY_all[pt][y]   = (TH1F*) in->Get(Form("Mass_pt_%i_y%i_all",pt,y));
		setTH1(MassPtY_all[pt][y],"M_{#gamma e^{+}e^{-}} [MeV/c^{2}]","#counts");
		MassPtY_all[pt][y]->SetLineColor(1);
		MassPtY_all[pt][y]->SetLineWidth(2);
              //  MassPtY_all[pt][y]->Sumw2();

                binMin = MassPtY_all[pt][y]->FindBin(250);
                binMax = MassPtY_all[pt][y]->FindBin(600);

                Float_t int_all = MassPtY_all[pt][y]->Integral(binMin,binMax);
		Float_t int_mix = MassPtY_mix[pt][y]->Integral(binMin,binMax);

		MassPtY_mix[pt][y]->Scale(int_all/int_mix);


		MassPtY_sig[pt][y] = (TH1F*)MassPtY_all[pt][y]->Clone();
		MassPtY_sig[pt][y]->SetNameTitle(Form("MassSig_pt_%i_y%i_all",pt,y),Form("MassSig_pt_%i_y%i_all",pt,y));
		MassPtY_sig[pt][y]->Add(MassPtY_mix[pt][y],-1);
		MassPtY_sig[pt][y]->SetLineColor(2);
		MassPtY_sig[pt][y]->SetLineWidth(2);
                MassPtY_sig[pt][y]->SetMarkerColor(kRed);
                MassPtY_sig[pt][y]->Sumw2();


	    }
	}

     Int_t binMinmore;
     Int_t binMaxmore;
     Int_t binMinless;
     Int_t binMaxless;

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
            //    MassPtY[c][pt][y]->Sumw2();


                MassPtY_More2_all[c][pt][y] = (TH1F*) in->Get(Form("Mass_More_pt%i_y%i_c%i",pt,y,c));
                MassPtY_More2_all[c][pt][y]->SetLineColor(2);
                MassPtY_More2_all[c][pt][y]->SetLineWidth(2);
                MassPtY_Less2_all[c][pt][y] = (TH1F*) in->Get(Form("Mass_Less_pt%i_y%i_c%i",pt,y,c));
                MassPtY_Less2_all[c][pt][y]->SetLineColor(3);
		MassPtY_Less2_all[c][pt][y]->SetLineWidth(2);
             //   MassPtY_Less2_all[c][pt][y]->Sumw2();

		MassPtY_More2_mix[c][pt][y] = (TH1F*) inMIX->Get(Form("Mass_More_pt%i_y%i_c%i",pt,y,c));
                MassPtY_More2_mix[c][pt][y]->SetLineColor(2);
                MassPtY_More2_mix[c][pt][y]->SetLineWidth(2);
                MassPtY_Less2_mix[c][pt][y] = (TH1F*) inMIX->Get(Form("Mass_Less_pt%i_y%i_c%i",pt,y,c));
                MassPtY_Less2_mix[c][pt][y]->SetLineColor(3);
		MassPtY_Less2_mix[c][pt][y]->SetLineWidth(2);
             //   MassPtY_Less2_mix[c][pt][y]->Sumw2();

                binMinmore =  MassPtY_More2_all[c][pt][y]->FindBin(250);
                binMaxmore =  MassPtY_More2_all[c][pt][y]->FindBin(600);
                binMinless =  MassPtY_Less2_all[c][pt][y]->FindBin(250);
		binMinless =  MassPtY_Less2_all[c][pt][y]->FindBin(600);

		Float_t int_more_all = MassPtY_More2_all[c][pt][y]->Integral(binMinmore,binMaxmore);
		Float_t int_more_mix = MassPtY_More2_mix[c][pt][y]->Integral(binMinmore,binMaxmore);
		Float_t int_less_all = MassPtY_Less2_all[c][pt][y]->Integral(binMinless,binMaxless);
		Float_t int_less_mix = MassPtY_Less2_mix[c][pt][y]->Integral(binMinless,binMaxless);

		MassPtY_More2_mix[c][pt][y]->Scale(int_more_all/int_more_mix);
		MassPtY_Less2_mix[c][pt][y]->Scale(int_less_all/int_less_mix);

		MassPtY_More2_sig[c][pt][y] = (TH1F*)MassPtY_More2_all[c][pt][y]->Clone();
		MassPtY_More2_sig[c][pt][y]->SetNameTitle(Form("MassSig_more_pt_%i_y%i_all",pt,y),Form("MassSig_more_pt_%i_y%i_all",pt,y));
		MassPtY_More2_sig[c][pt][y]->Add(MassPtY_More2_mix[c][pt][y],-1);
		MassPtY_More2_sig[c][pt][y]->SetLineColor(2);
		MassPtY_More2_sig[c][pt][y]->SetLineWidth(2);
	        MassPtY_More2_sig[c][pt][y]->Sumw2();
 		MassPtY_Less2_sig[c][pt][y] = (TH1F*)MassPtY_Less2_all[c][pt][y]->Clone();
		MassPtY_Less2_sig[c][pt][y]->SetNameTitle(Form("MassSig_less_pt_%i_y%i_all",pt,y),Form("MassSig_less_pt_%i_y%i_all",pt,y));
		MassPtY_Less2_sig[c][pt][y]->Add(MassPtY_Less2_mix[c][pt][y],-1);
		MassPtY_Less2_sig[c][pt][y]->SetLineColor(2);
		MassPtY_Less2_sig[c][pt][y]->SetLineWidth(2);
                MassPtY_Less2_sig[c][pt][y]->Sumw2();


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


    TCanvas *canPtY[2][8];

    for(Int_t c=0; c<2; c++)
    {

        for(Int_t i=0; i<8; i++)
        {

            canPtY[c][i] = new TCanvas(Form("canc%iPtY%i",c,i),Form("can%iPtY%i",c,i),10,10,1000,1000);
            canPtY[c][i]->Divide(5,4);

            for(Int_t j=0; j<20; j++)
            {
                canPtY[c][i]->cd(j+1);
                MassPtY[c][j][i]->GetXaxis()->SetRangeUser(0,600);
                MassPtY[c][j][i]->SetLineColor(kBlack);
                MassPtY[c][j][i]->Draw("hist");

                MassPtY_More2_all[c][j][i]->Draw("same");

		MassPtY_Less2_all[c][j][i]->Draw("same");

                setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.7,0.9,0.999,1,0.2);

            }
        }
    }


    Double_t sig_binMin=0.;
    Double_t sig_binMax=0.;
    Double_t sig_count=0.;
    Double_t sig_peak=0;

    TCanvas *canPtY_all[8];
    for(Int_t i=0; i<8; i++) // rapidity
    {

	canPtY_all[i] = new TCanvas(Form("cancPtY%i",i),Form("canPtY%i",i),10,10,1000,1000);
	canPtY_all[i]->Divide(5,4);
	for(Int_t j=0; j<20; j++) // pt
	{

	    canPtY_all[i]->cd(j+1);
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

	    sig_count =   MassPtY_sig[j][i]->Integral(sig_peak-3,sig_peak+3);

	    if (sig_count>0)
	    {
		PtY_sig->SetBinContent(j+1,i+1,sig_count);
	    }
	    else             PtY_sig->SetBinContent(j+1,i+1,0);
	}
    }

    Double_t sig_more_binMin=0.;
    Double_t sig_more_binMax=0.;
    Double_t sig_more_count=0.;
    Double_t sig_more_peak=0;

    Double_t sig_less_binMin=0.;
    Double_t sig_less_binMax=0.;
    Double_t sig_less_count=0.;
    Double_t sig_less_peak=0;


    for (Int_t c=0; c<2; c++)
    {
	for(Int_t i=0; i<8; i++) // rapidity
	{

	    for(Int_t j=0; j<20; j++) // pt
	    {

		sig_more_peak = MassPtY_More2_sig[c][j][i]->GetMaximumBin();
		sig_more_binMin = MassPtY_More2_sig[c][j][i]->FindBin(100);
		sig_more_binMax = MassPtY_More2_sig[c][j][i]->FindBin(200);
		sig_more_count =   MassPtY_More2_sig[c][j][i]->Integral(sig_more_peak-3,sig_more_peak+3);

	        if (sig_more_count>0)
		{
		    PtY_more_sig[c]->SetBinContent(j+1,i+1,sig_more_count);
		}
		else             PtY_more_sig[c]->SetBinContent(j+1,i+1,0);


		sig_less_peak = MassPtY_Less2_sig[c][j][i]->GetMaximumBin();
		sig_less_binMin = MassPtY_Less2_sig[c][j][i]->FindBin(100);
		sig_less_binMax = MassPtY_Less2_sig[c][j][i]->FindBin(200);
		sig_less_count =   MassPtY_Less2_sig[c][j][i]->Integral(sig_less_peak-3,sig_less_peak+3);

	        if (sig_less_count>0)
		{
		    PtY_less_sig[c]->SetBinContent(j+1,i+1,sig_less_count);
		}
		else             PtY_less_sig[c]->SetBinContent(j+1,i+1,0);
	    }
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
