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
#include "TText.h"
#include <string>
#include "/u/sudol/unigen/v2.1/drawStyle.C"

using namespace std;

void daraw_sig_stats(TH2F* PtY_sig,Int_t sig_count[8][10],Int_t sig_gaus_peak[8][10],Int_t sig_gaus_sigma[8][10],char can_name[8])
{
    TCanvas * sig_can = new TCanvas(can_name,can_name);
    sig_can->cd();
    PtY_sig->Draw("colztexte");
    TText* te = new TText();
    te->SetTextFont(43);
    te->SetTextSize(18);
    for(Int_t i=0; i<8; i++)
    {
	for(Int_t j=0; j< 10; j++)
	{
            
	    int p =  (sig_gaus_peak[i][j]);
	    int s =  (sig_gaus_sigma[i][j]);

	    char val_p[3];
	    sprintf(val_p,"(%d , %d)",p,s);

	    if ( (p > 130 && p < 160) && sig_count[i][j] >0)
	    {

		te->DrawText(PtY_sig->GetXaxis()->GetBinCenter(j)+60,PtY_sig->GetYaxis()->GetBinCenter(i)+0.12,val_p);
	    }

	}

    }

    delete te;
    sig_can->Write();
}

void plot_comparison2()
{
   drawStyle();

    gStyle->SetOptTitle(0);

    TFile *inMIX = new TFile("diphoton_massEXP_gen2_opAng5_C10_50_MIX.root","read");
    TFile *in    = new TFile("diphoton_massEXP_gen2_opAng5_C10_50.root","read");

    TFile *out    = new TFile("pi0sec_2_4_opAng5_C10_50_sig.root","RECREATE");
    out->cd();

    // ------>  some 2D histograms which are showing phps of gamma e+e-

    //  ------>   histograms for mass spectra in pt-y bins  <--------------- //

    Int_t ptbins = 10;
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


    TH2F *GEPMassPt[2][ptbins][ybins];
    TH2F *GEPMassPt_Less2[2][ptbins][ybins];
    TH2F *GEPMassPt_More2[2][ptbins][ybins];
    TH2F *GEPMassY[2][ptbins][ybins];
    TH2F *GEPMassY_Less2[2][ptbins][ybins];
    TH2F *GEPMassY_More2[2][ptbins][ybins];

////////////////////////////////////////////////
    TH1F *MassPtY_all[ptbins][ybins];
    TH1F *MassPtY_mix[ptbins][ybins];
    TH1F *MassPtY_sig[ptbins][ybins];
    TH2F *PtY_sig =  new TH2F("PtY_sig",";p_{t} [MeV/c];y",ptbins,ptmin,ptmax,ybins,ymin,ymax);
//////////////////////////////////////////////
   // MassPtY_More2

    Int_t binMin;
    Int_t binMax;

     for (Int_t pt = 0; pt < ptbins; pt++)
        {
            for (Int_t y = 0; y < ybins; y++)
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
        for (Int_t pt = 0; pt < ptbins; pt++)
        {
            for (Int_t y = 0; y < ybins; y++)
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

     TCanvas *canPtY_all[ybins];

    for(Int_t i=0; i<ybins; i++)
        {

	    canPtY_all[i] = new TCanvas(Form("cancPtY%i",i),Form("canPtY%i",i),10,10,1000,1000);
	    canPtY_all[i]->Divide(3,3);

	    for(Int_t j=1; j<ptbins; j++)
            {
                canPtY_all[i]->cd(j);
		MassPtY_all[j][i]->GetXaxis()->SetRangeUser(0,600);
		MassPtY_all[j][i]->GetYaxis()->SetRangeUser(-100, MassPtY_all[j][i]->GetMaximum()+100);

		MassPtY_all[j][i]->Draw("");
		setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);

		MassPtY_mix[j][i]->GetXaxis()->SetRangeUser(0,600);
		MassPtY_mix[j][i]->GetYaxis()->SetRangeUser(-100, MassPtY_all[j][i]->GetMaximum()+100);

		MassPtY_mix[j][i]->Draw("same");
		setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);


		MassPtY_sig[j][i]->GetXaxis()->SetRangeUser(0,600);
		MassPtY_sig[j][i]->GetYaxis()->SetRangeUser(-100, MassPtY_all[j][i]->GetMaximum()+100);

		MassPtY_sig[j][i]->Draw("same");
		setOPT_text(Form("#it{p_{t} %.0f-%.0f }",(j*PtBinSize),(((j+1)*PtBinSize))) ,0.5,0.9,0.9,0.999,1,0.2);

            }
	}




    TCanvas *canPtY[2][ybins];

    for(Int_t c=0; c<2; c++)
    {

        for(Int_t i=0; i<ybins; i++)
        {

            canPtY[c][i] = new TCanvas(Form("canc%iPtY%i",c,i),Form("can%iPtY%i",c,i),10,10,1000,1000);
            canPtY[c][i]->Divide(3,3);

            for(Int_t j=1; j<ptbins; j++)
            {
                canPtY[c][i]->cd(j);
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
    Double_t sig_gaus_3sig_min=0;
    Double_t sig_gaus_3sig_max=0;


    Int_t sig_gaus_peak[8][10];
    Int_t sig_gaus_sigma[8][10];
    Int_t sig_gaus_sigcount[8][10];

    TF1* fh = new TF1("fh", "gaus",  100, 200);

    for(Int_t i=0; i<ybins; i++) // rapidity
    {
	for(Int_t j=0; j<ptbins; j++) // pt
	{

	    MassPtY_sig[j][i]->Fit(fh,"RQ");
	    sig_gaus_3sig_min = fh->GetParameter(1) - (3 * fh->GetParameter(2)) ;
	    sig_gaus_3sig_max = fh->GetParameter(1) + (3 * fh->GetParameter(2)) ;

	    sig_peak = MassPtY_sig[j][i]->GetMaximumBin();


	    sig_binMin = MassPtY_sig[j][i]->FindBin(sig_gaus_3sig_min);
	    sig_binMax = MassPtY_sig[j][i]->FindBin(sig_gaus_3sig_max);
            sig_gaus_peak[i][j]  = fh->GetParameter(1)   ;
	    sig_gaus_sigma[i][j] = fh->GetParameter(2)   ;


	    if (  fh->GetParameter(1) > 130 &&  fh->GetParameter(1) < 160 )
	    {
		sig_count =   MassPtY_sig[j][i]->Integral(sig_binMin,sig_binMax);

	    }
	    else sig_count = 0;
            sig_gaus_sigcount[i][j] = sig_count;

	    if (sig_count>0)
	    {
		PtY_sig->SetBinContent(j+1,i+1,sig_count);
	    }
	    else             PtY_sig->SetBinContent(j+1,i+1,0);
	}
    }


    char sig_can_name[ybins] = "sig_all";
    daraw_sig_stats(PtY_sig,sig_gaus_sigcount,sig_gaus_peak,sig_gaus_sigma,sig_can_name);

   /*
    Double_t sig_more_binMin=0.;
    Double_t sig_more_binMax=0.;
    Double_t sig_more_peak=0;

    Double_t sig_less_binMin=0.;
    Double_t sig_less_binMax=0.;
    Double_t sig_less_peak=0;

    Double_t sig_gaus_3sig_min_less =0;
    Double_t sig_gaus_3sig_max_less =0;
    Double_t sig_gaus_3sig_min_more =0;
    Double_t sig_gaus_3sig_max_more =0;      */
    Double_t sig_more_count=0.;
    Double_t sig_less_count=0.;
    Int_t sig_gaus_peak_more[2][8][10];
    Int_t sig_gaus_sigma_more[2][8][10];
    Int_t sig_gaus_peak_less[2][8][10];
    Int_t sig_gaus_sigma_less[2][8][10];
    Int_t sig_gaus_sigma_sigcount_less[2][8][10];
    Int_t sig_gaus_sigma_sigcount_more[2][8][10];

    for (Int_t c=0; c<2; c++)
    {
	for(Int_t i=0; i<ybins; i++) // rapidity
	{

	    for(Int_t j=0; j<ptbins; j++) // pt
	    {

		MassPtY_More2_sig[c][j][i]->Fit(fh,"RQ");
		sig_gaus_3sig_min = fh->GetParameter(1) - (3 * fh->GetParameter(2)) ;
		sig_gaus_3sig_max = fh->GetParameter(1) + (3 * fh->GetParameter(2)) ;
		sig_binMin = MassPtY_More2_sig[c][j][i]->FindBin(sig_gaus_3sig_min);
		sig_binMax = MassPtY_More2_sig[c][j][i]->FindBin(sig_gaus_3sig_max);
		sig_gaus_peak_more[c][i][j]  = fh->GetParameter(1)   ;
		sig_gaus_sigma_more[c][i][j] = fh->GetParameter(2)   ;
		if (  fh->GetParameter(1) > 130 &&  fh->GetParameter(1) < 160 )
		{
		    sig_more_count =   MassPtY_More2_sig[c][j][i]->Integral(sig_binMin,sig_binMax);
		}
                else sig_more_count = 0;
		sig_gaus_sigma_sigcount_more[c][i][j] = sig_more_count;
	        if (sig_more_count>0)
		{
		    PtY_more_sig[c]->SetBinContent(j+1,i+1,sig_more_count);
		}
		else             PtY_more_sig[c]->SetBinContent(j+1,i+1,0);





		MassPtY_Less2_sig[c][j][i]->Fit(fh,"RQ");
		sig_gaus_3sig_min = fh->GetParameter(1) - (3 * fh->GetParameter(2)) ;
		sig_gaus_3sig_max = fh->GetParameter(1) + (3 * fh->GetParameter(2)) ;
		sig_binMin =  MassPtY_Less2_sig[c][j][i]->FindBin(sig_gaus_3sig_min);
		sig_binMax =  MassPtY_Less2_sig[c][j][i]->FindBin(sig_gaus_3sig_max);
		sig_gaus_peak_less[c][i][j]  = fh->GetParameter(1)   ;
		sig_gaus_sigma_less[c][i][j] = fh->GetParameter(2)   ;
		if (  fh->GetParameter(1) > 130 &&  fh->GetParameter(1) < 160 )
		{

		    sig_less_count =   MassPtY_Less2_sig[c][j][i]->Integral(sig_binMin,sig_binMax);

		}
		else sig_less_count = 0;
		sig_gaus_sigma_sigcount_less[c][i][j] = sig_less_count;
	        if (sig_less_count>0)
		{
		    PtY_less_sig[c]->SetBinContent(j+1,i+1,sig_less_count);
		}
		else             PtY_less_sig[c]->SetBinContent(j+1,i+1,0);
	    }
	}
    }

    char sig_can_name_more[ybins] = "sig_mor";
    daraw_sig_stats(PtY_more_sig[0],sig_gaus_sigma_sigcount_more[0],sig_gaus_peak_more[0],sig_gaus_sigma_more[0],sig_can_name_more);
    char sig_can_name_less[ybins] = "sig_les";
    daraw_sig_stats(PtY_less_sig[0],sig_gaus_sigma_sigcount_less[0],sig_gaus_peak_less[0],sig_gaus_sigma_less[0],sig_can_name_less);

    TH1F * proj_sig;
    TCanvas* can_proj_sig= new TCanvas("pi_sig_proj","pi_sig_proj") ;

    can_proj_sig->cd();
 //   proj_sig = (TH1F*)PtY_sig->ProjectionX(Form("Bin0",0),0,0+1);
 //   proj_sig->Draw("APL");
    const EColor colours[] = {kBlue,kRed,kGreen,kBlack,kTeal,kOrange,kPink,kYellow};

    TLegend* leg = new TLegend ( 0.5,0.7,0.9,0.9 );
    leg->SetHeader ( "Rapidity Y" );
    leg->SetFillColor ( 1 );

    Double_t Y_firstbin = 0.4;

    for (Int_t a=0; a<ybins; a++)
    {
	proj_sig = (TH1F*)PtY_sig->ProjectionX(Form("Bin%i",a),a,a+1);
        proj_sig->GetXaxis()->SetRangeUser(0,8000);
	proj_sig->Draw("PL,same");
	proj_sig->SetLineColor(colours[a]);
	proj_sig->SetLineWidth(2);
	leg->AddEntry (proj_sig,Form("Y %.2f-%.2f",Y_firstbin,Y_firstbin+0.2),"lep" );
        Y_firstbin=Y_firstbin+0.2;

    }
    leg->SetFillStyle ( 0 );
    leg->Draw();
  //  can_proj_sig->GetXaxis()->SetRangeUser(0,8000);

    gStyle->SetPaintTextFormat("4.0f");


    out->cd();
    can_proj_sig->Write();
    out->GetList()->Write();

    for (Int_t a=0; a<2; a++){
        for(Int_t b=0; b<ybins; b++){
            canPtY[a][b]->Write();
           // canPtY_all[a][b]->Write();
       //     canPtYless[a][b]->Write();
        }
    }
    for (Int_t a =0; a<ybins; a++)
    {
	canPtY_all[a]->Write();
    }

    out->Save();

}
