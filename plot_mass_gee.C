#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>
#include <map>
#include <stdio.h>

#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"


typedef struct {Float_t mult_tofrpc,mult_emc,sectorP,betaP,thetaP,phiP,chargeP,richQaP,momP,
sectorE,betaE,thetaE,phiE,chargeE,richQaE,momE,
sectorG,betaG,cellG,energyG,timeG,sizeG,thetaG,phiG,matchCellG,matchPartG,
massLL,angleLL,ptLL,yLL,thetaLL,phiLL,momLL,
massGEP,angleGEP,ptGEP,yGEP,thetaGEP,phiGEP,momGEP;} PID;


const Double_t d2r = 0.01745329;
void plot_mass_gee(TString inFile, TString outputFile, Int_t nev=1000., Bool_t useMix=kFALSE)

    //void plot_mass_gee(TString inFile, TString outputFile, Int_t nev=100000000., Bool_t useMix=kFALSE)
{

    //data

    TChain nt = TChain("nt");
    nt.AddFile(inFile);

    PID pid;
    nt.SetBranchAddress("pid" , &pid.mult_tofrpc);
    Long64_t nentries        = nt.GetEntries();

    cout<<"nentries:  "<<nentries<<endl;

    TFile *out    = new TFile("diphoton_massEXP_gen2_opAng5_C10_50_MIX.root","RECREATE");
    out->cd();

    // histogram definition

    TH1F* hAll_BetaGamma    = new TH1F("hAll_BetaGamma",";#beta #gamma",500,0.8,1.3);
    TH1F* hAll_massGEP    = new TH1F("hAll_massGEP",";M_{#gamma e^{+}e^{-} [MeV/c]}",1000,0,1000);
    TH1F* hAll_EEAng    = new TH1F("hAll_EEAng",";#alpha _{e^{+}e^{-}}",150,0,150);
    TH2F* hAll_GEPmassAngll = new TH2F("hAll_GEPmassAngll",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];#alpha _{e^{+}e^{-}}",1000,0,1000.0,150,0,150);
    TH2F* hAll_GEPmassMultmulttofrpc = new TH2F("hAll_GEPmassMultmulttofrpc",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];Mult^{TOF+RPC}",1000,0,1000.0,160,0,160);
    TH2F* hAll_GEPmassMultmulttofrpc0 = new TH2F("hAll_GEPmassMultmulttofrpc0",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];Mult^{TOF+RPC}",1000,0,1000.0,160,0,160);

    TH2F* hAll_GEPmassyGEP = new TH2F("hAll_GEPmassyGEP",";M_{#gamma e^{+}e^{-}} [MeV/c^];y",1000,0,1000.0,300,0,3);
    TH2F* hAll_GEPmassptGEP = new TH2F("hAll_GEPmassptGEP",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];p_{t}",1000,0,1000.0,1000,0,1000);
    TH2F* hAll_GEPmassAng_en120_mult = new TH2F("hAll_GEPmassAng_en120_mult",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];#alpha _{#gamma e^{+} e^{-}}",300,0,600.0,150,0,150);
    TH2F* hAll_ptGEPyGEP = new TH2F("hAll_ptGEPyGEP",";p_{t} [MeV/c];y",20,0,1000.0,8,0.4,2);
    TH2F *beta = new TH2F("beta","beta;pxq/|q| (MeV/c);#beta",400,-1000,1000,100,0.8,1.4);

    TH2F* hAll_AngGEP_AngLL = new TH2F("hAll_AngGEP_AngLL",";#alpha _{#gamma e^{+} e^{-}};#alpha _{e^{+} e^{-}}",140,0,140,140,0,140);

    TH2F*  hAll_GEPmass_en_mult[5];
    TH2F*  hAll_ptGEP_GEPmass[5];
    TH2F*  hAll_yGEP_GEPmass[5];
    // Int_t  angl_limit = 2;
    for (int a =0; a<5; a++)
    {
	hAll_GEPmass_en_mult[a] = new TH2F ( Form ( "hAll_GEPmass_en_mult_%d", a ) , Form ( "hAll_GEPmass_en_mult_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];Mult^{TOF+RPC}", a ), 1000, 0,1000,160,0,160 );
	hAll_ptGEP_GEPmass[a] = new TH2F ( Form ( "hAll_ptGEP_GEPmass_%d", a ) , Form ( "hAll_ptGEP_GEPmass_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];p_{t}", a ), 1000, 0,1000,16,0,800 );
	hAll_yGEP_GEPmass[a] = new TH2F ( Form ( "hAll_yGEP_GEPmass_%d", a ) , Form ( "hAll_yGEP_GEPmass_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];y", a ), 1000, 0,1000,10,0,2 );
	//   angl_limit = angl_limit+2;


    }
    TH2F*  hAll_GEPmass_en_mult_noOp = new TH2F ("hAll_GEPmass_en_mult_noOp","hAll_GEPmass_en_mult_noOp;M_{#gamma e^{+}e^{-}} [MeV/c^];Mult^{TOF+RPC}",1000, 0,1000,160,0,160);
    TH2F*  hAll_ptGEP_GEPmass_noOp =  new TH2F ( "hAll_ptGEP_GEPmass_noOp" , "hAll_ptGEP_GEPmass_noOp;M_{#gamma e^{+}e^{-}} [MeV/c^];p_{t}", 1000, 0,1000,16,0,800 );
    TH2F*  hAll_yGEP_GEPmass_noOp   = new TH2F ( "hAll_yGEP_GEPmass_noOp" , "hAll_yGEP_GEPmass_noOp;M_{#gamma e^{+}e^{-}} [MeV/c^];y", 1000, 0,1000,10,0,2 );



    Int_t ptbins = 20.; // for pi0 reduced
    Int_t ptmin  = 0;
    Int_t ptmax  = 1000;

    //Int_t ybins    = 15.; // for eta
    Int_t ybins    = 8.; // for pi0
    Double_t ymin  = 0.4;  //0.015
    Double_t ymax  = 2.;  //2.015
    Double_t deltay = 0.2;

    Double_t PtBinSize = (ptmax-ptmin)/ptbins; // MeV/bin
    Double_t YBinSize  = (ymax-ymin)/ybins;  // rapidity/bin

    Int_t massbins = 60;
    Int_t massmin  = 0;
    Int_t massmax  = 600;

    TH2F*  hAll_ptGEP_GEPmass_Less2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_ptGEP_GEPmass_More2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_yGEP_GEPmass_Less2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_yGEP_GEPmass_More2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_ptGEP_yGEP_Less2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_ptGEP_yGEP_More2Deg[2];// [0] : 30-50%, [1] : 10-30%

    TH2F*  hAll_massGEP_yGEP_Less2Deg[2];// [0] : 30-50%, [1] : 10-30%
    TH2F*  hAll_massGEP_yGEP_More2Deg[2];// [0] : 30-50%, [1] : 10-30%

    for(Int_t i=0; i<2; i++){
	hAll_ptGEP_GEPmass_Less2Deg[i] =  new TH2F ( Form("hAll_ptGEP_GEPmass_Less2Deg_cent%i",i) , Form("hAll_ptGEP_GEPmass_Less2Deg_cent%i;M_{#gamma e^{+}e^{-}} [MeV/c^];p_{t}",i),massbins,massmin,massmax,ptbins,ptmin,ptmax);
	hAll_ptGEP_GEPmass_More2Deg[i] =  new TH2F ( Form("hAll_ptGEP_GEPmass_More2Deg_cent%i",i) , Form("hAll_ptGEP_GEPmass_More2Deg_cent%i;M_{#gamma e^{+}e^{-}} [MeV/c^];p_{t}",i),massbins,massmin,massmax,ptbins,ptmin,ptmax);
	hAll_yGEP_GEPmass_Less2Deg[i] =  new TH2F ( Form("hAll_yGEP_GEPmass_Less2Deg_cent%i",i) , Form("hAll_yGEP_GEPmass_Less2Deg_cent%i;M_{#gamma e^{+}e^{-}} [MeV/c^];y",i),massbins,massmin,massmax,ybins,ymin,ymax);
	hAll_yGEP_GEPmass_More2Deg[i] =  new TH2F ( Form("hAll_yGEP_GEPmass_More2Deg_cent%i",i) , Form("hAll_yGEP_GEPmass_More2Deg_cent%i;M_{#gamma e^{+}e^{-}} [MeV/c^];y",i),massbins,massmin,massmax,ybins,ymin,ymax);
	hAll_ptGEP_yGEP_Less2Deg[i] =  new TH2F ( Form("hAll_ptGEP_yGEP_Less2Deg_cent%i",i) , Form("hAll_ptGEP_yGEP_Less2Deg_cent%i;p_{t};y",i),ptbins,ptmin,ptmax,ybins,ymin,ymax);
	hAll_ptGEP_yGEP_More2Deg[i] =  new TH2F ( Form("hAll_ptGEP_yGEP_More2Deg_cent%i",i) , Form("hAll_ptGEP_yGEP_More2Deg_cent%i;p_{t};y",i),ptbins,ptmin,ptmax,ybins,ymin,ymax);
    }


    Int_t YBin=0;
    Int_t PtBin=0;
    TString HistName;
    TH1F *MassPtY[2][20][13];
    TH1F *MassPtY_More2[2][20][13];
    TH1F *MassPtY_Less2[2][20][13];

    TH2F *GEPMass_Pt[2][20][13];
    TH2F *GEPMass_PtMore2[2][20][13];
    TH2F *GEPMass_PtLess2[2][20][13];

    TH2F *GEPMass_Y[2][20][13];
    TH2F *GEPMass_YMore2[2][20][13];
    TH2F *GEPMass_YLess2[2][20][13];



   TH1F*  MassPtY_all[ptbins][ybins] ;
 

    for (Int_t ptindex = 0; ptindex < ptbins; ptindex++)
    {
	for (Int_t yindex = 0; yindex < ybins; yindex++)
	{
	    HistName = "Mass_pt_";
	    HistName += ptindex;
	    HistName += "_y";
	    HistName += yindex;
	    HistName += "_all";
	    MassPtY_all[ptindex][yindex] = new TH1F(HistName.Data(),HistName.Data(),massbins,massmin,massmax);

	}
    }

    for (Int_t cindex = 0; cindex < 2; cindex++)
    {
	for (Int_t ptindex = 0; ptindex < ptbins; ptindex++)
	{
	    for (Int_t yindex = 0; yindex < ybins; yindex++)
	    {
		HistName = "Mass_pt";
		HistName += ptindex;
                HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		MassPtY[cindex][ptindex][yindex] = new TH1F(HistName.Data(),HistName.Data(),massbins,massmin,massmax);

		HistName = "Mass_More_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		MassPtY_More2[cindex][ptindex][yindex] = new TH1F(HistName.Data(),HistName.Data(),massbins,massmin,massmax);

		HistName = "Mass_Less_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		MassPtY_Less2[cindex][ptindex][yindex] = new TH1F(HistName.Data(),HistName.Data(),massbins,massmin,massmax);

                HistName = "MassGEP_Pt_pt";
		HistName += ptindex;
                HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_Pt[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ptbins,ptmax,ptmin);

		HistName = "MoreMassGEP_Pt_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_PtMore2[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ptbins,ptmax,ptmin);

		HistName = "LessMassGEP_Pt_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_PtLess2[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ptbins,ptmax,ptmin);

                HistName = "MassGEP_Y_pt";
		HistName += ptindex;
                HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_Y[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ybins,ymax,ymin);

		HistName = "MoreMassGEP_Y_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_YMore2[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ybins,ymax,ymin);

		HistName = "LessMassGEP_Y_pt";
		HistName += ptindex;
		HistName += "_y";
		HistName += yindex;
		HistName += "_c";
		HistName += cindex;
		GEPMass_YLess2[cindex][ptindex][yindex] = new TH2F(HistName.Data(),HistName.Data(),massbins,massmin,massmax,ybins,ymax,ymin);

	    }
	}
    }


    Int_t centrality2[2][2]=  {{26,56},{56,102}};

    for (Long64_t ii = 1; ii <= nentries ; ii++)
    {
	if (ii != 0  &&  ii%5000000 == 0)    cout << "." << flush;
	if (ii != 0  &&  ii%50000000 == 0) {
	    Double_t event_percent = 100.0*((Double_t)ii)/((Double_t)nentries);
	    cout << " reco loop " << ii << " (" << event_percent << "%) " << "\n" << "==> Processing data" << flush;
	}

	nt.GetEntry(ii);

	if(pid.energyG < 200)                    continue;  // in the tree already 100MeV is applied
	if(pid.betaG < 0.95 )                    continue;  // to have good selection on gamma
	if(pid.betaG > 1.2 )                     continue;  // to have good selection on gamma
	if(pid.sizeG !=1)                        continue;  // still only cluster size 1 are calibrated, higer size are not calibrated

	// in ECAL we have only 1,2,4,5;  2 & 4 were well performing
	// if(!(pid.sectorG ==2)) continue;  // only those 2 sectors are good in ECAL
	if(!(pid.sectorG ==2 || pid.sectorG ==4)) continue;  // only those 2 sectors are good in ECAL

	hAll_GEPmass_en_mult_noOp->Fill(pid.massGEP,pid.mult_tofrpc);
	hAll_ptGEP_GEPmass_noOp->Fill(pid.massGEP,pid.ptGEP);
	hAll_yGEP_GEPmass_noOp->Fill(pid.massGEP,pid.yGEP);


	if(pid.angleGEP< 5  ) continue;
//	if(pid.angleGEP< 5 || pid.angleGEP> 60 ) continue;
	// applying higher energy condition

	hAll_GEPmass_en_mult[0]->Fill(pid.massGEP,pid.mult_tofrpc);
	hAll_ptGEP_GEPmass[0]->Fill(pid.massGEP,pid.ptGEP);
	hAll_yGEP_GEPmass[0]->Fill(pid.massGEP,pid.yGEP);

	if(pid.angleLL<2)
	{
	    hAll_GEPmass_en_mult[1]->Fill(pid.massGEP,pid.mult_tofrpc);
	    hAll_ptGEP_GEPmass[1]->Fill(pid.massGEP,pid.ptGEP);
	    hAll_yGEP_GEPmass[1]->Fill(pid.massGEP,pid.yGEP);
	    if(pid.mult_tofrpc > 26 && pid.mult_tofrpc <56)
	    {

		hAll_ptGEP_GEPmass_Less2Deg[0]->Fill(pid.massGEP,pid.ptGEP);
		hAll_yGEP_GEPmass_Less2Deg[0]->Fill(pid.massGEP,pid.yGEP);
		hAll_ptGEP_yGEP_Less2Deg[0]->Fill(pid.ptGEP,pid.yGEP);


	    }
	    else if(pid.mult_tofrpc > 26 && pid.mult_tofrpc <102)
	    {

		hAll_ptGEP_GEPmass_Less2Deg[1]->Fill(pid.massGEP,pid.ptGEP);
		hAll_yGEP_GEPmass_Less2Deg[1]->Fill(pid.massGEP,pid.yGEP);
		hAll_ptGEP_yGEP_Less2Deg[1]->Fill(pid.ptGEP,pid.yGEP);


	    }


	}
	if(pid.angleLL>2)
	{
	    hAll_GEPmass_en_mult[2]->Fill(pid.massGEP,pid.mult_tofrpc);
	    hAll_ptGEP_GEPmass[2]->Fill(pid.massGEP,pid.ptGEP);
	    hAll_yGEP_GEPmass[2]->Fill(pid.massGEP,pid.yGEP);

	    if(pid.mult_tofrpc > 26 && pid.mult_tofrpc <56)
	    {

		hAll_ptGEP_GEPmass_More2Deg[0]->Fill(pid.massGEP,pid.ptGEP);
		hAll_yGEP_GEPmass_More2Deg[0]->Fill(pid.massGEP,pid.yGEP);
		hAll_ptGEP_yGEP_More2Deg[0]->Fill(pid.ptGEP,pid.yGEP);


	    }
	    else if(pid.mult_tofrpc > 56 && pid.mult_tofrpc <102)
	    {

		hAll_ptGEP_GEPmass_More2Deg[1]->Fill(pid.massGEP,pid.ptGEP);
		hAll_yGEP_GEPmass_More2Deg[1]->Fill(pid.massGEP,pid.yGEP);
		hAll_ptGEP_yGEP_More2Deg[1]->Fill(pid.ptGEP,pid.yGEP);


	    }

	}
	if (pid.angleLL>4 )
	{
	    hAll_GEPmass_en_mult[3]->Fill(pid.massGEP,pid.mult_tofrpc);
	    hAll_ptGEP_GEPmass[3]->Fill(pid.massGEP,pid.ptGEP);
	    hAll_yGEP_GEPmass[3]->Fill(pid.massGEP,pid.yGEP);

	}
	if (pid.angleLL>6 )
	{
	    hAll_GEPmass_en_mult[4]->Fill(pid.massGEP,pid.mult_tofrpc);
	    hAll_ptGEP_GEPmass[4]->Fill(pid.massGEP,pid.ptGEP);
	    hAll_yGEP_GEPmass[4]->Fill(pid.massGEP,pid.yGEP);

	}

        if(pid.mult_tofrpc < 26 || pid.mult_tofrpc > 102) continue;
	//        if(pid.mult_tofrpc < 26 || pid.mult_tofrpc > 102) continue;

	// here you fill histograms

	// hAll_GEPmassAng_en150_mult->Fill(pid.massGEP,pid.angleGEP);
	beta->Fill(pid.momP,pid.betaP);
	beta->Fill(pid.momE*pid.chargeE,pid.betaE);
	hAll_BetaGamma->Fill(pid.betaG);
	hAll_massGEP->Fill(pid.massGEP);
	hAll_EEAng->Fill(pid.angleLL);
	hAll_GEPmassAngll->Fill(pid.massGEP,pid.angleLL);
	hAll_GEPmassMultmulttofrpc->Fill(pid.massGEP,pid.mult_tofrpc);
	hAll_ptGEPyGEP->Fill(pid.ptGEP,pid.yGEP);
	hAll_AngGEP_AngLL->Fill(pid.angleGEP,pid.angleLL);


	//           if ( (pid.mult_tofrpc > 20 && pid.mult_tofrpc <120) && (pid.angleGEP>5 && pid.angleGEP<60) && (pid.massLL<140) )
	//           {
	hAll_GEPmassyGEP->Fill(pid.massGEP,pid.yGEP);
	hAll_GEPmassptGEP->Fill(pid.massGEP,pid.ptGEP);
	hAll_GEPmassAng_en120_mult->Fill(pid.massGEP,pid.angleGEP);

	//	    }

	YBin  = (Int_t)((pid.yGEP-ymin)/YBinSize);
	PtBin = (Int_t)((pid.ptGEP-ptmin)/PtBinSize);

	if(PtBin >= 0 && YBin>=0 && PtBin < ptbins && YBin < ybins) {

	    MassPtY_all[PtBin][YBin]->Fill(pid.massGEP);
	    for(Int_t c=0; c<2; c++)
	    {

		MassPtY[c][PtBin][YBin]->Fill(pid.massGEP);
		GEPMass_Pt[c][PtBin][YBin]->Fill(pid.massGEP,pid.ptGEP);
		GEPMass_Y[c][PtBin][YBin]->Fill(pid.massGEP,pid.yGEP);


	     //   if(pid.mult_tofrpc > centrality2[c][0] && pid.mult_tofrpc < centrality2[c][1])
	     //   {
		    if(pid.angleLL < 2 )
		    {
			MassPtY_Less2[c][PtBin][YBin]->Fill(pid.massGEP);
			GEPMass_PtLess2[c][PtBin][YBin]->Fill(pid.massGEP,pid.ptGEP);
                        GEPMass_YLess2[c][PtBin][YBin]->Fill(pid.massGEP,pid.ptGEP);
		    }
		    else
		    {
			MassPtY_More2[c][PtBin][YBin]->Fill(pid.massGEP);
			GEPMass_PtMore2[c][PtBin][YBin]->Fill(pid.massGEP,pid.ptGEP);
                        GEPMass_YMore2[c][PtBin][YBin]->Fill(pid.massGEP,pid.ptGEP);
		    }

	      //  }
	    }
	}
    }

    TCanvas * Ct[6];

    for (int ct =0 ; ct< 6 ; ct++)
    {
	Ct[ct] = new TCanvas(Form("Ct%d",ct),Form("Ct%d",ct));

    }
    TH1F * projh2X;
    TH1F * proj_noOp;
    TH2F * dum;

    Int_t  centrality_arr[5][2]={{26,39},{ 39,56},{56,77},{77,102},{102,169}};

    //TString color[5]={"kBlue","kRed","kGrey","kBlack","kGreen"} ;
    const EColor colours[] = {kBlue,kRed,kGreen,kBlack,kTeal};
    for (int b=0; b<5; b++)
    { // loop over centrality
	Ct[b]->cd(); // canvas


	TLegend* leg = new TLegend ( 0.5,0.7,0.9,0.9 );
	leg->SetHeader ( "Opening angle condition" );
	leg->SetFillColor ( 1 );

	proj_noOp = (TH1F*) hAll_GEPmass_en_mult_noOp->ProjectionX(Form("cent%i_Nocut",b),centrality_arr[b][0],centrality_arr[b][1]);
       // proj_noOp->Draw();
	leg->AddEntry (proj_noOp,"No opening angl cndt","lep" );


	for (int a =0; a<5; a++)
	{ // loop over cuts

	    dum =    hAll_GEPmass_en_mult[a];
	    projh2X = (TH1F*) dum->ProjectionX(Form("cent%i_cut%i",b,a),centrality_arr[b][0],centrality_arr[b][1]);
	    projh2X->Draw("same");
	    projh2X->SetLineColor(colours[a]);
	    leg->AddEntry (projh2X,Form("Condition%d",a),"lep" );


	}



	leg->SetFillStyle ( 0 );
	leg->Draw();

	//        Ct[b]->GetYaxis()->SetRangeUser(0,proj_noOp);
	Ct[b]->Write();
    }


    out->cd();
    out->GetList()->Write();
    out->Save();

    // 1. beta of gamma                       1D
    // 2. beta vs. momentum for e+ and e- for 2D
    // 3. invaraint mass of GEP for           1D
    // 4. opening angle of e+e-               1D
    // 5. opening angle of e+e- vs. mass GEP  2D

    // 5a. opening angle of e+e- vs. op. ang. GEP

    // 6.  apply opening angle between e+e- > 2 and plot
    // 6a. GEPmass vs. mult
    // 6b. ptGEP vs. GEPmass
    // 6c. yGEP vs. GEPmass

    // 6.  apply opening angle between e+e- > 4 and plot
    // 6a. GEPmass vs. mult
    // 6b. ptGEP vs. GEPmass
    // 6c. yGEP vs. GEPmass

    // 6.  apply opening angle between e+e- > 6 and plot
    // 6a. GEPmass vs. mult
    // 6b. ptGEP vs. GEPmass
    // 6c. yGEP vs. GEPmass




}