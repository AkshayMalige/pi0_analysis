

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


    typedef struct {Float_t mult_tofrpc,mult_emc,sectorP,betaP,thetaP,phiP,chargeP,richQaP,momP,
    sectorE,betaE,thetaE,phiE,chargeE,richQaE,momE,
    sectorG,betaG,cellG,energyG,timeG,sizeG,thetaG,phiG,matchCellG,matchPartG,
    massLL,angleLL,ptLL,yLL,thetaLL,phiLL,momLL,
    massGEP,angleGEP,ptGEP,yGEP,thetaGEP,phiGEP,momGEP;} PID;


const Double_t d2r = 0.01745329;

void plot_mass_gee(TString inFile, TString outputFile, Int_t nev=100000000., Bool_t useMix=kFALSE)
{

    //data

    TChain nt = TChain("nt");
    nt.AddFile(inFile);

    PID pid;
    nt.SetBranchAddress("pid" , &pid.mult_tofrpc);
    Long64_t nentries        = nt.GetEntries();

    cout<<"nentries:  "<<nentries<<endl;

    TFile *out    = new TFile("diphoton_massEXP_gen2.root","RECREATE");
    out->cd();

    // histogram definition

    TH2F* hAll_GEPmassAng_en150_mult = new TH2F("hAll_GEPmassAng_en150_mult",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];#alpha _{e^{+} e^{-}}",300,0,1000.0,150,0,150);
    TH2F* hAll_GEPmass_en150_mult    = new TH2F("hAll_GEPmass_en150_mult",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}]",300,0,1000.0,150,0,150);
    TH1F* hAll_BetaGamma    = new TH1F("hAll_BetaGamma",";#beta #gamma",500,0.8,1.3);
    TH1F* hAll_massGEP    = new TH1F("hAll_massGEP",";M_{#gamma e^{+}e^{-} [MeV/c]}",1000,0,1000);
    TH1F* hAll_EEAng    = new TH1F("hAll_EEAng",";#alpha _{e^{+}e^{-}}",150,0,150);
    TH2F* hAll_GEPmassAngll = new TH2F("hAll_GEPmassAngll",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];#alpha _{e^{+}e^{-}}",1000,0,1000.0,150,0,150);
    TH2F* hAll_GEPmassMultmulttofrpc = new TH2F("hAll_GEPmassMultmulttofrpc",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];Mult^{TOF+RPC}",1000,0,1000.0,160,0,160);
    TH2F* hAll_GEPmassMultmulttofrpc0 = new TH2F("hAll_GEPmassMultmulttofrpc0",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];Mult^{TOF+RPC}",1000,0,1000.0,160,0,160);

    TH2F* hAll_GEPmassyGEP = new TH2F("hAll_GEPmassyGEP",";M_{#gamma e^{+}e^{-}} [MeV/c^];y",1000,0,1000.0,300,0,3);
    TH2F* hAll_GEPmassptGEP = new TH2F("hAll_GEPmassptGEP",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];p_{t}",1000,0,1000.0,1000,0,1000);
    TH2F* hAll_GEPmassAng_en120_mult = new TH2F("hAll_GEPmassAng_en120_mult",";M_{#gamma e^{+}e^{-}} [MeV/c^{2}];#alpha _{#gamma e^{+} e^{-}}",300,0,600.0,150,0,150);
    TH2F* hAll_ptGEPyGEP = new TH2F("hAll_ptGEPyGEP",";p_{t} [MeV/c];y",600,200,800.0,8,0.4,2);
    TH2F *beta = new TH2F("beta","beta;pxq/|q| (MeV/c);#beta",400,-1000,1000,100,0.8,1.4);

    TH2F* hAll_AngGEP_AngLL = new TH2F("hAll_AngGEP_AngLL",";#alpha _{#gamma e^{+} e^{-}};#alpha _{e^{+} e^{-}}",140,0,140,140,0,140);

    TH2F*  hAll_GEPmass_en_mult[3];
    TH2F*  hAll_ptGEP_GEPmass[3];
    TH2F*  hAll_yGEP_GEPmass[3];
    Int_t  angl_limit = 2;
    for (int a =0; a<3; a++)
    {
	hAll_GEPmass_en_mult[a] = new TH2F ( Form ( "hAll_GEPmass_en_mult_%d", angl_limit ) , Form ( "hAll_GEPmass_en_mult_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];Mult^{TOF+RPC}", angl_limit ), 1000, 0,1000,160,0,160 );
	hAll_ptGEP_GEPmass[a] = new TH2F ( Form ( "hAll_ptGEP_GEPmass_%d", angl_limit ) , Form ( "hAll_ptGEP_GEPmass_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];p_{t}", angl_limit ), 1000, 0,1000,16,0,800 );
	hAll_yGEP_GEPmass[a] = new TH2F ( Form ( "hAll_yGEP_GEPmass_%d", angl_limit ) , Form ( "hAll_yGEP_GEPmass_%d;M_{#gamma e^{+}e^{-}} [MeV/c^];y", angl_limit ), 1000, 0,1000,10,0,2 );
        angl_limit = angl_limit+2;


    }


    Int_t ptbins = 20.; // for pi0 reduced
    Int_t ptmin  = 0;
    Int_t ptmax  = 800;

    //Int_t ybins    = 15.; // for eta
    Int_t ybins    = 20.; // for pi0
    Double_t ymin  = 0.;  //0.015
    Double_t ymax  = 3.;  //2.015
    Double_t deltay = 0.2;

    Double_t PtBinSize = (ptmax-ptmin)/ptbins; // MeV/bin
    Double_t YBinSize  = (ymax-ymin)/ybins;  // rapidity/bin

    Int_t massbins = 40;
    Int_t massmin  = 0;
    Int_t massmax  = 800;

    Int_t YBin=0;
    Int_t PtBin=0;
    TString HistName;
    TH1F *MassPtY[20][15];

    for (Int_t ptindex = 0; ptindex < ptbins; ptindex++)
    {
	for (Int_t yindex = 0; yindex < ybins; yindex++)
	{
	    HistName = "Mass_pt";
	    HistName += ptindex;
	    HistName += "_y";
	    HistName += yindex;
	    MassPtY[ptindex][yindex] = new TH1F(HistName.Data(),HistName.Data(),massbins,massmin,massmax);

	}
    }



    for (Long64_t ii = 1; ii <= nentries ; ii++)
    {
	if (ii != 0  &&  ii%5000000 == 0)    cout << "." << flush;
	if (ii != 0  &&  ii%50000000 == 0) {
	    Double_t event_percent = 100.0*((Double_t)ii)/((Double_t)nentries);
	    cout << " reco loop " << ii << " (" << event_percent << "%) " << "\n" << "==> Processing data" << flush;
	}

	nt.GetEntry(ii);

	if(pid.energyG < 100)                     continue;  // in the tree already 100MeV is applied
	if(pid.betaG < 0.95 )                     continue;  // to have good selection on gamma
	if(pid.sizeG !=1)                         continue;  // still only cluster size 1 are calibrated, higer size are not calibrated

	// in ECAL we have only 1,2,4,5;  2 & 4 were well performing
	if(!(pid.sectorG ==2 || pid.sectorG ==4)) continue;  // only those 2 sectors are good in ECAL

	// applying higher energy condition
	if(pid.energyG > 150) {

	    hAll_GEPmass_en150_mult->Fill(pid.massGEP,pid.mult_tofrpc);
	   // hAll_GEPmassMultmulttofrpc0->Fill(pid.massGEP,pid.mult_tofrpc);


	    if(pid.mult_tofrpc < 30 || pid.mult_tofrpc>120) continue;

            // here you fill histograms

	    hAll_GEPmassAng_en150_mult->Fill(pid.massGEP,pid.angleGEP);
            beta->Fill(pid.momP,pid.betaP);
	    beta->Fill(pid.momE*pid.chargeE,pid.betaE);
	    hAll_BetaGamma->Fill(pid.betaG);
	    hAll_massGEP->Fill(pid.massGEP);
	    hAll_EEAng->Fill(pid.angleLL);
	    hAll_GEPmassAngll->Fill(pid.massGEP,pid.angleLL);
	    hAll_GEPmassMultmulttofrpc->Fill(pid.massGEP,pid.mult_tofrpc);
	    hAll_ptGEPyGEP->Fill(pid.ptGEP,pid.yGEP);
	    hAll_AngGEP_AngLL->Fill(pid.angleGEP,pid.angleLL);


	    if(pid.angleLL<2)
	    {
		hAll_GEPmass_en_mult[0]->Fill(pid.massGEP,pid.mult_tofrpc);
		hAll_ptGEP_GEPmass[0]->Fill(pid.massGEP,pid.ptGEP);
                hAll_yGEP_GEPmass[0]->Fill(pid.massGEP,pid.yGEP);
	    }
	    else if (pid.angleLL>2 && pid.angleLL<4)
	    {
		hAll_GEPmass_en_mult[1]->Fill(pid.massGEP,pid.mult_tofrpc);
		hAll_ptGEP_GEPmass[1]->Fill(pid.massGEP,pid.ptGEP);
                hAll_yGEP_GEPmass[1]->Fill(pid.massGEP,pid.yGEP);

	    }
	    else if (pid.angleLL>4 && pid.angleLL<6)
	    {
		hAll_GEPmass_en_mult[2]->Fill(pid.massGEP,pid.mult_tofrpc);
		hAll_ptGEP_GEPmass[2]->Fill(pid.massGEP,pid.ptGEP);
                hAll_yGEP_GEPmass[2]->Fill(pid.massGEP,pid.yGEP);

	    }




 //           if ( (pid.mult_tofrpc > 20 && pid.mult_tofrpc <120) && (pid.angleGEP>5 && pid.angleGEP<60) && (pid.massLL<140) )
 //           {
               hAll_GEPmassyGEP->Fill(pid.massGEP,pid.yGEP);
	       hAll_GEPmassptGEP->Fill(pid.massGEP,pid.ptGEP);
               hAll_GEPmassAng_en120_mult->Fill(pid.massGEP,pid.angleGEP);

//	    }

            YBin  = (Int_t)((pid.yGEP-ymin)/YBinSize);
	    PtBin = (Int_t)((pid.ptGEP-ptmin)/PtBinSize);

	    if(PtBin >= 0 && YBin>=0 && PtBin < ptbins && YBin < ybins) MassPtY[PtBin][YBin]->Fill(pid.massGEP);

	}

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