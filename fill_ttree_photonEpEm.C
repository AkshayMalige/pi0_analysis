
//#ifndef __CINT__

#include "hades.h"
#include "htool.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hgeantkine.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "hemcclustersim.h"
#include "TLorentzVector.h"
#include "hcategorymanager.h"
#include "TMacro.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticleevtinfo.h"
//#include "hgeantkine.h"

#include "hstart2cal.h"
#include "hstart2hit.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "simulationdef.h"
#include "tofdef.h"
#include "walldef.h"
#include "hemccluster.h"
#include "emcdef.h"
#include "hrpchit.h"
#include "htofhit.h"

#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLegendEntry.h"
//#include "PID_Cuts_gosia.h"


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
#include "hparticlebooker.h"
#include "htime.h"
#include "hparticleanglecor.h"


void doTrackCorr(HParticleAngleCor &trackcor, HCategory *candCat) {
    Int_t size = candCat->getEntries();
    for(Int_t j = 0; j < size; j ++){
	HParticleCand *cand;
	cand = HCategoryManager::getObject(cand,candCat,j);
	trackcor.recalcSetEmission(cand);  // changes values inside candidate
	trackcor.realignRichRing(cand);    //
    }
}


//#include "/u/harabasz/mar19/loopDST_workingcopy/selectFunctions.h"
static TH1F *hcutMassSys0 = NULL;
static TH1F *hcutMassSys1 = NULL;
static TH1F *hcutRichQaSys0 = NULL;
static TH1F *hcutRichQaSys1 = NULL;

static Bool_t selectLeptonsBetaLoose(HParticleCand* pcand){

    //  selection function for lepton candidates.

    Bool_t selectEpEm = kFALSE;

    // if(pcand->isFakeRejected()) return kFALSE;
    if(pcand->isFlagAND(5,
			Particle::kIsAcceptedHitRICH,
			Particle::kIsAcceptedHitInnerMDC,
			Particle::kIsAcceptedHitOuterMDC,
			Particle::kIsAcceptedHitMETA,
			Particle::kIsAcceptedRK)
       &&
       pcand->getInnerSegmentChi2() > 0
       &&
       pcand->getOuterSegmentChi2() > 0
       &&
       pcand->getMetaMatchQuality() < 4
       &&
       pcand->getMetaMatchQualityEmc() < 4
       &&
       HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
       &&
       pcand->getChi2() < 1000
       &&
       pcand->getBeta() > 0.9
       &&
       pcand->getMomentum() > 0.
       &&
       pcand->getMomentum() < 1000.
      ) selectEpEm = kTRUE;

    return selectEpEm;
}
void setupMassCuts(TString filename) {
    TFile *fileCutMass = new TFile(filename);
    if (fileCutMass) {
	hcutMassSys0 = (TH1F*)fileCutMass->Get("m2_cut_sys0");
	hcutMassSys1 = (TH1F*)fileCutMass->Get("m2_cut_sys1");
	if (!hcutMassSys0 || !hcutMassSys1) {
	    cout << "No histogram with mass cuts found in file. Exitting..." << endl;
	    exit(-1);
	}
    }
    else {
	cout << "No file with mass cuts found. Exiting..." << endl;
	exit(-1);
    }
    //   hcutMassSys1->Scale(0.25);
}


typedef struct { Float_t sector,beta,cell,energy,time,size,theta,phi,matchCell,matchPart;
    TLorentzVector v;} photon;

typedef struct { Float_t sector, beta, theta, phi, charge, richQa, mom;
    TLorentzVector e;} lep;


void fill_ttree_diphoton(TString inFile, TString outputFile, Int_t nev=1000000., Bool_t isMix=kFALSE)
{

    TROOT dst_production( "DstProductionApp", "Produce a DST File" );

    Bool_t isSimulation          = kFALSE;
    Bool_t createHades = kTRUE; //kFALSE;  // kTRUE = create HADES new
    HLoop      *loop = new HLoop(createHades); // create HADES (needed if one wants to use gHades)

    cout<<"plik inputowy : "<<gSystem->BaseName(inFile.Data())<<endl;

    loop->addMultFiles(inFile);
    //    loop->addFiles(files);
    if(!loop->setInput(""))  exit(1);
    HEventHeader *evHeader = loop->getEventHeader();

    HCategory *catEmcCal  = loop->getCategory("HEmcCal");
    HCategory *catEmcClus = loop->getCategory("HEmcCluster");
    HCategory *catRpcClus = loop->getCategory("HRpcCluster");
    HCategory *catRpcHit  = loop->getCategory("HRpcHit");
    HCategory *catTofHit  = loop->getCategory("HTofHit");
    HCategory *catStartHit= loop->getCategory("HStart2Hit");
    HCategory *candCat    = loop->getCategory("HParticleCand");
    HCategory *catGeantKine = loop->getCategory("HGeantKine");
    HParticleCand* cand1;

    HGeantKine*       kine;
    //HLinearCategory* kineCat = (HLinearCategory*)HCategoryManager::getCategory(catGeantKine);

    if(catGeantKine) {isSimulation=kTRUE; Info("init()","GeantKine found - is Simulation");}
    else        {isSimulation=kFALSE;Info("init()","GeantKine not found - is not Simulation");}

    HGeantHeader* fSubHeader = NULL;
    HPartialEvent* fSimul = NULL;

    //old calibration

//    TString outputFileName = fileWoPath+"_"+outputFile;

    TFile* out = new TFile(outputFile,"RECREATE");
    out->cd();

    typedef struct {Float_t mult_tofrpc,mult_emc,sectorP,betaP,thetaP,phiP,chargeP,richQaP,momP,
    sectorE,betaE,thetaE,phiE,chargeE,richQaE,momE,
    sectorG,betaG,cellG,energyG,timeG,sizeG,thetaG,phiG,matchCellG,matchPartG,
    massLL,angleLL,ptLL,yLL,thetaLL,phiLL,momLL,
    massGEP,angleGEP,ptGEP,yGEP,thetaGEP,phiGEP,momGEP;} PID;

    static PID pid;

    TTree *nt = new TTree("nt","pid");
    nt->Branch("pid",&pid,"mult_tofrpc:mult_emc:sectorP:betaP:thetaP:phiP:chargeP:richQaP:momP:sectorE:betaE:thetaE:phiE:chargeE:richQaE:momE:sectorG:betaG:cellG:energyG:timeG:sizeG:thetaG:phiG:matchCellG:matchPartS:massLL:opangleLL:ptLL:yLL:thetaLL:phiLL:momLL:massGEP:opangleGEP:ptGEP:yGEP:thetaGEP:phiGEP:momGEP");

    Int_t cellMap[163] = {
	6,   7,   8,   9,  10,
	23,  24,  25,  26,  27,
	39,  40,  41,  42,  43,  44,  45,
	56,  57,  58,  59,  60,  61,  62,
	72,  73,  74,  75,  76,  77,  78,  79,  80,
	89,  90,  91,  92,  93,  94,  95,  96,  97,
	106, 107, 108, 109, 110, 111, 112, 113, 114,
	122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
	139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
	155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
	172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184,
	188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202,
	205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
	221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
	238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254
    };

    Int_t systemMap[255];
    for (int i = 0; i < 255; ++i) systemMap[i] = -1;
    for (int i = 0; i < 163; ++i) systemMap[cellMap[i]] = i;

    Int_t sectorMap[6];
    Int_t secMap[4]={1,2,4,5};

    for (int i = 0; i < 6; ++i) sectorMap[i] = -1;
    for (int i = 0; i < 4; ++i) sectorMap[secMap[i]] = i;

    Float_t emc_ModMult[6][163];

    setupMassCuts("m2_cut.root");

    //--------------------------CONFIGURATION---------------------------------------------------
    //At begin of the program (outside the event loop)
    HParticleTrackSorter sorter;
    sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
    sorter.init();                                                  // get catgegory pointers etc...

    HParticleAngleCor trackcor;
    trackcor.setDefaults("apr12");


    //=======================================================================
    Int_t nEntries = loop->getEntries();
    printf("nEntries = %i\n",nEntries);
    if(nEntries>nev && nev>0) nEntries = nev;
    for (Int_t i=0; i<nEntries ; i++) if(loop->nextEvent(i) > 0)  {

	//            cout<<"kuku0"<<endl;

	if(loop->nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
	HTool::printProgress(i,nEntries,1,"Analyze pairs :");

        if (!catStartHit || catStartHit->getEntries()<1)  continue;

	if(gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ() < -70) continue;

	//	cout<<"i : "<<i<<endl;

	HParticleEvtInfo* evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);

	Int_t multMeta =  evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
	pid.mult_tofrpc=multMeta;

	pid.mult_emc=catEmcClus->getEntries();


	// different flags for leptons
	Int_t sizeCand = candCat->getEntries();
	Bool_t isUsed[sizeCand];
	Bool_t isUsedPion[sizeCand];
	Bool_t richQaFlag[sizeCand];
	Bool_t betaFlag[sizeCand];
	Bool_t betaFlagTof[sizeCand];
	Bool_t betaFlagRpc[sizeCand];
	Bool_t showerFlag[sizeCand];
	Bool_t leptonFlag[sizeCand];
	Bool_t isUsed_noRich[sizeCand];

	Int_t mult_photon =0;

	for(Int_t i=0; i<6; i++){
	    for(Int_t j=0; j<163; j++){
		emc_ModMult[i][j] = 0.;
	    }
	}

	//////////////////////////////////////////////////////////////////////////////////
	//      this part is for leptons
	//
	//////////////////////////////////////////////////////////////////////////////////

      //  doTrackCorr(trackcor, candCat); // multiple scatering correction for leptons

	sorter.cleanUp();
	sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
        sorter.fill(selectLeptonsBetaLoose);   // fill only good leptons
	sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);
//	sorter.setIgnoreInnerMDC();                            // do not reject Double_t inner MDC hits
//	sorter.init();                                         // get catgegory pointers etc...

	vector<lep> electrons;
	vector<lep> positrons;
	vector<photon> gammas;
	vector<photon> dileptons;

	Int_t sizeCatCand = candCat->getEntries();

	for(Int_t j = 0; j < sizeCand; j ++){ // loop over charge particle candidate
	    richQaFlag[j]                     = kFALSE;
	    betaFlag[j]                       = kFALSE;
	    betaFlagTof[j]                    = kFALSE;
	    betaFlagRpc[j]                    = kFALSE;
	    leptonFlag[j]                     = kFALSE;
            isUsed[j]                         = kFALSE;

	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

	    if(cand1->getSystemUsed()==-1)  continue;

	    lep ll;

	    Float_t richQa  = cand1->getRichMatchingQuality();
	    Float_t _mom    = cand1->getMomentum();
	    Float_t _beta   = cand1->getBeta();
	    Float_t m2      = _mom*_mom*(1-_beta*_beta)/(_beta*_beta);
	    Int_t charge    = cand1->getCharge();
	    Int_t sys       = cand1->getSystemUsed();
	    Int_t sec       = cand1->getSector();
	    Int_t th        = cand1->getTheta();
	    Int_t ph        = cand1->getPhi();

	    ll.sector = sec;
	    ll.charge = charge;
	    ll.beta   = _beta;
	    ll.richQa = richQa;
	    ll.theta  = th;
	    ll.phi    = ph;
	    ll.mom    = _mom;

	    cand1->calc4vectorProperties(HPhysicsConstants::mass(3));
            ll.e      = *cand1;
	    //            if (!cand1->isFlagBit(Particle::kIsLepton)) continue;

	    ///////////////////////////////////////////////////////////////////////////////////
	    //
	    // Hard cut flags
	    //
	    ///////////////////////////////////////////////////////////////////////////////////
	    // Ring matching

            if(!selectLeptonsBetaLoose(cand1)) continue;
	    isUsed[j] = cand1->isFlagBit(Particle::kIsLepton);

	    if (richQa < 1. && richQa > -1) richQaFlag[j]=kTRUE;

	    // beta
	    Int_t charge_index;
	    if (charge < 0) {
		charge_index = 0;
	    }
	    else {
		charge_index = 1;
	    }
  //          cout<<"sys1 : "<<hcutMassSys1<<"  sys0  "<<hcutMassSys0<<endl;

	    if (sys==1) {
		if (m2 < hcutMassSys1->Interpolate(_mom*charge)) betaFlagTof[j] = kTRUE;
	    }
	    if (sys==0) {
		if (m2 < hcutMassSys0->Interpolate(_mom*charge)) betaFlagRpc[j] = kTRUE;
	    }
	    if (betaFlagTof[j] || betaFlagRpc[j]) betaFlag[j] = kTRUE;

//	    cout<<"flagi  "<<isUsed[j] <<"  "<< richQaFlag[j]<<" richQa " <<richQa <<"  "<< betaFlag[j]<<endl;


	    if (isUsed[j] && richQaFlag[j] && betaFlag[j]) {
		if (cand1->getCharge() == -1) {
		    electrons.push_back(ll);

		}
		if (cand1->getCharge() == 1) {
		    positrons.push_back(ll);
		}
	    }
	} // end of loop over ParticleCand


	//////////////////////////////////////////////////////////////////////////////////
	//
	//   this part is for gamma
	//
	//////////////////////////////////////////////////////////////////////////////////


	sorter.cleanUp();
	sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
	//        sorter.fill(selectLeptonsBetaLoose);   // fill only good leptons
	sorter.fill(HParticleTrackSorter::selectHadrons);
	sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

	Int_t size = candCat->getEntries();

	for(Int_t j = 0; j < size; j ++){ // loop over Particlecand
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

	    if(!cand1->isFlagBit(Particle::kIsUsed)) continue;

	    Int_t emc_index = cand1->getEmcInd();
	    Double_t emc_mq = cand1->getMetaMatchQualityEmc();
	    if(emc_index >=0 && emc_mq < 2.5) {
		HEmcCluster* emcClst = (HEmcCluster*)catEmcClus->getObject(emc_index);
		Int_t sec = cand1->getSector();
		Int_t cell= emcClst->getCell();

		emc_ModMult[sec][systemMap[cell]]++ ;
	    }
	}

	Int_t newEvent = 1;
	Int_t nEmcCluster = catEmcClus->getEntries();

	for(Int_t e=0; e<nEmcCluster-1; e++) { // loop over EMC cluster for gamma
	    HEmcCluster* cls = (HEmcCluster*)catEmcClus->getObject(e);
	    HEmcClusterSim* clsSim = dynamic_cast<HEmcClusterSim*>(cls);

	    photon g;

	    g.sector     = cls->getSector();
	    g.size       = cls->getNCells(); // = cluster size
	    g.phi        = cls->getPhi();
	    g.theta      = cls->getTheta();
	    g.cell       = systemMap[cls->getCell()];
	    g.matchCell  = cls->getNMatchedCells();
	    g.matchPart  = emc_ModMult[cls->getSector()][systemMap[cls->getCell()]];
	    g.time          = cls->getTime();
	    g.energy        = cls->getMaxEnergy();

	    // Calculate Beta:
	    HGeomVector trackVec(cls->getXLab(),cls->getYLab(),cls->getZLab()  +30.);
	    Double_t trackLength = trackVec.length();
	    trackVec  *= (g.energy/trackLength);

	    Double_t Beta        = (trackLength/1000.)/(g.time*1.e-9) / TMath::C();
	    g.beta=Beta;
            g.v.SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);

            if(g.beta > 0.9 && g.beta < 1.4 && g.energy > 100 && g.matchCell == 0 && g.matchPart==0) gammas.push_back(g);
	} // end of loop over EMC cluster


        Int_t sizeEle = electrons.size();
        Int_t sizePos = positrons.size();
	Int_t sizeGam = gammas.size();

        //cout<<sizeEle<<"  "<<sizePos<<"  "<<sizeGam<<endl;

	for(Int_t i=0; i<sizeEle; i++){
	    for(Int_t j=0; j<sizePos; j++){

		TLorentzVector pairLL = electrons[i].e + positrons[j].e;

		pid.sectorP = positrons[j].sector;
		pid.betaP   = positrons[j].beta;
		pid.thetaP  = positrons[j].theta;
		pid.phiP    = positrons[j].phi;
		pid.chargeP = positrons[j].charge;
		pid.richQaP = positrons[j].richQa;
		pid.momP    = positrons[j].mom;
		pid.sectorE = electrons[j].sector;
		pid.betaE   = electrons[j].beta;
		pid.thetaE  = electrons[j].theta;
		pid.phiE    = electrons[j].phi;
		pid.chargeE = electrons[j].charge;
		pid.richQaE = electrons[j].richQa;
		pid.momE    = electrons[j].mom;


		Double_t angle = electrons[i].e.Angle(positrons[j].e.Vect())*TMath::RadToDeg();
		Double_t mass  = pairLL.M();
		Float_t pt     = pairLL.Pt();
		Float_t y      = pairLL.Rapidity();

                pid.massLL  = mass;
                pid.angleLL = angle;
                pid.ptLL    = pt;
		pid.yLL     = y;
		pid.thetaLL = pairLL.Theta();
		pid.phiLL   = pairLL.Phi();
		pid.momLL   = pairLL.P();


		for(Int_t k=0; k<sizeGam; k++){

                    TLorentzVector pairGEP = pairLL + gammas[k].v;

                    pid.sectorG     = gammas[k].sector;
                    pid.betaG       = gammas[k].beta;
                    pid.cellG       = gammas[k].cell;
                    pid.sizeG       = gammas[k].size;
                    pid.energyG     = gammas[k].energy;
                    pid.timeG       = gammas[k].time;
                    pid.thetaG      = gammas[k].theta;
                    pid.phiG        = gammas[k].phi;
                    pid.matchCellG  = gammas[k].matchCell;
                    pid.matchPartG  = gammas[k].matchPart;

		    Double_t angle  = pairLL.Angle(gammas[k].v.Vect())*TMath::RadToDeg();
		    Double_t mass   = pairGEP.M();
		    Double_t pt     = pairGEP.Pt();
		    Double_t y      = pairGEP.Rapidity();

		    //cout<<"mass : "<<mass<<endl;

		    pid.massGEP  = mass;
		    pid.angleGEP = angle;
		    pid.ptGEP    = pt;
		    pid.yGEP     = y;
		    pid.thetaGEP = pairGEP.Theta();
		    pid.phiGEP   = pairGEP.Phi();
		    pid.momGEP   = pairGEP.P();
		    nt->Fill();

		} // end of loop over gamma
	    } // end of loop over positron
	} // end of look over electron
    }  // end of event loop

    out->cd();
    out->GetList()->Write();
    TMacro *m = new TMacro("fill_ttree_photonEpEm.C");
    TTree macros("macros","macros");
    macros.GetUserInfo()->Add(m);
    macros.Write();
    out->Save();
    out->Close();

    return ;

}