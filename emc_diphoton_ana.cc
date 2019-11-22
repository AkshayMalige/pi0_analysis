
#include "fill_ttree_photonEpEm.C"
//#include "emc_recalibration.C"
#include "hades.h"
#include "hloop.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include <iostream>
using namespace std;

#ifndef __CINT__


int main(int argc, char **argv)
{
    TROOT Analysis_Tree("Analysis","compiled analysis macro");

    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1
    cout<<argc<<" arguments "<<endl;
    if(argc>1) cout<<"arg1 ="<<argv[1]<<endl;
    if(argc>2) cout<<"arg2 ="<<argv[2]<<endl;
    if(argc>3) cout<<"arg3 ="<<argv[3]<<endl;
    if(argc>4) cout<<"arg4 ="<<argv[4]<<endl;

    TString number ;
    TString nevts ;
    bool useMix;
    if(TString(argv[4]).Contains("yes")) useMix = true;
    else useMix=false;

    switch (argc)
    {
    case 5:       // just inputfile name + nEvents
	nevts  = argv[3];
	fill_ttree_diphoton(TString(argv[1]),TString(argv[2]),nevts.Atoi(),useMix);
	break;
    default: cerr<<"ERROR: analysis() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="",  nevents=1000, bool"<<endl;

	return 1; // fail
    }
}
#endif


