#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "constants.h"

using namespace std;
using namespace constants;

void calculate_targets_flux() {

    double density = 1.3836; //g/cm^3
    double volume = 5.82515e7; //cm^3
    double m_mol = 39.95; //g/mol
    double NA = 6.022e23; // N/mol

    double Nt = density * volume * NA / m_mol;
    cout << "Nt = " << Nt << endl; // 1.2149e+30 argon targets   
    
	TFile* FluxFile = TFile::Open("../super_unified/mySTVAnalysis/mcc9_10/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));    
	double DataPOT = 1.2e20;						
	double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);  
    cout << "IntegratedFlux = " << IntegratedFlux << endl;
    cout << "IntegratedFlux * Nt = " << IntegratedFlux * Nt << endl << endl;

    double FVx = 256., FVy = 232, FVz = 1037., borderx = 25., bordery = 25., borderz = 25.; // cm
    double mod_volume = (FVx - 2* borderx) * (FVy - 2*bordery) * (FVz - 2* borderz); // cm^3

    double mod_Nt = density * mod_volume * NA / m_mol;
    cout << "mod Nt = " << mod_Nt << endl; // 1.2149e+30 argon targets

}