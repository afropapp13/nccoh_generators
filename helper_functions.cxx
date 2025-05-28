#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <vector>

#include "constants.h"

using namespace std;
using namespace constants;

//----------------------------------------//

double find_bin_max_value(TH1D* h){

	int NBins = h->GetXaxis()->GetNbins();
	double HistoMax = -999.;	
	int bin_max = -1;

	for (int ibin = 1; ibin<= NBins; ibin++) {

		double LocalMax = h->GetBinContent(ibin);
		if (LocalMax > HistoMax) { HistoMax = LocalMax; bin_max = ibin; }

	}

	return HistoMax;

}
//----------------------------------------//

double PeLEE_ReturnBeamOnRunPOT(TString Run) {

	double DataPOT = -99.;

	if ( Run.Contains("Run1") ) { DataPOT = Fulltor860_wcut_Run1 ; }
	if ( Run.Contains("Run1A_open_trigger") ) { DataPOT = Fulltor860_wcut_Run1A_open_trigger ; }
	if ( Run.Contains("Run1B_open_trigger") ) { DataPOT = Fulltor860_wcut_Run1B_open_trigger ; }
	if ( Run.Contains("Run2") ) { DataPOT = Fulltor860_wcut_Run2 ; }
	if ( Run.Contains("Run3") ) { DataPOT = Fulltor860_wcut_Run3 ; }
	if ( Run.Contains("Run4a") ) { DataPOT = Fulltor860_wcut_Run4a ; }	
	if ( Run.Contains("Run4b") ) { DataPOT = Fulltor860_wcut_Run4b ; }
	if ( Run.Contains("Run4b_standalone") ) { DataPOT = Fulltor860_wcut_mcc9_10_Run4b_standalone; }	
	if ( Run.Contains("Run4b_unified") ) { DataPOT = Fulltor860_wcut_mcc9_10_Run4b_unified; }			
	if ( Run.Contains("Run4c") ) { DataPOT = Fulltor860_wcut_Run4c ; }	
	if ( Run.Contains("Run4d") ) { DataPOT = Fulltor860_wcut_Run4d ; }		
	if ( Run.Contains("Run5") ) { DataPOT = Fulltor860_wcut_Run5 ; }
	if ( Run.Contains("Combined") ) { DataPOT = Fulltor860_wcut_Combined ; }

	return DataPOT;

}
//----------------------------------------//

double get_median(const TH1D * h1) { 

	int n = h1->GetXaxis()->GetNbins();  
	std::vector<double>  x(n);
 	h1->GetXaxis()->GetCenter( &x[0] );
	const double * y = h1->GetArray(); 
	// exclude underflow/overflows from bin content array y
	return TMath::Median(n, &x[0], &y[1]); 
 
}

//----------------------------------------//

//Function to multiply by the bin width
void rm_bin_width(TH1D* h, double SF = 1.) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = SF * CurrentEntry * h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = SF * CurrentError * h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}

//----------------------------------------//

//Function to divide by the bin width and to get xsecs
void divide_bin_width(TH1D* h, double SF = 1.) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = SF * CurrentEntry / h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = SF * CurrentError / h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}

//----------------------------------------//

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
}

//----------------------------------------//

double round(double var,int acc = 0) 
{ 
    double value = (int)(var * TMath::Power(10.,acc) + .5); 
    return (double)value / TMath::Power(10.,acc); 
} 

//----------------------------------------//


