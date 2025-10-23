#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

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

void normalize_rows(TH2D* hist) {

    if (!hist) {
        std::cerr << "Error: null histogram pointer.\n";
        return;
    }

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();

    // Loop over Y bins (rows)
    for (int iy = 1; iy <= nBinsY; ++iy) {
        double rowSum = 0.0;

        // Sum over X bins
        for (int ix = 1; ix <= nBinsX; ++ix) {
            rowSum += hist->GetBinContent(ix, iy);
        }

        // Avoid division by zero
        if (rowSum == 0) continue;

        // Normalize this row
        for (int ix = 1; ix <= nBinsX; ++ix) {
            double content = hist->GetBinContent(ix, iy);
            double error   = hist->GetBinError(ix, iy);

            hist->SetBinContent(ix, iy, content / rowSum);
            hist->SetBinError(ix, iy, error / rowSum);
        }
    }
}


//----------------------------------------//

void normalize_columns(TH2D* hist) {

    if (!hist) {
        std::cerr << "Error: null histogram pointer.\n";
        return;
    }

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();

    // Loop over X bins (columns)
    for (int ix = 1; ix <= nBinsX; ++ix) {
        double colSum = 0.0;

        // Sum over Y bins
        for (int iy = 1; iy <= nBinsY; ++iy) {
            colSum += hist->GetBinContent(ix, iy);
        }

        // Avoid division by zero
        if (colSum == 0) continue;

        // Normalize this column
        for (int iy = 1; iy <= nBinsY; ++iy) {
            double content = hist->GetBinContent(ix, iy);
            double error   = hist->GetBinError(ix, iy);

            hist->SetBinContent(ix, iy, content / colSum);
            hist->SetBinError(ix, iy, error / colSum);
        }
    }
}

//----------------------------------------//

void calc_chi2(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval, double &sigma) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {			

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	//TMatrixD inverse_cov_m = cov_m.Invert();

  	TMatrixD inverse_cov_m = cov_m;

	TDecompSVD svd(cov_m);
	inverse_cov_m = TMatrixDSym(cov_m.GetNrows(), svd.Invert().GetMatrixArray());

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;
		}

	}

	ndof = h_data_clone->GetNbinsX();
	pval = TMath::Prob(chi, ndof);
	sigma = TMath::Sqrt( TMath::ChisquareQuantile( 1-pval, 1 ) ); 

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

//----------------------------------------//

int locate_bin_with_value(TH1D* h, double Value) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 1; i <= NBins; i++) {

    double CurrentEntry = h->GetBinContent(i);
    if (CurrentEntry == Value) { return i; } 

  }

  return -99;

}

//----------------------------------------//

void bin_number_x_title(TH2D* h) {

	TString hname = h->GetName();
	if (string(hname).find("Serial") != std::string::npos) {	

		TString axis_title = h->GetXaxis()->GetTitle();
		axis_title.ReplaceAll("deg","bin #");
		axis_title.ReplaceAll("(GeV/c)","bin #");
		axis_title.ReplaceAll("GeV","bin #");				
		h->GetXaxis()->SetTitle(axis_title);

	}

}

//----------------------------------------//

void bin_number_y_title(TH2D* h) {

	TString hname = h->GetName();
	if (string(hname).find("Serial") != std::string::npos) {	

		TString axis_title = h->GetYaxis()->GetTitle();
		axis_title.ReplaceAll("deg","bin #");
		axis_title.ReplaceAll("(GeV/c)","bin #");
		axis_title.ReplaceAll("GeV","bin #");				
		h->GetYaxis()->SetTitle(axis_title);

	}

}
//----------------------------------------//

void bin_number_x_title(TH1D* h) {

	TString hname = h->GetName();
	if (string(hname).find("Serial") != std::string::npos) {	

		TString axis_title = h->GetXaxis()->GetTitle();
		axis_title.ReplaceAll("deg","bin #");
		axis_title.ReplaceAll("(GeV/c)","bin #");
		axis_title.ReplaceAll("GeV","bin #");				
		h->GetXaxis()->SetTitle(axis_title);

	}

}

//----------------------------------------//

void bin_number_y_title(TH1D* h) {

	TString hname = h->GetName();
	if (string(hname).find("Serial") != std::string::npos) {	

		TString axis_title = h->GetYaxis()->GetTitle();
		axis_title.ReplaceAll("deg","bin #");
		axis_title.ReplaceAll("(GeV/c)","bin #");
		axis_title.ReplaceAll("GeV","bin #");				
		h->GetYaxis()->SetTitle(axis_title);

	}

}

//----------------------------------------//

void TV2TH(const TVectorD vec, TH1D* histo)
{
    // Fill vector to histogram,
    for(Int_t i=0; i<vec.GetNrows(); i++)
    {
        histo->SetBinContent(i+1, vec(i));
    }
}

//----------------------------------------//

void TH2TM(const TH2D* histo, TMatrixD& mat, bool rowcolumn) {

    // Fill 2D histogram into matrix
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE

    for (Int_t i=0; i<histo->GetNbinsX(); i++) {

        for (Int_t j=0; j<histo->GetNbinsY(); j++) {

            if (rowcolumn) { mat(i, j) = histo->GetBinContent(i+1, j+1); }
            else { mat(j, i) = histo->GetBinContent(i+1, j+1); }

        }

    }

}

//----------------------------------------//

void TH2TV(const TH1D* histo, TVectorD& vec)
{
    // Fill 1D histogram into matrix
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        vec(i) = histo->GetBinContent(i+1);
    }
}

//----------------------------------------//

double FindTwoDimHistoMaxValue(TH2D* h){

	int NXBins = h->GetXaxis()->GetNbins();
	int NYBins = h->GetYaxis()->GetNbins();	
	double StartHistoMax = -9999.;
	double HistoMax = StartHistoMax;	

	for (int xbin = 1; xbin <= NXBins; xbin++) {

		for (int ybin = 1; ybin <= NYBins; ybin++) {		

			double LocalMax = h->GetBinContent(xbin,ybin);
			if (LocalMax > HistoMax && !isinf(LocalMax) ) { HistoMax = LocalMax; }

		}

	}

	if (HistoMax == StartHistoMax) { cout << "HistoMax = " << HistoMax << endl; }

	return HistoMax;

}

//----------------------------------------//

double FindTwoDimHistoMinValue(TH2D* h){

	int NXBins = h->GetXaxis()->GetNbins();
	int NYBins = h->GetYaxis()->GetNbins();	
	double StartHistoMin = 999999.;
	double HistoMin = StartHistoMin;	

	for (int xbin = 1; xbin<= NXBins; xbin++) {

		for (int ybin = 1; ybin<= NYBins; ybin++) {		

			double LocalMin = h->GetBinContent(xbin,ybin);
			if (LocalMin < HistoMin) { HistoMin = LocalMin; }

		}

	}

	if (HistoMin == StartHistoMin) { cout << "HistoMin = " << HistoMin << endl; }
	return HistoMin;

}

//----------------------------------------//                                                                                               

TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {

  TH1D* TrueClone = (TH1D*)(True->Clone());

  int XBins = SmearMatrix->GetXaxis()->GetNbins();
  int YBins = SmearMatrix->GetYaxis()->GetNbins();

  if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

  TVectorD signal(XBins);
  TMatrixD response(XBins,XBins);

  TH2TV(TrueClone, signal);
  TH2TM(SmearMatrix, response, kTRUE);

  TVectorD RecoSpace = response * signal;
  TV2TH(RecoSpace, TrueClone);

  return TrueClone;

}

//----------------------------------------//

double back_proj_dist (TVector3 shower_v) {

	double backwards_projected_dist = -99999.0;

    double reco_shower_momentum_perp = shower_v.Pt();
    double shower_theta = shower_v.Theta() * (180. / TMath::Pi());
    double shower_phis = shower_v.Phi() * (180. / TMath::Pi());
    double shower_momentum_total_3d = shower_v.Mag();
	
    std::vector<double> shower_unit_vector_3d = {shower_v.X() / shower_momentum_total_3d,
                                 				 shower_v.Y() / shower_momentum_total_3d,
                                 				 shower_v.Z() / shower_momentum_total_3d};
    double center_x = 130.;
    double center_y = 0.;
    double center_z = 525.;
        /*double towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<double> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        double inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        double shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<double> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<double> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        double inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);
*/

    double min_backwards_projected_dist = 1e9;
/*
        //projecting to x walls
        if (shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
              min_backwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2];
        }
        if (isinf(min_backwards_projected_dist)) min_backwards_projected_dist = -99999.0;
*/
    backwards_projected_dist = min_backwards_projected_dist;

	return backwards_projected_dist;

}

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

TString pot_sci_notation(TString run, int decimals = 2) {

    double pot = PeLEE_ReturnBeamOnRunPOT(run);

    if (pot <= 0) return "Invalid POT";

    // Extract exponent
    int exponent = (int)TMath::Floor(TMath::Log10(pot));

    // Extract coefficient (mantissa)
    double coefficient = pot / TMath::Power(10, exponent);

    // Format coefficient with desired precision
    TString coeffFormat = TString::Format("%%.%df", decimals);
    TString coeffStr = TString::Format(coeffFormat, coefficient);

    // Return LaTeX-formatted string: e.g., "6.60 \\times 10^{20}"
    TString result = TString::Format("%s\\times 10^{%d}", coeffStr.Data(), exponent);

    return result;
	
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

double get_median_vector(std::vector<double> v) { 

	size_t n = v.size();

	if (n == 0) { return 0.; }

	std::nth_element(v.begin(), v.begin() + n / 2, v.end());
	return v[n / 2];

}

//----------------------------------------//

double cos_alpha(const TVector3& vertex1, const TVector3& direction, const TVector3& vertex2) {

    // Compute the vector from vertex1 to vertex2
    TVector3 displacement = vertex2 - vertex1;

    // Normalize the direction vector
    TVector3 dirNorm = direction.Unit();

    // Compute the cosine of the angle between the direction and displacement vectors
    double costheta = dirNorm.Dot(displacement.Unit());

    return costheta;

}

//----------------------------------------//

//Backtracked proton-like blips 
bool IsBackTrackedBlip( double d, double cosine_angle){
 // angle w.r.t shower direction and L~24cm radiation/convertion length
 // d < 70 & cosine_angle < -0.9975    , ~Cos(176°) & 3L 
 // d < 48 & cosine_angle < -0.9900    , ~Cos(172°) & 2L
 // d < 30 & cosine_angle < -0.9825    , ~Cos(169°) & 1L

  return ( (d < 70 && cosine_angle < -0.9975 ) ||   (d < 48 && cosine_angle < -0.9900 ) || (d < 30 && cosine_angle < -0.9825  ) );

}

//----------------------------------------//

// electron-like blips
bool IsWithinSphereOutsideConic(const TVector3& vertex1, const TVector3& direction, const TVector3& vertex2, double R) {

    // Compute the vector from vertex1 to vertex2
    TVector3 displacement = vertex2 - vertex1;

    // Compute the squared distance (more efficient than sqrt)
    double distance = displacement.Mag();
    double distanceSquared = displacement.Mag2();

    // Check if within the sphere of radius R
    if (distanceSquared > R * R) return false;

    // Normalize the direction vector
    TVector3 dirNorm = direction.Unit();

    // Compute the cosine of the angle between the direction and displacement vectors
    double cosTheta = dirNorm.Dot(displacement.Unit());

    // The cosine of 45 degrees is sqrt(2)/2 ≈ 0.707
    const double cos45 = std::cos(TMath::Pi() / 4);

    if ( IsBackTrackedBlip(distance, cosTheta) ) { return 0; }

    // Check if the point is inside the sphere but outside the conic region or proton-like backtracked region
    return cosTheta < cos45;

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


