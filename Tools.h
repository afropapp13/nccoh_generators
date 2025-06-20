#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TString.h"
#include "TMath.h"
#include <TVector3.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TVectorD.h>
#include <TMatrixD.h>

class Tools {

//	private:


	public:

		// Default constructor
		Tools(){ MuonMass = 106; ProtonMass = 938.272; MuonPdg = 13; ProtonPdg = 2212;
			 FVx = 256.; FVy = 232; FVz = 1037.; borderx = 25.; bordery = 25.; borderz = 25.;}

		// Default destructor
		~Tools(){}

		int ConcatRunSubRunEvent(int run, int subrun, int event, int univ);
		double PoissonRandomNumber(int seed);
		bool is_meson_or_antimeson(int pdg);
		bool IsContained(TVector3 TrackStart, TVector3 TrackEnd);
		bool inFV(double x, double y, double z);
		bool inFVVector(TVector3 vector);
		bool loose_inFVVector(TVector3 vector);		
		double PToKE(int pdg, double momentum);
		double KEToP(int pdg, double ke);
		TString to_string_with_precision(double a_value, const int n);
		TString ConvertToString(double value);
		int ReturnIndex(double value, std::vector<double> vec);
		void Reweight(TH1D* h, double SF);	
		void Reweight2D(TH2D* h, double SF);	
		std::vector<TMatrixD> MatrixDecomp(int nbins,TVectorD matrix_pred,TMatrixD matrix_syst);	
		
		int Return2DNBins(std::vector< std::vector<double> > BinVector);
		std::vector<double> Return2DBinIndices(std::vector< std::vector<double> > BinVector);			
		int ReturnIndexIn2DList(std::vector< std::vector<double> > BinEdgeVector, int SliceIndex, double ValueInSlice);
		
		int Return3DNBins(std::vector< std::vector< std::vector<double> > > BinVector);
		std::vector<double> Return3DBinIndices(std::vector< std::vector< std::vector<double> > > BinVector);			
		int ReturnIndexIn3DList(std::vector< std::vector< std::vector<double> > > BinEdgeVector, int FirstSliceIndex, int SecondSliceIndex, double ValueInSlice);
		
		TH1D* GetHistoBins(TH1D* h,int LowBin,int HighBin,double ScaleFactor,std::vector<double> Binning, TString Name);	
		TH2D* Get2DHistoBins(TH2D* h,int LowBin,int HighBin,double ScaleFactor,std::vector<double> Binning, bool Scale);	
		std::vector< std::vector<double> > CollapseMatrixIntoArray(std::vector< std::vector< std::vector<double> > > Matrix);			

		double CosAlpha(TVector3 blip, TVector3 vertex, TVector3 shower);

		double MuonMass; // MeV
		double ProtonMass; // MeV
		int MuonPdg;
		int ProtonPdg;
		double FVx;
		double FVy;
		double FVz;
		double borderx;
		double bordery;
		double borderz;

};

#endif
