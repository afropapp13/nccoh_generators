#define analyzer_cxx
#include "analyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>

#include "helper_functions.cxx"
#include "constants.h"

using namespace std;
using namespace constants;

//----------------------------------------//

void analyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//----------------------------------------//	

    // output root files

	TString FileNameAndPath = "output_files/" + fweights + "analyzer_ncpi0_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	// Plot declaration

	// Post FSI

	TH1D* TrueSingleBinPlot[NInte];	
	TH1D* TruePi0CosThetaPlot[NInte];
	TH1D* TruePi0MomentumPlot[NInte];
	TH1D* TruePi0InvMassPlot[NInte];

	//----------------------------------------//

	// Pre FSI

	TH1D* NoFSITrueSingleBinPlot[NInte];
	TH1D* NoFSITruePi0CosThetaPlot[NInte];
	TH1D* NoFSITruePi0MomentumPlot[NInte];
	TH1D* NoFSITruePi0InvMassPlot[NInte];

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//--------------------------------------------------//

		// Post FSI

		TrueSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSingleBinPlot", LabelXAxisSingleBin,NBinsSingleBin,ArrayNBinsSingleBin);	
		TruePi0CosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePi0CosThetaPlot", LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TruePi0MomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePi0MomentumPlot", LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TruePi0InvMassPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePi0InvMassPlot", LabelXAxisPi0InvMass,NBinsPi0InvMass,ArrayNBinsPi0InvMass);

		//--------------------------------------------------//

		// Pre FSI

		NoFSITrueSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueSingleBinPlot",";cos#theta_{#mu}",NBinsSingleBin,ArrayNBinsSingleBin);
		NoFSITruePi0CosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITruePi0CosThetaPlot",";cos#theta_{p}",NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		NoFSITruePi0MomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITruePi0MomentumPlot",";#pi^{0} momentum [GeV/c]",NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		NoFSITruePi0InvMassPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITruePi0InvMassPlot", LabelXAxisPi0InvMass,NBinsPi0InvMass,ArrayNBinsPi0InvMass);

		//--------------------------------------------------//

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;

	//----------------------------------------//

	TFile* fweights_file = nullptr;
	TTree* tweights = nullptr;
	float cv_weight = -99.;

	if (fweights == "Weights") {

		if (fOutputFile == "GENIE_v3_0_6") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_Nominal.root"); }

		tweights = (TTree*)fweights_file->Get("GenericVectors__VARS");
		tweights->SetBranchAddress("Weight", &cv_weight);

	}

	//----------------------------------------//
	
	// Loop over the events

	cout << "nentries = " << nentries << endl;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//----------------------------------------//	
	
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%100000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		//----------------------------------------//	

		double t2kweight = 1.;

		if (fweights == "Weights") {

			tweights->GetEntry(jentry); t2kweight = cv_weight;

		}

		//----------------------------------------//			

		double weight = fScaleFactor * Units * A * Weight * t2kweight;	

		if (fOutputFile == "GiBUU_2023") { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 105 samples	
		if (fOutputFile == "GiBUU_2023_medium") { weight = weight/150.; } // To increase the stats, the GiBUU sample has been produced in 73 samples	
		if (fOutputFile == "ACHILLES") { weight = weight*1000./(40./12.); } // ACHILLES scaling still under discussion

		//----------------------------------------//	

		// Signal definition

		// make sure that we don't have a lepton in the final state
		if ( PDGLep == 11 || PDGLep == 13 || PDGLep == 15) { continue; }
		// ACHILLES doesn't know how to handle the cc/nc branch yet
		if (fOutputFile != "ACHILLES") {
	    
			// make sure that we have only NC interactions
			if (cc != 0) { continue; }		
	  
		}

		int ProtonTagging = 0, ChargedPionTagging = 0, Pi0Tagging = 0;
		int heavy_meason_tagging = 0, SigmaTagging = 0, LambdaTagging = 0;
		int PhotonTagging = 0, LeptonTagging = 0 , cluster_tagging = 0;
		int neutron_tagging = 0;
		vector <int> Pi0ID; Pi0ID.clear();		

		int NoFSIProtonTagging = 0, NoFSIChargedPionTagging = 0, NoFSIPi0Tagging = 0;
		int NoFSIheavy_meason_tagging = 0, NoFSISigmaTagging = 0, NoFSILambdaTagging = 0;
		int NoFSIPhotonTagging = 0, NoFSILeptonTagging = 0, nofsi_cluster_tagging = 0;
		int nofsi_neutron_tagging = 0;		
		vector <int> NoFSIPi0ID; NoFSIPi0ID.clear();			

		//----------------------------------------//	

		// Loop over the final state particles / post FSI

	  	for (int i = 0; i < nfsp; i++) {
				
			if (pdg[i] == ProtonPdg ) {

				double ke = E[i] - ProtonMass_GeV;

				// proton kinetic energy threshold
				if ( ke > proton_ke_thres ) {

					ProtonTagging ++;

				}

			}

			else if ( fabs(pdg[i]) == AbsChargedPionPdg )  {

			ChargedPionTagging ++;

			}

			else if ( fabs(pdg[i]) == NeutralPionPdg)  {

			Pi0Tagging ++;
			Pi0ID.push_back(i);

			}

			else if ( fabs(pdg[i]) == KaonPdg || fabs(pdg[i]) == NeutralKaonPdg 
			       || fabs(pdg[i]) == rho_pdg || fabs(pdg[i]) == charged_rho_pdg || fabs(pdg[i]) == eta_pdg)  {

				heavy_meason_tagging ++;
	
			}	
			
			else if ( fabs(pdg[i]) == SigmaPlusPdg || fabs(pdg[i]) == SigmaMinusPdg || fabs(pdg[i]) == NeutralSigmaPdg)  {

				SigmaTagging ++;
	
			}		
			
			else if ( fabs(pdg[i]) == LambdaPdg)  {

				LambdaTagging ++;
	
			}
			
			else if ( fabs(pdg[i]) == PhotonPdg)  {

				PhotonTagging ++;
	
			}
			
			else if ( fabs(pdg[i]) == ElectronPdg || fabs(pdg[i]) == MuonPdg)  {

				LeptonTagging ++;
	
			}			

			else if ( fabs(pdg[i]) == hydrogen_cluster_pdg || fabs(pdg[i]) == nucleon_pair
			       || fabs(pdg[i]) == ArgonPdg || fabs(pdg[i]) == neutron_pair || fabs(pdg[i]) == proton_pair) {

				// Ignore neutrons, numus, nues

			}

			else if ( fabs(pdg[i]) == NuMuPdg || fabs(pdg[i]) == nue_pdg) {

				// Ignore numus, nues

			}

			else if ( fabs(pdg[i]) == NeutronPdg) {

				double ke = E[i] - NeutronMass_GeV;

				// neutron kinetic energy threshold
				if ( ke > neutron_ke_thres ) {

					neutron_tagging ++;

				}

			}			

			else { cout << "post fsi pdg " << pdg[i] << endl; }

	  	} // End of the loop over the final state particles / post FSI

		//----------------------------------------//	

		// Loop over final state particles / pre FSI

		for (int i = 0; i < nvertp; i++) {
			
			if (pdg_vert[i] == ProtonPdg ) {

				double ke = E_vert[i] - ProtonMass_GeV;

				// proton kinetic energy threshold
				if ( ke > proton_ke_thres ) {

					NoFSIProtonTagging ++;

				}

			}

			else if ( fabs(pdg_vert[i]) == AbsChargedPionPdg )  {

				NoFSIChargedPionTagging ++;

			}

			else if ( fabs(pdg_vert[i]) == NeutralPionPdg)  {

				NoFSIPi0Tagging ++;
				NoFSIPi0ID.push_back(i);

			}

			else if ( fabs(pdg_vert[i]) == KaonPdg || fabs(pdg_vert[i]) == NeutralKaonPdg 
			       || fabs(pdg_vert[i]) == rho_pdg || fabs(pdg_vert[i]) == charged_rho_pdg || fabs(pdg_vert[i]) == eta_pdg)  {

				NoFSIheavy_meason_tagging ++;
		
			}	
				
			else if ( fabs(pdg_vert[i]) == SigmaPlusPdg || fabs(pdg_vert[i]) == SigmaMinusPdg || fabs(pdg_vert[i]) == NeutralSigmaPdg)  {

				NoFSISigmaTagging ++;
		
			}		
				
			else if ( fabs(pdg_vert[i]) == LambdaPdg)  {

				NoFSILambdaTagging ++;
		
			}
				
			else if ( fabs(pdg_vert[i]) == PhotonPdg)  {

				NoFSIPhotonTagging ++;
		
			}
				
			else if ( fabs(pdg_vert[i]) == ElectronPdg || fabs(pdg_vert[i]) == MuonPdg)  {

				NoFSILeptonTagging ++;
		
			}			

			else if ( fabs(pdg_vert[i]) == hydrogen_cluster_pdg || fabs(pdg_vert[i]) == nucleon_pair
			       || fabs(pdg_vert[i]) == ArgonPdg || fabs(pdg_vert[i]) == neutron_pair || fabs(pdg_vert[i]) == proton_pair) {

				// Ignore neutrons, numus, nues

			}

			else if ( fabs(pdg_vert[i]) == NuMuPdg 
			       || fabs(pdg_vert[i]) == nue_pdg) {

				// Ignore neutrons, numus, nues

			}

			else if ( fabs(pdg_vert[i]) == NeutronPdg) {

				double ke = E_vert[i] - NeutronMass_GeV;

				// neutron kinetic energy threshold
				if ( ke > neutron_ke_thres ) {

					nofsi_neutron_tagging ++;

				}

			}	

			else { cout << "pre fsi pdg_vert " << pdg_vert[i] << endl; }

		} // End of the loop over the final state particles / pre FSI

	  //----------------------------------------//	

	  // Classify the events based on the interaction type

	  // https://arxiv.org/pdf/2106.15809.pdf

	  int genie_mode = -1.;

	  if (fOutputFile ==  "ACHILLES") {

		// ACHILLES has limited interactions
	    genie_mode = 1;

	  } else {

	    if (TMath::Abs(Mode) == 51 ||TMath::Abs(Mode) == 52) { genie_mode = 1; } // QE
	    else if (TMath::Abs(Mode) == 2) { genie_mode = 2; } // MEC
	    else if ( TMath::Abs(Mode) == 31 || TMath::Abs(Mode) == 32 
			   || TMath::Abs(Mode) == 33 || TMath::Abs(Mode) == 34 
			   || TMath::Abs(Mode) == 42 || TMath::Abs(Mode) == 43
			   || TMath::Abs(Mode) == 44 || TMath::Abs(Mode) == 45
		   ) { genie_mode = 3; } // RES
	    else if ( TMath::Abs(Mode) == 41 || TMath::Abs(Mode) == 46) { genie_mode = 4; } // DIS
	    else if (TMath::Abs(Mode) == 36) { genie_mode = 5;} // COH
	    else { cout << "new interaction " << Mode << endl;/*continue;*/ }  

	  }

	  //----------------------------------------//	

	  // If the signal definition post-FSI  is satisfied
	  if ( Pi0Tagging == 1 && ProtonTagging == 0 && ChargedPionTagging == 0 && 
		   heavy_meason_tagging == 0 && LambdaTagging == 0 && SigmaTagging == 0 &&
		   PhotonTagging == 0 && LeptonTagging == 0 && cluster_tagging == 0 && neutron_tagging == 0
		) { 

	    CounterEventsPassedSelection++;

	    // Kinematics of eutral pion in the final state

	    TLorentzVector Pi04Vector(px[Pi0ID[0]], py[Pi0ID[0]], pz[Pi0ID[0]], E[Pi0ID[0]]);
	    double Pi0InvMass = Pi04Vector.M();

	    //----------------------------------------//

	    // Variables of interest
	    // Assign twice to keep track of the old values as well

	    double Pi0Momentum = Pi04Vector.Rho();
	    double Pi0CosTheta = Pi04Vector.CosTheta();

		if (Pi0CosTheta < pi0_costheta_thres) { continue; }
	  
		//----------------------------------------//	

	    // Underflow / overflow
            
		if (Pi0Momentum < ArrayNBinsPi0Momentum[0]) { Pi0Momentum = (ArrayNBinsPi0Momentum[0] + ArrayNBinsPi0Momentum[1])/2.; }
        if (Pi0Momentum > ArrayNBinsPi0Momentum[NBinsPi0Momentum]) { Pi0Momentum = (ArrayNBinsPi0Momentum[NBinsPi0Momentum] + ArrayNBinsPi0Momentum[NBinsPi0Momentum-1])/2.; }

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    TrueSingleBinPlot[0]->Fill(0.5,weight);	
	    TruePi0CosThetaPlot[0]->Fill(Pi0CosTheta,weight);
	    TruePi0MomentumPlot[0]->Fill(Pi0Momentum,weight);	
	    TruePi0InvMassPlot[0]->Fill(Pi0InvMass,weight);	

	    //----------------------------------------//		

	    // filling in the histo based on the interaction mode

	    TrueSingleBinPlot[genie_mode]->Fill(0.5,weight);	
	    TruePi0CosThetaPlot[genie_mode]->Fill(Pi0CosTheta,weight);
	    TruePi0MomentumPlot[genie_mode]->Fill(Pi0Momentum,weight);
	    TruePi0InvMassPlot[genie_mode]->Fill(Pi0InvMass,weight);

	  } // End of the post-FSI selection

	  //----------------------------------------//

	  // If the signal definition pre-FSI is satisfied
	  if ( NoFSIPi0Tagging == 1 && NoFSIProtonTagging == 0 && NoFSIChargedPionTagging == 0 && 
		   NoFSIheavy_meason_tagging == 0 && NoFSILambdaTagging == 0 && NoFSISigmaTagging == 0 &&
		   NoFSIPhotonTagging == 0 && NoFSILeptonTagging == 0 && nofsi_cluster_tagging == 0 && nofsi_neutron_tagging == 0
	    ) { 

	    // Kinematics of neutral pion in the final state pre FSI

	    TLorentzVector Pi04Vector(px_vert[NoFSIPi0ID[0]], py_vert[NoFSIPi0ID[0]], pz_vert[NoFSIPi0ID[0]], E_vert[NoFSIPi0ID[0]]);
	    double Pi0InvMass = Pi04Vector.M();

	    //----------------------------------------//

	    double Pi0CosTheta = Pi04Vector.CosTheta();
	    double Pi0Momentum = Pi04Vector.Rho();

		if (Pi0CosTheta < pi0_costheta_thres) { continue; }		
	
	    //----------------------------------------//	

	    // Underflow / overflow

		if (Pi0Momentum < ArrayNBinsPi0Momentum[0]) { Pi0Momentum = (ArrayNBinsPi0Momentum[0] + ArrayNBinsPi0Momentum[1])/2.; }
        if (Pi0Momentum > ArrayNBinsPi0Momentum[NBinsPi0Momentum]) { Pi0Momentum = (ArrayNBinsPi0Momentum[NBinsPi0Momentum] + ArrayNBinsPi0Momentum[NBinsPi0Momentum-1])/2.; }

	    //----------------------------------------//	

	    // filling in the histo regardless of interaction mode

	    NoFSITrueSingleBinPlot[0]->Fill(0.5,weight);
	    NoFSITruePi0CosThetaPlot[0]->Fill(Pi0CosTheta,weight);
	    NoFSITruePi0MomentumPlot[0]->Fill(Pi0Momentum,weight);							
	    NoFSITruePi0InvMassPlot[0]->Fill(Pi0InvMass,weight);							

	    NoFSITrueSingleBinPlot[genie_mode]->Fill(0.5,weight);
	    NoFSITruePi0CosThetaPlot[genie_mode]->Fill(Pi0CosTheta,weight);
	    NoFSITruePi0MomentumPlot[genie_mode]->Fill(Pi0Momentum,weight);							
	    NoFSITruePi0InvMassPlot[genie_mode]->Fill(Pi0InvMass,weight);							



	  } // End of the post-FSI selection

	  //----------------------------------------//
	
	} // End of the loop over the events

	//----------------------------------------//	

	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	//----------------------------------------//	
	//----------------------------------------//	

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		divide_bin_width(NoFSITruePi0CosThetaPlot[inte]);
		divide_bin_width(NoFSITruePi0MomentumPlot[inte]);

		divide_bin_width(TruePi0CosThetaPlot[inte]);
		divide_bin_width(TruePi0MomentumPlot[inte]);		

	} // End of the loop over the interaction processes		

	//----------------------------------------//		

	file->cd();
	file->Write();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created" << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	//----------------------------------------//		

} // End of the program
