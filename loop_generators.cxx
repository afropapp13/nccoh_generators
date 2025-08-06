{
	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
			
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/geniev3/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6");
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/geniev3/14_1000180400_NC_v3_6_0_AR23_20i_00_000.flat.root"); WhichName.push_back("GENIE_v3_6_0_AR23");	
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/neut/bnb.ub.num.neut_5_4_0_1.flat.root"); WhichName.push_back("NEUT_5_4_0_1");
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/neut/bnb.ub.num.neut_5_6_0.flat.root"); WhichName.push_back("NEUT_5_6_0");	
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/nuwro/nuwro_numu_bnb_nc_25_03_1.flat.root"); WhichName.push_back("NuWro_25_03_1");
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/gibuu/gibuu2025_nc_numu_1000runs.root"); WhichName.push_back("GiBUU_2025");			
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/geniev3_reinsehgal/14_1000180400_NC_v3_6_0_G18_10a_02_11a.flat.root"); WhichName.push_back("GENIE_v3_6_0_RS");
	WhichSample.push_back("/exp/uboone/app/users/apapadop/test/BuildEventGenerators/jobcards/samples/nuwro_delta_res_NC_BNB_Ar/NuWro.flat.root"); WhichName.push_back("NuWro_25_03_1_delta");
	WhichSample.push_back("/exp/uboone/app/users/apapadop/test/tweaked_coh_nuwro/BuildEventGenerators/jobcards/samples/NuWro_MicroBooNE_BNB_numu_NC_tweakedCOH_A_scaling.flat.root"); WhichName.push_back("NuWro_25_03_1_tweaked_coh");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/achilles_samples/nccoh/achilles_nccoh.flat.root"); WhichName.push_back("achilles");	
	
	//----------------------------------------//

	gROOT->ProcessLine(".L analyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("analyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}

};
