{
	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
			
	WhichSample.push_back("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/geniev3/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6");

	//----------------------------------------//

	gROOT->ProcessLine(".L analyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("analyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}

};
