
	// ROOT stuff
	#include <TApplication.h>
	#include <TCanvas.h>
	#include <TFile.h>
	#include <TRandom3.h>
	#include <TH1F.h>
	#include <TH2F.h>
	#include <TGraph.h>
	#include <TTree.h>
	#include <TROOT.h>
	#include <TMinuit.h>
	#include <TNtuple.h>
	#include <TComplex.h>
	#include <TFile.h>
	#include <TStyle.h>
	#include <TH2Poly.h>
	#include <TGraphErrors.h>

	//Minuit
	#include <Minuit2/MnStrategy.h>
	#include <Minuit2/Minuit2Minimizer.h>

	// System stuff
	#include <CLI/Timer.hpp>
	#include <fstream>
	#include <string>
	#include <time.h>

	// GooFit stuff
	#include <goofit/Application.h>
	#include <goofit/BinnedDataSet.h>
	#include <goofit/FitManager.h>
	#include <goofit/PDFs/GooPdf.h>
	#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
	#include <goofit/PDFs/basic/PolynomialPdf.h>
	#include <goofit/PDFs/combine/AddPdf.h>
	#include <goofit/PDFs/combine/ProdPdf.h>
	#include <goofit/PDFs/physics/Amp3Body.h>
	#include <goofit/PDFs/physics/DalitzVetoPdf.h>
	#include <goofit/PDFs/physics/ResonancePdf.h>
	#include <goofit/UnbinnedDataSet.h>
	#include <goofit/Variable.h>
	#include <goofit/PDFs/physics/DalitzPlotter.h>
	#include <goofit/FunctorWriter.h>

	//Matrix 
	#include <Eigen/Core>
	#include <Eigen/LU>

	#define torad(x)(x*M_PI/180.)

	using namespace std;
	using namespace GooFit;
	using namespace ROOT;

	//D_MASS (PDG)
	double D_MASS   = 1.86966;
	double pi_MASS = 0.13957018;

	double Decay_MASS = D_MASS;
	double d1_MASS  = pi_MASS; 
	double d2_MASS  = pi_MASS;
	double d3_MASS  = pi_MASS;

	Variable Mother_Mass("DecayParticle_Mass",Decay_MASS);
	Variable Daughter1_Mass("DecayProduct_1_Mass",d1_MASS);
	Variable Daughter2_Mass("DecayProduct_2_Mass",d2_MASS);
	Variable Daughter3_Mass("DecayProduct_3_Mass",d3_MASS);

	//Bins for grid normalization
	const int bins = 500;

	//Dalitz Limits
	const fptype s12_min = (d1_MASS  + d2_MASS)*(d1_MASS  + d2_MASS);
	const fptype s12_max = (Decay_MASS   - d3_MASS)*(Decay_MASS   - d3_MASS);
	const fptype s13_min = (d1_MASS  + d3_MASS)*(d1_MASS  + d3_MASS);
	const fptype s13_max = (Decay_MASS   - d2_MASS)*(Decay_MASS   - d2_MASS);
	const fptype s23_min = (d2_MASS  + d3_MASS)*(d2_MASS  + d3_MASS);
	const fptype s23_max = (Decay_MASS   - d1_MASS)*(Decay_MASS   - d1_MASS);

	//Observables
	Observable s12("s12",s12_min,s12_max); 
	Observable s13("s13",s13_min,s13_max);
	Observable s23("s23",s23_min,s23_max);
	EventNumber eventNumber("eventNumber"); 

	GooPdf *polyEff( Observable s12 , Observable s13){
	    vector<Variable> offsets;
	    vector<Observable> observables;
	    vector<Variable> coefficients;
	    Variable constantOne("c1", 1);
	    Variable constantZero("c0", 0);
	    observables.push_back(s12);
	    observables.push_back(s13);
	    offsets.push_back(constantZero);
	    offsets.push_back(constantZero);
	    coefficients.push_back(constantOne);
	    PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); //No efficiency

	    return eff;
	}

	Amp3Body* makesignalpdf( Observable s12, Observable s13, EventNumber eventNumber, GooPdf* eff = 0){

	    //set up the decay channel
	    DecayInfo3 dtoppp;
	    dtoppp.motherMass   = D_MASS;
	    dtoppp.daug1Mass    = d1_MASS;
	    dtoppp.daug2Mass    = d2_MASS;
	    dtoppp.daug3Mass    = d3_MASS;
	    dtoppp.meson_radius = 1.5; //GeV^-1
        dtoppp.mother_meson_radius = 5.0; //GeV^-1

	   //Mass and width

		//parameters from Laura++
		double f0_980_MASS    = 0.965;
		double f0_980_GPP     = 0.165;
		double f0_980_GKK     = 4.21;
		double f0_980_amp     = 3.2;;
		double f0_980_img     = -1.35;;

		Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS);
		Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP);
		Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK);
		Variable v_f0_980_real("f0_980_REAL",f0_980_amp,0.01,0,0);
		Variable v_f0_980_img("f0_980_IMAG",f0_980_img,0.01,0,0);

		auto f0_980 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

		double f0_1370_MASS   = 1.259;
		double f0_1370_WIDTH  = 0.298;
		double f0_1370_amp    = 0.04408;
		double f0_1370_img  = -0.414616;

		Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS);
		Variable v_f0_1370_Width("f0_1370_WIDTH",f0_1370_WIDTH);
		Variable v_f0_1370_real("f0_1370_REAL",f0_1370_amp,0.01,0,0);
		Variable v_f0_1370_img("f0_1370_IMAG",f0_1370_img,0.01,0,0);

		auto f0_1370 = new Resonances::RBW("f0_1370",v_f0_1370_real,v_f0_1370_img,v_f0_1370_Mass,v_f0_1370_Width,0,PAIR_12,true,true);


		double f0_1500_MASS   = 1.506;
		double f0_1500_WIDTH  = 0.112;
		double f0_1500_amp    = 1.97753;
		double f0_1500_img  = 0.207441;

		Variable v_f0_1500_Mass("f0_1500_MASS",f0_1500_MASS);
		Variable v_f0_1500_Width("f0_1500_WIDTH",f0_1500_WIDTH);
		Variable v_f0_1500_real("f0_1500_REAL",f0_1500_amp,0.01,0,0);
		Variable v_f0_1500_img("f0_1500_IMAG",f0_1500_img,0.01,0,0);

		auto f0_1500 = new Resonances::RBW("f0_1500",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,0,PAIR_12,true,true);

		double f0_1710_MASS   = 1.704;
		double f0_1710_WIDTH  = 0.123;
		double f0_1710_amp    = -0.161453;
		double f0_1710_img  = 0.8019;

		Variable v_f0_1710_Mass("f0_1710_MASS",f0_1710_MASS);
		Variable v_f0_1710_Width("f0_1710_WIDTH",f0_1710_WIDTH);
		Variable v_f0_1710_real("f0_1710_REAL",f0_1710_amp,0.01,0,0);
		Variable v_f0_1710_img("f0_1710_IMAG",f0_1710_img,0.01,0,0);

		auto f0_1710 = new Resonances::RBW("f0_1710",v_f0_1710_real,v_f0_1710_img,v_f0_1710_Mass,v_f0_1710_Width,0,PAIR_12,true,true);


		double sigma_MASS   = 0.450;
		double sigma_WIDTH  = 0.275;
		double sigma_amp    = 3.48;;
		double sigma_img  = 2.0;;

		Variable v_sigma_Mass("sigma_MASS",sigma_MASS);
		Variable v_sigma_Width("sigma_WIDTH",sigma_WIDTH);
		Variable v_sigma_real("sigma_REAL",sigma_amp,0.01,0,0);
		Variable v_sigma_img("sigma_IMAG",sigma_img,0.01,0,0);

		auto sigma = new Resonances::RBW("sigma",v_sigma_real,v_sigma_img,v_sigma_Mass,v_sigma_Width,0,PAIR_12,true,true);


		double omega_MASS   = 0.78265;
		double omega_WIDTH  = 0.00849;
		double omega_amp    = 0.00130426;
		double omega_img  = -0.020938;

		Variable v_omega_Mass("omega_MASS",omega_MASS);
		Variable v_omega_Width("omega_WIDTH",omega_WIDTH);
		Variable v_omega_real("omega_REAL",omega_amp,0.01,0,0);
		Variable v_omega_img("omega_IMAG",omega_img,0.01,0,0);

		auto omega = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true,true);

		double rho770_MASS   = 0.77549;
		double rho770_WIDTH  = 0.1491;
		double rho770_amp    = 1.6;
		double rho770_img  =  0.;

		Variable v_rho770_Mass("rho770_MASS",rho770_MASS);
		Variable v_rho770_Width("rho770_WIDTH",rho770_WIDTH);
		Variable v_rho770_real("rho770_REAL",rho770_amp,0.01,0,0);
		Variable v_rho770_img("rho770_IMAG",rho770_img,0.01,0,0);

		auto rho770 = new Resonances::GS("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true,true);

		double rho1450_MASS = 1.465;
		double rho1450_WIDTH = 0.4;
		double rho1450_amp    = 2.3;
		double rho1450_img  =  1.7;

		Variable v_rho1450_Mass("rho1450_MASS",rho1450_MASS);
		Variable v_rho1450_Width("rho1450_WIDTH",rho1450_WIDTH);
		Variable v_rho1450_real("rho1450_REAL",rho1450_amp,0.01,0,0);
		Variable v_rho1450_img("rho1450_IMAG",rho1450_img,0.01,0,0);

		auto rho1450 = new Resonances::GS("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true,true);

		double f2_1270_MASS     = 1.2751; 
		double f2_1270_WIDTH    = 0.1851;
		double f2_1270_amp      = 2.262;
		double f2_1270_img    = -2.82;

		Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS);
		Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH);
		Variable v_f2_1270_real("f2_1270_REAL",f2_1270_amp);
		Variable v_f2_1270_img("f2_1270_IMAG",f2_1270_img);

		auto f2_1270 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true,true);

		//If you want include a resonance in your model, just push into the vector 'vec_resonances'
		std::vector<ResonancePdf *> vec_resonances;

		vec_resonances.push_back(sigma); 
		vec_resonances.push_back(f0_980);
		vec_resonances.push_back(f0_1370);
		vec_resonances.push_back(f0_1500);
		vec_resonances.push_back(f0_1710);
		vec_resonances.push_back(omega);
		vec_resonances.push_back(rho770); 
		vec_resonances.push_back(rho1450);
		vec_resonances.push_back(f2_1270);

		dtoppp.resonances = vec_resonances;

		if(!eff)
			eff = polyEff(s12,s13);

	return new Amp3Body("signalPDF", s12, s13, eventNumber, dtoppp, eff);
}

void getData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool toy) {
    //load data in a GooFit::dataset

    toyFileName = app.get_filename(toyFileName, "MC/");

    auto obs               = data.getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    auto openRoot   = new TFile(toyFileName.c_str());
    auto tree       = (TTree*)openRoot->Get("DecayTree");
    auto s12_val(0.);
    auto s13_val(0.);

    if(toy){
        tree->SetBranchAddress("s12",&s12_val);
        tree->SetBranchAddress("s13",&s13_val);
    }else{
        
        tree->SetBranchAddress("s12_pipi_DTF",&s12_val);
        tree->SetBranchAddress("s13_pipi_DTF",&s13_val);
    }

    size_t j = 0;
    for(size_t i = 0 ; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        s12.setValue(s12_val);
        s13.setValue(s13_val);
        eventNumber.setValue(data.getNumEvents());
        if((s12.getValue()<s12.getUpperLimit())
            &&(s13.getValue()<s13.getUpperLimit())
            &&(s12.getValue()>s12.getLowerLimit())
            &&(s13.getValue()>s13.getLowerLimit()))
        {
            data.addEvent();
            if(j<10) printf("[%d] = (%f , %f)\n",i,s12.getValue(),s13.getValue());
            j++;
        }
    }


}

void to_root(UnbinnedDataSet& toyMC , std::string name ){
//save GooFit::Dataset in a root file
    auto obs               = toyMC.getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    double _s12, _s13,_s23;
    auto f = new TFile(name.c_str(),"recreate");
    auto t = new TTree("DecayTree","toyMC");
    auto b_s12 = t->Branch("s12",&_s12,"s12/D");
    auto b_s13 = t->Branch("s13",&_s13,"s13/D");
    auto b_s23 = t->Branch("s23",&_s23,"s23/D");

   
    for(int i = 0; i < toyMC.getNumEvents(); i++){
		toyMC.loadEvent(i);
		t->GetEntry(i);
		_s12 = s12.getValue();
		_s13 = s13.getValue();
		_s23 = POW2(Decay_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s12.getValue() - s13.getValue() ;
		t->Fill();
    }
	t->Write("",TObject::kOverwrite);
	f->Write();
	f->Close();
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

}

Amp3Body* runFit(GooPdf *totalPdf, Amp3Body *signal, UnbinnedDataSet *data, std::string name) {
//This function run the data fit

    //Setting data and EventNumber
    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());

	auto output = fmt::format("Fit/{0}/fit_result_before.txt",name.c_str());
    writeToFile(totalPdf, output.c_str());

    //Fitter (it uses ROOT::FunctionMinimum API)
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(2);
    datapdf.setMaxCalls(200000);

    //Start fit
    auto func_min = datapdf.fit();


	output = fmt::format("Fit/{0}/fit_result_fitted.txt",name.c_str());
    writeToFile(totalPdf, output.c_str());

    return signal;
}


int main(int argc, char **argv){
    
    GooFit::Application app{"Genfit",argc,argv};
    
    std::string input_data_name = "input.root";
    std::string fit_name = "Fit";
	std::string acc_file = "acc_hist_0_Smoothed.root";
    std::string bkg_file = "bkg_hist_0_Smoothed.root";
    std::string toyName = "MC.root";
    bool save_toy = false;
    bool is_toy = false;
    bool no_acc_and_bkg = true;
    size_t Nevents=1000000;

    auto fit = app.add_subcommand("fit","fit data");
    fit->add_option("-f,--file",input_data_name,"name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit") ;
    fit->add_option("-s,--saveToy",save_toy,"save toy in root file");
    fit->add_option("-n,--fitName",fit_name,"name of this fit(useful to save results)")->required(true);
	fit->add_option("-a,--acc",acc_file,"name of acc file");
    fit->add_option("-b,--bkg",bkg_file,"name of bkg file"); 
    fit->add_option("-d,--disable-acc-bkg",no_acc_and_bkg,"disable-acc-bkg"); 
    
    auto makeToy = app.add_subcommand("makeToy","make a toy");
    makeToy->add_option("-e,--nevents",Nevents,"number of events");
    makeToy->add_option("-n,--name",toyName,"output_toy_name.root");
    makeToy->add_option("-s,--saveToy",save_toy,"save toy in root file");	
    makeToy->add_option("-a,--acc",acc_file,"name of acc file");
    makeToy->add_option("-b,--bkg",bkg_file,"name of bkg file");
    makeToy->add_option("-d,--disable-acc-bkg",no_acc_and_bkg,"disable-acc-bkg"); 

    app.require_subcommand();

    GOOFIT_PARSE(app);

    // Make the MC directory if it does not exist
    std::string command = "mkdir -p MC";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    // Make the Fit directory if it does not exist
    command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");

    
    s12.setNumBins(bins);
    s13.setNumBins(bins);
    
    auto efficiency = polyEff(s12,s13);
    auto background = nullptr;
    GooPdf* totalpdf = nullptr;
    
    // to be implemented...
    /*if(!no_acc_and_bkg){
         efficiency = makeHistogramPdf(efffile,effhist,s12,s13,true,false,false);
         background = makeHistogramPdf(bkgfile,bkghist,s12,s13,false,false,false);
         auto signal = makesignalpdf(s12, s13, eventNumber,efficiency);
         totalpdf = new AddPdf("totalpdf", Variable("frac",0.93), signal, background) ;
    
    }*/
   
    auto signal = makesignalpdf(s12, s13, eventNumber,efficiency);
    totalpdf = new ProdPdf("totalpdf", {signal}) ;
   

    if(*makeToy) {
        DalitzPlotter dplotter{totalpdf, signal};
        UnbinnedDataSet data({s12, s13, eventNumber});
        dplotter.fillDataSetMC(data, Nevents);
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout <<   data.getNumEvents() << " events was generated!"  << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        if(save_toy) {
		    auto fullName= fmt::format("MC/{0}",toyName);
        	to_root(data,fullName); 
		    std::cout <<   toyName << " root file was saved in MC folder" << std::endl;
        	std::cout << "----------------------------------------------------------" << std::endl;
        }
		
		std::cout << "Fit Fractions Interference" << '\n';
		
		//signal->normalize();
		auto frac = signal->fit_fractions(true);


		return 0;
    }   


	/* it uses accept and reject method (very inefficient for B phase space)
	if(*makeToy) {
		auto tuple = signal->GenerateSig(Nevents,200323);
		if(save_toy) {
		    auto name= fmt::format("MC/{0}",toyName);
			double _s12, _s13,_s23;

			auto f = new TFile(name.c_str(),"recreate");
			auto t = new TTree("DecayTree","toyMC");
			auto b_s12 = t->Branch("s12",&_s12,"s12/D");
			auto b_s13 = t->Branch("s13",&_s13,"s13/D");
			auto b_s23 = t->Branch("s23",&_s23,"s23/D");

			auto variables = std::get<1>(tuple);
			auto weights   = std::get<2>(tuple);
			auto flags     = std::get<3>(tuple);

			for(int i = 0; i < weights.size(); i++){
				if(flags[i] == 1) {
					t->GetEntry(i);
					_s12 = (*(variables[0]))[i];
					_s13 = (*(variables[2]))[i];
					_s23 = (*(variables[1]))[i];
					t->Fill();
				}
			}
			t->Write("",TObject::kOverwrite);
			f->Write();
			f->Close();
			std::cout << "------------------------------------------" << std::endl;
			std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
			std::cout << "------------------------------------------" << std::endl;
		}
		return 0;
    }*/

    if(*fit){

	    auto command = fmt::format("mkdir -p Fit/{0}",fit_name);
    	if(system(command.c_str()) != 0)
        	throw GooFit::GeneralError("Making directory failed");

        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Reading file --> " << input_data_name       << std::endl; 
        std::cout << "------------------------------------------" << std::endl;
	    UnbinnedDataSet data({s12, s13, eventNumber});
        getData(input_data_name, app, data,is_toy);	
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Num Entries Loaded =  " << data.getNumEvents()       << std::endl; 
        std::cout << "------------------------------------------" << std::endl;
        auto output_signal = runFit(totalpdf,signal, &data,fit_name);

		std::cout << "norm = " << output_signal->normalize() << std::endl;
    }
        
}