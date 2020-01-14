#include <goofit/Python.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator_TD.h>
#include <thrust/count.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_Amp4Body_TD(py::module &m) {
    py::class_<Amp4Body_TD, Amp4BodyBase> cls(m, "Amp4Body_TD");
    cls.def(py::init<std::string,
                     std::vector<Observable>,
                     DecayInfo4t,
                     MixingTimeResolution *,
                     GooPdf *,
                     Observable *,
                     unsigned int>(),
            "n"_a,
            "observables"_a,
            "decay"_a,
            "r"_a,
            "eff"_a,
            "mistag"_a       = nullptr,
            "MCeventsNorm"_a = 5e6,
            py::keep_alive<1, 4>(),
            py::keep_alive<1, 5>(),
            py::keep_alive<1, 6>(),
            py::keep_alive<1, 7>())

        .def("set_special_integral", &Amp4Body_TD::set_special_integral,"special"_a)
        .def("setDataSize", &Amp4Body_TD::setDataSize, "dataSize"_a, "evtSize"_a)
        .def("setGenerationOffset", &Amp4Body_TD::setGenerationOffset, "off"_a)
        .def("populateArrays",&Amp4Body_TD::populateArrays)
        .def("normalize",&Amp4Body_TD::normalize)
        .def("setForceIntegrals",&Amp4Body_TD::setForceIntegrals)
        .def("getMCevents",&Amp4Body_TD::getMCevents)
        .def("setMaxWeight",&Amp4Body_TD::setMaxWeight,"wmax"_a)
        .def("set_norm_dtime",[](Amp4Body_TD &self,py::array_t<fptype> pydtime){
	    mcbooster::RealVector_h norm_dtime_h = self.get_norm_dtime();
	    for(int i =0; i < norm_dtime_h.size();i++){
	      norm_dtime_h[i] = pydtime.mutable_at(i);
	    }
	    //copy back over to device
	    self.set_norm_dtime(norm_dtime_h);
	  })
        .def("set_norm_pdf_weights",[](Amp4Body_TD &self,py::array_t<fptype> pyweight){
	    mcbooster::RealVector_h norm_weight_h = self.get_norm_pdf_weights();
	  for(int i =0; i < norm_weight_h.size();i++){
	    norm_weight_h[i] = pyweight.mutable_at(i);
	  }
	  //copy back over to device                                                                                               
	  self.set_norm_pdf_weights(norm_weight_h);
	})
        .def("set_norm_eff",[](Amp4Body_TD &self, py::array_t<fptype> pyeff){
	    mcbooster::RealVector_h norm_eff_h = self.get_norm_eff();
	    for(int i = 0; i < norm_eff_h.size();i++){
	      norm_eff_h[i] = pyeff.mutable_at(i);
	    }
	    self.set_norm_eff(norm_eff_h);
	  })
        .def("set_norm_importance_weights",[](Amp4Body_TD &self,py::array_t<fptype> pyweight){
	    mcbooster::RealVector_h norm_weight_h = self.get_norm_importance_weights();
	    for(int i =0; i < norm_weight_h.size();i++){
	      norm_weight_h[i] = pyweight.mutable_at(i);
	    }
	    //copy back over to device                                                                                                                                                                        
	    self.set_norm_importance_weights(norm_weight_h);
	  })

        .def("get_norm_eff",[](Amp4Body_TD &self){
	    mcbooster::RealVector_h norm_eff = self.get_norm_eff();
	    py::array_t<fptype> pyeff{norm_eff.size()};
	    for(int i = 0 ;i < norm_eff.size();i++){
	      pyeff.mutable_at(i) = norm_eff[i];
	    }
	    return pyeff;
	  })
        .def("get_norm_m12",
	     [](Amp4Body_TD &self){
	       mcbooster::RealVector_h m12 = self.get_norm_m12();
	       py::array_t<fptype> pym12{m12.size()};
	       for(int i = 0; i < m12.size();i++){
		 pym12.mutable_at(i) = m12[i];
	       }
	       return pym12;
	     })
         .def("get_norm_m34",
	     [](Amp4Body_TD &self){
	       mcbooster::RealVector_h m34 = self.get_norm_m34();
	       py::array_t<fptype> pym34{m34.size()};
	       for(int i = 0; i < m34.size();i++){
		 pym34.mutable_at(i) = m34[i];
	       } 
	        return pym34;
	     })

         .def("get_norm_c34",
	      [](Amp4Body_TD &self){
		mcbooster::RealVector_h c34 = self.get_norm_c34();
		py::array_t<fptype> pyc34{c34.size()};
		for(int i = 0; i < c34.size();i++){
		  pyc34.mutable_at(i) = c34[i];
		}  
		return pyc34;
	      })
         .def("get_norm_c12",
	      [](Amp4Body_TD &self){
		mcbooster::RealVector_h c12 = self.get_norm_c12();
		py::array_t<fptype> pyc12{c12.size()};
		for(int i = 0; i < c12.size();i++){
		  pyc12.mutable_at(i) = c12[i];
		}  
		return pyc12;
	      })

         .def("get_norm_phi",
	      [](Amp4Body_TD &self){
		mcbooster::RealVector_h phi = self.get_norm_phi();
		py::array_t<fptype> pyphi{phi.size()};
		for(int i = 0; i < phi.size();i++){
		  pyphi.mutable_at(i) = phi[i];
		}  
		return pyphi;
	      })

          .def("get_norm_dtime",
	       [](Amp4Body_TD &self){
		 mcbooster::RealVector_h dtime = self.get_norm_dtime();
		 py::array_t<fptype> pydtime{dtime.size()};
		 for(int i = 0; i < dtime.size();i++){
		   pydtime.mutable_at(i) = dtime[i];
		 }
		 return pydtime;
	       })
      .def("get_norm_pdf_weights",
	   [](Amp4Body_TD &self){
	     mcbooster::RealVector_h weight = self.get_norm_pdf_weights();
	     py::array_t<fptype> pyweight{weight.size()};
	     for(int i = 0; i < weight.size();i++){
	       pyweight.mutable_at(i) = weight[i];
	     }
	     return pyweight;
	   })
      .def("get_norm_importance_weights",
           [](Amp4Body_TD &self){
	     mcbooster::RealVector_h weight = self.get_norm_importance_weights();
	     py::array_t<fptype> pyweight{weight.size()};
             for(int i = 0; i < weight.size();i++){
               pyweight.mutable_at(i) = weight[i];
             }
             return pyweight;
           })

        .def("GenerateSig",
             [](Amp4Body_TD &self, size_t numEvents) {
                 mcbooster::ParticlesSet_h particles; // vector of pointers to vectors of 4R
                 mcbooster::VariableSet_h variables;  // vector of pointers to vectors of Grealt
                 mcbooster::RealVector_h weights;     // vector of greal t
                 mcbooster::BoolVector_h flags;       // vector of gboolt

                 std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);
		 py::array_t<fptype> pyparticles{{(size_t)4 * 4, numEvents}};
                 py::array_t<fptype> pyvariables{{(size_t)6, numEvents}};
                 py::array_t<fptype> pyweights{numEvents};
                 py::array_t<fptype> pyflags{numEvents};
		 int accepted = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
		 std::cout << "nAcc (from length of weights vector is " << weights.size() << std::endl;
		 std::cout << "Number of accepted flags: " << accepted << std::endl;
                 for(int i = 0; i < 4; i++) {
                     for(int j = 0; j < weights.size(); j++) {
                         pyparticles.mutable_at(i * 4, j)     = (*(particles[i]))[j].get(0);
                         pyparticles.mutable_at(i * 4 + 1, j) = (*(particles[i]))[j].get(1);
                         pyparticles.mutable_at(i * 4 + 2, j) = (*(particles[i]))[j].get(2);
                         pyparticles.mutable_at(i * 4 + 3, j) = (*(particles[i]))[j].get(3);
                     }
                 }
		 
		 
                 for(int i = 0; i < 6; i++) {
                     for(int j = 0; j < weights.size(); j++) {
                         pyvariables.mutable_at(i, j) = (*(variables[i]))[j];
                     }
                 }

                 for(int i = 0; i < weights.size(); i++) {
                     pyweights.mutable_at(i) = weights[i];
                 }

                 for(int i = 0; i < weights.size(); i++) {
                     pyflags.mutable_at(i) = flags[i];
                 }
		 std::cout << "Check to see what weights have been initialized with" << std::endl;
		 for(int i = weights.size() -1; i < weights.size() + 10; i++){
		   std::cout << "Flag value: " << flags[i] << std::endl;
		 }
                 delete variables[0];
                 delete variables[1];
                 delete variables[2];
                 delete variables[3];
                 delete variables[4];
                 delete variables[5];

                 delete particles[0];
                 delete particles[1];
                 delete particles[2];
                 delete particles[3];

                 return std::make_tuple(pyparticles, pyvariables, pyweights, pyflags);
             })
            
        ;

    m.attr("TDDP4") = cls;
}
