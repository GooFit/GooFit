#pragma once

#include <cfloat>
#include <string>

#include <goofit/Log.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

namespace GooFit 
{
    class TestUtils final {
        public:
            static fptype getPDFVal(fpcomplex amp)
            {
                return ( amp*thrust::conj(amp) ).real();
            }

            static void testCompareComplex(
                fpcomplex testVal,
                fpcomplex expectedVal,
                fptype EPS, 
                bool printMe,
                const std::string& valName)
            {
                testCompareFloat(testVal.real(), expectedVal.real(), EPS, printMe, valName+".real()");
                testCompareFloat(testVal.imag(), expectedVal.imag(), EPS, printMe, valName+".imag()");
            }

            static void testCompareFloat(
                fptype testVal, 
                fptype expectedVal,
                fptype EPS, 
                bool printMe,
                const std::string& valName)
            {
                if (printMe)
                {
                    GOOFIT_INFO("{}: Expected {}, actual {}", valName, expectedVal, testVal);
                }

                fptype diff = std::abs(expectedVal - testVal);
                if (testVal == 0.0 || expectedVal == 0.0 || diff < FLT_MIN)
                {
                    REQUIRE_THAT( testVal, Catch::Matchers::WithinAbs(expectedVal, EPS*FLT_MIN) ); 
                }
                else
                {
                    REQUIRE_THAT( testVal, Catch::Matchers::WithinRel(expectedVal, EPS) );
                }
            }
    }; // end class
} // end namespace GooFit