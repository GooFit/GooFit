#include <goofit/Catch.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

using namespace GooFit;
using Catch::Matchers::Contains;

TEST_CASE("Adding values to unbinned dataset", "[simple][dataset]") {
    // Independent variable.
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    xvar.setValue(1);
    yvar.setValue(2);
    data.addEvent();

    xvar.setValue(3);
    yvar.setValue(4);
    data.addEvent();

    xvar.setValue(5);
    yvar.setValue(6);
    data.addEvent();

    data.loadEvent(0);

    CHECK(xvar.getValue() == 1);
    CHECK(yvar.getValue() == 2);

    data.loadEvent(1);
    CHECK(xvar.getValue() == 3_a);
    CHECK(yvar.getValue() == 4_a);

    CHECK(data.getValue(xvar, 0) == 1_a);
    CHECK(data.getValue(yvar, 0) == 2_a);
    CHECK(data.getValue(xvar, 1) == 3_a);
    CHECK(data.getValue(yvar, 1) == 4_a);
    CHECK(data.getValue(xvar, 2) == 5_a);
    CHECK(data.getValue(yvar, 2) == 6_a);
}
TEST_CASE("Variable setting and getting", "[simple][variable]") {
    // Independent variable.
    Observable var{"var", 0, 10};

    var = 1.0;

    fptype val = var;

    CHECK(val == 1.0);
    CHECK(var.getValue() == 1.0);
}

TEST_CASE("Fancy add event", "[simple][dataset]") {
    // Independent variable.
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    data.addEvent(1, 2);
    data.addEvent(3, 4);

    CHECK(data.getNumEvents() == 2);

    CHECK(data.getValue(xvar, 0) == 1_a);
    CHECK(data.getValue(yvar, 0) == 2_a);
    CHECK(data.getValue(xvar, 1) == 3_a);
    CHECK(data.getValue(yvar, 1) == 4_a);

    CHECK_THROWS_AS(data.addEvent(1), GooFit::GeneralError);
    CHECK_THROWS_AS(data.addEvent(1, 2, 3), GooFit::GeneralError);
}

TEST_CASE("Make a grid", "[simple][grid][dataset]") {
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    UnbinnedDataSet data{{xvar, yvar}};
    data.fillWithGrid();

    REQUIRE(data.getNumEvents() == 100 * 100);

    data.loadEvent(0);
    CHECK(xvar.getValue() == 0.05_a);
    CHECK(yvar.getValue() == 0.05_a);

    data.loadEvent(1);
    CHECK(xvar.getValue() == 0.15_a);
    CHECK(yvar.getValue() == 0.05_a);

    data.loadEvent(100);
    CHECK(xvar.getValue() == 0.05_a);
    CHECK(yvar.getValue() == 0.15_a);
}

TEST_CASE("Make a grid with event number", "[simple][grid][dataset]") {
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};
    EventNumber eventNumber{"eventNumber"};

    UnbinnedDataSet data{{xvar, yvar, eventNumber}};
    data.fillWithGrid();

    REQUIRE(data.getNumEvents() == 100 * 100);

    data.loadEvent(0);
    CHECK(xvar.getValue() == 0.05_a);
    CHECK(yvar.getValue() == 0.05_a);
    CHECK(eventNumber.getValue() == 0.0);

    data.loadEvent(1);
    CHECK(xvar.getValue() == 0.15_a);
    CHECK(yvar.getValue() == 0.05_a);
    CHECK(eventNumber.getValue() == 1.0);

    data.loadEvent(100);
    CHECK(xvar.getValue() == 0.05_a);
    CHECK(yvar.getValue() == 0.15_a);
    CHECK(eventNumber.getValue() == 100.0);
}

TEST_CASE("Argus variable output", "[output]") {
    Observable xvar{"xvar", 0, 1};
    Variable c{"c", xvar.getUpperLimit()};
    Variable xi{"xi", 0.25};
    Variable power{"power", .75};

    ArgusPdf argus{"argus", xvar, c, xi, true, power};

    std::stringstream str;
    str << argus;

    std::string output = "GooPdf::ArgusPdf(\"argus\") :\n"
                         "  Device function: ptr_to_Argus_Upper\n"
                         "  Observable: xvar\n"
                         "  Parameters: c, xi, power\n";

    CHECK(str.str() == output);

    CHECK_THAT(str.str(), Contains("ptr_to_Argus_Upper"));
}
