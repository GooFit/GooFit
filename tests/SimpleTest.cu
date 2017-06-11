#include "goofit/Variable.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/Error.h"

#include "catch.hpp"

using namespace GooFit;

TEST_CASE("Simple Unbinned Adding", "[simple]") {
    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}};

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

    CHECK(1 == xvar.getValue());
    CHECK(2 == yvar.getValue());

    data.loadEvent(1);
    CHECK(3 == xvar.getValue());
    CHECK(4 == yvar.getValue());

    CHECK(1 == data.getValue(&xvar, 0));
    CHECK(2 == data.getValue(&yvar, 0));
    CHECK(3 == data.getValue(&xvar, 1));
    CHECK(4 == data.getValue(&yvar, 1));
    CHECK(5 == data.getValue(&xvar, 2));
    CHECK(6 == data.getValue(&yvar, 2));
}
TEST_CASE("Simple Setting And Getting", "[simple]") {
    // Independent variable.
    Variable var{"var", 0, 10};

    var = 1.0;

    fptype val = var;

    CHECK(1.0 == val);
    CHECK(1.0 == var.getValue());
}

TEST_CASE("Fancy Add Event", "[simple]") {
    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}};

    data.addEvent(1, 2);
    data.addEvent(3, 4);

    CHECK(2 == data.getNumEvents());

    CHECK(data.getValue(&xvar, 0) == Approx(1));
    CHECK(data.getValue(&yvar, 0) == Approx(2));
    CHECK(data.getValue(&xvar, 1) == Approx(3));
    CHECK(data.getValue(&yvar, 1) == Approx(4));

    CHECK_THROWS_AS(data.addEvent(1), GeneralError);
    CHECK_THROWS_AS(data.addEvent(1, 2, 3), GeneralError);
}
