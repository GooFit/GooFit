#include <gtest/gtest.h>

#include "goofit/Variable.h"
#include "goofit/UnbinnedDataSet.h"

TEST(Simple, UnbinnedAdding) {

    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data {{&xvar, &yvar}};

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
    
    EXPECT_FLOAT_EQ(1, xvar.getValue());
    EXPECT_FLOAT_EQ(2, yvar.getValue());
    
    data.loadEvent(1);
    EXPECT_FLOAT_EQ(3, xvar.getValue());
    EXPECT_FLOAT_EQ(4, yvar.getValue());
    
    EXPECT_FLOAT_EQ(1, data.getValue(&xvar, 0));
    EXPECT_FLOAT_EQ(2, data.getValue(&yvar, 0));
    EXPECT_FLOAT_EQ(3, data.getValue(&xvar, 1));
    EXPECT_FLOAT_EQ(4, data.getValue(&yvar, 1));
    EXPECT_FLOAT_EQ(5, data.getValue(&xvar, 2));
    EXPECT_FLOAT_EQ(6, data.getValue(&yvar, 2));
}
