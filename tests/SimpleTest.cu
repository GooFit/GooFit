#include <gtest/gtest.h>

#include "goofit/Variable.h"
#include "goofit/UnbinnedDataSet.h"

TEST(Simple, UnbinnedAdding) {

    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data {{&xvar, &yvar}};

    xvar.value = 1;
    yvar.value = 2;
    data.addEvent();
    
    xvar.value = 3;
    yvar.value = 4;
    data.addEvent();
    
    xvar.value = 5;
    yvar.value = 6;
    data.addEvent();
    
    data.loadEvent(0);
    
    EXPECT_FLOAT_EQ(1, xvar.value);
    EXPECT_FLOAT_EQ(2, yvar.value);
    
    data.loadEvent(1);
    EXPECT_FLOAT_EQ(3, xvar.value);
    EXPECT_FLOAT_EQ(4, yvar.value);
    
    EXPECT_FLOAT_EQ(1, data.GetValue(&xvar, 0));
    EXPECT_FLOAT_EQ(2, data.GetValue(&yvar, 0));
    EXPECT_FLOAT_EQ(3, data.GetValue(&xvar, 1));
    EXPECT_FLOAT_EQ(4, data.GetValue(&yvar, 1));
    EXPECT_FLOAT_EQ(5, data.GetValue(&xvar, 2));
    EXPECT_FLOAT_EQ(6, data.GetValue(&yvar, 2));
}
