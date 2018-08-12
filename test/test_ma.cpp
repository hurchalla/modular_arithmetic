
#include "hurchalla/modular_arithmetic/ma.h"
//#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "gtest/gtest.h"

#include <iostream>


#define USE_GOOGLE_TESTS 1


#ifndef USE_GOOGLE_TESTS
 
int main(int argc, char *argv[])
{
   std::cout << "***Test MA***\n";
   ma();
   return 0;
}

#else
    
namespace {
    TEST(MAFunction1Test, TestName1) {
        EXPECT_TRUE(true);
    }
    TEST(MAFunction1Test, TestName2) {
        EXPECT_TRUE(true);
    }
    TEST(MAFunction2Test, TestName1) {
        EXPECT_FALSE(false);
    }
}

#endif
