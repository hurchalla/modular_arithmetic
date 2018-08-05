
#include "hurchalla/modular_arithmetic/ma.h"
//#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "gtest/gtest.h"

#include <iostream>


#define USE_GOOGLE_TESTS 1


#ifndef USE_GOOGLE_TESTS

int main(int argc, char *argv[])
{
   std::cout << "***Test MA2***\n";
   ma();
   return 0;
}

#else
    
namespace {
	TEST(MA2Function1Test, TestNameA) {
		EXPECT_TRUE(true);
	}
	TEST(MA2Function1Test, TestNameB) {
		EXPECT_TRUE(true);
	}
	TEST(MA2Function2Test, TestNameA) {
		EXPECT_FALSE(false);
	}
}

#endif
