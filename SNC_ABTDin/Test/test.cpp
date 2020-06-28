#include "pch.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include "../SNC_ABTDin/src/abt_compression.h"
using namespace std;

TEST(TestCase1, DeltaCode)
{
    std::vector<int> original, result;
    BitString code;
    srand((unsigned)time(NULL));
    for (int i = 0; i < 10000; ++i)
        original.push_back(i * i);
    DeltaCode delta;
    delta.Encode(original, &code);
    delta.Decode(code, &result);
    EXPECT_EQ(result, original);
}