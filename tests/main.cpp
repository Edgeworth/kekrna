#include "gtest/gtest.h"
#include "base.h"

using namespace kekrna;

int main(int argc, char** argv) {
  LoadEnergyModelFromDataDir("data");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
