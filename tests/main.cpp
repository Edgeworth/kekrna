#include "gtest/gtest.h"
#include "base.h"
#include "energy/energy_model.h"

using namespace kekrna;

int main(int argc, char** argv) {
  energy::LoadEnergyModelFromDataDir("data");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
