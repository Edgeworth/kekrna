#ifndef KEKRNA_COMMON_TEST_H
#define KEKRNA_COMMON_TEST_H

#include <iostream>
#include <cstdint>
#include "energy/energy_model.h"
#include "gtest/gtest.h"

#define ONLY_FOR_THIS_MODEL(em_, hash_) \
  do { \
    auto our_hash = em_->Checksum(); \
    if (our_hash != hash_) { \
      printf("Skipping energy model specific tests: %#010x != " #hash_ " (%#010x).\n", our_hash, hash_); \
      return; \
    } \
  } while(0)

namespace kekrna {

const uint32_t T04_MODEL_HASH = 0x03b94db8;
extern energy::EnergyModelPtr g_em;
extern std::vector<energy::EnergyModelPtr> g_ems;

std::ostream& operator<<(std::ostream& os, const secondary_t& s);
std::ostream& operator<<(std::ostream& os, const computed_t& computed);

}

#endif //KEKRNA_COMMON_TEST_H
