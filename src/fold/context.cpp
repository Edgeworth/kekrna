#include <stack>
#include <climits>
#include "parsing.h"
#include "fold/context.h"
#include "fold/suboptimal0.h"

namespace kekrna {
namespace fold {

using namespace constants;
using namespace energy;

constexpr context_options_t::TableAlg context_options_t::TABLE_ALGS[];
constexpr context_options_t::SuboptimalAlg context_options_t::SUBOPTIMAL_ALGS[];

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse) {
  context_options_t options;
  auto opt = argparse.GetOption("alg");
  if (opt == "0")
    options.table_alg = context_options_t::TableAlg::ZERO;
  else if (opt == "1")
    options.table_alg = context_options_t::TableAlg::ONE;
  else if (opt == "2")
    options.table_alg = context_options_t::TableAlg::TWO;
  else if (opt == "3")
    options.table_alg = context_options_t::TableAlg::THREE;
  else
    verify_expr(false, "unknown fold option");
  options.subopt_delta = atoi(argparse.GetOption("subopt-delta").c_str());
  options.subopt_num = atoi(argparse.GetOption("subopt-num").c_str());
  return options;
}

void Context::ComputeTables() {
  internal::SetGlobalState(r, *em);
  switch (options.table_alg) {
    case context_options_t::TableAlg::ZERO:
      internal::ComputeTables0();
      break;
    case context_options_t::TableAlg::ONE:
      internal::ComputeTables1();
      break;
    case context_options_t::TableAlg::TWO:
      internal::ComputeTables2();
      break;
    case context_options_t::TableAlg::THREE:
      internal::ComputeTables3();
      break;
  }
  internal::ComputeExterior();
}

computed_t Context::Fold() {
  ComputeTables();
  internal::Traceback();
  return {{internal::gr, internal::gp}, internal::gctd, internal::genergy};
}

std::vector<computed_t> Context::Suboptimal() {
  ComputeTables();
  energy_t max_energy = internal::gext[0][internal::EXT] + options.subopt_delta;
  int max_structures = options.subopt_num;
  if (options.subopt_delta < 0) max_energy = constants::CAP_E;
  if (options.subopt_num < 0) max_structures = INT_MAX;
  switch (options.suboptimal_alg) {
    case context_options_t::SuboptimalAlg::ZERO:
      return internal::Suboptimal0(max_energy, max_structures).Run();
    default:
      verify_expr(false, "bug");
  }
}

}
}