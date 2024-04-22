module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_probabilities_cross_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import archive;
import double3;

export struct MCMoveProbabilitiesCrossSystem
{
  uint64_t versionNumber{1};

  bool operator==(MCMoveProbabilitiesCrossSystem const &) const = default;

  double probabilityGibbsVolumeMove{0.0};
  double probabilityParallelTemperingSwap{0.0};

  void optimizeAcceptance();

  std::string printStatus() const;
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesCrossSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesCrossSystem &p);
};
