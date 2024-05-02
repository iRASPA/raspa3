module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <sstream>
#include <string>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <map>
#include <source_location>
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module mc_moves_probabilities_cross_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <algorithm>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import archive;
import double3;
import stringutils;

std::string MCMoveProbabilitiesCrossSystem::printStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Cross system move probabilities\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Gibbs Volume Move probability:         {}\n", probabilityGibbsVolumeMove);
  std::print(stream, "Parallel Tempering Swap probability:   {}\n", probabilityParallelTemperingSwap);

  std::print(stream, "\n\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesCrossSystem &p)
{
  archive << p.versionNumber;

  archive << p.probabilityGibbsVolumeMove;
  archive << p.probabilityParallelTemperingSwap;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesCrossSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'MCMoveProbabilitiesCrossSystem' at line {} in file {}\n", location.line(),
                    location.file_name()));
  }

  archive >> p.probabilityGibbsVolumeMove;
  archive >> p.probabilityParallelTemperingSwap;

  return archive;
}
