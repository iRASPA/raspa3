module mc_moves_statistics_system;

import <string>;
import <sstream>;
import <fstream>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import archive;
import double3;
import move_statistics;
import stringutils;


void MCMoveStatisticsSystem::clear()
{
  volumeMove.clear();
  GibbsVolumeMove.clear();
}

void MCMoveStatisticsSystem::optimizeAcceptance()
{
  volumeMove.optimizeAcceptance(0.01, 1.5);
  GibbsVolumeMove.optimizeAcceptance(0.01, 1.5);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsSystem &p)
{
  archive << p.versionNumber;

  archive << p.volumeMove;
  archive << p.GibbsVolumeMove;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.volumeMove;
  archive >> p.GibbsVolumeMove;

  return archive;
}