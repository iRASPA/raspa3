module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <tuple>
#include <optional>
#include <span>
#include <fstream>
#endif

export module mc_moves;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <fstream>;
#endif

import archive;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;
import running_energy;
import mc_moves_translation;
import mc_moves_probabilities_cross_system;

export namespace MC_Moves
{
void performRandomMove(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                       size_t selectedComponent, size_t& fractionalMoleculeSystem);

void performRandomMoveProduction(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                                 size_t selectedComponent, size_t& fractionalMoleculeSystem, size_t currentBlock);
void performRandomCrossSystemMove(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                                  MCMoveProbabilitiesCrossSystem& mc_moves_cs);
void performRandomCrossSystemMoveProduction(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                                            MCMoveProbabilitiesCrossSystem& mc_moves_cs);
};
