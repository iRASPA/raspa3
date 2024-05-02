module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <ios>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module parallel_tempering;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <algorithm>;
import <numeric>;
import <ranges>;
import <chrono>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <span>;
import <string>;
import <optional>;
import <fstream>;
import <sstream>;
import <filesystem>;
import <tuple>;
import <ios>;
import <complex>;
import <exception>;
import <source_location>;
import <future>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#if defined(__has_include) && __has_include(<mdspan>)
import <mdspan>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

import stringutils;
import hardware_info;
import archive;
import system;
import randomnumbers;
import mc_moves;
import input_reader;
import component;
import averages;
import loadings;
import units;
import enthalpy_of_adsorption;
import simulationbox;
import forcefield;
import sample_movies;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import atom;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import property_widom;
import property_simulationbox;
import property_energy;
import property_loading;
import property_enthalpy;
import mc_moves_probabilities_particles;
import mc_moves_cputime;
import mc_moves_count;
import property_pressure;
import transition_matrix;
import interactions_ewald;
import equation_of_states;
import mc_moves_probabilities_cross_system;

ParallelTempering::ParallelTempering() : random(std::nullopt){};

ParallelTempering::ParallelTempering(InputReader& reader) noexcept
    : numberOfCycles(reader.numberOfCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      mc_moves_probabilities_cross_system(reader.mc_moves_probabilities_cross_system),
      systems(std::move(reader.systems)),
      random(reader.randomSeed),
      estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
}

System& ParallelTempering::randomSystem()
{
  return systems[size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void ParallelTempering::run()
{
  switch (simulationStage)
  {
    case SimulationStage::Initialization:
      goto continueInitializationStage;
    case SimulationStage::Equilibration:
      goto continueEquilibrationStage;
    case SimulationStage::Production:
      goto continueProductionStage;
    default:
      break;
  }

continueInitializationStage:
  initialize();
continueEquilibrationStage:
  equilibrate();
continueProductionStage:
  production();

  output();
}

void ParallelTempering::createOutputFiles()
{
  std::filesystem::create_directories("output");
  for (System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.data", system.temperature, system.input_pressure, system.systemId);
    streams.emplace_back(fileNameString, std::ios::out);
  }
}

size_t ParallelTempering::runSystemCycleInitialize(System& system)
{
  size_t totalNumberOfMolecules{0uz};
  size_t totalNumberOfComponents{0uz};
  size_t numberOfStepsPerCycle{0uz};

  totalNumberOfMolecules = system.numberOfMolecules();
  totalNumberOfComponents = system.numerOfAdsorbateComponents();
  numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;
  for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
  {
    // move to 'slide' when implemented in llvm
    // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
    // std::ranges::views::slide(s, 2uz);

    // now using system and system for two selected components
    // system swaps and gibbs volume are not allowed anyway.

    size_t selectedComponent = system.randomComponent(random);

    MC_Moves::performRandomMove(random, system, system, selectedComponent, fractionalMoleculeSystem);

    for (Component& component : system.components)
    {
      component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
    }
  }

  return 1;
}

void ParallelTempering::initialize()
{
  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  createOutputFiles();

  for (System& system : systems)
  {
    // switch the fractional molecule on in the first system, and off in all others
    if (system.systemId == 0uz)
      system.containsTheFractionalMolecule = true;
    else
      system.containsTheFractionalMolecule = false;
  }

  for (const System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());
    std::print(stream, "{}", system.writeSystemStatus());
    std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(stream, "{}", system.forceField.printForceFieldStatus());
    std::print(stream, "{}", system.writeComponentStatus());
    std::print(stream, "{}", mc_moves_probabilities_cross_system.printStatus());
    std::print(stream, "{}", system.reactions.printStatus());
  }

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.recomputeTotalEnergies();

    std::ostream stream(streams[system.systemId].rdbuf());
    system.runningEnergies.print(stream, "Recomputed from scratch");
  };

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
    // auto task = [this](System& system) { runSystemCycleInitializate(system); };
    for (System& system : systems)
    {
      runSystemCycleInitialize(system);
      // threadPool.enqueue(task, system);
    }
    for (System& system : systems)
    {
      if (currentCycle % optimizeMCMovesEvery == 0uz)
      {
        system.optimizeMCMoves();
      }
    }

    for (size_t i = 0; i < systems.size(); ++i)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];
      MC_Moves::performRandomCrossSystemMove(random, selectedSystem, selectedSecondSystem,
                                             mc_moves_probabilities_cross_system);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::flush(stream);
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

  continueInitializationStage:;
  }
}

size_t ParallelTempering::runSystemCycleEquilibrate(System& system)
{
  size_t totalNumberOfMolecules{0uz};
  size_t totalNumberOfComponents{0uz};
  size_t numberOfStepsPerCycle{0uz};

  totalNumberOfMolecules = system.numberOfMolecules();
  totalNumberOfComponents = system.numerOfAdsorbateComponents();
  numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;
  for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
  {
    // move to 'slide' when implemented in llvm
    // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
    // std::ranges::views::slide(s, 2uz);

    // now using system and system for two selected components
    // system swaps and gibbs volume are not allowed anyway.

    size_t selectedComponent = system.randomComponent(random);

    MC_Moves::performRandomMove(random, system, system, selectedComponent, fractionalMoleculeSystem);

    system.components[selectedComponent].lambdaGC.WangLandauIteration(
        PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, system.containsTheFractionalMolecule);

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
  }

  if (currentCycle % rescaleWangLandauEvery == 0uz)
  {
    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
                                             system.containsTheFractionalMolecule);
    }
  }

  if (currentCycle % optimizeMCMovesEvery == 0uz)
  {
    system.optimizeMCMoves();
  }
  return 1;
}

void ParallelTempering::equilibrate()
{
  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.recomputeTotalEnergies();
    system.runningEnergies.print(stream, "Recomputed from scratch");

    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    // auto task = [this](System& system) { runSystemCycleEquilibration(system); };
    for (System& system : systems)
    {
      runSystemCycleEquilibrate(system);
      // threadPool2.enqueue(task, system);
    }
    for (size_t i = 0; i < systems.size(); ++i)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];
      MC_Moves::performRandomCrossSystemMove(random, selectedSystem, selectedSecondSystem,
                                             mc_moves_probabilities_cross_system);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::flush(stream);
      }
    }

    if (currentCycle % printEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }
  continueEquilibrationStage:;
  }
}

size_t ParallelTempering::runSystemCycleProduction(System& system)
{
  size_t totalNumberOfMolecules{0uz};
  size_t totalNumberOfComponents{0uz};
  size_t numberOfStepsPerCycle{0uz};

  totalNumberOfMolecules = system.numberOfMolecules();
  totalNumberOfComponents = system.numerOfAdsorbateComponents();
  numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;
  for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
  {
    // move to 'slide' when implemented in llvm
    // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
    // std::ranges::views::slide(s, 2uz);

    // now using system and system for two selected components
    // system swaps and gibbs volume are not allowed anyway.

    size_t selectedComponent = system.randomComponent(random);

    MC_Moves::performRandomMoveProduction(random, system, system, selectedComponent, fractionalMoleculeSystem,
                                          estimation.currentBin);
    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
  }

  // add the sample energy to the averages
  if (currentCycle % 10uz == 0uz || currentCycle % printEvery == 0uz)
  {
    std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
    std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
    system.currentEnergyStatus = molecularPressure.first;
    system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
    std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
    system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
    system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
  }
  if (currentCycle % optimizeMCMovesEvery == 0uz)
  {
    system.optimizeMCMoves();
  }

  system.sampleProperties(estimation.currentBin, currentCycle);

  numberOfSteps += numberOfStepsPerCycle;
  return 1;
}

void ParallelTempering::production()
{
  std::chrono::system_clock::time_point t1, t2;
  double minBias{0.0};

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.recomputeTotalEnergies();
    system.runningEnergies.print(stream, "Recomputed from scratch");

    system.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();
    system.mc_moves_count.clearCountStatistics();

    for (Component& component : system.components)
    {
      component.mc_moves_statistics.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();
      component.mc_moves_count.clearCountStatistics();

      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  minBias = std::numeric_limits<double>::max();
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      double currentMinBias =
          *std::min_element(component.lambdaGC.biasFactor.cbegin(), component.lambdaGC.biasFactor.cend());
      minBias = currentMinBias < minBias ? currentMinBias : minBias;
    }
  }
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      component.lambdaGC.normalize(minBias);
    }
  }

  ThreadPool& pool = ThreadPool::instance();

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfCycles; ++currentCycle)
  {
    std::cout << currentCycle << " ";
    t1 = std::chrono::system_clock::now();
    estimation.setCurrentSample(currentCycle);

    auto task = [this](System& system) -> size_t { return runSystemCycleProduction(system); };
    std::vector<std::future<size_t>> threads;
    for (size_t systemId = 0; systemId < systems.size(); systemId++)
    {
      threads[systemId] = pool.enqueue(task, systems[systemId]);
    }
    size_t total = 0;
    for (auto thread : threads)
    {
      total += thread.get();
    }

    for (size_t i = 0; i < systems.size(); ++i)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];
      MC_Moves::performRandomCrossSystemMoveProduction(random, selectedSystem, selectedSecondSystem,
                                                       mc_moves_probabilities_cross_system);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        std::print(stream, "{}", system.writeProductionStatusReport(currentCycle, numberOfCycles));
        std::flush(stream);
      }
    }

    // output properties to files
    for (System& system : systems)
    {
      if (system.propertyConventionalRadialDistributionFunction.has_value())
      {
        system.propertyConventionalRadialDistributionFunction->writeOutput(
            system.forceField, system.systemId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms,
            currentCycle);
      }

      if (system.propertyRadialDistributionFunction.has_value())
      {
        system.propertyRadialDistributionFunction->writeOutput(system.systemId, currentCycle);
      }
      if (system.propertyDensityGrid.has_value())
      {
        system.propertyDensityGrid->writeOutput(system.systemId, system.simulationBox, system.forceField,
                                                system.frameworkComponents, system.components, currentCycle);
      }
    }

    // write binary-restart file
    if (currentCycle % printEvery == 0uz)
    {
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }
    t2 = std::chrono::system_clock::now();
    totalSimulationTime += (t2 - t1);
  continueProductionStage:;
  }
}

void ParallelTempering::output()
{
  MCMoveCpuTime total;
  MCMoveCount countTotal;
  for (const System& system : systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_count;
  }

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies.print(stream, "Running energies");

    RunningEnergy recomputedEnergies = system.computeTotalEnergies();
    recomputedEnergies.print(stream, "Recomputed from scratch");

    RunningEnergy drift = system.runningEnergies - recomputedEnergies;
    drift.print(stream, "Monte-Carlo energy drift");

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", system.writeMCMoveStatistics());

    std::print(stream, "Production run counting of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for (const Component& component : system.components)
    {
      std::print(
          stream, "{}",
          component.mc_moves_count.writeComponentStatistics(numberOfSteps, component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_count.writeSystemStatistics(numberOfSteps));

    std::print(stream, "Production run counting of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", countTotal.writeAllSystemStatistics(numberOfSteps));

    std::print(stream, "\n\n");

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for (const Component& component : system.components)
    {
      std::print(stream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());

    std::print(stream, "Production run CPU timings of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalSimulationTime));
    std::print(stream, "\n\n");

    std::print(stream, "{}",
               system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.frameworkComponents,
                                                              system.components));
    std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());
    std::print(
        stream, "{}",
        system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swapableComponents, system.components));
    std::print(stream, "{}", system.averageLoadings.writeAveragesStatistics(system.components, system.frameworkMass));
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ParallelTempering& pt)
{
  archive << pt.versionNumber;

  archive << pt.numberOfCycles;
  archive << pt.numberOfSteps;
  archive << pt.numberOfInitializationCycles;
  archive << pt.numberOfEquilibrationCycles;

  archive << pt.printEvery;
  archive << pt.writeBinaryRestartEvery;
  archive << pt.rescaleWangLandauEvery;
  archive << pt.optimizeMCMovesEvery;

  archive << pt.currentCycle;
  archive << pt.simulationStage;

  archive << pt.systems;
  archive << pt.random;

  archive << pt.fractionalMoleculeSystem;

  archive << pt.estimation;

  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ParallelTempering& pt)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > pt.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ParallelTempering' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> pt.numberOfCycles;
  archive >> pt.numberOfSteps;
  archive >> pt.numberOfInitializationCycles;
  archive >> pt.numberOfEquilibrationCycles;

  archive >> pt.printEvery;
  archive >> pt.writeBinaryRestartEvery;
  archive >> pt.rescaleWangLandauEvery;
  archive >> pt.optimizeMCMovesEvery;

  archive >> pt.currentCycle;
  archive >> pt.simulationStage;

  archive >> pt.systems;
  archive >> pt.random;

  archive >> pt.fractionalMoleculeSystem;

  archive >> pt.estimation;

  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
  }
  std::cout << std::format("Magic number read correctly: {} vs {}\n", magicNumber, static_cast<uint64_t>(0x6f6b6179));
  return archive;
}
