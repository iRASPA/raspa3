module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <span>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <chrono>
#include <type_traits>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module breakthrough_simulation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <span>;
import <cmath>;
import <string>;
import <iostream>;
import <fstream>;
import <limits>;
import <filesystem>;
import <algorithm>;
import <numeric>;
import <sstream>;
import <chrono>;
import <type_traits>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import stringutils;
import hardware_info;
import input_reader;
import component;
import system;
import simulationbox;
import mixture_prediction;
import breakthrough;


BreakthroughSimulation::BreakthroughSimulation(InputReader &inputReader):
    systems(std::move(inputReader.systems))
{
  for(System &system: systems)
  {
    std::string directoryNameString = std::format("output/system_{}/", system.systemId);
    std::filesystem::path directoryName{ directoryNameString };
    std::filesystem::create_directories(directoryName);
  }
}

void BreakthroughSimulation::run()
{
  for(System &system: systems)
  {
    Breakthrough breakthrough(system);

    std::string fileNameString = std::format("output/system_{}/output_{}_{}.data",
        system.systemId, system.temperature, system.input_pressure);
    std::ofstream fstream(fileNameString, std::ios::out );
    std::ostream stream(fstream.rdbuf());
    //std::ostream stream(std::cout.rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "{}", HardwareInfo::writeInfo());
    std::print(stream, "{}", breakthrough.writeHeader());

    breakthrough.createPlotScript();
    breakthrough.createMovieScripts();

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    breakthrough.run(stream);

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    std::chrono::duration<double> totalSimulationTime = (t2 - t1);
    std::print(stream, "Breakthrough simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count());
  }
}
