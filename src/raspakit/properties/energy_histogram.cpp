module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <format>
#include <ranges>
#include <source_location>
#include <vector>
#endif

module property_energy_histogram;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
import <map>;
import <array>;
import <vector>;
import <algorithm>;
import <print>;
import <format>;
#endif

import double4;
import archive;
import units;


void PropertyEnergyHistogram::addSample(size_t blockIndex, double4 energy, const double &weight)
{
  size_t bin;

  bin = static_cast<size_t>((energy.x - range.first) * static_cast<double>(numberOfBins) / std::fabs(range.second - range.first));
  if(bin >=0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].x += weight;
  }
  bin = static_cast<size_t>((energy.y - range.first) * static_cast<double>(numberOfBins) / std::fabs(range.second - range.first));
  if(bin >=0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].y += weight;
  }
  bin = static_cast<size_t>((energy.z - range.first) * static_cast<double>(numberOfBins) / std::fabs(range.second - range.first));
  if(bin >=0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].z += weight;
  }
  bin = static_cast<size_t>((energy.w - range.first) * static_cast<double>(numberOfBins) / std::fabs(range.second - range.first));
  if(bin >=0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].w += weight;
  }

  numberOfCounts[blockIndex] += weight;
  totalNumberOfCounts += weight;
}

std::vector<double4> PropertyEnergyHistogram::averagedProbabilityHistogram(size_t blockIndex) const
{
  std::vector<double4> averagedData(numberOfBins);
  std::transform(bookKeepingEnergyHistogram[blockIndex].begin(),
                 bookKeepingEnergyHistogram[blockIndex].end(),
                 averagedData.begin(),
                 [&](const double4 &sample) { return sample / numberOfCounts[blockIndex]; });
  return averagedData;
}

std::vector<double4> PropertyEnergyHistogram::averagedProbabilityHistogram() const
{
  std::vector<double4> summedBlocks(numberOfBins);
  for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(),
                   bookKeepingEnergyHistogram[blockIndex].begin(),
                   summedBlocks.begin(), [](const double4 &a, const double4 &b) { return a + b; });
  }
  std::vector<double4> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double4 &sample) { return sample / totalNumberOfCounts; });

  return average;
}


std::pair<std::vector<double4>, std::vector<double4>> PropertyEnergyHistogram::averageProbabilityHistogram() const
{
  size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  std::vector<double4> average = averagedProbabilityHistogram();

  std::vector<double4> sumOfSquares(numberOfBins);
  for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<double4> blockAverage = averagedProbabilityHistogram(blockIndex);
    for (size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
    {
      double4 value = blockAverage[binIndex] - average[binIndex];
      sumOfSquares[binIndex] += value * value;
    }
  }
  std::vector<double4> standardDeviation(numberOfBins);
  std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(), [&](const double4 &sumofsquares)
                 { return sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

  std::vector<double4> standardError(numberOfBins);
  std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                 [&](const double4 &sigma) { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

  std::vector<double4> confidenceIntervalError(numberOfBins);
  std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                 [&](const double4 &error) { return intermediateStandardNormalDeviate * error; });

  return std::make_pair(average, confidenceIntervalError);
}


void PropertyEnergyHistogram::writeOutput(int systemId, size_t currentCycle)
{
   if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("energy_histogram");

  std::ofstream stream_output(std::format("energy_histogram/energy_histogram.s{}.data", systemId));

  stream_output << std::format("# energy_histogram, number of counts: {}\n", totalNumberOfCounts);
  stream_output << "# column 1: energy [K]\n";
  stream_output << "# column 2: total energy histogram [-]\n";
  stream_output << "# column 3: total energy histogram error [-]\n";
  stream_output << "# column 4: VDW energy histogram [-]\n";
  stream_output << "# column 5: VDW energy histogram error [-]\n";
  stream_output << "# column 6: Coulombic energy histogram [-]\n";
  stream_output << "# column 7: Coulombic energy histogram error [-]\n";
  stream_output << "# column 8: Polarization energy histogram [-]\n";
  stream_output << "# column 9: Polarization energy histogram error [-]\n";

  auto [average, error] = averageProbabilityHistogram();

  for (size_t bin = 0; bin != numberOfBins; ++bin)
  {
    double energy = static_cast<double>(bin) * std::fabs(range.second - range.first) / static_cast<double>(numberOfBins) + range.first;
    if(average[bin].x > 0.0)
    {
      stream_output << std::format("{} {} {} {} {} {} {} {} {}\n", 
          energy * Units::EnergyToKelvin,
          average[bin].x, error[bin].x,
          average[bin].y, error[bin].y,
          average[bin].z, error[bin].z,
          average[bin].w, error[bin].w);
    }
  }

}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &temp)
{
  archive << temp.versionNumber;
  archive << temp.numberOfBlocks;
  //archive << temp.bookKeepingTemperature;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &temp)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > temp.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyHistogram' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> temp.numberOfBlocks;
  //archive >> temp.bookKeepingTemperature;

  return archive;
}
