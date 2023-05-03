module;

module property_loading;

import <string>;
import <iostream>;
import <sstream>;
import <optional>;
import <vector>;

import print;

import units;
import loadings;
import component;

std::string PropertyLoading::writeAveragesStatistics(std::vector<Component> components, std::optional<double> frameworkMass) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    const double toMolePerKg = 1000.0 / frameworkMass.value();

    std::pair<Loadings, Loadings> loadingAverage = averageLoading();

    std::print(stream, "Loadings\n");
    std::print(stream, "===============================================================================\n\n");


    for(size_t i = 0; i < components.size(); ++i)
    {
      if(components[i].type != Component::Type::Framework)
      {
        const double toMgPerG = 1000.0 * components[i].mass / frameworkMass.value();
        
        std::print(stream, "Component {} ({})\n", components[i].componentId, components[i].name);
        
        for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
        {
          Loadings blockAverage = averagedLoading(j);
          std::print(stream, "    Block[ {:2d}] {: .6e}\n", j, blockAverage.numberOfMolecules[i]);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                           loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
          toMolePerKg * loadingAverage.first.numberOfMolecules[i], toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
          toMgPerG * loadingAverage.first.numberOfMolecules[i], toMgPerG * loadingAverage.second.numberOfMolecules[i]);
        
        std::print(stream, "\n");

        for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
        {
          Loadings blockAverage = averagedLoading(j);
          std::print(stream, "    Block[ {:2d}] {: .6e}\n", j, blockAverage.numberOfMolecules[i] - components[i].amountOfExcessMolecules);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
          loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules, loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
          toMolePerKg * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules), toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
          toMgPerG * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules), toMgPerG * loadingAverage.second.numberOfMolecules[i]);

        std::print(stream, "\n\n");
      }
    }
  }
  else
  {
    const double densityConversionFactor = 1.0 / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::pair<Loadings, Loadings> loadingAverage = averageLoading();

    std::print(stream, "Densities\n");
    std::print(stream, "===============================================================================\n\n");

    for (size_t i = 0; i < components.size(); ++i)
    {
      std::print(stream, "Component {} ({})\n", components[i].componentId, components[i].name);

      for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
      {
        Loadings blockAverage = averagedLoading(j);
        std::print(stream, "    Block[ {:2d}] {}\n", j, blockAverage.numberOfMolecules[i]);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molecules]\n",
        loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molec/A^3]\n",
        loadingAverage.first.numberDensities[i], loadingAverage.second.numberDensities[i]);
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [kg/m^3]\n",
        densityConversionFactor * components[i].mass * loadingAverage.first.numberDensities[i], densityConversionFactor * components[i].mass * loadingAverage.second.numberDensities[i]);
    }
  }

  std::print(stream, "\n");

  return stream.str();
}

