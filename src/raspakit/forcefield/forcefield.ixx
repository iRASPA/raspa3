module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <fstream>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>
#endif

export module forcefield;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <string>;
import <algorithm>;
import <iostream>;
import <ostream>;
import <fstream>;
import <optional>;
#endif

import archive;
import double4;
import double3;
import int3;
import pseudo_atom;
import vdwparameters;
import json;

export struct ForceField
{
  enum class ChargeMethod : int
  {
    Ewald = 0,
    Coulomb = 1,
    Wolf = 2,
    ModifiedWolf = 3
  };

  enum class MixingRule : int
  {
    Lorentz_Berthelot = 0
  };

  uint64_t versionNumber{1};

  // 2D-vector, size numberOfPseudoAtoms squared
  std::vector<VDWParameters> data{};
  std::vector<bool> shiftPotentials{};
  std::vector<bool> tailCorrections{};
  MixingRule mixingRule{MixingRule::Lorentz_Berthelot};
  double cutOffVDW{12.0};
  double cutOffCoulomb{12.0};
  double dualCutOff{6.0};

  size_t numberOfPseudoAtoms{0};
  std::vector<PseudoAtom> pseudoAtoms{};

  ChargeMethod chargeMethod{ChargeMethod::Ewald};

  double overlapCriteria{1e5};

  double EwaldPrecision{1e-6};
  double EwaldAlpha{0.265058};
  int3 numberOfWaveVectors{8, 8, 8};
  bool automaticEwald{true};

  bool useCharge{true};
  bool omitEwaldFourier{false};

  double minimumRosenbluthFactor{1e-150};
  double energyOverlapCriteria{1e6};

  bool useDualCutOff{false};
  bool computePolarization{false};

  ForceField() noexcept = default;
  ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> parameters, MixingRule mixingRule,
             double cutOff, bool shifted, bool tailCorrections, bool useCharge) noexcept(false);
  ForceField(std::string filePath) noexcept(false);

  VDWParameters &operator()(size_t row, size_t col) { return data[row * numberOfPseudoAtoms + col]; }
  const VDWParameters &operator()(size_t row, size_t col) const { return data[row * numberOfPseudoAtoms + col]; }
  bool operator==(const ForceField &other) const;

  void applyMixingRule();
  void preComputePotentialShift();
  void preComputeTailCorrection();

  static std::optional<ForceField> readForceField(std::optional<std::string> directoryName,
                                                  std::string forceFieldFileName) noexcept(false);

  std::string printPseudoAtomStatus() const;
  std::string printForceFieldStatus() const;
  std::vector<nlohmann::json> jsonPseudoAtomStatus() const;
  nlohmann::json jsonForceFieldStatus() const;

  std::optional<size_t> findPseudoAtom(const std::string &name) const;
  static std::optional<size_t> findPseudoAtom(const std::vector<PseudoAtom> pseudoAtoms, const std::string &name);

  void initializeEwaldParameters(double3 perpendicularWidths);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceField &f);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceField &f);

  std::string repr() const { return printPseudoAtomStatus() + "\n" + printForceFieldStatus(); }
};
