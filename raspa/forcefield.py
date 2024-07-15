from typing import Literal

import raspalib
from .base import RaspaBase
from .utils import RASPA_DIR, SHARE_DIR
import os

class PseudoAtom(RaspaBase):
    def __init__(
        self,
        name: str = "C",
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomicNumber: int = 8,
        printToPDB: bool = False,
        source: str = "-"
    ):
        super().__init__()
        self.__name = name
        self.__mass = mass
        self.__charge = charge
        self.__polarizability = polarizability
        self.__atomicNumber = atomicNumber
        self.__printToPDB = printToPDB
        self.__source = source

        self._cpp_obj = raspalib.PseudoAtom(name, mass, charge, polarizability, atomicNumber, printToPDB, source)


class VDWParameter(RaspaBase):
    def __init__(self, epsilon: float, sigma: float):
        super().__init__()
        self.__epsilon = epsilon
        self.__sigma = sigma

        self._cpp_obj = raspalib.VDWParameters(epsilon, sigma)


class ForceField(RaspaBase):
    def __init__(
        self,
        pseudoAtoms: list[PseudoAtom],
        vdwParameters: list[VDWParameter],
        mixingRule: Literal["Lorentz-Berthelot"] = "Lorentz-Berthelot",
        cutOff: float = 12.0,
        shifted: bool = False,
        tailCorrections: bool = False,
    ):
        super().__init__()
        self.__pseudoAtoms = pseudoAtoms
        self.__vdwParameters = vdwParameters
        mixingRule = {"Lorentz-Berthelot": raspalib.ForceField.MixingRule.Lorentz_Berthelot}[mixingRule]
        self.__mixingRule = mixingRule
        self.__cutOff = cutOff
        self.__shifted = shifted
        self.__tailCorrection = tailCorrections

        self._cpp_obj = raspalib.ForceField(
            [pA._cpp_obj for pA in pseudoAtoms],
            [vP._cpp_obj for vP in vdwParameters],
            mixingRule,
            cutOff,
            shifted,
            tailCorrections,
        )

    @classmethod
    def from_json(cls, ff_path: str):
        ff = cls([PseudoAtom()], [VDWParameter(1.0, 1.0)])
        ff._cpp_obj = raspalib.readForceField(*os.path.split(ff_path))
        ff.__pseudoAtoms = ff._cpp_obj.pseudoAtoms
        ff.__vdwParameters = ff._cpp_obj.vdwParameters
        ff.__mixingRule = ff._cpp_obj.MixingRule
        return ff


exampleMoleculeForceField = ForceField.from_json(
    os.path.join(SHARE_DIR, "forcefields", "example_molecule_forcefield", "force_field.json")
)
