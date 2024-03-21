import RaspaKit
from .base import RaspaBase

import numpy as np

class PseudoAtom(RaspaBase):
    def __init__(self, name: str, mass: float, charge: float, atomicNumber: int, printToPDB: bool):
        self.name = name
        self.mass = mass
        self.charge = charge
        self.atomicNumber = atomicNumber
        self.printToPDB = printToPDB

        self._cpp_obj = RaspaKit.PseudoAtom(self.name, self.mass, self.charge, self.atomicNumber, self.printToPDB)

class VDWParameters:
    def __init__(self, epsilon: float, sigma: float):
        self.epsilon = epsilon
        self.sigma = sigma

        self._cpp_obj = RaspaKit.VDWParameters(self.epsilon, self.sigma)

class ForceField:
    def __init__(self, pseudoAtoms: list[PseudoAtom], parameters: list[VDWParameters], mixingRule: str, cutOff: float, shifted: bool, tailCorrection: bool):
        self.pseudoAtoms = pseudoAtoms
        self.parameters = parameters
        self.mixingRule = getattr(RaspaKit.ForceField.MixingRule, mixingRule)
        self.cutOff = cutOff
        self.shifted = shifted
        self.tailCorrection = tailCorrection

        self._cpp_obj = RaspaKit.ForceField(self.pseudoAtoms, self.parameters, self.mixingRule, self.cutOff, self.shifted, self.tailCorrection)