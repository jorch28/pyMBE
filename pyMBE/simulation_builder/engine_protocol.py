from typing import Protocol,runtime_checkable

class EspressoParticleProtocol(Protocol):
    def add():
        return
    def by_id():
        return
class EspressoBondedInterProtocol(Protocol):
    def add():
        return
    
@runtime_checkable
class EspressoSystemProtocol(Protocol):
    """ Class that emulates the structure of the methods employed by Pymbe from the espressomd.System class
        . The decorator @runtime_checkable allows to only check for the structure not the types"""
    part: EspressoParticleProtocol
    bonded_inter: EspressoBondedInterProtocol
