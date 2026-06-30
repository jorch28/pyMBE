#
# Copyright (C) 2023-2026 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from typing import Protocol,runtime_checkable

class EspressoParticleProtocol(Protocol):
    """Class that emulates the estructure of the Espresso Particle class"""
    def add():
        return
    def by_id():
        return
class EspressoBondedInterProtocol(Protocol):
    """Class that emulates the structure of the EspressoBondedInterProtocol"""
    def add():
        return
    
@runtime_checkable
class EspressoSystemProtocol(Protocol):
    """ Class that emulates the structure of the methods employed by Pymbe from the espressomd.System class
        . The decorator @runtime_checkable allows to only check for the structure not the types"""
    part: EspressoParticleProtocol
    bonded_inter: EspressoBondedInterProtocol

@runtime_checkable
class LammpsProtocol(Protocol):
    """ Class that emulates the structure of the methods employed by Pymbe from the Lammps class
        . The decorator @runtime_checkable allows to only check for the structure not the types"""
    
