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
class EspressoSystemProtocolversion422(Protocol):
    """ Class that emulates the structure of the methods employed by Pymbe from the espressomd.System class
        . The decorator @runtime_checkable allows to only check for the structure not the types"""
    part: EspressoParticleProtocol
    bonded_inter: EspressoBondedInterProtocol


@runtime_checkable  
class EspressoSystemProtocolversion501(Protocol):
    def change_volume_and_rescale_particles(self, d_new, dir="xyz"):
        return
    def volume(self):
        return
    def distance(self, p1, p2):
        return
    def distance_vec(self, p1, p2):
        return
    def velocity_difference(self, p1, p2):
        return
    def auto_exclusions(self, distance):
        pass

@runtime_checkable
class LammpsProtocol(Protocol):
    """ Class that emulates the structure of the methods employed by Pymbe from the Lammps class
        . The decorator @runtime_checkable allows to only check for the structure not the types"""
    def display_commands(self):
        return
    def setup_lj_interactions(self):
        return
    def setup_langevin(self):
        return
    def setup_cpH(self):
        return
    def run_simulation(self):
        return

def is_engine_available(simulation_engine):
    if simulation_engine == EspressoSystemProtocolversion501:
        import contextlib
        with contextlib.suppress(ImportError):
            import espressomd
            if espressomd.version.friendly() == "4.2":
                return False
            version = espressomd.version.version()
            return version >= (5, 0, 0) and version < (5, 1, 0)
        return False
    if simulation_engine == EspressoSystemProtocolversion422:
        import contextlib
        with contextlib.suppress(ImportError):
            import espressomd
            version = espressomd.version.friendly()
            return version == "4.2"
        return False
    if simulation_engine == LammpsProtocol:
        import subprocess
        try:
            help_text = subprocess.check_output(["lmp", "-h"])
        except:
            return False
        help_text = "\n".join(help_text.decode().split("\n", 4)[:-1])
        return "Large-scale Atomic/Molecular Massively Parallel Simulator" in help_text
    raise NotImplementedError(f'Engine "{simulation_engine}" is not supported')
