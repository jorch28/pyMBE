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

from pyMBE.simulation_builder.base_engine import SimulationEngine


class LammpsSimulation(SimulationEngine):
    def __init__(self):
        pass
    def __getattr__(self, attr):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    def _add_angle(self):
        return
    def _check_bond_inputs(self):
        return
    def _create_bond_instance(self):
        return
    def _get_bond_instance(self):
        return
    def add_instances_to_engine(self):
        return
    
    