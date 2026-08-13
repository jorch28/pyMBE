
import unittest as ut
from pyMBE.simulation_builder.engine_protocol import EspressoSystemProtocolversion422,EspressoSystemProtocolversion501
from pyMBE import pymbe_library
from pyMBE.simulation_builder.lammps_engine import LammpsSimulation
import espressomd

box_l=[10]*3
espresso_system=espressomd.System(box_l=box_l)
pmb=pymbe_library(seed=32)
lammps_engine=LammpsSimulation()

class Test(ut.TestCase):
    def test_espresso_engine_protocols(self):
        self.assertEqual(True,(isinstance(espresso_system,EspressoSystemProtocolversion501) or isinstance(espresso_system,EspressoSystemProtocolversion422)))
    def test_dummy_engine(self):
        with self.assertRaises(RuntimeError):
            pmb.setup_lj_interactions()
    def test_missing_simulation_engine(self):
        with self.assertRaises(ValueError):
            pmb.set_simulation_engine(simulation_engine=None,box_l=box_l)
    def test_lammps_engine_not_implemented(self):
        with self.assertRaises(NotImplementedError):
            lammps_engine.setup_lj_interactions()
if __name__=='__main__':
    ut.main()
