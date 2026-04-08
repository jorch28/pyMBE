from espresso_engine import EspressoSimulation
from lammps_engine import LammpsSimulation

class EngineFactory:
    @staticmethod
    def get_simulation_engine(engine_name,box_l,db):
        if engine_name=='espresso':
            return EspressoSimulation(Box_L=box_l,db=db)
        elif engine_name=='lammps':
            return LammpsSimulation()
        else:
            raise ValueError('Only engines "espresso" and "lammps" have been implemented')