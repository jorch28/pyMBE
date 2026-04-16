from abc import ABC,abstractmethod

class SimulationEngine(ABC):
    def __init__(self):
        pass
    @abstractmethod
    def _check_bond_inputs(self):
        return
    @abstractmethod
    def _create_bond_instance(self):
        return
    @abstractmethod
    def _get_bond_instance(self):
        return
    @abstractmethod
    def add_instances_to_engine(self):
        return