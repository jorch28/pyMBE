from base_classes import SimulationBuilder
import espressomd



class EspressoSimulation(SimulationBuilder):
    def __init__(self,box_l,db):
        self.db=db
        self.box_l=box_l
        self.espresso_system=espressomd.System(box_l = [self.box_l.to('reduced_length').magnitude]*3)
        pass

    def _check_bond_inputs(self, bond_type, bond_parameters):
        """
        Checks that the input bond parameters are valid within the current pyMBE implementation.

        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.
            
            bond_parameters ('dict'): 
                parameters of the potential of the bond.
        """
        valid_bond_types   = ["harmonic", "FENE"] 
        if bond_type not in valid_bond_types:
            raise NotImplementedError(f"Bond type '{bond_type}' currently not implemented in pyMBE, accepted types are {valid_bond_types}")
        required_parameters = {"harmonic": ["r_0","k"],
                                "FENE": ["r_0","k","d_r_max"]}
        for required_parameter in required_parameters[bond_type]:
            if required_parameter not in bond_parameters.keys():
                raise ValueError(f"Missing required parameter {required_parameter} for {bond_type} bond")
        
    def _create_bond_instance(self, bond_type, bond_parameters):
        """
        Creates an ESPResSo bond instance.

        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.

            bond_parameters ('dict'): 
                parameters of the potential of the bond.

        Notes:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain:
                - k ('Pint.Quantity')      : Magnitude of the bond. It should have units of energy/length**2 
                using the 'pmb.units' UnitRegistry.
                - r_0 ('Pint.Quantity')    : Equilibrium bond length. It should have units of length using 
                the 'pmb.units' UnitRegistry.
           
            For a FENE bond the dictionary must additionally contain:
                - d_r_max ('Pint.Quantity'): Maximal stretching length for FENE. It should have 
                units of length using the 'pmb.units' UnitRegistry. Default 'None'.

        Returns:
            ('espressomd.interactions'): instance of an ESPResSo bond object
        """
        
        self._check_bond_inputs(bond_parameters=bond_parameters,
                                bond_type=bond_type)
        if bond_type == 'harmonic':
            bond_instance = espressomd.interactions.HarmonicBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                      r_0 = bond_parameters["r_0"].m_as("reduced_length"))
        elif bond_type == 'FENE':
            bond_instance    = espressomd.interactions.FeneBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                      r_0 = bond_parameters["r_0"].m_as("reduced_length"),
                                                      d_r_max = bond_parameters["d_r_max"].m_as("reduced_length"))    
        return bond_instance
    
    def _get_bond_instance(self, bond_template, espresso_system):
        """
        Retrieve or create a bond instance in an ESPResSo system for a given pair of particle names.

        Args:
            bond_template ('BondTemplate'): 
                BondTemplate object from the pyMBE database.
            espresso_system ('espressomd.system.System'): 
                An ESPResSo system object where the bond will be added or retrieved.

        Returns:
            ('espressomd.interactions.BondedInteraction'): 
                The ESPResSo bond instance object.

        Notes:
            When a new bond instance is created, it is not added to the ESPResSo system.
        """
        if bond_template.name in self.db.espresso_bond_instances.keys():
            bond_inst = self.db.espresso_bond_instances[bond_template.name]
        else:   
            # Create an instance of the bond 
            bond_inst = self._create_bond_instance(bond_type=bond_template.bond_type,
                                                            bond_parameters=bond_template.get_parameters(self.units))
            self.db.espresso_bond_instances[bond_template.name]= bond_inst
            espresso_system.bonded_inter.add(bond_inst)
        return bond_inst
    def get_box_side_length(self):
        return self.Box_L
    def save_molecule(self):
        return 
