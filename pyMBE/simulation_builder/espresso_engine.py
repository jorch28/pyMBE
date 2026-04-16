from pyMBE.simulation_builder.base_engine import SimulationEngine
import espressomd



class EspressoSimulation(SimulationEngine):
    def __init__(self,box_l,db,espresso_system,units):
        self.db=db
        self.box_l=box_l
        self.espresso_system=espresso_system
        self.units=units
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
    
    def _get_bond_instance(self, bond_template):
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
            self.espresso_system.bonded_inter.add(bond_inst)
        return bond_inst
    
    def _add_bond(self,particle_id1,particle_id2,bond_inst):
        bond_tpl=self.db.get_template(name=bond_inst.name, 
                                        pmb_type="bond")
        espresso_bond_inst=self._get_bond_instance(bond_template=bond_tpl)
        self.espresso_system.part.by_id(particle_id1).add_bond((espresso_bond_inst, particle_id2))
    
    def _add_particle(self,particle_id):
        particle_instance=self.db.get_instance(pmb_type='particle',
                                 instance_id=particle_id)
        # part_tpl = self.db.get_template(pmb_type="particle",
        #                                 name=particle_instance.name)
        part_state = self.db.get_template(pmb_type="particle_state",
                                        name=particle_instance.initial_state)
        
        if particle_instance.fix:
            kwargs = dict(id=particle_id, pos=particle_instance.position, type=part_state.es_type, q=part_state.z,fix=particle_instance.fix)
        else:
            kwargs = dict(id=particle_id, pos=particle_instance.position, type=part_state.es_type, q=part_state.z)
        
        self.espresso_system.part.add(**kwargs)

    def get_box_side_length(self):
        return self.box_l
    
    
    def add_instances_to_engine(self):
        ### test the method
        bond_instances=self.db._get_instances_df(pmb_type='bond')
        particle_instances=self.db._get_instances_df(pmb_type='particle')
        if bond_instances.index.size>0:

            particles_added=set()
            for i in range(bond_instances.index.size):
                bond_instance=self.db.get_instance(pmb_type='bond',
                                    instance_id=i)
                particle_id1=bond_instance.particle_id1
                particle_id2=bond_instance.particle_id2

                if particle_id1 not in particles_added:
                    self._add_particle(particle_id1)
                    particles_added.add(particle_id1)

                if particle_id2 not in particles_added:
                    self._add_particle(particle_id2)
                    particles_added.add(particle_id2)

                self._add_bond(particle_id1,particle_id2,bond_instance)

        elif particle_instances.index.size>0:
            for i in range(particle_instances.index.size):
                self._add_particle(i)
        else:
            raise RuntimeError('No particles, residues or molecules have been created so far')
                    
            
            
        return 
