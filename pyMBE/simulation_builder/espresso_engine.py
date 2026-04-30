from pyMBE.simulation_builder.base_engine import SimulationEngine
import espressomd
import warnings
from  typing import List,Set



class EspressoSimulation(SimulationEngine):
    def __init__(self,box_l,db,espresso_system,units):
        self.db=db
        self.box_l: List[float]=box_l
        self.espresso_system=espresso_system
        self.units=units
        pass
    
    def _add_bond(self,particle_id1,particle_id2,bond_inst):
        bond_tpl=self.db.get_template(name=bond_inst.name, 
                                        pmb_type="bond")
        espresso_bond_inst=self._get_bond_instance(bond_template=bond_tpl)
        self.espresso_system.part.by_id(particle_id1).add_bond((espresso_bond_inst, particle_id2))
        self.db._update_instance(instance_id=bond_inst.bond_id,
                                    pmb_type='bond',
                                    attribute='added_to_engine',
                                    value=True)
    
    def _add_particle(self,particle_id):
        particle_instance=self.db.get_instance(pmb_type='particle',
                                 instance_id=particle_id)
        part_state = self.db.get_template(pmb_type="particle_state",
                                        name=particle_instance.initial_state)
        kwargs = dict(id=particle_id, 
                        pos=particle_instance.position, 
                        type=part_state.es_type, 
                        q=part_state.z,
                        fix=particle_instance.fix)        
        self.espresso_system.part.add(**kwargs)
        self.db._update_instance(instance_id=particle_id,
                                    pmb_type='particle',
                                    attribute='added_to_engine',
                                    value=True)

            
    def _check_particle_exists_in_espresso(self,particle_id):
        particle_exists=self.espresso_system.exists(particle_id)
        return particle_exists
    
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
            bond_instance = espressomd.interactions.FeneBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                            r_0 = bond_parameters["r_0"].m_as("reduced_length"),
                                                            d_r_max = bond_parameters["d_r_max"].m_as("reduced_length"))    
        return bond_instance
    
    def _delete_particles(self, particle_ids):
        """
        Remove a list of particles from an ESPResSo simulation system.

        Args:
            particle_ids  ('Iterable[int]'):
                A list (or other iterable) of ESPResSo particle IDs to remove.

            espresso_system ('espressomd.system.System'):
                The ESPResSo simulation system from which the particles
                will be removed.

        Notess:
            - This method removes particles only from the ESPResSo simulation,
            **not** from the pyMBE database. Database cleanup must be handled
            separately by the caller.
            - Attempting to remove a non-existent particle ID will raise
            an ESPResSo error.
        """
        for pid in particle_ids:
            self.espresso_system.part.by_id(pid).remove()
            self.db._update_instance(instance_id=pid,
                                     pmb_type='particle',
                                     attribute='added_to_engine',
                                     value=False)

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
    
    def _get_particle_ids_in_espresso(self):
        espresso_particles=self.espresso_system.part.all()
        return espresso_particles.id    
    
    def _get_last_particle_id_in_espresso(self):
        espresso_particles=self.espresso_system.part.all()
        if len(espresso_particles)==0:
            return None
        last_particle_id=espresso_particles.id
        return last_particle_id[-1]
    
    def _get_particle_pos_espresso(self,id):
        return self.espresso_system.part.by_id(id).pos
    
    def get_box_side_length(self):
        return self.box_l
    
    def change_volume_and_rescale_particles(self, d_new, dir="xyz"):
        """
        Change the volume for a particular dimension into the espresso system.
        args:
            d_new(float): new value for the dimension
            dir(Literal[x,y,z]): coordinate in which to set the new dimension. 
        
        """
        self.espresso_system.change_volume_and_rescale_particles(d_new=d_new,
                                                                 dir=dir)
        return
    
    def update_particle_id(self, old_pid, new_pid):
        """
        Method to be called if particles have been previously added to a simulation engine without using pyMBE.

        Args:
            old_pid (int): old particle id to be replaced.
            new_pid (int): new particle id to be assigned to instance in the pyMBE database. 
        """
        
        if self.espresso_system == None:
            raise ValueError('No simulation engine has been set to pymbe')
        
        self.db._update_instance(instance_id=old_pid,
                                 pmb_type='particle',
                                 attribute='particle_id',
                                 value=new_pid)
    
        particle_id1_instances_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                           attribute='particle_id1', 
                                                                           value=old_pid)
        for bond_id in particle_id1_instances_ids:
            self.db._update_instance(instance_id=bond_id,
                                     pmb_type='bond',
                                     attribute='particle_id1',
                                     value=new_pid)
        particle_id2_instances_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                attribute='particle_id2', 
                                                                value=old_pid)
        for bond_id in particle_id2_instances_ids:
            self.db._update_instance(instance_id=bond_id,
                                     pmb_type='bond',
                                     attribute='particle_id2',
                                     value=new_pid)

        

    def update_instances_ids_according_to_engine_particles_ids(self):
        """
        Method to be called if particles have been previously added to a simulation engine without using pymbe.
        """
        
        if self.espresso_system == None:
            raise ValueError('No simulation engine has been set to pymbe')
        
        last_id=self._get_last_particle_id_in_espresso()

        if last_id == None:
            raise ValueError('This method is intended to be used if particles have been previously added to a simulation engine')
        
        particle_instances=self.db.get_instances(pmb_type='particle')
        for i in range(particle_instances.index.size):
            self.db._update_instance(instance_id=i,
                                     pmb_type='particle',
                                     attribute='particle_id',
                                     value=last_id+i)
            
        bond_instances=self.db.get_instances(pmb_type='bond')
        for i in range(bond_instances.index.size):
            bond_instance=self.db.get_instance(pmb_type='bond',
                                 instance_id=i)
            self.db._update_instance(instance_id=i,
                                     pmb_type='bond',
                                     attribute='particle_id1',
                                     value=last_id+bond_instance.particle_id1)
            self.db._update_instance(instance_id=i,
                                     pmb_type='bond',
                                     attribute='particle_id2',
                                     value=last_id+bond_instance.particle_id2)
        warnings.warn("""Please review your setup, you have previously added a set of particles to ESPResSo. The particle ids of the pyMBE database have been updated taking into account the id of the last particle from ESPResSo. """
        )
    def add_instance_to_engine(self,pmb_type,instance_id):

        if pmb_type=='particle':
            self._add_particle(instance_id)
        elif pmb_type=='bond':

            bond_instance=self.db.get_instance(pmb_type='bond',
                                     instance_id=instance_id)
            particle_id1=bond_instance.particle_id1
            particle_id2=bond_instance.particle_id2
            self._add_bond(particle_id1,particle_id2,bond_instance)

        elif pmb_type=='residue':
            raise TypeError('Use add instances to engine to add particles and bonds not residue instances')
        elif pmb_type=='molecule':
            raise TypeError('Use add instances to engine to add particles and bonds that correspond to a molecule instance')
        elif pmb_type=='protein':
            raise TypeError('Use add instances to engine to add particles and bonds not protein instances')
        elif pmb_type=='assembly':
            raise TypeError('Use add instances to engine to add particles and bonds not assembly instances')
    
    def add_instances_to_engine(self):
        """
            This method adds the set of particles instances and bond instances 
            that are not present in the pymbe data base 
        """
        
        missing_particle_ids=self.db._find_instance_ids_by_attribute(pmb_type='particle', 
                                                                     attribute='added_to_engine', 
                                                                     value=False)
        if not missing_particle_ids:
            raise RuntimeError('No particles instances in the pyMBE database set to be added to the simulation engine')
        
        added_particle_ids=self.db._find_instance_ids_by_attribute(pmb_type='particle',
                                                                   attribute='added_to_engine', 
                                                                   value=True)
        ids_in_espresso = list(self._get_particle_ids_in_espresso())
        
        overlapping_ids = set(ids_in_espresso).intersection(set(missing_particle_ids))
        for overlapping_id in overlapping_ids:
            missing_particle_ids.remove(overlapping_id)
            new_id = max(ids_in_espresso+added_particle_ids+missing_particle_ids)+1
            self.update_particle_id(old_pid=overlapping_id, 
                                    new_pid=new_id)
            missing_particle_ids.append(new_id)

        if overlapping_ids:
            warnings.warn("""Please review your setup, you have previously added a set of particles to ESPResSo. 
                          The particle ids of the pyMBE database have been updated taking into account the id of 
                          the last particle from ESPResSo. The following particle ids were updated: {}. """.format(overlapping_ids))
        
        for id in missing_particle_ids:
            self._add_particle(id)
        
        missing_bond_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                 attribute='added_to_engine', 
                                                                 value=False)
        
        if len(missing_bond_ids)>0:
            for id in missing_bond_ids:
                bond_instance=self.db.get_instance(pmb_type='bond',
                                                   instance_id=id)
                particle_id1=bond_instance.particle_id1
                particle_id2=bond_instance.particle_id2
                self._add_bond(particle_id1,
                               particle_id2,
                               bond_instance)
       
        return 
