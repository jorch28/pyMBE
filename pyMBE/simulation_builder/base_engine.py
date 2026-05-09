from abc import ABC,abstractmethod
import numpy as np

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
    def calculate_center_of_mass(self, instance_id, pmb_type):
        """
        Calculates the center of mass of a pyMBE object instance in an ESPResSo system.

        Args:
            instance_id ('int'):
                pyMBE instance ID of the object whose center of mass is calculated.

            pmb_type ('str'):
                Type of the pyMBE object. Must correspond to a particle-aggregating
                template type (e.g. '"molecule"', '"residue"', '"peptide"', '"protein"').

            espresso_system ('espressomd.system.System'):
                ESPResSo system containing the particle instances.

        Returns:
            ('numpy.ndarray'):
                Array of shape '(3,)' containing the Cartesian coordinates of the
                center of mass.

        Notes:
            - This method assumes equal mass for all particles.
            - Periodic boundary conditions are *not* unfolded; positions are taken
            directly from ESPResSo particle coordinates.
        """
        center_of_mass = np.zeros(3)
        axis_list = [0,1,2]
        inst = self.db.get_instance(pmb_type=pmb_type,
                                    instance_id=instance_id)
        particle_id_list = self.db.get_particle_id_map(object_name=inst.name)["all"]
        for pid in particle_id_list:
            for axis in axis_list:
                center_of_mass[axis] += self.db.get_instance(pmb_type='particle',
                                    instance_id=pid).position[axis]
        center_of_mass = center_of_mass / len(particle_id_list)
        return center_of_mass
    def determine_reservoir_concentrations(self, pH_res, c_salt_res, activity_coefficient_monovalent_pair, max_number_sc_runs=200):
        """
        Determines ionic concentrations in the reservoir at fixed pH and salt concentration.

        Args:
            pH_res ('float'):
                Target pH value in the reservoir.

            c_salt_res ('pint.Quantity'):
                Concentration of monovalent salt (e.g., NaCl) in the reservoir.

            activity_coefficient_monovalent_pair ('callable'):
                Function returning the activity coefficient of a monovalent ion pair
                as a function of ionic strength:
                'gamma = activity_coefficient_monovalent_pair(I)'.

            max_number_sc_runs ('int', optional):
                Maximum number of self-consistent iterations allowed before
                convergence is enforced. Defaults to 200.

        Returns:
            tuple:
                (cH_res, cOH_res, cNa_res, cCl_res)
                - cH_res ('pint.Quantity'): Concentration of H⁺ ions.
                - cOH_res ('pint.Quantity'): Concentration of OH⁻ ions.
                - cNa_res ('pint.Quantity'): Concentration of Na⁺ ions.
                - cCl_res ('pint.Quantity'): Concentration of Cl⁻ ions.

        Notess:
            - The algorithm enforces electroneutrality in the reservoir.
            - Water autodissociation is included via the equilibrium constant 'Kw'.
            - Non-ideal effects enter through activity coefficients depending on
            ionic strength.
            - The implementation follows the self-consistent scheme described in
            Landsgesell (PhD thesis, Sec. 5.3, doi:10.18419/opus-10831), adapted
            from the original code (doi:10.18419/darus-2237).
        """
        def determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res):
            """
            Iteratively determines reservoir ion concentrations self-consistently.

            Args:
                cH_res ('pint.Quantity'):
                    Current estimate of the H⁺ concentration.
                c_salt_res ('pint.Quantity'):
                    Concentration of monovalent salt in the reservoir.

            Returns:
                'tuple':
                    (cH_res, cOH_res, cNa_res, cCl_res)
            """
            # Initial ideal estimate
            cOH_res = self.Kw / cH_res
            if cOH_res >= cH_res:
                cNa_res = c_salt_res + (cOH_res - cH_res)
                cCl_res = c_salt_res
            else:
                cCl_res = c_salt_res + (cH_res - cOH_res)
                cNa_res = c_salt_res
            # Self-consistent iteration
            for _ in range(max_number_sc_runs):
                ionic_strength_res = 0.5 * (cNa_res + cCl_res + cOH_res + cH_res)
                cOH_new = self.Kw / (cH_res * activity_coefficient_monovalent_pair(ionic_strength_res))
                if cOH_new >= cH_res:
                    cNa_new = c_salt_res + (cOH_new - cH_res)
                    cCl_new = c_salt_res
                else:
                    cCl_new = c_salt_res + (cH_res - cOH_new)
                    cNa_new = c_salt_res
                # Update values
                cOH_res = cOH_new
                cNa_res = cNa_new
                cCl_res = cCl_new
            return cH_res, cOH_res, cNa_res, cCl_res
        # Initial guess for H+ concentration from target pH
        cH_res = 10 ** (-pH_res) * self.units.mol / self.units.l
        # First self-consistent solve
        cH_res, cOH_res, cNa_res, cCl_res = (determine_reservoir_concentrations_selfconsistently(cH_res, 
                                                                                                 c_salt_res))
        ionic_strength_res = 0.5 * (cNa_res + cCl_res + cOH_res + cH_res)
        determined_pH = -np.log10(cH_res.to("mol/L").magnitude* np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))
        # Outer loop to enforce target pH
        while abs(determined_pH - pH_res) > 1e-6:
            if determined_pH > pH_res:
                cH_res *= 1.005
            else:
                cH_res /= 1.003
            cH_res, cOH_res, cNa_res, cCl_res = (determine_reservoir_concentrations_selfconsistently(cH_res, 
                                                                                                     c_salt_res))
            ionic_strength_res = 0.5 * (cNa_res + cCl_res + cOH_res + cH_res)
            determined_pH = -np.log10(cH_res.to("mol/L").magnitude * np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))
        return cH_res, cOH_res, cNa_res, cCl_res