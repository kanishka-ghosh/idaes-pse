##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Example property package for the combustion of methane in air using
Gibbs energy minimisation.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, log, Param, \
                          PositiveReals, Reals, Set, value, Var
from pyomo.opt import SolverFactory, TerminationCondition

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBase,
                        StateBlockDataBase,
                        StateBlockBase)
from idaes.core.util.initialization import solve_indexed_blocks

# Some more inforation about this module
__author__ = "Andrew Lee, Jinliang Ma"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("PhysicalParameterBlock")
class PhysicalParameterData(PhysicalParameterBase):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self._make_params()

    def _make_params(self):
        '''
        This section is for parameters needed for the property models.
        '''
        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2', 'N2', 'O2', 'CH4',
                                              'CO', 'CO2', 'H2O', 'NH3'])

        # List of components in each phase (optional)
        self.phase_comp = {"Vap": self.component_list}

        # List of all chemical elements that constitute the chemical species
        self.element_list = Set(initialize=['H', 'N', 'O', 'C'])

        # Elemental composition of all species
        self.element_comp = {'H2': {'H': 2, 'N': 0, 'O': 0, 'C': 0},
                             'N2': {'H': 0, 'N': 2, 'O': 0, 'C': 0},
                             'O2': {'H': 0, 'N': 0, 'O': 2, 'C': 0},
                             'CH4': {'H': 4, 'N': 0, 'O': 0, 'C': 1},
                             'CO': {'H': 0, 'N': 0, 'O': 1, 'C': 1},
                             'CO2': {'H': 0, 'N': 0, 'O': 2, 'C': 1},
                             'H2O': {'H': 2, 'N': 0, 'O': 1, 'C': 0},
                             'NH3': {'H': 3, 'N': 1, 'O': 0, 'C': 0}}

        # Heat capacity parameters - Shomate eqn (from NIST webbook)
        cp_param_dict = {'H2': {1: 18.563083, 2: 12.257357, 3: -2.859786,
                                4: 0.268238, 5: 1.977990, 6: -1.147438,
                                7: 156.288133, 8: 0.0},
                         'N2': {1: 19.50583, 2: 19.88705, 3: -8.598535,
                                4: 1.369784, 5: 0.527601, 6: -4.935202,
                                7: 212.3900, 8: 0},
                         'O2': {1: 30.03235, 2: 8.772972, 3: -3.988133,
                                4: 0.788313, 5: -0.741599, 6: -11.32468,
                                7: 236.1663, 8: 0},
                         'CH4': {1: 85.81217, 2: 11.26467, 3: 2.114146,
                                 4: 0.138190, 5: -26.42221, 6: -153.5327,
                                 7: 224.4143, 8: -74.87310},
                         'CO': {1: 35.15070, 2: 1.300095, 3: -0.205921,
                                4: 0.013550, 5: -3.282780, 6: -127.8375,
                                7: 231.7120, 8: -110.5271},
                         'CO2': {1: 58.16639, 2: 2.720074, 3: -0.492289,
                                 4: 0.038844, 5: -6.447293, 6: -425.9186,
                                 7: 263.6125, 8: -393.5224},
                         'H2O': {1: 41.96426, 2: 8.622053, 3: -1.499780,
                                 4: 0.098119, 5: -11.15764, 6: -272.1797,
                                 7: 219.7809, 8: -241.8264},
                         'NH3': {1: 52.02427, 2: 18.48801, 3: -3.765128,
                                 4: 0.248541, 5: -12.45799, 6: -85.53895,
                                 7: 223.8022, 8: -45.89806}}
        self.cp_params = Param(self.component_list,
                               range(1, 9),
                               mutable=False,
                               initialize=cp_param_dict,
                               doc="Shomate equation heat capacity parameters")

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas Constant [J/mol.K]')

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325.0,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'flow_mol_comp': {'method': None, 'units': 'mol/s'},
                'pressure': {'method': None, 'units': 'Pa'},
                'temperature': {'method': None, 'units': 'K'},
                'mole_frac': {'method': None, 'units': None},
                'cp_mol': {'method': '_cp_mol', 'units': 'J/mol.K'},
                'cp_mol_comp': {'method': '_cp_mol_comp',
                                'units': 'J/mol.K'},
                'dens_mol_phase': {'method': '_dens_mol_phase',
                                   'units': 'mol/m^3'},
                'enth_mol': {'method': '_enth_mol', 'units': 'J/mol'},
                'enth_mol_comp': {'method': '_enth_mol_comp',
                                  'units': 'J/mol'},
                'entr_mol': {'method': '_entr_mol', 'units': 'J/mol'},
                'entr_mol_comp': {'method': '_entr_mol_comp',
                                  'units': 'J/mol'},
                'flow_mol': {'method': '_flow_mol', 'units': 'mol/s'},
                'gibbs_mol': {'method': '_gibbs_mol', 'units': 'J/mol'},
                'gibbs_mol_comp': {'method': '_gibbs_mol_comp',
                                   'units': 'J/mol'}})
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class StateBlock(StateBlockBase):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, flow_mol_comp=None, temperature=None, pressure=None,
                   hold_state=False, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-8}):
        '''
        Initialisation routine for property package.

        Keyword Arguments:
            flow_mol_comp : value at which to initialize component flows
                             (default=None)
            pressure : value at which to initialize pressure (default=None)
            temperature : value at which to initialize temperature
                          (default=None)
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        '''
        # Fix state variables if not already fixed
        Fcflag = {}
        Pflag = {}
        Tflag = {}

        for k in blk.keys():
            for j in blk[k].component_list:
                if blk[k].flow_mol_comp[j].fixed is True:
                    Fcflag[k, j] = True
                else:
                    Fcflag[k, j] = False
                    if flow_mol_comp is None:
                        blk[k].flow_mol_comp[j].fix(1.0)
                    else:
                        blk[k].flow_mol_comp[j].fix(flow_mol_comp[j])

            if blk[k].pressure.fixed is True:
                Pflag[k] = True
            else:
                Pflag[k] = False
                if pressure is None:
                    blk[k].pressure.fix(101325.0)
                else:
                    blk[k].pressure.fix(pressure)

            if blk[k].temperature.fixed is True:
                Tflag[k] = True
            else:
                Tflag[k] = False
                if temperature is None:
                    blk[k].temperature.fix(1500.0)
                else:
                    blk[k].temperature.fix(temperature)

            for j in blk[k].component_list:
                blk[k].mole_frac[j] = (value(blk[k].flow_mol_comp[j]) /
                                       sum(value(blk[k].flow_mol_comp[i])
                                           for i in blk[k].component_list))

        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            for j in blk[k].component_list:

                if hasattr(blk[k], "cp_shomate_eqn"):
                    blk[k].cp_mol_comp[j] = value(
                            blk[k].cp_params[j][1] +
                            blk[k].cp_params[j][2]*blk[k].temperature +
                            blk[k].cp_params[j][3]*blk[k].temperature**2 +
                            blk[k].cp_params[j][4]*blk[k].temperature**3 +
                            blk[k].cp_params[j][5]/blk[k].temperature**2)

                if hasattr(blk[k], "enthalpy_shomate_eqn"):
                    blk[k].enth_mol_comp[j] = value(
                        blk[k].cp_params[j][1]*blk[k].temperature +
                        blk[k].cp_params[j][2]*blk[k].temperature**2/2 +
                        blk[k].cp_params[j][3]*blk[k].temperature**3/3 +
                        blk[k].cp_params[j][4]*blk[k].temperature**4/4 -
                        blk[k].cp_params[j][5]/blk[k].temperature +
                        blk[k].cp_params[j][6] -
                        blk[k].cp_params[j][8])

                if hasattr(blk[k], "entropy_shomate_eqn"):
                    blk[k].entr_mol_comp[j] = value(
                            blk[k].cp_params[j][1]*log(blk[k].temperature) +
                            blk[k].cp_params[j][2]*blk[k].temperature +
                            blk[k].cp_params[j][3]*blk[k].temperature**2/2 +
                            blk[k].cp_params[j][4]*blk[k].temperature**3/3 -
                            blk[k].cp_params[j][5]/(2*blk[k].temperature**2) +
                            blk[k].cp_params[j][7] -
                            blk[k].gas_const*log(blk[k].mole_frac[j]))

                if hasattr(blk[k], "partial_gibbs_energy_eqn"):
                    blk[k].gibbs_mol_comp[j] = value(
                            blk[k].enth_mol_comp[j] -
                            blk[k].temperature*blk[k].entr_mol_comp[j])

            if hasattr(blk[k], "ideal_gas"):
                blk[k].dens_mol_phase["Vap"] = value(
                        blk[k].pressure/(blk[k].gas_const*blk[k].temperature))

            if hasattr(blk[k], "mixture_heat_capacity_eqn"):
                blk[k].cp_mol = value(sum(
                        blk[k].cp_mol_comp[j]*blk[k].mole_frac[k]
                        for j in blk[k].component_list))

            if hasattr(blk[k], "mixture_enthalpy_eqn"):
                blk[k].enth_mol = value(sum(blk[k].mole_frac[j] *
                                            blk[k].enth_mol_comp[j]
                                            for j in blk[k].component_list))

            if hasattr(blk[k], "mixture_entropy_eqn"):
                blk[k].entr_mol = value(sum(blk[k].mole_frac[j] *
                                            blk[k].entr_mol_comp[j]
                                            for j in blk[k].component_list))

            if hasattr(blk[k], "total_flow_eqn"):
                blk[k].flow_mol = value(sum(blk[k].flow_mol_comp[j]
                                            for j in blk[k].component_list))

            if hasattr(blk[k], "mixture_gibbs_eqn"):
                blk[k].gibbs_mol = value(blk[k].enth_mol -
                                         blk[k].temperature*blk[k].entr_mol)

        results = solve_indexed_blocks(opt, blk, tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info('{} Initialisation Step 1 Complete.'
                          .format(blk.name))
            else:
                _log.warning('{} Initialisation Step 1 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        flags = {"Fcflag": Fcflag, "Pflag": Pflag,
                 "Tflag": Tflag}

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} Initialisation Complete.'.format(blk.name))

        if hold_state is True:
            return flags
        else:
            blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        # Unfix state variables
        for k in blk.keys():
            for j in blk[k].component_list:
                if flags['Fcflag'][k, j] is False:
                    blk[k].flow_mol_comp[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))


@declare_process_block_class("StateBlock",
                             block_class=StateBlockBase)
class StateBlockData(StateBlockDataBase):
    """
    An example property package for ideal gas properties with Gibbs energy
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(StateBlockData, self).build()

        self._make_params()
        self._make_vars()
        self._make_constraints()
        self._make_balance_terms()

    def _make_params(self):
        # Create package parameters
        ''' This section is for parameters needed for the property models.'''

        # List of valid phases in property package
        object.__setattr__(self,
                           "phase_list",
                           self.config.parameters.phase_list)

        # Component list - a list of component identifiers
        object.__setattr__(self,
                           "component_list",
                           self.config.parameters.component_list)

        # List of all chemical elements that constitute the chemical species
        object.__setattr__(self,
                           "element_list",
                           self.config.parameters.element_list)

        # Elemental composition of all species
        object.__setattr__(self,
                           "element_comp",
                           self.config.parameters.element_comp)

        # Heat capacity correlation parameters
        object.__setattr__(self,
                           "cp_params",
                           self.config.parameters.cp_params)

        # Gas constant
        object.__setattr__(self, "gas_const",
                           self.config.parameters.gas_const)

        # Thermodynamic reference state
        object.__setattr__(self, "pressure_ref",
                           self.config.parameters.pressure_ref)
        object.__setattr__(self, "temperature_ref",
                           self.config.parameters.temperature_ref)

    def _make_vars(self):
        # Create state variables
        ''' This section contains the state variables to be used to calucate
            the properties, and needs to be linked to the state properties in
            the unit model.'''

        self.flow_mol_comp = Var(self.component_list,
                                 initialize=1.0,
                                 doc='Component molar flowrate [mol/s]')
        self.pressure = Var(domain=Reals,
                            initialize=101325.0,
                            bounds=(1e3, 1e6),
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=1500,
                               bounds=(1500, 3000),
                               doc='State temperature [K]')
        self.mole_frac = Var(self.component_list,
                             domain=Reals,
                             initialize=0.0,
                             doc='State component mole fractions [-]')

    def _make_constraints(self):
        # Create standard constraints
        ''' This section creates the necessary constraints for calculating
            the standard property values.
            All calcuations assume ideal gas behaviour
        '''
        # Mole fractions
        def mole_frac_constraint(b, j):
            return b.flow_mol_comp[j] == (
                       b.mole_frac[j] *
                       sum(b.flow_mol_comp[k] for k in b.component_list))
        self.mole_frac_constraint = Constraint(self.component_list,
                                               rule=mole_frac_constraint)

    def _dens_mol_phase(self):
        # Molar density
        self.dens_mol_phase = Var(self.phase_list,
                                  doc="Molar density")

        def ideal_gas(b, p):
            return (b.dens_mol_phase[p]*b.gas_const*b.temperature ==
                    b.pressure)
        try:
            # Try to build constraint
            self.ideal_gas = Constraint(self.phase_list, rule=ideal_gas)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mol_phase)
            self.del_component(self.ideal_gas)
            raise

    def _cp_mol_comp(self):
        # Pure component vapour heat capacities
        self.cp_mol_comp = Var(self.component_list,
                               domain=Reals,
                               initialize=1.0,
                               doc="Pure component vapour heat capacities "
                               "[J/mol.K]")

        def pure_component_cp_mol(b, j):
            return b.cp_mol_comp[j] == (
                        b.cp_params[j][1] +
                        b.cp_params[j][2]*b.temperature +
                        b.cp_params[j][3]*b.temperature**2 +
                        b.cp_params[j][4]*b.temperature**3 +
                        b.cp_params[j][5]/b.temperature**2)
        try:
            # Try to build constraint
            self.cp_shomate_eqn = Constraint(self.component_list,
                                             rule=pure_component_cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_shomate_eqn)
            raise

    def _cp_mol(self):
        # Mixture heat capacities
        self.cp_mol = Var(domain=Reals,
                          initialize=1.0,
                          doc="Mixture heat capacity [J/mol.K]")

        def cp_mol(b, j):
            return b.cp_mol == sum(b.cp_mol_comp[j]*b.mole_frac[j]
                                   for j in b.component_list)
        try:
            # Try to build constraint
            self.mixture_heat_capacity_eqn = Constraint(
                                            rule=cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol)
            self.del_component(self.mixture_heat_capacity_eqn)
            raise

    def _enth_mol_comp(self):
        # Pure component vapour enthalpies
        self.enth_mol_comp = Var(
                self.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Pure component vapour enthalpies [J/mol]")

        def pure_comp_enthalpy(b, j):
            return b.enth_mol_comp[j] == (
                    b.cp_params[j][1]*b.temperature +
                    b.cp_params[j][2]*b.temperature**2/2 +
                    b.cp_params[j][3]*b.temperature**3/3 +
                    b.cp_params[j][4]*b.temperature**4/4 -
                    b.cp_params[j][5]/b.temperature +
                    b.cp_params[j][6] -
                    b.cp_params[j][8])
        try:
            # Try to build constraint
            self.enthalpy_shomate_eqn = Constraint(self.component_list,
                                                   rule=pure_comp_enthalpy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_comp)
            self.del_component(self.enthalpy_shomate_eqn)
            raise

    def _enth_mol(self):
        # Mixture molar enthalpy
        self.enth_mol = Var(domain=Reals,
                            initialize=0.0,
                            doc='Mixture specific enthalpy [J/mol]')
        try:
            # Try to build constraint
            self.mixture_enthalpy_eqn = Constraint(expr=(
                        self.enth_mol == sum(self.mole_frac[j] *
                                             self.enth_mol_comp[j]
                                             for j in self.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol)
            self.del_component(self.mixture_enthalpy_eqn)
            raise

    def _entr_mol_comp(self):
        # Partial component vapour entropies
        self.entr_mol_comp = Var(
                self.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Partial component vapour entropies [J/mol.K]")

        def pure_comp_entropy(b, j):
            return b.entr_mol_comp[j] == (
                    b.cp_params[j][1]*log(b.temperature) +
                    b.cp_params[j][2]*b.temperature +
                    b.cp_params[j][3]*b.temperature**2/2 +
                    b.cp_params[j][4]*b.temperature**3/3 -
                    b.cp_params[j][5]/(2*b.temperature**2) +
                    b.cp_params[j][7] -
                    b.gas_const*log(b.mole_frac[j]))
        try:
            # Try to build constraint
            self.entropy_shomate_eqn = Constraint(self.component_list,
                                                  rule=pure_comp_entropy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.entr_mol_comp)
            self.del_component(self.entropy_shomate_eqn)
            raise

    def _entr_mol(self):
        # Mixture molar entropy
        self.entr_mol = Var(domain=Reals,
                            initialize=0.0,
                            doc='Mixture specific entropy [J/mol.K]')
        try:
            # Try to build constraint
            self.mixture_entropy_eqn = Constraint(expr=(
                        self.entr_mol == sum(self.mole_frac[j] *
                                             self.entr_mol_comp[j]
                                             for j in self.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.entr_mol)
            self.del_component(self.mixture_entropy_eqn)
            raise

    def _flow_mol(self):
        # Total molar material flow
        self.flow_mol = Var(domain=Reals,
                            initialize=1.0,
                            doc="Total mixture molar flowrate [mol/s]")

        try:
            # Try to build constraint
            self.total_flow_eqn = Constraint(expr=(
                        self.flow_mol == sum(self.flow_mol_comp[j]
                                             for j in self.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.flow_mol)
            self.del_component(self.total_flow_eqn)
            raise

    def _gibbs_mol_comp(self):
        # Partial component Gibbs free energy
        self.gibbs_mol_comp = Var(
                            self.component_list,
                            domain=Reals,
                            initialize=0.0,
                            doc="Partial component Gibbs free energy [J/mol]")

        # Assume constant cp_mol for simplicity and vapour phase only
        def comp_gibbs_energy_equation(b, j):
            return b.gibbs_mol_comp[j] == (b.enth_mol_comp[j] -
                                           b.temperature*b.entr_mol_comp[j])
        try:
            # Try to build constraint
            self.partial_gibbs_energy_eqn = Constraint(
                                            self.component_list,
                                            rule=comp_gibbs_energy_equation)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.gibbs_mol_comp)
            self.del_component(self.partial_gibbs_energy_eqn)
            raise

    def _gibbs_mol(self):
        # Mixture Gibbs energy
        self.gibbs_mol = Var(domain=Reals,
                             initialize=0.0,
                             doc='Mixture Gibbs energy [J/mol]')
        try:
            # Try to build constraint
            self.mixture_gibbs_eqn = Constraint(expr=(
                        self.gibbs_mol == (self.enth_mol -
                                           self.temperature*self.entr_mol)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.gibbs_mol)
            self.del_component(self.mixture_gibbs_eqn)
            raise

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return b.flow_mol*b.enth_mol

    def get_material_density_terms(b, p, j):
        return b.dens_mol_phase[p]*b.mole_frac[j]

    def get_enthalpy_density_terms(b, p):
        return b.dens_mol_phase[p]*b.enth_mol

    def define_state_vars(b):
        return {"flow_mol_comp": b.flow_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}

    def model_check(blk):
        """
        Model checks for property block
        """
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error('{} Temperature set below lower bound.'
                       .format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error('{} Temperature set above upper bound.'
                       .format(blk.name))

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error('{} Pressure set below lower bound.'.format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error('{} Pressure set above upper bound.'.format(blk.name))
