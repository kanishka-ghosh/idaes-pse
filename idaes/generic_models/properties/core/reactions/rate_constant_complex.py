#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Methods for calculating rate constants
"""
from pyomo.environ import exp, Var, log, value, units as pyunits

from idaes.core import MaterialFlowBasis
from idaes.generic_models.properties.core.generic.utility import \
    ConcentrationForm
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant dh_rxn
class arrhenius_complex():

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        rbasis = parent.config.reaction_basis
        if rbasis == MaterialFlowBasis.molar:
            r_base = units["amount"]
        elif rbasis == MaterialFlowBasis.mass:
            r_base = units["mass"]
        else:
            raise BurntToast(
                "{} for unexpected reaction basis {}. This should not happen "
                "so please contact the IDAES developers with this bug."
                .format(rblock.name, rbasis))

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction or
              c_form == ConcentrationForm.activity):
            r_units = r_base*units["volume"]**-1*units["time"]**-1
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += -rblock.reaction_order[p, j].value

            if c_form == ConcentrationForm.molarity:
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                c_units = units["amount"]*units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units["pressure"]
            else:
                raise BurntToast(
                    "{} received unrecognised ConcentrationForm ({}). "
                    "This should not happen - please contact the IDAES "
                    "developers with this bug."
                    .format(rblock.name, c_form))
            
            r_units = (r_base *
                       units["length"]**-3 *
                       units["time"]**-1 *
                       c_units**order)
        
        rblock.alpha_olig = Var(
            doc="Oligomerization pre-exponential")
            # units=r_units_alpha_olig)
        set_param_from_config(rblock, param="alpha_olig",config=config)
        
        rblock.alpha_crack = Var(
            doc="Cracking pre-exponential")
            # units=r_units_alpha_crack)
        set_param_from_config(rblock, param="alpha_crack",config=config)
        
        rblock.gamma = Var(
            doc="Chain-length coefficient for heat of formation",
            units=units["energy_mole"])
        set_param_from_config(rblock, param="gamma",config=config)
        
        rblock.delta = Var(
            doc="Standard heat of formation for olefins",
            units=units["energy_mole"])
        set_param_from_config(rblock, param="delta",config=config)
        
        rblock.alpha_ads = Var(
            doc="Dispersive van der Waals interaction parameter",
            units=units["energy_mole"])
        set_param_from_config(rblock, param="alpha_ads",config=config)
        
        rblock.beta_ads = Var(
            doc="Parameter for local interaction between olefin and acid site",
            units=units["energy_mole"])
        set_param_from_config(rblock, param="beta_ads",config=config)
        
        rblock.kappa_olig = Var(
            doc="Transfer coefficient for oligomerization",
            units=None)
        set_param_from_config(rblock, param="kappa_olig",config=config)
    
        rblock.kappa_crack = Var(
            doc="Transfer coefficient for cracking",
            units=None)
        set_param_from_config(rblock, param="kappa_crack",config=config)
        
        rblock.E0 = Var(
            doc="Intrinsic energy barrier",
            units=units["energy_mole"])
        set_param_from_config(rblock, param="E0",config=config)
        
        rblock.C_n = Var(
            doc="Carbon number of adsorbed reactant",
            units=None)
        set_param_from_config(rblock, param="C_n",config=config)
            
        rblock.C_m = Var(
            doc="Carbon number of gas-phase reactant reactant",
            units=None)
        set_param_from_config(rblock, param="C_m",config=config)
        
        # rblock.r_type = Var(
        #    doc="Reaction type: oligomerization (1) or cracking (2)",
        #    units=None)
        # integer boolean to indicate reaction type:
        # 1: oligomerization
        # 2: cracking
        # TODO: add validation for r_type
        rblock.r_type = config.parameter_data["r_type"]
        # set_param_from_config(rblock, param="r_type",config=config)
        
    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units
        
        lamda = log((9.9/(5.6*10**7))*1.01325*10**5)
        delH_ads = rblock.alpha_ads + rblock.C_n * rblock.beta_ads
        
        delH_formation_n = rblock.gamma * rblock.C_n + rblock.delta
        delH_formation_m = rblock.gamma * rblock.C_m + rblock.delta
        delH_formation_nm = rblock.gamma * (rblock.C_n + rblock.C_m) + rblock.delta
        
        # delH_reaction = None
        
        Ea = None
        
        log_k_nm = None
        
        
        if rblock.r_type == 1:
            delH_reaction = delH_formation_nm - delH_formation_m - delH_formation_n
            if value(delH_reaction) <= 0.0:
                Ea = rblock.E0 + rblock.kappa_olig * delH_reaction
            else:
                Ea = rblock.E0 + (1 - rblock.kappa_olig) * delH_reaction
            
            log_k_nm = rblock.alpha_olig + log(T) - log(298.15) + lamda + (delH_ads-Ea)/(
                                                                           pyunits.convert(c.gas_constant,
                                                                           to_units=units["gas_constant"])*T)
        elif rblock.r_type == 2:
            delH_reaction = delH_formation_m + delH_formation_n - delH_formation_nm
            if value(delH_reaction) <= 0.0:
                Ea = rblock.E0 + rblock.kappa_crack * delH_reaction
            else:
                Ea = rblock.E0 + (1 - rblock.kappa_crack) * delH_reaction
            
            log_k_nm = rblock.alpha_crack + log(T) - log(298.15) + lamda + (delH_ads-Ea)/(
                                                                            pyunits.convert(c.gas_constant,
                                                                            to_units=units["gas_constant"])*T)
                
        return exp(log_k_nm)
