# Methods to deal with solids pic tab, slip off from solids_handler.py

"""Particle in Cell Model Task Pane Window: (requires PIC solver)"""

from __future__ import print_function, absolute_import, unicode_literals, division

from tools.general import get_combobox_item, set_item_enabled

class SolidsPIC(object):
    def init_solids_pic(self):
        s = self.ui.solids

    def setup_pic_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        s = self.ui.solids

        #Specify void fraction at close pack (required)
        # Sets keyword EP_STAR [0.42]
        #Define parcel CFL number
        # Sets keyword CFL_PIC
        # DEFAULT value of 0.1
        for (key, default) in [('ep_star', 0.42),
                               ('cfl_pic', 0.1)]:
            val = self.project.get_value(key)
            if val is None:
                self.update_keyword(key, default)

        #Select solid stress model
        # Selection always available
        # Available selection
        # Snider, 2001 [DEFAULT]
        # Selection always available
        # Sets keyword MPPIC_SOLID_STRESS_SNIDER=.T.
        self.update_keyword('mppic_solid_stress_snider', True)
        # TODO FIXME how does this ever get cleared?

        #Option to enable implicit treatment of drag force
        # Sets keyword: MPPIC_PDRAG_IMPLICIT
        # Disabled [DEFAULT]
        pass

        #Define particle-particle momentum retention
        # Sets keyword MPPIC_COEFF_EN1
        # DEFAULT value of 0.4
         #Define wall normal momentum retention
        # Sets keyword MPPIC_COEFF_EN_WALL
        # DEFAULT value of 0.3
         #Define wall tangential momentum retention
        # Sets keyword MPPIC_COEFF_ET_WALL
        # DEFAULT value of 0.99
        for (key, default) in [('mppic_coeff_en1', 0.4),
                               ('mppic_coeff_en_wall', 0.3),
                               ('mppic_coeff_et_wall', 0.99)]:
            val = self.project.get_value(key)
            if val is None:
                self.project.update_keyword(key, default)

        #Select gas-solids coupling scheme:
        # Selection unavailable if fluid model is disabled
        # Available selections:
        # One-way Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED true
        # Fully Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED false

        #Select interpolation framework:
        # Selection always available
        # Available selections:
        # field-to-particle and particle-to-field [DEFAULT]
        # Sets keyword DES_INTERP_ON to true
        # Sets keyword DES_INTERP_MEAN_FIELDS to true
        # field-to-particle only
        # Sets keyword DES_INTERP_ON to true
        # Sets keyword DES_INTERP_MEAN_FIELDS to false
        # particle-to-field only
        # Sets keyword DES_INTERP_ON to false
        # Sets keyword DES_INTERP_MEAN_FIELDS to true
        # no-interpolation
        # Sets keyword DES_INTERP_ON to false
        # Sets keyword DES_INTERP_MEAN_FIELDS to false

        #Select interpolation scheme:
        # Selection available except when no-interpolation framework is selected
        # Available selections:
        # None [locked default for no-interpolation framework]
        # Selection not available
        # Sets keyword DES_INTERP_SCHEME='SQUARE_DPVM'
        # Garg 2012
        # Selection not available with explicit coupling enabled
        # Sets keyword DES_INTERP_SCHEME='GARG_2012'
        # Square DPVM
        # Selection always available
        # Requires an interpolation width, DES_INTERP_WIDTH
        # Sets keyword DES_INTERP_SCHEME='SQUARE_DPVM'

        #Define interpolation width (DPVM only) (required)
        # Specification only available with SQUARE_DPVM interpolation scheme
        # Sets keyword DES_INTERP_WIDTH

        #Define solids stress model parameter: pressure constant
        # Sets keyword PSFAC_FRIC_PIC
        # DEFAULT value of 10.0
        #Define solids stress model parameter: volume fraction exponent
        # Sets keyword FRIC_EXP_PIC
        # DEFAULT value of 3.0
        #Define solids stress model parameter: non-singularity factor
        # Sets keyword FRIC_NON_SING_FAC
        # DEFAULT value of 1.0E-8
