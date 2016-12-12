# Methods to deal with solids tfm tab, slip off from solids_handler.py
from __future__ import print_function, absolute_import, unicode_literals, division

from mfixgui.tools.general import get_combobox_item, set_item_enabled

from mfixgui.constants import *

class SolidsTFM(object):
    def init_solids_tfm(self):
        ui = self.ui.solids
        key = 'kt_type'
        cb = ui.combobox_kt_type
        cb.activated.connect(self.set_kt_type)
        self.add_tooltip(cb, key)
        for (i, v) in enumerate(KT_TYPES):
            self.add_tooltip(get_combobox_item(cb, i), key, value=v)

        key = 'friction_model'
        cb = ui.combobox_friction_model
        self.add_tooltip(cb, key)
        for (i, v) in enumerate(FRICTION_MODELS):
            self.add_tooltip(get_combobox_item(cb, i), key, value=v)
        ui.combobox_friction_model.activated.connect(self.set_friction_model)

        ui.combobox_rdf_type.activated.connect(self.set_rdf_type)
        ui.combobox_blending_function.activated.connect(self.set_blending_function)
        ui.combobox_max_packing_correlation.activated.connect(self.set_max_packing_correlation)

    def setup_tfm_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        ui = self.ui.solids

        # Select Viscous Stress Model (KTGS):
        # Selection is unavailable for constant solids viscosity (MU_S0 defined)
        # FIXME This is not right, solids viscosity model is phase-dependent
        #enabled = (self.solids_viscosity_model != CONSTANT) # SRS p18
        #for item in (ui.label_kt_type, ui.combobox_kt_type,
        #             ui.label_friction_model, ui.combobox_friction_model):
        #    item.setEnabled(enabled)

        # SRS p18 - enable/disable menu items in viscous stress model
        kt_type = self.project.get_value('kt_type', default=DEFAULT_KT_TYPE)
        cb = ui.combobox_kt_type
        if kt_type:
            if kt_type not in KT_TYPES:
                self.warn("Invalid kt_type %s" % kt_type)
                self.unset_keyword('kt_type')
            else:
                cb.setCurrentIndex(KT_TYPES.index(kt_type))
        else:
            pass # TODO:  can kt_type be unset?
        mmax = self.project.get_value('mmax', default=1)
        k_e = (self.project.get_value('turbulence_model') == 'K_EPSILON')
        added_mass = self.project.get_value('m_am') or self.project.get_value('added_mass')
        drag_type = self.project.get_value('drag_type')
        friction_model = self.project.get_value('friction_model', default='NONE')
        enabled = [True, #ALGEBRAIC
                   True, #LUN_1984
                   True, #IA_NONEP
                   k_e,  #SIMONIN
                   k_e,  #AHMADI
                   mmax==1, #GD99
                   mmax==1, #GTSH
                   mmax<=2 and (not added_mass) and drag_type in ('WEN_YU', 'HYS')] # GHD
        #assert len(enabled) == len(KT_TYPES)
        for (i,e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)

        # Select Frictional Stress Model
        cb = ui.combobox_friction_model
        # Srivastava and Sundaresan, 2003
        # Unavailable for Algebraic Formulation viscous stress model
        set_item_enabled(get_combobox_item(cb,1), (kt_type!='ALGEBRAIC'))
        if kt_type == 'ALGEBRAIC':
            if friction_model == 'SRIVASTAVA': # Forbidden
                cb.setCurrentIndex(2) # None
                friction_model = FRICTION_MODELS[2]
                self.update_keyword('friction_model', friction_model)

        if friction_model not in FRICTION_MODELS:
            self.warn("Unrecogized friction_model %s" % friction_model)
            cb.setCurrentIndex(2) # None
            friction_model = FRICTION_MODELS[2]
            self.update_keyword('friction_model', friction_model)
        else:
            cb.setCurrentIndex(FRICTION_MODELS.index(friction_model))

        # Specify solids volume fraction at onset of friction
        enabled = (friction_model == 'SRIVASTAVA')
        for item in (ui.label_eps_f_min, ui.lineedit_keyword_eps_f_min):
            item.setEnabled(enabled)

        #Specify particle-particle restitution coefficient
        # Specification available only when required
        # Required for MMAX >=2
        # Required for viscous stress models except GHD and algebraic formulation
        # Sets keyword C_E
        enabled = (mmax>=2) or (kt_type not in ('GHD', 'ALGEBRAIC'))
        for item in (ui.label_c_e, ui.lineedit_keyword_c_e):
            item.setEnabled(enabled)

        #  Garzo, Hrenya and Dufty, 2007
        #    Selection not available for MMAX > 2
        #    Selection not available with added mass force
        #    Sets keyword KT_TYPE to GHD
        #    Requires WEN_YU or HYS drag model
        #    Specify coefficient of restitution; R_p (optional)
        # note R_P *replaces* C_E for kt_type==GHD

        ghd = (kt_type=='GHD')
        if ghd:
            names = list(self.solids.keys())

            if names:
                ui.label_r_p_1_1.setText(
                    "%s restitution coeff." % names[0])
            if len(names) > 1:
                ui.label_r_p_1_2.setText(
                    "%s-%s restitution coeff." % (names[0], names[1]))
                ui.label_r_p_2_2.setText(
                    "%s restitution coeff." % names[1])

        for item in (ui.label_c_e,
                     ui.lineedit_keyword_c_e):
            item.hide() if ghd else item.show()

        for item in (ui.label_r_p_1_1,
                     ui.lineedit_keyword_r_p_args_1_1):
            item.show() if ghd else item.hide()

        for item in (ui.label_r_p_1_2,
                     ui.label_r_p_2_2,
                     ui.lineedit_keyword_r_p_args_1_2,
                     ui.lineedit_keyword_r_p_args_2_2):
            item.show() if (ghd and len(names)>1) else item.hide()


        #Specify interphase friction coefficient
        # Specification available only when required
        # Required for MMAX >= 2
        # Sets keyword C_F
        enabled = (mmax>=2)
        for item in (ui.label_c_f, ui.lineedit_keyword_c_f):
            item.setEnabled(enabled)

        #Specify angle of particle-particle friction
        # Specification available only when required
        # Required for FRICTION_MODEL=SCHAEFFER
        # Required for FRICTION_MODEL=SRIVASTAVA
        # Sets keyword PHI
        enabled = friction_model in ('SCHAEFFER', 'SRIVASTAVA')
        for item in (ui.label_phi, ui.lineedit_keyword_phi):
            item.setEnabled(enabled)

        ### Advanced
        # Select radial distribution function
        rdf_type = self.project.get_value('rdf_type', default=DEFAULT_RDF_TYPE)
        if rdf_type not in RDF_TYPES:
            self.warn('Invalid rdf_type %s' % rdf_type)
            rdf_type = 'LEBOWITZ'
            self.update_keyword('rdf_type', rdf_type)

        if rdf_type == 'LEBOWITZ':
            if mmax == 1:
                index = 0
            else:
                index = 1
        else:
            index = RDF_TYPES.index(rdf_type)

        cb = ui.combobox_rdf_type
        cb.setCurrentIndex(index)

        enabled = [mmax==1] + 4*[mmax>1]

        for (i,e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)


        # Select stress blending model
        # Selection only available with FRICTION_MODEL=SCHAEFFER
        blending_function = self.project.get_value('blending_function',
                                                   default=DEFAULT_BLENDING_FUNCTION)
        if blending_function not in BLENDING_FUNCTIONS:
            self.warn('Invalid blending_function %s' % blending_function)
            blending_function = DEFAULT_BLENDING_FUNCTION
            self.update_keyword('blending_function', blending_function)

        ui.combobox_blending_function.setCurrentIndex(BLENDING_FUNCTIONS.index(blending_function))
        enabled = (friction_model=='SCHAEFFER')
        for item in (ui.label_blending_function, ui.combobox_blending_function):
                    item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('blending_function') #
        else:
            v = ui.combobox_blending_function.currentIndex()
            self.update_keyword('blending_function', BLENDING_FUNCTIONS[v])


        # Specify the segregation slope coefficient
        #  Only available for MMAX > 1 in conjunction with the following viscous stress
        # algebraic formulation; Lun. 1984; Simonin, 1996; Ahmadi, 1995
        enabled = (mmax>1) and kt_type in ['ALGEBRAIC', 'LUN_1984', 'SIMONIN', 'AHMADI']
        for item in (ui.label_segregation_slope_coefficient,
                     ui.lineedit_keyword_segregation_slope_coefficient):
                    item.setEnabled(enabled)

        # Select maximum packing correlation
        # Selection only available with FRICTION_MODEL=SCHAEFFER and MMAX >1
        enabled = (friction_model=='SCHAEFFER') and (mmax>1)
        for item in (ui.label_max_packing_correlation, ui.combobox_max_packing_correlation):
            item.setEnabled(enabled)
        # Constant [DEFAULT]
        # Selection always available
        # Yu & Standish
        # Selection only available for MMAX = 2
        cb = ui.combobox_max_packing_correlation
        set_item_enabled(get_combobox_item(cb, 1), mmax==2)
        yu_standis = self.project.get_value('yu_standis')
        fedors_landel = self.project.get_value('fedors_landel')
        if yu_standis and fedors_landel:
            self.warn("YU_STANDIS and FEDORS_LANDEL both set")
            self.unset_keyword('yu_standis')
            self.unset_keyword('fedors_landel')
            cb.setCurrentIndex(0) # Constant
        if yu_standis and (mmax != 2):
            self.warn("YU_STANDIS only valid for MMAX=2")
            self.unset_keyword('yu_standis')
            cb.setCurrentIndex(0) # Constant
        if yu_standis:
            cb.setCurrentIndex(1)
        elif fedors_landel:
            cb.setCurrentIndex(2)
        else:
            cb.setCurrentIndex(0)

        # Specify excluded volume in Boyle-Massoudi stress (optional)
        # Only available with algebraic formulation of viscous stress model
        enabled = (kt_type=='ALGEBRAIC')
        for item in (ui.label_v_ex, ui.lineedit_keyword_v_ex):
            item.setEnabled(enabled)


    def set_kt_type(self, val):
        kt_type = KT_TYPES[val]
        self.update_keyword('kt_type', kt_type)
        self.setup_tfm_tab()


    def set_friction_model(self, val):
        self.update_keyword('friction_model',
                            FRICTION_MODELS[val])

        self.setup_tfm_tab()

    def set_rdf_type(self, val):
        self.update_keyword('rdf_type', RDF_TYPES[val])
        self.setup_tfm_tab()

    def set_blending_function(self, val):
        self.update_keyword('blending_function', BLENDING_FUNCTIONS[val])
        self.setup_tfm_tab()

    def set_max_packing_correlation(self, val):
        self.update_keyword('yu_standis', val==1)
        self.update_keyword('fedors_landel', val==2)
        self.setup_tfm_tab()
