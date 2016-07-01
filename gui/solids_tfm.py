# Methods to deal with solids tfm tab, slip off from solids_handler.py
from __future__ import print_function, absolute_import, unicode_literals, division
import logging
log = logging.getLogger(__name__)

from tools.general import get_combobox_item, set_item_enabled

kt_types = ['ALGEBRAIC', 'LUN_1984', 'IA_NONEP', 'SIMONIN',
            'AHMADI', 'GD_99', 'GTSH', 'GHD']

friction_models = ['SCHAEFFER', 'SRIVASTAVA', 'NONE']

rdf_types = ['LEBOWITZ', 'LEBOWITZ', #sic
             'MANSOORI', 'MODIFIED_LEBOWITZ', 'MODIFIED_MANSOORI']

blending_functions = ['NONE', 'TANH_BLEND', 'SIGM_BLEND']

class SolidsTFM(object):
    def init_solids_tfm(self):
        s = self.ui.solids
        s.combobox_kt_type.currentIndexChanged.connect(self.set_kt_type)
        s.combobox_friction_model.currentIndexChanged.connect(self.set_friction_model)
        s.combobox_rdf_type.currentIndexChanged.connect(self.set_rdf_type)
        s.combobox_blending_function.currentIndexChanged.connect(self.set_blending_function)
        s.combobox_max_packing_correlation.currentIndexChanged.connect(self.set_max_packing_correlation)

    def setup_tfm_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        s = self.ui.solids

        # Select Viscous Stress Model (KTGS):
        # Selection is unavailable for constant solids viscosity (MU_S0 defined)
        # FIXME This is not right, solids visc. model is phase-dependent
        #enabled = (self.solids_viscosity_model != CONSTANT) # SRS p18
        #for item in (s.label_kt_type, s.combobox_kt_type,
        #             s.label_friction_model, s.combobox_friction_model):
        #    item.setEnabled(enabled)

        # SRS p18 - enable/disable menu items in viscous stress model
        kt_type = self.project.get_value('kt_type')
        cb = s.combobox_kt_type
        if kt_type:
            if kt_type not in kt_types:
                log.warn("Unrecogized kt_type %s" % kt_type)
                self.unset_keyword('kt_type')
            else:
                cb.setCurrentIndex(kt_types.index(kt_type))
        else:
            pass # TODO:  can kt_type be unset?
        mmax = self.project.get_value('mmax', default=1)
        k_e = (self.project.get_value('turbulence_model') == 'K_EPSILON')
        added_mass = self.project.get_value('m_am') or self.project.get_value('added_mass')
        drag_type = self.project.get_value('drag_type')
        friction_model = self.project.get_value('friction_model')
        enabled = [True, #Algebraic
                   True, #Lun
                   True, #Iddir
                   k_e,  #Simonin
                   k_e,  #Cao
                   mmax==1, #Garzo99
                   mmax==1, #Garzo12
                   mmax<=2 and (not added_mass) and drag_type in ('WEN_YU', 'HYS')] # Garzo07
        for (i,e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)

        # Select Frictional Stress Model
        cb = s.combobox_friction_model
        # Srivastava and Sundaresan, 2003
        # Unavailable for Algebraic Formulation viscous stress model
        set_item_enabled(get_combobox_item(cb,1), (kt_type!='ALGEBRAIC'))
        if kt_type == 'ALGEBRAIC':
            if friction_model == 'SRIVASTAVA': # Forbidden
                self.set_friction_model(2) # None

        if friction_model not in friction_models:
            log.warn("Unrecogized friction_model %s" % friction_model)
            self.set_friction_model(2) # None
        else:
            cb.setCurrentIndex(friction_models.index(friction_model))

        # Specify solids volume fraction at onset of friction
        enabled = (self.project.get_value("friction_model") == 'SRIVASTAVA')
        for item in (s.label_eps_f_min, s.lineedit_keyword_eps_f_min):
            item.setEnabled(enabled)

        #Specify particle-particle restitution coefficient
        # Specification available only when required
        # Required for MMAX >=2
        # Required for viscous stress models except GHD and algebraic formulation
        # Sets keyword C_E
        enabled = (mmax>=2) or (kt_type not in ('GHD', 'ALGEBRAIC'))
        for item in (s.label_c_e, s.lineedit_keyword_c_e):
            item.setEnabled(enabled)

        #Specify interphase friction coefficient
        # Specification available only when required
        # Required for MMAX >= 2
        # Sets keyword C_F
        enabled = (mmax>=2)
        for item in (s.label_c_f, s.lineedit_keyword_c_f):
            item.setEnabled(enabled)

        #Specify angle of particle-particle friction
        # Specification available only when required
        # Required for FRICTION_MODEL=SCHAEFFER
        # Required for FRICTION_MODEL=SRIVASTAVA
        # Sets keyword PHI
        enabled = friction_model in ('SCHAEFFER', 'SRIVASTAVA')
        for item in (s.label_phi, s.lineedit_keyword_phi):
            item.setEnabled(enabled)

        ### Advanced
        # Select radial distribution function
        rdf_type = self.project.get_value('rdf_type', 'LEBOWITZ') #default??
        if rdf_type not in rdf_types:
            log.warn('Invalid rdf_type %s' % rdf_type)
            self.update_keyword('rdf_type', 'LEBOWITZ')
            rdf_type = 'LEBOWITZ'
        if rdf_type == 'LEBOWITZ':
            if mmax == 1:
                index = 0
            else:
                index = 1
        else:
            index = rdf_types.index(rdf_type)

        cb = s.combobox_rdf_type
        cb.setCurrentIndex(index)

        enabled = [mmax==1] + 4*[mmax>1]

        for (i,e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)


        # Select stress blending model
        # Selection only available with FRICTION_MODEL=SCHAEFFER
        blending_function = self.project.get_value('blending_function', 'NONE')
        if blending_function not in blending_functions:
            log.warn('Invalid blending_function %s' % blending_function)
            self.unset_keyword('blending_function')
            blending_function = 'NONE'
        s.combobox_blending_function.setCurrentIndex(blending_functions.index(blending_function))
        enabled = (friction_model=='SCHAEFFER')
        for item in (s.label_blending_function, s.combobox_blending_function):
                    item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('blending_function') #
        else:
            v = s.combobox_blending_function.currentIndex()
            self.update_keyword('blending_function', blending_function[v])


        # Specify the segregation slope coefficient
        #  Only available for MMAX > 1 in conjunction with the following viscous stress
        # algebraic formulation; Lun. 1984; Simonin, 1996; Ahmadi, 1995
        enabled = (mmax>1) and kt_type in ['ALGEBRAIC', 'LUN_1984', 'SIMONIN', 'AHMADI']
        for item in (s.label_segregation_slope_coefficient,
                     s.lineedit_keyword_segregation_slope_coefficient):
                    item.setEnabled(enabled)

        # Select maximum packing correlation
        # Selection only available with FRICTION_MODEL=SCHAEFFER and MMAX >1
        enabled = (friction_model=='SCHAEFFER') and (mmax>1)
        for item in (s.label_max_packing_correlation, s.combobox_max_packing_correlation):
            item.setEnabled(enabled)
        # Constant [DEFAULT]
        # Selection always available
        # Yu & Standish
        # Selection only available for MMAX = 2
        cb = s.combobox_max_packing_correlation
        set_item_enabled(get_combobox_item(cb, 1), mmax==2)
        yu_standis = self.project.get_value('yu_standis')
        fedors_landel = self.project.get_value('fedors_landel')
        if yu_standis and fedors_landel:
            log.warn("YU_STANDIS and FEDORS_LANDEL both set")
            self.unset_keyword('yu_standis')
            self.unset_keyword('fedors_landel')
            cb.setCurrentIndex(0) # Constant
        if yu_standis and (mmax != 2):
            log.warn("YU_STANDIS only valid for MMAX=2")
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
        for item in (s.label_v_ex, s.lineedit_keyword_v_ex):
            item.setEnabled(enabled)


    def set_kt_type(self, val):
        s = self.ui.solids
        cb = s.combobox_kt_type
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        kt_type = kt_types[val]
        self.update_keyword('kt_type', kt_type)
        self.setup_tfm_tab()

        #set_item_enabled(get_combobox_item(self.ui.solids.combobox_friction_model,1),
        #                 enabled = (kt_type!='ALGEBRAIC')) # Algebraic model forbids Srivastava
        #friction_model = self.project.get_value('friction_model')
        #if val==0 and friction_model=='SRIVASTAVA':
        #    self.set_friction_model(2) # None

    def set_friction_model(self, val):
        s = self.ui.solids
        cb = s.combobox_friction_model
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return

        self.update_keyword('friction_model',
                            friction_models[val])

        self.setup_tfm_tab()
        # Specify solids volume fraction at onset of friction
        # Only available with FRICTION_MODEL=SRIVASTAVA
        #enabled = (val==1) # Srivastava
        #for item in (s.label_eps_f_min, s.lineedit_keyword_eps_f_min):
        #    item.setEnabled(enabled)

        #Specify angle of particle-particle friction
        # Specification available only when required
        # Required for FRICTION_MODEL=SCHAEFFER
        # Required for FRICTION_MODEL=SRIVASTAVA
        #enabled = (val in (0,1))
        #for item in (s.label_phi, s.lineedit_keyword_phi):
        #    item.setEnabled(enabled)

    def set_rdf_type(self, val):
        s = self.ui.solids
        cb = s.combobox_rdf_type
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        self.update_keyword('rdf_type', rdf_types[val])
        self.setup_tfm_tab()

    def set_blending_function(self, val):
        s = self.ui.solids
        cb = s.combobox_blending_function
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        self.update_keyword('blending_function', blending_functions[val])
        self.setup_tfm_tab()

    def set_max_packing_correlation(self, val):
        s = self.ui.solids
        cb = s.combobox_max_packing_correlation
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        self.update_keyword('yu_standis', val==1)
        self.update_keyword('fedors_landel', val==2)
        self.setup_tfm_tab()
