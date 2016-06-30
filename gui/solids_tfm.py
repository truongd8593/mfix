from tools.general import get_combobox_item, set_item_enabled

kt_types = ['ALGEBRAIC', 'LUN_1984', 'IA_NONEP', 'SIMONIN',
            'AHMADI', 'GD_99', 'GTSH', 'GHD']

friction_models = ['SCHAEFFER', 'SRIVASTAVA', 'NONE']

class SolidsTFM(object):
    def init_solids_tfm(self):

        s = self.ui.solids
        s.combobox_kt_type.currentIndexChanged.connect(self.set_kt_type)
        s.combobox_friction_model.currentIndexChanged.connect(self.set_friction_model)

    def setup_tfm_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        s = self.ui.solids

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

        mmax = self.project.get_value('mmax',0)
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


    def set_kt_type(self, val):
        self.update_keyword('kt_type', kt_types[val])
        set_item_enabled(get_combobox_item(self.ui.solids.combobox_friction_model,1),
                         enabled = (val!=0)) # Algebraic model forbids Srivastava
        friction_model = self.project.get_value('friction_model')
        if val==0 and friction_model=='SRIVASTAVA':
            self.set_friction_model(2) # None

    def set_friction_model(self, val):
        s = self.ui.solids
        cb = s.combobox_friction_model
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return

        self.update_keyword('friction_model',
                            friction_models[val])
        enabled = (val==1) # Srivastava
        # Specify solids volume fraction at onset of friction

        for item in (s.label_eps_f_min, s.lineedit_keyword_eps_f_min):
            item.setEnabled(enabled)
