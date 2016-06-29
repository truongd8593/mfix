kt_types = ['ALGEBRAIC', 'LUN_1984', 'IA_NONEP', 'SIMONIN',
            'AHMADI', 'GD_99', 'GTSH', 'GHD']


class SolidsTFM(object):
    def init_solids_tfm(self):

        s = self.ui.solids
        s.combobox_kt_type.currentIndexChanged.connect(self.set_kt_type)
        s.combobox_friction_model.currentIndexChanged.connect(self.set_friction_model)

    def setup_tfm_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        s = self.ui.solids

        # FIXME This is not right, solids visc. model is phase-dependent
        #enabled = (self.solids_viscosity_model != CONSTANT) # SRS p18
        #for item in (s.label_kt_type, s.combobox_kt_type,
        #             s.label_friction_model, s.combobox_friction_model):
        #    item.setEnabled(enabled)

        enabled = (self.project.get_value("friction_model") == 'SRIVASTAVA')
        for item in (s.label_eps_f_min, s.lineedit_keyword_eps_f_min):
            item.setEnabled(enabled)

        kt_type = self.project.get_value('kt_type')
        if kt_type:
            if kt_type not in kt_types:
                log.warn("Unrecogized kt_type %s" % kt_type)
                self.unset_keyword('kt_type')
            else:
                s.combobox_kt_type.setCurrentIndex(kt_types.index(kt_type))

        mmax = self.project.get_value('mmax')

    def set_kt_type(self, val):
        self.update_keyword('kt_type', kt_types[val])

    def set_friction_model(self, val):
        print(val)
