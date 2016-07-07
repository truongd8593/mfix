# Methods to deal with solids tfm tab, slip off from solids_handler.py
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtWidgets

from widgets.base import LineEdit
from tools.general import (get_combobox_item, set_item_enabled, set_item_noedit)



"""Discrete Element Model Task Pane Window: (requires DEM solver)"""

# In-line comments from MFIX-UI SRS as of 2016-07-01
#  Please update comments along with SRS/code changes!


des_intg_methods = ['EULER', 'ADAMS_BASHFORTH']
des_coll_models = ['LSD', 'HERTZIAN']
des_interp_schemes = ['SQUARE_DPVM', 'GARG_2012', 'SQUARE_DPVM']


class SolidsDEM(object):
    # After any item in this tab is edited, we call setup_dem_tab again to make
    # sure all constraints (disabled/enabled inputs) are enforced
    # Therefore, 'setup_dem_tab' must not call any of the setters

    def init_solids_dem(self):
        s = self.ui.solids
        s.checkbox_keyword_gener_part_config.clicked.connect(self.set_gener_part_config)
        s.combobox_des_intg_method.activated.connect(self.set_des_intg_method)
        s.combobox_des_coll_model.activated.connect(self.set_des_coll_model)
        s.combobox_coupling_method.activated.connect(self.set_coupling_method)
        s.checkbox_keyword_des_explicitly_coupled.clicked.connect(self.setup_dem_tab)
        s.combobox_des_interp.activated.connect(self.set_des_interp)
        s.combobox_des_interp_scheme.activated.connect(self.set_des_interp_scheme)
        s.checkbox_enable_des_diffuse_width.clicked.connect(self.enable_des_diffuse_width)


    def set_gener_part_config(self, val):
        s = self.ui.solids
        self.update_keyword('gener_part_config', val)
        if val:
            self.unset_keyword("particles")
        else:
            self.update_keyword("particles", s.lineedit_particles.text()) # Restore prev. value
        enabled = not val
        for item in (s.label_particles, s.lineedit_particles):
            item.setEnabled(enabled)

    def set_des_intg_method(self, val):
        self.update_keyword('des_intg_method', des_intg_methods[val])
        self.setup_dem_tab()

    def set_des_coll_model(self, val):
        self.update_keyword('des_coll_model', des_coll_models[val])
        self.setup_dem_tab()

    def set_coupling_method(self, val):
        self.update_keyword('des_oneway_coupled', [True,False][val])
        self.setup_dem_tab()

    def set_des_interp(self, val):
        des_interp_on = not bool(val>>1)
        des_interp_mean_fields = not bool(val%2)
        self.update_keyword('des_interp_on', des_interp_on)
        self.update_keyword('des_interp_mean_fields', des_interp_mean_fields)
        self.setup_dem_tab()

    def set_des_interp_scheme(self, val):
        des_interp_scheme = des_interp_schemes[val]
        self.update_keyword('des_interp_scheme', des_interp_scheme)
        self.setup_dem_tab()

    def enable_des_diffuse_width(self, val):
        s = self.ui.solids
        enabled = val
        for item in (s.label_des_diffuse_width, s.lineedit_keyword_des_diffuse_width,
                     s.label_des_diffuse_width_units):
            item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('des_diffuse_width')
        else: #Restore value
            self.update_keyword('des_diffuse_width',
                                s.lineedit_keyword_des_diffuse_width.text())

    def setup_dem_tab(self):
        # Ensures all constraints (items enabled/disabled) are set
        # called by each 'set_' function, so don't call those here

        s = self.ui.solids
        # Inline comments from MFIX-UI SRS as of 2016-07-01
        #  Pleas update as needed!

        #MFIX-UI SRS
        #Enable automatic particle generation
        # Enabled sets keyword GENER_PART_CONFIG to true
        # Disabled enables the user to specify number of entries in particle input file
        # Default value is 0
        # Sets keyword PARTICLES
        gener_part_config = self.project.get_value('gener_part_config')
        particles = self.project.get_value('particles', default=0)
        if particles < 0:
            self.warning("Invalid particles %s" % particles)
            self.update_keyword('particles', 0)
            particles = 0
        if particles > 0 and gener_part_config:
            self.warning("gener_part_config set, particles=%s" % particles)
            self.update_keyword('gener_part_config', False)
            gener_part_config = False
        enabled = not gener_part_config
        for item in (s.label_particles, s.lineedit_particles):
            item.setEnabled(enabled)

        #Select numerical integration method
        # Selection always available
        # Available selections
        #Euler [DEFAULT]
        # Selection always available
        # Sets keyword DES_INTG_METHOD to 'EULER'
        #Adams-Bashforth
        # Selection always available
        # Sets keyword DES_INTG_METHOD to 'ADAMS_BASHFORTH'1
        des_intg_method = self.project.get_value('des_intg_method', default='EULER')
        if des_intg_method not in des_intg_methods:
            self.warn("Invalid des_intg_method %s" % des_intg_method)
            des_intg_method = 'EULER'
        s.combobox_des_intg_method.setCurrentIndex(des_intg_methods.index(des_intg_method))

        #Selection collision model
        # Selection always available
        # Available selections
        #Linear Spring-Dashpot [DEFAULT]
        # Selection always available
        # Sets keyword DES_COLL_MODEL to 'LSD'
        #Hertzian
        # Selection always available
        # Sets keyword DES_COLL_MODEL to 'HERTZIAN'
        des_coll_model = self.project.get_value('des_coll_model', default='LSD')
        if des_coll_model not in des_coll_models:
            self.warn("Invalid des_coll_model %s" % des_coll_model)
            des_coll_model = 'LSD'
        s.combobox_des_coll_model.setCurrentIndex(des_coll_models.index(des_coll_model))

        #Select gas-solids coupling scheme:
        # Selection unavailable if fluid model is disabled
        # Available selections:
        # One-way Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED true
        # Fully Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED false
        enabled = not self.fluid_solver_disabled
        for item in (s.label_coupling_method, s.combobox_coupling_method):
            item.setEnabled(enabled)
        des_oneway_coupled = self.project.get_value('des_oneway_coupled', False) #?default?
        if des_oneway_coupled not in (True, False):
            self.warn("Invalid des_oneway_coupled %s" % des_oneway_coupled)
            des_oneway_coupled = False
            self.update_keyword('des_oneway_coupled', des_oneway_coupled)
        s.combobox_coupling_method.setCurrentIndex(0 if des_oneway_coupled else 1)

        #Optional to enable explicitly coupled simulation
        # Unavailable for GARG_2012 interpolation
        des_interp_scheme = self.project.get_value('des_interp_scheme')
        enabled = (des_interp_scheme!='GARG_2012')
        s.checkbox_keyword_des_explicitly_coupled.setEnabled(enabled)

        #Select interpolation framework:
        # Selection always available
        # Available selections:
        # Field-to-Particle and Particle-to-Field [DEFAULT]
        # Sets keyword DES_INTERP_ON to true
        # Sets keyword DES_INTERP_MEAN_FIELDS to true
        # Field-to-Particle only
        # Sets keyword DES_INTERP_ON to true
        # Sets keyword DES_INTERP_MEAN_FIELDS to false
        # Particle-to-Field only
        # Sets keyword DES_INTERP_ON to false
        # Sets keyword DES_INTERP_MEAN_FIELDS to true
        # No Interpolation
        # Sets keyword DES_INTERP_ON to false
        # Sets keyword DES_INTERP_MEAN_FIELDS to false
        des_interp_on = self.project.get_value('des_interp_on', True)
        if des_interp_on not in (True, False):
            self.warn("Invalid des_interp_on %s" % des_interp_on)
            des_interp_on = True
            self.update_keyword('des_interp_on', des_interp_on)

        des_interp_mean_fields = self.project.get_value('des_interp_mean_fields', True)
        if des_interp_mean_fields not in (True, False):
            self.warn("Invalid des_interp_mean_fields %s" % des_interp_mean_fields)
            des_interp_mean_fields = True
            self.update_keyword('des_interp_mean_fields', des_interp_mean_fields)

        index = 2*(1-des_interp_on) + (1-des_interp_mean_fields)
        s.combobox_des_interp.setCurrentIndex(index)

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
        # Sets keyword DES_INTERP_SCHEME='SQUARE_DPM'

        cb = s.combobox_des_interp_scheme
        des_interp_scheme = self.project.get_value('des_interp_scheme')
        des_explicity_coupled = self.project.get_value('des_explicity_coupled')
        interp_enabled = des_interp_on or des_interp_mean_fields # not no-interp
        cb.setEnabled(interp_enabled)
        if not interp_enabled:
            des_interp_scheme = 'SQUARE_DPVM'
            self.update_keyword('des_interp_scheme', des_interp_scheme)
            cb.setCurrentIndex(0) # Must be 'None'
        # per-item enable flags
        enabled = (not interp_enabled, not des_explicity_coupled, True)
        for (i,e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)
        # Make sure we don't leave the combobox on an invalid item!
        if not enabled[cb.currentIndex()]:
            index = enabled.index(True) # 2 is always True so this is safe
            cb.setCurrentIndex(index)
            des_interp_scheme = des_interp_schemes[index]
            self.update_keyword('des_interp_scheme', des_interp_scheme)

        #Define interpolation width (DPVM only) (required)
        # Specification only available with SQUARE_DPVM interpolation scheme
        # Sets keyword DES_INTERP_WIDTH
        enabled = interp_enabled and (des_interp_scheme=='SQUARE_DPVM') #?
        for item in (s.label_des_interp_width, s.lineedit_keyword_des_interp_width,
                     s.label_des_interp_width_units):
            item.setEnabled(enabled)

        #Option to enable diffusion of particle data
        # Selection unavailable with GARG_2012 interpolation scheme
        # No keyword is set by this option
        # Enables the user to specify a diffusion width
        # Sets keyword DES_DIFFUSE_WIDTH
        enabled = (des_interp_scheme!='GARG_2012') # and interp_enabled #?
        s.checkbox_enable_des_diffuse_width.setEnabled(enabled)
        if not enabled:
            s.checkbox_enable_des_diffuse_width.setChecked(False)
            self.unset_keyword('des_diffuse_width')
            s.lineedit_keyword_des_diffuse_width.clear() # ??? FIXME
        enabled = s.checkbox_enable_des_diffuse_width.isChecked()
        for item in (s.label_des_diffuse_width, s.lineedit_keyword_des_diffuse_width,
                     s.label_des_diffuse_width_units):
            item.setEnabled(enabled)


        #Specify friction coefficient
        # Specification always required
        # Sets keyword MEW (MEW_W)
        pass

        #Specify normal spring constant
        # Only available for LSD collision model
        # Sets keyword KN (KN_W)
        #Specify tangential spring constant factor
        # Only available for LSD collision model
        # Sets keyword KT_FAC (KT_W_FAC)
        # Default values of 2.0/7.0
        #Specify tangential damping coefficient factor
        # Only available for LSD collision model
        # Sets keyword DES_ETAT_FAC (DES_ETAT_W_FAC)
        # Default values of 0.5
        enabled = (des_coll_model=='LSD')
        for item in (s.label_kn,
                     s.lineedit_keyword_kn,
                     s.lineedit_keyword_kn_w,
                     s.label_kt_fac,
                     s.lineedit_keyword_kt_fac,
                     s.lineedit_keyword_kt_w_fac,
                     s.label_des_etat_fac,
                     s.lineedit_keyword_des_etat_fac,
                     s.lineedit_keyword_des_etat_w_fac):
            item.setEnabled(enabled)

        if enabled: # TODO set these defaults at load-time, not when this tab is shown
            for (key, default) in [('kt_fac', 2.0), ('kt_w_fac', 7.0),
                                   ('des_etat_fac', 0.5), ('des_etat_w_fac', 0.5)]:
                if self.project.get_value(key) is None:
                    self.update_keyword(key, default)
        # Unset keywords if not enabled? TODO

        #Specify Young's modulus
        # Only available for Hertzian collision model
        # Sets keyword E_YOUNG (EW_YOUNG)
        #Specify Poisson ratio:
        # Only available for Hertzian collision model
        # Sets keyword V_POISSON (VW_POISSON)
        enabled = (des_coll_model=='HERTZIAN')
        for item in (s.label_e_young,
                     s.lineedit_keyword_e_young,
                     s.lineedit_keyword_ew_young,
                     s.label_v_poisson,
                     s.lineedit_keyword_v_poisson,
                     s.lineedit_keyword_vw_poisson):
            item.setEnabled(enabled)

        #Specify normal restitution coefficient
        # Specification always required
        # Sets keyword DES_EN_INPUT (DES_EN_WALL_INPUT)
        # Input given as an upper triangular matrix
        mmax = self.project.get_value('mmax', 1)
        tw = s.tablewidget_des_en_input
        names = list(self.solids.keys())
        # Really only need to do this if a name changes
        tw.setHorizontalHeaderLabels(names)
        tw.setVerticalHeaderLabels(names + ['Wall'])
        # Table size changed
        if tw.rowCount() != mmax+1 or tw.columnCount() != mmax:
            # Clear out old lineedit widgets
            for row in range(tw.rowCount()):
                for col in range(tw.columnCount()):
                    w = tw.cellWidget(row, col)
                    if w:
                        self.project.unregister_widget(w)
                        del w
            tw.clearContents()

            # Make a new batch
            tw.setRowCount(mmax+1) # extra row for "Wall"
            tw.setColumnCount(mmax)

            def make_item(str):
                item = QtWidgets.QTableWidgetItem(str)
                set_item_noedit(item)
                set_item_enabled(item, False)
                return item
            arg = 1 # One-based
            key = 'des_en_input'
            for row in range(mmax):
                for col in range(mmax):
                    if col < row:
                        tw.setItem(row, col, make_item('--'))
                    else:
                        le = LineEdit()#FIXME.  lineedit in table works but looks a bit odd
                        le.setMaximumWidth(150)
                        le.key = key
                        le.args = [arg]
                        le.setdtype('d')
                        tw.setCellWidget(row, col, le)
                        # make sure DES_EN_INPUT is a vector
                        # (Some files have a single element)
                        val = None
                        if mmax==1:
                            val = self.project.get_value(key)
                            if val is not None:
                                self.unset_keyword(key)
                                self.update_keyword(key, val, args=[arg])
                        if val is None:
                            val = self.project.get_value(key, args=[arg])
                        if val is not None:
                            le.updateValue(key, val)
                        self.project.register_widget(le, keys=[key], args=[arg])
                        arg += 1
            key = 'des_en_wall_input'
            row = mmax
            arg = 1
            for col in range(mmax):
                le = LineEdit()
                le.setMaximumWidth(150)
                le.key = key
                le.args = [arg]
                le.setdtype('d')
                tw.setCellWidget(row, col, le)
                val = None
                if mmax==1: # vectorize single scalar
                    val = self.project.get_value(key)
                    if val is not None:
                        self.unset_keyword(key)
                        self.update_keyword(key, val, args=[arg])
                if val is None:
                    val = self.project.get_value(key, args=[arg])
                if val is not None:
                    le.updateValue(key, val)
                self.project.register_widget(le, keys=[key], args=[arg])
                arg += 1

        self.fixup_solids_table(tw, stretch_column=mmax-1)

        #Specify tangential restitution coefficient
        # Specification available for Hertzian collision model
        # Sets keyword DES_ET_INPUT (DES_ET_WALL_INPUT)
        # Input given as an upper triangular matrix
        enabled = (des_coll_model=='HERTZIAN')
        s.label_des_et_input.setEnabled(enabled)
        s.tablewidget_des_et_input.setVisible(enabled) # ?


"""
Select cohesion model
o Selection always available
o Available selections
 None [DEFAULT]
 Selection always available
 Sets keyword USE_COHESION to false
 Sets keyword VAN_DER_WAALS to false
 Van der Waals
 Selection always available
 Sets keyword USE_COHESION to true
 Sets keyword VAN_DER_WAALS to true

Specify Hamaker constant
o Specification only available for Van der Waals cohesion model
o Sets keyword HAMAKER_CONSTANT (WALL_HAMAKER_CONSTANT)

Specify outer cutoff;
o Specification only available for Van der Waals cohesion model
o Sets keyword VDW_OUTTER_CUTOFF (WALL_OUTTER_CUTOFF)

Specify inner cutoff
o Specification only available for Van der Waals cohesion model
o Sets keyword VDW_INNER_CUTOFF (WALL_INNER_CUTOFF)

Specify asperities
o Specification only available for Van der Waals cohesion model
o Sets keyword ASPERITIES

List the following options under an 'Advanced' section header.

Select Neighbor Search Method
o Selection always available
o Available selection
 Grid-based [DEFAULT]
  Selection always available
  Sets keyword DES_NEIGHBOR_SEARCH 4
 N-Square
  Selection always available
  Sets keyword DES_NEIGHBOR_SEARCH 1

Specify maximum steps between neighbor search
 Specification always available
 Sets keyword NEIGHBOR_SEARCH_N

Specify factor defining particle neighborhood
 Specification always available
 Sets keyword FACTOR_RLM

Specify neighborhood search radius ratio
 Specification always available
 Sets keyword NEIGHBOR_SEARCH_RAD_RATIO

Specify search grid partitions (optional)
 Specification always available
 Sets keyword DESGRIDSEARCH_IMAX
 Sets keyword DESGRIDSEARCH_JMAX
 Sets keyword DESGRIDSEARCH_KMAX

Enable user scalar tracking
 Selection always available
 Does not directly set any keywords
 Enables specification of number of user scalars
 Sets keyword DES_USR_VAR_SIZE

Define minimum distance for contact conduction (optional)
 Unavailable if not solving energy equations

Define fluid lens proportion constant (optional)
 Unavailable if not solving energy equations
"""
