# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
from copy import deepcopy

from qtpy import QtWidgets, PYQT5
from qtpy.QtWidgets import QCheckBox, QComboBox

from qtpy.QtGui import QValidator

from mfixgui.tools.general import (set_item_noedit, set_item_enabled,
                           widget_iter,
                           get_selected_row, get_combobox_item)

from mfixgui.widgets.base import (LineEdit, ComboBox)

from mfixgui.reaction_parser import ReactionParser

from json import JSONDecoder, JSONEncoder

class ExtLineEdit(LineEdit):
    # Track focus of lineedits, to emulate selection behavior
    def focusInEvent(self, ev):
        self.tw.current_row = self.row
        if self.side == 'reactants':
            self.parent.chemistry_handle_reactant_selection(row=self.row)
        else:
            self.parent.chemistry_handle_product_selection(row=self.row)
        LineEdit.focusInEvent(self, ev)
        #self.selectAll()

#Column numbers
COL_ENABLE, COL_RXN_NAME, COL_CHEM_EQ = (0,1,2)
COL_PHASE, COL_SPECIES, COL_COEFF = (0,1,2)

class Chemistry(object):
    #Chemistry Task Pane Window: This section allows a user to define chemical reaction input.

    def init_chemistry(self):
        ui = self.ui.chemistry

        # data members
        self.current_reaction_name = None
        self.reaction_edited = False
        self.working_reaction = None
        self.reaction_mass_totals = [None, None]
        self.disabled_reactions = {}

        # Toolbuttons
        ui.toolbutton_add_reaction.clicked.connect(self.chemistry_add_reaction)
        ui.toolbutton_delete_reaction.clicked.connect(lambda checked: self.chemistry_delete_reaction())
        ui.toolbutton_add_reactant.clicked.connect(self.chemistry_add_reactant)
        ui.toolbutton_delete_reactant.clicked.connect(self.chemistry_delete_reactant)
        ui.toolbutton_add_product.clicked.connect(self.chemistry_add_product)
        ui.toolbutton_delete_product.clicked.connect(self.chemistry_delete_product)
        ui.toolbutton_ok.clicked.connect(self.chemistry_apply_changes)
        ui.toolbutton_cancel.clicked.connect(self.chemistry_revert_changes)

        ui.toolbutton_delete_reaction.setEnabled(False) # Need a selection
        ui.toolbutton_delete_reactant.setEnabled(False)
        ui.toolbutton_delete_product.setEnabled(False)

        ui.toolbutton_ok.setEnabled(False)
        ui.toolbutton_cancel.setEnabled(False)

        # Tablewidgets
        ui.tablewidget_reactions.itemSelectionChanged.connect(self.chemistry_handle_selection)

        # Note, since tw_reactants and tw_products are populated with cell widgets, selection
        # callbacks do not work http://www.qtcentre.org/threads/56547-QTableWidget-cellwidget-selection
        #ui.tablewidget_reactants.itemSelectionChanged.connect(self.chemistry_handle_reactant_selection)
        #ui.tablewidget_products.itemSelectionChanged.connect(self.chemistry_handle_product_selection)
        ui.tablewidget_reactants.current_row = None
        ui.tablewidget_products.current_row = None

        # Set reaction name
        ui.lineedit_reaction_name.editingFinished.connect(self.set_reaction_name)
        class RxnIdValidator(QValidator):
            #  Alphanumeric combinations (no special characters excluding underscores)
            #  Limited to 32 characters
            #  First character must be a letter
            #  No blank spaces
            def __init__(self, parent=None):
                super(RxnIdValidator, self).__init__()
                self.parent = parent

            def validate(self, text, pos):
                #self.parent.len = len(text)
                #if len(text) == 0:
                #    # How to reset the lineedit after user blanks out input?
                #    return (QValidator.Intermediate, text, pos)
                if len(text) == 0:
                    return (QValidator.Acceptable, text, pos)
                elif 1 <= len(text) <= 32 and text[0].isalpha() and all(c.isalnum() or c=='_' for c in text):
                    if text.lower() in self.parent.keyword_doc: # cannot use keywords as reaction names!
                        return (QValidator.Intermediate, text, pos)
                    else:
                        return (QValidator.Acceptable, text, pos)
                else:
                    return (QValidator.Invalid, text, pos)
        ui.lineedit_reaction_name.setValidator(RxnIdValidator(parent=self))

        ui.groupbox_dh.clicked.connect(self.handle_dh_checkbox)

        ui.lineedit_dh.key = 'dh'
        ui.lineedit_dh.min = 0
        ui.lineedit_dh.dtype = float
        for le in ui.lineedit_fracdh_1, ui.lineedit_fracdh_2:
            le.key = 'fracdh'
            le.dtype = float
            le.min = 0
            le.max = 1
            le.value_updated.connect(self.handle_fracdh)
        ui.lineedit_fracdh_1.args = [0]
        ui.lineedit_fracdh_2.args = [1]

    def handle_dh_checkbox(self, enabled):
        ui = self.ui.chemistry
        for w in (ui.label_dh, ui.lineedit_dh, ui.label_dh_units,
                  ui.label_fracdh_1, ui.lineedit_fracdh_1,
                  ui.label_fracdh_2, ui.lineedit_fracdh_2):
            w.setEnabled(enabled)
        reaction = self.working_reaction
        if reaction is None:
            return
        if enabled:
            phases = self.chemistry_reaction_phases(reaction)
            num_phases = len(phases)

            dh = reaction.get('dh')
            if dh is None:
                dh = ui.lineedit_dh.value #restore
                if dh=='':
                    dh = None
                if dh is not None:
                    reaction['dh'] = dh
            if dh is None:
                dh = 0.0 #Default
                reaction['dh'] = dh
            ui.lineedit_dh.updateValue('dh', dh)

            if num_phases == 0:
                reaction.pop('fracdh', None)
            elif num_phases == 1:
                reaction['fracdh'] = {phases[0], 1.0}
            elif num_phases == 2:
                fracdh1 = ui.lineedit_fracdh_1.value
                fracdh2 = ui.lineedit_fracdh_2.value
                if fracdh1 == '': fracdh1 = None
                if fracdh2 == '': fracdh2 = None
                if fracdh1 is None and fracdh2 is None:
                    fracdh1 = fracdh2 = 0.5
                elif fracdh1 is None:
                    fracdh1 = fracdh2 - 1.0
                else:
                    fracdh2 = 1.0 - fracdh1 # ensure they sum to 1
                ui.lineedit_fracdh_1.updateValue('fracdh', fracdh1)
                ui.lineedit_fracdh_2.updateValue('fracdh', fracdh2)
                reaction['fracdh'] = {phases[0]: fracdh1, phases[1]: fracdh2}
            else:
                self.error("num_phases > 2")
                return
        else:
            reaction.pop('dh', None)
            reaction.pop('fracdh', None)
        self.reaction_edited = True
        self.chemistry_update_buttons()


    def handle_dh(self, widget, val, args):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        reaction['dh'] = val['dh']
        self.reaction_edited = True
        self.chemistry_update_buttons()


    def handle_fracdh(self, widget, vals, args):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        phases = self.chemistry_reaction_phases(reaction)
        num_phases = len(phases)
        if num_phases != 2: # What are we doing here?
            return
        val = vals.get('fracdh')
        if val in (None, ''):
            val = 0.0
        if args==[0]:
            reaction['fracdh'] = {phases[0]: val, phases[1]: 1.0-val}
            ui.lineedit_fracdh_2.updateValue('fracdh', 1.0-val)
        elif args==[1]:
            reaction['fracdh'] = {phases[0]: 1.0-val, phases[1]: val}
            ui.lineedit_fracdh_1.updateValue('fracdh', 1.0-val)

        self.reaction_edited = True
        self.chemistry_update_buttons()

    def set_reaction_name(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        if row is None:
            return
        name = ui.lineedit_reaction_name.value
        if len(name) == 0: # Reset name if empty input
            name = tw.item(row, COL_RXN_NAME).text()
            ui.lineedit_reaction_name.setText(name)
            return
        tw.item(row, COL_RXN_NAME).setText(name)
        self.current_reaction_name = name # We can only edit the current reaction
        # update ordered dict, keeping order
        keys = list(self.project.reactions.keys())
        keys[row] = name
        self.project.reactions = OrderedDict(zip(keys, self.project.reactions.values()))
        self.set_unsaved_flag()


    def chemistry_restrict_phases(self):
        # Set up comboboxes to only allow homogeneous or 2-phase heterogenous reactions
        ui = self.ui.chemistry
        # Collect all current phase info (reactants and products)
        phases = {}
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            for row in range(tw.rowCount()-1): # Skip 'total'
                cb = tw.cellWidget(row, COL_PHASE)
                phase = cb.currentIndex()
                phases[(tw, row)] = phase
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            # We need to determine, for each combobox item, how many phases would be
            # referenced if the combobox were set to that item
            for row in range(tw.rowCount()-1):
                cb = tw.cellWidget(row, COL_PHASE)
                orig_phase = cb.currentIndex()
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    phases[(tw, row)] = i # consider what setting this combobox would do...
                    if i == orig_phase:
                        enabled = True
                    else:
                        enabled = (len(set(phases.values())) <= 2 # Allow at most 2-phase reactions
                                   # and don't allow phases for which there are no available species
                                   and self.chemistry_find_available_species(side, match_phase=i))
                    set_item_enabled(item, enabled)

                phases[(tw, row)] = orig_phase # ..and set it back


    def chemistry_restrict_species(self):
        # Set up comboboxes to ensure that no species is duplicated as a reactant or product
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            # Collect species info
            species = {}
            for row in range(tw.rowCount()-1): # Skip 'total'
                cb = tw.cellWidget(row, COL_SPECIES)
                species[row] = cb.currentText()
            # For each combobox item, determine whether setting the combobox to that value
            # results in duplicated species
            n_rows = tw.rowCount() - 1 # Skip 'total'
            for row in range(n_rows):
                cb = tw.cellWidget(row, COL_SPECIES)
                orig_species = cb.currentText()
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    species[row] = item.text() # consider what setting this combobox would do...
                    enabled = (len(set(species.values())) == n_rows) # should be as many species as rows
                    set_item_enabled(item, enabled)
                species[row] = orig_species # ... and set it back


    def chemistry_handle_selection(self):
        # selection callback for main table
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        enabled = (row is not None)
        self.reaction_mass_totals = [None, None]
        if enabled:
            self.current_reaction_name = tw.item(row, COL_RXN_NAME).text()
            self.working_reaction = deepcopy(self.project.reactions[self.current_reaction_name])
        else:
            self.current_reaction_name = None
            self.working_reaction = None
        self.reaction_edited = False
        self.chemistry_update_detail_pane()
        ui.scrollarea_detail.ensureVisible(0, 0)
    def chemistry_update_detail_pane(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        reaction = self.working_reaction
        enabled = (row is not None) # and tw.item(row,COL_CHEM_EQ).text() != '') # chem eq is blank when defining new reaction
        ui.toolbutton_delete_reaction.setEnabled(enabled)

        ui.toolbutton_add_reactant.setEnabled(bool(self.chemistry_find_available_species('reactants')))
        ui.toolbutton_add_product.setEnabled(bool(self.chemistry_find_available_species('products')))

        for widget in (ui.label_reaction_name,
                       ui.lineedit_reaction_name,
                       ui.groupbox_reactants,
                       ui.groupbox_products,
                       ui.groupbox_dh):
            widget.setEnabled(enabled)

        # Note, leave checkbox_keyword_stiff_chemistry enabled
        # even if no selection
        #ui.bottom_frame.setEnabled(enabled)

        if not enabled:
            self.chemistry_clear_tables()
            ui.lineedit_reaction_name.clear()
            ui.groupbox_dh.setChecked(False)
            for widget in widget_iter(ui.groupbox_dh):
                if isinstance(widget, LineEdit):
                    widget.clear()
            self.current_reaction_name = None
            return

        tw = ui.tablewidget_reactions
        name = tw.item(row, COL_RXN_NAME).text()
        ui.lineedit_reaction_name.setText(name)
        self.current_reaction_name = name

        def handle_phase(tw, cb, row, idx):
            ui = self.ui.chemistry
            # We have to replace the species widget
            old_species_cb = tw.cellWidget(row,  COL_SPECIES)
            if self.working_reaction is None:
                return
            reaction = self.working_reaction
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            species = self.chemistry_find_available_species(side, idx)
            reaction[side][row][0] = species
            item = make_species_item(tw, row, idx, species)
            tw.setCellWidget(row, COL_SPECIES, item)
            if old_species_cb: # necessary?
                try:
                    old_species_cb.activated.disconnect()
                except:
                    pass
                try:
                    old_species_cb.currentIndexChanged.disconnect()
                except:
                    pass
                old_species_cb.deleteLater()
            self.reaction_edited = True
            self.chemistry_restrict_phases()
            self.chemistry_restrict_species()
            self.chemistry_check_fracdh(reaction) # sets 'phases'
            self.chemistry_update_detail_pane()


        def make_phase_item(tw, row, phase):
            ui = self.ui.chemistry
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            cb = ComboBox()
            phases = [self.fluid_phase_name]
            for name in self.solids.keys():
                phases.append(name)
            for p in phases:
                cb.addItem(p)
            if phase is not None:
                cb.setCurrentIndex(phase)
            cb.currentIndexChanged.connect(lambda idx, tw=tw, cb=cb, row=row: handle_phase(tw, cb, row, idx))
            if side=='reactants': # additional callbacks for pseudo-selection
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_reactant_selection(row))
            else:
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_product_selection(row))
            return cb

        def handle_species(tw, cb, row, idx):
            ui = self.ui.chemistry
            #tw.cellWidget(row, COL_SPECIES).setText('1.0')
            self.reaction_edited = True
            species = tw.cellWidget(row, COL_SPECIES).currentText()
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            #reaction = self.project.reactions[self.current_reaction_name]
            reaction = self.working_reaction
            reaction[side][row][0] = species
            self.chemistry_restrict_phases()
            self.chemistry_restrict_species()
            self.chemistry_update_totals()

        def make_species_item(tw, row, phase, species):
            cb = QComboBox()
            idx = 0
            for s in self.species_of_phase(phase):
                cb.addItem(s)
                if s == species:
                    cb.setCurrentIndex(idx)
                idx += 1
            cb.currentIndexChanged.connect(lambda idx, tw=tw, cb=cb, row=row:
                                           handle_species(tw, cb, row, idx))
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            if side=='reactants': # additional callbacks for pseudo-selection
                cb.activated.connect(lambda idx, row=row:
                                     self.chemistry_handle_reactant_selection(row))
            else:
                cb.activated.connect(lambda idx, row=row:
                                     self.chemistry_handle_product_selection(row))
            return cb

        def handle_coeff(widget, val, args):
            ui = self.ui.chemistry
            #assert not widget.zombie
            tw = ui.tablewidget_reactants if widget.side == 'reactants' else ui.tablewidget_products
            #reaction = self.project.reactions[self.current_reaction_name]
            if not self.working_reaction:
                return
            reaction = self.working_reaction
            val = widget.value
            row = widget.row
            tw.current_row = row
            if val in (None, ''):
                val = 1.0
                widget.setText('1.0')
            if val != reaction[widget.side][widget.row][1]:
                self.reaction_edited = True
                reaction[widget.side][widget.row][1] = val
            self.chemistry_update_totals()
            tw.setCurrentCell(widget.row, COL_COEFF)
            widget.selectAll() # simulate selection

        def make_coeff_item(tw, row, val):
            le = ExtLineEdit()
            #le.zombie = False
            le.parent = self
            le.setMaximumWidth(80) #?
            le.dtype = float
            le.min = 0
            le.setToolTip("Stoichometric coefficient")
            le.key = ''
            le.updateValue('', val)
            le.side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            le.tw = tw
            le.row = row
            le.value_updated.connect(handle_coeff)
            return le

        self.chemistry_clear_tables()
        for side in 'reactants', 'products':
            data = reaction.get(side,[])
            #data = self.project.reactions[name].get(side, [])
            # Add a "total" row, only if there is data
            tw = ui.tablewidget_reactants if side=='reactants' else ui.tablewidget_products
            tw.setRowCount(len(data)+1 if data else 0)
            for (row, (species, coeff)) in enumerate(data):
                phase = self.find_species_phase(species)
                if phase is None:
                    self.error("Species %s not found in any phase" % species)
                tw.setCellWidget(row, COL_PHASE, make_phase_item(tw, row, phase))
                tw.setCellWidget(row, COL_SPECIES, make_species_item(tw, row, phase, species))
                tw.setCellWidget(row, COL_COEFF, make_coeff_item(tw, row, coeff))
            if data:
                item = QtWidgets.QTableWidgetItem('Total mol. weight')
                set_item_noedit(item)
                tw.setItem(row+1, COL_SPECIES, item)
                item = QtWidgets.QTableWidgetItem('0.0')
                set_item_noedit(item)
                font=item.font()
                font.setBold(True)
                item.setFont(font)
                tw.setItem(row+1, COL_COEFF, item)
            self.fixup_chemistry_table(tw)
            row = 0 if len(data)==1 else None # Auto-select single row
            if side == 'reactants':
                self.chemistry_handle_reactant_selection(row)
            else:
                self.chemistry_handle_product_selection(row)

        self.fixup_chemistry_table(ui.tablewidget_reactions, stretch_column=2)#chem eq
        self.chemistry_restrict_phases()
        self.chemistry_restrict_species()
        self.chemistry_update_totals()

        dh = reaction.get('dh')
        ui.lineedit_dh.updateValue('dh', dh)
        enabled = (dh is not None)
        ui.groupbox_dh.setChecked(enabled)
        for w in (ui.label_dh, ui.lineedit_dh, ui.label_dh_units,
                  ui.label_fracdh_1, ui.lineedit_fracdh_1,
                  ui.label_fracdh_2, ui.lineedit_fracdh_2):
            w.setEnabled(enabled)

        phases = self.chemistry_reaction_phases(reaction)
        num_phases = len(phases)
        def phase_name(phase):
            if phase == 0:
                return self.fluid_phase_name
            else:
                return list(self.solids.keys())[phase-1]
        if num_phases == 0:
            for w in (ui.label_fracdh_1, ui.lineedit_fracdh_1,
                      ui.label_fracdh_2, ui.lineedit_fracdh_2):
                w.hide()
        elif num_phases == 1:
            if phases[0] is not None:
                ui.label_fracdh_1.setText('Fraction assigned to %s' % phase_name(phases[0]))
                ui.label_fracdh_1.show()
                ui.lineedit_fracdh_1.show()
                ui.lineedit_fracdh_1.updateValue('fracdh', 1.0)
                ui.lineedit_fracdh_1.setReadOnly(True)
                ui.lineedit_fracdh_1.setEnabled(False)
                for w in ui.label_fracdh_2, ui.lineedit_fracdh_2:
                    w.hide()
        elif num_phases == 2:
            fracdh = reaction.get('fracdh', {})
            if phases[0] is not None:
                fracdh1 = fracdh.get(phases[0])
                ui.label_fracdh_1.setText('Fraction assigned to %s' % phase_name(phases[0]))
                ui.label_fracdh_1.show()
                ui.lineedit_fracdh_1.show()
                ui.lineedit_fracdh_1.setReadOnly(False)
                ui.lineedit_fracdh_1.setEnabled(ui.groupbox_dh.isChecked())
                ui.lineedit_fracdh_1.updateValue('fracdh', fracdh1)
            if phases[1] is not None:
                fracdh2 = fracdh.get(phases[1])
                ui.label_fracdh_2.setText('Fraction assigned to %s' % phase_name(phases[1]))
                ui.label_fracdh_2.show()
                ui.lineedit_fracdh_2.show()
                ui.lineedit_fracdh_2.updateValue('fracdh', fracdh2)
        else:
            self.error("num_phases = %s" % num_phases)


    def format_chem_eq(self, reactants, products):
        tmp = []
        for side in reactants, products:
            tmp.append(' + '.join(species if coeff==1.0 else '%g*%s' % (coeff, species)
                                  for (species, coeff) in side))
        chem_eq = ' --> '.join(tmp)
        return chem_eq


    def chemistry_update_chem_eq(self, name):
        # Allow updating non-selected reaction, used when renaming species
        ui = self.ui.chemistry
        #reaction = self.working_reaction
        reaction = self.project.reactions.get(name)
        if reaction is None:
            return
        chem_eq = self.format_chem_eq(reaction.get('reactants',[]), reaction.get('products',[]))
        display_text = chem_eq.replace('-->', '→')
        if reaction['chem_eq'].upper() == 'NONE': # Update disabled reaction
            self.disabled_reactions[name] = chem_eq
        else: # Update active reaction
            reaction['chem_eq'] = chem_eq
        tw = ui.tablewidget_reactions
        for row in range(tw.rowCount()):
            if tw.item(row, COL_RXN_NAME).text() == name:
                ui.tablewidget_reactions.item(row, COL_CHEM_EQ).setText(display_text)
                break


    def chemistry_update_totals(self):
        ui = self.ui.chemistry
        self.reaction_mass_totals = [None, None]
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            nrows = tw.rowCount()
            if nrows < 2: #
                continue
            tot = 0.0
            for row in range(nrows-1):
                species = tw.cellWidget(row, COL_SPECIES).currentText()
                m_w = self.species_mol_weight(species)
                if m_w is None: # Undefined mol. wt - species_mol_wt printed a warning
                    self.warning("Molecular weight for '%s' not found" % species)
                    continue
                coeff = tw.cellWidget(row, COL_COEFF).value
                if coeff in (None, ''):
                    continue # Empty input field
                tot += m_w * float(coeff)
            self.reaction_mass_totals[tw==ui.tablewidget_products] = tot
            tot = round(tot, 6)
            tw.item(nrows-1, COL_COEFF).setText(str(tot))
        self.chemistry_update_buttons()


    def chemistry_check_reaction_balance(self):
        ui = self.ui.chemistry
        if not self.working_reaction:
            return False, "No reaction"
        reaction = self.working_reaction
        reactants = reaction.get('reactants',[])
        products = reaction.get('products',[])
        if not reactants:
            return False, "No reactants defined"
        if not products:
            return False, "No products defined"
        if sorted(reactants) == sorted(products):
            return False, "Reactants = products"
        totals = self.reaction_mass_totals
        if any (t is None for t in totals):
            return False, "Total mol. weight unavailable"
        mass_reactants, mass_products = totals
        if mass_reactants == 0.0:
            return False, "Reactant mass total = 0"
        if mass_products == 0.0:
            return False, "Product mass total = 0"
        if abs(mass_products/mass_reactants - 1.0) < 1e-4 :
            return True, "Reaction is balanced"
        else:
            return False, "Reaction is unbalanced"


    def chemistry_handle_reactant_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        tw.current_row = row
        enabled = (row is not None)
        ui.toolbutton_delete_reactant.setEnabled(enabled)
        for r in range(tw.rowCount()-1):
            tw.cellWidget(r, COL_COEFF).deselect()
        if enabled:
            tw.setCurrentCell(row, COL_COEFF)
            tw.cellWidget(row, COL_COEFF).selectAll()


    def chemistry_handle_product_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        tw.current_row = row
        enabled = (row is not None)
        ui.toolbutton_delete_product.setEnabled(enabled)
        for r in range(tw.rowCount()-1):
            tw.cellWidget(r, COL_COEFF).deselect()
        if enabled:
            tw.setCurrentCell(row, COL_COEFF)
            tw.cellWidget(row, COL_COEFF).selectAll()


    def chemistry_reaction_phases(self, reaction):
        """return sorted list of phase indices involved in specified reaction"""
        alias_list = [k[0] for k in reaction.get('reactants',[]) + reaction.get('products',[])]
        phases = list(set(map(self.find_species_phase, alias_list)))
        phases.sort()
        return phases


    def chemistry_update_enabled(self):
        #Chemistry pane is disabled if any solids are specified as PIC.
        disabled = False
        # Fewer than 2 species, no chemistry possible
        if len(self.fluid_species) + sum(map(len, self.solids_species.values())) < 2:
            disabled = True
        if any(self.project.get_value('solids_model', args=[i])=='PIC'
               for (i,s) in enumerate(self.solids, 1)):
            disabled = True
        if self.project.reactions: # Don't ever disable pane if reactions are defined
            disabled = False
        self.find_navigation_tree_item("Chemistry").setDisabled(disabled)


    def fixup_chemistry_table(self, tw, stretch_column=1): # species column, for reactant/product tables
        ui = self.ui.chemistry
        hv = QtWidgets.QHeaderView
        if PYQT5:
            resize = tw.horizontalHeader().setSectionResizeMode
        else:
            resize = tw.horizontalHeader().setResizeMode
        ncols = tw.columnCount()

        for n in range(0, ncols):
            resize(n, hv.Stretch if n==stretch_column else hv.ResizeToContents)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # TODO FIXME scrollbar handling is not right - scrollbar status can change
        # outside of this function.  We need to call this everytime window geometry changes
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()
        if nrows==0:
            row_height = 0
            height = header_height+scrollbar_height
        else:
            row_height = tw.rowHeight(0)
            height =  (header_height+scrollbar_height
                       + nrows*row_height + 4) # extra to avoid unneeded scrollbar

        if tw == ui.tablewidget_reactions:
            ui.top_frame.setMaximumHeight(height+24)
            ui.top_frame.setMinimumHeight(header_height+24+row_height*min(nrows,3))
            ui.top_frame.updateGeometry()
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
        else:
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(height)
        tw.updateGeometry() #? needed?


    def chemistry_add_reaction(self):
        ui = self.ui.chemistry
        aliases = list(self.species_all_aliases())
        if not aliases:
            ui.toolbutton_add_reaction.setEnabled(False)
            return

        count = 1
        while 'Reaction_%s' % count in self.project.reactions:
            count += 1
        name = 'Reaction_%s' % count

        self.working_reaction =  {'reactants': [],
                                  'products': [],
                                  'chem_eq': ''}

        self.project.reactions[name] = self.working_reaction
        self.setup_chemistry() # adds row to table
        tw = ui.tablewidget_reactions
        row = tw.rowCount()-1
        tw.cellWidget(row, COL_ENABLE).setChecked(True) #Enable new reaction
        tw.setCurrentCell(row, COL_RXN_NAME) # and select it
        self.reaction_edited = True
        self.chemistry_update_buttons()


    def chemistry_clear_tables(self):
        ui = self.ui.chemistry
        # Work from end, since removing down-shifts entries
        for tw in ui.tablewidget_reactants, ui.tablewidget_products:
            tw.current_row = None
            for row in range(tw.rowCount()-1, -1, -1):
                for col in (2, 1, 0):
                    w = tw.cellWidget(row, col)
                    if w:
                        #w.zombie = True
                        try:
                            w.value_updated.disconnect()
                        except:
                            pass
                        try:
                            w.activated.disconnect()
                        except:
                            pass
                        try:
                            w.currentIndexChanged.disconnect()
                        except:
                            pass
                        tw.removeCellWidget(row, col)
                        w.deleteLater()
            tw.clearContents()
            tw.setRowCount(0)
            self.fixup_chemistry_table(tw)


    def chemistry_delete_reaction(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        if row is None:
            row = get_selected_row(tw)
        if row is None:
            return
        self.chemistry_clear_tables()
        name = tw.item(row, COL_RXN_NAME).text()
        tw.removeRow(row)
        self.fixup_chemistry_table(tw, stretch_column=1)
        self.print_internal(name, font='strikeout')
        del self.project.reactions[name]
        self.set_unsaved_flag()
        #self.setup_chemistry() # handled by selection change


    def chemistry_update_buttons(self):
        # Update the apply/revert/delete buttons
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        balanced, message = self.chemistry_check_reaction_balance()
        if not self.reaction_edited:
            ui.label_status.setText('') #"Reaction unmodified")
        else:
            ui.label_status.setText(message)
        ui.toolbutton_ok.setEnabled(self.reaction_edited and balanced)
        ui.toolbutton_cancel.setEnabled(self.reaction_edited)
        ui.toolbutton_add_reaction.setEnabled(not self.reaction_edited)
        ui.toolbutton_delete_reaction.setEnabled(get_selected_row(tw) is not None)
        tw.setEnabled(not self.reaction_edited)


    def chemistry_apply_changes(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        ui.toolbutton_add_reaction.setEnabled(True)
        ui.toolbutton_delete_reaction.setEnabled(row is not None)
        self.reaction_edited = False
        self.chemistry_update_buttons()
        if row is None or self.working_reaction is None:
            return

        self.project.reactions[self.current_reaction_name] = deepcopy(self.working_reaction)
        self.chemistry_update_chem_eq(self.current_reaction_name)
        self.set_unsaved_flag()


    def chemistry_revert_changes(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        ui.toolbutton_add_reaction.setEnabled(True)
        if row is not None:
            chem_eq = tw.item(row, COL_CHEM_EQ)
            if not (chem_eq and chem_eq.text()):
                self.chemistry_delete_reaction()
        self.reaction_edited = False
        self.chemistry_update_buttons()
        self.chemistry_handle_selection()


    def chemistry_find_available_species(self, side, match_phase=None):
        # side is 'reactants' or 'products'
        # if match_phase is passed, species must belong to that phase
        if not self.current_reaction_name:
            return
        #reaction = self.project.reactions.get(self.current_reaction_name)
        reaction = self.working_reaction
        if reaction is None:
            return

        # Collect phase info
        if match_phase is None:
            phases = self.chemistry_reaction_phases(reaction)
        else:
            phases = []
        used = set(k[0] for k in reaction.get(side, []))

        for alias in self.species_all_aliases():
            # This is a bit inefficient - we should already know the phase
            # (but the species list should not be too long)
            species_phase = self.find_species_phase(alias)
            if match_phase is not None and species_phase != match_phase:
                continue
            if len(phases) > 1 and species_phase not in phases:
                continue
            if alias not in used:
                return alias

    def chemistry_add_reactant(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('reactants')
        if not alias:
            return
        reaction['reactants'].append([alias, 1.0])
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        self.chemistry_handle_reactant_selection(row=len(reaction['reactants'])-1)
        if not self.chemistry_find_available_species('reactants'):
            ui.toolbutton_add_reactant.setEnabled(False)


    def chemistry_delete_reactant(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        row = tw.current_row
        if row is None:
            return
        reaction = self.working_reaction
        # reaction = self.project.reactions[self.current_reaction_name]
        del reaction['reactants'][row]
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_reactants = len(reaction['reactants'])
        if n_reactants == 0:
            row = None
        elif row > n_reactants-1:
            row = n_reactants-1
        if row is not None:
            tw.setCurrentCell(row, COL_COEFF)
        tw.current_row = row
        self.chemistry_handle_reactant_selection(row=row)
        if self.chemistry_find_available_species('reactants'):
            ui.toolbutton_add_reactant.setEnabled(True)
        self.chemistry_update_buttons()


    def chemistry_add_product(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('products')
        if not alias:
            return
        reaction['products'].append([alias, 1.0])
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        self.chemistry_handle_product_selection(row=len(reaction['products'])-1)
        if not self.chemistry_find_available_species('products'):
            ui.toolbutton_add_product.setEnabled(False)


    def chemistry_delete_product(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        row = tw.current_row
        if row is None:
            return
        reaction = self.working_reaction
        del reaction['products'][row]
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_products = len(reaction['products'])
        if n_products == 0:
            row = None
        elif row > n_products-1:
            row = n_products-1
        if row is not None:
            tw.setCurrentCell(row, COL_COEFF)
        tw.current_row = row
        self.chemistry_handle_product_selection(row=row)
        if self.chemistry_find_available_species('products'):
            ui.toolbutton_add_product.setEnabled(True)


    def chemistry_check_species_in_use(self, alias):
        for (rxn_name, reaction) in self.project.reactions.items():
            for side in 'reactants', 'products':
                if any(v[0]==alias for v in reaction.get(side,[])):
                    return rxn_name
        reaction = self.working_reaction
        if reaction:
            for side in 'reactants', 'products':
                if any(v[0]==alias for v in reaction.get(side,[])):
                    return self.current_reaction_name
        return False


    def chemistry_rename_species(self, old_name, new_name):
        reaction = self.working_reaction
        if reaction is not None:
            for side in 'reactants', 'products':
                for v in reaction.get(side, []):
                    if v[0]==old_name:
                        v[0] = new_name

        for (rxn_name, reaction) in self.project.reactions.items():
            changed = False
            for side in 'reactants', 'products':
                for v in reaction.get(side, []):
                    if v[0]==old_name:
                        v[0] = new_name
                        changed = True
            if changed:
                if rxn_name is not None:
                    self.chemistry_update_chem_eq(rxn_name)


    def setup_chemistry(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions

        # Note, because we clear and reconstruct this tab each time
        #  we lose the current selection
        old_selection = get_selected_row(tw)

        def make_item(sval):
            item = QtWidgets.QTableWidgetItem(sval)
            set_item_noedit(item)
            return item

        tw.setRowCount(len(self.project.reactions))
        for row, (name, data) in enumerate(self.project.reactions.items()):
            chem_eq = data.get('chem_eq', '')
            enabled = bool(chem_eq and chem_eq.upper()!='NONE')
            item = QCheckBox()
            item.setChecked(enabled)
            item.clicked.connect(lambda enabled, row=row:
                                 self.chemistry_toggle_reaction(row, enabled))
            tw.setCellWidget(row, COL_ENABLE, item)

            item = make_item(name)
            tw.setItem(row, COL_RXN_NAME, item)

            if not enabled:
                chem_eq = self.disabled_reactions.get(name, '')

            text = chem_eq.replace('==', '→')
            text = text.replace('-->', '→')
            item = make_item(text)
            tw.setItem(row, COL_CHEM_EQ, item)

        # Autoselect if only 1 row
        if tw.rowCount() == 1:
            tw.setCurrentCell(0, COL_RXN_NAME)
        elif old_selection is not None and old_selection < tw.rowCount():
            tw.setCurrentCell(old_selection, COL_RXN_NAME)

        self.chemistry_update_detail_pane()
        self.fixup_chemistry_table(tw, stretch_column=2)#chem eq
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            self.fixup_chemistry_table(tw)


    def chemistry_extract_info(self):
        """extract additional chemistry info after loading project file"""
        for reaction in self.project.reactions.values():
            self.chemistry_check_fracdh(reaction) # sets 'phases' field


    def chemistry_toggle_reaction(self, row, enabled):
        name = list(self.project.reactions.keys())[row]
        reaction = self.project.reactions[name]
        if enabled:
            reaction['chem_eq'] = self.disabled_reactions.pop(name, 'NONE')
        else:
            self.disabled_reactions[name] = reaction.get('chem_eq', 'NONE')
            reaction['chem_eq'] = 'NONE'
        self.set_unsaved_flag()


    def chemistry_check_fracdh(self, reaction):
        phases = self.chemistry_reaction_phases(reaction)
        reaction['phases'] = phases
        num_phases = len(phases)
        dh = reaction.get('dh')
        fracdh = reaction.get('fracdh')
        if dh is None:
            reaction.pop('fracdh', None)
        else:
            if num_phases == 0:
                reaction.pop('fracdh', None)
            elif len(phases) == 1:
                reaction['fracdh'] = {phases[0]: 1.0}
            elif len(phases) == 2:
                if fracdh is None:
                    reaction['fracdh'] = {phases[0]:0.5, phases[1]:0.5}
                elif sorted(fracdh.keys()) == phases:
                    pass # Nothing to do
                elif phases[0] in fracdh:
                    val = fracdh[phases[0]]
                    reaction['fracdh'] = {phases[0]:val, phases[1]:1.0-val}
                elif phases[1] in fracdh:
                    val = fracdh[phases[1]]
                    reaction['fracdh'] = {phases[0]:1.0-val, phases[1]:val}
                else:
                    reaction['fracdh'] = {phases[0]:0.5, phases[1]:0.5}
            else:
                self.error('num_phases = %s' % num_phases)


    def chemistry_to_str(self):
        if self.disabled_reactions:
            data = {'disabled_reactions': self.disabled_reactions}
            return JSONEncoder().encode(data)
        return ''


    def chemistry_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)

        if data:
            RP = ReactionParser()
            val = data.get('disabled_reactions')
            if val:
                self.disabled_reactions = val
                for (name, eq) in self.disabled_reactions.items():
                    if name in self.project.reactions:
                        reaction = self.project.reactions[name]
                        reaction['reactants'], reaction['products'] = RP.parse_chem_eq(eq)
                        self.chemistry_check_fracdh(reaction)


    def reset_chemistry(self):
        ui = self.ui.chemistry
        self.current_reaction_name = None
        self.reaction_edited = False
        self.reaction_mass_totals = [None, None]
        self.working_reaction = None
        self.project.reactions.clear() # done in project.reset()
        self.disabled_reactions.clear()
        self.chemistry_clear_tables()
        tw = ui.tablewidget_reactions
        tw.clearSelection()
        tw.clearContents()
        tw.setRowCount(0)
        self.fixup_chemistry_table(tw)
        tw.current_row = None
        tw.setEnabled(True)
        ui.toolbutton_add_reaction.setEnabled(True)
        ui.toolbutton_delete_reaction.setEnabled(False)

        for tb in (ui.toolbutton_ok, ui.toolbutton_cancel):
            tb.setEnabled(False)


# Documentation/spec

#Enable the stiff chemistry solver
# Selection always available
# Sets keyword STIFF_CHEMISTRY to .TRUE.

#Chemical reaction input is handled different that keyword pair inputs. All homogeneous gas phase
#chemical reactions and all heterogeneous gas-tfm solids reactions are specified between @(RXNS)
#and @(END) reaction block. All heterogeneous gas-dem and gas-pic reactions are specified
#between @(DES_RXNS) and @(DES_END) reaction block.

#Users use the globally unique species aliases to specify the chemical reaction equation. Each
#reaction is specified with a unique reaction identify.

#Specify the reaction identifier (Name)
# Specification always available
# DEFAULT value reactionN (for the Nth reaction)
# Reaction identifiers must be "Fortran compilable"
#  Alphanumeric combinations (no special characters excluding underscores)
#  Limited to 32 characters
#  First character must be a letter
#  No blank spaces

#Specify chemical reaction reactants (table format)
# Use +/- buttons to add/remove reactants
# Column 1 - Select phase for reactant
#  Use drop down list of user defined phase names
#  Reactions are limited to homogeneous or two-phase heterogeneous reactions
#(e.g., species from three separate phases cannot be referenced by any single chemical reaction)
# Column 2 - Select reactant species
#  Use drop down list to show species in selected phase
#  A species can only appear once as a reactant in the same reaction
# Column 3 - Enter stoichiometric coefficient
#  Numerical value (integer or float)
#  Value must be non-negative
#Specify chemical reaction products (table format)
# Use +/- buttons to add/remove products
# Column 1 - Select phase for product
#  Use drop down list of user defined phase names
#Reactions are limited to homogeneous or two-phase heterogeneous reactions
#(e.g., species from three separate phases cannot be referenced by any single
#chemical reaction)
# Column 2 - Select product species
#  Use drop down list to show species in selected phase
#  A species can only appear once as a product in the same reaction
# Column 3 - Enter stoichiometric coefficient
#  Numerical value (integer or float)
#  Value must be non-negative

#Reactant/Product information is combined with stoichiometric coefficients to define the
#chemical reaction as a string.
# Sets reaction construct keyword CHEM_EQ
# Example: CHEM_EQ = "rcoeff1*reactant1 + rcoeff2*reactant2 --> pcoeff1*product1"

# Error check: Mass of the reactants equal mass of the products (within a tolerance, 1.0e-4).
#  abs(mass_products/mass_reactants - 1.0) < 1e-4



#Enable user-defined heat of reaction
# Selection always available
# DEFAULT disabled

#Specify heat of reaction
# Only available if user-defined heat of reaction is enabled
# DEFAULT value 0.0
# Sets reaction construct keyword DH

#Specify HoR fraction assigned to phase
# Only available if user-defined heat of reaction is enabled
# Homogeneous chemical reactions
#  Specification is not available
#  Set reaction construct keyword fracDH(#) to 1.0 where # is the phase index
# Heterogeneous chemical reactions
#  Entry for each phase referenced by the reaction
#  DEFAULT value 0.5 for both entries
#  Sets reaction construct keyword fracDH(#) for each referenced phase

# The user cannot 'save' the reaction if there are errors. After
# saving (adding?) the reaction, the reaction identifier (name) and
# chemical equation are shown in the summary box at the top. A
# chemical reaction is activated/deactivated by checking/unchecking the box. If the user 'deactivates'
# the chemical equation, the CHEM_EQ reaction construct keyword should get set to "NONE."

# NB: user's guide says: Aliases cannot conflict with existing MFIX variable names (e.g., a
#  species alias of MU_g will cause an error when compiling)
