# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
import copy

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox, QComboBox)

from qtpy.QtGui import QValidator

from mfixgui.tools.general import (set_item_noedit, set_item_enabled,
                           widget_iter,
                           get_selected_row, get_combobox_item)

from mfixgui.widgets.base import (LineEdit, ComboBox)

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

class Chemistry(object):
    #Chemistry Task Pane Window: This section allows a user to define chemical reaction input.

    def init_chemistry(self):
        ui = self.ui.chemistry

        # data members
        self.current_reaction_name = None
        self.reaction_edited = False
        self.working_reaction = None
        self.reaction_mass_totals = [None, None]

        # Toolbuttons
        ui.toolbutton_add_reaction.clicked.connect(self.chemistry_add_reaction)
        ui.toolbutton_delete_reaction.clicked.connect(lambda checked: self.chemistry_delete_reaction())
        ui.toolbutton_add_reactant.clicked.connect(self.chemistry_add_reactant)
        ui.toolbutton_delete_reactant.clicked.connect(self.chemistry_delete_reactant)
        ui.toolbutton_add_product.clicked.connect(self.chemistry_add_product)
        ui.toolbutton_delete_product.clicked.connect(self.chemistry_delete_product)
        ui.toolbutton_save.clicked.connect(self.chemistry_save_reaction)
        ui.toolbutton_cancel.clicked.connect(self.chemistry_cancel)

        ui.toolbutton_delete_reaction.setEnabled(False) # Need a selection
        ui.toolbutton_delete_reactant.setEnabled(False)
        ui.toolbutton_delete_product.setEnabled(False)

        ui.toolbutton_save.setEnabled(False)
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
                    return (QValidator.Acceptable, text, pos)
                else:
                    return (QValidator.Invalid, text, pos)
        ui.lineedit_reaction_name.setValidator(RxnIdValidator(parent=self))

    def set_reaction_name(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        if row is None:
            return
        name = ui.lineedit_reaction_name.value
        if len(name) == 0: # Reset name if empty input
            name = tw.item(row,0).text()
            ui.lineedit_reaction_name.setText(name)
            return
        tw.item(row,0).setText(name)
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
                cb = tw.cellWidget(row, 0)
                phase = cb.currentIndex()
                phases[(tw, row)] = phase
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            # We need to determine, for each combobox item, how many phases would be
            # referenced if the combobox were set to that item
            for row in range(tw.rowCount()-1):
                cb = tw.cellWidget(row, 0)
                orig_phase = cb.currentIndex()
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    phases[(tw, row)] = i # consider what setting this combobox would do...
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
                cb = tw.cellWidget(row, 1)
                species[row] = cb.currentText()
            # For each combobox item, determine whether setting the combobox to that value
            # results in duplicated species
            n_rows = tw.rowCount() - 1 # Skip 'total'
            for row in range(n_rows):
                cb = tw.cellWidget(row, 1)
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
            self.current_reaction_name = tw.item(row,0).text()
            self.working_reaction = copy.deepcopy(self.project.reactions[self.current_reaction_name])
        else:
            self.working_reaction = None
        # trim partly-defined reactions.  note, this causes runtime crashes
        #try:
        #    tw.itemSelectionChanged.disconnect()
        #    for row2 in range(tw.rowCount()-1, -1, -1): # We navigated away from a partly-defined reaction
        #        if row2 != row and tw.item(row2, 1).text() == '':
        #            self.chemistry_delete_reaction(row2)
        #finally:
        #    tw.itemSelectionChanged.connect(self.chemistry_handle_selection)

        self.reaction_edited = False
        self.chemistry_update_detail_pane()


    def chemistry_update_detail_pane(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        enabled = (row is not None) # and tw.item(row,1).text() != '') # chem eq is blank when defining new reaction
        ui.toolbutton_delete_reaction.setEnabled(enabled)

        for widget in (ui.label_reaction_name,
                       ui.lineedit_reaction_name,
                       ui.groupbox_reactants,
                       ui.groupbox_products,
                       ui.groupbox_heat_of_reaction):
            widget.setEnabled(enabled)

        # Note, leave checkbox_keyword_stiff_chemistry enabled
        # even if no selection
        #ui.bottom_frame.setEnabled(enabled)

        if not enabled:
            for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
                for row in range(tw.rowCount()):
                    for col in (0,1,2):
                        widget = tw.cellWidget(row, col)
                        if isinstance(widget, ExtLineEdit):
                            widget.value_updated.disconnect()

                tw.clearContents()
                tw.setRowCount(0)
                self.fixup_chemistry_table(tw)
            ui.lineedit_reaction_name.clear()
            ui.groupbox_heat_of_reaction.setChecked(False)
            for widget in widget_iter(ui.groupbox_heat_of_reaction):
                if isinstance(widget, LineEdit):
                    widget.clear()
            self.current_reaction_name = None
            return

        tw = ui.tablewidget_reactions
        name = tw.item(row,0).text()
        ui.lineedit_reaction_name.setText(name)
        self.current_reaction_name = name

        def handle_phase(tw, row, idx):
            ui = self.ui.chemistry
            if not self.working_reaction:
                return
            old_item = tw.cellWidget(row, 1)
            reaction = self.working_reaction[0]
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            species = self.chemistry_find_available_species(side, idx)
            reaction[side][row][0] = species
            item = make_species_item(tw, row, idx, species)
            tw.setCellWidget(row, 1, item)

            #w = tw.cellWidget(row,2)
            #if w:
            #    w.setText('1.0') # Default
            if old_item:
                old_item.deleteLater()
            self.reaction_edited = True
            self.chemistry_update_totals()
            self.chemistry_restrict_phases()
            self.chemistry_restrict_species()


        def make_phase_item(tw, row, phase):
            ui = self.ui.chemistry
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            cb = ComboBox()
            phases = [self.fluid_phase_name]
            for name in self.solids.keys():
                phases.append(name)
            for p in phases:
                cb.addItem(p)
            cb.setCurrentIndex(phase)
            cb.currentIndexChanged.connect(lambda idx, tw=tw, row=row: handle_phase(tw, row, idx))
            if side=='reactants': # additional callbacks for pseudo-selection
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_reactant_selection(row))
            else:
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_product_selection(row))
            return cb

        def handle_species(tw, row, idx):
            ui = self.ui.chemistry
            if not self.current_reaction_name:
                return
            #tw.cellWidget(row, 2).setText('1.0')
            self.reaction_edited = True
            species = tw.cellWidget(row, 1).currentText()
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            #reaction = self.project.reactions[self.current_reaction_name][0]
            reaction = self.working_reaction[0]
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
            cb.currentIndexChanged.connect(lambda idx, tw=tw, row=row: handle_species(tw, row, idx))
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            if side=='reactants': # additional callbacks for pseudo-selection
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_reactant_selection(row))
            else:
                cb.activated.connect(lambda idx, row=row: self.chemistry_handle_product_selection(row))
            return cb

        def handle_coeff(widget, val, args):
            ui = self.ui.chemistry
            tw = ui.tablewidget_reactants if widget.side == 'reactants' else ui.tablewidget_products
            #reaction_data = self.project.reactions[self.current_reaction_name][0]
            if not self.working_reaction:
                return
            reaction_data = self.working_reaction[0]
            val = widget.value
            row = widget.row
            tw.current_row = row
            if val in (None, ''):
                val = 1.0
                widget.setText('1.0')
            if val == 0.0:
                if widget.side == 'reactants':
                    self.chemistry_delete_reactant()
                else:
                    self.chemistry_delete_product()
                self.reaction_edited = True
            else:
                if val != reaction_data[widget.side][widget.row][1]:
                    self.reaction_edited = True
                reaction_data[widget.side][widget.row][1] = val

            self.chemistry_update_totals()
            tw.setCurrentCell(widget.row, 2)
            widget.selectAll() # simulate selection


        def make_coeff_item(tw, row, val):
            le = ExtLineEdit()
            le.setObjectName('LE-%s-%s' % (tw.objectName(), row))
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

        for side in 'reactants', 'products':
            tw = ui.tablewidget_reactants if side=='reactants' else ui.tablewidget_products
            for row in range(tw.rowCount()):
                for col in (0,1,2):
                    widget = tw.cellWidget(row, col)
                    if isinstance(widget, ExtLineEdit):
                        widget.value_updated.disconnect()

            data = self.working_reaction[0].get(side,[])
            #data = self.project.reactions[name][0].get(side, [])

            tw.clearContents()
            # Add a "total" row, only if there is data
            tw.setRowCount(len(data)+1 if data else 0)
            for (row, (species, coeff)) in enumerate(data):
                phase = self.find_species_phase(species)
                if phase is None:
                    self.error("Species %s not found in any phase" % species)
                    continue
                tw.setCellWidget(row, 0, make_phase_item(tw, row, phase))
                tw.setCellWidget(row, 1, make_species_item(tw, row, phase, species))
                tw.setCellWidget(row, 2, make_coeff_item(tw, row, coeff))
            if data:
                item = QtWidgets.QTableWidgetItem('Total mol. weight')
                set_item_noedit(item)
                tw.setItem(row+1, 1, item)
                item = QtWidgets.QTableWidgetItem('0.0')
                set_item_noedit(item)
                font=item.font()
                font.setBold(True)
                item.setFont(font)
                tw.setItem(row+1, 2, item)
            self.fixup_chemistry_table(tw)
            if len(data) == 1: # autoselect
                if side == 'reactants':
                    self.chemistry_handle_reactant_selection(row=0)
                else:
                    self.chemistry_handle_product_selection(row=0)

        self.fixup_chemistry_table(ui.tablewidget_reactions, stretch_column=1)
        self.chemistry_restrict_phases()
        self.chemistry_restrict_species()
        self.chemistry_update_totals()

    def chemistry_update_chem_eq(self):
        ui = self.ui.chemistry
        if not self.current_reaction_name:
            return
        #reaction = self.project.reactions[self.current_reaction_name][0]
        reaction = self.working_reaction[0]
        fmt = {}
        for side in 'reactants', 'products':
            fmt[side] = ' + '.join(species if coeff==1.0 else '%s*%s' % (coeff, species)
                                   for (species, coeff) in reaction[side] if coeff)
        chem_eq = "%s --> %s" % (fmt['reactants'], fmt['products'])
        display_text = chem_eq.replace('-->', '→')
        reaction['chem_eq'] = chem_eq
        alias_list = [k[0] for k in reaction.get('reactants',[])] + [k[0] for k in reaction.get('products',[])]
        reaction['num_phases'] = self.chemistry_num_phases(alias_list)
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        if row is not None:
            ui.tablewidget_reactions.item(row, 1).setText(display_text)
        #self.setup_chemistry()


    def chemistry_update_totals(self):
        ui = self.ui.chemistry
        self.reaction_mass_totals = [None, None]
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            nrows = tw.rowCount()
            if nrows < 2: #
                continue
            tot = 0.0
            for row in range(nrows-1):
                species = tw.cellWidget(row,1).currentText()
                m_w = self.species_mol_weight(species)
                if m_w is None: # Undefined mol. wt - species_mol_wt sprinted a warning
                    self.warning("Molecular weight for %s not found" % species)
                    continue
                coeff = tw.cellWidget(row,2).value
                if coeff in (None, ''):
                    continue # Empty input field
                tot += m_w * float(coeff)
            self.reaction_mass_totals[tw==ui.tablewidget_products] = tot
            tot = round(tot, 6)
            tw.item(nrows-1, 2).setText(str(tot))
        self.chemistry_update_save_cancel()


    def chemistry_check_reaction_balance(self):
        ui = self.ui.chemistry
        if not self.working_reaction:
            return False
        reaction = self.working_reaction[0]
        reactants = reaction.get('reactants',[])
        products = reaction.get('products',[])
        if sorted(reactants) == sorted(products):
            return False # Reject trivial reaction

        totals = self.reaction_mass_totals
        if any (t is None or t==0.0 for t in totals):
            return False

        mass_reactants, mass_products = totals
        balanced = bool(abs(mass_products/mass_reactants - 1.0) < 1e-4)
        return balanced


    def chemistry_handle_reactant_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        tw.current_row = row
        enabled = (row is not None)
        ui.toolbutton_delete_reactant.setEnabled(enabled)
        for r in range(tw.rowCount()-1):
            tw.cellWidget(r, 2).deselect()
        if enabled:
            tw.setCurrentCell(row, 2)
            tw.cellWidget(row, 2).selectAll()


    def chemistry_handle_product_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        tw.current_row = row
        enabled = (row is not None)
        ui.toolbutton_delete_product.setEnabled(enabled)
        for r in range(tw.rowCount()-1):
            tw.cellWidget(r, 2).deselect()
        if enabled:
            tw.setCurrentCell(row, 2)
            tw.cellWidget(row, 2).selectAll()


    def chemistry_num_phases(self, alias_list):
        """determine minimum number of phases required to find all listed species,
        i.e. if this returns 1, reaction is homogeneous"""
        # This should go away after species/alias integration!
        if not alias_list:
            return 0
        all_aliases = [(0, list(data.get('alias', name) for (name, data) in self.fluid_species.items()))]

        for (num, phase) in enumerate(self.solids_species.values(), 1):
            all_aliases.append((num, list(data.get('alias', name) for (name, data) in phase.items())))

        def recurse(alias_list):
            if len(alias_list) == 0:
                yield []
            elif len(alias_list) == 1:
                alias = alias_list[0]
                for (num, sublist) in all_aliases:
                    if alias in sublist:
                        yield [num]
            else:
                head, tail = alias_list[:1], alias_list[1:]
                for r in recurse(head):
                    for s in recurse(tail):
                        yield r+s

        r = [(len(set(indices)), indices) # Sort by # of distinct elements]
             for indices in recurse(alias_list)]
        return min(r)[0] if r else 0


    def chemistry_update_enabled(self):
        #Chemistry pane is disabled if any solids are specified as PIC.
        disabled = False
        if any(self.project.get_value('solids_model', args=[i])=='PIC'
               for (i,s) in enumerate(self.solids, 1)):
            disabled = True
        if self.project.reactions: # Don't disable panes if reactions are defined (?)
            disabled = False
        self.find_navigation_tree_item("Chemistry").setDisabled(disabled)


    def fixup_chemistry_table(self, tw, stretch_column=1):
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
            height = header_height+scrollbar_height
        else:
            height =  (header_height+scrollbar_height
                       + nrows*tw.rowHeight(0) + 4) # extra to avoid unneeded scrollbar

        if tw == ui.tablewidget_reactions:
            ui.top_frame.setMaximumHeight(height+40)
            ui.top_frame.setMinimumHeight(header_height+40)
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
            return

        count = 1
        while 'Reaction_%s' % count in self.project.reactions:
            count += 1
        name = 'Reaction_%s' % count

        self.working_reaction =  ({'reactants': [],
                                   'products': [],
                                   'chem_eq': ''},
                                  {})
        self.project.reactions[name] = self.working_reaction

        self.reaction_edited = True
        self.setup_chemistry()
        self.chemistry_update_save_cancel()
        # Auto-select new reaction
        tw = ui.tablewidget_reactions
        tw.setCurrentCell(tw.rowCount()-1, 0)

    def chemistry_delete_reaction(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        if row is None:
            row = get_selected_row(tw)
        if row is None:
            return
        name = tw.item(row,0).text()
        tw.removeRow(row)
        self.fixup_chemistry_table(tw, stretch_column=1)
        self.print_internal(name, font='strikeout')
        del self.project.reactions[name]
        self.set_unsaved_flag()
        #self.setup_chemistry() # handled by selection change


    def chemistry_update_save_cancel(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        ui.toolbutton_save.setEnabled(self.reaction_edited and self.chemistry_check_reaction_balance())
        ui.toolbutton_cancel.setEnabled(self.reaction_edited)
        if self.reaction_edited:
            for tb in (ui.toolbutton_add_reaction, ui.toolbutton_delete_reaction):
                tb.setEnabled(False)
        else:
            ui.toolbutton_add_reaction.setEnabled(True)
            ui.toolbutton_delete_reaction.setEnabled(get_selected_row(tw) is not None)

    def chemistry_save_reaction(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        ui.toolbutton_add_reaction.setEnabled(True)
        ui.toolbutton_delete_reaction.setEnabled(row is not None)
        self.reaction_edited = False
        self.chemistry_update_save_cancel()
        if row is None or self.working_reaction is None:
            return
        self.chemistry_update_chem_eq()
        if self.project.reactions[self.current_reaction_name] != self.working_reaction:
            self.project.reactions[self.current_reaction_name] = copy.deepcopy(self.working_reaction)
            self.set_unsaved_flag()


    def chemistry_cancel(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        ui.toolbutton_add_reaction.setEnabled(True)
        if row is not None:
            chem_eq = tw.item(row, 1).text()
            if not chem_eq:  # user cancelled an add
                self.chemistry_delete_reaction()

        self.chemistry_handle_selection()


    def chemistry_find_available_species(self, side, match_phase=None):
        # side is 'reactants' or 'products'
        # if match_phase is passed, species must belong to that phase
        if not self.current_reaction_name:
            return
        #reaction = self.project.reactions.get(self.current_reaction_name)
        reaction = self.working_reaction
        if reaction is None:
            self.error("reaction %s undefined" % self.current_reaction_name)
            return
        reaction_data = reaction[0]

        # Collect phase info
        if match_phase is None:
            phases = set(self.find_species_phase(s)
                         for (s,c) in
                         reaction_data.get('reactants',[])+reaction_data.get('products',[]))
        else:
            phases = set()
        data = reaction_data.get(side, [])
        used = set(s for (s,c) in data)

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
        tw = ui.tablewidget_reactants
        #reaction = self.project.reactions.get(self.current_reaction_name)
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('reactants')
        if not alias:
            return
        reaction[0]['reactants'].append([alias, 1.0])
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        self.chemistry_handle_reactant_selection(row=len(reaction[0]['reactants'])-1)

    def chemistry_delete_reactant(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        row = tw.current_row
        #tw.removeRow(row)
        reaction_data = self.working_reaction[0]
        # reaction_data = self.project.reactions[self.current_reaction_name][0]
        del reaction_data['reactants'][row]
        self.reaction_edited = True
        #self.chemistry_update_totals()
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_reactants = len(reaction_data['reactants'])
        if n_reactants == 0:
            row = None
        elif row > n_reactants-1:
            row = n_reactants-1
        if row is not None:
            tw.setCurrentCell(row, 0)
        tw.current_row = row
        self.chemistry_handle_reactant_selection(row=row)


    def chemistry_add_product(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        #reaction = self.project.reactions.get(self.current_reaction_name)
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('products')
        if not alias:
            return
        reaction[0]['products'].append([alias, 1.0])
        self.reaction_edited = True
        self.chemistry_update_detail_pane()
        self.chemistry_handle_product_selection(row=len(reaction[0]['products'])-1)


    def chemistry_delete_product(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        row = tw.current_row
        #tw.removeRow(row)
        reaction_data = self.working_reaction[0]
        # reaction_data = self.project.reactions[self.current_reaction_name][0]
        del reaction_data['products'][row]
        self.reaction_edited = True
        #self.chemistry_update_totals()
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_products = len(reaction_data['products'])
        if n_products == 0:
            row = None
        elif row > n_products-1:
            row = n_products-1
        if row is not None:
            tw.setCurrentCell(row, 0)
        tw.current_row = row
        self.chemistry_handle_product_selection(row=row)


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
        for row, rowdata in enumerate(self.project.reactions.items()):
            name, (data, indices) = rowdata
            item = make_item(name)
            tw.setItem(row, 0, item)
            text = data.get('chem_eq')
            if text is None:
                continue
            text = text.replace('==', '→')
            text = text.replace('-->', '→')
            item = make_item(text)
            tw.setItem(row, 1, item)

        # Autoselect if only 1 row
        if tw.rowCount() == 1:
            tw.setCurrentCell(0, 0)
        elif old_selection is not None and old_selection < tw.rowCount():
            tw.setCurrentCell(old_selection, 0)

        self.chemistry_update_detail_pane()
        self.fixup_chemistry_table(tw, stretch_column=1)
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            self.fixup_chemistry_table(tw)


    def chemistry_extract_phases(self):
        """extract additional chemistry info after loading project file"""
        for (name, (data, indices)) in self.project.reactions.items():
            alias_list = [k[0] for k in data.get('reactants',[])] + [k[0] for k in data.get('products',[])]
            data['num_phases'] = self.chemistry_num_phases(alias_list)


    def reset_chemistry(self):
        ui = self.ui.chemistry
        self.current_reaction_name = None
        self.reaction_edited = False
        self.reaction_mass_totals = [None, None]
        self.working_reaction = None
        self.project.reactions.clear() # done in project.reset()
        for tw in (ui.tablewidget_reactions, ui.tablewidget_reactants, ui.tablewidget_products):
            tw.clearSelection()
            tw.clearContents()
            tw.setRowCount(0)
            self.fixup_chemistry_table(tw)
            tw.current_row = None
        ui.toolbutton_add_reaction.setEnabled(True)
        ui.toolbutton_delete_reaction.setEnabled(False)

        for tb in (ui.toolbutton_save, ui.toolbutton_cancel):
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
    #  No black spaces

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

    # Error check: Mass of the reactants equal mass of the products (within a tolerance, 1.0e-6).
    # abs[(rcoeff1*MW(reactant1) + rcoeff2*MW(reactant2) - (pcoeff1*product1)] < 1.0e-6

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
    #  DEFAULT value 0.5 for both entires
    #  Sets reaction construct keyword fracDH(#) for each referenced phase

    # The user cannot 'save' the reaction if there are errors. After
    # saving (adding?) the reaction, the reaction identifier (name) and
    # chemical equation are shown in the summary box at the top. A
    # chemical reaction is activated/deactivated by checking/unchecking the box. If the user 'deactivates'
    # the chemical equation, the CHEM_EQ reaction construct keyword should get set to "NONE."

    # NB: user's guide says: Aliases cannot conflict with existing MFIX variable names (e.g., a species alias of MU_g will cause an error when compiling

"""
Evaporation{ chem_eq = "Vapor --> Liquid"}

Ash_to_Ash { chem_eq = 'Ash2 --> FlyAsh' } ! Cold --> Hot


Charring {

chem_eq = 'Biomass --> 8.3334 * Char'

DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}

Pyrolysis {

chem_eq = 'Biomass --> ' &
'0.9639 * CO + 0.8771 * CO2 + 0.3491 * CH4 + ' &
'1.6276 * H2 + 1.4210 * H2O'

DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}

"""

# reaction ID:  32 chars max, alphanumeric + underscore
