# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox, QComboBox)

from qtpy.QtGui import QValidator

from mfixgui.tools.general import (set_item_noedit, set_item_enabled,
                           widget_iter,
                           get_selected_row, get_combobox_item)

from mfixgui.widgets.base import (LineEdit, ComboBox)

class Chemistry(object):
    #Chemistry Task Pane Window: This section allows a user to define chemical reaction input.

    def init_chemistry(self):
        self.chemistry_current_reaction = None
        ui = self.ui.chemistry
        ui.toolbutton_add_reaction.clicked.connect(self.chemistry_add_reaction)
        ui.toolbutton_delete_reaction.clicked.connect(self.chemistry_delete_reaction)
        ui.toolbutton_add_reactant.clicked.connect(self.chemistry_add_reactant)
        ui.toolbutton_add_product.clicked.connect(self.chemistry_add_product)

        # Toolbuttons
        ui.toolbutton_delete_reaction.setEnabled(False) # Need a selection
        ui.toolbutton_delete_reactant.setEnabled(False)
        ui.toolbutton_delete_product.setEnabled(False)
        ui.tablewidget_reactions.itemSelectionChanged.connect(self.chemistry_handle_selection)
        ui.tablewidget_reactants.itemSelectionChanged.connect(self.chemistry_handle_reactant_selection)
        ui.tablewidget_products.itemSelectionChanged.connect(self.chemistry_handle_product_selection)

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
                self.parent.Xlen = len(text)
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
                phase = cb.currentText()
                phases[(tw, row)] = phase
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            # We need to determine, for each combobox item, how many phases would be
            # referenced if the combobox were set to that item
            for row in range(tw.rowCount()-1):
                cb = tw.cellWidget(row, 0)
                orig_phase = cb.currentText()
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    phases[(tw, row)] = item.text() # consider what setting this combobox would do...
                    enabled = len(set(phases.values())) <= 2 # Allow at most 2-phase reactions
                    set_item_enabled(item, enabled)
                    if not enabled:
                        item.setToolTip('Only homogeneous or 2-phase heterogenous reactions supported')
                    else:
                        item.setToolTip(None)
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
                tw.clearContents()
                tw.setRowCount(0)
                self.fixup_chemistry_table(tw)
            ui.lineedit_reaction_name.clear()
            ui.groupbox_heat_of_reaction.setChecked(False)
            for widget in widget_iter(ui.groupbox_heat_of_reaction):
                if isinstance(widget, LineEdit):
                    widget.clear()
            self.chemistry_current_reaction = None
            return

        tw = ui.tablewidget_reactions
        name = tw.item(row,0).text()
        ui.lineedit_reaction_name.setText(name)
        self.chemistry_current_reaction = name

        def handle_phase(tw, row, idx):
            old_item = tw.cellWidget(row, 1)
            item = make_species_item(tw, row, idx, None)
            tw.setCellWidget(row, 1, item)
            w = tw.cellWidget(row,2)
            if w:
                w.setText('1.0') # Default
            if old_item:
                old_item.deleteLater()
            self.chemistry_update_totals()
            self.chemistry_restrict_phases()
            self.set_unsaved_flag()

        def make_phase_item(tw, row, phase):
            cb = ComboBox()
            phases = []
            if not self.fluid_solver_disabled:
                phases.append(self.fluid_phase_name)
            for name in self.solids.keys():
                phases.append(name)
            for p in phases:
                cb.addItem(p)
            cb.setCurrentIndex(phase)
            cb.activated.connect(lambda idx, tw=tw, row=row: handle_phase(tw, row, idx))
            return cb

        def handle_species(tw, row, idx):
            ui = self.ui.chemistry
            if not self.chemistry_current_reaction:
                return
            tw.cellWidget(row,2).setText('1.0')
            species = tw.cellWidget(row, 1).currentText()
            side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            reaction = self.project.reactions[self.chemistry_current_reaction][0]
            reaction[side][row][0] = species
            self.chemistry_restrict_species()
            self.chemistry_update_totals()
            self.chemistry_update_chem_eq()
            self.set_unsaved_flag()

        def make_species_item(tw, row, phase, species):
            cb = QComboBox()
            idx = 0
            for s in self.species_of_phase(phase):
                cb.addItem(s)
                if s == species:
                    cb.setCurrentIndex(idx)
                idx += 1
            cb.activated.connect(lambda idx, tw=tw, row=row: handle_species(tw, row, idx))
            return cb

        def handle_coeff(widget, val, args):
            ui = self.ui.chemistry
            val = widget.value
            if val in (None, ''):
                val = 1.0
                widget.setText('1.0')
            # if val == 0.0: delete row #?
            self.chemistry_update_totals()
            reaction_data = self.project.reactions[self.chemistry_current_reaction][0]
            reaction_data[widget.side][widget.row][1] = val
            self.chemistry_update_chem_eq()
            self.set_unsaved_flag()

        def make_coeff_item(tw, row, val):
            le = LineEdit()
            le.setMaximumWidth(80) #?
            le.dtype = float
            le.min = 0
            le.setToolTip("Stoichometric coefficient")
            le.key = ''
            le.updateValue('', val)
            le.side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            le.row = row
            le.value_updated.connect(handle_coeff)
            return le

        for (side, tw) in (('reactants', ui.tablewidget_reactants),
                           ('products', ui.tablewidget_products)):
            data = self.project.reactions[name][0].get(side, [])
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
        self.fixup_chemistry_table(ui.tablewidget_reactions, stretch_column=1)
        self.chemistry_restrict_phases()
        self.chemistry_restrict_species()
        self.chemistry_update_totals()

    def chemistry_update_chem_eq(self):
        ui = self.ui.chemistry
        if not self.chemistry_current_reaction:
            return
        reaction = self.project.reactions[self.chemistry_current_reaction][0]
        data = {}
        for side in 'reactants', 'products':
            data[side] = ' + '.join(species if coeff==1.0 else '%.4g * %s' % (coeff, species)
                                    for (species, coeff) in reaction[side] if coeff)
        chem_eq = "%s --> %s" % (data['reactants'], data['products'])
        reaction['chem_eq'] = chem_eq
        self.setup_chemistry()


    def chemistry_update_totals(self):
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            nrows = tw.rowCount()
            if nrows < 2: # Should not happen - either 0 rows or a 'total' at end
                return
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
            tot = round(tot, 6)
            tw.item(nrows-1, 2).setText(str(tot))

    def chemistry_handle_reactant_selection(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        row = get_selected_row(tw)
        enabled = (row is not None)
        ui.toolbutton_delete_reactant.setEnabled(enabled)

    def chemistry_handle_product_selection(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        row = get_selected_row(tw)
        enabled = (row is not None)
        ui.toolbutton_delete_product.setEnabled(enabled)

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
        if self.project.reactions: # Don't disable panes if reactions are defined (?)
            disabled = False
        if any(self.project.get_value('solids_model', args=[i])=='PIC'
               for (i,s) in enumerate(self.solids, 1)):
            disabled = True
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

    def chemistry_available_phases(self):
        ui = self.ui.chemistry
        phases = set()
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            n_rows = tw.rowCount() - 1 # skip 'total'
            for row in range(n_rows):
                cb = tw.cellWidget(row, 0)
                text = cb.currentText()
                if text not in phases: # list won't be very long
                    phases.add(text)
        if len(phases) < 2: # all phases are valid
            phases = []
            if not self.fluid_solver_disabled:
                phases.append(self.fluid_phase_name)
            for name in self.solids.keys():
                phases.append(name)
            return phases
        if len(phases) == 2:
            return phases
        else:
            self.error("more than 2 phases selected")
            return []


    def chemistry_available_species(self, side):
        # side is 'reactants' or 'products'
        ui = self.ui.chemistry
        in_use = set()
        tw = ui.tablewidget_reactants if side=='reactants' else ui.tablewidget_products
        n_rows = tw.rowCount()-1
        for row in range(n_rows):
            in_use.add(tw.cellWidget(row,1).currentText())







    def chemistry_add_reaction(self):
        ui = self.ui.chemistry
        aliases = list(self.species_all_aliases())
        if not aliases:
            return
        count = 1
        while 'Reaction_%s' % count in self.project.reactions:
            count += 1
        name = 'Reaction_%s' % count
        # TODO find first unused
        alias = aliases[0]
        chem_eq = '%s --> %s' % (alias, alias)
        self.project.reactions[name] = (
            {'reactants': [[alias, 1.0]],
             'products': [[alias, 1.0]],
             'chem_eq': chem_eq},
            {})

        self.set_unsaved_flag()
        self.setup_chemistry()
        # Auto-select new reaction
        tw = ui.tablewidget_reactions
        tw.setCurrentCell(tw.rowCount()-1, 0)

    def chemistry_delete_reaction(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
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


    def chemistry_add_reactant(self):
        ui = self.ui.chemistry
        aliases = self.species_all_aliases()
        aliases = list(aliases)
        if not aliases:
            return
        #TODO Find first unused
        alias = aliases[0]
        self.project.reactions[self.chemistry_current_reaction][0]['reactants'].append([alias, 1.0])
        self.setup_chemistry()
        self.chemistry_handle_selection() # force update


    def chemistry_delete_reactant(self):
        ui = self.chemistry.ui
        tw = ui.tablewidget_reactants
        row = get_selected_row(tw)
        if row is None:
            return
        tw.deleteRow(row)
        reaction_data =  self.project.reactions[self.chemistry_current_reaction][0]
        del reaction_data['reactants'][row]
        self.update_chem_eq()
        self.chemistry_update_totals()
        #self.setup_chemistry() #?


    def chemistry_add_product(self):
        ui = self.ui.chemistry
        aliases = self.species_all_aliases()
        aliases = list(aliases)
        if not aliases:
            return
        #TODO Find first unused
        alias = aliases[0]
        self.project.reactions[self.chemistry_current_reaction][0]['products'].append([alias, 1.0])
        self.setup_chemistry()
        self.chemistry_handle_selection() # force update

    def chemistry_delete_product(self):
        pass

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

        self.chemistry_handle_selection() # update bottom pane

        self.fixup_chemistry_table(tw, stretch_column=1)
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            self.fixup_chemistry_table(tw)


    def chemistry_extract_phases(self):
        """extract additional chemistry info after loading project file"""
        for (name, (data, indices)) in self.project.reactions.items():
            alias_list = [k[0] for k in data.get('reactants',[])] + [k[0] for k in data.get('products',[])]
            data['num_phases'] = self.chemistry_num_phases(alias_list)


    def reset_chemistry(self):
        self.project.reactions.clear() # done in project.reset()
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_reactions, ui.tablewidget_reactants, ui.tablewidget_products):
            tw.clearContents()
            tw.setRowCount(0)
            self.fixup_chemistry_table(tw)

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
