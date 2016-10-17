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
        ui = self.ui.chemistry
        ui.toolbutton_add.clicked.connect(self.chemistry_add)
        ui.toolbutton_delete.clicked.connect(self.chemistry_delete)
        # TODO implement 'duplicate' (what does this do?)
        ui.toolbutton_delete.setEnabled(False) # Need a selection
        ui.tablewidget_chemistry.itemSelectionChanged.connect(self.handle_chemistry_selection)
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
        tw = ui.tablewidget_chemistry
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


    def handle_chemistry_selection(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_chemistry
        row = get_selected_row(tw)
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)

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
            return

        tw = ui.tablewidget_chemistry
        name = tw.item(row,0).text()
        ui.lineedit_reaction_name.setText(name)

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
            tw.cellWidget(row,2).setText('1.0')
            self.chemistry_update_totals()

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
            if widget.text() == '':
                widget.setText('1.0')
            self.chemistry_update_totals()

        def make_coeff_item(tw, row, val):
            le = LineEdit()
            le.setMaximumWidth(80) #?
            le.dtype = float
            le.min = 0
            le.setToolTip("Stoichometric coefficient")
            le.key = ''
            le.updateValue('', val)
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
        self.fixup_chemistry_table(ui.tablewidget_chemistry, stretch_column=1)
        self.chemistry_update_totals()

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


    def chemistry_num_phases(self, alias_list):
        """determine minimum number of phases required to find all listed species,
        i.e. if this returns 1, reaction is homogeneous"""
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

        if tw == ui.tablewidget_chemistry:
            ui.top_frame.setMaximumHeight(height+40)
            ui.top_frame.setMinimumHeight(header_height+40)
            ui.top_frame.updateGeometry()
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
        else:
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(height)
        tw.updateGeometry() #? needed?


    def chemistry_add(self):
        pass


    def chemistry_delete(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_chemistry
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


    def setup_chemistry(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_chemistry

        # Note, because we clear and reconstruct this tab each time
        #  we lose the current selection
        tw.clearContents()

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

        self.fixup_chemistry_table(tw, stretch_column=1)
        # Autoselect if only 1 row
        if tw.rowCount() == 1:
            tw.setCurrentCell(0, 0)

        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            self.fixup_chemistry_table(tw)
        self.handle_chemistry_selection() # enable/disable inputs


    def chemistry_extract_phases(self):
        """extract additional chemistry info after loading project file"""
        for (name, (data, indices)) in self.project.reactions.items():
            alias_list = [k[0] for k in data.get('reactants',[])] + [k[0] for k in data.get('products',[])]
            data['num_phases'] = self.chemistry_num_phases(alias_list)


    def reset_chemistry(self):
        self.project.reactions.clear() # done in project.reset()
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_chemistry, ui.tablewidget_reactants, ui.tablewidget_products):
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

    #Specify HoR fraction assigned to -phase
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

    #  user's guide says: Aliases cannot conflict with existing MFIX variable names (e.g., a species alias of MU_g will cause an error when compiling

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