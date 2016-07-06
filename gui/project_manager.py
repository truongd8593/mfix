"""ProjectManager handles interaction between gui and Project objects,
updating gui widgets when keywords are changed, and vice-versa.

It is a subclass of Project, and the .project member in the main MfixGui
object is actually a ProjectManager, not a Project.

Widgets get associated with keywords via the 'register' method.  Multiple
widgets may map to the same keyword.  (Can a widget register multiple keywords?)
Widgets must support updateValue method and emit value_updated signal (see
widgets in widgets/base.py)

It is important that widgets emit the 'value_updated" signal only on user
interaction, not a programmatic setValue/setText/setChecked/etc, otherwise
a notification loop is possible """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

from collections import OrderedDict
import sys
import traceback
import logging
log = logging.getLogger(__name__)
import warnings

from project import Project, Keyword
from constants import *

from widgets.base import LineEdit # a little special handling needed

from tools.general import (format_key_with_args, unformat_key_with_args,
                           plural, to_text_string)
from tools import read_burcat


class ProjectManager(Project):
    """handles interaction between gui and mfix project"""
    def __init__(self, gui=None, keyword_doc=None):
        Project.__init__(self, keyword_doc=keyword_doc)
        self.gui = gui
        self.keyword_and_args_to_widget = {}
        self.registered_keywords = set()
        self.solver = SINGLE  # default

    def submit_change(self, widget, newValueDict, args=None):

        """Submit a value change, for example
        submitChange(lineEdit, {'run_name':'new run name'}, args)"""

        # Note, this may be a callback from Qt (in which case 'widget' is
        # the widget that the user activated), or from the initial mfix
        # loading (in which case widget == None)

        if args is None:
            args = []
        if isinstance(args, int):
            args = [args]

        # Special argument handling!
        # Currently only "S" = solids phase is defined
        args = list(args) # Copy since we're modifying
        for (i, arg) in enumerate(args):
            if arg == 'S': #selected solid phase
                phase = self.gui.solids_current_phase
                if phase is None:
                    return
                args[i] = phase

        for (key, newValue) in newValueDict.items():
            if isinstance(newValue, dict):
                if args:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=args+[ind])

                else:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=[ind])
            else:
                self._change(widget, key, newValue, args=args)

        self._cleanDeletedItems()


    def _change(self, widget, key, newValue, args=None):
        key = key.lower()
        updatedValue = None
        assert not isinstance(newValue, Keyword)
        if args is None:
            args = []

        try:
            updatedKeyword = self.updateKeyword(key, newValue, args)
            updatedValue = updatedKeyword.value
        except Exception as e:
            self.gui.print_internal("Warning: %s: %s" %
                                       (format_key_with_args(key, args), e),
                                       color='red')
            traceback.print_exception(*sys.exc_info())
            return

        keytuple = tuple([key]+args)
        widgets_to_update = self.keyword_and_args_to_widget.get(keytuple)
        keytuple_star = tuple([key]+['*'])
        widgets_star = self.keyword_and_args_to_widget.get(keytuple_star)
        widgets_S = None
        if args and self.gui.solids_current_phase == args[0]:
            keytuple_S = tuple([key]+['S'])
            widgets_S = self.keyword_and_args_to_widget.get(keytuple_S)

        if widgets_to_update == None:
            widgets_to_update = []
        if widgets_star:
            widgets_to_update.extend(widgets_star)
        if widgets_S:
            widgets_to_update.extend(widgets_S)

        for w in widgets_to_update:
            try:
                w.updateValue(key, updatedValue, args)
                if isinstance(w, LineEdit):
                    w.text_changed_flag = False # Needed?
            except Exception as e:
                ka = format_key_with_args(key, args)
                #log.warn("%s: %s" % (e, ka))
                msg = "Cannot set %s = %s: %s" % (ka, updatedValue, e)
                self.gui.print_internal(msg, color='red')
                traceback.print_exception(*sys.exc_info())
                raise ValueError(msg)


        if updatedValue is None or updatedValue=='': # why are we getting ''?  has to do with
                                                     # validation and "saved_value"
            self.gui.unset_keyword(key, args) # prints msg in window.
        else:
            val_str = to_text_string(updatedValue) # Just used for log message
            #if isinstance(updatedValue, bool):
            #    val_str = '.%s.' % val_str
            self.gui.set_unsaved_flag()
            self.gui.print_internal("%s = %s" % (format_key_with_args(key, args), val_str),
                                    font="Monospace")


    def guess_solver(self):
        """ Attempt to derive solver type, after reading mfix file"""
        mmax = self.get_value('mmax', default=1)
        if mmax == 0:
            return SINGLE
        solids_models = set(self.get_value(['solids_model', n], default='TFM').upper()
                           for n in range(1, mmax+1))
        if solids_models == set(["TFM"]):
            return TFM
        elif solids_models == set(["DEM"]):
            return DEM
        elif solids_models == set(["PIC"]):
            return PIC
        elif solids_models == set(["TFM", "DEM"]):
            return HYBRID
        # mfix settings are inconsistent, warn user.  (Popup here?)
        msg = "Warning, cannot deduce solver type"
        self.gui.print_internal(msg, color='red')
        #log.warn(msg)
        #default
        return SINGLE


    def load_project_file(self, project_file):
        """Load an MFiX project file."""
        # See also gui.open_project
        n_errs = 0
        errlist = []
        with warnings.catch_warnings(record=True) as ws:
            self.parsemfixdat(fname=project_file)
            # emit loaded keys
            # some of these changes may cause new keywords to be instantiated,
            # so iterate over a copy of the list, which may change
            kwlist = list(self.keywordItems())
            # Let's guess the solver type from the file
            self.solver = self.guess_solver()
            # Now put the GUI into the correct state before setting up interface
            self.gui.set_solver(self.solver)

            # deal with species, since they can cause other widgets to instantiate
            nmax_g = self.get_value('nmax_g', default=0)
            # TODO: handle cases which use rrates.f, like tutorials/reactor1b
            # which has 'nmax(0)' instead of 'nmax_g'
            if len(self.gasSpecies) != nmax_g:
                warnings.warn("nmax_g = %d, %d gas species defined" %
                              (nmax_g, len(self.gasSpecies)))

            # Make sure they are sorted by index before inserting into gui
            self.gasSpecies.sort(key=lambda a: a.ind) # override 'sort' in class Project?
            self.solids.sort(key=lambda a:a.ind)

            # TODO: integrate project.gasSpecies with gui.fluid_species
            db = self.gui.species_popup.db

            # Note that parsemfixdat does not modify the THERMO DATA section into
            # Species objects
            user_species = {}
            if self.thermo_data is not None:
                for (species, lines) in self.thermo_data.items():
                    section = ['User defined', 'loaded from %s' % project_file] + lines
                    data = read_burcat.parse_section(section)
                    for (species, phase, tmin, tmax, mol_weight, coeffs, comment) in data:
                        user_species[(species, phase)] = (tmin, tmax, mol_weight, coeffs, comment)


            for g in self.gasSpecies:
                # First look for definition in THERMO DATA section
                phase = g.phase.upper() # phase and species_g are guaranteed to be set
                species = g.get('species_g')

                source = "User Defined"
                if species is None:
                    species = 'Gas %s' % g.ind
                    source = "Auto"
                    #warnings.warn("no species_g for gas %d" % g.ind)
                alias = g.get('species_alias_g', species)
                mw_g = g.get('mw_g', None)
                # Note, we're going to unset mw_g and migrate it into THERMO DATA

                user_def = user_species.get((species, phase))
                # Look for mismatched phase
                if not user_def:
                    for ((s,p),v) in user_species.items():
                        if s == species:
                            # This is all-too-common.  existing mfix files all have 'S' for
                            # phase in thermo data.
                            #warnings.warn("species '%s' defined as phase '%s', expected '%s'"
                            #              % (species, p, phase))
                            user_def = v
                            break
                if user_def:
                    (tmin, tmax, mol_weight, coeffs, comment) = user_def
                    if mw_g is not None:
                        mol_weight = mw_g # mw_g overrides values in THERMO DATA
                    species_data = {'source': source,
                                    'phase': phase,
                                    'alias': alias,
                                    'mol_weight': mol_weight,
                                    'h_f': coeffs[14],
                                    'tmin': tmin,
                                    'tmax': tmax,
                                    'a_low': coeffs[:7],
                                    'a_high': coeffs[7:14]}

                else:
                    # get this from the species popup so we don't have to load
                    # another copy of the database.  the database is owned by
                    # the species popup.
                    species_data = self.gui.species_popup.get_species_data(species, phase)
                    if not species_data:
                        # Look for mismatched phase definition
                        # This is not really 'mismatched'.  Fluid phase may contain solids and v.v.
                        for p in 'GLSC':
                            if p == phase:
                                continue
                            species_data = self.gui.species_popup.get_species_data(species, p)
                            if species_data:
                                #warnings.warn("species '%s' defined as phase '%s', expected '%s'"
                                #              % (species, p, phase))
                                log.info("species '%s' defined as phase '%s', expected '%s'"
                                         % (species, p, phase))
                                break
                    if species_data:
                        species_data['alias'] = alias
                        if mw_g is not None:
                            species_data['mol_weight'] = mw_g
                            source = 'Burcat*' # Modifed mol. weight

                if not species_data:
                    warnings.warn("no definition found for species '%s' phase '%s'" % (species, phase))
                    species_data = {
                        'alias' : alias,
                        'source': 'Auto',
                        'phase': phase,
                        'mol_weight': mw_g or 0,
                        'h_f': 0.0,
                        'tmin':  0.0,
                        'tmax': 0.0,
                        'a_low': [0.0]*7,
                        'a_high': [0.0]*7}

                if species in self.gui.fluid_species:
                    self.gui.print_internal("Copying %s to %s" % (species, alias))
                    species = alias # Create a new species, so we can override mol. weight, etc

                self.gui.fluid_species[species] = species_data

            self.update_thermo_data(self.gui.fluid_species)


            # Build the 'solids' dict that the gui expects.
            #  Note that this is derived from the project.solids SolidsCollection
            #  but is slightly different.
            default_solids_model = ("DEM" if self.solver==DEM
                                   else "TFM" if self.solver==TFM
                                   else "PIC" if self.solver==PIC
                                   else "TFM")


            for s in self.solids:
                name = s.name
                if not s.get('solids_model'):
                    # Make sure all solids have a defined solids model, even if its just the default
                    s['solids_model'] = default_solids_model
                    self.updateKeyword('solids_model', default_solids_model, args=[s.ind])

                self.gui.solids_species[s.ind] = OrderedDict()

                species = list(s.species)
                species.sort(key=lambda a:a.ind)

                for sp in species:
                    # First look for definition in THERMO DATA section
                    phase = sp.phase.upper() # phase and species_g are guaranteed to be set
                    species = sp.get('species_s')
                    source = "User Defined"
                    if species is None:
                        species = 'Solid %s' % sp.ind
                        source = "Auto"
                        #warnings.warn("no species_s for solid %d" % sp.ind)
                    alias = sp.get('species_alias_s', species)
                    mw_s = sp.get('mw_s', None)
                    density = sp.get('ro_xs0',None)
                    # Note, we're going to unset mw_s and migrate it into THERMO DATA

                    # TODO:  make sure alias is set & unique
                    user_def = user_species.get((species, phase))
                    # Hack, look for mismatched phase
                    if not user_def:
                        for ((s_,p),v) in user_species.items():
                            if s_ == species:
                                # This is all-too-common.  existing mfix files all have 'S' for
                                # phase in thermo data.
                                #warnings.warn("species '%s' defined as phase '%s', expected '%s'"
                                #              % (species, p, phase))
                                user_def = v
                                break
                    if user_def:
                        (tmin, tmax, mol_weight, coeffs, comment) = user_def
                        if mw_s is not None:
                            mol_weight = mw_s # mw_s overrides values in THERMO DATA
                        species_data = {'source': source,
                                        'phase': phase,
                                        'alias': alias,
                                        'mol_weight': mol_weight,
                                        'h_f': coeffs[14],
                                        'tmin': tmin,
                                        'tmax': tmax,
                                        'a_low': coeffs[:7],
                                        'a_high': coeffs[7:14]}

                    else:
                        # get this from the species popup so we don't have to load
                        # another copy of the database.  currently the database is
                        # owned by the species popup.
                        species_data = self.gui.species_popup.get_species_data(species, phase)
                        if not species_data:
                            # Look for  mismatched phase definition
                            # This is not really 'mismatched'.  Solids phase may contain fluids and v.v.
                            for p in 'SCLG':
                                if p == phase:
                                    continue
                                species_data = self.gui.species_popup.get_species_data(species, p)
                                if species_data:
                                    #warnings.warn("species '%s' defined as phase '%s', expected '%s'"
                                    #              % (species, p, phase))
                                    log.info("species '%s' defined as phase '%s', expected '%s'"
                                             % (species, p, phase))
                                    break
                        if species_data:
                            species_data['alias'] = alias
                            if mw_s is not None:
                                species_data['mol_weight'] = mw_s
                                source = 'Burcat*' # Modifed mol. weight

                    if not species_data:
                        warnings.warn("no definition found for species '%s' phase '%s'" % (species, phase))
                        species_data = {
                            'alias' : alias,
                            'source': 'Auto',
                            'phase': phase,
                            'mol_weight': mw_s or 0,
                            'h_f': 0.0,
                            'tmin':  0.0,
                            'tmax': 0.0,
                            'a_low': [0.0]*7,
                            'a_high': [0.0]*7}

                    species_data['density'] = density
                    self.gui.solids_species[s.ind][species] = species_data

                self.update_thermo_data(self.gui.solids_species[s.ind])



                solids_data =  {'model': s.get("solids_model"),
                               'diameter': s.get('d_p0'),
                               'density': s.get('ro_s0'),
                                #'species': species
                }
                self.gui.solids[name] = solids_data


            # Now submit all remaining keyword updates, except the ones we're skipping
            skipped_keys = set(['mw_g', 'mw_s'])
            for kw in kwlist:
                if kw.key in skipped_keys:
                    self.gui.print_internal("%s=%s moved to THERMO DATA section" % (
                        format_key_with_args(kw.key, kw.args),
                        kw.value))

                    self.gui.unset_keyword(kw.key, args=kw.args) # print msg in window
                    continue
                try:
                    self.submit_change(None, {kw.key: kw.value}, args=kw.args)
                except ValueError as e:
                    errlist.append(e)

            # report any errors (this should probably be in gui class)
            for w in errlist + ws:
                self.gui.print_internal("Warning: %s" % w.message, color='red')
            n_errs = len(errlist) + len(ws)
            if n_errs:
                self.gui.print_internal("Warning: %s loading %s" %
                                           (plural(n_errs, "error") , project_file),
                                           color='red')
            else:
                self.gui.print_internal("Loaded %s" % project_file, color='blue')

    def register_widget(self, widget, keys=None, args=None):
        """ Register a widget with the project manager. The widget must have a
        value_updated signal to connect to.  If args is not None, widget will
        be updated only when the keyword with matching args is updated.  If
        args=['*'], widget recieves updates regardless of args.

        Special args:  'S' will be substituted with the currently selected
                        solids phase (the one the user is editing)  """


        if args is None:
            args = []
        else:
            args = list(args)

        log.debug('ProjectManager: Registering {} with keys {}, args {}'.format(
            widget.objectName(),
            keys, args))

        # add widget to dictionary of widgets to update
        d = self.keyword_and_args_to_widget
        for key in keys:
            keytuple = tuple([key]+args)
            if keytuple not in d:
                d[keytuple] = []
            d[keytuple].append(widget)

        if hasattr(widget, 'value_updated'):
            widget.value_updated.connect(self.submit_change)
            #widget.value_updated.connect(self.gui.set_unsaved_flag)

        self.registered_keywords = self.registered_keywords.union(set(keys))

    def unregister_widget(self, widget):
        widget.disconnect()

        #for key in widget.keys:
        #    self.registered_keywords.discard(key) #possibly overzealous

        # This is inefficient
        for (k,v) in self.keyword_and_args_to_widget.items():
            if widget in v:
                v.remove(widget)


    def update_thermo_data(self, species_dict):
        """Update definitions in self.thermo_data based on data in species_dict.
        Unmatching entries are not modified"""

        changed = False
        update_dict =  dict((species, format_burcat(species, data))
                            for (species, data) in species_dict.items())
        for (k,v) in update_dict.items():
            if self.thermo_data.get(k) != v:
                changed = True
                break

        if changed:
            self.thermo_data.update(update_dict)
            self.gui.set_unsaved_flag()

    def update_parameters(self, params):
        """parameter values have changed, loop through and update"""
        for p in params:
            if p in self.parameter_key_map:
                key_args = self.parameter_key_map[p]
                for key_arg in key_args:
                    key, args = unformat_key_with_args(key_arg)
                    keyword_obj = self.keywordLookup(key, args)
                    self._change(self, key, keyword_obj.value, args=args)

    def objectName(self):
        return 'Project Manager'


def format_burcat(species, data):
    """Return a list of lines in BURCAT.THR format"""
    ## TODO: move this somewhere else
    lines = []


    calc_quality = 'B' # This appears in column 68, valid values are A-F.
                       # We'll just assign 'B' to fill the space, it's the
                       # most common value.

    # First line
    row = (species, 'User Defined', data['phase'].upper(),
           data['tmin'], data['tmax'],
           calc_quality, data['mol_weight'])
    lines.append('%-24s%-18s0.%c%10.3f%10.3f%3s%10.5f 1' % row)

    # remaining 3 lines
    fmt = '%15.8E' * 5
    a_low = data['a_low']
    a_high = data['a_high']
    row = tuple(a_low[:5])
    lines.append( (fmt%row) + '    2')
    row = tuple(a_low[5:7] + a_high[:3])
    lines.append( (fmt%row) + '    3')
    row = tuple(a_high[3:] + [data['h_f']])
    lines.append( (fmt%row) + '    4')
    lines.append('') # blank line
    return lines
