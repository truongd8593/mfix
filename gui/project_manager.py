# Project manager class, mediates between mfix gui and mfix project

import logging
log = logging.getLogger(__name__)
import warnings

from collections import OrderedDict

from project import Project, Keyword
from constants import *

from tools.general import (format_key_with_args, plural)
from tools import read_burcat

# --- Project Manager ---
class ProjectManager(Project):
    """handles interaction between gui and mfix project"""

    def __init__(self, gui=None, keyword_doc=None):
        Project.__init__(self, keyword_doc=keyword_doc)
        self.gui = gui
        self.keyword_and_args_to_widget = {}
        self.registered_keywords = set()
        self._widget_update_stack = [] # prevent circular updates
        self.solver = SINGLE # default


    def submit_change(self, widget, newValueDict, args=None,
                      forceUpdate=False): # forceUpdate unused (?)
        """Submit a value change, for example
        submitChange(lineEdit, {'run_name':'new run name'}, args)"""

        # Note, this may be a callback from Qt (in which case 'widget' is
        # the widget that the user activated), or from the initial mfix
        # loading (in which case widget == None)

        if isinstance(args, int):
            args = [args]

        for key, newValue in newValueDict.items():
            if isinstance(newValue, dict):
                if args:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=args+[ind],
                                     forceUpdate=forceUpdate)
                else:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=[ind],
                                     forceUpdate=forceUpdate)
            else:
                self._change(widget, key, newValue, args=args,
                             forceUpdate=forceUpdate)

        self._cleanDeletedItems()


    def _change(self, widget, key, newValue, args=None, forceUpdate=False):
        # prevent circular updates, from this widget or any higher in stack
        if widget in self._widget_update_stack:
            return

        key = key.lower()
        updatedValue = None
        if isinstance(newValue, Keyword): # why needed?
            keyword = newValue
            newValue = keyword.value
        else:
            keyword = None
        if args is None:
            args = []

        try:
            updatedValue = self.updateKeyword(key, newValue, args)
        except Exception as e:
            self.gui.print_internal("Warning: %s: %s" %
                                       (format_key_with_args(key, args), e),
                                       color='red')
            return


        keytuple = tuple([key]+args)
        widgets_to_update = self.keyword_and_args_to_widget.get(keytuple)
        keytuple_star = tuple([key]+['*'])
        widgets_star = self.keyword_and_args_to_widget.get(keytuple_star)

        if widgets_to_update == None:
            widgets_to_update = []
        if widgets_star:
            widgets_to_update.extend(widgets_star)

        # Are we using the 'all' mechanism?
        #widgets_to_update.extend(
        #    self.keyword_and_args_to_widget.get("all", []))

        for w in widgets_to_update:
            # Prevent circular updates
            self._widget_update_stack.append(w)
            try:
                w.updateValue(key, updatedValue, args)
            except Exception as e:
                ka = format_key_with_args(key, args)
                #log.warn("%s: %s" % (e, ka))
                msg = "Cannot set %s = %s: %s" % (ka, updatedValue, e)
                if widget: # We're in a callback, not loading
                    self.gui.print_internal(msg, color='red')
                raise ValueError(msg)

            finally:
                self._widget_update_stack.pop()

        self.gui.print_internal("%s = %s" % (format_key_with_args(key, args),
                                                updatedValue),
                                font="Monospace",
                                color='green' if widgets_to_update else None)


    def guess_solver(self):
        """ Attempt to derive solver type, after reading mfix file"""
        keys = self.keywordItems()
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
        log.warn(msg)
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

            # TODO:  make sure aliases are unique

            # Make sure they are sorted by index before inserting into gui
            self.gasSpecies.sort(key=lambda a: a.ind) # override 'sort' in class Project?
            self.solids.sort(key=lambda a:a.ind)

            # TODO: integrate project.gasSpecies with gui.fluid_species
            db = self.gui.species_popup.db

            # Note that parsemfixdat does not modify the THERMO DATA section into
            # Species objects
            user_species = {}
            if self.thermo_data is not None:
                thermo_data = self.thermo_data[:]
                thermo_data.append('') # slight hack, add blank line to force parsing last block
                section = []
                for line in thermo_data:
                    line = line.strip()
                    if not line: # sections separated by blank lines
                        if not section:
                            continue # multiple blank lines
                        # slight hack: read_burcat expects a CAS ID and a comment block
                        # Note, currently this comment is not exposed to the user ... only
                        #  comments from BURCAT.THR wind up in the gui.
                        section.insert(0, 'User defined,')
                        section.insert(1, 'loaded from %s' % project_file)
                        data = read_burcat.parse_section(section)
                        for (species, phase, tmin, tmax, mol_weight, coeffs, comment) in data:
                            user_species[(species, phase)] = (tmin, tmax, mol_weight, coeffs, comment)
                        section = []
                    else:
                        section.append(line)


            for g in self.gasSpecies:
                # First look for definition in THERMO DATA section
                phase = g.phase.upper() # phase and species are guaranteed to be set
                species = g.get('species_g')
                source = "User Defined"
                if species is None:
                    species = 'Gas %s' % g.ind
                    source = "Auto"
                    #warnings.warn("no species_g for gas %d" % g.ind)
                alias = g.get('species_alias_g', species)
                mw_g = g.get('mw_g', None)
                # Note, we're going to unset mw_g and migrate it into THERMO DATA

                # TODO:  make sure alias is set & unique
                user_def = user_species.get((species, phase))
                # Hack, look for mismatched phase
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
                    # another copy of the database.  currently the database is
                    # owned by the species popup.
                    species_data = self.gui.species_popup.get_species_data(species, phase)
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

                self.gui.fluid_species[species] = species_data

            self.update_thermo_data(self.gui.fluid_species)

            for s in self.solids:
                name = s.name
                species = list(s.species)
                species.sort(key=lambda a:a.ind)
                # FIXME just use the solid object instead of doing this translation
                solids_data =  {'model': s.get("solids_model"),
                               'diameter': s.get('d_p0'),
                               'density': s.get('ro_s0'),
                               'species': species}
                self.gui.solids[name] = solids_data

            # Now submit all remaining keyword updates, except the ones we're skipping
            skipped_keys = set(['mw_g'])
            for kw in kwlist:
                if kw.key in skipped_keys:
                    # is this a warn or info?
                    log.warn("%s=%s moved to THERMO DATA section",
                             format_key_with_args(kw.key, kw.args),
                             kw.value)

                    self.gui.unset_keyword(kw.key, args=kw.args)
                    continue
                try:
                    self.submit_change(None, {kw.key: kw.value},
                                       args=kw.args, forceUpdate=True)
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
        args=['*'], widget recieves updates regardless of args. """

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

        widget.value_updated.connect(self.submit_change)

        self.registered_keywords = self.registered_keywords.union(set(keys))


    def update_thermo_data(self, species_dict):
        """Update definitions in self.thermo_data based on data in species_dict.
        Unmatching entries are not modified"""

        new_thermo_data = []
        #species_to_save = set(k for (k,v) in species_dict.items() if v['source'] != 'BURCAT')
        # We're going to save all of them, even the ones from BURCAT. so that mfix does
        #  not need mfix.dat at runtime.  Note, that means "source" will be "User Decoolfined"
        #  next time this project is loaded
        species_to_save = set(species_dict.keys())
        # Keep sections of thermo_data not mentioned in species_dict, for now
        replace_entry = False
        for line in self.thermo_data:
            line = line.rstrip()
            if replace_entry:
                if line:
                    continue
                else:
                    skip = False
            elif line.endswith(' 1'):
                species = line[:18].strip()
                if species in species_to_save:
                    species_to_save.remove(species) # We've handled this one
                    replace_entry = True
                    continue
            # Avoid repeated blanks
            if not line:
                if new_thermo_data and not new_thermo_data[-1] :
                    continue
            # Keep the line
            new_thermo_data.append(line)
        self.thermo_data = new_thermo_data
        # Now append records for species not handled yet
        for species in species_to_save:
            data = species_dict[species]
            #if data['source'] != 'BURCAT':
            #    self.thermo_data.extend(format_burcat(species,data))
            self.thermo_data.extend(format_burcat(species,data))

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
