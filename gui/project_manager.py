# Project manager class, mediates between mfix gui and mfix project

import logging
log = logging.getLogger(__name__)
import warnings

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

        warn = False
        if widgets_to_update == None:
            widgets_to_update = []
        if widgets_star:
            widgets_to_update.extend(widgets_star)
        if not widgets_to_update:
            warn = True

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
                                   font="Monospace", color='red' if warn else None)


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
            nmax_g = self.get_value('nmax_g', 0)
            if len(self.gasSpecies) != nmax_g:
                warnings.warn("nmax_g = %d, %d gas species defined" %
                              (nmax_g, len(self.gasSpecies)))

            # TODO:  make sure aliases are unique

            # Make sure they are sorted by index before inserting into gui
            self.gasSpecies.sort(cmp=lambda a,b: cmp(a.ind, b.ind))
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
                phase = g.phase.upper()
                species = g.species_g
                alias = g.species_alias_g
                # TODO:  make sure alias is set & unique
                tmp = user_species.get((species, phase))
                if tmp:
                    (tmin, tmax, mol_weight, coeffs, comment) = tmp
                    species_data = {'source': 'User Defined',
                                    'phase': phase,
                                    'alias': alias,
                                    'molecular_weight': mol_weight,
                                    'heat_of_formation': coeffs[14],
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

                if species_data:
                    self.gui.fluid_species[species] = species_data
                else:
                    warnings.warn("%s: not defined" % g)

            # Now submit all remaining keyword updates
            for keyword in kwlist:
                try:
                    self.submit_change(None, {keyword.key: keyword.value},
                                       args=keyword.args, forceUpdate=True)
                except ValueError as e:
                    errlist.append(e)

            # report any errors
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
        '''
        Register a widget with the project manager. The widget must have a
        value_updated signal to connect to.
        '''
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
        new_thermo_data = []
        user_species = set(k for (k,v) in species_dict.items() if v['source'] != 'BURCAT')
        # Keep sections of thermo_data not mentioned in species_dict, for now
        skip = False
        for line in self.thermo_data:
            line = line.rstrip()
            if skip:
                if line:
                    continue
                else:
                    skip = False
            elif line.endswith(' 1'):
                species = line[:18].strip()
                if species in user_species:
                    skip = True
                    continue
            # Avoid repeated blanks
            if not line:
                if new_thermo_data and not new_thermo_data[-1] :
                    continue
            # Keep the line
            new_thermo_data.append(line)
        self.thermo_data = new_thermo_data
        # Now append records for all user species
        for (k, v) in species_dict.items():
            if v['source'] != 'BURCAT':
                self.thermo_data.extend(format_burcat(k,v))

    def objectName(self):
        return 'Project Manager'

def format_burcat(species, data):
    """Return a list of lines in BURCAT.THR format"""
    ## TODO: move this somewhere else
    lines = []

    # First line
    calc_quality = 'B' # appears in col 68, values A-F
    row = (species, 'User Defined', data['phase'].upper(),
           data['tmin'], data['tmax'],
           calc_quality, data['molecular_weight'])
    lines.append('%-24s%-18s0.%c%10.3f%10.3f%3s%10.5f 1' % row)

    # remaining 3 lines
    fmt = '%15.8E' * 5
    a_low = data['a_low']
    a_high = data['a_high']
    row = tuple(a_low[:5])
    lines.append( (fmt%row) + '    2')
    row = tuple(a_low[5:7] + a_high[:3])
    lines.append( (fmt%row) + '    3')
    row = tuple(a_high[3:] + [data['heat_of_formation']])
    lines.append( (fmt%row) + '    4')
    lines.append('') # blank line
    return lines
