"""classes to monitor MFIX output files & executables"""

import glob
import logging
import os
import subprocess

log = logging.getLogger(__name__)

class Monitor(object):
    """class for monitoring available MFIX executables and output files"""
    def __init__(self, parent):
        self.parent = parent
        self.cache = {}
        self.outputs = None
        self.exes = {}
        self.update_exes()

    def get_exes(self):
        """returns a dict mapping full [mfix|pymfix] paths
        to configuration options."""
        self.update_exes()
        return self.exes

    def update_exes(self):
        """update self.exes"""
        def mfix_print_flags(mfix_exe, cache=self.cache):
            """Determine mfix configuration by running mfix --print-flags.  Cache results"""
            try: # Possible race, file may have been deleted/renamed since isfile check!
                stat = os.stat(mfix_exe)
            except OSError as err:
                log.exception("could not run %s --print-flags", mfix_exe)
                return ''

            cached_stat, cached_flags = cache.get(mfix_exe, (None, None))
            if cached_stat and cached_stat == stat:
                return cached_flags

            exe_dir = os.path.dirname(mfix_exe)
            popen = subprocess.Popen(mfix_exe + " --print-flags",
                                     cwd=exe_dir,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True)
            (out, err) = popen.communicate()
            flags = '' if err else out.strip()
            cache[mfix_exe] = (stat, flags)
            return flags

        exes = {} # map exes -> flags
        for d in [self.parent.get_project_dir(),]:
            for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                exe = os.path.abspath(os.path.join(d, name))
                if os.path.isfile(exe):
                    log.debug("found %s executable in %s", name, d)
                    exes[exe] = str(mfix_print_flags(exe))
        self.exes = exes

    def get_res_files(self): # Why is this in monitor class?
        if not self.parent.get_project_dir():
            return
        pattern = os.path.join(self.parent.get_project_dir(), '*.RES')
        return glob.glob(pattern)

    def get_outputs(self, patterns=[]):
        project_dir = self.parent.get_project_dir()
        if project_dir is None:
            return
        if len(patterns) == 0:
            patterns = [
                '*.LOG', '*.OUT', '*.RES', '*.SP?', 'MFIX.STOP',
                '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        outputs = []
        for pat in patterns:
            outputs.extend(glob.glob(os.path.join(project_dir, pat)))
        return outputs
