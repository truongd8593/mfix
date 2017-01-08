"""classes to monitor MFIX output files & executables"""

import glob
import logging
import os
from mfixgui.constants import RESTART_FILES, SPX_FILES, VTK_FILES, OTHER_FILES

log = logging.getLogger(__name__)

class Monitor(object):
    """class for monitoring available MFIX executables and output files"""
    def __init__(self, parent):
        self.parent = parent
        self.outputs = None

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
            patterns = RESTART_FILES + SPX_FILES + VTK_FILES + OTHER_FILES
        outputs = []
        for pat in patterns:
            outputs.extend(glob.glob(os.path.join(project_dir, pat)))
        return outputs
