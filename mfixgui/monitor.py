"""classes to monitor MFiX output files & executables"""

import glob
import os
from mfixgui.constants import RESTART_FILES, SPX_FILES, VTK_FILES, OTHER_FILES

class Monitor(object):
    """class for monitoring available MFiX executables and output files"""
    def __init__(self, parent):
        self.parent = parent
        self.outputs = None

    def get_res_files(self): # Why is this in monitor class?
        """ get the residual files of an MFIX solver job """
        if not self.parent.get_project_dir():
            return
        pattern = os.path.join(self.parent.get_project_dir(), '*.RES')
        return glob.glob(pattern)

    def get_outputs(self, patterns=None):
        """ get the output files of an MFIX solver job """
        project_dir = self.parent.get_project_dir()
        if project_dir is None:
            return
        if patterns is None or len(patterns) == 0:
            patterns = RESTART_FILES + SPX_FILES + VTK_FILES + OTHER_FILES
        outputs = []
        for pat in patterns:
            outputs.extend(glob.glob(os.path.join(project_dir, pat)))
        return outputs
