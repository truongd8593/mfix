# Methods to deal with solids pic tab, slip off from solids_handler.py
from __future__ import print_function, absolute_import, unicode_literals, division

from tools.general import get_combobox_item, set_item_enabled

class SolidsPIC(object):
    def init_solids_pic(self):
        s = self.ui.solids

    def setup_pic_tab(self):
        # Note - we are doing this setup on first show of this tab, rather
        # than at init or project load time
        s = self.ui.solids
