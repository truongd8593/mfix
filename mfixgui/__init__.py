# -*- coding: utf-8 -*-

import subprocess
from os import path

def build_in_dir(dirpath):
    here = path.abspath(path.dirname(__file__))
    configure_mfix_path = path.join(path.dirname(here), 'configure_mfix')
    conf_out, conf_err = subprocess.Popen([configure_mfix_path, '--python'], cwd=dirpath ).communicate()
    make_out, make_err = subprocess.Popen(['make',], cwd=dirpath ).communicate()
    return str(conf_out)+str(conf_err), str(make_out)+str(make_err)
