"""This file is will collect the .mfixguiinfo files data and save the info
to a single file for the gui. To use, run:

> make thumbnails

from the mfixgui directory"""

import json
import os

info_dict = {}
for path, dirs, files in os.walk('../'):
    # make sure we don't look in the build dir
    if 'build' not in path and '.mfixguiinfo' in files:
        name = os.path.basename(path)
        with open(os.path.join(path, '.mfixguiinfo')) as f:
            info = f.readlines()[0].split(',')
            d = info_dict[name] = {}
            for k, v in zip(('solver', 'geometry', 'chemistry', 'description'), info):
                d[k] = v
        os.remove(os.path.join(path, '.mfixguiinfo'))

with open('./tools/template_data.json', 'w') as f:
    json.dump(info_dict, f)
