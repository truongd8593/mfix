# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
import subprocess
import argparse

from qtpy import QtWidgets, QtCore

class BuildPopup(QtWidgets.QProgressDialog):
    def __init__(self, parent=None, cwd='./'):
        QtWidgets.QProgressDialog.__init__(self, parent)

        self.error_list = []
        self.cwd = cwd

        self.setModal(False)

        self.check_timer = QtCore.QTimer()
        self.check_timer.timeout.connect(self.check_progress)
        self.check_timer.start(0)

        self.build()

    def build(self):
        self.line_count = 0
        self.build_proc = subprocess.Popen(
            'build_mfixsolver',
            cwd=self.cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=True)

    def check_progress(self):

        out = self.build_proc.stdout.readline()
        if 'error' in out:
            self.error_list.append(out)
        if not out:
            self.check_timer.stop()
            self.close()
        self.line_count += 1
        self.setValue(int(self.line_count/3000.0*100))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cwd', action='store', nargs='?', default='./',
                        help='')
    args = parser.parse_args()

    app = QtWidgets.QApplication([])

    popup = BuildPopup(None, args.cwd)

    popup.exec_()
