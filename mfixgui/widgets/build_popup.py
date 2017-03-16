# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
import subprocess
import argparse
import logging

from qtpy import QtWidgets, QtCore

log = logging.getLogger('mfix-gui' if __name__ == '__main__' else __name__)

class BuildPopup(QtWidgets.QDialog):
    def __init__(self, parent=None, cwd='./'):
        QtWidgets.QDialog.__init__(self, parent)

        self.error_list = []
        self.cwd = cwd
        self.setModal(False)
        self.setWindowTitle('Build Solver')

        # generate ui
        self.layout = QtWidgets.QGridLayout(self)
        self.layout.setSizeConstraint(self.layout.SetFixedSize)

        self.progressbar = QtWidgets.QProgressBar()
        self.layout.addWidget(self.progressbar, 10, 0, 1, -1)

        self.build_btn = QtWidgets.QPushButton('Build Solver')
        self.build_btn.clicked.connect(self.build)
        self.layout.addWidget(self.build_btn, 11, 1)

        self.cancel_btn = QtWidgets.QPushButton('Cancel')
        self.cancel_btn.clicked.connect(self.cancel)
        self.cancel_btn.setEnabled(False)
        self.layout.addWidget(self.cancel_btn, 11, 2)

        self.output_btn = QtWidgets.QPushButton('Show Output')
        self.output_btn.clicked.connect(self.show_output)
        self.layout.addWidget(self.output_btn, 11, 3)

        self.close_btn = QtWidgets.QPushButton('Close')
        self.close_btn.clicked.connect(self.close)
        self.layout.addWidget(self.close_btn, 11, 4)

        self.output = QtWidgets.QTextBrowser()
        self.layout.addWidget(self.output, 100, 0, 1, -1)
        self.output.hide()

    def cancel(self):
        self.build_proc.kill()

    def show_output(self):
        if not self.output.isVisible():
            self.output.show()
        else:
            self.output.hide()

    def build(self):
        self.line_count = 0
        self.build_proc = QtCore.QProcess()
        self.build_proc.setWorkingDirectory(self.cwd)
        self.build_proc.start('build_mfixsolver')
        self.build_proc.readyReadStandardOutput.connect(self.check_progress)
        self.build_proc.finished.connect(self.finished_building)
        self.build_proc.error.connect(self.error)

        self.cancel_btn.setEnabled(True)

    def finished_building(self):
        self.progressbar.setValue(100)

    def error(self, error):
        if error == QtCore.QProcess.FailedToStart:
            msg = "Process failed to start "+ 'build_mfixsolver'
        elif error == QtCore.QProcess.Crashed:
            msg = "Process exit " + 'build_mfixsolver'
        elif error == QtCore.QProcess.Timedout:
            msg = "Process timeout "+ 'build_mfixsolver'
        elif error in (QtCore.QProcess.WriteError, QtCore.QProcess.ReadError):
            msg = "Process communication error " + 'build_mfixsolver'
        else:
            msg = "Unknown error " + 'build_mfixsolver'
        log.warn(msg)
        # make the message print in red
        self.parent.stderr_signal.emit(msg)

    def check_progress(self):

        out = bytes(self.build_proc.readAllStandardOutput()).decode('utf-8')
        if 'error' in out:
            self.error_list.append(out)
        if not out:
            self.check_timer.stop()
        self.line_count += 1
        self.progressbar.setValue(int(self.line_count/3000.0*100))

        cursor = self.output.textCursor()
        cursor.movePosition(cursor.End)

        scrollbar = self.output.verticalScrollBar()
        scrolled_to_end = (scrollbar.value() == scrollbar.maximum())
        cursor.insertText(out)
        if scrolled_to_end:
            scrollbar.setValue(scrollbar.maximum())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cwd', action='store', nargs='?', default='./',
                        help='')
    args = parser.parse_args()

    app = QtWidgets.QApplication([])

    popup = BuildPopup(None, args.cwd)

    popup.exec_()
