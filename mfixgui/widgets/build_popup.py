# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
import argparse
import logging
import os

from qtpy import QtWidgets, QtCore, QtGui

log = logging.getLogger('mfix-gui' if __name__ == '__main__' else __name__)

BUILD_CMD = 'build_mfixsolver'
LINECOUNT = 3000

class BuildPopup(QtWidgets.QDialog):
    def __init__(self, parent=None, cwd='./'):
        QtWidgets.QDialog.__init__(self, parent)

        self.error_list = []
        self.cwd = cwd
        self.build_proc = None
        self.setModal(False)
        self.setWindowTitle('Build Solver')

        # generate ui
        self.layout = QtWidgets.QGridLayout(self)
        self.layout.setSizeConstraint(self.layout.SetFixedSize)

        # don't show options on windows
        visible = os.name != 'nt'
        self.dmp = QtWidgets.QCheckBox('Distributed memory parallel (DMP)')
        self.dmp.setVisible(visible)
        self.layout.addWidget(self.dmp , 2, 0, 1, -1)

        self.smp = QtWidgets.QCheckBox('Shared memory parallel (SMP)')
        self.smp.setVisible(visible)
        self.layout.addWidget(self.smp , 3, 0, 1, -1)

        label = QtWidgets.QLabel('FCFLAGS')
        label.setVisible(visible)
        self.layout.addWidget(label, 5, 0)

        self.fc_flags = QtWidgets.QLineEdit()
        self.fc_flags.setVisible(visible)
        self.layout.addWidget(self.fc_flags, 5, 1, 1, -1)

        label = QtWidgets.QLabel('Other Flags')
        label.setVisible(visible)
        self.layout.addWidget(label, 6, 0)

        self.other_flags = QtWidgets.QLineEdit()
        self.other_flags.setVisible(visible)
        self.layout.addWidget(self.other_flags, 6, 1, 1, -1)

        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setValue(0)
        self.layout.addWidget(self.progressbar, 10, 0, 1, -1)

        self.output_btn = QtWidgets.QPushButton('Show Output')
        self.output_btn.clicked.connect(self.toggle_output)
        self.layout.addWidget(self.output_btn, 11, 1)

        self.build_btn = QtWidgets.QPushButton('Build Solver')
        self.build_btn.clicked.connect(self.build)
        self.layout.addWidget(self.build_btn, 11, 2)

        self.cancel_btn = QtWidgets.QPushButton('Cancel')
        self.cancel_btn.clicked.connect(self.cancel)
        self.layout.addWidget(self.cancel_btn, 11, 3)

        self.output = QtWidgets.QTextBrowser()
        self.layout.addWidget(self.output, 100, 0, 1, -1)
        self.output.hide()

    def cancel(self):
        if self.build_proc is not None:
            self.build_proc.kill()
        self.close()

    def toggle_output(self):
        if not self.output.isVisible():
            self.output.show()
        else:
            self.output.hide()

    def build(self):
        self.line_count = 0
        self.build_proc = QtCore.QProcess()
        self.build_proc.setWorkingDirectory(self.cwd)

        cmd = BUILD_CMD

        fc = self.fc_flags.text()
        if fc:
            cmd += ' FCFLAGS="%s"' % fc

        of = self.other_flags.text()
        if of:
            cmd += ' %s' % of

        if self.smp.isChecked():
            cmd += ' --smp'

        if self.dmp.isChecked():
            cmd += ' --dmp'

        self.print_to_output('Command: %s' % cmd)
        self.build_proc.start(cmd)
        self.build_proc.readyReadStandardOutput.connect(self.check_progress)
        self.build_proc.readyReadStandardError.connect(self.read_err)
        self.build_proc.finished.connect(self.finished_building)
        self.build_proc.error.connect(self.error)
        self.cancel_btn.setText('Cancel')
        self.build_btn.setEnabled(False)

    def finished_building(self):
        self.progressbar.setValue(100)
        self.cancel_btn.setText('Close')

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

        self.print_to_output(msg, error=True)
        self.output.show()

    def read_err(self):
        err_str = bytes(self.build_proc.readAllStandardError()).decode('utf-8')
        self.print_to_output(err_str, error=True)
        self.output.show()

    def check_progress(self):

        out = bytes(self.build_proc.readAllStandardOutput()).decode('utf-8')
        error = True if 'error' in out else False

        if not out:
            self.check_timer.stop()
        self.line_count += 1
        self.progressbar.setValue(int(self.line_count/LINECOUNT*100))
        self.print_to_output(out, error)

    def print_to_output(self, msg, error = False):
        cursor = self.output.textCursor()
        cursor.movePosition(cursor.End)

        scrollbar = self.output.verticalScrollBar()
        scrolled_to_end = (scrollbar.value() == scrollbar.maximum())

        # formating

        char_format = QtGui.QTextCharFormat()
        if error:
            char_format.setForeground(QtGui.QColor('red'))
        cursor.setCharFormat(char_format)
        if not msg.endswith('\n'):
            msg += '\n'
        cursor.insertText(msg)
        if scrolled_to_end:
            scrollbar.setValue(scrollbar.maximum())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cwd', action='store', nargs='?', default='./', help='')
    args = parser.parse_args()

    app = QtWidgets.QApplication([])

    popup = BuildPopup(None, args.cwd)

    popup.exec_()
