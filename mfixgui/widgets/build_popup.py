# -*- coding: utf-8 -*-
'''
Dialog to build the mfixsolver from the GUI
On windows, the following paths need to be added to PATH:
  - ANACONDA_HOME/Library/mingw-w64/bin
  - ANACONDA_HOME/Library/usr/bin
'''
from __future__ import print_function, absolute_import, unicode_literals, division
import argparse
import logging
import os
import subprocess
from distutils import spawn

# FIXME: use six instead
try:
    # Python 3
    import urllib.request as urlparse
    import urllib.request as urllib
except ImportError:
    # Python 2
    import urlparse
    import urllib

from qtpy import QtWidgets, QtCore, QtGui

try:
    from mfixgui.tools.util import (
        SCRIPT_DIRECTORY,
    )
except ImportError:
    SCRIPT_DIRECTORY='./'

log = logging.getLogger('mfix-gui' if __name__ == '__main__' else __name__)

BUILD_CMD = 'build_mfixsolver'
LINECOUNT = 3000
WINDOWS_CONDA_PACKAGES = ['m2-base', 'm2-autoconf', 'm2-automake-wrapper',
    'm2-make', 'm2-tar', 'm2w64-gcc', 'm2w64-gcc-fortran']

def path2url(path):
    """Convert path to url."""
    return urlparse.urljoin(
        'file:', urllib.pathname2url(path))

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

        # don't show options on windows
        label = QtWidgets.QLabel('Compiler Options:')
        self.layout.addWidget(label, 0, 0, 1, -1)

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
        self.fc_flags.setToolTip(
            'Flags to be passed to the Fortran compiler.')
        self.fc_flags.setVisible(visible)
        self.layout.addWidget(self.fc_flags, 5, 1, 1, -1)

        label = QtWidgets.QLabel('Other Flags' if visible else 'Flags')
        self.layout.addWidget(label, 6, 0)

        self.other_flags = QtWidgets.QLineEdit()
        self.other_flags.setToolTip(
            'Flags to be passed to build_mfixsolver such as -j to build in parallel.')
        self.layout.addWidget(self.other_flags, 6, 1, 1, -1)

        spacer = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum,)
        self.layout.addItem(spacer, 8, 0)

        self.compile_label = QtWidgets.QLabel('Press "Build Solver" to compile.')
        self.layout.addWidget(self.compile_label, 9, 0, 1, -1)

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

    def show(self):
        # Check to see if compiling is possible
        # Check if BUILD_CMD in path
        possible = True
        if spawn.find_executable(BUILD_CMD) is None:
            possible = False
            message = 'The command: "{}" does not exist in your path.'.format(BUILD_CMD)
        else:
            # check platform dependencies
            if os.name == 'nt':
                possible, message = self.check_windows()
            else:
                possible, message = self.check_unix()

        if not possible:
            setup_guide = path2url(os.path.join(SCRIPT_DIRECTORY, 'doc', 'SETUP_GUIDE.html'))

            message = '\n'.join([
                'Can not build solver in current environment:'
                '<blockquote>' + message + '</blockquote>',
                'Please see the <a href="%s">Setup Guide</a> for more information.' % setup_guide])
            self.parent().warn(message, popup=True)

            self.cancel()
            self.finished.emit(1)
        else:
            QtWidgets.QDialog.show(self)

    def check_unix(self):
        return True, 'blahh' #TODO check unix

    def check_windows(self):
        """check for packages on windows"""

        try:
            packages = subprocess.Popen("conda list", stdout=subprocess.PIPE).stdout.read()
        except FileNotFoundError:
            # conda does not exists
            return False, 'The command "conda" is not avaliable, can not check dependencies.'
        packages = str(packages)

        found = [False]*len(WINDOWS_CONDA_PACKAGES)
        for i, package in enumerate(WINDOWS_CONDA_PACKAGES):
            if package in packages:
                found[i] = True

        missing = '\n'.join([package for i, package in enumerate(WINDOWS_CONDA_PACKAGES) if not found[i]])
        found = all(found)

        return found, '' if found else 'The following packages are missing:\n' + missing

    def cancel(self):
        if self.build_proc is not None:
            self.build_proc.kill()
        self.close()

    def toggle_output(self):
        if not self.output.isVisible():
            self.output.show()
            constraint = self.layout.SetMaximumSize
        else:
            self.output.hide()
            constraint = self.layout.SetFixedSize
        self.layout.setSizeConstraint(constraint)

    def get_environment(self):
        env = QtCore.QProcessEnvironment.systemEnvironment()

        if os.name == 'nt':
            # find conda home
            conda = spawn.find_executable('conda')
            if conda is None:
                self.parent().warn('Can not find "conda" to prepend PATH')
            else:
                anaconda_home = os.path.dirname(os.path.dirname(conda))
                path = env.value('PATH')

                # prepend
                #  - ANACONDA_HOME/Library/mingw-w64/bin
                #  - ANACONDA_HOME/Library/usr/bin
                new_path = os.pathsep.join([
                    os.path.join(anaconda_home, 'Library', 'mingw-w64', 'bin'),
                    os.path.join(anaconda_home, 'Library', 'usr', 'bin'),
                    path
                    ])
                env.insert("PATH", new_path)
        return env

    def build(self):
        self.compile_label.setText('Building Solver...')

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
        self.build_proc.setProcessEnvironment(self.get_environment())
        self.build_proc.readyReadStandardOutput.connect(self.check_progress)
        self.build_proc.readyReadStandardError.connect(self.read_err)
        self.build_proc.finished.connect(self.finished_building)
        self.build_proc.error.connect(self.error)
        self.cancel_btn.setText('Cancel')
        self.build_btn.setEnabled(False)
        self.build_proc.start(cmd)

    def finished_building(self):
        self.progressbar.setValue(100)
        self.compile_label.setText('Finished Building.')
        self.cancel_btn.setText('Close')

    def error(self, error):
        if error == QtCore.QProcess.FailedToStart:
            msg = "Process failed to start "+ BUILD_CMD
        elif error == QtCore.QProcess.Crashed:
            msg = "Process exit " + BUILD_CMD
        elif error == QtCore.QProcess.Timedout:
            msg = "Process timeout "+ BUILD_CMD
        elif error in (QtCore.QProcess.WriteError, QtCore.QProcess.ReadError):
            msg = "Process communication error " + BUILD_CMD
        else:
            msg = "Unknown error " + BUILD_CMD

        self.compile_label.setText('Process Error')
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
