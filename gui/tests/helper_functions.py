from qtpy import QtWidgets, QtTest, PYQT5
import unittest
from datetime import datetime as datetime_, timedelta

# Note, qWaitForWindowShown is deprecated in qt5, which provides a better
# version, including a timeout.
if PYQT5:
    waitForWindow = QtTest.QTest.qWaitForWindowActive
else:
    waitForWindow = QtTest.QTest.qWaitForWindowShown

# holds a global QApplication instance created in the qapp fixture; keeping
# this reference alive avoids it being garbage collected too early
QTAPP = None

def qapp():
    """
    instantiates a global QApplication instance that will be used by the tests.
    """
#    app = QtGui.QApplication.instance()
    global QTAPP
    if QTAPP is None:
        QTAPP = QtWidgets.QApplication([])

    return QTAPP


class TestQApplication(unittest.TestCase):
    '''
    A base class that handles creation and removal of a QApplication instance.
    '''
    def setUp(self):
        '''Creates the QApplication instance'''
        unittest.TestCase.setUp(self)
        self.qapp = qapp()

    def tearDown(self):
        '''Deletes the reference owned by self'''
        if hasattr(self, 'mfix'):
            self.mfix.close()
#       self.qapp.deleteLater() # JMW - causes segfault inbetween gui tests
        unittest.TestCase.tearDown(self)


def waitFor(t):
    end = datetime_.now() + timedelta(milliseconds=t)
    while datetime_.now() < end:
        QtWidgets.QApplication.processEvents()
