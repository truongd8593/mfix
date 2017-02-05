import os
from qtpy import QtCore, QtWidgets, QtGui, PYQT5
from qtpy.QtCore import Qt, QFileSystemWatcher, QSettings, Signal

try:
    from general import get_image_path
except ImportError:
    from mfixgui.tools.general import get_image_path

def create_thumbnail(name='./.thumbnail', model='tfm', geometry=False, chemistry=False, background=None):
    if background is not None and os.path.exists(background):
        img = QtGui.QImage(background)
        base = img.scaled(300, 300).scaled(128, 128, Qt.IgnoreAspectRatio)
    else:
        base = QtGui.QImage(QtCore.QSize(128, 128), QtGui.QImage.Format_ARGB32)
        base.fill(Qt.white)

    painter = QtGui.QPainter(base)

    # add images
    model = QtGui.QImage(get_image_path(model+'.png'))
    painter.drawImage(0, 128-24, model)

    if geometry:
        geo = QtGui.QImage(get_image_path('geometry.png'))
        painter.drawImage(24, 128-24, geo)

    if chemistry:
        geo = QtGui.QImage(get_image_path('chemistry.png'))
        painter.drawImage(24*2, 128-24, geo)

    base.save(name, "PNG")
    del painter

if __name__  == '__main__':
    create_thumbnail('.thumbnail', 'dem', True, True, 'setup.png')