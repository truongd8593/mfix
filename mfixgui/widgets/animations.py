# -*- coding: utf-8 -*-
import sys
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtCore import Qt
import random

class Base(QtWidgets.QWidget):
    timer_id = -1
    size_hint = (100, 100)
    speed = 80
    animating = False
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

    def start(self, speed=None):
        '''start the animation'''
        if speed is None:
            speed = self.speed
        else:
            self.speed = speed
        if self.timer_id == -1:
            self.timer_id = self.startTimer(speed)
        self.animating = True

    def stop(self):
        '''stop the animation'''
        if self.timer_id != -1:
            self.killTimer(self.timer_id)
        self.timer_id = -1
        self.animating = False
        self.update()

    def timerEvent(self, event):
        '''called when timer times out'''
        self.update()

    def sizeHint(self):
        '''return a size so that layouts can pack it'''
        return QtCore.QSize(*self.size_hint)


class BusyIndicator(Base):
    def __init__(self, parent=None):
        Base.__init__(self, parent)
        self.particles = []
        self.count = 20
        self.size_factor = 10
        self.particle_color = QtGui.QColor('#1974f4')

    def init_particles(self):
        '''initialize the particles'''
        width = self.width()
        height = self.height()
        ps = self.particles = []
        for i in range(self.count):
            m = int(min(width, height)/self.size_factor)
            x = random.randint(0, width)
            y = random.randint(0, height)
            vx = random.randint(-m, m)
            if vx == 0: vx=1
            vy = random.randint(-m, m)
            r = random.randint(m//2, m)
            ps.append([x, y, vx, vy, r])

    def advance(self):
        '''Advance the time'''
        ps = self.particles
        width = self.width()
        height = self.height()
        for p in ps:
            r = p[-1]
            nx = p[0] + p[2]
            if (nx+r > width and p[2] > 0) or (nx-r < 0 and p[2] < 0):
                p[2] = -p[2]
            ny = p[1] + p[3]
            if (ny+r > height and p[3] > 0) or (ny-r < 0 and p[3] < 0):
                p[3] = -p[3]

            p[0] = nx
            p[1] = ny

    def paintEvent (self, event):
        '''Paint the widget'''
        if not self.particles:
            self.init_particles()

        if self.animating:
            self.advance()

        # draw the particles
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)
        painter.setBrush(self.particle_color)

        for x, y, vx, vy, r in self.particles:
            painter.drawEllipse(x-r, y-r, 2*r, 2*r)


class StatusIndicator(Base):
    def __init__(self, parent=None):
        Base.__init__(self, parent)
        self.progress = 0
        self.progressbar_color = QtGui.QColor('#1974f4')
        self.text_color = QtGui.QColor(Qt.black)
        self.text = 'Status'
        self.setMinimumWidth(500)
        self.setMinimumHeight(30)

    def change_text(self, text='test'):
        self.text = text
        self.update()

    def set_progress(self, percent=0.0):
        '''set the progress [0-1], animation changes'''
        self.progress = percent

    def paintEvent (self, event):
        '''Paint the widget'''

        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)


        # background
        rect = self.rect()
        rect.setHeight(self.height()-8)
        rect.moveTop(4)
        rect.setWidth(self.width()-8)
        rect.moveLeft(4)
        gradient = QtGui.QLinearGradient(rect.topLeft(), rect.bottomLeft())
        gradient.setColorAt(0, QtGui.QColor(200,200,200,200))
        gradient.setColorAt(1, QtGui.QColor(200,200,200,255))
        painter.setBrush(gradient)
        painter.setPen(QtGui.QColor(0,0,0,80))
        painter.drawRoundedRect(rect, 4, 4)

        # progress bar
        if self.progress:
            rect.setHeight(rect.height()-4)
            rect.moveTop(6)
            rect.setWidth(rect.width()-4)
            rect.moveLeft(6)
            rect.setWidth(int(rect.width()*self.progress))
            gradient = QtGui.QLinearGradient(rect.topLeft(), rect.bottomLeft())
            gradient.setColorAt(0, self.progressbar_color.lighter())
            gradient.setColorAt(1, self.progressbar_color)
            painter.setBrush(gradient)
            painter.setPen(Qt.NoPen)
            painter.drawRoundedRect(rect, 4, 4)

        # draw the text
        font = painter.font()
        painter.setFont(font)
        painter.setPen(self.text_color)
        painter.drawText(self.rect(), Qt.AlignCenter|Qt.AlignVCenter,self.text)



if __name__  == '__main__':

    # create the QApplication
    qapp = QtWidgets.QApplication([])

#    bi = BusyIndicator()
#    bi.setStyleSheet('''QWidget{
#                    	   background-color: #E0E0E0;
#                        }''')
#    bi.show()
#    bi.start(50)

    si = StatusIndicator()
    si.setGeometry(QtCore.QRect(100, 0, 500, 25))
    si.setStyleSheet('''QWidget{
                    	   background-color: #E0E0E0;
                        }''')
    si.show()
    si.start(50)
    si.set_progress(.5)

    qapp.exec_()
    sys.exit()