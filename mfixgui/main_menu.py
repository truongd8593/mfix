import os
from qtpy import QtCore, QtWidgets, QtGui, PYQT5

from mfixgui.tools.general import SCRIPT_DIRECTORY
from mfixgui.tools.general import (get_icon, get_mfix_home, widget_iter,
                           is_text_string, is_unicode, get_image_path,
                           format_key_with_args, to_unicode_from_fs)

class MainMenu(object):
    main_menu_animation_speed = 150
    def init_main_menu(self):
        """build the main menu"""
        self.main_menu = QtWidgets.QWidget(self)
        self.main_menu.setObjectName('main_menu')
        self.main_menu.setStyleSheet('QWidget#main_menu{background-color: #E0E0E0;}')
        self.main_menu.hide()

        layout = self.ui.main_menu_layout = QtWidgets.QGridLayout(self.main_menu)
        layout.setContentsMargins(0, 0, 0, 0)

        # return
        btn = self.ui.main_menu_return = QtWidgets.QToolButton(self.main_menu)
        btn.setIcon(get_icon('left.png'))
        btn.setText('Back')
        btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        btn.setAutoRaise(True)
        btn.clicked.connect(self.handle_main_menu_hide)
        layout.addWidget(btn, 0, 0)

        # list widget
        lw = self.ui.main_menu_list = QtWidgets.QListWidget()
        lw.setMaximumWidth(200)
        lw.setFrameStyle(lw.NoFrame)
        lw.setStyleSheet('QListView{background-color: #E0E0E0;}')
        lw.selectionModel().selectionChanged.connect(self.handle_main_menu_selection_changed)
        for name, icon in zip(['Info', 'New', 'Open', 'Save', 'Save As', 'Export', 'Settings', 'Help', 'About', 'Close'],
                              ['infooutline', 'newfolder', 'openfolder', 'save', 'save', 'open_in_new', 'settings', 'help', 'infooutline', 'close']):
            li = QtWidgets.QListWidgetItem(get_icon(icon+'.png'), name)
            lw.addItem(li)
        layout.addWidget(lw, 1, 0)

        # stacked widget
        st = self.ui.main_menu_stackedwidget = QtWidgets.QStackedWidget()
        layout.addWidget(st, 0, 1, 2, 1)

        # blank widget
        bw = self.ui.main_menu_blank_widget = QtWidgets.QWidget()
        bw.setStyleSheet('QWidget{background-color: white;}')
        st.addWidget(bw)

        # build open
        ow = self.ui.main_menu_open_widget = QtWidgets.QWidget()
        ow.setObjectName('browse_stack')
        ow.setStyleSheet('QWidget#browse_stack{background-color: white;}')
        st.addWidget(ow)
        ow_layout = QtWidgets.QGridLayout(ow)

        browse = self.ui.main_menu_browse = QtWidgets.QToolButton()
        browse.setText('Browse')
        browse.setIcon(get_icon('openfolder.png'))
        browse.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        browse.clicked.connect(self.handle_open)
        ow_layout.addWidget(browse, 0, 0)

        lw = self.ui.main_menu_loc_lw = QtWidgets.QListWidget()
        lw.setMaximumWidth(200)
        lw.selectionModel().selectionChanged.connect(self.handle_main_menu_browse_loc_changes)
        loc = ['Recent']
        mfx_dir = os.path.dirname(SCRIPT_DIRECTORY)

        self.tutorial_paths = []
        for root, dirs, files in os.walk(os.path.join(mfx_dir, 'tutorials')):
            if any(f.endswith('.dat') for f in files):
                self.tutorial_paths.append(root)
        if self.tutorial_paths:
            loc += ['Tutorials']

        self.benchmark_paths = []
        for root, dirs, files in os.walk(os.path.join(mfx_dir, 'benchmarks')):
            if any(f.endswith('.dat') for f in files):
                self.benchmark_paths.append(root)
        if self.benchmark_paths:
            loc += ['Benchmarks']

        lw.addItems(loc)
        ow_layout.addWidget(lw, 1, 0)

        lw = self.ui.main_menu_file_lw = QtWidgets.QListWidget()
        lw.itemDoubleClicked.connect(self.handle_main_menu_open_project)
        ow_layout.addWidget(lw, 1, 1)

    def handle_main_menu_open_project(self, item):
        if self.unsaved_flag:
            confirm = self.message(text="Project not saved\nData will be lost!\nProceed?",
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm != 'yes':
                return
            self.clear_unsaved_flag()

        text = str(item.text())
        self.open_project(text)

    def handle_main_menu_browse_loc_changes(self, selected, deselected):
        if selected:
            text = str(self.ui.main_menu_loc_lw.item(selected.indexes()[0].row()).text()).lower()

            self.ui.main_menu_file_lw.clear()

            if text == 'tutorials':
                self.ui.main_menu_file_lw.addItems(self.tutorial_paths)
            elif text == 'benchmarks':
                self.ui.main_menu_file_lw.addItems(self.benchmark_paths)
            elif text == 'recent':
                prjs = self.settings.value('recent_projects')
                if prjs:
                    self.ui.main_menu_file_lw.addItems(prjs.split('|'))


    def handle_main_menu_selection_changed(self, selected, deselected):
        if selected:
            text = str(self.ui.main_menu_list.item(selected.indexes()[0].row()).text()).lower()

            if text == 'new':
                self.handle_main_menu_hide()
                self.new_project()
            elif text == 'save':
                self.handle_main_menu_hide()
                self.handle_save()
            elif text == 'save as':
                self.handle_main_menu_hide()
                self.handle_save_as()
            elif text == 'export':
                self.handle_main_menu_hide()
                self.handle_export()
            elif text == 'close':
                self.close()
            elif text == 'open':
                self.ui.main_menu_stackedwidget.setCurrentIndex(1)
            else:
                self.ui.main_menu_stackedwidget.setCurrentIndex(0)


    def handle_main_menu(self):
        """Show the main menu"""

        if self.get_project_file() is None:
            self.ui.main_menu_list.setCurrentRow(2)

        else:
            self.ui.main_menu_list.setCurrentRow(0)

        self.ui.main_menu_loc_lw.setCurrentRow(0)

        tw, th = self.ui.toolbutton_file.width(), self.ui.toolbutton_file.height()
        self.ui.main_menu_return.setMinimumWidth(tw)
        self.ui.main_menu_return.setMinimumHeight(th)

        # animate
        w, h = self.width(), self.height()
        self.main_menu.setGeometry(-w/2, 0, w, h)
        self.main_menu.show()
        self.main_menu.raise_()
        ani = self.main_menu_animation = self.create_main_menu_animation(self.main_menu, -w/4, 0, 0, 0)
        ani.finished.connect(self.main_menu_animation_finished)
        ani.start()

    def create_main_menu_animation(self, target, x_start, y_start, x_end,y_end):
        animation = QtCore.QPropertyAnimation(target, "pos".encode('utf-8'))
        animation.setDuration(self.main_menu_animation_speed)
        animation.setStartValue(QtCore.QPoint(x_start, y_start))
        animation.setEndValue(QtCore.QPoint(x_end,y_end))

        return animation

    def handle_main_menu_hide(self):
        """Show the main menu"""
        # animate
        w= self.width()
        ani = self.main_menu_animation = self.create_main_menu_animation(self.main_menu, 0, 0, -w/4, 0)
        ani.finished.connect(self.main_menu_animation_finished_hide)
        ani.start()

    def main_menu_animation_finished(self):
        w, h = self.width(), self.height()
        self.main_menu.setGeometry(0, 0, w, h)
        self.main_menu.show()
        self.main_menu.raise_()

    def main_menu_animation_finished_hide(self):
        self.main_menu.hide()