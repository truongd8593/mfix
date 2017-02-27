"""the main (file menu) bar."""
import json
import os
import subprocess
import sys

import qtpy
from mfixgui.tools.general import (SCRIPT_DIRECTORY, get_icon, get_mfix_home,
                                   get_pixmap, get_separator)
from mfixgui.version import __version__, get_git_revision_short_hash
from mfixgui.widgets.workflow import PYQTNODE_AVAILABLE
from qtpy import API_NAME, QtCore, QtGui, QtWidgets

try:
    from qtpy import QT_VERSION
except ImportError:
    QT_VERSION = 'Unknown'


try:
    import numpy as np
    numpy_version = np.__version__
except ImportError:
    numpy_version = 'Import Failed'

try:
    from vtk import vtkVersion
    vtkVersion = vtkVersion.GetVTKVersion()
except ImportError:
    vtkVersion = 'Import Failed'


# FIXME: should we use six.moves.urllib_parse.urljoin six.moves.urllib.request.pathname2url instead?
try:
    # Python 3
    import urllib.request as urlparse
    import urllib.request as urllib
except ImportError:
    # Python 2
    import urlparse
    import urllib


def path2url(path):
    """Convert path to url."""
    return urlparse.urljoin(
      'file:', urllib.pathname2url(path))

class MainMenu(object):
    """Main menu mixin for the gui."""

    main_menu_animation_speed = 150
    def init_main_menu(self):
        """Build the main menu."""

        self.default_paths = None
        self.tutorial_paths = None
        self.benchmark_paths = None

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
        lw.setStyleSheet('''QListView{background-color: #E0E0E0;}
                            QListView::item:selected{background: #64B5F6;
                                color: white;}
                            QListView::item:hover{background:#BBDEFB;}
                         ''')
        lw.selectionModel().selectionChanged.connect(self.handle_main_menu_selection_changed)

        names = ['Project Info', 'New', 'Open', 'Save', 'Save As', 'Export Project', 'sep', 'Settings', 'Help', 'About', 'Quit']
        icons = ['infooutline', 'newfolder', 'openfolder', 'save', 'save', 'open_in_new', '', 'settings', 'help', 'infooutline', 'close']
        if PYQTNODE_AVAILABLE:
            for n, i in reversed([('Export Workflow', 'open_in_new'), ('Import Workflow', 'import'),
                         ('sep', '')]):
                names.insert(7, n)
                icons.insert(7, i)

        for name, icon in zip(names, icons):
            if name == 'sep':
                li = QtWidgets.QListWidgetItem()
                li.setFlags(QtCore.Qt.NoItemFlags)
                li.setSizeHint(QtCore.QSize(0, 10))
                lw.addItem(li)
                sep = get_separator(vertical=False)
                sep.setEnabled(False)
                lw.setItemWidget(li, sep)
            else:
                li = QtWidgets.QListWidgetItem(get_icon(icon+'.png'), name)
                lw.addItem(li)
        layout.addWidget(lw, 1, 0)

        # logo
        l = QtWidgets.QLabel()
        pixmap = get_pixmap('mfix.png', 84, 84)
        l.setPixmap(pixmap)
        l.setAlignment(QtCore.Qt.AlignCenter)
        layout.addWidget(l, 2, 0)

        # stacked widget
        st = self.ui.main_menu_stackedwidget = QtWidgets.QStackedWidget()
        layout.addWidget(st, 0, 1, -1, 1)

        # blank widget
        bw = self.ui.main_menu_blank_widget = QtWidgets.QWidget()
        bw.setObjectName('default')
        bw.setStyleSheet('QWidget{background-color: white;}')
        st.addWidget(bw)

        # size policies
        label_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        label_policy_m = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        spacer = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum,)

        # --- build open ---
        ow = self.ui.main_menu_open_widget = QtWidgets.QWidget()
        ow.setObjectName('open')
        ow.setStyleSheet('QWidget#open{background-color: white;}')
        st.addWidget(ow)
        ow_layout = QtWidgets.QGridLayout(ow)
        ow_layout.setContentsMargins(9, 9, 0, 0)

        l = self.ui.main_menu_label = QtWidgets.QLabel('Open')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        ow_layout.addWidget(l, 0, 0, 1, -1)

        browse = self.ui.main_menu_browse = QtWidgets.QToolButton()
        browse.setText('Browse')
        browse.setIcon(get_icon('openfolder.png'))
        browse.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        browse.clicked.connect(self.handle_open)
        ow_layout.addWidget(browse, 1, 0)

        spacer_exp = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum,)
        ow_layout.addItem(spacer_exp, 1, 1)

        lw_f = self.ui.main_menu_file_lw = QtWidgets.QListWidget()
        lw_f.setFrameStyle(lw_f.NoFrame)
        lw_f.setIconSize(QtCore.QSize(128,128))
        lw_f.setUniformItemSizes(True)
        lw_f.setResizeMode(QtWidgets.QListWidget.Adjust)
        lw_f.itemDoubleClicked.connect(self.handle_main_menu_open_project)
        ow_layout.addWidget(lw_f, 2, 0, 1, -1)

        tb = QtWidgets.QToolButton()
        tb.setText('Clear Recent')
        tb.setToolTip('Clear list of recent projects')
        tb.setAutoRaise(True)
        tb.pressed.connect(self.handle_clear_recent)
        ow_layout.addWidget(tb, 1, 2)

        tb_list = QtWidgets.QToolButton()
        tb_list.setIcon(get_icon('list.png'))
        tb_list.setToolTip('View projects as list')
        tb_list.setAutoRaise(True)
        tb_list.setCheckable(True)
        ow_layout.addWidget(tb_list, 1, 3)

        tb_tile = QtWidgets.QToolButton()
        tb_tile.setIcon(get_icon('tile.png'))
        tb_tile.setToolTip('View projects as grid')
        tb_tile.setAutoRaise(True)
        tb_tile.setCheckable(True)
        ow_layout.addWidget(tb_tile, 1, 4)

        # apply previous state
        if self.settings.value('open_list_mode', 'icon') == 'icon':
            lw_f.setViewMode(QtWidgets.QListWidget.IconMode)
            tb_tile.setChecked(True)
        else:
            tb_list.setChecked(True)

        def callback_tile():
            lw_f.setViewMode(QtWidgets.QListWidget.IconMode)
            self.settings.setValue('open_list_mode', 'icon')
            tb_list.setChecked(False)
            if not tb_tile.isChecked():
                tb_tile.setChecked(True)
        tb_tile.clicked.connect(callback_tile)

        def callback_list():
            lw_f.setViewMode(QtWidgets.QListWidget.ListMode)
            self.settings.setValue('open_list_mode', 'list')
            tb_tile.setChecked(False)
            if not tb_list.isChecked():
                tb_list.setChecked(True)
        tb_list.clicked.connect(callback_list)

        #--- build new ---
        nw = self.ui.main_menu_open_widget = QtWidgets.QWidget()
        nw.setObjectName('new')
        nw.setStyleSheet('QWidget#new{background-color: white;}')
        st.addWidget(nw)
        nw_layout = QtWidgets.QGridLayout(nw)
        nw_layout.setContentsMargins(9, 9, 0, 0)

        l = self.ui.main_menu_label = QtWidgets.QLabel('New')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        nw_layout.addWidget(l, 0, 0, 1, -1)

        fw = QtWidgets.QWidget()
        nw_layout.addWidget(fw, 1, 0, 1, -1)
        fw_layout = QtWidgets.QHBoxLayout(fw)
        fw_layout.setContentsMargins(0, 0, 0, 0)

        self.main_menu_new_enable_list = [True]*7
        for i, f in enumerate(['Single', 'TFM', 'PIC', 'DEM', 'Hybrid', 'Cartesian', 'Chemistry']):
            cb = QtWidgets.QCheckBox()
            cb.setIcon(get_icon('geometry.png' if f == 'Cartesian' else f.lower()+'.png'))
            cb.setChecked(True)
            cb.setToolTip(f)
            cb.toggled.connect(lambda checked, idx=i: self.main_menu_filter_new(idx, checked))
            fw_layout.addWidget(cb)
        hspacer = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Maximum,)
        fw_layout.addSpacerItem(hspacer)

        lw = self.ui.main_menu_new_list = QtWidgets.QListWidget()
        lw.setFrameStyle(lw.NoFrame)
        lw.setIconSize(QtCore.QSize(128, 128))
        lw.setUniformItemSizes(True)
        lw.setResizeMode(QtWidgets.QListWidget.Adjust)
        lw.itemDoubleClicked.connect(self.handle_main_menu_new_proect)
        nw_layout.addWidget(lw, 2, 0, 1, -1)

        tb_l = QtWidgets.QToolButton()
        tb_l.setIcon(get_icon('list.png'))
        tb_l.setToolTip('View projects as list')
        tb_l.setAutoRaise(True)
        tb_l.setCheckable(True)
        fw_layout.addWidget(tb_l)

        tb_t = QtWidgets.QToolButton()
        tb_t.setIcon(get_icon('tile.png'))
        tb_t.setToolTip('View projects as grid')
        tb_t.setAutoRaise(True)
        tb_t.setCheckable(True)
        fw_layout.addWidget(tb_t)

        # apply previous state
        if self.settings.value('new_list_mode', 'icon') == 'icon':
            lw.setViewMode(QtWidgets.QListWidget.IconMode)
            tb_t.setChecked(True)
        else:
            tb_l.setChecked(True)

        def callback():
            lw.setViewMode(QtWidgets.QListWidget.IconMode)
            self.settings.setValue('new_list_mode', 'icon')
            tb_l.setChecked(False)
            if not tb_t.isChecked():
                tb_t.setChecked(True)
        tb_t.clicked.connect(callback)

        def callback2():
            lw.setViewMode(QtWidgets.QListWidget.ListMode)
            self.settings.setValue('new_list_mode', 'list')
            tb_t.setChecked(False)
            if not tb_l.isChecked():
                tb_l.setChecked(True)
        tb_l.clicked.connect(callback2)

        # --- build info ---
        iw = self.ui.main_menu_info_widget = QtWidgets.QWidget()
        iw.setObjectName('project_info')
        iw.setStyleSheet('QWidget#project_info{background-color: white;}')
        st.addWidget(iw)
        iw_layout = QtWidgets.QGridLayout(iw)

        l = QtWidgets.QLabel('Project Info')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 0, 0, 1, -1)

        iw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 1, 0)

        l = QtWidgets.QLabel('Project Version:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 2, 0)

        l = self.ui.main_menu_project_version = QtWidgets.QLabel('Unknown')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 2, 1)

        l = QtWidgets.QLabel('Created with MFiX GUI Version:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 3, 0)

        l = self.ui.main_menu_gui_version = QtWidgets.QLabel('Unknown')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 3, 1)

        iw_layout.addItem(spacer, 4, 0)

        l = QtWidgets.QLabel('Author:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 5, 0)

        l = self.ui.main_menu_created_by = QtWidgets.QLabel('Unknown')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 5, 1)

        l = QtWidgets.QLabel('Modified By:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 6, 0)

        l = self.ui.main_menu_modified_by = QtWidgets.QLabel('Unknown')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 6, 1)

        iw_layout.addItem(spacer, 7, 0)

        l = QtWidgets.QLabel('Last Modified:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 8, 0)

        l = self.ui.main_menu_modified_time = QtWidgets.QLabel('')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 8, 1)

        l = QtWidgets.QLabel('Created:')
        l.setSizePolicy(label_policy_m)
        iw_layout.addWidget(l, 9, 0)

        l = self.ui.main_menu_created_time = QtWidgets.QLabel('')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 9, 1)

        iw_layout.addItem(spacer, 10, 0)

        l = QtWidgets.QLabel('Notes:')
        l.setSizePolicy(label_policy)
        iw_layout.addWidget(l, 30, 0, 1, -1)

        nt = self.ui.main_menu_project_notes = QtWidgets.QTextEdit()
        iw_layout.addWidget(nt, 31, 0, 1, -1)

        iw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 100, 0)

        # --- build settings ---
        spacer = QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum,)

        sw = self.ui.main_menu_settings_widget = QtWidgets.QWidget()
        sw.setObjectName('settings')
        sw.setStyleSheet('QWidget#settings{background-color: white;}')
        l.setSizePolicy(label_policy)
        st.addWidget(sw)
        sw_layout = QtWidgets.QGridLayout(sw)

        l = QtWidgets.QLabel('Settings')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        sw_layout.addWidget(l, 0, 0, 1, -1)

        sw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 1, 0)

        # style
        l = QtWidgets.QLabel('Style:')
        l.setSizePolicy(label_policy_m)
        sw_layout.addWidget(l, 2, 0)

        sc = QtWidgets.QComboBox()
        sc.addItems([s.lower() for s in QtWidgets.QStyleFactory.keys()])
        cur_style = self.settings.value('app_style')
        sc.setCurrentIndex(sc.findText(cur_style))
        sc.currentIndexChanged.connect(lambda: self.change_app_style(sc.currentText()))
        sw_layout.addWidget(sc, 2, 1)

        sw_layout.addItem(spacer, 3, 0)

        # animation settings
        gb = QtWidgets.QGroupBox()
        gb.setCheckable(True)
        gb.setChecked(self.settings.value('animate', True))
        gb.setTitle('Enable Animations')
        sw_layout.addWidget(gb, 4, 0, 1, 2)
        gb_layout = QtWidgets.QGridLayout(gb)

        al = QtWidgets.QLabel('Animation Speed')
        gb_layout.addWidget(al, 0, 0)
        gb_layout.setContentsMargins(5, 5, 5, 5)
        asb = QtWidgets.QSpinBox()
        asb.valueChanged.connect(lambda n: self.settings.setValue('animation_speed', n))
        asb.setRange(0, 1000)
        asb.setValue(int(self.settings.value('animation_speed', 400)))
        gb_layout.addWidget(asb, 0, 1)

        au = QtWidgets.QLabel('milli seconds')
        gb_layout.addWidget(au, 0, 2)

        sw_layout.addItem(spacer, 5, 0)

        # developer mode
        cb = QtWidgets.QCheckBox()
        cb.setText('Enable developer tools.')
        cb.setChecked(int(self.settings.value('developer_mode', 0)))
        cb.stateChanged.connect(self.enable_developer_mode)
        sw_layout.addWidget(cb, 6, 0, 1, 2)

        sw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 100, 0)

        # --- build help ---
        hw = self.ui.main_menu_help_widget = QtWidgets.QWidget()
        hw.setObjectName('help')
        hw.setStyleSheet('QWidget#help{background-color: white;}')
        st.addWidget(hw)
        hw_layout = QtWidgets.QGridLayout(hw)

        l = QtWidgets.QLabel('Help')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        hw_layout.addWidget(l, 0, 0, 1, -1)

        hw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 1, 0)

        # link only works after running: python setup.py build_doc
        def open_user_guide(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))

        for i, help_text in enumerate([
            'See <a href="%s">User Guide</a> for documentation on using the GUI' % path2url(os.path.join(SCRIPT_DIRECTORY, 'doc', 'USER_GUIDE.html')),
            'See <a href="%s">Setup Guide</a> for documentation on building custom mfixsolvers' % path2url(os.path.join(SCRIPT_DIRECTORY, 'doc', 'SETUP_GUIDE.html')),
            'See <a href="%s">Tutorials</a> text based model setup tutorials' % path2url(os.path.join(SCRIPT_DIRECTORY, 'doc', 'TUTORIALS.html')),
            ]):

            help_info = QtWidgets.QLabel(help_text)
            help_info.setStyleSheet('background-color: white;')
            help_info.setWordWrap(True)
            help_info.linkActivated.connect(open_user_guide)
            help_info.setOpenExternalLinks(True)
            hw_layout.addWidget(help_info, 10+i, 0, 1, -1)

        hw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 100, 0)

        # --- build about ---
        aw = self.ui.main_menu_about_widget = QtWidgets.QWidget()
        aw.setObjectName('about')
        aw.setStyleSheet('QWidget#about{background-color: white;}')
        st.addWidget(aw)
        aw_layout = QtWidgets.QGridLayout(aw)

        l = QtWidgets.QLabel('About')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        aw_layout.addWidget(l, 0, 0, 1, -1)

        aw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 1, 0)

        il = QtWidgets.QLabel(' '.join(l.strip() for l in
        '''
        MFiX is an open-source multiphase flow solver and is free to download
        and use. MFiX provides a suite of models that treat the carrier phase
        (typically the gas phase) and disperse phase (typically the solids
        phase) differently. These models include the Two Fluid Model (TFM), the
        Particle in Cell (PIC) model, the Discrete Element Model (DEM), and the
        Eulerian-Lagrangian-Eulerian (Hybrid) model.
        '''.strip().split('\n')))
        il.setStyleSheet('background-color: white;')
        il.setWordWrap(True)
        aw_layout.addWidget(il, 2, 0, 1, -1)

        aw_layout.addItem(QtWidgets.QSpacerItem(100, 20, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Minimum,), 3, 0)

        git_des = None
        if int(self.settings.value('developer_mode', 0)):
            git_des = get_git_revision_short_hash()

        self.version_labels = [
            '<b>MFiX GUI version:</b> {}'.format(__version__),
            '<b>Git description:</b> {}'.format(git_des) if git_des is not None else None,
            '<b>Python version:</b> {}'.format(sys.version),
            '<b>Qt Wrapper:</b> {}'.format(API_NAME),
            '<b>Qt Version:</b> {}'.format(QT_VERSION),
            '<b>qtpy Version:</b> {}'.format(qtpy.__version__),
            '<b>Numpy Version:</b> {}'.format(numpy_version),
            '<b>VTK Version:</b> {}'.format(vtkVersion),
            ]

        for i, label in enumerate(self.version_labels):
            if label is None:
                continue
            ql = QtWidgets.QLabel(label)
            ql.setStyleSheet('background-color: white;')
            ql.setWordWrap(True)
            aw_layout.addWidget(ql, 10+i, 0, 1, -1)

        # copy to clipboard btn
        def copy_version_info():
            cp = self.app.clipboard()
            cp.setText('\n'.join(str(v).replace('<b>', '').replace('</b>', '') for v in self.version_labels))
        copy_btn = QtWidgets.QToolButton()
        copy_btn.setText('Copy Version Info')
        copy_btn.clicked.connect(copy_version_info)
        aw_layout.addWidget(copy_btn, 100, 0, 1, 2)

        aw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 1000, 0)

    def handle_main_menu_open_project(self, item):
        """Open the project of the selected item"""
        item = self.ui.main_menu_file_lw.currentItem()
        if not item:
            return

        if self.unsaved_flag:
            confirm = self.message(text="Project not saved\nData will be lost!\nProceed?",
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm != 'yes':
                return
            self.clear_unsaved_flag()

        project_path = item.full_path
        if os.path.exists(project_path):
            self.open_project(project_path)
        else:
            self.message(text="File does not exist: %s" % project_path)

    def handle_main_menu_new_proect(self, item):
        if self.unsaved_flag:
            confirm = self.message(text="Project not saved\nData will be lost!\nProceed?",
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm != 'yes':
                return
            self.clear_unsaved_flag()

        self.open_new_from_template(os.path.join(item.full_path, 'mfix.dat'))

    def handle_clear_recent(self):
        self.settings.setValue('recent_projects', '|'.join([]))
        self.ui.main_menu_file_lw.clear()

    def handle_main_menu_selection_changed(self, selected, deselected):
        if selected and selected.indexes():
            text = str(self.ui.main_menu_list.item(selected.indexes()[0].row()).text()).lower()
            sw = self.ui.main_menu_stackedwidget

            if text == 'save':
                self.handle_save()
            elif text == 'save as':
                self.handle_save_as()
            elif text == 'export project':
                self.handle_export()
            elif text == 'export workflow':
                self.ui.workflow_widget.handle_export()
            elif text == 'import workflow':
                self.ui.workflow_widget.handle_import()
            elif text == 'quit':
                self.close()
            else:
                index = 0

                matches = [i
                   for i in range(sw.count())
                   if sw.widget(i).objectName().replace('_', ' ') == text]

                if len(matches) == 1:
                    index = matches[0]
                sw.setCurrentIndex(index)

                if text == 'new' and self.default_paths is None:
                    self.collect_template_files()
                elif text == 'open':
                    self.populate_recent_projects()

    def handle_main_menu(self):
        """Show the main menu"""

        project_file = self.get_project_file()
        if project_file is None:
            if self.settings.value('recent_projects', ''):
                self.ui.main_menu_list.setCurrentRow(2)
            else:
                self.ui.main_menu_list.setCurrentRow(1)
            self.disable_main_menu_items(['project info', 'save', 'save as', 'export project', 'export workflow', 'import workflow'])
        else:
            self.ui.main_menu_modified_time.setText(
                self.project.modified_time
                )
            self.ui.main_menu_created_time.setText(
                self.project.mfix_gui_comments.get('created_date', 'Unknown'))
            self.ui.main_menu_project_version.setText(
                str(self.project.mfix_gui_comments.get('project_version', 'Unknown')))
            self.ui.main_menu_gui_version.setText(
                str(self.project.mfix_gui_comments.get('gui_version', 'Unknown')))
            self.ui.main_menu_created_by.setText(
                self.project.mfix_gui_comments.get('author', 'Unknown'))
            self.ui.main_menu_modified_by.setText(
                self.project.mfix_gui_comments.get('modified_by', '').replace('|', ', '))
            self.ui.main_menu_project_notes.setText(
                json.loads(self.project.mfix_gui_comments.get('project_notes', '""')))
            self.disable_main_menu_items([], ['save'])
            self.ui.main_menu_list.setCurrentRow(0)

        tw, th = self.ui.toolbutton_file.width(), self.ui.toolbutton_file.height()
        self.ui.main_menu_return.setMinimumWidth(tw)
        self.ui.main_menu_return.setMinimumHeight(th)

        # hide the current widget, issue #291, VTK widget overlay issues
        self.ui.tabWidgetGraphics.currentWidget().hide()

        # animate
        w, h = self.width(), self.height()
        self.main_menu.setGeometry(-w/2, 0, w, h)
        self.main_menu.show()
        self.main_menu.raise_()
        ani = self.main_menu_animation = self.create_main_menu_animation(self.main_menu, -w/4, 0, 0, 0)
        ani.finished.connect(self.main_menu_animation_finished)
        ani.start()

    def populate_recent_projects(self):
        projs = self.settings.value('recent_projects', '').split('|')

        if not projs:
            return

        lw = self.ui.main_menu_file_lw
        lw.clear()
        for proj in projs:
            if not os.path.exists(proj):
                continue
            name = os.path.basename(proj)
            dir_ = os.path.dirname(proj)
            thumb_path = os.path.join(dir_, '.thumbnail')

            if os.path.exists(thumb_path):
                pix = QtGui.QPixmap()
                img = QtGui.QImage(thumb_path, 'PNG')
                pix.convertFromImage(img)
                icon = QtGui.QIcon(pix)
            else:
                icon = get_icon('missing_thumbnail.png')

            # look for description
            description = ''
            with open(proj) as f:
                for line in f:
                    clean = line.lower().split('#')[0].split('!')[0].strip()
                    if 'description' in clean:
                        toks = [tok.strip() for tok in line.split('=')]
                        des_ind = toks.index('description')
                        description = toks[des_ind + 1].replace('"', '').replace("'", '')
                        break

            name = name.split('.')[0]
            text = '\n'.join([name, description, proj])
            item = QtWidgets.QListWidgetItem(icon, text)
            item.full_path = proj
            lw.addItem(item)

        sb = lw.verticalScrollBar()
        sb.setValue(0)

    def create_main_menu_animation(self, target, x_start, y_start, x_end, y_end):
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
        """callback when the show animation is finished"""
        w, h = self.width(), self.height()
        self.main_menu.setGeometry(0, 0, w, h)
        self.main_menu.show()
        self.main_menu.raise_()

    def main_menu_animation_finished_hide(self):
        """callback when the hide animation is finished"""
        # show the current widget, issue #291, VTK widget overlay issues
        self.ui.tabWidgetGraphics.currentWidget().show()
        self.main_menu.hide()

    def disable_main_menu_items(self, items, except_items=[]):
        for r in range(self.ui.main_menu_list.count()):
            i = self.ui.main_menu_list.item(r)
            text = str(i.text()).lower()
            if text in except_items:
                continue
            i.setFlags(i.flags() & ~QtCore.Qt.ItemIsEnabled if text in items else i.flags() | QtCore.Qt.ItemIsEnabled )

    def change_app_style(self, style):
        self.app.setStyle(style)
        self.settings.setValue('app_style', style)

    def enable_developer_mode(self, enable):
        self.change_mode('modeler')
        self.ui.pushButtonDeveloper.setVisible(enable)
        self.ui.pushButtonInterpreter.setVisible(enable)
        if enable:
            self.ui.tabWidgetGraphics.addTab(self.ui.api_response, 'api response')
        else:
            self.ui.tabWidgetGraphics.removeTab(self.ui.tabWidgetGraphics.indexOf(self.ui.api_response))

        self.settings.setValue('developer_mode', int(enable))

    def main_menu_filter_new(self, idx, checked):
        self.main_menu_new_enable_list[idx] = checked
        for r in range(self.ui.main_menu_new_list.count()):
            item = self.ui.main_menu_new_list.item(r)
            show = any(e==True and i==True for e, i in zip(self.main_menu_new_enable_list, item.enable_list))
            item.setHidden(not show)

    def collect_template_files(self):

        # --- look for template files ---
        mfx_dir = get_mfix_home()
        def set_paths(dirname, ui_name):
            path_var = []
            top_path = os.path.join(mfx_dir, dirname, '')
            for root, dirs, files in os.walk(top_path):
                if any(f.endswith('mfix.dat') for f in files):
                    path_var.append(root.replace(top_path,''))
            if path_var:
                path_var.sort(key=lambda y: y.lower())
            return path_var

        self.default_paths = set_paths('defaults', 'Defaults')
        self.tutorial_paths = set_paths('tutorials', 'Tutorials')
        self.benchmark_paths = set_paths('benchmarks', 'Benchmarks')

        # read info about the template files
        temp_info = {}
        temp_info_file = os.path.join(SCRIPT_DIRECTORY, 'tools', 'template_data.json')
        if os.path.exists(temp_info_file):
            with open(temp_info_file) as f:
                temp_info = json.load(f)

        lw = self.ui.main_menu_new_list
        # loop through template files, building listwidgetitems
        for dirname, paths in (('defaults', self.default_paths), ('tutorials', self.tutorial_paths), ('benchmarks',self.benchmark_paths)):
            for path in paths:
                name = os.path.basename(path)
                full_path = os.path.join(mfx_dir, dirname, path)
                thumb_path = os.path.join(full_path, '.thumbnail')

                if os.path.exists(thumb_path):
                    pix = QtGui.QPixmap()
                    img = QtGui.QImage(thumb_path, 'PNG')
                    pix.convertFromImage(img)
                    icon = QtGui.QIcon(pix)
                else:
                    icon = get_icon('missing_thumbnail.png')

                # extract info from the template info file
                info = temp_info.get(name, {})
                description = info.get('description', None)
                solver = info.get('solver', 'single')
                geo = info.get('geometry', 'false')
                chem = info.get('chemistry', 'false')

                # build a list of booleans for filtering templates
                # [single, tfm, pic, dem, hybrid, geometry, chemistry]
                enable_list = [solver == t for t in ['single', 'tfm', 'pic', 'dem', 'hybrid']]
                enable_list.extend([i.lower() == 'true' for i in [geo, chem]])

                # if there is a desciption, add it to the name
                if description is not None:
                    name = '\n'.join([name, description])

                item = QtWidgets.QListWidgetItem(icon, name)
                item.full_path = full_path
                item.enable_list = enable_list

                lw.addItem(item)
