import os
import json
from qtpy import QtCore, QtWidgets, QtGui

from mfixgui.tools.general import get_icon, get_mfix_home, get_separator, get_pixmap
from mfixgui.version import __version__
from mfixgui.widgets.workflow import PYQTNODE_AVAILABLE

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
                lw.setItemWidget(li, get_separator(vertical=False))
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

        label_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        label_policy_m = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        spacer = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum,)
        # --- build open ---
        ow = self.ui.main_menu_open_widget = QtWidgets.QWidget()
        ow.setObjectName('open')
        ow.setStyleSheet('QWidget#open{background-color: white;}')
        st.addWidget(ow)
        ow_layout = QtWidgets.QGridLayout(ow)

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

        lw = self.ui.main_menu_loc_lw = QtWidgets.QListWidget()
        lw.setMaximumWidth(200)
        lw.selectionModel().selectionChanged.connect(self.handle_main_menu_browse_loc_changes)
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

        loc = ['Recent']
        lw.addItems(loc)
        ow_layout.addWidget(lw, 2, 0)

        lw = self.ui.main_menu_file_lw = QtWidgets.QListWidget()
        lw.itemDoubleClicked.connect(self.handle_main_menu_open_project)
        ow_layout.addWidget(lw, 2, 1)

        #--- build new ---
        nw = self.ui.main_menu_open_widget = QtWidgets.QWidget()
        nw.setObjectName('new')
        nw.setStyleSheet('QWidget#new{background-color: white;}')
        st.addWidget(nw)
        nw_layout = QtWidgets.QGridLayout(nw)

        l = self.ui.main_menu_label = QtWidgets.QLabel('New')
        l.setStyleSheet('font-size: 24pt;background-color: white;')
        l.setSizePolicy(label_policy)
        nw_layout.addWidget(l, 0, 0, 1, -1)

        fw = QtWidgets.QWidget()
        nw_layout.addWidget(fw, 1, 0, 1, -1)
        fw_layout = QtWidgets.QHBoxLayout(fw)

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
        lw.setViewMode(QtWidgets.QListWidget.IconMode)
        lw.setIconSize(QtCore.QSize(128,128))
        lw.setUniformItemSizes(True)
        lw.setResizeMode(QtWidgets.QListWidget.Adjust)
        lw.itemDoubleClicked.connect(self.handle_main_menu_new_proect)
        nw_layout.addWidget(lw, 2, 0, 1, -1)

        mfx_dir = get_mfix_home()
        for dirname, paths in (('defaults', self.default_paths), ('tutorials', self.tutorial_paths), ('benchmarks',self.benchmark_paths)):
            for path in paths:
                name = os.path.basename(path)
                full_path = os.path.join(mfx_dir, dirname, path)
                thumb_path = os.path.join(full_path, '.thumbnail')
                info_path = os.path.join(full_path, '.mfixguiinfo')
                pix = QtGui.QPixmap()
                img = QtGui.QImage(thumb_path, 'PNG')
                pix.convertFromImage(img)
                item = QtWidgets.QListWidgetItem(QtGui.QIcon(pix), name)
                item.full_path = full_path
                item.enable_list = [False]*7
                info = None
                if os.path.exists(info_path):
                    with open(info_path) as f:
                        info = f.readlines()[0].split(',')
                if info is not None and len(info) == 3:
                    s = [info[0] == t for t in ['single', 'tfm', 'pic', 'dem', 'hybrid']]
                    s.extend([i.lower() == 'true' for i in info[-2:]])
                    item.enable_list = s
                lw.addItem(item)

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

        l = QtWidgets.QLabel('MFiX GUI Version:')
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

        vl = QtWidgets.QLabel('MFiX GUI version: {}'.format(__version__))
        vl.setStyleSheet('background-color: white;')
        aw_layout.addWidget(vl, 2, 0, 1, -1)

        il = QtWidgets.QLabel('''
        MFiX is an open-source multiphase flow solver and is free to download
        and use. MFiX provides a suite of models that treat the carrier phase
        (typically the gas phase) and disperse phase (typically the solids
        phase) differently. These models include the Two Fluid Model (TFM), the
        Particle in Cell (PIC) model, the Discrete Element Model (DEM), and the
        Eulerian-Lagrangian-Eulerian (Hybrid) model.
        ''')
        il.setStyleSheet('background-color: white;')
        il.setWordWrap(True)
        aw_layout.addWidget(il, 3, 0, 1, -1)


        aw_layout.addItem(QtWidgets.QSpacerItem(100, 100, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.MinimumExpanding,), 100, 0)

    def handle_main_menu_open_project(self, item):

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

        project_path = str(item.text())

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

    def handle_main_menu_browse_loc_changes(self, selected, deselected):
        if selected:
            self.set_file_listwidget(self.ui.main_menu_loc_lw.item(selected.indexes()[0].row()).text().lower())

    def set_file_listwidget(self, text):

        self.ui.main_menu_file_lw.clear()

        if text == 'tutorials':
            self.ui.main_menu_file_lw.addItems(self.tutorial_paths)
        elif text == 'benchmarks':
            self.ui.main_menu_file_lw.addItems(self.benchmark_paths)
        elif text == 'defaults':
            self.ui.main_menu_file_lw.addItems(self.default_paths)
        elif text == 'recent':
            projects = self.settings.value('recent_projects')
            if projects:
                existing_projects = [ project for project in projects.split('|') if os.path.exists(project)]
                self.ui.main_menu_file_lw.addItems(existing_projects)
        elif text == 'clear recent':
            self.settings.setValue('recent_projects', '|'.join([]))


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


    def handle_main_menu(self):
        """Show the main menu"""

        project_file = self.get_project_file()
        if project_file is None:
            self.ui.main_menu_list.setCurrentRow(2)
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
            hide = any( e==False and i==True for e, i in zip(self.main_menu_new_enable_list, item.enable_list))
            item.setHidden(hide)