# -*- mode: python -*-

block_cipher = None


a = Analysis(['gui.py'],
             pathex=['/Users/mwm126/gui/gui'],
             binaries=None,
             datas=[('widgets/burcat.pickle','widgets'),
                    ('pymfix','.'),
                    ('../mfix*.so','.'),
                    ('uifiles/regions.ui','uifiles'),
                    ('uifiles/run_popup.ui','uifiles'),
                    ('uifiles/species_popup.ui','uifiles'),
             ],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='gui',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='gui')
