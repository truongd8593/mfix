if not exist "%PREFIX%\Menu" mkdir "%PREFIX%\Menu"
copy "mfixgui\menu-windows.json" "%PREFIX%\Menu"
copy "mfixgui\icons\mfix.ico" "%PREFIX%\Menu"

"%PYTHON%" setup.py install
if errorlevel 1 exit 1

:: Add more build steps here, if they are necessary.

:: See
:: http://docs.continuum.io/conda/build.html
:: for a list of environment variables that are set during the build process.
