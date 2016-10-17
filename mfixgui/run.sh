#!/bin/sh
# Convenience script for developers: run the mfix gui from current directory, without installing

PYTHONPATH="`pwd`/..:$PYTHONPATH" python -m mfixgui.gui "$@"
