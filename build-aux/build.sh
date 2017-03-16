#!/bin/bash

ARCH="$(uname 2>/dev/null)"

LinuxInstallation() {

    $PYTHON setup.py install
    cp "mfixgui/menu-linux.json" "$PREFIX/Menu"
    cp "mfixgui/app.svg" "$PREFIX/Menu"
    return $?;
}

DarwinInstallation() {

    # on Mac, we want the homebrew-installed gcc-6 not the Xcode /usr/bin/gcc which is just a shim for clang
    env CC=gcc-6 "${PYTHON}" setup.py install
    cp "mfixgui/menu-osx.json" "$PREFIX/Menu"
    cp "mfixgui/icons/mfix.icns" "$PREFIX/Menu"
    return $?;
}

mkdir -p "$PREFIX/Menu"

case ${ARCH} in
    'Linux')
        LinuxInstallation || exit 1;
        ;;
    'Darwin')
        DarwinInstallation || exit 1;
  ;;
    *)
        echo -e "Unsupported machine type: ${ARCH}";
        exit 1;
        ;;
esac

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
