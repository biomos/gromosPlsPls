#!/bin/sh

echo "preparing local settings"

mkdir -p config

aclocal --force &&
libtoolize  --copy --force &&
autoconf --force &&
autoheader --force &&
automake --add-missing --copy --force
svn revert INSTALL

echo "run configure next"
echo "for a list of options run configure --help"
echo ""

