#!/bin/sh

echo "preparing local settings"

mkdir -p config

aclocal &&
libtoolize &&
autoconf &&
automake --add-missing

echo "run configure next"
echo "for a list of options run configure --help"
echo ""

