#!/bin/sh

echo "preparing local settings"

mkdir -p config

aclocal &&
libtoolize  --copy &&
autoconf &&
autoheader &&
automake --add-missing --copy

echo "run configure next"
echo "for a list of options run configure --help"
echo ""

