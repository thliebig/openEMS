#!/bin/sh

find . -type f -name \*.cpp -exec astyle --style=allman --indent=tab --indent-cases {} \;
find . -type f -name \*.h -exec astyle --style=allman --indent=tab --keep-one-line-blocks --indent-cases {} \;
