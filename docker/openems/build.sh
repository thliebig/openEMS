#!/bin/bash
# usage: build.sh [dev|runtime|all] [openEMS-Project dir]  -- build the images (default: all)
# The openEMS sources: the given directory, else an openEMS-Project tree next to this
# directory (../openEMS-Project, or the tree this directory is in), else the GPU branch on GitHub.
set -e
cd "$(dirname "$0")"
TARGET=${1:-all}
src_ok() { [ -d "$1/fparser" ] && [ -d "$1/CSXCAD" ] && [ -d "$1/openEMS" ] && [ -d "$1/.git/modules" ]; }
SRC=
for d in "$2" ../openEMS-Project ../../..; do
	if [ -n "$d" ] && src_ok "$d"; then SRC=$(cd "$d" && pwd); break; fi
done
if [ -n "$SRC" ]; then
	echo "build.sh: openEMS sources from $SRC"
	SOURCE=(--build-arg OPENEMS_SOURCE=local --build-context openems-src="$SRC")
else
	echo "build.sh: openEMS sources from GitHub"
	SOURCE=(--build-arg OPENEMS_SOURCE=github)
fi
if [ "$TARGET" = dev ] || [ "$TARGET" = all ]; then
	docker build --target dev -t seanmollet/openems-dev:latest .
fi
if [ "$TARGET" = runtime ] || [ "$TARGET" = all ]; then
	docker build --target runtime "${SOURCE[@]}" -t seanmollet/openems:latest .
fi
