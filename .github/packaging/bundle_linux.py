#!/usr/bin/env python3
"""Make an openEMS package directory relocatable on Linux.

usage: bundle_linux.py <package dir>

The package directory holds bin/ (the programs) and lib/ (the openEMS, CSXCAD,
fparser and nf2ff libraries), copied from an installation. The shared libraries
they load are copied into lib/deps, except those every system has and that must
match it: glibc, libstdc++, libgcc_s and the NVIDIA driver (libcuda, loaded at
run time anyway). The run paths become relative (patchelf):

- bin/:      $ORIGIN/../lib:$ORIGIN/../lib/deps
- lib/:      $ORIGIN:$ORIGIN/deps
- lib/deps/: $ORIGIN

lib/ then holds only our own libraries, so it can be put in LD_LIBRARY_PATH for
the Python modules without shadowing the system's libraries.
"""

import os
import re
import shutil
import subprocess
import sys

# provided by the system: the C and C++ runtimes and the NVIDIA driver
SYSTEM = re.compile(r'^(linux-vdso|ld-linux|libc|libm|libpthread|libdl|librt|libutil|libresolv|libmvec|'
                    r'libstdc\+\+|libgcc_s|libcuda|libnvidia-)[.-]')


def run(*args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def is_elf(path):
    if not os.path.isfile(path) or os.path.islink(path):
        return False
    with open(path, 'rb') as f:
        return f.read(4) == b'\x7fELF'


def main():
    pkg = os.path.abspath(sys.argv[1])
    bin_dir, lib_dir = os.path.join(pkg, 'bin'), os.path.join(pkg, 'lib')
    deps_dir = os.path.join(lib_dir, 'deps')
    os.makedirs(deps_dir, exist_ok=True)
    ours = {n for n in os.listdir(lib_dir) if '.so' in n}
    files = [os.path.join(bin_dir, n) for n in sorted(os.listdir(bin_dir))]
    files += [os.path.join(lib_dir, n) for n in sorted(ours)]
    files = [f for f in files if is_elf(f)]

    # ldd lists the whole closure: resolved with the installed libraries of the build machine
    env = dict(os.environ, LD_LIBRARY_PATH=lib_dir)
    copied = {}
    for f in files:
        out = subprocess.run(['ldd', f], capture_output=True, text=True, env=env).stdout
        for line in out.splitlines():
            m = re.match(r'\s*(\S+) => (\S+)', line)
            if not m:
                if 'not found' in line:
                    sys.exit(f'bundle_linux: {line.strip()} ({f})')
                continue
            name, path = m.groups()
            if name in ours or SYSTEM.match(name) or name in copied:
                continue
            if path == 'not':
                sys.exit(f'bundle_linux: {name} not found ({f})')
            shutil.copy2(os.path.realpath(path), os.path.join(deps_dir, name))
            os.chmod(os.path.join(deps_dir, name), 0o755)
            copied[name] = path

    for f in files:
        rpath = '$ORIGIN/../lib:$ORIGIN/../lib/deps' if f.startswith(bin_dir + os.sep) else '$ORIGIN:$ORIGIN/deps'
        run('patchelf', '--set-rpath', rpath, f)
    for name in copied:
        run('patchelf', '--set-rpath', '$ORIGIN', os.path.join(deps_dir, name))
    for name in sorted(copied):
        print(f'bundle_linux: lib/deps/{name}')
    print(f'bundle_linux: {len(copied)} libraries in lib/deps')


if __name__ == '__main__':
    main()
