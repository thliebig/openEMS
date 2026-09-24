#!/usr/bin/env python3
"""Make an openEMS package directory relocatable on macOS.

usage: bundle_macos.py <package dir>

The package directory holds bin/ (the programs) and lib/ (the openEMS, CSXCAD,
fparser and nf2ff libraries), copied from an installation. The libraries they
load from anywhere but the system (/usr/lib, /System), e.g. Homebrew, are copied
into lib/deps. All references are rewritten to be relative:

- our libraries:   @rpath/<name>, rpath @executable_path/../lib (programs) or
                   @loader_path (libraries)
- other libraries: @executable_path/../lib/deps/<name> (programs),
                   @loader_path/deps/<name> (lib/), @loader_path/<name> (lib/deps/)

The old rpaths are removed and every changed file is signed ad hoc again
(required on Apple silicon). lib/ then holds only our own libraries, so it can
be put in DYLD_LIBRARY_PATH for the Python modules without shadowing others.
"""

import os
import subprocess
import sys


def run(*args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def is_macho(path):
    return os.path.isfile(path) and not os.path.islink(path) and run('file', '-b', path).startswith('Mach-O')


def install_id(path):
    lines = run('otool', '-D', path).splitlines()
    return lines[1].strip() if len(lines) > 1 else None


def references(path):
    """the libraries a Mach-O file loads (without its own id)"""
    own = install_id(path)
    refs = []
    for line in run('otool', '-L', path).splitlines()[1:]:
        ref = line.strip().split(' (compatibility')[0]
        if ref and ref != own:
            refs.append(ref)
    return refs


def rpaths(path):
    out, result = run('otool', '-l', path).splitlines(), []
    for i, line in enumerate(out):
        if line.strip() == 'cmd LC_RPATH':
            result.append(out[i + 2].strip().split()[1])
    return result


def is_system(ref):
    return ref.startswith('/usr/lib/') or ref.startswith('/System/')


def resolve(ref, origin, executable_dir):
    """the file a reference of the file originally at <origin> loads"""
    def expand(p):
        return (p.replace('@loader_path', os.path.dirname(origin))
                 .replace('@executable_path', executable_dir))
    if ref.startswith('@rpath/'):
        for rp in rpaths(origin) + ['/opt/homebrew/lib', '/usr/local/lib']:
            cand = os.path.join(expand(rp), ref[len('@rpath/'):])
            if os.path.exists(cand):
                return os.path.realpath(cand)
        return None
    cand = expand(ref)
    return os.path.realpath(cand) if os.path.exists(cand) else None


def main():
    pkg = os.path.abspath(sys.argv[1])
    bin_dir, lib_dir = os.path.join(pkg, 'bin'), os.path.join(pkg, 'lib')
    deps_dir = os.path.join(lib_dir, 'deps')
    os.makedirs(deps_dir, exist_ok=True)
    ours = {n for n in os.listdir(lib_dir) if n.endswith('.dylib')}

    # (file in the package, where it was originally, kind)
    queue = [(os.path.join(bin_dir, n), os.path.join(bin_dir, n), 'bin') for n in sorted(os.listdir(bin_dir))]
    queue += [(os.path.join(lib_dir, n), os.path.join(lib_dir, n), 'lib') for n in sorted(ours)]
    queue = [q for q in queue if is_macho(q[0])]
    copied = {}   # dependency name -> original file
    done = set()
    while queue:
        path, origin, kind = queue.pop(0)
        if path in done:
            continue
        done.add(path)
        changes = []
        for ref in references(path):
            if is_system(ref):
                continue
            name = os.path.basename(ref)
            if name in ours:
                new = '@rpath/' + name
            else:
                if name not in copied:
                    src = resolve(ref, origin, bin_dir)
                    if src is None:
                        sys.exit(f'bundle_macos: cannot resolve {ref} of {path}')
                    dst = os.path.join(deps_dir, name)
                    run('cp', src, dst)
                    os.chmod(dst, 0o755)
                    copied[name] = src
                    queue.append((dst, src, 'deps'))
                new = {'bin': '@executable_path/../lib/deps/', 'lib': '@loader_path/deps/',
                       'deps': '@loader_path/'}[kind] + name
            if new != ref:
                changes += ['-change', ref, new]
        for rp in rpaths(path):
            changes += ['-delete_rpath', rp]
        if kind == 'bin':
            changes += ['-add_rpath', '@executable_path/../lib']
        elif kind == 'lib':
            changes += ['-id', '@rpath/' + os.path.basename(install_id(path) or path), '-add_rpath', '@loader_path']
        else:
            changes += ['-id', '@loader_path/' + os.path.basename(path)]
        run('install_name_tool', *changes, path)
        run('codesign', '--force', '--sign', '-', path)
        print(f'bundle_macos: {os.path.relpath(path, pkg)}')
    print(f'bundle_macos: {len(copied)} libraries in lib/deps')


if __name__ == '__main__':
    main()
