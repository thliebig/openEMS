# openEMS images

Two images from plain Ubuntu 24.04 (x86_64), one Dockerfile:

- **dev** (`seanmollet/openems-dev`, ~3.2 GB): the build tools, the CUDA 12.8
  compiler and the dependencies of fparser, CSXCAD and openEMS, a Python venv, and
  scripts to build and test. openEMS itself is built from synced sources.
- **runtime** (`seanmollet/openems`, ~380 MB): openEMS, nf2ff and the Python bindings,
  built in dev, with a venv of numpy, h5py and matplotlib. No compilers or headers: it
  has only the shared libraries openEMS loads (copied from the build, see
  `collect-runtime-libs`) and Python from apt.

The HIP engine is built for Turing to Blackwell (`75-real;80-real;86-real;89-real;
90-real;100-real;120`, the build argument `HIP_ARCHITECTURES`), with PTX of the newest
for later GPUs: one image for all NVIDIA GPUs. Both images have btop, built from source
with GPU monitoring.

Neither image contains an NVIDIA driver. The NVIDIA Container Toolkit of the host mounts
the host's driver libraries and `nvidia-smi` (`NVIDIA_DRIVER_CAPABILITIES=compute,utility`),
so they always match its kernel module. Built with CUDA 12.8: drivers 570 and newer run
it; older CUDA 12 drivers (525+) should through CUDA's minor version compatibility
(untested).

## Build

The build context is the openEMS-Project directory with the fparser, CSXCAD and openEMS
submodules and their git data (`.git/modules`), which the version numbers come from:

```
cd openEMS-Project
docker build -f openEMS/docker/openems/Dockerfile --target dev     -t seanmollet/openems-dev .
docker build -f openEMS/docker/openems/Dockerfile --target runtime -t seanmollet/openems .
```

## Use

Runtime:

```
docker run --rm --gpus all -v $PWD:/workspace seanmollet/openems python my_simulation.py
```

`python` is the venv's, and `openEMS` and `nf2ff` are on the path.

Dev:

1. Sync the sources to `/workspace/openEMS-Project` (fparser, CSXCAD, openEMS and
   `.git/modules/<name>` of each), e.g. with rsync.
2. `build-openems`: builds into `/opt/openEMS` for the local GPU (`HIP_ARCH=86`
   overrides it) and installs the Python bindings. Rebuilds are incremental.
3. `run-test <Test> [engine]` runs `openEMS/python/Tests/<Test>.py`, optionally with
   every simulation on the given engine (e.g. `gpu`). `run-bench <script.py> [engine]`
   records the run time and peak host and GPU memory of a script.
