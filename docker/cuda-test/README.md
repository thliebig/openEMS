# CUDA test image

A Docker image for building and testing the HIP GPU engine on rented NVIDIA machines
(e.g. Vast.ai). It contains everything except openEMS itself:
- the build dependencies of fparser, CSXCAD and openEMS,
- a Python venv for the bindings and the tests,
- the scripts below, which sync, build and run.

It is based on the Vast.ai base image (ssh and portal), pinned to CUDA 12.8.1 and Ubuntu
24.04. CUDA 12.8 covers everything from Turing (RTX 20xx) to Blackwell (RTX 50xx), with an
NVIDIA driver 570 or newer on the host.

## Build and push

On an x86_64 Linux machine with Docker:

```
docker build -t seanmollet/openems-cuda-test:latest openEMS/docker/cuda-test
docker push seanmollet/openems-cuda-test:latest
```

## Use

1. Start the instances with the image `seanmollet/openems-cuda-test:latest` and ssh launch
   mode.

2. Sync the sources. fparser, CSXCAD and openEMS are submodules, and the CSXCAD build needs
   git describe, so also copy their git data from the parent repository:

   ```
   cd openEMS-Project
   for d in fparser CSXCAD openEMS; do
       rsync -az --exclude python/build/ $d/ root@<host>:/workspace/openEMS-Project/$d/
       rsync -az .git/modules/$d/ root@<host>:/workspace/openEMS-Project/.git/modules/$d/
   done
   ```

   Use `ssh -p <port>` with rsync `-e` for the instance's ssh port. Vast.ai closes
   connections when many short ones arrive in quick succession, so reuse one connection
   per host (`ControlMaster auto` in `~/.ssh/config`).

3. Build: `build-openems`

   - It builds into `/opt/openEMS`, for the architecture of the local GPU (`HIP_ARCH=86`
     overrides it).
   - It installs the Python bindings.
   - Rebuilds after the next sync are incremental.
   - Logs go to `/workspace/logs/build_*.log`.

4. Run:
   - `run-test GPU_Engine` runs `openEMS/python/Tests/GPU_Engine.py`.
   - `run-test Coax gpu` runs a test with all of its simulations on the GPU engine.
   - `run-bench <script.py> [engine]` runs a simulation script (default engine `gpu`) and
     records the output, the time, the peak memory and the peak GPU memory above the idle
     baseline in `/workspace/logs/<script>_<engine>.*`.

   The engine is chosen by a `sitecustomize.py` (in `/opt/openEMS/tools/inject`). It makes
   every `openEMS.Run()` use the engine given in `OPENEMS_TEST_ENGINE`, without changing
   the scripts.

`btop` is included for watching the machine.

## Local test

The Vast.ai entry point starts its services (supervisord) and ignores a command, so bypass it:

```
docker run --rm --entrypoint bash -e HIP_ARCH=86 -v <sources>:/workspace/openEMS-Project \
    <image> -c "build-openems && run-test Stripline"
```

Add `--gpus all` on a machine with an NVIDIA GPU; without one, set `HIP_ARCH`.
