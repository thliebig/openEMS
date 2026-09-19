# openEMS images with the HIP GPU engine, from plain Ubuntu (see README.md). Targets:
#   dev:     the build tools, the CUDA compiler and the dependencies of openEMS; sync the
#            sources and run build-openems
#   runtime: openEMS and its Python bindings, built in dev for all supported GPU
#            architectures, with only the libraries they load; no build tools
# Build context: the openEMS-Project directory (fparser, CSXCAD, openEMS, .git/modules).
#
# No NVIDIA driver in the images: the NVIDIA Container Toolkit of the host mounts its driver
# libraries and nvidia-smi (NVIDIA_DRIVER_CAPABILITIES), matching the host's kernel module.

ARG UBUNTU=ubuntu:24.04
ARG HIP_VERSION=12-8
ARG BTOP_VERSION=v1.4.7
# machine code for Turing to Blackwell, and PTX of the newest for later GPUs
ARG HIP_ARCHITECTURES="75-real;80-real;86-real;89-real;90-real;100-real;120"

FROM ${UBUNTU} AS base
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

# btop from source, with GPU monitoring (it loads the NVML of the host driver). It needs C++23:
# GCC 14 (Ubuntu 24.04's libstdc++ is from GCC 14, so the runtime image has what it links).
FROM base AS btop
ARG BTOP_VERSION
RUN apt-get update && apt-get install -y --no-install-recommends build-essential g++-14 git ca-certificates \
    && git clone --depth 1 --branch ${BTOP_VERSION} https://github.com/aristocratos/btop /tmp/btop \
    && make -C /tmp/btop -j$(nproc) GPU_SUPPORT=true CXX=g++-14 \
    && make -C /tmp/btop install PREFIX=/usr/local DESTDIR=/out \
    && strip /out/usr/local/bin/btop

FROM base AS dev
ARG HIP_VERSION
# the CUDA compiler and runtime from NVIDIA's repository, the build dependencies of fparser,
# CSXCAD and openEMS
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
    && curl -fsSL -o /tmp/cuda-keyring.deb \
        https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb \
    && dpkg -i /tmp/cuda-keyring.deb && rm /tmp/cuda-keyring.deb \
    && apt-get update && apt-get install -y --no-install-recommends \
        cuda-hipcc-${HIP_VERSION} cuda-cudart-dev-${HIP_VERSION} \
        build-essential cmake git rsync time \
        libhdf5-dev libvtk9-dev libcgal-dev libtinyxml-dev libgmp-dev libmpfr-dev \
        libboost-program-options-dev libboost-thread-dev libboost-date-time-dev \
        libboost-serialization-dev libboost-chrono-dev libboost-system-dev \
        python3-dev python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Python environment for the openEMS/CSXCAD bindings and the tests
RUN python3 -m venv /opt/openEMS/venv \
    && /opt/openEMS/venv/bin/pip install --no-cache-dir --upgrade \
        pip setuptools setuptools_scm cython numpy h5py matplotlib scipy

# the synced submodules are owned by another user
RUN git config --system --add safe.directory '*'

COPY --from=btop /out/usr/local/ /usr/local/
COPY openEMS/docker/openems/build-openems openEMS/docker/openems/run-test openEMS/docker/openems/run-bench \
     openEMS/docker/openems/collect-runtime-libs /opt/openEMS/tools/
COPY openEMS/docker/openems/sitecustomize.py /opt/openEMS/tools/inject/

ENV PATH=/opt/openEMS/tools:/opt/openEMS/bin:/usr/local/cuda/bin:${PATH}
WORKDIR /workspace

# what the runtime image gets from apt: its libraries are not copied (see collect-runtime-libs)
FROM base AS runtime-base
RUN apt-get update && apt-get install -y --no-install-recommends python3 \
    && rm -rf /var/lib/apt/lists/*

FROM runtime-base AS runtime-packages
RUN dpkg-query -W -f='${Package}\n' > /runtime.packages

# openEMS for the runtime image: the libraries and programs, a venv with the bindings and
# their runtime packages only, and the shared libraries they load that apt does not provide there
FROM dev AS openems-build
ARG HIP_ARCHITECTURES
COPY fparser /src/fparser
COPY CSXCAD /src/CSXCAD
COPY openEMS /src/openEMS
COPY .git/modules /src/.git/modules
RUN HIP_ARCH="${HIP_ARCHITECTURES}" BUILD_PYTHON=0 build-openems /src
RUN set -e; P=/opt/openEMS; \
    for py in CSXCAD openEMS; do \
        cd /src/$py/python && rm -rf build; \
        CSXCAD_INSTALL_PATH=$P OPENEMS_INSTALL_PATH=$P $P/venv/bin/pip wheel --no-build-isolation --no-deps -w /wheels .; \
        $P/venv/bin/pip install --no-deps /wheels/$(echo $py | tr A-Z a-z)-*.whl; \
    done; \
    rm -rf $P/venv; \
    python3 -m venv $P/venv; \
    $P/venv/bin/pip install --no-cache-dir numpy h5py matplotlib /wheels/*.whl; \
    # only our files: stripping the libraries bundled in the wheels can break them (patchelf)
    for f in $(find $P/bin $P/lib $P/venv/lib/python3*/site-packages/CSXCAD $P/venv/lib/python3*/site-packages/openEMS \
               -type f \( -name '*.so*' -o -perm -u+x \)); do \
        strip --strip-unneeded $f 2>/dev/null || true; \
    done
COPY --from=runtime-packages /runtime.packages /tmp/
RUN collect-runtime-libs /tmp/runtime.packages /opt/openEMS/deps /opt/openEMS/bin /opt/openEMS/lib /opt/openEMS/venv

FROM runtime-base AS runtime
COPY --from=btop /out/usr/local/ /usr/local/
COPY --from=openems-build /opt/openEMS/bin /opt/openEMS/bin
COPY --from=openems-build /opt/openEMS/lib /opt/openEMS/lib
COPY --from=openems-build /opt/openEMS/deps /opt/openEMS/deps
COPY --from=openems-build /opt/openEMS/venv /opt/openEMS/venv
# every library resolves (the CUDA driver is loaded at run time), and the bindings import
RUN printf '/opt/openEMS/lib\n/opt/openEMS/deps\n' > /etc/ld.so.conf.d/openems.conf && ldconfig \
    && ! find /opt/openEMS -type f \( -name '*.so*' -o -perm -u+x \) -exec ldd {} + 2>/dev/null | grep "not found" \
    && cd /tmp && /opt/openEMS/venv/bin/python -c "import CSXCAD, openEMS, h5py, matplotlib"
ENV PATH=/opt/openEMS/venv/bin:/opt/openEMS/bin:${PATH}
WORKDIR /workspace
