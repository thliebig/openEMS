# Build and test environment for the HIP GPU engine, e.g. on Vast.ai instances.
# Everything but openEMS itself: sync the sources to /workspace/openEMS-Project and
# run build-openems (see README.md).

# Vast.ai base image (ssh, portal), pinned: the "auto" tag resolves to different Ubuntu versions
ARG BASE=vastai/base-image:cuda-12.8.1-cudnn-devel-ubuntu24.04-2026-09-07
FROM ${BASE}

ENV DEBIAN_FRONTEND=noninteractive

# build dependencies of fparser, CSXCAD and openEMS (never an NVIDIA driver: the host provides it), btop
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git rsync time ca-certificates btop \
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

COPY build-openems run-test run-bench /opt/openEMS/tools/
COPY sitecustomize.py /opt/openEMS/tools/inject/
RUN chmod +x /opt/openEMS/tools/build-openems /opt/openEMS/tools/run-test /opt/openEMS/tools/run-bench

ENV PATH=/opt/openEMS/tools:/opt/openEMS/bin:/usr/local/cuda/bin:${PATH}
WORKDIR /workspace
