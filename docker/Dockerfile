# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    g++ \
    gfortran \
    git \
    python3 \
    python3-pip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Note: For legacy Fortran builds (DELPHIPP_ENABLE_LEGACY_FORTRAN=ON),
# mount or install CERNLIB 2006 / nypatchy (e.g. via CERN CVMFS
# /cvmfs/delphi.cern.ch or Debian/Fedora cernlib packages).

WORKDIR /workspace/delphipp
