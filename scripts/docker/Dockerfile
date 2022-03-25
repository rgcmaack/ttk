FROM ubuntu:focal as base
LABEL maintainer="Christoph Garth <garth@cs.uni-kl.de>"

ENV DEBIAN_FRONTEND=noninteractive

# --------------------------------------------------------------------------

FROM base as build-base

ENV CMAKE_BUILD_TYPE=MinSizeRel \
    CMAKE_GENERATOR=Ninja

# bring in additional apt sources
RUN apt-get update \ 
 && apt-get install --no-install-recommends -yqq ca-certificates \
 && echo "deb [trusted=yes] https://apt.kitware.com/ubuntu/ focal main" > /etc/apt/sources.list.d/kitware.list \
 && apt-get update

# install base development env
RUN apt-get install --no-install-recommends -yqq \
    build-essential \
    ninja-build \
    cmake \
    dlocate \
    file \
    curl

# --------------------------------------------------------------------------

FROM build-base as builder

COPY scripts/docker/install-helper /usr/bin

# install OSPRay + dependencies
COPY scripts/docker/pkg/ispc.sh /tmp
RUN  install-helper /tmp/ispc.sh

COPY scripts/docker/pkg/tbb.sh /tmp
RUN  install-helper /tmp/tbb.sh

COPY scripts/docker/pkg/embree.sh /tmp
RUN  install-helper /tmp/embree.sh

COPY scripts/docker/pkg/rkcommon.sh /tmp
RUN  install-helper /tmp/rkcommon.sh

COPY scripts/docker/pkg/openimagedenoise.sh /tmp
RUN  install-helper /tmp/openimagedenoise.sh

COPY scripts/docker/pkg/openvkl.sh /tmp
RUN  install-helper /tmp/openvkl.sh

COPY scripts/docker/pkg/ospray.sh /tmp
RUN  install-helper /tmp/ospray.sh

# install OSMesa
COPY scripts/docker/pkg/mesa.sh /tmp
RUN  install-helper /tmp/mesa.sh

# install ZFP
COPY scripts/docker/pkg/zfp.sh /tmp
RUN  install-helper /tmp/zfp.sh

# install Spectra
COPY scripts/docker/pkg/spectra.sh /tmp
RUN  install-helper /tmp/spectra.sh

# install ParaView
ARG paraview=5.10.0
ENV PARAVIEW_VERSION=${paraview}

COPY scripts/docker/pkg/paraview.sh /tmp
RUN  install-helper /tmp/paraview.sh

# --------------------------------------------------------------------------

FROM builder as builder-ttk

COPY scripts/docker/pkg/ttk-local.sh /tmp
COPY . /build/ttk-local
RUN  install-helper /tmp/ttk-local.sh

# --------------------------------------------------------------------------

FROM base as ttk

COPY --from=builder-ttk /usr/local /usr/local

RUN apt-get update \
 && apt-get -yqq --no-install-recommends install \
        $(cat /usr/local/.pkgs) \
        libgoogle-perftools-dev \
        google-perftools \
        binutils \
 && apt-get clean \
 && rm -rf /var/cache/apt/lists

# run pvserver by default
CMD /usr/local/bin/pvserver
EXPOSE 11111

# --------------------------------------------------------------------------

FROM build-base as ttk-dev

COPY --from=builder /usr/local /usr/local

RUN apt-get update \
 && apt-get -yqq --no-install-recommends install $(cat /usr/local/.pkgs) \
 && apt-get clean \
 && rm -rf /var/cache/apt/lists

# run pvserver by default
CMD /usr/local/bin/pvserver
EXPOSE 11111
