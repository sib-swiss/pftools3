FROM debian:buster
ENV DEBIAN_VERSION=buster-20200327

WORKDIR "/var/lib/pftools"
RUN DEBIAN_FRONTEND=noninteractive \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
      build-essential \
      libpcre++-dev \
      gfortran \
      ca-certificates \
      git \
      cmake \
      wget \
    && rm -fr /var/lib/apt/lists/* \
    && update-ca-certificates \
    && git clone "https://github.com/sib-swiss/pftools3.git" \
    && cd pftools3 \
    && cmake -DCMAKE_INSTALL_PREFIX:PATH=/var/lib/pftools -H. -Bbuild \
    && cmake --build build -- -j $(nproc)  \
    && cmake --build build -- install  \
    && cd .. \
    && rm -rf pftools3 \
    && apt-get remove -y --autoremove \
      build-essential \
      gfortran \
      ca-certificates \
      git \
      cmake \
      wget \
      libpcre++-dev \
    && apt-get autoremove \
    && apt-get clean
ENV PATH="/var/lib/pftools/bin/:${PATH}"
CMD ["bash"]
