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
    && rm -fr /var/lib/apt/lists/* \
    && update-ca-certificates \
    && git clone "https://github.com/sib-swiss/pftools3.git" \
    && cd pftools3 \
    && cmake -H. -Bbuild \
    && cmake --build build -- -j $(nproc) \
    && ln -s /var/lib/pftools/pftools3/build/src/C/pfscanV3 /usr/bin \
    && ln -s /var/lib/pftools/pftools3/build/src/C/pfsearchV3 /usr/bin \
    && apt-get remove -y --autoremove \
      build-essential \
      gfortran \
      ca-certificates \
      git \
      cmake \
    && apt-get clean

CMD ["bash"]
