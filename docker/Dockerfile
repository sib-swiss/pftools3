FROM debian:buster-slim as compile
ENV DEBIAN_VERSION=buster-20200327

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
      libfile-slurp-perl \
    && rm -fr /var/lib/apt/lists/* \
    && update-ca-certificates 
WORKDIR "/var/lib/pftools"
RUN git clone --depth=1 "https://github.com/sib-swiss/pftools3.git" \
    && cd pftools3 \
    && cmake -DCMAKE_INSTALL_PREFIX:PATH=/var/lib/pftools -H. -Bbuild \
    && cmake --build build -- -j $(nproc)  \
    && cmake --build build -- install 
WORKDIR "/var/lib/pftools/pftools3"
RUN cmake --build build -- test
FROM debian:buster-slim AS binary
COPY --from=compile /var/lib/pftools/bin/* /usr/local/bin/
CMD ["bash"]
