FROM ubuntu:20.04 AS compile

ARG PFTOOLS_VERSION
ARG BUILD_DATE
ENV PFTOOLS_PATH=/var/lib/pftools

RUN echo '# OS update' \
	&& export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update -y \
	&& apt-get install -y --no-install-recommends \
      build-essential \
      libpcre2-dev \
      gfortran \
      libgfortran5 \
      ca-certificates \
      git \
      cmake \
      zlib1g-dev \
      libpng-dev \
      libfile-slurp-perl \
	&& echo '# Install Pftools V3' \
	&& git clone --depth=1 "https://github.com/sib-swiss/pftools3.git" \
	&& mkdir pftools3/build \
	&& cd pftools3/build \
	&& cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$PFTOOLS_PATH -DCMAKE_BUILD_TYPE=Release -DUSE_GRAPHICS=OFF -DUSE_PDF=ON \
	&& make \
	&& make install \
	&& make test \
	&& gzip -9 $PFTOOLS_PATH/share/man/man5/*.5 $PFTOOLS_PATH/share/man/man1/*.1


FROM ubuntu:20.04 AS binary

ARG PFTOOLS_VERSION
ARG BUILD_DATE
ENV PFTOOLS_PATH=/var/lib/pftools

# METADATA
# Use Open Containers Initiative (OCI)
# See https://rehansaeed.com/docker-labels-depth/
# See https://github.com/opencontainers/image-spec/blob/master/annotations.md
# Exist also another structured label initiative: http://label-schema.org/ superceded by OCI now
LABEL org.opencontainers.image.title="pftools3 - A suite of tools to build and search generalized profiles"
LABEL org.opencontainers.image.version=$PFTOOLS_VERSION
LABEL org.opencontainers.image.vendor="SIB Swiss Institute of Bioinformatics"
LABEL org.opencontainers.image.authors="sebastien.moretti@sib.swiss"
LABEL org.opencontainers.image.url="https://github.com/sib-swiss/pftools3"
LABEL org.opencontainers.image.documentation="https://github.com/sib-swiss/pftools3"
LABEL org.opencontainers.image.licenses="GPL-2.0"
LABEL org.opencontainers.image.description="pftools3 - A suite of tools to build and search generalized profiles"
LABEL org.opencontainers.image.created=$BUILD_DATE

COPY --from=compile $PFTOOLS_PATH/ $PFTOOLS_PATH/
RUN echo '# OS update' \
	&& export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update -y \
	&& apt-get install -y --no-install-recommends \
	   man-db \
	   ca-certificates \
	   libfile-slurp-perl \
	   libc6 \
	   libgcc1 \
	   libgfortran5 \
	   libpcre2-8-0 \
	   libquadmath0 \
	&& ln -s $PFTOOLS_PATH/bin/* /usr/local/bin/ \
	&& ln -s $PFTOOLS_PATH/share/man/man5/* /usr/share/man/man5/ \
	&& ln -s $PFTOOLS_PATH/share/man/man1/* /usr/share/man/man1/ \
	&& mandb \
	&& apt autoremove -y \
	&& apt-get clean \
	&& rm -fr /var/lib/apt/lists/* \
	&& update-ca-certificates


# SECURITY
## Control root access
ENV USER=pftools3
RUN groupadd -g 1001 ${USER} && useradd -r -u 1001 -g ${USER} ${USER}
USER ${USER}
ENV HOME=/home/${USER}
WORKDIR ${HOME}
##Rest of Dockerfile with this user
ENV LC_ALL="C"


#TODO man pages do not work???
# TEST
RUN for pl in $PFTOOLS_PATH/bin/*.pl; do perl -cw $pl; done \
	&& pfscanV3 -h \
	&& pfscan -h || true \
	&& echo \
	&& pfsearchV3 -h \
	&& pfsearch -h || true


HEALTHCHECK CMD pfscanV3 -h || exit 1


CMD ["bash"]

