FROM ubuntu:22.04

RUN apt-get update; \
    DEBIAN_FRONTEND=noninteractive; \
    apt-get install -y --no-install-recommends \
        bzip2=1.0.8-5build1 cmake=3.22.1-1ubuntu1.22.04.1 csh=20110502-7 \
        make=4.3-4.1build1 gcc=4:11.2.0-1ubuntu1 gfortran=4:11.2.0-1ubuntu1 \
        g++=4:11.2.0-1ubuntu1 flex=2.6.4-8build2 bison=2:3.8.2+dfsg-1build1 \
        patch=2.7.6-7build2 bc=1.07.1-3build1 libbz2-dev=1.0.8-5build1 \
        wget=1.21.2-2ubuntu1 openmpi-bin=4.1.2-2ubuntu1 \
        libopenmpi-dev=4.1.2-2ubuntu1 openssh-client=1:8.9p1-3 \
        ca-certificates=20211016ubuntu0.22.04.1 git=1:2.34.1-1ubuntu1.5

WORKDIR /usr/bin
COPY AmberTools23.tar.bz2 .
RUN tar xjvf AmberTools23.tar.bz2 && rm AmberTools23.tar.bz2

WORKDIR amber22_src
RUN git clone https://github.com/kzinovjev/string-amber
RUN cp string-amber/lapack/* AmberTools/src/lapack/
RUN cp string-amber/sander/* AmberTools/src/sander/
RUN cp string-amber/docker/run_cmake build/

WORKDIR build
RUN chmod +x run_cmake
RUN ./run_cmake
RUN make -j 4 install
RUN echo "source /usr/bin/amber22/amber.sh" >> /etc/bash.bashrc
SHELL ["/bin/bash", "-c"]
