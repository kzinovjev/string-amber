FROM ubuntu:21.04

RUN apt-get update; \
    DEBIAN_FRONTEND=noninteractive; \
    apt-get install -y --no-install-recommends \
        bzip2=1.0.8-4ubuntu3 cmake=3.18.4-2ubuntu1 csh=20110502-6 \
        make=4.3-4ubuntu1 gcc=4:10.3.0-1ubuntu1 gfortran=4:10.3.0-1ubuntu1 \
        g++=4:10.3.0-1ubuntu1 flex=2.6.4-8 bison=2:3.7.5+dfsg-1 \
        patch=2.7.6-7 bc=1.07.1-2build2 libbz2-dev=1.0.8-4ubuntu3 \
        wget=1.21-1ubuntu3 openmpi-bin=4.1.0-7ubuntu2 \
        libopenmpi-dev=4.1.0-7ubuntu2 openssh-client=1:8.4p1-5ubuntu1 \
        ca-certificates=20210119build1 git=1:2.30.2-1ubuntu1

WORKDIR /usr/bin
COPY AmberTools21.tar.bz2 .
RUN tar xjvf AmberTools21.tar.bz2 && rm AmberTools21.tar.bz2

WORKDIR amber20_src
RUN git clone https://github.com/kzinovjev/string-amber
RUN cp string-amber/lapack/* AmberTools/src/lapack/
RUN cp string-amber/sander/* AmberTools/src/sander/
RUN cp string-amber/docker/run_cmake build/

WORKDIR build
RUN chmod +x run_cmake
RUN ./run_cmake
RUN make -j 4 install
RUN echo "source /usr/bin/amber20/amber.sh" >> /etc/bash.bashrc
SHELL ["/bin/bash", "-c"]
