FROM debian:stable-20210721-slim

LABEL maintainer "Daniel Cordier (https://orcid.org/0000-0003-4515-6271)"
LABEL description "Fortran 2008 implementation of the PC-SAFT EoS for Titan"

RUN apt update && \
    apt install -y gfortran \
                   make

WORKDIR /opt/src/

# Copy source code
COPY mod_pcsaft.f08 .
COPY utils_dc.f08 .
COPY foul.f90 .

COPY binary_N2CH4_demo.f08 .
COPY pcsaft_demo.f08 .

COPY Makefile .

RUN make && \
    make clean && \
    mv pcsaft_demo /usr/bin/. && \
    mv binary_N2CH4_demo /usr/bin/.

WORKDIR /data
