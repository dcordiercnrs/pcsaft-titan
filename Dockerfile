# Dockerfile make an image dedicated to PCSAFT-TITAN.
# D. Cordier, 22 juillet 2021.
# http://orcid.org/0000-0003-4515-6271
# https://www.researchgate.net/profile/Daniel_Cordier
FROM debian:10-slim
RUN apt-get update && \
    apt-get upgrade --yes && \  
    apt-get install make --yes && \ 
    apt-get install git --yes && \
    apt-get install gfortran --yes
#
# PC-SAFT Titan installation:
WORKDIR /home
CMD ls
