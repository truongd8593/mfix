# modified from https://github.com/ContinuumIO/docker-images/blob/master/anaconda/Dockerfile

FROM continuumio/anaconda:4.2.0

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

# Install gfortran for building MFIX
RUN apt-get update --fix-missing && apt-get install -y \
    autoconf \
    automake \
    gfortran \
    openmpi-bin \
    make

# Install XVFB and dependencies for unit tests
RUN apt-get update --fix-missing && apt-get install -y \
    libgl1-mesa-glx \
    libglu1-mesa \
    libopenmpi-dev \
    libxt6 \
    x11-utils \
    xvfb

RUN conda install -y vtk
RUN pip install coverage
RUN pip install future
RUN pip install packaging
RUN pip install pyinstaller==3.1.1
RUN pip install qtpy
RUN pip install xvfbwrapper
