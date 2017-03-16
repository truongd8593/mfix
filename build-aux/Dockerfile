# Docker image for building MFIX on Linux

FROM continuumio/conda_builder_linux
ENV PATH="/opt/miniconda/bin:${PATH}"

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>


RUN conda install -y \
    automake \
    conda-build \
    nomkl \
    numpy==1.11.3 \
    pandoc \
    pyqt \
    qtpy


# RUN conda install -y vtk
# RUN pip install coverage
# RUN pip install future
# RUN pip install packaging
# RUN pip install pyinstaller==3.1.1
# RUN pip install qtpy
# RUN pip install xvfbwrapper
