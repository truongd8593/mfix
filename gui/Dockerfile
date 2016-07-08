# modified from https://github.com/ContinuumIO/docker-images/blob/master/anaconda/Dockerfile

FROM debian:8

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get update --fix-missing && apt-get install -y wget ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git gfortran automake autoconf make libgl1-mesa-glx libglu1-mesa libxt6 x11-utils xvfb openmpi-bin libopenmpi-dev
RUN apt-get install -y bzip2
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash /Anaconda2-4.0.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm /Anaconda2-4.0.0-Linux-x86_64.sh

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

RUN conda install -y vtk
RUN pip install coverage
RUN pip install future
RUN pip install pyinstaller
RUN pip install qtpy
RUN pip install xvfbwrapper

# http://bugs.python.org/issue19846
# > At the moment, setting "LANG=C" on a Linux system *fundamentally breaks Python 3*, and that's not OK.
ENV LANG C.UTF-8

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
