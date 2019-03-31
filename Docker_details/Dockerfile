#To build this docker file
#docker build -t pyg .

#To run this, with an interactive python2 temrinal, mounting your current host directory in the container directory at /workspace use:
#sudo docker run -it --rm -v`pwd`:/workspace pyg /bin/bash

#To run a jupyter notebook server with the python 2 conda environement
#sudo docker run -p 8888:8888 -it --rm -v`pwd`:/workspace pyg /bin/bash -c "source activate py2GEOL && jupyter notebook --allow-root --ip=0.0.0.0 --no-browser"

# Pull base image.
FROM ubuntu:trusty
MAINTAINER Nathaniel Butterworth

#Create some directories to work with
WORKDIR /build
RUN mkdir /workspace 

#Install ubuntu libraires and packages
RUN apt-get update -y
RUN apt-get install libglew-dev python-dev libboost-dev libboost-python-dev libboost-thread-dev libboost-program-options-dev libboost-test-dev libboost-system-dev libqt4-gui libqt4-dev libgdal1-dev libcgal-dev libproj-dev libqwt-dev libfreetype6-dev libfontconfig1-dev libxrender-dev libice-dev libsm-dev -y
RUN apt-get install -y wget

#Download pygplates and install it
RUN wget https://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/pygplates-ubuntu-trusty_1.5_1_amd64.deb
RUN dpkg -i pygplates-ubuntu-trusty_1.5_1_amd64.deb

#Download conda python because I like it
RUN wget -O ~/miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh  && \
     chmod +x ~/miniconda.sh && \
     ~/miniconda.sh -b -p /opt/conda && \
     rm ~/miniconda.sh 
     
RUN /opt/conda/bin/conda create -n py2GEOL python=2.7 scipy=1.2 scikit-learn=0.20 matplotlib=2.0 pyshp=1.2 numpy=1.15 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4

RUN /opt/conda/bin/conda create -n pyGEOL python=3.7 scipy=1.2 scikit-learn=0.20 matplotlib=3.0 pyshp=2.0 numpy=1.16 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4

#Make the workspace persistant
VOLUME /workspace
WORKDIR /workspace
EXPOSE 8888

#Set the path so python can find pygplates
ENV PYTHONPATH ${PYTHONPATH}:/usr/lib:/usr/lib/pygplates/revision12/
ENV PATH "/opt/conda/bin:$PATH"

RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh



