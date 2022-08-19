# syntax=docker/dockerfile:1.0.2-experimental
# Start first stage of build
#FROM registry.gitlab.com/jupiterintel/baseimages/ubuntu-minimal:0.1.0 AS base
FROM ubuntu-minimal

# From here on, this is now Callies docker file!

#FROM registry.gitlab.com/jupiterintel/baseimages/centos7:0.1.1
LABEL maintainer="Callie McNicholas <callie.mcnicholas@jupiterintel.com>"
ARG ROOT_PATH=/home/jupiter/model/wrf
ARG APP_PATH=${ROOT_PATH}/container

# This Dockerfile compiles WRF from source during "docker build" step
ENV WRF_VERSION=4.3.3 \
    WPS_VERSION=4.3.1 \
    CONDA_VERSION=latest
# Avoid writting pyc or pyo files to save space
ENV PYTHONDONTWRITEBYTECODE=true

# 
RUN apt-get update && \
    apt-get install -y curl file
# Get the license from NCAR
RUN curl -SL https://ral.ucar.edu/sites/default/files/public/projects/ncar-docker-wrf/ucar-bsd-3-clause-license.pdf > /UCAR-BSD-3-Clause-License.pdf
# Setup input/output directories for WRF as well as directories for input model data
RUN mkdir /home/jupiter && mkdir /home/jupiter/model && mkdir ${ROOT_PATH} && \
    mkdir ${ROOT_PATH}/static_data && mkdir ${ROOT_PATH}/wpsdata && \
    mkdir ${ROOT_PATH}/wpsdata/geogrid && mkdir ${ROOT_PATH}/wpsdata/intermediate && \
    mkdir ${ROOT_PATH}/wpsdata/metgrid && mkdir ${ROOT_PATH}/wrfdata && \
    mkdir ${ROOT_PATH}/wrfdata/real && mkdir ${ROOT_PATH}/wrfdata/wrf && \
    mkdir ${ROOT_PATH}/data_sources && mkdir ${ROOT_PATH}/data_sources/ERA5

#  Copy container files
#COPY configs ${APP_PATH}/configs
#COPY input_models ${APP_PATH}/input_models
#COPY main ${APP_PATH}/main
#COPY tests ${APP_PATH}/tests
#COPY utilities ${APP_PATH}/utilities
#COPY * ${APP_PATH}/
COPY conda_requirements.txt ${APP_PATH}/conda_requirements.txt
COPY docker-clean ${APP_PATH}/docker-clean



# Retrieve miniconda from URL
RUN curl -o ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh && \
   bash ~/miniconda.sh -b -p /miniconda && \
   export PATH="/miniconda/bin:$PATH" && \
   rm ~/miniconda.sh

ENV CONDA_DIR=/miniconda
ENV PATH=/miniconda/bin:$PATH
# Add condarc file and set PATH variable
COPY condarc /root/.condarc

# Install conda libraries and remove unnecessary files
RUN conda update -n base -c defaults conda && \
    conda install -y -c conda-forge --file ${APP_PATH}/conda_requirements.txt && \ 
    conda clean --all -yq && \
    find ${CONDA_DIR} -type f -name '*.pyc' -exec rm -f {} \; && \
    find ${CONDA_DIR} -type f -name '*.js.map*' -exec strip --strip-all {} \; && \
    rm -rf ${CONDA_DIR}/pkgs && \
    rm -rf ${CONDA_DIR}/lib/*_mc.so && \
    rm -rf ${CONDA_DIR}/lib/*_mc2.so && \
    rm -rf ${CONDA_DIR}/lib/*_mc3.so && \
    rm -rf ${CONDA_DIR}/lib/*_avx512* && \
    rm -rf ${CONDA_DIR}/lib/*_avx.* 

# Install pip libraries
#RUN mkdir -p -m 0600 ~/.ssh && ssh-keyscan gitlab.com >> ~/.ssh/known_hosts
#RUN --mount=type=ssh bash -c "mkdir -p -m 0600 ~/.ssh && ssh-keyscan gitlab.com >> ~/.ssh/known_hosts"
#RUN --mount=type=ssh pip install --ignore-installed --no-cache-dir -r ${APP_PATH}/pip_requirements.txt
#RUN --mount=type=ssh pip install --ignore-installed --no-cache-dir -e ${APP_PATH}

# Retrieve WRF and WPS from github
RUN curl -SL https://github.com/wrf-model/WRF/archive/refs/tags/v${WRF_VERSION}.tar.gz | tar zxC ${ROOT_PATH} \
    && mv ${ROOT_PATH}/WRF-${WRF_VERSION} ${ROOT_PATH}/WRF \
    && curl -SL https://github.com/wrf-model/WPS/archive/refs/tags/v${WPS_VERSION}.tar.gz | tar zxC ${ROOT_PATH} \
    && mv ${ROOT_PATH}/WPS-${WPS_VERSION} ${ROOT_PATH}/WPS \
    && curl -SL https://ral.ucar.edu/sites/default/files/public/projects/ncar-docker-wrf/ucar-bsd-3-clause-license.pdf > /UCAR-BSD-3-Clause-License.pdf

# Start second stage of build
#FROM registry.gitlab.com/jupiterintel/baseimages/ubuntu-minimal:0.1.0 AS refactor_batch
##FROM ubuntu:18.04 AS refactor_batch

# Install relevant libraries (added scp Aaron A.)
RUN apt-get update && apt-get install -y --no-install-recommends unzip csh make m4 file && \
    apt-get autoclean -y && \
    apt-get autoremove -y && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/archives/* && \
    rm -rf /var/tmp/* && \
    rm -rf /tmp/*

# Respecify args
ARG ROOT_PATH=/home/jupiter/model/wrf
ARG APP_PATH=${ROOT_PATH}/container

ENV PATH=/miniconda/bin:$PATH
# Add condarc file and set PATH variable
COPY condarc /root/.condarc

# Copy files from base
#COPY --from=base /miniconda /miniconda
#COPY --from=base /home/jupiter /home/jupiter

# On to compiling WRF
WORKDIR ${ROOT_PATH}

# Set those environment variables here for compile time
ENV PATH=/miniconda/bin:$PATH \
    CC=gcc \
    FC=gfortran \
    CXX=g++ \
    HDF5=/miniconda \
    NETCDF=/miniconda \
    JASPERLIB=/miniconda/lib \
    JASPERINC=/miniconda/include

# Set environment variables, compile WRF, compile WPS, and remove extras.
RUN export NETCDF=${NETCDF} \
    && export JASPERINC=${JASPERINC} \
    && export JASPERLIB=${JASPERLIB} \
    && export HDF5=${HDF5} \
    && export CC=${CC} \
    && export CXX=${CXX}} \
    && export FC=${FC} \
    && cd ${ROOT_PATH}/WRF \
    && echo '35\r5\r' > config.in \
    && ./configure < config.in \
    && sed -i -e 's?/lib/cpp?/miniconda/bin/cpp?' ./configure.wrf \
    && ./compile em_real \
    && cd ${ROOT_PATH}/WPS \
    && echo '1\r' > config.in \
    && ./configure < config.in \
    && sed -i -e 's?/usr/bin/cpp?/miniconda/bin/cpp?' ./configure.wps \
    && sed -i -e 's?LDFLAGS.*?LDFLAGS             = -lgomp?' ./configure.wps \ 
    && ./compile \ 
    && chmod +x ${APP_PATH}/docker-clean \
    && ${APP_PATH}/docker-clean

#    && mv ${APP_PATH}/tests/test_files/test_run_files/*.grib ${ROOT_PATH}/data_sources/ERA5/ \
#    && scp ${ROOT_PATH}/WPS/ungrib/Variable_Tables/Vtable.ERA-interim.pl ${ROOT_PATH}/WPS/ungrib/Variable_Tables/Vtable.ERA5 \
#    && rm -rf ${APP_PATH}/tests/test_files/test_run_files \
#    && chmod +x ${APP_PATH}/docker-clean \
#    && ${APP_PATH}/docker-clean

## Everything below has been added by Aaron A. 
# Adding a new "layer"
ARG NOAHMP=/home/jupiter/model/noahmp 

# create a new layer to download files into
RUN mkdir $NOAHMP \ 
    && mkdir -p /var/run/secrets/eks.amazonaws.com/serviceaccount \
    && pip install awscli
 
#Grab this from Aaron A.'s GITHUB \
#COPY ./token /var/run/secrets/eks.amazonaws.com/serviceaccount/token
RUN git clone https://github.com/GAaronAlexander/NOAH-MP_HUE.git ${NOAHMP} \
    && cd ${NOAHMP}/hrldas  \
    && rm user_build_options \
    && ln -s user_build_options.gfortran.cloud.parallel user_build_options \ 
    && make clean \
    && make 

COPY ./token /var/run/secrets/eks.amazonaws.com/serviceaccount/token
COPY ./geogrid-files/* ${NOAHMP}/geogrid-files/ 
COPY *.py ${NOAHMP}/

#you need to make sure you replace your 
#.cdsapric file to have updated permissions
COPY ./.cdsapirc /root/.cdsapirc

#control your physics options here
COPY namelist.hrldas.draft ${NOAHMP}/hrldas/run/


# Start bash
CMD ["/bin/bash"]

