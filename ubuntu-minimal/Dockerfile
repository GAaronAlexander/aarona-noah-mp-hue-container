# https://hub.docker.com/_/ubuntu/
# Use LTS version (currently 18.04)
FROM ubuntu:18.04
MAINTAINER @jupiterintel

# For more information about this base image, see
# https://blog.ubuntu.com/2018/07/09/minimal-ubuntu-released

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /opt

# Add Tini - https://github.com/krallin/tini/releases/download/v0.18.0/tini
COPY tini-0.18.0 /tini
RUN chmod +x /tini

# provide additional installation options that downstream
# containers can install, as required; these scripts should
# be on the PATH (i.e. copied to /usr/local/bin).
COPY scripts /usr/local/bin/

# install_base.sh enables the locale, set the ENV after that
# but before installing additional libs
RUN install_base.sh
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# ssl-certs
RUN install_openssh.sh
ENV CURL_CA_BUNDLE /etc/ssl/certs/ca-certificates.crt

# Install common HTTP tools
RUN install_web_tools.sh

COPY version /version

# interactive login shell
ENV DEBIAN_FRONTEND teletype
ENTRYPOINT ["/tini", "--"]
CMD ["/bin/bash", "-l"]
