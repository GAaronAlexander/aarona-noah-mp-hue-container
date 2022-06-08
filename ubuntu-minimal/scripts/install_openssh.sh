#!/bin/bash

apt-get -yq update
apt-get install -y --no-install-recommends apt-utils

# ssh is required for docker builds with ssh bind mounts and ssh-keyscan
# ssl libs may be required for cryptography and secure sockets etc.
apt-get install -yq --no-install-recommends \
    git \
    ca-certificates \
    openssh-client \
    openssl

update-ca-certificates
export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
git config --global http.sslVerify true
git config --global http.sslCAinfo /etc/ssl/certs/ca-certificates.crt
git config --global http.sslBackend "openssl"

# Clean up
apt-get autoclean -y
apt-get autoremove -y
apt-get clean -y -q
rm -rf /var/cache/apt/archives/*
rm -rf /var/tmp/*
rm -rf /tmp/*
