#!/bin/bash

apt-get -yq update

# Install common HTTP tools
apt-get install -yq --no-install-recommends \
    curl \
    wget

# Clean up
apt-get autoclean -y
apt-get autoremove -y
apt-get clean -y -q
rm -rf /var/cache/apt/archives/*
rm -rf /var/lib/apt/lists/*
rm -rf /var/tmp/*
rm -rf /tmp/*
