#!/bin/bash

apt-get -yq update
apt-get install -y --no-install-recommends apt-utils
apt-get -yq upgrade

# Set locale
apt-get install -yq --no-install-recommends locales

echo "LC_ALL=en_US.UTF-8" >> /etc/environment
echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
echo "LANG=en_US.UTF-8" > /etc/locale.conf
locale-gen en_US.UTF-8
dpkg-reconfigure --frontend noninteractive locales
# not sure about doing the following:
#echo "America/Los_Angeles" > /etc/timezone
#dpkg-reconfigure --frontend noninteractive tzdata

# git is commonly required to install dependencies, so include it here.
# The majority of the build tools are optional in this base image, to install
# them using a common pattern, use the /usr/src/scripts/install_build_tools.sh
apt-get install -yq --no-install-recommends \
    git

# Clean up
apt-get autoclean -y
apt-get autoremove -y
apt-get clean -y -q
rm -rf /var/cache/apt/archives/*
rm -rf /var/lib/apt/lists/*
rm -rf /var/tmp/*
rm -rf /tmp/*
