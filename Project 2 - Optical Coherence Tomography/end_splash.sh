#!/usr/bin/sh

SWAP_FILE=/tmp/swap_file

sudo swapoff -v ${SWAP_FILE}
rm -vf ${SWAP_FILE}
