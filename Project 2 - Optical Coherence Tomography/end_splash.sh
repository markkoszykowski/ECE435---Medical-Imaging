#!/usr/bin/env bash

SWAP_FILE=/tmp/swap_file

sudo swapoff -v ${SWAP_FILE}
rm -vf ${SWAP_FILE}
