#!/usr/bin/sh

SWAP_FILE=/tmp/swap_file
SIZE_MB=30000
TO_RUN="matlab -nosplash -nodesktop -r run('./project3.m');"

dd if=/dev/zero of=${SWAP_FILE} bs=1M count=${SIZE_MB}
mkswap ${SWAP_FILE}
chmod 0600 ${SWAP_FILE}
sudo swapon -v ${SWAP_FILE}
echo "Swap enabled."
echo ${TO_RUN}
${TO_RUN}
echo "Removing the swap."
sudo swapoff -v ${SWAP_FILE}
rm -vf ${SWAP_FILE}
