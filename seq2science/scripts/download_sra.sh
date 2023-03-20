#!/bin/bash/
set -o nounset -e
# Downloads an SRA file to the current folder.

RUN=$1
LOG=$2

# Remove the top level folder since prefetch will assume we are done otherwise
if [ -e $RUN ]; then
    rm -r $RUN
fi

# three attempts
prefetch --max-size 999999999999 --output-directory ./ --log-level debug --progress $RUN >> $LOG 2>&1

# TODO: this is the strangest bug, in that on some machines (ocimum) prefetch downloads
# to a different location. Not sure what causes this, but this should fix that. Could simply
# be a global setting that we haven't discovered yet...
# bug report: https://github.com/ncbi/sra-tools/issues/533
if [[ -f "$RUN.sra" ]]; then
    mkdir -p $RUN
    mv $RUN.sra $RUN/
fi

# If an sralite file was downloaded instead of a sra file, just rename it
if [[ -f "$RUN.sralite" ]]; then
    mkdir -p $RUN
    mv $RUN.sralite $RUN/$RUN.sra
fi
