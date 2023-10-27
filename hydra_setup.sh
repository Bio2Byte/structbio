#!/usr/bin/bash
#
# This script creates a couple of soft links to the repository location
# to facilitate working on hydra.

check_success() {
    ERR="$?"
    if [ $ERR -ne 0 ]; then
        echo "[ERROR] $*"
        exit $ERR
    fi
}

# Check that you are really working on hydra
test -d "$VSC_HOME"; check_success "VSC_HOME is not defined. Not working on hydra?"
test -d "$VSC_DATA"; check_success "VSC_DATA is not defined. Not working on hydra?"
test -d "$VSC_SCRATCH"; check_success "VSC_SCRATCH is not defined. Not working on hydra?"

# Get important file locations
SCRIPT_PATH="$(readlink -f $0)"
REPO_PATH="$(dirname $SCRIPT_PATH)"
REPO_NAME="$(basename $REPO_PATH)"

# Create softlinks to the repository from all entry points on hydra
ln -sv "$REPO_PATH" "$VSC_HOME/$REPO_NAME"
ln -sv "$REPO_PATH" "$VSC_DATA/$REPO_NAME"
