#!/bin/bash

REPO="holab-hku/vOMIX-MEGA" # Replace with your GitHub repository
DIR=$1

mkdir -p tmp
cd tmp

echo "Cloning repository $REPO..."
git clone "https://github.com/$REPO.git"

if [ $? -ne 0 ]; then
    echo "Error cloning repository $REPO."
    exit 1
fi

cd ${SOFTWARE}


