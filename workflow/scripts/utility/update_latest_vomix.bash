#!/bin/bash

SOFTWARE="vOMIX-MEGA"
REPO="holab-hku/${SOFTWARE}" # Replace with your GitHub repository

mkdir -p tmp
cd tmp

echo "Cloning repository $REPO..."
git clone "https://github.com/$REPO.git"

if [ $? -ne 0 ]; then
    echo "Error cloning repository $REPO."
    exit 1
fi


