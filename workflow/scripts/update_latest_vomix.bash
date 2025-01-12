#!/bin/bash

REPO="holab-hku/vOMIX-MEGA" # Replace with your GitHub repository
PWD=$(pwd)
DIR=$1
DIR_CLEAN="${PWD}/${DIR%/}"


echo "Cloning repository $REPO..."
git clone "https://github.com/$REPO.git" vomix-tmp

if [ $? -ne 0 ]; then
    echo "Error cloning repository $REPO."
    exit 1
fi

echo "[1/5] Updating Snakefile..."
oldfile="${DIR_CLEAN}/workflow/Snakefile"
newfile="vomix-tmp/workflow/Snakefile"
if cmp -s $oldfile $newfile; then : ; else echo "Updating Snakefile"; mv $oldfile $newfile; fi

echo "[2/5] Updating config.yml..."
oldfile="${DIR_CLEAN}/config/config.yml"
newfile="vomix-tmp/config/config.yml"
if [[ -e $oldfile; then cp $oldfile "${original_file%.yml}_old.yml"; else : ;fi
if cmp -s $oldfile $newfile; then : ; else echo "Updating config.yml"; mv $oldfile $newfile; fi

echo "[3/5] Updating rules..." 
for newfile in $(find vomix-tmp/workflow/rules/ -maxdepth 1 -type f); do BASE=$(basename ${newfile}); oldfile="${DIR_CLEAN}/workflow/rules/${BASE}"; if cmp -s "$oldfile" "$newfile"; then : ; else echo "Updating $BASE"; mv $oldfile $newfile; fi; done

echo "[4/5] Updating envs..." 
for newfile in $(find vomix-tmp/workflow/envs/ -maxdepth 1 -type f); do BASE=$(basename ${newfile}); oldfile="${DIR_CLEAN}/workflow/envs/${BASE}"; if cmp -s "$oldfile" "$newfile"; then : ; else echo "Updating $BASE"; mv $oldfile $newfile; fi;  done

echo "[5/5] Updating scripts..." 
for newfile in $(find vomix-tmp/workflow/scripts/ -maxdepth 1 -type f); do BASE=$(basename ${newfile}); oldfile="${DIR_CLEAN}/workflow/scripts/${BASE}"; if cmp -s "$oldfile" "$newfile"; then : ; else echo "Updating $BASE"; mv $oldfile $newfile; fi; done

echo "Done!" 

rm -r 
