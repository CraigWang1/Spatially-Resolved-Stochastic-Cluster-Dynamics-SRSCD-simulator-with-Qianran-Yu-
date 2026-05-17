#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ==========================================
# CONFIGURATION
# ==========================================
# 1. Choose your style: "define" (for #define VAR VAL) or "equals" (for VAR = VAL)
STYLE="equals"  # Options: "equals" or "define" for #define MYVAR {val} cases

# 2. File and variable configuration
FILE_TO_EDIT="SCDWrapper.cpp"
VAR_TO_EDIT="probFrenkel"      # The variable/macro name in your code

# 3. Base/source value and directory we copy FROM
BASE_VAL="5.0e-6"
BASE_DIR="retention_823K_backDesorb_explengths_dftDesorbAbsorbEsMaxSurfConc_houHVBindEs_preciseSAVcutoffs_HVsavE0.95eV_hbindEGB0.79eV_DFrenkelPairProb5e-6"

# 4. Target values you want to generate, compile, and submit
TARGET_VALUES=("7e-6" "1e-5" "1.5e-5" "2e-5")

# Ensure the base template directory actually exists before starting
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Base source directory $BASE_DIR does not exist!"
    exit 1
fi

# ==========================================
# LOOP THROUGH TARGET VALUES
# ==========================================
START_DIR=$(pwd)

for NEW_VAL in "${TARGET_VALUES[@]}"; do
    echo "------------------------------------------------"
    echo "Processing value: $NEW_VAL"
    echo "------------------------------------------------"
    
    # Define the new directory name for this specific iteration
    NEW_DIR="retention_823K_backDesorb_explengths_dftDesorbAbsorbEsMaxSurfConc_houHVBindEs_preciseSAVcutoffs_HVsavE0.95eV_hbindEGB0.79eV_DFrenkelPairProb${NEW_VAL}"

    echo "Step 1: Creating $NEW_DIR and copying contents from $BASE_DIR..."
    mkdir -p "$NEW_DIR"
    cp -a "$START_DIR/$BASE_DIR/." "$START_DIR/$NEW_DIR/"

    echo "Step 2: Entering $NEW_DIR..."
    cd "$START_DIR/$NEW_DIR"
    rm {output,joblog}.[0-9]*

    echo "Step 3: Updating $FILE_TO_EDIT (setting $VAR_TO_EDIT to $NEW_VAL)..."
    if [ ! -f "$FILE_TO_EDIT" ]; then
        echo "Error: $FILE_TO_EDIT not found in $NEW_DIR!"
        exit 1
    fi
    
    # DYNAMIC SEARCH AND REPLACE BASED ON STYLE
    if [ "$STYLE" = "define" ]; then
        # Matches: #define VAR_TO_EDIT [spaces] BASE_VAL
        # Replaces with: #define VAR_TO_EDIT NEW_VAL
        # Note: We escape the # symbol just in case, though usually not strictly necessary.
        sed -i "s/\#define[[:space:]]\+${VAR_TO_EDIT}[[:space:]]\+${BASE_VAL}/\#define ${VAR_TO_EDIT} ${NEW_VAL}/g" "$FILE_TO_EDIT"
        
    elif [ "$STYLE" = "equals" ]; then
        # Matches: VAR_TO_EDIT [spaces] = [spaces] BASE_VAL
        # Replaces with: VAR_TO_EDIT = NEW_VAL
        sed -i "s/${VAR_TO_EDIT}[[:space:]]*=[[:space:]]*${BASE_VAL}/${VAR_TO_EDIT} = ${NEW_VAL}/g" "$FILE_TO_EDIT"
    else
        echo "Error: Invalid STYLE selected. Choose 'equals' or 'define'."
        exit 1
    fi

    echo "Step 4: Cleaning and building..."
    make clean
    make

    echo "Step 5: Submitting the job..."
    if [ -f "MPI_SUBMIT.sh" ]; then
        qsub MPI_SUBMIT.sh
        echo "Job successfully submitted for $NEW_VAL!"
    else
        echo "Error: MPI_SUBMIT.sh not found in $NEW_DIR. Compilation succeeded, but job was not submitted."
        exit 1
    fi

    # Return to the starting directory to prepare for the next iteration
    cd "$START_DIR"
    echo "Done with $NEW_VAL."
    echo ""
done

echo "================================================"
echo "All jobs processed successfully!"
echo "================================================"
