#!/usr/bin/env bash
# Enable nullglob so the loop doesn't run with the literal pattern if no files match
shopt -s nullglob

if [ $# -eq 0 ]; then
    echo "Usage: $0 <file_pattern>" >&2
    exit 1
fi

# Print header
echo -e "Filename\tQubits\tOrbitals\tMaxMemoryMB\tLowestEnergy0\tLowestEnergy1\tLowestEnergy2\tSolveTime"

for file in "$@"; do
    [ -f "$file" ] || continue
    
    filename=$(basename "$file")
    
    # Extract fields
    qubits=$(grep "^Spin orbitals/qubits:" "$file" | awk '{print $NF}')
    orbitals=$(grep "^Spatial orbitals:" "$file" | awk '{print $NF}')
    memory=$(grep "^PySCF max_memory:" "$file" | awk '{print $(NF-1)}')
    
    # Extract the three lowest energies
    energy0=$(grep -A4 "^Lowest energies" "$file" | awk '$1=="0" {print $2}')
    energy1=$(grep -A4 "^Lowest energies" "$file" | awk '$1=="1" {print $2}')
    energy2=$(grep -A4 "^Lowest energies" "$file" | awk '$1=="2" {print $2}')
    
    # Extract solve time (0.945 is field 3: "Solve" $1, "time:" $2, "0.945" $3)
    solvetime=$(grep "^Solve time:" "$file" | awk '{print $3}')
    
    # Output formatted row
    echo -e "${filename}\t${qubits}\t${orbitals}\t${memory}\t${energy0}\t${energy1}\t${energy2}\t${solvetime}"
done

