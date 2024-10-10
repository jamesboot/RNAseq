# Script to find R1 and R2 for all samples in project and save to CSV along with sample ID
# Input must be directory to search in

# Find R1 and R2 files 
find -L $1 -name "*R1_001.fastq.gz" | sort >> r1.txt
find -L $1 -name "*R2_001.fastq.gz" | sort >> r2.txt

# Create CSV for results to go in
echo "sample,fastq_1,fastq_2" >> samplesheet.csv

# In for loop put all into csv
for i in $(seq 1 $(wc -l < r1.txt));
do
    # Get R1
    R1=$(sed -n "${i}p" r1.txt)
    # Get R2
    R2=$(sed -n "${i}p" r2.txt)
    # Get basename
    BASE=$(basename "$R1")
    # Get ID
    ID="${BASE%%_*}"
    # Write to CSV
    echo "$ID,$R1,$R2" >> samplesheet.csv
done

# Clean up
rm r1.txt
rm r2.txt