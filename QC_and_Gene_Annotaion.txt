#===========================================================================================================
# Codes for Quality Control and gene annotation using Compute Canada High-Performing Computing System (HPC)
#===========================================================================================================

# Request compute resources
salloc --time=3:00:00 --mem=16G --cpus-per-task=4

#------------------
# Getting the data
#------------------

# Move to the directory where the data is stored
Pwd 
cd ../..
cd scratch/lukens/

# Copy the SRA files to my directory
cp /scratch/lukens/deg_Seqs/* /scratch/yhezema/deg/sequences

# Copy the Reference Genome (RG) to my directory
cp /scratch/lukens/deg_Genome/* /scratch/yhezema/deg/genome

# Extract the .gz files
cd sequences
gunzip SRR10551660.fastq.gz
gunzip SRR10551661.fastq.gz

#---------------------
# Genome exploration
#---------------------

cd genome
grep ">" yeast_genome.fa
head genomic.gtf
wc -l genomic.gtf  # 38215 genomic.gtf
wc -l yeast_genome.fa #151990 yeast_genome.fa

#-----------------------------------------------------------
# Check the quality control (QC) of Raw Reads using fastp 
#-----------------------------------------------------------

# Create a script for trimming using fastp 
Module load fastp 

# =======Step 1-open a text editor=======
nano trim_reads.sh

# =======Step 2-Write the script and save it=======
# Define the output directory
output_dir="/home/yhezema/scratch/deg/trimmed/"
# Loop through all the .fastq files in the sequences directory
for srr_file in /home/yhezema/scratch/deg/sequences/*.fastq
do
    # Define the output file name by appending ".trimmed" to the input file name
    output_file="${output_dir}$(basename $srr_file .fastq)_trimmed.fastq"
    # Define the HTML and JSON report filenames
    html_report="${output_dir}$(basename $srr_file .fastq)_fastp_report.html"
    json_report="${output_dir}$(basename $srr_file .fastq)_fastp_report.json"
        # Run fastp to trim the reads and generate a QC report (HTML and JSON)
    fastp -i $srr_file \
          -o $output_file \
          --cut_front --cut_tail --length_required 36 \
          --html $html_report \
          --json $json_report \
          --thread 4  # Use 4 threads for faster processing
        # Check if fastp ran successfully
    if [ $? -eq 0 ]; then
        echo "Trimming successful for $srr_file. QC report saved as $html_report"
    else
        echo "Error occurred while trimming $srr_file. Check the log file for details."
    fi
done

# =======Step 3-Make the script executable=======
chmod +x trim_reads.sh

# =======Step 4-Run the script=======
./trim_reads.sh

# Review the trimmed output:

Use Control + Shift + m to toggle the tab key moving focus. Alternatively, use esc then tab to move to the next interactive element on the page.
