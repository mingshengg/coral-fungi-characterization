# Use conda activate qiime2-amplicon-2024.10
# Directory containing the input files
input_dir="./My_raw_reads"
output_dir="./ITSxpress"

# Loop through all forward read files
for r1 in "$input_dir"/*_trim_1.fastq.gz; do
	# Get the sample name by stripping directory and suffix
	sample=$(basename "$r1" "_trim_1.fastq.gz")

	# Define the reverse read file
	r2="$input_dir/${sample}_trim_2.fastq.gz"

	# Define the output files
	trimmed_r1="$output_dir/${sample}_trimmed_r1.fastq.gz"
	trimmed_r2="$output_dir/${sample}_trimmed_r2.fastq.gz"
	logfile="$output_dir/${sample}_logfile.txt"

	# Run ITSxpress
	itsxpress --fastq "$r1" --fastq2 "$r2" --region ITS2 --taxa All --log "$logfile" --outfile "$trimmed_r1" --outfile2 "$trimmed_r2" --threads 4

	echo "Processed $sample"

done
