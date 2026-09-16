#!/bin/bash

#Check if all the required files exist
if [[ -f combined.unassembled.forward.fastq && -f combined.unassembled.reverse.fastq && -f combined.discarded.fastq && -f poor_reads.csv && -f good_reads.csv ]]; then
	counts=($(wc -l < combined.unassembled.forward.fastq) $(wc -l < combined.discarded.fastq) $(wc -l < poor_reads.csv) $(wc -l < good_reads.csv))
	total_reads=$((counts[0]/4 + counts[1]/4 + counts[2] + counts[3]))
	percentage=$(awk "BEGIN {printf \"%.1f\", (${counts[3]} / $total_reads) * 100}")
	echo "total_reads,unassembled_reads,discarded_reads,poor_reads,good_reads" > "info.csv"
	echo "$total_reads,$((counts[0]/4)),$((counts[1]/4)),${counts[2]},${counts[3]}" >> "info.csv"
	echo -e "           $(printf "%'d" ${counts[3]}) / $(printf "%'d" $total_reads) reads passing all filters (${percentage}%)"
elif [[ -f poor_reads.csv && -f good_reads.csv && ! -f combined.unassembled.forward.fastq && ! -f combined.unassembled.reverse.fastq && ! -f combined.discarded.fastq ]]; then
	counts=($(wc -l < poor_reads.csv) $(wc -l < good_reads.csv))
	total_reads=$((counts[0] + counts[1]))
	percentage=$(awk "BEGIN {printf \"%.1f\", (${counts[1]} / $total_reads) * 100}")
	echo "total_reads,poor_reads,good_reads" > "info.csv"
	echo "$total_reads,${counts[0]},${counts[1]}" >> "info.csv"
	echo -e "           $(printf "%'d" ${counts[1]}) / $(printf "%'d" $total_reads) reads passing all filters (${percentage}%)"
else
	echo -e "\n$(date '+%I:%M%p') -- Error: One or more required files do not exist"
	exit 1
fi