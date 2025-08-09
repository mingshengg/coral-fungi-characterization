#repeated for fungi_only_ASVs.fa
cd /home/svu/hpcusers306/db-blast_nt
blastn -query /home/svu/e0321475/Chapter2/unassigned_ASVs.fa -db nt \
	-out /home/svu/e0321475/Chatper2/unassigned_ASVs.out \
	-outfmt "6 qseqid sseqid staxids ssciname scomnames evalue pident length mismatch gapopen" \
	-task megablast \
  	-max_target_seqs 5
