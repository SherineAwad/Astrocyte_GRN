
astrocyte.h5ad: 
	python src/preprocess_scRNA.py --project_name astrocyte --samples samples.csv


clustered_astrocyte.h5ad: 
	python src/cluster_scRNA.py \
	--project_name astrocyte.h5ad \
	--markers markers.txt \
	--output clustered_astrocyte.h5ad \
	--prefix astrocyte



scenicOuts/consensus_peak_calling/pseudobulk_bed_files/: 
	python src/pseudobulk.py  --samples_map samples_map.txt --out_dir scenicOuts --n_cpu 1 --temp_dir /tmp 
scenicOuts/:
	python src/peak_calling.py --samples_map samples_map.txt -o scenicOuts


mm10_tss.bed: 
	wget https://resources.aertslab.org/cistarget/regions/mm10-limited-upstream10000-tss-downstream10000-full-transcript.bed
	mv mm10-limited-upstream10000-tss-downstream10000-full-transcript.bed mm10_tss.bed 

