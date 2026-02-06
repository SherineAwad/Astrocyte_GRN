
astrocyte.h5ad: 
	python src/preprocess_scRNA.py --project_name astrocyte --samples samples.csv


clustered_astrocyte.h5ad: 
	python src/cluster_scRNA.py \
	--project_name astrocyte.h5ad \
	--markers markers.txt \
	--output clustered_astrocyte.h5ad \
	--prefix astrocyte



