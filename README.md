# Scenic+ pipeline using mouse Astrocyte data 


## scRNA preprocessing to get the barcodes 

#### Before and after filtering 

##### Before filtering 
![](figures/violin_QC.png?v=2)

##### After filtering 
![](figures/violin_AfterQC.png?v=2)


#### Samples UMAP
![](figures/umap_astrocyte_Samples.png?v=2)

#### Clustering UMAP

![](figures/umap_astrocyte_Clusters.png?v=3)

##### Dotplot 
![](figures/astrocyte_Dotplot.png?v=3)

#### Marker genes UMAP 

<img src="figures/umap_astrocyte_Aldh1l1.png?v=2" alt="aldh1l1" width="33%"><img src="figures/umap_astrocyte_Aldoc.png?v=2" alt="aldoc" width="33%"><img src="figures/umap_astrocyte_Apoe.png?v=2" alt="apoe" width="33%">

<img src="figures/umap_astrocyte_Aqp4.png?v=2" alt="aqp4" width="33%"><img src="figures/umap_astrocyte_Arx.png?v=2" alt="arx" width="33%"><img src="figures/umap_astrocyte_Ascl1.png?v=2" alt="ascl1" width="33%">

<img src="figures/umap_astrocyte_Bcl11a.png?v=2" alt="bcl11a" width="33%"><img src="figures/umap_astrocyte_Bcl11b.png?v=2" alt="bcl11b" width="33%"><img src="figures/umap_astrocyte_Calb1.png?v=2" alt="calb1" width="33%">

<img src="figures/umap_astrocyte_Calb2.png?v=2" alt="calb2" width="33%"><img src="figures/umap_astrocyte_Cdk1.png?v=2" alt="cdk1" width="33%"><img src="figures/umap_astrocyte_Clu.png?v=2" alt="clu" width="33%">

<img src="figures/umap_astrocyte_Dcx.png?v=2" alt="dcx" width="33%"><img src="figures/umap_astrocyte_Dlx1.png?v=2" alt="dlx1" width="33%"><img src="figures/umap_astrocyte_Dlx2.png?v=2" alt="dlx2" width="33%">

<img src="figures/umap_astrocyte_Dlx5.png?v=2" alt="dlx5" width="33%"><img src="figures/umap_astrocyte_Dlx6.png?v=2" alt="dlx6" width="33%"><img src="figures/umap_astrocyte_Elavl3.png?v=2" alt="elavl3" width="33%">

<img src="figures/umap_astrocyte_Elavl4.png?v=2" alt="elavl4" width="33%"><img src="figures/umap_astrocyte_Foxp2.png?v=2" alt="foxp2" width="33%"><img src="figures/umap_astrocyte_Gad1.png?v=2" alt="gad1" width="33%">

<img src="figures/umap_astrocyte_Gad2.png?v=2" alt="gad2" width="33%"><img src="figures/umap_astrocyte_Glul.png?v=2" alt="glul" width="33%"><img src="figures/umap_astrocyte_Hes5.png?v=2" alt="hes5" width="33%">

<img src="figures/umap_astrocyte_Malat1.png?v=2" alt="malat1" width="33%"><img src="figures/umap_astrocyte_Mbp.png?v=2" alt="mbp" width="33%"><img src="figures/umap_astrocyte_Mpeg1.png?v=2" alt="mpeg1" width="33%">

<img src="figures/umap_astrocyte_Necab1.png?v=2" alt="necab1" width="33%"><img src="figures/umap_astrocyte_Necab2.png?v=2" alt="necab2" width="33%"><img src="figures/umap_astrocyte_Neurog2.png?v=2" alt="neurog2" width="33%">

<img src="figures/umap_astrocyte_Npy.png?v=2" alt="npy" width="33%"><img src="figures/umap_astrocyte_Ntsr2.png?v=2" alt="ntsr2" width="33%"><img src="figures/umap_astrocyte_Olig1.png?v=2" alt="olig1" width="33%">

<img src="figures/umap_astrocyte_Olig2.png?v=2" alt="olig2" width="33%"><img src="figures/umap_astrocyte_Pcp4.png?v=2" alt="pcp4" width="33%"><img src="figures/umap_astrocyte_Pdgfra.png?v=2" alt="pdgfra" width="33%">

<img src="figures/umap_astrocyte_Plp1.png?v=2" alt="plp1" width="33%"><img src="figures/umap_astrocyte_Pou3f2.png?v=2" alt="pou3f2" width="33%"><img src="figures/umap_astrocyte_Prox1.png?v=2" alt="prox1" width="33%">

<img src="figures/umap_astrocyte_Reln.png?v=2" alt="reln" width="33%"><img src="figures/umap_astrocyte_S100b.png?v=2" alt="s100b" width="33%"><img src="figures/umap_astrocyte_Scrt1.png?v=2" alt="scrt1" width="33%">

<img src="figures/umap_astrocyte_Scrt2.png?v=2" alt="scrt2" width="33%"><img src="figures/umap_astrocyte_Slc17a6.png?v=2" alt="slc17a6" width="33%"><img src="figures/umap_astrocyte_Slc1a3.png?v=2" alt="slc1a3" width="33%">

<img src="figures/umap_astrocyte_Slc32a1.png?v=2" alt="slc32a1" width="33%"><img src="figures/umap_astrocyte_Sox10.png?v=2" alt="sox10" width="33%"><img src="figures/umap_astrocyte_Sox2.png?v=2" alt="sox2" width="33%">

<img src="figures/umap_astrocyte_Sox9.png?v=2" alt="sox9" width="33%"><img src="figures/umap_astrocyte_Sst.png?v=2" alt="sst" width="33%"><img src="figures/umap_astrocyte_Tcf7l2.png?v=2" alt="tcf7l2" width="33%">

<img src="figures/umap_astrocyte_Tubb3.png?v=2" alt="tubb3" width="33%">


## 🔍 Fragment Count QC Summary

we used `check.py` for this QC check 

## 📊 Fragment Distribution Analysis

---

### 📁 Sample: `13005-TH2_atac_fragments.tsv.gz`

**Overview**
* **Total fragments:** 143,673,298  
* **Total barcodes:** 507,783  
* **Mean fragments per barcode:** 282.94  
* **Median fragments per barcode:** 2.00  
* **Minimum fragments per barcode:** 1  
* **Maximum fragments per barcode:** 1,625,975  

**Percentiles (fragments per barcode)**
* 10th percentile: 1.0  
* 25th percentile: 1.0  
* 50th percentile (median): 2.0  
* 75th percentile: 6.0  
* 90th percentile: 106.0  
* 95th percentile: 205.0  
* 99th percentile: 6,339.2  

**Cells meeting fragment thresholds**
* ≥10 fragments: **94,766** cells (18.7%)  
* ≥50 fragments: **62,349** cells (12.3%)  
* ≥100 fragments: **52,419** cells (10.3%)  
* ≥500 fragments: **9,424** cells (1.9%)  
* ≥1000 fragments: **6,717** cells (1.3%)  

**Distribution summary**
* Top 10% of cells have **≥106.0** fragments  
* Bottom 10% of cells have **≤1.0** fragments  

---

### 📁 Sample: `13784-TH1_atac_fragments.tsv.gz`

**Overview**
* **Total fragments:** 344,732,502  
* **Total barcodes:** 620,928  
* **Mean fragments per barcode:** 555.19  
* **Median fragments per barcode:** 5.00  
* **Minimum fragments per barcode:** 1  
* **Maximum fragments per barcode:** 645,039  

**Percentiles (fragments per barcode)**
* 10th percentile: 1.0  
* 25th percentile: 2.0  
* 50th percentile (median): 5.0  
* 75th percentile: 12.0  
* 90th percentile: 312.0  
* 95th percentile: 1,501.0  
* 99th percentile: 14,666.3  

**Cells meeting fragment thresholds**
* ≥10 fragments: **189,411** cells (30.5%)  
* ≥50 fragments: **84,823** cells (13.7%)  
* ≥100 fragments: **72,069** cells (11.6%)  
* ≥500 fragments: **58,522** cells (9.4%)  
* ≥1000 fragments: **50,235** cells (8.1%)  

**Distribution summary**
* Top 10% of cells have **≥312.0** fragments  
* Bottom 10% of cells have **≤1.0** fragments  

---

### 📁 Sample: `13784-TH2_atac_fragments.tsv.gz`

**Overview**
* **Total fragments:** 361,488,946  
* **Total barcodes:** 624,185  
* **Mean fragments per barcode:** 579.14  
* **Median fragments per barcode:** 5.00  
* **Minimum fragments per barcode:** 1  
* **Maximum fragments per barcode:** 312,897  

**Percentiles (fragments per barcode)**
* 10th percentile: 1.0  
* 25th percentile: 2.0  
* 50th percentile (median): 5.0  
* 75th percentile: 12.0  
* 90th percentile: 239.0  
* 95th percentile: 1,082.0  
* 99th percentile: 18,388.3  

**Cells meeting fragment thresholds**
* ≥10 fragments: **190,804** cells (30.6%)  
* ≥50 fragments: **83,727** cells (13.4%)  
* ≥100 fragments: **70,674** cells (11.3%)  
* ≥500 fragments: **56,133** cells (9.0%)  
* ≥1000 fragments: **36,244** cells (5.8%)  

**Distribution summary**
* Top 10% of cells have **≥239.0** fragments  
* Bottom 10% of cells have **≤1.0** fragments  

---

