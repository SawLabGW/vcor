## process genbank files
for i in `ls gbks/ | sed 's/.gbk//g'`;do anvi-script-process-genbank -i gbks/${i}.gbk -O ${i} ; done
mv *.fa processed_gbks/
mv *.txt processed_gbks/


## create contig dbs and import functions
cd processed_gbks
for i in `ls *.fa | sed 's/-contigs.fa//g'`;do anvi-gen-contigs-database -f ${i}-contigs.fa -n ${i} --external-gene-calls ${i}-external-gene-calls.txt -o ${i}.db; done
for i in `ls *.db`;do anvi-run-hmms -T 4 -c ${i}; done
for i in `ls *.db | sed 's/.db//g'`;do
    anvi-import-functions -c ${i}.db -i ${i}-external-functions.txt
done

## create input file tabs
echo -e "name\tcontigs_db_path" > external_db.list
awk '{print $1 "\t" "processed_gbks/" $1 ".db"}' genomes.list >> external_db.list

## create storage db
anvi-gen-genomes-storage -e external_db.list -o GENOMES.db --gene-caller 'NCBI_PGAP'

anvi-pan-genome -g GENOMES.db -n Vcors_pan

cd Vcors_pan

anvi-display-pan -g ../GENOMES.db -p Vcors_pan-PAN.db

## list the bins manually selected and stored in a colection
anvi-get-sequences-for-gene-clusters --list-bins -g ../GENOMES.db -p Vcors_pan-PAN.db -C bins_found
mkdir gene_clusters_bins

## extract gene clusters only found in each strain
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b OCN008 -o gene_clusters_bins/OCN008_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b RE98 -o gene_clusters_bins/RE98_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b core -o gene_clusters_bins/core_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b H1 -o gene_clusters_bins/H1_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b BAA450 -o gene_clusters_bins/BAA450_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b OCN014 -o gene_clusters_bins/OCN014_clusters.fa
anvi-get-sequences-for-gene-clusters -g GENOMES.db -p Vcors_pan/Vcors_pan-PAN.db -C bins_found -b RE22 -o gene_clusters_bins/RE22_clusters.fa

## get gene cluster annotations
cd gene_clusters_bins
for i in `grep ">" RE98_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_RE98.txt
for i in `grep ">" OCN008_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_OCN008.txt
for i in `grep ">" H1_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_H1.txt
for i in `grep ">" BAA450_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_BAA450.txt
for i in `grep ">" OCN014_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_OCN014.txt
for i in `grep ">" RE22_clusters.fa | cut -d ':' -f2 | cut -d '|' -f1`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz ; done > cluster_annotations_RE22.txt

cat table_header.txt cluster_annotations_RE98.txt > tmp
mv tmp cluster_annotations_RE98.txt
cat table_header.txt cluster_annotations_OCN008.txt > tmp
mv tmp cluster_annotations_OCN008.txt
cat table_header.txt cluster_annotations_H1.txt > tmp
mv tmp cluster_annotations_H1.txt
cat table_header.txt cluster_annotations_BAA450.txt > tmp
mv tmp cluster_annotations_BAA450.txt
cat table_header.txt cluster_annotations_OCN014.txt > tmp
mv tmp cluster_annotations_OCN014.txt
cat table_header.txt cluster_annotations_RE22.txt > tmp
mv tmp cluster_annotations_RE22.txt

## summarize the pangenome
anvi-summarize -p Vcors_pan/Vcors_pan-PAN.db -g GENOMES.db -C bins_found -o summary_pangenome

## calculate ANI values
anvi-compute-genome-similarity --external-genomes external_db.list --program fastANI --output-dir fastANI --num-threads 4 --pan-db Vcors_pan/Vcors_pan-PAN.db

## extract annotation information of gene clusters
python extract_cluster_annotations.py -c cluster_annotations_OCN008.txt -g ../gbks/OCN008.gbk > OCN008_annotations.txt
python extract_cluster_annotations.py -c cluster_annotations_RE98.txt -g ../gbks/RE98.gbk > RE98_annotations.txt
python extract_cluster_annotations.py -c cluster_annotations_H1.txt -g ../gbks/H1.gbk > H1_annotations.txt
python extract_cluster_annotations.py -c cluster_annotations_BAA450.txt -g ../gbks/BAA450.gbk > BAA450_annotations.txt
python extract_cluster_annotations.py -c cluster_annotations_OCN014.txt -g ../gbks/OCN014.gbk > OCN014_annotations.txt
python extract_cluster_annotations.py -c cluster_annotations_RE22.txt -g ../gbks/RE22.gbk > RE22_annotations.txt

## get gene clusters present in all but H1
for i in `cat all_but_H1_clusters.list`;do zgrep ${i} ../Vcors_pan/SUMMARY_bins_found/Vcors_pan_gene_clusters_summary.txt.gz; done > all_but_H1_clusters_annotations.txt
cat table_header.txt all_but_H1_clusters_annotations.txt > tmp
mv tmp all_but_H1_clusters_annotations.txt
for i in RE22 RE98 OCN008 BAA450 OCN014;do
    python extract_cluster_annotations.py -c all_but_H1_clusters_annotations.txt -g ../gbks/${i}.gbk
done
for i in RE22 RE98 OCN008 BAA450 OCN014;do     python extract_cluster_annotations.py -c all_but_H1_clusters_annotations.txt -g ../gbks/${i}.gbk -s ${i}; done > all_but_H1.txt

## run phylogenomic pipeline to build guide tree
anvi-setup-ncbi-cogs --num-threads 4

cd processed_gbks
for i in `ls *.db`;do anvi-migrate --migrate-dbs-safely ${i}; done
for i in `ls *.db`;do anvi-run-ncbi-cogs -c ${i} --num-threads 4; done

anvi-get-sequences-for-hmm-hits --external-genomes external_db.list --hmm-source Bacteria_71 -o concatenated-proteins.fa --return-best-hit --get-aa-sequences --concatenate
grep ">" concatenated-proteins.fa | sed 's/>//g' > ids_concatted.txt
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa -o phylogenomic-tree.txt

head -1 concatenated-proteins.fa | awk '{print $2}' | awk 'BEGIN{FS=","}{for(i=1;i<=NF;i++) print $i}' > genes_list.txt
echo -e 'item_name\tdata_type\tdata_value' > misc_data.txt

## use rooted iqtree output as a guide (from phylophlan)
echo -e 'phylogenomic_tree\tnewick\t(((RE22:9.018449999999999E-4,H1:8.42373E-4):2.872600000000001E-4,OCN014:0.001863164):1.99884E-4,((BAA450:0.001037833,RE98:0.001087706):4.886060000000002E-4,OCN008:0.00201175):5.129799999999997E-5);' >> misc_data.txt

anvi-import-misc-data -p Vcors_pan/Vcors_pan-PAN.db -t layer_orders misc_data.txt
