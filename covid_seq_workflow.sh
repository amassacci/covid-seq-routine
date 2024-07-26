#!/bin/bash
set -e

#bs auth
export PATH="/home/amassacci/:$PATH"

run=$1    #(e.g. run_060522)
sampleDir="/data/sequence/covid-19/covid-seq/${run}"
mkdir -p $sampleDir
tmp=$sampleDir/tmp
mkdir -p $tmp

#bs list appsessions | grep "DRAGEN COVID Lineage"
id=$(bs list appsessions -f csv | grep "DRAGEN COVID Lineage" | grep $(echo $run | sed -r 's/run_//') | cut -d "," -f 2)

bs appsession download -i $id --extension=.consensus_hard_masked_sequence.fa -o $tmp

find $tmp -type f -not -name '*.json' -exec mv {} $tmp \;

for f in $(ls $tmp/*.consensus_hard_masked_sequence.fa); do
  base=$(basename $f)
  #[[ $base == PositiveControl*.consensus_hard_masked_sequence.fa ]] && continue
  #[[ $base == NTC*.consensus_hard_masked_sequence.fa ]] && continue
  #echo $f
  #echo $base
  cat $f >> $sampleDir/all_sequences.fasta
done

rm -r $tmp/*



#bs list appsessions | grep "DRAGEN RNA Pathogen Detection"
id=$(bs list appsessions -f csv | grep "DRAGEN RNA Pathogen Detection" | grep $(echo $run | sed -r 's/run_//') | cut -d "," -f 2)
coverageDir=$sampleDir/QC
mkdir -p $coverageDir

bs appsession download -i $id --extension=.pathogen_full_res.bed -o $tmp

find $tmp -type f -not -name '*.json' -exec mv {} $tmp \;

echo -e "Sample\tSpike_median_coverage\tSpike_min_coverage\tSpike_max_coverage\tRBD_median_coverage\tRBD_min_coverage\tRBD_max_coverage" > $sampleDir/coverage.tsv
for f in $(ls $tmp/*.pathogen_full_res.bed); do
  base=$(basename $f)
  #[[ $base == PositiveControl*.pathogen_full_res.bed ]] && continue
  #[[ $base == NTC*.pathogen_full_res.bed ]] && continue

  #NC_045512.2:21563-25384
  grep MN908947 $f | awk '{ if ($2 >= 21563 && $2 <= 25384) print }' > $tmp/${base%%.*}_spike_coverage.bed
  #NC_045512.2:22517-23185
  grep MN908947 $f | awk '{ if ($2 >= 22517 && $2 <= 23185) print }' > $tmp/${base%%.*}_rbd_coverage.bed
  echo -e "${base%%.*}\t$(datamash median 4 < $tmp/${base%%.*}_spike_coverage.bed)\t$(sort -n -k4 $tmp/${base%%.*}_spike_coverage.bed | head -1 | cut -f 4)\t$(sort -n -k4 $tmp/${base%%.*}_spike_coverage.bed | tail -1 | cut -f 4)\t$(datamash median 4 < $tmp/${base%%.*}_rbd_coverage.bed)\t$(sort -n -k4 $tmp/${base%%.*}_rbd_coverage.bed | head -1 | cut -f 4)\t$(sort -n -k4 $tmp/${base%%.*}_rbd_coverage.bed | tail -1 | cut -f 4)" >> $sampleDir/coverage.tsv
  #-------coverage locus VoCs
  bedtools intersect -a $f -b /data/sequence/covid-19/covid-seq/voc.bed -nonamecheck > $coverageDir/${base%%.*}_coverage_vocs.tsv
done

rm -r $tmp/*



bs appsession download -i $id --extension=.consensus_hard_masked_sequence.fa -o $tmp

find $tmp -type f -not -name '*.json' -exec mv {} $tmp \;

for f in $(ls $tmp/*.consensus_hard_masked_sequence.fa); do
  base=$(basename $f)
  cat $f >> $sampleDir/all_sequences_rna_pathogen.fasta
done
rm -r $tmp/*


##################################################################################################
Rscript --vanilla /data/sequence/covid-19/covid-seq/hotspot_QC.R $run


source /data/miniconda3/etc/profile.d/conda.sh
conda activate pangolin
pangolin --all-versions
pangolin --update

pangolin --all-versions > $sampleDir/pangolin_version_$(echo $run | sed -r 's/run_//').txt

pangolin $sampleDir/all_sequences.fasta --outdir $sampleDir/ --outfile $sampleDir/${run}_lineage_report.csv --tempdir /scratch/amassacci/testing -t 4 --max-ambig 0.5
pangolin $sampleDir/all_sequences_rna_pathogen.fasta --outdir $sampleDir/ --outfile $sampleDir/${run}_lineage_report_rna_pathogen.csv --tempdir /scratch/amassacci/testing -t 4 --max-ambig 0.5


python3 /data/sequence/covid-19/covid-seq/rename_pre_covid_miner.py ${sampleDir}

gzip ${sampleDir}/all_sequences_rna_pathogen.fa

# Run Covid-miner
#docker run --rm -it -v $sampleDir:/fasta_input -v $sampleDir:/covid-miner/results covid19:latest -f /fasta_input/all_sequences_rna_pathogen.fa.gz --merge no

