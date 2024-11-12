#!/bin/bash 
# This is for Docker pipe, default root dir is /atac-seq
# This is for Docker!!!


species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
 

# common tools:
adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"
fastq_dump_tool='fastq-dump.2.10.0'
preseq="preseq"
cutadapt="cutadapt"
fastqc="fastqc"
samtools="samtools"
bwa="bwa"
methylQA="/ref/twlab/software/methylQA"
macs2="macs2"

# genome specific resources:
if [[ $species == mm10 ]]; 
	then
	bwa_ref="/atac_seq/Resource/Genome/mm10/bwa_index_mm10/mm10.fa"
	chrom_size="/atac_seq/Resource/Genome/mm10/mm10.chrom.sizes"
	black_list="/atac_seq/Resource/Genome/mm10/mm10_black_list.bed"
	genome_size=2730871774
	promoter_file="/atac_seq/Resource/Genome/mm10/mm10_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/mm10/mm10_promoter_coding_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == mm9 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/mm9/bwa_index_mm9/mm9.fa"
	chrom_size="/atac_seq/Resource/Genome/mm9/mm9_chrom_sizes"
	black_list="/atac_seq/Resource/Genome/mm9/mm9-blacklist.bed"
	genome_size=2725765481
	promoter_file="/atac_seq/Resource/Genome/mm9/mm9_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/mm9/mm9_coding_promoter_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == hg38 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/hg38/bwa_index_hg38.25/hg38.25_chromsome.fa"
	chrom_size="/atac_seq/Resource/Genome/hg38/hg38.25_chromsome.sizes"
	black_list="/atac_seq/Resource/Genome/hg38/hg38_black_list.bed"
	genome_size=3209286105
	promoter_file="/atac_seq/Resource/Genome/hg38/hg38_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/hg38/hg38_coding_promoter_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == hg19 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/hg19/bwa_index_0.7.5/hg19.fa"
	chrom_size="/atac_seq/Resource/Genome/hg19/hg19_chromosome.size"
	black_list="/atac_seq/Resource/Genome/hg19/hg19_blacklist.bed"
	genome_size=3137161264
	promoter_file="/atac_seq/Resource/Genome/hg19/hg19_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/hg19/hg19_promoter_coding_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == danRer10 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/danRer10/bwa_index_denRer10/danRer10.fa"
	chrom_size="/atac_seq/Resource/Genome/danRer10/danRer10.chrom.sizes"
	touch pesudo_bl.txt
	black_list="pesudo_bl.txt"
	genome_size=1340447187
	promoter_file="/atac_seq/Resource/Genome/danRer10/promoter_region_danRer10_bistream_1k.bed"
	coding_promoter="/atac_seq/Resource/Genome/danRer10/danRer10_coding_promoter_bistream_1k.bed"
	macs2_genome='mm'
elif [[ $species == danRer11 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/danRer11/bwa_index_GTCz11/GRCz11.fa"
	chrom_size="/atac_seq/Resource/Genome/danRer11/GRCz11_chrom.size"
	touch pesudo_bl.txt
        black_list="pesudo_bl.txt"
        genome_size=1345118429
	promoter_file="/atac_seq/Resource/Genome/danRer11/GRCz11_promoter_region.bed"
	coding_promoter="/atac_seq/Resource/Genome/danRer11/pseudo_GRCz11_coding_promoter_region.bed"
	macs2_genome='mm'
elif [[ $species == dm6 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/dm6/bwa_index_dm6/dm6.fa"
	chrom_size="/atac_seq/Resource/Genome/dm6/d.mel.chrom.sizes"
	touch pesudo_bl.txt
        black_list="pesudo_bl.txt"
        genome_size=143726002
	promoter_file="/atac_seq/Resource/Genome/dm6/promoter_region_from_Dmel.bed"
	coding_promoter="/atac_seq/Resource/Genome/dm6/pseudo_coding_promoter_region.bed"
	macs2_genome="dm"
elif [[ $species == rn6 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/rn6/bwa_index_rn6/rn6.fa"
	chrom_size="/atac_seq/Resource/Genome/rn6/rn6.chrom.sizes"
	touch pesudo_bl.txt
        black_list="pesudo_bl.txt"
        genome_size=2870182909
	promoter_file="/atac_seq/Resource/Genome/rn6/promoter_region.bed"
	coding_promoter="/atac_seq/Resource/Genome/rn6/coding_promoter_region.bed"
	macs2_genome="mm"
elif [[ $species == hg38.LOCAL ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/hg38/bwa_index_hg38.25/hg38.25_chromsome.fa"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/hg38/hg38.25_chromsome.sizes"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/hg38/hg38_black_list.bed"
        genome_size=3209286105
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/hg38/hg38_promoter_bistream_1kb.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/hg38/hg38_coding_promoter_bistream_1kb.bed"
elif [[ $species == hg38.NoBlackList ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/hg38_without_blacklist/bwa_index_hg38.25/hg38.25_chromsome.fa"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/hg38_without_blacklist/hg38.25_chromsome.sizes"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/hg38_without_blacklist/hg38_black_list.bed"
        genome_size=3209286105
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/hg38_without_blacklist/hg38_promoter_bistream_1kb.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/hg38_without_blacklist/hg38_coding_promoter_bistream_1kb.bed"
elif [[ $species == personalize ]];
	then
	echo "please add all your preferred file as reference, please make sure that you are very clear of which file is for which"
	echo "remove exit 1 after adding your file"
	exit 1
	bwa_ref=" "
	chrom_size=" "
	black_list=" "
	genome_size=
	promoter_file=" "
	macs2_genome=" "
	coding_promoter=" "
elif [[ $species == HG00621.maternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.maternal/bwa_index_HG00621.maternal/HG00621_mat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.maternal/HG00621.maternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.maternal/HG00621_mat_BlackList.bed"
        genome_size=2998271117
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.maternal/HG00621_mat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.maternal/HG00621_mat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG00741.maternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.maternal/bwa_index_HG00741.maternal/HG00741_mat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.maternal/HG00741.maternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.maternal/HG00741_mat_BlackList.bed"
        genome_size=3026809976
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.maternal/HG00741_mat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.maternal/HG00741_mat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01952.maternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.maternal/bwa_index_HG01952.maternal/HG01952_mat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.maternal/HG01952.maternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.maternal/HG01952_mat_BlackList.bed"
        genome_size=3002018008
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.maternal/HG01952_mat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.maternal/HG01952_mat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01978.maternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.maternal/bwa_index_HG01978.maternal/HG01978_mat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.maternal/HG01978.maternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.maternal/HG01978_mat_BlackList.bed"
        genome_size=3034686331
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.maternal/HG01978_mat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.maternal/HG01978_mat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG03516.maternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.maternal/bwa_index_HG03516.maternal/HG03516_mat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.maternal/HG03516.maternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.maternal/HG03516_mat_BlackList.bed"
        genome_size=3006030321
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.maternal/HG03516_mat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.maternal/HG03516_mat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG00621.paternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal/bwa_index_HG00621.paternal/HG00621_pat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal/HG00621.paternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal/HG00621_pat_BlackList.bed"
        genome_size=2883970837
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal/HG00621_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal/HG00621_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG00741.paternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal/bwa_index_HG00741.paternal/HG00741_pat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal/HG00741.paternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal/HG00741_pat_BlackList.bed"
        genome_size=3012046167
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal/HG00741_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal/HG00741_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01952.paternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal/bwa_index_HG01952.paternal/HG01952_pat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal/HG01952.paternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal/HG01952_pat_BlackList.bed"
        genome_size=2888062891
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal/HG01952_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal/HG01952_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01978.paternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal/bwa_index_HG01978.paternal/HG01978_pat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal/HG01978.paternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal/HG01978_pat_BlackList.bed"
        genome_size=3024711997
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal/HG01978_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal/HG01978_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG03516.paternal ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal/bwa_index_HG03516.paternal/HG03516_pat_ragtag.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal/HG03516.paternal.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal/HG03516_pat_BlackList.bed"
        genome_size=3029039178
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal/HG03516_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal/HG03516_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == CHM13 ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/CHM13/bwa_index_CHM13/chm13v2.0.fa"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/CHM13/CHM13.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/CHM13/CHM13_BlackList.bed"
        genome_size=3117292070
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/CHM13/CHM13_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/CHM13/CHM13_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG00621.paternal.wMito ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal.wMito/bwa_index_HG00621.paternal.wMito/HG00621_pat_ragtag_wMito.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal.wMito/HG00621.paternal_wMito.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal.wMito/HG00621_pat_BlackList.bed"
        genome_size=2883987407
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal.wMito/HG00621_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00621.paternal.wMito/HG00621_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG00741.paternal.wMito ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal.wMito/bwa_index_HG00741.paternal.wMito/HG00741_pat_ragtag_wMito.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal.wMito/HG00741.paternal_wMito.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal.wMito/HG00741_pat_BlackList.bed"
        genome_size=3012062734
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal.wMito/HG00741_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG00741.paternal.wMito/HG00741_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01952.paternal.wMito ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal.wMito/bwa_index_HG01952.paternal.wMito/HG01952_pat_ragtag_wMito.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal.wMito/HG01952.paternal_wMito.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal.wMito/HG01952_pat_BlackList.bed"
        genome_size=2888079471
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal.wMito/HG01952_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01952.paternal.wMito/HG01952_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG01978.paternal.wMito ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal.wMito/bwa_index_HG01978.paternal.wMito/HG01978_pat_ragtag_wMito.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal.wMito/HG01978.paternal_wMito.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal.wMito/HG01978_pat_BlackList.bed"
        genome_size=3024728568
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal.wMito/HG01978_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG01978.paternal.wMito/HG01978_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == HG03516.paternal.wMito ]];
        then
        bwa_ref="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal.wMito/bwa_index_HG03516.paternal.wMito/HG03516_pat_ragtag_wMito.fasta"
        chrom_size="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal.wMito/HG03516.paternal_wMito.ChromSize.txt"
        black_list="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal.wMito/HG03516_pat_BlackList.bed"
        genome_size=3029055749
        promoter_file="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal.wMito/HG03516_pat_Promoter_RustyBAM.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/hllab/Juan/ATACseq/Genomes/HG03516.paternal.wMito/HG03516_pat_CodingPromoter_RustyBAM.bed"
elif [[ $species == hg38.twlab ]];
        then
        bwa_ref="/scratch/twlab/genome/hg38/AIAP_Reference_Files/bwa_index_hg38.25/hg38.25_chromsome.fa"
        chrom_size="/scratch/twlab/genome/hg38/AIAP_Reference_Files/hg38.25_chromsome.sizes"
        black_list="/scratch/twlab/genome/hg38/AIAP_Reference_Files/hg38_black_list.bed"
        genome_size=3209286105
        promoter_file="/scratch/twlab/genome/hg38/AIAP_Reference_Files/hg38_promoter_bistream_1kb.bed"
        macs2_genome="hs"
        coding_promoter="/scratch/twlab/genome/hg38/AIAP_Reference_Files/hg38_coding_promoter_bistream_1kb.bed"
fi




