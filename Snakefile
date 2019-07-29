import os
import glob
from snakemake.utils import R

workdir: os.getcwd()

REGIONS = ["exon-11-363", "exon-1-2", "exon-2-3", "exon-3-6", "exon-6-11", 'introns', 'promoters', 'exons']
rule all:
    input: "output/supp_fig3/supplfig3.pdf", "input/hg38.chromInfo", \
           "output/supp_fig2/supplfig2.pdf", \
           expand("output/supp_fig1_{peak}/supplfig1_{peak}.pdf", peak=["H3K4me3_ENCFF616DLO", "P300_ENCFF433PKW"]), \
           expand("output/bedtools_fisher/bedtools_fisher_{region}.txt", region=REGIONS), \
           expand("output/binomial_test/binomial_test_{region}.txt", region=REGIONS), \

#--------------------------------------------------------------
# Retrieve a set of peaks (GRCh38, H3K4me3, Homo sapiens K562,
# ENCFF616DLO). Keep only conventional chromosomes
#--------------------------------------------------------------

rule get_h3k4me3:
    output: "input/peaks/H3K4me3_ENCFF616DLO.bed"
    shell: '''
    wget https://www.encodeproject.org/files/ENCFF616DLO/@@download/ENCFF616DLO.bed.gz -O H3K4me3_ENCFF616DLO.bed.gz
    gunzip -f H3K4me3_ENCFF616DLO.bed.gz
    cat H3K4me3_ENCFF616DLO.bed | perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output}
    rm -f H3K4me3_ENCFF616DLO.bed
    '''

#--------------------------------------------------------------
# Retrieve a set of peaks (GRCh38, H3K4me3, Homo sapiens K562,
# ENCFF616DLO). Keep only conventional chromosomes
#--------------------------------------------------------------

rule get_p300:
    output: "input/peaks/P300_ENCFF433PKW.bed"
    shell: '''
    wget https://www.encodeproject.org/files/ENCFF433PKW/@@download/ENCFF433PKW.bed.gz -O P300_ENCFF433PKW.bed.gz
    gunzip -f P300_ENCFF433PKW.bed.gz
    cat P300_ENCFF433PKW.bed | perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output}
    rm -f P300_ENCFF433PKW.bed
    '''

#--------------------------------------------------------------
# Retrieve a GTF from ensembl (Homo sapiens, release 92)
# add 'chr' prefix (-C) and select conventional chromosomes
#--------------------------------------------------------------

rule get_gtf_from_ensembl:
    output: "input/Homo_sapiens_GRCh38_92_chr.gtf"
    shell: '''
    gtftk  retrieve -V 1 -Ccd -r 92 -s homo_sapiens | \
    gtftk select_by_regexp -V 1 -k chrom -r '^chr[0-9XY]+$'  -o {output}
    '''

#--------------------------------------------------------------
# Retrieve chromosome sizes subsequently used by Bedtools
#--------------------------------------------------------------

rule get_chrom_size:
    output: "input/hg38.chromInfo"
    shell: '''
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
       -e "select chrom, size from hg38.chromInfo" | \
       perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output}
    '''

#--------------------------------------------------------------
# Supplementary file 1: OLOGRAM on ChromHMM regions
# defined in K562
# As expected, EP300 is mostly enriched in enhancer-associated
# states, and H3K4me3 in promoter-associated ones.
#--------------------------------------------------------------

def get_chrom_hmm_file(wildcards):
    file_list = glob.glob('input/chromHMM_K562_hg38/*bed')
    return " ".join(file_list)

def get_chrom_hmm_label(wildcards):
    file_list = glob.glob('input/chromHMM_K562_hg38/*bed')
    lab_list = [os.path.splitext(os.path.basename(x))[0] for x in file_list]
    return ",".join(lab_list)

rule supp_fig1:
    input: peak="input/peaks/{peak}.bed"
    params: ref=get_chrom_hmm_file, labels=get_chrom_hmm_label
    output: "output/supp_fig1_{peak}/supplfig1_{peak}.pdf"
    shell: '''
    gtftk ologram -z -y  -c hg38  -p {input.peak} -V 1 --more-bed {params.ref}  \
       --more-bed-labels {params.labels}  -o output/supp_fig1_{wildcards.peak}  \
       -V 3 -k 8 -pf output/supp_fig1_{wildcards.peak}/supplfig1_{wildcards.peak}.pdf
    '''

#--------------------------------------------------------------
# Supplementary file 2:
# We calculate the significance of intersections between H3K4me3
# peaks and GTF-defined numbered exons. The peaks are much more
# present in the first exons, likely due to the broadness of
# H3K4me3 peaks.
#--------------------------------------------------------------

rule supp_fig2:
    input: gtf="input/Homo_sapiens_GRCh38_92_chr.gtf", peak="input/peaks/H3K4me3_ENCFF616DLO.bed"
    output: pdf="output/supp_fig2/supplfig2.pdf", \
            bed1="output/regions_benchmark/ologram_exon-11-363_pygtftk.bed", \
            bed2="output/regions_benchmark/ologram_exon-1-2_pygtftk.bed", \
            bed3="output/regions_benchmark/ologram_exon-2-3_pygtftk.bed", \
            bed4="output/regions_benchmark/ologram_exon-3-6_pygtftk.bed", \
            bed5="output/regions_benchmark/ologram_exon-6-11_pygtftk.bed"
    shell: """
       mkdir -p output/supp_fig2/tmp
       gtftk add_exon_nb -k exon_nbr -i {input.gtf} | \
       gtftk discretize_key -p -d exon_nbr_cat -n 5 -k exon_nbr | \
       gtftk ologram -p {input.peak} -c hg38 -D -n -y -m exon_nbr_cat \
       -pf {output.pdf}  -k 8 -V 3 \
       -j summed_bp_overlaps_true -k 8 -D -o output/supp_fig2 \
       -K output/supp_fig2/tmp
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__11_0_363_0_*.bed {output.bed1}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__1_0_2_0_*.bed {output.bed2}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__2_0_3_0_*.bed {output.bed3}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__3_0_6_0_*.bed {output.bed4}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__6_0_11_0_*.bed {output.bed5}
    """
    
#--------------------------------------------------------------
# Supplementary file 3: OLOGRAM on all genomic features of the GTF
# Only results related to introns are displayed in the finale
# figure. We compare the results for N (number of intersections)
# and S (total number of overlapping nucleotides) for a subset
# of GTF elements. For example, the peaks appear to be significantly
# enriched in introns based on N, but for S that it is not the case ;
# and vice-versa for intergenic regions. Hence S is an important
# statistic to consider : in this particular example it may mean that
# the overlaps of peaks and introns are frequent but short.
# Note that here, an "intersection" means having at least one
# nucleotide in common.
# NB: all tested regions are kept (-K) so that they will be
# subsequently used to test the results of other approaches
#--------------------------------------------------------------

rule supp_fig3:
    input: gtf="input/Homo_sapiens_GRCh38_92_chr.gtf", bed="input/peaks/H3K4me3_ENCFF616DLO.bed"
    output: pdf="output/supp_fig3/supplfig3.pdf", \
            intron='output/regions_benchmark/ologram_introns_pygtftk.bed', \
            exon='output/regions_benchmark/ologram_exons_pygtftk.bed', \
            promoter='output/regions_benchmark/ologram_promoters_pygtftk.bed'
    shell: '''
    mkdir -p output/regions_benchmark
    gtftk ologram -y -V 1 -c hg38 -p {input.bed} -k 8 -o output/supp_fig3 -D \
    -i {input.gtf} -u 1000 -d 1000 -K output/supp_fig3/tmp -pf {output.pdf}
    mv output/supp_fig3/tmp/ologram_introns_pygtftk_*.bed {output.intron}
    mv output/supp_fig3/tmp/ologram_exon_pygtftk_*.bed {output.exon}
    mv output/supp_fig3/tmp/ologram_promoters_pygtftk_*.bed {output.promoter}
    '''

#--------------------------------------------------------------
# Supplementary file 4: The results obtained with bedtools
# fisher when checking the significance of
# overlaps between H3K4me3 and promoter/introns                  
#--------------------------------------------------------------


rule bedtools_fisher:
    # pdf input won't be used but indicates that this rule
    # should be preformed after supplementary figure 3
    input: pdf="output/supp_fig3/supplfig3.pdf", \
           peak="input/peaks/H3K4me3_ENCFF616DLO.bed", \
           chrom="input/hg38.chromInfo", \
           region='output/regions_benchmark/ologram_{region}_pygtftk.bed'
    output: "output/bedtools_fisher/bedtools_fisher_{region}.txt"
    shell: """
    bedtools sort -i {input.region} | bedtools merge > {input.region}.tmp
    bedtools sort -i {input.peak} | bedtools merge > {input.peak}.tmp
    sort {input.chrom} > {input.chrom}.tmp
    time bedtools fisher -b {input.region}.tmp -a {input.peak}.tmp -g {input.chrom}.tmp > {output}
    rm -f {input.region}.tmp {input.peak}.tmp {input.chrom}.tmp
    """

#--------------------------------------------------------------
# Supplementary file 4: The results obtained with a binomial
# test (cf CEAS) when checking the significance of
# overlaps between H3K4me3 and promoter/introns
# NB: Given that CEAS does not provide an annotation
# database for hg38 we implemented this simple solution
# to emulate CEAS approach.
# NB: we will consider here than genome size is of
# hg38 is 2913022398 (see https://tinyurl.com/y6j367hv)
#--------------------------------------------------------------

# here, the number of succes is the number of time H3K4me3 intersect a region
rule compute_nb_success:
    input: pdf3="output/supp_fig3/supplfig3.pdf", \
           pdf2="output/supp_fig2/supplfig2.pdf", \
           peak="input/peaks/H3K4me3_ENCFF616DLO.bed", \
           chrom="input/hg38.chromInfo", \
           region='output/regions_benchmark/ologram_{region}_pygtftk.bed'
    output: "output/bedtools_intersect/bedtools_intersect_{region}.txt"
    shell: """
    bedtools intersect -a {input.peak} -b {input.region} -wa | sort | uniq > {output}
    """
def capitalize_region_name(wildcards):
    return wildcards.region.capitalize()


rule binomial_test:
    input: pdf3="output/supp_fig3/supplfig3.pdf", \
           pdf2="output/supp_fig2/supplfig2.pdf", \
           intersections="output/bedtools_intersect/bedtools_intersect_{region}.txt", \
           peak="input/peaks/H3K4me3_ENCFF616DLO.bed", \
           chrom="input/hg38.chromInfo", \
           region_bed='output/regions_benchmark/ologram_{region}_pygtftk.bed'
    params: region=capitalize_region_name, hg38_size=2913022398
    output: "output/binomial_test/binomial_test_{region}.txt"
    run: R('''
             nb_intersections <- nrow(read.table('{input.intersections}'))
             nb_trials <- nrow(read.table('{input.peak}'))
             region_bed <- read.table('{input.region_bed}', head=F)
             sum_labeled_nuc <- sum(region_bed[,3] - region_bed[,2])
             prob <- sum_labeled_nuc/{params.hg38_size}
             p_val <- pbinom(q=nb_intersections-1, 
                             size=nb_trials, 
                             prob=prob, 
                             lower.tail = FALSE)
             out_df <- data.frame(nb_success=nb_intersections, 
                                  sum_labeled_nuc=sum_labeled_nuc,genome_size={params.hg38_size},
                                  nb_trials=nb_trials,
                                  prob=prob,p_val=p_val, row.names='{params.region}')
             write.table(out_df, file='{output}', col.names=NA, quote=F)
             ''')

    
