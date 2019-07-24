import os
import glob

workdir: os.getcwd()

rule all:
    input: "output/supp_fig3/supplfig3.pdf", "input/hg38.chromInfo", \
           expand("output/supp_fig1_{peak}/supplfig1_{peak}.pdf", peak=["H3K4me3_ENCFF616DLO", "P300_ENCFF433PKW"])

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
    gtftk ologram -z  -c hg38  -p {input.peak} -V 1 --more-bed {params.ref}  \
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
    output: "output/supp_fig2/supplfig2.pdf"
    shell: """
    gtftk add_exon_nb -k exon_nbr -i {input.gtf} | \
       gtftk discretize_key -p -d exon_nbr_cat -n 5 -k exon_nbr | \
       gtftk ologram -p {input.peak} -c hg38 -D -n -m exon_nbr_cat \
       -pf output/supp_fig2/supplfig2.pdf  -k 8 -V 3 \
       -j summed_bp_overlaps_true -k 8 -D -o output/supp_fig2
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
    output: "output/supp_fig3/supplfig3.pdf"
    shell: '''
    gtftk ologram -V 1 -c hg38 -p {input.bed} -k 8 -o output/supp_fig3 -D \
    -i {input.gtf} -u 1000 -d 1000 -K output/supp_fig3/tmp -pf output/supp_fig3/supplfig3.pdf
    '''

#--------------------------------------------------------------
# Supplementary file 4: The results obtained with bedtools
# fisher when checking the significance of
# overlaps between H3K4me3 and promoter/introns                  
#--------------------------------------------------------------

def get_region(wildcard):
    file_list = glob.glob('output/supp_fig3/tmp/ologram*pygtftk_*.bed')
    file_list = [x for x in file_list if "_" + wildcards.region + "_" in x]
    return file_list[0]

rule bedtools_fisher:
    # pdf input won't be used but indicates that this rule
    # should be preformed after supplementary figure 3
    input: pdf="output/supp_fig3/supplfig3.pdf", \
           peak="input/peaks/H3K4me3_ENCFF616DLO.bed", \
           chrom="input/hg38.chromInfo"
    params: region=get_region
    output: "output/bedtools_fisher/bedtools_fisher_{region}.txt"
    shell: """
    bedtools fisher -b {params.region} -a {input.peak} -g {input.chrom} > {output}
    """

#--------------------------------------------------------------
# Supplementary file 4: The results obtained with a binomial
# test (cf CEAS) when checking the significance of
# overlaps between H3K4me3 and promoter/introns
# NB: Given that CEAS does not provide an annotation
# database for hg38 we implemented this simple solution
# to emulate CEAS approach.
#--------------------------------------------------------------

rule binomial_test:
    input: pdf="output/supp_fig3/supplfig3.pdf", \
           peak="input/peaks/H3K4me3_ENCFF616DLO.bed", \
           chrom="input/hg38.chromInfo"
        params: region=get_region
    output: "output/binomial_test/binomial_test_{region}.txt"
    shell: R('''
           
           ''')
