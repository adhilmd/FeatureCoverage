# FeatureCoverage
Horizontal and Vertical Coverage Calculation for Whole Exome Sequence data.

This tool(coveragecalc.py) calculates horizontal and vertical coverage for every feature of gtf file such as gene, transcript and exon. Vertical coverage for Gene and Transcript are calculated using the exons in those features by removing introns and overlapping regions. Vertical coverage is calculated for a given bam file also using the index file in the same folder of bam file. Horizontal coverage is calculated for a given interval file, which is used for library prepration or sequencing. Also cytoband and entrez gene ids are annotated for every feature.

# Required python libraries

pandas
re
pybedtools
threading
math
pickle
argparse
pysam
sqlalchemy
sqlite3
sys

#
coveragecalc.py [-h] [-bam <bamfilepath> -ifl <intervalfile> -sn <samplename> -mf <yes|no> -gtf <gtffilepath> -th <multithreads> -ty <grch37|grch38> -pk <picklefilepath> -dir <outputdir>]

---------------------- Horizontal and Vertical Coverage Calculation for Whole
Exome Sequence data (Works only with GRCh GTF files) -------------------- [Author: Mohamood Adhil], [Email:
adhil.md@gmail.com]

optional arguments:
  -h, --help            show this help message and exit
  -bam BAM, --bamfilepath BAM
                        Bam file path
  -ifl IFL, --intervalfile IFL
                        Interval file used for library prepration
  -sn SN, --samplename SN
                        Sample Name
  -gid GEI, --geneidpath GEI
                        gene to geneid map (pickle file)
  -cyt CYT, --cytopath CYT
                        cytoband file path (pickle file)
  -mf {yes,no}, --mergefeature {yes,no}
                        Merge features for exons value: [yes] or [no], if -mf
                        is [yes] a pickle file will be generated on the output
                        directory
  -gtf GTF, --gtffilepath GTF
                        (Required if -mf is [yes]) gtf file path containing 9
                        columns with tab seperated, 3rd column should contain
                        'exon' feature and 9th column should contain three
                        attributes 'gene_name, transcript_id and exon_number'
  -np NP, --numproc NP  Number of cores to use for processing (maxumum 30)
  -ty {hg19,grch37}, --annotype {hg19,grch37}
                        Annotation type value: [hg19] or [grch37]
  -pk PK, --picklefile PK
                        (Required if -mf is [no]) Pickle file contating the
                        merged intervals
  -dir DIR1, --outputpath DIR1
                        Output directory path for storing the results
