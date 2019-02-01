# FeatureCoverage
Horizontal and Vertical Coverage Calculation for Whole Exome Sequence data.

This tool(coveragecalc.py) calculates horizontal and vertical coverage for every feature of gtf file such as gene, transcript and exon. Vertical coverage for Gene and Transcript are calculated using the exons in those features by removing introns and overlapping regions. Vertical coverage is calculated for a given bam file also using the index file in the same folder of bam file. Horizontal coverage is calculated for a given interval file, which is used for library prepration or sequencing. Also cytoband and entrez gene ids are annotated for every feature.

# Required python libraries

pandas,
re,
pybedtools,
threading,
math,
pickle,
argparse,
pysam,
sqlalchemy,
sqlite3,
sys
