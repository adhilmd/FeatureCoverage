# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:14:32 2016

@author: adhil
"""
import pandas as pd
import re
from pybedtools import BedTool as bdt
from multiprocessing import Pool
import math
import pickle
import argparse
import pysam
import sqlalchemy
import sqlite3
import sys
import numpy
import copy

# finding duplicates in the list
def list_duplicates(seq):
  seen = set()
  seen_add = seen.add
  # adds all elements it doesn't know yet to seen and all other to seen_twice
  seen_twice = set( x for x in seq if x in seen or seen_add(x) )
  # turn the set into a list (as requested)
  return list(seen_twice)

# split a list into equal number of parts 
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

# merging overlapping intervals
def mergingoverlapintervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged

# removing fist element of features seperated by colon
def listproc(x):
    return ':'.join(x.split(':')[1:])

# Merging intervals for features
def horcov(genepd):
    feature = 'feature'
    genepd = genepd[['chrom','start','end',feature]].drop_duplicates()
    count = pd.DataFrame(genepd.groupby(feature).size())
    count[feature] = count.index
    count.columns = ['count',feature]
    genepd = genepd.merge(count, left_on=feature, right_on=feature, how='left')
    genepd_single = genepd[genepd['count'] == 1]
    genepd_duplicate = genepd[genepd['count'] > 1]
    geneli = list(set(genepd_duplicate[feature].tolist()))
    genepd_alt = pd.DataFrame(columns=['chrom','start','end',feature])    
    for items1 in geneli:
        genepd_temp = pd.DataFrame(columns=['chrom','start','end',feature])
        temp = genepd_duplicate[genepd_duplicate[feature] == items1][['chrom','start','end',feature]]
        temp = temp.sort(['start'])
        start = temp['start'].tolist()
        end = temp['end'].tolist()
        intervals = [(start[i],end[i]) for i,j in enumerate(start)]
        mergedint = mergingoverlapintervals(intervals)
        rstart = []
        rend = []
        for it1,it2 in mergedint:
            rstart.append(it1)
            rend.append(it2)
        chrom = re.sub(',$','',(items1.split(':')[0]+',')*len(rstart)).split(',')
        feat = re.sub(',$','',(items1+',')*len(rend)).split(',')
        genepd_temp['chrom'] = chrom
        genepd_temp['start'] = rstart
        genepd_temp['end'] = rend
        genepd_temp[feature] = feat
        genepd_alt = pd.concat([genepd_alt,genepd_temp])
    genepd_single = genepd_single.drop('count', 1)
    genepd_final = pd.concat([genepd_single,genepd_alt])
    return (genepd_final)

# remove chrom prefix tag
def chrpref(dframe,colm):
    chromp = ['Chromosome', 'chromosome', 'chrom', 'Chrom', 'chr', 'Chr']    
    for itemsc in chromp: 
        dframe[colm] = dframe[colm].replace(to_replace=itemsc, value='', regex=True)
    return (dframe)

# sorting features with chrom number and start position 
def chrsort(dframe,ch,st,ed,pr):
    dframe = dframe.sort([ch,st,ed])
    dframe.index = range(dframe.shape[0])
    if (pr == 'y'):
        dframe[ch] = 'chr'+dframe[ch].astype(str)
    return (dframe)    

# Processing gtffile for merging intervals
def dataproc(gtffile):
    #Reading files into pandas dataframe
    print ("Reading GTF file")
    gtforig = pd.read_csv(gtffile, sep = "\t", comment = "#", header = None, names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    chrrange = range(1, 23)+['X','Y','MT']
    gtf = gtforig[gtforig['seqname'].isin(chrrange)]
    exon = gtf[gtf['feature'] == 'exon']
    exon = chrpref(exon,'seqname')
    #Extracting exon attributes
    exonattrib = exon['attribute'].tolist()
    exonfeature = []
    for items in exonattrib:
        te = re.findall(r'gene_name "[^"]*"',items)
        te1 = re.findall(r'transcript_id "[^"]*"',items)
        te2 = re.findall(r'exon_number "[^"]*"',items)
        exonfeature.append(te[0].strip('gene_name ').strip('\"')+str('|||')+te1[0].strip('transcript_id ').strip('\"')+str('|||')+te1[0].strip('transcript_id ').strip('\"')+':'+te2[0].strip('exon_number ').strip('\"'))
    tempexon =  pd.DataFrame(data=exonfeature, columns=["reqattrib"])
    exon1 = exon.set_index([tempexon.index.values])
    exon1['reqattrib'] = tempexon
    exon1 = exon1[['seqname','start','end','reqattrib']]
    #Extracting required attributes
    id_data1 = [str(x).split('|||') for x in exon1['reqattrib'].tolist()]
    df_id1 = pd.DataFrame(id_data1, columns=["gene","transcript","exonnumber"], index=exon1.index, dtype=str)
    exon1 = exon1.join(df_id1)
    exon1['seqgene'] = exon1['seqname'].map(str)+':'+exon1['gene']
    exon1['seqtrans'] = exon1['seqname'].map(str)+':'+exon1['gene']+':'+exon1['transcript']
    exon1['seqexonid'] = exon1['seqname'].map(str)+':'+exon1['gene']+':'+exon1['exonnumber']    
    #Extracting hg19 exon attributes
    genes = exon1[['seqname','start','end','seqgene']]
    genes.columns = ['chrom','start','end','feature']
    trans = exon1[['seqname','start','end','seqtrans']]
    trans.columns = ['chrom','start','end','feature']
    exons = exon1[['seqname','start','end','seqexonid']]
    exons.columns = ['chrom','start','end','feature']
    return ([genes, trans, exons])


# removing fist element of features seperated by colon
def chromremove(dframe):
    tro = dframe['feature'].tolist()
    trifeat = [listproc(i) for i in tro]
    dframe['feature']=trifeat
    return (dframe)

# coverage calculation for features
def coverageproc(dfdata):    
    dfdata['intbases'] = dfdata['end']-dfdata['start']
    intoverlap = pd.DataFrame(dfdata.groupby('feature')['intbases'].sum())
    intoverlap['feature']=intoverlap.index
    intoverlap.columns = ['intoverlapsum','feature']
    dfdata = dfdata.merge(intoverlap, left_on='feature', right_on='feature', how='left')    
    coververt = pd.DataFrame(dfdata.groupby('feature')['coverage'].sum())
    coververt['feature']=coververt.index
    coververt.columns = ['coverval','feature']
    dfdata = dfdata.merge(coververt, left_on='feature', right_on='feature', how='left')    
    minimum = pd.DataFrame(dfdata.groupby('feature')['start'].min())
    minimum['feature']=minimum.index
    minimum.columns = ['startmin','feature']
    dfdata = dfdata.merge(minimum, left_on='feature', right_on='feature', how='left')    
    maximum = pd.DataFrame(dfdata.groupby('feature')['end'].max())
    maximum['feature']=maximum.index
    maximum.columns = ['endmax','feature']
    dfdata = dfdata.merge(maximum, left_on='feature', right_on='feature', how='left')
    dfdata['coverval'] = pd.to_numeric(dfdata['coverval'], errors='coerce')
    dfdata['intoverlapsum'] = pd.to_numeric(dfdata['intoverlapsum'], errors='coerce')
    dfdata['coverageavg'] = dfdata['coverval']/dfdata['intoverlapsum']
    dfdata['coverageavg'].fillna(0, inplace=True)
    dfdata = dfdata[['chrom','startmin','endmax','feature','coverageavg']].drop_duplicates() 
    return (dfdata)

# coverage with target calculation
def bedhor(fdata,intervalbed):
    fdata['start'] = fdata['start'].map(int)
    fdata['end'] = fdata['end'].map(int)
    fdata = fdata[['chrom','start','end','feature']]
    beddata = bdt.from_dataframe(fdata)
    beddata = beddata.intersect(intervalbed, wao=True)
    df = pd.read_table(beddata.fn, names=['chrom', 'start', 'end', 'feature', 'chrom2', 'start2', 'end2', 'feature2', 'coverage'])
    df = df[['chrom', 'start', 'end', 'feature', 'coverage']] 
    df['chrom'] = map(str, df['chrom'].tolist())  
    df = coverageproc(df)
    df.columns = ['chrom', 'start', 'end', 'feature', 'Coverage with Reference']
    df['Coverage with Reference'] = df['Coverage with Reference'].round(2)
    return (df)

# individual nucleotide coverage calculation
def bedvert(indrttup):
    featr = indrttup[0]
    liint = indrttup[1]
    sfile = indrttup[2]
    samfile = pysam.AlignmentFile(sfile, "rb" )
    featr.index = range(featr.shape[0])
    colsf = ['Chrom','Start','End','featurex','Depth of Coverage','Interval']+[str(items)+'x' for items in liint]
    f1 = pd.DataFrame(columns = colsf)
    tcoxv = [(str(itn)+'x',0) for itn in liint]
    fc = -1
    indli = featr.index.tolist()
    for items in indli:
        coxv = copy.copy(tcoxv)
        read = 0
        a1 = str(featr.iloc[items,0])
        b1 = featr.iloc[items,1]
        c1 = featr.iloc[items,2]
        d1 = str(featr.iloc[items,3])
        for pileupcolumn in samfile.pileup(a1,b1,c1):
            tr = (pileupcolumn.pos, pileupcolumn.n)
            if (b1<=tr[0]<=c1):
                x = tr[1]
                read=read+x
                for i,j in enumerate(liint):
                    if (x >= j):
                        coxv[i] = (coxv[i][0], coxv[i][1]+1)
                    else:
                        break
        vals = [a1,b1,c1,d1,read,c1-b1+1] + [im[1] for im in coxv]
        fc = fc + 1
        f1.loc[fc] = vals
    return (f1)

# cytoband annotation
def cytoanno(fdata,cytobed):
    fdata['start'] = fdata['start'].map(int)
    fdata['end'] = fdata['end'].map(int)
    fdata = fdata[['chrom','start','end','feature']]
    beddata = bdt.from_dataframe(fdata)
    beddata = beddata.intersect(cytobed, wao=True)
    df = pd.read_table(beddata.fn, names=['chrom', 'start', 'end', 'feature', 'chrom2', 'start2', 'end2', 'cytoband', 'strain','bases'])
    df['chrom'] = df['chrom'].astype(str)
    df['feature'] = df['chrom'] + ':' + df['chrom'] + df['cytoband'] + ':' + df['feature']
    df = df[['chrom', 'start', 'end', 'feature']]
    return (df)

def allcovcal(dfdata,liint):
    minimum = pd.DataFrame(dfdata.groupby('Feature')['Start'].min())
    minimum['Feature']=minimum.index
    maximum = pd.DataFrame(dfdata.groupby('Feature')['End'].max())
    maximum['Feature']=maximum.index
    chrm = dfdata[['Chrom','Feature']].drop_duplicates()
    allsum = pd.DataFrame(dfdata.groupby('Feature')[['Depth of Coverage', 'Interval']+[str(items)+'x' for items in liint]].sum())
    allsum = allsum.div(allsum['Interval'], axis = 0)
    allsum['Feature']=allsum.index
    finvt = chrm.merge(minimum, left_on='Feature', right_on='Feature', how='left')
    finvt = finvt.merge(maximum, left_on='Feature', right_on='Feature', how='left')
    finvt = finvt.merge(allsum, left_on='Feature', right_on='Feature', how='left')
    finvt['Depth of Coverage'] = finvt['Depth of Coverage'].round(2)
    itle = [str(items)+'x' for items in liint]
    for items4 in itle:
	finvt[items4] = finvt[items4].round(2)
    return (finvt)

# main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = '''\n\ncoveragecalc.py''', description='''------------- Horizontal and Vertical Coverage Calculation for Whole Exome Sequence data (Works only with GRCh GTF files)-----------		
    [Author: Mohamood Adhil], [Email: adhil.md@gmail.com], [Date: 4th August 2016], [help: python converagecalc.py -h]''', 
    usage = 'coveragecalc.py [-h] [-bam <bamfilepath> -ifl <intervalfile> -sn <samplename> -mf <yes|no> -gtf <gtffilepath> -th <multithreads> -ty <grch37|grch38> -pk <picklefilepath> -dir <outputdir>]')
    parser.add_argument('-bam','--bamfilepath', type=str, dest='bam', help="Bam file path", action = 'store', required = True)
    parser.add_argument('-ifl', '--intervalfile', type=str, dest='ifl', help= "Interval file used for library prepration", action = 'store', required = True)
    parser.add_argument('-sn', '--samplename', type=str, dest='sn', help= "Sample Name", action = 'store', required = True)    
    parser.add_argument('-gid','--geneidpath', type=str, dest='gei', help="gene to geneid map (pickle file)", action = 'store', required = True)
    parser.add_argument('-cyt','--cytopath', type=str, dest='cyt', help="cytoband file path (pickle file)", action = 'store', required = True)
    parser.add_argument('-mf', '--mergefeature', dest='mf', choices=['yes', 'no'], help= "Merge features for exons value: [yes] or [no], if -mf is [yes] a pickle file will be generated on the output directory", action = 'store', required = True)        
    parser.add_argument('-gtf', '--gtffilepath', type=str, dest='gtf', help= "(Required if -mf is [yes]) gtf file path containing 9 columns with tab seperated, 3rd column should contain 'exon' feature and 9th column should contain three attributes 'gene_name, transcript_id and exon_number'", action = 'store')
    parser.add_argument('-gname', '--gtfname', type=str, dest='gpickle', help= "(Required if -mf is [yes]) prefix name for gtf pickle file", action = 'store')
    parser.add_argument('-np', '--numproc', type=int, dest='np', help= "Number of cores to use for processing (maxumum 30)", action = 'store', required = True)
    parser.add_argument('-ty', '--annotype', dest='ty', choices=['hg19', 'grch37', 'grch38'], help= "Annotation type value: [hg19] or [grch37] or [grch38]", action = 'store', required = True)    
    parser.add_argument('-pk', '--picklefile', type=str, dest='pk', help= "(Required if -mf is [no]) Pickle file contating the merged intervals", action = 'store')
    parser.add_argument('-dir', '--outputpath', type=str, dest='dir1', help= "Output directory path for storing the results", action = 'store', required = True)
    args = parser.parse_args()
    sfile = args.bam    
    intfile = args.ifl    
    sn = args.sn    
    mf = args.mf
    numproc = args.np    
    ty = args.ty    
    odir = args.dir1
    gei = args.gei
    print ("Reading GeneID mapping file")
    with open(gei) as fe:
	geid = pickle.load(fe)
    print ("Reading Cytoband file")    
    cyto = args.cyt
    with open(cyto) as fe1:
	cytof = pickle.load(fe1)
    if (mf == 'yes'):
        required_together = ['gtf']
        con = [getattr(args,x) for x in required_together]
        if None in con:
            print("Please provide gtf and multithread")
            quit()
        gtffile = args.gtf
        annot = dataproc(gtffile)
        exons = annot[2]    
        exons = chromremove(exons)
        exons['feature'] = [':'.join(x.split(':')[:3]) for x in exons['feature'].tolist()]
        del annot[2]
        collapseall = pd.concat(annot)
        threads=[]
        finlist=[]
        drt = []
        featseq = list(set(collapseall['feature'].tolist()))
        number = int(math.ceil(len(featseq)/float(numproc)))
        splititems = list(chunks(featseq, number))
        dattemp = [collapseall[collapseall['feature'].isin(its)] for its in splititems]
        print ("Processing whole gtf file, this might take some time")
        p2 = Pool()
        finlist = p2.map(horcov, dattemp)
        #status = parallelprocessintervals(dattemp)
        allge = pd.concat(finlist) 
        allge = chromremove(allge)
        allge = pd.concat([allge,exons])
        allge.index = range(allge.shape[0])
        allge['start'] = allge['start'].map(int)
	allge['end'] = allge['end'].map(int)
        fpt = str(odir)+'/'+args.ty+'.pickle'
        with open(fpt, 'w') as f:
            pickle.dump(allge, f)
    if (mf == 'no'):
        required_together = ['pk']
        con = [getattr(args,x) for x in required_together]
        if None in con:
            print("Please provide pickle file")
            quit()
        fpt = args.pk
        print ("Reading pickle file")
        with open(fpt) as f:
            allge = pickle.load(f)
    allge['chrom'] = map(str, allge['chrom'].tolist())
    allge = chrsort(allge,'chrom','start','end','n')
    allge = allge[allge['chrom'] != 'MT']
    allge = allge[allge['chrom'] != 'plasmid']
    allge = allge.drop_duplicates()
    cytof = chrpref(cytof,'chrom')
    cytof = chrsort(cytof,'chrom','start','end','n')
    cytobed = bdt.from_dataframe(cytof)
    allge = cytoanno(allge, cytobed)
    allgebed = bdt.from_dataframe(allge)
    print ("Reading Interval file")
    intervalpd = pd.read_csv(intfile, sep = "\t", comment = "#", header = None, usecols=[0,1,2,3], names = ['chrom', 'start', 'end', 'feature'])
    intervalpd['chrom'] = map(str, intervalpd['chrom'].tolist())
    intervalpd = chrpref(intervalpd,'chrom')
    intervalpd = chrsort(intervalpd,'chrom','start','end','n')
    intervalpd['start'] = intervalpd['start'].map(int)
    intervalpd['end'] = intervalpd['end'].map(int)
    intervalbed = bdt.from_dataframe(intervalpd)
    allintbed = allgebed.intersect(intervalbed, wo=True)
    allt = pd.read_table(allintbed.fn, names=['chrom', 'startx', 'endx', 'feature', 'chromx', 'start', 'end', 'featurex','overlap'])
    allt = allt[allt['overlap'] >= 1]
    allht = allt[['chrom','startx','endx','feature']]
    allht.columns = ['chrom','start','end','feature']
    allt = allt[['chrom','start','end','feature']]
    intervalte = intervalpd[['chrom','start','end']].drop_duplicates()
    intefeat = list(numpy.arange(1, intervalte.shape[0]+1, 1))
    intefeat = ['Interval-'+str(i) for i in intefeat]
    intervalte['feature'] = intefeat
    intervalanno = cytoanno(intervalte, cytobed)
    allint = pd.concat([allt, intervalanno], axis=0)
    allint_temp = allint[['chrom','start','end']].drop_duplicates()
    inteft = list(numpy.arange(1, allint_temp.shape[0]+1, 1))
    inteft = ['comb-'+str(i) for i in inteft]
    allint_temp['featurex'] = inteft
    allintbed = bdt.from_dataframe(allint)
    allint_tempbed = bdt.from_dataframe(allint_temp)
    allintbedint = allintbed.intersect(allint_tempbed, wo=True,f=1,r=True)
    allint1 = pd.read_table(allintbedint.fn, names=['chrom', 'start', 'end', 'feature', 'chromx', 'startx', 'endx', 'featurex','overlap'])
    allint1 = allint1[['chrom', 'start', 'end', 'feature','featurex']]
    if numproc >= 30:
        numproc = 30
    featge = list(set(allint_temp['featurex'].tolist()))
    number1 = int(math.ceil(len(featge)/float(numproc)))
    splititems1 = list(chunks(featge, number1))
    indrt = [chrsort(allint_temp[allint_temp['featurex'].isin(its)],'chrom','start','end','n') for its in splititems1]
    liint = list(numpy.arange(10, 110, 10))
    indrtfin = [(i,liint,sfile) for i in indrt]
    print ("Calculating Vertical Coverage")
    p1 = Pool()
    drt = p1.map(bedvert, indrtfin)
    dfvert = pd.concat(drt)
    dfvert.index = range(dfvert.shape[0])
    allint1 = allint1.merge(dfvert, left_on='featurex', right_on='featurex', how='left')    
    allint1 = allint1.drop('featurex', 1)         
    allint1 = allint1.drop_duplicates()
    allint1=allint1.rename(columns = {'feature':'Feature'})
    finvt = allcovcal(allint1,liint)
    print ("Calculating Horizontal Coverage")
    dfhor = bedhor(allht,intervalbed)
    dfhor.index = range(dfhor.shape[0])
    intervalanno['Coverage with Reference'] = ['--']*intervalanno.shape[0]
    dfhor = pd.concat([dfhor,intervalanno],axis=0)
    print ("Merging Coverage")
    dfall = dfhor.merge(finvt, left_on='feature', right_on='Feature', how='left')
    dfall = dfall[pd.notnull(dfall['Chrom'])]
    dfall.index = range(dfall.shape[0])
    dfall = dfall[['chrom','start','end','feature','Coverage with Reference','Depth of Coverage']+[str(items)+'x' for items in liint]]
    dfall = chromremove(dfall)
    dftemp = pd.DataFrame(dfall.feature.str.split(':').tolist())
    dftemp.columns = ['Cytoband','Gene/Interval','Transcript','Exon']
    dfall = pd.concat([dfall,dftemp],axis=1)
    dfall.columns = ['Chrom','Start','End','feature','Coverage with Reference', 'Depth of Coverage']+[str(items)+'x' for items in liint]+['Cytoband','Gene/Interval', 'Transcript', 'Exon']
    dfall = dfall.merge(geid, left_on='Gene/Interval', right_on='Gene', how='left')
    dfall['GeneID'] = dfall['GeneID'].fillna(0)
    dfall['GeneID'] = dfall['GeneID'].astype(int)
    dfall['GeneID'] = dfall['GeneID'].replace(to_replace=0, value='--', regex=False)
    dfall['Transcript'] = dfall['Transcript'].fillna('--')
    dfall['Exon'] = dfall['Exon'].fillna('--')
    dfall = dfall[['Chrom','Start','End','Gene/Interval','GeneID','Transcript','Exon','Cytoband','Coverage with Reference', 'Depth of Coverage']+[str(items)+'x' for items in liint]]
    dfall = chrsort(dfall,'Chrom','Start','End','n')
    outfile = str(odir) + '/' +str(sn) +'_'+str(ty)+ '.csv'
    outdb = str(odir) + '/' +str(sn) +'_'+str(ty)+ '.db'
    tablename = str(sn)+'_'+str(ty)
    print ("Writting Coverage Outputfile .csv and .db")
    dfall.to_csv(outfile, sep = '\t', header=True, index=False)
    engine = sqlalchemy.create_engine('sqlite:///' + outdb)
    lic = ['col'+str(i) for i in range(0,len(dfall.columns.tolist()))]
    dfall.columns = lic
    dfall.to_sql(tablename,engine,if_exists='replace',flavor='sqlite',index=False)
