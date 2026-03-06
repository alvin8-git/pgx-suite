import pandas as pd
class VCF:
    def __init__(self):
        self.meta = []
        self.header = []
        self.data = []

def read_vcf_simple(file):
    f = gzip.open(file, 'rt') if '.gz' in file else open(file)
    vcf = VCF()
    for line in f:
        if '##' in line:
            vcf.meta.append(line)
            continue
        fields = line.strip().split('\t')
        if fields[0] == '#CHROM':
            vcf.header = fields
            continue
        chr = fields[0].replace('chr', '')
        vcf.data.append([chr] + fields[1:])
    f.close()
    return vcf

def partial(lst, query):
    for i,s in enumerate(lst):
         if query in s:
              return({'sample':i, 'gt':s})
         else:
              pass

#def check_impute_calls(phaseme_vcf, phased_vcf,genome_build):
#   starTable=pd.read_table("grc38_star.txt")
#   check=[]
#   for i,line in enumerate(phaseme_vcf.data):
#      x=partial(line,"./.")
#      if x is not None:
#          x['pos']=line[1]
#          x['ref']=line[3]
#          x['alt']=line[4]
#          check.append(x)
#      
#
#   check2=[]
#   phased_df=pd.DataFrame(phased_vcf.data, columns=phased_vcf.header)
#   for l in check:
#       gt_phased=phased_df.loc[phased_df['POS']==l['pos'],].values.flatten().tolist()[l['sample']]
#       if ("0|0" not in gt_phased):
#            l['sample']=phased_df.columns[l['sample']]
#            l['gt_beagle']=gt_phased
#            check2.append(l)
#
#   check3=[]
#   for l in check2:
##       print(l)
#       tmp=starTable.loc[starTable['pos']==int(l['pos']),].reset_index()
#       #print(tmp)
#       if tmp.shape[0]==0:
#           next
#       elif tmp.shape[0]==1 and tmp.loc[0,'ref']==l['ref'] and tmp.loc[0,'alt']==l['alt']:
#           l['star']=tmp['name'].to_list()
#           check3.append(l)
#       else:
#           tmp1=tmp[['pos','ref','alt']].drop_duplicates()
#           for i,row in tmp1.iterrows():
#               if row['ref']==l['ref'] and row['alt']==l['alt']:
#                  stars=tmp.loc[((tmp['ref']==row['ref']) & (tmp['alt']==row['alt'])),'name'].to_list()
#                  l['star']=stars
#                  check3.append(l)
#   return(check3)


def check_impute_calls(phaseme_vcf, phased_vcf,genome_build):
   starTableName=genome_build+"_star.txt"
   starTable=pd.read_table(starTableName)
   check=[]
   for i,line in enumerate(phaseme_vcf.data):
      x=partial(line,"./.")
      if x is not None:
          x['pos']=line[1]
          check.append(x)      

   check2=[]
   phased_df=pd.DataFrame(phased_vcf.data, columns=phased_vcf.header)
   for l in check:
       tmp=starTable.loc[starTable['pos']==int(l['pos']),].reset_index()
       if tmp.shape[0]==0:
           next
       else:
           gt_phased=phased_df.loc[phased_df['POS']==l['pos'],].values.flatten().tolist()[l['sample']]
           if "0|0" in gt_phased:
               next
           elif tmp,shape[0]==1 andtmp.loc[0,'ref']==l['ref'] and tmp.loc[0,'alt']==l['alt']:
               l['sample']=phased_df.columns[l['sample']]
               l['gt_beagle']=gt_phased
               l['star']=tmp['name'].to_list()
               check2.append(l)
           else:
               tmp1=tmp[['pos','ref','alt']].drop_duplicates()
               for i,row in tmp1.iterrows():
                  if row['ref']==l['ref'] and row['alt']==l['alt']:
                     l['sample']=phased_df.columns[l['sample']]
                     l['gt_beagle']=gt_phased
                     stars=tmp.loc[((tmp['ref']==row['ref']) & (tmp['alt']==row['alt'])),'name'].to_list()
                     l['star']=stars
                     check2.append(l)

   return(check2)

phaseme_vcf=read_vcf_simple("phaseme_1.vcf")
phased_vcf=read_vcf_simple("phased_1.vcf")
pd.DataFrame(check_impute_calls(phaseme_vcf,phased_vcf,genome_build="grc38")).to_csv("check_impute_report.txt",sep="\t",index=False)

