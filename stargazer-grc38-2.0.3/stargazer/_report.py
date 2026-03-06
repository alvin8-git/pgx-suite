import pandas as pd
def create_report(result_df, star_table_path,target_gene,impute_dict):
    star_table=pd.read_table(star_table_path)
    impute_dict_table=pd.DataFrame(impute_dict)
    report=[]
    def get_other_names(hap,genest):
        hap2=genest.loc[genest['name']==hap,"other_names"].to_string(index=False)
        if hap2=="." or "Series" in hap2:
           return hap
        else:
           return hap2
    df_gene=star_table.loc[star_table['gene']==target_gene,]
    for i,row in result_df.iterrows():
         hap1=get_other_names(row['hap1_main'],df_gene)
         hap2=get_other_names(row['hap2_main'],df_gene)
         dip=hap1+"/"+hap2
         imputed=[]
         if len(impute_dict_table)>0:
            tmp=impute_dict_table.loc[impute_dict_table['sample']==row['name'],]
            if tmp.shape[0]>0:
               for i,row in tmp.iterrows():
                  stars=tmp[['star']].tolist()
                  if hap1 in stars:
                      imputed.append(hap1)
                  if hap2 in stars:
                      imputed.append(hap2)
         else:
            imputed.append('none') 
         #print("Diplotype:",dip)
         hap1_cands=row['hap1_cand'].split(",")
         #for y in hap1_cands:
         hapCands_pre=[get_other_names(x,df_gene) for x in hap1_cands]
         #print(hapCands)
         hap2_cands=row['hap2_cand'].split(",")
         hapCands_pre.extend([get_other_names(x,df_gene) for x in hap2_cands])
         #print("before:",row['name'],hap1,hap2,hapCands_pre)
         hapCands_pre1 = list(filter((hap1).__ne__, hapCands_pre))
         #print("before1:",row['name'],hap1,hap2,hapCands_pre1)
         hapCands = list(filter((hap2).__ne__, hapCands_pre1))
         #hapCands.remove(hap1)
         #hapCands.remove(hap2)
         #print("after:",row['name'],hapCands)
         hapCands=list(set(hapCands))
         if target_gene=="g6pd" and len(hapCands)>1:
            mayAlsoBe=";".join(hapCands)
         elif len(hapCands)>1:
            mayAlsoBe=",".join(hapCands)
         elif len(hapCands)==1:
            mayAlsoBe=hapCands[0]
         else:
            mayAlsoBe=""
         #print("May also be:",hapCands)
         dip_cands=row['dip_cand'].split(",")
         dipCands=[get_other_names(x,df_gene) for x in dip_cands]
         #print(dipCands)
         if hap1 in dipCands:
             dipCands.remove(hap1)
         if hap2 in dipCands:
             dipCands.remove(hap2) 
         for y in hapCands:
              #print(y,hapCands)
              #print("up to now:",row['name'],dipCands,mayAlsoBe)
              if y in dipCands:
                   dipCands.remove(y)
         if target_gene=="g6pd" and len(dipCands)>1:
              alsoPossible=";".join(dipCands)
         elif len(dipCands)>1:
              alsoPossible=",".join(dipCands)
         elif len(dipCands)==1:
              alsoPossible=dipCands[0]
         else:
              alsoPossible=""
         #print("Also possible haplotype:",row['name'],dipCands,mayAlsoBe,alsoPossible)
         report.append({"Sample":row['name'],
                   "Gene":target_gene,
                   "Diplotype":dip,
                   "BEAGLE imputed":";".join(imputed),
                   "Phenotype":row['phenotype'], 
                   "Score":row['dip_score'],
                   "May also be":mayAlsoBe,
                   "Also possible haplotype":alsoPossible})
    return(report)
