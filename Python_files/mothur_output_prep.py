import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()

#-db DATABSE -u USERNAME -p PASSWORD -size 20
parser.add_argument("-i", "--input_dir", help="Full path to input files")
parser.add_argument("-ct", "--count", help="Count table")
parser.add_argument("-tx", "--tax", help="Taxonomy file (not summary)")
parser.add_argument("-t", "--top", help="Top n OTUs", type=int)

args = parser.parse_args()

count_df = pd.read_csv(args.input_dir+'/'+args.count, sep="\t", header=0, index_col=0) #Read in count file
count_df = count_df.drop(['total'], axis=1) #Delete total column

no_cols = len(count_df.columns)
no_rows = len(count_df.index)

print("\n Input count table has "+str(no_cols)+" samples and "+str(no_rows)+" OTUs. You have specified that you would like to find the top "+str(args.top)+" OTUs per sample")

cols = list(count_df) #get sample column headings
rel_perc_df = count_df[cols].div(count_df[cols].sum(axis=0), axis=1).multiply(100) #Create new df with rel percentages per sample
rel_perc_df.to_csv(args.input_dir+'/'+'Relative_abundance_table.tab', sep='\t', index=True) #Write out to file

tax_df = pd.read_csv(args.input_dir+'/'+args.tax, sep="\t", header=0, index_col=0) #Read in taxonomy info
tax_df[['Kingdom', 'Phylum', 'Class','Order','Family','Genus','Species']] = tax_df['Taxonomy'].str.split(';', expand=True) #Sep out tax data to columns
tax_df = tax_df.drop(['Size','Taxonomy'], axis=1) #Get rid of "Size" and original taxonomy column

#Add tax infor to rel abundance df
rel_perc_df.index.names = ['OTU'] #Change index name to match
taxrel_df = pd.merge(tax_df,rel_perc_df,on='OTU', how='inner') #Merge tables
taxrel_df.to_csv(args.input_dir+'/'+'Relative_abundance_table_with_taxonomy.tab', sep='\t', index=True) #print to file

list_of_topN_OTUs_lists = np.array([rel_perc_df[c].nlargest(args.top).index.values for c in rel_perc_df]) #Get top N OTUs per sample
flat_list = [item for sublist in list_of_topN_OTUs_lists for item in sublist] #Convert to regular list
final_OTU_list = list(set(flat_list)) #Remove duplicates

topN_df = rel_perc_df.loc[final_OTU_list]
topN_df = topN_df.sort_values(by=cols, ascending=False)
topN_df.to_csv(args.input_dir+'/'+'Top10_OTUs_per_sample.tab', sep='\t', index=True) 

taxrel_df_topN = pd.merge(tax_df,topN_df,on='OTU', how='inner') #Merge tables
taxrel_df_topN = taxrel_df_topN.sort_values(by=cols, ascending=False)
taxrel_df_topN.to_csv(args.input_dir+'/'+'Top10_OTUs_per_sample_with_taxonomy.tab', sep='\t', index=True) #print to file

print("Output successfully generated")