

import pandas as pd

###I class function
def is_ITopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) < 5:
        return ''
    ### check if sam-motif before cat motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
   
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]
    
###make table
def make_table(dt1):
    df = pd.DataFrame({'REBASE_name': dt1['REBASE_name'].unique()})
    df['I_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_ITopology(x)).unique())).strip('|'))
    return df

### input file paths
dt_path = input('Path to dataframe:')
dt = pd.read_csv(dt_path, sep = '\t')

#profiles groups
#path to saved table
path_save = input('where to save:')

#save table
make_table(dt).to_csv(path_save+'/'+'classes_I.csv')
