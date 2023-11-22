import pandas as pd

###E class function
def is_ETopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group E
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in E_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) < 5:
        return ''
    ###check if cat-motif from topology E:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in ['NPPY']:
        #print(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A)
        #if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A:
            return ''
    ### check if sam-motif before cat motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
   
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###F class function
def is_FTopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group F
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in F_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) < 1:
        return ''
    ###check if cat-motif from topology F:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in ['NPPF']:
        #print(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A)
        #if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A:
            return ''
    ### check if sam-motif before cat motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
   
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###G class function
def is_GTopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group G
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in G_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) < 1:
        return ''
    ###check if cat-motif from topology G:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in ['DPPW', 'EPPW']:
        #print(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A)
        #if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A:
            return ''
    ### check if sam-motif before cat motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) > int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
   
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###make table
def make_table(dt1):
    df = pd.DataFrame({'REBASE_name': dt1['REBASE_name'].unique()})
    df['E_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_ETopology(x)).unique())).strip('|'))
    df['F_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_FTopology(x)).unique())).strip('|'))
    df['G_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_GTopology(x)).unique())).strip('|'))
    return df

### input file paths
#/home/masha/RM_project/New/algorithm/rebnr_SMD_pE2_np.regaln.csv
dt_path = input('Path to dataframe:')
dt = pd.read_csv(dt_path, sep = '\t')

#profiles groups
E_profiles = ['Dam']
F_profiles = ['EcoRI_methylase']
G_profiles = ['MT-A70']

#path to saved table
path_save = input('where to save:')

#save table
make_table(dt).to_csv(path_save+'/'+'classes_Pfam.csv')