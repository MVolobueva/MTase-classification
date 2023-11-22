import pandas as pd

###A class function
def is_ATopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group A
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in A_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 8:
        return ''
    ###check if cat-motif from topology A:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in catmotif_A:
        if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_A:
            return ''
    ### check if sam-motif before cat motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return 'by motifs'
    
    ### check if Hd1 at the start of sam-domain
    if dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].str.replace(r'\.', '', regex=True).values[0][:6].count('-') > 1:
        return ''
   
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###B class function
def is_BTopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        #print(len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']))
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group B
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in B_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 8:
        return ''
    ###check if cat-motif from topology B:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in catmotif_B:
        if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_B:
            return ''
    ### check if cat-motif before sam-motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) > int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
     
    ### found two sam-subdomain regions
    if len(dt2[dt2['Region_name'] == 'sam_subdom']['Region_coords'].values[0].split(',')) < 2:
        return ''   
    
    
    ### check if S7-Hu3 at the start of sam-domain
    if dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',')[0].count('-') > 10 and dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',')[1].count('-') > 20:
        return ''
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###C class function 
def is_CTopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group C
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in C_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 8:
        return ''
    ###check if cat-motif from topology C:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in catmotif_C:
        if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_C:
            return ''
    ### check if cat-motif before sam-motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
    
    ### found two sam-subdomain regions
    if len(dt2[dt2['Region_name'] == 'sam_subdom']['Region_coords'].values[0].split(',')) < 2:
        return ''
    
    ### check if S7-Hu3 at the start of sam-domain
    if dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',').count('-') > 20 and dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',')[1].count('-') > 10:
        return ''
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###D class function
def is_DTopology(dt2):
    if len(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID']) == 0:
        return ''
    if len(dt2[dt2['Region_name'] == 'sam_motif']['Model_ID']) == 0:
        return ''
    ### check if profile from group D
    if str(dt2[dt2['Region_name'] == 'cat_motif']['Model_ID'].values[0]) not in D_profiles:
        return ''
    ### check if cat-motif == 4 letters
    if  len(dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 4:
        return ''
    ### check if sam-motif == 8 letters
    if  len(dt2[dt2['Region_name'] == 'sam_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0]) != 8:
        return ''
    ###check if cat-motif from topology D:
    if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] not in catmotif_D:
        if  dt2[dt2['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[0] in catmotif_others_D:
            return ''
    ### check if cat-motif before sam-motif:
    if  int(dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0].split('-')[0]) < int(dt2[dt2['Region_name'] == 'sam_motif']['Region_coords'].values[0].split('-')[0]):
        return ''
    
    ### found two sam-subdomain regions
    if len(dt2[dt2['Region_name'] == 'sam_subdom']['Region_coords'].values[0].split(',')) < 2:
        return ''
    
    ### check if S7-Hu3 at the start of sam-domain
    if dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',')[0].count('-') > 20 and dt2[dt2['Region_name'] == 'sam_subdom']['Alignment_frags'].values[0].split(',')[1].count('-') > 10:
        return ''
    
    return dt2[dt2['Region_name'] == 'cat_motif']['Region_coords'].values[0]

###K class function
def is_KTopology(dt3):
    #dt3 = dt3[(dt3['Model_ID'] == 46303)]
    dt3 = dt3[(dt3['Model_ID'] == 46303)]
    dt3 = dt3[~dt3['Region_coords'].isna()]
    dtsam =  dt3[dt3['Region_name'] == 'sam_subdom']
    dtcat =  dt3[dt3['Region_name'] == 'cat_subdom']
    dt4 = dt3[dt3['#:Hit_ID'].isin(list(dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']))]

    
    
    if len(dt4[(dt4['Region_name'] == 'sam_subdom')]) < 2:
        return ''
    if int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[0].split('-')[1]) - int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[0].split('-')[0]) > 30:
        return ''
    if len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))>0:
        for i in range(len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))):
            #print(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i])
            if i != '':
                if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] not in catmotif_C:
                    if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] in catmotif_others_C:    
                        return ''       
    
    return '|'.join(list(dt4[dt4['Region_name'] == 'cat_subdom']['Region_coords'].dropna()))


def is_LTopology(dt3):
    if 'Model_ID' not in dt3:
        return ''
    dt3 = dt3[(dt3['Model_ID'] == 45988)]
    dt3 = dt3[~dt3['Region_coords'].isna()]
    dtsam =  dt3[dt3['Region_name'] == 'sam_subdom']
    dtcat =  dt3[dt3['Region_name'] == 'cat_subdom']
    #dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']
    #filter hitID that contains two part of sam-subdom
    #dt4 = dt3[dt3['#:Hit_ID'].isin(list(dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']))]
    #filter 45988
    dt4 = dt3[dt3['#:Hit_ID'].isin(list(dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']))]
    #print(list(dtsam[dtsam['Region_coords'].str.contains(',')]))
    #print(dt4)
    
    
    if len(dt4[(dt4['Region_name'] == 'sam_subdom')]) < 2:
        #print(dt4)
        return ''
        
    if int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[1].split('-')[0]) - int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[0].split('-')[1]) >= 20:
        return ''
    #print(len(dt4[(dt4['Region_name'] == 'sam_subdom')]))
    ###check if cat-motif from topology B:
    #print(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))
    if len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))>0:
        for i in range(len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))):
            #print(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i])
            if i != '':
                if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] not in catmotif_B:
                    if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] in catmotif_others_B:    
                        return ''       
    if len(dt4[(dt4['Region_name'] == 'sam_subdom')]) == 4:
        if int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[3].split('-')[0]) - int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[2].split('-')[2]) > 20:
            return '|'.join(list(dt4[dt4['Region_name'] == 'cat_subdom']['Region_coords'][:2].dropna()))

    return '|'.join(list(dt4[dt4['Region_name'] == 'cat_motif']['Region_coords'].dropna()))


def is_MTopology(dt3):
    if 'Model_ID' not in dt3:
        return ''
    dt3 = dt3[(dt3['Model_ID'] == 45988)]
    dt3 = dt3[~dt3['Region_coords'].isna()]
    dtsam =  dt3[dt3['Region_name'] == 'sam_subdom']
    dtcat =  dt3[dt3['Region_name'] == 'cat_subdom']
    #dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']
    #filter hitID that contains two part of sam-subdom
    #dt4 = dt3[dt3['#:Hit_ID'].isin(list(dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']))]
    #filter 45988
    dt4 = dt3[dt3['#:Hit_ID'].isin(list(dtsam[~dtsam['Region_coords'].str.contains(',')]['#:Hit_ID']))]
    
    #print(list(dtsam[dtsam['Region_coords'].str.contains(',')]))
    #print(dt4)
    
    
    if len(dt4[(dt4['Region_name'] == 'sam_subdom')]) < 2:
        #print(dt4)
        return ''
        
    if int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[1].split('-')[0]) - int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[0].split('-')[1]) < 20:
        return ''
    #print(len(dt4[(dt4['Region_name'] == 'sam_subdom')]))
    ###check if cat-motif from topology B:
    #print(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))
    if len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))>0:
        for i in range(len(list(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True)))):
            #print(dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i])
            if i != '':
                if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] not in catmotif_B:
                    if  dt4[dt4['Region_name'] == 'cat_motif']['Alignment_frags'].str.replace(r'\W', '', regex=True).values[i] in catmotif_others_B:    
                        return ''       
    if len(dt4[(dt4['Region_name'] == 'sam_subdom')]) == 4:
        if int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[3].split('-')[0]) - int(dt4[(dt4['Region_name'] == 'sam_subdom')]['Region_coords'].values[2].split('-')[2]) > 20:
            return '|'.join(list(dt4[dt4['Region_name'] == 'cat_subdom']['Region_coords'][:2].dropna()))

    return '|'.join(list(dt4[dt4['Region_name'] == 'cat_motif']['Region_coords'].dropna()))
###function for making tables
def make_table(dt1):
    df = pd.DataFrame({'REBASE_name': dt1['REBASE_name'].unique()})
    df['A_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_ATopology(x)).unique())).strip('|'))
    df['B_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_BTopology(x)).unique())).strip('|'))
    df['C_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_CTopology(x)).unique())).strip('|'))
    df['D_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('#:Hit_ID').apply(lambda x: is_DTopology(x)).unique())).strip('|'))
    df['K_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('REBASE_name').apply(lambda x: is_KTopology(x)).unique())).strip('|'))
    df['M_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('REBASE_name').apply(lambda x: is_MTopology(x)).unique())).strip('|'))
    df['L_topology'] = df['REBASE_name'].apply(lambda x: '|'.join(list(dt1[dt1['REBASE_name'] == x].groupby('REBASE_name').apply(lambda x: is_LTopology(x)).unique())).strip('|'))

    return df

### input file paths
#/home/masha/RM_project/New/algorithm/rebnr_SMD_pE2_np.regaln.csv
dt_path = input('Path to dataframe:')
dt = pd.read_csv(dt_path, sep = '\t')

#Load A class motifs
catmotif_A_path = input('Path to cat-motif from class A:')
#'/home/masha/RM_project/New/motifs/A_class.csv'
catmotif_A = list(pd.read_csv(catmotif_A_path)['catMotif'])
#Load B class motifs
catmotif_B_path = input('Path to cat-motif from class B:')
#'/home/masha/RM_project/New/motifs/B_class.csv'
catmotif_B = list(pd.read_csv(catmotif_B_path)['catMotif'])
#Load C class motifs
catmotif_C_path = input('Path to cat-motif from class C:')
#'/home/masha/RM_project/New/motifs/C_class.csv'
catmotif_C = list(pd.read_csv(catmotif_C_path)['catMotif'])
#Load D class motifs
catmotif_D_path = input('Path to cat-motif from class D:')
#'/home/masha/RM_project/New/motifs/D_class.csv'
catmotif_D = list(pd.read_csv(catmotif_D_path)['catMotif'])

#profiles groups
A_profiles = ['51816', '52618', '53087', '54378']
B_profiles = ['45988', '37952', '36976']
C_profiles = ['45633', '46303', '46923']
D_profiles = ['48856', '52484']

#others motifs
catmotif_others_A = catmotif_D+catmotif_B+catmotif_C
catmotif_others_B = catmotif_A+catmotif_D+catmotif_C
catmotif_others_C = catmotif_A+catmotif_B+catmotif_D
catmotif_others_D = catmotif_A+catmotif_B+catmotif_C

#path to saved table
path_save = input('where to save:')

#save table
make_table(dt).to_csv(path_save+'/'+'classes.csv')