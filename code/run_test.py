import time
import os
from prody import parsePDB
import numpy as np
import pandas as pd
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

################################################### DEFINE FUNCTIONS ###########################################################

def get_centroids_df(pdb_dir):
    # parse pdb file and get atom info of all proteins (excludes water and other contents of file)
    atoms = parsePDB(pdb_dir).select('protein')
    
    # get amino acid names and coordinates of atoms
    #df = pd.DataFrame(atoms.getCoords(), columns=['x','y','z'])
    df = pd.DataFrame(atoms.getCoordsets(0), columns=['x','y','z']) # use first model only
    df['acid'] = atoms.getResnames()
    df['chain'] = atoms.getChids()
    df['resnum'] = atoms.getResnums()
    
    # get multiindex df with chain and amino acid resnum as indices
    df2 = df.set_index(['chain', 'resnum'])
    
    # group by amino acids and get centroids (mean coords)
    df3 = df2.groupby(level=['chain','resnum'], sort=False).mean()
    
    # get list original acid names in order without duplicates
    acid_names = [ df[(df['chain']==chain) & (df['resnum']==resnum)]['acid'].values[0] for (chain, resnum) in df3.index ]
    df3['acid'] = acid_names
    
    # final df
    return df3.reset_index()

def find_acid_sequences(df, k=2, theta=6, n=20):
    interaction_list = []

    for i, data_i in df.iterrows():
        chain_i = data_i['chain']
        resnum_i = data_i['resnum']
        chain_i_len = df[df['chain']==chain_i].shape[0]
        #print('\t\tchain_i_len =', chain_i_len)

        for j, data_j in df.iterrows():
            if i == j:
                continue # don't do calculation for the same acid
            else:
                dx = (data_i['x'] - data_j['x'])
                dy = (data_i['y'] - data_j['y'])
                dz = (data_i['z'] - data_j['z'])

                dist = np.sqrt( dx**2 + dy**2 + dz**2 )            
                chain_j = data_j['chain']
                resnum_j = data_j['resnum']


                if dist <= theta: # if 3D distance criteria satsified
                    #print ('\ttheta satisfied')

                    if ((chain_i != chain_j) | ((np.abs(resnum_i - resnum_j) >= n) & (chain_i == chain_j))): # if linear sequence distance criteria satisfied
                        #print ('\tn satisfied')

                        # if match found at begining of sequence, append next k
                        if ((resnum_i-k) < 0):
                            interaction_list.append(df.loc[(df['chain']==chain_i) & (df['resnum'] >= resnum_i) & (df['resnum'] < resnum_i+k), 'acid'].values)

                        # else if match found at end of sequence, append previous k
                        elif ((resnum_i+k) > chain_i_len):
                            interaction_list.append(df.loc[(df['chain']==chain_i) & (df['resnum'] >= resnum_i-k) & (df['resnum'] < resnum_i), 'acid'].values)

                        # else if match found in middle of sequence, append either side of i
                        else:
                            begin_index = int(resnum_i - (k/2))
                            end_index = begin_index + k
                            interaction_list.append(df.loc[(df['chain']==chain_i) & (df['resnum'] >= begin_index) & (df['resnum'] < end_index), 'acid'].values)
                            
    return np.array(interaction_list)

def list_all_file_paths(directory):
    all_files = []
    
    for root,dirs,files in os.walk(directory):
        for f in files:
            all_files.append(os.path.join(root, f))
            
    return all_files

def get_stats(comb): # calculate stats using combined results from all workers
    
    # convert results np array into pandas df
    col_names = list(range(comb.shape[1]-1))
    col_names.append('seq')
    comb_df = pd.DataFrame(comb, columns=col_names)
    
    # calculate amount of unique amino acid sequences involved in interactions
    unique_seqs = comb_df['seq'].unique()
    print('amount of unique amino acid sequences:',unique_seqs.shape[0])
    
    # count amount of times and percent interactions involved for each unique amino acid sequence
    stats_df = comb_df.groupby('seq', sort=False).agg(['count'])[0].sort_values(by=['count'], ascending=False)
    stats_df['pct'] = np.round(100*(stats_df['count'] / unique_seqs.shape[0]), 3)
    
    # count amount of times and percent interactions involved for each unique individual amino acid
    all_acids = []
    for c in list(range(comb.shape[1]-1)):
        all_acids.append(comb_df[c].values)
    all_acids = np.array(all_acids).reshape(-1,1)
    all_acids_df = pd.DataFrame(all_acids, columns=['acid'])
    indiv_acid_df = pd.DataFrame(all_acids_df['acid'].value_counts())
    indiv_acid_df.columns = ['count']
    indiv_acid_df['pct'] = np.round(100*(indiv_acid_df['count'] / indiv_acid_df['count'].sum()), 3)
    
    # save results to .csv files
    indiv_acid_df.to_csv('indiv_acids_{:d}_mers.csv'.format(k))
    print('SAVED: indiv_acids_{:d}_mers.csv'.format(k))
    stats_df.to_csv('acids_seqs_{:d}_mers.csv'.format(k))
    print('SAVED: acids_seqs_{:d}_mers.csv'.format(k))
    
    # plot results
    try:
        import matplotlib.pyplot as plt    
        top = 25

        # top sequnces interacting (percent)
        stats_df['pct'].head(top).plot(kind='bar', figsize=(10,5))
        plt.title('Amino Acid Sequence Involvement in Interactions')
        plt.xlabel('Amino Acid Sequence')
        plt.ylabel('%')
        plt.savefig('pct_involvement_top{:d}_{:d}_mers.png'.format(top, k), bbox_inches='tight', dpi=200)
        print('SAVED: pct_involvement_top{:d}_{:d}_mers.png'.format(top, k))
        plt.clf()

        # top sequnces interacting (count)
        stats_df['count'].head(top).plot(kind='bar', figsize=(10,5))
        plt.title('Amino Acid Sequence Involvement in Interactions')
        plt.xlabel('Amino Acid Sequence')
        plt.ylabel('Count')
        plt.savefig('count_involvement_top{:d}_{:d}_mers.png'.format(top, k), bbox_inches='tight', dpi=200)
        print('SAVED: count_involvement_top{:d}_{:d}_mers.png'.format(top, k))
        plt.clf()

        # top individual interacting (percent)
        indiv_acid_df['pct'].plot(kind='bar', figsize=(10,5))
        plt.title('Individual Amino Acid Involvement in Interactions')
        plt.xlabel('Amino Acid')
        plt.ylabel('%')
        plt.savefig('pct_involvement_indiv_{:d}_mers.png'.format(k), bbox_inches='tight', dpi=200)
        print('SAVED: pct_involvement_indiv_{:d}_mers.png'.format(k))
        plt.clf()

        # top individual interacting (count)
        indiv_acid_df['count'].plot(kind='bar', figsize=(10,5))
        plt.title('Individual Amino Acid Involvement in Interactions')
        plt.xlabel('Amino Acid')
        plt.ylabel('Count')
        plt.savefig('count_involvement_indiv_{:d}_mers.png'.format(k), bbox_inches='tight', dpi=200)
        print('SAVED: count_involvement_indiv_{:d}_mers.png'.format(k))
        plt.clf()
        
    except Exception as e:
        print('***** Stats plots failed *****')
        print('\t',e)

    return np.hstack((np.array(stats_df.index).reshape(-1,1), stats_df.values))

############################################### USE FUNCTIONS ON MASTER AND WORKER NODES ######################################################

### HEAD NODE:
if rank == 0:

    # start overall timer
    start_o = time.time()

    # get all files from input dir
    inp_path = sys.argv[1]
    all_files = list_all_file_paths(inp_path)
    print('files to scatter by head node:', len(all_files))

    # split file names into separate arrays that will be scattered among worker nodes
    split_files = np.array_split(all_files, size)


### WORKER NODE:
else:
    split_files = None


### PROCESSING:

# start timer for worker
start_w = time.time()

# get values for k, theta
k = int(sys.argv[2])
theta = float(sys.argv[3])

# scatter all file names in np array
files = comm.scatter(split_files, root=0)
print('node', rank, 'has files:', files)

# find interacting sequences
interacting_results = []
for f_dir in files:
    try:
        print('node', rank,'analysing file:',f_dir,':')
        print('\tcalculating centroids...')
        df = get_centroids_df(f_dir)

        # find interacting sequences of amino acids
        print('\tfinding interacting sequences...')
        seq = find_acid_sequences(df, k=k, theta=theta, n=20)
        assert seq.shape[1] == k # ensure correct shape of sequences
        
        seq_df = pd.DataFrame(seq)
        combine_acids = lambda row : '-'.join(row)
        seq_df['comb'] = seq_df.apply(combine_acids, axis = 1)        
        
        interacting_results.append(seq_df)
    except Exception as e:
        print('***** something went wrong with file:',f_dir,'on node', rank, '*****')
        print('\t',e)

# combine all results and convert to numpy array for gathering by head node
interacting_results_df = pd.concat(interacting_results, ignore_index=True)
worker_results = interacting_results_df.values

# gather answers from workers (np arrays)
col_data = comm.gather(worker_results, root=0)

# ead node vstack results from workers and calculate statistics
if rank == 0:
    # collect all results from workers, vsack array and print results
    col_data = np.array(col_data)
    col_data = np.vstack(col_data)

    stats = get_stats(col_data)
    
    print()
    print('***************')
    print('*** RESULTS ***')
    print('***************')
    print()
    print('SEQ | COUNT | PCT')
    print()
    print(stats)
    print()

    # stop overall timer
    end_o = time.time()

    # elapsed time overall
    elapsed_time_o = end_o - start_o
    print()
    print('OVERALL ELAPSED TIME:', elapsed_time_o, 'seconds')


# stop timer for worker
end_w = time.time()

# elapsed time for worker
elapsed_time_w = end_w - start_w
print('elapsed time for node ', rank, ':', elapsed_time_w, 'seconds')


############################################################# END ##############################################################
