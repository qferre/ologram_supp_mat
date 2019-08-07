#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt

# root_dir = './output/supp_fig4/'
# concat_file = root_dir + 'merged_ologram_test.tsv'
# res_file_1 = root_dir + 'S_mean_boxplot.pdf'
# res_file_2 = root_dir + 'S_var_boxplot.pdf'

concat_file = snakemake.input[0]
res_file_1 = snakemake.output[0]
res_file_2 = snakemake.output[1]


# Read result
# Skip even rows except for 0, they are the headers of other files
def logic(index):
    if (index % 2 == 0) and index != 0 : return True
    return False
df = pd.read_csv(concat_file, sep='\t', header = 0, index_col = None,
    skiprows = lambda x: logic(x))

# Iterate over each row and build a result object
result = list()
for i, row in df.iterrows():
    header = row['feature_type'].split('_')
    r = {
        'nb_minibatches':int(header[2]),
        'try_id':int(header[5]),
        'S_mean':row['summed_bp_overlaps_expectation_shuffled'],
        'S_var':row['summed_bp_overlaps_variance_shuffled'],
    }

    result += [r]


# Plot figures
result_df = pd.DataFrame(result)

MINIBATCH_SIZE = 10 # Hardcoded for now, keep it the same in the Snakefile
result_df['Number of shuffles'] = MINIBATCH_SIZE * result_df['nb_minibatches']

meanplot = result_df.boxplot(column=['S_mean'], by='Number of shuffles')
plt.savefig(res_file_1)

result_df.boxplot(column=['S_var'], by='Number of shuffles')
plt.savefig(res_file_2)
