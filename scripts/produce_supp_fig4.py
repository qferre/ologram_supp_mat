#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt

concat_file = snakemake.input[0]
res_file_1= snakemake.output[0]
res_file_2= snakemake.output[1]

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
result_df.boxplot(column=['S_mean'], by='nb_minibatches')
#plt.savefig('./S_mean_boxplot.pdf')
plt.savefig(res_file_1)
result_df.boxplot(column=['S_var'], by='nb_minibatches')
#plt.savefig('./S_var_boxplot.pdf')
plt.savefig(res_file_2)
