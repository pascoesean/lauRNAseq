import pandas as pd
import numpy as np
from scipy.linalg import svd
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("data/hvg_counts_fromDESEQ.csv")
data['donor'] = data['donor'].astype('category')
data['drug'] = data['drug'].astype('category')


genes = data.iloc[:, 1:-3]
metadata = data[['core_id', 'timepoint', 'donor', 'drug']].astype('category')


scaled_genes = scale(genes)

U, s, Vh = svd(scaled_genes)

scores = U[:, 0:(len(s))] @ np.diag(s) # here, scores @ Vh = scaled_palmer bc lossless bc nrows > ncols = npcs
eigenvals = (s**2)

scores_df = pd.DataFrame(scores, columns = [f"PC_{i + 1}" for i in range(len(s))])

# add back scores

scores_metadata = pd.merge(metadata, scores_df, left_index = True, right_index = True)


# visual inspection! 
plt.clf()
sns.scatterplot(data = scores_metadata, x = 'PC_3', y = 'PC_6', hue = 'donor', style = 'drug')
sns.swarmplot(data = scores_metadata, y = 'PC_11', hue = 'drug')

# PC3 maybe separates out 985 and 3013 a little bit? within donor

pc3_gene_loadings = pd.DataFrame(Vh[2, :], index = genes.columns, columns = ['weight'])

pc3_gene_loadings['abs_weight'] = np.abs(pc3_gene_loadings['weight'])


# maybe CMKLR1 is diff?

plt.clf()
sns.swarmplot(data = data, x = 'donor', y = 'SLC12A7', hue = 'drug')



# PCR
