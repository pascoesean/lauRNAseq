import pandas as pd
import numpy as np
from scipy.linalg import svd
from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder(handle_unknown='ignore', drop='if_binary')
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict, GroupKFold
from sklearn.metrics import mean_squared_error, r2_score

import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("data/processed/vst_seancounts.csv")
data['donor'] = data['donor'].astype('category')
#data['is_sixhour'] = 

# FILTERING : only dmso + jnki broad
data = data[data['drug_broad'] != 'ctrl']
data = data.reset_index(inplace=False).drop('index', axis=1)

X = data.iloc[:, 4:]
y = data[['drug_broad']].astype('category')
enc.fit(y)
y_oh = enc.transform(y).toarray()
donor = data['donor']

kfold_iterator = GroupKFold(n_splits=5)


# downsample jnkis to feature select with lasso and get tunable model

# it's problematic that the model is unable to generalize. need to make generalizable
test_donor = 1520
X_train = X.loc[donor != test_donor, :]
y_train = y_oh[donor != test_donor, :]
X_test = X.loc[donor == test_donor, :]
y_test = y_oh[donor == test_donor, :]

pls = PLSRegression(n_components=2)

pls.fit(X_train, y_train)

prescores = pd.DataFrame(pls.x_scores_, 
                         columns = [f'lv_{i+1}' for i in range(pls.x_scores_.shape[1])],
                         index = X_train.index)

postscores = pd.DataFrame(pls.transform(X_test), 
                         columns = [f'lv_{i+1}' for i in range(pls.x_scores_.shape[1])],
                         index = X_test.index)

scores_df = pd.concat([prescores, postscores], 
                       axis=0).merge(right = data[['timepoint', 'drug_broad', 'donor']], 
                                     left_index=True, right_index=True)


plt.clf()
sns.scatterplot(data=scores_df, x='lv_1', y='lv_2', hue='donor', style='drug_broad')

loadings_df = pd.DataFrame(pls.x_loadings_,
                           columns= [f'lv_{i+1}' for i in range(pls.x_loadings_.shape[1])],
                           index=data.columns[4:])

lv1s = loadings_df.sort_values('lv_1').index    


sns.scatterplot(data=data, x = 'RPL37AP1', y = 'IL1B', hue= 'donor', style = 'timepoint')

# retired?

preds = cross_val_predict(pls, X, y_oh,
                         groups=donor, cv=kfold_iterator)

pred_df = pd.DataFrame(data=preds, columns=['jnki']).reset_index()

sns.barplot(data=pred_df, x = 'index', y = 'jnki')

pls.fit(X_train, y_train)
preds = pls.predict(X_test)


# REDO PLSDA TO NOT INCLUDE CONTROL !! IT IS TOO HARD FOR THE MACHINE TO LEARN
