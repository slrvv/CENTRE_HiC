#import libraries

import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
import numpy as np # calculate the mean and standard deviation
import xgboost as xgb # XGBoost stuff
from xgboost import plot_importance
from sklearn.model_selection import train_test_split # split  data into training and testing sets
from sklearn.model_selection import GridSearchCV # cross validation
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import confusion_matrix # creates a confusion matrix
#from sklearn.metrics import plot_confusion_matrix # draws a confusion matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from matplotlib import pyplot
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedKFold, cross_val_score
from numpy import sort


#import traing set and process it

xgbost_df = pd.read_csv('/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.txt', header=0, sep=',')
print(xgbost_df.head())
com=xgbost_df.fillna(0)
com['distance'] =com['distance'].abs()
com=com.sort_values('pair19')
com=com.reset_index(drop=True)


#create CV folds for parameter optimization

cv_names=com["CV"].unique()
myCViterator = []
for i in range(len(cv_names)):
    trainIndices = com[ com['CV']!=cv_names[i] ].index.values.astype(int)
    testIndices =  com[ com['CV']==cv_names[i] ].index.values.astype(int)
    myCViterator.append( (trainIndices, testIndices) )


#run randomized search for optimal parameters

X_train= com.drop(['gene_id1','gene_id','symbol38','symbol19','pair','pair19','label','CV'], axis=1).copy()
y_train = com['label'].copy()
model = xgb.XGBClassifier(objective = "binary:logistic",scale_pos_weight=5,random_state=0)
param_grid = {
        'max_depth': [4, 5, 6,8,10,12],
        'learning_rate': [0.1, 0.05, 0.01],
        'gamma': [0, 0.25, 1.0],
        'reg_lambda': [0, 1.0, 10.0],
        'n_estimators': [100,200,300,400,500],
        'colsample_bytree': [0.5,0.6,0.7,0.9],
        'subsample': [0.7, 0.9]
    }
search = RandomizedSearchCV(estimator=model, param_distributions=param_grid,scoring='f1', cv=myCViterator, n_jobs=12, refit=True)
result = search.fit(X_train, y_train)
print('est=%.3f, cfg=%s' % (result.best_score_, result.best_params_))
