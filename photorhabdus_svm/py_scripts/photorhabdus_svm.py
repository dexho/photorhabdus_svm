import sklearn
from sklearn.svm import SVC
from sklearn.utils import shuffle
import pandas as pd

model = SVC()

neg_train = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/negative_training_genes.pkl')
neg_test = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/negative_testing_genes.pkl')
pos_train = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/positive_training_genes.pkl')
pos_test = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/positive_testing_genes.pkl')

frames = [neg_train, pos_train]
df = pd.concat(frames)
df = shuffle(df)
X_train = df.drop(['target'], axis = 'columns')
y_train = df.target

frames = [neg_test, pos_test]
df = pd.concat(frames)
df = shuffle(df)
X_test = df.drop(['target'], axis = 'columns')
y_test = df.target


print(X_train)
print(X_test)
print(y_train)
print(y_test)

model.fit(X_train, y_train)
print(model.score(X_test, y_test))
