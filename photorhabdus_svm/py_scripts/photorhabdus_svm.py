import sklearn
from sklearn.svm import SVC
from sklearn.utils import shuffle
import pandas as pd

model = SVC()

neg_train = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testSS/negative_training_ss.pkl')
neg_test = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testSS/negative_testing_ss.pkl')
pos_train = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testSS/positive_training_ss.pkl')
pos_test = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testSS/positive_testing_ss.pkl')

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


model.fit(X_train, y_train)
print(model.score(X_test, y_test))
