# Python code to illustrate 
# classification using data set 
#Importing the required library 

#https://www.geeksforgeeks.org/regression-classification-supervised-machine-learning/

import pandas as pd 
from sklearn.cross_validation import train_test_split 
from sklearn.ensemble import RandomForestClassifier 
from sklearn.preprocessing import LabelEncoder 
from sklearn.metrics import confusion_matrix 
from sklearn.metrics import accuracy_score 
from sklearn.metrics import classification_report 

#Importing the dataset 
dataset = pd.read_csv( 
		'https://archive.ics.uci.edu/ml/machine-learning-'+
		'databases/iris/iris.data',sep= ',', header= None) 
data = dataset.iloc[:, :] 

#checking for null values 
print("Sum of NULL values in each column. ") 
print(data.isnull().sum()) 

#seperating the predicting column from the whole dataset 
X = data.iloc[:, :-1].values 
y = dataset.iloc[:, 4].values 

#Encoding the predicting variable 
labelencoder_y = LabelEncoder() 
y = labelencoder_y.fit_transform(y) 

#Spliting the data into test and train dataset 
X_train, X_test, y_train, y_test = train_test_split( 
			X, y, test_size = 0.3, random_state = 0) 

#Using the random forest classifier for the prediction 
classifier=RandomForestClassifier() 
classifier=classifier.fit(X_train,y_train) 
predicted=classifier.predict(X_test) 

#printing the results 
print ('Confusion Matrix :') 
print(confusion_matrix(y_test, predicted)) 
print ('Accuracy Score :',accuracy_score(y_test, predicted)) 
print ('Report : ') 
print (classification_report(y_test, predicted)) 

