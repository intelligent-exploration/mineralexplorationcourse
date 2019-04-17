
import numpy as np 

import matplotlib.pyplot as plt


from sklearn import datasets

 
  
 
#https://medium.com/@martinpella/logistic-regression-from-scratch-in-python-124c5636b8ac


class LogisticRegression:
    def __init__(self, lr, num_iter):
        self.lr = lr
        self.num_iter = num_iter
        self.fit_intercept = True
        self.verbose = True
    
    def __add_intercept(self, X):
        intercept = np.ones((X.shape[0], 1))
        #print(intercept, 'intercept')
        return np.concatenate((intercept, X), axis=1)
    
    def __sigmoid(self, z):
        return 1 / (1 + np.exp(-z))
    def __loss(self, h, y):
        return (-y * np.log(h) - (1 - y) * np.log(1 - h)).mean()
    
    def fit(self, X, y):
        if self.fit_intercept:
            X = self.__add_intercept(X)
        
        # weights initialization
        self.theta = np.zeros(X.shape[1])
        
        for i in range(self.num_iter):
            z = np.dot(X, self.theta)
            h = self.__sigmoid(z)
            gradient = np.dot(X.T, (h - y)) / y.size
            self.theta -= self.lr * gradient
            
            if(self.verbose == True and i % 1000 == 0):
                z = np.dot(X, self.theta)
                h = self.__sigmoid(z)
                print( self.__loss(h, y))
    
    def predict_prob(self, X):
        if self.fit_intercept:
            X = self.__add_intercept(X)
    
        return self.__sigmoid(np.dot(X, self.theta))
    
    def predict(self, X, threshold):
        return self.predict_prob(X) >= threshold


 

 

def main(): 
 

    iris = datasets.load_iris()
    x_train = iris.data[:, :2]

    y_train_raw = iris.target

    print(y_train_raw)


    y_train = (y_train_raw!= 0) * 1

    print (y_train)

    num_iter=30000
    lr = 0.1

    model = LogisticRegression(lr, num_iter ) 
    model.fit(x_train, y_train) 
    preds = model.predict(x_train, 0.3) 


    score = (preds == y_train).mean()


    print(score, ' is fitness score - correct')


    print(preds)


    coefficients = model.theta



    print(coefficients, ' coefficients')






 






if __name__ == '__main__':
    main()

