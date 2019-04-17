
import numpy as np 

import matplotlib.pyplot as plt


from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


import random

#http://archive.ics.uci.edu/ml/machine-learning-databases/housing/

#http://archive.ics.uci.edu/ml/datasets/iris
#https://en.wikipedia.org/wiki/Iris_flower_data_set

#https://scikit-learn.org/stable/auto_examples/linear_model/plot_ols.html
#https://www.cs.toronto.edu/~frossard/post/linear_regression/



def get_data():

    #house_data = datasets.load_boston() #Scikit-learn provides a handy description of the dataset, and it can be easily viewed by:
    #print (data.DESCR) 
    #print (data)

    #dataset = np.loadtxt('raw_data/iris.data') # gives error 

    #dataset = np.genfromtxt('processed_data/iris_train.csv',delimiter=',')   # when you have .csv  comma sepeated data  in file
    
    

    #print(dataset, ' iris_data')


    # Load the diabetes dataset
    diabetes = datasets.load_diabetes() 

    data_input = diabetes.data[:, np.newaxis, 2] 

    x_train = data_input[:-20]
    x_test = data_input[-20:]

    # Split the targets into training/testing sets
    y_train = diabetes.target[:-20]
    y_test = diabetes.target[-20:]

    # Split the data into training/testing sets

    return x_train, x_test, y_train, y_test


def sigmoid(z):
    return 1. / (1 + np.exp(-z))


def z(theta, x):
    assert theta.shape[1] == 1
    assert theta.shape[0] == x.shape[1]  # Theta should have as many rows as x has features.
    return np.dot(x, theta)


def hypothesis(theta, x):
    return sigmoid(z(theta, x))


def cost(theta, x, y):
    assert x.shape[1] == theta.shape[0]  # x has a column for each feature, theta has a row for each feature.
    assert x.shape[0] == y.shape[0]  # One row per sample.
    assert y.shape[1] == 1
    assert theta.shape[1] == 1
    h = hypothesis(theta, x)
    one_case = np.matmul(-y.T, np.log(h))
    zero_case = np.matmul(-(1 - y).T, np.log(1 - h))
    return (one_case + zero_case) / len(x)


def gradient_descent(theta, x, y, learning_rate, regularization = 0):
    regularization = theta * regularization
    error = hypothesis(theta, x) - y
    n = (learning_rate / len(x)) * (np.matmul(x.T, error) + regularization)
    return theta - n


def minimize(theta, x, y, iterations, learning_rate, regularization = 0):
    costs = []
    for _ in range(iterations):
        theta = gradient_descent(theta, x, y, learning_rate, regularization)
        costs.append(cost(theta, x, y)[0][0])
    return theta, costs


def main(): 

    x_train, x_test, y_train, y_test = get_data()





    '''mushroom_data = pd.read_csv("../input/mushrooms.csv").dropna()
    mushroom_x = pd.get_dummies(mushroom_data.drop('class', axis=1))



    mushroom_x['bias'] = 1
    mushroom_x = mushroom_x.values
    mushroom_y = (np.atleast_2d(mushroom_data['class']).T == 'p').astype(int)

    x_train, x_test, y_train, y_test = train_test_split(mushroom_x, mushroom_y, train_size=0.85, test_size=0.15)'''
    
    candidate = np.atleast_2d([ np.random.uniform(-1, 1, 118) ]).T

    iterations = 1200
    learning_rate = 0.1
    regularization = 0.5

    theta, costs = minimize(candidate, x_train, y_train, iterations, learning_rate, regularization)
    plt.plot(range(len(costs)), costs)
    plt.show()
    print(costs[-1])

    predictions = x_test.dot(theta) > 0
    len(list(filter(lambda x: x[0] == x[1], np.dstack((predictions, y_test))[:,0]))) / len(predictions)


     






if __name__ == '__main__':
    main()

