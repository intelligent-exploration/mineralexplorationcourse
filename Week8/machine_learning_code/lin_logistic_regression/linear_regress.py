#import statsmodels.api as sm 

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

#https://en.wikipedia.org/wiki/Dot_product





def get_data():

    
    #house_data = datasets.load_boston() #Scikit-learn provides a handy description of the dataset, and it can be easily viewed by:
    #print (data.DESCR) 
    #print (data)

    #dataset = np.loadtxt('raw_data/iris.data') # gives error 

    dataset = np.genfromtxt('processed_data/iris_train.csv',delimiter=',')   # when you have .csv  comma sepeated data  in file
    
    

    print(dataset, ' iris_data')


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

def generate_syntheticdata():

    data_x = np.linspace(1.0, 10.0, 100)[:, np.newaxis]
    #print(data_x, ' ** ')
    data_y = np.sin(data_x) + 0.1 * np.power(data_x, 2) + 0.5 * np.random.randn(100, 1)
    #print(data_y, ' **** ')
    
    data_x /= np.max(data_x) 
    data_x = np.hstack((np.ones_like(data_x), data_x))

    order = np.random.permutation(len(data_x))
    portion = 20
    x_test = data_x[order[:portion]]
    y_test = data_y[order[:portion]]
    x_train = data_x[order[portion:]]
    y_train = data_y[order[portion:]]

    return x_train, x_test, y_train, y_test



def get_gradient(w, x, y):
    y_estimate = x.dot(w).flatten()
    error = (y.flatten() - y_estimate)
    gradient = -(1.0/len(x)) * error.dot(x)
    return gradient, np.power(error, 2)


def numpy_linear_mod(x_train, x_test, y_train, y_test):


    print(' running numpy linear model')

    w = np.random.randn(2)
    alpha = 0.5
    tolerance = 1e-5

    # Perform Gradient Descent
    iterations = 1
    while True:
        gradient, error = get_gradient(w, x_train, y_train)

        
        new_w = w - alpha * gradient
    
        # Stopping Condition
        if np.sum(abs(new_w - w)) < tolerance:
            print("Converged.")
            break
    
        # Print error every 50 iterations
        if iterations % 5 == 0:
            mean_squared_error =  np.sum(error)/error.shape
            print(iterations, mean_squared_error)
    
        iterations += 1
        w = new_w









    
def scipy_linear_mod(x_train, x_test, y_train, y_test):

    print(' running scipy linear model')

    regr = linear_model.LinearRegression()


    # Create linear regression object

    # Train the model using the training sets
    regr.fit(x_train, y_train)

    # Make predictions using the testing set
    y_pred = regr.predict(x_test)

    # The coefficients
    print('Coefficients: \n', regr.coef_)
    # The mean squared error
    print("Mean squared error: %.2f" % mean_squared_error(y_test, y_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y_test, y_pred))

    # Plot outputs
    plt.scatter(x_test, y_test,  color='black')
    plt.plot(x_test, y_pred, color='blue', linewidth=3)

    plt.xticks(())
    plt.yticks(())

    plt.savefig('resultlinear_reg.png')




def main(): 

    x_train, x_test, y_train, y_test = get_data()

    #print(x_train, ' x_train')
    #print(y_train, ' y_train')
    #print(x_test, ' x_test')

    scipy_linear_mod(x_train, x_test, y_train, y_test)


    x_train, x_test, y_train, y_test = generate_syntheticdata()

    #print(x_train, ' x_train')
    #print(y_train, ' y_train')
    #print(x_test, ' x_test')
    #print(y_test, ' x_test')


    numpy_linear_mod(x_train, x_test, y_train, y_test)







if __name__ == '__main__':
    main()