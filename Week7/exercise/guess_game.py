#!/usr/bin/env python

'''
A Python program  

'''
import numpy as np 

import matplotlib.pyplot as plt


import random

pi = 3.14 # global variable 

def plot_examples():
  x = [1, 2, 3, 4]
  #numpy_x = 
  plt.plot(x)
  plt.ylabel('some numbers')
  #plt.show()
  plt.savefig('x.png')

  plt.clf() # make sure you have this in all plots

def area_circle(radius):

  area = pi * radius * radius

  return area 

def numpy_lists():

 x = np.zeros(15)
 y = np.ones(15)
 z = x + y

 k = y * 10

 '''print(x)
 print(y)
 print (z)
 print(k)'''

 magic = np.random.rand(3,4)
  
 magic_three = np.random.rand(3,4,2)
 # homework, write a function for summing 3D magic
 magic_one = np.random.rand(10)

 #print(magic)
 #print(magic_one)
 print(magic_three)
 print(magic.shape)

 magic_sum = sum_numpy(magic)
 print( magic_sum, ' is magic sum')

 

def sum_3D(a):
  print('homework - sum elements of 3d numpy array')

  # use nested for loops

  #https://www.ict.social/python/basics/multidimensional-lists-in-python









def sum_numpy(a):

  sum = 0
  
  for x in range(a.shape[0]):
    for y in range(a.shape[1]):
        print(x, y, ' *' , a[x][y] )
        sum = sum + a[x][y]
  return sum 

def numpy_design(filename):

  #a = np.zeros((10,10))

  b = np.random.rand(10,10)
 
  for x in range(1, b.shape[0]):
    for y in range(1, b.shape[1], 2):
        b[x][y] = y*x 

  double_mat =   b * 2

 #https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html
  np.savetxt(filename, double_mat, delimiter = ' ',  fmt='%1.2f' )
  np.savetxt('mymagic.csv', double_mat, delimiter = ',',  fmt='%1.2f' )
  np.savetxt('mymagic_fmt.csv', double_mat, delimiter = ',',  fmt='%1.5f' )

def advanced_numpy():

  #https://docs.scipy.org/doc/numpy/reference/generated/numpy.vstack.html

  a = np.array([1, 2, 3])
  b = np.array([2, 3, 4])

  vertical_stack = np.vstack((a,b))
  print(vertical_stack, '  vertical_stack')

  horizontal_stack = np.hstack((a,b))
  print(horizontal_stack, '  horizontal_stack')

  #print(double_mat)

def load_file(filename):

 #https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html
  data_load = np.loadtxt(filename)

  print(data_load, ' loaded data')

  plt.imshow(data_load, cmap='hot', interpolation='nearest')
  plt.savefig('design_magic.png')
  plt.clf() # make sure you have this in all plots

  #https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html

  data_tras = np.transpose(data_load)


  plt.imshow(data_tras, cmap='hot', interpolation='nearest')
  plt.savefig('design_magic_traspose.png')
  plt.clf() # make sure you have this in all plots
 
  
        
 
    



def sum_lists(a):

  #a = [[1, 12, 3], [51, 6, 17], [7, 18, 9]]
  sum = 0
  
  for x in range(len(a)):
    for y in range(len(a[x])):
        print(x, y, ' *' , a[x][y] )
        sum = sum + a[x][y]
  return sum 
    



def speed_function_one(speed):

  if speed < 80:
    print " speed is ok"
  else:
    print " you. have to pay a fine"


def speed_function_two(speed):

   
        
    if speed < 80:
      print " speed is ok"

    elif (speed >= 80) and (speed < 100):
      base_fine = 200
      speed_diff = speed - 80
      extra_fine = speed_diff * 10
      total_fine = extra_fine + base_fine 
      print " fine of  ", total_fine
    elif (speed >= 100) and (speed < 140):
      print " fine of 600 AUD"
    else:
      print " time to go to prison. "

   









def guess_function(n):

    guess = int(raw_input("Enter an integer from 1 to 99: "))
    while n != "guess":
    
        if guess < n:
            print "guess is low"
            guess = int(raw_input("Enter an integer from 1 to 99: "))
        elif guess > n:
            print "guess is high"
            guess = int(raw_input("Enter an integer from 1 to 99: "))
        else:
            print "you guessed it!"
            break 


def main():
    #print("Guess game.")


    #n = random.randint(1, 99)

    #guess_function(n)

    #area = area_circle(5)

    #print(area)

    #speed = int(input("Enter speed "))
    

    #speed_function_two(speed)

 
    #x = [[1,2],[3,5]]

    #sum = sum_lists(x)

    numpy_lists()

    print(sum, ' is sum')

    filename = 'mymagic.txt' # filename to  output data and also read data later



    numpy_design(filename)

    plot_examples()

    load_file(filename)

    advanced_numpy()





if __name__ == '__main__':
    main()