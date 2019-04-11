#!/usr/bin/env python

'''
A Python program  

'''



import random

pi = 3.14 # global variable 

def area_circle(radius):

  area = pi * radius * radius

  return area 

def sum_lists():

  a = [[1, 12, 3], [51, 6, 17], [7, 18, 9]]
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

 

    sum = sum_lists()
    print(sum, ' is sum')





if __name__ == '__main__':
    main()