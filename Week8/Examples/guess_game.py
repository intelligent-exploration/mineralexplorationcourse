#!/usr/bin/env python

'''
A Python program  

'''



import random

pi = 3.14 # global variable 

def area_circle(radius):

  area = pi * radius * radius

  return area 


def speed_function_one(speed):

  if speed < 80:
    print " speed is ok"
  else:
    print " you. have to pay a fine"



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

    speed = int(input("Enter speed "))
    

    speed_function_one(speed)




if __name__ == '__main__':
    main()