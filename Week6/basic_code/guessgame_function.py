#!/usr/bin/env python

'''
A Python program  

'''



import random

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
    print("Guess game.")


    n = random.randint(1, 99)

    guess_function(n)




if __name__ == '__main__':
    main()
