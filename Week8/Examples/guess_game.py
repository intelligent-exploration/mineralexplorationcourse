'''
A Python program  

'''



import random

pi = 3.14 # global variable 

def area_circle(radius):

  area = pi * radius * radius

  return area 

def speed_function(speed):

    guess = int(raw_input("Enter an integer from 1 to 99: "))
    while n != "guess":
    
        if speed < n:
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


    n = random.randint(1, 99)

    #guess_function(n)

    area = area_circle(5)

    print(area)




if __name__ == '__main__':
    main()