import os
os.system("clear")

class Point:

    # constructor = function that gets called at the time of creating an object
    # init method to initialize the object of a class

    def __init__(self, x, y): # "self" = reference to the current object
        self.x = x  # like point.x = 10
        self.y = y

    def move(self):
        print("move")

    def draw(self):
        print("draw")


point = Point(10, 20) # x and y coordinates
point.x = 11
print(point.x)

