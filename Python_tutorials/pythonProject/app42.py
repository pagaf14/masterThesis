import os
os.system("clear")

# we use classes to define new types
# basic types: numbers, strings, booleans; lists, dictionaries

class Point:  # not use underscore but capitalize every word (ex: EmailClient) (naming convention)
    def move(self):
        print("move")
    
    def draw(self):
        print("draw")


point1 = Point() # this creates object
point1.draw()
point1.move()
point1.x = 10 # "x" = attrtibute of the object point1
point1.y = 20
print(point1.x)

point2 = Point()
print(point2.x) # error: attribute x not defined on object point2