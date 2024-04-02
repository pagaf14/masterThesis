import os
os.system("clear")

coordinates = (1, 2, 3)
# coordinates[0] * coordinates[1] * coordinates[2] ... too long
# x = coordinates[0]
# y =  coordinates[1]
# z = coordinates[2]
x, y, z, = coordinates # short hand to achieve the same result, we're unpacking the tuple in 3 items
print(y)