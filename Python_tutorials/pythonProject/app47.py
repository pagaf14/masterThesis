import os
os.system("clear")

from utils import find_max

list = [2, 46, 32, 58, 1, 0]
print(find_max(list))

# in python there is built-in max in matlab

print(max(list))

# but if

max = find_max(list)

print(max(list)) # ==> TypeError: 'int' object is not callable