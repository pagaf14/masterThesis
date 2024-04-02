# 2 dimensional lists

import os
os.system("clear")

# [
#     1 2 3
#     4 5 6
#     7 8 9
# ]

matrix = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
]

# matrix[0][1] = 20
# print(matrix[0][1]) # access second item of first list of items

for row in matrix:
    for item in row:
        print(item)