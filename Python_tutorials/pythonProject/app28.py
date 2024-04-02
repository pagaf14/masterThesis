import os
os.system("clear")

numbers = [3, 2, 4, 4, 6, 5, 1]
uniques = []

for number in numbers:
    if number not in uniques:
        uniques.append(number)
print(uniques)