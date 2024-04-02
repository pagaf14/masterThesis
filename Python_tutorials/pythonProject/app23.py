import os
os.system("clear")

numbers = [5, 2, 5, 2, 2]

#for item in numbers:
#    print("x"*item)

for x_count in numbers:
    output = ""
    for count in range(x_count):
        output += "x"
    print(output)
