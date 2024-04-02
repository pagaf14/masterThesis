import os
os.system("clear")

Weight = input("Weight: ")
Unit = input("(L)bs or (K)g: ")

if Unit.upper() == "L":
    print("You are " + str(Weight) + " pounds")
elif Unit.upper() == "K":
    Weight = int(Weight) * (1/0.45)
    #print("You are " + str(Weight) + " pounds")
    print(f"You are {Weight} pounds")
else:
    print("Invalid input for measurement unit")