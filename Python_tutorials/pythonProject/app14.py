import os
os.system('clear')

temperature = 30

if temperature >= 30: # if you use only =, you are re-assigning the variable, not boolean
    print("It's a hot day")
else:
    print("It's not a hot day")