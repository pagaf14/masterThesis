# how to handle errors

import os
os.system("clear")

try:
    age = int(input("Age: "))
    income = 20000
    risk = income / age # if age = 0 ==> ZeroDivisionError: division by zero
    print(age)
except ZeroDivisionError:
    print("\n\nYou must enter a non-zero value for your age.")


except ValueError:
    print("Invalid value")


# if input not a numer ==> ValueError: invalid literal for int() with base 10
# exit code 1 = crash, exit code 0 = success

