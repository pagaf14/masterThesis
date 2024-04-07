import os
os.system("clear")

# def lbs_to_kg(weight):
#     return weight * 0.45

# def kg_to_lbs(weight):
#     return weight / 0.45

# weight_kg =  float(input("\nEnter your weight in kilograms: "))
# print(f"\nYour weight in pounds is approximately: {kg_to_lbs(weight_kg)}")

import converters

print(converters.kg_to_lbs(78))


# or

from converters import kg_to_lbs # in this way don't need to call the function with prefix of the module name

print(kg_to_lbs(78))


