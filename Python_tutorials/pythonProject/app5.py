weight_kg = input("Weight (kg):\n")
# weight_lbs = weight_kg * (1/0.45) # TypeError: can't multiply sequence by non-int of type 'float'
weight_lbs = int(weight_kg) * (1/0.45)
print(weight_lbs)

