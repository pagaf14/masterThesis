
birth_year = input("What year were you born?\n")

current_year = input("What is the current year?\n")

age = int(current_year) - int(birth_year)
print("Your age is:\n" + str(age))

# float() to convert a string into a float number (decimal)
# bool() to convert a string into a boolean value (true or false)

print("The classe of the variable age is:\n")
print(type(age))
print("The classe of the variable birth_year is:\n")
print(type(birth_year))

boh = 2.9
print("The classe of the variable boh is:\n")
print(type(boh))

acqua = False
print("The classe of the variable acqua is:\n")
print(type(acqua))