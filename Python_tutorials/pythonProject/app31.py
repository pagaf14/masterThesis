# dictionaries

import os
os.system("clear")

customer = {
    "name": "Matteo Paganelli",
    "age": 24, # each key should be unique in a dictionary
    "is_verified": True # the value could be a string, a value, a boolean, a list, ... anything
}

print(customer["name"])
print(customer.get("email")) # if there is no email it returns None instead of raising an error
print(customer.get("birthdate", "Jan 1 1980")) # can give default value, if not defined
customer["name"] = "Franco"
print(customer["name"])