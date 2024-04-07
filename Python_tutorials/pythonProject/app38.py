import os
os.system('clear')

def greet_user():
    first_name = input("Insert your first name: ")
    last_name = input("Insert your last name: ")
    print(f"Hi {first_name} {last_name}!")
    print("Welcome aboard")


greet_user() # positional argument: the order matters! (not in this case because of the input build)


def insult(first_name, last_name):
    print(f"Oh, now that I remember: {first_name} {last_name}, you're a bitch!")


insult(last_name="Smith", first_name="John") # key ward arguments: doesn't matter the order!
# can use key ward argument after positional argument but not the opposite

# calc_cost(total = 50, shipping = 5, discount = 0.1) # example of improving readibility