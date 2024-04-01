# IF STATEMENTS
is_hot = False
is_cold = False

if is_hot:
    print("It's a hot day")
    print("Drink water")
# press shift and tab to exit this block (fuori indentazione)
elif is_cold:
    print("It's a cold day")
    print("Wear warm clothes")
else:
    print("It's a lovely day")


temperature = input("How is the temperature today?\n")

if temperature.lower() == "hot":
    print("Hot?? Drink water!")
elif temperature.lower() == "cold":
    print("Cold?? Wear warm clothes!")
else:
    print("Allora suca")