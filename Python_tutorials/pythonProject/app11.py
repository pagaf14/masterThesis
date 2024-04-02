house_price = 1e6
print("The house price is " + str(house_price) + "€ \n")

buyer_credit = input("How is the buyer income?\n")

if buyer_credit.lower() == "good":
    print("Allright! The buyer has to put down 10%, so " + str(house_price * 0.1) + "€")
else:
    print("Ouch! The buyer has to put down 20%, so " + str(house_price*0.2) + "€")

buyer_income = input("Now, has the buyer good credit?\n")

if buyer_income.lower() == "yes" and buyer_credit.lower() == "good":
    print("Daje Roma daje!")
else:
    print("Va bene lo stessooooo")