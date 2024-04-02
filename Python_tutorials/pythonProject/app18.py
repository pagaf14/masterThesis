import os
os.system("clear")

secret_number = 9
guess_count = 0
guess_limit = 3
while guess_count < guess_limit:
    guess = int(input("Guess: "))
    guess_count += 1 # equivalent to guess_count = guess_count + 1
    if guess == secret_number:
        print("You won! LESGOO")
        break
else:
    print("You failed! HAHA")