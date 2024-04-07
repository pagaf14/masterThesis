import os
os.system("clear")

import random

for i in range(3):
   print(random.random())
   
for i in range(3):
    print(random.randint(10, 20))

members = ["Matteo", "Marco", "Giovanni"]
leader = random.choice(members)
print(leader)