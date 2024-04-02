import os
os.system("clear")

numbers = [5, 2, 4, 1]

numbers.append(20) # adds the specified item at the end of the list

numbers.insert(0, 3) # inserts 3 in position 0

numbers.remove(4) # removes the specified item

# numbers.clear()

numbers.pop() 

print(numbers.index(5)) # check the existence of an item, giving as output its position in the list, otherwise in generates an error
print(7 in numbers) # same, but gives only a boolean value

numbers.insert(3, 5)

print(numbers.count(5))

print(numbers.sort()) # none
numbers.sort() # re-orders items in ascending order
print(numbers)

numbers.reverse() # re-orders items in descending order
print(numbers)

numbers2 = numbers.copy()
numbers.append(10)

print(numbers2)
print(numbers)