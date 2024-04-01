import math

print(10 + 3)
print(10 / 3)
print(10 // 3)
print(10 % 3)
print(10 ** 3)

x = 10
x = x + 3
print(x)

y = 10
y += 3
print(y)

z = 10
z -= 3
print(z)

x = 10 + 3 * 2 ** 2
# operations precedence (parenthesis > exponential > multiplication/division > summation/substraction
print(x)

x = 2.9
print(round(x))
print(abs(-x))

print(math.ceil(x)) # approssimazione per eccesso (ceiling)
print(math.floor(x)) # approssimazione per difetto (floor)
print(math.fabs(-x)) # absolute value