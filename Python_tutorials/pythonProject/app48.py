import os
os.system("clear")

from ecommerce.shipping import calc_shipping as cs

# otherwise import the whole module

from ecommerce import shipping



shipping.calc_shipping()
cs()