import os

os.system('clear')
pwd = os.getcwd()
print("La directory di lavoro corrente è:", pwd)

print("Ciao mondo, sono Bot! @.@")

messaggio_simpatico = "Oggi non mi sento Atzeco"
print(messaggio_simpatico)

cibo_prefe = "la fessa"
print("Il mio cibo preferito è " + cibo_prefe)

print("\nEhi tu chi cazzo sei?")

# La funzione input prende un'informazione dalla console e la dà ad una variabile
nome = input()

if  nome == "Rocco Siffredi":
    print("\nComplimenti per la ciolla!  Ti ammiro molto.")
elif nome == "Valentina Nappi":
    print("\nMazza che topa!")
else:
    print("\nChe cazzo vuoi " + nome + ", tornatene a casa")

lunghezza_nome = len(nome)
lunghezza_nome_stringa = str(lunghezza_nome)
print("Il tuo nome ha " + lunghezza_nome_stringa + " lettere")

# le variabili possono essere di tipo str se contengono una stringa 
# o di tipo int se contengono un numero intero

#print("\nIn che anno sei nato?")
#anno_di_nascita = input()
#Equivalente:
anno_di_nascita = input("\nIn che anno sei nato? ")
anno_corrente = input("\nIn che anno ci troviamo? ")
eta = int(anno_corrente) - int(anno_di_nascita)
print(nome + "! Hai esattamente " + str(eta) + " anni ...")

if(eta < 18):
    print("...beato te!!")
    anni_da_compiere = eta + 1 #Occhio all'indentazione!
    print("\nL'anno prossimo avrai " + str(anni_da_compiere) + " anni")
else:
    print("...sei proprio un vecchietto!")
    eta_Silente = 150
    print("\nTra " + str(eta_Silente - eta) + " anni, sarai vecchio come Albus Silente")

#fino a che l'istruzione tra parentesi sarà vera, le istruzioni indentate saranno eseguite,
# altrimenti il ciclo finisce. Se mettiamo True, l'istruzione saarà sempre vera
risposta = ""
while(risposta != "chi è?"): #!= significa diversa
    print("\ntoc toc...")
    risposta = input()
print("STO CAZZO")

animali = ["gatto", "cane", "topo", "ragno"]
print("\nQual è il tuo animale preferito?")
for animale in animali:
    print(animale)
print("...indicalo con un numero da 0 a 3") #gli elementi in python sono numerati da 0!
numero_animale = int(input())
if(numero_animale < 0 or numero_animale > 3):
    print("Non hai capito proprio una sega")
else:
    animale_scelto = animali[numero_animale]
    print(animale_scelto+ "...ottima scelta!")