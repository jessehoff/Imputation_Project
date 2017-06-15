import sys

targ= sys.argv[1]
Chr = sys.argv[2]
outn = sys.argv[3]

file = open(targ,'r')

pile = file.read()

jin = pile.split(',')

out = open(outn,'w')
out.write(Chr.join(jin))
out.close()
       
