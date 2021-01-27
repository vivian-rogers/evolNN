import math
import random


dx = 0.3


outfile = "outfile.csv"
print('Printing atoms to ' + outfile + "\n")
f = open(outfile,"w") 
f.write("\n")
for i in range(0,280):
   #f.write(str(i) + ',')
   data = []
   start = random.random()
   x = start*2*math.pi
   w = 1.2
   #w = random.random()/3 + 1
   for j in range(0,9):
      data.append(math.floor(8* (math.sin(w*x)+1)/2 ))
      x += dx
   print("w = " + str(w) + " " + str(data))
   line = str(data[0])
   for j in range(1,len(data)):
       line = line + ',' + str(data[j])
   line += "\n"
   f.write(line)
f.close() 

#   for j in range(1,3)
#      data.append(sin(x

