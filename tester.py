# https://www.w3schools.com/python/python_regex.asp
import re

s = 'PREDICTED: Danio rerio roundabout, axon guidance receptor, homolog 2 (Drosophila) (robo2), transcript variant X11, mRNA'

## this searches for PREDICTED.
x = re.search('PREDICTED',s)

if x:
    print('YES!')
else:
    print('NO')

match_str = x.group(0) #https://stackoverflow.com/questions/18493677/how-do-i-return-a-string-from-a-regex-match-in-python
#print(match_str)


## this searches for PREDICTED.
x = re.search('\([0-9A-Za-z]*\),',s)

if x:
    print('YES!')
else:
    print('NO')

match_str = x.group(0)
print(match_str)

name = match_str[1:len(match_str)-2]
print(name)