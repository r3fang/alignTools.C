import sys
MISMATCH = -1
MATCH = 1
GAP = -3
EXTENSION = -1

def encode(input_string):
    count = 1
    prev = ''
    lst = []
    for character in input_string:
        if character != prev:
            if prev:
                entry = (prev,count)
                lst.append(entry)
                #print lst
            count = 1
            prev = character
        else:
            count += 1
    else:
        entry = (character,count)
        lst.append(entry)
    return lst

lines = sys.stdin.readlines()
seq1 = lines[0].strip()
seq2 = lines[1].strip()
score = 0

for(x,y) in zip(seq1, seq2):
    if(x==y): 
        score += MATCH
    if(x!=y and x!= '-' and y!= '-'): 
        score+=MISMATCH
for(ch, num) in encode(seq1):
    if(ch=='-'):
        score = score + GAP + EXTENSION*num

for(ch, num) in encode(seq2):
    if(ch=='-'):
        score = score + GAP + EXTENSION*num

print score