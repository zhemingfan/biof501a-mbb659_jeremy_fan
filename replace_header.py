#!/usr/bin/env python

fasta= open('combined.fasta')
newfasta= open('new_combined.fasta', 'w')

i = 1
for line in fasta:
    if line.startswith('>'):
        newfasta.write('>'+ str(i) + '\n')
        i+=1
    else:
        newfasta.write(line)

fasta.close()
newfasta.close()