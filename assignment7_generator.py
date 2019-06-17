#! /usr/bin/env/python3

import random
import sys

introtext = '''
#! /bin/bash

# PART 1

# After mapping, we want to know how many reads come from which features in the genome. 
# We will use a tool called featureCounts from the subread package. 
# All of the following questions can be answered by going through the section about featureCounts
# in the manual, which you can find here: http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf 
# Please answer the questions by using `
# echo "ANSWER"

'''

PART2 = '''
# PART 2

# To speed up the rest of our course, please first install R,
# which you can find here: https://cran.r-project.org/

# Then download and install Rstudio from here:
# https://www.rstudio.com/products/rstudio/download/#download"

# If you want to be super prepared, please also take a look at SARtools,
# an R package, which we will use for differential gene expression analysis.
# Try to install it by following the instructions on their GitHub page:
# https://github.com/PF2-pasteur-fr/SARTools
# Their paper is also very informative:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157022

'''


class Question:

    def __init__(self, line):
        l = line.split(': ')
        self.q_num = l[0]
        self.q_text = l[1]
        
class Student:
    
    def __init__(self, matrikelnummer, Qs):
        self.MN = matrikelnummer.rstrip('\n')
        self.qset = random.sample(Qs, k=3)

questions = list()
students = list()

with open('../assignment7.txt', 'r') as qfile:
    for line in qfile:
        if line.startswith('Q'):
            questions.append(Question(line=line))

with open('../MN', 'r') as MN:
    for n in MN:
        students.append(Student(matrikelnummer=n, Qs=questions))

for s in students:
    with open(f'{s.MN}.txt', 'w') as outfile:
        outfile.write(introtext)
        
        for q in s.qset:
            outfile.write('# ' + q.q_text + '\n')
            
        outfile.write(PART2)
        
        q_ids = [q.q_num for q in s.qset]
        outfile.write('# ' + ' '.join(q_ids))



