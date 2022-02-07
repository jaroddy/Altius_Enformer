# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 11:34:05 2022

@author: jarod
"""


#Enformer Results Grapher

import pandas as pd
import matplotlib.pyplot as plt



contr = pd.read_csv('contributions_report_final.csv')

transcript_id = pd.DataFrame(index=range(len(contr)), columns=['Transcript ID'])

for x in range(len(contr)):
    
    temp_line = str(contr['Gff3 Line'][x])
    
    counter = 0
    
    start = 0
    
    end = 0
    
    for y in range(len(temp_line)):
        
        pos = str(temp_line[y])
        
        if pos ==';':
            
            counter = counter + 1
            
            if counter == 3 and [pos == '=']:
                
                start = y + 15
                
                end = y + 32
                
                
    transcript_id['Transcript ID'][x] = str(temp_line[start:end])

contr['Transcript ID'] = transcript_id
        
fig, ax = plt.subplots(figsize=(50,10))

plt.bar([i*2 for i in range(153)], contr['score'])

plt.xticks([i*2 for i in range(153)], contr['Transcript ID'])

plt.yticks(size=20)

for tick in ax.get_xticklabels():
    tick.set_rotation(90)
    
#plt.bar(contr['Name'],contr['score'], width=0.08)
plt.title('Contribution Scores for Areas 50kb around the Interval',size=30)
plt.xlabel('Interval',size=2)
plt.ylabel('Contribution Score', size=25)
plt.show()