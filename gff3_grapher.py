# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 11:34:05 2022

@author: jarod
"""


#Enformer Results Grapher


#Regions of interest (50kb around):

#chr17 4997480	4997720
#chr7	100693549	100694000
#chr1	150487200	150487440
#chr17	80254669	80255033
#chr12	9764723	9765172

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
    
    
def contr_score_per_transcript():
    
    fig, ax = plt.subplots(figsize=(50,10))

    plt.bar([i*2 for i in range(153)], contr['score'])

    plt.xticks([i*2 for i in range(153)], contr['Transcript ID'])

    plt.yticks(size=20)

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    
    #plt.bar(contr['Name'],contr['score'], width=0.08)
    plt.title('Contribution Scores for Areas 50kb around the Interval',size=30)
    plt.xlabel('Transcript ID',size=25)
    plt.ylabel('Contribution Score', size=25)
    plt.show()
    
def highlighted_transcript(full_contr, pos_report, high_start, high_end):
    
    index = np.arange(896)
    
    #mask1 = high_end >= index >= high_start
    #mask2 = high_end <= index <= high_start
    
    mask1 = index >= high_start
    mask3 = high_end >= index
    mask2 = index <= high_start
    mask4 = high_end <= index
    
    fig, ax = plt.subplots(figsize=(50,10))

    plt.bar(index[mask1], full_contr['Score'][mask1], color = 'red')
    plt.bar(index[mask3], full_contr['Score'][mask3], color = 'red')
    plt.bar(index[mask2], full_contr['Score'][mask2], color = 'blue')
    plt.bar(index[mask4], full_contr['Score'][mask4], color = 'blue')

    plt.yticks(size=30)
    plt.xticks(size=30)

    #for tick in ax.get_xticklabels():
     #   tick.set_rotation(90)

    plt.title('Pooled Contribution Scores for 114688 bp Window',size=40)
    plt.xlabel('Position',size=40)
    plt.ylabel('Contribution Score', size=40)
    plt.show()
    

loc = pd.read_csv('pos_report.csv')

full_contr_pass = pd.read_csv('contr_score_152.csv', dtype = str, header=None)

contr = pd.DataFrame(index = range(len(full_contr_pass)), columns=['Score'], dtype = np.float64)

for x in range(len(full_contr_pass)):
    
    contr['Score'][x] = float(full_contr_pass[0][x])
    
    
#transcript_id = pd.DataFrame(index=range(len(contr)), columns=['Transcript ID'])

# for x in range(len(contr)):
    
#     temp_line = str(contr['Gff3 Line'][x])
    
#     counter = 0
    
#     start = 0
    
#     end = 0
    
#     for y in range(len(temp_line)):
        
#         pos = str(temp_line[y])
        
#         if pos ==';':
            
#             counter = counter + 1
            
#             if counter == 3 and [pos == '=']:
                
#                 start = y + 15
                
#                 end = y + 32
                
                
#     transcript_id['Transcript ID'][x] = str(temp_line[start:end])



#contr['Transcript ID'] = transcript_id


start = ((80255033-80254669)/2) + 80254669 - (114688/2)

for x in [152]:
    
    start = ((80255033-80254669)/2) + 80254669 - (114688/2)
    
    begin = (loc['begin'][x]-start)/128
    
    end = (loc['end'][x]-start)/128
    
    highlighted_transcript(contr, loc, begin, end)
    
    
        
