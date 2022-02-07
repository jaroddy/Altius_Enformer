# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 15:16:21 2022

@author: jarod
"""


#GFF3 Parser

#Regions of interest (50kb around):

#chr17 4997480	4997720
#chr7	100693549	100694000
#chr1	150487200	150487440
#chr17	80254669	80255033
#chr12	9764723	9765172


import pandas as pd

def parser(gff3_str):
    
    tab_count = 0
    
    begin = int()
    
    end = int()
    
    chrom = str()
    
    transcript = str()
    
    for pos in range(len(gff3_str)):
        
        if gff3_str[pos] == '\t':
        
            tab_count = tab_count + 1
            
            if tab_count == 1:
                
                chrom = gff3_str[0:pos]
                
            elif tab_count == 2:
                
                trans_beg = pos
            
            elif tab_count == 3:
            
                start_begin = pos
                
                transcript = gff3_str[trans_beg+1:pos]
            
            elif tab_count == 4:
            
                end_begin = pos
            
                begin = gff3_str[start_begin+1:end_begin]
            
            elif tab_count == 5:
                
                end_end = pos
                
                end = gff3_str[end_begin+1:end_end]
    
    return begin, end, chrom, transcript
            

gff3 = pd.read_csv("v39.annotation.gff3",delimiter='/')

pos_report = pd.DataFrame(index=[0], columns=['begin','end', 'Chromosome', 'Type','Gff3 Line'])

chrom_list = ['chr17','chr7','chr1','chr17','chr12']

bounds_lower = [int(4997480), int(100693549), int(150487200), int(80254669), int(9764723)]

bounds_higher = [int(4997720), int(100694000), int(150487440), int(80255033), int(9765172)]

lower_bound = int()

upper_bound = int()


for x in range(len(gff3)):
    
    pos_report_temp = pd.DataFrame(index=[0], columns=['begin','end', 'Chromosome', 'Type','Gff3 Line'])
    
    temp_line = gff3.iloc[x]
    
    begin, end, chrom, transcript = parser(temp_line[0])
    
    begin = int(begin)
    
    end = int(end)
    
    for y in range(len(chrom_list)):
        
        lower_bound = bounds_lower[y] - 50000
        
        upper_bound = bounds_higher[y] + 50000
    
        if (chrom_list[y] == chrom) and (begin >= lower_bound) and (end <= upper_bound):
            
            print(x)
    
            pos_report_temp['begin'][0] = int(begin)
    
            pos_report_temp['end'][0] = int(end)
    
            pos_report_temp['Gff3 Line'][0] = temp_line[0]
        
            pos_report_temp['Chromosome'][0] = chrom
    
            pos_report_temp['Type'][0] = transcript
    
            pos_report = pos_report.append(pos_report_temp)
            
            
pos_report = pos_report[pos_report['Type'] == 'transcript']
    
pos_report.to_csv('pos_report.csv', index=False)
        
    
        
    
    
    
    
    
    
    
    



    
    
