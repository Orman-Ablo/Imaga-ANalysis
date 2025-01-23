import pandas as pd #pandas is a python library that is commonly used to store and index data
import matplotlib.pyplot as plt # the most frequently used plotting library, we are only using the pyplot section.
import numpy as np # numpy is an almost ubiquitous maths-based library, useful for certain formulas (not needed, but out of habit I added it)
#open all files
#if we collect more data, we only need to open it and add it to the list and 
one_b = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results 1b.csv")
one_c = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results 1c (removed weird slice 2).csv")
one_d = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results 1d (shaky).csv")
one_e = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results Cells 1e.csv")
one_f = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results Cells 1f.csv")
two_f = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results 2f.csv")
#initializing some lists to be placed in the final dataframe
Control = [two_f,one_f]
start_cell_count = 0
#turns controls into one sample
one_f.set_index('Label',inplace=True)
two_f.set_index('Label',inplace=True)
wpb_1 = one_f.filter(like='VWF:Cell',axis=0)
wpb_2 = two_f.filter(like='VWF:Cell',axis=0)
cells = one_f.filter(like='-', axis = 0)
cells_2f = two_f.filter(like='-', axis = 0)

total_cells = cells.shape[0]

wpb_2.reset_index(inplace=True)
cells_2f.reset_index(inplace=True)
one_f.reset_index('Label',inplace=True)




wpb_2.loc[:,'Number'] = wpb_2['Label'].str.extract(r'(\d+)$').astype(int) + total_cells

# Update df2's Labels with the new numbers
wpb_2.loc[:,'Label'] = wpb_2['Label'].str.replace(r'(\d+)$', '', regex=True) + wpb_2['Number'].astype(str)

Control = pd.concat([one_f,cells_2f,wpb_2], ignore_index=True)
print(Control)
