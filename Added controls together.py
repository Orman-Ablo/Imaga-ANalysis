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

#Turns the Controls into one control set
#This code predominantly served t oreplace the values of the WPB labels
#so the total amount of cells across both images is accounded for 
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



#Update number
wpb_2.loc[:,'Number'] = wpb_2['Label'].str.extract(r'(\d+)$').astype(int) + total_cells
#Add number
wpb_2.loc[:,'Label'] = wpb_2['Label'].str.replace(r'(\d+)$', '', regex=True) + wpb_2['Number'].astype(str)

Control = pd.concat([one_f,cells_2f,wpb_2], ignore_index=True)


list = [Control,one_e,one_d]
ratios = []
stdev_ratios = []
Feret = []
Feret_stdev = []
Circle = []
Circle_stdev = []
Normalize = []
stdev_norm=[]
names = ['Control','10\u03bcM:\nday 1','10\u03bcM 4:\nday 4']

# it will perform the same set of commands per image, to extract the same information
for condition in list:
    results = condition[['Label','Area','Feret','Circ.']]
    results1 = results
    results.set_index('Label',inplace=True)
    filt_results_WB = results.filter(like='VWF:Cell',axis=0)
    filt_results_Cell = results.filter(like='-', axis = 0)
    temp = float(np.mean(filt_results_WB['Feret']))
    Feret.append(temp)
    temp = filt_results_WB['Feret'].std()
    Feret_stdev.append(temp)
    temp = float(np.mean(filt_results_WB['Circ.']))
    Circle.append(temp)
    temp = float(np.std(filt_results_WB['Circ.']))
    Circle_stdev.append(temp)
    #count the number of cells
    cell_amount = filt_results_Cell.shape[0]
    #the list for theratio of area to WPB counts
    norm = []
    rat = []
    #checks the WPB counts per cell
    filt_results_Cell.reset_index(inplace=True)
    filt_results_WB.reset_index(inplace=True)
    label_counts = filt_results_WB['Label'].value_counts()
    
    #for each cell, it normalizes the data
    for i in range(cell_amount):
        filt_results_WB.set_index('Label',inplace=True)
        #thankfully, as the cells are analyzed first, they always sit in the beginning of the index
        # i takes the row index, 'Area' only checks the value for the 'Area column
        area = filt_results_Cell.loc[i,'Area' ]
        #Temporary solution: since cells have different sizes, this provides a simple ratio for cell area to WPB count
        ratio = label_counts['VWF:Cell '+str(i+1)]/area
        #add each ratio to the normalization list
        rat.append(ratio)
        relative =filt_results_WB.filter(like='VWF:Cell '+ str(i),axis=0)
        filt_results_WB.reset_index(inplace=True)
        Areas = relative['Area']
        norm.append((sum(Areas)/area)*100)

    #calculate the standard deviation and average   
    
    std_ratios = np.std(rat)
    
    avg_ratios= np.mean(rat)
    

    std_norm = np.std(norm)
    
    avg_norm= np.mean(norm)
    print(avg_norm)


    
    #add to the list
    ratios.append(avg_ratios)
    stdev_ratios.append(std_ratios)
    Normalize.append(avg_norm)
    stdev_norm.append(std_norm)
    
#Place all the lists into a single dataframe for easier/ more efficeint storage for later use
final = pd.DataFrame({'image':names,'ratio':ratios,'Standard Deviation': stdev_ratios, 'Feret':Feret, 'Feret stdev':Feret_stdev,'Circularity':Circle, 'Circularity stdev':Circle_stdev,'Normal':Normalize,'Normal deviation':stdev_norm})
#make the index the images taken
final.set_index('image')
#print(final)
#plotting our results: scatter plot, with the standard deviation as errorbars
extract = final[['ratio','Standard Deviation','Normal','Normal deviation']]
print(extract)
figure, axis = plt.subplots(2, 2)
axis[1,0].remove()
axis[1,1].remove()
axis[0,0].scatter(final['image'],final['ratio'])
axis[0,0].errorbar(final['image'],final['ratio'],yerr=final['Standard Deviation'],fmt='bo',ecolor='r',capsize=6)
axis[0,0].set_xlabel('images(a)')
axis[0,0].set_ylabel('Ratio WPB Counts per \u03bcm^2')
axis[0,0].set_title('Comparison of WPB Counts \n Across Set Conditions')

axis[0,1].scatter(final['image'],final['Normal'])
axis[0,1].errorbar(final['image'],final['Normal'],yerr=final['Normal deviation'],fmt='bo',ecolor='r',capsize=6)
axis[0,1].set_xlabel('images(b)')
axis[0,1].set_ylabel(u'Normalized areas (%)')
axis[0,1].set_title('Comparison of Normalized \n Area Coverage Across\n Set Conditions')

x = [1,1.5,2]
figure.suptitle('An analyisis on WPB count and Area Coverage')


plt.show()

   




