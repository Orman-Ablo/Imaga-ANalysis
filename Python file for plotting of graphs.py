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
list = [one_f,one_e,one_d]
norms = []
stdev = []
Feret = []
Feret_stdev = []
Circle = []
Circle_stdev = []
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
    #I didn't put the results of the cells into their own results. This boolean checks the WPB counts only
    #Bool_wb = results['Label'].str.startswith("VWF:Cell")
    #this one only considers the cells themselves to be True
    #cells = results['Label'].str.startswith("VWF:0")
    #count the number of cells
    cell_amount = filt_results_Cell.shape[0]
    #the list for theratio of area to WPB counts
    norm = []
    #checks the WPB counts per cell
    print(type(results1))
    filt_results_WB.reset_index(inplace=True)
    filt_results_Cell.reset_index(inplace=True)
 
    label_counts = filt_results_WB['Label'].value_counts()
    #for each cell, it normalizes the data
    for i in range(cell_amount):
        # print(label_counts['VWF:Cell '+str(i+1)])
        #thankfully, as the cells are analyzed first, they always sit in the beginning of the index
        # i takes the row index, 'Area' only checks the value for the 'Area column
        area = filt_results_Cell.loc[i,'Area' ]
        #Temporary solution: since cells have different sizes, this provides a simple ratio for cell area to WPB count
        normalize = label_counts['VWF:Cell '+str(i+1)]/area
        #add each ratio to the normalization list
        norm.append(normalize)
    #calculate the standard deviation and average    
    std_norm = np.std(norm)
    avg_norm= np.mean(norm)
    #add to the list
    norms.append(avg_norm)
    stdev.append(std_norm)
#Place all the lists into a single dataframe for easier/ more efficeint storage for later use
final = pd.DataFrame({'image':names,'normal':norms,'Standard Deviation': stdev, 'Feret':Feret, 'Feret stdev':Feret_stdev,'Circularity':Circle, 'Circularity stdev':Circle_stdev})
#make the index the images taken
final.set_index('image')
print(final)
#plotting our results: scatter plot, with the standard deviation as errorbars
figure, axis = plt.subplots(2, 2)
axis[0,0].scatter(final['image'],final['normal'])
axis[0,0].errorbar(final['image'],final['normal'],yerr=final['Standard Deviation'],fmt='bo',ecolor='r',capsize=6)
axis[0,0].set_xlabel('images(a)')
axis[0,0].set_ylabel('Normalized WPB Counts')
axis[0,0].set_title('Comparison of WPB Counts \n Across Set Conditions')
axis[1,0].remove()
axis[0,1].remove()
axis[1,1].remove()
# axis[0,1].scatter(final['image'],final['Feret'])
# axis[0,1].errorbar(final['image'],final['Feret'],yerr=final['Feret stdev'],fmt='bo',ecolor='r',capsize=6)
# axis[0,1].set_xlabel('images(b)')
# axis[0,1].set_ylabel(u'Average Feret Lengths (\u03bcm)')
# axis[0,1].set_title('Comparison of Averaged Feret \n Lengths Across Set Conditions')

# axis[1,0].scatter(final['image'],final['Circularity'])
# axis[1,0].errorbar(final['image'],final['Circularity'],yerr=final['Circularity stdev'],fmt='bo',ecolor='r',capsize=6)
# axis[1,0].set_xlabel('images(c)')
# axis[1,0].set_ylabel('Average Circularity')
# axis[1,0].set_title('Comparison of Averaged \n Circularities Across Set Conditions')
# figure.tight_layout(pad=2.0)

plt.show()

   




