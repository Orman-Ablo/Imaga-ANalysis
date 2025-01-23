import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#open all files
Results = pd.read_csv("C:/Users/Naam/Documents/CSBB Gene Pool Party/Results images/Results Area Background Filtering.csv")
Area = Results[['Label','Area','%Area']]
Area.iloc[1],Area.iloc[2] = Area.iloc[2],Area.iloc[1]

ax = plt.gca()
ax.set_ylim([0, 10])
fig = plt.gcf()
fig.set_size_inches(5,6)
# for row, column in Area_percent.items():
plt.bar(Area['Label'],Area['%Area'], width= 0.8)
plt.xlabel('Conditions')
plt.ylabel('Percentage of Area VWF')
plt.title('Percentage of image area occupied by VWF')
for i in range(len(Area['%Area'])):
    plt.text(Area['Label'][i], Area['%Area'][i],'%'+str(Area['%Area'][i]))
plt.show()
# for labels in Area['Label']:
