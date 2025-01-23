import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

results = pd.read_csv('C:/Users/Naam/Documents/CSBB Gene Pool Party/ELISA Values/ELISA Data.csv',sep=';')
print(results)


plt.bar(results['Sample'],results['nM'])

plt.xticks(rotation=30)
plt.xlabel('Conditions')
plt.ylabel('nM')
plt.title('ELISA Results VWF Concentration')
plt.tight_layout()
plt.show()