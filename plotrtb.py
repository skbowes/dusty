import matplotlib.pyplot as plt
import numpy as np
nmodels = 5
inputfile = 'sphere4'
# Read the data from the file
data = np.array([])
# print(data) 
with open(inputfile+'.rtb', 'r') as file:
    model = -1
    columns = {}
    for line in file:
        if line.strip().startswith('#') and line.strip().endswith('PROFILES'):
            model += 1
            columns[model] = []
        # Split the line into columns and convert to float
        if line.strip().startswith('#'):
            continue
        columns[model].append(list(map(float, line.split())))


tau0_values = []

with open(inputfile+'.spp', 'r') as file:
    for line in file:
        if line.strip().startswith('#'):
            continue
        taucol = line.split()
        tau0_values.append(float(taucol[1]))

# Now you can call on tau0_values
print(tau0_values)

data = []
lambda_cols = []

for i in range(nmodels):
    model_data = list(zip(*columns[i]))
    data.append(model_data)
    lambda_cols.append(model_data[0])


column_headers = ['y', 'eta', 't', 'tauF', 'epsilon', 'Td', 'rg']
column_descriptors = ['dimensionless radius', 'dimensionless normalized radial density profile', \
                      'radial profile of optical depth variation', 'radial profile of the flux-averaged optical depth', \
                        'fraction of grain heating due to the contribution of the envelope to the radiation field', \
                            'radial profile of dust temperature', 'radial profile of the ratio of P_rad to F_g']
# Plot the data
# plt.figure(figsize=(10, 6))
# # for i in range(nmodels):
# #     plt.plot(list(zip(*columns[i]))[0], list(zip(*columns[i]))[4], label=f'Model {i+1}')
# plt.plot(lambda_col1, fTot_col1, marker='o', markersize = 3, color='b', label='Model 1')
# plt.plot(lambda_col2, fTot_col2, marker='o',  markersize = 3,color='g', label='Model 2')

### Spectrum file

#plotting the dust temperature profile:
plt.figure(figsize=(10, 6))
for i in range(nmodels):
    plt.plot(lambda_cols[i], data[i][5], label=f'$\\tau$={tau0_values[i]}') # marker='o', markersize=3,
plt.xscale('log')
plt.yscale('log')
plt.xlabel('dimensionless radius')
plt.legend()
plt.ylabel(column_headers[5])
plt.title(column_descriptors[5])
plt.grid(True)
plt.show()

# for i in range(nmodels):
#     plt.plot



for j in range(1, len(data[0])):  
    plt.figure(figsize=(10, 6))
    for i in range(nmodels): 
        plt.plot(lambda_cols[i], data[i][j], label=f'$\\tau$={tau0_values[i]}') # marker='o', markersize=3,
    plt.xscale('log')
    plt.xlabel('dimensionless radius')
    plt.legend()
    plt.ylabel(column_headers[j])
    plt.title(column_descriptors[j])
    plt.grid(True)
    plt.show()



