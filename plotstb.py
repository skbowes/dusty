import matplotlib.pyplot as plt
import numpy as np
nmodels = 5
inputfile = 'sphere6'
# Read the data from the file
data = np.array([])
# print(data) 
with open(inputfile+'.stb', 'r') as file:
    model = -1
    columns = {}
    for line in file:
        if line.strip().startswith('#') and line.strip().endswith('SPECTRUM'):
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


column_headers = ['lambda', 'fTot', 'xAtt', 'xDs', 'xDe', 'fInp', 'tauT', 'albedo']
column_descriptors = ['$\\lamda$', '$\\lambda*F_{\\lambda}/int(F_{\\lambda}d\lambda)$', 'fractional contribution of the attenuated input radiation to fTot', \
                      ' fractional contribution of the scattered radiation to fTot', 'fractional contribution of the dust emission to fTot', \
                         'spectral shape of the input (unattenuated) radiation', 'overall optical depth at wavelength $\\lambda$', \
                             'albedo at wavelength lambda' ]
# Plot the data
# plt.figure(figsize=(10, 6))
# # for i in range(nmodels):
# #     plt.plot(list(zip(*columns[i]))[0], list(zip(*columns[i]))[4], label=f'Model {i+1}')
# plt.plot(lambda_col1, fTot_col1, marker='o', markersize = 3, color='b', label='Model 1')
# plt.plot(lambda_col2, fTot_col2, marker='o',  markersize = 3,color='g', label='Model 2')

### Spectrum file



# #plotting the dust temperature profile:
# plt.figure(figsize=(10, 6))  
# for i in range(nmodels):
#     plt.plot(lambda_cols[i], data[i][5], label=f'$\\tau$={tau0_values[i]}') # marker='o', markersize=3,
# plt.xscale('log')
# plt.xlabel('$\\lambda$')
# plt.legend()
# plt.ylabel(column_headers[5])
# plt.title(column_descriptors[5])
# plt.grid(True)
# plt.show()
# print(lambda_cols)
freq_cols = []
for lambda_col in lambda_cols:
    freq_col = [3e8 / lam for lam in lambda_col]
    freq_cols.append(freq_col)

    # Create a new object to store the result of the division
    divided_data = []

    # Divide all elements in data[:][1] by the corresponding element in lambda_cols
    for i in range(nmodels):
        divided_model_data = [d / lam for d, lam in zip(data[i][1], lambda_cols[i])]
        divided_data.append(divided_model_data)


print(data[:][1])


plt.figure(figsize=(10, 6))  
for i in range(nmodels):
    #data[i][1]
    plt.plot(freq_cols[i], divided_data[i], label=f'$\\tau$={tau0_values[i]}') # marker='o', markersize=3,
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\\nu$')
plt.legend()
plt.ylabel(column_headers[1])
plt.title('*F_{\\nu}/int(F_{\\lambda}d\lambda)$')
plt.grid(True)
plt.show()


# plt.figure(figsize=(10, 6))
# for j in range(1, len(data[0])):  
#     for i in range(nmodels):
#         plt.plot(lambda_cols[i], data[i][j], label=f'$\\tau$={tau0_values[i]}') # marker='o', markersize=3,
#     plt.xscale('log')
#     plt.xlabel('$\\lambda$')
#     plt.legend()
#     plt.ylabel(column_headers[j])
#     plt.title(column_descriptors[j])
#     plt.grid(True)
#     plt.show()



