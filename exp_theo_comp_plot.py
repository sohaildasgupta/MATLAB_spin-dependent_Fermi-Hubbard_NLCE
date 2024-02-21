#%%
# Define the methods 
import numpy as np
import os
import matplotlib.pyplot as plt

def find_matching_filename(relative_folder_path, strings_to_match):
    """
    Find the file that has all the strings_to_match in the filename
    """
    folder_path = os.path.join(os.getcwd(), relative_folder_path)
    for filename in os.listdir(folder_path):
        if all(string in filename for string in strings_to_match):
            return filename
    return None  # Return None if no matching filename is found


def plot_experimental_data(U,N=3):
    """
    Plot experimental data from .txt files.
    """
    # Load data from text file
    current_directory = os.getcwd()
    #Check for the file with correct U (DOES NOT CARE FOR THE ORDER!! NEED TO UPDATE!)
    filename = find_matching_filename("../experimental_data/",U)
    relative_path = "../experimental_data/%s"%filename
    file_path = os.path.join(current_directory,relative_path)
    data = np.genfromtxt(file_path, skip_header=1)  # Skip the header line

# Extract columns
    x = data[:, 0]
    y1 = data[:, 1]
    y2 = data[:, 2]
    y3 = data[:, 3]
    error1 = data[:, 4]
    error2 = data[:, 5]
    error3 = data[:, 6]

# Plot data with error bars
    for i in range(N): 
        plt.errorbar(x, data[:,i+1], yerr=data[:,i+1+N] , fmt='o-', capsize=5 ,label=r'$n_%i$'%(i+1))

# Add labels and legend
    plt.xlabel(r'Distance')
    plt.ylabel(r'$\langle n_\sigma\rangle$')
    plt.legend()
    plt.grid(True)
    U = [x.replace("p",".") for x in U]
    plt.title(r"$U_{12}= %s, U_{13}= %s, U_{23} = %s$"%tuple(U))
# Show plot
    plt.show()

#%%
# plots 
u = [87, 32.6, 15.5]
u = [str(num).replace(".","p") for num in u]
plot_experimental_data(u,N=3)
