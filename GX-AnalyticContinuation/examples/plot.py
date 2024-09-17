"""  
Copyright (C) 2020-2024 GreenX library
This file is distributed under the terms of the APACHE2 License.

This script is intendet to plot the output of pade_example.f90

usage: python plot.py <data_file> <output_image_path>

"""
import numpy as np 
import matplotlib.pyplot as plt 
from sys import argv

#------------------------------------------------------------------------------ 
# parse arguments
#------------------------------------------------------------------------------ 

data_file = argv[1]
img_path = argv[2]

#------------------------------------------------------------------------------ 
# load data
#------------------------------------------------------------------------------ 

data = np.loadtxt(data_file, comments="#")

x_re = data[:, 0]
x_im = data[:, 1]
y_re = data[:, 2]
y_im = data[:, 3]
y_ref_re = data[:, 4]
y_ref_im = data[:, 5]

#------------------------------------------------------------------------------ 
# plot
#------------------------------------------------------------------------------ 

plt.figure(figsize=(7, 7))

# real part
ax1 = plt.subplot(211)
ax1.axhline(0, color='k', linewidth=0.7)
ax1.axvline(0, color='k', linewidth=0.7)
ax1.plot(x_re, y_ref_re, '-', linewidth=3.5, color="gray", label="Re(y_ref)")
ax1.plot(x_re, y_re, '--', color="tab:orange", label="Re(y_interpolated)")
ax1.legend(fontsize=16, loc=4)
ax1.set_xlabel("Re(x)", fontsize=16)
ax1.set_ylabel("Re(y)", fontsize=16)


# imaginary part
ax2 = plt.subplot(212)
ax2.axhline(0, color='k', linewidth=0.7)
ax2.axvline(0, color='k', linewidth=0.7)
ax2.plot(x_re, y_ref_im, '-', linewidth=3.5, color="gray", label="Im(y_ref)")
ax2.plot(x_re, y_im, '--', color="tab:orange", label="Im(y_interpolated)")
ax2.legend(fontsize=16, loc=1)
ax2.set_xlabel("Re(x)", fontsize=16)
ax2.set_ylabel("Im(y)", fontsize=16)

plt.tight_layout()
plt.savefig(img_path, dpi=200)