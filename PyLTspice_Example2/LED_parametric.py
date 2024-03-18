# LED_parametric.py
# example of using the PyLTSpice library to loop through a spice model and plot the output
# Luminus Devices November 2022
# URL: https://luminusdevices.zendesk.com/hc/en-us/articles/10818558523661-Data-Analysis-Using-Python-to-run-LTspice-as-a-remote-process
# 
# Process:
# 1. Create a schematic template in LTspice that includes at least one paramater. Do not assign any default values in the .asc file.
# 2. Run the model parameters in a loop to create <projectname>_n.raw output files
# 3. extract desired output data from the raw files as lists
# 4. Plot the data  

# declarations
import os
import glob
import numpy as np
from matplotlib import pyplot as plt

# PyLTSpice can be found at https://github.com/nunobrum/PyLTSpice
from PyLTSpice.LTSpiceBatch import SimCommander
from PyLTSpice.LTSpice_RawRead import LTSpiceRawRead


# get this script's absolute path
meAbsPath = os.path.dirname(os.path.realpath(__file__))

# Initialize model variables defined in the .asc file
Rsh_list=[1.5e3, 1e4, 1e5, 1e6, 1.73937e8] # we will loop through this list of shunt resistances - the other variables are constants in this example but could also be lists if we modify the code to do multiple loops
IS=1e-19
RS=0.1
N=2

# Create the SimCommander object "LTC" for the .asc file
LTC = SimCommander(meAbsPath + "\\Voltage_Source_LED_sweep_parametric.asc")
print('LTC created')

# loop through the list of Rsh values and run the model for each one
# this will create a "raw" file for each run with _1, _2, _3, ..., _n appended to the base name of the model
# it also creates netlist files for each iteration with the same numbering scheme
for value in Rsh_list:
    LTC.set_parameters(Rsh_value = value) # modify the value of a named parameter in the asc file - this does not work if there are default statements in the asc file - you need to comment them out or you get multiple decimal points in the net file
    LTC.set_parameters(IS_value = IS)
    LTC.set_parameters(RS_value = RS)
    LTC.set_parameters(N_value = N)
    LTC.run()
    LTC.wait_completion() 
    LTC.reset_netlist() # skipping this reset step produces errors 

# Get the file names of the raw files that were created for this model run 
raw_file_list= glob.glob(meAbsPath+'\\Voltage_Source_LED_sweep_parametric_*.raw')

# read each raw file and extract the information needed to plot the curves
voltage_list=[]
LED_current_list=[]
shunt_current_list=[]
total_current_list=[]

for filename in raw_file_list: 
    LTR = LTSpiceRawRead(filename)
    
    # these are the three traces we want to plot - they are data structures at this point (not directly plottable)
    # you can programatically get these names, I know what they are named so I hard coded them
    Vin = LTR.get_trace('V(n001)')
    IR1 = LTR.get_trace("I(R1)") 
    ID1 = LTR.get_trace("I(D1)")
        
    voltage = Vin.get_wave(0) # <--------------- This is how to get one trace from the LTR object. If you are not running steps, 0 is the only trace.
    voltage.tolist() # <------------------------ Converting this to a list makes the plot routine happy 
    voltage_list.append(voltage) # <------------ This creates a list of lists
    
    LED_current=ID1.get_wave(0)
    LED_current.tolist()
    LED_current_list.append(LED_current)
    
    shunt_current=IR1.get_wave(0)
    shunt_current.tolist()
    shunt_current_list.append(shunt_current)
    
    total_current_list.append(shunt_current+LED_current)


# We need to delete some files to prevent latency. Deleting the numbered raw files should be enough for this example
# - comment this out if you want to look at the raw files or the timestamps
# - This could be made more bulletproof, it is mainly here to to deal with the case where the number of Rsh values is decreased between runs
# - This loop runs before the plot command because plot is the one that likes to crash when there are list length / data type mismatches. 
#   This may happen if the code is modified for a different case study
for filename in raw_file_list:
    print('deleting ' + filename)
    os.remove(filename)
    
# Plotting each trace in a loop
x=voltage # This is actually the last calculated voltage - since we are not changing the sweep in this analysis, this is OK
y=total_current_list # This is a list of lists
for i in range(len(y)):
    plt.plot(x,y[i],label = 'run {}: Rsh = {:,.2E} Ohms'.format(i+1,Rsh_list[i]))
plt.plot([0,2],[1e-5,1e-5], '--k', label='Reference line (10 uA)') # this is a reference line
plt.ylabel('Currrent (A)')
plt.xlabel('Voltage (V)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()



