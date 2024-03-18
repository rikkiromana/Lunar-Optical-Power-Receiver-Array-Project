# declarations
import os
import ltspice
import numpy as np
from matplotlib import pyplot as plt

# PyLTSpice can be found at https://github.com/nunobrum/PyLTSpice
from PyLTSpice import SimCommander

# get this script's absolute path
meAbsPath = os.path.dirname(os.path.realpath(__file__))

# Initialize model variables defined in the .asc file
illu = 700 # All cells are receiving 700 W/m^2 illumination as initial test

# Create the SimCommander object "LTC" for the .asc file
LTC = SimCommander(meAbsPath + "\\NineCell.asc")
print('LTC created')

# Setting values in the .asc file
LTC.set_parameters(illu1 = illu)
LTC.set_parameters(illu2 = illu)
LTC.set_parameters(illu3 = illu)
LTC.set_parameters(illu4 = illu)
LTC.set_parameters(illu5 = illu)
LTC.set_parameters(illu6 = illu)
LTC.set_parameters(illu7 = illu)
LTC.set_parameters(illu8 = illu)
LTC.set_parameters(illu9 = illu)
LTC.run()
LTC.wait_completion() 
LTC.reset_netlist() # skipping this reset step produces errors 

LTR = ltspice.Ltspice(meAbsPath + '\\NineCell_1.raw')
LTR.parse()

voltage = LTR.get_data('v(vout)')
print(voltage)
current = LTR.get_data('I(I1)')
print(current)

# Power Calculation
power = voltage * current

# Set Up Plot
plt.subplot(2,1,1)
plt.title('I-V Curve')
plt.ylabel('Current (A)')
plt.plot(voltage, current)
plt.subplot(2,1,2)
plt.title('P-V Curve')
plt.xlabel('Voltage (V)')
plt.ylabel('Power (W)')
plt.plot(voltage, power)
plt.show()