import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

# parameters-- change input depending which we will be using
cell_width_cm = 1.35  # width of the photovoltaic cell in cm
cell_height_cm = 1.6  # weight of the photovoltaic cell in cm
beam_center_x_cm = 0  # x-coordinate of the center of the Gaussian light beam in cm
beam_center_y_cm = 0  # y-coordinate of the center of the Gaussian light beam in cm
wavelength_nm = 1065  # wavelength of the laser in nm
average_power_w = 32  # average power of the laser in watts
optical_power_density_w_cm2 = 2  # optical power density in W/cm^2
output_nominal_voltage_v = 0.9  # output nominal voltage in volts
#output_nominal_current_a = 16  # output nominal current in amperes

# fefractive index of borosilicate glass (BK7?)
glass_refractive_index = 1.51446

# function to calculate Gaussian distribution
def gaussian(x, y):
    beam_waist_cm = cell_width_cm / 2  # assumption that beam waist is half of the cell width?
    sigma = beam_waist_cm / np.sqrt(2 * np.log(2))  # standard deviation for Gaussian beam
    return (1 / (2 * np.pi * sigma**2)) * np.exp(-((x - beam_center_x_cm)**2 + (y - beam_center_y_cm)**2) / (2 * sigma**2))

# function to calculate Fresnel coefficients for two reflections
def calculate_fresnel_coefficient(n1, n2, theta):
    cos_theta_i = np.cos(np.deg2rad(theta))
    cos_theta_t = np.sqrt(1 - ((n1 / n2) ** 2) * (1 - cos_theta_i ** 2))
    rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
    rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
    return (rs + rp) / 2

# function to integrate the light distribution over the cell area considering glass cover
def integrate_light_distribution_with_glass():
    result, _ = dblquad(lambda y, x: gaussian(x, y) * 
                        (1 - calculate_fresnel_coefficient(1, glass_refractive_index, 0)), 
                        -cell_width_cm/2, cell_width_cm/2,
                        lambda x: -cell_height_cm/2, lambda x: cell_height_cm/2)
    return result

# function to calculate the photocurrent generated in the cell
def calculate_photocurrent(light_intensity):
    return light_intensity * (wavelength_nm / 1000) * optical_power_density_w_cm2

def main():
    # calculating the distribution of light on the photovoltaic cell considering glass cover
    light_intensity = integrate_light_distribution_with_glass()
    print("Average Light Intensity on Cell with Glass Cover:", light_intensity)
    
    # calculating the photocurrent generated in the cell
    photocurrent = calculate_photocurrent(light_intensity)
    print("Photocurrent Generated in Cell:", photocurrent, "A")

    # plotting the Gaussian distribution of light
    x = np.linspace(-cell_width_cm/2, cell_width_cm/2, 100)
    y = np.linspace(-cell_height_cm/2, cell_height_cm/2, 100)
    X, Y = np.meshgrid(x, y)
    Z = gaussian(X, Y)
    plt.contourf(X, Y, Z, cmap='viridis')
    plt.colorbar(label='Intensity')
    plt.title('Gaussian Distribution of Light on Photovoltaic Cell')
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    plt.show()

if __name__ == "__main__":
    main()
