# Purpose: This code simulates a Gaussian laser beam's light distribution over a 3x3 solar cell array and calculates the total photocurrent generated.
# It helps in aligning solar panels to maximize energy capture from the beam.

import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

# defining cell and beam
cell_width_cm = 1.35
cell_height_cm = 1.6
beam_center_x_cm = 0
beam_center_y_cm = 0
beam_waist_cm = 0.7

# defining geometry of the array parameters
gap_size_cm = 0.1 #assumption
num_cells_x = 3
num_cells_y = 3
array_width_cm = num_cells_x * cell_width_cm + (num_cells_x - 1) * gap_size_cm
array_height_cm = num_cells_y * cell_height_cm + (num_cells_y - 1) * gap_size_cm


# material and optical properties
glass_refractive_index = 1.51446

# gaussian light distribution function
# this function models the spread of light intensity from the center of the beam, getting dimmer as it moves away from the center.
def gaussian(x, y, beam_center_x_cm=beam_center_x_cm, beam_center_y_cm=beam_center_y_cm, beam_waist_cm=beam_waist_cm):
    sigma = beam_waist_cm / np.sqrt(2 * np.log(2))
    return (1 / (2 * np.pi * sigma ** 2)) * np.exp(-((x - beam_center_x_cm) ** 2 + (y - beam_center_y_cm) ** 2) / (2 * sigma ** 2))

# fresnel losses calculation
# it calculates how much light is lost when passing from air into the glass of the solar cell, depending on the materials' properties.
def calculate_fresnel_losses(n1, n2, theta_i=0):
    cos_theta_i = np.cos(np.deg2rad(theta_i))
    sin_theta_t = n1 / n2 * np.sin(np.deg2rad(theta_i))
    cos_theta_t = np.sqrt(1 - sin_theta_t ** 2)
    rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
    rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
    fresnel_coefficient = (rs + rp) / 2
    total_loss = fresnel_coefficient ** 2
    return 1 - total_loss

# integration of light distribution over a cell, adjusted for Fresnel losses
# the function calculates how much light from the beam actually hits a cell, taking into account the Gaussian distribution and the Fresnel losses.
def integrate_light_distribution_with_glass(x_center, y_center):
    x_min = x_center - cell_width_cm / 2
    x_max = x_center + cell_width_cm / 2
    y_min = y_center - cell_height_cm / 2
    y_max = y_center + cell_height_cm / 2
    
    def integrand(x, y):
        light_intensity = gaussian(x, y)
        fresnel_adjustment = calculate_fresnel_losses(1, glass_refractive_index)
        return light_intensity * fresnel_adjustment
    
    result, _ = dblquad(integrand, x_min, x_max, lambda x: y_min, lambda x: y_max)
    return result

# placeholder I-V curve for a silicon photodiode
# this is a model that shows how much electrical current is generated for a given voltage in a silicon photodiode, which is a component that converts light into electrical current.
def iv_curve_silicon(voltage):
    Is = 1e-12  # Saturation current
    Vt = 0.026  # Thermal voltage at room temperature
    return Is * (np.exp(voltage / Vt) - 1)

# photocurrent calculation from light intensity, adjusted for silicon photodiodes
# this function uses the I-V curve to determine the electrical current generated by the amount of light that hits each cell.
def calculate_photocurrent(light_intensity):
    # Maximum voltage for a typical silicon photodiode
    max_voltage = 0.6
    photocurrent = iv_curve_silicon(max_voltage) * light_intensity
    return photocurrent

# simulating the entire array
# this loops over each cell in the array, calculates the light intensity for each, and then sums up the total current that would be generated by the entire array.
def simulate_array():
    total_photocurrent = 0
    for i in range(num_cells_x):
        for j in range(num_cells_y):
            x_center = (cell_width_cm + gap_size_cm) * i + cell_width_cm / 2 - ((num_cells_x * cell_width_cm + (num_cells_x - 1) * gap_size_cm) / 2)
            y_center = (cell_height_cm + gap_size_cm) * j + cell_height_cm / 2 - ((num_cells_y * cell_height_cm + (num_cells_y - 1) * gap_size_cm) / 2)
            light_intensity = integrate_light_distribution_with_glass(x_center, y_center)
            photocurrent = calculate_photocurrent(light_intensity)
            total_photocurrent += photocurrent
    print(f"Total Photocurrent for the Array: {total_photocurrent} A")

# visualization of Gaussian distribution across the array
# this generates a visual map of how light is spread across the array, showing the Gaussian distribution from the beam.
def visualize_gaussian_distribution():
    x = np.linspace(-array_width_cm / 2, array_width_cm / 2, 300)
    y = np.linspace(-array_height_cm / 2, array_height_cm / 2, 300)
    X, Y = np.meshgrid(x, y)
    Z = gaussian(X, Y)
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, Z, levels=50, cmap='viridis')
    plt.colorbar(label='Intensity')
    plt.title('Gaussian Light Distribution Over the 3x3 Solar Cell Array')
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    simulate_array()
    visualize_gaussian_distribution()

# the output image will be a color-coded visual created by the code that shows how the intensity of the light varies across the array
# the brightest area is at the center, where the beam is strongest, and it fades outwards
# this helps understand where the light is most concentrated and where it is not, which is useful for designing the array layout/aligning the light source properly