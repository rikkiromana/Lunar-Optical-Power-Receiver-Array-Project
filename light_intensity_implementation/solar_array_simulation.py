import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import os

# Defining cell and beam properties
cell_width_cm = 1.35
cell_height_cm = 1.6
beam_center_x_cm = 0
beam_center_y_cm = 0
beam_waist_cm = 0.7

# Defining geometry of the array parameters
gap_size_cm = 0.1  # Assumption
num_cells_x = 3
num_cells_y = 3

# Material and optical properties
glass_refractive_index = 1.51446

# Total laser power in Watts (W)
P_total = 500

def gaussian_intensity(x, y, P_total, beam_center_x_cm=beam_center_x_cm, beam_center_y_cm=beam_center_y_cm, beam_waist_cm=beam_waist_cm):
    sigma = beam_waist_cm / np.sqrt(2 * np.log(2))
    normalized_intensity = np.exp(-((x - beam_center_x_cm) ** 2 + (y - beam_center_y_cm) ** 2) / (2 * sigma ** 2))
    peak_intensity_W_per_cm2 = P_total / (np.pi * sigma**2)
    return normalized_intensity * peak_intensity_W_per_cm2

def calculate_fresnel_losses(n1, n2, theta_i=0):
    cos_theta_i = np.cos(np.deg2rad(theta_i))
    sin_theta_t = n1 / n2 * np.sin(np.deg2rad(theta_i))
    cos_theta_t = np.sqrt(1 - sin_theta_t ** 2)
    rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
    rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
    fresnel_coefficient = (rs + rp) / 2
    total_loss = fresnel_coefficient ** 2
    return 1 - total_loss

def integrate_light_distribution_with_glass(x_center, y_center):
    x_min = x_center - cell_width_cm / 2
    x_max = x_center + cell_width_cm / 2
    y_min = y_center - cell_height_cm / 2
    y_max = y_center + cell_height_cm / 2
    
    def integrand(x, y):
        light_intensity = gaussian_intensity(x, y, P_total)
        fresnel_adjustment = calculate_fresnel_losses(1, glass_refractive_index)
        return light_intensity * fresnel_adjustment
    
    result, _ = dblquad(integrand, x_min, x_max, lambda x: y_min, lambda x: y_max)
    return result

def iv_curve_silicon_realistic(voltage, light_intensity):
    Is = 1e-12
    n = 1.5
    Vt = 0.026
    Is_modified = Is * light_intensity
    return Is_modified * (np.exp(voltage / (n * Vt)) - 1)

def calculate_photocurrent(light_intensity):
    max_voltage = 0.6
    photocurrent = iv_curve_silicon_realistic(max_voltage, light_intensity)
    return photocurrent

def simulate_array_and_export():
    light_intensities = np.zeros((num_cells_x, num_cells_y))
    total_photocurrent = 0
    
    print("Light Intensities for Each Cell (W/cm^2):")
    for i in range(num_cells_x):
        for j in range(num_cells_y):
            x_center = (cell_width_cm + gap_size_cm) * i + cell_width_cm / 2 - ((num_cells_x * cell_width_cm + (num_cells_x - 1) * gap_size_cm) / 2)
            y_center = (cell_height_cm + gap_size_cm) * j + cell_height_cm / 2 - ((num_cells_y * cell_height_cm + (num_cells_y - 1) * gap_size_cm) / 2)
            light_intensity = integrate_light_distribution_with_glass(x_center, y_center)
            light_intensities[i, j] = light_intensity
            
            print(f"Cell ({i+1}, {j+1}): {light_intensity:.5f} W/cm^2")
            
            photocurrent = calculate_photocurrent(light_intensity)
            total_photocurrent += photocurrent

    # Exporting the light intensities to a file
    specific_path = r"C:\Users\rikki\Lunar-Optical-Power-Receiver-Array-Project\light_intensity_implementation"
    os.makedirs(specific_path, exist_ok=True)  # Ensure the directory exists
    filename = os.path.join(specific_path, "light_intensities.txt")
    np.savetxt(filename, light_intensities.flatten(), fmt="%0.5f")

    print(f"Total Photocurrent for the Array: {total_photocurrent} A")

def visualize_gaussian_distribution():
    # Recalculate array dimensions
    array_width_cm = num_cells_x * cell_width_cm + (num_cells_x - 1) * gap_size_cm
    array_height_cm = num_cells_y * cell_height_cm + (num_cells_y - 1) * gap_size_cm

    x = np.linspace(-array_width_cm / 2, array_width_cm / 2, 300)
    y = np.linspace(-array_height_cm / 2, array_height_cm / 2, 300)
    X, Y = np.meshgrid(x, y)
    # Updated to pass P_total to gaussian_intensity directly within X, Y calculations
    Z = gaussian_intensity(X, Y, P_total)
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, Z, levels=50, cmap='viridis')
    plt.colorbar(label='Intensity (W/cm^2)')
    plt.title('Gaussian Light Distribution Over the Solar Cell Array')
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    simulate_array_and_export()
    visualize_gaussian_distribution()
