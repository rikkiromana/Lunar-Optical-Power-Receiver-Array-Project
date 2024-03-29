import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import os

cell_width_m = 0.0135
cell_height_m = 0.016
beam_center_x_m = 0.0153
beam_center_y_m = 0.0188
beam_waist_m =  5e-3 # Initial beam waist (m)
gap_size_m = 0.001  # Assumption
distance_to_cells_m = 100
lambda_laser_m = 1520e-9 # Laser wavelength in meters, assuming a 1520nm laser

num_cells_x = 3
num_cells_y = 3

# material and optical properties
glass_refractive_index = 1.51446

P_total = 500

# gaussian intensity function using global variables
def gaussian_intensity(x, y):
    global P_total, distance_to_cells_m, lambda_laser_m, beam_center_x_m, beam_center_y_m, beam_waist_m
    w_z = beam_waist_m * np.sqrt(1 + (lambda_laser_m * distance_to_cells_m / (np.pi * beam_waist_m**2))**2)
    normalized_intensity = np.exp(-2 * ((x - beam_center_x_m) ** 2 + (y - beam_center_y_m) ** 2) / (w_z ** 2))
    peak_intensity_W_per_m2 = (2 * P_total) / (np.pi * w_z**2)
    return normalized_intensity * peak_intensity_W_per_m2

def calculate_fresnel_losses(n1, n2, theta_i=0):
    cos_theta_i = np.cos(np.deg2rad(theta_i))
    sin_theta_t = n1 / n2 * np.sin(np.deg2rad(theta_i))
    cos_theta_t = np.sqrt(1 - sin_theta_t ** 2)
    rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
    rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
    fresnel_coefficient = (rs + rp) / 2
    total_loss = fresnel_coefficient ** 2
    return 1 - total_loss

# function to integrate light distribution over a cell using only x_center and y_center as arguments
def integrate_light_distribution_with_glass(x_center, y_center):
    global cell_width_m, cell_height_m, P_total, lambda_laser_m, glass_refractive_index, beam_waist_m
    x_min = x_center - cell_width_m / 2
    x_max = x_center + cell_width_m / 2
    y_min = y_center - cell_height_m / 2
    y_max = y_center + cell_height_m / 2

    def integrand(y, x):
        light_intensity = gaussian_intensity(x, y)
        fresnel_adjustment = calculate_fresnel_losses(1, glass_refractive_index)
        return light_intensity * fresnel_adjustment

    result, error = dblquad(integrand, x_min, x_max, lambda x: y_min, lambda x: y_max)
    return result

def gaussian_intensity_with_fresnel(x, y):
    intensity = gaussian_intensity(x, y)
    fresnel_losses = calculate_fresnel_losses(1, glass_refractive_index, theta_i=0)
    intensity_with_fresnel = intensity * fresnel_losses
    return intensity_with_fresnel

def simulate_array_and_export():
    light_intensities = np.zeros((num_cells_y, num_cells_x))
    
    print("Light Intensities for Each Cell (W/m^2):")
    for j in range(num_cells_y - 1, -1, -1):  # Start from the top
        for i in range(num_cells_x):
            # calculating cell center, adjusting the formula if necessary
            x_center = (cell_width_m + gap_size_m) * i + cell_width_m / 2 - ((num_cells_x * cell_width_m + (num_cells_x - 1) * gap_size_m) / 2)
            y_center = (cell_height_m + gap_size_m) * j + cell_height_m / 2 - ((num_cells_y * cell_height_m + (num_cells_y - 1) * gap_size_m) / 2)
            
            light_intensity = integrate_light_distribution_with_glass(x_center, y_center)
            light_intensities[j, i] = light_intensity  # Adjust index for correct assignment
            
            # adjusting print statement to reflect correct cell numbering from top to bottom, left to right
            print(f"Cell ({num_cells_y - j},{i + 1}): {light_intensity:.12f} W/m^2")

    # exporting the light intensities to a file
    '''specific_path = r"C:\Users\rikki\Lunar-Optical-Power-Receiver-Array-Project\light_intensity_implementation"
    os.makedirs(specific_path, exist_ok=True)  # Ensure the directory exists
    filename = os.path.join(specific_path, "light_intensities.txt")
    np.savetxt(filename, light_intensities.flatten(), fmt="%0.12f")'''

    return light_intensities


# visualization function including cell numbering
def visualize_array_and_beam():
    global num_cells_x, num_cells_y, cell_width_m, cell_height_m, gap_size_m, beam_center_x_m, beam_center_y_m, P_total, lambda_laser_m, beam_waist_m, distance_to_cells_m, glass_refractive_index

    array_width_m = num_cells_x * cell_width_m + (num_cells_x - 1) * gap_size_m
    array_height_m = num_cells_y * cell_height_m + (num_cells_y - 1) * gap_size_m

    # creating a grid for the array
    x = np.linspace(-array_width_m / 2, array_width_m / 2, 500)
    y = np.linspace(-array_height_m / 2, array_height_m / 2, 500)
    X, Y = np.meshgrid(x, y)

    # calculating the intensity at each point on the grid
    Z = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Z[i, j] = gaussian_intensity_with_fresnel(X[i, j], Y[i, j])

    # setting up the figure and plot
    plt.figure(figsize=(10, 8))
    # for aspect
    plt.gca().set_aspect('equal', adjustable='box')

    min_intensity = Z.min()
    max_intensity = Z.max()
    levels = np.linspace(min_intensity, max_intensity, 1000)
    contour = plt.contourf(X, Y, Z, levels=levels, cmap='viridis', vmin=min_intensity, vmax=max_intensity)
    plt.colorbar(contour, label='Irradiance (W/m²)')
    plt.title('Gaussian Beam Distribution and Array Geometry')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')

    # visualizing the geometry of the array and number each cell
    for i in range(num_cells_x):
        for j in range(num_cells_y - 1, -1, -1):  # start from the top
            cell_x = (cell_width_m + gap_size_m) * i - array_width_m / 2
            cell_y = (cell_height_m + gap_size_m) * j - array_height_m / 2
            plt.gca().add_patch(plt.Rectangle((cell_x, cell_y), cell_width_m, cell_height_m, linewidth=1, edgecolor='r', facecolor='none'))
            cell_number = f"{num_cells_y - j},{i + 1}"
            plt.text(cell_x + cell_width_m / 2, cell_y + cell_height_m / 2, cell_number, ha='center', va='center', color='white')

    plt.grid(True)
    plt.xlim([-array_width_m / 2, array_width_m / 2])
    plt.ylim([-array_height_m / 2, array_height_m / 2])
    plt.show()


if __name__ == "__main__":
    simulate_array_and_export()
    visualize_array_and_beam()
