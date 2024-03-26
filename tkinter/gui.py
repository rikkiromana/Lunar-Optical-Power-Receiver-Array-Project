import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import StringVar, Text
import os

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
    intensity_data = "Light Intensities for Each Cell (W/m^2):\n"  # Ensure this line is at the beginning of the function
    
    for j in range(num_cells_y - 1, -1, -1):  # Start from the top
        for i in range(num_cells_x):
            x_center = (cell_width_m + gap_size_m) * i + cell_width_m / 2 - ((num_cells_x * cell_width_m + (num_cells_x - 1) * gap_size_m) / 2)
            y_center = (cell_height_m + gap_size_m) * j + cell_height_m / 2 - ((num_cells_y * cell_height_m + (num_cells_y - 1) * gap_size_m) / 2)
            
            light_intensity = integrate_light_distribution_with_glass(x_center, y_center)
            light_intensities[j, i] = light_intensity
            
            intensity_data += f"Cell ({num_cells_y - j},{i + 1}): {light_intensity:.12f} W/m^2\n"

    return intensity_data

def visualize_array_and_beam():
    array_width_m = num_cells_x * cell_width_m + (num_cells_x - 1) * gap_size_m
    array_height_m = num_cells_y * cell_height_m + (num_cells_y - 1) * gap_size_m

    # creating a grid for array
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
    # setting the aspect of the plot to be equal
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
        for j in range(num_cells_y - 1, -1, -1):  # starting from the top
            cell_x = (cell_width_m + gap_size_m) * i - array_width_m / 2
            cell_y = (cell_height_m + gap_size_m) * j - array_height_m / 2
            plt.gca().add_patch(plt.Rectangle((cell_x, cell_y), cell_width_m, cell_height_m, linewidth=1, edgecolor='r', facecolor='none'))
            cell_number = f"{num_cells_y - j},{i + 1}"
            plt.text(cell_x + cell_width_m / 2, cell_y + cell_height_m / 2, cell_number, ha='center', va='center', color='white')

    plt.grid(True)
    plt.xlim([-array_width_m / 2, array_width_m / 2])
    plt.ylim([-array_height_m / 2, array_height_m / 2])
    plt.show()

def simulate_array_and_export():
    global num_cells_x, num_cells_y, cell_width_m, cell_height_m
    light_intensities = np.zeros((num_cells_y, num_cells_x))
    intensity_data = "Light Intensities for Each Cell (W/m^2):\n"  # make sure this line is at the beginning of the function
    
    for j in range(num_cells_y - 1, -1, -1):  # starting from the top
        for i in range(num_cells_x):
            x_center = (cell_width_m + gap_size_m) * i + cell_width_m / 2 - ((num_cells_x * cell_width_m + (num_cells_x - 1) * gap_size_m) / 2)
            y_center = (cell_height_m + gap_size_m) * j + cell_height_m / 2 - ((num_cells_y * cell_height_m + (num_cells_y - 1) * gap_size_m) / 2)
            
            light_intensity = integrate_light_distribution_with_glass(x_center, y_center)
            light_intensities[j, i] = light_intensity
            
            intensity_data += f"Cell ({num_cells_y - j},{i + 1}): {light_intensity:.12f} W/m^2\n"
    
    # exporting the light intensities to a file, will comment out for now but useful for automation later
    #specific_path = r"C:\Users\rikki\Lunar-Optical-Power-Receiver-Array-Project\tkinter"
    #os.makedirs(specific_path, exist_ok=True)  # Ensure the directory exists
    #filename = os.path.join(specific_path, "light_intensities.txt")
    #np.savetxt(filename, light_intensities.flatten(), fmt="%0.12f")

    return intensity_data


# GUI Functions
def update_parameters_and_run():
    # updating global parameters with values from the GUI
    global cell_width_m, cell_height_m, beam_center_x_m, beam_center_y_m
    global beam_waist_m, gap_size_m, distance_to_cells_m, lambda_laser_m
    global num_cells_x, num_cells_y, glass_refractive_index, P_total

    try:
        # getting values from the GUI and convert to appropriate types
        cell_width_m = float(cell_width_var.get())
        cell_height_m = float(cell_height_var.get())
        beam_center_x_m = float(beam_center_x_var.get())
        beam_center_y_m = float(beam_center_y_var.get())
        beam_waist_m = float(beam_waist_var.get())
        gap_size_m = float(gap_size_var.get())
        distance_to_cells_m = float(distance_to_cells_var.get())
        lambda_laser_m = float(lambda_laser_var.get())
        num_cells_x = int(num_cells_x_var.get())
        num_cells_y = int(num_cells_y_var.get())
        glass_refractive_index = float(glass_refractive_index_var.get())
        P_total = float(P_total_var.get())

        # run simulation and visualization
        simulate_array_and_export()
        visualize_array_and_beam()
    # run simulation and get intensity data
        intensity_data = simulate_array_and_export()
        visualize_array_and_beam()

        # updating the text widget with the intensity data
        output_text.delete("1.0", "end")  # clear the existing output
        output_text.insert("end", intensity_data)  # insert the new data

    except ValueError as e:
        output_text.delete("1.0", "end")  # clear the existing output
        output_text.insert("end", f"Error in input: {e}")

# GUI setup
root = tk.Tk()
root.title("Lunar Optical Power Receiver Array Simulation")

# container for all input fields
inputs_frame = tk.Frame(root)
inputs_frame.grid(row=0, column=0, padx=10, pady=10)

# dictionary to keep track of input variables
input_vars = {
    "Cell Width (m)": StringVar(value='0.0135'),
    "Cell Height (m)": StringVar(value='0.016'),
    "Beam Center X (m)": StringVar(value='0.0'),
    "Beam Center Y (m)": StringVar(value='0.0'),
    "Beam Waist (m)": StringVar(value='5e-3'),
    "Gap Size (m)": StringVar(value='0.001'),
    "Distance to Cells (m)": StringVar(value='100'),
    "Laser Wavelength (m)": StringVar(value='1520e-9'),
    "Number of Cells X": StringVar(value='3'),
    "Number of Cells Y": StringVar(value='3'),
    "Glass Refractive Index": StringVar(value='1.51446'),
    "Total Laser Power (W)": StringVar(value='500')
}

# dynamically create labels and entry fields based on the input_vars
for i, (label_text, var) in enumerate(input_vars.items()):
    tk.Label(inputs_frame, text=label_text).grid(row=i, column=0, sticky="w")
    tk.Entry(inputs_frame, textvariable=var).grid(row=i, column=1)

# text widget for displaying results
output_text = Text(root, height=15, width=50)
output_text.grid(row=1, column=0, padx=10, pady=10)

def update_parameters_and_run():
    # update global parameters from GUI inputs
    global cell_width_m, cell_height_m, beam_center_x_m, beam_center_y_m
    global beam_waist_m, gap_size_m, distance_to_cells_m, lambda_laser_m
    global num_cells_x, num_cells_y, glass_refractive_index, P_total
    
    cell_width_m = float(input_vars["Cell Width (m)"].get())
    cell_height_m = float(input_vars["Cell Height (m)"].get())
    beam_center_x_m = float(input_vars["Beam Center X (m)"].get())
    beam_center_y_m = float(input_vars["Beam Center Y (m)"].get())
    beam_waist_m = float(input_vars["Beam Waist (m)"].get())
    gap_size_m = float(input_vars["Gap Size (m)"].get())
    distance_to_cells_m = float(input_vars["Distance to Cells (m)"].get())
    lambda_laser_m = float(input_vars["Laser Wavelength (m)"].get())
    num_cells_x = int(input_vars["Number of Cells X"].get())
    num_cells_y = int(input_vars["Number of Cells Y"].get())
    glass_refractive_index = float(input_vars["Glass Refractive Index"].get())
    P_total = float(input_vars["Total Laser Power (W)"].get())

    # run the simulation and display results
    intensity_data = simulate_array_and_export()
    visualize_array_and_beam()
    output_text.delete("1.0", "end")
    output_text.insert("end", intensity_data)

# button to run the simulation
run_button = tk.Button(root, text="Run Simulation", command=update_parameters_and_run)
run_button.grid(row=2, column=0, padx=10, pady=10)

root.mainloop()