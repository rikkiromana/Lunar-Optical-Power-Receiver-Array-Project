import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# photovoltaic cell parameters -- hange input depending which we will be using
cell_width_cm = 1.35
cell_height_cm = 1.6
beam_center_x_cm = 0
beam_center_y_cm = 0
wavelength_nm = 1065
average_power_w = 32
optical_power_density_w_cm2 = 2
output_nominal_voltage_v = 0.9

# glass cover properties (borosilicate glass similar to BK7?)
glass_refractive_index = 1.51446

# Gaussian light distribution function
def gaussian(x, y, beam_center_x_cm=0, beam_center_y_cm=0):
    beam_waist_cm = cell_width_cm / 2
    sigma = beam_waist_cm / np.sqrt(2 * np.log(2))
    return (1 / (2 * np.pi * sigma**2)) * np.exp(-((x - beam_center_x_cm)**2 + (y - beam_center_y_cm)**2) / (2 * sigma**2))

# adjusted Fresnel coefficient calculation for two reflections at normal incidence
def calculate_fresnel_losses(n1, n2, theta_i=0):
    cos_theta_i = np.cos(np.deg2rad(theta_i))
    sin_theta_t = n1 / n2 * np.sin(np.deg2rad(theta_i))  # Snell's law
    cos_theta_t = np.sqrt(1 - sin_theta_t**2)
    
    rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
    rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
    fresnel_coefficient = (rs + rp) / 2
    
    # applying the coefficient twice for two reflections
    total_loss = fresnel_coefficient ** 2
    return 1 - total_loss  # returning the remaining energy after losses

# light distribution integration with glass cover consideration, including two reflections
def integrate_light_distribution_with_glass(beam_center_x_cm=0, beam_center_y_cm=0):
    result, _ = dblquad(lambda y, x: gaussian(x, y, beam_center_x_cm, beam_center_y_cm) * 
                        calculate_fresnel_losses(1, glass_refractive_index), 
                        -cell_width_cm/2, cell_width_cm/2,
                        lambda x: -cell_height_cm/2, lambda x: cell_height_cm/2)
    return result

# this is a placeholder function for photocurrent calculation based on interpolated I-V curves-- WE WILL WORK ON THIS LATER :)
def calculate_photocurrent(light_intensity):
    voltage = np.array([0, 0.5, 0.9])  # example voltage values
    current = np.array([0, 10, 20])  # example current values for given voltages
    iv_curve = interp1d(voltage, current, kind='linear')
    photocurrent = iv_curve(output_nominal_voltage_v) * light_intensity
    return photocurrent

# function to simulate beam alignment -- WE WILL ALSO WORK ON THIS FUNCITON MORE LATER :)
def simulate_beam_alignment():
    alignment_range = np.linspace(-0.5, 0.5, 20)  # 20 positions along each axis
    max_photocurrent = 0
    optimal_position = (0, 0)
    
    for x_shift in alignment_range:
        for y_shift in alignment_range:
            light_intensity = integrate_light_distribution_with_glass(x_shift, y_shift)
            photocurrent = calculate_photocurrent(light_intensity)
            
            if photocurrent > max_photocurrent:
                max_photocurrent = photocurrent
                optimal_position = (x_shift, y_shift)

    return optimal_position, max_photocurrent

# main function-- to execute calculations and simulations
def main():
    # the initial light intensity and photocurrent calculations
    initial_light_intensity = integrate_light_distribution_with_glass()
    initial_photocurrent = calculate_photocurrent(initial_light_intensity)
    print("Initial Light Intensity on Cell with Glass Cover:", initial_light_intensity)
    print("Initial Photocurrent Generated in Cell:", initial_photocurrent, "A")
    
    # simulating beam alignment to find optimal position
    optimal_position, max_photocurrent = simulate_beam_alignment()
    print(f"Optimal Beam Position: {optimal_position}, Max Photocurrent: {max_photocurrent} A")
    
    # for displaying Gaussian light distribution
    x = np.linspace(-cell_width_cm/2, cell_width_cm/2, 100)
    y = np.linspace(-cell_height_cm/2, cell_height_cm/2, 100)
    X, Y = np.meshgrid(x, y)
    Z = gaussian(X, Y, optimal_position[0], optimal_position[1])
    plt.contourf(X, Y, Z, cmap='viridis')
    plt.colorbar(label='Intensity')
    plt.title('Gaussian Distribution of Light on Photovoltaic Cell at Optimal Beam Position')
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    plt.show()

if __name__ == "__main__":
    main()
