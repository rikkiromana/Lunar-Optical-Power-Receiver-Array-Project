import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

# Paramètres pour le faisceau gaussien, la cellule photovoltaïque et les pertes de Fresnel
# Parameters for the Gaussian beam, photovoltaic cell and Fresnel losses
cell_width_cm = 2
cell_height_cm = 2
beam_center_x_cm = 0
beam_center_y_cm = 0
beam_sigma_cm = 0.5
voltage = 1.5  # Tension appliquée à la cellule photovoltaïque # voltage applied to the photovoltaic cell in volts
n_air = 1  # Indice de réfraction de l'air # Air refraction index
n_glass = 1.5  # Indice de réfraction du verre # Refractive index of glass

# Fonction de distribution lumineuse gaussienne
def light_distribution(x, y, beam_center_x_cm=beam_center_x_cm, beam_center_y_cm=beam_center_y_cm, beam_sigma_cm=beam_sigma_cm):
    return np.exp(-((x - beam_center_x_cm)**2 + (y - beam_center_y_cm)**2) / (2 * beam_sigma_cm**2))

# Intégration de la distribution lumineuse sur la surface de la cellule
def integrate_light_distribution():
    result, _ = dblquad(light_distribution, -cell_width_cm / 2, cell_width_cm / 2, 
                        lambda x: -cell_height_cm / 2, lambda x: cell_height_cm / 2)
    return result

# Calculation of Fresnel losses for two reflections on the glass protector
# Calcul des pertes de Fresnel pour deux réflexions sur le protecteur en verre
def fresnel_losses(n1=n_air, n2=n_glass):
    R = ((n1 - n2) / (n1 + n2)) ** 2
    T_one_surface = 1 - R
    T_total = T_one_surface ** 2  # Transmission after two thoughts # Transmission après deux réflexions
    return T_total

# Calcul du photocourant (à affiner avec un modèle plus précis)
# Calculation of the photocurrent (to be refined with a more precise model)
def calculate_photocurrent(light_intensity, voltage=voltage):
    return light_intensity * voltage

# Fonction principale pour simuler, incluant les pertes de Fresnel
def main():
    T_total = fresnel_losses()
    
    # Calculation of light intensity without and with Fresnel losses
    # Calcul de l'intensité lumineuse sans et avec les pertes de Fresnel
    light_intensity_without_losses = integrate_light_distribution()
    light_intensity_with_losses = light_intensity_without_losses * T_total
    
    # Calcul du photocourant avec l'intensité lumineuse ajustée pour les pertes de Fresnel
    photocurrent = calculate_photocurrent(light_intensity_with_losses)
    print(f"Average light intensity on the cell (adjusted for Fresnel losses): {light_intensity_with_losses}")
    print(f"Photocurrent generated in the cell: {photocurrent} A")

    # Visualisation de la distribution de la lumière sur la cellule photovoltaïque, ajustée pour les pertes de Fresnel
    x = np.linspace(-cell_width_cm / 2, cell_width_cm / 2, 100)
    y = np.linspace(-cell_height_cm / 2, cell_height_cm / 2, 100)
    X, Y = np.meshgrid(x, y)
    Z = light_distribution(X, Y) * T_total  # Ajustement de la distribution lumineuse pour la visualisation
    plt.contourf(X, Y, Z, cmap='viridis')
    plt.colorbar(label='Intensité lumineuse (ajustée pour les pertes de Fresnel)')
    plt.title('Distribution de la lumière sur la cellule photovoltaïque (avec pertes de Fresnel)')
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    plt.show()

if __name__ == "__main__":
    main()
