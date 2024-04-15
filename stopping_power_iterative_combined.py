import numpy as np
import matplotlib.pyplot as plt

######################## 2H_in_C data set ########################
energy_levels1 = np.array([
    1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80,
    2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00,
    4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 9.00, 10.00
])

######################## Sum of electronic plus nuclear ########################

stopping_powers1 = np.array([
    0.3552004, 0.3350766, 0.3176564, 0.3022392, 0.2886242,
    0.2764111, 0.2653995, 0.2553892, 0.2462799, 0.2302640,
    0.2142479, 0.1988348, 0.1859239, 0.1747147, 0.1650069,
    0.1564001, 0.14879416, 0.14200892, 0.13028008, 0.12047291,
    0.11226697, 0.10516196, 0.09895768, 0.09355398, 0.08442790,
    0.07714310, 0.07103920
])

######################## 2H_in_CD2 data set ########################
energy_levels2 = np.array([
    1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80,
    2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00,
    4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 9.00, 10.00
])

######################## Sum of electronic plus nuclear ########################

stopping_powers2 = np.array([
    0.4005142, 0.376089, 0.3548678, 0.3362497, 0.320234,
    0.3056202, 0.292708, 0.2811972, 0.2706875, 0.2525708,
    0.2374539, 0.2222402, 0.2074888, 0.1947193, 0.1836111,
    0.173704, 0.16489776, 0.1569929, 0.1434387, 0.132176,
    0.1227941, 0.1146642, 0.107675, 0.101559, 0.091357, 
    0.083458, 0.0764554
])

def calculate_energy_loss(energy_levels, stopping_powers, initial_energy, rho, thickness, dx):
    total_energy_loss = 0.0
    x = 0.0
    x_steps = []
    energy_losses = []
    remaining_energies = []

    current_energy = initial_energy    
    while x < thickness:
        ####### Safety check energies outside data set
        if initial_energy <= energy_levels[0]: 
            S_total = stopping_powers[0]  
        elif initial_energy >= energy_levels[-1]:
            S_total = stopping_powers[-1]  
        else:
            S_total = np.interp(initial_energy, energy_levels, stopping_powers)  # linear interpolation

        ##############################################

        dE = S_total * rho * dx
        total_energy_loss += dE
        initial_energy -= dE
        x += dx

        x_steps.append(x)
        energy_losses.append(total_energy_loss)
        remaining_energies.append(initial_energy)

        if initial_energy <= 0:  # Safety check
            break

    return x_steps, energy_losses, remaining_energies, total_energy_loss

rho1 = 2253  # Density C in mg/cm3
d1 = 0.00019  # Thickness in cm
initial_energy1 = 2.00  # MeV

x_steps1, energy_losses1, remaining_energies1, total_energy_loss_1 = calculate_energy_loss(energy_levels1, stopping_powers1, initial_energy1, rho1, d1, 0.000000001)

rho2 = 1060  # Densityi CD2 in mg/cm3
d2 = 0.00056  # Thickness in cm
initial_energy2 = 2.00 

x_steps2, energy_losses2, remaining_energies2, total_energy_loss_2 = calculate_energy_loss(energy_levels2, stopping_powers2, initial_energy2, rho2, d2, 0.000000001)

# Convert x_steps to microns for better visualization and energy_losses to keV
x_steps1_microns = np.array(x_steps1) * 1e4
energy_losses1_keV = np.array(energy_losses1) * 1e3

x_steps2_microns = np.array(x_steps2) * 1e4
energy_losses2_keV = np.array(energy_losses2) * 1e3

print(f"Total energy loss for {initial_energy1} MeV initial energy over {d1} cm: {total_energy_loss_1 * 1e3} keV")
print(f"Total energy loss for {initial_energy2} MeV initial energy over {d2} cm: {total_energy_loss_2 * 1e3} keV")

# Plotting both datasets
plt.figure(figsize=(10, 6))
plt.plot(x_steps1_microns, energy_losses1_keV, label='Natural Carbon')
plt.plot(x_steps2_microns, energy_losses2_keV, label='Deuterated Polyethylene', linestyle='--')
plt.title('Comparison of Energy Loss vs. Thickness Step')
plt.xlabel('Thickness Step (microns)')
plt.ylabel('Energy Loss (keV)')
plt.legend()
plt.grid(True)
plt.show()

# Plotting both datasets
plt.figure(figsize=(10, 6))
plt.plot(x_steps2_microns, remaining_energies2, label='Deuterated Polyethylene')
plt.plot(x_steps1_microns, remaining_energies1, label='Natural Carbon', linestyle='--')
plt.title('Comparison of Remaining Energy vs. Thickness Step')
plt.xlabel('Thickness Step (microns)')
plt.ylabel('Energy Deuterium (MeV)')
plt.legend()
plt.grid(True)
plt.show()
