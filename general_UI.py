# Importing classes
from The_Main import RotorSizing, PowerAnalysis, SoundAnalysis, EnergyAnalysis

def get_configuration():
#user chosen configuration
    print("Choose a VTOL configuration:")
    print("1. Configuration 1: Co - Axial")
    print("2. Configuration 2: Quad - rotor")
    print("3. Configuration 3: Compound")
    config = int(input("Enter the configuration number (1, 2, or 3): "))
    
    #conf 1 is co-axial
    if config == 1:
        return {
            "MTOW": 718.89, 
            "n_blades": 2, 
            "n_rotor": 1, 
            "DL": 2.28 * ((718.89)**(1/3) - 2.34), #kg/m2
            "bank_angle": 30, #max bank angle
            "cto_fl": 0.11, #next 3 parameters defined depending on max advance ratio. same for all configurations
            "cto_turn": 0.14, 
            "cto_turb": 0.17,
            "co_ax": 2, # 1 if not coaxial, 2if coaxial 
            "d_fact": 0.04, 
            "max_v": 50,  #m/s
            "k_int": 1.28, #some coaxial thing
            "A_eq": 0.524, #equivalent flat plate area, assumed clean helicopter config
            "FM": 0.7 #figure of merit
        }
    
    #conf 2 is quad-rotor
    elif config == 2:
        return {
            "MTOW": 718.89, 
            "n_blades": 3, 
            "n_rotor": 4, 
            "DL": 14, # kg/m2
            "bank_angle": 30, #max bank angle
            "cto_fl": 0.11, #next 3 parameters defined depending on max advance ratio
            "cto_turn": 0.14, 
            "cto_turb": 0.17, 
            "co_ax": 1, # 1 if not coaxial, 2if coaxial 
            "d_fact": 0.05, 
            "max_v": 50,  #m/s
            "k_int": 1,
            "A_eq": 0.75, #equivalent flat plate area, assumed between clean and utility helicopter config
            "FM": 0.65 #figure of merit, low for multirotors
        }
    
    #conf 3 is compound
    elif config == 3:
        return {
            "MTOW": 718.89, 
            "n_blades": 5, 
            "n_rotor": 4, 
            "DL": 245.32, #kg/m2
            "bank_angle": 30, #max bank angle
            "cto_fl": 0.11, #next 3 parameters defined depending on max advance ratio
            "cto_turn": 0.14, 
            "cto_turb": 0.17, 
            "co_ax": 1, # 1 if not coaxial, 2if coaxial 
            "d_fact": 0.06, #assumption
            "max_v": 57,  #m/s
            "k_int": 1,
            "A_eq": 0.92, #equivalent flat plate area, utility helicopter assumed
            "FM": 0.6 #figure of merit
        }
    else:
        print("Invalid selection! Defaulting to Configuration 1.")
        return print('you wrong boi, try again')

def main():
    """
    Main function to execute the VTOL configuration and power analysis.
    """
    # Step 1: Prompt user for configuration
    config = get_configuration()

    # Step 2: Initialize RotorSizing with user-selected configuration
    rotor = RotorSizing(
        MTOW=config["MTOW"], 
        n_blades=config["n_blades"], 
        n_rotor=config["n_rotor"], 
        DL=config["DL"],
        bank_angle=config["bank_angle"], 
        cto_fl=config["cto_fl"], 
        cto_turn=config["cto_turn"], 
        cto_turb=config["cto_turb"], 
        co_ax=config["co_ax"],
        d_fact=config["d_fact"], 
        max_v=config["max_v"],
        k_int=config["k_int"],
        A_eq=config["A_eq"],
        FM=config["FM"]
    )
    
    print('----------------------------------------')
    print(f'Number of Blades = {rotor.n_blades}, Number of Rotors = {rotor.N_rotors} and MTOW = {rotor.MTOW}kg')
    print('----------------------------------------')
    rotor.display_parameters()

    # uncomment to plot blades vs aspect ratio graph
    rotor.visual_blade_vs_aspect_ratio()

    # Step 3: Initialize PowerAnalysis with the configured RotorSizing object
    power = PowerAnalysis(rotorsizing=rotor)  # Pass the rotor object directly to PowerAnalysis
    power.display_parameters()

    # Update power analysis with specific conditions if required
    power.iterate_design(new_rho=1.19011)  # Example: Change air density for cruise altitude
    power.iterate_design(new_gamma_CD=0)  # Example: Change climb/descent angle

    # Compute and display hover and forward flight power
    power.iterate_design(new_ROC_VCD=0)
    #print(f"HOGE power = {power.P_hoge / 1000:.2f} [kW]")
    #print(f"Vertical climb/descent power = {power.P_VCD / 1000:.2f} [kW]")

    # Step 4: Plot power components
    power.plot_power_components()
    power.final_power()
    power.get_highest_power()
    print('----------------------------------------')
    power.print_all_powers()
    print('----------------------------------------')
    # Print the noise outputs
    sound = SoundAnalysis()
    print('----------Rotational noise----------')
    sound.display_parameters_rotor()
    print('----------Vortex noise----------')
    sound.display_paramenters_vortex()


    # Instantiate the EnergyAnalysis class
    energy_analysis = EnergyAnalysis(power)

    # Calculate mission phase times
    times = energy_analysis.calculate_missionphase_time()

    # Calculate energies
    energies = energy_analysis.calculate_energy_required()

    # Calculate amps
    amps = energy_analysis.calculate_amps()

    print("Mission Phase Times:")
    print(times)
    total_time = times['total']/60
    print(f'Total Mission time = {total_time:.2f} [min]')

    print("\nMission Energies (Wh):")
    print(energies)
    total_energy = energies['total']
    loiter_energy = energies['Loiter']
    loiter_time = times['Loiter']/60
    print(f'Total Energy Consumption = {total_energy:.2f} [Wh]')

    print(f'Loiter power required = {power.min_power:.2f} [kW] at a speed of {power.min_power_velocity*3.6:.2f} [km/h]')
    print(f'Loiter time = {loiter_time} [min]')
    print(f'Energy Consumption during loiter = {loiter_energy:.2f} [Wh]')

    print("\nMission Amps:")
    #print(amps)
    max_amps = amps['max']
    print(f'Max amps = {max_amps:.2f} [A]')

    # plot the PEMFC power vs mission phase/time
    energy_analysis.visual_PEMFC_power()
    energy_analysis.visual_PEMFC_energy()
    
if __name__ == '__main__':
    main()

