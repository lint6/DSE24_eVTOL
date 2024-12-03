# Importing classes
from The_Main import RotorSizing, PowerAnalysis  

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
            "n_blades": 4, 
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
            "A_eq": 0.75, #equivalent flat plate area
            "FM": 0.7 #figure of merit
        }
    
    #conf 2 is quad-rotor
    elif config == 2:
        return {
            "MTOW": 718.89, 
            "n_blades": 4, 
            "n_rotor": 4, 
            "DL": 34.2, # kg/m2
            "bank_angle": 30, #max bank angle
            "cto_fl": 0.11, #next 3 parameters defined depending on max advance ratio
            "cto_turn": 0.14, 
            "cto_turb": 0.17, 
            "co_ax": 1, # 1 if not coaxial, 2if coaxial 
            "d_fact": 0.05, 
            "max_v": 50,  #m/s
            "k_int": 1,
            "A_eq": 0.48, #equivalent flat plate area
            "FM": 0.75 #figure of merit
        }
    
    #conf 3 is compound
    elif config == 3:
        return {
            "MTOW": 718.89, 
            "n_blades": 4, 
            "n_rotor": 6, 
            "DL": 100, #kg/m2
            "bank_angle": 30, #max bank angle
            "cto_fl": 0.11, #next 3 parameters defined depending on max advance ratio
            "cto_turn": 0.14, 
            "cto_turb": 0.17, 
            "co_ax": 1, # 1 if not coaxial, 2if coaxial 
            "d_fact": 0.06, 
            "max_v": 50,  #m/s
            "k_int": 1,
            "A_eq": 0.5, #equivalent flat plate area
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
        d_fact=config["d_fact"], 
        max_v=config["max_v"],
        k_int=config["k_int"],
        A_eq=config["A_eq"],
        FM=config["FM"]
    )
    
    print('----------------------------------------')
    print(f'Number of Blades = {rotor.n_blades} and MTOW = {rotor.MTOW}kg')
    print('----------------------------------------')
    rotor.display_parameters()
    rotor.visual_blade_vs_aspect_ratio()

    # Step 3: Initialize PowerAnalysis with the configured RotorSizing object
    power = PowerAnalysis(rotorsizing=rotor)  # Pass the rotor object directly to PowerAnalysis
    power.display_parameters()

    # Update power analysis with specific conditions if required
    power.iterate_design(new_rho=1.19011)  # Example: Change air density for cruise altitude
    power.iterate_design(new_gamma_CD=0)  # Example: Change climb/descent angle

    # Compute and display hover and forward flight power
    power.iterate_design(new_ROC_VCD=0)
    print(f"HOGE power = {power.P_hoge / 1000:.2f} [kW]")
    print(f"Vertical climb/descent power = {power.P_VCD / 1000:.2f} [kW]")

    # Step 4: Plot power components
    power.plot_power_components()

if __name__ == '__main__':
    main()
