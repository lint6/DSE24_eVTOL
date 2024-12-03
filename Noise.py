import numpy as np 
import matplotlib.pyplot as plt 

class WeightBalance:

    def __init__(self, MTOW=718.89, payload=185, airframe=170, batteries=80, fuelcell=120, h2tank=173, motor=75):
        self.MTOW = MTOW 
        
        # weights 
        self.m_payload = payload
        self.m_airframe = airframe
        self.m_batteries = batteries
        self.m_fuelcell = fuelcell
        self.m_h2tank = h2tank 
        self.m_motor = motor 

        # CG percentages
        self.cg_payload = 14.6 
        self.cg_airframe = 29 
        self.cg_batteries = 29 
        self.cg_fuel_cell = 14.6
        self.cg_h2tank = 40 
        self.cg_motor = 36 

        self.balancing()

    def balancing(self):
        self.total_weight = np.sum([self.m_payload, self.m_airframe, 
                              self.m_batteries, self.m_fuelcell, 
                              self.m_h2tank, self.m_motor])
        self.total_multiple = np.sum([self.m_payload * self.cg_payload,
                                self.m_airframe * self.cg_airframe,
                                self.m_batteries * self.cg_batteries,
                                self.m_fuelcell * self.cg_fuel_cell,
                                self.m_h2tank * self.cg_h2tank,
                                self.m_motor * self.cg_motor])
        self.final_cg = self.total_multiple / self.total_weight 

    def visual_cg_bar(self):
        components = ['Payload', 'Airframe', 'Batteries', 'Fuel Cell', 'H2 Tank', 'Motor']
        cg_values = [self.cg_payload, self.cg_airframe, self.cg_batteries, 
                     self.cg_fuel_cell, self.cg_h2tank, self.cg_motor]
        plt.figure(figsize=(10, 6))
        plt.bar(components, cg_values, color='skyblue', label='System CG Values')
        plt.axhline(y=self.final_cg, color='red', linestyle='--', label=f'Final CG = {self.final_cg:.2f}%')
        plt.axhline(y=29, color='green', linestyle = '-', label = 'Ideal CG = 29%')
        plt.title('System CG Values and Final CG', fontsize = 16)
        plt.xlabel('Components', fontsize=14)
        plt.ylabel('CG Position (%)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    def visual_cg_line(self):
        components = ['Payload', 'Airframe', 'Batteries', 'Fuel Cell', 'H2 Tank', 'Motor']
        cg_values = [self.cg_payload, self.cg_airframe, self.cg_batteries, 
                     self.cg_fuel_cell, self.cg_h2tank, self.cg_motor]
        weights = [self.m_payload, self.m_airframe, self.m_batteries, 
                   self.m_fuelcell, self.m_h2tank, self.m_motor]
        plt.figure(figsize=(12, 6))
        plt.plot([0, 100], [0, 0], color='black', linewidth=2)
        plt.scatter(cg_values, np.zeros_like(cg_values), s=np.array(weights) * 2, color='blue', label='Component CG')
        for i, label in enumerate(components):
            plt.text(cg_values[i], 0.05, label, ha='center', fontsize=10)
        plt.scatter(self.final_cg, 0, color='red', s=200, label=f'Final CG = {self.final_cg:.2f}%')
        plt.yticks([])
        plt.xticks(np.arange(0, 101, 10), fontsize=12)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()




if __name__ == '__main__':
    WB = WeightBalance()
    print(f"Final CG: {WB.final_cg:.2f}")
    WB.visual_cg_bar()
    WB.visual_cg_line()



