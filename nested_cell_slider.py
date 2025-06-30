import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve
import warnings
import tkinter as tk
from tkinter import Scale
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def L_matrix_4x4(f):
    """
    Returns the 4x4 lens matrix for focal length f, acting on [x, x', y, y'].
    """
    L_2x2 = np.array([[1, 0], [-1/f, 1]])
    L_4x4 = np.block([
        [L_2x2, np.zeros((2, 2))],
        [np.zeros((2, 2)), L_2x2]
    ])
    return L_4x4

def S_matrix_4x4(d):
    """
    Returns the 4x4 space (propagation) matrix for distance d, acting on [x, x', y, y'].
    """
    S_2x2 = np.array([[1, d], [0, 1]])
    S_4x4 = np.block([
        [S_2x2, np.zeros((2, 2))],
        [np.zeros((2, 2)), S_2x2]
    ])
    return S_4x4

def solve_system(d, R2, n, R1):
    """
    Solve the system of trigonometric equations for alpha and B/A
    
    Parameters:
    d, R2, n, theta: float values for the system parameters
    
    Returns:
    tuple: (alpha, B_over_A) or None if no solution found
    """
    theta = np.arccos(1- d/R1)
    def equations(vars):
        alpha, BA = vars
        
        # Pre-calculate common terms
        cos_n_half = np.cos((n - 0.5) * theta)
        sin_n_half = np.sin((n - 0.5) * theta)
        cos_n_3half = np.cos((n - 1.5) * theta)
        sin_n_3half = np.sin((n - 1.5) * theta)
        
        # Equation 1
        lhs1 = -cos_n_3half + 2*(1 - d/R2)*cos_n_half
        rhs1 = cos_n_half*np.cos(alpha) - BA*sin_n_half*np.sin(alpha)
        eq1 = lhs1 - rhs1
        
        # Equation 2  
        lhs2 = sin_n_3half - 2*(1 - d/R2)*sin_n_half
        rhs2 = (1/BA)*cos_n_half*np.sin(alpha) + sin_n_half*np.cos(alpha)
        eq2 = lhs2 - rhs2
        
        return [eq1, eq2]
    
    # Try multiple initial guesses to find solutions
    initial_guesses = [
        [2, 0.8],
        [np.pi - 0.01, 0.5], 
        [np.pi/2, 2.0],
        [np.pi/4, 0.3],
    ]
    
    for guess in initial_guesses:
        try:
            solution = fsolve(equations, guess, xtol=1e-12)
            alpha, BA = solution
            alpha = np.abs(np.pi - alpha % (2*np.pi))
            #BA = np.abs(BA)
            # Verify the solution by checking residuals
            residuals = equations(solution)
            if all(abs(r) < 1e-10 for r in residuals):
                return alpha, BA
                
        except:
            continue
    
    return None


R1 = 4000
R2 = 3355
hole_dist = 43.588
wavelength = 5.32e-4
w0 = 2.1 
f1 = R1/2
f2 = R2/2 
d = 535.89  # Initial d value
num_points = 366  # Initial number of points

class OpticalSimulationApp:
    def __init__(self, master):
        self.master = master
        master.title("Optical Simulation with Sliders")

        # Create sliders frame
        control_frame = tk.Frame(master)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        # Slider initialization (unchanged)
        self.d_slider = Scale(control_frame, from_=500, to=650, resolution=0.1,
                             orient='horizontal', label='d (distance)')
        self.d_slider.set(d)
        self.d_slider.pack(fill=tk.X)

        self.num_points_slider = Scale(control_frame, from_=0, to=700, resolution=1,
                                      orient='horizontal', label='Number of Points')
        self.num_points_slider.set(num_points)
        self.num_points_slider.pack(fill=tk.X)

        self.R1_slider = Scale(control_frame, from_=3950, to=4050, resolution=1,
                                      orient='horizontal', label='R1')
        self.R1_slider.set(R1)
        self.R1_slider.pack(fill=tk.X)

        self.R2_slider = Scale(control_frame, from_=3200, to=3500, resolution=1,
                                      orient='horizontal', label='R2')
        self.R2_slider.set(R2)
        self.R2_slider.pack(fill=tk.X)

        # Create figure with TWO subplots
        self.fig, (self.ax_blue, self.ax_red) = plt.subplots(1, 2, figsize=(16, 10))
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Initialize persistent plot objects
        self.scatter_blue = self.ax_blue.scatter([], [], color='blue', alpha=1, label='mirror 1', s=[])
        self.scatter_red = self.ax_red.scatter([], [], color='red', alpha=1, label='mirror 2', s=[])
        self.scatter_hole = self.ax_blue.scatter([], [], color='green', alpha=1, label='entry point', s=[])
        
        # Store circle references
        self.circles_blue = []
        self.circles_red = []

        # Label for info text
        self.info_label = tk.Label(master, text="", font=("Arial", 12), justify="left", anchor="w")
        self.info_label.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        # Bind sliders
        self.d_slider.config(command=lambda _: self.update_plot())
        self.num_points_slider.config(command=lambda _: self.update_plot())
        self.R1_slider.config(command=lambda _: self.update_plot())
        self.R2_slider.config(command=lambda _: self.update_plot())

        # Initial plot
        self.update_plot()

    def update_plot(self):
        # Get slider values (unchanged)
        current_d = self.d_slider.get()
        current_num_points = self.num_points_slider.get()
        current_R1 = self.R1_slider.get()
        current_R2 = self.R2_slider.get()

        # Simulation code (unchanged)
        theta = np.arccos(1 - current_d/current_R1)
        B_over_A = solve_system(current_d, current_R2, 3, current_R1)[1]
        A = np.sqrt(hole_dist**2/(np.cos(theta/2)**2 + B_over_A**2 * np.sin(theta/2)**2))
        B = B_over_A * A
        q = 1j * np.pi * w0**2 / wavelength

        x0 = A * np.cos(theta/2)
        y0 = B * np.sin(theta/2)
        x1 = x0 
        y1 = -y0 
        x0_grad = (x1 - x0)/current_d 
        y0_grad = (y1 - y0)/current_d 
        inner_mirror_radius = (np.sqrt(A**2 * np.sin(theta/2)**2 + B**2 * np.cos(theta/2)**2) + 
                              np.sqrt(A**2 * np.sin(3*theta/2)**2 + B**2 * np.cos(3*theta/2)**2))/2

        states = np.zeros((4, current_num_points + 1))
        states[:, 0] = [x0, x0_grad, y0, y0_grad]
        spotsize = np.zeros(current_num_points + 1)
        spotsize[0] = w0

        for i in range(1, current_num_points + 1):
            states[:,i] = S_matrix_4x4(current_d) @ states[:, i-1]
            q = q + np.linalg.norm([states[0, i] - states[0, i-1], states[2, i] - states[2, i-1], current_d])
            q_inv = 1/q 
            spotsize[i] = np.sqrt(-wavelength / (np.pi * q_inv.imag))
            if states[0, i]**2 + states[2, i]**2 > inner_mirror_radius ** 2:
                states[:,i] = L_matrix_4x4(current_R1/2) @ states[:, i]
                q = q / (1 - (2/current_R1)*q)
            else:
                states[:,i] = L_matrix_4x4(current_R2/2) @ states[:, i]
                q = q / (1 - (2/current_R2)*q)

        x_points = states[0, :] 
        y_points = states[2, :]
        radii = [hole_dist]

        # Update scatter plots (FAST update)
        blue_points = np.column_stack((x_points[2:current_num_points+1:2], y_points[2:current_num_points+1:2]))
        red_points = np.column_stack((x_points[1:current_num_points+1:2], y_points[1:current_num_points+1:2]))
        
        self.scatter_blue.set_offsets(blue_points)
        self.scatter_blue.set_sizes(7 * np.pi * spotsize[2:current_num_points+1:2]**2)
        
        self.scatter_red.set_offsets(red_points)
        self.scatter_red.set_sizes(7 * np.pi * spotsize[1:current_num_points+1:2]**2)
        
        self.scatter_hole.set_offsets([[x_points[0], y_points[0]]])
        self.scatter_hole.set_sizes([7 * np.pi * spotsize[0]**2])

        # Update circles (remove old, add new)
        for circle in self.circles_blue + self.circles_red:
            circle.remove()
        self.circles_blue = []
        self.circles_red = []
        
        for radius in radii:
            circle_blue = plt.Circle((0, 0), radius, color='green', fill=False, linestyle='--', alpha=0.7)
            circle_red = plt.Circle((0, 0), radius, color='green', fill=False, linestyle='--', alpha=0.7)
            self.ax_blue.add_patch(circle_blue)
            self.ax_red.add_patch(circle_red)
            self.circles_blue.append(circle_blue)
            self.circles_red.append(circle_red)

        # Update titles and labels
        self.ax_blue.set_title(f"R1 = {current_R1}, d={current_d:.2f}, N={current_num_points}")
        self.ax_red.set_title(f"R2 = {current_R2}, d={current_d:.2f}")
        
        self.ax_blue.axis('equal')
        self.ax_red.axis('equal')
        
        # Update info text
        entry_state = states[:, 0]
        exit_state = states[:, -1]
        info_text = (
            f"B/A = {B_over_A:.5f}\n"
            f"Entry state = [{entry_state[0]:.3f}, {entry_state[1]:.3f}, {entry_state[2]:.3f}, {entry_state[3]:.3f}]\n"
            f"Exit state = [{exit_state[0]:.3f}, {exit_state[1]:.3f}, {exit_state[2]:.3f}, {exit_state[3]:.3f}]\n"
            f"Time in cavity = {1000* current_d * current_num_points / 299793 }ns \n"
            f"final spot size = {spotsize[-1]}"
        )
        self.info_label.config(text=info_text)

        # Only redraw canvas (no full recreation)
        self.canvas.draw_idle()

root = tk.Tk()
app = OpticalSimulationApp(root)
root.mainloop()