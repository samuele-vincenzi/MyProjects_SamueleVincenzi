# CONSTANT TORQUE ATTITUDE PROPAGATION
# Task provided by EnduroSat

# %%  IMPORT LIBRARIES
from Utility_functions.QuaternionUtilities import quaternion_mult
from Utility_functions.AttitudeConversion import attitude_conversion
from Utility_functions.AttitudeProp import attitude_propagator
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"


# %% INITIALIZE THE PROBLEM

# Inertia parameters
I = np.array([
    [1.0, 0.1, 0.1],
    [0.1, 2.0, 0.1],
    [0.1, 0.1, 0.3]
])  # [kg m^2]

# Initial quaternion (inertial to body frame)
q_sensor = (0.936680, 0.216195, 0.124597, 0.245695)
q_sensor /= np.linalg.norm(q_sensor)  # Normalize
# Rotation matrix from sensor to body frame
R = np.array([
    [0, 0, -1],
    [0, -1, 0],
    [-1, 0, 0]
])
# quaternion from sensor to body frame
q_rot = attitude_conversion(R, 'DCM_2_quat')
q_in = quaternion_mult(q_sensor, q_rot)

# Initial angular velocity vector
omega_in = np.radians([1, 2, 0]).T  # [deg/s]

# Initial state
state_in = (*q_in, *omega_in)

# External torque
M = (0.001, 0.0015, 0.0)

# Time window for propagation
t_in = 0  # [s]
t_fin = 1000  # [s]
dt = 0.01  # [s]
n_s = int((t_fin-t_in)/dt)+1
tvec = np.linspace(t_in, t_fin, n_s)

# Numerical integration method
method = 'RK45'

# %% ATTITUDE PROPAGATION

# Propagation function
q_vec, eul_angles, DCM_vec, omega_vec = attitude_propagator(
    I, M, tvec, state_in, method, 'quaternion')

# Compute angular momentum and kinetic energy
h_vec = I@omega_vec
h_norm = np.linalg.norm(h_vec, axis=0)
T_vec = 0.5 * np.einsum('ij,ni,nj->n', I, omega_vec.T, omega_vec.T)

# %% KINETIC ENERGY AND ANGULAR MOMENTUM ELLIPSOIDS

# Disable interactive rotation
plt.ioff()

# Compute the initial kinetic energy and angular momentum
T0 = 0.5 * omega_in @ (I @ omega_in)
h0_norm = np.linalg.norm(I @ omega_in)

# Diagonalize inertia tensor
eigvals, eigvecs = np.linalg.eigh(I)
if np.linalg.det(eigvecs) < 0:
    eigvecs[:, 0] *= -1  # ensure right-handed basis

# Parametric variables
u = np.linspace(0, 2*np.pi, 120)
v = np.linspace(0, np.pi, 60)
u, v = np.meshgrid(u, v)

# Kinetic Energy Ellipsoid in ω-space (principal frame)
rx, ry, rz = np.sqrt(2*T0 / eigvals)
x = rx * np.cos(u) * np.sin(v)
y = ry * np.sin(u) * np.sin(v)
z = rz * np.cos(v)

# Kinetic Energy Ellipsoid in ω-space (body frame)
omega_body = eigvecs @ np.array([x.flatten(), y.flatten(), z.flatten()])
X, Y, Z = [a.reshape(x.shape) for a in omega_body]

# Map to L-space
H_ellipsoid = I @ np.array([X.flatten(), Y.flatten(), Z.flatten()])
H_X, H_Y, H_Z = [a.reshape(x.shape) for a in H_ellipsoid]

# Angular Momentum Sphere in L-space (principal frame)
Hx_diag = h0_norm * np.sin(v) * np.cos(u)
Hy_diag = h0_norm * np.sin(v) * np.sin(u)
Hz_diag = h0_norm * np.cos(v)

# Angular Momentum Sphere in L-space (body frame)
H_body = eigvecs @ np.array([Hx_diag.flatten(),
                            Hy_diag.flatten(), Hz_diag.flatten()])
Hx, Hy, Hz = [a.reshape(x.shape) for a in H_body]

# Create figure
fig = go.Figure()

# Plot the kinetic energy ellipsoid
fig.add_trace(go.Mesh3d(
    x=H_X.flatten(),
    y=H_Y.flatten(),
    z=H_Z.flatten(),
    opacity=1,       # semi-transparent
    color='blue',
    alphahull=0,
    name='Kinetic Energy Ellipsoid',
    showlegend=True
))

# Plot the angular momentum sphere
fig.add_trace(go.Mesh3d(
    x=Hx.flatten(),
    y=Hy.flatten(),
    z=Hz.flatten(),
    opacity=1,       # semi-transparent
    color='green',
    alphahull=0,
    name='Angular Momentum Sphere',
    showlegend=True
))

# Add time evolution of the angular momentum
fig.add_trace(go.Scatter3d(
    x=h_vec[0, :],
    y=h_vec[1, :],
    z=h_vec[2, :],
    mode='lines',
    line=dict(color='red', width=3),
    name='Propagated Angular Momentum (h)',
    showlegend=True
))

# Update layout
fig.update_layout(
    title=dict(
        # text="Kinetic Energy Ellipsoid & Angular Momentum Sphere",
        # x=0.5,          # center horizontally
        # y=0.95,         # position above the plot (1 is top of figure)
        # xanchor='center',
        # yanchor='top',
        # font=dict(size=24)  # title font size
    ),
    scene=dict(
        xaxis=dict(title="H₁"),  # Unicode subscript
        yaxis=dict(title="H₂"),
        zaxis=dict(title="H₃"),
        aspectmode='data'  # ensures equal scaling of axes
    ),
    scene_camera=dict(
        eye=dict(x=1.5, y=1.5, z=1.0)  # better 3D perspective
    ),
    legend=dict(
        x=0.50,
        y=0.75,
        font=dict(size=14),
        bgcolor='rgba(255,255,255,0.5)',
        bordercolor='black',
        borderwidth=1
    ),
    showlegend=True
)

fig.show()

# %% FINAL IMAGES GENERATION

# LaTeX-style fonts for better math rendering
plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "legend.fontsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "text.usetex": False,  # set True if you have LaTeX installed
    "figure.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

# Time evolution of the angular velocity
fig, (wx, wy, wz) = plt.subplots(
    1, 3, figsize=(15, 4), sharex=True, sharey=False)
fig.suptitle("Angular Velocity Components", fontsize=25)

labels = [r"$\omega_x$", r"$\omega_y$", r"$\omega_z$"]
axes = [wx, wy, wz]

for ax, comp, label in zip(axes, omega_vec, labels):
    ax.plot(tvec, comp, linewidth=2, label=label)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"Value $\omega$ (rad/s)")
    ax.legend(loc="upper right")
    ax.set_xlim(t_in, t_fin)
    ax.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Time evolution of the angular momentum
fig, (hx, hy, hz) = plt.subplots(
    1, 3, figsize=(15, 4), sharex=True, sharey=False)
fig.suptitle("Angular Momentum Components", fontsize=25)

labels_h = [r"$h_x$", r"$h_y$", r"$h_z$"]
axes_h = [hx, hy, hz]

for ax, comp, label in zip(axes_h, h_vec, labels_h):
    ax.plot(tvec, comp, linewidth=2, label=label, color="tab:orange")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"Value $h$ (N·m·s)")
    ax.legend(loc="upper right")
    ax.set_xlim(t_in, t_fin)
    ax.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Time evolution of the Euler angles
plt.figure(figsize=(10, 6))
plt.plot(tvec, eul_angles[0, :], label='Around z [deg]')
plt.plot(tvec, eul_angles[1, :], label='Around y [deg]')
plt.plot(tvec, eul_angles[2, :], label='Around x [deg]')
plt.xlabel('Time [s]', fontsize=20)
plt.ylabel('Angle [deg]', fontsize=20)
plt.title(
    'Euler Angles (1–2–3 sequence)',
    fontsize=23,
    pad=25,          # push title upward (default is ~6)
    loc='center',    # can be 'left', 'center', or 'right'
)
plt.legend(fontsize=18)
plt.grid(True)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=15)
plt.tight_layout()

# Time evolution of the quaternion components
fig, axes = plt.subplots(2, 2, figsize=(16, 6), sharex=True)
fig.suptitle("Quaternion Components", fontsize=30)

labels_quat = [r"$q_x$", r"$q_y$", r"$q_z$", r"$q_w$"]

# Flatten axes array for easy iteration
axes_flat = axes.flatten()

for ax, comp, label in zip(axes_flat, q_vec, labels_quat):
    ax.plot(tvec, comp, linewidth=2, label=label, color="tab:red")
    ax.set_xlabel("Time (s)", fontsize=25)
    ax.set_ylabel("Value (-)", fontsize=25)

    # Axis tick label size
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)

    ax.legend(loc="upper right", fontsize=20)
    ax.set_xlim(t_in, t_fin)
    ax.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()


# %% Demonstrate conservation of kinetic energy and norm of angular momentum
perc_diff_T = (max(T_vec)-min(T_vec))/max(T_vec)*100
perc_diff_h = (max(h_norm)-min(h_norm))/max(h_norm)*100
print(f"Percentage difference of the kinetic energy: {perc_diff_T:.2e} %")
print(
    f"Percentage difference of the angular momentum norm: {perc_diff_h:.2e} %")

# %%
