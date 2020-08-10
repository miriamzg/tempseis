import matplotlib.pylab as plt
import sys
import numpy as np
import matplotlib as mpl
import math


def rotate_vector_z(v, a_deg):
    a = np.deg2rad(a_deg)

    R = [[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]]

    v_rot = np.dot(v, R)
    return v_rot


def rotate_axis(x1, x2, angle):
    angle = math.radians(angle)
    x1_new = x1 * np.cos(angle) - x2 * np.sin(angle)
    x2_new = x1 * np.sin(angle) + x2 * np.cos(angle)
    return x1_new, x2_new


def full_rotation(x, y, z, strike, dip):
    x_new1, y_new1 = rotate_axis(x, y, strike + 90)
    z_new1 = z
    y_new2, z_new2 = rotate_axis(y_new1, z_new1, -dip + 90)
    x_new2 = x_new1
    return x_new2, y_new2, z_new2


result_folder = sys.argv[1]
result_file = result_folder + "/rfi_models"


fault_model = "./c000kn4n.fsp"
result_parameters = result_folder + "/best_parameters.txt"


ev_lon = 144.66
ev_lat = 37.17
ev_dep = 24.9

strike = 14.0
dip = 49.0


# load parameter results
lines = open(result_parameters).readlines()
rc_mod = float(lines[3].split()[2])
rc_ang = float(lines[4].split()[2])
Amax = float(lines[6].split()[1])
Amax_min = float(lines[6].split()[3])
Amax_max = float(lines[6].split()[5])
Amin = float(lines[7].split()[1])
Amin_min = float(lines[7].split()[3])
Amin_max = float(lines[7].split()[5])
phi = float(lines[8].split()[1])
phi_min = float(lines[8].split()[3])
phi_max = float(lines[8].split()[5])
v_abs = float(lines[9].split()[2])
v_ang = float(lines[10].split()[2])

# location of the point source used in the fault plane ref. system
r0_x = 0
r0_y = 0


v_abs_t = 4.0
v_ang_t = 0.0


vf = rotate_vector_z([v_abs, 0, 0], v_ang)
vf_t = rotate_vector_z([v_abs_t, 0, 0], v_ang_t)


rc_ang_r = np.deg2rad(rc_ang)
rc_x = r0_x + rc_mod * np.sin((np.pi / 2.0) - rc_ang_r)
rc_y = r0_y + rc_mod * np.cos((np.pi / 2.0) - rc_ang_r)

# plotting

ell = mpl.patches.Ellipse(
    xy=[rc_x, rc_y],
    width=Amax_min,
    height=Amin_min,
    angle=phi,
    edgecolor="red",
    lw=2,
    facecolor="None",
    linestyle="--",
)
ell2 = mpl.patches.Ellipse(
    xy=[rc_x, rc_y],
    width=Amax_max,
    height=Amin_max,
    angle=phi,
    edgecolor="red",
    lw=2,
    facecolor="None",
    linestyle="--",
)
ell3 = mpl.patches.Ellipse(
    xy=[rc_x, rc_y],
    width=Amax,
    height=Amin,
    angle=phi,
    edgecolor="red",
    lw=2,
    facecolor="None",
    linestyle="-",
)
fig, ax = plt.subplots()


plt.plot([rc_x], [rc_y], marker="o", markersize=5, color="red")
ax.add_patch(ell)
ax.add_patch(ell2)
ax.add_patch(ell3)
ax.set_aspect("equal")
ax.set_xlabel("Distance along strike (km)", fontsize=15)
ax.set_ylabel("Distance along dip (km)", fontsize=15)
ax.arrow(
    rc_x,
    rc_y,
    vf[0],
    vf[1],
    head_width=2,
    head_length=2,
    linewidth=2,
    color="black",
    zorder=2,
)
ax.autoscale()
# plt.show()
plt.savefig(result_folder + "/2D_finite_fault.png")
plt.close()
sys.exit()
