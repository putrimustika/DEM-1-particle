import math
import matplotlib.pyplot as plt

# particle properties
radius                = 0.01
density               = 1000
spring                = 10000
restitution           = 0.9
sliding_friction_coef = 0.3
rolling_friction_coef = 0.01
gravity               = -9.81

# Initial conditions
position = [[0.0, 0.01, 0.0]]
velocity = [[0.1, 0.0, 0.0]]
rotSpeed = [[0.0, 0.0, 0.0]]
rotAngle = [[0.0, 0.0, 0.0]]
displacement_t = [0.0, 0.0, 0.0]
particle = []

# ------------------------
mass    = 4/3 * math.pi * radius**3 * density
inertia = 2/5 * mass * radius**2
damping = -2 * math.log(restitution) * math.sqrt(mass*spring/(math.pi**2 + (math.log(restitution)**2)))

# simulation properties
dt = 0.000001
iteration = 1000000

# Time step check
def time_step_check():
    suggested_dt = math.pi/10 * math.sqrt(mass/spring)
    if dt > suggested_dt:
        text = 'time step should be <= ' + str(suggested_dt) + ' s.'
        return(print(text))
    else:
        pass

# Gravitational forces
def gravitational():
    return mass * gravity

# Normal contact force
def normal_contact(displacement_n, relative_vel_n):
    if displacement_n <= 0:
        return -spring * displacement_n - damping * relative_vel_n
    elif abs(displacement_n) < 1e-6:  # small gap: still treat as contact
        return 0.0
    else:
        return 0.0

# Tangential contact force
def tangential_contact(F_normal, displacement_t, relative_vel_t):
    tangential_force = []
    for displacement_component, relative_velocity_component in zip(displacement_t,relative_vel_t):
        total_force_component = -spring * displacement_component - damping * relative_velocity_component
        tangential_force.append(total_force_component)
    tangential_magnitude = math.sqrt(sum(f**2 for f in tangential_force))
    max_friction = sliding_friction_coef * abs (F_normal)
    
    if tangential_magnitude <= max_friction:
        return tangential_force, False
    else:
        scale = max_friction / tangential_magnitude
        limited_force = [f * scale for f in tangential_force] 
        return limited_force, True

# Run simulation
time_step_check()

for t in range(iteration):
    # Current state
    pos_now = position[-1]
    vel_now = velocity[-1]
    rotSpeed_now = rotSpeed[-1]

    # Update position
    new_pos = [p + v * dt for p, v in zip(pos_now, vel_now)]
    position.append(new_pos)

    # Compute normal force, assume flat floor at y=0
    displacement_y = new_pos[1] - radius
    F_normal = normal_contact(displacement_y, vel_now[1])

    # Compute relative velocity at contact point
    v_contact = [vel_now[i] + rotSpeed_now[2] * radius if i == 0 else 
                 vel_now[i] - rotSpeed_now[0] * radius if i == 2 else 
                 vel_now[i] for i in range(3)]
    rel_vel_t = [v_contact[i] if i != 1 else 0.0 for i in range(3)]  # Only tangential components (x and z)

    # Compute tengential displacement
    displacement_t = [d + v * dt for d, v in zip(displacement_t, rel_vel_t)]

    # Compute tangential force
    F_tangent, slipping = tangential_contact(F_normal, displacement_t, rel_vel_t)
    if slipping:
        displacement_t = [0.0, 0.0, 0.0]

    # Total force
    F_total = []
    F_gravity = gravitational()
    for i in range(3):
        if i == 1:
            F_total.append(F_tangent[i] + F_gravity)
        else:
            F_total.append(F_tangent[i])
    if F_normal != 0:
        F_total[1] += F_normal

    # Linear velocity update
    new_vel = [vel_now[i] + (F_total[i] / mass) * dt for i in range(3)]
    velocity.append(new_vel)

    # Compute rolling friction
    r_eff = radius
    if any(rotSpeed_now):
        mag = math.sqrt(sum(omega**2 for omega in rotSpeed_now))
        unit_vector = [omega / mag for omega in rotSpeed_now]
        Rolling_friction = [-rolling_friction_coef * F_normal * r_eff * d for d in unit_vector]
    else:
        Rolling_friction = [0.0, 0.0, 0.0]

    # Angular acceleration
    torque_contact = [-F_tangent[2] * radius, 0.0, F_tangent[0] * radius]
    torque = [torque_contact[i] + Rolling_friction[i] for i in range(3)]
    alpha = [torque[i] / inertia for i in range(3)]
    
    new_rotSpeed = [rotSpeed_now[i] + alpha[i] * dt for i in range(3)]
    rotSpeed.append(new_rotSpeed)

    new_rotAngle = [rotAngle[-1][i] + new_rotSpeed[i] * dt for i in range(3)]
    rotAngle.append(new_rotAngle)

# Plot x-direction velocity and rolling speed (ω·R)
time = [i * dt for i in range(len(velocity))]
v_x = [v[0] for v in velocity]
wR = [w[2] * radius for w in rotSpeed]  # Assuming rolling around z-axis

with open("C:/My Files/Events/1 particle project/motion_data.txt", "w") as f:
    f.write("step,time,radius,x,y,z,rx,ry,rz\n")
    for i in range(len(position)):
        x, y, z = position[i]
        rx, ry, rz = rotAngle[i]
        f.write(f"{i},{i*dt},{radius},{x},{y},{z},{rx},{ry},{rz}\n")

plt.figure(figsize=(10, 5))
plt.plot(time, v_x, label='TranslationalVelocity_X (m/s)')
plt.plot(time, wR, label='RollingSpeed_Z (rad/s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('3D Rolling Motion: Velocity vs Rolling Speed')
plt.legend()
plt.grid(True)
plt.show()