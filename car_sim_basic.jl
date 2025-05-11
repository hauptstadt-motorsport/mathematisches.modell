using DifferentialEquations
using Plots

# --- Physical Parameters ---
const m = 0.055            # car mass (kg)
const g = 9.81             # gravity (m/s²)
const ρ = 1.225            # air density (kg/m³)
const Cd = 0.32            # drag coefficient
const A = 0.002            # frontal area (m²)
const μ_static = 0.04      # static friction coefficient
const μ_kinetic = 0.015    # kinetic friction coefficient
const r = 0.015            # wheel radius (m)
const I_wheel = 0.5 * 0.002 * r^2  # wheel moment of inertia (approx)
const num_wheels = 4

# --- CO₂ Thrust Curve ---
function thrust(t)
    if t < 0.12
        return 12.5 * exp(-25 * t) + 2.0
    elseif t < 0.25
        return 2.0 * exp(-8 * (t - 0.12))
    else
        return 0.0
    end
end

# --- Reynolds-based Drag Coefficient Adjustment (optional) ---
function dynamic_Cd(v)
    Re = (ρ * v * 0.1) / (1.8e-5)
    return Cd + 0.01 * log1p(Re / 1e4)
end

# --- ODE System ---
# u[1] = v (m/s), u[2] = x (m), u[3] = a (m/s²), u[4] = ω (rad/s), u[5] = E (J), u[6] = F_net (N)
function dynamics!(du, u, p, t)
    v = max(u[1], 0.0)
    x = u[2]
    ω = u[4]

    F_thrust = thrust(t)
    Cd_eff = dynamic_Cd(v)
    F_drag = 0.5 * ρ * Cd_eff * A * v^2
    F_gravity = m * g
    F_friction = (v > 0.01) ? μ_kinetic * F_gravity : μ_static * F_gravity

    F_net = F_thrust - F_drag - F_friction

    # Convert thrust to torque at wheels: T = F * r
    T_wheel = F_thrust * r / num_wheels
    α = T_wheel / I_wheel  # angular acceleration

    a = F_net / m

    # Update derivatives
    du[1] = a                       # dv/dt = a
    du[2] = v                       # dx/dt = v
    du[3] = (a - u[3]) / 0.001      # smoothed jerk
    du[4] = α                       # dω/dt = α
    du[5] = 0.5 * m * v^2           # kinetic energy
    du[6] = F_net                   # force for plotting
end

# --- Initial Conditions ---
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [v, x, a, ω, E, F_net]
tspan = (0.0, 2.0)

prob = ODEProblem(dynamics!, u0, tspan)
sol = solve(prob, Vern9(), reltol=1e-10, abstol=1e-10)

# --- Extract Results ---
v = [u[1] for u in sol.u]
x = [u[2] for u in sol.u]
a = [u[3] for u in sol.u]
ω = [u[4] for u in sol.u]
E = [u[5] for u in sol.u]
F_net = [u[6] for u in sol.u]
t = sol.t

# --- Plotting ---
plot(t, v, label="Velocity (m/s)", lw=2)
plot!(t, x, label="Position (m)", lw=2)
plot!(t, a, label="Acceleration (m/s²)", lw=2)
plot!(t, F_net, label="Net Force (N)", lw=2)
plot!(t, E, label="Kinetic Energy (J)", lw=2)
xlabel!("Time (s)")
title!("Ultra-Complete STEM Racing Car Simulation")
legend!(:bottomright)
