# https://trinket.io/glowscript/eff0cc6b52

# global constants
g = 9.81
l = 2.0
dt = 0.01 
max_t = 20.0
initial_angle = 45.0 * pi / 180.0
initial_angle_velocity = 0.0


# simple pendulum (point mass)
def accel(theta):
  return -g / l * sin(theta)


# https://www.physics.udel.edu/~bnikolic/teaching/phys660/numerical_ode/node2.html
def update_state_euler(theta, theta_dot, dt):
  theta_dot_new = theta_dot + accel(theta) * dt
  theta_new = theta + theta_dot * dt
  return (theta_new, theta_dot_new)


# https://www.physics.udel.edu/~bnikolic/teaching/phys660/numerical_ode/node1.html
def update_state_euler_cromer(theta, theta_dot, dt):
  theta_dot_new = theta_dot + accel(theta) * dt
  theta_new = theta + theta_dot_new * dt
  return (theta_new, theta_dot_new)


# https://stackoverflow.com/questions/52985027/runge-kutta-4-and-pendulum-simulation-in-python
def update_state_runge_kutta(y, v, h):
    k1y = h*v
    k1v = h*accel(y)
    k2y = h*(v+0.5*k1v)
    k2v = h*accel(y+0.5*k1y)
    k3y = h*(v+0.5*k2v)
    k3v = h*accel(y+0.5*k2y)
    k4y = h*(v+k3v)
    k4v = h*accel(y+k3y)
    y_new = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0 
    v_new = v + (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
    return (y_new, v_new)

# ----------------------------------------------------------

def plot_approx(dt, max_t, initial_angle):
  approx_curve = gcurve(color=color.red, label="approx")
  t = 0.0
  while t < max_t:
    approx_curve.plot(t, initial_angle * cos(sqrt(g/l) * t))
    t = t + dt


def plot_euler(dt, max_t, initial_angle, initial_angle_velocity):
  euler_curve = gcurve(color=color.green, label="euler")
  theta = initial_angle
  theta_dot = initial_angle_velocity
  t = 0.0
  while t < max_t:
    theta, theta_dot = update_state_euler(theta, theta_dot, dt)
    euler_curve.plot(t, theta)
    t = t + dt


def plot_euler_cromer(dt, max_t, initial_angle, initial_angle_velocity):
  euler_cromer_curve = gcurve(color=color.blue, label="euler-cromer")
  theta = initial_angle
  theta_dot = initial_angle_velocity
  t = 0.0
  while t < max_t:
    theta, theta_dot = update_state_euler_cromer(theta, theta_dot, dt)
    euler_cromer_curve.plot(t, theta)
    t = t + dt


def plot_runge_kutta(dt, max_t, initial_angle, initial_angle_velocity):
  runge_kutta_curve = gcurve(color=color.black, label="runge-kutta")
  theta = initial_angle
  theta_dot = initial_angle_velocity
  t = 0.0
  while t < max_t:
    theta, theta_dot = update_state_euler_cromer(theta, theta_dot, dt)
    runge_kutta_curve.plot(t, theta)
    t = t + dt

#--------------------------------------------------------------

def main():
  graph(xtitle="t [s]", ytitle="theta [rad]")
  plot_approx(dt, max_t, initial_angle)
  plot_euler(dt, max_t, initial_angle, initial_angle_velocity)
  plot_euler_cromer(dt, max_t, initial_angle, initial_angle_velocity)
  # plot_runge_kutta(dt, max_t, initial_angle, initial_angle_velocity)
  
if __name__ == '__main__':
  main()
  