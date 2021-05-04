import numpy as np
import matplotlib.pyplot as plt
from casadi import *


T = 60
T_retirement = 45

beta = 0.94
gamma = 0.5
eta = 1

little_r = 0.1           
R = 1 + little_r
assets_min = -100
Lmax = 1

tC = 0.05
tL = 0.1
tK = 0

taxes = [tC, tL, tK]


def utility(consumption, labour, gamma, eta):
    # Standard additive CRRA utility function of consumption and labour
    # Input: Vector 1xT consumption, labour; 1x1 gamma, eta
    # Output: utility 1xT
    
    u_consumption = (consumption)**(1 - gamma) / (1 - gamma)
  
    u_labour = (labour)**(1 + eta) / (1 + eta)
    
    u = u_consumption - u_labour
    return u


def wage_fun(T, T_ret):
    # Return a 1xT vector of wages
    # input: T, and retirement age T_ret
    res = np.maximum(1.5,1/2 + np.array(range(1,T+1)) * (1 - np.array(range(1,T+1))/T)) / 16
    res[T_ret:] = 0
    return res

def tax(R, assets, consumption, wage, labour, taxes):
    # Calculate the combined tax on labour income, consumption and capital
    # income

    [tc, tl, tk] = taxes
  
    return tk * (R - 1) * assets + tl * wage * labour + tc * consumption


def asset_constraint(assets, wage, labour, consumption, R, taxes, T):
    temp = (R * assets[:T] + wage * labour[:T] - consumption[:T]
            - tax(R, assets[:T], consumption, wage, labour, taxes)  - assets[1:])
    return temp


consumption = SX.sym('consumption', T, 1)
labour = SX.sym('labour', T, 1)
assets = SX.sym('assets', T + 1, 1)

objective = -sum1(beta**(DM(range(T))) * utility(consumption, labour, gamma, eta))

lower_bound_C = DM.ones(T) * 1e-9    # lower bound on the consumption -> not binding anyway
lower_bound_L = DM.zeros(T)
lower_bound_A = vertcat([0], assets_min * DM.ones(T-1), [0])

upper_bound_C = DM.ones(T) * np.inf
upper_bound_L = vertcat(DM.ones(T_retirement) * Lmax, DM.zeros(T - T_retirement))
upper_bound_A = vertcat([0], DM.ones(T)*np.inf)


lb_x = vertcat(lower_bound_C, lower_bound_L, lower_bound_A)
ub_x = vertcat(upper_bound_C, upper_bound_L, upper_bound_A)



# Get wage profile from wage_fun function
wage = wage_fun(T, T_retirement)

# Define the start point
x_0 = vertcat(DM.ones(T), DM.zeros(T_retirement)+0.5, DM.zeros(T - T_retirement),DM.zeros(T+1))

nonlin_con = asset_constraint(assets, wage, labour, consumption, R, taxes, T)



nlp = {'x':vertcat(consumption,labour,assets), 'f':objective, 'g':nonlin_con}
solver = nlpsol('solver', 'ipopt', nlp)
solution = solver(x0=x_0,lbx=lb_x,ubx=ub_x,lbg=0,ubg=0)
sol = solution['x']


print('consumption', sol[:T])
print('labor', sol[T:2*T])
print('assets', sol[2*T:])

def plot_solution(solution,wage,T):
	plt.figure()
	plt.plot(solution[0:T],'.')
	plt.title('Consumption')
	plt.figure()
	plt.plot(solution[T:2*T],'.')
	plt.title('Labour')
	plt.figure()
	plt.plot(solution[2*T:],'.')
	plt.title('Assets')
	plt.figure()
	plt.plot(wage)
	plt.plot(solution[T:2*T]*wage, ".")
	plt.title("Wage and Earnings")
	plt.show()


plot_solution(sol, wage, T)