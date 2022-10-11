# %% Imports
import sys
import amici
import amici.plotting
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.optimize import dual_annealing as dh
import pypesto
import pypesto.optimize as optimize
import pypesto.visualize as visualize
sys.path.append('.')# for odes2py
from odes2py import odes2py

# %% Supress stderr
from contextlib import contextmanager, redirect_stderr
from os import devnull

@contextmanager
def silent_errors(): # Note: might no longer work
    with open(devnull,'w') as fnull:
        with redirect_stderr(fnull) as err:
            yield err


#%% Define functions for plotting the data
def plot_data(data): 
    plt.errorbar(data["time"], data["mean"], yerr=data["SEM"], fmt='o', capsize=3)
    plt.xlabel("time (s)")
    plt.ylabel("Response of Rp")


def plot_agreement(data, rdata, model):
    amici.plotting.plotStateTrajectories(rdata, model=model, state_indices=[1])
    plot_data(data)

# %% Determine if plots should be presented or not. 
plot = False

#%% Load the ODE file from MATLAB/SBtoolbox format, convert to yaml
fileName = "M1.txt"
imported_model = odes2py("M1.txt", 'M1.xml','sbml-yaml')
model_name = imported_model["name"]
observables = imported_model["observables"]
x0 = [p[1] for p in imported_model["parameters"]]


# %% Import the sbml file and convert/compile to AMICI
observables_tuple=observables.copy()
observables = {} 
for name, formula in observables_tuple:
    observables[name]={'name': '', 'formula': formula}
print(f"Observables: {observables}")
sbml_importer = amici.SbmlImporter(model_name+".xml")
sbml_importer.sbml2amici(model_name, model_name+"_amici", observables=observables, verbose=0)


# %% Import the AMICI model
model_module = amici.import_model_module(model_name, model_name+"_amici")
model = model_module.getModel() 


# %% Define the experimental data
with open("data.json",'r') as f:
    data = json.load(f)

edata = amici.ExpData(1,0,0,data["time"]) #specifies the size of the experimental data
edata.setObservedData(data["mean"])
edata.setObservedDataStdDev(data["SEM"])


#%% Simulate using the default parameters (should be bad). 
solver = model.getSolver()
model.setTimepoints(np.linspace(0, 2, 201)) 
rdata = amici.runAmiciSimulation(model, solver, edata)
print(f"Cost using the initial guess: {rdata['chi2']} (should be ~701)")
if plot:
  plot_agreement(data, rdata, model)


# %% Simulate with a set of known "optimal" parameters
with open("M1(13.316).json", "r") as f: 
    x_opt = json.load(f)

model.setParameters(x_opt)
rdata = amici.runAmiciSimulation(model, solver, edata)
print(f"Cost using optimal parameters: {rdata['chi2']} (should be around ~13.3)")
if plot:
  plot_agreement(data, rdata, model)


# %% Setup the model for optimization
model = model_module.getModel()
model.setTimepoints(data["time"])
x0 = np.array(model.getParameters())
solver = model.getSolver()

model.requireSensitivitiesForAllParameters()
solver.setSensitivityMethod(amici.SensitivityMethod_adjoint)
solver.setSensitivityOrder(amici.SensitivityOrder_first)

rdata = amici.runAmiciSimulation(model, solver, edata)
assert(np.floor(rdata["chi2"])==701, "The original cost does not correspond to the known cost of the initial parameter set")


#%% Optimization settings
objective = pypesto.AmiciObjective(model, solver, [edata])
options = {"disp": True}
optimizer = optimize.ScipyOptimizer(options=options)
# optimizer = optimize.PyswarmOptimizer()
# optimizer = optimize.IpoptOptimizer()
optimizer = optimize.FidesOptimizer()

lb = np.tile(1e-6, (1, len(model.getParameters())))
ub = np.tile(1e7, (1, len(model.getParameters())))
scales = ['log10']*5
problem = pypesto.Problem(objective=objective,
                        lb=lb, ub=ub, x_guesses=[x0], x_scales=scales)


# %% Optimize one time using the defined optimizer
with silent_errors():
    opts=pypesto.optimize.OptimizeOptions(allow_failed_starts=False)
    result = optimizer.minimize(problem, x0, "test", optimize_options=opts)


# %% Print the single run optimization results
print(result)
x = result["x"]
model.setParameters(x)
rdata = amici.runAmiciSimulation(model, solver, edata)
print(f"optimized cost: {rdata['chi2']} (optima ~13.3)")
if plot:
  plot_agreement(data, rdata, model)


# %% Multistart optimization, no start guess
with silent_errors():
    result = optimize.minimize(problem,optimizer=optimizer,n_starts=300) # 200
print(result.optimize_result.as_dataframe())
if plot:
    visualize.waterfall(result)
    visualize.parameters(result)


# %%  Print the results of the multistart optimization  
x = result.optimize_result.as_dataframe()["x"].values[0]

model.setParameters(x)
rdata = amici.runAmiciSimulation(model, solver, edata)
print(f"optimized cost: {rdata['chi2']} (optima ~13.3)")
if plot:
    plot_agreement(data, rdata, model)


# %% Optimize using a custom cost function and scipys dual anealing algorithm (unused)
# def cost(param, model, solver, edata):
#   try:
#     model.setParameters(param)
#     rdata = amici.runAmiciSimulation(model, solver, edata)
#     cost = rdata["chi2"]
#   except: 
#     cost=1e90
#   return cost

# iter=0
# def callback(x,f,c):
#     global iter
#     if iter%50==0:
#         print("Iter: {0:4d}, obj:{1:3.6f}".format(iter, f))
#     iter+=1

# bounds = np.array([[1.0e-5, 1.0e5]]*5)

# iter=0
# res = dh(func = cost,bounds = bounds ,args=(model,solver, edata), x0 = x0, callback = callback)

# print(f"Cost from the opimization: {res.fun}")
# model.setParameters(res.x)
# rdata = amici.runAmiciSimulation(model, solver, edata)
# print(f"Cost from the simulation of the optimized parameters: {rdata['chi2']} (should be the same)")
