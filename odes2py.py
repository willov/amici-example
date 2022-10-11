# Code written by William Lövfors.

from logging import error
import sys
import re

def import_odes(filename, do_print=False):
    """This function imports a model from the [IQM tools / SB toolbox] format. 

    Examples:
        Import a model from the [IQM tools / SB toolbox] format with file name 'model.txt'
            import_odes('model.txt')
        Import a model with file name 'model.txt', and print the imported model
            import_odes('model.txt')

    Args: 
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .py extension.
        do_print: 
            Set to True if you want the imported structure to be printed after importing. 
    
    Returns: 
         A dict containing the fields ["name", "states", "parameters", "variables", "observables", "reactions", events"]
    """

    modelName=re.search('([\w,-]+)\.',filename).group(1)
    states={}
    params=[]
    variables=[]
    reactions=[]
    events=[]
    observables=[]
    with open(filename) as f:
        for line in f:
            if len(line.strip())>1 and not line.strip()[0]=='%':
                line=line.rstrip()
                if line[0] == '*': #can most likely be updated to work with our python toolbox format as well
                    # switch list
                    if re.search('model name', line, re.IGNORECASE):
                        inputType='name'
                    elif re.search('model states', line, re.IGNORECASE):
                        inputType='states'
                    elif re.search('model parameters', line, re.IGNORECASE):
                        inputType='parameters'
                    elif re.search('model variables', line, re.IGNORECASE):
                        inputType='variables'
                    elif re.search('model reactions', line, re.IGNORECASE):
                        inputType='reactions'
                    elif re.search('model event', line, re.IGNORECASE):
                        inputType='events'
                    else:
                        inputType='none'
                else:
                    line=line.split('%')
                    if len(line)>1: comment = line[1], print('Warning: comments are not yet implemented. Removing comment.') #TODO: implement comments
                    line=line[0].strip()
                    if inputType=='name':
                        modelName=line
                    elif inputType=='states':
                        if re.search('d/dt\(',line,re.IGNORECASE):
                            match=re.search('d/dt\((\w+)\)\s*=\s*(.+)%*',line,re.IGNORECASE)
                            states[match.group(1)]=(match.group(2),)
                        elif re.search('(\w+)\(0\)',line,re.IGNORECASE):
                            match=re.search('(\w+)\(0\)\s*=\s*(.+)%*',line,re.IGNORECASE)
                            states[match.group(1)]+=(float(match.group(2)),)
                    elif inputType=='parameters':
                        match=re.search('(\w+)\s*=\s*(.+)%*',line,re.IGNORECASE)
                        params.append((match.group(1), float(match.group(2))))
                    elif inputType=='variables':
                        match=re.search('(\w+)\s*=\s*(.+)%*',line,re.IGNORECASE)
                        if line.startswith("y_"):
                            observables.append((match.group(1), match.group(2)))
                        else:
                            variables.append((match.group(1), match.group(2)))
                    elif inputType=='reactions':
                        match=re.search('(\w+)\s*=\s*(.+)%*',line,re.IGNORECASE)
                        reactions.append((match.group(1), match.group(2)))               
                    elif inputType=='events':
                        match=re.search("(\w+)\s*=\s*(\w+)\s*\((\w+)\s*,\s*([\w\.]+)\s*\)\s*,\s*(\w+)\s*,\s*([\w\.]+)\s*%*", line)
                        if match.group(2)=='eq':
                            sign='=='
                        elif match.group(2)=='lt':
                            sign='<'
                        elif match.group(2)=='gt':
                            sign='>'
                        elif match.group(2)=='le':
                            sign='<='
                        elif match.group(2)=='ge':
                            sign='>='
                        else:
                            print('Unknown event format: '+line)
                        events.append((match.group(1), match.group(3), sign, match.group(4), match.group(5), match.group(6)))


    states=[(k,*v) for k,v in states.items()]
    if events: 
        print("Warning, events are not fully supported. E.g. only one thing can be changed per event")

    if do_print:
        print('---ODEs---')
        for state, rhs,_ in states:
            print(f"    d/dt({state}) = {rhs}")
        print('---ICs---')
        for name, _, value in states:
            print(name,'=', value)
        print('---parameters---')
        for name, value in params:
            print(name, '=', value)
        print('---Variables---')
        for name, value in variables:
            print(name, '=', value)
        print('---reactions---')
        for name, value in reactions:
            print(name, '=', value)
        print('---observables---')
        for name, value in observables:
            print(name, '=', value)
        print('---events---')
        for name, x, cond, val,y,val2 in events:
            print(f"    {name}: {x} {cond} {val} -> {y} = {val2}")
    model=dict()
    model["name"]=modelName
    model["states"] = states
    model["parameters"] = params
    if variables: model["variables"]=variables
    if observables: model["observables"]=observables
    if reactions: model["reactions"]=reactions
    if events: model["events"]=events
    f.close()
    return model

def export_as_scipy(model, filename = None):
    """This function exports a model to the SciPy (odeint) format. 

    Examples:
        Exporting a model to SciPy with file name 'model.py'
            export_as_scipy('model.txt', 'model.py')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .py extension.
    
    """
    if filename:
        f=open(filename, "w")
    else:
        f=open(model["name"]+'.py','w')
    f.write("from math import exp as exp\n")
    f.write("from numba import jit\n\n")
    f.write("@jit\n")
    f.write("def {0}(state,t, param):\n".format(model["name"]))
    i=0
    f.write("\n#Defining starting values\n")
    state_names = [s[0] for s in model["states"]]
    for i,state in enumerate(state_names):
        f.write("  {0} = state[{1}]\n".format(state,i))

    f.write("\n#Defining parameter values\n")
    for i,param in enumerate([p[0] for p in model["parameters"]]):
        f.write("  {0} = param[{1}]\n".format(param,i))

    if any(key in model for key in ["variables", "observables", "reactions"]):
        f.write("\n#Defining variables\n")

    if "variables" in model:
        for var, val in model["variables"]:
            f.write(f"      {var} = {val}".replace('^','**')+"\n")

    if "observables" in model:
        print("Observables are not yet implemented, writing it as a variable.")
        for var, val in model["observables"]:
            f.write(f"      {var} = {val}".replace('^','**')+"\n")

    if "reactions" in model:
        f.write("\n#Defining reactions\n")
        for reaction, val in model["reactions"]:
            f.write(f"      {reaction} = {val}".replace('^','**')+"\n")

    if "events" in model:
        print("Events are not yet implemented")
     
    f.write("\n#Defining ODEs\n")
    for state, rhs,_ in model["states"]:
        f.write(f"      {state}_d = {rhs}".replace('^','**')+"\n")

    f.write("\n#Return ODE values\n")
    f.write("  return[")
    f.write(state_names[0]+"_d")

    for state in state_names[1:]:
        f.write(f"    , {state}_d")
    f.write("]")
    f.close()
    print(f'Converted {model["name"]} to SciPy model')

def export_as_yaml(model, filename=None):
    """This function exports a model to the yaml format. 

    Examples:
        Exporting a model to yaml with file name 'model.yaml'
            export_as_yaml('model.txt', 'model.yaml')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .yaml extension.
    
    """
    if filename:
        f=open(filename, "w")
    else:
        f=open(model["name"]+'.yml','w')
    f.write("time:  \n  variable: t\n\n") #not sure if t or time is the variable? 

    f.write("odes:  \n")
    for state, eq, ic in model["states"]:
        f.write(f"        - stateId: {state} \n")
        f.write(f"          rightHandSide: \"{eq}\" \n")
        f.write(f"          initialValue: {ic} \n\n")
        
    f.write("parameters:  \n")
    for param, val in model["parameters"]:
        f.write(f"        - parameterId: {param} \n")
        f.write(f"          nominalValue: {val} \n\n")
    
    if "variables" in model or "reactions" in model:
        f.write("assignments:  \n")
    if "variables" in model: 
        for var, val in model["variables"]:
            f.write(f"        - assignmentId: {var} \n")
            f.write(f"          formula: \"{val}\" \n\n")
    if "reactions" in model:
        for reaction, val in model["reactions"]:
            f.write(f"        - assignmentId: {reaction} \n")
            f.write(f"          formula: \"{val}\" \n\n")  

    if "observables" in model:
        print("Warning, noiseFormula is not implemented/readable from ODE file. Assuming it to be 0")
        f.write("observables:  \n")
        for obs, val in model["observables"]:
            f.write(f"        - observableId: {obs} \n")
            f.write(f"          observableFormula: \"{val}\" \n")
            f.write(f"          noiseFormula: 0.0 \n\n")
    f.close()
    print(f'Converted {model["name"]} to a YAML model')

def export_as_antimony(model, filename=None):
    """This function exports a model to the antimony format. 

    Examples:
        Exporting a model to antimony with file name 'model.txt'
            export_as_antimony(model, 'model_a.txt', 'yaml')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .txt extension.
    
    """
    if filename:
        f=open(filename, "w")
    else:
        f=open(model["name"]+'.txt','w')

    def reduce_equations(rhs, model):
        if "reactions" in model:
            reactions = dict(model["reactions"])
            for key, value in reactions.items():
                rhs = rhs.replace(key, value)
        if "variables" in model:
            variables = dict(model["variables"])
            for key, value in variables.items():
                rhs = rhs.replace(key, value)
        return rhs.replace(" ", "")
    
    outbound = []
    inbound = []
    for state,rhs,_ in model["states"]: 
        eq = reduce_equations(rhs, model)
        reactions =  re.split('(\+|\-)',eq)
        direction = "+"
        for reaction in reactions: 
            if reaction[0] == "-":
                direction = '-'
            elif reaction[0] == "+":
                direction = "+"
            else:
                if direction == "+":
                    inbound.append((state, reaction))
                elif direction =="-":
                    outbound.append((state, reaction))

    reactions=[]
    for state_out, eq_out in outbound:
        found = False
        for state_in, eq_in in inbound: 
            if eq_out == eq_in:
                reaction = (state_out, state_in, eq_out)
                found=True
                break
        if found:
            reactions.append(reaction)
        else:
            reactions.append((state_out, "", eq_out))
    eqs_out = [o[1] for o in outbound]
    for state_in, eq_in in inbound:
        if eq_in not in eqs_out:
            reactions.append(("", state_in, eq_in))
    for s1,s2,r in reactions:
        f.write(f"    {s1} -> {s2}; {r}\n")
    f.write("\n")
    for name, value in model["parameters"]:
        f.write(f"    {name} = {value};\n")
    f.write("\n")
    for name, _ ,value in model["states"]:
        f.write(f"    {name} = {value};\n")
    f.close()

def export_as_medigit(model, filename=None):
    """This function exports a model to the MeDigiT format. 

    For now assumes that all models are defined in seconds. 

    Examples:
        Exporting a model to medigit with file name 'model_mdt.txt'
            export_as_medigit(model, 'model_mdt.txt')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .txt extension.
    
    """
    if filename:
        f=open(filename, "w")
    else:
        f=open(model["name"]+'.txt','w')
    f.write("########## NAME\n")
    f.write(f"    {model['name']}\n")
    f.write("########## METADATA\n")
    print("Warning, time-unit is not imported. Assuming it is in seconds")
    f.write("timeunit = s\n")
    f.write("########## MACROS\n")
    if "macros" in model:
        for name, value in model["macros"]:
            f.write(f"    {name} = {value}\n")
    f.write("########## STATES\n")
    for state, rhs, _ in model["states"]:
        f.write(f"    d/dt({state}) = {rhs}\n")
    f.write("\n")
    for state, _, ic in model["states"]:
        f.write(f"    {state}(0) = {ic}\n")
    f.write("########## PARAMETERS\n")
    for name, value in model["parameters"]:
        f.write(f"    {name} = {value}\n")
    f.write("########## VARIABLES\n")
    if "variables" in model:
        for name, value in model["variables"]:
            f.write(f"    {name} = {value}\n")
    if "reactions" in model:
        for name, value in model["reactions"]:
            f.write(f"    {name} = {value}\n")
    f.write("########## FUNCTIONS\n")
    if "functions" in model:
        for name, value in model["functions"]:
            f.write(f"    {name} = {value}\n")
    f.write("########## EVENTS\n")
    if "events" in model:
        for name, value in model["events"]:
            f.write(f"    {name} = {value}")
    f.write("########## OUTPUTS\n")
    if "outputs" in model:
        for name, value in model["outputs"]:
            f.write(f"    {name} = {value}\n")
    f.write("########## INPUTS\n")
    if  "inputs" in model:
        for name, *value in model["inputs"]:
            f.write(f"    {name} = {value[0]} ")
            if len(value)>1: f.write(f"    @ {value[1]}")
            f.write("\n")
    f.write("########## FEATURES\n")
    if "observables" in model:
        for name, value in model["observables"]:
            f.write(f"    {name} = {value}\n")
    f.close()
    print(f'Converted {model["name"]} to a medigit model')

    #TODO: add optional time-step format
    #TODO: strip y_ for measurment variable? 

def export_as_latex(model, filename=None):
    """This function exports a model as a set of LaTeX equations. 
    Note that this function only calls the 'export_as_LaTeX' function. 
    """
    export_as_LaTeX(model, filename)

def export_as_LaTeX(model, filename=None):
    """This function exports a model as a set of LaTeX equations. 

    Examples:
        Exporting a model as LaTeX with file name 'model.tex'
            export_as_LaTeX(model, 'model.tex')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .tex extension.
    
    """
    if filename:
        f=open(filename, "w")
    else:
        f=open(model["name"]+'.tex','w')

    f.write("\\documentclass[12pt]{article}\n")
    f.write("\\usepackage[utf8]{inputenc}\n")
    f.write("\\usepackage{caption}\n")
    f.write("\\usepackage{amsmath}\n")
    f.write("\\usepackage{amssymb}\n")
    f.write("\\usepackage{float}\n")
    f.write("\\usepackage[capitalize]{cleveref}\n")

    f.write("\\begin{document}\n")

    f.write("% ODEs:\n")
    f.write("\section{System of ordinary different equations}\n")
    f.write("The ODEs in the model are given in \\cref{eq:odes}.\n")
    f.write("\\begin{equation}\n")
    f.write("    \\begin{alignedat}{4}\n")
    for state, rhs, _ in model["states"]:
        f.write(f"        &d/dt({state}) &&= &&{rhs} \\\\ \n".replace('_','\_').replace('*',' \cdot '))
    f.write("    \\label{eq:odes}\n")
    f.write("    \\end{alignedat}\n")
    f.write("\\end{equation}\n")

    f.write("\n% Initial values:\n")
    f.write("\section{Initial values}\n")
    f.write("The initial values of the model states are given in \\cref{eq:initialvalues} and \\cref{table:initialvalues}.\n")
    f.write("\\cref{eq:initialvalues} is given below.\n")

    f.write("\\begin{equation}\n")
    f.write("    \\begin{alignedat}{4}\n")
    for state, _, ic in model["states"]:
        ic="{:.4E}".format(ic)
        if "E+00" in ic:
            ic=ic.replace('E+00','')
        else:
           ic=ic.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"        &{state}(0) &&= &&{ic} \\\\ \n".replace('_','\_'))
    f.write("    \\label{eq:initialvalues}\n")
    f.write("    \\end{alignedat}\n")
    f.write("\\end{equation}\n")

    f.write("\\cref{table:initialvalues} is given below.\n")
    f.write("\\begin{table}[H]\n")
    f.write("\\centering\n")
    f.write("    \\begin{tabular}{ | l | l|} \\hline \n")
    f.write(f"        {'State':<12}   & Initial values               \\\\ \\hline \n".replace('_','\_'))
    for state, _, ic in model["states"]:
        ic="{:.4E}".format(ic)
        if "E+00" in ic:
            ic=ic.replace('E+00','')
        else:
           ic=ic.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"        ${state:<12}$ & ${ic:<25}$ \\\\ \\hline \n".replace('_','\_'))
    f.write("    \\end{tabular}\n")
    f.write("    \\caption{\\textbf{Initial values.}}\n")
    f.write("    \\label{table:initialvalues}\n")
    f.write("\\end{table}\n")

    f.write("The initial values are also given as a block of equations below in \\cref{eq:initialvalues_block}.\n")
    f.write("Note: You need to add line breaks yourself. Dont forget to add \\\\ before each new line. Also remember to remove the last trailing comma.\n")
    f.write("\\begin{equation}\n")
    f.write("    \\begin{alignedat}{2}\n")
    f.write(f"        &")

    for state, _, ic in model["states"]:
        ic="{:.4E}".format(ic)
        if "E+00" in ic:
            ic=ic.replace('E+00','')
        else:
           ic=ic.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"{state}(0)={ic}, \ ".replace('_','\_'))
    f.write("\n    \\label{eq:initialvalues_block}\n")
    f.write("    \\end{alignedat}\n")
    f.write("\\end{equation}\n")
    
    f.write("\n% Parameters:\n")
    f.write("\section{Parameter values}\n")
    f.write("The parameter values of the model states are given in \\cref{eq:parameters} and \\cref{table:parameters}.\n")
    f.write("\\cref{eq:parameters} is given below.\n")
    f.write("\\begin{equation}\n")
    f.write("    \\begin{alignedat}{4}\n")
    for name, value in model["parameters"]:
        value="{:.4E}".format(value)
        if "E+00" in value:
            value=value.replace('E+00','')
        else:
           value=value.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"        &{name} &&= &&{value} \\\\ \n".replace('_','\_'))
    f.write("    \\label{eq:parameters}\n")
    f.write("    \\end{alignedat}\n")
    f.write("\\end{equation}\n")
    
    f.write("\\cref{table:parameters} is given below.\n")
    f.write("\\begin{table}[H]\n")
    f.write("\\centering\n")
    f.write("    \\begin{tabular}{ | l | l|} \\hline \n")
    f.write(f"        {'Parameter':<12}   &  Value                   \\\\ \\hline \n")
    for name, value in model["parameters"]:
        value="{:.4E}".format(value)
        if "E+00" in value:
            value=value.replace('E+00','')
        else:
           value=value.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"        ${name:<12}$ & ${value:<25}$ \\\\ \\hline \n".replace('_','\_'))
    f.write("    \\end{tabular}\n")
    f.write("    \\caption{\\textbf{Parameter values.}}\n")
    f.write("    \\label{table:parameters}\n")
    f.write("\\end{table}\n")

    f.write("The parameter values are also given as a block of equations below in \\cref{eq:parameters_block}.\n")
    f.write("Note: You need to add line breaks yourself. Dont forget to add \\\\ before each new line. Also remember to remove the last trailing comma.\n")
    f.write("\\begin{equation}\n")
    f.write("    \\begin{alignedat}{2}\n")
    f.write(f"        &")

    for name, value in model["parameters"]:
        value="{:.4E}".format(value)
        if "E+00" in value:
            value=value.replace('E+00','')
        else:
           value=value.replace('E+0','\cdot 10^{').replace('E-0','\cdot 10^{-')+"}"
        f.write(f"{name}={ic}, \ ".replace('_','\_'))
    f.write("\n    \\label{eq:parameters_block}\n")
    f.write("    \\end{alignedat}\n")
    f.write("\\end{equation}\n")

    if "variables" in model:
        f.write("\n% Variable equations:\n")
        f.write("\section{Variables}\n")
        f.write("The variable equations are given in \\cref{eq:variables}.\n")

        f.write("\\begin{equation}\n")
        f.write("    \\begin{alignedat}{4}\n")
        for name, value in model["variables"]:
            f.write(f"        &{name} &&= &&{value} \\\\ \n".replace('_','\_').replace('*',' \cdot '))
        f.write("    \\label{eq:variables}\n")
        f.write("    \\end{alignedat}\n")
        f.write("\\end{equation}\n")

    if "reactions" in model:
        f.write("\n% Reaction equations:\n")
        f.write("\section{Reactions}\n")
        f.write("The reaction equations are given in \\cref{eq:reactions}.\n")

        f.write("\\begin{equation}\n")
        f.write("    \\begin{alignedat}{4}\n")
        for name, value in model["reactions"]:
            f.write(f"        &{name} &&= &&{value} \\\\ \n".replace('_','\_').replace('*',' \cdot '))
        f.write("    \\label{eq:reactions}\n")
        f.write("    \\end{alignedat}\n")
        f.write("\end{equation}\n")

    # if "functions" in model:
    # Not implemented

    # if "events" in model:       
    # Not implemented

    
    if "outputs" in model:
        f.write("\n% Output equations:\n")
        f.write("\section{Outputs}\n")
        f.write("The output equations are given in \\cref{eq:outputs}.\n")

        f.write("\\begin{equation}\n")
        f.write("    \\begin{alignedat}{4}\n")
        for name, value in model["outputs"]:
            f.write(f"        &{name} &&= &&{value} \\\\ \n".replace('_','\_'))
        f.write("    \\label{eq:outputs}\n")
        f.write("    \\end{alignedat}\n")
        f.write("\end{equation}\n")
    
    if  "inputs" in model:
        f.write("\n% Input equations:\n")
        f.write("\section{Inputs}\n")
        f.write("The input equations are given in \\cref{eq:inputs}.\n")

        f.write("\\begin{equation}\n")
        f.write("    \\begin{alignedat}{4}\n")
        for name, *value in model["inputs"]:
            f.write(f"        &{name} &&= &&{value[0]} ".replace('_','\_'))
            if len(value)>1: f.write(f"    @ {value[1]}")
            f.write(' \\\\ \n')
        f.write("    \\label{eq:inputs}\n")
        f.write("    \\end{alignedat}\n")
        f.write("\end{equation}\n")
    
    if "observables" in model:
        f.write("\n% Measurement equations:\n")
        f.write("\section{Measurement equations}\n")
        f.write("The measurement equations are given in \\cref{eq:observables}.\n")

        f.write("\\begin{equation}\n")
        f.write("    \\begin{alignedat}{4}\n")
        for name, value in model["observables"]:
            f.write(f"        &{name} &&= &&{value}\n".replace('_','\_'))
        f.write("    \\label{eq:observables}\n")
        f.write("    \\end{alignedat}\n")
        f.write("\end{equation}\n")

    f.write("\end{document}")
    f.close()
    print(f'Converted {model["name"]} to LaTeX equations')

def export_as_SBML(model, filename=None, route='yaml'): 
    """This function exports a model to SBML. 

    For SBML conversions, the model must first be converted via either yaml2sbml or Tellurium (te); see type argument.
    Assumes that all variables named y_* are observables. 

    Examples:
        Exporting a model to SBML with file name 'model.xml' using YAML 
            export_as_SBML(model, 'model.xml', 'yaml')
        Exporting a model to SBML with file name 'model.xml' using Tellurium 
            export_as_SBML(model, 'model.xml', 'te')

    Args: 
        model: 
            An model (typically imported) to be converted
        filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a .xml extension.
        route: 
            Desired route for the SBML conversion. Available options: {'yaml', 'te'}
    
    """
    import os
    if not filename:
        filename = model["name"]+'.xml'

    if route=='yaml':
        export_as_yaml(model, filename=model["name"]+"_temp.yml")
        import yaml2sbml
        yaml2sbml.validate_yaml(yaml_dir=model["name"]+"_temp.yml")
        yaml2sbml.yaml2sbml(model["name"]+"_temp.yml", filename)  
        if os.path.exists(model["name"]+"_temp.yml"): os.remove(model["name"]+"_temp.yml")

    elif route=='te':
        from warnings import warn
        warn("Note that there exists a compatibility issue when exporting to SBML using Tellurium and later importing into amici. Both things cannot be done in the same overall function call. Please convert the model first, and then in a new call compile the model in AMICI. See this page for reference: https://github.com/AMICI-dev/AMICI/issues/1483")
        export_as_antimony(model, filename=model["name"]+"_temp.txt")
        import tellurium as te
        r = te.loada(model["name"] + '_temp.txt',)
        r.exportToSBML(filename)
        if os.path.exists(model["name"]+"_temp.txt"): os.remove(model["name"]+"_temp.txt")
    else:
        error("Unknown option for how to create the SBML file. Acceptable options are 'yaml' for yaml+yaml2sbml or 'te' for antimony+Tellurium")
        
def odes2py(in_filename, out_filename=None, type = 'mdt', do_print=False):
    """This function can convert an textfile with ODEs in the IQM/SBtoolbox format to another type.

    For SBML conversions, the model must first be converted via either yaml2sbml or Tellurium (te); see type argument.
    Assumes that all variables named y_* are observables. 

    Examples:
        Converting a model in 'model.txt' to SBML with file name 'model.xml' using YAML: 
            odes2py('model.txt', 'model.xml', 'sbml-yaml')
        Converting a model in 'model.txt' to MeDigiT with default output file name:
            odes2py('model.txt', type = 'medigit')

    Args: 
        in_filename: 
            path to the model to be converted (including file extension)
        out_filename:
            Desired file name for the converted model (including file extension). If set to None then the name in the input file will be used, combined with a suitable extension.
        type: 
            Desired type for the converted model. Available options: ['scipy', 'medigit' | 'mdt', 'sbml-yaml', 'sbml-te', 'yaml', 'antimony', 'LaTeX'|'latex']
        do_print: 
            Set to True if you want the imported structure to be printed after importing. 
    Returns: 
        model (dict), containing the keys  ["name", "states", "parameters", "variables", "observables", "reactions", events"]
    
    """

    model = import_odes(in_filename, do_print)

    if type == "scipy":
        export_as_scipy(model, out_filename)
    elif type == "yaml":
        export_as_yaml(model, out_filename)
    elif type == "medigit":
        export_as_medigit(model, out_filename)
    elif type == "mdt":
        export_as_medigit(model, out_filename)
    elif type == "antimony" or type == "te" :
        export_as_antimony(model, out_filename)
    elif type == "sbml-yaml":
        export_as_SBML(model, out_filename)
    elif type == "sbml-te":
        export_as_SBML(model, out_filename, route = 'te')
    elif type == "latex" or type == "LaTeX" :
        export_as_LaTeX(model, out_filename)


    return model


if __name__ == '__main__':
    """ Command line usage for model conversions

    Use cases: 
        odes2py
        odes2py in_file.txt type
        odes2py in_file.txt type out_file.txt

    """
    from pathlib import Path
    
    argv = sys.argv[1:]
    if len(argv)==0:
        fileName=input('Enter model file name: ')
    else:
        fileName=argv[0]

    if len(argv)==1:
        type="mdt"
    else:
        type=argv[1]

    if len(argv)<3: 
        out_filename=f"    {Path(fileName).stem}-{type}{Path(fileName).suffix}"
    else:
        out_filename=argv[2]

    odes2py(fileName, out_filename, type)
