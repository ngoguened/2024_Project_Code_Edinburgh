#INSTALL LIBRARIES AND PAS DEPENDENCY
using TreeParzen
using DifferentialEquations
using Plots
include("Models/P-Aminostyrene_Synthesis.jl")
using .PAminostyreneSynthesis
using Peaks
using Statistics
using SplitApplyCombine
using BenchmarkTools
using CSV
using Tables
using DataFrames
using FlexiMaps

# FLAGS FOR CONSTRAINTS: BYPASSES FMIN RESTRICTIONS
global is_oscillatory_flag::Bool
global proteomic_constraint_weight::Float64
global proteomic_rate_constraint_weight::Float64

# FLAG FOR REGULARIZATION
global complexity_regularization_weight::Float64

# ODE PROBLEM PARAMETERS
TIME_SPAN = (0.0, 1.73E5)
LEN_SAVE_TIMES = 200
SAVE_TIMES = range(0., 1.73E5, LEN_SAVE_TIMES)
INITIAL_VALUES = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# LOSS WEIGHT PARAMETERS
MIN_PRODUCTION = 1E1
MAX_PRODUCTION = 1E4
MIN_BURDEN = 1E-1
MAX_BURDEN = 1E3

# BAYESOPT PARAMETERS
FMIN_THRESHOLD = .9
STATIC_YIELD = 1.1232746768211363e-13 #median static yield from Latin Hypercube Sampling of the 5D k space on the no control arch.

# RUN STATE
global run_state = []

# BAYESOPT SPACE
space = Dict(
    :architecture => [HP.Choice(:prom1, [[0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]), HP.Choice(:prom2, [[1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]), HP.Choice(:prom3, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]])],
    :theta_paf_prom1 => HP.Uniform(:theta_paf_prom1, 5E-6, 1E-4), 
    :theta_paf_prom2 => HP.Uniform(:theta_paf_prom2, 5E-6, 1E-4), 
    :theta_paf_prom3 => HP.Uniform(:theta_paf_prom3, 5E-6, 1E-4), 
    :theta_paca_prom1 => HP.Uniform(:theta_paca_prom1, 5E-6, 1E-4),
    :theta_paca_prom2 => HP.Uniform(:theta_paca_prom2, 5E-6, 1E-4),
    :theta_paca_prom3 => HP.Uniform(:theta_paca_prom3, 5E-6, 1E-4),
    
    :k_paf_papA => HP.Uniform(:k_paf_papA, 5E-6, 1E-4),
    :k_paf_papB => HP.Uniform(:k_paf_papB, 5E-6, 1E-4),
    :k_paf_papC => HP.Uniform(:k_paf_papC, 5E-6, 1E-2),
    :k_paf_prom2 => HP.Uniform(:k_paf_prom2, 5E-6, 1E-4),
    :k_paf_prom3 => HP.Uniform(:k_paf_prom3, 5E-6, 1E-4),
    :k_paca_papA => HP.Uniform(:k_paca_papA, 5E-6, 1E-4), 
    :k_paca_papB => HP.Uniform(:k_paca_papB, 5E-6, 1E-4), 
    :k_paca_papC => HP.Uniform(:k_paca_papC, 5E-6, 1E-4), 
    :k_paca_prom2 => HP.Uniform(:k_paca_prom2, 5E-6, 1E-4),
    :k_paca_prom3 => HP.Uniform(:k_paca_prom3, 5E-6, 1E-4)
    )

# ODE HELPER FUNCTIONS
function unpack(space)
    arch = space[:architecture]
    thetas = [space[:theta_paf_prom1], space[:theta_paf_prom2], space[:theta_paf_prom3], 
              space[:theta_paca_prom1], space[:theta_paca_prom2], space[:theta_paca_prom3]]
    ks = [space[:k_paf_papA], space[:k_paf_papB], space[:k_paf_papC], space[:k_paf_prom2], space[:k_paf_prom3], 
          space[:k_paca_papA], space[:k_paca_papB], space[:k_paca_papC], space[:k_paca_prom2], space[:k_paca_prom3]]
    return arch, thetas, ks
end

function solve_ode(space)
    arch, thetas, ks = unpack(space)
    problem = ODEProblem(p_aminostyrene, INITIAL_VALUES, TIME_SPAN, saveat=SAVE_TIMES, [arch, thetas, ks])
    return solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
end

function calculate_relative_yield(solution)
    chorismate = solution[1,:]
    paca = solution[6,:]
    pas = solution[7,:]
    yield = (sum(paca) + sum(pas))/sum(chorismate)
    return yield/STATIC_YIELD
end

#CONSTRAINT HELPER FUNCTIONS
function is_oscillatory(solution, weight)
    for i in range(1, length(solution))  
        p = findmaxima(solution[i])
        idxs = p[1]
        dists = [p[2][i+1] - p[2][i] for i in range(1, length(p[2])-1)]
        if length(idxs) > 3 && std(dists) < 1
            return weight
        end
    end
    return 0
end

function proteomic_constraint(solution, constraint, weight)
    #laao papA papB papC p_efflux
    #do we care about unfolded proteins?
    for i in range(1, LEN_SAVE_TIMES)
        total_protein = solution.u[i][22] + solution.u[i][11] + solution.u[i][14] + solution.u[i][17] + solution.u[i][26]
        if total_protein > constraint
            return weight
        end
    end
    return 0
end

function proteomic_rate_constraint(solution, constraint, weight)
    total_protein_rate = 0 # mrna values as our protein production rates
    for i in range(2, LEN_SAVE_TIMES)
        total_protein_rate += max((solution.u[i][9]-solution.u[i-1][9]),0) + max((solution.u[i][12]-solution.u[i-1][12]),0) + max((solution.u[i][15]-solution.u[i-1][15]),0) +
        max((solution.u[i][20]-solution.u[i-1][20]),0) + max((solution.u[i][24]-solution.u[i-1][24]),0)
        if total_protein_rate > constraint
            return weight
        end
    end
    return 0
end

# REGULARIZATION HELPER FUNCTIONS
function calculate_loops_and_sensors(arch)
    number_of_loops = 0
    for a in arch
        if a[5] != 1
            number_of_loops += 1
        end
    end
    sensors = [0,0,0,0]
    for a in arch
        for i in range(1,length(a)-1)
            if a[i] == 1
                sensors[i] = 1
            end
        end
    end
    number_of_sensors = sum(sensors)
    return number_of_loops, number_of_sensors
end

function complexity(arch, weight)
    loops, sensors = calculate_loops_and_sensors(arch)
    return weight * (loops + 2*sensors)
end


# BAYESOPT HELPER FUNCTIONS
function push_state(space, obj, yield)
    arch, thetas, ks = unpack(space)
    push!(run_state, [arch, thetas, ks, obj, yield])
end

function write_state(file, state, iterations)
    df = DataFrame("Iteration" => [i for i in range(1, iterations)], "Architecture" => [state[i][1] for i in range(1, iterations)], 
                            "Thetas" => [state[i][2] for i in range(1, iterations)], "Ks" => [state[i][3] for i in range(1, iterations)], 
                            "Objective" => [state[i][4] for i in range(1, iterations)], "Relative Yield" => [state[i][5] for i in range(1, iterations)],)

        CSV.write(file, df, append=isfile(file))
    
        
end

function loss(production, burden, space)
    alpha = 0.6
    j1 = (production - MIN_PRODUCTION)/(MAX_PRODUCTION - MIN_PRODUCTION)
    j2 = (burden - MIN_BURDEN)/(MAX_BURDEN - MIN_BURDEN)
    arch, _, _ = unpack(space)
    if complexity_regularization_weight != 0
        return alpha*j1 + (1-alpha)*j2 + complexity(arch, complexity_regularization_weight)
    end
    return alpha*j1 + (1-alpha)*j2
end

function objective(space)
    solution = solve_ode(space)
    production, burden = solution[end][7], solution[end][end]
    obj = loss(1/production, burden, space)
    if !iszero(proteomic_constraint_weight)
        obj += proteomic_constraint(solution, proteomic_constraint_weight, 1E8)
    end
    if !iszero(proteomic_rate_constraint_weight)
        obj += proteomic_rate_constraint(solution, proteomic_rate_constraint_weight, 1E8)
    end
    if is_oscillatory_flag
        obj += is_oscillatory(solution, 1E8)
    end
    yield = calculate_relative_yield(solution)
    push_state(space, obj, yield)
    return Float64(obj)
end

# BAYESOPT MAIN LOOP
function bayesopt_loop(iterations, space; proteomic_constraint=0, proteomic_rate_constraint=0, 
                       oscillatory_constraint::Bool=true, complexity_regularization=0, write_to_file::String="")
    global is_oscillatory_flag = oscillatory_constraint
    global proteomic_constraint_weight = proteomic_constraint
    global proteomic_rate_constraint_weight = proteomic_rate_constraint
    global complexity_regularization_weight = complexity_regularization

    best = fmin(objective, space, iterations, threshold=FMIN_THRESHOLD)

    if write_to_file != ""
        write_state(write_to_file, run_state, iterations)
    end

    global run_state = []

    return best
end


#TESTING
function main()
        for i in range(1,5)
            bayesopt_loop(1000, space, write_to_file="uhoh.csv")
        end
end
main()