module PAminostyreneSynthesis
export p_aminostyrene
using DifferentialEquations
using Plots
using BenchmarkTools

function activation(x, k, theta)
    return (k*(x/theta)^2)/(1+(x/theta)^2)
end
function repression(x, k, theta)
    return k/(1+(x/theta)^2)
end
function p_aminostyrene(du, u, p, t)
    # chorismate_production_rate = 1100 #1/s
    translation_initiation_rate = 2E-1 #1/s
    ribosome_elongation_rate = 20 #1/s
    chorismate_production_rate = 150 #1/s
    deaminase_production_rate = 1E1 #1/s
    pAF_loss = 1.4E-5 #1/s
    mrna_degradation_rate = 3E-3 #1/s
    protein_degradation_rate = 2E-4 #1/s
    protein_folding_rate = 2E0 #1/s
    dilution_rate = 5.79E-4 #1/s
    dna_duplication_rate = 5.78E-4 #1/s
    avogadro_number = 6.02214E23 #NA
    cell_volume = 2.5E-17 #L
    toxicity_constant = 5E-5 #M
    metabolite_induced_toxicity = 5E-4 #NA
    protein_induced_toxicity = 50 #NA
    enzyme_induced_toxicity = 50 #NA
    pap_operon_mrna_length = 3400 #nt
    efflux_pump_mrna_length = 2900 #nt
    laao_mrna_length = 1600 #nt
    deaminase_kcat = 5E0 #M/s
    deaminase_km = 1E-6 #M
    papA_kcat = 0.2975 #M/s
    papA_km = 0.056 #M
    papB_kcat = 39 #M/s
    papB_km = 0.38 #M
    papC_kcat = 20.44 #M/s
    papC_km = 0.555 #M
    laao_kcat = 1.29 #M/s
    laao_km = 10.82 #M
    # efflux_pump_rate = 275
    efflux_pump_rate = 5E1 #M

    chorismate, pA1, pA2, pA3, pAF, pACA, pAS,
    pr_1, mrna_papA, papA_uf, papA, 
    mrna_papB, papB_uf, papB,
    mrna_papC, papC_uf, papC, deaminase,
    pr_2, mrna_laao, laao_uf, laao,
    pr_3, mrna_p_efflux, p_efflux_uf, p_efflux, j1, j2 = u
    architecture, thetas, ks = p
    theta_1_pAF, theta_2_pAF, theta_3_pAF, theta_1_pACA, theta_2_pACA, theta_3_pACA = thetas
    k_1_pAF, k_2_pAF, k_3_pAF, k_4_pAF, k_5_pAF, k_1_pACA, k_2_pACA, k_3_pACA, k_4_pACA, k_5_pACA = ks

    toxicity_factor = 1 / 1 + (pACA/toxicity_constant/metabolite_induced_toxicity) + (p_efflux/toxicity_constant/protein_induced_toxicity) + (laao/toxicity_constant/enzyme_induced_toxicity)
    # toxicity_factor = toxicity_constant / (toxicity_constant + (pACA / metabolite_induced_toxicity) + (p_efflux / protein_induced_toxicity) + (laao / enzyme_induced_toxicity))
    external_efflux = p_efflux * ((pACA / avogadro_number) / cell_volume) * efflux_pump_rate * toxicity_factor

    # papA_catalyzed_biosynthesis = kinetic_reaction(papA, chorismate, papA_kcat, papA_km) * toxicity_factor
    # papB_catalyzed_biosynthesis = kinetic_reaction(papB, pA1, papB_kcat, papB_km) * toxicity_factor
    # papC_catalyzed_biosynthesis = kinetic_reaction(papC, pA2, papC_kcat, papC_km) * toxicity_factor
    papA_catalyzed_biosynthesis = papA_kcat * chorismate / (papA_km + chorismate)
    papB_catalyzed_biosynthesis = papB_kcat * papA / (papB_km + papA)
    papC_catalyzed_biosynthesis = papC_kcat * papB / (papC_km + papC)
    
    deaminase_catalyzed_biosynthesis = deaminase_kcat * ((pA3 / avogadro_number) / cell_volume) / (deaminase_km + ((pA3 / avogadro_number) / cell_volume)) * toxicity_factor
    laao_catalyzed_biosynthesis = laao_kcat * ((pAF/ avogadro_number) / cell_volume) / (laao_km + ((pAF / avogadro_number) / cell_volume)) * toxicity_factor
    chorismate_biosynthesis = chorismate_production_rate * toxicity_factor
    deaminase_biosynthesis = deaminase_production_rate * toxicity_factor
    
    papA_mrna_txn = sum(architecture[1].*[activation(pAF, k_1_pAF, theta_1_pAF), activation(pACA, k_1_pACA, theta_1_pACA), repression(pAF, k_1_pAF, theta_1_pAF), repression(pACA, k_1_pACA, theta_1_pACA), k_1_pAF])
    papB_mrna_txn = sum(architecture[1].*[activation(pAF, k_2_pAF, theta_1_pAF), activation(pACA, k_2_pACA, theta_1_pACA), repression(pAF, k_2_pAF, theta_1_pAF), repression(pACA, k_2_pACA, theta_1_pACA), k_2_pAF])
    papC_mrna_txn = sum(architecture[1].*[activation(pAF, k_3_pAF, theta_1_pAF), activation(pACA, k_3_pACA, theta_1_pACA), repression(pAF, k_3_pAF, theta_1_pAF), repression(pACA, k_3_pACA, theta_1_pACA), k_3_pAF])
    laao_mrna_txn = sum(architecture[2].*[activation(pAF, k_4_pAF, theta_2_pAF), activation(pACA, k_4_pACA, theta_2_pACA), repression(pAF, k_4_pAF, theta_2_pAF), repression(pACA, k_4_pACA, theta_2_pACA), k_4_pAF])
    eff_mrna_txn = sum(architecture[3].*[activation(pAF, k_5_pAF, theta_3_pAF), activation(pACA, k_5_pACA, theta_3_pACA), repression(pAF, k_5_pAF, theta_3_pAF), repression(pACA, k_5_pACA, theta_3_pACA), k_5_pAF])
    
    du[1] = chorismate_biosynthesis - papA_catalyzed_biosynthesis - dilution_rate * chorismate
    du[2] = papA_catalyzed_biosynthesis - papB_catalyzed_biosynthesis - dilution_rate * pA1
    du[3] = papB_catalyzed_biosynthesis - papC_catalyzed_biosynthesis - dilution_rate * pA2
    du[4] = papC_catalyzed_biosynthesis - deaminase_catalyzed_biosynthesis - dilution_rate * pA3
    du[5] = deaminase_catalyzed_biosynthesis - laao_catalyzed_biosynthesis - dilution_rate * pAF - pAF_loss * pAF
    du[6] = laao_catalyzed_biosynthesis - external_efflux - dilution_rate * pACA
    du[7] = external_efflux
    du[8] = pr_1 * dna_duplication_rate - dilution_rate * pr_1
    du[9] = papA_mrna_txn - dilution_rate * mrna_papA - toxicity_factor * mrna_degradation_rate * mrna_papA
    du[10] = mrna_papA / (translation_initiation_rate + (pap_operon_mrna_length / ribosome_elongation_rate)) * toxicity_factor - papA_uf * protein_folding_rate * toxicity_factor - papA_uf * dilution_rate - papA_uf * protein_degradation_rate  * toxicity_factor #papA_uf 
    du[11] = papA_uf * protein_folding_rate * toxicity_factor - papA * dilution_rate - papA * protein_degradation_rate * toxicity_factor  #papA 
    du[12] = papB_mrna_txn - dilution_rate * mrna_papB - toxicity_factor * mrna_degradation_rate * mrna_papB
    du[13] = mrna_papB / (translation_initiation_rate + (pap_operon_mrna_length / ribosome_elongation_rate)) * toxicity_factor - papB_uf * protein_folding_rate * toxicity_factor - papB_uf * dilution_rate - papB_uf * protein_degradation_rate  * toxicity_factor #papB_uf 
    du[14] = papB_uf * protein_folding_rate * toxicity_factor - papB * dilution_rate - papB * protein_degradation_rate * toxicity_factor  #papB
    du[15] = papC_mrna_txn - dilution_rate * mrna_papC - toxicity_factor * mrna_degradation_rate * mrna_papC #mrna_papC
    du[16] = mrna_papC / (translation_initiation_rate + (pap_operon_mrna_length / ribosome_elongation_rate)) * toxicity_factor - papC_uf * protein_folding_rate * toxicity_factor - papC_uf * dilution_rate - papC_uf * protein_degradation_rate  * toxicity_factor #papC_uf 
    du[17] = papC_uf * protein_folding_rate * toxicity_factor - papC * dilution_rate - papC * protein_degradation_rate * toxicity_factor  #papC
    du[18] = deaminase_biosynthesis - deaminase * dilution_rate #deaminase
    du[19] = pr_2 * dna_duplication_rate - dilution_rate * pr_2 #pr_2
    du[20] = laao_mrna_txn - dilution_rate * mrna_laao - toxicity_factor * mrna_degradation_rate * mrna_laao #mrna_laao
    du[21] = mrna_laao / (translation_initiation_rate + (laao_mrna_length / ribosome_elongation_rate)) * toxicity_factor - protein_folding_rate * toxicity_factor * laao_uf - dilution_rate * laao_uf - protein_degradation_rate * toxicity_factor * laao_uf #laao_uf
    du[22] = protein_folding_rate * toxicity_factor * laao_uf - dilution_rate * laao -  protein_degradation_rate * toxicity_factor * laao #laao
    du[23] = pr_3 * dna_duplication_rate - dilution_rate * pr_3
    du[24] = eff_mrna_txn - mrna_p_efflux * dilution_rate - mrna_p_efflux * mrna_degradation_rate * toxicity_factor
    du[25] = mrna_p_efflux / (translation_initiation_rate + (efflux_pump_mrna_length / ribosome_elongation_rate)) * toxicity_factor - protein_folding_rate * p_efflux_uf * toxicity_factor - dilution_rate * p_efflux_uf - protein_degradation_rate * toxicity_factor * p_efflux_uf
    du[26] = protein_folding_rate * toxicity_factor * p_efflux_uf - dilution_rate * p_efflux - protein_degradation_rate * toxicity_factor * p_efflux
    
    du[27] = (chorismate_biosynthesis - external_efflux)^2 #J1
    du[28] = papA_mrna_txn + papB_mrna_txn + papC_mrna_txn + laao_mrna_txn + eff_mrna_txn #J2
end

# time_span = (0.0, 50000)
# save_times = range(0., 50000, 200)
# initial_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# architecture = [[0, 0, 0, 0, 1], [0, 0, 0, 0, 1], [0, 0, 0, 0, 1]]
# thetas = [1E-4,1E-4,1E-4,1E-4,1E-4,1E-4]
# ks = [1E-5,1E-4,1E-4,1E-5,1E-7,1E-4,1E-4,1E-4,1E-4,1E-4]
# problem = ODEProblem(p_aminostyrene, initial_values, time_span, saveat=save_times, [architecture, thetas, ks])
# for i in range(1, 25)
#     solution1 = @btime solve(problem, FBDF(), reltol=1e-3, abstol=1e-6)
# end
# solution2 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
# solution3 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
# solution4 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
# solution5 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
# solution6 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6    )
# solution7 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)
# solution8 = @btime solve(problem, Rosenbrock23(), reltol=1e-3, abstol=1e-6)

end
