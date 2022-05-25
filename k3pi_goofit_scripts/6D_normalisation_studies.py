#script to run fits on a single toy to assess error on 6D normalisation method
import numpy as np
import pandas as pd
import subprocess as sb
import os
import argparse
from joblib import Parallel, delayed

executable_name = "WS_model_fit"
build_dir = "../build-cuda/k3pi_goofit_scripts/"
gen_executable_name = "WS_model_gen"
numDisjointSamples = 300
x = "0.0045" 
y = "0.0062"
nGenEvts = 70000 # will need to be changed for the full Run 2 WS sample

parser = argparse.ArgumentParser()
parser.parse_args()

def gen_disjoint_samples(gpu_id = 0,seeds=range(int(numDisjointSamples))):
    #generate samples with x = 0.0045 and y = 0.0062 but with different seeds to generate different events for each toy
    #We can then fit each toy (using the same normalisation, by default) to obtain a distribution of fitted x and y values over the ensemble of toys
    for i in seeds:
        batchNum = 1000000
        nGenEvts = 70000 # will need to be changed for the full Run 2 WS sample
        seed = i
        output_file = build_dir+f"ws_generated_mc_samples/signal_only_{nGenEvts}_sample_num_{i}_x_{x}_y_{y}.txt"
        command_arr = [build_dir+ gen_executable_name, str(batchNum), str(nGenEvts), str(seed), str(output_file), x, y, "--gpu-dev", str(gpu_id)]
        output = None
        try:
            output = sb.check_output(command_arr)
        except sb.CalledProcessError as e:
            print(f"failed to execute command with seed: {seed}")
            print(str(e.cmd))
        print(output)

def fit_disjoint_samples():
    pass

def fit_results_only_xy(gpu_id = 0,seeds=range(int(numDisjointSamples)),special_integral=True):
    #seeds = [i for i in range(100)]
    for seed in seeds:
        sampleNum = "70000"
        x_start = "0.0045"
        y_start = "0.0062"
        
        #input_file = "../build-cuda/k3pi_goofit_scripts/ws_generated_mc_samples/signal_only_70k_sample_x_0.0045_y_0.0062.txt"
        input_file = build_dir+f"ws_generated_mc_samples/signal_only_{nGenEvts}_sample_num_{seed}_x_{x}_y_{y}.txt"
        output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_{seed}_same_norm_seed.txt"
        #uncomment below if you want to change the normalisation seed when fitting
        #command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt","--genOff",str(seed), "--gpu-dev", str(gpu_id)]
        command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt", "--gpu-dev", str(gpu_id)]
        if not special_integral:
            output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_{seed}_same_norm_seed_default_integral.txt"
            command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file), "--gpu-dev", str(gpu_id)]
        print("executing command")
        output = None
        try:
            output = sb.check_output(command_arr)
        except sb.CalledProcessError as e:
            print(f"failed to execute command with seed: {seed}")
            print(str(e.cmd))
        print(output)
        
def fit_results_xy_floated_DCS(gpu_id = 0,seeds=range(int(numDisjointSamples)),special_integral=True):
    #seeds = [i for i in range(100)]
    for seed in seeds:
        sampleNum = "70000"
        x_start = "0.0045"
        y_start = "0.0062"
        input_file = build_dir+f"ws_generated_mc_samples/signal_only_{nGenEvts}_sample_num_{seed}_x_{x}_y_{y}.txt"
        output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_{seed}_same_norm_seed_floated_dcs.txt"

        #command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt","--genOff",str(seed),"--floatDCS", "--gpu-dev", str(gpu_id)]
        command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt","--floatDCS" , "--gpu-dev", str(gpu_id)]
        if not special_integral:
            output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_{seed}_same_norm_seed_floated_dcs_default_integral.txt"
            command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--floatDCS" , "--gpu-dev", str(gpu_id)]
        print("executing command")
        output = None
        try:
            output = sb.check_output(command_arr)
        except sb.CalledProcessError as e:
            print(f"failed to execute command with seed: {seed}")
            print(str(e.cmd))
        print(output)
#Parallel(n_jobs=2)(delayed(gen_disjoint_samples)(id,seeds) for (id,seeds) in [(0,range(int(numDisjointSamples/2))),(1,range(int(numDisjointSamples/2),numDisjointSamples))])
#fit toys using 6D numerical normalisation integral
#Parallel(n_jobs=2)(delayed(fit_results_only_xy)(id,seeds) for (id,seeds) in [(0,range(int(numDisjointSamples/2))),(1,range(int(numDisjointSamples/2),numDisjointSamples))])
#Parallel(n_jobs=2)(delayed(fit_results_xy_floated_DCS)(id,seeds) for (id,seeds) in [(0,range(int(numDisjointSamples/2))),(1,range(int(numDisjointSamples/2),numDisjointSamples))])
#fit toys using the default normalisation integral
Parallel(n_jobs=2)(delayed(fit_results_only_xy)(id,seeds,special_int) for (id,seeds,special_int) in [(0,range(int(numDisjointSamples/2)),False),(1,range(int(numDisjointSamples/2),numDisjointSamples),False)])
Parallel(n_jobs=2)(delayed(fit_results_xy_floated_DCS)(id,seeds,special_int) for (id,seeds,special_int) in [(0,range(int(numDisjointSamples/2)),False),(1,range(int(numDisjointSamples/2),numDisjointSamples),False)])
#fit_results_only_xy()
#fit_results_xy_floated_DCS()