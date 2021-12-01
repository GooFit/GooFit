#script to run fits on a single toy to assess error on 6D normalisation method
import numpy as np
import pandas as pd
import subprocess as sb
import os

executable_name = "WS_model_fit"
build_dir = "../build-cuda/k3pi_goofit_scripts/"

def fit_results_only_xy():
    seeds = [i for i in range(100)]
    for seed in seeds:
        sampleNum = "70000"
        x_start = "0.0045"
        y_start = "0.0062"
        input_file = "../build-cuda/k3pi_goofit_scripts/ws_generated_mc_samples/signal_only_70k_sample_x_0.0045_y_0.0062.txt"
        output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_seed_{seed}.txt"

        command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt","--genOff",str(seed)]
        print("executing command")
        output = None
        try:
            output = sb.check_output(command_arr)
        except sb.CalledProcessError as e:
            print(f"failed to execute command with seed: {seed}")
            print(str(e.cmd))
        print(output)
        
def fit_results_xy_floated_DCS():
    seeds = [i for i in range(100)]
    for seed in seeds:
        sampleNum = "70000"
        x_start = "0.0045"
        y_start = "0.0062"
        input_file = "../build-cuda/k3pi_goofit_scripts/ws_generated_mc_samples/signal_only_70k_sample_x_0.0045_y_0.0062.txt"
        output_file = f"../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_seed_{seed}_floated_DCS.txt"

        command_arr = [build_dir + executable_name,str(sampleNum),str(x_start),str(y_start),str(input_file),str(output_file),"--specInt","--genOff",str(seed),"--floatDCS"]
        print("executing command")
        output = None
        try:
            output = sb.check_output(command_arr)
        except sb.CalledProcessError as e:
            print(f"failed to execute command with seed: {seed}")
            print(str(e.cmd))
        print(output)
#fit_results_only_xy()
fit_results_xy_floated_DCS()