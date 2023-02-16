#!/usr/bin/env python

# This is a python script to test the output of the files.
# It requires plumbum (pip install plumbum or conda install -c conda-forge plumbum)


import re

try:
    from plumbum import TEE, cli, colors, local
except ImportError:
    print(
        "This file uses the plumbum library. Install with pip or conda (user directory or virtual environment OK)."
    )
    raise
import time


# Simple timer
class Timer:
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start


LOCAL_DIR = local.path(__file__).dirname


def test(filename, *args):
    command = local[LOCAL_DIR / filename / filename]
    command = command["--nosplash"]
    for arg in args:
        command = command[arg]
    colors.info.print("Running", command)
    with local.cwd(LOCAL_DIR / filename):
        with Timer() as t:
            code, stdout, stderr = command & TEE(retcode=None)
        if code == 0:
            colors.success.print(filename, "Successful")
        else:
            colors.fatal.print(filename, "Failed with status code:", code)
    return {
        "name": filename + " " + " ".join(map(str, args)),
        "code": code,
        "time": t.interval,
        "stdout": stdout,
        "stderr": stderr,
    }


def make_results(profile=False):
    results = [
        test("DP4"),
        test("SigGen"),
        test("2d_plot"),
        test("dalitz"),
        test("convolution"),
        test("addition"),
        test("exponential"),
        test("TDDP4", "-t", 100 if profile else 5),
        test("product"),
        test("simpleFit"),
        test("chisquare"),
    ]

    if (LOCAL_DIR / "zachFit/dataFiles/dstwidth_kpi_data.dat").exists():
        results.append(test("zachFit", 1, 0))
    if (LOCAL_DIR / "pipipi0DPFit/dataFiles/toyPipipi0/dalitz_toyMC_000.txt").exists():
        results.append(test("pipipi0DPFit", "toy", 0, 1))
        if profile:
            results.append(
                test(
                    "pipipi0DPFit",
                    "canonical",
                    "dataFiles/cocktail_pp_0.txt",
                    "--blindSeed=0",
                )
            )
    return results


MIN_TIME = re.compile(r"The minimization took: (.*?)$", re.MULTILINE)


class RunAll(cli.Application):
    profile = cli.Flag(["-p", "--profile"], help="Profile (longer) output")

    @cli.positional(int)
    def main(self, threads=None):
        env = local.env(OMP_NUM_THREADS=threads) if threads else local.env()
        with env:
            results = make_results(self.profile)
        failed = [result for result in results if result["code"] != 0]
        successes = [result for result in results if result["code"] == 0]

        if failed:
            colors.fatal.print("Failed:")
            for result in failed:
                colors.fatal.print(result["name"], result["code"])
        else:
            colors.success.print("All programs completed.")

        print()
        colors.info.print("{:20}:\tTotal time (s)\tFit times".format("Program"))
        for result in successes:
            fit = ", ".join(MIN_TIME.findall(result["stdout"]))
            print(
                (colors.success if result["code"] == 0 else colors.warn)
                | "{0[name]:20}:\t{0[time]}\t{1}".format(result, fit)
            )
        if threads:
            colors.info.print("OMP Threads:", threads)


if __name__ == "__main__":
    RunAll()
