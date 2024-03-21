import argparse
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


RUN_COUNT = 1


def run(cmd: list[str]) -> str:
    result = subprocess.run(cmd, stderr=subprocess.PIPE)
    return result.stderr.decode("utf-8")


# return (secs, joules)
def parse_output(output: str) -> tuple[float, float]:
    lines = output.splitlines()

    secs = float(lines[0].rstrip(" sec"))
    joules = float(lines[2])

    return secs, joules


# return average (secs, joules)
def measure(alg: str, k: int = 0, r: int = 0) -> tuple[float, float]:
    secs_sum = 0.
    joules_sum = 0.
    run_count = RUN_COUNT

    for _ in range(run_count):
        output = run(["sudo", "turbostat", "--Summary", "--quiet", "--Joules", "--show", "Pkg_J", "./bin/run_enc_dec", alg, str(k), str(r)])

        secs, joules = parse_output(output)
        secs_sum += secs
        joules_sum += joules

    secs_avg = secs_sum / run_count
    joules_avg = joules_sum / run_count

    return secs_avg, joules_avg


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", type=Path, help="directory where to save plots")

    return parser.parse_args()


def calc_watts(secs_avg: float, joules_avg: float, no_secs_avg: float, no_joules_avg: float) -> float:
    return (joules_avg - no_joules_avg) / (secs_avg - no_secs_avg)


def plot_cdf(ax, data: list[float], title: str, xlabel: str):
    ax.set_title(title)
    ax.ecdf(data)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("CDF", fontsize=12)
    ax.grid(True)


if __name__ == "__main__":
    args = parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    no_secs_avg, no_joules_avg = measure("NO")

    rs_secs_arr = []
    rs_joules_arr = []
    rs_watts_arr = []

    rlc_secs_arr = []
    rlc_joules_arr = []
    rlc_watts_arr = []

    k = 2000
    r = 40
    for i in range(100):
        if i == 0 or (i + 1) % 10 == 0:
            print("Running test", i + 1)

        rs_secs_avg, rs_joules_avg = measure("RS", k, r)
        rs_watts = calc_watts(rs_secs_avg, rs_joules_avg, no_secs_avg, no_joules_avg)

        rlc_secs_avg, rlc_joules_avg = measure("RLC", k, r)
        rlc_watts = calc_watts(rlc_secs_avg, rlc_joules_avg, no_secs_avg, no_joules_avg)

        rs_secs_arr.append(rs_secs_avg - no_secs_avg)
        rs_joules_arr.append(rs_joules_avg - no_joules_avg)
        rs_watts_arr.append(rs_watts)

        rlc_secs_arr.append(rlc_secs_avg - no_secs_avg)
        rlc_joules_arr.append(rlc_joules_avg - no_joules_avg)
        rlc_watts_arr.append(rlc_watts)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    fig.suptitle("Сравнение энергоэффективности (k=2000, r=40)")

    plot_cdf(axs[0],
             list(map(lambda t: t[0] / t[1], zip(rs_joules_arr, rlc_joules_arr))),
             "Отношение затраченной энергии",
             r"$\frac{E_{RS}}{E_{RLC}}$")

    plot_cdf(axs[1],
             list(map(lambda t: t[0] / t[1], zip(rs_watts_arr, rlc_watts_arr))),
             "Отношение потребляемой мощности",
             r"$\frac{P_{RS}}{P_{RLC}}$")

    fig.tight_layout()
    fig.savefig(args.output / "energy_efficiency.png")
