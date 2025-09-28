import os
import string
from io import StringIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

CONC_ABSORPTION_WAVELENGTH = 450  # wavelength of the absorption coefficient entered by the user

DIRECTORY_PATH = './data'

GRAPH_Y_MARGIN = 0.1  # Margin percentage of the y-Axis

FIT_ENABLE_BY_DEFAULT = True
FIT_GRAPH_DATAPOINT_N = 500
FIT_PEAK_PERCENTAGE = 0.95

NORMALIZE_ENABLE_BY_DEFAULT = True

WAVELENGTH_MIN = 350
WAVELENGTH_MAX = 800

INTERNAL_X_ID = 'Wavelength (nm)'
INTERNAL_Y_ID = 'Abs (a.u.)'


def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def get_sorted_csv_files(directory):
    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter for .csv files
    csv_files = [f for f in all_files if f.lower().endswith('.csv') and os.path.isfile(os.path.join(directory, f))]

    # Sort alphabetically
    csv_files.sort()

    return csv_files


def print_file_list(csv_files):
    print(f"There are {len(csv_files)} csv files available.")
    for n in range(len(csv_files)):
        print(f"({n}): {csv_files[n]}")


def read_csv(file_path):
    # Read lines until the first empty line
    with open(DIRECTORY_PATH + '/' + file_path, 'r') as f:
        lines = []
        for line in f:
            if not line.strip():  # empty line -> stop; below that, there is metadata which breaks the csv
                break
            lines.append(line)
    return pd.read_csv(StringIO(''.join(lines)))  # Create a file-like object which pandas can handle


def process_sample(sam):
    # Extract wavelength and absorbance columns
    wavelength = sam["df"].iloc[1:, sam["row_id"]].astype(float)
    absorbance = sam["df"].iloc[1:, sam["row_id"] + 1].astype(float)

    # Create a new DataFrame
    data = pd.DataFrame({
        INTERNAL_X_ID: wavelength,
        INTERNAL_Y_ID: absorbance
    })

    if sam["ext_coeff"]:  # Calculate concentration
        print(f"-- Calculating concentration for {sam["label"]}) --")
        conc_data_diff = (data[
                              INTERNAL_X_ID] - CONC_ABSORPTION_WAVELENGTH).abs()  # This seems rather inefficient but is good enough
        conc_idx = conc_data_diff.idxmin()
        print(f"Closest datapoint: {data[INTERNAL_X_ID][conc_idx]} nm | {data[INTERNAL_Y_ID][conc_idx]}")
        print(f"With extinction coefficient {sam["ext_coeff"]}")
        conc = data[INTERNAL_Y_ID][conc_idx] / sam["ext_coeff"]
        print(f"Concentration: {conc} M")

    # Normalize absorbance
    if sam["normalize"]:
        mask_normalization = (data[INTERNAL_X_ID] >= WAVELENGTH_MIN) & (data[INTERNAL_X_ID] <= WAVELENGTH_MAX)
        max_abs = data.loc[mask_normalization, INTERNAL_Y_ID].max()
        data[INTERNAL_Y_ID] = data[INTERNAL_Y_ID] / max_abs

    # Filter data between 350 and 800 nm for plotting
    mask_plot = (data[INTERNAL_X_ID] >= WAVELENGTH_MIN) & (data[INTERNAL_X_ID] <= WAVELENGTH_MAX)
    return data.loc[mask_plot]


# CSV files can contain one or multiple samples
# Baselines will be ignored
# Lets the user choose a sample (if there are multiple)
def select_sample(df):
    samples = []
    for k, column in enumerate(df.columns):
        if column.startswith("Unnamed") or column.startswith("Baseline"):
            continue
        samples.append((column, k))  # column: sample name

    if len(samples) == 1:
        print(f"There is only one sample in the file. Selecting \"{samples[0][0]}\"")
        return samples[0]
    else:
        while True:
            # list samples in the csv file
            print(f"There are {len(samples)} samples in the file.")
            for j, sample in enumerate(samples):
                print(f"({j}): {sample[0]}")
            # select one
            inp = input("Select a sample: ")
            try:
                input_number = int(inp)
                if input_number < len(samples):
                    return samples[input_number]
                else:
                    print(f"The entered number {input_number} is out of range.")
            except ValueError:
                print("Invalid input.")


# Fit a gaussian peak around the maximum peak in the data
# The fit will be applied between FIT_PEAK_PERCENTAGE * Y_MAX to the left and right of the maximum
def fit_peak(fit_data):
    # Find index of maximum absorption
    peak_idx = fit_data[INTERNAL_Y_ID].idxmax()
    peak_x = fit_data[INTERNAL_X_ID][peak_idx]
    peak_y = fit_data[INTERNAL_Y_ID][peak_idx]
    # Split left and right of peak
    left_df = fit_data.loc[:peak_idx].iloc[::-1]  # reversed to scan outward from peak
    right_df = fit_data.loc[peak_idx:]

    # Find first x where Abs < PEAK_PERCENTAGE on each side
    x_right = left_df[left_df[INTERNAL_Y_ID] < FIT_PEAK_PERCENTAGE * peak_y][INTERNAL_X_ID].iloc[0]
    x_left = right_df[right_df[INTERNAL_Y_ID] < FIT_PEAK_PERCENTAGE * peak_y][INTERNAL_X_ID].iloc[0]
    print(f"Fitting between {x_left} and {x_right}.")
    mask_fit = (fit_data[INTERNAL_X_ID] <= x_right) & (fit_data[INTERNAL_X_ID] >= x_left)
    popt_gauss, pcov_gauss = scipy.optimize.curve_fit(gaussian, fit_data.loc[mask_fit, INTERNAL_X_ID].values,
                                                      fit_data.loc[mask_fit, INTERNAL_Y_ID].values,
                                                      p0=[peak_y, peak_x, 10])
    print(popt_gauss)
    lin_space_x = np.linspace(x_left, x_right, FIT_GRAPH_DATAPOINT_N)
    return lin_space_x, gaussian(lin_space_x, *popt_gauss), popt_gauss[1]


sorted_csvs = get_sorted_csv_files(DIRECTORY_PATH)  # index csv files
selected_samples = []

user_input = 1  # while loop loops until input is an empty line

while user_input:
    print_file_list(sorted_csvs)  # numerical index

    if len(selected_samples) != 0:  # small letters as index
        print('Selected files: ')
        for i, sel in enumerate(selected_samples):
            print(f"({string.ascii_lowercase[i]}) {sel["path"]}: {sel["label"]}" + (
                " [N]" if sel["normalize"] else "") + (" [F]" if sel["fit"] else "") + (
                      " [C]" if sel["ext_coeff"] else ""))

    print("Enter a number to select a file")
    print("Enter [letter] [command] to modify options for a file (see documentation)")
    print("Press enter to draw the graph")

    user_input = input("Command: ")
    if user_input:  # false if empty input: Then draw the graph
        match user_input:
            case user_input if user_input.isdigit():  # open a csv file
                value = int(user_input)
                if value < 0 or value >= len(sorted_csvs):
                    print("Invalid input.")
                else:
                    selected_file_path = sorted_csvs[value]
                    csv_file = read_csv(selected_file_path)  # read file
                    # CSV files can contain one or multiple samples
                    name, row_id = select_sample(csv_file)
                    selected_samples.append(
                        {"df": csv_file, "path": selected_file_path, "row_id": row_id, "label": name,
                         "normalize": NORMALIZE_ENABLE_BY_DEFAULT, "fit": FIT_ENABLE_BY_DEFAULT, "ext_coeff": False})

            case user_input if len(user_input) > 2 and user_input[0].isalpha() and user_input[1] == ' ':
                target = user_input[0]
                cmd = user_input[2:]
                target_i = ord(target.lower()) - ord(
                    'a')  # this might break if there are unreasonably many samples selected
                if target_i >= len(selected_samples):
                    print(f"Invalid: input ({target}) is not set.")
                    continue
                match cmd:
                    case "f e" | "fit enable":
                        selected_samples[target_i]["fit"] = True
                        print(f"Enabled fit for {selected_samples[target_i]['path']}.")

                    case "f d" | "fit d":
                        selected_samples[target_i]["fit"] = False
                        print(f"Disabled fit for {selected_samples[target_i]['path']}.")

                    case "n e" | "normalize enable":
                        selected_samples[target_i]["normalize"] = True
                        print(f"Enabled normalization for {selected_samples[target_i]['path']}.")

                    case "n d" | "normalize d":
                        selected_samples[target_i]["normalize"] = False
                        print(f"Disabled normalization for {selected_samples[target_i]['path']}.")

                    case "rm" | "remove":
                        print(f"Removing {selected_samples[target_i]['path']} from selection.")
                        selected_samples.pop(target_i)
                        print("Removed.")

                    case label if label.startswith("label "):
                        label = label[6:]
                        selected_samples[target_i]["label"] = label
                        print(f"Set label for {selected_samples[target_i]['path']} as {label}.")

                    case input_conc if input_conc.startswith("calc-concentration"):
                        if len(input_conc) < 19:
                            print("You need to specify the extinction coefficient")
                        else:
                            try:
                                selected_samples[target_i]["ext_coeff"] = float(
                                    input_conc[18:])  # i.e.: 1.38e8 for 13 nm AuNPs
                                print(
                                    f"Set extinction coefficient for {selected_samples[target_i]['path']} as {selected_samples[target_i]['ext_coeff']}.")
                            except ValueError:
                                print("Invalid input for extinction coefficient.")

                    case input_conc if input_conc.startswith("cc"):
                        if len(input_conc) < 4:
                            print("You need to specify the extinction coefficient")
                        else:
                            try:
                                selected_samples[target_i]["ext_coeff"] = float(
                                    input_conc[3:])  # i.e.: 1.38e8 for 13 nm AuNPs
                                print(
                                    f"Set extinction coefficient for {selected_samples[target_i]['path']} as {selected_samples[target_i]['ext_coeff']}.")
                            except ValueError:
                                print("Invalid input for extinction coefficient.")

                    case _:
                        print("Invalid command.")

            case "":
                print("Generating graph")

            case _:
                print("Invalid input")

plot_data = [process_sample(sample) for sample in selected_samples]  # only plot in specified range, normalize, fit

plt.figure(figsize=(10, 5))

plot_y_min = 0  # absolute maximum
plot_y_max = 0  # absolute minimum

for i, plot in enumerate(plot_data):
    print(f"-- Fitting {selected_samples[i]["label"]} --")
    plt.plot(plot[INTERNAL_X_ID], plot[INTERNAL_Y_ID], label=selected_samples[i]["label"])

    plot_y_max = max(plot[INTERNAL_Y_ID].max(), plot_y_max)  # absolute maximum
    plot_y_min = min(plot[INTERNAL_Y_ID].min(), plot_y_min)  # absolute minimum

    if selected_samples[i]["fit"]:  # fit peak if activated
        x_fit, y_fit, peak = fit_peak(plot_data[i])
        plt.plot(x_fit, y_fit, 'r-', label=f"Peak: {peak: .2f} nm")

plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance (a.u.)')

plt.xlim(WAVELENGTH_MIN, WAVELENGTH_MAX)
plot_y_delta = plot_y_max - plot_y_min
# add plot_y_delta * GRAPH_Y_MARGIN to the Y max and min in the graph as margin
plt.ylim(plot_y_min - (plot_y_delta * GRAPH_Y_MARGIN), plot_y_max + (plot_y_delta * GRAPH_Y_MARGIN))

plt.title(f'UV-Vis {"spectrum" if len(selected_samples) == 1 else "spectra"}')
plt.grid(True)
plt.legend()
plt.show()
