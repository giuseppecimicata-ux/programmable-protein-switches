from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "/path/to/your/input_file.xlsx"

SHEET_NAME = 0
OUTPUT_PATH = "measured_vs_predicted_landscape.png"

X_LABEL = "Measured values"
Y_LABEL = "Calculated values"
TITLE = "Measured vs calculated switching landscape"

POINT_COLOR = "black"
POINT_SIZE = 35
LINE_COLOR = "red"
LINE_WIDTH = 2

# Set to True if the Excel file contains column headers.
HAS_HEADER = True

# Column indices containing measured and calculated values.
MEASURED_COLUMN = 0
CALCULATED_COLUMN = 1


# ============================================================================
# METRICS
# ============================================================================

def r2_score(y_true, y_pred):
    """
    Calculate the coefficient of determination R².
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    ss_res = np.sum(
        (y_true - y_pred) ** 2
    )

    ss_tot = np.sum(
        (y_true - np.mean(y_true)) ** 2
    )

    if ss_tot == 0:
        return np.nan

    return 1.0 - ss_res / ss_tot


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():

    file_path = Path(FILE_PATH)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {file_path}"
        )

    if HAS_HEADER:
        df = pd.read_excel(
            file_path,
            sheet_name=SHEET_NAME,
        )
    else:
        df = pd.read_excel(
            file_path,
            sheet_name=SHEET_NAME,
            header=None,
        )

    required_column = max(
        MEASURED_COLUMN,
        CALCULATED_COLUMN,
    )

    if df.shape[1] <= required_column:
        raise ValueError(
            f"The input file must contain at least "
            f"{required_column + 1} columns."
        )

    measured = pd.to_numeric(
        df.iloc[:, MEASURED_COLUMN],
        errors="coerce",
    ).to_numpy(dtype=float)

    calculated = pd.to_numeric(
        df.iloc[:, CALCULATED_COLUMN],
        errors="coerce",
    ).to_numpy(dtype=float)

    valid = (
        np.isfinite(measured)
        & np.isfinite(calculated)
    )

    measured = measured[valid]
    calculated = calculated[valid]

    if len(measured) == 0:
        raise ValueError(
            "No valid measured/calculated data pairs were found."
        )

    r2 = r2_score(
        measured,
        calculated,
    )

    # Use common limits for both axes.
    data_min = min(
        np.min(measured),
        np.min(calculated),
    )

    data_max = max(
        np.max(measured),
        np.max(calculated),
    )

    data_range = data_max - data_min

    padding = 0.05 * (
        data_range
        if data_range > 0
        else 1.0
    )

    lower_limit = data_min - padding
    upper_limit = data_max + padding


    # =========================================================================
    # PLOT
    # =========================================================================

    fig, ax = plt.subplots(
        figsize=(5.5, 5.5)
    )

    ax.scatter(
        measured,
        calculated,
        s=POINT_SIZE,
        color=POINT_COLOR,
        alpha=0.85,
    )

    # Ideal agreement line: calculated = measured.
    ax.plot(
        [lower_limit, upper_limit],
        [lower_limit, upper_limit],
        color=LINE_COLOR,
        linewidth=LINE_WIDTH,
        label="y = x",
    )

    ax.set_xlim(
        lower_limit,
        upper_limit,
    )

    ax.set_ylim(
        lower_limit,
        upper_limit,
    )

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)

    if TITLE:
        ax.set_title(TITLE)

    ax.legend(
        title=f"R² = {r2:.3f}"
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(False)

    output_path = Path(OUTPUT_PATH)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    print(f"R² = {r2:.6f}")
    print(f"Figure saved to: {output_path}")


if __name__ == "__main__":
    main()
