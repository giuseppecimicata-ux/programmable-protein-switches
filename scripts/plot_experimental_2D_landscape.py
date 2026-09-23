from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "path/to/your/input_file.xlsx"

USE_LOG_SCALE = True

SHEET_NAME = 0

X_LABEL = "Binder concentration (µM)"
Y_LABEL = "Switch-release concentration (µM)"
VALUE_LABEL = "Normalized activity (fraction ON)"
TITLE = "2D switching landscape"

COLORMAP = "cividis"
VMIN = 0.0
VMAX = 1.0


# ============================================================================
# AXIS EDGE CALCULATION
# ============================================================================

def log_edges(values):
    """
    Calculate cell edges for logarithmically spaced concentration values.

    The returned edges define the boundaries of the cells plotted by
    pcolormesh.
    """
    values = np.asarray(values, dtype=float)

    values = np.unique(values)
    values = np.sort(values)

    if len(values) == 0:
        raise ValueError("No concentration values were found.")

    if np.any(values <= 0):
        raise ValueError(
            "All concentration values must be greater than zero "
            "when using a logarithmic scale."
        )

    if len(values) == 1:
        value = values[0]
        return np.array([
            value / np.sqrt(10),
            value * np.sqrt(10),
        ])

    log_values = np.log(values)

    edges = np.empty(
        len(values) + 1,
        dtype=float,
    )

    edges[1:-1] = (
        log_values[:-1]
        + log_values[1:]
    ) / 2.0

    edges[0] = (
        log_values[0]
        - (edges[1] - log_values[0])
    )

    edges[-1] = (
        log_values[-1]
        + (log_values[-1] - edges[-2])
    )

    return np.exp(edges)


def linear_edges(values):
    """
    Calculate cell edges for linearly spaced concentration values.
    """
    values = np.asarray(values, dtype=float)

    values = np.unique(values)
    values = np.sort(values)

    if len(values) == 0:
        raise ValueError("No concentration values were found.")

    if len(values) == 1:
        value = values[0]
        return np.array([
            value - 0.5,
            value + 0.5,
        ])

    edges = np.empty(
        len(values) + 1,
        dtype=float,
    )

    edges[1:-1] = (
        values[:-1]
        + values[1:]
    ) / 2.0

    edges[0] = (
        values[0]
        - (edges[1] - values[0])
    )

    edges[-1] = (
        values[-1]
        + (values[-1] - edges[-2])
    )

    return edges


# ============================================================================
# LOAD AND PLOT EXPERIMENTAL LANDSCAPE
# ============================================================================

def plot_experimental_landscape(
    file_path,
    sheet_name=0,
    x_label="Binder concentration (µM)",
    y_label="Switch-release concentration (µM)",
    value_label="Normalized activity (fraction ON)",
    title="2D switching landscape",
    cmap="cividis",
    vmin=0.0,
    vmax=1.0,
):
    """
    Plot an experimental two-dimensional switching landscape.

    The input Excel sheet is expected to contain:
        - binder concentrations in the first row;
        - switch-release concentrations in the first column;
        - normalized activity values in the remaining matrix.
    """
    df = pd.read_excel(
        file_path,
        sheet_name=sheet_name,
        header=None,
    )

    # The current experimental matrix uses 8 x 8 concentration combinations,
    # plus one header row and one header column.
    if df.shape[0] < 9 or df.shape[1] < 9:
        raise ValueError(
            f"Expected at least a 9 x 9 worksheet, "
            f"found {df.shape[0]} x {df.shape[1]}."
        )

    x_values = pd.to_numeric(
        df.iloc[0, 1:9],
        errors="coerce",
    ).to_numpy(dtype=float)

    y_values = pd.to_numeric(
        df.iloc[1:9, 0],
        errors="coerce",
    ).to_numpy(dtype=float)

    z_values = (
        df.iloc[1:9, 1:9]
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy(dtype=float)
    )

    if np.any(~np.isfinite(x_values)):
        raise ValueError(
            "The binder concentration labels contain invalid values."
        )

    if np.any(~np.isfinite(y_values)):
        raise ValueError(
            "The switch-release concentration labels contain invalid values."
        )

    if np.any(~np.isfinite(z_values)):
        raise ValueError(
            "The activity matrix contains invalid values."
        )

    # Sort concentration axes and rearrange the activity matrix
    # accordingly.
    x_order = np.argsort(x_values)
    y_order = np.argsort(y_values)

    x_values = x_values[x_order]
    y_values = y_values[y_order]

    z_values = z_values[
        np.ix_(y_order, x_order)
    ]

    # Calculate cell boundaries according to the selected axis scale.
    if USE_LOG_SCALE:
        x_edges = log_edges(x_values)
        y_edges = log_edges(y_values)
    else:
        x_edges = linear_edges(x_values)
        y_edges = linear_edges(y_values)

    # ------------------------------------------------------------------------
    # Figure style
    # ------------------------------------------------------------------------

    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "text.color": "black",
        "font.size": 11,
    })

    fig, ax = plt.subplots(
        figsize=(6.6, 5.8)
    )

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        z_values,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
        edgecolors="none",
    )

    if USE_LOG_SCALE:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if title:
        ax.set_title(title)

    ax.set_xlim(
        x_edges[0],
        x_edges[-1],
    )

    ax.set_ylim(
        y_edges[0],
        y_edges[-1],
    )

    colorbar = fig.colorbar(
        mesh,
        ax=ax,
    )

    colorbar.set_label(value_label)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(False)

    fig.tight_layout()

    plt.show()


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":

    file_path = Path(FILE_PATH)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {file_path}"
        )

    plot_experimental_landscape(
        file_path=file_path,
        sheet_name=SHEET_NAME,
        x_label=X_LABEL,
        y_label=Y_LABEL,
        value_label=VALUE_LABEL,
        title=TITLE,
        cmap=COLORMAP,
        vmin=VMIN,
        vmax=VMAX,
    )
