from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "/path/to/your/input_file.xlsx"

SHEET_NAME = "Sheet1"

X_COLUMN = "Cycle"
Y_COLUMN = "Emission 529nm"
Y_ERROR_COLUMN = "stdv"

X_LABEL = "Cycle"
Y_LABEL = "FRET intensity at 529 nm (a.u.)"
TITLE = "Reversible switching cycles"

PLOT_COLOR = "#0080FE"

LINE_WIDTH = 2
MARKER_SIZE = 8
ERROR_CAP_SIZE = 4


# ============================================================================
# LOAD DATA
# ============================================================================

file_path = Path(FILE_PATH)

if not file_path.exists():
    raise FileNotFoundError(
        f"Input file not found: {file_path}"
    )

df = pd.read_excel(
    file_path,
    sheet_name=SHEET_NAME,
)


# ============================================================================
# VALIDATE INPUT
# ============================================================================

required_columns = [
    X_COLUMN,
    Y_COLUMN,
    Y_ERROR_COLUMN,
]

missing_columns = [
    column
    for column in required_columns
    if column not in df.columns
]

if missing_columns:
    raise ValueError(
        "The following required columns were not found: "
        + ", ".join(missing_columns)
    )

for column in required_columns:
    df[column] = pd.to_numeric(
        df[column],
        errors="coerce",
    )

valid = (
    df[X_COLUMN].notna()
    & df[Y_COLUMN].notna()
    & df[Y_ERROR_COLUMN].notna()
)

df = df.loc[valid].copy()

if df.empty:
    raise ValueError(
        "No valid numeric data points were found."
    )


# ============================================================================
# PLOT
# ============================================================================

fig, ax = plt.subplots(figsize=(6, 4))

ax.errorbar(
    df[X_COLUMN],
    df[Y_COLUMN],
    yerr=df[Y_ERROR_COLUMN],
    fmt="o-",
    color=PLOT_COLOR,
    markerfacecolor=PLOT_COLOR,
    linewidth=LINE_WIDTH,
    markersize=MARKER_SIZE,
    capsize=ERROR_CAP_SIZE,
)

ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)

if TITLE:
    ax.set_title(TITLE)

ax.grid(False)

fig.tight_layout()

plt.show()
