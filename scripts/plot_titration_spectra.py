from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "/path/to/your/input_file.xlsx"

X_COLUMN = "Wavelength"

SERIES = [
    ("0uM", "stdv1", "#1f77b4"),
    ("0.1uM", "stdv2", "#ff7f0e"),
    ("0.2uM", "stdv3", "#2ca02c"),
    ("0.5uM", "stdv4", "#d62728"),
    ("1uM", "stdv5", "#9467bd"),
    ("2uM", "stdv6", "#8c564b"),
    ("5uM", "stdv7", "#e377c2"),
    ("10uM", "stdv8", "#7f7f7f"),
    ("20uM", "stdv9", "#bcbd22"),
    ("50uM", "stdv10", "#17becf"),
    ("75uM", "stdv11", "#003f5c"),
    ("100uM", "stdv12", "#7a5195"),
    ("125uM", "stdv13", "#ef5675"),
    ("150uM", "stdv14", "#ffa600"),
]


# ============================================================================
# LOAD DATA
# ============================================================================

file_path = Path(FILE_PATH)

if not file_path.exists():
    raise FileNotFoundError(f"Input file not found: {file_path}")

df = pd.read_excel(file_path)

print("Available columns:")
print(df.columns.tolist())

print("\nFirst five rows:")
print(df.head())


# ============================================================================
# VALIDATE X-AXIS COLUMN
# ============================================================================

if X_COLUMN not in df.columns:
    raise ValueError(
        f"Required column '{X_COLUMN}' was not found in the input file."
    )

df[X_COLUMN] = pd.to_numeric(df[X_COLUMN], errors="coerce")

if df[X_COLUMN].isna().all():
    raise ValueError(
        f"Column '{X_COLUMN}' contains no valid numeric values."
    )


# ============================================================================
# VALIDATE DATA SERIES
# ============================================================================

valid_series = []
missing_series = []

for value_column, std_column, color in SERIES:
    if value_column not in df.columns or std_column not in df.columns:
        missing_series.append((value_column, std_column))
    else:
        valid_series.append((value_column, std_column, color))

if missing_series:
    print("\nMissing value/standard-deviation column pairs:")

    for value_column, std_column in missing_series:
        print(
            f"  - value: '{value_column}' "
            f"(found: {value_column in df.columns}), "
            f"standard deviation: '{std_column}' "
            f"(found: {std_column in df.columns})"
        )

if not valid_series:
    raise ValueError(
        "No valid data series were found in the input file."
    )


# ============================================================================
# CONVERT DATA TO NUMERIC
# ============================================================================

for value_column, std_column, _ in valid_series:

    df[value_column] = pd.to_numeric(
        df[value_column],
        errors="coerce",
    )

    df[std_column] = pd.to_numeric(
        df[std_column],
        errors="coerce",
    )

    if df[value_column].isna().all():
        print(
            f"Warning: '{value_column}' contains no valid numeric values."
        )

    if df[std_column].isna().all():
        print(
            f"Warning: '{std_column}' contains no valid numeric values."
        )


# ============================================================================
# PLOT
# ============================================================================

fig, ax = plt.subplots(figsize=(9, 6))

plotted_any = False

for value_column, std_column, color in valid_series:

    if df[value_column].dropna().empty:
        print(
            f"Skipping '{value_column}': "
            "the column contains no valid data."
        )
        continue

    ax.errorbar(
        df[X_COLUMN],
        df[value_column],
        yerr=df[std_column],
        fmt="o-",
        capsize=4,
        color=color,
        label=value_column,
        markersize=0.1,
        linewidth=4,
    )

    plotted_any = True

if not plotted_any:
    raise ValueError("No valid series could be plotted.")


# ============================================================================
# FIGURE FORMATTING
# ============================================================================

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Fluorescence")
ax.set_title("Switch Titration")

ax.legend(
    loc="best",
    fontsize="small",
)

ax.grid(
    False
)

ax.set_xlim(
    480,
    df[X_COLUMN].max(),
)

fig.tight_layout()

plt.show()
