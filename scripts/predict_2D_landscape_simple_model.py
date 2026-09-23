import csv
from pathlib import Path

import numpy as np


# ============================================================================
# USER INPUT
# ============================================================================

# Apparent Kd for the switch-binder interaction (µM)
KD_APP = 22.9

# Apparent Kd for the binder-switch-release interaction (µM)
KD_COMP = 2.05

# Minimum normalized FRET signal
FMIN = 0.114435327681147

# Total binder concentrations (µM)
BINDER_CONCENTRATIONS = [
    0.1,
    1,
    2,
    5,
    10,
    20,
    50,
    150,
]

# Total switch-release concentrations (µM)
SWITCH_RELEASE_CONCENTRATIONS = [
    0.1,
    1,
    2,
    5,
    10,
    20,
    50,
    150,
]

# Output directory
OUTPUT_DIR = "predictions"


# ============================================================================
# EQUILIBRIUM MODEL
# ============================================================================

def calculate_complex_concentration(
    a_total,
    b_total,
    kd,
):
    """
    Calculate the equilibrium concentration of a 1:1 complex.

    Equilibrium:
        A + B <-> AB

    Parameters
    ----------
    a_total : float
        Total concentration of species A.

    b_total : float
        Total concentration of species B.

    kd : float
        Dissociation constant.

    Returns
    -------
    float
        Equilibrium concentration of AB.
    """
    discriminant = (
        (a_total + b_total + kd) ** 2
        - 4.0 * a_total * b_total
    )

    # Protect against very small negative values caused by
    # floating-point numerical precision.
    discriminant = max(discriminant, 0.0)

    return (
        a_total
        + b_total
        + kd
        - np.sqrt(discriminant)
    ) / 2.0


def calculate_heatmap(
    kd_app,
    kd_comp,
    binder_concentrations,
    switch_release_concentrations,
    fmin,
):
    """
    Calculate a predicted 2D switching landscape using a simplified
    sequential equilibrium model.

    The model first calculates binder sequestration by switch-release
    and then uses the resulting free binder concentration to predict
    the normalized switching response.

    Parameters
    ----------
    kd_app : float
        Apparent dissociation constant for switch-binder interaction.

    kd_comp : float
        Apparent dissociation constant for binder-switch-release
        interaction.

    binder_concentrations : array-like
        Total binder concentrations.

    switch_release_concentrations : array-like
        Total switch-release concentrations.

    fmin : float
        Minimum normalized FRET signal.

    Returns
    -------
    heatmap : numpy.ndarray
        Predicted normalized FRET response.

    complex_matrix : numpy.ndarray
        Equilibrium concentration of binder-switch-release complex.

    free_binder_matrix : numpy.ndarray
        Equilibrium concentration of free binder.
    """
    binder_concentrations = np.asarray(
        binder_concentrations,
        dtype=float,
    )

    switch_release_concentrations = np.asarray(
        switch_release_concentrations,
        dtype=float,
    )

    heatmap = np.zeros(
        (
            len(switch_release_concentrations),
            len(binder_concentrations),
        )
    )

    complex_matrix = np.zeros_like(heatmap)
    free_binder_matrix = np.zeros_like(heatmap)

    for i, switch_release_total in enumerate(
        switch_release_concentrations
    ):

        for j, binder_total in enumerate(
            binder_concentrations
        ):

            # ----------------------------------------------------------------
            # Binder sequestration by switch-release
            # ----------------------------------------------------------------

            complex_concentration = (
                calculate_complex_concentration(
                    binder_total,
                    switch_release_total,
                    kd_comp,
                )
            )

            free_binder = (
                binder_total
                - complex_concentration
            )

            # ----------------------------------------------------------------
            # Predicted normalized switching response
            # ----------------------------------------------------------------

            predicted_response = (
                fmin
                + (1.0 - fmin)
                / (1.0 + free_binder / kd_app)
            )

            complex_matrix[i, j] = complex_concentration
            free_binder_matrix[i, j] = free_binder
            heatmap[i, j] = predicted_response

    return (
        heatmap,
        complex_matrix,
        free_binder_matrix,
    )


# ============================================================================
# OUTPUT
# ============================================================================

def save_matrix_as_csv(
    filename,
    matrix,
    row_values,
    column_values,
):
    """
    Save a matrix to CSV with switch-release concentrations as rows
    and binder concentrations as columns.
    """
    with open(
        filename,
        "w",
        newline="",
    ) as file:

        writer = csv.writer(file)

        writer.writerow(
            ["Switch-release / Binder"]
            + list(column_values)
        )

        for row_value, row in zip(
            row_values,
            matrix,
        ):
            writer.writerow(
                [row_value]
                + list(row)
            )


# ============================================================================
# MAIN
# ============================================================================

def main():

    (
        heatmap,
        complex_matrix,
        free_binder_matrix,
    ) = calculate_heatmap(
        kd_app=KD_APP,
        kd_comp=KD_COMP,
        binder_concentrations=BINDER_CONCENTRATIONS,
        switch_release_concentrations=SWITCH_RELEASE_CONCENTRATIONS,
        fmin=FMIN,
    )

    np.set_printoptions(
        precision=4,
        suppress=True,
    )

    print("\nPredicted normalized FRET landscape:")
    print(heatmap)

    print("\nBinder-switch-release complex [BR]:")
    print(complex_matrix)

    print("\nFree binder [B]:")
    print(free_binder_matrix)

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    save_matrix_as_csv(
        output_dir / "predicted_heatmap_simple_model.csv",
        heatmap,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "BR_matrix_simple_model.csv",
        complex_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "free_binder_matrix_simple_model.csv",
        free_binder_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    print("\nResults saved to:")
    print(
        f"  {output_dir / 'predicted_heatmap_simple_model.csv'}"
    )
    print(
        f"  {output_dir / 'BR_matrix_simple_model.csv'}"
    )
    print(
        f"  {output_dir / 'free_binder_matrix_simple_model.csv'}"
    )


if __name__ == "__main__":
    main()
