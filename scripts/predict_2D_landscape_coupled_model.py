import csv
from pathlib import Path

import numpy as np


# ============================================================================
# USER INPUT
# ============================================================================

# Apparent Kd for switch-binder interaction (µM)
KD_SB = 22.9

# Dissociation constant for switch-release interaction with the closed switch
# (µM)
KD_SR = 10.0

# Apparent Kd for binder-switch-release interaction (µM)
KD_BR = 2.05

# Minimum normalized FRET signal
FMIN = 0.114435327681147

# Total switch concentration (µM)
SWITCH_TOTAL = 1.0

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
# COUPLED EQUILIBRIUM MODEL
# ============================================================================

def calculate_free_ligands(
    binder_total,
    switch_release_total,
    kd_br,
):
    """
    Calculate free binder and free switch-release concentrations.

    The B-R equilibrium is:

        B + R <-> BR

    Under the approximation that the total switch concentration is much
    smaller than the total binder and switch-release concentrations, depletion
    of B and R by binding to the switch can be neglected in these mass
    balances.

    This leaves the B-R equilibrium to be solved from the total B and R
    concentrations.

    Returns
    -------
    b_free : float
        Free binder concentration.

    r_free : float
        Free switch-release concentration.

    br : float
        Equilibrium concentration of the B-R complex.
    """
    delta = switch_release_total - binder_total

    discriminant = (
        (delta + kd_br) ** 2
        + 4.0 * binder_total * kd_br
    )

    discriminant = max(discriminant, 0.0)

    b_free = (
        -(delta + kd_br)
        + np.sqrt(discriminant)
    ) / 2.0

    # Protect against tiny negative values caused by floating-point
    # numerical precision.
    b_free = max(b_free, 0.0)

    r_free = b_free + delta
    r_free = max(r_free, 0.0)

    br = binder_total - b_free
    br = max(br, 0.0)

    return b_free, r_free, br


def calculate_closed_fraction(
    b_free,
    r_free,
    kd_sb,
    kd_sr,
):
    """
    Calculate the equilibrium fraction of switch in the closed state.

    The closed switch participates in two ligand-binding equilibria:

        S* + B <-> SB
        S* + R <-> SR

    The closed-state fraction is:

        [S*] / [S]_tot
        = 1 / (1 + [B]/K_SB + [R]/K_SR)
    """
    return 1.0 / (
        1.0
        + b_free / kd_sb
        + r_free / kd_sr
    )


def calculate_heatmap(
    kd_sb,
    kd_sr,
    kd_br,
    binder_concentrations,
    switch_release_concentrations,
    fmin,
    switch_total,
):
    """
    Calculate the predicted 2D switching landscape using the coupled
    equilibrium model.

    Parameters
    ----------
    kd_sb : float
        Dissociation constant for switch-binder binding.

    kd_sr : float
        Dissociation constant for switch-release binding to the closed switch.

    kd_br : float
        Dissociation constant for binder-switch-release binding.

    binder_concentrations : array-like
        Total binder concentrations.

    switch_release_concentrations : array-like
        Total switch-release concentrations.

    fmin : float
        Minimum normalized FRET signal.

    switch_total : float
        Total switch concentration.

    Returns
    -------
    heatmap : numpy.ndarray
        Predicted normalized FRET response.

    br_matrix : numpy.ndarray
        Equilibrium concentration of the binder-switch-release complex.

    free_binder_matrix : numpy.ndarray
        Equilibrium concentration of free binder.

    free_switch_release_matrix : numpy.ndarray
        Equilibrium concentration of free switch-release.

    closed_fraction_matrix : numpy.ndarray
        Predicted fraction of switch in the closed state.
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

    br_matrix = np.zeros_like(heatmap)
    free_binder_matrix = np.zeros_like(heatmap)
    free_switch_release_matrix = np.zeros_like(heatmap)
    closed_fraction_matrix = np.zeros_like(heatmap)

    for i, switch_release_total in enumerate(
        switch_release_concentrations
    ):

        for j, binder_total in enumerate(
            binder_concentrations
        ):

            b_free, r_free, br = calculate_free_ligands(
                binder_total=binder_total,
                switch_release_total=switch_release_total,
                kd_br=kd_br,
            )

            closed_fraction = calculate_closed_fraction(
                b_free=b_free,
                r_free=r_free,
                kd_sb=kd_sb,
                kd_sr=kd_sr,
            )

            predicted_response = (
                fmin
                + (1.0 - fmin) * closed_fraction
            )

            br_matrix[i, j] = br
            free_binder_matrix[i, j] = b_free
            free_switch_release_matrix[i, j] = r_free
            closed_fraction_matrix[i, j] = closed_fraction
            heatmap[i, j] = predicted_response

    return (
        heatmap,
        br_matrix,
        free_binder_matrix,
        free_switch_release_matrix,
        closed_fraction_matrix,
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
                [row_value] + list(row)
            )


# ============================================================================
# MAIN
# ============================================================================

def main():

    (
        heatmap,
        br_matrix,
        free_binder_matrix,
        free_switch_release_matrix,
        closed_fraction_matrix,
    ) = calculate_heatmap(
        kd_sb=KD_SB,
        kd_sr=KD_SR,
        kd_br=KD_BR,
        binder_concentrations=BINDER_CONCENTRATIONS,
        switch_release_concentrations=SWITCH_RELEASE_CONCENTRATIONS,
        fmin=FMIN,
        switch_total=SWITCH_TOTAL,
    )

    np.set_printoptions(
        precision=4,
        suppress=True,
    )

    print("\nPredicted normalized FRET landscape:")
    print(heatmap)

    print("\nBinder-switch-release complex [BR]:")
    print(br_matrix)

    print("\nFree binder [B]:")
    print(free_binder_matrix)

    print("\nFree switch-release [R]:")
    print(free_switch_release_matrix)

    print("\nClosed switch fraction [S*]/[S]total:")
    print(closed_fraction_matrix)

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    save_matrix_as_csv(
        output_dir / "predicted_heatmap_coupled_model.csv",
        heatmap,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "BR_matrix_coupled_model.csv",
        br_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "free_binder_matrix_coupled_model.csv",
        free_binder_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "free_switch_release_matrix_coupled_model.csv",
        free_switch_release_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    save_matrix_as_csv(
        output_dir / "closed_fraction_coupled_model.csv",
        closed_fraction_matrix,
        SWITCH_RELEASE_CONCENTRATIONS,
        BINDER_CONCENTRATIONS,
    )

    print("\nResults saved to:")
    print(
        f"  {output_dir / 'predicted_heatmap_coupled_model.csv'}"
    )
    print(
        f"  {output_dir / 'BR_matrix_coupled_model.csv'}"
    )
    print(
        f"  {output_dir / 'free_binder_matrix_coupled_model.csv'}"
    )
    print(
        f"  {output_dir / 'free_switch_release_matrix_coupled_model.csv'}"
    )
    print(
        f"  {output_dir / 'closed_fraction_coupled_model.csv'}"
    )


if __name__ == "__main__":
    main()
