def debye_huckel_alpha(temperature):
    """
    Approximation of the temperature dependency of ionic strength effects
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    Returns
    -------
    float
    """
    return 1e-3*(9.20483 * temperature) - 1e-5 * (1.284668 * temperature ** 2) + 1e-8 * (4.95199 * temperature ** 3)


def debye_huckel(ionic_strength, temperature, beta=1.6):
    """
    Debye-Huckel

    Parameters
    ----------
    ionic_strength : float
    temperature : float
        Temperature in Kelvin

    Returns
    -------
    """
    return debye_huckel_alpha(temperature) * ionic_strength ** 0.5 / (1.0 + beta * ionic_strength ** 0.5)
