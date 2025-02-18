
from py_round_number import round_number

# Changing from Ellipticity to Axis ratio
def ell_to_axrat(ell):
    axrat = 1 - ell
    return round_number(axrat,3)

def axrat_to_ell(axrat):
    ell = 1 - axrat
    return round_number(ell,3)