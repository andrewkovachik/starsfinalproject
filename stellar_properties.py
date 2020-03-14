import numpy as np
import desolver as de

class Star:
    stellar_struct = {
        "density": de.DifferentialEquation("Density"),
        "temperature1": de.DifferentialEquation("Temperature1"),
        "temperature2": de.DifferentialEquation("Temperature2"),
        "mass": de.DifferentialEquation("Mass"),
        "luminosity": de.DifferentialEquation("Luminosity"),
        "opticaldepth": de.DifferentialEquation("Optical Depth")
        }

    def __init__(self):
        return

