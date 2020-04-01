"""
Takes in intial guesses for central temperature and pressure
as well as the core type and uses them to create a star and save
a text file.
"""
import stellar_properties as starprop

def make_star(central_temperature, central_density, core_type, name):

    star = starprop.Star(
            cent_density=central_density,
            cent_temperature=central_temperature,
            core=core_type,
            name=name)
