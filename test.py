import os
import shutil

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.synthetic_population import Population
from posydon.popsyn.synthetic_population import PopulationRunner


path_to_params = os.path.join(PATH_TO_POSYDON, "posydon/popsyn/population_params_default.ini")
shutil.copyfile(path_to_params, 'data/population_params.ini')

poprun = PopulationRunner('data/population_params.ini', verbose=True)

print('Number of binary populations:',len(poprun.binary_populations))
print('Metallicity:', poprun.binary_populations[0].metallicity)
print('Number of binaries:', poprun.binary_populations[0].number_of_binaries)

poprun.evolve()

