


# Let users know if they are missing some dependencies
hard_dependencies = ['pandas']
soft_dependencies = ['ROOT']
missing_dependencies = []

for dependency in hard_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(dependency)

for dependency in soft_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        print('ImportError: {} is missing. Useability is highly restricted without!'.format(dependency))


from .import mipcali
from .import phase2
