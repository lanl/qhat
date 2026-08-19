# _________________________________________________________________________________________________
# Example configuration demonstrating command-line parameter usage
#
# This config can be run with parameters like:
#   python3.11 driver.py examples/cmdline_params_demo_config.py -p distance=1.5 -p basis=sto-6g
#
# Parameters are accessed via the 'params' dictionary that is provided to all configuration files.
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
# General configuration

general.print_verbose()
general.file_stub = "Be-H_cmdline_demo"
general.file_format = "default"

# _________________________________________________________________________________________________
# Demonstrate different ways to use command-line parameters

# Method 1: Use parameter with a default value (recommended)
distance = params.get('distance', 1.30)  # Default to 1.30 if not provided
basis_set = params.get('basis', 'sto-3g')  # Default to 'sto-3g' if not provided

# Method 2: Require parameter (will raise KeyError if not provided)
# Uncomment this to require the parameter:
# distance = params['distance']  # Will raise KeyError if -p distance=... not provided

# Method 3: Use parameter with type conversion and validation
num_occupied = int(params.get('num_occupied', 3))
num_vacant = int(params.get('num_vacant', 3))

# _________________________________________________________________________________________________
# Describe the Hamiltonian using the parameters

hamiltonian.add_atom("Be", 0.00, 0, 0)
hamiltonian.add_atom("H", distance, 0, 0)  # Uses the distance parameter

hamiltonian.basis = basis_set  # Uses the basis parameter

hamiltonian.num_active_occupied = num_occupied  # Uses the num_occupied parameter
hamiltonian.num_active_vacant = num_vacant      # Uses the num_vacant parameter

# _________________________________________________________________________________________________
# Parameters can also be used in more complex expressions

# Example: You could pass a list of atoms
# python3.11 driver.py config.py -p atoms='[("Be",0.0),("H",1.3)]'
# Then use it like:
# if 'atoms' in params:
#     for atom_type, position in params['atoms']:
#         hamiltonian.add_atom(atom_type, position, 0, 0)
