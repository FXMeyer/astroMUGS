import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import FortranFile


def grid(radlist, nb_sizes, itime, species, chempath):

    columns_static = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a']
    grid_dict = {}  # Dictionary to store z arrays

    for radius in radlist[0]:
        static = pd.read_table(f'{chempath}{radius}AU/1D_static.dat', delimiter=r"\s+", comment='!', header=None, engine='python')
        static.columns = columns_static
        z = static['z'].values  # Extract the z column as a NumPy array
        nh = static['nH'].values  # Extract the nH column as a NumPy array
        tg = static['Tg'].values  # Extract the Tg column as a NumPy array

        if nb_sizes == 1:
            td = static['Td'].values
            ab_d = static['inv_ab'].values
            grain_size = static['a'].values
            
            abundance = np.loadtxt(f'{chempath}{radius}AU/ab/{species}.ab', comments='!')
            abundance = np.delete(abundance[itime], 0)
            chem_numdens = nh*abundance

            # Store z, nH, and Tg arrays in a nested dictionary for the current radius
            grid_dict[radius] = {
                'z': z,
                'nH': nh,
                'Tg': tg,
                'nd': ab_d,  # Store the extracted grain sizes
                'Td': td,  # Store the extracted grain sizes
                'numberdens_species' : chem_numdens
            }

        elif nb_sizes > 1:
            # Check if the '1D_grain_sizes.in' file exists
            grain_file_path = f'{chempath}{radius}AU/1D_grain_sizes.in'
            grain_sizes = None  # Default value if the file does not exist
            grain_densities = None
            grain_temperatures = None
            if os.path.exists(grain_file_path):
                # Read the file and extract the first nb_sizes columns
                grain_data = pd.read_table(grain_file_path, delimiter=r"\s+", comment='!', header=None, engine='python')
                grain_abundances = grain_data.iloc[:, nb_sizes:2*nb_sizes].values  # Extract the first nb_sizes columns as a NumPy array
                grain_temperatures = grain_data.iloc[:, 2*nb_sizes:3*nb_sizes].values

                grain_densities = nh / np.transpose(grain_abundances)  # Broadcasting division


            #---abundances---
            #ab_list =  network_species(chempath, radius)
            ab_folder = f'{chempath}{radius}AU/ab/'
            ab_list = [f for f in os.listdir(ab_folder) if f.endswith(".ab") and os.path.isfile(os.path.join(ab_folder, f))]

            ab_list = [
            os.path.splitext(f)[0]  # remove extension
            for f in os.listdir(ab_folder)
            if f.endswith(".ab") 
            and os.path.isfile(os.path.join(ab_folder, f))
            and 'space' not in f and 'l' not in f#get rid of all non-abundance names.
            ]

         
            abundance = np.loadtxt(f'{chempath}{radius}AU/ab/{species}.ab', comments='!')
            abundance = np.delete(abundance[itime], 0)
            chem_numdens = nh*abundance

            # Store z, nH, and Tg arrays in a nested dictionary for the current radius
            grid_dict[radius] = {
                'z': z,
                'nH': nh,
                'Tg': tg,
                'grain_numdens': grain_densities,  # Store the extracted grain sizes
                'grain_temperatures': np.transpose(grain_temperatures),  # Store the extracted grain sizes
                'numberdens_species' : chem_numdens
            }

    return grid_dict


def parameters(radlist, chempath):
    chem_parameters = {}  # Dictionary to store the parameters
    with open(f'{chempath}{radlist[0][0]}AU/parameters.in', 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if not line or line.startswith('!'):  # Skip blank lines or comments
                continue

            # Split the line into key and value at the '=' character
            if '=' in line:
                key, value = line.split('=', 1)  # Split only at the first '='
                key = key.strip()  # Remove whitespace around the key
                value = value.split('!')[0].strip()  # Remove comments after the value
                chem_parameters[key] = value  # Add the key-value pair to the dictionary

    return chem_parameters


def radii(chempath):
    radii = os.listdir (chempath) # get all files' and folders' names in the current directory
    radlist = []
    for radius in radii:
        if radius.endswith("AU"):
            rad_d = radius.replace('AU','')
            rad_d = int(rad_d)
            radlist.append(rad_d)

    radlist = sorted(radlist)
    return radlist

def network_species(chempath, radius):
    ab_folder = os.listdir (chempath+radius+"AU/ab/") # get all files' and folders' names in the current directory

    # List all .ab files in the ab/ directory for the given radius
    ab_list = [f for f in os.listdir(ab_folder) if f.endswith(".ab") and os.path.isfile(os.path.join(ab_folder, f))]


    ab_list = [
    os.path.splitext(f)[0]  # remove extension
    for f in os.listdir(ab_folder)
    if f.endswith(".ab") 
    and os.path.isfile(os.path.join(ab_folder, f))
    and 'space' not in f and 'l' not in f#get rid of all non-abundance names.
    ]

    return ab_list


def species_names(path):
    """Read the ordered species list from NMGC's species.out file.

    ``species.out`` is written by NMGC at the start of every run and
    records the exact order of the internal ``species_name`` array —
    gas species first (in ``gas_species.in`` order), then grain/surface
    species (in ``grain_species.in`` order). This order matches the
    second axis of the ``abundances.out`` binary output.

    Parameters
    ----------
    path : str
        Path to the directory containing ``species.out`` (i.e. the
        NMGC working / chemistry directory).

    Returns
    -------
    list of str
        Species names in the same order as the abundances array axis.
    """
    names = []
    with open(os.path.join(path, 'species.out'), 'r') as f:
        for line in f:
            # Format: (5(I4,")",1X,A11,1X)) → "   1) H              2) H+    ..."
            # Splitting on ")" gives: ["   1", " H         2", " H+    3", ..., " NAME  \n"]
            # parts[0] has only the first index (no name yet); parts[1:] each start with a name.
            parts = line.split(')')
            for part in parts[1:]:
                tokens = part.split()
                if tokens:
                    names.append(tokens[0])
    return names


def abundances_binary(path):
    """Read the NMGC binary abundance output file.

    Reads ``abundances.out`` produced by NMGC when compiled with the
    single-file append output format. All timesteps are stored as
    sequential Fortran unformatted records. ``spatial_resolution``,
    ``nb_species``, and ``nb_timesteps`` are all inferred automatically
    from the file contents and size — no user input required.

    Each timestep contains three records:

    1. ``current_time`` — scalar float64 [seconds]
    2. Physical state — ``(4 * spatial_resolution + 1)`` float64 values:
       gas temperature, dust temperature (grain 1), H number density,
       visual extinction (each of length ``spatial_resolution``), then
       the X-ray ionisation rate (scalar).
    3. Abundances — ``nb_species * spatial_resolution`` float64 values
       in Fortran column-major order.

    Parameters
    ----------
    path : str
        Path to the ``abundances.out`` binary file.

    Returns
    -------
    dict with keys:

    - ``time`` : ndarray, shape (nb_timesteps,) — time in seconds
    - ``gas_temperature`` : ndarray, shape (nb_timesteps, spatial_resolution) [K]
    - ``dust_temperature`` : ndarray, shape (nb_timesteps, spatial_resolution) [K]
    - ``H_number_density`` : ndarray, shape (nb_timesteps, spatial_resolution) [cm-3]
    - ``visual_extinction`` : ndarray, shape (nb_timesteps, spatial_resolution) [mag]
    - ``X_ionisation_rate`` : ndarray, shape (nb_timesteps,) [s-1]
    - ``abundances`` : ndarray, shape (nb_timesteps, nb_species, spatial_resolution)
    """
    # Auto-detect spatial_resolution, nb_species, and nb_timesteps
    # from the first timestep. Record 2 always has 4*R + 1 values by
    # construction, so R = (len(phys0) - 1) // 4 is unambiguous.
    with FortranFile(path, 'r') as f:
        f.read_reals(dtype=np.float64)          # time
        phys0 = f.read_reals(dtype=np.float64)  # physical state
        ab0   = f.read_reals(dtype=np.float64)  # abundances

    R          = (len(phys0) - 1) // 4
    nb_phys    = len(phys0)
    nb_species = len(ab0) // R

    # Bytes per timestep: each Fortran record = 4-byte marker + data + 4-byte marker
    bytes_per_ts = (
        (4 + 1           * 8 + 4) +   # time record
        (4 + nb_phys     * 8 + 4) +   # physical state record
        (4 + nb_species * R * 8 + 4)  # abundances record
    )
    nb_timesteps = os.path.getsize(path) // bytes_per_ts

    # Allocate output arrays
    times    = np.empty(nb_timesteps)
    phys_all = np.empty((nb_timesteps, nb_phys))
    ab_all   = np.empty((nb_timesteps, nb_species, R))

    with FortranFile(path, 'r') as f:
        for i in range(nb_timesteps):
            times[i]    = f.read_reals(dtype=np.float64)[0]
            phys_all[i] = f.read_reals(dtype=np.float64)
            ab_all[i]   = f.read_reals(dtype=np.float64).reshape(nb_species, R, order='F')

    return {
        'time':             times,
        'gas_temperature':  phys_all[:, :R],
        'dust_temperature': phys_all[:, R:2*R],
        'H_number_density': phys_all[:, 2*R:3*R],
        'visual_extinction':phys_all[:, 3*R:4*R],
        'X_ionisation_rate':phys_all[:, -1],
        'abundances':       ab_all,
    }


def rates_binary(path, spatial_resolution=1):
    """Read the NMGC binary reaction rates output file.

    Reads ``rates.out`` produced by NMGC. All timesteps are appended
    sequentially; ``nb_reactions`` and ``nb_timesteps`` are inferred
    from the file size.

    Each timestep contains one record: the ``reaction_rates_1D`` array
    of shape ``(spatial_resolution, nb_reactions)`` in Fortran
    column-major order.

    Parameters
    ----------
    path : str
        Path to the ``rates.out`` binary file.
    spatial_resolution : int, optional
        Number of spatial grid points (default 1 for 0D per-cell runs).

    Returns
    -------
    ndarray, shape (nb_timesteps, spatial_resolution, nb_reactions)
        Reaction rates for all timesteps.
    """
    R = spatial_resolution

    # Auto-detect nb_reactions from the first record
    with FortranFile(path, 'r') as f:
        rec0 = f.read_reals(dtype=np.float64)

    nb_reactions  = len(rec0) // R
    bytes_per_ts  = 4 + R * nb_reactions * 8 + 4
    nb_timesteps  = os.path.getsize(path) // bytes_per_ts

    rates = np.empty((nb_timesteps, R, nb_reactions))

    with FortranFile(path, 'r') as f:
        for i in range(nb_timesteps):
            rates[i] = f.read_reals(dtype=np.float64).reshape(R, nb_reactions, order='F')

    return rates


# def abundance(path, itime=47, species='CO'):
#     chemmodel = pd.read_table(path + 'disk_t{}.dat'.format(str(itime)), sep=" ", engine='python')
#     chemmodel.dropna(how='all',inplace=True)
#     chemmodel.reset_index(inplace=True)

#     nr = chemmodel['r'].nunique()
#     nz = int(len(chemmodel['r'])/nr)

#     r = chemmodel['r'].unique()
#     z = np.reshape(chemmodel['z'].values, (nr, nz))
#     z = np.fliplr(z)
#     z = np.transpose(z)
#     zz = np.flipud(z)

#     nH = np.reshape(chemmodel['nH'].values, (nr, nz))
#     nH = np.transpose(nH)
#     ab = np.reshape(chemmodel['ab({})'.format(species)].values, (nr, nz))
#     ab = np.transpose(ab)
#     dens_mol = ab*nH #number density cm-3

#     return r, zz, dens_mol