import sys
if sys.version_info.major > 2:
    from subprocess import STDOUT, run as sp_run
else:
    from subprocess32 import STDOUT, run as sp_run


def run(chempath, nmgc_cmd='nmgc', verbose=True, timelimit=None):
    """Run an NMGC gas-grain chemistry simulation.

    Calls ``nmgc run`` as a subprocess in the given chemistry directory,
    which must already contain all required input files (parameters.in,
    species files, reaction files, etc.). Typically called from
    ``Model.run_chemistry()`` rather than directly, so that ``chempath``
    and ``nmgc_cmd`` are taken from the model configuration.

    Parameters
    ----------
    chempath : str
        Path to the directory containing NMGC input files. NMGC writes
        its outputs (abundances.out, rates.out) in this directory.
    nmgc_cmd : str, optional
        Command used to invoke NMGC. Override to point to a specific
        binary path if ``nmgc`` is not on PATH.
    verbose : bool, optional
        If True, NMGC stdout/stderr stream to the terminal. If False,
        output is redirected to ``<chempath>/nmgc.out``.
    timelimit : float or None, optional
        Timeout in seconds. Raises ``subprocess.TimeoutExpired`` if
        exceeded. None means no timeout.

    Returns
    -------
    subprocess.CompletedProcess
        Result of the subprocess call. Check ``result.returncode == 0``
        to confirm success.

    Raises
    ------
    subprocess.TimeoutExpired
        If the simulation exceeds ``timelimit``.
    FileNotFoundError
        If ``nmgc_cmd`` is not found on PATH.
    RuntimeError
        If NMGC exits with a non-zero return code.
    """
    command = '{} run'.format(nmgc_cmd)

    if verbose:
        result = sp_run(command.split(), cwd=chempath, stderr=STDOUT, timeout=timelimit)
    else:
        with open('{}/nmgc.out'.format(chempath), 'w') as f:
            result = sp_run(command.split(), cwd=chempath, stdout=f, stderr=f, timeout=timelimit)

    if result.returncode != 0:
        raise RuntimeError(
            'NMGC exited with code {}. Check {}/nmgc.out for details.'.format(
                result.returncode, chempath
            )
        )

    return result
