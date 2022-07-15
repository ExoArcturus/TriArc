# TriArc
Bayesian atmospheric retrieval for assessing detectability of trace gases in exoplanet atmospheres.
author: A. B. Claringbold (ExoArcturus)

Built on petitRADTRANS (https://petitradtrans.readthedocs.io/en/latest/)
and utilising PandExo (https://natashabatalha.github.io/PandExo/)

Paper for citation will be available soon, in the meantime please contact me if you would like to use TriArc (abc38@cam.ac.uk).

Requires:
petitRADTRANS,
PandExo (PandArc.py only)

Tutorial to come, comments in code should provide assistance. Also see test_script.py and Claringbold2022ResultsCode.py. A super-quick walkthrough is as follows. There are:

Three basic high level functions: plot_spectrum, perform_retrieval, and assess_detectibility.
The back_atm argument of these must be provided using load_planet() taking the argument of the name imported data file, like Hycean.py (see test_script.py)

Two macro functions: run_band and run_atmosphere.
These require an addition input file as the band_input_data as an argument, band_input_data which must also be imported e.g. prebiosignatures.py (see Claringbold2022ResultsCode.py). These functions have load_planet() built in, so only need the data file as the argument for atmosphere_background_data.

To use JWST precision provided by PandExo, will also need simulate_JWST_noise from PandArc.
Simply use simulate_JWST_noise with a third input file as the argument detailing observation parameters e.g. GJ1132b.py (see Claringbold2022ResultsCode.py).
This will create the noise profiles in .txt files which will be read in by the run_band and run_atmosphere functions in TriArc with JWST_noise = True.
