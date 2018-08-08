"""
Modules here are used to compute phonon powder SQE (no Monte Carlo neutron ray tracing,
just simple modeling).

If FORCE_CONSTANTS and related files exist, it is better to use the use_phonopy module.
See tests/phonon/powderSQE/use_phonopy_TestCase.py

If phonon data only exists in DANSE IDF, "IDF" module is the only choice, but please 
keep in mind it may not be doing things correctly for anything other than cubic systems.

There are examples and tests here use "IDF" module but uses phonopy to calculate DANSE IDF
data and then calculate SQE. This is not the recommended way. Use "use_phonopy" module
instead.
"""
