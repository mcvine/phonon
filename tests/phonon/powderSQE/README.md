Tests and examples here are used to compute phonon powder SQE (no Monte Carlo neutron ray tracing,
just simple modeling).

If FORCE_CONSTANTS and related files exist, it is better to use the use_phonopy module.
See use_phonopy_TestCase.py

If phonon data only exists in DANSE IDF, "IDF" module have to be used then.

There are examples and tests here use "IDF" module but uses phonopy to calculate DANSE IDF
data and then calculate SQE. This is not the recommended way. Use "use_phonopy" module
instead.
