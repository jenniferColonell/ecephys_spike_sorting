CatGT Helper
==============
Python wrapper for CatGT, a C++ command application written by Bill Karsh for preprocessing data collected by SpikeGLX. CatGT can also scan the Neuropixels sync signal and auxiliary data to extract pulse edges for synchronization across data streams.  See the README for CatGT for details about parameters.

Dependencies
------------
[CatGT](https://billkarsh.github.io/SpikeGLX/#catgt)

Running
-------
```
python -m ecephys_spike_sorting.modules.catGT_helper--input_json <path to input json> --output_json <path to output json>
```
Two arguments must be included:
1. The location of an existing file in JSON format containing a list of paths and parameters.
2. The location to write a file in JSON format containing information generated by the module while it was run.

See the `_schemas.py` file for detailed information about the contents of the input JSON.

Input data
----------
- **SpikeGLX .bin files for ap and ni streams** : int16 binary files written by [SpikeGLX](https://github.com/billkarsh/spikeglx)

Output data
-----------
- **CatGT output files** : .bin files of concatenated, filtered data
- **CatGT edge files** : text files of edges found by CatGT scanning SYNC, XA and XD channels as specified in the CatGT command line.