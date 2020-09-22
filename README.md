# Raphael
PVT Simulator
16-09-2020

This is a PVT Simulator written in Python. It can be run using Jupyter Notebook or directly using the Python Command.

It requires Reservoir Fluid Composition, Reservoir Temperature and list of Pressures.
It uses the Peng Robinson 2 Parameter Equation of State to determine Vapour, Liquid quantities and compostions.
The final output is PVT tables in Eclipse Format PVTO, PVDG for Oils and PVTO, PVTG for Gas Condensates.

Caution:
Some Formatting may be necessary for the PVTO, PVDG and PVTG tables.
Viscosity uses LBC, but does not psudoize C7+

Future Work:
Improve output formatting
Include Pedersen Equation for calculating Viscosity

Please see the accompanying Presntation before use and let me have any comments.
