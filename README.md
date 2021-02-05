# Hydraulics2.jl

Hydraulic calculations for:

* Line sizing
* Pressure drop
* Networks

<!-- ABOUT THE PROJECT -->
## About The Project

This project provides functions in Julia for standard calculations for piping hydraulic calculations (line sizing and pressure drop calculation). The intention is to use the routines in a Jupyter Notebook file for documenting engineering work.  


### Built With

The code is written in Julia. The code is intended to be used in a Jupyter Notebook. I have not used the routines in a stand-alone Julia environment.



<!-- GETTING STARTED -->
## Getting Started

The following lines of code are needed in a Jupyter Notebook (Julia shell) to pull the package from GitHub and use the package.
~~~~
Pkg.add(PackageSpec(url="https://github.com/kevindorma/Hydraulics2.jl")
using Hydraulics
~~~~

### Prerequisites

The package requires the following packages
* DataFrames
* CSV
* ExcelFiles

<!-- TESTING -->
### Testing

The following code tests are available

* friction factor for a pipe

<!-- USAGE EXAMPLES -->
## Usage

Refer to the following files for working examples:

* Line Sizing
	* demoHydraulicsSizingRev1.ipynb
	* sizingHydraulicsSampleRev1.xls
* Simple hydraulics (pressure drop)
	* demoHydraulicsSimpleRev2.ipynb
	* simpleHydraulicsSampleRev2.xls
* Network hydraulics
	* demoHydraulicsNetworkRev1.ipynb
	* networkHydraulicsSampleRev1.xls

A spreadsheet is very helpful for tabulating the input information.

Refer to the ./docs/src/index.md for documentation on functions.



<!-- ROADMAP -->
## Roadmap

* Create piping network diagram from the information used to perform network calculations. This is useful for checking the layout of the network against an existing drawing.
* Implement gravity in the calculations.
* Permit negative mass flow rates.
* Implement different friction factors.
* Provide the ability to use equivalent length instead of fitting K value, useful for firewater networks.


<!-- CONTRIBUTING -->
## Contributing

Send me a note.



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Kevin Dorma - [@KevinDorma](https://twitter.com/KevinDorma) - kevin@kevindorma.ca

Project Link: [https://github.com/kevindorma/Hydraulics.jl](https://github.com/kevindorma/Hydraulics.jl)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

Not sure who to acknowledge.
