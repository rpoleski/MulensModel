# MulensModel

<dl>MulensModel is package for modeling microlensing (or &mu;-lensing) 
events. </dl>

[Latest release: 1.14.0](https://github.com/rpoleski/MulensModel/releases/latest) and we're working on further developing the code. The documentations is at [**https://rpoleski.github.io/MulensModel/**](https://rpoleski.github.io/MulensModel/).

MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best fit parameters.

If you want to learn more about microlensing, please visit [Microlensing Source website](http://microlensing-source.org/).

Currently, MulensModel supports:
* Lens Systems: point lens or binary lens.
* Source Stars: single source or binary source.
* Effects: finite source (1-parameter), parallax (satellite & annual), binary lens orbital motion, different parametrizations of microlensing models.

Need more? Open [an issue](https://github.com/rpoleski/MulensModel/issues) or send us an e-mail. 

Are you using MulensModel for scientific research? Please give us credit by citing the [paper published in "Astronomy and Computing"](https://ui.adsabs.harvard.edu/abs/2019A%26C....26...35P/abstract) and [ASCL reference](http://ascl.net/1803.006). For arXiv see [link](https://arxiv.org/abs/1803.01003).


## Examples and tutorials

We have more than a dozen of examples - starting from very simple ones (like plotting a model) to very advanced (like fitting a binary lens model with finite source effect). See [**this file**](documents/examples_list.md) for **a list of examples** and a few tutorials. Most importantly, we have file describing **all microlensing parameters used in MulensModel** [**here**](documents/parameter_names.pdf). Also we have a file describing **methods used to calculate magnification in MulensModel** [**here**](documents/magnification_methods.pdf). 
Finally, we have **high-level fitting example** [**here**](examples/example_16).

The full documentation of API is at [https://rpoleski.github.io/MulensModel/](https://rpoleski.github.io/MulensModel/).

## How to install?

Download the source code and run:
```
python setup.py install
```
MulensModel requires some standard packages plus [astropy package](http://www.astropy.org/). To make sure that you have everything that's needed, just run:
```
pip install -r requirements.txt
```
Alternatively, you can run makefiles: go to `source/VBBL/` and run `make`, then go to `source/AdaptiveContouring/` and do the same. Then and add the path `MulensModel/source` to your `PYTHONPATH`. If you have any problems, please contact the authors and we will try to help.

If you want to **install MulensModel on Windows**, please see notes [here](documents/windows_install.md).

If you have **problems with installing or running MulensModel on MacOS**, please see notes [here](documents/macos_install.md).

---
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

file revised Jun 2020

