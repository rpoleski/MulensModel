# MulensModel support for WFIRST Data Analysis Challenge

[WFIRST Data Analysis Challenge](http://microlensing-source.org/data-challenge/) is a great way to get involved in awesome microlensing survey.  The [WFIRST satellite](https://en.wikipedia.org/wiki/Wide_Field_Infrared_Survey_Telescope) will deliver thousands of space-quality images of galactic bulge. These images will contain photometry for a dozen or so thousand microlensing events with around 1000 of them showing planetary signals. Undoubtedly, analysing the data will be challenging and as a part of planning WFIRST, we need to find these challenges and also mitigate them. Hence, the WFIRST data challenge has been announced.

For the full description of the data challenge, see: [http://microlensing-source.org/data-challenge/](http://microlensing-source.org/data-challenge/)

Here, we present a set of examples and tools that may be helpful for data challenge participants. If you want to participate and you see that MulensModel is missing some important functionality, then please [open a GitHub issue](https://help.github.com/articles/creating-an-issue/) and we will try to help. We are always open to new ideas. 

Please note that MulensModel has not fitting capability itself, it provides a user-friendly method to calculate chi^2 and examples how to use it. You have to make a number of choices while fitting the model to the data: what should be the fitting algorithm, what is the staring point, when the process ends, have I explored all possible models, etc. And answering these questions is your task.

### Links

* [tutorial on basics of MulensModel](https://rpoleski.github.io/MulensModel/tutorial.html)
* [simple fitting](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html) and [more advanced fitting](https://rpoleski.github.io/MulensModel/tutorial_fit_pi_E.html)
* [fit a planetary model using grids of parameters](https://github.com/rpoleski/MulensModel/blob/master/examples/example_08_planet_grid_fitting.ipynb)

### Frequently Asked Questions

* How do I read [wfirst\_ephemeris\_W149.txt and wfirst\_ephemeris\_Z087.txt files](https://github.com/microlensing-data-challenge/data-challenge-1)?

   When you create instance of MulensData class, then simply add option `ephemerides_file=PATH_TO_THE_FILE`.

* Data challenge requires specifying code versions. How do I get these?

   To find which version of python you're using just type `import sys` and `print(sys.version)`.  For MulensModel use `print(MulensModel.__version__)` instead.  Similar command should work for numpy, scipy, matplotlib, ctypes, and astropy.

