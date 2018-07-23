# MulensModel support for WFIRST Data Analysis Challenge

The [WFIRST Data Analysis Challenge](http://microlensing-source.org/data-challenge/) is a great way to get involved in an awesome microlensing survey.  The [WFIRST satellite](https://en.wikipedia.org/wiki/Wide_Field_Infrared_Survey_Telescope) will deliver thousands of space-quality images of galactic bulge. These images will contain photometry for over ten thousand microlensing events with around 1,000 of them showing planetary signals. Undoubtedly, analyzing the data will be challenging and as a part of planning WFIRST, we need to find these challenges and also mitigate them. Hence, the WFIRST data challenge has been announced.

For the full description of the data challenge, see: [http://microlensing-source.org/data-challenge/](http://microlensing-source.org/data-challenge/).

Here, we present a set of examples and tools that may be helpful for data challenge participants. Please note that MulensModel does not have fitting capability itself. It provides a user-friendly method to calculate chi^2 and examples of how to use it. You have to make a number of choices while fitting the model to the data: what fitting algorithm to use, what is the starting point, when the process ends, have I explored all possible models?, etc. And answering these questions is your task.

### Tutorials

* [tutorial on basics of MulensModel](https://rpoleski.github.io/MulensModel/tutorial.html)
* [simple fitting of a single (star) lens](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html) and [more advanced fitting: including the parallax effect](https://rpoleski.github.io/MulensModel/tutorial_fit_pi_E.html)
* [fit a planetary model using grids of parameters](https://github.com/rpoleski/MulensModel/blob/master/examples/example_08_planet_grid_fitting.ipynb)

If you participate and you see that MulensModel is missing some important functionality, then please [open a GitHub issue](https://help.github.com/articles/creating-an-issue/) and we will try to help. We are always open to new ideas. 

### Frequently Asked Questions

* How do I read [wfirst\_ephemeris\_W149.txt and wfirst\_ephemeris\_Z087.txt files](https://github.com/microlensing-data-challenge/data-challenge-1)?

   When you create instance of MulensData class, simply set the keyword `ephemerides_file=PATH_TO_THE_FILE`.

* The data challenge requires specifying code versions. How do I get these?

   Convenient Method: Go to the directory where you installed MulensModel (the one above source/), then examples/, and run checks.py: `python checks.py`.

   General Method: To find which version of python you're using just type `import sys` and `print(sys.version)`.  For MulensModel use `print(MulensModel.__version__)` instead.  Similar commands should work for numpy, scipy, matplotlib, ctypes, and astropy.

