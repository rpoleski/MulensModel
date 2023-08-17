# How to get satellite ephemeris using JPL Horizons service?

The [JPL Horizons service](https://ssd.jpl.nasa.gov/horizons/) allows retrieving 
highly accurate positions of satellites. We need these positions to calculate 
how microlensing parallax affects event magnification seen by the satellite. 
This file explains how to access Horizons via e-mail.

The MulensModel package contains short template e-mail: 

`documents/horizons_short_template.txt`

This file contains all the basic settings that Horizons requires. There are only 
a few settings that one will modify on normal basis. These are:

### Satellite choice

This is set by keyword `COMMAND` and has a numerical value. There are many 
natural and many artificial Solar System objects, but we are interested in 
just a few. Let's have a look at 
[this list](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html)
Some particular choices:

* Hubble Space Telescope `-48`
* Spitzer Space Telescope `-79`
* New Horizons `-98`
* EPOXI `-140`
* Rosetta `-226`
* Kepler/K2 `-227`

Oops... Gaia is not there! And yes, Gaia is observing the microlensing events. 
In that case let's go to 
[another website](https://ssd.jpl.nasa.gov/horizons/app.html#/). 
By querying it, we can find:

* Gaia `-139479`
* Swift `-128485`
* Euclid `-680`

Note that all the values are __negative__ and as long as we don't start observing 
microlensing events from other planets or Moon, we will only use negative 
values for keyword `COMMAND`

---

### Epochs selection

You can select epochs for which you want calculations to be performed in two 
ways: 

1. by selecting `START_TIME`, `STOP_TIME`, and `STEP_SIZE`, or
2. by providing a list of epochs via `TLIST` keyword. 

The first one should be obvious. If you want to use the second one, then see 
the manual linked below. These epochs don't have to be exactly the same 
as the epochs in which you have observations. If they don't match, 
then they'll be interpolated. However, keep in mind that the time step 
used for interpolation cannot be too short -- if you're using data from 
satellite on a low Earth orbit, then the step size cannot be 1 day! 
It is a good practice to extend epochs range to be wider than your data. 
This allows one to make nice plots that contain both ground-based and 
space-based data.

---

### Expected output

The `QUANTITIES` keyword is set to `'1,20,21'` and this translates to, 
respectively:

* 'Astrometric RA & DEC'
* 'Obsrv range & rng rate'
* 'One-Way Light-Time'

These are the quantities required by MulensModel. If you need more quantities 
for whatever reason, then add the additional ones at the end of the list. 
There are currently 39 other quantities you may ask for. If you wish, there 
are even more options available if you change `TABLE_TYPE` keyword to 
something other than `'OBS'`, but we're not discussing it here as MulensModel 
will not read those. 

---

If you need more detailed description, see:

[ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.brief](ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.brief)

or

[ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.long](ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.long)

Once you are done with all the settings, you're ready to send an e-mail. 
Just open your e-mail program and: 

 * type `horizons@ssd.jpl.nasa.gov` in To field of the e-mail, 
 * type `JOB` in Subject field,
 * copy-paste the settings prepared above into the main body of the e-mail,
 * make sure that your e-mail is in __plain text mode__ (if you send html message, then the Horizons won't be able to process it!),
 * send an e-mail,
 * wait for reply and save it once it arrives,
 * if you're happy with what you got, switch your e-mail program back to html mode, rich text mode etc.

Congratulations! You now know how to access Horizons via e-mail.

If you send default or only slightly modified settings, then you can compare the results with mine, which are here:

`data/ephemeris_files/Spitzer_ephemeris_01.dat`

---

(C) R. Poleski, Mar 2017, last changes: Aug 2023

