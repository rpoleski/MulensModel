import MulensModel as mm


raise NotImplementedError('satellite keyword for MulensData not supported')

wfirst_data = mm.MulensData(
    file_name='interesting_event.dat', satellite='WFIRST')
print(wfirst_data.bandpass)
# automatically sets bandpass='W146'
