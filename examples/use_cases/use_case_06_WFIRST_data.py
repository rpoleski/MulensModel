import MulensModel

wfirst_data = MulensModel.MulensData(
    file_name='interesting_event.dat', satellite='WFIRST') 
print(wfirst_data.bandpass)
#automatically sets bandpass='W146'
