import MulensModel
import matplotlib.pyplot as pl

data = MulensModel.MulensData(file_name="my_data.dat")
model = MulensModel.Model(n_components=1)


#pandas style. Does this work without pandas???
data.bad[np.isnan(data.err)] = True


event = MulensModel.Event(datasets=data, model=model)
event.estimate_model_params() #aspirational
event.get_chi2() 
event.clean_data()
event.clean_data(sigma=3.) 
"""
if sigma==None, set sigma based on the number of data points. 
clean_data requires fitting a model...
"""
event.get_chi2() #should mask bad data

pl.scatter(data.time, data.mag, marker="o", facecolor=None)
pl.scatter(data.time[data.good], data.mag[not data.bad], marker="o", facecolor="black")
pl.scatter(data.time[data.bad], data.mag[data.bad], marker="x")

