## High level decisions we should make:
(starting from the most important ones)

1. How to handle changes in origin of the coordinate system? Internally we're working in the center of mass, but fitting is sometimes much faster if the system is relative to, e.g., planetary caustic or a cusp. Also think how this should be done for triple lenses. 
1. How to handle full reparametrization of the model? Cassan 2008 is particular example. 


## Use cases to be written 

* Anything from "high level decisions" above.
* Fitting PSPL with free blending and fixed no blending.
* Scaling of observed data to a scale of other dataset. We normally do it to transform follow-up data to survey magnitude scale so that they can be presented on a single plot. 
* Class Model should not allow accesing attributes that shouldn't be there, eg., q for single lens case.


