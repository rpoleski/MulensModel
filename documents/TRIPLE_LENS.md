## Triple lens

**First decide on frame used** - e.g. check:
* [ob07349](https://ui.adsabs.harvard.edu/#abs/2016AJ....152..125B/abstract) and make sure t\_0 and u\_0 can be predicted in center-of-mass frame.
* some other triple lens event - can we get 3L parameters from 2L solutions?

Jen also suggests primary at origin and axis pointing to second object.

My main idea - center of mass as origin and first object on X axis.

# To do:

* "reset" triple\_lens branch
* Binary lens positions - calculation and plotting (UC35)
* Use cases - we have 07 already, which does basic LC plotting
* Utils.\_triple\_lens\_positions\_to\_parameters() - easy
* Utils.\_parameters\_to\_triple\_lens\_positions() - already coded
* Utils: - (q\_21, q\_31) <=> (eps\_1, eps\_2, eps\_3) - easy
* unit tests
* Caustics - use the same class for triple lenses?
* MagnificationCurve.get\_triple\_lens\_magnification (include hooks to hexa, quad and PSPL)
* ModelParameters.\_\_repr\_\_ - finish
* ModelParameters.\_check\_valid\_combination\_1\_source - this needs some thinking
* ModelParameters - all the new parameters starting at s\_21
* Trajectory.get\_xy
* TripleLens.point\_source\_magnification
* TripleLens.hexadecapole\_magnification

