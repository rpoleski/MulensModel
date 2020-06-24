### Windows installation

We are currently working on updating MulensModel installation, so that it will work on either Linux, MacOS, or Windows. In the meantime, please follow the procedure below.

Go to `source/VBBL/`. Then replace the wrapper file:

```
mv VBBinaryLensingLibrary_wrapper_NEW.cpp VBBinaryLensingLibrary_wrapper.cpp
```

Change directory: `cd ../MulensModel/`. Now you have to edit `binarylens.py` file. Remove lines 45-76 and add:

```python
import VBBL
_vbbl_binary_mag_dark = VBBL.VBBinaryLensing_BinaryMagDark
_vbbl_SG12_5 = VBBL.VBBL_SG12_5
```

Save this file and go back: `cd ../../`. Finally, run the compilation:

```
python setup.py install
```

And most probably you still need to add path to MulensModel into your PYTHONPATH.

Much better installation should be available soon.

