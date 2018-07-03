### How to install?

1. Make sure you have python with [astropy package](http://www.astropy.org/) installed.
2. Download source code - either [recent release](https://github.com/rpoleski/MulensModel/releases) or current repository using green button above.
3. Unpack the archive.
4. Add the path to the unpack directory to the PYTHONPATH, e.g., if you've extracted the archive in your home directory (``/home/USER_NAME/``) in tcsh:
```
setenv PYTHONPATH /home/USER_NAME/MulensModel-1.3.0/source\:$PYTHONPATH
```
in bash:
```
export PYTHONPATH=/home/USER_NAME/MulensModel-1.3.0/source:$PYTHONPATH
```
In order to have this command invoked every time you open a terminal, please add this command to your startup file (``~/.cshrc``, ``~/.bashrc``, ``~/.profile`` or similar). If you didn't have ``PYTHONPATH`` defined before, then skip the last part of the above commands.

5. Go to subdirectory ```source/VBBL/``` and run ```make``` command. If it's not working and you're using Windows, then please run:
```
gcc -lm -lstdc++ -fPIC -c VBBinaryLensingLibrary.cpp
gcc -Wl,-soname,rapper -shared -o VBBinaryLensingLibrary_wrapper.so VBBinaryLensingLibrary_wrapper.cpp -lm -lstdc++ -fPIC VBBinaryLensingLibrary.o
```
6. Repeat above in ```source/AdaptiveContouring/```
7. Run ```py.test``` in ```source/MulensModel``` to check that all unit tests pass.
8. Congratulations! You have MulensModel installed fully.

