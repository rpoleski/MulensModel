MulensModel can be installed in a few ways:

### PIP - suggested

The easiest way is to run:
```
pip install MulensModel
```
which will download all files and also install all dependencies. After that you're ready to go!

### PIP in editable mode

If you want to modify the code and see the changes in installed version of MulensModel without need for re-installing, then run:
```
pip install -e .
```

### Using setup.py

One more option is to download the source (e.g., using `git clone https://github.com/rpoleski/MulensModel.git`) and run:
```
python setup.py install
```
MulensModel requires some standard packages plus [astropy package](http://www.astropy.org/). To make sure that you have everything that is needed, just run:
```
pip install -r requirements.txt
```

### Using PYTHONPATH

The last possibility is the most complicated. You should go to `source/VBBL/` and run `make`. Then do the same in `source/AdaptiveContouring/`. Finally, add the full path `/.../MulensModel/source` to your `PYTHONPATH`. And you probably want to do this everytime you login, to add `PYTHONPATH` changes to your `~/.bashrc`, `.cshrc`, `.login` etc. You also have to install all required packages (see above).

---
If you still have problems, then please contact code authors.
