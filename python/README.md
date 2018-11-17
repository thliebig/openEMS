# openEMS python interface

## Install
* Simple version:
```python
python setup.py install
```

* Extended options, e.g. for custom install path at */opt/openEMS*:
```python
python setup.py build_ext -I/opt/openEMS/include -L/opt/openEMS/lib -R/opt/openEMS/lib"
python setup.py install
```
**Note:** The install command may require root on Linux, or add --user to install to ~/.local
