# openEMS python interface

## Install
If openEMS was installed into `~/opt`, then install this package with:

```bash
python setup.py build_ext -I ~/opt/include -L ~/opt/lib -R ~/opt/lib
python setup.py install
```

Otherwise, replace `~/opt` with the path to the place where it was installed.
