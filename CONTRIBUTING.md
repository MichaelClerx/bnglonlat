# Contributing to bnglonlat

This is not encouraged!
It only exists because some users were having issues with `lonlat_bng`.
Hopefully these won't exist for future versions and this package can be dropped.

## Installing

Install for development with:

```
pip install -e .[dev]
```

## Testing

Make sure you've installed `lonlat_bng` and run

```
python -m unittest
```

## Uploading to PyPi

```
python setup.py sdist bdist_wheel
twine upload dist/*
rm build/ dist/ -rf
```

Check:

    https://pypi.org/project/bnglonlat/#history


