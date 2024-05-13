import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-complete',
    {
        'format': 'netcdf',
        'param': '129/130/131/132/157',
        'levtype': 'ml',
        'levelist': '1/to/137',
        'stream': 'oper',
        'type': 'an',
        'date': '2023-10-15/to/2023-10-16',
        'time': '00/to/23',
        'area': '-9/111/-45/160',
        'grid': '0.25/0.25',
    },
    '/g/data/w40/ab4502/IN2023_V06/data/era5/era5_ml.nc')

