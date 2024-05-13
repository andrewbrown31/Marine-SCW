import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '10m_wind_gust_since_previous_post_processing',
            '2m_dewpoint_temperature', '2m_temperature', 'convective_available_potential_energy',
            'convective_precipitation', 'geopotential', 'orography', 'surface_pressure',
            'total_precipitation',
        ],
        'year': '2023',
        'month': '10',
        'day': [
            '15','16',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            -9, 111, -45,
            160,
        ],
    },
    '/g/data/w40/ab4502/IN2023_V06/data/era5/era5_sfc.nc')

