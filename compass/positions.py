"""
Reference positions

"""

from astropy.coordinates import SkyCoord

reference_position = {
    "b335": SkyCoord("19h37m00.901s", "+07d34m09.44s", frame="icrs"),
    "bhr71": SkyCoord("12h01m36.476s", "-65d08m49.37s", frame="icrs"),  # Continuum peak
    "hops108": SkyCoord("05h35m27.084s", "-05d10m00.07s", frame="icrs"),
    "hops373": SkyCoord("05h46m30.905s", "-00d02m35.20s", frame="icrs"),
    "iras4a": SkyCoord("03h29m10.438s", "+31d13m32.05s", frame="icrs"),
    "l1551irs5": SkyCoord("04h31m34.164s", "+18d08m04.72s", frame="icrs"),
    "svs13": SkyCoord("03h29m03.748s", "+31d16m03.73s", frame="icrs"),
    "v883ori": SkyCoord("05h38m18.100s", "-07d02m26.00s", frame="icrs"),
    "serps18": SkyCoord("18h30m04.118s", "-02d03m02.55s", frame="icrs"),
    "serp11": SkyCoord("18h29m06.62s", "00d30m33.82s", frame="icrs"),
    "serp1a": SkyCoord("18h29m49.793s", "01d15m20.20s", frame="icrs"),
}
