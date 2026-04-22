import healpy as hp

hp_12 = 122654106
nside_12 = 4096

theta, phi = hp.pix2ang(nside_12, hp_12, nest=True)
import numpy as np
ra = np.degrees(phi)
dec = 90.0 - np.degrees(theta)
print(f"RA={ra}, Dec={dec}")
