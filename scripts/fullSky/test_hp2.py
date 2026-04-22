import healpy as hp

hp_12 = 122654106
hp_5 = hp_12 >> 14

print(f"Level 5 pixel: {hp_5}")

nside_5 = 32
theta, phi = hp.pix2ang(nside_5, hp_5, nest=True)
import numpy as np
ra = np.degrees(phi)
dec = 90.0 - np.degrees(theta)
print(f"Level 5 RA={ra}, Dec={dec}")
