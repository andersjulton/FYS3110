import numpy as np
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'



with open('Initial_values.txt', 'w') as outfile:
    names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
    outfile.write(str(len(names))+"\n")
    outfile.write("Mercury Venus Earth Mars Jupiter Saturn Uranus Neptune Pluto\n")
    for i, n in enumerate(names):
        id = str(i+1)+'99'
        obj = Horizons(id=id, id_type= 'majorbody', epochs={'start':'2018-10-16', 'stop':'2018-11-16','step':'1d'})
        vec = obj.vectors()
        outfile.write("%e   %e   %e   %e   %e   %e" % (vec['x'][0], vec['y'][0], vec['z'][0], vec['vx'][0], vec['vy'][0], vec['vz'][0]))
        outfile.write("\n")
