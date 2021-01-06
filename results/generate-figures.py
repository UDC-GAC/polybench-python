import os

for x in [1,4,5,7,9,10,11,12]:
    for v in ["paper","regenerated"]:
        os.system( "python fig%d.py %s" % (x,v) )
