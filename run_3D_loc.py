import misc_tools
from mtools import *
from numpy import arange,asarray
from antelope.stock import pfread, pfin
from antelope.datascope import closing, dbopen
params=pfin('eqloc3d.pf')
loc_params = params['location_parameters']
#These should go in parameter file.
#nr = int(loc_params['nr'])
#nlat = int(loc_params['nlat'])
#nlon = int(loc_params['nlon'])
#nx, ny, nz = nlon, nlat, nr
earth_rad=6371

#Load events
print 'Reading db'
with closing(dbopen('/Users/mcwhite/staging/dbs/anza_sub/anza')) as db:
    tbl_event = db.schema_tables['event']
    tbl_event = tbl_event.join('origin')
    tbl_event = tbl_event.subset('time >= _2013319 00:00:00_')
    tbl_event = tbl_event.separate('event')
    event_list = misc_tools.create_event_list(tbl_event, 'CSS3.0')
print 'Done reading db'

for ev in event_list:
    origin = ev.preferred_origin
    if origin.lon<-117 or origin.lon>-116 or origin.lat<33.0 or origin.lat>34.0:
        continue
    misc_tools.locate_eq(origin)
    sys.exit()
