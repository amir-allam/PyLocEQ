import sys
import ant_tools
import scec_tools
import core_tools
from antelope.datascope import closing, dbopen

#print 'Reading db'
#with closing(dbopen('/Users/mcwhite/staging/dbs/anza_sub/anza')) as db:
#    tbl_event = db.schema_tables['event']
#    tbl_event = tbl_event.join('origin')
#    tbl_event = tbl_event.subset('time >= _2013319 00:00:00_')
#    tbl_event = tbl_event.separate('event')
#    event_list = ant_tools.create_event_list(tbl_event)
#print 'Done reading db'

event_list = scec_tools.create_event_list('./example_files/pha.phase_210to9_july2011.dat')

locator = core_tools.Locator('pyloceq.cfg')

for ev in event_list:
    origin = ev.preferred_origin
    if origin.lon<-117 or origin.lon>-116 or origin.lat<33.0 or origin.lat>34.0:
        continue
    origin = locator.locate_eq(origin)
    print origin
    sys.exit()
