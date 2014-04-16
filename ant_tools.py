import sys
import os
sys.path.append('%s/data/python' % os.environ['ANTELOPE'])
from antelope.datascope import DbfindEnd,\
                               dbTABLE_FIELDS,\
                               dbTABLE_NAME
from antelope.stock import pfin
import antpy

def create_event_list(view):
    """
    Create and return a list of core_tools.Event objects.

    Arguments:
    view - A datascope database pointer to a (potentially subsetted) View
    of the Event table from the CSS3.0 database schema.

    Return Values:
    A list of Event objects.

    Behaviour:
    This method does NOT open or close the database passed in.
    """
    from temp_core_tools import Event, Phase
    import time as pytime
    import calendar
    event_list = []
    for record1 in view.iter_record():
        evid, evname, prefor, auth, commid, lddate = record1.getv('evid',
                                                                 'evname',
                                                                 'prefor',
                                                                 'auth',
                                                                 'commid',
                                                                 'lddate')
        event = Event(prefor,
                      evid=evid,
                      evname=evname,
                      auth=auth,
                      commid=commid,
                      lddate=lddate)
        view2 = view.subset('evid == %d' % evid)
        view2 = view2.join('origin')
        view2 = view2.separate('origin')
        for record2 in view2.iter_record():
            lat = record2.getv('lat')[0]
            lon = record2.getv('lon')[0]
            depth = record2.getv('depth')[0]
            time = record2.getv('time')[0]
            orid = record2.getv('orid')[0]
            evid = record2.getv('evid')[0]
            jdate = record2.getv('jdate')[0]
            nass = record2.getv('nass')[0]
            ndef = record2.getv('ndef')[0]
            ndp = record2.getv('ndp')[0]
            grn = record2.getv('grn')[0]
            srn = record2.getv('srn')[0]
            etype = record2.getv('etype')[0]
            review = record2.getv('review')[0]
            depdp = record2.getv('depdp')[0]
            dtype = record2.getv('dtype')[0]
            mb = record2.getv('mb')[0]
            mbid = record2.getv('mbid')[0]
            ms = record2.getv('ms')[0]
            msid = record2.getv('msid')[0]
            ml = record2.getv('ml')[0]
            mlid = record2.getv('mlid')[0]
            algorithm = record2.getv('algorithm')[0]
            auth = record2.getv('auth')[0]
            commid = record2.getv('commid')[0]
            lddate = record2.getv('lddate')[0]
            view3 = view2.subset('orid == %d' % orid)
            view3 = view3.join('assoc')
            view3 = view3.join('arrival')
            arrival_data = [record3.getv('sta',
                                         'arrival.time',
                                         'iphase', 'arid')\
                                         for record3 in view3.iter_record()]
            arrivals = [Phase(name, time, phase, arid=arid)
                        for name, time, phase, arid in arrival_data]
            event.add_origin(lat,
                             lon,
                             depth,
                             time,
                             auth,
                             arrivals=arrivals,
                             orid=orid,
                             evid=evid,
                             jdate=jdate,
                             nass=nass,
                             ndef=ndef,
                             ndp=ndp,
                             grn=grn,
                             srn=srn,
                             etype=etype,
                             review=review,
                             depdp=depdp,
                             dtype=dtype,
                             mb=mb,
                             mbid=mbid,
                             ms=ms,
                             msid=msid,
                             ml=ml,
                             mlid=mlid,
                             algorithm=algorithm,
                             commid=commid,
                             lddate=lddate)
        event.set_preferred_origin(event.prefor)
        event_list += [event]
    return event_list

def create_station_list(view):
    """
    Create and return a list of core_tools.Station objects.

    Arguments:
    view - A datascope database pointer to a (potentially subsetted) View
    of the Site table from the CSS3.0 database schema.

    Return Values:
    A list of Station objects.

    Behaviour:
    This method does NOT open or close the database passed in.
    """
    from core_tools import Station
    view = view.sort('sta', unique=True)
    station_list = []
    for record in view.iter_record():
        name, lat, lon, elev = record.getv('sta', 'lat', 'lon', 'elev')
        station_list += [Station(name, lat, lon, elev)]
    return station_list

def write_origin(origin, output):
    """
    Write an origin in a given output format.

    Arguments:
    origin - A core_tools.Origin object to be written out.
    output - A datascope database pointer to an open output database.

    Return Values:
    0 - Sucess
    -1 - Failure

    Behaviour:
    This method does NOT open or close the database passed in.

    Additional Comments:
    This method will assumes that the database being written out is the
    same as the input databse. Ie. NO arrival rows are created, they are
    assumed to already exist.
    """
    from time import time
    tbl_origin = output.schema_tables['origin']
    origin.orid = output.nextid('orid')
    origin = map_null_values(tbl_origin, origin)
    tbl_origin.record = tbl_origin.addnull()
    tbl_origin.putv(('lat', origin.lat),
                    ('lon', origin.lon),
                    ('depth', origin.depth),
                    ('time', origin.time),
                    ('orid', origin.orid),
                    ('evid', origin.evid),
                    ('auth', origin.auth),
                    ('jdate', origin.jdate),
                    ('nass', origin.nass),
                    ('ndef', origin.ndef),
                    ('ndp', origin.ndp),
                    ('grn', origin.grn),
                    ('srn', origin.srn),
                    ('etype', origin.etype),
                    ('review', origin.review),
                    ('depdp', origin.depdp),
                    ('dtype', origin.dtype),
                    ('mb', origin.mb),
                    ('mbid', origin.mbid),
                    ('ms', origin.ms),
                    ('msid', origin.msid),
                    ('ml', origin.ml),
                    ('mlid', origin.mlid),
                    ('algorithm', origin.algorithm),
                    ('commid', origin.commid))
                    #'lddate', origin.lddate)
    tbl_event = output.schema_tables['event']
    tbl_event.record = tbl_event.find('evid == %d' % origin.evid)
    tbl_event.putv(('prefor', origin.orid))
    tbl_assoc = output.schema_tables['assoc']
    tbl_predarr = output.schema_tables['predarr']
    tbl_site = output.schema_tables['site']
    for arrival in origin.arrivals:
#THIS IS A HACK!
#predicted arrival time should be calculated for EVERY arrival
        if arrival.tt_calc == None:
            continue
        tbl_assoc.record = tbl_assoc.addnull()
        tbl_assoc.putv(('arid', arrival.arid),
                       ('orid', origin.orid),
                       ('sta', arrival.sta),
                       ('phase', arrival.phase))
        view = tbl_site.subset('sta =~ /%s/ && ondate < _%f_ && (offdate == -1 || offdate > _%f_)' % (arrival.sta, time(), time()))
        view.record = 0
        stalat, stalon = view.getv('lat', 'lon')
        tbl_predarr.record = tbl_predarr.addnull()
        tbl_predarr.putv(('arid', arrival.arid),
                         ('orid', origin.orid),
                         ('time', origin.time + arrival.tt_calc),
                         ('slow', antpy.distance(stalat, stalon,
                                                 origin.lat, origin.lon)
                                                 / arrival.tt_calc),
                         ('seaz', antpy.azimuth(stalat,
                                                stalon,
                                                origin.lat,
                                                origin.lon)),
                         ('esaz', antpy.azimuth(origin.lat,
                                                origin.lon,
                                                stalat,
                                                stalon)))

    return 0

def map_null_values(table, obj):
    """
    Maps 'None' field values to appropriate CSS3.0 null field values.
    """
    from antpy import get_null_value
    for field in table.query(dbTABLE_FIELDS):
        if getattr(obj, field) == None:
           setattr(obj, field, get_null_value(table.query(dbTABLE_NAME), field))
    return obj

def pf_2_cfg(pf, config_file):
    import ConfigParser
    config_file = '%s.cfg' % config_file
    pf = '%s.pf' % pf
    if os.path.isfile(config_file):
        try:
            os.remove(config_file)
        except Exception:
            print 'Could not remove potentially stale configuration file - %s.\n'\
                    'Please remove and try again.' % config_file
            sys.exit(-1)
    config = ConfigParser.RawConfigParser()
    config.add_section('misc')
    pf = pfin('%s' % pf)
    for key1 in pf.keys():
        if isinstance(pf[key1], dict):
            config.add_section(key1)
            for key2 in pf[key1]:
                config.set(key1, key2, pf[key1][key2])
        else:
            config.set('misc', key1, pf[key1])
    with open(config_file, 'w') as config_file:
        config.write(config_file)
