def create_event_list(inp, fmt):
    """
    Create and return a list of Event objects.

    Arguments:
    inp - A (potentially subsetted) View of the Event table from the
    CSS3.0 database schema. Alternatively the path to a SCEDC format
    flat file.

    fmt - The format of the 'inp' argument; 'CSS3.0' or 'SCEDC'.

    Return Values:
    A list of Event objects.
    """
    from eqloc3d_classes import Event
    import time as pytime
    import calendar
    event_list = []
    if fmt == 'CSS3.0':
        for record1 in inp.iter_record():
            evid, evname, prefor, auth, commid, lddate = record1.getv('evid',
                                                                     'evname',
                                                                     'prefor',
                                                                     'auth',
                                                                     'commid',
                                                                     'lddate')
            event = Event(evid,
                          prefor,
                          evname=evname,
                          auth=auth,
                          commid=commid,
                          lddate=lddate)
            view2 = inp.subset('evid == %d' % evid)
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
                arrivals = [Phase(sta, time, phase, arid=arid)
                            for sta, time, phase, arid in arrival_data]
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
    elif fmt == 'SCEDC':
        infile = open(inp, 'r')
        event = None
        for line in infile:
            line = line.strip().split()
            #In SCEDC format, event lines begin with '#'
            if line[0] == '#':
                if event != None:
                    event_list += [event]
                year = int(line[1])
                month = int(line[2])
                day = int(line[3])
                hour = int(line[4])
                minute = int(line[5])
                second = float(line[6])
                lat = float(line[7])
                lon = float(line[8])
                depth = float(line[9])
                mag = float(line[10])
                orid = int(line[14])
                time = calendar.timegm(pytime.struct_time(list([year,
                                              month,
                                              day,
                                              hour,
                                              minute,
                                              second,
                                              0, 0, 0])))
                event = Event(evid_ctr,
                              orid,
                              auth='SCEDC',
                              lddate=pytime.time())
                event.add_origin(lat, lon, depth, time, orid, 'SCEDC')
                event.set_preferred_origin(orid)
                evid_ctr += 1
            else:
                sta = line[0]
                arrtime = float(line[1]) + time
                qual = float(line[2])
                iphase = line[3]
                event.preferred_origin.arrivals += [Phase(sta,
                                                          arrtime,
                                                          iphase)]

    else:
        raise Exception('Input format %s not recognized.' % fmt)
    return event_list

