#class StationList(list):
#    """
#    A container class for a list of Station objects.
#
#    This class replaces the deprecated Stalist class.
#    """
#    def __init__(self, sta, lat, lon, elev):
#        self.append(Station(sta, lat, lon, elev)
#
#    def __str__(self):
#        ret = 'StationList Object\n------------------\n'
#        for sta in self:
#             ret += '%s\n' % sta
##            ret += 'sta:\t\t%s\n' % sta.name
##            ret += 'lat:\t\t%s\n' % sta.lat
##            ret += 'lon:\t\t%s\n' % sta.lon
##            ret += 'elev:\t\t%s\n' % sta.elev
#        return ret
#
#    def _init_db(self, db):
#        """
#        Initialize station list using a CSS3.0 database as input.
#        """
#        with closing(dbopen(db, 'r')) as db:
#            tbl_site = db.schema_tables['site']
#            tbl_site = tbl_site.sort('sta', unique=True)
#            for record in tbl_site.iter_record():
#                sta, lat, lon, elev = record.getv('sta', 'lat', 'lon', 'elev')
#                self.append(Station(sta, lat, lon, elev))
#
#    def _init_scedc(self, infile):
#        """
#        Initialize station list using SCEDC format flat file as input.
#        """
#        infile = open(infile, 'r')
#        for line in infile:
#            line = line.strip().split() #stripping may be uneccessary
#            self.append(Station(line[0],
#                                 float(line[1]),
#                                 float(line[2]),
#                                 float(line[3])))

class Station:
    """
    A containter class for station location data.
    """
    def __init__(self, name, lat, lon, elev):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.elev = elev

    def __str__(self):
        ret = 'Station Object\n--------------\n'
        ret += 'name:\t\t%s' % name
        ret += 'lat:\t\t%s' % lat
        ret += 'lon:\t\t%s' % lon
        ret += 'elev:\t\t%s' % elev
        return ret

class Event():
    """
    A container class for earthquake event metadata.
    """
    #def __init__(self, time, lat, lon, depth, mag, magtype=None, evid=None):
    def __init__(self,
                 prefor,
                 evid=None,
                 evname=None,
                 auth=None,
                 commid=None,
                 lddate=None,
                 origins=None):
        """
        Initialize Event object using one of two possible inputs.
        """
        import time as pytime
        self.evid = evid
        self.evname = evname
        self.prefor = prefor
        self.auth = auth
        self.commid = commid
        self.lddate = lddate
        self.preferred_origin = None
        if origins == None: self.origins = []
        else: self.origins = origins

    def __str__(self):
        ret = 'Event Object\n------------\n'
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'evname:\t\t%s\n' % self.evname
        ret += 'prefor:\t\t%s\n' % self.prefor
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'origins:\n'
        if len(self.origins) == 0:
            ret += '\t\tNone\n'
        else:
            for i in range(len(self.origins)):
                for line in  ('%s' % self.origins[i]).split('\n'):
                    ret += '\t\t%s\n' % line
        return ret

    def set_preferred_origin(self, prefor):
        """
        Set self.preferred_origin to equal origin with orid == prefor.
        """
        for i in range(len(self.origins)):
            if self.origins[i].orid == prefor:
                self.preferred_origin = self.origins[i]
                return 0

    def add_origin(self,
                   lat,
                   lon,
                   depth,
                   time,
                   auth,
                   arrivals=[],
                   orid=None,
                   evid=None,
                   jdate=None,
                   nass=None,
                   ndef=None,
                   ndp=None,
                   grn=None,
                   srn=None,
                   etype=None,
                   review=None,
                   depdp=None,
                   dtype=None,
                   mb=None,
                   mbid=None,
                   ms=None,
                   msid=None,
                   ml=None,
                   mlid=None,
                   algorithm=None,
                   commid=None,
                   lddate=None):
        """
        Add an Origin object to the list of origins associated with this event.
        """
        self.origins += [Origin(lat,
                                lon,
                                depth,
                                time,
                                auth,
                                orid=orid,
                                evid=evid,
                                arrivals=arrivals,
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
                                lddate=lddate)]
    def write(self, out, fmt):
        """
        Write out newly authored data pertaining to event.

        Arguments:
        out - A datascope db pointer to an open CSS3.0 database for output.
        Alternatively the path to an output SCEDC format flat file.

        fmt - The format of the 'out' argument; 'CSS3.0' or 'SCEDC'.

        Return Values:
        0 - Success
        -1 - Failure
        """
        if fmt == 'CSS3.0':
            return self._write_CSS()
        elif fmt == 'SCEDC':
            #Do it the SCEDC way
            pass
        else:
            raise Exception('Output format %s not recognized' % fmt)

class Origin():
    """
    A container class for origin data.
    """
    def __init__(self,
                 lat,
                 lon,
                 depth,
                 time,
                 auth,
                 arrivals=[],
                 orid=None,
                 evid=None,
                 jdate=None,
                 nass=None,
                 ndef=None,
                 ndp=None,
                 grn=None,
                 srn=None,
                 etype=None,
                 review=None,
                 depdp=None,
                 dtype=None,
                 mb=None,
                 mbid=None,
                 ms=None,
                 msid=None,
                 ml=None,
                 mlid=None,
                 algorithm=None,
                 commid=None,
                 lddate=None):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.time = time
        self.orid = orid
        self.evid = evid
        self.auth = auth
        self.arrivals = arrivals
        self.jdate = jdate
        self.nass = nass
        self.ndef = ndef
        self.ndp = ndp
        self.grn = grn
        self.srn = srn
        self.etype = etype
        self.review = review
        self.depdp = depdp
        self.dtype = dtype
        self.mb = mb
        self.mbid = mbid
        self.ms = ms
        self.msid = msid
        self.ml = ml
        self.mlid = mlid
        self.algorithm = algorithm
        self.commid = commid
        self.lddate = lddate

    def __str__(self):
        """
        Return string representation of Origin object.
        """
        ret = 'Origin Object\n-------------\n'
        ret += 'lat:\t\t%s\n' % self.lat
        ret += 'lon:\t\t%s\n' % self.lon
        ret += 'depth:\t\t%s\n' % self.depth
        ret += 'time:\t\t%s\n' % self.time
        ret += 'orid:\t\t%s\n' % self.orid
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'jdate:\t\t%s\n' % self.jdate
        ret += 'nass:\t\t%s\n' % self.nass
        ret += 'ndef:\t\t%s\n' % self.ndef
        ret += 'ndp:\t\t%s\n' % self.ndp
        ret += 'grn:\t\t%s\n' % self.grn
        ret += 'srn:\t\t%s\n' % self.srn
        ret += 'etype:\t\t%s\n' % self.etype
        ret += 'review:\t\t%s\n' % self.review
        ret += 'depdp:\t\t%s\n' % self.depdp
        ret += 'dtype:\t\t%s\n' % self.dtype
        ret += 'mb:\t\t%s\n' % self.mb
        ret += 'mbid:\t\t%s\n' % self.mbid
        ret += 'ms:\t\t%s\n' % self.ms
        ret += 'msid:\t\t%s\n' % self.msid
        ret += 'ml:\t\t%s\n' % self.ml
        ret += 'mlid:\t\t%s\n' % self.mlid
        ret += 'algorithm:\t\t%s\n' % self.algorithm
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'arrivals:\n'
        for i in range(len(self.arrivals)):
            ret += '\t\t%s' % self.arrivals[i]
        return ret

class Phase():
    """
    A container class for phase data.
    """
    def __init__(self, sta, time, phase, chan=None, qual=None, arid=None):
        self.sta = sta
        self.time = time
        self.phase = phase
        self.chan = chan
        self.qual = qual
        self.arid = arid

    def __str__(self):
        ret = 'Arrival Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'time:\t\t%s\n' % self.time
        ret += 'phase:\t\t%s\n' % self.phase
        ret += 'qual:\t\t%s\n'  % self.qual
        return ret

