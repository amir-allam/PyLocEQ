import sys

def _main():
    import ant_tools
    from core_tools import Locator
    from antelope.datascope import closing, dbopen
    args = _parse_command_line()
    if args.pf: ant_tools.pf_2_cfg(args.pf, 'pyloceq')
    else: ant_tools.pf_2_cfg('pyloceq', 'pyloceq')
    locator = Locator('pyloceq.cfg')
    with closing(dbopen(args.db, 'r+')) as db:
        tbl_event = db.schema_tables['event']
        if args.subset:
            view = tbl_event.join('origin')
            view = view.subset(args.subset)
            tbl_event = view.separate('event')
        event_list = ant_tools.create_event_list(tbl_event)
        for event in event_list:
            origin = event.preferred_origin
            origin = locator.locate_eq(origin)
            if origin == None:
                continue
            ant_tools.write_origin(origin, db)
    return 0

def _parse_command_line():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('db', type=str, help='Input/output databse.')
    parser.add_argument('-s', '--subset', type=str, help='Subset expression.')
    parser.add_argument('-p', '--pf', type=str, help='Parameter file.')
    return parser.parse_args()

if __name__ == '__main__': sys.exit(_main())
else: raise ImportError
