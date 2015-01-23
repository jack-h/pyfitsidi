import yaml
import tables as tb
import pyfits as pf
import numpy as np

class FitsBuilder(object):
    def __init__(self, configfile, h5file)
        self.configfile = configfile
        self.h5file = h5file
        with open(configfile, 'r') as fh: self.config = yaml.load(fh.read())
        self.h5 = tb.openFile(h5file)
        (self.t_len, self.chan_len, self.bl_len, self.pol_len, self.ri_len) = self.h5.root.xeng_raw0.shape
        self.n_ants = len(self.config['antennas']) 
        self.start_freq = self.h5.attrs.start_freq
        self.chan_bandwidth = self.h5.attrs.chan_bandwidth
        self.primary_hdu = get_primary_table()
        self.n_bands = 1 #hardcoded for now

    def get_primary_table(self):
        hdu = pf.PrimaryHDU()
        for key, val in self.config['PRIMARY'].iteritems():
            hdu.header.set(key, val)

        hdu.verify()
        return hdu

    def get_array_geometry_table(self):
        try:
            ant_order = self.h5.attrs.ant_order
        except KeyError:
            ant_order = range(self.n_ants)

        ants = self.config['array']['antennas']
        n_orb = self.config['array']['NUMORB']
        if n_orb != 0:
            raise NotImplementedError('n_orb != 0 not currently supported!')

        cols = []
        cols.append(pf.Column(name='ANNAME', format='8A',
                              array=['ant%d'%n for n in xrange(self.n_ants)]))
        cols.append(pf.Column(name='STABXYZ', format='3D', unit='METERS',
                              array=[ant['STABXYZ'] for ant in ants]))
        cols.append(pf.Column(name='DERXYZ', format='3E', unit='METERS/SEC',
                              array=[ant['DERXYZ'] for ant in ants]))
        cols.append(pf.Column(name='ORBPARM', format='%dD'%n_orb,
                              array=[]))
        cols.append(pf.Column(name='NOSTA', format='1I',
                              array=ant_order))
        cols.append(pf.Column(name='MNTSTA', format='1J',
                              array=[ant['MNTSTA'] for ant in ants]))
        cols.append(pf.Column(name='STAXOF', format='3E', unit='METERS',
                              array=[ant['STAXOF'] for ant in ants]))
        cols.append(pf.Column(name='DIAMETER', format='1E', unit='METERS',
                              array=[ant['DIAMETER'] for ant in ants]))

        tblhdu = pf.BinTableHDU(from_columns(c))

        for key, val in self.config['ARRAY_GEOMETRY']:
            tblhdu.header.set(key, val)

        self._add_common_headers(tblhdu)

        tblhdu.verify()
        return tblhdu

    def _add_common_headers(self, tbl):
        tbl.header.set('TABREV', 1)
        tbl.header.set('NO_STKD', self.pol_len)
        tbl.header.set('STK_1', 1)
        tbl.header.set('NO_BAND', 1)
        tbl.header.set('NO_CHAN', self.chan_len)
        tbl.header.set('REF_FREQ', self.center_freq)
        tbl.header.set('CHAN_BW', self.chan_bandwidth)
        tbl.header.set('REF_PIXL', 0)

    def get_frequency_table(self):

        cols = []
        cols.append(pf.Column(name='FREQID', format='1J',
                              array=[self.n_bands])
        cols.append(pf.Column(name='BANDFREQ', format='%dD'%self.n_bands,
                              array=[0])
        cols.append(pf.Column(name='CH_WIDTH', format='%dE'%self.n_bands,
                              array=[self.chan_bandwidth])
        cols.append(pf.Column(name='TOTAL_BANDWIDTH', format='%dE'%self.n_bands,
                              array=[self.chan_len * self.chan_bandwidth])
        cols.append(pf.Column(name='SIDEBAND', format='%dE'%self.n_bands,
                              array=[1])

        tblhdu = pf.BinTableHDU(from_columns(c))
        tblhdu.header.set('EXTNAME', 'FREQUENCY')

        self._add_common_headers(tblhdu)

        tblhdu.verify()
        return tblhdu

    def get_source_table(self):
        sources = self.h5.sources

        col_fmt = {
            'SOURCE_ID':
              {'TYPE':'1J', 'UNITS':None},
            'SOURCE':
              {'TYPE':'16A', 'UNITS':None},
            'QUAL':
              {'TYPE':'1J', 'UNITS':None},
            'CALCODE':
              {'TYPE':'4A', 'UNITS':None},
            'FREQID':
              {'TYPE':'1J', 'UNITS':None},
            'IFLUX':
              {'TYPE':'%dE'%self.n_bands, 'UNITS':'Jy'},
            'QFLUX':
              {'TYPE':'%dE'%self.n_bands, 'UNITS':'Jy'},
            'UFLUX':
              {'TYPE':'%dE'%self.n_bands, 'UNITS':'Jy'},
            'VFLUX':
              {'TYPE':'%dE'%self.n_bands, 'UNITS':'Jy'},
            'FREQOFF':
              {'TYPE':'%dE'%self.n_bands, 'UNITS':'Hz'},
            'RAEPO':
              {'TYPE':'1D', 'UNITS':'degrees'},
            'DECEPO':
              {'TYPE':'1D', 'UNITS':'degrees'},
            'EQUINOX':
              {'TYPE':'8A', 'UNITS':None},
            'RAAPP':
              {'TYPE':'1D', 'UNITS':'degrees'},
            'DECAPP':
              {'TYPE':'1D', 'UNITS':'degrees'},
            'SYSVEL':
              {'TYPE':'%dD'%self.n_bands, 'UNITS':'meters/sec'},
            'VELTYP':
              {'TYPE':'8A', 'UNITS':None},
            'VELDEF':
              {'TYPE':'8A', 'UNITS':None},
            'RESTFREQ':
              {'TYPE':'%dD'%self.n_bands, 'UNITS':'Hz'},
            'PMRA':
              {'TYPE':'1D', 'UNITS':'degrees/day'},
            'PMDEC':
              {'TYPE':'1D', 'UNITS':'degrees/day'},
            'PARALLAX':
              {'TYPE':'1E', 'UNITS':'arcseconds'},
            'EPOCH':
              {'TYPE':'1D', 'UNITS':'years'},

        c = []
        for key, col in col_types.iteritems():
            c.append(pf.Column(name=key, format=col['TYPE'], unit=col['UNIT'],
                               array=get(self.h5.sources.attr, [0]))

        tblhdu = pf.BinTableHDU.from_columns(c)
        tblhdu.header.set('EXTNAME', 'SOURCE')

        self._add_common_headers(tblhdu)

        tblhdu.verify()
        return tblhdu

    def get_antenna_table(self):
        
        col_fmt = {
            'TIME':
              {'TYPE':'1D', 'UNITS':'days}',
            'TIME_INTERVAL':
              {'TYPE':'1D', 'UNITS':'days'},
            'ANNAME':
              {'TYPE':'8A', 'UNITS':None},
            'ANTENNA_NO':
              {'TYPE':'1J', 'UNITS':None},
            'ARRAY':
              {'TYPE':'1J', 'UNITS':None},
            'FREQID':
              {'TYPE':'1J', 'UNITS':None},
            'NO_LEVELS':
              {'TYPE':'1J', 'UNITS':None},
            'POLAA':
              {'TYPE':'%dE'self.n_bands, 'UNITS':'degrees'},
            'POLAB':
              {'TYPE':'%dE'self.n_bands, 'UNITS':'degrees'},
            'POLTYA':
              {'TYPE':'1A'},
            'POLTYB':
              {'TYPE':'1A'},
        }

        c = []
        for key, col in col_fmt.iteritems():
            if key in self.config['ANTENNA']['COLUMNS'].keys():
                c.append(pf.Column(name=key, format=col['TYPE'], unit=col['UNIT'],
                                   array=[self.config['ANTENNA']['COLUMNS'][key] for an in xrange(self.n_ants)])

        # antenna specific stuff
        c.append(pf.Column(name='ANNAME', format=col_fmt['ANNAME']['TYPE'],
                           array=[self.config['ARRAY_GEOMETRY']['antennas'][an]['ANNAME'] for an in xrange(self.n_ants)])
        c.append(pf.Column(name='ANTENNA_NO', format=col_fmt['ANTENNA_NO']['TYPE'],
                           array=[an for an in xrange(self.n_ants)])

        
        tblhdu = pf.BinTableHDU.from_columns(c)
        tblhdu.header.set('EXTNAME', 'ANTENNA')

        self._add_common_headers(tblhdu)

        tblhdu.verify()
        return tblhdu
        





        
        
        





        
        
        

