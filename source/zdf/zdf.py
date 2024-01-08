"""
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
"""

#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

class ZDF_Record:
    def __init__( self ):
        self.pos = 0
        self.id = 0
        self.name = ""
        self.len = 0
    def version(self):
        return( self.id & 0x0000FFFF )

    def typeid(self):
        return( self.id & 0xFFFF0000 )

class ZDF_Iteration:
    def __init__( self ):
        self.name = ""
        self.n = 0
        self.t = 0.
        self.tunits = ""

class ZDF_Grid_Axis:
    def __init__( self ):
        self.name = ''
        self.type = 0
        self.min = 0.
        self.max = 0.
        self.label = ''
        self.units = ''

class ZDF_Grid_Info:
    def __init__( self ):
        self.name = ''
        self.ndims = 0
        self.nx = []
        self.label = ''
        self.units = ''
        self.has_axis = 0
        self.axis = []

class ZDF_Part_Info:
    def __init__( self ):
        self.name = ''
        self.label = ''
        self.nquants = 0
        self.quants = []
        self.qlabels = dict()
        self.qunits = dict()
        self.nparts = 0

class ZDF_Tracks_Info:
    def __init__( self ):
        self.name = ''
        self.label = ''
        self.ntracks = 0
        self.ndump = 0
        self.niter = 0
        self.nquants = 0
        self.quants = []
        self.qlabels = []
        self.qunits = []
        
class ZDFfile:
    def __init__(self, file_name):
        self.__file = open(file_name, "rb")

        # Check magic number
        magic = self.__file.read(4)
        if (magic != b'ZDF1'):
            print('File is not a proper ZDF file, aborting', file=sys.stderr)
            self.__file.close

    def close(self):
        self.__file.close

    def __read_uint32(self):
        data = np.fromfile(self.__file,dtype='<u4',count=1)
        if ( data.size == 0 ):
            return False
        else:
            return data[0]

    def __read_int32(self):
        return np.fromfile(self.__file,dtype='<i4',count=1)[0]

    def __read_uint64(self):
        return np.fromfile(self.__file,dtype='<u8',count=1)[0]

    def __read_int64(self):
        return np.fromfile(self.__file,dtype='<i8',count=1)[0]

    def __read_float32(self):
        return np.fromfile(self.__file,dtype='<f4',count=1)[0]

    def __read_float64(self):
        return np.fromfile(self.__file,dtype='<f8',count=1)[0]

    def __read_string(self):

        length = self.__read_uint32()

        if (length > 0):
            data = self.__file.read(length)
            fstring = data.decode()

            # Data is always written in blocks of 4 byt
            pad = ((length - 1) // 4 + 1) * 4 - length
            if ( pad > 0 ):
                self.__file.seek(pad,1)
        else:
            fstring = ""

        return fstring

    def record_type(self, typeTag):
        version = typeTag & 0x0000FFFF
        typeID = typeTag & 0xFFFF0000

        types = {0x00010000: "int",
                 0x00020000: "double",
                 0x00030000: "string",
                 0x00100000: "dataset",
                 0x00110000: "cdset_start",
                 0x00120000: "cdset_chunk",
                 0x00130000: "cdset_end",
                 0x00200000: "iteration",
                 0x00210000: "grid_info",
                 0x00220000: "part_info",
                 0x00230000: "track_info", }

        if (typeID in types):
            return types[typeID]
        else:
            return "unknown"

# -----------------------------------------------------------------------------
# Array datatypes
# -----------------------------------------------------------------------------
    
    def __read_int32_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<i4',count=size)
        data.shape = np.flip(nx)
        return data

    def __read_uint32_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<u4',count=size)
        data.shape = np.flip(nx)
        return data

    def __read_int64_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<i8',count=size)
        data.shape = np.flip(nx)
        return data

    def __read_uint64_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<u8',count=size)
        data.shape = np.flip(nx)
        return data

    def __read_float32_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<f4',count=size)
        data.shape = np.flip(nx)
        return data

    def __read_float64_arr(self, nx):
        size = np.prod(nx)
        data = np.fromfile(self.__file,dtype='<f8',count=size)
        data.shape = np.flip(nx)
        return data
    
    def __read_arr( self, dtype, nx ):
        
        if ( dtype == 5 ):
            data = self.__read_int32_arr(nx)           
        elif ( dtype == 3 ):
            data = self.__read_uint32_arr(nx)           
        elif ( dtype == 7 ):
            data = self.__read_int64_arr(nx)           
        elif ( dtype == 8 ):
            data = self.__read_uint64_arr(nx)           
        elif ( dtype == 9 ):
            data = self.__read_float32_arr(nx)           
        elif ( dtype == 10 ):
            data = self.__read_float64_arr(nx)           
        else:
            print( '(*error*) ZDF: Data type not yet supported.' , file=sys.stderr)
            data = False
        
        return data
        
# -----------------------------------------------------------------------------
# Low level interfaces
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Read record header
# -----------------------------------------------------------------------------

    def read_record(self, skip=False):

        rec = ZDF_Record()
        rec.pos = self.__file.tell()

        # Read record id and check for EOF
        rec.id = self.__read_uint32()

        if (rec.id is False):
            # If end of file return false
            return False

        # Read name and length
        rec.name = self.__read_string()
        rec.len  = self.__read_uint64()

        if (skip):
            self.__file.seek(rec.len, 1)

        return rec
    
    def __record_skip( self, rec ):
        self.__file.seek(rec.len, 1)

# -----------------------------------------------------------------------------
# Read string
# -----------------------------------------------------------------------------

    def read_string(self, rec = False):
        if ( rec is False ):
            rec = self.read_record()
        fstring = self.__read_string()
        return fstring

# -----------------------------------------------------------------------------
# Read iteration
# -----------------------------------------------------------------------------

    def read_iteration(self, rec = False):
        if ( rec is False ):
            rec = self.read_record()
        iteration = ZDF_Iteration()
        iteration.name = rec.name
        iteration.n = self.__read_int32()
        iteration.t = self.__read_float64()
        iteration.tunits = self.__read_string()
        return iteration

# -----------------------------------------------------------------------------
# Read grid info
# -----------------------------------------------------------------------------

    def read_grid_info(self, rec = False):
        if ( rec is False ):
            rec = self.read_record()

        # Maximum supported version
        max_version = 0x00000001

        # Get version
        version = rec.version()
        if ( version > max_version ):
            print( '(*error*) ZDF: Grid info version is higher than supported.' , file=sys.stderr)
            print( '(*error*) ZDF: Please update the code to a newer version.' , file=sys.stderr)
            return False

        info = ZDF_Grid_Info()
        
        info.name = rec.name
        info.ndims = self.__read_uint32()
        info.nx = self.__read_uint64_arr( info.ndims )

        info.label = self.__read_string()
        info.units = self.__read_string()
        info.has_axis = self.__read_int32()

        if ( info.has_axis ):
            for i in range(info.ndims):
                ax = ZDF_Grid_Axis()
                if (version > 0):
                    ax.name = self.__read_string()
                else:
                    ax.name = 'axis_{}'.format(i)
                ax.type  = self.__read_int32()
                ax.min   = self.__read_float64()
                ax.max   = self.__read_float64()
                ax.label = self.__read_string()
                ax.units = self.__read_string()
                info.axis.append(ax)

        return info

# -----------------------------------------------------------------------------
# Read particle info
# -----------------------------------------------------------------------------
    
    def read_part_info(self, rec = False):
        if ( rec is False ):
            rec = self.read_record()

        # Maximum supported version
        max_version = 0x00000002

        # Get version
        version = rec.version()
        if ( version > max_version ):
            print( '(*error*) ZDF: Particles info version is higher than supported.' , file=sys.stderr)
            print( '(*error*) ZDF: Please update the code to a newer version.' , file=sys.stderr)
            return False
        
        info      = ZDF_Part_Info()
        info.name = rec.name
        info.label = self.__read_string()
        
        if ( version >= 1 ):
            # version 1
            info.nparts = self.__read_uint64()
            info.nquants = self.__read_uint32()
            
            for i in range(info.nquants):
                info.quants.append( self.__read_string() )
            for q in info.quants:
                info.qlabels[q] = self.__read_string()            
            for q in info.quants:
                info.qunits[q] = self.__read_string()
            
        else:
            # version 0
            info.nquants = self.__read_uint32()

            for i in range(info.nquants):
                info.quants.append( self.__read_string() )

            # version 0 does not have label information
            for q in info.quants:
                info.qlabels[q] = q

            for q in info.quants:
                info.qunits[q] = self.__read_string()

            info.nparts = self.__read_uint64()

        return info

# -----------------------------------------------------------------------------
# Read track info
# -----------------------------------------------------------------------------

    def read_track_info(self, rec = False):
        
        if ( rec is False ):
            rec = self.read_record()
        
        # Maximum supported version
        max_version = 0x00000001
        
        # Get version
        version = rec.version()
        if ( version > max_version ):
            print( '(*error*) ZDF: Tracks info version is higher than supported.' , file=sys.stderr)
            print( '(*error*) ZDF: Please update the code to a newer version.' , file=sys.stderr)
            return False
        
        info = ZDF_Tracks_Info()

        info.name = rec.name
        info.label    = self.__read_string()
        info.ntracks = self.__read_uint32()
        info.ndump   = self.__read_uint32()
        info.niter   = self.__read_uint32()
        info.nquants = self.__read_uint32()

        for i in range(info.nquants):
            info.quants.append( self.__read_string() )

        for i in range(info.nquants):
            info.qlabels.append( self.__read_string() )

        for i in range(info.nquants):
            info.qunits.append( self.__read_string() )
        
        # Iteration data is not supported so remove it from metadata
        info.nquants -= 1
        
        info.quants.pop(0)
        info.qlabels.pop(0)
        info.qunits.pop(0)

        return info

# -----------------------------------------------------------------------------
# Read dataset
# -----------------------------------------------------------------------------

    def read_dataset(self, rec = False):

        if ( rec is False ):
            rec = self.read_record()
            
        if ( self.record_type(rec.id) != 'dataset' ):
            print( '(*error*) ZDF: Expected dataset record but found {} instead.'.format(self.record_type(rec.id)),
                  file=sys.stderr)
            return False
            
        
        # Maximum supported version
        max_version = 0x00000002
        
        # Get version
        version = rec.version()
        if ( version > max_version ):
            print( '(*error*) ZDF: Dataset version is higher than supported.' , file=sys.stderr)
            print( '(*error*) ZDF: Please update the code to a newer version.' , file=sys.stderr)
            return False
        
        # Version 0x0001 includes id tag
        if ( version >= 1 ):
            id = self.__read_uint32()
        else:
            id = 0
        
        data_type = self.__read_int32()
        ndims     = self.__read_uint32()
        nx        = self.__read_uint64_arr(ndims)
        data      = self.__read_arr( data_type, nx )

        return data

# -----------------------------------------------------------------------------
# Read chunked dataset
# -----------------------------------------------------------------------------
    
    def read_cdset(self, rec = False, pos = 0 ):
                
        if ( rec is False ):
            rec = self.read_record()

        if ( self.record_type(rec.id) != 'cdset_start' ):
            print( '(*error*) ZDF: Expected cdset_start record but found {} instead.'.format(self.record_type(rec.id)),
                  file=sys.stderr)
            return False

        # Maximum supported version
        max_version = 0x00000001
        
        # Get version
        version = rec.id & 0x0000FFFF
        if ( version > max_version ):
            print( '(*error*) ZDF: Chunked dataset version is higher than supported.' , file=sys.stderr)
            print( '(*error*) ZDF: Please update the code to a newer version.' , file=sys.stderr)
            return -1
        
        id        = self.__read_uint32()
        data_type = self.__read_int32()
        ndims     = self.__read_uint32()
        nx        = self.__read_uint64_arr(ndims)
                
        size = np.prod(nx)
  
        # Create numpy array
        dt = {
            5 :'int32',
            3 :'uint32',
            7 :'int64',
            8 :'uint64',
            9 :'float32',
           10 :'float64',
        }       
        
        data = np.zeros( np.flip(nx), dtype = dt[data_type] )
        
        chunk_name = "{:#08}-chunk".format(id)
        end_name = "{:#08}-end".format(id)
        
        # Loop until the end record
        cdset_start_end = self.__file.tell()
        
        name = ""
        while ( name != end_name ):
            
            # Read next record
            rec = self.read_record()
            
            # Check if end of file reached
            if ( rec is False ):
                break
            
            name = rec.name
                        
            if ( name == chunk_name ):
                chunk_id = self.__read_uint32()
                
                count  = self.__read_int64_arr(ndims)
                start  = self.__read_int64_arr(ndims)
                stride = self.__read_int64_arr(ndims)
                                
                chunk  = self.__read_arr( data_type, count )
                                
                if ( ndims == 1 ) :
                    data[start[0]:start[0]+count[0]:stride[0]] = chunk
                elif ( ndims == 2 ) :
                    data[start[1]:start[1]+count[1]:stride[1],
                         start[0]:start[0]+count[0]:stride[0]] = chunk
                else :
                    data[start[2]:start[2]+count[2]:stride[2],
                         start[1]:start[1]+count[1]:stride[1],
                         start[0]:start[0]+count[0]:stride[0]] = chunk
            
            else:
                
                self.__record_skip(rec)
        
        # If requested, position file pointer at the end of the cdset_start record
        if ( pos == 1 ):
            self.__file.seek( cdset_start_end )

        return data    

# -----------------------------------------------------------------------------
# Read arbitrary ZDF element
# -----------------------------------------------------------------------------

    def read_element( self, rec = False, name = False, type_id = False):
        if ( rec is False ):
            rec = self.read_record()
        
        if ( name ):
            if (name != rec.name ):
                print("(*warning*) Requested name does not match record name", file=sys.stderr)
                print("(*warning*) expected '{}', found '{}".format(name, rec.name), file=sys.stderr)
                self.__record_skip(rec)
                return False
        else:
            name = rec.name
        
        if ( type_id ):
            if (type_id != self.record_type(rec.id) ):
                print("(*warning*) Requested type does not match record type", file=sys.stderr)
                print("(*warning*) expected '{}', found '{}".format(type_id, self.record_type(rec.id)), file=sys.stderr)
                self.__record_skip(rec)
                return False
        else:
            type_id = self.record_type(rec.id)
        
        if ( type_id == "int" ):
            data = self.__read_int32()
        elif( type_id == "double" ):
            data = self.__read_float64()
        elif( type_id == "string" ):
            data = self.__read_string()
        elif( type_id == "dataset" ):
            data = self.read_dataset( rec = rec )
        elif( type_id == "cdset_start" ):
            data = self.read_cdset( rec = rec )
        elif( type_id == "cdset_chunk" ):
            print("(*warning*) Dataset chunks are not meant to be read directly", file=sys.stderr)
            self.__record_skip(rec)
            data = False
        elif( type_id == "cdset_end" ):
            print("(*warning*) Dataset end marks have no data", file=sys.stderr)
            self.__record_skip(rec)
            data = False
        elif( type_id == "iteration" ):
            data = self.read_iteration( rec = rec )
        elif( type_id == "grid_info" ):
            data = self.read_grid_info( rec = rec )
        elif( type_id == "part_info" ):
            data = self.read_part_info( rec = rec )
        elif( type_id == "track_info" ):
            data = self.read_track_info( rec = rec )
        else:
            print("(*warning*) Unknown element type, skipping", file=sys.stderr)
            self.__record_skip(rec)
            data = False
        
        return data

# -----------------------------------------------------------------------------
# Read particle data
# -----------------------------------------------------------------------------
    
    def read_part_data( self, quants ):
        
        # Read all quantities
        data = dict()
        
        for q in quants:
            rec  = self.read_record()
            
            # Sanity check, this should never happen
            if ( rec.name != q ):
                print("(*error*) Expecting {} record, {} found".format(
                    q, rec.name ) )
                      
            # Particle data can be either a standard dataset or a chunked dataset
            type_id = self.record_type(rec.id)
            if ( type_id == "dataset" ):
                data[q] = self.read_dataset( rec = rec )
            elif ( type_id == "cdset_start" ):
                # Read chunked dataset and position file pointer after the
                # cdset_start record
                data[q] = self.read_cdset( rec = rec, pos = 1 )
            else:
                print("(*error*) Unable to read particle data, {} record found".format(type_id),
                      file=sys.stderr)
                data[q] = None
        
        return data

# -----------------------------------------------------------------------------
# Read track data
# -----------------------------------------------------------------------------
    
    def read_track_data( self, trackInfo ):
        
        rec  = self.read_record()
        if ( rec.name != 'itermap' ):
                print("(*error*) Expecting itermap record, {} found".format(rec.name ) ,
                      file=sys.stderr)
        itermap = self.read_cdset( rec = rec, pos = 1 )

        rec  = self.read_record()
        if ( rec.name != 'data' ):
                print("(*error*) Expecting data record, {} found".format(rec.name ) ,
                      file=sys.stderr)
        data = self.read_cdset( rec = rec )
        
        trackNp = np.zeros( trackInfo.ntracks, dtype = '<i8' )
        for i in range( itermap.shape[0] ):
            trackID = itermap[i,0]-1
            npoints = itermap[i,1]
            trackNp[trackID] += npoints
        
        trackData = [None] * trackInfo.ntracks
        for i in range( trackInfo.ntracks ):
            trackData[i] = np.zeros( [ trackNp[i], trackInfo.nquants ], dtype = '<f4' )
            trackNp[i] = 0
                    
        idx = 0
        for i in range( itermap.shape[0] ):
            trackID = itermap[i,0]-1
            npoints = itermap[i,1]
            (trackData[ trackID ])[ trackNp[trackID] : trackNp[trackID]+npoints, : ] = data[ idx : idx+npoints, : ]
            trackNp[trackID] += npoints
            idx += npoints

        return trackData

# -----------------------------------------------------------------------------
# Retrieve list of file contents / print file contents
# -----------------------------------------------------------------------------
    
    def list(self, printRec = False):
        # Position file after magic number
        self.__file.seek(4)

        rec_list = []
        while True:
            rec = self.read_record(skip=True)
            if (rec is False):
                break
            else:
                rec_list.append(rec)

        if (printRec and (len(rec_list) > 0)):
            print('Position     Size(bytes)  Type         Name')
            print('-----------------------------------------------------')
            for rec in rec_list:
                print('{:#010x}   {:#010x}   {:11}  {}'.format(
                    rec.pos, rec.len,
                    self.record_type(rec.id), rec.name)
                     )

        return rec_list
    
# -----------------------------------------------------------------------------
# High level interfaces
# -----------------------------------------------------------------------------

class ZDF_Info:
    """ZDF File information"""
    def __init__(self):
        self.type = ""
        self.grid = None
        self.particles = None
        self.tracks = None
        self.iteration = None


def info( file_name ):
    # Open file
    zdf = ZDFfile( file_name )

    # Check file type and read metadata
    info = ZDF_Info()
    info.type = zdf.read_string()
    if ( info.type == "grid" ):
        info.grid = zdf.read_grid_info()
        info.iteration = zdf.read_iteration()
    elif ( info.type == "particles" ):
        info.particles = zdf.read_part_info()
        info.iteration = zdf.read_iteration()
    elif ( info.type == "tracks-2" ):
        info.tracks = zdf.read_track_info()
    else:
        print("File is not a valid ZDF grid, particles or tracks file", file=sys.stderr)
        zdf.close()
        return False

    # Close file
    zdf.close()

    return info

def read( file_name ):
    # Open file
    zdf = ZDFfile( file_name )

    # Check file type and read metadata
    info = ZDF_Info()

    info.type = zdf.read_string()
    if ( info.type == "grid" ):
        info.grid = zdf.read_grid_info()
        info.iteration = zdf.read_iteration()
        data = zdf.read_element()
    elif ( info.type == "particles" ):
        info.particles = zdf.read_part_info()
        info.iteration = zdf.read_iteration()
        data = zdf.read_part_data( info.particles.quants )
    elif ( info.type == "tracks-2" ):
        info.tracks = zdf.read_track_info()
        data = zdf.read_track_data( info.tracks )
    else:
        print("File is not a valid ZDF grid, particles or tracks file", file=sys.stderr)
        zdf.close()
        return False

    # Close file
    zdf.close()

    return (data,info)

def list(file_name):
    zdf = ZDFfile(file_name)
    zdf.list(True)
    zdf.close()

