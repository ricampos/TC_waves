#!/usr/bin/env python3
import s3fs
import eccodes
import numpy as np

# Initialize s3fs for anonymous access
fs = s3fs.S3FileSystem(anon=True)

# GEFS files
grib_file = 'noaa-gefs-pds/gefs.20241106/00/wave/gridded/gefs.wave.t00z.c00.global.0p25.f009.grib2'
idx_file = 'noaa-gefs-pds/gefs.20241106/00/wave/gridded/gefs.wave.t00z.c00.global.0p25.f009.grib2.idx'

# Read the index file to find SWH message location
with fs.open(idx_file, 'r') as f:
    idx_lines = f.read().strip().split('\n')

# Find significant wave height message
for i, line in enumerate(idx_lines):
    parts = line.split(':')
    if len(parts) >= 4 and 'HTSGW' in line:
        byte_offset = int(parts[1])
        
        # Get next offset for message size
        if i + 1 < len(idx_lines):
            next_parts = idx_lines[i + 1].split(':')
            next_offset = int(next_parts[1])
            message_size = next_offset - byte_offset
        else:
            message_size = None
        
        # Read only the specific GRIB message
        with fs.open(grib_file, 'rb') as f:
            f.seek(byte_offset)
            grib_message = f.read(message_size) if message_size else f.read()
        
        # Decode GRIB message
        msg = eccodes.codes_new_from_message(grib_message)
        
        # Get missing value info
        missing_value = eccodes.codes_get(msg, 'missingValue')
        values = eccodes.codes_get_values(msg)
        
        # Replace missing values with NaN
        values[values == missing_value] = np.nan
        
        # Calculate max significant wave height
        max_swh = np.nanmax(values)
        
        print(f"Maximum significant wave height: {max_swh:.2f} m")
        print(f"Data shape: {values.shape}")
        print(f"Valid data points: {np.sum(~np.isnan(values))}")
        print(f"Downloaded only {len(grib_message)} bytes instead of full file!")
        
        eccodes.codes_release(msg)
        break

