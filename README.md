# ADCP_mooring_data_processing
Program for processing ADCP mooring data

---------------------

Short documentation:

Open template_get_adcp_data.m

Input parameters specific to the mooring:

- Raw file location (file with .000 extension)
- Output directory

- Cruise name
- Mooring name (ex: '0-10W')
- Mooring longitude and latitude ([deg-min])
- clock drift

- Adcp serial number (ex: 15258)
- Adcp type (ex: '150 khz Quartermaster')
- Adcp direction (ex: 'up') % upward-looking 'up', downward-looking 'dn'
- Adcp depth (ex: 150) % nominal instrument depth
- Adcp number if more than 1 Adcp

---------------------

Processing steps

- Data reading
- Correct clock drift (linear drift)
- Determine first and last indices when instrument was at depth
- Correction of magnetic deviation
- Percent good threshold
- Remove bad attitude (pitch-roll) data
- Calculate depth of each bin
- Depth correction from surface reflection (if ADCP is upward-looking) 
- Remove "shadow zone" data (if ADCP is upward-looking) 
- Grid data on a regular vertical grid
- Temporal interpolation to fill gaps
- Remove tide effect
- Interpolate data on a regular temporal grid (6hour)
- Save data as netcdf (version 1)
