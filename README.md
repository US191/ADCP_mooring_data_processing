# ADCP_mooring_data_processing
Program for processing ADCP mooring data

Open template_get_adcp_data.m

Input parameters specific to the mooring:

- Raw file location (file with .000 extension)
- Mooring name (ex: '0-10W')
- Mooring longitude and latitude (ex: 00+00/60)
- Adcp serial number (ex: 15258)
- Adcp type (ex: '150 khz Quartermaster')
- Adcp direction (ex: 'up') % upward-looking 'up', downward-looking 'dn'
- Adcp depth (ex: 150) % nominal instrument depth
- Mean of magnetic deviation at time of deployment and time of recovery if ADCP was not set up 
to correct for magnetic deviation internally ("EA0" code in configuration file). Use http://www.ngdc.noaa.gov/geomag-web/#declination 
(ex: rot=-(9+37/60+9+29/60)/2)
- First and last indices when instrument was at depth (you can do this by plotting 'raw.pressure'
- Range of bins which cover the surface reflection  (ex: sbins= [17:28])

Processing steps

- Data reading
- Correction of magnetic deviation
- Exclude data with percent good below 20%
- Calculate depth of each bin
- If ADCP is upward-looking a depth correction can be inferred from the surface reflection, which is done in adcp_surface_fit
- Remove bad data if ADCP is looking upward : depth below the surface which is contaminated by the surface reflection and and bad velocities closed to the surface

Saving data
