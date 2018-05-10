import time
from os.path import join, exists, basename, isdir, dirname
from os import remove, mkdir
from struct import error as StructError
import json
import pyasdf
from pyasdf import ASDFWarning
import warnings
from collections import Counter, defaultdict
from convert_logs.decode_datfile import decode_anulog

import glob

import numpy as np

from obspy import read_inventory, read, UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Site, Channel
from obspy.io.mseed.core import InternalMSEEDReadingWarning

import warnings

import sys
from query_input_yes_no import query_yes_no

warnings.filterwarnings("error")

code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_ANU'

# FDSN network identifier2
FDSNnetwork = '7W(2008-2010)'

# XML file input
XML_in = '/g/data1/ha3/Passive/_ANU/7W(2008-2010)/network_metadata/7W_prelim.xml'

# =========================================================================== #

XML_path_out = join(data_path, virt_net, FDSNnetwork, 'network_metadata')
path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_path, virt_net, FDSNnetwork, 'ASDF')

if not exists(ASDF_path_out):
    mkdir(ASDF_path_out)

# JSON filename for network
JSON_out = join(ASDF_path_out, FDSNnetwork + '_raw_dataDB.json')
# ASDF filename for network
ASDF_out = join(ASDF_path_out, FDSNnetwork + '.h5')
# Logfile output
ASDF_log_out = join(ASDF_path_out, FDSNnetwork + '.log')

keys_list = []
info_list = []
station_name_counter = Counter()
station_name_paras = {}

# remove log file if it exists
if exists(ASDF_log_out):
    remove(ASDF_log_out)

# query the user to overwrite JSON database file or not
if exists(JSON_out):
    delete_queary = query_yes_no("Remove Existing JSON database file?")
    if delete_queary == 'yes':
        # removing existing SQLdb
        remove(JSON_out)
    elif delete_queary == 'no':
        sys.exit(0)

# query the user to overwrite the ASDF database or not
if exists(ASDF_out):
    delete_queary = query_yes_no("Remove Existing ASDF File?")
    if delete_queary == 'yes':
        # removing existing ASDF
        remove(ASDF_out)
    elif delete_queary == 'no':
        sys.exit(0)

# create the log file
ASDF_log_file = open(ASDF_log_out, 'w')

# Create/open the ASDF file
ds = pyasdf.ASDFDataSet(ASDF_out, compression="gzip-3")

# read in station xml file
read_inv = read_inventory(XML_in)

# create empty inventory to add all inventories together
new_inv = Inventory(networks=[], source="Geoscience Australia AusArray")

# create the inventory object for the network
net_inv = Network(code=FDSNnetwork[:2])

# dictionary to keep end date/start date for each station
station_start_end_dict = {}

# dictionary to keep inventory for all stations (default dict)
station_inventory_dict = {}


# function to create the ASDF waveform ID tag
def make_ASDF_tag(tr, tag):
    # def make_ASDF_tag(ri, tag):
    data_name = "{net}.{sta}.{loc}.{cha}__{start}__{end}__{tag}".format(
        net=tr.stats.network,
        sta=tr.stats.station,
        loc=tr.stats.location,
        cha=tr.stats.channel,
        start=tr.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S"),
        end=tr.stats.endtime.strftime("%Y-%m-%dT%H:%M:%S"),
        tag=tag)
    return data_name


# function to make a number into a 4 digit string with leading zeros
def make_fourdig(a):
    if len(a) == 1:
        return '000' + a
    elif len(a) == 2:
        return '00' + a
    elif len(a) == 3:
        return '0' + a
    return a


waveforms_added = 0

# Get a list of service directories
service_dir_list = glob.glob(path_DATA + '*')

# counter for number of logfiles
logfile_counter = 0

# iterate through service directories
for service in service_dir_list:

    if not isdir(service):
        continue

    # ignore other temp directories
    if basename(service) == "temp":
        continue

    print '\r Processing: ', basename(service)


    stn_dir_list = glob.glob(service + '/*')

    for stn_dir in stn_dir_list:
        print '\r Processing Station: ', basename(stn_dir)


        day_dir_list = glob.glob(stn_dir + '/*')

        # iterate through station directories
        for day_path in day_dir_list:
            day_name = basename(day_path)

            # get miniseed files
            seed_files = glob.glob(join(day_path, '*.D*'))  # '*miniSEED/*.mseed*'))

            if len(seed_files) == 0:
                continue

            print '\r Working on Day: ', day_name

            # Iterate through the miniseed files, fix the header values and add waveforms
            for _i, filename in enumerate(seed_files):
                print "\r     Parsing miniseed file ", _i + 1, ' of ', len(seed_files), ' ....',
                sys.stdout.flush()

                try:
                    # Read the stream
                    st = read(filename)

                except (TypeError, StructError, InternalMSEEDReadingWarning) as e:
                    # the file is not miniseed
                    ASDF_log_file.write(filename + '\t' + "TypeError\n")
                    continue

                # iterate through traces in st (there will usually be only one trace per stream,
                # however if there are problems with the miniseed files - like much of the ANU data - there
                # can be more than one trace in a miniseed file (seperated by a large time gap)
                for tr in st:

                    if len(tr) == 0:
                        continue

                    waveforms_added += 1

                    # do some checks to make sure that the network, station, channel, location information is correct

                    # Network Code: the network code in the miniseed header is prone to user error
                    # (i.e. whatever the operator entered into the instrument in the field)
                    orig_net = tr.stats.network
                    # assume trace header is correct
                    new_net = orig_net

                    # # use the first two characters as network code. Temporary networks have start and end year as well
                    # new_net = FDSNnetwork[:2]
                    # # overwrite network code in miniseed header
                    # tr.stats.network = new_net

                    # Station Name: assume header is correct
                    orig_station = tr.stats.station
                    new_station = orig_station

                    # Channel assume header
                    orig_chan = tr.stats.channel
                    new_chan = orig_chan

                    # Location Code: use miniseed
                    orig_loc = tr.stats.location
                    new_loc = orig_loc

                    # ignore miniseed files if station name does not match directory
                    if not(new_station == basename(stn_dir)):
                        ASDF_log_file.write(filename + '\t' + "MismatchStnName\n")
                        continue


                    starttime = tr.stats.starttime.timestamp
                    endtime = tr.stats.endtime.timestamp

                    # see if station is already in start_end dict
                    if new_station in station_start_end_dict.keys():
                        # compare time to start and end times in dict and see if it is earlier/later
                        stored_starttime = station_start_end_dict[new_station][0]
                        stored_endtime = station_start_end_dict[new_station][1]
                        if starttime < stored_starttime:
                            station_start_end_dict[new_station][0] = starttime
                        elif endtime > stored_endtime:
                            station_start_end_dict[new_station][1] = endtime
                    else:
                        station_start_end_dict[new_station] = [starttime, endtime]

                    # The ASDF formatted waveform name [full_id, station_id, starttime, endtime, tag]
                    ASDF_tag = make_ASDF_tag(tr, "raw_recording").encode('ascii')

                    # make a dictionary for the trace that will then be appended to a larger dictionary for whole network
                    temp_dict = {"tr_starttime": starttime,
                                 "tr_endtime": endtime,
                                 "orig_network": str(orig_net),
                                 "new_network": str(new_net),
                                 "orig_station": str(orig_station),
                                 "new_station": str(new_station),
                                 "orig_channel": str(orig_chan),
                                 "new_channel": str(new_chan),
                                 "orig_location": str(orig_loc),
                                 "new_location": str(new_loc),
                                 "seed_path": str(dirname(filename)),
                                 "seed_filename": str(basename(filename)),
                                 "log_filename": ""}

                    # print(new_net,new_station,new_chan,new_loc)

                    # get the inventory object for the channel
                    select_inv = read_inv.select(network=new_net, station=new_station, channel=new_chan, location=new_loc)

                    # print(select_inv)

                    # see if station is already in the station inv dictionary

                    if not new_station in station_inventory_dict.keys():


                        # create 3* channels:

                        z_chan = Channel(code="BHZ", location_code="", depth=0, azimuth=0, dip=90,
                                        start_date=select_inv[0][0].start_date,
                                        end_date=select_inv[0][0].end_date,
                                        sample_rate=50,
                                        clock_drift_in_seconds_per_sample=0,
                                         latitude=select_inv[0][0].latitude,
                                         longitude=select_inv[0][0].longitude,
                                         elevation=select_inv[0][0].elevation)

                        n_chan = Channel(code="BHN", location_code="", depth=0, azimuth=0, dip=0,
                                         start_date=select_inv[0][0].start_date,
                                         end_date=select_inv[0][0].end_date,
                                         sample_rate=50,
                                         clock_drift_in_seconds_per_sample=0,
                                         latitude=select_inv[0][0].latitude,
                                         longitude=select_inv[0][0].longitude,
                                         elevation=select_inv[0][0].elevation)

                        e_chan = Channel(code="BHE", location_code="", depth=0, azimuth=90, dip=0,
                                         start_date=select_inv[0][0].start_date,
                                         end_date=select_inv[0][0].end_date,
                                         sample_rate=50,
                                         clock_drift_in_seconds_per_sample=0,
                                         latitude=select_inv[0][0].latitude,
                                         longitude=select_inv[0][0].longitude,
                                         elevation=select_inv[0][0].elevation)

                        # create the station inventory

                        sta_inv = Station(code=new_station,
                                          creation_date=select_inv[0][0].creation_date,
                                          start_date=select_inv[0][0].start_date,
                                          end_date=select_inv[0][0].end_date,
                                          latitude=select_inv[0][0].latitude,
                                          longitude=select_inv[0][0].longitude,
                                          elevation=select_inv[0][0].elevation,
                                          site=Site(new_station),
                                          channels=[z_chan, n_chan, e_chan])


                        # append it to the station inventory dict
                        station_inventory_dict[new_station] = sta_inv

                    try:
                        # Add waveform to the ASDF file
                        ds.add_waveforms(tr, tag="raw_recording", labels=[basename(service)])
                    except ASDFWarning:
                        # trace already exist in ASDF file!
                        ASDF_log_file.write(filename + '\t' + ASDF_tag + '\t' + "ASDFDuplicateError\n")
                        continue

                    keys_list.append(str(ASDF_tag))
                    info_list.append(temp_dict)

# go through the stations in the station inventory dict and append them to the network inventory
for station, sta_inv in station_inventory_dict.iteritems():
    start_date = UTCDateTime(station_start_end_dict[station][0])
    end_date = UTCDateTime(station_start_end_dict[station][1])

    # print(station)
    # print(start_date, end_date)
    # change the station start/end date
    # get the start/end dates from dict


    # print(sta_inv)

    channel_inventory_list = []
    # fix the channel inventory by fixing the start/end times
    for chan_inv in sta_inv:
        chan = chan_inv.code


        if 'Z' in chan:
            az = '0'
            dip = '90'
        elif 'N' in chan:
            az = '0'
            dip = '0'
        elif 'E' in chan:
            az = '90'
            dip = '0'

        channel_inventory_list.append(
            Channel(code=chan, location_code=chan_inv.location_code, depth=0, azimuth=az, dip=dip,
                    start_date=start_date,
                    end_date=end_date,
                    sample_rate=chan_inv.sample_rate,
                    clock_drift_in_seconds_per_sample=0,
                    latitude=sta_inv.latitude,
                    longitude=sta_inv.longitude,
                    elevation=sta_inv.elevation))

    site = Site(name=station)

    # make the station_level inventory
    new_sta_inv = Station(code=station, creation_date=start_date, termination_date=end_date,
                          start_date=start_date,
                          end_date=end_date,
                          site=site,
                          latitude=sta_inv.latitude,
                          longitude=sta_inv.longitude,
                          elevation=sta_inv.elevation,
                          vault="Transportable Array",
                          channels=channel_inventory_list)

    net_inv.stations.append(new_sta_inv)
    # net_inv.stations.append(sta_inv)

network_start_end = False
# go through station start/end date dict and get the overall start_end date
for key, (start, end) in station_start_end_dict.iteritems():
    if not network_start_end:
        network_start_end = [start, end]
        continue

    if start < network_start_end[0]:
        network_start_end[0] = start
    elif end > network_start_end[1]:
        network_start_end[1] = end

# now add the network start/end date
net_inv.start_date = UTCDateTime(network_start_end[0])
net_inv.end_date = UTCDateTime(network_start_end[1])

# print(net_inv)

# add the network inventory to the complete and updated inventory
new_inv.networks.append(net_inv)

XML_file = join(XML_path_out, FDSNnetwork + '_updated.xml')

if exists(XML_file):
    remove(XML_file)

# write the inventory into the default path
new_inv.write(path_or_file_object=XML_file, format='STATIONXML', validate=True)

# add it to ASDF file
ds.add_stationxml(new_inv)

big_dictionary = dict(zip(keys_list, info_list))

with open(JSON_out, 'w') as fp:
    json.dump(big_dictionary, fp)

del ds
print '\n'

exec_time = time.time() - code_start_time

exec_str = "--- Execution time: %s seconds ---" % exec_time
added_str = '--- Added ' + str(waveforms_added) + ' waveforms to ASDF and JSON database files ---'

print exec_str
print added_str

ASDF_log_file.write(exec_str + '\n')
ASDF_log_file.write(added_str + '\n')

ASDF_log_file.close()
