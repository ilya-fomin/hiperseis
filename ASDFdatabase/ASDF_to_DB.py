import json
import pyasdf
from os.path import join, exists
import sys
from os import remove
from query_input_yes_no import query_yes_no
import time
from obspy import UTCDateTime
code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_ANU'

# FDSN network identifier2
FDSNnetwork = '7W(2008-2010)'

# =========================================================================== #

ASDF_path = join(data_path, virt_net, FDSNnetwork, 'ASDF')
ASDF_filename = join(ASDF_path, FDSNnetwork+".h5")

# JSON filename for network
JSON_out = join(ASDF_path, FDSNnetwork + '_raw_dataDB.json')


# query the user to overwrite JSON database file or not
if exists(JSON_out):
    delete_queary = query_yes_no("Remove Existing JSON database file?")
    if delete_queary == 'yes':
        # removing existing SQLdb
        remove(JSON_out)
    elif delete_queary == 'no':
        sys.exit(0)


def make_three_char(a):
    if len(a) == 1:
        return '  ' + a
    elif len(a) == 2:
        return ' ' + a
    return a

#open the ASDF file
ds = pyasdf.ASDFDataSet(ASDF_filename)

# get list of stations
sta_list = ds.waveforms.list()

no_stns = len(sta_list)

waveforms_added = 0
keys_list = []
info_list = []

for _i, sta in enumerate(sta_list):
    #
    # if not sta == "7W.BL02":
    #     continue

    print("\nWorking on Station: " + sta + " (" + make_three_char(str(_i+1)) + " of " + str(no_stns) + ")")

    sta_accessor = ds.waveforms[sta]

    waveforms_list = sta_accessor.list()

    no_waveforms = len(waveforms_list)

    for _j, ASDF_tag in enumerate(waveforms_list):

        # ignore the station XML entry
        if ASDF_tag == "StationXML":
            continue

        print "\r  Parsing waveform " + make_three_char(str(_j + 1)) + ' of ' + str(no_waveforms-1) +' ....',
        sys.stdout.flush()

        ASDF_tag_split = ASDF_tag.split("__")
        id_split = ASDF_tag_split[0].split(".")

        network = id_split[0]
        station = id_split[1]
        channel = id_split[3]
        location = id_split[2]
        starttime = UTCDateTime(ASDF_tag_split[1])
        endtime = UTCDateTime(ASDF_tag_split[2])

        starttime = UTCDateTime(ASDF_tag_split[1]).timestamp
        endtime = UTCDateTime(ASDF_tag_split[2]).timestamp

        temp_dict = {"tr_starttime": starttime,
                     "tr_endtime": endtime,
                     "orig_network": str("__"),
                     "new_network": str(network),
                     "orig_station": str("__"),
                     "new_station": str(station),
                     "orig_channel": str("__"),
                     "new_channel": str(channel),
                     "orig_location": str("__"),
                     "new_location": str(location),
                     "seed_path": str("ASDF"),
                     "seed_filename": str("ASDF"),
                     "log_filename": ""}

        keys_list.append(str(ASDF_tag))
        info_list.append(temp_dict)

        waveforms_added += 1



big_dictionary = dict(zip(keys_list, info_list))

with open(JSON_out, 'w') as fp:
    json.dump(big_dictionary, fp)


del ds

exec_time = time.time() - code_start_time

print("\n")
exec_str = "--- Execution time: %s seconds ---" % exec_time
added_str = '--- Added ' + str(waveforms_added) + ' waveforms to ASDF and JSON database files ---'

print exec_str
print added_str




