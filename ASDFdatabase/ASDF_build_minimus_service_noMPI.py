import time
from os.path import join, exists, basename, dirname
from os import mkdir, remove
import warnings
import json
import random
import sys

import glob

import pyasdf
from obspy import read


code_start_time = time.time()


# =========================== User Input Required =========================== #

# Path to the data
data_in_path = '/g/data/ha3/Passive/'
data_out_path = "/g/data/ha3/Passive/"

# IRIS Virtual Network name
virt_net = '_AusArray'

# FDSN network identifier
FDSNnetwork = 'OA'

# =========================================================================== #

ASDF_path_out = join(data_out_path, virt_net, FDSNnetwork, 'ASDF_new')
JSON_process_filename = join(ASDF_path_out, FDSNnetwork+"_process_file.json")


def process(rank, process_path_dict):

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


    raw_DATA_path = process_path_dict["raw_DATA_path"]
    asdf_filename = process_path_dict["ASDF_filename"]
    json_filename = process_path_dict["ASDF_json_path"]
    asdf_log_filename = process_path_dict["ASDF_logfile_path"]

    waveforms_added = 0

    # get the station name from the directory name
    station_name = basename(raw_DATA_path).split("_")[1]

    # keys and info list for waveforms taht are added to ASDF
    waveform_keys_list = []
    waveform_info_list = []

    # create the log file
    ASDF_log_file = open(asdf_log_filename, 'w')
    # get miniseed files
    seed_files = glob.glob(join(raw_DATA_path, '*miniSEED*/*.mseed*'))

    print("Process {} recieved: {} - No. SEED files: {}".format(str(rank), basename(asdf_filename), len(seed_files)))



    if len(seed_files) == 0:
        return (basename(asdf_filename), 0)

    # now open ASDF and add new data
    ds = pyasdf.ASDFDataSet(asdf_filename)

    # Iterate through the miniseed files, fix the header values and add waveforms
    for _i, seed_filename in enumerate(seed_files):
        # print("\rParsing miniseed file ", _i + 1, ' of ', len(seed_files), ' ....',)
        # sys.stdout.flush()

        if not basename(seed_filename) in ["Sensor0TemprCRou_0000000010_00008.mseed",
                                           "Sensor0SeismoZSm_0000000200_00015.mseed",
                                           "Sensor0HumidBRou_0000000010_00006.mseed",
                                           "ClkDACVoltageRou_0000000010_00007.mseed",
                                           "Sensor0SeismoNSm_0000000200_00021.mseed",
                                           "Sensor0SeismoESm_0000000200_00014.mseed",
                                           "Sensor0SeismoZSm_0000000200_00033.mseed",
                                           "Sensor0VoltageRo_0000000010_00002.mseed",
                                           "Sensor0TemprCRou_0000000010_00004.mseed"]:
            continue

        print("Reading Stream")

        st_read_start = time.time()

        try:
            # Read the stream
            st = read(seed_filename)

        except Exception as e:
            # there is something wrong with the miniseed file
            ASDF_log_file.write(seed_filename + '\t' + str(type(e)) + "\t" + str(e.args) + "\n")
            continue

        st_read_exec_time = time.time() - st_read_start
        print("\tExec time: {} seconds".format(str(st_read_exec_time)))

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
            # use the first two characters as network code. Temporary networks have start and end year as well
            new_net = FDSNnetwork[:2]
            # overwrite network code in miniseed header
            tr.stats.network = new_net

            # Station Name: use directory name
            orig_station = tr.stats.station
            new_station = station_name
            # overwrite station code in miniseed header
            tr.stats.station = new_station

            # Channel use miniseed
            orig_chan = tr.stats.channel
            new_chan = orig_chan

            # Location Code: use miniseed
            orig_loc = tr.stats.location
            new_loc = orig_loc

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp

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
                         "seed_path": str(dirname(seed_filename)),
                         "seed_filename": str(basename(seed_filename)),
                         "log_filename": ""}

            print("\t\tTrace Length: {}".format(tr.stats.npts))

            print("\t\tAdding to ASDF")

            asdf_add_start = time.time()

            try:
                # Add waveform to the ASDF file
                # if station_name == "TEST3":
                #     print(tr)
                ds.add_waveforms(tr, tag="raw_recording")
            except Exception as e:
                # trace already exist in ASDF file!
                ASDF_log_file.write(seed_filename + '\t' + str(type(e)) + "\t" + str(e.args) + "\n")
                continue

            asdf_add_exec_time = time.time() - asdf_add_start
            print("\t\t\tExec time: {} seconds".format(str(asdf_add_exec_time)))


            waveform_keys_list.append(str(ASDF_tag))
            waveform_info_list.append(temp_dict)

            tr = None
        st = None

    print("\tFinished Iterating Traces")

    del ds

    waveform_dictionary = dict(zip(waveform_keys_list, waveform_info_list))

    print("\tDumping Waveform Dictionary")

    with open(json_filename, 'w') as fp:
        json.dump(waveform_dictionary, fp)

    return(basename(asdf_filename), waveforms_added)




# read in the process json
f = open(JSON_process_filename, 'r')
process_dict = json.load(f)


total_no_tasks = len(process_dict.keys())

task_index = 0

while task_index < total_no_tasks:
    result = process(0, process_dict[str(task_index)])

    print("{} Waveforms were added to File: {}".format(result[1], result[0]))

    task_index += 1



exec_time = time.time() - code_start_time

exec_str = "--- Execution time: %s seconds ---" % exec_time

print("--- All Done ---")
print(exec_str)