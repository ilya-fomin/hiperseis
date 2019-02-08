from os.path import join, exists, basename, dirname
from os import mkdir, remove
import glob
import json

import pyasdf


# =========================== User Input Required =========================== #

# Path to the data
data_in_path = '/Volumes/SeiOdyssey1/Passive/'
data_out_path = "/Users/ashbycooper/Desktop/new_Passive/"

# IRIS Virtual Network name
virt_net = '_AusArray'

# FDSN network identifier
FDSNnetwork = 'OA'

# =========================================================================== #

path_raw_DATA = join(data_in_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_out_path, virt_net, FDSNnetwork, 'ASDF')
JSON_process_filename = join(ASDF_path_out, FDSNnetwork+"_process_file.json")

if not exists(join(ASDF_path_out)):
    mkdir(join(ASDF_path_out))


# Get a list of service directories
service_dir_list = glob.glob(path_raw_DATA + '*')

process_dict = {}

service_station_index = 0

for service_dir in service_dir_list:
    service_name = basename(service_dir)
    if not exists(join(ASDF_path_out, service_name)):
        mkdir(join(ASDF_path_out, service_name))

    # get list of station dirs
    service_station_dirs = glob.glob(join(service_dir, "*"))

    for service_station_path in service_station_dirs:
        station_name = basename(service_station_path)
        if not exists(join(ASDF_path_out, service_name, station_name)):
            mkdir(join(ASDF_path_out, service_name, station_name))

        ASDF_filename = join(ASDF_path_out, service_name, station_name, station_name + "_" + service_name + ".h5")
        JSON_out = join(ASDF_path_out, service_name, station_name, station_name + "_" + service_name + '_raw_dataDB.json')
        ASDF_logfile_out = join(ASDF_path_out, service_name, station_name, station_name + "_" + service_name + '.log')

        if exists(ASDF_filename):
            remove(ASDF_filename)
        if exists(JSON_out):
            remove(JSON_out)
        if exists(ASDF_logfile_out):
            remove(ASDF_logfile_out)

        # make the ASDF file (empty)
        ds = pyasdf.ASDFDataSet(ASDF_filename)
        del ds

        process_dict[service_station_index] = {"ASDF_filename": ASDF_filename, "raw_DATA_path": service_station_path, "ASDF_logfile_path":ASDF_logfile_out ,"ASDF_json_path":JSON_out}

        service_station_index += 1



with open(JSON_process_filename, 'w') as fp:
    json.dump(process_dict, fp)