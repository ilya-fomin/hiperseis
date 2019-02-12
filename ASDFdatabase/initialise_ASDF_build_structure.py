from os.path import join, exists, basename, dirname
from os import mkdir, remove
import glob
import json

import pyasdf


# =========================== User Input Required =========================== #

# Path to the data
data_in_path = '/g/data/ha3/Passive/'
data_out_path = "/g/data/ha3/Passive/"

# IRIS Virtual Network name
virt_net = '_AusArray'

# FDSN network identifier
FDSNnetwork = 'OA'

process_service_list = ["AusArray_year1_pickup_Nov18"]

# processed_stations_list = ["OA_CC25", "OA_BW25", "OA_CC26", "OA_BZ26", "OA_CC24", "OA_CD25", "OA_CE25", "OA_CA20", "OA_CB22"]
processed_stations_list = ["OA_CC25", "OA_BW25", "OA_CC26", "OA_BZ26", "OA_CC24", "OA_CD25", "OA_CE25", "OA_CA20", "OA_CB22", "OA_CB23", "OA_BS26", "OA_BZ20", "OA_BZ21", "OA_BY22"]
# =========================================================================== #

path_raw_DATA = join(data_in_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_out_path, virt_net, FDSNnetwork, 'ASDF_new')
JSON_process_filename = join(ASDF_path_out, FDSNnetwork+"_process_file.json")

if not exists(join(ASDF_path_out)):
    mkdir(join(ASDF_path_out))


# Get a list of service directories
service_dir_list = glob.glob(path_raw_DATA + '*')

process_dict = {}
processing_time_stn_service_list = []
num_seed_files_list = []

service_station_index = 0

for service_dir in service_dir_list:
    service_name = basename(service_dir)

    # if not service_name in process_service_list:
    #     continue

    if not exists(join(ASDF_path_out, service_name)):
        mkdir(join(ASDF_path_out, service_name))

    # get list of station dirs
    service_station_dirs = glob.glob(join(service_dir, "*"))

    for service_station_path in service_station_dirs:
        station_name = basename(service_station_path)

        # if not station_name in processed_stations_list:
        #     continue

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

        # get the number of miniseed files for the service_station
        # get miniseed files
        seed_files = glob.glob(join(service_station_path, '*miniSEED*/*.mseed*'))
        # estimate prcessing time for station/service by assuming 12 seconds of time per seed file
        processing_time_stn_service_list.append(len(seed_files)*15)
        num_seed_files_list.append(len(seed_files))

        service_station_index += 1



with open(JSON_process_filename, 'w') as fp:
    json.dump(process_dict, fp)



# claculate processing time by summing proc times for each station/service and adding 10 percent then calculating average
est_proc_time = (max(processing_time_stn_service_list)*1.1)#/len(processing_time_stn_service_list)

print("Estimated Processing Time is: {} seconds".format(est_proc_time))
print("Estimated Wall Time is: {} seconds".format(est_proc_time*1.5))
print(len(processing_time_stn_service_list))
print(processing_time_stn_service_list)
print(sum(processing_time_stn_service_list))
print(num_seed_files_list)