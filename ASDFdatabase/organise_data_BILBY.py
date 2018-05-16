import os
from os.path import join, isdir, basename
import sys
import json
import glob
from obspy import UTCDateTime
# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework nsame
virt_net = '_ANU'

# FDSN network identifier2
FDSNnetwork = '7W(2008-2010)'

# =========================================================================== #

path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')
path_cleaned_out = join(path_DATA, 'raw_data_cleaned.json')




if os.path.exists(path_cleaned_out):
    os.remove(path_cleaned_out)

data_cleaned_dict = {}

data_cleaned_dict["ProblemFiles"] = []


# Get a list of service directories
service_dir_list = glob.glob(path_DATA + '*')

# iterate through service directories
for service in service_dir_list:

    if not isdir(service):
        continue

    # ignore other temp directories
    if basename(service) in ["temp"]:
        continue

    print '\r Processing: ', basename(service)


    stn_dir_list = glob.glob(service + '/*')

    for stn_dir in stn_dir_list:
        if not isdir(stn_dir):
            continue

        # if not basename(stn_dir) == "BL19":
        #     continue

        print '\r Processing Station: ', basename(stn_dir)


        day_dir_list = glob.glob(stn_dir + '/*')

        # print(day_dir_list)

        # iterate through station directories
        for j, day_path in enumerate(day_dir_list):
            if not isdir(day_path):
                continue
            day_name = basename(day_path)

            print '\r Working on Day: ', str(j+1), ' of ', len(day_dir_list),
            sys.stdout.flush()

            # get miniseed files
            seed_files = glob.glob(join(day_path, '*.D*'))  # '*miniSEED/*.mseed*'))

            if len(seed_files) == 0:
                continue

            # if not day_name in ["023","125"]:
            #     continue



            # Iterate through the miniseed files, fix the header values and add waveforms
            for _i, filename in enumerate(seed_files):
                # print "\r     Parsing miniseed file ", _i + 1, ' of ', len(seed_files), ' ....',
                # sys.stdout.flush()

                station = basename(filename)[0:4]
                year = basename(filename)[4:6]
                month = basename(filename)[6:8]
                day = basename(filename)[8:10]


                try:
                    file_date = UTCDateTime(year=int(year), month=int(month), day=int(day))
                except (ValueError) as e:
                    data_cleaned_dict["ProblemFiles"].append(filename + '\t' + "ValueError\n")
                    continue

                if not station == basename(stn_dir):
                    #could be because filename is for "XX station or something
                    data_cleaned_dict["ProblemFiles"].append(filename + '\t' + "MismatchStnName\n")
                    continue



                julian_day = str(file_date.julday)

                # print(station,year,month,day, str(file_date.julday))


                if station in data_cleaned_dict.keys():
                    if year in data_cleaned_dict[station].keys():
                        if julian_day in data_cleaned_dict[station][year].keys():
                            data_cleaned_dict[station][year][julian_day].append(filename)

                        else:
                            data_cleaned_dict[station][year][julian_day] = [filename]

                    else:
                        data_cleaned_dict[station][year] = {julian_day: [filename]}

                else:
                    data_cleaned_dict[station] = {year: {julian_day: [filename]}}


# print(data_cleaned_dict)

# for station in data_cleaned_dict.keys():
#     print("Data for Station:"+ station)
#     for year in data_cleaned_dict[station].keys():
#         print ("\tYear:"+ year)
#         for julian_day in data_cleaned_dict[station][year].keys():
#             print("\t\tJulian Day: "+ str(julian_day) + "Contains " + str(len(data_cleaned_dict[station][year][julian_day])) + " files")



with open(path_cleaned_out, 'w') as fp:
    json.dump(data_cleaned_dict, fp)
