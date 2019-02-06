import time
from os.path import join, exists
from os import mkdir
import warnings

import glob

from mpi4py import MPI

warnings.filterwarnings("error")

code_start_time = time.time()


# =========================== User Input Required =========================== #

# Path to the data
data_in_path = '/Volumes/SeiOdyssey1/Passive/'
data_out_path = "/Users/ashbycooper/Desktop/new_Passive"

# IRIS Virtual Network name
virt_net = '_AusArray'

# FDSN network identifier
FDSNnetwork = 'OA'

# =========================================================================== #


path_raw_DATA = join(data_in_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_out_path, virt_net, FDSNnetwork, 'ASDF')

if not exists(ASDF_path_out):
    mkdir(ASDF_path_out)


def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(npartitions)]


def process(r, workload):
    print 'On rank %d, processing folders: %s'%(r, str(workload))

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


comm = MPI.COMM_WORLD # communicator
nproc = comm.Get_size() # number for processors
rank = comm.Get_rank() # index of current processor


# Get a list of service directories
service_dir_list = glob.glob(path_raw_DATA + '*')


# now break list of service directories to process across the MPI cluster
workload = split_list(service_dir_list, nproc)



print("\n")
exec_time = time.time() - code_start_time

exec_str = "--- Execution time: %s seconds ---" % exec_time

print(exec_str)
