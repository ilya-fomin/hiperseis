# import time
# from os.path import join, exists
# from os import mkdir
# import warnings
#
# import glob
#
# from mpi4py import MPI
#
# import pyasdf
#
# # warnings.filterwarnings("errors")
#
# code_start_time = time.time()
#
#
# # =========================== User Input Required =========================== #
#
# # Path to the data
# data_in_path = '/g/data/ha3/Passive/'
# data_out_path = "/g/data/ha3/Passive/"
#
# # IRIS Virtual Network name
# virt_net = '_AusArray'
#
# # FDSN network identifier
# FDSNnetwork = 'OA'
#
# # =========================================================================== #
#
#
# path_raw_DATA = join(data_in_path, virt_net, FDSNnetwork, 'raw_DATA/')
# ASDF_path_out = join(data_out_path, virt_net, FDSNnetwork, 'ASDF_new')
#
# #
# # def split_list(lst, npartitions):
# #     k, m = divmod(len(lst), npartitions)
# #     return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# #
# #
# # def process(r, seed_path):
# #     print('On rank {}, processing number folders: {}'.format(r, seed_path))
# #
# #     # function to create the ASDF waveform ID tag
# #     def make_ASDF_tag(tr, tag):
# #         # def make_ASDF_tag(ri, tag):
# #         data_name = "{net}.{sta}.{loc}.{cha}__{start}__{end}__{tag}".format(
# #             net=tr.stats.network,
# #             sta=tr.stats.station,
# #             loc=tr.stats.location,
# #             cha=tr.stats.channel,
# #             start=tr.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S"),
# #             end=tr.stats.endtime.strftime("%Y-%m-%dT%H:%M:%S"),
# #             tag=tag)
# #         return data_name
# #
# #     # function to make a number into a 4 digit string with leading zeros
# #     def make_fourdig(a):
# #         if len(a) == 1:
# #             return '000' + a
# #         elif len(a) == 2:
# #             return '00' + a
# #         elif len(a) == 3:
# #             return '0' + a
# #         return a
# #
# #     return()
#
# #
# # comm = MPI.COMM_WORLD # communicator
# # nproc = comm.Get_size() # number for processors
# # rank = comm.Get_rank() # index of current processor
# #
# #
# # # Get a list of service directories
# # service_dir_list = glob.glob(path_raw_DATA + '*')
# #
# # station_dir_list = []
# #
# # for service_dir in service_dir_list:
# #     # get list of station dirs
# #     service_station_dirs = glob.glob(join(service_dir, "*"))
# #     [station_dir_list.append(x) for x in service_station_dirs]
# #
# # # print(station_dir_list)
# #
# #
# # # now break list of service directories to process across the MPI cluster
# # workload = split_list(station_dir_list, nproc)
# #
# #
# # process(rank, station_dir_list[rank])
#
#
# # Get a list of service directories
# service_dir_list = glob.glob(path_raw_DATA + '*')
#
# station_dir_list = []
#
# for service_dir in service_dir_list:
#     # get list of station dirs
#     service_station_dirs = glob.glob(join(service_dir, "*"))
#     [station_dir_list.append(x) for x in service_station_dirs]
#
# # print(station_dir_list)
#
# complete_station_dir_list = []
#
# comm = MPI.COMM_SELF # communicator
# nproc = comm.Get_size() # number for processors
# rank = comm.Get_rank() # index of current processor
#
#
# if rank == 0:
#     # master node
#     for
#
# else:
#     pass
#
#
#
# print("\n")
# exec_time = time.time() - code_start_time
#
# exec_str = "--- Execution time: {} seconds ---".format(exec_time)
#
# print(exec_str)

from mpi4py import MPI
import datetime, sys, numpy, time

################ Set up MPI variables ################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
status = MPI.Status()

################ Master code ################

icomms = []
script = 'ASDF_build_minimus_service_MPI_child.py'
for child_rank in range(size):
    if child_rank == 0:
        pass
    else:
        try:
           print('Trying to spawn child process {}...'.format(child_rank))
           icomm = MPI.COMM_SELF.Spawn(sys.executable, args=[script], maxprocs=1, root=0)
           icomm.send("I am a child_rank {}".format(child_rank), dest=0, tag=11)
           icomms.append(icomm)
           print('Spawned a child.')
        except: ValueError('Spawn failed to start.')

completed = False
while not completed and icomms:
    for icomm in icomms:
        if icomm.Iprobe(source=0, tag=MPI.ANY_TAG):
            print('A child responded...')
            solved\
                = icomm.recv(source=0, tag=MPI.ANY_TAG)
            icomm.Disconnect()
            icomms.remove(icomm)
            if solved: break
    if not solved:
        print('spawns doing some work...')
        time.sleep(1)
# make sure all pending sends get matched
for icomm in icomms:
    icomm.recv(source=0, tag=MPI.ANY_TAG)
    icomm.Disconnect()
print('received solution: %d' % solved)