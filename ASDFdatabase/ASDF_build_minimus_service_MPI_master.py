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

# from mpi4py import MPI
# import sys, time
#
# ################ Set up MPI variables ################
#
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
# name = MPI.Get_processor_name()
# status = MPI.Status()
#
# ################ Master code ################
#
# print("MPI settings")
#
# print(rank)
# print(size)
# print(name)
# print(status)
#
# print("----------")
#
# icomms_dict = {}
# script = 'ASDF_build_minimus_service_MPI_spawn.py'
#
# if rank == 0:
#     for _i in range(size-1):
#         child_rank = _i + 1
#         try:
#            print('sending info to child process {}...'.format(child_rank))
#            icomm = MPI.COMM_SELF.Spawn(sys.executable, args=[script], maxprocs=1, root=0)
#            icomm.send("I am a child_rank {}".format(child_rank), dest=0, tag=11)
#            icomms_dict[child_rank] = {"icomm": icomm, "status": True}
#            print('Spawned a child.')
#         except: ValueError('Spawn failed to start.')
#
# completed_dict = {}
# work_finished = False
#
# while not work_finished:
#     for mpi_rank, icomm_info in icomms_dict.items():
#         icomm = icomm_info["icomm"]
#
#         if icomm.Iprobe(source=0, tag=MPI.ANY_TAG):
#             print('{} child responded...'.format(mpi_rank))
#             rec_msg = icomm.recv(source=0, tag=MPI.ANY_TAG)
#             print(icomm.Get_rank())
#             icomm.Disconnect()
#             icomms_dict[mpi_rank]["status"] = False
#
#         else:
#             print('{} child doing some work...'.format(mpi_rank))
#
#     if all(value["status"] == False for value in icomms_dict.values()):
#         work_finished = True
#
#     else:
#         # sleep and wait another period of time to recheck for msg
#         time.sleep(10)
#
#
# print('work complete')

from mpi4py import MPI


def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')

# Initializations and preliminaries
comm = MPI.COMM_WORLD  # get MPI communicator object
size = comm.size  # total number of processes
rank = comm.rank  # rank of this process
status = MPI.Status()  # get MPI status object

if rank == 0:
    # Master process executes code below
    tasks = range(10 * size)
    task_index = 0
    num_workers = size - 1
    closed_workers = 0
    print("Master starting with %d workers" % num_workers)
    while closed_workers < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # Worker is ready, so send it a task
            if task_index < len(tasks):
                comm.send(tasks[task_index], dest=source, tag=tags.START)
                print("Sending task %d to worker %d" % (task_index, source))
                task_index += 1
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            results = data
            print("Got data from worker %d" % source)
        elif tag == tags.EXIT:
            print("Worker %d exited." % source)
            closed_workers += 1

    print("Master finishing")

else:
    # Worker processes execute code below
    name = MPI.Get_processor_name()
    print("I am a worker with rank %d on %s." % (rank, name))
    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        if tag == tags.START:
            # Do the work here
            result = task ** 2
            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            break

    comm.send(None, dest=0, tag=tags.EXIT)