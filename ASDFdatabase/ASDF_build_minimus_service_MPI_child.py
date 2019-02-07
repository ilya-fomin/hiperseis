from mpi4py import MPI
import datetime, sys, numpy

################ Set up MPI variables ################

icomm = MPI.Comm.Get_parent()
comm = MPI.COMM_WORLD
irank = comm.Get_rank()
rank = comm.Get_rank()

running = True
while running:
    data = None
    data = icomm.recv(source=0, tag=11)
    if data:
        print('Trying to send %s from worker rank %d to %d' % (data, rank, irank))
        icomm.send(data, dest=0, tag=22)
        break
print('Worker on rank %d done.' % rank)
icomm.Disconnect()