#script to collect waveforms corresponding to the automated S-wave picks and save them with an index to a directory
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','..'))
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import autopicks_approx as ap
from obspy.core import *
import random
import re

random.seed(0)

fds=FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True, logger=None)


import numpy as np

loadDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'
saveDir='/g/data/ha3/rlt118/neural-datasets/autopicks/Swaves/'
with open(loadDir+'ensemble.s.verified.txt','r') as f:
    picks=f.readlines()
pickCtr=len(picks)

for index in range(0,pickCtr):
    pick=picks[index].split()
    if pick[0][0]=='#': #this line is commented
        continue
    net=pick[6]
    st=pick[7]
    ch=pick[8]
    time=UTCDateTime(float(pick[9]))

    starttime=time-20
    endtime=time+40



    if ch=='00T':

        baz=float(pick[14])
        #get all the horizontal short-period channels at loc '00' or '' (assuming the 00 in 00T refers to
        #loc)
        try:

            inv=fds.get_stations(starttime,endtime,network=net,station=st)
            #find appropriate channels in the inventory (ASDF returns list, not inventory type)
            for chan in inv:
                if re.search('^[HBS].[N1]$',chan[3]):
                    break
            if not re.search('^[HBS].[N1]$',chan[3]):
                raise Exception('no appropriate horizontal channels found')

            stream=fds.get_waveforms(net,st,chan[2],chan[3],starttime,endtime)
            wf1=stream[0]
            wf2=fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'E',starttime,endtime)[0]
            #this is a hack and wrong. Hopefully 1 and 2 are only off from N and E by a few degrees max
            wf1.id=wf1.id[:-1]+'N'
            wf2.id=wf2.id[:-1]+'E'
            stream=Stream(traces=[wf1,wf2])
        except Exception as e:
            #handle data missing
            print >> sys.stderr, e
            stream=None
            
        wf=None
        if stream:
            try:
                stream=stream.rotate(method='NE->RT',back_azimuth=baz)
                for trywf in stream:
                    if trywf.stats['channel'][-1]=='T':
                        wf=trywf
                        break
            except Exception as e:
                print >>sys.stderr, e
                wf=None

    else:
        try:
            #just grab the channel that the pick was on
            inv=fds.get_stations(starttime,endtime,network=net,station=st,channel=ch)
            loc=inv[0][2]
            wf=fds.get_waveforms(net,st,loc,ch,starttime,endtime)[0]
        except Exception as e:
            print >> sys.stderr, e
            wf=None
    if wf and len(wf)>100:
        wf.resample(100)
        wf.detrend()
        wf.normalize()
        #save/plot the thing here
        wf.write(saveDir+str(index)+'.pkl',format="PICKLE")
        wf.plot(outfile=saveDir+str(index)+'.png')
    else:
        print >> sys.stderr, 'Waveform not found. Skipping pick...'
        continue


