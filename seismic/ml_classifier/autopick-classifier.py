from data_harvester.autopicks import pickLoader

import model
import numpy as np
import autopick_data

network=model.shakenet(pretrained_weights='shakenet-model.hdf5')

datagen=autopick_data.autoGenerator(10)

outfile='/g/data/ha3/rlt118/neural-datasets/autopicks/ensemble.s.verified-slow.txt'
infile='/g/data/ha3/rlt118/neural-datasets/autopicks/ensemble.s.txt'

corrCtr=0

with open(outfile,'w') as f, open(infile,'r') as g:
    #initialise the file with the same format as the input    
    f.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n")

    results=network.predict_generator(datagen)
    for i in range(0,len(results)):
        pick=g.readline()
        while pick[0]=='#':
            pick=g.readline()
        if results[i][0]>=0.5:
            corrCtr+=1
            f.write(pick)
            

    

