import numpy as np
import os
import pandas as pd
import sys
import time, re

from multiprocessing.dummy import Pool 


import xml.etree.ElementTree as ET
from xml.dom.minidom import parse
from contextlib import contextmanager

@contextmanager
def timer(title):
    t0 = time.time()
    yield
    print("{} - done in {:.2f}s".format(title, time.time() - t0))

def getValue(node):
    if node.childNodes == []:
        return None
    else:
        return node.childNodes[0].nodeValue

def proteomeDB_proteins(path):
    data = pd.DataFrame(columns=['PROTEINNAME','SEARCH_ENGINE','IDENTIFICATION_ID','PRECURSOR_MZ','PRECURSOR_CHARGE','RECALIBRATED_PRECURSOR_MZ',
                            'SEQUENCE','SCORE','DELTA_SCORE','THRESHOLD_SCORE','PEP','PEPTIDE_ID','Q_VALUE','LLD_FDR_CUTOFF',
                             'VARIABLE_MODIFICATION_STRING','FIXED_MODIFICATION_STRING','MODIFICATION_DELTA_MASS','EXPECTED_PEPTIDE_MASS',
                            'OBSERVED_PEPTIDE_MASS','MASS_ERROR_DA','MASS_ERROR_PPM','PROJECT_NAME','EXPERIMENT_NAME',
                            'EXPERIMENT_ID','FILE_NAME','SCAN_NUMBER'])

    

    PROTEINNAME = path.split('/')[-1].split('.')[0]
    # print(PROTEINNAME)


    DOMTree=parse(path)
    feed=DOMTree.documentElement

    ENTRYS = feed.getElementsByTagName('entry')
    for entry in ENTRYS:
        CONTENTS = entry.getElementsByTagName('content')
        for content in CONTENTS:
            PROPERTIES = content.getElementsByTagName('m:properties')
            for properties in PROPERTIES:
                SEARCH_ENGINE = getValue(properties.childNodes[0])
                IDENTIFICATION_ID = getValue(properties.childNodes[1])
                PRECURSOR_MZ = getValue(properties.childNodes[2])
                PRECURSOR_CHARGE = getValue(properties.childNodes[3])
                RECALIBRATED_PRECURSOR_MZ = getValue(properties.childNodes[4])
                SEQUENCE = getValue(properties.childNodes[5])
                SCORE = getValue(properties.childNodes[6])
                DELTA_SCORE = getValue(properties.childNodes[7])
                THRESHOLD_SCORE = getValue(properties.childNodes[8])
                PEP = getValue(properties.childNodes[9])
                PEPTIDE_ID = getValue(properties.childNodes[10])
                Q_VALUE = getValue(properties.childNodes[11])
                LLD_FDR_CUTOFF = getValue(properties.childNodes[12])
                VARIABLE_MODIFICATION_STRING = getValue(properties.childNodes[13])
                FIXED_MODIFICATION_STRING = getValue(properties.childNodes[14])
                MODIFICATION_DELTA_MASS = getValue(properties.childNodes[15])
                EXPECTED_PEPTIDE_MASS = getValue(properties.childNodes[16])
                OBSERVED_PEPTIDE_MASS = getValue(properties.childNodes[17])
                MASS_ERROR_DA = getValue(properties.childNodes[18])
                MASS_ERROR_PPM = getValue(properties.childNodes[19])
                PROJECT_NAME = getValue(properties.childNodes[20])
                EXPERIMENT_NAME = getValue(properties.childNodes[21])
                EXPERIMENT_ID = getValue(properties.childNodes[22])
                FILE_NAME = getValue(properties.childNodes[23])
                SCAN_NUMBER = getValue(properties.childNodes[24])
    #             print(SCAN_NUMBER)

                data.loc[len(data)] = [PROTEINNAME,SEARCH_ENGINE,IDENTIFICATION_ID,PRECURSOR_MZ,PRECURSOR_CHARGE,RECALIBRATED_PRECURSOR_MZ,
                                SEQUENCE,SCORE,DELTA_SCORE,THRESHOLD_SCORE,PEP,PEPTIDE_ID,Q_VALUE,LLD_FDR_CUTOFF,
                                 VARIABLE_MODIFICATION_STRING,FIXED_MODIFICATION_STRING,MODIFICATION_DELTA_MASS,EXPECTED_PEPTIDE_MASS,
                                OBSERVED_PEPTIDE_MASS,MASS_ERROR_DA,MASS_ERROR_PPM,PROJECT_NAME,EXPERIMENT_NAME,
                                EXPERIMENT_ID,FILE_NAME,SCAN_NUMBER]
    data.to_csv('K:/proteomeDB/'+PROTEINNAME+'.csv')
    return data


if __name__=='__main__':
    start = time.time()                                      

    rootdir = 'K:/proteomeDB_proteins/'
    file_list = os.listdir(rootdir)
    protein_list = [os.path.join(rootdir,file_list[i]) for i in range(len(file_list))]
    
    pool = Pool(8)
    result = pool.map(proteomeDB_proteins,protein_list)
    pool.close()
    pool.join()

    testData = pd.DataFrame(columns=['PROTEINNAME','SEARCH_ENGINE','IDENTIFICATION_ID','PRECURSOR_MZ','PRECURSOR_CHARGE','RECALIBRATED_PRECURSOR_MZ',
                            'SEQUENCE','SCORE','DELTA_SCORE','THRESHOLD_SCORE','PEP','PEPTIDE_ID','Q_VALUE','LLD_FDR_CUTOFF',
                             'VARIABLE_MODIFICATION_STRING','FIXED_MODIFICATION_STRING','MODIFICATION_DELTA_MASS','EXPECTED_PEPTIDE_MASS',
                            'OBSERVED_PEPTIDE_MASS','MASS_ERROR_DA','MASS_ERROR_PPM','PROJECT_NAME','EXPERIMENT_NAME',
                            'EXPERIMENT_ID','FILE_NAME','SCAN_NUMBER'])

	for r in result:
    	testData = testData.append(r)
    testData.to_csv('K:/proteomeDB/all.csv')

    end = time.time()
    print('Task %s runs %0.2f seconds.' % (name, (end - start)))
