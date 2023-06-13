import numpy as np
import pandas as pd
import struct 
import os
import re

try:
    from .file import File, WrongFormatError, BrokenReaderError
except:
    from file import File, WrongFormatError, BrokenReaderError
    
try:
    from .csv_file import CSVFile
except:
    from csv_file import CSVFile
# --------------------------------------------------------------------------------}
# --- OUT FILE 
# --------------------------------------------------------------------------------{
class QBladeOutputFile(File):
    """ 
    Read a QBlade ouput file (.outq).

    Main methods
    ------------
    - read, toDataFrame

    Examples
    --------

        # read an output file and convert it to pandas dataframe
        df = QBladeOutputFile('Output.outq').toDataFrame()
        time  = df['Time_[s]']
        Omega = df['RotSpeed_[rpm]']

    """

    @staticmethod
    def defaultExtensions():
        return ['.outq']

    @staticmethod
    def formatName():
        return 'QBlade output file'

    def _read(self):
        
        self.info={}
        self['binary']=False
        
        try:
                F=CSVFile(filename=self.filename, sep=',', commentLines=[0,2],colNamesLine=2)
                self.data = F.data
                del F
                self.info['attribute_units']=None
                self.info['attribute_names']=self.data.columns.values

        except MemoryError as e:    
            raise BrokenReaderError('QBlade Out File {}: Memory error encountered\n{}'.format(self.filename,e))
        except Exception as e:    
            raise WrongFormatError('QBlade Out File {}: {}'.format(self.filename,e.args))

        if self.info['attribute_units'] is not None:
            self.info['attribute_units'] = [re.sub(r'[()\[\]]','',u) for u in self.info['attribute_units']]


    def _write(self): 
        with open(self.filename,'w') as f:
            f.write(','.join(['{:>10s}'.format(c)         for c in self.info['attribute_names']])+'\n')
            #f.write('\t'.join(['{:>10s}'.format('('+u+')') for u in self.info['attribute_units']])+'\n')
            # TODO better..
            f.write('\n'.join([','.join(['{:10.4f}'.format(y[0])]+['{:10.3e}'.format(x) for x in y[1:]]) for y in self.data]))

    def _toDataFrame(self):
        if self.info['attribute_units'] is not None:
            cols=[n+'_['+u.replace('sec','s')+']' for n,u in zip(self.info['attribute_names'],self.info['attribute_units'])]
        else:
            cols=self.info['attribute_names']
        if isinstance(self.data, pd.DataFrame):
            df= self.data
            df.columns=cols
        else:
            df = pd.DataFrame(data=self.data,columns=cols)

        return df


# --------------------------------------------------------------------------------
# --- Helper low level functions 
# --------------------------------------------------------------------------------


if __name__ == "__main__":
    B=QBladeOutputFile('outputSim.outq')
    print(B.data)



