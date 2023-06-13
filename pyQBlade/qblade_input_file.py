# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 15:59:38 2022

@author: Guido Lazzerini e Giancarlo Troise
"""
from io import open
from builtins import range
from builtins import str
import re
import numpy as np

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    from file import File, WrongFormatError, BrokenFormatError

TABTYPE_NOT_A_TAB = 0
TABTYPE_NUM_WITH_HEADER    = 1
TABTYPE_NUM_NO_HEADER      = 4
TABTYPE_MIX_WITH_HEADER    = 5

class QBladeInputFile(File):
    """ 
    Read/write a QBlade input file. The object behaves like a dictionary.
    
    Only tested with .sim and main and substructure .str!!!
    
    Main methods
    ------------
    - read, write, keys, modTabElm, appendTabLine

    Main keys
    ---------
    The keys correspond to the keys used in the file. For instance for a .str file: 'TWR2SHFT'

    Examples
    --------

        filename = 'DTU10MW_Main.str'
        B=QBladeInputFile('./test_files/DTU10MW_Main.str')
        B['TWR2SHFT'] = 2.0
        B.write('MOD_Main.str')
        
    """
    @staticmethod
    def defaultExtensions():
        return ['sim','.dat','.str','.txt']

    @staticmethod
    def formatName():
        return 'QBlade input file'
    
    def __init__(self, filename=None, **kwargs):

        super(QBladeInputFile, self).__init__(filename=filename,**kwargs)
    
    def keys(self):
        self.labels = [ d['label'] for d in self.data if not d['isComment'] ]
        return self.labels

    def getID(self,label):
        i=self.getIDSafe(label)
        if i<0:
            raise KeyError('Variable `'+ label+'` not found in QBlade file:'+self.filename)
        else:
            return i
        
    def getIDs(self,label):
        I=[]
        # brute force search
        for i in range(len(self.data)):
            d = self.data[i]
            if d['label'].lower()==label.lower():
                I.append(i)
        if len(I)<0:
            raise KeyError('Variable `'+ label+'` not found in QBlade file:'+self.filename)
        else:
            return I

    def getIDSafe(self,label):
        # brute force search
        for i in range(len(self.data)):
            d = self.data[i]
            if d['label'].lower()==label.lower():
                return i
        return -1

    # Making object an iterator
    def __iter__(self):
        self.iCurrent=-1
        self.iMax=len(self.data)-1
        return self

    def __next__(self): # Python 2: def next(self)
        if self.iCurrent > self.iMax:
            raise StopIteration
        else:
            self.iCurrent += 1
            return self.data[self.iCurrent]

    # Making it behave like a dictionary
    def __setitem__(self,key,item):
        I = self.getIDs(key)
        for i in I: 
            self.data[i]['value'] = item

    def __getitem__(self,key):
        i = self.getID(key)
        return self.data[i]['value']

    def __repr__(self):
        s ='QBlade input file: {}\n'.format(self.filename)
        return s+'\n'.join(['{:15s}: {} (comment line: {}, description: {})'.format(d['label'],d['value'],d['isComment'],d['descr']) for i,d in enumerate(self.data)])

    def addKeyVal(self,key,val,descr=None):
        d=getDict()
        d['label']=key
        d['value']=val
        if descr is not None:
            d['descr']=descr
        self.data.append(d)
    
    def modTabElm(self,label,row,colName,newValue):
        i = self.getID(label)
        if type(row) == int:
            row = [row]
        for k in row:
                        
            if k <= len(self.data[i]['value']):
                
                if type(colName) == int:
                    try:
                        self.data[i]['value'][k][colName] = str(newValue)
                    except:
                        KeyError('Column `'+colName+'` not found in '+label)
                else:
                    colIndex =np.where(np.array(self.data[i]['tabColumnNames'])==colName)
                    
                    if len(colIndex[0])!=0:
                        #print('Modifying value at row '+str(k)+' and column '+colName+'.')
                        
                        if type(newValue) == float:
                            newValue = "%.4f" % (newValue)
                            
                        self.data[i]['value'][k][colIndex] = str(newValue)
                    else:
                        raise KeyError('Column `'+colName+'` not found in '+label)    
            else:
                raise KeyError('Row '+str(k)+' not found in '+label)
    
    def appendTabLine(self,label,line):
        i = self.getID(label)
        self.data[i]['value'] = np.append(self.data[i]['value'],[np.array(line)], axis=0)
                
    
    def _read(self,**kwargs):
        
        # DIM_VAR = 0 if uknown a priori
        # --- Tables that can be detected based on 'Label' (first entry on line)
        # Sub File
        NUMTAB_FROM_VAL_VARNAME = ['TP_INTERFACE_POS' , 'REF_COG_POS' , 'REF_HYDRO_POS' , 'SUB_MASS' , 'SUB_HYDROSTIFFNESS']
        NUMTAB_FROM_VAL_DIM_VAR = [1                  , 1             , 1               , 6          , 6                   ] 
        NUMTAB_FROM_VAL_NHEADER = [1                  , 1             , 1               , 0          , 1                   ]
        NUMTAB_FROM_VAL_TYPE    = ['num'              , 'num'         , 'num'           , 'num'      , 'num'               ]
        # # Sub File continues
        NUMTAB_FROM_VAL_VARNAME += ['SUB_HYDRODAMPING' , 'SUB_HYDROQUADDAMPING', 'SUB_HYDROADDEDMASS' , 'SUB_HYDROCONSTFORCE' , 'MARINEGROWTH' ]
        NUMTAB_FROM_VAL_DIM_VAR += [6                  , 6                     , 6                    , 1                     , 1              ]
        NUMTAB_FROM_VAL_NHEADER += [1                  , 1                     , 0                    , 1                     , 1              ]
        NUMTAB_FROM_VAL_TYPE    += ['num'              , 'num'                 , 'num'                , 'num'                 , 'num'          ]
        # # # Sub File continues
        NUMTAB_FROM_VAL_VARNAME += ['HYDROJOINTCOEFF'  ,'HYDROMEMBERCOEFF' , 'SUBELEMENTS', 'SUBELEMENTSRIGID' , 'JOINTOFFSET' ]
        NUMTAB_FROM_VAL_DIM_VAR += [0                  , 0                 , 0            ,0                   ,  1             ]
        NUMTAB_FROM_VAL_NHEADER += [1                  , 1                 , 1            ,1                   ,  0             ]
        NUMTAB_FROM_VAL_TYPE    += ['num'              , 'num'             , 'num'        ,'num'               , 'num'         ]
        # # # # Sub File continues
        NUMTAB_FROM_VAL_VARNAME += ['SUBJOINTS' , 'SUBCONSTRAINTS' , 'SUBMEMBERS' , 'MOORELEMENTS' , 'MOORMEMBERS' ]
        NUMTAB_FROM_VAL_DIM_VAR += [0           , 0                , 0            , 0              , 0             ]
        NUMTAB_FROM_VAL_NHEADER += [1           , 1                , 1            , 1              , 1             ]
        NUMTAB_FROM_VAL_TYPE    += ['num'       , 'num'            , 'mix'        , 'num'          , 'mix'         ]
        # # # # Sub File continues
        NUMTAB_FROM_VAL_VARNAME += ['TRANSITIONBLOCK' , 'TRANSITIONCYLINDER' , 'RGBCOLOR']
        NUMTAB_FROM_VAL_DIM_VAR += [1                 , 1                    , 1         ]
        NUMTAB_FROM_VAL_NHEADER += [1                 , 1                    , 1         ]
        NUMTAB_FROM_VAL_TYPE    += ['num'             , 'num'                , 'num'     ]
        
        # Using lower case to be more tolerant..
        NUMTAB_FROM_VAL_VARNAME_L = [s.lower() for s in NUMTAB_FROM_VAL_VARNAME]
        
        self.data   = []
        self.hasNodal=False
        self.module = None
        with open(self.filename, 'r', errors="surrogateescape") as f:
            lines=f.read().splitlines()

        # Parsing line by line, storing each line into a dictionary
        i=0    
        #nComments  = 0
        #nWrongLabels = 0
        allowSpaceSeparatedList=False
        while i<len(lines):
            # --- Parsing of standard lines: value(s) key comment
            line = lines[i]
            d = parseQBladeInputLine(line,i,allowSpaceSeparatedList)
            
            # --- Handling of tables
            if isStr(d['value']) and d['value'].lower() in NUMTAB_FROM_VAL_VARNAME_L:
                # Table with numerical values, 
                ii             = NUMTAB_FROM_VAL_VARNAME_L.index(d['value'].lower())
                tab_type       = NUMTAB_FROM_VAL_TYPE[ii]
                if tab_type=='num':
                    d['tabType']   = TABTYPE_NUM_WITH_HEADER
                else:
                    d['tabType']   = TABTYPE_MIX_WITH_HEADER
                d['label']     = NUMTAB_FROM_VAL_VARNAME[ii].upper()
                d['tabDimVar'] = NUMTAB_FROM_VAL_DIM_VAR[ii]
                nHeaders       = NUMTAB_FROM_VAL_NHEADER[ii]
                nTabLines=0
                if d['tabDimVar']!= 0:
                    nTabLines = d['tabDimVar']
                else:
                    iStart = i + nHeaders + 1
                    k = iStart
                    for k in range(iStart,len(lines)):
                        #print('k is now:'+str(k))
                        #print(lines[k])
                        try:
                            dummy = lines[k][0]
                        except:
                            nTabLines = k - iStart
                            d['tabDimVar'] = k-iStart
                            break
                    
                #print('Reading table {} Dimension {} (based on {})'.format(d['label'],d['tabDimVar'],d['tabDimVar']));
                d['value'], d['tabColumnNames'], d['tabUnits'] = parseQBladeNumTable(self.filename,lines[i:i+nTabLines+nHeaders+1], nTabLines, i, nHeaders, tableType=tab_type, varNumLines=d['tabDimVar'])
                i += nTabLines+nHeaders+1
                
                #GT ... 
                # check if d is a matrix and set to "not a comment"
                if d['tabType']!=TABTYPE_NOT_A_TAB:
                    d['isComment']=False
                #... GT 
                
            self.data.append(d)
            i += 1
            
    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString())
    
    def toString(self):
        
        s=''
        
        def toStringVLD(val,lab,descr):
            val='{}'.format(val)
            lab='{}'.format(lab)
            if len(val)<13:
                val='{:13s}'.format(val)
            if len(lab)<13:
                lab='{:13s}'.format(lab)
            return val+' '+lab+' - '+descr.strip().strip('-').strip()+'\n'
   
        for i in range(len(self.data)):
            d=self.data[i]

            if type(d['value'])==bool:
                d['value'] = str(d['value']).lower() # to avoid 'True' and 'False'
                
            if d['isComment']:
                s+='{}'.format(d['value'])
            elif d['tabType']==TABTYPE_NOT_A_TAB:
                if isinstance(d['value'], list):
                    sList=', '.join([str(x) for x in d['value']])
                    s+='{} {} {}'.format(sList,d['label'],d['descr'])
                else:
                    s+=toStringVLD(d['value'],d['label'],d['descr']).strip()
            elif d['tabType']==TABTYPE_NUM_WITH_HEADER:
                #GT... add a newline at the beginning, some lines were not separated
                # s+=d['label']+'\n'
                s+='\n'+d['label']+'\n'
                #...GT
                if d['tabColumnNames'] is not None:
                    s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                #s+=d['descr'] # Not ready for that
                    if d['tabUnits'] is not None:
                        s+='\n'
                        s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                    newline='\n'
                else:
                    newline=''
                if np.size(d['value'],0) > 0 :
                    s+=newline
                    s+='\n'.join('\t'.join( ('{:15.0f}'.format(x) if int(x)==x else '{:15.8e}'.format(x) )  for x in y) for y in d['value'])
            elif d['tabType']==TABTYPE_MIX_WITH_HEADER:
                #GT... write \n + table name
                s+='\n'+d['label']+'\n'
                #...GT
                s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                if d['tabUnits'] is not None:
                    s+='\n'
                    s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                if np.size(d['value'],0) > 0 :
                    s+='\n'
                    s+='\n'.join('\t'.join('{}'.format(x) for x in y) for y in d['value'])            
            else:
                raise Exception('Unknown table type for variable {}'.format(d))
            if i<len(self.data)-1:
                s+='\n'
        return s

    def _toDataFrame(self):
        print('Not implemented yet')
        pass

    def _fromDataFrame(self):
        print('Not implemented yet')
        pass
    
    def _fromDictionary(self):
        print('Not implemented yet')
        pass

    def _fromFile(self):
        print('Not implemented yet')
        pass
    
def isStr(s):
    # Python 2 and 3 compatible
    # Two options below
    # NOTE: all this avoided since we import str from builtins
    # --- Version 2
    #     isString = False;
    #     if(isinstance(s, str)):
    #         isString = True;
    #     try:
    #         if(isinstance(s, basestring)): # todo unicode as well
    #             isString = True;
    #     except NameError:
    #         pass; 
    #     return isString
    # --- Version 1
    #     try: 
    #        basestring # python 2
    #        return isinstance(s, basestring) or isinstance(s,unicode)
    #     except NameError:
    #          basestring=str #python 3
    #     return isinstance(s, str)
   return isinstance(s, str)
    
def strIsFloat(s):
    #return s.replace('.',',1').isdigit()
    try:
        float(s)
        return True
    except:
        return False

def strIsBool(s):
    return (s.lower() == 'true') or (s.lower() == 'false')

def strIsInt(s):
    s = str(s)
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()    

def hasSpecialChars(s):
    # fast allows for parenthesis
    # For now we allow for - but that's because of BeamDyn geometry members 
    return not re.match("^[\"\'a-zA-Z0-9_()-]*$", s)

def cleanLine(l):
    # makes a string single space separated
    l = l.replace('\t',' ')
    l = ' '.join(l.split())
    l = l.strip()
    return l

def cleanAfterChar(l,c):
    # remove whats after a character
    n = l.find(c);
    if n>0:
        return l[:n]
    else:
        return l
    
def detectUnits(s,nRef):
    nPOpen=s.count('(')
    nPClos=s.count(')')
    nBOpen=s.count('[')
    nBClos=s.count(']')

    sep='!#@#!'
    if (nPOpen == nPClos) and (nPOpen>=nRef):
        #that's pretty good
        Units=s.replace('(','').replace(')',sep).split(sep)[:-1]
    elif (nBOpen == nBClos) and (nBOpen>=nRef):
        Units=s.replace('[','').replace(']',sep).split(sep)[:-1]
    else:
        Units=s.split()
    return Units

def getDict():
    return {'value':None, 'label':'', 'isComment':False, 'descr':'', 'tabType':TABTYPE_NOT_A_TAB}

def parseQBladeInputLine(line_raw,i,allowSpaceSeparatedList=False):
    d = getDict()
    #print(line_raw)
    try:
        # preliminary cleaning (Note: loss of formatting)
        line = cleanLine(line_raw)
        # Generic Comment
        if any(line.startswith(c) for c in ['#','!','-','==','/']) or len(line)==0:
            d['isComment']=True
            d['value']=line_raw
            return d
        # .Sim Comments
        if any(line.startswith(c) for c in ['Generated','Archive','Time','Date','multiple','the']) or len(line)==0:
            d['isComment']=True
            d['value']=line_raw
            return d
        if line.lower().startswith('end'):
            sp =line.split()
            if len(sp)>2 and sp[1]=='of':
                d['isComment']=True
                d['value']=line_raw
        # TODO add detecting output lists of QBlade (now treated as comments)
        # Detecting lists
        List=[];
        iComma=line.find(',')
        if iComma>0 and iComma<30:
            fakeline=line.replace(' ',',')
            fakeline=re.sub(',+',',',fakeline)
            csplits=fakeline.split(',')
            # Splitting based on comma and looping while it's numbers or booleans
            ii=0
            s=csplits[ii]
            #print(csplits)
            while strIsFloat(s) or strIsBool(s) and ii<len(csplits):
                if strIsInt(s):
                    List.append(int(s))
                elif strIsFloat(s):
                    List.append(float(s))
                elif strIsBool(s):
                    List.append(bool(s))
                else:
                    raise WrongFormatError('Lists of strings not supported.')
                ii =ii+1
                if ii>=len(csplits):
                    raise WrongFormatError('Wrong number of list values')
                s = csplits[ii]
            #print('[INFO] Line {}: Found list: '.format(i),List)
        # Defining value and remaining splits
        if len(List)>=2:
            d['value']=List
            line_remaining=line
            # eating line, removing each values
            for iii in range(ii):
                sValue=csplits[iii]
                ipos=line_remaining.find(sValue)
                line_remaining = line_remaining[ipos+len(sValue):]
            splits=line_remaining.split()
            iNext=0
        else:
            # It's not a list, we just use space as separators
            line.replace('//',' ')
            splits=line.split(' ')
            s=splits[0]

            if strIsInt(s):
                d['value']=int(s)
                if allowSpaceSeparatedList and len(splits)>1:
                    if strIsInt(splits[1]):
                        d['value']=splits[0]+ ' '+splits[1]
            elif strIsFloat(s):
                d['value']=float(s)
            elif strIsBool(s):
                d['value']=bool(s)
            else:
                d['value']=s
            iNext=1
            #import pdb  ; pdb.set_trace();

        # Extracting label (TODO, for now only second split)
        bOK=False
        while (not bOK) and iNext<len(splits):
            # Nasty handling of !XXX: comments
            if splits[iNext][0]=='!' and splits[iNext][-1]==':': 
                iNext=iNext+2
                continue
            # Nasty handling of the fact that sometimes old values are repeated before the label
            if strIsFloat(splits[iNext]):
                iNext=iNext+1
                continue
            else:
                bOK=True
        if bOK:
            d['label']= splits[iNext].strip()
            iNext = iNext+1
        else:
            #print('[WARN] Line {}: No label found -> comment assumed'.format(i+1))
            d['isComment']=True
            d['value']=line_raw
            iNext = len(splits)+1
                
        # Recombining description
        if len(splits)>=iNext+1:
            d['descr']=' '.join(splits[(iNext):])
    except WrongFormatError as e:
        raise WrongFormatError('Line {}: '.format(i+1)+e.args[0])
    except Exception as e:
        raise Exception('Line {}: '.format(i+1)+e.args[0])

    return d

def parseQBladeNumTable(filename,lines,n,iStart,nHeaders=2,tableType='num',nOffset=1, varNumLines=''):
    """ 
    First lines of data starts at: nHeaders+nOffset
    
    """
    Tab = None
    ColNames = None
    Units = None
    

    if len(lines)!=n+nHeaders+nOffset:
        raise BrokenFormatError('Not enough lines in table: {} lines instead of {}\nFile:{}'.format(len(lines)-nHeaders,n,filename))
    try:
        if nHeaders==0:
            # Extract number of values from number of numerical values on first line
            numeric_const_pattern = r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
            rx = re.compile(numeric_const_pattern, re.VERBOSE)
            header = cleanAfterChar(lines[nOffset], '!')
            if tableType=='num':
                dat= np.array(rx.findall(header)).astype(float)
                ColNames=['C{}'.format(j) for j in range(len(dat))]
            else:
                raise NotImplementedError('Reading QBlade tables with no headers for type different than num not implemented yet')

        elif nHeaders>=1:
            # Extract column names
            i = nOffset
            sTmp = cleanLine(lines[i])
            # sTmp = cleanAfterChar(sTmp,'[')
            # sTmp = cleanAfterChar(sTmp,'(')
            sTmp = cleanAfterChar(sTmp,'!')
            sTmp = cleanAfterChar(sTmp,'#')
            if sTmp.startswith('!'):
                sTmp=sTmp[1:].strip()
            ColNames=sTmp.split()
        # if nHeaders>=2:
        #     # Extract units
        #     i = 1
        #     sTmp = cleanLine(lines[i])
        #     sTmp = cleanAfterChar(sTmp,'!')
        #     sTmp = cleanAfterChar(sTmp,'#')
        #     if sTmp.startswith('!'):
        #         sTmp=sTmp[1:].strip()

        #     Units = detectUnits(sTmp,len(ColNames))
        #     Units = ['({})'.format(u.strip()) for u in Units]
        #     # Forcing user to match number of units and column names
        #     if len(ColNames) != len(Units):
        #         print(ColNames)
        #         print(Units)
        #         print('[WARN] {}: Line {}: Number of column names different from number of units in table'.format(filename, iStart+i+1))

        nCols=len(ColNames)

        if tableType=='num':
            if n==0:
                Tab = np.zeros((n, nCols))
            for i in range(nHeaders+nOffset,n+nHeaders+nOffset):
                l = cleanAfterChar(lines[i].lower(),'!')
                l = cleanAfterChar(l,'#')
                l = cleanAfterChar(l,'/')
                v = l.split()
                if len(v) != nCols:
                    # Discarding SubDyn special cases
                    if ColNames[-1].lower() not in ['nodecnt']:
                        print('[WARN] {}: Line {}: number of data different from number of column names. ColumnNames: {}'.format(filename, iStart+i+1, ColNames))
                if i==nHeaders+nOffset:
                    # Node Cnt
                    if len(v) != nCols:
                        if ColNames[-1].lower()== 'nodecnt':
                            ColNames = ColNames+['Col']*(len(v)-nCols)
                            Units    = Units+['Col']*(len(v)-nCols)

                    nCols=len(v)
                    Tab = np.zeros((n, nCols))
                # Accounting for TRUE FALSE and converting to float
                v = [s.replace('true','1').replace('false','0').replace('noprint','0').replace('print','1') for s in v]
                v = [float(s) for s in v[0:nCols]]
                if len(v) < nCols:
                    raise Exception('Number of data is lower than number of column names')
                Tab[i-nHeaders-nOffset,:] = v
        elif tableType=='mix':
              # a mix table contains a mixed of strings and floats
              # For now, we are being a bit more relaxed about the number of columns
              if n==0:
                  Tab = np.zeros((n, nCols)).astype(object)
              for i in range(nHeaders+nOffset,n+nHeaders+nOffset):
                  l = lines[i]
                  l = cleanAfterChar(l,'!')
                  l = cleanAfterChar(l,'#')
                  v = l.split()
                  if l.startswith('---'):
                      raise BrokenFormatError('Error reading line {} while reading table. Is the variable `{}` set correctly?'.format(iStart+i+1, varNumLines))
                  if len(v) != nCols:
                      # Discarding SubDyn special cases
                      if ColNames[-1].lower() not in ['cosmid', 'ssifile']:
                          print('[WARN] {}: Line {}: Number of data is different than number of column names. Column Names: {}'.format(filename,iStart+1+i, ColNames))
                  if i==nHeaders+nOffset:
                      if len(v)>nCols:
                          ColNames = ColNames+['Col']*(len(v)-nCols)
                          Units    = Units+['Col']*(len(v)-nCols)
                      nCols=len(v)
                      Tab = np.zeros((n, nCols)).astype(object)
                  v=v[0:min(len(v),nCols)]
                  Tab[i-nHeaders-nOffset,0:len(v)] = v
              # If all values are float, we convert to float
              if all([strIsFloat(x) for x in Tab.ravel()]):
                  Tab=Tab.astype(float)      
        else:
            raise Exception('Unknown table type')

        ColNames = ColNames[0:nCols]
        if Units is not None:
            Units    = Units[0:nCols]
            Units    = ['('+u.replace('(','').replace(')','')+')' for u in Units]
        if nHeaders==0:
            ColNames=None
                
    except Exception as e:    
        raise BrokenFormatError('Line {}: {}'.format(iStart+i+1,e.args[0]))
    return Tab, ColNames, Units

if __name__ == "__main__":
    
    B=QBladeInputFile(r'C:/test_qbladedll_and_io/pyQBlade/test_files/DTU10MW_Blade.str')
    # B.modTabElm('MOORMEMBERS',0,'Len.[m]',910)
    #B.appendTabLine('SUBELEMENTS',[10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 ])
    #B.modTabElm('SUBELEMENTS',0,'Eiy_[N.m^2]',90000)
    #B=QBladeInputFile(r'D:\prove\test_qbladedll_and_io\simulation_files\NREL_5MW_AE_hydro_3moor_simple\Structure\OC3_Sparbuoy_Sub_LPMD.str')
    #B['TIMESTEP'] = 0.05
    #B['SUBMEMBERS'] = np.append(B['SUBMEMBERS'],['12', '11', '12', '0', '9', '1', '0', '0', '0', '5', 'UpperBuoy2'])
    print(B.keys())
    #B.addKeyVal('','MOO_3_1.0',)
    #B.write(r'prova5.str')

