# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:36:46 2022

@author: Guido Lazzerini
"""

class TurbModel(object):
    
    def __init__(self,filename):

        self.filename = filename
        self.data = []
        lines = self._read()
        i = 0
        
        while i<len(lines):
            # --- Parsing of standard lines: value - key 
            line_raw = lines[i]
            d = getDict()
            #print(line_raw)
            try:
                # Preliminary cleaning (Note:loss of formatting)
                line = cleanLine(line_raw)
                # Generic Comment
                if any(line.startswith(c) for c in ['#','!']) or len(line)==0:
                    d['isComment']=True
                    d['value']=line_raw
                    i +=1
                    self.data.append(d)
                else:
                    splits=line.split(' ')
                    s = splits[0]
                    s1 = splits[2]
                    d['value'] = s
                    d['label'] =  s1
                    d['isComment'] = False
                    i +=1
                    self.data.append(d)
            except:
                raise print('Something wrong happened @ line {}: '.format(i+1))
                
        # Check if QBlade or OpenFAST model file
        try: 
            a = self.getID('FSTFILENAME')
            d = getDict()
            d['value'] = 'OpenFAST'
            d['label'] = 'SIMSOFTWARE'
            d['isComment'] = False
            self.data.append(d)
        except:
            try:
                a = self.getID('SIMFILENAME')
                d = getDict()
                d['value'] = 'QBlade'
                d['label'] = 'SIMSOFTWARE'
                d['isComment'] = False
                self.data.append(d)                
            except KeyError:
                print('Definition file does not contain any .fst or .sim filenames, no simulation will be performed if I continue...')
    
    # Making it behave like a dictionary
    def __setitem__(self,key,item):
        I = self.getIDs(key)
        for i in I: 
            self.data[i]['value'] = item

    def __getitem__(self,key):
        i = self.getID(key)
        return self.data[i]['value']
    
    def getID(self,label):
        i=self.getIDSafe(label)
        if i<0:
            raise KeyError('Variable `'+ label+'` not found in file:'+self.filename)
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
            raise KeyError('Variable `'+ label+'` not found in file:'+self.filename)
        else:
            return I
    
    def getIDSafe(self,label):
        # Brute force search
        for i in range(len(self.data)):
            d = self.data[i]
            if d['label'].lower()==label.lower():
                return i
        return -1
    
    def addKeyVal(self,key,val):
        d=getDict()
        d['label']=key
        d['value']=val
        self.data.append(d)
    
    def fromQBladeModel(self,file):
        return
    
    def fromOpenFASTModel(self,file):
        return
    
    def _read(self):
        with open(self.filename, 'r', errors="surrogateescape") as f:
            lines=f.read().splitlines()
        return lines
    
def getDict():
    return {'value':None, 'label':'', 'isComment':False}

def cleanLine(l):
    # makes a string single space separated
    l = l.replace('\t',' ')
    l = ' '.join(l.split())
    l = l.strip()
    return l

if __name__ == "__main__":
    
    model = TurbModel(r'C:\FOWT_optim\FOWT_optim_test\sims\template_input_files\_DTU10MW3Spar_modeldefinition.dat')